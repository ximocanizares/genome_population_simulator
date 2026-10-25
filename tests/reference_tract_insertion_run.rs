//! F1.7.2: la inserción de tractos **conectada al producto** (REQ-REF-37,
//! DTC-098).
//!
//! `tests/reference_tract_insertion.rs` prueba el módulo que empalma; este
//! fichero prueba lo que lo rodea: la declaración en el TOML, la traducción de
//! los demás loci y regiones, la publicación de la referencia final con su
//! mapa, que la verdad, el VCF y la matriz hablen en sus coordenadas, los
//! rechazos de esta primera versión, y el juez independiente.
//!
//! Los oráculos son aritméticos y están escritos en cada prueba: «la base 2001
//! de partida acaba en la 2031 tras insertar 30 bases después de la 1000». El
//! genoma de partida es pseudoaleatorio a propósito: con un `ACGT` periódico,
//! un desplazamiento de cuatro bases sería invisible al comparar.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::{
    ExecutionAuthorizations, execute_validated_config, execute_validated_config_with_authorizations,
};
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::results::validate_published_resources;
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "00000000000000000000000000000000000000000000000000000000000000aa";
const TRACT_30: &str = "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG";

fn temporary(label: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "gps-tract-insertion-run-{label}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&path).expect("directorio temporal");
    path
}

/// Genoma de partida de 4000 bases: pseudoaleatorio salvo un tracto `(CA)20`
/// en 2001..2040, que es el locus de tracto existente de los escenarios.
fn source_bases() -> String {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut noisy = |length: usize| -> String {
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                char::from(b"ACGT"[(state % 4) as usize])
            })
            .collect()
    };
    let bases = format!("{}{}{}", noisy(2000), "CA".repeat(20), noisy(1960));
    assert_eq!(bases.len(), 4000);
    bases
}

/// FASTA en líneas de 60, distinto del canónico de 80: la copia suministrada
/// se publica byte a byte y la final en canónico, y la prueba no puede
/// confundirlas.
fn write_fasta(path: &Path, bases: &str) {
    let lines = bases
        .as_bytes()
        .chunks(60)
        .map(|chunk| String::from_utf8(chunk.to_vec()).expect("ASCII"))
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(path, format!(">chr1\n{lines}\n")).expect("FASTA de partida");
}

const MATRIX_GENERATE: &str = r#"
[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 3
spatial_mode = "weighted_regions"
outside_weight = 0.0
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.regions]]
contig = "chr1"
start = 3000
end = 3500
weight = 1.0

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
"#;

const POPULATION_AND_VCF: &str = r#"
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"

[vcf_output]
mode = "bgzf"
"#;

/// Locus que inserta `(CAG)10` tras la base 1000 de partida, con alelos
/// declarados para que el sorteo no dependa del generador.
const INSERT_AT_1000: &str = r#"
[[reference.repeat_loci]]
id = "STR-INS"
sequence = "chr1"
insert_after = 1000
motif = "CAG"
reference_units = [10]
[[reference.repeat_loci.alleles]]
units = 10
frequency = 0.5
[[reference.repeat_loci.alleles]]
units = 11
frequency = 0.5
"#;

/// Locus de tracto existente, declarado en coordenadas de partida.
const EXISTING_CA: &str = r#"
[[reference.repeat_loci]]
id = "STR-CA"
sequence = "chr1"
start = 2001
end = 2040
motif = "CA"
[[reference.repeat_loci.alleles]]
units = 20
frequency = 0.5
[[reference.repeat_loci.alleles]]
units = 21
frequency = 0.5
"#;

/// Documento general con referencia suministrada; `reference_extra` va dentro
/// de `[reference]` y `loci` después del catálogo.
fn document(reference_extra: &str, loci: &str, rest: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"
{MATRIX_GENERATE}
[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"
{reference_extra}

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 4000
{loci}
{rest}"#
    )
}

fn with_insertion(loci: &str) -> String {
    document(
        "inserted_reference_id = \"ref-main-str\"",
        loci,
        POPULATION_AND_VCF,
    )
}

struct Scenario {
    root: PathBuf,
    destination: PathBuf,
}

impl Scenario {
    fn read(&self, relative: &str) -> Vec<u8> {
        fs::read(self.destination.join(relative)).expect("artefacto publicado")
    }

    fn json(&self, relative: &str) -> Value {
        serde_json::from_slice(&self.read(relative)).expect("JSON publicado")
    }

    /// Bases de un FASTA publicado, sin cabecera ni saltos de línea.
    fn unfolded(&self, relative: &str) -> String {
        String::from_utf8(self.read(relative))
            .expect("FASTA ASCII")
            .lines()
            .filter(|line| !line.starts_with('>'))
            .collect()
    }
}

fn write_config(body: &str) -> (PathBuf, PathBuf) {
    let root = temporary("config");
    write_fasta(&root.join("reference.fa"), &source_bases());
    let config = root.join("general.toml");
    fs::write(&config, body).expect("configuración");
    (root, config)
}

fn publish(body: &str) -> Scenario {
    let (root, config) = write_config(body);
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Scenario { root, destination }
}

fn rejection(body: &str) -> String {
    let (root, config) = write_config(body);
    let error = load_and_validate_config(&config).expect_err("configuración rechazada");
    fs::remove_dir_all(root).ok();
    format!("{error:?}")
}

fn effective_loci(body: &str) -> Vec<(String, u64, u64)> {
    let (root, config) = write_config(body);
    let effective = load_and_validate_config(&config).expect("configuración válida");
    fs::remove_dir_all(root).ok();
    effective
        .repeat_loci()
        .iter()
        .map(|locus| {
            (
                locus.id().to_owned(),
                locus.start().expect("inicio"),
                locus.end().expect("final"),
            )
        })
        .collect()
}

/// Desplazamiento en base cero donde `grep -ob` encuentra el tracto en la
/// secuencia desplegada. Es un buscador **externo**, ajeno al productor.
fn external_offsets(sequence: &str, needle: &str) -> Vec<u64> {
    let root = temporary("grep");
    let path = root.join("unfolded.txt");
    fs::write(&path, sequence).expect("secuencia desplegada");
    let output = Command::new("grep")
        .arg("-ob")
        .arg(needle)
        .arg(&path)
        .output()
        .expect("grep");
    fs::remove_dir_all(root).ok();
    String::from_utf8(output.stdout)
        .expect("salida de grep")
        .lines()
        .map(|line| {
            line.split(':')
                .next()
                .expect("desplazamiento")
                .parse()
                .expect("número")
        })
        .collect()
}

fn external_sha256(path: &Path) -> String {
    let output = Command::new("sha256sum")
        .arg(path)
        .output()
        .expect("sha256sum");
    String::from_utf8(output.stdout)
        .expect("salida")
        .split_whitespace()
        .next()
        .expect("huella")
        .to_owned()
}

fn validates(document: &Value, schema: &str) -> bool {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let dir = temporary("schema");
    let path = dir.join("document.json");
    fs::write(
        &path,
        serde_json::to_vec_pretty(document).expect("serialización"),
    )
    .expect("documento");
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(&path)
        .output()
        .expect("validador externo de JSON Schema");
    fs::remove_dir_all(dir).ok();
    output.status.success()
}

// ---------------------------------------------------------------------------
// Coordenadas: la traducción al genoma final
// ---------------------------------------------------------------------------

/// **El caso de F1.7.2, de extremo a extremo.** Treinta bases tras la 1000: el
/// tracto aparece en el desplazamiento 1000 de la secuencia final, buscado con
/// `grep` y no con el productor; la base 2001 de partida está en la 2031 del
/// final; y el locus de tracto existente, declarado en 2001..2040, se publica
/// en 2031..2070.
///
/// *Una sola inserción no discrimina* que la segunda posición se lea en el
/// genoma intermedio: no hay ninguna anterior que la haya desplazado. Eso lo
/// caza `the_second_insertion_still_refers_to_the_original_genome`.
#[test]
fn thirty_bases_after_1000_move_base_2001_to_2031() {
    let scenario = publish(&with_insertion(&format!("{INSERT_AT_1000}{EXISTING_CA}")));
    let source = source_bases();
    let derived = scenario.unfolded("resources/references/ref-main-str.fa");
    assert_eq!(derived.len(), 4030);
    assert_eq!(external_offsets(&derived, TRACT_30), vec![1000]);
    // Base uno 2001 de partida = índice 2000; base uno 2031 final = índice 2030.
    assert_eq!(&derived[2030..2031], &source[2000..2001]);
    assert_eq!(&derived[..1000], &source[..1000]);
    assert_eq!(&derived[1030..], &source[1000..]);
    // La copia suministrada se publica tal cual, sin el tracto.
    assert_eq!(
        scenario.read("resources/references/ref-main.fa"),
        fs::read(scenario.root.join("reference.fa")).expect("original")
    );

    let effective = scenario.json("effective-config.json");
    let manifest = scenario.json("manifest.json");
    for (locus, start, end, declared) in [
        (
            "STR-INS",
            1001,
            1030,
            serde_json::json!({"insert_after": 1000}),
        ),
        (
            "STR-CA",
            2031,
            2070,
            serde_json::json!({"start": 2001, "end": 2040}),
        ),
    ] {
        let published = effective["repeat_loci"]
            .as_array()
            .expect("loci")
            .iter()
            .find(|candidate| candidate["id"] == locus)
            .expect("locus publicado");
        assert_eq!(published["start"], start, "{locus}");
        assert_eq!(published["end"], end, "{locus}");
        assert_eq!(published["declared"], declared, "{locus}");
        let in_manifest = manifest["repeat_loci"]["loci"]
            .as_array()
            .expect("loci del manifiesto")
            .iter()
            .find(|candidate| candidate["id"] == locus)
            .expect("locus en el manifiesto");
        assert_eq!(in_manifest["start"], start, "{locus}");
        assert_eq!(in_manifest["end"], end, "{locus}");
    }
    // Las bases del locus existente siguen siendo el (CA)20 en su sitio final.
    assert_eq!(&derived[2030..2070], "CA".repeat(20));
}

/// **La prueba que discrimina la regla de las coordenadas de partida.** Dos
/// inserciones, tras la 1000 (30 bases) y tras la 2000 (20 bases, declaradas en
/// el genoma original): el segundo tracto entra en el desplazamiento 2030 del
/// final, no en el 2000 que daría leer la segunda posición sobre el genoma ya
/// modificado. El locus que la declara se publica en 2031..2050.
#[test]
fn the_second_insertion_still_refers_to_the_original_genome() {
    let second = r#"
[[reference.repeat_loci]]
id = "STR-GATA"
sequence = "chr1"
insert_after = 2050
motif = "GATA"
reference_units = [5]
[[reference.repeat_loci.alleles]]
units = 5
frequency = 1.0
"#;
    let scenario = publish(&with_insertion(&format!("{second}{INSERT_AT_1000}")));
    let derived = scenario.unfolded("resources/references/ref-main-str.fa");
    let source = source_bases();
    assert_eq!(derived.len(), 4050);
    assert_eq!(external_offsets(&derived, TRACT_30), vec![1000]);
    assert_eq!(external_offsets(&derived, &"GATA".repeat(5)), vec![2080]);
    assert_eq!(&derived[2080 - 30..2080], &source[2050 - 30..2050]);
    assert_eq!(&derived[2100..], &source[2050..]);
    let loci = effective_loci(&with_insertion(&format!("{second}{INSERT_AT_1000}")));
    assert!(
        loci.contains(&("STR-GATA".to_owned(), 2081, 2100)),
        "{loci:?}"
    );
    assert!(
        loci.contains(&("STR-INS".to_owned(), 1001, 1030)),
        "{loci:?}"
    );
}

/// **El borde del punto de inserción.** Un tracto existente que **acaba** en la
/// base ancla no se mueve; uno que **empieza** justo después se desplaza.
///
/// Es la prueba que distingue «tractos insertados tras una base *anterior*» de
/// «tras una base anterior *o igual*»: con la segunda, el final del primer
/// locus saltaría treinta bases y lo partiría.
#[test]
fn a_tract_ending_at_the_anchor_stays_and_one_starting_after_it_moves() {
    let before = r#"
[[reference.repeat_loci]]
id = "STR-BEFORE"
sequence = "chr1"
start = 2001
end = 2040
motif = "CA"
[[reference.repeat_loci.alleles]]
units = 20
frequency = 1.0
"#;
    let inserting_at_2040 = r#"
[[reference.repeat_loci]]
id = "STR-INS"
sequence = "chr1"
insert_after = 2040
motif = "CAG"
reference_units = [10]
[[reference.repeat_loci.alleles]]
units = 10
frequency = 1.0
"#;
    let loci = effective_loci(&with_insertion(&format!("{before}{inserting_at_2040}")));
    assert!(
        loci.contains(&("STR-BEFORE".to_owned(), 2001, 2040)),
        "{loci:?}"
    );
    assert!(
        loci.contains(&("STR-INS".to_owned(), 2041, 2070)),
        "{loci:?}"
    );

    let inserting_at_2000 = inserting_at_2040.replace("insert_after = 2040", "insert_after = 2000");
    let loci = effective_loci(&with_insertion(&format!("{before}{inserting_at_2000}")));
    assert!(
        loci.contains(&("STR-BEFORE".to_owned(), 2031, 2070)),
        "{loci:?}"
    );
    assert!(
        loci.contains(&("STR-INS".to_owned(), 2001, 2030)),
        "{loci:?}"
    );
}

/// **Las regiones del generador de la matriz se traducen.** La región
/// 1001..1010 de partida, justo después del punto de inserción, pasa a
/// 1031..1040 del genoma final, y las tres variantes caen ahí.
///
/// El fixture está **saturado** a propósito: sin traducir, la región caería
/// encima del tracto insertado (1001..1030 del final), y la regla de que nada
/// toca un microsatélite rechazaría el escenario entero. No hay forma de que
/// una región sin traducir pase esta prueba.
#[test]
fn the_generation_regions_are_translated_to_the_final_genome() {
    let body = with_insertion(INSERT_AT_1000)
        .replace("start = 3000\nend = 3500", "start = 1001\nend = 1010");
    assert!(
        body.contains("start = 1001\nend = 1010"),
        "el fixture cambió la región"
    );
    let scenario = publish(&body);
    let matrix = String::from_utf8(scenario.read("population_snvs.tsv")).expect("TSV");
    let derived = scenario.unfolded("resources/references/ref-main-str.fa");
    let rows = matrix.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(rows.len(), 3, "{matrix}");
    for line in rows {
        let columns = line.split('\t').collect::<Vec<_>>();
        let position = columns[2].parse::<usize>().expect("posición");
        assert!((1031..=1040).contains(&position), "{line}");
        assert_eq!(
            &derived[position - 1..position],
            columns[3],
            "la REF de la matriz es la base del genoma final"
        );
    }
}

/// **Los sorteos se siembran con el identificador de la referencia final**
/// (DTC-098). El oráculo vuelve a sortear cada copia con la función pública
/// del sorteo y el identificador final, y cuenta las copias de cada alelo; el
/// manifiesto tiene que publicar esas cuentas.
///
/// El fixture se comprueba discriminante dentro de la prueba: con el
/// identificador de la suministrada las cuentas salen distintas. Si alguna
/// vez coincidieran, la prueba no podría distinguir los dos ámbitos y lo dice
/// fallando en vez de pasar en falso.
#[test]
fn the_repeat_draws_are_seeded_with_the_final_reference() {
    use genome_population_simulator::application::EffectiveInitialGenotypeModel;
    use genome_population_simulator::identity::{IndividualId, PopulationId};
    use genome_population_simulator::repeat_genotypes::draw_repeat_allele_index;

    let ids = (1..=20)
        .map(|index| format!("\"I{index:02}\""))
        .collect::<Vec<_>>();
    let population = format!(
        "\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 20\nindividual_ids = [{}]\nploidy = 2\ngenotype_model = \"hwe\"\n",
        ids.join(", ")
    );
    let body = document(
        "inserted_reference_id = \"ref-main-str\"",
        &format!("{INSERT_AT_1000}{EXISTING_CA}"),
        &population,
    );
    let (root, config) = write_config(&body);
    let effective = load_and_validate_config(&config).expect("configuración válida");
    let locus = effective
        .repeat_loci()
        .iter()
        .find(|locus| locus.id() == "STR-CA")
        .expect("locus");
    let counts = |reference_id: &str| {
        let mut counts = vec![0_u64; locus.alleles().len()];
        for index in 1..=20 {
            for copy in 0..2_u8 {
                let drawn = draw_repeat_allele_index(
                    effective.reproducibility().root_seed(),
                    reference_id,
                    &PopulationId::parse("P1", "test").expect("población"),
                    &IndividualId::parse(format!("I{index:02}"), "test").expect("individuo"),
                    EffectiveInitialGenotypeModel::HardyWeinberg,
                    copy,
                    locus,
                )
                .expect("sorteo");
                counts[drawn] += 1;
            }
        }
        counts
    };
    let expected = counts("ref-main-str");
    assert_ne!(
        expected,
        counts("ref-main"),
        "el fixture no distingue los dos ámbitos de azar"
    );
    let destination = root.join("published");
    execute_validated_config(&effective, &destination).expect("publicación");
    let manifest: Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifiesto"))
            .expect("JSON");
    let published = manifest["repeat_loci"]["loci"]
        .as_array()
        .expect("loci")
        .iter()
        .find(|candidate| candidate["id"] == "STR-CA")
        .expect("locus")["alleles"]
        .as_array()
        .expect("alelos")
        .iter()
        .map(|allele| allele["realized_copies"].as_u64().expect("copias"))
        .collect::<Vec<_>>();
    assert_eq!(published, expected);
}

// ---------------------------------------------------------------------------
// Identidad y publicación
// ---------------------------------------------------------------------------

/// La referencia final tiene identidad propia y todo lo publicado la nombra a
/// ella: la cabecera de la verdad con la huella de su FASTA, el `##contig` y el
/// `##reference` del VCF, la matriz, la sección de loci del manifiesto. El
/// manifiesto la registra derivada de la suministrada.
#[test]
fn every_published_output_names_the_final_reference() {
    let scenario = publish(&with_insertion(&format!("{INSERT_AT_1000}{EXISTING_CA}")));
    let manifest = scenario.json("manifest.json");
    let effective = scenario.json("effective-config.json");
    assert_eq!(
        manifest["references"],
        serde_json::json!([
            {"id": "ref-main", "fasta_artifact_id": "reference-fasta"},
            {"id": "ref-main-str", "fasta_artifact_id": "reference-fasta-ref-main-str"}
        ])
    );
    let insertion = &manifest["reference_insertion"];
    assert_eq!(insertion["id"], "ref-main-str");
    assert_eq!(insertion["source_reference_id"], "ref-main");
    assert_eq!(insertion["source_fasta_artifact_id"], "reference-fasta");
    assert_eq!(
        insertion["insertions"],
        serde_json::json!([{"locus_id": "STR-INS", "sequence_id": "chr1", "insert_after": 1000, "inserted_bases": 30}])
    );
    assert_eq!(effective["reference"]["id"], "ref-main");
    assert_eq!(effective["inserted_reference"]["id"], "ref-main-str");
    assert_eq!(
        effective["inserted_reference"]["sequences"][0]["length"],
        4030
    );
    assert_eq!(manifest["repeat_loci"]["reference_id"], "ref-main-str");

    let truth_header: Value = serde_json::from_str(
        String::from_utf8(scenario.read("outputs/truth.jsonl"))
            .expect("verdad")
            .lines()
            .next()
            .expect("cabecera"),
    )
    .expect("cabecera JSON");
    assert_eq!(truth_header["reference_id"], "ref-main-str");
    let final_fasta = scenario
        .destination
        .join("resources/references/ref-main-str.fa");
    let reference_digest = truth_header["source_digests"]
        .as_array()
        .expect("huellas")
        .iter()
        .find(|digest| digest["source"]["kind"] == "reference")
        .expect("huella de la referencia");
    assert_eq!(reference_digest["source"]["source_id"], "ref-main-str");
    assert_eq!(reference_digest["sha256"], external_sha256(&final_fasta));

    let vcf = Command::new("bgzip")
        .arg("-dc")
        .arg(scenario.destination.join("outputs/truth.vcf.gz"))
        .output()
        .expect("bgzip");
    let vcf = String::from_utf8(vcf.stdout).expect("VCF");
    assert!(vcf.contains("##contig=<ID=chr1,length=4030>"), "{vcf}");
    assert!(vcf.contains("##reference=ref-main-str"), "{vcf}");

    let metadata = scenario.json("population_snvs.metadata.json");
    assert_eq!(metadata["reference"]["reference_id"], "ref-main-str");
    assert_eq!(
        metadata["reference_fasta"]["sha256"],
        external_sha256(&final_fasta)
    );
}

/// Los tres documentos nuevos o completados validan contra sus esquemas con el
/// validador externo: la configuración efectiva v6, el manifiesto general v6 y
/// el mapa de coordenadas v2. El mapa además describe el empalme: tres
/// segmentos, el del medio insertado y sin intervalo de origen.
#[test]
fn the_published_documents_validate_against_their_schemas() {
    let scenario = publish(&with_insertion(INSERT_AT_1000));
    assert!(validates(
        &scenario.json("effective-config.json"),
        "effective-config-v6.schema.json"
    ));
    assert!(validates(
        &scenario.json("manifest.json"),
        "manifest-general-v6.schema.json"
    ));
    let map = scenario.json("resources/references/ref-main-str.map.json");
    assert!(validates(&map, "reference-coordinate-map-v2.schema.json"));
    assert_eq!(map["source_reference_id"], "ref-main");
    assert_eq!(map["derived_reference_id"], "ref-main-str");
    assert_eq!(
        map["segments"],
        serde_json::json!([
            {"materialization_order": 1, "segment_class": "retained", "sequence_id": "chr1",
             "sequence_class": "chromosome", "orientation": "forward",
             "source_interval": {"start": 1, "end": 1000}, "derived_interval": {"start": 1, "end": 1000}},
            {"materialization_order": 2, "segment_class": "inserted", "sequence_id": "chr1",
             "sequence_class": "chromosome", "orientation": "forward",
             "derived_interval": {"start": 1001, "end": 1030}},
            {"materialization_order": 3, "segment_class": "retained", "sequence_id": "chr1",
             "sequence_class": "chromosome", "orientation": "forward",
             "source_interval": {"start": 1001, "end": 4000}, "derived_interval": {"start": 1031, "end": 4030}}
        ])
    );
}

/// Con el generador de alelos por deslizamiento el alelo ancestral es el tracto
/// insertado, leído de la referencia final: `(CAG)10`. Si el generador leyera
/// la suministrada, en 1001..1030 encontraría bases al azar y no un tracto.
/// La segunda mitad publica con pasos de verdad para que el juez del abanico,
/// que lo rehace con la semilla y el identificador de la referencia final,
/// tenga algo que comparar.
#[test]
fn the_slippage_generator_starts_from_the_inserted_tract() {
    let generated = INSERT_AT_1000
        .split("[[reference.repeat_loci.alleles]]")
        .next()
        .expect("cabecera del locus")
        .to_owned()
        + "[reference.repeat_loci.slippage]\nmodel = \"stepwise\"\nmean_events = 0.0\n";
    let scenario = publish(&with_insertion(&generated));
    let effective = scenario.json("effective-config.json");
    let alleles = effective["repeat_loci"][0]["alleles"]
        .as_array()
        .expect("alelos");
    assert_eq!(alleles.len(), 1, "media cero: solo el ancestral");
    assert_eq!(alleles[0]["sequence"], TRACT_30);

    // Con pasos de verdad, el juez rehace el abanico sembrado con la
    // referencia final; un productor que sembrara con la suministrada
    // publicaría otro abanico y no pasaría de aquí.
    let stepping = generated.replace("mean_events = 0.0", "mean_events = 3.0");
    let scenario = publish(&with_insertion(&stepping));
    let effective = scenario.json("effective-config.json");
    assert!(
        effective["repeat_loci"][0]["slippage"]["applied_events"]
            .as_u64()
            .expect("pasos")
            > 0
    );
}

/// Con lecturas, el manifiesto del evento técnico enlaza el FASTA **final**, y
/// todas las lecturas salen de él: el tracto insertado aparece en ellas.
#[test]
fn the_technical_event_reads_the_final_reference() {
    let event = r#"
[technical_event]
event_id = "reads-main"
individual_ids = ["I1"]
target_coverage_wgs = 20
sequencing_errors = "disabled"
[[technical_event.libraries]]
individual_id = "I1"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#;
    let scenario = publish(&document(
        "inserted_reference_id = \"ref-main-str\"",
        INSERT_AT_1000,
        &format!("{POPULATION_AND_VCF}{event}"),
    ));
    let manifest = scenario.json("manifest.json");
    let event_manifest_path = manifest["artifacts"]
        .as_array()
        .expect("artefactos")
        .iter()
        .find(|artifact| artifact["kind"] == "event_manifest")
        .and_then(|artifact| artifact["path"].as_str())
        .expect("manifiesto del evento")
        .to_owned();
    let event_manifest = scenario.json(&event_manifest_path);
    let rendered = event_manifest.to_string();
    assert!(rendered.contains("ref-main-str"), "{rendered}");
    let fastq_bytes = manifest["artifacts"]
        .as_array()
        .expect("artefactos")
        .iter()
        .filter(|artifact| artifact["kind"] == "fastq")
        .map(|artifact| {
            let path = scenario
                .destination
                .join(artifact["path"].as_str().expect("ruta"));
            let output = Command::new("bgzip")
                .arg("-dc")
                .arg(&path)
                .output()
                .expect("bgzip");
            if output.status.success() {
                output.stdout
            } else {
                fs::read(&path).expect("FASTQ")
            }
        })
        .collect::<Vec<_>>()
        .concat();
    let fastq = String::from_utf8_lossy(&fastq_bytes);
    // Un trozo del flanco izquierdo pegado al principio del tracto: solo existe
    // en el genoma final.
    let source = source_bases();
    let junction = format!("{}{}", &source[985..1000], &TRACT_30[..15]);
    let reverse = junction
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect::<String>();
    assert!(
        fastq.contains(&junction) || fastq.contains(&reverse),
        "ninguna lectura cruza la unión del tracto insertado"
    );
}

/// **Sin inserción no cambia nada.** El mismo escenario sin `insert_after`
/// publica configuración efectiva y manifiesto sin ninguna de las secciones
/// nuevas, y el locus sin `declared`.
///
/// La identidad byte a byte frente a la base anterior se comprueba contra un
/// binario compilado aparte desde ella; dentro de la batería solo puede
/// comprobarse la ausencia de lo nuevo.
#[test]
fn without_insertion_nothing_new_is_published() {
    let scenario = publish(&document("", EXISTING_CA, POPULATION_AND_VCF));
    let effective = scenario.json("effective-config.json");
    let manifest = scenario.json("manifest.json");
    assert!(effective.get("inserted_reference").is_none());
    assert!(manifest.get("reference_insertion").is_none());
    assert!(effective["repeat_loci"][0].get("declared").is_none());
    assert_eq!(
        manifest["references"]
            .as_array()
            .expect("referencias")
            .len(),
        1
    );
    assert!(
        !scenario
            .destination
            .join("resources/references/ref-main-str.fa")
            .exists()
    );
}

// ---------------------------------------------------------------------------
// Dos secuencias: una inserción no mueve nada de la otra
// ---------------------------------------------------------------------------

/// Segunda secuencia del escenario de dos cromosomas: 1000 bases
/// pseudoaleatorias con un `(CA)10` en 501..520.
fn second_sequence_bases() -> String {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut noisy = |length: usize| -> String {
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                char::from(b"ACGT"[(state % 4) as usize])
            })
            .collect()
    };
    let bases = format!("{}{}{}", noisy(500), "CA".repeat(10), noisy(480));
    assert_eq!(bases.len(), 1000);
    bases
}

/// Escenario de dos secuencias: la región del generador está en `chr2`
/// 801..810 y `loci` decide qué se inserta y dónde.
fn two_sequence_document(loci: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 3
spatial_mode = "weighted_regions"
outside_weight = 0.0
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.regions]]
contig = "chr2"
start = 801
end = 810
weight = 1.0

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"
inserted_reference_id = "ref-main-str"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 4000

[[reference.sequences]]
id = "chr2"
class = "chromosome"
length = 1000
{loci}
{POPULATION_AND_VCF}"#
    )
}

fn publish_two_sequences(body: &str) -> Scenario {
    let root = temporary("two");
    let fasta = format!(
        ">chr1\n{}\n>chr2\n{}\n",
        source_bases(),
        second_sequence_bases()
    );
    fs::write(root.join("reference.fa"), fasta).expect("FASTA de dos secuencias");
    let config = root.join("general.toml");
    fs::write(&config, body).expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Scenario { root, destination }
}

/// Las secuencias del FASTA publicado, por identificador.
fn fasta_sequences(text: &str) -> Vec<(String, String)> {
    let mut sequences: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            sequences.push((
                name.split_whitespace().next().unwrap_or("").to_owned(),
                String::new(),
            ));
        } else if let Some((_, bases)) = sequences.last_mut() {
            bases.push_str(line);
        }
    }
    sequences
}

const CA_ON_CHR2: &str = r#"
[[reference.repeat_loci]]
id = "STR-CA2"
sequence = "chr2"
start = 501
end = 520
motif = "CA"
[[reference.repeat_loci.alleles]]
units = 10
frequency = 0.5
[[reference.repeat_loci.alleles]]
units = 11
frequency = 0.5
"#;

/// **Una inserción en `chr1` no mueve nada de `chr2`.** Treinta bases tras la
/// base 100 de `chr1`; en `chr2`, el locus 501..520 y la región 801..810 del
/// generador siguen exactamente donde estaban, y la secuencia de `chr2` del
/// FASTA final es la de partida, byte a byte.
///
/// Es la única prueba que discrimina el filtro por secuencia de la traducción:
/// sin él, `chr2` heredaría el desplazamiento de `chr1`, el locus caería fuera
/// de su tracto y la región, en 831..840.
#[test]
fn an_insertion_in_one_sequence_moves_nothing_in_another() {
    let insert_chr1 = INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 100");
    let scenario = publish_two_sequences(&two_sequence_document(&format!(
        "{insert_chr1}{CA_ON_CHR2}"
    )));
    let published = fasta_sequences(
        &String::from_utf8(scenario.read("resources/references/ref-main-str.fa")).expect("FASTA"),
    );
    assert_eq!(published.len(), 2);
    assert_eq!(published[0].1.len(), 4030);
    assert_eq!(external_offsets(&published[0].1, TRACT_30), vec![100]);
    assert_eq!(published[1].0, "chr2");
    assert_eq!(published[1].1, second_sequence_bases());

    let effective = scenario.json("effective-config.json");
    let ca2 = effective["repeat_loci"]
        .as_array()
        .expect("loci")
        .iter()
        .find(|locus| locus["id"] == "STR-CA2")
        .expect("locus de chr2");
    assert_eq!(
        (ca2["start"].clone(), ca2["end"].clone()),
        (Value::from(501), Value::from(520))
    );
    assert_eq!(
        ca2["declared"],
        serde_json::json!({"start": 501, "end": 520})
    );

    let matrix = String::from_utf8(scenario.read("population_snvs.tsv")).expect("TSV");
    let rows = matrix.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(rows.len(), 3, "{matrix}");
    for line in rows {
        let columns = line.split('\t').collect::<Vec<_>>();
        assert_eq!(columns[1], "chr2", "{line}");
        let position = columns[2].parse::<u64>().expect("posición");
        assert!((801..=810).contains(&position), "{line}");
    }
}

/// **El mismo punto en dos secuencias distintas no es el mismo punto.** Un
/// tracto tras la base 100 de `chr1` y otro tras la 100 de `chr2` se admiten,
/// y cada uno entra en su secuencia.
#[test]
fn the_same_point_in_two_sequences_is_admitted() {
    let insert_chr1 = INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 100");
    let insert_chr2 = insert_chr1
        .replace("STR-INS", "STR-INS-2")
        .replace("sequence = \"chr1\"", "sequence = \"chr2\"");
    let scenario = publish_two_sequences(&two_sequence_document(&format!(
        "{insert_chr1}{insert_chr2}"
    )));
    let published = fasta_sequences(
        &String::from_utf8(scenario.read("resources/references/ref-main-str.fa")).expect("FASTA"),
    );
    assert_eq!(external_offsets(&published[0].1, TRACT_30), vec![100]);
    assert_eq!(external_offsets(&published[1].1, TRACT_30), vec![100]);
    assert_eq!(published[1].1.len(), 1030);
}

// ---------------------------------------------------------------------------
// Rechazos de la configuración
// ---------------------------------------------------------------------------

/// Cada rechazo de la primera versión, con su control positivo al lado: la
/// configuración base se admite, así que cada rechazo se debe a la pieza que
/// cambia.
#[test]
fn each_unsupported_declaration_is_rejected_naming_what_collides() {
    let base = with_insertion(INSERT_AT_1000);
    let (root, config) = write_config(&base);
    load_and_validate_config(&config).expect("control positivo");
    fs::remove_dir_all(root).ok();

    let cases: Vec<(&str, String, &str)> = vec![
        (
            "insert_after junto a start",
            with_insertion(
                &INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 1000\nstart = 1001"),
            ),
            "son excluyentes",
        ),
        (
            "sin reference_units",
            with_insertion(&INSERT_AT_1000.replace("reference_units = [10]\n", "")),
            "no declara reference_units",
        ),
        (
            "unidades que no son una por bloque",
            with_insertion(
                &INSERT_AT_1000.replace("reference_units = [10]", "reference_units = [10, 2]"),
            ),
            "hace falta uno por bloque",
        ),
        (
            "un bloque sin unidades",
            with_insertion(
                &INSERT_AT_1000.replace("reference_units = [10]", "reference_units = [0]"),
            ),
            "cero unidades",
        ),
        (
            "un tracto que rebasa el tope",
            with_insertion(
                &INSERT_AT_1000.replace("reference_units = [10]", "reference_units = [201]"),
            ),
            "rebasa el tope de 200",
        ),
        (
            "antes de la primera base",
            with_insertion(&INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 0")),
            "no deja una base a cada lado",
        ),
        (
            "después de la última",
            with_insertion(&INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 4000")),
            "no deja una base a cada lado",
        ),
        (
            "sin inserted_reference_id",
            document("", INSERT_AT_1000, POPULATION_AND_VCF),
            "no declara inserted_reference_id",
        ),
        (
            "inserted_reference_id sin inserción",
            with_insertion(EXISTING_CA),
            "sin ningún locus que inserte",
        ),
        (
            "inserted_reference_id igual al de la suministrada",
            document(
                "inserted_reference_id = \"ref-main\"",
                INSERT_AT_1000,
                POPULATION_AND_VCF,
            ),
            "repite el identificador de la referencia suministrada",
        ),
        (
            "dos en el mismo punto",
            with_insertion(&format!(
                "{INSERT_AT_1000}{}",
                INSERT_AT_1000.replace("STR-INS", "STR-INS-2")
            )),
            "en el mismo punto que el locus STR-INS",
        ),
        (
            "dentro del tracto de otro locus",
            with_insertion(&format!(
                "{EXISTING_CA}{}",
                INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 2020")
            )),
            "lo partiría",
        ),
        (
            "una región del generador partida",
            with_insertion(&INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 3200")),
            "un tracto insertado la partiría",
        ),
        // Insertar tras la **primera** base de un tramo también lo parte: el
        // tracto quedaría entre la primera y la segunda.
        (
            "tras la primera base del tramo de otro locus",
            with_insertion(&format!(
                "{EXISTING_CA}{}",
                INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 2001")
            )),
            "lo partiría",
        ),
        (
            "tras la primera base de una región del generador",
            with_insertion(&INSERT_AT_1000.replace("insert_after = 1000", "insert_after = 3000")),
            "un tracto insertado la partiría",
        ),
        (
            "junto a una reducción",
            with_insertion(INSERT_AT_1000)
                + "\n[[reduced_references]]\nid = \"ref-small\"\nsource_reference_id = \"ref-main\"\nregion_size_bases = 100\ntarget_bases = 1000\n",
            "junto a una reducción de referencia",
        ),
        (
            "junto a cambios directos",
            with_insertion(INSERT_AT_1000)
                + "\n[[direct_events]]\nevent_id = \"E1\"\nindividual_id = \"I1\"\nphase = 0\nkind = \"snp\"\ncontig = \"chr1\"\nposition_1based = 10\nref = \"A\"\nalt = \"C\"\n",
            "junto a cambios directos",
        ),
    ];
    for (label, body, expected) in cases {
        let rendered = rejection(&body);
        assert!(
            rendered.contains(expected),
            "{label}: se esperaba «{expected}» en {rendered}"
        );
    }
}

/// La matriz de SNV suministrada o reutilizada se rechaza junto a una
/// inserción: está atada al genoma de partida.
#[test]
fn a_supplied_or_reused_matrix_is_rejected_with_an_insertion() {
    for (mode, key) in [
        ("supplied", "input_path = \"matrix.tsv\""),
        ("reuse", "metadata_path = \"population_snvs.metadata.json\""),
    ] {
        let body = with_insertion(INSERT_AT_1000).replace(
            MATRIX_GENERATE,
            &format!("\n[population_snvs]\nmode = \"{mode}\"\n{key}\n"),
        );
        let rendered = rejection(&body);
        assert!(
            rendered.contains(&format!(
                "matriz de SNV {}",
                if mode == "supplied" {
                    "suministrada"
                } else {
                    "reutilizada"
                }
            )),
            "{mode}: {rendered}"
        );
    }
}

/// Con referencia simulada no se inserta: el tracto se coloca en su tramo.
#[test]
fn an_insertion_on_a_simulated_reference_is_rejected() {
    let body = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[reference]
kind = "simulated"
id = "ref-sim"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 4000
{INSERT_AT_1000}"#
    );
    assert!(rejection(&body).contains("solo se admite con referencia suministrada"));
    let with_id = body.replace(
        "id = \"ref-sim\"",
        "id = \"ref-sim\"\ninserted_reference_id = \"x\"",
    );
    assert!(
        rejection(&with_id).contains("una referencia simulada no admite inserted_reference_id")
    );
}

// ---------------------------------------------------------------------------
// El esquema de la configuración general
// ---------------------------------------------------------------------------

/// El esquema general v4 admite el locus que inserta y rechaza cada forma mal
/// escrita. Con control positivo para que «rechazar todo» no pase.
#[test]
fn the_general_schema_accepts_an_insertion_and_rejects_its_malformed_forms() {
    let document = |reference_extra: Value, locus: Value| {
        let mut reference = serde_json::json!({
            "kind": "supplied",
            "id": "ref-main",
            "fasta": "reference.fa",
            "sequences": [{"id": "chr1", "class": "chromosome", "length": 4000}],
            "repeat_loci": [locus]
        });
        for (key, value) in reference_extra.as_object().expect("objeto") {
            reference[key] = value.clone();
        }
        serde_json::json!({
            "config_kind": "general",
            "schema_version": "4",
            "reference": reference
        })
    };
    let allele = serde_json::json!([{"units": 10, "frequency": 1.0}]);
    let inserting = serde_json::json!({
        "id": "STR-INS", "sequence": "chr1", "insert_after": 1000,
        "motif": "CAG", "reference_units": [10], "alleles": allele
    });
    let named = serde_json::json!({"inserted_reference_id": "ref-main-str"});
    let schema = "config-general-v4.schema.json";
    assert!(validates(
        &document(named.clone(), inserting.clone()),
        schema
    ));

    let mut with_start = inserting.clone();
    with_start["start"] = serde_json::json!(1001);
    let mut without_units = inserting.clone();
    without_units
        .as_object_mut()
        .expect("objeto")
        .remove("reference_units");
    let mut zero = inserting.clone();
    zero["insert_after"] = serde_json::json!(0);
    let existing = serde_json::json!({
        "id": "STR-CA", "sequence": "chr1", "start": 2001, "end": 2040,
        "motif": "CA", "alleles": allele
    });
    let mut existing_with_units = existing.clone();
    existing_with_units["reference_units"] = serde_json::json!([20]);
    for (label, candidate) in [
        (
            "insert_after con start",
            document(named.clone(), with_start),
        ),
        (
            "sin reference_units",
            document(named.clone(), without_units),
        ),
        ("insert_after cero", document(named.clone(), zero)),
        (
            "sin inserted_reference_id",
            document(serde_json::json!({}), inserting.clone()),
        ),
        (
            "inserted_reference_id sin inserción",
            document(named.clone(), existing.clone()),
        ),
        (
            "reference_units en tracto existente",
            document(serde_json::json!({}), existing_with_units),
        ),
    ] {
        assert!(
            !validates(&candidate, schema),
            "{label} tendría que rechazarse"
        );
    }
    // Y el locus de tracto existente sigue admitiéndose sin el nombre final.
    assert!(validates(
        &document(serde_json::json!({}), existing),
        schema
    ));
}

// ---------------------------------------------------------------------------
// El juez independiente
// ---------------------------------------------------------------------------

/// Reescribe un artefacto publicado y actualiza su huella y su tamaño en el
/// manifiesto, de modo que la corrupción llegue al juez semántico en vez de
/// caer en la comprobación de integridad.
fn republish(destination: &Path, relative: &str, bytes: &[u8]) {
    let path = destination.join(relative);
    fs::write(&path, bytes).expect("artefacto reescrito");
    let identity = identify_file(&path, relative).expect("huella");
    let manifest_path = destination.join("manifest.json");
    let mut manifest: Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("manifiesto")).expect("JSON");
    for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
        if artifact["path"] == relative {
            artifact["checksum"]["digest"] = Value::String(identity.digest().to_owned());
            artifact["size_bytes"] = Value::from(identity.size_bytes());
        }
    }
    fs::write(
        &manifest_path,
        serde_json::to_vec_pretty(&manifest).expect("manifiesto"),
    )
    .expect("manifiesto reescrito");
}

fn judge(destination: &Path) -> Result<(), String> {
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    validate_published_resources(
        destination,
        &effective_path,
        &fs::read(&effective_path).expect("efectiva"),
        &manifest_path,
        &fs::read(&manifest_path).expect("manifiesto"),
        None,
    )
    .map_err(|error| format!("{error:?}"))
}

/// Copia la publicación entera a un directorio nuevo para corromperla sin
/// tocar la original.
fn copy_tree(from: &Path, to: &Path) {
    let status = Command::new("cp")
        .arg("-r")
        .arg(from)
        .arg(to)
        .status()
        .expect("cp");
    assert!(status.success());
}

/// **Cada corrupción de lo publicado la rechaza su propia rama del juez.** Las
/// corrupciones mantienen la integridad de los ficheros —se rehace la huella
/// en el manifiesto—, de modo que lo que las rechaza es el contenido.
///
/// Control positivo: la publicación intacta pasa el juez.
#[test]
fn every_corruption_of_an_inserted_publication_is_rejected_by_its_own_branch() {
    let scenario = publish(&with_insertion(&format!("{INSERT_AT_1000}{EXISTING_CA}")));
    judge(&scenario.destination).expect("control positivo");

    type Corruption = fn(&mut Value, &mut Value, &mut String, &mut Value);
    let cases: Vec<(&str, Corruption, &str)> = vec![
        (
            "una base del FASTA final",
            |_, _, fasta, _| {
                let index = fasta.find('\n').expect("cabecera") + 1 + 5;
                let original = fasta.as_bytes()[index];
                let replacement = if original == b'A' { "C" } else { "A" };
                fasta.replace_range(index..=index, replacement);
            },
            "no coincide base a base",
        ),
        (
            "el segmento insertado del mapa una base más allá",
            |_, _, _, map| {
                map["segments"][1]["derived_interval"] =
                    serde_json::json!({"start": 1002, "end": 1031});
            },
            "no describe el empalme",
        ),
        (
            "el intervalo publicado de un locus",
            |effective, manifest, _, _| {
                effective["repeat_loci"][1]["start"] = serde_json::json!(2030);
                manifest["repeat_loci"]["loci"][1]["start"] = serde_json::json!(2030);
            },
            "no es la traducción al genoma final",
        ),
        (
            "las unidades del tracto insertado",
            |effective, _, _, _| {
                effective["repeat_loci"][0]["reference_units"] = serde_json::json!([11]);
            },
            "los cambios publicados no son los que declaran los loci",
        ),
        (
            "un locus que inserta sin unidades",
            |effective, _, _, _| {
                effective["repeat_loci"][0]
                    .as_object_mut()
                    .expect("locus")
                    .remove("reference_units");
            },
            "sin motivos y unidades por bloque",
        ),
        (
            "un tracto insertado por encima de su tope",
            |effective, _, _, _| {
                effective["repeat_loci"][0]["reference_units"] = serde_json::json!([201]);
            },
            "o con más unidades que su tope",
        ),
        // DT-025: el tope `max_units` sale del mismo documento que las unidades y
        // no acota nada; la cota que muerde es la física, que lo insertado
        // quepa en el FASTA final. El primer caso, de tres millones de unidades,
        // **termina** sin esa cota (reserva unos megabytes) y rechaza por otro
        // motivo, de modo que la prueba muere por el motivo y no cuelga la
        // batería; va delante del segundo, que sin la cota pediría 8 GiB.
        (
            "tres millones de unidades con un tope que las admite",
            |effective, _, _, _| {
                effective["repeat_loci"][0]["reference_units"] = serde_json::json!([3_000_000]);
                effective["repeat_loci"][0]["max_units"] = serde_json::json!(3_000_000);
            },
            "o con más unidades que su tope",
        ),
        (
            "u64::MAX / 3 unidades con un tope que las admite",
            |effective, _, _, _| {
                effective["repeat_loci"][0]["reference_units"] =
                    serde_json::json!([6_148_914_691_236_517_205_u64]);
                effective["repeat_loci"][0]["max_units"] =
                    serde_json::json!(6_148_914_691_236_517_205_u64);
            },
            "o con más unidades que su tope",
        ),
        (
            "el intervalo de un locus solo en el manifiesto",
            |_, manifest, _, _| {
                manifest["repeat_loci"]["loci"][1]["start"] = serde_json::json!(2032);
            },
            "no es la traducción al genoma final",
        ),
        (
            "el punto de un cambio solo en el manifiesto",
            |_, manifest, _, _| {
                manifest["reference_insertion"]["insertions"][0]["insert_after"] =
                    serde_json::json!(999);
            },
            "no coincide con la referencia declarada o con reference_insertion",
        ),
        (
            "un cambio publicado que partiría la región del generador",
            |effective, manifest, _, _| {
                effective["inserted_reference"]["insertions"][0]["insert_after"] =
                    serde_json::json!(3100);
                manifest["reference_insertion"]["insertions"][0]["insert_after"] =
                    serde_json::json!(3100);
            },
            "no son la traducción al genoma final de las del TOML",
        ),
        (
            "la sección de loci del manifiesto con la referencia suministrada",
            |_, manifest, _, _| {
                manifest["repeat_loci"]["reference_id"] = serde_json::json!("ref-main");
            },
            "la sección de loci del manifiesto no nombra la referencia de la simulación",
        ),
        (
            "lo declarado de un locus",
            |effective, _, _, _| {
                effective["repeat_loci"][1]["declared"] =
                    serde_json::json!({"start": 2002, "end": 2040});
            },
            "lo declarado de cada locus no coinciden con el TOML",
        ),
        (
            "la sección de la efectiva ausente",
            |effective, _, _, _| {
                effective
                    .as_object_mut()
                    .expect("efectiva")
                    .remove("inserted_reference");
            },
            "no coinciden en si hay referencia final",
        ),
        (
            "la sección de la efectiva sin forma de objeto",
            |effective, _, _, _| {
                effective["inserted_reference"] = serde_json::json!("ref-main-str");
            },
            "inserted_reference no es objeto",
        ),
        (
            "el mapa de la efectiva distinto del manifiesto",
            |effective, _, _, _| {
                effective["inserted_reference"]["map_artifact_id"] =
                    serde_json::json!("reference-fasta");
            },
            "no coincide con la referencia declarada o con reference_insertion",
        ),
        (
            "la sección del manifiesto sin forma de objeto",
            |_, manifest, _, _| {
                manifest["reference_insertion"] = serde_json::json!(1);
            },
            "reference_insertion no es objeto",
        ),
        (
            "la referencia final igual a la de origen",
            |_, manifest, _, _| {
                manifest["reference_insertion"]["id"] = serde_json::json!("ref-main");
            },
            "distinta de su origen",
        ),
        (
            "la versión del mapa",
            |_, manifest, _, _| {
                manifest["reference_insertion"]["map_schema_version"] = serde_json::json!("1");
            },
            "la clase y la versión del mapa",
        ),
        (
            "el estado de la inserción",
            |_, manifest, _, _| {
                manifest["reference_insertion"]["state"] = serde_json::json!("done");
            },
            "no tiene un estado válido",
        ),
        (
            "una inserción fallida en una ejecución completa",
            |_, manifest, _, _| {
                manifest["reference_insertion"]["state"] = serde_json::json!("failed");
            },
            "una ejecución complete exige",
        ),
        (
            "la lista de cambios vacía",
            |_, manifest, _, _| {
                manifest["reference_insertion"]["insertions"] = serde_json::json!([]);
            },
            "no publica ningún cambio",
        ),
        (
            "un cambio sin sus cuatro campos",
            |_, manifest, _, _| {
                manifest["reference_insertion"]["insertions"][0]
                    .as_object_mut()
                    .expect("cambio")
                    .remove("inserted_bases");
            },
            "no tiene exactamente locus, secuencia, punto y bases",
        ),
    ];
    for (label, corrupt, expected) in cases {
        let copy = temporary("corrupt").join("published");
        copy_tree(&scenario.destination, &copy);
        let mut effective = scenario.json("effective-config.json");
        let mut manifest = scenario.json("manifest.json");
        let mut fasta = String::from_utf8(scenario.read("resources/references/ref-main-str.fa"))
            .expect("FASTA");
        let mut map = scenario.json("resources/references/ref-main-str.map.json");
        corrupt(&mut effective, &mut manifest, &mut fasta, &mut map);
        republish(
            &copy,
            "resources/references/ref-main-str.fa",
            fasta.as_bytes(),
        );
        republish(
            &copy,
            "resources/references/ref-main-str.map.json",
            &serde_json::to_vec_pretty(&map).expect("mapa"),
        );
        // El manifiesto se reescribe con los cambios del caso **y** con las
        // huellas que acaban de rehacerse.
        let refreshed: Value =
            serde_json::from_slice(&fs::read(copy.join("manifest.json")).expect("manifiesto"))
                .expect("JSON");
        manifest["artifacts"] = refreshed["artifacts"].clone();
        fs::write(
            copy.join("manifest.json"),
            serde_json::to_vec_pretty(&manifest).expect("manifiesto"),
        )
        .expect("manifiesto");
        republish(
            &copy,
            "effective-config.json",
            &serde_json::to_vec_pretty(&effective).expect("efectiva"),
        );
        let rendered = judge(&copy).expect_err(label);
        assert!(
            rendered.contains(expected),
            "{label}: se esperaba «{expected}» en {rendered}"
        );
    }
}

/// **La cota física de lo insertado es acumulada entre loci** (DT-025): dos
/// tractos que caben cada uno en el FASTA final pero no los dos a la vez. Cada
/// tracto es el 60 % del tamaño del fichero; sin acumular, los dos pasan la cota
/// y el rechazo llega más tarde, por otro motivo (los cambios publicados no son
/// los declarados).
///
/// Control positivo: la publicación de dos inserciones, intacta, pasa el juez.
/// Debe fallar contra: una cota que no acumule. Se ejecutó.
#[test]
fn the_physical_bound_on_inserted_bases_adds_up_across_loci() {
    let second = INSERT_AT_1000
        .replace("STR-INS", "STR-INS-2")
        .replace("insert_after = 1000", "insert_after = 500");
    let scenario = publish(&with_insertion(&format!(
        "{INSERT_AT_1000}{second}{EXISTING_CA}"
    )));
    judge(&scenario.destination).expect("control positivo");
    let mut effective = scenario.json("effective-config.json");
    let size = scenario.read("resources/references/ref-main-str.fa").len();
    let units = (size * 6 / 10 / 3) as u64;
    let inserting = effective["repeat_loci"]
        .as_array_mut()
        .expect("loci")
        .iter_mut()
        .filter(|locus| locus["declared"].get("insert_after").is_some())
        .collect::<Vec<_>>();
    assert_eq!(inserting.len(), 2, "el fixture inserta en dos puntos");
    for locus in inserting {
        locus["reference_units"] = serde_json::json!([units]);
        locus["max_units"] = serde_json::json!(units);
    }
    let copy = temporary("corrupt").join("published");
    copy_tree(&scenario.destination, &copy);
    republish(
        &copy,
        "effective-config.json",
        &serde_json::to_vec_pretty(&effective).expect("efectiva"),
    );
    let rendered = judge(&copy).expect_err("dos tractos que no caben juntos");
    assert!(
        rendered.contains("o con más unidades que su tope"),
        "{rendered}"
    );
}

/// La autorización de matriz vacía también funciona con inserción: una
/// población sin una sola variante, sobre la referencia final.
#[test]
fn an_empty_generated_matrix_is_admitted_with_an_insertion() {
    let body = with_insertion(INSERT_AT_1000).replace("locus_count = 3", "locus_count = 0");
    let (root, config) = write_config(&body);
    let destination = root.join("published");
    execute_validated_config_with_authorizations(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
        ExecutionAuthorizations {
            allow_empty_snv_matrix: true,
            allow_incompatible_genotypes: false,
        },
    )
    .expect("publicación con matriz vacía");
    assert!(
        destination
            .join("resources/references/ref-main-str.fa")
            .exists()
    );
}

/// Descomprime el VCF BGZF publicado, aplica `edit` al texto y lo vuelve a
/// comprimir con `bgzip`, de modo que lo que cambia es el contenido y no el
/// formato.
fn rewrite_vcf(destination: &Path, edit: impl Fn(&str) -> String) {
    let path = destination.join("outputs/truth.vcf.gz");
    let plain = Command::new("bgzip")
        .arg("-dc")
        .arg(&path)
        .output()
        .expect("bgzip -dc");
    let edited = edit(&String::from_utf8(plain.stdout).expect("VCF"));
    let staging = destination.join("outputs/edited.vcf");
    fs::write(&staging, edited).expect("VCF editado");
    let compressed = Command::new("bgzip")
        .arg("-c")
        .arg(&staging)
        .output()
        .expect("bgzip -c");
    fs::remove_file(&staging).ok();
    republish(destination, "outputs/truth.vcf.gz", &compressed.stdout);
}

/// Reescribe la verdad aplicando `edit` a cada registro JSON.
fn rewrite_truth(destination: &Path, edit: impl Fn(usize, &mut Value)) {
    let relative = "outputs/truth.jsonl";
    let text = fs::read_to_string(destination.join(relative)).expect("verdad");
    let lines = text
        .lines()
        .enumerate()
        .map(|(index, line)| {
            let mut record: Value = serde_json::from_str(line).expect("registro");
            edit(index, &mut record);
            serde_json::to_string(&record).expect("registro")
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    republish(destination, relative, lines.as_bytes());
}

/// Cambia el primer `sequence_length` que encuentra en un registro.
fn first_sequence_length(value: &mut Value) -> Option<&mut Value> {
    match value {
        Value::Object(object) => {
            if object.contains_key("sequence_length") {
                return object.get_mut("sequence_length");
            }
            object.values_mut().find_map(first_sequence_length)
        }
        Value::Array(items) => items.iter_mut().find_map(first_sequence_length),
        _ => None,
    }
}

/// **La verdad y el VCF tienen que nombrar la referencia final** (DTC-098).
/// Cada corrupción hace que una salida nombre la suministrada o mida mal una
/// secuencia; la integridad del fichero se rehace en el manifiesto, de modo
/// que lo que la rechaza es el contenido, cada una en su propia rama.
///
/// Control positivo: la publicación intacta pasa el juez.
#[test]
fn every_output_that_names_the_wrong_reference_is_rejected() {
    let scenario = publish(&with_insertion(&format!("{INSERT_AT_1000}{EXISTING_CA}")));
    judge(&scenario.destination).expect("control positivo");
    type Corruption = fn(&Path);
    let cases: Vec<(&str, Corruption, &str)> = vec![
        (
            "la cabecera de la verdad con la suministrada",
            |destination| {
                rewrite_truth(destination, |index, record| {
                    if index == 0 {
                        record["reference_id"] = serde_json::json!("ref-main");
                    }
                });
            },
            "la cabecera de la verdad no nombra la referencia de la simulación",
        ),
        (
            "la huella del contenido de otra referencia",
            |destination| {
                rewrite_truth(destination, |index, record| {
                    if index == 0 {
                        record["reference_content_sha256"] = serde_json::json!("0".repeat(64));
                    }
                });
            },
            "la huella del contenido en la cabecera de la verdad",
        ),
        (
            "la huella del FASTA de la suministrada",
            |destination| {
                rewrite_truth(destination, |index, record| {
                    if index == 0 {
                        for digest in record["source_digests"].as_array_mut().expect("huellas") {
                            if digest["source"]["kind"] == "reference" {
                                digest["sha256"] = serde_json::json!("1".repeat(64));
                            }
                        }
                    }
                });
            },
            "la huella de referencia de la verdad no es la del FASTA",
        ),
        (
            "un fragmento con la longitud de la suministrada",
            |destination| {
                rewrite_truth(destination, |index, record| {
                    if index > 0
                        && let Some(length) = first_sequence_length(record)
                    {
                        *length = serde_json::json!(4000);
                    }
                });
            },
            "un fragmento de la verdad no está en",
        ),
        (
            "el ##reference del VCF",
            |destination| {
                rewrite_vcf(destination, |text| {
                    text.replace("##reference=ref-main-str", "##reference=ref-main")
                });
            },
            "la cabecera del VCF no declara ##reference=",
        ),
        (
            "el ##contig del VCF",
            |destination| {
                rewrite_vcf(destination, |text| {
                    text.replace("length=4030", "length=4000")
                });
            },
            "los ##contig del VCF no son el catálogo",
        ),
    ];
    for (label, corrupt, expected) in cases {
        let copy = temporary("names").join("published");
        copy_tree(&scenario.destination, &copy);
        corrupt(&copy);
        let rendered = judge(&copy).expect_err(label);
        assert!(
            rendered.contains(expected),
            "{label}: se esperaba «{expected}» en {rendered}"
        );
    }
}

/// **La referencia final no puede ser su propio origen.** El manifiesto se
/// corrompe para que la final repita el identificador y el FASTA de la
/// suministrada, de modo que las demás condiciones de esa comprobación se
/// cumplen y solo `id == source` puede rechazarlo. Se llama solo al lector
/// del manifiesto, que es donde vive esa comprobación.
#[test]
fn a_final_reference_equal_to_its_source_is_rejected_by_the_manifest_reader() {
    use genome_population_simulator::results::validate_general_manifest_bytes;
    let scenario = publish(&with_insertion(INSERT_AT_1000));
    let path = scenario.destination.join("manifest.json");
    let intact = scenario.read("manifest.json");
    validate_general_manifest_bytes(&path, &intact).expect("control positivo");
    let mut manifest = scenario.json("manifest.json");
    manifest["reference_insertion"]["id"] = serde_json::json!("ref-main");
    manifest["reference_insertion"]["fasta_artifact_id"] = serde_json::json!("reference-fasta");
    let corrupted = serde_json::to_vec_pretty(&manifest).expect("manifiesto");
    let rendered = format!(
        "{:?}",
        validate_general_manifest_bytes(&path, &corrupted).expect_err("final igual al origen")
    );
    assert!(rendered.contains("distinta de su origen"), "{rendered}");
}
