//! F6.4.3.3: la verdad publicada con alelos de repetición, versión **3**
//! hasta el 30/09/2026 y **5** desde entonces (la 5 es la 3 sin la copia
//! entera del haplotipo en los alelos estructurales).
//!
//! # Cómo están elegidos los fixtures
//!
//! Lo que estas pruebas afirman es que la verdad lleva **la secuencia** de cada
//! alelo y no un contador, y que la versión sube **por lo publicado** y no por
//! lo declarado. Para que puedan fallar:
//!
//! - la pareja de alelos que demuestra lo primero tiene **la misma longitud** y
//!   la interrupción en sitio distinto: con longitudes distintas, una
//!   proyección que publicara solo el número de unidades seguiría
//!   distinguiéndolos y la prueba no diría nada;
//! - la comparación de la versión se hace contra un escenario que **declara el
//!   locus** y cuyo único alelo es el de la referencia. Si se comparara con uno
//!   que no declarase nada, una promoción disparada por la declaración
//!   sobreviviría.

use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_fasta::{SequenceExpectation, validate_supplied_fasta};
use genome_population_simulator::truth_validation::validate_and_reconstruct_truth_jsonl;
use serde_json::Value;

#[path = "support/detection_section.rs"]
mod detection_section;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
/// Trescientas bases: cien de flanco, cien de tracto `(CA)50` y cien más.
const SEQUENCE_LENGTH: u64 = 300;

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f6-4-3-3-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

fn fasta_bases() -> String {
    let mut bases = "G".repeat(100);
    bases.push_str(&"CA".repeat(50));
    bases.push_str(&"G".repeat(100));
    assert_eq!(bases.len() as u64, SEQUENCE_LENGTH);
    bases
}

/// Escenario con referencia suministrada, una matriz mínima y dos individuos
/// diploides. `loci` se pega tal cual.
fn scenario(loci: &str) -> (PathBuf, PathBuf) {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    fs::write(
        root.join("reference.fa"),
        format!(">chr1\n{}\n", fasta_bases()),
    )
    .expect("FASTA");
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {SEQUENCE_LENGTH}
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"
"#
        )),
    )
    .expect("configuración");
    (root, config)
}

fn publish(loci: &str) -> (PathBuf, PathBuf) {
    let (root, config) = scenario(loci);
    let effective = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective, &destination).expect("publicación");
    (root, destination)
}

fn truth_lines(destination: &Path) -> Vec<Value> {
    fs::read_to_string(destination.join("outputs/truth.jsonl"))
        .expect("verdad publicada")
        .lines()
        .map(|line| serde_json::from_str(line).expect("línea JSON"))
        .collect()
}

/// Un locus polimórfico: dos alelos de **distinta** longitud.
const POLYMORPHIC: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 45
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 40
frequency = 0.5
"#;

/// Un locus cuyo **único alelo es el de la referencia**: cincuenta copias de
/// `CA`, que es exactamente lo que hay entre las bases 101 y 200.
const REFERENCE_ONLY: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 50
frequency = 1.0
"#;

/// Dos alelos de **igual longitud** con la interrupción en sitio distinto,
/// declarados por su secuencia exacta.
fn equal_length_pair() -> String {
    let left = format!("CATT{}", "CA".repeat(48));
    let right = format!("{}TTCA", "CA".repeat(48));
    assert_eq!(left.len(), 100);
    assert_eq!(right.len(), 100);
    format!(
        r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
min_purity = 0.9

[[reference.repeat_loci.alleles]]
sequence = "{left}"
frequency = 0.5

[[reference.repeat_loci.alleles]]
sequence = "{right}"
frequency = 0.5
"#
    )
}

/// **El registro de un diploide con alelos de repetición lleva cada uno con su
/// fase, su locus, su motivo y su número de unidades, y su secuencia.**
#[test]
fn the_published_truth_carries_the_locus_the_motif_the_units_and_the_sequence() {
    let (root, destination) = publish(POLYMORPHIC);
    let lines = truth_lines(&destination);
    let header = &lines[0];
    assert_eq!(header["record_type"], "header");
    assert_eq!(header["schema_version"], "5");

    let pieces: Vec<&Value> = lines
        .iter()
        .filter(|line| line["record_type"] == "haplotype_copy")
        .flat_map(|copy| copy["current_pieces"].as_array().expect("piezas"))
        .filter(|piece| piece["kind"] == "repeat_allele")
        .collect();
    assert!(
        !pieces.is_empty(),
        "el escenario tiene que publicar al menos un alelo de repetición"
    );
    for piece in &pieces {
        assert_eq!(piece["locus_id"], "STR-1");
        assert_eq!(piece["source"]["start"], 100);
        assert_eq!(piece["source"]["end"], 200);
        let bases = piece["bases"].as_str().expect("bases");
        let blocks = piece["blocks"].as_array().expect("bloques");
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0]["motif"], "CA");
        assert_eq!(blocks[0]["motif_length"], 2);
        // La correspondencia que este esquema no puede comprobar: las bases del
        // bloque son las del alelo y las unidades son su mitad.
        assert_eq!(blocks[0]["bases"], bases.len());
        assert_eq!(blocks[0]["units"], (bases.len() / 2) as f64);
        assert!(bases == "CA".repeat(45) || bases == "CA".repeat(40));
    }

    let cases: Vec<&Value> = lines
        .iter()
        .filter(|line| line["case_kind"] == "repeat")
        .collect();
    assert!(!cases.is_empty(), "hay al menos un caso de repetición");
    for case in cases {
        assert_eq!(case["locus"]["start"], 100);
        assert_eq!(case["locus"]["end"], 200);
        let alleles = case["alleles"].as_array().expect("alelos");
        assert_eq!(alleles[0]["kind"], "reference");
        for allele in &alleles[1..] {
            assert_eq!(allele["kind"], "repeat");
            assert_eq!(allele["locus_id"], "STR-1");
            assert!(allele["bases"].as_str().expect("bases").len() >= 80);
        }
        // Las dosis suman la ploidía y los portadores son las dos copias, con
        // su fase: es lo que hace faseado al registro.
        let dosages: u64 = case["allele_dosages"]
            .as_array()
            .expect("dosis")
            .iter()
            .map(|value| value.as_u64().expect("entero"))
            .sum();
        assert_eq!(dosages, 2);
        let carriers = case["carriers"].as_array().expect("portadores");
        assert_eq!(carriers.len(), 2);
        assert_eq!(carriers[0]["copy"]["copy_index"], 0);
        assert_eq!(carriers[1]["copy"]["copy_index"], 1);
    }
    fs::remove_dir_all(&root).ok();
}

/// **Dos alelos de igual longitud y distinta interrupción producen registros
/// distinguibles.**
///
/// *Qué hace que pueda fallar:* los dos miden cien bases y tienen el mismo
/// número de unidades, cincuenta. Una proyección que publicara solo el número
/// de unidades y no la secuencia —la variante defectuosa que el paso nombra—
/// haría los dos registros idénticos.
#[test]
fn two_alleles_of_equal_length_with_the_interruption_elsewhere_are_distinguishable() {
    let (root, destination) = publish(&equal_length_pair());
    let lines = truth_lines(&destination);
    let published: std::collections::BTreeSet<String> = lines
        .iter()
        .filter(|line| line["record_type"] == "haplotype_copy")
        .flat_map(|copy| copy["current_pieces"].as_array().expect("piezas"))
        .filter(|piece| piece["kind"] == "repeat_allele")
        .map(|piece| piece["bases"].as_str().expect("bases").to_owned())
        .collect();
    assert_eq!(
        published.len(),
        2,
        "las dos secuencias publicadas son distintas aunque midan lo mismo: {published:?}"
    );
    let lengths: std::collections::BTreeSet<usize> = published.iter().map(String::len).collect();
    assert_eq!(
        lengths.len(),
        1,
        "y miden lo mismo, que es lo que hace discriminante la comparación"
    );
    fs::remove_dir_all(&root).ok();
}

/// **Un escenario con un locus declarado cuyo único alelo es el de la
/// referencia publica la verdad en su versión anterior y byte a byte igual a la
/// de un escenario sin loci.**
///
/// *Qué hace que pueda fallar:* la comparación es contra un escenario que **sí
/// declara** el locus. Una promoción disparada por la declaración en vez de por
/// la publicación de un alelo cambiaría la cabecera y las dos aserciones
/// morirían.
///
/// *Sobre el número de versión:* el plan escribió «versión 2», pero la verdad
/// que publica este recorrido —el de la fase de población inicial— es la
/// **versión 4** (la 1 hasta el 30/09/2026); la 2 es la de la fase de población fundadora, que tiene su
/// propio publicador. Lo que el paso pide es que no se promueva, y eso es lo
/// que se comprueba.
///
/// *Y sobre el «byte a byte»:* tampoco es literal, y no puede serlo. El
/// identificador de ejecución se deriva de los bytes del fichero de
/// configuración, de modo que declarar el locus lo cambia necesariamente. Lo
/// que se compara es el resto del documento con ese identificador
/// normalizado, y que los dos identificadores sean **distintos** se comprueba
/// también, para que la normalización no tape una comparación vacía.
#[test]
fn a_locus_whose_only_allele_is_the_reference_one_does_not_promote_the_truth() {
    let (bare_root, bare) = publish("");
    let (locus_root, with_locus) = publish(REFERENCE_ONLY);
    let bare_bytes = fs::read_to_string(bare.join("outputs/truth.jsonl")).expect("verdad sin loci");
    let locus_bytes =
        fs::read_to_string(with_locus.join("outputs/truth.jsonl")).expect("verdad con locus");
    assert_eq!(
        truth_lines(&with_locus)[0]["schema_version"],
        "4",
        "sin alelo publicado no hay promoción"
    );
    // El identificador de ejecución se deriva de los bytes del fichero de
    // configuración, así que declarar el locus lo cambia por fuerza. Lo demás
    // —cabecera, copias, genotipos— tiene que ser idéntico.
    let normalize = |text: &str, id: &str| text.replace(id, "EXECUTION");
    let bare_id = truth_lines(&bare)[0]["execution_id"]
        .as_str()
        .expect("identificador")
        .to_owned();
    let locus_id = truth_lines(&with_locus)[0]["execution_id"]
        .as_str()
        .expect("identificador")
        .to_owned();
    assert_ne!(
        bare_id, locus_id,
        "los dos ficheros de entrada son distintos"
    );
    assert_eq!(
        normalize(&bare_bytes, &bare_id),
        normalize(&locus_bytes, &locus_id),
        "declarar un locus monomórfico no cambia nada más de la verdad"
    );
    fs::remove_dir_all(&bare_root).ok();
    fs::remove_dir_all(&locus_root).ok();
}

/// El **validador independiente** reconstruye la verdad con alelos de
/// repetición desde la referencia publicada.
#[test]
fn the_independent_validator_replays_a_truth_with_repeat_alleles() {
    let (root, destination) = publish(POLYMORPHIC);
    let reference = validate_supplied_fasta(
        &destination.join("resources/references/ref-main.fa"),
        "ref-main",
        &[SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: SEQUENCE_LENGTH,
        }],
    )
    .expect("referencia publicada");
    let file = fs::File::open(destination.join("outputs/truth.jsonl")).expect("verdad");
    validate_and_reconstruct_truth_jsonl(&mut BufReader::new(file), reference.reference())
        .expect("la verdad sigue siendo reconstruible de forma independiente");
    fs::remove_dir_all(&root).ok();
}

/// Cada línea publicada valida contra `truth-v5.schema.json` con un validador
/// **externo**, y la misma línea de cabecera **no** valida contra el esquema de
/// la versión 4 ni contra el histórico de la 3.
#[test]
fn every_published_line_validates_against_the_external_v5_schema() {
    let (root, destination) = publish(POLYMORPHIC);
    let lines = truth_lines(&destination);
    assert_external(&lines, "truth-v5.schema.json", true);
    // Control negativo: la cabecera declara la versión 5, que el esquema sin
    // repeticiones y el histórico rechazan. Sin él, un esquema que aceptara
    // cualquier cosa pasaría.
    assert_external(&lines[..1], "truth-v4.schema.json", false);
    assert_external(&lines[..1], "truth-v3.schema.json", false);
    fs::remove_dir_all(&root).ok();
}

/// **Guarda de lo que la proyección al VCF ya hace, que no es todo lo que
/// F6.5.1.1 pedirá.**
///
/// Hoy el alelo llega al VCF **siempre explícito**, con la secuencia completa
/// en `REF` y `ALT`, alineado a la izquierda, en **un solo registro
/// multialélico** y con el genotipo faseado. Lo que **todavía no hay** son los
/// campos de repetición en tándem del estándar —secuencia del motivo, su
/// longitud, número de copias, bases y número de secuencias repetidas—, que
/// son el trabajo de F6.5.1.1. Esta prueba fija lo que ya funciona para que el
/// paso siguiente no lo rompa sin enterarse; **no** da por hecho el paso.
#[test]
fn the_repeat_allele_already_reaches_the_vcf_as_an_explicit_multiallelic_record() {
    let (root, config) = scenario(POLYMORPHIC);
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            "{}\n[vcf_output]\nmode = \"text\"\n",
            fs::read_to_string(&config).expect("configuración")
        )),
    )
    .expect("configuración con VCF");
    let effective = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective, &destination).expect("publicación");
    let vcf = fs::read_to_string(destination.join("outputs/truth.vcf")).expect("VCF publicado");

    let records: Vec<&str> = vcf
        .lines()
        .filter(|line| line.starts_with("chr1\t"))
        .collect();
    // Ningún alelo simbólico, ni siquiera con un tracto de cien bases. Se mira
    // solo en los registros: la cabecera lleva `<` en cada declaración de
    // campo, y buscarlo en el fichero entero haría pasar la búsqueda siempre.
    assert!(
        records.iter().all(|line| !line.contains('<')),
        "la forma explícita aplica siempre: {records:?}"
    );
    assert!(!records.is_empty(), "control positivo: hay registros");
    let tract_records: Vec<&&str> = records
        .iter()
        .filter(|line| {
            line.split('\t')
                .nth(3)
                .is_some_and(|reference| reference.len() > 1)
        })
        .collect();
    assert_eq!(
        tract_records.len(),
        1,
        "un solo registro para el locus, no uno por alelo: {records:?}"
    );
    let fields: Vec<&str> = tract_records[0].split('\t').collect();
    let reference_allele = fields[3];
    let alternates: Vec<&str> = fields[4].split(',').collect();
    assert_eq!(alternates.len(), 2, "los dos alelos del locus");
    // El ancla de la izquierda más el tracto: el `REF` publicado mide lo que la
    // diferencia con el alelo más corto exige, y los dos `ALT` son más cortos.
    for alternate in &alternates {
        assert!(alternate.len() < reference_allele.len());
        assert!(alternate.starts_with(&reference_allele[..1]));
    }
    // Genotipo faseado en las dos muestras.
    for sample in &fields[9..] {
        assert!(sample.contains('|'), "genotipo faseado: {sample}");
    }
    fs::remove_dir_all(&root).ok();
}

/// **Los esquemas de verdad heredados no cambian ni un byte.**
#[test]
fn the_inherited_truth_schemas_are_byte_for_byte_unchanged() {
    let expected = [
        (
            "truth-v1.schema.json",
            "e3d49abd11447370d34e4acf29b755a288fb310641474218086ed60099233fd7",
        ),
        (
            "truth-v2.schema.json",
            "be53a6b32eb04689bcb4113731bb97f84e57f3566116e98fa2a2500313ba647b",
        ),
        (
            "truth-v3.schema.json",
            "09dcd19b6e764e83c72549fd6c8c0b3b9ff485fc249547a97a3bfc9af31b12bf",
        ),
    ];
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for (name, digest) in expected {
        let identity = genome_population_simulator::integrity::identify_file(
            &root.join("schemas").join(name),
            name,
        )
        .expect("huella del esquema");
        assert_eq!(identity.digest(), digest, "el esquema {name} ha cambiado");
    }
}

fn assert_external(lines: &[Value], schema: &str, expected: bool) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let bundle = temporary("schema");
    fs::create_dir_all(&bundle).expect("paquete aislado");
    fs::copy(root.join("schemas").join(schema), bundle.join(schema)).expect("copia del esquema");
    let cases: Vec<Value> = lines
        .iter()
        .map(|instance| serde_json::json!({"instance": instance, "expected": expected}))
        .collect();
    let case_path = bundle.join("cases.json");
    fs::write(&case_path, serde_json::to_vec(&cases).expect("casos")).expect("escritura");
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(bundle.join(schema))
        .arg(&case_path)
        .arg("--cases")
        .output()
        .expect("validador externo");
    fs::remove_dir_all(&bundle).ok();
    assert!(
        output.status.success(),
        "{schema}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Locus **compuesto y saturado**: tres bloques cuyos motivos, longitudes de
/// motivo, copias y bases son **todos distintos entre sí**.
///
/// La saturación es lo que hace discriminante la comprobación de las listas
/// paralelas. Con bloques de igual longitud de motivo o igual número de copias
/// —como los del microsatélite compuesto real `(TCTA)n(TCTG)n(TCTA)n`—
/// publicar el motivo del primer bloque en todos, o desalinear las listas, no
/// se notaría. Aquí: motivos `CAG`, `TA` y `GGGG`; longitudes 3, 2 y 4; copias
/// 7, 3 y 2; bases 21, 6 y 8.
const COMPOSITE_TRACT: &str = "CAGCAGCAGCAGCAGCAGCAGTATATAGGGGGGGG";

fn composite_locus() -> String {
    assert_eq!(COMPOSITE_TRACT.len(), 35);
    r#"
[[reference.repeat_loci]]
id = "STR-C"
sequence = "chr1"
start = 101
end = 135
motifs = ["CAG", "TA", "GGGG"]

[[reference.repeat_loci.alleles]]
units = [7, 3, 2]
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = [5, 3, 2]
frequency = 0.5
"#
    .to_owned()
}

/// Escenario con la referencia que lleva el tracto compuesto en su sitio.
fn publish_composite() -> (PathBuf, PathBuf) {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    // Flancos de `G` y `T` al azar, no tiras de `G`: entre dos tiras largas de
    // `G` el detector no ve el `(CAG)×7`, y el locus compuesto no coincidiría
    // con nada (ver `gt_background`).
    let bases = detection_section::gt_background(0x9E37_79B9_7F4A_7C15, 100)
        + COMPOSITE_TRACT
        + &detection_section::gt_background(0xD1B5_4A32_D192_ED03, 100);
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {}
{}
[[initial_populations]]
population_id = "P1"
individual_count = 4
individual_ids = ["I1", "I2", "I3", "I4"]
ploidy = 2
genotype_model = "hwe"

[vcf_output]
mode = "text"
"#,
            bases.len(),
            composite_locus()
        )),
    )
    .expect("configuración");
    let effective = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective, &destination).expect("publicación");
    (root, destination)
}

/// **La verdad publica un locus compuesto con sus tres bloques alineados.**
///
/// Es la reincidencia que la ronda anterior encontró en la configuración
/// efectiva: sin un fixture compuesto, publicar el motivo del primer bloque en
/// todos ellos sobrevivía a la batería entera.
#[test]
fn a_composite_locus_publishes_its_three_blocks_without_mixing_them() {
    let (root, destination) = publish_composite();
    let lines = truth_lines(&destination);
    let pieces: Vec<&Value> = lines
        .iter()
        .filter(|line| line["record_type"] == "haplotype_copy")
        .flat_map(|copy| copy["current_pieces"].as_array().expect("piezas"))
        .filter(|piece| piece["kind"] == "repeat_allele")
        .collect();
    assert!(!pieces.is_empty(), "el escenario publica algún alelo");
    for piece in &pieces {
        let blocks = piece["blocks"].as_array().expect("bloques");
        assert_eq!(blocks.len(), 3, "tres bloques");
        // Los motivos, en su orden y sin repetir el primero.
        let motifs: Vec<&str> = blocks
            .iter()
            .map(|block| block["motif"].as_str().expect("motivo"))
            .collect();
        assert_eq!(motifs, vec!["CAG", "TA", "GGGG"]);
        // Y cada lista paralela con el valor de SU bloque: las longitudes, las
        // bases y las copias son todas distintas entre sí, de modo que
        // cualquier desalineación cambia algún número.
        let lengths: Vec<u64> = blocks
            .iter()
            .map(|block| block["motif_length"].as_u64().expect("longitud"))
            .collect();
        assert_eq!(lengths, vec![3, 2, 4]);
        let bases: Vec<u64> = blocks
            .iter()
            .map(|block| block["bases"].as_u64().expect("bases"))
            .collect();
        let units: Vec<f64> = blocks
            .iter()
            .map(|block| block["units"].as_f64().expect("unidades"))
            .collect();
        // El primer bloque es el que cambia entre los dos alelos declarados
        // —siete copias o cinco—; los otros dos son los mismos en ambos, y por
        // eso sus valores están fijos y delatan cualquier cruce de listas.
        assert!(bases[0] == 21 || bases[0] == 15, "{bases:?}");
        assert_eq!(&bases[1..], &[6, 8]);
        assert_eq!(units[0], bases[0] as f64 / 3.0);
        assert_eq!(&units[1..], &[3.0, 2.0]);
        // El segundo bloque tiene 3 copias de un motivo de 2 y el tercero 2 de
        // uno de 4: si las dos listas se cruzaran, estos dos productos dejarían
        // de coincidir con sus bases.
        for index in 0..3 {
            assert_eq!(
                bases[index] as f64,
                units[index] * lengths[index] as f64,
                "el bloque {index} no cuadra"
            );
        }
        // Y la suma es el tracto entero, no el de otro alelo.
        assert_eq!(
            bases.iter().sum::<u64>() as usize,
            piece["bases"].as_str().expect("secuencia").len()
        );
    }
    fs::remove_dir_all(&root).ok();
}

/// **El caso heterocigoto referencia/alternativo existe y se publica bien.**
///
/// Es el genotipo más común en la realidad y ningún fixture lo producía: todos
/// eran homocigotos para el alternativo, de modo que la rama del lector
/// independiente que resuelve «esta copia conserva el tracto de la referencia»
/// no la ejercitaba nada.
#[test]
fn a_heterozygous_reference_over_alternative_case_is_published_and_replayed() {
    let (root, destination) = publish_composite();
    let cases: Vec<Value> = truth_lines(&destination)
        .into_iter()
        .filter(|line| line["case_kind"] == "repeat")
        .collect();
    let heterozygous: Vec<&Value> = cases
        .iter()
        .filter(|case| {
            case["allele_dosages"]
                .as_array()
                .expect("dosis")
                .iter()
                .map(|value| value.as_u64().expect("entero"))
                .collect::<Vec<_>>()
                == vec![1, 1]
        })
        .collect();
    assert!(
        !heterozygous.is_empty(),
        "el corpus tiene que producir al menos un heterocigoto referencia/alternativo: {cases:?}"
    );
    for case in heterozygous {
        let carriers = case["carriers"].as_array().expect("portadores");
        let indices: Vec<u64> = carriers
            .iter()
            .map(|carrier| carrier["allele_index"].as_u64().expect("índice"))
            .collect();
        assert_eq!(indices.len(), 2);
        assert!(indices.contains(&0), "una copia lleva el de la referencia");
        assert!(indices.contains(&1), "y la otra el alternativo");
    }
    // Y el lector independiente lo replica: es el camino que ningún fixture
    // recorría.
    let reference = validate_supplied_fasta(
        &destination.join("resources/references/ref-main.fa"),
        "ref-main",
        &[SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 235,
        }],
    )
    .expect("referencia publicada");
    let file = fs::File::open(destination.join("outputs/truth.jsonl")).expect("verdad");
    validate_and_reconstruct_truth_jsonl(&mut BufReader::new(file), reference.reference())
        .expect("la verdad heterocigota se reconstruye");
    fs::remove_dir_all(&root).ok();
}

/// **El oráculo semántico del VCF examina una verdad con repeticiones.**
///
/// Es código nuevo que hasta ahora no ejecutaba ninguna prueba: el oráculo no
/// corre en el recorrido de publicación, solo desde las pruebas, y ningún
/// documento con repeticiones llegaba a él.
#[test]
fn the_independent_vcf_oracle_accepts_a_truth_with_repeat_alleles() {
    let (root, destination) = publish_composite();
    let reference = validate_supplied_fasta(
        &destination.join("resources/references/ref-main.fa"),
        "ref-main",
        &[SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 235,
        }],
    )
    .expect("referencia publicada");
    let truth =
        genome_population_simulator::truth_validation::read_truth_jsonl(&mut BufReader::new(
            fs::File::open(destination.join("outputs/truth.jsonl")).expect("verdad"),
        ))
        .expect("verdad legible");
    assert!(
        truth.carries_repeat_alleles(),
        "control positivo: el documento lleva alelos de repetición"
    );
    genome_population_simulator::vcf_semantic_validation::validate_vcf_against_truth(
        reference.reference(),
        &truth,
        &mut BufReader::new(fs::File::open(destination.join("outputs/truth.vcf")).expect("VCF")),
    )
    .expect("el oráculo semántico independiente acepta el VCF");
    fs::remove_dir_all(&root).ok();
}

/// **Un alelo que acaba en unidad parcial se publica con su número decimal.**
///
/// Ningún fixture producía uno, de modo que la comparación con tolerancia del
/// lector independiente —la que comprueba que el decimal publicado es de
/// verdad las bases entre la longitud del motivo— nunca veía un valor no
/// entero y se ejercitaba sola contra números redondos.
#[test]
fn an_allele_ending_in_a_partial_unit_publishes_a_fractional_count() {
    let partial = format!("CA{}C", "CA".repeat(48));
    assert_eq!(partial.len(), 99, "noventa y nueve bases, motivo de dos");
    let loci = format!(
        r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 50
frequency = 0.5

[[reference.repeat_loci.alleles]]
sequence = "{partial}"
frequency = 0.5
"#
    );
    let (root, destination) = publish(&loci);
    let lines = truth_lines(&destination);
    let fractional: Vec<&Value> = lines
        .iter()
        .filter(|line| line["record_type"] == "haplotype_copy")
        .flat_map(|copy| copy["current_pieces"].as_array().expect("piezas"))
        .filter(|piece| piece["kind"] == "repeat_allele")
        .collect();
    assert!(
        !fractional.is_empty(),
        "el alelo parcial llega a alguna copia"
    );
    for piece in &fractional {
        let block = &piece["blocks"].as_array().expect("bloques")[0];
        assert_eq!(block["bases"], 99);
        assert_eq!(block["units"], 49.5, "las unidades son decimales");
        assert_eq!(piece["bases"].as_str().expect("secuencia").len(), 99);
    }
    // Y el lector independiente lo acepta: es su comparación con tolerancia la
    // que aquí se ejercita por primera vez con un valor no entero.
    let reference = validate_supplied_fasta(
        &destination.join("resources/references/ref-main.fa"),
        "ref-main",
        &[SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: SEQUENCE_LENGTH,
        }],
    )
    .expect("referencia publicada");
    let file = fs::File::open(destination.join("outputs/truth.jsonl")).expect("verdad");
    validate_and_reconstruct_truth_jsonl(&mut BufReader::new(file), reference.reference())
        .expect("la unidad parcial se reconstruye");
    fs::remove_dir_all(&root).ok();
}
