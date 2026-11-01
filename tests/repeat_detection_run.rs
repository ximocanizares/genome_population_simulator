//! Entrega 3 del plan de detección de microsatélites: **la detección en el
//! producto** (pasos F6.4.a.10, F6.4.a.10b, F6.4.a.11 y F6.4.a.12).
//!
//! # Los oráculos no llaman al detector (riesgo R-1)
//!
//! Cada prueba planta los tractos en posiciones **escritas aquí** sobre un fondo
//! sin repeticiones de periodo 1 ni 2 (cada base es distinta de las dos
//! anteriores), y espera exactamente esos loci, con su identificador escrito a
//! mano: `SSR-<índice de secuencia>-<inicio>-<fin>-<motivo>` en base uno. Lo que
//! se compara es lo que la ejecución **publica** en la configuración efectiva,
//! no una estructura interna.
//!
//! # Lo que estas pruebas no ven, dicho sin adornos
//!
//! La exhaustividad del detector sobre un fondo cualquiera no se prueba aquí:
//! el fondo está construido para no tener tractos, y un tracto absorbido en un
//! segmento impuro no se planta (riesgo R-7). Lo cubren las pruebas del propio
//! detector (`src/repeat_detection.rs`).

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::{
    ExecutionAuthorizations, execute_validated_config,
    preflight_validated_config_with_authorizations,
};
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED_A: &str = "00000000000000000000000000000000000000000000000000000000000000d1";
const SEED_B: &str = "7a0000000000000000000000000000000000000000000000000000000000ff02";

/// Deslizamiento común de todas las pruebas.
const SLIPPAGE: &str =
    "[reference.repeat_detection.slippage]\nmodel = \"stepwise\"\nmean_events = 1.0\n\n";

fn temporary(label: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "gps-detection-run-{label}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&path).expect("directorio temporal");
    path
}

/// Fondo sin repeticiones de periodo 1 ni 2: cada base es distinta de las dos
/// anteriores, elegida al azar (xorshift64) entre las que quedan.
fn background(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    let mut bases: Vec<u8> = Vec::with_capacity(length);
    while bases.len() < length {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let candidate = b"ACGT"[(state % 4) as usize];
        let n = bases.len();
        if (n >= 1 && bases[n - 1] == candidate) || (n >= 2 && bases[n - 2] == candidate) {
            continue;
        }
        bases.push(candidate);
    }
    bases
}

/// Planta `motif × units` desde la base `start` (base uno) y pone en cada
/// flanco una base que no prolonga el tracto: distinta de la primera y de la
/// última del motivo.
fn plant(bases: &mut [u8], start: usize, motif: &str, units: usize) -> (usize, usize) {
    let tract = motif.repeat(units);
    let begin = start - 1;
    let end = begin + tract.len();
    bases[begin..end].copy_from_slice(tract.as_bytes());
    let motif = motif.as_bytes();
    let flank = *b"ACGT"
        .iter()
        .find(|base| **base != motif[0] && **base != motif[motif.len() - 1])
        .expect("una base libre");
    bases[begin - 1] = flank;
    bases[end] = flank;
    (start, end)
}

fn write_fasta(path: &Path, sequences: &[(&str, Vec<u8>)]) {
    let mut text = String::new();
    for (id, bases) in sequences {
        text.push('>');
        text.push_str(id);
        text.push('\n');
        for chunk in bases.chunks(60) {
            text.push_str(std::str::from_utf8(chunk).expect("ASCII"));
            text.push('\n');
        }
    }
    fs::write(path, text).expect("FASTA");
}

/// Una fila de la matriz por cada `(secuencia, posición)`, con la base de
/// referencia leída de las secuencias.
fn write_matrix(path: &Path, sequences: &[(&str, Vec<u8>)], rows: &[(&str, usize)]) {
    let mut text = String::from("contig\tposition\tref\talt\tpopulation_id\talt_af\n");
    for (contig, position) in rows {
        let reference = sequences
            .iter()
            .find(|(id, _)| id == contig)
            .expect("secuencia")
            .1[position - 1];
        let alternate = *b"ACGT"
            .iter()
            .find(|base| **base != reference)
            .expect("alternativa");
        text.push_str(&format!(
            "{contig}\t{position}\t{}\t{}\tP1\t0.5\n",
            char::from(reference),
            char::from(alternate)
        ));
    }
    fs::write(path, text).expect("matriz");
}

/// Escribe un documento general con referencia suministrada, una población
/// de dos diploides, la matriz suministrada y lo que se pase en `reference`
/// (la sección de detección y los loci) y en `tail`.
fn write_config(
    root: &Path,
    seed: &str,
    sequences: &[(&str, Vec<u8>)],
    rows: &[(&str, usize)],
    reference: &str,
    tail: &str,
) -> PathBuf {
    write_fasta(&root.join("reference.fa"), sequences);
    write_matrix(&root.join("snvs.tsv"), sequences, rows);
    let declared: String = sequences
        .iter()
        .map(|(id, bases)| {
            format!(
                "[[reference.sequences]]\nid = \"{id}\"\nclass = \"chromosome\"\nlength = {}\n\n",
                bases.len()
            )
        })
        .collect();
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{seed}\"\n\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"snvs.tsv\"\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             {declared}{reference}\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n\n{tail}"
        ),
    )
    .expect("configuración");
    config
}

/// Ejecuta `run` por la biblioteca y devuelve la configuración efectiva
/// publicada, o el diagnóstico.
fn run(root: &Path, config: &Path) -> Result<Value, Diagnostic> {
    let effective = load_and_validate_config(config)?;
    let destination = root.join("published");
    execute_validated_config(&effective, &destination)?;
    Ok(serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("efectiva"),
    )
    .expect("JSON"))
}

/// `validate`: la misma frontera que `run`, sin publicar.
fn validate(config: &Path) -> Result<(), Diagnostic> {
    let effective = load_and_validate_config(config)?;
    preflight_validated_config_with_authorizations(&effective, ExecutionAuthorizations::default())
}

/// El catálogo publicado: `(id, origen, secuencia, inicio, fin, motivos)`.
fn catalog(effective: &Value) -> Vec<(String, String, String, u64, u64, String)> {
    effective["repeat_loci"]
        .as_array()
        .expect("repeat_loci")
        .iter()
        .map(|locus| {
            (
                locus["id"].as_str().expect("id").to_owned(),
                locus["origin"].as_str().expect("origin").to_owned(),
                locus["sequence_id"].as_str().expect("sequence").to_owned(),
                locus["start"].as_u64().expect("start"),
                locus["end"].as_u64().expect("end"),
                locus["motifs"]
                    .as_array()
                    .expect("motifs")
                    .iter()
                    .map(|motif| motif.as_str().expect("motivo"))
                    .collect::<Vec<_>>()
                    .join(","),
            )
        })
        .collect()
}

fn entry(
    id: &str,
    origin: &str,
    sequence: &str,
    start: u64,
    end: u64,
    motif: &str,
) -> (String, String, String, u64, u64, String) {
    (
        id.to_owned(),
        origin.to_owned(),
        sequence.to_owned(),
        start,
        end,
        motif.to_owned(),
    )
}

fn detail(diagnostic: &Diagnostic) -> String {
    diagnostic.render_human()
}

/// Seis tractos plantados de periodos 1 a 6 en el suelo por omisión o por
/// encima, y uno por debajo del suelo, en dos secuencias.
fn planted_reference() -> Vec<(&'static str, Vec<u8>)> {
    let mut chr1 = background(0x5DEE_CE66_D1CE_4E5B, 3000);
    plant(&mut chr1, 301, "A", 12);
    plant(&mut chr1, 801, "CA", 10);
    plant(&mut chr1, 1301, "GAT", 7);
    plant(&mut chr1, 1801, "AAAG", 6);
    // Por debajo del suelo de 6 copias de un dinucleótido: no se detecta.
    plant(&mut chr1, 2401, "TG", 5);
    let mut chr2 = background(0x2545_F491_4F6C_DD1D, 2000);
    plant(&mut chr2, 401, "ACGTT", 5);
    plant(&mut chr2, 1201, "AACGTT", 4);
    vec![("chr1", chr1), ("chr2", chr2)]
}

fn planted_catalog() -> Vec<(String, String, String, u64, u64, String)> {
    vec![
        entry("SSR-1-301-312-A", "detected", "chr1", 301, 312, "A"),
        entry("SSR-1-801-820-CA", "detected", "chr1", 801, 820, "CA"),
        entry("SSR-1-1301-1321-GAT", "detected", "chr1", 1301, 1321, "GAT"),
        entry(
            "SSR-1-1801-1824-AAAG",
            "detected",
            "chr1",
            1801,
            1824,
            "AAAG",
        ),
        entry("SSR-2-401-425-ACGTT", "detected", "chr2", 401, 425, "ACGTT"),
        entry(
            "SSR-2-1201-1224-AACGTT",
            "detected",
            "chr2",
            1201,
            1224,
            "AACGTT",
        ),
    ]
}

/// **Una referencia con tractos plantados da exactamente esos loci**, en
/// coordenadas y con identificador escritos aquí, con dos semillas: la
/// detección no consume azar.
#[test]
fn planted_tracts_are_detected_exactly_where_they_were_written() {
    for seed in [SEED_A, SEED_B] {
        let root = temporary("planted");
        let sequences = planted_reference();
        let config = write_config(
            &root,
            seed,
            &sequences,
            &[("chr1", 2000)],
            &format!("[reference.repeat_detection]\n\n{SLIPPAGE}"),
            "",
        );
        let effective = run(&root, &config).expect("la ejecución con detección");
        assert_eq!(catalog(&effective), planted_catalog(), "semilla {seed}");
        // Todos reciben alelos del deslizamiento común.
        for locus in effective["repeat_loci"].as_array().expect("loci") {
            assert_eq!(locus["slippage"]["mean_events"], 1.0, "{locus}");
            assert!(!locus["alleles"].as_array().expect("alelos").is_empty());
        }
        assert_eq!(
            effective["repeat_detection"]["min_copies"],
            serde_json::json!([8, 6, 6, 6, 5, 4])
        );
        assert_eq!(effective["repeat_detection"]["max_units"], 200);
        assert_eq!(effective["repeat_detection"]["regions"], Value::Null);
        assert_eq!(
            effective["repeat_detection"]["algorithm_version"],
            "ssr-detection-v1"
        );
        let _ = fs::remove_dir_all(&root);
    }
}

/// **`validate` y `run` deciden lo mismo**: los dos aceptan el escenario
/// plantado, y los dos rechazan con el mismo diagnóstico, palabra por palabra,
/// un locus declarado que no coincide con ninguno detectado. El diagnóstico
/// nombra el detectado más cercano, de modo que los dos han detectado al
/// menos eso igual; el catálogo entero que ve `validate` no se publica y esta
/// prueba **no** lo compara (el de `run` lo fija
/// `planted_tracts_are_detected_exactly_where_they_were_written`).
#[test]
fn validate_and_run_accept_and_reject_the_same_scenarios() {
    let root = temporary("validate");
    let sequences = planted_reference();
    let reference = format!(
        "[reference.repeat_detection]\n\n{SLIPPAGE}\
         [[reference.repeat_loci]]\nid = \"STR-NONE\"\nsequence = \"chr1\"\nstart = 1001\n\
         end = 1020\nmotif = \"CA\"\n\n"
    );
    let config = write_config(&root, SEED_A, &sequences, &[("chr1", 2000)], &reference, "");
    let validated = validate(&config).expect_err("validate rechaza");
    let ran = run(&root, &config).expect_err("run rechaza");
    assert_eq!(detail(&validated), detail(&ran));
    assert!(
        detail(&ran).contains("STR-NONE"),
        "nombra el locus: {}",
        detail(&ran)
    );
    let accepted = temporary("validate-ok");
    let config = write_config(
        &accepted,
        SEED_A,
        &sequences,
        &[("chr1", 2000)],
        &format!("[reference.repeat_detection]\n\n{SLIPPAGE}"),
        "",
    );
    validate(&config).expect("validate acepta el plantado");
    run(&accepted, &config).expect("run acepta el plantado");
    let _ = fs::remove_dir_all(&root);
    let _ = fs::remove_dir_all(&accepted);
}

/// **Con inserción se detecta sobre la referencia final.** El tracto
/// insertado tras la base 1000 desplaza 30 bases el `(CA)10` que había en
/// 2001..2020 de partida: el detectado está en 2031..2050 y se llama así. Si
/// la detección corriera sobre la de partida, saldría en 2001..2020.
#[test]
fn with_insertion_the_catalog_is_in_final_coordinates() {
    let root = temporary("insertion");
    let mut chr1 = background(0x1234_5678_9ABC_DEF1, 4000);
    plant(&mut chr1, 2001, "CA", 10);
    let sequences = vec![("chr1", chr1)];
    let reference = format!(
        "inserted_reference_id = \"ref-final\"\n\n[reference.repeat_detection]\n\n{SLIPPAGE}\
         [[reference.repeat_loci]]\nid = \"STR-INS\"\nsequence = \"chr1\"\ninsert_after = 1000\n\
         motif = \"CAG\"\nreference_units = [10]\n\n"
    );
    // La matriz suministrada no se admite con inserción: se genera.
    write_fasta(&root.join("reference.fa"), &sequences);
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED_A}\"\n\n\
             [population_snvs]\nmode = \"generate\"\n\n[population_snvs.generate]\n\
             count_mode = \"exact\"\nlocus_count = 3\nspatial_mode = \"uniform\"\n\
             substitution_profile = \"uniform_snv12_v1\"\n\n\
             [[population_snvs.generate.populations]]\npopulation_id = \"P1\"\n\
             frequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\
             {reference}\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 4000\n\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n"
        ),
    )
    .expect("configuración");
    let effective = run(&root, &config).expect("la ejecución con inserción");
    assert_eq!(
        catalog(&effective),
        vec![
            entry("STR-INS", "declared", "chr1", 1001, 1030, "CAG"),
            entry("SSR-1-2031-2050-CA", "detected", "chr1", 2031, 2050, "CA"),
        ]
    );
    let _ = fs::remove_dir_all(&root);
}

/// **Con referencia simulada, el tracto plantado se detecta** y es el
/// declarado quien queda en el catálogo, con su identificador. La región lo
/// rodea para que el relleno al azar no aporte otros.
#[test]
fn with_a_simulated_reference_the_planted_tract_is_detected() {
    let root = temporary("simulated");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED_A}\"\n\n\
             [population_snvs]\nmode = \"generate\"\n\n[population_snvs.generate]\n\
             count_mode = \"exact\"\nlocus_count = 2\nspatial_mode = \"uniform\"\n\
             substitution_profile = \"uniform_snv12_v1\"\n\n\
             [[population_snvs.generate.populations]]\npopulation_id = \"P1\"\n\
             frequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n\n\
             [reference]\nkind = \"simulated\"\nid = \"ref-sim\"\n\n\
             [reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 3000\n\n\
             [reference.repeat_detection]\n\n{SLIPPAGE}\
             [[reference.repeat_detection.regions]]\nsequence = \"chr1\"\nstart = 1490\n\
             end = 1540\n\n\
             [[reference.repeat_loci]]\nid = \"STR-SIM\"\nsequence = \"chr1\"\nstart = 1501\n\
             end = 1530\nmotif = \"TCA\"\n\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n"
        ),
    )
    .expect("configuración");
    let effective = run(&root, &config).expect("la ejecución simulada");
    let found = catalog(&effective);
    assert_eq!(
        found,
        vec![entry("STR-SIM", "declared", "chr1", 1501, 1530, "TCA")],
        "el plantado coincide con un detectado y conserva su nombre y sus fronteras"
    );
    let _ = fs::remove_dir_all(&root);
}

/// **Una `N` suelta dentro de un tracto largo** (apartado 11.6 bis): lo más
/// simple coherente con los criterios es que cuente como discrepancia de
/// pureza, como en la admisión. `(CA)20 N (CA)20` da un solo locus de 81
/// bases (pureza 80/81), y la ejecución publica su VCF, con la `N` dentro del
/// `REF` y de los `ALT`, y lo acepta su juez semántico. Confirmado por el
/// responsable el 30/10/2026. **Debe fallar contra** una `N` que corte el
/// tracto: saldrían dos `(CA)×20`.
#[test]
fn a_lone_n_inside_a_long_tract_counts_as_a_purity_mismatch() {
    let root = temporary("lone-n");
    let mut chr1 = background(0x0F0F_1234_A5A5_9876, 3000);
    plant(&mut chr1, 1001, "CA", 20);
    let (_, end) = plant(&mut chr1, 1042, "CA", 20);
    chr1[1040] = b'N';
    assert_eq!(end, 1081);
    let sequences = vec![("chr1", chr1)];
    let config = write_config(
        &root,
        SEED_A,
        &sequences,
        &[("chr1", 2500)],
        &format!("[reference.repeat_detection]\n\n{SLIPPAGE}"),
        "[vcf_output]\nmode = \"text\"\n",
    );
    let effective = run(&root, &config).expect("la ejecución con una N dentro del tracto");
    assert_eq!(
        catalog(&effective),
        vec![entry(
            "SSR-1-1001-1081-CA",
            "detected",
            "chr1",
            1001,
            1081,
            "CA"
        )]
    );
    let vcf = fs::read_to_string(root.join("published/outputs/truth.vcf")).expect("VCF");
    // El registro del tracto empieza en su primera base y su REF lleva la N.
    assert!(
        vcf.lines().any(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            fields.len() > 7
                && fields[0] == "chr1"
                && fields[1] == "1001"
                && fields[3].len() == 81
                && fields[3].as_bytes()[40] == b'N'
                && fields[7].contains("RN=")
        }),
        "el VCF lleva el registro del tracto con la N:\n{vcf}"
    );
    let _ = fs::remove_dir_all(&root);
}

/// **Regiones (DP-SSR-18):** un tramo se detecta solo si cabe **entero** en
/// una región. Una región que contiene el `(CA)10` de 801..820 da solo ese;
/// una que lo corta por la mitad no da nada; sin lista, todos (lo comprueba
/// `planted_tracts_are_detected_exactly_where_they_were_written`). **Debe
/// fallar contra** conservar un tramo que solo solapa la región: el caso
/// «cuts» daría el `(CA)10`. Decisión del responsable del 30/10/2026.
#[test]
fn regions_keep_only_tracts_that_fit_entirely_inside() {
    let sequences = planted_reference();
    let cases = [
        (
            "contains",
            "[[reference.repeat_detection.regions]]\nsequence = \"chr1\"\nstart = 790\nend = 830\n\n",
            vec![entry(
                "SSR-1-801-820-CA",
                "detected",
                "chr1",
                801,
                820,
                "CA",
            )],
        ),
        (
            "exact",
            "[[reference.repeat_detection.regions]]\nsequence = \"chr1\"\nstart = 801\nend = 820\n\n",
            vec![entry(
                "SSR-1-801-820-CA",
                "detected",
                "chr1",
                801,
                820,
                "CA",
            )],
        ),
        (
            "cuts",
            "[[reference.repeat_detection.regions]]\nsequence = \"chr1\"\nstart = 790\nend = 810\n\n",
            vec![],
        ),
    ];
    for (label, regions, expected) in cases {
        let root = temporary(label);
        let config = write_config(
            &root,
            SEED_A,
            &sequences,
            &[("chr1", 2000)],
            &format!("[reference.repeat_detection]\n\n{SLIPPAGE}{regions}"),
            "",
        );
        let effective = run(&root, &config).unwrap_or_else(|error| panic!("{label}: {error:?}"));
        assert_eq!(catalog(&effective), expected, "{label}");
        let _ = fs::remove_dir_all(&root);
    }
}

/// **Coincidencia de los declarados (REQ-REF-40).**
///
/// - Uno con fronteras distintas pero solapadas coincide y **mandan las
///   suyas**; conserva su nombre y sustituye al detectado.
/// - Uno de motivo `AC` sobre el detectado `CA` coincide: se empareja por
///   rotación, no por igualdad de cadena.
/// - Un compuesto declarado sobre **dos tramos detectados** los sustituye a los
///   dos. Los tramos son un `(CA)×10` en 2601..2620 y un `(GAT)×7` en
///   2623..2643, separados por dos bases de flanco: así el detector da las dos
///   piezas (lo fija el control, sin el declarado), y no se pierde ninguna por
///   la deuda de la frontera entre repeticiones pegadas (DT-061, cerrada en su
///   parte (a) el 01/11/2026 con el recorte de la regla de solapes), que es lo
///   que le pasaba al `(TC)×8 (TCTG)×6` de antes.
///
/// **Debe fallar contra** emparejar por igualdad de cadena del motivo, por
/// igualdad exacta de fronteras, y contra sustituir solo al primer detectado
/// que solapa el compuesto (el `(GAT)×7` seguiría en el catálogo y el
/// catálogo tendría dos loci solapados).
#[test]
fn declared_loci_match_by_overlap_and_rotation_and_keep_their_boundaries() {
    let mut sequences = planted_reference();
    plant(&mut sequences[0].1, 2601, "CA", 10);
    plant(&mut sequences[0].1, 2623, "GAT", 7);
    let section = format!("[reference.repeat_detection]\n\n{SLIPPAGE}");

    // Control positivo: sin declarados salen las dos piezas.
    let root = temporary("matching-control");
    let config = write_config(&root, SEED_A, &sequences, &[("chr1", 2000)], &section, "");
    let effective = run(&root, &config).expect("sin declarados");
    let mut expected = planted_catalog();
    expected.insert(
        4,
        entry("SSR-1-2601-2620-CA", "detected", "chr1", 2601, 2620, "CA"),
    );
    expected.insert(
        5,
        entry("SSR-1-2623-2643-GAT", "detected", "chr1", 2623, 2643, "GAT"),
    );
    assert_eq!(catalog(&effective), expected, "control: las dos piezas");
    let _ = fs::remove_dir_all(&root);

    let root = temporary("matching");
    let reference = format!(
        "{section}\
         [[reference.repeat_loci]]\nid = \"STR-WIDE\"\nsequence = \"chr1\"\nstart = 1299\n\
         end = 1321\nmotif = \"GAT\"\nmax_units = 30\n\n\
         [[reference.repeat_loci]]\nid = \"STR-ROT\"\nsequence = \"chr1\"\nstart = 802\n\
         end = 820\nmotif = \"AC\"\n\n\
         [[reference.repeat_loci]]\nid = \"STR-COMP\"\nsequence = \"chr1\"\nstart = 2601\n\
         end = 2643\nmotifs = [\"CA\", \"GAT\"]\nmin_purity = 0.5\nmin_longest_run_units = 2\n\n\
         [[reference.repeat_loci.alleles]]\nunits = [10, 7]\nfrequency = 0.5\n\n\
         [[reference.repeat_loci.alleles]]\nunits = [10, 8]\nfrequency = 0.5\n\n"
    );
    let config = write_config(&root, SEED_A, &sequences, &[("chr1", 2000)], &reference, "");
    let effective = run(&root, &config).expect("los declarados coinciden");
    let mut expected = planted_catalog();
    expected[1] = entry("STR-ROT", "declared", "chr1", 802, 820, "AC");
    expected[2] = entry("STR-WIDE", "declared", "chr1", 1299, 1321, "GAT");
    expected.insert(
        4,
        entry("STR-COMP", "declared", "chr1", 2601, 2643, "CA,GAT"),
    );
    assert_eq!(catalog(&effective), expected);
    let _ = fs::remove_dir_all(&root);
}

/// Escribe `text` desde la base `start` (base uno) con una `G` en cada flanco
/// (que no prolonga ninguno de los tractos de las pruebas que lo usan).
fn put(bases: &mut [u8], start: usize, text: &str) {
    let begin = start - 1;
    let end = begin + text.len();
    bases[begin..end].copy_from_slice(text.as_bytes());
    bases[begin - 1] = b'G';
    bases[end] = b'G';
}

/// **Tractos pegados: el recorte de la regla de solapes en el producto**
/// (DP-SSR-05, reescrita el 01/11/2026), con el identificador y las
/// coordenadas escritos a mano.
///
/// - `(CA)×10 A×12` en 2201: el buscador da `CA` y una tira de `A` de 13 bases
///   que comparten la última `A` del dinucleótido; la tira se recorta a 12 y
///   son **dos loci**, `SSR-1-2201-2220-CA` y `SSR-1-2221-2232-A`.
/// - `(TCTA)×8 (TCTG)×8 (TCTA)×8` en 2701 (2701..2796): tres trozos que se
///   pisan 3 bases con su vecino, y son **tres loci**. El del centro y el
///   tercero se recortan por la izquierda, y el motivo se lee al inicio del
///   tramo recortado, así que salen **rotaciones**: `GTCT` (de `TCTG`) y
///   `ATCT` (de `TCTA`). El identificador lleva el motivo rotado.
///
/// Con el compuesto **declarado** `(TCTA)(TCTG)(TCTA)`: sigue coincidiendo (su
/// primer motivo es rotación de un detectado, y el del centro también) y
/// **sustituye a los tres** detectados; manda su identificador y su tramo. El
/// control de la parte sin declarar es el que dice que sin él salen tres.
///
/// **Debe fallar contra** la regla de antes (descartar entero): sin recorte
/// faltaría la tira de `A` y el `(TCTG)` central; y contra sustituir solo a
/// los que el compuesto solapa de uno en uno (quedaría algún detectado).
#[test]
fn adjacent_tracts_are_trimmed_and_a_declared_compound_replaces_all_three() {
    let mut sequences = planted_reference();
    put(
        &mut sequences[0].1,
        2201,
        &["CA".repeat(10), "A".repeat(12)].concat(),
    );
    put(
        &mut sequences[0].1,
        2701,
        &["TCTA".repeat(8), "TCTG".repeat(8), "TCTA".repeat(8)].concat(),
    );
    let section = format!("[reference.repeat_detection]\n\n{SLIPPAGE}");
    let mut with_tracts = planted_catalog();
    // chr1 va primero: tras el `AAAG` (índice 3) entran los cinco nuevos.
    for (offset, (id, start, end, motif)) in [
        ("SSR-1-2201-2220-CA", 2201, 2220, "CA"),
        ("SSR-1-2221-2232-A", 2221, 2232, "A"),
        ("SSR-1-2701-2735-TCTA", 2701, 2735, "TCTA"),
        ("SSR-1-2736-2767-GTCT", 2736, 2767, "GTCT"),
        ("SSR-1-2768-2796-ATCT", 2768, 2796, "ATCT"),
    ]
    .into_iter()
    .enumerate()
    {
        with_tracts.insert(4 + offset, entry(id, "detected", "chr1", start, end, motif));
    }

    let root = temporary("trim-control");
    let config = write_config(&root, SEED_A, &sequences, &[("chr1", 2000)], &section, "");
    let effective = run(&root, &config).expect("sin declarados");
    assert_eq!(
        catalog(&effective),
        with_tracts,
        "control: cinco detectados"
    );
    let _ = fs::remove_dir_all(&root);

    let root = temporary("trim-compound");
    let reference = format!(
        "{section}\
         [[reference.repeat_loci]]\nid = \"STR-COMP3\"\nsequence = \"chr1\"\nstart = 2701\n\
         end = 2796\nmotifs = [\"TCTA\", \"TCTG\", \"TCTA\"]\nmin_purity = 0.5\n\
         min_longest_run_units = 2\n\n\
         [[reference.repeat_loci.alleles]]\nunits = [8, 8, 8]\nfrequency = 1.0\n\n"
    );
    let config = write_config(&root, SEED_A, &sequences, &[("chr1", 2000)], &reference, "");
    let effective = run(&root, &config).expect("el compuesto coincide con lo detectado");
    let mut expected = with_tracts.clone();
    expected.splice(
        6..9,
        [entry(
            "STR-COMP3",
            "declared",
            "chr1",
            2701,
            2796,
            "TCTA,TCTG,TCTA",
        )],
    );
    assert_eq!(
        catalog(&effective),
        expected,
        "el compuesto sustituye a los tres"
    );
    let _ = fs::remove_dir_all(&root);
}

/// **Un declarado empareja por rotación con un detectado recortado.** El
/// `(TCTG)×8` central del compuesto `(TCTA)×8 (TCTG)×8 (TCTA)×8` se detecta
/// recortado por la izquierda, con motivo `GTCT` (rotación de `TCTG`) en
/// 2736..2767. Un declarado `TCTG` en 2740..2760 solo solapa a ese detectado y
/// coincide por rotación: conserva su nombre y sus fronteras y lo sustituye;
/// los dos de los lados siguen detectados.
///
/// **Debe fallar contra** emparejar por igualdad de cadena (`GTCT` no es
/// `TCTG`): el locus se rechazaría por no coincidir con ningún detectado.
#[test]
fn a_declared_locus_matches_a_trimmed_detected_one_by_rotation() {
    let mut sequences = planted_reference();
    put(
        &mut sequences[0].1,
        2701,
        &["TCTA".repeat(8), "TCTG".repeat(8), "TCTA".repeat(8)].concat(),
    );
    let root = temporary("trim-rotation");
    let reference = format!(
        "[reference.repeat_detection]\n\n{SLIPPAGE}\
         [[reference.repeat_loci]]\nid = \"STR-CENTER\"\nsequence = \"chr1\"\nstart = 2740\n\
         end = 2760\nmotif = \"TCTG\"\nmax_units = 30\n\n"
    );
    let config = write_config(&root, SEED_A, &sequences, &[("chr1", 2000)], &reference, "");
    let effective = run(&root, &config).expect("el declarado coincide con el recortado");
    let mut expected = planted_catalog();
    expected.insert(
        4,
        entry(
            "SSR-1-2701-2735-TCTA",
            "detected",
            "chr1",
            2701,
            2735,
            "TCTA",
        ),
    );
    expected.insert(
        5,
        entry("STR-CENTER", "declared", "chr1", 2740, 2760, "TCTG"),
    );
    expected.insert(
        6,
        entry(
            "SSR-1-2768-2796-ATCT",
            "detected",
            "chr1",
            2768,
            2796,
            "ATCT",
        ),
    );
    assert_eq!(catalog(&effective), expected);
    let _ = fs::remove_dir_all(&root);
}

/// **Un declarado sin tracto se rechaza** nombrando el locus y lo detectado
/// más cercano; uno por debajo del suelo nombra la etapa que no supera.
#[test]
fn unmatched_declared_loci_are_rejected_naming_the_nearest_and_the_stage() {
    let sequences = planted_reference();
    let cases = [
        (
            "nothing",
            "[[reference.repeat_loci]]\nid = \"STR-NONE\"\nsequence = \"chr1\"\nstart = 1001\n\
             end = 1020\nmotif = \"CA\"\n\n",
            vec!["STR-NONE", "801..820"],
        ),
        (
            "below-floor",
            "[[reference.repeat_loci]]\nid = \"STR-SHORT\"\nsequence = \"chr1\"\nstart = 2401\n\
             end = 2410\nmotif = \"TG\"\n\n",
            vec!["STR-SHORT", "min_copies"],
        ),
    ];
    for (label, locus, needles) in cases {
        let root = temporary(label);
        let config = write_config(
            &root,
            SEED_A,
            &sequences,
            &[("chr1", 2000)],
            &format!("[reference.repeat_detection]\n\n{SLIPPAGE}{locus}"),
            "",
        );
        let error = run(&root, &config).expect_err("se rechaza");
        let text = detail(&error);
        for needle in needles {
            assert!(text.contains(needle), "{label}: falta «{needle}» en {text}");
        }
        let _ = fs::remove_dir_all(&root);
    }
}

/// **Un declarado con motivo de 7 pb queda exento** (REQ-CONF-60): la
/// detección no lo ve, pero tiene que existir en su tramo.
#[test]
fn a_declared_locus_with_a_seven_base_motif_is_exempt() {
    let root = temporary("seven");
    let mut sequences = planted_reference();
    plant(&mut sequences[0].1, 2601, "ACGTTGC", 5);
    let reference = format!(
        "[reference.repeat_detection]\n\n{SLIPPAGE}\
         [[reference.repeat_loci]]\nid = \"STR-SEVEN\"\nsequence = \"chr1\"\nstart = 2601\n\
         end = 2635\nmotif = \"ACGTTGC\"\n\n"
    );
    let config = write_config(&root, SEED_A, &sequences, &[("chr1", 2000)], &reference, "");
    let effective = run(&root, &config).expect("el locus de 7 pb pasa sin detector");
    assert!(
        catalog(&effective).contains(&entry(
            "STR-SEVEN",
            "declared",
            "chr1",
            2601,
            2635,
            "ACGTTGC"
        )),
        "{:?}",
        catalog(&effective)
    );
    let _ = fs::remove_dir_all(&root);
}

/// **Una fila suministrada dentro de un detectado lo retira** y la causa se
/// publica (REQ-CONF-143); la misma fila dentro de un **declarado** se
/// rechaza como siempre (S-32).
#[test]
fn a_supplied_row_inside_a_detected_tract_retires_it_and_inside_a_declared_one_is_rejected() {
    let sequences = planted_reference();
    let root = temporary("retire");
    let config = write_config(
        &root,
        SEED_A,
        &sequences,
        &[("chr1", 2000), ("chr1", 810)],
        &format!("[reference.repeat_detection]\n\n{SLIPPAGE}"),
        "",
    );
    let effective = run(&root, &config).expect("la fila retira el tracto");
    let mut expected = planted_catalog();
    expected.remove(1);
    assert_eq!(catalog(&effective), expected);
    let retired = &effective["repeat_detection"]["retired_loci"];
    assert_eq!(retired.as_array().map(Vec::len), Some(1), "{retired}");
    assert_eq!(retired[0]["id"], "SSR-1-801-820-CA");
    assert_eq!(retired[0]["start"], 801);
    assert_eq!(retired[0]["end"], 820);
    assert!(
        retired[0]["cause"]
            .as_str()
            .is_some_and(|cause| cause.contains("chr1:810")),
        "{retired}"
    );
    let _ = fs::remove_dir_all(&root);

    let root = temporary("retire-declared");
    let config = write_config(
        &root,
        SEED_A,
        &sequences,
        &[("chr1", 2000), ("chr1", 810)],
        &format!(
            "[reference.repeat_detection]\n\n{SLIPPAGE}\
             [[reference.repeat_loci]]\nid = \"STR-CA\"\nsequence = \"chr1\"\nstart = 801\n\
             end = 820\nmotif = \"CA\"\n\n"
        ),
        "",
    );
    let error = run(&root, &config).expect_err("dentro de un declarado se rechaza");
    assert!(
        detail(&error).contains("STR-CA") && detail(&error).contains("chr1:810"),
        "{}",
        detail(&error)
    );
    let _ = fs::remove_dir_all(&root);
}

/// Bases de una referencia con un 0,7 % de tractos plantados: 70 tractos
/// `(CA)10` de 20 bases en 200 000 bases, uno cada 2 857.
/// Secuencias con nombre y tramos plantados `(inicio, fin)` en base uno.
type PlantedReference = (Vec<(&'static str, Vec<u8>)>, Vec<(u64, u64)>);

fn dense_reference() -> PlantedReference {
    let mut chr1 = background(0x7777_1111_3333_9999, 200_000);
    let mut tracts = Vec::new();
    for index in 0..70 {
        let start = 1000 + index * 2857;
        let (start, end) = plant(&mut chr1, start, "CA", 10);
        tracts.push((start as u64, end as u64));
    }
    (vec![("chr1", chr1)], tracts)
}

/// **Una matriz generada de 1000 SNV nunca cae dentro de un locus**, con dos
/// semillas (riesgo R-5): el generador excluye los tractos del catálogo, de
/// modo que nada se retira y nada depende de la semilla.
#[test]
fn a_generated_matrix_never_falls_inside_a_locus_with_two_seeds() {
    let (sequences, tracts) = dense_reference();
    for seed in [SEED_A, SEED_B] {
        let root = temporary("dense");
        write_fasta(&root.join("reference.fa"), &sequences);
        let config = root.join("general.toml");
        fs::write(
            &config,
            format!(
                "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{seed}\"\n\n\
                 [population_snvs]\nmode = \"generate\"\n\n[population_snvs.generate]\n\
                 count_mode = \"exact\"\nlocus_count = 1000\nspatial_mode = \"uniform\"\n\
                 substitution_profile = \"uniform_snv12_v1\"\n\n\
                 [[population_snvs.generate.populations]]\npopulation_id = \"P1\"\n\
                 frequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n\n\
                 [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
                 [reference.repeat_detection]\n\n{SLIPPAGE}\
                 [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 200000\n\n\
                 [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\n\
                 ploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n"
            ),
        )
        .expect("configuración");
        let effective = run(&root, &config).expect("la matriz generada se acepta");
        // Los setenta plantados quedan en el catálogo. El fondo de 200 000
        // bases puede aportar algún tracto más de periodo 3 o mayor: no se
        // cuenta, pero también se excluye abajo.
        let loci: Vec<(u64, u64)> = catalog(&effective)
            .into_iter()
            .map(|(_, _, _, start, end, _)| (start, end))
            .collect();
        for tract in &tracts {
            assert!(loci.contains(tract), "{seed}: falta el plantado {tract:?}");
        }
        assert_eq!(
            effective["repeat_detection"]["retired_loci"],
            serde_json::json!([]),
            "{seed}"
        );
        let matrix =
            fs::read_to_string(root.join("published/population_snvs.tsv")).expect("matriz");
        let positions: Vec<u64> = matrix
            .lines()
            .skip(1)
            .map(|line| {
                line.split('\t')
                    .nth(2)
                    .expect("posición")
                    .parse()
                    .expect("u64")
            })
            .collect();
        assert_eq!(positions.len(), 1000, "{seed}");
        for position in positions {
            assert!(
                !loci
                    .iter()
                    .any(|(start, end)| (*start..=*end).contains(&position)),
                "{seed}: la fila de {position} cae dentro de un locus del catálogo"
            );
        }
        let _ = fs::remove_dir_all(&root);
    }
}

/// **La exclusión empieza en la primera base del tracto.** El generador solo
/// puede elegir dentro de la región 801..821, que es el `(CA)10` de 801..820
/// más una base: queda **una** posición elegible. Con un locus, la fila cae en
/// la 821; con dos, no hay sitio y se rechaza. Si la exclusión dejara fuera la
/// primera base del tracto, habría dos posiciones y el segundo caso pasaría.
#[test]
fn the_generator_exclusion_starts_at_the_first_base_of_the_tract() {
    let sequences = planted_reference();
    for (count, expected) in [(1, Some(821_u64)), (2, None)] {
        let root = temporary("first-base");
        write_fasta(&root.join("reference.fa"), &sequences);
        let config = root.join("general.toml");
        fs::write(
            &config,
            format!(
                "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED_A}\"\n\n\
                 [population_snvs]\nmode = \"generate\"\n\n[population_snvs.generate]\n\
                 count_mode = \"exact\"\nlocus_count = {count}\n\
                 spatial_mode = \"weighted_regions\"\noutside_weight = 0.0\n\
                 substitution_profile = \"uniform_snv12_v1\"\n\n\
                 [[population_snvs.generate.regions]]\ncontig = \"chr1\"\nstart = 801\n\
                 end = 821\nweight = 1.0\n\n\
                 [[population_snvs.generate.populations]]\npopulation_id = \"P1\"\n\
                 frequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n\n\
                 [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
                 [reference.repeat_detection]\n\n{SLIPPAGE}\
                 [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 3000\n\n\
                 [[reference.sequences]]\nid = \"chr2\"\nclass = \"chromosome\"\nlength = 2000\n\n\
                 [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\n\
                 ploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n"
            ),
        )
        .expect("configuración");
        match (run(&root, &config), expected) {
            (Ok(_), Some(position)) => {
                let matrix =
                    fs::read_to_string(root.join("published/population_snvs.tsv")).expect("matriz");
                let positions: Vec<&str> = matrix
                    .lines()
                    .skip(1)
                    .map(|line| line.split('\t').nth(2).expect("posición"))
                    .collect();
                assert_eq!(positions, vec![position.to_string().as_str()]);
            }
            (Err(error), None) => assert!(
                detail(&error).contains("population_snvs"),
                "{}",
                detail(&error)
            ),
            (outcome, expected) => panic!(
                "con {count} locus se esperaba {expected:?} y salió {:?}",
                outcome.map(|_| ())
            ),
        }
        let _ = fs::remove_dir_all(&root);
    }
}

/// **Sin individuos** (REQ-REF-37: «sin matriz» es sin individuos). Un
/// escenario con la sección y un solo locus declarado con sus alelos se
/// publica, y el juez físico, que detecta otra vez sobre la referencia
/// publicada, necesita sus bases aunque no haya población ni evento: sin
/// ellas la reconstrucción se rompía con un defecto interno (cazado al
/// reescribir `tests/repeat_admission_guard.rs`).
///
/// **Sin individuos no se generan alelos de los loci detectados** (decisión
/// del responsable del 30/10/2026): sin copias haplotípicas el generador no
/// tiene linajes (REQ-CONF-140), y un detectado no tiene otros alelos, de modo
/// que **no entra en el catálogo** y la ejecución se publica. La segunda mitad
/// lo fija sobre la referencia entera, que con individuos da los loci de
/// `planted_catalog` (`planted_tracts_are_detected_exactly_where_they_were_written`):
/// sin ellos el catálogo publicado queda vacío y la sección sí se publica.
/// **Debe fallar contra** intentar generar los alelos de los detectados, que
/// era lo de antes: el generador rechaza el escenario («no declara
/// ninguna»).
#[test]
fn without_individuals_a_declared_locus_publishes_and_detected_loci_are_left_out() {
    let sequences = planted_reference();
    let header = |reference: &str| {
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED_A}\"\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 3000\n\n\
             [[reference.sequences]]\nid = \"chr2\"\nclass = \"chromosome\"\nlength = 2000\n\n\
             [reference.repeat_detection]\n\n{SLIPPAGE}{reference}"
        )
    };
    let root = temporary("no-individuals");
    write_fasta(&root.join("reference.fa"), &sequences);
    let config = root.join("general.toml");
    fs::write(
        &config,
        header(
            "[[reference.repeat_detection.regions]]\nsequence = \"chr1\"\nstart = 795\nend = 826\n\n\
             [[reference.repeat_loci]]\nid = \"STR-CA\"\nsequence = \"chr1\"\nstart = 801\n\
             end = 820\nmotif = \"CA\"\n\n[[reference.repeat_loci.alleles]]\nunits = 10\n\
             frequency = 1.0\n\n",
        ),
    )
    .expect("configuración");
    let effective = run(&root, &config).expect("sin individuos, con el locus declarado");
    assert_eq!(
        catalog(&effective),
        vec![entry("STR-CA", "declared", "chr1", 801, 820, "CA")]
    );
    let _ = fs::remove_dir_all(&root);

    let root = temporary("no-individuals-detected");
    write_fasta(&root.join("reference.fa"), &sequences);
    let config = root.join("general.toml");
    fs::write(&config, header("")).expect("configuración");
    let effective = run(&root, &config)
        .unwrap_or_else(|error| panic!("sin individuos, sin loci: {}", detail(&error)));
    assert_eq!(catalog(&effective), vec![]);
    assert!(
        effective["repeat_detection"].is_object(),
        "la sección de microsatélites se publica: {effective}"
    );
    let _ = fs::remove_dir_all(&root);
}

/// **Con PCR, el tartamudeo de base solo se exige si el catálogo tiene algún
/// locus** (decisión del responsable del 30/10/2026; REQ-CONF-138 con «loci
/// declarados» como los del catálogo, REQ-CONF-142). La misma biblioteca
/// amplificada sin tartamudeo de base se admite con una región en la que no
/// se detecta nada y se rechaza, nombrando el primer locus, con una región que
/// contiene el `(CA)×10` plantado. Las dos por `validate`, que detecta igual
/// que `run`.
///
/// **Debe fallar contra** exigir el tartamudeo por la mera presencia de la
/// sección (la primera mitad se rechazaría) y contra no volver a comprobarlo
/// después de detectar (la segunda se admitiría).
#[test]
fn with_pcr_base_stutter_is_required_only_when_the_catalog_has_a_locus() {
    let sequences = planted_reference();
    let amplified = "[technical_event]\nevent_id = \"reads-main\"\n\
         individual_ids = [\"I1\", \"I2\"]\ntarget_coverage_wgs = 2\n\
         library_preparation = \"pcr\"\n\n\
         [technical_event.pcr]\nduplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n\n\
         [[technical_event.libraries]]\nindividual_id = \"I1\"\ni7 = \"ACGTACGT\"\n\
         i5 = \"TGCATGCA\"\n\n\
         [[technical_event.libraries]]\nindividual_id = \"I2\"\ni7 = \"CGTACGTA\"\n\
         i5 = \"GCATGCAT\"\n";
    let region = |start: u64, end: u64| {
        format!(
            "[reference.repeat_detection]\n\n{SLIPPAGE}[[reference.repeat_detection.regions]]\n\
             sequence = \"chr1\"\nstart = {start}\nend = {end}\n\n"
        )
    };

    let root = temporary("pcr-empty");
    let config = write_config(
        &root,
        SEED_A,
        &sequences,
        &[("chr1", 2500)],
        &region(2, 100),
        amplified,
    );
    validate(&config).unwrap_or_else(|error| {
        panic!(
            "sin ningún locus, la PCR sin tartamudeo se admite: {}",
            detail(&error)
        )
    });
    let _ = fs::remove_dir_all(&root);

    let root = temporary("pcr-locus");
    let config = write_config(
        &root,
        SEED_A,
        &sequences,
        &[("chr1", 2500)],
        &region(790, 830),
        amplified,
    );
    let error = validate(&config).expect_err("con un locus detectado la PCR exige tartamudeo");
    let rendered = detail(&error);
    assert!(
        rendered.contains("exige technical_event.stutter.enabled = true")
            && rendered.contains("SSR-1-801-820-CA"),
        "{rendered}"
    );
    let _ = fs::remove_dir_all(&root);
}

/// La biblioteca amplificada sin tartamudeo de base de las pruebas de la PCR.
const AMPLIFIED: &str = "[technical_event]\nevent_id = \"reads-main\"\n\
     individual_ids = [\"I1\", \"I2\"]\ntarget_coverage_wgs = 2\n\
     library_preparation = \"pcr\"\n\n\
     [technical_event.pcr]\nduplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n\n\
     [[technical_event.libraries]]\nindividual_id = \"I1\"\ni7 = \"ACGTACGT\"\n\
     i5 = \"TGCATGCA\"\n\n\
     [[technical_event.libraries]]\nindividual_id = \"I2\"\ni7 = \"CGTACGTA\"\n\
     i5 = \"GCATGCAT\"\n";

/// Los ficheros publicados de un directorio, con sus bytes, por ruta.
fn published_files(root: &Path) -> std::collections::BTreeMap<String, Vec<u8>> {
    fn walk(root: &Path, directory: &Path, out: &mut std::collections::BTreeMap<String, Vec<u8>>) {
        for entry in fs::read_dir(directory).expect("directorio publicado") {
            let path = entry.expect("entrada").path();
            if path.is_dir() {
                walk(root, &path, out);
            } else {
                out.insert(
                    path.strip_prefix(root)
                        .expect("relativa")
                        .to_string_lossy()
                        .into_owned(),
                    fs::read(&path).expect("fichero"),
                );
            }
        }
    }
    let mut out = std::collections::BTreeMap::new();
    walk(root, root, &mut out);
    out
}

/// **Con PCR y el catálogo final vacío, `run` publica lo mismo que sin la
/// sección** (DTC-092; decisión del responsable del 30/10/2026). Dos maneras
/// de quedarse sin loci: una región en la que no se detecta nada, y una fila
/// de la matriz que retira el único detectado, el `(CA)×10` de 801..820
/// (REQ-CONF-143). En las dos la ejecución se publica, el evento técnico lleva
/// el catálogo de azar `f4-technical-v5` de la amplificación sin
/// microsatélites, y los FASTQ y el manifiesto de evento son, byte a byte, los
/// de la misma ejecución sin la sección.
///
/// **Debe fallar contra** dejar conectado el tartamudeo clonal que se conectó
/// al leer la sección (el evento saldría con `f4-technical-v6`, el manifiesto
/// general con la versión 4, y la publicación moría con GPS-FMT-006), y
/// contra aplicar la regla de la PCR antes de retirar (la segunda ejecución se
/// rechazaría por falta de tartamudeo).
#[test]
fn with_pcr_an_empty_final_catalog_publishes_what_a_run_without_the_section_does() {
    let sequences = planted_reference();
    let section = |start: u64, end: u64| {
        format!(
            "[reference.repeat_detection]\n\n{SLIPPAGE}[[reference.repeat_detection.regions]]\n\
             sequence = \"chr1\"\nstart = {start}\nend = {end}\n\n"
        )
    };
    for (label, region, row) in [
        ("empty-region", (2, 100), ("chr1", 2500)),
        ("retired", (790, 830), ("chr1", 810)),
    ] {
        let without = temporary(&format!("pcr-{label}-without"));
        let config = write_config(&without, SEED_A, &sequences, &[row], "", AMPLIFIED);
        run(&without, &config)
            .unwrap_or_else(|error| panic!("{label}, sin la sección: {}", detail(&error)));
        let with = temporary(&format!("pcr-{label}-with"));
        let config = write_config(
            &with,
            SEED_A,
            &sequences,
            &[row],
            &section(region.0, region.1),
            AMPLIFIED,
        );
        let effective = run(&with, &config)
            .unwrap_or_else(|error| panic!("{label}, con la sección: {}", detail(&error)));
        assert_eq!(catalog(&effective), vec![], "{label}");
        let retired = effective["repeat_detection"]["retired_loci"]
            .as_array()
            .expect("retirados")
            .len();
        assert_eq!(retired, usize::from(label == "retired"), "{label}");
        let event: Value = serde_json::from_slice(
            &fs::read(with.join("published/events/reads-main/manifest.json")).expect("evento"),
        )
        .expect("JSON");
        assert_eq!(event["rng_catalog"], "f4-technical-v5", "{label}");
        assert_eq!(
            event["artifacts"]["pcr"]["rng_catalog"], "f4-technical-v5",
            "{label}"
        );
        let reference = published_files(&without.join("published"));
        let observed = published_files(&with.join("published"));
        let compared: Vec<&String> = reference
            .keys()
            .filter(|path| path.ends_with(".fastq.gz") || path.starts_with("events/"))
            .collect();
        assert!(
            compared
                .iter()
                .filter(|path| path.ends_with(".fastq.gz"))
                .count()
                == 4,
            "{label}: los cuatro FASTQ: {compared:?}"
        );
        for path in compared {
            assert!(
                observed.get(path) == reference.get(path),
                "{label}: {path} cambia con la sección"
            );
        }
        let _ = fs::remove_dir_all(&without);
        let _ = fs::remove_dir_all(&with);
    }
}

/// **El juez físico rechaza una efectiva sin su sección de microsatélites por
/// la sección, y no con un defecto interno.** Un paquete sin población (la
/// referencia se lee entonces sin bases si nada de la efectiva las pide) al
/// que se le quita `repeat_detection` de la efectiva, rehaciendo su huella en
/// el manifiesto. La reconstrucción detecta según el TOML copiado, que sí
/// lleva la sección, y necesita las bases.
///
/// **Debe fallar contra** decidir si se leen las bases solo por la efectiva:
/// la reconstrucción se rompía con GPS-INTERNAL-001 («bases no
/// materializadas»).
#[test]
fn the_physical_judge_rejects_an_effective_config_without_its_detection_section() {
    use genome_population_simulator::integrity::identify_file;
    use genome_population_simulator::results::{json_document_bytes, validate_published_resources};
    let sequences = planted_reference();
    let root = temporary("judge-without-section");
    write_fasta(&root.join("reference.fa"), &sequences);
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED_A}\"\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 3000\n\n\
             [[reference.sequences]]\nid = \"chr2\"\nclass = \"chromosome\"\nlength = 2000\n\n\
             [reference.repeat_detection]\n\n{SLIPPAGE}\
             [[reference.repeat_detection.regions]]\nsequence = \"chr1\"\nstart = 795\nend = 826\n\n\
             [[reference.repeat_loci]]\nid = \"STR-CA\"\nsequence = \"chr1\"\nstart = 801\n\
             end = 820\nmotif = \"CA\"\n\n[[reference.repeat_loci.alleles]]\nunits = 10\n\
             frequency = 1.0\n\n"
        ),
    )
    .expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let judge = |effective: &Value, manifest: &Value| {
        let effective_bytes = json_document_bytes(effective).expect("JSON");
        fs::write(&effective_path, &effective_bytes).expect("efectiva");
        let identity = identify_file(&effective_path, "effective-config.json").expect("huella");
        let mut manifest = manifest.clone();
        for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
            if artifact["path"] == "effective-config.json" {
                artifact["checksum"]["digest"] = Value::String(identity.digest().to_owned());
                artifact["size_bytes"] = Value::from(identity.size_bytes());
            }
        }
        let manifest_bytes = json_document_bytes(&manifest).expect("JSON");
        fs::write(&manifest_path, &manifest_bytes).expect("manifiesto");
        validate_published_resources(
            &destination,
            &effective_path,
            &effective_bytes,
            &manifest_path,
            &manifest_bytes,
            None,
        )
    };
    let original: Value =
        serde_json::from_slice(&fs::read(&effective_path).expect("efectiva")).expect("JSON");
    let manifest: Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("manifiesto")).expect("JSON");
    judge(&original, &manifest).expect("control positivo: el paquete rehecho igual se acepta");
    let mut without = original.clone();
    without
        .as_object_mut()
        .expect("objeto")
        .remove("repeat_detection")
        .expect("la efectiva publicaba la sección");
    let error = judge(&without, &manifest).expect_err("sin la sección se rechaza");
    let rendered = detail(&error);
    assert!(
        error.code() != "GPS-INTERNAL-001"
            && rendered.contains(
                "la sección de microsatélites efectiva no es la que el TOML, la referencia y la matriz producen"
            ),
        "{rendered}"
    );
    let _ = fs::remove_dir_all(&root);
}

/// **Un declarado no puede llevar el identificador de un detectado que
/// queda** (REQ-REF-40, DTC-096): `SSR-1-301-312-A` declarado sobre el
/// `(CA)×10` de 801..820 choca con el detectado `A×12` de 301..312, que sigue
/// en el catálogo, y se rechaza como identificador repetido. **Debe fallar
/// contra** no comprobar el choque: saldrían dos loci con el mismo nombre.
#[test]
fn a_declared_locus_named_like_a_remaining_detected_one_is_rejected() {
    let root = temporary("duplicate-id");
    let sequences = planted_reference();
    let reference = format!(
        "[reference.repeat_detection]\n\n{SLIPPAGE}\
         [[reference.repeat_loci]]\nid = \"SSR-1-301-312-A\"\nsequence = \"chr1\"\nstart = 801\n\
         end = 820\nmotif = \"CA\"\n\n"
    );
    let config = write_config(&root, SEED_A, &sequences, &[("chr1", 2000)], &reference, "");
    let error = run(&root, &config).expect_err("identificador repetido");
    let rendered = detail(&error);
    assert!(
        error.code() == "GPS-CONF-008"
            && rendered.contains("identificador de locus repetido: SSR-1-301-312-A"),
        "{rendered}"
    );
    let _ = fs::remove_dir_all(&root);
}

/// Una ejecución con un `(CAG)×10` insertado tras la base 1000 (30 bases) y un
/// `(CA)×10` plantado en 2001..2020 de partida, que en la referencia final
/// ocupa 2031..2050, con las regiones dadas.
fn insertion_with_regions(regions: &str) -> Result<Value, Diagnostic> {
    let root = temporary("insertion-regions");
    let mut chr1 = background(0x1234_5678_9ABC_DEF1, 4000);
    plant(&mut chr1, 2001, "CA", 10);
    write_fasta(&root.join("reference.fa"), &[("chr1", chr1)]);
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED_A}\"\n\n\
             [population_snvs]\nmode = \"generate\"\n\n[population_snvs.generate]\n\
             count_mode = \"exact\"\nlocus_count = 3\nspatial_mode = \"uniform\"\n\
             substitution_profile = \"uniform_snv12_v1\"\n\n\
             [[population_snvs.generate.populations]]\npopulation_id = \"P1\"\n\
             frequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\
             inserted_reference_id = \"ref-final\"\n\n\
             [reference.repeat_detection]\n\n{SLIPPAGE}{regions}\
             [[reference.repeat_loci]]\nid = \"STR-INS\"\nsequence = \"chr1\"\ninsert_after = 1000\n\
             motif = \"CAG\"\nreference_units = [10]\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 4000\n\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n"
        ),
    )
    .expect("configuración");
    let outcome = run(&root, &config);
    let _ = fs::remove_dir_all(&root);
    outcome
}

/// **Con inserción, las regiones van en coordenadas de la referencia final**
/// (DP-SSR-18, paso F6.4.a.10b). Una región en 2025..4030 contiene el
/// `(CA)×10`, que en la final ocupa 2031..2050, y llega hasta la última base
/// de la final, que mide 4030; la misma región en las coordenadas de partida
/// del tracto (1995..2025) no lo contiene. **Debe fallar contra** resolver las
/// regiones en coordenadas de partida: la primera se rechazaría por pasarse
/// de la secuencia de 4000 bases, y la segunda daría el `(CA)×10`.
#[test]
fn with_insertion_the_regions_are_in_final_coordinates() {
    let region = |start: u64, end: u64| {
        format!(
            "[[reference.repeat_detection.regions]]\nsequence = \"chr1\"\nstart = {start}\n\
             end = {end}\n\n"
        )
    };
    let inserted = region(990, 1040);
    let effective = insertion_with_regions(&format!("{inserted}{}", region(2025, 4030)))
        .unwrap_or_else(|error| panic!("regiones finales: {}", detail(&error)));
    assert_eq!(
        catalog(&effective),
        vec![
            entry("STR-INS", "declared", "chr1", 1001, 1030, "CAG"),
            entry("SSR-1-2031-2050-CA", "detected", "chr1", 2031, 2050, "CA"),
        ]
    );
    let effective = insertion_with_regions(&format!("{inserted}{}", region(1995, 2025)))
        .unwrap_or_else(|error| panic!("regiones de partida: {}", detail(&error)));
    assert_eq!(
        catalog(&effective),
        vec![entry("STR-INS", "declared", "chr1", 1001, 1030, "CAG")]
    );
}

/// **Un cambio de base escrito a mano dentro de un detectado lo retira**
/// (REQ-CONF-143), igual que una fila de la matriz: el SNV `snp-direct` en
/// chr1:810 cae en el `(CA)×10` de 801..820, que sale del catálogo y se
/// publica entre los retirados con su causa. **Debe fallar contra** mirar
/// solo la matriz: el locus quedaría con una variante dentro y la ejecución se
/// rechazaría.
#[test]
fn a_direct_base_change_inside_a_detected_tract_retires_it() {
    let root = temporary("direct-retire");
    let sequences = planted_reference();
    let reference_base = char::from(sequences[0].1[809]);
    let alternate = if reference_base == 'T' { 'G' } else { 'T' };
    let tail = format!(
        "[[direct_events]]\nevent_id = \"snp-direct\"\nindividual_id = \"I1\"\nphase = 0\n\
         kind = \"snp\"\ncontig = \"chr1\"\nposition_1based = 810\nref = \"{reference_base}\"\n\
         alt = \"{alternate}\"\n"
    );
    let config = write_config(
        &root,
        SEED_A,
        &sequences,
        &[("chr1", 2000)],
        &format!("[reference.repeat_detection]\n\n{SLIPPAGE}"),
        &tail,
    );
    let effective = run(&root, &config)
        .unwrap_or_else(|error| panic!("el cambio directo retira: {}", detail(&error)));
    let mut expected = planted_catalog();
    expected.remove(1);
    assert_eq!(catalog(&effective), expected);
    let retired = &effective["repeat_detection"]["retired_loci"];
    assert_eq!(retired.as_array().map(Vec::len), Some(1), "{retired}");
    assert_eq!(retired[0]["id"], "SSR-1-801-820-CA", "{retired}");
    assert!(
        retired[0]["cause"]
            .as_str()
            .is_some_and(|cause| cause.contains("snp-direct") && cause.contains("chr1:810")),
        "{retired}"
    );
    let _ = fs::remove_dir_all(&root);
}
