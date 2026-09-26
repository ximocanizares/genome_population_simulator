//! **Un diagnóstico cita el esquema del documento que rechaza, no uno fijo.**
//!
//! # Qué arregla, y por qué merecía una prueba
//!
//! Dos diagnósticos de formato citaban **siempre** el mismo fichero de esquema
//! aunque el documento fuera de otra versión: el del manifiesto de evento decía
//! siempre `manifest-event-v2.schema.json` y el de la configuración efectiva
//! siempre `effective-config-v4.schema.json`. Es cosmético en el sentido de que
//! no cambia ningún resultado, y **no** lo es en el sentido de que manda a quien
//! lee el rechazo a abrir el fichero equivocado. En un proyecto cuya regla es
//! que los documentos digan la verdad, eso chirría.
//!
//! # Por qué las dos mitades
//!
//! Cada prueba comprueba **las dos versiones**, no solo la nueva. Comprobar
//! solo la v3 y la v6 dejaría pasar una implementación que citara siempre la
//! versión alta, que es el mismo defecto con el error al otro lado.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::event_manifest::validate_event_manifest_bytes;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::validate_effective_config_against_manifest_bytes;
use serde_json::Value;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-diag-esquema-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

fn bases(length: usize) -> String {
    let mut state: u64 = 0x5DEE_CE66_D1B2_7F3C;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            match (state >> 33) % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            }
        })
        .collect()
}

struct Published {
    root: PathBuf,
    destination: PathBuf,
}

impl Drop for Published {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

impl Published {
    fn read(&self, relative: &str) -> Value {
        serde_json::from_slice(&fs::read(self.destination.join(relative)).expect("documento"))
            .expect("JSON")
    }

    fn path(&self, relative: &str) -> PathBuf {
        self.destination.join(relative)
    }
}

/// Publica un escenario con o sin preparación amplificada.
///
/// El escenario amplificado sube el documento general a la versión que conoce
/// las repeticiones, y con él la configuración efectiva a la 6 y el manifiesto
/// de evento a la 3. El sencillo se queda en la 4 y la 2.
fn run(amplified: bool) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let sequence = bases(8_000);
    let reference_base = sequence.as_bytes()[19] as char;
    let alternate = if reference_base == 'T' { 'G' } else { 'T' };
    fs::write(root.join("reference.fa"), format!(">chr1\n{sequence}\n")).expect("FASTA");
    fs::write(
        root.join("population-snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\n\
             chr1\t20\t{reference_base}\t{alternate}\tP1\t0.5\n"
        ),
    )
    .expect("matriz");
    let (version, extra) = if amplified {
        (
            "4",
            "library_preparation = \"pcr\"\n\n[technical_event.pcr]\n\
             duplicate_fraction = 0.2\ncycles = 10\nefficiency = 0.9\n",
        )
    } else {
        ("2", "")
    };
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "{version}"
seed = "{SEED}"
ploidy = 1

[population_snvs]
mode = "supplied"
input_path = "population-snvs.tsv"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 8000

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 4
{extra}
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        ),
    )
    .expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Published { root, destination }
}

/// Rompe un documento de forma que su validador tenga que rechazarlo.
fn corrupted(mut document: Value) -> Vec<u8> {
    document
        .as_object_mut()
        .expect("objeto")
        .insert("campo_inventado".to_owned(), Value::from(1));
    serde_json::to_vec(&document).expect("JSON")
}

fn cited_schema(error: &genome_population_simulator::diagnostics::Diagnostic) -> String {
    assert_eq!(error.code(), "GPS-FMT-006", "{error:?}");
    error
        .context()
        .get("schema")
        .cloned()
        .expect("el diagnóstico cita un esquema")
}

/// **El manifiesto de evento cita el esquema de su versión.**
#[test]
fn an_event_manifest_diagnostic_cites_the_schema_of_its_declared_version() {
    let relative = Path::new("events/reads-main/manifest.json");
    for (amplified, version, schema) in [
        (false, "2", "manifest-event-v2.schema.json"),
        (true, "3", "manifest-event-v3.schema.json"),
    ] {
        let published = run(amplified);
        let document = published.read("events/reads-main/manifest.json");
        // Control positivo: los dos escenarios producen versiones distintas de
        // verdad. Sin esto, la comparación de abajo podría estar mirando dos
        // veces el mismo documento.
        assert_eq!(
            document["schema_version"].as_str(),
            Some(version),
            "el escenario amplificado={amplified} tenía que publicar la versión {version}"
        );
        let error = validate_event_manifest_bytes(relative, &corrupted(document))
            .expect_err("un campo inventado tiene que rechazarse");
        assert_eq!(
            cited_schema(&error),
            schema,
            "un manifiesto de evento de la versión {version} tiene que citar {schema}"
        );
    }
}

/// **La configuración efectiva cita el esquema de su versión.**
#[test]
fn an_effective_config_diagnostic_cites_the_schema_of_its_declared_version() {
    for (amplified, version, schema) in [
        (false, "4", "effective-config-v4.schema.json"),
        (true, "6", "effective-config-v6.schema.json"),
    ] {
        let published = run(amplified);
        let effective = published.read("effective-config.json");
        assert_eq!(
            effective["schema_version"].as_str(),
            Some(version),
            "el escenario amplificado={amplified} tenía que publicar la versión {version}"
        );
        let manifest_bytes = fs::read(published.path("manifest.json")).expect("manifiesto");
        let error = validate_effective_config_against_manifest_bytes(
            &published.path("effective-config.json"),
            &corrupted(effective),
            &published.path("manifest.json"),
            &manifest_bytes,
        )
        .expect_err("un campo inventado tiene que rechazarse");
        assert_eq!(
            cited_schema(&error),
            schema,
            "una configuración efectiva de la versión {version} tiene que citar {schema}"
        );
    }
}

/// **Un rechazo del manifiesto general no se retiqueta con el esquema de la
/// configuración efectiva.**
///
/// Es la mitad que el retiquetado podría estropear: las dos familias de
/// diagnóstico viajan por el mismo resultado, y reescribir por el código solo
/// cambiaría un destino equivocado por otro.
#[test]
fn a_general_manifest_rejection_keeps_its_own_schema() {
    let published = run(true);
    let effective_bytes = fs::read(published.path("effective-config.json")).expect("efectiva");
    let manifest = published.read("manifest.json");
    let error = validate_effective_config_against_manifest_bytes(
        &published.path("effective-config.json"),
        &effective_bytes,
        &published.path("manifest.json"),
        &corrupted(manifest),
    )
    .expect_err("un campo inventado en el manifiesto tiene que rechazarse");
    let schema = cited_schema(&error);
    assert!(
        schema.starts_with("manifest-general-"),
        "un rechazo del manifiesto general tiene que citar su propio esquema, y cita {schema}"
    );
}
