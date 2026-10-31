//! DTC-099: el multiplicador configurable de la PCR en los homopolímeros.
//!
//! # Qué hay aquí
//!
//! 1. **La huella congelada sin el campo** (REQ-CONF-141, ACE-120): un
//!    escenario amplificado con un homopolímero y un dinucleótido, sin
//!    multiplicador, publica exactamente lo que publicaba antes de que el campo
//!    existiera.
//! 2. **La publicación con el campo** (REQ-MAN-97): los dos documentos lo
//!    llevan, la fila del homopolímero lleva la fracción resultante y nada más
//!    se mueve; sin él, la clave no aparece.
//! 3. **El juez** (REQ-VAL-61): las manipulaciones coherentes de fracción y
//!    probabilidad se rechazan en la configuración efectiva y en el manifiesto
//!    de evento.
//!
//! Los vectores a mano están en `src/clonal_stutter.rs`, la composición sobre
//! los FASTQ en `tests/clonal_composition.rs` y los rechazos de admisión en
//! `tests/pcr_admission.rs`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::clonal_stutter::{
    HOMOPOLYMER_PCR_FOLD_MAXIMUM, HOMOPOLYMER_PCR_FOLD_MINIMUM, per_cycle_slip_probability,
};
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::validate_published_resources;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

#[path = "support/detection_section.rs"]
mod detection_section;

static NONCE: AtomicU64 = AtomicU64::new(0);

/// El escenario es **el mismo fichero** que la evidencia de la comparación
/// entre binarios, `evidence/dt-009-multiplicador-congelacion/`: así la prueba
/// y la evidencia no pueden divergir.
const SCENARIO_TOML: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/general.toml"
);
const SCENARIO_FASTA: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/reference.fa"
);
const SCENARIO_SNVS: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/population-snvs.tsv"
);

// ---------------------------------------------------------------------------
// Las huellas congeladas.
//
// **Medidas el 08/10/2026 en el commit de partida, `fc2b205`, ANTES de tocar
// una línea del código**, con una copia sacada por `git archive` y compilada
// con su propio directorio de compilación. Coinciden con las de la ejecución
// por línea de órdenes archivada en `linea-base.sha256` de la evidencia.
//
// Si se midieran **después** del cambio, la prueba sería circular: congelaría
// lo que el código nuevo produce y no podría fallar por el defecto que dice
// cubrir. No se actualizan para que callen: un cambio aquí significa que el
// valor por omisión ha dejado de ser la calibración tal cual.
//
// **Excepción, 13/10/2026 (DT-034):** el manifiesto de evento y el bloque del
// tartamudeo clonal de la efectiva cambian de huella solo porque el texto de la
// fuente de la calibración ahora atribuye las cifras a Gymrek (2016) además de
// Raz (2019). Comprobado: con el texto antiguo restituido, los bytes dan las
// huellas de partida. Los FASTQ y la verdad no se mueven.
// ---------------------------------------------------------------------------

const FROZEN: [(&str, &str); 6] = [
    (
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz",
        "fe3024e5da040cef13f9372c1db3eb23e592f97e5283fbd66f7c831fb5960fe8",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz",
        "0df5d5c66a2d958dd623d176c1a77290bf883613cc8b564e11e0bac6cc201264",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R1.fastq.gz",
        "d5176776afce409db804d4ef6925261ac1d02a3169bfd6d912546ba41dd3a1be",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R2.fastq.gz",
        "a29364c401393696809b92c23ab1ca3089f478b0122baf0ebf7d63c3ccd67c5f",
    ),
    // Cambia a propósito el 30/10/2026 (entrega 3 del plan de detección): era
    // 0ee74f98b072b70e88ed344656e250b3a978def3ac97aa04e20f8a3de2e2d7ae. El TOML
    // de la evidencia lleva desde entonces la sección de microsatélites, y la
    // verdad cita el identificador de ejecución y la huella del TOML; quitados
    // esos dos, es la misma verdad que publica el binario de `dcae406` con el
    // TOML de antes (ver `tests/technical_catalog_contract.rs`).
    (
        "outputs/truth.jsonl",
        "fa398a3d952478bc4318fb37e3ac0c0246475cb91ed90c2b2c84e47d2ec11154",
    ),
    (
        "events/reads-main/manifest.json",
        "bc59706c711f5bdcabefddc8f39b00ac4ddcda3cf998c4f26aabfd198ec30576",
    ),
];

/// SHA-256 de los bytes canónicos —`serde_json::to_vec`, claves ordenadas— del
/// bloque `technical_event.artifacts.pcr.stutter` de la configuración
/// efectiva. La efectiva entera no se congela porque registra la ruta del
/// TOML, que aquí es un directorio temporal; el bloque no la lleva.
const FROZEN_EFFECTIVE_BLOCK: &str =
    "91c01cb95f5b63d920660c6bb843b7ccbdc08bbecb4dc9bfd233f516f6016539";

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-dtc-099-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
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
    fn read(&self, relative: &str) -> Vec<u8> {
        fs::read(self.destination.join(relative)).expect("artefacto publicado")
    }

    fn json(&self, relative: &str) -> Value {
        serde_json::from_slice(&self.read(relative)).expect("JSON")
    }

    fn effective(&self) -> Value {
        self.json("effective-config.json")
    }

    fn event(&self) -> Value {
        self.json("events/reads-main/manifest.json")
    }
}

/// El TOML del escenario con la tabla del tartamudeo clonal indicada, o sin
/// ella si `clonal` está vacío.
fn scenario(clonal: &str) -> String {
    if clonal.is_empty() {
        return SCENARIO_TOML.to_owned();
    }
    let anchor = "\n[technical_event.stutter]\n";
    assert_eq!(SCENARIO_TOML.matches(anchor).count(), 1, "fixture");
    SCENARIO_TOML.replacen(
        anchor,
        &format!("\n[technical_event.pcr.stutter]\n{clonal}{anchor}"),
        1,
    )
}

fn publish(toml: &str) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    fs::write(
        root.join("general.toml"),
        detection_section::with_detection_section(toml),
    )
    .expect("TOML");
    fs::write(root.join("reference.fa"), SCENARIO_FASTA).expect("FASTA");
    fs::write(root.join("population-snvs.tsv"), SCENARIO_SNVS).expect("matriz");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&root.join("general.toml")).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Published { root, destination }
}

fn effective_block(document: &Value) -> &Value {
    &document["technical_event"]["artifacts"]["pcr"]["stutter"]
}

fn event_block(document: &Value) -> &Value {
    &document["artifacts"]["pcr"]["stutter"]
}

/// **A. Sin el campo, el escenario amplificado con microsatélites publica los
/// bytes del commit de partida.**
///
/// **Debe fallar contra** un valor por omisión implementado como
/// `6,30 × 0,027` (0,1701 en vez de 0,17: cambian la probabilidad publicada y
/// los FASTQ), contra el multiplicador aplicado a todas las longitudes y
/// contra un campo publicado siempre (cambian el bloque y el manifiesto de
/// evento). Las huellas se midieron antes del cambio; ver la nota de arriba.
#[test]
fn the_amplified_scenario_without_the_fold_keeps_its_frozen_digests() {
    let published = publish(&scenario(""));
    for (relative, expected) in FROZEN {
        assert_eq!(
            digest(&published.read(relative)),
            expected,
            "{relative} ha cambiado de bytes sin declarar el multiplicador; NO se actualiza el \
             literal: el valor por omisión tiene que ser la calibración tal cual"
        );
    }
    let effective = published.effective();
    let block = effective_block(&effective);
    assert_eq!(
        digest(&serde_json::to_vec(block).expect("JSON")),
        FROZEN_EFFECTIVE_BLOCK,
        "el bloque del tartamudeo clonal de la configuración efectiva ha cambiado"
    );
    // Control positivo: el bloque congelado existe y es el del tartamudeo
    // clonal, no un `null` cuya huella también sería estable.
    assert_eq!(block["algorithm"], "clonal-stutter-v2");
    assert_eq!(block["per_cycle_slip"][0]["with_pcr_fraction"], 0.17);
}

/// **E. Con el multiplicador, los dos documentos lo publican, la fila del
/// homopolímero lleva la fracción resultante y nada más del bloque se mueve;
/// sin él, la clave no aparece en ninguno.**
///
/// **Debe fallar contra** un esquema olvidado en uno de los dos documentos,
/// contra publicarlo solo en la efectiva y contra publicarlo siempre. La
/// búsqueda negativa —sin el campo, la clave no está— lleva su control
/// positivo: con el campo sí está, en los dos documentos.
#[test]
fn the_fold_is_published_in_both_documents_only_when_declared() {
    let plain = publish(&scenario(""));
    let folded = publish(&scenario("homopolymer_pcr_fold = 2.5\n"));
    let (plain_effective, plain_event) = (plain.effective(), plain.event());
    let (folded_effective, folded_event) = (folded.effective(), folded.event());
    for (what, block) in [
        (
            "efectiva sin multiplicador",
            effective_block(&plain_effective),
        ),
        ("evento sin multiplicador", event_block(&plain_event)),
    ] {
        assert!(block.is_object(), "{what}: el bloque existe");
        assert!(
            block.get("homopolymer_pcr_fold").is_none(),
            "{what}: sin declarar, el campo no se publica"
        );
    }
    for (what, block, plain_block) in [
        (
            "efectiva",
            effective_block(&folded_effective),
            effective_block(&plain_effective),
        ),
        (
            "evento",
            event_block(&folded_event),
            event_block(&plain_event),
        ),
    ] {
        assert_eq!(block["homopolymer_pcr_fold"], 2.5, "{what}");
        let row = &block["per_cycle_slip"][0];
        assert_eq!(row["with_pcr_fraction"], 0.0675, "{what}");
        assert_eq!(row["without_pcr_fraction"], 0.027, "{what}");
        assert_eq!(row["source"], "derived", "{what}");
        let probability = row["probability"].as_f64().expect("probabilidad");
        assert!(
            (probability - 0.004_242_5).abs() < 1e-7,
            "{what}: la probabilidad publicada con 2,5 es 0,0042425: {probability}"
        );
        // Las filas de 2 a 6 y el resto del bloque, idénticos.
        for index in 1..6 {
            assert_eq!(
                block["per_cycle_slip"][index], plain_block["per_cycle_slip"][index],
                "{what}: la fila {index} no puede moverse"
            );
        }
        for key in [
            "algorithm",
            "cycles",
            "formula",
            "coarse",
            "coarse_reason",
            "calibration",
        ] {
            assert_eq!(block[key], plain_block[key], "{what}: {key}");
        }
    }
    assert_eq!(
        effective_block(&folded_effective),
        event_block(&folded_event),
        "los dos documentos publican el mismo bloque"
    );
    // La verdad no depende del tartamudeo; los FASTQ, sí. La cabecera de la
    // verdad lleva el identificador de ejecución y la huella de la
    // configuración, que cambian con los bytes del TOML; se comparan los
    // registros, que son lo que dice el genotipo.
    let records = |published: &Published| {
        String::from_utf8(published.read("outputs/truth.jsonl"))
            .expect("UTF-8")
            .lines()
            .skip(1)
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    let plain_records = records(&plain);
    assert!(plain_records.len() > 3, "la verdad tiene registros");
    assert_eq!(records(&folded), plain_records);
    assert_ne!(
        folded.read("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
        plain.read("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
        "control positivo: el multiplicador cambia lo que se secuencia"
    );
    // Los dos documentos con el campo pasan el validador externo de esquemas.
    assert_external_schema(
        &folded.destination.join("events/reads-main/manifest.json"),
        "manifest-event-v3.schema.json",
        true,
    );
    // **Y la efectiva y el manifiesto general también, directamente** (12/10/2026,
    // DT-040 y DT-041). Hasta ese día la efectiva publicaba
    // `technical_event.rng_catalog = "f4-technical-v5"` con tartamudeo clonal y
    // su esquema exige `f4-technical-v6`; aquí había un rodeo que comprobaba que
    // ese era el único motivo de rechazo. Con el catálogo corregido en el
    // productor el rodeo sobra: las dos efectivas —con y sin multiplicador— y
    // sus manifiestos generales pasan su esquema tal como se publican.
    for (what, published) in [
        ("con multiplicador", &folded),
        ("sin multiplicador", &plain),
    ] {
        assert_eq!(
            published.json("effective-config.json")["technical_event"]["rng_catalog"],
            "f4-technical-v6",
            "{what}"
        );
        assert_external_schema(
            &published.destination.join("effective-config.json"),
            "effective-config-v6.schema.json",
            true,
        );
        assert_external_schema(
            &published.destination.join("manifest.json"),
            "manifest-general-v6.schema.json",
            true,
        );
    }
    // Control negativo del esquema: la efectiva con el catálogo de antes del
    // arreglo, el v5, se rechaza; es el defecto de DT-040 visto desde fuera.
    let mut patched = folded_effective.clone();
    patched["technical_event"]["rng_catalog"] = json!("f4-technical-v5");
    let patched_path = folded.root.join("effective-patched.json");
    fs::write(
        &patched_path,
        serde_json::to_vec_pretty(&patched).expect("JSON"),
    )
    .expect("copia");
    assert_external_schema(&patched_path, "effective-config-v6.schema.json", false);
    patched["technical_event"]["rng_catalog"] = json!("f4-technical-v6");
    fs::write(
        &patched_path,
        serde_json::to_vec_pretty(&patched).expect("JSON"),
    )
    .expect("copia");
    assert_external_schema(&patched_path, "effective-config-v6.schema.json", true);
    // Control negativo del esquema: un multiplicador fuera de dominio en esa
    // misma copia se rechaza, de modo que el esquema mira el campo.
    patched["technical_event"]["artifacts"]["pcr"]["stutter"]["homopolymer_pcr_fold"] = json!(40.0);
    fs::write(
        &patched_path,
        serde_json::to_vec_pretty(&patched).expect("JSON"),
    )
    .expect("copia");
    assert_external_schema(&patched_path, "effective-config-v6.schema.json", false);
}

/// **El documento de entrada con el campo pasa su esquema, y con 40 no.**
///
/// El esquema declara el dominio para la interfaz y el validador externo; la
/// comprobación que manda es la de Rust (`tests/pcr_admission.rs`).
#[test]
fn the_input_schema_admits_the_fold_and_bounds_it() {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    for (literal, valid) in [
        ("2.5", true),
        ("37", true),
        ("1", true),
        ("40", false),
        ("0.5", false),
    ] {
        let toml = root.join(format!("general-{literal}.toml"));
        fs::write(
            &toml,
            scenario(&format!("homopolymer_pcr_fold = {literal}\n")),
        )
        .expect("TOML");
        let json = root.join(format!("general-{literal}.json"));
        let converted = Command::new("python3")
            .arg("-c")
            .arg(
                "import json, sys, tomllib; \
                 json.dump(tomllib.load(open(sys.argv[1], 'rb')), open(sys.argv[2], 'w'))",
            )
            .arg(&toml)
            .arg(&json)
            .status()
            .expect("python3");
        assert!(converted.success(), "conversión del TOML");
        assert_external_schema(&json, "config-general-v4.schema.json", valid);
    }
    let _ = fs::remove_dir_all(&root);
}

/// **Los topes de los tres esquemas son las constantes de Rust.**
///
/// Debe fallar contra un esquema con otro tope, o sin él, en cualquiera de los
/// tres documentos: se comprueban por separado.
#[test]
fn the_schema_bounds_of_the_fold_match_the_rust_constants() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for (schema, pointer) in [
        (
            "config-general-v4.schema.json",
            "/$defs/technical_event/properties/pcr/properties/stutter/properties/homopolymer_pcr_fold",
        ),
        (
            "effective-config-v6.schema.json",
            "/$defs/technical_artifacts/properties/pcr/oneOf/1/properties/stutter/properties/homopolymer_pcr_fold",
        ),
        (
            "manifest-event-v3.schema.json",
            "/$defs/technical_artifacts/properties/pcr/oneOf/1/properties/stutter/properties/homopolymer_pcr_fold",
        ),
    ] {
        let document: Value =
            serde_json::from_slice(&fs::read(root.join("schemas").join(schema)).expect("esquema"))
                .expect("JSON");
        let property = document
            .pointer(pointer)
            .unwrap_or_else(|| panic!("{schema} no declara el multiplicador en {pointer}"));
        assert_eq!(
            property["minimum"].as_f64(),
            Some(HOMOPOLYMER_PCR_FOLD_MINIMUM),
            "{schema}"
        );
        assert_eq!(
            property["maximum"].as_f64(),
            Some(HOMOPOLYMER_PCR_FOLD_MAXIMUM),
            "{schema}"
        );
        // Y el bloque que lo contiene sigue cerrado y sin exigirlo.
        let stutter = document
            .pointer(pointer.trim_end_matches("/properties/homopolymer_pcr_fold"))
            .expect("bloque");
        assert_eq!(stutter["additionalProperties"], false, "{schema}");
        assert!(
            !stutter["required"]
                .as_array()
                .is_some_and(|required| required.contains(&json!("homopolymer_pcr_fold"))),
            "{schema}: el campo es opcional"
        );
    }
}

/// **Con multiplicador 1, ningún resbalón de linaje en el homopolímero.**
///
/// La probabilidad por ciclo vale cero y `bernoulli-v1` no sortea. Control
/// positivo: el dinucleótido sí resbala, y sin multiplicador el homopolímero
/// también; sin él, un recorrido que no sorteara nada pasaría.
#[test]
fn a_unit_fold_leaves_the_homopolymer_without_lineage_slips() {
    let slips = |published: &Published| {
        published.event()["artifacts"]["pcr"]["lineage_slips_by_locus"]
            .as_array()
            .expect("recuentos")
            .iter()
            .map(|entry| {
                (
                    entry["locus_id"].as_str().expect("locus").to_owned(),
                    entry["slips"].as_u64().expect("recuento"),
                )
            })
            .collect::<std::collections::BTreeMap<_, _>>()
    };
    let unit = slips(&publish(&scenario("homopolymer_pcr_fold = 1\n")));
    assert_eq!(unit["STR-HOMO"], 0, "con 1 la PCR no añade tartamudeo");
    assert!(unit["STR-DI"] > 0, "el dinucleótido sigue resbalando");
    let plain = slips(&publish(&scenario("")));
    assert!(
        plain["STR-HOMO"] > 0,
        "sin multiplicador, el homopolímero resbala"
    );
}

// ---------------------------------------------------------------------------
// F. El juez.
// ---------------------------------------------------------------------------

/// Pasa el juez de los recursos publicados con la efectiva dada y el
/// manifiesto general dado.
fn judge(published: &Published, effective: &Value, manifest: &[u8]) -> Result<(), Diagnostic> {
    validate_published_resources(
        &published.destination,
        &published.destination.join("effective-config.json"),
        &serde_json::to_vec_pretty(effective).expect("JSON"),
        &published.destination.join("manifest.json"),
        manifest,
        None,
    )
}

fn detail(result: Result<(), Diagnostic>) -> String {
    let diagnostic = result.expect_err("el juez tenía que rechazar");
    diagnostic
        .context()
        .get("detail")
        .cloned()
        .unwrap_or_default()
}

/// Reescribe en disco el manifiesto de evento y **rehace sus huellas** en el
/// manifiesto general, que es lo que haría quien quisiera colarlo. Devuelve
/// los bytes del manifiesto general resultante.
fn forge_event_manifest(published: &Published, event: &Value) -> Vec<u8> {
    fn rehash(value: &mut Value, digest: &str, size: usize) {
        match value {
            Value::Object(map) => {
                if map.get("path").and_then(Value::as_str)
                    == Some("events/reads-main/manifest.json")
                {
                    map.insert(
                        "checksum".to_owned(),
                        json!({"algorithm": "sha256", "digest": digest}),
                    );
                    map.insert("size_bytes".to_owned(), json!(size));
                }
                for nested in map.values_mut() {
                    rehash(nested, digest, size);
                }
            }
            Value::Array(items) => {
                for item in items {
                    rehash(item, digest, size);
                }
            }
            _ => {}
        }
    }
    let bytes = serde_json::to_vec_pretty(event).expect("JSON");
    fs::write(
        published
            .destination
            .join("events/reads-main/manifest.json"),
        &bytes,
    )
    .expect("manifiesto de evento");
    let mut manifest = published.json("manifest.json");
    rehash(&mut manifest, &digest(&bytes), bytes.len());
    serde_json::to_vec_pretty(&manifest).expect("JSON")
}

const TOML_MISMATCH: &str = "technical_event efectivo no coincide exactamente con el TOML resuelto";
const EVENT_MISMATCH: &str =
    "el tartamudeo clonal del manifiesto técnico no es el de la configuración efectiva";

/// Cambia la fracción con PCR de una fila y recalcula su probabilidad con la
/// función del productor, de modo que la fila siga siendo coherente consigo
/// misma.
fn coherent_with_pcr(row: &mut Value, with_pcr: f64, cycles: u32) {
    let without = row["without_pcr_fraction"].as_f64().expect("sin PCR");
    row["with_pcr_fraction"] = json!(with_pcr);
    row["probability"] = json!(per_cycle_slip_probability(with_pcr, without, cycles));
}

/// **F. Sin multiplicador, la efectiva con fracciones cambiadas junto con su
/// probabilidad se rechaza.**
///
/// Antes de escribir nada se reprodujo, sobre el commit de partida, la
/// sospecha de la propuesta: que el juez solo comprobara que la probabilidad
/// sale de las fracciones publicadas. **En la efectiva no se reprodujo**: la
/// rechaza la comparación con lo que resuelve el TOML copiado, y es esa
/// comparación la que esta prueba fija. Con esa guarda borrada, las dos
/// manipulaciones pasan (se verificó mutándola).
#[test]
fn the_judge_rejects_coherently_altered_fractions_in_the_effective_config() {
    let published = publish(&scenario(""));
    let manifest = published.read("manifest.json");
    let original = published.effective();
    judge(&published, &original, &manifest).expect("control positivo: sin tocar se acepta");

    let mut homopolymer = original.clone();
    coherent_with_pcr(
        &mut homopolymer["technical_event"]["artifacts"]["pcr"]["stutter"]["per_cycle_slip"][0],
        0.0675,
        10,
    );
    assert_eq!(
        detail(judge(&published, &homopolymer, &manifest)),
        TOML_MISMATCH
    );

    let mut trinucleotide = original.clone();
    let row =
        &mut trinucleotide["technical_event"]["artifacts"]["pcr"]["stutter"]["per_cycle_slip"][2];
    row["without_pcr_fraction"] = json!(0.005);
    row["probability"] = json!(per_cycle_slip_probability(0.011, 0.005, 10));
    assert_eq!(
        detail(judge(&published, &trinucleotide, &manifest)),
        TOML_MISMATCH
    );
}

/// **F. Con multiplicador, cambiarlo, quitarlo, llevarlo fuera de dominio o
/// declararlo junto a una fila del homopolímero marcada como suya se rechaza.**
///
/// Las cuatro manipulaciones dejan la fila coherente consigo misma y pasan la
/// comprobación de la probabilidad; las rechaza la comparación con el TOML.
#[test]
fn the_judge_rejects_every_manipulation_of_a_declared_fold() {
    let published = publish(&scenario("homopolymer_pcr_fold = 2.5\n"));
    let manifest = published.read("manifest.json");
    let original = published.effective();
    judge(&published, &original, &manifest).expect("control positivo: sin tocar se acepta");
    type Mutation = fn(&mut Value);
    let cases: [(&str, Mutation); 4] = [
        ("otro multiplicador", |block| {
            block["homopolymer_pcr_fold"] = json!(3.0);
            coherent_with_pcr(&mut block["per_cycle_slip"][0], 0.081, 10);
        }),
        ("quitado, dejando 0,0675", |block| {
            block
                .as_object_mut()
                .expect("bloque")
                .remove("homopolymer_pcr_fold");
        }),
        ("fuera de dominio", |block| {
            block["homopolymer_pcr_fold"] = json!(40.0);
        }),
        ("con la fila marcada como declarada", |block| {
            block["per_cycle_slip"][0]["source"] = json!("declared");
        }),
    ];
    for (what, mutate) in cases {
        let mut altered = original.clone();
        mutate(&mut altered["technical_event"]["artifacts"]["pcr"]["stutter"]);
        assert_ne!(altered, original, "{what}: la mutación muta");
        assert_eq!(
            detail(judge(&published, &altered, &manifest)),
            TOML_MISMATCH,
            "{what}"
        );
    }
}

/// **F. El bloque del manifiesto de evento tiene que ser el de la efectiva.**
///
/// Este es el hueco que **sí** se reprodujo el 08/10/2026 sobre el commit de
/// partida: con la fracción con PCR del homopolímero y su probabilidad
/// cambiadas a la vez en el manifiesto de evento, y sus huellas rehechas en el
/// manifiesto general, el juez lo aceptaba; el contenido del bloque solo lo
/// comprobaba el esquema. **Debe fallar contra** el juez sin la igualdad entre
/// los dos bloques (se verificó borrándola) y contra una igualdad **parcial**
/// que solo comparara la tabla por ciclo y el multiplicador: por eso se
/// alteran también, uno a uno, los ciclos del bloque, la fuente de la
/// calibración, la fórmula y el motivo de la tosquedad (se verificó con esa
/// comparación parcial). Control positivo: el manifiesto general rehecho sin
/// cambiar nada se acepta, de modo que lo que rechaza es el bloque y no el
/// rehacer las huellas.
///
/// **Lo que esta prueba no cubre:** el resto de la sección de amplificación
/// del manifiesto de evento, fuera del bloque —`artifacts.pcr.cycles`, los
/// recuentos de linaje—, que el juez sigue sin atar a nada (DT-043).
#[test]
fn the_judge_rejects_an_event_manifest_whose_clonal_block_differs() {
    for clonal in ["", "homopolymer_pcr_fold = 2.5\n"] {
        let published = publish(&scenario(clonal));
        let effective = published.effective();
        let original = published.event();
        let untouched = forge_event_manifest(&published, &original);
        judge(&published, &effective, &untouched).expect("control positivo: rehecho sin cambios");

        let mut fraction = original.clone();
        coherent_with_pcr(
            &mut fraction["artifacts"]["pcr"]["stutter"]["per_cycle_slip"][0],
            0.05,
            10,
        );
        let manifest = forge_event_manifest(&published, &fraction);
        assert_eq!(
            detail(judge(&published, &effective, &manifest)),
            EVENT_MISMATCH,
            "fracción y probabilidad cambiadas a la vez ({clonal:?})"
        );

        let mut fold = original.clone();
        let block = fold["artifacts"]["pcr"]["stutter"]
            .as_object_mut()
            .expect("bloque");
        if block.remove("homopolymer_pcr_fold").is_none() {
            block.insert("homopolymer_pcr_fold".to_owned(), json!(6.3));
        }
        let manifest = forge_event_manifest(&published, &fold);
        assert_eq!(
            detail(judge(&published, &effective, &manifest)),
            EVENT_MISMATCH,
            "multiplicador quitado o añadido ({clonal:?})"
        );

        type Mutation = fn(&mut Value);
        let outside_the_table: [(&str, Mutation); 4] = [
            ("ciclos del bloque", |block| block["cycles"] = json!(12)),
            ("fuente de la calibración", |block| {
                block["calibration"]["source"] = json!("Gymrek (2016)");
            }),
            ("fórmula", |block| {
                block["formula"] = json!("r_pcr = r_with_pcr - r_without_pcr");
            }),
            ("motivo de la tosquedad", |block| {
                block["coarse_reason"] = json!("measured");
            }),
        ];
        for (what, mutate) in outside_the_table {
            let mut altered = original.clone();
            mutate(&mut altered["artifacts"]["pcr"]["stutter"]);
            assert_ne!(altered, original, "{what}: la mutación muta");
            let manifest = forge_event_manifest(&published, &altered);
            assert_eq!(
                detail(judge(&published, &effective, &manifest)),
                EVENT_MISMATCH,
                "{what} ({clonal:?})"
            );
        }
    }
}

fn external_schema(document: &Path, schema: &str) -> std::process::Output {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(document)
        .output()
        .expect("el validador externo se ejecuta")
}

fn assert_external_schema(document: &Path, schema: &str, valid: bool) {
    let output = external_schema(document, schema);
    assert_eq!(
        output.status.success(),
        valid,
        "{schema} sobre {document:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}
