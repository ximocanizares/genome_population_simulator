//! Cada rechazo del juez del manifiesto de evento, con su motivo exacto.
//!
//! DT-001, tanda 4. El juez (`validate_event_manifest_bytes`) tenía 147
//! comprobaciones de las que 117 no las disparaba ninguna prueba, y las que sí
//! se disparaban casi siempre bajo un `.is_err()` que no exige **cuál** de
//! ellas. Aquí cada caso parte de un manifiesto de evento **publicado de
//! verdad** por el programa (cuatro escenarios: sin amplificar con dos
//! bibliotecas, contaminación y salto de índice; con amplificación; con
//! amplificación y microsatélites; con tartamudeo sin amplificar), lo altera en
//! **un solo punto** y exige el motivo exacto del rechazo.
//!
//! # Cómo leer una tabla
//!
//! Cada tabla es una lista de [`Case`]: una etiqueta, la alteración y el
//! texto exacto que el juez tiene que dar en `detail`. El motivo exacto **es**
//! la discriminación: dos guardas con la misma causa tipada
//! (`GPS-FMT-006`) solo se distinguen por él. Además cada ejecución comprueba
//! que **la alteración alteró** el documento (si no, el caso no probaría nada)
//! y que el documento intacto se acepta (control positivo: sin él, un juez que
//! rechazara todo satisfaría la tabla entera). Los casos de «lado aceptado»
//! fijan el **lado aceptado** de las comparaciones (límites inclusivos y
//! dominios): un juez que endureciera de más los rechazaría.
//!
//! # Lo que este fichero no es
//!
//! Los documentos de evento **fallido** que aquí se juzgan se derivan del
//! completo alterándolo: el programa no publica hoy un manifiesto de evento
//! fallido por la vía del ensamblado (DT-036), de modo que el contrato de ese
//! estado se ejercita sobre documentos fabricados. Y el juez de los recursos
//! publicados (`src/results.rs`), que ata este manifiesto a la configuración
//! efectiva, queda fuera: sus huecos (DT-043) no se cubren ni se fijan aquí.
//!
//! # Las guardas inalcanzables
//!
//! Las que ningún manifiesto publicado y manipulado puede alcanzar están
//! anotadas **junto a la guarda** en `src/event_manifest.rs`, con el porqué.
//! Este fichero ejecuta el intento en `unreachable_guards_are_shadowed_by_other_checks`
//! y lo dice: esa prueba **no discrimina** esas guardas, solo muestra qué otra
//! comprobación gana la carrera.

use std::fs;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{
    EffectiveTechnicalEvent, build_effective_config, load_and_validate_config,
};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::error_quality::SequencingErrorsMode;
use genome_population_simulator::event_manifest::{
    EventManifest, EventManifestArtifact, EventManifestCoverage, EventManifestDemultiplexing,
    EventManifestErrorQuality, EventManifestFinalFastqSummary, EventManifestReference,
    event_manifest_document, validate_event_manifest_bytes,
};
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::fastq_validation::F4CoverageReport;
use genome_population_simulator::identity::{EventId, IndividualId, LibraryId, ReferenceId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::results::{ArtifactRecord, UnitState};
use genome_population_simulator::technical_identity::{
    TechnicalIdentityPlan, TechnicalMate, derive_technical_identity_plan,
};
use serde_json::{Value, json};

#[path = "support/detection_section.rs"]
mod detection_section;

static NONCE: AtomicU64 = AtomicU64::new(0);
const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
const DOCUMENT_PATH: &str = "events/reads-main/manifest.json";

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-event-guards-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

/// Publica una ejecución real y devuelve su manifiesto de evento.
fn publish(fasta: &str, matrix: &str, config: &str) -> Value {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    fs::write(root.join("reference.fa"), fasta).expect("FASTA");
    fs::write(root.join("population-snvs.tsv"), matrix).expect("matriz");
    let path = root.join("general.toml");
    fs::write(&path, detection_section::with_detection_section(config)).expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&path).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    let document = serde_json::from_slice(
        &fs::read(destination.join(DOCUMENT_PATH)).expect("manifiesto de evento"),
    )
    .expect("JSON");
    fs::remove_dir_all(root).expect("limpieza");
    document
}

// ---------------------------------------------------------------------------
// Los cuatro escenarios publicados
// ---------------------------------------------------------------------------

/// Dos bibliotecas, sin amplificar, con errores de secuenciación, un 50 % de
/// contaminación biológica de A con B y salto de índice a uno.
fn plain() -> Value {
    static CELL: OnceLock<Value> = OnceLock::new();
    CELL.get_or_init(|| {
        let config = format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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
length = 1000

[[initial_populations]]
population_id = "P1"
individual_count = 2
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A", "sample-B"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A", "sample-B"]
target_coverage_wgs = 1
contamination_min_percent = 50
contamination_max_percent = 50
contamination_recipient_ids = ["sample-A"]
index_hopping_rate = 1

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"

[[technical_event.libraries]]
individual_id = "sample-B"
i7 = "CGTACGTA"
i5 = "GCATGCAT"
"#
        );
        publish(
            &format!(">chr1\n{}\n", "A".repeat(1000)),
            "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t500\tA\tT\tP1\t0.5\n",
            &config,
        )
    })
    .clone()
}

/// Una biblioteca amplificada (30 % de duplicados, 10 ciclos) sin
/// microsatélites: la sección de amplificación sin linaje, versión 3 y
/// catálogo aleatorio v5.
fn pcr() -> Value {
    static CELL: OnceLock<Value> = OnceLock::new();
    CELL.get_or_init(|| {
        let config = format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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
length = 1000

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
sequencing_errors = "disabled"
library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.3
cycles = 10
efficiency = 0.9

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        );
        publish(
            &format!(">chr1\n{}\n", &"ACGTTGCAAGCT".repeat(84)[..1000]),
            "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tA\tT\tP1\t0.5\n",
            &config,
        )
    })
    .clone()
}

/// Una biblioteca amplificada con dos microsatélites y tartamudeo: la sección
/// de amplificación con su bloque de tartamudeo clonal y los tres recuentos de
/// linaje, más la sección de tartamudeo (versión 3, catálogo aleatorio v6).
fn clonal() -> Value {
    static CELL: OnceLock<Value> = OnceLock::new();
    CELL.get_or_init(|| {
        let short = "CA".repeat(10);
        let long = "CA".repeat(17);
        let head = format!("{}{short}GGGG{long}", "G".repeat(60));
        let bases = format!("{head}{}", "G".repeat(300 - head.len()));
        let config = format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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
length = 300

[[reference.repeat_loci]]
id = "STR-SHORT"
sequence = "chr1"
start = 61
end = 80
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{short}"
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-LONG"
sequence = "chr1"
start = 85
end = 118
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{long}"
frequency = 1.0

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 400
sequencing_errors = "disabled"
library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.3
cycles = 10
efficiency = 0.9

[technical_event.stutter]
enabled = true

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
rate = 0.3

[[technical_event.stutter.steps]]
units = -3
weight = 0.05

[[technical_event.stutter.steps]]
units = -2
weight = 0.15

[[technical_event.stutter.steps]]
units = -1
weight = 0.6

[[technical_event.stutter.steps]]
units = 1
weight = 0.2

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        );
        publish(
            &format!(">chr1\n{bases}\n"),
            "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tG\tT\tP1\t0.5\n",
            &config,
        )
    })
    .clone()
}

/// Tartamudeo sin amplificar, con dos microsatélites de dinucleótidos de
/// longitudes distintas (versión 3, catálogo aleatorio v4).
fn stutter_only() -> Value {
    static CELL: OnceLock<Value> = OnceLock::new();
    CELL.get_or_init(|| {
        let short = "CA".repeat(5);
        let long = "CA".repeat(17);
        let head = format!("{}{short}GGGG{long}", "G".repeat(60));
        let bases = format!("{head}{}", "G".repeat(300 - head.len()));
        let config = format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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
length = 300

[[reference.repeat_loci]]
id = "STR-SHORT"
sequence = "chr1"
start = 61
end = 70
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{short}"
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-LONG"
sequence = "chr1"
start = 75
end = 108
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{long}"
frequency = 1.0

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 20
sequencing_errors = "disabled"

[technical_event.stutter]
enabled = true

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
rate = 1.0

[[technical_event.stutter.steps]]
units = -3
weight = 0.05

[[technical_event.stutter.steps]]
units = -2
weight = 0.15

[[technical_event.stutter.steps]]
units = -1
weight = 0.6

[[technical_event.stutter.steps]]
units = 1
weight = 0.2

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        );
        publish(
            &format!(">chr1\n{bases}\n"),
            "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tG\tT\tP1\t0.5\n",
            &config,
        )
    })
    .clone()
}

// ---------------------------------------------------------------------------
// Infraestructura de las tablas
// ---------------------------------------------------------------------------

type Mutation = Box<dyn Fn(&mut Value)>;

/// Un rechazo esperado: qué se altera y qué motivo exacto tiene que dar.
struct Case {
    label: String,
    mutate: Mutation,
    detail: String,
}

fn case(
    label: impl Into<String>,
    mutate: impl Fn(&mut Value) + 'static,
    detail: impl Into<String>,
) -> Case {
    Case {
        label: label.into(),
        mutate: Box::new(mutate),
        detail: detail.into(),
    }
}

/// Un documento que, alterado, **tiene** que seguir aceptándose.
struct Accepted {
    label: String,
    mutate: Mutation,
}

fn accepted(label: impl Into<String>, mutate: impl Fn(&mut Value) + 'static) -> Accepted {
    Accepted {
        label: label.into(),
        mutate: Box::new(mutate),
    }
}

fn at<'a>(document: &'a mut Value, pointer: &str) -> &'a mut Value {
    document
        .pointer_mut(pointer)
        .unwrap_or_else(|| panic!("el documento no tiene {pointer}"))
}

/// Pone un valor; si la clave no existe en su objeto padre, la crea.
fn set(document: &mut Value, pointer: &str, value: Value) {
    if document.pointer(pointer).is_some() {
        *at(document, pointer) = value;
        return;
    }
    let (parent, key) = pointer.rsplit_once('/').expect("puntero con padre");
    at(document, parent)
        .as_object_mut()
        .unwrap_or_else(|| panic!("{parent} no es objeto"))
        .insert(key.to_owned(), value);
}

fn del(document: &mut Value, pointer: &str) {
    let (parent, key) = pointer.rsplit_once('/').expect("puntero con padre");
    match at(document, parent) {
        Value::Object(map) => {
            map.remove(key).unwrap_or_else(|| panic!("sin {pointer}"));
        }
        Value::Array(list) => {
            list.remove(key.parse::<usize>().expect("índice"));
        }
        other => panic!("{parent} no es contenedor: {other}"),
    }
}

fn push(document: &mut Value, pointer: &str, value: Value) {
    at(document, pointer)
        .as_array_mut()
        .unwrap_or_else(|| panic!("{pointer} no es lista"))
        .push(value);
}

fn text(document: &Value, pointer: &str) -> String {
    document
        .pointer(pointer)
        .and_then(Value::as_str)
        .unwrap_or_else(|| panic!("{pointer} no es texto"))
        .to_owned()
}

fn number(document: &Value, pointer: &str) -> u64 {
    document
        .pointer(pointer)
        .and_then(Value::as_u64)
        .unwrap_or_else(|| panic!("{pointer} no es u64"))
}

fn verdict(document: &Value) -> Result<(), Diagnostic> {
    validate_event_manifest_bytes(
        Path::new(DOCUMENT_PATH),
        &serde_json::to_vec(document).expect("JSON"),
    )
}

/// Ejecuta el juez sin dejar que un pánico (resta que desborda en un mutante)
/// tire la tabla entera: un pánico es un fallo del caso, no del programa.
fn judged(document: &Value) -> Result<Result<(), Diagnostic>, String> {
    catch_unwind(AssertUnwindSafe(|| verdict(document))).map_err(|payload| {
        payload
            .downcast_ref::<&str>()
            .map(|message| (*message).to_owned())
            .or_else(|| payload.downcast_ref::<String>().cloned())
            .unwrap_or_else(|| "pánico".to_owned())
    })
}

fn detail_of(error: &Diagnostic) -> String {
    assert_eq!(error.code(), "GPS-FMT-006", "{error:?}");
    assert_eq!(
        error.context().get("path").map(String::as_str),
        Some(DOCUMENT_PATH),
        "el diagnóstico nombra el fichero"
    );
    error
        .context()
        .get("detail")
        .cloned()
        .expect("el diagnóstico lleva detalle")
}

/// Los fallos de una prueba. Cada uno se **imprime al anotarlo**: una prueba
/// que luego entrara en pánico (un mutante puede romper lo que la prueba da por
/// supuesto) no perdería los que ya había visto.
#[derive(Default)]
struct Failures(Vec<String>);

impl Failures {
    fn push(&mut self, message: String) {
        eprintln!("{message}");
        self.0.push(message);
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }
}

/// Acumula los fallos de una tabla; `finish` los vuelca y falla una sola vez.
#[derive(Default)]
struct Ledger {
    failures: Failures,
}

impl Ledger {
    fn control(&mut self, name: &str, base: &Value) {
        match judged(base) {
            Ok(Ok(())) => {}
            other => self.failures.push(format!(
                "CASE-FAIL {name}/control: el documento intacto tiene que aceptarse: {other:?}"
            )),
        }
    }

    fn rejected(&mut self, name: &str, base: &Value, cases: Vec<Case>) {
        self.control(name, base);
        for case in cases {
            let mut document = base.clone();
            (case.mutate)(&mut document);
            let label = format!("{name}/{}", case.label);
            if &document == base {
                self.failures
                    .push(format!("CASE-FAIL {label}: la alteración no alteró nada"));
                continue;
            }
            match judged(&document) {
                Err(message) => self.failures.push(format!(
                    "CASE-FAIL {label}: el juez entró en pánico: {message} (esperado {:?})",
                    case.detail
                )),
                Ok(Ok(())) => self.failures.push(format!(
                    "CASE-FAIL {label}: se aceptó, y tenía que rechazarse (esperado {:?})",
                    case.detail
                )),
                Ok(Err(error)) => {
                    let detail = detail_of(&error);
                    if detail != case.detail {
                        self.failures.push(format!(
                            "CASE-FAIL {label}: motivo {detail:?}, esperado {:?}",
                            case.detail
                        ));
                    }
                }
            }
        }
    }

    fn accepted(&mut self, name: &str, base: &Value, cases: Vec<Accepted>) {
        self.control(name, base);
        for case in cases {
            let mut document = base.clone();
            (case.mutate)(&mut document);
            let label = format!("{name}/{}", case.label);
            if &document == base {
                self.failures
                    .push(format!("CASE-FAIL {label}: la alteración no alteró nada"));
                continue;
            }
            match judged(&document) {
                Ok(Ok(())) => {}
                other => self.failures.push(format!(
                    "CASE-FAIL {label}: tenía que aceptarse y dio {:?}",
                    other.map(|result| result.map_err(|error| detail_of(&error)))
                )),
            }
        }
    }

    fn finish(self) {
        assert!(
            self.failures.is_empty(),
            "{} casos fallan; ver arriba",
            self.failures.len()
        );
    }
}

// ---------------------------------------------------------------------------
// Alteraciones compartidas por varias tablas
// ---------------------------------------------------------------------------

/// Las dos bibliotecas del escenario, en el orden de sus identificadores.
fn libraries(document: &Value) -> (String, String) {
    (
        text(document, "/selected_individuals/0/library_id"),
        text(document, "/selected_individuals/1/library_id"),
    )
}

/// Un evento sin reasignaciones de índice (para poder quitar parejas finales).
fn neutral_demultiplexing(document: &mut Value) {
    set(document, "/demultiplexing/reassigned_pairs", json!(0));
    set(document, "/demultiplexing/reassigned_by_sample", json!([]));
}

/// Recompone los resúmenes FASTQ y las coberturas para una longitud de lectura
/// y unos recuentos de parejas, todo propio y sin N: el documento queda
/// coherente consigo mismo y cambia solo lo que se pide.
fn reshape_summary(document: &mut Value, read_length: u64, pairs: &[u64]) {
    set(document, "/profile/read_length", json!(read_length));
    for (index, pair_count) in pairs.iter().enumerate() {
        let bases = 2 * pair_count * read_length;
        let summary = format!("/final_fastq_summary/{index}");
        for (field, value) in [
            ("pair_count", *pair_count),
            ("read_count", 2 * pair_count),
            ("base_count", bases),
            ("own_pair_count", *pair_count),
            ("contaminant_pair_count", 0),
            ("own_base_count", bases),
            ("contaminant_base_count", 0),
            ("n_base_count", 0),
        ] {
            set(document, &format!("{summary}/{field}"), json!(value));
        }
        let coverage = format!("/coverage/{index}");
        let denominator = number(document, &format!("{coverage}/denominator_bases"));
        set(document, &format!("{coverage}/n_bases"), json!(0));
        set(
            document,
            &format!("{coverage}/final_x/own_x"),
            json!(bases as f64 / denominator as f64),
        );
        set(
            document,
            &format!("{coverage}/final_x/contaminant_x"),
            json!(0.0),
        );
    }
}

/// Añade una tercera biblioteca (`sample-C`) a un documento de dos.
fn with_third_library(document: &mut Value) {
    let mut selected = document["selected_individuals"][1].clone();
    selected["individual_id"] = json!("sample-C");
    selected["library_id"] = json!("LIB_C");
    selected["read_group_id"] = json!("RG_C");
    push(document, "/selected_individuals", selected);
    let mut coverage = document["coverage"][1].clone();
    coverage["individual_id"] = json!("sample-C");
    coverage["library_id"] = json!("LIB_C");
    push(document, "/coverage", coverage);
    let mut summary = document["final_fastq_summary"][1].clone();
    summary["individual_id"] = json!("sample-C");
    summary["library_id"] = json!("LIB_C");
    push(document, "/final_fastq_summary", summary);
    for (mate, source) in [(1, 2), (2, 3)] {
        let mut output = document["outputs"][source].clone();
        let id = format!("fastq-LIB_C-r{mate}");
        output["id"] = json!(id);
        output["individual_id"] = json!("sample-C");
        output["library_id"] = json!("LIB_C");
        output["artifact"]["id"] = json!(id);
        output["artifact"]["path"] = json!(format!(
            "outputs/fastq/reads-main/sample-C_R{mate}.fastq.gz"
        ));
        push(document, "/outputs", output);
    }
}

fn format_diagnostic() -> Value {
    json!({
        "code": "GPS-FMT-006",
        "message": "Un documento JSON producido no satisface su esquema publicado.",
    })
}

/// El mismo evento con un FASTQ fallido: sin coberturas ni resúmenes, con el
/// diagnóstico que lo explica y sin reasignaciones (no hay parejas finales).
fn failed_of(complete: &Value) -> Value {
    let mut failed = complete.clone();
    failed["state"] = json!("failed");
    failed["outputs"][0]["state"] = json!("failed");
    failed["outputs"][0]["artifact"] = Value::Null;
    failed["outputs"][0]["diagnostic_code"] = json!("GPS-FMT-006");
    failed["coverage"] = json!([]);
    failed["final_fastq_summary"] = json!([]);
    failed["diagnostics"] = json!([format_diagnostic()]);
    neutral_demultiplexing(&mut failed);
    failed
}

/// El mismo evento con descartes por UDI no reconocido, que el programa no
/// produce hoy pero que el contrato admite: el juez tiene que aceptarlo.
fn with_discards(complete: &Value) -> Value {
    let (first, second) = libraries(complete);
    let mut document = complete.clone();
    let entries = json!([
        {"source_library_id": first, "cause": "unrecognized_udi", "pair_count": 2},
        {"source_library_id": second, "cause": "unrecognized_udi", "pair_count": 1},
    ]);
    set(&mut document, "/demultiplexing/discarded_pairs", json!(3));
    set(&mut document, "/demultiplexing/unassigned_pairs", json!(3));
    set(
        &mut document,
        "/demultiplexing/discarded_by_sample",
        entries.clone(),
    );
    set(
        &mut document,
        "/demultiplexing/unassigned_by_sample",
        entries,
    );
    document
}

/// Alteraciones de forma: un objeto que deja de serlo o una lista que deja de
/// serlo, con el contexto que el juez nombra.
fn shape_cases(pointers: &[(&str, &str, bool)]) -> Vec<Case> {
    pointers
        .iter()
        .map(|(pointer, context, is_object)| {
            let kind = if *is_object { "objeto" } else { "lista" };
            let pointer = (*pointer).to_owned();
            case(
                format!("{pointer} deja de ser {kind}"),
                move |document| set(document, &pointer, json!("x")),
                format!("{context} no es {kind}"),
            )
        })
        .collect()
}

/// Alteraciones de campos cerrados: una clave de más y una de menos.
fn closed_cases(pointers: &[(&str, &str)]) -> Vec<Case> {
    let mut cases = Vec::new();
    for (pointer, context) in pointers {
        let (extra, missing, renamed) = (
            (*pointer).to_owned(),
            (*pointer).to_owned(),
            (*pointer).to_owned(),
        );
        let detail = format!("campos inválidos en {context}");
        cases.push(case(
            format!("{pointer}: una clave de más"),
            move |document| {
                at(document, &extra)
                    .as_object_mut()
                    .expect("objeto")
                    .insert("zz_extra".to_owned(), json!(1));
            },
            detail.clone(),
        ));
        cases.push(case(
            format!("{pointer}: una clave de menos"),
            move |document| {
                let first = at(document, &missing)
                    .as_object()
                    .expect("objeto")
                    .keys()
                    .next()
                    .expect("alguna clave")
                    .clone();
                at(document, &missing)
                    .as_object_mut()
                    .expect("objeto")
                    .remove(&first);
            },
            detail.clone(),
        ));
        // Mismo número de claves, una con otro nombre: sin esta, la mitad
        // «todas las claves pertenecen a la lista» de `validate_closed` no
        // la distingue nadie (la longitud ya caza las claves de más y de menos).
        cases.push(case(
            format!("{pointer}: una clave con otro nombre"),
            move |document| {
                let object = at(document, &renamed).as_object_mut().expect("objeto");
                let first = object.keys().next().expect("alguna clave").clone();
                let value = object.remove(&first).expect("valor");
                object.insert("zz_renamed".to_owned(), value);
            },
            detail,
        ));
    }
    cases
}

// ---------------------------------------------------------------------------
// Forma: objetos y listas que dejan de serlo, y campos cerrados
// ---------------------------------------------------------------------------

/// **Cada objeto y cada lista del manifiesto se comprueba por su contexto.**
///
/// `required_object`, `required_array` y `validate_closed` son ayudantes
/// compartidos: la puerta de cobertura los da por disparados con un solo uso,
/// y quitar la comprobación de **un** llamante no lo notaría ninguna prueba que
/// no nombre ese contexto. Por eso cada llamante tiene su caso, con el
/// contexto exacto en el motivo.
#[test]
fn every_object_list_and_closed_object_is_judged_by_its_own_context() {
    let mut ledger = Ledger::default();
    let plain = plain();
    ledger.rejected(
        "forma/plain",
        &plain,
        [
            shape_cases(&[
                ("", "raíz", true),
                ("/id_catalog", "id_catalog", true),
                ("/id_catalog/derivations", "id_catalog.derivations", true),
                ("/id_catalog/derivations/v1", "v1", false),
                ("/id_catalog/derivations/v2", "v2", false),
                ("/profile", "profile", true),
                ("/profile/quality", "profile.quality", true),
                ("/selected_individuals", "selected_individuals", false),
                ("/selected_individuals/0", "selected_individuals[]", true),
                ("/reference", "reference", true),
                ("/reference/fasta", "artifact", true),
                ("/reference/fasta/checksum", "checksum", true),
                ("/coverage", "coverage", false),
                ("/coverage/0", "coverage[]", true),
                ("/coverage/0/final_x", "coverage[].final_x", true),
                ("/final_fastq_summary", "final_fastq_summary", false),
                ("/final_fastq_summary/0", "final_fastq_summary[]", true),
                ("/error_quality", "error_quality", true),
                ("/error_quality/r1", "error_quality mate", true),
                ("/demultiplexing", "demultiplexing", true),
                ("/artifacts", "artifacts", true),
                ("/artifacts/chimera", "artifacts.chimera", true),
                ("/artifacts/pcr", "artifacts.pcr", true),
                (
                    "/artifacts/optical_duplicates",
                    "artifacts.optical_duplicates",
                    true,
                ),
                ("/artifacts/gc", "artifacts.gc", true),
                ("/artifacts/adapters", "artifacts.adapters", true),
                (
                    "/artifacts/biological_contamination",
                    "artifacts.biological_contamination",
                    false,
                ),
                (
                    "/artifacts/biological_contamination/0",
                    "artifacts.biological_contamination[]",
                    true,
                ),
                (
                    "/artifacts/biological_contamination/0/percent_range",
                    "percent_range",
                    true,
                ),
                (
                    "/artifacts/biological_contamination/0/subseeds",
                    "subseeds",
                    true,
                ),
                (
                    "/artifacts/biological_contamination/0/donor_contributions",
                    "donor_contributions",
                    false,
                ),
                (
                    "/artifacts/biological_contamination/0/donor_contributions/0",
                    "donor_contributions[]",
                    true,
                ),
                ("/outputs", "outputs", false),
                ("/outputs/0", "outputs[]", true),
                ("/outputs/0/artifact/checksum", "checksum", true),
            ]),
            closed_cases(&[
                ("", "raíz"),
                ("/id_catalog", "id_catalog"),
                ("/id_catalog/derivations", "id_catalog.derivations"),
                ("/profile", "profile"),
                ("/profile/quality", "profile.quality"),
                ("/selected_individuals/0", "selected_individuals[]"),
                ("/reference", "reference"),
                ("/reference/fasta", "artifact"),
                ("/reference/fasta/checksum", "checksum"),
                ("/coverage/0", "coverage[]"),
                ("/coverage/0/final_x", "coverage[].final_x"),
                ("/final_fastq_summary/0", "final_fastq_summary[]"),
                ("/error_quality", "error_quality"),
                ("/error_quality/r1", "error_quality mate"),
                ("/demultiplexing", "demultiplexing"),
                (
                    "/demultiplexing/reassigned_by_sample/0",
                    "reassigned_by_sample",
                ),
                ("/artifacts", "artifacts"),
                ("/artifacts/chimera", "artifacts.chimera"),
                ("/artifacts/pcr", "artifacts.pcr"),
                (
                    "/artifacts/optical_duplicates",
                    "artifacts.optical_duplicates",
                ),
                ("/artifacts/gc", "artifacts.gc"),
                ("/artifacts/adapters", "artifacts.adapters"),
                (
                    "/artifacts/biological_contamination/0",
                    "artifacts.biological_contamination[]",
                ),
                (
                    "/artifacts/biological_contamination/0/percent_range",
                    "percent_range",
                ),
                ("/artifacts/biological_contamination/0/subseeds", "subseeds"),
                (
                    "/artifacts/biological_contamination/0/donor_contributions/0",
                    "donor_contributions[]",
                ),
                ("/outputs/0", "outputs[]"),
                ("/outputs/0/artifact", "artifact"),
                ("/outputs/0/artifact/checksum", "checksum"),
            ]),
        ]
        .concat(),
    );
    let failed = failed_of(&plain);
    ledger.rejected(
        "forma/failed",
        &failed,
        [
            shape_cases(&[
                ("/diagnostics", "diagnostics", false),
                ("/diagnostics/0", "diagnostics[]", true),
            ]),
            closed_cases(&[("/diagnostics/0", "diagnostics[]")]),
        ]
        .concat(),
    );
    let pcr = pcr();
    ledger.rejected(
        "forma/pcr",
        &pcr,
        [
            shape_cases(&[
                (
                    "/artifacts/pcr/duplicate_fraction_inversion",
                    "artifacts.pcr.duplicate_fraction_inversion",
                    true,
                ),
                (
                    "/artifacts/pcr/resource_bound",
                    "artifacts.pcr.resource_bound",
                    true,
                ),
                (
                    "/artifacts/pcr/family_size_histogram",
                    "artifacts.pcr.family_size_histogram",
                    false,
                ),
                (
                    "/artifacts/pcr/family_size_histogram/0",
                    "artifacts.pcr.family_size_histogram[]",
                    true,
                ),
            ]),
            closed_cases(&[
                ("/artifacts/pcr", "artifacts.pcr"),
                (
                    "/artifacts/pcr/duplicate_fraction_inversion",
                    "artifacts.pcr.duplicate_fraction_inversion",
                ),
                (
                    "/artifacts/pcr/resource_bound",
                    "artifacts.pcr.resource_bound",
                ),
                (
                    "/artifacts/pcr/family_size_histogram/0",
                    "artifacts.pcr.family_size_histogram[]",
                ),
            ]),
        ]
        .concat(),
    );
    let clonal = clonal();
    let mut lineage_shapes = Vec::new();
    let mut lineage_closed = Vec::new();
    for field in [
        "lineage_slips_by_locus",
        "lineage_inherited_slips_by_locus",
        "lineage_applied_steps_by_locus",
    ] {
        let leaked: &'static str = Box::leak(format!("/artifacts/pcr/{field}").into_boxed_str());
        let context: &'static str = Box::leak(format!("artifacts.pcr.{field}").into_boxed_str());
        let entry: &'static str = Box::leak(format!("{leaked}/0").into_boxed_str());
        let entry_context: &'static str = Box::leak(format!("{context}[]").into_boxed_str());
        lineage_shapes.push((leaked, context, false));
        lineage_shapes.push((entry, entry_context, true));
        lineage_closed.push((entry, entry_context));
    }
    ledger.rejected(
        "forma/clonal",
        &clonal,
        [shape_cases(&lineage_shapes), closed_cases(&lineage_closed)].concat(),
    );
    let stutter = stutter_only();
    let dependence = "/artifacts/stutter/tract_length_dependence";
    let by_motif = "artifacts.stutter.tract_length_dependence.by_motif_length";
    ledger.rejected(
        "forma/stutter",
        &stutter,
        [
            shape_cases(&[
                ("/artifacts/stutter", "artifacts.stutter", true),
                (
                    "/artifacts/stutter/run_length_coefficient",
                    "artifacts.stutter.run_length_coefficient",
                    true,
                ),
                (
                    "/artifacts/stutter/calibration",
                    "artifacts.stutter.calibration",
                    true,
                ),
                (
                    dependence,
                    "artifacts.stutter.tract_length_dependence",
                    true,
                ),
                (
                    "/artifacts/stutter/tract_length_dependence/by_motif_length",
                    by_motif,
                    false,
                ),
                (
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0",
                    by_motif,
                    true,
                ),
                (
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/factors",
                    "artifacts.stutter.tract_length_dependence.by_motif_length.factors",
                    false,
                ),
                (
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/factors/0",
                    "artifacts.stutter.tract_length_dependence.by_motif_length.factors",
                    true,
                ),
                (
                    "/artifacts/stutter/step_weights",
                    "artifacts.stutter.step_weights",
                    false,
                ),
                (
                    "/artifacts/stutter/step_weights/0",
                    "artifacts.stutter.step_weights",
                    true,
                ),
                (
                    "/artifacts/stutter/steps_by_locus",
                    "artifacts.stutter.steps_by_locus",
                    false,
                ),
                (
                    "/artifacts/stutter/steps_by_locus/0",
                    "artifacts.stutter.steps_by_locus",
                    true,
                ),
                (
                    "/artifacts/stutter/steps_by_locus/0/steps",
                    "artifacts.stutter.steps_by_locus.steps",
                    false,
                ),
                (
                    "/artifacts/stutter/steps_by_locus/0/steps/0",
                    "artifacts.stutter.steps_by_locus.steps",
                    true,
                ),
            ]),
            closed_cases(&[
                ("/artifacts/stutter", "artifacts.stutter"),
                (
                    "/artifacts/stutter/run_length_coefficient",
                    "artifacts.stutter.run_length_coefficient",
                ),
                (
                    "/artifacts/stutter/calibration",
                    "artifacts.stutter.calibration",
                ),
                (dependence, "artifacts.stutter.tract_length_dependence"),
                (
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0",
                    by_motif,
                ),
                (
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/factors/0",
                    "artifacts.stutter.tract_length_dependence.by_motif_length.factors",
                ),
                (
                    "/artifacts/stutter/step_weights/0",
                    "artifacts.stutter.step_weights",
                ),
                (
                    "/artifacts/stutter/steps_by_locus/0",
                    "artifacts.stutter.steps_by_locus",
                ),
                (
                    "/artifacts/stutter/steps_by_locus/0/steps/0",
                    "artifacts.stutter.steps_by_locus.steps",
                ),
            ]),
        ]
        .concat(),
    );
    ledger.finish();
}

/// Une varias tablas de casos (`Case` no se puede clonar, así que no sirve
/// `concat` de la biblioteca estándar).
trait JoinCases {
    fn concat(self) -> Vec<Case>;
}

impl<const N: usize> JoinCases for [Vec<Case>; N] {
    fn concat(self) -> Vec<Case> {
        self.into_iter().flatten().collect()
    }
}

// ---------------------------------------------------------------------------
// Raíz, versión, catálogos, perfil y referencia
// ---------------------------------------------------------------------------

fn swap_items(document: &mut Value, pointer: &str, first: usize, second: usize) {
    at(document, pointer)
        .as_array_mut()
        .expect("lista")
        .swap(first, second);
}

/// **La raíz, la versión que va con las etapas, los catálogos, el perfil y la
/// referencia.** Cada literal que `require_text` compara tiene su caso:
/// quitar la comprobación de **un** llamante no la cazaría ninguna otra.
#[test]
fn the_root_the_version_the_catalogs_the_profile_and_the_reference_are_judged_exactly() {
    let mut ledger = Ledger::default();
    let plain = plain();
    ledger.rejected(
        "raíz/plain",
        &plain,
        vec![
            case(
                "tipo de manifiesto",
                |d| set(d, "/manifest_kind", json!("otro")),
                "manifest_kind incompatible",
            ),
            case(
                "tipo de evento",
                |d| set(d, "/event_type", json!("otro")),
                "event_type incompatible",
            ),
            case(
                "tecnología",
                |d| set(d, "/technology", json!("otra")),
                "technology incompatible",
            ),
            case(
                "versión desconocida",
                |d| set(d, "/schema_version", json!("9")),
                "schema_version debe ser 2 o 3, se encontró 9",
            ),
            case(
                "versión que no es texto",
                |d| set(d, "/schema_version", json!(2)),
                "schema_version debe ser texto no vacío",
            ),
            case(
                "versión vacía",
                |d| set(d, "/schema_version", json!("")),
                "schema_version debe ser texto no vacío",
            ),
            case(
                "v3 sin ninguna etapa nueva",
                |d| set(d, "/schema_version", json!("3")),
                "la versión del documento no corresponde a las etapas que declara",
            ),
            case(
                "v2 que declara amplificación",
                |d| set(d, "/artifacts/pcr/applicable", json!(true)),
                "la versión del documento no corresponde a las etapas que declara",
            ),
            case(
                "v2 con sección de tartamudeo",
                |d| set(d, "/artifacts/stutter", Value::Null),
                "la versión del documento no corresponde a las etapas que declara",
            ),
            case(
                "catálogo aleatorio de otra etapa",
                |d| set(d, "/rng_catalog", json!("f4-technical-v4")),
                "rng_catalog incompatible",
            ),
            case(
                "catálogo aleatorio de la amplificación",
                |d| set(d, "/rng_catalog", json!("f4-technical-v5")),
                "rng_catalog incompatible",
            ),
            case(
                "catálogo aleatorio vacío",
                |d| set(d, "/rng_catalog", json!("")),
                "rng_catalog debe ser texto no vacío",
            ),
            case(
                "identificador de evento inválido",
                |d| set(d, "/event_id", json!("no valido!")),
                "ID event inválido",
            ),
            case(
                "identificador de evento vacío",
                |d| set(d, "/event_id", json!("")),
                "event_id debe ser texto no vacío",
            ),
            case(
                "identificador de evento que no coincide con las rutas",
                |d| set(d, "/event_id", json!("reads-otro")),
                "artefacto FASTQ completo incoherente",
            ),
            case(
                "versión del catálogo de identificadores",
                |d| set(d, "/id_catalog/version", json!("9")),
                "version incompatible",
            ),
            case(
                "derivaciones v1 en otro orden",
                |d| swap_items(d, "/id_catalog/derivations/v1", 0, 1),
                "id_catalog incompatible",
            ),
            case(
                "derivaciones v1 con una de más",
                |d| push(d, "/id_catalog/derivations/v1", json!("extra")),
                "id_catalog incompatible",
            ),
            case(
                "derivaciones v2 en otro orden",
                |d| swap_items(d, "/id_catalog/derivations/v2", 0, 1),
                "id_catalog incompatible",
            ),
            case(
                "derivaciones v2 con una de menos",
                |d| del(d, "/id_catalog/derivations/v2/1"),
                "id_catalog incompatible",
            ),
            case(
                "estado inválido",
                |d| set(d, "/state", json!("quizá")),
                "estado inválido",
            ),
            case(
                "identificador de tanda inválido",
                |d| set(d, "/technical_batch_id", json!("no valido!")),
                "ID batch inválido",
            ),
            case(
                "subsemilla del evento en mayúsculas",
                |d| {
                    let upper = text(d, "/event_subseed").to_uppercase();
                    set(d, "/event_subseed", json!(upper));
                },
                "event_subseed debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "subsemilla del evento larga",
                |d| {
                    let long = format!("{}a", text(d, "/event_subseed"));
                    set(d, "/event_subseed", json!(long));
                },
                "event_subseed debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "huella de la referencia larga",
                |d| set(d, "/reference/fasta/checksum/digest", json!("a".repeat(65))),
                "checksum.digest debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "subsemilla del evento corta",
                |d| {
                    let short = text(d, "/event_subseed")[1..].to_owned();
                    set(d, "/event_subseed", json!(short));
                },
                "event_subseed debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "subsemilla del evento con un carácter que no es hexadecimal",
                |d| {
                    let mut seed = text(d, "/event_subseed");
                    seed.replace_range(0..1, "g");
                    set(d, "/event_subseed", json!(seed));
                },
                "event_subseed debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "perfil desconocido",
                |d| set(d, "/profile/id", json!("otro-perfil")),
                "id incompatible",
            ),
            case(
                "identificador del perfil con amplificación y recurso sin ella",
                |d| set(d, "/profile/id", json!("gps-hiseqx-pcr-pe150-v1")),
                "resource_path incompatible",
            ),
            case(
                "ruta del recurso del perfil",
                |d| set(d, "/profile/resource_path", json!("resources/f4/otro.json")),
                "resource_path incompatible",
            ),
            case(
                "huella del recurso del perfil",
                |d| set(d, "/profile/sha256", json!("0".repeat(64))),
                "sha256 incompatible",
            ),
            case(
                "preparación del perfil",
                |d| set(d, "/profile/preparation", json!("pcr")),
                "preparation incompatible",
            ),
            case(
                "longitud de lectura cero",
                |d| set(d, "/profile/read_length", json!(0)),
                "longitud de lectura de perfil inválida",
            ),
            case(
                "longitud de lectura por encima de 150",
                |d| set(d, "/profile/read_length", json!(151)),
                "longitud de lectura de perfil inválida",
            ),
            case(
                "longitud de lectura que es texto",
                |d| set(d, "/profile/read_length", json!("150")),
                "longitud de lectura de perfil inválida",
            ),
            case(
                "longitud de lectura nula",
                |d| set(d, "/profile/read_length", Value::Null),
                "longitud de lectura de perfil inválida",
            ),
            case(
                "calidad de las bases ACGT",
                |d| set(d, "/profile/quality/acgt_phred", json!(39)),
                "calidad de perfil inválida",
            ),
            case(
                "calidad de las bases N",
                |d| set(d, "/profile/quality/n_phred", json!(1)),
                "calidad de perfil inválida",
            ),
            case(
                "calidad ACGT que es texto",
                |d| set(d, "/profile/quality/acgt_phred", json!("40")),
                "calidad de perfil inválida",
            ),
            case(
                "identificador de la referencia inválido",
                |d| set(d, "/reference/reference_id", json!("no valido!")),
                "ID reference inválido",
            ),
            case(
                "la referencia no enlaza un FASTA",
                |d| set(d, "/reference/fasta/kind", json!("fastq")),
                "la referencia no enlaza un FASTA",
            ),
            case(
                "ruta de la referencia absoluta",
                |d| set(d, "/reference/fasta/path", json!("/etc/ref.fa")),
                "ruta de artefacto no canónica",
            ),
            case(
                "ruta de la referencia con ..",
                |d| set(d, "/reference/fasta/path", json!("a/../ref.fa")),
                "ruta de artefacto no canónica",
            ),
            case(
                "ruta de una salida absoluta",
                |d| set(d, "/outputs/0/artifact/path", json!("/etc/x.fastq.gz")),
                "ruta de artefacto no canónica",
            ),
            case(
                "tamaño de la referencia que es texto",
                |d| set(d, "/reference/fasta/size_bytes", json!("1")),
                "tamaño de artefacto inválido",
            ),
            case(
                "tamaño de una salida negativo",
                |d| set(d, "/outputs/0/artifact/size_bytes", json!(-1)),
                "tamaño de artefacto inválido",
            ),
            case(
                "algoritmo de la huella de la referencia",
                |d| set(d, "/reference/fasta/checksum/algorithm", json!("md5")),
                "algorithm incompatible",
            ),
            case(
                "huella de la referencia en mayúsculas",
                |d| set(d, "/reference/fasta/checksum/digest", json!("A".repeat(64))),
                "checksum.digest debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "huella de la referencia corta",
                |d| set(d, "/reference/fasta/checksum/digest", json!("a".repeat(63))),
                "checksum.digest debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "huella de la referencia con un carácter que no es hexadecimal",
                |d| set(d, "/reference/fasta/checksum/digest", json!("g".repeat(64))),
                "checksum.digest debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "huella de una salida en mayúsculas",
                |d| {
                    set(
                        d,
                        "/outputs/0/artifact/checksum/digest",
                        json!("A".repeat(64)),
                    )
                },
                "checksum.digest debe ser SHA-256 hexadecimal minúsculo",
            ),
            case(
                "huella de una salida con un carácter que no es hexadecimal",
                |d| {
                    set(
                        d,
                        "/outputs/0/artifact/checksum/digest",
                        json!("g".repeat(64)),
                    )
                },
                "checksum.digest debe ser SHA-256 hexadecimal minúsculo",
            ),
        ],
    );
    ledger.accepted(
        "raíz/plain (lado aceptado)",
        &plain,
        vec![
            accepted("longitud de lectura uno, recontada", |d| {
                reshape_summary(d, 1, &[3, 3]);
            }),
            accepted("longitud de lectura 150, recontada", |d| {
                reshape_summary(d, 150, &[3, 3]);
                // Control del control: el documento recompuesto y el original
                // coinciden salvo en la mezcla de contaminación, que no se
                // juzga aquí.
                set(d, "/final_fastq_summary/1/own_pair_count", json!(3));
            }),
        ],
    );
    // El perfil con amplificación: la terna va junta en el otro sentido.
    let pcr = pcr();
    ledger.rejected(
        "raíz/pcr",
        &pcr,
        vec![
            case(
                "perfil sin amplificación con el recurso de la amplificación",
                |d| set(d, "/profile/id", json!("gps-hiseqx-pcrfree-pe150-v1")),
                "resource_path incompatible",
            ),
            case(
                "recurso del perfil con amplificación cambiado",
                |d| set(d, "/profile/resource_path", json!("resources/f4/otro.json")),
                "resource_path incompatible",
            ),
            case(
                "huella del perfil con amplificación",
                |d| set(d, "/profile/sha256", json!("0".repeat(64))),
                "sha256 incompatible",
            ),
            case(
                "preparación del perfil con amplificación",
                |d| set(d, "/profile/preparation", json!("pcr_free")),
                "preparation incompatible",
            ),
            case(
                "v2 con amplificación",
                |d| set(d, "/schema_version", json!("2")),
                "la versión del documento no corresponde a las etapas que declara",
            ),
            case(
                "catálogo aleatorio del tartamudeo con amplificación",
                |d| set(d, "/rng_catalog", json!("f4-technical-v4")),
                "rng_catalog incompatible",
            ),
            case(
                "catálogo aleatorio del linaje sin linaje",
                |d| set(d, "/rng_catalog", json!("f4-technical-v6")),
                "rng_catalog incompatible",
            ),
        ],
    );
    let stutter = stutter_only();
    ledger.rejected(
        "raíz/stutter",
        &stutter,
        vec![
            case(
                "v2 con tartamudeo",
                |d| set(d, "/schema_version", json!("2")),
                "la versión del documento no corresponde a las etapas que declara",
            ),
            case(
                "catálogo aleatorio sin etapas",
                |d| set(d, "/rng_catalog", json!("f4-technical-v3")),
                "rng_catalog incompatible",
            ),
            case(
                "catálogo aleatorio de la amplificación sin amplificar",
                |d| set(d, "/rng_catalog", json!("f4-technical-v5")),
                "rng_catalog incompatible",
            ),
            case(
                "sin la sección de tartamudeo en v3",
                |d| del(d, "/artifacts/stutter"),
                "la versión del documento no corresponde a las etapas que declara",
            ),
        ],
    );
    let clonal = clonal();
    ledger.rejected(
        "raíz/clonal",
        &clonal,
        vec![
            case(
                "catálogo aleatorio de la amplificación con linaje",
                |d| set(d, "/rng_catalog", json!("f4-technical-v5")),
                "rng_catalog incompatible",
            ),
            case(
                "catálogo aleatorio del tartamudeo con linaje",
                |d| set(d, "/rng_catalog", json!("f4-technical-v4")),
                "rng_catalog incompatible",
            ),
        ],
    );
    // Un JSON que no se lee no pasa por la tabla (no es un documento alterado
    // en un punto): se juzga sobre los bytes, con el texto exacto del error
    // del lector de JSON.
    for (label, bytes) in [
        ("bytes sin JSON", b"{no es json".to_vec()),
        ("bytes vacíos", Vec::new()),
        ("JSON truncado", b"{\"manifest_kind\": ".to_vec()),
    ] {
        let expected = format!(
            "JSON no legible: {}",
            serde_json::from_slice::<Value>(&bytes).expect_err("no es JSON")
        );
        match validate_event_manifest_bytes(Path::new(DOCUMENT_PATH), &bytes) {
            Err(error) => {
                let detail = detail_of(&error);
                if detail != expected {
                    ledger.failures.push(format!(
                        "CASE-FAIL {label}: motivo {detail:?}, esperado {expected:?}"
                    ));
                }
            }
            Ok(()) => ledger.failures.push(format!(
                "CASE-FAIL {label}: se aceptó, y tenía que rechazarse (esperado {expected:?})"
            )),
        }
    }
    ledger.finish();
}

// ---------------------------------------------------------------------------
// Individuos seleccionados, coberturas y resúmenes FASTQ
// ---------------------------------------------------------------------------

const SELECTED_ORDER: &str = "individuos, bibliotecas o grupos no son únicos y ordenados";
const COVERAGE_ORDER: &str =
    "coberturas no son únicas, ordenadas o compatibles con las bibliotecas";
const SUMMARY_ORDER: &str = "resúmenes FASTQ no son únicos, ordenados o compatibles";

/// **Los individuos seleccionados, las coberturas y los resúmenes FASTQ.**
///
/// Las tres listas comparten la forma de su guarda de unicidad y orden
/// (`selected.get(..) != Some(..) || previous >= actual || !insert(..)`), y
/// cada una tiene la suya: un caso por mitad **y por lista**.
#[test]
fn the_selected_individuals_the_coverage_and_the_summaries_are_judged_exactly() {
    let mut ledger = Ledger::default();
    let plain = plain();
    let mut cases = vec![
        case(
            "ningún individuo seleccionado",
            |d| set(d, "/selected_individuals", json!([])),
            "el evento no selecciona individuos",
        ),
        case(
            "identificador de individuo inválido",
            |d| {
                set(
                    d,
                    "/selected_individuals/0/individual_id",
                    json!("no valido!"),
                )
            },
            "ID individual inválido",
        ),
        case(
            "identificador de biblioteca inválido",
            |d| set(d, "/selected_individuals/0/library_id", json!("no valido!")),
            "ID library inválido",
        ),
        case(
            "identificador de grupo de lectura inválido",
            |d| {
                set(
                    d,
                    "/selected_individuals/0/read_group_id",
                    json!("no valido!"),
                )
            },
            "ID read_group inválido",
        ),
        case(
            "individuos en orden inverso",
            |d| swap_items(d, "/selected_individuals", 0, 1),
            SELECTED_ORDER,
        ),
        case(
            "individuo repetido",
            |d| {
                let first = text(d, "/selected_individuals/0/individual_id");
                set(d, "/selected_individuals/1/individual_id", json!(first));
            },
            SELECTED_ORDER,
        ),
        case(
            "grupo de lectura repetido",
            |d| {
                let first = text(d, "/selected_individuals/0/read_group_id");
                set(d, "/selected_individuals/1/read_group_id", json!(first));
            },
            SELECTED_ORDER,
        ),
    ];
    // Los dos UDI, cada uno con cada mitad de su guarda (longitud y alfabeto).
    for field in ["i7", "i5"] {
        for (label, value) in [
            ("de nueve bases", "ACGTACGTA"),
            ("de siete bases", "ACGTACG"),
            ("con una N", "ACGTACGN"),
            ("en minúsculas", "acgtacgt"),
        ] {
            cases.push(case(
                format!("UDI {field} {label}"),
                move |d| set(d, &format!("/selected_individuals/0/{field}"), json!(value)),
                "UDI de biblioteca inválido",
            ));
        }
    }
    cases.extend([
        case(
            "cobertura de un individuo desconocido",
            |d| set(d, "/coverage/0/individual_id", json!("sample-Z")),
            COVERAGE_ORDER,
        ),
        case(
            "cobertura con la biblioteca de otro",
            |d| {
                let other = text(d, "/coverage/1/library_id");
                set(d, "/coverage/0/library_id", json!(other));
            },
            COVERAGE_ORDER,
        ),
        case(
            "coberturas en orden inverso",
            |d| swap_items(d, "/coverage", 0, 1),
            COVERAGE_ORDER,
        ),
        case(
            "cobertura repetida",
            |d| {
                let first = d["coverage"][0].clone();
                set(d, "/coverage/1", first);
            },
            COVERAGE_ORDER,
        ),
        case(
            "objetivo de cobertura negativo",
            |d| set(d, "/coverage/0/target_x", json!(-1.0)),
            "target_x de cobertura inválida",
        ),
        case(
            "objetivo de cobertura que es texto",
            |d| set(d, "/coverage/0/target_x", json!("1")),
            "target_x de cobertura inválida",
        ),
        case(
            "cobertura básica negativa",
            |d| set(d, "/coverage/0/basic_x", json!(-0.5)),
            "basic_x de cobertura inválida",
        ),
        case(
            "cobertura propia negativa",
            |d| set(d, "/coverage/0/final_x/own_x", json!(-0.5)),
            "own_x de cobertura inválida",
        ),
        case(
            "cobertura contaminante negativa",
            |d| set(d, "/coverage/0/final_x/contaminant_x", json!(-0.5)),
            "contaminant_x de cobertura inválida",
        ),
        case(
            "denominador de cobertura cero",
            |d| set(d, "/coverage/0/denominator_bases", json!(0)),
            "denominador de cobertura inválido",
        ),
        case(
            "denominador de cobertura negativo",
            |d| set(d, "/coverage/0/denominator_bases", json!(-1)),
            "denominador de cobertura inválido",
        ),
        case(
            "denominador de cobertura que es texto",
            |d| set(d, "/coverage/0/denominator_bases", json!("1000")),
            "denominador de cobertura inválido",
        ),
        case(
            "recuento de N que es texto",
            |d| set(d, "/coverage/0/n_bases", json!("0")),
            "recuento N de cobertura inválido",
        ),
        case(
            "recuento de N negativo",
            |d| set(d, "/coverage/0/n_bases", json!(-1)),
            "recuento N de cobertura inválido",
        ),
        case(
            "resumen de un individuo desconocido",
            |d| set(d, "/final_fastq_summary/0/individual_id", json!("sample-Z")),
            SUMMARY_ORDER,
        ),
        case(
            "resumen con la biblioteca de otro",
            |d| {
                let other = text(d, "/final_fastq_summary/1/library_id");
                set(d, "/final_fastq_summary/0/library_id", json!(other));
            },
            SUMMARY_ORDER,
        ),
        case(
            "resúmenes en orden inverso",
            |d| swap_items(d, "/final_fastq_summary", 0, 1),
            SUMMARY_ORDER,
        ),
        case(
            "resumen repetido",
            |d| {
                let first = d["final_fastq_summary"][0].clone();
                set(d, "/final_fastq_summary/1", first);
            },
            SUMMARY_ORDER,
        ),
    ]);
    for field in [
        "pair_count",
        "read_count",
        "base_count",
        "own_pair_count",
        "contaminant_pair_count",
        "own_base_count",
        "contaminant_base_count",
        "n_base_count",
    ] {
        cases.push(case(
            format!("{field} que es texto"),
            move |d| set(d, &format!("/final_fastq_summary/0/{field}"), json!("1")),
            format!("{field} no es entero sin signo"),
        ));
    }
    cases.extend([
        case(
            "bases propias y contaminantes que desbordan",
            |d| {
                set(d, "/final_fastq_summary/0/own_base_count", json!(u64::MAX));
                set(d, "/final_fastq_summary/0/contaminant_base_count", json!(1));
            },
            "contadores agregados FASTQ desbordan",
        ),
        case(
            "parejas cuyo doble desborda",
            |d| set(d, "/final_fastq_summary/0/pair_count", json!(1_u64 << 63)),
            "contador de lecturas FASTQ desborda",
        ),
        case(
            "lecturas cuyo producto por la longitud desborda",
            |d| set(d, "/final_fastq_summary/0/pair_count", json!(1_u64 << 62)),
            "contador de bases FASTQ desborda",
        ),
        case(
            "propias más contaminantes no suman las parejas",
            |d| set(d, "/final_fastq_summary/0/own_pair_count", json!(4)),
            "contadores agregados FASTQ incoherentes",
        ),
        case(
            "parejas propias que desbordan al sumar",
            |d| {
                set(d, "/final_fastq_summary/0/own_pair_count", json!(u64::MAX));
                set(d, "/final_fastq_summary/0/contaminant_pair_count", json!(4));
            },
            "contadores agregados FASTQ incoherentes",
        ),
        case(
            "lecturas distintas del doble de las parejas",
            |d| set(d, "/final_fastq_summary/0/read_count", json!(7)),
            "contadores agregados FASTQ incoherentes",
        ),
        case(
            "bases distintas de lecturas por longitud",
            |d| set(d, "/final_fastq_summary/0/base_count", json!(901)),
            "contadores agregados FASTQ incoherentes",
        ),
        case(
            "bases propias y contaminantes distintas del total",
            |d| set(d, "/final_fastq_summary/0/own_base_count", json!(901)),
            "contadores agregados FASTQ incoherentes",
        ),
        case(
            "más bases N que bases",
            |d| set(d, "/final_fastq_summary/0/n_base_count", json!(901)),
            "contadores agregados FASTQ incoherentes",
        ),
        case(
            "resumen sin cobertura asociada",
            |d| del(d, "/coverage/1"),
            "resumen FASTQ no tiene cobertura asociada",
        ),
        case(
            "N del resumen distintas de las de la cobertura",
            |d| set(d, "/final_fastq_summary/0/n_base_count", json!(1)),
            "resumen FASTQ no coincide con su cobertura final",
        ),
        case(
            "cobertura propia que no sale de las bases propias",
            |d| set(d, "/coverage/0/final_x/own_x", json!(0.95)),
            "resumen FASTQ no coincide con su cobertura final",
        ),
        case(
            "cobertura contaminante que no sale de las bases contaminantes",
            |d| set(d, "/coverage/1/final_x/contaminant_x", json!(0.5)),
            "resumen FASTQ no coincide con su cobertura final",
        ),
        case(
            "cobertura propia con cero negativo",
            |d| {
                // B sin parejas propias: su cobertura propia exacta es +0.0.
                for (field, value) in [
                    ("own_pair_count", 0),
                    ("contaminant_pair_count", 3),
                    ("own_base_count", 0),
                    ("contaminant_base_count", 900),
                ] {
                    set(d, &format!("/final_fastq_summary/1/{field}"), json!(value));
                }
                set(d, "/coverage/1/final_x/own_x", json!(-0.0));
                set(d, "/coverage/1/final_x/contaminant_x", json!(0.9));
            },
            "resumen FASTQ no coincide con su cobertura final",
        ),
        case(
            "cobertura contaminante con cero negativo",
            |d| set(d, "/coverage/0/final_x/contaminant_x", json!(-0.0)),
            "resumen FASTQ no coincide con su cobertura final",
        ),
        case(
            "pares finales que desbordan al sumar",
            |d| {
                with_third_library(d);
                reshape_summary(d, 1, &[(1 << 63) - 1, (1 << 63) - 1, 2]);
            },
            "pares FASTQ finales desbordan",
        ),
        case(
            "resúmenes que no cubren las mismas bibliotecas que las coberturas",
            |d| {
                neutral_demultiplexing(d);
                del(d, "/final_fastq_summary/1");
            },
            "resúmenes FASTQ y coberturas no cubren las mismas bibliotecas",
        ),
        case(
            "falta la cobertura de una biblioteca completa",
            |d| {
                neutral_demultiplexing(d);
                del(d, "/final_fastq_summary/1");
                del(d, "/coverage/1");
            },
            "falta cobertura de una biblioteca completa",
        ),
    ]);
    ledger.rejected("individuos/plain", &plain, cases);
    ledger.accepted(
        "individuos/plain (lado aceptado)",
        &plain,
        vec![
            accepted("objetivo y cobertura básica en cero", |d| {
                set(d, "/coverage/0/target_x", json!(0.0));
                set(d, "/coverage/0/basic_x", json!(0.0));
            }),
            accepted("cobertura propia con cero positivo exacto", |d| {
                for (field, value) in [
                    ("own_pair_count", 0),
                    ("contaminant_pair_count", 3),
                    ("own_base_count", 0),
                    ("contaminant_base_count", 900),
                ] {
                    set(d, &format!("/final_fastq_summary/1/{field}"), json!(value));
                }
                set(d, "/coverage/1/final_x/own_x", json!(0.0));
                set(d, "/coverage/1/final_x/contaminant_x", json!(0.9));
            }),
            accepted("tantas N como bases", |d| {
                set(d, "/final_fastq_summary/0/n_base_count", json!(900));
                set(d, "/coverage/0/n_bases", json!(900));
            }),
            accepted("tres bibliotecas con parejas cuya suma no desborda", |d| {
                with_third_library(d);
                reshape_summary(d, 1, &[(1 << 63) - 1, (1 << 63) - 1, 0]);
                neutral_demultiplexing(d);
            }),
        ],
    );
    ledger.finish();
}

// ---------------------------------------------------------------------------
// Demultiplexado y calidad de error
// ---------------------------------------------------------------------------

const DEMULTIPLEXING_SUMMARY: &str = "resumen de demultiplexado inválido o publica procedencia";
const AGGREGATES_ORDER: &str = "agregados de demultiplexado desordenados";
const AGGREGATES_TOTAL: &str = "agregados de demultiplexado no suman el total";
const FIXED_ADJUSTMENTS: &str =
    "F4 no admite ajustes de error/calidad distintos de quality_shift=0 y error_scale=1";

/// **El resumen de demultiplexado y sus tres listas de agregados.**
///
/// `validate_demultiplexing_aggregates` es **un** ayudante con tres
/// llamantes (reasignadas, descartadas y sin asignar) y dos formas de fila
/// (con y sin destino): un caso por llamante y por forma. El escenario sin
/// amplificar trae reasignaciones reales por salto de índice; los descartes
/// por UDI no reconocido no los produce el programa hoy pero el contrato los
/// admite, así que se parte de un documento con ellos que el juez acepta.
#[test]
fn the_demultiplexing_summary_and_its_aggregates_are_judged_exactly() {
    let mut ledger = Ledger::default();
    let plain = plain();
    let mut cases = vec![
        case(
            "tasa de salto por encima de uno",
            |d| set(d, "/demultiplexing/index_hopping_rate", json!(1.5)),
            DEMULTIPLEXING_SUMMARY,
        ),
        case(
            "tasa de salto negativa",
            |d| set(d, "/demultiplexing/index_hopping_rate", json!(-0.5)),
            DEMULTIPLEXING_SUMMARY,
        ),
        case(
            "tasa de salto que es texto",
            |d| set(d, "/demultiplexing/index_hopping_rate", json!("1")),
            DEMULTIPLEXING_SUMMARY,
        ),
        case(
            "reasignadas que son texto",
            |d| set(d, "/demultiplexing/reassigned_pairs", json!("6")),
            DEMULTIPLEXING_SUMMARY,
        ),
        case(
            "descartadas que son texto",
            |d| set(d, "/demultiplexing/discarded_pairs", json!("0")),
            DEMULTIPLEXING_SUMMARY,
        ),
        case(
            "sin asignar que son texto",
            |d| set(d, "/demultiplexing/unassigned_pairs", json!("0")),
            DEMULTIPLEXING_SUMMARY,
        ),
        case(
            "publica el FASTQ limpio",
            |d| set(d, "/demultiplexing/clean_fastq_published", json!(true)),
            DEMULTIPLEXING_SUMMARY,
        ),
        case(
            "publica la procedencia individual",
            |d| {
                set(
                    d,
                    "/demultiplexing/individual_provenance_published",
                    json!(true),
                );
            },
            DEMULTIPLEXING_SUMMARY,
        ),
        case(
            "más sin asignar que descartadas",
            |d| set(d, "/demultiplexing/unassigned_pairs", json!(1)),
            "parejas sin asignar exceden los descartes",
        ),
        case(
            "más reasignaciones que parejas finales",
            |d| set(d, "/demultiplexing/reassigned_pairs", json!(7)),
            "reasignaciones exceden las parejas FASTQ finales",
        ),
        case(
            "tasa de salto nula con reasignaciones",
            |d| set(d, "/demultiplexing/index_hopping_rate", json!(0.0)),
            "un índice de salto nulo no puede reasignar parejas",
        ),
        case(
            "reasignadas que no son lista",
            |d| set(d, "/demultiplexing/reassigned_by_sample", json!("x")),
            "agregados de demultiplexado inválidos",
        ),
        case(
            "descartadas que no son lista",
            |d| set(d, "/demultiplexing/discarded_by_sample", json!("x")),
            "agregados de demultiplexado inválidos",
        ),
        case(
            "sin asignar que no son lista",
            |d| set(d, "/demultiplexing/unassigned_by_sample", json!("x")),
            "agregados de demultiplexado inválidos",
        ),
        case(
            "fila de reasignadas que no es objeto",
            |d| set(d, "/demultiplexing/reassigned_by_sample/0", json!("x")),
            "agregado de demultiplexado inválido",
        ),
        case(
            "fila de descartadas que no es objeto",
            |d| push(d, "/demultiplexing/discarded_by_sample", json!(1)),
            "agregado de demultiplexado inválido",
        ),
        case(
            "fila de sin asignar que no es objeto",
            |d| push(d, "/demultiplexing/unassigned_by_sample", json!(1)),
            "agregado de demultiplexado inválido",
        ),
        case(
            "reasignada sin destino",
            |d| {
                del(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/destination_library_id",
                )
            },
            "campos inválidos en reassigned_by_sample",
        ),
        case(
            "origen de una reasignación desconocido",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/source_library_id",
                    json!("LIB_otra"),
                );
            },
            "origen de demultiplexado inválido",
        ),
        case(
            "origen de una reasignación que no es texto",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/source_library_id",
                    json!(7),
                )
            },
            "origen de demultiplexado inválido",
        ),
        case(
            "origen de una reasignación con sintaxis inválida",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/source_library_id",
                    json!("no valido!"),
                );
            },
            "origen de demultiplexado inválido",
        ),
        case(
            "destino de una reasignación desconocido",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/destination_library_id",
                    json!("LIB_otra"),
                );
            },
            "destino de demultiplexado inválido",
        ),
        case(
            "destino de una reasignación que no es texto",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/destination_library_id",
                    json!(7),
                );
            },
            "destino de demultiplexado inválido",
        ),
        case(
            "reasignación circular",
            |d| {
                let source = text(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/source_library_id",
                );
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/destination_library_id",
                    json!(source),
                );
            },
            "reasignación de demultiplexado circular",
        ),
        case(
            "causa de una reasignación que no es el salto de índice",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/cause",
                    json!("unrecognized_udi"),
                );
            },
            "causa de demultiplexado inválida",
        ),
        case(
            "fila de reasignadas vacía",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/pair_count",
                    json!(0),
                )
            },
            "agregado de demultiplexado vacío",
        ),
        case(
            "recuento de una reasignación que es texto",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/pair_count",
                    json!("3"),
                )
            },
            "pair_count no es entero sin signo",
        ),
        case(
            "recuentos de reasignadas que desbordan",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/pair_count",
                    json!(u64::MAX),
                );
            },
            "recuento de demultiplexado desbordado",
        ),
        case(
            "reasignadas en orden inverso",
            |d| swap_items(d, "/demultiplexing/reassigned_by_sample", 0, 1),
            AGGREGATES_ORDER,
        ),
        case(
            "reasignada repetida",
            |d| {
                let first = d["demultiplexing"]["reassigned_by_sample"][0].clone();
                set(d, "/demultiplexing/reassigned_by_sample/1", first);
            },
            AGGREGATES_ORDER,
        ),
        case(
            "reasignadas que no suman el total",
            |d| {
                set(
                    d,
                    "/demultiplexing/reassigned_by_sample/0/pair_count",
                    json!(4),
                )
            },
            AGGREGATES_TOTAL,
        ),
    ];
    // Los mismos ejes sobre las dos listas sin destino.
    let discards = with_discards(&plain);
    let mut discard_cases = Vec::new();
    for (field, pairs_field) in [
        ("discarded_by_sample", "discarded_pairs"),
        ("unassigned_by_sample", "unassigned_pairs"),
    ] {
        let list = format!("/demultiplexing/{field}");
        let pairs = format!("/demultiplexing/{pairs_field}");
        let (a, b, c, e, f, g, h) = (
            list.clone(),
            list.clone(),
            list.clone(),
            list.clone(),
            list.clone(),
            list.clone(),
            list.clone(),
        );
        discard_cases.push(case(
            format!("{field}: fila con destino"),
            move |d| {
                let library = text(d, &format!("{a}/0/source_library_id"));
                set(d, &format!("{a}/0/destination_library_id"), json!(library));
            },
            format!("campos inválidos en {field}"),
        ));
        discard_cases.push(case(
            format!("{field}: origen desconocido"),
            move |d| set(d, &format!("{b}/0/source_library_id"), json!("LIB_otra")),
            "origen de demultiplexado inválido",
        ));
        discard_cases.push(case(
            format!("{field}: causa que no es UDI no reconocido"),
            move |d| set(d, &format!("{c}/0/cause"), json!("index_hopping")),
            "causa de demultiplexado inválida",
        ));
        discard_cases.push(case(
            format!("{field}: fila vacía"),
            move |d| set(d, &format!("{e}/0/pair_count"), json!(0)),
            "agregado de demultiplexado vacío",
        ));
        discard_cases.push(case(
            format!("{field}: filas en orden inverso"),
            move |d| swap_items(d, &f, 0, 1),
            AGGREGATES_ORDER,
        ));
        discard_cases.push(case(
            format!("{field}: fila repetida"),
            move |d| {
                let first = d.pointer(&format!("{g}/0")).expect("fila").clone();
                set(d, &format!("{g}/1"), first);
            },
            AGGREGATES_ORDER,
        ));
        discard_cases.push(case(
            format!("{field}: recuentos que desbordan"),
            move |d| set(d, &format!("{h}/0/pair_count"), json!(u64::MAX)),
            "recuento de demultiplexado desbordado",
        ));
        discard_cases.push(case(
            format!("{field}: total distinto del declarado"),
            move |d| set(d, &pairs, json!(2)),
            if field == "discarded_by_sample" {
                // Descartadas = 2 y sin asignar = 3 exceden antes que la suma.
                "parejas sin asignar exceden los descartes"
            } else {
                AGGREGATES_TOTAL
            },
        ));
    }
    discard_cases.push(case(
        "descartadas con un total declarado mayor que la suma",
        |d| set(d, "/demultiplexing/discarded_pairs", json!(4)),
        AGGREGATES_TOTAL,
    ));
    discard_cases.push(case(
        "descartes y sin asignar repartidos de otra forma",
        |d| {
            set(
                d,
                "/demultiplexing/unassigned_by_sample/0/pair_count",
                json!(1),
            );
            set(
                d,
                "/demultiplexing/unassigned_by_sample/1/pair_count",
                json!(2),
            );
        },
        "descartes y parejas sin asignar no coinciden por biblioteca",
    ));
    cases.shrink_to_fit();
    ledger.rejected("demultiplexado/plain", &plain, cases);
    ledger.rejected("demultiplexado/descartes", &discards, discard_cases);
    // El lado aceptado de los límites ya lo fijan los controles de las dos
    // tablas: el escenario real tiene la tasa de salto en 1,0 (extremo
    // inclusivo del intervalo) y tantas reasignaciones (6) como parejas
    // finales, y el de los descartes acepta `sin asignar == descartadas`.
    ledger.finish();
}

/// **El modo de errores, la marca de calidad y los recuentos por mate.**
#[test]
fn the_error_quality_section_is_judged_exactly() {
    let mut ledger = Ledger::default();
    let plain = plain();
    let mut cases = vec![
        case(
            "modo pedido desconocido",
            |d| set(d, "/error_quality/requested_mode", json!("quizá")),
            "modo de errores de secuenciación inválido",
        ),
        case(
            "modo efectivo desconocido",
            |d| set(d, "/error_quality/effective_mode", json!("quizá")),
            "modo de errores de secuenciación inválido",
        ),
        case(
            "marca de calidad que no es booleana",
            |d| {
                set(
                    d,
                    "/error_quality/qualities_represent_realized_error_probability",
                    json!("sí"),
                );
            },
            "marca de calidad de error inválida",
        ),
        case(
            "modo pedido distinto del efectivo",
            |d| set(d, "/error_quality/requested_mode", json!("disabled")),
            "modo de error y marca de calidad incoherentes",
        ),
        case(
            "marca de calidad que no corresponde al modo",
            |d| {
                set(
                    d,
                    "/error_quality/qualities_represent_realized_error_probability",
                    json!(false),
                );
            },
            "modo de error y marca de calidad incoherentes",
        ),
        case(
            "oportunidades que son texto",
            |d| set(d, "/error_quality/r1/opportunities", json!("900")),
            "opportunities no es entero sin signo",
        ),
        case(
            "errores que son texto",
            |d| set(d, "/error_quality/r2/errors", json!("0")),
            "errors no es entero sin signo",
        ),
        case(
            "errores de R1 por encima de las oportunidades",
            |d| set(d, "/error_quality/r1/errors", json!(901)),
            "errores realizados exceden oportunidades",
        ),
        case(
            "errores de R2 por encima de las oportunidades",
            |d| set(d, "/error_quality/r2/errors", json!(901)),
            "errores realizados exceden oportunidades",
        ),
    ];
    for (label, mate, field, value) in [
        (
            "desplazamiento de calidad de R1",
            "r1",
            "quality_shift",
            json!(1),
        ),
        ("escala de error de R1", "r1", "error_scale", json!(2.0)),
        (
            "desplazamiento de calidad de R2",
            "r2",
            "quality_shift",
            json!(-1),
        ),
        ("escala de error de R2", "r2", "error_scale", json!(0.5)),
    ] {
        cases.push(case(
            label,
            move |d| set(d, &format!("/error_quality/{mate}/{field}"), value.clone()),
            FIXED_ADJUSTMENTS,
        ));
    }
    ledger.rejected("calidad/plain", &plain, cases);
    ledger.accepted(
        "calidad/plain (lado aceptado)",
        &plain,
        vec![
            accepted("tantos errores como oportunidades en R1", |d| {
                set(d, "/error_quality/r1/errors", json!(900));
            }),
            accepted("tantos errores como oportunidades en R2", |d| {
                set(d, "/error_quality/r2/errors", json!(900));
            }),
        ],
    );
    // Con los errores desactivados, ningún mate puede traer errores.
    let mut disabled = plain.clone();
    set(
        &mut disabled,
        "/error_quality/requested_mode",
        json!("disabled"),
    );
    set(
        &mut disabled,
        "/error_quality/effective_mode",
        json!("disabled"),
    );
    set(
        &mut disabled,
        "/error_quality/qualities_represent_realized_error_probability",
        json!(false),
    );
    ledger.rejected(
        "calidad/desactivada",
        &disabled,
        vec![
            case(
                "errores en R1 con los errores desactivados",
                |d| set(d, "/error_quality/r1/errors", json!(1)),
                "sequencing_errors disabled exige cero errores por mate",
            ),
            case(
                "errores en R2 con los errores desactivados",
                |d| set(d, "/error_quality/r2/errors", json!(1)),
                "sequencing_errors disabled exige cero errores por mate",
            ),
            case(
                "errores en los dos mates con los errores desactivados",
                |d| {
                    set(d, "/error_quality/r1/errors", json!(1));
                    set(d, "/error_quality/r2/errors", json!(2));
                },
                "sequencing_errors disabled exige cero errores por mate",
            ),
        ],
    );
    ledger.finish();
}

// ---------------------------------------------------------------------------
// Artefactos técnicos fijos, amplificación y linaje
// ---------------------------------------------------------------------------

const FIXED_ARTIFACTS: &str = "artefactos técnicos fijos inválidos";
const CHIMERA: &str = "artefacto de quimera inválido";
const FAMILIES: &str = "los recuentos de familias y duplicados no cuadran";
const HISTOGRAM_ORDER: &str = "el histograma de tamaños de familia no está en orden estricto";
const HISTOGRAM_OVERFLOW: &str = "el histograma de familias desborda";
const HISTOGRAM_SUMS: &str =
    "el histograma de tamaños de familia no suma las familias y parejas publicadas";
const PCR_COUNTS: &str =
    "los recuentos de intentos de composición o de familias truncadas no cuadran";
const LINEAGE_TOGETHER: &str =
    "la sección de tartamudeo clonal y sus tres recuentos de linaje van juntos o no van";

/// **Los artefactos fijos (quimera, PCR ausente, ópticos, GC y adaptadores).**
#[test]
fn the_fixed_technical_artifacts_are_judged_exactly() {
    let mut ledger = Ledger::default();
    let plain = plain();
    ledger.rejected(
        "fijos/plain",
        &plain,
        vec![
            case(
                "bloque de tartamudeo clonal en la sección de una biblioteca sin amplificar",
                |d| set(d, "/artifacts/pcr/stutter", json!({})),
                "campos inválidos en artifacts.pcr",
            ),
            case(
                "tasa de quimera por encima de uno",
                |d| set(d, "/artifacts/chimera/rate", json!(1.1)),
                CHIMERA,
            ),
            case(
                "tasa de quimera negativa",
                |d| set(d, "/artifacts/chimera/rate", json!(-0.1)),
                CHIMERA,
            ),
            case(
                "tasa de quimera que es texto",
                |d| set(d, "/artifacts/chimera/rate", json!("0")),
                CHIMERA,
            ),
            case(
                "parejas quiméricas que son texto",
                |d| set(d, "/artifacts/chimera/realized_pairs", json!("1")),
                CHIMERA,
            ),
            case(
                "parejas quiméricas que son texto con tasa positiva",
                |d| {
                    set(d, "/artifacts/chimera/rate", json!(0.5));
                    set(d, "/artifacts/chimera/realized_pairs", json!("1"));
                },
                CHIMERA,
            ),
            case(
                "parejas quiméricas con tasa cero",
                |d| set(d, "/artifacts/chimera/realized_pairs", json!(1)),
                CHIMERA,
            ),
            case(
                "amplificación declarada con un valor que no es booleano",
                |d| set(d, "/artifacts/pcr/applicable", json!(0)),
                FIXED_ARTIFACTS,
            ),
            case(
                "amplificación declarada como texto",
                |d| set(d, "/artifacts/pcr/applicable", json!("false")),
                FIXED_ARTIFACTS,
            ),
            case(
                "duplicados ópticos realizados",
                |d| set(d, "/artifacts/optical_duplicates/realized_pairs", json!(1)),
                FIXED_ARTIFACTS,
            ),
            case(
                "modo de GC con curva sin amplificación",
                |d| set(d, "/artifacts/gc/mode", json!("curve")),
                FIXED_ARTIFACTS,
            ),
            case(
                "modo de GC desconocido",
                |d| set(d, "/artifacts/gc/mode", json!("otro")),
                FIXED_ARTIFACTS,
            ),
            case(
                "modo de GC que no es texto",
                |d| set(d, "/artifacts/gc/mode", json!(1)),
                FIXED_ARTIFACTS,
            ),
            case(
                "adaptadores activados",
                |d| set(d, "/artifacts/adapters/enabled", json!(true)),
                FIXED_ARTIFACTS,
            ),
        ],
    );
    ledger.accepted(
        "fijos/plain (lado aceptado)",
        &plain,
        vec![
            accepted("tasa de quimera uno sin parejas quiméricas", |d| {
                set(d, "/artifacts/chimera/rate", json!(1.0));
            }),
            accepted("tasa de quimera positiva con parejas quiméricas", |d| {
                set(d, "/artifacts/chimera/rate", json!(0.5));
                set(d, "/artifacts/chimera/realized_pairs", json!(7));
            }),
        ],
    );
    let pcr = pcr();
    ledger.accepted(
        "fijos/pcr (lado aceptado)",
        &pcr,
        vec![accepted("modo de GC con curva y amplificación", |d| {
            set(d, "/artifacts/gc/mode", json!("curve"));
        })],
    );
    ledger.finish();
}

/// **La sección de amplificación.** Los recuentos de familias, el histograma
/// ponderado, la fracción realizada y los contadores de intentos son las
/// relaciones que un esquema JSON no puede expresar, y cada una tiene sus
/// casos, uno por mitad de cada disyunción y con su límite aceptado.
#[test]
fn the_pcr_section_is_judged_exactly() {
    let mut ledger = Ledger::default();
    let pcr = pcr();
    let mut cases = vec![
        case(
            "algoritmo de la amplificación",
            |d| set(d, "/artifacts/pcr/algorithm", json!("otro")),
            "algorithm incompatible",
        ),
        case(
            "catálogo aleatorio del linaje en la amplificación sin linaje",
            |d| set(d, "/artifacts/pcr/rng_catalog", json!("f4-technical-v6")),
            "rng_catalog incompatible",
        ),
        case(
            "un recuento de linaje suelto",
            |d| set(d, "/artifacts/pcr/lineage_slips_by_locus", json!([])),
            LINEAGE_TOGETHER,
        ),
        case(
            "solo el bloque de tartamudeo clonal",
            |d| {
                set(d, "/rng_catalog", json!("f4-technical-v6"));
                set(d, "/artifacts/pcr/stutter", json!({}));
            },
            LINEAGE_TOGETHER,
        ),
        case(
            "método de inversión de la fracción de duplicados",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/duplicate_fraction_inversion/method",
                    json!("otro"),
                )
            },
            "method incompatible",
        ),
        case(
            "cota de recursos publicada como medida",
            |d| set(d, "/artifacts/pcr/resource_bound/measured", json!(true)),
            "la cota de recursos se publica calculada, no medida",
        ),
        case(
            "cota de recursos con la marca que no es booleana",
            |d| set(d, "/artifacts/pcr/resource_bound/measured", json!("no")),
            "la cota de recursos se publica calculada, no medida",
        ),
        case(
            "más familias que parejas",
            |d| set(d, "/artifacts/pcr/families", json!(14)),
            FAMILIES,
        ),
        case(
            "parejas duplicadas que no son parejas menos familias",
            |d| set(d, "/artifacts/pcr/duplicate_pairs", json!(7)),
            FAMILIES,
        ),
        case(
            "familias que son texto",
            |d| set(d, "/artifacts/pcr/families", json!("7")),
            "families no es entero sin signo",
        ),
        case(
            "parejas que son texto",
            |d| set(d, "/artifacts/pcr/pairs", json!("13")),
            "pairs no es entero sin signo",
        ),
        case(
            "parejas duplicadas que son texto",
            |d| set(d, "/artifacts/pcr/duplicate_pairs", json!("6")),
            "duplicate_pairs no es entero sin signo",
        ),
        case(
            "tamaño de familia cero",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/family_size_histogram/0/family_size",
                    json!(0),
                )
            },
            HISTOGRAM_ORDER,
        ),
        case(
            "celda sin familias",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/family_size_histogram/0/families",
                    json!(0),
                )
            },
            HISTOGRAM_ORDER,
        ),
        case(
            "histograma en orden inverso",
            |d| swap_items(d, "/artifacts/pcr/family_size_histogram", 0, 1),
            HISTOGRAM_ORDER,
        ),
        case(
            "tamaño de familia repetido",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/family_size_histogram/1/family_size",
                    json!(1),
                )
            },
            HISTOGRAM_ORDER,
        ),
        case(
            "tamaño de familia que es texto",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/family_size_histogram/0/family_size",
                    json!("1"),
                )
            },
            "family_size no es entero sin signo",
        ),
        case(
            "familias de una celda que son texto",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/family_size_histogram/0/families",
                    json!("2"),
                )
            },
            "families no es entero sin signo",
        ),
        case(
            "familias del histograma que desbordan al sumar",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/family_size_histogram",
                    json!([
                        {"family_size": 1, "families": u64::MAX},
                        {"family_size": 2, "families": 1},
                    ]),
                );
            },
            HISTOGRAM_OVERFLOW,
        ),
        case(
            "producto tamaño por familias que desborda",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/family_size_histogram",
                    json!([{"family_size": 2, "families": 1_u64 << 63}]),
                );
            },
            HISTOGRAM_OVERFLOW,
        ),
        case(
            "parejas del histograma que desbordan al sumar",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/family_size_histogram",
                    json!([
                        {"family_size": 2, "families": 1_u64 << 62},
                        {"family_size": 4, "families": 1_u64 << 61},
                    ]),
                );
            },
            HISTOGRAM_OVERFLOW,
        ),
        case(
            "familias del histograma distintas de las publicadas",
            |d| {
                set(d, "/artifacts/pcr/families", json!(8));
                set(d, "/artifacts/pcr/duplicate_pairs", json!(5));
            },
            HISTOGRAM_SUMS,
        ),
        case(
            "parejas del histograma distintas de las publicadas",
            |d| {
                set(d, "/artifacts/pcr/pairs", json!(14));
                set(d, "/artifacts/pcr/duplicate_pairs", json!(7));
            },
            HISTOGRAM_SUMS,
        ),
        case(
            "fracción realizada ausente",
            |d| del(d, "/artifacts/pcr/realized_duplicate_fraction"),
            "campos inválidos en artifacts.pcr",
        ),
        case(
            "fracción realizada que es texto",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/realized_duplicate_fraction",
                    json!("0.46"),
                )
            },
            "falta la fracción realizada de duplicados",
        ),
        case(
            "fracción realizada con cero negativo",
            |d| {
                for (field, value) in [
                    ("pairs", json!(0)),
                    ("families", json!(0)),
                    ("duplicate_pairs", json!(0)),
                    ("family_size_histogram", json!([])),
                    ("realized_duplicate_fraction", json!(-0.0)),
                    ("gc_attempts_exhausted_families", json!(0)),
                    ("gc_attempts_total", json!(0)),
                    ("families_truncated_at_end", json!(0)),
                ] {
                    set(d, &format!("/artifacts/pcr/{field}"), value);
                }
            },
            "la fracción realizada de duplicados no sale de los recuentos publicados",
        ),
        case(
            "fracción realizada que no sale de los recuentos",
            |d| set(d, "/artifacts/pcr/realized_duplicate_fraction", json!(0.5)),
            "la fracción realizada de duplicados no sale de los recuentos publicados",
        ),
        case(
            "familias agotadas por encima de las familias",
            |d| set(d, "/artifacts/pcr/gc_attempts_exhausted_families", json!(8)),
            PCR_COUNTS,
        ),
        case(
            "menos intentos de composición que familias",
            |d| set(d, "/artifacts/pcr/gc_attempts_total", json!(6)),
            PCR_COUNTS,
        ),
        case(
            "familias truncadas por encima de las familias",
            |d| set(d, "/artifacts/pcr/families_truncated_at_end", json!(8)),
            PCR_COUNTS,
        ),
        case(
            "familias agotadas que son texto",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/gc_attempts_exhausted_families",
                    json!("0"),
                )
            },
            "gc_attempts_exhausted_families no es entero sin signo",
        ),
        case(
            "intentos de composición que son texto",
            |d| set(d, "/artifacts/pcr/gc_attempts_total", json!("7")),
            "gc_attempts_total no es entero sin signo",
        ),
        case(
            "familias truncadas que son texto",
            |d| set(d, "/artifacts/pcr/families_truncated_at_end", json!("1")),
            "families_truncated_at_end no es entero sin signo",
        ),
    ];
    cases.shrink_to_fit();
    ledger.rejected("pcr/pcr", &pcr, cases);
    ledger.accepted(
        "pcr/pcr (lado aceptado)",
        &pcr,
        vec![
            accepted("los tres contadores en su límite inclusivo", |d| {
                set(d, "/artifacts/pcr/gc_attempts_exhausted_families", json!(7));
                set(d, "/artifacts/pcr/gc_attempts_total", json!(7));
                set(d, "/artifacts/pcr/families_truncated_at_end", json!(7));
            }),
            accepted("una amplificación sin parejas", |d| {
                set(d, "/artifacts/pcr/pairs", json!(0));
                set(d, "/artifacts/pcr/families", json!(0));
                set(d, "/artifacts/pcr/duplicate_pairs", json!(0));
                set(d, "/artifacts/pcr/family_size_histogram", json!([]));
                set(d, "/artifacts/pcr/realized_duplicate_fraction", json!(0.0));
                set(d, "/artifacts/pcr/gc_attempts_exhausted_families", json!(0));
                set(d, "/artifacts/pcr/gc_attempts_total", json!(0));
                set(d, "/artifacts/pcr/families_truncated_at_end", json!(0));
            }),
        ],
    );
    ledger.finish();
}

/// **Los tres recuentos de linaje.** Un solo ayudante (`locus_slip_counts`)
/// sirve a las tres listas, con dos nombres de contador distintos: un caso por
/// lista y por nombre.
#[test]
fn the_lineage_counts_are_judged_exactly() {
    let mut ledger = Ledger::default();
    let clonal = clonal();
    const SAME_LOCI: &str = "los tres recuentos de linaje no nombran los mismos loci";
    const TWICE: &str = "un locus aparece dos veces en un recuento de linaje";
    let mut cases = vec![
        case(
            "catálogo aleatorio de la amplificación sin linaje",
            |d| set(d, "/artifacts/pcr/rng_catalog", json!("f4-technical-v5")),
            "rng_catalog incompatible",
        ),
        case(
            "sin los resbalones por locus",
            |d| del(d, "/artifacts/pcr/lineage_slips_by_locus"),
            LINEAGE_TOGETHER,
        ),
        case(
            "sin los resbalones heredados",
            |d| del(d, "/artifacts/pcr/lineage_inherited_slips_by_locus"),
            LINEAGE_TOGETHER,
        ),
        case(
            "sin los pasos aplicados",
            |d| del(d, "/artifacts/pcr/lineage_applied_steps_by_locus"),
            LINEAGE_TOGETHER,
        ),
        case(
            "sin el bloque de tartamudeo clonal",
            |d| {
                set(d, "/rng_catalog", json!("f4-technical-v5"));
                del(d, "/artifacts/pcr/stutter");
            },
            LINEAGE_TOGETHER,
        ),
        case(
            "parejas que son texto",
            |d| set(d, "/artifacts/pcr/pairs", json!("400")),
            "pairs no es entero sin signo",
        ),
        case(
            "los heredados nombran menos loci",
            |d| del(d, "/artifacts/pcr/lineage_inherited_slips_by_locus/1"),
            SAME_LOCI,
        ),
        case(
            "los pasos aplicados nombran menos loci",
            |d| del(d, "/artifacts/pcr/lineage_applied_steps_by_locus/1"),
            SAME_LOCI,
        ),
        case(
            "los resbalones nombran un locus de más",
            |d| {
                push(
                    d,
                    "/artifacts/pcr/lineage_slips_by_locus",
                    json!({"locus_id": "STR-ZZZ", "slips": 0}),
                );
            },
            SAME_LOCI,
        ),
        case(
            "más heredados que resbalones",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/lineage_inherited_slips_by_locus/0/slips",
                    json!(11),
                )
            },
            "hay más resbalones heredados que resbalones en un locus",
        ),
        case(
            "más pasos escritos que parejas",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/lineage_applied_steps_by_locus/0/steps",
                    json!(401),
                )
            },
            "hay más pasos de linaje escritos que parejas en un locus",
        ),
        case(
            "locus repetido en los resbalones",
            |d| {
                let first = d["artifacts"]["pcr"]["lineage_slips_by_locus"][0].clone();
                push(d, "/artifacts/pcr/lineage_slips_by_locus", first);
            },
            TWICE,
        ),
        case(
            "locus repetido en los heredados",
            |d| {
                let first = d["artifacts"]["pcr"]["lineage_inherited_slips_by_locus"][0].clone();
                push(d, "/artifacts/pcr/lineage_inherited_slips_by_locus", first);
            },
            TWICE,
        ),
        case(
            "locus repetido en los pasos aplicados",
            |d| {
                let first = d["artifacts"]["pcr"]["lineage_applied_steps_by_locus"][0].clone();
                push(d, "/artifacts/pcr/lineage_applied_steps_by_locus", first);
            },
            TWICE,
        ),
        case(
            "los pasos aplicados con el nombre de contador de los resbalones",
            |d| {
                let steps =
                    d["artifacts"]["pcr"]["lineage_applied_steps_by_locus"][0]["steps"].clone();
                del(d, "/artifacts/pcr/lineage_applied_steps_by_locus/0/steps");
                set(
                    d,
                    "/artifacts/pcr/lineage_applied_steps_by_locus/0/slips",
                    steps,
                );
            },
            "campos inválidos en artifacts.pcr.lineage_applied_steps_by_locus[]",
        ),
        case(
            "los resbalones con el nombre de contador de los pasos",
            |d| {
                let slips = d["artifacts"]["pcr"]["lineage_slips_by_locus"][0]["slips"].clone();
                del(d, "/artifacts/pcr/lineage_slips_by_locus/0/slips");
                set(d, "/artifacts/pcr/lineage_slips_by_locus/0/steps", slips);
            },
            "campos inválidos en artifacts.pcr.lineage_slips_by_locus[]",
        ),
        case(
            "locus que no es texto",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/lineage_slips_by_locus/0/locus_id",
                    json!(7),
                )
            },
            "locus_id debe ser texto no vacío",
        ),
        case(
            "resbalones que son texto",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/lineage_slips_by_locus/0/slips",
                    json!("10"),
                )
            },
            "slips no es entero sin signo",
        ),
        case(
            "pasos aplicados que son texto",
            |d| {
                set(
                    d,
                    "/artifacts/pcr/lineage_applied_steps_by_locus/0/steps",
                    json!("10"),
                )
            },
            "steps no es entero sin signo",
        ),
    ];
    cases.shrink_to_fit();
    ledger.rejected("linaje/clonal", &clonal, cases);
    ledger.accepted(
        "linaje/clonal (lado aceptado)",
        &clonal,
        vec![
            accepted("tantos heredados como resbalones", |d| {
                set(
                    d,
                    "/artifacts/pcr/lineage_inherited_slips_by_locus/0/slips",
                    json!(10),
                );
            }),
            accepted("tantos pasos escritos como parejas", |d| {
                set(
                    d,
                    "/artifacts/pcr/lineage_applied_steps_by_locus/0/steps",
                    json!(400),
                );
            }),
        ],
    );
    ledger.finish();
}

// ---------------------------------------------------------------------------
// Contaminación biológica
// ---------------------------------------------------------------------------

const CONTAMINATION: &str = "/artifacts/biological_contamination/0";
const CONSERVATION: &str = "los contadores realizados de contaminación no conservan las parejas";
const PERCENT: &str = "porcentaje de contaminación biológica inválido";
const DONOR: &str = "aporte de donante inválido";

fn contamination(field: &str) -> String {
    format!("{CONTAMINATION}/{field}")
}

/// **El registro de contaminación biológica de la receptora.**
///
/// El escenario trae un registro real: 3 parejas limpias, un 50 %, objetivo 2,
/// una retenida, un donante con 2 aportes. Las cinco igualdades de
/// conservación son **mutuamente redundantes** (cada una se deduce de las
/// otras cuatro), de modo que ninguna mutación **sola** de una de ellas puede
/// morir; lo que sí tiene que morir es quitarlas **por pares**, y para eso hay
/// un caso que viola cada par posible sin violar las otras tres.
#[test]
fn the_biological_contamination_record_is_judged_exactly() {
    let mut ledger = Ledger::default();
    let plain = plain();
    let mut cases = vec![
        case(
            "receptora con identificador inválido",
            |d| set(d, &contamination("recipient_id"), json!("no valido!")),
            "ID individual inválido",
        ),
        case(
            "biblioteca de la receptora con identificador inválido",
            |d| {
                set(
                    d,
                    &contamination("recipient_library_id"),
                    json!("no valido!"),
                )
            },
            "ID library inválido",
        ),
        case(
            "receptora repetida",
            |d| {
                let first = d["artifacts"]["biological_contamination"][0].clone();
                push(d, "/artifacts/biological_contamination", first);
            },
            "receptora de contaminación biológica repetida",
        ),
        case(
            "mínimo por encima del porcentaje muestreado",
            |d| set(d, &contamination("percent_range/minimum"), json!(60.0)),
            PERCENT,
        ),
        case(
            "porcentaje muestreado por encima del máximo",
            |d| set(d, &contamination("percent_range/maximum"), json!(40.0)),
            PERCENT,
        ),
        case(
            "mínimo negativo",
            |d| set(d, &contamination("percent_range/minimum"), json!(-1.0)),
            PERCENT,
        ),
        case(
            "máximo por encima de cien",
            |d| set(d, &contamination("percent_range/maximum"), json!(101.0)),
            PERCENT,
        ),
        case(
            "mínimo que es texto",
            |d| set(d, &contamination("percent_range/minimum"), json!("50")),
            PERCENT,
        ),
        case(
            "máximo que es texto",
            |d| set(d, &contamination("percent_range/maximum"), json!("50")),
            PERCENT,
        ),
        case(
            "porcentaje muestreado que es texto",
            |d| set(d, &contamination("sampled_percent"), json!("50")),
            PERCENT,
        ),
        case(
            "parejas limpias que son texto",
            |d| set(d, &contamination("clean_pair_count"), json!("3")),
            "recuento limpio de contaminación inválido",
        ),
        case(
            "objetivo que es texto",
            |d| set(d, &contamination("target_pair_count"), json!("2")),
            "objetivo de contaminación inválido",
        ),
        case(
            "retiradas que son texto",
            |d| set(d, &contamination("removed_own_pair_count"), json!("2")),
            "recuento retirado de contaminación inválido",
        ),
        case(
            "insertadas que son texto",
            |d| {
                set(
                    d,
                    &contamination("inserted_contaminant_pair_count"),
                    json!("2"),
                )
            },
            "recuento insertado de contaminación inválido",
        ),
        case(
            "retenidas que son texto",
            |d| set(d, &contamination("retained_own_pair_count"), json!("1")),
            "recuento retenido de contaminación inválido",
        ),
        case(
            "finales que son texto",
            |d| set(d, &contamination("final_pair_count"), json!("3")),
            "recuento final de contaminación inválido",
        ),
        case(
            "objetivo por encima de las parejas limpias",
            |d| set(d, &contamination("target_pair_count"), json!(4)),
            "objetivo de contaminación supera parejas limpias",
        ),
        case(
            "objetivo que no se puede recomputar",
            |d| {
                set(
                    d,
                    &contamination("clean_pair_count"),
                    json!(16_487_797_971_023_415_208_u64),
                );
                set(d, &contamination("target_pair_count"), json!(0));
                for field in ["minimum", "maximum"] {
                    set(
                        d,
                        &contamination(&format!("percent_range/{field}")),
                        json!(100.0),
                    );
                }
                set(d, &contamination("sampled_percent"), json!(100.0));
            },
            "objetivo de contaminación no es recomputable",
        ),
        case(
            "objetivo que no aplica el redondeo normativo",
            |d| set(d, &contamination("target_pair_count"), json!(1)),
            "objetivo de contaminación no aplica el redondeo normativo",
        ),
        case(
            "retiradas y insertadas distintas del objetivo",
            |d| {
                set(d, &contamination("removed_own_pair_count"), json!(3));
                set(
                    d,
                    &contamination("inserted_contaminant_pair_count"),
                    json!(3),
                );
                set(d, &contamination("retained_own_pair_count"), json!(0));
            },
            CONSERVATION,
        ),
        case(
            "retenidas de más",
            |d| set(d, &contamination("retained_own_pair_count"), json!(2)),
            CONSERVATION,
        ),
        case(
            "finales de más",
            |d| set(d, &contamination("final_pair_count"), json!(4)),
            CONSERVATION,
        ),
        case(
            "limpias y finales de más con el resto coherente",
            |d| {
                set(d, &contamination("clean_pair_count"), json!(4));
                set(d, &contamination("final_pair_count"), json!(3));
            },
            CONSERVATION,
        ),
        case(
            "retiradas distintas del objetivo y finales que no suman lo retenido",
            |d| {
                set(d, &contamination("removed_own_pair_count"), json!(3));
                set(d, &contamination("retained_own_pair_count"), json!(0));
            },
            CONSERVATION,
        ),
        case(
            "insertadas distintas del objetivo y retiradas que no suman lo limpio",
            |d| {
                set(
                    d,
                    &contamination("inserted_contaminant_pair_count"),
                    json!(3),
                );
                set(d, &contamination("retained_own_pair_count"), json!(0));
            },
            CONSERVATION,
        ),
        case(
            "retiradas distintas del objetivo y finales distintos de los limpios",
            |d| {
                set(d, &contamination("removed_own_pair_count"), json!(3));
                set(d, &contamination("retained_own_pair_count"), json!(0));
                set(d, &contamination("final_pair_count"), json!(2));
            },
            CONSERVATION,
        ),
        case(
            "insertadas distintas del objetivo y finales distintos de los limpios",
            |d| {
                set(
                    d,
                    &contamination("inserted_contaminant_pair_count"),
                    json!(3),
                );
                set(d, &contamination("final_pair_count"), json!(4));
            },
            CONSERVATION,
        ),
        case(
            "solo las retiradas distintas del objetivo",
            |d| set(d, &contamination("removed_own_pair_count"), json!(3)),
            CONSERVATION,
        ),
        case(
            "solo las insertadas distintas del objetivo",
            |d| {
                set(
                    d,
                    &contamination("inserted_contaminant_pair_count"),
                    json!(3),
                )
            },
            CONSERVATION,
        ),
        case(
            "retenidas que desbordan al sumar",
            |d| {
                set(
                    d,
                    &contamination("retained_own_pair_count"),
                    json!(u64::MAX),
                )
            },
            CONSERVATION,
        ),
        case(
            "aportes que no son lista",
            |d| set(d, &contamination("donor_contributions"), json!("x")),
            "donor_contributions no es lista",
        ),
        case(
            "donante con identificador inválido",
            |d| {
                set(
                    d,
                    &contamination("donor_contributions/0/individual_id"),
                    json!("no valido!"),
                )
            },
            "ID individual inválido",
        ),
        case(
            "aporte que es texto",
            |d| set(d, &contamination("donor_contributions/0/pairs"), json!("2")),
            DONOR,
        ),
        case(
            "la receptora como donante",
            |d| {
                set(
                    d,
                    &contamination("donor_contributions/0/individual_id"),
                    json!("sample-A"),
                )
            },
            DONOR,
        ),
        case(
            "aporte vacío",
            |d| set(d, &contamination("donor_contributions/0/pairs"), json!(0)),
            DONOR,
        ),
        case(
            "donante repetido",
            |d| {
                set(
                    d,
                    &contamination("donor_contributions"),
                    json!([
                        {"individual_id": "sample-B", "pairs": 1},
                        {"individual_id": "sample-B", "pairs": 1},
                    ]),
                );
            },
            DONOR,
        ),
        case(
            "aportes que desbordan al sumar",
            |d| {
                set(
                    d,
                    &contamination("donor_contributions"),
                    json!([
                        {"individual_id": "sample-B", "pairs": u64::MAX},
                        {"individual_id": "sample-C", "pairs": 1},
                    ]),
                );
            },
            "desbordamiento de aportes",
        ),
        case(
            "aportes que no suman el objetivo",
            |d| set(d, &contamination("donor_contributions/0/pairs"), json!(1)),
            "los aportes no suman el objetivo de contaminación",
        ),
    ];
    for name in ["percent", "donors", "removal"] {
        for (label, mutate) in [
            (
                "de 63 caracteres",
                (|seed: String| seed[1..].to_owned()) as fn(String) -> String,
            ),
            ("con un carácter que no es hexadecimal", |seed| {
                format!("g{}", &seed[1..])
            }),
            ("de 65 caracteres", |seed| format!("{seed}a")),
            ("en mayúsculas", |seed| seed.to_uppercase()),
        ] {
            cases.push(case(
                format!("subsemilla {name} {label}"),
                move |d| {
                    let pointer = contamination(&format!("subseeds/{name}"));
                    let altered = mutate(text(d, &pointer));
                    set(d, &pointer, json!(altered));
                },
                format!("{name} no es una subsemilla hexadecimal de 64 caracteres"),
            ));
        }
    }
    ledger.rejected("contaminación/plain", &plain, cases);
    ledger.accepted(
        "contaminación/plain (lado aceptado)",
        &plain,
        vec![
            accepted("mínimo en cero", |d| {
                set(d, &contamination("percent_range/minimum"), json!(0.0));
            }),
            accepted("máximo en cien", |d| {
                set(d, &contamination("percent_range/maximum"), json!(100.0));
            }),
            accepted("dos donantes distintos", |d| {
                set(
                    d,
                    &contamination("donor_contributions"),
                    json!([
                        {"individual_id": "sample-B", "pairs": 1},
                        {"individual_id": "sample-C", "pairs": 1},
                    ]),
                );
            }),
            accepted("contaminación del cien por cien", |d| {
                // El límite inclusivo de `target > clean`: objetivo igual a las
                // parejas limpias, y el intervalo cerrado en cien.
                for field in ["minimum", "maximum"] {
                    set(
                        d,
                        &contamination(&format!("percent_range/{field}")),
                        json!(100.0),
                    );
                }
                set(d, &contamination("sampled_percent"), json!(100.0));
                for field in [
                    "target_pair_count",
                    "removed_own_pair_count",
                    "inserted_contaminant_pair_count",
                ] {
                    set(d, &contamination(field), json!(3));
                }
                set(d, &contamination("retained_own_pair_count"), json!(0));
                set(d, &contamination("donor_contributions/0/pairs"), json!(3));
            }),
            accepted("sin ninguna receptora", |d| {
                set(d, "/artifacts/biological_contamination", json!([]));
            }),
        ],
    );
    ledger.finish();
}

// ---------------------------------------------------------------------------
// Resultados FASTQ, estado del evento y diagnósticos
// ---------------------------------------------------------------------------

const OUTPUTS_ORDER: &str = "artefactos FASTQ no son únicos, completos por biblioteca u ordenados";
const OUTPUT_INCOHERENT: &str = "artefacto FASTQ completo incoherente";
const NOT_COMPLETE_COMPLETE: &str =
    "un evento no completo exige diagnóstico y un FASTQ no completo";
const DIAGNOSTIC_REJECTED: &str = "diagnóstico desconocido, alterado o repetido";

/// **Las salidas FASTQ.** Un caso por mitad de la guarda de unicidad y orden
/// y por cada una de las siete mitades de la coherencia de un artefacto
/// completo.
#[test]
fn the_fastq_outputs_are_judged_exactly() {
    let mut ledger = Ledger::default();
    let plain = plain();
    ledger.rejected(
        "salidas/plain",
        &plain,
        vec![
            case(
                "ninguna salida",
                |d| set(d, "/outputs", json!([])),
                "el evento no tiene resultados FASTQ",
            ),
            case(
                "salida que no es un FASTQ",
                |d| set(d, "/outputs/0/kind", json!("reference_fasta")),
                "F3 solo admite artefactos FASTQ",
            ),
            case(
                "estado de una salida inválido",
                |d| set(d, "/outputs/0/state", json!("quizá")),
                "estado inválido",
            ),
            case(
                "mate tres",
                |d| set(d, "/outputs/0/mate", json!(3)),
                "mate FASTQ inválido",
            ),
            case(
                "mate cero",
                |d| set(d, "/outputs/0/mate", json!(0)),
                "mate FASTQ inválido",
            ),
            case(
                "mate que es texto",
                |d| set(d, "/outputs/0/mate", json!("1")),
                "mate FASTQ inválido",
            ),
            case(
                "salida de un individuo desconocido",
                |d| set(d, "/outputs/0/individual_id", json!("sample-Z")),
                OUTPUTS_ORDER,
            ),
            case(
                "salida con la biblioteca de otro",
                |d| {
                    let other = text(d, "/outputs/2/library_id");
                    set(d, "/outputs/0/library_id", json!(other));
                },
                OUTPUTS_ORDER,
            ),
            case(
                "identificador de salida repetido",
                |d| {
                    let first = text(d, "/outputs/0/id");
                    set(d, "/outputs/1/id", json!(first));
                },
                OUTPUTS_ORDER,
            ),
            case(
                "mate repetido para un individuo",
                |d| set(d, "/outputs/1/mate", json!(1)),
                OUTPUTS_ORDER,
            ),
            case(
                "salidas fuera de orden",
                |d| swap_items(d, "/outputs", 0, 1),
                OUTPUTS_ORDER,
            ),
            case(
                "falta el resultado de un mate",
                |d| del(d, "/outputs/3"),
                "falta el resultado FASTQ de una biblioteca seleccionada",
            ),
            case(
                "artefacto de una salida completa nulo",
                |d| set(d, "/outputs/0/artifact", Value::Null),
                "artifact no puede ser nulo",
            ),
            case(
                "artefacto de una salida completa que es texto",
                |d| set(d, "/outputs/0/artifact", json!("x")),
                "artifact no puede ser nulo",
            ),
            case(
                "identificador del artefacto distinto del de la salida",
                |d| set(d, "/outputs/0/artifact/id", json!("otro")),
                OUTPUT_INCOHERENT,
            ),
            case(
                "clase del artefacto distinta de la de la salida",
                |d| set(d, "/outputs/0/artifact/kind", json!("otro")),
                OUTPUT_INCOHERENT,
            ),
            case(
                "ruta del artefacto que no es la canónica",
                |d| {
                    set(
                        d,
                        "/outputs/0/artifact/path",
                        json!("outputs/fastq/reads-main/otro_R1.fastq.gz"),
                    );
                },
                OUTPUT_INCOHERENT,
            ),
            case(
                "artefacto con el identificador de la referencia",
                |d| {
                    let id = text(d, "/outputs/0/artifact/id");
                    set(d, "/reference/fasta/id", json!(id));
                },
                OUTPUT_INCOHERENT,
            ),
            case(
                "artefacto con la ruta de la referencia",
                |d| {
                    let path = text(d, "/outputs/0/artifact/path");
                    set(d, "/reference/fasta/path", json!(path));
                },
                OUTPUT_INCOHERENT,
            ),
            case(
                "salida completa con código de diagnóstico",
                |d| set(d, "/outputs/0/diagnostic_code", json!("GPS-FMT-006")),
                OUTPUT_INCOHERENT,
            ),
            case(
                "evento completo con un diagnóstico",
                |d| set(d, "/diagnostics", json!([format_diagnostic()])),
                "un evento completo exige FASTQ completos y ningún diagnóstico",
            ),
        ],
    );
    ledger.finish();
}

/// **El estado del evento, sus salidas no completas y sus diagnósticos.**
#[test]
fn the_event_state_the_incomplete_outputs_and_the_diagnostics_are_judged_exactly() {
    let mut ledger = Ledger::default();
    let plain = plain();
    let failed = failed_of(&plain);
    ledger.rejected(
        "estado/failed",
        &failed,
        vec![
            case(
                "salida no completa con identidad física",
                |d| {
                    let record = d["reference"]["fasta"].clone();
                    set(d, "/outputs/0/artifact", record);
                },
                "artefacto no completo no puede tener identidad física",
            ),
            case(
                "salida no completa sin código de diagnóstico",
                |d| set(d, "/outputs/0/diagnostic_code", Value::Null),
                "diagnostic_code debe ser texto no vacío",
            ),
            case(
                "salida no completa con un código que no está en los diagnósticos",
                |d| set(d, "/outputs/0/diagnostic_code", json!("GPS-IO-001")),
                "diagnóstico de artefacto ausente",
            ),
            case(
                "evento fallido con todas las salidas completas",
                |d| {
                    let complete = d["outputs"][1].clone();
                    let mut first = complete;
                    first["mate"] = json!(1);
                    first["id"] = d["outputs"][0]["id"].clone();
                    first["artifact"]["id"] = d["outputs"][0]["id"].clone();
                    first["artifact"]["path"] =
                        json!("outputs/fastq/reads-main/sample-A_R1.fastq.gz");
                    set(d, "/outputs/0", first);
                },
                NOT_COMPLETE_COMPLETE,
            ),
            case(
                "evento no completo con cobertura",
                |d| {
                    let coverage = plain_coverage();
                    set(d, "/coverage", coverage);
                },
                "un evento no completo no publica cobertura F4",
            ),
            case(
                "código de diagnóstico desconocido",
                |d| {
                    set(d, "/diagnostics/0/code", json!("GPS-BOGUS-001"));
                    set(d, "/outputs/0/diagnostic_code", json!("GPS-BOGUS-001"));
                },
                DIAGNOSTIC_REJECTED,
            ),
            case(
                "mensaje de un diagnóstico alterado",
                |d| set(d, "/diagnostics/0/message", json!("otro mensaje")),
                DIAGNOSTIC_REJECTED,
            ),
            case(
                "diagnósticos en orden inverso",
                |d| {
                    push(
                        d,
                        "/diagnostics",
                        json!({
                            "code": "GPS-IO-001",
                            "message": "Una operación de lectura o escritura falló.",
                        }),
                    );
                    swap_items(d, "/diagnostics", 0, 1);
                },
                DIAGNOSTIC_REJECTED,
            ),
            case(
                "diagnóstico repetido",
                |d| push(d, "/diagnostics", format_diagnostic()),
                DIAGNOSTIC_REJECTED,
            ),
            case(
                "código de diagnóstico que no es texto",
                |d| set(d, "/diagnostics/0/code", json!(7)),
                "code debe ser texto no vacío",
            ),
            case(
                "mensaje de diagnóstico vacío",
                |d| set(d, "/diagnostics/0/message", json!("")),
                "message debe ser texto no vacío",
            ),
        ],
    );
    ledger.accepted(
        "estado/failed (lado aceptado)",
        &failed,
        vec![
            accepted("evento incompleto en vez de fallido", |d| {
                set(d, "/state", json!("incomplete"));
            }),
            accepted("salida incompleta en vez de fallida", |d| {
                set(d, "/outputs/0/state", json!("incomplete"));
            }),
            accepted("dos diagnósticos en orden", |d| {
                push(
                    d,
                    "/diagnostics",
                    json!({
                        "code": "GPS-IO-001",
                        "message": "Una operación de lectura o escritura falló.",
                    }),
                );
            }),
        ],
    );
    ledger.finish();
}

/// La cobertura del escenario completo, para pegársela a un evento fallido.
fn plain_coverage() -> Value {
    plain()["coverage"].clone()
}

// ---------------------------------------------------------------------------
// Tartamudeo
// ---------------------------------------------------------------------------

const REFERENCE_TRACT: &str = "una longitud de tracto de referencia no es la del algoritmo";
const MEASURED_RANGE: &str = "el tramo medido publicado no es el del algoritmo";
const CURVE_FACTOR: &str = "un factor de la curva no es el del algoritmo";
const STEP_WEIGHTS: &str = "los pesos publicados no cubren exactamente la estructura del algoritmo";

/// **La sección del tartamudeo.** El escenario de tartamudeo sin amplificar
/// trae la sección real, con 40 tractos evaluados y 20 que tartamudean; la
/// curva de dependencia con la longitud del tracto se compara contra el
/// algoritmo, fila por fila.
#[test]
fn the_stutter_section_is_judged_exactly() {
    let mut ledger = Ledger::default();
    let stutter = stutter_only();
    let cases = vec![
        case(
            "sección que se publica desactivada",
            |d| set(d, "/artifacts/stutter/enabled", json!(false)),
            "un manifiesto que publica la sección de tartamudeo la publica activa",
        ),
        case(
            "marca de activación que es texto",
            |d| set(d, "/artifacts/stutter/enabled", json!("true")),
            "un manifiesto que publica la sección de tartamudeo la publica activa",
        ),
        case(
            "algoritmo del tartamudeo",
            |d| set(d, "/artifacts/stutter/algorithm", json!("otro")),
            "algorithm incompatible",
        ),
        case(
            "catálogo aleatorio del tartamudeo",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/rng_catalog",
                    json!("f4-technical-v5"),
                )
            },
            "rng_catalog incompatible",
        ),
        case(
            "coeficiente de la tirada publicado como calibrado",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/run_length_coefficient/calibrated",
                    json!(true),
                )
            },
            "el coeficiente de la tirada tiene que publicarse como no calibrado",
        ),
        case(
            "marca de calibración que es texto",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/run_length_coefficient/calibrated",
                    json!("no"),
                )
            },
            "el coeficiente de la tirada tiene que publicarse como no calibrado",
        ),
        case(
            "discrepancia conocida de la calibración",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/calibration/known_mismatch",
                    json!("otra"),
                )
            },
            "known_mismatch incompatible",
        ),
        case(
            "procedencia de la curva",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/source",
                    json!("otra"),
                )
            },
            "source incompatible",
        ),
        case(
            "motivo de la tosquedad de la curva",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/coarse_reason",
                    json!("otro"),
                )
            },
            "coarse_reason incompatible",
        ),
        case(
            "curva publicada como no tosca",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/coarse",
                    json!(false),
                )
            },
            "la dependencia con la longitud del tracto se publica marcada como tosca",
        ),
        case(
            "curva extrapolada",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/extrapolated",
                    json!(true),
                )
            },
            "fuera del tramo medido se mantiene el extremo y no se extrapola",
        ),
        case(
            "marca de extrapolación que es texto",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/extrapolated",
                    json!("no"),
                )
            },
            "fuera del tramo medido se mantiene el extremo y no se extrapola",
        ),
        case(
            "tabla de longitudes con una fila de menos",
            |d| {
                del(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/5",
                )
            },
            "la tabla de longitudes de referencia no cubre el tramo calibrado",
        ),
        case(
            "tabla de longitudes con una fila de más",
            |d| {
                let extra = d["artifacts"]["stutter"]["tract_length_dependence"]["by_motif_length"]
                    [5]
                .clone();
                push(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length",
                    extra,
                );
            },
            "la tabla de longitudes de referencia no cubre el tramo calibrado",
        ),
        case(
            "longitud de referencia de un motivo medido cambiada",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/reference_tract_bases",
                    json!(13),
                );
            },
            REFERENCE_TRACT,
        ),
        case(
            "longitud de referencia de un motivo medido nula",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/reference_tract_bases",
                    Value::Null,
                );
            },
            REFERENCE_TRACT,
        ),
        case(
            "longitud de referencia de un motivo medido que es texto",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/reference_tract_bases",
                    json!("15"),
                );
            },
            REFERENCE_TRACT,
        ),
        case(
            "longitud de referencia de un motivo no medido",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/2/reference_tract_bases",
                    json!(5),
                );
            },
            REFERENCE_TRACT,
        ),
        case(
            "motivo medido publicado como no medido",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/measured",
                    json!(false),
                );
            },
            REFERENCE_TRACT,
        ),
        case(
            "motivo no medido publicado como medido",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/2/measured",
                    json!(true),
                );
            },
            REFERENCE_TRACT,
        ),
        case(
            "mínimo del tramo medido cambiado",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/measured_units_min",
                    json!(6),
                );
            },
            MEASURED_RANGE,
        ),
        case(
            "máximo del tramo medido cambiado",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/measured_units_max",
                    json!(14),
                );
            },
            MEASURED_RANGE,
        ),
        case(
            "mínimo del tramo medido nulo",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/measured_units_min",
                    Value::Null,
                );
            },
            MEASURED_RANGE,
        ),
        case(
            "tramo medido en un motivo que no lo tiene (mínimo)",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/2/measured_units_min",
                    json!(5),
                );
            },
            MEASURED_RANGE,
        ),
        case(
            "tramo medido en un motivo que no lo tiene (máximo)",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/2/measured_units_max",
                    json!(5),
                );
            },
            MEASURED_RANGE,
        ),
        case(
            "tabla de factores con una fila de menos",
            |d| {
                let pointer =
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/factors";
                let last = d
                    .pointer(pointer)
                    .and_then(Value::as_array)
                    .expect("lista")
                    .len()
                    - 1;
                del(d, &format!("{pointer}/{last}"));
            },
            "la tabla de factores no cubre el tramo medido",
        ),
        case(
            "tabla de factores en un motivo que no tiene tramo",
            |d| {
                push(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/2/factors",
                    json!({"units": 5, "factor": 1.0}),
                );
            },
            "la tabla de factores no cubre el tramo medido",
        ),
        case(
            "unidades de un factor cambiadas",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/factors/0/units",
                    json!(6),
                );
            },
            CURVE_FACTOR,
        ),
        case(
            "unidades de un factor posterior cambiadas",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/1/factors/2/units",
                    json!(9),
                );
            },
            CURVE_FACTOR,
        ),
        case(
            "valor de un factor cambiado",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/0/factors/3/factor",
                    json!(0.5),
                );
            },
            CURVE_FACTOR,
        ),
        case(
            "valor del factor de otro motivo cambiado",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/tract_length_dependence/by_motif_length/1/factors/0/factor",
                    json!(0.5),
                );
            },
            CURVE_FACTOR,
        ),
        case(
            "paso publicado sin unidades",
            |d| set(d, "/artifacts/stutter/step_weights/0/units", json!("x")),
            "un paso publicado no declara sus unidades",
        ),
        case(
            "pesos con un paso de menos",
            |d| del(d, "/artifacts/stutter/step_weights/3"),
            STEP_WEIGHTS,
        ),
        case(
            "pesos con los pasos en otro orden",
            |d| swap_items(d, "/artifacts/stutter/step_weights", 0, 1),
            STEP_WEIGHTS,
        ),
        case(
            "pesos con un paso de más",
            |d| {
                push(
                    d,
                    "/artifacts/stutter/step_weights",
                    json!({"units": 2, "weight": 0.1}),
                )
            },
            STEP_WEIGHTS,
        ),
        case(
            "más tartamudeos que tractos evaluados",
            |d| set(d, "/artifacts/stutter/stuttered_tracts", json!(41)),
            "no pueden tartamudear más moléculas de las evaluadas",
        ),
        case(
            "tractos evaluados que son texto",
            |d| set(d, "/artifacts/stutter/evaluated_tracts", json!("40")),
            "evaluated_tracts no es entero sin signo",
        ),
        case(
            "tractos que tartamudean que son texto",
            |d| set(d, "/artifacts/stutter/stuttered_tracts", json!("20")),
            "stuttered_tracts no es entero sin signo",
        ),
        case(
            "paso realizado sin unidades",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/steps_by_locus/0/steps/0/units",
                    json!("x"),
                )
            },
            "un paso realizado no declara sus unidades",
        ),
        case(
            "paso realizado que no es de la estructura",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/steps_by_locus/0/steps/0/units",
                    json!(2),
                )
            },
            "el paso realizado 2 no pertenece a la estructura",
        ),
        case(
            "paso realizado negativo que no es de la estructura",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/steps_by_locus/0/steps/0/units",
                    json!(-4),
                )
            },
            "el paso realizado -4 no pertenece a la estructura",
        ),
        case(
            "recuento de un paso que es texto",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/steps_by_locus/0/steps/0/count",
                    json!("1"),
                )
            },
            "count no es entero sin signo",
        ),
        case(
            "distribución de pasos que desborda",
            |d| {
                set(
                    d,
                    "/artifacts/stutter/steps_by_locus/0/steps",
                    json!([{"units": -3, "count": u64::MAX}, {"units": -2, "count": 1}]),
                );
            },
            "la distribución de pasos desborda",
        ),
        case(
            "pasos que suman más que los tartamudeos",
            |d| set(d, "/artifacts/stutter/stuttered_tracts", json!(19)),
            "la distribución de pasos suma 20 y las moléculas que tartamudearon son 19",
        ),
        case(
            "pasos que suman menos que los tartamudeos",
            |d| set(d, "/artifacts/stutter/stuttered_tracts", json!(21)),
            "la distribución de pasos suma 20 y las moléculas que tartamudearon son 21",
        ),
    ];
    ledger.rejected("tartamudeo/stutter", &stutter, cases);
    ledger.accepted(
        "tartamudeo/stutter (lado aceptado)",
        &stutter,
        vec![
            accepted("tantos tartamudeos como tractos evaluados", |d| {
                set(d, "/artifacts/stutter/evaluated_tracts", json!(20));
            }),
            accepted("ningún tartamudeo", |d| {
                set(d, "/artifacts/stutter/stuttered_tracts", json!(0));
                set(d, "/artifacts/stutter/steps_by_locus", json!([]));
            }),
        ],
    );
    // La misma sección dentro de un documento amplificado y con linaje.
    let clonal = clonal();
    ledger.rejected(
        "tartamudeo/clonal",
        &clonal,
        vec![
            case(
                "sección publicada desactivada",
                |d| set(d, "/artifacts/stutter/enabled", json!(false)),
                "un manifiesto que publica la sección de tartamudeo la publica activa",
            ),
            case(
                "pasos que no suman los tartamudeos",
                |d| {
                    let stuttered = number(d, "/artifacts/stutter/stuttered_tracts");
                    set(d, "/artifacts/stutter/stuttered_tracts", json!(stuttered + 1));
                },
                format!(
                    "la distribución de pasos suma {total} y las moléculas que tartamudearon son {next}",
                    total = number(&clonal, "/artifacts/stutter/stuttered_tracts"),
                    next = number(&clonal, "/artifacts/stutter/stuttered_tracts") + 1
                )
                .as_str(),
            ),
        ],
    );
    ledger.finish();
}

// ---------------------------------------------------------------------------
// Los hechos del escritor (`validate_writer_facts`)
// ---------------------------------------------------------------------------

/// Las piezas de un evento mínimo completo (una biblioteca, PCR-free) con las
/// que armar un `EventManifest` y estropearlo en un solo hecho.
struct Parts {
    event: &'static EffectiveTechnicalEvent,
    individual: &'static IndividualId,
    library: &'static LibraryId,
    reference_id: &'static ReferenceId,
    fasta: &'static ArtifactRecord,
    identities: &'static TechnicalIdentityPlan,
    records: [&'static ArtifactRecord; 2],
    ids: [&'static str; 2],
}

fn leak<T>(value: T) -> &'static T {
    Box::leak(Box::new(value))
}

fn record(id: &str, kind: &str, path: &str, digest: char) -> ArtifactRecord {
    ArtifactRecord {
        id: id.to_owned(),
        kind: kind.to_owned(),
        path: path.to_owned(),
        size_bytes: 1,
        digest: digest.to_string().repeat(64),
    }
}

fn plan(
    event_id: &str,
    reference: &ReferenceId,
    individuals: &[IndividualId],
) -> TechnicalIdentityPlan {
    derive_technical_identity_plan(
        RootSeed::parse("root", SEED).expect("semilla"),
        reference,
        EventId::parse(event_id, "test").expect("evento"),
        individuals,
    )
    .expect("identidades técnicas")
}

fn parts() -> Parts {
    let text = b"config_kind = \"general\"\nschema_version = \"2\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1000\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\"]\n\n[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 30\n\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n";
    let declared = parse_declared_config(
        admit_toml_header(
            Path::new("writer-facts.toml"),
            text.to_vec(),
            ExpectedConfigKind::General,
        )
        .expect("cabecera"),
    )
    .expect("configuración declarada");
    let effective = leak(build_effective_config(&declared).expect("efectiva"));
    let event = effective.technical_event().expect("evento técnico");
    let reference_id = leak(ReferenceId::parse("reference-main", "test").expect("referencia"));
    let identities = leak(plan("reads-main", reference_id, event.individual_ids()));
    let individual = event.individual_ids().first().expect("individuo");
    let library = identities.library(individual).expect("biblioteca");
    let first: &'static str =
        Box::leak(format!("fastq-{}", library.library_id().as_str()).into_boxed_str());
    let second: &'static str = Box::leak(format!("{first}-r2").into_boxed_str());
    Parts {
        event,
        individual,
        library: library.library_id(),
        reference_id,
        fasta: leak(record(
            "reference-main-fasta",
            "reference_fasta",
            "resources/references/reference-main.fa",
            'a',
        )),
        identities,
        records: [
            leak(record(
                first,
                "fastq",
                "outputs/fastq/reads-main/sample-A_R1.fastq.gz",
                'b',
            )),
            leak(record(
                second,
                "fastq",
                "outputs/fastq/reads-main/sample-A_R2.fastq.gz",
                'a',
            )),
        ],
        ids: [first, second],
    }
}

fn writer_artifact(parts: &Parts, mate: usize) -> EventManifestArtifact<'static> {
    EventManifestArtifact {
        id: parts.ids[mate],
        kind: "fastq",
        state: UnitState::Complete,
        individual_id: parts.individual,
        library_id: parts.library,
        mate: if mate == 0 {
            TechnicalMate::R1
        } else {
            TechnicalMate::R2
        },
        artifact: Some(parts.records[mate]),
        diagnostic: None,
    }
}

fn writer_coverage(parts: &Parts) -> EventManifestCoverage<'static> {
    EventManifestCoverage {
        individual_id: parts.individual,
        library_id: parts.library,
        report: F4CoverageReport::from_base_counts(30.0, 30_000, 30_000, 0, 1_000, 0)
            .expect("cobertura"),
    }
}

fn writer_summary(parts: &Parts) -> EventManifestFinalFastqSummary<'static> {
    EventManifestFinalFastqSummary {
        individual_id: parts.individual,
        library_id: parts.library,
        pair_count: 100,
        read_count: 200,
        base_count: 30_000,
        own_pair_count: 100,
        contaminant_pair_count: 0,
        own_base_count: 30_000,
        contaminant_base_count: 0,
        n_base_count: 0,
    }
}

/// Un resumen cuyas parejas propias y contaminantes no suman las parejas: el
/// escritor no lo mira y el validador lo rechaza antes que casi todo lo demás.
fn incoherent_summary(parts: &Parts) -> EventManifestFinalFastqSummary<'static> {
    EventManifestFinalFastqSummary {
        own_pair_count: 99,
        ..writer_summary(parts)
    }
}

/// Un diagnóstico `GPS-FMT-006` de verdad: lo da el propio juez ante un
/// documento que no es JSON (los constructores de `Diagnostic` no son
/// públicos).
fn writer_diagnostic() -> Diagnostic {
    let error = validate_event_manifest_bytes(Path::new(DOCUMENT_PATH), b"no es json")
        .expect_err("un documento que no es JSON se rechaza");
    assert_eq!(error.code(), "GPS-FMT-006");
    error
}

/// Un diagnóstico distinto del anterior, para el artefacto cuyo código el
/// manifiesto no lleva: el de una configuración que no existe.
fn other_diagnostic() -> Diagnostic {
    let error = load_and_validate_config(Path::new("/no/existe/general.toml"))
        .expect_err("una configuración que no existe se rechaza");
    assert_ne!(error.code(), "GPS-FMT-006");
    error
}

fn writer_manifest(parts: &Parts) -> EventManifest<'static> {
    EventManifest {
        event: parts.event,
        state: UnitState::Complete,
        identities: parts.identities,
        reference: EventManifestReference {
            reference_id: parts.reference_id,
            fasta: parts.fasta,
        },
        artifacts: leak([writer_artifact(parts, 0), writer_artifact(parts, 1)]),
        coverage: leak([writer_coverage(parts)]),
        final_fastq_summary: leak([writer_summary(parts)]),
        error_quality: EventManifestErrorQuality::neutral(),
        demultiplexing: EventManifestDemultiplexing::disabled(),
        chimera_rate: parts.event.chimera_rate(),
        realized_chimera_pairs: 0,
        biological_contamination: &[],
        stutter: None,
        pcr: None,
        diagnostics: &[],
    }
}

/// Un evento fallido: los dos FASTQ fallan con el mismo diagnóstico y no hay
/// coberturas ni resúmenes.
fn failed_manifest(parts: &Parts) -> EventManifest<'static> {
    let diagnostic = leak(writer_diagnostic());
    let failed = |mate: usize| EventManifestArtifact {
        state: UnitState::Failed,
        artifact: None,
        diagnostic: Some(diagnostic),
        ..writer_artifact(parts, mate)
    };
    EventManifest {
        state: UnitState::Failed,
        artifacts: leak([failed(0), failed(1)]),
        coverage: &[],
        final_fastq_summary: &[],
        diagnostics: std::slice::from_ref(diagnostic),
        ..writer_manifest(parts)
    }
}

type Tweak = Box<dyn Fn(&mut EventManifest<'static>, &Parts)>;

struct WriterCase {
    label: &'static str,
    failed: bool,
    tweak: Tweak,
    detail: &'static str,
}

fn writer_case(
    label: &'static str,
    detail: &'static str,
    tweak: impl Fn(&mut EventManifest<'static>, &Parts) + 'static,
) -> WriterCase {
    WriterCase {
        label,
        failed: false,
        tweak: Box::new(tweak),
        detail,
    }
}

fn failed_writer_case(
    label: &'static str,
    detail: &'static str,
    tweak: impl Fn(&mut EventManifest<'static>, &Parts) + 'static,
) -> WriterCase {
    WriterCase {
        failed: true,
        ..writer_case(label, detail, tweak)
    }
}

/// **Los hechos que el escritor comprueba antes de serializar.**
///
/// Estas guardas (`validate_writer_facts`) juzgan el `EventManifest` en
/// memoria, no un fichero publicado: un manifiesto publicado y manipulado no
/// pasa por ellas, pero sí un productor que arme los hechos mal, y la API
/// pública permite armarlos. Cada caso estropea **un** hecho.
#[test]
fn the_writer_rejects_each_incoherent_fact_with_its_own_reason() {
    let mut failures = Failures::default();
    let control = parts();
    for (name, manifest) in [
        ("completo", writer_manifest(&control)),
        ("fallido", failed_manifest(&control)),
    ] {
        if let Err(error) = event_manifest_document(&manifest) {
            failures.push(format!(
                "CASE-FAIL control/{name}: tenía que aceptarse: {error:?}"
            ));
        }
    }
    let cases = vec![
        writer_case(
            "modo pedido distinto del efectivo del evento",
            "el modo de errores no coincide con la configuración efectiva",
            |m, _| m.error_quality.requested_mode = SequencingErrorsMode::Disabled,
        ),
        writer_case(
            "modo efectivo distinto del del evento",
            "el modo de errores no coincide con la configuración efectiva",
            |m, _| m.error_quality.effective_mode = SequencingErrorsMode::Disabled,
        ),
        writer_case(
            "marca de calidad distinta de la del modo",
            "el modo de errores no coincide con la configuración efectiva",
            |m, _| {
                m.error_quality
                    .qualities_represent_realized_error_probability = false
            },
        ),
        writer_case(
            "desplazamiento de calidad de R1",
            FIXED_ADJUSTMENTS,
            |m, _| m.error_quality.r1.quality_shift = 1,
        ),
        writer_case("escala de error de R1", FIXED_ADJUSTMENTS, |m, _| {
            m.error_quality.r1.error_scale = 2.0;
        }),
        writer_case(
            "desplazamiento de calidad de R2",
            FIXED_ADJUSTMENTS,
            |m, _| m.error_quality.r2.quality_shift = -1,
        ),
        writer_case("escala de error de R2", FIXED_ADJUSTMENTS, |m, _| {
            m.error_quality.r2.error_scale = 0.5;
        }),
        writer_case(
            "tasa de quimera distinta de la del evento",
            "la tasa de quimera no coincide con la configuración efectiva",
            |m, _| m.chimera_rate = 0.5,
        ),
        writer_case(
            "parejas quiméricas con tasa cero",
            "la tasa de quimera no coincide con la configuración efectiva",
            |m, _| m.realized_chimera_pairs = 1,
        ),
        writer_case(
            "plan técnico de otro evento",
            "el plan técnico no coincide con el evento o referencia",
            |m, parts| {
                m.identities = leak(plan(
                    "reads-otro",
                    parts.reference_id,
                    parts.event.individual_ids(),
                ));
            },
        ),
        writer_case(
            "plan técnico de otra referencia",
            "el plan técnico no coincide con el evento o referencia",
            |m, parts| {
                m.identities = leak(plan(
                    "reads-main",
                    &ReferenceId::parse("reference-otra", "test").expect("referencia"),
                    parts.event.individual_ids(),
                ));
            },
        ),
        writer_case(
            "la referencia no enlaza un FASTA",
            "el enlace de referencia no apunta a un FASTA",
            |m, _| {
                m.reference.fasta = leak(record(
                    "reference-main-fasta",
                    "fastq",
                    "resources/references/reference-main.fa",
                    'a',
                ));
            },
        ),
        writer_case(
            "plan con una biblioteca que el evento no tiene",
            "individuos técnicos y bibliotecas no coinciden",
            |m, parts| {
                let mut individuals = parts.event.individual_ids().to_vec();
                individuals.push(IndividualId::parse("sample-B", "test").expect("individuo"));
                m.identities = leak(plan("reads-main", parts.reference_id, &individuals));
            },
        ),
        writer_case(
            "diagnósticos repetidos",
            "diagnósticos repetidos",
            |m, _| {
                m.diagnostics = leak([writer_diagnostic(), writer_diagnostic()]);
            },
        ),
        writer_case(
            "artefacto de un individuo no seleccionado",
            "un artefacto refiere un individuo no seleccionado",
            |m, parts| {
                let stranger = leak(IndividualId::parse("sample-Z", "test").expect("individuo"));
                m.artifacts = leak([
                    EventManifestArtifact {
                        individual_id: stranger,
                        ..writer_artifact(parts, 0)
                    },
                    writer_artifact(parts, 1),
                ]);
            },
        ),
        writer_case(
            "artefacto con la biblioteca de otro",
            "un artefacto refiere una biblioteca ajena",
            |m, parts| {
                let foreign = leak(LibraryId::parse("LIB_ajena", "test").expect("biblioteca"));
                m.artifacts = leak([
                    EventManifestArtifact {
                        library_id: foreign,
                        ..writer_artifact(parts, 0)
                    },
                    writer_artifact(parts, 1),
                ]);
            },
        ),
        writer_case(
            "mate repetido",
            "una biblioteca repite un mate FASTQ",
            |m, parts| m.artifacts = leak([writer_artifact(parts, 0), writer_artifact(parts, 0)]),
        ),
        writer_case(
            "artefacto completo sin identidad física",
            "un artefacto completo no tiene identidad física",
            |m, parts| {
                m.artifacts = leak([
                    EventManifestArtifact {
                        artifact: None,
                        ..writer_artifact(parts, 0)
                    },
                    writer_artifact(parts, 1),
                ]);
            },
        ),
        writer_case(
            "artefacto completo con diagnóstico",
            "artefacto completo incoherente",
            |m, parts| {
                m.artifacts = leak([
                    EventManifestArtifact {
                        diagnostic: Some(leak(writer_diagnostic())),
                        ..writer_artifact(parts, 0)
                    },
                    writer_artifact(parts, 1),
                ]);
            },
        ),
        writer_case(
            "identificador del artefacto distinto del de su identidad",
            "artefacto completo incoherente",
            |m, parts| {
                m.artifacts = leak([
                    EventManifestArtifact {
                        id: "otro",
                        ..writer_artifact(parts, 0)
                    },
                    writer_artifact(parts, 1),
                ]);
            },
        ),
        writer_case(
            "clase del artefacto distinta de la de su identidad",
            "artefacto completo incoherente",
            |m, parts| {
                m.artifacts = leak([
                    EventManifestArtifact {
                        kind: "otra",
                        ..writer_artifact(parts, 0)
                    },
                    writer_artifact(parts, 1),
                ]);
            },
        ),
        writer_case(
            "ruta del artefacto que no es la canónica",
            "artefacto completo incoherente",
            |m, parts| {
                m.artifacts = leak([
                    EventManifestArtifact {
                        artifact: Some(leak(record(
                            parts.ids[0],
                            "fastq",
                            "outputs/fastq/reads-main/otro_R1.fastq.gz",
                            'b',
                        ))),
                        ..writer_artifact(parts, 0)
                    },
                    writer_artifact(parts, 1),
                ]);
            },
        ),
        writer_case(
            "falta el resultado de un mate",
            "falta el resultado FASTQ de una biblioteca seleccionada",
            |m, parts| m.artifacts = leak([writer_artifact(parts, 0)]),
        ),
        writer_case(
            "cobertura repetida",
            "cobertura por biblioteca incompleta o repetida",
            |m, parts| m.coverage = leak([writer_coverage(parts), writer_coverage(parts)]),
        ),
        writer_case(
            "evento completo sin cobertura",
            "cobertura por biblioteca incompleta o repetida",
            |m, _| m.coverage = &[],
        ),
        writer_case(
            "resumen repetido",
            "resumen FASTQ por biblioteca incompleto o repetido",
            |m, parts| {
                m.final_fastq_summary = leak([writer_summary(parts), writer_summary(parts)]);
            },
        ),
        writer_case(
            "evento completo sin resumen",
            "resumen FASTQ por biblioteca incompleto o repetido",
            |m, _| m.final_fastq_summary = &[],
        ),
        // Los tres rechazos del escritor cuyo texto **repite** el del validador
        // que corre al final (`validate_event_manifest_bytes`): sin una segunda
        // incoherencia que el validador vea **antes**, quitar la guarda del
        // escritor daría el mismo motivo por otra vía y ninguna prueba lo
        // notaría. El resumen con las parejas propias y contaminantes que no
        // suman las parejas lo rechaza `validate_final_fastq_summary`, que corre
        // antes que el resto, y el escritor no lo mira.
        writer_case(
            "ajuste de calidad de R1 con un resumen incoherente detrás",
            FIXED_ADJUSTMENTS,
            |m, parts| {
                m.error_quality.r1.quality_shift = 1;
                m.final_fastq_summary = leak([incoherent_summary(parts)]);
            },
        ),
        writer_case(
            "escala de error de R2 con un resumen incoherente detrás",
            FIXED_ADJUSTMENTS,
            |m, parts| {
                m.error_quality.r2.error_scale = 0.5;
                m.final_fastq_summary = leak([incoherent_summary(parts)]);
            },
        ),
        writer_case(
            "falta un mate con un resumen incoherente detrás",
            "falta el resultado FASTQ de una biblioteca seleccionada",
            |m, parts| {
                m.artifacts = leak([writer_artifact(parts, 0)]);
                m.final_fastq_summary = leak([incoherent_summary(parts)]);
            },
        ),
        failed_writer_case(
            "evento no completo con cobertura y un resumen incoherente detrás",
            "un evento no completo no publica cobertura F4",
            |m, parts| {
                m.coverage = leak([writer_coverage(parts)]);
                m.final_fastq_summary = leak([incoherent_summary(parts)]);
            },
        ),
        failed_writer_case(
            "artefacto fallido sin diagnóstico",
            "artefacto no completo sin diagnóstico",
            |m, parts| {
                m.artifacts = leak([
                    EventManifestArtifact {
                        diagnostic: None,
                        ..failed_artifact(parts, 0, m.diagnostics)
                    },
                    failed_artifact(parts, 1, m.diagnostics),
                ]);
            },
        ),
        failed_writer_case(
            "artefacto incompleto sin diagnóstico",
            "artefacto no completo sin diagnóstico",
            |m, parts| {
                m.artifacts = leak([
                    EventManifestArtifact {
                        state: UnitState::Incomplete,
                        diagnostic: None,
                        ..failed_artifact(parts, 0, m.diagnostics)
                    },
                    failed_artifact(parts, 1, m.diagnostics),
                ]);
            },
        ),
        failed_writer_case(
            "artefacto fallido con identidad física",
            "artefacto no completo incoherente",
            |m, parts| {
                m.artifacts = leak([
                    EventManifestArtifact {
                        artifact: Some(parts.records[0]),
                        ..failed_artifact(parts, 0, m.diagnostics)
                    },
                    failed_artifact(parts, 1, m.diagnostics),
                ]);
            },
        ),
        failed_writer_case(
            "artefacto fallido con un diagnóstico que el manifiesto no lleva",
            "artefacto no completo incoherente",
            |m, parts| {
                let stranger = leak(other_diagnostic());
                m.artifacts = leak([
                    EventManifestArtifact {
                        diagnostic: Some(stranger),
                        ..failed_artifact(parts, 0, m.diagnostics)
                    },
                    failed_artifact(parts, 1, m.diagnostics),
                ]);
            },
        ),
        failed_writer_case(
            "evento no completo con cobertura",
            "un evento no completo no publica cobertura F4",
            |m, parts| m.coverage = leak([writer_coverage(parts)]),
        ),
    ];
    for case in cases {
        let parts = parts();
        let mut manifest = if case.failed {
            failed_manifest(&parts)
        } else {
            writer_manifest(&parts)
        };
        (case.tweak)(&mut manifest, &parts);
        match catch_unwind(AssertUnwindSafe(|| event_manifest_document(&manifest))) {
            Err(_) => failures.push(format!(
                "CASE-FAIL escritor/{}: pánico (esperado {:?})",
                case.label, case.detail
            )),
            Ok(Ok(_)) => failures.push(format!(
                "CASE-FAIL escritor/{}: se aceptó, y tenía que rechazarse (esperado {:?})",
                case.label, case.detail
            )),
            Ok(Err(error)) => {
                let detail = detail_of(&error);
                if detail != case.detail {
                    failures.push(format!(
                        "CASE-FAIL escritor/{}: motivo {detail:?}, esperado {:?}",
                        case.label, case.detail
                    ));
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} casos fallan; ver arriba",
        failures.len()
    );
}

fn failed_artifact(
    parts: &Parts,
    mate: usize,
    diagnostics: &'static [Diagnostic],
) -> EventManifestArtifact<'static> {
    EventManifestArtifact {
        state: UnitState::Failed,
        artifact: None,
        diagnostic: diagnostics.first(),
        ..writer_artifact(parts, mate)
    }
}

// ---------------------------------------------------------------------------
// Las guardas que ningún manifiesto publicado y manipulado alcanza
// ---------------------------------------------------------------------------

/// **Ejecuta el intento sobre las guardas inalcanzables y muestra qué otra
/// comprobación gana la carrera. NO DISCRIMINA esas guardas**: sus pruebas no
/// pueden fallar por ellas, porque ningún documento llega a ejercitarlas. Está
/// para que el porqué anotado junto a cada una en `src/event_manifest.rs` sea
/// una afirmación ejecutada y no solo escrita.
///
/// - «falta resumen FASTQ de una biblioteca completa»: con el estado completo,
///   coberturas y resúmenes pueden ser cada uno cualquier subconjunto de las
///   dos bibliotecas (16 combinaciones); el motivo nunca es ese.
/// - «un evento no completo no publica resúmenes FASTQ»: un resumen sin
///   cobertura lo rechaza antes `resumen FASTQ no tiene cobertura asociada`, y
///   con cobertura lo rechaza antes la cobertura del evento no completo.
/// - `cobertura propia/contaminante inválida` del cruce resumen-cobertura: una
///   cobertura con esos valores rota ya la rechazó `validate_coverage`.
/// - «descarte de demultiplexado con destino»: la clave de más la rechaza
///   antes `validate_closed`.
/// - «porcentaje de contaminación biológica ausente»: lo rechaza antes la
///   comprobación del porcentaje (`Option` desempaquetado dentro de un `matches!`
///   que ya exigió `Some`).
/// - `falta {campo}` de `required_value`: todo llamante lo precede un
///   `validate_closed` sobre el mismo objeto, y un objeto cerrado trae todas
///   sus claves; quitar cualquiera da `campos inválidos`.
/// - Ruta de salida repetida (`paths.insert`): un artefacto con la ruta de otro
///   incumple antes `record.path != paired_path`.
/// - Evento completo con una salida no completa y sin diagnósticos, y evento
///   no completo sin diagnósticos con una salida no completa: lo rechaza antes
///   `diagnóstico de artefacto ausente`.
#[test]
fn unreachable_guards_are_shadowed_by_other_checks() {
    let plain = plain();
    let detail = |document: &Value| -> String {
        match verdict(document) {
            Ok(()) => "ACEPTADO".to_owned(),
            Err(error) => detail_of(&error),
        }
    };
    let mut seen = Vec::new();

    // Coberturas y resúmenes, todos los subconjuntos, con el evento completo.
    for coverage_mask in 0..4_u8 {
        for summary_mask in 0..4_u8 {
            let mut document = plain.clone();
            neutral_demultiplexing(&mut document);
            for (list, mask) in [
                ("coverage", coverage_mask),
                ("final_fastq_summary", summary_mask),
            ] {
                let kept: Vec<Value> = (0..2)
                    .filter(|index| mask & (1 << index) != 0)
                    .map(|index| plain[list][index].clone())
                    .collect();
                set(&mut document, &format!("/{list}"), Value::Array(kept));
            }
            seen.push(detail(&document));
        }
    }
    assert!(
        !seen
            .iter()
            .any(|d| d == "falta resumen FASTQ de una biblioteca completa"),
        "esa guarda no es inalcanzable: {seen:?}"
    );

    // Evento no completo con resúmenes: coberturas y resúmenes cada uno
    // vacío o completo.
    let failed = failed_of(&plain);
    for coverage in [false, true] {
        for summary in [false, true] {
            let mut document = failed.clone();
            if coverage {
                set(&mut document, "/coverage", plain["coverage"].clone());
            }
            if summary {
                set(
                    &mut document,
                    "/final_fastq_summary",
                    plain["final_fastq_summary"].clone(),
                );
            }
            let got = detail(&document);
            assert_ne!(
                got, "un evento no completo no publica resúmenes FASTQ",
                "esa guarda no es inalcanzable (cobertura {coverage}, resumen {summary})"
            );
            if summary && !coverage {
                assert_eq!(got, "resumen FASTQ no tiene cobertura asociada");
            }
            if coverage {
                assert_eq!(got, "un evento no completo no publica cobertura F4");
            }
        }
    }

    // Cobertura propia o contaminante rota: la rechaza `validate_coverage`.
    for (pointer, expected) in [
        ("/coverage/0/final_x/own_x", "own_x de cobertura inválida"),
        (
            "/coverage/0/final_x/contaminant_x",
            "contaminant_x de cobertura inválida",
        ),
    ] {
        for value in [json!("x"), json!(null), json!(-1.0)] {
            let mut document = plain.clone();
            set(&mut document, pointer, value);
            assert_eq!(detail(&document), expected);
        }
    }

    // Descarte con destino: lo rechaza `validate_closed`.
    let discards = with_discards(&plain);
    let mut document = discards.clone();
    let library = text(
        &document,
        "/demultiplexing/discarded_by_sample/0/source_library_id",
    );
    set(
        &mut document,
        "/demultiplexing/discarded_by_sample/0/destination_library_id",
        json!(library),
    );
    assert_eq!(detail(&document), "campos inválidos en discarded_by_sample");

    // Porcentaje muestreado ausente o roto: lo rechaza el control del porcentaje.
    let mut document = plain.clone();
    set(&mut document, &contamination("sampled_percent"), json!("x"));
    assert_eq!(detail(&document), PERCENT);
    let mut document = plain.clone();
    del(&mut document, &contamination("sampled_percent"));
    assert_eq!(
        detail(&document),
        "campos inválidos en artifacts.biological_contamination[]"
    );

    // Una ruta de salida repetida: antes incumple la canónica.
    let mut document = plain.clone();
    let first = text(&document, "/outputs/0/artifact/path");
    set(&mut document, "/outputs/2/artifact/path", json!(first));
    assert_eq!(detail(&document), OUTPUT_INCOHERENT);

    // Salida no completa sin diagnósticos, con el evento completo o fallido.
    let mut document = plain.clone();
    set(&mut document, "/outputs/0/state", json!("failed"));
    set(&mut document, "/outputs/0/artifact", Value::Null);
    set(
        &mut document,
        "/outputs/0/diagnostic_code",
        json!("GPS-FMT-006"),
    );
    assert_eq!(detail(&document), "diagnóstico de artefacto ausente");
    set(&mut document, "/state", json!("failed"));
    assert_eq!(detail(&document), "diagnóstico de artefacto ausente");
}
