//! Cada rechazo del juez de los recursos publicados en lo que toca a la
//! configuración efectiva, con su motivo exacto.
//!
//! DT-001, tanda 7 (segunda parte de `src/results.rs`). El juez tenía 787
//! comprobaciones y 497 no las disparaba ninguna prueba; esta tanda cubre las
//! de la configuración efectiva: `validate_effective_config_v3_against_manifest_bytes`
//! (la efectiva frente al manifiesto), el enlace del evento técnico, las
//! poblaciones iniciales, los eventos directos, el evento técnico sin
//! amplificación, la referencia, las reducciones, los catálogos derivados, el
//! texto reversible y `validate_effective_targets_against_configuration` (la
//! efectiva atada al TOML copiado).
//!
//! # Cómo se parte de un paquete real
//!
//! Cada caso parte de un paquete **publicado de verdad** por el programa y
//! altera **un solo punto**; el motivo exacto del rechazo es la discriminación
//! (casi todas las guardas comparten la causa tipada `GPS-FMT-006`). Escenarios:
//!
//! - **`rich`**: referencia suministrada con dos secuencias (una con bases
//!   desconocidas), una reducción, dos cohortes (una con `founder_keys` y
//!   endogamia), tres eventos directos (SNP, inserción y deleción), un evento
//!   técnico con contaminación biológica, VCF de texto, ploidía global 3 y
//!   tolerancia propia (para que los valores esperados no coincidan con los
//!   de por defecto).
//! - **`simulated`**: referencia simulada con composición y bases
//!   desconocidas, sin poblaciones ni eventos, con dos secuencias.
//! - **`simulated` sin semilla declarada**: la semilla sale del sistema.
//! - **`inserting`**: referencia final por inserción de dos tractos.
//!
//! El juez de la efectiva frente al manifiesto
//! (`validate_effective_config_against_manifest_bytes`) recibe la efectiva y el
//! manifiesto alterados; el juez completo (`validate_published_resources`)
//! además lee el TOML copiado y los artefactos. Cada tabla lleva su control
//! positivo (el paquete intacto se acepta: sin él, un juez que rechazara todo
//! satisfaría la tabla entera), comprueba que **la alteración alteró**, y los
//! casos «aceptado» fijan el lado aceptado de los límites y de las listas
//! cerradas. Varios casos construyen un manifiesto **válido** distinto (sin
//! eventos, con dos eventos, fallido, sin reducción, sin VCF): las guardas
//! que el lector del manifiesto no sombrea solo se alcanzan así.
//!
//! # Lo que este fichero no es
//!
//! Las guardas que ningún paquete publicado y manipulado puede alcanzar están
//! anotadas **junto a la guarda** en `src/results.rs`, con el porqué real. La
//! prueba que solo «ejecuta el intento» sobre ellas lo dice en su docstring:
//! **no discrimina** esas guardas. El código muerto de la versión 1 del juez
//! (sin ningún llamante) se borró el 11/10/2026 y no se prueba.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::derive_execution_id;
use genome_population_simulator::integrity::configuration_digest;
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::results::{
    validate_effective_config_against_manifest_bytes, validate_published_resources,
};
use serde_json::{Value, json};

#[path = "support/detection_section.rs"]
mod detection_section;

static NONCE: AtomicU64 = AtomicU64::new(0);
const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

struct Package {
    root: PathBuf,
    destination: PathBuf,
}

impl Drop for Package {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

impl Package {
    fn json(&self, relative: &str) -> Value {
        serde_json::from_slice(&self.bytes(relative)).expect("JSON publicado")
    }
    fn bytes(&self, relative: &str) -> Vec<u8> {
        fs::read(self.destination.join(relative)).expect("artefacto publicado")
    }
    fn manifest_path(&self) -> PathBuf {
        self.destination.join("manifest.json")
    }
    fn effective_path(&self) -> PathBuf {
        self.destination.join("effective-config.json")
    }
}

fn publish(files: &[(&str, &str)]) -> Package {
    let root = std::env::temp_dir().join(format!(
        "gps-results-effective-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("directorio");
    for (name, content) in files {
        fs::write(
            root.join(name),
            detection_section::with_detection_section(content),
        )
        .expect("fichero de entrada");
    }
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&root.join("general.toml")).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Package { root, destination }
}

fn rich_reference() -> String {
    let mut chr1 = "ACGT".repeat(250);
    chr1.replace_range(10..20, "NNNNNNNNNN");
    format!(
        ">chr1 primary\n{chr1}\n>scaffold-1\n{}\n",
        "CGTA".repeat(100)
    )
}

fn rich_toml(technical: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"
ploidy = 3
normalization_tolerance = 1e-7

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
start = 101
end = 300
weight = 1.0

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[[population_snvs.generate.populations]]
population_id = "P2"
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
length = 1000

[[reference.sequences]]
id = "scaffold-1"
class = "scaffold"
length = 400

[[reduced_references]]
id = "ref-small"
source_reference_id = "ref-main"
region_size_bases = 100
target_bases = 700

[[initial_populations]]
population_id = "P1"
individual_count = 3
individual_ids = ["I001", "I002", "I003"]
ploidy = 2
genotype_model = "hwe"

[[initial_populations]]
population_id = "P2"
individual_count = 2
founder_keys = ["fk-a", "fk-b"]
ploidy = 2
genotype_model = "f_diploid"
inbreeding_coefficient = 0.5

[[direct_events]]
event_id = "ev-snp"
individual_id = "I001"
phase = 0
kind = "snp"
contig = "chr1"
position_1based = 500
ref = "T"
alt = "G"

[[direct_events]]
event_id = "ev-ins"
individual_id = "I002"
phase = 1
kind = "insertion"
contig = "chr1"
point = 600
bases = "AC"

[[direct_events]]
event_id = "ev-del"
individual_id = "I003"
phase = 0
kind = "deletion"
contig = "chr1"
start = 700
end = 710

[vcf_output]
mode = "text"
{technical}"#
    )
}

fn rich_technical(event_id: &str) -> String {
    format!(
        r#"
[technical_event]
event_id = "{event_id}"
individual_ids = ["I001", "I002"]
target_coverage_wgs = 1
contamination_min_percent = 10
contamination_max_percent = 20
contamination_recipient_ids = ["I002"]

[[technical_event.libraries]]
individual_id = "I001"
i7 = "ACGTACGT"
i5 = "TGCATGCA"

[[technical_event.libraries]]
individual_id = "I002"
i7 = "TTGGCCAA"
i5 = "GGAATTCC"
"#
    )
}

fn publish_rich() -> Package {
    publish_rich_event("reads-main")
}

fn publish_rich_event(event_id: &str) -> Package {
    publish(&[
        ("general.toml", &rich_toml(&rich_technical(event_id))),
        ("reference.fa", &rich_reference()),
    ])
}

// ---------------------------------------------------------------------------
// Casos
// ---------------------------------------------------------------------------

fn detail(result: Result<(), Diagnostic>) -> String {
    result
        .expect_err("el juez tenía que rechazar")
        .context()
        .get("detail")
        .cloned()
        .unwrap_or_default()
}

/// Qué juez recibe los documentos alterados.
#[derive(Clone, Copy)]
enum Judge {
    /// `validate_effective_config_against_manifest_bytes`: la efectiva frente
    /// al manifiesto, sin leer nada más del paquete.
    Pair,
    /// `validate_published_resources`: el juez completo, que además ata la
    /// efectiva al TOML copiado y lee los artefactos publicados.
    Full,
}

fn judge(
    kind: Judge,
    package: &Package,
    effective: &Value,
    manifest: &Value,
) -> Result<(), Diagnostic> {
    let effective_bytes = serde_json::to_vec_pretty(effective).expect("JSON");
    let manifest_bytes = serde_json::to_vec_pretty(manifest).expect("JSON");
    match kind {
        Judge::Pair => validate_effective_config_against_manifest_bytes(
            &package.effective_path(),
            &effective_bytes,
            &package.manifest_path(),
            &manifest_bytes,
        ),
        Judge::Full => validate_published_resources(
            &package.destination,
            &package.effective_path(),
            &effective_bytes,
            &package.manifest_path(),
            &manifest_bytes,
            None,
        ),
    }
}

type Mutate = Box<dyn Fn(&mut Value, &mut Value)>;

/// Un caso: una etiqueta, el texto **exacto** del rechazo y la alteración de
/// la configuración efectiva y del manifiesto general.
struct Case {
    label: String,
    expected: String,
    mutate: Mutate,
}

fn case(
    label: impl Into<String>,
    expected: impl Into<String>,
    mutate: impl Fn(&mut Value, &mut Value) + 'static,
) -> Case {
    Case {
        label: label.into(),
        expected: expected.into(),
        mutate: Box::new(mutate),
    }
}

/// Un caso que altera **solo** un punto de la efectiva.
fn eff(
    label: impl Into<String>,
    expected: impl Into<String>,
    pointer: impl Into<String>,
    value: Value,
) -> Case {
    let pointer = pointer.into();
    case(label, expected, move |effective, _| {
        *effective
            .pointer_mut(&pointer)
            .unwrap_or_else(|| panic!("el fixture no tiene {pointer}")) = value.clone();
    })
}

/// Un caso que quita una clave de la efectiva (el puntero termina en ella).
fn eff_remove(
    label: impl Into<String>,
    expected: impl Into<String>,
    pointer: impl Into<String>,
) -> Case {
    let pointer = pointer.into();
    case(label, expected, move |effective, _| {
        remove(effective, &pointer);
    })
}

fn remove(document: &mut Value, pointer: &str) {
    let (parent, key) = pointer.rsplit_once('/').expect("puntero");
    let parent = document
        .pointer_mut(parent)
        .unwrap_or_else(|| panic!("el fixture no tiene {parent}"));
    match parent {
        Value::Object(map) => {
            assert!(map.remove(key).is_some(), "no había {pointer}");
        }
        Value::Array(items) => {
            items.remove(key.parse::<usize>().expect("índice"));
        }
        _ => panic!("{pointer} no cuelga de un contenedor"),
    }
}

fn assert_cases(kind: Judge, package: &Package, cases: &[Case]) {
    let effective = package.json("effective-config.json");
    let manifest = package.json("manifest.json");
    judge(kind, package, &effective, &manifest)
        .expect("control positivo: el paquete intacto se acepta");
    for case in cases {
        let mut altered_effective = effective.clone();
        let mut altered_manifest = manifest.clone();
        (case.mutate)(&mut altered_effective, &mut altered_manifest);
        assert!(
            altered_effective != effective || altered_manifest != manifest,
            "«{}»: la alteración no alteró nada",
            case.label
        );
        assert_eq!(
            detail(judge(kind, package, &altered_effective, &altered_manifest)),
            case.expected,
            "«{}» se rechazó por otra razón",
            case.label
        );
    }
}

/// El lado aceptado de un límite o de una lista cerrada.
fn assert_accepted(kind: Judge, package: &Package, cases: &[Case]) {
    let effective = package.json("effective-config.json");
    let manifest = package.json("manifest.json");
    for case in cases {
        let mut altered_effective = effective.clone();
        let mut altered_manifest = manifest.clone();
        (case.mutate)(&mut altered_effective, &mut altered_manifest);
        assert!(
            altered_effective != effective || altered_manifest != manifest,
            "«{}»: la alteración no alteró nada",
            case.label
        );
        judge(kind, package, &altered_effective, &altered_manifest)
            .unwrap_or_else(|error| panic!("«{}» se tenía que aceptar: {error:?}", case.label));
    }
}

/// Un caso «aceptado» (sin texto de rechazo).
fn ok(label: impl Into<String>, pointer: impl Into<String>, value: Value) -> Case {
    eff(label, "", pointer, value)
}

// ---------------------------------------------------------------------------
// Poblaciones iniciales
// ---------------------------------------------------------------------------

#[test]
fn the_initial_populations_are_judged() {
    let package = publish_rich();
    let cases = [
        eff(
            "una cohorte que no es objeto",
            "initial_populations[] inválida",
            "/initial_populations/0",
            json!(5),
        ),
        eff(
            "dos cohortes con el mismo identificador",
            "population_id inicial repetido",
            "/initial_populations/1/population_id",
            json!("P1"),
        ),
        eff(
            "ploidía 11 en una cohorte",
            "ploidy inicial fuera de rango",
            "/initial_populations/0/ploidy",
            json!(11),
        ),
        eff(
            "un identificador de individuo que no es texto",
            "individual_id inicial inválido",
            "/initial_populations/0/individual_ids/0",
            json!(5),
        ),
        eff(
            "un identificador de individuo vacío",
            "individual_id inicial inválido",
            "/initial_populations/0/individual_ids/0",
            json!(""),
        ),
        eff(
            "unas founder_keys que no son una lista",
            "founder_keys inicial inválida",
            "/initial_populations/1/founder_keys",
            json!("fk-a"),
        ),
        eff(
            "menos founder_keys que individuos",
            "founder_keys inicial no coincide con individual_ids",
            "/initial_populations/1/founder_keys",
            json!(["fk-a"]),
        ),
        eff(
            "más founder_keys que individuos",
            "founder_keys inicial no coincide con individual_ids",
            "/initial_populations/1/founder_keys",
            json!(["fk-a", "fk-b", "fk-c"]),
        ),
        eff(
            "una founder_key que no es texto",
            "founder_key inicial inválida",
            "/initial_populations/1/founder_keys/0",
            json!(5),
        ),
        eff(
            "una founder_key con sintaxis inválida",
            "founder_key inicial inválida",
            "/initial_populations/1/founder_keys/0",
            json!("bad key!"),
        ),
        eff(
            "una founder_key repetida",
            "founder_key inicial repetida",
            "/initial_populations/1/founder_keys/1",
            json!("fk-a"),
        ),
        eff(
            "coeficiente de endogamia mayor que 1",
            "inbreeding_coefficient inicial inválido",
            "/initial_populations/1/genotype_model/inbreeding_coefficient",
            json!(1.5),
        ),
        eff(
            "coeficiente de endogamia negativo",
            "inbreeding_coefficient inicial inválido",
            "/initial_populations/1/genotype_model/inbreeding_coefficient",
            json!(-0.1),
        ),
        eff(
            "coeficiente de endogamia que no es número",
            "inbreeding_coefficient inicial inválido",
            "/initial_populations/1/genotype_model/inbreeding_coefficient",
            json!("0.5"),
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);
    assert_accepted(
        Judge::Pair,
        &package,
        &[
            ok(
                "ploidía 10 en una cohorte",
                "/initial_populations/0/ploidy",
                json!(10),
            ),
            ok(
                "coeficiente de endogamia 0",
                "/initial_populations/1/genotype_model/inbreeding_coefficient",
                json!(0.0),
            ),
            ok(
                "coeficiente de endogamia 1",
                "/initial_populations/1/genotype_model/inbreeding_coefficient",
                json!(1.0),
            ),
        ],
    );
}

// ---------------------------------------------------------------------------
// Eventos directos
// ---------------------------------------------------------------------------

#[test]
fn the_direct_events_are_judged() {
    let package = publish_rich();
    let cases = [
        eff(
            "un evento directo que no es objeto",
            "direct_events[] inválido",
            "/direct_events/0",
            json!(5),
        ),
        eff(
            "una clase de evento desconocida",
            "direct_events[].kind inválido",
            "/direct_events/0/kind",
            json!("inversion"),
        ),
        eff(
            "un event_id con sintaxis inválida",
            "event_id directo inválido",
            "/direct_events/0/event_id",
            json!("bad id!"),
        ),
        eff(
            "un event_id repetido",
            "event_id directo repetido",
            "/direct_events/1/event_id",
            json!("ev-snp"),
        ),
        eff(
            "un individual_id con sintaxis inválida",
            "individual_id directo inválido",
            "/direct_events/0/individual_id",
            json!("bad id!"),
        ),
        eff(
            "un contig con sintaxis inválida",
            "contig directo inválido",
            "/direct_events/0/contig",
            json!("bad id!"),
        ),
        eff(
            "fase 10",
            "fase directa fuera de rango",
            "/direct_events/0/phase",
            json!(10),
        ),
        eff(
            "un SNP con la base de referencia N",
            "alelo directo inválido",
            "/direct_events/0/ref",
            json!("N"),
        ),
        eff(
            "un SNP con la base de referencia en minúscula",
            "alelo directo inválido",
            "/direct_events/0/ref",
            json!("t"),
        ),
        eff(
            "un SNP con dos bases de referencia",
            "alelo directo inválido",
            "/direct_events/0/ref",
            json!("TG"),
        ),
        eff(
            "un SNP con la base alterna N",
            "alelo directo inválido",
            "/direct_events/0/alt",
            json!("N"),
        ),
        eff(
            "un SNP con dos bases alternas",
            "alelo directo inválido",
            "/direct_events/0/alt",
            json!("GG"),
        ),
        eff(
            "una inserción con una base que no es ACGT",
            "bases de inserción inválidas",
            "/direct_events/1/bases",
            json!("AN"),
        ),
        eff(
            "una inserción con bases en minúscula",
            "bases de inserción inválidas",
            "/direct_events/1/bases",
            json!("ac"),
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);
    // Ninguna base en minúscula, una por una: la lista cerrada es ACGT en
    // mayúscula, y admitir una sola minúscula más la rompería.
    let lowercase = ["a", "c", "g", "t"].map(|base| {
        [
            eff(
                format!("un SNP con la base de referencia {base} en minúscula"),
                "alelo directo inválido",
                "/direct_events/0/ref",
                json!(base),
            ),
            eff(
                format!("una inserción de la base {base} en minúscula"),
                "bases de inserción inválidas",
                "/direct_events/1/bases",
                json!(base),
            ),
        ]
    });
    for pair in &lowercase {
        assert_cases(Judge::Pair, &package, pair);
    }
    assert_accepted(
        Judge::Pair,
        &package,
        &[
            ok("fase 9", "/direct_events/0/phase", json!(9)),
            ok(
                "una inserción larga",
                "/direct_events/1/bases",
                json!("ACGTACGT"),
            ),
        ],
    );
}

// ---------------------------------------------------------------------------
// El evento técnico (sin amplificación)
// ---------------------------------------------------------------------------

/// Convierte el documento en uno de la versión 4: sin las tres secciones de
/// repeticiones, que solo existen en la 6.
fn as_version_4(effective: &mut Value) {
    let object = effective.as_object_mut().expect("objeto");
    object.insert("schema_version".to_owned(), json!("4"));
    for key in ["repeat_loci", "repeat_admission", "repeat_output_budget"] {
        assert!(object.remove(key).is_some(), "el fixture no tiene {key}");
    }
    // La versión 4 del `vcf_output` no lleva la versión de la normalización.
    assert!(
        effective["vcf_output"]
            .as_object_mut()
            .expect("vcf_output")
            .remove("normalization_version")
            .is_some()
    );
}

const NO_PCR_FREE: &str = "artefactos de duplicación del perfil PCR-free inválidos";
const GC_OR_ADAPTERS: &str = "artefactos GC o adaptadores del perfil inicial inválidos";

#[test]
fn the_technical_event_is_judged() {
    let package = publish_rich();
    let cases = [
        case(
            "tartamudeo en un documento de la versión 4",
            "la sección de tartamudeo no pertenece a esta versión de configuración efectiva",
            |effective, _| {
                as_version_4(effective);
                effective["technical_event"]["stutter"] = json!({});
            },
        ),
        eff(
            "un event_id técnico con sintaxis inválida",
            "event_id técnico inválido",
            "/technical_event/event_id",
            json!("bad id!"),
        ),
        eff(
            "un event_id técnico igual al de un evento directo",
            "event_id técnico colisiona con un evento directo",
            "/technical_event/event_id",
            json!("ev-snp"),
        ),
        eff(
            "un catálogo de identificadores que no es el del programa",
            "catálogo técnico F4 inválido",
            "/technical_event/id_catalog/version",
            json!("99"),
        ),
        eff(
            "un perfil que no es el del programa",
            "perfil técnico mínimo inválido",
            "/technical_event/profile/sha256",
            json!("0".repeat(64)),
        ),
        eff(
            "longitud de lectura 0",
            "perfil técnico mínimo inválido",
            "/technical_event/read_length",
            json!(0),
        ),
        eff(
            "longitud de lectura 151",
            "perfil técnico mínimo inválido",
            "/technical_event/read_length",
            json!(151),
        ),
        eff(
            "rango de contaminación con el mínimo sobre el máximo",
            "rango de contaminación biológica inválido",
            "/technical_event/biological_contamination/percent_range/minimum",
            json!(25.0),
        ),
        eff(
            "rango de contaminación con el mínimo negativo",
            "rango de contaminación biológica inválido",
            "/technical_event/biological_contamination/percent_range/minimum",
            json!(-1.0),
        ),
        eff(
            "rango de contaminación con el máximo sobre 100",
            "rango de contaminación biológica inválido",
            "/technical_event/biological_contamination/percent_range/maximum",
            json!(100.5),
        ),
        eff(
            "rango de contaminación con el máximo negativo",
            "rango de contaminación biológica inválido",
            "/technical_event/biological_contamination/percent_range/maximum",
            json!(-0.5),
        ),
        eff(
            "rango de contaminación con el mínimo sobre 100",
            "rango de contaminación biológica inválido",
            "/technical_event/biological_contamination/percent_range",
            json!({"minimum": 100.5, "maximum": 100.5}),
        ),
        eff(
            "rango de contaminación que no es número",
            "rango de contaminación biológica inválido",
            "/technical_event/biological_contamination/percent_range/maximum",
            json!("20"),
        ),
        eff(
            "calidad de ACGT distinta de 40",
            "calidad técnica mínima inválida",
            "/technical_event/quality/acgt_phred",
            json!(41),
        ),
        eff(
            "calidad de N distinta de 0",
            "calidad técnica mínima inválida",
            "/technical_event/quality/n_phred",
            json!(1),
        ),
        eff(
            "modo efectivo de errores distinto del pedido",
            "modo de errores de secuenciación efectivo inválido",
            "/technical_event/sequencing_errors/effective_mode",
            json!("disabled"),
        ),
        case(
            "un modo de errores fuera de la lista cerrada",
            "modo de errores de secuenciación efectivo inválido",
            |effective, _| {
                let errors = &mut effective["technical_event"]["sequencing_errors"];
                errors["requested_mode"] = json!("auto");
                errors["effective_mode"] = json!("auto");
                errors["qualities_represent_realized_error_probability"] = json!(false);
            },
        ),
        eff(
            "la marca de calidades realizadas contraria al modo",
            "modo de errores de secuenciación efectivo inválido",
            "/technical_event/sequencing_errors/qualities_represent_realized_error_probability",
            json!(false),
        ),
        eff(
            "tasa de quimeras negativa",
            "artefacto de quimera efectivo inválido",
            "/technical_event/artifacts/chimera/rate",
            json!(-0.1),
        ),
        eff(
            "tasa de quimeras sobre 1",
            "artefacto de quimera efectivo inválido",
            "/technical_event/artifacts/chimera/rate",
            json!(1.1),
        ),
        eff(
            "tasa de quimeras que no es número",
            "artefacto de quimera efectivo inválido",
            "/technical_event/artifacts/chimera/rate",
            json!("0"),
        ),
        eff(
            "tasa de index hopping negativa",
            "artefacto de index hopping efectivo inválido",
            "/technical_event/artifacts/index_hopping/rate",
            json!(-0.1),
        ),
        eff(
            "tasa de index hopping sobre 1",
            "artefacto de index hopping efectivo inválido",
            "/technical_event/artifacts/index_hopping/rate",
            json!(1.1),
        ),
        eff(
            "tasa de index hopping que no es número",
            "artefacto de index hopping efectivo inválido",
            "/technical_event/artifacts/index_hopping/rate",
            json!("0"),
        ),
        eff(
            "PCR no aplicable que no es falso",
            NO_PCR_FREE,
            "/technical_event/artifacts/pcr/applicable",
            json!("no"),
        ),
        eff(
            "duplicados ópticos realizados",
            NO_PCR_FREE,
            "/technical_event/artifacts/optical_duplicates/realized_pairs",
            json!(1),
        ),
        eff(
            "curva de GC sin amplificación",
            GC_OR_ADAPTERS,
            "/technical_event/artifacts/gc/mode",
            json!("curve"),
        ),
        eff(
            "un modo de GC desconocido",
            GC_OR_ADAPTERS,
            "/technical_event/artifacts/gc/mode",
            json!("flat"),
        ),
        eff(
            "adaptadores activados",
            GC_OR_ADAPTERS,
            "/technical_event/artifacts/adapters/enabled",
            json!(true),
        ),
        eff(
            "un evento sin individuos",
            "technical_event sin individuos",
            "/technical_event/individual_ids",
            json!([]),
        ),
        eff(
            "un individuo técnico que no es texto",
            "individual_id técnico inválido",
            "/technical_event/individual_ids/0",
            json!(5),
        ),
        eff(
            "un individuo técnico con sintaxis inválida",
            "individual_id técnico inválido",
            "/technical_event/individual_ids/0",
            json!("bad id!"),
        ),
        eff(
            "un individuo técnico repetido",
            "individual_id técnico repetido",
            "/technical_event/individual_ids/1",
            json!("I001"),
        ),
        eff(
            "un individuo técnico que no está en las poblaciones iniciales",
            "individual_id técnico fuera de las poblaciones iniciales",
            "/technical_event/individual_ids/1",
            json!("I009"),
        ),
        eff(
            "contaminación sin receptoras",
            "contaminación biológica sin receptoras",
            "/technical_event/biological_contamination/recipient_ids",
            json!([]),
        ),
        eff(
            "una receptora que no es texto",
            "receptora de contaminación biológica inválida",
            "/technical_event/biological_contamination/recipient_ids/0",
            json!(5),
        ),
        eff(
            "una receptora con sintaxis inválida",
            "receptora de contaminación biológica inválida",
            "/technical_event/biological_contamination/recipient_ids/0",
            json!("bad id!"),
        ),
        eff(
            "una receptora que no es individuo del evento",
            "receptoras de contaminación inválidas o repetidas",
            "/technical_event/biological_contamination/recipient_ids/0",
            json!("I003"),
        ),
        eff(
            "una receptora repetida",
            "receptoras de contaminación inválidas o repetidas",
            "/technical_event/biological_contamination/recipient_ids",
            json!(["I002", "I002"]),
        ),
        eff_remove(
            "menos bibliotecas que individuos",
            "libraries no coincide con individual_ids",
            "/technical_event/libraries/1",
        ),
        case(
            "una biblioteca de más",
            "libraries no coincide con individual_ids",
            |effective, _| {
                let copy = effective["technical_event"]["libraries"][0].clone();
                effective["technical_event"]["libraries"]
                    .as_array_mut()
                    .expect("bibliotecas")
                    .push(copy);
            },
        ),
        eff(
            "una biblioteca que no es objeto",
            "library inválida",
            "/technical_event/libraries/0",
            json!(5),
        ),
        eff(
            "una biblioteca con un individuo de sintaxis inválida",
            "library individual_id inválido",
            "/technical_event/libraries/0/individual_id",
            json!("bad id!"),
        ),
        eff(
            "una biblioteca de un individuo que no está en el evento",
            "library técnica no derivable",
            "/technical_event/libraries/0/individual_id",
            json!("I003"),
        ),
        eff(
            "dos bibliotecas del mismo individuo, solo con eso mal",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/1",
            json!({
                "individual_id": "I001",
                "library_id": "LIB_JCNGZ52VSMZ7MMFPVCDDYABCBQ",
                "i7": "TTGGCCAA",
                "i5": "GGAATTCC",
                "is_contamination_recipient": false,
            }),
        ),
        eff(
            "dos bibliotecas del mismo individuo",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/1/individual_id",
            json!("I001"),
        ),
        eff(
            "un library_id que no es el derivado",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/library_id",
            json!("LIB_OTRA"),
        ),
        eff(
            "la marca de receptora contraria a la lista",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/is_contamination_recipient",
            json!(true),
        ),
        eff(
            "un i7 de siete bases",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i7",
            json!("ACGTACG"),
        ),
        eff(
            "un i7 con una base que no es ACGT",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i7",
            json!("ACGTACGN"),
        ),
        eff(
            "un i7 en minúscula",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i7",
            json!("acgtacgt"),
        ),
        eff(
            "un i7 con la primera base a en minúscula",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i7",
            json!("aCGTACGT"),
        ),
        eff(
            "un i7 con la primera base c en minúscula",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i7",
            json!("cCGTACGT"),
        ),
        eff(
            "un i7 con la primera base g en minúscula",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i7",
            json!("gCGTACGT"),
        ),
        eff(
            "un i7 con la primera base t en minúscula",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i7",
            json!("tCGTACGT"),
        ),
        eff(
            "un i5 de nueve bases",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i5",
            json!("TGCATGCAT"),
        ),
        eff(
            "un i5 con una base que no es ACGT",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/0/i5",
            json!("TGCATGCX"),
        ),
        eff(
            "un i7 repetido",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/1/i7",
            json!("ACGTACGT"),
        ),
        eff(
            "un i5 repetido",
            "biblioteca UDI inválida o repetida",
            "/technical_event/libraries/1/i5",
            json!("TGCATGCA"),
        ),
        eff(
            "cobertura WGS negativa",
            "cobertura WGS técnica inválida",
            "/technical_event/target/coverage_wgs",
            json!(-1.0),
        ),
        eff(
            "cobertura WGS que no es número",
            "cobertura WGS técnica inválida",
            "/technical_event/target/coverage_wgs",
            json!("1"),
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);
    assert_accepted(
        Judge::Pair,
        &package,
        &[
            case(
                "un documento de la versión 4 sin tartamudeo",
                "",
                |effective, _| as_version_4(effective),
            ),
            ok(
                "longitud de lectura 1",
                "/technical_event/read_length",
                json!(1),
            ),
            ok(
                "rango de contaminación con mínimo igual al máximo",
                "/technical_event/biological_contamination/percent_range",
                json!({"minimum": 15.0, "maximum": 15.0}),
            ),
            ok(
                "rango de contaminación de 100 a 100",
                "/technical_event/biological_contamination/percent_range",
                json!({"minimum": 100.0, "maximum": 100.0}),
            ),
            ok(
                "rango de contaminación de 0 a 0",
                "/technical_event/biological_contamination/percent_range",
                json!({"minimum": 0.0, "maximum": 0.0}),
            ),
            ok(
                "rango de contaminación de 0 a 100",
                "/technical_event/biological_contamination/percent_range",
                json!({"minimum": 0.0, "maximum": 100.0}),
            ),
            case(
                "errores de secuenciación desactivados",
                "",
                |effective, _| {
                    let errors = &mut effective["technical_event"]["sequencing_errors"];
                    errors["requested_mode"] = json!("disabled");
                    errors["effective_mode"] = json!("disabled");
                    errors["qualities_represent_realized_error_probability"] = json!(false);
                },
            ),
            ok(
                "quimeras en 1",
                "/technical_event/artifacts/chimera/rate",
                json!(1.0),
            ),
            ok(
                "index hopping en 1",
                "/technical_event/artifacts/index_hopping/rate",
                json!(1.0),
            ),
            ok(
                "cobertura WGS 0",
                "/technical_event/target/coverage_wgs",
                json!(0.0),
            ),
        ],
    );
}

// ---------------------------------------------------------------------------
// La configuración efectiva frente al manifiesto, la referencia y las
// reducciones
// ---------------------------------------------------------------------------

/// Una referencia **simulada**, con una fracción de bases desconocidas, y sin
/// poblaciones ni eventos.
fn publish_simulated() -> Package {
    publish_simulated_with_seed(&format!("seed = \"{SEED}\"\n"))
}

/// Lo mismo **sin semilla declarada**: la raíz sale del sistema.
fn publish_simulated_with_a_system_seed() -> Package {
    publish_simulated_with_seed("")
}

fn publish_simulated_with_seed(seed_line: &str) -> Package {
    let toml = format!(
        r#"config_kind = "general"
schema_version = "4"
{seed_line}
[reference]
kind = "simulated"
id = "ref-sim"

[reference.composition]
A = 0.4
C = 0.1
G = 0.1
T = 0.4
N = 0.05

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1000

[[reference.sequences]]
id = "scaffold-1"
class = "scaffold"
length = 400
"#
    );
    publish(&[("general.toml", &toml)])
}

fn publish_inserting() -> Package {
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
    let reference = format!(">chr1\n{bases}\n>scaffold-1\n{}\n", noisy(300));
    let config = format!(
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
contig = "chr1"
start = 3000
end = 3500
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
id = "scaffold-1"
class = "scaffold"
length = 300

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

[[reference.repeat_loci]]
id = "STR-INS2"
sequence = "scaffold-1"
insert_after = 100
motif = "CAG"
reference_units = [10]

[[reference.repeat_loci.alleles]]
units = 10
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 11
frequency = 0.5

[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"

[vcf_output]
mode = "bgzf"
"#
    );
    publish(&[("general.toml", &config), ("reference.fa", &reference)])
}

// ---------------------------------------------------------------------------
// El enlace del evento técnico con el manifiesto
// ---------------------------------------------------------------------------

/// Un manifiesto **válido** sin ningún evento: la lista `events` falta, la de
/// `technical_events` está vacía y el catálogo no lleva ni el manifiesto de
/// evento ni los FASTQ.
fn without_events(manifest: &mut Value) {
    let object = manifest.as_object_mut().expect("objeto");
    object.remove("events");
    without_event_artifacts(manifest);
}

/// Lo mismo, con la lista `events` presente y vacía.
fn with_an_empty_events_list(manifest: &mut Value) {
    manifest["events"] = json!([]);
    without_event_artifacts(manifest);
}

fn without_event_artifacts(manifest: &mut Value) {
    let object = manifest.as_object_mut().expect("objeto");
    object["technical_events"] = json!([]);
    for field in ["minimum", "requested"] {
        manifest["output_categories"][field]
            .as_array_mut()
            .expect("categoría")
            .retain(|id| {
                let id = id.as_str().expect("id");
                !(id.starts_with("fastq-") || id.starts_with("event-manifest-"))
            });
    }
    manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
        .retain(|artifact| artifact["kind"] != "fastq" && artifact["kind"] != "event_manifest");
}

/// Un manifiesto **válido** con un segundo evento técnico, `reads-zzz`.
fn with_a_second_event(manifest: &mut Value) {
    let mut link = manifest["events"][0].clone();
    link["event_id"] = json!("reads-zzz");
    link["manifest_artifact_id"] = json!("event-manifest-reads-zzz");
    manifest["events"]
        .as_array_mut()
        .expect("events")
        .push(link);
    let mut record = manifest["artifacts"]
        .as_array()
        .expect("catálogo")
        .iter()
        .find(|artifact| artifact["id"] == "event-manifest-reads-main")
        .expect("manifiesto de evento")
        .clone();
    record["id"] = json!("event-manifest-reads-zzz");
    record["path"] = json!("events/reads-zzz/manifest.json");
    manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
        .push(record.clone());
    manifest["output_categories"]["minimum"]
        .as_array_mut()
        .expect("mínimo")
        .push(json!("event-manifest-reads-zzz"));
    let mut technical = manifest["technical_events"][0].clone();
    technical["event_id"] = json!("reads-zzz");
    technical["event_manifest"] = json!({
        "id": "event-manifest-reads-zzz",
        "kind": "event_manifest",
        "path": "events/reads-zzz/manifest.json",
        "checksum": record["checksum"].clone(),
        "size_bytes": record["size_bytes"].clone(),
    });
    manifest["technical_events"]
        .as_array_mut()
        .expect("technical_events")
        .push(technical);
}

/// Una ejecución **fallida** y válida: el manifiesto dice `failed` con un
/// diagnóstico del catálogo, y no publica ni la verdad ni el VCF.
fn as_a_failed_run(effective: &mut Value, manifest: &mut Value) {
    let line = include_str!("../catalogs/diagnostics-v3.tsv")
        .lines()
        .find(|line| !line.is_empty() && !line.starts_with('#'))
        .expect("un diagnóstico del catálogo");
    let fields = line.split('|').collect::<Vec<_>>();
    effective["vcf_output"] = Value::Null;
    manifest["state"] = json!("failed");
    manifest["diagnostics"] = json!([{"code": fields[0], "message": fields[2]}]);
    for field in ["minimum", "requested"] {
        manifest["output_categories"][field]
            .as_array_mut()
            .expect("categoría")
            .retain(|id| !id.as_str().expect("id").starts_with("truth-"));
    }
    manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
        .retain(|artifact| !artifact["id"].as_str().expect("id").starts_with("truth-"));
    manifest["vcf_outputs"] = json!([]);
}

/// Una ejecución fallida con el evento fallido.
fn with_a_failed_event(effective: &mut Value, manifest: &mut Value) {
    as_a_failed_run(effective, manifest);
    manifest["events"][0]["state"] = json!("failed");
}

/// Una ejecución fallida que **no llegó a publicar** la configuración
/// efectiva: sin enlace a ella ni artefacto, y sin eventos.
fn without_the_effective_link(effective: &mut Value, manifest: &mut Value) {
    without_events(manifest);
    as_a_failed_run(effective, manifest);
    manifest["effective_config_artifact_id"] = Value::Null;
    for field in ["minimum", "requested"] {
        manifest["output_categories"][field]
            .as_array_mut()
            .expect("categoría")
            .retain(|id| id != "effective-config");
    }
    manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
        .retain(|artifact| artifact["id"] != "effective-config");
}

#[test]
fn the_technical_event_link_is_judged() {
    let package = publish_rich();
    let other_event =
        publish_rich_event("reads-other").json("effective-config.json")["technical_event"].clone();
    let cases = [
        eff_remove(
            "un manifiesto con eventos y una efectiva sin technical_event",
            "un manifiesto con eventos exige technical_event en la configuración efectiva",
            "/technical_event",
        ),
        case(
            "un technical_event y un manifiesto sin eventos",
            "technical_event no tiene enlace en el manifiesto general",
            |_, manifest| without_events(manifest),
        ),
        case(
            "un technical_event y un manifiesto con la lista de eventos vacía",
            "F3 exige exactamente un enlace de evento técnico",
            |_, manifest| with_an_empty_events_list(manifest),
        ),
        case(
            "un technical_event y un manifiesto con dos eventos",
            "F3 exige exactamente un enlace de evento técnico",
            |_, manifest| with_a_second_event(manifest),
        ),
        case(
            "un technical_event con otro identificador que el enlace",
            "el enlace del evento no coincide con technical_event efectivo",
            move |effective, _| effective["technical_event"] = other_event.clone(),
        ),
        case(
            "un enlace de evento que no está complete",
            "el enlace del evento no coincide con technical_event efectivo",
            with_a_failed_event,
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);

    // El lado aceptado: sin evento técnico, el manifiesto no trae `events` o
    // lo trae vacío.
    let package = publish_simulated();
    assert_accepted(
        Judge::Pair,
        &package,
        &[case("events vacío", "", |_, manifest| {
            manifest["events"] = json!([])
        })],
    );
    judge(
        Judge::Pair,
        &package,
        &package.json("effective-config.json"),
        &package.json("manifest.json"),
    )
    .expect("sin events");
}

// ---------------------------------------------------------------------------
// La configuración efectiva frente al manifiesto: referencia y reducciones
// ---------------------------------------------------------------------------

/// Un manifiesto **válido** con una segunda reducción, `ref-small-2`, que
/// reutiliza los artefactos de la primera.
fn with_a_second_reduction(manifest: &mut Value) {
    let mut reduction = manifest["reference_reductions"][0].clone();
    reduction["id"] = json!("ref-small-2");
    manifest["reference_reductions"]
        .as_array_mut()
        .expect("reducciones")
        .push(reduction);
    manifest["references"]
        .as_array_mut()
        .expect("referencias")
        .push(json!({"id": "ref-small-2", "fasta_artifact_id": "reference-fasta-ref-small"}));
}

#[test]
fn the_effective_config_is_judged_against_the_manifest() {
    let package = publish_rich();
    let cases = [
        case(
            "un manifiesto fallido sin enlace a la configuración efectiva",
            "el manifiesto no enlaza la configuración efectiva",
            without_the_effective_link,
        ),
        eff(
            "ploidía global 11",
            "ploidy fuera de rango",
            "/ploidy",
            json!(11),
        ),
        eff(
            "un source_path que no es la procedencia de la copia",
            "source_path no coincide con la procedencia de la copia TOML",
            "/source_path",
            json!("/otro/general.toml"),
        ),
        eff(
            "otro identificador de referencia",
            "la referencia efectiva no coincide con el enlace del manifiesto",
            "/reference/id",
            json!("otra"),
        ),
        eff(
            "otro FASTA de referencia",
            "la referencia efectiva no coincide con el enlace del manifiesto",
            "/reference/fasta_artifact_id",
            json!("reference-fasta-ref-small"),
        ),
        eff_remove(
            "una reducción de menos en la efectiva",
            "las reducciones efectivas no coinciden con el manifiesto",
            "/reduced_references/0",
        ),
        case(
            "una reducción de más en la efectiva",
            "las reducciones efectivas no coinciden con el manifiesto",
            |effective, _| {
                let copy = effective["reduced_references"][0].clone();
                effective["reduced_references"]
                    .as_array_mut()
                    .expect("reducciones")
                    .push(copy);
            },
        ),
        case(
            "dos reducciones efectivas iguales frente a dos del manifiesto",
            "ID de reducción efectiva repetido",
            |effective, manifest| {
                with_a_second_reduction(manifest);
                let copy = effective["reduced_references"][0].clone();
                effective["reduced_references"]
                    .as_array_mut()
                    .expect("reducciones")
                    .push(copy);
            },
        ),
        case(
            "una reducción que reutiliza el identificador de la referencia",
            "una reducción reutiliza un ID de referencia",
            |effective, manifest| {
                for reduction in [
                    &mut effective["reduced_references"][0],
                    &mut manifest["reference_reductions"][0],
                ] {
                    reduction["id"] = json!("ref-main");
                    reduction["fasta_artifact_id"] = json!("reference-fasta");
                    reduction["source_reference_id"] = json!("ref-small");
                    reduction["source_fasta_artifact_id"] = json!("reference-fasta-ref-small");
                }
                // El FASTA de la referencia pasa a ser el de la reducción
                // (rol `requested`) y el de la antigua reducción, el del
                // mínimo.
                let minimum = manifest["output_categories"]["minimum"]
                    .as_array_mut()
                    .expect("mínimo");
                minimum.retain(|id| id != "reference-fasta");
                minimum.push(json!("reference-fasta-ref-small"));
                let requested = manifest["output_categories"]["requested"]
                    .as_array_mut()
                    .expect("solicitados");
                requested.retain(|id| id != "reference-fasta-ref-small");
                requested.push(json!("reference-fasta"));
            },
        ),
        eff(
            "una reducción efectiva que el manifiesto no tiene",
            "la reducción otra no existe en el manifiesto",
            "/reduced_references/0/id",
            json!("otra"),
        ),
        eff(
            "una reducción efectiva que no es objeto",
            "una reducción no es objeto",
            "/reduced_references/0",
            json!(5),
        ),
        eff(
            "un modo de objetivo desconocido",
            "target.mode inválido",
            "/reduced_references/0/target/mode",
            json!("percent"),
        ),
        eff(
            "una selección vacía",
            "selection_order inválido",
            "/reduced_references/0/selection_order",
            json!([]),
        ),
        eff(
            "una selección con un cero",
            "selection_order inválido",
            "/reduced_references/0/selection_order/2",
            json!(0),
        ),
        eff(
            "una selección con un valor que no es entero",
            "selection_order inválido",
            "/reduced_references/0/selection_order/2",
            json!("3"),
        ),
        eff(
            "una reducción efectiva que no está complete",
            "estado de reducción efectiva inválido",
            "/reduced_references/0/state",
            json!("incomplete"),
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);
    assert_accepted(
        Judge::Pair,
        &package,
        &[
            ok("ploidía global 10", "/ploidy", json!(10)),
            case(
                "una selección distinta, en las dos mitades",
                "",
                |effective, manifest| {
                    effective["reduced_references"][0]["selection_order"] = json!([1, 2, 3]);
                    manifest["reference_reductions"][0]["selection_order"] = json!([1, 2, 3]);
                },
            ),
        ],
    );
}

// ---------------------------------------------------------------------------
// La referencia efectiva
// ---------------------------------------------------------------------------

const BAD_ORIGIN: &str = "reference.composition no corresponde al origen declarado";
const BAD_PATH: &str = "reference.original_fasta_path no corresponde al origen declarado";
const BAD_FORMAT: &str = "reference.canonical_format no corresponde al origen declarado";

fn acgt(a: f64, c: f64, g: f64, t: f64) -> Value {
    json!({"A": a, "C": c, "G": g, "T": t})
}

/// La referencia **suministrada**: su composición es nula, su ruta es texto y
/// no tiene formato canónico.
#[test]
fn a_supplied_reference_is_judged() {
    let package = publish_rich();
    let cases = [
        eff(
            "un origen desconocido",
            "reference.source inválido",
            "/reference/source",
            json!("downloaded"),
        ),
        eff(
            "una composición con una clave desconocida",
            "campos desconocidos en reference.composition",
            "/reference/composition",
            json!({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25, "Z": 0.0}),
        ),
        eff(
            "una composición con la A que no es número",
            "reference.composition.A inválida",
            "/reference/composition",
            json!({"A": "0.25", "C": 0.25, "G": 0.25, "T": 0.5}),
        ),
        eff(
            "una composición sin la T",
            "reference.composition.T inválida",
            "/reference/composition",
            json!({"A": 0.25, "C": 0.25, "G": 0.5}),
        ),
        eff(
            "una composición con la C sobre 1",
            "reference.composition.C fuera de rango",
            "/reference/composition",
            acgt(0.0, 1.5, -0.25, -0.25),
        ),
        eff(
            "una composición con la G negativa",
            "reference.composition.G fuera de rango",
            "/reference/composition",
            acgt(0.5, 0.5, -0.25, 0.25),
        ),
        eff(
            "una composición con la N que no es número",
            "reference.composition.N inválida",
            "/reference/composition",
            json!({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25, "N": "0.1"}),
        ),
        eff(
            "una composición con la N sobre 1",
            "reference.composition.N fuera de rango",
            "/reference/composition",
            json!({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25, "N": 1.5}),
        ),
        eff(
            "una composición con la N negativa",
            "reference.composition.N fuera de rango",
            "/reference/composition",
            json!({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25, "N": -0.1}),
        ),
        eff(
            "una composición que no es objeto ni nula",
            "reference.composition inválida",
            "/reference/composition",
            json!("uniforme"),
        ),
        eff(
            "una composición válida en una referencia suministrada",
            BAD_ORIGIN,
            "/reference/composition",
            acgt(0.25, 0.25, 0.25, 0.25),
        ),
        eff(
            "sin ruta del FASTA original",
            BAD_PATH,
            "/reference/original_fasta_path",
            Value::Null,
        ),
        eff(
            "con la ruta del FASTA original vacía",
            BAD_PATH,
            "/reference/original_fasta_path",
            json!(""),
        ),
        eff(
            "con la ruta del FASTA original que no es texto",
            BAD_PATH,
            "/reference/original_fasta_path",
            json!(5),
        ),
        eff(
            "un formato canónico en una referencia suministrada",
            BAD_FORMAT,
            "/reference/canonical_format",
            json!({"version": "1", "line_base_count": 80, "line_ending": "LF", "base_case": "uppercase"}),
        ),
        eff(
            "una referencia sin secuencias",
            "reference.sequences está vacío",
            "/reference/sequences",
            json!([]),
        ),
        eff(
            "una secuencia que no es objeto",
            "una secuencia no es objeto",
            "/reference/sequences/1",
            json!(5),
        ),
        eff(
            "un identificador de secuencia repetido",
            "ID de secuencia repetido",
            "/reference/sequences/1/id",
            json!("chr1"),
        ),
        eff(
            "una clase de secuencia desconocida",
            "clase de secuencia inválida",
            "/reference/sequences/1/class",
            json!("plasmid"),
        ),
        case(
            "unas longitudes cuya suma desborda",
            "la suma de longitudes desborda u64",
            |effective, _| {
                effective["reference"]["sequences"][0]["length"] = json!(u64::MAX);
            },
        ),
        eff(
            "unas bases desconocidas que no son número",
            "unknown_base_count inválido",
            "/reference/sequences/0/unknown_base_count",
            json!("10"),
        ),
        eff(
            "unas bases desconocidas negativas",
            "unknown_base_count inválido",
            "/reference/sequences/0/unknown_base_count",
            json!(-1),
        ),
        eff(
            "más bases desconocidas que la longitud",
            "unknown_base_count excede length",
            "/reference/sequences/1/unknown_base_count",
            json!(401),
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);
    assert_accepted(
        Judge::Pair,
        &package,
        &[
            ok(
                "todas las bases desconocidas",
                "/reference/sequences/1/unknown_base_count",
                json!(400),
            ),
            ok(
                "una secuencia de clase chromosome",
                "/reference/sequences/1/class",
                json!("chromosome"),
            ),
        ],
    );
}

/// La referencia **simulada**: su composición es un objeto, su ruta es nula y
/// su formato canónico es el del programa.
#[test]
fn a_simulated_reference_is_judged() {
    let package = publish_simulated();
    let cases = [
        eff(
            "una composición con la A fuera de rango",
            "reference.composition.A fuera de rango",
            "/reference/composition",
            json!({"A": 1.5, "C": 0.1, "G": 0.1, "T": -0.7}),
        ),
        eff(
            "una composición sin la A",
            "reference.composition.A inválida",
            "/reference/composition",
            json!({"C": 0.1, "G": 0.1, "T": 0.8}),
        ),
        eff(
            "una composición que no suma uno",
            "reference.composition no suma 1 dentro de la tolerancia",
            "/reference/composition",
            json!({"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.3, "N": 0.05}),
        ),
        eff(
            "una composición nula en una referencia simulada",
            BAD_ORIGIN,
            "/reference/composition",
            Value::Null,
        ),
        eff(
            "una ruta de FASTA original en una referencia simulada",
            BAD_PATH,
            "/reference/original_fasta_path",
            json!("reference.fa"),
        ),
        eff(
            "sin formato canónico en una referencia simulada",
            BAD_FORMAT,
            "/reference/canonical_format",
            Value::Null,
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);
    assert_accepted(
        Judge::Pair,
        &package,
        &[
            ok(
                "una composición sin N",
                "/reference/composition",
                acgt(0.4, 0.1, 0.1, 0.4),
            ),
            ok(
                "una N de 0",
                "/reference/composition",
                json!({"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4, "N": 0.0}),
            ),
            ok(
                "una N de 1",
                "/reference/composition",
                json!({"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4, "N": 1.0}),
            ),
            ok(
                "una A de 1 y las demás de 0",
                "/reference/composition",
                acgt(1.0, 0.0, 0.0, 0.0),
            ),
        ],
    );
}

// ---------------------------------------------------------------------------
// Texto reversible y catálogos derivados
// ---------------------------------------------------------------------------

/// Los casos de un texto reversible (`{encoding, value}`) en `pointer`, que
/// la configuración efectiva llama `field`. Lo usan la cabecera y la
/// descripción de cada secuencia, de la referencia y de los catálogos
/// derivados: cada llamante tiene sus casos.
fn reversible_text_cases(field: &str, pointer: &str) -> (Vec<Case>, Vec<Case>) {
    let rejected = vec![
        eff(
            format!("{field}: no es objeto"),
            format!("{field} no es objeto"),
            pointer,
            json!("texto"),
        ),
        eff(
            format!("{field}: valor que no es texto"),
            format!("{field}.value inválido"),
            pointer,
            json!({"encoding": "utf-8", "value": 5}),
        ),
        eff(
            format!("{field}: valor vacío"),
            format!("{field}.value está vacío"),
            pointer,
            json!({"encoding": "utf-8", "value": ""}),
        ),
        eff(
            format!("{field}: codificación desconocida"),
            format!("{field}.encoding inválido"),
            pointer,
            json!({"encoding": "base64", "value": "AAAA"}),
        ),
        eff(
            format!("{field}: bytes-hex de longitud impar"),
            format!("{field}.encoding inválido"),
            pointer,
            json!({"encoding": "bytes-hex", "value": "abc"}),
        ),
        eff(
            format!("{field}: bytes-hex en mayúsculas"),
            format!("{field}.encoding inválido"),
            pointer,
            json!({"encoding": "bytes-hex", "value": "AB"}),
        ),
        eff(
            format!("{field}: bytes-hex con un carácter que no es hexadecimal"),
            format!("{field}.encoding inválido"),
            pointer,
            json!({"encoding": "bytes-hex", "value": "zz"}),
        ),
    ];
    let accepted = vec![
        ok(
            format!("{field}: bytes-hex válido"),
            pointer,
            json!({"encoding": "bytes-hex", "value": "00ff09af"}),
        ),
        ok(
            format!("{field}: utf-8"),
            pointer,
            json!({"encoding": "utf-8", "value": "texto"}),
        ),
    ];
    (rejected, accepted)
}

/// Los casos de un catálogo derivado de secuencias, en `prefix` (la lista
/// `sequences` del documento). El primer elemento es el cromosoma, el segundo
/// el andamio.
fn derived_catalog_cases(prefix: &str) -> (Vec<Case>, Vec<Case>) {
    let at = |tail: &str| format!("{prefix}{tail}");
    let rejected = vec![
        eff(
            "un catálogo derivado vacío",
            "catálogo derivado vacío",
            at(""),
            json!([]),
        ),
        eff(
            "una secuencia derivada que no es objeto",
            "una secuencia derivada no es objeto",
            at("/1"),
            json!(5),
        ),
        eff(
            "un identificador de secuencia derivada repetido",
            "ID de secuencia derivada repetido",
            at("/1/id"),
            json!("chr1"),
        ),
        eff(
            "una clase de secuencia derivada desconocida",
            "clase de secuencia derivada inválida",
            at("/1/class"),
            json!("plasmid"),
        ),
        eff(
            "unas longitudes derivadas cuya suma desborda",
            "la longitud derivada desborda u64",
            at("/0/length"),
            json!(u64::MAX),
        ),
        eff(
            "unas N derivadas que no son número",
            "n_count derivado inválido",
            at("/1/n_count"),
            json!("0"),
        ),
        eff(
            "unas N derivadas negativas",
            "n_count derivado inválido",
            at("/1/n_count"),
            json!(-1),
        ),
        eff(
            "más N derivadas que la longitud",
            "n_count derivado excede length",
            at("/1/n_count"),
            json!(100_000),
        ),
        eff(
            "un alfabeto derivado desconocido",
            "alphabet_validation derivada inválida",
            at("/1/alphabet_validation"),
            json!("ACGT"),
        ),
    ];
    (rejected, Vec::new())
}

fn derived_length_case(prefix: &str, length: u64) -> Case {
    ok(
        "todas las bases derivadas son N",
        format!("{prefix}/1/n_count"),
        json!(length),
    )
}

#[test]
fn the_reversible_text_and_the_derived_catalogs_of_a_reduction_are_judged() {
    let package = publish_rich();
    let mut cases = Vec::new();
    let mut accepted = Vec::new();
    for (field, pointer) in [
        (
            "reference.sequences[].header",
            "/reference/sequences/0/header",
        ),
        (
            "reference.sequences[].description",
            "/reference/sequences/0/description",
        ),
        (
            "reduced_references[].sequences[].header",
            "/reduced_references/0/sequences/0/header",
        ),
        (
            "reduced_references[].sequences[].description",
            "/reduced_references/0/sequences/0/description",
        ),
    ] {
        let (rejected, ok_cases) = reversible_text_cases(field, pointer);
        cases.extend(rejected);
        accepted.extend(ok_cases);
    }
    let (rejected, _) = derived_catalog_cases("/reduced_references/0/sequences");
    cases.extend(rejected);
    assert_cases(Judge::Pair, &package, &cases);
    accepted.push(derived_length_case("/reduced_references/0/sequences", 200));
    assert_accepted(Judge::Pair, &package, &accepted);
}

#[test]
fn the_derived_catalog_of_an_inserted_reference_is_judged() {
    let package = publish_inserting();
    let mut cases = Vec::new();
    let mut accepted = Vec::new();
    for (field, pointer) in [
        (
            "reduced_references[].sequences[].header",
            "/inserted_reference/sequences/0/header",
        ),
        (
            "reduced_references[].sequences[].description",
            "/inserted_reference/sequences/0/description",
        ),
    ] {
        let (rejected, ok_cases) = reversible_text_cases(field, pointer);
        cases.extend(rejected);
        accepted.extend(ok_cases);
    }
    let (rejected, _) = derived_catalog_cases("/inserted_reference/sequences");
    cases.extend(rejected);
    assert_cases(Judge::Pair, &package, &cases);
    accepted.push(derived_length_case("/inserted_reference/sequences", 330));
    assert_accepted(Judge::Pair, &package, &accepted);
}

// ---------------------------------------------------------------------------
// La configuración efectiva atada al TOML copiado
// ---------------------------------------------------------------------------

/// La huella del TOML copiado, la que entra en el identificador de ejecución.
fn toml_digest(package: &Package) -> String {
    let manifest = package.json("manifest.json");
    let relative = manifest["artifacts"]
        .as_array()
        .expect("catálogo")
        .iter()
        .find(|artifact| artifact["id"] == "config-general")
        .and_then(|artifact| artifact["path"].as_str())
        .expect("TOML copiado");
    configuration_digest(&package.bytes(relative)).expect("huella del TOML")
}

fn next_up(value: f64) -> f64 {
    assert!(value > 0.0 && value.is_finite());
    f64::from_bits(value.to_bits() + 1)
}

/// Lo que dice el juez completo sobre el TOML copiado cuando se sustituye por
/// `rewrite(texto)`, y se restaura después.
fn judge_with_toml(package: &Package, rewrite: impl Fn(&str) -> String) -> String {
    let manifest = package.json("manifest.json");
    let relative = manifest["artifacts"]
        .as_array()
        .expect("catálogo")
        .iter()
        .find(|artifact| artifact["id"] == "config-general")
        .and_then(|artifact| artifact["path"].as_str())
        .expect("TOML copiado")
        .to_owned();
    let path = package.destination.join(&relative);
    let original = fs::read_to_string(&path).expect("TOML copiado");
    let altered = rewrite(&original);
    assert_ne!(altered, original, "la alteración del TOML no alteró nada");
    fs::write(&path, detection_section::with_detection_section(&altered)).expect("TOML alterado");
    let result = judge(
        Judge::Full,
        package,
        &package.json("effective-config.json"),
        &manifest,
    );
    fs::write(&path, detection_section::with_detection_section(&original))
        .expect("TOML restaurado");
    detail(result)
}

/// El rechazo de una autoridad física (el TOML copiado), con el código y el
/// significado que el catálogo de diagnósticos le da.
fn physical(operation: &str, code: &str, meaning: &str) -> String {
    format!("autoridad física rechazada en {operation}: {code} {meaning}")
}

/// **El TOML copiado tiene que ser admisible, declarable y reconstruible.**
///
/// El juez completo relee el TOML copiado y rehace la configuración efectiva
/// desde él; si el TOML publicado no pasa la admisión, la lectura o la
/// reconstrucción, el paquete se rechaza con el diagnóstico del paso que
/// falla. Cada alteración se aplica sobre el fichero copiado y se restaura.
#[test]
fn the_copied_toml_must_be_admitted_declared_and_rebuilt() {
    let package = publish_rich();
    judge(
        Judge::Full,
        &package,
        &package.json("effective-config.json"),
        &package.json("manifest.json"),
    )
    .expect("control positivo: el paquete intacto se acepta");
    for (label, expected, from, to) in [
        (
            "una sintaxis TOML inválida",
            physical(
                "admit_configuration",
                "GPS-CONF-001",
                "El documento de configuración no tiene sintaxis TOML válida.",
            ),
            "schema_version = \"4\"",
            "schema_version = = \"4\"",
        ),
        (
            "un campo desconocido",
            physical(
                "parse_configuration",
                "GPS-CONF-006",
                "La configuración contiene un campo desconocido.",
            ),
            "schema_version = \"4\"",
            "schema_version = \"4\"\nbogus = 1",
        ),
        (
            "un objetivo de reducción que no cabe",
            physical(
                "rebuild_effective_config",
                "GPS-CONF-021",
                "El objetivo de la referencia reducida no es menor que el espacio elegible.",
            ),
            "target_bases = 700",
            "target_bases = 5000",
        ),
    ] {
        assert_eq!(
            judge_with_toml(&package, |text| text.replacen(from, to, 1)),
            expected,
            "«{label}» se rechazó por otra razón"
        );
    }
}

/// Quita del manifiesto y de la efectiva la reducción, con sus artefactos.
fn without_the_reduction(effective: &mut Value, manifest: &mut Value) {
    effective["reduced_references"] = json!([]);
    manifest["reference_reductions"] = json!([]);
    manifest["references"]
        .as_array_mut()
        .expect("referencias")
        .retain(|reference| reference["id"] != "ref-small");
    for field in ["minimum", "requested"] {
        manifest["output_categories"][field]
            .as_array_mut()
            .expect("categoría")
            .retain(|id| !id.as_str().expect("id").ends_with("-ref-small"));
    }
    manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
        .retain(|artifact| !artifact["id"].as_str().expect("id").ends_with("-ref-small"));
}

/// Quita del manifiesto y de la efectiva el VCF.
fn without_the_vcf(effective: &mut Value, manifest: &mut Value) {
    effective["vcf_output"] = Value::Null;
    manifest["vcf_outputs"] = json!([]);
    for field in ["minimum", "requested"] {
        manifest["output_categories"][field]
            .as_array_mut()
            .expect("categoría")
            .retain(|id| id.as_str().expect("id") != "truth-vcf");
    }
    manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
        .retain(|artifact| artifact["id"] != "truth-vcf");
}

/// El borde de `normalization.tolerance` en el juez de la pareja: el máximo
/// 1e-6 se acepta, y el primer valor por encima y el primer negativo por debajo
/// de 0 se rechazan con su motivo exacto.
#[test]
fn the_effective_normalization_tolerance_is_judged_at_its_edges() {
    let package = publish_rich();
    assert_cases(
        Judge::Pair,
        &package,
        &[
            eff(
                "tolerancia justo por encima de 1e-6",
                "normalization.tolerance inválida",
                "/normalization/tolerance",
                json!(next_up(1e-6)),
            ),
            eff(
                "tolerancia negativa mínima",
                "normalization.tolerance inválida",
                "/normalization/tolerance",
                json!(-f64::from_bits(1)),
            ),
        ],
    );
    assert_accepted(
        Judge::Pair,
        &package,
        &[
            ok("tolerancia 1e-6", "/normalization/tolerance", json!(1e-6)),
            ok("tolerancia 0", "/normalization/tolerance", json!(0.0)),
        ],
    );
}

#[test]
fn the_effective_config_must_be_the_one_the_toml_resolves() {
    let package = publish_rich();
    let tolerance = package.json("effective-config.json")["normalization"]["tolerance"]
        .as_f64()
        .expect("tolerancia");
    let cases = [
        eff(
            "un origen de la semilla que el TOML no declara",
            "la semilla efectiva no coincide con la declaración TOML",
            "/seed/origin",
            json!("system_csprng"),
        ),
        eff(
            "otra ploidía global",
            "la ploidía efectiva no coincide con la declaración TOML",
            "/ploidy",
            json!(4),
        ),
        eff(
            "una tolerancia un ulp mayor",
            "la tolerancia efectiva no coincide con la declaración TOML",
            "/normalization/tolerance",
            json!(next_up(tolerance)),
        ),
        case(
            "una reducción de más, coherente con el manifiesto",
            "las reducciones efectivas no coinciden con las declaradas en TOML",
            |effective, manifest| {
                with_a_second_reduction(manifest);
                let mut copy = effective["reduced_references"][0].clone();
                copy["id"] = json!("ref-small-2");
                effective["reduced_references"]
                    .as_array_mut()
                    .expect("reducciones")
                    .push(copy);
            },
        ),
        case(
            "sin la reducción que el TOML declara",
            "las reducciones efectivas no coinciden con las declaradas en TOML",
            without_the_reduction,
        ),
        case(
            "sin el VCF que el TOML declara",
            "vcf_output efectivo no coincide exactamente con el TOML resuelto",
            without_the_vcf,
        ),
    ];
    assert_cases(Judge::Full, &package, &cases);

    let package = publish_simulated();
    let cases = [
        case(
            "una secuencia de menos con la misma suma",
            "el número de secuencias efectivo no coincide con el TOML",
            |effective, _| {
                effective["reference"]["sequences"]
                    .as_array_mut()
                    .expect("secuencias")
                    .pop();
                effective["reference"]["sequences"][0]["length"] = json!(1400);
            },
        ),
        case(
            "una secuencia de más con la misma suma",
            "el número de secuencias efectivo no coincide con el TOML",
            |effective, _| {
                let mut extra = effective["reference"]["sequences"][1].clone();
                extra["id"] = json!("chr3");
                extra["length"] = json!(10);
                extra["unknown_base_count"] = json!(0);
                effective["reference"]["sequences"]
                    .as_array_mut()
                    .expect("secuencias")
                    .push(extra);
                effective["reference"]["sequences"][0]["length"] = json!(990);
            },
        ),
        case(
            "otro identificador de secuencia",
            "orden, ID, clase o longitud de secuencia no coincide con el TOML",
            |effective, _| effective["reference"]["sequences"][0]["id"] = json!("chrX"),
        ),
        case(
            "otra clase de secuencia",
            "orden, ID, clase o longitud de secuencia no coincide con el TOML",
            |effective, _| effective["reference"]["sequences"][0]["class"] = json!("scaffold"),
        ),
        case(
            "otro identificador de referencia, en las dos mitades",
            "la referencia efectiva no coincide con la declaración TOML",
            |effective, manifest| {
                effective["reference"]["id"] = json!("ref-otra");
                manifest["references"][0]["id"] = json!("ref-otra");
            },
        ),
        case(
            "otro tamaño de genoma con la suma de secuencias coherente",
            "la referencia efectiva no coincide con la declaración TOML",
            |effective, _| {
                effective["reference"]["sequences"][0]["length"] = json!(1001);
                effective["reference"]["total_genome_size"] = json!(1401);
            },
        ),
        case(
            "una referencia simulada presentada como suministrada",
            "la referencia efectiva no coincide con la declaración TOML",
            |effective, _| {
                effective["reference"]["source"] = json!("supplied");
                effective["reference"]["composition"] = Value::Null;
                effective["reference"]["original_fasta_path"] = json!("x.fa");
                effective["reference"]["canonical_format"] = Value::Null;
            },
        ),
        case(
            "otra longitud de secuencia con la misma suma",
            "orden, ID, clase o longitud de secuencia no coincide con el TOML",
            |effective, _| {
                effective["reference"]["sequences"][0]["length"] = json!(1001);
                effective["reference"]["sequences"][1]["length"] = json!(399);
            },
        ),
    ];
    assert_cases(Judge::Full, &package, &cases);

    // La raíz cambiada **con el identificador de ejecución rehecho** (es
    // público): solo la comparación de la raíz con la semilla que el TOML
    // declara la rechaza, y es la única defensa.
    for package in [publish_simulated(), publish_inserting()] {
        let cases = [case(
            "una raíz de semilla distinta con el execution_id rehecho",
            "la semilla efectiva no coincide con la declaración TOML",
            {
                let toml = toml_digest(&package);
                move |effective, manifest| {
                    let other = "1".repeat(64);
                    manifest["execution_id"] = json!(
                        derive_execution_id(
                            RootSeed::parse("semilla", &other).expect("semilla"),
                            &toml
                        )
                        .expect("execution_id")
                        .as_str()
                    );
                    effective["seed"]["root"] = json!(other);
                }
            },
        )];
        assert_cases(Judge::Full, &package, &cases);
    }

    // Sin semilla en el TOML, el origen es el sistema; el juez rehace el
    // origen desde el TOML y no admite que la efectiva diga que fue declarada.
    let package = publish_simulated_with_a_system_seed();
    assert_eq!(
        package.json("effective-config.json")["seed"]["origin"],
        "system_csprng"
    );
    let cases = [eff(
        "una semilla que se dice declarada y no lo es",
        "la semilla efectiva no coincide con la declaración TOML",
        "/seed/origin",
        json!("declared"),
    )];
    assert_cases(Judge::Full, &package, &cases);
}

// ---------------------------------------------------------------------------
// Lo que ningún paquete publicado y manipulado puede alcanzar
// ---------------------------------------------------------------------------

/// Un número JSON que no cabe en un `f64` finito (`1e999`).
fn overflowing_decimal() -> Value {
    serde_json::from_str("1e999").expect("un número JSON válido")
}

/// Lo que dice el juez de la efectiva frente al manifiesto con estos bytes.
fn judge_raw(package: &Package, effective: &[u8], manifest: &[u8]) -> String {
    detail(validate_effective_config_against_manifest_bytes(
        &package.effective_path(),
        effective,
        &package.manifest_path(),
        manifest,
    ))
}

/// **Intento ejecutado de alcanzar las guardas inalcanzables.** Cada caso
/// intenta llegar a una guarda anotada en `src/results.rs` como sombreada por
/// una comprobación anterior y comprueba que el rechazo viene de **esa otra**
/// comprobación, con su texto. **No discrimina las guardas anotadas**: si se
/// borrara una de ellas la prueba seguiría pasando, porque otra comprobación
/// ya rechaza antes; solo demuestra que la comprobación anterior existe y
/// rechaza lo que la guarda no vería.
///
/// Fija además la premisa de las mitades `is_finite` que se dan por
/// equivalentes: un decimal no finito no llega a la comprobación.
#[test]
fn unreachable_effective_guards_are_shadowed_by_other_checks() {
    let package = publish_rich();
    // Un texto o una descripción de secuencia ausentes: la lista cerrada de
    // campos de la secuencia ya exige que estén.
    let cases = [
        eff_remove(
            "sin la cabecera de una secuencia de la referencia",
            "campos ausentes o desconocidos en reference.sequences[]",
            "/reference/sequences/0/header",
        ),
        eff_remove(
            "sin la descripción de una secuencia de la referencia",
            "campos ausentes o desconocidos en reference.sequences[]",
            "/reference/sequences/0/description",
        ),
        eff_remove(
            "sin la cabecera de una secuencia derivada",
            "campos ausentes o desconocidos en reduced_references[].sequences[]",
            "/reduced_references/0/sequences/0/header",
        ),
        eff_remove(
            "sin la descripción de una secuencia derivada",
            "campos ausentes o desconocidos en reduced_references[].sequences[]",
            "/reduced_references/0/sequences/0/description",
        ),
        // Una raíz de semilla inválida: la valida `validate_effective_seed`
        // antes de que el TOML la reconstruya.
        eff(
            "una raíz de semilla en mayúsculas",
            "seed.root SHA-256 inválido",
            "/seed/root",
            json!("ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789"),
        ),
        // Una secuencia o una reducción que no son objeto: las rechaza la
        // comprobación de la referencia o de la reducción, antes de atarlas
        // al TOML.
        eff(
            "una secuencia de la referencia que no es objeto",
            "una secuencia no es objeto",
            "/reference/sequences/0",
            json!(5),
        ),
        eff(
            "una reducción que no es objeto",
            "una reducción no es objeto",
            "/reduced_references/0",
            json!(5),
        ),
        // La configuración del manifiesto: la lee y la valida el lector del
        // manifiesto antes que nadie.
        case(
            "un artifact_id de la copia TOML que no es texto",
            "falta artifact_id",
            |_, manifest| manifest["configuration_input"]["artifact_id"] = json!(5),
        ),
        case(
            "un artifact_id de la copia TOML que el catálogo no tiene",
            "el enlace nada no existe",
            |_, manifest| manifest["configuration_input"]["artifact_id"] = json!("nada"),
        ),
        case(
            "un original_path vacío",
            "original_path está vacío",
            |_, manifest| manifest["configuration_input"]["original_path"] = json!(""),
        ),
        case(
            "un manifiesto con otra versión",
            "versión de manifiesto no admitida: Some(\"99\"); se admiten exclusivamente \"4\" (manifest-general-v4.schema.json) y \"6\" (manifest-general-v6.schema.json)",
            |_, manifest| manifest["schema_version"] = json!("99"),
        ),
        case(
            "una reducción del manifiesto que no es objeto",
            "una reducción no es objeto",
            |_, manifest| manifest["reference_reductions"][0] = json!(5),
        ),
        case(
            "un enlace de evento que no es objeto",
            "un enlace de evento no es objeto",
            |_, manifest| manifest["events"][0] = json!(5),
        ),
        case(
            "un enlace de evento sin su manifiesto",
            "falta manifest_artifact_id",
            |_, manifest| manifest["events"][0]["manifest_artifact_id"] = Value::Null,
        ),
        case(
            "un enlace de evento a un artefacto que no es manifiesto de evento",
            "el enlace reference-fasta exige event_manifest, pero el catálogo declara reference_fasta",
            |_, manifest| manifest["events"][0]["manifest_artifact_id"] = json!("reference-fasta"),
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);

    let package = publish_inserting();
    let cases = [
        eff_remove(
            "sin la cabecera de una secuencia de la referencia final",
            "campos ausentes o desconocidos en reduced_references[].sequences[]",
            "/inserted_reference/sequences/0/header",
        ),
        eff_remove(
            "sin la descripción de una secuencia de la referencia final",
            "campos ausentes o desconocidos en reduced_references[].sequences[]",
            "/inserted_reference/sequences/0/description",
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);

    // Lo que no es JSON, o no es objeto: lo rechazan antes las lecturas del
    // manifiesto y de la efectiva.
    let intact_effective = package.bytes("effective-config.json");
    let intact_manifest = package.bytes("manifest.json");
    for raw in [&b"{"[..], &b"[]"[..]] {
        let manifest_says = judge_raw(&package, &intact_effective, raw);
        assert!(
            manifest_says.starts_with("JSON no legible: ")
                || manifest_says == "la raíz no es un objeto",
            "un manifiesto ilegible: {manifest_says}"
        );
        let effective_says = judge_raw(&package, raw, &intact_manifest);
        assert!(
            effective_says.starts_with("JSON no legible: ")
                || effective_says == "la raíz no es un objeto",
            "una efectiva ilegible: {effective_says}"
        );
    }

    // Las mitades `is_finite` de los decimales: con `serde_json` y
    // `arbitrary_precision`, `as_f64` devuelve `None` para un número que no
    // cabe en un `f64` finito, así que un decimal no finito nunca llega a la
    // comprobación `is_finite` (la rechaza la lectura del decimal, con el mismo
    // motivo). Se comprueba la premisa y cada sitio.
    assert!(
        serde_json::from_str::<Value>("1e999")
            .expect("un número JSON válido")
            .as_f64()
            .is_none()
    );
    let package = publish_rich();
    let cases = [
        eff(
            "una endogamia que no cabe en un f64",
            "inbreeding_coefficient inicial inválido",
            "/initial_populations/1/genotype_model/inbreeding_coefficient",
            overflowing_decimal(),
        ),
        eff(
            "un mínimo de contaminación que no cabe en un f64",
            "rango de contaminación biológica inválido",
            "/technical_event/biological_contamination/percent_range/minimum",
            overflowing_decimal(),
        ),
        eff(
            "un máximo de contaminación que no cabe en un f64",
            "rango de contaminación biológica inválido",
            "/technical_event/biological_contamination/percent_range/maximum",
            overflowing_decimal(),
        ),
        eff(
            "una tasa de quimeras que no cabe en un f64",
            "artefacto de quimera efectivo inválido",
            "/technical_event/artifacts/chimera/rate",
            overflowing_decimal(),
        ),
        eff(
            "una tasa de index hopping que no cabe en un f64",
            "artefacto de index hopping efectivo inválido",
            "/technical_event/artifacts/index_hopping/rate",
            overflowing_decimal(),
        ),
        eff(
            "una cobertura que no cabe en un f64",
            "cobertura WGS técnica inválida",
            "/technical_event/target/coverage_wgs",
            overflowing_decimal(),
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);
    let package = publish_simulated();
    let cases = [eff(
        "una N de composición que no cabe en un f64",
        "reference.composition.N inválida",
        "/reference/composition",
        json!({"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4, "N": overflowing_decimal()}),
    )];
    assert_cases(Judge::Pair, &package, &cases);

    // Una raíz de semilla distinta y válida, sin rehacer el identificador de
    // ejecución: lo rechaza este, que sale de la raíz y de los bytes del TOML.
    let package = publish_simulated();
    let cases = [eff(
        "una raíz de semilla válida y distinta",
        "execution_id no coincide con la semilla y los bytes TOML",
        "/seed/root",
        json!("1".repeat(64)),
    )];
    assert_cases(Judge::Full, &package, &cases);
}

// ---------------------------------------------------------------------------
// Guardas que ya disparaban otras pruebas, ahora con su motivo exacto
// ---------------------------------------------------------------------------

const AMPLIFIED_TOML: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/general.toml"
);
const AMPLIFIED_FASTA: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/reference.fa"
);
const AMPLIFIED_SNVS: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/population-snvs.tsv"
);

/// Las guardas de las mismas funciones que ya disparaban otras pruebas, **sin
/// exigir su motivo**: la medición de esta tanda mostró que ninguna prueba de
/// la batería cambiaba si se alteraba su texto. Se fijan aquí con el motivo
/// exacto, sobre el escenario `rich`.
#[test]
fn the_guards_other_tests_only_triggered_have_their_exact_motive() {
    let package = publish_rich();
    let cases = [
        eff(
            "otra clase de documento",
            "clase o versión incompatible",
            "/effective_config_kind",
            json!("otra"),
        ),
        eff(
            "unos loci declarados que el manifiesto no registra",
            "la configuración efectiva declara loci de repetición y el manifiesto no los registra",
            "/repeat_loci",
            json!([{}]),
        ),
        case(
            "una referencia de más en el manifiesto",
            "el conjunto de referencias del manifiesto no coincide exactamente con la configuración efectiva",
            |_, manifest| {
                manifest["references"]
                    .as_array_mut()
                    .expect("referencias")
                    .push(json!({"id": "ref-extra", "fasta_artifact_id": "reference-fasta"}));
            },
        ),
        case(
            "una reducción que el manifiesto cuenta con otras bases elegibles",
            "la reducción efectiva ref-small no coincide con el manifiesto",
            |_, manifest| manifest["reference_reductions"][0]["eligible_bases"] = json!(1399),
        ),
        eff(
            "un identificador de población con sintaxis inválida",
            "population_id inicial inválido",
            "/initial_populations/0/population_id",
            json!("bad id!"),
        ),
        eff(
            "un recuento de individuos que no coincide con la lista",
            "individual_count inicial no coincide con individual_ids",
            "/initial_populations/0/individual_count",
            json!(4),
        ),
        eff(
            "un identificador de individuo con sintaxis inválida",
            "individual_id inicial inválido",
            "/initial_populations/0/individual_ids/0",
            json!("bad id!"),
        ),
        eff(
            "un identificador de individuo repetido",
            "individual_id inicial repetido",
            "/initial_populations/0/individual_ids/1",
            json!("I001"),
        ),
        eff(
            "un modelo de genotipos desconocido",
            "genotype_model inicial inválido",
            "/initial_populations/0/genotype_model",
            json!({"kind": "otro"}),
        ),
        eff(
            "un modelo hwe con un campo de más",
            "genotype_model inicial inválido",
            "/initial_populations/0/genotype_model",
            json!({"kind": "hwe", "inbreeding_coefficient": 0.5}),
        ),
        eff(
            "un modelo f_diploid con ploidía distinta de 2",
            "genotype_model inicial inválido",
            "/initial_populations/1/ploidy",
            json!(3),
        ),
        case(
            "cohortes sin estadísticas iniciales",
            "faltan estadísticas iniciales para cohortes iniciales",
            |effective, manifest| {
                for document in [effective, manifest] {
                    assert!(
                        document["population_snvs"]
                            .as_object_mut()
                            .expect("matriz")
                            .remove("initial_population_statistics")
                            .is_some()
                    );
                }
            },
        ),
        eff(
            "un catálogo de azar que no es el del evento",
            "catálogo técnico F4 inválido",
            "/technical_event/rng_catalog",
            json!("f4-technical-v9"),
        ),
        eff(
            "un alfabeto de secuencia desconocido",
            "alphabet_validation inválida",
            "/reference/sequences/0/alphabet_validation",
            json!("ACGT"),
        ),
        eff(
            "un tamaño de genoma que no es la suma de las secuencias",
            "reference.total_genome_size no coincide con la suma de secuencias",
            "/reference/total_genome_size",
            json!(1500),
        ),
        eff(
            "un objetivo de reducción que no suma lo del catálogo derivado",
            "la suma del catálogo derivado no coincide con target.bases",
            "/reduced_references/0/target/bases",
            json!(701),
        ),
    ];
    assert_cases(Judge::Pair, &package, &cases);

    // Frente al TOML copiado.
    let cases = [
        case(
            "una reducción efectiva distinta del TOML, en las dos mitades",
            "la reducción efectiva ref-small no coincide completamente con el TOML",
            |effective, manifest| {
                effective["reduced_references"][0]["region_size_bases"] = json!(101);
                manifest["reference_reductions"][0]["region_size_bases"] = json!(101);
            },
        ),
        eff(
            "una procedencia del FASTA distinta de la del TOML",
            "la procedencia FASTA efectiva no coincide con la declaración TOML",
            "/reference/original_fasta_path",
            json!("/otro/reference.fa"),
        ),
        eff(
            "un evento directo distinto del TOML",
            "direct_events efectivo no coincide exactamente con el TOML resuelto",
            "/direct_events/0/alt",
            json!("C"),
        ),
        eff(
            "un evento técnico distinto del TOML",
            "technical_event efectivo no coincide exactamente con el TOML resuelto",
            "/technical_event/target/coverage_wgs",
            json!(2.0),
        ),
    ];
    assert_cases(Judge::Full, &package, &cases);

    let package = publish_simulated();
    let cases = [eff(
        "una composición distinta de la del TOML con la suma correcta",
        "la composición efectiva no coincide con la declaración TOML",
        "/reference/composition",
        json!({"A": 0.3, "C": 0.2, "G": 0.1, "T": 0.4, "N": 0.05}),
    )];
    assert_cases(Judge::Full, &package, &cases);

    // Las cohortes frente al TOML, sin evento técnico que las nombre.
    let package = publish_inserting();
    let cases = [eff(
        "un individuo inicial distinto del TOML",
        "las cohortes iniciales efectivas no coinciden exactamente con el TOML",
        "/initial_populations/0/individual_ids/0",
        json!("J1"),
    )];
    assert_cases(Judge::Full, &package, &cases);

    // La referencia final por inserción.
    let cases = [
        eff_remove(
            "una efectiva sin la referencia final que el manifiesto declara",
            "la configuración efectiva y el manifiesto no coinciden en si hay referencia final por inserción",
            "/inserted_reference",
        ),
        eff(
            "una referencia final que no es objeto",
            "inserted_reference no es objeto",
            "/inserted_reference",
            json!(5),
        ),
        eff(
            "otro identificador de la referencia final",
            "inserted_reference no coincide con la referencia declarada o con reference_insertion del manifiesto",
            "/inserted_reference/id",
            json!("ref-otra"),
        ),
        eff(
            "otro mapa de la referencia final",
            "inserted_reference no coincide con la referencia declarada o con reference_insertion del manifiesto",
            "/inserted_reference/map_artifact_id",
            json!("reference-fasta"),
        ),
        eff(
            "otros cambios de la referencia final",
            "inserted_reference no coincide con la referencia declarada o con reference_insertion del manifiesto",
            "/inserted_reference/insertions/0/locus_id",
            json!("STR-OTRO"),
        ),
        case(
            "una base insertada de más, con el catálogo y el manifiesto coherentes",
            "las regiones del generador de la matriz no son la traducción al genoma final de las del TOML",
            |effective, manifest| {
                effective["inserted_reference"]["insertions"][0]["inserted_bases"] = json!(31);
                manifest["reference_insertion"]["insertions"][0]["inserted_bases"] = json!(31);
                effective["inserted_reference"]["sequences"][0]["length"] = json!(4031);
            },
        ),
        case(
            "unos cambios de la referencia final distintos del TOML, en las dos mitades",
            "la referencia final, sus cambios o lo declarado de cada locus no coinciden con el TOML",
            |effective, manifest| {
                effective["inserted_reference"]["insertions"][0]["locus_id"] = json!("STR-OTRO");
                manifest["reference_insertion"]["insertions"][0]["locus_id"] = json!("STR-OTRO");
            },
        ),
    ];
    assert_cases(Judge::Full, &package, &cases);

    // El abanico del manifiesto con loci frente a una efectiva sin ellos.
    let package = publish(&[
        ("general.toml", AMPLIFIED_TOML),
        ("reference.fa", AMPLIFIED_FASTA),
        ("population-snvs.tsv", AMPLIFIED_SNVS),
    ]);
    let cases = [eff(
        "un manifiesto con loci y una efectiva sin ellos",
        "el manifiesto registra loci de repetición que la configuración efectiva no declara",
        "/repeat_loci",
        json!([]),
    )];
    assert_cases(Judge::Pair, &package, &cases);
}
