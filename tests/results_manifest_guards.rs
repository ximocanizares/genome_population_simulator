//! Cada rechazo del juez de los recursos publicados en lo que toca al
//! manifiesto general y a los recursos físicos, con su motivo exacto.
//!
//! DT-001, tanda 9 (tercera parte de `src/results.rs`). Cubre lo que las
//! tandas 6 y 7 dejaron: el lector del manifiesto general
//! (`validate_general_manifest_bytes`), el juez de los manifiestos de evento
//! (`validate_published_event_manifests`) y los jueces físicos de
//! `validate_published_resources` (evento técnico, contaminación biológica,
//! VCF, referencias y reducciones).
//!
//! # Cómo se parte de un paquete real
//!
//! Cada caso parte de un paquete **publicado de verdad** por el programa y
//! altera **un solo punto** (salvo donde el docstring del caso diga otra
//! cosa); el motivo exacto del rechazo es la discriminación. Cada tabla lleva
//! su control positivo (el paquete intacto se acepta: sin él, un juez que
//! rechazara todo satisfaría la tabla entera), comprueba que **la alteración
//! alteró**, y los casos «aceptado» fijan el lado aceptado de los límites y de
//! las listas cerradas.
//!
//! # Lo que este fichero no es
//!
//! Las guardas que ningún paquete publicado y manipulado puede alcanzar están
//! anotadas **junto a la guarda** en `src/results.rs`, con el porqué real. Las
//! pruebas que solo «ejecutan el intento» sobre ellas lo dicen en su
//! docstring: **no discriminan** esas guardas.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::{
    validate_general_manifest_bytes, validate_published_event_manifests,
    validate_published_resources,
};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

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
        "gps-results-manifest-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("directorio");
    // El paquete se crea **antes** de publicar: si algo falla y la prueba
    // entra en pánico, su `Drop` borra el directorio temporal.
    let destination = root.join("published");
    let package = Package { root, destination };
    for (name, content) in files {
        fs::write(
            package.root.join(name),
            detection_section::with_detection_section(content),
        )
        .expect("fichero de entrada");
    }
    execute_validated_config(
        &load_and_validate_config(&package.root.join("general.toml"))
            .expect("configuración válida"),
        &package.destination,
    )
    .expect("publicación");
    package
}

fn rich_reference() -> String {
    let mut chr1 = "ACGT".repeat(250);
    chr1.replace_range(10..20, "NNNNNNNNNN");
    format!(
        ">chr1 primary\n{chr1}\n>scaffold-1\n{}\n",
        "CGTA".repeat(100)
    )
}

fn rich_toml(vcf: &str, technical: &str) -> String {
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
mode = "{vcf}"
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

/// El paquete rico: referencia suministrada, una reducción, dos cohortes, tres
/// eventos directos, un evento técnico con contaminación biológica y VCF de
/// texto.
fn publish_rich() -> Package {
    publish(&[
        (
            "general.toml",
            &rich_toml("text", &rich_technical("reads-main")),
        ),
        ("reference.fa", &rich_reference()),
    ])
}

/// Lo mismo con el VCF **comprimido** (BGZF con su índice).
fn publish_rich_bgzf() -> Package {
    publish(&[
        (
            "general.toml",
            &rich_toml("bgzf", &rich_technical("reads-main")),
        ),
        ("reference.fa", &rich_reference()),
    ])
}

/// Una referencia **simulada**, sin poblaciones ni eventos.
fn publish_simulated() -> Package {
    let toml = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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

// ---------------------------------------------------------------------------
// Ayudantes
// ---------------------------------------------------------------------------

fn detail(result: Result<(), Diagnostic>) -> String {
    result
        .expect_err("el juez tenía que rechazar")
        .context()
        .get("detail")
        .cloned()
        .unwrap_or_default()
}

fn bytes_of(document: &Value) -> Vec<u8> {
    serde_json::to_vec_pretty(document).expect("JSON")
}

/// Lo que dice el lector del manifiesto general de este documento.
fn parse(manifest: &Value) -> Result<(), Diagnostic> {
    validate_general_manifest_bytes(Path::new("manifest.json"), &bytes_of(manifest))
}

/// Una alteración del manifiesto general.
type Mutate = Box<dyn Fn(&mut Value)>;

/// Un intento: etiqueta, alteración y el texto del rechazo que sí se produce.
type Attempt = (&'static str, Mutate, &'static str);

/// Un caso: una etiqueta, el texto **exacto** del rechazo (vacío si se
/// acepta) y la alteración del manifiesto.
struct Case {
    label: String,
    expected: String,
    mutate: Mutate,
}

fn case(
    label: impl Into<String>,
    expected: impl Into<String>,
    mutate: impl Fn(&mut Value) + 'static,
) -> Case {
    Case {
        label: label.into(),
        expected: expected.into(),
        mutate: Box::new(mutate),
    }
}

/// Un caso que cambia **un** punto del manifiesto.
fn set(
    label: impl Into<String>,
    expected: impl Into<String>,
    pointer: impl Into<String>,
    value: Value,
) -> Case {
    let pointer = pointer.into();
    case(label, expected, move |manifest| {
        *manifest
            .pointer_mut(&pointer)
            .unwrap_or_else(|| panic!("el fixture no tiene {pointer}")) = value.clone();
    })
}

/// Un caso que **quita** una clave (o un elemento) del manifiesto.
fn remove_case(
    label: impl Into<String>,
    expected: impl Into<String>,
    pointer: impl Into<String>,
) -> Case {
    let pointer = pointer.into();
    case(label, expected, move |manifest| remove(manifest, &pointer))
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

/// Aplica cada caso al manifiesto del paquete y exige el motivo exacto; los
/// de texto vacío se exigen **aceptados**.
fn assert_parse_cases(package: &Package, cases: &[Case]) {
    let manifest = package.json("manifest.json");
    parse(&manifest).expect("control positivo: el paquete intacto se acepta");
    for case in cases {
        let mut altered = manifest.clone();
        (case.mutate)(&mut altered);
        assert!(
            altered != manifest,
            "«{}»: la alteración no alteró nada",
            case.label
        );
        if case.expected.is_empty() {
            parse(&altered)
                .unwrap_or_else(|error| panic!("«{}» se tenía que aceptar: {error:?}", case.label));
        } else {
            assert_eq!(
                detail(parse(&altered)),
                case.expected,
                "«{}» se rechazó por otra razón",
                case.label
            );
        }
    }
}

/// Quita un artefacto del catálogo y de las cuatro categorías de salida.
fn drop_artifact(manifest: &mut Value, id: &str) {
    manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
        .retain(|artifact| artifact["id"] != id);
    for field in ["minimum", "defaulted", "requested", "disabled"] {
        manifest["output_categories"][field]
            .as_array_mut()
            .expect("categoría")
            .retain(|listed| listed != id);
    }
}

fn artifact_of(manifest: &Value, kind: &str) -> Value {
    manifest["artifacts"]
        .as_array()
        .expect("catálogo")
        .iter()
        .find(|artifact| artifact["kind"] == kind)
        .unwrap_or_else(|| panic!("el fixture no tiene un artefacto {kind}"))
        .clone()
}

fn catalog_path(manifest: &Value, id: &str) -> String {
    manifest["artifacts"]
        .as_array()
        .expect("catálogo")
        .iter()
        .find(|artifact| artifact["id"] == id)
        .and_then(|artifact| artifact["path"].as_str())
        .expect("ruta del artefacto")
        .to_owned()
}

/// Rehace, en todo el documento, la identidad de lo que está en `path`.
fn rehash_artifact(manifest: &mut Value, path: &str, bytes: &[u8]) {
    fn walk(value: &mut Value, path: &str, digest: &str, size: usize) {
        match value {
            Value::Object(map) => {
                if map.get("path").and_then(Value::as_str) == Some(path) {
                    if map.contains_key("checksum") {
                        map.insert(
                            "checksum".to_owned(),
                            json!({"algorithm": "sha256", "digest": digest}),
                        );
                    }
                    if map.contains_key("sha256") {
                        map.insert("sha256".to_owned(), json!(digest));
                    }
                    if map.contains_key("size_bytes") {
                        map.insert("size_bytes".to_owned(), json!(size));
                    }
                }
                for nested in map.values_mut() {
                    walk(nested, path, digest, size);
                }
            }
            Value::Array(items) => {
                for item in items {
                    walk(item, path, digest, size);
                }
            }
            _ => {}
        }
    }
    walk(
        manifest,
        path,
        &format!("{:x}", Sha256::digest(bytes)),
        bytes.len(),
    );
}

#[test]
fn the_fixtures_publish_and_are_accepted() {
    for package in [publish_rich(), publish_rich_bgzf(), publish_simulated()] {
        let manifest = package.json("manifest.json");
        let effective = package.json("effective-config.json");
        parse(&manifest).expect("el lector acepta el paquete");
        validate_published_resources(
            &package.destination,
            &package.effective_path(),
            &bytes_of(&effective),
            &package.manifest_path(),
            &bytes_of(&manifest),
            None,
        )
        .expect("el juez completo acepta el paquete");
        validate_published_event_manifests(
            &package.destination,
            &package.manifest_path(),
            &bytes_of(&manifest),
        )
        .expect("el juez de los manifiestos de evento acepta el paquete");
    }
}

// ---------------------------------------------------------------------------
// El lector del manifiesto general: la raíz, la entrada, las autorizaciones
// ---------------------------------------------------------------------------

const COMPLETE_NEEDS: &str =
    "una ejecución complete exige entrada, configuración efectiva, referencia y ningún diagnóstico";

/// Un diagnóstico del catálogo de la versión 3 de la categoría pedida.
fn catalog_diagnostic_of(category: &str) -> Value {
    let line = include_str!("../catalogs/diagnostics-v3.tsv")
        .lines()
        .find(|line| line.starts_with(&format!("GPS-{category}-")))
        .unwrap_or_else(|| panic!("el catálogo no tiene un diagnóstico {category}"));
    let fields = line.split('|').collect::<Vec<_>>();
    json!({"code": fields[0], "message": fields[2]})
}

/// **La raíz del manifiesto, su entrada y sus listas se juzgan.**
///
/// La clase del documento, la presencia y la forma de `configuration_input`,
/// el identificador de la configuración efectiva, el recuento de los
/// artefactos, el tipo de `artifacts`, de `output_categories` y de los
/// diagnósticos, la forma cerrada de las autorizaciones y la huella de los
/// metadatos de la matriz de SNV frente a su artefacto.
#[test]
fn the_manifest_root_is_judged() {
    let package = publish_rich();
    let mismatch = "population_snvs.metadata_sha256 no coincide con el artefacto";
    let bad_digest = "population_snvs.metadata_sha256 inválido";
    let cases = [
        set(
            "otra clase de documento",
            "clase o versión incompatible",
            "/manifest_kind",
            json!("gps_other_manifest"),
        ),
        remove_case(
            "sin configuration_input",
            "falta configuration_input",
            "/configuration_input",
        ),
        set(
            "configuration_input que no es objeto ni nulo",
            "configuration_input inválido",
            "/configuration_input",
            json!(5),
        ),
        set(
            "configuration_input de otra clase",
            "configuration_input.kind inválido",
            "/configuration_input/kind",
            json!("event"),
        ),
        set(
            "effective_config_artifact_id que no es texto",
            "effective_config_artifact_id inválido",
            "/effective_config_artifact_id",
            json!(5),
        ),
        set(
            "effective_config_artifact_id vacío",
            "effective_config_artifact_id inválido",
            "/effective_config_artifact_id",
            json!(""),
        ),
        set(
            "artifacts que no es una lista",
            "artifacts no es lista",
            "/artifacts",
            json!(5),
        ),
        set(
            "un artefacto que no es objeto",
            "un artefacto no es objeto",
            "/artifacts/0",
            json!(5),
        ),
        set(
            "una huella que no es sha256",
            "checksum no usa sha256",
            "/artifacts/0/checksum/algorithm",
            json!("md5"),
        ),
        set(
            "un tamaño que no es número",
            "size_bytes inválido",
            "/artifacts/0/size_bytes",
            json!("12"),
        ),
        set(
            "un tamaño negativo",
            "size_bytes inválido",
            "/artifacts/0/size_bytes",
            json!(-1),
        ),
        remove_case(
            "un artefacto sin tamaño",
            "size_bytes inválido",
            "/artifacts/0/size_bytes",
        ),
        case(
            "dos artefactos con el mismo identificador",
            "ID o ruta de artefacto repetidos",
            |manifest| {
                let first = manifest["artifacts"][0].clone();
                manifest["artifacts"][1]["id"] = first["id"].clone();
            },
        ),
        case(
            "dos artefactos con la misma ruta",
            "ID o ruta de artefacto repetidos",
            |manifest| {
                let first = manifest["artifacts"][0].clone();
                manifest["artifacts"][1]["path"] = first["path"].clone();
            },
        ),
        set(
            "output_categories que no es objeto",
            "output_categories no es objeto",
            "/output_categories",
            json!(5),
        ),
        set(
            "diagnostics con un valor que no es objeto",
            "un diagnóstico no es objeto",
            "/diagnostics",
            json!([5]),
        ),
        remove_case(
            "sin authorizations",
            "falta authorizations",
            "/authorizations",
        ),
        set(
            "authorizations que no es objeto",
            "authorizations inválido",
            "/authorizations",
            json!(5),
        ),
        case(
            "authorizations con una clave de más",
            "authorizations inválido",
            |manifest| {
                manifest["authorizations"]["extra"] = json!(true);
            },
        ),
        remove_case(
            "authorizations sin allow_empty_snv_matrix",
            "authorizations inválido",
            "/authorizations/allow_empty_snv_matrix",
        ),
        set(
            "allow_empty_snv_matrix que no es booleano",
            "authorizations inválido",
            "/authorizations/allow_empty_snv_matrix",
            json!("false"),
        ),
        set(
            "allow_incompatible_genotypes que no es booleano",
            "authorizations inválido",
            "/authorizations/allow_incompatible_genotypes",
            json!(0),
        ),
        // La huella de los metadatos: forma (64 hexadecimales en minúscula) y
        // coincidencia con el artefacto.
        set(
            "huella de los metadatos con una mayúscula",
            bad_digest,
            "/population_snvs/metadata_sha256",
            json!("A".repeat(64)),
        ),
        set(
            "huella de los metadatos con un carácter que no es hexadecimal",
            bad_digest,
            "/population_snvs/metadata_sha256",
            json!("g".repeat(64)),
        ),
        set(
            "huella de los metadatos de 63 caracteres",
            bad_digest,
            "/population_snvs/metadata_sha256",
            json!("a".repeat(63)),
        ),
        set(
            "huella de los metadatos de 65 caracteres",
            bad_digest,
            "/population_snvs/metadata_sha256",
            json!("a".repeat(65)),
        ),
        set(
            "huella de los metadatos con forma válida pero distinta (letras)",
            mismatch,
            "/population_snvs/metadata_sha256",
            json!("f".repeat(64)),
        ),
        set(
            "huella de los metadatos con forma válida pero distinta (dígitos)",
            mismatch,
            "/population_snvs/metadata_sha256",
            json!("0".repeat(64)),
        ),
        // La forma de una huella de artefacto: 64 hexadecimales en minúscula;
        // los dos extremos de la letra (`f` pasa, `g` no) y la longitud por los
        // dos lados.
        set(
            "una huella con un carácter después de la f",
            "digest SHA-256 inválido",
            "/artifacts/0/checksum/digest",
            json!("g".repeat(64)),
        ),
        set(
            "una huella de 63 caracteres",
            "digest SHA-256 inválido",
            "/artifacts/0/checksum/digest",
            json!("a".repeat(63)),
        ),
        set(
            "una huella de 65 caracteres",
            "digest SHA-256 inválido",
            "/artifacts/0/checksum/digest",
            json!("a".repeat(65)),
        ),
        set(
            "una huella con una mayúscula",
            "digest SHA-256 inválido",
            "/artifacts/0/checksum/digest",
            json!("A".repeat(64)),
        ),
        set(
            "una huella con la letra más baja y la más alta de la lista",
            "",
            "/artifacts/0/checksum/digest",
            json!(format!("{}{}", "a".repeat(32), "f".repeat(32))),
        ),
        set(
            "una huella con un dígito",
            "",
            "/artifacts/0/checksum/digest",
            json!("0".repeat(64)),
        ),
        set(
            "una semilla de reducción de 63 caracteres",
            "root_seed SHA-256 inválido",
            "/reference_reductions/0/root_seed",
            json!("a".repeat(63)),
        ),
        set(
            "una subsemilla de selección con una letra fuera de la lista",
            "selection_subseed SHA-256 inválido",
            "/reference_reductions/0/selection_subseed",
            json!("g".repeat(64)),
        ),
        set(
            "population_snvs que no es objeto",
            "population_snvs inválido",
            "/population_snvs",
            json!(5),
        ),
        set(
            "un modo de matriz desconocido",
            "population_snvs.mode inválido",
            "/population_snvs/mode",
            json!("imported"),
        ),
        set(
            "metadatos de la matriz generada sin identificador de artefacto",
            "population_snvs.metadata_artifact_id inválido",
            "/population_snvs/metadata_artifact_id",
            Value::Null,
        ),
    ];
    assert_parse_cases(&package, &cases);
}

/// **Un código de diagnóstico tiene la forma `GPS-<categoría>-<NNN>`.**
///
/// Cada alteración rompe una sola pieza del código: el prefijo, la categoría,
/// la longitud o el carácter del número, o sobra o falta un tramo. El lado
/// aceptado fija la lista cerrada de categorías: un diagnóstico del catálogo
/// de cada una de las ocho pasa el código (y se rechaza después por otra cosa:
/// una ejecución `complete` no puede llevar diagnósticos).
#[test]
fn a_diagnostic_code_has_the_shape_gps_category_number() {
    let package = publish_rich();
    let invalid = "código de diagnóstico inválido";
    let with_code = |code: &'static str| {
        case(format!("código {code}"), invalid, move |manifest| {
            manifest["diagnostics"] = json!([{"code": code, "message": "x"}]);
        })
    };
    let mut cases = vec![
        with_code("XYZ-CONF-001"),
        with_code("GPS-FOO-001"),
        with_code("GPS-CONF-01"),
        with_code("GPS-CONF-0001"),
        with_code("GPS-CONF-0a1"),
        with_code("GPS-CONF-001-X"),
        with_code("GPS-CONF"),
        with_code("GPS"),
        with_code("GPS-001"),
    ];
    for category in [
        "CONF", "INPUT", "FMT", "DOMAIN", "RES", "TOOL", "IO", "INTERNAL",
    ] {
        cases.push(case(
            format!("un diagnóstico {category} del catálogo se acepta como código"),
            COMPLETE_NEEDS,
            move |manifest| {
                manifest["diagnostics"] = json!([catalog_diagnostic_of(category)]);
            },
        ));
    }
    assert_parse_cases(&package, &cases);
}

/// **Una ejecución `complete` exige entrada, configuración efectiva,
/// referencia y ningún diagnóstico; una `failed`, al menos uno.**
///
/// Cada término de la condición se rompe por separado: la entrada nula, la
/// configuración efectiva nula (con el evento técnico que la embebe también
/// nulo), las referencias vacías, una reducción que no está completa
/// (`failed` e `incomplete`), un diagnóstico y un evento fallido. La entrada
/// y la configuración efectiva nulas se acompañan de quitar el artefacto del
/// catálogo y de las categorías, porque si no las rechazaría antes el
/// catálogo; es la ejecución que no llegó a publicarlos.
#[test]
fn a_complete_run_needs_its_input_its_effective_config_its_reference_and_no_diagnostics() {
    let package = publish_rich();
    let cases = [
        case("sin entrada", COMPLETE_NEEDS, |manifest| {
            manifest["configuration_input"] = Value::Null;
            drop_artifact(manifest, "config-general");
        }),
        case("sin configuración efectiva", COMPLETE_NEEDS, |manifest| {
            manifest["effective_config_artifact_id"] = Value::Null;
            drop_artifact(manifest, "effective-config");
            manifest["technical_events"][0]["effective_config"] = Value::Null;
        }),
        set(
            "una reducción fallida",
            COMPLETE_NEEDS,
            "/reference_reductions/0/state",
            json!("failed"),
        ),
        set(
            "una reducción incompleta",
            COMPLETE_NEEDS,
            "/reference_reductions/0/state",
            json!("incomplete"),
        ),
        case("un diagnóstico", COMPLETE_NEEDS, |manifest| {
            manifest["diagnostics"] = json!([catalog_diagnostic_of("IO")]);
        }),
        set(
            "un evento fallido",
            COMPLETE_NEEDS,
            "/events/0/state",
            json!("failed"),
        ),
    ];
    assert_parse_cases(&package, &cases);

    // Sin referencias (en un paquete que no tiene nada más que la referencia).
    let simulated = publish_simulated();
    assert_parse_cases(
        &simulated,
        &[case("sin referencias", COMPLETE_NEEDS, |manifest| {
            manifest["references"] = json!([]);
            drop_artifact(manifest, "reference-fasta");
        })],
    );

    // Una ejecución fallida sin ningún diagnóstico se rechaza; con uno, se
    // acepta. No publica la verdad: se quita del catálogo y de las categorías.
    fn failed(manifest: &mut Value) {
        manifest["state"] = json!("failed");
        drop_artifact(manifest, "truth-jsonl");
    }
    fn failed_with_a_diagnostic(manifest: &mut Value) {
        failed(manifest);
        manifest["diagnostics"] = json!([catalog_diagnostic_of("IO")]);
    }
    assert_parse_cases(
        &simulated,
        &[
            case(
                "fallida sin diagnósticos",
                "una ejecución failed exige al menos un diagnóstico",
                failed,
            ),
            case("fallida con un diagnóstico", "", failed_with_a_diagnostic),
        ],
    );
}

// ---------------------------------------------------------------------------
// Las categorías de salida
// ---------------------------------------------------------------------------

fn push(manifest: &mut Value, pointer: &str, value: Value) {
    manifest
        .pointer_mut(pointer)
        .unwrap_or_else(|| panic!("el fixture no tiene {pointer}"))
        .as_array_mut()
        .expect("lista")
        .push(value);
}

/// Mueve un identificador de una categoría de salida a otra.
fn move_category(manifest: &mut Value, id: &str, from: &str, to: &str) {
    manifest["output_categories"][from]
        .as_array_mut()
        .expect("categoría")
        .retain(|listed| listed != id);
    push(manifest, &format!("/output_categories/{to}"), json!(id));
}

/// **Las categorías de salida se juzgan.**
///
/// Cada identificador listado tiene que ser un texto no vacío, estar en el
/// catálogo y no repetirse ni dentro de una categoría ni entre categorías;
/// `defaulted` y `disabled` están vacías en F1; `minimum` y `requested` son
/// exactamente los roles normativos; y el catálogo no lleva nada sin rol.
///
/// El caso de `minimum` y `requested` a la vez construye, **en tres puntos**,
/// un manifiesto cuyo FASTA derivado es el mismo artefacto que el de la
/// referencia de origen (una reducción fallida sin enlace de origen): es la
/// única forma de que los roles esperados se solapen, porque el lector impone
/// antes todo lo demás.
#[test]
fn the_output_categories_are_judged() {
    let package = publish_rich();
    let cases = [
        set(
            "un identificador que no es texto en minimum",
            "minimum contiene un ID inválido",
            "/output_categories/minimum/0",
            json!(5),
        ),
        set(
            "un identificador vacío en minimum",
            "minimum contiene un ID inválido",
            "/output_categories/minimum/0",
            json!(""),
        ),
        case(
            "un identificador que no es texto en defaulted",
            "defaulted contiene un ID inválido",
            |manifest| push(manifest, "/output_categories/defaulted", json!(5)),
        ),
        set(
            "un identificador vacío en requested",
            "requested contiene un ID inválido",
            "/output_categories/requested/0",
            json!(""),
        ),
        case(
            "un identificador vacío en disabled",
            "disabled contiene un ID inválido",
            |manifest| push(manifest, "/output_categories/disabled", json!("")),
        ),
        case(
            "un artefacto ausente en requested",
            "la categoría requested enlaza un artefacto ausente: ghost",
            |manifest| push(manifest, "/output_categories/requested", json!("ghost")),
        ),
        case(
            "un identificador repetido en minimum",
            "el artefacto config-general se repite en minimum",
            |manifest| {
                push(
                    manifest,
                    "/output_categories/minimum",
                    json!("config-general"),
                )
            },
        ),
        case(
            "un artefacto en dos categorías",
            "el artefacto truth-vcf aparece en más de una categoría",
            |manifest| push(manifest, "/output_categories/minimum", json!("truth-vcf")),
        ),
        case(
            "un artefacto mínimo en defaulted",
            "defaulted y disabled deben estar vacías en F1",
            |manifest| move_category(manifest, "config-general", "minimum", "defaulted"),
        ),
        case(
            "un artefacto mínimo en disabled",
            "defaulted y disabled deben estar vacías en F1",
            |manifest| move_category(manifest, "config-general", "minimum", "disabled"),
        ),
        case(
            "un artefacto pedido en minimum",
            "minimum/requested no coinciden con los roles normativos de F1",
            |manifest| move_category(manifest, "truth-vcf", "requested", "minimum"),
        ),
        case(
            "un artefacto mínimo en requested",
            "minimum/requested no coinciden con los roles normativos de F1",
            |manifest| move_category(manifest, "config-general", "minimum", "requested"),
        ),
        case(
            "un artefacto en el catálogo sin categoría",
            "el catálogo contiene un artefacto sin rol de salida permitido en F1",
            |manifest| {
                let mut extra = artifact_of(manifest, "truth_vcf");
                extra["id"] = json!("extra-artifact");
                extra["path"] = json!("outputs/extra.txt");
                push(manifest, "/artifacts", extra);
            },
        ),
        case(
            "un mismo FASTA con rol de mínimo y de pedido",
            "un mismo artefacto no puede tener simultáneamente rol minimum y requested",
            |manifest| {
                manifest["references"][1]["fasta_artifact_id"] = json!("reference-fasta");
                manifest["reference_reductions"][0]["fasta_artifact_id"] = json!("reference-fasta");
                manifest["reference_reductions"][0]["source_fasta_artifact_id"] = Value::Null;
                manifest["reference_reductions"][0]["state"] = json!("failed");
            },
        ),
    ];
    assert_parse_cases(&package, &cases);
}

// ---------------------------------------------------------------------------
// Los enlaces de evento y el inventario técnico
// ---------------------------------------------------------------------------

/// **Los enlaces de evento y su ruta canónica se juzgan.**
///
/// La lista y cada enlace tienen la forma debida, el identificador es un
/// `EventId`, los enlaces van en orden **estrictamente** creciente (un
/// repetido y uno descendente se rechazan; uno creciente, con el
/// manifiesto de evento de otro, llega a la comprobación de la ruta), el
/// estado es `complete` o `failed` y la ruta del manifiesto de evento es la
/// canónica de su evento.
#[test]
fn the_event_links_are_judged() {
    let package = publish_rich();
    let invalid_id =
        genome_population_simulator::identity::EventId::parse("bad id!", "manifest.events")
            .expect_err("un identificador con espacio no es válido")
            .render_human();
    let link = |event_id: &'static str| {
        move |manifest: &mut Value| {
            let mut link = manifest["events"][0].clone();
            link["event_id"] = json!(event_id);
            push(manifest, "/events", link);
        }
    };
    let cases = [
        set(
            "events que no es una lista",
            "events no es una lista",
            "/events",
            json!(5),
        ),
        set(
            "un enlace que no es objeto",
            "un enlace de evento no es objeto",
            "/events/0",
            json!(5),
        ),
        set(
            "un event_id con sintaxis inválida",
            invalid_id.clone(),
            "/events/0/event_id",
            json!("bad id!"),
        ),
        // La primera lectura del identificador va antes que la del estado: con
        // los dos inválidos manda el identificador.
        case(
            "un event_id inválido con un estado inválido",
            invalid_id,
            |manifest| {
                manifest["events"][0]["event_id"] = json!("bad id!");
                manifest["events"][0]["state"] = json!("pending");
            },
        ),
        case(
            "un enlace repetido",
            "events no está ordenado por event_id",
            link("reads-main"),
        ),
        case(
            "un enlace descendente",
            "events no está ordenado por event_id",
            link("reads-aaa"),
        ),
        case(
            "un enlace creciente con el manifiesto de evento de otro",
            "la ruta del manifiesto de evento no es canónica",
            link("reads-zzz"),
        ),
        set(
            "un estado de evento desconocido",
            "estado de evento inválido",
            "/events/0/state",
            json!("pending"),
        ),
        case(
            "un manifiesto de evento fuera de su ruta canónica",
            "la ruta del manifiesto de evento no es canónica",
            |manifest| {
                let artifacts = manifest["artifacts"].as_array_mut().expect("catálogo");
                let record = artifacts
                    .iter_mut()
                    .find(|artifact| artifact["id"] == "event-manifest-reads-main")
                    .expect("manifiesto de evento");
                record["path"] = json!("events/reads-other/manifest.json");
            },
        ),
    ];
    assert_parse_cases(&package, &cases);
}

/// **El inventario técnico agregado se juzga.**
///
/// Cada entrada de `technical_events` es un objeto del evento que enlaza, con
/// el perfil técnico, el catálogo de azar y el catálogo de identificadores
/// admitidos; sus artefactos agregados y el resumen de quimeras existen; y
/// la configuración efectiva y el manifiesto de evento que embebe son los del
/// catálogo. **No hay casos «aceptados» para el otro perfil técnico ni para
/// los catálogos v4 y v5 en un manifiesto de versión 4**: el esquema v4 los
/// prohíbe y ningún productor los escribe (el juez completo los rechaza, ver
/// `the_general_manifest_carries_the_profile_and_the_catalog_of_the_effective_config`).
/// El catálogo v4, que el esquema v6 sí admite, se acepta en un paquete de
/// versión 6 (`a_version_6_manifest_accepts_the_random_catalogs_its_schema_admits`).
#[test]
fn the_technical_event_inventory_is_judged() {
    use genome_population_simulator::technical_profile::GPS_HISEQX_PCR_PE150_V1;
    let package = publish_rich();
    let incompatible = "inventario técnico F4 incompatible";
    let embedded = "artefacto técnico embebido no coincide con el catálogo";
    let cases = [
        set(
            "una entrada que no es objeto",
            "technical_event inválido",
            "/technical_events/0",
            json!(5),
        ),
        set(
            "otro identificador de evento",
            incompatible,
            "/technical_events/0/event_id",
            json!("reads-other"),
        ),
        set(
            "otro perfil técnico",
            incompatible,
            "/technical_events/0/profile/id",
            json!("gps-hiseqx-other-pe150-v1"),
        ),
        set(
            "la huella de otro perfil técnico",
            incompatible,
            "/technical_events/0/profile/sha256",
            json!("0".repeat(64)),
        ),
        remove_case(
            "un perfil técnico sin huella",
            incompatible,
            "/technical_events/0/profile/sha256",
        ),
        set(
            "el perfil PCR con la huella del otro",
            incompatible,
            "/technical_events/0/profile/id",
            json!(GPS_HISEQX_PCR_PE150_V1),
        ),
        set(
            "un catálogo de azar desconocido",
            incompatible,
            "/technical_events/0/rng_catalog",
            json!("f4-technical-v2"),
        ),
        set(
            "un catálogo de azar que no es texto",
            incompatible,
            "/technical_events/0/rng_catalog",
            json!(3),
        ),
        set(
            "otra versión del catálogo de identificadores",
            incompatible,
            "/technical_events/0/id_catalog/version",
            json!("4"),
        ),
        set(
            "artefactos agregados que no son objeto",
            "artefactos técnicos agregados inválidos",
            "/technical_events/0/artifacts",
            Value::Null,
        ),
        set(
            "un resumen de quimeras que no es objeto",
            "resumen de quimera inválido",
            "/technical_events/0/artifacts/chimera",
            Value::Null,
        ),
        set(
            "la configuración efectiva embebida ausente",
            embedded,
            "/technical_events/0/effective_config",
            Value::Null,
        ),
        set(
            "la configuración efectiva embebida con otro tamaño",
            embedded,
            "/technical_events/0/effective_config/size_bytes",
            json!(1),
        ),
        set(
            "el manifiesto de evento embebido con otra huella",
            embedded,
            "/technical_events/0/event_manifest/checksum/digest",
            json!("0".repeat(64)),
        ),
    ];
    assert_parse_cases(&package, &cases);
}

// ---------------------------------------------------------------------------
// Las referencias, las reducciones y los VCF
// ---------------------------------------------------------------------------

/// **Las referencias y las reducciones del manifiesto se juzgan.**
///
/// Cada referencia es un objeto con su identificador único; cada reducción
/// tiene identificador único, un modo conocido, un orden de selección de
/// enteros positivos, un estado conocido y los enlaces de sus artefactos
/// coherentes con las referencias: una reducción `complete` exige los tres, y
/// el FASTA derivado y el de origen tienen que ser los que el manifiesto
/// registra para cada referencia. Los enteros positivos (`eligible_bases`,
/// `target.bases`, `region_size_bases`) se rechazan en 0, con otro tipo y
/// ausentes; los modos y estados conocidos se aceptan.
#[test]
fn the_references_and_the_reductions_of_the_manifest_are_judged() {
    let package = publish_rich();
    let reduction = |field: &str| format!("/reference_reductions/0/{field}");
    let three_links = "una reducción complete exige sus tres enlaces de artefacto";
    let positive = |field: &'static str, pointer: String| {
        [
            set(
                format!("{field} a cero"),
                format!("{field} inválido"),
                pointer.clone(),
                json!(0),
            ),
            set(
                format!("{field} que no es número"),
                format!("{field} inválido"),
                pointer.clone(),
                json!("7"),
            ),
            remove_case(
                format!("{field} ausente"),
                format!("{field} inválido"),
                pointer,
            ),
        ]
    };
    let mut cases = vec![
        set(
            "una referencia que no es objeto",
            "una referencia no es objeto",
            "/references/1",
            json!(5),
        ),
        set(
            "una reducción que no es objeto",
            "una reducción no es objeto",
            "/reference_reductions/0",
            json!(5),
        ),
        case(
            "un identificador de referencia repetido",
            "ID de referencia repetido",
            |manifest| {
                manifest["references"][1]["id"] = json!("ref-main");
            },
        ),
        case(
            "un identificador de reducción repetido",
            "ID de reducción repetido",
            |manifest| {
                let copy = manifest["reference_reductions"][0].clone();
                push(manifest, "/reference_reductions", copy);
            },
        ),
        set(
            "un modo de reducción desconocido",
            "modo de reducción inválido",
            reduction("target/mode"),
            json!("reduction_fraction"),
        ),
        set(
            "el modo de reducción por porcentaje",
            "",
            reduction("target/mode"),
            json!("reduction_percent"),
        ),
        set(
            "un orden de selección vacío",
            "selection_order está vacío",
            reduction("selection_order"),
            json!([]),
        ),
        set(
            "un orden de selección con un cero",
            "selection_order inválido",
            reduction("selection_order"),
            json!([1, 0, 3]),
        ),
        set(
            "un orden de selección con un texto",
            "selection_order inválido",
            reduction("selection_order"),
            json!([1, "2", 3]),
        ),
        set(
            "un orden de selección con un negativo",
            "selection_order inválido",
            reduction("selection_order"),
            json!([1, -2, 3]),
        ),
        set(
            "un orden de selección que empieza en 1",
            "",
            reduction("selection_order"),
            json!([1]),
        ),
        set(
            "un estado de reducción desconocido",
            "estado de reducción inválido",
            reduction("state"),
            json!("pending"),
        ),
        set(
            "una reducción complete sin el FASTA derivado",
            three_links,
            reduction("fasta_artifact_id"),
            Value::Null,
        ),
        set(
            "una reducción complete sin el mapa",
            three_links,
            reduction("map_artifact_id"),
            Value::Null,
        ),
        set(
            "una reducción complete sin el FASTA de origen",
            three_links,
            reduction("source_fasta_artifact_id"),
            Value::Null,
        ),
        set(
            "una reducción con el mismo FASTA de origen y derivado",
            "una reducción no puede reutilizar el mismo FASTA como fuente y derivado",
            reduction("source_fasta_artifact_id"),
            json!("reference-fasta-ref-small"),
        ),
        set(
            "una reducción cuyo FASTA no es el de su referencia derivada",
            "la reducción no coincide con el enlace de su referencia derivada",
            "/references/1/fasta_artifact_id",
            json!("reference-fasta"),
        ),
        set(
            "una reducción cuyo FASTA de origen no es el de su referencia fuente",
            "la reducción no coincide con el enlace de su referencia fuente",
            "/references/0/fasta_artifact_id",
            json!("reference-fasta-ref-small"),
        ),
        set(
            "un objetivo que no es objeto",
            "target no es objeto",
            reduction("target"),
            json!(5),
        ),
    ];
    cases.extend(positive("eligible_bases", reduction("eligible_bases")));
    cases.extend(positive(
        "region_size_bases",
        reduction("region_size_bases"),
    ));
    cases.extend(positive("bases", reduction("target/bases")));
    assert_parse_cases(&package, &cases);
}

/// **Un VCF publicado se juzga.**
///
/// Como mucho uno, objeto, `complete`, de una referencia que el manifiesto
/// conoce, con muestras de texto no vacías y distintas, y su modo:
/// texto, sin índice (ni `index_kind` ni `index_artifact_id`, **ni presentes
/// con otro valor ni ausentes**), o BGZF, con el tipo de índice y el
/// artefacto de índice. El caso del VCF con una referencia duplicada no
/// existe porque el límite de uno lo rechaza antes: esa mitad de la guarda
/// es inalcanzable (anotado junto a ella).
#[test]
fn a_published_vcf_is_judged() {
    let text = publish_rich();
    let bad_text_index = "VCF de texto no admite índice";
    let cases = [
        case(
            "dos VCF",
            "solo se admite un VCF final por ejecución",
            |manifest| {
                let copy = manifest["vcf_outputs"][0].clone();
                push(manifest, "/vcf_outputs", copy);
            },
        ),
        set(
            "un VCF que no es objeto",
            "vcf_outputs contiene un valor no objeto",
            "/vcf_outputs/0",
            json!(5),
        ),
        set(
            "un VCF que no está complete",
            "un VCF publicado debe estar complete",
            "/vcf_outputs/0/state",
            json!("failed"),
        ),
        set(
            "un VCF de una referencia desconocida",
            "VCF con referencia ausente o duplicada",
            "/vcf_outputs/0/reference_id",
            json!("ref-unknown"),
        ),
        set(
            "una muestra que no es texto",
            "muestra VCF inválida",
            "/vcf_outputs/0/samples/0",
            json!(5),
        ),
        set(
            "una muestra vacía",
            "muestra VCF inválida",
            "/vcf_outputs/0/samples/0",
            json!(""),
        ),
        set(
            "una muestra repetida",
            "muestra VCF duplicada",
            "/vcf_outputs/0/samples/1",
            json!("I001"),
        ),
        set(
            "un VCF de texto con tipo de índice",
            bad_text_index,
            "/vcf_outputs/0/index_kind",
            json!("tbi"),
        ),
        set(
            "un VCF de texto con artefacto de índice",
            bad_text_index,
            "/vcf_outputs/0/index_artifact_id",
            json!("truth-vcf"),
        ),
        remove_case(
            "un VCF de texto sin la clave del tipo de índice",
            bad_text_index,
            "/vcf_outputs/0/index_kind",
        ),
        remove_case(
            "un VCF de texto sin la clave del artefacto de índice",
            bad_text_index,
            "/vcf_outputs/0/index_artifact_id",
        ),
        set(
            "un modo de VCF desconocido",
            "modo VCF inválido",
            "/vcf_outputs/0/mode",
            json!("zip"),
        ),
    ];
    assert_parse_cases(&text, &cases);

    let bgzf = publish_rich_bgzf();
    let cases = [
        set(
            "un VCF comprimido sin tipo de índice",
            "VCF BGZF exige índice TBI o CSI",
            "/vcf_outputs/0/index_kind",
            Value::Null,
        ),
        set(
            "un VCF comprimido con un tipo de índice desconocido",
            "VCF BGZF exige índice TBI o CSI",
            "/vcf_outputs/0/index_kind",
            json!("bai"),
        ),
        set(
            "un VCF comprimido sin artefacto de índice",
            "VCF BGZF exige artefacto de índice",
            "/vcf_outputs/0/index_artifact_id",
            Value::Null,
        ),
        // El lado aceptado de la lista cerrada de tipos de índice: `csi` con
        // un artefacto de su clase (el lector del manifiesto no mira la
        // longitud de las secuencias; eso lo juzga la efectiva).
        case("un VCF comprimido con índice csi", "", |manifest| {
            manifest["vcf_outputs"][0]["index_kind"] = json!("csi");
            let index = manifest["artifacts"]
                .as_array_mut()
                .expect("catálogo")
                .iter_mut()
                .find(|artifact| artifact["id"] == "truth-vcf-index")
                .expect("índice");
            index["kind"] = json!("vcf_csi_index");
        }),
        case(
            "un VCF comprimido con índice csi y un artefacto tbi",
            "el enlace truth-vcf-index exige vcf_csi_index, pero el catálogo declara vcf_tbi_index",
            |manifest| {
                manifest["vcf_outputs"][0]["index_kind"] = json!("csi");
            },
        ),
    ];
    assert_parse_cases(&bgzf, &cases);
}

// ---------------------------------------------------------------------------
// Los manifiestos de evento contra el catálogo general
// ---------------------------------------------------------------------------

/// Lo que dice `validate_published_event_manifests` de este manifiesto, con
/// los ficheros que haya ahora en el paquete.
fn judge_events(package: &Package, manifest: &Value) -> Result<(), Diagnostic> {
    validate_published_event_manifests(
        &package.destination,
        &package.manifest_path(),
        &bytes_of(manifest),
    )
}

/// Altera el manifiesto de evento del paquete, rehace su identidad en todo el
/// manifiesto general, juzga y **restaura** el fichero.
fn with_a_forged_event(package: &Package, edit: &dyn Fn(&mut Value)) -> Result<(), Diagnostic> {
    let manifest = package.json("manifest.json");
    let event_path = catalog_path(&manifest, "event-manifest-reads-main");
    let original = package.bytes(&event_path);
    let mut event: Value = serde_json::from_slice(&original).expect("evento");
    edit(&mut event);
    let bytes = bytes_of(&event);
    fs::write(package.destination.join(&event_path), &bytes).expect("evento alterado");
    let mut forged = manifest;
    rehash_artifact(&mut forged, &event_path, &bytes);
    let result = judge_events(package, &forged);
    fs::write(package.destination.join(&event_path), original).expect("restaurar");
    result
}

/// Pone otra huella o otro tamaño al artefacto que está en `path`, en todo el
/// manifiesto (catálogo y artefactos embebidos).
fn set_identity(manifest: &mut Value, path: &str, digest: Option<&str>, size: Option<u64>) {
    fn walk(value: &mut Value, path: &str, digest: Option<&str>, size: Option<u64>) {
        match value {
            Value::Object(map) => {
                if map.get("path").and_then(Value::as_str) == Some(path) {
                    if let (Some(digest), true) = (digest, map.contains_key("checksum")) {
                        map.insert(
                            "checksum".to_owned(),
                            json!({"algorithm": "sha256", "digest": digest}),
                        );
                    }
                    if let (Some(size), true) = (size, map.contains_key("size_bytes")) {
                        map.insert("size_bytes".to_owned(), json!(size));
                    }
                }
                for nested in map.values_mut() {
                    walk(nested, path, digest, size);
                }
            }
            Value::Array(items) => {
                for item in items {
                    walk(item, path, digest, size);
                }
            }
            _ => {}
        }
    }
    walk(manifest, path, digest, size);
}

/// **El juez de los manifiestos de evento ata cada uno a su identidad
/// publicada y al catálogo general.**
///
/// El manifiesto de evento tiene que ser, byte a byte, el que el catálogo
/// declara (tamaño y huella, **cada uno por separado**); su identificador y
/// su estado, los del enlace general; y el FASTA de referencia y cada salida
/// completa que nombra tienen que existir en el catálogo y ser **idénticos** a
/// su registro. Los casos de identidad cambian el registro del catálogo (y el
/// embebido) sin tocar el fichero; los demás cambian el fichero y rehacen su
/// identidad.
#[test]
fn the_event_manifests_are_judged_against_the_general_catalog() {
    let package = publish_rich();
    let manifest = package.json("manifest.json");
    judge_events(&package, &manifest).expect("control positivo");
    let event_path = catalog_path(&manifest, "event-manifest-reads-main");
    let recorded = artifact_of(&manifest, "event_manifest");
    let size = recorded["size_bytes"].as_u64().expect("tamaño");
    let identity = "el manifiesto de evento no coincide con su identidad publicada";

    // La identidad: el tamaño solo, y la huella sola.
    let mut altered = manifest.clone();
    set_identity(&mut altered, &event_path, None, Some(size + 1));
    assert_ne!(altered, manifest);
    assert_eq!(
        detail(judge_events(&package, &altered)),
        identity,
        "otro tamaño"
    );
    let mut altered = manifest.clone();
    set_identity(&mut altered, &event_path, Some(&"0".repeat(64)), None);
    assert_ne!(altered, manifest);
    assert_eq!(
        detail(judge_events(&package, &altered)),
        identity,
        "otra huella"
    );

    // El enlace: un evento que el general da por fallido y el fichero por
    // completo (la ejecución es fallida y válida, sin verdad ni VCF).
    let mut altered = manifest.clone();
    altered["state"] = json!("failed");
    altered["diagnostics"] = json!([catalog_diagnostic_of("IO")]);
    drop_artifact(&mut altered, "truth-jsonl");
    drop_artifact(&mut altered, "truth-vcf");
    altered["vcf_outputs"] = json!([]);
    altered["events"][0]["state"] = json!("failed");
    parse(&altered).expect("la ejecución fallida es un manifiesto válido");
    assert_eq!(
        detail(judge_events(&package, &altered)),
        "el enlace general no coincide con el manifiesto de evento",
        "estado del enlace distinto del del fichero"
    );

    // El enlace con otro identificador de evento: el manifiesto de evento se
    // copia a la ruta canónica del nuevo (para que el lector del general lo
    // acepte) y su catálogo y su artefacto embebido se actualizan.
    let moved = "events/reads-zzz/manifest.json";
    fs::create_dir_all(package.destination.join("events/reads-zzz")).expect("directorio");
    fs::write(package.destination.join(moved), package.bytes(&event_path)).expect("copia");
    let mut altered = manifest.clone();
    altered["events"][0]["event_id"] = json!("reads-zzz");
    altered["technical_events"][0]["event_id"] = json!("reads-zzz");
    for pointer in [
        "/artifacts/12/path",
        "/technical_events/0/event_manifest/path",
    ] {
        *altered.pointer_mut(pointer).expect("ruta") = json!(moved);
    }
    assert_eq!(altered["artifacts"][12]["id"], "event-manifest-reads-main");
    parse(&altered).expect("el enlace movido es un manifiesto válido");
    assert_eq!(
        detail(judge_events(&package, &altered)),
        "el enlace general no coincide con el manifiesto de evento",
        "identificador del enlace distinto del del fichero"
    );

    // El FASTA de referencia y las salidas: registro distinto o ausente.
    assert_eq!(
        detail(with_a_forged_event(&package, &|event| {
            event["reference"]["fasta"]["size_bytes"] = json!(1);
        })),
        "artefacto de evento no coincide con el catálogo físico general",
        "FASTA con otro tamaño"
    );
    assert_eq!(
        detail(with_a_forged_event(&package, &|event| {
            event["reference"]["fasta"]["id"] = json!("reference-fasta-unknown");
        })),
        "artefacto de evento no catalogado",
        "FASTA fuera del catálogo"
    );
    assert_eq!(
        detail(with_a_forged_event(&package, &|event| {
            event["outputs"][0]["artifact"]["size_bytes"] = json!(1);
        })),
        "artefacto de evento no coincide con el catálogo físico general",
        "salida completa con otro tamaño"
    );
    assert_eq!(
        detail(with_a_forged_event(&package, &|event| {
            event["outputs"][0]["artifact"]["id"] = json!("fastq-unknown");
        })),
        "artefacto FASTQ completo incoherente",
        "salida completa fuera del catálogo: la rechaza antes el juez del evento"
    );
}

/// **Intento ejecutado de alcanzar las guardas del juez de los manifiestos de
/// evento que el juez propio del manifiesto de evento sombrea.** Cada caso
/// altera el fichero del evento en el punto que la guarda comprueba (y rehace
/// su identidad) y comprueba que el rechazo viene del juez del manifiesto de
/// evento, con su texto. **No discrimina las guardas anotadas**: si se borrara
/// una de ellas la prueba seguiría pasando, porque el juez del evento rechaza
/// antes; solo demuestra que esa comprobación anterior existe y rechaza lo
/// que la guarda no vería. Las que quedan sin caso aquí (la lectura del
/// manifiesto general, que ya hizo `parse_general_manifest_bytes`, y la lista
/// `events`, que ya validó `validate_manifest_events`) no tienen ni siquiera
/// intento: el lector del general las rechaza antes de que este juez llegue.
#[test]
fn unreachable_event_catalog_guards_are_shadowed_by_the_event_judge() {
    let package = publish_rich();
    let attempts: Vec<Attempt> = vec![
        (
            "un manifiesto de evento sin referencia",
            Box::new(|event| {
                event.as_object_mut().expect("objeto").remove("reference");
            }),
            "campos inválidos en raíz",
        ),
        (
            "una referencia sin FASTA",
            Box::new(|event| {
                event["reference"]
                    .as_object_mut()
                    .expect("objeto")
                    .remove("fasta");
            }),
            "campos inválidos en reference",
        ),
        (
            "un FASTA que no es objeto",
            Box::new(|event| event["reference"]["fasta"] = json!(5)),
            "artifact no es objeto",
        ),
        (
            "un FASTA sin identificador",
            Box::new(|event| {
                event["reference"]["fasta"]
                    .as_object_mut()
                    .expect("objeto")
                    .remove("id");
            }),
            "campos inválidos en artifact",
        ),
        (
            "un manifiesto de evento sin salidas",
            Box::new(|event| {
                event.as_object_mut().expect("objeto").remove("outputs");
            }),
            "campos inválidos en raíz",
        ),
        (
            "una salida completa sin artefacto",
            Box::new(|event| {
                event["outputs"][0]
                    .as_object_mut()
                    .expect("objeto")
                    .remove("artifact");
            }),
            "campos inválidos en outputs[]",
        ),
        (
            "una salida completa con un artefacto sin identificador",
            Box::new(|event| {
                event["outputs"][0]["artifact"]
                    .as_object_mut()
                    .expect("objeto")
                    .remove("id");
            }),
            "campos inválidos en artifact",
        ),
        (
            "una salida completa con un artefacto que no es objeto",
            Box::new(|event| event["outputs"][0]["artifact"] = json!(5)),
            "artifact no puede ser nulo",
        ),
    ];
    for (label, edit, expected) in attempts {
        assert_eq!(
            detail(with_a_forged_event(&package, &*edit)),
            expected,
            "«{label}» se rechazó por otra razón"
        );
    }
}

// ---------------------------------------------------------------------------
// resolve_manifest_reference
// ---------------------------------------------------------------------------

/// **Resolver una referencia desde su manifiesto.**
///
/// Control positivo con la referencia suministrada y con la reducida (el
/// juez completo del paquete corre dentro). Una referencia que el manifiesto
/// no tiene se rechaza con su nombre; y una ejecución fallida que no llegó a
/// publicar su configuración efectiva (sin enlace a ella) se rechaza **después
/// de resolver el FASTA**, porque la identidad del paquete no se puede juzgar
/// sin ella. Las demás guardas de la función (sin directorio padre, JSON,
/// `references`, FASTA y artefacto de configuración efectiva de clase
/// distinta) las rechaza antes el lector del manifiesto: están anotadas junto
/// a ellas.
#[test]
fn a_reference_is_resolved_from_its_manifest() {
    use genome_population_simulator::identity::ReferenceId;
    use genome_population_simulator::results::resolve_manifest_reference;
    let package = publish_rich();
    let id = |text: &str| ReferenceId::parse(text, "test").expect("identificador");
    for reference in ["ref-main", "ref-small"] {
        let resolved = resolve_manifest_reference(&package.manifest_path(), &id(reference), None)
            .unwrap_or_else(|error| panic!("{reference} se tenía que resolver: {error:?}"));
        assert!(
            resolved.fasta_path().starts_with(&package.destination),
            "{reference}"
        );
    }
    assert_eq!(
        detail_of(resolve_manifest_reference(
            &package.manifest_path(),
            &id("ref-unknown"),
            None
        )),
        "la referencia solicitada no existe: ref-unknown"
    );

    // Una ejecución fallida sin enlace a su configuración efectiva.
    let simulated = publish_simulated();
    let mut failed = simulated.json("manifest.json");
    failed["state"] = json!("failed");
    failed["diagnostics"] = json!([catalog_diagnostic_of("IO")]);
    failed["effective_config_artifact_id"] = Value::Null;
    drop_artifact(&mut failed, "truth-jsonl");
    drop_artifact(&mut failed, "effective-config");
    parse(&failed).expect("la ejecución fallida es un manifiesto válido");
    fs::write(simulated.manifest_path(), bytes_of(&failed)).expect("manifiesto fallido");
    assert_eq!(
        detail_of(resolve_manifest_reference(
            &simulated.manifest_path(),
            &id("ref-sim"),
            None
        )),
        "falta effective_config_artifact_id"
    );
}

/// El rechazo de una función que devuelve un valor, como texto.
fn detail_of<T>(result: Result<T, Diagnostic>) -> String {
    match result {
        Ok(_) => panic!("el juez tenía que rechazar"),
        Err(error) => error.context().get("detail").cloned().unwrap_or_default(),
    }
}

// ---------------------------------------------------------------------------
// Ayudantes de los jueces físicos: sustituir ficheros del paquete y juzgar
// ---------------------------------------------------------------------------

type Edit = Box<dyn Fn(&mut Value)>;

fn full_judge(package: &Package, effective: &Value, manifest: &Value) -> Result<(), Diagnostic> {
    validate_published_resources(
        &package.destination,
        &package.effective_path(),
        &bytes_of(effective),
        &package.manifest_path(),
        &bytes_of(manifest),
        None,
    )
}

/// Lo que dice el juez completo con la efectiva, el manifiesto de evento y el
/// manifiesto general alterados; el manifiesto de evento se restaura después.
fn probe(
    package: &Package,
    effective_edit: &dyn Fn(&mut Value),
    event_edit: &dyn Fn(&mut Value),
    manifest_edit: &dyn Fn(&mut Value),
) -> String {
    let mut manifest = package.json("manifest.json");
    let mut effective = package.json("effective-config.json");
    let event_path = catalog_path(&manifest, "event-manifest-reads-main");
    let original = package.bytes(&event_path);
    let mut event: Value = serde_json::from_slice(&original).expect("evento");
    effective_edit(&mut effective);
    event_edit(&mut event);
    let bytes = bytes_of(&event);
    fs::write(package.destination.join(&event_path), &bytes).expect("evento");
    rehash_artifact(&mut manifest, &event_path, &bytes);
    manifest_edit(&mut manifest);
    let result = full_judge(package, &effective, &manifest);
    fs::write(package.destination.join(&event_path), &original).expect("restaurar");
    match result {
        Ok(()) => "OK".to_owned(),
        Err(error) => error
            .context()
            .get("detail")
            .cloned()
            .unwrap_or_else(|| format!("{error:?}")),
    }
}

/// Como `probe`, pero con cualquier fichero del paquete sustituido (y, si se
/// pide, su identidad rehecha en el manifiesto). Los ficheros se restauran.
fn probe_files(
    package: &Package,
    effective_edit: &dyn Fn(&mut Value),
    manifest_edit: &dyn Fn(&mut Value),
    files: &[(String, Vec<u8>)],
    rehash: bool,
) -> Result<(), Diagnostic> {
    let mut manifest = package.json("manifest.json");
    let mut effective = package.json("effective-config.json");
    let mut originals = Vec::new();
    for (relative, bytes) in files {
        originals.push((relative.clone(), package.bytes(relative)));
        fs::write(package.destination.join(relative), bytes).expect("fichero alterado");
        if rehash {
            rehash_artifact(&mut manifest, relative, bytes);
        }
    }
    effective_edit(&mut effective);
    manifest_edit(&mut manifest);
    let result = full_judge(package, &effective, &manifest);
    for (relative, bytes) in originals {
        fs::write(package.destination.join(relative), bytes).expect("restaurar");
    }
    result
}

fn describe(result: Result<(), Diagnostic>) -> String {
    match result {
        Ok(()) => "OK".to_owned(),
        Err(error) => format!(
            "{} {}",
            error.code(),
            error
                .context()
                .get("detail")
                .cloned()
                .unwrap_or_else(|| format!("{:?}", error.context()))
        ),
    }
}

/// Dos receptoras de contaminación (I001 e I002, de la cohorte P1) y un tercer
/// individuo técnico de otra cohorte (J001); I003, de P1, no se secuencia.
fn publish_contaminating() -> Package {
    let toml = rich_toml("text", "")
        .replace(
            "founder_keys = [\"fk-a\", \"fk-b\"]\nploidy = 2\ngenotype_model = \"f_diploid\"\ninbreeding_coefficient = 0.5",
            "individual_ids = [\"J001\", \"J002\"]\nploidy = 2\ngenotype_model = \"hwe\"",
        )
        + r#"
[technical_event]
event_id = "reads-main"
individual_ids = ["I001", "I002", "J001"]
target_coverage_wgs = 1
contamination_min_percent = 10
contamination_max_percent = 20
contamination_recipient_ids = ["I001", "I002"]

[[technical_event.libraries]]
individual_id = "I001"
i7 = "ACGTACGT"
i5 = "TGCATGCA"

[[technical_event.libraries]]
individual_id = "I002"
i7 = "TTGGCCAA"
i5 = "GGAATTCC"

[[technical_event.libraries]]
individual_id = "J001"
i7 = "CCGGTTAA"
i5 = "AATTGGCC"
"#;
    publish(&[("general.toml", &toml), ("reference.fa", &rich_reference())])
}

/// Una referencia suministrada **sin nada más**: el juez solo comprueba su
/// ficha y su alfabeto, sin cargar las bases.
fn publish_summary() -> Package {
    let toml = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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
"#
    );
    publish(&[("general.toml", &toml), ("reference.fa", &rich_reference())])
}

/// El rechazo de una autoridad física, tal como lo cita el juez.
fn physical(operation: &str, code: &str, meaning: &str) -> String {
    format!("autoridad física rechazada en {operation}: {code} {meaning}")
}

/// El rechazo con el que `artifact_identity_error` marca un artefacto cuyo
/// tamaño o huella no son los del catálogo, con el identificador que cita.
fn identity_rejection(result: Result<(), Diagnostic>) -> (String, String, String) {
    let error = result.expect_err("el juez tenía que rechazar");
    (
        error.code().to_owned(),
        error
            .context()
            .get("artifact_id")
            .cloned()
            .unwrap_or_default(),
        error.context().get("path").cloned().unwrap_or_default(),
    )
}

/// **Los recursos que el juez completo rehace por su cuenta.**
///
/// - El orden de selección de una reducción tiene que ser el que sale de la
///   subsemilla (la efectiva y el manifiesto, con el mismo orden permutado,
///   pasan todo lo demás).
/// - El mapa de reducción tiene que ser JSON legible.
/// - El catálogo físico derivado tiene que tener las secuencias del plan
///   (una reducción con sus dos secuencias fundidas en una, con la efectiva y
///   el FASTA a juego, tiene la misma suma de bases).
/// - Cada artefacto VCF publicado tiene el tamaño y la huella del catálogo
///   (**cada uno por separado**).
/// - Una referencia que solo se comprueba por su ficha se lee del FASTA
///   publicado.
#[test]
fn the_physical_resources_are_rebuilt_and_their_identity_is_kept() {
    let package = publish_rich();
    let none: &dyn Fn(&mut Value) = &|_| {};
    let manifest = package.json("manifest.json");
    let map_path = catalog_path(&manifest, "reference-map-ref-small");
    let small_path = catalog_path(&manifest, "reference-fasta-ref-small");
    let vcf_path = catalog_path(&manifest, "truth-vcf");

    assert_eq!(
        describe(probe_files(&package, none, none, &[], false)),
        "OK",
        "control positivo"
    );

    // El orden de selección.
    let permuted = json!([2, 1, 3, 4, 5, 6, 7]);
    assert_ne!(
        package.json("effective-config.json")["reduced_references"][0]["selection_order"],
        permuted
    );
    assert_eq!(
        describe(probe_files(
            &package,
            &|effective| effective["reduced_references"][0]["selection_order"] = permuted.clone(),
            &|manifest| manifest["reference_reductions"][0]["selection_order"] = permuted.clone(),
            &[],
            false,
        )),
        "GPS-FMT-006 selection_order de ref-small no coincide con la subsemilla"
    );

    // El mapa que no es JSON.
    let not_json = b"{not json".to_vec();
    let serde_error = serde_json::from_slice::<Value>(&not_json).expect_err("no es JSON");
    assert_eq!(
        describe(probe_files(
            &package,
            none,
            none,
            &[(map_path, not_json)],
            true
        )),
        format!("GPS-FMT-006 mapa reference-map-ref-small no es JSON legible: {serde_error}")
    );

    // El catálogo derivado con otras secuencias que las del plan.
    let merged = format!(">chr1\n{}\n", "ACGT".repeat(175));
    assert_eq!(
        describe(probe_files(
            &package,
            &|effective| {
                let sequences = effective["reduced_references"][0]["sequences"]
                    .as_array_mut()
                    .expect("secuencias");
                sequences.pop();
                sequences[0]["length"] = json!(700);
            },
            none,
            &[(small_path, merged.into_bytes())],
            true,
        )),
        "GPS-FMT-006 el catálogo físico derivado no coincide con el plan"
    );

    // La identidad del VCF de texto: otro tamaño, y misma longitud con otra huella.
    let original = package.bytes(&vcf_path);
    let mut longer = original.clone();
    longer.extend_from_slice(b"x\n");
    let mut other_digest = original.clone();
    let last = other_digest.len() - 2;
    other_digest[last] = if other_digest[last] == b'0' {
        b'1'
    } else {
        b'0'
    };
    assert_eq!(other_digest.len(), original.len());
    for (label, bytes) in [("otro tamaño", longer), ("otra huella", other_digest)] {
        assert_eq!(
            identity_rejection(probe_files(
                &package,
                none,
                none,
                &[(vcf_path.clone(), bytes)],
                false
            )),
            (
                "GPS-IO-002".to_owned(),
                "truth-vcf".to_owned(),
                package.destination.join(&vcf_path).display().to_string()
            ),
            "VCF de texto con {label}"
        );
    }

    // La del índice del VCF comprimido.
    let bgzf = publish_rich_bgzf();
    let bgzf_manifest = bgzf.json("manifest.json");
    let index_path = catalog_path(&bgzf_manifest, "truth-vcf-index");
    let original = bgzf.bytes(&index_path);
    let mut longer = original.clone();
    longer.push(0);
    let mut other_digest = original.clone();
    other_digest[0] ^= 0xff;
    for (label, bytes) in [("otro tamaño", longer), ("otra huella", other_digest)] {
        assert_eq!(
            identity_rejection(probe_files(
                &bgzf,
                none,
                none,
                &[(index_path.clone(), bytes)],
                false
            )),
            (
                "GPS-IO-002".to_owned(),
                "truth-vcf-index".to_owned(),
                bgzf.destination.join(&index_path).display().to_string()
            ),
            "índice con {label}"
        );
    }

    // El tamaño declarado solo: el fichero es el publicado (su huella es la del
    // catálogo) y es el catálogo el que declara otro tamaño.
    for (package, id, path) in [
        (&package, "truth-vcf", &vcf_path),
        (&bgzf, "truth-vcf-index", &index_path),
    ] {
        let size = package.bytes(path).len() as u64;
        assert_eq!(
            identity_rejection(probe_files(
                package,
                none,
                &|manifest| set_identity(manifest, path, None, Some(size + 1)),
                &[],
                false
            )),
            (
                "GPS-IO-002".to_owned(),
                id.to_owned(),
                package.destination.join(path).display().to_string()
            ),
            "{id} con otro tamaño declarado"
        );
    }

    // La ruta canónica: el catálogo enlaza el VCF (y su índice) a **otro
    // fichero**, con los mismos bytes y, por tanto, la misma identidad; el
    // canónico queda intacto. El lector del manifiesto no exige la ruta y la
    // identidad se comprueba sobre la canónica, así que solo esta guarda lo ve.
    for (package, id, path) in [
        (&package, "truth-vcf", &vcf_path),
        (&bgzf, "truth-vcf-index", &index_path),
    ] {
        let moved = format!("outputs/moved-{id}");
        fs::write(package.destination.join(&moved), package.bytes(path)).expect("copia");
        let result = probe_files(
            package,
            none,
            &|manifest| {
                for artifact in manifest["artifacts"].as_array_mut().expect("catálogo") {
                    if artifact["id"] == id {
                        artifact["path"] = json!(moved);
                    }
                }
            },
            &[],
            false,
        );
        fs::remove_file(package.destination.join(&moved)).expect("borrar la copia");
        assert_eq!(
            describe(result),
            format!("GPS-FMT-006 el artefacto {id} no tiene clase o ruta VCF canónica"),
            "{id} enlazado a otro fichero"
        );
    }

    // Una referencia que solo se comprueba por su ficha.
    let summary = publish_summary();
    let summary_manifest = summary.json("manifest.json");
    let fasta_path = catalog_path(&summary_manifest, "reference-fasta");
    assert_eq!(
        describe(probe_files(&summary, none, none, &[], false)),
        "OK",
        "control positivo de la referencia por su ficha"
    );
    assert_eq!(
        describe(probe_files(
            &summary,
            none,
            none,
            &[(fasta_path, b">chr1 primary\nACGT\n".to_vec())],
            false
        )),
        format!(
            "GPS-FMT-006 {}",
            physical(
                "validate_reference_fasta",
                "GPS-CONF-008",
                "La configuración contiene una combinación de campos incompatible."
            )
        )
    );
}

/// **La cobertura del evento es la que la configuración efectiva pide.**
///
/// Cada entrada de cobertura del manifiesto de evento tiene el objetivo de
/// cobertura de la efectiva **y** como denominador la longitud total de la
/// referencia, cada uno por separado. El denominador se altera en las dos
/// entradas, con las razones de bases propias y contaminantes rehechas a
/// juego (si no, lo rechaza antes el juez del manifiesto de evento, que ata
/// el resumen de FASTQ a la cobertura final).
#[test]
fn the_event_coverage_is_the_one_the_effective_config_asks_for() {
    let package = publish_rich();
    let none: &dyn Fn(&mut Value) = &|_| {};
    assert_eq!(probe(&package, none, none, none), "OK", "control positivo");
    assert_eq!(
        probe(
            &package,
            none,
            &|event| event["coverage"][0]["target_x"] = json!(2.0),
            none
        ),
        "la cobertura F4 no coincide con el objetivo o denominador estructural efectivo",
        "otro objetivo de cobertura"
    );
    assert_eq!(
        probe(
            &package,
            none,
            &|event| {
                for entry in 0..2 {
                    let own = event["final_fastq_summary"][entry]["own_base_count"]
                        .as_u64()
                        .expect("bases propias") as f64;
                    let contaminant = event["final_fastq_summary"][entry]["contaminant_base_count"]
                        .as_u64()
                        .expect("bases contaminantes") as f64;
                    event["coverage"][entry]["denominator_bases"] = json!(1401);
                    event["coverage"][entry]["final_x"]["own_x"] = json!(own / 1401.0);
                    event["coverage"][entry]["final_x"]["contaminant_x"] =
                        json!(contaminant / 1401.0);
                }
            },
            none
        ),
        "la cobertura F4 no coincide con el objetivo o denominador estructural efectivo",
        "otro denominador"
    );
}

/// Cambia un punto del registro de contaminación de una receptora.
fn contamination(index: usize, tail: &'static str, value: Value) -> Box<dyn Fn(&mut Value)> {
    Box::new(move |event| {
        *event
            .pointer_mut(&format!(
                "/artifacts/biological_contamination/{index}/{tail}"
            ))
            .unwrap_or_else(|| panic!("el fixture no tiene {tail}")) = value.clone();
    })
}

/// **El resumen de contaminación biológica se juzga contra la configuración
/// efectiva.** Paquete con dos receptoras (I001 e I002, de la cohorte P1) y un
/// tercer individuo técnico de otra cohorte (J001); I003, de P1, no se
/// secuencia.
///
/// Un rango positivo exige un resumen por receptora; el rango de cada
/// resumen es el de la efectiva (mínimo y máximo **por separado**); un
/// donante tiene que estar secuenciado **y** ser de la población de la
/// receptora (los dos casos por separado: I003 es de la población pero no
/// está secuenciado; J001 está secuenciado pero es de otra población); y el
/// resumen cubre **exactamente** las receptoras efectivas.
#[test]
fn the_contamination_summary_is_judged_against_the_effective_config() {
    let package = publish_contaminating();
    let none: &dyn Fn(&mut Value) = &|_| {};
    assert_eq!(probe(&package, none, none, none), "OK", "control positivo");
    let outside = "donante contaminante fuera de la población seleccionada";
    let range = "rango del resumen de contaminación no coincide con la configuración efectiva";
    assert_eq!(
        probe(
            &package,
            none,
            &|event| event["artifacts"]["biological_contamination"] = json!([]),
            none
        ),
        "un rango positivo de contaminación exige resumen por receptora",
        "ningún resumen con un rango positivo"
    );
    // El mínimo de la segunda receptora (su porcentaje muestreado, 13,49, cabe).
    assert_eq!(
        probe(
            &package,
            none,
            &*contamination(1, "percent_range/minimum", json!(11.0)),
            none
        ),
        range,
        "otro mínimo"
    );
    assert_eq!(
        probe(
            &package,
            none,
            &*contamination(0, "percent_range/maximum", json!(21.0)),
            none
        ),
        range,
        "otro máximo"
    );
    assert_eq!(
        probe(
            &package,
            none,
            &*contamination(0, "donor_contributions/0/individual_id", json!("I003")),
            none
        ),
        outside,
        "un donante de la población que no se secuencia"
    );
    assert_eq!(
        probe(
            &package,
            none,
            &*contamination(0, "donor_contributions/0/individual_id", json!("J001")),
            none
        ),
        outside,
        "un donante secuenciado de otra población"
    );
    assert_eq!(
        probe(
            &package,
            none,
            &|event| {
                event["artifacts"]["biological_contamination"]
                    .as_array_mut()
                    .expect("resúmenes")
                    .pop();
            },
            none
        ),
        "el resumen de contaminación no cubre exactamente las receptoras efectivas",
        "una receptora sin resumen"
    );
}

/// **El enlace del VCF y las autorizaciones de la configuración efectiva se
/// juzgan contra el manifiesto.**
///
/// El `vcf_output` efectivo es un objeto con un modo conocido; si pide un VCF
/// el manifiesto lo enlaza **una vez** (un manifiesto válido sin VCF, con la
/// efectiva pidiéndolo, se rechaza), y un VCF comprimido enlaza el índice del
/// tipo y con el identificador que corresponden a la referencia (cada uno por
/// separado: el tipo `csi` con su artefacto, y el artefacto con otro
/// identificador); las autorizaciones son las mismas en los dos documentos.
///
/// **No hay caso para la tercera mitad de la comprobación del índice**
/// (`catalog.get(VCF_INDEX_ARTIFACT_ID) != Some(&artifact_kind)`): el lector
/// del manifiesto ya exige que la clase del artefacto de índice sea la del
/// `index_kind`, y el tipo se compara antes; es implicada por las otras dos.
#[test]
fn the_effective_vcf_link_and_the_authorizations_are_judged_against_the_manifest() {
    let text = publish_rich();
    let none: &dyn Fn(&mut Value) = &|_| {};
    assert_eq!(probe(&text, none, none, none), "OK", "control positivo");
    assert_eq!(
        probe(&text, &|e| e["vcf_output"] = json!(5), none, none),
        "vcf_output inválido"
    );
    assert_eq!(
        probe(
            &text,
            &|e| e["vcf_output"]["mode"] = json!("zip"),
            none,
            none
        ),
        "vcf_output.mode inválido"
    );
    assert_eq!(
        probe(&text, none, none, &|manifest| {
            manifest["vcf_outputs"] = json!([]);
            drop_artifact(manifest, "truth-vcf");
        }),
        "la configuración efectiva VCF exige exactamente un enlace en el manifiesto",
        "la efectiva pide un VCF y el manifiesto, válido, no lo enlaza"
    );
    assert_eq!(
        probe(&text, none, none, &|manifest| {
            let copy = manifest["vcf_outputs"][0].clone();
            push(manifest, "/vcf_outputs", copy);
        }),
        "solo se admite un VCF final por ejecución",
        "dos enlaces: lo rechaza antes el lector del manifiesto"
    );
    assert_eq!(
        probe(&text, &|e| e["vcf_output"] = Value::Null, none, none),
        "el manifiesto declara VCF aunque la configuración efectiva lo desactiva",
        "la efectiva no pide VCF y el manifiesto lo enlaza"
    );
    assert_eq!(
        probe(
            &text,
            &|e| e["authorizations"]["allow_empty_snv_matrix"] = json!(true),
            none,
            none
        ),
        "authorizations no coincide con el manifiesto",
        "autorizaciones distintas"
    );

    let bgzf = publish_rich_bgzf();
    assert_eq!(
        probe(&bgzf, none, none, none),
        "OK",
        "control positivo (BGZF)"
    );
    let index_mismatch =
        "índice VCF del manifiesto no coincide con referencia y configuración efectiva";
    assert_eq!(
        probe(&bgzf, none, none, &|manifest| {
            manifest["vcf_outputs"][0]["index_kind"] = json!("csi");
            let index = manifest["artifacts"]
                .as_array_mut()
                .expect("catálogo")
                .iter_mut()
                .find(|artifact| artifact["id"] == "truth-vcf-index")
                .expect("índice");
            index["kind"] = json!("vcf_csi_index");
        }),
        index_mismatch,
        "un índice csi para una referencia que cabe en tbi"
    );
    assert_eq!(
        probe(&bgzf, none, none, &|manifest| {
            manifest["vcf_outputs"][0]["index_artifact_id"] = json!("other-index");
            for pointer in ["/artifacts", "/output_categories/requested"] {
                for entry in manifest
                    .pointer_mut(pointer)
                    .and_then(Value::as_array_mut)
                    .expect("lista")
                {
                    if entry == "truth-vcf-index" {
                        *entry = json!("other-index");
                    } else if entry["id"] == "truth-vcf-index" {
                        entry["id"] = json!("other-index");
                    }
                }
            }
        }),
        index_mismatch,
        "el índice con otro identificador"
    );
}

fn probe_effective(package: &Package, edit: &dyn Fn(&mut Value)) -> String {
    let manifest = package.json("manifest.json");
    let mut effective = package.json("effective-config.json");
    edit(&mut effective);
    match full_judge(package, &effective, &manifest) {
        Ok(()) => "OK".to_owned(),
        Err(error) => error.context().get("detail").cloned().unwrap_or_default(),
    }
}

/// **Intento ejecutado de alcanzar las guardas físicas que otras comprobaciones
/// sombrean.** Cada caso altera lo que la guarda comprueba y exige que el
/// rechazo venga de **esa otra** comprobación, con su texto (distinto, salvo
/// donde se dice). **No discrimina las guardas anotadas**: si se borrara una
/// de ellas la prueba seguiría pasando; solo demuestra que la comprobación
/// anterior existe y rechaza lo que la guarda no vería. Lo que dice cada
/// grupo, junto a la guarda en `src/results.rs`:
///
/// - Lo que la configuración efectiva declara del evento técnico, de la
///   contaminación, de las poblaciones y de las reducciones lo ha validado ya
///   `validate_effective_config_v3_against_manifest_bytes`.
/// - El manifiesto de evento publicado lo ha validado ya su propio juez
///   (`validate_event_manifest_bytes`), con las mismas cuentas: tipos de los
///   recuentos, el redondeo del objetivo, la conservación de las parejas, los
///   aportes de cada donante y su suma.
/// - El manifiesto general lo ha leído ya `parse_general_manifest_bytes`.
#[test]
fn unreachable_physical_guards_are_shadowed_by_other_checks() {
    type Triple = (&'static str, Edit, Edit, Edit);
    let package = publish_rich();
    let none = || -> Edit { Box::new(|_| {}) };
    let contamination = |tail: &'static str, value: Value| -> Edit {
        Box::new(move |event: &mut Value| {
            *event
                .pointer_mut(&format!("/artifacts/biological_contamination/0/{tail}"))
                .expect("punto") = value.clone();
        })
    };
    let effective = |pointer: &'static str, value: Value| -> Edit {
        Box::new(move |document: &mut Value| {
            *document.pointer_mut(pointer).expect("punto") = value.clone();
        })
    };
    let attempts: Vec<(Triple, &str)> = vec![
        // El evento técnico efectivo.
        (
            (
                "technical_event que no es objeto",
                effective("/technical_event", json!(5)),
                none(),
                none(),
            ),
            "technical_event no es objeto",
        ),
        (
            (
                "event_id técnico inválido (mismo texto que la guarda)",
                effective("/technical_event/event_id", json!("bad id!")),
                none(),
                none(),
            ),
            "event_id técnico inválido",
        ),
        (
            (
                "individual_id técnico inválido (mismo texto que la guarda)",
                effective("/technical_event/individual_ids/0", json!(5)),
                none(),
                none(),
            ),
            "individual_id técnico inválido",
        ),
        (
            (
                "una biblioteca efectiva menos que individuos",
                Box::new(|e: &mut Value| {
                    e["technical_event"]["libraries"]
                        .as_array_mut()
                        .expect("bibliotecas")
                        .pop();
                }),
                none(),
                none(),
            ),
            "libraries no coincide con individual_ids",
        ),
        (
            (
                "sin enlace de evento en el manifiesto",
                none(),
                none(),
                Box::new(|m: &mut Value| m["events"] = json!([])),
            ),
            "technical_events no coincide con events",
        ),
        (
            (
                "un enlace de evento a un artefacto que no existe",
                none(),
                none(),
                Box::new(|m: &mut Value| m["events"][0]["manifest_artifact_id"] = json!("x")),
            ),
            "el enlace x no existe",
        ),
        // La contaminación efectiva y las poblaciones.
        (
            (
                "mínimo de contaminación que no es número",
                effective(
                    "/technical_event/biological_contamination/percent_range/minimum",
                    json!("x"),
                ),
                none(),
                none(),
            ),
            "rango de contaminación biológica inválido",
        ),
        (
            (
                "máximo de contaminación que no es número",
                effective(
                    "/technical_event/biological_contamination/percent_range/maximum",
                    json!("x"),
                ),
                none(),
                none(),
            ),
            "rango de contaminación biológica inválido",
        ),
        (
            (
                "mínimo de contaminación mayor que el máximo",
                effective(
                    "/technical_event/biological_contamination/percent_range/minimum",
                    json!(30.0),
                ),
                none(),
                none(),
            ),
            "rango de contaminación biológica inválido",
        ),
        (
            (
                "receptora que no es texto",
                effective(
                    "/technical_event/biological_contamination/recipient_ids/0",
                    json!(5),
                ),
                none(),
                none(),
            ),
            "receptora de contaminación biológica inválida",
        ),
        (
            (
                "biblioteca técnica que no es objeto",
                effective("/technical_event/libraries/0", json!(5)),
                none(),
                none(),
            ),
            "library inválida",
        ),
        (
            (
                "cohorte que no es objeto",
                effective("/initial_populations/0", json!(5)),
                none(),
                none(),
            ),
            "initial_populations[] inválida",
        ),
        (
            (
                "individuo de una cohorte que no es texto",
                effective("/initial_populations/0/individual_ids/0", json!(5)),
                none(),
                none(),
            ),
            "individual_id inicial inválido",
        ),
        // El resumen de contaminación publicado.
        (
            (
                "sin resumen de contaminación",
                none(),
                Box::new(|v: &mut Value| {
                    v["artifacts"]
                        .as_object_mut()
                        .expect("artefactos")
                        .remove("biological_contamination");
                }),
                none(),
            ),
            "campos inválidos en artifacts",
        ),
        (
            (
                "un resumen que no es objeto",
                none(),
                Box::new(|v: &mut Value| v["artifacts"]["biological_contamination"][0] = json!(5)),
                none(),
            ),
            "artifacts.biological_contamination[] no es objeto",
        ),
        (
            (
                "un porcentaje muestreado que no es número",
                none(),
                contamination("sampled_percent", json!("x")),
                none(),
            ),
            "porcentaje de contaminación biológica inválido",
        ),
        (
            (
                "un recuento limpio que no es entero",
                none(),
                contamination("clean_pair_count", json!("x")),
                none(),
            ),
            "recuento limpio de contaminación inválido",
        ),
        (
            (
                "un objetivo que no es entero",
                none(),
                contamination("target_pair_count", json!("x")),
                none(),
            ),
            "objetivo de contaminación inválido",
        ),
        (
            (
                "un recuento retirado que no es entero",
                none(),
                contamination("removed_own_pair_count", json!("x")),
                none(),
            ),
            "recuento retirado de contaminación inválido",
        ),
        (
            (
                "un recuento insertado que no es entero",
                none(),
                contamination("inserted_contaminant_pair_count", json!("x")),
                none(),
            ),
            "recuento insertado de contaminación inválido",
        ),
        (
            (
                "un recuento retenido que no es entero",
                none(),
                contamination("retained_own_pair_count", json!("x")),
                none(),
            ),
            "recuento retenido de contaminación inválido",
        ),
        (
            (
                "un recuento final que no es entero",
                none(),
                contamination("final_pair_count", json!("x")),
                none(),
            ),
            "recuento final de contaminación inválido",
        ),
        (
            (
                "una receptora con identificador inválido",
                none(),
                contamination("recipient_id", json!("bad id!")),
                none(),
            ),
            "ID individual inválido",
        ),
        (
            (
                "una biblioteca con identificador inválido",
                none(),
                contamination("recipient_library_id", json!("bad id!")),
                none(),
            ),
            "ID library inválido",
        ),
        (
            (
                "un objetivo que no es el redondeo normativo",
                none(),
                contamination("target_pair_count", json!(2)),
                none(),
            ),
            "objetivo de contaminación no aplica el redondeo normativo",
        ),
        (
            (
                "contadores que no conservan las parejas",
                none(),
                contamination("final_pair_count", json!(6)),
                none(),
            ),
            "los contadores realizados de contaminación no conservan las parejas",
        ),
        (
            (
                "un aporte que no es objeto",
                none(),
                contamination("donor_contributions/0", json!(5)),
                none(),
            ),
            "donor_contributions[] no es objeto",
        ),
        (
            (
                "unos pares de aporte que no son entero",
                none(),
                contamination("donor_contributions/0/pairs", json!("x")),
                none(),
            ),
            "aporte de donante inválido",
        ),
        (
            (
                "unos aportes que no suman el objetivo",
                none(),
                contamination("donor_contributions/0/pairs", json!(2)),
                none(),
            ),
            "los aportes no suman el objetivo de contaminación",
        ),
        (
            (
                "la receptora como donante de sí misma",
                none(),
                contamination("donor_contributions/0/individual_id", json!("I002")),
                none(),
            ),
            "aporte de donante inválido",
        ),
        // La cobertura del evento.
        (
            (
                "sin cobertura en el evento",
                none(),
                Box::new(|v: &mut Value| {
                    v.as_object_mut().expect("evento").remove("coverage");
                }),
                none(),
            ),
            "campos inválidos en raíz",
        ),
        (
            (
                "una entrada de cobertura que no es objeto",
                none(),
                Box::new(|v: &mut Value| v["coverage"][0] = json!(5)),
                none(),
            ),
            "coverage[] no es objeto",
        ),
        (
            (
                "una cobertura vacía",
                none(),
                Box::new(|v: &mut Value| v["coverage"] = json!([])),
                none(),
            ),
            "resumen FASTQ no tiene cobertura asociada",
        ),
        (
            (
                "un objetivo de cobertura negativo",
                effective("/technical_event/target/coverage_wgs", json!(-1.0)),
                none(),
                none(),
            ),
            "cobertura WGS técnica inválida",
        ),
        // Las reducciones.
        (
            (
                "una reducción que no es objeto",
                effective("/reduced_references/0", json!(5)),
                none(),
                none(),
            ),
            "una reducción no es objeto",
        ),
        (
            (
                "una reducción de otra referencia fuente",
                effective("/reduced_references/0/source_reference_id", json!("ref-x")),
                none(),
                none(),
            ),
            "la reducción efectiva ref-small no coincide con el manifiesto",
        ),
        (
            (
                "otra subsemilla de selección",
                effective(
                    "/reduced_references/0/selection_subseed",
                    json!("0".repeat(64)),
                ),
                none(),
                none(),
            ),
            "la reducción efectiva ref-small no coincide con el manifiesto",
        ),
        (
            (
                "un orden de selección con un texto",
                effective("/reduced_references/0/selection_order/0", json!("x")),
                none(),
                none(),
            ),
            "selection_order inválido",
        ),
        (
            (
                "un mapa de reducción que no está en el catálogo",
                effective(
                    "/reduced_references/0/map_artifact_id",
                    json!("reference-map-x"),
                ),
                none(),
                none(),
            ),
            "el enlace reference-map-x no existe",
        ),
        // El manifiesto, las secuencias y el VCF.
        (
            (
                "semilla que no es texto",
                effective("/seed/root", json!(5)),
                none(),
                none(),
            ),
            "root falta o está vacío",
        ),
        (
            (
                "autorizaciones que no son objeto",
                effective("/authorizations", json!(5)),
                none(),
                none(),
            ),
            "authorizations no coincide con el manifiesto",
        ),
        (
            (
                "sin autorizaciones en la efectiva",
                Box::new(|e: &mut Value| {
                    e.as_object_mut()
                        .expect("efectiva")
                        .remove("authorizations");
                }),
                none(),
                none(),
            ),
            "campos ausentes o desconocidos en effective-config.json",
        ),
        (
            (
                "sin enlace de VCF y sin artefactos de VCF",
                none(),
                none(),
                Box::new(|m: &mut Value| m["vcf_outputs"] = json!([])),
            ),
            "minimum/requested no coinciden con los roles normativos de F1",
        ),
        (
            (
                "VCF efectivo comprimido con enlace de texto",
                effective("/vcf_output/mode", json!("bgzf")),
                none(),
                none(),
            ),
            "referencia, modo o artefacto VCF del manifiesto no coincide con la configuración efectiva",
        ),
        (
            (
                "una secuencia que no es objeto",
                effective("/reference/sequences/0", json!(5)),
                none(),
                none(),
            ),
            "una secuencia no es objeto",
        ),
        (
            (
                "un FASTA que no está en el catálogo",
                effective("/reference/fasta_artifact_id", json!("reference-fasta-x")),
                none(),
                none(),
            ),
            "el enlace reference-fasta-x no existe",
        ),
        (
            (
                "un recuento de bases desconocidas que no es entero",
                effective("/reference/sequences/0/unknown_base_count", json!("x")),
                none(),
                none(),
            ),
            "unknown_base_count inválido",
        ),
        (
            (
                "una clase de secuencia desconocida",
                effective("/reference/sequences/0/class", json!("plasmid")),
                none(),
                none(),
            ),
            "clase de secuencia inválida",
        ),
        (
            (
                "sin matriz con poblaciones iniciales, coherente en los dos documentos",
                effective("/population_snvs", Value::Null),
                none(),
                Box::new(|m: &mut Value| {
                    m["population_snvs"] = Value::Null;
                    drop_artifact(m, "population-snvs-tsv");
                    drop_artifact(m, "population-snvs-metadata");
                }),
            ),
            "faltan estadísticas iniciales para cohortes iniciales",
        ),
        (
            (
                "allow_incompatible_genotypes que no es booleano en los dos documentos",
                effective("/authorizations/allow_incompatible_genotypes", json!("x")),
                none(),
                Box::new(|m: &mut Value| {
                    m["authorizations"]["allow_incompatible_genotypes"] = json!("x");
                }),
            ),
            "authorizations inválido",
        ),
        (
            (
                "un individuo técnico fuera de las poblaciones",
                Box::new(|e: &mut Value| {
                    e["initial_populations"][0]["individual_ids"] = json!(["I001", "I003"]);
                    e["initial_populations"][0]["individual_count"] = json!(2);
                }),
                none(),
                none(),
            ),
            "individual_id técnico fuera de las poblaciones iniciales",
        ),
        (
            (
                "un objetivo de reducción mayor que la referencia",
                effective("/reduced_references/0/target/bases", json!(1_u64 << 62)),
                none(),
                none(),
            ),
            "la suma del catálogo derivado no coincide con target.bases",
        ),
        (
            (
                "un artefacto VCF de más en el catálogo",
                none(),
                none(),
                Box::new(|m: &mut Value| {
                    let mut extra = artifact_of(m, "truth_vcf");
                    extra["id"] = json!("extra-vcf");
                    extra["path"] = json!("outputs/extra.vcf");
                    push(m, "/artifacts", extra);
                    push(m, "/output_categories/requested", json!("extra-vcf"));
                }),
            ),
            "minimum/requested no coinciden con los roles normativos de F1",
        ),
        (
            (
                "el artefacto VCF enlazado fuera del catálogo",
                none(),
                none(),
                Box::new(|m: &mut Value| drop_artifact(m, "truth-vcf")),
            ),
            "el enlace truth-vcf no existe",
        ),
    ];
    for ((label, effective_edit, event_edit, manifest_edit), expected) in attempts {
        assert_eq!(
            probe(&package, &*effective_edit, &*event_edit, &*manifest_edit),
            expected,
            "«{label}» se rechazó por otra razón"
        );
    }
    let manifest = package.json("manifest.json");
    let event_path = catalog_path(&manifest, "event-manifest-reads-main");
    assert_eq!(
        describe(probe_files(
            &package,
            &|_| {},
            &|_| {},
            &[(event_path, b"{no".to_vec())],
            true
        )),
        "GPS-FMT-006 JSON no legible: key must be a string at line 1 column 2",
        "un manifiesto de evento que no es JSON: lo rechaza el juez del evento"
    );

    // La reconstrucción de las estadísticas iniciales no puede fallar con una
    // matriz que ya pasó el par físico: un TSV alterado con su huella rehecha
    // lo rechaza antes la comprobación de la pareja de ficheros.
    let tsv = catalog_path(&manifest, "population-snvs-tsv");
    let mut changed = package.bytes(&tsv);
    let last = changed.len() - 2;
    changed[last] = if changed[last] == b'1' { b'2' } else { b'1' };
    assert_eq!(
        describe(probe_files(
            &package,
            &|_| {},
            &|_| {},
            &[(tsv, changed)],
            true
        )),
        "GPS-FMT-006 par SNV físico inválido: AuthoritativeFasta { detail: \"SHA-256 del TSV reutilizado no coincide con metadata\" }",
        "un TSV con la huella rehecha: lo rechaza el par físico de la matriz"
    );

    // La referencia simulada: la composición, la fracción de bases
    // desconocidas, las secuencias y los loci los ha validado ya el juez de la
    // efectiva (las comprobaciones de la regeneración repiten sus textos).
    let simulated = publish_simulated();
    let attempts: Vec<Attempt> = vec![
        (
            "una composición que no es número",
            Box::new(|e| e["reference"]["composition"]["A"] = json!("x")),
            "reference.composition.A inválida",
        ),
        (
            "una fracción de bases desconocidas negativa",
            Box::new(|e| e["reference"]["composition"]["N"] = json!(-0.1)),
            "reference.composition.N fuera de rango",
        ),
        (
            "una secuencia que no es objeto",
            Box::new(|e| e["reference"]["sequences"][0] = json!(5)),
            "una secuencia no es objeto",
        ),
        (
            "unos loci que no son una lista",
            Box::new(|e| e["repeat_loci"] = json!(5)),
            "tractos de la referencia simulada: repeat_loci no es una lista",
        ),
        (
            "una clase de secuencia desconocida",
            Box::new(|e| e["reference"]["sequences"][0]["class"] = json!("plasmid")),
            "clase de secuencia inválida",
        ),
    ];
    assert_eq!(
        probe_effective(&simulated, &|_| {}),
        "OK",
        "control positivo"
    );
    for (label, edit, expected) in attempts {
        assert_eq!(
            probe_effective(&simulated, &*edit),
            expected,
            "«{label}» se rechazó por otra razón"
        );
    }
}

/// **Restos de la configuración efectiva** que la tanda 7 no ejercitó: el
/// origen de la semilla, la tolerancia de normalización (sus dos extremos y
/// el valor siguiente), los metadatos de la normalización y los tres
/// ayudantes compartidos con que la efectiva lee una lista, un entero positivo
/// y un entero no negativo (una lista de individuos que no es lista, una
/// longitud a cero y una fase que no es un entero).
#[test]
fn the_effective_seed_the_normalization_and_the_shared_readers_are_judged() {
    let package = publish_rich();
    let none: &dyn Fn(&mut Value) = &|_| {};
    let on = |edit: fn(&mut Value)| probe(&package, &edit, none, none);
    assert_eq!(
        on(|e| e["seed"]["origin"] = json!("random")),
        "seed.origin inválido"
    );
    // La tolerancia admitida es [0, 1e-6]. Con la del TOML copiado la efectiva
    // diverge, pero antes la rechaza la validación propia del documento.
    assert_eq!(
        on(|e| e["normalization"]["tolerance"] = json!("x")),
        "normalization.tolerance inválida"
    );
    assert_eq!(
        on(|e| e["normalization"]["tolerance"] = json!(-1e-12)),
        "normalization.tolerance inválida"
    );
    assert_eq!(
        on(|e| e["normalization"]["tolerance"] = json!(f64::from_bits(1e-6_f64.to_bits() + 1))),
        "normalization.tolerance inválida"
    );
    // El lado aceptado de los dos extremos: pasa la validación del documento y
    // la rechaza después la comparación con el TOML copiado.
    for tolerance in [0.0, 1e-6] {
        let altered = probe(
            &package,
            &|e| e["normalization"]["tolerance"] = json!(tolerance),
            none,
            none,
        );
        assert_eq!(
            altered,
            "la tolerancia efectiva no coincide con la declaración TOML"
        );
    }
    for edit in [
        (|e: &mut Value| e["normalization"]["canonical_order"] = json!("population_id,contig"))
            as fn(&mut Value),
        |e| e["normalization"]["sum_algorithm"] = json!("kahan_sum_v1"),
    ] {
        assert_eq!(on(edit), "metadatos de normalización incompatibles");
    }
    assert_eq!(
        on(|e| e["initial_populations"][0]["individual_ids"] = json!(5)),
        "individual_ids no es lista"
    );
    assert_eq!(
        on(|e| e["reference"]["sequences"][1]["length"] = json!(0)),
        "length inválido"
    );
    assert_eq!(
        on(|e| e["direct_events"][0]["phase"] = json!("x")),
        "phase inválido"
    );
}

/// Un paquete de la **versión 6** (con sección de repeticiones): una
/// referencia suministrada con un microsatélite declarado a mano, cuatro
/// individuos diploides de una cohorte y VCF de texto.
fn publish_repeats() -> Package {
    publish_repeats_with("")
}

/// Lo mismo con un bloque TOML añadido al final (el evento técnico).
fn publish_repeats_with(tail: &str) -> Package {
    let bases = format!("{}{}{}", "G".repeat(100), "CA".repeat(50), "G".repeat(100));
    let toml = format!(
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
length = 300

[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 45
frequency = 0.7

[[reference.repeat_loci.alleles]]
units = 40
frequency = 0.3

[[initial_populations]]
population_id = "P1"
individual_count = 4
individual_ids = ["I001", "I002", "I003", "I004"]
ploidy = 2
genotype_model = "hwe"

[vcf_output]
mode = "text"
{tail}"#
    );
    publish(&[
        ("general.toml", &toml),
        ("reference.fa", &format!(">chr1\n{bases}\n")),
    ])
}

/// La matriz de SNV generada y la referencia del paquete rico, **sin
/// poblaciones iniciales, eventos ni VCF**.
fn publish_matrix_only() -> Package {
    let full = rich_toml("text", "");
    let toml = full
        .split("[[initial_populations]]")
        .next()
        .expect("TOML")
        .to_owned();
    publish(&[("general.toml", &toml), ("reference.fa", &rich_reference())])
}

/// Lo que dice `validate_published_tract_output_accounting` con estos
/// documentos y ficheros.
fn tract_judge(
    package: &Package,
    effective: &[u8],
    manifest: &[u8],
    truth: &Path,
    vcf: Option<&Path>,
) -> String {
    describe(
        genome_population_simulator::results::validate_published_tract_output_accounting(
            &package.effective_path(),
            effective,
            &package.manifest_path(),
            manifest,
            truth,
            vcf,
        ),
    )
}

/// El VCF con cada registro de repetición (los que llevan `RN=`) alterado por
/// `edit` sobre sus campos.
fn with_the_repeat_records_edited(vcf: &str, edit: impl Fn(&mut Vec<String>)) -> String {
    vcf.lines()
        .map(|line| {
            if !line.starts_with('#') && line.contains("RN=") {
                let mut fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
                edit(&mut fields);
                fields.join("\t")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

/// **El juez de la contabilidad de bases de tracto relee los bytes
/// publicados.** Sobre un paquete de la versión 6, con las dos puertas: la
/// pública (`validate_published_tract_output_accounting`, a la que se pasan
/// documentos y ficheros alterados) y el juez completo.
///
/// Cada caso altera un solo punto: los bytes del manifiesto o de la efectiva
/// que no son JSON, la contabilidad de tracto que falta o con una cifra que no
/// es un entero, la cota efectiva ausente, la verdad que falta, con una línea
/// que no es JSON o con un byte que no es UTF-8, el VCF que falta o con un
/// byte que no es UTF-8, y un registro de repetición que ya no contiene el
/// intervalo de su locus (con la huella rehecha, por el juez completo). Los
/// textos de los errores del sistema se calculan con el mismo error: lo que
/// se comprueba es el motivo y que el detalle lo lleva.
#[test]
fn the_tract_accounting_judge_rereads_the_published_bytes() {
    let package = publish_repeats();
    let manifest = package.json("manifest.json");
    let effective = package.json("effective-config.json");
    let (effective_bytes, manifest_bytes) = (bytes_of(&effective), bytes_of(&manifest));
    let truth = package
        .destination
        .join(catalog_path(&manifest, "truth-jsonl"));
    let vcf_relative = catalog_path(&manifest, "truth-vcf");
    let vcf = package.destination.join(&vcf_relative);
    let judge = |effective: &[u8], manifest: &[u8], truth: &Path, vcf: Option<&Path>| {
        tract_judge(&package, effective, manifest, truth, vcf)
    };
    let fmt = |detail: String| format!("GPS-FMT-006 {detail}");
    assert_eq!(
        judge(&effective_bytes, &manifest_bytes, &truth, Some(&vcf)),
        "OK",
        "control positivo"
    );

    // Los documentos.
    let not_json = serde_json::from_slice::<Value>(b"{no").expect_err("no es JSON");
    assert_eq!(
        judge(&effective_bytes, b"{no", &truth, Some(&vcf)),
        fmt(format!("JSON no legible: {not_json}")),
        "un manifiesto que no es JSON"
    );
    assert_eq!(
        judge(b"{no", &manifest_bytes, &truth, Some(&vcf)),
        fmt(format!("JSON no legible: {not_json}")),
        "una efectiva que no es JSON"
    );
    let mut altered = manifest.clone();
    altered["repeat_loci"]
        .as_object_mut()
        .expect("sección")
        .remove("tract_output");
    assert_eq!(
        judge(&effective_bytes, &bytes_of(&altered), &truth, Some(&vcf)),
        fmt("falta la contabilidad de tracto".to_owned()),
        "sin contabilidad de tracto"
    );
    let mut altered = manifest.clone();
    altered["repeat_loci"]["tract_output"]["total_tract_bases_bound"] = json!("x");
    assert_eq!(
        judge(&effective_bytes, &bytes_of(&altered), &truth, Some(&vcf)),
        fmt("una cifra de tracto no es un entero".to_owned()),
        "una cifra que no es un entero"
    );
    let mut altered = effective.clone();
    altered["repeat_output_budget"]
        .as_object_mut()
        .expect("presupuesto")
        .remove("total_tract_bases_bound");
    assert_eq!(
        judge(&bytes_of(&altered), &manifest_bytes, &truth, Some(&vcf)),
        fmt("la configuración efectiva no publica la cota de bases de tracto".to_owned()),
        "sin la cota efectiva"
    );

    // La verdad.
    let missing = package.destination.join("no-such-truth.jsonl");
    let open_error = fs::File::open(&missing).expect_err("no existe");
    assert_eq!(
        judge(&effective_bytes, &manifest_bytes, &missing, Some(&vcf)),
        fmt(format!("verdad no legible: {open_error}")),
        "una verdad que falta"
    );
    let bad_truth = package.destination.join("bad-truth.jsonl");
    fs::write(&bad_truth, b"{not json\n").expect("verdad alterada");
    let line_error = serde_json::from_slice::<Value>(b"{not json").expect_err("no es JSON");
    assert_eq!(
        judge(&effective_bytes, &manifest_bytes, &bad_truth, Some(&vcf)),
        fmt(format!("línea de verdad ilegible: {line_error}")),
        "una línea de la verdad que no es JSON"
    );
    fs::write(&bad_truth, b"\xff\xfe\n").expect("verdad alterada");
    let utf8_error = std::io::BufRead::lines(&b"\xff\xfe\n"[..])
        .next()
        .expect("una línea")
        .expect_err("no es UTF-8");
    assert_eq!(
        judge(&effective_bytes, &manifest_bytes, &bad_truth, Some(&vcf)),
        fmt(format!("verdad no legible: {utf8_error}")),
        "una verdad con un byte que no es UTF-8"
    );

    // El VCF.
    let missing_vcf = package.destination.join("no-such.vcf");
    let open_error = fs::File::open(&missing_vcf).expect_err("no existe");
    assert_eq!(
        judge(
            &effective_bytes,
            &manifest_bytes,
            &truth,
            Some(&missing_vcf)
        ),
        fmt(format!("VCF no legible: {open_error}")),
        "un VCF que falta"
    );
    let original = String::from_utf8(package.bytes(&vcf_relative)).expect("VCF");
    let bad_vcf = package.destination.join("bad.vcf");
    fs::write(&bad_vcf, b"##fileformat=VCFv4.2\n\xff\n").expect("VCF alterado");
    assert_eq!(
        judge(&effective_bytes, &manifest_bytes, &truth, Some(&bad_vcf)),
        fmt("VCF no es UTF-8".to_owned()),
        "un VCF con un byte que no es UTF-8"
    );
    let interval = "un registro de repetición del VCF no contiene el intervalo de ningún locus";
    // Cada condición del intervalo por separado: el inicio (POS una base más
    // allá), el final (un REF más corto, que acaba antes del fin del locus) y
    // la secuencia (otro CHROM).
    let shifted = with_the_repeat_records_edited(&original, |fields| {
        let position: u64 = fields[1].parse().expect("POS");
        fields[1] = (position + 1).to_string();
    });
    let shorter = with_the_repeat_records_edited(&original, |fields| {
        fields[3].pop();
    });
    let other_sequence = with_the_repeat_records_edited(&original, |fields| {
        fields[0] = "chr2".to_owned();
    });
    for (label, altered) in [
        ("un registro de repetición desplazado", &shifted),
        ("un registro de repetición con un REF más corto", &shorter),
        (
            "un registro de repetición de otra secuencia",
            &other_sequence,
        ),
    ] {
        assert_ne!(altered, &original, "la alteración alteró");
        fs::write(&bad_vcf, altered).expect("VCF alterado");
        assert_eq!(
            judge(&effective_bytes, &manifest_bytes, &truth, Some(&bad_vcf)),
            fmt(interval.to_owned()),
            "{label}"
        );
    }

    // El juez completo, con la huella del VCF rehecha en el manifiesto.
    let none: &dyn Fn(&mut Value) = &|_| {};
    assert_eq!(
        describe(probe_files(
            &package,
            none,
            none,
            &[(vcf_relative.clone(), shifted.into_bytes())],
            true
        )),
        fmt(interval.to_owned()),
        "juez completo: un registro de repetición desplazado"
    );
    let mut invalid = original.into_bytes();
    let header_end = invalid
        .windows(2)
        .position(|pair| pair == b"\nc")
        .expect("primera línea de datos");
    invalid.insert(header_end + 2, 0xff);
    assert_eq!(
        describe(probe_files(
            &package,
            none,
            none,
            &[(vcf_relative, invalid)],
            true
        )),
        fmt("VCF no es UTF-8".to_owned()),
        "juez completo: un byte que no es UTF-8 en una línea de datos"
    );

    // Con un paquete de la versión 4 el VCF se lee **solo la cabecera** (para
    // el `execution_id`): un byte que no es UTF-8 en ella se rechaza igual.
    let rich = publish_rich();
    let rich_vcf = catalog_path(&rich.json("manifest.json"), "truth-vcf");
    let mut header = rich.bytes(&rich_vcf);
    let first = header
        .windows(2)
        .position(|pair| pair == b"##")
        .expect("cabecera");
    header.insert(first + 5, 0xff);
    assert_eq!(
        describe(probe_files(&rich, none, none, &[(rich_vcf, header)], true)),
        fmt("VCF no es UTF-8".to_owned()),
        "versión 4: un byte que no es UTF-8 en la cabecera"
    );
}

/// **Los artefactos de la matriz de SNV conservan su identidad publicada**,
/// cada uno con el tamaño y la huella **por separado**: el TSV y los
/// metadatos. Un byte cambiado sin rehacer la huella (mismo tamaño) y un
/// tamaño declarado distinto con la huella verdadera.
#[test]
fn the_snv_matrix_artifacts_keep_their_published_identity() {
    let package = publish_rich();
    let none: &dyn Fn(&mut Value) = &|_| {};
    let manifest = package.json("manifest.json");
    let metadata_id = manifest["population_snvs"]["metadata_artifact_id"]
        .as_str()
        .expect("metadatos")
        .to_owned();
    for id in ["population-snvs-tsv", metadata_id.as_str()] {
        let path = catalog_path(&manifest, id);
        let expected = format!("GPS-FMT-006 artefacto SNV físico {id} no coincide con su SHA-256");
        let mut changed = package.bytes(&path);
        let position = changed
            .iter()
            .position(|byte| *byte == b' ' || byte.is_ascii_digit())
            .expect("un byte que cambiar");
        changed[position] = if changed[position] == b'0' {
            b'1'
        } else {
            b'0'
        };
        assert_ne!(changed, package.bytes(&path));
        assert_eq!(
            describe(probe_files(
                &package,
                none,
                none,
                &[(path.clone(), changed)],
                false
            )),
            expected,
            "{id} con otra huella"
        );
        let size = package.bytes(&path).len() as u64;
        assert_eq!(
            describe(probe_files(
                &package,
                none,
                &|manifest| set_identity(manifest, &path, None, Some(size + 1)),
                &[],
                false
            )),
            expected,
            "{id} con otro tamaño declarado"
        );
    }
}

/// **Unas estadísticas iniciales sin poblaciones iniciales se rechazan.** Un
/// paquete con la matriz y sin poblaciones (control positivo: se acepta); las
/// estadísticas de otro paquete, puestas en la efectiva y en el manifiesto a
/// la vez, se rechazan en el juez físico de las estadísticas.
#[test]
fn statistics_without_initial_populations_are_rejected() {
    let stats = publish_rich().json("effective-config.json")["population_snvs"]
        ["initial_population_statistics"]
        .clone();
    assert!(stats.is_object(), "el fixture lleva estadísticas");
    let package = publish_matrix_only();
    let none: &dyn Fn(&mut Value) = &|_| {};
    assert_eq!(
        describe(probe_files(&package, none, none, &[], false)),
        "OK",
        "control positivo"
    );
    assert_eq!(
        describe(probe_files(
            &package,
            &|effective| effective["population_snvs"]["initial_population_statistics"] =
                stats.clone(),
            &|manifest| manifest["population_snvs"]["initial_population_statistics"] =
                stats.clone(),
            &[],
            false
        )),
        "GPS-FMT-006 hay estadísticas iniciales sin poblaciones iniciales"
    );
}

/// **El manifiesto general lleva el perfil técnico y el catálogo de azar de la
/// configuración efectiva** (DT-056). En el paquete rico (sin PCR: la efectiva
/// dice `f4-technical-v3` y el perfil sin PCR) se cambia en el manifiesto
/// general, de uno en uno, el catálogo de azar (v4 y v5) y el perfil entero
/// (el de PCR, con su huella). Hasta el 12/10/2026 el juez completo
/// **aceptaba** los tres. El identificador y la huella del perfil **no se
/// prueban por separado**: el lector del manifiesto exige el par entero de
/// uno de los dos perfiles, así que cambiar solo uno lo rechaza él antes
/// («inventario técnico F4 incompatible»); las dos mitades de la comparación
/// del juez completo se implican.
///
/// **Desde el 12/10/2026 (DT-041) el lector ata además el catálogo al perfil y
/// a la versión del manifiesto**, de modo que los tres cambios sueltos los
/// rechaza él antes de llegar a la comparación con la efectiva: en un
/// manifiesto v4 sin PCR solo cabe el v3. Para llegar a la comparación del
/// perfil se cambian a la vez el perfil y el catálogo (PCR con v5, coherentes
/// entre sí). La comparación del catálogo, con perfil y catálogo coherentes y
/// distintos de la efectiva, la ejercita
/// `the_judges_and_the_schemas_reject_a_catalog_that_contradicts_the_run`
/// (`tests/technical_catalog_contract.rs`) sobre la ejecución sin PCR con loci
/// (enlace en v3 frente a la efectiva en v4).
#[test]
fn the_general_manifest_carries_the_profile_and_the_catalog_of_the_effective_config() {
    use genome_population_simulator::technical_profile::{
        GPS_HISEQX_PCR_PE150_SHA256, GPS_HISEQX_PCR_PE150_V1,
    };
    let package = publish_rich();
    let none: &dyn Fn(&mut Value) = &|_| {};
    assert_eq!(probe(&package, none, none, none), "OK", "control positivo");
    let catalog = "el catálogo de azar del inventario técnico no corresponde a su perfil o a la versión del manifiesto";
    let profile = "el perfil técnico del manifiesto general no es el de la configuración efectiva";
    for (label, edit, expected) in [
        (
            "el catálogo v4",
            Box::new(|m: &mut Value| {
                m["technical_events"][0]["rng_catalog"] = json!("f4-technical-v4")
            }) as Mutate,
            catalog,
        ),
        (
            "el catálogo v5",
            Box::new(|m: &mut Value| {
                m["technical_events"][0]["rng_catalog"] = json!("f4-technical-v5")
            }),
            catalog,
        ),
        (
            "el perfil de PCR con su huella",
            Box::new(|m: &mut Value| {
                m["technical_events"][0]["profile"] = json!({
                    "id": GPS_HISEQX_PCR_PE150_V1,
                    "sha256": GPS_HISEQX_PCR_PE150_SHA256,
                });
            }),
            catalog,
        ),
        (
            "el perfil de PCR con su huella y el catálogo v5",
            Box::new(|m: &mut Value| {
                m["technical_events"][0]["profile"] = json!({
                    "id": GPS_HISEQX_PCR_PE150_V1,
                    "sha256": GPS_HISEQX_PCR_PE150_SHA256,
                });
                m["technical_events"][0]["rng_catalog"] = json!("f4-technical-v5");
            }),
            profile,
        ),
    ] {
        assert_eq!(probe(&package, none, none, &*edit), expected, "{label}");
    }
}

/// **Un manifiesto de versión 6 acepta los catálogos de azar que su esquema
/// admite.** El esquema v6 admite `f4-technical-v3` y `f4-technical-v4`: el
/// lector los acepta en un paquete de esa versión con un evento técnico (el
/// de repeticiones del fixture más lecturas), y sigue rechazando un catálogo
/// desconocido. El juez completo acepta el paquete intacto.
#[test]
fn a_version_6_manifest_accepts_the_random_catalogs_its_schema_admits() {
    let technical = r#"
[technical_event]
event_id = "reads-main"
individual_ids = ["I001", "I002"]
target_coverage_wgs = 1
"#;
    let toml_tail = format!(
        "{technical}\n[[technical_event.libraries]]\nindividual_id = \"I001\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n\n[[technical_event.libraries]]\nindividual_id = \"I002\"\ni7 = \"TTGGCCAA\"\ni5 = \"GGAATTCC\"\n"
    );
    let package = publish_repeats_with(&toml_tail);
    let manifest = package.json("manifest.json");
    assert_eq!(manifest["schema_version"], "6");
    let cases = [
        set(
            "el catálogo de azar v4",
            "",
            "/technical_events/0/rng_catalog",
            json!("f4-technical-v4"),
        ),
        set(
            "un catálogo de azar desconocido",
            "inventario técnico F4 incompatible",
            "/technical_events/0/rng_catalog",
            json!("f4-technical-v2"),
        ),
    ];
    assert_parse_cases(&package, &cases);
    let effective = package.json("effective-config.json");
    validate_published_resources(
        &package.destination,
        &package.effective_path(),
        &bytes_of(&effective),
        &package.manifest_path(),
        &bytes_of(&manifest),
        None,
    )
    .expect("el juez completo acepta el paquete de la versión 6 con evento");
}
