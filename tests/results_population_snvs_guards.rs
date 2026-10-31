//! Cada rechazo del juez de los recursos publicados en lo que toca a las
//! estadísticas iniciales de las poblaciones y a la auditoría de la
//! normalización, con su motivo exacto.
//!
//! DT-001, tanda 8 (tercera parte de `src/results.rs`). Cubre
//! `validate_initial_population_statistics` y sus ayudantes (procedencia,
//! población, modelo genotípico y locus), `validate_serialized_normalization_audit`
//! y la parte del juez físico que lee la matriz de SNV
//! (`validate_physical_population_snvs`).
//!
//! # Cómo se parte de un paquete real
//!
//! Cada caso parte de un paquete **publicado de verdad** por el programa
//! (referencia simulada de dos secuencias, matriz de SNV generada, tres
//! poblaciones de dos ploidías —2, 2 y 3— con los dos modelos genotípicos) y
//! altera **un solo punto** del manifiesto general. El juez lee el manifiesto
//! antes que nada, de modo que el rechazo sale del lector del manifiesto y su
//! texto es el motivo exacto; la efectiva queda intacta, y por eso, si una
//! guarda dejara de rechazar, el caso no pasaría en silencio: caería después en
//! la comparación efectiva/manifiesto, con **otro** texto.
//!
//! Las alteraciones que **deben aceptarse** (lado aceptado de cada límite, de
//! cada lista cerrada y de las ramas que el paquete publicado no produce, como
//! un locus de tres alelos o una publicación sin copias observadas) se aplican
//! a los dos documentos por igual, de modo que el juez de la pareja las vea
//! coherentes. Cada tabla lleva su control positivo y comprueba que la
//! alteración alteró.
//!
//! # Locus y poblaciones sintéticos
//!
//! La matriz generada solo trae loci de dos alelos y sin copias ausentes. Para
//! las ramas que el generador no produce, la base de algunas tablas añade a los
//! dos documentos un locus (o una población) escrito a mano: tres alelos, o
//! ninguna copia observada, o una ploidía distinta. Es información sintética en
//! un documento publicado de verdad, y lo dice el nombre de cada función que la
//! construye.
//!
//! **Cumplen el contrato de la arquitectura**, no solo la lista de claves: los
//! genotipos realizados cuentan solo observaciones completas, así que por alelo
//! la suma de dosis por recuento no supera el recuento realizado y las
//! observaciones completas por la ploidía no superan las copias observadas. El
//! juez de la pareja no comprueba esas dos desigualdades (DT-055; las caza
//! después el juez completo), de modo que nada del juez obligaría a cumplirlas:
//! lo hace `assert_architecture_inequalities` sobre la base y sobre cada caso
//! «aceptado», para no consagrar un documento que el contrato no admite.
//!
//! # Lo que este fichero no es
//!
//! Las guardas que ningún paquete publicado y manipulado puede alcanzar están
//! anotadas **junto a la guarda** en `src/results.rs` con el porqué real; la
//! prueba que solo «ejecuta el intento» sobre ellas lo dice en su docstring y
//! **no discrimina** esas guardas.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::derive_execution_id;
use genome_population_simulator::integrity::{configuration_digest, identify_file};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::results::{
    validate_effective_config_against_manifest_bytes, validate_published_resources,
};
use serde_json::{Value, json};

#[path = "support/detection_section.rs"]
mod detection_section;

static NONCE: AtomicU64 = AtomicU64::new(0);

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
        serde_json::from_slice(&fs::read(self.destination.join(relative)).expect("artefacto"))
            .expect("JSON publicado")
    }
}

fn publish(toml: &str) -> Package {
    publish_files(&[("general.toml", toml)])
}

fn publish_files(files: &[(&str, &str)]) -> Package {
    let root = std::env::temp_dir().join(format!(
        "gps-results-population-snvs-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    // La guarda que borra el directorio existe **antes** de crearlo: si la
    // publicación falla, el directorio no queda.
    let package = Package {
        destination: root.join("published"),
        root,
    };
    fs::create_dir_all(&package.root).expect("directorio");
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

/// Tres poblaciones —ploidía 2 con `hwe`, ploidía 2 con `f_diploid` y ploidía
/// 3 con `hwe`—, una matriz generada de tres loci en dos secuencias y una
/// referencia simulada.
const STATISTICS_TOML: &str = r#"config_kind = "general"
schema_version = "4"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 3
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"

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

[[population_snvs.generate.populations]]
population_id = "P3"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "simulated"
id = "ref-main"

[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 400

[[reference.sequences]]
id = "scaffold-1"
class = "scaffold"
length = 300

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

[[initial_populations]]
population_id = "P3"
individual_count = 2
individual_ids = ["J001", "J002"]
ploidy = 3
genotype_model = "hwe"
"#;

// ---------------------------------------------------------------------------
// Maquinaria de casos
// ---------------------------------------------------------------------------

fn detail(result: Result<(), Diagnostic>) -> String {
    result
        .expect_err("el juez tenía que rechazar")
        .context()
        .get("detail")
        .cloned()
        .unwrap_or_default()
}

/// El juez de la pareja efectiva/manifiesto: lee el manifiesto antes que la
/// efectiva y ata la una a la otra.
fn judge(package: &Package, effective: &Value, manifest: &Value) -> Result<(), Diagnostic> {
    validate_effective_config_against_manifest_bytes(
        &package.destination.join("effective-config.json"),
        &serde_json::to_vec_pretty(effective).expect("JSON"),
        &package.destination.join("manifest.json"),
        &serde_json::to_vec_pretty(manifest).expect("JSON"),
    )
}

/// Qué se altera: solo el manifiesto (un rechazo) o los dos documentos por
/// igual (un caso que debe aceptarse).
struct Case {
    label: String,
    expected: String,
    mutate: Box<dyn Fn(&mut Value)>,
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

/// Un rechazo: cambia un valor del manifiesto.
fn set(
    label: impl Into<String>,
    expected: impl Into<String>,
    pointer: impl Into<String>,
    value: Value,
) -> Case {
    let pointer = pointer.into();
    case(label, expected, move |document| {
        *document
            .pointer_mut(&pointer)
            .unwrap_or_else(|| panic!("el fixture no tiene {pointer}")) = value.clone();
    })
}

/// Un rechazo: quita una clave.
fn remove(
    label: impl Into<String>,
    expected: impl Into<String>,
    pointer: impl Into<String>,
) -> Case {
    let pointer = pointer.into();
    case(label, expected, move |document| {
        remove_key(document, &pointer)
    })
}

fn remove_key(document: &mut Value, pointer: &str) {
    let (parent, key) = pointer.rsplit_once('/').expect("puntero");
    match document
        .pointer_mut(parent)
        .unwrap_or_else(|| panic!("el fixture no tiene {parent}"))
    {
        Value::Object(map) => {
            assert!(map.remove(key).is_some(), "no había {pointer}");
        }
        Value::Array(items) => {
            items.remove(key.parse::<usize>().expect("índice"));
        }
        _ => panic!("{pointer} no cuelga de un contenedor"),
    }
}

/// Un rechazo: añade una clave que el contrato no conoce.
fn add_key(
    label: impl Into<String>,
    expected: impl Into<String>,
    pointer: impl Into<String>,
) -> Case {
    let pointer = pointer.into();
    case(label, expected, move |document| {
        document
            .pointer_mut(&pointer)
            .and_then(Value::as_object_mut)
            .unwrap_or_else(|| panic!("{pointer} no es un objeto"))
            .insert("unexpected".to_owned(), json!(1));
    })
}

/// Un rechazo: renombra una clave (el objeto conserva su número de claves).
fn rename(
    label: impl Into<String>,
    expected: impl Into<String>,
    pointer: impl Into<String>,
    key: &'static str,
) -> Case {
    let pointer = pointer.into();
    case(label, expected, move |document| {
        let object = document
            .pointer_mut(&pointer)
            .and_then(Value::as_object_mut)
            .unwrap_or_else(|| panic!("{pointer} no es un objeto"));
        let value = object
            .remove(key)
            .unwrap_or_else(|| panic!("no había {key}"));
        object.insert(format!("{key}_renamed"), value);
    })
}

/// Un caso que debe aceptarse: lo aplica a la efectiva **y** al manifiesto.
fn ok(label: impl Into<String>, pointer: impl Into<String>, value: Value) -> Case {
    set(label, "", pointer, value)
}

fn ok_with(label: impl Into<String>, mutate: impl Fn(&mut Value) + 'static) -> Case {
    case(label, "", mutate)
}

type Base = fn(&mut Value);

fn no_base(_: &mut Value) {}

/// Las dos desigualdades de la arquitectura que el juez de la pareja **no**
/// comprueba (el juez completo las caza después): los genotipos realizados
/// cuentan solo observaciones completas, así que por alelo la suma de dosis por
/// recuento no supera su recuento realizado, y las observaciones completas por
/// la ploidía no superan las copias observadas. Todo documento que una prueba
/// de este fichero da por **aceptado** (la base y cada caso «aceptado») las
/// cumple: no se consagra un documento que el contrato no admite.
fn assert_architecture_inequalities(document: &Value, label: &str) {
    let Some(populations) = document
        .pointer(&stats("/populations"))
        .and_then(Value::as_array)
    else {
        return;
    };
    for population in populations {
        let ploidy = population["ploidy"].as_u64().expect("ploidía");
        for locus in population["loci"].as_array().expect("loci") {
            let complete = locus["complete_genotype_observations"]
                .as_u64()
                .expect("completas");
            let observed = locus["observed_copies"].as_u64().expect("observadas");
            assert!(
                complete * ploidy <= observed,
                "«{label}»: {complete} observaciones completas por ploidía {ploidy} superan las {observed} copias observadas"
            );
            let alleles = locus["alleles"].as_array().expect("alelos");
            let mut used = vec![0_u64; alleles.len()];
            for genotype in locus["genotypes"].as_array().expect("genotipos") {
                let realized = genotype["realized_count"].as_u64().expect("recuento");
                for (index, dosage) in genotype["allele_dosages"]
                    .as_array()
                    .expect("dosis")
                    .iter()
                    .enumerate()
                {
                    used[index] += dosage.as_u64().expect("dosis") * realized;
                }
            }
            for (index, allele) in alleles.iter().enumerate() {
                let count = allele["realized_count"].as_u64().expect("recuento alélico");
                assert!(
                    used[index] <= count,
                    "«{label}»: los genotipos usan {} copias del alelo {index} y solo hay {count}",
                    used[index]
                );
            }
        }
    }
}

/// Los rechazos: parten de `base` aplicada a los dos documentos y alteran solo
/// el manifiesto.
fn assert_rejected(package: &Package, base: Base, cases: &[Case]) {
    let mut effective = package.json("effective-config.json");
    let mut manifest = package.json("manifest.json");
    base(&mut effective);
    base(&mut manifest);
    judge(package, &effective, &manifest)
        .expect("control positivo: el paquete (con su base) se acepta");
    assert_architecture_inequalities(&manifest, "base");
    for case in cases {
        let mut altered = manifest.clone();
        (case.mutate)(&mut altered);
        assert_ne!(
            serde_json::to_string(&altered).expect("JSON"),
            serde_json::to_string(&manifest).expect("JSON"),
            "«{}»: la alteración no alteró nada",
            case.label
        );
        assert_eq!(
            detail(judge(package, &effective, &altered)),
            case.expected,
            "«{}» se rechazó por otra razón",
            case.label
        );
    }
}

/// Lo aceptado: la alteración se aplica a los dos documentos por igual.
fn assert_accepted(package: &Package, base: Base, cases: &[Case]) {
    let mut effective = package.json("effective-config.json");
    let mut manifest = package.json("manifest.json");
    base(&mut effective);
    base(&mut manifest);
    for case in cases {
        let (mut altered_effective, mut altered_manifest) = (effective.clone(), manifest.clone());
        (case.mutate)(&mut altered_effective);
        (case.mutate)(&mut altered_manifest);
        assert_ne!(
            serde_json::to_string(&altered_manifest).expect("JSON"),
            serde_json::to_string(&manifest).expect("JSON"),
            "«{}»: la alteración no alteró nada",
            case.label
        );
        assert_architecture_inequalities(&altered_manifest, &case.label);
        judge(package, &altered_effective, &altered_manifest)
            .unwrap_or_else(|error| panic!("«{}» se tenía que aceptar: {error:?}", case.label));
    }
}

fn statistics_package() -> Package {
    publish(STATISTICS_TOML)
}

// Raíces de los punteros.
const STATS: &str = "/population_snvs/initial_population_statistics";
const AUDIT: &str = "/population_snvs/normalization_audit";

fn stats(suffix: &str) -> String {
    format!("{STATS}{suffix}")
}

fn audit(suffix: &str) -> String {
    format!("{AUDIT}{suffix}")
}

fn append(document: &mut Value, pointer: &str, item: Value) {
    document
        .pointer_mut(pointer)
        .and_then(Value::as_array_mut)
        .unwrap_or_else(|| panic!("{pointer} no es una lista"))
        .push(item);
}

fn replace_at(document: &mut Value, pointer: &str, value: Value) {
    *document
        .pointer_mut(pointer)
        .unwrap_or_else(|| panic!("el fixture no tiene {pointer}")) = value;
}

fn pointer_value(document: &Value, pointer: &str) -> Value {
    document
        .pointer(pointer)
        .unwrap_or_else(|| panic!("el fixture no tiene {pointer}"))
        .clone()
}

// ---------------------------------------------------------------------------
// Las estadísticas iniciales: la raíz y las poblaciones
// ---------------------------------------------------------------------------

const STATS_INVALID: &str = "estadísticas iniciales inválidas";

#[test]
fn the_statistics_root_is_judged() {
    let package = statistics_package();
    let cases = [
        set(
            "estadísticas que no son un objeto",
            STATS_INVALID,
            STATS,
            json!(5),
        ),
        add_key("una clave de más en la raíz", STATS_INVALID, STATS),
        remove(
            "una clave de menos en la raíz",
            STATS_INVALID,
            stats("/version"),
        ),
        rename(
            "una clave desconocida con el mismo número de claves",
            STATS_INVALID,
            STATS,
            "version",
        ),
        // La procedencia ausente o renombrada tiene su propio motivo si la
        // lista cerrada de claves no la parase antes: es lo que distingue las
        // dos mitades de esa guarda de las comprobaciones que vienen después
        // (quitar `version` o renombrarla da el mismo texto en las dos).
        remove(
            "sin procedencia: lo para la lista cerrada de claves",
            STATS_INVALID,
            stats("/provenance"),
        ),
        rename(
            "procedencia con otro nombre: lo para la lista cerrada de claves",
            STATS_INVALID,
            STATS,
            "provenance",
        ),
        set(
            "otro tipo de estadística",
            STATS_INVALID,
            stats("/kind"),
            json!("other"),
        ),
        set(
            "otra versión de estadística",
            STATS_INVALID,
            stats("/version"),
            json!("2"),
        ),
        set(
            "versión numérica",
            STATS_INVALID,
            stats("/version"),
            json!(1),
        ),
        set(
            "procedencia que no es un objeto",
            "procedencia de estadísticas iniciales inválida",
            stats("/provenance"),
            json!([]),
        ),
        set(
            "poblaciones que no son una lista",
            STATS_INVALID,
            stats("/populations"),
            json!({}),
        ),
        set(
            "ninguna población",
            "estadísticas iniciales sin poblaciones",
            stats("/populations"),
            json!([]),
        ),
        set(
            "dos poblaciones con el mismo identificador",
            "población estadística repetida",
            stats("/populations/1/population_id"),
            json!("P1"),
        ),
    ];
    assert_rejected(&package, no_base, &cases);
}

const POPULATION_INVALID: &str = "población estadística inválida";
const MODEL_INVALID: &str = "modelo genotípico estadístico inválido";

/// Una población sintética: ploidía, número de individuos y cuántos individuos
/// hay de cada genotipo (`realized[k]` con `k` copias alternas, `k = 0..=ploidy`,
/// que suman los individuos). Los recuentos de los alelos y la frecuencia
/// realizada **salen de los genotipos** (cada individuo completo aporta
/// `ploidía` copias), de modo que el documento cumple el contrato de la
/// arquitectura: los genotipos realizados solo cuentan observaciones completas,
/// y aquí todas las copias lo son. Un locus de dos alelos.
fn synthetic_population(id: &str, ploidy: u64, realized: &[u64], model: Value) -> Value {
    assert_eq!(realized.len() as u64, ploidy + 1);
    let individuals: u64 = realized.iter().sum();
    let total = individuals * ploidy;
    let alternate_count: u64 = realized.iter().zip(0..).map(|(count, k)| count * k).sum();
    let reference_count = total - alternate_count;
    let genotypes: Vec<Value> = (0..=ploidy)
        .map(|alternate| {
            json!({
                "allele_dosages": [ploidy - alternate, alternate],
                "expected_probability": 1.0 / (ploidy + 1) as f64,
                "expected_count": individuals as f64 / (ploidy + 1) as f64,
                "realized_count": realized[alternate as usize],
            })
        })
        .collect();
    json!({
        "population_id": id,
        "genotype_model": model,
        "individual_count": individuals,
        "ploidy": ploidy,
        "loci": [{
            "contig": "chr1",
            "position_1based": 8,
            "ref": "A",
            "total_copies": total,
            "observed_copies": total,
            "missing_by_deletion": 0,
            "alleles": [
                {"allele": "A", "kind": "ref", "requested_frequency": 0.5,
                 "realized_count": reference_count,
                 "realized_frequency": reference_count as f64 / total as f64},
                {"allele": "G", "kind": "alt", "variant_id": "VAR_SYNTHETIC",
                 "source_id": "supplied", "requested_frequency": 0.5,
                 "realized_count": alternate_count,
                 "realized_frequency": alternate_count as f64 / total as f64},
            ],
            "complete_genotype_observations": individuals,
            "genotypes": genotypes,
        }],
    })
}

fn with_synthetic_populations(document: &mut Value) {
    // Ploidía 1 con tres individuos (uno de referencia y dos alternos) y
    // ploidía 10 con dos (uno homocigoto de referencia y otro de alterno).
    let mut ploidy_ten = vec![0; 11];
    (ploidy_ten[0], ploidy_ten[10]) = (1, 1);
    for (id, ploidy, realized) in [("P8", 1, vec![1, 2]), ("P9", 10, ploidy_ten)] {
        append(
            document,
            &stats("/populations"),
            synthetic_population(id, ploidy, &realized, json!({"kind": "hwe"})),
        );
    }
}

#[test]
fn the_statistics_populations_are_judged() {
    let package = statistics_package();
    let populations = stats("/populations");
    let p0 = |suffix: &str| format!("{populations}/0{suffix}");
    let cases = [
        set(
            "una población que no es un objeto",
            POPULATION_INVALID,
            p0(""),
            json!(5),
        ),
        add_key(
            "una clave de más en la población",
            POPULATION_INVALID,
            p0(""),
        ),
        remove(
            "una clave de menos en la población",
            POPULATION_INVALID,
            p0("/ploidy"),
        ),
        rename(
            "una clave desconocida en la población",
            POPULATION_INVALID,
            p0(""),
            "ploidy",
        ),
        remove(
            "sin modelo genotípico: lo para la lista cerrada de claves",
            POPULATION_INVALID,
            p0("/genotype_model"),
        ),
        rename(
            "los loci con otro nombre: lo para la lista cerrada de claves",
            POPULATION_INVALID,
            p0(""),
            "loci",
        ),
        set(
            "identificador de población vacío",
            POPULATION_INVALID,
            p0("/population_id"),
            json!(""),
        ),
        set(
            "identificador de población numérico",
            POPULATION_INVALID,
            p0("/population_id"),
            json!(5),
        ),
        set(
            "ningún individuo",
            POPULATION_INVALID,
            p0("/individual_count"),
            json!(0),
        ),
        set(
            "individuos que no son un número",
            POPULATION_INVALID,
            p0("/individual_count"),
            json!("3"),
        ),
        set(
            "individuos negativos",
            POPULATION_INVALID,
            p0("/individual_count"),
            json!(-3),
        ),
        set("ploidía 0", POPULATION_INVALID, p0("/ploidy"), json!(0)),
        set("ploidía 11", POPULATION_INVALID, p0("/ploidy"), json!(11)),
        set(
            "ploidía que no es un número",
            POPULATION_INVALID,
            p0("/ploidy"),
            json!("2"),
        ),
        set(
            "modelo que no es un objeto",
            MODEL_INVALID,
            p0("/genotype_model"),
            json!("hwe"),
        ),
        set(
            "modelo hwe con una clave de más",
            MODEL_INVALID,
            p0("/genotype_model"),
            json!({"kind": "hwe", "inbreeding_coefficient": 0.5}),
        ),
        set(
            "modelo de otra clase",
            MODEL_INVALID,
            p0("/genotype_model"),
            json!({"kind": "other"}),
        ),
        set(
            "modelo sin clase",
            MODEL_INVALID,
            p0("/genotype_model"),
            json!({"inbreeding_coefficient": 0.5}),
        ),
        set(
            "f_diploid con ploidía 3",
            MODEL_INVALID,
            format!("{populations}/2/genotype_model"),
            json!({"kind": "f_diploid", "inbreeding_coefficient": 0.5}),
        ),
        case("f_diploid con ploidía 1", MODEL_INVALID, {
            let populations = populations.clone();
            move |document| {
                replace_at(document, &format!("{populations}/0/ploidy"), json!(1));
                replace_at(
                    document,
                    &format!("{populations}/0/genotype_model"),
                    json!({"kind": "f_diploid", "inbreeding_coefficient": 0.5}),
                );
            }
        }),
        remove(
            "f_diploid sin coeficiente",
            MODEL_INVALID,
            format!("{populations}/1/genotype_model/inbreeding_coefficient"),
        ),
        add_key(
            "f_diploid con una clave de más",
            MODEL_INVALID,
            format!("{populations}/1/genotype_model"),
        ),
        set(
            "coeficiente menor que 0",
            MODEL_INVALID,
            format!("{populations}/1/genotype_model/inbreeding_coefficient"),
            json!(-0.1),
        ),
        set(
            "coeficiente mayor que 1",
            MODEL_INVALID,
            format!("{populations}/1/genotype_model/inbreeding_coefficient"),
            json!(1.5),
        ),
        set(
            "coeficiente textual",
            MODEL_INVALID,
            format!("{populations}/1/genotype_model/inbreeding_coefficient"),
            json!("0.5"),
        ),
        set(
            "coeficiente nulo",
            MODEL_INVALID,
            format!("{populations}/1/genotype_model/inbreeding_coefficient"),
            json!(null),
        ),
        set(
            "loci que no son una lista",
            "loci estadísticos inválidos",
            p0("/loci"),
            json!({}),
        ),
        case("un locus repetido", "locus estadístico repetido", {
            let populations = populations.clone();
            move |document| {
                let locus = pointer_value(document, &format!("{populations}/0/loci/0"));
                append(document, &format!("{populations}/0/loci"), locus);
            }
        }),
    ];
    let mut cases = Vec::from(cases);
    let coefficient = format!("{populations}/1/genotype_model/inbreeding_coefficient");
    cases.push(set(
        "coeficiente justo por encima de 1",
        MODEL_INVALID,
        coefficient.clone(),
        json!(next_up(1.0)),
    ));
    cases.push(set(
        "coeficiente negativo mínimo",
        MODEL_INVALID,
        coefficient,
        json!(tiny_negative()),
    ));
    assert_rejected(&package, no_base, &cases);

    // El lado aceptado: ploidías 1 y 10, el coeficiente en sus dos extremos y
    // un locus que solo difiere de otro en una de las tres componentes de su
    // clave (secuencia, posición, base de referencia).
    let near_duplicate = |label: &str, edit: fn(&mut Value)| {
        let populations = populations.clone();
        ok_with(label.to_owned(), move |document| {
            let mut locus = pointer_value(document, &format!("{populations}/0/loci/0"));
            edit(&mut locus);
            append(document, &format!("{populations}/0/loci"), locus);
        })
    };
    let accepted = [
        ok(
            "coeficiente 0",
            format!("{populations}/1/genotype_model/inbreeding_coefficient"),
            json!(0.0),
        ),
        ok(
            "coeficiente 1",
            format!("{populations}/1/genotype_model/inbreeding_coefficient"),
            json!(1.0),
        ),
        near_duplicate("otra secuencia", |locus| {
            locus["contig"] = json!("scaffold-1")
        }),
        near_duplicate("otra posición", |locus| {
            locus["position_1based"] = json!(9)
        }),
        near_duplicate("otra base de referencia", |locus| {
            locus["ref"] = json!("C");
            locus["alleles"][0]["allele"] = json!("C");
        }),
    ];
    assert_accepted(&package, no_base, &accepted);
    // Las ploidías extremas son poblaciones sintéticas (la base las añade).
    let mut cases = Vec::new();
    for (index, label) in [(3, "ploidía 1"), (4, "ploidía 10")] {
        cases.push(set(
            format!("{label}: la base se acepta y se altera la clase del modelo"),
            MODEL_INVALID,
            format!("{populations}/{index}/genotype_model"),
            json!({"kind": "f_diploid", "inbreeding_coefficient": 0.5}),
        ));
    }
    assert_rejected(&package, with_synthetic_populations, &cases);
}

// ---------------------------------------------------------------------------
// La procedencia de las estadísticas
// ---------------------------------------------------------------------------

const PROVENANCE_INVALID: &str = "procedencia de estadísticas iniciales inválida";
const ALGORITHMS_INVALID: &str = "algoritmos de genotipo inválidos";
const CONSTRAINT_INVALID: &str = "restricción SNV explícita inválida";

fn provenance(suffix: &str) -> String {
    format!("{STATS}/provenance{suffix}")
}

#[test]
fn the_statistics_provenance_is_judged() {
    let package = statistics_package();
    let manifest = package.json("manifest.json");
    let mut cases = vec![
        add_key(
            "una clave de más en la procedencia",
            PROVENANCE_INVALID,
            provenance(""),
        ),
        remove(
            "una clave de menos en la procedencia",
            PROVENANCE_INVALID,
            provenance("/explicit_snp_constraints"),
        ),
        rename(
            "una clave desconocida en la procedencia",
            PROVENANCE_INVALID,
            provenance(""),
            "missing_locus_policy",
        ),
        rename(
            "los algoritmos con otro nombre: lo para la lista cerrada de claves",
            PROVENANCE_INVALID,
            provenance(""),
            "genotype_algorithms",
        ),
        set(
            "referencia sin nombre",
            PROVENANCE_INVALID,
            provenance("/matrix_reference_id"),
            json!(""),
        ),
        set(
            "referencia con nombre numérico",
            PROVENANCE_INVALID,
            provenance("/matrix_reference_id"),
            json!(5),
        ),
        set(
            "política de ausentes de otra clase",
            PROVENANCE_INVALID,
            provenance("/missing_locus_policy"),
            json!("other"),
        ),
        set(
            "política de ausentes numérica",
            PROVENANCE_INVALID,
            provenance("/missing_locus_policy"),
            json!(5),
        ),
        set(
            "algoritmos que no son una lista",
            ALGORITHMS_INVALID,
            provenance("/genotype_algorithms"),
            json!("hwe-composition-v1"),
        ),
        set(
            "ningún algoritmo",
            ALGORITHMS_INVALID,
            provenance("/genotype_algorithms"),
            json!([]),
        ),
        set(
            "un algoritmo de más",
            ALGORITHMS_INVALID,
            provenance("/genotype_algorithms"),
            json!([
                "hwe-composition-v1",
                "f-diploid-genotype-v1",
                "phase-permutation-v1",
                "other"
            ]),
        ),
        set(
            "el primer algoritmo cambiado",
            ALGORITHMS_INVALID,
            provenance("/genotype_algorithms/0"),
            json!("other"),
        ),
        set(
            "el segundo algoritmo cambiado",
            ALGORITHMS_INVALID,
            provenance("/genotype_algorithms/1"),
            json!("other"),
        ),
        set(
            "el tercer algoritmo cambiado",
            ALGORITHMS_INVALID,
            provenance("/genotype_algorithms/2"),
            json!("other"),
        ),
        set(
            "restricciones que no son una lista",
            "restricciones SNV explícitas inválidas",
            provenance("/explicit_snp_constraints"),
            json!({}),
        ),
        set(
            "una restricción que no es un objeto",
            CONSTRAINT_INVALID,
            provenance("/explicit_snp_constraints"),
            json!([5]),
        ),
    ];
    // La huella de la referencia y la de la matriz pasan por el mismo ayudante
    // y se prueban las dos: una por llamante.
    for field in [
        "matrix_reference_fingerprint",
        "matrix_canonical_tsv_sha256",
    ] {
        let digest = pointer_value(&manifest, &provenance(&format!("/{field}")));
        let digest = digest.as_str().expect("huella").to_owned();
        let pointer = provenance(&format!("/{field}"));
        cases.push(set(
            format!("{field} de 63 caracteres"),
            PROVENANCE_INVALID,
            pointer.clone(),
            json!(digest[..63]),
        ));
        cases.push(set(
            format!("{field} de 65 caracteres"),
            PROVENANCE_INVALID,
            pointer.clone(),
            json!(format!("{digest}0")),
        ));
        cases.push(set(
            format!("{field} con un carácter que no es hexadecimal"),
            PROVENANCE_INVALID,
            pointer.clone(),
            json!("g".repeat(64)),
        ));
        assert_ne!(
            digest.to_uppercase(),
            digest,
            "el fixture tiene que llevar letras"
        );
        cases.push(set(
            format!("{field} en mayúsculas"),
            PROVENANCE_INVALID,
            pointer.clone(),
            json!(digest.to_uppercase()),
        ));
        cases.push(set(
            format!("{field} numérica"),
            PROVENANCE_INVALID,
            pointer,
            json!(5),
        ));
    }
    assert_rejected(&package, no_base, &cases);
}

/// Una restricción explícita coherente, que el paquete generado no produce
/// (no declara cambios directos sobre loci de la matriz).
fn constraint() -> Value {
    json!({
        "population_id": "P1",
        "individual_id": "I001",
        "event_id": "ev-snp",
        "phase": 0,
        "matrix_variant_id": "VAR_C7V4CZRTKPH4LK3OYGCOEHDVTI",
        "sampled_alternate": "C",
        "explicit_alternate": "G",
        "incompatible": false,
    })
}

fn with_one_constraint(document: &mut Value) {
    append(
        document,
        &provenance("/explicit_snp_constraints"),
        constraint(),
    );
}

fn constraint_pointer(suffix: &str) -> String {
    provenance(&format!("/explicit_snp_constraints/0{suffix}"))
}

#[test]
fn the_explicit_snp_constraints_are_judged() {
    let package = statistics_package();
    let pointer = |suffix: &str| constraint_pointer(suffix);
    let mut cases = vec![
        add_key(
            "una clave de más en la restricción",
            CONSTRAINT_INVALID,
            pointer(""),
        ),
        remove(
            "una clave de menos en la restricción",
            CONSTRAINT_INVALID,
            pointer("/incompatible"),
        ),
        rename(
            "una clave desconocida en la restricción",
            CONSTRAINT_INVALID,
            pointer(""),
            "incompatible",
        ),
        set("fase 10", CONSTRAINT_INVALID, pointer("/phase"), json!(10)),
        set(
            "fase negativa",
            CONSTRAINT_INVALID,
            pointer("/phase"),
            json!(-1),
        ),
        set(
            "fase textual",
            CONSTRAINT_INVALID,
            pointer("/phase"),
            json!("0"),
        ),
        set(
            "fase fraccionaria",
            CONSTRAINT_INVALID,
            pointer("/phase"),
            json!(1.5),
        ),
        set(
            "alterno muestreado fuera de ACGT",
            CONSTRAINT_INVALID,
            pointer("/sampled_alternate"),
            json!("N"),
        ),
        set(
            "alterno muestreado numérico",
            CONSTRAINT_INVALID,
            pointer("/sampled_alternate"),
            json!(5),
        ),
        set(
            "alterno muestreado en minúscula",
            CONSTRAINT_INVALID,
            pointer("/sampled_alternate"),
            json!("c"),
        ),
        set(
            "alterno explícito nulo",
            CONSTRAINT_INVALID,
            pointer("/explicit_alternate"),
            json!(null),
        ),
        set(
            "alterno explícito fuera de ACGT",
            CONSTRAINT_INVALID,
            pointer("/explicit_alternate"),
            json!("N"),
        ),
        set(
            "alterno explícito de dos bases",
            CONSTRAINT_INVALID,
            pointer("/explicit_alternate"),
            json!("GG"),
        ),
        set(
            "alterno explícito numérico",
            CONSTRAINT_INVALID,
            pointer("/explicit_alternate"),
            json!(5),
        ),
        set(
            "incompatible textual",
            CONSTRAINT_INVALID,
            pointer("/incompatible"),
            json!("false"),
        ),
        set(
            "incompatible nulo",
            CONSTRAINT_INVALID,
            pointer("/incompatible"),
            json!(null),
        ),
    ];
    // Los cuatro textos no vacíos: vacío y no textual, uno por uno.
    for field in [
        "population_id",
        "individual_id",
        "event_id",
        "matrix_variant_id",
    ] {
        cases.push(set(
            format!("{field} vacío"),
            CONSTRAINT_INVALID,
            pointer(&format!("/{field}")),
            json!(""),
        ));
        cases.push(set(
            format!("{field} numérico"),
            CONSTRAINT_INVALID,
            pointer(&format!("/{field}")),
            json!(5),
        ));
    }
    // La clave de unicidad es (población, individuo, evento): repetirla se
    // rechaza aunque lo demás difiera, y difiriendo en una sola componente se
    // acepta.
    cases.push(case(
        "la misma restricción dos veces",
        "restricción SNV explícita repetida",
        |document| {
            append(
                document,
                &provenance("/explicit_snp_constraints"),
                constraint(),
            )
        },
    ));
    cases.push(case(
        "la misma clave con otra fase y otra variante",
        "restricción SNV explícita repetida",
        |document| {
            let mut other = constraint();
            other["phase"] = json!(3);
            other["matrix_variant_id"] = json!("VAR_OTHER");
            other["incompatible"] = json!(true);
            append(document, &provenance("/explicit_snp_constraints"), other);
        },
    ));
    assert_rejected(&package, with_one_constraint, &cases);

    let different_in = |label: &'static str, field: &'static str| {
        ok_with(label, move |document| {
            let mut other = constraint();
            other[field] = json!("OTHER");
            append(document, &provenance("/explicit_snp_constraints"), other);
        })
    };
    let accepted = [
        ok("fase 9", pointer("/phase"), json!(9)),
        ok(
            "alterno muestreado nulo",
            pointer("/sampled_alternate"),
            json!(null),
        ),
        ok(
            "alterno muestreado T",
            pointer("/sampled_alternate"),
            json!("T"),
        ),
        ok(
            "incompatible verdadero",
            pointer("/incompatible"),
            json!(true),
        ),
        different_in("otra población", "population_id"),
        different_in("otro individuo", "individual_id"),
        different_in("otro evento", "event_id"),
    ];
    assert_accepted(&package, with_one_constraint, &accepted);
}

// ---------------------------------------------------------------------------
// Los loci de las estadísticas
// ---------------------------------------------------------------------------

const LOCUS_INVALID: &str = "locus estadístico inválido";
const TOTAL_INVALID: &str = "total de copias inválido";
const COPIES_SUM: &str = "las copias observadas y ausentes no suman el total";
const ALLELES_INVALID: &str = "alelos estadísticos inválidos";
const ALLELE_INVALID: &str = "alelo estadístico inválido";
const ALLELE_ORDER: &str = "clase u orden de alelo estadístico inválido";
const ALLELE_DUPLICATED: &str = "alelo estadístico duplicado o referencia incoherente";
const ALLELE_PROVENANCE: &str = "procedencia de alelo alternativo inválida";
const FREQUENCY_RANGE: &str = "frecuencia estadística no finita o fuera de rango";
const FREQUENCY_BITS: &str = "frecuencia realizada no coincide con recuento y denominador";
const REQUESTED_SUM: &str = "las frecuencias objetivo no suman uno";
const GENOTYPE_INVALID: &str = "genotipo estadístico inválido";
const DOSAGES_INVALID: &str = "dosis genotípicas inválidas";
const DOSAGES_BAD: &str = "dosis genotípicas inválidas o repetidas";
const PROBABILITY_SUM: &str = "probabilidades genotípicas no suman uno";

/// El primer locus de la primera población (ploidía 2, tres individuos).
fn l0(suffix: &str) -> String {
    stats(&format!("/populations/0/loci/0{suffix}"))
}

/// El locus sintético que añaden las bases (`with_*_locus`), tras los tres que
/// publica el programa.
fn l3(suffix: &str) -> String {
    stats(&format!("/populations/0/loci/3{suffix}"))
}

/// Pone los recuentos realizados de los genotipos de un locus de dos alelos
/// (por dosis) y deja a 0 los demás.
fn set_genotype_counts(document: &mut Value, locus: &str, counts: &[([u64; 2], u64)]) {
    let genotypes = document
        .pointer_mut(&format!("{locus}/genotypes"))
        .and_then(Value::as_array_mut)
        .expect("genotipos");
    for genotype in genotypes {
        let dosages = genotype["allele_dosages"].clone();
        let count = counts
            .iter()
            .find(|(candidate, _)| json!(candidate) == dosages)
            .map_or(0, |(_, count)| *count);
        genotype["realized_count"] = json!(count);
    }
}

/// El negativo más pequeño representable en valor absoluto (`-5e-324`): el
/// primer valor fuera por debajo de un límite inferior que vale 0.
fn tiny_negative() -> f64 {
    -f64::from_bits(1)
}

fn next_up(value: f64) -> f64 {
    assert!(value > 0.0 && value.is_finite());
    f64::from_bits(value.to_bits() + 1)
}

/// Un locus de tres alelos que cumple el contrato (ploidía 2, tres
/// individuos, todos completos). **Sintético**: la matriz generada solo trae
/// loci de dos alelos. Los genotipos realizados cuadran con los recuentos de
/// los alelos (A = 2, G = 3, C = 1): un individuo A/G, otro A/C y otro G/G, de
/// modo que la suma de dosis por recuento de cada alelo es exactamente su
/// recuento y las observaciones completas por la ploidía son las seis copias
/// observadas.
fn synthetic_three_allele_locus() -> Value {
    let genotypes = [
        ([2, 0, 0], 0.25, 0),
        ([1, 1, 0], 0.25, 1),
        ([1, 0, 1], 0.25, 1),
        ([0, 2, 0], 0.125, 1),
        ([0, 1, 1], 0.0625, 0),
        ([0, 0, 2], 0.0625, 0),
    ]
    .map(|(dosages, probability, realized)| {
        json!({
            "allele_dosages": dosages,
            "expected_probability": probability,
            "expected_count": probability * 3.0,
            "realized_count": realized,
        })
    });
    json!({
        "contig": "chr1",
        "position_1based": 300,
        "ref": "A",
        "total_copies": 6,
        "observed_copies": 6,
        "missing_by_deletion": 0,
        "alleles": [
            {"allele": "A", "kind": "ref", "requested_frequency": 0.3,
             "realized_count": 2, "realized_frequency": 2.0 / 6.0},
            {"allele": "G", "kind": "alt", "variant_id": "VAR_SYNTHETIC_1",
             "source_id": "supplied", "requested_frequency": 0.5,
             "realized_count": 3, "realized_frequency": 0.5},
            {"allele": "C", "kind": "alt", "variant_id": "VAR_SYNTHETIC_2",
             "source_id": "supplied", "requested_frequency": 0.2,
             "realized_count": 1, "realized_frequency": 1.0 / 6.0},
        ],
        "complete_genotype_observations": 3,
        "genotypes": genotypes,
    })
}

fn with_three_allele_locus(document: &mut Value) {
    append(
        document,
        &stats("/populations/0/loci"),
        synthetic_three_allele_locus(),
    );
}

/// Un locus sin ninguna copia observada: las seis copias se han perdido por
/// deleción. **Sintético**: lo que ocurre cuando una deleción directa se come
/// el sitio de todos los individuos.
fn synthetic_unobserved_locus() -> Value {
    json!({
        "contig": "chr1",
        "position_1based": 300,
        "ref": "A",
        "total_copies": 6,
        "observed_copies": 0,
        "missing_by_deletion": 6,
        "alleles": [
            {"allele": "A", "kind": "ref", "requested_frequency": 0.5,
             "realized_count": 0, "realized_frequency": 0.0},
            {"allele": "G", "kind": "alt", "variant_id": "VAR_SYNTHETIC_1",
             "source_id": "supplied", "requested_frequency": 0.5,
             "realized_count": 0, "realized_frequency": 0.0},
        ],
        "complete_genotype_observations": 0,
        "genotypes": [
            {"allele_dosages": [2, 0], "expected_probability": 0.25, "expected_count": 0.75, "realized_count": 0},
            {"allele_dosages": [1, 1], "expected_probability": 0.5, "expected_count": 1.5, "realized_count": 0},
            {"allele_dosages": [0, 2], "expected_probability": 0.25, "expected_count": 0.75, "realized_count": 0},
        ],
    })
}

fn with_unobserved_locus(document: &mut Value) {
    append(
        document,
        &stats("/populations/0/loci"),
        synthetic_unobserved_locus(),
    );
}

#[test]
fn the_statistics_loci_are_judged() {
    let package = statistics_package();
    let manifest = package.json("manifest.json");
    // El control positivo ya cubre el límite «tantas observaciones genotípicas
    // como individuos»: el locus publicado las tiene completas.
    assert_eq!(
        pointer_value(&manifest, &l0("/complete_genotype_observations")),
        pointer_value(&manifest, &stats("/populations/0/individual_count")),
    );
    let realized_frequency = |index: usize| {
        pointer_value(
            &manifest,
            &l0(&format!("/alleles/{index}/realized_frequency")),
        )
        .as_f64()
        .expect("frecuencia")
    };
    let requested = |index: usize| {
        pointer_value(
            &manifest,
            &l0(&format!("/alleles/{index}/requested_frequency")),
        )
        .as_f64()
        .expect("frecuencia")
    };
    let probability = pointer_value(&manifest, &l0("/genotypes/0/expected_probability"))
        .as_f64()
        .expect("probabilidad");
    let (frequency_1, requested_0) = (realized_frequency(1), requested(0));
    let cases = [
        set(
            "un locus que no es un objeto",
            LOCUS_INVALID,
            l0(""),
            json!(5),
        ),
        add_key("una clave de más en el locus", LOCUS_INVALID, l0("")),
        remove(
            "una clave de menos en el locus",
            LOCUS_INVALID,
            l0("/genotypes"),
        ),
        rename(
            "una clave desconocida en el locus",
            LOCUS_INVALID,
            l0(""),
            "genotypes",
        ),
        set("secuencia vacía", LOCUS_INVALID, l0("/contig"), json!("")),
        set("secuencia numérica", LOCUS_INVALID, l0("/contig"), json!(5)),
        set(
            "posición 0",
            LOCUS_INVALID,
            l0("/position_1based"),
            json!(0),
        ),
        set(
            "posición textual",
            LOCUS_INVALID,
            l0("/position_1based"),
            json!("8"),
        ),
        set(
            "posición negativa",
            LOCUS_INVALID,
            l0("/position_1based"),
            json!(-1),
        ),
        set(
            "base de referencia fuera de ACGT",
            LOCUS_INVALID,
            l0("/ref"),
            json!("N"),
        ),
        set(
            "base de referencia en minúscula",
            LOCUS_INVALID,
            l0("/ref"),
            json!("a"),
        ),
        set(
            "base de referencia de dos letras",
            LOCUS_INVALID,
            l0("/ref"),
            json!("AC"),
        ),
        set(
            "base de referencia numérica",
            LOCUS_INVALID,
            l0("/ref"),
            json!(5),
        ),
        set(
            "total de copias que no es un número",
            TOTAL_INVALID,
            l0("/total_copies"),
            json!("6"),
        ),
        set(
            "copias observadas que no son un número",
            "copias observadas inválidas",
            l0("/observed_copies"),
            json!("6"),
        ),
        set(
            "copias ausentes que no son un número",
            "copias ausentes inválidas",
            l0("/missing_by_deletion"),
            json!("0"),
        ),
        case(
            "un total que no es individuos por ploidía",
            COPIES_SUM,
            |document| {
                replace_at(document, &l0("/total_copies"), json!(7));
                replace_at(document, &l0("/observed_copies"), json!(7));
            },
        ),
        case(
            "observadas y ausentes que suman el total solo al desbordar",
            COPIES_SUM,
            |document| {
                replace_at(document, &l0("/observed_copies"), json!(u64::MAX));
                replace_at(document, &l0("/missing_by_deletion"), json!(7));
            },
        ),
        set(
            "observadas y ausentes que no suman el total",
            COPIES_SUM,
            l0("/observed_copies"),
            json!(5),
        ),
        case(
            "individuos por ploidía que desbordan",
            TOTAL_INVALID,
            |document| {
                replace_at(
                    document,
                    &stats("/populations/0/individual_count"),
                    json!(u64::MAX),
                );
            },
        ),
        set(
            "alelos que no son una lista",
            ALLELES_INVALID,
            l0("/alleles"),
            json!({}),
        ),
        remove("un solo alelo", ALLELES_INVALID, l0("/alleles/1")),
        set("ningún alelo", ALLELES_INVALID, l0("/alleles"), json!([])),
        set(
            "un alelo que no es un objeto",
            ALLELE_INVALID,
            l0("/alleles/0"),
            json!(5),
        ),
        set(
            "un alelo de referencia que declara ser alterno",
            ALLELE_ORDER,
            l0("/alleles/0/kind"),
            json!("alt"),
        ),
        set(
            "un alelo alterno que declara ser de referencia",
            ALLELE_ORDER,
            l0("/alleles/1/kind"),
            json!("ref"),
        ),
        set(
            "una clase de alelo numérica",
            ALLELE_ORDER,
            l0("/alleles/0/kind"),
            json!(5),
        ),
        set(
            "una clase de alelo desconocida",
            ALLELE_ORDER,
            l0("/alleles/1/kind"),
            json!("other"),
        ),
        add_key(
            "una clave de más en el alelo de referencia",
            ALLELE_INVALID,
            l0("/alleles/0"),
        ),
        add_key(
            "una clave de más en el alelo alterno",
            ALLELE_INVALID,
            l0("/alleles/1"),
        ),
        remove(
            "un alelo alterno sin source_id",
            ALLELE_INVALID,
            l0("/alleles/1/source_id"),
        ),
        remove(
            "un alelo de referencia sin frecuencia objetivo",
            ALLELE_INVALID,
            l0("/alleles/0/requested_frequency"),
        ),
        rename(
            "una clave desconocida en el alelo alterno",
            ALLELE_INVALID,
            l0("/alleles/1"),
            "variant_id",
        ),
        rename(
            "una clave desconocida en el alelo de referencia",
            ALLELE_INVALID,
            l0("/alleles/0"),
            "realized_frequency",
        ),
        set(
            "un alelo alterno fuera de ACGT",
            ALLELE_INVALID,
            l0("/alleles/1/allele"),
            json!("N"),
        ),
        set(
            "un alelo alterno de dos letras",
            ALLELE_INVALID,
            l0("/alleles/1/allele"),
            json!("GG"),
        ),
        set(
            "un alelo alterno en minúscula",
            ALLELE_INVALID,
            l0("/alleles/1/allele"),
            json!("g"),
        ),
        set(
            "un alelo alterno numérico",
            ALLELE_INVALID,
            l0("/alleles/1/allele"),
            json!(5),
        ),
        set(
            "un alelo alterno igual a la referencia",
            ALLELE_DUPLICATED,
            l0("/alleles/1/allele"),
            json!("A"),
        ),
        set(
            "un alelo de referencia que no es la base de referencia",
            ALLELE_DUPLICATED,
            l0("/alleles/0/allele"),
            json!("C"),
        ),
        set(
            "un alelo alterno sin variante",
            ALLELE_PROVENANCE,
            l0("/alleles/1/variant_id"),
            json!(""),
        ),
        set(
            "un alelo alterno con variante numérica",
            ALLELE_PROVENANCE,
            l0("/alleles/1/variant_id"),
            json!(5),
        ),
        set(
            "un alelo alterno con fuente numérica",
            ALLELE_PROVENANCE,
            l0("/alleles/1/source_id"),
            json!(5),
        ),
        set(
            "un alelo alterno con fuente nula",
            ALLELE_PROVENANCE,
            l0("/alleles/1/source_id"),
            json!(null),
        ),
        set(
            "un recuento alélico que no es un número",
            "recuento alélico inválido",
            l0("/alleles/0/realized_count"),
            json!("2"),
        ),
        set(
            "una frecuencia objetivo que no es un número",
            "frecuencia objetivo inválida",
            l0("/alleles/1/requested_frequency"),
            json!("0.5"),
        ),
        set(
            "una frecuencia realizada que no es un número",
            "frecuencia realizada inválida",
            l0("/alleles/0/realized_frequency"),
            json!("0.5"),
        ),
        set(
            "frecuencia objetivo mayor que 1",
            FREQUENCY_RANGE,
            l0("/alleles/0/requested_frequency"),
            json!(1.5),
        ),
        set(
            "frecuencia objetivo negativa",
            FREQUENCY_RANGE,
            l0("/alleles/1/requested_frequency"),
            json!(-0.1),
        ),
        set(
            "frecuencia realizada mayor que 1",
            FREQUENCY_RANGE,
            l0("/alleles/0/realized_frequency"),
            json!(1.5),
        ),
        set(
            "frecuencia realizada negativa",
            FREQUENCY_RANGE,
            l0("/alleles/1/realized_frequency"),
            json!(-0.1),
        ),
        set(
            "frecuencia realizada un bit por encima de recuento entre observadas",
            FREQUENCY_BITS,
            l0("/alleles/1/realized_frequency"),
            json!(next_up(frequency_1)),
        ),
        set(
            "frecuencia realizada que no es recuento entre observadas",
            FREQUENCY_BITS,
            l0("/alleles/0/realized_frequency"),
            json!(0.99),
        ),
        case(
            "recuentos que no suman las copias observadas",
            "los recuentos alélicos no suman las copias observadas",
            |document| {
                for index in 0..2 {
                    replace_at(
                        document,
                        &l0(&format!("/alleles/{index}/realized_count")),
                        json!(1),
                    );
                    replace_at(
                        document,
                        &l0(&format!("/alleles/{index}/realized_frequency")),
                        json!(1.0 / 6.0),
                    );
                }
            },
        ),
        case(
            "recuentos alélicos que desbordan",
            "recuento alélico desbordado",
            |document| {
                let total = (u64::MAX / 2) * 2;
                replace_at(
                    document,
                    &stats("/populations/0/individual_count"),
                    json!(u64::MAX / 2),
                );
                replace_at(document, &l0("/total_copies"), json!(total));
                replace_at(document, &l0("/observed_copies"), json!(total));
                for index in 0..2 {
                    replace_at(
                        document,
                        &l0(&format!("/alleles/{index}/realized_count")),
                        json!(1_u64 << 63),
                    );
                    replace_at(
                        document,
                        &l0(&format!("/alleles/{index}/realized_frequency")),
                        json!(0.5),
                    );
                    replace_at(
                        document,
                        &l0(&format!("/alleles/{index}/requested_frequency")),
                        json!(0.5),
                    );
                }
            },
        ),
        set(
            "frecuencias objetivo que suman más de uno",
            REQUESTED_SUM,
            l0("/alleles/0/requested_frequency"),
            json!(requested_0 + 1.5e-12),
        ),
        set(
            "frecuencias objetivo que suman menos de uno",
            REQUESTED_SUM,
            l0("/alleles/0/requested_frequency"),
            json!(requested_0 - 1.5e-12),
        ),
        set(
            "observaciones genotípicas que no son un número",
            "observaciones genotípicas inválidas",
            l0("/complete_genotype_observations"),
            json!("2"),
        ),
        set(
            "más observaciones genotípicas que individuos",
            "observaciones genotípicas fuera de rango",
            l0("/complete_genotype_observations"),
            json!(4),
        ),
        set(
            "genotipos que no son una lista",
            "genotipos estadísticos inválidos",
            l0("/genotypes"),
            json!({}),
        ),
        set(
            "ningún genotipo",
            "genotipos estadísticos inválidos",
            l0("/genotypes"),
            json!([]),
        ),
        set(
            "un genotipo que no es un objeto",
            GENOTYPE_INVALID,
            l0("/genotypes/0"),
            json!(5),
        ),
        add_key(
            "una clave de más en el genotipo",
            GENOTYPE_INVALID,
            l0("/genotypes/0"),
        ),
        remove(
            "una clave de menos en el genotipo",
            GENOTYPE_INVALID,
            l0("/genotypes/0/realized_count"),
        ),
        rename(
            "una clave desconocida en el genotipo",
            GENOTYPE_INVALID,
            l0("/genotypes/0"),
            "expected_count",
        ),
        set(
            "dosis que no son una lista",
            DOSAGES_INVALID,
            l0("/genotypes/0/allele_dosages"),
            json!("2,0"),
        ),
        set(
            "una dosis que no es un número",
            DOSAGES_INVALID,
            l0("/genotypes/0/allele_dosages"),
            json!([1, "1"]),
        ),
        set(
            "dosis que desbordan al sumarse",
            "dosis genotípicas desbordadas",
            l0("/genotypes/0/allele_dosages"),
            json!([u64::MAX, 1]),
        ),
        set(
            "dosis de menos",
            DOSAGES_BAD,
            l0("/genotypes/0/allele_dosages"),
            json!([2]),
        ),
        set(
            "dosis de más con la ploidía bien sumada",
            DOSAGES_BAD,
            l0("/genotypes/0/allele_dosages"),
            json!([1, 1, 0]),
        ),
        set(
            "dosis que suman más que la ploidía",
            DOSAGES_BAD,
            l0("/genotypes/0/allele_dosages"),
            json!([1, 2]),
        ),
        set(
            "dosis que suman menos que la ploidía",
            DOSAGES_BAD,
            l0("/genotypes/0/allele_dosages"),
            json!([1, 0]),
        ),
        case(
            "dos genotipos con las mismas dosis",
            DOSAGES_BAD,
            |document| {
                let dosages = pointer_value(document, &l0("/genotypes/0/allele_dosages"));
                replace_at(document, &l0("/genotypes/1/allele_dosages"), dosages);
            },
        ),
        set(
            "probabilidad que no es un número",
            "probabilidad genotípica inválida",
            l0("/genotypes/0/expected_probability"),
            json!("0.1"),
        ),
        set(
            "probabilidad mayor que 1",
            "probabilidad genotípica inválida",
            l0("/genotypes/0/expected_probability"),
            json!(1.5),
        ),
        set(
            "probabilidad negativa",
            "probabilidad genotípica inválida",
            l0("/genotypes/0/expected_probability"),
            json!(-0.1),
        ),
        set(
            "recuento esperado que no es un número",
            "recuento genotípico esperado inválido",
            l0("/genotypes/0/expected_count"),
            json!("0.1"),
        ),
        set(
            "recuento esperado negativo",
            "recuento genotípico esperado inválido",
            l0("/genotypes/0/expected_count"),
            json!(-1.0),
        ),
        set(
            "recuento esperado nulo",
            "recuento genotípico esperado inválido",
            l0("/genotypes/0/expected_count"),
            json!(null),
        ),
        set(
            "recuento realizado que no es un número",
            "recuento genotípico realizado inválido",
            l0("/genotypes/0/realized_count"),
            json!("0"),
        ),
        set(
            "recuento realizado negativo",
            "recuento genotípico realizado inválido",
            l0("/genotypes/0/realized_count"),
            json!(-1),
        ),
        case(
            "recuentos genotípicos realizados que desbordan",
            "recuento genotípico realizado desbordado",
            |document| {
                for index in 0..2 {
                    replace_at(
                        document,
                        &l0(&format!("/genotypes/{index}/realized_count")),
                        json!(u64::MAX),
                    );
                }
            },
        ),
        set(
            "probabilidades que suman más de uno",
            PROBABILITY_SUM,
            l0("/genotypes/0/expected_probability"),
            json!(probability + 2e-12),
        ),
        set(
            "probabilidades que suman menos de uno",
            PROBABILITY_SUM,
            l0("/genotypes/0/expected_probability"),
            json!(probability - 2e-12),
        ),
        case(
            "recuentos genotípicos que no suman las observaciones completas",
            "recuentos genotípicos realizados incoherentes",
            |document| {
                let realized = pointer_value(document, &l0("/genotypes/0/realized_count"))
                    .as_u64()
                    .expect("recuento");
                replace_at(
                    document,
                    &l0("/genotypes/0/realized_count"),
                    json!(realized + 1),
                );
            },
        ),
    ];
    // Los límites, justo fuera de cada rango y de cada tolerancia (el primer
    // valor por encima de 1, el primer negativo por debajo de 0 y la primera
    // desviación representable por encima de 1e-12: 4504 epsilon, frente a
    // 4503 epsilon, que se acepta).
    let mut cases = Vec::from(cases);
    for (index, name) in [(0, "referencia"), (1, "alterno")] {
        for (label, field) in [
            ("objetivo", "requested_frequency"),
            ("realizada", "realized_frequency"),
        ] {
            let pointer = l0(&format!("/alleles/{index}/{field}"));
            cases.push(set(
                format!("frecuencia {label} del alelo de {name} justo por encima de 1"),
                FREQUENCY_RANGE,
                pointer.clone(),
                json!(next_up(1.0)),
            ));
            cases.push(set(
                format!("frecuencia {label} del alelo de {name} negativa mínima"),
                FREQUENCY_RANGE,
                pointer,
                json!(tiny_negative()),
            ));
        }
    }
    cases.push(set(
        "probabilidad justo por encima de 1",
        "probabilidad genotípica inválida",
        l0("/genotypes/0/expected_probability"),
        json!(next_up(1.0)),
    ));
    cases.push(set(
        "probabilidad negativa mínima",
        "probabilidad genotípica inválida",
        l0("/genotypes/0/expected_probability"),
        json!(tiny_negative()),
    ));
    cases.push(set(
        "recuento esperado negativo mínimo",
        "recuento genotípico esperado inválido",
        l0("/genotypes/0/expected_count"),
        json!(tiny_negative()),
    ));
    cases.push(case(
        "frecuencias objetivo con la primera desviación por encima de 1e-12",
        REQUESTED_SUM,
        |document| {
            replace_at(document, &l0("/alleles/0/requested_frequency"), json!(0.5));
            replace_at(
                document,
                &l0("/alleles/1/requested_frequency"),
                json!(0.5 + 4504.0 * f64::EPSILON),
            );
        },
    ));
    cases.push(case(
        "probabilidades con la primera desviación por encima de 1e-12",
        PROBABILITY_SUM,
        |document| {
            for (index, probability) in [0.5, 0.5 + 4504.0 * f64::EPSILON, 0.0]
                .into_iter()
                .enumerate()
            {
                replace_at(
                    document,
                    &l0(&format!("/genotypes/{index}/expected_probability")),
                    json!(probability),
                );
            }
        },
    ));
    assert_rejected(&package, no_base, &cases);

    // El lado aceptado. Con dos copias ausentes (las observadas pasan a 4, con
    // sus recuentos y frecuencias), con las frecuencias y las probabilidades
    // en sus extremos y con una desviación bajo la tolerancia de 1e-12.
    let accepted = [
        // Un individuo pierde sus dos copias por deleción: quedan cuatro copias
        // observadas (A = 1, G = 3) y dos individuos completos, uno A/G y otro
        // G/G, que usan exactamente esas cuatro copias.
        ok_with("dos copias ausentes por deleción", |document| {
            replace_at(document, &l0("/observed_copies"), json!(4));
            replace_at(document, &l0("/missing_by_deletion"), json!(2));
            replace_at(document, &l0("/alleles/0/realized_count"), json!(1));
            replace_at(document, &l0("/alleles/0/realized_frequency"), json!(0.25));
            replace_at(document, &l0("/alleles/1/realized_count"), json!(3));
            replace_at(document, &l0("/alleles/1/realized_frequency"), json!(0.75));
            replace_at(document, &l0("/complete_genotype_observations"), json!(2));
            set_genotype_counts(document, &l0(""), &[([1, 1], 1), ([0, 2], 1)]);
        }),
        // Tres individuos G/G: ninguna copia de referencia, seis de alterno.
        ok_with(
            "frecuencias objetivo y realizadas en 0 y en 1",
            |document| {
                replace_at(document, &l0("/alleles/0/requested_frequency"), json!(0.0));
                replace_at(document, &l0("/alleles/0/realized_count"), json!(0));
                replace_at(document, &l0("/alleles/0/realized_frequency"), json!(0.0));
                replace_at(document, &l0("/alleles/1/requested_frequency"), json!(1.0));
                replace_at(document, &l0("/alleles/1/realized_count"), json!(6));
                replace_at(document, &l0("/alleles/1/realized_frequency"), json!(1.0));
                replace_at(document, &l0("/complete_genotype_observations"), json!(3));
                set_genotype_counts(document, &l0(""), &[([0, 2], 3)]);
            },
        ),
        ok_with("probabilidades 0 y 1 y recuentos esperados 0", |document| {
            for (index, probability) in [0.0, 1.0, 0.0].into_iter().enumerate() {
                replace_at(
                    document,
                    &l0(&format!("/genotypes/{index}/expected_probability")),
                    json!(probability),
                );
                replace_at(
                    document,
                    &l0(&format!("/genotypes/{index}/expected_count")),
                    json!(0.0),
                );
            }
        }),
        ok(
            "frecuencias objetivo con una desviación bajo la tolerancia",
            l0("/alleles/0/requested_frequency"),
            json!(requested_0 + 5e-13),
        ),
        ok(
            "probabilidades con una desviación bajo la tolerancia",
            l0("/genotypes/0/expected_probability"),
            json!(probability + 5e-13),
        ),
        ok_with(
            "frecuencias objetivo con la última desviación bajo 1e-12",
            |document| {
                replace_at(document, &l0("/alleles/0/requested_frequency"), json!(0.5));
                replace_at(
                    document,
                    &l0("/alleles/1/requested_frequency"),
                    json!(0.5 + 4503.0 * f64::EPSILON),
                );
            },
        ),
        ok_with(
            "probabilidades con la última desviación bajo 1e-12",
            |document| {
                for (index, probability) in [0.5, 0.5 + 4503.0 * f64::EPSILON, 0.0]
                    .into_iter()
                    .enumerate()
                {
                    replace_at(
                        document,
                        &l0(&format!("/genotypes/{index}/expected_probability")),
                        json!(probability),
                    );
                }
            },
        ),
        ok(
            "una fuente de alelo alterno vacía",
            l0("/alleles/1/source_id"),
            json!(""),
        ),
    ];
    assert_accepted(&package, no_base, &accepted);
}

/// Las ramas que solo un locus de tres alelos alcanza.
#[test]
fn the_statistics_three_allele_locus_is_judged() {
    let package = statistics_package();
    let cases = [
        set(
            "dos alelos alternos iguales",
            ALLELE_DUPLICATED,
            l3("/alleles/2/allele"),
            json!("G"),
        ),
        set(
            "un alelo alterno igual a la referencia",
            ALLELE_DUPLICATED,
            l3("/alleles/2/allele"),
            json!("A"),
        ),
        set(
            "dosis de menos que alelos",
            DOSAGES_BAD,
            l3("/genotypes/0/allele_dosages"),
            json!([1, 1]),
        ),
        set(
            "dosis de más que alelos",
            DOSAGES_BAD,
            l3("/genotypes/0/allele_dosages"),
            json!([1, 1, 0, 0]),
        ),
        set(
            "un tercer alelo con una frecuencia objetivo que rompe la suma",
            REQUESTED_SUM,
            l3("/alleles/2/requested_frequency"),
            json!(0.25),
        ),
        set(
            "probabilidades que no suman uno en los seis genotipos",
            PROBABILITY_SUM,
            l3("/genotypes/5/expected_probability"),
            json!(0.125),
        ),
    ];
    assert_rejected(&package, with_three_allele_locus, &cases);
}

/// La rama sin copias observadas: el denominador de la frecuencia es 0 y la
/// frecuencia realizada tiene que ser **exactamente +0.0**.
#[test]
fn the_statistics_locus_without_observed_copies_is_judged() {
    let package = statistics_package();
    let cases = [
        set(
            "frecuencia realizada -0.0 sin copias observadas",
            FREQUENCY_BITS,
            l3("/alleles/1/realized_frequency"),
            json!(-0.0),
        ),
        set(
            "frecuencia realizada 0.5 sin copias observadas",
            FREQUENCY_BITS,
            l3("/alleles/1/realized_frequency"),
            json!(0.5),
        ),
    ];
    assert_rejected(&package, with_unobserved_locus, &cases);
}

// ---------------------------------------------------------------------------
// La auditoría de la normalización
// ---------------------------------------------------------------------------

/// Todos los rechazos de la auditoría comparten este motivo: lo que
/// discrimina cada caso es que la alteración es **única** y que, si la guarda
/// no rechazara, el caso caería en la comparación efectiva/manifiesto con otro
/// texto.
const AUDIT_INVALID: &str = "population_snvs.normalization_audit inválida";

fn g(index: usize, suffix: &str) -> String {
    audit(&format!("/{index}{suffix}"))
}

fn swap_groups(first: usize, second: usize) -> impl Fn(&mut Value) {
    move |document| {
        document
            .pointer_mut(AUDIT)
            .and_then(Value::as_array_mut)
            .expect("auditoría")
            .swap(first, second);
    }
}

fn audit_value(alt: &str, variant: &str, source: &str) -> Value {
    json!({"alt": alt, "variant_id": variant, "source_id": source, "original": 0.1, "effective": 0.1})
}

fn set_values(values: Vec<Value>) -> impl Fn(&mut Value) {
    move |document| replace_at(document, &g(0, "/values"), Value::Array(values.clone()))
}

/// Inserta en `index` una copia del primer grupo con la base de referencia
/// cambiada.
fn insert_group_with_ref(index: usize, reference: &'static str) -> impl Fn(&mut Value) {
    move |document| {
        let mut group = pointer_value(document, &g(0, ""));
        group["ref"] = json!(reference);
        document
            .pointer_mut(AUDIT)
            .and_then(Value::as_array_mut)
            .expect("auditoría")
            .insert(index, group);
    }
}

#[test]
fn the_normalization_audit_groups_are_judged() {
    let package = statistics_package();
    let mut cases = vec![
        set("una auditoría nula", AUDIT_INVALID, AUDIT, json!(null)),
        set(
            "una auditoría que es un objeto",
            AUDIT_INVALID,
            AUDIT,
            json!({}),
        ),
        set("una auditoría textual", AUDIT_INVALID, AUDIT, json!("none")),
        set(
            "un grupo que no es un objeto",
            AUDIT_INVALID,
            g(0, ""),
            json!(5),
        ),
        add_key("una clave de más en el grupo", AUDIT_INVALID, g(0, "")),
        remove(
            "una clave de menos en el grupo",
            AUDIT_INVALID,
            g(0, "/decision"),
        ),
        rename(
            "una clave desconocida en el grupo",
            AUDIT_INVALID,
            g(0, ""),
            "decision",
        ),
        set(
            "población vacía",
            AUDIT_INVALID,
            g(0, "/population_id"),
            json!(""),
        ),
        set(
            "población numérica",
            AUDIT_INVALID,
            g(0, "/population_id"),
            json!(5),
        ),
        set("secuencia vacía", AUDIT_INVALID, g(0, "/contig"), json!("")),
        set(
            "secuencia numérica",
            AUDIT_INVALID,
            g(0, "/contig"),
            json!(5),
        ),
        set(
            "posición numérica",
            AUDIT_INVALID,
            g(0, "/position"),
            json!(8),
        ),
        set("posición 0", AUDIT_INVALID, g(0, "/position"), json!("0")),
        set(
            "posición que no es un número",
            AUDIT_INVALID,
            g(0, "/position"),
            json!("abc"),
        ),
        set(
            "posición negativa",
            AUDIT_INVALID,
            g(0, "/position"),
            json!("-1"),
        ),
        set(
            "posición vacía",
            AUDIT_INVALID,
            g(0, "/position"),
            json!(""),
        ),
        set(
            "posición que no cabe en 64 bits",
            AUDIT_INVALID,
            g(0, "/position"),
            json!("18446744073709551616"),
        ),
        set(
            "base de referencia fuera de ACGT",
            AUDIT_INVALID,
            g(0, "/ref"),
            json!("N"),
        ),
        set(
            "base de referencia de dos letras",
            AUDIT_INVALID,
            g(0, "/ref"),
            json!("AC"),
        ),
        set(
            "base de referencia en minúscula",
            AUDIT_INVALID,
            g(0, "/ref"),
            json!("a"),
        ),
        set(
            "base de referencia numérica",
            AUDIT_INVALID,
            g(0, "/ref"),
            json!(5),
        ),
        set(
            "tolerancia textual",
            AUDIT_INVALID,
            g(0, "/tolerance"),
            json!("1e-9"),
        ),
        set(
            "tolerancia negativa",
            AUDIT_INVALID,
            g(0, "/tolerance"),
            json!(-1e-9),
        ),
        set(
            "tolerancia mayor que 1e-6",
            AUDIT_INVALID,
            g(0, "/tolerance"),
            json!(2e-6),
        ),
        set(
            "tolerancia nula",
            AUDIT_INVALID,
            g(0, "/tolerance"),
            json!(null),
        ),
        set(
            "orden canónico distinto",
            AUDIT_INVALID,
            g(0, "/canonical_order"),
            json!("population_id,contig"),
        ),
        set(
            "orden canónico numérico",
            AUDIT_INVALID,
            g(0, "/canonical_order"),
            json!(5),
        ),
        set(
            "algoritmo de suma distinto",
            AUDIT_INVALID,
            g(0, "/sum_algorithm"),
            json!("kahan_sum_v1"),
        ),
        set(
            "algoritmo de suma numérico",
            AUDIT_INVALID,
            g(0, "/sum_algorithm"),
            json!(5),
        ),
        set(
            "frecuencia de referencia mayor que 1",
            AUDIT_INVALID,
            g(0, "/reference_frequency"),
            json!(1.5),
        ),
        set(
            "frecuencia de referencia negativa",
            AUDIT_INVALID,
            g(0, "/reference_frequency"),
            json!(-0.1),
        ),
        set(
            "frecuencia de referencia textual",
            AUDIT_INVALID,
            g(0, "/reference_frequency"),
            json!("0.5"),
        ),
        set(
            "decisión desconocida",
            AUDIT_INVALID,
            g(0, "/decision"),
            json!("other"),
        ),
        set(
            "decisión numérica",
            AUDIT_INVALID,
            g(0, "/decision"),
            json!(5),
        ),
        set(
            "valores que no son una lista",
            AUDIT_INVALID,
            g(0, "/values"),
            json!({}),
        ),
        set("ningún valor", AUDIT_INVALID, g(0, "/values"), json!([])),
        // El orden de los grupos: población, secuencia, posición (numérica) y
        // base de referencia, y sin repeticiones.
        case(
            "dos grupos en orden de posición descendente",
            AUDIT_INVALID,
            swap_groups(0, 1),
        ),
        case(
            "dos grupos en orden de secuencia descendente",
            AUDIT_INVALID,
            swap_groups(1, 2),
        ),
        case(
            "dos grupos en orden de población descendente",
            AUDIT_INVALID,
            swap_groups(2, 3),
        ),
        case(
            "dos grupos en orden de base de referencia descendente",
            AUDIT_INVALID,
            insert_group_with_ref(0, "C"),
        ),
        case("un grupo repetido", AUDIT_INVALID, |document| {
            let group = pointer_value(document, &g(0, ""));
            document
                .pointer_mut(AUDIT)
                .and_then(Value::as_array_mut)
                .expect("auditoría")
                .insert(1, group);
        }),
    ];
    // Los cuatro números no negativos, cada uno con su rechazo.
    for field in [
        "original_sum",
        "original_delta",
        "effective_sum",
        "effective_delta",
    ] {
        cases.push(set(
            format!("{field} textual"),
            AUDIT_INVALID,
            g(0, &format!("/{field}")),
            json!("0.5"),
        ));
        cases.push(set(
            format!("{field} negativo"),
            AUDIT_INVALID,
            g(0, &format!("/{field}")),
            json!(-0.5),
        ));
    }
    // La posición se escribe como el contrato (`^[1-9][0-9]*$`): ni cero a la
    // izquierda ni signo. (`str::parse::<u64>` admite las dos.)
    cases.push(set(
        "posición con un cero a la izquierda",
        AUDIT_INVALID,
        g(0, "/position"),
        json!("08"),
    ));
    cases.push(set(
        "posición con signo",
        AUDIT_INVALID,
        g(0, "/position"),
        json!("+8"),
    ));
    // Los límites, justo fuera: el primer valor por encima del máximo y el
    // primer negativo por debajo de 0.
    cases.push(set(
        "tolerancia justo por encima de 1e-6",
        AUDIT_INVALID,
        g(0, "/tolerance"),
        json!(next_up(1e-6)),
    ));
    cases.push(set(
        "tolerancia negativa mínima",
        AUDIT_INVALID,
        g(0, "/tolerance"),
        json!(tiny_negative()),
    ));
    cases.push(set(
        "frecuencia de referencia justo por encima de 1",
        AUDIT_INVALID,
        g(0, "/reference_frequency"),
        json!(next_up(1.0)),
    ));
    cases.push(set(
        "frecuencia de referencia negativa mínima",
        AUDIT_INVALID,
        g(0, "/reference_frequency"),
        json!(tiny_negative()),
    ));
    for field in [
        "original_sum",
        "original_delta",
        "effective_sum",
        "effective_delta",
    ] {
        cases.push(set(
            format!("{field} negativo mínimo"),
            AUDIT_INVALID,
            g(0, &format!("/{field}")),
            json!(tiny_negative()),
        ));
    }
    assert_rejected(&package, no_base, &cases);

    let mut accepted = vec![
        ok("posición 1", g(0, "/position"), json!("1")),
        ok("tolerancia 0", g(0, "/tolerance"), json!(0.0)),
        ok("tolerancia 1e-6", g(0, "/tolerance"), json!(1e-6)),
        ok(
            "frecuencia de referencia 0",
            g(0, "/reference_frequency"),
            json!(0.0),
        ),
        ok(
            "frecuencia de referencia 1",
            g(0, "/reference_frequency"),
            json!(1.0),
        ),
        ok(
            "decisión normalizada",
            g(0, "/decision"),
            json!("normalized"),
        ),
        ok_with(
            "dos grupos que solo difieren en la base de referencia, en orden",
            insert_group_with_ref(1, "C"),
        ),
    ];
    for field in [
        "original_sum",
        "original_delta",
        "effective_sum",
        "effective_delta",
    ] {
        accepted.push(ok(
            format!("{field} 0"),
            g(0, &format!("/{field}")),
            json!(0.0),
        ));
        accepted.push(ok(
            format!("{field} mayor que 1"),
            g(0, &format!("/{field}")),
            json!(1.5),
        ));
    }
    assert_accepted(&package, no_base, &accepted);
}

#[test]
fn the_normalization_audit_values_are_judged() {
    let package = statistics_package();
    let cases = [
        set(
            "un valor que no es un objeto",
            AUDIT_INVALID,
            g(0, "/values/0"),
            json!(5),
        ),
        add_key(
            "una clave de más en el valor",
            AUDIT_INVALID,
            g(0, "/values/0"),
        ),
        remove(
            "una clave de menos en el valor",
            AUDIT_INVALID,
            g(0, "/values/0/effective"),
        ),
        rename(
            "una clave desconocida en el valor",
            AUDIT_INVALID,
            g(0, "/values/0"),
            "original",
        ),
        set(
            "un alelo alterno fuera de ACGT",
            AUDIT_INVALID,
            g(0, "/values/0/alt"),
            json!("N"),
        ),
        set(
            "un alelo alterno de dos letras",
            AUDIT_INVALID,
            g(0, "/values/0/alt"),
            json!("GG"),
        ),
        set(
            "un alelo alterno numérico",
            AUDIT_INVALID,
            g(0, "/values/0/alt"),
            json!(5),
        ),
        set(
            "un alelo alterno igual a la referencia",
            AUDIT_INVALID,
            g(0, "/values/0/alt"),
            json!("A"),
        ),
        set(
            "una variante vacía",
            AUDIT_INVALID,
            g(0, "/values/0/variant_id"),
            json!(""),
        ),
        set(
            "una variante numérica",
            AUDIT_INVALID,
            g(0, "/values/0/variant_id"),
            json!(5),
        ),
        set(
            "una fuente numérica",
            AUDIT_INVALID,
            g(0, "/values/0/source_id"),
            json!(5),
        ),
        set(
            "una fuente nula",
            AUDIT_INVALID,
            g(0, "/values/0/source_id"),
            json!(null),
        ),
        set(
            "frecuencia original mayor que 1",
            AUDIT_INVALID,
            g(0, "/values/0/original"),
            json!(1.5),
        ),
        set(
            "frecuencia original negativa",
            AUDIT_INVALID,
            g(0, "/values/0/original"),
            json!(-0.1),
        ),
        set(
            "frecuencia original textual",
            AUDIT_INVALID,
            g(0, "/values/0/original"),
            json!("0.5"),
        ),
        set(
            "frecuencia efectiva mayor que 1",
            AUDIT_INVALID,
            g(0, "/values/0/effective"),
            json!(1.5),
        ),
        set(
            "frecuencia efectiva negativa",
            AUDIT_INVALID,
            g(0, "/values/0/effective"),
            json!(-0.1),
        ),
        set(
            "frecuencia efectiva textual",
            AUDIT_INVALID,
            g(0, "/values/0/effective"),
            json!("0.5"),
        ),
        case(
            "dos valores en orden de alelo descendente",
            AUDIT_INVALID,
            set_values(vec![
                audit_value("G", "VAR_X", "s"),
                audit_value("C", "VAR_X", "s"),
            ]),
        ),
        case(
            "dos valores en orden de variante descendente",
            AUDIT_INVALID,
            set_values(vec![
                audit_value("G", "VAR_B", "s"),
                audit_value("G", "VAR_A", "s"),
            ]),
        ),
        case(
            "dos valores en orden de fuente descendente",
            AUDIT_INVALID,
            set_values(vec![
                audit_value("G", "VAR_X", "b"),
                audit_value("G", "VAR_X", "a"),
            ]),
        ),
        case(
            "un valor repetido",
            AUDIT_INVALID,
            set_values(vec![
                audit_value("G", "VAR_X", "s"),
                audit_value("G", "VAR_X", "s"),
            ]),
        ),
    ];
    let mut cases = Vec::from(cases);
    for field in ["original", "effective"] {
        let pointer = g(0, &format!("/values/0/{field}"));
        cases.push(set(
            format!("frecuencia {field} justo por encima de 1"),
            AUDIT_INVALID,
            pointer.clone(),
            json!(next_up(1.0)),
        ));
        cases.push(set(
            format!("frecuencia {field} negativa mínima"),
            AUDIT_INVALID,
            pointer,
            json!(tiny_negative()),
        ));
    }
    assert_rejected(&package, no_base, &cases);

    let accepted = [
        ok(
            "frecuencias original y efectiva 0",
            g(0, "/values/0/original"),
            json!(0.0),
        ),
        ok(
            "frecuencia efectiva 1",
            g(0, "/values/0/effective"),
            json!(1.0),
        ),
        ok("una fuente vacía", g(0, "/values/0/source_id"), json!("")),
        ok_with(
            "dos valores en orden de alelo",
            set_values(vec![
                audit_value("C", "VAR_X", "s"),
                audit_value("G", "VAR_X", "s"),
            ]),
        ),
        ok_with(
            "dos valores en orden de variante",
            set_values(vec![
                audit_value("G", "VAR_A", "s"),
                audit_value("G", "VAR_B", "s"),
            ]),
        ),
        ok_with(
            "dos valores en orden de fuente",
            set_values(vec![
                audit_value("G", "VAR_X", "a"),
                audit_value("G", "VAR_X", "b"),
            ]),
        ),
    ];
    assert_accepted(&package, no_base, &accepted);
}

fn inserting() -> Package {
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
    let config = r#"config_kind = "general"
schema_version = "4"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

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
"#;
    publish_files(&[("general.toml", config), ("reference.fa", &reference)])
}

// ---------------------------------------------------------------------------
// El juez físico de la matriz de SNV
// ---------------------------------------------------------------------------

/// La misma configuración sin poblaciones iniciales: solo la matriz.
fn matrix_only_toml() -> String {
    let end = STATISTICS_TOML
        .find("[[initial_populations]]")
        .expect("poblaciones iniciales");
    STATISTICS_TOML[..end].to_owned()
}

/// Lo que dice el juez completo con estos dos documentos, o `ACCEPTED`.
fn full(package: &Package, effective: &Value, manifest: &Value) -> String {
    match validate_published_resources(
        &package.destination,
        &package.destination.join("effective-config.json"),
        &serde_json::to_vec_pretty(effective).expect("JSON"),
        &package.destination.join("manifest.json"),
        &serde_json::to_vec_pretty(manifest).expect("JSON"),
        None,
    ) {
        Ok(()) => "ACCEPTED".to_owned(),
        Err(error) => error.context().get("detail").cloned().unwrap_or_default(),
    }
}

/// Presenta la matriz como **reutilizada** en los dos documentos: modo
/// `reuse`, sin artefacto de metadatos, y fuera del catálogo y de las
/// categorías de salida la pareja metadatos/TSV, que una reutilización no
/// publica.
fn present_as_reused(effective: &mut Value, manifest: &mut Value) {
    for document in [&mut *effective, &mut *manifest] {
        document["population_snvs"]["mode"] = json!("reuse");
        document["population_snvs"]["metadata_artifact_id"] = json!(null);
    }
    let is_pair = |id: &str| id == "population-snvs-metadata" || id == "population-snvs-tsv";
    for field in ["minimum", "requested"] {
        manifest["output_categories"][field]
            .as_array_mut()
            .expect("categoría")
            .retain(|id| !is_pair(id.as_str().expect("id")));
    }
    manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
        .retain(|artifact| !is_pair(artifact["id"].as_str().expect("id")));
}

/// **El juez físico lee la matriz que el TOML copiado declara**, no la que los
/// documentos dicen: una efectiva que presenta como reutilizada una matriz que
/// el TOML genera, o que el TOML no declara, se rechaza antes de leer nada.
///
/// El segundo caso rehace por completo lo que lo acompaña: el TOML copiado
/// sin su sección de matriz, el identificador de ejecución del manifiesto
/// (la derivación es pública) y la presentación como reutilizada.
#[test]
fn the_physical_judge_reads_the_matrix_source_the_toml_declares() {
    let package = publish(&matrix_only_toml());
    let effective = package.json("effective-config.json");
    let manifest = package.json("manifest.json");
    assert_eq!(
        full(&package, &effective, &manifest),
        "ACCEPTED",
        "control positivo: el paquete intacto se acepta"
    );

    // El TOML genera la matriz; la efectiva y el manifiesto dicen que se
    // reutiliza.
    let (mut reused_effective, mut reused_manifest) = (effective.clone(), manifest.clone());
    present_as_reused(&mut reused_effective, &mut reused_manifest);
    assert_eq!(
        full(&package, &reused_effective, &reused_manifest),
        "modo SNV efectivo no coincide con reutilización"
    );

    // El TOML copiado no declara ninguna matriz.
    let toml_path = package.destination.join(
        manifest["artifacts"]
            .as_array()
            .expect("catálogo")
            .iter()
            .find(|artifact| artifact["id"] == "config-general")
            .and_then(|artifact| artifact["path"].as_str())
            .expect("TOML copiado"),
    );
    let original = fs::read_to_string(&toml_path).expect("TOML copiado");
    let start = original
        .find("[population_snvs]")
        .expect("sección de matriz");
    let end = original.find("[reference]").expect("referencia");
    let without_matrix = format!("{}{}", &original[..start], &original[end..]);
    assert_ne!(without_matrix, original, "la alteración no alteró nada");
    fs::write(
        &toml_path,
        detection_section::with_detection_section(&without_matrix),
    )
    .expect("TOML alterado");
    let digest = configuration_digest(without_matrix.as_bytes()).expect("huella del TOML");
    let seed = effective["seed"]["root"].as_str().expect("semilla");
    reused_manifest["execution_id"] = json!(
        derive_execution_id(RootSeed::parse("semilla", seed).expect("semilla"), &digest)
            .expect("execution_id")
            .as_str()
    );
    let result = full(&package, &reused_effective, &reused_manifest);
    fs::write(
        &toml_path,
        detection_section::with_detection_section(&original),
    )
    .expect("TOML restaurado");
    assert_eq!(result, "reutilización SNV sin fuente física");
}

/// **El juez físico relee el par publicado (metadatos y TSV) y ata a él el
/// enlace de los documentos.** Dos rechazos que antes solo disparaban otras
/// pruebas, sin afirmar su texto: el par que no supera su propia validación
/// (un TSV que ya no es el que los metadatos declaran, con el catálogo
/// rehecho para que la identidad del fichero cuadre) y el enlace de los
/// documentos (auditoría de normalización) que no coincide con el par.
#[test]
fn the_physical_judge_ties_the_documents_to_the_published_pair() {
    let package = publish(&matrix_only_toml());
    let effective = package.json("effective-config.json");
    let manifest = package.json("manifest.json");

    // El enlace: la auditoría dice otra cosa que el par, en los dos documentos.
    for (label, pointer, value) in [
        (
            "otra decisión de normalización",
            g(0, "/decision"),
            json!("normalized"),
        ),
        (
            "un origen declarado distinto del que dicen los metadatos",
            "/population_snvs/mode".to_owned(),
            json!("supplied"),
        ),
    ] {
        let (mut altered_effective, mut altered_manifest) = (effective.clone(), manifest.clone());
        for document in [&mut altered_effective, &mut altered_manifest] {
            replace_at(document, &pointer, value.clone());
        }
        assert_eq!(
            full(&package, &altered_effective, &altered_manifest),
            "enlace population_snvs no coincide con el par físico",
            "«{label}»"
        );
    }

    // El par: un TSV con otro contenido y su identidad rehecha en el catálogo.
    let tsv_path = package.destination.join("population_snvs.tsv");
    let original = fs::read_to_string(&tsv_path).expect("TSV");
    let altered = original.replacen("0.6267319433564043", "0.6267319433564044", 1);
    assert_ne!(altered, original, "la alteración no alteró nada");
    fs::write(&tsv_path, &altered).expect("TSV alterado");
    let identity = identify_file(&tsv_path, PathBuf::from("population_snvs.tsv")).expect("huella");
    let mut altered_manifest = manifest.clone();
    for artifact in altered_manifest["artifacts"]
        .as_array_mut()
        .expect("catálogo")
    {
        if artifact["id"] == "population-snvs-tsv" {
            artifact["checksum"]["digest"] = json!(identity.digest());
            artifact["size_bytes"] = json!(identity.size_bytes());
        }
    }
    let result = full(&package, &effective, &altered_manifest);
    fs::write(&tsv_path, original).expect("TSV restaurado");
    assert_eq!(
        result,
        "par SNV físico inválido: AuthoritativeFasta { detail: \"SHA-256 del TSV reutilizado no coincide con metadata\" }",
        "el par con otro TSV"
    );
}

/// **Intento ejecutado de alcanzar las guardas inalcanzables del juez físico
/// de la matriz y de la raíz de las estadísticas.** Cada caso intenta llegar a
/// una guarda anotada en `src/results.rs` como sombreada por una comprobación
/// anterior y comprueba que el rechazo viene de **esa otra** comprobación, con
/// su texto. **No discrimina las guardas anotadas**: si se borrara una de ellas
/// la prueba seguiría pasando, porque otra comprobación ya rechaza antes.
#[test]
fn unreachable_physical_snv_guards_are_shadowed_by_other_checks() {
    type Alter = fn(&mut Value, &mut Value);
    let shadowed: [(&str, Alter, &str); 8] = [
        (
            "reutilizada y sin artefacto FASTA en la referencia",
            |effective, manifest| {
                present_as_reused(effective, manifest);
                effective["reference"]
                    .as_object_mut()
                    .expect("referencia")
                    .remove("fasta_artifact_id");
            },
            "campos ausentes o desconocidos en reference",
        ),
        (
            "reutilizada y con un artefacto FASTA que el catálogo no tiene",
            |effective, manifest| {
                present_as_reused(effective, manifest);
                effective["reference"]["fasta_artifact_id"] = json!("nope");
            },
            "el enlace nope no existe",
        ),
        (
            "publicada con los metadatos sin identificador",
            |effective, manifest| {
                for document in [effective, manifest] {
                    document["population_snvs"]["metadata_artifact_id"] = json!(null);
                }
            },
            "population_snvs.metadata_artifact_id inválido",
        ),
        (
            "publicada con metadatos que el catálogo no tiene",
            |effective, manifest| {
                for document in [effective, manifest] {
                    document["population_snvs"]["metadata_artifact_id"] = json!("nope");
                }
            },
            "el enlace nope no existe",
        ),
        (
            "publicada sin el TSV en el catálogo",
            |_, manifest| {
                manifest["artifacts"]
                    .as_array_mut()
                    .expect("catálogo")
                    .retain(|artifact| artifact["id"] != "population-snvs-tsv");
            },
            "el enlace population-snvs-tsv no existe",
        ),
        (
            "sin artefacto FASTA en la referencia",
            |effective, _| {
                effective["reference"]
                    .as_object_mut()
                    .expect("referencia")
                    .remove("fasta_artifact_id");
            },
            "campos ausentes o desconocidos en reference",
        ),
        (
            "con un artefacto FASTA que el catálogo no tiene",
            |effective, _| effective["reference"]["fasta_artifact_id"] = json!("nope"),
            "el enlace nope no existe",
        ),
        (
            "sin tolerancia de normalización",
            |effective, _| {
                effective["normalization"]
                    .as_object_mut()
                    .expect("normalización")
                    .remove("tolerance");
            },
            "campos ausentes o desconocidos en normalization",
        ),
    ];
    let package = publish(&matrix_only_toml());
    for (label, alter, expected) in shadowed {
        let (mut effective, mut manifest) = (
            package.json("effective-config.json"),
            package.json("manifest.json"),
        );
        alter(&mut effective, &mut manifest);
        assert_eq!(full(&package, &effective, &manifest), expected, "«{label}»");
    }

    // La misma pareja con la referencia **final por inserción**: la referencia
    // que el juez físico mira es la insertada, y su artefacto lo ata antes el
    // juez de la referencia final.
    let package = inserting();
    let inserted: [(&str, Alter, &str); 4] = [
        (
            "reutilizada y sin artefacto FASTA en la referencia final",
            |effective, manifest| {
                present_as_reused(effective, manifest);
                effective["inserted_reference"]
                    .as_object_mut()
                    .expect("insertada")
                    .remove("fasta_artifact_id");
            },
            "campos ausentes o desconocidos en inserted_reference",
        ),
        (
            "reutilizada y con un artefacto FASTA final que el catálogo no tiene",
            |effective, manifest| {
                present_as_reused(effective, manifest);
                effective["inserted_reference"]["fasta_artifact_id"] = json!("nope");
            },
            "inserted_reference no coincide con la referencia declarada o con reference_insertion del manifiesto",
        ),
        (
            "sin artefacto FASTA en la referencia final",
            |effective, _| {
                effective["inserted_reference"]
                    .as_object_mut()
                    .expect("insertada")
                    .remove("fasta_artifact_id");
            },
            "campos ausentes o desconocidos en inserted_reference",
        ),
        (
            "con un artefacto FASTA final que el catálogo no tiene",
            |effective, _| effective["inserted_reference"]["fasta_artifact_id"] = json!("nope"),
            "inserted_reference no coincide con la referencia declarada o con reference_insertion del manifiesto",
        ),
    ];
    for (label, alter, expected) in inserted {
        let (mut effective, mut manifest) = (
            package.json("effective-config.json"),
            package.json("manifest.json"),
        );
        alter(&mut effective, &mut manifest);
        assert_eq!(full(&package, &effective, &manifest), expected, "«{label}»");
    }

    // El identificador de una población estadística lo comprueba antes la
    // validación de la propia población.
    let package = statistics_package();
    let cases = [
        remove(
            "una población estadística sin identificador",
            POPULATION_INVALID,
            stats("/populations/0/population_id"),
        ),
        set(
            "una población estadística con identificador numérico",
            POPULATION_INVALID,
            stats("/populations/0/population_id"),
            json!(5),
        ),
    ];
    // Un identificador que falta deja a la población con una clave de menos.
    assert_rejected(&package, no_base, &cases);
}

/// **Fija la premisa de las mitades `is_finite` que se dan por equivalentes.**
/// Un decimal que no cabe en un `f64` finito (`1e999`) se lee con `as_f64` como
/// `None`, de modo que ni una frecuencia, ni una probabilidad, ni una
/// tolerancia, ni un coeficiente llegan a la comprobación `is_finite` de
/// `validate_initial_population_statistics_locus`, de `..._model` y de
/// `validate_serialized_normalization_audit`: los rechaza antes la lectura del
/// número, con su propio motivo. **No discrimina esas mitades**: si se
/// borraran la prueba seguiría pasando.
#[test]
fn non_finite_decimals_never_reach_the_finite_halves() {
    const MARK: &str = "__OVERFLOW__";
    let package = statistics_package();
    let effective = package.json("effective-config.json");
    let manifest = package.json("manifest.json");
    let cases = [
        (
            l0("/alleles/0/requested_frequency"),
            "frecuencia objetivo inválida",
        ),
        (
            l0("/alleles/0/realized_frequency"),
            "frecuencia realizada inválida",
        ),
        (
            l0("/genotypes/0/expected_probability"),
            "probabilidad genotípica inválida",
        ),
        (
            l0("/genotypes/0/expected_count"),
            "recuento genotípico esperado inválido",
        ),
        (
            stats("/populations/1/genotype_model/inbreeding_coefficient"),
            MODEL_INVALID,
        ),
        (g(0, "/tolerance"), AUDIT_INVALID),
        (g(0, "/reference_frequency"), AUDIT_INVALID),
        (g(0, "/original_sum"), AUDIT_INVALID),
        (g(0, "/values/0/original"), AUDIT_INVALID),
    ];
    for (pointer, expected) in cases {
        let mut altered = manifest.clone();
        replace_at(&mut altered, &pointer, json!(MARK));
        let bytes = serde_json::to_string_pretty(&altered)
            .expect("JSON")
            .replace(&format!("\"{MARK}\""), "1e999");
        assert!(
            bytes.contains("1e999"),
            "la alteración no alteró nada: {pointer}"
        );
        let result = validate_effective_config_against_manifest_bytes(
            &package.destination.join("effective-config.json"),
            &serde_json::to_vec_pretty(&effective).expect("JSON"),
            &package.destination.join("manifest.json"),
            bytes.as_bytes(),
        );
        assert_eq!(detail(result), expected, "{pointer}");
    }
}
