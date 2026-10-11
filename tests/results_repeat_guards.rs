//! Cada rechazo del juez de los recursos publicados en lo que toca a los
//! microsatélites, el tartamudeo y la PCR, con su motivo exacto.
//!
//! DT-001, tanda 6 (primera parte de `src/results.rs`). El juez tiene 779
//! comprobaciones y 618 no las disparaba ninguna prueba; esta tanda cubre las
//! de diez funciones: la sección de loci del manifiesto general, el abanico de
//! alelos generado por deslizamiento (el recorrido que rehace el juez y su
//! comparación), el bloque del generador, y las secciones de tartamudeo, de
//! dependencia con la longitud del tracto, de tartamudeo clonal y de
//! amplificación de la configuración efectiva.
//!
//! # Cómo se parte de un paquete real
//!
//! Cada caso parte de un paquete **publicado de verdad** por el programa y lo
//! altera **en un solo punto**; el motivo exacto de cada rechazo es la
//! discriminación (varias guardas comparten la causa tipada `GPS-FMT-006`). Se
//! usan dos escenarios:
//!
//! - **`mixed`**: referencia suministrada con dos loci, uno **declarado a mano**
//!   (con un alelo interrumpido y otro puro) y otro **generado** por
//!   deslizamiento, cinco individuos diploides. Lo juzga el lector del
//!   manifiesto (`validate_general_manifest_bytes`) y el juez completo
//!   (`validate_published_resources`) para el abanico generado.
//! - **`amplified`**: el escenario de la evidencia de DTC-099 (dos loci, PCR y
//!   tartamudeo de base y clonal) y su variante sin loci, que amplifica sin
//!   microsatélites. Lo juzga `validate_effective_config_against_manifest_bytes`
//!   con la configuración efectiva alterada y el manifiesto intacto: las
//!   comprobaciones de estructura del evento técnico corren antes que todo
//!   lo que ata la efectiva al TOML.
//!
//! Cada tabla lleva su control positivo (el documento intacto se acepta: sin
//! él, un juez que rechazara todo satisfaría la tabla entera), comprueba que
//! **la alteración alteró**, y los casos «aceptado» fijan el lado aceptado de
//! los límites y de las listas cerradas.
//!
//! # DT-051: los loci publicados, atados al TOML y a la verdad
//!
//! La última sección comprueba que el juez completo ata los loci publicados
//! —los de la configuración efectiva, sus mínimos de admisión y la sección de
//! loci del manifiesto general— a lo que el TOML copiado produce, y las copias
//! realizadas de cada alelo a lo que la verdad lleva. Son las alteraciones que
//! la revisión adversarial de la tanda 6 reprodujo, aceptadas todas hasta el
//! 10/10/2026, el trasplante de un abanico con sus parámetros y el lado
//! aceptado de dos alelos que la verdad no distingue.
//!
//! # DT-052 y DT-043: la contabilidad de tracto y la amplificación
//!
//! Las dos últimas secciones comprueban que el juez completo recuenta al
//! publicar las bases de tracto y ata su cota, en el manifiesto y en la
//! configuración efectiva, al TOML copiado (DT-052), y que cada parámetro de la
//! amplificación del manifiesto de evento es el de la configuración efectiva
//! (DT-043), con el lado aceptado sobre paquetes publicados de verdad.
//!
//! # Lo que este fichero no es
//!
//! Las guardas que ningún paquete publicado y manipulado puede alcanzar están
//! anotadas **junto a la guarda** en `src/results.rs`, con el porqué real.
//! Las pruebas de este fichero que solo «ejecutan el intento» sobre ellas lo
//! dicen en su docstring: **no discriminan** esas guardas.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::{
    validate_effective_config_against_manifest_bytes, validate_general_manifest_bytes,
    validate_published_resources,
};
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);
const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

const AMPLIFIED_TOML: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/general.toml"
);
const AMPLIFIED_FASTA: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/reference.fa"
);
const AMPLIFIED_SNVS: &str = include_str!(
    "../evidence/dt-009-multiplicador-congelacion/inputs/amplificado-loci/population-snvs.tsv"
);

// ---------------------------------------------------------------------------
// Paquetes publicados
// ---------------------------------------------------------------------------

/// Un paquete publicado, que se borra al soltarlo.
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
        "gps-results-repeat-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("directorio");
    for (name, content) in files {
        fs::write(root.join(name), content).expect("fichero de entrada");
    }
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&root.join("general.toml")).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Package { root, destination }
}

/// Doscientas bases de flanco, cincuenta de `(GT)25` y cien de `(CA)50`, y
/// cien de flanco: los dos loci del escenario `mixed`.
fn mixed_reference() -> String {
    format!(
        ">chr1\n{}{}{}{}\n",
        "G".repeat(50),
        "GT".repeat(25),
        "CA".repeat(50),
        "G".repeat(100)
    )
}

/// El primer alelo del locus declarado: cincuenta bases `GTAA` + `GT` × 23, con
/// una interrupción `AA` tras la primera copia.
const INTERRUPTED_ALLELE: &str = "GTAAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT";

/// Escenario `mixed`: un locus declarado (alelo interrumpido y alelo puro) y un
/// locus generado por deslizamiento, con `individuals` individuos diploides.
fn publish_mixed(individuals: usize) -> Package {
    publish_mixed_with_ploidy(individuals, 2)
}

/// Lo mismo con la ploidía de la población que se pida (los linajes del
/// generador son individuos por ploidía).
fn publish_mixed_with_ploidy(individuals: usize, ploidy: u32) -> Package {
    publish_mixed_with(individuals, ploidy, "3.0")
}

/// Lo mismo con la media de pasos del generador que se pida (DT-051: el
/// trasplante de un abanico con otros parámetros).
fn publish_mixed_with(individuals: usize, ploidy: u32, mean_events: &str) -> Package {
    publish_mixed_with_tail(individuals, ploidy, mean_events, "")
}

/// Lo mismo, con un VCF de texto (DT-052: la parte VCF de la contabilidad de
/// tracto).
fn publish_mixed_with_vcf(individuals: usize) -> Package {
    publish_mixed_with_tail(individuals, 2, "3.0", "\n[vcf_output]\nmode = \"text\"\n")
}

/// El escenario `mixed` con lo que se pida añadido al final del TOML.
fn publish_mixed_with_tail(
    individuals: usize,
    ploidy: u32,
    mean_events: &str,
    tail: &str,
) -> Package {
    let ids = (1..=individuals)
        .map(|index| format!("\"I{index:03}\""))
        .collect::<Vec<_>>()
        .join(", ");
    let reference = mixed_reference();
    let length = reference.trim_end().lines().nth(1).expect("bases").len();
    let config = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
outside_weight = 0.0
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.regions]]
contig = "chr1"
start = 1
end = 50
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

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {length}

[[reference.repeat_loci]]
id = "STR-LEFT"
sequence = "chr1"
start = 51
end = 100
motif = "GT"

[[reference.repeat_loci.alleles]]
sequence = "{INTERRUPTED_ALLELE}"
frequency = 0.6

[[reference.repeat_loci.alleles]]
units = 24
frequency = 0.4

[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = {mean_events}
expansion_bias = 0.5

[[initial_populations]]
population_id = "P1"
individual_count = {individuals}
individual_ids = [{ids}]
ploidy = {ploidy}
genotype_model = "hwe"
{tail}"#
    );
    publish(&[("general.toml", &config), ("reference.fa", &reference)])
}

/// Escenario `edges`: cuatro loci generados cuyo recorrido **llega a los
/// bordes** del modelo, para que el juez del abanico tenga que reproducir lo
/// que ocurre al detenerse un linaje.
///
/// - `EDGE-CONTRACT`: `(CA)3`, sesgo 0 (solo contrae): baja a dos unidades, a
///   una, y se detiene.
/// - `EDGE-CAP`: `(CA)10` con tope 10: no puede expandirse, sí contraerse.
/// - `EDGE-INT`: `CACACA T CACACA`, un tracto interrumpido: sus alelos
///   ordenados por longitud no coinciden con el orden de sus bases.
/// - `EDGE-UP`: `(CA)4` con tope 5 y sesgo 1 (solo expande): sube a cinco y se
///   detiene.
fn publish_edges() -> Package {
    let reference = format!(
        ">chr1\n{}{}{}{}{}{}{}{}{}\n",
        "G".repeat(20),
        "CA".repeat(3),
        "G".repeat(20),
        "CA".repeat(10),
        "G".repeat(20),
        "CACACATCACACA",
        "G".repeat(20),
        "CA".repeat(4),
        "G".repeat(20),
    );
    let locus = |id: &str, start: u32, end: u32, cap: u32, bias: &str| {
        format!(
            r#"
[[reference.repeat_loci]]
id = "{id}"
sequence = "chr1"
start = {start}
end = {end}
motif = "CA"
max_units = {cap}

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 4.0
expansion_bias = {bias}
"#
        )
    };
    let loci = [
        locus("EDGE-CONTRACT", 21, 26, 90, "0.0"),
        locus("EDGE-CAP", 47, 66, 10, "0.5"),
        locus("EDGE-INT", 87, 99, 90, "0.5"),
        locus("EDGE-UP", 120, 127, 5, "1.0"),
    ]
    .concat();
    let config = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
outside_weight = 0.0
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.regions]]
contig = "chr1"
start = 1
end = 20
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

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 147
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 12
individual_ids = ["I01","I02","I03","I04","I05","I06","I07","I08","I09","I10","I11","I12"]
ploidy = 2
genotype_model = "hwe"
"#
    );
    publish(&[("general.toml", &config), ("reference.fa", &reference)])
}

/// Escenario `inserting`: un tracto existente generado por deslizamiento
/// (`(CA)20`) y otro locus que **inserta** un tracto `(CAG)10` en la
/// referencia: la referencia sobre la que se simula es la final, con otro
/// identificador, y el ámbito de azar del abanico se siembra con él.
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
    let reference = format!(">chr1\n{bases}\n");
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
id = "STR-CA"
sequence = "chr1"
start = 2001
end = 2040
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5

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

/// El escenario de la evidencia de DTC-099: amplificado, con dos loci
/// (homopolímero y dinucleótido), tartamudeo de base y tartamudeo clonal.
/// `clonal` se inserta en `[technical_event.pcr.stutter]` si no está vacío.
fn publish_amplified(clonal: &str) -> Package {
    let toml = if clonal.is_empty() {
        AMPLIFIED_TOML.to_owned()
    } else {
        let anchor = "\n[technical_event.stutter]\n";
        assert_eq!(AMPLIFIED_TOML.matches(anchor).count(), 1, "fixture");
        AMPLIFIED_TOML.replacen(
            anchor,
            &format!("\n[technical_event.pcr.stutter]\n{clonal}{anchor}"),
            1,
        )
    };
    publish(&[
        ("general.toml", &toml),
        ("reference.fa", AMPLIFIED_FASTA),
        ("population-snvs.tsv", AMPLIFIED_SNVS),
    ])
}

/// Amplifica **sin** microsatélites: sin loci no hay sección de tartamudeo
/// clonal en la amplificación, y su catálogo es el v5. `fraction` es la
/// fracción de duplicados declarada.
fn publish_amplified_without_loci(fraction: &str) -> Package {
    let toml = format!(
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

[[initial_populations]]
population_id = "P1"
individual_count = 2
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A", "sample-B"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A", "sample-B"]
target_coverage_wgs = 20
library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = {fraction}
cycles = 10
efficiency = 0.9

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"

[[technical_event.libraries]]
individual_id = "sample-B"
i7 = "TTGGCCAA"
i5 = "GGAATTCC"
"#
    );
    publish(&[
        ("general.toml", &toml),
        ("reference.fa", AMPLIFIED_FASTA),
        ("population-snvs.tsv", AMPLIFIED_SNVS),
    ])
}

fn detail(result: Result<(), Diagnostic>) -> String {
    result
        .expect_err("el juez tenía que rechazar")
        .context()
        .get("detail")
        .cloned()
        .unwrap_or_default()
}

/// Un caso: una etiqueta, el texto **exacto** del rechazo y la alteración.
struct Case {
    label: String,
    expected: &'static str,
    mutate: Box<dyn Fn(&mut Value)>,
}

fn case(
    label: impl Into<String>,
    expected: &'static str,
    mutate: impl Fn(&mut Value) + 'static,
) -> Case {
    Case {
        label: label.into(),
        expected,
        mutate: Box::new(mutate),
    }
}

/// Un caso que el juez **tiene que aceptar**: el lado aceptado de un límite o
/// de una lista cerrada.
struct Accepted {
    label: String,
    mutate: Box<dyn Fn(&mut Value)>,
}

fn accepted(label: impl Into<String>, mutate: impl Fn(&mut Value) + 'static) -> Accepted {
    Accepted {
        label: label.into(),
        mutate: Box::new(mutate),
    }
}

/// Aplica la alteración a una copia del documento y comprueba que **alteró**.
fn altered(original: &Value, label: &str, mutate: &dyn Fn(&mut Value)) -> Value {
    let mut document = original.clone();
    mutate(&mut document);
    assert_ne!(
        &document, original,
        "«{label}»: la alteración no alteró nada"
    );
    document
}

/// Un número JSON que no cabe en un `f64` finito (`1e999`).
///
/// Con `serde_json` y su característica `arbitrary_precision`, el número se
/// conserva como texto y `as_f64` devuelve `None` para él: **ningún documento
/// publicado puede aportar un decimal no finito**, de modo que las mitades
/// `!x.is_finite()` de las comprobaciones de decimales no deciden nunca; lo
/// que rechaza este valor es la lectura del decimal, con su motivo.
fn overflowing_decimal() -> Value {
    serde_json::from_str("1e999").expect("un número JSON válido")
}

/// El primer `f64` mayor que `value` (positivo).
fn next_up(value: f64) -> f64 {
    assert!(value > 0.0 && value.is_finite());
    f64::from_bits(value.to_bits() + 1)
}

// ---------------------------------------------------------------------------
// La sección de loci del manifiesto general
// ---------------------------------------------------------------------------

fn judge_manifest(package: &Package, document: &Value) -> Result<(), Diagnostic> {
    validate_general_manifest_bytes(
        &package.manifest_path(),
        &serde_json::to_vec(document).expect("JSON"),
    )
}

/// Control positivo, y después cada caso: la alteración altera, y el rechazo
/// es el exacto.
fn assert_manifest_cases(package: &Package, cases: &[Case]) {
    let original = package.json("manifest.json");
    judge_manifest(package, &original).expect("control positivo: el manifiesto intacto se acepta");
    for case in cases {
        let document = altered(&original, &case.label, case.mutate.as_ref());
        assert_eq!(
            detail(judge_manifest(package, &document)),
            case.expected,
            "«{}» se rechazó por otra razón",
            case.label
        );
    }
}

fn assert_manifest_accepted(package: &Package, cases: &[Accepted]) {
    let original = package.json("manifest.json");
    for case in cases {
        let document = altered(&original, &case.label, case.mutate.as_ref());
        judge_manifest(package, &document)
            .unwrap_or_else(|error| panic!("«{}» se tenía que aceptar: {error:?}", case.label));
    }
}

fn section(document: &mut Value) -> &mut Value {
    &mut document["repeat_loci"]
}

/// El locus declarado a mano (`STR-LEFT`).
fn declared(document: &mut Value) -> &mut Value {
    &mut document["repeat_loci"]["loci"][0]
}

/// El locus generado por deslizamiento (`STR-RIGHT`).
fn generated(document: &mut Value) -> &mut Value {
    &mut document["repeat_loci"]["loci"][1]
}

/// El alelo interrumpido del locus declarado.
fn allele0(document: &mut Value) -> &mut Value {
    &mut document["repeat_loci"]["loci"][0]["alleles"][0]
}

/// El alelo puro del locus declarado.
fn allele1(document: &mut Value) -> &mut Value {
    &mut document["repeat_loci"]["loci"][0]["alleles"][1]
}

const ADMISSION_NOT_OBJECT: &str = "unos mínimos de admisión no son un objeto";
const ADMISSION_FIELDS: &str = "unos mínimos de admisión no tienen exactamente sus campos";
const PURITY_NOT_DECIMAL: &str = "la pureza mínima no es un decimal";
const PURITY_RANGE: &str = "la pureza mínima no es finita ni pertenece a [0,1]";
const RUN_NOT_POSITIVE: &str = "la tirada mínima no es un entero positivo";

/// **Los mínimos de admisión se juzgan con los mismos motivos desde sus dos
/// llamantes**: los de omisión de la sección y los de cada locus.
///
/// Un ayudante compartido se prueba por llamante: si un día los dos dejaran
/// de usarlo, una tabla sobre uno solo no se enteraría.
#[test]
fn the_admission_minimums_are_judged_from_both_callers() {
    let package = publish_mixed(5);
    type Pointer = fn(&mut Value) -> &mut Value;
    let callers: [(&str, Pointer); 2] = [
        (
            "de omisión",
            |d| &mut d["repeat_loci"]["admission_defaults"],
        ),
        (
            "del locus",
            |d| &mut d["repeat_loci"]["loci"][0]["admission"],
        ),
    ];
    let mut cases = Vec::new();
    let mut accept = Vec::new();
    for (who, at) in callers {
        cases.push(case(
            format!("mínimos {who} que no son objeto"),
            ADMISSION_NOT_OBJECT,
            move |d| *at(d) = json!(5),
        ));
        cases.push(case(
            format!("mínimos {who} con un campo de más"),
            ADMISSION_FIELDS,
            move |d| at(d)["sobra"] = json!(1),
        ));
        cases.push(case(
            format!("mínimos {who} con un campo de menos"),
            ADMISSION_FIELDS,
            move |d| {
                at(d).as_object_mut().expect("mínimos").remove("min_purity");
            },
        ));
        cases.push(case(
            format!("mínimos {who} con un campo cambiado por otro"),
            ADMISSION_FIELDS,
            move |d| {
                let object = at(d).as_object_mut().expect("mínimos");
                let value = object.remove("min_purity").expect("pureza");
                object.insert("pureza".to_owned(), value);
            },
        ));
        cases.push(case(
            format!("pureza {who} que no es un decimal"),
            PURITY_NOT_DECIMAL,
            move |d| at(d)["min_purity"] = json!("x"),
        ));
        cases.push(case(
            format!("pureza {who} que no cabe en un decimal finito"),
            PURITY_NOT_DECIMAL,
            move |d| at(d)["min_purity"] = overflowing_decimal(),
        ));
        for (what, value) in [("por encima de uno", 1.5), ("negativa", -0.1)] {
            cases.push(case(
                format!("pureza {who} {what}"),
                PURITY_RANGE,
                move |d| at(d)["min_purity"] = json!(value),
            ));
        }
        cases.push(case(
            format!("tirada mínima {who} igual a cero"),
            RUN_NOT_POSITIVE,
            move |d| at(d)["min_longest_run_units"] = json!(0),
        ));
        cases.push(case(
            format!("tirada mínima {who} que no es un entero"),
            RUN_NOT_POSITIVE,
            move |d| at(d)["min_longest_run_units"] = json!("x"),
        ));
        // El lado aceptado: los dos extremos cerrados de la pureza y la
        // tirada mínima de una unidad.
        accept.push(accepted(format!("pureza {who} igual a 0"), move |d| {
            at(d)["min_purity"] = json!(0.0);
        }));
        accept.push(accepted(format!("pureza {who} igual a 1"), move |d| {
            at(d)["min_purity"] = json!(1.0);
        }));
        accept.push(accepted(format!("tirada mínima {who} de 1"), move |d| {
            at(d)["min_longest_run_units"] = json!(1);
        }));
    }
    assert_manifest_cases(&package, &cases);
    assert_manifest_accepted(&package, &accept);
}

/// **La contabilidad de bases de tracto: forma y aritmética.**
///
/// Las seis cifras tienen que existir y ser enteros; cada total, la suma de
/// sus dos términos —con aritmética comprobada: un desborde **rechaza**, no
/// da la vuelta—; y lo realizado no pasa de su cota. El lado aceptado es la
/// igualdad (realizado igual a la cota).
#[test]
fn the_tract_output_accounting_is_judged_for_shape_and_arithmetic() {
    let package = publish_mixed(5);
    const NOT_OBJECT: &str = "la contabilidad de tracto no es un objeto";
    const FIELDS: &str = "la contabilidad de tracto no tiene exactamente sus campos";
    const NOT_INTEGER: &str = "una cifra de la contabilidad de tracto no es un entero";
    const BOUND_SUM: &str = "la cota total no es la suma de sus dos términos";
    const REALIZED_SUM: &str = "el total realizado no es la suma de sus dos términos";
    fn accounting(d: &mut Value) -> &mut Value {
        &mut d["repeat_loci"]["tract_output"]
    }
    fn set(d: &mut Value, values: [u64; 6]) {
        let names = [
            "truth_tract_bases_bound",
            "vcf_tract_bases_bound",
            "total_tract_bases_bound",
            "realized_truth_tract_bases",
            "realized_vcf_tract_bases",
            "realized_total_tract_bases",
        ];
        for (name, value) in names.iter().zip(values) {
            accounting(d)[*name] = json!(value);
        }
    }
    let mut cases = vec![
        case("contabilidad que no es objeto", NOT_OBJECT, |d| {
            *accounting(d) = json!(5);
        }),
        case("contabilidad con un campo de más", FIELDS, |d| {
            accounting(d)["sobra"] = json!(1);
        }),
        case("contabilidad con un campo de menos", FIELDS, |d| {
            accounting(d)
                .as_object_mut()
                .expect("contabilidad")
                .remove("realized_vcf_tract_bases");
        }),
        case(
            "lo realizado de la verdad por encima de su cota",
            "las bases de tracto realizadas superan su cota",
            |d| set(d, [10, 20, 30, 11, 0, 11]),
        ),
        case(
            "lo realizado del VCF por encima de su cota",
            "las bases de tracto realizadas superan su cota",
            |d| set(d, [10, 20, 30, 0, 21, 21]),
        ),
        case("contabilidad con un campo cambiado por otro", FIELDS, |d| {
            let object = accounting(d).as_object_mut().expect("contabilidad");
            let value = object.remove("realized_vcf_tract_bases").expect("cifra");
            object.insert("realizado_vcf".to_owned(), value);
        }),
        case("cota total que no es la suma", BOUND_SUM, |d| {
            set(d, [10, 20, 31, 10, 20, 30]);
        }),
        case("cota total cuya suma desborda", BOUND_SUM, |d| {
            // La suma desborda: con aritmética que diera la vuelta, 0 sería
            // la cota «correcta».
            set(d, [u64::MAX, 1, 0, 0, 0, 0]);
        }),
        case("total realizado que no es la suma", REALIZED_SUM, |d| {
            set(d, [10, 20, 30, 10, 20, 31]);
        }),
        case("total realizado cuya suma desborda", REALIZED_SUM, |d| {
            set(d, [10, 20, 30, u64::MAX, 1, 0]);
        }),
    ];
    for name in [
        "truth_tract_bases_bound",
        "vcf_tract_bases_bound",
        "total_tract_bases_bound",
        "realized_truth_tract_bases",
        "realized_vcf_tract_bases",
        "realized_total_tract_bases",
    ] {
        cases.push(case(
            format!("{name} que no es un entero"),
            NOT_INTEGER,
            move |d| accounting(d)[name] = json!("x"),
        ));
    }
    assert_manifest_cases(&package, &cases);
    assert_manifest_accepted(
        &package,
        &[
            accepted("lo realizado igual a su cota", |d| {
                set(d, [10, 20, 30, 10, 20, 30]);
            }),
            accepted("sin nada realizado", |d| set(d, [10, 20, 30, 0, 0, 0])),
        ],
    );
}

/// **Los campos de la sección.**
#[test]
fn the_repeat_section_fields_are_judged() {
    let package = publish_mixed(5);
    assert_manifest_cases(
        &package,
        &[
            case(
                "sección que no es un objeto",
                "la sección repeat_loci no es un objeto",
                |d| *section(d) = json!(5),
            ),
            case(
                "referencia con nombre vacío",
                "la sección repeat_loci no nombra su referencia",
                |d| section(d)["reference_id"] = json!(""),
            ),
            case(
                "referencia cuyo nombre no es texto",
                "la sección repeat_loci no nombra su referencia",
                |d| section(d)["reference_id"] = json!(5),
            ),
            case(
                "sección con un campo cambiado por otro",
                "la sección repeat_loci no tiene exactamente sus campos",
                |d| {
                    let object = section(d).as_object_mut().expect("sección");
                    let value = object.remove("tract_output").expect("contabilidad");
                    object.insert("contabilidad".to_owned(), value);
                },
            ),
            case(
                "sección sin los mínimos de omisión y con un campo ajeno",
                "la sección repeat_loci no tiene exactamente sus campos",
                |d| {
                    let object = section(d).as_object_mut().expect("sección");
                    object.remove("admission_defaults");
                    object.insert("sobra".to_owned(), json!(1));
                },
            ),
            case(
                "sección sin la contabilidad de tracto y con un campo ajeno",
                "la sección repeat_loci no tiene exactamente sus campos",
                |d| {
                    let object = section(d).as_object_mut().expect("sección");
                    object.remove("tract_output");
                    object.insert("sobra".to_owned(), json!(1));
                },
            ),
            case(
                "lista de loci que no es una lista",
                "la sección repeat_loci no publica ningún locus",
                |d| section(d)["loci"] = json!(5),
            ),
            case(
                "lista de loci vacía",
                "la sección repeat_loci no publica ningún locus",
                |d| section(d)["loci"] = json!([]),
            ),
            case(
                "una procedencia de la referencia inventada",
                "la procedencia de la referencia no es simulada ni suministrada",
                |d| section(d)["reference_source"] = json!("derivada"),
            ),
            case(
                "una procedencia de la referencia que no es texto",
                "la procedencia de la referencia no es simulada ni suministrada",
                |d| section(d)["reference_source"] = json!(5),
            ),
            case(
                "dos loci con el mismo identificador",
                "dos loci publicados comparten identificador",
                |d| {
                    let first = declared(d).clone();
                    document_loci(d).push(first);
                },
            ),
        ],
    );
    assert_manifest_accepted(
        &package,
        &[accepted("una referencia simulada", |d| {
            section(d)["reference_source"] = json!("simulated");
        })],
    );
}

const LOCUS_FIELDS: &str = "un locus publicado no tiene exactamente sus campos";

/// **Cada campo de un locus publicado.**
///
/// Incluye los límites de las coordenadas (el inicio de base uno, el final que
/// llega al inicio, el tope positivo) con su lado aceptado.
#[test]
fn the_locus_fields_are_judged() {
    let package = publish_mixed(5);
    let mut cases = vec![
        case(
            "un locus que no es un objeto",
            "un locus publicado no es un objeto",
            |d| document_loci(d)[0] = json!(5),
        ),
        case("un locus con un campo de menos", LOCUS_FIELDS, |d| {
            declared(d).as_object_mut().expect("locus").remove("end");
        }),
        case("un locus con un campo de más", LOCUS_FIELDS, |d| {
            declared(d)["sobra"] = json!(1);
        }),
        case(
            "un locus con un campo cambiado por otro",
            LOCUS_FIELDS,
            |d| {
                let object = declared(d).as_object_mut().expect("locus");
                let value = object.remove("end").expect("final");
                object.insert("fin".to_owned(), value);
            },
        ),
        case(
            "un locus sin los mínimos y con un campo ajeno",
            LOCUS_FIELDS,
            |d| {
                let object = declared(d).as_object_mut().expect("locus");
                object.remove("admission");
                object.insert("sobra".to_owned(), json!(1));
            },
        ),
        case(
            "un locus con identificador vacío",
            "un locus publicado no tiene identificador",
            |d| declared(d)["id"] = json!(""),
        ),
        case(
            "un locus cuyo identificador no es texto",
            "un locus publicado no tiene identificador",
            |d| declared(d)["id"] = json!(5),
        ),
        case(
            "un locus con secuencia vacía",
            "un locus publicado no nombra su secuencia",
            |d| declared(d)["sequence_id"] = json!(""),
        ),
        case(
            "un locus cuya secuencia no es texto",
            "un locus publicado no nombra su secuencia",
            |d| declared(d)["sequence_id"] = json!(5),
        ),
        case(
            "un inicio igual a cero",
            "el inicio del locus no es una coordenada de base uno",
            |d| declared(d)["start"] = json!(0),
        ),
        case(
            "un inicio que no es un entero",
            "el inicio del locus no es una coordenada de base uno",
            |d| declared(d)["start"] = json!("x"),
        ),
        case(
            "un final justo por debajo del inicio",
            "el final del locus no llega a su inicio",
            |d| declared(d)["end"] = json!(50),
        ),
        case(
            "un final que no es un entero",
            "el final del locus no llega a su inicio",
            |d| declared(d)["end"] = json!("x"),
        ),
        case(
            "un tope igual a cero",
            "el tope del locus no es un entero positivo",
            |d| declared(d)["max_units"] = json!(0),
        ),
        case(
            "un tope que no es un entero",
            "el tope del locus no es un entero positivo",
            |d| declared(d)["max_units"] = json!("x"),
        ),
        case(
            "copias sorteadas que no son un entero",
            "las copias sorteadas del locus no son un entero",
            |d| declared(d)["drawn_copies"] = json!("x"),
        ),
        case(
            "un locus sin alelos",
            "el locus publicado no tiene alelos",
            |d| declared(d)["alleles"] = json!([]),
        ),
        case(
            "alelos que no son una lista",
            "el locus publicado no tiene alelos",
            |d| declared(d)["alleles"] = json!(5),
        ),
        case(
            "un alelo menos que lo sorteado: no suman",
            "las copias realizadas no suman las copias sorteadas del locus",
            |d| {
                declared(d)["alleles"].as_array_mut().expect("alelos").pop();
            },
        ),
    ];
    // Las longitudes de motivo contra los motivos.
    cases.push(case(
        "más longitudes de motivo que motivos",
        "las longitudes de motivo del locus no son las de sus motivos",
        |d| declared(d)["motif_lengths"] = json!([2, 2]),
    ));
    cases.push(case(
        "menos longitudes de motivo que motivos",
        "las longitudes de motivo del locus no son las de sus motivos",
        |d| {
            declared(d)["motifs"] = json!(["GT", "GT"]);
            declared(d)["motif_lengths"] = json!([2]);
        },
    ));
    cases.push(case(
        "una longitud de motivo por debajo de la de su motivo",
        "las longitudes de motivo del locus no son las de sus motivos",
        |d| declared(d)["motif_lengths"] = json!([1]),
    ));
    cases.push(case(
        "una longitud de motivo que no es la de su motivo",
        "las longitudes de motivo del locus no son las de sus motivos",
        |d| declared(d)["motif_lengths"] = json!([3]),
    ));
    // Las listas del locus, por llamante del ayudante de listas.
    for (field, list_message, item_message) in [
        (
            "motifs",
            "una lista de motivos está ausente o vacía",
            "un motivo publicado no es texto",
        ),
        (
            "motif_lengths",
            "una lista de enteros está ausente o vacía",
            "un entero publicado no lo es",
        ),
    ] {
        cases.push(case(
            format!("el locus con {field} vacío"),
            list_message,
            move |d| declared(d)[field] = json!([]),
        ));
        cases.push(case(
            format!("el locus con {field} que no es una lista"),
            list_message,
            move |d| declared(d)[field] = json!(5),
        ));
        cases.push(case(
            format!("el locus con un elemento de {field} de otro tipo"),
            item_message,
            move |d| {
                declared(d)[field] = json!([if field == "motifs" {
                    json!(5)
                } else {
                    json!("x")
                }])
            },
        ));
    }
    cases.push(case(
        "el locus con un motivo vacío",
        "un motivo publicado no son bases A, C, G o T",
        |d| {
            declared(d)["motifs"] = json!([""]);
            declared(d)["motif_lengths"] = json!([0]);
        },
    ));
    cases.push(case(
        "el locus con un motivo que no son bases",
        "un motivo publicado no son bases A, C, G o T",
        |d| declared(d)["motifs"] = json!(["GN"]),
    ));
    assert_manifest_cases(&package, &cases);
    assert_manifest_accepted(
        &package,
        &[
            accepted("un inicio igual a uno", |d| {
                declared(d)["start"] = json!(1);
            }),
            accepted("un final igual al inicio", |d| {
                declared(d)["end"] = json!(51);
            }),
            accepted("un tope de una unidad", |d| {
                declared(d)["max_units"] = json!(1);
            }),
        ],
    );
}

fn document_loci(document: &mut Value) -> &mut Vec<Value> {
    document["repeat_loci"]["loci"]
        .as_array_mut()
        .expect("loci")
}

/// **Las copias realizadas de un locus no desbordan.**
///
/// Dos alelos con todas las copias sorteadas (`u64::MAX`) y su frecuencia
/// realizada exacta (1): cada alelo es válido por sí mismo y su suma desborda.
#[test]
fn the_realized_copies_of_a_locus_cannot_overflow() {
    let package = publish_mixed(5);
    assert_manifest_cases(
        &package,
        &[case(
            "copias realizadas que desbordan al sumarse",
            "las copias realizadas del locus desbordan",
            |d| {
                declared(d)["drawn_copies"] = json!(u64::MAX);
                for allele in [allele0 as fn(&mut Value) -> &mut Value, allele1] {
                    allele(d)["realized_copies"] = json!(u64::MAX);
                    allele(d)["realized_frequency"] = json!(1.0);
                }
            },
        )],
    );
}

const ALLELE_FIELDS: &str = "un alelo publicado no tiene exactamente sus campos";

/// **Cada campo de un alelo publicado.**
#[test]
fn the_allele_fields_are_judged() {
    let package = publish_mixed(5);
    let mut cases = vec![
        case(
            "un alelo que no es un objeto",
            "un alelo publicado no es un objeto",
            |d| declared(d)["alleles"][0] = json!(5),
        ),
        case("un alelo con un campo de menos", ALLELE_FIELDS, |d| {
            allele0(d).as_object_mut().expect("alelo").remove("purity");
        }),
        case("un alelo con un campo de más", ALLELE_FIELDS, |d| {
            allele0(d)["sobra"] = json!(1);
        }),
        case(
            "un alelo con un campo cambiado por otro",
            ALLELE_FIELDS,
            |d| {
                let object = allele0(d).as_object_mut().expect("alelo");
                let value = object.remove("purity").expect("pureza");
                object.insert("pureza".to_owned(), value);
            },
        ),
        case(
            "una secuencia vacía",
            "la secuencia del alelo no son bases",
            |d| allele0(d)["sequence"] = json!(""),
        ),
        case(
            "una secuencia con una letra que no es base",
            "la secuencia del alelo no son bases",
            |d| {
                allele0(d)["sequence"] = json!("GTAAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGX")
            },
        ),
        case(
            "una secuencia en minúsculas",
            "la secuencia del alelo no son bases",
            |d| {
                allele0(d)["sequence"] = json!("gtaagtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgt")
            },
        ),
        case(
            "una secuencia que no es texto",
            "la secuencia del alelo no son bases",
            |d| allele0(d)["sequence"] = json!(5),
        ),
        case(
            "un origen inventado",
            "el origen del alelo no es declarado ni generado",
            |d| allele0(d)["origin"] = json!("inventado"),
        ),
        case(
            "unas unidades que no son un decimal",
            "las unidades del alelo no son un decimal",
            |d| allele0(d)["units"] = json!("x"),
        ),
        case(
            "unas unidades que no caben en un decimal finito",
            "las unidades del alelo no son un decimal",
            |d| allele0(d)["units"] = overflowing_decimal(),
        ),
        case(
            "una frecuencia realizada que no cabe en un decimal finito",
            "la frecuencia realizada del alelo no es un decimal",
            |d| allele0(d)["realized_frequency"] = overflowing_decimal(),
        ),
        case(
            "unas unidades iguales a cero",
            "las unidades del alelo no son un decimal positivo",
            |d| allele0(d)["units"] = json!(0.0),
        ),
        case(
            "unas unidades negativas",
            "las unidades del alelo no son un decimal positivo",
            |d| allele0(d)["units"] = json!(-1.0),
        ),
        case(
            "una tirada más larga que no es un entero",
            "la tirada más larga del alelo no es un entero",
            |d| allele0(d)["longest_run_units"] = json!("x"),
        ),
        case(
            "una tirada más larga negativa",
            "la tirada más larga del alelo no es un entero",
            |d| allele0(d)["longest_run_units"] = json!(-1),
        ),
        case(
            "unas copias realizadas que no son un entero",
            "las copias realizadas del alelo no son un entero",
            |d| allele0(d)["realized_copies"] = json!("x"),
        ),
        case(
            "una frecuencia realizada que no es un decimal",
            "la frecuencia realizada del alelo no es un decimal",
            |d| allele0(d)["realized_frequency"] = json!("x"),
        ),
        case(
            "una frecuencia realizada nula con copias sorteadas",
            "la frecuencia realizada del alelo no es un decimal",
            |d| allele0(d)["realized_frequency"] = Value::Null,
        ),
        case(
            "una frecuencia realizada de cero negativo con cero copias",
            "la frecuencia realizada no es las copias realizadas entre las sorteadas",
            |d| {
                let drawn = declared(d)["drawn_copies"].as_u64().expect("sorteadas");
                allele0(d)["realized_copies"] = json!(0);
                allele0(d)["realized_frequency"] = json!(-0.0);
                allele1(d)["realized_copies"] = json!(drawn);
                allele1(d)["realized_frequency"] = json!(1.0);
            },
        ),
        case(
            "un recuento de secuencias repetidas igual a cero",
            "el recuento de secuencias repetidas no es un entero positivo",
            |d| allele0(d)["repeat_sequence_count"] = json!(0),
        ),
        case(
            "un recuento de secuencias repetidas que no es un entero",
            "el recuento de secuencias repetidas no es un entero positivo",
            |d| allele0(d)["repeat_sequence_count"] = json!("x"),
        ),
        case(
            "copias por bloque que no son una lista",
            "las copias por bloque no son una lista",
            |d| allele0(d)["block_unit_counts"] = json!(5),
        ),
        case(
            "una cifra de copias por bloque que no es un decimal",
            "una cifra de copias por bloque no es un decimal",
            |d| allele0(d)["block_unit_counts"] = json!(["x"]),
        ),
        case(
            "un bloque sin bases",
            "un bloque sin bases no es un bloque",
            |d| allele0(d)["block_bases"] = json!([0]),
        ),
        case(
            "copias de un bloque que no son sus bases entre el motivo",
            "las copias de un bloque no son sus bases entre la longitud del motivo",
            |d| allele0(d)["block_unit_counts"] = json!([26.0]),
        ),
        case(
            "copias de un bloque justo fuera de la tolerancia",
            "las copias de un bloque no son sus bases entre la longitud del motivo",
            |d| allele0(d)["block_unit_counts"] = json!([25.000_000_1]),
        ),
        case(
            "copias de un bloque justo fuera de la tolerancia, por debajo",
            "las copias de un bloque no son sus bases entre la longitud del motivo",
            |d| allele0(d)["block_unit_counts"] = json!([24.999_999_9]),
        ),
        case(
            "copias de un bloque pequeño fuera del suelo de la tolerancia",
            "las copias de un bloque no son sus bases entre la longitud del motivo",
            |d| {
                allele1(d)["sequence"] = json!("G");
                allele1(d)["block_bases"] = json!([1]);
                allele1(d)["block_unit_counts"] = json!([0.5 + 2e-9]);
            },
        ),
        case(
            "bases de bloque que desbordan al sumarse",
            "las bases de los bloques desbordan",
            |d| {
                let half = u64::MAX as f64 / 2.0;
                let allele = allele0(d);
                allele["repeat_sequence_count"] = json!(2);
                allele["block_motifs"] = json!(["GT", "GT"]);
                allele["block_motif_lengths"] = json!([2, 2]);
                allele["block_bases"] = json!([u64::MAX, 1]);
                allele["block_unit_counts"] = json!([half, 0.5]);
            },
        ),
        case(
            "un alelo con más copias que las sorteadas",
            "un alelo recibió más copias que las sorteadas del locus",
            |d| {
                let drawn = declared(d)["drawn_copies"].as_u64().expect("sorteadas");
                allele0(d)["realized_copies"] = json!(drawn + 1);
            },
        ),
        case(
            "una frecuencia realizada no nula sin copias sorteadas",
            "sin copias sorteadas la frecuencia realizada tiene que ser nula",
            |d| {
                declared(d)["drawn_copies"] = json!(0);
                for allele in [allele0 as fn(&mut Value) -> &mut Value, allele1] {
                    allele(d)["realized_copies"] = json!(0);
                    allele(d)["realized_frequency"] = Value::Null;
                }
                allele0(d)["realized_frequency"] = json!(0.5);
            },
        ),
        case(
            "una longitud de motivo de bloque por encima de la de su motivo",
            "la longitud de un motivo de bloque no es la de su secuencia",
            |d| allele0(d)["block_motif_lengths"] = json!([3]),
        ),
        case(
            "una longitud de motivo de bloque por debajo de la de su motivo",
            "la longitud de un motivo de bloque no es la de su secuencia",
            |d| allele0(d)["block_motif_lengths"] = json!([1]),
        ),
        case(
            "bloques que suman más bases que la secuencia",
            "las bases de los bloques no suman las de la secuencia del alelo",
            |d| {
                allele0(d)["block_bases"] = json!([52]);
                allele0(d)["block_unit_counts"] = json!([26.0]);
            },
        ),
        case(
            "bloques que suman menos bases que la secuencia",
            "las bases de los bloques no suman las de la secuencia del alelo",
            |d| {
                allele0(d)["block_bases"] = json!([48]);
                allele0(d)["block_unit_counts"] = json!([24.0]);
            },
        ),
        case(
            "copias de un bloque entre uno y dos fuera de la tolerancia del suelo uno",
            "las copias de un bloque no son sus bases entre la longitud del motivo",
            |d| {
                // 3 bases de un motivo de 2: 1,5 copias. La tolerancia es
                // 1e-9 × 1,5; con 1,8e-9 de diferencia se rechaza, y con un
                // suelo de 2 (tolerancia 2e-9) se aceptaría.
                allele1(d)["sequence"] = json!("GTG");
                allele1(d)["block_bases"] = json!([3]);
                allele1(d)["block_unit_counts"] = json!([1.5 + 1.8e-9]);
            },
        ),
        case(
            "interrupciones que no son una lista",
            "las interrupciones del alelo no son una lista",
            |d| allele0(d)["interruptions"] = json!(5),
        ),
        case(
            "una interrupción que no es un objeto",
            "una interrupción publicada no es un objeto",
            |d| allele0(d)["interruptions"] = json!([5]),
        ),
        case(
            "una interrupción con un campo de menos",
            "una interrupción publicada no tiene exactamente sus campos",
            |d| allele0(d)["interruptions"] = json!([{"offset": 2}]),
        ),
        case(
            "una interrupción con un campo de más",
            "una interrupción publicada no tiene exactamente sus campos",
            |d| allele0(d)["interruptions"] = json!([{"offset": 2, "length": 2, "sobra": 1}]),
        ),
        case(
            "una interrupción con un campo cambiado por otro",
            "una interrupción publicada no tiene exactamente sus campos",
            |d| allele0(d)["interruptions"] = json!([{"offset": 2, "largo": 2}]),
        ),
        case(
            "una interrupción con un desplazamiento que no es un entero",
            "una interrupción publicada no tiene desplazamiento y longitud válidos",
            |d| allele0(d)["interruptions"] = json!([{"offset": "x", "length": 2}]),
        ),
        case(
            "una interrupción de longitud cero",
            "una interrupción publicada no tiene desplazamiento y longitud válidos",
            |d| allele0(d)["interruptions"] = json!([{"offset": 2, "length": 0}]),
        ),
        case(
            "una interrupción con una longitud que no es un entero",
            "una interrupción publicada no tiene desplazamiento y longitud válidos",
            |d| allele0(d)["interruptions"] = json!([{"offset": 2, "length": "x"}]),
        ),
    ];
    // Las proporciones del alelo, campo a campo.
    for field in ["declared_frequency", "effective_frequency", "purity"] {
        cases.push(case(
            format!("{field} que no es un decimal"),
            "una proporción del alelo no es un decimal",
            move |d| allele0(d)[field] = json!("x"),
        ));
        cases.push(case(
            format!("{field} que no cabe en un decimal finito"),
            "una proporción del alelo no es un decimal",
            move |d| allele0(d)[field] = overflowing_decimal(),
        ));
        for (what, value) in [("por encima de uno", 1.5), ("negativa", -0.1)] {
            cases.push(case(
                format!("{field} {what}"),
                "una proporción del alelo no es finita ni pertenece a [0,1]",
                move |d| allele0(d)[field] = json!(value),
            ));
        }
    }
    // Las listas del alelo, por llamante del ayudante de listas.
    for (field, list_message, item_message) in [
        (
            "block_motifs",
            "una lista de motivos está ausente o vacía",
            "un motivo publicado no es texto",
        ),
        (
            "block_motif_lengths",
            "una lista de enteros está ausente o vacía",
            "un entero publicado no lo es",
        ),
        (
            "block_bases",
            "una lista de enteros está ausente o vacía",
            "un entero publicado no lo es",
        ),
    ] {
        cases.push(case(format!("{field} vacío"), list_message, move |d| {
            allele0(d)[field] = json!([])
        }));
        cases.push(case(
            format!("{field} que no es una lista"),
            list_message,
            move |d| allele0(d)[field] = json!(5),
        ));
        cases.push(case(
            format!("un elemento de {field} de otro tipo"),
            item_message,
            move |d| {
                allele0(d)[field] = json!([if field == "block_motifs" {
                    json!(5)
                } else {
                    json!("x")
                }])
            },
        ));
    }
    for list in [
        "block_motifs",
        "block_motif_lengths",
        "block_bases",
        "block_unit_counts",
    ] {
        cases.push(case(
            format!("{list} con un elemento de más que el recuento"),
            "las listas de bloque no miden lo que el recuento declara",
            move |d| {
                let item = allele0(d)[list][0].clone();
                allele0(d)[list] = json!([item, item]);
            },
        ));
        cases.push(case(
            format!("{list} con un elemento de menos que el recuento"),
            "las listas de bloque no miden lo que el recuento declara",
            move |d| {
                allele0(d)["repeat_sequence_count"] = json!(2);
                for other in [
                    "block_motifs",
                    "block_motif_lengths",
                    "block_bases",
                    "block_unit_counts",
                ] {
                    if other != list {
                        let item = allele0(d)[other][0].clone();
                        allele0(d)[other] = json!([item, item]);
                    }
                }
            },
        ));
    }
    cases.push(case(
        "un motivo de bloque vacío",
        "un motivo publicado no son bases A, C, G o T",
        |d| {
            allele0(d)["block_motifs"] = json!([""]);
            allele0(d)["block_motif_lengths"] = json!([0]);
        },
    ));
    assert_manifest_cases(&package, &cases);

    let mut accept = vec![
        accepted("una secuencia con una N", |d| {
            allele0(d)["sequence"] = json!("GTAAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGN");
        }),
        accepted(
            "una interrupción en el desplazamiento cero y longitud uno",
            |d| allele0(d)["interruptions"] = json!([{"offset": 0, "length": 1}]),
        ),
        accepted("copias de un bloque justo dentro de la tolerancia", |d| {
            allele0(d)["block_unit_counts"] = json!([25.000_000_01]);
        }),
        accepted(
            "copias de un bloque pequeño dentro del suelo de la tolerancia",
            |d| {
                allele1(d)["sequence"] = json!("G");
                allele1(d)["block_bases"] = json!([1]);
                allele1(d)["block_unit_counts"] = json!([0.5 + 8e-10]);
            },
        ),
        accepted(
            "un alelo sin ninguna copia con su frecuencia realizada exacta",
            |d| {
                let drawn = declared(d)["drawn_copies"].as_u64().expect("sorteadas");
                allele0(d)["realized_copies"] = json!(0);
                allele0(d)["realized_frequency"] = json!(0.0);
                allele1(d)["realized_copies"] = json!(drawn);
                allele1(d)["realized_frequency"] = json!(1.0);
            },
        ),
        accepted("una sola copia sorteada con su frecuencia realizada", |d| {
            declared(d)["drawn_copies"] = json!(1);
            allele0(d)["realized_copies"] = json!(1);
            allele0(d)["realized_frequency"] = json!(1.0);
            allele1(d)["realized_copies"] = json!(0);
            allele1(d)["realized_frequency"] = json!(0.0);
        }),
        accepted("ninguna copia sorteada y las frecuencias nulas", |d| {
            declared(d)["drawn_copies"] = json!(0);
            for allele in [allele0 as fn(&mut Value) -> &mut Value, allele1] {
                allele(d)["realized_copies"] = json!(0);
                allele(d)["realized_frequency"] = Value::Null;
            }
        }),
        accepted("sin interrupciones", |d| {
            allele0(d)["interruptions"] = json!([]);
        }),
    ];
    for field in ["declared_frequency", "effective_frequency", "purity"] {
        for value in [0.0, 1.0] {
            accept.push(accepted(format!("{field} igual a {value}"), move |d| {
                allele0(d)[field] = json!(value)
            }));
        }
    }
    assert_manifest_accepted(&package, &accept);
}

/// **El bloque del generador, campo a campo.**
#[test]
fn the_generator_block_is_judged_field_by_field() {
    let package = publish_mixed(5);
    fn slippage(d: &mut Value) -> &mut Value {
        &mut generated(d)["slippage"]
    }
    let mut cases = vec![
        case(
            "un bloque del generador que no es un objeto",
            "el bloque del generador no es un objeto",
            |d| *slippage(d) = json!(5),
        ),
        case(
            "una media de pasos que no cabe en un decimal finito",
            "la media de pasos no es un decimal",
            |d| slippage(d)["mean_events"] = overflowing_decimal(),
        ),
        case(
            "un sesgo de expansión que no cabe en un decimal finito",
            "el sesgo de expansión no es un decimal",
            |d| slippage(d)["expansion_bias"] = overflowing_decimal(),
        ),
        case(
            "un bloque del generador con un campo de más",
            "el bloque del generador no tiene exactamente sus campos",
            |d| slippage(d)["sobra"] = json!(1),
        ),
        case(
            "un bloque del generador con un campo de menos",
            "el bloque del generador no tiene exactamente sus campos",
            |d| {
                slippage(d)
                    .as_object_mut()
                    .expect("bloque")
                    .remove("lineages");
            },
        ),
        case(
            "un bloque del generador con un campo cambiado por otro",
            "el bloque del generador no tiene exactamente sus campos",
            |d| {
                let object = slippage(d).as_object_mut().expect("bloque");
                let value = object.remove("lineages").expect("linajes");
                object.insert("linajes".to_owned(), value);
            },
        ),
        case(
            "un modelo del generador que no es el escalonado",
            "el modelo del generador no es el escalonado",
            |d| slippage(d)["model"] = json!("infinite_alleles"),
        ),
        case(
            "un modelo del generador que no es texto",
            "el modelo del generador no es el escalonado",
            |d| slippage(d)["model"] = json!(5),
        ),
        case(
            "un paso del generador de dos unidades",
            "el paso del generador no es de una unidad completa",
            |d| slippage(d)["step_units"] = json!(2),
        ),
        case(
            "un paso del generador que no es un entero",
            "el paso del generador no es de una unidad completa",
            |d| slippage(d)["step_units"] = json!("x"),
        ),
        case(
            "un locus con generador cuyos alelos se dicen declarados",
            "el origen del alelo contradice al del locus",
            |d| generated(d)["alleles"][0]["origin"] = json!("declared"),
        ),
        case(
            "un locus declarado cuyos alelos se dicen generados",
            "el origen del alelo contradice al del locus",
            |d| allele0(d)["origin"] = json!("generated"),
        ),
        case(
            "un locus sin bloque del generador con alelos generados",
            "el origen del alelo contradice al del locus",
            |d| {
                generated(d)
                    .as_object_mut()
                    .expect("locus")
                    .remove("slippage");
            },
        ),
        case(
            "una media de pasos que no es un decimal",
            "la media de pasos no es un decimal",
            |d| slippage(d)["mean_events"] = json!("x"),
        ),
        case(
            "un sesgo de expansión que no es un decimal",
            "el sesgo de expansión no es un decimal",
            |d| slippage(d)["expansion_bias"] = json!("x"),
        ),
        case(
            "unas limitaciones que no son una lista",
            "las limitaciones del generador no son una lista",
            |d| slippage(d)["limitations"] = json!("x"),
        ),
        case(
            "unos linajes que no son un entero",
            "los linajes del generador no son un entero positivo",
            |d| slippage(d)["lineages"] = json!("x"),
        ),
        case(
            "una media de pasos justo por debajo de cero",
            "la media de pasos no es un decimal no negativo",
            |d| slippage(d)["mean_events"] = json!(-1e-300),
        ),
        case(
            "un sesgo justo por encima de uno",
            "el sesgo de expansión no pertenece a [0,1]",
            |d| slippage(d)["expansion_bias"] = json!(1.000_000_000_000_000_2),
        ),
        case(
            "un sesgo justo por debajo de cero",
            "el sesgo de expansión no pertenece a [0,1]",
            |d| slippage(d)["expansion_bias"] = json!(-1e-300),
        ),
        case(
            "un linaje menos que alelos",
            "el abanico publica más alelos que linajes simulados",
            |d| {
                let alleles = generated(d)["alleles"].as_array().expect("alelos").len();
                slippage(d)["lineages"] = json!(alleles - 1);
            },
        ),
        case(
            "una limitación de menos",
            "las limitaciones del generador no son las tres declaradas",
            |d| {
                slippage(d)["limitations"]
                    .as_array_mut()
                    .expect("limitaciones")
                    .pop();
            },
        ),
        case(
            "cero linajes",
            "los linajes del generador no son un entero positivo",
            |d| slippage(d)["lineages"] = json!(0),
        ),
        case(
            "un solo linaje con más alelos que linajes",
            "el abanico publica más alelos que linajes simulados",
            |d| slippage(d)["lineages"] = json!(1),
        ),
        case(
            "una limitación de más",
            "las limitaciones del generador no son las tres declaradas",
            |d| {
                slippage(d)["limitations"]
                    .as_array_mut()
                    .expect("limitaciones")
                    .push(json!("una cuarta"));
            },
        ),
        case(
            "una limitación que no es texto",
            "las limitaciones del generador no son las tres declaradas",
            |d| slippage(d)["limitations"][2] = json!(5),
        ),
        case(
            "una limitación cambiada de orden",
            "las limitaciones del generador no son las tres declaradas",
            |d| {
                slippage(d)["limitations"]
                    .as_array_mut()
                    .expect("limitaciones")
                    .swap(0, 2);
            },
        ),
    ];
    for (field, message) in [
        ("drawn_events", "los pasos sorteados no son un entero"),
        ("applied_events", "los pasos aplicados no son un entero"),
        (
            "unapplied_events",
            "los pasos no aplicados no son un entero",
        ),
    ] {
        cases.push(case(
            format!("{field} que no es un entero"),
            message,
            move |d| slippage(d)[field] = json!("x"),
        ));
    }
    for field in [
        "algorithm",
        "count_algorithm",
        "bias_algorithm",
        "rng_scope",
        "rng_scope_catalog_version",
    ] {
        cases.push(case(
            format!("{field} cambiado"),
            "una versión o un ámbito del generador no es el normativo",
            move |d| slippage(d)[field] = json!("otro"),
        ));
        cases.push(case(
            format!("{field} que no es texto"),
            "una versión o un ámbito del generador no es el normativo",
            move |d| slippage(d)[field] = json!(5),
        ));
    }
    cases.push(case(
        "una contabilidad con un paso aplicado de menos",
        "los pasos aplicados y los no aplicados no suman los sorteados",
        |d| {
            let applied = slippage(d)["applied_events"].as_u64().expect("aplicados");
            slippage(d)["applied_events"] = json!(applied - 1);
        },
    ));
    cases.push(case(
        "una contabilidad cuya suma desborda",
        "los pasos aplicados y los no aplicados no suman los sorteados",
        |d| {
            slippage(d)["applied_events"] = json!(u64::MAX);
            slippage(d)["unapplied_events"] = json!(1);
            slippage(d)["drawn_events"] = json!(0);
        },
    ));
    assert_manifest_cases(&package, &cases);
    assert_manifest_accepted(
        &package,
        &[
            accepted("una media de pasos de cero", |d| {
                slippage(d)["mean_events"] = json!(0.0);
            }),
            accepted("un sesgo de cero", |d| {
                slippage(d)["expansion_bias"] = json!(0.0);
            }),
            accepted("un sesgo de uno", |d| {
                slippage(d)["expansion_bias"] = json!(1.0);
            }),
            accepted("tantos linajes como alelos", |d| {
                // Cada linaje termina en un alelo: con tantos linajes como
                // alelos publicados la relación se cumple justo.
                let alleles = generated(d)["alleles"].as_array().expect("alelos").len();
                slippage(d)["lineages"] = json!(alleles);
            }),
        ],
    );
}

/// **El intento sobre las guardas de la sección que ningún manifiesto
/// publicado y manipulado puede alcanzar.** No discrimina esas guardas: solo
/// muestra qué otra comprobación gana la carrera.
///
/// - Los mínimos de omisión, la contabilidad de tracto, los mínimos del locus
///   y la frecuencia realizada de un alelo se leen con `get(..)` **después** de
///   comprobar que el objeto tiene exactamente sus campos: ninguno puede faltar
///   sin que ese cierre rechace antes.
/// - El origen del alelo lo exige ya `validate_repeat_allele`, que corre antes
///   que la lectura del locus.
#[test]
fn unreachable_guards_are_shadowed_by_other_checks() {
    let package = publish_mixed(5);
    assert_manifest_cases(
        &package,
        &[
            case(
                "sin mínimos de omisión",
                "la sección repeat_loci no tiene exactamente sus campos",
                |d| {
                    section(d)
                        .as_object_mut()
                        .expect("sección")
                        .remove("admission_defaults");
                },
            ),
            case(
                "sin contabilidad de tracto",
                "la sección repeat_loci no tiene exactamente sus campos",
                |d| {
                    section(d)
                        .as_object_mut()
                        .expect("sección")
                        .remove("tract_output");
                },
            ),
            case("un locus sin mínimos", LOCUS_FIELDS, |d| {
                declared(d)
                    .as_object_mut()
                    .expect("locus")
                    .remove("admission");
            }),
            case("un alelo sin frecuencia realizada", ALLELE_FIELDS, |d| {
                allele0(d)
                    .as_object_mut()
                    .expect("alelo")
                    .remove("realized_frequency");
            }),
            case("un alelo sin origen", ALLELE_FIELDS, |d| {
                allele0(d).as_object_mut().expect("alelo").remove("origin");
            }),
        ],
    );
}

// ---------------------------------------------------------------------------
// El abanico generado por deslizamiento: el juez que lo rehace por su cuenta
// ---------------------------------------------------------------------------

fn judge_resources(package: &Package, effective: &Value) -> Result<(), Diagnostic> {
    validate_published_resources(
        &package.destination,
        &package.effective_path(),
        &serde_json::to_vec_pretty(effective).expect("JSON"),
        &package.manifest_path(),
        &package.bytes("manifest.json"),
        None,
    )
}

/// El locus generado de la configuración efectiva (`STR-RIGHT`, el segundo).
fn effective_generated(document: &mut Value) -> &mut Value {
    &mut document["repeat_loci"][1]
}

fn effective_allele(document: &mut Value, index: usize) -> &mut Value {
    &mut document["repeat_loci"][1]["alleles"][index]
}

/// El rechazo del juez del abanico del locus `STR-RIGHT`: el juez antepone
/// siempre el identificador del locus al motivo.
fn fan(motive: &str) -> String {
    format!("abanico generado de STR-RIGHT: {motive}")
}

/// Lo mismo, para un motivo literal.
macro_rules! fan {
    ($motive:literal) => {
        concat!("abanico generado de STR-RIGHT: ", $motive)
    };
}

fn assert_effective_cases(package: &Package, cases: &[Case]) {
    let original = package.json("effective-config.json");
    judge_resources(package, &original).expect("control positivo: el paquete intacto se acepta");
    for case in cases {
        let document = altered(&original, &case.label, case.mutate.as_ref());
        assert_eq!(
            detail(judge_resources(package, &document)),
            case.expected,
            "«{}» se rechazó por otra razón",
            case.label
        );
    }
}

/// **Cada dato que el juez del abanico lee del locus generado de la
/// configuración efectiva.**
///
/// La configuración efectiva se altera en un punto y el juez completo
/// (`validate_published_resources`) tiene que rechazarla por el motivo exacto
/// de ese punto: el recorrido del abanico no llega a ejecutarse.
#[test]
fn the_fan_judge_rejects_a_generated_locus_it_cannot_read() {
    let package = publish_mixed(5);
    let cases = [
        case(
            "sin identificador de locus generado",
            "abanico generado de locus sin id: el locus no tiene identificador",
            |d| {
                effective_generated(d)
                    .as_object_mut()
                    .expect("locus")
                    .remove("id");
            },
        ),
        case(
            "un identificador de locus que no es texto",
            "abanico generado de locus sin id: el locus no tiene identificador",
            |d| effective_generated(d)["id"] = json!(5),
        ),
        case(
            "un tope que no es un entero",
            fan!("el tope del locus no es un entero"),
            |d| effective_generated(d)["max_units"] = json!("x"),
        ),
        case(
            "un tope que falta",
            fan!("el tope del locus no es un entero"),
            |d| {
                effective_generated(d)
                    .as_object_mut()
                    .expect("locus")
                    .remove("max_units");
            },
        ),
        case(
            "una media de pasos que no es un decimal",
            fan!("la media de pasos no es un decimal"),
            |d| effective_generated(d)["slippage"]["mean_events"] = json!("x"),
        ),
        case(
            "un sesgo de expansión que no es un decimal",
            fan!("el sesgo de expansión no es un decimal"),
            |d| effective_generated(d)["slippage"]["expansion_bias"] = json!("x"),
        ),
        case(
            "unos linajes que no son un entero",
            fan!("los linajes no son un entero"),
            |d| effective_generated(d)["slippage"]["lineages"] = json!("x"),
        ),
        case(
            "unos linajes que no son las copias haplotípicas publicadas",
            fan!("los linajes no son las copias haplotípicas de las poblaciones iniciales"),
            |d| effective_generated(d)["slippage"]["lineages"] = json!(11),
        ),
        case(
            "unos linajes que son un entero distinto de las copias, por debajo",
            fan!("los linajes no son las copias haplotípicas de las poblaciones iniciales"),
            |d| effective_generated(d)["slippage"]["lineages"] = json!(9),
        ),
        case(
            "una media de pasos justo por encima del tope",
            fan!("la media de pasos rebasa el tope que el generador admite"),
            |d| effective_generated(d)["slippage"]["mean_events"] = json!(next_up(1000.0)),
        ),
        case(
            "la media máxima admitida: pasa la cota y cae en el abanico",
            fan!("el número de alelos publicados no es el que el modelo produce"),
            |d| effective_generated(d)["slippage"]["mean_events"] = json!(1000.0),
        ),
        case(
            "una media de cero: pasa la cota y cae en el abanico",
            fan!("el número de alelos publicados no es el que el modelo produce"),
            |d| effective_generated(d)["slippage"]["mean_events"] = json!(0.0),
        ),
        case(
            "una media de pasos negativa",
            fan!("la media de pasos rebasa el tope que el generador admite"),
            |d| effective_generated(d)["slippage"]["mean_events"] = json!(-1e-300),
        ),
        case(
            "un sesgo fuera de [0,1] distinto de los extremos",
            fan!("el ensayo del sesgo falló: InvalidProbability"),
            |d| effective_generated(d)["slippage"]["expansion_bias"] = json!(2.0),
        ),
        case(
            "un locus sin alelos publicados",
            fan!("el locus no publica alelos"),
            |d| {
                effective_generated(d)
                    .as_object_mut()
                    .expect("locus")
                    .remove("alleles");
            },
        ),
        case(
            "unos alelos que no son una lista",
            fan!("el locus no publica alelos"),
            |d| effective_generated(d)["alleles"] = json!(5),
        ),
        case(
            "un alelo publicado de menos",
            fan!("el número de alelos publicados no es el que el modelo produce"),
            |d| {
                effective_generated(d)["alleles"]
                    .as_array_mut()
                    .expect("alelos")
                    .pop();
            },
        ),
        case(
            "un alelo publicado de más",
            fan!("el número de alelos publicados no es el que el modelo produce"),
            |d| {
                let extra = effective_allele(d, 0).clone();
                effective_generated(d)["alleles"]
                    .as_array_mut()
                    .expect("alelos")
                    .push(extra);
            },
        ),
        case(
            "un alelo sin secuencia",
            fan!("un alelo publicado no tiene secuencia"),
            |d| {
                effective_allele(d, 1)
                    .as_object_mut()
                    .expect("alelo")
                    .remove("sequence");
            },
        ),
        case(
            "un alelo cuya secuencia no es texto",
            fan!("un alelo publicado no tiene secuencia"),
            |d| effective_allele(d, 1)["sequence"] = json!(5),
        ),
        case(
            "dos alelos en otro orden que el canónico",
            fan!(
                "la secuencia de un alelo publicado no es la que el modelo produce, o el orden canónico no es el suyo"
            ),
            |d| {
                effective_generated(d)["alleles"]
                    .as_array_mut()
                    .expect("alelos")
                    .swap(0, 1);
            },
        ),
        case(
            "un alelo con una base cambiada",
            fan!(
                "la secuencia de un alelo publicado no es la que el modelo produce, o el orden canónico no es el suyo"
            ),
            |d| {
                let sequence = effective_allele(d, 0)["sequence"]
                    .as_str()
                    .expect("secuencia")
                    .to_owned();
                let mut bytes = sequence.into_bytes();
                bytes[0] = b'T';
                effective_allele(d, 0)["sequence"] =
                    json!(String::from_utf8(bytes).expect("UTF-8"));
            },
        ),
    ];
    assert_effective_cases(&package, &cases);
}

/// **El tramo ancestral que el juez saca de la referencia publicada.**
///
/// El intervalo del locus se lee de la configuración efectiva y se corta de las
/// bases del FASTA publicado: cada dato ausente, cada intervalo que no cabe y
/// cada motivo que no puede estructurar ese tramo se rechaza con su motivo.
#[test]
fn the_fan_judge_rejects_an_ancestral_tract_it_cannot_cut_from_the_reference() {
    let package = publish_mixed(5);
    let cases = [
        case(
            "un locus sin secuencia",
            fan!("el locus no nombra su secuencia"),
            |d| {
                effective_generated(d)
                    .as_object_mut()
                    .expect("locus")
                    .remove("sequence_id");
            },
        ),
        case(
            "un locus cuya secuencia no es texto",
            fan!("el locus no nombra su secuencia"),
            |d| effective_generated(d)["sequence_id"] = json!(5),
        ),
        case(
            "una secuencia que la referencia no contiene",
            fan!("la referencia publicada no contiene esa secuencia"),
            |d| effective_generated(d)["sequence_id"] = json!("chrX"),
        ),
        case(
            "un inicio igual a cero",
            fan!("el inicio no es una coordenada de base uno"),
            |d| effective_generated(d)["start"] = json!(0),
        ),
        case(
            "un inicio que no es un entero",
            fan!("el inicio no es una coordenada de base uno"),
            |d| effective_generated(d)["start"] = json!("x"),
        ),
        case(
            "un inicio igual a uno se lee como base uno y el tramo de una G no es el abanico",
            fan!("el número de alelos publicados no es el que el modelo produce"),
            |d| {
                effective_generated(d)["start"] = json!(1);
                effective_generated(d)["end"] = json!(1);
            },
        ),
        case(
            "un final que no es un entero",
            fan!("el final no es una coordenada"),
            |d| effective_generated(d)["end"] = json!("x"),
        ),
        case(
            "un final que falta",
            fan!("el final no es una coordenada"),
            |d| {
                effective_generated(d)
                    .as_object_mut()
                    .expect("locus")
                    .remove("end");
            },
        ),
        case(
            "un final justo más allá de la secuencia",
            fan!("el intervalo rebasa la secuencia publicada"),
            |d| effective_generated(d)["end"] = json!(301),
        ),
        case(
            "un inicio posterior al final",
            fan!("el intervalo rebasa la secuencia publicada"),
            |d| effective_generated(d)["start"] = json!(250),
        ),
        case(
            "sin motivos publicados",
            fan!("el locus no publica sus motivos"),
            |d| {
                effective_generated(d)
                    .as_object_mut()
                    .expect("locus")
                    .remove("motifs");
            },
        ),
        case(
            "unos motivos que no son una lista",
            fan!("el locus no publica sus motivos"),
            |d| effective_generated(d)["motifs"] = json!("CA"),
        ),
        case(
            "una lista de motivos vacía",
            fan!("un locus generado tiene que tener un solo motivo"),
            |d| effective_generated(d)["motifs"] = json!([]),
        ),
        case(
            "dos motivos",
            fan!("un locus generado tiene que tener un solo motivo"),
            |d| effective_generated(d)["motifs"] = json!(["CA", "GT"]),
        ),
        case(
            "un motivo que no es texto",
            fan!("el motivo publicado no es texto"),
            |d| effective_generated(d)["motifs"] = json!([5]),
        ),
        case("un motivo vacío", fan!("the repeat motif is empty"), |d| {
            effective_generated(d)["motifs"] = json!([""])
        }),
        case(
            "un motivo con una base que no es A, C, G ni T",
            fan!("the repeat motif contains 0x4e, which is not A, C, G or T"),
            |d| effective_generated(d)["motifs"] = json!(["N"]),
        ),
        case(
            "un motivo en minúsculas",
            fan!("the repeat motif contains 0x63, which is not A, C, G or T"),
            |d| effective_generated(d)["motifs"] = json!(["ca"]),
        ),
        case(
            "un tramo vacío: el final justo antes del inicio",
            fan!("the repeat tract is empty"),
            |d| {
                effective_generated(d)["start"] = json!(101);
                effective_generated(d)["end"] = json!(100);
            },
        ),
    ];
    assert_effective_cases(&package, &cases);
}

/// **La contabilidad y la frecuencia de cada alelo se comparan con las que el
/// juez rehace.**
///
/// La frecuencia se compara **bit a bit** (`f64 != f64` sobre lo publicado):
/// el siguiente decimal representable de cada una se rechaza, alelo a alelo,
/// y el diagnóstico lleva las dos cifras. Los tres contadores de eventos se
/// alteran de uno en uno y, a la vez, en la contabilidad **que sí cuadra
/// consigo misma** del manifiesto: este rechazo es el del juez, no el de la
/// suma del lector.
#[test]
fn the_fan_judge_rejects_a_frequency_or_an_accounting_it_did_not_reproduce() {
    let package = publish_mixed(5);
    let original = package.json("effective-config.json");
    judge_resources(&package, &original).expect("control positivo");
    let lineages = original["repeat_loci"][1]["slippage"]["lineages"]
        .as_u64()
        .expect("linajes");
    let alleles = original["repeat_loci"][1]["alleles"]
        .as_array()
        .expect("alelos");
    assert!(
        alleles.len() >= 3,
        "el fixture publica un abanico de varios alelos"
    );
    for (index, allele) in alleles.iter().enumerate() {
        let frequency = allele["declared_frequency"].as_f64().expect("frecuencia");
        let share = (frequency * lineages as f64).round() as u64;
        assert_eq!(
            share as f64 / lineages as f64,
            frequency,
            "el fixture: la frecuencia publicada es una proporción de linajes"
        );
        for (what, wrong) in [
            ("el siguiente decimal", next_up(frequency)),
            ("cero", 0.0),
            ("uno", 1.0),
        ] {
            let mut document = original.clone();
            effective_allele(&mut document, index)["declared_frequency"] = json!(wrong);
            assert_ne!(document, original);
            assert_eq!(
                detail(judge_resources(&package, &document)),
                fan(&format!(
                    "la frecuencia de un alelo no es la proporción de linajes que terminan en él (alelo {index}: publicada Some({wrong:?}), rehecha {frequency} = {share}/{lineages})"
                )),
                "alelo {index} con {what}"
            );
        }
        // Sin frecuencia, el rechazo es el mismo y la publicada es `None`.
        let mut document = original.clone();
        effective_allele(&mut document, index)
            .as_object_mut()
            .expect("alelo")
            .remove("declared_frequency");
        assert_eq!(
            detail(judge_resources(&package, &document)),
            fan(&format!(
                "la frecuencia de un alelo no es la proporción de linajes que terminan en él (alelo {index}: publicada None, rehecha {frequency} = {share}/{lineages})"
            )),
            "alelo {index} sin frecuencia"
        );
    }
    for field in ["drawn_events", "applied_events", "unapplied_events"] {
        for (what, delta) in [("uno más", 1_i64), ("uno menos", -1)] {
            let current = original["repeat_loci"][1]["slippage"][field]
                .as_u64()
                .expect("contador");
            let Some(wrong) = current.checked_add_signed(delta) else {
                continue;
            };
            let mut document = original.clone();
            effective_generated(&mut document)["slippage"][field] = json!(wrong);
            assert_eq!(
                detail(judge_resources(&package, &document)),
                fan("la contabilidad publicada del generador no es la que el modelo produce"),
                "{field} con {what}"
            );
        }
        let mut document = original.clone();
        effective_generated(&mut document)["slippage"][field] = json!("x");
        assert_eq!(
            detail(judge_resources(&package, &document)),
            fan("la contabilidad publicada del generador no es la que el modelo produce"),
            "{field} que no es un entero"
        );
    }
}

/// **El intento sobre las guardas del juez del abanico que ningún paquete
/// publicado y manipulado puede alcanzar.** No discrimina esas guardas: solo
/// muestra qué otra comprobación gana la carrera.
///
/// - El identificador de la referencia lo exige no vacío `read_physical_reference`
///   antes de que el juez del abanico lo lea, y el del documento de inserción
///   `judge_inserted_reference`.
/// - Las poblaciones iniciales y la matriz de un paquete publicado las ata la
///   comparación con el manifiesto y con el TOML antes de leer la referencia
///   (las bases no materializadas **sí** se alcanzan, con un paquete rehecho
///   entero: `a_package_forged_without_populations_reaches_the_unmaterialized_reference`).
#[test]
fn unreachable_fan_guards_are_shadowed_by_other_checks() {
    let package = publish_mixed(5);
    let original = package.json("effective-config.json");
    for (what, mutate, expected) in [
        (
            "referencia sin identificador",
            (|d: &mut Value| {
                d["reference"]
                    .as_object_mut()
                    .expect("referencia")
                    .remove("id");
            }) as fn(&mut Value),
            "campos ausentes o desconocidos en reference",
        ),
        (
            "referencia con identificador vacío",
            |d: &mut Value| d["reference"]["id"] = json!(""),
            "id falta o está vacío",
        ),
        (
            "sin poblaciones iniciales ni matriz",
            |d: &mut Value| {
                d["initial_populations"] = json!([]);
                d["population_snvs"] = Value::Null;
            },
            "population_snvs no coincide con el manifiesto",
        ),
        (
            "sin poblaciones iniciales",
            |d: &mut Value| d["initial_populations"] = json!([]),
            "las cohortes iniciales efectivas no coinciden exactamente con el TOML",
        ),
    ] {
        let mut document = original.clone();
        mutate(&mut document);
        assert_ne!(document, original, "{what}");
        assert_eq!(
            detail(judge_resources(&package, &document)),
            expected,
            "{what}"
        );
    }
}

/// **Un locus generado sin poblaciones iniciales no llega a publicarse**: el
/// programa nunca produce el paquete que no materializa la referencia. No
/// discrimina ninguna guarda del juez.
#[test]
fn a_generated_locus_without_populations_is_not_even_published() {
    let root = std::env::temp_dir().join(format!(
        "gps-results-repeat-nopop-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("directorio");
    fs::write(root.join("reference.fa"), mixed_reference()).expect("FASTA");
    let config = format!(
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
length = 300

[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5
"#
    );
    fs::write(root.join("general.toml"), config).expect("configuración");
    let validated = load_and_validate_config(&root.join("general.toml")).expect("configuración");
    let error = execute_validated_config(&validated, &root.join("published"))
        .expect_err("un generador sin copias haplotípicas se rechaza al publicar");
    assert_eq!(error.code(), "GPS-CONF-008");
    assert_eq!(
        error.context().get("detail").map(String::as_str),
        Some(
            "locus STR-RIGHT: el generador simula un linaje por copia haplotípica y el escenario no declara ninguna"
        )
    );
    fs::remove_dir_all(root).ok();
}

// ---------------------------------------------------------------------------
// La configuración efectiva: tartamudeo, curva de longitud de tracto,
// tartamudeo clonal y amplificación
// ---------------------------------------------------------------------------

fn judge_effective(package: &Package, effective: &Value) -> Result<(), Diagnostic> {
    validate_effective_config_against_manifest_bytes(
        &package.effective_path(),
        &serde_json::to_vec_pretty(effective).expect("JSON"),
        &package.manifest_path(),
        &package.bytes("manifest.json"),
    )
}

fn assert_technical_cases(package: &Package, cases: &[Case]) {
    let original = package.json("effective-config.json");
    judge_effective(package, &original)
        .expect("control positivo: la configuración efectiva intacta se acepta");
    for case in cases {
        let document = altered(&original, &case.label, case.mutate.as_ref());
        assert_eq!(
            detail(judge_effective(package, &document)),
            case.expected,
            "«{}» se rechazó por otra razón",
            case.label
        );
    }
}

fn assert_technical_accepted(package: &Package, cases: &[Accepted]) {
    let original = package.json("effective-config.json");
    for case in cases {
        let document = altered(&original, &case.label, case.mutate.as_ref());
        judge_effective(package, &document)
            .unwrap_or_else(|error| panic!("«{}» se tenía que aceptar: {error:?}", case.label));
    }
}

/// La sección del tartamudeo de base (`technical_event.stutter`).
fn stutter(document: &mut Value) -> &mut Value {
    &mut document["technical_event"]["stutter"]
}

/// La curva de longitud de tracto del tartamudeo de base.
fn dependence(document: &mut Value) -> &mut Value {
    &mut document["technical_event"]["stutter"]["tract_length_dependence"]
}

fn pcr(document: &mut Value) -> &mut Value {
    &mut document["technical_event"]["artifacts"]["pcr"]
}

/// El tartamudeo clonal de la amplificación.
fn clonal(document: &mut Value) -> &mut Value {
    &mut document["technical_event"]["artifacts"]["pcr"]["stutter"]
}

/// **El tartamudeo de base: algoritmo, procedencia, coeficiente, tasas y
/// pesos.**
#[test]
fn the_effective_stutter_section_is_judged() {
    let package = publish_amplified("");
    let mut cases = vec![
        case(
            "un algoritmo de tartamudeo que no es el publicado",
            "algoritmo de tartamudeo inválido",
            |d| stutter(d)["algorithm"] = json!("stutter-v1"),
        ),
        case(
            "un catálogo de azar que no es el del tartamudeo",
            "algoritmo de tartamudeo inválido",
            |d| stutter(d)["rng_catalog"] = json!("f4-technical-v5"),
        ),
        case(
            "un coeficiente de tirada publicado como calibrado",
            "el coeficiente de la tirada tiene que publicarse como no calibrado",
            |d| stutter(d)["run_length_coefficient"]["calibrated"] = json!(true),
        ),
        case(
            "un coeficiente de tirada con la marca en otro tipo",
            "el coeficiente de la tirada tiene que publicarse como no calibrado",
            |d| stutter(d)["run_length_coefficient"]["calibrated"] = json!("false"),
        ),
        case(
            "unas tasas por longitud de motivo que no son una lista",
            "faltan las tasas por longitud de motivo",
            |d| stutter(d)["motif_length_rates"] = json!(5),
        ),
        case(
            "una tasa por longitud de motivo de menos",
            "la calibración publica exactamente seis tasas por longitud de motivo",
            |d| {
                stutter(d)["motif_length_rates"]
                    .as_array_mut()
                    .expect("tasas")
                    .pop();
            },
        ),
        case(
            "una tasa por longitud de motivo de más",
            "la calibración publica exactamente seis tasas por longitud de motivo",
            |d| {
                let extra = stutter(d)["motif_length_rates"][5].clone();
                stutter(d)["motif_length_rates"]
                    .as_array_mut()
                    .expect("tasas")
                    .push(extra);
            },
        ),
        case(
            "unos pesos de los pasos que no son una lista",
            "faltan los pesos de los pasos",
            |d| stutter(d)["step_weights"] = json!(5),
        ),
        case(
            "un peso sin sus unidades",
            "un peso publicado no declara sus unidades",
            |d| {
                stutter(d)["step_weights"][0]
                    .as_object_mut()
                    .expect("peso")
                    .remove("units");
            },
        ),
        case(
            "un peso con unidades que no son un entero",
            "un peso publicado no declara sus unidades",
            |d| stutter(d)["step_weights"][0]["units"] = json!("x"),
        ),
        case(
            "un peso que no es un objeto",
            "un peso publicado no declara sus unidades",
            |d| stutter(d)["step_weights"][0] = json!(5),
        ),
        case(
            "un paso de menos",
            "los pesos publicados no cubren exactamente la estructura del algoritmo",
            |d| {
                stutter(d)["step_weights"]
                    .as_array_mut()
                    .expect("pesos")
                    .pop();
            },
        ),
        case(
            "un paso de más",
            "los pesos publicados no cubren exactamente la estructura del algoritmo",
            |d| {
                stutter(d)["step_weights"]
                    .as_array_mut()
                    .expect("pesos")
                    .push(json!({"units": 2, "weight": 0.0}));
            },
        ),
        case(
            "los pasos en otro orden",
            "los pesos publicados no cubren exactamente la estructura del algoritmo",
            |d| {
                stutter(d)["step_weights"]
                    .as_array_mut()
                    .expect("pesos")
                    .swap(0, 3);
            },
        ),
        case(
            "un paso con otras unidades",
            "los pesos publicados no cubren exactamente la estructura del algoritmo",
            |d| stutter(d)["step_weights"][3]["units"] = json!(2),
        ),
    ];
    for field in ["source", "species", "technology", "column"] {
        cases.push(case(
            format!("calibración con {field} cambiado"),
            "la procedencia de la calibración del tartamudeo no es la publicada",
            move |d| stutter(d)["calibration"][field] = json!("otra"),
        ));
    }
    assert_technical_cases(&package, &cases);
}

/// **La curva de longitud de tracto del tartamudeo de base.**
///
/// Procedencia, tosquedad, motivo de la tosquedad, no extrapolación y la tabla
/// por longitud de motivo: seis filas en orden, cada una con su marca de curva
/// medida, su longitud de referencia, su tramo y su tabla de factores.
#[test]
fn the_tract_length_dependence_is_judged() {
    let package = publish_amplified("");
    let mut cases = vec![
        case(
            "una procedencia de la curva que no es la publicada",
            "la procedencia de la curva de longitud de tracto no es la publicada",
            |d| dependence(d)["source"] = json!("otra"),
        ),
        case(
            "una curva que no se marca como tosca",
            "la dependencia con la longitud del tracto se publica marcada como tosca",
            |d| dependence(d)["coarse"] = json!(false),
        ),
        case(
            "una marca de tosquedad que no es un booleano",
            "la dependencia con la longitud del tracto se publica marcada como tosca",
            |d| dependence(d)["coarse"] = json!("true"),
        ),
        case(
            "una marca de extrapolación que no es un booleano",
            "fuera del tramo medido se mantiene el extremo y no se extrapola",
            |d| dependence(d)["extrapolated"] = json!("false"),
        ),
        case(
            "un motivo de la tosquedad cambiado",
            "el motivo de la tosquedad de la curva no es el publicado",
            |d| dependence(d)["coarse_reason"] = json!("otro"),
        ),
        case(
            "una curva que se extrapola",
            "fuera del tramo medido se mantiene el extremo y no se extrapola",
            |d| dependence(d)["extrapolated"] = json!(true),
        ),
        case(
            "una tabla de longitudes que no es una lista",
            "falta la tabla de longitudes de tracto de referencia",
            |d| dependence(d)["by_motif_length"] = json!(5),
        ),
        case(
            "una fila de menos en la tabla de longitudes",
            "la tabla de longitudes de referencia no cubre el tramo calibrado",
            |d| {
                dependence(d)["by_motif_length"]
                    .as_array_mut()
                    .expect("tabla")
                    .pop();
            },
        ),
        case(
            "una fila de más en la tabla de longitudes",
            "la tabla de longitudes de referencia no cubre el tramo calibrado",
            |d| {
                let extra = dependence(d)["by_motif_length"][5].clone();
                dependence(d)["by_motif_length"]
                    .as_array_mut()
                    .expect("tabla")
                    .push(extra);
            },
        ),
        case(
            "una fila que no es un objeto",
            "una fila de la tabla de longitudes no es objeto",
            |d| dependence(d)["by_motif_length"][2] = json!(5),
        ),
        case(
            "dos filas intercambiadas",
            "la tabla de longitudes no va en orden de longitud de motivo",
            |d| {
                dependence(d)["by_motif_length"]
                    .as_array_mut()
                    .expect("tabla")
                    .swap(0, 1);
            },
        ),
        case(
            "una fila sin la marca de curva medida",
            "una fila no declara si tiene curva medida",
            |d| dependence(d)["by_motif_length"][0]["measured"] = json!("x"),
        ),
        case(
            "un mononucleótido que dice no tener curva medida",
            "una fila dice tener curva medida y no la tiene, o al revés",
            |d| dependence(d)["by_motif_length"][0]["measured"] = json!(false),
        ),
        case(
            "un trinucleótido que dice tener curva medida",
            "una fila dice tener curva medida y no la tiene, o al revés",
            |d| dependence(d)["by_motif_length"][2]["measured"] = json!(true),
        ),
        case(
            "una longitud de referencia que no es la del algoritmo",
            "una longitud de tracto de referencia no es la del algoritmo",
            |d| dependence(d)["by_motif_length"][1]["reference_tract_bases"] = json!(33),
        ),
        case(
            "una longitud de referencia donde no debe haber",
            "una longitud de tracto de referencia no es la del algoritmo",
            |d| dependence(d)["by_motif_length"][2]["reference_tract_bases"] = json!(30),
        ),
        case(
            "un tramo medido que no es el del algoritmo",
            "el tramo medido publicado no es el del algoritmo",
            |d| dependence(d)["by_motif_length"][0]["measured_units_min"] = json!(4),
        ),
        case(
            "un tramo medido en un motivo sin curva",
            "el tramo medido publicado no es el del algoritmo",
            |d| dependence(d)["by_motif_length"][2]["measured_units_max"] = json!(15),
        ),
        case(
            "una tabla de factores que no es una lista",
            "falta la tabla de factores",
            |d| dependence(d)["by_motif_length"][0]["factors"] = json!(5),
        ),
        case(
            "una tabla de factores que no es una lista, sin curva",
            "falta la tabla de factores",
            |d| dependence(d)["by_motif_length"][2]["factors"] = json!(5),
        ),
        case(
            "un factor de menos",
            "la tabla de factores no cubre el tramo medido",
            |d| {
                dependence(d)["by_motif_length"][0]["factors"]
                    .as_array_mut()
                    .expect("factores")
                    .pop();
            },
        ),
        case(
            "un factor de más en un motivo sin curva",
            "la tabla de factores no cubre el tramo medido",
            |d| {
                dependence(d)["by_motif_length"][2]["factors"] =
                    json!([{"units": 5, "factor": 1.0}]);
            },
        ),
        case(
            "una fila de factores que no es un objeto",
            "una fila de la tabla de factores no es objeto",
            |d| dependence(d)["by_motif_length"][0]["factors"][3] = json!(5),
        ),
        case(
            "un factor que no es un decimal",
            "falta un factor de la tabla",
            |d| dependence(d)["by_motif_length"][0]["factors"][3]["factor"] = json!("x"),
        ),
        case(
            "unas unidades de un factor que no son las del tramo",
            "un factor de la curva no es el del algoritmo",
            |d| dependence(d)["by_motif_length"][1]["factors"][3]["units"] = json!(9),
        ),
        case(
            "una fila de factores desplazada: unidades y factor coherentes entre sí pero no los esperados",
            "un factor de la curva no es el del algoritmo",
            |d| {
                // Las unidades 9 con **su** factor, en el lugar de las 8: el
                // factor es el correcto para las unidades publicadas, así que
                // solo las unidades delatan la fila.
                let shifted = dependence(d)["by_motif_length"][1]["factors"][4].clone();
                dependence(d)["by_motif_length"][1]["factors"][3] = shifted;
            },
        ),
        case(
            "un factor que difiere en el último decimal",
            "un factor de la curva no es el del algoritmo",
            |d| {
                let factor = dependence(d)["by_motif_length"][1]["factors"][3]["factor"]
                    .as_f64()
                    .expect("factor");
                dependence(d)["by_motif_length"][1]["factors"][3]["factor"] =
                    json!(next_up(factor));
            },
        ),
    ];
    for index in 0..6 {
        cases.push(case(
            format!("la fila {index} con la longitud de motivo cambiada"),
            "la tabla de longitudes no va en orden de longitud de motivo",
            move |d| {
                dependence(d)["by_motif_length"][index]["motif_length_bases"] = json!(index + 2)
            },
        ));
    }
    assert_technical_cases(&package, &cases);
}

/// **El tartamudeo clonal de la amplificación.**
///
/// La probabilidad por ciclo de una fila **derivada** tiene que salir, al bit,
/// de sus dos fracciones y de los ciclos publicados; una fila **declarada** no
/// se deriva y la procedencia solo admite esas dos palabras. Con el
/// multiplicador de los homopolímeros publicado, el bloque lo admite y sigue
/// cerrado.
#[test]
fn the_clonal_stutter_section_is_judged() {
    let package = publish_amplified("");
    let mut cases = vec![
        case(
            "un algoritmo clonal que no es el publicado",
            "algoritmo o fórmula del tartamudeo clonal inválidos",
            |d| clonal(d)["algorithm"] = json!("clonal-stutter-v1"),
        ),
        case(
            "una fórmula clonal que no es la publicada",
            "algoritmo o fórmula del tartamudeo clonal inválidos",
            |d| clonal(d)["formula"] = json!("r_pcr = r_with_pcr - r_without_pcr"),
        ),
        case(
            "unas tasas por ciclo que no se marcan como toscas",
            "las tasas por ciclo se publican marcadas como toscas",
            |d| clonal(d)["coarse"] = json!(false),
        ),
        case(
            "una marca de tosquedad clonal que no es un booleano",
            "las tasas por ciclo se publican marcadas como toscas",
            |d| clonal(d)["coarse"] = json!("true"),
        ),
        case(
            "unos ciclos que no caben en treinta y dos bits",
            "ciclos fuera de dominio",
            |d| clonal(d)["cycles"] = json!(u64::from(u32::MAX) + 1),
        ),
        case(
            "una tabla por ciclo que no es una lista",
            "falta la tabla de probabilidades por ciclo",
            |d| clonal(d)["per_cycle_slip"] = json!(5),
        ),
        case(
            "una fila de menos en la tabla por ciclo",
            "la tabla de probabilidades por ciclo no cubre el tramo calibrado",
            |d| {
                clonal(d)["per_cycle_slip"]
                    .as_array_mut()
                    .expect("tabla")
                    .pop();
            },
        ),
        case(
            "una fila de más en la tabla por ciclo",
            "la tabla de probabilidades por ciclo no cubre el tramo calibrado",
            |d| {
                let extra = clonal(d)["per_cycle_slip"][5].clone();
                clonal(d)["per_cycle_slip"]
                    .as_array_mut()
                    .expect("tabla")
                    .push(extra);
            },
        ),
        case(
            "una fila de la tabla por ciclo que no es un objeto",
            "una fila de la tabla por ciclo no es objeto",
            |d| clonal(d)["per_cycle_slip"][1] = json!(5),
        ),
        case(
            "dos filas de la tabla por ciclo intercambiadas",
            "la tabla por ciclo no va en orden de longitud de motivo",
            |d| {
                clonal(d)["per_cycle_slip"]
                    .as_array_mut()
                    .expect("tabla")
                    .swap(0, 1);
            },
        ),
        case(
            "una fracción con PCR que no es un decimal",
            "falta la fracción con PCR",
            |d| clonal(d)["per_cycle_slip"][0]["with_pcr_fraction"] = json!("x"),
        ),
        case(
            "una fracción sin PCR que no es un decimal",
            "falta la fracción sin PCR",
            |d| clonal(d)["per_cycle_slip"][0]["without_pcr_fraction"] = json!("x"),
        ),
        case(
            "una probabilidad por ciclo que no es un decimal",
            "falta la probabilidad por ciclo",
            |d| clonal(d)["per_cycle_slip"][0]["probability"] = json!("x"),
        ),
        case(
            "una probabilidad derivada que difiere en el último decimal",
            "una probabilidad por ciclo no sale de las fracciones publicadas",
            |d| {
                let probability = clonal(d)["per_cycle_slip"][1]["probability"]
                    .as_f64()
                    .expect("probabilidad");
                clonal(d)["per_cycle_slip"][1]["probability"] = json!(next_up(probability));
            },
        ),
        case(
            "una fracción con PCR cambiada sin recalcular la probabilidad",
            "una probabilidad por ciclo no sale de las fracciones publicadas",
            |d| clonal(d)["per_cycle_slip"][2]["with_pcr_fraction"] = json!(0.012),
        ),
        case(
            "unos ciclos cambiados sin recalcular las probabilidades",
            "una probabilidad por ciclo no sale de las fracciones publicadas",
            |d| clonal(d)["cycles"] = json!(11),
        ),
        case(
            "una procedencia de probabilidad que no es ni derivada ni declarada",
            "la procedencia de una probabilidad por ciclo no es ni derivada ni declarada",
            |d| clonal(d)["per_cycle_slip"][0]["source"] = json!("measured"),
        ),
        case(
            "un campo ajeno junto al multiplicador de los homopolímeros",
            "campos ausentes o desconocidos en technical_event.artifacts.pcr.stutter",
            |d| clonal(d)["sobra"] = json!(1),
        ),
    ];
    // Una fila marcada como derivada pero con una fracción hecha a medida: la
    // procedencia tiene que ser una de las dos palabras, en cada fila.
    for index in 0..6 {
        cases.push(case(
            format!("la fila {index} con una procedencia inventada"),
            "la procedencia de una probabilidad por ciclo no es ni derivada ni declarada",
            move |d| clonal(d)["per_cycle_slip"][index]["source"] = json!("inventada"),
        ));
    }
    assert_technical_cases(&package, &cases);
    assert_technical_accepted(
        &package,
        &[
            accepted(
                "una fila declarada no se deriva: su probabilidad es la que diga",
                |d| {
                    clonal(d)["per_cycle_slip"][0]["source"] = json!("declared");
                    clonal(d)["per_cycle_slip"][0]["probability"] = json!(0.5);
                },
            ),
            accepted(
                "el multiplicador de los homopolímeros, que es opcional, se admite \
                 aunque su valor no se juzgue en este paso",
                |d| clonal(d)["homopolymer_pcr_fold"] = json!(2.5),
            ),
            accepted(
                "unos ciclos de treinta y dos bits justos, con todas las filas declaradas",
                |d| {
                    clonal(d)["cycles"] = json!(u64::from(u32::MAX));
                    for index in 0..6 {
                        clonal(d)["per_cycle_slip"][index]["source"] = json!("declared");
                    }
                },
            ),
        ],
    );
}

/// **Con el multiplicador de los homopolímeros declarado, el bloque clonal lo
/// admite y sigue cerrado.**
///
/// El campo es opcional: sin él, el bloque no lo lleva (el resto de la tabla
/// cubre ese caso). Con él, la lista de campos admitidos crece en ese y solo
/// en ese, de modo que un campo ajeno se rechaza igual.
#[test]
fn the_clonal_block_admits_the_declared_homopolymer_fold_and_stays_closed() {
    let package = publish_amplified("homopolymer_pcr_fold = 2.5\n");
    let original = package.json("effective-config.json");
    assert_eq!(
        original["technical_event"]["artifacts"]["pcr"]["stutter"]["homopolymer_pcr_fold"], 2.5,
        "el fixture publica el multiplicador"
    );
    assert_technical_cases(
        &package,
        &[
            case(
                "un campo ajeno junto al multiplicador",
                "campos ausentes o desconocidos en technical_event.artifacts.pcr.stutter",
                |d| clonal(d)["sobra"] = json!(1),
            ),
            case(
                "el multiplicador cambiado de nombre",
                "campos ausentes o desconocidos en technical_event.artifacts.pcr.stutter",
                |d| {
                    let object = clonal(d).as_object_mut().expect("bloque");
                    let value = object.remove("homopolymer_pcr_fold").expect("campo");
                    object.insert("multiplicador".to_owned(), value);
                },
            ),
        ],
    );
}

/// **La amplificación: catálogo, inversión de la fracción de duplicados,
/// fracción esperada, clases de eficiencia y curva de composición.**
///
/// Con microsatélites la amplificación abre el tartamudeo clonal y publica el
/// catálogo v6; sin ellos, el v5. El catálogo se juzga contra **esa
/// derivación** en los dos sentidos.
#[test]
fn the_effective_pcr_section_is_judged() {
    let package = publish_amplified("");
    let original = package.json("effective-config.json");
    assert_eq!(
        original["technical_event"]["artifacts"]["pcr"]["rng_catalog"], "f4-technical-v6",
        "con microsatélites, la amplificación publica el catálogo v6"
    );
    let mut cases = vec![
        case(
            "un algoritmo de amplificación que no es el publicado",
            "algoritmo de amplificación inválido",
            |d| pcr(d)["algorithm"] = json!("pcr-amplification-v0"),
        ),
        case(
            "un catálogo v5 en una amplificación que abre el tartamudeo clonal",
            "algoritmo de amplificación inválido",
            |d| pcr(d)["rng_catalog"] = json!("f4-technical-v5"),
        ),
        case(
            "un método de inversión que no es el publicado",
            "método de inversión de la fracción de duplicados inválido",
            |d| pcr(d)["duplicate_fraction_inversion"]["method"] = json!("newton"),
        ),
        case(
            "una fracción de duplicados declarada que no es un decimal",
            "falta la fracción de duplicados declarada",
            |d| pcr(d)["declared_duplicate_fraction"] = json!("x"),
        ),
        case(
            "una m derivada que no es un decimal",
            "falta la m derivada",
            |d| pcr(d)["mean_pairs_per_molecule"] = json!("x"),
        ),
        case(
            "una iteración de menos en la inversión",
            "la inversión publicada no corresponde a la fracción declarada",
            |d| pcr(d)["duplicate_fraction_inversion"]["iterations"] = json!(63),
        ),
        case(
            "una iteración de más en la inversión",
            "la inversión publicada no corresponde a la fracción declarada",
            |d| pcr(d)["duplicate_fraction_inversion"]["iterations"] = json!(65),
        ),
        case(
            "ninguna iteración con una fracción declarada distinta de cero",
            "la inversión publicada no corresponde a la fracción declarada",
            |d| pcr(d)["duplicate_fraction_inversion"]["iterations"] = json!(0),
        ),
        case(
            "una fracción esperada que no es un decimal",
            "falta la fracción esperada",
            |d| pcr(d)["expected_duplicate_fraction"] = json!("x"),
        ),
        case(
            "una fracción esperada que difiere en el último decimal",
            "la fracción esperada publicada no sale de la m publicada",
            |d| {
                let expected = pcr(d)["expected_duplicate_fraction"]
                    .as_f64()
                    .expect("esperada");
                pcr(d)["expected_duplicate_fraction"] = json!(next_up(expected));
            },
        ),
        case(
            "una m cambiada sin recalcular la fracción esperada",
            "la fracción esperada publicada no sale de la m publicada",
            |d| pcr(d)["mean_pairs_per_molecule"] = json!(0.8),
        ),
        case(
            "unas clases de eficiencia que no son una lista",
            "faltan las clases de eficiencia",
            |d| pcr(d)["efficiency_classes"] = json!(5),
        ),
        case(
            "ninguna clase de eficiencia",
            "las clases de eficiencia publicadas están fuera del contrato",
            |d| pcr(d)["efficiency_classes"] = json!([]),
        ),
        case(
            "una clase de eficiencia más que el máximo",
            "las clases de eficiencia publicadas están fuera del contrato",
            |d| {
                let class = pcr(d)["efficiency_classes"][0].clone();
                pcr(d)["efficiency_classes"] = json!([class, class, class, class, class]);
            },
        ),
        case(
            "una curva de composición que no es una lista",
            "falta la curva de composición",
            |d| pcr(d)["gc_relative_efficiency"] = json!(5),
        ),
        case(
            "una curva de composición con un solo nudo",
            "la curva de composición publicada necesita al menos dos nudos",
            |d| {
                pcr(d)["gc_relative_efficiency"]
                    .as_array_mut()
                    .expect("curva")
                    .pop();
            },
        ),
        case(
            "una curva de composición sin nudos",
            "la curva de composición publicada necesita al menos dos nudos",
            |d| pcr(d)["gc_relative_efficiency"] = json!([]),
        ),
    ];
    cases.push(case(
        "un campo ajeno en la amplificación",
        "campos ausentes o desconocidos en technical_event.artifacts.pcr",
        |d| pcr(d)["sobra"] = json!(1),
    ));
    assert_technical_cases(&package, &cases);
    assert_technical_accepted(
        &package,
        &[
            accepted("cuatro clases de eficiencia, el máximo", |d| {
                let class = pcr(d)["efficiency_classes"][0].clone();
                pcr(d)["efficiency_classes"] = json!([class, class, class, class]);
            }),
            accepted("una curva de composición de tres nudos", |d| {
                let knot = pcr(d)["gc_relative_efficiency"][1].clone();
                pcr(d)["gc_relative_efficiency"]
                    .as_array_mut()
                    .expect("curva")
                    .push(knot);
            }),
        ],
    );
}

/// **Sin microsatélites la amplificación publica el catálogo v5 y no lleva
/// sección de tartamudeo clonal.**
///
/// Es el otro sentido de la derivación del catálogo: el v6 sin la sección se
/// rechaza igual que la sección sin el v6, y una sección clonal añadida a mano
/// obliga al v6.
#[test]
fn the_pcr_catalog_follows_whether_the_clonal_section_is_present() {
    let package = publish_amplified_without_loci("0.3");
    let original = package.json("effective-config.json");
    assert_eq!(
        original["technical_event"]["artifacts"]["pcr"]["rng_catalog"], "f4-technical-v5",
        "sin microsatélites, la amplificación publica el catálogo v5"
    );
    assert!(
        original["technical_event"]["artifacts"]["pcr"]
            .get("stutter")
            .is_none(),
        "sin microsatélites no hay tartamudeo clonal"
    );
    assert_technical_cases(
        &package,
        &[
            case(
                "un catálogo v6 sin tartamudeo clonal",
                "algoritmo de amplificación inválido",
                |d| pcr(d)["rng_catalog"] = json!("f4-technical-v6"),
            ),
            case(
                "un catálogo v5 con un tartamudeo clonal añadido",
                "algoritmo de amplificación inválido",
                |d| {
                    pcr(d)["stutter"] = json!({});
                },
            ),
            case(
                "un algoritmo de amplificación que no es el publicado",
                "algoritmo de amplificación inválido",
                |d| pcr(d)["algorithm"] = json!("otro"),
            ),
        ],
    );
}

/// **Con fracción de duplicados cero la inversión no itera y la `m` es cero.**
///
/// Solo juzga la estructura del evento técnico: que la `m` publicada sea la de
/// la fracción declarada lo ata la comparación con el TOML, que corre después
/// y que este fichero no ejercita.
///
/// La inversión devuelve cero por retorno anticipado: ni una iteración, ni un
/// residuo. Cada una de las dos mitades de esa condición se rechaza por
/// separado (una iteración publicada; una `m` distinta de cero con cero
/// iteraciones), y el lado contrario —una fracción positiva mínima exige las
/// 64 iteraciones— se acepta.
#[test]
fn a_zero_duplicate_fraction_publishes_no_iterations_and_a_zero_m() {
    let package = publish_amplified_without_loci("0.0");
    let original = package.json("effective-config.json");
    let published = &original["technical_event"]["artifacts"]["pcr"];
    assert_eq!(published["declared_duplicate_fraction"], 0.0);
    assert_eq!(published["duplicate_fraction_inversion"]["iterations"], 0);
    assert_eq!(published["mean_pairs_per_molecule"], 0.0);
    assert_technical_cases(
        &package,
        &[
            case(
                "una iteración publicada con fracción cero",
                "la inversión publicada no corresponde a la fracción declarada",
                |d| pcr(d)["duplicate_fraction_inversion"]["iterations"] = json!(64),
            ),
            case(
                "una fracción negativa con cero iteraciones",
                "la inversión publicada no corresponde a la fracción declarada",
                |d| pcr(d)["declared_duplicate_fraction"] = json!(-0.5),
            ),
            case(
                "una m negativa con fracción cero y cero iteraciones",
                "la inversión publicada no corresponde a la fracción declarada",
                |d| pcr(d)["mean_pairs_per_molecule"] = json!(-0.5),
            ),
            case(
                "una m distinta de cero con fracción cero y cero iteraciones",
                "la inversión publicada no corresponde a la fracción declarada",
                |d| pcr(d)["mean_pairs_per_molecule"] = json!(0.5),
            ),
        ],
    );
    assert_technical_accepted(
        &package,
        &[accepted(
            "una fracción positiva mínima con las 64 iteraciones",
            |d| {
                pcr(d)["declared_duplicate_fraction"] = json!(f64::MIN_POSITIVE);
                pcr(d)["duplicate_fraction_inversion"]["iterations"] = json!(64);
            },
        )],
    );
}

/// **El juez del abanico reproduce lo que ocurre en los bordes del modelo.**
///
/// Cuatro loci cuyos linajes se detienen —al agotar la tirada, contra el tope,
/// por una expansión que ya no cabe— y uno interrumpido, cuyos alelos ordenados
/// por longitud no lo están por sus bases. El control es que el paquete
/// publicado por el programa **se acepta**: cualquier regla del recorrido
/// (viabilidad, sentido del paso, recuento de pasos aplicados y no aplicados,
/// orden canónico) que el juez escribiera distinta de la del productor haría
/// divergir el abanico y rechazaría un paquete bueno.
#[test]
fn the_fan_judge_reproduces_the_edges_of_the_model() {
    let package = publish_edges();
    let effective = package.json("effective-config.json");
    judge_resources(&package, &effective).expect("control positivo: el paquete intacto se acepta");
    for index in [0, 1, 3] {
        let locus = &effective["repeat_loci"][index];
        assert!(
            locus["slippage"]["unapplied_events"]
                .as_u64()
                .expect("no aplicados")
                > 0,
            "el fixture: {} llega a un borde",
            locus["id"]
        );
    }
    // El tracto interrumpido: el alelo más largo no es el último por sus bases.
    let interrupted = effective["repeat_loci"][2]["alleles"]
        .as_array()
        .expect("alelos");
    let sequences: Vec<&str> = interrupted
        .iter()
        .map(|allele| allele["sequence"].as_str().expect("secuencia"))
        .collect();
    let mut by_bases = sequences.clone();
    by_bases.sort_unstable();
    assert_ne!(
        sequences, by_bases,
        "el fixture: el orden por longitud no es el orden por bases"
    );
}

/// **Los linajes son individuos por ploidía.** Con una población tetraploide,
/// el juez tiene que sumar cuatro copias por individuo y no dos.
#[test]
fn the_fan_judge_counts_the_ploidy_of_the_populations() {
    let package = publish_mixed_with_ploidy(3, 4);
    let effective = package.json("effective-config.json");
    assert_eq!(effective["repeat_loci"][1]["slippage"]["lineages"], 12);
    judge_resources(&package, &effective).expect("tres individuos tetraploides se aceptan");
}

/// **Con inserción, el abanico se siembra con la referencia final.**
///
/// El juez lee el identificador de `inserted_reference` y no el de la
/// referencia suministrada; con el otro, el abanico rehecho es otro y un
/// paquete bueno se rechazaría.
#[test]
fn the_fan_judge_seeds_the_fan_with_the_final_reference_when_there_is_an_insertion() {
    let package = publish_inserting();
    let effective = package.json("effective-config.json");
    assert_eq!(effective["inserted_reference"]["id"], "ref-main-str");
    judge_resources(&package, &effective).expect("control positivo: el paquete intacto se acepta");
}

/// **Un paquete manipulado sin poblaciones ni matriz de SNV sí llega a las
/// bases no materializadas de la referencia.**
///
/// La referencia se lee con las bases solo si el documento declara matriz,
/// poblaciones, evento técnico, VCF o eventos directos
/// (`effective_needs_reference_bases`). Quitar todo eso de un paquete
/// publicado de verdad exige rehacer a la vez la copia del TOML, su huella y
/// el identificador de ejecución (que depende de ella), la matriz y las
/// poblaciones de la efectiva, la matriz del manifiesto, sus dos roles y sus
/// artefactos: con todo coherente, los demás jueces no encuentran nada que
/// rechazar y el del abanico llega a cortar el tramo ancestral de una
/// referencia sin bases. Cada paso previo se rechaza por su propio motivo.
#[test]
fn a_package_forged_without_populations_reaches_the_unmaterialized_reference() {
    let package = publish_mixed(5);
    let toml_path = package.destination.join("inputs/config/general.toml");
    let toml = fs::read_to_string(&toml_path).expect("copia del TOML");
    let mut forged = String::new();
    let mut skip = false;
    for line in toml.lines() {
        if line.starts_with('[') {
            skip = line.contains("population_snvs") || line.contains("initial_populations");
        }
        if !skip {
            forged.push_str(line);
            forged.push('\n');
        }
    }
    assert_ne!(forged, toml, "la copia del TOML perdió las tablas");
    fs::write(&toml_path, &forged).expect("TOML manipulado");
    let digest = genome_population_simulator::integrity::configuration_digest(forged.as_bytes())
        .expect("huella");
    let mut effective = package.json("effective-config.json");
    let mut manifest = package.json("manifest.json");
    let root = genome_population_simulator::randomness::RootSeed::parse(
        "effective-config.seed.root",
        effective["seed"]["root"].as_str().expect("semilla"),
    )
    .expect("semilla");
    let execution = genome_population_simulator::identity::derive_execution_id(root, &digest)
        .expect("identificador de ejecución");
    manifest["execution_id"] = json!(execution.as_str());
    for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
        if artifact["id"] == "config-general" {
            artifact["checksum"]["digest"] = json!(digest.clone());
            artifact["size_bytes"] = json!(forged.len());
        }
    }
    effective["initial_populations"] = json!([]);
    effective["population_snvs"] = Value::Null;
    manifest["population_snvs"] = Value::Null;
    let judge = |effective: &Value, manifest: &Value| {
        validate_published_resources(
            &package.destination,
            &package.effective_path(),
            &serde_json::to_vec_pretty(effective).expect("JSON"),
            &package.manifest_path(),
            &serde_json::to_vec_pretty(manifest).expect("JSON"),
            None,
        )
    };
    // Primero el paquete sin matriz en el manifiesto pero con sus roles: otro
    // juez lo rechaza antes, de modo que el intento no basta.
    assert_ne!(
        detail(judge(&effective, &manifest)),
        fan("las bases de la referencia publicada no están materializadas"),
        "con los roles de la matriz todavía en el manifiesto, gana otro rechazo"
    );
    let roles = ["population-snvs-tsv", "population-snvs-metadata"];
    manifest["output_categories"]["minimum"]
        .as_array_mut()
        .expect("mínimos")
        .retain(|role| !roles.contains(&role.as_str().expect("rol")));
    manifest["artifacts"]
        .as_array_mut()
        .expect("artefactos")
        .retain(|artifact| !roles.contains(&artifact["id"].as_str().expect("id")));
    assert_eq!(
        detail(judge(&effective, &manifest)),
        fan("las bases de la referencia publicada no están materializadas")
    );
}

// ---------------------------------------------------------------------------
// DT-051: los loci publicados, atados al TOML copiado y a la verdad
// ---------------------------------------------------------------------------

const EFFECTIVE_LOCI: &str =
    "los loci de repetición efectivos no son los que el TOML y la referencia producen";
const EFFECTIVE_ADMISSION: &str = "los mínimos de admisión efectivos no son los del TOML";
const MANIFEST_SECTION: &str =
    "la sección de loci del manifiesto no es la que el TOML y la referencia producen";
const REALIZED_LEFT: &str =
    "las copias realizadas de los alelos de STR-LEFT no son las que la verdad lleva";
const REALIZED_RIGHT: &str =
    "las copias realizadas de los alelos de STR-RIGHT no son las que la verdad lleva";

/// El juez completo con la configuración efectiva y el manifiesto que se den.
fn judge_both(package: &Package, effective: &Value, manifest: &Value) -> Result<(), Diagnostic> {
    validate_published_resources(
        &package.destination,
        &package.effective_path(),
        &serde_json::to_vec_pretty(effective).expect("JSON"),
        &package.manifest_path(),
        &serde_json::to_vec_pretty(manifest).expect("JSON"),
        None,
    )
}

/// El juez completo con el manifiesto alterado y la efectiva intacta.
fn assert_manifest_cases_against_the_toml(package: &Package, cases: &[Case]) {
    let effective = package.json("effective-config.json");
    let original = package.json("manifest.json");
    judge_both(package, &effective, &original)
        .expect("control positivo: el paquete intacto se acepta");
    for case in cases {
        let document = altered(&original, &case.label, case.mutate.as_ref());
        assert_eq!(
            detail(judge_both(package, &effective, &document)),
            case.expected,
            "«{}» se rechazó por otra razón",
            case.label
        );
    }
}

/// **Cada campo de los loci de la configuración efectiva, y sus mínimos de
/// admisión, tienen que ser los que el TOML produce** (DT-051).
///
/// Son las dieciocho alteraciones que la revisión adversarial de la tanda 6
/// reprodujo, más la tirada mínima de admisión (la segunda revisión de DT-051
/// comprobó que comparar solo la pureza sobrevivía a toda la batería): hasta el 10/10/2026 el juez completo las **aceptaba** todas
/// (comprobado ejecutando contra el código anterior). En el locus generado,
/// las que tocan lo que el juez del abanico rehace (la media y el sesgo
/// cambiados en un ulp, las unidades) no cambian el abanico rehecho, así que
/// ese juez no las ve y llegan a la comparación con el TOML.
#[test]
fn the_effective_repeat_loci_are_the_ones_the_toml_produces() {
    let package = publish_mixed(5);
    let declared_sequence = "GTGTAAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT";
    assert_eq!(declared_sequence.len(), INTERRUPTED_ALLELE.len());
    assert_effective_cases(
        &package,
        &[
            case("alelo declarado: secuencia", EFFECTIVE_LOCI, move |d| {
                d["repeat_loci"][0]["alleles"][0]["sequence"] = json!(declared_sequence);
            }),
            case(
                "alelo declarado: frecuencia declarada",
                EFFECTIVE_LOCI,
                |d| {
                    d["repeat_loci"][0]["alleles"][0]["declared_frequency"] = json!(0.5);
                },
            ),
            case(
                "alelo declarado: frecuencia efectiva",
                EFFECTIVE_LOCI,
                |d| {
                    d["repeat_loci"][0]["alleles"][0]["effective_frequency"] = json!(0.5);
                },
            ),
            case("alelo declarado: pureza", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][0]["alleles"][0]["purity"] = json!(0.5);
            }),
            case("alelo declarado: unidades", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][0]["alleles"][0]["units"] = json!(30.0);
            }),
            case("locus declarado: inicio", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][0]["start"] = json!(52);
            }),
            case("locus declarado: motivos", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][0]["motifs"] = json!(["TG"]);
            }),
            case("locus declarado: tope", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][0]["max_units"] = json!(3);
            }),
            case("locus declarado: pureza mínima", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][0]["admission"]["min_purity"] = json!(0.1);
            }),
            case("locus declarado: identificador", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][0]["id"] = json!("STR-FAKE");
            }),
            case(
                "locus generado: tope que el abanico no alcanza",
                EFFECTIVE_LOCI,
                |d| d["repeat_loci"][1]["max_units"] = json!(60),
            ),
            case("locus generado: media un ulp más", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][1]["slippage"]["mean_events"] = json!(next_up(3.0));
            }),
            case("locus generado: sesgo un ulp más", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][1]["slippage"]["expansion_bias"] = json!(next_up(0.5));
            }),
            case("alelo generado: pureza", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][1]["alleles"][0]["purity"] = json!(0.9);
            }),
            case("alelo generado: frecuencia efectiva", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][1]["alleles"][0]["effective_frequency"] = json!(0.2);
            }),
            case("alelo generado: unidades", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][1]["alleles"][0]["units"] = json!(46.0);
            }),
            case("alelo generado: tirada más larga", EFFECTIVE_LOCI, |d| {
                d["repeat_loci"][1]["alleles"][0]["longest_run_units"] = json!(3);
            }),
            case("mínimos de admisión: pureza", EFFECTIVE_ADMISSION, |d| {
                d["repeat_admission"]["min_purity"] = json!(0.1);
            }),
            case("mínimos de admisión: tirada", EFFECTIVE_ADMISSION, |d| {
                let run = d["repeat_admission"]["min_longest_run_units"]
                    .as_u64()
                    .expect("tirada mínima");
                d["repeat_admission"]["min_longest_run_units"] = json!(run + 1);
            }),
        ],
    );
}

/// **La sección de loci del manifiesto general tiene que ser la que el TOML
/// produce** (DT-051), salvo lo que sale del sorteo, que se juzga contra la
/// verdad en la prueba siguiente.
///
/// Son las alteraciones del manifiesto que la revisión adversarial de la tanda
/// 6 reprodujo, todas **aceptadas** hasta el 10/10/2026 (comprobado
/// ejecutando contra el código anterior): la procedencia de la referencia, los
/// mínimos de omisión, campos de un alelo y de un locus declarados, la
/// frecuencia declarada y los parámetros de un generado, y quitar un alelo
/// generado sin copias dándole su frecuencia a otro (punto (c): el manifiesto y
/// la efectiva no se comparaban entre sí).
#[test]
fn the_repeat_section_of_the_manifest_is_the_one_the_toml_produces() {
    let package = publish_mixed(5);
    let declared_sequence = "GTGTAAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT";
    assert_manifest_cases_against_the_toml(
        &package,
        &[
            case("procedencia de la referencia", MANIFEST_SECTION, |d| {
                section(d)["reference_source"] = json!("simulated");
            }),
            case("mínimos de omisión: pureza", MANIFEST_SECTION, |d| {
                section(d)["admission_defaults"]["min_purity"] = json!(0.1);
            }),
            case("alelo declarado: secuencia", MANIFEST_SECTION, move |d| {
                allele0(d)["sequence"] = json!(declared_sequence);
            }),
            case("alelo declarado: pureza", MANIFEST_SECTION, |d| {
                allele0(d)["purity"] = json!(0.5);
            }),
            case(
                "alelo declarado: frecuencia declarada",
                MANIFEST_SECTION,
                |d| {
                    allele0(d)["declared_frequency"] = json!(0.9);
                },
            ),
            case(
                "alelos generados: frecuencias declaradas intercambiadas",
                MANIFEST_SECTION,
                |d| {
                    let alleles = &mut generated(d)["alleles"];
                    let first = alleles[0]["declared_frequency"].clone();
                    let fourth = alleles[3]["declared_frequency"].clone();
                    assert_ne!(first, fourth, "el fixture: frecuencias distintas");
                    alleles[0]["declared_frequency"] = fourth;
                    alleles[3]["declared_frequency"] = first;
                },
            ),
            case("locus generado: media", MANIFEST_SECTION, |d| {
                generated(d)["slippage"]["mean_events"] = json!(1.0);
            }),
            case("locus generado: tope", MANIFEST_SECTION, |d| {
                generated(d)["max_units"] = json!(1000);
            }),
            case("locus declarado: inicio", MANIFEST_SECTION, |d| {
                declared(d)["start"] = json!(1);
            }),
            case("locus declarado: identificador", MANIFEST_SECTION, |d| {
                declared(d)["id"] = json!("STR-FAKE");
            }),
            case("alelo generado: otra secuencia", MANIFEST_SECTION, |d| {
                let allele = &mut generated(d)["alleles"][0];
                allele["sequence"] = json!("CA".repeat(46));
                allele["units"] = json!(46.0);
                allele["block_bases"] = json!([92]);
                allele["block_unit_counts"] = json!([46.0]);
                allele["longest_run_units"] = json!(46);
            }),
            case(
                "alelo generado sin copias quitado y su frecuencia dada a otro",
                MANIFEST_SECTION,
                |d| {
                    let alleles = generated(d)["alleles"].as_array_mut().expect("alelos");
                    assert_eq!(alleles[2]["realized_copies"], 0, "el fixture");
                    alleles.remove(2);
                    alleles[1]["declared_frequency"] = json!(0.6);
                    alleles[1]["effective_frequency"] = json!(0.6);
                },
            ),
        ],
    );
}

/// **Las copias realizadas de cada alelo son las que la verdad lleva**
/// (DT-051, punto (b)).
///
/// Intercambiar las copias realizadas (y con ellas la frecuencia realizada)
/// entre dos alelos conserva la suma y la aritmética que el lector del
/// manifiesto comprueba; hasta el 10/10/2026 el juez completo lo **aceptaba**
/// (comprobado ejecutando contra el código anterior). El tercer caso mueve la
/// copia que lleva el alelo **igual a la referencia** (cien bases de `CA`, sin
/// caso en la verdad) a otro alelo: discrimina que una copia sin caso cuente
/// como el alelo de la referencia.
#[test]
fn the_realized_copies_are_the_ones_the_truth_carries() {
    let package = publish_mixed(5);
    let manifest = package.json("manifest.json");
    let copies = |locus: usize| -> Vec<u64> {
        manifest["repeat_loci"]["loci"][locus]["alleles"]
            .as_array()
            .expect("alelos")
            .iter()
            .map(|allele| allele["realized_copies"].as_u64().expect("copias"))
            .collect()
    };
    assert_eq!(copies(0), [4, 6], "el fixture del locus declarado");
    assert_eq!(copies(1), [0, 7, 0, 1, 2], "el fixture del locus generado");
    assert_eq!(
        manifest["repeat_loci"]["loci"][1]["alleles"][3]["sequence"],
        json!("CA".repeat(50)),
        "el fixture: el cuarto alelo generado es el tracto de la referencia"
    );
    let set = |allele: &mut Value, copies: u64, drawn: u64| {
        allele["realized_copies"] = json!(copies);
        allele["realized_frequency"] = json!(copies as f64 / drawn as f64);
    };
    assert_manifest_cases_against_the_toml(
        &package,
        &[
            case(
                "alelos declarados: copias intercambiadas",
                REALIZED_LEFT,
                move |d| {
                    set(allele0(d), 6, 10);
                    set(allele1(d), 4, 10);
                },
            ),
            case(
                "alelos generados: copias intercambiadas",
                REALIZED_RIGHT,
                move |d| {
                    set(&mut generated(d)["alleles"][1], 0, 10);
                    set(&mut generated(d)["alleles"][2], 7, 10);
                },
            ),
            case(
                "la copia del alelo de la referencia, dada a otro",
                REALIZED_RIGHT,
                move |d| {
                    set(&mut generated(d)["alleles"][2], 1, 10);
                    set(&mut generated(d)["alleles"][3], 0, 10);
                },
            ),
        ],
    );
}

/// **El trasplante de un abanico con sus parámetros se rechaza** (DT-051,
/// punto (a)).
///
/// Se publica el mismo escenario con `mean_events = 1.0` y su locus generado
/// (parámetros y abanico) se lleva al paquete original (3.0). El juez del
/// abanico lo rehace desde los parámetros **publicados** y lo acepta, como
/// debe; hasta el 10/10/2026 nada más lo miraba y el juez completo aceptaba
/// los tres trasplantes (comprobado ejecutando contra el código anterior).
/// Ahora la comparación con el TOML, que corre después, los rechaza.
#[test]
fn a_fan_transplanted_with_its_parameters_is_rejected() {
    let original = publish_mixed(5);
    let donor = publish_mixed_with(5, 2, "1.0");
    let effective = original.json("effective-config.json");
    let manifest = original.json("manifest.json");
    let mut forged_effective = effective.clone();
    forged_effective["repeat_loci"][1] =
        donor.json("effective-config.json")["repeat_loci"][1].clone();
    let mut forged_manifest = manifest.clone();
    forged_manifest["repeat_loci"]["loci"][1] =
        donor.json("manifest.json")["repeat_loci"]["loci"][1].clone();
    assert_ne!(
        forged_effective, effective,
        "el trasplante alteró la efectiva"
    );
    assert_ne!(
        forged_manifest, manifest,
        "el trasplante alteró el manifiesto"
    );
    assert_eq!(
        forged_effective["repeat_loci"][1]["slippage"]["mean_events"],
        json!(1.0)
    );
    judge_both(&original, &effective, &manifest).expect("control positivo");
    judge_both(
        &donor,
        &donor.json("effective-config.json"),
        &donor.json("manifest.json"),
    )
    .expect("control positivo del donante");
    assert_eq!(
        detail(judge_both(&original, &forged_effective, &manifest)),
        EFFECTIVE_LOCI,
        "solo la efectiva"
    );
    assert_eq!(
        detail(judge_both(&original, &effective, &forged_manifest)),
        MANIFEST_SECTION,
        "solo el manifiesto"
    );
    assert_eq!(
        detail(judge_both(&original, &forged_effective, &forged_manifest)),
        EFFECTIVE_LOCI,
        "los dos"
    );
}

/// Un locus declarado con **dos alelos iguales** (`(GT)24` dos veces, mitad y
/// mitad) y uno distinto, sobre la referencia del escenario `mixed`.
fn publish_twin_alleles() -> Package {
    let reference = mixed_reference();
    let length = reference.trim_end().lines().nth(1).expect("bases").len();
    let config = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
outside_weight = 0.0
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.regions]]
contig = "chr1"
start = 1
end = 50
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

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {length}

[[reference.repeat_loci]]
id = "STR-LEFT"
sequence = "chr1"
start = 51
end = 100
motif = "GT"

[[reference.repeat_loci.alleles]]
units = 24
frequency = 0.4

[[reference.repeat_loci.alleles]]
units = 24
frequency = 0.4

[[reference.repeat_loci.alleles]]
units = 26
frequency = 0.2

[[initial_populations]]
population_id = "P1"
individual_count = 6
individual_ids = ["I1", "I2", "I3", "I4", "I5", "I6"]
ploidy = 2
genotype_model = "hwe"
"#
    );
    publish(&[("general.toml", &config), ("reference.fa", &reference)])
}

/// **Dos alelos que la verdad no distingue se comparan por su suma** (DT-051,
/// punto (b)).
///
/// La verdad publica de cada copia las bases y los bloques del tracto, no el
/// índice del alelo que le tocó: con dos alelos iguales, repartir sus copias
/// de otra forma es indistinguible, y el juez no puede exigir más de lo que la
/// verdad dice. Lado aceptado: el reparto cambiado se acepta. Lado rechazado:
/// darle una copia de los gemelos al alelo distinto se rechaza.
#[test]
fn alleles_the_truth_cannot_tell_apart_are_compared_by_their_sum() {
    let package = publish_twin_alleles();
    let effective = package.json("effective-config.json");
    let manifest = package.json("manifest.json");
    judge_both(&package, &effective, &manifest).expect("control positivo");
    let alleles = &manifest["repeat_loci"]["loci"][0]["alleles"];
    let copies: [u64; 3] =
        std::array::from_fn(|index| alleles[index]["realized_copies"].as_u64().expect("copias"));
    assert_eq!(
        alleles[0]["sequence"], alleles[1]["sequence"],
        "el fixture: gemelos"
    );
    assert!(
        copies[0] > 0 && copies[1] > 0,
        "el fixture: los dos gemelos tienen copias: {copies:?}"
    );
    let set = |allele: &mut Value, copies: u64| {
        allele["realized_copies"] = json!(copies);
        allele["realized_frequency"] = json!(copies as f64 / 12.0);
    };
    let shifted = altered(&manifest, "reparto entre gemelos", &move |d| {
        set(allele0(d), copies[0] - 1);
        set(allele1(d), copies[1] + 1);
    });
    judge_both(&package, &effective, &shifted)
        .expect("un reparto distinto entre alelos indistinguibles se acepta");
    let leaked = altered(
        &manifest,
        "copia de un gemelo dada al distinto",
        &move |d| {
            set(allele0(d), copies[0] - 1);
            set(
                &mut d["repeat_loci"]["loci"][0]["alleles"][2],
                copies[2] + 1,
            );
        },
    );
    assert_eq!(
        detail(judge_both(&package, &effective, &leaked)),
        REALIZED_LEFT
    );
}

// ---------------------------------------------------------------------------
// DT-052: la contabilidad de bases de tracto, juzgada al publicar y atada al TOML
// ---------------------------------------------------------------------------

const BOUND_MISMATCH: &str = "la cota del manifiesto no es la de la configuración efectiva";
const REALIZED_TRACT: &str =
    "las bases de tracto realizadas no son las que los documentos publicados llevan";
const MANIFEST_TRACT: &str =
    "la contabilidad de tracto del manifiesto no es la que resuelven el TOML y lo medido";
const EFFECTIVE_TRACT: &str =
    "la cota de bases de tracto efectiva no es la que resuelven el TOML y lo que el VCF lleva";

/// Suma `delta` a un entero de un documento (con signo, para poder rebajar).
fn shift(value: &mut Value, delta: i64) {
    let current = value.as_u64().expect("un entero");
    *value = json!(current.checked_add_signed(delta).expect("sin desbordar"));
}

/// Una alteración de la contabilidad de tracto: la etiqueta, el motivo
/// exacto, y lo que se cambia en el manifiesto y en la configuración
/// efectiva (cualquiera de los dos puede quedar intacto).
type TractCase = (&'static str, &'static str, fn(&mut Value), fn(&mut Value));

fn untouched(_: &mut Value) {}

/// **La contabilidad de bases de tracto se juzga al publicar, y su cota se ata
/// al TOML copiado** (DT-052).
///
/// Son las alteraciones que la segunda revisión adversarial de DT-051
/// reprodujo, más la versión coherente en los dos documentos y la del tope
/// declarado. Hasta el 10/10/2026 el juez completo las **aceptaba** todas: el
/// juez que recuenta las bases (`validate_published_tract_output_accounting`)
/// no lo llamaba nadie al publicar, y nada ataba la cota al TOML. Se comprobó
/// ejecutando esta prueba contra una copia del código anterior.
///
/// Cada caso altera las cifras de forma que la aritmética que el lector del
/// manifiesto comprueba (cada total es la suma de sus términos, lo realizado
/// no pasa de su cota) siga cuadrando: si no, el rechazo vendría de ese lector
/// y la prueba no diría nada del juez nuevo. Se ejecuta con y sin VCF; la
/// alteración de las bases del VCF solo tiene sentido con él.
///
/// Qué guarda caza cada uno, en el orden en que corren:
///
/// - el juez que recuenta, si un documento dice otra cota que el otro o unas
///   bases realizadas que los ficheros no llevan;
/// - la comparación del manifiesto con el TOML, si la cota de la verdad se
///   infla **a la vez** en los dos documentos (solo así se llega a ella);
/// - la comparación de la efectiva con el TOML, si cambia una cifra que solo
///   ella publica: la cota de un locus, el tope declarado, las copias de la
///   verdad por copia de haplotipo, o la cota de la verdad sin tocar su total.
#[test]
fn the_tract_accounting_is_judged_when_publishing_and_bound_to_the_toml() {
    fn manifest_truth_bound(d: &mut Value) {
        let accounting = &mut d["repeat_loci"]["tract_output"];
        shift(&mut accounting["truth_tract_bases_bound"], 1000);
        shift(&mut accounting["total_tract_bases_bound"], 1000);
    }
    fn effective_truth_bound(d: &mut Value) {
        shift(
            &mut d["repeat_output_budget"]["truth_tract_bases_bound"],
            1000,
        );
    }
    fn effective_truth_bound_and_total(d: &mut Value) {
        effective_truth_bound(d);
        shift(
            &mut d["repeat_output_budget"]["total_tract_bases_bound"],
            1000,
        );
    }
    let common: [TractCase; 8] = [
        (
            "manifiesto: cota de la verdad +1000, con su total",
            BOUND_MISMATCH,
            manifest_truth_bound,
            untouched,
        ),
        (
            "manifiesto: bases realizadas de la verdad -2, con su total",
            REALIZED_TRACT,
            |d| {
                let accounting = &mut d["repeat_loci"]["tract_output"];
                shift(&mut accounting["realized_truth_tract_bases"], -2);
                shift(&mut accounting["realized_total_tract_bases"], -2);
            },
            untouched,
        ),
        (
            "efectiva: cota de la verdad +1000",
            EFFECTIVE_TRACT,
            untouched,
            effective_truth_bound,
        ),
        (
            "efectiva: cota de la verdad +1000, con su total",
            BOUND_MISMATCH,
            untouched,
            effective_truth_bound_and_total,
        ),
        (
            "efectiva: cota de un locus +5",
            EFFECTIVE_TRACT,
            untouched,
            |d| {
                shift(
                    &mut d["repeat_output_budget"]["per_locus"][0]["truth_tract_bases_bound"],
                    5,
                )
            },
        ),
        (
            "efectiva: tope declarado 7",
            EFFECTIVE_TRACT,
            untouched,
            |d| d["repeat_output_budget"]["declared_max_tract_output_bases"] = json!(7),
        ),
        (
            "efectiva: copias de la verdad por copia de haplotipo 4 (publica 3)",
            EFFECTIVE_TRACT,
            untouched,
            |d| d["repeat_output_budget"]["truth_copies_per_haplotype_copy"] = json!(4),
        ),
        (
            "los dos: cota de la verdad +1000, con su total",
            MANIFEST_TRACT,
            manifest_truth_bound,
            effective_truth_bound_and_total,
        ),
    ];
    let with_vcf: [TractCase; 1] = [(
        "manifiesto: bases realizadas del VCF -1, con su total",
        REALIZED_TRACT,
        |d| {
            let accounting = &mut d["repeat_loci"]["tract_output"];
            shift(&mut accounting["realized_vcf_tract_bases"], -1);
            shift(&mut accounting["realized_total_tract_bases"], -1);
        },
        untouched,
    )];
    for (package, vcf) in [(publish_mixed(5), false), (publish_mixed_with_vcf(5), true)] {
        let effective = package.json("effective-config.json");
        let manifest = package.json("manifest.json");
        judge_both(&package, &effective, &manifest)
            .expect("control positivo: el paquete intacto se acepta");
        let accounting = &manifest["repeat_loci"]["tract_output"];
        assert_eq!(
            accounting["realized_vcf_tract_bases"]
                .as_u64()
                .expect("VCF")
                > 0,
            vcf,
            "el fixture: el VCF lleva bases de tracto si y solo si se publica"
        );
        assert!(
            accounting["realized_truth_tract_bases"].as_u64()
                < accounting["truth_tract_bases_bound"].as_u64(),
            "el fixture: lo realizado queda por debajo de la cota, así que rebajarlo cuadra"
        );
        let cases = common
            .iter()
            .chain(vcf.then_some(&with_vcf).into_iter().flatten());
        for (label, expected, mutate_manifest, mutate_effective) in cases {
            let mut forged_manifest = manifest.clone();
            mutate_manifest(&mut forged_manifest);
            let mut forged_effective = effective.clone();
            mutate_effective(&mut forged_effective);
            assert!(
                forged_manifest != manifest || forged_effective != effective,
                "«{label}»: la alteración no alteró nada"
            );
            // La aritmética del manifiesto sigue cuadrando: el lector lo acepta.
            judge_manifest(&package, &forged_manifest).unwrap_or_else(|error| {
                panic!("«{label}»: el lector del manifiesto no debía rechazarlo: {error:?}")
            });
            assert_eq!(
                detail(judge_both(&package, &forged_effective, &forged_manifest)),
                *expected,
                "«{label}» (VCF: {vcf}) se rechazó por otra razón"
            );
        }
    }
}

/// **Con VCF, los loci se atan al TOML antes de recontar el tracto** (DT-052).
///
/// El juez que recuenta las bases atribuye cada registro de repetición del VCF
/// a un locus **del manifiesto**. Si corriera antes de la comparación de los
/// loci con el TOML (DT-051), una sección de loci alterada se rechazaría por un
/// síntoma del recuento —«un registro de repetición del VCF no contiene el
/// intervalo de ningún locus», o las bases de un locus que no cuadran— en vez
/// de por su causa. Las pruebas de DT-051 son sin VCF y no veían ese orden; la
/// revisión adversarial de DT-052 comprobó que adelantar el recuento
/// sobrevivía a la batería entera. Estos dos casos, con VCF, exigen el motivo
/// de la comparación con el TOML.
#[test]
fn with_a_vcf_the_loci_are_bound_to_the_toml_before_the_tract_recount() {
    let package = publish_mixed_with_vcf(5);
    assert_manifest_cases_against_the_toml(
        &package,
        &[
            case("secuencia de un locus cambiada", MANIFEST_SECTION, |d| {
                declared(d)["sequence_id"] = json!("chr9")
            }),
            case(
                "identificadores de los dos loci intercambiados",
                MANIFEST_SECTION,
                |d| {
                    let left = declared(d)["id"].clone();
                    let right = generated(d)["id"].clone();
                    declared(d)["id"] = right;
                    generated(d)["id"] = left;
                },
            ),
        ],
    );
}

/// **Sin loci, la cota publicada es cero entera** (DT-052).
///
/// Un documento v4 publica siempre `repeat_output_budget`, aunque no declare
/// loci, y el manifiesto no publica sección de loci. El juez que recuenta no
/// tiene nada que mirar sin esa sección, de modo que antes nada impedía
/// publicar en la efectiva unas bases de VCF inventadas: hasta el 10/10/2026 se
/// **aceptaba** (comprobado contra el código anterior). Ahora la efectiva tiene
/// que ser la que resuelve el TOML con lo medido, que sin loci es cero.
#[test]
fn without_loci_the_effective_tract_bound_is_all_zero() {
    let package = publish_amplified_without_loci("0.3");
    let effective = package.json("effective-config.json");
    let manifest = package.json("manifest.json");
    assert!(
        manifest.get("repeat_loci").is_none(),
        "el fixture: sin loci"
    );
    assert_eq!(
        effective["repeat_output_budget"]["total_tract_bases_bound"],
        json!(0),
        "el fixture: la cota publicada es cero"
    );
    judge_both(&package, &effective, &manifest).expect("control positivo");
    let forged = altered(&effective, "bases de VCF inventadas", &|d| {
        let budget = &mut d["repeat_output_budget"];
        budget["vcf_tract_bases_bound"] = json!(50);
        budget["total_tract_bases_bound"] = json!(50);
    });
    assert_eq!(
        detail(judge_both(&package, &forged, &manifest)),
        EFFECTIVE_TRACT
    );
}

/// El escenario amplificado con loci y un VCF comprimido.
fn publish_amplified_with_vcf() -> Package {
    let toml = format!("{AMPLIFIED_TOML}\n[vcf_output]\nmode = \"bgzf\"\n");
    publish(&[
        ("general.toml", &toml),
        ("reference.fa", AMPLIFIED_FASTA),
        ("population-snvs.tsv", AMPLIFIED_SNVS),
    ])
}

/// **El lado aceptado de DT-052 y DT-043 con paquetes publicados de verdad.**
///
/// Cada paquete que el programa publica pasa por el juez completo al
/// publicarse, de modo que la batería entera ya es el lado aceptado; esta
/// prueba lo hace explícito para las combinaciones que tocan las dos
/// comparaciones nuevas: loci con y sin VCF (de texto y comprimido), un tracto
/// insertado, amplificación con loci (con y sin multiplicador de los
/// homopolímeros, con VCF) y sin loci. Que haya bases de tracto en el VCF
/// donde se espera se comprueba, para que el lado aceptado no sea trivial.
#[test]
fn real_packages_pass_the_tract_and_amplification_comparisons() {
    let packages = [
        ("loci sin VCF", publish_mixed(5), false),
        ("loci con VCF de texto", publish_mixed_with_vcf(5), true),
        (
            "tracto insertado con VCF comprimido",
            publish_inserting(),
            true,
        ),
        ("amplificado con loci", publish_amplified(""), false),
        (
            "amplificado con multiplicador",
            publish_amplified("homopolymer_pcr_fold = 2.5\n"),
            false,
        ),
        (
            "amplificado con loci y VCF",
            publish_amplified_with_vcf(),
            true,
        ),
        (
            "amplificado sin loci",
            publish_amplified_without_loci("0.3"),
            false,
        ),
    ];
    for (label, package, vcf_bases) in &packages {
        let manifest = package.json("manifest.json");
        let realized_vcf = manifest
            .pointer("/repeat_loci/tract_output/realized_vcf_tract_bases")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        assert_eq!(realized_vcf > 0, *vcf_bases, "«{label}»: el fixture");
        judge_both(package, &package.json("effective-config.json"), &manifest)
            .unwrap_or_else(|error| panic!("«{label}» se tenía que aceptar: {error:?}"));
    }
}

// ---------------------------------------------------------------------------
// DT-043: la amplificación del manifiesto de evento, atada a la efectiva
// ---------------------------------------------------------------------------

const EVENT_MANIFEST: &str = "events/reads-main/manifest.json";

/// Reescribe en disco el manifiesto de evento y **rehace sus huellas** en el
/// manifiesto general, que es lo que haría quien quisiera colarlo. Devuelve el
/// manifiesto general resultante.
fn forge_event_manifest(package: &Package, event: &Value) -> Value {
    use sha2::{Digest, Sha256};
    fn rehash(value: &mut Value, digest: &str, size: usize) {
        match value {
            Value::Object(map) => {
                if map.get("path").and_then(Value::as_str) == Some(EVENT_MANIFEST) {
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
    fs::write(package.destination.join(EVENT_MANIFEST), &bytes).expect("manifiesto de evento");
    let mut manifest = package.json("manifest.json");
    rehash(
        &mut manifest,
        &format!("{:x}", Sha256::digest(&bytes)),
        bytes.len(),
    );
    manifest
}

fn amplification_mismatch(key: &str) -> String {
    format!(
        "el campo {key} de la amplificación del manifiesto técnico no es el de la configuración efectiva"
    )
}

/// **Cada parámetro de la amplificación del manifiesto de evento tiene que ser
/// el de la configuración efectiva** (DT-043).
///
/// La efectiva ya está atada exactamente al TOML copiado; el manifiesto de
/// evento publica el mismo bloque más los recuentos realizados. Cada caso
/// cambia **una** clave del bloque en el manifiesto de evento, rehace sus
/// huellas en el manifiesto general y exige el motivo exacto, que nombra la
/// clave. Hasta el 10/10/2026 el juez completo **aceptaba** todos (comprobado
/// contra el código anterior): solo se comparaba el sub-bloque del tartamudeo
/// clonal. El primero es el caso reproducido por la revisión de DT-009:
/// `cycles` a 12 con el tartamudeo clonal en 10.
///
/// Se ejecuta sobre los tres paquetes amplificados: con loci (con tartamudeo
/// clonal, catálogo v6), con el multiplicador de los homopolímeros declarado y
/// sin loci (sin tartamudeo clonal, catálogo v5).
///
/// **Cinco claves del bloque no las alcanza esta comparación**, y no se
/// finge: `applicable`, `algorithm`, `rng_catalog` y el método de
/// `duplicate_fraction_inversion` los fija ya el juez propio del manifiesto de
/// evento como literales (y `applicable` decide la forma entera del bloque), y
/// el sub-bloque `stutter` lo compara antes la igualdad de REQ-VAL-61, con su
/// motivo. Sus casos comprueban que se rechazan por **otra** guarda que corre
/// antes; no discriminan la comparación nueva.
#[test]
fn the_amplification_of_the_event_manifest_is_the_effective_one() {
    type Mutation = fn(&mut Value);
    let reachable: [(&str, Mutation); 11] = [
        ("cycles", |pcr| pcr["cycles"] = json!(12)),
        ("efficiency", |pcr| pcr["efficiency"] = json!(0.8)),
        ("early_cycles", |pcr| shift(&mut pcr["early_cycles"], 1)),
        ("declared_duplicate_fraction", |pcr| {
            pcr["declared_duplicate_fraction"] = json!(0.25);
        }),
        ("expected_duplicate_fraction", |pcr| {
            pcr["expected_duplicate_fraction"] = json!(0.25);
        }),
        ("mean_pairs_per_molecule", |pcr| {
            pcr["mean_pairs_per_molecule"] = json!(1.5);
        }),
        ("duplicate_fraction_inversion", |pcr| {
            shift(&mut pcr["duplicate_fraction_inversion"]["iterations"], 1);
        }),
        ("max_family_size", |pcr| {
            shift(&mut pcr["max_family_size"], 1)
        }),
        ("gc_max_attempts", |pcr| {
            shift(&mut pcr["gc_max_attempts"], 1)
        }),
        ("gc_relative_efficiency", |pcr| {
            pcr["gc_relative_efficiency"][0]["relative_efficiency"] = json!(0.123);
        }),
        ("efficiency_classes", |pcr| {
            pcr["efficiency_classes"][0]["efficiency"] = json!(0.123);
        }),
    ];
    let shadowed: [(&str, Mutation); 5] = [
        ("applicable", |pcr| pcr["applicable"] = json!(false)),
        ("algorithm", |pcr| pcr["algorithm"] = json!("otro-v1")),
        ("rng_catalog", |pcr| {
            pcr["rng_catalog"] = json!("f4-technical-v3")
        }),
        ("duplicate_fraction_inversion", |pcr| {
            pcr["duplicate_fraction_inversion"]["method"] = json!("otro");
        }),
        ("stutter", |pcr| {
            if let Some(block) = pcr.get_mut("stutter") {
                block["cycles"] = json!(12);
            }
        }),
    ];
    for (label, package, clonal) in [
        ("con loci", publish_amplified(""), true),
        (
            "con multiplicador",
            publish_amplified("homopolymer_pcr_fold = 2.5\n"),
            true,
        ),
        ("sin loci", publish_amplified_without_loci("0.3"), false),
    ] {
        let effective = package.json("effective-config.json");
        let original = package.json(EVENT_MANIFEST);
        let effective_pcr = &effective["technical_event"]["artifacts"]["pcr"];
        let keys = effective_pcr.as_object().expect("bloque").len();
        assert_eq!(
            keys,
            if clonal { 15 } else { 14 },
            "«{label}»: el fixture, el bloque de la efectiva tiene sus claves"
        );
        let untouched = forge_event_manifest(&package, &original);
        judge_both(&package, &effective, &untouched)
            .unwrap_or_else(|error| panic!("«{label}»: control positivo rehecho: {error:?}"));
        for (key, mutate) in reachable {
            let mut event = original.clone();
            mutate(&mut event["artifacts"]["pcr"]);
            assert_ne!(event, original, "«{label}», {key}: la alteración altera");
            let manifest = forge_event_manifest(&package, &event);
            assert_eq!(
                detail(judge_both(&package, &effective, &manifest)),
                amplification_mismatch(key),
                "«{label}», {key}"
            );
        }
        for (key, mutate) in shadowed {
            let mut event = original.clone();
            mutate(&mut event["artifacts"]["pcr"]);
            if event == original {
                // Sin tartamudeo clonal no hay sub-bloque que alterar.
                assert!(!clonal, "«{label}», {key}: la alteración altera");
                continue;
            }
            let manifest = forge_event_manifest(&package, &event);
            let motive = detail(judge_both(&package, &effective, &manifest));
            assert_ne!(
                motive,
                amplification_mismatch(key),
                "«{label}», {key}: lo caza otra guarda que corre antes"
            );
        }
        // Se deja el manifiesto de evento como estaba.
        forge_event_manifest(&package, &original);
    }
}

// ---------------------------------------------------------------------------
// DT-054: la verdad y el VCF son de esta ejecución
// ---------------------------------------------------------------------------

const TRUTH_EXECUTION: &str =
    "la cabecera de la verdad publicada no lleva el execution_id del manifiesto general";
const VCF_EXECUTION: &str =
    "la cabecera del VCF publicado no lleva el execution_id del manifiesto general";

/// Rehace, en todo el manifiesto general, la huella y el tamaño del artefacto
/// con esa ruta.
fn rehash_artifact(manifest: &mut Value, path: &str, bytes: &[u8]) {
    use sha2::{Digest, Sha256};
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

/// Cambia el TOML copiado de forma **coherente**: le añade un tope de bases de
/// tracto, rehace su huella en el manifiesto, pone el mismo tope en la
/// efectiva y sustituye el `execution_id` del manifiesto por el que sale del
/// TOML nuevo y la semilla. Devuelve la efectiva, el manifiesto y los dos
/// identificadores de ejecución, el viejo y el nuevo.
fn forge_configuration(package: &Package) -> (Value, Value, String, String) {
    let manifest = package.json("manifest.json");
    let mut effective = package.json("effective-config.json");
    let configuration = catalog_path(
        &manifest,
        manifest["configuration_input"]["artifact_id"]
            .as_str()
            .expect("TOML copiado"),
    );
    let text = String::from_utf8(package.bytes(&configuration)).expect("TOML");
    let anchor = "fasta = \"reference.fa\"\n";
    assert_eq!(text.matches(anchor).count(), 1, "el fixture");
    let altered = text.replacen(
        anchor,
        &format!("{anchor}max_tract_output_bases = 999999\n"),
        1,
    );
    fs::write(package.destination.join(&configuration), &altered).expect("TOML alterado");
    let mut forged = manifest.clone();
    rehash_artifact(&mut forged, &configuration, altered.as_bytes());
    effective["repeat_output_budget"]["declared_max_tract_output_bases"] = json!(999_999);
    let root_seed = genome_population_simulator::randomness::RootSeed::parse(
        "semilla",
        effective["seed"]["root"].as_str().expect("semilla"),
    )
    .expect("semilla");
    let digest = genome_population_simulator::integrity::configuration_digest(altered.as_bytes())
        .expect("huella del TOML");
    let new = genome_population_simulator::identity::derive_execution_id(root_seed, &digest)
        .expect("execution_id")
        .as_str()
        .to_owned();
    let old = manifest["execution_id"]
        .as_str()
        .expect("execution_id")
        .to_owned();
    assert_ne!(old, new, "el fixture: el TOML nuevo da otro execution_id");
    forged["execution_id"] = json!(new);
    (effective, forged, old, new)
}

/// **El `execution_id` de la cabecera de la verdad y del VCF tiene que ser el
/// del manifiesto general** (DT-054).
///
/// Es la reproducción de la revisión adversarial de DT-052: un TOML copiado
/// cambiado de forma coherente (tope de tracto añadido, su huella, el tope de
/// la efectiva y el `execution_id` del manifiesto rehechos) lo **aceptaba** el
/// juez completo hasta el 11/10/2026, aunque la cabecera de la verdad seguía
/// con el identificador de la ejecución original (comprobado contra el código
/// anterior). Sin VCF lo rechaza la comparación de la verdad. Con VCF, si se
/// rehace también la cabecera de la verdad (con su huella), lo rechaza la
/// comparación del VCF; con la verdad sin rehacer, la de la verdad, que va
/// antes. El manifiesto de evento no publica `execution_id`.
#[test]
fn the_truth_and_the_vcf_carry_the_execution_id_of_the_manifest() {
    let package = publish_mixed(5);
    let (effective, forged, _, _) = forge_configuration(&package);
    assert_eq!(
        detail(judge_both(&package, &effective, &forged)),
        TRUTH_EXECUTION,
        "sin VCF"
    );

    let package = publish_mixed_with_vcf(5);
    judge_both(
        &package,
        &package.json("effective-config.json"),
        &package.json("manifest.json"),
    )
    .expect("control positivo");
    let (effective, mut forged, old, new) = forge_configuration(&package);
    assert_eq!(
        detail(judge_both(&package, &effective, &forged)),
        TRUTH_EXECUTION,
        "con VCF y la verdad sin rehacer"
    );
    let truth = catalog_path(&forged, "truth-jsonl");
    let text = String::from_utf8(package.bytes(&truth)).expect("verdad");
    let (header, rest) = text.split_once('\n').expect("cabecera de la verdad");
    let field = |id: &str| format!("\"execution_id\":\"{id}\"");
    assert_eq!(header.matches(&field(&old)).count(), 1, "el fixture");
    let rewritten = format!("{}\n{rest}", header.replace(&field(&old), &field(&new)));
    fs::write(package.destination.join(&truth), &rewritten).expect("verdad rehecha");
    rehash_artifact(&mut forged, &truth, rewritten.as_bytes());
    assert_eq!(
        detail(judge_both(&package, &effective, &forged)),
        VCF_EXECUTION,
        "con VCF y la verdad rehecha"
    );
}
