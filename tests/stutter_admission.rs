//! F4.8.1.1 y F4.8.3: la guarda de admisión y los rechazos de configuración.
//!
//! # Qué guarda la guarda, en una frase
//!
//! Un locus que **vaya a tartamudear** necesita una tirada mínima efectiva de
//! al menos dos unidades. El motivo es que el paso de unidades no permite que
//! una contracción vacíe su tirada, de modo que un tracto cuya tirada más larga
//! sea de una sola copia **no podría contraerse nunca mientras que expandirse
//! sí**: una clase entera de tractos quedaría congelada frente a la contracción
//! y el modelo saldría sesgado hacia el alargamiento.
//!
//! # Por qué las dos configuraciones, y no solo el locus
//!
//! El microsatélite se declara junto a la referencia y el tartamudeo en el
//! evento técnico: son dos sitios distintos del fichero. Una guarda que mirara
//! **solo el locus** rechazaría la tirada mínima 1 aunque no hubiera
//! tartamudeo, y eso sería un rechazo injustificado. La segunda comprobación de
//! `a_minimum_run_of_one_is_admitted_without_stutter` es la que lo caza.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;

#[path = "support/detection_section.rs"]
mod detection_section;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-8-1-1-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

const FOUR_STEPS: &str = r#"
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
"#;

struct Attempt {
    root: PathBuf,
    result: Result<(), genome_population_simulator::diagnostics::Diagnostic>,
}

impl Drop for Attempt {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

impl Attempt {
    fn rejection(&self) -> &genome_population_simulator::diagnostics::Diagnostic {
        self.result
            .as_ref()
            .expect_err("la configuración tenía que rechazarse")
    }

    fn accepted(&self) -> bool {
        self.result.is_ok()
    }
}

/// Valida una configuración con el locus y el bloque de tartamudeo indicados.
///
/// `locus_extra` se pega dentro de la tabla del locus, que es donde vive su
/// tirada mínima propia.
fn validate(locus_extra: &str, stutter: &str, motif: &str, tract: &str) -> Attempt {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let bases = format!("{}{tract}{}", "G".repeat(100), "G".repeat(100));
    let length = bases.len();
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let start = 101;
    let end = 100 + tract.len();
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
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
length = {length}

[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = {start}
end = {end}
motif = "{motif}"
{locus_extra}
[[reference.repeat_loci.alleles]]
sequence = "{tract}"
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
target_coverage_wgs = 1
{stutter}
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        )),
    )
    .expect("configuración");
    let result = load_and_validate_config(&config).map(|_| ());
    Attempt { root, result }
}

fn enabled_stutter() -> String {
    format!("\n[technical_event.stutter]\nenabled = true\n{FOUR_STEPS}")
}

/// Un tracto de treinta bases que **sí** tiene una tirada larga, para que el
/// rechazo no pueda venir de la admisión ordinaria del locus.
fn pure_tract() -> String {
    "CA".repeat(15)
}

/// **La guarda rechaza la tirada mínima 1 cuando hay tartamudeo, nombrando el
/// locus y los dos valores.**
#[test]
fn a_minimum_run_of_one_is_rejected_when_the_locus_will_stutter() {
    let attempt = validate(
        "min_longest_run_units = 1\n",
        &enabled_stutter(),
        "CA",
        &pure_tract(),
    );
    let rejection = attempt.rejection();
    let rendered = format!("{rejection:?}");
    // **Nombrar el código y el campo, no conformarse con que aparezcan las
    // cifras sueltas.** Casi cualquier rechazo del sistema contiene un `1` y un
    // `2` en alguna parte; eso no distinguía este rechazo de ningún otro.
    assert_eq!(rejection.code(), "GPS-CONF-008");
    assert!(
        rendered.contains("reference.repeat_loci.min_longest_run_units,technical_event.stutter"),
        "nombra los dos campos que se cruzan: {rendered}"
    );
    assert!(rendered.contains("STR-1"), "nombra el locus: {rendered}");
    assert!(
        rendered.contains("tirada mínima efectiva de 1 unidad(es)"),
        "nombra el valor declarado: {rendered}"
    );
    assert!(
        rendered.contains("exige al menos 2"),
        "nombra el mínimo exigido: {rendered}"
    );
}

/// **El mismo locus sin tartamudeo se admite.**
///
/// Es la comprobación que demuestra que la guarda mira **las dos
/// configuraciones** y no solo el locus. Una guarda escrita junto a la lectura
/// del locus rechazaría también este caso, y es exactamente la variante
/// defectuosa que el plan manda cazar.
#[test]
fn a_minimum_run_of_one_is_admitted_without_stutter() {
    let attempt = validate("min_longest_run_units = 1\n", "", "CA", &pure_tract());
    assert!(
        attempt.accepted(),
        "sin tartamudeo la tirada mínima 1 es legítima: {:?}",
        attempt.result
    );
}

/// Con tirada mínima 2 y tartamudeo activo se admite: la guarda pide dos, no
/// tres.
#[test]
fn a_minimum_run_of_two_is_admitted_with_stutter() {
    let attempt = validate(
        "min_longest_run_units = 2\n",
        &enabled_stutter(),
        "CA",
        &pure_tract(),
    );
    assert!(attempt.accepted(), "{:?}", attempt.result);
}

/// **Un locus que no declara mínimo propio hereda el general.**
///
/// Sin esta comprobación, una guarda que mirara el mínimo **declarado** en vez
/// del **efectivo** dejaría pasar un locus sin mínimo propio bajo un general
/// de uno. Aquí el general se baja a uno a propósito y el locus calla: el
/// mínimo efectivo es uno y tiene que rechazarse igual.
#[test]
fn the_guard_looks_at_the_effective_minimum_and_not_at_the_declared_one() {
    let attempt = validate(
        "",
        &format!(
            "{}\n[reference.repeat_admission]\nmin_longest_run_units = 1\n",
            enabled_stutter()
        ),
        "CA",
        &pure_tract(),
    );
    // El bloque de admisión general se escribe **después** del evento técnico,
    // que es donde el TOML lo admite: la tabla de la referencia ya está cerrada
    // arriba. Si el orden cambiara, la configuración no sería válida y la
    // prueba dejaría de comprobar lo que dice.
    let rendered = format!("{:?}", attempt.rejection());
    assert!(
        rendered.contains("STR-1") && rendered.contains("tartamudeo"),
        "el mínimo efectivo heredado de 1 tiene que rechazarse por la guarda, no por otra cosa: \
         {rendered}"
    );
}

/// Un locus que no declara mínimo propio y hereda el general de 3 entra sin
/// tocar nada.
#[test]
fn a_locus_without_its_own_minimum_inherits_the_general_three_and_is_admitted() {
    let attempt = validate("", &enabled_stutter(), "CA", &pure_tract());
    assert!(attempt.accepted(), "{:?}", attempt.result);
}

/// **F4.8.3. El tartamudeo activo sin pesos declarados se rechaza nombrando el
/// parámetro.**
#[test]
fn enabled_stutter_without_declared_weights_is_rejected() {
    let attempt = validate(
        "",
        "\n[technical_event.stutter]\nenabled = true\n",
        "CA",
        &pure_tract(),
    );
    let rendered = format!("{:?}", attempt.rejection());
    assert!(
        rendered.contains("technical_event.stutter.steps"),
        "nombra el parámetro: {rendered}"
    );

    // **El rechazo tiene que traer el ejemplo.** Quien tropieza aquí no puede
    // seguir sin cuatro números que nadie ha medido para este uso, y sin el
    // ejemplo no sabe por dónde empezar. Se exigen los cuatro pesos **y** la
    // advertencia, porque dar la tabla sin decir que no es una calibración
    // sería peor que no darla: alguien la copiaría creyéndola medida.
    // **Las parejas enteras, no fragmentos.** «peso 1» y «peso 10» son
    // subcadenas de «peso 100»: comprobarlas sueltas dejaba borrar la mitad
    // del ejemplo sin que la prueba se enterase. Lo destapó la revisión de
    // cierre, y es la forma «aserción implicada por una condición anterior».
    for fragmento in [
        "-1 con peso 100",
        "+1 con peso 20",
        "-2 con peso 10",
        "-3 con peso 1",
        "ILUSTRATIVO",
    ] {
        assert!(
            rendered.contains(fragmento),
            "el rechazo tiene que traer el ejemplo y su advertencia, y le falta {fragmento:?}: {rendered}"
        );
    }
}

/// Un paso `+2` no pertenece a la estructura y se rechaza antes de simular.
#[test]
fn a_step_outside_the_admitted_structure_is_rejected() {
    let steps = FOUR_STEPS.replace("units = 1\nweight = 0.2", "units = 2\nweight = 0.2");
    let attempt = validate(
        "",
        &format!("\n[technical_event.stutter]\nenabled = true\n{steps}"),
        "CA",
        &pure_tract(),
    );
    let rejection = attempt.rejection();
    let rendered = format!("{rejection:?}");
    assert_eq!(rejection.code(), "GPS-CONF-008");
    assert!(
        rendered.contains("technical_event.stutter.steps.units"),
        "nombra el campo: {rendered}"
    );
    assert!(
        rendered.contains("el paso 2 no pertenece a la estructura del algoritmo"),
        "nombra el paso y por qué se rechaza: {rendered}"
    );
}

/// Un peso negativo se rechaza.
#[test]
fn a_negative_step_weight_is_rejected() {
    let steps = FOUR_STEPS.replace("weight = 0.6", "weight = -0.6");
    let attempt = validate(
        "",
        &format!("\n[technical_event.stutter]\nenabled = true\n{steps}"),
        "CA",
        &pure_tract(),
    );
    let rendered = format!("{:?}", attempt.rejection());
    assert!(
        rendered.contains("el peso -0.6 del paso -1 debe ser finito y no negativo"),
        "nombra el peso y el paso: {rendered}"
    );
}

/// Una tasa fuera de `[0,1]` se rechaza.
#[test]
fn a_rate_outside_the_unit_interval_is_rejected() {
    let attempt = validate(
        "",
        &format!(
            "\n[technical_event.stutter]\nenabled = true\n\n[[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = 1.5\n{FOUR_STEPS}"
        ),
        "CA",
        &pure_tract(),
    );
    let rendered = format!("{:?}", attempt.rejection());
    assert!(
        rendered.contains("la tasa 1.5 del motivo de 2 pb"),
        "nombra la tasa y la longitud de motivo: {rendered}"
    );
}

/// **Un paso expresado en bases no existe en el contrato.**
///
/// La tabla del paso no admite claves desconocidas: `bases` se rechaza al leer
/// el TOML, antes de que nadie pueda interpretarla como unidades.
#[test]
fn a_step_expressed_in_base_pairs_is_rejected() {
    let steps = FOUR_STEPS.replace("units = -1", "units = -1\nbases = -2");
    let attempt = validate(
        "",
        &format!("\n[technical_event.stutter]\nenabled = true\n{steps}"),
        "CA",
        &pure_tract(),
    );
    let rejection = attempt.rejection();
    let rendered = format!("{rejection:?}");
    // El rechazo es **de forma**, no semántico: la clave es desconocida para la
    // tabla del paso, y por eso el código es el del campo no admitido y no el
    // de la combinación incompatible. El diagnóstico nombra la clave y su
    // línea, que es lo que el lector de TOML publica para un campo desconocido;
    // no nombra la tabla padre, y exigirlo sería exigir algo que este
    // diagnóstico no da en ningún sitio del proyecto.
    assert_eq!(rejection.code(), "GPS-CONF-006");
    assert!(
        rendered.contains("\"field\": \"bases\""),
        "nombra exactamente la clave rechazada: {rendered}"
    );
    // El número de línea exacto no se fija: cambiaría al tocar la plantilla del
    // fixture y la prueba moriría por un motivo que no tiene nada que ver.
    assert!(
        rendered.contains("\"location\": \"line "),
        "y dice dónde está, para poder corregirla: {rendered}"
    );
}

/// Una longitud de motivo fuera del tramo calibrado no se extrapola: se
/// rechaza nombrando el locus y la longitud.
#[test]
fn a_motif_longer_than_the_calibration_is_rejected_instead_of_extrapolated() {
    // Motivo de siete bases, repetido cinco veces: 35 pb de tracto puro.
    let attempt = validate("", &enabled_stutter(), "ACGTACG", &"ACGTACG".repeat(5));
    let rendered = format!("{:?}", attempt.rejection());
    assert!(rendered.contains("STR-1"), "nombra el locus: {rendered}");
    assert!(
        rendered.contains("declara un motivo de 7 pb y la calibración del tartamudeo solo cubre"),
        "dice qué longitud y por qué no se extrapola: {rendered}"
    );
}
