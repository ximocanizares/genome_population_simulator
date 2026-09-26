//! F4.9.3 y F4.9.11: los rechazos de la preparación amplificada.
//!
//! # Qué se rechaza y por qué importa el control positivo
//!
//! Cada rechazo de este fichero va acompañado de un **control positivo**: el
//! mismo escenario sin el defecto tiene que admitirse. Sin él, una guarda que
//! rechazara todo satisfaría la mitad de las pruebas sin discriminar nada.
//!
//! # La combinación con microsatélites, antes y ahora
//!
//! Hasta la entrega 4, una biblioteca amplificada sobre un escenario con loci
//! declarados se rechazaba **en bloque**, porque el tartamudeo clonal
//! heredable por linaje no existía. Ya existe, y ese rechazo temporal está
//! levantado.
//!
//! Lo que queda en su sitio es una condición permanente y más estrecha: el
//! tartamudeo clonal solo despeja el **incremento** de la amplificación, de
//! modo que sin el tartamudeo de base la biblioteca quedaría por debajo de lo
//! medido —entre un 14 % y un 28 % en motivos de uno a cinco pares de base—.
//! La configuración exige entonces que el tartamudeo de base esté activo.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-9-3-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

struct Attempt {
    root: PathBuf,
    result: Result<(), Diagnostic>,
}

impl Drop for Attempt {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

impl Attempt {
    fn rejection(&self) -> &Diagnostic {
        self.result
            .as_ref()
            .expect_err("la configuración tenía que rechazarse")
    }

    fn accepted(&self) -> bool {
        self.result.is_ok()
    }
}

/// Valida una configuración con el bloque del evento técnico indicado.
///
/// `locus` se pega junto a la referencia, y por omisión está vacío: **ningún
/// paso de este hito cita un locus de repetición** salvo el rechazo temporal.
fn validate(event_extra: &str, locus: &str) -> Attempt {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let tract = "CA".repeat(15);
    let bases = format!("{}{tract}{}", "G".repeat(400), "G".repeat(400));
    let length = bases.len();
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let locus_block = if locus.is_empty() {
        String::new()
    } else {
        format!(
            r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 401
end = {}
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{tract}"
frequency = 1.0
"#,
            400 + tract.len()
        )
    };
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
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
{locus_block}
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
{event_extra}
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        ),
    )
    .expect("configuración");
    let result = load_and_validate_config(&config).map(|_| ());
    Attempt { root, result }
}

/// El bloque amplificado admisible, del que parten todas las variantes.
fn amplified(extra: &str) -> String {
    format!(
        r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.05
cycles = 10
efficiency = 0.9
{extra}"#
    )
}

/// Comprueba un rechazo: su código, su **campo exacto** y su explicación.
///
/// # Por qué el campo se compara por igualdad y no por subcadena
///
/// Un campo es **prefijo de todos sus subcampos**: `technical_event.pcr` está
/// contenido en `technical_event.pcr.cycles`, en `…efficiency` y en todos los
/// demás. Una comprobación por subcadena admitiría un rechazo que nombrara el
/// campo equivocado con tal de que empezara igual. Es exactamente la forma que
/// tenía el defecto de «peso 1» dentro de «peso 100», que ya mordió a este
/// proyecto esta semana.
fn assert_rejects(attempt: &Attempt, field: &str, needle: &str) {
    let rejection = attempt.rejection();
    let rendered = format!("{rejection:?}");
    assert_eq!(rejection.code(), "GPS-CONF-008", "{rendered}");
    assert_eq!(
        rejection.context().get("fields").map(String::as_str),
        Some(field),
        "el rechazo tiene que nombrar exactamente {field}: {rendered}"
    );
    let detail = rejection
        .context()
        .get("detail")
        .map(String::as_str)
        .unwrap_or_default();
    assert!(
        detail.contains(needle),
        "el rechazo tiene que explicar {needle}: {rendered}"
    );
}

/// **Control positivo de toda la batería.** El escenario base se admite.
///
/// Sin esta prueba, una guarda que rechazara cualquier preparación amplificada
/// haría pasar todas las demás de este fichero.
#[test]
fn the_baseline_amplified_scenario_is_admitted() {
    let attempt = validate(&amplified(""), "");
    assert!(
        attempt.accepted(),
        "el escenario amplificado base tiene que admitirse: {:?}",
        attempt.result
    );
}

/// **La tabla de amplificación con preparación sin PCR se rechaza.**
#[test]
fn an_amplification_table_without_the_amplified_preparation_is_rejected() {
    let attempt = validate(
        r#"
[technical_event.pcr]
duplicate_fraction = 0.05
cycles = 10
efficiency = 0.9
"#,
        "",
    );
    assert_rejects(&attempt, "technical_event.pcr", "solo se admite con");
}

/// **La preparación amplificada sin su tabla se rechaza.**
#[test]
fn the_amplified_preparation_without_its_table_is_rejected() {
    let attempt = validate("library_preparation = \"pcr\"\n", "");
    assert_rejects(&attempt, "technical_event.pcr", "obligatorios");
}

/// **Sin fracción de duplicados se rechaza, y sin valor por omisión.**
#[test]
fn an_amplification_without_the_duplicate_fraction_is_rejected() {
    let attempt = validate(
        r#"library_preparation = "pcr"

[technical_event.pcr]
cycles = 10
efficiency = 0.9
"#,
        "",
    );
    assert!(
        attempt.result.is_err(),
        "la fracción de duplicados es obligatoria y no tiene valor por omisión"
    );
}

/// **Sin ciclos se rechaza.**
#[test]
fn an_amplification_without_cycles_is_rejected() {
    let attempt = validate(
        r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.05
efficiency = 0.9
"#,
        "",
    );
    assert!(attempt.result.is_err(), "los ciclos son obligatorios");
}

/// **Treinta y un ciclos se rechazan nombrando el máximo.**
#[test]
fn thirty_one_cycles_are_rejected_naming_the_maximum() {
    let attempt = validate(
        r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.05
cycles = 31
efficiency = 0.9
"#,
        "",
    );
    assert_rejects(&attempt, "technical_event.pcr.cycles", "entre 1 y 30");
    // Control positivo: treinta ciclos **sí** se admiten, de modo que el
    // rechazo viene del tope y no de que los ciclos molesten.
    let accepted = validate(
        r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.05
cycles = 30
efficiency = 0.9
"#,
        "",
    );
    assert!(
        accepted.accepted(),
        "treinta ciclos tienen que admitirse: {:?}",
        accepted.result
    );
}

/// **Eficiencia cero se rechaza.**
#[test]
fn a_zero_efficiency_is_rejected() {
    let attempt = validate(
        r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.05
cycles = 10
efficiency = 0
"#,
        "",
    );
    assert_rejects(&attempt, "technical_event.pcr.efficiency", "(0, 1]");
}

/// **Una fracción de duplicados negativa o igual o mayor que 1 se rechaza.**
#[test]
fn a_duplicate_fraction_outside_its_domain_is_rejected() {
    for value in ["-0.1", "1", "1.5"] {
        let attempt = validate(
            &format!(
                r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = {value}
cycles = 10
efficiency = 0.9
"#
            ),
            "",
        );
        assert!(
            attempt.result.is_err(),
            "la fracción {value} está fuera de [0, 1) y tenía que rechazarse"
        );
    }
    // Control positivo de los dos extremos admisibles.
    for value in ["0", "0.9"] {
        let attempt = validate(
            &format!(
                r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = {value}
cycles = 10
efficiency = 0.9
"#
            ),
            "",
        );
        assert!(
            attempt.accepted(),
            "la fracción {value} es legítima: {:?}",
            attempt.result
        );
    }
}

/// **Una fracción inalcanzable se rechaza nombrando el máximo alcanzable, en
/// vez de pinzarse en silencio.**
#[test]
fn an_unreachable_duplicate_fraction_is_rejected_instead_of_being_clamped() {
    let attempt = validate(
        r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.99
cycles = 10
efficiency = 0.9
"#,
        "",
    );
    assert_rejects(
        &attempt,
        "technical_event.pcr.duplicate_fraction",
        "máxima alcanzable",
    );
}

/// **Una curva con eficiencia relativa cero se rechaza.**
#[test]
fn a_composition_curve_with_zero_efficiency_is_rejected() {
    let attempt = validate(
        &amplified(
            r#"
[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.0
relative_efficiency = 0.0

[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 1.0
relative_efficiency = 1.0
"#,
        ),
        "",
    );
    assert_rejects(
        &attempt,
        "technical_event.pcr.gc_relative_efficiency",
        "fuera de dominio",
    );
}

/// **Una curva con nudos decrecientes se rechaza.**
#[test]
fn a_composition_curve_with_decreasing_knots_is_rejected() {
    let attempt = validate(
        &amplified(
            r#"
[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.7
relative_efficiency = 1.0

[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.3
relative_efficiency = 0.5
"#,
        ),
        "",
    );
    assert_rejects(
        &attempt,
        "technical_event.pcr.gc_relative_efficiency",
        "estrictamente crecientes",
    );
    // Control positivo: los mismos dos nudos en orden creciente se admiten, de
    // modo que el rechazo viene del orden y no de la curva en sí.
    let accepted = validate(
        &amplified(
            r#"
[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.3
relative_efficiency = 0.5

[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.7
relative_efficiency = 1.0
"#,
        ),
        "",
    );
    assert!(
        accepted.accepted(),
        "una curva creciente es legítima: {:?}",
        accepted.result
    );
}

/// **Unas clases de eficiencia cuyos pesos no suman 1 se rechazan.**
#[test]
fn efficiency_classes_whose_weights_do_not_sum_to_one_are_rejected() {
    let attempt = validate(
        &amplified(
            r#"
[[technical_event.pcr.efficiency_classes]]
efficiency = 0.95
weight = 0.7

[[technical_event.pcr.efficiency_classes]]
efficiency = 0.6
weight = 0.5
"#,
        ),
        "",
    );
    assert_rejects(
        &attempt,
        "technical_event.pcr.efficiency_classes",
        "dentro de la tolerancia de normalización",
    );
    // Control positivo: las mismas dos clases con pesos que sí suman 1.
    let accepted = validate(
        &amplified(
            r#"
[[technical_event.pcr.efficiency_classes]]
efficiency = 0.95
weight = 0.7

[[technical_event.pcr.efficiency_classes]]
efficiency = 0.6
weight = 0.3
"#,
        ),
        "",
    );
    assert!(
        accepted.accepted(),
        "unas clases normalizadas son legítimas: {:?}",
        accepted.result
    );
}

/// **Una clave desconocida dentro de la tabla se rechaza.**
///
/// Es lo que falla contra `additionalProperties: true` por el lado del lector
/// de Rust: el esquema lo impide por su lado y esta prueba por el otro.
#[test]
fn an_unknown_key_inside_the_amplification_table_is_rejected() {
    let attempt = validate(&amplified("mean_pairs_per_molecule = 0.1\n"), "");
    let rejection = attempt.rejection();
    assert_eq!(
        rejection.code(),
        "GPS-CONF-006",
        "una clave desconocida tiene su propio código: {rejection:?}"
    );
}

/// **Una preparación inventada se rechaza nombrando las admitidas.**
#[test]
fn an_invented_library_preparation_is_rejected() {
    let attempt = validate("library_preparation = \"nanopore\"\n", "");
    assert_rejects(&attempt, "technical_event.library_preparation", "no existe");
}

/// **F4.10.6. El rechazo temporal está levantado: la combinación se acepta.**
///
/// Hasta la entrega 4 esta combinación se rechazaba en bloque, porque el
/// tartamudeo clonal heredable por linaje no existía. Ya existe, y lo que
/// queda es la condición permanente que comprueba la prueba siguiente.
///
/// **Debe fallar contra** un levantamiento a medias que deje el rechazo en pie.
///
/// # Lo que esta prueba **no** comprueba, dicho sin adornos
///
/// Solo mira que la configuración se **acepte**. No comprueba que el
/// levantamiento venga acompañado del bloque de tartamudeo clonal: una
/// implementación que aceptara la combinación y no abriera ninguna decisión de
/// linaje pasaría por aquí sin inmutarse. Quien caza eso son otras tres, y se
/// nombran para que nadie crea que esta las cubre:
/// `the_amplified_run_reproduces_the_measured_with_pcr_fractions` y
/// `the_amplified_manifest_publishes_the_derivation_and_the_lineage_counts`,
/// las dos en `tests/clonal_composition.rs`, y
/// `the_two_parts_add_up_on_every_published_read` en `tests/clonal_lineage.rs`.
#[test]
fn the_amplified_preparation_with_declared_repeat_loci_is_now_accepted() {
    let attempt = validate(&amplified(&stutter_block()), "STR-1");
    assert!(
        attempt.accepted(),
        "la combinación ya no se rechaza: {:?}",
        attempt.result
    );
}

/// **Lo que queda en pie: con loci declarados, la amplificación exige el
/// tartamudeo de base.**
///
/// El tartamudeo clonal solo despeja el **incremento** de la amplificación.
/// Sin la parte de base la aritmética no cierra y la biblioteca quedaría entre
/// un 14 % y un 28 % por debajo de lo medido. El rechazo lo dice con esas
/// cifras.
///
/// **Debe fallar contra** un levantamiento que acepte la combinación sin
/// comprobar que el tartamudeo de base está activo.
#[test]
fn amplification_over_declared_loci_demands_the_base_stutter() {
    let attempt = validate(&amplified(""), "STR-1");
    let rejection = attempt.rejection();
    let rendered = format!("{rejection:?}");
    assert_eq!(rejection.code(), "GPS-CONF-008", "{rendered}");
    assert!(
        rendered.contains("technical_event.library_preparation,technical_event.stutter.enabled"),
        "nombra los dos campos que se cruzan: {rendered}"
    );
    assert!(
        rendered.contains("STR-1"),
        "nombra el locus declarado: {rendered}"
    );
    assert!(
        rendered.contains("14 %") && rendered.contains("28 %"),
        "explica cuánto se subestimaría: {rendered}"
    );
    // Control positivo del **mismo** escenario: con el interruptor apagado
    // explícitamente también se rechaza, de modo que la guarda mira el valor y
    // no la mera presencia de la tabla.
    let switched_off = validate(
        &amplified("\n[technical_event.stutter]\nenabled = false\n"),
        "STR-1",
    );
    assert!(
        switched_off.result.is_err(),
        "una tabla presente con el interruptor apagado no basta: {:?}",
        switched_off.result
    );
}

/// **La sustitución directa de la probabilidad por ciclo se valida.**
///
/// Las dos comprobaciones miran cosas distintas: una longitud de motivo fuera
/// del tramo calibrado y una probabilidad fuera de dominio. Cada una lleva su
/// control positivo con el mismo escenario, que es lo que impide que una
/// guarda que rechazara toda la tabla las satisficiera a las dos.
#[test]
fn a_declared_per_cycle_slip_is_checked_against_its_domain() {
    let out_of_range = validate(
        &amplified(&format!(
            "{}\n[[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 7\nprobability = 0.01\n",
            stutter_block()
        )),
        "STR-1",
    );
    assert_rejects(
        &out_of_range,
        "technical_event.pcr.stutter.per_cycle_slip.motif_length_bases",
        "fuera del tramo calibrado",
    );
    let out_of_domain = validate(
        &amplified(&format!(
            "{}\n[[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 2\nprobability = 1.5\n",
            stutter_block()
        )),
        "STR-1",
    );
    assert_rejects(
        &out_of_domain,
        "technical_event.pcr.stutter.per_cycle_slip.probability",
        "es la probabilidad de UN ciclo",
    );
    let duplicated = validate(
        &amplified(&format!(
            "{}\n[[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 2\nprobability = 0.01\n\
             \n[[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 2\nprobability = 0.02\n",
            stutter_block()
        )),
        "STR-1",
    );
    assert_rejects(
        &duplicated,
        "technical_event.pcr.stutter.per_cycle_slip.motif_length_bases",
        "se declara dos veces",
    );
    // Control positivo: una sustitución legítima se admite.
    let accepted = validate(
        &amplified(&format!(
            "{}\n[[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 2\nprobability = 0.02\n",
            stutter_block()
        )),
        "STR-1",
    );
    assert!(
        accepted.accepted(),
        "una sustitución dentro de dominio es legítima: {:?}",
        accepted.result
    );
}

/// Bloque de tartamudeo de base mínimo y admisible: interruptor y pesos, que
/// son obligatorios cuando está activo.
fn stutter_block() -> String {
    r#"
[technical_event.stutter]
enabled = true

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
"#
    .to_owned()
}

/// **El mismo locus sin preparación amplificada se admite**, y la misma
/// preparación amplificada sin loci también.
///
/// Es el par de controles positivos que demuestra que la guarda mira **las dos
/// configuraciones** y no una sola. Sin ellos, una guarda que rechazara
/// cualquier locus, o cualquier amplificación, pasaría la prueba de arriba.
#[test]
fn either_half_alone_is_admitted() {
    let only_locus = validate("", "STR-1");
    assert!(
        only_locus.accepted(),
        "un locus sin amplificación es legítimo: {:?}",
        only_locus.result
    );
    let only_amplification = validate(&amplified(""), "");
    assert!(
        only_amplification.accepted(),
        "una amplificación sin loci es legítima: {:?}",
        only_amplification.result
    );
}

/// **Los ciclos explícitos no pueden superar a los ciclos.**
#[test]
fn early_cycles_above_the_declared_cycles_are_rejected() {
    let attempt = validate(
        r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.05
cycles = 2
efficiency = 0.9
early_cycles = 5
"#,
        "",
    );
    assert_rejects(
        &attempt,
        "technical_event.pcr.early_cycles",
        "no superar los 2 ciclos declarados",
    );
    // Control positivo: con ciclos explícitos iguales a los ciclos se admite.
    let accepted = validate(
        r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.05
cycles = 2
efficiency = 0.9
early_cycles = 2
"#,
        "",
    );
    assert!(
        accepted.accepted(),
        "dos ciclos explícitos sobre dos ciclos son legítimos: {:?}",
        accepted.result
    );
}

/// Los dos topes opcionales se comprueban contra su constante de contrato.
#[test]
fn the_optional_caps_are_checked_against_their_contract_maximum() {
    let too_big = validate(&amplified("max_family_size = 257\n"), "");
    assert_rejects(
        &too_big,
        "technical_event.pcr.max_family_size",
        "entre 1 y 256",
    );
    let too_many = validate(&amplified("gc_max_attempts = 17\n"), "");
    assert_rejects(
        &too_many,
        "technical_event.pcr.gc_max_attempts",
        "entre 1 y 16",
    );
    // Control positivo de los dos, justo en su tope.
    assert!(validate(&amplified("max_family_size = 256\n"), "").accepted());
    assert!(validate(&amplified("gc_max_attempts = 16\n"), "").accepted());
}

/// **Cinco clases de eficiencia se rechazan; cuatro se admiten.**
#[test]
fn more_efficiency_classes_than_the_contract_allows_are_rejected() {
    let five = (0..5)
        .map(|_| {
            "\n[[technical_event.pcr.efficiency_classes]]\nefficiency = 0.9\nweight = 0.2\n"
                .to_owned()
        })
        .collect::<String>();
    let attempt = validate(&amplified(&five), "");
    assert_rejects(
        &attempt,
        "technical_event.pcr.efficiency_classes",
        "el máximo es 4",
    );
    let four = (0..4)
        .map(|_| {
            "\n[[technical_event.pcr.efficiency_classes]]\nefficiency = 0.9\nweight = 0.25\n"
                .to_owned()
        })
        .collect::<String>();
    let accepted = validate(&amplified(&four), "");
    assert!(
        accepted.accepted(),
        "cuatro clases son legítimas: {:?}",
        accepted.result
    );
}
