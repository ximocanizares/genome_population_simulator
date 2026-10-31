//! F6.4.2: guarda de admisión de un locus declarado sobre una **referencia
//! real**.
//!
//! Comprueba que en la posición declarada hay de verdad un tracto de los
//! motivos declarados. Sin ella se anotaría como microsatélite un trozo
//! cualquiera del genoma.
//!
//! # Cómo están elegidos los fixtures, que es lo que hace o no discriminante a
//! este fichero
//!
//! El criterio son **dos** mínimos —pureza y unidades de la tirada más larga— y
//! el peligro evidente es que todos los casos rechazados incumplan los dos a la
//! vez: entonces quitar cualquiera de las dos comprobaciones no mataría nada.
//! Por eso hay una pareja construida al revés:
//!
//! - un tramo que **cumple la pureza y no la tirada**, y
//! - un tramo que **cumple la tirada y no la pureza**.
//!
//! Y la pareja del umbral —justo por debajo y justo por encima— está construida
//! de modo que **lo único que cambia entre las dos es la pureza**: las dos
//! tienen exactamente cien bases y exactamente la misma tirada más larga, cinco,
//! de modo que la tirada no puede ser lo que las separa.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;

#[path = "support/detection_section.rs"]
mod detection_section;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

/// Cien bases con `runs` tiradas de `units` copias de `CA`, separadas por `TT`,
/// y relleno de `T` hasta completar.
///
/// El separador es `TT` a propósito: pegado a un `CA` no forma ninguna copia
/// nueva, de modo que el número de copias es exactamente el que se pide.
fn tract_of_100(runs: usize, units: usize, short_run_units: usize) -> String {
    let mut bases = String::new();
    for index in 0..runs {
        let this = if index == 0 { short_run_units } else { units };
        bases.push_str(&"CA".repeat(this));
        bases.push_str("TT");
    }
    while bases.len() < 100 {
        bases.push('T');
    }
    assert_eq!(bases.len(), 100, "el tramo mide siempre cien bases");
    bases
}

fn scenario(tract: &str, admission: &str) -> (PathBuf, String) {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    let fasta = root.join("reference.fa");
    // El tracto ocupa las bases 101 a 200 de una secuencia de 300.
    let sequence = "G".repeat(100) + tract + &"G".repeat(100);
    assert_eq!(sequence.len(), 300);
    fs::write(&fasta, format!(">chr1\n{sequence}\n")).expect("FASTA");
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[reference]
kind = "supplied"
id = "ref-real"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 300

[[reference.repeat_loci]]
id = "STR-R"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
{admission}

[[reference.repeat_loci.alleles]]
units = 20
frequency = 1.0
"#
        )),
    )
    .expect("configuración");
    (root, config.display().to_string())
}

fn admit(tract: &str, admission: &str) -> Result<(), Diagnostic> {
    let (root, config) = scenario(tract, admission);
    let result = load_and_validate_config(Path::new(&config));
    // La guarda vive en la preparación de la referencia, que es donde las bases
    // están materializadas; se llega a ella ejecutando.
    let outcome = match result {
        Err(diagnostic) => Err(diagnostic),
        Ok(effective) => genome_population_simulator::execution::execute_validated_config(
            &effective,
            &PathBuf::from(&root).join("published"),
        )
        .map(|_| ()),
    };
    fs::remove_dir_all(&root).ok();
    outcome
}

/// Un tramo **puro** se admite.
///
/// *Cuidado:* es el caso de referencia y **no discrimina casi nada**. Un tracto
/// puro pasa cualquier umbral de pureza y cualquier mínimo de tirada razonable;
/// está aquí para que se vea que la guarda no rechaza lo que debe admitir, no
/// como guarda de nada.
#[test]
fn a_pure_tract_in_the_reference_is_admitted() {
    // Cincuenta copias seguidas de `CA`: cien bases, pureza 1, tirada 50.
    admit(&"CA".repeat(50), "").expect("un tramo puro se admite");
}

/// **Un tracto real impuro por encima del mínimo se admite.**
///
/// Es el caso que una guarda que exigiera pureza 1 rechazaría por error, y por
/// eso el fixture es impuro a propósito: cuarenta copias sobre cien bases dan
/// pureza 0,80, por encima del mínimo de 0,75 que declara.
#[test]
fn a_real_impure_tract_above_the_minimum_is_admitted() {
    let tract = tract_of_100(8, 5, 5);
    admit(&tract, "min_purity = 0.75\nmin_longest_run_units = 5")
        .expect("ochenta de cien bases en copias completas");
}

/// **Una región sin ninguna repetición del motivo se rechaza**, nombrando el
/// locus.
///
/// Desde la entrega 3 del plan de detección el rechazo llega **antes** de esta
/// guarda: un locus declarado tiene que coincidir con un microsatélite
/// detectado (REQ-REF-40), y en cien bases de `G` no hay ninguno de su motivo,
/// de modo que el diagnóstico es el de la coincidencia y no nombra la pureza.
/// Que la guarda siga viva lo demuestran las parejas de abajo, que pasan la
/// coincidencia y mueren si la guarda desaparece.
#[test]
fn a_region_without_a_single_repeat_is_rejected() {
    let error = admit(&"G".repeat(100), "").expect_err("cien bases sin ni una copia");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-R"), "{rendered}");
    assert!(
        rendered.contains("no coincide con ningún microsatélite detectado"),
        "{rendered}"
    );
}

/// **Justo por debajo del mínimo se rechaza y justo por encima se admite**, y lo
/// único que cambia entre los dos tramos es la pureza.
///
/// Los dos miden cien bases y los dos tienen la misma tirada más larga, cinco
/// copias: el que se rechaza tiene una tirada de cuatro en vez de cinco en el
/// primer bloque, de modo que pierde una copia —pureza 0,78 frente a 0,80— sin
/// tocar la tirada máxima. **Si la tirada fuera lo que los separa, esta prueba no
/// diría nada sobre la pureza.**
#[test]
fn just_below_the_purity_minimum_is_rejected_and_just_above_is_admitted() {
    let above = tract_of_100(8, 5, 5);
    let below = tract_of_100(8, 5, 4);
    let minimum = "min_purity = 0.80\nmin_longest_run_units = 5";

    admit(&above, minimum).expect("pureza 0,80 con el mínimo en 0,80");
    let error = admit(&below, minimum).expect_err("pureza 0,78 con el mínimo en 0,80");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-R"), "{rendered}");
    assert!(rendered.contains("pureza 0.78"), "{rendered}");

    // **Bajar el mínimo convierte el rechazo en admisión**: el mismo tramo, otra
    // configuración. Es lo que demuestra que el umbral se lee y no está
    // codificado en el programa.
    admit(&below, "min_purity = 0.75\nmin_longest_run_units = 5")
        .expect("el mismo tramo con el mínimo en 0,75");
}

/// **La pareja que impide quitar cualquiera de las dos comprobaciones.**
///
/// El primer tramo **cumple la pureza y no la tirada**; el segundo **cumple la
/// tirada y no la pureza**. Quitar la comprobación de tirada deja pasar el
/// primero; quitar la de pureza deja pasar el segundo. Sin esta pareja, un
/// rechazo que incumpliera las dos cosas a la vez dejaría sobrevivir a las dos
/// mutaciones.
#[test]
fn each_of_the_two_minima_rejects_something_the_other_admits() {
    let minimum = "min_purity = 0.80\nmin_longest_run_units = 5";

    // Diez tiradas de cuatro copias: cuarenta copias, pureza 0,80 —cumple—, y
    // tirada más larga cuatro —no cumple—.
    let short_runs = tract_of_100(10, 4, 4);
    let error = admit(&short_runs, minimum).expect_err("tiradas demasiado cortas");
    let rendered = format!("{error:?}");
    assert!(
        rendered.contains("tirada más larga de 4"),
        "lo que rechaza es la tirada, no la pureza: {rendered}"
    );
    assert!(rendered.contains("el mínimo es 5"), "{rendered}");

    // Y el de la prueba anterior: tirada cinco —cumple— y pureza 0,78 —no—.
    let impure = tract_of_100(8, 5, 4);
    let error = admit(&impure, minimum).expect_err("pureza insuficiente");
    assert!(
        format!("{error:?}").contains("pureza 0.78"),
        "lo que rechaza es la pureza, no la tirada: {error:?}"
    );
}

/// Un locus compuesto se observa contra **todos** sus motivos.
///
/// Con `motifs = ["CAG", "CA"]` sobre un tramo `(CAG)n(CA)m`, una guarda que
/// mirase solo el primer motivo vería la mitad del tramo como interrupción y
/// rechazaría un locus legítimo.
#[test]
fn a_composite_locus_is_observed_against_every_one_of_its_motifs() {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    let tract = "CAG".repeat(20) + &"CA".repeat(20); // 60 + 40 = 100 bases puras
    assert_eq!(tract.len(), 100);
    let fasta = root.join("reference.fa");
    fs::write(
        &fasta,
        format!(">chr1\n{}\n", "G".repeat(100) + &tract + &"G".repeat(100)),
    )
    .expect("FASTA");
    let config = root.join("general.toml");
    let body = |motifs: &str| {
        format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[reference]
kind = "supplied"
id = "ref-real"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 300

[[reference.repeat_loci]]
id = "STR-R"
sequence = "chr1"
start = 101
end = 200
{motifs}

[[reference.repeat_loci.alleles]]
units = {}
frequency = 1.0
"#,
            if motifs.contains("motifs") {
                "[10, 10]"
            } else {
                "20"
            }
        )
    };

    fs::write(
        &config,
        detection_section::with_detection_section(&body("motifs = [\"CAG\", \"CA\"]")),
    )
    .expect("configuración");
    let effective = load_and_validate_config(&config).expect("locus compuesto");
    genome_population_simulator::execution::execute_validated_config(
        &effective,
        &root.join("published-composite"),
    )
    .expect("el tramo es puro contra los dos motivos");

    // Y con un solo motivo el mismo tramo no llega al mínimo, porque las cuarenta
    // bases del segundo bloque dejan de contar.
    fs::write(
        &config,
        detection_section::with_detection_section(&body("motif = \"CAG\"")),
    )
    .expect("configuración");
    let outcome = load_and_validate_config(&config).and_then(|effective| {
        genome_population_simulator::execution::execute_validated_config(
            &effective,
            &root.join("published-single"),
        )
        .map(|_| ())
    });
    let error = outcome.expect_err("leído con un solo motivo no llega al mínimo");
    assert!(format!("{error:?}").contains("pureza 0.6"), "{error:?}");
    fs::remove_dir_all(&root).ok();
}

/// Con **referencia simulada** la guarda no se aplica, y es deliberado: allí el
/// tracto no existe hasta que se crea, así que no hay nada observado contra lo
/// que contrastar.
#[test]
fn the_guard_does_not_apply_to_a_simulated_reference() {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[reference]
kind = "simulated"
id = "ref-main"

[reference.composition]
A = 1.0
C = 0.0
G = 0.0
T = 0.0

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 300

[[reference.repeat_loci]]
id = "STR-R"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 20
frequency = 1.0
"#
        )),
    )
    .expect("configuración");
    // La referencia generada son trescientas `A`: no contiene ni una copia de
    // `CA`. Si la guarda se aplicara aquí, esto se rechazaría.
    let effective = load_and_validate_config(&config).expect("referencia simulada");
    genome_population_simulator::execution::execute_validated_config(
        &effective,
        &root.join("published"),
    )
    .expect("con referencia simulada el tracto se escribe, no se lee");
    fs::remove_dir_all(&root).ok();
}

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f6-4-2-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}
