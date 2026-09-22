//! Decisión S-32: un cambio pequeño **dentro** del intervalo de un locus de
//! repetición rechaza el escenario, y lo rechaza **antes de simular**.
//!
//! # Cómo están elegidos los fixtures, que es donde esta prueba se juega
//!
//! El defecto que estas pruebas cubren tenía una trampa: la guarda del estado
//! solo se consultaba **cuando había un alelo que aplicar**, y el alelo se
//! sortea. Con el alelo de la referencia entre los declarados, la misma
//! configuración se aceptaba o se rechazaba **según la semilla**.
//!
//! De ahí las dos condiciones que todo fixture de este fichero cumple, y sin
//! las cuales no discriminaría nada:
//!
//! 1. **El alelo de la referencia está entre los declarados.** Sin él, toda
//!    copia recibe siempre un alelo que hay que aplicar, se consulta siempre la
//!    guarda del estado, y el código defectuoso pasaba la prueba igual.
//! 2. **El rechazo se exige con varias semillas.** Una sola no distingue
//!    «rechaza siempre» de «rechazó con esta».
//!
//! Y cada rechazo lleva su **control positivo**: el mismo escenario sin la
//! modificación invasora se admite. Sin él, una guarda que rechazara todo
//! pasaría.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;

static NONCE: AtomicU64 = AtomicU64::new(0);

/// Cuatro semillas distintas. El rechazo no puede depender de ninguna.
const SEEDS: [&str; 4] = [
    "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
    "1f1e1d1c1b1a191817161514131211100f0e0d0c0b0a09080706050403020100",
    "a0a1a2a3a4a5a6a7a8a9aaabacadaeafb0b1b2b3b4b5b6b7b8b9babbbcbdbebf",
    "5555555555555555555555555555555555555555555555555555555555555555",
];

const SEQUENCE_LENGTH: u64 = 300;

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-s32-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Cien `G`, cien bases de `(CA)50` entre la 101 y la 200, y cien `G`.
fn fasta_bases() -> String {
    let bases = "G".repeat(100) + &"CA".repeat(50) + &"G".repeat(100);
    assert_eq!(bases.len() as u64, SEQUENCE_LENGTH);
    bases
}

/// **El alelo de la referencia está entre los declarados**, que es lo que hace
/// que el sorteo pueda no aplicar nada y, con ello, lo que hacía que el
/// escenario se aceptara o no según la semilla.
const ALLELES_INCLUDING_THE_REFERENCE_ONE: &str = r#"
[[reference.repeat_loci.alleles]]
units = 50
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 45
frequency = 0.5
"#;

fn locus(start: u64, end: u64) -> String {
    format!(
        r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = {start}
end = {end}
motif = "CA"
{ALLELES_INCLUDING_THE_REFERENCE_ONE}"#
    )
}

/// Ejecuta el escenario entero: si algo se rechaza, se rechaza antes de
/// publicar nada.
fn run(seed: &str, matrix: &str, loci: &str, events: &str) -> Result<(), Diagnostic> {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    fs::write(
        root.join("reference.fa"),
        format!(">chr1\n{}\n", fasta_bases()),
    )
    .expect("FASTA");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{seed}"
{matrix}
[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {SEQUENCE_LENGTH}
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"
{events}
"#
        ),
    )
    .expect("configuración");
    let outcome = load_and_validate_config(&config).and_then(|effective| {
        execute_validated_config(&effective, &root.join("published")).map(|_| ())
    });
    fs::remove_dir_all(&root).ok();
    outcome
}

/// Matriz suministrada con una fila en la posición indicada.
fn supplied_matrix(position_1based: u64) -> String {
    let path = temporary("tsv");
    fs::write(
        &path,
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\tsource_id\n\
             chr1\t{position_1based}\t{}\tT\tP1\t0.5\tfixed\n",
            base_at(position_1based)
        ),
    )
    .expect("matriz");
    format!(
        "\n[population_snvs]\nmode = \"supplied\"\ninput_path = \"{}\"\n",
        path.display()
    )
}

/// Matriz generada con `count` variantes repartidas por la secuencia.
fn generated_matrix(count: u64) -> String {
    format!(
        r#"
[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = {count}
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
"#
    )
}

/// La base que la referencia tiene de verdad en esa posición.
///
/// Se lee de la secuencia y no se deduce de una regla escrita a mano: un
/// cambio de base declarado con la base equivocada lo rechaza otra guarda
/// distinta, y entonces esta prueba estaría midiendo esa otra.
fn base_at(position_1based: u64) -> char {
    fasta_bases()
        .chars()
        .nth((position_1based - 1) as usize)
        .expect("la posición cae dentro de la secuencia")
}

fn direct_snp(event_id: &str, position_1based: u64) -> String {
    let reference_base = base_at(position_1based);
    format!(
        r#"
[[direct_events]]
event_id = "{event_id}"
individual_id = "I1"
phase = 0
kind = "snp"
contig = "chr1"
position_1based = {position_1based}
ref = "{reference_base}"
alt = "T"
"#
    )
}

fn rendered(error: &Diagnostic) -> String {
    error.render_human()
}

/// **Un cambio de base declarado dentro del tracto se rechaza, con cualquier
/// semilla**, y el mismo escenario sin él se admite.
#[test]
fn a_direct_base_change_inside_the_tract_is_rejected_with_any_seed() {
    for seed in SEEDS {
        // Control positivo: sin el cambio, el escenario corre entero.
        run(seed, &supplied_matrix(50), &locus(101, 200), "").unwrap_or_else(|error| {
            panic!(
                "el escenario sin invasión tiene que correr: {}",
                rendered(&error)
            )
        });

        let error = run(
            seed,
            &supplied_matrix(50),
            &locus(101, 200),
            &direct_snp("snp-inside", 150),
        )
        .expect_err("un cambio de base dentro del tracto rechaza el escenario");
        let text = rendered(&error);
        assert!(text.contains("STR-1"), "{text}");
        assert!(text.contains("snp-inside"), "{text}");
        assert!(text.contains("chr1:150"), "{text}");
    }
}

/// **El diagnóstico es legible.**
///
/// Ni comillas escapadas ni nombres de variantes de Rust: era literalmente lo
/// que el usuario recibía ante el error de configuración más probable de esta
/// funcionalidad.
#[test]
fn the_rejection_reads_like_spanish_and_not_like_a_debug_dump() {
    let error = run(
        SEEDS[0],
        &supplied_matrix(50),
        &locus(101, 200),
        &direct_snp("snp-inside", 150),
    )
    .expect_err("rechazo");
    let text = rendered(&error);
    assert!(!text.contains("\\\""), "comillas escapadas: {text}");
    assert!(
        !text.contains("RepeatIntervalAlreadyModified") && !text.contains("Diagnostic {"),
        "nombre de variante de Rust en el mensaje: {text}"
    );
    assert!(
        text.contains("cae dentro del tracto declarado"),
        "el mensaje dice qué pasa: {text}"
    );
    assert!(
        text.contains("interrupciones"),
        "y dice qué hacer en su lugar: {text}"
    );
}

/// Un cambio **justo fuera** del tracto se admite, y **justo dentro** se
/// rechaza.
///
/// La pareja está pegada al borde —la última base antes del tracto y la
/// primera del tracto— para que lo único que las separe sea el borde. Sin la
/// mitad de fuera, una guarda que rechazara cualquier cambio en la secuencia
/// pasaría igual.
#[test]
fn the_boundary_of_the_tract_is_where_the_rejection_starts() {
    for seed in SEEDS {
        run(
            seed,
            &supplied_matrix(50),
            &locus(101, 200),
            &direct_snp("snp-outside", 100),
        )
        .expect("la base 100 es la última de fuera");

        let error = run(
            seed,
            &supplied_matrix(50),
            &locus(101, 200),
            &direct_snp("snp-first", 101),
        )
        .expect_err("la base 101 es la primera del tracto");
        assert!(rendered(&error).contains("chr1:101"), "{error:?}");
    }
}

/// **Una fila de la matriz suministrada dentro del tracto se rechaza.**
///
/// Es un camino distinto del de los cambios directos: la matriz no pasa por la
/// lista de eventos declarados. Su control positivo es la misma fila fuera.
#[test]
fn a_supplied_matrix_row_inside_the_tract_is_rejected() {
    run(SEEDS[0], &supplied_matrix(50), &locus(101, 200), "")
        .expect("la fila en la base 50 está fuera del tracto");

    let error = run(SEEDS[0], &supplied_matrix(150), &locus(101, 200), "")
        .expect_err("la fila en la base 150 está dentro");
    let text = rendered(&error);
    assert!(text.contains("STR-1"), "{text}");
    // El identificador de la variante lo deriva el programa —la tabla del
    // usuario no lo lleva—, de modo que lo que el mensaje tiene que dar es la
    // procedencia y la coordenada, que es lo que el usuario puede buscar.
    assert!(text.contains("de la matriz en chr1:150"), "{text}");
}

/// **Las variantes generadas tampoco pueden caer dentro de un tracto.**
///
/// El locus cubre la secuencia entera, de modo que **cualquier** posición que
/// el generador elija cae dentro: el rechazo no depende de la suerte y la
/// prueba no es probabilística. El control positivo es el mismo escenario sin
/// locus declarado, que corre entero.
#[test]
fn generated_variants_inside_a_declared_tract_are_rejected() {
    let whole_sequence = format!(
        r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
{ALLELES_INCLUDING_THE_REFERENCE_ONE}"#
    );
    for seed in SEEDS {
        // Sin locus declarado, treinta variantes generadas corren.
        run(seed, &generated_matrix(30), "", "")
            .unwrap_or_else(|error| panic!("sin locus el escenario corre: {}", rendered(&error)));
    }
    // Con el locus y suficientes variantes, alguna cae dentro con cualquiera de
    // las semillas: el tracto ocupa un tercio de la secuencia y se piden
    // ochenta posiciones distintas.
    for seed in SEEDS {
        let error = run(seed, &generated_matrix(80), &whole_sequence, "")
            .expect_err("con ochenta variantes sobre trescientas bases alguna cae dentro");
        let text = rendered(&error);
        assert!(text.contains("STR-1"), "{text}");
        assert!(text.contains("de la matriz"), "{text}");
    }
}

/// **Un cambio directo llamado igual que un locus se rechaza antes de
/// simular.**
///
/// El identificador del evento que coloca un alelo de repetición **es el del
/// propio locus**, de modo que el choque haría ambiguo el origen de la
/// modificación. Antes, la ejecución corría entera y moría al final con un
/// error de formato sobre un documento producido, en vez de con un error de
/// configuración sobre lo que el usuario escribió.
#[test]
fn a_direct_event_named_like_the_locus_is_rejected_before_simulating() {
    // Control positivo: el mismo cambio, fuera del tracto y con otro nombre.
    run(
        SEEDS[0],
        &supplied_matrix(50),
        &locus(101, 200),
        &direct_snp("otro-nombre", 80),
    )
    .expect("un cambio fuera del tracto y con otro nombre se admite");

    let error = run(
        SEEDS[0],
        &supplied_matrix(50),
        &locus(101, 200),
        &direct_snp("STR-1", 80),
    )
    .expect_err("el cambio directo lleva el nombre del locus");
    let text = rendered(&error);
    assert_eq!(error.code(), "GPS-CONF-008", "{text}");
    assert!(text.contains("mismo identificador que el locus"), "{text}");
}

// ---------------------------------------------------------------------------
// La deriva de coordenadas: una inserción y una deleción declaran su punto en
// **coordenadas vigentes**, no en las de la referencia
// ---------------------------------------------------------------------------
//
// # Por qué estas tres pruebas van juntas y ninguna sobra
//
// Un cambio de base y una fila de la matriz declaran su posición en
// coordenadas de la **referencia**, el mismo espacio del intervalo del locus.
// Una inserción o una deleción, no: su punto es el de la secuencia **después**
// de los cambios anteriores de esa misma copia. Con un solo indel los dos
// espacios coinciden y no se nota; encadenando dos se separan, y entonces
// compararlos contra el intervalo se equivoca **en las dos direcciones**.
//
// Las tres pruebas cubren las tres situaciones, y la del medio —la que se
// admite— es la que impide «arreglarlo» rechazando de más:
//
// 1. sin deriva, un indel dentro del tracto se rechaza;
// 2. con deriva hacia fuera, un indel que **no** toca el tracto se admite;
// 3. con deriva hacia dentro, un indel que **sí** lo toca se rechaza.
//
// Las tres usan la misma inserción previa de sesenta bases, de modo que lo
// único que las separa es dónde acaba cayendo el segundo indel.

/// Sesenta bases insertadas en el punto 50, que es el desplazamiento que
/// separa las coordenadas vigentes de las de la referencia.
const DRIFT: u64 = 60;

fn insertion(event_id: &str, point: u64, length: usize) -> String {
    format!(
        r#"
[[direct_events]]
event_id = "{event_id}"
individual_id = "I1"
phase = 0
kind = "insertion"
contig = "chr1"
point = {point}
bases = "{}"
"#,
        "T".repeat(length)
    )
}

/// **Sin deriva, una inserción dentro del tracto se rechaza.**
///
/// Es el caso corriente y el control positivo de los otros dos: sin él, una
/// implementación que no mirara nunca los indeles pasaría las tres.
#[test]
fn an_insertion_inside_the_tract_is_rejected_when_no_previous_indel_shifts_it() {
    for seed in SEEDS {
        let error = run(
            seed,
            &supplied_matrix(50),
            &locus(101, 200),
            &insertion("ins-plain", 150, 4),
        )
        .expect_err("el punto 150 cae dentro del tracto");
        let text = rendered(&error);
        assert!(text.contains("STR-1"), "{text}");
        assert!(text.contains("no llega intacto"), "{text}");
        assert!(
            text.contains("ins-plain"),
            "el mensaje nombra el cambio: {text}"
        );
    }
}

/// **Con deriva hacia dentro, la invasión se caza igual.**
///
/// El punto vigente 210 son las 150 de la referencia, dentro del tracto. Una
/// comprobación que comparase el punto declarado contra el intervalo vería
/// 210, lo daría por fuera y **dejaría escapar la invasión**: publicaría una
/// verdad que dice que esa copia conserva el tracto de la referencia cuando no
/// lo conserva.
#[test]
fn an_insertion_that_drifts_into_the_tract_is_still_rejected() {
    for seed in SEEDS {
        let error = run(
            seed,
            &supplied_matrix(50),
            &locus(101, 200),
            &format!(
                "{}{}",
                insertion("ins-before", 50, DRIFT as usize),
                insertion("ins-drifted-in", 200 + DRIFT / 6, 4)
            ),
        )
        .expect_err("el punto vigente 210 son las 150 de la referencia");
        let text = rendered(&error);
        assert!(text.contains("STR-1"), "{text}");
        assert!(text.contains("no llega intacto"), "{text}");
        assert!(text.contains("ins-drifted-in"), "{text}");
    }
}

/// **Con deriva hacia fuera, el escenario se admite.**
///
/// El punto vigente 150 son las 90 de la referencia, once bases **antes** del
/// tracto. Es un escenario legítimo —un indel antes del microsatélite y otro
/// después— y rechazarlo diciendo que «cae dentro del tracto» sería afirmar
/// algo falso y dar un consejo imposible de seguir.
///
/// *Qué hace que pueda fallar:* la prueba anterior usa exactamente la misma
/// inserción previa y el mismo locus, y sí se rechaza. Lo único que cambia
/// entre las dos es el punto del segundo indel.
#[test]
fn an_insertion_that_drifts_out_of_the_tract_is_admitted() {
    for seed in SEEDS {
        run(
            seed,
            &supplied_matrix(50),
            &locus(101, 200),
            &format!(
                "{}{}",
                insertion("ins-before", 50, DRIFT as usize),
                insertion("ins-drifted-out", 150, 4)
            ),
        )
        .unwrap_or_else(|error| {
            panic!(
                "el punto vigente 150 son las 90 de la referencia, fuera del tracto: {}",
                rendered(&error)
            )
        });
    }
}
