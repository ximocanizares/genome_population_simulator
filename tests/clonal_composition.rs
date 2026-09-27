//! F4.10.5: la comprobación aritmética de la composición.
//!
//! # Qué valida exactamente, dicho sin adornos
//!
//! La recomposición `1 − (1−r_sin)·(1−r_pcr) = r_con` es **exacta por
//! construcción**: la probabilidad por ciclo se derivó despejando justamente
//! eso, y acumularla sobre los ciclos devuelve `r_pcr` por definición. Por
//! tanto esta prueba **no comprueba que el modelo sea correcto**. Es
//! **tautológica respecto del modelo**.
//!
//! Lo que sí comprueba, y es su única razón de ser, es que **la
//! implementación compone como el diseño dice**. Caza los dos errores de
//! implementación que el orden de todo el plan existe para evitar:
//!
//! 1. una preparación amplificada entregada **sin la parte de base**, que
//!    daría 0,147 en el homopolímero y 0,0273 en el dinucleótido en vez de
//!    0,17 y 0,038;
//! 2. una que **sume** las dos columnas en vez de restarlas, que daría por
//!    encima de lo medido.
//!
//! Y **no** comprueba la genealogía: el fixture usa fracción de duplicados
//! cero, de modo que cada familia tiene una sola copia y no hay de quién
//! heredar. Se hace así a propósito, porque la aritmética de arriba supone
//! observaciones **independientes** y las copias de una familia no lo son.
//! Quien comprueba la herencia es `tests/clonal_lineage.rs`.
//!
//! # El tamaño muestral está calculado, no puesto a ojo
//!
//! Para distinguir `p₀` de la `p₁` que produciría una entrega sin la parte de
//! base hacen falta
//! `n ≈ (z_{α/2}·√(p₀(1−p₀)) + z_β·√(p₁(1−p₁)))² / (p₀−p₁)²` moléculas que
//! cubran el tracto. Con **α = 0,05 a dos colas y potencia del 90 %**: para el
//! homopolímero, `p₀ = 0,17` frente a `p₁ = 0,147`, **n ≈ 2700**; para el
//! dinucleótido, `p₀ = 0,038` frente a `p₁ = 0,0273`, **n ≈ 3000**.
//!
//! El fixture usa **n = 6000**, por encima de los dos, y de ahí sale la
//! tolerancia `z_{α/2}·√(p₀(1−p₀)/n)`, que vale **0,0095** para el
//! homopolímero y **0,0048** para el dinucleótido. Las dos son cómodamente
//! menores que las separaciones que hay que distinguir, 0,023 y 0,0107.
//!
//! # El sesgo por cancelación, que hay que declarar en vez de esconder
//!
//! La fracción que el oráculo cuenta es la de lecturas **que muestran otro
//! número de unidades**, no la de moléculas que resbalaron alguna vez. Dos
//! resbalones que se cancelan —un `+1` y un `−1`— dejan la lectura indistinguible
//! de una que no resbaló nunca. El efecto es pequeño y va **siempre hacia
//! abajo**: con la estructura de pasos de este fixture está por debajo de
//! 0,004 en el homopolímero y es despreciable en el dinucleótido. La
//! tolerancia lo absorbe, y por eso el intervalo se comprueba a dos lados en
//! vez de exigir la igualdad.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use noodles_bgzf as bgzf;
use serde_json::Value;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

/// Moléculas que cubren cada tracto. Ver la nota de tamaño muestral de arriba.
const MOLECULES: u64 = 6_000;
/// Unidades verdaderas del homopolímero y del dinucleótido.
///
/// **Son exactamente las longitudes de tracto de referencia**, y eso es lo que
/// convierte esta prueba en la demostración de que el modelo **no cuenta dos
/// veces** la longitud del tracto. A la referencia el factor de modulación vale
/// exactamente `1.0`, de modo que las fracciones que las lecturas muestran
/// tienen que salir **idénticas a las medidas** —0,17 y 0,038 con PCR, 0,027 y
/// 0,011 sin ella—. Si alguien multiplicara la curva encima de la tasa en vez
/// de normalizarla, o la sumara, estas cuatro cifras se moverían.
///
/// Quince unidades en el homopolímero son 15 pares de base; diecisiete en el
/// dinucleótido son 34. Las dos son el tope del tramo medido de su curva.
const TRUE_HOMOPOLYMER_UNITS: u64 = 15;
const TRUE_DINUCLEOTIDE_UNITS: u64 = 17;
/// `z` de la cola del 2,5 %.
const Z_TWO_SIDED: f64 = 1.959_964;

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-10-5-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
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
    fn event(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("events/reads-main/manifest.json"))
                .expect("manifiesto de evento"),
        )
        .expect("JSON")
    }

    /// Lecturas R1 publicadas, que son las que el oráculo mira.
    fn first_reads(&self) -> Vec<Vec<u8>> {
        read_sequences(
            &self
                .destination
                .join("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
        )
    }
}

fn read_sequences(path: &Path) -> Vec<Vec<u8>> {
    let mut text = Vec::new();
    bgzf::io::Reader::new(fs::File::open(path).expect("FASTQ"))
        .read_to_end(&mut text)
        .expect("legible");
    String::from_utf8(text)
        .expect("ASCII")
        .lines()
        .skip(1)
        .step_by(4)
        .map(|line| line.as_bytes().to_vec())
        .collect()
}

/// **El oráculo: copias completas contiguas del motivo en la tirada más larga
/// de una lectura.**
///
/// Está escrito **aquí**, en la prueba, y no reutiliza una sola línea del
/// productor: comparar la salida del productor con su propio contador sería un
/// oráculo circular. Mira la lectura entera, sin que le digan dónde está el
/// tracto, y eso es posible porque los flancos del fixture son todo `G`: el
/// único sitio donde el motivo puede repetirse es el tracto.
fn longest_run_units(sequence: &[u8], motif: &[u8]) -> u64 {
    let mut best = 0_u64;
    let mut current = 0_u64;
    let mut offset = 0_usize;
    while offset < sequence.len() {
        if offset + motif.len() <= sequence.len()
            && &sequence[offset..offset + motif.len()] == motif
        {
            current += 1;
            best = best.max(current);
            offset += motif.len();
        } else {
            current = 0;
            offset += 1;
        }
    }
    best
}

/// Fracción de lecturas cuyo tracto muestra otro número de unidades.
///
/// El denominador son las lecturas que **cubren el tracto entero**, que en
/// este fixture son todas: la geometría lo garantiza y la prueba lo comprueba
/// exigiendo que el denominador valga exactamente el número de moléculas.
fn stuttered_fraction(reads: &[Vec<u8>], motif: &[u8], truth: u64) -> (u64, u64) {
    let mut covered = 0_u64;
    let mut stuttered = 0_u64;
    for read in reads {
        let observed = longest_run_units(read, motif);
        if observed == 0 {
            continue;
        }
        covered += 1;
        if observed != truth {
            stuttered += 1;
        }
    }
    (stuttered, covered)
}

/// La referencia del fixture: dos tractos dentro de los primeros 150 pb de
/// cualquier molécula, y todo lo demás `G`.
///
/// # La aritmética que garantiza que R1 los ve enteros, siempre
///
/// La molécula es `copia[s, s+I)` con inserto `I ∈ [250, 300]` y comienzo
/// `s ∈ [0, 300−I]`, de modo que `s ≤ 50`. La lectura R1 es
/// `copia[s, s+150)`. Un tracto que ocupa `[a, a+T)` cabe entera en R1 para
/// **todo** comienzo posible exactamente cuando `a ≥ 50` y `a + T ≤ 150`.
/// Los dos tractos empiezan en 60 y en 79 y acaban, con la mayor expansión
/// posible —una unidad más—, en 76 y en 115. Ninguno de los dos números es
/// redondo por casualidad.
fn reference_bases() -> String {
    let homopolymer = "A".repeat(TRUE_HOMOPOLYMER_UNITS as usize);
    let dinucleotide = "CA".repeat(TRUE_DINUCLEOTIDE_UNITS as usize);
    let prefix = "G".repeat(60);
    let separator = "GGGG";
    let head = format!("{prefix}{homopolymer}{separator}{dinucleotide}");
    let bases = format!("{head}{}", "G".repeat(300 - head.len()));
    assert_eq!(bases.len(), 300);
    bases
}

/// Bloque del tartamudeo de base, con los pesos del ejemplo trabajado.
fn base_stutter() -> String {
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

/// Ejecuta el escenario con el bloque de preparación indicado.
///
/// Todo lo demás es idéntico entre las dos ejecuciones, que es lo que permite
/// atribuir la diferencia a la preparación y no a otra cosa.
fn run(preparation: &str) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let bases = reference_bases();
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    // El sitio de la matriz cae en el flanco, fuera de los dos tractos: dentro
    // invadiría un locus y la admisión lo rechazaría por un motivo ajeno.
    fs::write(
        root.join("population-snvs.tsv"),
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tG\tT\tP1\t0.5\n",
    )
    .expect("matriz");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
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
id = "STR-HOMO"
sequence = "chr1"
start = 61
end = 75
motif = "A"

[[reference.repeat_loci.alleles]]
sequence = "{homopolymer}"
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-DI"
sequence = "chr1"
start = 80
end = 113
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{dinucleotide}"
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
target_coverage_wgs = {coverage}
sequencing_errors = "disabled"
{preparation}{stutter}
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
            homopolymer = "A".repeat(TRUE_HOMOPOLYMER_UNITS as usize),
            dinucleotide = "CA".repeat(TRUE_DINUCLEOTIDE_UNITS as usize),
            coverage = MOLECULES,
            stutter = base_stutter(),
        ),
    )
    .expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Published { root, destination }
}

/// El bloque de la preparación amplificada, con fracción de duplicados cero.
///
/// **Cero no significa «sin PCR»**: la biblioteca sigue pasando por sus diez
/// ciclos y su tartamudeo clonal. Lo que la fracción cero quita son los
/// duplicados, y con ellos la correlación entre copias de una misma familia
/// que estropearía la aritmética de la muestra.
fn amplified() -> String {
    r#"library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0
cycles = 10
efficiency = 0.9
"#
    .to_owned()
}

/// Tolerancia de una fracción esperada con el tamaño muestral del fixture.
fn tolerance(expected: f64, sample: u64) -> f64 {
    Z_TWO_SIDED * (expected * (1.0 - expected) / (sample as f64)).sqrt()
}

fn assert_within(observed: f64, expected: f64, sample: u64, what: &str) {
    let tolerance = tolerance(expected, sample);
    assert!(
        (observed - expected).abs() <= tolerance,
        "{what}: observado {observed:.5}, esperado {expected} ± {tolerance:.5} con n = {sample}"
    );
}

/// **Con las dos partes activas, la fracción observada reproduce la fracción
/// medida con PCR.**
///
/// **Debe fallar contra** una entrega de la preparación amplificada sin la
/// parte de base, que daría 0,147 y 0,0273 en vez de 0,17 y 0,038; y contra
/// una que sume las dos columnas, que daría por encima de lo medido.
#[test]
fn the_amplified_run_reproduces_the_measured_with_pcr_fractions() {
    let published = run(&amplified());
    let reads = published.first_reads();
    assert_eq!(
        reads.len() as u64,
        MOLECULES,
        "el plan de cobertura tiene que dar exactamente las moléculas calculadas"
    );
    let (homopolymer, covered_homopolymer) =
        stuttered_fraction(&reads, b"A", TRUE_HOMOPOLYMER_UNITS);
    let (dinucleotide, covered_dinucleotide) =
        stuttered_fraction(&reads, b"CA", TRUE_DINUCLEOTIDE_UNITS);
    assert_eq!(
        covered_homopolymer, MOLECULES,
        "todas las lecturas tienen que cubrir el homopolímero entero"
    );
    assert_eq!(
        covered_dinucleotide, MOLECULES,
        "todas las lecturas tienen que cubrir el dinucleótido entero"
    );
    let observed_homopolymer = (homopolymer as f64) / (covered_homopolymer as f64);
    let observed_dinucleotide = (dinucleotide as f64) / (covered_dinucleotide as f64);
    assert_within(
        observed_homopolymer,
        0.17,
        MOLECULES,
        "homopolímero con PCR",
    );
    assert_within(
        observed_dinucleotide,
        0.038,
        MOLECULES,
        "dinucleótido con PCR",
    );
    // **La discriminación, escrita como aserción y no como confianza.** La
    // fracción que produciría una entrega sin la parte de base queda fuera del
    // intervalo; si algún día dejara de quedar fuera, la prueba habría dejado
    // de distinguir y hay que saberlo.
    assert!(
        (0.147_f64 - 0.17).abs() > tolerance(0.17, MOLECULES),
        "la variante sin la parte de base tiene que quedar fuera del intervalo"
    );
    assert!(
        (0.027_3_f64 - 0.038).abs() > tolerance(0.038, MOLECULES),
        "la variante sin la parte de base tiene que quedar fuera del intervalo"
    );
}

/// **Sin amplificar, las mismas dos fracciones caen en la columna sin PCR.**
///
/// Es el control que cierra la aritmética por el otro lado: sin él, unas
/// fracciones que estuvieran altas por cualquier motivo pasarían la prueba de
/// arriba y nadie vería que la parte de base sola ya daba 0,17.
#[test]
fn the_unamplified_run_reproduces_the_measured_pcr_free_fractions() {
    let published = run("");
    let reads = published.first_reads();
    assert_eq!(reads.len() as u64, MOLECULES);
    let (homopolymer, covered) = stuttered_fraction(&reads, b"A", TRUE_HOMOPOLYMER_UNITS);
    assert_eq!(covered, MOLECULES);
    assert_within(
        (homopolymer as f64) / (covered as f64),
        0.027,
        MOLECULES,
        "homopolímero sin PCR",
    );
    let (dinucleotide, covered) = stuttered_fraction(&reads, b"CA", TRUE_DINUCLEOTIDE_UNITS);
    assert_eq!(covered, MOLECULES);
    assert_within(
        (dinucleotide as f64) / (covered as f64),
        0.011,
        MOLECULES,
        "dinucleótido sin PCR",
    );
}

/// **El manifiesto de la ejecución amplificada publica el catálogo v6, la
/// derivación con su marca de tosquedad y los tres recuentos de linaje.**
///
/// Y los heredados valen **exactamente cero**, porque con fracción de
/// duplicados cero todas las familias tienen una sola copia y no hay de quién
/// heredar. **Debe fallar contra** un recuento de heredados que cuente todos
/// los resbalones.
#[test]
fn the_amplified_manifest_publishes_the_derivation_and_the_lineage_counts() {
    let published = run(&amplified());
    let event = published.event();
    assert_eq!(event["rng_catalog"], "f4-technical-v6");
    let pcr = &event["artifacts"]["pcr"];
    assert_eq!(pcr["rng_catalog"], "f4-technical-v6");
    let clonal = &pcr["stutter"];
    assert_eq!(clonal["algorithm"], "clonal-stutter-v2");
    assert_eq!(clonal["cycles"], 10);
    assert_eq!(
        clonal["coarse"], true,
        "la marca de tosquedad viaja siempre"
    );
    assert!(
        clonal["coarse_reason"]
            .as_str()
            .expect("motivo")
            .contains("different individuals"),
        "el motivo dice de dónde viene la tosquedad"
    );
    let entries = clonal["per_cycle_slip"].as_array().expect("tabla");
    assert_eq!(entries.len(), 6);
    // Las dos fracciones de entrada se publican junto a la derivada, que es lo
    // que permite rehacer la resta sin este programa.
    assert_eq!(entries[0]["with_pcr_fraction"], 0.17);
    assert_eq!(entries[0]["without_pcr_fraction"], 0.027);
    assert_eq!(entries[0]["source"], "derived");
    let probability = entries[0]["probability"].as_f64().expect("probabilidad");
    assert!(
        (probability - 0.0158).abs() < 5e-5,
        "la probabilidad por ciclo del homopolímero, calculada a mano, es 0,0158: {probability}"
    );

    let slips = pcr["lineage_slips_by_locus"]
        .as_array()
        .expect("resbalones");
    let inherited = pcr["lineage_inherited_slips_by_locus"]
        .as_array()
        .expect("heredados");
    assert_eq!(slips.len(), 2, "los dos loci aparecen aunque no resbalen");
    assert_eq!(inherited.len(), 2);
    let total_slips: u64 = slips
        .iter()
        .map(|entry| entry["slips"].as_u64().expect("recuento"))
        .sum();
    assert!(
        total_slips > 0,
        "control positivo: con diez ciclos tiene que haber resbalones"
    );
    for entry in inherited {
        assert_eq!(
            entry["slips"], 0,
            "con familias de tamaño 1 no hay de quién heredar: {entry}"
        );
    }
}
