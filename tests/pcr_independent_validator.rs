//! F4.9.12: validador independiente de la amplificación.
//!
//! # Qué lo hace independiente, y por qué importa
//!
//! **No llama al productor.** Rederiva las decisiones desde la semilla y la
//! configuración publicada volviendo a escribir aquí, a mano, las piezas
//! normativas que necesita:
//!
//! 1. el material canónico de derivación de subsemilla con su ámbito
//!    (DTC-037 versión 2),
//! 2. la elección ponderada `categorical-f64-v1` (DTC-072) con su suma
//!    compensada (DTC-049),
//! 3. la programación dinámica de los primeros ciclos y su sesgo por tamaño,
//! 4. la Poisson condicionada y **pinzada** del tamaño de familia.
//!
//! Las primitivas criptográficas y el generador —BLAKE3 y ChaCha12— sí son las
//! mismas bibliotecas, y la exponencial es la misma `libm`: **son la
//! especificación, no el productor**. Lo que no se reutiliza es una sola línea
//! de `src/pcr_amplification.rs` ni de `src/wgs_sampling.rs`.
//!
//! Si el validador leyera el histograma que el productor publica y lo comparara
//! consigo mismo sería un espejo y no podría fallar nunca. La prueba
//! `a_mutated_manifest_histogram_is_rejected` lo demuestra al revés: **muta el
//! histograma publicado y exige que el validador lo rechace**.
//!
//! # Hasta dónde vale, dicho antes de que nadie lo estire
//!
//! **Solo para curva de composición neutra y una sola clase de eficiencia.**
//! Con una curva que rechace, el factor de composición depende de la secuencia
//! de cada molécula y habría que rederivar también el sorteo del inserto, la
//! copia y el inicio; con varias clases, la eficiencia de cada familia sale de
//! un sorteo más. Los dos fixtures de aquí usan la curva neutra y una clase, y
//! el validador comprueba en el documento que así es antes de empezar.
//!
//! Y **solo para escenarios sin contaminación biológica**: con contaminación,
//! las parejas del fichero publicado dejan de ser las moléculas muestreadas y
//! la correspondencia de uno a uno que este validador usa se rompe.

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use noodles_bgzf as bgzf;
use rand_chacha::ChaCha12Rng;
use rand_chacha::rand_core::{Rng, SeedableRng};
use serde_json::Value;

const SEED_HEX: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
/// Contexto de derivación del proyecto, escrito aquí a mano.
const DERIVATION_CONTEXT: &str = "genome_population_simulator/v1/subseed";
/// Versión de derivación v2, la que lleva ámbito.
const DERIVATION_VERSION_V2: u32 = 2;
/// Catálogo de la etapa `amplification`.
const CATALOG_V5: &str = "f4-technical-v5";

static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-9-12-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

// ---------------------------------------------------------------------------
// Las piezas normativas, reescritas aquí.
// ---------------------------------------------------------------------------

/// Material canónico de derivación v2: raíz, versión, ámbito y componentes con
/// longitud por delante.
fn canonical_material(root: &[u8; 32], scope: &str, components: &[(&str, String)]) -> Vec<u8> {
    let mut material = Vec::new();
    material.extend_from_slice(root);
    material.extend_from_slice(&DERIVATION_VERSION_V2.to_be_bytes());
    let scope = scope.as_bytes();
    material.extend_from_slice(&(scope.len() as u32).to_be_bytes());
    material.extend_from_slice(scope);
    material.extend_from_slice(&(components.len() as u32).to_be_bytes());
    for (name, value) in components {
        material.extend_from_slice(&(name.len() as u32).to_be_bytes());
        material.extend_from_slice(name.as_bytes());
        material.extend_from_slice(&(value.len() as u32).to_be_bytes());
        material.extend_from_slice(value.as_bytes());
    }
    material
}

/// Clave de la etapa `amplification`, con sus cinco componentes en el orden del
/// catálogo. El **ordinal de familia** es el de plantilla de su primera pareja.
fn amplification_rng(
    root: &[u8; 32],
    event_id: &str,
    library_id: &str,
    family_ordinal: u64,
    decision: &str,
) -> ChaCha12Rng {
    let components = [
        ("event_id", event_id.to_owned()),
        ("library_id", library_id.to_owned()),
        ("stage", "amplification".to_owned()),
        ("family_ordinal", family_ordinal.to_string()),
        ("decision", decision.to_owned()),
    ];
    let seed = blake3::derive_key(
        DERIVATION_CONTEXT,
        &canonical_material(root, CATALOG_V5, &components),
    );
    ChaCha12Rng::from_seed(seed)
}

/// Suma compensada de Neumaier, el algoritmo `neumaier_sum_v1` de DTC-049.
fn neumaier_sum(values: &[f64]) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        let next = sum + value;
        if sum.abs() >= value.abs() {
            correction += (sum - next) + value;
        } else {
            correction += (value - next) + sum;
        }
        sum = next;
    }
    sum + correction
}

/// `categorical-f64-v1`. Con un solo peso positivo devuelve sin sortear.
fn categorical(rng: &mut ChaCha12Rng, weights: &[f64]) -> usize {
    let total = neumaier_sum(weights);
    let positives = weights.iter().filter(|weight| **weight > 0.0).count();
    if positives == 1 {
        return weights
            .iter()
            .position(|weight| *weight > 0.0)
            .expect("hay uno positivo");
    }
    let target = (((rng.next_u64() >> 12) as f64) + 0.5) * (1.0 / 4_503_599_627_370_496.0) * total;
    let mut cumulative = 0.0;
    let mut last_positive = 0;
    for (index, weight) in weights.iter().enumerate() {
        if *weight > 0.0 {
            last_positive = index;
        }
        cumulative += weight;
        if cumulative > target {
            return index;
        }
    }
    last_positive
}

/// Distribución exacta de copias tras `c_early` ciclos, por el triángulo de
/// Pascal ponderado. El índice `i` corresponde a `z = i + 1`.
fn early_copies(early_cycles: u32, efficiency: f64) -> Vec<f64> {
    let states = 1_usize << early_cycles;
    let mut distribution = vec![0.0_f64; states];
    distribution[0] = 1.0;
    if early_cycles == 0 {
        return distribution;
    }
    let failure = 1.0 - efficiency;
    let mut rows: Vec<Vec<f64>> = vec![vec![1.0]];
    for trials in 1..=states {
        let previous = rows[trials - 1].clone();
        let mut row = vec![0.0_f64; trials + 1];
        for (successes, slot) in row.iter_mut().enumerate() {
            let stay = if successes < previous.len() {
                previous[successes] * failure
            } else {
                0.0
            };
            let grow = if successes > 0 {
                previous[successes - 1] * efficiency
            } else {
                0.0
            };
            *slot = stay + grow;
        }
        rows.push(row);
    }
    for _ in 0..early_cycles {
        let mut next = vec![0.0_f64; states];
        for (index, probability) in distribution.iter().enumerate() {
            if *probability == 0.0 {
                continue;
            }
            let copies = index + 1;
            for (successes, weight) in rows[copies].iter().enumerate() {
                next[copies + successes - 1] += probability * weight;
            }
        }
        distribution = next;
    }
    distribution
}

/// `p*(z) ∝ z · p(z)`: una familia observada es más probable cuantas más copias
/// tenga su molécula.
fn size_biased(distribution: &[f64]) -> Vec<f64> {
    let weighted: Vec<f64> = distribution
        .iter()
        .enumerate()
        .map(|(index, probability)| ((index + 1) as f64) * probability)
        .collect();
    let total = neumaier_sum(&weighted);
    weighted.into_iter().map(|weight| weight / total).collect()
}

/// Poisson de media `λ` condicionada a `K ≥ 1` y **pinzada** en `k_max`: el
/// último peso recoge toda la masa restante.
fn family_size_weights(lambda: f64, max_family_size: usize) -> Vec<f64> {
    let mut weights = vec![0.0_f64; max_family_size];
    if lambda <= 0.0 || max_family_size == 1 {
        weights[0] = 1.0;
        return weights;
    }
    let exponential = libm::exp(-lambda);
    let denominator = 1.0 - exponential;
    let mut term = lambda * exponential / denominator;
    weights[0] = term;
    for (size, slot) in weights
        .iter_mut()
        .enumerate()
        .take(max_family_size - 1)
        .skip(1)
    {
        term = term * lambda / ((size + 1) as f64);
        *slot = term;
    }
    let head = neumaier_sum(&weights[..max_family_size - 1]);
    weights[max_family_size - 1] = (1.0 - head).max(0.0);
    weights
}

// ---------------------------------------------------------------------------
// El validador.
// ---------------------------------------------------------------------------

/// Lo que el validador deduce por su cuenta de una ejecución publicada.
#[derive(Debug, Eq, PartialEq)]
struct Recomputed {
    families: u64,
    pairs: u64,
    histogram: BTreeMap<u64, u64>,
}

/// Rederiva las familias de una ejecución.
///
/// Toma de los documentos publicados **los hechos de la configuración** —el
/// evento, la biblioteca, la `m` derivada, los ciclos explícitos, la eficiencia
/// y cuántas parejas hay—, nunca los recuentos realizados, que son justo lo que
/// va a comprobar.
fn recompute(effective: &Value, event: &Value) -> Recomputed {
    let root = decode_hex(
        effective["seed"]["root"]
            .as_str()
            .expect("semilla raíz publicada"),
    );
    let event_id = event["event_id"].as_str().expect("evento");
    let library_id = event["selected_individuals"][0]["library_id"]
        .as_str()
        .expect("biblioteca");
    let pcr = &event["artifacts"]["pcr"];
    // El alcance del validador, comprobado en el documento antes de empezar.
    assert_eq!(
        pcr["efficiency_classes"].as_array().map(Vec::len),
        Some(1),
        "este validador solo cubre una clase de eficiencia"
    );
    assert!(
        pcr["gc_relative_efficiency"]
            .as_array()
            .expect("curva")
            .iter()
            .all(|knot| knot["relative_efficiency"].as_f64() == Some(1.0)),
        "este validador solo cubre la curva neutra"
    );
    let efficiency = pcr["efficiency_classes"][0]["efficiency"]
        .as_f64()
        .expect("eficiencia");
    let early_cycles = u32::try_from(pcr["early_cycles"].as_u64().expect("ciclos explícitos"))
        .expect("ciclos explícitos");
    let mean_pairs = pcr["mean_pairs_per_molecule"].as_f64().expect("m derivada");
    let max_family_size =
        usize::try_from(pcr["max_family_size"].as_u64().expect("tope")).expect("tope");
    // Las parejas publicadas del fichero, que son un hecho del FASTQ y no un
    // recuento de la amplificación.
    let pairs = event["final_fastq_summary"][0]["pair_count"]
        .as_u64()
        .expect("parejas publicadas");

    // `E[Z] = (1+e)^c_early`, por multiplicación repetida.
    let mut mean_early = 1.0_f64;
    for _ in 0..early_cycles {
        mean_early *= 1.0 + efficiency;
    }
    let biased = size_biased(&early_copies(early_cycles, efficiency));

    let mut histogram: BTreeMap<u64, u64> = BTreeMap::new();
    let mut families = 0_u64;
    let mut ordinal = 0_u64;
    while ordinal < pairs {
        // Con una sola clase la eficiencia no se sortea.
        let copies = if early_cycles == 0 || mean_pairs == 0.0 {
            1_u64
        } else {
            let mut rng = amplification_rng(&root, event_id, library_id, ordinal, "early_copies");
            (categorical(&mut rng, &biased) as u64) + 1
        };
        // Con curva neutra el factor de composición vale exactamente 1 y el
        // ensayo de aceptación no se abre.
        let lambda = mean_pairs * ((copies as f64) / mean_early);
        let drawn = if lambda == 0.0 {
            1_u64
        } else {
            let mut rng = amplification_rng(&root, event_id, library_id, ordinal, "family_size");
            let weights = family_size_weights(lambda, max_family_size);
            (categorical(&mut rng, &weights) as u64) + 1
        };
        let emitted = drawn.min(pairs - ordinal);
        *histogram.entry(emitted).or_insert(0) += 1;
        families += 1;
        ordinal += emitted;
    }
    Recomputed {
        families,
        pairs,
        histogram,
    }
}

/// Lo que el manifiesto dice haber realizado.
fn published(event: &Value) -> Recomputed {
    let pcr = &event["artifacts"]["pcr"];
    let histogram = pcr["family_size_histogram"]
        .as_array()
        .expect("histograma")
        .iter()
        .map(|entry| {
            (
                entry["family_size"].as_u64().expect("tamaño"),
                entry["families"].as_u64().expect("familias"),
            )
        })
        .collect();
    Recomputed {
        families: pcr["families"].as_u64().expect("familias"),
        pairs: pcr["pairs"].as_u64().expect("parejas"),
        histogram,
    }
}

fn decode_hex(text: &str) -> [u8; 32] {
    let mut bytes = [0_u8; 32];
    for (index, slot) in bytes.iter_mut().enumerate() {
        *slot = u8::from_str_radix(&text[index * 2..index * 2 + 2], 16).expect("hexadecimal");
    }
    bytes
}

// ---------------------------------------------------------------------------
// El escenario y las pruebas.
// ---------------------------------------------------------------------------

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
    fn effective(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("effective-config.json")).expect("configuración"),
        )
        .expect("JSON")
    }

    fn event(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("events/reads-main/manifest.json"))
                .expect("manifiesto de evento"),
        )
        .expect("JSON")
    }

    /// Tamaños de los grupos de parejas idénticas en el FASTQ publicado.
    fn written_group_sizes(&self) -> BTreeMap<u64, u64> {
        let first = read_sequences(
            &self
                .destination
                .join("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
        );
        let second = read_sequences(
            &self
                .destination
                .join("outputs/fastq/reads-main/sample-A_R2.fastq.gz"),
        );
        let mut groups: BTreeMap<(Vec<u8>, Vec<u8>), u64> = BTreeMap::new();
        for pair in first.into_iter().zip(second) {
            *groups.entry(pair).or_insert(0) += 1;
        }
        let mut sizes: BTreeMap<u64, u64> = BTreeMap::new();
        for size in groups.values() {
            *sizes.entry(*size).or_insert(0) += 1;
        }
        sizes
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

/// Bases equilibradas y sin periodo corto, para que dos posiciones distintas no
/// produzcan la misma lectura.
fn bases(length: usize) -> String {
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            match (state >> 33) % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            }
        })
        .collect()
}

fn run(duplicate_fraction: &str, early_cycles: u32) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let sequence = bases(20_000);
    let reference_base = sequence.as_bytes()[19] as char;
    let alternate = if reference_base == 'T' { 'G' } else { 'T' };
    fs::write(root.join("reference.fa"), format!(">chr1\n{sequence}\n")).expect("FASTA");
    fs::write(
        root.join("population-snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\n\
             chr1\t20\t{reference_base}\t{alternate}\tP1\t0.5\n"
        ),
    )
    .expect("matriz");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED_HEX}"
ploidy = 1

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
length = 20000

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 8
sequencing_errors = "disabled"
library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = {duplicate_fraction}
cycles = 10
efficiency = 0.9
early_cycles = {early_cycles}

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
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

/// **La rederivación coincide con los recuentos publicados, en dos fixtures.**
#[test]
fn the_independent_recomputation_matches_the_published_aggregates() {
    for (fraction, early) in [("0.3", 3_u32), ("0.1", 0_u32)] {
        let published_run = run(fraction, early);
        let event = published_run.event();
        let recomputed = recompute(&published_run.effective(), &event);
        assert!(
            recomputed.families > 0,
            "tiene que haber familias, o la comparación no compara nada"
        );
        assert!(
            recomputed.pairs > recomputed.families,
            "con fracción {fraction} tiene que haber duplicados de verdad"
        );
        assert_eq!(recomputed, published(&event), "fixture {fraction}/{early}");
        // Control positivo de la distribución: si el histograma solo tuviera
        // un tamaño, la igualdad de arriba sería mucho más débil de lo que
        // parece.
        assert!(
            recomputed.histogram.len() > 1,
            "el histograma tiene que contener varios tamaños: {:?}",
            recomputed.histogram
        );
    }
}

/// **Un histograma mutado del manifiesto se rechaza.**
///
/// Es la prueba de que el validador **no es un espejo**: si leyera lo publicado
/// y lo comparara consigo mismo, mutar el manifiesto cambiaría los dos lados a
/// la vez y la comparación seguiría cuadrando.
#[test]
fn a_mutated_manifest_histogram_is_rejected() {
    let published_run = run("0.3", 3);
    let mut event = published_run.event();
    let original = published(&event);
    // Se mueve una familia de un tamaño a otro: el número de familias no
    // cambia, de modo que el rechazo tiene que venir del histograma y no de un
    // total que deje de cuadrar.
    let inflated = original
        .histogram
        .keys()
        .next()
        .copied()
        .expect("algún tamaño");
    event["artifacts"]["pcr"]["family_size_histogram"][0]["families"] =
        Value::from(original.histogram[&inflated] + 1);
    let recomputed = recompute(&published_run.effective(), &event);
    assert_ne!(
        recomputed,
        published(&event),
        "el validador tiene que rechazar un histograma que el productor no realizó"
    );
    // Y sigue coincidiendo con el recuento verdadero, para que el rechazo no
    // pueda venir de que el validador se haya roto.
    assert_eq!(recomputed, original);
}

/// **Las decisiones rederivadas están escritas en los bytes publicados.**
///
/// # Qué agujero cierra
///
/// Las pruebas anteriores rederivan el flujo de decisiones y lo comparan con
/// los recuentos del manifiesto. Pero ninguna mira lo que se **escribió**: un
/// productor que tomara exactamente las mismas decisiones y emitiera una sola
/// pareja por familia las satisfaría todas mientras el manifiesto siguiera
/// diciendo lo contrario. Aquí se cierra el lazo: los tamaños de los grupos de
/// parejas idénticas del FASTQ tienen que ser los tamaños de familia
/// rederivados.
///
/// *Vale para escenarios sin contaminación biológica, sin errores de
/// secuenciación y sin loci declarados*, y este lo es. Con errores activados
/// los duplicados dejan de ser idénticos y el recuento por bytes no
/// discriminaría nada.
#[test]
fn the_recomputed_family_sizes_match_the_groups_written_in_the_reads() {
    let published_run = run("0.3", 3);
    let recomputed = recompute(&published_run.effective(), &published_run.event());
    let written = published_run.written_group_sizes();
    assert_eq!(
        written, recomputed.histogram,
        "cada tamaño de familia rederivado tiene que aparecer escrito en el mismo número de \
         grupos de lecturas idénticas"
    );
    // Control positivo: el corpus tiene que contener grupos de más de una
    // pareja, o la igualdad se cumpliría con un histograma trivial.
    assert!(
        written.keys().any(|size| *size > 1),
        "el corpus tiene que mezclar familias de varios tamaños: {written:?}"
    );
}
