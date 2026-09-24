//! F4.8.8: validador independiente del tartamudeo.
//!
//! # Qué lo hace independiente, y por qué importa
//!
//! **No llama al productor.** Rederiva las decisiones desde la semilla y la
//! configuración publicada volviendo a escribir aquí, a mano, las tres piezas
//! normativas que necesita:
//!
//! 1. el material canónico de derivación de subsemilla con su ámbito
//!    (DTC-037 versión 2),
//! 2. el ensayo de Bernoulli `bernoulli-v1` (DTC-088),
//! 3. la elección ponderada `categorical-f64-v1` (DTC-072).
//!
//! Las primitivas criptográficas y el generador —BLAKE3 y ChaCha12— sí son las
//! mismas bibliotecas: **son la especificación, no el productor**. Lo que no se
//! reutiliza es una sola línea de `src/stutter.rs` ni de `src/randomness.rs`.
//!
//! Si el validador leyera los recuentos que el productor publica y los
//! comparara consigo mismos sería un espejo y no podría fallar nunca. La prueba
//! `a_mutated_manifest_count_is_rejected` lo demuestra al revés: **muta un
//! recuento del manifiesto y exige que el validador lo rechace**.
//!
//! # Hasta dónde vale, dicho antes de que nadie lo estire
//!
//! **Solo para escenarios sin contaminación biológica.** Los recuentos del
//! manifiesto cuentan las moléculas **muestreadas**; con contaminación, una
//! receptora retira parejas propias y publica parejas de una donante bajo su
//! biblioteca, de modo que las parejas del fichero publicado dejan de ser las
//! moléculas que se evaluaron y la correspondencia de uno a uno que este
//! validador usa se rompe. Los escenarios de aquí no declaran contaminación.
//!
//! # El fixture, y por qué permite rederivar sin simular
//!
//! La copia mide 280 pb y el tracto de diez bases está colocado de modo que
//! **toda** molécula lo contiene entero (la aritmética está en
//! `tests/stutter_pipeline.rs`). Por tanto la tirada más larga es la misma en
//! todas, 5 unidades, y con ella la probabilidad de resbalón y el conjunto de
//! pasos viables. Eso es lo que permite recorrer los ordinales de plantilla sin
//! volver a muestrear un solo fragmento.

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
const CATALOG_V4: &str = "f4-technical-v4";
const TWO_POW_53: u64 = 1 << 53;
/// Los cuatro pasos, en el orden canónico del catálogo.
const ADMITTED_STEPS: [i64; 4] = [-3, -2, -1, 1];

static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-8-8-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

// ---------------------------------------------------------------------------
// Las tres piezas normativas, reescritas aquí.
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

/// Clave de la etapa `stutter`, con sus seis componentes en el orden del
/// catálogo.
fn stutter_rng(
    root: &[u8; 32],
    event_id: &str,
    library_id: &str,
    ordinal: u64,
    locus_id: &str,
    decision: &str,
) -> ChaCha12Rng {
    let components = [
        ("event_id", event_id.to_owned()),
        ("library_id", library_id.to_owned()),
        ("stage", "stutter".to_owned()),
        ("template_ordinal", ordinal.to_string()),
        ("locus_id", locus_id.to_owned()),
        ("decision", decision.to_owned()),
    ];
    let seed = blake3::derive_key(
        DERIVATION_CONTEXT,
        &canonical_material(root, CATALOG_V4, &components),
    );
    ChaCha12Rng::from_seed(seed)
}

/// `uniform-index-v1` sobre `2^53`, con su rechazo de la cola alta.
fn uniform_index_53(rng: &mut ChaCha12Rng) -> u64 {
    let n = TWO_POW_53;
    let remainder = (u64::MAX % n + 1) % n;
    let limit = if remainder == 0 {
        None
    } else {
        Some(u64::MAX - (remainder - 1))
    };
    loop {
        let value: u64 = rng.next_u64();
        if limit.is_none_or(|bound| value < bound) {
            return value % n;
        }
    }
}

/// `bernoulli-v1`. Las probabilidades exactas 0 y 1 **no consumen** nada.
fn bernoulli(rng: &mut ChaCha12Rng, probability: f64) -> bool {
    if probability == 0.0 {
        return false;
    }
    if probability == 1.0 {
        return true;
    }
    (uniform_index_53(rng) as f64 / TWO_POW_53 as f64) < probability
}

/// `categorical-f64-v1`. Con un solo peso positivo devuelve sin sortear.
fn categorical(rng: &mut ChaCha12Rng, weights: &[f64]) -> usize {
    let total = weights.iter().sum::<f64>();
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

// ---------------------------------------------------------------------------
// El validador.
// ---------------------------------------------------------------------------

/// Lo que el validador deduce por su cuenta de una ejecución publicada.
#[derive(Debug, Eq, PartialEq)]
struct Recomputed {
    stuttered: u64,
    steps: std::collections::BTreeMap<i64, u64>,
}

/// Rederiva las decisiones de tartamudeo de una ejecución.
///
/// Toma de los documentos publicados **los hechos de la configuración** —el
/// evento, la biblioteca, el locus, la tasa efectiva, los pesos y cuántas
/// parejas hay—, nunca los recuentos realizados, que son justo lo que va a
/// comprobar.
fn recompute(effective: &Value, event: &Value, longest_run_units: u64) -> Recomputed {
    let root = decode_hex(
        effective["seed"]["root"]
            .as_str()
            .expect("semilla raíz publicada"),
    );
    let event_id = event["event_id"].as_str().expect("evento");
    let library_id = event["selected_individuals"][0]["library_id"]
        .as_str()
        .expect("biblioteca");
    let locus_id = event["artifacts"]["stutter"]["locus_parameters"][0]["locus_id"]
        .as_str()
        .expect("locus");
    let probability = event["artifacts"]["stutter"]["locus_parameters"][0]["base_rates"][0]["rate"]
        .as_f64()
        .expect("tasa base");
    let coefficient = event["artifacts"]["stutter"]["run_length_coefficient"]["value"]
        .as_f64()
        .expect("coeficiente");
    // La forma del modelo, reescrita aquí: dependencia **lineal por ciclo** de
    // la tirada ininterrumpida más larga, acotada a [0,1].
    let probability = (probability + coefficient * longest_run_units as f64).clamp(0.0, 1.0);
    let weights: Vec<f64> = ADMITTED_STEPS
        .iter()
        .map(|step| {
            event["artifacts"]["stutter"]["step_weights"]
                .as_array()
                .expect("pesos")
                .iter()
                .find(|entry| entry["units"].as_i64() == Some(*step))
                .and_then(|entry| entry["weight"].as_f64())
                .expect("peso del paso")
        })
        .collect();
    let pairs = event["final_fastq_summary"][0]["pair_count"]
        .as_u64()
        .expect("parejas publicadas");
    // Filtrado previo, reescrito aquí: un paso entra si deja al menos una
    // unidad completa en la tirada. La molécula nunca baja de la lectura en
    // este fixture, de modo que esa mitad del filtro no recorta nada.
    let candidates: Vec<i64> = ADMITTED_STEPS
        .into_iter()
        .filter(|step| *step > 0 || step.unsigned_abs() < longest_run_units)
        .collect();
    let candidate_weights: Vec<f64> = candidates
        .iter()
        .map(|step| {
            weights[ADMITTED_STEPS
                .iter()
                .position(|admitted| admitted == step)
                .expect("paso admitido")]
        })
        .collect();
    let mut stuttered = 0;
    let mut steps = std::collections::BTreeMap::new();
    for ordinal in 0..pairs {
        let mut trial = stutter_rng(&root, event_id, library_id, ordinal, locus_id, "bernoulli");
        if !bernoulli(&mut trial, probability) {
            continue;
        }
        let mut choice = stutter_rng(&root, event_id, library_id, ordinal, locus_id, "step");
        let step = candidates[categorical(&mut choice, &candidate_weights)];
        stuttered += 1;
        *steps.entry(step).or_insert(0) += 1;
    }
    Recomputed { stuttered, steps }
}

/// Lo que el manifiesto dice haber realizado.
fn published(event: &Value) -> Recomputed {
    let stutter = &event["artifacts"]["stutter"];
    let mut steps = std::collections::BTreeMap::new();
    for entry in stutter["steps_by_locus"].as_array().expect("distribución") {
        for step in entry["steps"].as_array().expect("pasos") {
            *steps
                .entry(step["units"].as_i64().expect("unidades"))
                .or_insert(0) += step["count"].as_u64().expect("recuento");
        }
    }
    Recomputed {
        stuttered: stutter["stuttered_tracts"].as_u64().expect("realizadas"),
        steps,
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

    fn mates(&self) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
        (
            read_sequences(
                &self
                    .destination
                    .join("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
            ),
            read_sequences(
                &self
                    .destination
                    .join("outputs/fastq/reads-main/sample-A_R2.fastq.gz"),
            ),
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

fn run(rate: f64) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let tract = "CA".repeat(5);
    let bases = format!("{}{tract}{}", "G".repeat(135), "G".repeat(135));
    assert_eq!(bases.len(), 280);
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
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
seed = "{SEED_HEX}"

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
length = 280

[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 136
end = 145
motif = "CA"

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
target_coverage_wgs = 100
sequencing_errors = "disabled"

[technical_event.stutter]
enabled = true

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
rate = {rate}

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

/// Copias completas contiguas del motivo, contadas **aquí**.
fn units(sequence: &[u8], motif: &[u8]) -> u64 {
    let mut best = 0;
    let mut current = 0;
    let mut offset = 0;
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

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect()
}

/// **El validador coincide con los recuentos agregados del manifiesto.**
#[test]
fn the_independent_recomputation_matches_the_published_aggregates() {
    let published_run = run(0.5);
    let event = published_run.event();
    let recomputed = recompute(&published_run.effective(), &event, 5);
    assert!(
        recomputed.stuttered > 0,
        "con la tasa a 0,5 tiene que haber tartamudeos, o la comparación no compara nada"
    );
    assert_eq!(recomputed, published(&event));
    // Control positivo de la distribución: si el corpus solo tuviera un paso,
    // la igualdad de arriba sería mucho más débil de lo que parece.
    assert!(
        recomputed.steps.len() > 1,
        "la distribución tiene que contener varios pasos: {:?}",
        recomputed.steps
    );
}

/// **Un recuento mutado del manifiesto se rechaza.**
///
/// Es la prueba de que el validador **no es un espejo**: si leyera lo publicado
/// y lo comparara consigo mismo, mutar el manifiesto cambiaría los dos lados a
/// la vez y la comparación seguiría cuadrando.
#[test]
fn a_mutated_manifest_count_is_rejected() {
    let published_run = run(0.5);
    let mut event = published_run.event();
    let original = published(&event);
    let inflated = original.stuttered + 1;
    event["artifacts"]["stutter"]["stuttered_tracts"] = Value::from(inflated);
    let recomputed = recompute(&published_run.effective(), &event, 5);
    assert_ne!(
        recomputed,
        published(&event),
        "el validador tiene que rechazar un recuento que el productor no realizó"
    );
    // Y sigue coincidiendo con el recuento verdadero, para que el rechazo no
    // pueda venir de que el validador se haya roto.
    assert_eq!(recomputed, original);
}

/// **Sobre los FASTQ finales, los dos mates de cada pareja muestran el mismo
/// número de unidades.**
#[test]
fn both_mates_of_every_published_pair_agree_on_the_unit_count() {
    let published_run = run(0.5);
    let (r1, r2) = published_run.mates();
    assert!(!r1.is_empty());
    let mut counts = std::collections::BTreeSet::new();
    for (first, second) in r1.iter().zip(&r2) {
        let forward = units(first, b"CA");
        let reverse = units(&reverse_complement(second), b"CA");
        assert_eq!(forward, reverse);
        counts.insert(forward);
    }
    assert!(
        counts.len() > 1,
        "con la tasa a 0,5 el corpus mezcla tractos tartamudeados y limpios: {counts:?}"
    );
}

/// Con la tasa a cero el validador **no consume aleatoriedad** y no encuentra
/// ningún tartamudeo, que es lo mismo que publica el manifiesto.
#[test]
fn a_zero_rate_recomputes_to_no_stutter_at_all() {
    let published_run = run(0.0);
    let event = published_run.event();
    let recomputed = recompute(&published_run.effective(), &event, 5);
    assert_eq!(recomputed.stuttered, 0);
    assert_eq!(recomputed, published(&event));
}

/// **La comprobación que ata las decisiones a los bytes.**
///
/// # Qué agujero cierra
///
/// Las demás pruebas de este fichero rederivan el **flujo de decisiones** y lo
/// comparan con los recuentos del manifiesto, que es lo que el plan pide. Pero
/// ninguna mira lo que se **escribió**: un productor que tomara exactamente las
/// mismas decisiones y luego aplicara el paso mal —por ejemplo dos veces por
/// molécula— las satisfaría todas. Se ha medido: esa edición sobrevive a las
/// otras cuatro pruebas.
///
/// Aquí se cierra el lazo. El tracto tiene **cinco unidades** y el paso se
/// elige entre `{−3,−2,−1,+1}`, de modo que una molécula publicada solo puede
/// mostrar 2, 3, 4, 5 o 6 unidades, y cada tamaño corresponde a un paso y solo
/// a uno. Se cuenta cada tamaño en los ficheros de lecturas y se compara con la
/// distribución de pasos **rederivada**, no con la publicada.
///
/// *Vale para escenarios sin contaminación biológica*, y este lo es: con
/// contaminación, las parejas publicadas dejan de ser las moléculas
/// muestreadas y la correspondencia de uno a uno se rompe.
#[test]
fn the_recomputed_decisions_match_the_units_written_in_the_reads() {
    const REFERENCE_UNITS: i64 = 5;
    let published_run = run(0.5);
    let recomputed = recompute(&published_run.effective(), &published_run.event(), 5);
    let (r1, _) = published_run.mates();
    assert!(!r1.is_empty());

    // Lo que de verdad hay escrito en las lecturas, contado aquí.
    let mut observed: std::collections::BTreeMap<i64, u64> = std::collections::BTreeMap::new();
    for read in &r1 {
        let seen = i64::try_from(units(read, b"CA")).expect("unidades");
        *observed.entry(seen - REFERENCE_UNITS).or_insert(0) += 1;
    }

    // Lo que las decisiones rederivadas predicen que debería haber, con el
    // paso cero para las moléculas que no resbalaron.
    let mut expected = recomputed.steps.clone();
    let stuttered: u64 = expected.values().sum();
    let quiet = r1.len() as u64 - stuttered;
    if quiet > 0 {
        expected.insert(0, quiet);
    }

    assert_eq!(
        observed, expected,
        "cada paso rederivado tiene que aparecer escrito en el mismo número de lecturas"
    );
    // Control positivo: el corpus tiene que contener más de un tamaño y algún
    // paso realizado. Sin esto, un corpus sin tartamudeo satisfaría la
    // igualdad con los dos lados vacíos salvo el cero.
    assert!(
        observed.len() > 2,
        "el corpus tiene que mezclar varios pasos: {observed:?}"
    );
    assert!(stuttered > 0);
}
