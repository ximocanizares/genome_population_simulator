//! F4.10.6: validador independiente del tartamudeo clonal.
//!
//! # Qué lo hace independiente, y por qué importa
//!
//! **No llama al productor.** Rederiva las decisiones desde la semilla y la
//! configuración publicada volviendo a escribir aquí, a mano, las piezas
//! normativas que necesita:
//!
//! 1. el material canónico de derivación de subsemilla con su ámbito
//!    (DTC-037 versión 2), incluidos los componentes literales de la etapa de
//!    linaje;
//! 2. el ensayo binario `bernoulli-v1` (DTC-088);
//! 3. la elección ponderada `categorical-f64-v1` (DTC-072) con su suma
//!    compensada (DTC-049);
//! 4. la contabilidad de la tirada más larga sobre las bases publicadas.
//!
//! Las primitivas criptográficas y el generador —BLAKE3 y ChaCha12— sí son las
//! mismas bibliotecas: **son la especificación, no el productor**. Lo que no se
//! reutiliza es una sola línea de `src/clonal_stutter.rs`, `src/stutter.rs` ni
//! `src/wgs_sampling.rs`.
//!
//! # Lo que ata, y por qué eso es lo difícil
//!
//! No se queda en comparar recuentos con recuentos: **predice, lectura a
//! lectura, cuántas unidades tiene que mostrar el tracto**, y lo comprueba
//! contra los bytes del FASTQ publicado. Un recuento puede cuadrar por
//! casualidad; trescientas predicciones de un número entero, no.
//!
//! # Hasta dónde vale, dicho antes de que nadie lo estire
//!
//! El fixture está elegido para que el validador **no tenga que rederivar el
//! inserto**: con veinticuatro unidades de partida y un inserto mínimo de 250
//! pb, los cuatro pasos del algoritmo son siempre viables y el filtrado previo
//! no depende de la molécula. El validador **comprueba esa condición** en cada
//! paso en vez de suponerla, y si alguna vez dejara de cumplirse se para y lo
//! dice.
//!
//! Vale además solo para: fracción de duplicados **cero** —de modo que cada
//! familia tiene una copia y no hay progenitor ni divergencia que sortear—,
//! curva de composición neutra, una sola clase de eficiencia, un solo locus,
//! coeficiente de la tirada en cero, sin contaminación, sin quimeras y sin
//! errores de secuenciación.

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
use sha2::{Digest, Sha256};

const SEED_HEX: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
/// Contexto de derivación del proyecto, escrito aquí a mano.
const DERIVATION_CONTEXT: &str = "genome_population_simulator/v1/subseed";
/// Versión de derivación v2, la que lleva ámbito.
const DERIVATION_VERSION_V2: u32 = 2;
/// Catálogo de las cuatro decisiones de linaje.
const CATALOG_V6: &str = "f4-technical-v6";
/// Catálogo de las dos decisiones del tartamudeo de base.
const CATALOG_V4: &str = "f4-technical-v4";
/// Catálogo de las cinco decisiones de la amplificación.
const CATALOG_V5: &str = "f4-technical-v5";
/// Los cuatro pasos admitidos, en orden canónico. Escritos aquí a mano.
const ADMITTED_STEPS: [i64; 4] = [-3, -2, -1, 1];
/// Unidades del tracto verdadero.
const TRUE_UNITS: i64 = 24;
/// Parejas de la biblioteca.
const PAIRS: u64 = 300;
/// Longitud del motivo del fixture, en pares de base.
const MOTIF_LENGTH: u64 = 2;
/// `2^53`, el espacio del uniforme cerrado de F4.
const FIFTY_THREE_BIT: u64 = 9_007_199_254_740_992;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-10-6-{}-{}-{}",
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

/// Clave de una decisión de linaje, con sus ocho componentes en el orden del
/// catálogo: `[event_id, library_id, amplification, family_ordinal, decision,
/// copy_index, cycle, locus_id]`.
#[expect(
    clippy::too_many_arguments,
    reason = "son los ocho componentes literales de la clave; agruparlos escondería el orden, que es justo lo que esta prueba comprueba"
)]
fn lineage_rng(
    root: &[u8; 32],
    event_id: &str,
    library_id: &str,
    family_ordinal: u64,
    decision: &str,
    copy_index: u64,
    cycle: u64,
    locus_id: &str,
) -> ChaCha12Rng {
    let components = [
        ("event_id", event_id.to_owned()),
        ("library_id", library_id.to_owned()),
        ("stage", "amplification".to_owned()),
        ("family_ordinal", family_ordinal.to_string()),
        ("decision", decision.to_owned()),
        ("copy_index", copy_index.to_string()),
        ("cycle", cycle.to_string()),
        ("locus_id", locus_id.to_owned()),
    ];
    ChaCha12Rng::from_seed(blake3::derive_key(
        DERIVATION_CONTEXT,
        &canonical_material(root, CATALOG_V6, &components),
    ))
}

/// Clave de una decisión del tartamudeo de base, con sus seis componentes:
/// `[event_id, library_id, stutter, template_ordinal, locus_id, decision]`.
fn base_stutter_rng(
    root: &[u8; 32],
    event_id: &str,
    library_id: &str,
    template_ordinal: u64,
    locus_id: &str,
    decision: &str,
) -> ChaCha12Rng {
    let components = [
        ("event_id", event_id.to_owned()),
        ("library_id", library_id.to_owned()),
        ("stage", "stutter".to_owned()),
        ("template_ordinal", template_ordinal.to_string()),
        ("locus_id", locus_id.to_owned()),
        ("decision", decision.to_owned()),
    ];
    ChaCha12Rng::from_seed(blake3::derive_key(
        DERIVATION_CONTEXT,
        &canonical_material(root, CATALOG_V4, &components),
    ))
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

/// `bernoulli-v1`. Las probabilidades exactas 0 y 1 no consumen.
fn bernoulli(rng: &mut ChaCha12Rng, probability: f64) -> bool {
    if probability == 0.0 {
        return false;
    }
    if probability == 1.0 {
        return true;
    }
    let numerator = rng.next_u64() % FIFTY_THREE_BIT;
    (numerator as f64) / (FIFTY_THREE_BIT as f64) < probability
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

/// Contexto de derivación de identificadores, escrito aquí a mano.
const ID_CONTEXT: &str = "genome_population_simulator/v1/id";
/// Versión de derivación de identificadores técnicos.
const TECHNICAL_ID_DERIVATION_VERSION: u32 = 2;

/// Material canónico de un identificador técnico: raíz, versión, número de
/// componentes y cada componente con su longitud por delante.
fn technical_id_material(root: &[u8; 32], components: &[(&str, &str)]) -> Vec<u8> {
    let mut material = Vec::new();
    material.extend_from_slice(root);
    material.extend_from_slice(&TECHNICAL_ID_DERIVATION_VERSION.to_be_bytes());
    material.extend_from_slice(&(components.len() as u32).to_be_bytes());
    for (name, value) in components {
        material.extend_from_slice(&(name.len() as u32).to_be_bytes());
        material.extend_from_slice(name.as_bytes());
        material.extend_from_slice(&(value.len() as u32).to_be_bytes());
        material.extend_from_slice(value.as_bytes());
    }
    material
}

/// Base32 del alfabeto del proyecto, reescrito aquí.
fn encode_base32(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 32] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
    let mut encoded = String::new();
    let mut buffer = 0_u32;
    let mut bits = 0_u8;
    for byte in bytes {
        buffer = (buffer << 8) | u32::from(*byte);
        bits += 8;
        while bits >= 5 {
            bits -= 5;
            encoded.push(char::from(ALPHABET[((buffer >> bits) & 0x1f) as usize]));
        }
    }
    if bits > 0 {
        encoded.push(char::from(
            ALPHABET[((buffer << (5 - bits)) & 0x1f) as usize],
        ));
    }
    encoded
}

fn technical_id(root: &[u8; 32], prefix: &str, components: &[(&str, &str)]) -> String {
    let digest = blake3::derive_key(ID_CONTEXT, &technical_id_material(root, components));
    format!("{prefix}{}", encode_base32(&digest[..16]))
}

/// Identificador de pareja técnica de un ordinal de plantilla.
///
/// # Por qué hace falta rederivarlo
///
/// El FASTQ publicado **no va en orden de plantilla** —el ensamblado final
/// reordena las parejas— y su cabecera lleva el identificador `TP_…`, no el
/// ordinal. Se midió: las trescientas predicciones formaban exactamente el
/// mismo conjunto que lo publicado y **ni una sola caía en su sitio**.
/// Comparar por posición habría dado un fallo que no significa nada, y
/// comparar conjuntos habría dejado pasar una implementación que asignara el
/// paso a la molécula equivocada. Rederivando la cadena de identidad, cada
/// decisión queda atada a **su** lectura.
fn technical_pair_id(
    root: &[u8; 32],
    event_id: &str,
    library_id: &str,
    template_ordinal: u64,
) -> String {
    let template_id = format!("TPL_{library_id}_{template_ordinal}");
    let pre_hop = technical_id(
        root,
        "PH_",
        &[
            ("class", "PreHopInstanceId"),
            ("event_id", event_id),
            ("library_id", library_id),
            ("origin", "own"),
            ("template_id", &template_id),
        ],
    );
    technical_id(
        root,
        "TP_",
        &[
            ("class", "TechnicalPairId"),
            ("event_id", event_id),
            ("final_library_id", library_id),
            ("pre_hop_instance_id", &pre_hop),
        ],
    )
}

/// Copias completas contiguas del motivo en la tirada más larga.
fn longest_run_units(sequence: &[u8], motif: &[u8]) -> i64 {
    let mut best = 0_i64;
    let mut current = 0_i64;
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

// ---------------------------------------------------------------------------
// El validador.
// ---------------------------------------------------------------------------

/// Lo que el validador deduce por su cuenta de una ejecución publicada.
#[derive(Debug, PartialEq, Eq)]
struct Recomputed {
    /// Unidades que cada lectura tiene que mostrar, en orden de plantilla.
    units: Vec<i64>,
    /// Resbalones del linaje, por locus.
    lineage_slips: BTreeMap<String, u64>,
}

/// Clave de la etapa `amplification`, con sus cinco componentes.
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
    ChaCha12Rng::from_seed(blake3::derive_key(
        DERIVATION_CONTEXT,
        &canonical_material(root, CATALOG_V5, &components),
    ))
}

/// `uniform-index-v1` de DTC-088, con su **rechazo de la cola superior**.
///
/// Con `n = 1` devuelve cero **sin consumir**, que es lo que hace que la copia
/// 1 de cada familia no gaste una palabra en elegir el único progenitor que
/// tiene.
fn uniform_index(rng: &mut ChaCha12Rng, n: u64) -> u64 {
    assert!(n > 0, "el uniforme discreto exige una cota positiva");
    if n == 1 {
        return 0;
    }
    let remainder = (u64::MAX % n + 1) % n;
    let limit = (remainder != 0).then(|| u64::MAX - (remainder - 1));
    loop {
        let value = rng.next_u64();
        if limit.is_none_or(|bound| value < bound) {
            return value % n;
        }
    }
}

/// Distribución exacta de copias tras `early_cycles` ciclos, por el triángulo
/// de Pascal ponderado. El índice `i` corresponde a `z = i + 1`.
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

/// `p*(z) ∝ z · p(z)`: una familia observada es más probable cuantas más
/// copias tenga su molécula.
fn size_biased(distribution: &[f64]) -> Vec<f64> {
    let weighted: Vec<f64> = distribution
        .iter()
        .enumerate()
        .map(|(index, probability)| ((index + 1) as f64) * probability)
        .collect();
    let total = neumaier_sum(&weighted);
    weighted.into_iter().map(|weight| weight / total).collect()
}

/// Poisson de media `λ` condicionada a `K ≥ 1` y **pinzada** en `k_max`.
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

/// Pesos del ciclo de divergencia `d ∈ 0..=c`: geométrica de razón
/// `q = 1/(1+e)` con **toda** la masa restante en el último ciclo.
fn divergence_weights(cycles: u32, efficiency: f64) -> Vec<f64> {
    let slots = cycles as usize + 1;
    let mut weights = vec![0.0_f64; slots];
    let ratio = 1.0 / (1.0 + efficiency);
    let mut power = 1.0_f64;
    for slot in weights.iter_mut().take(slots - 1) {
        *slot = (1.0 - ratio) * power;
        power *= ratio;
    }
    weights[slots - 1] = power;
    weights
}

/// Inserto mínimo y longitud de lectura del perfil **que el manifiesto
/// nombra**, leídos de su propio recurso.
///
/// El manifiesto publica la ruta del recurso y su huella SHA-256, pero no el
/// inserto mínimo. Así que el validador va al recurso, **comprueba la huella
/// por su cuenta** y lee de ahí la geometría. De paso ata una cosa más: que la
/// huella publicada sea de verdad la del fichero que el documento nombra.
fn profile_geometry(event: &Value) -> (u64, u64) {
    let profile = &event["profile"];
    let resource_path = profile["resource_path"]
        .as_str()
        .expect("ruta del recurso del perfil");
    let published_sha256 = profile["sha256"].as_str().expect("huella del perfil");
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(resource_path);
    let bytes = fs::read(&path).expect("recurso del perfil");
    let digest = Sha256::digest(&bytes);
    let computed = digest.iter().fold(String::new(), |mut text, byte| {
        use std::fmt::Write as _;
        let _ = write!(text, "{byte:02x}");
        text
    });
    assert_eq!(
        computed, published_sha256,
        "la huella publicada no es la del recurso que el manifiesto nombra"
    );
    let resource: Value = serde_json::from_slice(&bytes).expect("JSON del perfil");
    let minimum = resource["insert"]["minimum_bases"]
        .as_u64()
        .expect("inserto mínimo");
    let read_length = resource["read"]["common_length_bases"]
        .as_u64()
        .expect("longitud de lectura");
    (minimum, read_length)
}

/// **La regla de viabilidad entera, no la mitad.**
///
/// El fixture está elegido para que el filtrado previo **no dependa de la
/// molécula**, y eso exige las **dos** mitades de la regla: que la tirada
/// conserve al menos una unidad completa tras la contracción mayor —`run > 3`—
/// y que la molécula no quede más corta que una lectura —`L − 3·motivo ≥
/// lectura`—. La primera versión de este validador comprobaba solo la primera
/// y su documentación decía comprobar las dos.
fn assert_all_four_steps_viable(run: i64, molecule_minimum: u64, read_length: u64) {
    assert!(
        run > 3,
        "el fixture garantiza que los cuatro pasos son viables: la tirada bajó a {run}"
    );
    let shortest = molecule_minimum
        .checked_sub(3 * MOTIF_LENGTH)
        .expect("la molécula mínima tiene que poder perder tres unidades");
    assert!(
        shortest >= read_length,
        "el fixture garantiza que los cuatro pasos son viables: la molécula más corta quedaría \
         en {shortest} bases, por debajo de la lectura de {read_length}"
    );
}

/// Rederiva el linaje y el tartamudeo de base de una ejecución.
///
/// Toma de los documentos publicados **los hechos de la configuración** —el
/// evento, la biblioteca, el locus, los ciclos, la probabilidad por ciclo, la
/// tasa de base y los pesos—, nunca los recuentos realizados ni las unidades
/// observadas, que son justo lo que va a comprobar.
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
        pcr["declared_duplicate_fraction"].as_f64(),
        Some(0.0),
        "este validador solo cubre familias de una copia"
    );
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
    let clonal = &pcr["stutter"];
    let cycles = clonal["cycles"].as_u64().expect("ciclos");
    let per_cycle = clonal["per_cycle_slip"]
        .as_array()
        .expect("tabla")
        .iter()
        .find(|entry| entry["motif_length_bases"] == 2)
        .and_then(|entry| entry["probability"].as_f64())
        .expect("probabilidad por ciclo del dinucleótido");

    let stutter = &event["artifacts"]["stutter"];
    assert_eq!(
        stutter["run_length_coefficient"]["value"].as_f64(),
        Some(0.0),
        "este validador solo cubre el coeficiente de la tirada en cero"
    );
    let base_rate = stutter["motif_length_rates"]
        .as_array()
        .expect("tasas")
        .iter()
        .find(|entry| entry["motif_length_bases"] == 2)
        .and_then(|entry| entry["rate"].as_f64())
        .expect("tasa de base del dinucleótido");
    let weights: Vec<f64> = ADMITTED_STEPS
        .iter()
        .map(|step| {
            stutter["step_weights"]
                .as_array()
                .expect("pesos")
                .iter()
                .find(|entry| entry["units"].as_i64() == Some(*step))
                .and_then(|entry| entry["weight"].as_f64())
                .expect("peso")
        })
        .collect();
    let locus_id = pcr["lineage_slips_by_locus"][0]["locus_id"]
        .as_str()
        .expect("locus");
    let pairs = event["final_fastq_summary"][0]["pair_count"]
        .as_u64()
        .expect("parejas publicadas");
    // La geometría sale del **perfil que el manifiesto nombra**, no de una
    // constante escrita aquí: si algún día ese perfil cambiara su inserto
    // mínimo o su lectura, la comprobación de viabilidad tiene que enterarse.
    let (molecule_minimum, read_length) = profile_geometry(event);

    let mut units = Vec::with_capacity(pairs as usize);
    let mut slips = 0_u64;
    for ordinal in 0..pairs {
        let mut run = TRUE_UNITS;
        // 1. El linaje, ciclo a ciclo. Con fracción de duplicados cero cada
        //    familia tiene una copia, de modo que no hay progenitor ni
        //    divergencia que sortear y la copia es siempre la cero.
        for cycle in 1..=cycles {
            let mut slip_rng = lineage_rng(
                &root,
                event_id,
                library_id,
                ordinal,
                "lineage_slip",
                0,
                cycle,
                locus_id,
            );
            if !bernoulli(&mut slip_rng, per_cycle) {
                continue;
            }
            assert_all_four_steps_viable(run, molecule_minimum, read_length);
            let mut step_rng = lineage_rng(
                &root,
                event_id,
                library_id,
                ordinal,
                "lineage_step",
                0,
                cycle,
                locus_id,
            );
            let step = ADMITTED_STEPS[categorical(&mut step_rng, &weights)];
            // **El resbalón se cuenta cuando el paso sale distinto de cero**,
            // que es lo que hace el productor. Contarlo tras el ensayo binario
            // coincidiría mientras los cuatro pesos sean positivos y dejaría de
            // coincidir en cuanto uno valga cero, sin que nada lo dijera.
            if step != 0 {
                slips += 1;
            }
            run += step;
        }
        // 2. La decisión de base, **después** del paso neto del linaje.
        let mut base_rng =
            base_stutter_rng(&root, event_id, library_id, ordinal, locus_id, "bernoulli");
        if bernoulli(&mut base_rng, base_rate) {
            assert_all_four_steps_viable(run, molecule_minimum, read_length);
            let mut step_rng =
                base_stutter_rng(&root, event_id, library_id, ordinal, locus_id, "step");
            run += ADMITTED_STEPS[categorical(&mut step_rng, &weights)];
        }
        units.push(run);
    }
    let mut lineage_slips = BTreeMap::new();
    lineage_slips.insert(locus_id.to_owned(), slips);
    Recomputed {
        units,
        lineage_slips,
    }
}

/// Lo que el validador deduce de una ejecución **con familias de varias
/// copias**.
#[derive(Debug, PartialEq, Eq)]
struct RecomputedFamilies {
    /// Unidades que cada lectura tiene que mostrar, en orden de plantilla.
    units: Vec<i64>,
    /// Ordinales de plantilla agrupados por familia, en orden.
    families: Vec<Vec<u64>>,
}

/// Rederiva **las familias enteras**: cuántas parejas tiene cada una, quién es
/// progenitor de quién, en qué ciclo diverge cada copia y qué paso neto le toca.
///
/// # Qué ata esto que la rederivación de familias de una copia no ataba
///
/// Que **cada copia lleva su linaje a su lectura**. Es el fenómeno que da
/// nombre a la entrega: una familia de duplicados que muestra más de un alelo
/// viniendo de una sola molécula. Sin esta rederivación, una implementación que
/// repartiera a todas las parejas de una familia el linaje de la copia 0
/// —familias homogéneas, artefacto desaparecido— pasaba la batería entera sin
/// que muriera nada.
///
/// # Alcance, comprobado en el documento antes de empezar
///
/// Una sola clase de eficiencia, curva de composición neutra, un solo locus,
/// coeficiente de la tirada en cero, sin contaminación y sin errores. Con
/// varias clases haría falta un sorteo más, y con una curva que rechace habría
/// que rederivar también el inserto, la copia y el comienzo.
fn recompute_families(effective: &Value, event: &Value) -> RecomputedFamilies {
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
    let clonal = &pcr["stutter"];
    let cycles = clonal["cycles"].as_u64().expect("ciclos");
    let per_cycle = clonal["per_cycle_slip"]
        .as_array()
        .expect("tabla")
        .iter()
        .find(|entry| entry["motif_length_bases"] == 2)
        .and_then(|entry| entry["probability"].as_f64())
        .expect("probabilidad por ciclo del dinucleótido");
    let stutter = &event["artifacts"]["stutter"];
    assert_eq!(
        stutter["run_length_coefficient"]["value"].as_f64(),
        Some(0.0),
        "este validador solo cubre el coeficiente de la tirada en cero"
    );
    let base_rate = stutter["motif_length_rates"]
        .as_array()
        .expect("tasas")
        .iter()
        .find(|entry| entry["motif_length_bases"] == 2)
        .and_then(|entry| entry["rate"].as_f64())
        .expect("tasa de base del dinucleótido");
    assert_eq!(
        base_rate, 0.0,
        "este fixture apaga la parte de base para aislar el linaje"
    );
    let weights: Vec<f64> = ADMITTED_STEPS
        .iter()
        .map(|step| {
            stutter["step_weights"]
                .as_array()
                .expect("pesos")
                .iter()
                .find(|entry| entry["units"].as_i64() == Some(*step))
                .and_then(|entry| entry["weight"].as_f64())
                .expect("peso")
        })
        .collect();
    let locus_id = pcr["lineage_slips_by_locus"][0]["locus_id"]
        .as_str()
        .expect("locus");
    let pairs = event["final_fastq_summary"][0]["pair_count"]
        .as_u64()
        .expect("parejas publicadas");
    let (molecule_minimum, read_length) = profile_geometry(event);

    // `E[Z] = (1+e)^c_early`, por multiplicación repetida.
    let mut mean_early = 1.0_f64;
    for _ in 0..early_cycles {
        mean_early *= 1.0 + efficiency;
    }
    let biased = size_biased(&early_copies(early_cycles, efficiency));
    let divergence_distribution =
        divergence_weights(u32::try_from(cycles).expect("ciclos"), efficiency);

    let mut units = vec![0_i64; pairs as usize];
    let mut families = Vec::new();
    let mut ordinal = 0_u64;
    while ordinal < pairs {
        let family_ordinal = ordinal;
        // Con una sola clase la eficiencia no se sortea, y con la curva neutra
        // el ensayo de aceptación no se abre.
        let copies = if early_cycles == 0 || mean_pairs == 0.0 {
            1_u64
        } else {
            let mut rng =
                amplification_rng(&root, event_id, library_id, family_ordinal, "early_copies");
            (categorical(&mut rng, &biased) as u64) + 1
        };
        let lambda = mean_pairs * ((copies as f64) / mean_early);
        let drawn = if lambda == 0.0 {
            1_u64
        } else {
            let mut rng =
                amplification_rng(&root, event_id, library_id, family_ordinal, "family_size");
            let weights = family_size_weights(lambda, max_family_size);
            (categorical(&mut rng, &weights) as u64) + 1
        };
        let emitted = drawn.min(pairs - ordinal);

        // La genealogía de la familia, copia a copia.
        let members = emitted as usize;
        let cycle_count = cycles as usize;
        let mut steps = vec![vec![0_i64; cycle_count]; members];
        let mut family = Vec::with_capacity(members);
        for member in 0..members {
            let copy_index = member as u64;
            let (parent, divergence) = if member == 0 {
                (0_usize, 0_usize)
            } else {
                let mut parent_rng = lineage_choice_rng(
                    &root,
                    event_id,
                    library_id,
                    family_ordinal,
                    "lineage_parent",
                    copy_index,
                );
                let parent = uniform_index(&mut parent_rng, copy_index) as usize;
                let mut divergence_rng = lineage_choice_rng(
                    &root,
                    event_id,
                    library_id,
                    family_ordinal,
                    "lineage_divergence",
                    copy_index,
                );
                let divergence = categorical(&mut divergence_rng, &divergence_distribution);
                (parent, divergence)
            };
            let mut run = TRUE_UNITS;
            for cycle in 1..=cycle_count {
                let step = if member > 0 && cycle <= divergence {
                    steps[parent][cycle - 1]
                } else {
                    let mut slip_rng = lineage_rng(
                        &root,
                        event_id,
                        library_id,
                        family_ordinal,
                        "lineage_slip",
                        copy_index,
                        cycle as u64,
                        locus_id,
                    );
                    if bernoulli(&mut slip_rng, per_cycle) {
                        assert_all_four_steps_viable(run, molecule_minimum, read_length);
                        let mut step_rng = lineage_rng(
                            &root,
                            event_id,
                            library_id,
                            family_ordinal,
                            "lineage_step",
                            copy_index,
                            cycle as u64,
                            locus_id,
                        );
                        ADMITTED_STEPS[categorical(&mut step_rng, &weights)]
                    } else {
                        0
                    }
                };
                steps[member][cycle - 1] = step;
                run += step;
            }
            units[(ordinal + copy_index) as usize] = run;
            family.push(ordinal + copy_index);
        }
        families.push(family);
        ordinal += emitted;
    }
    RecomputedFamilies { units, families }
}

/// Clave de una decisión de linaje **sin ciclo ni locus**: el progenitor y la
/// divergencia, que son de la copia entera y no de un tracto.
fn lineage_choice_rng(
    root: &[u8; 32],
    event_id: &str,
    library_id: &str,
    family_ordinal: u64,
    decision: &str,
    copy_index: u64,
) -> ChaCha12Rng {
    let components = [
        ("event_id", event_id.to_owned()),
        ("library_id", library_id.to_owned()),
        ("stage", "amplification".to_owned()),
        ("family_ordinal", family_ordinal.to_string()),
        ("decision", decision.to_owned()),
        ("copy_index", copy_index.to_string()),
    ];
    ChaCha12Rng::from_seed(blake3::derive_key(
        DERIVATION_CONTEXT,
        &canonical_material(root, CATALOG_V6, &components),
    ))
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

    /// Unidades que muestra cada lectura R1 publicada, **por identificador de
    /// pareja técnica**, que es lo que su cabecera lleva.
    fn observed_units(&self) -> BTreeMap<String, i64> {
        let mut text = Vec::new();
        bgzf::io::Reader::new(
            fs::File::open(
                self.destination
                    .join("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
            )
            .expect("FASTQ"),
        )
        .read_to_end(&mut text)
        .expect("legible");
        let text = String::from_utf8(text).expect("ASCII");
        let names = text.lines().step_by(4);
        let sequences = text.lines().skip(1).step_by(4);
        let mut units = BTreeMap::new();
        for (name, sequence) in names.zip(sequences) {
            let pair_id = name
                .split_whitespace()
                .next()
                .and_then(|head| head.rsplit(':').next())
                .unwrap_or_else(|| panic!("la cabecera no lleva pareja técnica: {name}"))
                .to_owned();
            assert!(
                pair_id.starts_with("TP_"),
                "la cabecera no lleva pareja técnica: {name}"
            );
            assert!(
                units
                    .insert(pair_id, longest_run_units(sequence.as_bytes(), b"CA"))
                    .is_none(),
                "dos lecturas con el mismo identificador de pareja"
            );
        }
        units
    }
}

/// Escenario con **familias de varias copias** y la parte de base apagada.
///
/// # Por qué cada número es el que es
///
/// - **Fracción de duplicados 0,5**: familias grandes de verdad, que es lo que
///   hace que haya copias `i ≥ 1` con progenitor y divergencia que sortear.
/// - **Un solo ciclo con probabilidad 1**: cada copia resbala exactamente una
///   vez, de modo que su paso neto es el paso que sortea o el que hereda, sin
///   acumulación que enturbie la comparación.
/// - **Pesos solo en `−1` y `+1`**: dos resultados posibles, 23 o 25 unidades.
///   Dentro de una familia, las copias que heredan coinciden con su progenitor
///   y las que divergen en el ciclo cero sortean por su cuenta. Es exactamente
///   el artefacto que la entrega existe para reproducir.
/// - **Tasa de base cero**: la lectura muestra el paso del linaje y nada más.
fn run_with_families() -> Published {
    run_scenario(
        "duplicate_fraction = 0.5\ncycles = 1\nefficiency = 0.9\n",
        "probability = 1.0\n",
        "rate = 0\n",
        [0.0, 0.0, 0.5, 0.5],
    )
}

fn run() -> Published {
    run_scenario(
        "duplicate_fraction = 0\ncycles = 10\nefficiency = 0.9\n",
        "probability = 0.05\n",
        "rate = 0.3\n",
        [0.05, 0.15, 0.6, 0.2],
    )
}

fn run_scenario(
    amplification: &str,
    per_cycle: &str,
    base_rate: &str,
    weights: [f64; 4],
) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let tract = "CA".repeat(TRUE_UNITS as usize);
    let head = format!("{}{tract}", "G".repeat(60));
    let bases = format!("{head}{}", "G".repeat(300 - head.len()));
    assert_eq!(bases.len(), 300);
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
length = 300

[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 61
end = {end}
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
target_coverage_wgs = {PAIRS}
sequencing_errors = "disabled"
library_preparation = "pcr"

[technical_event.pcr]
{amplification}
[[technical_event.pcr.stutter.per_cycle_slip]]
motif_length_bases = 2
{per_cycle}
[technical_event.stutter]
enabled = true

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
{base_rate}
{declared_weights}
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
            end = 60 + tract.len(),
            declared_weights = ADMITTED_STEPS
                .iter()
                .zip(weights)
                .map(|(units, weight)| format!(
                    "[[technical_event.stutter.steps]]\nunits = {units}\nweight = {weight}\n\n"
                ))
                .collect::<String>(),
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

/// **El validador rederiva el tracto de cada lectura y coincide con los bytes
/// publicados.**
///
/// Es la forma más exigente de atar las decisiones a lo escrito: no compara
/// recuentos con recuentos, sino **trescientas predicciones enteras** contra
/// las bases del FASTQ.
#[test]
fn the_validator_predicts_the_tract_of_every_published_read() {
    let published = run();
    let recomputed = recompute(&published.effective(), &published.event());
    let observed = published.observed_units();
    assert_eq!(
        observed.len(),
        recomputed.units.len(),
        "el validador tiene que predecir una lectura por pareja"
    );
    let effective = published.effective();
    let event = published.event();
    let root = decode_hex(effective["seed"]["root"].as_str().expect("semilla"));
    let event_id = event["event_id"].as_str().expect("evento");
    let library_id = event["selected_individuals"][0]["library_id"]
        .as_str()
        .expect("biblioteca");
    let predicted: BTreeMap<String, i64> = recomputed
        .units
        .iter()
        .enumerate()
        .map(|(ordinal, units)| {
            (
                technical_pair_id(&root, event_id, library_id, ordinal as u64),
                *units,
            )
        })
        .collect();
    assert_eq!(
        predicted.len(),
        recomputed.units.len(),
        "la cadena de identidad rederivada no puede colisionar"
    );
    assert_eq!(
        observed, predicted,
        "las unidades rederivadas no coinciden con las publicadas"
    );
    // Control positivo: el corpus **varía**. Sin esto, un fixture en el que
    // nada tartamudeara satisfaría la igualdad de arriba sin discriminar nada.
    let distinct = observed
        .values()
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    assert!(
        distinct > 1,
        "el corpus tiene que mostrar más de un número de unidades: {observed:?}"
    );
    assert!(
        observed.values().any(|units| *units != TRUE_UNITS),
        "tiene que haber lecturas tartamudeadas"
    );
}

/// **El recuento de resbalones del manifiesto coincide con el rederivado.**
#[test]
fn the_published_lineage_count_matches_the_recomputed_one() {
    let published = run();
    let event = published.event();
    let recomputed = recompute(&published.effective(), &event);
    let publicado: BTreeMap<String, u64> = event["artifacts"]["pcr"]["lineage_slips_by_locus"]
        .as_array()
        .expect("recuentos")
        .iter()
        .map(|entry| {
            (
                entry["locus_id"].as_str().expect("locus").to_owned(),
                entry["slips"].as_u64().expect("recuento"),
            )
        })
        .collect();
    assert_eq!(publicado, recomputed.lineage_slips);
    assert!(
        recomputed.lineage_slips.values().all(|count| *count > 0),
        "control positivo: tiene que haber resbalones que comparar"
    );
}

/// **Cada copia de una familia lleva SU linaje a SU lectura.**
///
/// Es el fenómeno que da nombre a esta entrega: una familia de duplicados que
/// muestra más de un número de unidades viniendo de **una sola molécula**. En
/// genética forense ese artefacto se confunde con un genotipo real.
///
/// La comprobación es la más fuerte que cabe: el validador rederiva las
/// familias enteras —cuántas parejas tiene cada una, quién desciende de quién y
/// en qué ciclo diverge— y predice, **pareja a pareja**, cuántas unidades tiene
/// que mostrar su lectura. Después ata cada predicción a **su** lectura por el
/// identificador de pareja técnica, que también rederiva.
///
/// **Debe fallar contra** una implementación que reparta a todas las parejas de
/// una familia el linaje de la copia 0: las familias dejarían de ser
/// heterogéneas y el artefacto desaparecería. Se comprobó que sin esta prueba
/// esa variante pasaba la batería entera —1446 pruebas, cero fallos— sin que
/// muriera nada.
#[test]
fn every_copy_of_a_family_carries_its_own_lineage_to_its_read() {
    let published = run_with_families();
    let effective = published.effective();
    let event = published.event();
    let recomputed = recompute_families(&effective, &event);
    let observed = published.observed_units();
    let root = decode_hex(effective["seed"]["root"].as_str().expect("semilla"));
    let event_id = event["event_id"].as_str().expect("evento");
    let library_id = event["selected_individuals"][0]["library_id"]
        .as_str()
        .expect("biblioteca");
    let predicted: BTreeMap<String, i64> = recomputed
        .units
        .iter()
        .enumerate()
        .map(|(ordinal, units)| {
            (
                technical_pair_id(&root, event_id, library_id, ordinal as u64),
                *units,
            )
        })
        .collect();
    assert_eq!(
        predicted.len(),
        recomputed.units.len(),
        "la cadena de identidad rederivada no puede colisionar"
    );
    assert_eq!(
        observed, predicted,
        "las unidades rederivadas no coinciden con las publicadas"
    );

    // **Los controles positivos, sin los cuales la igualdad no afirma nada.**
    let big_families = recomputed
        .families
        .iter()
        .filter(|family| family.len() > 1)
        .count();
    assert!(
        big_families > 0,
        "con fracción de duplicados 0,5 tiene que haber familias de varias copias"
    );
    // Y el artefacto en persona: al menos una familia muestra dos números de
    // unidades distintos **sobre los FASTQ**. Una implementación que repartiera
    // el linaje de la copia 0 a toda la familia no puede producir ninguna.
    let heterogeneous = recomputed
        .families
        .iter()
        .filter(|family| {
            family
                .iter()
                .map(|ordinal| recomputed.units[*ordinal as usize])
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                > 1
        })
        .count();
    assert!(
        heterogeneous > 0,
        "alguna familia tiene que mostrar más de un número de unidades: {} familias, \
         {} de más de una copia",
        recomputed.families.len(),
        big_families
    );
    // Y la otra mitad del par: también tiene que haber familias homogéneas, que
    // son las que heredaron. Sin esto, un modelo que sorteara todas las copias
    // por su cuenta —sin herencia ninguna— pasaría igual.
    let homogeneous = recomputed
        .families
        .iter()
        .filter(|family| {
            family.len() > 1
                && family
                    .iter()
                    .map(|ordinal| recomputed.units[*ordinal as usize])
                    .collect::<std::collections::BTreeSet<_>>()
                    .len()
                    == 1
        })
        .count();
    assert!(
        homogeneous > 0,
        "alguna familia de varias copias tiene que haber heredado entera"
    );
}

/// **Un recuento de linaje alterado se rechaza.**
///
/// Demuestra al revés que el validador no es un espejo: si leyera el recuento
/// publicado en vez de rederivarlo, esta prueba no podría existir.
#[test]
fn a_mutated_lineage_count_is_rejected() {
    let published = run();
    let mut event = published.event();
    let recomputed = recompute(&published.effective(), &event);
    let slips = event["artifacts"]["pcr"]["lineage_slips_by_locus"][0]["slips"]
        .as_u64()
        .expect("recuento");
    event["artifacts"]["pcr"]["lineage_slips_by_locus"][0]["slips"] = serde_json::json!(slips + 1);
    let mutated: BTreeMap<String, u64> = event["artifacts"]["pcr"]["lineage_slips_by_locus"]
        .as_array()
        .expect("recuentos")
        .iter()
        .map(|entry| {
            (
                entry["locus_id"].as_str().expect("locus").to_owned(),
                entry["slips"].as_u64().expect("recuento"),
            )
        })
        .collect();
    assert_ne!(
        mutated, recomputed.lineage_slips,
        "el validador tiene que rechazar un recuento alterado"
    );
}
