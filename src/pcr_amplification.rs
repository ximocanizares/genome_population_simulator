//! Preparación amplificada por PCR (F4.9).
//!
//! # Qué cambia, en una frase
//!
//! Hasta aquí el recorrido técnico producía **una pareja por fragmento**. Con
//! preparación amplificada produce **familias**: grupos de parejas que
//! descienden de la misma molécula original. De ahí salen los duplicados, la
//! sobrerrepresentación de unas moléculas frente a otras y el sesgo por
//! composición.
//!
//! # No se simula la población de moléculas
//!
//! Esa población crece como `2^ciclos` y no cabe en ninguna máquina. Se
//! muestrea **hacia atrás** desde las parejas que se van a publicar: como el
//! número de parejas `N` de una biblioteca ya es determinista desde la
//! cobertura, se generan familias hasta cubrirlas y cada familia cuesta un
//! número de operaciones acotado por las constantes de contrato de
//! [`crate::technical_profile`].
//!
//! # El usuario declara la fracción de duplicados, no el parámetro del modelo
//!
//! Un biólogo piensa en el porcentaje de duplicados que ve en su informe de
//! calidad, no en «parejas por molécula». El parámetro interno `m` se deriva
//! por **inversión numérica** de `f(m) = 1 − (1 − e^(−m))/m`, y esa inversión
//! tiene que dar el mismo resultado en cualquier máquina. Ver
//! [`invert_duplicate_fraction`].
//!
//! # Lo que la fracción declarada significa, y lo que no
//!
//! Es la fracción **esperada para una biblioteca homogénea**, no la realizada.
//! La realizada es **mayor** por tres motivos acumulativos: las clases de
//! eficiencia, la curva de composición y —el que se olvida— la **desigualdad
//! de los primeros ciclos**, que concentra las parejas en menos familias
//! aunque no haya ni clases ni curva. Basta con que los ciclos explícitos sean
//! más de cero, que es el caso por omisión. El manifiesto publica la realizada
//! aparte.
//!
//! Y **no** es la fracción que marca una herramienta de análisis contando
//! parejas que caen en el mismo sitio: esa incluye además las copias de la
//! química de agrupación, que no son PCR. Quien compare las dos sin saberlo no
//! las verá cuadrar nunca.

use crate::randomness::neumaier_sum;
use crate::technical_profile::{
    PCR_MAXIMUM_CYCLES, PCR_MAXIMUM_EARLY_CYCLES, PCR_MAXIMUM_EFFICIENCY_CLASSES,
    PCR_MAXIMUM_FAMILY_SIZE, PCR_MAXIMUM_GC_ATTEMPTS, PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE,
};

/// Versión del algoritmo, publicada en la configuración efectiva y en el
/// manifiesto. Cambiar cualquiera de sus piezas obliga a estrenar versión.
pub const PCR_AMPLIFICATION_ALGORITHM: &str = "pcr-amplification-v1";

/// Método de inversión de la fracción de duplicados, publicado con la `m`.
pub const PCR_DUPLICATE_FRACTION_INVERSION_METHOD: &str = "bisection-64-fixed-v1";

/// Iteraciones **fijas** de la bisección. No hay criterio de parada por
/// tolerancia: el número de pasos es el mismo para toda entrada distinta de
/// cero, y eso es lo que hace que el **coste** no dependa del objetivo.
pub const PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS: u32 = 64;

/// Valores por omisión del contrato, para los parámetros que sí tienen uno.
pub const PCR_DEFAULT_EARLY_CYCLES: u32 = 3;
pub const PCR_DEFAULT_MAX_FAMILY_SIZE: u32 = 64;
pub const PCR_DEFAULT_GC_ATTEMPTS: u32 = 8;

/// Fallo estable de la configuración o del modelo de amplificación.
///
/// Cada variante conserva el parámetro y el valor que lo provocó: es lo que
/// permite corregir la configuración sin adivinar.
#[derive(Clone, Debug, PartialEq)]
pub enum PcrAmplificationError {
    InvalidDuplicateFraction {
        value: f64,
    },
    /// La fracción pedida está por encima de la que `M_MAX` puede producir.
    /// **No se pinza en silencio**: se nombra el valor y el máximo alcanzable.
    UnreachableDuplicateFraction {
        value: f64,
        maximum_reachable: f64,
        maximum_mean_pairs_per_molecule: f64,
    },
    InvalidCycles {
        value: u32,
        maximum: u32,
    },
    InvalidEfficiency {
        value: f64,
    },
    InvalidEarlyCycles {
        value: u32,
        maximum: u32,
        cycles: u32,
    },
    InvalidMaxFamilySize {
        value: u32,
        maximum: u32,
    },
    InvalidGcAttempts {
        value: u32,
        maximum: u32,
    },
    InvalidGcCurveKnot {
        gc_fraction: f64,
        relative_efficiency: f64,
    },
    NonIncreasingGcCurveKnots {
        previous: f64,
        found: f64,
    },
    ShortGcCurve {
        knots: usize,
    },
    TooManyEfficiencyClasses {
        value: usize,
        maximum: u32,
    },
    InvalidEfficiencyClass {
        efficiency: f64,
        weight: f64,
    },
    EfficiencyClassWeightsNotNormalized {
        sum: f64,
        tolerance: f64,
    },
    NumericFailure,
}

// ---------------------------------------------------------------------------
// F4.9.4 — inversión numérica de la fracción de duplicados.
// ---------------------------------------------------------------------------

/// Fracción **esperada** de parejas duplicadas de una biblioteca homogénea.
///
/// `f(m) = 1 − (1 − e^(−m))/m`, continua y **estrictamente creciente** en `m`,
/// con `f(0) = 0` y límite 1. Esa monotonía es la precondición que hace única
/// la inversa, y una prueba la comprueba sobre una malla de cien puntos.
///
/// La exponencial procede **exclusivamente** de `libm = =0.2.16`, nunca de la
/// matemática estándar: es política vigente del proyecto (DTC-067, repetida en
/// DTC-072) y aquí además es lo que hace reproducible la bisección.
#[must_use]
pub fn expected_duplicate_fraction(mean_pairs_per_molecule: f64) -> f64 {
    if mean_pairs_per_molecule <= 0.0 {
        return 0.0;
    }
    1.0 - (1.0 - libm::exp(-mean_pairs_per_molecule)) / mean_pairs_per_molecule
}

/// Resultado de invertir la fracción declarada.
///
/// La reproducibilidad de una ejecución se ancla en la `m` **publicada**, no
/// en volver a derivarla: por eso el número de iteraciones viaja con ella.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuplicateFractionInversion {
    mean_pairs_per_molecule: f64,
    iterations: u32,
}

impl DuplicateFractionInversion {
    #[must_use]
    pub const fn mean_pairs_per_molecule(self) -> f64 {
        self.mean_pairs_per_molecule
    }

    #[must_use]
    pub const fn iterations(self) -> u32 {
        self.iterations
    }

    #[must_use]
    pub const fn method(self) -> &'static str {
        PCR_DUPLICATE_FRACTION_INVERSION_METHOD
    }
}

/// Deriva `m` desde la fracción de duplicados que el usuario declara.
///
/// # El método, con sus dos mitades
///
/// **Bisección sobre `[0, M_MAX]` con exactamente
/// [`PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS`] iteraciones**, sin criterio
/// de parada por tolerancia. Agotar la resolución de un `f64` en ese intervalo
/// cuesta entre 53 y 68 pasos según el tope, así que 64 bastan y sobran. No se
/// usa Newton, cuyo número de pasos depende del punto de partida.
///
/// Y la exponencial de [`expected_duplicate_fraction`] es la de `libm`, no la
/// del sistema. **Las dos mitades hacen falta**: un número fijo de iteraciones
/// garantiza el mismo **coste**, no el mismo **resultado**, porque la bisección
/// compara `f(mid)` con el objetivo y, si el último bit de la exponencial varía
/// entre plataformas, la comparación cae del otro lado y la bisección se va por
/// otra rama aunque el número de pasos sea idéntico. La revisión lo midió:
/// 194 875 de dos millones de exponenciales difieren en bits entre
/// implementaciones, y 2 591 de 20 000 fracciones objetivo acaban dando una `m`
/// distinta con 64 iteraciones fijas.
///
/// # Los dos casos límite
///
/// `f = 0` devuelve `m = 0` **por retorno anticipado, sin iterar**. Sin ese
/// retorno la bisección devuelve un residuo del orden de `10⁻¹⁸` en vez de
/// cero, y ese residuo **rompe la identidad byte a byte** de la reducción al
/// caso sin amplificar, porque un `m` que no es exactamente cero hace que el
/// tamaño de familia deje de ser trivial y se consuma aleatoriedad. Es la clase
/// de defecto que no da error y cambia todos los FASTQ.
///
/// Una `f` por encima de la que `M_MAX` puede producir se **rechaza**,
/// nombrando el valor y el máximo alcanzable.
///
/// # Cuidado al probar esto
///
/// Una prueba de reproducibilidad ejecutada **dos veces en la misma máquina**
/// no discrimina nada sobre el determinismo entre plataformas, porque la
/// exponencial del sistema también es determinista dentro de un sistema. Lo que
/// protege el contrato es la búsqueda negativa que comprueba que este módulo no
/// llama a la exponencial de la matemática estándar.
pub fn invert_duplicate_fraction(
    target: f64,
) -> Result<DuplicateFractionInversion, PcrAmplificationError> {
    if !target.is_finite() || !(0.0..1.0).contains(&target) {
        return Err(PcrAmplificationError::InvalidDuplicateFraction { value: target });
    }
    // Retorno anticipado exacto. No es un adorno: ver la sección de casos
    // límite de esta documentación.
    if target == 0.0 {
        return Ok(DuplicateFractionInversion {
            mean_pairs_per_molecule: 0.0,
            iterations: 0,
        });
    }
    let maximum_reachable = expected_duplicate_fraction(PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE);
    if target > maximum_reachable {
        return Err(PcrAmplificationError::UnreachableDuplicateFraction {
            value: target,
            maximum_reachable,
            maximum_mean_pairs_per_molecule: PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE,
        });
    }
    let mut low = 0.0_f64;
    let mut high = PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE;
    let mut iterations = 0_u32;
    while iterations < PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS {
        let midpoint = 0.5 * (low + high);
        if expected_duplicate_fraction(midpoint) < target {
            low = midpoint;
        } else {
            high = midpoint;
        }
        iterations += 1;
    }
    // Se devuelve el extremo superior: la invariante del bucle es
    // `f(high) >= target`, comprobada antes de entrar, así que la `m` devuelta
    // nunca se queda por debajo de lo pedido.
    Ok(DuplicateFractionInversion {
        mean_pairs_per_molecule: high,
        iterations,
    })
}

// ---------------------------------------------------------------------------
// F4.9.5 — copias tras los primeros ciclos, por programación dinámica exacta.
// ---------------------------------------------------------------------------

/// Distribución exacta del número de copias tras `early_cycles` ciclos.
///
/// Con eficiencia `e` y `z` copias presentes, el ciclo siguiente produce
/// `z + Binomial(z, e)` copias. Partiendo de una copia, la distribución tras
/// `c_early` ciclos se calcula **exactamente** por programación dinámica sobre
/// los estados `1..=2^c_early`. El índice `i` del vector devuelto corresponde a
/// `z = i + 1`.
///
/// Las probabilidades binomiales se acumulan por multiplicaciones sucesivas de
/// `e` y `1−e` en orden canónico —el triángulo de Pascal ponderado—, **sin**
/// `libm`: no hay ni factoriales ni potencias.
#[must_use]
pub fn early_copy_distribution(early_cycles: u32, efficiency: f64) -> Vec<f64> {
    EARLY_COPY_BUILDS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let states = 1_usize << early_cycles;
    let mut distribution = vec![0.0_f64; states];
    distribution[0] = 1.0;
    if early_cycles == 0 {
        return distribution;
    }
    let rows = binomial_rows(states, efficiency);
    for _ in 0..early_cycles {
        let mut next = vec![0.0_f64; states];
        for (index, probability) in distribution.iter().enumerate() {
            if *probability == 0.0 {
                continue;
            }
            let copies = index + 1;
            for (successes, weight) in rows[copies].iter().enumerate() {
                // `copies + successes` nunca supera `2^c_early`, porque cada
                // ciclo a lo sumo duplica y el vector se dimensionó para eso.
                next[copies + successes - 1] += probability * weight;
            }
        }
        distribution = next;
    }
    distribution
}

/// Cuántas veces se ha construido la distribución de los primeros ciclos.
///
/// # Por qué un contador de verdad y no una búsqueda en el texto
///
/// La primera guarda de esta propiedad era **léxica**: comprobaba que el
/// muestreador no nombrara estas funciones. Se midió que **se evade importando
/// la función con otro nombre**, con las dos baterías en verde; es decir, el
/// defecto que la guarda existía para impedir podía volver sin que muriera
/// nada. Un contador no se evade renombrando: cuenta construcciones, se llamen
/// como se llamen.
///
/// La tabla depende **solo** de los ciclos explícitos y de la eficiencia, y la
/// eficiencia toma a lo sumo `Q_MAX` valores, de modo que una tanda entera no
/// puede construirla más veces que clases tenga.
static EARLY_COPY_BUILDS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Lectura del contador de construcciones. Vive para que una prueba pueda
/// comprobar que la tabla se construye una vez por clase y no una por familia.
#[must_use]
pub fn early_copy_build_count() -> u64 {
    EARLY_COPY_BUILDS.load(std::sync::atomic::Ordering::Relaxed)
}

/// Triángulo de Pascal ponderado: `rows[n][k] = C(n,k)·e^k·(1−e)^(n−k)`.
///
/// Se construye fila a fila por la recurrencia
/// `rows[n][k] = rows[n−1][k]·(1−e) + rows[n−1][k−1]·e`, que es exactamente
/// «multiplicaciones sucesivas de `e` y `1−e` en orden canónico». No aparece
/// ningún coeficiente binomial explícito, que a partir de `C(64,32)` ya no
/// cabría exacto en un `f64`.
fn binomial_rows(max_trials: usize, efficiency: f64) -> Vec<Vec<f64>> {
    let failure = 1.0 - efficiency;
    let mut rows = Vec::with_capacity(max_trials + 1);
    rows.push(vec![1.0_f64]);
    for trials in 1..=max_trials {
        let previous = &rows[trials - 1];
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
    rows
}

/// Versión **sesgada por tamaño** de una distribución de copias.
///
/// Es la que se sortea, porque una familia observada es más probable cuantas
/// más copias tenga su molécula: `p*(z) ∝ z · p(z)`. Es también la razón por la
/// que `Z` se sortea **antes** que el tamaño de familia.
#[must_use]
pub fn size_biased(distribution: &[f64]) -> Vec<f64> {
    let weighted: Vec<f64> = distribution
        .iter()
        .enumerate()
        .map(|(index, probability)| ((index + 1) as f64) * probability)
        .collect();
    let total = neumaier_sum(weighted.iter().copied());
    if !total.is_finite() || total <= 0.0 {
        return weighted;
    }
    weighted.into_iter().map(|weight| weight / total).collect()
}

/// `E[Z] = (1+e)^c_early`, el factor determinista por el que se normaliza
/// `λ`. Se calcula por multiplicación repetida, sin `libm`.
#[must_use]
pub fn mean_early_copies(early_cycles: u32, efficiency: f64) -> f64 {
    let mut value = 1.0_f64;
    for _ in 0..early_cycles {
        value *= 1.0 + efficiency;
    }
    value
}

// ---------------------------------------------------------------------------
// F4.9.6 — tamaño de familia.
// ---------------------------------------------------------------------------

/// Pesos del tamaño de familia: Poisson de media `λ` **condicionada a `K ≥ 1`**
/// y **pinzada** en `max_family_size`.
///
/// El índice `i` corresponde a `k = i + 1`. Para `1 ≤ k < k_max`,
/// `p_k = (λ^k/k!)·e^(−λ)/(1 − e^(−λ))`, calculado por multiplicación sucesiva
/// (`p_1 = λ·e^(−λ)/(1−e^(−λ))`, `p_(k+1) = p_k·λ/(k+1)`) con la exponencial de
/// `libm`; y `p_(k_max)` recoge **toda** la masa restante de la cola, sumada
/// con la suma compensada de DTC-049.
///
/// # Por qué no `poisson-count-v1`
///
/// Su muestreo por rechazo no acota el número de sorteos y su cola no está
/// pinzada. El proyecto ya fijó que un recuento se **pinza** en vez de
/// rechazarse, «porque el rechazo no acota el número de sorteos y el
/// presupuesto de operaciones sí». Que conste aquí para que nadie lea DTC-067
/// como una prohibición: aquí no hay recuento de sucesos sobre oportunidades,
/// sino una **elección entre alternativas**, cuyo consumidor normativo es
/// `categorical-f64-v1`.
///
/// Con `λ = 0` el vector es `[1, 0, …]`, de modo que el sorteo categórico
/// devolvería sin consumir; el muestreador además ni lo llama.
#[must_use]
pub fn family_size_weights(lambda: f64, max_family_size: u32) -> Vec<f64> {
    let slots = max_family_size.max(1) as usize;
    let mut weights = vec![0.0_f64; slots];
    if !(lambda.is_finite() && lambda > 0.0) {
        weights[0] = 1.0;
        return weights;
    }
    let exponential = libm::exp(-lambda);
    let denominator = 1.0 - exponential;
    if !(denominator.is_finite() && denominator > 0.0) {
        weights[0] = 1.0;
        return weights;
    }
    if slots == 1 {
        weights[0] = 1.0;
        return weights;
    }
    let mut term = lambda * exponential / denominator;
    weights[0] = term;
    for (size, slot) in weights.iter_mut().enumerate().take(slots - 1).skip(1) {
        term = term * lambda / ((size + 1) as f64);
        *slot = term;
    }
    // El último recoge **toda** la masa de la cola, no el término que le
    // tocaría. Sin este pinzamiento la distribución perdería masa y el tamaño
    // máximo de familia dejaría de ser un tope.
    let head = neumaier_sum(weights[..slots - 1].iter().copied());
    weights[slots - 1] = (1.0 - head).max(0.0);
    weights
}

// ---------------------------------------------------------------------------
// F4.9.7 — composición de bases.
// ---------------------------------------------------------------------------

/// Curva de eficiencia relativa por composición, en tramos lineales.
#[derive(Clone, Debug, PartialEq)]
pub struct GcCurve {
    knots: Vec<(f64, f64)>,
}

impl Default for GcCurve {
    /// La curva neutra: eficiencia relativa 1 en todo el intervalo. Con ella
    /// el factor de composición vale **exactamente** `1.0` y el ensayo de
    /// aceptación no consume ni una palabra de generador.
    fn default() -> Self {
        Self {
            knots: vec![(0.0, 1.0), (1.0, 1.0)],
        }
    }
}

impl GcCurve {
    /// Construye la curva comprobando sus dominios: al menos dos nudos,
    /// fracciones **estrictamente crecientes** en `[0,1]` y eficiencias
    /// relativas en `(0,1]`.
    pub fn new(knots: Vec<(f64, f64)>) -> Result<Self, PcrAmplificationError> {
        if knots.len() < 2 {
            return Err(PcrAmplificationError::ShortGcCurve { knots: knots.len() });
        }
        let mut previous: Option<f64> = None;
        for (gc_fraction, relative_efficiency) in &knots {
            if !gc_fraction.is_finite()
                || !(0.0..=1.0).contains(gc_fraction)
                || !relative_efficiency.is_finite()
                || *relative_efficiency <= 0.0
                || *relative_efficiency > 1.0
            {
                return Err(PcrAmplificationError::InvalidGcCurveKnot {
                    gc_fraction: *gc_fraction,
                    relative_efficiency: *relative_efficiency,
                });
            }
            if let Some(previous) = previous
                && *gc_fraction <= previous
            {
                return Err(PcrAmplificationError::NonIncreasingGcCurveKnots {
                    previous,
                    found: *gc_fraction,
                });
            }
            previous = Some(*gc_fraction);
        }
        Ok(Self { knots })
    }

    #[must_use]
    pub fn knots(&self) -> &[(f64, f64)] {
        &self.knots
    }

    /// `true` cuando la curva no distingue ninguna composición.
    #[must_use]
    pub fn is_neutral(&self) -> bool {
        self.knots.iter().all(|(_, efficiency)| *efficiency == 1.0)
    }

    /// Eficiencia relativa en una fracción de GC, por interpolación lineal y
    /// pinzada fuera del intervalo de los nudos.
    #[must_use]
    pub fn evaluate(&self, gc_fraction: f64) -> f64 {
        let first = self.knots[0];
        if gc_fraction <= first.0 {
            return first.1;
        }
        let last = self.knots[self.knots.len() - 1];
        if gc_fraction >= last.0 {
            return last.1;
        }
        for window in self.knots.windows(2) {
            let (left_x, left_y) = window[0];
            let (right_x, right_y) = window[1];
            if gc_fraction <= right_x {
                let span = right_x - left_x;
                if span <= 0.0 {
                    return left_y;
                }
                return left_y + (right_y - left_y) * ((gc_fraction - left_x) / span);
            }
        }
        last.1
    }
}

/// Fracción de GC de un fragmento, **sobre A/C/G/T**.
///
/// Las `N` no cuentan ni en el numerador ni en el denominador; un fragmento sin
/// ninguna base A/C/G/T recibe `1.0`, que con cualquier curva es el valor que
/// no distorsiona nada.
#[must_use]
pub fn gc_fraction(sequence: &[u8]) -> f64 {
    let mut acgt = 0_u64;
    let mut gc = 0_u64;
    for base in sequence {
        match base {
            b'G' | b'C' | b'g' | b'c' => {
                acgt += 1;
                gc += 1;
            }
            b'A' | b'T' | b'a' | b't' => acgt += 1,
            _ => {}
        }
    }
    if acgt == 0 {
        return 1.0;
    }
    (gc as f64) / (acgt as f64)
}

/// Factor de composición de una molécula: `G = ((1 + e·g(GC)) / (1 + e))^c`.
///
/// Con la curva neutra `g = 1`, el cociente es **exactamente** `1.0` y la
/// potencia también, de modo que la aceptación no consume aleatoriedad. La
/// potencia se calcula por multiplicación repetida —el exponente es un entero
/// acotado por `C_MAX`—, así que tampoco aquí hace falta `libm`.
#[must_use]
pub fn composition_factor(efficiency: f64, relative_efficiency: f64, cycles: u32) -> f64 {
    let ratio = (1.0 + efficiency * relative_efficiency) / (1.0 + efficiency);
    let mut value = 1.0_f64;
    for _ in 0..cycles {
        value *= ratio;
    }
    value.clamp(0.0, 1.0)
}

// ---------------------------------------------------------------------------
// Parámetros efectivos.
// ---------------------------------------------------------------------------

/// Una clase de eficiencia declarada.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EfficiencyClass {
    efficiency: f64,
    weight: f64,
}

impl EfficiencyClass {
    #[must_use]
    pub const fn efficiency(self) -> f64 {
        self.efficiency
    }

    #[must_use]
    pub const fn weight(self) -> f64 {
        self.weight
    }
}

/// Distribución de copias tras los primeros ciclos de **una** clase de
/// eficiencia, ya sesgada por tamaño, y su media.
///
/// # Por qué está precalculada
///
/// La distribución depende **solo** de los ciclos explícitos y de la
/// eficiencia, y la eficiencia toma a lo sumo `Q_MAX` valores distintos. La
/// primera versión de este módulo la reconstruía **dentro de cada familia**, y
/// con ella el triángulo de Pascal de 65 filas: una ejecución la ejecutó 13 741
/// veces para publicar una cota que declaraba haberla calculado una sola vez
/// por clase. Extrapolado a un genoma humano a 30× eran del orden de cincuenta
/// minutos de reloj **y** una discrepancia entre el código y la fórmula
/// publicada que ninguna prueba podía ver.
///
/// Calcularla aquí no cambia un solo bit del resultado —es una función pura de
/// dos números— y hace que el código diga lo mismo que la cota.
#[derive(Clone, Debug, PartialEq)]
pub struct EarlyCopyTable {
    size_biased_weights: Vec<f64>,
    mean_copies: f64,
}

impl EarlyCopyTable {
    /// Pesos de la distribución **sesgada por tamaño**, que es la que se
    /// sortea: el índice `i` corresponde a `z = i + 1`.
    #[must_use]
    pub fn size_biased_weights(&self) -> &[f64] {
        &self.size_biased_weights
    }

    /// `E[Z] = (1+e)^c_early`, el factor determinista por el que se normaliza.
    #[must_use]
    pub const fn mean_copies(self_: &Self) -> f64 {
        self_.mean_copies
    }
}

/// Parámetros efectivos de la preparación amplificada.
///
/// Todos los dominios están ya comprobados contra las constantes de contrato:
/// quien recibe estos ajustes no vuelve a interpretarlos.
#[derive(Clone, Debug, PartialEq)]
pub struct PcrSettings {
    declared_duplicate_fraction: f64,
    inversion: DuplicateFractionInversion,
    cycles: u32,
    efficiency: f64,
    early_cycles: u32,
    max_family_size: u32,
    gc_max_attempts: u32,
    gc_curve: GcCurve,
    efficiency_classes: Vec<EfficiencyClass>,
    /// Una tabla por clase de eficiencia, en el mismo orden que las clases.
    early_copy_tables: Vec<EarlyCopyTable>,
    /// Unidades de operación que costó construir esas tablas. Es el término
    /// `Q × c_early × 4^c_early` de la cota, y se contabiliza **una vez**,
    /// porque una vez es lo que se paga.
    early_copy_build_units: u64,
    /// Tartamudeo clonal heredable por linaje, si el escenario declara loci de
    /// repetición.
    ///
    /// Va **dentro** de los ajustes de la amplificación y no aparte porque es
    /// la amplificación la que lo produce: sin ciclos de PCR no hay linaje en
    /// el que heredar nada. Sin loci declarados vale `None` y una biblioteca
    /// amplificada publica exactamente los mismos bytes que antes de F4.10.
    clonal: Option<crate::clonal_stutter::ClonalStutterSettings>,
}

/// Lo que el autor del escenario declara, ya leído del TOML pero sin validar.
#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredPcrParameters {
    pub duplicate_fraction: f64,
    pub cycles: u32,
    pub efficiency: f64,
    pub early_cycles: Option<u32>,
    pub max_family_size: Option<u32>,
    pub gc_max_attempts: Option<u32>,
    pub gc_curve: Option<Vec<(f64, f64)>>,
    pub efficiency_classes: Vec<(f64, f64)>,
}

impl PcrSettings {
    /// Valida los dominios y deriva `m`.
    ///
    /// `tolerance` es la tolerancia de normalización efectiva del escenario
    /// (DTC-049), la misma con la que se juzga cualquier otro reparto de pesos.
    pub fn new(
        declared: &DeclaredPcrParameters,
        tolerance: f64,
    ) -> Result<Self, PcrAmplificationError> {
        if declared.cycles == 0 || declared.cycles > PCR_MAXIMUM_CYCLES {
            return Err(PcrAmplificationError::InvalidCycles {
                value: declared.cycles,
                maximum: PCR_MAXIMUM_CYCLES,
            });
        }
        if !declared.efficiency.is_finite()
            || declared.efficiency <= 0.0
            || declared.efficiency > 1.0
        {
            return Err(PcrAmplificationError::InvalidEfficiency {
                value: declared.efficiency,
            });
        }
        let early_cycles = declared
            .early_cycles
            .unwrap_or_else(|| PCR_DEFAULT_EARLY_CYCLES.min(declared.cycles));
        if early_cycles > PCR_MAXIMUM_EARLY_CYCLES || early_cycles > declared.cycles {
            return Err(PcrAmplificationError::InvalidEarlyCycles {
                value: early_cycles,
                maximum: PCR_MAXIMUM_EARLY_CYCLES,
                cycles: declared.cycles,
            });
        }
        let max_family_size = declared
            .max_family_size
            .unwrap_or(PCR_DEFAULT_MAX_FAMILY_SIZE);
        if max_family_size == 0 || max_family_size > PCR_MAXIMUM_FAMILY_SIZE {
            return Err(PcrAmplificationError::InvalidMaxFamilySize {
                value: max_family_size,
                maximum: PCR_MAXIMUM_FAMILY_SIZE,
            });
        }
        let gc_max_attempts = declared.gc_max_attempts.unwrap_or(PCR_DEFAULT_GC_ATTEMPTS);
        if gc_max_attempts == 0 || gc_max_attempts > PCR_MAXIMUM_GC_ATTEMPTS {
            return Err(PcrAmplificationError::InvalidGcAttempts {
                value: gc_max_attempts,
                maximum: PCR_MAXIMUM_GC_ATTEMPTS,
            });
        }
        let gc_curve = match declared.gc_curve.clone() {
            Some(knots) => GcCurve::new(knots)?,
            None => GcCurve::default(),
        };
        let efficiency_classes = resolve_efficiency_classes(declared, tolerance)?;
        let inversion = invert_duplicate_fraction(declared.duplicate_fraction)?;
        // **Una sola vez por clase**, aquí y no dentro de cada familia.
        let mut early_copy_build_units = 0_u64;
        let early_copy_tables = efficiency_classes
            .iter()
            .map(|class| {
                let distribution = early_copy_distribution(early_cycles, class.efficiency());
                // El coste real de construirla: el triángulo de Pascal y la
                // propagación ciclo a ciclo. Se cuenta donde se paga.
                let states = 1_u64 << early_cycles;
                early_copy_build_units = early_copy_build_units.saturating_add(
                    u64::from(early_cycles).saturating_mul(states.saturating_mul(states)),
                );
                EarlyCopyTable {
                    size_biased_weights: size_biased(&distribution),
                    mean_copies: mean_early_copies(early_cycles, class.efficiency()),
                }
            })
            .collect();
        Ok(Self {
            declared_duplicate_fraction: declared.duplicate_fraction,
            inversion,
            cycles: declared.cycles,
            efficiency: declared.efficiency,
            early_cycles,
            max_family_size,
            gc_max_attempts,
            gc_curve,
            efficiency_classes,
            early_copy_tables,
            early_copy_build_units,
            clonal: None,
        })
    }

    /// Conecta el tartamudeo clonal, que solo existe cuando el escenario
    /// declara loci de repetición.
    ///
    /// Es un paso aparte de [`Self::new`] porque los dos hechos que lo
    /// deciden viven en sitios distintos del fichero: la amplificación se
    /// declara en el evento técnico y los microsatélites junto a la
    /// referencia. Solo cuando las dos configuraciones están resueltas puede
    /// saberse si hay linaje que heredar.
    #[must_use]
    pub fn with_clonal_stutter(
        mut self,
        clonal: crate::clonal_stutter::ClonalStutterSettings,
    ) -> Self {
        self.clonal = Some(clonal);
        self
    }

    /// Parámetros del tartamudeo clonal de este evento, si los hay.
    #[must_use]
    pub const fn clonal_stutter(&self) -> Option<&crate::clonal_stutter::ClonalStutterSettings> {
        self.clonal.as_ref()
    }

    /// Catálogo de ámbitos aleatorios del que salen las decisiones de **esta**
    /// amplificación: el v6 en cuanto abre las decisiones de linaje, y el v5
    /// mientras no lo haga.
    #[must_use]
    pub const fn rng_catalog(&self) -> &'static str {
        if self.clonal.is_some() {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V6
        } else {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V5
        }
    }

    /// Tabla precalculada de la clase de eficiencia con este índice.
    #[must_use]
    pub fn early_copy_table(&self, class_index: usize) -> Option<&EarlyCopyTable> {
        self.early_copy_tables.get(class_index)
    }

    /// Unidades que costó construir las tablas, para el contador de la cota.
    #[must_use]
    pub const fn early_copy_build_units(&self) -> u64 {
        self.early_copy_build_units
    }

    #[must_use]
    pub const fn declared_duplicate_fraction(&self) -> f64 {
        self.declared_duplicate_fraction
    }

    #[must_use]
    pub const fn inversion(&self) -> DuplicateFractionInversion {
        self.inversion
    }

    #[must_use]
    pub const fn mean_pairs_per_molecule(&self) -> f64 {
        self.inversion.mean_pairs_per_molecule()
    }

    #[must_use]
    pub const fn cycles(&self) -> u32 {
        self.cycles
    }

    #[must_use]
    pub const fn efficiency(&self) -> f64 {
        self.efficiency
    }

    #[must_use]
    pub const fn early_cycles(&self) -> u32 {
        self.early_cycles
    }

    #[must_use]
    pub const fn max_family_size(&self) -> u32 {
        self.max_family_size
    }

    #[must_use]
    pub const fn gc_max_attempts(&self) -> u32 {
        self.gc_max_attempts
    }

    #[must_use]
    pub const fn gc_curve(&self) -> &GcCurve {
        &self.gc_curve
    }

    #[must_use]
    pub fn efficiency_classes(&self) -> &[EfficiencyClass] {
        &self.efficiency_classes
    }

    /// Modo que el manifiesto publica para la composición.
    #[must_use]
    pub fn gc_mode(&self) -> &'static str {
        if self.gc_curve.is_neutral() {
            "neutral"
        } else {
            "curve"
        }
    }
}

/// Resuelve las clases de eficiencia declaradas.
///
/// # La consecuencia que hay que decir en voz alta
///
/// Sin clases declaradas hay **una sola**, con la eficiencia escalar de la
/// tabla. Con clases declaradas, esa escalar **deja de intervenir**: cada
/// familia usa la eficiencia de su clase, tanto en el factor de composición
/// como en la tabla de los primeros ciclos. La escalar sigue siendo
/// obligatoria —lo exige el contrato— y se publica, pero quien declare
/// eficiencia 0,9 y clases 0,95 y 0,5 verá publicado un 0,9 que no cambia
/// nada. Una prueba de extremo a extremo lo demuestra cambiando solo la
/// escalar y comprobando que los FASTQ no se mueven.
fn resolve_efficiency_classes(
    declared: &DeclaredPcrParameters,
    tolerance: f64,
) -> Result<Vec<EfficiencyClass>, PcrAmplificationError> {
    if declared.efficiency_classes.is_empty() {
        return Ok(vec![EfficiencyClass {
            efficiency: declared.efficiency,
            weight: 1.0,
        }]);
    }
    if declared.efficiency_classes.len() > PCR_MAXIMUM_EFFICIENCY_CLASSES as usize {
        return Err(PcrAmplificationError::TooManyEfficiencyClasses {
            value: declared.efficiency_classes.len(),
            maximum: PCR_MAXIMUM_EFFICIENCY_CLASSES,
        });
    }
    let mut classes = Vec::with_capacity(declared.efficiency_classes.len());
    for (efficiency, weight) in &declared.efficiency_classes {
        if !efficiency.is_finite()
            || *efficiency <= 0.0
            || *efficiency > 1.0
            || !weight.is_finite()
            || *weight < 0.0
        {
            return Err(PcrAmplificationError::InvalidEfficiencyClass {
                efficiency: *efficiency,
                weight: *weight,
            });
        }
        classes.push(EfficiencyClass {
            efficiency: *efficiency,
            weight: *weight,
        });
    }
    // Los pesos se juzgan con **la misma** suma compensada y la misma
    // tolerancia que cualquier otro reparto del proyecto (DTC-049).
    let sum = neumaier_sum(classes.iter().map(|class| class.weight()));
    if !sum.is_finite() || (sum - 1.0).abs() > tolerance {
        return Err(PcrAmplificationError::EfficiencyClassWeightsNotNormalized { sum, tolerance });
    }
    Ok(classes)
}

/// Plan de amplificación que se presta al muestreador PE.
///
/// El plan por omisión está **apagado**, de modo que un muestreador construido
/// sin él se comporta exactamente como antes de que la amplificación existiera.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PcrPlan {
    settings: Option<PcrSettings>,
}

impl PcrPlan {
    #[must_use]
    pub const fn new(settings: PcrSettings) -> Self {
        Self {
            settings: Some(settings),
        }
    }

    #[must_use]
    pub const fn settings(&self) -> Option<&PcrSettings> {
        self.settings.as_ref()
    }

    #[must_use]
    pub const fn is_active(&self) -> bool {
        self.settings.is_some()
    }
}

// ---------------------------------------------------------------------------
// Recuentos agregados (F4.9.9) y cota publicada (F4.9.10).
// ---------------------------------------------------------------------------

/// Recuentos agregados de la amplificación de una tanda.
///
/// **Nada por pareja.** Lo más fino que aparece es un tamaño de familia, que es
/// un número, no una identidad.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PcrReport {
    families: u64,
    pairs: u64,
    family_size_histogram: std::collections::BTreeMap<u32, u64>,
    families_truncated_at_end: u64,
    gc_attempts_total: u64,
    gc_attempts_exhausted_families: u64,
    /// Resbalones del linaje por locus: **todos** los ciclos con paso distinto
    /// de cero, de todas las copias de todas las familias.
    lineage_slips_by_locus: std::collections::BTreeMap<crate::identity::RepeatLocusId, u64>,
    /// De los anteriores, los que una copia **recibió de su progenitor**.
    ///
    /// Con familias de tamaño 1 vale exactamente cero, y esa es la prueba que
    /// distingue este recuento de uno que contara todos los resbalones.
    lineage_inherited_slips_by_locus:
        std::collections::BTreeMap<crate::identity::RepeatLocusId, u64>,
    /// Pasos netos de linaje que de verdad **llegaron a escribirse** en una
    /// molécula, por locus.
    ///
    /// # Por qué hace falta un tercer recuento y no basta con los dos de
    /// arriba
    ///
    /// Los dos de arriba cuentan lo que la genealogía **sorteó**. Entre el
    /// sorteo y la escritura hay dos filtros: un paso neto que no quepa en la
    /// molécula que se va a emitir —cosa que puede pasar con una quimera, que
    /// deja de ser la molécula sobre la que se sorteó— y un paso neto que salga
    /// cero por cancelación. Sin este recuento, un paso descartado al aplicarlo
    /// es **invisible**: los dos primeros seguirían diciendo que hubo
    /// resbalones y nadie podría saber cuántos llegaron a la lectura.
    lineage_applied_steps_by_locus: std::collections::BTreeMap<crate::identity::RepeatLocusId, u64>,
    /// Unidades de operación realmente gastadas por la amplificación, por
    /// encima del coste que F4 ya tenía. Existe para que una prueba pueda
    /// comparar el gasto real con la cota publicada; no se publica.
    operation_units: u64,
}

impl PcrReport {
    #[must_use]
    pub const fn families(&self) -> u64 {
        self.families
    }

    #[must_use]
    pub const fn pairs(&self) -> u64 {
        self.pairs
    }

    /// Parejas que son duplicado de otra de su familia: `N − familias`.
    #[must_use]
    pub const fn duplicate_pairs(&self) -> u64 {
        self.pairs.saturating_sub(self.families)
    }

    /// Fracción **realizada**, `1 − familias/parejas`. Es la que el manifiesto
    /// publica, y es **mayor** que la declarada siempre que haya desigualdad.
    #[must_use]
    pub fn realized_duplicate_fraction(&self) -> f64 {
        if self.pairs == 0 {
            return 0.0;
        }
        1.0 - (self.families as f64) / (self.pairs as f64)
    }

    #[must_use]
    pub const fn family_size_histogram(&self) -> &std::collections::BTreeMap<u32, u64> {
        &self.family_size_histogram
    }

    #[must_use]
    pub const fn families_truncated_at_end(&self) -> u64 {
        self.families_truncated_at_end
    }

    #[must_use]
    pub const fn gc_attempts_total(&self) -> u64 {
        self.gc_attempts_total
    }

    #[must_use]
    pub const fn gc_attempts_exhausted_families(&self) -> u64 {
        self.gc_attempts_exhausted_families
    }

    #[must_use]
    pub const fn lineage_slips_by_locus(
        &self,
    ) -> &std::collections::BTreeMap<crate::identity::RepeatLocusId, u64> {
        &self.lineage_slips_by_locus
    }

    #[must_use]
    pub const fn lineage_inherited_slips_by_locus(
        &self,
    ) -> &std::collections::BTreeMap<crate::identity::RepeatLocusId, u64> {
        &self.lineage_inherited_slips_by_locus
    }

    #[must_use]
    pub const fn lineage_applied_steps_by_locus(
        &self,
    ) -> &std::collections::BTreeMap<crate::identity::RepeatLocusId, u64> {
        &self.lineage_applied_steps_by_locus
    }

    #[must_use]
    pub const fn operation_units(&self) -> u64 {
        self.operation_units
    }
}

/// Acumulador de los recuentos de una biblioteca, fusionable entre bibliotecas.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PcrAccumulator(PcrReport);

impl PcrAccumulator {
    /// Registra una familia con el tamaño **realmente emitido**, que es el que
    /// hace que la suma del histograma ponderada por su índice valga
    /// exactamente el número de parejas.
    pub fn record_family(&mut self, emitted_pairs: u32, truncated: bool) {
        self.0.families += 1;
        self.0.pairs += u64::from(emitted_pairs);
        *self
            .0
            .family_size_histogram
            .entry(emitted_pairs)
            .or_insert(0) += 1;
        if truncated {
            self.0.families_truncated_at_end += 1;
        }
    }

    /// Registra los intentos de aceptación de una familia. `exhausted` indica
    /// que se agotó el tope y **se aceptó el último intento**.
    pub fn record_gc_attempts(&mut self, attempts: u64, exhausted: bool) {
        self.0.gc_attempts_total += attempts;
        if exhausted {
            self.0.gc_attempts_exhausted_families += 1;
        }
    }

    pub fn record_operation_units(&mut self, units: u64) {
        self.0.operation_units = self.0.operation_units.saturating_add(units);
    }

    /// Registra los resbalones de la genealogía de una familia.
    ///
    /// Los dos recuentos van **por separado** a propósito: el segundo es un
    /// subconjunto del primero, y publicarlos así es lo que permite ver si la
    /// herencia está ocurriendo de verdad. Un solo número no lo diría.
    pub fn record_lineage(&mut self, lineage: &crate::clonal_stutter::FamilyLineage) {
        for (locus_id, count) in lineage.slips_by_locus() {
            *self
                .0
                .lineage_slips_by_locus
                .entry(locus_id.clone())
                .or_insert(0) += count;
        }
        for (locus_id, count) in lineage.inherited_slips_by_locus() {
            *self
                .0
                .lineage_inherited_slips_by_locus
                .entry(locus_id.clone())
                .or_insert(0) += count;
        }
        // Los loci que nunca llegaron a escribir un paso aparecen con cero, no
        // faltan: un recuento ausente y un cero se leen igual de mal, y solo
        // uno de los dos dice que el locus estuvo y no llegó nada.
        for locus_id in lineage.slips_by_locus().keys() {
            self.0
                .lineage_applied_steps_by_locus
                .entry(locus_id.clone())
                .or_insert(0);
        }
    }

    /// Registra los pasos netos de linaje que **se escribieron** en una
    /// molécula. Los recibe ya filtrados por quien los aplica.
    pub fn record_applied_lineage(&mut self, applications: &[crate::stutter::StutterApplication]) {
        for application in applications {
            *self
                .0
                .lineage_applied_steps_by_locus
                .entry(application.locus_id().clone())
                .or_insert(0) += 1;
        }
    }

    pub fn merge(&mut self, report: &PcrReport) {
        self.0.families += report.families;
        self.0.pairs += report.pairs;
        for (size, count) in &report.family_size_histogram {
            *self.0.family_size_histogram.entry(*size).or_insert(0) += count;
        }
        self.0.families_truncated_at_end += report.families_truncated_at_end;
        self.0.gc_attempts_total += report.gc_attempts_total;
        self.0.gc_attempts_exhausted_families += report.gc_attempts_exhausted_families;
        for (locus_id, count) in &report.lineage_slips_by_locus {
            *self
                .0
                .lineage_slips_by_locus
                .entry(locus_id.clone())
                .or_insert(0) += count;
        }
        for (locus_id, count) in &report.lineage_inherited_slips_by_locus {
            *self
                .0
                .lineage_inherited_slips_by_locus
                .entry(locus_id.clone())
                .or_insert(0) += count;
        }
        for (locus_id, count) in &report.lineage_applied_steps_by_locus {
            *self
                .0
                .lineage_applied_steps_by_locus
                .entry(locus_id.clone())
                .or_insert(0) += count;
        }
        self.0.operation_units = self
            .0
            .operation_units
            .saturating_add(report.operation_units);
    }

    #[must_use]
    pub fn report(self) -> PcrReport {
        self.0
    }

    #[must_use]
    pub const fn view(&self) -> &PcrReport {
        &self.0
    }
}

/// Cota de recursos **publicada, no impuesta**.
///
/// # Por qué no se impone
///
/// La fase de secuenciación **nunca ha tenido** control de recursos —los
/// límites viven solo en la vía de haplotipos—, de modo que estrenar aquí un
/// rechazo por presupuesto sería añadir un contrato que no le corresponde
/// (decisión S-27). Se publica y se deja la imposición para quien tenga esa
/// responsabilidad.
///
/// # Lo que **no** acota, dicho expresamente
///
/// Los sorteos por rechazo que el recorrido ya usa —la normal truncada del
/// inserto y el índice uniforme— no tienen cota superior de palabras de
/// generador en el peor caso, y el proyecto ya lo reconoce. Este trabajo no
/// añade ningún sorteo por rechazo nuevo; sí multiplica por el tope de intentos
/// los sorteos de inserto de una familia.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PcrResourceBound {
    operation_units: u64,
    memory_bytes: u64,
    regeneration_passes: u32,
}

impl PcrResourceBound {
    #[must_use]
    pub const fn operation_units(self) -> u64 {
        self.operation_units
    }

    #[must_use]
    pub const fn memory_bytes(self) -> u64 {
        self.memory_bytes
    }

    /// El recorrido regenera cada biblioteca **dos veces**. La cota lo declara
    /// en vez de suponer una sola pasada.
    #[must_use]
    pub const fn regeneration_passes(self) -> u32 {
        self.regeneration_passes
    }

    /// Cota de la tanda a partir de las cotas de sus bibliotecas.
    ///
    /// Las operaciones se **suman**, porque las bibliotecas se recorren una
    /// tras otra; la memoria se toma como el **máximo**, porque el estado de
    /// familia de una biblioteca se libera antes de abrir la siguiente.
    #[must_use]
    pub const fn combine(self, other: Self) -> Self {
        Self {
            operation_units: self.operation_units.saturating_add(other.operation_units),
            memory_bytes: if self.memory_bytes > other.memory_bytes {
                self.memory_bytes
            } else {
                other.memory_bytes
            },
            regeneration_passes: PCR_REGENERATION_PASSES,
        }
    }

    /// Cota neutra para empezar a combinar.
    #[must_use]
    pub const fn empty() -> Self {
        Self {
            operation_units: 0,
            memory_bytes: 0,
            regeneration_passes: PCR_REGENERATION_PASSES,
        }
    }
}

/// Número de pasadas que el recorrido F4 hace sobre cada biblioteca.
pub const PCR_REGENERATION_PASSES: u32 = 2;

/// Calcula la cota **sin ejecutar nada**, desde la configuración y la
/// geometría conocidas antes de abrir el primer ámbito aleatorio.
///
/// ```text
/// ops ≤ N × [ Q + 2^c_early + A × (P + 3·L_max) + k_max + 2·c·T + c + 1 + k_max ]
///       + Q × c_early × 4^c_early
/// ```
///
/// donde `N` son las parejas de la biblioteca, `P` sus copias haplotípicas,
/// `L_max` el inserto máximo del perfil y `T` los tractos declarados que caben
/// en un inserto.
///
/// El término por pareja `2·c·T + c + 1` es el del **linaje**: dos decisiones
/// por ciclo y tracto —el resbalón y, como mucho, su paso—, más los `c + 1`
/// pesos de la divergencia y la elección del progenitor.
///
/// # Dos holguras declaradas, para que nadie las descubra por sorpresa
///
/// - El término de la programación dinámica queda **fuera** del producto por
///   `N` porque las tablas se construyen **una vez por clase de eficiencia** al
///   resolver el plan, no una vez por familia. El código hace exactamente eso
///   desde que se corrigió; antes no, y la cota mentía.
/// - El factor de composición se recalcula en **cada intento**, mientras la
///   fórmula lo carga una vez por pareja. Son a lo sumo `C_MAX` multiplicaciones
///   por intento, absorbidas de sobra por el término `A × (P + 3·L_max)` del
///   mismo intento, que vale al menos `A × 3 × L_max`. La cota sigue siendo
///   válida; queda dicho para que no parezca un descuido.
/// - Con un objetivo por **número de parejas** en vez de por cobertura no hay
///   plan de cobertura y la cota se calcularía con cero parejas. No es
///   alcanzable desde configuración —F4 exige cobertura— y queda anotado por si
///   alguna vez se reabre ese camino.
#[must_use]
pub fn resource_bound(
    settings: &PcrSettings,
    pair_count: u64,
    haplotype_copies: u64,
    maximum_insert_bases: u64,
    tracts_per_insert: u64,
) -> PcrResourceBound {
    let classes = settings.efficiency_classes().len() as u64;
    let early_states = 1_u64 << settings.early_cycles();
    let attempts = u64::from(settings.gc_max_attempts());
    let family_slots = u64::from(settings.max_family_size());
    let cycles = u64::from(settings.cycles());
    let per_pair = classes
        .saturating_add(early_states)
        .saturating_add(attempts.saturating_mul(
            haplotype_copies.saturating_add(3_u64.saturating_mul(maximum_insert_bases)),
        ))
        .saturating_add(family_slots)
        .saturating_add(
            2_u64
                .saturating_mul(cycles)
                .saturating_mul(tracts_per_insert),
        )
        .saturating_add(cycles)
        .saturating_add(1)
        .saturating_add(family_slots);
    let dynamic_programming = classes
        .saturating_mul(u64::from(settings.early_cycles()))
        .saturating_mul(early_states.saturating_mul(early_states));
    let operation_units = pair_count
        .saturating_mul(per_pair)
        .saturating_add(dynamic_programming);
    // Memoria añadida por encima del pico actual, independiente de `N` y del
    // genoma: la molécula, las tablas de pesos, los pasos netos de la familia
    // en curso y el histograma.
    let memory_bytes = maximum_insert_bases
        .saturating_add(
            8_u64.saturating_mul(
                classes
                    .saturating_mul(early_states)
                    .saturating_add(family_slots)
                    .saturating_add(cycles)
                    .saturating_add(1),
            ),
        )
        // Los pasos del linaje de la familia en curso. **Lleva un factor
        // `c + 1` que la primera redacción de esta cota no tenía**: para que
        // una copia herede el prefijo de su progenitor hacen falta los pasos
        // de este **ciclo a ciclo**, no solo su suma, y eso son `c` valores
        // por copia y tracto más el neto que sobrevive al cálculo. Declararlo
        // con un solo valor por copia y tracto habría dejado la cota por
        // debajo del gasto real.
        .saturating_add(
            u64::from(crate::technical_profile::PCR_MAXIMUM_FAMILY_SIZE)
                .saturating_mul(tracts_per_insert)
                .saturating_mul(cycles.saturating_add(1))
                .saturating_mul(8),
        )
        .saturating_add(8_u64.saturating_mul(family_slots.saturating_add(1)))
        // Los dos recuentos de linaje por locus, los de la familia en curso y
        // los del evento entero. Son **cuatro** mapas con a lo sumo `T`
        // entradas cada uno, y su clave es un `RepeatLocusId`, que comparte la
        // cadena con el locus declarado en vez de copiarla: veinticuatro bytes
        // por entrada y no la longitud del nombre que el autor del escenario
        // haya elegido. Se cuentan aquí para que la frase «el estado añadido no
        // depende de N ni del genoma» siga siendo cierta término a término y no
        // por absorción.
        .saturating_add(4_u64.saturating_mul(tracts_per_insert).saturating_mul(24));
    PcrResourceBound {
        operation_units,
        memory_bytes,
        regeneration_passes: PCR_REGENERATION_PASSES,
    }
}

// ---------------------------------------------------------------------------
// Documentos publicados.
// ---------------------------------------------------------------------------

/// Bloque de parámetros que publica la **configuración efectiva**.
///
/// Publica la fracción **declarada**, la `m` derivada con precisión de ida y
/// vuelta, el número de iteraciones y el nombre del método. La
/// reproducibilidad de una ejecución se ancla en la `m` publicada, **no** en
/// volver a derivarla: quien quiera repetir la ejecución en otra máquina lee
/// este número y no depende de que su exponencial coincida.
///
/// `expected_duplicate_fraction` es `f(m)` recalculada desde la `m` derivada, y
/// está aquí para que la ida y vuelta se vea sin tener que rehacerla.
#[must_use]
pub fn effective_settings_document(settings: &PcrSettings) -> serde_json::Value {
    let mut document = serde_json::json!({
        "applicable": true,
        "algorithm": PCR_AMPLIFICATION_ALGORITHM,
        "rng_catalog": settings.rng_catalog(),
        "cycles": settings.cycles(),
        "efficiency": settings.efficiency(),
        "early_cycles": settings.early_cycles(),
        "declared_duplicate_fraction": settings.declared_duplicate_fraction(),
        "expected_duplicate_fraction":
            expected_duplicate_fraction(settings.mean_pairs_per_molecule()),
        "mean_pairs_per_molecule": settings.mean_pairs_per_molecule(),
        "duplicate_fraction_inversion": {
            "method": settings.inversion().method(),
            "iterations": settings.inversion().iterations(),
        },
        "max_family_size": settings.max_family_size(),
        "gc_max_attempts": settings.gc_max_attempts(),
        "gc_relative_efficiency": settings
            .gc_curve()
            .knots()
            .iter()
            .map(|(gc_fraction, relative_efficiency)| serde_json::json!({
                "gc_fraction": gc_fraction,
                "relative_efficiency": relative_efficiency,
            }))
            .collect::<Vec<_>>(),
        "efficiency_classes": settings
            .efficiency_classes()
            .iter()
            .map(|class| serde_json::json!({
                "efficiency": class.efficiency(),
                "weight": class.weight(),
            }))
            .collect::<Vec<_>>(),
    });
    if let Some(clonal) = settings.clonal_stutter()
        && let Some(object) = document.as_object_mut()
    {
        object.insert(
            "stutter".to_owned(),
            crate::clonal_stutter::effective_settings_document(clonal),
        );
    }
    document
}

/// Bloque que publica el **manifiesto de evento**: los mismos parámetros más
/// los recuentos realizados y la cota calculada.
///
/// `unique_molecule_coverage_x` es la cobertura de **moléculas distintas**, que
/// es la que un analista compara tras eliminar duplicados: la cobertura básica
/// multiplicada por la proporción de familias sobre parejas. Se deriva de la
/// fracción realizada, no se mide aparte.
///
/// # Por qué el agregado de varias bibliotecas es exacto y no una media
///
/// La duda es legítima: la fórmula usa la proporción **agregada**
/// `familias_total / parejas_total`, y lo correcto sería sumar la de cada
/// biblioteca. Coinciden, y no por casualidad. La cobertura básica de una
/// biblioteca es `parejas_lib × 2 × longitud_de_lectura / denominador`, y el
/// **denominador es el mismo para todas**: hay una sola referencia por evento.
/// Llamando `c` a esa constante,
///
/// ```text
/// Σ_lib c·parejas_lib × (familias_lib/parejas_lib) = c · Σ_lib familias_lib
/// (Σ_lib c·parejas_lib) × (familias_total/parejas_total) = c · familias_total
/// ```
///
/// y las dos expresiones son la misma. Si algún día un evento admitiera
/// bibliotecas con longitudes de lectura distintas, dejarían de serlo.
#[must_use]
pub fn manifest_document(
    settings: &PcrSettings,
    report: &PcrReport,
    bound: PcrResourceBound,
    basic_wgs_coverage_x: f64,
) -> serde_json::Value {
    let mut document = effective_settings_document(settings);
    let unique_molecule_coverage_x = if report.pairs() == 0 {
        0.0
    } else {
        basic_wgs_coverage_x * (report.families() as f64) / (report.pairs() as f64)
    };
    if let Some(object) = document.as_object_mut() {
        object.insert("families".to_owned(), serde_json::json!(report.families()));
        object.insert("pairs".to_owned(), serde_json::json!(report.pairs()));
        object.insert(
            "duplicate_pairs".to_owned(),
            serde_json::json!(report.duplicate_pairs()),
        );
        object.insert(
            "realized_duplicate_fraction".to_owned(),
            serde_json::json!(report.realized_duplicate_fraction()),
        );
        object.insert(
            "family_size_histogram".to_owned(),
            serde_json::Value::Array(
                report
                    .family_size_histogram()
                    .iter()
                    .map(|(size, families)| {
                        serde_json::json!({"family_size": size, "families": families})
                    })
                    .collect(),
            ),
        );
        object.insert(
            "families_truncated_at_end".to_owned(),
            serde_json::json!(report.families_truncated_at_end()),
        );
        object.insert(
            "gc_attempts_total".to_owned(),
            serde_json::json!(report.gc_attempts_total()),
        );
        object.insert(
            "gc_attempts_exhausted_families".to_owned(),
            serde_json::json!(report.gc_attempts_exhausted_families()),
        );
        object.insert(
            "unique_molecule_coverage_x".to_owned(),
            serde_json::json!(unique_molecule_coverage_x),
        );
        // Los dos recuentos del linaje solo aparecen cuando hay linaje. Una
        // biblioteca amplificada sin loci de repetición publica el mismo
        // documento que antes de F4.10, byte a byte.
        if settings.clonal_stutter().is_some() {
            object.insert(
                "lineage_slips_by_locus".to_owned(),
                locus_counts_document(report.lineage_slips_by_locus(), "slips"),
            );
            object.insert(
                "lineage_inherited_slips_by_locus".to_owned(),
                locus_counts_document(report.lineage_inherited_slips_by_locus(), "slips"),
            );
            object.insert(
                "lineage_applied_steps_by_locus".to_owned(),
                locus_counts_document(report.lineage_applied_steps_by_locus(), "steps"),
            );
        }
        object.insert(
            "resource_bound".to_owned(),
            serde_json::json!({
                "operation_units": bound.operation_units(),
                "memory_bytes": bound.memory_bytes(),
                "regeneration_passes": bound.regeneration_passes(),
                "measured": false,
            }),
        );
    }
    document
}

/// Un recuento por locus, como lista ordenada de objetos con su nombre.
///
/// Se publica como lista y no como objeto con el locus de clave porque un
/// identificador de locus lo elige el autor del escenario, y un esquema con
/// claves libres no podría cerrarse con `additionalProperties: false`.
fn locus_counts_document(
    counts: &std::collections::BTreeMap<crate::identity::RepeatLocusId, u64>,
    field: &str,
) -> serde_json::Value {
    serde_json::Value::Array(
        counts
            .iter()
            .map(|(locus_id, count)| {
                serde_json::json!({"locus_id": locus_id.as_str(), field: count})
            })
            .collect(),
    )
}
