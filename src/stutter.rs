//! Tartamudeo de base de la polimerasa sobre un microsatélite (F4.8).
//!
//! # Qué es, en una frase
//!
//! Una molécula que contiene un tracto declarado puede **mostrar** un número de
//! unidades distinto del que el individuo realmente tiene, porque la polimerasa
//! resbala mientras la copia. Es un artefacto de la lectura, no un cambio
//! biológico: **la verdad del individuo no cambia**. El tartamudeo vive
//! exclusivamente en las moléculas que se secuencian, y por eso este módulo
//! solo toca [`PairedEndPair`] y jamás una copia haplotípica.
//!
//! # Por qué cabe en una fase que declara «sin PCR»
//!
//! Una biblioteca sin PCR evita la amplificación **en la preparación**, pero no
//! evita toda polimerasa: la formación de agrupaciones en la celda de flujo es
//! amplificación en puente y la secuenciación por síntesis es copia con
//! polimerasa. El deslizamiento ocurre en las dos. Por eso este artefacto **no
//! es una preparación amplificada**, **no introduce duplicados** y el perfil se
//! sigue registrando como `pcr_free` (decisión S-12).
//!
//! # El predictor es la tirada más larga, no el número total de unidades
//!
//! En genética forense está establecido que la **tirada ininterrumpida más
//! larga** predice el tartamudeo mejor que la longitud total del alelo, y es
//! coherente con que las interrupciones **estabilicen** el tracto. El valor lo
//! calcula la observación de F6.4.1 ([`observe_tract`]) y aquí se **consume**;
//! no se recalcula con otra convención.
//!
//! # La forma lineal no es una aproximación tosca, aunque lo parezca
//!
//! La literatura describe la fracción observada de tartamudeo de extremo a
//! extremo con una curva que **satura**. Quien conozca esa curva leerá la
//! dependencia lineal de este módulo como una simplificación grosera, y no lo
//! es: son **dos magnitudes distintas**. La dependencia lineal es la de la
//! probabilidad de resbalón **por ciclo**; la curva con meseta sale de
//! **acumular** esa probabilidad por ciclo sobre muchos ciclos. Las dos
//! afirmaciones son ciertas a la vez y no se contradicen.
//!
//! # La longitud del tracto modula, no suma
//!
//! El efecto más fuerte que hay sobre el tartamudeo es la longitud del tracto:
//! a lo largo del tramo medido de los mononucleótidos la fracción de lecturas
//! erróneas sube cerca de **19 000 veces**. Una recta no puede con eso, y por
//! eso esta dependencia es una **tabla interpolada en escala logarítmica** y
//! no un coeficiente.
//!
//! Lo delicado no es la tabla sino **dónde se ancla**. Las tasas por longitud
//! de motivo que este módulo ya usaba son promedios sobre los tractos de
//! aquella muestra: **ya llevan dentro** el efecto de la longitud. Si encima
//! se multiplicaran por la curva se contaría dos veces el mismo fenómeno. Por
//! eso la curva está **normalizada**: a la longitud de referencia el factor
//! vale exactamente `1.0` y la tasa medida sale intacta. La referencia es el
//! **tope del tramo medido**, de modo que la modulación solo actúa hacia
//! abajo y ningún tracto tartamudea más de lo que la tabla medida autoriza.
//!
//! Tres cosas que hay que saber antes de usarla, y que están dichas también en
//! el documento que publica cada ejecución:
//!
//! - **Solo hay curva medida para mononucleótidos y dinucleótidos.** Los
//!   motivos de tres a seis pares de base **no tienen dependencia ninguna**
//!   con la longitud del tracto. No se extrapola: una extrapolación aquí sería
//!   una cifra inventada con apariencia de medida.
//! - **Fuera del tramo medido se mantiene el valor del extremo.** Como el
//!   ancla es el extremo superior, un tracto de 30 unidades **no tartamudea
//!   más** que uno de 15. Es la consecuencia honesta de no extrapolar.
//! - **La sustitución del autor del escenario prevalece y no se modula.**
//!   Quien declara la tasa de un locus la ha elegido para ese locus, con su
//!   longitud dentro.
//!
//! # Por qué DTC-067 no aplica aquí, dicho por escrito
//!
//! DTC-067 obliga a que los **recuentos** aleatorios salgan de
//! `poisson-count-v1`. Aquí no hay ningún recuento de sucesos sobre
//! oportunidades: hay un **ensayo binario** —¿resbala o no?— y una **elección
//! entre alternativas** —¿cuánto y en qué dirección?—. Los consumidores
//! normativos que corresponden ya existen: `bernoulli-v1` (DTC-088) y
//! `categorical-f64-v1` (DTC-072). Decirlo evita que alguien lea la regla al
//! revés y bloquee un uso que nunca estuvo prohibido.
//!
//! # Sin reintentos
//!
//! El proyecto rechaza los reintentos en la selección aleatoria: los candidatos
//! inviables se excluyen **antes** de sortear, nunca se sortea y se descarta.
//! Aquí eso significa que un paso solo entra en el sorteo categórico si deja al
//! menos una unidad completa en la tirada y deja la molécula al menos tan larga
//! como una lectura. Si solo queda un candidato, `categorical-f64-v1` lo
//! devuelve **sin consumir** ni una palabra, y si no queda ninguno no se llega
//! ni al ensayo de Bernoulli.

use std::collections::BTreeMap;
use std::sync::Arc;

use crate::coordinates::GenomicInterval;
use crate::haplotype::{HaplotypeCopy, Individual};
use crate::identity::{EventId, RepeatLocusId, SequenceId};
use crate::paired_end::{PairedEndFragment, PairedEndPair, emit_paired_end};
use crate::randomness::{
    F4TechnicalDecision, RandomSamplingError, RootSeed, categorical_f64_with_source,
    derive_f4_technical_subseed, f4_bernoulli,
};
use crate::segments::CurrentIndelPiece;
use crate::tandem_repeat::{RepeatMotif, observe_tract};

/// Versión del algoritmo, publicada en la configuración efectiva y en el
/// manifiesto. Cambiar la forma de la probabilidad o la estructura del paso
/// obliga a estrenar versión.
pub const STUTTER_ALGORITHM: &str = "stutter-v2";

/// Longitud de motivo mayor que la calibración cubre, en pares de base.
pub const STUTTER_MOTIF_LENGTH_MAX: usize = 6;

/// **Los cuatro pasos admitidos, en orden canónico.**
///
/// La estructura **no se configura**: la fija el algoritmo. Contracciones de
/// hasta tres unidades, expansión de una sola. Lo que sí se configura son los
/// pesos relativos dentro de esta estructura, que es donde se expresa el sesgo
/// hacia la contracción. Una configuración que declare `+2` se rechaza.
pub const STUTTER_ADMITTED_STEPS: [i64; 4] = [-3, -2, -1, 1];

/// Fracción de lecturas afectadas por longitud de motivo, **columna sin PCR**.
///
/// El índice es la longitud del motivo menos uno: la posición 0 es el
/// mononucleótido. Con el coeficiente de la tirada en cero estas cifras **son**
/// directamente la probabilidad de resbalón, sin más aritmética.
///
/// El último valor, el del motivo de 6 pb, es **el menos fiable**: el propio
/// autor advierte de que procede del ruido de analizar pocos loci.
pub const STUTTER_DEFAULT_MOTIF_LENGTH_RATES: [f64; STUTTER_MOTIF_LENGTH_MAX] =
    [0.027, 0.011, 0.0029, 0.0016, 0.0015, 0.0010];

/// Fracción de lecturas afectadas por longitud de motivo, **columna con PCR**.
///
/// Mismo índice y misma fuente que la columna sin PCR de arriba: es la otra
/// mitad de la misma tabla. Vive aquí, junto a su hermana, porque separarlas
/// invitaría a componerlas mal.
///
/// # Las dos columnas son **totales**, no componentes
///
/// Cada conjunto pasó por una celda de flujo con su amplificación en puente y
/// por la síntesis de la secuenciación, de modo que la columna sin PCR **ya
/// contiene** el mosaicismo somático, el deslizamiento del puente y el de la
/// síntesis, y esta contiene todo eso **más** la PCR de biblioteca. Por eso lo
/// correcto es **restar** para despejar el incremento propio de la
/// amplificación, nunca sumar: sumar contaría dos veces el mismo fenómeno.
/// Quien la use es [`crate::clonal_stutter`].
pub const STUTTER_WITH_PCR_MOTIF_LENGTH_RATES: [f64; STUTTER_MOTIF_LENGTH_MAX] =
    [0.17, 0.038, 0.011, 0.0069, 0.0074, 0.0140];

/// Procedencia de esas cifras. Viaja con los valores a la configuración
/// efectiva y al manifiesto: unos números embebidos sin su fuente no se
/// distinguen de unos inventados.
pub const STUTTER_CALIBRATION_SOURCE: &str = "Rates: Gymrek (2016), bioRxiv 043448, table 1, CC BY 4.0; model form: Raz et al. (2019), Nucleic Acids Research 47(5):2436-2445, CC BY 4.0";
/// Especie de la calibración.
pub const STUTTER_CALIBRATION_SPECIES: &str = "human";
/// Tecnología de la calibración.
pub const STUTTER_CALIBRATION_TECHNOLOGY: &str = "illumina";
/// Columna concreta de la que salen las seis cifras.
pub const STUTTER_CALIBRATION_COLUMN: &str = "pcr_free";

/// **Coeficiente de la tirada: arranca en cero, a propósito.**
///
/// La curva que lo calibraría está publicada pero **sin digitalizar**, y el
/// responsable no ha decidido hacerlo (pendiente P-1). Con el coeficiente en
/// cero la probabilidad se reduce a la fracción medida por longitud de motivo,
/// que **sí** está calibrada. Un valor por omisión distinto de cero afirmaría
/// una calibración que nadie ha hecho, y por eso el manifiesto lo publica
/// marcado como no calibrado en vez de callarlo.
pub const STUTTER_DEFAULT_RUN_LENGTH_COEFFICIENT: f64 = 0.0;

/// Lo que el manifiesto publica sobre la calibración del coeficiente.
pub const STUTTER_RUN_LENGTH_COEFFICIENT_CALIBRATED: bool = false;

/// El desajuste conocido entre la fuente del modelo y la calibración adoptada.
///
/// La fuente distingue el **tipo** de motivo; la calibración del proyecto va
/// por **longitud** de motivo. No es un error: es una diferencia que el lector
/// del manifiesto tiene que conocer para no atribuir a la fuente una precisión
/// que aquí no se conserva.
pub const STUTTER_KNOWN_CALIBRATION_MISMATCH: &str =
    "the source distinguishes motif type while this calibration is indexed by motif length";

// ---------------------------------------------------------------------------
// Dependencia con la longitud del tracto
// ---------------------------------------------------------------------------

/// **De dónde sale la curva de la longitud del tracto.**
pub const STUTTER_TRACT_LENGTH_SOURCE: &str =
    "Fungtammasan et al. (2015), Genome Research 25(5):736-749, figure 1A, CC BY-NC";

/// **La dependencia con la longitud del tracto es TOSCA**, y se publica así.
pub const STUTTER_TRACT_LENGTH_COARSE: bool = true;

/// Por qué es tosca, dicho entero para que nadie tenga que ir a buscarlo.
///
/// Son cuatro cosas distintas y las cuatro importan: cómo se obtuvo la curva,
/// dónde se ancló, qué se hace fuera del tramo medido y para qué motivos no
/// hay curva ninguna.
pub const STUTTER_TRACT_LENGTH_COARSE_REASON: &str = "the tract-length curve was read by eye from a figure plotted in logarithmic scale, with an \
     uncertainty of about +/-0.15 log10 units (about +/-40 percent in linear scale); it is \
     normalised at the top of the measured range, so it only modulates the measured motif-length \
     rate downwards and never above it; outside the measured range the end value is held and \
     nothing is extrapolated; and only mononucleotides and dinucleotides have a measured curve, so \
     motifs of three to six bases carry no tract-length dependence at all";

/// **La longitud de tracto de referencia, en pares de base, por longitud de
/// motivo.** `None` significa que esa longitud de motivo no tiene curva medida
/// y por tanto **no tiene dependencia ninguna** con la longitud del tracto.
///
/// # De dónde sale, porque es la decisión entera
///
/// Es el **tope del tramo medido** de cada curva: 15 pares de base en
/// mononucleótidos y 34 en dinucleótidos. El responsable lo adoptó el 27 de
/// septiembre de 2026, y las tres razones son estas:
///
/// 1. Es la única elección que **no satura en ningún punto**. La curva sube
///    cerca de 19 000 veces a lo largo del tramo medido de los
///    mononucleótidos, y la tasa medida con PCR ya vale 0,17: cualquier ancla
///    por debajo de unos 11,7 pares de base llevaría la probabilidad a 1 —
///    resbalón seguro en todas las moléculas— dentro del propio tramo medido.
/// 2. Es el punto que el documento de calibración **ya usó** para contrastar
///    la curva con la tabla de tasas por longitud de motivo.
/// 3. Deja la modulación actuando **solo hacia abajo**: ningún tracto
///    tartamudea más de lo que la tabla medida ya autoriza.
///
/// # Su defecto, dicho con la misma claridad
///
/// **Un tracto de 30 unidades no tartamudea más que uno de 15.** Por encima
/// del tramo medido se mantiene el valor del extremo, que es el ancla, de modo
/// que la modulación vale exactamente 1 para todo lo que sea más largo. Es la
/// consecuencia honesta de no extrapolar, no un error del modelo, y el
/// responsable la aceptó sabiéndolo.
pub const STUTTER_TRACT_LENGTH_REFERENCE_BASES: [Option<u64>; STUTTER_MOTIF_LENGTH_MAX] =
    [Some(15), Some(34), None, None, None, None];

/// Primera entrada tabulada de las dos curvas, en **unidades**.
///
/// Las dos empiezan en cinco unidades porque ahí empieza el tramo medido:
/// 5 pares de base en mononucleótidos y 10 en dinucleótidos.
const TRACT_FACTOR_FIRST_UNITS: u64 = 5;

/// Factores de los **mononucleótidos**, de 5 a 15 unidades.
///
/// Uno por unidad, ya normalizados: el último vale exactamente `1.0` porque es
/// el ancla. Salen de interpolar en **escala logarítmica** la media de las dos
/// columnas digitalizadas y restarles el valor del ancla. La prueba
/// `the_tract_length_table_is_the_published_curve_normalised_at_the_anchor`
/// los rederiva desde los valores publicados sin mirar esta tabla.
const MONONUCLEOTIDE_TRACT_FACTORS: [f64; 11] = [
    5.30884444230989e-05,  // 5 unidades = 5 pb
    0.0002660725059798814, // 6 unidades = 6 pb
    0.001333521432163324,  // 7 unidades = 7 pb
    0.0047315125896148025, // 8 unidades = 8 pb
    0.016788040181225598,  // 9 unidades = 9 pb
    0.04466835921509633,   // 10 unidades = 10 pb
    0.11885022274370183,   // 11 unidades = 11 pb
    0.2238721138568339,    // 12 unidades = 12 pb
    0.4216965034285822,    // 13 unidades = 13 pb
    0.6493816315762113,    // 14 unidades = 14 pb
    1.0,                   // 15 unidades = 15 pb, el ancla
];

/// Factores de los **dinucleótidos**, de 5 a 17 unidades, es decir de 10 a 34
/// pares de base. Mismo origen y misma normalización que los de arriba.
const DINUCLEOTIDE_TRACT_FACTORS: [f64; 13] = [
    0.003758374042884443, // 5 unidades = 10 pb
    0.008175230379436493, // 6 unidades = 12 pb
    0.01778279410038922,  // 7 unidades = 14 pb
    0.03868120546330523,  // 8 unidades = 16 pb
    0.08413951416451952,  // 9 unidades = 18 pb
    0.12409377607517195,  // 10 unidades = 20 pb
    0.18302061063110558,  // 11 unidades = 22 pb
    0.26992928231542934,  // 12 unidades = 24 pb
    0.39810717055349715,  // 13 unidades = 26 pb
    0.5011872336272722,   // 14 unidades = 28 pb
    0.6309573444801934,   // 15 unidades = 30 pb
    0.7943282347242814,   // 16 unidades = 32 pb
    1.0,                  // 17 unidades = 34 pb, el ancla
];

/// La tabla de una longitud de motivo, o `None` si no hay curva medida.
#[must_use]
const fn tract_factors(motif_length: u64) -> Option<&'static [f64]> {
    match motif_length {
        1 => Some(&MONONUCLEOTIDE_TRACT_FACTORS),
        2 => Some(&DINUCLEOTIDE_TRACT_FACTORS),
        _ => None,
    }
}

/// **Cuánto multiplica la longitud del tracto a la tasa por longitud de
/// motivo.**
///
/// Devuelve exactamente `1.0` —ninguna dependencia— cuando esa longitud de
/// motivo no tiene curva medida, que es el caso de los motivos de tres a seis
/// pares de base. **No se extrapola**: por debajo y por encima del tramo
/// medido se mantiene el valor del extremo, y como el ancla es el extremo
/// superior, el factor nunca pasa de `1.0`.
///
/// Es una función pura y no consume aleatoriedad. Vive aparte para que una
/// prueba pueda comparar dos tractos sin montar una ejecución.
#[must_use]
pub fn tract_length_factor(motif_length: u64, run_units: u64) -> f64 {
    let Some(factors) = tract_factors(motif_length) else {
        return 1.0;
    };
    let index = run_units.saturating_sub(TRACT_FACTOR_FIRST_UNITS);
    let last = factors.len() - 1;
    let index = usize::try_from(index).unwrap_or(last).min(last);
    factors[index]
}

/// El tramo medido de una longitud de motivo, **en unidades**, o `None` si no
/// hay curva. Lo publica el documento efectivo para que quien lea un número
/// sepa dónde deja de estar medido.
#[must_use]
pub fn tract_length_measured_units(motif_length: u64) -> Option<(u64, u64)> {
    tract_factors(motif_length).map(|factors| {
        (
            TRACT_FACTOR_FIRST_UNITS,
            TRACT_FACTOR_FIRST_UNITS + factors.len() as u64 - 1,
        )
    })
}

/// **Probabilidad de resbalón por ciclo**, acotada a `[0,1]`.
///
/// `base_rate` es la tasa **ya modulada por la longitud del tracto** que
/// devuelve [`StutterSettings::modulated_rate_for`]; `coefficient` es el
/// coeficiente de la tirada (cero mientras no esté calibrado) y
/// `longest_run_units` es el predictor.
///
/// **La modulación no se hace aquí a propósito.** Se hace un paso antes,
/// porque depende de si el autor del escenario sustituyó la tasa de ese locus
/// —en cuyo caso su sustitución prevalece sin modular—, y esa es una decisión
/// de los parámetros, no de la aritmética.
///
/// Es una función pura y no consume aleatoriedad: existe aparte para que una
/// prueba pueda comparar dos tractos sin montar una ejecución.
#[must_use]
pub fn slip_probability(base_rate: f64, coefficient: f64, longest_run_units: u64) -> f64 {
    let value = coefficient.mul_add(longest_run_units as f64, base_rate);
    if !value.is_finite() {
        return 0.0;
    }
    value.clamp(0.0, 1.0)
}

/// Pesos relativos de los cuatro pasos admitidos, en el orden de
/// [`STUTTER_ADMITTED_STEPS`].
///
/// **Son obligatorios y no tienen valor por omisión** cuando el tartamudeo está
/// activo. Se buscó una fuente adoptable y no la hay: juntar cifras de trabajos
/// distintos sería síntesis propia disfrazada de fuente (pendiente P-3). La
/// documentación trae un ejemplo trabajado, marcado como ilustrativo, y ese
/// ejemplo **no** es un valor por omisión.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StutterStepWeights {
    weights: [f64; STUTTER_ADMITTED_STEPS.len()],
}

impl StutterStepWeights {
    /// Construye los pesos ya validados: finitos, no negativos y con al menos
    /// uno positivo. Sin ese último, el sorteo categórico no tendría ninguna
    /// alternativa que elegir.
    pub fn new(weights: [f64; STUTTER_ADMITTED_STEPS.len()]) -> Result<Self, StutterError> {
        for (index, weight) in weights.iter().enumerate() {
            if !weight.is_finite() || *weight < 0.0 {
                return Err(StutterError::InvalidStepWeight {
                    step: STUTTER_ADMITTED_STEPS[index],
                    weight: *weight,
                });
            }
        }
        if !weights.iter().any(|weight| *weight > 0.0) {
            return Err(StutterError::NoPositiveStepWeight);
        }
        Ok(Self { weights })
    }

    #[must_use]
    pub const fn weights(&self) -> &[f64; STUTTER_ADMITTED_STEPS.len()] {
        &self.weights
    }

    /// Peso del paso pedido, o `None` si ese paso no pertenece a la estructura.
    #[must_use]
    pub fn weight_of(&self, step: i64) -> Option<f64> {
        STUTTER_ADMITTED_STEPS
            .iter()
            .position(|admitted| *admitted == step)
            .map(|index| self.weights[index])
    }
}

/// Tasas efectivas por longitud de motivo, de 1 a [`STUTTER_MOTIF_LENGTH_MAX`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StutterMotifLengthRates {
    rates: [f64; STUTTER_MOTIF_LENGTH_MAX],
}

impl Default for StutterMotifLengthRates {
    fn default() -> Self {
        Self {
            rates: STUTTER_DEFAULT_MOTIF_LENGTH_RATES,
        }
    }
}

impl StutterMotifLengthRates {
    pub fn new(rates: [f64; STUTTER_MOTIF_LENGTH_MAX]) -> Result<Self, StutterError> {
        for (index, rate) in rates.iter().enumerate() {
            if !rate.is_finite() || !(0.0..=1.0).contains(rate) {
                return Err(StutterError::InvalidRate {
                    motif_length: index as u64 + 1,
                    rate: *rate,
                });
            }
        }
        Ok(Self { rates })
    }

    #[must_use]
    pub const fn rates(&self) -> &[f64; STUTTER_MOTIF_LENGTH_MAX] {
        &self.rates
    }

    /// Tasa de una longitud de motivo concreta, en pares de base.
    ///
    /// Devuelve `None` fuera del tramo calibrado. No se extrapola: una
    /// extrapolación aquí sería una cifra inventada con apariencia de medida.
    #[must_use]
    pub fn rate_for_motif_length(&self, motif_length: u64) -> Option<f64> {
        usize::try_from(motif_length)
            .ok()
            .filter(|length| (1..=STUTTER_MOTIF_LENGTH_MAX).contains(length))
            .map(|length| self.rates[length - 1])
    }
}

/// Parámetros efectivos del tartamudeo de un evento técnico.
#[derive(Clone, Debug, PartialEq)]
pub struct StutterSettings {
    enabled: bool,
    motif_length_rates: StutterMotifLengthRates,
    run_length_coefficient: f64,
    step_weights: Option<StutterStepWeights>,
    locus_rates: BTreeMap<String, f64>,
}

impl Default for StutterSettings {
    /// Desactivado. Es lo que hace que una ejecución anterior a este trabajo
    /// publique exactamente los mismos bytes: con el tartamudeo apagado no se
    /// abre la etapa `stutter` y no se consume ni una palabra aleatoria.
    fn default() -> Self {
        Self {
            enabled: false,
            motif_length_rates: StutterMotifLengthRates::default(),
            run_length_coefficient: STUTTER_DEFAULT_RUN_LENGTH_COEFFICIENT,
            step_weights: None,
            locus_rates: BTreeMap::new(),
        }
    }
}

impl StutterSettings {
    /// Construye los parámetros ya resueltos.
    ///
    /// Con el tartamudeo activo los pesos son obligatorios: sin ellos la
    /// configuración se rechaza nombrando el parámetro, en vez de repartir por
    /// omisión un sesgo que nadie ha declarado.
    pub fn new(
        enabled: bool,
        motif_length_rates: StutterMotifLengthRates,
        run_length_coefficient: f64,
        step_weights: Option<StutterStepWeights>,
        locus_rates: BTreeMap<String, f64>,
    ) -> Result<Self, StutterError> {
        if !run_length_coefficient.is_finite() || run_length_coefficient < 0.0 {
            return Err(StutterError::InvalidRunLengthCoefficient {
                coefficient: run_length_coefficient,
            });
        }
        if enabled && step_weights.is_none() {
            return Err(StutterError::MissingStepWeights);
        }
        for (locus_id, rate) in &locus_rates {
            if !rate.is_finite() || !(0.0..=1.0).contains(rate) {
                return Err(StutterError::InvalidLocusRate {
                    locus_id: locus_id.clone(),
                    rate: *rate,
                });
            }
        }
        Ok(Self {
            enabled,
            motif_length_rates,
            run_length_coefficient,
            step_weights,
            locus_rates,
        })
    }

    #[must_use]
    pub const fn enabled(&self) -> bool {
        self.enabled
    }

    #[must_use]
    pub const fn motif_length_rates(&self) -> &StutterMotifLengthRates {
        &self.motif_length_rates
    }

    #[must_use]
    pub const fn run_length_coefficient(&self) -> f64 {
        self.run_length_coefficient
    }

    #[must_use]
    pub const fn step_weights(&self) -> Option<&StutterStepWeights> {
        self.step_weights.as_ref()
    }

    /// Sustituciones de la tasa declaradas locus a locus.
    #[must_use]
    pub const fn locus_rates(&self) -> &BTreeMap<String, f64> {
        &self.locus_rates
    }

    /// Tasa base efectiva de un locus para una longitud de motivo dada: la
    /// sustitución del locus si la declara, y si no la de la calibración.
    ///
    /// **Sin modular por la longitud del tracto.** Quien quiera la tasa que el
    /// modelo usa de verdad tiene [`Self::modulated_rate_for`]; esta existe
    /// para publicar la tabla y para que la modulación se vea aplicarse en un
    /// solo sitio.
    #[must_use]
    pub fn base_rate_for(&self, locus_id: &str, motif_length: u64) -> Option<f64> {
        self.locus_rates.get(locus_id).copied().map_or_else(
            || self.motif_length_rates.rate_for_motif_length(motif_length),
            Some,
        )
    }

    /// **La tasa que el modelo usa**: la de la calibración modulada por la
    /// longitud del tracto, o la del autor del escenario tal cual.
    ///
    /// # La sustitución del autor prevalece, y no se modula
    ///
    /// Es la misma regla que la probabilidad por ciclo del tartamudeo clonal,
    /// y por el mismo motivo. Quien declara la tasa de un locus concreto la
    /// ha medido —o elegido— **para ese locus**, de modo que la longitud de
    /// su tracto ya está dentro. Multiplicarla encima por la curva contaría
    /// dos veces el mismo efecto, que es exactamente el error contra el que
    /// existe esta función.
    ///
    /// # Por qué modular y no sumar
    ///
    /// Las tasas por longitud de motivo de la calibración son **promedios
    /// sobre los tractos de aquella muestra**: ya llevan dentro el efecto de
    /// la longitud del tracto, promediado. La curva solo puede repartir ese
    /// promedio entre tractos largos y cortos, nunca añadirse a él. Por eso
    /// está normalizada: a la longitud de referencia el factor vale
    /// exactamente `1.0` y la tasa sale intacta.
    #[must_use]
    pub fn modulated_rate_for(
        &self,
        locus_id: &str,
        motif_length: u64,
        run_units: u64,
    ) -> Option<f64> {
        if let Some(rate) = self.locus_rates.get(locus_id) {
            return Some(*rate);
        }
        let base = self
            .motif_length_rates
            .rate_for_motif_length(motif_length)?;
        Some(base * tract_length_factor(motif_length, run_units))
    }
}

/// **Lo que se publica de la dependencia con la longitud del tracto.**
///
/// Sale de aquí y no de dos sitios para que la configuración efectiva y el
/// manifiesto no puedan discrepar. Lleva la marca de tosquedad con su motivo
/// entero, la fuente, y **por cada longitud de motivo** su longitud de
/// referencia y su tramo medido — o la declaración explícita de que no tiene
/// curva y por tanto no tiene dependencia ninguna.
#[must_use]
pub fn tract_length_dependence_document() -> serde_json::Value {
    let by_motif_length = (1..=STUTTER_MOTIF_LENGTH_MAX)
        .map(|motif_length| {
            let motif_length = motif_length as u64;
            let measured = tract_length_measured_units(motif_length);
            // **La tabla entera se publica**, no solo sus extremos. Es lo que
            // permite rehacer la modulación sin este programa, exactamente
            // igual que la derivación del tartamudeo clonal publica sus dos
            // fracciones de entrada y no solo su resultado. Un validador
            // independiente que tuviera que llevar su propia copia de estas
            // cifras dejaría de ser independiente el día que una de las dos
            // se moviera.
            let factors = measured.map_or_else(Vec::new, |(first, last)| {
                (first..=last)
                    .map(|units| {
                        serde_json::json!({
                            "units": units,
                            "factor": tract_length_factor(motif_length, units),
                        })
                    })
                    .collect()
            });
            serde_json::json!({
                "motif_length_bases": motif_length,
                "measured": measured.is_some(),
                "reference_tract_bases":
                    STUTTER_TRACT_LENGTH_REFERENCE_BASES[motif_length as usize - 1],
                "measured_units_min": measured.map(|(first, _)| first),
                "measured_units_max": measured.map(|(_, last)| last),
                "factors": factors,
            })
        })
        .collect::<Vec<_>>();
    serde_json::json!({
        "source": STUTTER_TRACT_LENGTH_SOURCE,
        "coarse": STUTTER_TRACT_LENGTH_COARSE,
        "coarse_reason": STUTTER_TRACT_LENGTH_COARSE_REASON,
        "extrapolated": false,
        "by_motif_length": by_motif_length,
    })
}

/// Documento común de los parámetros efectivos del tartamudeo.
///
/// Lo publican **la configuración efectiva y el manifiesto**, y sale de aquí
/// para que no puedan discrepar: dos constructores separados divergirían el día
/// que alguien añadiera un parámetro a uno solo.
#[must_use]
pub fn effective_settings_document(settings: &StutterSettings) -> serde_json::Value {
    let motif_length_rates = settings
        .motif_length_rates()
        .rates()
        .iter()
        .enumerate()
        .map(|(index, rate)| serde_json::json!({"motif_length_bases": index + 1, "rate": rate}))
        .collect::<Vec<_>>();
    let step_weights = STUTTER_ADMITTED_STEPS
        .iter()
        .map(|step| {
            serde_json::json!({
                "units": step,
                "weight": settings
                    .step_weights()
                    .and_then(|weights| weights.weight_of(*step))
                    .unwrap_or(0.0),
            })
        })
        .collect::<Vec<_>>();
    let locus_rates = settings
        .locus_rates()
        .iter()
        .map(|(locus_id, rate)| serde_json::json!({"locus_id": locus_id, "rate": rate}))
        .collect::<Vec<_>>();
    serde_json::json!({
        "enabled": settings.enabled(),
        "algorithm": STUTTER_ALGORITHM,
        "rng_catalog": crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4,
        "calibration": {
            "source": STUTTER_CALIBRATION_SOURCE,
            "species": STUTTER_CALIBRATION_SPECIES,
            "technology": STUTTER_CALIBRATION_TECHNOLOGY,
            "column": STUTTER_CALIBRATION_COLUMN,
            "known_mismatch": STUTTER_KNOWN_CALIBRATION_MISMATCH,
        },
        "motif_length_rates": motif_length_rates,
        "tract_length_dependence": tract_length_dependence_document(),
        "run_length_coefficient": {
            "value": settings.run_length_coefficient(),
            "calibrated": STUTTER_RUN_LENGTH_COEFFICIENT_CALIBRATED,
        },
        "step_weights": step_weights,
        "locus_rates": locus_rates,
    })
}

/// Un locus que este evento puede hacer tartamudear, ya resuelto.
#[derive(Clone, Debug, PartialEq)]
pub struct StutterLocus {
    id: RepeatLocusId,
    interval: GenomicInterval,
    motifs: Arc<[RepeatMotif]>,
}

impl StutterLocus {
    #[must_use]
    pub fn new(id: RepeatLocusId, interval: GenomicInterval, motifs: Vec<RepeatMotif>) -> Self {
        Self {
            id,
            interval,
            motifs: Arc::from(motifs),
        }
    }

    #[must_use]
    pub const fn id(&self) -> &RepeatLocusId {
        &self.id
    }

    #[must_use]
    pub const fn interval(&self) -> &GenomicInterval {
        &self.interval
    }

    #[must_use]
    pub fn motifs(&self) -> &[RepeatMotif] {
        &self.motifs
    }

    #[must_use]
    pub fn sequence_id(&self) -> &SequenceId {
        self.interval.sequence_id()
    }
}

/// Lo que un evento técnico necesita para tartamudear: sus parámetros y los
/// loci declarados.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StutterPlan {
    settings: StutterSettings,
    loci: Vec<StutterLocus>,
}

impl StutterPlan {
    #[must_use]
    pub fn new(settings: StutterSettings, loci: Vec<StutterLocus>) -> Self {
        Self { settings, loci }
    }

    #[must_use]
    pub const fn settings(&self) -> &StutterSettings {
        &self.settings
    }

    #[must_use]
    pub fn loci(&self) -> &[StutterLocus] {
        &self.loci
    }

    /// Un plan apagado no abre la etapa `stutter` ni consume aleatoriedad.
    ///
    /// Sin loci declarados tampoco hay nada que hacer, aunque el interruptor
    /// esté encendido: ese es exactamente el caso de la congelación byte a
    /// byte de F4.8.6.
    #[must_use]
    pub fn is_active(&self) -> bool {
        self.settings.enabled && !self.loci.is_empty()
    }

    /// Catálogo de ámbitos aleatorios que publica una ejecución con este plan.
    #[must_use]
    pub fn rng_catalog(&self) -> &'static str {
        if self.is_active() {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
        } else {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
        }
    }
}

/// Un paso realizado sobre una molécula concreta.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StutterApplication {
    locus_id: RepeatLocusId,
    step_units: i64,
}

impl StutterApplication {
    #[must_use]
    pub const fn locus_id(&self) -> &RepeatLocusId {
        &self.locus_id
    }

    #[must_use]
    pub const fn step_units(&self) -> i64 {
        self.step_units
    }
}

/// Resultado por pareja. La procedencia es transitoria, como en la quimera:
/// no pertenece al FASTQ ni a la verdad.
#[derive(Clone, Debug, PartialEq)]
pub struct StutterTransformation {
    pair: PairedEndPair,
    evaluated: u64,
    applied: Vec<StutterApplication>,
    /// Pasos netos del linaje de PCR aplicados a esta molécula, **aparte** de
    /// los de la parte de base.
    ///
    /// Van separados porque cuentan cosas distintas y un recuento agregado que
    /// los sumara mentiría en las dos direcciones: la fracción realizada del
    /// tartamudeo de base dejaría de compararse con su esperada, y el número
    /// de resbalones del linaje no es el número de moléculas que acabaron
    /// desplazadas —varios resbalones de una misma copia dan **un** paso neto,
    /// y dos que se cancelan dan **ninguno**—.
    lineage_applied: Vec<StutterApplication>,
}

impl StutterTransformation {
    #[must_use]
    pub const fn pair(&self) -> &PairedEndPair {
        &self.pair
    }

    #[must_use]
    pub fn applied(&self) -> &[StutterApplication] {
        &self.applied
    }

    /// Pasos netos de linaje que de verdad se escribieron en esta molécula.
    #[must_use]
    pub fn lineage_applied(&self) -> &[StutterApplication] {
        &self.lineage_applied
    }

    /// Tractos declarados de esta molécula que llegaron al ensayo de Bernoulli.
    #[must_use]
    pub const fn evaluated(&self) -> u64 {
        self.evaluated
    }

    #[must_use]
    pub const fn transformed(&self) -> bool {
        !self.applied.is_empty()
    }

    #[must_use]
    pub fn into_pair(self) -> PairedEndPair {
        self.pair
    }
}

/// Recuentos **agregados** del evento. Nunca guarda una lectura ni una pareja.
///
/// # Qué cuenta exactamente cada cifra
///
/// Un **tracto evaluado** es el par de una molécula y un tracto declarado, no
/// una molécula: una molécula que contuviera dos tractos declarados aporta dos.
/// Los campos se llaman así porque es lo que cuentan; llamarlos moléculas sería
/// cómodo y falso en cuanto un escenario declare dos loci solapados, cosa que
/// la configuración ya permite. Contar por tracto es lo que hace que la suma de
/// la distribución de pasos valga **exactamente** [`Self::stuttered_tracts`], y
/// no solo al menos.
///
/// # Qué moléculas, que es la otra mitad de la pregunta
///
/// Las **muestreadas**, no las parejas publicadas. Con contaminación biológica
/// los dos conjuntos dejan de coincidir: una receptora retira parejas propias
/// —ya contadas aquí— e incorpora parejas de una donante, que se publican bajo
/// su biblioteca. Es **la misma convención que ya usa el recuento de
/// quimeras**, de modo que no es una incoherencia nueva; pero hay que decirla,
/// porque quien sume las parejas del fichero publicado no va a encontrar esta
/// cifra.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct StutterReport {
    evaluated_tracts: u64,
    stuttered_tracts: u64,
    steps_by_locus: BTreeMap<String, BTreeMap<i64, u64>>,
}

impl StutterReport {
    #[must_use]
    pub const fn evaluated_tracts(&self) -> u64 {
        self.evaluated_tracts
    }

    #[must_use]
    pub const fn stuttered_tracts(&self) -> u64 {
        self.stuttered_tracts
    }

    /// Distribución de los pasos aplicados, por locus y por paso.
    #[must_use]
    pub const fn steps_by_locus(&self) -> &BTreeMap<String, BTreeMap<i64, u64>> {
        &self.steps_by_locus
    }
}

/// Acumula solo cifras agregadas, igual que el de la quimera.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct StutterAccumulator(StutterReport);

impl StutterAccumulator {
    pub fn record(&mut self, transformation: &StutterTransformation) -> Result<(), StutterError> {
        self.0.evaluated_tracts = self
            .0
            .evaluated_tracts
            .checked_add(transformation.evaluated)
            .ok_or(StutterError::CounterOverflow)?;
        for application in &transformation.applied {
            self.0.stuttered_tracts = self
                .0
                .stuttered_tracts
                .checked_add(1)
                .ok_or(StutterError::CounterOverflow)?;
            let entry = self
                .0
                .steps_by_locus
                .entry(application.locus_id.as_str().to_owned())
                .or_default()
                .entry(application.step_units)
                .or_insert(0);
            *entry = entry.checked_add(1).ok_or(StutterError::CounterOverflow)?;
        }
        Ok(())
    }

    /// Suma el informe de una biblioteca al del evento.
    ///
    /// Existe porque el evento recorre una biblioteca por muestreador y cada
    /// muestreador lleva su propio acumulador: sin esto habría que exponer las
    /// parejas, que es justo lo que este tipo evita.
    pub fn merge(&mut self, report: &StutterReport) -> Result<(), StutterError> {
        self.0.evaluated_tracts = self
            .0
            .evaluated_tracts
            .checked_add(report.evaluated_tracts)
            .ok_or(StutterError::CounterOverflow)?;
        self.0.stuttered_tracts = self
            .0
            .stuttered_tracts
            .checked_add(report.stuttered_tracts)
            .ok_or(StutterError::CounterOverflow)?;
        for (locus_id, steps) in &report.steps_by_locus {
            let destination = self.0.steps_by_locus.entry(locus_id.clone()).or_default();
            for (step, count) in steps {
                let entry = destination.entry(*step).or_insert(0);
                *entry = entry
                    .checked_add(*count)
                    .ok_or(StutterError::CounterOverflow)?;
            }
        }
        Ok(())
    }

    #[must_use]
    pub fn report(self) -> StutterReport {
        self.0
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum StutterError {
    InvalidRate {
        motif_length: u64,
        rate: f64,
    },
    InvalidLocusRate {
        locus_id: String,
        rate: f64,
    },
    InvalidStepWeight {
        step: i64,
        weight: f64,
    },
    InvalidRunLengthCoefficient {
        coefficient: f64,
    },
    NoPositiveStepWeight,
    MissingStepWeights,
    UnknownStep {
        step: i64,
    },
    /// La copia haplotípica de la que salió la pareja no está en el individuo.
    /// No es un caso legítimo: la pareja acaba de salir de él.
    CopyNotFound {
        sequence_id: String,
        copy_index: u8,
    },
    LengthOverflow,
    CounterOverflow,
    Random(RandomSamplingError),
    Pair(crate::diagnostics::Diagnostic),
}

/// Aplica el tartamudeo a una pareja ya formada.
///
/// # Dónde encaja, y por qué no puede ir antes
///
/// Se llama **después de la quimera y antes de calidad y errores**. La quimera
/// vuelve a materializar el prefijo del fragmento desde la copia haplotípica y
/// lo compara **byte a byte**; un fragmento ya tartamudeado rompería esa
/// comprobación y devolvería `ChimeraError::PairSourceMismatch`. Después de la
/// quimera la pareja todavía lleva su fragmento completo, así que la molécula
/// sigue disponible y se vuelve a emitir R1 y R2 exactamente como hace la
/// quimera.
///
/// La decisión es **una por molécula y tracto**, y se aplica a la **molécula
/// completa**: nunca por lectura. Es lo que hace que R1 y el complemento
/// inverso de R2 muestren el mismo número de unidades cuando los dos cubren el
/// tracto.
///
/// # Moléculas quiméricas
///
/// Una molécula quimérica se compone de **dos tramos de copia**: el prefijo,
/// de la copia y el inicio del fragmento, y el sufijo, del segundo origen que
/// eligió la quimera. Los tractos se sitúan **en cada tramo con su propia
/// copia** (`chimera`, el origen que entrega la quimera), de modo que se
/// evalúa cada tracto declarado contenido en la molécula y ninguna región que
/// no lo sea (REQ-DAT-74 a 76). El paso heredado del linaje de PCR solo cae
/// sobre los tractos del prefijo, que es la parte que procede de la molécula
/// original de la familia; el sufijo no tiene historia de linaje.
///
/// Hasta el 27 de octubre de 2026 la ventana de la primera copia se aplicaba a
/// la molécula entera: más allá del punto de corte leía bases del segundo
/// origen, que podían tartamudear sin ser un tracto, y un tracto del segundo
/// origen no se evaluaba nunca (DT-060).
#[expect(
    clippy::too_many_arguments,
    reason = "el origen de la quimera es el octavo hecho de la pareja; agruparlo con la pareja movería la lista sin acortarla"
)]
pub fn transform_paired_end(
    root: RootSeed,
    event_id: &EventId,
    individual: &Individual,
    plan: &StutterPlan,
    lineage: Option<&crate::clonal_stutter::LineageNetSteps>,
    chimera: Option<&crate::chimera::ChimeraOrigin>,
    pair: PairedEndPair,
    read_length: u16,
) -> Result<StutterTransformation, StutterError> {
    if !plan.is_active() {
        return Ok(StutterTransformation {
            pair,
            evaluated: 0,
            applied: Vec::new(),
            lineage_applied: Vec::new(),
        });
    }
    let template = pair.fragment().template();
    let ordinal = template.template_ordinal().get();
    let library_id = template.library_id().clone();
    let copy_id = pair.fragment().haplotype_copy_id().clone();

    // La secuencia solo se copia **cuando de verdad hay que reescribirla**. Con
    // la tasa calibrada del dinucleótido, 0,011, la inmensa mayoría de las
    // moléculas no tartamudea, y copiar sus cuatrocientas bases por si acaso
    // sería una asignación por molécula a cambio de nada.
    let mut rewritten: Option<Vec<u8>> = None;
    let mut evaluated = 0_u64;
    let mut applied = Vec::new();
    let mut lineage_applied = Vec::new();
    let original_length = u64::try_from(pair.fragment().sequence().len())
        .map_err(|_| StutterError::LengthOverflow)?;
    // Los tractos se recorren en **orden canónico de inicio sobre la
    // molécula**, con desempate por identificador de locus. El orden importa
    // porque un paso desplaza todo lo que queda a su derecha.
    let placements = match chimera {
        None => place_tracts(
            individual,
            plan,
            &copy_id,
            pair.fragment().start(),
            original_length,
        )?,
        Some(origin) => {
            let breakpoint = origin.breakpoint().get();
            let second = origin.second();
            place_tracts_in_segments(
                individual,
                plan,
                &[
                    MoleculeSegment {
                        copy_id: origin.first().copy_id(),
                        copy_start: origin.first().start(),
                        molecule_start: 0,
                        length: breakpoint,
                    },
                    MoleculeSegment {
                        copy_id: second.copy_id(),
                        copy_start: second.start(),
                        molecule_start: breakpoint,
                        length: second.end() - second.start(),
                    },
                ],
            )?
        }
    };

    // Desplazamiento acumulado por los pasos ya aplicados. Todos caen **a la
    // izquierda** del tracto siguiente, porque los tractos se recorren en
    // orden, de modo que basta con sumarlo a la ventana.
    let mut shift = 0_i64;
    for placed in placements {
        let locus = placed.locus();
        let (raw_start, raw_end) = placed.window();
        let window_start = usize::try_from(
            raw_start
                .checked_add_signed(shift)
                .ok_or(StutterError::LengthOverflow)?,
        )
        .map_err(|_| StutterError::LengthOverflow)?;
        let mut window_end = usize::try_from(
            raw_end
                .checked_add_signed(shift)
                .ok_or(StutterError::LengthOverflow)?,
        )
        .map_err(|_| StutterError::LengthOverflow)?;

        // ------------------------------------------------------------------
        // 1. **Primero** el paso neto heredado del linaje de PCR.
        // ------------------------------------------------------------------
        let net = if placed.from_first_segment {
            lineage.map_or(0, |steps| steps.step_for(locus.id().as_str()))
        } else {
            0
        };
        if net != 0 {
            let sequence: &[u8] = rewritten
                .as_deref()
                .unwrap_or_else(|| pair.fragment().sequence());
            let length = u64::try_from(sequence.len()).map_err(|_| StutterError::LengthOverflow)?;
            if window_end > sequence.len() {
                return Err(StutterError::LengthOverflow);
            }
            let observed = observe_tract(&sequence[window_start..window_end], locus.motifs());
            let run_units = observed.longest_run_units();
            let motif_length = observed.longest_run_motif_length();
            // El paso se sorteó sobre la molécula original de la familia. Aquí
            // se vuelve a comprobar contra la molécula **que se va a emitir**,
            // que con una quimera puede no ser la misma. Si no cabe, no se
            // aplica; no hay nada que devolver al generador porque en este
            // punto ya no se consume aleatoriedad.
            if run_units > 0
                && step_is_viable(run_units, motif_length, length, u64::from(read_length), net)
            {
                let run_start = window_start
                    + usize::try_from(observed.longest_run_start())
                        .map_err(|_| StutterError::LengthOverflow)?;
                let target = match rewritten {
                    Some(ref mut sequence) => sequence,
                    None => rewritten.insert(pair.fragment().sequence().to_vec()),
                };
                apply_step_to_molecule(target, run_start, run_units, motif_length, net)?;
                let delta = net
                    .checked_mul(
                        i64::try_from(motif_length).map_err(|_| StutterError::LengthOverflow)?,
                    )
                    .ok_or(StutterError::LengthOverflow)?;
                // La ventana de **este** tracto crece o mengua con su propio
                // paso, además de desplazar todo lo que queda a su derecha.
                window_end = usize::try_from(
                    u64::try_from(window_end)
                        .map_err(|_| StutterError::LengthOverflow)?
                        .checked_add_signed(delta)
                        .ok_or(StutterError::LengthOverflow)?,
                )
                .map_err(|_| StutterError::LengthOverflow)?;
                shift = shift
                    .checked_add(delta)
                    .ok_or(StutterError::LengthOverflow)?;
                lineage_applied.push(StutterApplication {
                    locus_id: locus.id().clone(),
                    step_units: net,
                });
            }
        }

        // ------------------------------------------------------------------
        // 2. **Después** la decisión de base por pareja, sobre la molécula ya
        //    actualizada. La tirada más larga se vuelve a leer aquí: el paso
        //    de linaje puede haberla movido, y con un tracto impuro puede
        //    incluso haber dejado otra tirada como la más larga.
        // ------------------------------------------------------------------
        let sequence: &[u8] = rewritten
            .as_deref()
            .unwrap_or_else(|| pair.fragment().sequence());
        let length = u64::try_from(sequence.len()).map_err(|_| StutterError::LengthOverflow)?;
        if window_end > sequence.len() {
            return Err(StutterError::LengthOverflow);
        }
        let observed = observe_tract(&sequence[window_start..window_end], locus.motifs());
        let run_units = observed.longest_run_units();
        if run_units == 0 {
            // No hay ninguna copia completa dentro de la molécula: no hay
            // tirada donde resbalar y no se consume aleatoriedad.
            continue;
        }
        let motif_length = observed.longest_run_motif_length();
        let Some(base_rate) =
            plan.settings()
                .modulated_rate_for(locus.id().as_str(), motif_length, run_units)
        else {
            // Longitud de motivo fuera del tramo calibrado. No se extrapola y
            // no se consume aleatoriedad; la admisión ya lo rechaza antes de
            // llegar aquí, y esta rama es su red.
            continue;
        };
        let run_start = window_start
            + usize::try_from(observed.longest_run_start())
                .map_err(|_| StutterError::LengthOverflow)?;
        // La longitud que cuenta es la **actual** de la molécula, no la que
        // tenía antes de un paso anterior: la guarda de la lectura mira lo que
        // se va a emitir.
        let candidates = viable_steps(run_units, motif_length, length, u64::from(read_length));
        if candidates.is_empty() {
            // **No puede ocurrir con una tirada de al menos una unidad**: la
            // expansión de una unidad ni vacía la tirada ni acorta la
            // molécula, de modo que siempre sobrevive al filtrado. La guarda
            // queda porque el coste es nulo y porque quien lea el filtrado
            // tiene derecho a ver qué pasaría si algún día dejara de ser
            // cierto. Ninguna prueba puede ejercitarla.
            continue;
        }
        evaluated += 1;
        let probability = slip_probability(
            base_rate,
            plan.settings().run_length_coefficient(),
            run_units,
        );
        let mut bernoulli_rng = derive_f4_technical_subseed(
            root,
            F4TechnicalDecision::stutter_bernoulli(event_id, &library_id, ordinal, locus.id()),
        )
        .map_err(StutterError::Pair)?
        .local_rng();
        if !f4_bernoulli(&mut bernoulli_rng, probability).map_err(StutterError::Random)? {
            continue;
        }
        let weights = plan
            .settings()
            .step_weights()
            .ok_or(StutterError::MissingStepWeights)?;
        let mut step_rng = derive_f4_technical_subseed(
            root,
            F4TechnicalDecision::stutter_step(event_id, &library_id, ordinal, locus.id()),
        )
        .map_err(StutterError::Pair)?
        .local_rng();
        let Some(step) = choose_step_with_source(&mut step_rng, weights, &candidates)? else {
            continue;
        };
        let target = match rewritten {
            Some(ref mut sequence) => sequence,
            None => rewritten.insert(pair.fragment().sequence().to_vec()),
        };
        apply_step_to_molecule(target, run_start, run_units, motif_length, step)?;
        shift = shift
            .checked_add(step.saturating_mul(
                i64::try_from(motif_length).map_err(|_| StutterError::LengthOverflow)?,
            ))
            .ok_or(StutterError::LengthOverflow)?;
        applied.push(StutterApplication {
            locus_id: locus.id().clone(),
            step_units: step,
        });
    }

    let Some(sequence) = rewritten else {
        return Ok(StutterTransformation {
            pair,
            evaluated,
            applied,
            lineage_applied,
        });
    };
    let fragment =
        PairedEndFragment::new(template.clone(), copy_id, pair.fragment().start(), sequence)
            .map_err(StutterError::Pair)?;
    let transformed = emit_paired_end(fragment, pair.pre_hop_identity().clone(), read_length)
        .map_err(StutterError::Pair)?;
    Ok(StutterTransformation {
        pair: transformed,
        evaluated,
        applied,
        lineage_applied,
    })
}

/// Probabilidad de resbalón **derivada de las bases** de un tramo de molécula.
///
/// Es el camino real del transformador, expuesto aparte para que una prueba
/// pueda comparar dos tractos sin montar una ejecución entera. Devuelve `None`
/// cuando no hay ninguna copia completa o cuando la longitud de motivo cae
/// fuera del tramo calibrado.
///
/// **El predictor sale de aquí**: es la tirada ininterrumpida más larga que
/// [`observe_tract`] encuentra, no el número total de unidades.
#[must_use]
pub fn observed_slip_probability(
    settings: &StutterSettings,
    locus_id: &str,
    bases: &[u8],
    motifs: &[RepeatMotif],
) -> Option<f64> {
    let observed = observe_tract(bases, motifs);
    if observed.longest_run_units() == 0 {
        return None;
    }
    let base_rate = settings.modulated_rate_for(
        locus_id,
        observed.longest_run_motif_length(),
        observed.longest_run_units(),
    )?;
    Some(slip_probability(
        base_rate,
        settings.run_length_coefficient(),
        observed.longest_run_units(),
    ))
}

/// Elige el paso entre los candidatos **ya filtrados**.
///
/// Vive aparte de [`transform_paired_end`] para que una prueba pueda **contar**
/// las palabras que consume: con un solo candidato tienen que ser cero, y sin
/// poder contarlas esa afirmación no sería comprobable. Es el mismo motivo por
/// el que el sorteo del alelo de repetición tiene su propia versión
/// instrumentable.
pub(crate) fn choose_step_with_source(
    rng: &mut impl crate::randomness::RandomWordSource,
    weights: &StutterStepWeights,
    candidates: &[i64],
) -> Result<Option<i64>, StutterError> {
    let candidate_weights = candidates
        .iter()
        .map(|step| {
            weights
                .weight_of(*step)
                .ok_or(StutterError::UnknownStep { step: *step })
        })
        .collect::<Result<Vec<_>, _>>()?;
    if !candidate_weights.iter().any(|weight| *weight > 0.0) {
        // Todos los candidatos viables tienen peso cero: el usuario ha
        // declarado que esos pasos no ocurren. No hay nada que elegir y no se
        // sortea.
        return Ok(None);
    }
    let chosen =
        categorical_f64_with_source(rng, &candidate_weights).map_err(StutterError::Random)?;
    Ok(Some(candidates[chosen]))
}

/// Los pasos que sobreviven al filtrado previo, en el orden canónico.
///
/// Un paso entra si deja al menos una unidad completa en la tirada —la misma
/// regla que `TandemRepeatTract::apply_unit_step`— y si deja la molécula al
/// menos tan larga como una lectura, porque emitir la pareja exige que el
/// inserto no sea más corto que el mate.
#[must_use]
pub fn viable_steps(
    run_units: u64,
    motif_length: u64,
    molecule_length: u64,
    read_length: u64,
) -> Vec<i64> {
    STUTTER_ADMITTED_STEPS
        .into_iter()
        .filter(|step| step_is_viable(run_units, motif_length, molecule_length, read_length, *step))
        .collect()
}

/// La regla de viabilidad, para un paso de **cualquier** tamaño.
///
/// [`viable_steps`] la aplica a los cuatro pasos de un ciclo; el tartamudeo
/// clonal la aplica al **paso neto acumulado** de un linaje, que puede valer
/// cualquier entero. **Es una sola función a propósito**: con dos copias, un
/// día una admitiría lo que la otra prohíbe y el modelo tendría dos reglas
/// distintas para el mismo fenómeno.
#[must_use]
pub fn step_is_viable(
    run_units: u64,
    motif_length: u64,
    molecule_length: u64,
    read_length: u64,
    step: i64,
) -> bool {
    if step >= 0 {
        return true;
    }
    let removed = step.unsigned_abs();
    if removed >= run_units {
        return false;
    }
    molecule_length
        .checked_sub(removed.saturating_mul(motif_length))
        .is_some_and(|length| length >= read_length)
}

/// Un tracto declarado ya situado dentro de una molécula concreta.
///
/// Los desplazamientos son los de la molécula **antes** de aplicar ningún
/// paso, y el orden de la lista es el canónico de inicio con desempate por
/// identificador de locus.
#[derive(Clone, Copy, Debug)]
pub struct PlacedTract<'a> {
    locus: &'a StutterLocus,
    start: u64,
    end: u64,
    /// El tracto sale del primer tramo de la molécula, el que procede de la
    /// molécula original de la familia. En una molécula no quimérica, siempre.
    from_first_segment: bool,
}

impl<'a> PlacedTract<'a> {
    #[must_use]
    pub const fn locus(&self) -> &'a StutterLocus {
        self.locus
    }

    /// Ventana `[inicio, fin)` dentro de la molécula.
    #[must_use]
    pub const fn window(&self) -> (u64, u64) {
        (self.start, self.end)
    }
}

/// Sitúa en una molécula los tractos declarados que contiene, en orden
/// canónico.
///
/// Vive aquí y no duplicado en el tartamudeo clonal porque las dos partes
/// tienen que ver **los mismos tractos en el mismo orden**: si el linaje
/// ordenara de otra manera, el paso neto y la decisión de base caerían sobre
/// tractos distintos en cuanto hubiera dos.
pub fn place_tracts<'a>(
    individual: &Individual,
    plan: &'a StutterPlan,
    copy_id: &crate::haplotype::HaplotypeCopyId,
    fragment_start: u64,
    fragment_length: u64,
) -> Result<Vec<PlacedTract<'a>>, StutterError> {
    place_tracts_in_segments(
        individual,
        plan,
        &[MoleculeSegment {
            copy_id,
            copy_start: fragment_start,
            molecule_start: 0,
            length: fragment_length,
        }],
    )
}

/// Un tramo contiguo de copia dentro de una molécula.
struct MoleculeSegment<'b> {
    copy_id: &'b crate::haplotype::HaplotypeCopyId,
    copy_start: u64,
    molecule_start: u64,
    length: u64,
}

/// Sitúa los tractos de cada tramo con su propia copia y los ordena en el
/// orden canónico de la molécula entera.
fn place_tracts_in_segments<'a>(
    individual: &Individual,
    plan: &'a StutterPlan,
    segments: &[MoleculeSegment<'_>],
) -> Result<Vec<PlacedTract<'a>>, StutterError> {
    let mut placements = Vec::new();
    for (index, segment) in segments.iter().enumerate() {
        let copy_id = segment.copy_id;
        let Some(copy) = individual
            .copies_for(copy_id.sequence_id())
            .and_then(|copies| copies.iter().find(|candidate| *candidate.id() == *copy_id))
        else {
            return Err(StutterError::CopyNotFound {
                sequence_id: copy_id.sequence_id().as_str().to_owned(),
                copy_index: copy_id.copy_index(),
            });
        };
        for locus in plan.loci() {
            if locus.sequence_id() != copy_id.sequence_id() {
                continue;
            }
            let Some(placement) = locate_tract_in_copy(copy, locus)? else {
                continue;
            };
            let Some((start, end)) =
                clip_to_molecule(placement, segment.copy_start, segment.length)
            else {
                continue;
            };
            placements.push(PlacedTract {
                locus,
                start: start + segment.molecule_start,
                end: end + segment.molecule_start,
                from_first_segment: index == 0,
            });
        }
    }
    placements.sort_by(|left, right| {
        (left.start, left.locus.id().as_str()).cmp(&(right.start, right.locus.id().as_str()))
    });
    Ok(placements)
}

/// Dónde cae el tracto dentro de la copia haplotípica.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct TractPlacement {
    start: u64,
    end: u64,
}

fn clip_to_molecule(
    placement: TractPlacement,
    fragment_start: u64,
    fragment_length: u64,
) -> Option<(u64, u64)> {
    let fragment_end = fragment_start.checked_add(fragment_length)?;
    let start = placement.start.max(fragment_start);
    let end = placement.end.min(fragment_end);
    if start >= end {
        return None;
    }
    Some((start - fragment_start, end - fragment_start))
}

/// Traduce el intervalo de referencia del locus a coordenadas **vigentes** de
/// la copia, recorriendo sus piezas.
///
/// Hay dos casos y los dos son normales:
///
/// - la copia lleva un alelo de repetición aplicado, y entonces la pieza es
///   suya y su longitud es la del tracto realizado;
/// - el alelo sorteado coincidía byte a byte con el de la referencia, así que
///   no se aplicó nada y el tracto vive dentro de una pieza de fuente. **Este
///   segundo caso es justo el de la prueba de congelación acotada**, de modo
///   que ignorarlo dejaría sin tartamudear precisamente el escenario que el
///   plan usa para comprobar que el tartamudeo ocurre.
fn locate_tract_in_copy(
    copy: &HaplotypeCopy,
    locus: &StutterLocus,
) -> Result<Option<TractPlacement>, StutterError> {
    let interval = locus.interval();
    let mut offset = 0_u64;
    for piece in copy.segments().current_modifications() {
        match piece {
            CurrentIndelPiece::RepeatAllele(allele) if allele.locus_id() == locus.id() => {
                let length = allele.tract().length();
                let end = offset
                    .checked_add(length)
                    .ok_or(StutterError::LengthOverflow)?;
                return Ok(Some(TractPlacement { start: offset, end }));
            }
            CurrentIndelPiece::Source(source)
                if source.source().sequence_id() == interval.sequence_id()
                    && source.source().start() <= interval.start()
                    && interval.end() <= source.source().end() =>
            {
                let start = offset
                    .checked_add(interval.start() - source.source().start())
                    .ok_or(StutterError::LengthOverflow)?;
                let end = start
                    .checked_add(interval.len())
                    .ok_or(StutterError::LengthOverflow)?;
                return Ok(Some(TractPlacement { start, end }));
            }
            _ => {}
        }
        offset = offset
            .checked_add(piece_length(piece)?)
            .ok_or(StutterError::LengthOverflow)?;
    }
    Ok(None)
}

fn piece_length(piece: &CurrentIndelPiece) -> Result<u64, StutterError> {
    match piece {
        CurrentIndelPiece::Source(source) => Ok(source.source().len()),
        CurrentIndelPiece::Inserted(inserted) => {
            u64::try_from(inserted.bases().len()).map_err(|_| StutterError::LengthOverflow)
        }
        CurrentIndelPiece::RepeatAllele(allele) => Ok(allele.tract().length()),
    }
}

/// Aplica el paso **en la tirada más larga**, con la misma convención que
/// `TandemRepeatTract::apply_unit_step`: las copias se añaden o se quitan por
/// el **extremo derecho de esa tirada**, no por el extremo derecho del tracto.
///
/// La molécula cambia de longitud en `paso × longitud del motivo`, y con ella
/// el inserto observado. **No se compensa ni se descarta la molécula**: es lo
/// que ocurre físicamente, porque la molécula ya estaba seleccionada por tamaño
/// cuando la polimerasa la copió.
fn apply_step_to_molecule(
    sequence: &mut Vec<u8>,
    run_start: usize,
    run_units: u64,
    motif_length: u64,
    step: i64,
) -> Result<(), StutterError> {
    let motif_length = usize::try_from(motif_length).map_err(|_| StutterError::LengthOverflow)?;
    let run_units = usize::try_from(run_units).map_err(|_| StutterError::LengthOverflow)?;
    let run_end = run_start
        .checked_add(
            run_units
                .checked_mul(motif_length)
                .ok_or(StutterError::LengthOverflow)?,
        )
        .ok_or(StutterError::LengthOverflow)?;
    if run_end > sequence.len() {
        return Err(StutterError::LengthOverflow);
    }
    if step > 0 {
        let added = usize::try_from(step).map_err(|_| StutterError::LengthOverflow)?;
        let motif = sequence[run_start..run_start + motif_length].to_vec();
        let mut inserted = Vec::with_capacity(added * motif_length);
        for _ in 0..added {
            inserted.extend_from_slice(&motif);
        }
        sequence.splice(run_end..run_end, inserted);
    } else {
        let removed =
            usize::try_from(step.unsigned_abs()).map_err(|_| StutterError::LengthOverflow)?;
        let cut = run_end
            .checked_sub(
                removed
                    .checked_mul(motif_length)
                    .ok_or(StutterError::LengthOverflow)?,
            )
            .ok_or(StutterError::LengthOverflow)?;
        sequence.drain(cut..run_end);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        STUTTER_ADMITTED_STEPS, STUTTER_DEFAULT_MOTIF_LENGTH_RATES, StutterMotifLengthRates,
        StutterSettings, StutterStepWeights, apply_step_to_molecule, choose_step_with_source,
        observed_slip_probability, viable_steps,
    };
    use crate::randomness::RandomWordSource;
    use crate::tandem_repeat::RepeatMotif;

    /// DT-034: la licencia CC BY 4.0 de Gymrek (2016) exige atribuir las cifras,
    /// y la forma del modelo es de Raz (2019); el texto publicado cita a ambos.
    #[test]
    fn calibration_source_credits_the_rates_and_the_model_form() {
        let source = super::STUTTER_CALIBRATION_SOURCE;
        assert!(source.contains("Gymrek (2016)"), "{source}");
        assert!(source.contains("table 1"), "{source}");
        assert!(source.contains("Raz et al. (2019)"), "{source}");
        assert_eq!(source.matches("CC BY 4.0").count(), 2, "{source}");
    }

    /// Fuente instrumentada: cuenta las palabras que se le piden.
    ///
    /// Sin ella, «no consume aleatoriedad» sería una afirmación que ninguna
    /// prueba podría distinguir de «consume y da el mismo resultado».
    struct Counting {
        words: u64,
    }

    impl RandomWordSource for Counting {
        fn next_word(&mut self) -> u64 {
            self.words += 1;
            0
        }
    }

    fn motifs(text: &str) -> Vec<RepeatMotif> {
        vec![RepeatMotif::new(text.as_bytes()).expect("motivo")]
    }

    fn settings(coefficient: f64) -> StutterSettings {
        StutterSettings::new(
            true,
            StutterMotifLengthRates::default(),
            coefficient,
            Some(StutterStepWeights::new([0.05, 0.15, 0.6, 0.2]).expect("pesos")),
            std::collections::BTreeMap::new(),
        )
        .expect("parámetros")
    }

    /// **(a) El predictor es la tirada más larga, no el total de unidades.**
    ///
    /// Los dos tractos tienen **diez unidades completas** de `CA`. El primero
    /// las lleva en una sola tirada; el segundo, partidas en seis y cuatro por
    /// una interrupción `TT`. Un modelo cuyo predictor fuera el número total
    /// de unidades daría **la misma** probabilidad a los dos, y esta es la
    /// única prueba que lo caza.
    ///
    /// **Cuidado, y por eso el fixture es impuro:** con tractos puros el
    /// número total de unidades y la tirada más larga **coinciden**, de modo
    /// que un fixture puro no discriminaría entre los dos predictores. El
    /// total se ha igualado a propósito.
    ///
    /// El coeficiente es positivo porque con el coeficiente en cero la
    /// probabilidad no depende de la tirada en absoluto y la aserción no
    /// podría fallar.
    #[test]
    fn longest_run_predicts_more_slippage_than_a_split_tract_of_equal_total_units() {
        let settings = settings(0.001);
        let motifs = motifs("CA");
        let single_run = "CA".repeat(10);
        let split = format!("{}TT{}", "CA".repeat(6), "CA".repeat(4));
        let long = observed_slip_probability(&settings, "STR-1", single_run.as_bytes(), &motifs)
            .expect("tirada larga");
        let short = observed_slip_probability(&settings, "STR-1", split.as_bytes(), &motifs)
            .expect("tirada partida");
        // El total es el mismo en los dos: diez copias completas.
        assert_eq!(
            crate::tandem_repeat::observe_tract(single_run.as_bytes(), &motifs).covered_bases(),
            crate::tandem_repeat::observe_tract(split.as_bytes(), &motifs).covered_bases(),
        );
        assert!(
            long > short,
            "la tirada de diez tiene que resbalar más que la de seis: {long} frente a {short}"
        );
    }

    /// **(b) La longitud del motivo manda sobre los pares de base.**
    ///
    /// Veinte pares de base en las dos, mononucleótido contra
    /// tetranucleótido. Un modelo que dependiera solo de los pares de base
    /// daría la misma cifra. Con el coeficiente en cero la probabilidad
    /// **reproduce exactamente** la tabla calibrada, y eso también se
    /// comprueba aquí.
    #[test]
    fn motif_length_drives_the_rate_beyond_the_tract_length_in_base_pairs() {
        let settings = settings(0.0);
        let mononucleotide =
            observed_slip_probability(&settings, "STR-1", "A".repeat(20).as_bytes(), &motifs("A"))
                .expect("mononucleótido");
        let tetranucleotide = observed_slip_probability(
            &settings,
            "STR-1",
            "ACGT".repeat(5).as_bytes(),
            &motifs("ACGT"),
        )
        .expect("tetranucleótido");
        assert!(mononucleotide > tetranucleotide);
        assert_eq!(mononucleotide, STUTTER_DEFAULT_MOTIF_LENGTH_RATES[0]);
        assert_eq!(tetranucleotide, STUTTER_DEFAULT_MOTIF_LENGTH_RATES[3]);
    }

    /// **(c) La estructura del paso es exactamente `{−3,−2,−1,+1}`.**
    #[test]
    fn the_admitted_step_structure_is_asymmetric_and_closed() {
        assert_eq!(STUTTER_ADMITTED_STEPS, [-3, -2, -1, 1]);
        let weights = StutterStepWeights::new([0.05, 0.15, 0.6, 0.2]).expect("pesos");
        assert_eq!(weights.weight_of(2), None);
        assert_eq!(weights.weight_of(-4), None);
    }

    /// **(d) Con un solo candidato viable el sorteo no consume aleatoriedad.**
    ///
    /// Una tirada de **una sola unidad completa** deja fuera las tres
    /// contracciones, porque cualquiera de ellas la vaciaría, y solo sobrevive
    /// la expansión de una unidad. Una implementación que sorteara entre los
    /// cuatro pasos y **descartara** los inviables consumiría al menos una
    /// palabra: esta prueba cuenta el consumo, no el resultado, y es la única
    /// que distingue el filtrado previo del descarte posterior.
    #[test]
    fn a_single_viable_step_is_chosen_without_consuming_randomness() {
        let candidates = viable_steps(1, 2, 400, 150);
        assert_eq!(candidates, vec![1]);
        let weights = StutterStepWeights::new([0.05, 0.15, 0.6, 0.2]).expect("pesos");
        let mut counter = Counting { words: 0 };
        assert_eq!(
            choose_step_with_source(&mut counter, &weights, &candidates).expect("paso"),
            Some(1)
        );
        assert_eq!(counter.words, 0, "el candidato único no se sortea");
    }

    /// Con varios candidatos sí se consume: es el control positivo de la
    /// prueba anterior. Sin él, un sorteo que nunca consumiera nada la
    /// satisfaría igual.
    #[test]
    fn several_viable_steps_do_consume_randomness() {
        let candidates = viable_steps(4, 2, 400, 150);
        assert_eq!(candidates, vec![-3, -2, -1, 1]);
        let weights = StutterStepWeights::new([0.05, 0.15, 0.6, 0.2]).expect("pesos");
        let mut counter = Counting { words: 0 };
        let chosen = choose_step_with_source(&mut counter, &weights, &candidates).expect("paso");
        assert!(chosen.is_some());
        assert!(counter.words > 0, "cuatro candidatos exigen sortear");
    }

    /// El filtrado previo excluye la contracción que vaciaría la tirada y la
    /// que dejaría la molécula más corta que una lectura.
    #[test]
    fn viable_steps_exclude_contractions_that_empty_the_run_or_shorten_the_molecule() {
        // Tirada de tres: quitar tres la vaciaría.
        assert_eq!(viable_steps(3, 2, 400, 150), vec![-2, -1, 1]);
        // Molécula de 155 pb con motivo de 4 pb: quitar una la deja en 151 y
        // cabe; quitar dos la dejaría en 147, por debajo de la lectura de 150.
        assert_eq!(viable_steps(10, 4, 155, 150), vec![-1, 1]);
    }

    /// **F4.8.6 (b). Cuándo se abre la etapa `stutter`, y cuándo no.**
    ///
    /// El transformador sale por la primera línea cuando el plan no está
    /// activo, de modo que «cero consumo en la etapa» y «plan inactivo» son la
    /// misma afirmación. Un escenario **sin loci declarados** deja el plan
    /// inactivo aunque el interruptor esté encendido, y publica el catálogo
    /// vigente: eso es lo que mantiene byte a byte lo ya publicado.
    #[test]
    fn the_stage_only_opens_with_the_switch_on_and_a_declared_locus() {
        let locus = super::StutterLocus::new(
            crate::identity::RepeatLocusId::parse("STR-1", "test").expect("locus"),
            crate::coordinates::GenomicInterval::new(
                crate::identity::ReferenceId::parse("ref-main", "test").expect("referencia"),
                crate::identity::SequenceId::parse("chr1", "test").expect("secuencia"),
                100,
                130,
                1000,
            )
            .expect("intervalo"),
            motifs("CA"),
        );
        let off = super::StutterPlan::default();
        assert!(!off.is_active());
        assert_eq!(
            off.rng_catalog(),
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
        );

        let on_without_loci = super::StutterPlan::new(settings(0.0), Vec::new());
        assert!(!on_without_loci.is_active());
        assert_eq!(
            on_without_loci.rng_catalog(),
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
        );

        let active = super::StutterPlan::new(settings(0.0), vec![locus]);
        assert!(active.is_active());
        assert_eq!(
            active.rng_catalog(),
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
        );
    }

    /// **F4.8.5. La molécula cambia de longitud y no se compensa.**
    ///
    /// Una molécula que pierde **una unidad de un motivo de 4 pb** produce un
    /// fragmento exactamente 4 pb más corto. La aserción mira **la molécula**:
    /// observar las lecturas no discriminaría nada, porque siguen midiendo lo
    /// mismo pase lo que pase. Una implementación que recortara o rellenara
    /// para conservar la longitud sorteada muere aquí.
    #[test]
    fn a_contraction_shortens_the_molecule_by_the_motif_length() {
        let mut sequence = format!("GGGG{}GGGG", "ACGT".repeat(6)).into_bytes();
        let before = sequence.len();
        apply_step_to_molecule(&mut sequence, 4, 6, 4, -1).expect("contracción");
        assert_eq!(before - sequence.len(), 4);
        assert_eq!(
            String::from_utf8(sequence.clone()).expect("texto"),
            format!("GGGG{}GGGG", "ACGT".repeat(5))
        );
        // Y la expansión la alarga en la misma cantidad, no en otra.
        let before = sequence.len();
        apply_step_to_molecule(&mut sequence, 4, 5, 4, 1).expect("expansión");
        assert_eq!(sequence.len() - before, 4);
    }

    // -----------------------------------------------------------------
    // Dependencia con la longitud del tracto
    // -----------------------------------------------------------------

    /// Los valores digitalizados de la figura, **copiados aquí a mano** desde
    /// la sección 9 del documento de calibración.
    ///
    /// Son el oráculo independiente de la tabla del módulo: viven aquí como
    /// literales, no se importan de `super`, y por eso la prueba que los usa
    /// puede distinguir una tabla mal generada. Si se importaran, la
    /// comparación sería circular y no cazaría nada.
    const PUBLISHED_MONONUCLEOTIDE_BASES: [u64; 6] = [5, 7, 9, 11, 13, 15];
    const PUBLISHED_MONONUCLEOTIDE_WITHOUT_PCR: [f64; 6] = [-5.4, -4.1, -3.0, -2.2, -1.65, -1.3];
    const PUBLISHED_MONONUCLEOTIDE_WITH_PCR: [f64; 6] = [-5.3, -3.8, -2.7, -1.8, -1.25, -0.85];
    const PUBLISHED_DINUCLEOTIDE_BASES: [u64; 4] = [10, 18, 26, 34];
    const PUBLISHED_DINUCLEOTIDE_WITHOUT_PCR: [f64; 4] = [-3.6, -2.3, -1.6, -1.25];
    const PUBLISHED_DINUCLEOTIDE_WITH_PCR: [f64; 4] = [-3.4, -2.0, -1.35, -0.9];

    /// Interpola en escala logarítmica, escrito aquí aparte del módulo.
    fn published_log10(grid: &[u64], values: &[f64], bases: u64) -> f64 {
        let bases = bases as f64;
        if bases <= grid[0] as f64 {
            return values[0];
        }
        if bases >= grid[grid.len() - 1] as f64 {
            return values[values.len() - 1];
        }
        for window in 0..grid.len() - 1 {
            let (low, high) = (grid[window] as f64, grid[window + 1] as f64);
            if low <= bases && bases <= high {
                let fraction = (bases - low) / (high - low);
                return fraction.mul_add(values[window + 1] - values[window], values[window]);
            }
        }
        unreachable!("la rejilla cubre todo el tramo")
    }

    /// **La tabla del módulo es la curva publicada, normalizada en el ancla.**
    ///
    /// Rederiva **cada** entrada desde los valores de la figura escritos a
    /// mano arriba: media logarítmica de las dos columnas, interpolación
    /// lineal en escala logarítmica y resta del valor del ancla.
    ///
    /// **Debe fallar contra** una tabla generada con otra ancla, con una sola
    /// columna, con interpolación lineal en escala lineal, o con una entrada
    /// cambiada a mano.
    #[test]
    fn the_tract_length_table_is_the_published_curve_normalised_at_the_anchor() {
        for (motif_length, grid, without_pcr, with_pcr, anchor_bases) in [
            (
                1_u64,
                &PUBLISHED_MONONUCLEOTIDE_BASES[..],
                &PUBLISHED_MONONUCLEOTIDE_WITHOUT_PCR[..],
                &PUBLISHED_MONONUCLEOTIDE_WITH_PCR[..],
                15_u64,
            ),
            (
                2,
                &PUBLISHED_DINUCLEOTIDE_BASES[..],
                &PUBLISHED_DINUCLEOTIDE_WITHOUT_PCR[..],
                &PUBLISHED_DINUCLEOTIDE_WITH_PCR[..],
                34,
            ),
        ] {
            assert_eq!(
                super::STUTTER_TRACT_LENGTH_REFERENCE_BASES[motif_length as usize - 1],
                Some(anchor_bases),
                "el ancla publicada tiene que ser el tope del tramo medido"
            );
            let mean = |bases: u64| {
                (published_log10(grid, without_pcr, bases) + published_log10(grid, with_pcr, bases))
                    / 2.0
            };
            let anchor = mean(anchor_bases);
            let (first, last) =
                super::tract_length_measured_units(motif_length).expect("tramo medido");
            for units in first..=last {
                let expected = 10_f64.powf(mean(units * motif_length) - anchor);
                let observed = super::tract_length_factor(motif_length, units);
                assert!(
                    (observed - expected).abs() <= expected * 1e-12,
                    "motivo de {motif_length} pb, {units} unidades: \
                     la tabla dice {observed} y la curva publicada {expected}"
                );
            }
        }
    }

    /// **La modulación toma la TIRADA ININTERRUMPIDA MÁS LARGA, no el total de
    /// unidades del tracto.**
    ///
    /// # Por qué es esa y no la otra, y por qué importa
    ///
    /// Porque es lo que la fuente midió. Fungtammasan y colaboradores dicen
    /// en sus métodos que **«solo se consideraron STR ininterrumpidos
    /// (perfectos)»**, de modo que el eje de su figura es la longitud de una
    /// tirada sin interrupciones. Aplicar su curva al total de unidades de un
    /// tracto impuro le atribuiría una medida que no hizo. Y coincide con el
    /// predictor que DTC-090 ya había fijado para la parte lineal, de modo que
    /// el modelo no tiene dos predictores distintos para el mismo fenómeno.
    ///
    /// # El fixture es IMPURO a propósito
    ///
    /// Con un tracto puro el total de unidades y la tirada más larga
    /// **coinciden**, y ninguna prueba con tractos puros puede distinguir las
    /// dos lecturas del modelo. Aquí los dos fixtures tienen **las mismas diez
    /// unidades completas y las mismas veinte bases cubiertas**; lo único que
    /// los separa es que en uno van seguidas y en el otro partidas en cinco y
    /// cinco.
    ///
    /// **El coeficiente va en cero**, que es su valor por omisión, y ahí está
    /// lo interesante: antes de la modulación la probabilidad no dependía de
    /// la longitud en absoluto con el coeficiente en cero. Que ahora se separen
    /// treinta y tres veces es la modulación y solo la modulación.
    ///
    /// **Debe fallar contra** una entrega que module por el total de unidades
    /// del tracto o por sus bases cubiertas: las dos darían a los dos fixtures
    /// exactamente el mismo número.
    #[test]
    fn the_modulation_takes_the_longest_run_and_not_the_total_units() {
        let settings = settings(0.0);
        let motifs = motifs("CA");
        let single_run = "CA".repeat(10);
        let split = format!("{}TT{}", "CA".repeat(5), "CA".repeat(5));
        // **Los dos totales son idénticos**, que es lo que hace impuro al
        // fixture y lo que le da poder de discriminación.
        let observed_single = crate::tandem_repeat::observe_tract(single_run.as_bytes(), &motifs);
        let observed_split = crate::tandem_repeat::observe_tract(split.as_bytes(), &motifs);
        assert_eq!(
            observed_single.covered_bases(),
            observed_split.covered_bases(),
            "el fixture solo discrimina si los dos cubren las mismas bases"
        );
        assert_ne!(
            observed_single.longest_run_units(),
            observed_split.longest_run_units(),
            "y solo discrimina si las tiradas más largas difieren"
        );
        let single = observed_slip_probability(&settings, "STR-1", single_run.as_bytes(), &motifs)
            .expect("tirada seguida");
        let partida = observed_slip_probability(&settings, "STR-1", split.as_bytes(), &motifs)
            .expect("tirada partida");
        assert!(
            single > partida * 30.0,
            "diez unidades seguidas frente a cinco y cinco: la curva las separa \
             treinta y tres veces, y {single} frente a {partida} no lo hace"
        );
    }

    /// **A la longitud de referencia se reproduce EXACTAMENTE la tasa medida.**
    ///
    /// Es la prueba de que el modelo **no cuenta dos veces** la longitud del
    /// tracto. La igualdad es exacta, bit a bit, y no aproximada: el factor
    /// vale `1.0` y multiplicar por uno no puede mover un `f64`.
    ///
    /// **Debe fallar contra** una normalización anclada en cualquier otro
    /// punto, contra una curva que se sume en vez de multiplicar, y contra una
    /// que multiplique sin normalizar.
    #[test]
    fn at_the_reference_length_the_measured_motif_rate_is_reproduced_exactly() {
        let settings = settings(0.0);
        for (motif_length, reference_bases) in [(1_u64, 15_u64), (2, 34)] {
            let units = reference_bases / motif_length;
            assert_eq!(
                super::tract_length_factor(motif_length, units),
                1.0,
                "en el ancla el factor vale exactamente uno"
            );
            let modulated = settings
                .modulated_rate_for("STR-sin-sustituir", motif_length, units)
                .expect("tasa modulada");
            assert_eq!(
                modulated,
                STUTTER_DEFAULT_MOTIF_LENGTH_RATES[motif_length as usize - 1],
                "a la longitud de referencia la tasa medida tiene que salir intacta"
            );
        }
    }

    /// **Un tracto largo tartamudea mucho más que uno corto del mismo motivo.**
    ///
    /// No basta con «más»: se exige el orden de magnitud, porque un modelo que
    /// dependiera de la longitud de forma **lineal** también daría «más» y
    /// pasaría una aserción de mera desigualdad. A lo largo del tramo medido
    /// de los mononucleótidos la curva sube más de diez mil veces, y una recta
    /// solo daría tres.
    ///
    /// **Debe fallar contra** el modelo lineal anterior y contra una tabla
    /// plana.
    #[test]
    fn a_long_tract_slips_far_more_than_a_short_one_of_the_same_motif() {
        let settings = settings(0.0);
        let short = settings
            .modulated_rate_for("STR-1", 1, 5)
            .expect("tracto corto");
        let long = settings
            .modulated_rate_for("STR-1", 1, 15)
            .expect("tracto largo");
        assert!(long > short, "el largo tiene que tartamudear más");
        assert!(
            long / short > 10_000.0,
            "de 5 a 15 pb la curva sube cuatro órdenes de magnitud, no un factor 3: {}",
            long / short
        );
        // Y el corto no es cero: si lo fuera, el cociente sería infinito y la
        // aserción de arriba pasaría sin decir nada.
        assert!(short > 0.0, "el extremo corto tiene que ser positivo");
    }

    /// **Fuera del tramo medido se mantiene el valor del extremo.**
    ///
    /// Por debajo y por encima. Es lo que la decisión llama «no extrapolar», y
    /// tiene una consecuencia que conviene ver escrita: **un tracto de 30
    /// unidades no tartamudea más que uno de 15**.
    ///
    /// **Debe fallar contra** una implementación que prolongue la recta
    /// logarítmica más allá del último punto medido, que daría al tracto de 30
    /// unidades un factor muy superior a uno.
    #[test]
    fn outside_the_measured_range_the_end_value_is_held() {
        // Por encima del ancla, que es el extremo superior.
        for units in [16_u64, 20, 30, 100] {
            assert_eq!(
                super::tract_length_factor(1, units),
                1.0,
                "{units} unidades están por encima del tramo medido: se mantiene el extremo"
            );
        }
        // Por debajo del primer punto medido.
        let floor = super::tract_length_factor(1, 5);
        for units in [1_u64, 2, 3, 4] {
            assert_eq!(
                super::tract_length_factor(1, units),
                floor,
                "{units} unidades están por debajo del tramo medido"
            );
        }
        // Y el extremo inferior **no** es el mismo número que el superior: sin
        // esto, una tabla constante satisfaría las dos aserciones anteriores.
        assert!(floor < 1.0, "los dos extremos tienen que ser distintos");
    }

    /// **Un motivo sin curva medida no tiene dependencia ninguna.**
    ///
    /// Los motivos de tres a seis pares de base no aparecen en la figura, y no
    /// se extrapola nada desde las otras dos curvas. Su factor vale
    /// exactamente uno para **cualquier** longitud de tracto.
    ///
    /// **Debe fallar contra** una implementación que reutilizara la curva de
    /// los dinucleótidos para los demás motivos, o que interpolara entre
    /// longitudes de motivo.
    #[test]
    fn a_motif_without_a_measured_curve_has_no_tract_length_dependence() {
        let settings = settings(0.0);
        for motif_length in 3_u64..=6 {
            assert_eq!(
                super::tract_length_measured_units(motif_length),
                None,
                "el motivo de {motif_length} pb no tiene curva medida"
            );
            assert_eq!(
                super::STUTTER_TRACT_LENGTH_REFERENCE_BASES[motif_length as usize - 1],
                None,
                "y por tanto tampoco tiene longitud de referencia"
            );
            let table = STUTTER_DEFAULT_MOTIF_LENGTH_RATES[motif_length as usize - 1];
            for units in [1_u64, 5, 12, 40] {
                assert_eq!(super::tract_length_factor(motif_length, units), 1.0);
                assert_eq!(
                    settings
                        .modulated_rate_for("STR-1", motif_length, units)
                        .expect("tasa"),
                    table,
                    "la tasa de un motivo sin curva no cambia con la longitud"
                );
            }
        }
        // Control: los dos motivos que **sí** tienen curva no salen ilesos de
        // la misma comprobación. Sin él, un `tract_length_factor` que
        // devolviera siempre uno pasaría el bucle de arriba.
        for motif_length in 1_u64..=2 {
            assert!(super::tract_length_factor(motif_length, 5) < 1.0);
        }
    }

    /// **La sustitución del autor del escenario gana, y no se modula.**
    ///
    /// Quien declara la tasa de un locus la ha elegido para ese locus, con su
    /// longitud dentro. Modularla encima contaría dos veces lo mismo.
    ///
    /// **Debe fallar contra** una implementación que aplicara la curva también
    /// a la tasa sustituida: con cinco unidades de dinucleótido el factor vale
    /// 0,0038 y la tasa declarada se quedaría en 0,0019 en vez de en 0,5.
    #[test]
    fn a_declared_locus_rate_prevails_and_is_not_modulated() {
        let mut declared = std::collections::BTreeMap::new();
        declared.insert("STR-SUYO".to_owned(), 0.5);
        let settings = StutterSettings::new(
            true,
            StutterMotifLengthRates::default(),
            0.0,
            Some(StutterStepWeights::new([0.05, 0.15, 0.6, 0.2]).expect("pesos")),
            declared,
        )
        .expect("parámetros");
        for units in [1_u64, 5, 13, 17, 40] {
            assert_eq!(
                settings
                    .modulated_rate_for("STR-SUYO", 2, units)
                    .expect("tasa declarada"),
                0.5,
                "la sustitución vale lo mismo mida lo que mida el tracto"
            );
        }
        // Control positivo: el locus **no** sustituido sí se modula, y con
        // cinco unidades queda muy por debajo de su tasa de tabla. Sin este
        // control, una implementación que no modulara nada pasaría el bucle.
        let table = STUTTER_DEFAULT_MOTIF_LENGTH_RATES[1];
        let other = settings
            .modulated_rate_for("STR-AJENO", 2, 5)
            .expect("tasa de tabla");
        assert!(
            other < table / 100.0,
            "el locus sin sustituir sí tiene que modularse: {other} frente a {table}"
        );
    }

    /// El paso se aplica **en la tirada más larga**, por su extremo derecho, y
    /// deja intacto todo lo demás de la molécula.
    #[test]
    fn the_step_is_applied_at_the_longest_run_and_leaves_the_rest_untouched() {
        let mut sequence = format!("GGG{}TT{}GGG", "CA".repeat(5), "CA".repeat(3)).into_bytes();
        // La tirada más larga es la de cinco, que empieza en el desplazamiento 3.
        apply_step_to_molecule(&mut sequence, 3, 5, 2, -1).expect("contracción");
        assert_eq!(
            String::from_utf8(sequence).expect("texto"),
            format!("GGG{}TT{}GGG", "CA".repeat(4), "CA".repeat(3))
        );
    }
}
