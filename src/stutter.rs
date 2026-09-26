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
pub const STUTTER_ALGORITHM: &str = "stutter-v1";

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
pub const STUTTER_CALIBRATION_SOURCE: &str =
    "Raz et al. (2019), Nucleic Acids Research 47(5):2436-2445, CC BY 4.0";
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

/// **Probabilidad de resbalón por ciclo**, acotada a `[0,1]`.
///
/// `base_rate` es la fracción calibrada de la longitud de motivo de la tirada;
/// `coefficient` es el coeficiente de la tirada (cero mientras no esté
/// calibrado) y `longest_run_units` es el predictor.
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
    #[must_use]
    pub fn base_rate_for(&self, locus_id: &str, motif_length: u64) -> Option<f64> {
        self.locus_rates.get(locus_id).copied().map_or_else(
            || self.motif_length_rates.rate_for_motif_length(motif_length),
            Some,
        )
    }
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
/// # Una limitación conocida con las quimeras, dicha en voz alta
///
/// El tracto se sitúa traduciendo el intervalo del locus a coordenadas de la
/// copia y recortándolo contra el intervalo de la molécula. En una molécula
/// **quimérica** ese intervalo deja de describir un tramo contiguo de la copia
/// más allá del punto de rotura, y la procedencia de la segunda mitad es
/// transitoria: la quimera no la conserva. En consecuencia, sobre una molécula
/// quimérica la ventana puede caer donde el tracto ya no está.
///
/// **Lo que eso produce, y lo que no.** La ventana se vuelve a leer con
/// [`observe_tract`] antes de decidir nada: si dentro no hay ninguna copia
/// completa del motivo, no hay tirada, no se evalúa y **no se consume
/// aleatoriedad**. Nunca se escribe un paso sobre bases que no repiten. Lo que
/// sí puede pasar es que una molécula quimérica que sigue conteniendo el tracto
/// no llegue a tartamudear. Con la tasa de quimeras del perfil vigente, cero,
/// el caso no existe.
pub fn transform_paired_end(
    root: RootSeed,
    event_id: &EventId,
    individual: &Individual,
    plan: &StutterPlan,
    lineage: Option<&crate::clonal_stutter::LineageNetSteps>,
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
    let placements = place_tracts(
        individual,
        plan,
        &copy_id,
        pair.fragment().start(),
        original_length,
    )?;

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
        let net = lineage.map_or(0, |steps| steps.step_for(locus.id().as_str()));
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
        let Some(base_rate) = plan
            .settings()
            .base_rate_for(locus.id().as_str(), motif_length)
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
    let base_rate = settings.base_rate_for(locus_id, observed.longest_run_motif_length())?;
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
    let Some(copy) = individual
        .copies_for(copy_id.sequence_id())
        .and_then(|copies| copies.iter().find(|candidate| *candidate.id() == *copy_id))
    else {
        return Err(StutterError::CopyNotFound {
            sequence_id: copy_id.sequence_id().as_str().to_owned(),
            copy_index: copy_id.copy_index(),
        });
    };
    let mut placements = Vec::new();
    for locus in plan.loci() {
        if locus.sequence_id() != copy_id.sequence_id() {
            continue;
        }
        let Some(placement) = locate_tract_in_copy(copy, locus)? else {
            continue;
        };
        let Some((start, end)) = clip_to_molecule(placement, fragment_start, fragment_length)
        else {
            continue;
        };
        placements.push(PlacedTract { locus, start, end });
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
