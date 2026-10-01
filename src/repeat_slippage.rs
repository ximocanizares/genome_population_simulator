//! F6.4.4 y F6.4.11: **deslizamiento biológico** de un microsatélite.
//!
//! # Qué decide este módulo, y qué no
//!
//! Hasta aquí los alelos de un locus de repetición y sus frecuencias los
//! **escribía el usuario a mano**. Este módulo los **genera**: parte del alelo
//! que el genoma de referencia lleva en ese locus, aplica pasos de una unidad
//! completa y produce el abanico de alelos con sus frecuencias, que después el
//! sorteo de `repeat_genotypes` reparte entre las copias haplotípicas.
//!
//! Es **biología, no artefacto**, y por eso va en la **verdad**. No se
//! confunde con el tartamudeo de `stutter`, que comparte el mecanismo físico
//! pero ensucia **la medida** y vive en las lecturas. La prueba de que no son
//! intercambiables: si los alelos verdaderos salieran del modelo de
//! tartamudeo, cambiar el número de ciclos de PCR cambiaría el genotipo real
//! de los individuos, es decir, un parámetro de laboratorio decidiendo la
//! biología.
//!
//! # Es **una sola tirada de inicialización**
//!
//! No es una simulación de generaciones. Cada copia haplotípica de la
//! población recibe un número de pasos y acumula su deriva de una vez. El
//! modelo escalonado aplicado generación tras generación pertenece a la fase
//! de evolución y este trabajo no la invade.
//!
//! # El algoritmo, adoptado por el responsable el 20 de septiembre de 2026
//!
//! Se simulan **tantos linajes independientes como copias haplotípicas** tenga
//! la población. Cada linaje parte del alelo de la referencia y recibe un
//! **número de pasos sorteado de una Poisson de media declarada**, con
//! `poisson-count-v1`, que es el único recuento normativo del proyecto; cada
//! paso lo da el servicio escalonado de este mismo módulo con su sesgo. La
//! **frecuencia de cada alelo es la proporción de linajes que terminan en él**.
//!
//! El alelo de la referencia entra en el abanico **si algún linaje termina en
//! él**, y su presencia **no se fuerza**: que desaparezca depende de los
//! parámetros que se den, no del modelo. Con una media de pasos baja la mayoría
//! de los linajes no se mueven y el alelo está; con una media alta puede
//! perderse, que es lo que ocurre en una población real.
//!
//! # La regla de agotamiento, igual que la del tartamudeo y por la misma razón
//!
//! Un paso que dejaría **menos de una unidad completa en la tirada más larga**
//! no es viable: se excluye **antes** de sortear y **no se reintenta**. Si no
//! queda ningún paso viable —o ninguno con peso positivo—, el evento **no se
//! aplica**, **no consume aleatoriedad** y se **cuenta aparte** como evento no
//! aplicado. No es un error: un microsatélite que ha llegado a su mínimo es un
//! resultado biológico legítimo. Un linaje que rebasaría el tope de unidades
//! del locus se detiene ahí por la misma vía.
//!
//! # Las tres limitaciones declaradas, que son biología y no detalles
//!
//! Van en [`SLIPPAGE_LIMITATIONS`], se publican en el manifiesto y se repiten
//! aquí porque quien lea el código tiene que verlas antes que nadie:
//!
//! 1. **No hay ancestría compartida.** Cada copia deriva por su cuenta desde
//!    el ancestro, de modo que dos individuos emparentados no comparten
//!    historia en este locus. El espectro de frecuencias sale **más plano** que
//!    el de una población real, donde la deriva concentra unos pocos alelos.
//!    Importa si alguien usa el conjunto para probar herramientas de estructura
//!    poblacional.
//! 2. **El número de pasos no es un proceso en el tiempo.** Cada copia acumula
//!    mutaciones sin decir en cuántas generaciones.
//! 3. **El genoma de referencia puede llevar un alelo que ningún individuo
//!    tiene.** Ocurre de verdad, pero conviene saberlo antes de encadenar un
//!    análisis que lo dé por imposible.
//!
//! # El ámbito de azar, y por qué no puede ser el del tartamudeo
//!
//! El deslizamiento biológico y el tartamudeo **no pueden compartir ámbito de
//! azar**, porque uno decide la verdad y el otro las lecturas. El ámbito de
//! este módulo se identifica por el **evento** —que es el nombre del propio
//! ámbito, `initial_repeat_slippage`— y por el **identificador de locus**,
//! nunca por la posición del tracto ni por el orden en que los loci se
//! declararon. La consecuencia observable es la que importa: **insertar un
//! locus que ordene antes que los existentes no cambia el abanico de los
//! previos**.
//!
//! **Lo que *no* demuestra la separación**, y conviene decirlo porque parece
//! que sí: «cambiar los ciclos de PCR no cambia ni un alelo de la verdad» se
//! cumple sola. DTC-081 exige que cada decisión arranque su propio generador
//! desde su clave y que no haya consumo acumulado, de modo que un consumidor no
//! desplaza al otro; y además la verdad se genera **antes** que las lecturas.
//! Lo que sí caza un ámbito compartido es **comparar las dos subsemillas**.

use std::collections::BTreeMap;

use crate::diagnostics::Diagnostic;
use crate::randomness::{
    F4_BERNOULLI_ALGORITHM, POISSON_COUNT_ALGORITHM, RandomSamplingError, RandomWordSource,
    RootSeed, RuntimeSubseedScope, SUBSEED_CATALOG_VERSION, derive_subseed,
    f4_bernoulli_with_source, neumaier_sum, poisson_count_with_source,
};
use crate::tandem_repeat::{TandemRepeatError, TandemRepeatTract};

/// Ámbito de subsemilla del generador, tal y como lo publica el catálogo.
pub const SLIPPAGE_SCOPE: &str = "initial_repeat_slippage";

/// Versión del catálogo de ámbitos en la que nace este ámbito.
pub const SLIPPAGE_SCOPE_CATALOG_VERSION: &str = SUBSEED_CATALOG_VERSION;

/// Escenario de los ámbitos de azar de la población inicial, el mismo que usan
/// la composición genotípica, la fase y el sorteo del alelo de repetición.
pub const SLIPPAGE_SCENARIO: &str = "direct";

/// Versión del algoritmo del generador. Se publica con cada abanico.
pub const SLIPPAGE_ALGORITHM: &str = "biological-slippage-v1";

/// El único modelo admitido. Cualquier otro se rechaza nombrando parámetro y
/// valor, antes de simular.
pub const SLIPPAGE_MODEL_STEPWISE: &str = "stepwise";

/// El único tamaño de paso admitido, en unidades completas. DTF-015 mantiene
/// fuera el salto de varias unidades.
pub const SLIPPAGE_STEP_UNITS: u64 = 1;

/// Los dos pasos del modelo, en orden canónico.
pub const SLIPPAGE_ADMITTED_STEPS: [i64; 2] = [-1, 1];

/// Sesgo de expansión por omisión: simétrico, ni expansivo ni contractivo.
pub const SLIPPAGE_DEFAULT_EXPANSION_BIAS: f64 = 0.5;

/// Techo de la media de pasos por linaje.
///
/// # Por qué hay techo, y por qué este
///
/// `poisson-count-v1` admite medias hasta `2^52`, que es su dominio como
/// primitiva. Aquí el parámetro lo **escribe el usuario**, y sin techo se van de
/// las manos dos cosas: el **coste**, porque cada paso clona el tracto y el
/// trabajo por linaje crece con la media; y el **sentido**, porque un tracto
/// tiene a lo sumo `max_units` unidades y un linaje se detiene en cuanto toca un
/// borde, de modo que por encima de unos pocos cientos de pasos la media deja de
/// cambiar el resultado y solo cambia lo que se tarda.
///
/// Mil está holgadamente por encima de cualquier tracto que este simulador
/// admita, así que no recorta ningún escenario con sentido biológico, y acota el
/// trabajo. **Solo rechaza; no cambia el resultado de ningún escenario
/// aceptado.**
pub const SLIPPAGE_MAX_MEAN_EVENTS: f64 = 1000.0;

/// Recuento normativo del número de pasos por linaje (DTC-067).
pub const SLIPPAGE_COUNT_ALGORITHM: &str = POISSON_COUNT_ALGORITHM;

/// Ensayo binario que elige entre expandir y contraer cuando los dos pasos son
/// viables. Las probabilidades exactas 0 y 1 **no consumen aleatoriedad**.
pub const SLIPPAGE_BIAS_ALGORITHM: &str = F4_BERNOULLI_ALGORITHM;

/// Las tres limitaciones biológicas del modelo, tal y como se publican.
///
/// Van en el manifiesto **palabra por palabra**. Son decisiones del modelo y no
/// detalles de implementación: quien encadene un análisis sobre estos datos
/// tiene que leerlas sin abrir el código.
///
/// # Su texto está anclado FUERA de este módulo, y hace falta que lo esté
///
/// El juez del manifiesto compara el texto publicado contra **esta misma
/// constante**, de modo que por sí solo es un espejo: cambiar aquí lo que el
/// simulador le cuenta al usuario sobre las limitaciones de su propio modelo
/// pasaría sin que nada fallara. Se midió, y pasaba.
///
/// Lo que lo impide es `tests/repeat_slippage.rs`, que lleva los tres textos
/// **escritos a mano** y los compara contra los **bytes publicados**. Si se
/// toca esta constante hay que tocar aquella prueba, y entonces el cambio deja
/// de ser invisible. La tercera copia, en castellano, vive en DTC-095.
pub const SLIPPAGE_LIMITATIONS: [&str; 3] = [
    "no shared ancestry: every haplotype copy derives independently from the reference allele, so \
related individuals share no history at this locus and the resulting frequency spectrum is \
flatter than a real population's, where drift concentrates a few alleles; it matters if the \
dataset is used to test population-structure tools",
    "the step count is not a process in time: each copy accumulates a random number of mutations \
without stating over how many generations, which is consistent with this being a single \
initialisation draw and not a simulation of generations",
    "the reference genome may carry an allele that no individual carries: it does happen in \
reality, but it must be known before chaining an analysis that rules it out",
];

/// Parámetros del generador, ya validados.
///
/// Se construye por [`SlippageSettings::new`], que es quien rechaza un modelo
/// distinto del admitido y un salto de varias unidades **nombrando el parámetro
/// y su valor**.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlippageSettings {
    mean_events: f64,
    expansion_bias: f64,
}

impl SlippageSettings {
    /// Valida el modelo, el tamaño de paso, la media de pasos y el sesgo.
    ///
    /// El rechazo llega **antes de simular** y nombra el parámetro y el valor,
    /// que es lo que hace accionable el diagnóstico.
    pub fn new(
        model: &str,
        step_units: u64,
        mean_events: f64,
        expansion_bias: f64,
    ) -> Result<Self, SlippageError> {
        if model != SLIPPAGE_MODEL_STEPWISE {
            return Err(SlippageError::UnknownModel {
                value: model.to_owned(),
            });
        }
        if step_units != SLIPPAGE_STEP_UNITS {
            return Err(SlippageError::UnsupportedStepUnits { value: step_units });
        }
        if !mean_events.is_finite() || !(0.0..=SLIPPAGE_MAX_MEAN_EVENTS).contains(&mean_events) {
            return Err(SlippageError::InvalidMeanEvents { value: mean_events });
        }
        if !expansion_bias.is_finite() || !(0.0..=1.0).contains(&expansion_bias) {
            return Err(SlippageError::InvalidExpansionBias {
                value: expansion_bias,
            });
        }
        Ok(Self {
            mean_events,
            expansion_bias,
        })
    }

    /// Media de la Poisson que decide cuántos pasos recibe cada linaje.
    #[must_use]
    pub const fn mean_events(self) -> f64 {
        self.mean_events
    }

    /// Probabilidad de que un paso viable sea una expansión.
    #[must_use]
    pub const fn expansion_bias(self) -> f64 {
        self.expansion_bias
    }

    /// Peso del paso, que es lo que decide si un paso viable llega a aplicarse.
    ///
    /// Vive aparte porque es la pieza que hace que **el recuento de eventos no
    /// aplicados** signifique algo: un paso viable con peso cero no se aplica y
    /// no consume aleatoriedad, igual que uno inviable.
    #[must_use]
    fn weight_of(self, step: i64) -> f64 {
        if step > 0 {
            self.expansion_bias
        } else {
            1.0 - self.expansion_bias
        }
    }
}

/// Lo que un evento hizo, que es lo que permite contarlos por separado.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SlippageEvent {
    /// El evento aplicó un paso de una unidad completa.
    Applied(i64),
    /// El linaje había llegado a un **borde** y se detiene ahí. **No consumió
    /// aleatoriedad** y se cuenta aparte.
    Stopped,
}

/// El resultado de aplicar una tanda de eventos a un tracto.
#[derive(Clone, Debug, PartialEq)]
pub struct SlippageRun {
    tract: TandemRepeatTract,
    applied_events: u64,
    unapplied_events: u64,
}

impl SlippageRun {
    #[must_use]
    pub const fn tract(&self) -> &TandemRepeatTract {
        &self.tract
    }

    #[must_use]
    pub const fn applied_events(&self) -> u64 {
        self.applied_events
    }

    /// Eventos que **no se aplicaron** al agotar el tracto o al topar.
    #[must_use]
    pub const fn unapplied_events(&self) -> u64 {
        self.unapplied_events
    }

    #[must_use]
    pub fn into_tract(self) -> TandemRepeatTract {
        self.tract
    }
}

/// Un alelo del abanico generado.
#[derive(Clone, Debug, PartialEq)]
pub struct SlippageAllele {
    tract: TandemRepeatTract,
    lineages: u64,
    frequency: f64,
}

impl SlippageAllele {
    #[must_use]
    pub const fn tract(&self) -> &TandemRepeatTract {
        &self.tract
    }

    /// Linajes que terminaron en este alelo: el **numerador** de la frecuencia.
    #[must_use]
    pub const fn lineages(&self) -> u64 {
        self.lineages
    }

    /// Proporción de linajes que terminaron en este alelo.
    #[must_use]
    pub const fn frequency(&self) -> f64 {
        self.frequency
    }

    #[must_use]
    pub fn into_parts(self) -> (TandemRepeatTract, u64, f64) {
        (self.tract, self.lineages, self.frequency)
    }
}

/// El abanico completo de un locus, con la contabilidad de su generación.
#[derive(Clone, Debug, PartialEq)]
pub struct SlippageFan {
    alleles: Vec<SlippageAllele>,
    lineages: u64,
    drawn_events: u64,
    applied_events: u64,
    unapplied_events: u64,
}

impl SlippageFan {
    /// Los alelos en orden canónico: **por longitud del tracto y, a igual
    /// longitud, por sus bases**.
    ///
    /// El orden no depende de en qué linaje apareció cada alelo por primera
    /// vez, que sería un orden de recorrido y no una propiedad del abanico.
    #[must_use]
    pub fn alleles(&self) -> &[SlippageAllele] {
        &self.alleles
    }

    /// Linajes simulados, que son las copias haplotípicas de la población.
    #[must_use]
    pub const fn lineages(&self) -> u64 {
        self.lineages
    }

    /// Pasos **sorteados** por la Poisson, sumados sobre todos los linajes.
    #[must_use]
    pub const fn drawn_events(&self) -> u64 {
        self.drawn_events
    }

    /// Pasos que llegaron a cambiar el tracto.
    #[must_use]
    pub const fn applied_events(&self) -> u64 {
        self.applied_events
    }

    /// Pasos que **no se aplicaron** al agotar el tracto o al topar.
    ///
    /// Es la cifra que distingue un servicio que ignora en silencio los eventos
    /// sobrantes de uno que los contabiliza. Sin ella los dos producirían el
    /// mismo tracto final y nadie notaría la diferencia.
    #[must_use]
    pub const fn unapplied_events(&self) -> u64 {
        self.unapplied_events
    }

    #[must_use]
    pub fn into_alleles(self) -> Vec<SlippageAllele> {
        self.alleles
    }
}

/// Los pasos que sobreviven al filtrado previo, en el orden canónico.
///
/// Un paso entra si deja **al menos una unidad completa en la tirada más
/// larga** —la misma regla que `TandemRepeatTract::apply_unit_step` y que el
/// tartamudeo— y si no rebasa el tope de unidades del locus.
///
/// La contracción mira la **tirada más larga** y no el tracto entero, igual que
/// el paso del dominio: es donde la unidad se gana o se pierde.
///
/// # Que falte uno significa que el linaje se detiene, no que se le empuje
///
/// Cuando esta lista **no contiene los dos pasos admitidos**, el linaje ha
/// llegado a un borde —el suelo de una unidad completa, o el tope del locus— y
/// [`apply_slippage_events_with_source`] lo detiene ahí. **No se aplica el paso
/// que queda.**
///
/// Es una decisión del responsable, tomada el 27 de septiembre de 2026 con las
/// dos opciones delante, y descartó expresamente la otra: dejar que el borde
/// empujara. Con un solo candidato viable el servicio lo aplicaba de forma
/// **determinista**, de modo que un linaje que tocaba el tope se iba contrayendo
/// paso a paso, todos los eventos se contaban como aplicados, ninguno como no
/// aplicado, y el documento publicaba un sesgo intermedio cuando no había habido
/// ni una expansión. El borde decidía la biología en lugar del sesgo.
#[must_use]
pub fn viable_unit_steps(tract: &TandemRepeatTract, max_units: u64) -> Vec<i64> {
    let annotation = tract.annotation();
    let longest_run_units = annotation.longest_run_units();
    let units = annotation.units_as_f64();
    SLIPPAGE_ADMITTED_STEPS
        .into_iter()
        .filter(|step| {
            if *step > 0 {
                units + 1.0 <= max_units as f64
            } else {
                // Una contracción de una unidad exige que la tirada más larga
                // tenga al menos dos: quitar la única copia la vaciaría.
                longest_run_units >= 2
            }
        })
        .collect()
}

/// Elige el paso entre los candidatos **ya filtrados**, o no elige ninguno.
///
/// Vive aparte del recorrido para que una prueba pueda **contar** las palabras
/// que consume: en los extremos del sesgo tienen que ser cero, y sin poder
/// contarlas esa afirmación no sería comprobable. Es el mismo motivo por el que
/// el tartamudeo y el sorteo del alelo tienen su propia versión instrumentable.
///
/// Solo se llama con **los dos pasos admitidos** viables: en un borde el linaje
/// se detiene y no llega aquí. Con los dos viables, **alguno tiene peso
/// positivo** —los dos pesos suman uno y el sesgo se validó finito en `[0,1]`—,
/// de modo que siempre hay un paso que elegir.
///
/// Con el sesgo en un extremo solo queda un candidato con peso y **no se
/// sortea**: aplicarlo no puede desplazar la secuencia de azar de nada.
pub(crate) fn choose_unit_step_with_source(
    rng: &mut impl RandomWordSource,
    settings: SlippageSettings,
    candidates: &[i64],
) -> Result<i64, SlippageError> {
    let positive: Vec<i64> = candidates
        .iter()
        .copied()
        .filter(|step| settings.weight_of(*step) > 0.0)
        .collect();
    match positive.as_slice() {
        // Los dos pesos suman uno: que los dos sean cero es imposible con un
        // sesgo finito en `[0,1]`, que es el único que la validación admite.
        [] => Err(SlippageError::NoWeightedStep),
        [only] => Ok(*only),
        _ => {
            let expands = f4_bernoulli_with_source(rng, settings.expansion_bias())
                .map_err(SlippageError::Random)?;
            Ok(if expands { 1 } else { -1 })
        }
    }
}

/// Aplica **un** evento al tracto, con la regla del borde completa.
///
/// Si el linaje está en un borde —le falta alguno de los dos pasos admitidos—,
/// **se detiene**: el tracto no cambia, no se consume aleatoriedad y el evento
/// se cuenta aparte.
pub(crate) fn apply_slippage_event_with_source(
    rng: &mut impl RandomWordSource,
    tract: &TandemRepeatTract,
    max_units: u64,
    settings: SlippageSettings,
) -> Result<(TandemRepeatTract, SlippageEvent), SlippageError> {
    let candidates = viable_unit_steps(tract, max_units);
    if candidates.len() != SLIPPAGE_ADMITTED_STEPS.len() {
        return Ok((tract.clone(), SlippageEvent::Stopped));
    }
    let step = choose_unit_step_with_source(rng, settings, &candidates)?;
    let stepped = tract.apply_unit_step(step).map_err(SlippageError::Tract)?;
    Ok((stepped, SlippageEvent::Applied(step)))
}

/// Aplica una tanda de eventos y devuelve el tracto final con su contabilidad.
pub(crate) fn apply_slippage_events_with_source(
    rng: &mut impl RandomWordSource,
    tract: &TandemRepeatTract,
    max_units: u64,
    settings: SlippageSettings,
    events: u64,
) -> Result<SlippageRun, SlippageError> {
    let mut current = tract.clone();
    let mut applied_events = 0_u64;
    let mut unapplied_events = 0_u64;
    // **Una vez detenido, el linaje no vuelve.** No es un atajo: detenerse y
    // volver a mirar en el evento siguiente serían modelos distintos, porque el
    // tracto no ha cambiado y la respuesta sería la misma. Se escribe así para
    // que se vea que el borde detiene el linaje y no solo el evento.
    let mut stopped = false;
    for _ in 0..events {
        if stopped {
            unapplied_events += 1;
            continue;
        }
        let (next, outcome) = apply_slippage_event_with_source(rng, &current, max_units, settings)?;
        current = next;
        match outcome {
            SlippageEvent::Applied(_) => applied_events += 1,
            SlippageEvent::Stopped => {
                stopped = true;
                unapplied_events += 1;
            }
        }
    }
    Ok(SlippageRun {
        tract: current,
        applied_events,
        unapplied_events,
    })
}

/// El servicio escalonado completo contra un generador local, para quien no
/// necesite instrumentar el consumo.
pub fn apply_slippage_events(
    rng: &mut crate::randomness::LocalRng,
    tract: &TandemRepeatTract,
    max_units: u64,
    settings: SlippageSettings,
    events: u64,
) -> Result<SlippageRun, SlippageError> {
    apply_slippage_events_with_source(rng, tract, max_units, settings, events)
}

/// Genera el abanico de alelos de un locus a partir del alelo de la referencia.
///
/// El ámbito de azar se deriva del escenario, la referencia y el
/// **identificador del locus**. Ni la posición del tracto ni el orden de
/// declaración entran, de modo que insertar un locus que ordene antes que los
/// existentes no cambia el abanico de los previos.
pub fn generate_allele_fan(
    root_seed: RootSeed,
    reference_id: &str,
    locus_id: &str,
    ancestral: &TandemRepeatTract,
    max_units: u64,
    settings: SlippageSettings,
    lineages: u64,
) -> Result<SlippageFan, SlippageError> {
    let mut rng = slippage_subseed(root_seed, reference_id, locus_id)?.local_rng();
    generate_allele_fan_with_source(&mut rng, ancestral, max_units, settings, lineages)
}

/// La subsemilla **que el generador usa de verdad**.
///
/// Vive expuesta a propósito: es lo único que permite comprobar desde fuera que
/// el deslizamiento biológico **no comparte ámbito con el tartamudeo**.
/// Comparar las dos subsemillas es la prueba de la separación; «cambiar los
/// ciclos de PCR no cambia ni un alelo de la verdad» **no lo es**, porque se
/// cumple sola.
pub fn slippage_subseed(
    root_seed: RootSeed,
    reference_id: &str,
    locus_id: &str,
) -> Result<crate::randomness::Subseed, SlippageError> {
    derive_subseed(
        root_seed,
        RuntimeSubseedScope::InitialRepeatSlippage {
            scenario_id: SLIPPAGE_SCENARIO,
            reference_id,
            locus_id,
        },
    )
    .map_err(|diagnostic| SlippageError::Subseed {
        detail: format!("{diagnostic:?}"),
    })
}

/// El mismo generador contra una fuente de palabras cualquiera.
///
/// Existe para que una prueba pueda **contar** cuántas palabras consume: con
/// una media de cero pasos tienen que ser cero, y sin poder contarlas esa
/// afirmación no sería comprobable.
pub(crate) fn generate_allele_fan_with_source(
    rng: &mut impl RandomWordSource,
    ancestral: &TandemRepeatTract,
    max_units: u64,
    settings: SlippageSettings,
    lineages: u64,
) -> Result<SlippageFan, SlippageError> {
    if lineages == 0 {
        return Err(SlippageError::NoLineages);
    }
    let annotation = ancestral.annotation();
    if annotation.longest_run_units() == 0 {
        return Err(SlippageError::AncestralWithoutCompleteUnit);
    }
    // **El tope se comprueba sobre el ancestro, no solo sobre los pasos.**
    // Mirándolo solo al expandir, el tope se convertía en un «no crezcas más» y
    // dejaba pasar un tramo de referencia que ya lo rebasaba: el mismo alelo
    // declarado a mano se rechazaba y generado se admitía. Y el generador se
    // volvía **unidireccional**, porque desde encima del tope la expansión
    // nunca es viable y el sesgo declarado dejaba de significar nada.
    let units = annotation.units_as_f64();
    if units > max_units as f64 {
        return Err(SlippageError::AncestralAboveCap { units, max_units });
    }
    // Clave canónica: longitud del tracto y, a igual longitud, sus bases. No es
    // el orden de aparición, que sería un orden de recorrido.
    let mut counts: BTreeMap<(usize, Vec<u8>), (TandemRepeatTract, u64)> = BTreeMap::new();
    let mut drawn_events = 0_u64;
    let mut applied_events = 0_u64;
    let mut unapplied_events = 0_u64;
    for _ in 0..lineages {
        let events = poisson_count_with_source(rng, settings.mean_events(), None)
            .map_err(SlippageError::Random)?;
        drawn_events = drawn_events
            .checked_add(events)
            .ok_or(SlippageError::EventCountOverflow)?;
        let run = apply_slippage_events_with_source(rng, ancestral, max_units, settings, events)?;
        applied_events = applied_events
            .checked_add(run.applied_events())
            .ok_or(SlippageError::EventCountOverflow)?;
        unapplied_events = unapplied_events
            .checked_add(run.unapplied_events())
            .ok_or(SlippageError::EventCountOverflow)?;
        let tract = run.into_tract();
        let key = (tract.bases().len(), tract.bases().to_vec());
        counts
            .entry(key)
            .and_modify(|(_, count)| *count += 1)
            .or_insert((tract, 1));
    }
    let alleles: Vec<SlippageAllele> = counts
        .into_values()
        .map(|(tract, count)| SlippageAllele {
            tract,
            lineages: count,
            frequency: count as f64 / lineages as f64,
        })
        .collect();
    // La suma se recomputa con el mismo algoritmo compensado que consume el
    // sorteo categórico (DTC-049), en vez de darla por uno porque los
    // numeradores sumen el denominador.
    let sum = neumaier_sum(alleles.iter().map(SlippageAllele::frequency));
    if !sum.is_finite() || sum <= 0.0 {
        return Err(SlippageError::DegenerateFrequencies { sum });
    }
    Ok(SlippageFan {
        alleles,
        lineages,
        drawn_events,
        applied_events,
        unapplied_events,
    })
}

/// Los rechazos del generador. Todos llegan **antes de simular** salvo los dos
/// numéricos, que son invariantes internas.
#[derive(Clone, Debug, PartialEq)]
pub enum SlippageError {
    /// Un modelo distinto del escalonado. DTF-015 mantiene fuera los demás.
    UnknownModel {
        value: String,
    },
    /// Un salto de más de una unidad completa.
    UnsupportedStepUnits {
        value: u64,
    },
    InvalidMeanEvents {
        value: f64,
    },
    InvalidExpansionBias {
        value: f64,
    },
    /// El generador necesita al menos una copia haplotípica que simular.
    NoLineages,
    /// El tramo de la referencia no contiene ni una copia completa del motivo,
    /// de modo que no hay alelo ancestral del que partir.
    ///
    /// **No se alcanza desde una configuración, y se dice en vez de fingir que
    /// sí.** Con referencia suministrada, la guarda de admisión de F6.4.2 exige
    /// una tirada mínima —tres unidades por omisión— y caza antes el tramo que
    /// no es un microsatélite. Con referencia simulada, desde el 01/10/2026, el
    /// simulador coloca el tracto y la misma guarda lo juzga después (S-06
    /// corregida, DT-011). Se conserva
    /// porque es el contrato del generador y no puede depender de que otra
    /// guarda siga existiendo, y lo ejercita la prueba del módulo, que llama al
    /// generador directamente. Anotado como DT-011.
    AncestralWithoutCompleteUnit,
    /// El tramo de la referencia ya rebasa el tope de unidades del locus.
    AncestralAboveCap {
        units: f64,
        max_units: u64,
    },
    /// Invariante interna: con los dos pasos viables, los dos pesos suman uno y
    /// alguno tiene que ser positivo.
    NoWeightedStep,
    /// El generador solo admite un locus de **un** motivo: sin fronteras de
    /// bloque declaradas no hay forma de partir el tracto de la referencia.
    CompositeLocus {
        motifs: usize,
    },
    EventCountOverflow,
    /// **No se alcanza desde una configuración**, y se dice: las frecuencias son
    /// proporciones de linajes con denominador positivo, de modo que su suma
    /// solo podría no ser finita si el recuento desbordara, y eso lo caza antes
    /// `EventCountOverflow`. Se conserva como red de la invariante de DTC-049.
    /// Anotado como DT-011.
    DegenerateFrequencies {
        sum: f64,
    },
    Subseed {
        detail: String,
    },
    Random(RandomSamplingError),
    Tract(TandemRepeatError),
}

impl SlippageError {
    /// El parámetro al que apunta el rechazo, para que el diagnóstico sea
    /// accionable sin leer el código.
    #[must_use]
    pub const fn parameter(&self) -> &'static str {
        match self {
            Self::UnknownModel { .. } => "reference.repeat_loci.slippage.model",
            Self::UnsupportedStepUnits { .. } => "reference.repeat_loci.slippage.step_units",
            Self::InvalidMeanEvents { .. } => "reference.repeat_loci.slippage.mean_events",
            Self::InvalidExpansionBias { .. } => "reference.repeat_loci.slippage.expansion_bias",
            Self::NoLineages => "initial_populations",
            Self::AncestralWithoutCompleteUnit | Self::CompositeLocus { .. } => {
                "reference.repeat_loci"
            }
            Self::AncestralAboveCap { .. } => "reference.repeat_loci.max_units",
            Self::NoWeightedStep => "reference.repeat_loci.slippage.expansion_bias",
            Self::EventCountOverflow
            | Self::DegenerateFrequencies { .. }
            | Self::Subseed { .. }
            | Self::Random(_)
            | Self::Tract(_) => "reference.repeat_loci.slippage",
        }
    }

    /// El detalle, que nombra el valor rechazado siempre que haya uno.
    #[must_use]
    pub fn detail(&self) -> String {
        match self {
            Self::UnknownModel { value } => format!(
                "el modelo «{value}» no existe: el único admitido es «{SLIPPAGE_MODEL_STEPWISE}»"
            ),
            Self::UnsupportedStepUnits { value } => format!(
                "un paso de {value} unidades no está admitido: el modelo escalonado da pasos de \
{SLIPPAGE_STEP_UNITS} unidad completa"
            ),
            Self::InvalidMeanEvents { value } => {
                format!("la media de pasos {value} está fuera de [0, {SLIPPAGE_MAX_MEAN_EVENTS}]")
            }
            Self::InvalidExpansionBias { value } => {
                format!("el sesgo de expansión {value} está fuera de [0, 1]")
            }
            Self::NoLineages => {
                "el generador simula un linaje por copia haplotípica y el escenario no declara \
ninguna"
                    .to_owned()
            }
            Self::AncestralWithoutCompleteUnit => {
                "el tramo de la referencia no contiene ni una copia completa del motivo, de modo \
que no hay alelo ancestral del que partir"
                    .to_owned()
            }
            Self::CompositeLocus { motifs } => format!(
                "el locus declara {motifs} motivos y el generador solo admite uno: sin fronteras \
de bloque declaradas no hay forma de partir el tracto de la referencia"
            ),
            Self::AncestralAboveCap { units, max_units } => format!(
                "el tramo de la referencia mide {units} unidades y rebasa el tope de {max_units} \
del locus; el generador parte de ese tramo, de modo que el tope tiene que caberlo"
            ),
            Self::NoWeightedStep => {
                "ningún paso viable tiene peso positivo con un sesgo finito en [0,1]".to_owned()
            }
            Self::EventCountOverflow => "el recuento de eventos desborda".to_owned(),
            Self::DegenerateFrequencies { sum } => {
                format!("las frecuencias generadas suman {sum}")
            }
            Self::Subseed { detail } => detail.clone(),
            Self::Random(error) => format!("{error:?}"),
            Self::Tract(error) => error.detail(),
        }
    }

    /// El rechazo ya como diagnóstico de configuración, nombrando el locus.
    ///
    /// Usa el mismo código que el resto de los rechazos de esta frontera, con
    /// el parámetro exacto en `fields` y el locus dentro del detalle: es lo que
    /// permite corregir el documento sin abrir el código.
    #[must_use]
    pub fn into_diagnostic(self, path: &str, locus_id: &str) -> Diagnostic {
        Diagnostic::configuration(
            "GPS-CONF-008",
            "La configuración contiene una combinación de campos incompatible.",
        )
        .with("path", path)
        .with("fields", self.parameter())
        .with("detail", format!("locus {locus_id}: {}", self.detail()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tandem_repeat::RepeatMotif;

    /// Fuente de palabras que **cuenta** lo que se le pide.
    struct Counting {
        words: Vec<u64>,
        cursor: usize,
        consumed: usize,
    }

    impl Counting {
        fn new(words: &[u64]) -> Self {
            Self {
                words: words.to_vec(),
                cursor: 0,
                consumed: 0,
            }
        }
    }

    impl RandomWordSource for Counting {
        fn next_word(&mut self) -> u64 {
            self.consumed += 1;
            let word = self.words.get(self.cursor).copied().unwrap_or(0);
            self.cursor += 1;
            word
        }
    }

    fn motif() -> RepeatMotif {
        RepeatMotif::new("CA").expect("motivo válido")
    }

    fn tract(units: u64) -> TandemRepeatTract {
        TandemRepeatTract::pure(&motif(), units).expect("tracto puro")
    }

    fn settings(mean: f64, bias: f64) -> SlippageSettings {
        SlippageSettings::new(SLIPPAGE_MODEL_STEPWISE, 1, mean, bias).expect("parámetros válidos")
    }

    /// **Mil eventos con sesgo 1 son todos expansiones, y ninguno consume
    /// aleatoriedad.**
    ///
    /// El fixture está elegido para que la propiedad **se vea**: el tope se
    /// pone en 1005 unidades, mil por encima del tracto de partida, de modo que
    /// los mil eventos caben sin topar. Con un tope alcanzable, la prueba
    /// mediría el tope y no el sesgo.
    #[test]
    fn a_thousand_events_with_full_expansion_bias_only_expand_and_consume_nothing() {
        let mut rng = Counting::new(&[]);
        let run =
            apply_slippage_events_with_source(&mut rng, &tract(5), 1005, settings(0.0, 1.0), 1000)
                .expect("tanda válida");
        assert_eq!(run.applied_events(), 1000);
        assert_eq!(run.unapplied_events(), 0);
        assert_eq!(run.tract().annotation().complete_copies(), 1005);
        assert_eq!(rng.consumed, 0, "el sesgo extremo no sortea nada");
    }

    /// **Mil eventos con sesgo 0 sobre un tracto de cinco unidades lo dejan en
    /// una, sin error, y el recuento de no aplicados vale exactamente 996.**
    ///
    /// La cifra exacta es lo que discrimina: sin ella, un servicio que ignorase
    /// los eventos sobrantes y otro que los contase darían el mismo tracto
    /// final. Y el tope se pone muy por encima para que lo que detiene la
    /// contracción sea el **agotamiento del tracto** y no el tope.
    #[test]
    fn a_thousand_contractions_stop_at_one_unit_and_count_the_rest_apart() {
        let mut rng = Counting::new(&[]);
        let run =
            apply_slippage_events_with_source(&mut rng, &tract(5), 1000, settings(0.0, 0.0), 1000)
                .expect("tanda válida");
        assert_eq!(run.applied_events(), 4);
        assert_eq!(run.unapplied_events(), 996);
        assert_eq!(run.tract().annotation().complete_copies(), 1);
        assert_eq!(rng.consumed, 0, "el sesgo extremo no sortea nada");
    }

    /// **Cada evento aplicado cambia la longitud en exactamente una unidad.**
    ///
    /// Con sesgo intermedio los dos pasos son viables y el ensayo binario sí
    /// consume, que es la mitad que la prueba anterior no puede ejercitar.
    #[test]
    fn every_applied_event_moves_the_tract_by_exactly_one_unit() {
        let motif_length = motif().length();
        let mut current = tract(20);
        // Palabras alternas: la mitad alta hace verdadero el ensayo binario.
        let mut rng = Counting::new(&[u64::MAX, 0, u64::MAX, 0, u64::MAX, 0]);
        for _ in 0..6 {
            let before = current.length();
            let (next, outcome) =
                apply_slippage_event_with_source(&mut rng, &current, 100, settings(1.0, 0.5))
                    .expect("evento válido");
            let SlippageEvent::Applied(step) = outcome else {
                panic!("con dos pasos viables el evento se aplica");
            };
            assert_eq!(step.abs(), 1, "el paso es de una unidad");
            assert_eq!(
                next.length().abs_diff(before),
                motif_length,
                "la longitud cambia en una unidad completa"
            );
            current = next;
        }
        assert_eq!(rng.consumed, 6, "cada elección real consume una palabra");
    }

    /// **El tope detiene el linaje, y el evento se cuenta como no aplicado.**
    ///
    /// El fixture pone el tope **dos unidades por encima** del tracto de
    /// partida: así se ven los dos lados, los eventos que caben y los que no.
    #[test]
    fn the_unit_cap_stops_expansion_and_counts_the_event_apart() {
        let mut rng = Counting::new(&[]);
        let run = apply_slippage_events_with_source(&mut rng, &tract(5), 7, settings(0.0, 1.0), 10)
            .expect("tanda válida");
        assert_eq!(run.applied_events(), 2);
        assert_eq!(run.unapplied_events(), 8);
        assert_eq!(run.tract().annotation().complete_copies(), 7);
    }

    /// **Un modelo distinto y un salto de dos unidades se rechazan antes de
    /// simular, nombrando el parámetro y el valor.**
    #[test]
    fn another_model_and_a_two_unit_jump_are_rejected_naming_parameter_and_value() {
        let rejected = SlippageSettings::new("infinite_alleles", 1, 1.0, 0.5)
            .expect_err("un modelo distinto se rechaza");
        assert_eq!(rejected.parameter(), "reference.repeat_loci.slippage.model");
        assert!(rejected.detail().contains("infinite_alleles"));
        let rejected =
            SlippageSettings::new(SLIPPAGE_MODEL_STEPWISE, 2, 1.0, 0.5).expect_err("dos unidades");
        assert_eq!(
            rejected.parameter(),
            "reference.repeat_loci.slippage.step_units"
        );
        assert!(rejected.detail().contains('2'));
    }

    /// **Con media cero el abanico es solo el alelo de la referencia y no se
    /// consume ni una palabra.**
    #[test]
    fn a_zero_mean_returns_only_the_reference_allele_without_consuming_randomness() {
        let ancestral = tract(12);
        let mut rng = Counting::new(&[]);
        let fan =
            generate_allele_fan_with_source(&mut rng, &ancestral, 100, settings(0.0, 0.5), 40)
                .expect("abanico válido");
        assert_eq!(fan.alleles().len(), 1);
        assert_eq!(fan.alleles()[0].tract().bases(), ancestral.bases());
        assert_eq!(fan.alleles()[0].frequency(), 1.0);
        assert_eq!(fan.lineages(), 40);
        assert_eq!(fan.drawn_events(), 0);
        assert_eq!(fan.applied_events(), 0);
        assert_eq!(fan.unapplied_events(), 0);
        assert_eq!(rng.consumed, 0);
    }

    /// **Ningún alelo del abanico tiene fracciones de unidad**, y el fixture
    /// está elegido para que el abanico **tenga varios alelos**: con un solo
    /// alelo la comprobación se cumpliría sola.
    #[test]
    fn no_generated_allele_carries_a_fractional_unit() {
        let root = RootSeed::parse(
            "root",
            "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
        )
        .expect("semilla de 256 bits");
        let fan = generate_allele_fan(
            root,
            "ref-main",
            "STR-1",
            &tract(30),
            60,
            settings(3.0, 0.5),
            200,
        )
        .expect("abanico válido");
        assert!(
            fan.alleles().len() > 1,
            "el fixture tiene que producir varios alelos para que la propiedad se vea"
        );
        for allele in fan.alleles() {
            let units = allele.tract().annotation().units_as_f64();
            assert_eq!(units.fract(), 0.0, "alelo con fracción de unidad: {units}");
        }
        let total: u64 = fan.alleles().iter().map(SlippageAllele::lineages).sum();
        assert_eq!(
            total,
            fan.lineages(),
            "los linajes se reparten entre los alelos"
        );
    }

    /// **La misma semilla produce el mismo abanico.**
    #[test]
    fn the_same_seed_produces_the_same_fan() {
        let root = RootSeed::parse(
            "root",
            "101112131415161718191a1b1c1d1e1f000102030405060708090a0b0c0d0e0f",
        )
        .expect("semilla de 256 bits");
        let first = generate_allele_fan(
            root,
            "ref-main",
            "STR-1",
            &tract(25),
            60,
            settings(2.0, 0.4),
            120,
        )
        .expect("abanico");
        let second = generate_allele_fan(
            root,
            "ref-main",
            "STR-1",
            &tract(25),
            60,
            settings(2.0, 0.4),
            120,
        )
        .expect("abanico");
        assert_eq!(first, second);
    }

    /// **El generador rechaza un tramo sin ni una copia completa del motivo.**
    #[test]
    fn an_ancestral_without_a_complete_unit_is_rejected() {
        let ancestral =
            TandemRepeatTract::new(motif(), "CGTGT").expect("tracto sin copia completa de CA");
        let mut rng = Counting::new(&[]);
        let rejected =
            generate_allele_fan_with_source(&mut rng, &ancestral, 100, settings(1.0, 0.5), 10)
                .expect_err("sin unidad completa no hay ancestro");
        assert_eq!(rejected, SlippageError::AncestralWithoutCompleteUnit);
    }

    /// **En el tope el linaje se DETIENE, no se le empuja hacia abajo.**
    ///
    /// # El fixture usa sesgo intermedio, y esa es toda la prueba
    ///
    /// Con sesgo exactamente 0 o exactamente 1 —los dos que usan las pruebas
    /// vecinas— esta contradicción es **invisible**, porque el paso que el borde
    /// deja fuera es justo el que el sesgo ya anulaba. Con sesgo 0,5 los dos
    /// pasos tienen peso, de modo que la versión anterior aplicaba la
    /// contracción de forma determinista y el linaje se iba hacia abajo sin que
    /// el sesgo pintara nada. Aquí eso muere.
    ///
    /// **Control positivo**: el mismo fixture una unidad por debajo del tope sí
    /// se mueve y sí consume, de modo que el resultado no es «aquí nunca pasa
    /// nada».
    #[test]
    fn at_the_cap_the_lineage_stops_instead_of_being_pushed_down() {
        let mut rng = Counting::new(&[u64::MAX, u64::MAX, u64::MAX, u64::MAX, u64::MAX]);
        let run =
            apply_slippage_events_with_source(&mut rng, &tract(10), 10, settings(1.0, 0.5), 5)
                .expect("tanda válida");
        assert_eq!(
            run.applied_events(),
            0,
            "el tope no puede empujar al linaje"
        );
        assert_eq!(run.unapplied_events(), 5);
        assert_eq!(run.tract().annotation().complete_copies(), 10);
        assert_eq!(rng.consumed, 0, "un linaje detenido no sortea");

        // Control positivo: una unidad por debajo del tope sí hay dos pasos
        // viables, el sorteo ocurre y el tracto se mueve.
        let mut rng = Counting::new(&[u64::MAX, 0, u64::MAX, 0, u64::MAX]);
        let run = apply_slippage_events_with_source(&mut rng, &tract(9), 10, settings(1.0, 0.5), 1)
            .expect("tanda válida");
        assert_eq!(run.applied_events(), 1);
        assert_eq!(run.unapplied_events(), 0);
        assert_eq!(rng.consumed, 1);
    }

    /// **En el suelo el linaje se DETIENE, no se le empuja hacia arriba.**
    ///
    /// Es la cara simétrica de la anterior y hace falta igual: la regla del
    /// borde tiene dos lados y una sola de las dos pruebas dejaría el otro sin
    /// cubrir. Mismo fixture de sesgo intermedio, misma razón.
    #[test]
    fn at_the_floor_the_lineage_stops_instead_of_being_pushed_up() {
        let mut rng = Counting::new(&[0, 0, 0, 0, 0]);
        let run = apply_slippage_events_with_source(&mut rng, &tract(1), 20, settings(1.0, 0.5), 5)
            .expect("tanda válida");
        assert_eq!(
            run.applied_events(),
            0,
            "el suelo no puede empujar al linaje"
        );
        assert_eq!(run.unapplied_events(), 5);
        assert_eq!(run.tract().annotation().complete_copies(), 1);
        assert_eq!(rng.consumed, 0);

        // Control positivo: con dos unidades los dos pasos son viables.
        let mut rng = Counting::new(&[u64::MAX]);
        let run = apply_slippage_events_with_source(&mut rng, &tract(2), 20, settings(1.0, 0.5), 1)
            .expect("tanda válida");
        assert_eq!(run.applied_events(), 1);
        assert_eq!(rng.consumed, 1);
    }

    /// **Un tramo de referencia que ya rebasa el tope del locus se rechaza.**
    ///
    /// El mismo alelo escrito a mano se rechaza desde la entrega 1; generado se
    /// admitía, y encima dejaba el generador **unidireccional**, porque desde
    /// encima del tope la expansión nunca es viable y el sesgo declarado deja de
    /// significar nada.
    ///
    /// **Control positivo**: con el tope justo en las unidades del ancestro se
    /// acepta, de modo que el rechazo es por rebasar y no por rozar.
    #[test]
    fn an_ancestral_tract_above_the_locus_cap_is_rejected() {
        let mut rng = Counting::new(&[]);
        let rejected =
            generate_allele_fan_with_source(&mut rng, &tract(10), 5, settings(1.0, 0.5), 4)
                .expect_err("el ancestro rebasa el tope");
        assert_eq!(
            rejected,
            SlippageError::AncestralAboveCap {
                units: 10.0,
                max_units: 5
            }
        );
        assert_eq!(rejected.parameter(), "reference.repeat_loci.max_units");
        assert!(rejected.detail().contains("10") && rejected.detail().contains('5'));

        let mut rng = Counting::new(&[]);
        generate_allele_fan_with_source(&mut rng, &tract(10), 10, settings(0.0, 0.5), 4)
            .expect("con el tope justo en el ancestro se acepta");
    }

    /// **El abanico exacto de un fixture calculado a mano, palabra a palabra.**
    ///
    /// # Por qué hace falta, dicho sin rodeos
    ///
    /// La afirmación central del modelo es que **la frecuencia de cada alelo es
    /// la proporción de linajes que terminan en él**. Comprobar que hay varios
    /// alelos, que las frecuencias suman uno y que dos ejecuciones con la misma
    /// semilla coinciden **no la comprueba**: un abanico con el espectro
    /// emparejado al revés cumple las tres. Aquí se fija **qué alelo recibe qué
    /// linaje**, que es justo lo que aquello no mira.
    ///
    /// # La aritmética, para que se pueda seguir sin ejecutar nada
    ///
    /// La fuente entrega palabras leídas de una lista. Con media 1 el recuento
    /// de Poisson multiplica uniformes hasta bajar del umbral `e⁻¹ ≈ 0,3679`:
    /// una palabra de todo unos da un uniforme de casi exactamente 1 y **no**
    /// baja del umbral, y una palabra de ceros da un uniforme de `1,1·10⁻¹⁶` y
    /// lo cruza de golpe. Y el ensayo del sesgo compara `palabra mód 2⁵³` con la
    /// mitad: **cero expande y todo unos contrae**, que es al revés de lo que
    /// parece y por eso se dice aquí.
    ///
    /// - **Linaje 1**, palabras `[0]`: cero pasos. Termina en 10 unidades.
    /// - **Linaje 2**, palabras `[MAX, 0]` y luego `[0]`: un paso, expansión.
    ///   Termina en 11 unidades.
    /// - **Linaje 3**, palabras `[MAX, MAX, 0]` y luego `[MAX, MAX]`: dos pasos,
    ///   las dos contracciones. Termina en 8 unidades.
    ///
    /// Nueve palabras en total, tres alelos con un linaje cada uno.
    #[test]
    fn a_hand_computed_fixture_produces_exactly_the_expected_fan() {
        const MAX: u64 = u64::MAX;
        let mut rng = Counting::new(&[
            0, // linaje 1: cero pasos
            MAX, 0, // linaje 2: un paso
            0, // linaje 2: el paso expande
            MAX, MAX, 0, // linaje 3: dos pasos
            MAX, MAX, // linaje 3: los dos contraen
        ]);
        let fan = generate_allele_fan_with_source(&mut rng, &tract(10), 20, settings(1.0, 0.5), 3)
            .expect("abanico válido");
        assert_eq!(rng.consumed, 9, "la aritmética de la docstring no cuadra");
        let observed: Vec<(u64, u64)> = fan
            .alleles()
            .iter()
            .map(|allele| {
                (
                    allele.tract().annotation().complete_copies(),
                    allele.lineages(),
                )
            })
            .collect();
        // Orden canónico: por longitud del tracto. 8, 10 y 11 unidades.
        assert_eq!(observed, vec![(8, 1), (10, 1), (11, 1)]);
        for allele in fan.alleles() {
            assert_eq!(allele.frequency(), 1.0 / 3.0);
        }
        assert_eq!(fan.drawn_events(), 3);
        assert_eq!(fan.applied_events(), 3);
        assert_eq!(fan.unapplied_events(), 0);
    }

    /// **Una media de pasos por encima del techo se rechaza antes de simular.**
    ///
    /// **Control positivo**: el techo exacto se acepta, de modo que el rechazo
    /// es por rebasarlo y no por rozarlo.
    #[test]
    fn a_mean_above_the_ceiling_is_rejected_and_the_ceiling_itself_is_accepted() {
        let rejected = SlippageSettings::new(
            SLIPPAGE_MODEL_STEPWISE,
            1,
            SLIPPAGE_MAX_MEAN_EVENTS + 1.0,
            0.5,
        )
        .expect_err("por encima del techo");
        assert_eq!(
            rejected.parameter(),
            "reference.repeat_loci.slippage.mean_events"
        );
        assert!(rejected.detail().contains("1000"));
        SlippageSettings::new(SLIPPAGE_MODEL_STEPWISE, 1, SLIPPAGE_MAX_MEAN_EVENTS, 0.5)
            .expect("el techo exacto se acepta");
    }

    /// **Sin copias haplotípicas no hay abanico que generar.**
    #[test]
    fn a_population_without_haplotype_copies_is_rejected() {
        let mut rng = Counting::new(&[]);
        let rejected =
            generate_allele_fan_with_source(&mut rng, &tract(10), 100, settings(1.0, 0.5), 0)
                .expect_err("cero linajes");
        assert_eq!(rejected, SlippageError::NoLineages);
    }
}
