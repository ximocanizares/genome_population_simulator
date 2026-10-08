//! Tartamudeo clonal heredable por linaje (F4.10).
//!
//! # Qué es, en una frase
//!
//! Con preparación amplificada, el deslizamiento de la polimerasa puede
//! ocurrir en **cualquier ciclo** de la PCR, y una copia que resbala transmite
//! su número de unidades a **toda su descendencia**. Un resbalón en el primer
//! ciclo alcanza a media familia; uno en el último, a una sola copia.
//!
//! # El fenómeno que de verdad importa
//!
//! Una familia de duplicados puede mostrar **dos o tres alelos en proporciones
//! que imitan a un heterocigoto**, cuando procede de una sola molécula. En
//! genética forense ese artefacto se confunde con un genotipo real, y es
//! exactamente lo que este módulo reproduce.
//!
//! # Por qué esto no puede existir sin el tartamudeo de base
//!
//! Las dos columnas de la calibración —con PCR y sin PCR— son **totales, no
//! componentes**. Este módulo **despeja** el incremento propio de la
//! amplificación, `r_pcr = 1 − (1−r_con)/(1−r_sin)`, y lo compone después con
//! la parte de base, de modo que `1 − (1−r_sin)(1−r_pcr) = r_con`, exactamente
//! la fracción medida. Una preparación amplificada entregada **sin** la parte
//! de base se quedaría en el incremento despejado y subestimaría el tartamudeo
//! —con la calibración por omisión, entre un 14 % y un 28 % en motivos de uno a
//! cinco pares de base; con el multiplicador de los homopolímeros de DTC-099
//! la cifra del homopolímero cambia—. La admisión de la configuración lo
//! impide.
//!
//! # Honestidad sobre lo que la resta puede y no puede afirmar
//!
//! Las dos columnas proceden de **dos individuos distintos**, de proyectos
//! distintos y con coberturas distintas. Que cada columna sea un total sigue
//! en pie, y con ello la resta; lo que se debilita es atribuir la diferencia
//! **solo** a la PCR, porque puede llevar dentro efectos de individuo y de
//! lote. En consecuencia **las tasas por ciclo derivadas son toscas**, y el
//! manifiesto las publica marcadas como tales en vez de callarlo.
//!
//! # La genealogía es una aproximación declarada
//!
//! No es el coalescente exacto del proceso de ramificación, del que no se
//! afirma tener forma cerrada. La pareja `0` de cada familia tiene linaje
//! propio completo; la pareja `i ≥ 1` elige un progenitor uniforme entre las
//! anteriores y un ciclo de divergencia geométrico, hereda los resbalones del
//! progenitor **anteriores** a la divergencia y sortea los posteriores.

use std::collections::BTreeMap;

use crate::identity::{EventId, LibraryId, RepeatLocusId};
use crate::randomness::{
    F4TechnicalDecision, RootSeed, categorical_f64, derive_f4_technical_subseed, f4_bernoulli,
    f4_uniform_index,
};
use crate::stutter::{
    STUTTER_CALIBRATION_SOURCE, STUTTER_CALIBRATION_SPECIES, STUTTER_CALIBRATION_TECHNOLOGY,
    STUTTER_DEFAULT_MOTIF_LENGTH_RATES, STUTTER_MOTIF_LENGTH_MAX,
    STUTTER_WITH_PCR_MOTIF_LENGTH_RATES, StutterError, StutterStepWeights, choose_step_with_source,
    viable_steps,
};

/// Versión del algoritmo, publicada en la configuración efectiva y en el
/// manifiesto. Cambiar la genealogía o la derivación obliga a estrenar
/// versión.
pub const CLONAL_STUTTER_ALGORITHM: &str = "clonal-stutter-v2";

/// La derivación, escrita tal cual para que viaje con los números.
pub const CLONAL_STUTTER_FORMULA: &str =
    "r_pcr = 1 - (1 - r_with_pcr) / (1 - r_without_pcr); p_cycle = 1 - (1 - r_pcr)^(1/cycles)";

/// **Marca de tosquedad**, publicada siempre que estas tasas se usen.
pub const CLONAL_STUTTER_COARSE: bool = true;

/// Por qué son toscas, dicho entero para que nadie tenga que buscarlo.
pub const CLONAL_STUTTER_COARSE_REASON: &str = "the two calibration columns come from two different individuals, from different projects and \
     at different coverages, so the subtraction may carry individual and batch effects besides \
     the amplification: these per-cycle rates are coarse and are published as such";

/// Menor multiplicador de la PCR admitido en los homopolímeros (DTC-099).
///
/// Con 1 la fracción con PCR iguala a la sin PCR: la amplificación no añade
/// tartamudeo, la probabilidad por ciclo vale cero y no se consume
/// aleatoriedad. Por debajo diría que la PCR **quita** tartamudeo, que la
/// derivación no representa: solo despeja un incremento no negativo.
pub const HOMOPOLYMER_PCR_FOLD_MINIMUM: f64 = 1.0;

/// Mayor multiplicador de la PCR admitido en los homopolímeros (DTC-099).
///
/// Es el mayor entero que deja la fracción con PCR por debajo de 1:
/// `37 × 0,027 = 0,999`. **Se deriva de la constante sin PCR del
/// homopolímero**, `STUTTER_DEFAULT_MOTIF_LENGTH_RATES[0]`: si esa constante
/// cambiara, este tope hay que recalcularlo. Una prueba lo comprueba.
pub const HOMOPOLYMER_PCR_FOLD_MAXIMUM: f64 = 37.0;

/// Fallo estable de los parámetros del tartamudeo clonal que declara el
/// escenario.
#[derive(Clone, Debug, PartialEq)]
pub enum ClonalStutterError {
    /// Longitud de motivo fuera del tramo calibrado.
    UnknownMotifLength { motif_length: u64 },
    /// Probabilidad fuera de `[0, 1]` o no finita.
    InvalidProbability { motif_length: u64, value: f64 },
    /// Multiplicador de los homopolímeros no finito o fuera de
    /// `[HOMOPOLYMER_PCR_FOLD_MINIMUM, HOMOPOLYMER_PCR_FOLD_MAXIMUM]`.
    InvalidHomopolymerFold { value: f64 },
    /// Multiplicador declarado a la vez que una sustitución directa de la
    /// probabilidad por ciclo de los homopolímeros: dos valores para lo mismo.
    HomopolymerFoldWithSubstitution { value: f64 },
}

/// Comprueba el dominio del multiplicador de los homopolímeros.
///
/// Rechaza en vez de pinzar: un valor por debajo de 1 llegaría a
/// [`per_cycle_slip_probability`] como una fracción con PCR menor que la sin
/// PCR, y esa función devolvería cero sin avisar. El tramo cerrado rechaza
/// también `NaN` y los infinitos: ninguno de los tres está contenido en él.
pub fn validate_homopolymer_pcr_fold(value: f64) -> Result<f64, ClonalStutterError> {
    if (HOMOPOLYMER_PCR_FOLD_MINIMUM..=HOMOPOLYMER_PCR_FOLD_MAXIMUM).contains(&value) {
        Ok(value)
    } else {
        Err(ClonalStutterError::InvalidHomopolymerFold { value })
    }
}

/// Probabilidad de resbalón **por ciclo** despejada de las dos columnas.
///
/// # Las dos mitades de la fórmula
///
/// 1. **Despejar el incremento de la PCR**, `r_pcr = 1 − (1−r_con)/(1−r_sin)`.
///    Es una resta, no una suma: cada columna es un total, y sumarlas contaría
///    dos veces el mosaicismo somático, la amplificación en puente y la
///    síntesis, que ya están dentro de las dos.
/// 2. **Repartirlo entre los ciclos**, `p = 1 − (1 − r_pcr)^(1/c)`, que es la
///    tasa por paso de linaje cuya acumulación en `c` ciclos reproduce
///    `r_pcr`.
///
/// La potencia sale de `libm`, no de la matemática estándar: es la política
/// vigente del proyecto y aquí además es lo que la hace reproducible entre
/// máquinas.
///
/// # Los casos que devuelven cero
///
/// Si la fracción con PCR **no supera** a la de sin PCR, la probabilidad vale
/// cero y no se consume aleatoriedad. Es lo que corresponde: no hay incremento
/// que repartir, y fabricar uno negativo sería inventar. Lo mismo con una
/// columna sin PCR que valiera 1, donde el cociente no está definido.
#[must_use]
pub fn per_cycle_slip_probability(with_pcr: f64, without_pcr: f64, cycles: u32) -> f64 {
    if cycles == 0
        || !with_pcr.is_finite()
        || !without_pcr.is_finite()
        || !(0.0..=1.0).contains(&with_pcr)
        || !(0.0..1.0).contains(&without_pcr)
        || with_pcr <= without_pcr
    {
        return 0.0;
    }
    let increment = 1.0 - (1.0 - with_pcr) / (1.0 - without_pcr);
    if !increment.is_finite() || increment <= 0.0 {
        return 0.0;
    }
    if increment >= 1.0 {
        return 1.0;
    }
    let probability = 1.0 - libm::pow(1.0 - increment, 1.0 / f64::from(cycles));
    if !probability.is_finite() {
        return 0.0;
    }
    probability.clamp(0.0, 1.0)
}

/// De dónde sale la probabilidad por ciclo de una longitud de motivo.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PerCycleSlipSource {
    /// Despejada de las dos columnas de la calibración y del número de ciclos.
    Derived,
    /// Sustituida directamente por el autor del escenario, que prevalece.
    Declared,
}

impl PerCycleSlipSource {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Derived => "derived",
            Self::Declared => "declared",
        }
    }
}

/// Parámetros efectivos del tartamudeo clonal de un evento.
#[derive(Clone, Debug, PartialEq)]
pub struct ClonalStutterSettings {
    cycles: u32,
    /// La columna **con PCR efectiva**: la de la calibración, salvo el
    /// homopolímero cuando el escenario declara su multiplicador. La leen los
    /// tres sitios que la necesitan —la derivación publicada, la del sorteo y
    /// el documento publicado— para que no puedan discrepar.
    with_pcr: [f64; STUTTER_MOTIF_LENGTH_MAX],
    homopolymer_pcr_fold: Option<f64>,
    per_cycle: [f64; STUTTER_MOTIF_LENGTH_MAX],
    source: [PerCycleSlipSource; STUTTER_MOTIF_LENGTH_MAX],
}

impl ClonalStutterSettings {
    /// Deriva las seis probabilidades por ciclo y aplica las sustituciones.
    ///
    /// **Las dos columnas de entrada son las constantes de la calibración**, no
    /// las tasas efectivas que el evento pueda haber sustituido. Es
    /// deliberado: `r_pcr` es lo que **la amplificación añade**, una propiedad
    /// del paso de PCR que no cambia porque el autor declare que su biblioteca
    /// tartamudea menos en la agrupación y la síntesis. Quien baje la tasa de
    /// base verá bajar la fracción de extremo a extremo, que es lo que ha
    /// pedido; quien quiera mover la parte de PCR tiene la sustitución directa.
    ///
    /// # El multiplicador de los homopolímeros (DTC-099)
    ///
    /// Con `homopolymer_pcr_fold = Some(k)`, la fracción con PCR del
    /// homopolímero pasa a ser `k × 0,027`, **la constante** sin PCR, y de ahí
    /// se despeja como siempre. Sin él se usa la constante con PCR, 0,17,
    /// **tal cual**: no se calcula `6,30 × 0,027`, que daría 0,1701 y movería
    /// los bytes publicados. Las longitudes de dos a seis no cambian nunca.
    pub fn new(
        cycles: u32,
        substitutions: &BTreeMap<u64, f64>,
        homopolymer_pcr_fold: Option<f64>,
    ) -> Result<Self, ClonalStutterError> {
        let mut with_pcr = STUTTER_WITH_PCR_MOTIF_LENGTH_RATES;
        if let Some(fold) = homopolymer_pcr_fold {
            let fold = validate_homopolymer_pcr_fold(fold)?;
            if substitutions.contains_key(&1) {
                return Err(ClonalStutterError::HomopolymerFoldWithSubstitution { value: fold });
            }
            with_pcr[0] = fold * STUTTER_DEFAULT_MOTIF_LENGTH_RATES[0];
        }
        let mut per_cycle = [0.0_f64; STUTTER_MOTIF_LENGTH_MAX];
        let mut source = [PerCycleSlipSource::Derived; STUTTER_MOTIF_LENGTH_MAX];
        for (index, slot) in per_cycle.iter_mut().enumerate() {
            *slot = per_cycle_slip_probability(
                with_pcr[index],
                STUTTER_DEFAULT_MOTIF_LENGTH_RATES[index],
                cycles,
            );
        }
        for (motif_length, value) in substitutions {
            let index = usize::try_from(*motif_length)
                .ok()
                .filter(|length| (1..=STUTTER_MOTIF_LENGTH_MAX).contains(length))
                .ok_or(ClonalStutterError::UnknownMotifLength {
                    motif_length: *motif_length,
                })?
                - 1;
            if !value.is_finite() || !(0.0..=1.0).contains(value) {
                return Err(ClonalStutterError::InvalidProbability {
                    motif_length: *motif_length,
                    value: *value,
                });
            }
            per_cycle[index] = *value;
            source[index] = PerCycleSlipSource::Declared;
        }
        Ok(Self {
            cycles,
            with_pcr,
            homopolymer_pcr_fold,
            per_cycle,
            source,
        })
    }

    #[must_use]
    pub const fn cycles(&self) -> u32 {
        self.cycles
    }

    /// Probabilidad por ciclo de una longitud de motivo **a la longitud de
    /// tracto de referencia**, o `None` fuera del tramo calibrado.
    ///
    /// Es la que publica el manifiesto, y la que se usaba antes de que el
    /// modelo dependiera de la longitud del tracto. El sorteo usa
    /// [`Self::probability_for`], que es esta modulada.
    #[must_use]
    pub fn probability_for_motif_length(&self, motif_length: u64) -> Option<f64> {
        usize::try_from(motif_length)
            .ok()
            .filter(|length| (1..=STUTTER_MOTIF_LENGTH_MAX).contains(length))
            .map(|length| self.per_cycle[length - 1])
    }

    /// **La probabilidad por ciclo que el sorteo usa de verdad**, ya modulada
    /// por la longitud del tracto.
    ///
    /// # Por qué se modulan las DOS columnas con el MISMO factor
    ///
    /// Porque de la diferencia entre ellas es de donde sale esta
    /// probabilidad, y esa derivación es un modelo ya cerrado. Modular las dos
    /// por el mismo factor `m` deja el cociente entre columnas **independiente
    /// de la longitud**, que es lo que la derivación supone, y hace que la
    /// composición siga cerrando **exactamente** a cualquier longitud:
    ///
    /// ```text
    /// 1 − (1 − r_sin·m)(1 − r_pcr(m)) = r_con·m
    /// ```
    ///
    /// Modular solo la parte de base rompería esa identidad en todas las
    /// longitudes salvo la de referencia.
    ///
    /// # La sustitución del autor prevalece, y no se modula
    ///
    /// Misma regla que la tasa de base: quien sustituye la probabilidad por
    /// ciclo la ha elegido él, y la curva no tiene nada que corregirle.
    #[must_use]
    pub fn probability_for(&self, motif_length: u64, run_units: u64) -> Option<f64> {
        let index = usize::try_from(motif_length)
            .ok()
            .filter(|length| (1..=STUTTER_MOTIF_LENGTH_MAX).contains(length))?
            - 1;
        if self.source[index] == PerCycleSlipSource::Declared {
            return Some(self.per_cycle[index]);
        }
        let factor = crate::stutter::tract_length_factor(motif_length, run_units);
        Some(per_cycle_slip_probability(
            self.with_pcr[index] * factor,
            STUTTER_DEFAULT_MOTIF_LENGTH_RATES[index] * factor,
            self.cycles,
        ))
    }
}

/// Documento común de los parámetros efectivos, que publican **la
/// configuración efectiva y el manifiesto**.
///
/// Sale de aquí para que no puedan discrepar, y lleva **las dos fracciones de
/// entrada, el número de ciclos, la fórmula y la marca de tosquedad**: unos
/// números sin su procedencia no se distinguen de unos inventados, y una tasa
/// tosca sin su marca se leería como una medida.
///
/// El multiplicador de los homopolímeros se publica **solo si se declaró**
/// (REQ-MAN-97): publicarlo siempre movería los bytes de toda ejecución
/// amplificada con microsatélites.
#[must_use]
pub fn effective_settings_document(settings: &ClonalStutterSettings) -> serde_json::Value {
    let per_cycle_slip = (1..=STUTTER_MOTIF_LENGTH_MAX)
        .map(|motif_length| {
            let index = motif_length - 1;
            serde_json::json!({
                "motif_length_bases": motif_length,
                "with_pcr_fraction": settings.with_pcr[index],
                "without_pcr_fraction": STUTTER_DEFAULT_MOTIF_LENGTH_RATES[index],
                "probability": settings.per_cycle[index],
                "source": settings.source[index].as_str(),
            })
        })
        .collect::<Vec<_>>();
    let mut document = serde_json::json!({
        "algorithm": CLONAL_STUTTER_ALGORITHM,
        "cycles": settings.cycles(),
        "formula": CLONAL_STUTTER_FORMULA,
        "coarse": CLONAL_STUTTER_COARSE,
        "coarse_reason": CLONAL_STUTTER_COARSE_REASON,
        "calibration": {
            "source": STUTTER_CALIBRATION_SOURCE,
            "species": STUTTER_CALIBRATION_SPECIES,
            "technology": STUTTER_CALIBRATION_TECHNOLOGY,
        },
        "per_cycle_slip": per_cycle_slip,
    });
    if let (Some(fold), Some(object)) = (settings.homopolymer_pcr_fold, document.as_object_mut()) {
        object.insert("homopolymer_pcr_fold".to_owned(), serde_json::json!(fold));
    }
    document
}

// ---------------------------------------------------------------------------
// F4.10.2 — la genealogía mínima.
// ---------------------------------------------------------------------------

/// Pesos del **ciclo de divergencia** `d ∈ 0..=c`, en su orden natural.
///
/// `P(d = j) = (1−q)·q^j` para `j < c`, y **toda** la masa restante en
/// `d = c`, con `q = 1/(1+e)`. El pinzamiento del último no es un detalle: sin
/// él la distribución perdería masa y el número de ciclos dejaría de ser un
/// tope.
///
/// Con `d = 0` la copia no comparte nada con su progenitor. Es el caso más
/// frecuente con eficiencia alta —con `e = 0,9`, `P(d=0) = e/(1+e) ≈ 0,47`—, y
/// corresponde a dos copias que descienden de productos distintos del primer
/// ciclo.
///
/// La potencia se calcula por multiplicación repetida: el exponente es un
/// entero acotado por el número de ciclos, así que no hace falta `libm`.
#[must_use]
pub fn divergence_weights(cycles: u32, efficiency: f64) -> Vec<f64> {
    let slots = cycles as usize + 1;
    let mut weights = vec![0.0_f64; slots];
    if !efficiency.is_finite() || efficiency <= 0.0 {
        // Sin eficiencia no hay ramificación: todo diverge en el último ciclo.
        weights[slots - 1] = 1.0;
        return weights;
    }
    let ratio = 1.0 / (1.0 + efficiency);
    let mut power = 1.0_f64;
    for slot in weights.iter_mut().take(slots - 1) {
        *slot = (1.0 - ratio) * power;
        power *= ratio;
    }
    weights[slots - 1] = power;
    weights
}

/// Pasos netos de linaje de **una** copia, uno por tracto.
///
/// Es lo único que sale de la genealogía hacia el recorrido: ni el árbol, ni
/// los ciclos, ni quién es progenitor de quién. El manifiesto tampoco los
/// publica por pareja (REQ-DAT-56); un validador independiente los rederiva
/// desde la semilla.
///
/// La lista de loci se **comparte** entre todas las copias de la familia en
/// vez de repetirse en cada una: con un identificador por copia y tracto, la
/// memoria de una familia dependería de cómo de largos sean los nombres que el
/// autor del escenario elija, que es justo lo que una cota publicada no puede
/// permitirse.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LineageNetSteps {
    loci: Vec<RepeatLocusId>,
    net: Vec<i64>,
}

impl LineageNetSteps {
    /// Paso neto de un locus, o cero si esta copia no lo desplazó.
    #[must_use]
    pub fn step_for(&self, locus_id: &str) -> i64 {
        self.loci
            .iter()
            .position(|candidate| candidate.as_str() == locus_id)
            .and_then(|index| self.net.get(index).copied())
            .unwrap_or(0)
    }

    /// Los pasos, en el orden canónico de los tractos de la familia.
    #[must_use]
    pub fn steps(&self) -> &[i64] {
        &self.net
    }
}

/// Un tracto de la molécula original de la familia, ya observado.
#[derive(Clone, Copy, Debug)]
pub struct LineageTract<'a> {
    locus_id: &'a RepeatLocusId,
    run_units: u64,
    motif_length: u64,
    probability: f64,
}

impl<'a> LineageTract<'a> {
    #[must_use]
    pub const fn new(
        locus_id: &'a RepeatLocusId,
        run_units: u64,
        motif_length: u64,
        probability: f64,
    ) -> Self {
        Self {
            locus_id,
            run_units,
            motif_length,
            probability,
        }
    }
}

/// Lo que la genealogía de una familia deja para el recorrido y para el
/// manifiesto.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FamilyLineage {
    loci: Vec<RepeatLocusId>,
    /// `net[copia × tractos + tracto]`, aplanado. Es exactamente el término
    /// `K_MAX × T × 8` que la cota de recursos declara.
    net: Vec<i64>,
    slips_by_locus: BTreeMap<RepeatLocusId, u64>,
    inherited_slips_by_locus: BTreeMap<RepeatLocusId, u64>,
}

impl FamilyLineage {
    /// Pasos netos de la copia `index` de la familia.
    #[must_use]
    pub fn member(&self, index: usize) -> Option<LineageNetSteps> {
        if self.loci.is_empty() {
            return None;
        }
        let width = self.loci.len();
        let start = index.checked_mul(width)?;
        let slice = self.net.get(start..start.checked_add(width)?)?;
        Some(LineageNetSteps {
            loci: self.loci.clone(),
            net: slice.to_vec(),
        })
    }

    /// Cuántos tractos entraron en la genealogía de esta familia.
    #[must_use]
    pub fn loci_count(&self) -> usize {
        self.loci.len()
    }

    /// Resbalones por locus: **todos** los ciclos con paso distinto de cero,
    /// heredados o propios.
    #[must_use]
    pub const fn slips_by_locus(&self) -> &BTreeMap<RepeatLocusId, u64> {
        &self.slips_by_locus
    }

    /// De los anteriores, los que la copia **recibió de su progenitor** en vez
    /// de sortearlos.
    ///
    /// Con familias de tamaño 1 vale **exactamente cero**, porque ahí no hay
    /// de quién heredar. Es lo que distingue este recuento de uno que contara
    /// todos los resbalones.
    #[must_use]
    pub const fn inherited_slips_by_locus(&self) -> &BTreeMap<RepeatLocusId, u64> {
        &self.inherited_slips_by_locus
    }
}

/// Identidad de la familia cuyos linajes se resuelven.
#[derive(Clone, Copy, Debug)]
pub struct LineageKey<'a> {
    pub root: RootSeed,
    pub event_id: &'a EventId,
    pub library_id: &'a LibraryId,
    pub family_ordinal: u64,
}

/// Geometría de la molécula sobre la que se sortea el linaje.
#[derive(Clone, Copy, Debug)]
pub struct LineageMolecule {
    pub length_bases: u64,
    pub read_length_bases: u64,
}

/// Resuelve la genealogía de una familia entera y devuelve el paso neto de
/// cada copia y cada tracto.
///
/// # El orden del sorteo, que es el contrato
///
/// Para cada copia, **ciclo a ciclo** y dentro de cada ciclo en el orden
/// canónico de los tractos. Recorrerlo por tractos y luego por ciclos daría
/// las mismas claves pero llevaría mal la cuenta de la longitud de la
/// molécula, que es lo que la guarda de la lectura mira.
///
/// # Las guardas se aplican al paso neto acumulado
///
/// Un resbalón entra en el sorteo solo si deja al menos una unidad completa en
/// la tirada **tal como está después de los resbalones anteriores de esta
/// misma copia**, y si deja la molécula al menos tan larga que una lectura.
/// Se excluye **antes** de sortear y no se reintenta. Aplicarlas a cada paso
/// por separado dejaría pasar linajes que vacían el tracto en dos mitades.
///
/// # Lo que esta contabilidad supone, dicho sin adornos
///
/// La tirada después de `k` unidades de paso se calcula **sumando** `k` a la
/// tirada observada, no volviendo a leer las bases. Para un tracto **puro** es
/// exacto. Para uno impuro puede no serlo, porque una contracción bastante
/// grande podría dejar otra tirada como la más larga; en ese caso la guarda
/// del sorteo y la del momento de aplicar discrepan, y la de aplicar manda:
/// [`crate::stutter::transform_paired_end`] vuelve a leer las bases y no
/// escribe un paso que no quepa.
///
/// # Por qué los pesos de la divergencia entran como argumento
///
/// Porque así el fixture saturado puede **forzarlos** sin tocar el código que
/// prueba: con toda la masa en el último ciclo, todas las copias heredan el
/// linaje entero; con toda la masa en el cero, ninguna hereda nada. Sin ese
/// asidero, la prueba de la herencia tendría que confiar en que un sorteo con
/// eficiencia extrema caiga siempre del mismo lado, que es precisamente la
/// clase de aserción que no discrimina. El recorrido real los calcula con
/// [`divergence_weights`].
pub fn resolve_family_lineage(
    key: LineageKey<'_>,
    family_size: u64,
    divergence_distribution: &[f64],
    cycles: u32,
    tracts: &[LineageTract<'_>],
    weights: &StutterStepWeights,
    molecule: LineageMolecule,
) -> Result<FamilyLineage, StutterError> {
    let members = usize::try_from(family_size).map_err(|_| StutterError::LengthOverflow)?;
    if members == 0 || tracts.is_empty() || cycles == 0 {
        return Ok(FamilyLineage::default());
    }
    let cycle_count = cycles as usize;
    // `steps[copia][tracto][ciclo − 1]`. Es el estado que la cota de recursos
    // declara como «los pasos netos de la familia en curso»: no crece con el
    // número de parejas de la biblioteca ni con el genoma.
    let mut steps = vec![vec![vec![0_i64; cycle_count]; tracts.len()]; members];
    if divergence_distribution.len() != cycle_count + 1 {
        return Err(StutterError::LengthOverflow);
    }
    let mut lineage = FamilyLineage {
        loci: tracts
            .iter()
            .map(|tract| tract.locus_id.clone())
            .collect::<Vec<_>>(),
        net: Vec::with_capacity(members.saturating_mul(tracts.len())),
        slips_by_locus: BTreeMap::new(),
        inherited_slips_by_locus: BTreeMap::new(),
    };
    for member in 0..members {
        let copy_index = member as u64;
        // La copia 0 tiene linaje propio completo. La copia `i ≥ 1` elige un
        // progenitor uniforme entre las anteriores; con `i = 1` hay un solo
        // candidato y `uniform-index-v1` devuelve sin consumir.
        let (parent, divergence) = if member == 0 {
            (0_usize, 0_usize)
        } else {
            let mut parent_rng = derive_f4_technical_subseed(
                key.root,
                F4TechnicalDecision::lineage_parent(
                    key.event_id,
                    key.library_id,
                    key.family_ordinal,
                    copy_index,
                ),
            )
            .map_err(StutterError::Pair)?
            .local_rng();
            let parent = f4_uniform_index(&mut parent_rng, copy_index)
                .map_err(StutterError::Random)? as usize;
            let mut divergence_rng = derive_f4_technical_subseed(
                key.root,
                F4TechnicalDecision::lineage_divergence(
                    key.event_id,
                    key.library_id,
                    key.family_ordinal,
                    copy_index,
                ),
            )
            .map_err(StutterError::Pair)?
            .local_rng();
            let divergence = categorical_f64(&mut divergence_rng, divergence_distribution)
                .map_err(StutterError::Random)?;
            (parent, divergence)
        };
        let mut net = vec![0_i64; tracts.len()];
        let mut length =
            i64::try_from(molecule.length_bases).map_err(|_| StutterError::LengthOverflow)?;
        for cycle in 1..=cycle_count {
            for (index, tract) in tracts.iter().enumerate() {
                let inherited = member > 0 && cycle <= divergence;
                let step = if inherited {
                    steps[parent][index][cycle - 1]
                } else {
                    draw_cycle_step(
                        key,
                        copy_index,
                        cycle as u64,
                        tract,
                        net[index],
                        length,
                        molecule.read_length_bases,
                        weights,
                    )?
                };
                if step == 0 {
                    continue;
                }
                steps[member][index][cycle - 1] = step;
                net[index] = net[index]
                    .checked_add(step)
                    .ok_or(StutterError::LengthOverflow)?;
                length = length
                    .checked_add(
                        step.saturating_mul(
                            i64::try_from(tract.motif_length)
                                .map_err(|_| StutterError::LengthOverflow)?,
                        ),
                    )
                    .ok_or(StutterError::LengthOverflow)?;
                // La clave es el **identificador**, no su texto copiado: un
                // `RepeatLocusId` es un `Arc<str>` y clonarlo es un contador,
                // no una asignación nueva. Con `String` la memoria de estos
                // recuentos crecería con la longitud de los nombres que el
                // autor del escenario elija, que es justo lo que una cota
                // publicada no puede permitirse.
                *lineage
                    .slips_by_locus
                    .entry(tract.locus_id.clone())
                    .or_insert(0) += 1;
                if inherited {
                    *lineage
                        .inherited_slips_by_locus
                        .entry(tract.locus_id.clone())
                        .or_insert(0) += 1;
                }
            }
        }
        lineage.net.extend_from_slice(&net);
    }
    // Los loci que nunca resbalaron se publican con cero en vez de faltar: un
    // recuento ausente y un cero se leen igual de mal, y solo uno de los dos
    // dice que el locus estuvo y no pasó nada.
    for tract in tracts {
        lineage
            .slips_by_locus
            .entry(tract.locus_id.clone())
            .or_insert(0);
        lineage
            .inherited_slips_by_locus
            .entry(tract.locus_id.clone())
            .or_insert(0);
    }
    Ok(lineage)
}

/// Un ciclo de una copia sobre un tracto: ¿resbala?, y si sí, ¿cuánto?
///
/// Con probabilidad exactamente cero, `bernoulli-v1` devuelve sin consumir ni
/// una palabra; con un solo candidato viable, el categórico tampoco consume.
#[expect(
    clippy::too_many_arguments,
    reason = "son los siete datos del sorteo de un ciclo; agruparlos en un tipo solo movería la lista"
)]
fn draw_cycle_step(
    key: LineageKey<'_>,
    copy_index: u64,
    cycle: u64,
    tract: &LineageTract<'_>,
    accumulated: i64,
    length: i64,
    read_length: u64,
    weights: &StutterStepWeights,
) -> Result<i64, StutterError> {
    let mut slip_rng = derive_f4_technical_subseed(
        key.root,
        F4TechnicalDecision::lineage_slip(
            key.event_id,
            key.library_id,
            key.family_ordinal,
            copy_index,
            cycle,
            tract.locus_id,
        ),
    )
    .map_err(StutterError::Pair)?
    .local_rng();
    if !f4_bernoulli(&mut slip_rng, tract.probability).map_err(StutterError::Random)? {
        return Ok(0);
    }
    let run_units = i64::try_from(tract.run_units)
        .map_err(|_| StutterError::LengthOverflow)?
        .checked_add(accumulated)
        .ok_or(StutterError::LengthOverflow)?;
    if run_units <= 0 || length <= 0 {
        // **Inalcanzable con las guardas de arriba**: ningún paso admitido
        // puede vaciar la tirada ni dejar la molécula sin bases, porque
        // `step_is_viable` los excluye antes. La rama queda porque el coste es
        // nulo y porque la aritmética con signo la exige; ninguna prueba puede
        // ejercitarla.
        return Ok(0);
    }
    let candidates = viable_steps(
        run_units.unsigned_abs(),
        tract.motif_length,
        length.unsigned_abs(),
        read_length,
    );
    if candidates.is_empty() {
        // Tampoco es alcanzable con una tirada de al menos una unidad: la
        // expansión de una unidad sobrevive siempre al filtrado.
        return Ok(0);
    }
    let mut step_rng = derive_f4_technical_subseed(
        key.root,
        F4TechnicalDecision::lineage_step(
            key.event_id,
            key.library_id,
            key.family_ordinal,
            copy_index,
            cycle,
            tract.locus_id,
        ),
    )
    .map_err(StutterError::Pair)?
    .local_rng();
    Ok(choose_step_with_source(&mut step_rng, weights, &candidates)?.unwrap_or(0))
}

#[cfg(test)]
mod tests {
    use super::{
        ClonalStutterError, ClonalStutterSettings, HOMOPOLYMER_PCR_FOLD_MAXIMUM,
        HOMOPOLYMER_PCR_FOLD_MINIMUM, STUTTER_DEFAULT_MOTIF_LENGTH_RATES,
        STUTTER_WITH_PCR_MOTIF_LENGTH_RATES, divergence_weights, effective_settings_document,
        per_cycle_slip_probability, validate_homopolymer_pcr_fold,
    };
    use std::collections::BTreeMap;

    /// **DTC-099. El multiplicador de los homopolímeros llega a los DOS
    /// cálculos —el publicado y el del sorteo—, modulado, y no toca nada más.**
    ///
    /// Las cifras esperadas están escritas a mano desde las dos fracciones
    /// literales, no pedidas a otra llamada del módulo:
    ///
    /// - en el ancla (15 unidades), `r_pcr = 1 − 0,9325/0,973` y
    ///   `p = 1 − (1 − r_pcr)^(1/10) ≈ 0,0042425`;
    /// - en un tracto de 11 unidades, las dos fracciones multiplicadas por
    ///   `m(11) = 0,11885022274370183` dan `p ≈ 0,00048395`.
    ///
    /// **Debe fallar contra** (a) una entrega que trate el multiplicador como
    /// una sustitución sin modular: en el tracto corto daría 0,0042425, nueve
    /// veces más; (b) una que actualice la derivación publicada y olvide la del
    /// sorteo: en el tracto corto daría 0,00171825, la de 0,17; (c) una que
    /// aplique el multiplicador a todas las longitudes: las filas de 2 a 6 ya
    /// no serían bit a bit las de sin multiplicador.
    ///
    /// **En el ancla el factor de longitud vale 1** y la variante (a) no se
    /// distingue de la buena; por eso hace falta el tracto corto, y por eso se
    /// comprueba que `m(11)` es de verdad el literal y no 1.
    #[test]
    fn a_homopolymer_fold_reaches_both_derivations_modulated_and_nothing_else() {
        let none =
            ClonalStutterSettings::new(10, &BTreeMap::new(), None).expect("sin multiplicador");
        let fold = ClonalStutterSettings::new(10, &BTreeMap::new(), Some(2.5)).expect("con 2,5");
        assert_eq!(
            0.027_f64 * 2.5,
            0.0675,
            "la fracción con PCR resultante es exacta"
        );
        let anchor = fold.probability_for_motif_length(1).expect("homopolímero");
        assert!(
            (anchor - 0.004_242_5).abs() < 1e-7,
            "probabilidad publicada del homopolímero con 2,5: {anchor}"
        );
        assert_eq!(
            fold.probability_for(1, 15),
            Some(anchor),
            "en el ancla el sorteo usa la publicada"
        );
        let factor = crate::stutter::tract_length_factor(1, 11);
        assert_eq!(
            factor, 0.118_850_222_743_701_83,
            "el fixture del tracto corto"
        );
        let short = fold.probability_for(1, 11).expect("tracto corto");
        assert!(
            (short - 0.000_483_95).abs() < 1e-8,
            "probabilidad del sorteo en 11 unidades con 2,5: {short}"
        );
        for motif_length in 2..=6_u64 {
            assert_eq!(
                fold.probability_for_motif_length(motif_length)
                    .map(f64::to_bits),
                none.probability_for_motif_length(motif_length)
                    .map(f64::to_bits),
                "la longitud de motivo {motif_length} no puede moverse"
            );
            for units in [3_u64, 8, 17] {
                assert_eq!(
                    fold.probability_for(motif_length, units).map(f64::to_bits),
                    none.probability_for(motif_length, units).map(f64::to_bits),
                    "la longitud de motivo {motif_length} con {units} unidades no puede moverse"
                );
            }
        }
        // Control positivo: el homopolímero sí se mueve, en los dos caminos.
        assert!(anchor < none.probability_for_motif_length(1).expect("0,0158") / 3.0);
        assert!(short < none.probability_for(1, 11).expect("sin multiplicador") / 3.0);
    }

    /// **Con un multiplicador cuyo producto tiene muchos decimales, la
    /// fracción publicada es el producto en coma flotante, bit a bit.**
    ///
    /// `k = 2,82` es la separación que la figura de Fungtammasan y col. da a
    /// 15 pb, el ancla. `2,82 × 0,027 = 0,07614`, con cinco decimales: **debe
    /// fallar contra** una entrega que redondee la fracción con PCR —a cuatro
    /// decimales daría 0,0761—, que con 2,5, 1 o 37 no se distingue porque sus
    /// productos son exactos con pocos decimales. La probabilidad esperada está
    /// calculada a mano: `r_pcr = 1 − (1 − 0,07614)/0,973 ≈ 0,0505036` y
    /// `p = 1 − (1 − r_pcr)^(1/10) ≈ 0,00516895`.
    #[test]
    fn a_fold_with_a_long_product_publishes_the_exact_product() {
        let settings =
            ClonalStutterSettings::new(10, &BTreeMap::new(), Some(2.82)).expect("con 2,82");
        let document = effective_settings_document(&settings);
        let published = document["per_cycle_slip"][0]["with_pcr_fraction"]
            .as_f64()
            .expect("fracción");
        assert_eq!(published.to_bits(), (2.82_f64 * 0.027).to_bits());
        assert_ne!(published, 0.0761, "la variante redondeada es distinguible");
        let probability = settings
            .probability_for_motif_length(1)
            .expect("homopolímero");
        assert!(
            (probability - 0.005_168_95).abs() < 1e-8,
            "probabilidad publicada con 2,82: {probability}"
        );
        assert_eq!(document["per_cycle_slip"][0]["probability"], probability);
        assert_eq!(document["homopolymer_pcr_fold"], 2.82);
    }

    /// **Sin multiplicador se usa la constante 0,17 tal cual, no 6,30 × 0,027.**
    ///
    /// Debe fallar contra una entrega que implemente la ausencia como el
    /// multiplicador de Gymrek redondeado: `6,30 × 0,027 = 0,1701`, que mueve la
    /// probabilidad publicada y con ella los bytes.
    #[test]
    fn without_a_fold_the_with_pcr_constant_is_used_verbatim() {
        let none =
            ClonalStutterSettings::new(10, &BTreeMap::new(), None).expect("sin multiplicador");
        let document = effective_settings_document(&none);
        assert_eq!(document["per_cycle_slip"][0]["with_pcr_fraction"], 0.17);
        assert_ne!(
            6.30_f64 * 0.027,
            0.17,
            "la variante defectuosa es distinguible"
        );
        assert_eq!(
            none.probability_for_motif_length(1).map(f64::to_bits),
            Some(per_cycle_slip_probability(0.17, 0.027, 10).to_bits())
        );
        assert!(
            document.get("homopolymer_pcr_fold").is_none(),
            "sin declarar, el campo no se publica"
        );
        // Control positivo: declarado, sí se publica.
        let fold = ClonalStutterSettings::new(10, &BTreeMap::new(), Some(2.5)).expect("con 2,5");
        let document = effective_settings_document(&fold);
        assert_eq!(document["homopolymer_pcr_fold"], 2.5);
        assert_eq!(document["per_cycle_slip"][0]["with_pcr_fraction"], 0.0675);
        assert_eq!(document["per_cycle_slip"][0]["without_pcr_fraction"], 0.027);
        assert_eq!(document["per_cycle_slip"][0]["source"], "derived");
        assert_eq!(document["per_cycle_slip"][1]["with_pcr_fraction"], 0.038);
    }

    /// **El dominio es `[1, 37]`, cerrado por los dos lados, y finito.**
    ///
    /// Debe fallar contra una cota inferior estricta (rechazaría 1), contra
    /// una superior ausente (admitiría 37,01) y contra una comprobación escrita
    /// con comparaciones negadas, `!(k < 1 || k > 37)`, que admitiría `NaN`.
    #[test]
    fn the_fold_domain_is_closed_between_one_and_thirty_seven() {
        for admitted in [1.0, 2.5, 37.0] {
            assert_eq!(validate_homopolymer_pcr_fold(admitted), Ok(admitted));
        }
        for rejected in [0.99, 0.0, -2.5, 37.01, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                validate_homopolymer_pcr_fold(rejected),
                Err(ClonalStutterError::InvalidHomopolymerFold { value: rejected }),
                "{rejected} tiene que rechazarse"
            );
        }
        assert!(matches!(
            validate_homopolymer_pcr_fold(f64::NAN),
            Err(ClonalStutterError::InvalidHomopolymerFold { .. })
        ));
        // Y el constructor lo aplica: no se puede saltar la guarda.
        assert!(ClonalStutterSettings::new(10, &BTreeMap::new(), Some(37.5)).is_err());
    }

    /// **El tope 37 se deriva de la constante sin PCR del homopolímero.**
    ///
    /// Es el mayor entero `k` con `k × 0,027 < 1`. Si la constante cambiara,
    /// esta prueba obliga a recalcular el tope en vez de dejarlo desfasado.
    #[test]
    fn the_fold_maximum_is_the_largest_integer_keeping_the_fraction_below_one() {
        let without = STUTTER_DEFAULT_MOTIF_LENGTH_RATES[0];
        assert!(HOMOPOLYMER_PCR_FOLD_MAXIMUM * without < 1.0);
        assert!((HOMOPOLYMER_PCR_FOLD_MAXIMUM + 1.0) * without >= 1.0);
        assert_eq!(HOMOPOLYMER_PCR_FOLD_MAXIMUM.fract(), 0.0);
        assert_eq!(HOMOPOLYMER_PCR_FOLD_MINIMUM, 1.0);
    }

    /// **Con multiplicador 1 la PCR no añade tartamudeo al homopolímero.**
    ///
    /// La probabilidad vale exactamente cero en los dos caminos, y por tanto el
    /// sorteo no consume (`bernoulli-v1` con cero no lee el generador).
    #[test]
    fn a_unit_fold_gives_a_zero_homopolymer_probability() {
        let unit = ClonalStutterSettings::new(10, &BTreeMap::new(), Some(1.0)).expect("con 1");
        assert_eq!(unit.probability_for_motif_length(1), Some(0.0));
        assert_eq!(unit.probability_for(1, 11), Some(0.0));
        // Control positivo: el dinucleótido sigue tartamudeando.
        assert!(unit.probability_for_motif_length(2).expect("dinucleótido") > 0.0);
    }

    /// **El multiplicador y la sustitución directa del homopolímero se
    /// excluyen; con la de otra longitud conviven.**
    ///
    /// Debe fallar contra una comprobación que mire cualquier sustitución y no
    /// solo la de longitud 1.
    #[test]
    fn the_fold_excludes_only_a_homopolymer_substitution() {
        let mut homopolymer = BTreeMap::new();
        homopolymer.insert(1_u64, 0.01);
        assert_eq!(
            ClonalStutterSettings::new(10, &homopolymer, Some(2.5)),
            Err(ClonalStutterError::HomopolymerFoldWithSubstitution { value: 2.5 })
        );
        let mut dinucleotide = BTreeMap::new();
        dinucleotide.insert(2_u64, 0.01);
        let settings = ClonalStutterSettings::new(10, &dinucleotide, Some(2.5)).expect("conviven");
        assert_eq!(settings.probability_for_motif_length(2), Some(0.01));
        let document = effective_settings_document(&settings);
        assert_eq!(document["per_cycle_slip"][1]["source"], "declared");
        assert_eq!(document["per_cycle_slip"][0]["with_pcr_fraction"], 0.0675);
    }

    /// **La composición sigue cerrando EXACTAMENTE a cualquier longitud de
    /// tracto, no solo en la de referencia.**
    ///
    /// # Qué propiedad se comprueba y por qué importa
    ///
    /// El modelo se apoya en que la parte de base y la de PCR se compongan
    /// devolviendo la fracción medida: `1 − (1−r_sin)(1−r_pcr) = r_con`. Al
    /// añadir la dependencia con la longitud del tracto esa identidad podía
    /// romperse en todas las longitudes salvo una. No se rompe **porque las
    /// dos columnas se modulan con el mismo factor**, y eso es exactamente lo
    /// que esta prueba fija.
    ///
    /// **Debe fallar contra** una entrega que module solo la columna de base,
    /// solo la de con PCR, o que use factores distintos para cada una: en los
    /// tres casos la recomposición se separa de `r_con · m`.
    ///
    /// La tolerancia es de **redondeo puro**, `1e-9` relativo, y no de modelo:
    /// la identidad pasa dos veces por la potencia de `libm` —una al repartir
    /// el incremento entre los ciclos y otra al volver a acumularlo—, y ese
    /// viaje de ida y vuelta deja un residuo medido de unas `5e-11`. No es
    /// holgura para el modelo: las entregas que esta prueba caza se separan por
    /// órdenes de magnitud, no por decimales.
    #[test]
    fn the_composition_closes_at_every_tract_length_not_only_at_the_reference() {
        let cycles = 10_u32;
        let settings =
            ClonalStutterSettings::new(cycles, &BTreeMap::new(), None).expect("parámetros");
        let mut checked = 0_u32;
        let mut moved = 0_u32;
        let mut worst = 0.0_f64;
        for (motif_length, units_range) in [(1_u64, 3..=18_u64), (2, 3..=20)] {
            for units in units_range {
                let index = motif_length as usize - 1;
                let factor = crate::stutter::tract_length_factor(motif_length, units);
                let modulated_with = STUTTER_WITH_PCR_MOTIF_LENGTH_RATES[index] * factor;
                let modulated_without = STUTTER_DEFAULT_MOTIF_LENGTH_RATES[index] * factor;
                let per_cycle = settings
                    .probability_for(motif_length, units)
                    .expect("probabilidad por ciclo");
                // Acumular la probabilidad por ciclo sobre los ciclos devuelve
                // el incremento propio de la PCR...
                let accumulated = 1.0 - libm::pow(1.0 - per_cycle, f64::from(cycles));
                // ...y componerlo con la parte de base tiene que devolver la
                // fracción medida, ya modulada.
                let recomposed = 1.0 - (1.0 - modulated_without) * (1.0 - accumulated);
                let relative = (recomposed - modulated_with).abs() / modulated_with;
                assert!(
                    relative <= 1e-9,
                    "motivo de {motif_length} pb, {units} unidades: recompone {recomposed} \
                     y lo medido modulado es {modulated_with}, error relativo {relative}"
                );
                worst = worst.max(relative);
                checked += 1;
                if factor < 1.0 {
                    moved += 1;
                }
            }
        }
        // **Control positivo.** Si la modulación no moviera nada, la identidad
        // se cumpliría sola y esta prueba no distinguiría ninguna entrega.
        assert!(checked >= 30, "el barrido tiene que ser amplio: {checked}");
        assert!(
            moved >= 20,
            "la mayoría de las longitudes tienen que estar moduladas de verdad: {moved}"
        );
        // Y la tolerancia no está tapando nada: el peor error medido queda muy
        // por debajo de ella, y las entregas que esta prueba tiene que cazar se
        // separan por **factores**, no por decimales. Modular solo la columna
        // con PCR, por ejemplo, la deja por debajo de la de sin PCR en los
        // tractos cortos, con lo que el incremento se anula y la recomposición
        // devuelve **1,43·10⁻⁶** donde tendría que devolver **9,03·10⁻⁶**: un
        // **factor 6,3**, siete órdenes de magnitud por encima de la
        // tolerancia. (Las cifras de este comentario están medidas; una
        // versión anterior decía «0,027 frente a 0,0000090, un factor 3000» y
        // era falsa, porque olvidaba que la recomposición usa la columna sin
        // PCR **ya modulada**.)
        assert!(
            worst < 1e-9,
            "el peor error relativo medido es {worst}, que ya no es redondeo"
        );
    }

    /// **La sustitución de la probabilidad por ciclo PREVALECE y NO se
    /// modula.**
    ///
    /// Es la gemela clonal de la prueba del camino de base. Existía la del
    /// camino de base y esta no, de modo que la prevalencia del lado clonal
    /// solo moría **de rebote**, y una entrega que modulara la sustitución
    /// declarada podía pasar sin que ninguna prueba la nombrara.
    ///
    /// **Debe fallar contra** una entrega que aplique la curva también a la
    /// probabilidad declarada: con cinco unidades de dinucleótido el factor
    /// vale 0,0038 y el 0,5 declarado se quedaría en 0,0019.
    #[test]
    fn a_declared_per_cycle_probability_prevails_and_is_not_modulated() {
        let mut declared = BTreeMap::new();
        declared.insert(2_u64, 0.5);
        let settings = ClonalStutterSettings::new(10, &declared, None).expect("parámetros");
        for units in [1_u64, 5, 13, 17, 40] {
            assert_eq!(
                settings.probability_for(2, units).expect("probabilidad"),
                0.5,
                "la sustitución vale lo mismo mida lo que mida el tracto"
            );
        }
        // **Control positivo**: la longitud de motivo **no** sustituida sí se
        // modula, y a cinco unidades queda muy por debajo de su valor en el
        // ancla. Sin este control, una entrega que no modulara nada pasaría el
        // bucle de arriba.
        let anchor = settings.probability_for(1, 15).expect("en el ancla");
        let short = settings.probability_for(1, 5).expect("tracto corto");
        assert!(
            short < anchor / 100.0,
            "la longitud de motivo sin sustituir sí tiene que modularse:              {short} frente a {anchor}"
        );
    }

    /// **F4.10.3 (a). Las probabilidades por ciclo coinciden a tres cifras con
    /// el cálculo hecho a mano desde las dos columnas.**
    ///
    /// Las seis cifras de la derecha están calculadas **fuera** de este código,
    /// con la fórmula escrita a mano, y son las que el diseño publica. Compararlas
    /// con lo que devolviera otra llamada a la misma función sería un espejo.
    #[test]
    fn the_per_cycle_probabilities_match_the_hand_computation_to_three_figures() {
        let expected = [0.0158, 0.00276, 0.000815, 0.000532, 0.000593, 0.00131];
        let with_pcr = [0.17, 0.038, 0.011, 0.0069, 0.0074, 0.0140];
        let without_pcr = [0.027, 0.011, 0.0029, 0.0016, 0.0015, 0.0010];
        for index in 0..6 {
            let derived = per_cycle_slip_probability(with_pcr[index], without_pcr[index], 10);
            let relative = (derived - expected[index]).abs() / expected[index];
            assert!(
                relative < 5e-3,
                "motivo de {} pb: derivado {derived}, a mano {}",
                index + 1,
                expected[index]
            );
        }
    }

    /// **F4.10.3 (b). La composición devuelve exactamente la fracción medida.**
    ///
    /// Es la restricción aritmética que ordena el plan entero: acumular `p`
    /// sobre los ciclos da `r_pcr`, y componerlo con la parte de base devuelve
    /// `r_con`. Se comprueba a partir de la `p` derivada, **no** volviendo a
    /// despejarla, que sería la identidad escrita al revés.
    #[test]
    fn composing_the_two_parts_returns_the_measured_with_pcr_fraction() {
        let with_pcr = [0.17, 0.038, 0.011, 0.0069, 0.0074, 0.0140];
        let without_pcr = [0.027, 0.011, 0.0029, 0.0016, 0.0015, 0.0010];
        for cycles in [1_u32, 7, 10, 30] {
            for index in 0..6 {
                let probability =
                    per_cycle_slip_probability(with_pcr[index], without_pcr[index], cycles);
                let accumulated = 1.0 - (1.0 - probability).powi(cycles as i32);
                let composed = 1.0 - (1.0 - without_pcr[index]) * (1.0 - accumulated);
                assert!(
                    (composed - with_pcr[index]).abs() < 1e-9,
                    "ciclos {cycles}, motivo de {} pb: compuesto {composed}, medido {}",
                    index + 1,
                    with_pcr[index]
                );
            }
        }
    }

    /// **F4.10.3 (c). Sumar las dos columnas da un resultado distinto.**
    ///
    /// Es el control que distingue la resta de la suma. Sin él, «se compone
    /// bien» no significaría nada: la prueba anterior pasaría igual si la
    /// derivación fuera cualquier otra cosa que compusiera consigo misma.
    #[test]
    fn adding_the_two_columns_instead_of_subtracting_them_overshoots_the_measurement() {
        for index in 0..6 {
            let with_pcr = [0.17, 0.038, 0.011, 0.0069, 0.0074, 0.0140][index];
            let without_pcr = [0.027, 0.011, 0.0029, 0.0016, 0.0015, 0.0010][index];
            let subtracted = 1.0 - (1.0 - with_pcr) / (1.0 - without_pcr);
            let added = with_pcr;
            let composed_from_added = 1.0 - (1.0 - without_pcr) * (1.0 - added);
            assert!(subtracted < added);
            assert!(
                composed_from_added > with_pcr,
                "usar la columna con PCR como si fuera el incremento sobrepasa lo medido"
            );
        }
    }

    /// Con la fracción con PCR igual o menor que la de sin PCR, la
    /// probabilidad vale cero.
    #[test]
    fn a_with_pcr_fraction_that_does_not_exceed_the_other_gives_zero() {
        assert_eq!(per_cycle_slip_probability(0.02, 0.02, 10), 0.0);
        assert_eq!(per_cycle_slip_probability(0.01, 0.02, 10), 0.0);
        // Y el control positivo: cuando sí la supera, no vale cero.
        assert!(per_cycle_slip_probability(0.17, 0.027, 10) > 0.0);
    }

    /// **La sustitución directa del usuario prevalece y se marca como suya.**
    #[test]
    fn a_declared_substitution_wins_and_is_published_as_the_users() {
        let mut substitutions = BTreeMap::new();
        substitutions.insert(1_u64, 0.5_f64);
        let settings = ClonalStutterSettings::new(10, &substitutions, None).expect("parámetros");
        assert_eq!(settings.probability_for_motif_length(1), Some(0.5));
        let document = effective_settings_document(&settings);
        let entries = document["per_cycle_slip"].as_array().expect("tabla");
        assert_eq!(entries[0]["source"], "declared");
        assert_eq!(entries[1]["source"], "derived");
        // Y la derivada no se ha movido: sustituir una no toca las demás.
        let derived = ClonalStutterSettings::new(10, &BTreeMap::new(), None).expect("parámetros");
        assert_eq!(
            settings.probability_for_motif_length(2),
            derived.probability_for_motif_length(2)
        );
        // La marca de tosquedad viaja siempre.
        assert_eq!(document["coarse"], true);
    }

    /// **Los pesos de la divergencia son la geométrica pinzada en el último
    /// ciclo, y suman uno.**
    ///
    /// La masa del último no es el término que le tocaría, sino **toda** la
    /// cola: `q^c`. Una implementación que pusiera ahí `(1−q)q^c` perdería
    /// masa, y esta prueba lo caza por la suma.
    #[test]
    fn the_divergence_weights_are_geometric_and_clamped_at_the_last_cycle() {
        let weights = divergence_weights(10, 0.9);
        assert_eq!(weights.len(), 11);
        let total: f64 = weights.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "suman {total}");
        // `P(d=0) = e/(1+e)`, escrito a mano.
        assert!((weights[0] - 0.9 / 1.9).abs() < 1e-12);
        // Y decrecen: `q = 1/1,9 < 1`.
        for pair in weights[..10].windows(2) {
            assert!(pair[1] < pair[0]);
        }
        // La cola pinzada vale `q^c`, no el término geométrico.
        let clamped = (1.0_f64 / 1.9).powi(10);
        assert!((weights[10] - clamped).abs() < 1e-12);
    }
}
