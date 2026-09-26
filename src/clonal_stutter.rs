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
//! entre un 14 % y un 28 % en motivos de uno a cinco pares de base. La
//! admisión de la configuración lo impide.
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
pub const CLONAL_STUTTER_ALGORITHM: &str = "clonal-stutter-v1";

/// La derivación, escrita tal cual para que viaje con los números.
pub const CLONAL_STUTTER_FORMULA: &str =
    "r_pcr = 1 - (1 - r_with_pcr) / (1 - r_without_pcr); p_cycle = 1 - (1 - r_pcr)^(1/cycles)";

/// **Marca de tosquedad**, publicada siempre que estas tasas se usen.
pub const CLONAL_STUTTER_COARSE: bool = true;

/// Por qué son toscas, dicho entero para que nadie tenga que buscarlo.
pub const CLONAL_STUTTER_COARSE_REASON: &str = "the two calibration columns come from two different individuals, from different projects and \
     at different coverages, so the subtraction may carry individual and batch effects besides \
     the amplification: these per-cycle rates are coarse and are published as such";

/// Fallo estable de la sustitución directa de las tasas por ciclo.
#[derive(Clone, Debug, PartialEq)]
pub enum ClonalStutterError {
    /// Longitud de motivo fuera del tramo calibrado.
    UnknownMotifLength { motif_length: u64 },
    /// Probabilidad fuera de `[0, 1]` o no finita.
    InvalidProbability { motif_length: u64, value: f64 },
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
    pub fn new(
        cycles: u32,
        substitutions: &BTreeMap<u64, f64>,
    ) -> Result<Self, ClonalStutterError> {
        let mut per_cycle = [0.0_f64; STUTTER_MOTIF_LENGTH_MAX];
        let mut source = [PerCycleSlipSource::Derived; STUTTER_MOTIF_LENGTH_MAX];
        for (index, slot) in per_cycle.iter_mut().enumerate() {
            *slot = per_cycle_slip_probability(
                STUTTER_WITH_PCR_MOTIF_LENGTH_RATES[index],
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
            per_cycle,
            source,
        })
    }

    #[must_use]
    pub const fn cycles(&self) -> u32 {
        self.cycles
    }

    /// Probabilidad por ciclo de una longitud de motivo, o `None` fuera del
    /// tramo calibrado. No se extrapola.
    #[must_use]
    pub fn probability_for_motif_length(&self, motif_length: u64) -> Option<f64> {
        usize::try_from(motif_length)
            .ok()
            .filter(|length| (1..=STUTTER_MOTIF_LENGTH_MAX).contains(length))
            .map(|length| self.per_cycle[length - 1])
    }
}

/// Documento común de los parámetros efectivos, que publican **la
/// configuración efectiva y el manifiesto**.
///
/// Sale de aquí para que no puedan discrepar, y lleva **las dos fracciones de
/// entrada, el número de ciclos, la fórmula y la marca de tosquedad**: unos
/// números sin su procedencia no se distinguen de unos inventados, y una tasa
/// tosca sin su marca se leería como una medida.
#[must_use]
pub fn effective_settings_document(settings: &ClonalStutterSettings) -> serde_json::Value {
    let per_cycle_slip = (1..=STUTTER_MOTIF_LENGTH_MAX)
        .map(|motif_length| {
            let index = motif_length - 1;
            serde_json::json!({
                "motif_length_bases": motif_length,
                "with_pcr_fraction": STUTTER_WITH_PCR_MOTIF_LENGTH_RATES[index],
                "without_pcr_fraction": STUTTER_DEFAULT_MOTIF_LENGTH_RATES[index],
                "probability": settings.per_cycle[index],
                "source": settings.source[index].as_str(),
            })
        })
        .collect::<Vec<_>>();
    serde_json::json!({
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
    })
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
        ClonalStutterSettings, divergence_weights, effective_settings_document,
        per_cycle_slip_probability,
    };
    use std::collections::BTreeMap;

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
        let settings = ClonalStutterSettings::new(10, &substitutions).expect("parámetros");
        assert_eq!(settings.probability_for_motif_length(1), Some(0.5));
        let document = effective_settings_document(&settings);
        let entries = document["per_cycle_slip"].as_array().expect("tabla");
        assert_eq!(entries[0]["source"], "declared");
        assert_eq!(entries[1]["source"], "derived");
        // Y la derivada no se ha movido: sustituir una no toca las demás.
        let derived = ClonalStutterSettings::new(10, &BTreeMap::new()).expect("parámetros");
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
