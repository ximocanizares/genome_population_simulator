//! F5.3.1 y F5.3.2. Modo por bloque×subpoblación, cuotas controladas,
//! muestreo probabilístico de las copias fundadoras, F diploide y F realizado.
//!
//! DTC-089 «Fundadores, F y factibilidad controlada» fija cinco cosas que este
//! módulo implementa y ninguna más:
//!
//! 1. **Resolución de modo.** Cada bloque×subpoblación resuelve su modo
//!    solicitado y efectivo —el predeterminado de la subpoblación o la
//!    sustitución única por `population_id`×`block_index`— y, en controlado, su
//!    tolerancia declarada y efectiva, cuyo predeterminado es `0.5/C`.
//! 2. **Modo probabilístico.** Una fila por copia muestreada con `q`, cada una
//!    en su propio ámbito `founder_sample`.
//! 3. **Modo controlado.** Cuotas enteras por mayores restos con comprobación
//!    de factibilidad `O(H)` **antes** de abrir ningún RNG, desempate único por
//!    Fisher–Yates del grupo fronterizo y permutación Fisher–Yates de las `C`
//!    identidades de copia que asigna el multiconjunto ordenado por
//!    `haplotype_rank`.
//! 4. **F diploide.** En `P=2` la copia 1 reutiliza con probabilidad `F`
//!    exactamente la fila de la copia 0 del mismo individuo; `founder_f` se
//!    consume **antes** que `founder_sample` de la copia 1 y, si decide
//!    reutilizar, ese segundo flujo no se abre. `F=0` y `F=1` no consumen
//!    `founder_f`. `F` declarado con `P=1` o `P>2` es error.
//! 5. **F realizado.** Por población×bloque se publica `F_realized`,
//!    `not_applicable` fuera de `P=2` y cuando `H_expected=0`.
//!
//! El Ne reproductivo de DTC-089 es `not_applicable` en F5 **por
//! construcción**: F5 no reproduce, aquí no se calcula ninguno y ningún
//! esquema publicado tiene dónde escribirlo. Una búsqueda negativa sobre el
//! módulo impide que la fórmula por contribuciones gaméticas de F7.3 aparezca.
//!
//! F5.3.3 añade la sexta pieza, que no decide nada y solo materializa: las
//! identidades de fundador y su pertenencia única —`IndividualId`, `FCP_` por
//! secuencia y copia, `HAP_` por fila haplotípica y `FHA_` por bloque como
//! asignación separada— derivadas con DTC-039 **antes de toda salida**.
//!
//! Lo que **no** hace, por pertenecer a otros hitos: la fórmula de Ne por
//! contribuciones gaméticas es de F7.3, y escribir los siete TSV y la
//! proyección de verdad es de F5.4.

use std::collections::{BTreeMap, BTreeSet};

use crate::application::EffectiveF5Config;
use crate::configuration::{DeclaredF5SamplingOverride, DeclaredF5Subpopulation};
use crate::diagnostics::Diagnostic;
use crate::f5_haplotypes::{F5AlleleMatrix, F5BlockHaplotypes};
use crate::f5_layout::F5PhysicalLayout;
use crate::f5_marginals::NeumaierSum;
use crate::f5_preflight::{F5OtherResources, SealedF5Admission};
use crate::f5_search::F5HaplotypeSearch;
use crate::identity::{
    AutomaticIdRequest, ExecutionId, F5BlockHaplotypeAssignmentScope, F5BlockId,
    F5FounderCopyScope, F5HaplotypeId, F5IdentityContext, F5IdentityRegistry,
    F5PopulationIdentityContext, FounderCopyId, FounderHaplotypeAssignmentId, FounderKey, IdClass,
    IdentityRegistry, IndividualId, PopulationId, ReferenceId, ScenarioId, SequenceId,
};
use crate::randomness::{
    F5RandomDecision, LocalRng, RandomWordSource, RootSeed, derive_f5_subseed,
    f5_bernoulli_with_source, f5_fisher_yates_with_source, f5_uniform_53,
};

// ---------------------------------------------------------------------------
// Tipos publicados
// ---------------------------------------------------------------------------

/// Modo de muestreo de una unidad bloque×subpoblación.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum F5SamplingMode {
    Probabilistic,
    Controlled,
}

impl F5SamplingMode {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Probabilistic => "probabilistic",
            Self::Controlled => "controlled",
        }
    }

    /// Los dos únicos literales que el esquema de escenario admite.
    fn parse(text: &str) -> Option<Self> {
        match text {
            "probabilistic" => Some(Self::Probabilistic),
            "controlled" => Some(Self::Controlled),
            _ => None,
        }
    }
}

/// Origen del modo de una unidad: el predeterminado de su subpoblación o la
/// sustitución única de esa pareja `population_id`×`block_index`.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum F5SamplingSource {
    SubpopulationDefault,
    BlockOverride,
}

impl F5SamplingSource {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::SubpopulationDefault => "subpopulation_default",
            Self::BlockOverride => "block_override",
        }
    }
}

/// Una fila exacta de la lista `sampling` del efectivo y del manifiesto.
///
/// `requested_mode` es el modo que el documento pide **para esta unidad** —el
/// de la sustitución si existe, y si no el predeterminado de la
/// subpoblación—; `source` dice de cuál de los dos vino. `effective_mode`
/// coincide siempre con él: DTC-089 prohíbe corregir o relajar una cuota, de
/// modo que una unidad inviable se rechaza y no se degrada a probabilístico.
#[derive(Clone, Debug, PartialEq)]
pub struct F5SamplingResolution {
    population_id: PopulationId,
    block_id: F5BlockId,
    block_index: u64,
    requested_mode: F5SamplingMode,
    effective_mode: F5SamplingMode,
    source: F5SamplingSource,
    controlled_max_abs_error_declared: Option<f64>,
    controlled_max_abs_error_effective: Option<f64>,
    /// La misma tolerancia efectiva medida **en copias**, que es donde la
    /// comparación es exacta: `0.5` exacto por omisión y `t * C` si se declaró.
    controlled_max_abs_error_copies: Option<f64>,
    copies: u64,
}

impl F5SamplingResolution {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    #[must_use]
    pub const fn block_id(&self) -> &F5BlockId {
        &self.block_id
    }
    #[must_use]
    pub const fn block_index(&self) -> u64 {
        self.block_index
    }
    #[must_use]
    pub const fn requested_mode(&self) -> F5SamplingMode {
        self.requested_mode
    }
    #[must_use]
    pub const fn effective_mode(&self) -> F5SamplingMode {
        self.effective_mode
    }
    #[must_use]
    pub const fn source(&self) -> F5SamplingSource {
        self.source
    }
    /// Tolerancia que la sustitución declaró; `None` si no la declaró o si la
    /// unidad es probabilística, donde el contrato la prohíbe.
    #[must_use]
    pub const fn controlled_max_abs_error_declared(&self) -> Option<f64> {
        self.controlled_max_abs_error_declared
    }
    /// Tolerancia aplicada: la declarada o `0.5/C`; `None` en probabilístico.
    #[must_use]
    pub const fn controlled_max_abs_error_effective(&self) -> Option<f64> {
        self.controlled_max_abs_error_effective
    }
    /// `C = individual_count * ploidy` de la subpoblación de esta fila.
    #[must_use]
    pub const fn copies(&self) -> u64 {
        self.copies
    }
}

/// Un individuo fundador declarado: su identidad y, si vino de una clave, la
/// clave que la derivó.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5FounderIndividual {
    individual_id: IndividualId,
    founder_key: Option<FounderKey>,
}

impl F5FounderIndividual {
    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }
    #[must_use]
    pub const fn founder_key(&self) -> Option<&FounderKey> {
        self.founder_key.as_ref()
    }
}

/// Censo fundador de una subpoblación, **ordenado por `individual_id`**.
///
/// Es el orden que DTC-089 exige a la lista `founders` del efectivo, y el que
/// F5.3.3 y F5.4 heredan al reutilizar este censo. No cambia ningún sorteo:
/// las identidades de copia se ordenan igualmente antes de permutarse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5FounderRoster {
    population_id: PopulationId,
    individual_count: u64,
    copies: u64,
    individuals: Vec<F5FounderIndividual>,
}

impl F5FounderRoster {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    #[must_use]
    pub const fn individual_count(&self) -> u64 {
        self.individual_count
    }
    /// `C = individual_count * ploidy`.
    #[must_use]
    pub const fn copies(&self) -> u64 {
        self.copies
    }
    #[must_use]
    pub fn individuals(&self) -> &[F5FounderIndividual] {
        &self.individuals
    }
}

/// Asignación de una copia: `haplotype_rank` de la fila que recibe.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5CopyAssignment {
    individual_id: IndividualId,
    copy_index: u32,
    haplotype_rank: u64,
}

impl F5CopyAssignment {
    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }
    #[must_use]
    pub const fn copy_index(&self) -> u32 {
        self.copy_index
    }
    #[must_use]
    pub const fn haplotype_rank(&self) -> u64 {
        self.haplotype_rank
    }
}

/// Frontera empatada de una cuota controlada.
///
/// Es el objeto de `quota_ties` del manifiesto: `candidate_ranks` junto a
/// `select_count` representan **todas** las combinaciones óptimas sin
/// enumerarlas, y `selected_ranks` conserva el prefijo en el orden sorteado.
#[derive(Clone, Debug, PartialEq)]
pub struct F5QuotaTie {
    remainder: f64,
    candidate_ranks: Vec<u64>,
    select_count: u64,
    selected_ranks: Vec<u64>,
}

impl F5QuotaTie {
    #[must_use]
    pub const fn remainder(&self) -> f64 {
        self.remainder
    }
    /// Rangos empatados, en orden ascendente.
    #[must_use]
    pub fn candidate_ranks(&self) -> &[u64] {
        &self.candidate_ranks
    }
    #[must_use]
    pub const fn select_count(&self) -> u64 {
        self.select_count
    }
    /// Prefijo elegido, en el orden que produjo la Fisher–Yates del grupo.
    #[must_use]
    pub fn selected_ranks(&self) -> &[u64] {
        &self.selected_ranks
    }
}

/// Cuota controlada de una unidad: recuentos por rango y su frontera.
#[derive(Clone, Debug, PartialEq)]
pub struct F5ControlledQuota {
    copies: u64,
    counts: Vec<u64>,
    max_abs_error_effective: f64,
    tie: Option<F5QuotaTie>,
}

impl F5ControlledQuota {
    #[must_use]
    pub const fn copies(&self) -> u64 {
        self.copies
    }
    /// `n_h` por `haplotype_rank`; su suma es exactamente `C`.
    #[must_use]
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }
    #[must_use]
    pub const fn max_abs_error_effective(&self) -> f64 {
        self.max_abs_error_effective
    }
    /// `None` cuando ningún grupo de restos iguales cruzó la frontera y, por
    /// tanto, el desempate no consumió RNG.
    #[must_use]
    pub const fn tie(&self) -> Option<&F5QuotaTie> {
        self.tie.as_ref()
    }
}

/// Estado de un `F_realized` publicado: evaluado o no aplicable.
///
/// No existe un tercer estado y `not_applicable` **nunca** se sustituye por
/// cero: REQ-EVO-52 lo prohíbe expresamente.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum F5RealizedFStatus {
    Evaluated,
    NotApplicable,
}

impl F5RealizedFStatus {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Evaluated => "evaluated",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// `F_realized` de una población×bloque, con la forma exacta que el manifiesto
/// publica: `target_f`, `status` y `value`.
///
/// DTC-089: `F_realized = 1 - H_observed/H_expected`, con
/// `H_observed = heterozygous_individuals/individual_count`,
/// `q_realized_h = n_h/C` y `H_expected = 1 - Σ_h q_realized_h²`. La
/// frecuencia es la **realizada**, no el objetivo `q` del perfil. Un valor
/// negativo es un exceso de heterocigotos legítimo y se conserva (REQ-EVO-52).
#[derive(Clone, Debug, PartialEq)]
pub struct F5RealizedF {
    population_id: PopulationId,
    block_id: F5BlockId,
    target_f: Option<f64>,
    status: F5RealizedFStatus,
    value: Option<f64>,
}

impl F5RealizedF {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    #[must_use]
    pub const fn block_id(&self) -> &F5BlockId {
        &self.block_id
    }
    /// F efectivo de la subpoblación —el declarado o `0.0`— en `P=2`, y `None`
    /// en cualquier otra ploidía, donde declarar F es error.
    #[must_use]
    pub const fn target_f(&self) -> Option<f64> {
        self.target_f
    }
    #[must_use]
    pub const fn status(&self) -> F5RealizedFStatus {
        self.status
    }
    /// `Some` si y solo si el estado es `evaluated`.
    #[must_use]
    pub const fn value(&self) -> Option<f64> {
        self.value
    }
}

/// Unidad bloque×subpoblación resuelta y asignada.
#[derive(Clone, Debug, PartialEq)]
pub struct F5FounderUnit {
    population_id: PopulationId,
    block_id: F5BlockId,
    block_index: u64,
    mode: F5SamplingMode,
    assignments: Vec<F5CopyAssignment>,
    quota: Option<F5ControlledQuota>,
    operations: u64,
}

impl F5FounderUnit {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    #[must_use]
    pub const fn block_id(&self) -> &F5BlockId {
        &self.block_id
    }
    #[must_use]
    pub const fn block_index(&self) -> u64 {
        self.block_index
    }
    #[must_use]
    pub const fn mode(&self) -> F5SamplingMode {
        self.mode
    }
    /// Asignaciones ordenadas por `individual_id` y `copy_index`, que es el
    /// orden inicial de las identidades de copia antes de la permutación.
    #[must_use]
    pub fn assignments(&self) -> &[F5CopyAssignment] {
        &self.assignments
    }
    /// `Some` solo en controlado.
    #[must_use]
    pub const fn quota(&self) -> Option<&F5ControlledQuota> {
        self.quota.as_ref()
    }

    /// Visitas elementales que costó resolver esta unidad, contadas **dentro**
    /// de cada bucle: los rangos recorridos, las comparaciones de la selección
    /// parcial, los sorteos de las dos permutaciones y las copias escritas.
    ///
    /// Es la medida que acredita el coste lineal en `C`: ni el reparto ni la
    /// asignación ordenan nada, de modo que el total no puede crecer como
    /// `C log C`.
    #[must_use]
    pub const fn operations(&self) -> u64 {
        self.operations
    }
}

/// Resultado completo de F5.3.1.
#[derive(Clone, Debug, PartialEq)]
pub struct F5FounderAssignment {
    rosters: Vec<F5FounderRoster>,
    sampling: Vec<F5SamplingResolution>,
    units: Vec<F5FounderUnit>,
    f_realized: Vec<F5RealizedF>,
}

impl F5FounderAssignment {
    /// Censos por subpoblación, ordenados por `population_id`.
    #[must_use]
    pub fn rosters(&self) -> &[F5FounderRoster] {
        &self.rosters
    }
    /// Lista `sampling` ordenada por `population_id` y `block_id`, con una
    /// fila exacta por bloque×subpoblación.
    #[must_use]
    pub fn sampling(&self) -> &[F5SamplingResolution] {
        &self.sampling
    }
    /// Unidades en el mismo orden que [`Self::sampling`].
    #[must_use]
    pub fn units(&self) -> &[F5FounderUnit] {
        &self.units
    }
    /// `F_realized` por población×bloque, en el mismo orden que
    /// [`Self::sampling`] —población y después `block_id`— y con una fila
    /// exacta por unidad, también fuera de `P=2`.
    #[must_use]
    pub fn f_realized(&self) -> &[F5RealizedF] {
        &self.f_realized
    }
    /// Suma de las visitas de todas las unidades.
    #[must_use]
    pub fn operations(&self) -> u64 {
        self.units
            .iter()
            .fold(0_u64, |total, unit| total.saturating_add(unit.operations))
    }

    #[must_use]
    pub fn unit(&self, population_id: &PopulationId, block_index: u64) -> Option<&F5FounderUnit> {
        self.units
            .iter()
            .find(|unit| unit.population_id == *population_id && unit.block_index == block_index)
    }
}

// ---------------------------------------------------------------------------
// Frontera pública
// ---------------------------------------------------------------------------

/// Resuelve modos y tolerancias, asigna una fila haplotípica a cada copia,
/// aplica F en diploides y publica `F_realized`.
///
/// Consume la matriz ganadora de cada bloque×subpoblación que produjo F5.2.3:
/// los rangos que asigna son filas de esa matriz y ningún peso se recalcula
/// aquí. Exige la admisión sellada de la misma configuración porque abre
/// ámbitos RNG.
///
/// # Errors
///
/// - `GPS-CONF-101` si la sustitución nombra una población desconocida, repite
///   una pareja, sale de `0..B-1`, si el modo declarado no es un literal del
///   contrato, si la búsqueda no cubre exactamente las unidades del escenario,
///   si una subpoblación controlada declara `F>0` —combinación que DTC-089
///   declara incompatible— o si una subpoblación declara F con `P=1` o `P>2`.
/// - `GPS-DOMAIN-103` si la cuota controlada no puede satisfacerse dentro de
///   su tolerancia efectiva; siempre **antes** de abrir ningún RNG de cuotas.
/// - `GPS-INTERNAL-101` ante una aritmética no finita o un desbordamiento
///   posterior a las cotas ya demostradas por el preflight.
pub fn assign_f5_founders(
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
    search: &F5HaplotypeSearch,
) -> Result<F5FounderAssignment, Diagnostic> {
    if !sealed.matches_admission(config, F5OtherResources::default()) {
        return Err(configuration(
            config,
            "preflight.binding",
            "la admisión sellada no corresponde a esta configuración",
        ));
    }
    let scenario_id = ScenarioId::parse(config.scenario_id(), "scenario_id")
        .map_err(|_| configuration(config, "scenario_id", "no es una identidad tipada"))?;
    let root_seed = RootSeed::parse("seed", config.seed())
        .map_err(|_| configuration(config, "seed", "la semilla efectiva no es canónica"))?;

    let weights = config.profile().haplotype_weights();
    if weights.is_empty() {
        return Err(configuration(
            config,
            "profile.haplotype_weights",
            "el perfil admitido no conserva ningún peso haplotípico",
        ));
    }
    let haplotype_count = weights.len() as u64;
    let prefix = HaplotypePrefix::new(weights)
        .map_err(|detail| configuration(config, "profile.haplotype_weights", detail))?;

    let rosters = resolve_f5_founder_roster(config)?;
    let blocks = index_search_blocks(config, sealed, search, &rosters, weights.len())?;
    let sampling = resolve_f5_sampling(config, &rosters, &blocks)?;

    // Las `C` identidades de copia de una subpoblación son las mismas en todos
    // sus bloques: se ordenan una vez, no una vez por unidad.
    let ploidy = config.ploidy();
    let mut seats_by_population: BTreeMap<PopulationId, Vec<CopySeat>> = BTreeMap::new();
    // F efectivo de cada subpoblación y su objetivo publicable, resueltos una
    // vez por subpoblación: no dependen del bloque. Aquí se rechaza F con
    // `P=1` o `P>2`, que REQ-EVO-48 exige identificando ámbito y ploidía.
    let mut inbreeding_by_population: BTreeMap<PopulationId, (FounderInbreeding, Option<f64>)> =
        BTreeMap::new();
    for roster in &rosters {
        seats_by_population.insert(roster.population_id.clone(), copy_seats(roster, ploidy));
        let subpopulation = declared_subpopulation(config, &roster.population_id)?;
        let declared = effective_inbreeding_f(config, subpopulation)?;
        let inbreeding = resolve_founder_inbreeding(ploidy, declared).map_err(|detail| {
            configuration(
                config,
                "subpopulations.inbreeding_f,ploidy",
                format!(
                    "la subpoblación {}: {detail}",
                    roster.population_id.as_str()
                ),
            )
        })?;
        // DTC-089 publica `target_f` solo en diploides; en otra ploidía el
        // objetivo es `null`, nunca cero.
        let target_f = (ploidy == 2).then(|| declared.unwrap_or(0.0));
        inbreeding_by_population.insert(roster.population_id.clone(), (inbreeding, target_f));
    }

    let mut units = Vec::with_capacity(sampling.len());
    let mut f_realized = Vec::with_capacity(sampling.len());
    for row in &sampling {
        let seats = seats_by_population.get(&row.population_id).ok_or_else(|| {
            internal(
                &scenario_id,
                "la fila de muestreo nombra una subpoblación que el censo no contiene",
            )
        })?;
        if seats.len() as u64 != row.copies {
            return Err(internal(
                &scenario_id,
                format!(
                    "la subpoblación {} produjo {} identidades de copia y su C declarado es {}",
                    row.population_id.as_str(),
                    seats.len(),
                    row.copies
                ),
            ));
        }
        let (inbreeding, target_f) = inbreeding_by_population
            .get(&row.population_id)
            .copied()
            .ok_or_else(|| {
                internal(
                    &scenario_id,
                    "la fila de muestreo nombra una subpoblación sin F resuelto",
                )
            })?;
        let (assignments, quota, operations) = match row.effective_mode {
            F5SamplingMode::Probabilistic => {
                let (assignments, operations) = sample_probabilistic_unit(
                    &scenario_id,
                    root_seed,
                    row,
                    seats,
                    &prefix,
                    inbreeding,
                )?;
                (assignments, None, operations)
            }
            F5SamplingMode::Controlled => {
                let (assignments, quota, operations) =
                    sample_controlled_unit(&scenario_id, root_seed, row, seats, weights)?;
                (assignments, Some(quota), operations)
            }
        };
        if assignments
            .iter()
            .any(|assignment| assignment.haplotype_rank >= haplotype_count)
        {
            return Err(internal(
                &scenario_id,
                "una asignación nombra un rango haplotípico que la matriz no contiene",
            ));
        }
        let (realized, realized_operations) = realized_f_row(
            &scenario_id,
            row,
            ploidy,
            target_f,
            &assignments,
            weights.len(),
        )?;
        let operations = operations.saturating_add(realized_operations);
        f_realized.push(realized);
        units.push(F5FounderUnit {
            population_id: row.population_id.clone(),
            block_id: row.block_id.clone(),
            block_index: row.block_index,
            mode: row.effective_mode,
            assignments,
            quota,
            operations,
        });
    }

    Ok(F5FounderAssignment {
        rosters,
        sampling,
        units,
        f_realized,
    })
}

/// Censo fundador de cada subpoblación, ordenado por `population_id`.
///
/// Un `founder_id` aportado se usa literalmente como `IndividualId`; una
/// `founder_key` deriva el suyo con el contrato población×clave de REQ-IND-25,
/// el mismo que ya usa la vía F1/F2, sin ordinal físico.
///
/// DTC-089 exige que `IndividualId` identifique al fundador **una sola vez en
/// toda la ejecución**. La admisión ya rechaza un `founder_id` literal repetido
/// globalmente y el registro de identidades rechaza dos claves que deriven el
/// mismo identificador, pero nadie cruzaba las dos mitades: un `founder_id`
/// literal con el texto exacto de un `IND_` derivado de una `founder_key` de
/// otra subpoblación es sintácticamente válido y pasaba la admisión. Aquí se
/// cierra, que es el primer punto donde ambos conjuntos existen a la vez: la
/// admisión trabaja sobre el documento declarado, donde los identificadores
/// derivados todavía no se han construido, y derivarlos allí obligaría a la
/// máquina de admisión a abrir la derivación de identidades con la semilla
/// raíz, cosa que hoy no hace.
///
/// # Errors
///
/// `GPS-CONF-101` si la subpoblación no conserva exactamente un conjunto de
/// fundadores, si un identificador no es válido, si `C` desborda `u64`, si la
/// derivación de una clave colisiona o si un mismo `IndividualId` aparece en
/// más de una subpoblación.
pub fn resolve_f5_founder_roster(
    config: &EffectiveF5Config,
) -> Result<Vec<F5FounderRoster>, Diagnostic> {
    let root_seed = RootSeed::parse("seed", config.seed())
        .map_err(|_| configuration(config, "seed", "la semilla efectiva no es canónica"))?;
    let ploidy = u64::from(config.ploidy());
    let mut registry = IdentityRegistry::default();
    let mut rosters = Vec::with_capacity(config.declared_scenario().subpopulations.len());
    let mut declared_individuals: BTreeSet<IndividualId> = BTreeSet::new();
    for subpopulation in &config.declared_scenario().subpopulations {
        let population_id =
            PopulationId::parse(subpopulation.population_id.value().clone(), "population_id")
                .map_err(|_| {
                    configuration(
                        config,
                        "subpopulations.population_id",
                        "no es una identidad tipada",
                    )
                })?;
        let individual_count = *subpopulation.individual_count.value();
        let copies = individual_count.checked_mul(ploidy).ok_or_else(|| {
            configuration(
                config,
                "subpopulations.individual_count,ploidy",
                "el número de copias desborda u64",
            )
        })?;
        let mut individuals = resolve_founders(
            config,
            root_seed,
            &population_id,
            subpopulation,
            &mut registry,
        )?;
        for individual in &individuals {
            if !declared_individuals.insert(individual.individual_id.clone()) {
                return Err(configuration(
                    config,
                    "subpopulations.founder_ids,subpopulations.founder_keys",
                    format!(
                        "el fundador {} aparece en más de una subpoblación, y DTC-089 lo identifica una sola vez en toda la ejecución",
                        individual.individual_id.as_str()
                    ),
                ));
            }
        }
        individuals.sort_by(|left, right| left.individual_id.cmp(&right.individual_id));
        if individuals.len() as u64 != individual_count {
            return Err(configuration(
                config,
                "subpopulations.founder_ids,subpopulations.founder_keys",
                "la cardinalidad fundadora no coincide con individual_count",
            ));
        }
        rosters.push(F5FounderRoster {
            population_id,
            individual_count,
            copies,
            individuals,
        });
    }
    rosters.sort_by(|left, right| left.population_id.cmp(&right.population_id));
    if rosters
        .windows(2)
        .any(|pair| pair[0].population_id == pair[1].population_id)
    {
        return Err(configuration(
            config,
            "subpopulations.population_id",
            "una subpoblación aparece más de una vez",
        ));
    }
    Ok(rosters)
}

fn resolve_founders(
    config: &EffectiveF5Config,
    root_seed: RootSeed,
    population_id: &PopulationId,
    subpopulation: &DeclaredF5Subpopulation,
    registry: &mut IdentityRegistry,
) -> Result<Vec<F5FounderIndividual>, Diagnostic> {
    match (&subpopulation.founder_ids, &subpopulation.founder_keys) {
        (Some(ids), None) => {
            ids.iter()
                .map(|id| {
                    Ok(F5FounderIndividual {
                        individual_id: IndividualId::parse(id.value().clone(), "founder_ids")
                            .map_err(|_| {
                                configuration(
                                    config,
                                    "subpopulations.founder_ids",
                                    "no es una identidad tipada",
                                )
                            })?,
                        founder_key: None,
                    })
                })
                .collect()
        }
        (None, Some(keys)) => keys
            .iter()
            .map(|key| {
                let founder_key =
                    FounderKey::parse(key.value().clone(), "founder_keys").map_err(|_| {
                        configuration(
                            config,
                            "subpopulations.founder_keys",
                            "no es una clave fundadora válida",
                        )
                    })?;
                let request = AutomaticIdRequest::new(
                    IdClass::Individual,
                    population_id.as_str(),
                    &[("founder_key", founder_key.as_str())],
                )
                .map_err(Diagnostic::into_f5_configuration)?;
                let individual_id = registry
                    .derive_individual(root_seed, &request)
                    .map_err(Diagnostic::into_f5_configuration)?;
                Ok(F5FounderIndividual {
                    individual_id,
                    founder_key: Some(founder_key),
                })
            })
            .collect(),
        _ => Err(configuration(
            config,
            "subpopulations.founder_ids,subpopulations.founder_keys",
            "la subpoblación no conserva exactamente un conjunto de fundadores",
        )),
    }
}

// ---------------------------------------------------------------------------
// Resolución de modo y tolerancia
// ---------------------------------------------------------------------------

/// Bloque de la búsqueda ya validado contra el escenario.
#[derive(Clone, Debug)]
struct UnitBlock {
    population_id: PopulationId,
    block_id: F5BlockId,
    block_index: u64,
}

/// Índice de los bloques ganadores por población, con la comprobación de que
/// la búsqueda cubre exactamente las `B` unidades de cada subpoblación.
fn index_search_blocks(
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
    search: &F5HaplotypeSearch,
    rosters: &[F5FounderRoster],
    haplotype_count: usize,
) -> Result<Vec<UnitBlock>, Diagnostic> {
    let block_count = sealed.block_count();
    let mut seen: BTreeSet<(PopulationId, u64)> = BTreeSet::new();
    let mut blocks = Vec::with_capacity(search.blocks().len());
    for block in search.blocks() {
        if block.matrix().haplotype_count() != haplotype_count {
            return Err(configuration(
                config,
                "search.matrix",
                "la matriz ganadora no tiene tantas filas como pesos el perfil",
            ));
        }
        if block.block_index() >= block_count {
            return Err(configuration(
                config,
                "search.block_index",
                "la búsqueda nombra un bloque fuera del layout sellado",
            ));
        }
        if !rosters
            .iter()
            .any(|roster| roster.population_id == *block.population_id())
        {
            return Err(configuration(
                config,
                "search.population_id",
                "la búsqueda nombra una subpoblación que el escenario no declara",
            ));
        }
        if !seen.insert((block.population_id().clone(), block.block_index())) {
            return Err(configuration(
                config,
                "search.block_index",
                "la búsqueda repite un bloque dentro de la misma subpoblación",
            ));
        }
        blocks.push(UnitBlock {
            population_id: block.population_id().clone(),
            block_id: block.block_id().clone(),
            block_index: block.block_index(),
        });
    }
    let expected = u64::try_from(rosters.len())
        .ok()
        .and_then(|populations| populations.checked_mul(block_count));
    if expected != u64::try_from(blocks.len()).ok() {
        return Err(configuration(
            config,
            "search.blocks",
            "la búsqueda no cubre exactamente cada bloque de cada subpoblación",
        ));
    }
    Ok(blocks)
}

/// Una fila exacta por bloque×subpoblación, ordenada por `population_id` y
/// `block_id` como exige la lista `sampling` del efectivo.
fn resolve_f5_sampling(
    config: &EffectiveF5Config,
    rosters: &[F5FounderRoster],
    blocks: &[UnitBlock],
) -> Result<Vec<F5SamplingResolution>, Diagnostic> {
    let overrides = index_overrides(config)?;
    let mut consumed: BTreeSet<(PopulationId, u64)> = BTreeSet::new();
    let mut rows = Vec::with_capacity(blocks.len());
    for block in blocks {
        let roster = rosters
            .iter()
            .find(|roster| roster.population_id == block.population_id)
            .ok_or_else(|| {
                configuration(
                    config,
                    "search.population_id",
                    "la búsqueda nombra una subpoblación que el escenario no declara",
                )
            })?;
        let subpopulation = declared_subpopulation(config, &block.population_id)?;
        let default_mode = match subpopulation.sampling.as_ref() {
            // `probabilistic` es el modo predeterminado del contrato.
            None => F5SamplingMode::Probabilistic,
            Some(mode) => F5SamplingMode::parse(mode.value()).ok_or_else(|| {
                configuration(
                    config,
                    "subpopulations.sampling",
                    "el modo declarado no es probabilistic ni controlled",
                )
            })?,
        };
        let key = (block.population_id.clone(), block.block_index);
        let substitution = overrides.get(&key);
        if substitution.is_some() {
            consumed.insert(key);
        }
        let (requested_mode, source, declared) = match substitution {
            Some(row) => {
                let mode = F5SamplingMode::parse(row.mode.value()).ok_or_else(|| {
                    configuration(
                        config,
                        "sampling_overrides.mode",
                        "el modo declarado no es probabilistic ni controlled",
                    )
                })?;
                let declared = row
                    .controlled_max_abs_error
                    .as_ref()
                    .map(|number| {
                        number.literal().parse::<f64>().map_err(|_| {
                            configuration(
                                config,
                                "sampling_overrides.controlled_max_abs_error",
                                "la tolerancia declarada no es un decimal recuperable",
                            )
                        })
                    })
                    .transpose()?;
                (mode, F5SamplingSource::BlockOverride, declared)
            }
            None => (default_mode, F5SamplingSource::SubpopulationDefault, None),
        };

        let default_tolerance = 0.5 / roster.copies as f64;
        let (declared, effective, copies_tolerance) = match requested_mode {
            F5SamplingMode::Probabilistic => {
                if declared.is_some() {
                    return Err(configuration(
                        config,
                        "sampling_overrides.controlled_max_abs_error",
                        "el modo probabilístico prohíbe declarar una tolerancia controlada",
                    ));
                }
                (None, None, None)
            }
            F5SamplingMode::Controlled => {
                // DTC-089: la combinación `F>0` con controlado es
                // incompatible y se rechaza antes de abrir ningún RNG.
                let inbreeding = effective_inbreeding_f(config, subpopulation)?.unwrap_or(0.0);
                if inbreeding != 0.0 {
                    return Err(configuration(
                        config,
                        "subpopulations.inbreeding_f,sampling",
                        format!(
                            "la unidad controlada {}×{} exige F efectivo igual a cero y su F es {inbreeding}",
                            block.population_id.as_str(),
                            block.block_id.as_str()
                        ),
                    ));
                }
                let effective = declared.unwrap_or(default_tolerance);
                if !effective.is_finite() || effective < 0.0 || effective > default_tolerance {
                    return Err(configuration(
                        config,
                        "sampling_overrides.controlled_max_abs_error",
                        format!(
                            "la tolerancia efectiva {effective} no pertenece a [0,{default_tolerance}]"
                        ),
                    ));
                }
                // En copias. El predeterminado es `0.5` **exacto**, no
                // `0.5/C` reescalado, que no vuelve a valer 0,5 para todo `C`.
                let copies_tolerance = match declared {
                    None => 0.5,
                    Some(value) => value * roster.copies as f64,
                };
                if !copies_tolerance.is_finite() || copies_tolerance < 0.0 {
                    return Err(configuration(
                        config,
                        "sampling_overrides.controlled_max_abs_error",
                        format!("la tolerancia {copies_tolerance} en copias no es finita"),
                    ));
                }
                (declared, Some(effective), Some(copies_tolerance))
            }
        };

        rows.push(F5SamplingResolution {
            population_id: block.population_id.clone(),
            block_id: block.block_id.clone(),
            block_index: block.block_index,
            requested_mode,
            // No hay degradación posible: una unidad inviable se rechaza.
            effective_mode: requested_mode,
            source,
            controlled_max_abs_error_declared: declared,
            controlled_max_abs_error_effective: effective,
            controlled_max_abs_error_copies: copies_tolerance,
            copies: roster.copies,
        });
    }
    // Cada sustitución debe corresponder a una unidad real de esta ejecución.
    // La población desconocida y el índice fuera de `0..B-1` ya los rechazan la
    // admisión y el preflight —que es quien resuelve `B`—; esta comprobación
    // impide que una sustitución quede sin aplicar sin que nadie lo note.
    if consumed.len() != overrides.len() {
        let orphan = overrides
            .keys()
            .find(|key| !consumed.contains(*key))
            .map(|(population_id, block_index)| format!("{}×{block_index}", population_id.as_str()))
            .unwrap_or_default();
        return Err(configuration(
            config,
            "sampling_overrides",
            format!("la sustitución {orphan} no corresponde a ninguna unidad de esta ejecución"),
        ));
    }
    rows.sort_by(|left, right| {
        left.population_id
            .cmp(&right.population_id)
            .then_with(|| left.block_id.cmp(&right.block_id))
    });
    Ok(rows)
}

fn declared_subpopulation<'a>(
    config: &'a EffectiveF5Config,
    population_id: &PopulationId,
) -> Result<&'a DeclaredF5Subpopulation, Diagnostic> {
    config
        .declared_scenario()
        .subpopulations
        .iter()
        .find(|subpopulation| subpopulation.population_id.value() == population_id.as_str())
        .ok_or_else(|| {
            configuration(
                config,
                "subpopulations.population_id",
                "la subpoblación resuelta no existe en el escenario declarado",
            )
        })
}

/// F declarado para una subpoblación: el suyo o el de `[f5]`.
///
/// Devuelve `None` cuando **ningún** documento lo declara, que no es lo mismo
/// que declarar `0.0`: REQ-EVO-48 rechaza *configurar* F fuera de `P=2` y esa
/// distinción es la que permite hacerlo.
fn effective_inbreeding_f(
    config: &EffectiveF5Config,
    subpopulation: &DeclaredF5Subpopulation,
) -> Result<Option<f64>, Diagnostic> {
    let declared = subpopulation
        .inbreeding_f
        .as_ref()
        .or(config.declared_scenario().inbreeding_f.as_ref());
    let Some(number) = declared else {
        return Ok(None);
    };
    number.literal().parse::<f64>().map(Some).map_err(|_| {
        configuration(
            config,
            "inbreeding_f",
            "el F declarado no es un decimal recuperable",
        )
    })
}

/// Sustituciones indexadas por `population_id`×`block_index`.
///
/// El duplicado ya lo rechaza la admisión; aquí se vuelve a comprobar porque
/// el índice no puede representarlo: sin esta guarda, dos filas para la misma
/// pareja se sobrescribirían en silencio.
fn index_overrides(
    config: &EffectiveF5Config,
) -> Result<BTreeMap<(PopulationId, u64), &DeclaredF5SamplingOverride>, Diagnostic> {
    let mut indexed = BTreeMap::new();
    for row in &config.declared_scenario().sampling_overrides {
        let population_id = PopulationId::parse(
            row.population_id.value().clone(),
            "sampling_overrides.population_id",
        )
        .map_err(|_| {
            configuration(
                config,
                "sampling_overrides.population_id",
                "no es una identidad tipada",
            )
        })?;
        if indexed
            .insert((population_id, *row.block_index.value()), row)
            .is_some()
        {
            return Err(configuration(
                config,
                "sampling_overrides",
                "dos sustituciones nombran la misma pareja población×bloque",
            ));
        }
    }
    Ok(indexed)
}

// ---------------------------------------------------------------------------
// Identidades de copia
// ---------------------------------------------------------------------------

/// Identidad de copia dentro de un bloque: `(individual_id, copy_index)`.
///
/// El bloque fija la secuencia, de modo que la tupla completa de DTC-089
/// —referencia, individuo, secuencia y copia— queda determinada por la unidad
/// y no hace falta repetirla aquí.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CopySeat {
    individual_id: IndividualId,
    copy_index: u32,
}

/// Las `C` identidades de copia ordenadas por `individual_id` y `copy_index`,
/// que es el orden inicial que exige la permutación controlada.
fn copy_seats(roster: &F5FounderRoster, ploidy: u8) -> Vec<CopySeat> {
    let mut seats = Vec::with_capacity(roster.individuals.len() * usize::from(ploidy));
    for individual in &roster.individuals {
        for copy_index in 0..u32::from(ploidy) {
            seats.push(CopySeat {
                individual_id: individual.individual_id.clone(),
                copy_index,
            });
        }
    }
    seats.sort();
    seats
}

// ---------------------------------------------------------------------------
// Modo probabilístico
// ---------------------------------------------------------------------------

/// F aplicable a las copias de una unidad probabilística.
///
/// Las tres variantes son exactamente los tres consumos que DTC-089 fija para
/// `founder_f`: ninguno con `F=0`, ninguno con `F=1` y un U53 con `0<F<1`.
#[derive(Clone, Copy, Debug, PartialEq)]
enum FounderInbreeding {
    /// `P≠2`, F sin declarar o `F=0`: cada copia se muestrea por su cuenta y
    /// `founder_f` no se abre.
    Independent,
    /// `P=2` con `0<F<1`: la copia 1 consulta `founder_f` antes de su propio
    /// `founder_sample`.
    Bernoulli(f64),
    /// `P=2` con `F=1`: la copia 1 reutiliza siempre y `founder_f` no se abre,
    /// conforme a la regla Bernoulli del catálogo.
    AlwaysReuse,
}

/// Resuelve el F aplicable a una subpoblación a partir de la ploidía efectiva
/// y del F **declarado**.
///
/// `None` significa que ningún documento declaró F. REQ-EVO-48 rechaza
/// configurar F con `P=1` o `P>2` —no solo aplicarlo—, de modo que aquí basta
/// con que exista una declaración para rechazarla fuera de `P=2`, y el
/// diagnóstico cita ámbito y ploidía.
fn resolve_founder_inbreeding(
    ploidy: u8,
    declared: Option<f64>,
) -> Result<FounderInbreeding, String> {
    let Some(inbreeding) = declared else {
        return Ok(FounderInbreeding::Independent);
    };
    if ploidy != 2 {
        return Err(format!(
            "declara inbreeding_f={inbreeding} y F solo se admite con ploidía 2; la ploidía efectiva es {ploidy}"
        ));
    }
    if !inbreeding.is_finite() || !(0.0..=1.0).contains(&inbreeding) {
        return Err(format!(
            "declara inbreeding_f={inbreeding}, que no pertenece a [0,1]"
        ));
    }
    if inbreeding == 0.0 {
        return Ok(FounderInbreeding::Independent);
    }
    if inbreeding == 1.0 {
        return Ok(FounderInbreeding::AlwaysReuse);
    }
    Ok(FounderInbreeding::Bernoulli(inbreeding))
}

fn sample_probabilistic_unit(
    scenario_id: &ScenarioId,
    root_seed: RootSeed,
    row: &F5SamplingResolution,
    seats: &[CopySeat],
    prefix: &HaplotypePrefix,
    inbreeding: FounderInbreeding,
) -> Result<(Vec<F5CopyAssignment>, u64), Diagnostic> {
    let mut assignments = Vec::with_capacity(seats.len());
    let mut operations = 0_u64;
    // Las copias de un individuo se recorren juntas y en orden de
    // `copy_index`, de modo que la fila de la copia 0 sigue disponible cuando
    // le toca a la copia 1 y la reutilización no necesita recorrer nada.
    let mut previous: Option<(&IndividualId, u64)> = None;
    for seat in seats {
        operations = operations.saturating_add(1);
        // Solo la copia 1 **del mismo individuo** puede reutilizar una fila.
        let copy_zero_rank = match previous {
            Some((individual, rank))
                if seat.copy_index == 1 && *individual == seat.individual_id =>
            {
                Some(rank)
            }
            _ => None,
        };
        let (haplotype_rank, visits) = draw_copy_rank(
            scenario_id,
            inbreeding,
            copy_zero_rank,
            prefix,
            || {
                Ok(derive_f5_subseed(
                    root_seed,
                    F5RandomDecision::founder_f(
                        scenario_id,
                        &row.population_id,
                        &row.block_id,
                        &seat.individual_id,
                    ),
                )?
                .local_rng())
            },
            || {
                Ok(derive_f5_subseed(
                    root_seed,
                    F5RandomDecision::founder_sample(
                        scenario_id,
                        &row.population_id,
                        &row.block_id,
                        &seat.individual_id,
                        seat.copy_index,
                    ),
                )?
                .local_rng())
            },
        )?;
        operations = operations.saturating_add(visits);
        previous = Some((&seat.individual_id, haplotype_rank));
        assignments.push(F5CopyAssignment {
            individual_id: seat.individual_id.clone(),
            copy_index: seat.copy_index,
            haplotype_rank,
        });
    }
    Ok((assignments, operations))
}

/// Fila de una copia y **orden exacto en que abre sus ámbitos**.
///
/// `copy_zero_rank` es la fila de la copia 0 del mismo individuo, y solo la
/// trae la copia 1 de un diploide. El orden que fija DTC-089 —`founder_f`
/// antes que el `founder_sample` de la copia 1, y ese segundo flujo sin abrir
/// cuando se reutiliza— no se observa en ninguna salida: por eso los dos
/// ámbitos entran como cierres y las unitarias pueden comprobar cuál se abre,
/// en qué orden y cuál no llega a invocarse.
///
/// Devuelve el rango y las visitas que costó, `founder_f` incluido.
fn draw_copy_rank<R: RandomWordSource>(
    scenario_id: &ScenarioId,
    inbreeding: FounderInbreeding,
    copy_zero_rank: Option<u64>,
    prefix: &HaplotypePrefix,
    open_f: impl FnOnce() -> Result<R, Diagnostic>,
    open_sample: impl FnOnce() -> Result<R, Diagnostic>,
) -> Result<(u64, u64), Diagnostic> {
    let mut operations = 0_u64;
    if let Some(rank) = copy_zero_rank {
        let reuse = match inbreeding {
            FounderInbreeding::Independent => false,
            FounderInbreeding::AlwaysReuse => true,
            FounderInbreeding::Bernoulli(probability) => {
                let mut rng = open_f()?;
                operations = operations.saturating_add(1);
                f5_bernoulli_with_source(&mut rng, probability).map_err(|_| {
                    internal(
                        scenario_id,
                        format!("el F efectivo {probability} no es una probabilidad"),
                    )
                })?
            }
        };
        if reuse {
            // El segundo flujo no se abre: la copia 1 recibe **exactamente**
            // la fila de la copia 0.
            return Ok((rank, operations));
        }
    }
    let mut rng = open_sample()?;
    let (rank, visits) = prefix.choose(&mut rng);
    Ok((rank, operations.saturating_add(visits)))
}

/// Prefijos compensados de `q`, calculados **una vez por unidad**.
///
/// DTC-089 elige fila recorriendo `q` en orden de rango y quedándose con el
/// primer prefijo acumulado estrictamente mayor que el U53. Ese prefijo no
/// depende de la copia: es el mismo vector para las `C` copias de la unidad,
/// así que se compone una vez —con la misma suma de Neumaier y en el mismo
/// orden— y cada copia lo consulta. Los prefijos son no decrecientes, de modo
/// que la primera posición que supera el umbral se localiza por partición
/// binaria y el resultado es exactamente el del recorrido secuencial.
struct HaplotypePrefix {
    /// Prefijos compensados, uno por rango.
    prefixes: Vec<f64>,
    /// Último rango con peso positivo, que absorbe el residuo de redondeo.
    last_positive: u64,
    /// `true` cuando solo hay un haplotipo: no hay elección ni RNG.
    single: bool,
}

impl HaplotypePrefix {
    /// # Errors
    ///
    /// Devuelve el detalle si algún peso no es finito y no negativo, si
    /// ninguno es positivo o si la suma compensada deja de ser finita.
    fn new(weights: &[f64]) -> Result<Self, &'static str> {
        let mut prefixes = Vec::with_capacity(weights.len());
        let mut compensated = NeumaierSum::default();
        let mut last_positive: Option<u64> = None;
        for (rank, weight) in weights.iter().enumerate() {
            if !weight.is_finite() || *weight < 0.0 {
                return Err("un peso haplotípico no es finito y no negativo");
            }
            if *weight > 0.0 {
                last_positive = Some(rank as u64);
            }
            prefixes.push(compensated.add(*weight)?);
        }
        let last_positive = last_positive.ok_or("ningún haplotipo tiene peso positivo")?;
        Ok(Self {
            single: weights.len() == 1,
            prefixes,
            last_positive,
        })
    }

    /// Rango elegido y visitas que costó localizarlo.
    ///
    /// Una elección con una sola alternativa no consume RNG. En cualquier otro
    /// caso consume exactamente un U53, y el residuo de redondeo —cuando el
    /// umbral no llega a ningún prefijo— lo absorbe el último rango con peso
    /// positivo, no la última fila: un `q_h` nulo no puede recibir copias.
    fn choose(&self, rng: &mut impl RandomWordSource) -> (u64, u64) {
        if self.single {
            return (0, 0);
        }
        let threshold = f5_uniform_53(rng);
        let mut visits = 0_u64;
        let position = self.prefixes.partition_point(|prefix| {
            visits = visits.saturating_add(1);
            *prefix <= threshold
        });
        let rank = u64::try_from(position).unwrap_or(u64::MAX);
        (rank.min(self.last_positive), visits)
    }
}

// ---------------------------------------------------------------------------
// Modo controlado
// ---------------------------------------------------------------------------

/// Estado de un rango dentro del reparto por mayores restos.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RankState {
    /// Grupo excluido por completo: solo el piso es aplicable.
    Floor,
    /// Grupo incluido por completo: solo el piso más uno es aplicable.
    FloorPlusOne,
    /// Grupo fronterizo: los dos estados son aplicables.
    Boundary,
}

/// Reparto por mayores restos ya comprobado, antes del desempate.
struct QuotaPlan {
    floors: Vec<u64>,
    states: Vec<RankState>,
    /// Rangos empatados en la frontera, en orden ascendente.
    boundary_ranks: Vec<u64>,
    boundary_remainder: f64,
    select_count: u64,
    /// Visitas y comparaciones realmente ejecutadas, contadas dentro de cada
    /// bucle. No es una identidad de coste: crece con el trabajo hecho.
    operations: u64,
}

/// Causa exacta por la que una cuota controlada no puede satisfacerse.
enum QuotaError {
    /// Aritmética indefinida o desbordamiento: invariante interna.
    Internal(String),
    /// Cota incumplida o suma `C` inalcanzable: `GPS-DOMAIN-103`.
    Infeasible(String),
}

/// Reparto por mayores restos con comprobación previa de factibilidad.
///
/// La comprobación es `O(H)` y demuestra la cota para **todas** las
/// combinaciones óptimas sin enumerar ninguna: cada rango fijo tiene un único
/// estado posible y los del grupo fronterizo comparten resto, de modo que
/// exigir la cota en sus dos estados —piso y piso más uno— la exige en
/// cualquier prefijo que el desempate elija.
///
/// La cota llega **en copias**: `|n_h - C*q_h| <= t*C`. Es la desigualdad del
/// contrato multiplicada por `C`, y se evalúa así porque `C*q_h` y su
/// distancia a un entero contiguo son exactos en binary64 con pesos diádicos,
/// mientras que `n_h/C` introduce un redondeo que rechazaba cuotas admisibles
/// —entre ellas la del resto `0,5`, que es justo la que genera el empate
/// fronterizo que el contrato regula—.
///
/// El umbral de los mayores restos se obtiene con una selección parcial en
/// `O(H)`, no ordenando: solo hace falta saber qué resto ocupa la posición
/// `residual-1` del orden completo y quiénes lo igualan.
fn plan_controlled_quota(
    weights: &[f64],
    copies: u64,
    tolerance_copies: f64,
) -> Result<QuotaPlan, QuotaError> {
    if copies == 0 {
        return Err(QuotaError::Internal(
            "una subpoblación con C=0 no puede repartir cuotas".to_owned(),
        ));
    }
    if !tolerance_copies.is_finite() || tolerance_copies < 0.0 {
        return Err(QuotaError::Internal(format!(
            "la tolerancia efectiva {tolerance_copies} en copias no es finita y no negativa"
        )));
    }
    let copies_f64 = copies as f64;
    let mut operations = 0_u64;
    let mut floors = Vec::with_capacity(weights.len());
    let mut remainders = Vec::with_capacity(weights.len());
    let mut base = 0_u64;
    for weight in weights {
        operations = operations.saturating_add(1);
        let ideal = copies_f64 * *weight;
        if !ideal.is_finite() || ideal < 0.0 {
            return Err(QuotaError::Internal(format!(
                "la cuota ideal {ideal} de un peso {weight} no es finita y no negativa"
            )));
        }
        let floor = ideal.floor();
        if floor > copies_f64 {
            return Err(QuotaError::Infeasible(format!(
                "el piso {floor} de un peso {weight} supera C={copies}"
            )));
        }
        let floor_u64 = floor as u64;
        base = base
            .checked_add(floor_u64)
            .ok_or_else(|| QuotaError::Internal("la suma de pisos desborda u64".to_owned()))?;
        floors.push(floor_u64);
        remainders.push(ideal - floor);
    }
    let residual = copies.checked_sub(base).ok_or_else(|| {
        QuotaError::Infeasible(format!(
            "los pisos suman {base} y superan C={copies}: ninguna cuota entera suma C"
        ))
    })?;
    if residual > weights.len() as u64 {
        return Err(QuotaError::Infeasible(format!(
            "el residual {residual} excede los {} rangos disponibles",
            weights.len()
        )));
    }
    let residual_position = usize::try_from(residual)
        .map_err(|_| QuotaError::Internal("el residual no cabe en usize".to_owned()))?;

    let mut states = vec![RankState::Floor; weights.len()];
    let mut boundary_ranks = Vec::new();
    let mut boundary_remainder = 0.0_f64;
    let mut select_count = 0_u64;
    if residual_position > 0 {
        // Selección parcial en `O(H)`: basta el resto que ocuparía la posición
        // `residual-1` del orden decreciente. Ordenar los `H` restos costaría
        // `H log H` y no aporta nada, porque dentro del grupo empatado el
        // orden lo decide el desempate y fuera de él solo importa el lado.
        let mut order: Vec<usize> = (0..weights.len()).collect();
        operations = operations.saturating_add(order.len() as u64);
        let pivot = {
            let (_, pivot, _) =
                order.select_nth_unstable_by(residual_position - 1, |left, right| {
                    operations = operations.saturating_add(1);
                    remainders[*right]
                        .total_cmp(&remainders[*left])
                        .then_with(|| left.cmp(right))
                });
            *pivot
        };
        let threshold = remainders[pivot];
        let mut greater = 0_u64;
        let mut group = Vec::new();
        for (rank, remainder) in remainders.iter().enumerate() {
            operations = operations.saturating_add(1);
            // Igualdad significa el mismo valor finito, sin tolerancia alguna.
            if *remainder > threshold {
                states[rank] = RankState::FloorPlusOne;
                greater = greater.saturating_add(1);
            } else if *remainder == threshold {
                group.push(rank as u64);
            }
        }
        let selected = residual.checked_sub(greater).ok_or_else(|| {
            QuotaError::Internal("los restos mayores que el umbral superan el residual".to_owned())
        })?;
        if selected > group.len() as u64 {
            return Err(QuotaError::Internal(
                "el prefijo empatado excede el tamaño de su grupo".to_owned(),
            ));
        }
        if selected == group.len() as u64 {
            // El corte cae justo en el borde del grupo: entra entero y no hay
            // frontera que sortear, de modo que el prefijo empatado es cero.
            for rank in &group {
                states[*rank as usize] = RankState::FloorPlusOne;
            }
        } else {
            select_count = selected;
            for rank in &group {
                states[*rank as usize] = RankState::Boundary;
            }
            // `group` se construye recorriendo los rangos en orden creciente,
            // que es el orden que el desempate exige antes de permutar.
            boundary_ranks = group;
            boundary_remainder = threshold;
        }
    }

    // Factibilidad previa: ningún estado aplicable puede salirse de la cota.
    for (rank, weight) in weights.iter().enumerate() {
        operations = operations.saturating_add(1);
        let floor = floors[rank];
        let ideal = copies_f64 * *weight;
        let applicable: [Option<u64>; 2] = match states[rank] {
            RankState::Floor => [Some(floor), None],
            RankState::FloorPlusOne => [Some(floor + 1), None],
            RankState::Boundary => [Some(floor), Some(floor + 1)],
        };
        for count in applicable.into_iter().flatten() {
            let error = (count as f64 - ideal).abs();
            if !error.is_finite() {
                return Err(QuotaError::Internal(format!(
                    "el error de cuota del rango {rank} no es finito"
                )));
            }
            if error > tolerance_copies {
                return Err(QuotaError::Infeasible(format!(
                    "rango {rank}: |{count}-{ideal}|={error} copias supera la tolerancia efectiva {tolerance_copies} copias (q={weight}, C={copies})"
                )));
            }
        }
    }

    Ok(QuotaPlan {
        floors,
        states,
        boundary_ranks,
        boundary_remainder,
        select_count,
        operations,
    })
}

/// Error del desempate: el del ámbito RNG o una invariante interna propia.
enum TieError {
    Scope(Diagnostic),
    Internal(&'static str),
}

impl TieError {
    fn into_diagnostic(self, scenario_id: &ScenarioId) -> Diagnostic {
        match self {
            Self::Scope(diagnostic) => diagnostic,
            Self::Internal(detail) => internal(scenario_id, detail),
        }
    }
}

impl From<Diagnostic> for TieError {
    fn from(diagnostic: Diagnostic) -> Self {
        Self::Scope(diagnostic)
    }
}

/// Prefijo sorteado del grupo fronterizo, o `None` si ningún grupo cruza la
/// frontera.
///
/// `open_scope` **solo se invoca cuando existe esa frontera**. Es la forma de
/// hacer comprobable la regla de DTC-089 —«grupos incluidos o excluidos por
/// completo no consumen RNG»—: sin frontera no se deriva la subsemilla de
/// `controlled_quota_tie`, no solo no se sortea con ella.
fn select_tie_prefix(
    plan: &QuotaPlan,
    open_scope: impl FnOnce() -> Result<LocalRng, Diagnostic>,
) -> Result<(Option<F5QuotaTie>, u64), TieError> {
    if plan.boundary_ranks.is_empty() {
        return Ok((None, 0));
    }
    let mut candidates = plan.boundary_ranks.clone();
    let mut rng = open_scope()?;
    f5_fisher_yates_with_source(&mut rng, &mut candidates)
        .map_err(|_| TieError::Internal("la permutación del grupo empatado no pudo completarse"))?;
    let selected = usize::try_from(plan.select_count)
        .map_err(|_| TieError::Internal("el prefijo empatado no cabe en usize"))?;
    let selected_ranks = candidates
        .get(..selected)
        .ok_or(TieError::Internal(
            "el prefijo empatado excede el tamaño de su grupo",
        ))?
        .to_vec();
    let operations = candidates.len().saturating_sub(1) as u64;
    Ok((
        Some(F5QuotaTie {
            remainder: plan.boundary_remainder,
            candidate_ranks: plan.boundary_ranks.clone(),
            select_count: plan.select_count,
            selected_ranks,
        }),
        operations,
    ))
}

fn sample_controlled_unit(
    scenario_id: &ScenarioId,
    root_seed: RootSeed,
    row: &F5SamplingResolution,
    seats: &[CopySeat],
    weights: &[f64],
) -> Result<(Vec<F5CopyAssignment>, F5ControlledQuota, u64), Diagnostic> {
    let tolerance = row.controlled_max_abs_error_effective.ok_or_else(|| {
        internal(
            scenario_id,
            "una unidad controlada llegó sin tolerancia efectiva",
        )
    })?;
    let tolerance_copies = row.controlled_max_abs_error_copies.ok_or_else(|| {
        internal(
            scenario_id,
            "una unidad controlada llegó sin tolerancia efectiva en copias",
        )
    })?;
    let plan =
        plan_controlled_quota(weights, row.copies, tolerance_copies).map_err(
            |error| match error {
                QuotaError::Internal(detail) => internal(scenario_id, detail),
                QuotaError::Infeasible(detail) => Diagnostic::f5_controlled_quota_infeasible(
                    scenario_id.as_str(),
                    row.population_id.as_str(),
                    format!(
                        "block_id={};tolerance={tolerance};{detail}",
                        row.block_id.as_str()
                    ),
                ),
            },
        )?;

    let mut operations = plan.operations;
    let mut counts = plan.floors.clone();
    for (rank, state) in plan.states.iter().enumerate() {
        operations = operations.saturating_add(1);
        if *state == RankState::FloorPlusOne {
            counts[rank] += 1;
        }
    }
    // Un grupo incluido o excluido por completo no consume RNG: solo el que
    // cruza la frontera abre el ámbito `controlled_quota_tie`.
    let (tie, tie_operations) = select_tie_prefix(&plan, || {
        Ok(derive_f5_subseed(
            root_seed,
            F5RandomDecision::controlled_quota_tie(scenario_id, &row.population_id, &row.block_id),
        )?
        .local_rng())
    })
    .map_err(|error| error.into_diagnostic(scenario_id))?;
    operations = operations.saturating_add(tie_operations);
    if let Some(tie) = &tie {
        for rank in tie.selected_ranks() {
            let slot = usize::try_from(*rank)
                .map_err(|_| internal(scenario_id, "un rango empatado no cabe en usize"))?;
            let count = counts
                .get_mut(slot)
                .ok_or_else(|| internal(scenario_id, "un rango empatado no existe en la cuota"))?;
            *count += 1;
        }
    }

    let total = counts
        .iter()
        .try_fold(0_u64, |total, count| total.checked_add(*count))
        .ok_or_else(|| internal(scenario_id, "la suma de cuotas desborda u64"))?;
    if total != row.copies {
        return Err(internal(
            scenario_id,
            format!("las cuotas suman {total} y C vale {}", row.copies),
        ));
    }

    // El multiconjunto se ordena por `haplotype_rank` y la permutación de las
    // identidades de copia decide qué copia recibe cada posición.
    let mut multiset = Vec::with_capacity(seats.len());
    for (rank, count) in counts.iter().enumerate() {
        for _ in 0..*count {
            operations = operations.saturating_add(1);
            multiset.push(rank as u64);
        }
    }
    // Se permutan los **índices** de las identidades de copia, no las
    // identidades. Fisher–Yates intercambia posiciones, de modo que la
    // secuencia RNG y la permutación resultante son idénticas a permutar los
    // asientos, pero el resultado se escribe ya en el orden canónico
    // `individual_id`+`copy_index` y no hay que reordenar `C` asignaciones
    // después: el coste baja de `C log C` a `C`.
    let mut permutation: Vec<u32> = (0..u32::try_from(seats.len())
        .map_err(|_| internal(scenario_id, "C no cabe en la permutación de copias"))?)
        .collect();
    let mut rng = derive_f5_subseed(
        root_seed,
        F5RandomDecision::founder_controlled_permutation(
            scenario_id,
            &row.population_id,
            &row.block_id,
        ),
    )?
    .local_rng();
    f5_fisher_yates_with_source(&mut rng, &mut permutation).map_err(|_| {
        internal(
            scenario_id,
            "la permutación de identidades de copia no pudo completarse",
        )
    })?;
    operations = operations.saturating_add(permutation.len().saturating_sub(1) as u64);

    let mut ranks = vec![0_u64; seats.len()];
    for (position, haplotype_rank) in multiset.into_iter().enumerate() {
        operations = operations.saturating_add(1);
        let seat = *permutation
            .get(position)
            .ok_or_else(|| internal(scenario_id, "la permutación no cubre el multiconjunto"))?;
        let slot = ranks
            .get_mut(seat as usize)
            .ok_or_else(|| internal(scenario_id, "la permutación nombra una copia inexistente"))?;
        *slot = haplotype_rank;
    }
    let assignments = seats
        .iter()
        .zip(ranks)
        .map(|(seat, haplotype_rank)| F5CopyAssignment {
            individual_id: seat.individual_id.clone(),
            copy_index: seat.copy_index,
            haplotype_rank,
        })
        .collect::<Vec<_>>();

    Ok((
        assignments,
        F5ControlledQuota {
            copies: row.copies,
            counts,
            max_abs_error_effective: tolerance,
            tie,
        },
        operations,
    ))
}

// ---------------------------------------------------------------------------
// F realizado
// ---------------------------------------------------------------------------

/// Fila `f_realized` de una unidad, en los dos modos y en toda ploidía, junto
/// a las visitas que costó contarla.
fn realized_f_row(
    scenario_id: &ScenarioId,
    row: &F5SamplingResolution,
    ploidy: u8,
    target_f: Option<f64>,
    assignments: &[F5CopyAssignment],
    haplotype_count: usize,
) -> Result<(F5RealizedF, u64), Diagnostic> {
    // Fuera de `P=2` no hay heterocigoto que contar: la fila se publica igual,
    // con objetivo y valor nulos y estado `not_applicable`, nunca con un cero
    // ni con una estimación adaptada (REQ-EVO-52).
    let (value, operations) = if ploidy == 2 {
        realized_diploid_f(assignments, haplotype_count)
            .map_err(|detail| internal(scenario_id, detail))?
    } else {
        (None, 0)
    };
    Ok((
        F5RealizedF {
            population_id: row.population_id.clone(),
            block_id: row.block_id.clone(),
            target_f,
            status: if value.is_some() {
                F5RealizedFStatus::Evaluated
            } else {
                F5RealizedFStatus::NotApplicable
            },
            value,
        },
        operations,
    ))
}

/// `F_realized = 1 - H_observed/H_expected` de una unidad **diploide**, o
/// `None` cuando el denominador no lo permite, junto a las visitas que costó.
///
/// `H_observed = heterozygous_individuals/individual_count`,
/// `q_realized_h = n_h/C` y `H_expected = 1 - Σ_h q_realized_h²`. La suma
/// recorre los rangos en orden ascendente y se evalúa en el orden escrito, que
/// es la regla binary64 de DTC-089; no interviene aquí la suma compensada de
/// DTC-049, que gobierna las normalizaciones del perfil y no esta fórmula.
///
/// Devuelve `None` —`not_applicable`, nunca cero— cuando `H_expected` no es
/// positivo. El caso del contrato es `H_expected = 0`: con una sola categoría
/// realizada `q` vale 1 exacto y la resta da cero **exacto**, sin depender del
/// redondeo. Las otras dos formas de no ser positivo —negativo o `NaN`— son
/// defensa **inalcanzable**: con más de una categoría el mínimo exacto de
/// `C² - Σ n_h²` es `2C - 2`, de modo que `H_expected >= 2/C - 2/C²` domina
/// cualquier error de redondeo, y una búsqueda dirigida hasta `C = 2^51` no
/// encontró ninguno. Se conservan porque un cociente sin denominador no es un
/// resultado, no porque el contrato las prevea.
///
/// Un valor negativo es un exceso de heterocigotos y se conserva tal cual
/// (REQ-EVO-52).
///
/// # Errors
///
/// Devuelve el detalle si las asignaciones no vienen emparejadas por individuo
/// y copia, o si nombran un rango que la matriz no contiene: en un diploide
/// ambas cosas son invariantes internas ya garantizadas por el muestreo.
fn realized_diploid_f(
    assignments: &[F5CopyAssignment],
    haplotype_count: usize,
) -> Result<(Option<f64>, u64), String> {
    let mut counts = vec![0_u64; haplotype_count];
    let mut individuals = 0_u64;
    let mut heterozygous = 0_u64;
    let mut index = 0_usize;
    // El recuento es `O(C + H)` y se contabiliza **dentro** de sus dos bucles,
    // igual que el resto del coste de la unidad.
    let mut operations = 0_u64;
    while index < assignments.len() {
        operations = operations.saturating_add(1);
        let first = &assignments[index];
        let second = assignments.get(index + 1).ok_or_else(|| {
            "una unidad diploide deja una copia sin pareja en su individuo".to_owned()
        })?;
        if first.copy_index != 0
            || second.copy_index != 1
            || first.individual_id != second.individual_id
        {
            return Err(format!(
                "el individuo {} no aporta exactamente las copias 0 y 1 consecutivas",
                first.individual_id.as_str()
            ));
        }
        for assignment in [first, second] {
            let slot = usize::try_from(assignment.haplotype_rank)
                .ok()
                .and_then(|rank| counts.get_mut(rank))
                .ok_or_else(|| {
                    "una asignación nombra un rango que la matriz no contiene".to_owned()
                })?;
            *slot += 1;
        }
        individuals += 1;
        if first.haplotype_rank != second.haplotype_rank {
            heterozygous += 1;
        }
        index += 2;
    }
    // `C=0` no tiene frecuencia realizada. El censo admitido exige al menos un
    // individuo por subpoblación, de modo que la vía productiva no llega aquí;
    // la rama existe para que la función sea total y la fijan las unitarias.
    if individuals == 0 {
        return Ok((None, operations));
    }
    let copies = assignments.len() as f64;
    let h_observed = heterozygous as f64 / individuals as f64;
    let mut homozygosity = 0.0_f64;
    for count in &counts {
        operations = operations.saturating_add(1);
        let frequency = *count as f64 / copies;
        homozygosity += frequency * frequency;
    }
    let h_expected = 1.0 - homozygosity;
    // `== 0.0` es el caso del contrato; `< 0.0` y `NaN` son inalcanzables y se
    // conservan solo para no dividir nunca sin denominador.
    if h_expected.is_nan() || h_expected <= 0.0 {
        return Ok((None, operations));
    }
    let value = 1.0 - h_observed / h_expected;
    if !value.is_finite() {
        return Err(format!(
            "F realizado no finito con H_observed={h_observed} y H_expected={h_expected}"
        ));
    }
    Ok((Some(value), operations))
}

// ---------------------------------------------------------------------------
// F5.3.3 — identidades fundadoras, pertenencia y enlaces `FCP_`/`FHA_`→`HAP_`
// ---------------------------------------------------------------------------

/// Copia fundadora persistente: su `FCP_` y la subpoblación a la que
/// pertenece.
///
/// La tupla publicada es exactamente la que DTC-089 da a la proyección de
/// verdad —`(reference_id, individual_id, sequence_id, copy_index)`— más la
/// pertenencia, de modo que un lector rederiva el `FCP_` con DTC-039 sin
/// preguntar al productor. La referencia es única por ejecución F5 y vive en la
/// configuración efectiva, no en cada fila.
///
/// La copia se identifica por su índice, no por una fase: en F5 no hay meiosis
/// y el identificador de fase es aportado por la usuaria y sin constructor F5
/// (F5-R04).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5FounderCopy {
    population_id: PopulationId,
    individual_id: IndividualId,
    sequence_id: SequenceId,
    copy_index: u32,
    founder_copy_id: FounderCopyId,
}

impl F5FounderCopy {
    /// La pertenencia de la copia: la única subpoblación que la contiene.
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }
    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }
    #[must_use]
    pub const fn copy_index(&self) -> u32 {
        self.copy_index
    }
    #[must_use]
    pub const fn founder_copy_id(&self) -> &FounderCopyId {
        &self.founder_copy_id
    }
}

/// Fila haplotípica materializada de un bloque×subpoblación: su rango, sus
/// alelos y su `HAP_`.
///
/// `alleles` es la fila de la matriz ganadora tal cual: exactamente `S`
/// caracteres `0`/`1`. La polaridad de DTC-089 es una propiedad publicada del
/// sitio en `f5-frequencies.tsv` y no reescribe estos bits.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5HaplotypeIdentity {
    population_id: PopulationId,
    block_id: F5BlockId,
    block_index: u64,
    haplotype_rank: u64,
    alleles: String,
    haplotype_id: F5HaplotypeId,
}

impl F5HaplotypeIdentity {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    #[must_use]
    pub const fn block_id(&self) -> &F5BlockId {
        &self.block_id
    }
    /// El índice que DTC-039 usa como componente del `HAP_`; el orden canónico
    /// del layout lo recupera desde `f5-block-layout.tsv`.
    #[must_use]
    pub const fn block_index(&self) -> u64 {
        self.block_index
    }
    #[must_use]
    pub const fn haplotype_rank(&self) -> u64 {
        self.haplotype_rank
    }
    #[must_use]
    pub fn alleles(&self) -> &str {
        &self.alleles
    }
    #[must_use]
    pub const fn haplotype_id(&self) -> &F5HaplotypeId {
        &self.haplotype_id
    }
}

/// Asignación por bloque de una copia a una fila haplotípica.
///
/// DTC-089 es explícita: `FHA_` **no es una copia**, es la asignación por
/// bloque a un `HAP_` y nunca sustituye la identidad persistente de la copia.
/// Por eso esta fila conserva las dos identidades y no deriva una de la otra.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5BlockAssignment {
    population_id: PopulationId,
    individual_id: IndividualId,
    sequence_id: SequenceId,
    copy_index: u32,
    block_id: F5BlockId,
    block_index: u64,
    founder_copy_id: FounderCopyId,
    block_assignment_id: FounderHaplotypeAssignmentId,
    haplotype_id: F5HaplotypeId,
    haplotype_rank: u64,
}

impl F5BlockAssignment {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }
    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }
    #[must_use]
    pub const fn copy_index(&self) -> u32 {
        self.copy_index
    }
    #[must_use]
    pub const fn block_id(&self) -> &F5BlockId {
        &self.block_id
    }
    #[must_use]
    pub const fn block_index(&self) -> u64 {
        self.block_index
    }
    #[must_use]
    pub const fn founder_copy_id(&self) -> &FounderCopyId {
        &self.founder_copy_id
    }
    #[must_use]
    pub const fn block_assignment_id(&self) -> &FounderHaplotypeAssignmentId {
        &self.block_assignment_id
    }
    #[must_use]
    pub const fn haplotype_id(&self) -> &F5HaplotypeId {
        &self.haplotype_id
    }
    #[must_use]
    pub const fn haplotype_rank(&self) -> u64 {
        self.haplotype_rank
    }
}

/// Repositorio de identidades fundadoras materializado **antes de toda
/// salida**.
///
/// Contiene las tres poblaciones de identidades que F5.4 escribe y ninguna
/// más: las copias `FCP_` con su pertenencia, las filas `HAP_` de cada
/// bloque×subpoblación y las asignaciones `FHA_` que enlazan unas con otras.
/// No abre ningún RNG ni decide ninguna asignación: consume la que F5.3.1 y
/// F5.3.2 ya resolvieron.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5FounderIdentities {
    copies: Vec<F5FounderCopy>,
    haplotypes: Vec<F5HaplotypeIdentity>,
    assignments: Vec<F5BlockAssignment>,
    membership: BTreeMap<FounderCopyId, usize>,
}

impl F5FounderIdentities {
    /// Copias ordenadas por `population_id`, `individual_id`, `sequence_id` y
    /// `copy_index`.
    #[must_use]
    pub fn copies(&self) -> &[F5FounderCopy] {
        &self.copies
    }
    /// Filas haplotípicas ordenadas por `population_id`, `block_id` y
    /// `haplotype_rank`, que es la clave de `f5-haplotypes.tsv`.
    #[must_use]
    pub fn haplotypes(&self) -> &[F5HaplotypeIdentity] {
        &self.haplotypes
    }
    /// Asignaciones ordenadas por `population_id`, `individual_id`,
    /// `sequence_id`, `copy_index` y `block_id`, que es la clave de
    /// `f5-founders.tsv`.
    #[must_use]
    pub fn assignments(&self) -> &[F5BlockAssignment] {
        &self.assignments
    }

    /// Resuelve la pertenencia de una copia desde su identidad publicada.
    ///
    /// Es la operación que la proyección de verdad necesita para escribir
    /// `founder_membership` sin inferir la subpoblación del identificador, y
    /// es total e inyectiva por construcción: cada `FCP_` se deriva una sola
    /// vez y el registro de identidades rechaza cualquier repetición con
    /// `GPS-DOMAIN-002`.
    #[must_use]
    pub fn copy(&self, founder_copy_id: &FounderCopyId) -> Option<&F5FounderCopy> {
        self.membership
            .get(founder_copy_id)
            .and_then(|position| self.copies.get(*position))
    }
}

/// Lo único que la identidad necesita saber de un bloque colocado.
#[derive(Clone, Debug)]
struct BlockGeometry {
    block_id: F5BlockId,
    sequence_id: SequenceId,
    site_count: usize,
}

/// Materializa las identidades de fundador y su pertenencia única antes de
/// toda salida.
///
/// Deriva con DTC-039/DTC-037 —nunca por concatenación textual— el `FCP_` de
/// cada copia por secuencia, el `HAP_` de cada fila haplotípica y el `FHA_` de
/// cada asignación por bloque, y comprueba que ningún fundador pertenece a más
/// de una subpoblación. Las copias existen por cada secuencia que el layout
/// puebla: DTC-089 identifica la copia «por secuencia … a través de todos los
/// bloques de esa secuencia», y una secuencia sin bloques no aporta ninguna
/// fila a `f5-founders.tsv` ni ningún genotipo a la proyección de verdad.
///
/// Exige la admisión sellada y la identidad de ejecución porque el layout que
/// aporta las secuencias debe ser el de esta ejecución: cada `block_id` se
/// rederiva y cada bloque debe caber dentro de un tramo elegible de **esta**
/// admisión, con su mismo ordinal de secuencia y sus mismas coordenadas.
///
/// # Errors
///
/// - `GPS-CONF-101` si la admisión sellada no es la de esta configuración, si
///   el escenario, la referencia o la semilla efectiva no son canónicos, si el
///   layout recibido no lo produjo esta ejecución o no cabe en sus tramos
///   elegibles, si la asignación recibida nombra bloques de otra ejecución, si
///   la matriz ganadora no tiene los sitios de su bloque o si un
///   `individual_id` aparece en más de una subpoblación.
/// - `GPS-DOMAIN-002` si la derivación de un `FCP_`, un `HAP_` o un `FHA_`
///   repite un identificador ya registrado.
/// - `GPS-INTERNAL-101` si una unidad asignada nombra un bloque, una matriz,
///   una copia o un rango haplotípico que sus productores no publicaron.
pub fn materialize_f5_founder_identities(
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
    execution_id: &ExecutionId,
    layout: &F5PhysicalLayout,
    search: &F5HaplotypeSearch,
    assignment: &F5FounderAssignment,
) -> Result<F5FounderIdentities, Diagnostic> {
    if !sealed.matches_admission(config, F5OtherResources::default()) {
        return Err(configuration(
            config,
            "preflight.binding",
            "la admisión sellada no corresponde a esta configuración",
        ));
    }
    let scenario_id = ScenarioId::parse(config.scenario_id(), "scenario_id")
        .map_err(|_| configuration(config, "scenario_id", "no es una identidad tipada"))?;
    let reference_id = ReferenceId::parse(config.reference_id(), "reference_id")
        .map_err(|_| configuration(config, "reference_id", "no es una identidad tipada"))?;
    let root_seed = RootSeed::parse("seed", config.seed())
        .map_err(|_| configuration(config, "seed", "la semilla efectiva no es canónica"))?;
    let context = F5IdentityContext::new(root_seed, execution_id, &scenario_id, &reference_id);
    let geometry = index_layout_blocks(config, context, layout, sealed)?;
    // Las secuencias pobladas son el segundo eje de la identidad de copia; el
    // conjunto ordenado las recorre una sola vez por individuo.
    let sequences: BTreeSet<SequenceId> = geometry
        .values()
        .map(|block| block.sequence_id.clone())
        .collect();

    let mut registry = F5IdentityRegistry::default();
    if let Some((individual_id, previous, current)) = conflicting_membership(assignment.rosters()) {
        return Err(configuration(
            config,
            "subpopulations.founder_ids,subpopulations.founder_keys",
            format!(
                "el fundador {} pertenece a {} y a {}, y DTC-089 lo identifica una sola vez en toda la ejecución",
                individual_id.as_str(),
                previous.as_str(),
                current.as_str()
            ),
        ));
    }

    // Los tres bucles recorren el censo en su orden canónico —subpoblaciones
    // por `population_id`, individuos por `individual_id`, secuencias por
    // `sequence_id` y copias por índice—, de modo que `copies` queda ordenado
    // por la clave con la que después se busca cada pertenencia.
    let ploidy = u32::from(config.ploidy());
    let mut copies = Vec::new();
    let mut membership: BTreeMap<FounderCopyId, usize> = BTreeMap::new();
    for roster in assignment.rosters() {
        let population = F5PopulationIdentityContext::new(context, roster.population_id());
        for individual in roster.individuals() {
            for sequence_id in &sequences {
                for copy_index in 0..ploidy {
                    let founder_copy_id = registry.derive_founder_copy(
                        population,
                        F5FounderCopyScope::new(
                            individual.individual_id(),
                            sequence_id,
                            copy_index,
                        ),
                    )?;
                    let position = copies.len();
                    membership.insert(founder_copy_id.clone(), position);
                    copies.push(F5FounderCopy {
                        population_id: roster.population_id().clone(),
                        individual_id: individual.individual_id().clone(),
                        sequence_id: sequence_id.clone(),
                        copy_index,
                        founder_copy_id,
                    });
                }
            }
        }
    }

    // Un índice por bloque×subpoblación; barrer la lista de la búsqueda una vez
    // por unidad recorrería `B*P` filas `B*P` veces.
    let winners: BTreeMap<(&PopulationId, u64), &F5BlockHaplotypes> = search
        .blocks()
        .iter()
        .map(|block| ((block.population_id(), block.block_index()), block))
        .collect();

    let mut haplotypes = Vec::new();
    let mut assignments = Vec::new();
    for unit in assignment.units() {
        let block = geometry.get(&unit.block_index()).ok_or_else(|| {
            internal(
                &scenario_id,
                "una unidad asignada nombra un bloque que el layout no colocó",
            )
        })?;
        // `BLK_` se deriva del índice y del contexto, así que dos piezas de la
        // misma ejecución no pueden discrepar aquí; una asignación construida
        // sobre la búsqueda de **otra** ejecución sí.
        if &block.block_id != unit.block_id() {
            return Err(configuration(
                config,
                "assignment.identity",
                "la asignación recibida no la produjo esta ejecución: su identificador de bloque no coincide con el del layout",
            ));
        }
        let population = F5PopulationIdentityContext::new(context, unit.population_id());
        let matrix = winners
            .get(&(unit.population_id(), unit.block_index()))
            .ok_or_else(|| {
                internal(
                    &scenario_id,
                    "una unidad asignada nombra un bloque×subpoblación que la búsqueda no resolvió",
                )
            })?
            .matrix();
        // La matriz ganadora tiene una columna por sitio de **su** bloque: sin
        // esta comprobación, `alleles` podría medir los sitios de otro y nadie
        // lo notaría, porque el `HAP_` se deriva de la cadena publicada.
        if matrix.site_count() != block.site_count {
            return Err(configuration(
                config,
                "search.geometry",
                format!(
                    "la matriz ganadora del bloque {} tiene {} sitios y su bloque colocado tiene {}",
                    unit.block_index(),
                    matrix.site_count(),
                    block.site_count
                ),
            ));
        }
        let mut rows = Vec::with_capacity(matrix.haplotype_count());
        for rank in 0..matrix.haplotype_count() {
            let alleles = allele_string(matrix, rank).ok_or_else(|| {
                internal(
                    &scenario_id,
                    "la matriz ganadora no conserva la fila de un rango haplotípico",
                )
            })?;
            let haplotype_rank = rank as u64;
            let haplotype_id = registry.derive_haplotype(
                population,
                unit.block_index(),
                haplotype_rank,
                &alleles,
            )?;
            rows.push(haplotype_id.clone());
            haplotypes.push(F5HaplotypeIdentity {
                population_id: unit.population_id().clone(),
                block_id: block.block_id.clone(),
                block_index: unit.block_index(),
                haplotype_rank,
                alleles,
                haplotype_id,
            });
        }
        for seat in unit.assignments() {
            // `copies` ya está en la clave que buscamos —población, individuo,
            // secuencia y copia—, de modo que la pertenencia se localiza sin
            // construir ningún índice paralelo ni clonar la clave.
            let position = copies
                .binary_search_by(|copy| {
                    (
                        &copy.population_id,
                        &copy.individual_id,
                        &copy.sequence_id,
                        copy.copy_index,
                    )
                        .cmp(&(
                            unit.population_id(),
                            seat.individual_id(),
                            &block.sequence_id,
                            seat.copy_index(),
                        ))
                })
                .map_err(|_| {
                    internal(
                        &scenario_id,
                        "una asignación nombra una copia sin pertenencia materializada",
                    )
                })?;
            let founder_copy_id = copies[position].founder_copy_id.clone();
            let haplotype_id = usize::try_from(seat.haplotype_rank())
                .ok()
                .and_then(|rank| rows.get(rank))
                .ok_or_else(|| {
                    internal(
                        &scenario_id,
                        "una asignación nombra un rango haplotípico que la matriz no contiene",
                    )
                })?
                .clone();
            let block_assignment_id = registry.derive_block_haplotype_assignment(
                population,
                F5BlockHaplotypeAssignmentScope::new(
                    unit.block_index(),
                    &founder_copy_id,
                    &haplotype_id,
                ),
            )?;
            assignments.push(F5BlockAssignment {
                population_id: unit.population_id().clone(),
                individual_id: seat.individual_id().clone(),
                sequence_id: block.sequence_id.clone(),
                copy_index: seat.copy_index(),
                block_id: block.block_id.clone(),
                block_index: unit.block_index(),
                founder_copy_id,
                block_assignment_id,
                haplotype_id,
                haplotype_rank: seat.haplotype_rank(),
            });
        }
    }
    // Las unidades llegan por `population_id`+`block_id`; `f5-founders.tsv` se
    // indexa por población, individuo, secuencia, copia y bloque. La clave es
    // única en cada fila, de modo que el orden es total.
    assignments.sort_by(|left, right| {
        (
            &left.population_id,
            &left.individual_id,
            &left.sequence_id,
            left.copy_index,
            &left.block_id,
        )
            .cmp(&(
                &right.population_id,
                &right.individual_id,
                &right.sequence_id,
                right.copy_index,
                &right.block_id,
            ))
    });

    Ok(F5FounderIdentities {
        copies,
        haplotypes,
        assignments,
        membership,
    })
}

/// Índice de los bloques colocados, rederivando cada `block_id`.
///
/// Un índice repetido no necesita guarda propia: la segunda derivación produce
/// el mismo identificador y el registro la rechaza como colisión.
fn index_layout_blocks(
    config: &EffectiveF5Config,
    context: F5IdentityContext<'_>,
    layout: &F5PhysicalLayout,
    sealed: &SealedF5Admission,
) -> Result<BTreeMap<u64, BlockGeometry>, Diagnostic> {
    if u64::try_from(layout.blocks().len()).ok() != Some(sealed.block_count()) {
        return Err(configuration(
            config,
            "preflight.block_count",
            "el layout recibido no pertenece a la admisión sellada de esta configuración",
        ));
    }
    // Tramos elegibles de **esta** admisión, indexados por ordinal de
    // secuencia. El layout coloca cada bloque dentro de uno de ellos y solo los
    // parte al consumirlos, de modo que la contención es exacta: un bloque de
    // otra referencia nombra un ordinal que no existe o unas coordenadas que
    // ningún tramo propio contiene.
    let mut eligible: BTreeMap<u64, Vec<(u64, u64)>> = BTreeMap::new();
    for tract in sealed.eligible_tracts() {
        eligible
            .entry(tract.sequence_ordinal())
            .or_default()
            .push((tract.start_1based(), tract.end_1based()));
    }
    let mut registry = F5IdentityRegistry::default();
    let mut blocks = BTreeMap::new();
    for block in layout.blocks() {
        let expected = registry
            .derive_block(context, block.block_index())
            .map_err(|_| {
                configuration(
                    config,
                    "layout.identity",
                    "no se pudo rederivar el identificador de un bloque colocado",
                )
            })?;
        if &expected != block.block_id() {
            return Err(configuration(
                config,
                "layout.identity",
                "el layout recibido no lo produjo esta ejecución: su identificador de bloque no coincide",
            ));
        }
        let contained = eligible
            .get(&block.sequence_ordinal())
            .is_some_and(|tracts| {
                tracts.iter().any(|(start, end)| {
                    *start <= block.start_1based() && block.end_1based() <= *end
                })
            });
        if !contained {
            return Err(configuration(
                config,
                "layout.geometry",
                format!(
                    "el bloque {} ocupa {}:{}-{}, que no cabe en ningún tramo elegible de esta referencia",
                    block.block_index(),
                    block.sequence_id().as_str(),
                    block.start_1based(),
                    block.end_1based()
                ),
            ));
        }
        blocks.insert(
            block.block_index(),
            BlockGeometry {
                block_id: block.block_id().clone(),
                sequence_id: block.sequence_id().clone(),
                site_count: block.sites().len(),
            },
        );
    }
    Ok(blocks)
}

/// El primer fundador que pertenece a dos subpoblaciones, si lo hay.
///
/// DTC-089: «`IndividualId` identifica al fundador una sola vez en toda la
/// ejecución». El rechazo vive en [`resolve_f5_founder_roster`], que es quien
/// puede emitirlo pronto y con las tres vías cruzadas —literal↔literal,
/// derivado↔derivado y literal↔derivado—; esta función es la red que protege
/// al repositorio de identidades de un censo que no viniera de allí. Sigue
/// existiendo porque una pertenencia duplicada **no** produciría colisión de
/// `FCP_` —`population_id` es componente del material, de modo que los dos
/// identificadores serían distintos— y la proyección de verdad publicaría un
/// portador con dos pertenencias. Como ningún censo admitido puede llegar roto,
/// se prueba de forma unitaria, sobre censos construidos a mano.
fn conflicting_membership(
    rosters: &[F5FounderRoster],
) -> Option<(&IndividualId, &PopulationId, &PopulationId)> {
    let mut seen: BTreeMap<&IndividualId, &PopulationId> = BTreeMap::new();
    for roster in rosters {
        for individual in roster.individuals() {
            if let Some(previous) = seen.insert(individual.individual_id(), &roster.population_id) {
                return Some((individual.individual_id(), previous, &roster.population_id));
            }
        }
    }
    None
}

/// Fila de la matriz ganadora como los `S` caracteres `0`/`1` que DTC-089 hace
/// componente del `HAP_`.
fn allele_string(matrix: &F5AlleleMatrix, haplotype: usize) -> Option<String> {
    let mut alleles = String::with_capacity(matrix.site_count());
    for site in 0..matrix.site_count() {
        alleles.push(if matrix.allele(haplotype, site)? {
            '1'
        } else {
            '0'
        });
    }
    Some(alleles)
}

// ---------------------------------------------------------------------------
// Diagnósticos
// ---------------------------------------------------------------------------

fn configuration(
    config: &EffectiveF5Config,
    fields: impl Into<String>,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::f5_configuration(config.scenario_path().display().to_string(), fields, detail)
}

fn internal(scenario_id: &ScenarioId, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::f5_internal(scenario_id.as_str(), detail)
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;

    use super::{
        CopySeat, F5CopyAssignment, F5FounderIndividual, F5FounderRoster, F5RandomDecision,
        F5SamplingMode, F5SamplingResolution, F5SamplingSource, FounderInbreeding, HaplotypePrefix,
        NeumaierSum, QuotaError, RankState, RootSeed, conflicting_membership, derive_f5_subseed,
        draw_copy_rank, plan_controlled_quota, realized_diploid_f, resolve_founder_inbreeding,
        sample_probabilistic_unit, select_tie_prefix,
    };
    use crate::diagnostics::Diagnostic;
    use crate::identity::{F5BlockId, IndividualId, PopulationId, ScenarioId};
    use crate::randomness::{LocalRng, RandomWordSource};

    /// Flujo real de `controlled_quota_tie` para las unitarias del desempate.
    fn tie_rng() -> LocalRng {
        let root = RootSeed::parse(
            "seed",
            "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8",
        )
        .expect("semilla");
        let scenario = ScenarioId::parse("scenario", "test").expect("escenario");
        let population = PopulationId::parse("P1", "test").expect("población");
        let block = F5BlockId::parse("BLK_TEST", "test").expect("bloque");
        derive_f5_subseed(
            root,
            F5RandomDecision::controlled_quota_tie(&scenario, &population, &block),
        )
        .expect("subsemilla")
        .local_rng()
    }

    /// Fuente que entrega palabras fijadas y cuenta cuántas se consumieron.
    struct SpySource {
        words: Vec<u64>,
        consumed: usize,
    }

    impl SpySource {
        fn new(words: &[u64]) -> Self {
            Self {
                words: words.to_vec(),
                consumed: 0,
            }
        }
    }

    impl RandomWordSource for SpySource {
        fn next_word(&mut self) -> u64 {
            let word = self.words.get(self.consumed).copied().unwrap_or(0);
            self.consumed += 1;
            word
        }
    }

    // -----------------------------------------------------------------------
    // Oráculo entero exacto del reparto por mayores restos
    // -----------------------------------------------------------------------

    /// Estado esperado de cada rango y grupo fronterizo, calculado con enteros
    /// de 128 bits y sin llamar a ninguna función del productor.
    struct OraclePlan {
        floors: Vec<u64>,
        states: Vec<RankState>,
        boundary_ranks: Vec<u64>,
        select_count: u64,
        /// `true` si **todos** los estados aplicables respetan la cota.
        feasible: bool,
    }

    /// `numerators/denominator` son los pesos exactos y `tolerance` el
    /// racional `tn/td`. La cota se comprueba como
    /// `|n*D - C*w| * td <= tn * C * D`, sin dividir nunca.
    fn oracle_plan(
        numerators: &[u128],
        denominator: u128,
        copies: u128,
        tolerance: (u128, u128),
    ) -> OraclePlan {
        let mut floors = Vec::with_capacity(numerators.len());
        let mut remainders = Vec::with_capacity(numerators.len());
        let mut base = 0_u128;
        for numerator in numerators {
            let ideal = copies * numerator;
            floors.push((ideal / denominator) as u64);
            remainders.push(ideal % denominator);
            base += ideal / denominator;
        }
        let residual = copies - base;
        let mut order: Vec<usize> = (0..numerators.len()).collect();
        order.sort_by(|left, right| {
            remainders[*right]
                .cmp(&remainders[*left])
                .then_with(|| left.cmp(right))
        });
        let mut states = vec![RankState::Floor; numerators.len()];
        let mut boundary_ranks = Vec::new();
        let mut select_count = 0_u64;
        let mut start = 0_usize;
        while start < order.len() {
            let remainder = remainders[order[start]];
            let mut end = start + 1;
            while end < order.len() && remainders[order[end]] == remainder {
                end += 1;
            }
            if (end as u128) <= residual {
                for rank in &order[start..end] {
                    states[*rank] = RankState::FloorPlusOne;
                }
            } else if (start as u128) < residual {
                for rank in &order[start..end] {
                    states[*rank] = RankState::Boundary;
                }
                boundary_ranks = order[start..end].iter().map(|rank| *rank as u64).collect();
                boundary_ranks.sort_unstable();
                select_count = (residual - start as u128) as u64;
            }
            start = end;
        }
        let (tolerance_numerator, tolerance_denominator) = tolerance;
        let mut feasible = true;
        for (rank, numerator) in numerators.iter().enumerate() {
            let counts: Vec<u128> = match states[rank] {
                RankState::Floor => vec![u128::from(floors[rank])],
                RankState::FloorPlusOne => vec![u128::from(floors[rank]) + 1],
                RankState::Boundary => {
                    vec![u128::from(floors[rank]), u128::from(floors[rank]) + 1]
                }
            };
            for count in counts {
                let deviation = (count * denominator) as i128 - (copies * numerator) as i128;
                if deviation.unsigned_abs() * tolerance_denominator
                    > tolerance_numerator * copies * denominator
                {
                    feasible = false;
                }
            }
        }
        OraclePlan {
            floors,
            states,
            boundary_ranks,
            select_count,
            feasible,
        }
    }

    /// Los cuatro conjuntos de pesos diádicos de la batería, como pesos
    /// binary64 exactos y como racionales exactos sobre un denominador común.
    fn weight_sets() -> Vec<(Vec<f64>, Vec<u128>, u128)> {
        vec![
            (vec![0.5, 0.5], vec![1, 1], 2),
            (vec![0.75, 0.25], vec![3, 1], 4),
            (vec![0.25, 0.25, 0.25, 0.25], vec![1, 1, 1, 1], 4),
            (vec![0.5, 0.25, 0.125, 0.125], vec![4, 2, 1, 1], 8),
            (
                vec![0.5, 0.125, 0.125, 0.125, 0.125],
                vec![4, 1, 1, 1, 1],
                8,
            ),
            (vec![0.125; 8], vec![1, 1, 1, 1, 1, 1, 1, 1], 8),
        ]
    }

    /// El reparto, sus estados, su grupo fronterizo y su **veredicto** de
    /// factibilidad coinciden con el oráculo entero para todo `C` de la
    /// barrida, sin saltarse ninguno.
    ///
    /// Con pesos diádicos `C*q_h` es exacto en binary64 y la cota se compara
    /// en copias, donde `n - C*q_h` también lo es: la igualdad con el racional
    /// `|n*D - C*w| * 2C <= C * D` es exacta para todo `C`, no solo para las
    /// potencias de dos. La versión anterior de esta prueba dejaba los `C` no
    /// potencia de dos dentro de un `if let Ok(...)`, de modo que un rechazo
    /// no comprobaba nada; el recuento `compared` lo fija ahora.
    #[test]
    fn the_largest_remainder_split_matches_the_integer_oracle() {
        let mut compared = 0_u64;
        let mut with_boundary = 0_u64;
        let mut rejected = 0_u64;
        for (weights, numerators, denominator) in weight_sets() {
            for copies in 1_u64..=32 {
                let observed = plan_controlled_quota(&weights, copies, 0.5);
                let expected = oracle_plan(
                    &numerators,
                    denominator,
                    u128::from(copies),
                    (1, u128::from(2 * copies)),
                );
                match observed {
                    Ok(plan) => {
                        assert!(
                            expected.feasible,
                            "C={copies}: el oráculo la declara inviable"
                        );
                        assert_eq!(plan.floors, expected.floors, "C={copies}");
                        assert_eq!(plan.states, expected.states, "C={copies}");
                        assert_eq!(plan.boundary_ranks, expected.boundary_ranks, "C={copies}");
                        assert_eq!(plan.select_count, expected.select_count, "C={copies}");
                        let total: u64 = plan
                            .floors
                            .iter()
                            .zip(&plan.states)
                            .map(|(floor, state)| {
                                floor + u64::from(*state == RankState::FloorPlusOne)
                            })
                            .sum();
                        assert_eq!(
                            total + plan.select_count,
                            copies,
                            "C={copies}: los fijos más el prefijo deben sumar C"
                        );
                        if !plan.boundary_ranks.is_empty() {
                            with_boundary += 1;
                            assert!(
                                plan.select_count > 0
                                    && plan.select_count < plan.boundary_ranks.len() as u64,
                                "C={copies}: un grupo fronterizo cruza la frontera"
                            );
                        }
                    }
                    Err(QuotaError::Infeasible(detail)) => {
                        assert!(!expected.feasible, "C={copies}: {detail}");
                        rejected += 1;
                    }
                    Err(QuotaError::Internal(detail)) => {
                        panic!("C={copies}: invariante interna rota: {detail}")
                    }
                }
                compared += 1;
            }
        }
        assert_eq!(
            compared,
            weight_sets().len() as u64 * 32,
            "toda la barrida debe comparar veredicto y estructura, sin saltarse ningún C"
        );
        assert!(
            rejected >= 5,
            "la barrida debe incluir cuotas que el oráculo declara inviables"
        );
        assert!(
            with_boundary >= 5,
            "la barrida debe incluir grupos fronterizos"
        );
    }

    /// Todos los restos iguales y un grupo que abarca el residual entero.
    #[test]
    fn a_single_group_can_span_the_whole_residual() {
        let plan = plan_controlled_quota(&[0.25, 0.25, 0.25, 0.25], 2, 0.5)
            .unwrap_or_else(|_| panic!("cuota factible"));
        assert_eq!(plan.floors, vec![0, 0, 0, 0]);
        assert_eq!(plan.states, vec![RankState::Boundary; 4]);
        assert_eq!(plan.boundary_ranks, vec![0, 1, 2, 3]);
        assert_eq!(plan.select_count, 2, "el residual entero cae en el grupo");
        assert_eq!(plan.boundary_remainder, 0.5);
    }

    /// Un grupo incluido o excluido por completo no deja frontera, y sin
    /// frontera el desempate no abre su ámbito RNG.
    #[test]
    fn groups_included_or_excluded_in_full_leave_no_boundary() {
        // `C=8`: cada ideal es entero y el residual es cero.
        let plan = plan_controlled_quota(&[0.5, 0.25, 0.125, 0.125], 8, 0.5)
            .unwrap_or_else(|_| panic!("cuota factible"));
        assert_eq!(plan.floors, vec![4, 2, 1, 1]);
        assert!(
            plan.boundary_ranks.is_empty(),
            "no hay residual que repartir"
        );
        assert_eq!(plan.select_count, 0);
        assert!(plan.states.iter().all(|state| *state == RankState::Floor));

        // `C=2`: el residual lo absorbe entero un grupo de un solo rango.
        let plan = plan_controlled_quota(&[0.5, 0.25, 0.125, 0.125], 2, 0.5)
            .unwrap_or_else(|_| panic!("cuota factible"));
        assert_eq!(plan.floors, vec![1, 0, 0, 0]);
        assert_eq!(
            plan.states,
            vec![
                RankState::Floor,
                RankState::FloorPlusOne,
                RankState::Floor,
                RankState::Floor
            ]
        );
        assert!(
            plan.boundary_ranks.is_empty(),
            "un grupo incluido por completo no cruza la frontera"
        );
    }

    /// La factibilidad exige los **dos** estados del grupo fronterizo: basta
    /// que falle el piso más uno para rechazar, aunque el piso quepa.
    #[test]
    fn the_boundary_group_is_rejected_when_only_one_of_its_two_states_fits() {
        let weights = [0.25, 0.25, 0.25, 0.25];
        // `C=5`: cada ideal vale 1,25 y el residual alcanza a un solo rango.
        let plan = plan_controlled_quota(&weights, 5, 0.5);
        let detail = match plan {
            Err(QuotaError::Infeasible(detail)) => detail,
            Err(QuotaError::Internal(detail)) => panic!("invariante interna: {detail}"),
            Ok(_) => panic!("el estado piso más uno no cabe en la cota"),
        };
        assert!(
            detail.contains("|2-1.25|"),
            "el rechazo debe citar el estado que no cabe: {detail}"
        );
        // Control positivo: con una tolerancia que admite los dos estados
        // —`0,75` copias cubre `|2-1,25|`— la misma cuota es factible y deja
        // su frontera empatada. Comprobar solo el piso aceptaría las dos.
        let plan = plan_controlled_quota(&weights, 5, 0.75)
            .unwrap_or_else(|_| panic!("con 0,75 copias los dos estados caben"));
        assert_eq!(plan.boundary_ranks, vec![0, 1, 2, 3]);
        assert_eq!(plan.select_count, 1);
    }

    /// Un resto nulo nunca recibe el incremento: el residual se agota antes
    /// de llegar a él, aunque haya una frontera empatada por delante.
    #[test]
    fn a_null_weight_keeps_a_zero_quota() {
        let plan = plan_controlled_quota(&[0.5, 0.25, 0.25, 0.0], 2, 0.5)
            .unwrap_or_else(|_| panic!("cuota factible"));
        assert_eq!(plan.floors, vec![1, 0, 0, 0]);
        assert_eq!(plan.states[3], RankState::Floor, "el peso nulo no compite");
        assert_eq!(plan.boundary_ranks, vec![1, 2], "la frontera no lo incluye");
        assert_eq!(plan.select_count, 1);
    }

    /// `C=0` no es una cuota que repartir, sino una invariante rota.
    #[test]
    fn zero_copies_are_an_internal_invariant() {
        assert!(matches!(
            plan_controlled_quota(&[0.5, 0.5], 0, 0.5),
            Err(QuotaError::Internal(_))
        ));
    }

    /// La factibilidad no puede consumir RNG: `plan_controlled_quota` no
    /// recibe ninguna fuente y el rechazo es, por firma, anterior a todo
    /// sorteo. La prueba fija la firma —cualquier parámetro aleatorio futuro
    /// rompería esta llamada— y el rechazo que ocurre sin ella.
    #[test]
    fn the_feasibility_check_rejects_without_any_random_source() {
        let outcome = plan_controlled_quota(&[0.25, 0.25, 0.25, 0.25], 5, 0.5);
        assert!(matches!(outcome, Err(QuotaError::Infeasible(_))));
    }

    /// Un grupo que no cruza la frontera **no abre** el ámbito
    /// `controlled_quota_tie`: el cierre que derivaría su subsemilla no llega
    /// a invocarse, que es la regla exacta del contrato.
    #[test]
    fn a_group_without_a_boundary_never_opens_its_scope() {
        let plan = plan_controlled_quota(&[0.5, 0.25, 0.125, 0.125], 8, 0.5)
            .unwrap_or_else(|_| panic!("cuota factible"));
        assert!(plan.boundary_ranks.is_empty());
        let (tie, operations) = select_tie_prefix(&plan, || {
            panic!("el ámbito controlled_quota_tie no puede abrirse sin frontera")
        })
        .unwrap_or_else(|_| panic!("sin frontera no hay desempate"));
        assert!(tie.is_none());
        assert_eq!(operations, 0, "sin frontera no se consume ningún entero");

        // Con frontera sí se abre, y consume `g-1` enteros.
        let plan = plan_controlled_quota(&[0.5, 0.25, 0.125, 0.125], 4, 0.5)
            .unwrap_or_else(|_| panic!("cuota factible"));
        assert_eq!(plan.boundary_ranks, vec![2, 3]);
        let mut opened = false;
        let (tie, operations) = select_tie_prefix(&plan, || {
            opened = true;
            Ok(tie_rng())
        })
        .unwrap_or_else(|_| panic!("desempate completo"));
        assert!(opened, "la frontera empatada sí abre su ámbito");
        assert_eq!(operations, 1, "un grupo de dos consume un entero");
        let tie = tie.expect("frontera empatada");
        assert_eq!(tie.candidate_ranks(), &[2, 3]);
        assert_eq!(tie.selected_ranks().len(), 1);
    }

    // -----------------------------------------------------------------------
    // F diploide
    // -----------------------------------------------------------------------

    /// Escenario de las unitarias de F, que solo necesitan un identificador.
    fn scenario_id() -> ScenarioId {
        ScenarioId::parse("scenario", "test").expect("escenario")
    }

    /// Individuo de las unitarias de `F_realized`.
    fn individual(label: &str) -> IndividualId {
        IndividualId::parse(label.to_owned(), "test").expect("individuo")
    }

    /// Visitas que cuesta el paseo de prefijos con una palabra dada, para
    /// separar del coste total el U53 que añade `founder_f`.
    fn walk_cost(prefix: &HaplotypePrefix, word: u64) -> u64 {
        let mut spy = SpySource::new(&[word]);
        prefix.choose(&mut spy).1
    }

    /// Asignación suelta, para construir censos sintéticos.
    fn assignment(label: &str, copy_index: u32, haplotype_rank: u64) -> F5CopyAssignment {
        F5CopyAssignment {
            individual_id: individual(label),
            copy_index,
            haplotype_rank,
        }
    }

    /// Declarar F con `P=1` o `P>2` es error, y no declararlo nunca lo es.
    ///
    /// REQ-EVO-48 rechaza *configurar* F fuera de los diploides, no solo
    /// aplicarlo: por eso la entrada es el F **declarado** y no un cero que no
    /// distinguiría «sin declarar» de «declarado 0».
    #[test]
    fn declaring_f_outside_a_diploid_is_an_error() {
        for ploidy in [1_u8, 3, 4, 10] {
            for declared in [0.0, 0.25, 1.0] {
                let detail = resolve_founder_inbreeding(ploidy, Some(declared))
                    .expect_err("F fuera de P=2 debe rechazarse");
                assert!(
                    detail.contains("ploidía 2") && detail.contains(&ploidy.to_string()),
                    "el rechazo debe citar ámbito y ploidía: {detail}"
                );
            }
            assert_eq!(
                resolve_founder_inbreeding(ploidy, None),
                Ok(FounderInbreeding::Independent),
                "sin F declarado, cualquier ploidía muestrea copias independientes"
            );
        }
    }

    /// Los tres consumos del catálogo: `F=0` y `F=1` no abren `founder_f`, y
    /// solo `0<F<1` lo consulta.
    #[test]
    fn the_diploid_f_resolves_its_three_consumption_regimes() {
        assert_eq!(
            resolve_founder_inbreeding(2, Some(0.0)),
            Ok(FounderInbreeding::Independent)
        );
        assert_eq!(
            resolve_founder_inbreeding(2, Some(1.0)),
            Ok(FounderInbreeding::AlwaysReuse)
        );
        assert_eq!(
            resolve_founder_inbreeding(2, Some(0.25)),
            Ok(FounderInbreeding::Bernoulli(0.25))
        );
        assert_eq!(
            resolve_founder_inbreeding(2, None),
            Ok(FounderInbreeding::Independent)
        );
        for invalid in [-0.25_f64, 1.5, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                resolve_founder_inbreeding(2, Some(invalid)).is_err(),
                "F={invalid} no pertenece a [0,1]"
            );
        }
    }

    /// La copia 0 nunca consulta `founder_f`, cualquiera que sea F.
    #[test]
    fn the_first_copy_never_opens_the_f_scope() {
        let prefix = HaplotypePrefix::new(&[0.5, 0.5]).expect("prefijos");
        for inbreeding in [
            FounderInbreeding::Independent,
            FounderInbreeding::Bernoulli(0.5),
            FounderInbreeding::AlwaysReuse,
        ] {
            let (rank, operations) = draw_copy_rank(
                &scenario_id(),
                inbreeding,
                None,
                &prefix,
                || -> Result<SpySource, Diagnostic> {
                    panic!("la copia 0 no puede abrir founder_f")
                },
                || Ok(SpySource::new(&[u64::MAX])),
            )
            .expect("copia 0 muestreada");
            assert_eq!(rank, 1, "la palabra máxima elige el último rango");
            assert_eq!(
                operations,
                walk_cost(&prefix, u64::MAX),
                "solo el paseo de prefijos, sin el U53 de founder_f"
            );
        }
    }

    /// `0<F<1`: `founder_f` se abre **antes** que el `founder_sample` de la
    /// copia 1 y, cuando decide reutilizar, ese segundo flujo no se abre.
    ///
    /// El orden no se ve en ninguna salida —los dos ámbitos derivan
    /// subsemillas independientes—, de modo que solo esta observación directa
    /// lo fija.
    #[test]
    fn the_f_scope_is_consumed_before_the_second_copy_sample() {
        let prefix = HaplotypePrefix::new(&[0.5, 0.5]).expect("prefijos");
        let threshold = (0.5_f64 * (1_u64 << 53) as f64).floor() as u64;

        // Bits53 por debajo del umbral: reutiliza y `founder_sample` no se
        // abre. El rango es el de la copia 0, no el que el prefijo elegiría.
        let opened = RefCell::new(Vec::new());
        let (rank, operations) = draw_copy_rank(
            &scenario_id(),
            FounderInbreeding::Bernoulli(0.5),
            Some(0),
            &prefix,
            || {
                opened.borrow_mut().push("founder_f");
                Ok(SpySource::new(&[(threshold - 1) << 11]))
            },
            || -> Result<SpySource, Diagnostic> {
                opened.borrow_mut().push("founder_sample");
                panic!("reutilizar no puede abrir el segundo flujo")
            },
        )
        .expect("copia 1 reutilizada");
        assert_eq!(rank, 0, "la copia 1 reutiliza exactamente la fila de la 0");
        assert_eq!(operations, 1, "un U53 de founder_f y ningún paseo");
        assert_eq!(*opened.borrow(), vec!["founder_f"]);

        // Bits53 en el umbral: no reutiliza y abre su propio muestreo,
        // siempre después de `founder_f`.
        let opened = RefCell::new(Vec::new());
        let (rank, operations) = draw_copy_rank(
            &scenario_id(),
            FounderInbreeding::Bernoulli(0.5),
            Some(0),
            &prefix,
            || {
                opened.borrow_mut().push("founder_f");
                Ok(SpySource::new(&[threshold << 11]))
            },
            || {
                opened.borrow_mut().push("founder_sample");
                Ok(SpySource::new(&[u64::MAX]))
            },
        )
        .expect("copia 1 muestreada");
        assert_eq!(rank, 1, "el muestreo independiente eligió su propia fila");
        assert_eq!(
            operations,
            1 + walk_cost(&prefix, u64::MAX),
            "el U53 de F y el paseo de prefijos"
        );
        assert_eq!(*opened.borrow(), vec!["founder_f", "founder_sample"]);
    }

    /// `F=0` y `F=1` no abren `founder_f`: el cierre que derivaría su
    /// subsemilla no llega a invocarse, que es la regla Bernoulli literal del
    /// catálogo.
    #[test]
    fn the_extreme_f_values_never_open_the_f_scope() {
        let prefix = HaplotypePrefix::new(&[0.5, 0.5]).expect("prefijos");

        // `F=0`: muestreo independiente de la copia 1.
        let (rank, operations) = draw_copy_rank(
            &scenario_id(),
            FounderInbreeding::Independent,
            Some(0),
            &prefix,
            || -> Result<SpySource, Diagnostic> { panic!("F=0 no consume founder_f") },
            || Ok(SpySource::new(&[u64::MAX])),
        )
        .expect("copia 1 muestreada");
        assert_eq!(rank, 1, "con F=0 la copia 1 se muestrea por su cuenta");
        assert_eq!(
            operations,
            walk_cost(&prefix, u64::MAX),
            "F=0 no añade el U53 de founder_f"
        );

        // `F=1`: reutilización completa sin abrir ninguno de los dos flujos.
        let (rank, operations) = draw_copy_rank(
            &scenario_id(),
            FounderInbreeding::AlwaysReuse,
            Some(1),
            &prefix,
            || -> Result<SpySource, Diagnostic> { panic!("F=1 no consume founder_f") },
            || -> Result<SpySource, Diagnostic> { panic!("F=1 no abre el muestreo de la copia 1") },
        )
        .expect("copia 1 reutilizada");
        assert_eq!(rank, 1, "homocigosis completa");
        assert_eq!(operations, 0, "ningún sorteo");
    }

    /// Vector de la Bernoulli de F: la decisión es exactamente
    /// `bits53 < floor(F*2^53)`, con el borde inferior y el superior.
    #[test]
    fn the_f_bernoulli_follows_the_catalogue_vector() {
        let prefix = HaplotypePrefix::new(&[0.5, 0.5]).expect("prefijos");
        let mut compared = 0_u64;
        let mut reused = 0_u64;
        for probability in [0.25_f64, 0.5, 0.75, 1.0 / 3.0] {
            let threshold = (probability * (1_u64 << 53) as f64).floor() as u64;
            for bits in [
                0,
                1,
                threshold - 1,
                threshold,
                threshold + 1,
                (1_u64 << 53) - 1,
            ] {
                let expected_reuse = bits < threshold;
                let (rank, operations) = draw_copy_rank(
                    &scenario_id(),
                    FounderInbreeding::Bernoulli(probability),
                    Some(0),
                    &prefix,
                    || Ok(SpySource::new(&[bits << 11])),
                    // La palabra máxima elige el rango 1, distinto del 0 que
                    // reutilizaría: la fila distingue las dos decisiones.
                    || Ok(SpySource::new(&[u64::MAX])),
                )
                .expect("copia 1 resuelta");
                assert_eq!(
                    rank == 0,
                    expected_reuse,
                    "F={probability} con bits53={bits} debe {} reutilizar",
                    if expected_reuse { "" } else { "no" }
                );
                assert_eq!(
                    operations,
                    if expected_reuse {
                        1
                    } else {
                        1 + walk_cost(&prefix, u64::MAX)
                    },
                    "F={probability} con bits53={bits}: el U53 de F se consume siempre"
                );
                reused += u64::from(expected_reuse);
                compared += 1;
            }
        }
        assert_eq!(compared, 24, "la rejilla completa debe compararse");
        assert!(
            reused > 0 && reused < compared,
            "la rejilla cubre ambas decisiones"
        );
    }

    /// Semilla raíz de las unitarias que llaman al muestreo completo.
    fn unit_root_seed() -> RootSeed {
        RootSeed::parse(
            "seed",
            "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8",
        )
        .expect("semilla")
    }

    /// Fila probabilística sintética de `copies` copias.
    fn probabilistic_row(copies: u64) -> F5SamplingResolution {
        F5SamplingResolution {
            population_id: PopulationId::parse("P1", "test").expect("población"),
            block_id: F5BlockId::parse("BLK_TEST", "test").expect("bloque"),
            block_index: 0,
            requested_mode: F5SamplingMode::Probabilistic,
            effective_mode: F5SamplingMode::Probabilistic,
            source: F5SamplingSource::SubpopulationDefault,
            controlled_max_abs_error_declared: None,
            controlled_max_abs_error_effective: None,
            controlled_max_abs_error_copies: None,
            copies,
        }
    }

    /// Fila que una copia sortearía por su cuenta, reconstruida desde su
    /// ámbito con el recorrido secuencial y sin llamar al productor.
    fn own_rank(row: &F5SamplingResolution, weights: &[f64], label: &str, copy_index: u32) -> u64 {
        let mut rng = derive_f5_subseed(
            unit_root_seed(),
            F5RandomDecision::founder_sample(
                &scenario_id(),
                &row.population_id,
                &row.block_id,
                &individual(label),
                copy_index,
            ),
        )
        .expect("subsemilla")
        .local_rng();
        let threshold = (rng.next_u64() >> 11) as f64 / (1_u64 << 53) as f64;
        sequential_walk(weights, threshold)
    }

    /// Solo la copia 1 **del mismo individuo** puede reutilizar una fila.
    ///
    /// El recorrido las entrega juntas y en orden, de modo que la comprobación
    /// de identidad parece redundante; lo que esta prueba fija es que, si el
    /// orden cambiara, una copia 1 no heredaría la fila de otro individuo.
    #[test]
    fn only_the_second_copy_of_the_same_individual_can_reuse_a_row() {
        let weights = [0.5_f64, 0.5];
        let prefix = HaplotypePrefix::new(&weights).expect("prefijos");
        let root_seed = unit_root_seed();
        let scenario = scenario_id();
        let row = probabilistic_row(2);

        // Pareja ajena cuyas filas independientes difieren: si la copia 1
        // heredase la del asiento anterior, la diferencia se vería.
        let first = own_rank(&row, &weights, "I0", 0);
        let (label, expected) = ["I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8"]
            .into_iter()
            .map(|label| (label, own_rank(&row, &weights, label, 1)))
            .find(|(_, rank)| *rank != first)
            .expect("alguna copia 1 ajena sortea otra fila");
        let seats = vec![
            CopySeat {
                individual_id: individual("I0"),
                copy_index: 0,
            },
            CopySeat {
                individual_id: individual(label),
                copy_index: 1,
            },
        ];
        // `F=1` reutilizaría siempre si la identidad no se comprobara.
        let (assignments, _) = sample_probabilistic_unit(
            &scenario,
            root_seed,
            &row,
            &seats,
            &prefix,
            FounderInbreeding::AlwaysReuse,
        )
        .expect("unidad muestreada");
        assert_eq!(assignments[0].haplotype_rank, first);
        assert_eq!(
            assignments[1].haplotype_rank, expected,
            "la copia 1 de otro individuo se muestrea por su cuenta"
        );
    }

    /// Solo la copia **1** hereda: con `P=3`, la copia 2 se muestrea por su
    /// cuenta aunque el individuo sea el mismo y la reutilización esté activa.
    ///
    /// La vía productiva no llega aquí —F exige `P=2`—, pero sin la condición
    /// sobre `copy_index` la reutilización se propagaría copia a copia en
    /// cuanto otra ploidía la admitiera.
    #[test]
    fn only_the_second_copy_index_can_inherit_a_row() {
        let weights = [0.5_f64, 0.5];
        let prefix = HaplotypePrefix::new(&weights).expect("prefijos");
        let scenario = scenario_id();
        let row = probabilistic_row(3);

        // Individuo cuya copia 2 sortearía una fila distinta de la copia 0, que
        // es la que heredaría si la condición no existiera.
        let (label, expected) = ["I0", "I1", "I2", "I3", "I4", "I5", "I6", "I7"]
            .into_iter()
            .map(|label| (label, own_rank(&row, &weights, label, 2)))
            .find(|(label, rank)| *rank != own_rank(&row, &weights, label, 0))
            .expect("alguna copia 2 sortea otra fila que su copia 0");
        let seats = (0..3)
            .map(|copy_index| CopySeat {
                individual_id: individual(label),
                copy_index,
            })
            .collect::<Vec<_>>();
        let (assignments, _) = sample_probabilistic_unit(
            &scenario,
            unit_root_seed(),
            &row,
            &seats,
            &prefix,
            FounderInbreeding::AlwaysReuse,
        )
        .expect("unidad muestreada");
        assert_eq!(
            assignments[1].haplotype_rank, assignments[0].haplotype_rank,
            "la copia 1 sí hereda"
        );
        assert_eq!(
            assignments[2].haplotype_rank, expected,
            "la copia 2 se muestrea por su cuenta"
        );
    }

    // -----------------------------------------------------------------------
    // F realizado
    // -----------------------------------------------------------------------

    /// `F_realized` bialélico contrastado con su racional exacto, con los tres
    /// regímenes: homocigosis completa, exceso de heterocigotos y mezcla.
    ///
    /// El oráculo evalúa `1 - het*C²/(N*(C² - Σ n_h²))`, una expresión
    /// distinta de la del productor y exacta en enteros para estos censos.
    #[test]
    fn the_realized_f_matches_its_rational_oracle() {
        let cases: [(&[(u64, u64)], f64); 4] = [
            // Dos individuos homocigotos de rangos distintos: H_obs=0.
            (&[(0, 0), (1, 1)], 1.0),
            // Dos heterocigotos: H_obs=1 y H_exp=1/2, F=-1.
            (&[(0, 1), (0, 1)], -1.0),
            // Mitad y mitad: H_obs=1/2, H_exp=1/2, F=0.
            (&[(0, 0), (0, 1), (1, 1), (0, 1)], 0.0),
            // Cuatro individuos con tres rangos: F=40/168.
            (&[(0, 0), (1, 2), (0, 1), (2, 2)], 1.0 - 128.0 / 168.0),
        ];
        for (individuals, expected) in cases {
            let mut assignments = Vec::new();
            for (index, (left, right)) in individuals.iter().enumerate() {
                let label = format!("I{index}");
                assignments.push(assignment(&label, 0, *left));
                assignments.push(assignment(&label, 1, *right));
            }
            let (value, operations) =
                realized_diploid_f(&assignments, 3).expect("censo diploide bien formado");
            let value = value.expect("H_expected positivo");
            assert_eq!(
                operations,
                individuals.len() as u64 + 3,
                "el recuento visita cada individuo y cada rango, y lo contabiliza"
            );
            assert!(
                (value - expected).abs() <= 1e-12,
                "F realizado {value} frente al racional {expected}"
            );
        }
    }

    /// Una sola categoría realizada deja `H_expected=0` y el estado es no
    /// aplicable, **no** cero: el cociente no llega a evaluarse.
    #[test]
    fn a_null_expected_heterozygosity_is_not_applicable() {
        let assignments = vec![
            assignment("I0", 0, 2),
            assignment("I0", 1, 2),
            assignment("I1", 0, 2),
            assignment("I1", 1, 2),
        ];
        assert_eq!(
            realized_diploid_f(&assignments, 4)
                .expect("censo bien formado")
                .0,
            None,
            "con una sola categoría no hay heterocigosis esperada"
        );
        // Control positivo: basta mover una copia a otro rango para que la
        // misma fórmula sí sea evaluable.
        let mut mixed = assignments.clone();
        mixed[3] = assignment("I1", 1, 0);
        assert!(
            realized_diploid_f(&mixed, 4)
                .expect("censo bien formado")
                .0
                .is_some(),
            "con dos categorías el denominador deja de ser nulo"
        );
    }

    /// Sin copias no hay frecuencia realizada: la función es total y devuelve
    /// no aplicable en lugar de dividir por cero.
    #[test]
    fn an_empty_unit_is_not_applicable() {
        assert_eq!(
            realized_diploid_f(&[], 2).expect("censo vacío"),
            (None, 0),
            "sin copias no se visita nada: el `NaN` del cociente vacío daría el \
             mismo `None` pero recorrería los rangos"
        );
    }

    /// Un censo diploide mal emparejado es una invariante rota, no un cálculo:
    /// sin esta comprobación, dos copias 0 de individuos distintos contarían
    /// como un individuo heterocigoto inexistente.
    #[test]
    fn a_broken_diploid_pairing_is_rejected_instead_of_counted() {
        assert!(
            realized_diploid_f(&[assignment("I0", 0, 0), assignment("I1", 0, 1)], 2).is_err(),
            "dos copias 0 de individuos distintos no son un individuo"
        );
        assert!(
            realized_diploid_f(&[assignment("I0", 1, 0), assignment("I0", 0, 1)], 2).is_err(),
            "las copias deben venir en orden 0 y 1"
        );
        assert!(
            realized_diploid_f(&[assignment("I0", 0, 0)], 2).is_err(),
            "una copia suelta no completa un individuo"
        );
        assert!(
            realized_diploid_f(&[assignment("I0", 0, 0), assignment("I0", 1, 7)], 2).is_err(),
            "un rango fuera de la matriz no puede contarse"
        );
    }

    // -----------------------------------------------------------------------
    // Muestreo probabilístico
    // -----------------------------------------------------------------------

    /// Recorrido secuencial literal de DTC-089, escrito desde su enunciado:
    /// suma compensada en orden de rango y primer prefijo estrictamente mayor
    /// que el umbral, con el último rango positivo como residuo.
    ///
    /// Es el oráculo de [`HaplotypePrefix`], que hace lo mismo pero con la
    /// tabla calculada una vez por unidad y una partición binaria por copia.
    fn sequential_walk(weights: &[f64], threshold: f64) -> u64 {
        let mut compensated = NeumaierSum::default();
        let mut last_positive = 0_u64;
        for (rank, weight) in weights.iter().enumerate() {
            if *weight > 0.0 {
                last_positive = rank as u64;
            }
            if threshold < compensated.add(*weight).expect("suma finita") {
                return rank as u64;
            }
        }
        last_positive
    }

    /// La tabla de prefijos elige exactamente el mismo rango que el recorrido
    /// secuencial, para toda la rejilla de umbrales y todos los perfiles.
    #[test]
    fn the_prefix_table_agrees_with_the_sequential_walk() {
        let sets = [
            vec![0.5, 0.5],
            vec![0.75, 0.25],
            vec![0.25, 0.25, 0.25, 0.25],
            vec![0.5, 0.25, 0.125, 0.125],
            vec![0.5, 0.0, 0.5],
            vec![0.2, 0.2, 0.2, 0.2, 0.2],
            vec![0.1, 0.2, 0.3, 0.4],
        ];
        let mut compared = 0_u64;
        for weights in &sets {
            let prefix = HaplotypePrefix::new(weights).expect("prefijos");
            // El U53 pertenece a `[0,1)`: el umbral nunca alcanza 1.
            for step in 0..512_u64 {
                let threshold = step as f64 / 512.0;
                let bits = (threshold * (1_u64 << 53) as f64) as u64;
                let mut source = SpySource::new(&[bits << 11]);
                let (observed, _) = prefix.choose(&mut source);
                assert_eq!(
                    observed,
                    sequential_walk(weights, threshold),
                    "umbral {threshold} con {weights:?}"
                );
                compared += 1;
            }
        }
        assert_eq!(compared, 7 * 512, "la rejilla completa debe compararse");
    }

    /// Una sola alternativa no consume RNG; con dos o más, exactamente un U53.
    #[test]
    fn the_prefix_walk_consumes_one_word_per_copy() {
        let prefix = HaplotypePrefix::new(&[1.0]).expect("prefijos");
        let mut spy = SpySource::new(&[u64::MAX]);
        assert_eq!(prefix.choose(&mut spy).0, 0);
        assert_eq!(spy.consumed, 0, "una sola alternativa no sortea");

        let prefix = HaplotypePrefix::new(&[0.5, 0.5]).expect("prefijos");
        let mut spy = SpySource::new(&[0]);
        assert_eq!(prefix.choose(&mut spy).0, 0);
        assert_eq!(spy.consumed, 1, "un U53 por muestreo independiente");

        let mut spy = SpySource::new(&[u64::MAX]);
        assert_eq!(prefix.choose(&mut spy).0, 1);
        assert_eq!(spy.consumed, 1);
    }

    /// El prefijo compensado decide la frontera exacta entre dos rangos.
    #[test]
    fn the_prefix_walk_switches_rank_at_the_exact_threshold() {
        let prefix = HaplotypePrefix::new(&[0.5, 0.5]).expect("prefijos");
        // `U = bits53 / 2^53`; con `bits53 = 2^52` el umbral vale 0,5 y ya no
        // es estrictamente menor que el primer prefijo.
        let mut spy = SpySource::new(&[((1_u64 << 52) - 1) << 11]);
        assert_eq!(prefix.choose(&mut spy).0, 0);
        let mut spy = SpySource::new(&[(1_u64 << 52) << 11]);
        assert_eq!(
            prefix.choose(&mut spy).0,
            1,
            "el borde exacto pertenece al rango siguiente"
        );
    }

    /// Un `q_h` nulo es inalcanzable, incluso con la palabra máxima, y el
    /// residuo de redondeo lo absorbe el último rango con peso positivo.
    #[test]
    fn a_null_weight_is_unreachable_by_the_prefix_walk() {
        let prefix = HaplotypePrefix::new(&[0.5, 0.0, 0.5]).expect("prefijos");
        for word in [0, 1, 1 << 32, u64::MAX] {
            let mut spy = SpySource::new(&[word]);
            assert_ne!(
                prefix.choose(&mut spy).0,
                1,
                "el peso nulo no puede recibir una copia"
            );
        }
        // Pesos que no alcanzan 1: el residuo cae en el último positivo, no en
        // la última fila.
        let prefix = HaplotypePrefix::new(&[0.25, 0.25, 0.0]).expect("prefijos");
        let mut spy = SpySource::new(&[u64::MAX]);
        assert_eq!(prefix.choose(&mut spy).0, 1);
    }

    /// Un peso no finito o negativo, o un perfil sin ningún peso positivo, es
    /// una invariante rota y no una elección.
    #[test]
    fn an_invalid_weight_is_rejected_instead_of_sampled() {
        assert!(HaplotypePrefix::new(&[f64::NAN, 1.0]).is_err());
        assert!(HaplotypePrefix::new(&[-0.5, 1.5]).is_err());
        assert!(HaplotypePrefix::new(&[0.0, 0.0]).is_err());
    }

    /// La tabla se compone una sola vez por unidad: su coste es `H`, no `C*H`.
    #[test]
    fn the_prefix_table_costs_at_most_a_logarithm_per_copy() {
        let weights = vec![0.0625; 16];
        let prefix = HaplotypePrefix::new(&weights).expect("prefijos");
        for word in [0, 1 << 40, u64::MAX] {
            let mut spy = SpySource::new(&[word]);
            let (_, visits) = prefix.choose(&mut spy);
            assert!(
                visits <= 5,
                "16 rangos exigen como mucho 5 comparaciones, no 16: {visits}"
            );
        }
    }

    /// Censo de prueba: la pertenencia se comprueba sobre censos construidos a
    /// mano porque ninguna configuración admitida puede producir uno roto.
    fn roster(population: &str, individuals: &[&str]) -> F5FounderRoster {
        let population_id = PopulationId::parse(population, "test").expect("población");
        let individuals = individuals
            .iter()
            .map(|individual| F5FounderIndividual {
                individual_id: IndividualId::parse(*individual, "test").expect("individuo"),
                founder_key: None,
            })
            .collect::<Vec<_>>();
        F5FounderRoster {
            population_id,
            individual_count: individuals.len() as u64,
            copies: individuals.len() as u64 * 2,
            individuals,
        }
    }

    /// Un fundador en dos subpoblaciones es pertenencia duplicada: sus `FCP_`
    /// no colisionarían, porque `population_id` es componente del material.
    #[test]
    fn a_founder_in_two_subpopulations_is_a_duplicated_membership() {
        let rosters = [
            roster("P1", &["alpha", "beta"]),
            roster("P2", &["beta", "gamma"]),
        ];
        let (individual_id, previous, current) =
            conflicting_membership(&rosters).expect("la pertenencia duplicada se detecta");
        assert_eq!(individual_id.as_str(), "beta");
        assert_eq!(previous.as_str(), "P1");
        assert_eq!(current.as_str(), "P2");
    }

    /// Un fundador repetido dentro de un mismo censo también es pertenencia
    /// duplicada, y ahí las dos subpoblaciones nombradas coinciden.
    #[test]
    fn a_founder_repeated_inside_one_roster_is_also_rejected() {
        let rosters = [roster("P1", &["alpha", "alpha"])];
        let (individual_id, previous, current) =
            conflicting_membership(&rosters).expect("la repetición local se detecta");
        assert_eq!(individual_id.as_str(), "alpha");
        assert_eq!(previous.as_str(), "P1");
        assert_eq!(current.as_str(), "P1");
    }

    /// Censos disjuntos no producen ningún conflicto.
    #[test]
    fn disjoint_rosters_have_no_membership_conflict() {
        let rosters = [
            roster("P1", &["alpha", "beta"]),
            roster("P2", &["gamma", "delta"]),
        ];
        assert!(conflicting_membership(&rosters).is_none());
    }
}
