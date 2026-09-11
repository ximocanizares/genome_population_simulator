//! F5.3.1. Modo por bloque×subpoblación, cuotas controladas y muestreo
//! probabilístico de las copias fundadoras.
//!
//! DTC-089 «Fundadores, F y factibilidad controlada» fija tres cosas que este
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
//!
//! Lo que **no** hace, por pertenecer a otros hitos: F diploide y `F_realized`
//! son de F5.3.2 —aquí solo se rechaza la combinación `F>0`×controlado, que
//! DTC-089 declara incompatible—, y las identidades `FCP_`/`FHA_` y la
//! pertenencia publicada son de F5.3.3.

use std::collections::{BTreeMap, BTreeSet};

use crate::application::EffectiveF5Config;
use crate::configuration::{DeclaredF5SamplingOverride, DeclaredF5Subpopulation};
use crate::diagnostics::Diagnostic;
use crate::f5_marginals::NeumaierSum;
use crate::f5_preflight::{F5OtherResources, SealedF5Admission};
use crate::f5_search::F5HaplotypeSearch;
use crate::identity::{
    AutomaticIdRequest, F5BlockId, FounderKey, IdClass, IdentityRegistry, IndividualId,
    PopulationId, ScenarioId,
};
use crate::randomness::{
    F5RandomDecision, LocalRng, RandomWordSource, RootSeed, derive_f5_subseed,
    f5_fisher_yates_with_source, f5_uniform_53,
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

/// Resuelve modos y tolerancias y asigna una fila haplotípica a cada copia.
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
///   contrato, si la búsqueda no cubre exactamente las unidades del escenario
///   o si una subpoblación controlada declara `F>0`, combinación que DTC-089
///   declara incompatible.
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
    let mut seats_by_population: BTreeMap<PopulationId, Vec<CopySeat>> = BTreeMap::new();
    for roster in &rosters {
        seats_by_population.insert(
            roster.population_id.clone(),
            copy_seats(roster, config.ploidy()),
        );
    }

    let mut units = Vec::with_capacity(sampling.len());
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
        let (assignments, quota, operations) = match row.effective_mode {
            F5SamplingMode::Probabilistic => {
                let (assignments, operations) =
                    sample_probabilistic_unit(&scenario_id, root_seed, row, seats, &prefix)?;
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
    })
}

/// Censo fundador de cada subpoblación, ordenado por `population_id`.
///
/// Un `founder_id` aportado se usa literalmente como `IndividualId`; una
/// `founder_key` deriva el suyo con el contrato población×clave de REQ-IND-25,
/// el mismo que ya usa la vía F1/F2, sin ordinal físico.
///
/// # Errors
///
/// `GPS-CONF-101` si la subpoblación no conserva exactamente un conjunto de
/// fundadores, si un identificador no es válido, si `C` desborda `u64` o si la
/// derivación de una clave colisiona.
pub fn resolve_f5_founder_roster(
    config: &EffectiveF5Config,
) -> Result<Vec<F5FounderRoster>, Diagnostic> {
    let root_seed = RootSeed::parse("seed", config.seed())
        .map_err(|_| configuration(config, "seed", "la semilla efectiva no es canónica"))?;
    let ploidy = u64::from(config.ploidy());
    let mut registry = IdentityRegistry::default();
    let mut rosters = Vec::with_capacity(config.declared_scenario().subpopulations.len());
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
                let inbreeding = effective_inbreeding_f(config, subpopulation)?;
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

/// F efectivo de una subpoblación: el suyo, el de `[f5]` o `0.0`.
///
/// Solo se usa para rechazar la combinación incompatible con el modo
/// controlado; aplicar F y publicar `F_realized` pertenece a F5.3.2.
fn effective_inbreeding_f(
    config: &EffectiveF5Config,
    subpopulation: &DeclaredF5Subpopulation,
) -> Result<f64, Diagnostic> {
    let declared = subpopulation
        .inbreeding_f
        .as_ref()
        .or(config.declared_scenario().inbreeding_f.as_ref());
    let Some(number) = declared else {
        return Ok(0.0);
    };
    number.literal().parse::<f64>().map_err(|_| {
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

fn sample_probabilistic_unit(
    scenario_id: &ScenarioId,
    root_seed: RootSeed,
    row: &F5SamplingResolution,
    seats: &[CopySeat],
    prefix: &HaplotypePrefix,
) -> Result<(Vec<F5CopyAssignment>, u64), Diagnostic> {
    let mut assignments = Vec::with_capacity(seats.len());
    let mut operations = 0_u64;
    // Las copias de un individuo se recorren juntas y en orden de
    // `copy_index`: F5.3.2 añadirá aquí la reutilización de la fila de la
    // copia 0 con probabilidad F sin reordenar ningún consumo.
    for seat in seats {
        operations = operations.saturating_add(1);
        let mut rng = derive_f5_subseed(
            root_seed,
            F5RandomDecision::founder_sample(
                scenario_id,
                &row.population_id,
                &row.block_id,
                &seat.individual_id,
                seat.copy_index,
            ),
        )?
        .local_rng();
        let (haplotype_rank, visits) = prefix.choose(&mut rng);
        operations = operations.saturating_add(visits);
        assignments.push(F5CopyAssignment {
            individual_id: seat.individual_id.clone(),
            copy_index: seat.copy_index,
            haplotype_rank,
        });
    }
    Ok((assignments, operations))
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
    use super::{
        F5RandomDecision, HaplotypePrefix, NeumaierSum, QuotaError, RankState, RootSeed,
        derive_f5_subseed, plan_controlled_quota, select_tie_prefix,
    };
    use crate::identity::{F5BlockId, PopulationId, ScenarioId};
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
}
