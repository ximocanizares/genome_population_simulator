//! Orquestación de casos de uso mediante contratos del dominio.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::configuration::{
    DeclaredConfig, DeclaredGeneralConfig, DeclaredReductionTarget, DeclaredReferenceSource,
    DeclaredScenarioConfig, ExpectedConfigKind, FrequencyMode, Route, SequenceClass,
};
use crate::diagnostics::Diagnostic;
use crate::exact_decimal::ExactDecimal;
use crate::identity::{
    BlockId, CrossId, EventId, GroupId, IdentityRegistry, IndividualId, PhaseId, PopulationId,
    ReferenceId, ScenarioId, SequenceId, SubpopulationId,
};
use crate::numeric::{NormalizedCollection, normalize_declared};
use crate::randomness::{
    EntropySource, ReproducibilityMetadata, Subseed, SubseedRequest, SystemEntropy, derive_subseed,
    resolve_root_seed,
};
use crate::rules::{
    AxisKind, ParameterKey, ResolutionContext, Rule, RuleProvenance, RuleSet, RuleValue,
    SelectorIdentityCatalog, SelectorSet, SpatialLevel, metadata, validate_rules,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ValidatedGeneralSummary {
    pub total_genome_size: u64,
    pub explicit_event_ids: Vec<String>,
    pub individual_ids: Vec<String>,
    pub sequence_classes: BTreeMap<String, SequenceClass>,
}

/// Secuencia completamente validada y poseída por la configuración efectiva.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveSequence {
    id: String,
    class: SequenceClass,
    length: u64,
    ploidy: u8,
}

impl EffectiveSequence {
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub const fn class(&self) -> SequenceClass {
        self.class
    }

    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }

    #[must_use]
    pub const fn ploidy(&self) -> u8 {
        self.ploidy
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveCount {
    id: String,
    category: String,
    count: u64,
}

impl EffectiveCount {
    #[must_use]
    pub fn identity(&self) -> (&str, &str) {
        (&self.id, &self.category)
    }

    #[must_use]
    pub const fn count(&self) -> u64 {
        self.count
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveControlledSet {
    id: String,
    counts: Vec<EffectiveCount>,
}

/// Origen de referencia ya resuelto respecto al fichero general admitido.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EffectiveReferenceSource {
    Simulated,
    Supplied { original_path: PathBuf },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EffectiveReductionTarget {
    Bases { declared: u64 },
    ReductionPercent { literal: String },
}

impl EffectiveReductionTarget {
    #[must_use]
    pub const fn mode(&self) -> &'static str {
        match self {
            Self::Bases { .. } => "target_bases",
            Self::ReductionPercent { .. } => "reduction_percent",
        }
    }

    #[must_use]
    pub fn declared_value(&self) -> String {
        match self {
            Self::Bases { declared } => declared.to_string(),
            Self::ReductionPercent { literal } => literal.clone(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveReducedReference {
    id: String,
    source_reference_id: String,
    region_size_bases: u64,
    target: EffectiveReductionTarget,
    target_bases: u64,
    selection_subseed: Subseed,
}

impl EffectiveReducedReference {
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn source_reference_id(&self) -> &str {
        &self.source_reference_id
    }

    #[must_use]
    pub const fn region_size_bases(&self) -> u64 {
        self.region_size_bases
    }

    #[must_use]
    pub const fn target(&self) -> &EffectiveReductionTarget {
        &self.target
    }

    #[must_use]
    pub const fn target_bases(&self) -> u64 {
        self.target_bases
    }

    #[must_use]
    pub const fn selection_subseed(&self) -> Subseed {
        self.selection_subseed
    }
}

/// Entrada de configuración ya admitida, conservada byte a byte para que la
/// ejecución no vuelva a leer una fuente mutable después de validarla.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveConfigurationInput {
    kind: &'static str,
    id: String,
    original_path: String,
    bytes: Arc<[u8]>,
}

impl EffectiveConfigurationInput {
    #[must_use]
    pub const fn kind(&self) -> &'static str {
        self.kind
    }

    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn original_path(&self) -> &str {
        &self.original_path
    }

    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

impl EffectiveControlledSet {
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn counts(&self) -> &[EffectiveCount] {
        &self.counts
    }
}

/// Objeto y contexto cuya configuración debe resolverse antes de consumirla.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ResolutionRequest {
    pub object_id: String,
    pub parameter: ParameterKey,
    pub context: ResolutionContext,
}

/// Relación población-subpoblación validada por la fase que posee ese dominio.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PopulationMembership {
    pub population: String,
    pub subpopulation: String,
}

/// Intervalo cerrado de generaciones validado por la fase evolutiva propietaria.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PhaseGenerationInterval {
    pub phase: String,
    pub start: u32,
    pub end: u32,
}

/// Locus y coordenada validados por la fase que posee el catálogo de loci.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DomainLocus {
    pub id: String,
    pub reference: String,
    pub sequence: String,
    pub position: u64,
}

/// Prueba inmutable de relaciones que F1 no puede inferir desde listas de IDs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ResolutionDomain {
    provider: String,
    population_memberships: BTreeMap<String, String>,
    phase_generation_intervals: BTreeMap<String, (u32, u32)>,
    loci: BTreeMap<String, DomainLocus>,
}

impl ResolutionDomain {
    #[must_use]
    pub fn identity_only() -> Self {
        Self {
            provider: "f1-identity-only".to_owned(),
            population_memberships: BTreeMap::new(),
            phase_generation_intervals: BTreeMap::new(),
            loci: BTreeMap::new(),
        }
    }

    pub fn validated(
        provider: impl Into<String>,
        populations: impl IntoIterator<Item = PopulationMembership>,
        phases: impl IntoIterator<Item = PhaseGenerationInterval>,
        loci: impl IntoIterator<Item = DomainLocus>,
    ) -> Result<Self, String> {
        let provider = provider.into();
        if provider.is_empty() {
            return Err("el proveedor del dominio de resolución está vacío".to_owned());
        }
        let mut population_memberships = BTreeMap::new();
        for membership in populations {
            if membership.population.is_empty()
                || membership.subpopulation.is_empty()
                || population_memberships
                    .insert(membership.subpopulation, membership.population)
                    .is_some()
            {
                return Err("relación población-subpoblación vacía o repetida".to_owned());
            }
        }
        let mut phase_generation_intervals = BTreeMap::new();
        for interval in phases {
            if interval.phase.is_empty()
                || interval.start > interval.end
                || phase_generation_intervals
                    .insert(interval.phase, (interval.start, interval.end))
                    .is_some()
            {
                return Err("intervalo de fase vacío, invertido o repetido".to_owned());
            }
        }
        let mut locus_catalog = BTreeMap::new();
        for locus in loci {
            if locus.id.is_empty()
                || locus.reference.is_empty()
                || locus.sequence.is_empty()
                || locus.position == 0
                || locus_catalog.insert(locus.id.clone(), locus).is_some()
            {
                return Err("locus vacío, inválido o repetido".to_owned());
            }
        }
        Ok(Self {
            provider,
            population_memberships,
            phase_generation_intervals,
            loci: locus_catalog,
        })
    }

    #[must_use]
    pub fn provider(&self) -> &str {
        &self.provider
    }

    #[must_use]
    pub const fn population_memberships(&self) -> &BTreeMap<String, String> {
        &self.population_memberships
    }

    #[must_use]
    pub const fn phase_generation_intervals(&self) -> &BTreeMap<String, (u32, u32)> {
        &self.phase_generation_intervals
    }

    #[must_use]
    pub const fn loci(&self) -> &BTreeMap<String, DomainLocus> {
        &self.loci
    }
}

/// Valor efectivo con la identidad y procedencia completas de la regla ganadora.
#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveValue {
    object_id: String,
    parameter: ParameterKey,
    resolution_context: ResolutionContext,
    value: RuleValue,
    winning_rule_id: String,
    selectors: SelectorSet,
    provenance: RuleProvenance,
}

impl EffectiveValue {
    #[must_use]
    pub fn object_id(&self) -> &str {
        &self.object_id
    }

    #[must_use]
    pub const fn parameter(&self) -> ParameterKey {
        self.parameter
    }

    #[must_use]
    pub const fn resolution_context(&self) -> &ResolutionContext {
        &self.resolution_context
    }

    #[must_use]
    pub const fn value(&self) -> &RuleValue {
        &self.value
    }

    #[must_use]
    pub fn winning_rule_id(&self) -> &str {
        &self.winning_rule_id
    }

    #[must_use]
    pub const fn selectors(&self) -> &SelectorSet {
        &self.selectors
    }

    #[must_use]
    pub const fn provenance(&self) -> &RuleProvenance {
        &self.provenance
    }
}

/// Instantánea poseída, validada y resuelta. No expone ninguna mutación.
#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveConfig {
    schema_version: &'static str,
    source_path: String,
    route: Route,
    reference_id: String,
    reference_source: EffectiveReferenceSource,
    reduced_references: Vec<EffectiveReducedReference>,
    total_genome_size: u64,
    sequences: Vec<EffectiveSequence>,
    scenario_ids: Vec<String>,
    configuration_inputs: Vec<EffectiveConfigurationInput>,
    explicit_event_ids: Vec<String>,
    individual_ids: Vec<String>,
    normalized: Vec<NormalizedCollection>,
    controlled: Vec<EffectiveControlledSet>,
    values: Vec<EffectiveValue>,
    resolution_domain: ResolutionDomain,
    resolution_requests: Vec<ResolutionRequest>,
    reproducibility: ReproducibilityMetadata,
}

impl EffectiveConfig {
    #[must_use]
    pub const fn schema_version(&self) -> &'static str {
        self.schema_version
    }

    #[must_use]
    pub fn source_path(&self) -> &str {
        &self.source_path
    }

    #[must_use]
    pub const fn route(&self) -> Route {
        self.route
    }

    #[must_use]
    pub fn reference_id(&self) -> &str {
        &self.reference_id
    }

    #[must_use]
    pub const fn reference_source(&self) -> &EffectiveReferenceSource {
        &self.reference_source
    }

    #[must_use]
    pub fn reduced_references(&self) -> &[EffectiveReducedReference] {
        &self.reduced_references
    }

    #[must_use]
    pub const fn total_genome_size(&self) -> u64 {
        self.total_genome_size
    }

    #[must_use]
    pub fn sequences(&self) -> &[EffectiveSequence] {
        &self.sequences
    }

    #[must_use]
    pub fn scenario_ids(&self) -> &[String] {
        &self.scenario_ids
    }

    #[must_use]
    pub fn configuration_inputs(&self) -> &[EffectiveConfigurationInput] {
        &self.configuration_inputs
    }

    #[must_use]
    pub fn explicit_event_ids(&self) -> &[String] {
        &self.explicit_event_ids
    }

    #[must_use]
    pub fn individual_ids(&self) -> &[String] {
        &self.individual_ids
    }

    #[must_use]
    pub fn normalized_collections(&self) -> &[NormalizedCollection] {
        &self.normalized
    }

    #[must_use]
    pub fn controlled_sets(&self) -> &[EffectiveControlledSet] {
        &self.controlled
    }

    #[must_use]
    pub fn values(&self) -> &[EffectiveValue] {
        &self.values
    }

    #[must_use]
    pub const fn resolution_domain(&self) -> &ResolutionDomain {
        &self.resolution_domain
    }

    /// Solicitudes autoritativas que originaron los valores efectivos.
    #[must_use]
    pub fn resolution_requests(&self) -> &[ResolutionRequest] {
        &self.resolution_requests
    }

    #[must_use]
    pub const fn reproducibility(&self) -> &ReproducibilityMetadata {
        &self.reproducibility
    }
}

/// Frontera de solo lectura usada por el subcomando `validate`. No crea ni
/// publica un destino, incluso cuando el llamador ha indicado uno.
pub fn load_and_validate_config(path: &Path) -> Result<EffectiveConfig, Diagnostic> {
    load_and_validate_config_with_entropy(path, &mut SystemEntropy)
}

pub fn load_and_validate_config_with_entropy(
    path: &Path,
    entropy: &mut impl EntropySource,
) -> Result<EffectiveConfig, Diagnostic> {
    let bytes = std::fs::read(path).map_err(|error| input_error(path, error))?;
    let admitted = admit_toml_header(path, bytes, ExpectedConfigKind::General)?;
    let DeclaredConfig::General(general) = parse_declared_config(admitted)? else {
        return Err(config_error_at(
            path,
            "config_kind",
            "la entrada principal no es una configuración general",
        ));
    };

    let mut scenarios = Vec::with_capacity(general.scenarios().len());
    for link in general.scenarios() {
        let scenario_path = resolve_link(path, link.path.value());
        let bytes =
            std::fs::read(&scenario_path).map_err(|error| input_error(&scenario_path, error))?;
        let admitted = admit_toml_header(&scenario_path, bytes, ExpectedConfigKind::Scenario)?;
        let DeclaredConfig::Scenario(scenario) = parse_declared_config(admitted)? else {
            return Err(config_error_at(
                &scenario_path,
                "config_kind",
                "una entrada de escenario no es configuración de escenario",
            ));
        };
        scenarios.push(scenario);
    }
    let scenario_refs = scenarios.iter().map(Box::as_ref).collect::<Vec<_>>();
    build_effective_config_with_entropy(&general, &scenario_refs, &[], entropy)
}

fn resolve_link(general_path: &Path, declared: &str) -> PathBuf {
    let declared = Path::new(declared);
    if declared.is_absolute() {
        declared.to_path_buf()
    } else {
        general_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(declared)
    }
}

fn input_error(path: &Path, error: std::io::Error) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-001",
        "No se puede abrir o leer una entrada declarada.",
    )
    .with("path", path.display().to_string())
    .with("detail", error.to_string())
}

/// Valida relaciones cruzadas sin iniciar simulación ni crear un destino.
pub fn validate_general_config(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
) -> Result<ValidatedGeneralSummary, Diagnostic> {
    validate_route_and_scenarios(general, scenarios)?;

    let _reference_id = ReferenceId::parse(general.reference().id.value(), "execution")?;
    let mut reference_ids = BTreeSet::from([general.reference().id.value().clone()]);
    for reduction in general.reduced_references() {
        let _ = ReferenceId::parse(reduction.id.value(), "execution")?;
        if !reference_ids.insert(reduction.id.value().clone()) {
            return Err(config_error(
                general,
                "reduced_references.id",
                format!(
                    "identificador de referencia repetido: {}",
                    reduction.id.value()
                ),
            ));
        }
        if reduction.source_reference_id.value() != general.reference().id.value() {
            return Err(config_error(
                general,
                "reduced_references.source_reference_id",
                format!(
                    "la referencia fuente {} no coincide con la referencia suministrada {}",
                    reduction.source_reference_id.value(),
                    general.reference().id.value()
                ),
            ));
        }
        if !matches!(
            general.reference().source,
            DeclaredReferenceSource::Supplied { .. }
        ) {
            return Err(config_error(
                general,
                "reduced_references.source_reference_id",
                "una reducción solo puede usar una referencia FASTA suministrada",
            ));
        }
        if *reduction.region_size_bases.value() == 0 {
            return Err(config_error(
                general,
                "reduced_references.region_size_bases",
                "region_size_bases debe ser mayor que cero",
            ));
        }
        if matches!(&reduction.target, DeclaredReductionTarget::Bases(value) if *value.value() == 0)
        {
            return Err(config_error(
                general,
                "reduced_references.target_bases",
                "target_bases debe ser mayor que cero",
            ));
        }
    }
    for sequence in &general.reference().sequences {
        let _ = SequenceId::parse(sequence.id.value(), general.reference().id.value())?;
    }
    for scenario in general.scenarios() {
        let _ = ScenarioId::parse(scenario.id.value(), "execution")?;
    }
    let mut phase_owners = BTreeMap::<String, (String, String)>::new();
    for phase in general.phases() {
        let _ = PhaseId::parse(phase.id().value(), "execution")?;
        phase_owners.insert(
            phase.id().value().clone(),
            (
                "general".to_owned(),
                general.source().path().display().to_string(),
            ),
        );
    }
    for scenario in scenarios {
        let scenario_id = ScenarioId::parse(scenario.id().value(), "execution")?;
        for phase in scenario.phases() {
            let _ = PhaseId::parse(phase.id().value(), scenario_id.as_str())?;
            if let Some((first_owner, first_path)) = phase_owners.insert(
                phase.id().value().clone(),
                (
                    scenario_id.as_str().to_owned(),
                    scenario.source().path().display().to_string(),
                ),
            ) {
                return Err(config_error(
                    general,
                    "phases.id",
                    "los identificadores de fase son globales en v1",
                )
                .with("phase", phase.id().value())
                .with("first_scenario", first_owner)
                .with("first_path", first_path)
                .with("second_scenario", scenario_id.as_str())
                .with(
                    "second_path",
                    scenario.source().path().display().to_string(),
                ));
            }
        }
        for cross in scenario.crosses() {
            let _ = CrossId::parse(cross.id().value(), scenario_id.as_str())?;
        }
    }
    for population in general.populations() {
        let _ = PopulationId::parse(population.id().value(), "execution")?;
    }
    for subpopulation in general.subpopulations() {
        let _ = SubpopulationId::parse(subpopulation.id().value(), "execution")?;
    }
    for group in general.groups() {
        let _ = GroupId::parse(group.id().value(), "execution")?;
    }
    for block in general.blocks() {
        let _ = BlockId::parse(block.id().value(), "execution")?;
    }
    let mut identities = IdentityRegistry::default();
    let mut individual_ids = Vec::with_capacity(general.individuals().len());
    for individual in general.individuals() {
        let id = IndividualId::parse(individual.id().value(), "execution")?;
        individual_ids.push(identities.register_individual("execution", id)?.to_string());
    }
    individual_ids.sort();

    let declared_events = general
        .events()
        .iter()
        .chain(scenarios.iter().flat_map(|scenario| scenario.events()));
    let mut explicit_event_ids = Vec::new();
    for event in declared_events {
        let id = EventId::parse(event.id().value(), "execution")?;
        explicit_event_ids.push(identities.register_event("execution", id)?.to_string());
    }
    explicit_event_ids.sort();

    let mut total_genome_size = 0_u64;
    let mut sequence_classes = BTreeMap::new();
    let default_ploidy = match &general.reference().source {
        DeclaredReferenceSource::Simulated { default_ploidy, .. } => Some(*default_ploidy.value()),
        DeclaredReferenceSource::Supplied { .. } => None,
    };
    if general.reference().sequences.is_empty() {
        return Err(config_error(
            general,
            "reference.sequences",
            "la referencia debe declarar al menos una secuencia",
        ));
    }
    for sequence in &general.reference().sequences {
        if *sequence.length.value() == 0 {
            return Err(config_error(
                general,
                "reference.sequences.length",
                format!("la secuencia {} tiene longitud cero", sequence.id.value()),
            ));
        }
        total_genome_size = total_genome_size
            .checked_add(*sequence.length.value())
            .ok_or_else(|| {
                config_error(
                    general,
                    "reference.sequences.length",
                    "el tamaño total del genoma desborda u64",
                )
            })?;
        if sequence_classes
            .insert(sequence.id.value().clone(), *sequence.class.value())
            .is_some()
        {
            return Err(config_error(
                general,
                "reference.sequences.id",
                format!(
                    "identificador de secuencia repetido: {}",
                    sequence.id.value()
                ),
            ));
        }
        let effective_ploidy = sequence
            .ploidy
            .as_ref()
            .map(|value| *value.value())
            .or(default_ploidy);
        if !matches!(effective_ploidy, Some(1..=10)) {
            return Err(Diagnostic::configuration(
                "GPS-CONF-014",
                "Una secuencia carece de ploidía efectiva válida.",
            )
            .with("reference", general.reference().id.value())
            .with("sequence", sequence.id.value())
            .with(
                "found",
                effective_ploidy.map_or_else(|| "missing".to_owned(), |value| value.to_string()),
            )
            .with("accepted", "1..=10"));
        }
    }

    validate_declared_rules(general, scenarios)?;

    Ok(ValidatedGeneralSummary {
        total_genome_size,
        explicit_event_ids,
        individual_ids,
        sequence_classes,
    })
}

/// Construye una instantánea efectiva sin conservar referencias al TOML ni a
/// la configuración declarada. Las solicitudes representan objetos realmente
/// evaluados; para cada parámetro declarado sin solicitud se registra además
/// su valor base global.
pub fn build_effective_config(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
    requests: &[ResolutionRequest],
) -> Result<EffectiveConfig, Diagnostic> {
    build_effective_config_with_domain(
        general,
        scenarios,
        requests,
        &ResolutionDomain::identity_only(),
    )
}

/// Resuelve solicitudes contra relaciones previamente validadas por su fase
/// propietaria. F1 no construye ni infiere esas relaciones.
pub fn build_effective_config_with_domain(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
    requests: &[ResolutionRequest],
    domain: &ResolutionDomain,
) -> Result<EffectiveConfig, Diagnostic> {
    build_effective_config_with_entropy_and_domain(
        general,
        scenarios,
        requests,
        &mut SystemEntropy,
        domain,
    )
}

pub fn build_effective_config_with_entropy(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
    requests: &[ResolutionRequest],
    entropy: &mut impl EntropySource,
) -> Result<EffectiveConfig, Diagnostic> {
    build_effective_config_with_entropy_and_domain(
        general,
        scenarios,
        requests,
        entropy,
        &ResolutionDomain::identity_only(),
    )
}

pub(crate) fn build_effective_config_with_entropy_and_domain(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
    requests: &[ResolutionRequest],
    entropy: &mut impl EntropySource,
    domain: &ResolutionDomain,
) -> Result<EffectiveConfig, Diagnostic> {
    let summary = validate_general_config(general, scenarios)?;
    validate_resolution_domain(general, scenarios, domain)?;
    let (root_seed, seed_origin) = resolve_root_seed(
        general.seed().map(|declared| declared.value().as_str()),
        entropy,
    )
    .map_err(|diagnostic| {
        let mut diagnostic = diagnostic.with("path", general.source().path().display().to_string());
        if let Some(seed) = general.seed() {
            diagnostic = diagnostic.with("location", seed.location().to_string());
        }
        diagnostic
    })?;
    let default_ploidy = match &general.reference().source {
        DeclaredReferenceSource::Simulated { default_ploidy, .. } => Some(*default_ploidy.value()),
        DeclaredReferenceSource::Supplied { .. } => None,
    };
    let mut sequences = general
        .reference()
        .sequences
        .iter()
        .map(|sequence| {
            let ploidy = sequence
                .ploidy
                .as_ref()
                .map(|value| *value.value())
                .or(default_ploidy)
                .ok_or_else(|| {
                    config_error(
                        general,
                        "reference.sequences.ploidy",
                        "la secuencia validada carece de ploidía efectiva",
                    )
                })?;
            Ok(EffectiveSequence {
                id: sequence.id.value().clone(),
                class: *sequence.class.value(),
                length: *sequence.length.value(),
                ploidy,
            })
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;

    let mut normalized = Vec::new();
    if let DeclaredReferenceSource::Simulated { composition, .. } = &general.reference().source {
        normalized.push(normalize_reference_composition(general, composition)?);
    }
    let mut controlled = Vec::new();
    for set in general.frequency_sets() {
        match *set.mode.value() {
            FrequencyMode::Probabilistic => normalized.push(normalize_declared(
                &format!("frequency_sets.{}", set.id.value()),
                &set.values,
                set.tolerance
                    .as_ref()
                    .or_else(|| general.normalization_tolerance()),
            )?),
            FrequencyMode::Controlled => {
                let mut counts = set
                    .counts
                    .iter()
                    .map(|count| EffectiveCount {
                        id: count.id.value().clone(),
                        category: count.category.value().clone(),
                        count: *count.count.value(),
                    })
                    .collect::<Vec<_>>();
                counts.sort_by(|first, second| {
                    (&first.id, &first.category).cmp(&(&second.id, &second.category))
                });
                controlled.push(EffectiveControlledSet {
                    id: set.id.value().clone(),
                    counts,
                });
            }
        }
    }
    normalized.sort_by(|first, second| first.parameter.cmp(&second.parameter));
    controlled.sort_by(|first, second| first.id.cmp(&second.id));

    let grouped = group_rules(general.rules());
    let rule_sets = grouped
        .into_iter()
        .map(|(parameter, rules)| {
            RuleSet::new_with_fallback(parameter, rules, parameter == ParameterKey::Ploidy)
                .map(|set| (parameter, set))
                .map_err(|error| error.with("path", general.source().path().display().to_string()))
        })
        .collect::<Result<BTreeMap<_, _>, _>>()?;
    let mut resolution_requests = requests.to_vec();
    resolution_requests.sort_by(|first, second| {
        (&first.object_id, first.parameter.canonical_rank())
            .cmp(&(&second.object_id, second.parameter.canonical_rank()))
    });
    let mut effective_requests = resolution_requests.clone();
    if let Some(ploidy_rules) = rule_sets.get(&ParameterKey::Ploidy) {
        for sequence in &sequences {
            let context = ResolutionContext {
                reference: Some(general.reference().id.value().clone()),
                sequence: Some(sequence.id.clone()),
                ..ResolutionContext::default()
            };
            if !ploidy_rules.matching_rules(&context).is_empty() {
                effective_requests.push(ResolutionRequest {
                    object_id: sequence.id.clone(),
                    parameter: ParameterKey::Ploidy,
                    context,
                });
            }
        }
    }
    for (parameter, set) in &rule_sets {
        if *parameter != ParameterKey::Ploidy
            && set.has_base_rule()
            && !effective_requests
                .iter()
                .any(|request| request.parameter == *parameter)
        {
            effective_requests.push(ResolutionRequest {
                object_id: "global".to_owned(),
                parameter: *parameter,
                context: ResolutionContext::default(),
            });
        }
    }
    effective_requests.sort_by(|first, second| {
        (&first.object_id, first.parameter.canonical_rank())
            .cmp(&(&second.object_id, second.parameter.canonical_rank()))
    });
    let mut seen_requests = BTreeSet::new();
    let mut values = Vec::with_capacity(effective_requests.len());
    for request in effective_requests {
        if !seen_requests.insert((request.object_id.clone(), request.parameter)) {
            return Err(config_error(
                general,
                "resolution_requests",
                format!(
                    "solicitud repetida para {}:{}",
                    request.object_id,
                    request.parameter.as_str()
                ),
            ));
        }
        validate_resolution_request(general, scenarios, &request, domain)?;
        let set = rule_sets.get(&request.parameter).ok_or_else(|| {
            config_error(
                general,
                "resolution_requests.parameter",
                format!(
                    "no hay reglas declaradas para {}",
                    request.parameter.as_str()
                ),
            )
        })?;
        let winner = set
            .resolve(&request.context)
            .map_err(|error| error.with("path", general.source().path().display().to_string()))?;
        if request.parameter == ParameterKey::Ploidy {
            let RuleValue::Integer(ploidy @ 1..=10) = winner.value() else {
                return Err(Diagnostic::configuration(
                    "GPS-CONF-014",
                    "Una secuencia carece de ploidía efectiva válida.",
                )
                .with("reference", general.reference().id.value())
                .with("sequence", &request.object_id)
                .with("found", format!("{:?}", winner.value()))
                .with("accepted", "integer 1..=10"));
            };
            let sequence = sequences
                .iter_mut()
                .find(|sequence| sequence.id == request.object_id)
                .ok_or_else(|| {
                    config_error(
                        general,
                        "resolution_requests.object_id",
                        format!("secuencia inexistente: {}", request.object_id),
                    )
                })?;
            sequence.ploidy = u8::try_from(*ploidy).map_err(|_| {
                config_error(general, "rules.value", "la ploidía ganadora no cabe en u8")
            })?;
        }
        values.push(EffectiveValue {
            object_id: request.object_id,
            parameter: request.parameter,
            resolution_context: request.context,
            value: winner.value().clone(),
            winning_rule_id: winner.id().to_owned(),
            selectors: winner.selectors().clone(),
            provenance: winner.provenance().clone(),
        });
    }
    values.sort_by(|first, second| {
        (&first.object_id, first.parameter.canonical_rank())
            .cmp(&(&second.object_id, second.parameter.canonical_rank()))
    });
    let reproducibility = ReproducibilityMetadata::new(
        root_seed,
        seed_origin,
        general.reference().id.value(),
        &summary.explicit_event_ids,
    )?;
    let mut reduced_references = Vec::with_capacity(general.reduced_references().len());
    for reduction in general.reduced_references() {
        let (target, target_bases) = match &reduction.target {
            DeclaredReductionTarget::Bases(value) => (
                EffectiveReductionTarget::Bases {
                    declared: *value.value(),
                },
                *value.value(),
            ),
            DeclaredReductionTarget::ReductionPercent(value) => {
                let exact = ExactDecimal::parse("reduction_percent", value.literal()).map_err(
                    |diagnostic| {
                        diagnostic
                            .with("path", general.source().path().display().to_string())
                            .with("location", value.location().to_string())
                    },
                )?;
                (
                    EffectiveReductionTarget::ReductionPercent {
                        literal: value.literal().to_owned(),
                    },
                    exact
                        .retained_after_percent(summary.total_genome_size)
                        .map_err(|diagnostic| {
                            diagnostic
                                .with("path", general.source().path().display().to_string())
                                .with("location", value.location().to_string())
                        })?,
                )
            }
        };
        if target_bases == 0 || target_bases >= summary.total_genome_size {
            return Err(Diagnostic::configuration(
                "GPS-CONF-021",
                "El objetivo de la referencia reducida no es menor que el espacio elegible.",
            )
            .with("path", general.source().path().display().to_string())
            .with("reference", reduction.id.value())
            .with("target_bases", target_bases.to_string())
            .with("eligible_bases", summary.total_genome_size.to_string()));
        }
        let request = SubseedRequest::from_scope(
            "reference_reduction",
            &[
                ("source_reference_id", reduction.source_reference_id.value()),
                ("derived_reference_id", reduction.id.value()),
                ("algorithm", "reduced_reference_selection_v1"),
            ],
        )?;
        reduced_references.push(EffectiveReducedReference {
            id: reduction.id.value().clone(),
            source_reference_id: reduction.source_reference_id.value().clone(),
            region_size_bases: *reduction.region_size_bases.value(),
            target,
            target_bases,
            selection_subseed: derive_subseed(root_seed, &request)?,
        });
    }
    let reference_source = match &general.reference().source {
        DeclaredReferenceSource::Simulated { .. } => EffectiveReferenceSource::Simulated,
        DeclaredReferenceSource::Supplied { fasta } => EffectiveReferenceSource::Supplied {
            original_path: resolve_link(general.source().path(), fasta.value()),
        },
    };

    let mut configuration_inputs = vec![EffectiveConfigurationInput {
        kind: "general",
        id: "general".to_owned(),
        original_path: general.source().path().display().to_string(),
        bytes: Arc::from(general.source().bytes()),
    }];
    configuration_inputs.extend(
        scenarios
            .iter()
            .map(|scenario| EffectiveConfigurationInput {
                kind: "scenario",
                id: scenario.id().value().clone(),
                original_path: scenario.source().path().display().to_string(),
                bytes: Arc::from(scenario.source().bytes()),
            }),
    );
    configuration_inputs
        .sort_by(|first, second| (first.kind, &first.id).cmp(&(second.kind, &second.id)));
    let mut scenario_ids = general
        .scenarios()
        .iter()
        .map(|scenario| scenario.id.value().clone())
        .collect::<Vec<_>>();
    scenario_ids.sort();

    Ok(EffectiveConfig {
        schema_version: "1",
        source_path: general.source().path().display().to_string(),
        route: *general.route().value(),
        reference_id: general.reference().id.value().clone(),
        reference_source,
        reduced_references,
        total_genome_size: summary.total_genome_size,
        sequences,
        scenario_ids,
        configuration_inputs,
        explicit_event_ids: summary.explicit_event_ids,
        individual_ids: summary.individual_ids,
        normalized,
        controlled,
        values,
        resolution_domain: domain.clone(),
        resolution_requests,
        reproducibility,
    })
}

fn validate_resolution_domain(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
    domain: &ResolutionDomain,
) -> Result<(), Diagnostic> {
    let contains = |entities: &[crate::configuration::DeclaredEntity], id: &str| {
        entities.iter().any(|entity| entity.id().value() == id)
    };
    let invalid = |detail: String| {
        config_error(general, "resolution_domain", detail).with("provider", domain.provider())
    };
    for (subpopulation, population) in domain.population_memberships() {
        if !contains(general.populations(), population)
            || !contains(general.subpopulations(), subpopulation)
        {
            return Err(invalid(format!(
                "relación con identidad no declarada: {population}/{subpopulation}"
            )));
        }
    }
    for phase in domain.phase_generation_intervals().keys() {
        let declared = contains(general.phases(), phase)
            || scenarios
                .iter()
                .any(|scenario| contains(scenario.phases(), phase));
        if !declared {
            return Err(invalid(format!("intervalo de fase no declarada: {phase}")));
        }
    }
    for locus in domain.loci().values() {
        let sequence = general
            .reference()
            .sequences
            .iter()
            .find(|sequence| sequence.id.value() == &locus.sequence);
        if locus.reference != *general.reference().id.value()
            || sequence.is_none()
            || locus.position > *sequence.expect("checked sequence").length.value()
        {
            return Err(invalid(format!(
                "locus fuera de la referencia declarada: {}",
                locus.id
            )));
        }
    }
    Ok(())
}

fn normalize_reference_composition(
    general: &DeclaredGeneralConfig,
    composition: &[crate::configuration::DeclaredFrequencyValue],
) -> Result<NormalizedCollection, Diagnostic> {
    let mut normalized = normalize_declared(
        "reference.composition",
        composition,
        general.normalization_tolerance(),
    )?;
    let mut present = BTreeSet::new();
    for value in &normalized.values {
        if !matches!(value.id.as_str(), "A" | "C" | "G" | "T")
            || value.category != "base"
            || !present.insert(value.id.clone())
        {
            return Err(config_error(
                general,
                "reference.composition",
                format!(
                    "cada elemento debe usar una vez un ID exacto A/C/G/T y category=base; encontrado {}/{}",
                    value.id, value.category
                ),
            ));
        }
    }
    for id in ["A", "C", "G", "T"] {
        if !present.contains(id) {
            normalized.values.push(crate::numeric::NormalizedValue {
                id: id.to_owned(),
                category: "base".to_owned(),
                before: 0.0,
                after: 0.0,
            });
        }
    }
    normalized
        .values
        .sort_by(|first, second| (&first.id, &first.category).cmp(&(&second.id, &second.category)));
    Ok(normalized)
}

fn validate_declared_rules(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
) -> Result<(), Diagnostic> {
    let mut identities = SelectorIdentityCatalog::default();
    identities.add_reference(general.reference().id.value());
    for sequence in &general.reference().sequences {
        identities.add_sequence(general.reference().id.value(), sequence.id.value());
    }
    for entity in general.individuals() {
        identities.add_individual(entity.id().value());
    }
    for entity in general.populations() {
        identities.add_population(entity.id().value());
    }
    for entity in general.subpopulations() {
        identities.add_subpopulation(entity.id().value());
    }
    for entity in general.phases() {
        identities.add_phase(entity.id().value());
    }
    for scenario in scenarios {
        for entity in scenario.phases() {
            identities.add_phase(entity.id().value());
        }
    }
    for entity in general.groups() {
        identities.add_group(entity.id().value());
    }
    for entity in general.blocks() {
        identities.add_block(entity.id().value());
    }
    for entity in general.events() {
        identities.add_technical_event(entity.id().value());
    }
    validate_rules(general.rules(), &identities)?;
    for (parameter, rules) in group_rules(general.rules()) {
        let set = RuleSet::new_with_fallback(parameter, rules, parameter == ParameterKey::Ploidy)
            .map_err(|error| {
            error.with("path", general.source().path().display().to_string())
        })?;
        set.validate_declared_intersections()
            .map_err(|error| error.with("path", general.source().path().display().to_string()))?;
    }
    Ok(())
}

fn group_rules(rules: &[Rule]) -> BTreeMap<ParameterKey, Vec<Rule>> {
    let mut grouped = BTreeMap::new();
    for rule in rules {
        grouped
            .entry(rule.parameter())
            .or_insert_with(Vec::new)
            .push(rule.clone());
    }
    grouped
}

pub(crate) fn validate_resolution_request(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
    request: &ResolutionRequest,
    domain: &ResolutionDomain,
) -> Result<(), Diagnostic> {
    let context = &request.context;
    let contains = |entities: &[crate::configuration::DeclaredEntity], id: &str| {
        entities.iter().any(|entity| entity.id().value() == id)
    };
    let invalid = |detail: String| {
        config_error(general, "resolution_requests", detail)
            .with("object_id", &request.object_id)
            .with("parameter", request.parameter.as_str())
    };
    let present_axes = [
        (context.individual.is_some(), AxisKind::Individual),
        (
            context.population.is_some() || context.subpopulation.is_some(),
            AxisKind::Population,
        ),
        (
            context.phase.is_some() || context.generation.is_some(),
            AxisKind::Temporal,
        ),
        (
            context.reference.is_some()
                || context.sequence.is_some()
                || context.position.is_some()
                || context.locus.is_some(),
            AxisKind::Spatial,
        ),
        (context.group.is_some(), AxisKind::Group),
        (context.block.is_some(), AxisKind::Block),
        (context.technical_event.is_some(), AxisKind::TechnicalEvent),
    ];
    for (present, axis) in present_axes {
        if present && !metadata(request.parameter).allowed_axes.contains(&axis) {
            return Err(Diagnostic::configuration(
                "GPS-CONF-012",
                "La familia del parámetro no admite uno de los selectores declarados.",
            )
            .with("parameter", request.parameter.as_str())
            .with("selector", axis.as_str())
            .with("scope", "resolution_request")
            .with("object_id", &request.object_id)
            .with("path", general.source().path().display().to_string()));
        }
    }
    let spatial_level = if context.locus.is_some() {
        Some(SpatialLevel::Locus)
    } else if context.position.is_some() {
        Some(SpatialLevel::Region)
    } else if context.sequence.is_some() || context.reference.is_some() {
        Some(SpatialLevel::Sequence)
    } else {
        None
    };
    if let Some(level) = spatial_level
        && !metadata(request.parameter)
            .allowed_spatial_levels
            .contains(&level)
    {
        return Err(Diagnostic::configuration(
            "GPS-CONF-012",
            "La familia del parámetro no admite el nivel espacial solicitado.",
        )
        .with("parameter", request.parameter.as_str())
        .with("spatial_level", level.as_str())
        .with("scope", "resolution_request")
        .with("object_id", &request.object_id)
        .with("path", general.source().path().display().to_string()));
    }
    if context.subpopulation.is_some() && context.population.is_none() {
        return Err(invalid(
            "una subpoblación exige su población en el contexto".to_owned(),
        ));
    }
    if context.generation.is_some() && context.phase.is_none() {
        return Err(invalid(
            "una generación exige su fase en el contexto".to_owned(),
        ));
    }
    if context.reference.is_some() != context.sequence.is_some() {
        return Err(invalid(
            "una solicitud separa referencia y secuencia".to_owned(),
        ));
    }
    if (context.position.is_some() || context.locus.is_some())
        && (context.reference.is_none() || context.sequence.is_none())
    {
        return Err(invalid(
            "posición o locus exigen referencia y secuencia".to_owned(),
        ));
    }
    if context.locus.is_some() && context.position.is_none() {
        return Err(invalid(
            "un locus evaluado exige su posición en el contexto".to_owned(),
        ));
    }
    if let (Some(population), Some(subpopulation)) = (&context.population, &context.subpopulation)
        && domain.population_memberships().get(subpopulation) != Some(population)
    {
        return Err(invalid(
            "la relación población-subpoblación carece de prueba autoritativa".to_owned(),
        ));
    }
    if let (Some(phase), Some(generation)) = (&context.phase, context.generation) {
        let proven = domain
            .phase_generation_intervals()
            .get(phase)
            .is_some_and(|(start, end)| (*start..=*end).contains(&generation));
        if !proven {
            return Err(invalid(
                "la relación fase-generación carece de prueba autoritativa".to_owned(),
            ));
        }
    }
    if let (Some(locus), Some(reference), Some(sequence), Some(position)) = (
        &context.locus,
        &context.reference,
        &context.sequence,
        context.position,
    ) {
        let proven = domain.loci().get(locus).is_some_and(|declared| {
            declared.reference == *reference
                && declared.sequence == *sequence
                && declared.position == position
        });
        if !proven {
            return Err(invalid(
                "el locus y su coordenada carecen de prueba autoritativa".to_owned(),
            ));
        }
    }
    let mut candidates = BTreeSet::new();
    if let Some(id) = &context.individual {
        if !contains(general.individuals(), id) {
            return Err(invalid(format!("individuo inexistente: {id}")));
        }
        candidates.insert(id.as_str());
    }
    if let Some(id) = &context.population {
        if !contains(general.populations(), id) {
            return Err(invalid(format!("población inexistente: {id}")));
        }
        candidates.insert(id.as_str());
    }
    if let Some(id) = &context.subpopulation {
        if !contains(general.subpopulations(), id) {
            return Err(invalid(format!("subpoblación inexistente: {id}")));
        }
        candidates.insert(id.as_str());
    }
    if let Some(id) = &context.phase {
        let declared = contains(general.phases(), id)
            || scenarios
                .iter()
                .any(|scenario| contains(scenario.phases(), id));
        if !declared {
            return Err(invalid(format!("fase inexistente: {id}")));
        }
        candidates.insert(id.as_str());
    }
    if let (Some(reference), Some(sequence)) = (&context.reference, &context.sequence) {
        let declared_sequence = general
            .reference()
            .sequences
            .iter()
            .find(|declared| declared.id.value() == sequence);
        if reference != general.reference().id.value() || declared_sequence.is_none() {
            return Err(invalid(format!(
                "secuencia inexistente: {reference}:{sequence}"
            )));
        }
        if let Some(position) = context.position
            && (position == 0
                || position
                    > *declared_sequence
                        .expect("checked declared sequence")
                        .length
                        .value())
        {
            return Err(invalid(format!("posición fuera de secuencia: {position}")));
        }
        if let Some(locus) = &context.locus {
            candidates.insert(locus.as_str());
        } else {
            candidates.insert(sequence.as_str());
        }
    }
    for (id, entities, kind) in [
        (context.group.as_ref(), general.groups(), "grupo"),
        (context.block.as_ref(), general.blocks(), "bloque"),
        (
            context.technical_event.as_ref(),
            general.events(),
            "evento técnico",
        ),
    ] {
        if let Some(id) = id {
            if !contains(entities, id) {
                return Err(invalid(format!("{kind} inexistente: {id}")));
            }
            candidates.insert(id.as_str());
        }
    }
    if *context == ResolutionContext::default() {
        if request.object_id != "global" || request.parameter == ParameterKey::Ploidy {
            return Err(invalid(
                "object_id no identifica correctamente el ámbito global".to_owned(),
            ));
        }
    } else if !candidates.contains(request.object_id.as_str()) {
        return Err(invalid(
            "object_id no identifica ningún objeto declarado de su contexto".to_owned(),
        ));
    }
    if request.parameter == ParameterKey::Ploidy
        && context.sequence.as_deref() != Some(request.object_id.as_str())
    {
        return Err(invalid(
            "object_id no identifica la secuencia cuya ploidía se resuelve".to_owned(),
        ));
    }
    Ok(())
}

fn validate_route_and_scenarios(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
) -> Result<(), Diagnostic> {
    let route = *general.route().value();
    match route {
        Route::Direct if general.scenarios().is_empty() && scenarios.is_empty() => return Ok(()),
        Route::Direct => {
            return Err(config_error(
                general,
                "route,scenarios",
                "la vía directa no admite escenarios poblacionales",
            ));
        }
        Route::Evolution | Route::Families if general.scenarios().is_empty() => {
            return Err(config_error(
                general,
                "route,scenarios",
                "la vía declarada exige al menos un escenario",
            ));
        }
        Route::Evolution | Route::Families => {}
    }

    let by_id: BTreeMap<_, _> = scenarios
        .iter()
        .map(|scenario| (scenario.id().value().as_str(), *scenario))
        .collect();
    if by_id.len() != scenarios.len() || scenarios.len() != general.scenarios().len() {
        return Err(config_error(
            general,
            "scenarios.id",
            "los escenarios cargados no coinciden uno a uno con los declarados",
        ));
    }
    let mut declared_ids = BTreeSet::new();
    for link in general.scenarios() {
        if !declared_ids.insert(link.id.value()) {
            return Err(config_error(
                general,
                "scenarios.id",
                "identificador de escenario repetido",
            ));
        }
        let Some(scenario) = by_id.get(link.id.value().as_str()) else {
            return Err(config_error(
                general,
                "scenarios.id",
                format!("no se cargó el escenario {}", link.id.value()),
            ));
        };
        if *link.route.value() != route || *scenario.route().value() != route {
            return Err(config_error(
                general,
                "route,scenarios.route",
                format!(
                    "la vía de {} no coincide con {}",
                    link.id.value(),
                    route.as_str()
                ),
            ));
        }
    }
    Ok(())
}

fn config_error(
    general: &DeclaredGeneralConfig,
    fields: impl Into<String>,
    detail: impl Into<String>,
) -> Diagnostic {
    config_error_at(general.source().path(), fields, detail)
}

fn config_error_at(
    path: &Path,
    fields: impl Into<String>,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", path.display().to_string())
    .with("fields", fields.into())
    .with("detail", detail.into())
}
