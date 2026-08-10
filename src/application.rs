//! Orquestación de casos de uso mediante contratos del dominio.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::configuration::{
    DeclaredConfig, DeclaredGeneralConfig, DeclaredReferenceSource, DeclaredScenarioConfig,
    ExpectedConfigKind, FrequencyMode, Route, SequenceClass,
};
use crate::diagnostics::Diagnostic;
use crate::identity::{
    BlockId, CrossId, EventId, GroupId, IdentityRegistry, IndividualId, PhaseId, PopulationId,
    ReferenceId, ScenarioId, SequenceId, SubpopulationId,
};
use crate::numeric::{NormalizedCollection, normalize_declared};
use crate::randomness::{EntropySource, ReproducibilityMetadata, SystemEntropy, resolve_root_seed};
use crate::rules::{
    ParameterKey, ResolutionContext, Rule, RuleProvenance, RuleSet, RuleValue,
    SelectorIdentityCatalog, SelectorSet, validate_rules,
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

/// Valor efectivo con la identidad y procedencia completas de la regla ganadora.
#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveValue {
    object_id: String,
    parameter: ParameterKey,
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
    total_genome_size: u64,
    sequences: Vec<EffectiveSequence>,
    scenario_ids: Vec<String>,
    configuration_inputs: Vec<EffectiveConfigurationInput>,
    explicit_event_ids: Vec<String>,
    individual_ids: Vec<String>,
    normalized: Vec<NormalizedCollection>,
    controlled: Vec<EffectiveControlledSet>,
    values: Vec<EffectiveValue>,
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
    for sequence in &general.reference().sequences {
        let _ = SequenceId::parse(sequence.id.value(), general.reference().id.value())?;
    }
    for scenario in general.scenarios() {
        let _ = ScenarioId::parse(scenario.id.value(), "execution")?;
    }
    for scenario in scenarios {
        let scenario_id = ScenarioId::parse(scenario.id().value(), "execution")?;
        for phase in scenario.phases() {
            let _ = PhaseId::parse(phase.id().value(), scenario_id.as_str())?;
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
    for phase in general.phases() {
        let _ = PhaseId::parse(phase.id().value(), "execution")?;
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
    build_effective_config_with_entropy(general, scenarios, requests, &mut SystemEntropy)
}

pub fn build_effective_config_with_entropy(
    general: &DeclaredGeneralConfig,
    scenarios: &[&DeclaredScenarioConfig],
    requests: &[ResolutionRequest],
    entropy: &mut impl EntropySource,
) -> Result<EffectiveConfig, Diagnostic> {
    let summary = validate_general_config(general, scenarios)?;
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
        normalized.push(normalize_declared(
            "reference.composition",
            composition,
            general.normalization_tolerance(),
        )?);
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
    let mut effective_requests = requests.to_vec();
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
            value: winner.value().clone(),
            winning_rule_id: winner.id().to_owned(),
            selectors: winner.selectors().clone(),
            provenance: winner.provenance().clone(),
        });
    }
    values.sort_by(|first, second| {
        (&first.object_id, first.parameter).cmp(&(&second.object_id, second.parameter))
    });
    let reproducibility = ReproducibilityMetadata::new(
        root_seed,
        seed_origin,
        general.reference().id.value(),
        &summary.explicit_event_ids,
    )?;

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
        total_genome_size: summary.total_genome_size,
        sequences,
        scenario_ids,
        configuration_inputs,
        explicit_event_ids: summary.explicit_event_ids,
        individual_ids: summary.individual_ids,
        normalized,
        controlled,
        values,
        reproducibility,
    })
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
