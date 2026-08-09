//! Catálogo normativo y tipos del resolvedor de reglas de configuración.

use std::collections::{BTreeMap, BTreeSet};

use crate::diagnostics::Diagnostic;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum AxisKind {
    Individual,
    Population,
    Temporal,
    Spatial,
    Group,
    Block,
    TechnicalEvent,
}

impl AxisKind {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Individual => "individual",
            Self::Population => "population",
            Self::Temporal => "temporal",
            Self::Spatial => "spatial",
            Self::Group => "group",
            Self::Block => "block",
            Self::TechnicalEvent => "technical_event",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum ParameterKey {
    Ploidy,
    VariantRate,
    SelfingRate,
    RecombinationRate,
    Coverage,
    ErrorRate,
    ContaminationRate,
    OutputEnabled,
}

impl ParameterKey {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Ploidy => "ploidy",
            Self::VariantRate => "variant_rate",
            Self::SelfingRate => "selfing_rate",
            Self::RecombinationRate => "recombination_rate",
            Self::Coverage => "coverage",
            Self::ErrorRate => "error_rate",
            Self::ContaminationRate => "contamination_rate",
            Self::OutputEnabled => "output_enabled",
        }
    }

    #[must_use]
    pub fn parse(value: &str) -> Option<Self> {
        match value {
            "ploidy" => Some(Self::Ploidy),
            "variant_rate" => Some(Self::VariantRate),
            "selfing_rate" => Some(Self::SelfingRate),
            "recombination_rate" => Some(Self::RecombinationRate),
            "coverage" => Some(Self::Coverage),
            "error_rate" => Some(Self::ErrorRate),
            "contamination_rate" => Some(Self::ContaminationRate),
            "output_enabled" => Some(Self::OutputEnabled),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum ParameterFamily {
    StructureAndPloidy,
    CreationFrequencyAndMutation,
    DemographyReproductionAndSelection,
    Recombination,
    WgsCoverageProfileAndError,
    ContaminationAndDemultiplexing,
    OutputsAndObservation,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SpatialLevel {
    Sequence,
    Region,
    Locus,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParameterValueDomain {
    Ploidy,
    Probability,
    NonNegativeReal,
    Boolean,
}

impl ParameterValueDomain {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Ploidy => "integer_1_10",
            Self::Probability => "probability_f64_0_1",
            Self::NonNegativeReal => "finite_non_negative_real",
            Self::Boolean => "boolean",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ParameterMetadata {
    pub parameter: ParameterKey,
    pub family: ParameterFamily,
    pub is_biological: bool,
    pub allowed_axes: &'static [AxisKind],
    pub allowed_spatial_levels: &'static [SpatialLevel],
    pub value_domain: ParameterValueDomain,
}

const INDIVIDUAL: AxisKind = AxisKind::Individual;
const POPULATION: AxisKind = AxisKind::Population;
const TEMPORAL: AxisKind = AxisKind::Temporal;
const SPATIAL: AxisKind = AxisKind::Spatial;
const GROUP: AxisKind = AxisKind::Group;
const BLOCK: AxisKind = AxisKind::Block;
const TECHNICAL_EVENT: AxisKind = AxisKind::TechnicalEvent;
const SEQUENCE: SpatialLevel = SpatialLevel::Sequence;
const REGION: SpatialLevel = SpatialLevel::Region;
const LOCUS: SpatialLevel = SpatialLevel::Locus;

pub const PARAMETER_CATALOG: [ParameterMetadata; 8] = [
    ParameterMetadata {
        parameter: ParameterKey::Ploidy,
        family: ParameterFamily::StructureAndPloidy,
        is_biological: true,
        allowed_axes: &[SPATIAL],
        allowed_spatial_levels: &[SEQUENCE],
        value_domain: ParameterValueDomain::Ploidy,
    },
    ParameterMetadata {
        parameter: ParameterKey::VariantRate,
        family: ParameterFamily::CreationFrequencyAndMutation,
        is_biological: true,
        allowed_axes: &[INDIVIDUAL, POPULATION, TEMPORAL, SPATIAL, GROUP, BLOCK],
        allowed_spatial_levels: &[SEQUENCE, REGION, LOCUS],
        value_domain: ParameterValueDomain::NonNegativeReal,
    },
    ParameterMetadata {
        parameter: ParameterKey::SelfingRate,
        family: ParameterFamily::DemographyReproductionAndSelection,
        is_biological: true,
        allowed_axes: &[POPULATION, TEMPORAL],
        allowed_spatial_levels: &[],
        value_domain: ParameterValueDomain::Probability,
    },
    ParameterMetadata {
        parameter: ParameterKey::RecombinationRate,
        family: ParameterFamily::Recombination,
        is_biological: true,
        allowed_axes: &[SPATIAL, GROUP],
        allowed_spatial_levels: &[SEQUENCE],
        value_domain: ParameterValueDomain::NonNegativeReal,
    },
    ParameterMetadata {
        parameter: ParameterKey::Coverage,
        family: ParameterFamily::WgsCoverageProfileAndError,
        is_biological: false,
        allowed_axes: &[TECHNICAL_EVENT, INDIVIDUAL, SPATIAL],
        allowed_spatial_levels: &[SEQUENCE, REGION],
        value_domain: ParameterValueDomain::NonNegativeReal,
    },
    ParameterMetadata {
        parameter: ParameterKey::ErrorRate,
        family: ParameterFamily::WgsCoverageProfileAndError,
        is_biological: false,
        allowed_axes: &[TECHNICAL_EVENT, INDIVIDUAL, SPATIAL],
        allowed_spatial_levels: &[SEQUENCE, REGION],
        value_domain: ParameterValueDomain::Probability,
    },
    ParameterMetadata {
        parameter: ParameterKey::ContaminationRate,
        family: ParameterFamily::ContaminationAndDemultiplexing,
        is_biological: false,
        allowed_axes: &[TECHNICAL_EVENT, INDIVIDUAL],
        allowed_spatial_levels: &[],
        value_domain: ParameterValueDomain::Probability,
    },
    ParameterMetadata {
        parameter: ParameterKey::OutputEnabled,
        family: ParameterFamily::OutputsAndObservation,
        is_biological: false,
        allowed_axes: &[TECHNICAL_EVENT, INDIVIDUAL],
        allowed_spatial_levels: &[],
        value_domain: ParameterValueDomain::Boolean,
    },
];

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum PopulationSelector {
    Population(String),
    Subpopulation(String),
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum TemporalSelector {
    Phase(String),
    Generation(u32),
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct RegionSelector {
    pub reference: String,
    pub sequence: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum SpatialSelector {
    Sequence {
        reference: String,
        sequence: String,
    },
    Region(RegionSelector),
    Locus {
        reference: String,
        sequence: String,
        locus: String,
    },
}

impl SpatialSelector {
    #[must_use]
    pub const fn level(&self) -> SpatialLevel {
        match self {
            Self::Sequence { .. } => SpatialLevel::Sequence,
            Self::Region(_) => SpatialLevel::Region,
            Self::Locus { .. } => SpatialLevel::Locus,
        }
    }
}

#[derive(Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct SelectorSet {
    pub individual: Option<String>,
    pub population: Option<PopulationSelector>,
    pub temporal: Option<TemporalSelector>,
    pub spatial: Option<SpatialSelector>,
    pub group: Option<String>,
    pub block: Option<String>,
    pub technical_event: Option<String>,
}

impl SelectorSet {
    #[must_use]
    pub fn axes(&self) -> Vec<AxisKind> {
        let mut result = Vec::with_capacity(7);
        if self.individual.is_some() {
            result.push(AxisKind::Individual);
        }
        if self.population.is_some() {
            result.push(AxisKind::Population);
        }
        if self.temporal.is_some() {
            result.push(AxisKind::Temporal);
        }
        if self.spatial.is_some() {
            result.push(AxisKind::Spatial);
        }
        if self.group.is_some() {
            result.push(AxisKind::Group);
        }
        if self.block.is_some() {
            result.push(AxisKind::Block);
        }
        if self.technical_event.is_some() {
            result.push(AxisKind::TechnicalEvent);
        }
        result
    }
}

#[must_use]
pub fn metadata(parameter: ParameterKey) -> &'static ParameterMetadata {
    match parameter {
        ParameterKey::Ploidy => &PARAMETER_CATALOG[0],
        ParameterKey::VariantRate => &PARAMETER_CATALOG[1],
        ParameterKey::SelfingRate => &PARAMETER_CATALOG[2],
        ParameterKey::RecombinationRate => &PARAMETER_CATALOG[3],
        ParameterKey::Coverage => &PARAMETER_CATALOG[4],
        ParameterKey::ErrorRate => &PARAMETER_CATALOG[5],
        ParameterKey::ContaminationRate => &PARAMETER_CATALOG[6],
        ParameterKey::OutputEnabled => &PARAMETER_CATALOG[7],
    }
}

pub fn validate_selector_set(
    parameter: ParameterKey,
    selectors: &SelectorSet,
) -> Result<(), Diagnostic> {
    let metadata = metadata(parameter);
    for axis in selectors.axes() {
        if !metadata.allowed_axes.contains(&axis) {
            return Err(selector_error(parameter, axis.as_str()));
        }
    }
    if let Some(spatial) = &selectors.spatial
        && !metadata.allowed_spatial_levels.contains(&spatial.level())
    {
        return Err(selector_error(parameter, "spatial_level"));
    }
    if metadata.is_biological && selectors.technical_event.is_some() {
        return Err(selector_error(parameter, "technical_event"));
    }
    Ok(())
}

fn selector_error(parameter: ParameterKey, selector: &str) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-012",
        "La familia del parámetro no admite uno de los selectores declarados.",
    )
    .with("parameter", parameter.as_str())
    .with("selector", selector)
    .with("scope", "declared_rule")
}

#[derive(Clone, Debug, PartialEq)]
pub enum RuleValue {
    Float(f64),
    Integer(i64),
    Boolean(bool),
    Text(String),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RuleProvenance {
    pub path: String,
    pub location: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Rule {
    id: String,
    parameter: ParameterKey,
    value: RuleValue,
    selectors: SelectorSet,
    provenance: RuleProvenance,
}

impl Rule {
    pub fn new(
        id: impl Into<String>,
        parameter: ParameterKey,
        value: RuleValue,
        selectors: SelectorSet,
        provenance: RuleProvenance,
    ) -> Result<Self, Diagnostic> {
        let id = id.into();
        validate_selector_set(parameter, &selectors)?;
        validate_rule_value(&id, parameter, &value, &selectors, &provenance)?;
        Ok(Self {
            id,
            parameter,
            value,
            selectors,
            provenance,
        })
    }

    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
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
    pub const fn selectors(&self) -> &SelectorSet {
        &self.selectors
    }

    #[must_use]
    pub const fn provenance(&self) -> &RuleProvenance {
        &self.provenance
    }

    #[must_use]
    pub fn matches(&self, context: &ResolutionContext) -> bool {
        selector_matches(
            self.selectors.individual.as_deref(),
            context.individual.as_deref(),
        ) && population_matches(self.selectors.population.as_ref(), context)
            && temporal_matches(self.selectors.temporal.as_ref(), context)
            && spatial_matches(self.selectors.spatial.as_ref(), context)
            && selector_matches(self.selectors.group.as_deref(), context.group.as_deref())
            && selector_matches(self.selectors.block.as_deref(), context.block.as_deref())
            && selector_matches(
                self.selectors.technical_event.as_deref(),
                context.technical_event.as_deref(),
            )
    }
}

fn validate_rule_value(
    id: &str,
    parameter: ParameterKey,
    value: &RuleValue,
    selectors: &SelectorSet,
    provenance: &RuleProvenance,
) -> Result<(), Diagnostic> {
    let domain = metadata(parameter).value_domain;
    let is_valid = match (domain, value) {
        (ParameterValueDomain::Ploidy, RuleValue::Integer(value)) => (1..=10).contains(value),
        (ParameterValueDomain::Probability, RuleValue::Float(value)) => {
            value.is_finite() && (0.0..=1.0).contains(value)
        }
        (ParameterValueDomain::NonNegativeReal, RuleValue::Float(value)) => {
            value.is_finite() && *value >= 0.0
        }
        (ParameterValueDomain::Boolean, RuleValue::Boolean(_)) => true,
        _ => false,
    };
    if is_valid {
        return Ok(());
    }
    let found = render_rule_value(value);
    let accepted = match domain {
        ParameterValueDomain::Ploidy => "integer 1..=10",
        ParameterValueDomain::Probability => "finite f64 [0,1]",
        ParameterValueDomain::NonNegativeReal => "finite f64 >= 0",
        ParameterValueDomain::Boolean => "boolean",
    };
    let diagnostic = match (domain, value) {
        (ParameterValueDomain::Probability, RuleValue::Float(_)) => Diagnostic::configuration(
            "GPS-CONF-015",
            "Una probabilidad, frecuencia o tolerancia no es finita o está fuera de su intervalo.",
        )
        .with("parameter", parameter.as_str())
        .with("identity", id)
        .with("found", found)
        .with("accepted", accepted),
        (ParameterValueDomain::Ploidy, _) => {
            let (reference, sequence) = ploidy_scope(selectors);
            Diagnostic::configuration(
                "GPS-CONF-014",
                "Una secuencia carece de ploidía efectiva válida.",
            )
            .with("reference", reference)
            .with("sequence", sequence)
            .with("found", found)
            .with("accepted", accepted)
        }
        (
            ParameterValueDomain::Probability
            | ParameterValueDomain::NonNegativeReal
            | ParameterValueDomain::Boolean,
            _,
        ) => Diagnostic::configuration(
            "GPS-CONF-008",
            "La configuración contiene una combinación de campos incompatible.",
        )
        .with("fields", "rules.value")
        .with(
            "detail",
            format!(
                "{} exige {accepted}; se encontró {found}",
                parameter.as_str()
            ),
        )
        .with("parameter", parameter.as_str())
        .with("found", found)
        .with("accepted", accepted),
    };
    Err(diagnostic
        .with("path", &provenance.path)
        .with("location", &provenance.location))
}

fn ploidy_scope(selectors: &SelectorSet) -> (&str, &str) {
    match selectors.spatial.as_ref() {
        Some(SpatialSelector::Sequence {
            reference,
            sequence,
        }) => (reference, sequence),
        _ => ("<general>", "<all>"),
    }
}

fn render_rule_value(value: &RuleValue) -> String {
    match value {
        RuleValue::Float(value) => value.to_string(),
        RuleValue::Integer(value) => value.to_string(),
        RuleValue::Boolean(value) => value.to_string(),
        RuleValue::Text(value) => format!("string:{value}"),
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ResolutionContext {
    pub individual: Option<String>,
    pub population: Option<String>,
    pub subpopulation: Option<String>,
    pub phase: Option<String>,
    pub generation: Option<u32>,
    pub reference: Option<String>,
    pub sequence: Option<String>,
    pub position: Option<u64>,
    pub locus: Option<String>,
    pub group: Option<String>,
    pub block: Option<String>,
    pub technical_event: Option<String>,
}

#[derive(Clone, Debug)]
pub struct RuleSet {
    parameter: ParameterKey,
    rules: Vec<Rule>,
}

impl RuleSet {
    pub fn new(parameter: ParameterKey, rules: Vec<Rule>) -> Result<Self, Diagnostic> {
        Self::new_with_fallback(parameter, rules, false)
    }

    pub fn new_with_fallback(
        parameter: ParameterKey,
        rules: Vec<Rule>,
        fallback_available: bool,
    ) -> Result<Self, Diagnostic> {
        if rules.iter().any(|rule| rule.parameter != parameter) {
            return Err(Diagnostic::configuration(
                "GPS-CONF-008",
                "La configuración contiene una combinación de campos incompatible.",
            )
            .with("fields", "rules.parameter")
            .with("detail", "un conjunto contiene reglas de otro parámetro"));
        }
        if !fallback_available
            && !rules
                .iter()
                .any(|rule| rule.selectors == SelectorSet::default())
        {
            return Err(Diagnostic::configuration(
                "GPS-CONF-008",
                "La configuración contiene una combinación de campos incompatible.",
            )
            .with("fields", parameter.as_str())
            .with("detail", "falta la regla base obligatoria"));
        }
        Ok(Self { parameter, rules })
    }

    #[must_use]
    pub fn has_base_rule(&self) -> bool {
        self.rules
            .iter()
            .any(|rule| rule.selectors == SelectorSet::default())
    }

    #[must_use]
    pub const fn parameter(&self) -> ParameterKey {
        self.parameter
    }

    #[must_use]
    pub fn rules(&self) -> &[Rule] {
        &self.rules
    }

    #[must_use]
    pub fn matching_rules<'a>(&'a self, context: &ResolutionContext) -> Vec<&'a Rule> {
        self.rules
            .iter()
            .filter(|rule| rule.matches(context))
            .collect()
    }

    pub fn resolve<'a>(&'a self, context: &ResolutionContext) -> Result<&'a Rule, Diagnostic> {
        let matching = self.matching_rules(context);
        let maximal = matching
            .iter()
            .copied()
            .filter(|candidate| {
                !matching
                    .iter()
                    .any(|other| dominates(&other.selectors, &candidate.selectors))
            })
            .collect::<Vec<_>>();
        match maximal.as_slice() {
            [winner] => Ok(*winner),
            [] => Err(Diagnostic::configuration(
                "GPS-CONF-008",
                "La configuración contiene una combinación de campos incompatible.",
            )
            .with("fields", self.parameter.as_str())
            .with("detail", "ninguna regla resulta aplicable")),
            conflicts => Err(Diagnostic::configuration(
                "GPS-CONF-009",
                "Varias reglas máximas incomparables impiden obtener un valor efectivo único.",
            )
            .with("parameter", self.parameter.as_str())
            .with("scope", format!("{context:?}"))
            .with("rules", sorted_rule_ids(conflicts))),
        }
    }

    /// Comprueba todos los cruces compatibles entre ámbitos declarados. Cada
    /// ambigüedad real contiene al menos un par de máximas incomparables; la
    /// intersección mínima de ese par permite resolver contra el conjunto
    /// completo y respeta una posible regla combinada dominante.
    pub fn validate_declared_intersections(&self) -> Result<(), Diagnostic> {
        let mut rules = self.rules.iter().collect::<Vec<_>>();
        rules.sort_by(|first, second| first.id.cmp(&second.id));
        for (index, first) in rules.iter().enumerate() {
            for second in &rules[index + 1..] {
                if dominates(&first.selectors, &second.selectors)
                    || dominates(&second.selectors, &first.selectors)
                {
                    continue;
                }
                let Some(context) = intersection_context(&first.selectors, &second.selectors)
                else {
                    continue;
                };
                if let Err(error) = self.resolve(&context)
                    && error.code() == "GPS-CONF-009"
                {
                    return Err(error);
                }
            }
        }
        Ok(())
    }
}

fn intersection_context(first: &SelectorSet, second: &SelectorSet) -> Option<ResolutionContext> {
    let individual = merge_exact(&first.individual, &second.individual)?;
    let (population, subpopulation) =
        merge_population(first.population.as_ref(), second.population.as_ref())?;
    let (phase, generation) = merge_temporal(first.temporal.as_ref(), second.temporal.as_ref())?;
    let (reference, sequence, position, locus) =
        merge_spatial(first.spatial.as_ref(), second.spatial.as_ref())?;
    Some(ResolutionContext {
        individual,
        population,
        subpopulation,
        phase,
        generation,
        reference,
        sequence,
        position,
        locus,
        group: merge_exact(&first.group, &second.group)?,
        block: merge_exact(&first.block, &second.block)?,
        technical_event: merge_exact(&first.technical_event, &second.technical_event)?,
    })
}

fn merge_exact(first: &Option<String>, second: &Option<String>) -> Option<Option<String>> {
    match (first, second) {
        (Some(first), Some(second)) if first != second => None,
        (Some(value), _) | (_, Some(value)) => Some(Some(value.clone())),
        (None, None) => Some(None),
    }
}

fn merge_population(
    first: Option<&PopulationSelector>,
    second: Option<&PopulationSelector>,
) -> Option<(Option<String>, Option<String>)> {
    match (first, second) {
        (None, None) => Some((None, None)),
        (Some(PopulationSelector::Population(id)), None)
        | (None, Some(PopulationSelector::Population(id))) => Some((Some(id.clone()), None)),
        (Some(PopulationSelector::Subpopulation(id)), None)
        | (None, Some(PopulationSelector::Subpopulation(id))) => Some((None, Some(id.clone()))),
        (
            Some(PopulationSelector::Population(population)),
            Some(PopulationSelector::Subpopulation(subpopulation)),
        )
        | (
            Some(PopulationSelector::Subpopulation(subpopulation)),
            Some(PopulationSelector::Population(population)),
        ) => Some((Some(population.clone()), Some(subpopulation.clone()))),
        (
            Some(PopulationSelector::Population(first)),
            Some(PopulationSelector::Population(second)),
        ) if first == second => Some((Some(first.clone()), None)),
        (
            Some(PopulationSelector::Subpopulation(first)),
            Some(PopulationSelector::Subpopulation(second)),
        ) if first == second => Some((None, Some(first.clone()))),
        _ => None,
    }
}

fn merge_temporal(
    first: Option<&TemporalSelector>,
    second: Option<&TemporalSelector>,
) -> Option<(Option<String>, Option<u32>)> {
    match (first, second) {
        (None, None) => Some((None, None)),
        (Some(TemporalSelector::Phase(id)), None) | (None, Some(TemporalSelector::Phase(id))) => {
            Some((Some(id.clone()), None))
        }
        (Some(TemporalSelector::Generation(value)), None)
        | (None, Some(TemporalSelector::Generation(value))) => Some((None, Some(*value))),
        (Some(TemporalSelector::Phase(phase)), Some(TemporalSelector::Generation(generation)))
        | (Some(TemporalSelector::Generation(generation)), Some(TemporalSelector::Phase(phase))) => {
            Some((Some(phase.clone()), Some(*generation)))
        }
        (Some(TemporalSelector::Phase(first)), Some(TemporalSelector::Phase(second)))
            if first == second =>
        {
            Some((Some(first.clone()), None))
        }
        (Some(TemporalSelector::Generation(first)), Some(TemporalSelector::Generation(second)))
            if first == second =>
        {
            Some((None, Some(*first)))
        }
        _ => None,
    }
}

type SpatialContext = (Option<String>, Option<String>, Option<u64>, Option<String>);

fn merge_spatial(
    first: Option<&SpatialSelector>,
    second: Option<&SpatialSelector>,
) -> Option<SpatialContext> {
    match (first, second) {
        (None, None) => Some((None, None, None, None)),
        (Some(selector), None) | (None, Some(selector)) => Some(context_for_spatial(selector)),
        (Some(first), Some(second)) => {
            let (first_reference, first_sequence) = spatial_reference_sequence(first);
            let (second_reference, second_sequence) = spatial_reference_sequence(second);
            if first_reference != second_reference || first_sequence != second_sequence {
                return None;
            }
            let (position, locus) = match (first, second) {
                (SpatialSelector::Sequence { .. }, selector)
                | (selector, SpatialSelector::Sequence { .. }) => {
                    let (_, _, position, locus) = context_for_spatial(selector);
                    (position, locus)
                }
                (SpatialSelector::Region(first), SpatialSelector::Region(second)) => {
                    let start = first.start.max(second.start);
                    (start <= first.end.min(second.end)).then_some((Some(start), None))?
                }
                (SpatialSelector::Region(region), SpatialSelector::Locus { locus, .. })
                | (SpatialSelector::Locus { locus, .. }, SpatialSelector::Region(region)) => {
                    (Some(region.start), Some(locus.clone()))
                }
                (
                    SpatialSelector::Locus { locus: first, .. },
                    SpatialSelector::Locus { locus: second, .. },
                ) if first == second => (None, Some(first.clone())),
                (SpatialSelector::Locus { .. }, SpatialSelector::Locus { .. }) => return None,
            };
            Some((
                Some(first_reference.to_owned()),
                Some(first_sequence.to_owned()),
                position,
                locus,
            ))
        }
    }
}

fn context_for_spatial(selector: &SpatialSelector) -> SpatialContext {
    match selector {
        SpatialSelector::Sequence {
            reference,
            sequence,
        } => (Some(reference.clone()), Some(sequence.clone()), None, None),
        SpatialSelector::Region(region) => (
            Some(region.reference.clone()),
            Some(region.sequence.clone()),
            Some(region.start),
            None,
        ),
        SpatialSelector::Locus {
            reference,
            sequence,
            locus,
        } => (
            Some(reference.clone()),
            Some(sequence.clone()),
            None,
            Some(locus.clone()),
        ),
    }
}

fn spatial_reference_sequence(selector: &SpatialSelector) -> (&str, &str) {
    match selector {
        SpatialSelector::Sequence {
            reference,
            sequence,
        }
        | SpatialSelector::Locus {
            reference,
            sequence,
            ..
        } => (reference, sequence),
        SpatialSelector::Region(region) => (&region.reference, &region.sequence),
    }
}

fn sorted_rule_ids(rules: &[&Rule]) -> String {
    let mut ids = rules
        .iter()
        .map(|rule| rule.id.as_str())
        .collect::<Vec<_>>();
    ids.sort_unstable();
    ids.join(",")
}

/// Indica dominancia estricta entre dos reglas que ya son aplicables al mismo
/// objeto. No establece prioridad entre ejes distintos.
#[must_use]
pub fn dominates(candidate: &SelectorSet, other: &SelectorSet) -> bool {
    let candidate_ranks = specificity(candidate);
    let other_ranks = specificity(other);
    candidate_ranks
        .iter()
        .zip(other_ranks.iter())
        .all(|(candidate, other)| candidate >= other)
        && candidate_ranks
            .iter()
            .zip(other_ranks.iter())
            .any(|(candidate, other)| candidate > other)
}

fn specificity(selectors: &SelectorSet) -> [u8; 7] {
    [
        u8::from(selectors.individual.is_some()),
        match selectors.population {
            None => 0,
            Some(PopulationSelector::Population(_)) => 1,
            Some(PopulationSelector::Subpopulation(_)) => 2,
        },
        match selectors.temporal {
            None => 0,
            Some(TemporalSelector::Phase(_)) => 1,
            Some(TemporalSelector::Generation(_)) => 2,
        },
        match selectors.spatial {
            None => 0,
            Some(SpatialSelector::Sequence { .. }) => 1,
            Some(SpatialSelector::Region(_)) => 2,
            Some(SpatialSelector::Locus { .. }) => 3,
        },
        u8::from(selectors.group.is_some()),
        u8::from(selectors.block.is_some()),
        u8::from(selectors.technical_event.is_some()),
    ]
}

fn selector_matches(declared: Option<&str>, actual: Option<&str>) -> bool {
    declared.is_none_or(|declared| actual == Some(declared))
}

fn population_matches(declared: Option<&PopulationSelector>, context: &ResolutionContext) -> bool {
    match declared {
        None => true,
        Some(PopulationSelector::Population(id)) => context.population.as_ref() == Some(id),
        Some(PopulationSelector::Subpopulation(id)) => context.subpopulation.as_ref() == Some(id),
    }
}

fn temporal_matches(declared: Option<&TemporalSelector>, context: &ResolutionContext) -> bool {
    match declared {
        None => true,
        Some(TemporalSelector::Phase(id)) => context.phase.as_ref() == Some(id),
        Some(TemporalSelector::Generation(generation)) => context.generation == Some(*generation),
    }
}

fn spatial_matches(declared: Option<&SpatialSelector>, context: &ResolutionContext) -> bool {
    match declared {
        None => true,
        Some(SpatialSelector::Sequence {
            reference,
            sequence,
        }) => {
            context.reference.as_ref() == Some(reference)
                && context.sequence.as_ref() == Some(sequence)
        }
        Some(SpatialSelector::Region(region)) => {
            context.reference.as_ref() == Some(&region.reference)
                && context.sequence.as_ref() == Some(&region.sequence)
                && context
                    .position
                    .is_some_and(|position| region.start <= position && position <= region.end)
        }
        Some(SpatialSelector::Locus {
            reference,
            sequence,
            locus,
        }) => {
            context.reference.as_ref() == Some(reference)
                && context.sequence.as_ref() == Some(sequence)
                && context.locus.as_ref() == Some(locus)
        }
    }
}

/// Identidades ya declaradas y disponibles durante la validación previa.
#[derive(Clone, Debug, Default)]
pub struct SelectorIdentityCatalog {
    references: BTreeSet<String>,
    sequences: BTreeSet<(String, String)>,
    individuals: BTreeSet<String>,
    populations: BTreeSet<String>,
    subpopulations: BTreeSet<String>,
    phases: BTreeSet<String>,
    groups: BTreeSet<String>,
    blocks: BTreeSet<String>,
    technical_events: BTreeSet<String>,
}

impl SelectorIdentityCatalog {
    pub fn add_reference(&mut self, id: impl Into<String>) {
        self.references.insert(id.into());
    }

    pub fn add_sequence(&mut self, reference: impl Into<String>, sequence: impl Into<String>) {
        self.sequences.insert((reference.into(), sequence.into()));
    }

    pub fn add_individual(&mut self, id: impl Into<String>) {
        self.individuals.insert(id.into());
    }

    pub fn add_population(&mut self, id: impl Into<String>) {
        self.populations.insert(id.into());
    }

    pub fn add_subpopulation(&mut self, id: impl Into<String>) {
        self.subpopulations.insert(id.into());
    }

    pub fn add_phase(&mut self, id: impl Into<String>) {
        self.phases.insert(id.into());
    }

    pub fn add_group(&mut self, id: impl Into<String>) {
        self.groups.insert(id.into());
    }

    pub fn add_block(&mut self, id: impl Into<String>) {
        self.blocks.insert(id.into());
    }

    pub fn add_technical_event(&mut self, id: impl Into<String>) {
        self.technical_events.insert(id.into());
    }

    fn category_of(&self, id: &str) -> Option<&'static str> {
        if self.references.contains(id) {
            Some("reference")
        } else if self.individuals.contains(id) {
            Some("individual")
        } else if self.populations.contains(id) {
            Some("population")
        } else if self.subpopulations.contains(id) {
            Some("subpopulation")
        } else if self.phases.contains(id) {
            Some("phase")
        } else if self.groups.contains(id) {
            Some("group")
        } else if self.blocks.contains(id) {
            Some("block")
        } else if self.technical_events.contains(id) {
            Some("technical_event")
        } else {
            None
        }
    }
}

pub fn validate_rules(
    rules: &[Rule],
    identities: &SelectorIdentityCatalog,
) -> Result<(), Diagnostic> {
    validate_exact_duplicates(rules)?;
    validate_regions(rules)?;
    for rule in rules {
        validate_rule_identities(rule, identities)?;
    }
    Ok(())
}

fn validate_exact_duplicates(rules: &[Rule]) -> Result<(), Diagnostic> {
    let mut seen: BTreeMap<(ParameterKey, &SelectorSet), &Rule> = BTreeMap::new();
    for rule in rules {
        if let Some(previous) = seen.insert((rule.parameter, &rule.selectors), rule) {
            return Err(Diagnostic::configuration(
                "GPS-CONF-010",
                "Dos reglas del mismo parámetro tienen exactamente los mismos selectores.",
            )
            .with("parameter", rule.parameter.as_str())
            .with("selectors", format!("{:?}", rule.selectors))
            .with("scope", format!("{:?}", rule.selectors))
            .with(
                "locations",
                format!(
                    "{}:{}; {}:{}",
                    previous.provenance.path,
                    previous.provenance.location,
                    rule.provenance.path,
                    rule.provenance.location
                ),
            ));
        }
    }
    Ok(())
}

fn validate_regions(rules: &[Rule]) -> Result<(), Diagnostic> {
    for rule in rules {
        if let Some(SpatialSelector::Region(region)) = &rule.selectors.spatial
            && (region.start == 0 || region.start > region.end)
        {
            return Err(Diagnostic::configuration(
                "GPS-CONF-008",
                "La configuración contiene una combinación de campos incompatible.",
            )
            .with("fields", "region_start,region_end")
            .with(
                "detail",
                format!(
                    "intervalo cerrado inválido: {}..{}",
                    region.start, region.end
                ),
            ));
        }
    }
    for (index, first) in rules.iter().enumerate() {
        let Some(SpatialSelector::Region(first_region)) = &first.selectors.spatial else {
            continue;
        };
        for second in &rules[index + 1..] {
            let Some(SpatialSelector::Region(second_region)) = &second.selectors.spatial else {
                continue;
            };
            if first.parameter == second.parameter
                && same_non_spatial_selectors(&first.selectors, &second.selectors)
                && first_region.reference == second_region.reference
                && first_region.sequence == second_region.sequence
                && first_region.start <= second_region.end
                && second_region.start <= first_region.end
            {
                return Err(Diagnostic::configuration(
                    "GPS-CONF-011",
                    "Dos regiones del mismo parámetro y selectores no espaciales se solapan.",
                )
                .with("parameter", first.parameter.as_str())
                .with("reference", &first_region.reference)
                .with("sequence", &first_region.sequence)
                .with("scope", format!("{:?}", first.selectors))
                .with(
                    "regions",
                    format!(
                        "{}..{};{}..{}",
                        first_region.start,
                        first_region.end,
                        second_region.start,
                        second_region.end
                    ),
                ));
            }
        }
    }
    Ok(())
}

fn same_non_spatial_selectors(first: &SelectorSet, second: &SelectorSet) -> bool {
    first.individual == second.individual
        && first.population == second.population
        && first.temporal == second.temporal
        && first.group == second.group
        && first.block == second.block
        && first.technical_event == second.technical_event
}

fn validate_rule_identities(
    rule: &Rule,
    identities: &SelectorIdentityCatalog,
) -> Result<(), Diagnostic> {
    if let Some(id) = &rule.selectors.individual {
        require_identity(identities, "individual", id, &identities.individuals)?;
    }
    if let Some(population) = &rule.selectors.population {
        match population {
            PopulationSelector::Population(id) => {
                require_identity(identities, "population", id, &identities.populations)?;
            }
            PopulationSelector::Subpopulation(id) => {
                require_identity(identities, "subpopulation", id, &identities.subpopulations)?
            }
        }
    }
    if let Some(TemporalSelector::Phase(id)) = &rule.selectors.temporal {
        require_identity(identities, "phase", id, &identities.phases)?;
    }
    if let Some(spatial) = &rule.selectors.spatial {
        let (reference, sequence) = match spatial {
            SpatialSelector::Sequence {
                reference,
                sequence,
            }
            | SpatialSelector::Locus {
                reference,
                sequence,
                ..
            } => (reference, sequence),
            SpatialSelector::Region(region) => (&region.reference, &region.sequence),
        };
        require_identity(identities, "reference", reference, &identities.references)?;
        if !identities
            .sequences
            .contains(&(reference.clone(), sequence.clone()))
        {
            return Err(identity_error(
                identities,
                "sequence",
                &format!("{reference}:{sequence}"),
            ));
        }
    }
    if let Some(id) = &rule.selectors.group {
        require_identity(identities, "group", id, &identities.groups)?;
    }
    if let Some(id) = &rule.selectors.block {
        require_identity(identities, "block", id, &identities.blocks)?;
    }
    if let Some(id) = &rule.selectors.technical_event {
        require_identity(
            identities,
            "technical_event",
            id,
            &identities.technical_events,
        )?;
    }
    Ok(())
}

fn require_identity(
    identities: &SelectorIdentityCatalog,
    expected: &'static str,
    id: &str,
    category: &BTreeSet<String>,
) -> Result<(), Diagnostic> {
    if category.contains(id) {
        Ok(())
    } else {
        Err(identity_error(identities, expected, id))
    }
}

fn identity_error(
    identities: &SelectorIdentityCatalog,
    expected: &'static str,
    id: &str,
) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-013",
        "Una regla referencia una identidad inexistente o de otra categoría.",
    )
    .with("field", expected)
    .with("identifier", id)
    .with("expected_category", expected)
    .with(
        "found_category",
        identities.category_of(id).unwrap_or("missing"),
    )
}
