//! Contratos propios para configuración declarada y efectiva.

use std::fmt;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::diagnostics::Diagnostic;
pub use crate::reference::SequenceClass;
use crate::rules::Rule;

/// Clase de un documento TOML propio.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConfigKind {
    General,
    Scenario,
}

impl ConfigKind {
    pub(crate) fn parse(value: &str) -> Option<Self> {
        match value {
            "general" => Some(Self::General),
            "scenario" => Some(Self::Scenario),
            _ => None,
        }
    }

    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::General => "general",
            Self::Scenario => "scenario",
        }
    }
}

/// Clase exigida por la frontera que recibe un documento.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ExpectedConfigKind {
    General,
    Scenario,
    Either,
}

impl ExpectedConfigKind {
    pub(crate) const fn accepts(self, actual: ConfigKind) -> bool {
        matches!(self, Self::Either)
            || matches!(
                (self, actual),
                (Self::General, ConfigKind::General) | (Self::Scenario, ConfigKind::Scenario)
            )
    }

    pub(crate) const fn as_str(self) -> &'static str {
        match self {
            Self::General => "general",
            Self::Scenario => "scenario",
            Self::Either => "general|scenario",
        }
    }
}

/// Documento cuya sintaxis, clase y versión ya superaron la frontera inicial.
#[derive(Clone, Debug)]
pub struct AdmittedToml {
    path: PathBuf,
    bytes: Arc<[u8]>,
    kind: ConfigKind,
}

impl AdmittedToml {
    pub(crate) fn new(path: PathBuf, bytes: Arc<[u8]>, kind: ConfigKind) -> Self {
        Self { path, bytes, kind }
    }

    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    #[must_use]
    pub const fn kind(&self) -> ConfigKind {
        self.kind
    }
}

/// Ubicación de un valor o campo sobre los bytes originales.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SourceLocation {
    pub byte_start: usize,
    pub byte_end: usize,
    pub line: usize,
    pub column: usize,
}

impl SourceLocation {
    #[must_use]
    pub fn from_span(bytes: &[u8], span: Range<usize>) -> Self {
        let mut location = Self::from_offset(bytes, span.start);
        location.byte_end = span.end;
        location
    }

    #[must_use]
    pub fn from_offset(bytes: &[u8], offset: usize) -> Self {
        let safe_offset = offset.min(bytes.len());
        let line_start = bytes[..safe_offset]
            .iter()
            .rposition(|byte| *byte == b'\n')
            .map_or(0, |index| index + 1);
        Self {
            byte_start: safe_offset,
            byte_end: safe_offset,
            line: bytes[..safe_offset]
                .iter()
                .filter(|byte| **byte == b'\n')
                .count()
                + 1,
            column: String::from_utf8_lossy(&bytes[line_start..safe_offset])
                .chars()
                .count()
                + 1,
        }
    }
}

impl fmt::Display for SourceLocation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "line {}, column {}, bytes {}..{}",
            self.line, self.column, self.byte_start, self.byte_end
        )
    }
}

/// Valor declarado junto con su ubicación exacta en el TOML.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Located<T> {
    value: T,
    location: SourceLocation,
}

impl<T> Located<T> {
    pub(crate) const fn new(value: T, location: SourceLocation) -> Self {
        Self { value, location }
    }

    #[must_use]
    pub const fn value(&self) -> &T {
        &self.value
    }

    #[must_use]
    pub const fn location(&self) -> &SourceLocation {
        &self.location
    }
}

/// Vía biológica declarada.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Route {
    Direct,
    Evolution,
    Families,
}

impl Route {
    pub(crate) fn parse(value: &str) -> Option<Self> {
        match value {
            "direct" => Some(Self::Direct),
            "evolution" => Some(Self::Evolution),
            "families" => Some(Self::Families),
            _ => None,
        }
    }

    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Direct => "direct",
            Self::Evolution => "evolution",
            Self::Families => "families",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredEntity {
    id: Located<String>,
}

impl DeclaredEntity {
    pub(crate) const fn new(id: Located<String>) -> Self {
        Self { id }
    }

    #[must_use]
    pub const fn id(&self) -> &Located<String> {
        &self.id
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredSequence {
    pub id: Located<String>,
    pub class: Located<SequenceClass>,
    pub length: Located<u64>,
    pub ploidy: Option<Located<u8>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DeclaredReferenceSource {
    Simulated {
        default_ploidy: Located<u8>,
        composition: Vec<DeclaredFrequencyValue>,
    },
    Supplied {
        fasta: Located<String>,
    },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredNumber {
    literal: String,
    location: SourceLocation,
}

impl DeclaredNumber {
    pub(crate) const fn new(literal: String, location: SourceLocation) -> Self {
        Self { literal, location }
    }

    #[must_use]
    pub fn literal(&self) -> &str {
        &self.literal
    }

    #[must_use]
    pub const fn location(&self) -> &SourceLocation {
        &self.location
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredFrequencyValue {
    pub id: Located<String>,
    pub category: Located<String>,
    pub value: DeclaredNumber,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FrequencyMode {
    Probabilistic,
    Controlled,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredFrequencyCount {
    pub id: Located<String>,
    pub category: Located<String>,
    pub count: Located<u64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredFrequencySet {
    pub id: Located<String>,
    pub mode: Located<FrequencyMode>,
    pub tolerance: Option<DeclaredNumber>,
    pub values: Vec<DeclaredFrequencyValue>,
    pub counts: Vec<DeclaredFrequencyCount>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredReference {
    pub id: Located<String>,
    pub source: DeclaredReferenceSource,
    pub sequences: Vec<DeclaredSequence>,
}

/// Objetivo declarado de una referencia reducida. Las variantes son
/// excluyentes por construcción después de atravesar el adaptador TOML.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DeclaredReductionTarget {
    Bases(Located<u64>),
    ReductionPercent(DeclaredNumber),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredReducedReference {
    pub id: Located<String>,
    pub source_reference_id: Located<String>,
    pub region_size_bases: Located<u64>,
    pub target: DeclaredReductionTarget,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ScenarioLink {
    pub id: Located<String>,
    pub path: Located<String>,
    pub route: Located<Route>,
}

#[derive(Clone, Debug)]
pub struct DeclaredGeneralConfig {
    source: AdmittedToml,
    route: Located<Route>,
    seed: Option<Located<String>>,
    normalization_tolerance: Option<DeclaredNumber>,
    reference: DeclaredReference,
    reduced_references: Vec<DeclaredReducedReference>,
    scenarios: Vec<ScenarioLink>,
    events: Vec<DeclaredEntity>,
    individuals: Vec<DeclaredEntity>,
    populations: Vec<DeclaredEntity>,
    subpopulations: Vec<DeclaredEntity>,
    phases: Vec<DeclaredEntity>,
    groups: Vec<DeclaredEntity>,
    blocks: Vec<DeclaredEntity>,
    rules: Vec<Rule>,
    frequency_sets: Vec<DeclaredFrequencySet>,
}

pub(crate) struct DeclaredGeneralParts {
    pub source: AdmittedToml,
    pub route: Located<Route>,
    pub seed: Option<Located<String>>,
    pub normalization_tolerance: Option<DeclaredNumber>,
    pub reference: DeclaredReference,
    pub reduced_references: Vec<DeclaredReducedReference>,
    pub scenarios: Vec<ScenarioLink>,
    pub events: Vec<DeclaredEntity>,
    pub individuals: Vec<DeclaredEntity>,
    pub populations: Vec<DeclaredEntity>,
    pub subpopulations: Vec<DeclaredEntity>,
    pub phases: Vec<DeclaredEntity>,
    pub groups: Vec<DeclaredEntity>,
    pub blocks: Vec<DeclaredEntity>,
    pub rules: Vec<Rule>,
    pub frequency_sets: Vec<DeclaredFrequencySet>,
}

impl DeclaredGeneralConfig {
    pub(crate) fn new(parts: DeclaredGeneralParts) -> Self {
        let DeclaredGeneralParts {
            source,
            route,
            seed,
            normalization_tolerance,
            reference,
            reduced_references,
            scenarios,
            events,
            individuals,
            populations,
            subpopulations,
            phases,
            groups,
            blocks,
            rules,
            frequency_sets,
        } = parts;
        Self {
            source,
            route,
            seed,
            normalization_tolerance,
            reference,
            reduced_references,
            scenarios,
            events,
            individuals,
            populations,
            subpopulations,
            phases,
            groups,
            blocks,
            rules,
            frequency_sets,
        }
    }

    #[must_use]
    pub const fn source(&self) -> &AdmittedToml {
        &self.source
    }

    #[must_use]
    pub const fn route(&self) -> &Located<Route> {
        &self.route
    }

    #[must_use]
    pub const fn seed(&self) -> Option<&Located<String>> {
        self.seed.as_ref()
    }

    #[must_use]
    pub const fn normalization_tolerance(&self) -> Option<&DeclaredNumber> {
        self.normalization_tolerance.as_ref()
    }

    #[must_use]
    pub const fn reference(&self) -> &DeclaredReference {
        &self.reference
    }

    #[must_use]
    pub fn reduced_references(&self) -> &[DeclaredReducedReference] {
        &self.reduced_references
    }

    #[must_use]
    pub fn scenarios(&self) -> &[ScenarioLink] {
        &self.scenarios
    }

    #[must_use]
    pub fn events(&self) -> &[DeclaredEntity] {
        &self.events
    }

    #[must_use]
    pub fn individuals(&self) -> &[DeclaredEntity] {
        &self.individuals
    }

    #[must_use]
    pub fn populations(&self) -> &[DeclaredEntity] {
        &self.populations
    }

    #[must_use]
    pub fn subpopulations(&self) -> &[DeclaredEntity] {
        &self.subpopulations
    }

    #[must_use]
    pub fn phases(&self) -> &[DeclaredEntity] {
        &self.phases
    }

    #[must_use]
    pub fn groups(&self) -> &[DeclaredEntity] {
        &self.groups
    }

    #[must_use]
    pub fn blocks(&self) -> &[DeclaredEntity] {
        &self.blocks
    }

    #[must_use]
    pub fn rules(&self) -> &[Rule] {
        &self.rules
    }

    #[must_use]
    pub fn frequency_sets(&self) -> &[DeclaredFrequencySet] {
        &self.frequency_sets
    }
}

#[derive(Clone, Debug)]
pub struct DeclaredScenarioConfig {
    source: AdmittedToml,
    id: Located<String>,
    route: Located<Route>,
    phases: Vec<DeclaredEntity>,
    events: Vec<DeclaredEntity>,
    crosses: Vec<DeclaredEntity>,
}

impl DeclaredScenarioConfig {
    pub(crate) fn new(
        source: AdmittedToml,
        id: Located<String>,
        route: Located<Route>,
        phases: Vec<DeclaredEntity>,
        events: Vec<DeclaredEntity>,
        crosses: Vec<DeclaredEntity>,
    ) -> Self {
        Self {
            source,
            id,
            route,
            phases,
            events,
            crosses,
        }
    }

    #[must_use]
    pub const fn source(&self) -> &AdmittedToml {
        &self.source
    }

    #[must_use]
    pub const fn id(&self) -> &Located<String> {
        &self.id
    }

    #[must_use]
    pub const fn route(&self) -> &Located<Route> {
        &self.route
    }

    #[must_use]
    pub fn phases(&self) -> &[DeclaredEntity] {
        &self.phases
    }

    #[must_use]
    pub fn events(&self) -> &[DeclaredEntity] {
        &self.events
    }

    #[must_use]
    pub fn crosses(&self) -> &[DeclaredEntity] {
        &self.crosses
    }
}

#[derive(Clone, Debug)]
pub enum DeclaredConfig {
    General(Box<DeclaredGeneralConfig>),
    Scenario(Box<DeclaredScenarioConfig>),
}

impl DeclaredConfig {
    #[must_use]
    pub fn source(&self) -> &AdmittedToml {
        match self {
            Self::General(config) => config.source(),
            Self::Scenario(config) => config.source(),
        }
    }
}

pub(crate) fn semantic_error(
    source: &AdmittedToml,
    fields: impl Into<String>,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", source.path().display().to_string())
    .with("fields", fields.into())
    .with("detail", detail.into())
}
