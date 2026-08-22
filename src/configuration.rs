//! Contrato TOML mínimo de F1 para referencias iniciales y reducidas.

use std::fmt;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::diagnostics::Diagnostic;
pub use crate::reference::SequenceClass;

/// El recorrido vigente solo admite el documento general.
/// Los escenarios se incorporarán con su primer consumidor biológico real.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConfigKind {
    General,
}

impl ConfigKind {
    pub(crate) fn parse(value: &str) -> Option<Self> {
        (value == "general").then_some(Self::General)
    }

    #[must_use]
    pub const fn as_str(self) -> &'static str {
        "general"
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ExpectedConfigKind {
    General,
}

impl ExpectedConfigKind {
    pub(crate) const fn accepts(self, actual: ConfigKind) -> bool {
        matches!((self, actual), (Self::General, ConfigKind::General))
    }

    pub(crate) const fn as_str(self) -> &'static str {
        "general"
    }
}

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
pub struct DeclaredCompositionValue {
    pub base: Located<String>,
    pub value: DeclaredNumber,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredSequence {
    pub id: Located<String>,
    pub class: Located<SequenceClass>,
    pub length: Located<u64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DeclaredReferenceSource {
    Simulated {
        composition: Vec<DeclaredCompositionValue>,
    },
    Supplied {
        fasta: Located<String>,
    },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredReference {
    pub id: Located<String>,
    pub source: DeclaredReferenceSource,
    pub sequences: Vec<DeclaredSequence>,
}

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

/// Selección declarada del recurso común de frecuencias SNV.
///
/// Solo describe la frontera de autoría. La preparación, la reutilización
/// verificable y los parámetros de generación tienen sus propios pasos F2.2.
#[derive(Clone, Debug, PartialEq)]
pub enum DeclaredPopulationSnvs {
    Supplied {
        input_path: Located<String>,
    },
    Reuse {
        metadata_path: Located<String>,
    },
    Generate {
        parameters: Box<DeclaredPopulationSnvsGeneration>,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredPopulationSnvsGeneration {
    pub count_mode: Located<String>,
    pub locus_count: Option<Located<u64>>,
    pub mean_loci: Option<DeclaredNumber>,
    pub mean_density_per_base: Option<DeclaredNumber>,
    pub spatial_mode: Located<String>,
    pub substitution_profile: Located<String>,
    pub outside_weight: Option<DeclaredNumber>,
    pub regions: Vec<DeclaredPopulationSnvsRegion>,
    pub populations: Vec<DeclaredPopulationSnvsPopulation>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredPopulationSnvsRegion {
    pub contig: Located<String>,
    pub start: Located<u64>,
    pub end: Located<u64>,
    pub weight: DeclaredNumber,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredPopulationSnvsPopulation {
    pub population_id: Located<String>,
    pub frequency_model: Located<String>,
    pub frequency_mode: Option<Located<String>>,
    pub min_alt_af: Option<DeclaredNumber>,
    pub max_alt_af: Option<DeclaredNumber>,
    pub spectrum: Vec<DeclaredFrequencySpectrumComponent>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredFrequencySpectrumComponent {
    pub kind: Located<String>,
    pub fixed_alt_af: Option<DeclaredNumber>,
    pub min_alt_af: Option<DeclaredNumber>,
    pub max_alt_af: Option<DeclaredNumber>,
    pub weight: DeclaredNumber,
    pub sampling: Option<Located<String>>,
}

#[derive(Clone, Debug)]
pub struct DeclaredGeneralConfig {
    source: AdmittedToml,
    seed: Option<Located<String>>,
    ploidy: Option<Located<u64>>,
    normalization_tolerance: Option<DeclaredNumber>,
    reference: DeclaredReference,
    reduced_references: Vec<DeclaredReducedReference>,
    population_snvs: Option<DeclaredPopulationSnvs>,
}

impl DeclaredGeneralConfig {
    pub(crate) fn new(
        source: AdmittedToml,
        seed: Option<Located<String>>,
        ploidy: Option<Located<u64>>,
        normalization_tolerance: Option<DeclaredNumber>,
        reference: DeclaredReference,
        reduced_references: Vec<DeclaredReducedReference>,
        population_snvs: Option<DeclaredPopulationSnvs>,
    ) -> Self {
        Self {
            source,
            seed,
            ploidy,
            normalization_tolerance,
            reference,
            reduced_references,
            population_snvs,
        }
    }

    #[must_use]
    pub const fn source(&self) -> &AdmittedToml {
        &self.source
    }

    #[must_use]
    pub const fn seed(&self) -> Option<&Located<String>> {
        self.seed.as_ref()
    }

    /// Ploidía global solicitada; ausente conserva el valor diploide compatible con F1.
    #[must_use]
    pub const fn ploidy(&self) -> Option<&Located<u64>> {
        self.ploidy.as_ref()
    }

    /// Tolerancia general DTC-049 solicitada para las frecuencias de F2.
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

    /// Fuente de matriz SNV solicitada, si la configuración la declara.
    #[must_use]
    pub const fn population_snvs(&self) -> Option<&DeclaredPopulationSnvs> {
        self.population_snvs.as_ref()
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
