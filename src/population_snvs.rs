//! Formatos estructurales del recurso canónico de frecuencias SNV.
//!
//! Este módulo solo fija la envoltura TSV+JSON de F2.2.1. La validación contra
//! FASTA, la semántica de `alt_af`, sumas, IDs canónicos y hashes operativos se
//! incorpora en los pasos consumidores posteriores.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;

use serde_json::Value;

/// Cabecera literal y ordenada del TSV canónico de matriz SNV.
pub const TSV_HEADER: [&str; 8] = [
    "variant_id",
    "contig",
    "position",
    "ref",
    "alt",
    "population_id",
    "alt_af",
    "source_id",
];

/// Orden reproducible usado para sumar los ALT de un locus multialélico.
pub const ALT_FREQUENCY_CANONICAL_ORDER: &str =
    "population_id,contig,position,ref,alt,variant_id,source_id";
/// Algoritmo y versión de suma auditados por DTC-049.
pub const ALT_FREQUENCY_SUM_ALGORITHM: &str = "neumaier_sum_v1";

type AltGroupKey<'a> = (&'a str, &'a str, &'a str, &'a str);
type ParsedAltRow<'a> = (usize, &'a PopulationSnvRow, f64);

/// Un registro tabular aún no interpretado semánticamente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PopulationSnvRow {
    fields: [String; 8],
}

impl PopulationSnvRow {
    #[must_use]
    pub const fn fields(&self) -> &[String; 8] {
        &self.fields
    }

    #[must_use]
    pub fn field(&self, index: usize) -> &str {
        &self.fields[index]
    }
}

/// TSV canónico con su cabecera comprobada y filas de ocho columnas.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PopulationSnvsTsv {
    rows: Vec<PopulationSnvRow>,
}

impl PopulationSnvsTsv {
    #[must_use]
    pub fn rows(&self) -> &[PopulationSnvRow] {
        &self.rows
    }
}

/// Error explícito de formato anterior a cualquier interpretación biológica.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PopulationSnvsFormatError {
    InvalidUtf8,
    VcfInput,
    MissingHeader,
    IncorrectHeader {
        found: String,
    },
    IncorrectColumnCount {
        line: usize,
        found: usize,
    },
    InvalidMetadataJson,
    MetadataMustBeObject,
    MissingMetadataField {
        field: String,
    },
    UnexpectedMetadataField {
        field: String,
    },
    InvalidMetadataField {
        field: String,
        expected: &'static str,
    },
    InvalidMetadataLiteral {
        field: String,
        expected: &'static str,
    },
}

impl fmt::Display for PopulationSnvsFormatError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidUtf8 => write!(formatter, "el recurso no es UTF-8"),
            Self::VcfInput => write!(formatter, "VCF no es una entrada de matriz SNV"),
            Self::MissingHeader => write!(formatter, "falta la cabecera TSV canónica"),
            Self::IncorrectHeader { found } => {
                write!(formatter, "cabecera TSV no canónica: {found}")
            }
            Self::IncorrectColumnCount { line, found } => {
                write!(
                    formatter,
                    "la fila {line} contiene {found} columnas; se esperaban 8"
                )
            }
            Self::InvalidMetadataJson => write!(formatter, "metadatos JSON inválidos"),
            Self::MetadataMustBeObject => {
                write!(formatter, "los metadatos deben ser un objeto JSON")
            }
            Self::MissingMetadataField { field } => {
                write!(formatter, "falta el campo de metadatos {field}")
            }
            Self::UnexpectedMetadataField { field } => {
                write!(formatter, "campo de metadatos no admitido: {field}")
            }
            Self::InvalidMetadataField { field, expected } => {
                write!(formatter, "{field} debe ser {expected}")
            }
            Self::InvalidMetadataLiteral { field, expected } => {
                write!(formatter, "{field} debe ser {expected}")
            }
        }
    }
}

impl std::error::Error for PopulationSnvsFormatError {}

/// Decisión aplicada a un grupo población×locus.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AltFrequencyDecision {
    ReferenceRemainder,
    Normalized,
}

impl AltFrequencyDecision {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ReferenceRemainder => "reference_remainder",
            Self::Normalized => "normalized",
        }
    }
}

/// Valor ALT conservado en el orden canónico, antes y después de la decisión.
#[derive(Clone, Debug, PartialEq)]
pub struct AuditedAltFrequency {
    row_index: usize,
    alt: String,
    variant_id: String,
    source_id: String,
    original: f64,
    effective: f64,
}

impl AuditedAltFrequency {
    #[must_use]
    pub const fn row_index(&self) -> usize {
        self.row_index
    }
    #[must_use]
    pub fn alt(&self) -> &str {
        &self.alt
    }
    #[must_use]
    pub fn variant_id(&self) -> &str {
        &self.variant_id
    }
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }
    #[must_use]
    pub const fn original(&self) -> f64 {
        self.original
    }
    #[must_use]
    pub const fn effective(&self) -> f64 {
        self.effective
    }
}

/// Auditoría completa de una suma ALT población×locus.
#[derive(Clone, Debug, PartialEq)]
pub struct AltFrequencyNormalization {
    population_id: String,
    contig: String,
    position: String,
    reference: String,
    tolerance: f64,
    canonical_order: &'static str,
    sum_algorithm: &'static str,
    original_sum: f64,
    original_delta: f64,
    effective_sum: f64,
    effective_delta: f64,
    reference_frequency: f64,
    decision: AltFrequencyDecision,
    values: Vec<AuditedAltFrequency>,
}

impl AltFrequencyNormalization {
    #[must_use]
    pub fn population_id(&self) -> &str {
        &self.population_id
    }
    #[must_use]
    pub fn contig(&self) -> &str {
        &self.contig
    }
    #[must_use]
    pub fn position(&self) -> &str {
        &self.position
    }
    #[must_use]
    pub fn reference(&self) -> &str {
        &self.reference
    }
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }
    #[must_use]
    pub const fn canonical_order(&self) -> &'static str {
        self.canonical_order
    }
    #[must_use]
    pub const fn sum_algorithm(&self) -> &'static str {
        self.sum_algorithm
    }
    #[must_use]
    pub const fn original_sum(&self) -> f64 {
        self.original_sum
    }
    #[must_use]
    pub const fn original_delta(&self) -> f64 {
        self.original_delta
    }
    #[must_use]
    pub const fn effective_sum(&self) -> f64 {
        self.effective_sum
    }
    #[must_use]
    pub const fn effective_delta(&self) -> f64 {
        self.effective_delta
    }
    #[must_use]
    pub const fn reference_frequency(&self) -> f64 {
        self.reference_frequency
    }
    #[must_use]
    pub const fn decision(&self) -> AltFrequencyDecision {
        self.decision
    }
    #[must_use]
    pub fn values(&self) -> &[AuditedAltFrequency] {
        &self.values
    }
}

/// Resultado normalizado, sin afirmar todavía que la tabla sea un recurso preparado.
#[derive(Clone, Debug, PartialEq)]
pub struct NormalizedPopulationSnvs {
    groups: Vec<AltFrequencyNormalization>,
}

impl NormalizedPopulationSnvs {
    #[must_use]
    pub fn groups(&self) -> &[AltFrequencyNormalization] {
        &self.groups
    }
}

/// Error de valores y sumas ALT, separado de la forma TSV de F2.2.1.
#[derive(Clone, Debug, PartialEq)]
pub enum AltFrequencyNormalizationError {
    InvalidTolerance {
        value: f64,
    },
    InvalidFrequency {
        row_index: usize,
        literal: String,
    },
    ExcessAltSum {
        population_id: String,
        contig: String,
        position: String,
        reference: String,
        sum: f64,
        delta: f64,
        tolerance: f64,
    },
    InvalidNormalizedSum {
        population_id: String,
        contig: String,
        position: String,
        reference: String,
    },
}

impl fmt::Display for AltFrequencyNormalizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidTolerance { value } => {
                write!(formatter, "tolerancia general inválida: {value}")
            }
            Self::InvalidFrequency { row_index, literal } => write!(
                formatter,
                "alt_af inválido en fila {}: {literal}",
                row_index + 2
            ),
            Self::ExcessAltSum {
                population_id,
                contig,
                position,
                reference,
                sum,
                delta,
                tolerance,
            } => write!(
                formatter,
                "suma ALT excesiva para {population_id}/{contig}:{position}:{reference}: suma={sum}, delta={delta}, tolerancia={tolerance}"
            ),
            Self::InvalidNormalizedSum {
                population_id,
                contig,
                position,
                reference,
            } => write!(
                formatter,
                "la normalización ALT no produjo una suma válida para {population_id}/{contig}:{position}:{reference}"
            ),
        }
    }
}

impl std::error::Error for AltFrequencyNormalizationError {}

/// Valida y normaliza las sumas ALT multialélicas conforme a DTC-049/DTC-054.
///
/// Las filas se agrupan por población y locus textual, y sus ALT se recorren en
/// orden canónico estable. No interpreta aún coordenadas, alelos, IDs, FASTA ni
/// enlaces de recurso: esas autoridades pertenecen a los pasos posteriores.
pub fn normalize_alt_frequencies(
    table: &PopulationSnvsTsv,
    tolerance: f64,
) -> Result<NormalizedPopulationSnvs, AltFrequencyNormalizationError> {
    if !tolerance.is_finite() || !(0.0..=1e-6).contains(&tolerance) {
        return Err(AltFrequencyNormalizationError::InvalidTolerance { value: tolerance });
    }
    let mut groups: BTreeMap<AltGroupKey<'_>, Vec<ParsedAltRow<'_>>> = BTreeMap::new();
    for (row_index, row) in table.rows.iter().enumerate() {
        let literal = row.field(6);
        let value = literal.parse::<f64>().map_err(|_| {
            AltFrequencyNormalizationError::InvalidFrequency {
                row_index,
                literal: literal.to_owned(),
            }
        })?;
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(AltFrequencyNormalizationError::InvalidFrequency {
                row_index,
                literal: literal.to_owned(),
            });
        }
        groups
            .entry((row.field(5), row.field(1), row.field(2), row.field(3)))
            .or_default()
            .push((row_index, row, value));
    }

    let mut audited = Vec::with_capacity(groups.len());
    for ((population_id, contig, position, reference), mut rows) in groups {
        rows.sort_by(|left, right| canonical_alt_order(left.1, right.1));
        let original_sum = neumaier(rows.iter().map(|(_, _, value)| *value));
        let original_delta = (original_sum - 1.0).abs();
        let decision = if original_sum <= 1.0 {
            AltFrequencyDecision::ReferenceRemainder
        } else if original_sum.is_finite() && original_sum > 0.0 && original_delta <= tolerance {
            AltFrequencyDecision::Normalized
        } else {
            return Err(AltFrequencyNormalizationError::ExcessAltSum {
                population_id: population_id.to_owned(),
                contig: contig.to_owned(),
                position: position.to_owned(),
                reference: reference.to_owned(),
                sum: original_sum,
                delta: original_delta,
                tolerance,
            });
        };
        let values = rows
            .into_iter()
            .map(|(row_index, row, original)| {
                let effective = match decision {
                    AltFrequencyDecision::ReferenceRemainder => original,
                    AltFrequencyDecision::Normalized => original / original_sum,
                };
                AuditedAltFrequency {
                    row_index,
                    alt: row.field(4).to_owned(),
                    variant_id: row.field(0).to_owned(),
                    source_id: row.field(7).to_owned(),
                    original,
                    effective,
                }
            })
            .collect::<Vec<_>>();
        if values
            .iter()
            .any(|value| !value.effective.is_finite() || !(0.0..=1.0).contains(&value.effective))
        {
            return Err(AltFrequencyNormalizationError::InvalidNormalizedSum {
                population_id: population_id.to_owned(),
                contig: contig.to_owned(),
                position: position.to_owned(),
                reference: reference.to_owned(),
            });
        }
        let effective_sum = neumaier(values.iter().map(AuditedAltFrequency::effective));
        let effective_delta = (effective_sum - 1.0).abs();
        if decision == AltFrequencyDecision::Normalized
            && (!effective_sum.is_finite() || effective_delta > tolerance)
        {
            return Err(AltFrequencyNormalizationError::InvalidNormalizedSum {
                population_id: population_id.to_owned(),
                contig: contig.to_owned(),
                position: position.to_owned(),
                reference: reference.to_owned(),
            });
        }
        audited.push(AltFrequencyNormalization {
            population_id: population_id.to_owned(),
            contig: contig.to_owned(),
            position: position.to_owned(),
            reference: reference.to_owned(),
            tolerance,
            canonical_order: ALT_FREQUENCY_CANONICAL_ORDER,
            sum_algorithm: ALT_FREQUENCY_SUM_ALGORITHM,
            original_sum,
            original_delta,
            effective_sum,
            effective_delta,
            reference_frequency: match decision {
                AltFrequencyDecision::ReferenceRemainder => 1.0 - effective_sum,
                // Tras normalizar, la semántica DTC-054 fija REF exactamente a cero;
                // recomputar 1-suma puede dejar un residuo negativo por redondeo.
                AltFrequencyDecision::Normalized => 0.0,
            },
            decision,
            values,
        });
    }
    Ok(NormalizedPopulationSnvs { groups: audited })
}

fn canonical_alt_order(left: &PopulationSnvRow, right: &PopulationSnvRow) -> Ordering {
    left.field(4)
        .cmp(right.field(4))
        .then_with(|| left.field(0).cmp(right.field(0)))
        .then_with(|| left.field(7).cmp(right.field(7)))
}

fn neumaier(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

/// Comprueba la forma del TSV canónico; no interpreta aún sus valores.
pub fn parse_canonical_tsv(bytes: &[u8]) -> Result<PopulationSnvsTsv, PopulationSnvsFormatError> {
    let text = std::str::from_utf8(bytes).map_err(|_| PopulationSnvsFormatError::InvalidUtf8)?;
    if text.starts_with("##fileformat=VCF") || text.starts_with("#CHROM\t") {
        return Err(PopulationSnvsFormatError::VcfInput);
    }
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or(PopulationSnvsFormatError::MissingHeader)?;
    let expected = TSV_HEADER.join("\t");
    if header != expected {
        return Err(PopulationSnvsFormatError::IncorrectHeader {
            found: header.to_owned(),
        });
    }
    let mut rows = Vec::new();
    for (index, line) in lines.enumerate() {
        if line.is_empty() {
            return Err(PopulationSnvsFormatError::IncorrectColumnCount {
                line: index + 2,
                found: 0,
            });
        }
        let fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
        let found = fields.len();
        let fields =
            fields
                .try_into()
                .map_err(|_| PopulationSnvsFormatError::IncorrectColumnCount {
                    line: index + 2,
                    found,
                })?;
        rows.push(PopulationSnvRow { fields });
    }
    Ok(PopulationSnvsTsv { rows })
}

/// Metadatos canónicos, preservados como JSON tras verificar su forma cerrada.
#[derive(Clone, Debug, PartialEq)]
pub struct PopulationSnvsMetadata {
    document: Value,
}

impl PopulationSnvsMetadata {
    #[must_use]
    pub const fn document(&self) -> &Value {
        &self.document
    }
}

/// Comprueba la envoltura JSON v1 sin verificar aún hashes o contenido TSV.
pub fn parse_metadata_json(
    bytes: &[u8],
) -> Result<PopulationSnvsMetadata, PopulationSnvsFormatError> {
    let value: Value = serde_json::from_slice(bytes)
        .map_err(|_| PopulationSnvsFormatError::InvalidMetadataJson)?;
    let object = value
        .as_object()
        .ok_or(PopulationSnvsFormatError::MetadataMustBeObject)?;
    closed_object(
        object,
        &[
            "document_kind",
            "schema_version",
            "coordinate_system",
            "reference",
            "matrix",
            "populations",
            "provenance",
        ],
        "",
    )?;
    required_literal(object, "document_kind", "population_snvs_metadata")?;
    required_literal(object, "schema_version", "1")?;
    required_literal(object, "coordinate_system", "base1_closed")?;
    let reference = required_object(object, "reference")?;
    closed_object(reference, &["reference_id", "sha256"], "reference")?;
    required_string(reference, "reference.reference_id")?;
    required_string(reference, "reference.sha256")?;
    let matrix = required_object(object, "matrix")?;
    closed_object(matrix, &["path", "sha256"], "matrix")?;
    required_literal(matrix, "path", "population_snvs.tsv")?;
    required_string(matrix, "matrix.sha256")?;
    let populations = object
        .get("populations")
        .ok_or_else(|| missing("populations"))?
        .as_array()
        .ok_or_else(|| invalid("populations", "una lista"))?;
    for population in populations {
        let population = population
            .as_object()
            .ok_or_else(|| invalid("populations[]", "un objeto"))?;
        closed_object(
            population,
            &[
                "population_id",
                "description",
                "sample_size",
                "allele_count",
                "source",
                "filters",
                "data_version",
            ],
            "populations[]",
        )?;
        required_string(population, "populations[].population_id")?;
    }
    required_provenance(object)?;
    Ok(PopulationSnvsMetadata { document: value })
}

fn closed_object(
    object: &serde_json::Map<String, Value>,
    accepted: &[&str],
    prefix: &str,
) -> Result<(), PopulationSnvsFormatError> {
    for field in object.keys() {
        if !accepted.contains(&field.as_str()) {
            let field = if prefix.is_empty() {
                field.clone()
            } else {
                format!("{prefix}.{field}")
            };
            return Err(PopulationSnvsFormatError::UnexpectedMetadataField { field });
        }
    }
    Ok(())
}

fn required_object<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
) -> Result<&'a serde_json::Map<String, Value>, PopulationSnvsFormatError> {
    object
        .get(field)
        .ok_or_else(|| missing(field))?
        .as_object()
        .ok_or_else(|| invalid(field, "un objeto"))
}

fn required_string(
    object: &serde_json::Map<String, Value>,
    field: &str,
) -> Result<(), PopulationSnvsFormatError> {
    let leaf = field.rsplit('.').next().unwrap_or(field);
    object
        .get(leaf)
        .ok_or_else(|| missing(field))?
        .as_str()
        .ok_or_else(|| invalid(field, "una cadena"))?;
    Ok(())
}

fn required_literal(
    object: &serde_json::Map<String, Value>,
    field: &str,
    expected: &'static str,
) -> Result<(), PopulationSnvsFormatError> {
    object
        .get(field)
        .ok_or_else(|| missing(field))?
        .as_str()
        .filter(|value| *value == expected)
        .ok_or_else(|| PopulationSnvsFormatError::InvalidMetadataLiteral {
            field: field.to_owned(),
            expected,
        })?;
    Ok(())
}

fn required_provenance(
    object: &serde_json::Map<String, Value>,
) -> Result<(), PopulationSnvsFormatError> {
    let provenance = object
        .get("provenance")
        .ok_or_else(|| missing("provenance"))?;
    let kind = provenance
        .as_str()
        .ok_or_else(|| invalid("provenance", "\"supplied\" o \"generated\""))?;
    if !matches!(kind, "supplied" | "generated") {
        return Err(PopulationSnvsFormatError::InvalidMetadataLiteral {
            field: "provenance".to_owned(),
            expected: "\"supplied\" o \"generated\"",
        });
    }
    Ok(())
}

fn missing(field: &str) -> PopulationSnvsFormatError {
    PopulationSnvsFormatError::MissingMetadataField {
        field: field.to_owned(),
    }
}

fn invalid(field: &str, expected: &'static str) -> PopulationSnvsFormatError {
    PopulationSnvsFormatError::InvalidMetadataField {
        field: field.to_owned(),
        expected,
    }
}
