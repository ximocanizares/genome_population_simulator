//! Formatos estructurales del recurso canónico de frecuencias SNV.
//!
//! Este módulo solo fija la envoltura TSV+JSON de F2.2.1. La validación contra
//! FASTA, la semántica de `alt_af`, sumas, IDs canónicos y hashes operativos se
//! incorpora en los pasos consumidores posteriores.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use crate::identity::{PopulationId, SequenceId, VariantId, derive_matrix_variant_id};
use crate::randomness::{
    RandomSamplingError, RootSeed, RuntimeSubseedScope, categorical_f64, derive_subseed,
    fisher_yates, linear_uniform, log_uniform, poisson_count, weighted_without_replacement,
};
use crate::reference::Reference;
use crate::results::ResolvedManifestReference;

use serde_json::{Value, json};
use sha2::{Digest, Sha256};

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

/// Fila activa, canónica y separada de la tabla suministrada original.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedPopulationSnv {
    variant_id: VariantId,
    contig: SequenceId,
    position_1based: u64,
    reference: u8,
    alternate: u8,
    population_id: PopulationId,
    alt_af: f64,
    source_id: String,
}

impl PreparedPopulationSnv {
    #[must_use]
    pub const fn variant_id(&self) -> &VariantId {
        &self.variant_id
    }
    #[must_use]
    pub const fn contig(&self) -> &SequenceId {
        &self.contig
    }
    #[must_use]
    pub const fn position_1based(&self) -> u64 {
        self.position_1based
    }
    #[must_use]
    pub const fn reference(&self) -> u8 {
        self.reference
    }
    #[must_use]
    pub const fn alternate(&self) -> u8 {
        self.alternate
    }
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    #[must_use]
    pub const fn alt_af(&self) -> f64 {
        self.alt_af
    }
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }
}

/// Recurso preparado en memoria; no publica ni altera la tabla suministrada.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedPopulationSnvs {
    rows: Vec<PreparedPopulationSnv>,
    discarded_globally_absent: Vec<VariantId>,
    empty_authorized: bool,
}

/// Procedencia declarada del par canónico publicado.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PopulationSnvsProvenance {
    Supplied,
    Generated,
}

impl PopulationSnvsProvenance {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Supplied => "supplied",
            Self::Generated => "generated",
        }
    }

    fn parse(value: &str) -> Option<Self> {
        match value {
            "supplied" => Some(Self::Supplied),
            "generated" => Some(Self::Generated),
            _ => None,
        }
    }
}

impl PreparedPopulationSnvs {
    #[must_use]
    pub fn rows(&self) -> &[PreparedPopulationSnv] {
        &self.rows
    }
    #[must_use]
    pub fn discarded_globally_absent(&self) -> &[VariantId] {
        &self.discarded_globally_absent
    }
    #[must_use]
    pub const fn empty_authorized(&self) -> bool {
        self.empty_authorized
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum PopulationSnvsPreparationError {
    AuthoritativeFasta {
        detail: String,
    },
    InvalidPopulation {
        row_index: usize,
        value: String,
    },
    InvalidContig {
        row_index: usize,
        value: String,
    },
    InvalidPosition {
        row_index: usize,
        value: String,
    },
    InvalidAllele {
        row_index: usize,
        field: &'static str,
        value: String,
    },
    ReferenceMismatch {
        row_index: usize,
        contig: String,
        position: u64,
        expected: u8,
        found: u8,
    },
    DuplicateAltPopulation {
        contig: String,
        position: u64,
        reference: u8,
        alternate: u8,
        population: String,
    },
    IncompletePopulationCoverage {
        contig: String,
        position: u64,
        reference: u8,
        alternate: u8,
    },
    EmptyActiveMatrixRequiresAuthorization,
    Frequencies(AltFrequencyNormalizationError),
}

impl fmt::Display for PopulationSnvsPreparationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for PopulationSnvsPreparationError {}

/// Valida y prepara una matriz sin modificar sus bytes de entrada.
///
/// La referencia solo puede proceder de `resolve_manifest_reference` (DTC-060).
/// La publicación y los hashes del nuevo par siguen siendo responsabilidad de
/// F2.2.7; este paso devuelve exclusivamente la vista canónica en memoria.
pub fn prepare_supplied_population_snvs(
    table: &PopulationSnvsTsv,
    reference: &ResolvedManifestReference,
    tolerance: f64,
    allow_empty_matrix: bool,
) -> Result<PreparedPopulationSnvs, PopulationSnvsPreparationError> {
    let fasta = load_authoritative_fasta(reference)?;
    prepare_population_snvs(
        table,
        reference.reference_id(),
        &fasta,
        tolerance,
        allow_empty_matrix,
    )
}

/// Variante integrada sobre la referencia ya verificada y materializada por la
/// ejecución. No acepta rutas libres ni altera la tabla de entrada.
pub fn prepare_population_snvs_from_reference(
    table: &PopulationSnvsTsv,
    reference: &Reference,
    tolerance: f64,
    allow_empty_matrix: bool,
) -> Result<PreparedPopulationSnvs, PopulationSnvsPreparationError> {
    let fasta = AuthoritativeFasta {
        sequences: reference
            .sequences()
            .iter()
            .map(|sequence| {
                (
                    sequence.id().as_str().to_owned(),
                    sequence.sequence().as_bytes().to_vec(),
                )
            })
            .collect(),
        order: reference
            .sequences()
            .iter()
            .enumerate()
            .map(|(index, sequence)| (sequence.id().as_str().to_owned(), index))
            .collect(),
    };
    prepare_population_snvs(table, reference.id(), &fasta, tolerance, allow_empty_matrix)
}

fn prepare_population_snvs(
    table: &PopulationSnvsTsv,
    reference_id: &crate::identity::ReferenceId,
    fasta: &AuthoritativeFasta,
    tolerance: f64,
    allow_empty_matrix: bool,
) -> Result<PreparedPopulationSnvs, PopulationSnvsPreparationError> {
    let normalized = normalize_alt_frequencies(table, tolerance)
        .map_err(PopulationSnvsPreparationError::Frequencies)?;
    let mut effective = BTreeMap::new();
    for group in normalized.groups() {
        for value in group.values() {
            effective.insert(value.row_index(), value.effective());
        }
    }
    let populations = table
        .rows()
        .iter()
        .enumerate()
        .map(|(i, row)| {
            PopulationId::parse(row.field(5), "population_snvs").map_err(|_| {
                PopulationSnvsPreparationError::InvalidPopulation {
                    row_index: i,
                    value: row.field(5).to_owned(),
                }
            })
        })
        .collect::<Result<BTreeSet<_>, _>>()?;
    let mut seen = BTreeSet::new();
    let mut alt_populations: BTreeMap<(String, u64, u8, u8), BTreeSet<PopulationId>> =
        BTreeMap::new();
    let mut parsed = Vec::with_capacity(table.rows().len());
    for (row_index, row) in table.rows().iter().enumerate() {
        let contig = SequenceId::parse(row.field(1), "population_snvs").map_err(|_| {
            PopulationSnvsPreparationError::InvalidContig {
                row_index,
                value: row.field(1).to_owned(),
            }
        })?;
        let position = row
            .field(2)
            .parse::<u64>()
            .ok()
            .filter(|value| *value > 0)
            .ok_or_else(|| PopulationSnvsPreparationError::InvalidPosition {
                row_index,
                value: row.field(2).to_owned(),
            })?;
        let reference_base = parse_snv_base(row.field(3), row_index, "ref")?;
        let alternate = parse_snv_base(row.field(4), row_index, "alt")?;
        if reference_base == alternate {
            return Err(PopulationSnvsPreparationError::InvalidAllele {
                row_index,
                field: "alt",
                value: row.field(4).to_owned(),
            });
        }
        let population = PopulationId::parse(row.field(5), "population_snvs").map_err(|_| {
            PopulationSnvsPreparationError::InvalidPopulation {
                row_index,
                value: row.field(5).to_owned(),
            }
        })?;
        let sequence = fasta.sequences.get(contig.as_str()).ok_or_else(|| {
            PopulationSnvsPreparationError::InvalidContig {
                row_index,
                value: contig.as_str().to_owned(),
            }
        })?;
        let offset = usize::try_from(position - 1).ok();
        let expected = offset
            .and_then(|index| sequence.get(index).copied())
            .ok_or_else(|| PopulationSnvsPreparationError::ReferenceMismatch {
                row_index,
                contig: contig.as_str().to_owned(),
                position,
                expected: b'?',
                found: reference_base,
            })?;
        if expected != reference_base {
            return Err(PopulationSnvsPreparationError::ReferenceMismatch {
                row_index,
                contig: contig.as_str().to_owned(),
                position,
                expected,
                found: reference_base,
            });
        }
        let key = (
            contig.as_str().to_owned(),
            position,
            reference_base,
            alternate,
            population.as_str().to_owned(),
        );
        if !seen.insert(key) {
            return Err(PopulationSnvsPreparationError::DuplicateAltPopulation {
                contig: contig.as_str().to_owned(),
                position,
                reference: reference_base,
                alternate,
                population: population.as_str().to_owned(),
            });
        }
        alt_populations
            .entry((
                contig.as_str().to_owned(),
                position,
                reference_base,
                alternate,
            ))
            .or_default()
            .insert(population.clone());
        parsed.push((
            row_index,
            contig,
            position,
            reference_base,
            alternate,
            population,
            row.field(7).to_owned(),
        ));
    }
    if alt_populations
        .values()
        .any(|covered| covered != &populations)
    {
        let Some(((contig, position, reference_base, alternate), _)) = alt_populations
            .iter()
            .find(|(_, covered)| *covered != &populations)
        else {
            return Err(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "invariante de cobertura de población".to_owned(),
            });
        };
        return Err(
            PopulationSnvsPreparationError::IncompletePopulationCoverage {
                contig: contig.clone(),
                position: *position,
                reference: *reference_base,
                alternate: *alternate,
            },
        );
    }
    let mut discarded = Vec::new();
    let inactive: BTreeSet<(String, u64, u8, u8)> = alt_populations
        .iter()
        .filter_map(|((contig, position, reference_base, alternate), _)| {
            let all_zero = parsed
                .iter()
                .filter(|(_, c, p, r, a, _, _)| {
                    c.as_str() == contig && p == position && r == reference_base && a == alternate
                })
                .all(|(index, ..)| effective[index] == 0.0);
            all_zero.then_some((contig.clone(), *position, *reference_base, *alternate))
        })
        .collect::<BTreeSet<_>>();
    let mut rows = Vec::new();
    for (index, contig, position, reference_base, alternate, population, source_id) in parsed {
        let variant_id =
            derive_matrix_variant_id(reference_id, &contig, position, reference_base, alternate)
                .map_err(|error| PopulationSnvsPreparationError::AuthoritativeFasta {
                    detail: error.render_human(),
                })?;
        if inactive.contains(&(
            contig.as_str().to_owned(),
            position,
            reference_base,
            alternate,
        )) {
            if !discarded.contains(&variant_id) {
                discarded.push(variant_id);
            }
            continue;
        }
        rows.push(PreparedPopulationSnv {
            variant_id,
            contig,
            position_1based: position,
            reference: reference_base,
            alternate,
            population_id: population,
            alt_af: effective[&index],
            source_id,
        });
    }
    rows.sort_by(|left, right| {
        fasta
            .ordinal(left.contig().as_str())
            .cmp(&fasta.ordinal(right.contig().as_str()))
            .then(left.position_1based().cmp(&right.position_1based()))
            .then(left.reference().cmp(&right.reference()))
            .then(left.alternate().cmp(&right.alternate()))
            .then(
                left.population_id()
                    .as_str()
                    .cmp(right.population_id().as_str()),
            )
    });
    if rows.is_empty() && !allow_empty_matrix {
        return Err(PopulationSnvsPreparationError::EmptyActiveMatrixRequiresAuthorization);
    }
    let empty_authorized = rows.is_empty();
    Ok(PreparedPopulationSnvs {
        rows,
        discarded_globally_absent: discarded,
        empty_authorized,
    })
}

/// Serialización canónica, estable y sin formato numérico dependiente de la
/// plataforma del recurso TSV F2.2.7.
#[must_use]
pub fn canonical_tsv_bytes(prepared: &PreparedPopulationSnvs) -> Vec<u8> {
    let mut text = TSV_HEADER.join("\t");
    text.push('\n');
    for row in prepared.rows() {
        let _ = writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.variant_id().as_str(),
            row.contig().as_str(),
            row.position_1based(),
            char::from(row.reference()),
            char::from(row.alternate()),
            row.population_id().as_str(),
            row.alt_af(),
            row.source_id(),
        );
    }
    text.into_bytes()
}

/// Construye el JSON cerrado del par canónico después de calcular el SHA-256
/// del TSV. Las poblaciones se derivan de las filas activas y quedan ordenadas.
#[must_use]
pub fn canonical_metadata_document(
    prepared: &PreparedPopulationSnvs,
    reference_id: &str,
    reference_sha256: &str,
    matrix_sha256: &str,
    provenance: PopulationSnvsProvenance,
) -> Value {
    let populations = prepared
        .rows()
        .iter()
        .map(|row| row.population_id().as_str())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(|population_id| json!({"population_id": population_id}))
        .collect::<Vec<_>>();
    json!({
        "document_kind": "population_snvs_metadata",
        "schema_version": "1",
        "coordinate_system": "base1_closed",
        "reference": {"reference_id": reference_id, "sha256": reference_sha256},
        "matrix": {"path": "population_snvs.tsv", "sha256": matrix_sha256},
        "populations": populations,
        "provenance": provenance.as_str(),
    })
}

/// Valida un par reutilizado y devuelve su vista semántica. El TSV se exige en
/// bytes canónicos para evitar que una entrada equivalente reintroduzca IDs,
/// orden o representación numérica no canónicos.
pub fn validate_reusable_population_snvs(
    metadata_path: &Path,
    reference: &Reference,
    reference_sha256: &str,
    tolerance: f64,
    allow_empty_matrix: bool,
) -> Result<
    (PreparedPopulationSnvs, PopulationSnvsProvenance, String),
    PopulationSnvsPreparationError,
> {
    validate_reusable_population_snvs_with_hooks(
        metadata_path,
        reference,
        reference_sha256,
        tolerance,
        allow_empty_matrix,
        || {},
        || {},
    )
}

/// Igual que [`validate_reusable_population_snvs`], con puntos de prueba entre
/// las dos lecturas. Ninguna comprobación posterior vuelve a abrir ninguno de
/// los dos recursos: hashes y semántica proceden de estas instantáneas.
pub fn validate_reusable_population_snvs_with_hooks(
    metadata_path: &Path,
    reference: &Reference,
    reference_sha256: &str,
    tolerance: f64,
    allow_empty_matrix: bool,
    after_metadata_snapshot: impl FnOnce(),
    after_matrix_snapshot: impl FnOnce(),
) -> Result<
    (PreparedPopulationSnvs, PopulationSnvsProvenance, String),
    PopulationSnvsPreparationError,
> {
    let metadata_bytes = fs::read(metadata_path).map_err(|error| {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: format!("no se puede leer metadata_path: {error}"),
        }
    })?;
    after_metadata_snapshot();
    let matrix_path = metadata_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("population_snvs.tsv");
    let matrix_bytes = fs::read(&matrix_path).map_err(|error| {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: error.to_string(),
        }
    })?;
    after_matrix_snapshot();
    validate_reusable_population_snvs_bytes(
        &metadata_bytes,
        &matrix_bytes,
        reference,
        reference_sha256,
        tolerance,
        allow_empty_matrix,
    )
}

/// Valida las instantáneas exactas que forman un par reutilizable. Esta
/// frontera también la reutiliza la validación física previa a publicar.
pub fn validate_reusable_population_snvs_bytes(
    metadata_bytes: &[u8],
    matrix_bytes: &[u8],
    reference: &Reference,
    reference_sha256: &str,
    tolerance: f64,
    allow_empty_matrix: bool,
) -> Result<
    (PreparedPopulationSnvs, PopulationSnvsProvenance, String),
    PopulationSnvsPreparationError,
> {
    let metadata = parse_metadata_json(metadata_bytes).map_err(|error| {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: error.to_string(),
        }
    })?;
    let document = metadata.document();
    let object =
        document
            .as_object()
            .ok_or_else(|| PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata no es un objeto".to_owned(),
            })?;
    let reference_object = object
        .get("reference")
        .and_then(Value::as_object)
        .ok_or_else(|| PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "metadata.reference inválida".to_owned(),
        })?;
    let declared_reference = reference_object.get("reference_id").and_then(Value::as_str);
    let declared_sha = reference_object.get("sha256").and_then(Value::as_str);
    if declared_reference != Some(reference.id().as_str()) || declared_sha != Some(reference_sha256)
    {
        return Err(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "la matriz reutilizada no está ligada a la referencia autoritativa".to_owned(),
        });
    }
    let matrix = object
        .get("matrix")
        .and_then(Value::as_object)
        .ok_or_else(|| PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "metadata.matrix inválida".to_owned(),
        })?;
    let declared_digest = matrix
        .get("sha256")
        .and_then(Value::as_str)
        .ok_or_else(|| PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "metadata.matrix.sha256 inválido".to_owned(),
        })?;
    if sha256_bytes(matrix_bytes) != declared_digest {
        return Err(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "SHA-256 del TSV reutilizado no coincide con metadata".to_owned(),
        });
    }
    let table = parse_canonical_tsv(matrix_bytes).map_err(|error| {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: error.to_string(),
        }
    })?;
    let prepared =
        prepare_population_snvs_from_reference(&table, reference, tolerance, allow_empty_matrix)?;
    if canonical_tsv_bytes(&prepared) != matrix_bytes {
        return Err(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "el TSV reutilizado no está en orden y representación canónicos".to_owned(),
        });
    }
    let populations = object
        .get("populations")
        .and_then(Value::as_array)
        .ok_or_else(|| PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "metadata.populations inválida".to_owned(),
        })?;
    let expected = prepared
        .rows()
        .iter()
        .map(|row| row.population_id().as_str())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let observed = populations
        .iter()
        .filter_map(|value| value.get("population_id").and_then(Value::as_str))
        .collect::<Vec<_>>();
    if observed != expected {
        return Err(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "metadata.populations no coincide con el TSV canónico".to_owned(),
        });
    }
    let provenance = object
        .get("provenance")
        .and_then(Value::as_str)
        .and_then(PopulationSnvsProvenance::parse)
        .ok_or_else(|| PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "metadata.provenance inválida".to_owned(),
        })?;
    Ok((prepared, provenance, sha256_bytes(metadata_bytes)))
}

fn sha256_bytes(bytes: &[u8]) -> String {
    let mut digest = Sha256::new();
    digest.update(bytes);
    format!("{:x}", digest.finalize())
}

fn parse_snv_base(
    value: &str,
    row_index: usize,
    field: &'static str,
) -> Result<u8, PopulationSnvsPreparationError> {
    let bytes = value.as_bytes();
    if bytes.len() == 1 && matches!(bytes[0], b'A' | b'C' | b'G' | b'T') {
        Ok(bytes[0])
    } else {
        Err(PopulationSnvsPreparationError::InvalidAllele {
            row_index,
            field,
            value: value.to_owned(),
        })
    }
}

struct AuthoritativeFasta {
    sequences: BTreeMap<String, Vec<u8>>,
    order: BTreeMap<String, usize>,
}
impl AuthoritativeFasta {
    fn ordinal(&self, id: &str) -> usize {
        self.order.get(id).copied().unwrap_or(usize::MAX)
    }
}

fn load_authoritative_fasta(
    reference: &ResolvedManifestReference,
) -> Result<AuthoritativeFasta, PopulationSnvsPreparationError> {
    let bytes = fs::read(reference.fasta_path()).map_err(|error| {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: error.to_string(),
        }
    })?;
    let mut records = BTreeMap::new();
    let mut order = BTreeMap::new();
    let mut name = None::<String>;
    let mut bases = Vec::new();
    for raw in bytes.split(|byte| *byte == b'\n') {
        let line = raw.strip_suffix(b"\r").unwrap_or(raw);
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(previous) = name.replace(
                header
                    .split(u8::is_ascii_whitespace)
                    .next()
                    .unwrap_or_default()
                    .iter()
                    .copied()
                    .map(char::from)
                    .collect(),
            ) {
                if previous.is_empty()
                    || records
                        .insert(previous.clone(), std::mem::take(&mut bases))
                        .is_some()
                {
                    return Err(PopulationSnvsPreparationError::AuthoritativeFasta {
                        detail: "cabecera FASTA duplicada o vacía".to_owned(),
                    });
                }
                order.insert(previous, order.len());
            }
        } else {
            bases.extend(line.iter().map(|base| base.to_ascii_uppercase()));
        }
    }
    let Some(last) = name else {
        return Err(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "FASTA vacío".to_owned(),
        });
    };
    if last.is_empty() || records.insert(last.clone(), bases).is_some() {
        return Err(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "cabecera FASTA duplicada o vacía".to_owned(),
        });
    }
    order.insert(last, order.len());
    Ok(AuthoritativeFasta {
        sequences: records,
        order,
    })
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

/// Analiza la tabla mínima suministrada (seis columnas y `source_id` opcional)
/// sin confundirla con el recurso canónico publicado de ocho columnas.
pub fn parse_supplied_tsv(bytes: &[u8]) -> Result<PopulationSnvsTsv, PopulationSnvsFormatError> {
    let text = std::str::from_utf8(bytes).map_err(|_| PopulationSnvsFormatError::InvalidUtf8)?;
    if text.starts_with("##fileformat=VCF") || text.starts_with("#CHROM\t") {
        return Err(PopulationSnvsFormatError::VcfInput);
    }
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or(PopulationSnvsFormatError::MissingHeader)?;
    let accepted = [
        "contig\tposition\tref\talt\tpopulation_id\talt_af",
        "contig\tposition\tref\talt\tpopulation_id\talt_af\tsource_id",
    ];
    if !accepted.contains(&header) {
        return Err(PopulationSnvsFormatError::IncorrectHeader {
            found: header.to_owned(),
        });
    }
    let has_source = header.ends_with("\tsource_id");
    let mut rows = Vec::new();
    for (index, line) in lines.enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != if has_source { 7 } else { 6 } {
            return Err(PopulationSnvsFormatError::IncorrectColumnCount {
                line: index + 2,
                found: fields.len(),
            });
        }
        rows.push(PopulationSnvRow {
            fields: [
                String::new(),
                fields[0].to_owned(),
                fields[1].to_owned(),
                fields[2].to_owned(),
                fields[3].to_owned(),
                fields[4].to_owned(),
                fields[5].to_owned(),
                if has_source {
                    fields[6].to_owned()
                } else {
                    String::new()
                },
            ],
        });
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

/// Recuento de candidatos de una matriz generada. `Poisson` conserva la media,
/// mientras que `Exact` no abre el ámbito RNG de recuento (DTC-067).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GeneratedLocusCount {
    Exact(u64),
    PoissonMean(f64),
    PoissonDensityPerBase(f64),
}

/// Región base 1 cerrada para el espacio ponderado de la matriz.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneratedWeightedRegion {
    pub contig: String,
    pub start_1based: u64,
    pub end_1based: u64,
    pub weight: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum GeneratedSpatialMode {
    Uniform,
    WeightedRegions {
        outside_weight: f64,
        regions: Vec<GeneratedWeightedRegion>,
    },
}

/// Método con que un componente de espectro produce una frecuencia ALT.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FrequencySampling {
    LinearUniform,
    LogUniform,
}

/// Componente ordenado de un espectro de frecuencias ALT.
#[derive(Clone, Debug, PartialEq)]
pub enum FrequencySpectrumComponent {
    Fixed {
        value: f64,
        weight: f64,
    },
    Interval {
        min: f64,
        max: f64,
        weight: f64,
        sampling: FrequencySampling,
    },
}

impl FrequencySpectrumComponent {
    fn weight(&self) -> f64 {
        match self {
            Self::Fixed { weight, .. } | Self::Interval { weight, .. } => *weight,
        }
    }
}

/// Política de asociación de componentes a los loci de una población.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FrequencyAssignmentMode {
    Stochastic,
    Exact,
}

/// Espectro explícito de frecuencias ALT, en orden estable de componentes.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneratedFrequencySpectrum {
    pub mode: FrequencyAssignmentMode,
    pub components: Vec<FrequencySpectrumComponent>,
}

/// Perfil de pesos para los doce cambios dirigidos de SNV.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GeneratedSubstitutionProfile {
    UniformSnv12V1,
    HumanHarris2017Snv12ApproxV1,
}

/// Compatibilidad F2.2.5: cuando `spectrum` es `None`, `min_alt_af` y
/// `max_alt_af` representan el único intervalo lineal estocástico con peso 1.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneratedPopulationFrequency {
    pub population_id: String,
    pub min_alt_af: f64,
    pub max_alt_af: f64,
    pub spectrum: Option<GeneratedFrequencySpectrum>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PopulationSnvsGenerationConfig {
    pub count: GeneratedLocusCount,
    pub spatial: GeneratedSpatialMode,
    pub substitution_profile: GeneratedSubstitutionProfile,
    pub populations: Vec<GeneratedPopulationFrequency>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PopulationSnvsGenerationError {
    InvalidCount,
    InvalidPopulation(String),
    DuplicatePopulation(String),
    InvalidFrequencyInterval { population_id: String },
    InvalidFrequencySpectrum { population_id: String },
    InvalidRegion { contig: String },
    OverlappingRegions { contig: String },
    InvalidWeight,
    EmptyEligibleSpace,
    InsufficientEligibleSpace { requested: u64, available: u64 },
    Random(RandomSamplingError),
    Reference(String),
}

impl fmt::Display for PopulationSnvsGenerationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for PopulationSnvsGenerationError {}

/// Vista mínima de una referencia ya validada: identidad, orden FASTA y bases.
/// No reinterpreta rutas ni crea una identidad alternativa.
struct GenerationReference {
    id: crate::identity::ReferenceId,
    sequences: Vec<(SequenceId, Vec<u8>)>,
}

/// Genera la vista activa de una matriz nueva desde una referencia F1 resuelta.
pub fn generate_population_snvs(
    config: &PopulationSnvsGenerationConfig,
    reference: &ResolvedManifestReference,
    root_seed: RootSeed,
) -> Result<PreparedPopulationSnvs, PopulationSnvsGenerationError> {
    let fasta = load_authoritative_fasta(reference)
        .map_err(|error| PopulationSnvsGenerationError::Reference(error.to_string()))?;
    let sequence_order = fasta.order.clone();
    let mut ordered = fasta.sequences.into_iter().collect::<Vec<_>>();
    ordered.sort_by_key(|(id, _)| sequence_order.get(id).copied());
    generate_population_snvs_view(
        config,
        GenerationReference {
            id: reference.reference_id().clone(),
            sequences: ordered
                .into_iter()
                .map(|(id, bases)| {
                    Ok((
                        SequenceId::parse(id, "generated_population_snvs").map_err(|_| {
                            PopulationSnvsGenerationError::Reference("invalid sequence".into())
                        })?,
                        bases,
                    ))
                })
                .collect::<Result<_, _>>()?,
        },
        root_seed,
    )
}

/// Misma operación en el recorrido directo, sobre la referencia ya validada en
/// memoria. No publica la matriz ni sustituye el par canónico de F2.2.7.
pub fn generate_population_snvs_from_reference(
    config: &PopulationSnvsGenerationConfig,
    reference: &Reference,
    root_seed: RootSeed,
) -> Result<PreparedPopulationSnvs, PopulationSnvsGenerationError> {
    generate_population_snvs_view(
        config,
        GenerationReference {
            id: reference.id().clone(),
            sequences: reference
                .sequences()
                .iter()
                .map(|sequence| {
                    (
                        sequence.id().clone(),
                        sequence.sequence().as_bytes().to_vec(),
                    )
                })
                .collect(),
        },
        root_seed,
    )
}

fn generate_population_snvs_view(
    config: &PopulationSnvsGenerationConfig,
    reference: GenerationReference,
    root_seed: RootSeed,
) -> Result<PreparedPopulationSnvs, PopulationSnvsGenerationError> {
    validate_generation_config(config)?;
    validate_generation_regions(config, &reference.sequences)?;
    let mut candidates = Vec::<(SequenceId, u64, u8)>::new();
    let mut weights = Vec::new();
    let order = reference
        .sequences
        .iter()
        .enumerate()
        .map(|(index, (id, _))| (id.as_str().to_owned(), index))
        .collect::<BTreeMap<_, _>>();
    for (id, bases) in &reference.sequences {
        let contig = id.as_str();
        for (offset, base) in bases.iter().copied().enumerate() {
            if !matches!(base, b'A' | b'C' | b'G' | b'T') {
                continue;
            }
            let position =
                u64::try_from(offset).map_err(|_| PopulationSnvsGenerationError::InvalidCount)? + 1;
            let weight = candidate_weight(&config.spatial, contig, position)?;
            candidates.push((id.clone(), position, base));
            weights.push(weight);
        }
    }
    if weights.iter().all(|weight| *weight == 0.0) {
        return Err(PopulationSnvsGenerationError::EmptyEligibleSpace);
    }
    let requested = match config.count {
        GeneratedLocusCount::Exact(value) => value,
        GeneratedLocusCount::PoissonMean(mean) => {
            sample_candidate_count(root_seed, &reference.id, mean)?
        }
        GeneratedLocusCount::PoissonDensityPerBase(density) => {
            let eligible = weights.iter().filter(|weight| **weight > 0.0).count() as f64;
            sample_candidate_count(root_seed, &reference.id, density * eligible)?
        }
    };
    let available = weights.iter().filter(|weight| **weight > 0.0).count() as u64;
    if requested > available {
        return Err(PopulationSnvsGenerationError::InsufficientEligibleSpace {
            requested,
            available,
        });
    }
    let requested_usize =
        usize::try_from(requested).map_err(|_| PopulationSnvsGenerationError::InvalidCount)?;
    let mut selector = derive_subseed(
        root_seed,
        RuntimeSubseedScope::MatrixLocusSelection {
            scenario_id: "direct",
            reference_id: reference.id.as_str(),
        },
    )
    .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?
    .local_rng();
    let selected = weighted_without_replacement(&mut selector, &weights, requested_usize)
        .map_err(PopulationSnvsGenerationError::Random)?;
    let mut loci = Vec::with_capacity(selected.len());
    for index in selected {
        let (contig, position, reference_base) = &candidates[index];
        let mut alt_rng = derive_subseed(
            root_seed,
            RuntimeSubseedScope::MatrixAlternateAllele {
                scenario_id: "direct",
                reference_id: reference.id.as_str(),
                sequence_id: contig.as_str(),
                position_1based: *position,
                reference_base: std::str::from_utf8(&[*reference_base]).unwrap_or("A"),
            },
        )
        .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?
        .local_rng();
        let alternatives: Vec<u8> = (*b"ACGT")
            .into_iter()
            .filter(|base| base != reference_base)
            .collect();
        let alternate = alternatives[categorical_f64(
            &mut alt_rng,
            &substitution_weights(config.substitution_profile, *reference_base, &alternatives),
        )
        .map_err(PopulationSnvsGenerationError::Random)?];
        let variant_id =
            derive_matrix_variant_id(&reference.id, contig, *position, *reference_base, alternate)
                .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?;
        loci.push((
            contig.clone(),
            *position,
            *reference_base,
            alternate,
            variant_id,
        ));
    }
    loci.sort_by(|left, right| {
        order
            .get(left.0.as_str())
            .cmp(&order.get(right.0.as_str()))
            .then(left.1.cmp(&right.1))
    });
    let mut rows = Vec::with_capacity(loci.len() * config.populations.len());
    for population in &config.populations {
        let population_id = PopulationId::parse(
            &population.population_id,
            "generated_population_snvs",
        )
        .map_err(|_| {
            PopulationSnvsGenerationError::InvalidPopulation(population.population_id.clone())
        })?;
        let spectrum = population_spectrum(population);
        let assignments = assign_frequency_components(
            &spectrum,
            &loci,
            root_seed,
            &reference.id,
            population_id.as_str(),
        )?;
        for ((contig, position, reference_base, alternate, variant_id), component_index) in
            loci.iter().zip(assignments)
        {
            let component = &spectrum.components[component_index];
            let alt_af = sample_frequency_component(
                component,
                root_seed,
                &reference.id,
                population_id.as_str(),
                variant_id.as_str(),
            )?;
            rows.push(PreparedPopulationSnv {
                variant_id: variant_id.clone(),
                contig: contig.clone(),
                position_1based: *position,
                reference: *reference_base,
                alternate: *alternate,
                population_id: population_id.clone(),
                alt_af,
                source_id: "generated".to_owned(),
            });
        }
    }
    rows.sort_by(|left, right| {
        order
            .get(left.contig().as_str())
            .cmp(&order.get(right.contig().as_str()))
            .then(left.position_1based().cmp(&right.position_1based()))
            .then(left.alternate().cmp(&right.alternate()))
            .then(
                left.population_id()
                    .as_str()
                    .cmp(right.population_id().as_str()),
            )
    });
    Ok(PreparedPopulationSnvs {
        rows,
        discarded_globally_absent: Vec::new(),
        empty_authorized: false,
    })
}

fn substitution_weights(
    profile: GeneratedSubstitutionProfile,
    reference: u8,
    alternatives: &[u8],
) -> Vec<f64> {
    let weights = match profile {
        GeneratedSubstitutionProfile::UniformSnv12V1 => [1.0; 3],
        GeneratedSubstitutionProfile::HumanHarris2017Snv12ApproxV1 => {
            harris_substitution_weights(reference, alternatives)
        }
    };
    weights.to_vec()
}

fn harris_substitution_weights(reference: u8, alternatives: &[u8]) -> [f64; 3] {
    const PROFILE: &str = include_str!("../resources/f2_2_6/human_harris2017_snv12_approx_v1.tsv");
    let mut weights = BTreeMap::new();
    for line in PROFILE.lines().skip(1) {
        let mut fields = line.split('\t');
        let Some(change) = fields.next() else {
            continue;
        };
        let Some(weight) = fields.next().and_then(|value| value.parse::<f64>().ok()) else {
            continue;
        };
        weights.insert(change, weight);
    }
    let mut result = [0.0; 3];
    for (index, alternate) in alternatives.iter().copied().enumerate() {
        let change = format!("{}>{}", char::from(reference), char::from(alternate));
        result[index] = weights.get(change.as_str()).copied().unwrap_or(0.0);
    }
    result
}

fn population_spectrum(population: &GeneratedPopulationFrequency) -> GeneratedFrequencySpectrum {
    population
        .spectrum
        .clone()
        .unwrap_or(GeneratedFrequencySpectrum {
            mode: FrequencyAssignmentMode::Stochastic,
            components: vec![FrequencySpectrumComponent::Interval {
                min: population.min_alt_af,
                max: population.max_alt_af,
                weight: 1.0,
                sampling: FrequencySampling::LinearUniform,
            }],
        })
}

fn assign_frequency_components(
    spectrum: &GeneratedFrequencySpectrum,
    loci: &[(SequenceId, u64, u8, u8, VariantId)],
    root_seed: RootSeed,
    reference_id: &crate::identity::ReferenceId,
    population_id: &str,
) -> Result<Vec<usize>, PopulationSnvsGenerationError> {
    let weights = spectrum
        .components
        .iter()
        .map(FrequencySpectrumComponent::weight)
        .collect::<Vec<_>>();
    match spectrum.mode {
        FrequencyAssignmentMode::Stochastic => loci
            .iter()
            .map(|(_, _, _, _, variant_id)| {
                let mut rng = derive_subseed(
                    root_seed,
                    RuntimeSubseedScope::MatrixFrequencyComponent {
                        scenario_id: "direct",
                        reference_id: reference_id.as_str(),
                        population_id,
                        variant_id: variant_id.as_str(),
                    },
                )
                .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?
                .local_rng();
                categorical_f64(&mut rng, &weights).map_err(PopulationSnvsGenerationError::Random)
            })
            .collect(),
        FrequencyAssignmentMode::Exact => {
            let counts = largest_remainder_counts(loci.len(), &weights)?;
            let mut permuted = (0..loci.len()).collect::<Vec<_>>();
            permuted.sort_by(|left, right| loci[*left].4.as_str().cmp(loci[*right].4.as_str()));
            let mut rng = derive_subseed(
                root_seed,
                RuntimeSubseedScope::MatrixExactAssignment {
                    scenario_id: "direct",
                    reference_id: reference_id.as_str(),
                    population_id,
                },
            )
            .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?
            .local_rng();
            fisher_yates(&mut rng, &mut permuted).map_err(PopulationSnvsGenerationError::Random)?;
            let mut assigned = vec![0; loci.len()];
            let mut offset = 0;
            for (component, count) in counts.into_iter().enumerate() {
                for locus in &permuted[offset..offset + count] {
                    assigned[*locus] = component;
                }
                offset += count;
            }
            Ok(assigned)
        }
    }
}

fn largest_remainder_counts(
    count: usize,
    weights: &[f64],
) -> Result<Vec<usize>, PopulationSnvsGenerationError> {
    let total = neumaier(weights.iter().copied());
    if !total.is_finite() || total <= 0.0 {
        return Err(PopulationSnvsGenerationError::InvalidFrequencySpectrum {
            population_id: "unknown".to_owned(),
        });
    }
    let mut quotas = weights
        .iter()
        .enumerate()
        .map(|(index, weight)| {
            let quota = (count as f64) * *weight / total;
            let floor = quota.floor() as usize;
            (index, floor, quota - floor as f64)
        })
        .collect::<Vec<_>>();
    let assigned = quotas.iter().map(|(_, floor, _)| *floor).sum::<usize>();
    quotas.sort_by(|left, right| {
        right
            .2
            .partial_cmp(&left.2)
            .unwrap_or(Ordering::Equal)
            .then(left.0.cmp(&right.0))
    });
    for (_, floor, _) in quotas.iter_mut().take(count - assigned) {
        *floor += 1;
    }
    quotas.sort_by_key(|(index, _, _)| *index);
    Ok(quotas.into_iter().map(|(_, floor, _)| floor).collect())
}

fn sample_frequency_component(
    component: &FrequencySpectrumComponent,
    root_seed: RootSeed,
    reference_id: &crate::identity::ReferenceId,
    population_id: &str,
    variant_id: &str,
) -> Result<f64, PopulationSnvsGenerationError> {
    match component {
        FrequencySpectrumComponent::Fixed { value, .. } => Ok(*value),
        FrequencySpectrumComponent::Interval {
            min, max, sampling, ..
        } => {
            let mut rng = derive_subseed(
                root_seed,
                RuntimeSubseedScope::MatrixFrequencyValue {
                    scenario_id: "direct",
                    reference_id: reference_id.as_str(),
                    population_id,
                    variant_id,
                },
            )
            .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?
            .local_rng();
            match sampling {
                FrequencySampling::LinearUniform => linear_uniform(&mut rng, *min, *max),
                FrequencySampling::LogUniform => log_uniform(&mut rng, *min, *max),
            }
            .map_err(PopulationSnvsGenerationError::Random)
        }
    }
}

fn sample_candidate_count(
    root: RootSeed,
    reference: &crate::identity::ReferenceId,
    mean: f64,
) -> Result<u64, PopulationSnvsGenerationError> {
    let mut rng = derive_subseed(
        root,
        RuntimeSubseedScope::MatrixCandidateCount {
            scenario_id: "direct",
            reference_id: reference.as_str(),
        },
    )
    .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?
    .local_rng();
    poisson_count(&mut rng, mean).map_err(PopulationSnvsGenerationError::Random)
}

fn candidate_weight(
    spatial: &GeneratedSpatialMode,
    contig: &str,
    position: u64,
) -> Result<f64, PopulationSnvsGenerationError> {
    match spatial {
        GeneratedSpatialMode::Uniform => Ok(1.0),
        GeneratedSpatialMode::WeightedRegions {
            outside_weight,
            regions,
        } => {
            let mut weight = *outside_weight;
            for region in regions {
                if region.contig == contig
                    && (region.start_1based..=region.end_1based).contains(&position)
                {
                    weight = region.weight;
                    break;
                }
            }
            Ok(weight)
        }
    }
}

pub(crate) fn validate_generation_config(
    config: &PopulationSnvsGenerationConfig,
) -> Result<(), PopulationSnvsGenerationError> {
    match config.count {
        GeneratedLocusCount::Exact(_) => {}
        GeneratedLocusCount::PoissonMean(v) | GeneratedLocusCount::PoissonDensityPerBase(v)
            if v.is_finite() && v >= 0.0 => {}
        _ => return Err(PopulationSnvsGenerationError::InvalidCount),
    }
    if config.populations.is_empty() {
        return Err(PopulationSnvsGenerationError::InvalidPopulation(
            "missing".to_owned(),
        ));
    }
    let mut populations = BTreeSet::new();
    for population in &config.populations {
        PopulationId::parse(&population.population_id, "generated_population_snvs").map_err(
            |_| PopulationSnvsGenerationError::InvalidPopulation(population.population_id.clone()),
        )?;
        if !populations.insert(&population.population_id) {
            return Err(PopulationSnvsGenerationError::DuplicatePopulation(
                population.population_id.clone(),
            ));
        }
        if let Some(spectrum) = &population.spectrum {
            validate_frequency_spectrum(spectrum, &population.population_id)?;
        } else if !population.min_alt_af.is_finite()
            || !population.max_alt_af.is_finite()
            || population.min_alt_af < 0.0
            || population.min_alt_af >= population.max_alt_af
            || population.max_alt_af > 1.0
        {
            return Err(PopulationSnvsGenerationError::InvalidFrequencyInterval {
                population_id: population.population_id.clone(),
            });
        }
    }
    if let GeneratedSpatialMode::WeightedRegions {
        outside_weight,
        regions,
    } = &config.spatial
    {
        if !outside_weight.is_finite() || *outside_weight < 0.0 {
            return Err(PopulationSnvsGenerationError::InvalidWeight);
        }
        let mut sorted = regions.clone();
        sorted.sort_by(|a, b| {
            a.contig
                .cmp(&b.contig)
                .then(a.start_1based.cmp(&b.start_1based))
        });
        for region in &sorted {
            if region.start_1based == 0 || region.start_1based > region.end_1based {
                return Err(PopulationSnvsGenerationError::InvalidRegion {
                    contig: region.contig.clone(),
                });
            }
            if !region.weight.is_finite() || region.weight < 0.0 {
                return Err(PopulationSnvsGenerationError::InvalidWeight);
            }
        }
        for pair in sorted.windows(2) {
            if pair[0].contig == pair[1].contig && pair[1].start_1based <= pair[0].end_1based {
                return Err(PopulationSnvsGenerationError::OverlappingRegions {
                    contig: pair[0].contig.clone(),
                });
            }
        }
    }
    Ok(())
}

fn validate_frequency_spectrum(
    spectrum: &GeneratedFrequencySpectrum,
    population_id: &str,
) -> Result<(), PopulationSnvsGenerationError> {
    if spectrum.components.is_empty() {
        return Err(PopulationSnvsGenerationError::InvalidFrequencySpectrum {
            population_id: population_id.to_owned(),
        });
    }
    let mut intervals = Vec::new();
    let mut fixed = Vec::new();
    let mut positive_weight = false;
    for component in &spectrum.components {
        match component {
            FrequencySpectrumComponent::Fixed { value, weight } => {
                if !value.is_finite()
                    || !(0.0..=1.0).contains(value)
                    || !weight.is_finite()
                    || *weight < 0.0
                {
                    return Err(PopulationSnvsGenerationError::InvalidFrequencySpectrum {
                        population_id: population_id.to_owned(),
                    });
                }
                fixed.push(*value);
                positive_weight |= *weight > 0.0;
            }
            FrequencySpectrumComponent::Interval {
                min,
                max,
                weight,
                sampling,
            } => {
                if !min.is_finite()
                    || !max.is_finite()
                    || *min < 0.0
                    || *min >= *max
                    || *max > 1.0
                    || !weight.is_finite()
                    || *weight < 0.0
                    || (*sampling == FrequencySampling::LogUniform && *min <= 0.0)
                {
                    return Err(PopulationSnvsGenerationError::InvalidFrequencySpectrum {
                        population_id: population_id.to_owned(),
                    });
                }
                intervals.push((*min, *max));
                positive_weight |= *weight > 0.0;
            }
        }
    }
    intervals.sort_by(|left, right| left.0.partial_cmp(&right.0).unwrap_or(Ordering::Equal));
    fixed.sort_by(|left, right| left.partial_cmp(right).unwrap_or(Ordering::Equal));
    if !positive_weight
        || intervals.windows(2).any(|pair| pair[1].0 < pair[0].1)
        || fixed.windows(2).any(|pair| pair[0] == pair[1])
        || fixed.iter().any(|value| {
            intervals
                .iter()
                .any(|(min, max)| *min <= *value && *value < *max)
        })
    {
        return Err(PopulationSnvsGenerationError::InvalidFrequencySpectrum {
            population_id: population_id.to_owned(),
        });
    }
    Ok(())
}

pub(crate) fn validate_generation_regions(
    config: &PopulationSnvsGenerationConfig,
    sequences: &[(SequenceId, Vec<u8>)],
) -> Result<(), PopulationSnvsGenerationError> {
    let lengths = sequences
        .iter()
        .map(|(id, bases)| (id.as_str(), u64::try_from(bases.len()).unwrap_or(u64::MAX)))
        .collect::<BTreeMap<_, _>>();
    if let GeneratedSpatialMode::WeightedRegions { regions, .. } = &config.spatial {
        for region in regions {
            let Some(length) = lengths.get(region.contig.as_str()) else {
                return Err(PopulationSnvsGenerationError::InvalidRegion {
                    contig: region.contig.clone(),
                });
            };
            if region.start_1based == 0 || region.end_1based > *length {
                return Err(PopulationSnvsGenerationError::InvalidRegion {
                    contig: region.contig.clone(),
                });
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod profile_resource_tests {
    use super::{GeneratedSubstitutionProfile, substitution_weights};

    #[test]
    fn harris_profile_weights_read_by_rust_match_the_independent_versioned_resource() {
        let expected = include_str!("../resources/f2_2_6/human_harris2017_snv12_approx_v1.tsv")
            .lines()
            .skip(1)
            .map(|line| {
                let (change, weight) = line
                    .split_once('\t')
                    .expect("every versioned Harris row has its weight");
                (
                    change,
                    weight.parse::<f64>().expect("Harris weight is finite"),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(expected.len(), 12);
        for reference in *b"ACGT" {
            let alternatives = (*b"ACGT")
                .into_iter()
                .filter(|alternate| *alternate != reference)
                .collect::<Vec<_>>();
            for (alternate, observed) in alternatives.iter().copied().zip(substitution_weights(
                GeneratedSubstitutionProfile::HumanHarris2017Snv12ApproxV1,
                reference,
                &alternatives,
            )) {
                let change = format!("{}>{}", char::from(reference), char::from(alternate));
                let (_, expected_weight) = expected
                    .iter()
                    .find(|(expected_change, _)| *expected_change == change)
                    .expect("the twelve directed substitutions are exhaustive");
                assert_eq!(observed.to_bits(), (*expected_weight).to_bits(), "{change}");
            }
        }
    }
}
