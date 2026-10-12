//! Formatos estructurales del recurso canónico de frecuencias SNV.
//!
//! Este módulo solo fija la envoltura TSV+JSON de F2.2.1. La validación contra
//! FASTA, la semántica de `alt_af`, sumas, IDs canónicos y hashes operativos se
//! incorpora en los pasos consumidores posteriores.

use crate::guard_trace::rejected;
use std::borrow::Cow;
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};
use std::fmt;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use crate::identity::{PopulationId, ReferenceId, SequenceId, VariantId, derive_matrix_variant_id};
use crate::randomness::{
    RandomSamplingError, RootSeed, RuntimeSubseedScope, categorical_f64, derive_subseed,
    fisher_yates, linear_uniform, log_uniform, open_unit_f64, poisson_count,
};
use crate::reference::{Reference, ReferenceError};
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

/// Huella SHA-256 de la vista canónica `id`+secuencias+bases de una referencia.
///
/// No es el SHA-256 físico del FASTA: ese hash autentica un artefacto publicado
/// y se conserva separado en sus metadatos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CanonicalReferenceFingerprint(String);

impl CanonicalReferenceFingerprint {
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Identidad tipada que enlaza una matriz preparada con la referencia exacta.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PreparedReferenceIdentity {
    reference_id: ReferenceId,
    canonical_fingerprint: CanonicalReferenceFingerprint,
}

impl PreparedReferenceIdentity {
    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub const fn canonical_fingerprint(&self) -> &CanonicalReferenceFingerprint {
        &self.canonical_fingerprint
    }

    /// Huella canónica de una referencia materializada.
    ///
    /// Falla si la referencia se validó sin conservar sus bases: la huella se
    /// calcula sobre ellas y no puede derivarse de la ficha.
    pub fn from_reference(reference: &Reference) -> Result<Self, ReferenceError> {
        let mut sequences = Vec::with_capacity(reference.sequences().len());
        for sequence in reference.sequences() {
            sequences.push((sequence.id().as_str(), sequence.sequence().as_bytes()?));
        }
        Ok(Self::from_sequences(reference.id().clone(), sequences))
    }

    fn from_authoritative_fasta(reference_id: ReferenceId, fasta: &AuthoritativeFasta<'_>) -> Self {
        let mut sequences = fasta.sequences.iter().collect::<Vec<_>>();
        sequences.sort_by_key(|(id, _)| fasta.ordinal(id));
        Self::from_sequences(
            reference_id,
            sequences
                .into_iter()
                .map(|(id, bases)| (id.as_str(), bases.as_ref())),
        )
    }

    fn from_sequences<'a>(
        reference_id: ReferenceId,
        sequences: impl IntoIterator<Item = (&'a str, &'a [u8])>,
    ) -> Self {
        let mut digest = Sha256::new();
        digest.update(b"genome_population_simulator/v1/prepared-reference\0");
        append_digest_component(&mut digest, reference_id.as_str().as_bytes());
        for (sequence_id, bases) in sequences {
            append_digest_component(&mut digest, sequence_id.as_bytes());
            append_digest_component(&mut digest, bases);
        }
        Self {
            reference_id,
            canonical_fingerprint: CanonicalReferenceFingerprint(format!(
                "{:x}",
                digest.finalize()
            )),
        }
    }
}

/// Recurso preparado en memoria; no publica ni altera la tabla suministrada.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedPopulationSnvs {
    prepared_reference: PreparedReferenceIdentity,
    rows: Vec<PreparedPopulationSnv>,
    population_ids: Vec<PopulationId>,
    discarded_globally_absent: Vec<VariantId>,
    empty_authorized: bool,
    normalization_audit: Vec<AltFrequencyNormalization>,
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
    /// Identidad canónica de la referencia validada al preparar esta matriz.
    #[must_use]
    pub const fn prepared_reference(&self) -> &PreparedReferenceIdentity {
        &self.prepared_reference
    }

    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        self.prepared_reference.reference_id()
    }

    #[must_use]
    pub fn reference_fingerprint(&self) -> &str {
        self.prepared_reference.canonical_fingerprint().as_str()
    }

    #[must_use]
    pub fn rows(&self) -> &[PreparedPopulationSnv] {
        &self.rows
    }
    #[must_use]
    pub fn population_ids(&self) -> &[PopulationId] {
        &self.population_ids
    }
    #[must_use]
    pub fn discarded_globally_absent(&self) -> &[VariantId] {
        &self.discarded_globally_absent
    }
    #[must_use]
    pub const fn empty_authorized(&self) -> bool {
        self.empty_authorized
    }
    #[must_use]
    pub fn normalization_audit(&self) -> &[AltFrequencyNormalization] {
        &self.normalization_audit
    }

    /// Frecuencias alélicas preparadas para una única población×locus.
    ///
    /// La frecuencia REF procede de la decisión de normalización ya retenida:
    /// en particular, una decisión `normalized` conserva REF exactamente a
    /// cero y ningún consumidor vuelve a calcular `1 - sum(ALT)`.
    #[must_use]
    pub fn allele_frequencies_for_rows(&self, rows: &[&PreparedPopulationSnv]) -> Option<Vec<f64>> {
        let first = *rows.first()?;
        if rows.iter().any(|row| {
            row.population_id() != first.population_id()
                || row.contig() != first.contig()
                || row.position_1based() != first.position_1based()
                || row.reference() != first.reference()
        }) {
            return None;
        }
        let position = first.position_1based().to_string();
        let reference = char::from(first.reference()).to_string();
        let group = self.normalization_audit.iter().find(|group| {
            group.population_id() == first.population_id().as_str()
                && group.contig() == first.contig().as_str()
                && group.position() == position
                && group.reference() == reference
        })?;
        let mut expected = rows.to_vec();
        expected.sort_by(|left, right| {
            left.alternate()
                .cmp(&right.alternate())
                .then(left.variant_id().cmp(right.variant_id()))
        });
        let active_variants = expected
            .iter()
            .map(|row| (row.variant_id().as_str(), (row.alternate(), row.alt_af())))
            .collect::<BTreeMap<_, _>>();
        let discarded = self
            .discarded_globally_absent
            .iter()
            .map(VariantId::as_str)
            .collect::<BTreeSet<_>>();
        let mut seen_active = BTreeSet::new();
        for value in group.values() {
            let variant_id = value.variant_id();
            if let Some((alternate, frequency)) = active_variants.get(variant_id) {
                if value.alt().as_bytes() != [*alternate]
                    || value.effective().to_bits() != frequency.to_bits()
                    || !seen_active.insert(variant_id)
                {
                    return None;
                }
            } else if value.effective() != 0.0 || !discarded.contains(variant_id) {
                return None;
            }
        }
        if seen_active.len() != expected.len() {
            return None;
        }
        let mut frequencies = Vec::with_capacity(expected.len() + 1);
        frequencies.push(group.reference_frequency());
        frequencies.extend(expected.iter().map(|row| row.alt_af()));
        Some(frequencies)
    }

    /// Huella SHA-256 de la proyección canónica exacta consumida por F2.3.
    #[must_use]
    pub fn canonical_fingerprint(&self) -> String {
        use sha2::{Digest, Sha256};

        format!("{:x}", Sha256::digest(canonical_tsv_bytes(self)))
    }

    /// Exige la misma confirmación explícita para una matriz vacía generada
    /// que para una matriz suministrada o reutilizada.
    pub fn authorize_empty_matrix(
        mut self,
        allow_empty_matrix: bool,
    ) -> Result<Self, PopulationSnvsPreparationError> {
        if self.rows.is_empty() {
            if !allow_empty_matrix {
                return Err(rejected(
                    PopulationSnvsPreparationError::EmptyActiveMatrixRequiresAuthorization,
                ));
            }
            self.empty_authorized = true;
        }
        Ok(self)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum PopulationSnvsPreparationError {
    /// Se pidieron bases de una referencia validada sin conservarlas.
    BasesNotMaterialized,
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
    let prepared_reference = PreparedReferenceIdentity::from_authoritative_fasta(
        reference.reference_id().clone(),
        &fasta,
    );
    prepare_population_snvs(
        table,
        &prepared_reference,
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
    let prepared_reference = PreparedReferenceIdentity::from_reference(reference)
        .map_err(|_| rejected(PopulationSnvsPreparationError::BasesNotMaterialized))?;
    let mut fasta_sequences = BTreeMap::new();
    for sequence in reference.sequences() {
        fasta_sequences.insert(
            sequence.id().as_str().to_owned(),
            Cow::Borrowed(
                sequence
                    .sequence()
                    // Inalcanzable porque `PreparedReferenceIdentity::from_reference`, que corre justo
                    // antes, ya pide `as_bytes()` de cada secuencia y devuelve `BasesNotMaterialized`
                    // (línea anterior) si alguna no está materializada: aquí todas lo están.
                    .as_bytes()
                    .map_err(|_| rejected(PopulationSnvsPreparationError::BasesNotMaterialized))?,
            ),
        );
    }
    let fasta = AuthoritativeFasta {
        sequences: fasta_sequences,
        order: reference
            .sequences()
            .iter()
            .enumerate()
            .map(|(index, sequence)| (sequence.id().as_str().to_owned(), index))
            .collect(),
    };
    prepare_population_snvs(
        table,
        &prepared_reference,
        &fasta,
        tolerance,
        allow_empty_matrix,
    )
}

fn prepare_population_snvs(
    table: &PopulationSnvsTsv,
    prepared_reference: &PreparedReferenceIdentity,
    fasta: &AuthoritativeFasta<'_>,
    tolerance: f64,
    allow_empty_matrix: bool,
) -> Result<PreparedPopulationSnvs, PopulationSnvsPreparationError> {
    let mut normalized = normalize_alt_frequencies(table, tolerance)
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
                rejected(PopulationSnvsPreparationError::InvalidPopulation {
                    row_index: i,
                    value: row.field(5).to_owned(),
                })
            })
        })
        .collect::<Result<BTreeSet<_>, _>>()?;
    let mut seen = BTreeSet::new();
    let mut alt_populations: BTreeMap<(String, u64, u8, u8), BTreeSet<PopulationId>> =
        BTreeMap::new();
    let mut parsed = Vec::with_capacity(table.rows().len());
    for (row_index, row) in table.rows().iter().enumerate() {
        let contig = SequenceId::parse(row.field(1), "population_snvs").map_err(|_| {
            rejected(PopulationSnvsPreparationError::InvalidContig {
                row_index,
                value: row.field(1).to_owned(),
            })
        })?;
        let position = parse_canonical_position(row.field(2)).ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::InvalidPosition {
                row_index,
                value: row.field(2).to_owned(),
            })
        })?;
        let reference_base = parse_snv_base(row.field(3), row_index, "ref")?;
        let alternate = parse_snv_base(row.field(4), row_index, "alt")?;
        if reference_base == alternate {
            return Err(rejected(PopulationSnvsPreparationError::InvalidAllele {
                row_index,
                field: "alt",
                value: row.field(4).to_owned(),
            }));
        }
        // Inalcanzable porque el primer recorrido de `prepare_population_snvs` (el que
        // reúne `populations`) ya llamó a `PopulationId::parse(row.field(5), ..)` sobre
        // cada fila y devolvió `InvalidPopulation` si fallaba: es el mismo texto.
        let population = PopulationId::parse(row.field(5), "population_snvs").map_err(|_| {
            rejected(PopulationSnvsPreparationError::InvalidPopulation {
                row_index,
                value: row.field(5).to_owned(),
            })
        })?;
        let sequence = fasta.sequences.get(contig.as_str()).ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::InvalidContig {
                row_index,
                value: contig.as_str().to_owned(),
            })
        })?;
        let offset = usize::try_from(position - 1).ok();
        let expected = offset
            .and_then(|index| sequence.get(index).copied())
            .ok_or_else(|| {
                rejected(PopulationSnvsPreparationError::ReferenceMismatch {
                    row_index,
                    contig: contig.as_str().to_owned(),
                    position,
                    expected: b'?',
                    found: reference_base,
                })
            })?;
        if expected != reference_base {
            return Err(rejected(
                PopulationSnvsPreparationError::ReferenceMismatch {
                    row_index,
                    contig: contig.as_str().to_owned(),
                    position,
                    expected,
                    found: reference_base,
                },
            ));
        }
        let key = (
            contig.as_str().to_owned(),
            position,
            reference_base,
            alternate,
            population.as_str().to_owned(),
        );
        if !seen.insert(key) {
            return Err(rejected(
                PopulationSnvsPreparationError::DuplicateAltPopulation {
                    contig: contig.as_str().to_owned(),
                    position,
                    reference: reference_base,
                    alternate,
                    population: population.as_str().to_owned(),
                },
            ));
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
    let variant_ids = parsed
        .iter()
        .map(|(_, contig, position, reference_base, alternate, ..)| {
            derive_matrix_variant_id(
                prepared_reference.reference_id(),
                contig,
                *position,
                *reference_base,
                *alternate,
            )
            .map(|variant_id| {
                (
                    (
                        contig.as_str().to_owned(),
                        *position,
                        *reference_base,
                        *alternate,
                    ),
                    variant_id,
                )
            })
            // Inalcanzable porque `derive_matrix_variant_id` solo falla si un componente supera
            // `u32::MAX` bytes (`append_length_prefixed`) y sus componentes son el ID de la
            // referencia y el de la secuencia (64 bytes como máximo, ya validados), la
            // posición y dos bases.
            .map_err(|error| {
                rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                    detail: error.render_human(),
                })
            })
        })
        .collect::<Result<BTreeMap<_, _>, _>>()?;
    for group in &mut normalized.groups {
        let contig = group.contig().to_owned();
        let position_text = group.position().to_owned();
        let reference = group.reference().to_owned();
        for value in &mut group.values {
            let reference_base = parse_snv_base(&reference, value.row_index(), "ref")?;
            let alternate = parse_snv_base(value.alt(), value.row_index(), "alt")?;
            // Inalcanzable porque `position_text` es el `row.field(2)` de una fila que el
            // bucle anterior ya analizó como `u64` positivo: si no, devolvió `InvalidPosition`
            // antes de llegar aquí.
            let position = position_text.parse::<u64>().map_err(|_| {
                rejected(PopulationSnvsPreparationError::InvalidPosition {
                    row_index: value.row_index(),
                    value: position_text.clone(),
                })
            })?;
            value.variant_id = variant_ids
                .get(&(contig.clone(), position, reference_base, alternate))
                // Inalcanzable porque cada valor de cada grupo procede de una fila de `parsed`: su
                // contig, posición, REF y ALT son los mismos textos con los que se rellenó
                // `variant_ids`, de modo que la clave siempre existe.
                .ok_or_else(|| {
                    rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                        detail: "invariante de identidad de auditoría".to_owned(),
                    })
                })?
                .as_str()
                .to_owned();
        }
    }
    if alt_populations
        .values()
        .any(|covered| covered != &populations)
    {
        let Some(((contig, position, reference_base, alternate), _)) = alt_populations
            // Inalcanzable porque este `find` usa el mismo predicado que el `any` que acaba
            // de ser cierto: siempre encuentra el grupo.
            .iter()
            .find(|(_, covered)| *covered != &populations)
        else {
            return Err(rejected(
                PopulationSnvsPreparationError::AuthoritativeFasta {
                    detail: "invariante de cobertura de población".to_owned(),
                },
            ));
        };
        return Err(rejected(
            PopulationSnvsPreparationError::IncompletePopulationCoverage {
                contig: contig.clone(),
                position: *position,
                reference: *reference_base,
                alternate: *alternate,
            },
        ));
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
        let variant_id = derive_matrix_variant_id(
            prepared_reference.reference_id(),
            &contig,
            position,
            reference_base,
            alternate,
        )
        // Inalcanzable por la misma razón que el `derive_matrix_variant_id` de
        // `variant_ids` (arriba): componentes de tamaño acotado, y es la misma llamada
        // con los mismos argumentos, que ya salió bien.
        .map_err(|error| {
            rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: error.render_human(),
            })
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
        return Err(rejected(
            PopulationSnvsPreparationError::EmptyActiveMatrixRequiresAuthorization,
        ));
    }
    let empty_authorized = rows.is_empty();
    Ok(PreparedPopulationSnvs {
        prepared_reference: prepared_reference.clone(),
        rows,
        population_ids: populations.into_iter().collect(),
        discarded_globally_absent: discarded,
        empty_authorized,
        normalization_audit: normalized.groups,
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
    reference_fasta_sha256: &str,
    matrix_sha256: &str,
    provenance: PopulationSnvsProvenance,
) -> Value {
    let populations = prepared
        .population_ids()
        .iter()
        .map(PopulationId::as_str)
        .map(|population_id| json!({"population_id": population_id}))
        .collect::<Vec<_>>();
    let mut discarded_variant_pairs = prepared
        .population_ids()
        .iter()
        .flat_map(|population_id| {
            prepared
                .discarded_globally_absent()
                .iter()
                .map(move |variant_id| {
                    json!({
                        "population_id": population_id.as_str(),
                        "variant_id": variant_id.as_str(),
                    })
                })
        })
        .collect::<Vec<_>>();
    discarded_variant_pairs.sort_by(|left, right| {
        (left["population_id"].as_str(), left["variant_id"].as_str()).cmp(&(
            right["population_id"].as_str(),
            right["variant_id"].as_str(),
        ))
    });
    json!({
        "document_kind": "population_snvs_metadata",
        "schema_version": "1",
        "coordinate_system": "base1_closed",
        "reference": {
            "reference_id": prepared.reference_id().as_str(),
            "canonical_fingerprint": prepared.reference_fingerprint(),
        },
        "reference_fasta": {"sha256": reference_fasta_sha256},
        "matrix": {"path": "population_snvs.tsv", "sha256": matrix_sha256},
        "populations": populations,
        "provenance": provenance.as_str(),
        "discarded_variant_pairs": discarded_variant_pairs,
        "normalization_audit": normalization_audit_document(prepared),
    })
}

/// Auditoría completa de normalización que acompaña al par canónico hasta el
/// resultado. Sus valores `effective` son comprobables contra el TSV; los
/// `original` conservan la decisión previa a cualquier normalización.
#[must_use]
pub fn normalization_audit_document(prepared: &PreparedPopulationSnvs) -> Value {
    let mut groups = prepared.normalization_audit().iter().collect::<Vec<_>>();
    groups.sort_by(|left, right| {
        left.population_id()
            .cmp(right.population_id())
            .then(left.contig().cmp(right.contig()))
            .then_with(|| {
                left.position()
                    .parse::<u64>()
                    .unwrap_or(u64::MAX)
                    .cmp(&right.position().parse::<u64>().unwrap_or(u64::MAX))
            })
            .then(left.reference().cmp(right.reference()))
    });
    Value::Array(
        groups
            .into_iter()
            .map(|group| {
                json!({
                    "population_id": group.population_id(),
                    "contig": group.contig(),
                    "position": group.position(),
                    "ref": group.reference(),
                    "tolerance": group.tolerance(),
                    "canonical_order": group.canonical_order(),
                    "sum_algorithm": group.sum_algorithm(),
                    "original_sum": group.original_sum(),
                    "original_delta": group.original_delta(),
                    "effective_sum": group.effective_sum(),
                    "effective_delta": group.effective_delta(),
                    "reference_frequency": group.reference_frequency(),
                    "decision": group.decision().as_str(),
                    "values": group.values().iter().map(|value| json!({
                        "alt": value.alt(),
                        "variant_id": value.variant_id(),
                        "source_id": value.source_id(),
                        "original": value.original(),
                        "effective": value.effective(),
                    })).collect::<Vec<_>>(),
                })
            })
            .collect(),
    )
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
    (
        PreparedPopulationSnvs,
        PopulationSnvsProvenance,
        String,
        Value,
    ),
    PopulationSnvsPreparationError,
> {
    validate_reusable_population_snvs_with_matrix_sha256(
        metadata_path,
        reference,
        reference_sha256,
        tolerance,
        allow_empty_matrix,
    )
    .map(|(prepared, provenance, metadata_sha256, audit, _)| {
        (prepared, provenance, metadata_sha256, audit)
    })
}

/// Igual que [`validate_reusable_population_snvs`], pero conserva además el
/// SHA-256 de los bytes TSV que se validaron. Los consumidores que publican
/// procedencia de verdad deben usar esta identidad del recurso, nunca el hash
/// del documento metadata que lo describe.
pub fn validate_reusable_population_snvs_with_matrix_sha256(
    metadata_path: &Path,
    reference: &Reference,
    reference_sha256: &str,
    tolerance: f64,
    allow_empty_matrix: bool,
) -> Result<
    (
        PreparedPopulationSnvs,
        PopulationSnvsProvenance,
        String,
        Value,
        String,
    ),
    PopulationSnvsPreparationError,
> {
    validate_reusable_population_snvs_snapshots(
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
    (
        PreparedPopulationSnvs,
        PopulationSnvsProvenance,
        String,
        Value,
    ),
    PopulationSnvsPreparationError,
> {
    validate_reusable_population_snvs_snapshots(
        metadata_path,
        reference,
        reference_sha256,
        tolerance,
        allow_empty_matrix,
        after_metadata_snapshot,
        after_matrix_snapshot,
    )
    .map(|(prepared, provenance, metadata_sha256, audit, _)| {
        (prepared, provenance, metadata_sha256, audit)
    })
}

fn validate_reusable_population_snvs_snapshots(
    metadata_path: &Path,
    reference: &Reference,
    reference_sha256: &str,
    tolerance: f64,
    allow_empty_matrix: bool,
    after_metadata_snapshot: impl FnOnce(),
    after_matrix_snapshot: impl FnOnce(),
) -> Result<
    (
        PreparedPopulationSnvs,
        PopulationSnvsProvenance,
        String,
        Value,
        String,
    ),
    PopulationSnvsPreparationError,
> {
    let metadata_bytes = fs::read(metadata_path).map_err(|error| {
        rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: format!("no se puede leer metadata_path: {error}"),
        })
    })?;
    after_metadata_snapshot();
    let matrix_path = metadata_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("population_snvs.tsv");
    let matrix_bytes = fs::read(&matrix_path).map_err(|error| {
        rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: format!(
                "no se puede leer la matriz {}, que debe estar junto a los metadatos: {error}",
                matrix_path.display()
            ),
        })
    })?;
    after_matrix_snapshot();
    let (prepared, provenance, metadata_sha256, audit) = validate_reusable_population_snvs_bytes(
        &metadata_bytes,
        &matrix_bytes,
        reference,
        reference_sha256,
        tolerance,
        allow_empty_matrix,
    )?;
    Ok((
        prepared,
        provenance,
        metadata_sha256,
        audit,
        sha256_bytes(&matrix_bytes),
    ))
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
    (
        PreparedPopulationSnvs,
        PopulationSnvsProvenance,
        String,
        Value,
    ),
    PopulationSnvsPreparationError,
> {
    let metadata = parse_metadata_json(metadata_bytes).map_err(|error| {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: error.to_string(),
        }
    })?;
    let document = metadata.document();
    // Inalcanzable porque `parse_metadata_json`, que acaba de aceptar estos mismos
    // bytes, exige un objeto en la raíz (`MetadataMustBeObject`).
    let object = document.as_object().ok_or_else(|| {
        rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "metadata no es un objeto".to_owned(),
        })
    })?;
    let reference_object = object
        .get("reference")
        .and_then(Value::as_object)
        // Inalcanzable porque `parse_metadata_json` exige `reference` como objeto
        // (`required_object`).
        .ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.reference inválida".to_owned(),
            })
        })?;
    let expected_reference = PreparedReferenceIdentity::from_reference(reference)
        .map_err(|_| rejected(PopulationSnvsPreparationError::BasesNotMaterialized))?;
    let declared_reference = reference_object.get("reference_id").and_then(Value::as_str);
    let declared_fingerprint = reference_object
        .get("canonical_fingerprint")
        .and_then(Value::as_str);
    let fasta_object = object
        .get("reference_fasta")
        .and_then(Value::as_object)
        // Inalcanzable porque `parse_metadata_json` exige `reference_fasta` como objeto
        // (`required_object`).
        .ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.reference_fasta inválida".to_owned(),
            })
        })?;
    let declared_fasta_sha = fasta_object.get("sha256").and_then(Value::as_str);
    if declared_reference != Some(expected_reference.reference_id().as_str())
        || declared_fingerprint != Some(expected_reference.canonical_fingerprint().as_str())
        || declared_fasta_sha != Some(reference_sha256)
    {
        return Err(rejected(
            PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "la matriz reutilizada no está ligada a la referencia autoritativa"
                    .to_owned(),
            },
        ));
    }
    let matrix = object
        .get("matrix")
        .and_then(Value::as_object)
        // Inalcanzable porque `parse_metadata_json` exige `matrix` como objeto
        // (`required_object`).
        .ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.matrix inválida".to_owned(),
            })
        })?;
    let declared_digest = matrix
        .get("sha256")
        .and_then(Value::as_str)
        // Inalcanzable porque `parse_metadata_json` exige `matrix.sha256` como cadena
        // (`required_string`).
        .ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.matrix.sha256 inválido".to_owned(),
            })
        })?;
    if sha256_bytes(matrix_bytes) != declared_digest {
        return Err(rejected(
            PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "SHA-256 del TSV reutilizado no coincide con metadata".to_owned(),
            },
        ));
    }
    // Inalcanzable porque `parse_metadata_json` exige que `normalization_audit` exista
    // y sea una lista (su última comprobación).
    let audit = object.get("normalization_audit").ok_or_else(|| {
        rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: "metadata.normalization_audit ausente".to_owned(),
        })
    })?;
    let table = parse_canonical_tsv(matrix_bytes).map_err(|error| {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: error.to_string(),
        }
    })?;
    let mut prepared =
        prepare_population_snvs_from_reference(&table, reference, tolerance, allow_empty_matrix)?;
    restore_published_allele_frequencies(&mut prepared, audit)?;
    if canonical_tsv_bytes(&prepared) != matrix_bytes {
        return Err(rejected(
            PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "el TSV reutilizado no está en orden y representación canónicos".to_owned(),
            },
        ));
    }
    let populations = object
        .get("populations")
        .and_then(Value::as_array)
        // Inalcanzable porque `parse_metadata_json` exige `populations` como lista.
        .ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.populations inválida".to_owned(),
            })
        })?;
    let observed = populations
        .iter()
        .map(|value| {
            value
                .get("population_id")
                .and_then(Value::as_str)
                // Inalcanzable porque `parse_metadata_json` exige que cada población tenga
                // `population_id` como cadena (`required_string`). Su sintaxis sí la juzga la
                // guarda siguiente, que sí es alcanzable.
                .ok_or_else(|| {
                    rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                        detail: "metadata.populations inválida".to_owned(),
                    })
                })
                .and_then(|value| {
                    PopulationId::parse(value, "population_snvs_metadata").map_err(|_| {
                        rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                            detail: "metadata.populations inválida".to_owned(),
                        })
                    })
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    if observed.windows(2).any(|window| window[0] >= window[1]) {
        return Err(rejected(
            PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.populations no está ordenada canónicamente".to_owned(),
            },
        ));
    }
    if prepared.population_ids().is_empty() {
        prepared.population_ids = observed;
    } else if observed != prepared.population_ids() {
        return Err(rejected(
            PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.populations no coincide con el TSV canónico".to_owned(),
            },
        ));
    }
    let provenance = object
        .get("provenance")
        .and_then(Value::as_str)
        .and_then(PopulationSnvsProvenance::parse)
        // Inalcanzable porque `parse_metadata_json` exige `provenance` como cadena y
        // solo admite `supplied` o `generated`, que son exactamente lo que acepta
        // `PopulationSnvsProvenance::parse` (`required_provenance`).
        .ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.provenance inválida".to_owned(),
            })
        })?;
    let discarded_variant_pairs = object
        .get("discarded_variant_pairs")
        .and_then(Value::as_array)
        // Inalcanzable porque `parse_metadata_json` exige `discarded_variant_pairs` como lista.
        .ok_or_else(|| {
            rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "metadata.discarded_variant_pairs inválido".to_owned(),
            })
        })?
        .iter()
        .map(|value| {
            // Inalcanzable porque `parse_metadata_json` exige que cada par sea un objeto.
            let pair = value.as_object().ok_or_else(|| {
                rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                    detail: "metadata.discarded_variant_pairs inválido".to_owned(),
                })
            })?;
            Ok((
                pair.get("population_id")
                    .and_then(Value::as_str)
                    .map(str::to_owned)
                    // Inalcanzable porque `parse_metadata_json` exige `population_id` como cadena en
                    // cada par descartado.
                    .ok_or_else(|| {
                        rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                            detail: "metadata.discarded_variant_pairs inválido".to_owned(),
                        })
                    })?,
                pair.get("variant_id")
                    .and_then(Value::as_str)
                    .map(str::to_owned)
                    // Inalcanzable porque `parse_metadata_json` exige `variant_id` como cadena en
                    // cada par descartado.
                    .ok_or_else(|| {
                        rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
                            detail: "metadata.discarded_variant_pairs inválido".to_owned(),
                        })
                    })?,
            ))
        })
        .collect::<Result<BTreeSet<_>, _>>()?;
    validate_normalization_audit_document(
        audit,
        &prepared,
        reference,
        tolerance,
        &discarded_variant_pairs,
    )?;
    Ok((
        prepared,
        provenance,
        sha256_bytes(metadata_bytes),
        audit.clone(),
    ))
}

/// Restaura la semántica de frecuencias que acompaña al TSV canónico.
///
/// El TSV conserva los ALT efectivos, pero una normalización que declaró REF
/// exactamente cero no puede reconstruirse con seguridad como `1 - sum(ALT)`
/// tras una segunda lectura de `f64`. La auditoría autenticada es por tanto la
/// representación preparada autoritativa para REF y para cada ALT.
fn restore_published_allele_frequencies(
    prepared: &mut PreparedPopulationSnvs,
    audit: &Value,
) -> Result<(), PopulationSnvsPreparationError> {
    type AlleleKey = (String, String, String, String, String, String, String);
    type GroupKey = (String, String, String, String);
    let groups = audit
        // Inalcanzable por el camino público: esta función es privada y solo la llama
        // `validate_reusable_population_snvs_bytes`, tras `parse_metadata_json`, que exige
        // que `normalization_audit` sea una lista. La misma comprobación de
        // `validate_normalization_audit_document` (pública) sí se prueba.
        .as_array()
        .ok_or_else(|| audit_error("debe ser una lista"))?;
    let mut values = BTreeMap::<AlleleKey, f64>::new();
    let mut references = BTreeMap::<GroupKey, f64>::new();
    for group in groups {
        let group = group
            .as_object()
            .ok_or_else(|| audit_error("grupo inválido"))?;
        let population = group
            .get("population_id")
            .and_then(Value::as_str)
            .ok_or_else(|| audit_error("population_id"))?;
        let contig = group
            .get("contig")
            .and_then(Value::as_str)
            .ok_or_else(|| audit_error("contig"))?;
        let position = group
            .get("position")
            .and_then(Value::as_str)
            .ok_or_else(|| audit_error("position"))?;
        let reference = group
            .get("ref")
            .and_then(Value::as_str)
            .ok_or_else(|| audit_error("ref"))?;
        let reference_frequency = group
            .get("reference_frequency")
            .and_then(Value::as_f64)
            // `is_finite` no discrimina (anotado, no quitado): JSON no puede
            // representar `NaN` ni el infinito, y `1e999` se lee como número pero
            // `as_f64()` devuelve `None` (no un infinito), así que un `f64` leído de
            // aquí siempre es finito; solo `>= 0.0` actúa. Lo comprueba la prueba
            // `an_out_of_range_json_number_is_not_a_finite_float_for_the_audit_judge`.
            .filter(|value| value.is_finite() && *value >= 0.0)
            .ok_or_else(|| audit_error("reference_frequency"))?;
        if references
            .insert(
                (
                    population.to_owned(),
                    contig.to_owned(),
                    position.to_owned(),
                    reference.to_owned(),
                ),
                reference_frequency,
            )
            .is_some()
        {
            return Err(audit_error("grupo duplicado"));
        }
        for value in group
            .get("values")
            .and_then(Value::as_array)
            .ok_or_else(|| audit_error("values"))?
        {
            let value = value
                .as_object()
                .ok_or_else(|| audit_error("valor inválido"))?;
            let alt = value
                .get("alt")
                .and_then(Value::as_str)
                .ok_or_else(|| audit_error("alt"))?;
            let variant_id = value
                .get("variant_id")
                .and_then(Value::as_str)
                .ok_or_else(|| audit_error("variant_id"))?;
            let source_id = value
                .get("source_id")
                .and_then(Value::as_str)
                .ok_or_else(|| audit_error("source_id"))?;
            let effective = value
                .get("effective")
                .and_then(Value::as_f64)
                // Igual que `reference_frequency`: `is_finite` no discrimina.
                .filter(|value| value.is_finite() && *value >= 0.0)
                .ok_or_else(|| audit_error("effective"))?;
            if values
                .insert(
                    (
                        population.to_owned(),
                        contig.to_owned(),
                        position.to_owned(),
                        reference.to_owned(),
                        alt.to_owned(),
                        variant_id.to_owned(),
                        source_id.to_owned(),
                    ),
                    effective,
                )
                .is_some()
            {
                return Err(audit_error("valor duplicado"));
            }
        }
    }
    for row in &mut prepared.rows {
        let key = (
            row.population_id().as_str().to_owned(),
            row.contig().as_str().to_owned(),
            row.position_1based().to_string(),
            char::from(row.reference()).to_string(),
            char::from(row.alternate()).to_string(),
            row.variant_id().as_str().to_owned(),
            row.source_id().to_owned(),
        );
        row.alt_af = *values
            .get(&key)
            .ok_or_else(|| audit_error("efectivo activo ausente"))?;
    }
    for group in &mut prepared.normalization_audit {
        let key = (
            group.population_id.clone(),
            group.contig.clone(),
            group.position.clone(),
            group.reference.clone(),
        );
        group.reference_frequency = *references
            .get(&key)
            .ok_or_else(|| audit_error("referencia ausente"))?;
        for value in &mut group.values {
            let key = (
                group.population_id.clone(),
                group.contig.clone(),
                group.position.clone(),
                group.reference.clone(),
                value.alt.clone(),
                value.variant_id.clone(),
                value.source_id.clone(),
            );
            value.effective = *values
                .get(&key)
                .ok_or_else(|| audit_error("efectivo ausente"))?;
        }
    }
    Ok(())
}

/// Comprueba que la auditoría publicada conserva la decisión original y que
/// sus valores efectivos siguen correspondiendo exactamente al TSV canónico.
pub fn validate_normalization_audit_document(
    audit: &Value,
    prepared: &PreparedPopulationSnvs,
    reference: &Reference,
    tolerance: f64,
    declared_discarded_variant_pairs: &BTreeSet<(String, String)>,
) -> Result<(), PopulationSnvsPreparationError> {
    let groups = audit
        .as_array()
        .ok_or_else(|| audit_error("debe ser una lista"))?;
    let mut active = prepared
        .rows()
        .iter()
        .map(|row| {
            (
                (
                    row.population_id().as_str().to_owned(),
                    row.contig().as_str().to_owned(),
                    row.position_1based().to_string(),
                    char::from(row.reference()).to_string(),
                    char::from(row.alternate()).to_string(),
                    row.variant_id().as_str().to_owned(),
                    row.source_id().to_owned(),
                ),
                row.alt_af(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let mut seen = BTreeSet::new();
    let mut previous_group = None;
    let mut absent_variant_pairs = BTreeSet::new();
    for group in groups {
        let object = group
            .as_object()
            .ok_or_else(|| audit_error("grupo inválido"))?;
        const FIELDS: [&str; 14] = [
            "population_id",
            "contig",
            "position",
            "ref",
            "tolerance",
            "canonical_order",
            "sum_algorithm",
            "original_sum",
            "original_delta",
            "effective_sum",
            "effective_delta",
            "reference_frequency",
            "decision",
            "values",
        ];
        if object.len() != FIELDS.len()
            || object.keys().any(|field| !FIELDS.contains(&field.as_str()))
        {
            return Err(audit_error("campos de grupo inválidos"));
        }
        let text = |field| {
            object
                .get(field)
                .and_then(Value::as_str)
                .ok_or_else(|| audit_error(field))
        };
        let number = |field| {
            object
                .get(field)
                .and_then(Value::as_f64)
                .filter(|value| value.is_finite())
                .ok_or_else(|| audit_error(field))
        };
        let population = text("population_id")?;
        let contig = text("contig")?;
        let position_text = text("position")?;
        let reference_text = text("ref")?;
        PopulationId::parse(population, "normalization_audit")
            .map_err(|_| audit_error("population_id"))?;
        if !prepared
            .population_ids()
            .iter()
            .any(|known| known.as_str() == population)
        {
            return Err(audit_error("population_id no pertenece a la matriz"));
        }
        let position =
            parse_canonical_position(position_text).ok_or_else(|| audit_error("position"))?;
        let reference_base =
            parse_snv_base(reference_text, 0, "ref").map_err(|_| audit_error("ref"))?;
        let sequence_id =
            SequenceId::parse(contig, "normalization_audit").map_err(|_| audit_error("contig"))?;
        let sequence = reference
            .sequences()
            .iter()
            .find(|sequence| sequence.id() == &sequence_id)
            .ok_or_else(|| audit_error("contig ausente de referencia"))?;
        if sequence
            .sequence()
            .as_bytes()
            .map_err(|_| rejected(PopulationSnvsPreparationError::BasesNotMaterialized))?
            // Inalcanzable porque `usize::try_from(u64)` no falla en 64 bits (el único
            // destino que se compila y se publica) y `position - 1` no desborda: se filtró
            // `position > 0`.
            .get(usize::try_from(position - 1).map_err(|_| audit_error("position"))?)
            .copied()
            != Some(reference_base)
        {
            return Err(audit_error("ref no coincide con referencia"));
        }
        let group_key = (population, contig, position, reference_text);
        if number("tolerance")?.to_bits() != tolerance.to_bits()
            || text("canonical_order")? != ALT_FREQUENCY_CANONICAL_ORDER
            || text("sum_algorithm")? != ALT_FREQUENCY_SUM_ALGORITHM
            // Las dos últimas mitades se solapan (anotado, no quitado): el orden
            // estricto con el grupo anterior ya impide repetir una clave, porque la
            // secuencia aceptada es estrictamente creciente. `!seen.insert(..)` solo
            // actúa si falla antes `last >= key` y, por tanto, no cambia ningún
            // resultado; por la misma razón `>=` y `>` dan el mismo veredicto (con
            // `>` una clave igual la caza `seen`). El orden descendente sí lo
            // caza únicamente la primera de las dos.
            || previous_group
                .as_ref()
                .is_some_and(|last| last >= &group_key)
            || !seen.insert(group_key)
        {
            return Err(audit_error("identidad o algoritmo de grupo inválido"));
        }
        previous_group = Some(group_key);
        let values = object
            .get("values")
            .and_then(Value::as_array)
            .ok_or_else(|| audit_error("values"))?;
        if values.is_empty() {
            return Err(audit_error("values vacíos"));
        }
        let mut originals = Vec::with_capacity(values.len());
        let mut effectives = Vec::with_capacity(values.len());
        let mut previous_value = None;
        for value in values {
            let value = value
                .as_object()
                .ok_or_else(|| audit_error("valor inválido"))?;
            const VALUE_FIELDS: [&str; 5] =
                ["alt", "variant_id", "source_id", "original", "effective"];
            if value.len() != VALUE_FIELDS.len()
                || value
                    .keys()
                    .any(|field| !VALUE_FIELDS.contains(&field.as_str()))
            {
                return Err(audit_error("campos de valor inválidos"));
            }
            let alt = value
                .get("alt")
                .and_then(Value::as_str)
                .ok_or_else(|| audit_error("alt"))?;
            let source = value
                .get("source_id")
                .and_then(Value::as_str)
                .ok_or_else(|| audit_error("source_id"))?;
            let variant_id = value
                .get("variant_id")
                .and_then(Value::as_str)
                // `!value.is_empty()` no discrimina (anotado, no quitado): un texto vacío
                // lo rechazaría igual `VariantId::parse` en la línea siguiente (exige 1 a
                // 64 bytes) con el mismo motivo, `variant_id`.
                .filter(|value| !value.is_empty())
                .ok_or_else(|| audit_error("variant_id"))?;
            VariantId::parse(variant_id, "normalization_audit")
                .map_err(|_| audit_error("variant_id"))?;
            let alternate = parse_snv_base(alt, 0, "alt").map_err(|_| audit_error("alt"))?;
            if alternate == reference_base {
                return Err(audit_error("alt igual a ref"));
            }
            let expected_variant_id = derive_matrix_variant_id(
                reference.id(),
                &sequence_id,
                position,
                reference_base,
                alternate,
            )
            // Inalcanzable porque `derive_matrix_variant_id` solo falla si un componente supera
            // `u32::MAX` bytes; aquí son IDs de 64 bytes como máximo, una posición y dos bases.
            .map_err(|_| audit_error("variant_id"))?;
            if variant_id != expected_variant_id.as_str() {
                return Err(audit_error("variant_id no canónico"));
            }
            let value_key = (alt, variant_id, source);
            if previous_value
                .as_ref()
                .is_some_and(|last| last >= &value_key)
            {
                return Err(audit_error("orden canónico de values inválido"));
            }
            previous_value = Some(value_key);
            let original = value
                .get("original")
                .and_then(Value::as_f64)
                // `is_finite` no discrimina: ver `reference_frequency` en la restauración;
                // el intervalo cerrado ya excluye `NaN` y el infinito.
                .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
                .ok_or_else(|| audit_error("original"))?;
            let effective = value
                .get("effective")
                .and_then(Value::as_f64)
                // Igual que `original`: `is_finite` no discrimina.
                .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
                .ok_or_else(|| audit_error("effective"))?;
            let key = (
                population.to_owned(),
                contig.to_owned(),
                position.to_string(),
                reference_text.to_owned(),
                alt.to_owned(),
                variant_id.to_owned(),
                source.to_owned(),
            );
            if let Some(expected) = active.remove(&key) {
                if expected.to_bits() != effective.to_bits() {
                    return Err(audit_error("effective no coincide con TSV"));
                }
            } else if effective != 0.0 {
                return Err(audit_error("valor efectivo activo ausente del TSV"));
            } else {
                absent_variant_pairs.insert((population.to_owned(), variant_id.to_owned()));
            }
            originals.push(original);
            effectives.push(effective);
        }
        let original_sum = neumaier(originals.iter().copied());
        let effective_sum = neumaier(effectives.iter().copied());
        let original_delta = (original_sum - 1.0).abs();
        let effective_delta = (effective_sum - 1.0).abs();
        if number("original_sum")?.to_bits() != original_sum.to_bits()
            || number("original_delta")?.to_bits() != original_delta.to_bits()
            || number("effective_sum")?.to_bits() != effective_sum.to_bits()
            || number("effective_delta")?.to_bits() != effective_delta.to_bits()
        {
            return Err(audit_error("sumas no coinciden"));
        }
        match text("decision")? {
            "reference_remainder" if original_sum <= 1.0 => {
                if originals
                    .iter()
                    .map(|value| value.to_bits())
                    .ne(effectives.iter().map(|value| value.to_bits()))
                    || number("reference_frequency")?.to_bits() != (1.0 - effective_sum).to_bits()
                {
                    return Err(audit_error("resto de referencia inválido"));
                }
            }
            "normalized" if original_sum > 1.0 && original_delta <= tolerance => {
                if effectives
                    .iter()
                    .zip(&originals)
                    .any(|(effective, original)| {
                        effective.to_bits() != (original / original_sum).to_bits()
                    })
                    || effective_delta > tolerance
                    || number("reference_frequency")?.to_bits() != 0.0_f64.to_bits()
                {
                    return Err(audit_error("normalización inválida"));
                }
            }
            _ => return Err(audit_error("decisión inválida")),
        }
    }
    if active.is_empty() {
        let discarded_variant_ids = declared_discarded_variant_pairs
            .iter()
            .map(|(_, variant_id)| variant_id)
            .collect::<BTreeSet<_>>();
        let expected_discarded_pairs = prepared
            .population_ids()
            .iter()
            .flat_map(|population_id| {
                discarded_variant_ids.iter().map(move |variant_id| {
                    (population_id.as_str().to_owned(), (*variant_id).to_owned())
                })
            })
            .collect::<BTreeSet<_>>();
        if absent_variant_pairs == *declared_discarded_variant_pairs
            && expected_discarded_pairs == *declared_discarded_variant_pairs
        {
            Ok(())
        } else {
            Err(audit_error(
                "ceros ausentes no coinciden con pares población×variante descartados",
            ))
        }
    } else {
        Err(audit_error("faltan valores efectivos del TSV"))
    }
}

#[track_caller]
fn audit_error(detail: impl Into<String>) -> PopulationSnvsPreparationError {
    crate::guard_trace::record();
    PopulationSnvsPreparationError::AuthoritativeFasta {
        detail: format!("auditoría de normalización inválida: {}", detail.into()),
    }
}

fn sha256_bytes(bytes: &[u8]) -> String {
    let mut digest = Sha256::new();
    digest.update(bytes);
    format!("{:x}", digest.finalize())
}

/// Una posición de 1 base tal como la escribe el contrato (`^[1-9][0-9]*$`):
/// sin signo, sin ceros a la izquierda y mayor que cero. `str::parse::<u64>`
/// admite `"+8"` y `"08"`, que el esquema prohíbe, y por eso no basta.
pub(crate) fn parse_canonical_position(text: &str) -> Option<u64> {
    // Con el primer carácter en 1..=9 solo queda un signo o un dígito de más,
    // y `parse` rechaza lo que no sea un número: no hace falta mirar el resto.
    matches!(text.as_bytes().first(), Some(b'1'..=b'9'))
        .then(|| text.parse::<u64>().ok())
        .flatten()
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
        Err(rejected(PopulationSnvsPreparationError::InvalidAllele {
            row_index,
            field,
            value: value.to_owned(),
        }))
    }
}

/// Bases del FASTA autoritativo indexadas por secuencia.
///
/// Se prestan cuando ya están en memoria dentro de una referencia y se poseen
/// cuando proceden de leer un fichero. Copiarlas siempre duplicaría el genoma:
/// con un cromosoma grande, ese duplicado marca el pico de la ejecución.
struct AuthoritativeFasta<'a> {
    sequences: BTreeMap<String, Cow<'a, [u8]>>,
    order: BTreeMap<String, usize>,
}
impl AuthoritativeFasta<'_> {
    fn ordinal(&self, id: &str) -> usize {
        self.order.get(id).copied().unwrap_or(usize::MAX)
    }
}

fn append_digest_component(digest: &mut Sha256, bytes: &[u8]) {
    digest.update(u64::try_from(bytes.len()).unwrap_or(u64::MAX).to_be_bytes());
    digest.update(bytes);
}

fn load_authoritative_fasta(
    reference: &ResolvedManifestReference,
) -> Result<AuthoritativeFasta<'static>, PopulationSnvsPreparationError> {
    let bytes = fs::read(reference.fasta_path()).map_err(|error| {
        rejected(PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: error.to_string(),
        })
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
                        .insert(previous.clone(), Cow::Owned(std::mem::take(&mut bases)))
                        .is_some()
                {
                    return Err(rejected(
                        PopulationSnvsPreparationError::AuthoritativeFasta {
                            detail: "cabecera FASTA duplicada o vacía".to_owned(),
                        },
                    ));
                }
                order.insert(previous, order.len());
            }
        } else {
            bases.extend(line.iter().map(|base| base.to_ascii_uppercase()));
        }
    }
    let Some(last) = name else {
        return Err(rejected(
            PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "FASTA vacío".to_owned(),
            },
        ));
    };
    if last.is_empty() || records.insert(last.clone(), Cow::Owned(bases)).is_some() {
        return Err(rejected(
            PopulationSnvsPreparationError::AuthoritativeFasta {
                detail: "cabecera FASTA duplicada o vacía".to_owned(),
            },
        ));
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
    // Equivalente (anotado, no quitado): `NaN` e infinito tampoco están en
    // `0.0..=1e-6`, de modo que la mitad `!tolerance.is_finite()` rechaza solo lo que
    // la otra ya rechazaría, con el mismo error. Lo comprueba la prueba
    // `alt_frequencies_are_rejected_for_the_exact_value_that_breaks`.
    if !tolerance.is_finite() || !(0.0..=1e-6).contains(&tolerance) {
        return Err(rejected(AltFrequencyNormalizationError::InvalidTolerance {
            value: tolerance,
        }));
    }
    let mut groups: BTreeMap<AltGroupKey<'_>, Vec<ParsedAltRow<'_>>> = BTreeMap::new();
    for (row_index, row) in table.rows.iter().enumerate() {
        let literal = row.field(6);
        let value = literal.parse::<f64>().map_err(|_| {
            rejected(AltFrequencyNormalizationError::InvalidFrequency {
                row_index,
                literal: literal.to_owned(),
            })
        })?;
        // La mitad `!value.is_finite()` no discrimina (anotado, no quitado): `NaN` y el
        // infinito tampoco pertenecen a `0.0..=1.0`, que rechaza lo mismo con el mismo
        // error (casos `NaN`, `inf` y `-inf` de la prueba de frecuencias).
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(rejected(AltFrequencyNormalizationError::InvalidFrequency {
                row_index,
                literal: literal.to_owned(),
            }));
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
            return Err(rejected(AltFrequencyNormalizationError::ExcessAltSum {
                population_id: population_id.to_owned(),
                contig: contig.to_owned(),
                position: position.to_owned(),
                reference: reference.to_owned(),
                sum: original_sum,
                delta: original_delta,
                tolerance,
            }));
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
        // Inalcanzable porque `effective` es `original` (ya validado en [0, 1]) o
        // `original / original_sum` con `original_sum > 1` y finito, que cae en
        // (0, 1]: el cociente de un término entre la suma que lo incluye no pasa
        // de 1. Se probaron 1,19 millones de grupos normalizados al azar (hasta
        // 40 términos y sumas de 1 a 1 + 5 ulp) sin ningún efectivo fuera de rango.
        if values
            .iter()
            .any(|value| !value.effective.is_finite() || !(0.0..=1.0).contains(&value.effective))
        {
            return Err(rejected(
                AltFrequencyNormalizationError::InvalidNormalizedSum {
                    population_id: population_id.to_owned(),
                    contig: contig.to_owned(),
                    position: position.to_owned(),
                    reference: reference.to_owned(),
                },
            ));
        }
        let effective_sum = neumaier(values.iter().map(AuditedAltFrequency::effective));
        let effective_delta = (effective_sum - 1.0).abs();
        // Inalcanzable porque tras normalizar la suma de los efectivos se desvía de 1
        // como mucho un ulp de 1 (2,2e-16) —el error de redondeo de las divisiones y
        // de la suma—, y `original_delta <= tolerance` con `original_sum > 1` obliga a
        // `tolerance >= 2,2e-16`. La misma búsqueda al azar (1,19 millones de
        // grupos, tolerancia mínima admisible en cada uno) no encontró ninguno con
        // `effective_delta > tolerance`; el máximo observado fue exactamente 1 ulp.
        // La mitad `!effective_sum.is_finite()` tampoco actúa: la suma de valores en
        // [0, 1] es finita.
        if decision == AltFrequencyDecision::Normalized
            && (!effective_sum.is_finite() || effective_delta > tolerance)
        {
            return Err(rejected(
                AltFrequencyNormalizationError::InvalidNormalizedSum {
                    population_id: population_id.to_owned(),
                    contig: contig.to_owned(),
                    position: position.to_owned(),
                    reference: reference.to_owned(),
                },
            ));
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
    let text =
        std::str::from_utf8(bytes).map_err(|_| rejected(PopulationSnvsFormatError::InvalidUtf8))?;
    if text.starts_with("##fileformat=VCF") || text.starts_with("#CHROM\t") {
        return Err(rejected(PopulationSnvsFormatError::VcfInput));
    }
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| rejected(PopulationSnvsFormatError::MissingHeader))?;
    let expected = TSV_HEADER.join("\t");
    if header != expected {
        return Err(rejected(PopulationSnvsFormatError::IncorrectHeader {
            found: header.to_owned(),
        }));
    }
    let mut rows = Vec::new();
    for (index, line) in lines.enumerate() {
        if line.is_empty() {
            return Err(rejected(PopulationSnvsFormatError::IncorrectColumnCount {
                line: index + 2,
                found: 0,
            }));
        }
        let fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
        let found = fields.len();
        let fields = fields.try_into().map_err(|_| {
            rejected(PopulationSnvsFormatError::IncorrectColumnCount {
                line: index + 2,
                found,
            })
        })?;
        rows.push(PopulationSnvRow { fields });
    }
    Ok(PopulationSnvsTsv { rows })
}

/// Analiza la tabla mínima suministrada (seis columnas y `source_id` opcional)
/// sin confundirla con el recurso canónico publicado de ocho columnas.
pub fn parse_supplied_tsv(bytes: &[u8]) -> Result<PopulationSnvsTsv, PopulationSnvsFormatError> {
    let text =
        std::str::from_utf8(bytes).map_err(|_| rejected(PopulationSnvsFormatError::InvalidUtf8))?;
    if text.starts_with("##fileformat=VCF") || text.starts_with("#CHROM\t") {
        return Err(rejected(PopulationSnvsFormatError::VcfInput));
    }
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| rejected(PopulationSnvsFormatError::MissingHeader))?;
    let accepted = [
        "contig\tposition\tref\talt\tpopulation_id\talt_af",
        "contig\tposition\tref\talt\tpopulation_id\talt_af\tsource_id",
    ];
    if !accepted.contains(&header) {
        return Err(rejected(PopulationSnvsFormatError::IncorrectHeader {
            found: header.to_owned(),
        }));
    }
    let has_source = header.ends_with("\tsource_id");
    let mut rows = Vec::new();
    for (index, line) in lines.enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != if has_source { 7 } else { 6 } {
            return Err(rejected(PopulationSnvsFormatError::IncorrectColumnCount {
                line: index + 2,
                found: fields.len(),
            }));
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
        .map_err(|_| rejected(PopulationSnvsFormatError::InvalidMetadataJson))?;
    let object = value
        .as_object()
        .ok_or_else(|| rejected(PopulationSnvsFormatError::MetadataMustBeObject))?;
    closed_object(
        object,
        &[
            "document_kind",
            "schema_version",
            "coordinate_system",
            "reference",
            "reference_fasta",
            "matrix",
            "populations",
            "provenance",
            "discarded_variant_pairs",
            "normalization_audit",
        ],
        "",
    )?;
    required_literal(object, "document_kind", "population_snvs_metadata")?;
    required_literal(object, "schema_version", "1")?;
    required_literal(object, "coordinate_system", "base1_closed")?;
    let reference = required_object(object, "reference")?;
    closed_object(
        reference,
        &["reference_id", "canonical_fingerprint"],
        "reference",
    )?;
    required_string(reference, "reference.reference_id")?;
    required_string(reference, "reference.canonical_fingerprint")?;
    let reference_fasta = required_object(object, "reference_fasta")?;
    closed_object(reference_fasta, &["sha256"], "reference_fasta")?;
    required_string(reference_fasta, "reference_fasta.sha256")?;
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
    let discarded = object
        .get("discarded_variant_pairs")
        .ok_or_else(|| missing("discarded_variant_pairs"))?
        .as_array()
        .ok_or_else(|| invalid("discarded_variant_pairs", "una lista"))?;
    let mut previous = None;
    for value in discarded {
        let pair = value
            .as_object()
            .ok_or_else(|| invalid("discarded_variant_pairs[]", "un objeto"))?;
        closed_object(
            pair,
            &["population_id", "variant_id"],
            "discarded_variant_pairs[]",
        )?;
        let population_id = pair
            .get("population_id")
            .and_then(Value::as_str)
            .ok_or_else(|| invalid("discarded_variant_pairs[].population_id", "un PopulationId"))?;
        PopulationId::parse(population_id, "population_snvs_metadata")
            .map_err(|_| invalid("discarded_variant_pairs[].population_id", "un PopulationId"))?;
        let variant_id = pair
            .get("variant_id")
            .and_then(Value::as_str)
            .ok_or_else(|| invalid("discarded_variant_pairs[].variant_id", "un VariantId"))?;
        VariantId::parse(variant_id, "population_snvs_metadata")
            .map_err(|_| invalid("discarded_variant_pairs[].variant_id", "un VariantId"))?;
        let key = (population_id, variant_id);
        if previous.is_some_and(|last| last >= key) {
            return Err(invalid(
                "discarded_variant_pairs",
                "pares población×VariantId ordenados y sin duplicados",
            ));
        }
        previous = Some(key);
    }
    if !matches!(object.get("normalization_audit"), Some(Value::Array(_))) {
        return Err(invalid("normalization_audit", "una lista"));
    }
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
            return Err(rejected(
                PopulationSnvsFormatError::UnexpectedMetadataField { field },
            ));
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
        .ok_or_else(|| {
            rejected(PopulationSnvsFormatError::InvalidMetadataLiteral {
                field: field.to_owned(),
                expected,
            })
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
        return Err(rejected(
            PopulationSnvsFormatError::InvalidMetadataLiteral {
                field: "provenance".to_owned(),
                expected: "\"supplied\" o \"generated\"",
            },
        ));
    }
    Ok(())
}

#[track_caller]
fn missing(field: &str) -> PopulationSnvsFormatError {
    crate::guard_trace::record();
    PopulationSnvsFormatError::MissingMetadataField {
        field: field.to_owned(),
    }
}

#[track_caller]
fn invalid(field: &str, expected: &'static str) -> PopulationSnvsFormatError {
    crate::guard_trace::record();
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
    /// Se pidieron bases de una referencia validada sin conservarlas.
    BasesNotMaterialized,
    InvalidCount,
    InvalidPopulation(String),
    DuplicatePopulation(String),
    InvalidFrequencyInterval {
        population_id: String,
    },
    InvalidFrequencySpectrum {
        population_id: String,
    },
    InvalidRegion {
        contig: String,
    },
    OverlappingRegions {
        contig: String,
    },
    InvalidWeight,
    EmptyEligibleSpace,
    InsufficientEligibleSpace {
        requested: u64,
        available: u64,
    },
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
/// Referencia desde la que se genera una matriz nueva.
///
/// Las bases se prestan cuando ya están en memoria: copiarlas duplicaría el
/// genoma y marcaría el pico de la ejecución.
struct GenerationReference<'a> {
    id: ReferenceId,
    prepared_identity: PreparedReferenceIdentity,
    sequences: Vec<(SequenceId, Cow<'a, [u8]>)>,
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
    let prepared_identity = PreparedReferenceIdentity::from_authoritative_fasta(
        reference.reference_id().clone(),
        &fasta,
    );
    let mut ordered = fasta.sequences.into_iter().collect::<Vec<_>>();
    ordered.sort_by_key(|(id, _)| sequence_order.get(id).copied());
    generate_population_snvs_view(
        config,
        GenerationReference {
            id: reference.reference_id().clone(),
            prepared_identity,
            sequences: ordered
                .into_iter()
                .map(|(id, bases)| {
                    Ok((
                        SequenceId::parse(id, "generated_population_snvs").map_err(|_| {
                            PopulationSnvsGenerationError::Reference("invalid sequence".into())
                        })?,
                        Cow::Owned(bases.into_owned()),
                    ))
                })
                .collect::<Result<_, _>>()?,
        },
        root_seed,
        1e-9,
    )
}

/// Misma operación en el recorrido directo, sobre la referencia ya validada en
/// memoria. No publica la matriz ni sustituye el par canónico de F2.2.7.
pub fn generate_population_snvs_from_reference(
    config: &PopulationSnvsGenerationConfig,
    reference: &Reference,
    root_seed: RootSeed,
) -> Result<PreparedPopulationSnvs, PopulationSnvsGenerationError> {
    generate_population_snvs_from_reference_with_tolerance(config, reference, root_seed, 1e-9)
}

/// Variante del recorrido directo que registra la misma tolerancia efectiva
/// que se publicará con el resultado. Mantiene la API histórica para los
/// consumidores que usan la tolerancia por defecto de F2.2.
pub fn generate_population_snvs_from_reference_with_tolerance(
    config: &PopulationSnvsGenerationConfig,
    reference: &Reference,
    root_seed: RootSeed,
    normalization_tolerance: f64,
) -> Result<PreparedPopulationSnvs, PopulationSnvsGenerationError> {
    generate_population_snvs_view(
        config,
        GenerationReference {
            id: reference.id().clone(),
            prepared_identity: PreparedReferenceIdentity::from_reference(reference)
                .map_err(|_| PopulationSnvsGenerationError::BasesNotMaterialized)?,
            sequences: reference
                .sequences()
                .iter()
                .map(|sequence| {
                    sequence
                        .sequence()
                        .as_bytes()
                        .map(|bases| (sequence.id().clone(), Cow::Borrowed(bases)))
                        .map_err(|_| PopulationSnvsGenerationError::BasesNotMaterialized)
                })
                .collect::<Result<_, _>>()?,
        },
        root_seed,
        normalization_tolerance,
    )
}

fn generate_population_snvs_view(
    config: &PopulationSnvsGenerationConfig,
    reference: GenerationReference,
    root_seed: RootSeed,
    normalization_tolerance: f64,
) -> Result<PreparedPopulationSnvs, PopulationSnvsGenerationError> {
    validate_generation_config(config)?;
    validate_generation_regions(config, &reference.sequences)?;
    let population_ids = config
        .populations
        .iter()
        .map(|population| {
            PopulationId::parse(&population.population_id, "generated_population_snvs").map_err(
                |_| {
                    PopulationSnvsGenerationError::InvalidPopulation(
                        population.population_id.clone(),
                    )
                },
            )
        })
        .collect::<Result<BTreeSet<_>, _>>()?
        .into_iter()
        .collect();
    let order = reference
        .sequences
        .iter()
        .enumerate()
        .map(|(index, (id, _))| (id.as_str().to_owned(), index))
        .collect::<BTreeMap<_, _>>();
    let available = generated_locus_available_count(config, &reference.sequences)?;
    if available == 0 {
        return Err(PopulationSnvsGenerationError::EmptyEligibleSpace);
    }
    let requested = match config.count {
        GeneratedLocusCount::Exact(value) => value,
        GeneratedLocusCount::PoissonMean(mean) => {
            sample_candidate_count(root_seed, &reference.id, mean)?
        }
        GeneratedLocusCount::PoissonDensityPerBase(density) => {
            sample_candidate_count(root_seed, &reference.id, density * available as f64)?
        }
    };
    if requested > available {
        return Err(PopulationSnvsGenerationError::InsufficientEligibleSpace {
            requested,
            available,
        });
    }
    let requested_usize =
        usize::try_from(requested).map_err(|_| PopulationSnvsGenerationError::InvalidCount)?;
    let selected = select_generated_loci(
        config,
        &reference.sequences,
        root_seed,
        reference.id.as_str(),
        requested_usize,
    )?;
    let mut loci = Vec::with_capacity(selected.len());
    for selected in selected {
        let (contig, _) = &reference.sequences[selected.sequence_index];
        let position = selected.position_1based;
        let reference_base = selected.reference_base;
        let mut alt_rng = derive_subseed(
            root_seed,
            RuntimeSubseedScope::MatrixAlternateAllele {
                scenario_id: "direct",
                reference_id: reference.id.as_str(),
                sequence_id: contig.as_str(),
                position_1based: position,
                reference_base: std::str::from_utf8(&[reference_base]).unwrap_or("A"),
            },
        )
        .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?
        .local_rng();
        let alternatives: Vec<u8> = (*b"ACGT")
            .into_iter()
            .filter(|base| base != &reference_base)
            .collect();
        let alternate = alternatives[categorical_f64(
            &mut alt_rng,
            &substitution_weights(config.substitution_profile, reference_base, &alternatives),
        )
        .map_err(PopulationSnvsGenerationError::Random)?];
        let variant_id =
            derive_matrix_variant_id(&reference.id, contig, position, reference_base, alternate)
                .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?;
        loci.push((
            contig.clone(),
            position,
            reference_base,
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
    let mut prepared = PreparedPopulationSnvs {
        prepared_reference: reference.prepared_identity,
        rows,
        population_ids,
        discarded_globally_absent: Vec::new(),
        empty_authorized: false,
        normalization_audit: Vec::new(),
    };
    let generated_table = parse_canonical_tsv(&canonical_tsv_bytes(&prepared))
        .map_err(|error| PopulationSnvsGenerationError::Reference(error.to_string()))?;
    prepared.normalization_audit =
        normalize_alt_frequencies(&generated_table, normalization_tolerance)
            .map_err(|error| PopulationSnvsGenerationError::Reference(error.to_string()))?
            .groups;
    Ok(prepared)
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

/// Candidato retenido por el selector de loci generado.
///
/// El orden invierte el índice canónico en los empates para que, al usar
/// `Reverse` como montón mínimo, se descarte primero el candidato posterior.
/// Los empates de puntuación son extraordinarios, pero esta regla mantiene el
/// resultado definido incluso entonces.
#[derive(Debug, PartialEq)]
struct SelectedGeneratedLocus {
    selection_score: f64,
    canonical_index: u64,
    sequence_index: usize,
    position_1based: u64,
    reference_base: u8,
}

impl Eq for SelectedGeneratedLocus {}

impl Ord for SelectedGeneratedLocus {
    fn cmp(&self, other: &Self) -> Ordering {
        self.selection_score
            .total_cmp(&other.selection_score)
            .then_with(|| other.canonical_index.cmp(&self.canonical_index))
    }
}

impl PartialOrd for SelectedGeneratedLocus {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn generated_locus_available_count(
    config: &PopulationSnvsGenerationConfig,
    sequences: &[(SequenceId, Cow<'_, [u8]>)],
) -> Result<u64, PopulationSnvsGenerationError> {
    let mut available = 0_u64;
    for (id, bases) in sequences {
        for (offset, base) in bases.iter().copied().enumerate() {
            if !matches!(base, b'A' | b'C' | b'G' | b'T') {
                continue;
            }
            let position =
                u64::try_from(offset).map_err(|_| PopulationSnvsGenerationError::InvalidCount)? + 1;
            if candidate_weight(&config.spatial, id.as_str(), position)? > 0.0 {
                available = available
                    .checked_add(1)
                    .ok_or(PopulationSnvsGenerationError::InvalidCount)?;
            }
        }
    }
    Ok(available)
}

/// Selecciona una muestra ponderada sin reemplazo sin conservar el espacio
/// completo de candidatos. Gumbel-top-k induce la misma distribución de
/// muestreo secuencial proporcional a los pesos que el selector histórico.
///
/// A diferencia del selector histórico, esta implementación consume una
/// palabra aleatoria por candidato de peso positivo (en vez de una por locus
/// seleccionado). Se conserva el mismo ámbito de subsemilla; por ello la
/// semilla sigue determinando el resultado, pero las salidas pueden diferir de
/// versiones anteriores para una misma semilla.
fn select_generated_loci(
    config: &PopulationSnvsGenerationConfig,
    sequences: &[(SequenceId, Cow<'_, [u8]>)],
    root_seed: RootSeed,
    reference_id: &str,
    requested: usize,
) -> Result<Vec<SelectedGeneratedLocus>, PopulationSnvsGenerationError> {
    if requested == 0 {
        return Ok(Vec::new());
    }
    let mut rng = derive_subseed(
        root_seed,
        RuntimeSubseedScope::MatrixLocusSelection {
            scenario_id: "direct",
            reference_id,
        },
    )
    .map_err(|error| PopulationSnvsGenerationError::Reference(error.render_human()))?
    .local_rng();
    let mut selected = BinaryHeap::with_capacity(requested);
    let mut canonical_index = 0_u64;

    for (sequence_index, (id, bases)) in sequences.iter().enumerate() {
        for (offset, reference_base) in bases.iter().copied().enumerate() {
            if !matches!(reference_base, b'A' | b'C' | b'G' | b'T') {
                continue;
            }
            let position_1based =
                u64::try_from(offset).map_err(|_| PopulationSnvsGenerationError::InvalidCount)? + 1;
            let weight = candidate_weight(&config.spatial, id.as_str(), position_1based)?;
            let candidate_index = canonical_index;
            canonical_index = canonical_index
                .checked_add(1)
                .ok_or(PopulationSnvsGenerationError::InvalidCount)?;
            if weight == 0.0 {
                continue;
            }
            let candidate = SelectedGeneratedLocus {
                selection_score: libm::log(weight) - libm::log(-libm::log(open_unit_f64(&mut rng))),
                canonical_index: candidate_index,
                sequence_index,
                position_1based,
                reference_base,
            };
            if selected.len() < requested {
                selected.push(Reverse(candidate));
            } else if let Some(lowest) = selected.peek()
                && candidate > lowest.0
            {
                selected.pop();
                selected.push(Reverse(candidate));
            }
        }
    }
    Ok(selected
        .into_iter()
        .map(|Reverse(candidate)| candidate)
        .collect())
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
    sequences: &[(SequenceId, Cow<'_, [u8]>)],
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

#[cfg(test)]
mod guard_tests {
    //! Pruebas de las guardas del juez de la matriz de SNV (DT-001, tanda 3).
    //!
    //! Cada caso parte de **una matriz y unos metadatos reales** (preparados con
    //! el código de producción) y los altera **en un solo punto**, de modo que el
    //! único rechazo posible es el de la guarda que se quiere probar; la prueba
    //! exige el motivo exacto (la variante y su `detail`). Las comprobaciones que
    //! ningún dato publicado puede alcanzar quedan anotadas junto a su guarda.
    //!
    //! No se usan aquí `rejected`, `audit_error`, `missing` ni `invalid`: los
    //! errores esperados se escriben con la variante literal, para que la prueba
    //! no dependa de la constructora que vigila y para no anotar disparos falsos.

    use super::*;
    use crate::reference::{DnaSequence, FastaHeader, ReferenceOrigin, ReferenceSequence};
    use crate::reference::{Reference as DomainReference, SequenceClass};

    const TOLERANCE: f64 = 1e-9;
    const FASTA_SHA: &str = "abababababababababababababababababababababababababababababababab";

    /// Seis columnas, dos poblaciones y cinco loci:
    /// - `chr1:1` A con C y G: el resto de referencia (suma menor que 1);
    /// - `chr1:2` C>T: en P2 vale cero pero la fila sigue activa porque P1 no es cero;
    /// - `chr2:1` G>A: cero en ambas, el locus entero se descarta;
    /// - `chr2:2` G con A y C: P1 suma 1.0000000005 (se normaliza) y P2 suma 1 exacto;
    /// - `chr2:4` A con C y G: C vale cero en ambas (se descarta) y G queda activa.
    const SUPPLIED: &str = "contig\tposition\tref\talt\tpopulation_id\talt_af
chr1\t1\tA\tC\tP1\t0.3
chr1\t1\tA\tG\tP1\t0.2
chr1\t1\tA\tC\tP2\t0.1
chr1\t1\tA\tG\tP2\t0.1
chr1\t2\tC\tT\tP1\t0.5
chr1\t2\tC\tT\tP2\t0.0
chr2\t1\tG\tA\tP1\t0
chr2\t1\tG\tA\tP2\t0
chr2\t2\tG\tA\tP1\t0.5000000005
chr2\t2\tG\tC\tP1\t0.5
chr2\t2\tG\tA\tP2\t0.5
chr2\t2\tG\tC\tP2\t0.5
chr2\t4\tA\tC\tP1\t0
chr2\t4\tA\tG\tP1\t0.2
chr2\t4\tA\tC\tP2\t0
chr2\t4\tA\tG\tP2\t0.1
";

    fn build_reference(materialized: bool) -> DomainReference {
        let sequences = [
            ("chr1", b"ACGTACGTAC".as_slice()),
            ("chr2", b"GGCATTCAGT".as_slice()),
        ]
        .iter()
        .map(|(id, bases)| {
            let id = SequenceId::parse(*id, "test").expect("sequence id");
            let sequence = if materialized {
                DnaSequence::try_from_bytes(bases).expect("DNA")
            } else {
                DnaSequence::summary(u64::try_from(bases.len()).expect("length"), 0)
                    .expect("summary")
            };
            ReferenceSequence::new(
                id.clone(),
                FastaHeader::parse(id.as_str().as_bytes()).expect("header"),
                SequenceClass::Chromosome,
                &sequence,
            )
            .expect("reference sequence")
        })
        .collect();
        DomainReference::new(
            ReferenceId::parse("ref-a", "test").expect("reference id"),
            ReferenceOrigin::Generated,
            sequences,
        )
        .expect("reference")
    }

    struct Pair {
        reference: DomainReference,
        prepared: PreparedPopulationSnvs,
        tsv: Vec<u8>,
        metadata: Value,
    }

    fn pair() -> Pair {
        let reference = build_reference(true);
        let table = parse_supplied_tsv(SUPPLIED.as_bytes()).expect("supplied table");
        let prepared = prepare_population_snvs_from_reference(&table, &reference, TOLERANCE, false)
            .expect("prepared");
        let tsv = canonical_tsv_bytes(&prepared);
        let metadata = canonical_metadata_document(
            &prepared,
            FASTA_SHA,
            &sha256_bytes(&tsv),
            PopulationSnvsProvenance::Supplied,
        );
        Pair {
            reference,
            prepared,
            tsv,
            metadata,
        }
    }

    fn validate_bytes(
        pair: &Pair,
        metadata: &[u8],
        tsv: &[u8],
    ) -> Result<(), PopulationSnvsPreparationError> {
        validate_reusable_population_snvs_bytes(
            metadata,
            tsv,
            &pair.reference,
            FASTA_SHA,
            TOLERANCE,
            false,
        )
        .map(|_| ())
    }

    fn validate_metadata(
        pair: &Pair,
        metadata: &Value,
    ) -> Result<(), PopulationSnvsPreparationError> {
        validate_bytes(
            pair,
            &serde_json::to_vec(metadata).expect("json"),
            &pair.tsv,
        )
    }

    fn declared_pairs(metadata: &Value) -> BTreeSet<(String, String)> {
        metadata["discarded_variant_pairs"]
            .as_array()
            .expect("pairs")
            .iter()
            .map(|pair| {
                (
                    pair["population_id"]
                        .as_str()
                        .expect("population")
                        .to_owned(),
                    pair["variant_id"].as_str().expect("variant").to_owned(),
                )
            })
            .collect()
    }

    fn check_audit(
        pair: &Pair,
        audit: &Value,
        tolerance: f64,
    ) -> Result<(), PopulationSnvsPreparationError> {
        validate_normalization_audit_document(
            audit,
            &pair.prepared,
            &pair.reference,
            tolerance,
            &declared_pairs(&pair.metadata),
        )
    }

    fn fasta_detail(detail: &str) -> PopulationSnvsPreparationError {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: detail.to_owned(),
        }
    }

    fn audit_reason(reason: &str) -> PopulationSnvsPreparationError {
        PopulationSnvsPreparationError::AuthoritativeFasta {
            detail: format!("auditoría de normalización inválida: {reason}"),
        }
    }

    fn group<'a>(
        audit: &'a mut Value,
        population: &str,
        contig: &str,
        position: &str,
    ) -> &'a mut Value {
        audit
            .as_array_mut()
            .expect("audit list")
            .iter_mut()
            .find(|group| {
                group["population_id"] == population
                    && group["contig"] == contig
                    && group["position"] == position
            })
            .expect("group")
    }

    fn group_index(audit: &Value, population: &str, contig: &str, position: &str) -> usize {
        audit
            .as_array()
            .expect("audit list")
            .iter()
            .position(|group| {
                group["population_id"] == population
                    && group["contig"] == contig
                    && group["position"] == position
            })
            .expect("group")
    }

    fn next_up(value: f64) -> f64 {
        f64::from_bits(value.to_bits() + 1)
    }

    /// Recalcula `original_sum` y `original_delta` con la misma suma que el
    /// juez, para alterar un `original` sin que lo cace antes `sumas no coinciden`.
    fn refresh_original_sums(group: &mut Value) {
        let originals = group["values"]
            .as_array()
            .expect("values")
            .iter()
            .map(|value| value["original"].as_f64().expect("original"))
            .collect::<Vec<_>>();
        let sum = neumaier(originals);
        group["original_sum"] = json!(sum);
        group["original_delta"] = json!((sum - 1.0).abs());
    }

    type Edit = Box<dyn Fn(&mut Value)>;

    /// Lo que las demás pruebas dan por supuesto: la pareja sin alterar es válida
    /// por las tres entradas (bytes, auditoría directa y restauración).
    #[test]
    fn the_unaltered_pair_is_accepted_and_exercises_every_decision() {
        let pair = pair();
        validate_metadata(&pair, &pair.metadata).expect("bytes entry");
        check_audit(&pair, &pair.metadata["normalization_audit"], TOLERANCE).expect("audit entry");
        let mut prepared = pair.prepared.clone();
        restore_published_allele_frequencies(&mut prepared, &pair.metadata["normalization_audit"])
            .expect("restore entry");
        assert_eq!(prepared, pair.prepared);
        let audit = pair.metadata["normalization_audit"]
            .as_array()
            .expect("list");
        let decisions = audit
            .iter()
            .map(|group| group["decision"].as_str().expect("decision"))
            .collect::<BTreeSet<_>>();
        assert_eq!(
            decisions,
            BTreeSet::from(["normalized", "reference_remainder"])
        );
        assert_eq!(
            pair.metadata["discarded_variant_pairs"]
                .as_array()
                .unwrap()
                .len(),
            4
        );
        assert!(
            pair.tsv.starts_with(b"variant_id\tcontig"),
            "the canonical matrix is the one under test"
        );
    }

    // ----------------------------------------------------------------------
    // `restore_published_allele_frequencies`, llamada directa sobre la
    // auditoría alterada en un solo punto. El juez de bytes la corre antes que
    // `validate_normalization_audit_document`, que repite varios motivos con el
    // mismo texto: por eso estas pruebas la llaman sola, y así una guarda
    // borrada no queda tapada por la segunda comprobación.
    // ----------------------------------------------------------------------

    #[test]
    fn restore_rejects_each_malformed_published_audit_point_for_its_own_reason() {
        let pair = pair();
        let cases: Vec<(&str, Edit, &str)> = vec![
            (
                "group is not an object",
                Box::new(|audit| audit[0] = json!(7)),
                "grupo inválido",
            ),
            (
                "population_id missing",
                Box::new(|audit| {
                    audit[0].as_object_mut().unwrap().remove("population_id");
                }),
                "population_id",
            ),
            (
                "population_id not a string",
                Box::new(|audit| audit[0]["population_id"] = json!(5)),
                "population_id",
            ),
            (
                "contig missing",
                Box::new(|audit| {
                    audit[0].as_object_mut().unwrap().remove("contig");
                }),
                "contig",
            ),
            (
                "position not a string",
                Box::new(|audit| audit[0]["position"] = json!(1)),
                "position",
            ),
            (
                "ref missing",
                Box::new(|audit| {
                    audit[0].as_object_mut().unwrap().remove("ref");
                }),
                "ref",
            ),
            (
                "reference_frequency missing",
                Box::new(|audit| {
                    audit[0]
                        .as_object_mut()
                        .unwrap()
                        .remove("reference_frequency");
                }),
                "reference_frequency",
            ),
            (
                "reference_frequency not a number",
                Box::new(|audit| audit[0]["reference_frequency"] = json!("0.5")),
                "reference_frequency",
            ),
            (
                "reference_frequency negative",
                Box::new(|audit| audit[0]["reference_frequency"] = json!(-0.25)),
                "reference_frequency",
            ),
            (
                "values missing",
                Box::new(|audit| {
                    audit[0].as_object_mut().unwrap().remove("values");
                }),
                "values",
            ),
            (
                "values not a list",
                Box::new(|audit| audit[0]["values"] = json!({})),
                "values",
            ),
            (
                "value is not an object",
                Box::new(|audit| audit[0]["values"][0] = json!("x")),
                "valor inválido",
            ),
            (
                "alt missing",
                Box::new(|audit| {
                    audit[0]["values"][0].as_object_mut().unwrap().remove("alt");
                }),
                "alt",
            ),
            (
                "variant_id not a string",
                Box::new(|audit| audit[0]["values"][0]["variant_id"] = json!(3)),
                "variant_id",
            ),
            (
                "source_id missing",
                Box::new(|audit| {
                    audit[0]["values"][0]
                        .as_object_mut()
                        .unwrap()
                        .remove("source_id");
                }),
                "source_id",
            ),
            (
                "effective missing",
                Box::new(|audit| {
                    audit[0]["values"][0]
                        .as_object_mut()
                        .unwrap()
                        .remove("effective");
                }),
                "effective",
            ),
            (
                "effective negative",
                Box::new(|audit| audit[0]["values"][0]["effective"] = json!(-0.125)),
                "effective",
            ),
            (
                "duplicate group",
                Box::new(|audit| {
                    let copy = audit[0].clone();
                    audit.as_array_mut().unwrap().push(copy);
                }),
                "grupo duplicado",
            ),
            (
                "duplicate value inside a group",
                Box::new(|audit| {
                    let copy = audit[0]["values"][0].clone();
                    audit[0]["values"].as_array_mut().unwrap().push(copy);
                }),
                "valor duplicado",
            ),
        ];
        for (name, edit, reason) in cases {
            let mut audit = pair.metadata["normalization_audit"].clone();
            edit(&mut audit);
            assert_ne!(
                audit, pair.metadata["normalization_audit"],
                "{name}: no-op edit"
            );
            let mut prepared = pair.prepared.clone();
            assert_eq!(
                restore_published_allele_frequencies(&mut prepared, &audit),
                Err(audit_reason(reason)),
                "{name}"
            );
        }
    }

    /// Las tres negativas de la segunda mitad de la restauración: una fila
    /// activa sin su valor efectivo, un grupo de la matriz que la auditoría no
    /// trae y un valor de un grupo que sí trae pero sin ese valor.
    #[test]
    fn restore_rejects_an_audit_that_lacks_what_the_matrix_needs() {
        let pair = pair();
        let base = pair.metadata["normalization_audit"].clone();

        let mut audit = base.clone();
        audit[0]["values"].as_array_mut().unwrap().remove(0);
        // El primer grupo del documento publicado es P1 chr1:1, que tiene dos
        // valores activos: quitar uno deja una fila activa sin su efectivo.
        assert_eq!(audit[0]["population_id"], "P1");
        assert_eq!(audit[0]["position"], "1");
        assert_eq!(
            restore_published_allele_frequencies(&mut pair.prepared.clone(), &audit),
            Err(audit_reason("efectivo activo ausente"))
        );

        // El locus `chr2:1` está descartado entero: no tiene filas activas, de
        // modo que solo lo echa de menos el bucle de grupos.
        let mut audit = base.clone();
        let index = group_index(&audit, "P1", "chr2", "1");
        audit.as_array_mut().unwrap().remove(index);
        assert_eq!(
            restore_published_allele_frequencies(&mut pair.prepared.clone(), &audit),
            Err(audit_reason("referencia ausente"))
        );

        // En `chr2:4` la C está descartada y la G activa: quitar la C deja el
        // grupo, deja la fila activa y solo falta el valor del descartado.
        let mut audit = base;
        let discarded = group(&mut audit, "P1", "chr2", "4");
        let values = discarded["values"].as_array_mut().unwrap();
        assert_eq!(values[0]["alt"], "C");
        values.remove(0);
        assert_eq!(
            restore_published_allele_frequencies(&mut pair.prepared.clone(), &audit),
            Err(audit_reason("efectivo ausente"))
        );
    }

    /// Un locus multialélico normalizado conserva REF a cero tras releerlo del
    /// JSON: la restauración reescribe la frecuencia de las filas y del grupo.
    #[test]
    fn restore_overwrites_rows_and_groups_with_the_published_values() {
        let pair = pair();
        let mut audit = pair.metadata["normalization_audit"].clone();
        let tweaked = next_up(0.3);
        let first = group(&mut audit, "P1", "chr1", "1");
        first["values"][0]["effective"] = json!(tweaked);
        first["reference_frequency"] = json!(0.123);
        let mut prepared = pair.prepared.clone();
        restore_published_allele_frequencies(&mut prepared, &audit).expect("restored");
        let row = prepared
            .rows()
            .iter()
            .find(|row| {
                row.population_id().as_str() == "P1"
                    && row.contig().as_str() == "chr1"
                    && row.position_1based() == 1
                    && row.alternate() == b'C'
            })
            .expect("row");
        assert_eq!(row.alt_af().to_bits(), tweaked.to_bits());
        let restored = prepared
            .normalization_audit()
            .iter()
            .find(|group| {
                group.population_id() == "P1" && group.contig() == "chr1" && group.position() == "1"
            })
            .expect("group");
        assert_eq!(
            restored.reference_frequency().to_bits(),
            0.123_f64.to_bits()
        );
        assert_eq!(
            restored.values()[0].effective().to_bits(),
            tweaked.to_bits()
        );
    }

    // ----------------------------------------------------------------------
    // `validate_normalization_audit_document`
    // ----------------------------------------------------------------------

    #[test]
    fn the_audit_judge_rejects_each_malformed_group_field_for_its_own_reason() {
        let pair = pair();
        let remove = |field: &'static str| -> Edit {
            Box::new(move |audit| {
                audit[0].as_object_mut().unwrap().remove(field);
            })
        };
        let set = |field: &'static str, value: Value| -> Edit {
            Box::new(move |audit| audit[0][field] = value.clone())
        };
        let cases: Vec<(&str, Edit, &str)> = vec![
            (
                "audit is not a list",
                Box::new(|audit| *audit = json!({})),
                "debe ser una lista",
            ),
            (
                "group is not an object",
                Box::new(|audit| audit[0] = json!(1)),
                "grupo inválido",
            ),
            (
                "one field missing, none extra",
                remove("tolerance"),
                "campos de grupo inválidos",
            ),
            (
                "one extra field",
                Box::new(|audit| audit[0]["extra"] = json!(1)),
                "campos de grupo inválidos",
            ),
            (
                "one field renamed, same count",
                Box::new(|audit| {
                    let value = audit[0]
                        .as_object_mut()
                        .unwrap()
                        .remove("tolerance")
                        .unwrap();
                    audit[0]["tolerancia"] = value;
                }),
                "campos de grupo inválidos",
            ),
            // Cada llamante de la clausura `text` y de la clausura `number`
            // cuenta como una guarda para la puerta, así que lleva su caso.
            (
                "population_id as number",
                set("population_id", json!(5)),
                "population_id",
            ),
            ("contig as number", set("contig", json!(5)), "contig"),
            ("position as number", set("position", json!(1)), "position"),
            ("ref as number", set("ref", json!(1)), "ref"),
            (
                "tolerance as string",
                set("tolerance", json!("1e-9")),
                "tolerance",
            ),
            (
                "canonical_order as number",
                set("canonical_order", json!(1)),
                "canonical_order",
            ),
            (
                "sum_algorithm as number",
                set("sum_algorithm", json!(1)),
                "sum_algorithm",
            ),
            (
                "original_sum as string",
                set("original_sum", json!("0.5")),
                "original_sum",
            ),
            (
                "original_delta missing",
                remove("original_delta"),
                "campos de grupo inválidos",
            ),
            (
                "original_delta as string",
                set("original_delta", json!("x")),
                "original_delta",
            ),
            (
                "effective_sum as string",
                set("effective_sum", json!("x")),
                "effective_sum",
            ),
            (
                "effective_delta as string",
                set("effective_delta", json!("x")),
                "effective_delta",
            ),
            ("decision as number", set("decision", json!(5)), "decision"),
            (
                "reference_frequency as string, reference_remainder",
                set("reference_frequency", json!("x")),
                "reference_frequency",
            ),
            ("values as object", set("values", json!({})), "values"),
            ("values empty", set("values", json!([])), "values vacíos"),
        ];
        for (name, edit, reason) in cases {
            let mut audit = pair.metadata["normalization_audit"].clone();
            assert_eq!(audit[0]["decision"], "reference_remainder", "case group");
            edit(&mut audit);
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason(reason)),
                "{name}"
            );
        }

        // `reference_frequency` también la lee la rama `normalized`.
        let mut audit = pair.metadata["normalization_audit"].clone();
        let normalized = group(&mut audit, "P1", "chr2", "2");
        assert_eq!(normalized["decision"], "normalized");
        normalized["reference_frequency"] = json!("x");
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("reference_frequency")),
            "reference_frequency as string, normalized"
        );
    }

    #[test]
    fn the_audit_judge_rejects_an_identity_that_does_not_fit_the_matrix_or_the_reference() {
        let pair = pair();
        let check = |edit: &dyn Fn(&mut Value), reason: &str, name: &str| {
            let mut audit = pair.metadata["normalization_audit"].clone();
            edit(&mut audit);
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason(reason)),
                "{name}"
            );
        };
        check(
            &|a| a[0]["population_id"] = json!("bad id"),
            "population_id",
            "invalid population id",
        );
        check(
            &|a| a[0]["population_id"] = json!("P9"),
            "population_id no pertenece a la matriz",
            "population outside the matrix",
        );
        check(
            &|a| a[0]["position"] = json!("abc"),
            "position",
            "position not a number",
        );
        check(
            &|a| a[0]["position"] = json!("0"),
            "position",
            "position zero",
        );
        check(
            &|a| a[0]["position"] = json!("-1"),
            "position",
            "position negative",
        );
        check(
            &|a| a[0]["position"] = json!("08"),
            "position",
            "position with a leading zero",
        );
        check(
            &|a| a[0]["position"] = json!("+8"),
            "position",
            "position with a plus sign",
        );
        check(&|a| a[0]["ref"] = json!("N"), "ref", "ref not ACGT");
        check(&|a| a[0]["ref"] = json!("AC"), "ref", "ref two bases");
        check(
            &|a| a[0]["contig"] = json!("bad contig"),
            "contig",
            "invalid contig id",
        );
        check(
            &|a| a[0]["contig"] = json!("chr9"),
            "contig ausente de referencia",
            "contig outside the reference",
        );
        check(
            &|a| a[0]["ref"] = json!("C"),
            "ref no coincide con referencia",
            "ref differs from the reference base",
        );
        check(
            &|a| a[0]["position"] = json!("11"),
            "ref no coincide con referencia",
            "position past the end of the contig",
        );
    }

    /// Una referencia sin bases materializadas no puede juzgar la auditoría:
    /// es un error propio, no un «ref no coincide».
    #[test]
    fn the_audit_judge_reports_an_unmaterialized_reference_as_such() {
        let pair = pair();
        let unmaterialized = build_reference(false);
        let result = validate_normalization_audit_document(
            &pair.metadata["normalization_audit"],
            &pair.prepared,
            &unmaterialized,
            TOLERANCE,
            &declared_pairs(&pair.metadata),
        );
        assert_eq!(
            result,
            Err(PopulationSnvsPreparationError::BasesNotMaterialized)
        );
    }

    #[test]
    fn the_audit_judge_rejects_each_half_of_the_group_identity_and_order_check() {
        let pair = pair();
        let cases: Vec<(&str, Edit)> = vec![
            (
                "tolerance differs",
                Box::new(|a| a[0]["tolerance"] = json!(2e-9)),
            ),
            (
                "tolerance differs in the last bit",
                Box::new(|a| a[0]["tolerance"] = json!(next_up(TOLERANCE))),
            ),
            (
                "canonical_order differs",
                Box::new(|a| a[0]["canonical_order"] = json!("otro")),
            ),
            (
                "sum_algorithm differs",
                Box::new(|a| a[0]["sum_algorithm"] = json!("otro")),
            ),
            (
                "groups out of order",
                Box::new(|a| {
                    a.as_array_mut().unwrap().swap(0, 1);
                }),
            ),
            (
                "group repeated right after itself",
                Box::new(|a| {
                    let copy = a[0].clone();
                    a.as_array_mut().unwrap().insert(1, copy);
                }),
            ),
            (
                "group repeated at the end",
                Box::new(|a| {
                    let copy = a[0].clone();
                    a.as_array_mut().unwrap().push(copy);
                }),
            ),
        ];
        for (name, edit) in cases {
            let mut audit = pair.metadata["normalization_audit"].clone();
            edit(&mut audit);
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason("identidad o algoritmo de grupo inválido")),
                "{name}"
            );
        }
    }

    #[test]
    fn the_audit_judge_rejects_each_malformed_value_for_its_own_reason() {
        let pair = pair();
        // Premisa de la guarda `!value.is_empty()` (equivalente, anotada): un texto
        // vacío tampoco es un `VariantId`, así que lo rechazaría igual la línea siguiente.
        assert!(VariantId::parse("", "test").is_err());
        let value_edit = |edit: fn(&mut Value)| -> Edit {
            Box::new(move |audit| edit(&mut audit[0]["values"][0]))
        };
        let cases: Vec<(&str, Edit, &str)> = vec![
            (
                "value not an object",
                Box::new(|a| a[0]["values"][0] = json!(1)),
                "valor inválido",
            ),
            (
                "one value field missing",
                value_edit(|v| {
                    v.as_object_mut().unwrap().remove("source_id");
                }),
                "campos de valor inválidos",
            ),
            (
                "one extra value field",
                value_edit(|v| v["extra"] = json!(1)),
                "campos de valor inválidos",
            ),
            (
                "one value field renamed",
                value_edit(|v| {
                    let original = v.as_object_mut().unwrap().remove("original").unwrap();
                    v["origin"] = original;
                }),
                "campos de valor inválidos",
            ),
            ("alt as number", value_edit(|v| v["alt"] = json!(1)), "alt"),
            (
                "source_id as number",
                value_edit(|v| v["source_id"] = json!(1)),
                "source_id",
            ),
            (
                "variant_id as number",
                value_edit(|v| v["variant_id"] = json!(1)),
                "variant_id",
            ),
            (
                "variant_id empty",
                value_edit(|v| v["variant_id"] = json!("")),
                "variant_id",
            ),
            (
                "variant_id with an invalid character",
                value_edit(|v| v["variant_id"] = json!("bad id")),
                "variant_id",
            ),
            ("alt not ACGT", value_edit(|v| v["alt"] = json!("N")), "alt"),
            (
                "alt two bases",
                value_edit(|v| v["alt"] = json!("CC")),
                "alt",
            ),
            (
                "alt equal to ref",
                value_edit(|v| v["alt"] = json!("A")),
                "alt igual a ref",
            ),
            (
                "variant_id valid but not the derived one",
                value_edit(|v| v["variant_id"] = json!("VAR_OTRO")),
                "variant_id no canónico",
            ),
            (
                "original missing (the field set no longer fits)",
                value_edit(|v| {
                    v.as_object_mut().unwrap().remove("original");
                }),
                "campos de valor inválidos",
            ),
            (
                "original as string",
                value_edit(|v| v["original"] = json!("0.3")),
                "original",
            ),
            (
                "original negative",
                value_edit(|v| v["original"] = json!(-0.1)),
                "original",
            ),
            (
                "original above one",
                value_edit(|v| v["original"] = json!(1.5)),
                "original",
            ),
            (
                "effective as string",
                value_edit(|v| v["effective"] = json!("0.3")),
                "effective",
            ),
            (
                "effective negative",
                value_edit(|v| v["effective"] = json!(-0.1)),
                "effective",
            ),
            (
                "effective above one",
                value_edit(|v| v["effective"] = json!(1.5)),
                "effective",
            ),
        ];
        for (name, edit, reason) in cases {
            let mut audit = pair.metadata["normalization_audit"].clone();
            assert_eq!(audit[0]["values"][0]["alt"], "C", "case value");
            edit(&mut audit);
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason(reason)),
                "{name}"
            );
        }
    }

    #[test]
    fn the_audit_judge_rejects_values_out_of_canonical_order_or_repeated() {
        let pair = pair();
        let reason = audit_reason("orden canónico de values inválido");
        // P1 chr1:1 trae C y luego G.
        let mut audit = pair.metadata["normalization_audit"].clone();
        audit[0]["values"].as_array_mut().unwrap().swap(0, 1);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(reason.clone()),
            "descending alt"
        );

        let mut audit = pair.metadata["normalization_audit"].clone();
        let copy = audit[0]["values"][0].clone();
        audit[0]["values"].as_array_mut().unwrap().insert(1, copy);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(reason.clone()),
            "same value twice"
        );

        // Mismo alt y mismo variant_id, source_id descendente: el valor cero de
        // `chr2:4` (la C descartada) no está en la matriz y se acepta como ausente.
        let mut audit = pair.metadata["normalization_audit"].clone();
        let mixed = group(&mut audit, "P1", "chr2", "4");
        let discarded = mixed["values"][0].clone();
        assert_eq!(discarded["alt"], "C");
        let mut second = discarded.clone();
        second["source_id"] = json!("a");
        mixed["values"] = json!([discarded, second]);
        mixed["values"][0]["source_id"] = json!("b");
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(reason),
            "same alt and variant, descending source"
        );
    }

    #[test]
    fn the_audit_judge_compares_effective_values_with_the_matrix_bit_by_bit() {
        let pair = pair();
        let mut audit = pair.metadata["normalization_audit"].clone();
        let tweaked = next_up(0.3);
        audit[0]["values"][0]["effective"] = json!(tweaked);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("effective no coincide con TSV"))
        );

        // P2 chr1:2 vale cero y la fila sigue activa: `-0.0` es igual a `0.0` como
        // número pero no como bits, y el juez exige los bits.
        let mut audit = pair.metadata["normalization_audit"].clone();
        let zero = group(&mut audit, "P2", "chr1", "2");
        assert_eq!(
            zero["values"][0]["effective"].as_f64().unwrap().to_bits(),
            0
        );
        zero["values"][0]["effective"] = json!(-0.0_f64);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("effective no coincide con TSV"))
        );

        // Un valor que no está en la matriz solo puede valer cero.
        let mut audit = pair.metadata["normalization_audit"].clone();
        let mixed = group(&mut audit, "P1", "chr2", "4");
        assert_eq!(mixed["values"][0]["alt"], "C");
        mixed["values"][0]["effective"] = json!(0.1);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("valor efectivo activo ausente del TSV"))
        );

        // Una auditoría a la que le falta un grupo con filas activas.
        let mut audit = pair.metadata["normalization_audit"].clone();
        let index = group_index(&audit, "P1", "chr1", "2");
        audit.as_array_mut().unwrap().remove(index);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("faltan valores efectivos del TSV"))
        );
    }

    #[test]
    fn the_audit_judge_rejects_each_sum_that_does_not_match_its_values() {
        let pair = pair();
        for field in [
            "original_sum",
            "original_delta",
            "effective_sum",
            "effective_delta",
        ] {
            let mut audit = pair.metadata["normalization_audit"].clone();
            let current = audit[0][field].as_f64().expect("number");
            audit[0][field] = json!(next_up(current));
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason("sumas no coinciden")),
                "{field}"
            );
        }
    }

    #[test]
    fn the_audit_judge_rejects_a_decision_that_the_values_do_not_support() {
        let pair = pair();

        // --- reference_remainder (P1 chr1:1) ---
        let mut audit = pair.metadata["normalization_audit"].clone();
        let remainder = group(&mut audit, "P1", "chr1", "1");
        remainder["values"][0]["original"] = json!(0.25);
        refresh_original_sums(remainder);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("resto de referencia inválido")),
            "original differs from effective"
        );
        let mut audit = pair.metadata["normalization_audit"].clone();
        let remainder = group(&mut audit, "P1", "chr1", "1");
        let frequency = remainder["reference_frequency"].as_f64().unwrap();
        remainder["reference_frequency"] = json!(next_up(frequency));
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("resto de referencia inválido")),
            "reference frequency off by one bit"
        );

        // --- reference_remainder: los bits, no el valor, de `-0.0` ---
        // P2 chr2:2 suma exactamente 1 y publica REF = +0.0.
        let mut audit = pair.metadata["normalization_audit"].clone();
        let exact = group(&mut audit, "P2", "chr2", "2");
        assert_eq!(exact["decision"], "reference_remainder");
        assert_eq!(exact["reference_frequency"].as_f64().unwrap().to_bits(), 0);
        exact["reference_frequency"] = json!(-0.0_f64);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("resto de referencia inválido")),
            "REF -0.0 published for +0.0 (direct)"
        );
        let mut metadata = pair.metadata.clone();
        metadata["normalization_audit"] = audit;
        assert_eq!(
            validate_metadata(&pair, &metadata),
            Err(audit_reason("resto de referencia inválido")),
            "REF -0.0 published for +0.0 (bytes)"
        );
        // Un `original` de `-0.0` frente a un efectivo `+0.0` (P2 chr1:2 vale cero).
        let mut audit = pair.metadata["normalization_audit"].clone();
        let zero = group(&mut audit, "P2", "chr1", "2");
        assert_eq!(
            zero["values"][0]["effective"].as_f64().unwrap().to_bits(),
            0
        );
        zero["values"][0]["original"] = json!(-0.0_f64);
        refresh_original_sums(zero);
        assert_eq!(
            zero["original_sum"].as_f64().unwrap().to_bits(),
            0,
            "premise"
        );
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("resto de referencia inválido")),
            "original -0.0 against effective +0.0 (direct)"
        );
        let mut metadata = pair.metadata.clone();
        metadata["normalization_audit"] = audit;
        assert_eq!(
            validate_metadata(&pair, &metadata),
            Err(audit_reason("resto de referencia inválido")),
            "original -0.0 against effective +0.0 (bytes)"
        );

        // --- normalized (P1 chr2:2) ---
        let mut audit = pair.metadata["normalization_audit"].clone();
        let normalized = group(&mut audit, "P1", "chr2", "2");
        assert_eq!(normalized["decision"], "normalized");
        normalized["values"][0]["original"] = json!(0.5000000004);
        refresh_original_sums(normalized);
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("normalización inválida")),
            "effective is not original / sum"
        );
        for (name, frequency) in [("positive", 0.1_f64), ("negative zero", -0.0_f64)] {
            let mut audit = pair.metadata["normalization_audit"].clone();
            group(&mut audit, "P1", "chr2", "2")["reference_frequency"] = json!(frequency);
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason("normalización inválida")),
                "reference frequency {name}"
            );
        }

        // --- decisión inválida ---
        let decision_cases: Vec<(&str, &str, &str, &str, &str)> = vec![
            ("unknown decision", "P1", "1", "chr1", "otra"),
            (
                "remainder declared on a sum above one",
                "P1",
                "2",
                "chr2",
                "reference_remainder",
            ),
            (
                "normalized declared on a sum below one",
                "P1",
                "1",
                "chr1",
                "normalized",
            ),
            (
                "normalized declared on a sum of exactly one",
                "P2",
                "2",
                "chr2",
                "normalized",
            ),
        ];
        for (name, population, position, contig, decision) in decision_cases {
            let mut audit = pair.metadata["normalization_audit"].clone();
            group(&mut audit, population, contig, position)["decision"] = json!(decision);
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason("decisión inválida")),
                "{name}"
            );
        }
        // `normalized` exige además que el delta original no pase la tolerancia:
        // con la tolerancia publicada y la llamante rebajadas a 1e-12, la suma de
        // P1 chr2:2 (delta 5e-10) ya no se puede normalizar.
        let mut audit = pair.metadata["normalization_audit"].clone();
        for group_value in audit.as_array_mut().unwrap() {
            group_value["tolerance"] = json!(1e-12);
        }
        assert_eq!(
            check_audit(&pair, &audit, 1e-12),
            Err(audit_reason("decisión inválida")),
            "normalized with a delta above the tolerance"
        );
    }

    #[test]
    fn the_audit_judge_rejects_discarded_pairs_that_do_not_match_the_zero_values() {
        let pair = pair();
        let reason =
            audit_reason("ceros ausentes no coinciden con pares población×variante descartados");
        let audit = pair.metadata["normalization_audit"].clone();

        // Declarar menos de lo que la auditoría pone a cero.
        assert_eq!(
            validate_normalization_audit_document(
                &audit,
                &pair.prepared,
                &pair.reference,
                TOLERANCE,
                &BTreeSet::new()
            ),
            Err(reason.clone()),
            "nothing declared"
        );

        // La auditoría omite el grupo de P2 del locus descartado y la declaración
        // lo calla igual: coinciden entre sí (primera mitad) pero la matriz tiene
        // dos poblaciones y el par de P2 debería estar (segunda mitad).
        let mut shortened = audit;
        let index = group_index(&shortened, "P2", "chr2", "1");
        shortened.as_array_mut().unwrap().remove(index);
        let index = group_index(&shortened, "P2", "chr2", "4");
        shortened[index]["values"].as_array_mut().unwrap().remove(0);
        let mut declared = declared_pairs(&pair.metadata);
        declared.retain(|(population, _)| population == "P1");
        assert_eq!(declared.len(), 2);
        assert_eq!(
            validate_normalization_audit_document(
                &shortened,
                &pair.prepared,
                &pair.reference,
                TOLERANCE,
                &declared
            ),
            Err(reason),
            "declaration and audit agree but the matrix has another population"
        );
    }

    // ----------------------------------------------------------------------
    // `validate_reusable_population_snvs_bytes` y `..._snapshots`
    // ----------------------------------------------------------------------

    #[test]
    fn a_reused_pair_is_rejected_for_the_exact_binding_that_fails() {
        let pair = pair();
        let binding =
            fasta_detail("la matriz reutilizada no está ligada a la referencia autoritativa");
        let cases: Vec<(&str, Edit)> = vec![
            (
                "reference id",
                Box::new(|m| m["reference"]["reference_id"] = json!("otra")),
            ),
            (
                "canonical fingerprint",
                Box::new(|m| m["reference"]["canonical_fingerprint"] = json!("0".repeat(64))),
            ),
            (
                "fasta sha256",
                Box::new(|m| m["reference_fasta"]["sha256"] = json!("0".repeat(64))),
            ),
        ];
        for (name, edit) in cases {
            let mut metadata = pair.metadata.clone();
            edit(&mut metadata);
            assert_eq!(
                validate_metadata(&pair, &metadata),
                Err(binding.clone()),
                "{name}"
            );
        }

        let mut metadata = pair.metadata.clone();
        metadata["matrix"]["sha256"] = json!("0".repeat(64));
        assert_eq!(
            validate_metadata(&pair, &metadata),
            Err(fasta_detail(
                "SHA-256 del TSV reutilizado no coincide con metadata"
            ))
        );
    }

    #[test]
    fn a_reused_matrix_not_in_canonical_order_is_rejected_even_with_a_matching_digest() {
        let pair = pair();
        let text = String::from_utf8(pair.tsv.clone()).expect("utf8");
        let mut lines = text.lines().map(str::to_owned).collect::<Vec<_>>();
        lines.swap(1, 2);
        let reordered = format!("{}\n", lines.join("\n"));
        let mut metadata = pair.metadata.clone();
        metadata["matrix"]["sha256"] = json!(sha256_bytes(reordered.as_bytes()));
        assert_eq!(
            validate_bytes(
                &pair,
                &serde_json::to_vec(&metadata).unwrap(),
                reordered.as_bytes()
            ),
            Err(fasta_detail(
                "el TSV reutilizado no está en orden y representación canónicos"
            ))
        );
    }

    #[test]
    fn a_reused_pair_with_inconsistent_populations_is_rejected_for_the_exact_reason() {
        let pair = pair();
        let cases: Vec<(&str, Value, &str)> = vec![
            (
                "invalid population id",
                json!([{"population_id": "bad id"}, {"population_id": "P2"}]),
                "metadata.populations inválida",
            ),
            (
                "descending",
                json!([{"population_id": "P2"}, {"population_id": "P1"}]),
                "metadata.populations no está ordenada canónicamente",
            ),
            (
                "repeated",
                json!([{"population_id": "P1"}, {"population_id": "P1"}]),
                "metadata.populations no está ordenada canónicamente",
            ),
            (
                "another population",
                json!([{"population_id": "P1"}, {"population_id": "P3"}]),
                "metadata.populations no coincide con el TSV canónico",
            ),
            (
                "one population fewer",
                json!([{"population_id": "P1"}]),
                "metadata.populations no coincide con el TSV canónico",
            ),
        ];
        for (name, populations, detail) in cases {
            let mut metadata = pair.metadata.clone();
            metadata["populations"] = populations;
            assert_eq!(
                validate_metadata(&pair, &metadata),
                Err(fasta_detail(detail)),
                "{name}"
            );
        }
    }

    #[test]
    fn a_reused_pair_needs_a_materialized_reference() {
        let pair = pair();
        let unmaterialized = build_reference(false);
        let result = validate_reusable_population_snvs_bytes(
            &serde_json::to_vec(&pair.metadata).unwrap(),
            &pair.tsv,
            &unmaterialized,
            FASTA_SHA,
            TOLERANCE,
            false,
        );
        assert_eq!(
            result.map(|_| ()),
            Err(PopulationSnvsPreparationError::BasesNotMaterialized)
        );
    }

    /// El par reutilizado se lee de dos ficheros hermanos: si falta la matriz
    /// junto a los metadatos, el motivo nombra la ruta que se buscó.
    #[test]
    fn a_reused_pair_without_its_sibling_matrix_names_the_path_it_looked_for() {
        let pair = pair();
        let directory = std::env::temp_dir().join(format!(
            "gps-population-snvs-no-sibling-matrix-{}",
            std::process::id()
        ));
        fs::create_dir_all(&directory).expect("directory");
        let metadata_path = directory.join("population_snvs.metadata.json");
        fs::write(&metadata_path, serde_json::to_vec(&pair.metadata).unwrap()).expect("metadata");
        let expected_matrix = directory.join("population_snvs.tsv");
        assert!(!expected_matrix.exists());
        let error = validate_reusable_population_snvs(
            &metadata_path,
            &pair.reference,
            FASTA_SHA,
            TOLERANCE,
            false,
        )
        .expect_err("missing matrix");
        fs::remove_dir_all(&directory).expect("cleanup");
        let PopulationSnvsPreparationError::AuthoritativeFasta { detail } = error else {
            panic!("expected AuthoritativeFasta");
        };
        let prefix = format!(
            "no se puede leer la matriz {}, que debe estar junto a los metadatos: ",
            expected_matrix.display()
        );
        assert!(detail.starts_with(&prefix), "unexpected detail: {detail}");
    }

    #[test]
    fn a_reused_pair_whose_metadata_file_cannot_be_read_is_rejected() {
        let pair = pair();
        let path = std::env::temp_dir().join(format!(
            "gps-population-snvs-no-such-metadata-{}.json",
            std::process::id()
        ));
        assert!(!path.exists(), "the fixture path must not exist");
        let error =
            validate_reusable_population_snvs(&path, &pair.reference, FASTA_SHA, TOLERANCE, false)
                .expect_err("unreadable metadata");
        let PopulationSnvsPreparationError::AuthoritativeFasta { detail } = error else {
            panic!("expected AuthoritativeFasta");
        };
        assert!(
            detail.starts_with("no se puede leer metadata_path: "),
            "unexpected detail: {detail}"
        );
    }

    // ----------------------------------------------------------------------
    // Preparación desde una referencia en memoria
    // ----------------------------------------------------------------------

    fn supplied(rows: &[&str]) -> PopulationSnvsTsv {
        let mut text = String::from("contig\tposition\tref\talt\tpopulation_id\talt_af\n");
        for row in rows {
            text.push_str(row);
            text.push('\n');
        }
        parse_supplied_tsv(text.as_bytes()).expect("supplied shape")
    }

    fn prepare(rows: &[&str]) -> Result<PreparedPopulationSnvs, PopulationSnvsPreparationError> {
        prepare_population_snvs_from_reference(
            &supplied(rows),
            &build_reference(true),
            TOLERANCE,
            false,
        )
    }

    #[test]
    fn a_row_is_rejected_for_the_exact_field_that_breaks() {
        let invalid_allele = |row_index: usize, field: &'static str, value: &str| {
            PopulationSnvsPreparationError::InvalidAllele {
                row_index,
                field,
                value: value.to_owned(),
            }
        };
        let cases: Vec<(&str, Vec<&str>, PopulationSnvsPreparationError)> = vec![
            (
                "population of the second row",
                vec!["chr1\t2\tC\tT\tP1\t0.5", "chr1\t2\tC\tT\tbad id\t0.5"],
                PopulationSnvsPreparationError::InvalidPopulation {
                    row_index: 1,
                    value: "bad id".to_owned(),
                },
            ),
            (
                "contig syntax",
                vec!["bad contig\t2\tC\tT\tP1\t0.5"],
                PopulationSnvsPreparationError::InvalidContig {
                    row_index: 0,
                    value: "bad contig".to_owned(),
                },
            ),
            (
                "contig not in the reference",
                vec!["chr9\t2\tC\tT\tP1\t0.5"],
                PopulationSnvsPreparationError::InvalidContig {
                    row_index: 0,
                    value: "chr9".to_owned(),
                },
            ),
            (
                "position not a number",
                vec!["chr1\tx\tC\tT\tP1\t0.5"],
                PopulationSnvsPreparationError::InvalidPosition {
                    row_index: 0,
                    value: "x".to_owned(),
                },
            ),
            (
                "position zero",
                vec!["chr1\t0\tC\tT\tP1\t0.5"],
                PopulationSnvsPreparationError::InvalidPosition {
                    row_index: 0,
                    value: "0".to_owned(),
                },
            ),
            (
                "position with a leading zero",
                vec!["chr1\t08\tC\tT\tP1\t0.5"],
                PopulationSnvsPreparationError::InvalidPosition {
                    row_index: 0,
                    value: "08".to_owned(),
                },
            ),
            (
                "position with a plus sign",
                vec!["chr1\t+8\tC\tT\tP1\t0.5"],
                PopulationSnvsPreparationError::InvalidPosition {
                    row_index: 0,
                    value: "+8".to_owned(),
                },
            ),
            (
                "position one past the contig",
                vec!["chr1\t11\tA\tT\tP1\t0.5"],
                PopulationSnvsPreparationError::ReferenceMismatch {
                    row_index: 0,
                    contig: "chr1".to_owned(),
                    position: 11,
                    expected: b'?',
                    found: b'A',
                },
            ),
            (
                "ref differs from the reference base",
                vec!["chr1\t2\tA\tT\tP1\t0.5"],
                PopulationSnvsPreparationError::ReferenceMismatch {
                    row_index: 0,
                    contig: "chr1".to_owned(),
                    position: 2,
                    expected: b'C',
                    found: b'A',
                },
            ),
            (
                "same ALT and population twice",
                vec!["chr1\t2\tC\tT\tP1\t0.5", "chr1\t2\tC\tT\tP1\t0.25"],
                PopulationSnvsPreparationError::DuplicateAltPopulation {
                    contig: "chr1".to_owned(),
                    position: 2,
                    reference: b'C',
                    alternate: b'T',
                    population: "P1".to_owned(),
                },
            ),
            (
                "an ALT that one population lacks",
                vec!["chr1\t2\tC\tT\tP1\t0.5", "chr1\t3\tG\tT\tP2\t0.5"],
                PopulationSnvsPreparationError::IncompletePopulationCoverage {
                    contig: "chr1".to_owned(),
                    position: 2,
                    reference: b'C',
                    alternate: b'T',
                },
            ),
            (
                "alt equal to ref",
                vec!["chr1\t2\tC\tC\tP1\t0.5"],
                invalid_allele(0, "alt", "C"),
            ),
            (
                "ref is N",
                vec!["chr1\t2\tN\tT\tP1\t0.5"],
                invalid_allele(0, "ref", "N"),
            ),
            (
                "ref has two bases",
                vec!["chr1\t2\tCC\tT\tP1\t0.5"],
                invalid_allele(0, "ref", "CC"),
            ),
            (
                "ref is lower case",
                vec!["chr1\t2\tc\tT\tP1\t0.5"],
                invalid_allele(0, "ref", "c"),
            ),
            (
                "alt is N",
                vec!["chr1\t2\tC\tN\tP1\t0.5"],
                invalid_allele(0, "alt", "N"),
            ),
            (
                "alt has two bases",
                vec!["chr1\t2\tC\tTT\tP1\t0.5"],
                invalid_allele(0, "alt", "TT"),
            ),
            (
                "alt is empty",
                vec!["chr1\t2\tC\t\tP1\t0.5"],
                invalid_allele(0, "alt", ""),
            ),
        ];
        for (name, rows, expected) in cases {
            assert_eq!(prepare(&rows).map(|_| ()), Err(expected), "{name}");
        }
        // El último nucleótido del contig sí se acepta: el límite exacto de
        // «posición fuera del contig» es longitud + 1.
        prepare(&["chr1\t10\tC\tT\tP1\t0.5"]).expect("the last base is inside the contig");
    }

    #[test]
    fn preparing_from_an_unmaterialized_reference_is_its_own_error() {
        let table = supplied(&["chr1\t2\tC\tT\tP1\t0.5"]);
        let result = prepare_population_snvs_from_reference(
            &table,
            &build_reference(false),
            TOLERANCE,
            false,
        );
        assert_eq!(
            result.map(|_| ()),
            Err(PopulationSnvsPreparationError::BasesNotMaterialized)
        );
    }

    #[test]
    fn an_empty_active_matrix_needs_authorization_each_time_it_is_asked() {
        let table = supplied(&["chr1\t2\tC\tT\tP1\t0", "chr1\t2\tC\tT\tP2\t0"]);
        let reference = build_reference(true);
        assert_eq!(
            prepare_population_snvs_from_reference(&table, &reference, TOLERANCE, false)
                .map(|_| ()),
            Err(PopulationSnvsPreparationError::EmptyActiveMatrixRequiresAuthorization)
        );
        let prepared = prepare_population_snvs_from_reference(&table, &reference, TOLERANCE, true)
            .expect("authorized empty matrix");
        assert!(prepared.rows().is_empty());
        assert!(prepared.empty_authorized());
        assert_eq!(
            prepared.clone().authorize_empty_matrix(false).map(|_| ()),
            Err(PopulationSnvsPreparationError::EmptyActiveMatrixRequiresAuthorization)
        );
        assert!(
            prepared
                .authorize_empty_matrix(true)
                .expect("authorized")
                .empty_authorized()
        );
    }

    #[test]
    fn alt_frequencies_are_rejected_for_the_exact_value_that_breaks() {
        let table = supplied(&["chr1\t2\tC\tT\tP1\t0.5", "chr1\t3\tG\tT\tP1\tabc"]);
        assert_eq!(
            normalize_alt_frequencies(&table, TOLERANCE).map(|_| ()),
            Err(AltFrequencyNormalizationError::InvalidFrequency {
                row_index: 1,
                literal: "abc".to_owned(),
            })
        );
        let table = supplied(&["chr1\t2\tC\tT\tP1\t"]);
        assert_eq!(
            normalize_alt_frequencies(&table, TOLERANCE).map(|_| ()),
            Err(AltFrequencyNormalizationError::InvalidFrequency {
                row_index: 0,
                literal: String::new(),
            })
        );
        // Cada tolerancia fuera de [0, 1e-6]. `NaN` y el infinito los rechazaría
        // también el intervalo, de modo que la mitad `!is_finite()` es redundante
        // (ver el comentario de la guarda): aquí se comprueba el resultado.
        let table = supplied(&["chr1\t2\tC\tT\tP1\t0.5"]);
        for tolerance in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e-9, 2e-6, 1.0] {
            let result = normalize_alt_frequencies(&table, tolerance);
            match result {
                Err(AltFrequencyNormalizationError::InvalidTolerance { value }) => {
                    assert_eq!(value.to_bits(), tolerance.to_bits(), "{tolerance}");
                }
                other => panic!("tolerance {tolerance}: {:?}", other.map(|_| ())),
            }
        }
        for tolerance in [0.0, 1e-6] {
            normalize_alt_frequencies(&table, tolerance).expect("boundary tolerance");
        }
    }

    #[test]
    fn alt_frequencies_outside_zero_one_or_summing_too_high_are_rejected_with_their_numbers() {
        for literal in ["1.5", "-0.1", "NaN", "inf", "-inf"] {
            let table = supplied(&[&format!("chr1\t2\tC\tT\tP1\t{literal}")]);
            assert_eq!(
                normalize_alt_frequencies(&table, TOLERANCE).map(|_| ()),
                Err(AltFrequencyNormalizationError::InvalidFrequency {
                    row_index: 0,
                    literal: literal.to_owned(),
                }),
                "{literal}"
            );
        }
        // Dos ALT del mismo locus que suman 1,2: exceso, con la suma, el delta y
        // la tolerancia exactos.
        let table = supplied(&["chr1\t2\tC\tA\tP1\t0.6", "chr1\t2\tC\tT\tP1\t0.6"]);
        assert_eq!(
            normalize_alt_frequencies(&table, TOLERANCE).map(|_| ()),
            Err(AltFrequencyNormalizationError::ExcessAltSum {
                population_id: "P1".to_owned(),
                contig: "chr1".to_owned(),
                position: "2".to_owned(),
                reference: "C".to_owned(),
                sum: 1.2,
                delta: (1.2_f64 - 1.0).abs(),
                tolerance: TOLERANCE,
            })
        );
        // El límite es inclusivo: con la tolerancia igual al delta exacto, se normaliza.
        let upper = 0.25 + f64::EPSILON;
        let table = supplied(&[
            "chr1\t2\tC\tA\tP1\t0.75",
            &format!("chr1\t2\tC\tT\tP1\t{upper}"),
        ]);
        let normalized = normalize_alt_frequencies(&table, f64::EPSILON).expect("inclusive limit");
        assert_eq!(
            normalized.groups()[0].decision(),
            AltFrequencyDecision::Normalized
        );
        assert_eq!(
            normalized.groups()[0].original_delta().to_bits(),
            f64::EPSILON.to_bits()
        );
        // Un ulp menos de tolerancia ya no lo admite.
        assert!(matches!(
            normalize_alt_frequencies(&table, f64::EPSILON / 2.0),
            Err(AltFrequencyNormalizationError::ExcessAltSum { .. })
        ));
    }

    // ----------------------------------------------------------------------
    // Formatos TSV y metadatos
    // ----------------------------------------------------------------------

    #[test]
    fn the_canonical_tsv_parser_rejects_each_malformed_shape_for_its_own_reason() {
        let header = TSV_HEADER.join("\t");
        let row = "VAR_X\tchr1\t1\tA\tC\tP1\t0.5\t";
        let cases: Vec<(&str, Vec<u8>, PopulationSnvsFormatError)> = vec![
            (
                "not UTF-8",
                vec![0xff, 0xfe],
                PopulationSnvsFormatError::InvalidUtf8,
            ),
            (
                "VCF file format line",
                b"##fileformat=VCFv4.3\n".to_vec(),
                PopulationSnvsFormatError::VcfInput,
            ),
            (
                "VCF column header",
                b"#CHROM\tPOS\n".to_vec(),
                PopulationSnvsFormatError::VcfInput,
            ),
            (
                "empty input",
                Vec::new(),
                PopulationSnvsFormatError::MissingHeader,
            ),
            (
                "wrong header",
                b"a\tb\n".to_vec(),
                PopulationSnvsFormatError::IncorrectHeader {
                    found: "a\tb".to_owned(),
                },
            ),
            (
                "empty line after the header",
                format!("{header}\n\n").into_bytes(),
                PopulationSnvsFormatError::IncorrectColumnCount { line: 2, found: 0 },
            ),
            (
                "empty line in the middle",
                format!("{header}\n{row}\n\n").into_bytes(),
                PopulationSnvsFormatError::IncorrectColumnCount { line: 3, found: 0 },
            ),
            (
                "one column too many",
                format!("{header}\n{row}\textra\n").into_bytes(),
                PopulationSnvsFormatError::IncorrectColumnCount { line: 2, found: 9 },
            ),
            (
                "one column too few",
                format!("{header}\n{row}\n{}\n", "VAR_X\tchr1").into_bytes(),
                PopulationSnvsFormatError::IncorrectColumnCount { line: 3, found: 2 },
            ),
        ];
        for (name, bytes, expected) in cases {
            assert_eq!(
                parse_canonical_tsv(&bytes).map(|_| ()),
                Err(expected),
                "{name}"
            );
        }
    }

    #[test]
    fn the_supplied_tsv_parser_rejects_each_malformed_shape_for_its_own_reason() {
        let six = "contig\tposition\tref\talt\tpopulation_id\talt_af";
        let seven = format!("{six}\tsource_id");
        let cases: Vec<(&str, Vec<u8>, PopulationSnvsFormatError)> = vec![
            (
                "not UTF-8",
                vec![0xff],
                PopulationSnvsFormatError::InvalidUtf8,
            ),
            (
                "VCF file format line",
                b"##fileformat=VCFv4.3\n".to_vec(),
                PopulationSnvsFormatError::VcfInput,
            ),
            (
                "VCF column header",
                b"#CHROM\tPOS\n".to_vec(),
                PopulationSnvsFormatError::VcfInput,
            ),
            (
                "empty input",
                Vec::new(),
                PopulationSnvsFormatError::MissingHeader,
            ),
            (
                "wrong header",
                b"contig\tposition\n".to_vec(),
                PopulationSnvsFormatError::IncorrectHeader {
                    found: "contig\tposition".to_owned(),
                },
            ),
            (
                "six-column header with a seven-column row",
                format!("{six}\nchr1\t1\tA\tC\tP1\t0.5\tsrc\n").into_bytes(),
                PopulationSnvsFormatError::IncorrectColumnCount { line: 2, found: 7 },
            ),
            (
                "seven-column header with a six-column row",
                format!("{seven}\nchr1\t1\tA\tC\tP1\t0.5\n").into_bytes(),
                PopulationSnvsFormatError::IncorrectColumnCount { line: 2, found: 6 },
            ),
            (
                "six-column header with a short row",
                format!("{six}\nchr1\t1\n").into_bytes(),
                PopulationSnvsFormatError::IncorrectColumnCount { line: 2, found: 2 },
            ),
        ];
        for (name, bytes, expected) in cases {
            assert_eq!(
                parse_supplied_tsv(&bytes).map(|_| ()),
                Err(expected),
                "{name}"
            );
        }
        parse_supplied_tsv(format!("{seven}\nchr1\t1\tA\tC\tP1\t0.5\tsrc\n").as_bytes())
            .expect("seven columns with a source");
    }

    fn missing_field(field: &str) -> PopulationSnvsFormatError {
        PopulationSnvsFormatError::MissingMetadataField {
            field: field.to_owned(),
        }
    }

    fn invalid_field(field: &str, expected: &'static str) -> PopulationSnvsFormatError {
        PopulationSnvsFormatError::InvalidMetadataField {
            field: field.to_owned(),
            expected,
        }
    }

    fn unexpected_field(field: &str) -> PopulationSnvsFormatError {
        PopulationSnvsFormatError::UnexpectedMetadataField {
            field: field.to_owned(),
        }
    }

    fn invalid_literal(field: &str, expected: &'static str) -> PopulationSnvsFormatError {
        PopulationSnvsFormatError::InvalidMetadataLiteral {
            field: field.to_owned(),
            expected,
        }
    }

    fn remove_at(value: &mut Value, pointer: &str, key: &str) {
        value
            .pointer_mut(pointer)
            .expect("pointer")
            .as_object_mut()
            .expect("object")
            .remove(key)
            .expect("key");
    }

    /// Las cuatro sumas y la tolerancia se comparan por bits, no como números:
    /// `-0.0` es igual a `0.0` pero no lo es el patrón de bits, y el juez exige
    /// los bits que escribió el productor.
    #[test]
    fn the_audit_judge_compares_sums_and_tolerance_bit_by_bit() {
        let pair = pair();
        // P1 chr2:1 está descartado entero (sumas 0), P2 chr2:2 suma exactamente 1
        // (deltas 0).
        for (field, population, contig, position) in [
            ("original_sum", "P1", "chr2", "1"),
            ("effective_sum", "P1", "chr2", "1"),
            ("original_delta", "P2", "chr2", "2"),
            ("effective_delta", "P2", "chr2", "2"),
        ] {
            let mut audit = pair.metadata["normalization_audit"].clone();
            let target = group(&mut audit, population, contig, position);
            assert_eq!(
                target[field].as_f64().unwrap().to_bits(),
                0,
                "{field}: premise"
            );
            target[field] = json!(-0.0_f64);
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason("sumas no coinciden")),
                "{field}"
            );
        }
        // Una tolerancia de cero positivo frente a un `-0.0` publicado.
        let mut audit = pair.metadata["normalization_audit"].clone();
        for group_value in audit.as_array_mut().unwrap() {
            group_value["tolerance"] = json!(-0.0_f64);
        }
        assert_eq!(
            check_audit(&pair, &audit, 0.0),
            Err(audit_reason("identidad o algoritmo de grupo inválido"))
        );
    }

    /// Los cuatro campos de texto de cabecera de un grupo se leen todos antes de
    /// parsear ninguno: con dos fallos a la vez, el motivo es el del **primero** de
    /// la lista. Distingue cada lectura de una que cayera a «vacío» y fallara más
    /// tarde por otro motivo.
    #[test]
    fn the_audit_judge_names_the_first_group_text_field_that_is_not_a_string() {
        let pair = pair();
        let cases: Vec<(&str, Edit, &str)> = vec![
            (
                "population_id before contig",
                Box::new(|a| {
                    a[0]["population_id"] = json!(5);
                    a[0]["contig"] = json!(5);
                }),
                "population_id",
            ),
            (
                "contig before position",
                Box::new(|a| {
                    a[0]["contig"] = json!(5);
                    a[0]["position"] = json!(5);
                }),
                "contig",
            ),
            (
                "position before ref",
                Box::new(|a| {
                    a[0]["position"] = json!(5);
                    a[0]["ref"] = json!(5);
                }),
                "position",
            ),
            (
                "ref before the population id syntax",
                Box::new(|a| {
                    a[0]["ref"] = json!(5);
                    a[0]["population_id"] = json!("bad id");
                }),
                "ref",
            ),
        ];
        for (name, edit, reason) in cases {
            let mut audit = pair.metadata["normalization_audit"].clone();
            edit(&mut audit);
            assert_eq!(
                check_audit(&pair, &audit, TOLERANCE),
                Err(audit_reason(reason)),
                "{name}"
            );
        }
    }

    /// `1e999` es un número JSON legible pero no cabe en un `f64`: `as_f64()` da
    /// `None`, no un infinito. Por eso las mitades `is_finite()` de los filtros de
    /// la auditoría no actúan (están anotadas como equivalentes) y el rechazo lo
    /// produce el propio `as_f64`.
    #[test]
    fn an_out_of_range_json_number_is_not_a_finite_float_for_the_audit_judge() {
        let huge: Value = serde_json::from_str("1e999").expect("a legal JSON number");
        assert_eq!(huge.as_f64(), None);
        let pair = pair();
        let mut audit = pair.metadata["normalization_audit"].clone();
        audit[0]["tolerance"] = huge.clone();
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("tolerance"))
        );
        let mut audit = pair.metadata["normalization_audit"].clone();
        audit[0]["values"][0]["original"] = huge.clone();
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("original"))
        );
        let mut audit = pair.metadata["normalization_audit"].clone();
        audit[0]["values"][0]["effective"] = huge.clone();
        assert_eq!(
            check_audit(&pair, &audit, TOLERANCE),
            Err(audit_reason("effective"))
        );
        let mut audit = pair.metadata["normalization_audit"].clone();
        audit[0]["reference_frequency"] = huge.clone();
        assert_eq!(
            restore_published_allele_frequencies(&mut pair.prepared.clone(), &audit),
            Err(audit_reason("reference_frequency"))
        );
        let mut audit = pair.metadata["normalization_audit"].clone();
        audit[0]["values"][0]["effective"] = huge;
        assert_eq!(
            restore_published_allele_frequencies(&mut pair.prepared.clone(), &audit),
            Err(audit_reason("effective"))
        );
    }

    /// Intento ejecutado sobre las guardas **inalcanzables** de
    /// `validate_reusable_population_snvs_bytes` y de la restauración (las que
    /// repiten una forma que `parse_metadata_json` ya exige). No discrimina ninguna
    /// de ellas —por definición no se disparan—: demuestra la premisa de su
    /// anotación, que es que cada metadato roto de esa forma lo rechaza antes el
    /// analizador, con el texto del analizador y no con el de la guarda.
    #[test]
    fn the_shape_guards_of_the_reuse_judge_are_shadowed_by_the_metadata_parser() {
        let pair = pair();
        let provenance = "\"supplied\" o \"generated\"";
        let cases: Vec<(&str, Edit, PopulationSnvsFormatError)> = vec![
            (
                "reference not an object",
                Box::new(|m| m["reference"] = json!(1)),
                invalid_field("reference", "un objeto"),
            ),
            (
                "reference_fasta not an object",
                Box::new(|m| m["reference_fasta"] = json!(1)),
                invalid_field("reference_fasta", "un objeto"),
            ),
            (
                "matrix not an object",
                Box::new(|m| m["matrix"] = json!(1)),
                invalid_field("matrix", "un objeto"),
            ),
            (
                "matrix.sha256 not a string",
                Box::new(|m| m["matrix"]["sha256"] = json!(1)),
                invalid_field("matrix.sha256", "una cadena"),
            ),
            (
                "normalization_audit absent",
                Box::new(|m| {
                    m.as_object_mut().unwrap().remove("normalization_audit");
                }),
                invalid_field("normalization_audit", "una lista"),
            ),
            (
                "normalization_audit not a list",
                Box::new(|m| m["normalization_audit"] = json!({})),
                invalid_field("normalization_audit", "una lista"),
            ),
            (
                "populations not a list",
                Box::new(|m| m["populations"] = json!({})),
                invalid_field("populations", "una lista"),
            ),
            (
                "population_id not a string",
                Box::new(|m| m["populations"][0]["population_id"] = json!(1)),
                invalid_field("populations[].population_id", "una cadena"),
            ),
            (
                "provenance not a string",
                Box::new(|m| m["provenance"] = json!(1)),
                invalid_field("provenance", provenance),
            ),
            (
                "provenance outside the two literals",
                Box::new(|m| m["provenance"] = json!("otra")),
                invalid_literal("provenance", provenance),
            ),
            (
                "discarded pairs not a list",
                Box::new(|m| m["discarded_variant_pairs"] = json!({})),
                invalid_field("discarded_variant_pairs", "una lista"),
            ),
            (
                "discarded pair not an object",
                Box::new(|m| m["discarded_variant_pairs"][0] = json!(1)),
                invalid_field("discarded_variant_pairs[]", "un objeto"),
            ),
            (
                "discarded pair population_id not a string",
                Box::new(|m| m["discarded_variant_pairs"][0]["population_id"] = json!(1)),
                invalid_field("discarded_variant_pairs[].population_id", "un PopulationId"),
            ),
            (
                "discarded pair variant_id not a string",
                Box::new(|m| m["discarded_variant_pairs"][0]["variant_id"] = json!(1)),
                invalid_field("discarded_variant_pairs[].variant_id", "un VariantId"),
            ),
        ];
        for (name, edit, parser_error) in cases {
            let mut metadata = pair.metadata.clone();
            edit(&mut metadata);
            assert_eq!(
                validate_metadata(&pair, &metadata),
                Err(fasta_detail(&parser_error.to_string())),
                "{name}"
            );
        }
        assert_eq!(
            validate_bytes(&pair, b"[]", &pair.tsv),
            Err(fasta_detail(
                &PopulationSnvsFormatError::MetadataMustBeObject.to_string()
            ))
        );
    }

    /// El lado «aceptado» de la lista cerrada de campos de cada población
    /// (especificaciones: `description`, `sample_size`, `allele_count`, `source`,
    /// `filters` y `data_version`, además de `population_id`).
    #[test]
    fn a_population_with_every_admitted_field_is_accepted() {
        let pair = pair();
        let mut metadata = pair.metadata.clone();
        for population in metadata["populations"].as_array_mut().unwrap() {
            population["description"] = json!("d");
            population["sample_size"] = json!(10);
            population["allele_count"] = json!(20);
            population["source"] = json!("s");
            population["filters"] = json!(["PASS"]);
            population["data_version"] = json!("v1");
        }
        parse_metadata_json(&serde_json::to_vec(&metadata).unwrap()).expect("parser accepts");
        validate_metadata(&pair, &metadata).expect("the reuse judge accepts");
        for field in [
            "description",
            "sample_size",
            "allele_count",
            "source",
            "filters",
            "data_version",
        ] {
            let mut single = pair.metadata.clone();
            single["populations"][0][field] = metadata["populations"][0][field].clone();
            parse_metadata_json(&serde_json::to_vec(&single).unwrap())
                .unwrap_or_else(|error| panic!("{field}: {error}"));
        }
    }

    #[test]
    fn the_metadata_parser_rejects_each_malformed_point_for_its_own_reason() {
        let pair = pair();
        let provenance = "\"supplied\" o \"generated\"";
        let mut cases: Vec<(String, Edit, PopulationSnvsFormatError)> = Vec::new();
        let mut add = |name: &str, edit: Edit, expected: PopulationSnvsFormatError| {
            cases.push((name.to_owned(), edit, expected));
        };
        // `required_object`: tres llamantes.
        for key in ["reference", "reference_fasta", "matrix"] {
            let key_owned = key.to_owned();
            add(
                &format!("{key} missing"),
                Box::new(move |m| {
                    m.as_object_mut().unwrap().remove(&key_owned);
                }),
                missing_field(key),
            );
            let key_owned = key.to_owned();
            add(
                &format!("{key} not an object"),
                Box::new(move |m| m[key_owned.as_str()] = json!([])),
                invalid_field(key, "un objeto"),
            );
        }
        // `required_string`: cinco llamantes.
        for (pointer, leaf, full) in [
            ("/reference", "reference_id", "reference.reference_id"),
            (
                "/reference",
                "canonical_fingerprint",
                "reference.canonical_fingerprint",
            ),
            ("/reference_fasta", "sha256", "reference_fasta.sha256"),
            ("/matrix", "sha256", "matrix.sha256"),
            (
                "/populations/0",
                "population_id",
                "populations[].population_id",
            ),
        ] {
            add(
                &format!("{full} missing"),
                Box::new(move |m| remove_at(m, pointer, leaf)),
                missing_field(full),
            );
            add(
                &format!("{full} not a string"),
                Box::new(move |m| *m.pointer_mut(&format!("{pointer}/{leaf}")).unwrap() = json!(1)),
                invalid_field(full, "una cadena"),
            );
        }
        // `required_literal`: cuatro llamantes.
        for (pointer, key, expected) in [
            ("", "document_kind", "population_snvs_metadata"),
            ("", "schema_version", "1"),
            ("", "coordinate_system", "base1_closed"),
            ("/matrix", "path", "population_snvs.tsv"),
        ] {
            add(
                &format!("{key} missing"),
                Box::new(move |m| remove_at(m, pointer, key)),
                missing_field(key),
            );
            add(
                &format!("{key} wrong literal"),
                Box::new(move |m| {
                    *m.pointer_mut(&format!("{pointer}/{key}")).unwrap() = json!("otro")
                }),
                invalid_literal(key, expected),
            );
        }
        // `closed_object`: seis llamantes.
        for (pointer, prefix) in [
            ("", ""),
            ("/reference", "reference."),
            ("/reference_fasta", "reference_fasta."),
            ("/matrix", "matrix."),
            ("/populations/0", "populations[]."),
            ("/discarded_variant_pairs/0", "discarded_variant_pairs[]."),
        ] {
            add(
                &format!("extra field under '{prefix}'"),
                Box::new(move |m| {
                    let target = if pointer.is_empty() {
                        &mut *m
                    } else {
                        m.pointer_mut(pointer).unwrap()
                    };
                    target["extra"] = json!(1);
                }),
                unexpected_field(&format!("{prefix}extra")),
            );
        }
        // Poblaciones.
        add(
            "populations missing",
            Box::new(|m| {
                m.as_object_mut().unwrap().remove("populations");
            }),
            missing_field("populations"),
        );
        add(
            "populations not a list",
            Box::new(|m| m["populations"] = json!({})),
            invalid_field("populations", "una lista"),
        );
        add(
            "population not an object",
            Box::new(|m| m["populations"][0] = json!("P1")),
            invalid_field("populations[]", "un objeto"),
        );
        // Procedencia.
        add(
            "provenance missing",
            Box::new(|m| {
                m.as_object_mut().unwrap().remove("provenance");
            }),
            missing_field("provenance"),
        );
        add(
            "provenance not a string",
            Box::new(|m| m["provenance"] = json!(1)),
            invalid_field("provenance", provenance),
        );
        add(
            "provenance unknown literal",
            Box::new(|m| m["provenance"] = json!("otra")),
            invalid_literal("provenance", provenance),
        );
        // Pares descartados.
        add(
            "pairs missing",
            Box::new(|m| {
                m.as_object_mut().unwrap().remove("discarded_variant_pairs");
            }),
            missing_field("discarded_variant_pairs"),
        );
        add(
            "pairs not a list",
            Box::new(|m| m["discarded_variant_pairs"] = json!({})),
            invalid_field("discarded_variant_pairs", "una lista"),
        );
        add(
            "pair not an object",
            Box::new(|m| m["discarded_variant_pairs"][0] = json!(1)),
            invalid_field("discarded_variant_pairs[]", "un objeto"),
        );
        add(
            "pair population_id missing",
            Box::new(|m| remove_at(m, "/discarded_variant_pairs/0", "population_id")),
            invalid_field("discarded_variant_pairs[].population_id", "un PopulationId"),
        );
        add(
            "pair population_id not a valid id",
            Box::new(|m| m["discarded_variant_pairs"][0]["population_id"] = json!("bad id")),
            invalid_field("discarded_variant_pairs[].population_id", "un PopulationId"),
        );
        add(
            "pair variant_id missing",
            Box::new(|m| remove_at(m, "/discarded_variant_pairs/0", "variant_id")),
            invalid_field("discarded_variant_pairs[].variant_id", "un VariantId"),
        );
        add(
            "pair variant_id not a valid id",
            Box::new(|m| m["discarded_variant_pairs"][0]["variant_id"] = json!("bad id")),
            invalid_field("discarded_variant_pairs[].variant_id", "un VariantId"),
        );
        let ordering = invalid_field(
            "discarded_variant_pairs",
            "pares población×VariantId ordenados y sin duplicados",
        );
        add(
            "pairs descending",
            Box::new(|m| {
                m["discarded_variant_pairs"]
                    .as_array_mut()
                    .unwrap()
                    .swap(0, 1)
            }),
            ordering.clone(),
        );
        add(
            "pairs repeated",
            Box::new(|m| {
                let copy = m["discarded_variant_pairs"][0].clone();
                m["discarded_variant_pairs"]
                    .as_array_mut()
                    .unwrap()
                    .insert(1, copy);
            }),
            ordering,
        );
        // Auditoría.
        add(
            "audit not a list",
            Box::new(|m| m["normalization_audit"] = json!({})),
            invalid_field("normalization_audit", "una lista"),
        );
        add(
            "audit missing",
            Box::new(|m| {
                m.as_object_mut().unwrap().remove("normalization_audit");
            }),
            invalid_field("normalization_audit", "una lista"),
        );

        for (name, edit, expected) in cases {
            let mut metadata = pair.metadata.clone();
            edit(&mut metadata);
            assert_ne!(metadata, pair.metadata, "{name}: no-op edit");
            let bytes = serde_json::to_vec(&metadata).unwrap();
            assert_eq!(
                parse_metadata_json(&bytes).map(|_| ()),
                Err(expected),
                "{name}"
            );
        }
        parse_metadata_json(&serde_json::to_vec(&pair.metadata).unwrap()).expect("valid metadata");
        assert_eq!(
            parse_metadata_json(b"{").map(|_| ()),
            Err(PopulationSnvsFormatError::InvalidMetadataJson)
        );
        assert_eq!(
            parse_metadata_json(b"[]").map(|_| ()),
            Err(PopulationSnvsFormatError::MetadataMustBeObject)
        );
    }
}
