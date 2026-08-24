//! Contrato externo v1 de `truth.jsonl`, independiente del dominio Rust.
//!
//! Solo comprueba estructura en memoria. Proyección, E/S y reconstrucción son
//! responsabilidades posteriores (F2.4.2 y F2.4.3).

use std::collections::{BTreeMap, BTreeSet};

pub const TRUTH_SCHEMA_VERSION: &str = "1";
pub const TRUTH_DOCUMENT_KIND: &str = "gps_truth";
pub const TRUTH_COORDINATE_CONVENTION: &str = "zero_based_half_open_v1";

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TruthHeader {
    pub document_kind: String,
    pub schema_version: String,
    pub execution_id: String,
    pub reference_id: String,
    pub coordinate_convention: String,
    pub copy_count: u64,
    pub genotype_count: u64,
    pub source_digests: Vec<NamedDigest>,
    pub group_digests: Vec<NamedDigest>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NamedDigest {
    pub source: TruthSourceRef,
    pub sha256: String,
}
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct TruthSourceRef {
    pub kind: TruthSourceKind,
    pub source_id: String,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum TruthSourceKind {
    Reference,
    Founder,
    Explicit,
    Random,
    Matrix,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TruthDocument {
    pub header: TruthHeader,
    pub records: Vec<TruthRecord>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TruthRecord {
    HaplotypeCopy(HaplotypeCopyRecord),
    Genotype(GenotypeRecord),
}
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct CopyIdentity {
    pub reference_id: String,
    pub individual_id: String,
    pub sequence_id: String,
    pub copy_index: u8,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HaplotypeCopyRecord {
    pub identity: CopyIdentity,
    pub phase: u8,
    pub provenance: CopyProvenance,
    pub segments: Vec<TruthSegment>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CopyProvenance {
    pub origin: String,
    pub sources: Vec<TruthSourceRef>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TruthSegment {
    pub source: SegmentSource,
    pub modifications: Vec<TruthModification>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SegmentSource {
    pub source: TruthSourceRef,
    pub reference_id: String,
    pub sequence_id: String,
    pub start: u64,
    pub end: u64,
    pub sequence_length: u64,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TruthModification {
    Snp {
        event_id: String,
        position: u64,
        reference: String,
        alternate: String,
        provenance: ModificationProvenance,
    },
    Insertion {
        event_id: String,
        point: u64,
        bases: String,
        provenance: ModificationProvenance,
    },
    Deletion {
        event_id: String,
        start: u64,
        end: u64,
        provenance: ModificationProvenance,
    },
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModificationProvenance {
    pub kind: ModificationProvenanceKind,
    pub version: String,
    pub source: TruthSourceRef,
    pub matrix: Option<MatrixProvenance>,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ModificationProvenanceKind {
    Explicit,
    Random,
    Matrix,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MatrixProvenance {
    pub matrix_source_id: String,
    pub population_id: String,
    pub variant_id: String,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GenotypeRecord {
    pub case_id: String,
    pub reference_id: String,
    pub individual_id: String,
    pub locus: TruthLocus,
    pub alleles: Vec<String>,
    pub allele_dosages: Vec<u8>,
    pub carriers: Vec<AlleleCarrier>,
}
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct TruthLocus {
    pub sequence_id: String,
    pub start: u64,
    pub end: u64,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlleleCarrier {
    pub allele_index: usize,
    pub copy: CopyIdentity,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TruthContractError {
    UnsupportedDocumentKind,
    UnsupportedVersion,
    UnsupportedCoordinateConvention,
    EmptySourceDigests,
    DuplicateDigest {
        source: TruthSourceRef,
    },
    EmptySegment {
        identity: CopyIdentity,
    },
    InvalidInsertion {
        event_id: String,
    },
    InvalidDeletion {
        event_id: String,
    },
    MissingDigest {
        source: TruthSourceRef,
    },
    InvalidCopyPhase {
        identity: CopyIdentity,
        phase: u8,
    },
    CopyReferenceMismatch {
        identity: CopyIdentity,
    },
    DuplicateCopy {
        identity: CopyIdentity,
    },
    RecordsOutOfOrder,
    InvalidSegmentInterval {
        identity: CopyIdentity,
    },
    SegmentSourceMismatch {
        identity: CopyIdentity,
    },
    ModificationOutsideSegment {
        event_id: String,
    },
    ModificationsNotOrdered {
        event_id: String,
    },
    OverlappingModifications {
        event_id: String,
    },
    SnpDoesNotChangeBase {
        event_id: String,
    },
    InvalidModificationProvenance {
        event_id: String,
    },
    CopyCountMismatch {
        declared: u64,
        actual: u64,
    },
    GenotypeCountMismatch {
        declared: u64,
        actual: u64,
    },
    EmptyGenotype {
        case_id: String,
    },
    GenotypeReferenceMismatch {
        case_id: String,
    },
    DuplicateGenotype {
        case_id: String,
    },
    UnknownCarrierCopy {
        case_id: String,
        copy: CopyIdentity,
    },
    CarrierOutsideGenotype {
        case_id: String,
        copy: CopyIdentity,
    },
    DuplicateCarrier {
        case_id: String,
        copy: CopyIdentity,
    },
    CarrierAlleleOutOfRange {
        case_id: String,
        allele_index: usize,
    },
    DosageMismatch {
        case_id: String,
        allele_index: usize,
        declared: u8,
        actual: u8,
    },
}

impl TruthDocument {
    pub fn validate_structure(&self) -> Result<(), TruthContractError> {
        let h = &self.header;
        if h.document_kind != TRUTH_DOCUMENT_KIND {
            return Err(TruthContractError::UnsupportedDocumentKind);
        }
        if h.schema_version != TRUTH_SCHEMA_VERSION {
            return Err(TruthContractError::UnsupportedVersion);
        }
        if h.coordinate_convention != TRUTH_COORDINATE_CONVENTION {
            return Err(TruthContractError::UnsupportedCoordinateConvention);
        }
        let digests = source_digest_index(&h.source_digests)?;
        let groups = digest_index(&h.group_digests)?;
        require_digest(
            &digests,
            &TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: h.reference_id.clone(),
            },
        )?;
        let (mut previous, mut copies, mut genotypes) = (None, BTreeSet::new(), BTreeSet::new());
        let (mut copy_count, mut genotype_count) = (0, 0);
        for record in &self.records {
            let key = RecordKey::of(record);
            if previous.as_ref().is_some_and(|prior| prior > &key) {
                return Err(TruthContractError::RecordsOutOfOrder);
            }
            previous = Some(key);
            match record {
                TruthRecord::HaplotypeCopy(copy) => {
                    validate_copy(copy, h, &digests, &groups)?;
                    if !copies.insert(copy.identity.clone()) {
                        return Err(TruthContractError::DuplicateCopy {
                            identity: copy.identity.clone(),
                        });
                    }
                    copy_count += 1;
                }
                TruthRecord::Genotype(genotype) => {
                    validate_genotype(genotype, h, &copies)?;
                    if !genotypes.insert((
                        genotype.reference_id.clone(),
                        genotype.individual_id.clone(),
                        genotype.locus.clone(),
                        genotype.case_id.clone(),
                    )) {
                        return Err(TruthContractError::DuplicateGenotype {
                            case_id: genotype.case_id.clone(),
                        });
                    }
                    genotype_count += 1;
                }
            }
        }
        if h.copy_count != copy_count {
            return Err(TruthContractError::CopyCountMismatch {
                declared: h.copy_count,
                actual: copy_count,
            });
        }
        if h.genotype_count != genotype_count {
            return Err(TruthContractError::GenotypeCountMismatch {
                declared: h.genotype_count,
                actual: genotype_count,
            });
        }
        Ok(())
    }
}

fn source_digest_index(
    digests: &[NamedDigest],
) -> Result<BTreeMap<TruthSourceRef, &NamedDigest>, TruthContractError> {
    if digests.is_empty() {
        return Err(TruthContractError::EmptySourceDigests);
    }
    digest_index(digests)
}
fn digest_index(
    digests: &[NamedDigest],
) -> Result<BTreeMap<TruthSourceRef, &NamedDigest>, TruthContractError> {
    let mut indexed = BTreeMap::new();
    for digest in digests {
        if indexed.insert(digest.source.clone(), digest).is_some() {
            return Err(TruthContractError::DuplicateDigest {
                source: digest.source.clone(),
            });
        }
    }
    Ok(indexed)
}
fn require_digest(
    digests: &BTreeMap<TruthSourceRef, &NamedDigest>,
    source: &TruthSourceRef,
) -> Result<(), TruthContractError> {
    if digests.contains_key(source) {
        Ok(())
    } else {
        Err(TruthContractError::MissingDigest {
            source: source.clone(),
        })
    }
}
fn validate_copy(
    copy: &HaplotypeCopyRecord,
    header: &TruthHeader,
    digests: &BTreeMap<TruthSourceRef, &NamedDigest>,
    groups: &BTreeMap<TruthSourceRef, &NamedDigest>,
) -> Result<(), TruthContractError> {
    if copy.identity.reference_id != header.reference_id {
        return Err(TruthContractError::CopyReferenceMismatch {
            identity: copy.identity.clone(),
        });
    }
    if copy.phase != copy.identity.copy_index {
        return Err(TruthContractError::InvalidCopyPhase {
            identity: copy.identity.clone(),
            phase: copy.phase,
        });
    }
    for source in &copy.provenance.sources {
        require_digest(digests, source)?;
    }
    if copy.provenance.sources.is_empty() || copy.segments.is_empty() {
        return Err(TruthContractError::EmptySegment {
            identity: copy.identity.clone(),
        });
    }
    for segment in &copy.segments {
        let mut previous: Option<ModificationExtent> = None;
        if segment.source.start >= segment.source.end
            || segment.source.end > segment.source.sequence_length
        {
            return Err(TruthContractError::InvalidSegmentInterval {
                identity: copy.identity.clone(),
            });
        }
        if segment.source.reference_id != copy.identity.reference_id
            || segment.source.sequence_id != copy.identity.sequence_id
            || !matches!(
                segment.source.source.kind,
                TruthSourceKind::Reference | TruthSourceKind::Founder
            )
        {
            return Err(TruthContractError::SegmentSourceMismatch {
                identity: copy.identity.clone(),
            });
        }
        require_digest(digests, &segment.source.source)?;
        for modification in &segment.modifications {
            let extent = ModificationExtent::from(modification);
            if !extent.contained(&segment.source) {
                return Err(TruthContractError::ModificationOutsideSegment {
                    event_id: extent.event_id.clone(),
                });
            }
            if let Some(prior) = previous {
                if prior.sort_key() >= extent.sort_key() {
                    return Err(TruthContractError::ModificationsNotOrdered {
                        event_id: extent.event_id.clone(),
                    });
                }
                if prior.overlaps(extent.clone()) {
                    return Err(TruthContractError::OverlappingModifications {
                        event_id: extent.event_id.clone(),
                    });
                }
            }
            validate_modification(
                modification,
                digests,
                groups,
                segment.source.sequence_length,
            )?;
            previous = Some(extent);
        }
    }
    Ok(())
}
fn validate_modification(
    modification: &TruthModification,
    digests: &BTreeMap<TruthSourceRef, &NamedDigest>,
    groups: &BTreeMap<TruthSourceRef, &NamedDigest>,
    sequence_length: u64,
) -> Result<(), TruthContractError> {
    let (event_id, provenance) = match modification {
        TruthModification::Snp {
            event_id,
            reference,
            alternate,
            provenance,
            ..
        } => {
            if reference == alternate {
                return Err(TruthContractError::SnpDoesNotChangeBase {
                    event_id: event_id.clone(),
                });
            }
            (event_id, provenance)
        }
        TruthModification::Insertion {
            event_id,
            point,
            bases,
            provenance,
        } => {
            if *point == 0
                || *point >= sequence_length
                || bases.is_empty()
                || !bases
                    .bytes()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
            {
                return Err(TruthContractError::InvalidInsertion {
                    event_id: event_id.clone(),
                });
            }
            (event_id, provenance)
        }
        TruthModification::Deletion {
            event_id,
            start,
            end,
            provenance,
        } => {
            if *start == 0 || *start >= *end || *end >= sequence_length {
                return Err(TruthContractError::InvalidDeletion {
                    event_id: event_id.clone(),
                });
            }
            (event_id, provenance)
        }
    };
    let expected = match provenance.kind {
        ModificationProvenanceKind::Explicit => TruthSourceKind::Explicit,
        ModificationProvenanceKind::Random => TruthSourceKind::Random,
        ModificationProvenanceKind::Matrix => TruthSourceKind::Matrix,
    };
    if provenance.version != "1"
        || provenance.source.kind != expected
        || (provenance.kind == ModificationProvenanceKind::Matrix) != provenance.matrix.is_some()
    {
        return Err(TruthContractError::InvalidModificationProvenance {
            event_id: event_id.clone(),
        });
    }
    if let Some(matrix) = &provenance.matrix
        && (matrix.matrix_source_id != provenance.source.source_id
            || matrix.population_id.is_empty()
            || matrix.variant_id.is_empty())
    {
        return Err(TruthContractError::InvalidModificationProvenance {
            event_id: event_id.clone(),
        });
    }
    if let Some(matrix) = &provenance.matrix {
        require_digest(
            groups,
            &TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: matrix.population_id.clone(),
            },
        )?;
    }
    require_digest(digests, &provenance.source)
}
#[derive(Clone)]
struct ModificationExtent {
    event_id: String,
    start: u64,
    end: u64,
    point: bool,
    deletion: bool,
}
impl ModificationExtent {
    fn from(modification: &TruthModification) -> Self {
        match modification {
            TruthModification::Snp {
                event_id, position, ..
            } => Self {
                event_id: event_id.clone(),
                start: *position,
                end: position.saturating_add(1),
                point: false,
                deletion: false,
            },
            TruthModification::Insertion {
                event_id, point, ..
            } => Self {
                event_id: event_id.clone(),
                start: *point,
                end: *point,
                point: true,
                deletion: false,
            },
            TruthModification::Deletion {
                event_id,
                start,
                end,
                ..
            } => Self {
                event_id: event_id.clone(),
                start: *start,
                end: *end,
                point: false,
                deletion: true,
            },
        }
    }
    fn contained(&self, source: &SegmentSource) -> bool {
        (!self.point
            && self.start < self.end
            && self.start >= source.start
            && self.end <= source.end)
            || (self.point
                && self.start == self.end
                && self.start >= source.start
                && self.start <= source.end)
    }
    const fn sort_key(&self) -> (u64, u64, bool) {
        (self.start, self.end, self.point)
    }
    fn overlaps(&self, other: Self) -> bool {
        if self.point && other.point {
            self.start == other.start
        } else if self.point && other.deletion {
            other.start < self.start && self.start <= other.end
        } else if other.point && self.deletion {
            self.start < other.start && other.start <= self.end
        } else if self.point || other.point {
            false
        } else {
            self.start < other.end && other.start < self.end
        }
    }
}
fn validate_genotype(
    g: &GenotypeRecord,
    header: &TruthHeader,
    copies: &BTreeSet<CopyIdentity>,
) -> Result<(), TruthContractError> {
    if g.alleles.is_empty() || g.alleles.len() != g.allele_dosages.len() {
        return Err(TruthContractError::EmptyGenotype {
            case_id: g.case_id.clone(),
        });
    }
    if g.reference_id != header.reference_id {
        return Err(TruthContractError::GenotypeReferenceMismatch {
            case_id: g.case_id.clone(),
        });
    }
    let mut seen = BTreeSet::new();
    let mut observed = vec![0_u8; g.alleles.len()];
    for c in &g.carriers {
        if c.allele_index >= g.alleles.len() {
            return Err(TruthContractError::CarrierAlleleOutOfRange {
                case_id: g.case_id.clone(),
                allele_index: c.allele_index,
            });
        }
        if c.copy.reference_id != g.reference_id
            || c.copy.individual_id != g.individual_id
            || c.copy.sequence_id != g.locus.sequence_id
        {
            return Err(TruthContractError::CarrierOutsideGenotype {
                case_id: g.case_id.clone(),
                copy: c.copy.clone(),
            });
        }
        if !copies.contains(&c.copy) {
            return Err(TruthContractError::UnknownCarrierCopy {
                case_id: g.case_id.clone(),
                copy: c.copy.clone(),
            });
        }
        if !seen.insert(c.copy.clone()) {
            return Err(TruthContractError::DuplicateCarrier {
                case_id: g.case_id.clone(),
                copy: c.copy.clone(),
            });
        }
        observed[c.allele_index] += 1;
    }
    for (index, (&declared, &actual)) in g.allele_dosages.iter().zip(&observed).enumerate() {
        if declared != actual {
            return Err(TruthContractError::DosageMismatch {
                case_id: g.case_id.clone(),
                allele_index: index,
                declared,
                actual,
            });
        }
    }
    Ok(())
}
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum RecordKey {
    Copy(CopyIdentity),
    Genotype(String, String, TruthLocus, String),
}
impl RecordKey {
    fn of(record: &TruthRecord) -> Self {
        match record {
            TruthRecord::HaplotypeCopy(c) => Self::Copy(c.identity.clone()),
            TruthRecord::Genotype(g) => Self::Genotype(
                g.reference_id.clone(),
                g.individual_id.clone(),
                g.locus.clone(),
                g.case_id.clone(),
            ),
        }
    }
}
