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
    /// Huella canónica del contenido semántico de `Reference`, distinta del
    /// SHA-256 del artefacto FASTA declarado entre `source_digests`.
    pub reference_content_sha256: String,
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
    /// Vista autoritativa de la secuencia vigente, que puede contener bases
    /// insertadas sin coordenada de fuente tras composiciones de indels.
    pub current_pieces: Vec<TruthCurrentPiece>,
    /// Historial cronológico de cambios aplicados en coordenadas vigentes.
    pub history: Vec<TruthHistoryEvent>,
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
pub enum TruthCurrentPiece {
    Source {
        source: SegmentSource,
        snps: Vec<TruthCurrentSnp>,
    },
    Inserted {
        bases: String,
        provenance: ModificationProvenance,
    },
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TruthCurrentSnp {
    pub event_id: String,
    pub position: u64,
    pub reference: String,
    pub alternate: String,
    pub provenance: ModificationProvenance,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TruthHistoryEvent {
    Snp {
        event_id: String,
        position: u64,
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
    /// `Snp` uses a one-base locus; `Indel` names the canonical source
    /// interval (or the zero-width insertion boundary) affected by the event.
    pub case_kind: GenotypeCaseKind,
    pub locus: TruthLocus,
    /// Estructura tipada y autoritativa del alelo. No se codifican estados
    /// estructurales como cadenas centinela dentro de un caso SNP.
    pub alleles: Vec<TruthAllele>,
    pub allele_dosages: Vec<u8>,
    pub carriers: Vec<AlleleCarrier>,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GenotypeCaseKind {
    Snp,
    Indel,
}
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum TruthAllele {
    Reference,
    Snp {
        base: String,
    },
    Insertion {
        event_id: String,
        bases: String,
        final_haplotype_bases: String,
        source_event: TruthStructuralEvent,
        piece_links: Vec<TruthPieceLink>,
    },
    Deletion {
        event_id: String,
        start: u64,
        end: u64,
        final_haplotype_bases: String,
        source_event: TruthStructuralEvent,
        piece_links: Vec<TruthPieceLink>,
    },
}
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum TruthStructuralEvent {
    Insertion {
        event_id: String,
        point: u64,
        bases: String,
    },
    Deletion {
        event_id: String,
        start: u64,
        end: u64,
    },
}
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct TruthPieceLink {
    pub piece_index: u64,
    pub start: u64,
    pub end: u64,
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
    InvalidLexicalField {
        field: &'static str,
        value: String,
    },
    EmptySourceDigests,
    DuplicateDigest {
        source: TruthSourceRef,
    },
    EmptySegment {
        identity: CopyIdentity,
    },
    EmptyCurrentPieces {
        identity: CopyIdentity,
    },
    InvalidCurrentPiece {
        identity: CopyIdentity,
    },
    DuplicateHistoryEvent {
        event_id: String,
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
    DuplicateAllele {
        case_id: String,
    },
    InvalidGenotypeAllele {
        case_id: String,
    },
    ArithmeticOverflow {
        field: &'static str,
    },
}

impl TruthDocument {
    pub fn validate_structure(&self) -> Result<(), TruthContractError> {
        let h = &self.header;
        validate_lexical_document(self)?;
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
        let (mut copy_count, mut genotype_count) = (0_u64, 0_u64);
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
                    copy_count = copy_count.checked_add(1).ok_or(
                        TruthContractError::ArithmeticOverflow {
                            field: "copy_count",
                        },
                    )?;
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
                    genotype_count = genotype_count.checked_add(1).ok_or(
                        TruthContractError::ArithmeticOverflow {
                            field: "genotype_count",
                        },
                    )?;
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

fn validate_lexical_document(document: &TruthDocument) -> Result<(), TruthContractError> {
    let h = &document.header;
    valid_id("execution_id", &h.execution_id)?;
    valid_id("reference_id", &h.reference_id)?;
    valid_sha(&h.reference_content_sha256)?;
    for digest in h.source_digests.iter().chain(&h.group_digests) {
        valid_source(&digest.source)?;
        valid_sha(&digest.sha256)?;
    }
    for record in &document.records {
        match record {
            TruthRecord::HaplotypeCopy(copy) => {
                for value in [
                    &copy.identity.reference_id,
                    &copy.identity.individual_id,
                    &copy.identity.sequence_id,
                ] {
                    valid_id("copy_identity", value)?;
                }
                if copy.identity.copy_index > 9 || copy.phase > 9 {
                    return Err(TruthContractError::InvalidLexicalField {
                        field: "copy_index_or_phase",
                        value: copy.identity.copy_index.to_string(),
                    });
                }
                if copy.provenance.origin.is_empty() {
                    return Err(TruthContractError::InvalidLexicalField {
                        field: "copy_provenance_origin",
                        value: copy.provenance.origin.clone(),
                    });
                }
                for source in &copy.provenance.sources {
                    valid_source(source)?;
                }
                for segment in &copy.segments {
                    valid_segment_source(&segment.source)?;
                    for modification in &segment.modifications {
                        valid_modification_lexical(modification)?;
                    }
                }
                for piece in &copy.current_pieces {
                    match piece {
                        TruthCurrentPiece::Source { source, snps } => {
                            valid_segment_source(source)?;
                            for snp in snps {
                                valid_id("event_id", &snp.event_id)?;
                                valid_base("current_snp_reference", &snp.reference)?;
                                valid_base("current_snp_alternate", &snp.alternate)?;
                                valid_provenance(&snp.provenance)?;
                            }
                        }
                        TruthCurrentPiece::Inserted { bases, provenance } => {
                            valid_bases("current_inserted_bases", bases)?;
                            valid_provenance(provenance)?
                        }
                    }
                }
                for event in &copy.history {
                    match event {
                        TruthHistoryEvent::Snp {
                            event_id,
                            provenance,
                            ..
                        }
                        | TruthHistoryEvent::Insertion {
                            event_id,
                            provenance,
                            ..
                        }
                        | TruthHistoryEvent::Deletion {
                            event_id,
                            provenance,
                            ..
                        } => {
                            valid_id("event_id", event_id)?;
                            valid_provenance(provenance)?;
                        }
                    }
                    if let TruthHistoryEvent::Insertion { bases, .. } = event {
                        valid_bases("history_insertion_bases", bases)?;
                    }
                }
            }
            TruthRecord::Genotype(genotype) => {
                for value in [
                    &genotype.case_id,
                    &genotype.reference_id,
                    &genotype.individual_id,
                    &genotype.locus.sequence_id,
                ] {
                    valid_id("genotype", value)?;
                }
                for allele in &genotype.alleles {
                    valid_allele_lexical(allele)?;
                }
                for carrier in &genotype.carriers {
                    for value in [
                        &carrier.copy.reference_id,
                        &carrier.copy.individual_id,
                        &carrier.copy.sequence_id,
                    ] {
                        valid_id("carrier_copy", value)?;
                    }
                    if carrier.copy.copy_index > 9 {
                        return Err(TruthContractError::InvalidLexicalField {
                            field: "carrier_copy_index",
                            value: carrier.copy.copy_index.to_string(),
                        });
                    }
                }
            }
        }
    }
    Ok(())
}
fn valid_id(field: &'static str, value: &str) -> Result<(), TruthContractError> {
    if value.len() <= 64
        && !value.is_empty()
        && value.as_bytes()[0].is_ascii_alphanumeric()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
    {
        Ok(())
    } else {
        Err(TruthContractError::InvalidLexicalField {
            field,
            value: value.to_owned(),
        })
    }
}
fn valid_sha(value: &str) -> Result<(), TruthContractError> {
    if value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'a'..=b'f'))
    {
        Ok(())
    } else {
        Err(TruthContractError::InvalidLexicalField {
            field: "sha256",
            value: value.to_owned(),
        })
    }
}
fn valid_source(source: &TruthSourceRef) -> Result<(), TruthContractError> {
    valid_id("source_id", &source.source_id)
}
fn valid_segment_source(source: &SegmentSource) -> Result<(), TruthContractError> {
    valid_source(&source.source)?;
    valid_id("segment_reference_id", &source.reference_id)?;
    valid_id("segment_sequence_id", &source.sequence_id)
}
fn valid_provenance(provenance: &ModificationProvenance) -> Result<(), TruthContractError> {
    valid_source(&provenance.source)?;
    if let Some(matrix) = &provenance.matrix {
        valid_id("matrix_source_id", &matrix.matrix_source_id)?;
        valid_id("population_id", &matrix.population_id)?;
        valid_id("variant_id", &matrix.variant_id)?;
    }
    Ok(())
}
fn valid_modification_lexical(modification: &TruthModification) -> Result<(), TruthContractError> {
    match modification {
        TruthModification::Snp {
            event_id,
            reference,
            alternate,
            provenance,
            ..
        } => {
            valid_id("event_id", event_id)?;
            valid_base("historical_snp_reference", reference)?;
            valid_base("historical_snp_alternate", alternate)?;
            valid_provenance(provenance)
        }
        TruthModification::Insertion {
            event_id,
            bases,
            provenance,
            ..
        } => {
            valid_id("event_id", event_id)?;
            valid_bases("historical_insertion_bases", bases)?;
            valid_provenance(provenance)
        }
        TruthModification::Deletion {
            event_id,
            provenance,
            ..
        } => {
            valid_id("event_id", event_id)?;
            valid_provenance(provenance)
        }
    }
}
fn valid_allele_lexical(allele: &TruthAllele) -> Result<(), TruthContractError> {
    match allele {
        TruthAllele::Reference => Ok(()),
        TruthAllele::Snp { base } => valid_base("snp_allele_base", base),
        TruthAllele::Insertion {
            event_id,
            bases,
            final_haplotype_bases,
            source_event,
            piece_links,
        } => {
            valid_id("indel_event_id", event_id)?;
            valid_bases("insertion_allele_bases", bases)?;
            valid_bases("final_haplotype_bases", final_haplotype_bases)?;
            let final_length = u64::try_from(final_haplotype_bases.len()).map_err(|_| {
                TruthContractError::ArithmeticOverflow {
                    field: "final_haplotype_bases_length",
                }
            })?;
            if !matches!(source_event, TruthStructuralEvent::Insertion { point, bases: event_bases, .. } if *point > 0 && event_bases == bases)
                || piece_links.is_empty()
                || piece_links.iter().any(|link| link.end > final_length)
            {
                return Err(TruthContractError::InvalidGenotypeAllele {
                    case_id: event_id.clone(),
                });
            }
            validate_structural_event(event_id, source_event, piece_links)
        }
        TruthAllele::Deletion {
            event_id,
            start,
            end,
            final_haplotype_bases,
            source_event,
            piece_links,
            ..
        } => {
            valid_id("indel_event_id", event_id)?;
            valid_bases("final_haplotype_bases", final_haplotype_bases)?;
            let final_length = u64::try_from(final_haplotype_bases.len()).map_err(|_| {
                TruthContractError::ArithmeticOverflow {
                    field: "final_haplotype_bases_length",
                }
            })?;
            if *start >= 1 && *end >= 2 && start < end {
                if piece_links.iter().any(|link| link.end > final_length) {
                    return Err(TruthContractError::InvalidGenotypeAllele {
                        case_id: event_id.clone(),
                    });
                }
                if !matches!(source_event, TruthStructuralEvent::Deletion { start: source_start, end: source_end, .. } if source_start == start && source_end == end && *source_start >= 1 && *source_end >= 2)
                {
                    return Err(TruthContractError::InvalidGenotypeAllele {
                        case_id: event_id.clone(),
                    });
                }
                validate_structural_event(event_id, source_event, piece_links)
            } else {
                Err(TruthContractError::InvalidLexicalField {
                    field: "deletion_allele_interval",
                    value: format!("{start}:{end}"),
                })
            }
        }
    }
}
fn validate_structural_event(
    allele_event_id: &str,
    source_event: &TruthStructuralEvent,
    piece_links: &[TruthPieceLink],
) -> Result<(), TruthContractError> {
    let event_id = match source_event {
        TruthStructuralEvent::Insertion {
            event_id,
            point,
            bases,
        } => {
            if *point == 0 {
                return Err(TruthContractError::InvalidLexicalField {
                    field: "source_event_insertion_point",
                    value: point.to_string(),
                });
            }
            valid_bases("source_event_insertion_bases", bases)?;
            event_id
        }
        TruthStructuralEvent::Deletion {
            event_id,
            start,
            end,
        } => {
            if *start == 0 || *end < 2 || start >= end {
                return Err(TruthContractError::InvalidLexicalField {
                    field: "source_event_deletion_interval",
                    value: format!("{start}:{end}"),
                });
            }
            event_id
        }
    };
    valid_id("source_event_id", event_id)?;
    if event_id != allele_event_id {
        return Err(TruthContractError::InvalidLexicalField {
            field: "source_event_id",
            value: event_id.clone(),
        });
    }
    let mut prior_end = 0;
    let mut piece_indices = BTreeSet::new();
    for link in piece_links {
        if link.start >= link.end
            || link.start < prior_end
            || !piece_indices.insert(link.piece_index)
        {
            return Err(TruthContractError::InvalidLexicalField {
                field: "piece_links",
                value: format!("{}:{}:{}", link.piece_index, link.start, link.end),
            });
        }
        prior_end = link.end;
    }
    Ok(())
}
fn valid_base(field: &'static str, value: &str) -> Result<(), TruthContractError> {
    if matches!(value, "A" | "C" | "G" | "T") {
        Ok(())
    } else {
        Err(TruthContractError::InvalidLexicalField {
            field,
            value: value.to_owned(),
        })
    }
}
fn valid_bases(field: &'static str, value: &str) -> Result<(), TruthContractError> {
    if !value.is_empty()
        && value
            .bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
    {
        Ok(())
    } else {
        Err(TruthContractError::InvalidLexicalField {
            field,
            value: value.to_owned(),
        })
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
    let mut provenance_sources = BTreeSet::new();
    for source in &copy.provenance.sources {
        if !provenance_sources.insert(source.clone()) {
            return Err(TruthContractError::InvalidLexicalField {
                field: "copy_provenance_sources",
                value: source.source_id.clone(),
            });
        }
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
            let extent = ModificationExtent::from(modification)?;
            if !extent.contained(&segment.source) {
                return Err(TruthContractError::ModificationOutsideSegment {
                    event_id: extent.event_id.clone(),
                });
            }
            if let Some(prior) = previous {
                let deletion_covers_snp = prior.deletion != extent.deletion
                    && ((prior.deletion && prior.start <= extent.start && extent.end <= prior.end)
                        || (extent.deletion
                            && extent.start <= prior.start
                            && prior.end <= extent.end));
                if prior.sort_key() >= extent.sort_key() && !deletion_covers_snp {
                    return Err(TruthContractError::ModificationsNotOrdered {
                        event_id: extent.event_id.clone(),
                    });
                }
                if prior.overlaps(extent.clone()) && !deletion_covers_snp {
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
    validate_current_state(copy, digests, groups)?;
    Ok(())
}

fn validate_current_state(
    copy: &HaplotypeCopyRecord,
    digests: &BTreeMap<TruthSourceRef, &NamedDigest>,
    groups: &BTreeMap<TruthSourceRef, &NamedDigest>,
) -> Result<(), TruthContractError> {
    if copy.current_pieces.is_empty() {
        return Err(TruthContractError::EmptyCurrentPieces {
            identity: copy.identity.clone(),
        });
    }
    for piece in &copy.current_pieces {
        match piece {
            TruthCurrentPiece::Source { source, snps } => {
                if source.start >= source.end
                    || source.end > source.sequence_length
                    || source.reference_id != copy.identity.reference_id
                    || source.sequence_id != copy.identity.sequence_id
                    || !matches!(
                        source.source.kind,
                        TruthSourceKind::Reference | TruthSourceKind::Founder
                    )
                {
                    return Err(TruthContractError::InvalidCurrentPiece {
                        identity: copy.identity.clone(),
                    });
                }
                require_digest(digests, &source.source)?;
                let mut prior = None;
                for snp in snps {
                    if snp.position < source.start
                        || snp.position >= source.end
                        || snp.reference == snp.alternate
                    {
                        return Err(TruthContractError::InvalidCurrentPiece {
                            identity: copy.identity.clone(),
                        });
                    }
                    if prior.is_some_and(|position| position >= snp.position) {
                        return Err(TruthContractError::InvalidCurrentPiece {
                            identity: copy.identity.clone(),
                        });
                    }
                    validate_modification_provenance(
                        &snp.event_id,
                        &snp.provenance,
                        digests,
                        groups,
                    )?;
                    prior = Some(snp.position);
                }
            }
            TruthCurrentPiece::Inserted { bases, provenance } => {
                if bases.is_empty()
                    || !bases
                        .bytes()
                        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
                {
                    return Err(TruthContractError::InvalidCurrentPiece {
                        identity: copy.identity.clone(),
                    });
                }
                validate_modification_provenance("inserted_piece", provenance, digests, groups)?;
            }
        }
    }
    let mut history = BTreeSet::new();
    for event in &copy.history {
        let (event_id, valid, provenance) = match event {
            TruthHistoryEvent::Snp {
                event_id,
                provenance,
                ..
            } => (event_id, true, provenance),
            TruthHistoryEvent::Insertion {
                event_id,
                point,
                bases,
                provenance,
            } => (
                event_id,
                *point > 0
                    && !bases.is_empty()
                    && bases
                        .bytes()
                        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N')),
                provenance,
            ),
            TruthHistoryEvent::Deletion {
                event_id,
                start,
                end,
                provenance,
            } => (event_id, *start > 0 && start < end, provenance),
        };
        if !valid {
            return Err(TruthContractError::InvalidCurrentPiece {
                identity: copy.identity.clone(),
            });
        }
        if !history.insert(event_id) {
            return Err(TruthContractError::DuplicateHistoryEvent {
                event_id: event_id.clone(),
            });
        }
        validate_modification_provenance(event_id, provenance, digests, groups)?;
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
    validate_modification_provenance(event_id, provenance, digests, groups)
}

fn validate_modification_provenance(
    event_id: &str,
    provenance: &ModificationProvenance,
    digests: &BTreeMap<TruthSourceRef, &NamedDigest>,
    groups: &BTreeMap<TruthSourceRef, &NamedDigest>,
) -> Result<(), TruthContractError> {
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
            event_id: event_id.to_owned(),
        });
    }
    if let Some(matrix) = &provenance.matrix
        && (matrix.matrix_source_id != provenance.source.source_id
            || matrix.population_id.is_empty()
            || matrix.variant_id.is_empty())
    {
        return Err(TruthContractError::InvalidModificationProvenance {
            event_id: event_id.to_owned(),
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
    fn from(modification: &TruthModification) -> Result<Self, TruthContractError> {
        Ok(match modification {
            TruthModification::Snp {
                event_id, position, ..
            } => Self {
                event_id: event_id.clone(),
                start: *position,
                end: position
                    .checked_add(1)
                    .ok_or(TruthContractError::ArithmeticOverflow {
                        field: "snp_position_plus_one",
                    })?,
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
        })
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
pub(crate) fn validate_genotype(
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
    if g.locus.start > g.locus.end
        || (matches!(g.case_kind, GenotypeCaseKind::Snp)
            && g.locus.end
                != g.locus
                    .start
                    .checked_add(1)
                    .ok_or(TruthContractError::ArithmeticOverflow {
                        field: "snp_locus_start_plus_one",
                    })?)
        || (matches!(g.case_kind, GenotypeCaseKind::Indel)
            && !g.alleles.iter().any(|allele| {
                matches!(
                    allele,
                    TruthAllele::Insertion { .. } | TruthAllele::Deletion { .. }
                )
            }))
    {
        return Err(TruthContractError::InvalidGenotypeAllele {
            case_id: g.case_id.clone(),
        });
    }
    let mut unique_alleles = BTreeSet::new();
    for allele in &g.alleles {
        if !unique_alleles.insert(allele) {
            return Err(TruthContractError::DuplicateAllele {
                case_id: g.case_id.clone(),
            });
        }
        match (g.case_kind, allele) {
            (GenotypeCaseKind::Snp, TruthAllele::Reference | TruthAllele::Snp { .. })
            | (
                GenotypeCaseKind::Indel,
                TruthAllele::Reference
                | TruthAllele::Insertion { .. }
                | TruthAllele::Deletion { .. },
            ) => {}
            _ => {
                return Err(TruthContractError::InvalidGenotypeAllele {
                    case_id: g.case_id.clone(),
                });
            }
        }
        match allele {
            TruthAllele::Insertion {
                source_event: TruthStructuralEvent::Insertion { point, .. },
                ..
            } if g.locus.start != *point || g.locus.end != *point => {
                return Err(TruthContractError::InvalidGenotypeAllele {
                    case_id: g.case_id.clone(),
                });
            }
            TruthAllele::Deletion {
                source_event: TruthStructuralEvent::Deletion { start, end, .. },
                ..
            } if g.locus.start != *start || g.locus.end != *end => {
                return Err(TruthContractError::InvalidGenotypeAllele {
                    case_id: g.case_id.clone(),
                });
            }
            _ => {}
        }
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
        observed[c.allele_index] = observed[c.allele_index].checked_add(1).ok_or(
            TruthContractError::ArithmeticOverflow {
                field: "observed_allele_dosage",
            },
        )?;
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
