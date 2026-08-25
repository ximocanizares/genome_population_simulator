//! Agrupación multialélica y genotipos faseados de la vista VCF F2.5.3.
//!
//! Recibe los alelos ya normalizados por `vcf_projection`; no escribe VCF ni
//! decide el tratamiento de `N` o del complemento, que corresponden a los
//! pasos siguientes. La salida conserva los índices de alelo VCF (`0` para
//! REF, `1..` para ALT), la dosis por muestra y la relación ALT↔casos.

use std::collections::{BTreeMap, BTreeSet};

use crate::reference::Reference;
use crate::truth_schema::{CopyIdentity, TruthDocument, TruthRecord};
use crate::truth_validation::{TruthValidationError, validate_and_reconstruct_truth};
use crate::vcf_projection::NormalizedAllele;

/// Genotipo real de una muestra en el orden de fase de sus copias.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PhasedGenotype {
    /// Índices VCF: `0` es REF y `1..` son los ALT del registro.
    pub alleles: Vec<u32>,
    /// Número de copias de REF seguido de cada ALT del registro.
    pub allele_dosages: Vec<u8>,
}

/// Registro lógico listo para que F2.5.5 lo materialice como una línea VCF.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GroupedVcfRecord {
    /// Coordenada interna zero-based; el escritor será el único que use POS
    /// one-based.
    pub sequence_id: String,
    pub position: u64,
    pub reference: String,
    pub alternates: Vec<String>,
    /// Muestras en orden canónico por identificador estable.
    pub samples: Vec<(String, PhasedGenotype)>,
    /// `GPS_CASE`: casos fuente únicos en orden canónico.
    pub gps_cases: Vec<String>,
    /// Relación explícita, paralela a `alternates`, ALT↔casos fuente.
    pub alternate_cases: Vec<Vec<String>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum VcfGroupingError {
    Truth(TruthValidationError),
    UnknownCopy {
        copy: CopyIdentity,
    },
    CopySequenceMismatch {
        copy: CopyIdentity,
        sequence_id: String,
    },
    NonContiguousCopyPhases {
        individual_id: String,
        sequence_id: String,
        observed: Vec<u8>,
    },
    DuplicateCopyAtLocus {
        copy: CopyIdentity,
        sequence_id: String,
        position: u64,
    },
    ReferenceOutOfBounds {
        sequence_id: String,
        position: u64,
    },
    ReferenceMismatch {
        sequence_id: String,
        position: u64,
    },
    EmptyAllele {
        sequence_id: String,
        position: u64,
    },
    InvalidSourceSpan {
        sequence_id: String,
        position: u64,
        start: u64,
        end: u64,
    },
    ReferenceEqualsAlternate {
        sequence_id: String,
        position: u64,
    },
    IncompatibleReferences {
        sequence_id: String,
        position: u64,
    },
    CaseAlleleNotProjected {
        case_id: String,
        copy: CopyIdentity,
    },
    AlternateWithoutCase {
        sequence_id: String,
        position: u64,
        alternate: String,
    },
    ArithmeticOverflow,
}

impl std::fmt::Display for VcfGroupingError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for VcfGroupingError {}

/// Agrupa los alelos normalizados en un único registro F2 por `CHROM+POS`.
///
/// El orden físico de entrada no es semántico: se ordenan registros, ALT,
/// muestras y casos canónicamente. Antes de producir nada se vuelve a validar
/// `truth.jsonl`, para que una copia o un caso desconectado nunca se convierta
/// silenciosamente en un `GT` o `GPS_CASE` falso.
pub fn group_normalized_alleles(
    truth: &TruthDocument,
    reference: &Reference,
    normalized: &[NormalizedAllele],
) -> Result<Vec<GroupedVcfRecord>, VcfGroupingError> {
    validate_and_reconstruct_truth(truth.clone(), reference).map_err(VcfGroupingError::Truth)?;

    let sequence_order = reference
        .sequences()
        .iter()
        .enumerate()
        .map(|(index, sequence)| (sequence.id().as_str().to_owned(), index))
        .collect::<BTreeMap<_, _>>();
    let copies = truth
        .records
        .iter()
        .filter_map(|record| match record {
            TruthRecord::HaplotypeCopy(copy) => Some((copy.identity.clone(), copy.phase)),
            TruthRecord::Genotype(_) => None,
        })
        .collect::<BTreeMap<_, _>>();
    validate_contiguous_copy_phases(&copies)?;

    let mut grouped = BTreeMap::<(usize, String, u64), Vec<&NormalizedAllele>>::new();
    for allele in normalized {
        let sequence_index = *sequence_order.get(&allele.sequence_id).ok_or_else(|| {
            VcfGroupingError::ReferenceOutOfBounds {
                sequence_id: allele.sequence_id.clone(),
                position: allele.position,
            }
        })?;
        let Some(_) = copies.get(&allele.copy) else {
            return Err(VcfGroupingError::UnknownCopy {
                copy: allele.copy.clone(),
            });
        };
        if allele.copy.sequence_id != allele.sequence_id {
            return Err(VcfGroupingError::CopySequenceMismatch {
                copy: allele.copy.clone(),
                sequence_id: allele.sequence_id.clone(),
            });
        }
        validate_normalized_allele(allele, reference)?;
        grouped
            .entry((sequence_index, allele.sequence_id.clone(), allele.position))
            .or_default()
            .push(allele);
    }

    let mut records = Vec::with_capacity(grouped.len());
    let mut carriers_by_record = Vec::<BTreeMap<CopyIdentity, ProjectedCarrier>>::new();
    for ((_sequence_index, sequence_id, position), alleles) in grouped {
        let (record, carriers) = build_record(sequence_id.clone(), position, alleles, &copies)?;
        carriers_by_record.push(carriers);
        records.push(record);
    }

    // La asociación se deduce de los intervalos fuente finales de la copia,
    // emitidos por F2.5.2 después de componer sus cambios. Nunca se infiere
    // renormalizando un `TruthAllele` aislado: eso perdería los casos de los
    // cambios que forman un ALT compuesto.
    for genotype in truth.records.iter().filter_map(|record| match record {
        TruthRecord::Genotype(genotype) => Some(genotype),
        TruthRecord::HaplotypeCopy(_) => None,
    }) {
        for carrier in genotype
            .carriers
            .iter()
            .filter(|carrier| carrier.allele_index != 0)
        {
            let mut mapped = false;
            for (record_index, carriers) in carriers_by_record.iter().enumerate() {
                let Some(projected) = carriers.get(&carrier.copy) else {
                    continue;
                };
                if records[record_index].sequence_id == genotype.locus.sequence_id
                    && source_span_contains_locus(
                        projected,
                        genotype.locus.start,
                        genotype.locus.end,
                    )
                {
                    records[record_index].alternate_cases[projected.alt_index]
                        .push(genotype.case_id.clone());
                    mapped = true;
                }
            }
            if !mapped {
                return Err(VcfGroupingError::CaseAlleleNotProjected {
                    case_id: genotype.case_id.clone(),
                    copy: carrier.copy.clone(),
                });
            }
        }
    }

    for record in &mut records {
        for cases in &mut record.alternate_cases {
            cases.sort();
            cases.dedup();
        }
        record.gps_cases = record
            .alternate_cases
            .iter()
            .flatten()
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        for (alternate, cases) in record.alternates.iter().zip(&record.alternate_cases) {
            if cases.is_empty() {
                return Err(VcfGroupingError::AlternateWithoutCase {
                    sequence_id: record.sequence_id.clone(),
                    position: record.position,
                    alternate: alternate.clone(),
                });
            }
        }
    }
    Ok(records)
}

fn validate_normalized_allele(
    allele: &NormalizedAllele,
    reference: &Reference,
) -> Result<(), VcfGroupingError> {
    if allele.reference.is_empty() || allele.alternate.is_empty() {
        return Err(VcfGroupingError::EmptyAllele {
            sequence_id: allele.sequence_id.clone(),
            position: allele.position,
        });
    }
    if allele.source_start > allele.source_end {
        return Err(VcfGroupingError::InvalidSourceSpan {
            sequence_id: allele.sequence_id.clone(),
            position: allele.position,
            start: allele.source_start,
            end: allele.source_end,
        });
    }
    if allele.reference == allele.alternate {
        return Err(VcfGroupingError::ReferenceEqualsAlternate {
            sequence_id: allele.sequence_id.clone(),
            position: allele.position,
        });
    }
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id().as_str() == allele.sequence_id)
        .ok_or_else(|| VcfGroupingError::ReferenceOutOfBounds {
            sequence_id: allele.sequence_id.clone(),
            position: allele.position,
        })?;
    let end = allele
        .position
        .checked_add(
            u64::try_from(allele.reference.len())
                .map_err(|_| VcfGroupingError::ArithmeticOverflow)?,
        )
        .ok_or(VcfGroupingError::ArithmeticOverflow)?;
    let start = usize::try_from(allele.position).ok();
    let end_usize = usize::try_from(end).ok();
    let Some((start, end)) = start
        .zip(end_usize)
        .filter(|(_, end)| *end <= sequence.sequence().as_bytes().len())
    else {
        return Err(VcfGroupingError::ReferenceOutOfBounds {
            sequence_id: allele.sequence_id.clone(),
            position: allele.position,
        });
    };
    if sequence.sequence().as_bytes()[start..end] != *allele.reference.as_bytes() {
        return Err(VcfGroupingError::ReferenceMismatch {
            sequence_id: allele.sequence_id.clone(),
            position: allele.position,
        });
    }
    Ok(())
}

#[derive(Clone, Copy, Debug)]
struct ProjectedCarrier {
    alt_index: usize,
    source_start: u64,
    source_end: u64,
}

fn validate_contiguous_copy_phases(
    copies: &BTreeMap<CopyIdentity, u8>,
) -> Result<(), VcfGroupingError> {
    let mut phases = BTreeMap::<(String, String), Vec<u8>>::new();
    for (copy, phase) in copies {
        phases
            .entry((copy.individual_id.clone(), copy.sequence_id.clone()))
            .or_default()
            .push(*phase);
    }
    for ((individual_id, sequence_id), mut observed) in phases {
        observed.sort_unstable();
        let valid = (1..=10).contains(&observed.len())
            && observed
                .iter()
                .enumerate()
                .all(|(expected, phase)| usize::from(*phase) == expected);
        if !valid {
            return Err(VcfGroupingError::NonContiguousCopyPhases {
                individual_id,
                sequence_id,
                observed,
            });
        }
    }
    Ok(())
}

fn source_span_contains_locus(projected: &ProjectedCarrier, start: u64, end: u64) -> bool {
    if start == end {
        projected.source_start <= start && start <= projected.source_end
    } else {
        projected.source_start < end && start < projected.source_end
    }
}

fn build_record(
    sequence_id: String,
    position: u64,
    alleles: Vec<&NormalizedAllele>,
    copies: &BTreeMap<CopyIdentity, u8>,
) -> Result<(GroupedVcfRecord, BTreeMap<CopyIdentity, ProjectedCarrier>), VcfGroupingError> {
    let common_reference = alleles
        .iter()
        .map(|allele| allele.reference.as_str())
        .max_by_key(|reference| reference.len())
        .expect("a grouping key always has an allele");
    if alleles
        .iter()
        .any(|allele| !common_reference.starts_with(&allele.reference))
    {
        return Err(VcfGroupingError::IncompatibleReferences {
            sequence_id,
            position,
        });
    }
    let alternate_by_copy = alleles
        .iter()
        .map(|allele| {
            let suffix = &common_reference[allele.reference.len()..];
            (allele.copy.clone(), format!("{}{suffix}", allele.alternate))
        })
        .collect::<BTreeMap<_, _>>();
    if alternate_by_copy.len() != alleles.len() {
        let mut seen = BTreeSet::new();
        let duplicate = alleles
            .iter()
            .map(|allele| allele.copy.clone())
            .find(|copy| !seen.insert(copy.clone()))
            .expect("a duplicate exists");
        return Err(VcfGroupingError::DuplicateCopyAtLocus {
            copy: duplicate,
            sequence_id,
            position,
        });
    }
    let alternates = alternate_by_copy
        .values()
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let alt_index = alternates
        .iter()
        .enumerate()
        .map(|(index, alternate)| (alternate.clone(), index))
        .collect::<BTreeMap<_, _>>();
    let copy_to_alt = alternate_by_copy
        .iter()
        .map(|(copy, alternate)| {
            let source = alleles
                .iter()
                .find(|allele| allele.copy == *copy)
                .expect("each copy is present in the normalized input");
            (
                copy.clone(),
                ProjectedCarrier {
                    alt_index: alt_index[alternate],
                    source_start: source.source_start,
                    source_end: source.source_end,
                },
            )
        })
        .collect::<BTreeMap<_, _>>();

    let mut sample_copies = BTreeMap::<String, Vec<(u8, CopyIdentity)>>::new();
    for (copy, phase) in copies {
        if copy.sequence_id == sequence_id {
            sample_copies
                .entry(copy.individual_id.clone())
                .or_default()
                .push((*phase, copy.clone()));
        }
    }
    let samples = sample_copies
        .into_iter()
        .map(|(sample, mut sample_copies)| {
            sample_copies.sort_by_key(|(phase, copy)| (*phase, copy.clone()));
            let mut alleles = Vec::with_capacity(sample_copies.len());
            let mut dosage = vec![0_u8; alternates.len() + 1];
            for (_, copy) in sample_copies {
                let allele = copy_to_alt
                    .get(&copy)
                    .map_or(0_usize, |projected| projected.alt_index + 1);
                dosage[allele] = dosage[allele]
                    .checked_add(1)
                    .ok_or(VcfGroupingError::ArithmeticOverflow)?;
                alleles
                    .push(u32::try_from(allele).map_err(|_| VcfGroupingError::ArithmeticOverflow)?);
            }
            Ok((
                sample,
                PhasedGenotype {
                    alleles,
                    allele_dosages: dosage,
                },
            ))
        })
        .collect::<Result<Vec<_>, VcfGroupingError>>()?;
    let alternate_cases = vec![Vec::new(); alternates.len()];
    Ok((
        GroupedVcfRecord {
            sequence_id,
            position,
            reference: common_reference.to_owned(),
            alternates,
            samples,
            gps_cases: Vec::new(),
            alternate_cases,
        },
        copy_to_alt,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::identity::{ReferenceId, SequenceId};
    use crate::reference::{
        DnaSequence, FastaHeader, ReferenceOrigin, ReferenceSequence, SequenceClass,
    };
    use crate::truth_schema::{
        AlleleCarrier, CopyProvenance, GenotypeCaseKind, GenotypeRecord, HaplotypeCopyRecord,
        ModificationProvenance, ModificationProvenanceKind, NamedDigest, SegmentSource,
        TRUTH_COORDINATE_CONVENTION, TRUTH_DOCUMENT_KIND, TRUTH_SCHEMA_VERSION, TruthAllele,
        TruthCurrentPiece, TruthCurrentSnp, TruthHeader, TruthHistoryEvent, TruthLocus,
        TruthModification, TruthPieceLink, TruthSourceKind, TruthSourceRef, TruthStructuralEvent,
    };
    use crate::truth_validation::reference_content_sha256;

    fn reference() -> Reference {
        let sequence = ReferenceSequence::new(
            SequenceId::parse("chr1", "test").expect("sequence"),
            FastaHeader::parse(b"chr1").expect("header"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(b"ACGT").expect("DNA"),
        )
        .expect("reference sequence");
        Reference::new(
            ReferenceId::parse("REF_main", "test").expect("reference"),
            ReferenceOrigin::Generated,
            vec![sequence],
        )
        .expect("reference")
    }

    fn long_reference() -> Reference {
        let sequence = ReferenceSequence::new(
            SequenceId::parse("chr1", "test").expect("sequence"),
            FastaHeader::parse(b"chr1").expect("header"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(b"CAAAAAAA").expect("DNA"),
        )
        .expect("reference sequence");
        Reference::new(
            ReferenceId::parse("REF_main", "test").expect("reference"),
            ReferenceOrigin::Generated,
            vec![sequence],
        )
        .expect("reference")
    }

    fn identity(individual: &str, index: u8) -> CopyIdentity {
        CopyIdentity {
            reference_id: "REF_main".into(),
            individual_id: individual.into(),
            sequence_id: "chr1".into(),
            copy_index: index,
        }
    }

    fn provenance() -> ModificationProvenance {
        ModificationProvenance {
            kind: ModificationProvenanceKind::Explicit,
            version: "1".into(),
            source: TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "REF_main".into(),
            },
            matrix: None,
        }
    }

    fn source() -> SegmentSource {
        SegmentSource {
            source: TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "REF_main".into(),
            },
            reference_id: "REF_main".into(),
            sequence_id: "chr1".into(),
            start: 0,
            end: 4,
            sequence_length: 4,
        }
    }

    fn long_source(start: u64, end: u64) -> SegmentSource {
        SegmentSource {
            source: TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "REF_main".into(),
            },
            reference_id: "REF_main".into(),
            sequence_id: "chr1".into(),
            start,
            end,
            sequence_length: 8,
        }
    }

    fn copy(individual: &str, index: u8, alternate: Option<&str>) -> HaplotypeCopyRecord {
        let snps = alternate
            .map(|alternate| {
                vec![TruthCurrentSnp {
                    event_id: format!("SNP_{individual}"),
                    position: 1,
                    reference: "C".into(),
                    alternate: alternate.into(),
                    provenance: provenance(),
                }]
            })
            .unwrap_or_default();
        let modifications = alternate
            .map(|alternate| {
                vec![TruthModification::Snp {
                    event_id: format!("SNP_{individual}"),
                    position: 1,
                    reference: "C".into(),
                    alternate: alternate.into(),
                    provenance: provenance(),
                }]
            })
            .unwrap_or_default();
        let history = alternate
            .map(|_| {
                vec![crate::truth_schema::TruthHistoryEvent::Snp {
                    event_id: format!("SNP_{individual}"),
                    position: 1,
                    provenance: provenance(),
                }]
            })
            .unwrap_or_default();
        HaplotypeCopyRecord {
            identity: identity(individual, index),
            phase: index,
            provenance: CopyProvenance {
                origin: "test".into(),
                sources: vec![TruthSourceRef {
                    kind: TruthSourceKind::Reference,
                    source_id: "REF_main".into(),
                }],
            },
            segments: vec![crate::truth_schema::TruthSegment {
                source: source(),
                modifications,
            }],
            current_pieces: vec![TruthCurrentPiece::Source {
                source: source(),
                snps,
            }],
            history,
        }
    }

    fn document(reference: &Reference, alternate_b: &str) -> TruthDocument {
        let a = copy("IND_A", 0, Some("T"));
        let b = copy("IND_B", 0, Some(alternate_b));
        let genotype = |individual: &str, alternate: &str, case_id: &str| GenotypeRecord {
            case_id: case_id.into(),
            reference_id: "REF_main".into(),
            individual_id: individual.into(),
            case_kind: GenotypeCaseKind::Snp,
            locus: TruthLocus {
                sequence_id: "chr1".into(),
                start: 1,
                end: 2,
            },
            alleles: vec![
                TruthAllele::Reference,
                TruthAllele::Snp {
                    base: alternate.into(),
                },
            ],
            allele_dosages: vec![0, 1],
            carriers: vec![AlleleCarrier {
                allele_index: 1,
                copy: identity(individual, 0),
            }],
        };
        TruthDocument {
            header: TruthHeader {
                document_kind: TRUTH_DOCUMENT_KIND.into(),
                schema_version: TRUTH_SCHEMA_VERSION.into(),
                execution_id: "RUN_01".into(),
                reference_id: "REF_main".into(),
                reference_content_sha256: reference_content_sha256(reference),
                coordinate_convention: TRUTH_COORDINATE_CONVENTION.into(),
                copy_count: 2,
                genotype_count: 2,
                source_digests: vec![NamedDigest {
                    source: TruthSourceRef {
                        kind: TruthSourceKind::Reference,
                        source_id: "REF_main".into(),
                    },
                    sha256: "a".repeat(64),
                }],
                group_digests: vec![],
            },
            records: vec![
                TruthRecord::HaplotypeCopy(a),
                TruthRecord::HaplotypeCopy(b),
                TruthRecord::Genotype(genotype("IND_A", "T", "CASE_A")),
                TruthRecord::Genotype(genotype("IND_B", alternate_b, "CASE_B")),
            ],
        }
    }

    #[test]
    fn groups_multiple_alts_cases_and_reordered_projection_canonically() {
        let reference = reference();
        let truth = document(&reference, "G");
        let projected = vec![
            NormalizedAllele {
                sequence_id: "chr1".into(),
                position: 1,
                reference: "C".into(),
                alternate: "G".into(),
                source_start: 1,
                source_end: 2,
                copy: identity("IND_B", 0),
            },
            NormalizedAllele {
                sequence_id: "chr1".into(),
                position: 1,
                reference: "C".into(),
                alternate: "T".into(),
                source_start: 1,
                source_end: 2,
                copy: identity("IND_A", 0),
            },
        ];

        let records = group_normalized_alleles(&truth, &reference, &projected).expect("grouping");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].alternates, ["G", "T"]);
        assert_eq!(records[0].gps_cases, ["CASE_A", "CASE_B"]);
        assert_eq!(records[0].alternate_cases, [vec!["CASE_B"], vec!["CASE_A"]]);
        assert_eq!(records[0].samples[0].1.alleles, [2]);
        assert_eq!(records[0].samples[1].1.alleles, [1]);
    }

    #[test]
    fn keeps_all_source_cases_when_one_alt_has_multiple_cases() {
        let reference = reference();
        let truth = document(&reference, "T");
        let projected = vec![
            NormalizedAllele {
                sequence_id: "chr1".into(),
                position: 1,
                reference: "C".into(),
                alternate: "T".into(),
                source_start: 1,
                source_end: 2,
                copy: identity("IND_A", 0),
            },
            NormalizedAllele {
                sequence_id: "chr1".into(),
                position: 1,
                reference: "C".into(),
                alternate: "T".into(),
                source_start: 1,
                source_end: 2,
                copy: identity("IND_B", 0),
            },
        ];

        let records = group_normalized_alleles(&truth, &reference, &projected).expect("grouping");
        assert_eq!(records[0].alternates, ["T"]);
        assert_eq!(records[0].alternate_cases, [vec!["CASE_A", "CASE_B"]]);
        assert_eq!(records[0].gps_cases, ["CASE_A", "CASE_B"]);
    }

    #[test]
    fn keeps_all_cases_of_a_composed_final_allele() {
        let reference = long_reference();
        let copy = HaplotypeCopyRecord {
            identity: identity("IND", 0),
            phase: 0,
            provenance: CopyProvenance {
                origin: "test".into(),
                sources: vec![TruthSourceRef {
                    kind: TruthSourceKind::Reference,
                    source_id: "REF_main".into(),
                }],
            },
            segments: vec![crate::truth_schema::TruthSegment {
                source: long_source(0, 8),
                modifications: vec![TruthModification::Snp {
                    event_id: "SNP_03".into(),
                    position: 3,
                    reference: "A".into(),
                    alternate: "T".into(),
                    provenance: provenance(),
                }],
            }],
            current_pieces: vec![
                TruthCurrentPiece::Source {
                    source: long_source(0, 6),
                    snps: vec![TruthCurrentSnp {
                        event_id: "SNP_03".into(),
                        position: 3,
                        reference: "A".into(),
                        alternate: "T".into(),
                        provenance: provenance(),
                    }],
                },
                TruthCurrentPiece::Source {
                    source: long_source(7, 8),
                    snps: vec![],
                },
            ],
            history: vec![
                TruthHistoryEvent::Snp {
                    event_id: "SNP_03".into(),
                    position: 3,
                    provenance: provenance(),
                },
                TruthHistoryEvent::Deletion {
                    event_id: "DEL_67".into(),
                    start: 6,
                    end: 7,
                    provenance: provenance(),
                },
            ],
        };
        let truth = TruthDocument {
            header: TruthHeader {
                document_kind: TRUTH_DOCUMENT_KIND.into(),
                schema_version: TRUTH_SCHEMA_VERSION.into(),
                execution_id: "RUN_01".into(),
                reference_id: "REF_main".into(),
                reference_content_sha256: reference_content_sha256(&reference),
                coordinate_convention: TRUTH_COORDINATE_CONVENTION.into(),
                copy_count: 1,
                genotype_count: 2,
                source_digests: vec![NamedDigest {
                    source: TruthSourceRef {
                        kind: TruthSourceKind::Reference,
                        source_id: "REF_main".into(),
                    },
                    sha256: "a".repeat(64),
                }],
                group_digests: vec![],
            },
            records: vec![
                TruthRecord::HaplotypeCopy(copy),
                TruthRecord::Genotype(GenotypeRecord {
                    case_id: "CASE_SNP".into(),
                    reference_id: "REF_main".into(),
                    individual_id: "IND".into(),
                    case_kind: GenotypeCaseKind::Snp,
                    locus: TruthLocus {
                        sequence_id: "chr1".into(),
                        start: 3,
                        end: 4,
                    },
                    alleles: vec![
                        TruthAllele::Reference,
                        TruthAllele::Snp { base: "T".into() },
                    ],
                    allele_dosages: vec![0, 1],
                    carriers: vec![AlleleCarrier {
                        allele_index: 1,
                        copy: identity("IND", 0),
                    }],
                }),
                TruthRecord::Genotype(GenotypeRecord {
                    case_id: "CASE_DEL".into(),
                    reference_id: "REF_main".into(),
                    individual_id: "IND".into(),
                    case_kind: GenotypeCaseKind::Indel,
                    locus: TruthLocus {
                        sequence_id: "chr1".into(),
                        start: 6,
                        end: 7,
                    },
                    alleles: vec![
                        TruthAllele::Reference,
                        TruthAllele::Deletion {
                            event_id: "DEL_67".into(),
                            start: 6,
                            end: 7,
                            final_haplotype_bases: "CAATAAA".into(),
                            source_event: TruthStructuralEvent::Deletion {
                                event_id: "DEL_67".into(),
                                start: 6,
                                end: 7,
                            },
                            piece_links: Vec::<TruthPieceLink>::new(),
                        },
                    ],
                    allele_dosages: vec![0, 1],
                    carriers: vec![AlleleCarrier {
                        allele_index: 1,
                        copy: identity("IND", 0),
                    }],
                }),
            ],
        };

        let normalized = crate::vcf_projection::project_normalized_alleles(&truth, &reference)
            .expect("projection");
        assert_eq!(normalized.len(), 1);
        assert_eq!(
            (
                normalized[0].position,
                normalized[0].reference.as_str(),
                normalized[0].alternate.as_str(),
                normalized[0].source_start,
                normalized[0].source_end
            ),
            (3, "AA", "T", 3, 7)
        );
        let records = group_normalized_alleles(&truth, &reference, &normalized).expect("grouping");
        assert_eq!(records[0].alternate_cases, [vec!["CASE_DEL", "CASE_SNP"]]);
        assert_eq!(records[0].gps_cases, ["CASE_DEL", "CASE_SNP"]);
    }

    #[test]
    fn rejects_noncontiguous_copy_phases_before_gt_construction() {
        for observed in [vec![0, 2], vec![1]] {
            let copies = observed
                .iter()
                .map(|phase| (identity("IND", *phase), *phase))
                .collect::<BTreeMap<_, _>>();
            let error = validate_contiguous_copy_phases(&copies)
                .expect_err("phases must be exactly zero through P minus one");
            assert!(matches!(
                error,
                VcfGroupingError::NonContiguousCopyPhases { observed: actual, .. } if actual == observed
            ));
        }
    }

    #[test]
    fn expands_to_a_common_ref_and_preserves_phase_and_dosage() {
        let copies = BTreeMap::from([(identity("IND", 0), 0), (identity("IND", 1), 1)]);
        let (record, _) = build_record(
            "chr1".into(),
            0,
            vec![
                &NormalizedAllele {
                    sequence_id: "chr1".into(),
                    position: 0,
                    reference: "C".into(),
                    alternate: "CA".into(),
                    source_start: 0,
                    source_end: 1,
                    copy: identity("IND", 0),
                },
                &NormalizedAllele {
                    sequence_id: "chr1".into(),
                    position: 0,
                    reference: "CA".into(),
                    alternate: "C".into(),
                    source_start: 0,
                    source_end: 2,
                    copy: identity("IND", 1),
                },
            ],
            &copies,
        )
        .expect("compatible alleles");
        assert_eq!(record.reference, "CA");
        assert_eq!(record.alternates, ["C", "CAA"]);
        assert_eq!(
            record.samples,
            vec![(
                "IND".into(),
                PhasedGenotype {
                    alleles: vec![2, 1],
                    allele_dosages: vec![0, 1, 1]
                }
            )]
        );
    }

    #[test]
    fn supports_haploid_diploid_and_decaploid_gt_dosages() {
        for ploidy in [1_u8, 2, 10] {
            let copies = (0..ploidy)
                .map(|index| (identity("IND", index), index))
                .collect::<BTreeMap<_, _>>();
            let (record, _) = build_record(
                "chr1".into(),
                1,
                vec![&NormalizedAllele {
                    sequence_id: "chr1".into(),
                    position: 1,
                    reference: "C".into(),
                    alternate: "T".into(),
                    source_start: 1,
                    source_end: 2,
                    copy: identity("IND", 0),
                }],
                &copies,
            )
            .expect("GT");
            let gt = &record.samples[0].1;
            assert_eq!(gt.alleles.len(), usize::from(ploidy));
            assert_eq!(gt.alleles[0], 1);
            assert!(gt.alleles[1..].iter().all(|allele| *allele == 0));
            assert_eq!(gt.allele_dosages, [ploidy - 1, 1]);
        }
    }

    #[test]
    fn rejects_a_duplicate_copy_at_one_locus() {
        let copy = identity("IND", 0);
        let copies = BTreeMap::from([(copy.clone(), 0)]);
        let error = build_record(
            "chr1".into(),
            1,
            vec![
                &NormalizedAllele {
                    sequence_id: "chr1".into(),
                    position: 1,
                    reference: "C".into(),
                    alternate: "T".into(),
                    source_start: 1,
                    source_end: 2,
                    copy: copy.clone(),
                },
                &NormalizedAllele {
                    sequence_id: "chr1".into(),
                    position: 1,
                    reference: "C".into(),
                    alternate: "G".into(),
                    source_start: 1,
                    source_end: 2,
                    copy: copy.clone(),
                },
            ],
            &copies,
        )
        .expect_err("a copy has only one final allele at a locus");
        assert!(matches!(
            error,
            VcfGroupingError::DuplicateCopyAtLocus { .. }
        ));
    }
}
