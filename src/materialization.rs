//! Reconstrucción determinista bajo demanda y consultas SNP desde el estado canónico.

use crate::coordinates::BasePosition;
use crate::haplotype::{HaplotypeCopy, HaplotypePhase, Individual};
use crate::identity::SequenceId;
use crate::reference::{Reference, ReferenceSequence};
use crate::segments::{CanonicalBase, SmallModification};

/// Error al reconstruir una copia desde sus segmentos y una referencia concreta.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum MaterializationError {
    CopyReferenceMismatch,
    SourceReferenceMismatch,
    SourceSequenceMismatch,
    SourceLengthMismatch {
        expected: u64,
        found: u64,
    },
    MissingReferenceSequence,
    IndexOverflow,
    SnpReferenceMismatch {
        position: u64,
        expected: CanonicalBase,
        found: u8,
    },
    ModificationOrderMismatch,
}

/// SNP cuyo alelo se consulta sobre las coordenadas de su fuente canónica.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SnpLocus {
    position: BasePosition,
    reference: CanonicalBase,
    alternate: CanonicalBase,
}

impl SnpLocus {
    pub fn new(
        position: BasePosition,
        reference: CanonicalBase,
        alternate: CanonicalBase,
    ) -> Result<Self, SnpLocusError> {
        if reference == alternate {
            return Err(SnpLocusError::AllelesEqual);
        }
        Ok(Self {
            position,
            reference,
            alternate,
        })
    }

    #[must_use]
    pub const fn position(&self) -> &BasePosition {
        &self.position
    }

    #[must_use]
    pub const fn reference(&self) -> CanonicalBase {
        self.reference
    }

    #[must_use]
    pub const fn alternate(&self) -> CanonicalBase {
        self.alternate
    }
}

/// Rechazo de un locus que no representa una sustitución SNP real.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SnpLocusError {
    AllelesEqual,
}

/// Alelo observado en una copia y su fase real.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PhasedAllele {
    phase: HaplotypePhase,
    allele: CanonicalBase,
}

impl PhasedAllele {
    #[must_use]
    pub const fn phase(self) -> HaplotypePhase {
        self.phase
    }

    #[must_use]
    pub const fn allele(self) -> CanonicalBase {
        self.allele
    }
}

/// Genotipo SNP derivado de las copias, sin consultar ninguna salida serializada.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PhasedSnpGenotype {
    alleles: Vec<PhasedAllele>,
    alternate_dosage: u8,
}

impl PhasedSnpGenotype {
    #[must_use]
    pub fn alleles(&self) -> &[PhasedAllele] {
        &self.alleles
    }

    #[must_use]
    pub const fn alternate_dosage(&self) -> u8 {
        self.alternate_dosage
    }
}

/// Error al derivar un genotipo desde copias canónicas.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GenotypeQueryError {
    IndividualReferenceMismatch,
    LocusReferenceMismatch,
    MissingSequence,
    LocusLengthMismatch { expected: u64, found: u64 },
    ReferenceBaseMismatch { expected: CanonicalBase, found: u8 },
    Materialization(MaterializationError),
    LocusAbsentFromCopy { phase: HaplotypePhase },
    AmbiguousLocusInCopy { phase: HaplotypePhase },
    LocusDeletedFromCopy { phase: HaplotypePhase },
    IncompatibleSnpInCopy { phase: HaplotypePhase },
}

/// Materializa una copia concreta, sin almacenar el resultado en el estado haplotípico.
pub fn materialize_copy(
    copy: &HaplotypeCopy,
    reference: &Reference,
) -> Result<Vec<u8>, MaterializationError> {
    if copy.id().reference_id() != reference.id() {
        return Err(MaterializationError::CopyReferenceMismatch);
    }

    let mut materialized = Vec::new();
    for segment in copy.segments().segments() {
        let source = segment.source();
        if source.reference_id() != reference.id() {
            return Err(MaterializationError::SourceReferenceMismatch);
        }
        if source.sequence_id() != copy.id().sequence_id() {
            return Err(MaterializationError::SourceSequenceMismatch);
        }
        let sequence = reference_sequence(reference, source.sequence_id())?;
        if source.sequence_length() != sequence.length() {
            return Err(MaterializationError::SourceLengthMismatch {
                expected: source.sequence_length(),
                found: sequence.length(),
            });
        }
        materialize_segment(&mut materialized, segment, sequence)?;
    }
    Ok(materialized)
}

/// Deriva los alelos faseados y la dosis ALT de un SNP desde las copias del individuo.
pub fn derive_snp_genotype(
    individual: &Individual,
    reference: &Reference,
    locus: &SnpLocus,
) -> Result<PhasedSnpGenotype, GenotypeQueryError> {
    if individual.reference_id() != reference.id() {
        return Err(GenotypeQueryError::IndividualReferenceMismatch);
    }
    if locus.position.reference_id() != reference.id() {
        return Err(GenotypeQueryError::LocusReferenceMismatch);
    }
    let sequence = reference_sequence(reference, locus.position.sequence_id())
        .map_err(GenotypeQueryError::Materialization)?;
    if locus.position.sequence_length() != sequence.length() {
        return Err(GenotypeQueryError::LocusLengthMismatch {
            expected: sequence.length(),
            found: locus.position.sequence_length(),
        });
    }
    let reference_base = byte_at(sequence.sequence().as_bytes(), locus.position.offset())
        .map_err(GenotypeQueryError::Materialization)?;
    if reference_base != locus.reference.as_byte() {
        return Err(GenotypeQueryError::ReferenceBaseMismatch {
            expected: locus.reference,
            found: reference_base,
        });
    }
    let copies = individual
        .copies_for(locus.position.sequence_id())
        .ok_or(GenotypeQueryError::MissingSequence)?;
    let mut alleles = Vec::with_capacity(copies.len());
    for copy in copies {
        let allele = allele_in_copy(copy, locus)?;
        alleles.push(PhasedAllele {
            phase: copy.phase(),
            allele,
        });
    }
    alleles.sort_by_key(|allele| allele.phase);
    let alternate_dosage = u8::try_from(
        alleles
            .iter()
            .filter(|allele| allele.allele == locus.alternate)
            .count(),
    )
    .map_err(|_| GenotypeQueryError::Materialization(MaterializationError::IndexOverflow))?;
    Ok(PhasedSnpGenotype {
        alleles,
        alternate_dosage,
    })
}

fn materialize_segment(
    destination: &mut Vec<u8>,
    segment: &crate::segments::Segment,
    sequence: &ReferenceSequence,
) -> Result<(), MaterializationError> {
    let source = segment.source();
    let bases = sequence.sequence().as_bytes();
    let mut cursor = source.start();
    for modification in segment.modifications() {
        match modification {
            SmallModification::Snp(snp) => {
                append_source(destination, bases, cursor, snp.position().offset())?;
                let found = byte_at(bases, snp.position().offset())?;
                if found != snp.reference().as_byte() {
                    return Err(MaterializationError::SnpReferenceMismatch {
                        position: snp.position().offset(),
                        expected: snp.reference(),
                        found,
                    });
                }
                destination.push(snp.alternate().as_byte());
                cursor = snp
                    .position()
                    .offset()
                    .checked_add(1)
                    .ok_or(MaterializationError::IndexOverflow)?;
            }
            SmallModification::Insertion(insertion) => {
                let point = insertion.anchor().point().offset();
                append_source(destination, bases, cursor, point)?;
                destination.extend_from_slice(insertion.bases());
                cursor = point;
            }
            SmallModification::Deletion(deletion) => {
                let interval = deletion.deletion().interval();
                append_source(destination, bases, cursor, interval.start())?;
                cursor = interval.end();
            }
        }
    }
    append_source(destination, bases, cursor, source.end())
}

fn allele_in_copy(
    copy: &HaplotypeCopy,
    locus: &SnpLocus,
) -> Result<CanonicalBase, GenotypeQueryError> {
    let mut matching_segment = None;
    for segment in copy.segments().segments() {
        let source = segment.source();
        if source.reference_id() != locus.position.reference_id()
            || source.sequence_id() != locus.position.sequence_id()
        {
            return Err(GenotypeQueryError::Materialization(
                MaterializationError::SourceSequenceMismatch,
            ));
        }
        if source.sequence_length() != locus.position.sequence_length() {
            return Err(GenotypeQueryError::Materialization(
                MaterializationError::SourceLengthMismatch {
                    expected: locus.position.sequence_length(),
                    found: source.sequence_length(),
                },
            ));
        }
        if locus.position.offset() >= source.start()
            && locus.position.offset() < source.end()
            && matching_segment.replace(segment).is_some()
        {
            return Err(GenotypeQueryError::AmbiguousLocusInCopy {
                phase: copy.phase(),
            });
        }
    }
    let segment = matching_segment.ok_or(GenotypeQueryError::LocusAbsentFromCopy {
        phase: copy.phase(),
    })?;
    for modification in segment.modifications() {
        match modification {
            SmallModification::Deletion(deletion)
                if locus.position.offset() >= deletion.deletion().interval().start()
                    && locus.position.offset() < deletion.deletion().interval().end() =>
            {
                return Err(GenotypeQueryError::LocusDeletedFromCopy {
                    phase: copy.phase(),
                });
            }
            SmallModification::Snp(snp) if snp.position() == locus.position() => {
                if snp.reference() != locus.reference || snp.alternate() != locus.alternate {
                    return Err(GenotypeQueryError::IncompatibleSnpInCopy {
                        phase: copy.phase(),
                    });
                }
                return Ok(locus.alternate);
            }
            _ => {}
        }
    }
    Ok(locus.reference)
}

fn reference_sequence<'a>(
    reference: &'a Reference,
    sequence_id: &SequenceId,
) -> Result<&'a ReferenceSequence, MaterializationError> {
    reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id() == sequence_id)
        .ok_or(MaterializationError::MissingReferenceSequence)
}

fn append_source(
    destination: &mut Vec<u8>,
    source: &[u8],
    start: u64,
    end: u64,
) -> Result<(), MaterializationError> {
    if start > end {
        return Err(MaterializationError::ModificationOrderMismatch);
    }
    let start = usize::try_from(start).map_err(|_| MaterializationError::IndexOverflow)?;
    let end = usize::try_from(end).map_err(|_| MaterializationError::IndexOverflow)?;
    let bases = source
        .get(start..end)
        .ok_or(MaterializationError::IndexOverflow)?;
    destination.extend_from_slice(bases);
    Ok(())
}

fn byte_at(source: &[u8], offset: u64) -> Result<u8, MaterializationError> {
    let offset = usize::try_from(offset).map_err(|_| MaterializationError::IndexOverflow)?;
    source
        .get(offset)
        .copied()
        .ok_or(MaterializationError::IndexOverflow)
}
