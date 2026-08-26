//! Estado mínimo de individuos y copias haplotípicas, sin materializar bases.

use crate::coordinates::GenomicInterval;
use crate::identity::{IndividualId, ReferenceId, SequenceId};
use crate::reference::Reference;
use crate::segments::SegmentChain;

/// Ploidía uniforme de toda una ejecución, limitada por el contrato F2.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Ploidy(u8);

impl Ploidy {
    pub const MIN: u8 = 1;
    pub const MAX: u8 = 10;
    pub const DIPLOID: Self = Self(2);

    pub fn new(value: u8) -> Result<Self, PloidyError> {
        if !(Self::MIN..=Self::MAX).contains(&value) {
            return Err(PloidyError::OutOfRange { value });
        }
        Ok(Self(value))
    }

    #[must_use]
    pub const fn copy_count(self) -> u8 {
        self.0
    }
}

/// Rechazo local de una ploidía que no pertenece al dominio F2.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PloidyError {
    OutOfRange { value: u8 },
}

/// Posición faseada de una copia dentro de una secuencia del individuo.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct HaplotypePhase(u8);

impl HaplotypePhase {
    pub fn from_copy_index(copy_index: u8) -> Result<Self, HaplotypeStateError> {
        if copy_index >= Ploidy::MAX {
            return Err(HaplotypeStateError::PhaseOutOfRange { copy_index });
        }
        Ok(Self(copy_index))
    }

    #[must_use]
    pub const fn copy_index(self) -> u8 {
        self.0
    }
}

/// Identidad de una copia haplotípica, aun cuando varias copias tengan igual contenido futuro.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct HaplotypeCopyId {
    reference_id: ReferenceId,
    individual_id: IndividualId,
    sequence_id: SequenceId,
    copy_index: u8,
}

impl HaplotypeCopyId {
    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }

    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub const fn copy_index(&self) -> u8 {
        self.copy_index
    }
}

/// Copia haplotípica con una cadena compartible de segmentos no materializados.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HaplotypeCopy {
    id: HaplotypeCopyId,
    phase: HaplotypePhase,
    segments: SegmentChain,
}

impl HaplotypeCopy {
    /// Devuelve una nueva versión lógica de la copia con otra cadena ya validada.
    ///
    /// La sustitución conserva identidad y fase, y solo admite segmentos de la
    /// misma referencia, secuencia y longitud de fuente que la copia.
    pub fn with_segments(&self, segments: SegmentChain) -> Result<Self, HaplotypeStateError> {
        if segments.segments().iter().any(|segment| {
            segment.source().reference_id() != self.id.reference_id()
                || segment.source().sequence_id() != self.id.sequence_id()
        }) {
            return Err(HaplotypeStateError::SegmentSourceMismatch);
        }
        let expected_length = self.segments.segments()[0].source().sequence_length();
        if segments
            .segments()
            .iter()
            .any(|segment| segment.source().sequence_length() != expected_length)
        {
            return Err(HaplotypeStateError::SegmentSourceLengthMismatch {
                expected: expected_length,
                found: segments.segments()[0].source().sequence_length(),
            });
        }
        Ok(Self {
            id: self.id.clone(),
            phase: self.phase,
            segments,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &HaplotypeCopyId {
        &self.id
    }

    #[must_use]
    pub const fn phase(&self) -> HaplotypePhase {
        self.phase
    }

    #[must_use]
    pub const fn segments(&self) -> &SegmentChain {
        &self.segments
    }
}

/// Incompatibilidad entre una copia y los segmentos que pretende contener.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HaplotypeStateError {
    SegmentSourceMismatch,
    SegmentSourceLengthMismatch { expected: u64, found: u64 },
    MissingSequence,
    MissingCopyPhase,
    PhaseOutOfRange { copy_index: u8 },
}

/// Copias de una sola secuencia de referencia dentro de un individuo.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceGenotype {
    sequence_id: SequenceId,
    copies: Vec<HaplotypeCopy>,
}

impl SequenceGenotype {
    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub fn copies(&self) -> &[HaplotypeCopy] {
        &self.copies
    }
}

/// Estado genotípico mínimo asociado a una única referencia y ploidía uniforme.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Individual {
    id: IndividualId,
    reference_id: ReferenceId,
    ploidy: Ploidy,
    sequences: Vec<SequenceGenotype>,
}

impl Individual {
    /// Crea una copia lógica por índice y secuencia, sin clonar bases de la referencia.
    #[must_use]
    pub fn from_reference(id: IndividualId, reference: &Reference, ploidy: Ploidy) -> Self {
        let sequences = reference
            .sequences()
            .iter()
            .map(|sequence| {
                let sequence_id = sequence.id().clone();
                let segments = SegmentChain::intact(
                    GenomicInterval::new(
                        reference.id().clone(),
                        sequence_id.clone(),
                        0,
                        sequence.sequence().len(),
                        sequence.sequence().len(),
                    )
                    // Una `ReferenceSequence` válida nunca está vacía, de modo que
                    // el intervalo completo es una invariante del dominio.
                    .expect("a validated reference sequence has a complete interval"),
                );
                let copies = (0..ploidy.copy_count())
                    .map(|copy_index| HaplotypeCopy {
                        id: HaplotypeCopyId {
                            reference_id: reference.id().clone(),
                            individual_id: id.clone(),
                            sequence_id: sequence_id.clone(),
                            copy_index,
                        },
                        phase: HaplotypePhase(copy_index),
                        segments: segments.clone(),
                    })
                    .collect();
                SequenceGenotype {
                    sequence_id,
                    copies,
                }
            })
            .collect();
        Self {
            id,
            reference_id: reference.id().clone(),
            ploidy,
            sequences,
        }
    }

    #[must_use]
    pub const fn id(&self) -> &IndividualId {
        &self.id
    }

    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub const fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    #[must_use]
    pub fn sequences(&self) -> &[SequenceGenotype] {
        &self.sequences
    }

    #[must_use]
    pub fn copies_for(&self, sequence_id: &SequenceId) -> Option<&[HaplotypeCopy]> {
        self.sequences
            .iter()
            .find(|sequence| sequence.sequence_id() == sequence_id)
            .map(SequenceGenotype::copies)
    }

    /// Devuelve un nuevo individuo con una sola copia sustituida de forma persistente.
    pub fn with_copy_segments(
        &self,
        sequence_id: &SequenceId,
        phase: HaplotypePhase,
        segments: SegmentChain,
    ) -> Result<Self, HaplotypeStateError> {
        let mut updated = self.clone();
        updated.replace_copy_segments(sequence_id, phase, segments)?;
        Ok(updated)
    }

    /// Sustituye internamente una sola copia ya localizada. Esta operación no
    /// reconstruye el individuo: conserva todas las secuencias y copias no
    /// objetivo tal como están compartidas en la versión poseída.
    pub(crate) fn replace_copy_segments(
        &mut self,
        sequence_id: &SequenceId,
        phase: HaplotypePhase,
        segments: SegmentChain,
    ) -> Result<(), HaplotypeStateError> {
        let sequence = self
            .sequences
            .iter_mut()
            .find(|sequence| sequence.sequence_id() == sequence_id)
            .ok_or(HaplotypeStateError::MissingSequence)?;
        let copy = sequence
            .copies
            .iter_mut()
            .find(|copy| copy.phase() == phase)
            .ok_or(HaplotypeStateError::MissingCopyPhase)?;
        *copy = copy.with_segments(segments)?;
        Ok(())
    }
}
