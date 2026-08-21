//! Segmentos haplotípicos y modificaciones pequeñas sin materializar referencias.
//!
//! Una [`SegmentChain`] describe la secuencia lógica de una copia haplotípica.
//! Sus segmentos ordinarios conservan solamente coordenadas de una fuente; las
//! únicas bases que posee el estado son las de una inserción nueva.

use std::sync::Arc;

use crate::coordinates::{BasePosition, GenomicInterval, InsertionAnchor, InternalDeletion};
use crate::identity::EventId;

/// Base concreta admitida como alelo de un SNP.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CanonicalBase {
    A,
    C,
    G,
    T,
}

impl CanonicalBase {
    pub fn from_byte(value: u8) -> Result<Self, SegmentError> {
        match value.to_ascii_uppercase() {
            b'A' => Ok(Self::A),
            b'C' => Ok(Self::C),
            b'G' => Ok(Self::G),
            b'T' => Ok(Self::T),
            _ => Err(SegmentError::NonCanonicalSnpBase { value }),
        }
    }

    #[must_use]
    pub const fn as_byte(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
        }
    }
}

/// Procedencia explícita de una modificación pequeña.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModificationProvenance {
    event_id: EventId,
}

impl ModificationProvenance {
    #[must_use]
    pub const fn new(event_id: EventId) -> Self {
        Self { event_id }
    }

    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        &self.event_id
    }
}

/// Sustitución de una base de fuente por otra base canónica.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SnpModification {
    position: BasePosition,
    reference: CanonicalBase,
    alternate: CanonicalBase,
    provenance: ModificationProvenance,
}

impl SnpModification {
    pub fn new(
        position: BasePosition,
        reference: CanonicalBase,
        alternate: CanonicalBase,
        provenance: ModificationProvenance,
    ) -> Result<Self, SegmentError> {
        if reference == alternate {
            return Err(SegmentError::SnpDoesNotChangeBase);
        }
        Ok(Self {
            position,
            reference,
            alternate,
            provenance,
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

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }
}

/// Inserción con bases propias, pues no procede de un intervalo de fuente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InsertionModification {
    anchor: InsertionAnchor,
    bases: Arc<[u8]>,
    provenance: ModificationProvenance,
}

impl InsertionModification {
    pub fn new(
        anchor: InsertionAnchor,
        bases: impl Into<Arc<[u8]>>,
        provenance: ModificationProvenance,
    ) -> Result<Self, SegmentError> {
        let bases = bases.into();
        if bases.is_empty() {
            return Err(SegmentError::EmptyInsertion);
        }
        for &base in bases.iter() {
            if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                return Err(SegmentError::InvalidInsertedBase { base });
            }
        }
        Ok(Self {
            anchor,
            bases,
            provenance,
        })
    }

    #[must_use]
    pub const fn anchor(&self) -> &InsertionAnchor {
        &self.anchor
    }

    /// Las bases devueltas son el único contenido secuencial poseído por el estado.
    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }
}

/// Deleción pequeña de un intervalo interno de fuente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeletionModification {
    deletion: InternalDeletion,
    provenance: ModificationProvenance,
}

impl DeletionModification {
    #[must_use]
    pub const fn new(deletion: InternalDeletion, provenance: ModificationProvenance) -> Self {
        Self {
            deletion,
            provenance,
        }
    }

    #[must_use]
    pub const fn deletion(&self) -> &InternalDeletion {
        &self.deletion
    }

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }
}

/// Una modificación pequeña tipada sobre un segmento de fuente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SmallModification {
    Snp(SnpModification),
    Insertion(InsertionModification),
    Deletion(DeletionModification),
}

/// Segmento ordinario que puede reconstruirse desde su intervalo de fuente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Segment {
    source: GenomicInterval,
    modifications: Vec<SmallModification>,
}

impl Segment {
    pub fn new(
        source: GenomicInterval,
        modifications: Vec<SmallModification>,
    ) -> Result<Self, SegmentError> {
        let mut previous: Option<ModificationExtent> = None;
        for modification in &modifications {
            let extent = ModificationExtent::from_modification(modification);
            if !modification_matches_source(modification, &source) {
                return Err(SegmentError::ModificationSourceMismatch);
            }
            validate_extent_in_source(&extent, &source)?;
            if let Some(previous) = previous {
                if previous.sort_key() >= extent.sort_key() {
                    return Err(SegmentError::ModificationsNotStrictlyOrdered);
                }
                if previous.overlaps(extent) {
                    return Err(SegmentError::OverlappingModifications);
                }
            }
            previous = Some(extent);
        }
        Ok(Self {
            source,
            modifications,
        })
    }

    #[must_use]
    pub const fn source(&self) -> &GenomicInterval {
        &self.source
    }

    #[must_use]
    pub fn modifications(&self) -> &[SmallModification] {
        &self.modifications
    }
}

/// Cadena ordenada, persistente y compartible de segmentos de una copia.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SegmentChain {
    segments: Arc<[Segment]>,
}

impl SegmentChain {
    /// Crea la cadena intacta de una copia desde un intervalo ya validado.
    #[must_use]
    pub fn intact(source: GenomicInterval) -> Self {
        Self {
            segments: Arc::from(vec![Segment {
                source,
                modifications: Vec::new(),
            }]),
        }
    }

    pub fn new(segments: Vec<Segment>) -> Result<Self, SegmentError> {
        if segments.is_empty() {
            return Err(SegmentError::EmptyChain);
        }
        let first = segments[0].source();
        if segments.iter().skip(1).any(|segment| {
            let source = segment.source();
            source.reference_id() != first.reference_id()
                || source.sequence_id() != first.sequence_id()
                || source.sequence_length() != first.sequence_length()
        }) {
            return Err(SegmentError::MixedSourceSpace);
        }
        Ok(Self {
            segments: Arc::from(segments),
        })
    }

    #[must_use]
    pub fn segments(&self) -> &[Segment] {
        &self.segments
    }

    /// Indica si dos cadenas comparten realmente la misma asignación persistente.
    #[must_use]
    pub fn is_shared_with(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.segments, &other.segments)
    }
}

#[derive(Clone, Copy)]
struct ModificationExtent {
    start: u64,
    end: u64,
    is_insertion: bool,
    is_deletion: bool,
}

impl ModificationExtent {
    fn from_modification(modification: &SmallModification) -> Self {
        match modification {
            SmallModification::Snp(snp) => Self {
                start: snp.position.offset(),
                end: snp.position.offset() + 1,
                is_insertion: false,
                is_deletion: false,
            },
            SmallModification::Insertion(insertion) => Self {
                start: insertion.anchor.point().offset(),
                end: insertion.anchor.point().offset(),
                is_insertion: true,
                is_deletion: false,
            },
            SmallModification::Deletion(deletion) => Self {
                start: deletion.deletion.interval().start(),
                end: deletion.deletion.interval().end(),
                is_insertion: false,
                is_deletion: true,
            },
        }
    }

    const fn sort_key(self) -> (u64, u64, bool) {
        (self.start, self.end, self.is_insertion)
    }

    const fn overlaps(self, other: Self) -> bool {
        if self.is_insertion && other.is_insertion {
            self.start == other.start
        } else if self.is_insertion && other.is_deletion {
            other.start < self.start && self.start <= other.end
        } else if other.is_insertion && self.is_deletion {
            self.start < other.start && other.start <= self.end
        } else if self.is_insertion || other.is_insertion {
            false
        } else {
            self.start < other.end && other.start < self.end
        }
    }
}

fn validate_extent_in_source(
    extent: &ModificationExtent,
    source: &GenomicInterval,
) -> Result<(), SegmentError> {
    if extent.start < source.start() || extent.end > source.end() {
        return Err(SegmentError::ModificationOutsideSource);
    }
    Ok(())
}

fn modification_matches_source(modification: &SmallModification, source: &GenomicInterval) -> bool {
    let (reference_id, sequence_id, sequence_length) = match modification {
        SmallModification::Snp(value) => (
            value.position().reference_id(),
            value.position().sequence_id(),
            value.position().sequence_length(),
        ),
        SmallModification::Insertion(value) => (
            value.anchor().point().reference_id(),
            value.anchor().point().sequence_id(),
            value.anchor().point().sequence_length(),
        ),
        SmallModification::Deletion(value) => (
            value.deletion().interval().reference_id(),
            value.deletion().interval().sequence_id(),
            value.deletion().interval().sequence_length(),
        ),
    };
    reference_id == source.reference_id()
        && sequence_id == source.sequence_id()
        && sequence_length == source.sequence_length()
}

/// Invariantes que evitan una cadena o una composición ambigua.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SegmentError {
    NonCanonicalSnpBase { value: u8 },
    SnpDoesNotChangeBase,
    EmptyInsertion,
    InvalidInsertedBase { base: u8 },
    EmptyChain,
    MixedSourceSpace,
    ModificationSourceMismatch,
    ModificationOutsideSource,
    ModificationsNotStrictlyOrdered,
    OverlappingModifications,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::{
        GenomicInterval, InsertionAnchor, InterbasePosition, InternalDeletion,
    };
    use crate::identity::{EventId, ReferenceId, SequenceId};

    fn interval(start: u64, end: u64) -> GenomicInterval {
        GenomicInterval::new(
            ReferenceId::parse("ref", "test").expect("reference"),
            SequenceId::parse("chr1", "test").expect("sequence"),
            start,
            end,
            10,
        )
        .expect("interval")
    }

    fn provenance() -> ModificationProvenance {
        ModificationProvenance::new(EventId::parse("event-1", "test").expect("event"))
    }

    #[test]
    fn retains_an_intact_source_segment_without_reference_bases() {
        let segment = Segment::new(interval(0, 10), vec![]).expect("intact segment");
        assert_eq!(segment.source().len(), 10);
        assert!(segment.modifications().is_empty());
    }

    #[test]
    fn keeps_each_small_edit_typed_and_its_provenance() {
        let snp = SmallModification::Snp(
            SnpModification::new(
                BasePosition::new(
                    ReferenceId::parse("ref", "test").expect("reference"),
                    SequenceId::parse("chr1", "test").expect("sequence"),
                    2,
                    10,
                )
                .expect("position"),
                CanonicalBase::G,
                CanonicalBase::T,
                provenance(),
            )
            .expect("SNP"),
        );
        let insertion = SmallModification::Insertion(
            InsertionModification::new(
                InsertionAnchor::new(
                    InterbasePosition::new(
                        ReferenceId::parse("ref", "test").expect("reference"),
                        SequenceId::parse("chr1", "test").expect("sequence"),
                        4,
                        10,
                    )
                    .expect("point"),
                )
                .expect("anchor"),
                Arc::<[u8]>::from(&b"NN"[..]),
                provenance(),
            )
            .expect("insertion"),
        );
        let deletion = SmallModification::Deletion(DeletionModification::new(
            InternalDeletion::new(interval(6, 8)).expect("deletion"),
            provenance(),
        ));
        let segment =
            Segment::new(interval(0, 10), vec![snp, insertion, deletion]).expect("ordered edits");
        assert_eq!(segment.modifications().len(), 3);
        assert!(matches!(
            segment.modifications()[0],
            SmallModification::Snp(_)
        ));
        assert!(matches!(
            segment.modifications()[1],
            SmallModification::Insertion(_)
        ));
        assert!(matches!(
            segment.modifications()[2],
            SmallModification::Deletion(_)
        ));
        for modification in segment.modifications() {
            let observed = match modification {
                SmallModification::Snp(value) => value.provenance(),
                SmallModification::Insertion(value) => value.provenance(),
                SmallModification::Deletion(value) => value.provenance(),
            };
            assert_eq!(observed.event_id().as_str(), "event-1");
        }
    }

    #[test]
    fn cloned_chains_share_segments_instead_of_reference_or_segment_copies() {
        let chain = SegmentChain::new(vec![
            Segment::new(interval(0, 10), vec![]).expect("segment"),
        ])
        .expect("chain");
        let shared = chain.clone();
        assert!(chain.is_shared_with(&shared));
        assert_eq!(chain.segments()[0].source().reference_id().as_str(), "ref");
    }

    #[test]
    fn preserves_the_declared_order_of_a_short_known_composition() {
        let chain = SegmentChain::new(vec![
            Segment::new(interval(0, 3), vec![]).expect("first segment"),
            Segment::new(interval(7, 10), vec![]).expect("second segment"),
        ])
        .expect("composition");
        let observed: Vec<_> = chain
            .segments()
            .iter()
            .map(|segment| (segment.source().start(), segment.source().end()))
            .collect();
        assert_eq!(observed, vec![(0, 3), (7, 10)]);
    }

    #[test]
    fn rejects_ambiguous_or_overlapping_modifications() {
        let first = SmallModification::Snp(
            SnpModification::new(
                BasePosition::new(
                    ReferenceId::parse("ref", "test").expect("reference"),
                    SequenceId::parse("chr1", "test").expect("sequence"),
                    2,
                    10,
                )
                .expect("position"),
                CanonicalBase::A,
                CanonicalBase::C,
                provenance(),
            )
            .expect("SNP"),
        );
        let repeated = first.clone();
        assert_eq!(
            Segment::new(interval(0, 10), vec![first, repeated]),
            Err(SegmentError::ModificationsNotStrictlyOrdered)
        );

        let outside = SmallModification::Snp(
            SnpModification::new(
                BasePosition::new(
                    ReferenceId::parse("ref", "test").expect("reference"),
                    SequenceId::parse("chr1", "test").expect("sequence"),
                    7,
                    10,
                )
                .expect("position"),
                CanonicalBase::A,
                CanonicalBase::C,
                provenance(),
            )
            .expect("SNP"),
        );
        assert_eq!(
            Segment::new(interval(0, 4), vec![outside]),
            Err(SegmentError::ModificationOutsideSource)
        );

        let deletion = SmallModification::Deletion(DeletionModification::new(
            InternalDeletion::new(interval(2, 5)).expect("deletion"),
            provenance(),
        ));
        let inside_deleted_bases = SmallModification::Insertion(
            InsertionModification::new(
                InsertionAnchor::new(
                    InterbasePosition::new(
                        ReferenceId::parse("ref", "test").expect("reference"),
                        SequenceId::parse("chr1", "test").expect("sequence"),
                        4,
                        10,
                    )
                    .expect("point"),
                )
                .expect("anchor"),
                Arc::<[u8]>::from(&b"A"[..]),
                provenance(),
            )
            .expect("insertion"),
        );
        assert_eq!(
            Segment::new(interval(0, 10), vec![deletion, inside_deleted_bases]),
            Err(SegmentError::OverlappingModifications)
        );

        let other_space = SmallModification::Snp(
            SnpModification::new(
                BasePosition::new(
                    ReferenceId::parse("other-ref", "test").expect("reference"),
                    SequenceId::parse("chr1", "test").expect("sequence"),
                    2,
                    10,
                )
                .expect("position"),
                CanonicalBase::A,
                CanonicalBase::C,
                provenance(),
            )
            .expect("SNP"),
        );
        assert_eq!(
            Segment::new(interval(0, 10), vec![other_space]),
            Err(SegmentError::ModificationSourceMismatch)
        );
    }
}
