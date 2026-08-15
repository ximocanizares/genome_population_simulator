//! Coordenadas genómicas tipadas y conversiones explícitas de frontera.

use crate::identity::{ReferenceId, SequenceId};

/// Errores de construcción y conversión de coordenadas.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CoordinateError {
    ZeroBasedOneCoordinate,
    EmptyInterval,
    BaseOutOfBounds {
        position: u64,
        sequence_length: u64,
    },
    InterbaseOutOfBounds {
        position: u64,
        sequence_length: u64,
    },
    IntervalOutOfBounds {
        start: u64,
        end: u64,
        sequence_length: u64,
    },
    ReferenceMismatch,
    SequenceMismatch,
    SequenceLengthMismatch {
        start_length: u64,
        end_length: u64,
    },
    ArithmeticOverflow,
}

/// Posición de una base en convenio interno base 0.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BasePosition {
    reference_id: ReferenceId,
    sequence_id: SequenceId,
    offset: u64,
    sequence_length: u64,
}

impl BasePosition {
    pub fn new(
        reference_id: ReferenceId,
        sequence_id: SequenceId,
        offset: u64,
        sequence_length: u64,
    ) -> Result<Self, CoordinateError> {
        if offset >= sequence_length {
            return Err(CoordinateError::BaseOutOfBounds {
                position: offset,
                sequence_length,
            });
        }
        Ok(Self {
            reference_id,
            sequence_id,
            offset,
            sequence_length,
        })
    }

    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub const fn offset(&self) -> u64 {
        self.offset
    }

    pub fn to_base_1(&self) -> Result<u64, CoordinateError> {
        self.offset
            .checked_add(1)
            .ok_or(CoordinateError::ArithmeticOverflow)
    }
}

/// Punto entre bases en convenio interno base 0.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InterbasePosition {
    reference_id: ReferenceId,
    sequence_id: SequenceId,
    offset: u64,
    sequence_length: u64,
}

impl InterbasePosition {
    pub fn new(
        reference_id: ReferenceId,
        sequence_id: SequenceId,
        offset: u64,
        sequence_length: u64,
    ) -> Result<Self, CoordinateError> {
        if offset > sequence_length {
            return Err(CoordinateError::InterbaseOutOfBounds {
                position: offset,
                sequence_length,
            });
        }
        Ok(Self {
            reference_id,
            sequence_id,
            offset,
            sequence_length,
        })
    }

    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub const fn offset(&self) -> u64 {
        self.offset
    }
}

/// DTO de frontera para un intervalo base 1 cerrado `[start, end]`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ClosedInterval1 {
    start: u64,
    end: u64,
}

impl ClosedInterval1 {
    pub const fn new(start: u64, end: u64) -> Result<Self, CoordinateError> {
        if start == 0 || end == 0 {
            return Err(CoordinateError::ZeroBasedOneCoordinate);
        }
        if start > end {
            return Err(CoordinateError::EmptyInterval);
        }
        Ok(Self { start, end })
    }

    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn end(self) -> u64 {
        self.end
    }
}

/// Intervalo genómico interno base 0 semiabierto `[start, end)`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GenomicInterval {
    reference_id: ReferenceId,
    sequence_id: SequenceId,
    start: u64,
    end: u64,
}

impl GenomicInterval {
    pub fn new(
        reference_id: ReferenceId,
        sequence_id: SequenceId,
        start: u64,
        end: u64,
        sequence_length: u64,
    ) -> Result<Self, CoordinateError> {
        if start >= end {
            return Err(CoordinateError::EmptyInterval);
        }
        if end > sequence_length {
            return Err(CoordinateError::IntervalOutOfBounds {
                start,
                end,
                sequence_length,
            });
        }
        Ok(Self {
            reference_id,
            sequence_id,
            start,
            end,
        })
    }

    pub fn from_positions(
        start: BasePosition,
        end: InterbasePosition,
    ) -> Result<Self, CoordinateError> {
        if start.reference_id != end.reference_id {
            return Err(CoordinateError::ReferenceMismatch);
        }
        if start.sequence_id != end.sequence_id {
            return Err(CoordinateError::SequenceMismatch);
        }
        if start.sequence_length != end.sequence_length {
            return Err(CoordinateError::SequenceLengthMismatch {
                start_length: start.sequence_length,
                end_length: end.sequence_length,
            });
        }
        if start.offset >= end.offset {
            return Err(CoordinateError::EmptyInterval);
        }
        Ok(Self {
            reference_id: start.reference_id,
            sequence_id: start.sequence_id,
            start: start.offset,
            end: end.offset,
        })
    }

    pub fn from_closed_1(
        reference_id: ReferenceId,
        sequence_id: SequenceId,
        interval: ClosedInterval1,
        sequence_length: u64,
    ) -> Result<Self, CoordinateError> {
        let start = interval
            .start
            .checked_sub(1)
            .ok_or(CoordinateError::ArithmeticOverflow)?;
        Self::new(
            reference_id,
            sequence_id,
            start,
            interval.end,
            sequence_length,
        )
    }

    pub fn to_closed_1(&self) -> Result<ClosedInterval1, CoordinateError> {
        let start = self
            .start
            .checked_add(1)
            .ok_or(CoordinateError::ArithmeticOverflow)?;
        ClosedInterval1::new(start, self.end)
    }

    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn end(&self) -> u64 {
        self.end
    }

    #[must_use]
    pub const fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Los intervalos de dominio nunca están vacíos por construcción.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reference_id(value: &str) -> ReferenceId {
        ReferenceId::parse(value, "coordinate-test").expect("valid reference ID")
    }

    fn sequence_id(value: &str) -> SequenceId {
        SequenceId::parse(value, "coordinate-test").expect("valid sequence ID")
    }

    #[test]
    fn converts_first_base_and_complete_sequence_exactly() {
        let first = GenomicInterval::from_closed_1(
            reference_id("ref"),
            sequence_id("chr1"),
            ClosedInterval1::new(1, 1).expect("closed interval"),
            10,
        )
        .expect("internal interval");
        assert_eq!((first.start(), first.end()), (0, 1));
        assert_eq!(
            first.to_closed_1().expect("external interval"),
            ClosedInterval1::new(1, 1).expect("closed interval")
        );

        let complete = GenomicInterval::from_closed_1(
            reference_id("ref"),
            sequence_id("chr1"),
            ClosedInterval1::new(1, 10).expect("closed interval"),
            10,
        )
        .expect("internal interval");
        assert_eq!(
            (complete.start(), complete.end(), complete.len()),
            (0, 10, 10)
        );
        assert_eq!(complete.to_closed_1().expect("round trip").end(), 10);
    }

    #[test]
    fn accepts_both_interbase_boundaries_but_not_an_outside_point() {
        let last_base =
            BasePosition::new(reference_id("ref"), sequence_id("chr1"), 9, 10).expect("last base");
        assert_eq!(last_base.to_base_1().expect("base-1 position"), 10);
        assert_eq!(
            BasePosition::new(reference_id("ref"), sequence_id("chr1"), 10, 10),
            Err(CoordinateError::BaseOutOfBounds {
                position: 10,
                sequence_length: 10
            })
        );
        assert!(InterbasePosition::new(reference_id("ref"), sequence_id("chr1"), 0, 10).is_ok());
        assert!(InterbasePosition::new(reference_id("ref"), sequence_id("chr1"), 10, 10).is_ok());
        assert_eq!(
            InterbasePosition::new(reference_id("ref"), sequence_id("chr1"), 11, 10),
            Err(CoordinateError::InterbaseOutOfBounds {
                position: 11,
                sequence_length: 10
            })
        );
    }

    #[test]
    fn rejects_empty_out_of_bounds_and_zero_based_one_intervals() {
        assert_eq!(
            GenomicInterval::new(reference_id("ref"), sequence_id("chr1"), 4, 4, 10),
            Err(CoordinateError::EmptyInterval)
        );
        assert!(matches!(
            GenomicInterval::new(reference_id("ref"), sequence_id("chr1"), 9, 11, 10),
            Err(CoordinateError::IntervalOutOfBounds { .. })
        ));
        assert_eq!(
            ClosedInterval1::new(0, 1),
            Err(CoordinateError::ZeroBasedOneCoordinate)
        );
    }

    #[test]
    fn positions_from_different_spaces_cannot_form_an_interval() {
        let start =
            BasePosition::new(reference_id("ref-a"), sequence_id("chr1"), 0, 10).expect("base");
        let end = InterbasePosition::new(reference_id("ref-b"), sequence_id("chr1"), 1, 10)
            .expect("interbase");
        assert_eq!(
            GenomicInterval::from_positions(start, end),
            Err(CoordinateError::ReferenceMismatch)
        );

        let start =
            BasePosition::new(reference_id("ref"), sequence_id("chr1"), 0, 10).expect("base");
        let end = InterbasePosition::new(reference_id("ref"), sequence_id("chr2"), 1, 10)
            .expect("interbase");
        assert_eq!(
            GenomicInterval::from_positions(start, end),
            Err(CoordinateError::SequenceMismatch)
        );

        let start =
            BasePosition::new(reference_id("ref"), sequence_id("chr1"), 0, 10).expect("base");
        let end = InterbasePosition::new(reference_id("ref"), sequence_id("chr1"), 1, 11)
            .expect("interbase");
        assert_eq!(
            GenomicInterval::from_positions(start, end),
            Err(CoordinateError::SequenceLengthMismatch {
                start_length: 10,
                end_length: 11
            })
        );
    }
}
