//! Planificación reproducible de referencias reducidas.

use std::collections::{BTreeMap, BTreeSet};

use crate::coordinates::{CoordinateError, GenomicInterval};
use crate::identity::{ReferenceId, SequenceId};
use crate::randomness::LocalRng;
use crate::reference::SequenceClass;

pub const REDUCED_REFERENCE_SELECTION_ALGORITHM: &str = "reduced_reference_selection_v1";

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReductionError {
    EmptyEligibleSpace,
    ZeroLengthSequence(String),
    DuplicateSequenceId(String),
    ZeroRegionSize,
    ZeroTarget,
    TargetNotSmallerThanEligible { target: u64, eligible: u64 },
    InfeasibleLayout,
    ArithmeticOverflow,
    TooManySegments,
    AllocationFailed,
    Coordinate(CoordinateError),
    SelectionDoesNotMatchCatalog,
}

impl From<CoordinateError> for ReductionError {
    fn from(value: CoordinateError) -> Self {
        Self::Coordinate(value)
    }
}

/// Secuencia elegible en el orden del FASTA autoritativo.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EligibleSequence {
    id: SequenceId,
    class: SequenceClass,
    length: u64,
}

impl EligibleSequence {
    pub fn new(id: SequenceId, class: SequenceClass, length: u64) -> Result<Self, ReductionError> {
        if length == 0 {
            return Err(ReductionError::ZeroLengthSequence(id.as_str().to_owned()));
        }
        Ok(Self { id, class, length })
    }

    #[must_use]
    pub const fn id(&self) -> &SequenceId {
        &self.id
    }

    #[must_use]
    pub const fn class(&self) -> SequenceClass {
        self.class
    }

    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReductionSegmentClass {
    Full,
    Residual,
}

/// Segmento en el orden en que fue sorteado.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SelectedSegment {
    selection_order: u64,
    class: ReductionSegmentClass,
    sequence_class: SequenceClass,
    source_interval: GenomicInterval,
}

impl SelectedSegment {
    #[must_use]
    pub const fn selection_order(&self) -> u64 {
        self.selection_order
    }

    #[must_use]
    pub const fn class(&self) -> ReductionSegmentClass {
        self.class
    }

    #[must_use]
    pub const fn sequence_class(&self) -> SequenceClass {
        self.sequence_class
    }

    #[must_use]
    pub const fn source_interval(&self) -> &GenomicInterval {
        &self.source_interval
    }
}

/// Resultado completo e inmutable del sorteo.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SelectionPlan {
    source_reference_id: ReferenceId,
    target_bases: u64,
    region_size_bases: u64,
    segments: Vec<SelectedSegment>,
}

impl SelectionPlan {
    #[must_use]
    pub const fn source_reference_id(&self) -> &ReferenceId {
        &self.source_reference_id
    }

    #[must_use]
    pub const fn target_bases(&self) -> u64 {
        self.target_bases
    }

    #[must_use]
    pub const fn region_size_bases(&self) -> u64 {
        self.region_size_bases
    }

    #[must_use]
    pub fn segments_in_selection_order(&self) -> &[SelectedSegment] {
        &self.segments
    }
}

/// Segmento ordenado por origen y proyectado sobre la referencia derivada.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MaterializedSegment {
    materialization_order: u64,
    selection_order: u64,
    class: ReductionSegmentClass,
    sequence_class: SequenceClass,
    source_interval: GenomicInterval,
    derived_interval: GenomicInterval,
}

impl MaterializedSegment {
    #[must_use]
    pub const fn materialization_order(&self) -> u64 {
        self.materialization_order
    }

    #[must_use]
    pub const fn selection_order(&self) -> u64 {
        self.selection_order
    }

    #[must_use]
    pub const fn class(&self) -> ReductionSegmentClass {
        self.class
    }

    #[must_use]
    pub const fn sequence_class(&self) -> SequenceClass {
        self.sequence_class
    }

    #[must_use]
    pub const fn source_interval(&self) -> &GenomicInterval {
        &self.source_interval
    }

    #[must_use]
    pub const fn derived_interval(&self) -> &GenomicInterval {
        &self.derived_interval
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MaterializedPlan {
    source_reference_id: ReferenceId,
    derived_reference_id: ReferenceId,
    segments: Vec<MaterializedSegment>,
}

impl MaterializedPlan {
    #[must_use]
    pub const fn source_reference_id(&self) -> &ReferenceId {
        &self.source_reference_id
    }

    #[must_use]
    pub const fn derived_reference_id(&self) -> &ReferenceId {
        &self.derived_reference_id
    }

    #[must_use]
    pub fn segments(&self) -> &[MaterializedSegment] {
        &self.segments
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct FreeGap {
    sequence_index: usize,
    start: u64,
    end: u64,
}

impl FreeGap {
    const fn len(self) -> u64 {
        self.end - self.start
    }
}

/// Ejecuta el preflight completo y, solo si es viable, consume el RNG local.
pub fn select_reduced_reference(
    source_reference_id: &ReferenceId,
    sequences: &[EligibleSequence],
    target_bases: u64,
    region_size_bases: u64,
    rng: &mut LocalRng,
) -> Result<SelectionPlan, ReductionError> {
    preflight_reduced_reference(sequences, target_bases, region_size_bases)?;
    let full_count = target_bases / region_size_bases;
    let residual = target_bases % region_size_bases;
    let mut gaps = sequences
        .iter()
        .enumerate()
        .map(|(sequence_index, sequence)| FreeGap {
            sequence_index,
            start: 0,
            end: sequence.length,
        })
        .collect::<Vec<_>>();

    let segment_count = full_count
        .checked_add(u64::from(residual > 0))
        .ok_or(ReductionError::ArithmeticOverflow)?;
    let capacity = usize::try_from(segment_count).map_err(|_| ReductionError::TooManySegments)?;
    let mut segments = Vec::new();
    segments
        .try_reserve_exact(capacity)
        .map_err(|_| ReductionError::AllocationFailed)?;
    gaps.try_reserve_exact(capacity)
        .map_err(|_| ReductionError::AllocationFailed)?;

    for selected_full in 0..full_count {
        let remaining_full = full_count - selected_full - 1;
        let (gap_index, start) =
            choose_safe_full_start(&gaps, remaining_full, residual, region_size_bases, rng)?;
        append_selected_segment(
            source_reference_id,
            sequences,
            &mut gaps,
            &mut segments,
            gap_index,
            start,
            region_size_bases,
            ReductionSegmentClass::Full,
        )?;
    }

    if residual > 0 {
        let (gap_index, start) = choose_any_start(&gaps, residual, rng)?;
        append_selected_segment(
            source_reference_id,
            sequences,
            &mut gaps,
            &mut segments,
            gap_index,
            start,
            residual,
            ReductionSegmentClass::Residual,
        )?;
    }

    Ok(SelectionPlan {
        source_reference_id: source_reference_id.clone(),
        target_bases,
        region_size_bases,
        segments,
    })
}

/// Demuestra la viabilidad completa sin recibir ni consumir una fuente aleatoria.
pub fn preflight_reduced_reference(
    sequences: &[EligibleSequence],
    target_bases: u64,
    region_size_bases: u64,
) -> Result<(), ReductionError> {
    let eligible = validate_request(sequences, target_bases, region_size_bases)?;
    if target_bases >= eligible {
        return Err(ReductionError::TargetNotSmallerThanEligible {
            target: target_bases,
            eligible,
        });
    }

    let full_count = target_bases / region_size_bases;
    let residual = target_bases % region_size_bases;
    let gaps = sequences
        .iter()
        .enumerate()
        .map(|(sequence_index, sequence)| FreeGap {
            sequence_index,
            start: 0,
            end: sequence.length,
        })
        .collect::<Vec<_>>();

    if !layout_feasible(&gaps, full_count, residual, region_size_bases)? {
        return Err(ReductionError::InfeasibleLayout);
    }
    Ok(())
}

fn validate_request(
    sequences: &[EligibleSequence],
    target_bases: u64,
    region_size_bases: u64,
) -> Result<u64, ReductionError> {
    if sequences.is_empty() {
        return Err(ReductionError::EmptyEligibleSpace);
    }
    if target_bases == 0 {
        return Err(ReductionError::ZeroTarget);
    }
    if region_size_bases == 0 {
        return Err(ReductionError::ZeroRegionSize);
    }
    let mut ids = BTreeSet::new();
    let mut eligible = 0_u64;
    for sequence in sequences {
        if sequence.length == 0 {
            return Err(ReductionError::ZeroLengthSequence(
                sequence.id.as_str().to_owned(),
            ));
        }
        if !ids.insert(sequence.id.as_str()) {
            return Err(ReductionError::DuplicateSequenceId(
                sequence.id.as_str().to_owned(),
            ));
        }
        eligible = eligible
            .checked_add(sequence.length)
            .ok_or(ReductionError::ArithmeticOverflow)?;
    }
    Ok(eligible)
}

fn layout_feasible(
    gaps: &[FreeGap],
    full_count: u64,
    residual: u64,
    region_size: u64,
) -> Result<bool, ReductionError> {
    let mut capacity = 0_u64;
    let mut minimum_residual_loss = None;
    for gap in gaps {
        let length = gap.len();
        let full_capacity = length / region_size;
        capacity = capacity
            .checked_add(full_capacity)
            .ok_or(ReductionError::ArithmeticOverflow)?;
        if residual > 0 && length >= residual {
            let capacity_with_residual = (length - residual) / region_size;
            let loss = full_capacity - capacity_with_residual;
            minimum_residual_loss =
                Some(minimum_residual_loss.map_or(loss, |old: u64| old.min(loss)));
        }
    }
    if residual == 0 {
        return Ok(capacity >= full_count);
    }
    Ok(minimum_residual_loss
        .and_then(|loss| capacity.checked_sub(loss))
        .is_some_and(|available| available >= full_count))
}

#[derive(Clone, Copy)]
enum SplitRequirement {
    Any,
    HostResidual,
    ZeroLossResidual,
}

#[derive(Clone, Copy, Debug)]
struct ResidueRange {
    start: u64,
    end: u64,
}

#[derive(Clone, Copy, Debug)]
struct CandidatePiece {
    low: u64,
    high: u64,
    residue: ResidueRange,
}

#[derive(Debug)]
struct CandidateGroup {
    gap_index: usize,
    maximum_offset: u64,
    modulus: u64,
    pieces: Vec<CandidatePiece>,
    count: u64,
}

impl CandidateGroup {
    fn new(
        gap_index: usize,
        maximum_offset: u64,
        modulus: u64,
        pieces: Vec<CandidatePiece>,
    ) -> Result<Option<Self>, ReductionError> {
        let mut count = 0_u64;
        for piece in &pieces {
            count = count
                .checked_add(count_piece(*piece, modulus)?)
                .ok_or(ReductionError::ArithmeticOverflow)?;
        }
        Ok((count > 0).then_some(Self {
            gap_index,
            maximum_offset,
            modulus,
            pieces,
            count,
        }))
    }

    fn nth(&self, index: u64) -> Result<u64, ReductionError> {
        if index >= self.count {
            return Err(ReductionError::InfeasibleLayout);
        }
        let mut low = 0_u64;
        let mut high = self.maximum_offset;
        while low < high {
            let middle = low + (high - low) / 2;
            if self.count_through(middle)? > index {
                high = middle;
            } else {
                low = middle
                    .checked_add(1)
                    .ok_or(ReductionError::ArithmeticOverflow)?;
            }
        }
        Ok(low)
    }

    fn count_through(&self, maximum: u64) -> Result<u64, ReductionError> {
        let mut count = 0_u64;
        for piece in &self.pieces {
            if piece.low > maximum {
                continue;
            }
            let clipped = CandidatePiece {
                low: piece.low,
                high: piece.high.min(maximum),
                residue: piece.residue,
            };
            if clipped.low <= clipped.high {
                count = count
                    .checked_add(count_piece(clipped, self.modulus)?)
                    .ok_or(ReductionError::ArithmeticOverflow)?;
            }
        }
        Ok(count)
    }
}

fn choose_safe_full_start(
    gaps: &[FreeGap],
    remaining_full: u64,
    residual: u64,
    region_size: u64,
    rng: &mut LocalRng,
) -> Result<(usize, u64), ReductionError> {
    let total_capacity = gaps.iter().try_fold(0_u64, |total, gap| {
        total
            .checked_add(gap.len() / region_size)
            .ok_or(ReductionError::ArithmeticOverflow)
    })?;
    let mut groups = Vec::new();
    let mut total_candidates = 0_u64;
    for (gap_index, gap) in gaps.iter().copied().enumerate() {
        if gap.len() < region_size {
            continue;
        }
        if let Some(group) = safe_candidates_for_gap(
            gaps,
            gap_index,
            total_capacity,
            remaining_full,
            residual,
            region_size,
        )? {
            total_candidates = total_candidates
                .checked_add(group.count)
                .ok_or(ReductionError::ArithmeticOverflow)?;
            groups.push(group);
        }
    }
    choose_from_groups(&groups, total_candidates, gaps, rng)
}

fn choose_any_start(
    gaps: &[FreeGap],
    segment_length: u64,
    rng: &mut LocalRng,
) -> Result<(usize, u64), ReductionError> {
    let mut groups = Vec::new();
    let mut total_candidates = 0_u64;
    for (gap_index, gap) in gaps.iter().copied().enumerate() {
        if gap.len() < segment_length {
            continue;
        }
        let maximum_offset = gap.len() - segment_length;
        let piece = CandidatePiece {
            low: 0,
            high: maximum_offset,
            residue: ResidueRange { start: 0, end: 0 },
        };
        if let Some(group) = CandidateGroup::new(gap_index, maximum_offset, 1, vec![piece])? {
            total_candidates = total_candidates
                .checked_add(group.count)
                .ok_or(ReductionError::ArithmeticOverflow)?;
            groups.push(group);
        }
    }
    choose_from_groups(&groups, total_candidates, gaps, rng)
}

fn choose_from_groups(
    groups: &[CandidateGroup],
    total_candidates: u64,
    gaps: &[FreeGap],
    rng: &mut LocalRng,
) -> Result<(usize, u64), ReductionError> {
    if total_candidates == 0 {
        return Err(ReductionError::InfeasibleLayout);
    }
    let mut selected = uniform_index(rng, total_candidates);
    for group in groups {
        if selected < group.count {
            let local_offset = group.nth(selected)?;
            let start = gaps[group.gap_index]
                .start
                .checked_add(local_offset)
                .ok_or(ReductionError::ArithmeticOverflow)?;
            return Ok((group.gap_index, start));
        }
        selected -= group.count;
    }
    Err(ReductionError::InfeasibleLayout)
}

fn safe_candidates_for_gap(
    gaps: &[FreeGap],
    gap_index: usize,
    total_capacity: u64,
    remaining_full: u64,
    residual: u64,
    region_size: u64,
) -> Result<Option<CandidateGroup>, ReductionError> {
    let gap = gaps[gap_index];
    let length = gap.len();
    let maximum_offset = length - region_size;
    let other_capacity = total_capacity - length / region_size;
    let other_residual_loss = minimum_residual_loss(gaps, gap_index, residual, region_size);
    let high_split_capacity = maximum_offset / region_size;
    let boundary = maximum_offset % region_size;
    let mut pieces = Vec::new();

    add_capacity_class(
        &mut pieces,
        ResidueRange {
            start: 0,
            end: boundary,
        },
        other_capacity
            .checked_add(high_split_capacity)
            .ok_or(ReductionError::ArithmeticOverflow)?,
        maximum_offset,
        region_size,
        remaining_full,
        residual,
        other_residual_loss,
    );

    if boundary + 1 < region_size && high_split_capacity > 0 {
        add_capacity_class(
            &mut pieces,
            ResidueRange {
                start: boundary + 1,
                end: region_size - 1,
            },
            other_capacity
                .checked_add(high_split_capacity - 1)
                .ok_or(ReductionError::ArithmeticOverflow)?,
            maximum_offset,
            region_size,
            remaining_full,
            residual,
            other_residual_loss,
        );
    }

    CandidateGroup::new(gap_index, maximum_offset, region_size, pieces)
}

fn minimum_residual_loss(
    gaps: &[FreeGap],
    excluded: usize,
    residual: u64,
    region_size: u64,
) -> Option<u64> {
    if residual == 0 {
        return Some(0);
    }
    gaps.iter()
        .enumerate()
        .filter(|(index, _)| *index != excluded)
        .filter_map(|(_, gap)| {
            let length = gap.len();
            (length >= residual).then(|| length / region_size - (length - residual) / region_size)
        })
        .min()
}

#[allow(clippy::too_many_arguments)]
fn add_capacity_class(
    pieces: &mut Vec<CandidatePiece>,
    capacity_residues: ResidueRange,
    capacity: u64,
    maximum_offset: u64,
    modulus: u64,
    remaining_full: u64,
    residual: u64,
    other_residual_loss: Option<u64>,
) {
    if capacity < remaining_full {
        return;
    }
    let requirement = if residual == 0 || other_residual_loss == Some(0) {
        SplitRequirement::Any
    } else if capacity > remaining_full {
        if other_residual_loss.is_some() {
            SplitRequirement::Any
        } else {
            SplitRequirement::HostResidual
        }
    } else {
        SplitRequirement::ZeroLossResidual
    };

    let absolute_ranges = match requirement {
        SplitRequirement::HostResidual if maximum_offset < residual.saturating_mul(2) => {
            let mut ranges = Vec::new();
            if maximum_offset >= residual {
                ranges.push((0, maximum_offset - residual));
                ranges.push((residual, maximum_offset));
            }
            ranges
        }
        _ => vec![(0, maximum_offset)],
    };
    let residue_ranges = match requirement {
        SplitRequirement::ZeroLossResidual => intersect_ranges(
            &[capacity_residues],
            &zero_loss_residues(maximum_offset % modulus, residual, modulus),
        ),
        _ => vec![capacity_residues],
    };
    for (low, high) in absolute_ranges {
        for residue in &residue_ranges {
            pieces.push(CandidatePiece {
                low,
                high,
                residue: *residue,
            });
        }
    }
}

fn zero_loss_residues(offset_remainder: u64, residual: u64, modulus: u64) -> Vec<ResidueRange> {
    if residual == 0 {
        return vec![ResidueRange {
            start: 0,
            end: modulus - 1,
        }];
    }
    let first_forbidden = ResidueRange {
        start: 0,
        end: residual - 1,
    };
    let backwards = (residual - 1) % modulus;
    let circular_start = if offset_remainder >= backwards {
        offset_remainder - backwards
    } else {
        modulus - (backwards - offset_remainder)
    };
    let circular_forbidden = if circular_start <= offset_remainder {
        vec![ResidueRange {
            start: circular_start,
            end: offset_remainder,
        }]
    } else {
        vec![
            ResidueRange {
                start: 0,
                end: offset_remainder,
            },
            ResidueRange {
                start: circular_start,
                end: modulus - 1,
            },
        ]
    };
    let forbidden = intersect_ranges(&[first_forbidden], &circular_forbidden);
    subtract_ranges(
        &[ResidueRange {
            start: 0,
            end: modulus - 1,
        }],
        &forbidden,
    )
}

fn intersect_ranges(left: &[ResidueRange], right: &[ResidueRange]) -> Vec<ResidueRange> {
    let mut intersections = Vec::new();
    for a in left {
        for b in right {
            let start = a.start.max(b.start);
            let end = a.end.min(b.end);
            if start <= end {
                intersections.push(ResidueRange { start, end });
            }
        }
    }
    intersections
}

fn subtract_ranges(source: &[ResidueRange], removed: &[ResidueRange]) -> Vec<ResidueRange> {
    let mut result = source.to_vec();
    for cut in removed {
        let mut next = Vec::new();
        for range in result {
            if cut.end < range.start || cut.start > range.end {
                next.push(range);
                continue;
            }
            if cut.start > range.start {
                next.push(ResidueRange {
                    start: range.start,
                    end: cut.start - 1,
                });
            }
            if cut.end < range.end {
                next.push(ResidueRange {
                    start: cut.end + 1,
                    end: range.end,
                });
            }
        }
        result = next;
    }
    result
}

fn count_piece(piece: CandidatePiece, modulus: u64) -> Result<u64, ReductionError> {
    let through_high = count_residues_through(piece.high, piece.residue, modulus)?;
    let before_low = if piece.low == 0 {
        0
    } else {
        count_residues_through(piece.low - 1, piece.residue, modulus)?
    };
    through_high
        .checked_sub(before_low)
        .ok_or(ReductionError::ArithmeticOverflow)
}

fn count_residues_through(
    maximum: u64,
    residues: ResidueRange,
    modulus: u64,
) -> Result<u64, ReductionError> {
    let width = residues.end - residues.start + 1;
    let complete_cycles = maximum / modulus;
    let remainder = maximum % modulus;
    let tail = if remainder < residues.start {
        0
    } else {
        remainder.min(residues.end) - residues.start + 1
    };
    complete_cycles
        .checked_mul(width)
        .and_then(|count| count.checked_add(tail))
        .ok_or(ReductionError::ArithmeticOverflow)
}

fn uniform_index(rng: &mut LocalRng, upper_exclusive: u64) -> u64 {
    debug_assert!(upper_exclusive > 0);
    let rejection_threshold = upper_exclusive.wrapping_neg() % upper_exclusive;
    loop {
        let value = rng.next_u64();
        if value >= rejection_threshold {
            return value % upper_exclusive;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn append_selected_segment(
    source_reference_id: &ReferenceId,
    sequences: &[EligibleSequence],
    gaps: &mut Vec<FreeGap>,
    segments: &mut Vec<SelectedSegment>,
    gap_index: usize,
    start: u64,
    length: u64,
    class: ReductionSegmentClass,
) -> Result<(), ReductionError> {
    let gap = gaps[gap_index];
    let end = start
        .checked_add(length)
        .ok_or(ReductionError::ArithmeticOverflow)?;
    if start < gap.start || end > gap.end {
        return Err(ReductionError::InfeasibleLayout);
    }
    let sequence = &sequences[gap.sequence_index];
    let selection_order = u64::try_from(segments.len())
        .map_err(|_| ReductionError::TooManySegments)?
        .checked_add(1)
        .ok_or(ReductionError::ArithmeticOverflow)?;
    segments.push(SelectedSegment {
        selection_order,
        class,
        sequence_class: sequence.class,
        source_interval: GenomicInterval::new(
            source_reference_id.clone(),
            sequence.id.clone(),
            start,
            end,
            sequence.length,
        )?,
    });

    gaps.remove(gap_index);
    if gap.start < start {
        gaps.push(FreeGap {
            sequence_index: gap.sequence_index,
            start: gap.start,
            end: start,
        });
    }
    if end < gap.end {
        gaps.push(FreeGap {
            sequence_index: gap.sequence_index,
            start: end,
            end: gap.end,
        });
    }
    gaps.sort_by_key(|free| (free.sequence_index, free.start));
    Ok(())
}

/// Ordena por FASTA/coordenada y asigna coordenadas acumuladas por secuencia.
pub fn materialize_selection(
    selection: &SelectionPlan,
    derived_reference_id: &ReferenceId,
    sequences: &[EligibleSequence],
) -> Result<MaterializedPlan, ReductionError> {
    let mut catalog = BTreeMap::new();
    for (index, sequence) in sequences.iter().enumerate() {
        if catalog
            .insert(sequence.id.as_str(), (index, sequence.class))
            .is_some()
        {
            return Err(ReductionError::DuplicateSequenceId(
                sequence.id.as_str().to_owned(),
            ));
        }
    }
    let mut ordered = selection.segments.iter().collect::<Vec<_>>();
    ordered.sort_by_key(|segment| {
        let sequence_order = catalog
            .get(segment.source_interval.sequence_id().as_str())
            .map_or(usize::MAX, |(index, _)| *index);
        (sequence_order, segment.source_interval.start())
    });
    if ordered.iter().any(|segment| {
        catalog
            .get(segment.source_interval.sequence_id().as_str())
            .is_none_or(|(_, class)| *class != segment.sequence_class)
    }) {
        return Err(ReductionError::SelectionDoesNotMatchCatalog);
    }

    let mut derived_lengths = vec![0_u64; sequences.len()];
    for segment in &ordered {
        let index = catalog
            .get(segment.source_interval.sequence_id().as_str())
            .map(|(index, _)| *index)
            .ok_or(ReductionError::SelectionDoesNotMatchCatalog)?;
        derived_lengths[index] = derived_lengths[index]
            .checked_add(segment.source_interval.len())
            .ok_or(ReductionError::ArithmeticOverflow)?;
    }
    let mut offsets = vec![0_u64; sequences.len()];
    let mut materialized = Vec::new();
    materialized
        .try_reserve_exact(ordered.len())
        .map_err(|_| ReductionError::AllocationFailed)?;
    for (index, segment) in ordered.into_iter().enumerate() {
        let sequence_index = catalog
            .get(segment.source_interval.sequence_id().as_str())
            .map(|(index, _)| *index)
            .ok_or(ReductionError::SelectionDoesNotMatchCatalog)?;
        let start = offsets[sequence_index];
        let end = start
            .checked_add(segment.source_interval.len())
            .ok_or(ReductionError::ArithmeticOverflow)?;
        offsets[sequence_index] = end;
        materialized.push(MaterializedSegment {
            materialization_order: u64::try_from(index)
                .map_err(|_| ReductionError::TooManySegments)?
                .checked_add(1)
                .ok_or(ReductionError::ArithmeticOverflow)?,
            selection_order: segment.selection_order,
            class: segment.class,
            sequence_class: segment.sequence_class,
            source_interval: segment.source_interval.clone(),
            derived_interval: GenomicInterval::new(
                derived_reference_id.clone(),
                segment.source_interval.sequence_id().clone(),
                start,
                end,
                derived_lengths[sequence_index],
            )?,
        });
    }

    Ok(MaterializedPlan {
        source_reference_id: selection.source_reference_id.clone(),
        derived_reference_id: derived_reference_id.clone(),
        segments: materialized,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::randomness::{RootSeed, SubseedRequest, derive_subseed};

    fn reference_id(value: &str) -> ReferenceId {
        ReferenceId::parse(value, "reduction-test").expect("reference ID")
    }

    fn sequence(id: &str, class: SequenceClass, length: u64) -> EligibleSequence {
        EligibleSequence::new(
            SequenceId::parse(id, "ref-source").expect("sequence ID"),
            class,
            length,
        )
        .expect("eligible sequence")
    }

    fn rng(byte: u8) -> LocalRng {
        let root = RootSeed::from_bytes([byte; 32]);
        let request = SubseedRequest::from_scope(
            "reference_reduction",
            &[
                ("source_reference_id", "ref-source"),
                ("derived_reference_id", "ref-derived"),
                ("algorithm", REDUCED_REFERENCE_SELECTION_ALGORITHM),
            ],
        )
        .expect("subseed request");
        derive_subseed(root, &request).expect("subseed").local_rng()
    }

    fn rng_from_number(value: u64) -> LocalRng {
        let mut bytes = [0_u8; 32];
        bytes[..8].copy_from_slice(&value.to_be_bytes());
        let root = RootSeed::from_bytes(bytes);
        let request = SubseedRequest::from_scope(
            "reference_reduction",
            &[
                ("source_reference_id", "ref-source"),
                ("derived_reference_id", "ref-derived"),
                ("algorithm", REDUCED_REFERENCE_SELECTION_ALGORITHM),
            ],
        )
        .expect("subseed request");
        derive_subseed(root, &request).expect("subseed").local_rng()
    }

    #[test]
    fn creates_exact_non_overlapping_plan_with_one_residual() {
        let sequences = vec![
            sequence("chr1", SequenceClass::Chromosome, 11),
            sequence("scaf1", SequenceClass::Scaffold, 8),
        ];
        let mut local_rng = rng(7);
        let plan = select_reduced_reference(
            &reference_id("ref-source"),
            &sequences,
            13,
            5,
            &mut local_rng,
        )
        .expect("selection plan");
        assert_eq!(plan.segments.len(), 3);
        assert_eq!(
            plan.segments
                .iter()
                .map(|segment| segment.source_interval.len())
                .sum::<u64>(),
            13
        );
        assert_eq!(
            plan.segments
                .iter()
                .filter(|segment| segment.class == ReductionSegmentClass::Residual)
                .count(),
            1
        );
        for (index, first) in plan.segments.iter().enumerate() {
            for second in &plan.segments[index + 1..] {
                if first.source_interval.sequence_id() == second.source_interval.sequence_id() {
                    assert!(
                        first.source_interval.end() <= second.source_interval.start()
                            || second.source_interval.end() <= first.source_interval.start()
                    );
                }
            }
        }
    }

    #[test]
    fn exact_multiple_has_no_residual_segment() {
        let sequences = vec![sequence("chr1", SequenceClass::Chromosome, 20)];
        let mut rng = rng(19);
        let plan = select_reduced_reference(&reference_id("ref"), &sequences, 10, 5, &mut rng)
            .expect("exact multiple plan");
        assert_eq!(plan.segments_in_selection_order().len(), 2);
        assert!(
            plan.segments_in_selection_order()
                .iter()
                .all(|segment| segment.class() == ReductionSegmentClass::Full)
        );
    }

    #[test]
    fn preflight_failure_does_not_consume_rng() {
        let sequences = vec![
            sequence("chr1", SequenceClass::Chromosome, 3),
            sequence("chr2", SequenceClass::Chromosome, 3),
        ];
        let mut after_failure = rng(9);
        assert_eq!(
            select_reduced_reference(
                &reference_id("ref-source"),
                &sequences,
                5,
                4,
                &mut after_failure,
            ),
            Err(ReductionError::InfeasibleLayout)
        );
        let mut untouched = rng(9);
        let first = select_reduced_reference(
            &reference_id("ref-source"),
            &sequences,
            2,
            2,
            &mut after_failure,
        )
        .expect("valid plan");
        let second = select_reduced_reference(
            &reference_id("ref-source"),
            &sequences,
            2,
            2,
            &mut untouched,
        )
        .expect("same valid plan");
        assert_eq!(first, second);
    }

    #[test]
    fn materialization_uses_fasta_order_and_accumulates_per_sequence() {
        let sequences = vec![
            sequence("chr1", SequenceClass::Chromosome, 12),
            sequence("scaf1", SequenceClass::Scaffold, 8),
        ];
        let mut local_rng = rng(17);
        let selection = select_reduced_reference(
            &reference_id("ref-source"),
            &sequences,
            13,
            5,
            &mut local_rng,
        )
        .expect("selection");
        let materialized =
            materialize_selection(&selection, &reference_id("ref-derived"), &sequences)
                .expect("materialized plan");
        assert!(materialized.segments.windows(2).all(|pair| {
            let first_sequence = sequences
                .iter()
                .position(|sequence| sequence.id() == pair[0].source_interval().sequence_id())
                .expect("first sequence");
            let second_sequence = sequences
                .iter()
                .position(|sequence| sequence.id() == pair[1].source_interval().sequence_id())
                .expect("second sequence");
            first_sequence < second_sequence
                || (first_sequence == second_sequence
                    && pair[0].source_interval().start() < pair[1].source_interval().start())
        }));
        for window in materialized.segments.windows(2) {
            if window[0].derived_interval().sequence_id()
                == window[1].derived_interval().sequence_id()
            {
                assert_eq!(
                    window[0].derived_interval().end(),
                    window[1].derived_interval().start()
                );
            }
        }
    }

    #[test]
    fn safe_candidate_groups_match_exhaustive_oracle_for_small_gaps() {
        for length in 1_u64..=14 {
            for region_size in 1_u64..=length {
                let gap = FreeGap {
                    sequence_index: 0,
                    start: 0,
                    end: length,
                };
                for remaining_full in 0..=3 {
                    for residual in 0..region_size {
                        let total_capacity = length / region_size;
                        let group = safe_candidates_for_gap(
                            &[gap],
                            0,
                            total_capacity,
                            remaining_full,
                            residual,
                            region_size,
                        )
                        .expect("candidate group");
                        let expected = (0..=length - region_size)
                            .filter(|start| {
                                let after = [
                                    FreeGap {
                                        sequence_index: 0,
                                        start: 0,
                                        end: *start,
                                    },
                                    FreeGap {
                                        sequence_index: 0,
                                        start: start + region_size,
                                        end: length,
                                    },
                                ];
                                layout_feasible(&after, remaining_full, residual, region_size)
                                    .expect("oracle")
                            })
                            .collect::<Vec<_>>();
                        let actual = group
                            .as_ref()
                            .map(|candidates| {
                                (0..candidates.count)
                                    .map(|index| candidates.nth(index).expect("candidate"))
                                    .collect::<Vec<_>>()
                            })
                            .unwrap_or_default();
                        assert_eq!(
                            actual, expected,
                            "L={length}, R={region_size}, k={remaining_full}, r={residual}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn safe_candidate_groups_match_exhaustive_oracle_with_other_gaps() {
        for first_length in 1_u64..=8 {
            for second_length in 1_u64..=8 {
                for region_size in 1_u64..=first_length.max(second_length) {
                    let gaps = [
                        FreeGap {
                            sequence_index: 0,
                            start: 0,
                            end: first_length,
                        },
                        FreeGap {
                            sequence_index: 1,
                            start: 0,
                            end: second_length,
                        },
                    ];
                    let total_capacity =
                        gaps.iter().map(|gap| gap.len() / region_size).sum::<u64>();
                    for gap_index in 0..2 {
                        if gaps[gap_index].len() < region_size {
                            continue;
                        }
                        for remaining_full in 0..=3 {
                            for residual in 0..region_size {
                                let group = safe_candidates_for_gap(
                                    &gaps,
                                    gap_index,
                                    total_capacity,
                                    remaining_full,
                                    residual,
                                    region_size,
                                )
                                .expect("candidate group");
                                let selected_gap = gaps[gap_index];
                                let expected = (0..=selected_gap.len() - region_size)
                                    .filter(|offset| {
                                        let start = selected_gap.start + offset;
                                        let mut after = gaps.to_vec();
                                        after.remove(gap_index);
                                        if selected_gap.start < start {
                                            after.push(FreeGap {
                                                sequence_index: selected_gap.sequence_index,
                                                start: selected_gap.start,
                                                end: start,
                                            });
                                        }
                                        let end = start + region_size;
                                        if end < selected_gap.end {
                                            after.push(FreeGap {
                                                sequence_index: selected_gap.sequence_index,
                                                start: end,
                                                end: selected_gap.end,
                                            });
                                        }
                                        layout_feasible(
                                            &after,
                                            remaining_full,
                                            residual,
                                            region_size,
                                        )
                                        .expect("oracle")
                                    })
                                    .collect::<Vec<_>>();
                                let actual = group
                                    .as_ref()
                                    .map(|candidates| {
                                        (0..candidates.count)
                                            .map(|index| candidates.nth(index).expect("candidate"))
                                            .collect::<Vec<_>>()
                                    })
                                    .unwrap_or_default();
                                assert_eq!(
                                    actual, expected,
                                    "L1={first_length}, L2={second_length}, gap={gap_index}, R={region_size}, k={remaining_full}, r={residual}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn starts_are_uniform_over_bases_not_over_sequences() {
        const REPETITIONS: u64 = 7_000;
        const EXPECTED_PER_START: u64 = 1_000;
        const MAXIMUM_ABSOLUTE_DEVIATION: u64 = 150;
        let sequences = vec![
            sequence("chr1", SequenceClass::Chromosome, 3),
            sequence("scaf1", SequenceClass::Scaffold, 4),
        ];
        let mut counts = [0_u64; 7];
        for seed in 0..REPETITIONS {
            let plan = select_reduced_reference(
                &reference_id("ref-source"),
                &sequences,
                1,
                1,
                &mut rng_from_number(seed),
            )
            .expect("one-base plan");
            let segment = &plan.segments[0];
            let offset = usize::try_from(segment.source_interval.start()).expect("small offset");
            let index = if segment.source_interval.sequence_id().as_str() == "chr1" {
                offset
            } else {
                3 + offset
            };
            counts[index] += 1;
        }
        for count in counts {
            assert!(count.abs_diff(EXPECTED_PER_START) <= MAXIMUM_ABSOLUTE_DEVIATION);
        }
    }

    #[test]
    fn request_boundaries_are_rejected() {
        let sequences = vec![sequence("chr1", SequenceClass::Chromosome, 10)];
        assert_eq!(
            select_reduced_reference(&reference_id("ref-source"), &sequences, 10, 3, &mut rng(1),),
            Err(ReductionError::TargetNotSmallerThanEligible {
                target: 10,
                eligible: 10
            })
        );
        assert_eq!(
            select_reduced_reference(&reference_id("ref-source"), &sequences, 1, 0, &mut rng(1),),
            Err(ReductionError::ZeroRegionSize)
        );
    }

    #[test]
    fn uniform_integer_mapping_has_a_frozen_vector() {
        let mut local_rng = rng(23);
        let observed = (0..8)
            .map(|_| uniform_index(&mut local_rng, 7))
            .collect::<Vec<_>>();
        assert_eq!(observed, vec![5, 1, 3, 3, 3, 4, 1, 3]);
    }

    #[test]
    fn complete_selection_has_a_frozen_vector() {
        let sequences = vec![sequence("chr1", SequenceClass::Chromosome, 6)];
        let plan =
            select_reduced_reference(&reference_id("ref-source"), &sequences, 5, 3, &mut rng(9))
                .expect("selection");
        let observed = plan
            .segments
            .iter()
            .map(|segment| {
                (
                    segment.class,
                    segment.source_interval.start(),
                    segment.source_interval.end(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            observed,
            vec![
                (ReductionSegmentClass::Full, 0, 3),
                (ReductionSegmentClass::Residual, 3, 5),
            ]
        );
    }
}
