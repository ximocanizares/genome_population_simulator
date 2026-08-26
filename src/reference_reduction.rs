//! Planificación reproducible de referencias reducidas.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};

use crate::coordinates::{CoordinateError, GenomicInterval};
use crate::identity::{ReferenceId, SequenceId};
use crate::randomness::LocalRng;
use crate::reference::SequenceClass;

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
    let full = target_bases / region_size_bases;
    let residual = target_bases % region_size_bases;
    let capacity = total_capacity(sequences, region_size_bases)?;
    if residual == 0 {
        return (capacity >= full)
            .then_some(())
            .ok_or(ReductionError::InfeasibleLayout);
    }
    sequences
        .iter()
        .any(|sequence| {
            sequence.length >= residual
                && capacity - sequence.length / region_size_bases
                    + (sequence.length - residual) / region_size_bases
                    >= full
        })
        .then_some(())
        .ok_or(ReductionError::InfeasibleLayout)
}

/// Selector directo proporcional: asigna bloques completos por longitud con
/// límites de capacidad y después materializa posiciones sin solapes.
pub fn select_reduced_reference(
    source_reference_id: &ReferenceId,
    sequences: &[EligibleSequence],
    target_bases: u64,
    region_size_bases: u64,
    rng: &mut LocalRng,
) -> Result<SelectionPlan, ReductionError> {
    preflight_reduced_reference(sequences, target_bases, region_size_bases)?;
    let full = target_bases / region_size_bases;
    let residual = target_bases % region_size_bases;
    let host = if residual == 0 {
        None
    } else {
        Some(choose_residual_host(
            sequences,
            full,
            residual,
            region_size_bases,
            rng,
        )?)
    };
    let capacities = sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| {
            let reserved = if Some(index) == host { residual } else { 0 };
            (sequence.length - reserved) / region_size_bases
        })
        .collect::<Vec<_>>();
    let allocation = allocate_proportional(sequences, &capacities, full)?;
    let mut segments = Vec::new();
    segments
        .try_reserve_exact(
            usize::try_from(full + u64::from(residual > 0))
                .map_err(|_| ReductionError::TooManySegments)?,
        )
        .map_err(|_| ReductionError::AllocationFailed)?;
    for (index, sequence) in sequences.iter().enumerate() {
        append_layout(
            source_reference_id,
            sequence,
            allocation[index],
            (Some(index) == host).then_some(residual),
            region_size_bases,
            &mut segments,
            rng,
        )?;
    }
    let mut order = 1_u64;
    for segment in &mut segments {
        if segment.class == ReductionSegmentClass::Full {
            segment.selection_order = order;
            order += 1;
        }
    }
    for segment in &mut segments {
        if segment.class == ReductionSegmentClass::Residual {
            segment.selection_order = order;
        }
    }
    segments.sort_by_key(SelectedSegment::selection_order);
    Ok(SelectionPlan {
        source_reference_id: source_reference_id.clone(),
        target_bases,
        region_size_bases,
        segments,
    })
}

fn validate_request(
    sequences: &[EligibleSequence],
    target: u64,
    region: u64,
) -> Result<u64, ReductionError> {
    if sequences.is_empty() {
        return Err(ReductionError::EmptyEligibleSpace);
    }
    if target == 0 {
        return Err(ReductionError::ZeroTarget);
    }
    if region == 0 {
        return Err(ReductionError::ZeroRegionSize);
    }
    let mut ids = BTreeSet::new();
    sequences.iter().try_fold(0_u64, |sum, sequence| {
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
        sum.checked_add(sequence.length)
            .ok_or(ReductionError::ArithmeticOverflow)
    })
}

fn total_capacity(sequences: &[EligibleSequence], region: u64) -> Result<u64, ReductionError> {
    sequences.iter().try_fold(0_u64, |sum, sequence| {
        sum.checked_add(sequence.length / region)
            .ok_or(ReductionError::ArithmeticOverflow)
    })
}

fn choose_residual_host(
    sequences: &[EligibleSequence],
    full: u64,
    residual: u64,
    region: u64,
    rng: &mut LocalRng,
) -> Result<usize, ReductionError> {
    let total = total_capacity(sequences, region)?;
    let candidates = sequences
        .iter()
        .enumerate()
        .filter_map(|(index, sequence)| {
            (sequence.length >= residual
                && total - sequence.length / region + (sequence.length - residual) / region >= full)
                .then_some((index, sequence.length - residual + 1))
        })
        .collect::<Vec<_>>();
    let weight = candidates.iter().try_fold(0_u64, |sum, (_, value)| {
        sum.checked_add(*value)
            .ok_or(ReductionError::ArithmeticOverflow)
    })?;
    let mut selected = uniform_index(rng, weight);
    for (index, count) in candidates {
        if selected < count {
            return Ok(index);
        }
        selected -= count;
    }
    Err(ReductionError::InfeasibleLayout)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Priority {
    deficit: i128,
    index: usize,
}
impl Ord for Priority {
    fn cmp(&self, other: &Self) -> Ordering {
        self.deficit
            .cmp(&other.deficit)
            .then_with(|| other.index.cmp(&self.index))
    }
}
impl PartialOrd for Priority {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn allocate_proportional(
    sequences: &[EligibleSequence],
    capacities: &[u64],
    full: u64,
) -> Result<Vec<u64>, ReductionError> {
    if sequences.len() != capacities.len() {
        return Err(ReductionError::SelectionDoesNotMatchCatalog);
    }
    let total = sequences.iter().try_fold(0_u64, |sum, sequence| {
        sum.checked_add(sequence.length)
            .ok_or(ReductionError::ArithmeticOverflow)
    })?;
    let mut allocated = Vec::with_capacity(sequences.len());
    let mut heap = BinaryHeap::new();
    let mut assigned = 0_u64;
    for (index, (sequence, capacity)) in sequences.iter().zip(capacities).enumerate() {
        let requested = u128::from(full) * u128::from(sequence.length);
        let base = u64::try_from(requested / u128::from(total))
            .map_err(|_| ReductionError::ArithmeticOverflow)?
            .min(*capacity);
        assigned += base;
        allocated.push(base);
        if base < *capacity {
            heap.push(Priority {
                deficit: i128::try_from(requested - u128::from(base) * u128::from(total))
                    .map_err(|_| ReductionError::ArithmeticOverflow)?,
                index,
            });
        }
    }
    while assigned < full {
        let priority = heap.pop().ok_or(ReductionError::InfeasibleLayout)?;
        allocated[priority.index] += 1;
        assigned += 1;
        if allocated[priority.index] < capacities[priority.index] {
            heap.push(Priority {
                deficit: priority.deficit - i128::from(total),
                index: priority.index,
            });
        }
    }
    Ok(allocated)
}

fn append_layout(
    source: &ReferenceId,
    sequence: &EligibleSequence,
    full: u64,
    residual: Option<u64>,
    region: u64,
    segments: &mut Vec<SelectedSegment>,
    rng: &mut LocalRng,
) -> Result<(), ReductionError> {
    let residual_length = residual.unwrap_or(0);
    let count = full + u64::from(residual.is_some());
    if count == 0 {
        return Ok(());
    }
    let occupied = full
        .checked_mul(region)
        .and_then(|value| value.checked_add(residual_length))
        .ok_or(ReductionError::ArithmeticOverflow)?;
    let slack = sequence
        .length
        .checked_sub(occupied)
        .ok_or(ReductionError::InfeasibleLayout)?;
    let separators = sample_sorted_distinct(
        rng,
        slack
            .checked_add(count)
            .ok_or(ReductionError::ArithmeticOverflow)?,
        count,
    )?;
    let gaps = gaps_from_separators(slack, count, &separators)?;
    let residual_rank = residual
        .map(|_| uniform_index(rng, full + 1))
        .unwrap_or(u64::MAX);
    let mut start = gaps[0];
    for rank in 0..count {
        let (length, class) = if rank == residual_rank {
            (residual_length, ReductionSegmentClass::Residual)
        } else {
            (region, ReductionSegmentClass::Full)
        };
        let end = start
            .checked_add(length)
            .ok_or(ReductionError::ArithmeticOverflow)?;
        segments.push(SelectedSegment {
            selection_order: 0,
            class,
            sequence_class: sequence.class,
            source_interval: GenomicInterval::new(
                source.clone(),
                sequence.id.clone(),
                start,
                end,
                sequence.length,
            )?,
        });
        start = end
            .checked_add(
                gaps[usize::try_from(rank + 1).map_err(|_| ReductionError::TooManySegments)?],
            )
            .ok_or(ReductionError::ArithmeticOverflow)?;
    }
    Ok(())
}

fn uniform_index(rng: &mut LocalRng, upper: u64) -> u64 {
    let threshold = upper.wrapping_neg() % upper;
    loop {
        let value = rng.next_u64();
        if value >= threshold {
            return value % upper;
        }
    }
}

/// Floyd sampling draws `count` separators without iterating through slack.
fn sample_sorted_distinct(
    rng: &mut LocalRng,
    universe: u64,
    count: u64,
) -> Result<Vec<u64>, ReductionError> {
    if count > universe {
        return Err(ReductionError::InfeasibleLayout);
    }
    let mut selected = BTreeSet::new();
    for value in (universe - count)..universe {
        let candidate = uniform_index(rng, value + 1);
        if !selected.insert(candidate) {
            selected.insert(value);
        }
    }
    Ok(selected.into_iter().collect())
}

fn gaps_from_separators(
    slack: u64,
    count: u64,
    separators: &[u64],
) -> Result<Vec<u64>, ReductionError> {
    if separators.len() != usize::try_from(count).map_err(|_| ReductionError::TooManySegments)? {
        return Err(ReductionError::InfeasibleLayout);
    }
    if separators.is_empty() {
        return Ok(vec![slack]);
    }
    let mut gaps = Vec::with_capacity(separators.len() + 1);
    gaps.push(separators[0]);
    for pair in separators.windows(2) {
        gaps.push(pair[1] - pair[0] - 1);
    }
    let last = slack
        .checked_add(count)
        .and_then(|value| value.checked_sub(1))
        .ok_or(ReductionError::ArithmeticOverflow)?;
    gaps.push(last - separators[separators.len() - 1]);
    Ok(gaps)
}

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
        (
            catalog
                .get(segment.source_interval.sequence_id().as_str())
                .map_or(usize::MAX, |(index, _)| *index),
            segment.source_interval.start(),
        )
    });
    let mut lengths = vec![0_u64; sequences.len()];
    for segment in &ordered {
        let index = catalog
            .get(segment.source_interval.sequence_id().as_str())
            .ok_or(ReductionError::SelectionDoesNotMatchCatalog)?
            .0;
        lengths[index] = lengths[index]
            .checked_add(segment.source_interval.len())
            .ok_or(ReductionError::ArithmeticOverflow)?;
    }
    let mut offsets = vec![0_u64; sequences.len()];
    let mut materialized = Vec::with_capacity(ordered.len());
    for (index, segment) in ordered.into_iter().enumerate() {
        let sequence_index = catalog
            .get(segment.source_interval.sequence_id().as_str())
            .ok_or(ReductionError::SelectionDoesNotMatchCatalog)?
            .0;
        let start = offsets[sequence_index];
        let end = start
            .checked_add(segment.source_interval.len())
            .ok_or(ReductionError::ArithmeticOverflow)?;
        offsets[sequence_index] = end;
        materialized.push(MaterializedSegment {
            materialization_order: u64::try_from(index)
                .map_err(|_| ReductionError::TooManySegments)?
                + 1,
            selection_order: segment.selection_order,
            class: segment.class,
            sequence_class: segment.sequence_class,
            source_interval: segment.source_interval.clone(),
            derived_interval: GenomicInterval::new(
                derived_reference_id.clone(),
                segment.source_interval.sequence_id().clone(),
                start,
                end,
                lengths[sequence_index],
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
    use crate::randomness::{RootSeed, RuntimeSubseedScope, derive_subseed};
    fn reference_id(value: &str) -> ReferenceId {
        ReferenceId::parse(value, "test").unwrap()
    }
    fn sequence(id: &str, class: SequenceClass, length: u64) -> EligibleSequence {
        EligibleSequence::new(SequenceId::parse(id, "test").unwrap(), class, length).unwrap()
    }
    fn rng(seed: u8) -> LocalRng {
        derive_subseed(
            RootSeed::from_bytes([seed; 32]),
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: "source",
                derived_reference_id: "derived",
            },
        )
        .unwrap()
        .local_rng()
    }
    fn rng_from_number(value: u64) -> LocalRng {
        let mut bytes = [0_u8; 32];
        bytes[..8].copy_from_slice(&value.to_be_bytes());
        derive_subseed(
            RootSeed::from_bytes(bytes),
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: "ref-source",
                derived_reference_id: "ref-derived",
            },
        )
        .unwrap()
        .local_rng()
    }
    #[test]
    fn canonical_selector_reproduces_the_frozen_direct_vector() {
        let sequences = vec![
            sequence("chr-a", SequenceClass::Chromosome, 13),
            sequence("chr-b", SequenceClass::Chromosome, 13),
            sequence("scaf-c", SequenceClass::Scaffold, 6),
        ];
        let mut local_rng = rng_from_number(42);
        let plan = select_reduced_reference(
            &reference_id("ref-source"),
            &sequences,
            14,
            5,
            &mut local_rng,
        )
        .unwrap();
        let observed = plan
            .segments_in_selection_order()
            .iter()
            .map(|segment| {
                (
                    segment.selection_order(),
                    segment.class(),
                    segment.source_interval().sequence_id().as_str().to_owned(),
                    segment.source_interval().start(),
                    segment.source_interval().end(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            observed,
            vec![
                (1, ReductionSegmentClass::Full, "chr-a".to_owned(), 0, 5),
                (2, ReductionSegmentClass::Full, "chr-b".to_owned(), 7, 12),
                (3, ReductionSegmentClass::Residual, "chr-b".to_owned(), 3, 7),
            ]
        );
    }
    #[test]
    fn direct_selector_is_deterministic_non_overlapping_and_has_one_residual() {
        let sequences = vec![
            sequence("chr1", SequenceClass::Chromosome, 31),
            sequence("chr2", SequenceClass::Chromosome, 17),
            sequence("scaf", SequenceClass::Scaffold, 4),
        ];
        let mut first_rng = rng(41);
        let first =
            select_reduced_reference(&reference_id("source"), &sequences, 29, 5, &mut first_rng)
                .unwrap();
        let mut second_rng = rng(41);
        assert_eq!(
            first,
            select_reduced_reference(&reference_id("source"), &sequences, 29, 5, &mut second_rng)
                .unwrap()
        );
        assert_eq!(
            first
                .segments
                .iter()
                .map(|segment| segment.source_interval.len())
                .sum::<u64>(),
            29
        );
        assert_eq!(
            first
                .segments
                .iter()
                .filter(|segment| segment.class == ReductionSegmentClass::Residual)
                .count(),
            1
        );
        for (index, segment) in first.segments.iter().enumerate() {
            for other in &first.segments[index + 1..] {
                if segment.source_interval.sequence_id() == other.source_interval.sequence_id() {
                    assert!(
                        segment.source_interval.end() <= other.source_interval.start()
                            || other.source_interval.end() <= segment.source_interval.start()
                    );
                }
            }
        }
    }
    #[test]
    fn exact_multiple_has_no_residual() {
        let sequences = vec![sequence("chr1", SequenceClass::Chromosome, 20)];
        let mut local_rng = rng(19);
        let plan =
            select_reduced_reference(&reference_id("source"), &sequences, 10, 5, &mut local_rng)
                .unwrap();
        assert!(
            plan.segments
                .iter()
                .all(|segment| segment.class == ReductionSegmentClass::Full)
        );
    }

    #[test]
    fn direct_layout_handles_hundred_megabases_with_ten_blocks_without_per_base_work() {
        let sequences = vec![sequence("chr1", SequenceClass::Chromosome, 100_000_000)];
        let mut local_rng = rng(3);
        let plan = select_reduced_reference(
            &reference_id("source"),
            &sequences,
            10_000_000,
            1_000_000,
            &mut local_rng,
        )
        .unwrap();
        assert_eq!(plan.segments_in_selection_order().len(), 10);
    }
    #[test]
    fn boundaries_are_rejected_without_consuming_rng() {
        let sequences = vec![
            sequence("chr1", SequenceClass::Chromosome, 3),
            sequence("chr2", SequenceClass::Chromosome, 3),
        ];
        let mut after = rng(9);
        assert_eq!(
            select_reduced_reference(&reference_id("source"), &sequences, 5, 4, &mut after),
            Err(ReductionError::InfeasibleLayout)
        );
        let mut untouched = rng(9);
        assert_eq!(
            select_reduced_reference(&reference_id("source"), &sequences, 2, 2, &mut after)
                .unwrap(),
            select_reduced_reference(&reference_id("source"), &sequences, 2, 2, &mut untouched)
                .unwrap()
        );
    }
}
