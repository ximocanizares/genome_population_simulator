//! Quimeras transitorias de parejas F4.
//!
//! Esta transformación ocurre entre el fragmento biológico y la orientación
//! *paired-end*. Conserva la procedencia de los dos segmentos únicamente en su
//! resultado temporal; el llamador puede emitir la pareja y retener solo el
//! contador agregado.

use std::num::NonZeroU64;

use crate::haplotype::{HaplotypeCopy, HaplotypeCopyId, Individual};
use crate::identity::EventId;
use crate::materialization::{MaterializationError, materialize_copy_interval};
use crate::paired_end::{PairedEndFragment, PairedEndPair, emit_paired_end};
use crate::randomness::{
    F4TechnicalDecision, RandomSamplingError, RootSeed, derive_f4_technical_subseed, f4_bernoulli,
    f4_uniform_index,
};
use crate::reference::Reference;

/// Frecuencia efectiva de quimeras para una biblioteca. El perfil inicial usa
/// el neutro `0`; una fase de configuración posterior puede construir este
/// tipo con una tasa distinta sin cambiar el algoritmo.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChimeraSettings {
    rate: f64,
}

impl Default for ChimeraSettings {
    fn default() -> Self {
        Self { rate: 0.0 }
    }
}

impl ChimeraSettings {
    pub fn new(rate: f64) -> Result<Self, ChimeraError> {
        if !rate.is_finite() || !(0.0..=1.0).contains(&rate) {
            return Err(ChimeraError::InvalidRate);
        }
        Ok(Self { rate })
    }

    #[must_use]
    pub const fn rate(self) -> f64 {
        self.rate
    }
}

/// Un intervalo de origen, en orientación directa de su copia.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ChimeraSegment {
    copy_id: HaplotypeCopyId,
    start: u64,
    end: u64,
}

impl ChimeraSegment {
    #[must_use]
    pub const fn copy_id(&self) -> &HaplotypeCopyId {
        &self.copy_id
    }

    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn end(&self) -> u64 {
        self.end
    }
}

/// Procedencia efímera de una quimera. No pertenece a FASTQ ni al manifiesto.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ChimeraOrigin {
    first: ChimeraSegment,
    second: ChimeraSegment,
    breakpoint: NonZeroU64,
}

impl ChimeraOrigin {
    #[must_use]
    pub const fn first(&self) -> &ChimeraSegment {
        &self.first
    }

    #[must_use]
    pub const fn second(&self) -> &ChimeraSegment {
        &self.second
    }

    #[must_use]
    pub const fn breakpoint(&self) -> NonZeroU64 {
        self.breakpoint
    }
}

/// Resultado por pareja: la procedencia, si existe, sigue siendo transitoria.
#[derive(Clone, Debug, PartialEq)]
pub struct ChimeraTransformation {
    pair: PairedEndPair,
    origin: Option<ChimeraOrigin>,
}

impl ChimeraTransformation {
    #[must_use]
    pub const fn pair(&self) -> &PairedEndPair {
        &self.pair
    }

    #[must_use]
    pub const fn origin(&self) -> Option<&ChimeraOrigin> {
        self.origin.as_ref()
    }

    #[must_use]
    pub const fn transformed(&self) -> bool {
        self.origin.is_some()
    }

    #[must_use]
    pub fn into_pair(self) -> PairedEndPair {
        self.pair
    }
}

/// Recuento agregado, apto para el resumen de un evento.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ChimeraReport {
    realized_pairs: u64,
}

impl ChimeraReport {
    #[must_use]
    pub const fn realized_pairs(self) -> u64 {
        self.realized_pairs
    }
}

/// Acumula solamente el recuento realizado; nunca guarda parejas ni orígenes.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ChimeraAccumulator(ChimeraReport);

impl ChimeraAccumulator {
    pub fn record(&mut self, transformation: &ChimeraTransformation) -> Result<(), ChimeraError> {
        if transformation.transformed() {
            self.0.realized_pairs = self
                .0
                .realized_pairs
                .checked_add(1)
                .ok_or(ChimeraError::CounterOverflow)?;
        }
        Ok(())
    }

    #[must_use]
    pub const fn report(self) -> ChimeraReport {
        self.0
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ChimeraError {
    InvalidRate,
    PairIndividualMismatch,
    PairSourceMissing,
    PairSourceMismatch,
    NoSecondSource,
    /// El índice sorteado se sale de las candidatas que se acaban de contar:
    /// una invariante interna, distinta de que no haya candidatas (DT-031).
    SelectedSourceOutOfRange,
    LengthOverflow,
    CounterOverflow,
    Random(RandomSamplingError),
    Materialization(MaterializationError),
    Pair(crate::diagnostics::Diagnostic),
}

/// Aplica el Bernoulli F4 a una pareja limpia y, si procede, compone el
/// inserto quimérico antes de emitir de nuevo su orientación PE normal.
pub fn transform_paired_end(
    root: RootSeed,
    event_id: &EventId,
    reference: &Reference,
    individual: &Individual,
    pair: PairedEndPair,
    settings: ChimeraSettings,
) -> Result<ChimeraTransformation, ChimeraError> {
    let template = pair.fragment().template();
    let ordinal = template.template_ordinal().get();
    let library_id = template.library_id();
    let mut selected_rng = derive_f4_technical_subseed(
        root,
        F4TechnicalDecision::chimera_bernoulli(event_id, library_id, ordinal),
    )
    .map_err(ChimeraError::Pair)?
    .local_rng();
    if !f4_bernoulli(&mut selected_rng, settings.rate).map_err(ChimeraError::Random)? {
        return Ok(ChimeraTransformation { pair, origin: None });
    }

    let fragment = pair.fragment();
    if fragment.haplotype_copy_id().individual_id() != individual.id() {
        return Err(ChimeraError::PairIndividualMismatch);
    }
    let insert_length =
        u64::try_from(fragment.sequence().len()).map_err(|_| ChimeraError::LengthOverflow)?;
    if insert_length < 2 {
        return Err(ChimeraError::NoSecondSource);
    }
    let mut breakpoint_rng = derive_f4_technical_subseed(
        root,
        F4TechnicalDecision::chimera_breakpoint(event_id, library_id, ordinal),
    )
    .map_err(ChimeraError::Pair)?
    .local_rng();
    let breakpoint = f4_uniform_index(&mut breakpoint_rng, insert_length - 1)
        .map_err(ChimeraError::Random)?
        .checked_add(1)
        .and_then(NonZeroU64::new)
        .ok_or(ChimeraError::LengthOverflow)?;
    let suffix_length = insert_length
        .checked_sub(breakpoint.get())
        .ok_or(ChimeraError::LengthOverflow)?;

    let first_copy = find_copy(individual, fragment.haplotype_copy_id())
        .ok_or(ChimeraError::PairSourceMissing)?;
    let first_end = fragment
        .start()
        .checked_add(breakpoint.get())
        .ok_or(ChimeraError::LengthOverflow)?;
    let first_sequence =
        materialize_copy_interval(first_copy, reference, fragment.start(), first_end)
            .map_err(ChimeraError::Materialization)?;
    if first_sequence.as_slice()
        != &fragment.sequence()
            [..usize::try_from(breakpoint.get()).map_err(|_| ChimeraError::LengthOverflow)?]
    {
        return Err(ChimeraError::PairSourceMismatch);
    }

    let candidates =
        canonical_second_source_ranges(individual, suffix_length, fragment, first_end)?;
    let candidate_count = candidates.iter().try_fold(0_u64, |total, candidate| {
        total
            .checked_add(candidate.count)
            .ok_or(ChimeraError::LengthOverflow)
    })?;
    if candidate_count == 0 {
        return Err(ChimeraError::NoSecondSource);
    }
    let mut source_rng = derive_f4_technical_subseed(
        root,
        F4TechnicalDecision::chimera_second_source(event_id, library_id, ordinal, breakpoint),
    )
    .map_err(ChimeraError::Pair)?
    .local_rng();
    let selected =
        f4_uniform_index(&mut source_rng, candidate_count).map_err(ChimeraError::Random)?;
    let (second_copy, second_start) = select_second_source(&candidates, selected)?;
    let second_end = second_start
        .checked_add(suffix_length)
        .ok_or(ChimeraError::LengthOverflow)?;
    let second_sequence =
        materialize_copy_interval(second_copy, reference, second_start, second_end)
            .map_err(ChimeraError::Materialization)?;

    let mut sequence = first_sequence;
    sequence.extend_from_slice(&second_sequence);
    let fragment = PairedEndFragment::new(
        template.clone(),
        fragment.haplotype_copy_id().clone(),
        fragment.start(),
        sequence,
    )
    .map_err(ChimeraError::Pair)?;
    let transformed_pair = emit_paired_end(
        fragment,
        pair.pre_hop_identity().clone(),
        pair.r1().sequence().len() as u16,
    )
    .map_err(ChimeraError::Pair)?;
    Ok(ChimeraTransformation {
        pair: transformed_pair,
        origin: Some(ChimeraOrigin {
            first: ChimeraSegment {
                copy_id: first_copy.id().clone(),
                start: pair.fragment().start(),
                end: first_end,
            },
            second: ChimeraSegment {
                copy_id: second_copy.id().clone(),
                start: second_start,
                end: second_end,
            },
            breakpoint,
        }),
    })
}

struct CandidateRange<'a> {
    copy: &'a HaplotypeCopy,
    start: u64,
    count: u64,
}

fn canonical_second_source_ranges<'a>(
    individual: &'a Individual,
    suffix_length: u64,
    fragment: &PairedEndFragment,
    first_end: u64,
) -> Result<Vec<CandidateRange<'a>>, ChimeraError> {
    let mut copies = individual
        .sequences()
        .iter()
        .flat_map(|genotype| genotype.copies().iter())
        .collect::<Vec<_>>();
    copies.sort_by(|left, right| left.id().cmp(right.id()));

    let mut ranges = Vec::with_capacity(copies.len() * 2);
    for copy in copies {
        let length = copy
            .segments()
            .current_length()
            .map_err(|_| ChimeraError::LengthOverflow)?;
        let Some(maximum_start) = length.checked_sub(suffix_length) else {
            continue;
        };
        if copy.id() != fragment.haplotype_copy_id() {
            push_range(&mut ranges, copy, 0, maximum_start)?;
            continue;
        }

        // Para intervalos semiabiertos, solapar o tocar el prefijo equivale a
        // `start <= first_end && start + suffix_length >= first_start`.
        // Los dos rangos restantes preservan el orden creciente de `start`.
        if let Some(left_last) = fragment.start().checked_sub(suffix_length + 1) {
            push_range(&mut ranges, copy, 0, left_last.min(maximum_start))?;
        }
        if let Some(right_first) = first_end.checked_add(1)
            && right_first <= maximum_start
        {
            push_range(&mut ranges, copy, right_first, maximum_start)?;
        }
    }
    Ok(ranges)
}

fn push_range<'a>(
    ranges: &mut Vec<CandidateRange<'a>>,
    copy: &'a HaplotypeCopy,
    first: u64,
    last: u64,
) -> Result<(), ChimeraError> {
    let count = last
        .checked_sub(first)
        .and_then(|value| value.checked_add(1))
        .ok_or(ChimeraError::LengthOverflow)?;
    ranges.push(CandidateRange {
        copy,
        start: first,
        count,
    });
    Ok(())
}

fn select_second_source<'a>(
    candidates: &'a [CandidateRange<'a>],
    selected: u64,
) -> Result<(&'a HaplotypeCopy, u64), ChimeraError> {
    let mut cursor = 0_u64;
    for candidate in candidates {
        let next = cursor
            .checked_add(candidate.count)
            .ok_or(ChimeraError::LengthOverflow)?;
        if selected < next {
            let offset = selected
                .checked_sub(cursor)
                .ok_or(ChimeraError::LengthOverflow)?;
            let start = candidate
                .start
                .checked_add(offset)
                .ok_or(ChimeraError::LengthOverflow)?;
            return Ok((candidate.copy, start));
        }
        cursor = next;
    }
    Err(ChimeraError::SelectedSourceOutOfRange)
}

fn find_copy<'a>(
    individual: &'a Individual,
    copy_id: &HaplotypeCopyId,
) -> Option<&'a HaplotypeCopy> {
    individual
        .sequences()
        .iter()
        .flat_map(|genotype| genotype.copies())
        .find(|copy| copy.id() == copy_id)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// DT-031: un índice sorteado fuera de las candidatas contadas es una
    /// invariante interna, y no debe confundirse con que no haya candidatas
    /// (`NoSecondSource`), que es lo que el ejecutor presenta como fallo
    /// previsto. Con la lista vacía no hay ningún rango que cubra el índice.
    #[test]
    fn a_selected_index_outside_the_counted_candidates_is_not_a_missing_source() {
        assert!(matches!(
            select_second_source(&[], 0),
            Err(ChimeraError::SelectedSourceOutOfRange)
        ));
    }
}
