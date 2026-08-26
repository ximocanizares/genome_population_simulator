//! Composición transaccional de SNP e indels pequeños en orden declarado.
//!
//! Esta frontera no interpreta TOML ni crea individuos. Recibe una copia ya
//! identificada y los eventos que le corresponden; cada evento observa el
//! estado dejado por el anterior. El preflight recorre exactamente esa misma
//! secuencia sobre estados persistentes, de modo que ningún prefijo se publica
//! si un conflicto o una entrada incompatible aparece después.

use std::collections::BTreeSet;
use std::fmt;

use crate::haplotype::Individual;
use crate::identity::{ScenarioId, SequenceId};
use crate::indel_application::{
    IndelApplicationError, RandomSmallIndelRequest, SmallIndelApplicationResult, SmallIndelLengths,
    apply_random_small_indels, apply_sampled_small_indels,
};
use crate::randomness::RootSeed;
use crate::reference::Reference;
use crate::segments::{ModificationProvenance, SegmentError};
use crate::snp_application::{
    SnpApplicationError, SnpOrigin, TargetedSnp, apply_preflighted_snps_to_owned, apply_snps,
    preflight_snps,
};

/// Petición de indels pequeños que posee sus identificadores y conserva el
/// contexto de RNG de F2.2.9 sin depender de la vida de la configuración.
#[derive(Clone, Debug)]
pub struct DeclaredSmallIndelRequest {
    scenario_id: ScenarioId,
    sequence_id: SequenceId,
    phase: crate::haplotype::HaplotypePhase,
    spatial_scope_id: String,
    insertion_rate_per_opportunity: f64,
    deletion_rate_per_opportunity: f64,
}

impl DeclaredSmallIndelRequest {
    #[must_use]
    pub fn new(
        scenario_id: ScenarioId,
        sequence_id: SequenceId,
        phase: crate::haplotype::HaplotypePhase,
        spatial_scope_id: String,
        insertion_rate_per_opportunity: f64,
        deletion_rate_per_opportunity: f64,
    ) -> Self {
        Self {
            scenario_id,
            sequence_id,
            phase,
            spatial_scope_id,
            insertion_rate_per_opportunity,
            deletion_rate_per_opportunity,
        }
    }

    fn as_request(&self) -> RandomSmallIndelRequest<'_> {
        RandomSmallIndelRequest {
            scenario_id: &self.scenario_id,
            sequence_id: &self.sequence_id,
            phase: self.phase,
            spatial_scope_id: &self.spatial_scope_id,
            insertion_rate_per_opportunity: self.insertion_rate_per_opportunity,
            deletion_rate_per_opportunity: self.deletion_rate_per_opportunity,
        }
    }
}

/// Recuento de indels que ya ha sido sorteado para la frontera de pruebas.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SampledSmallIndelCounts {
    insertion_count: u64,
    deletion_count: u64,
}

impl SampledSmallIndelCounts {
    #[must_use]
    pub const fn new(insertion_count: u64, deletion_count: u64) -> Self {
        Self {
            insertion_count,
            deletion_count,
        }
    }
}

/// Inserción explícita en la coordenada de la secuencia vigente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TargetedCurrentInsertion {
    sequence_id: SequenceId,
    phase: crate::haplotype::HaplotypePhase,
    point: u64,
    bases: std::sync::Arc<[u8]>,
    provenance: ModificationProvenance,
}

impl TargetedCurrentInsertion {
    #[must_use]
    pub fn new(
        sequence_id: SequenceId,
        phase: crate::haplotype::HaplotypePhase,
        point: u64,
        bases: impl Into<std::sync::Arc<[u8]>>,
        provenance: ModificationProvenance,
    ) -> Self {
        Self {
            sequence_id,
            phase,
            point,
            bases: bases.into(),
            provenance,
        }
    }

    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub const fn phase(&self) -> crate::haplotype::HaplotypePhase {
        self.phase
    }

    #[must_use]
    pub const fn point(&self) -> u64 {
        self.point
    }

    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }
}

/// Deleción explícita en las coordenadas de la secuencia vigente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TargetedCurrentDeletion {
    sequence_id: SequenceId,
    phase: crate::haplotype::HaplotypePhase,
    start: u64,
    end: u64,
    provenance: ModificationProvenance,
}

impl TargetedCurrentDeletion {
    #[must_use]
    pub const fn new(
        sequence_id: SequenceId,
        phase: crate::haplotype::HaplotypePhase,
        start: u64,
        end: u64,
        provenance: ModificationProvenance,
    ) -> Self {
        Self {
            sequence_id,
            phase,
            start,
            end,
            provenance,
        }
    }

    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub const fn phase(&self) -> crate::haplotype::HaplotypePhase {
        self.phase
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
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }
}

/// Evento explícito, aplicado estrictamente en el orden declarado.
#[derive(Clone, Debug)]
pub enum DeclaredSmallEvent {
    Snp(TargetedSnp),
    Insertion(TargetedCurrentInsertion),
    Deletion(TargetedCurrentDeletion),
}

/// Bloque aleatorio posterior a todos los eventos explícitos.
#[derive(Clone, Debug)]
pub enum RandomSmallEventBlock {
    RandomIndels {
        lengths: SmallIndelLengths,
        request: DeclaredSmallIndelRequest,
    },
    SampledIndels {
        lengths: SmallIndelLengths,
        request: DeclaredSmallIndelRequest,
        counts: SampledSmallIndelCounts,
    },
}

/// Clase que identifica el elemento declarado que falló, sin depender del
/// texto de sus errores internos.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DeclaredSmallEventKind {
    Snp,
    Insertion,
    Deletion,
    RandomIndels,
    SampledIndels,
}

impl DeclaredSmallEvent {
    #[must_use]
    pub const fn kind(&self) -> DeclaredSmallEventKind {
        match self {
            Self::Snp(_) => DeclaredSmallEventKind::Snp,
            Self::Insertion(_) => DeclaredSmallEventKind::Insertion,
            Self::Deletion(_) => DeclaredSmallEventKind::Deletion,
        }
    }
}

impl RandomSmallEventBlock {
    #[must_use]
    pub const fn kind(&self) -> DeclaredSmallEventKind {
        match self {
            Self::RandomIndels { .. } => DeclaredSmallEventKind::RandomIndels,
            Self::SampledIndels { .. } => DeclaredSmallEventKind::SampledIndels,
        }
    }
}

/// Resumen de cada evento realmente aplicado al estado final.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AppliedSmallEvent {
    Snp,
    Insertion,
    Deletion,
    Indels {
        insertion_sampled: u64,
        insertion_accepted: u64,
        deletion_sampled: u64,
        deletion_accepted: u64,
    },
}

/// Estado final y trazabilidad mínima de la composición declarada.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SmallEventCompositionResult {
    individual: Individual,
    applied_events: Vec<AppliedSmallEvent>,
}

impl SmallEventCompositionResult {
    #[must_use]
    pub const fn individual(&self) -> &Individual {
        &self.individual
    }

    #[must_use]
    pub fn applied_events(&self) -> &[AppliedSmallEvent] {
        &self.applied_events
    }
}

/// Rechazo contextualizado de un elemento de la lista declarada.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SmallEventCompositionError {
    Event {
        declared_index: usize,
        kind: DeclaredSmallEventKind,
        source: SmallEventSourceError,
    },
}

/// Causa tipada conservada por el diagnóstico de composición.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SmallEventSourceError {
    Snp(SnpApplicationError),
    Indel(IndelApplicationError),
    Segment(SegmentError),
    IndividualReferenceMismatch,
    MissingTargetCopy {
        phase: crate::haplotype::HaplotypePhase,
    },
}

impl fmt::Display for SmallEventCompositionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for SmallEventCompositionError {}

/// Valida toda la composición sobre estados candidatos y no cambia el
/// individuo aportado. La aplicación usa el mismo recorrido para preservar la
/// atomicidad incluso si el llamador omite este preflight explícito.
pub fn preflight_declared_small_events(
    individual: &Individual,
    reference: &Reference,
    root_seed: RootSeed,
    explicit_events: &[DeclaredSmallEvent],
    random_blocks: &[RandomSmallEventBlock],
) -> Result<(), SmallEventCompositionError> {
    compose_declared_small_events(
        individual,
        reference,
        root_seed,
        explicit_events,
        random_blocks,
    )
    .map(|_| ())
}

/// Aplica SNP e indels pequeños en el orden declarado y devuelve únicamente
/// el resultado completo. Los indels de un mismo elemento mantienen el orden
/// local inserción→deleción de F2.2.9.
pub fn apply_declared_small_events(
    individual: &Individual,
    reference: &Reference,
    root_seed: RootSeed,
    explicit_events: &[DeclaredSmallEvent],
    random_blocks: &[RandomSmallEventBlock],
) -> Result<SmallEventCompositionResult, SmallEventCompositionError> {
    compose_declared_small_events(
        individual,
        reference,
        root_seed,
        explicit_events,
        random_blocks,
    )
}

fn compose_declared_small_events(
    individual: &Individual,
    reference: &Reference,
    root_seed: RootSeed,
    explicit_events: &[DeclaredSmallEvent],
    random_blocks: &[RandomSmallEventBlock],
) -> Result<SmallEventCompositionResult, SmallEventCompositionError> {
    if individual.reference_id() != reference.id()
        && let Some(kind) = explicit_events
            .first()
            .map(DeclaredSmallEvent::kind)
            .or_else(|| random_blocks.first().map(RandomSmallEventBlock::kind))
    {
        return Err(SmallEventCompositionError::Event {
            declared_index: 0,
            kind,
            source: SmallEventSourceError::IndividualReferenceMismatch,
        });
    }
    let mut updated = individual.clone();
    let mut applied_events = Vec::with_capacity(explicit_events.len() + random_blocks.len());
    let mut declared_index = 0;
    while declared_index < explicit_events.len() {
        if let Some(batch_end) = consecutive_matrix_snp_batch_end(explicit_events, declared_index) {
            let batch = &explicit_events[declared_index..batch_end];
            for (offset, event) in batch.iter().enumerate() {
                let DeclaredSmallEvent::Snp(snp) = event else {
                    unreachable!("a matrix SNP batch contains only SNP events");
                };
                preflight_snps(&updated, reference, std::slice::from_ref(snp)).map_err(
                    |source| SmallEventCompositionError::Event {
                        declared_index: declared_index + offset,
                        kind: event.kind(),
                        source: SmallEventSourceError::Snp(source),
                    },
                )?;
            }
            updated = apply_preflighted_snps_to_owned(
                updated,
                batch.iter().filter_map(|event| match event {
                    DeclaredSmallEvent::Snp(snp) => Some(snp),
                    DeclaredSmallEvent::Insertion(_) | DeclaredSmallEvent::Deletion(_) => None,
                }),
            )
            .map_err(|source| SmallEventCompositionError::Event {
                declared_index,
                kind: DeclaredSmallEventKind::Snp,
                source: SmallEventSourceError::Snp(source),
            })?;
            applied_events.extend(std::iter::repeat_n(AppliedSmallEvent::Snp, batch.len()));
            declared_index = batch_end;
            continue;
        }

        let event = &explicit_events[declared_index];
        let result = match event {
            DeclaredSmallEvent::Snp(snp) => {
                apply_snps(&updated, reference, std::slice::from_ref(snp))
                    .map(|individual| (individual, AppliedSmallEvent::Snp))
                    .map_err(SmallEventSourceError::Snp)
            }
            DeclaredSmallEvent::Insertion(insertion) => {
                apply_current_insertion(&updated, insertion)
                    .map(|individual| (individual, AppliedSmallEvent::Insertion))
            }
            DeclaredSmallEvent::Deletion(deletion) => apply_current_deletion(&updated, deletion)
                .map(|individual| (individual, AppliedSmallEvent::Deletion)),
        };
        let (next, applied) = result.map_err(|source| SmallEventCompositionError::Event {
            declared_index,
            kind: event.kind(),
            source,
        })?;
        updated = next;
        applied_events.push(applied);
        declared_index += 1;
    }
    for (block_index, block) in random_blocks.iter().enumerate() {
        let declared_index = explicit_events.len() + block_index;
        let result = match block {
            RandomSmallEventBlock::RandomIndels { lengths, request } => apply_random_small_indels(
                &updated,
                reference,
                root_seed,
                *lengths,
                &request.as_request(),
            )
            .map(indel_result)
            .map_err(SmallEventSourceError::Indel),
            RandomSmallEventBlock::SampledIndels {
                lengths,
                request,
                counts,
            } => apply_sampled_small_indels(
                &updated,
                reference,
                root_seed,
                *lengths,
                &request.as_request(),
                counts.insertion_count,
                counts.deletion_count,
            )
            .map(indel_result)
            .map_err(SmallEventSourceError::Indel),
        };
        let (next, applied) = result.map_err(|source| SmallEventCompositionError::Event {
            declared_index,
            kind: block.kind(),
            source,
        })?;
        updated = next;
        applied_events.push(applied);
    }
    Ok(SmallEventCompositionResult {
        individual: updated,
        applied_events,
    })
}

/// Devuelve el final exclusivo de un lote de SNP de matriz contiguos. Los
/// eventos explícitos y los indels nunca se atraviesan: conservan el orden
/// declarado y ven exactamente el estado dejado por el lote anterior.
fn consecutive_matrix_snp_batch_end(events: &[DeclaredSmallEvent], start: usize) -> Option<usize> {
    let DeclaredSmallEvent::Snp(first) = events.get(start)? else {
        return None;
    };
    if !matches!(first.snp().origin(), SnpOrigin::PopulationMatrix { .. }) {
        return None;
    }
    let end = events[start..]
        .iter()
        .take_while(|event| {
            matches!(
                event,
                DeclaredSmallEvent::Snp(snp)
                    if matches!(snp.snp().origin(), SnpOrigin::PopulationMatrix { .. })
            )
        })
        .count()
        + start;
    if end <= start + 1 {
        return None;
    }
    let mut targets = BTreeSet::new();
    events[start..end]
        .iter()
        .all(|event| {
            let DeclaredSmallEvent::Snp(snp) = event else {
                return false;
            };
            let position = snp.snp().locus().position();
            targets.insert((position.sequence_id(), snp.phase(), position.offset()))
        })
        .then_some(end)
}

fn apply_current_insertion(
    individual: &Individual,
    insertion: &TargetedCurrentInsertion,
) -> Result<Individual, SmallEventSourceError> {
    let copy = individual
        .copies_for(&insertion.sequence_id)
        .and_then(|copies| copies.iter().find(|copy| copy.phase() == insertion.phase))
        .ok_or(SmallEventSourceError::MissingTargetCopy {
            phase: insertion.phase,
        })?;
    let segments = copy
        .segments()
        .with_current_insertion(
            insertion.point,
            insertion.bases.clone(),
            insertion.provenance.clone(),
        )
        .map_err(SmallEventSourceError::Segment)?;
    individual
        .with_copy_segments(&insertion.sequence_id, insertion.phase, segments)
        .map_err(|_| SmallEventSourceError::MissingTargetCopy {
            phase: insertion.phase,
        })
}

fn apply_current_deletion(
    individual: &Individual,
    deletion: &TargetedCurrentDeletion,
) -> Result<Individual, SmallEventSourceError> {
    let copy = individual
        .copies_for(&deletion.sequence_id)
        .and_then(|copies| copies.iter().find(|copy| copy.phase() == deletion.phase))
        .ok_or(SmallEventSourceError::MissingTargetCopy {
            phase: deletion.phase,
        })?;
    let segments = copy
        .segments()
        .with_current_deletion(deletion.start, deletion.end, deletion.provenance.clone())
        .map_err(SmallEventSourceError::Segment)?;
    individual
        .with_copy_segments(&deletion.sequence_id, deletion.phase, segments)
        .map_err(|_| SmallEventSourceError::MissingTargetCopy {
            phase: deletion.phase,
        })
}

fn indel_result(result: SmallIndelApplicationResult) -> (Individual, AppliedSmallEvent) {
    let applied = AppliedSmallEvent::Indels {
        insertion_sampled: result.insertion().sampled_count(),
        insertion_accepted: result.insertion().accepted_count(),
        deletion_sampled: result.deletion().sampled_count(),
        deletion_accepted: result.deletion().accepted_count(),
    };
    (result.individual().clone(), applied)
}
