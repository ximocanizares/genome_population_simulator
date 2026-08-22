//! Aplicación dinámica y transaccional de indels pequeños aleatorios.
//!
//! F2.2.9 opera una clase cada vez: las inserciones se completan antes que las
//! deleciones y, dentro de cada clase, los candidatos se reenumeran después de
//! cada aceptación. La composición con SNP y la resolución global permanecen
//! en F2.2.10.

use std::fmt;
use std::sync::Arc;

use crate::haplotype::{HaplotypePhase, Individual};
use crate::identity::{
    AutomaticIdRequest, EventId, IdClass, IdentityRegistry, ScenarioId, SequenceId,
};
use crate::materialization::{MaterializationError, materialize_copy};
use crate::randomness::{
    LocalRng, RandomSamplingError, RootSeed, RuntimeSubseedScope, categorical_f64, derive_subseed,
    poisson_count, uniform_u64_below,
};
use crate::reference::Reference;
use crate::segments::{ModificationProvenance, SegmentError};

pub const DEFAULT_STRUCTURAL_INDEL_THRESHOLD: u64 = 50;
pub const DYNAMIC_INDEL_SELECTION_ALGORITHM: &str = "dynamic-indel-selection-v1";

/// Longitudes exactas de las dos clases de indels aleatorios.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SmallIndelLengths {
    insertion_length: u64,
    deletion_length: u64,
    effective_structural_indel_threshold: u64,
}

impl Default for SmallIndelLengths {
    fn default() -> Self {
        // 50 y 1 pertenecen al dominio publicado y no requieren validación.
        Self {
            insertion_length: 1,
            deletion_length: 1,
            effective_structural_indel_threshold: DEFAULT_STRUCTURAL_INDEL_THRESHOLD,
        }
    }
}

impl SmallIndelLengths {
    pub fn new(
        insertion_length: u64,
        deletion_length: u64,
        effective_structural_indel_threshold: u64,
    ) -> Result<Self, IndelApplicationError> {
        let value = Self {
            insertion_length,
            deletion_length,
            effective_structural_indel_threshold,
        };
        value.validate_length(insertion_length, IndelClass::Insertion)?;
        value.validate_length(deletion_length, IndelClass::Deletion)?;
        Ok(value)
    }

    #[must_use]
    pub const fn insertion_length(self) -> u64 {
        self.insertion_length
    }

    #[must_use]
    pub const fn deletion_length(self) -> u64 {
        self.deletion_length
    }

    #[must_use]
    pub const fn effective_structural_indel_threshold(self) -> u64 {
        self.effective_structural_indel_threshold
    }

    fn validate_length(self, length: u64, class: IndelClass) -> Result<(), IndelApplicationError> {
        if length == 0 || length >= self.effective_structural_indel_threshold {
            return Err(IndelApplicationError::InvalidLength {
                class,
                length,
                effective_structural_indel_threshold: self.effective_structural_indel_threshold,
            });
        }
        Ok(())
    }

    fn length_for(self, class: IndelClass) -> u64 {
        match class {
            IndelClass::Insertion => self.insertion_length,
            IndelClass::Deletion => self.deletion_length,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IndelClass {
    Insertion,
    Deletion,
}

impl IndelClass {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Insertion => "insertion",
            Self::Deletion => "deletion",
        }
    }
}

/// Contexto estable de los tres ámbitos RNG de una clase y copia.
#[derive(Clone, Debug)]
pub struct RandomIndelRequest<'a> {
    pub scenario_id: &'a ScenarioId,
    pub sequence_id: &'a SequenceId,
    pub phase: HaplotypePhase,
    pub spatial_scope_id: &'a str,
    pub class: IndelClass,
    pub rate_per_opportunity: f64,
}

/// Petición conjunta de F2.2.9. El orden es siempre inserciones y deleciones.
#[derive(Clone, Debug)]
pub struct RandomSmallIndelRequest<'a> {
    pub scenario_id: &'a ScenarioId,
    pub sequence_id: &'a SequenceId,
    pub phase: HaplotypePhase,
    pub spatial_scope_id: &'a str,
    pub insertion_rate_per_opportunity: f64,
    pub deletion_rate_per_opportunity: f64,
}

/// Resultado de la aplicación conjunta, con las dos clases separadas.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SmallIndelApplicationResult {
    individual: Individual,
    insertion: IndelApplicationResult,
    deletion: IndelApplicationResult,
}

impl SmallIndelApplicationResult {
    #[must_use]
    pub const fn individual(&self) -> &Individual {
        &self.individual
    }

    #[must_use]
    pub const fn insertion(&self) -> &IndelApplicationResult {
        &self.insertion
    }

    #[must_use]
    pub const fn deletion(&self) -> &IndelApplicationResult {
        &self.deletion
    }
}

/// Resultado de una clase aleatoria publicada de forma atómica.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IndelApplicationResult {
    individual: Individual,
    initial_opportunities: u64,
    sampled_count: u64,
    accepted_count: u64,
}

impl IndelApplicationResult {
    #[must_use]
    pub const fn individual(&self) -> &Individual {
        &self.individual
    }

    #[must_use]
    pub const fn initial_opportunities(&self) -> u64 {
        self.initial_opportunities
    }

    #[must_use]
    pub const fn sampled_count(&self) -> u64 {
        self.sampled_count
    }

    #[must_use]
    pub const fn accepted_count(&self) -> u64 {
        self.accepted_count
    }
}

/// Rechazo estable de configuración, selección o aplicación de una clase.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum IndelApplicationError {
    InvalidLength {
        class: IndelClass,
        length: u64,
        effective_structural_indel_threshold: u64,
    },
    InvalidRate,
    IndividualReferenceMismatch,
    MissingTargetCopy {
        phase: HaplotypePhase,
    },
    Materialization(MaterializationError),
    Segment(SegmentError),
    Random(RandomSamplingError),
    Randomness {
        detail: String,
    },
    ArithmeticOverflow,
    ExhaustedOpportunities {
        sampled_count: u64,
        accepted_before_failure: u64,
        active_before_rollback: u64,
        opportunities_at_failure: u64,
    },
}

impl fmt::Display for IndelApplicationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for IndelApplicationError {}

/// Cuenta una vez con Poisson y aplica una clase de indels sin publicar prefijos.
pub fn apply_random_indels(
    individual: &Individual,
    reference: &Reference,
    root_seed: RootSeed,
    lengths: SmallIndelLengths,
    request: &RandomIndelRequest<'_>,
) -> Result<IndelApplicationResult, IndelApplicationError> {
    if individual.reference_id() != reference.id() {
        return Err(IndelApplicationError::IndividualReferenceMismatch);
    }
    if !request.rate_per_opportunity.is_finite() || request.rate_per_opportunity < 0.0 {
        return Err(IndelApplicationError::InvalidRate);
    }
    lengths.validate_length(lengths.length_for(request.class), request.class)?;
    let initial_opportunities = enumerate_opportunities(
        current_bases(individual, reference, request.sequence_id, request.phase)?,
        request.class,
        lengths.length_for(request.class),
    )
    .len()
    .try_into()
    .map_err(|_| IndelApplicationError::ArithmeticOverflow)?;

    // Una tasa cero no deriva ni avanza ningún ámbito aleatorio.
    if request.rate_per_opportunity == 0.0 {
        return Ok(IndelApplicationResult {
            individual: individual.clone(),
            initial_opportunities,
            sampled_count: 0,
            accepted_count: 0,
        });
    }
    let lambda = request.rate_per_opportunity * initial_opportunities as f64;
    if !lambda.is_finite() {
        return Err(IndelApplicationError::InvalidRate);
    }
    let mut count_rng = indel_rng(root_seed, individual, reference, request, true)?;
    let sampled_count =
        poisson_count(&mut count_rng, lambda).map_err(IndelApplicationError::Random)?;
    apply_sampled_indels(
        individual,
        reference,
        root_seed,
        lengths,
        request,
        sampled_count,
    )
}

/// Aplica inserciones y después deleciones sobre el estado producido por ellas.
///
/// Esta prioridad local solo ordena las dos clases de F2.2.9. No compone SNP
/// con indels ni decide conflictos globales, responsabilidades de F2.2.10.
pub fn apply_random_small_indels(
    individual: &Individual,
    reference: &Reference,
    root_seed: RootSeed,
    lengths: SmallIndelLengths,
    request: &RandomSmallIndelRequest<'_>,
) -> Result<SmallIndelApplicationResult, IndelApplicationError> {
    let insertion = apply_random_indels(
        individual,
        reference,
        root_seed,
        lengths,
        &RandomIndelRequest {
            scenario_id: request.scenario_id,
            sequence_id: request.sequence_id,
            phase: request.phase,
            spatial_scope_id: request.spatial_scope_id,
            class: IndelClass::Insertion,
            rate_per_opportunity: request.insertion_rate_per_opportunity,
        },
    )?;
    let deletion = apply_random_indels(
        insertion.individual(),
        reference,
        root_seed,
        lengths,
        &RandomIndelRequest {
            scenario_id: request.scenario_id,
            sequence_id: request.sequence_id,
            phase: request.phase,
            spatial_scope_id: request.spatial_scope_id,
            class: IndelClass::Deletion,
            rate_per_opportunity: request.deletion_rate_per_opportunity,
        },
    )?;
    Ok(SmallIndelApplicationResult {
        individual: deletion.individual().clone(),
        insertion,
        deletion,
    })
}

/// Aplica recuentos ya decididos en el orden local inserción→deleción.
///
/// Esta frontera permite comprobar el selector dinámico sin sustituir el único
/// productor Poisson de [`apply_random_small_indels`].
pub fn apply_sampled_small_indels(
    individual: &Individual,
    reference: &Reference,
    root_seed: RootSeed,
    lengths: SmallIndelLengths,
    request: &RandomSmallIndelRequest<'_>,
    insertion_count: u64,
    deletion_count: u64,
) -> Result<SmallIndelApplicationResult, IndelApplicationError> {
    let insertion = apply_sampled_indels(
        individual,
        reference,
        root_seed,
        lengths,
        &RandomIndelRequest {
            scenario_id: request.scenario_id,
            sequence_id: request.sequence_id,
            phase: request.phase,
            spatial_scope_id: request.spatial_scope_id,
            class: IndelClass::Insertion,
            rate_per_opportunity: request.insertion_rate_per_opportunity,
        },
        insertion_count,
    )?;
    let deletion = apply_sampled_indels(
        insertion.individual(),
        reference,
        root_seed,
        lengths,
        &RandomIndelRequest {
            scenario_id: request.scenario_id,
            sequence_id: request.sequence_id,
            phase: request.phase,
            spatial_scope_id: request.spatial_scope_id,
            class: IndelClass::Deletion,
            rate_per_opportunity: request.deletion_rate_per_opportunity,
        },
        deletion_count,
    )?;
    Ok(SmallIndelApplicationResult {
        individual: deletion.individual().clone(),
        insertion,
        deletion,
    })
}

/// Aplica un recuento ya decidido. Es la frontera de prueba del selector dinámico.
pub fn apply_sampled_indels(
    individual: &Individual,
    reference: &Reference,
    root_seed: RootSeed,
    lengths: SmallIndelLengths,
    request: &RandomIndelRequest<'_>,
    sampled_count: u64,
) -> Result<IndelApplicationResult, IndelApplicationError> {
    if individual.reference_id() != reference.id() {
        return Err(IndelApplicationError::IndividualReferenceMismatch);
    }
    lengths.validate_length(lengths.length_for(request.class), request.class)?;
    let initial_opportunities = enumerate_opportunities(
        current_bases(individual, reference, request.sequence_id, request.phase)?,
        request.class,
        lengths.length_for(request.class),
    )
    .len()
    .try_into()
    .map_err(|_| IndelApplicationError::ArithmeticOverflow)?;
    if sampled_count == 0 {
        return Ok(IndelApplicationResult {
            individual: individual.clone(),
            initial_opportunities,
            sampled_count,
            accepted_count: 0,
        });
    }
    let mut position_rng = indel_rng(root_seed, individual, reference, request, false)?;
    let mut updated = individual.clone();
    let length = lengths.length_for(request.class);
    for event_index in 0..sampled_count {
        let candidates = enumerate_opportunities(
            current_bases(&updated, reference, request.sequence_id, request.phase)?,
            request.class,
            length,
        );
        if candidates.is_empty() {
            return Err(IndelApplicationError::ExhaustedOpportunities {
                sampled_count,
                accepted_before_failure: event_index,
                active_before_rollback: event_index,
                opportunities_at_failure: 0,
            });
        }
        let weights = vec![1.0; candidates.len()];
        let selected =
            categorical_f64(&mut position_rng, &weights).map_err(IndelApplicationError::Random)?;
        let copy = updated
            .copies_for(request.sequence_id)
            .and_then(|copies| copies.iter().find(|copy| copy.phase() == request.phase))
            .ok_or(IndelApplicationError::MissingTargetCopy {
                phase: request.phase,
            })?;
        let event_id = derive_indel_event_id(
            root_seed,
            reference,
            &updated,
            request,
            event_index,
            u64::try_from(copy.segments().current_indel_events().len())
                .map_err(|_| IndelApplicationError::ArithmeticOverflow)?,
        )?;
        let provenance = ModificationProvenance::new(event_id);
        let segments = match candidates[selected] {
            Opportunity::Insertion { point } => {
                let mut content_rng = derive_subseed(
                    root_seed,
                    RuntimeSubseedScope::InsertionContent {
                        scenario_id: request.scenario_id.as_str(),
                        reference_id: reference.id().as_str(),
                        individual_id: updated.id().as_str(),
                        sequence_id: request.sequence_id.as_str(),
                        spatial_scope_id: request.spatial_scope_id,
                        event_index,
                    },
                )
                .map_err(|error| IndelApplicationError::Randomness {
                    detail: error.render_human(),
                })?
                .local_rng();
                let bases = random_insertion_bases(&mut content_rng, length)?;
                copy.segments()
                    .with_current_insertion(point, Arc::from(bases), provenance)
                    .map_err(IndelApplicationError::Segment)?
            }
            Opportunity::Deletion { start, end } => copy
                .segments()
                .with_current_deletion(start, end, provenance)
                .map_err(IndelApplicationError::Segment)?,
        };
        updated = updated
            .with_copy_segments(request.sequence_id, request.phase, segments)
            .map_err(|_| IndelApplicationError::MissingTargetCopy {
                phase: request.phase,
            })?;
    }
    Ok(IndelApplicationResult {
        individual: updated,
        initial_opportunities,
        sampled_count,
        accepted_count: sampled_count,
    })
}

fn indel_rng(
    root_seed: RootSeed,
    individual: &Individual,
    reference: &Reference,
    request: &RandomIndelRequest<'_>,
    count: bool,
) -> Result<LocalRng, IndelApplicationError> {
    let scope = if count {
        RuntimeSubseedScope::IndelCount {
            scenario_id: request.scenario_id.as_str(),
            reference_id: reference.id().as_str(),
            individual_id: individual.id().as_str(),
            sequence_id: request.sequence_id.as_str(),
            variant_class: request.class.as_str(),
            spatial_scope_id: request.spatial_scope_id,
        }
    } else {
        RuntimeSubseedScope::IndelPosition {
            scenario_id: request.scenario_id.as_str(),
            reference_id: reference.id().as_str(),
            individual_id: individual.id().as_str(),
            sequence_id: request.sequence_id.as_str(),
            variant_class: request.class.as_str(),
            spatial_scope_id: request.spatial_scope_id,
        }
    };
    derive_subseed(root_seed, scope)
        .map_err(|error| IndelApplicationError::Randomness {
            detail: error.render_human(),
        })
        .map(crate::randomness::Subseed::local_rng)
}

fn derive_indel_event_id(
    root_seed: RootSeed,
    reference: &Reference,
    individual: &Individual,
    request: &RandomIndelRequest<'_>,
    event_index: u64,
    history_index: u64,
) -> Result<EventId, IndelApplicationError> {
    let event_key = format!(
        "scenario={};reference={};individual={};sequence={};phase={};class={};scope={};ordinal={event_index};history={history_index}",
        request.scenario_id.as_str(),
        reference.id().as_str(),
        individual.id().as_str(),
        request.sequence_id.as_str(),
        request.phase.copy_index(),
        request.class.as_str(),
        request.spatial_scope_id,
    );
    let request = AutomaticIdRequest::new(
        IdClass::Event,
        "f2.2.9",
        &[("event_key", event_key.as_str())],
    )
    .map_err(|error| IndelApplicationError::Randomness {
        detail: error.render_human(),
    })?;
    IdentityRegistry::default()
        .derive_event(root_seed, &request)
        .map_err(|error| IndelApplicationError::Randomness {
            detail: error.render_human(),
        })
}

fn current_bases(
    individual: &Individual,
    reference: &Reference,
    sequence_id: &SequenceId,
    phase: HaplotypePhase,
) -> Result<Vec<u8>, IndelApplicationError> {
    let copy = individual
        .copies_for(sequence_id)
        .and_then(|copies| copies.iter().find(|copy| copy.phase() == phase))
        .ok_or(IndelApplicationError::MissingTargetCopy { phase })?;
    materialize_copy(copy, reference).map_err(IndelApplicationError::Materialization)
}

#[derive(Clone, Copy)]
enum Opportunity {
    Insertion { point: u64 },
    Deletion { start: u64, end: u64 },
}

fn enumerate_opportunities(bases: Vec<u8>, class: IndelClass, length: u64) -> Vec<Opportunity> {
    let Ok(sequence_length) = u64::try_from(bases.len()) else {
        return Vec::new();
    };
    match class {
        IndelClass::Insertion => (1..sequence_length)
            .filter(|point| canonical_base_at(&bases, point - 1))
            .map(|point| Opportunity::Insertion { point })
            .collect(),
        IndelClass::Deletion => {
            let Some(last_start) = sequence_length.checked_sub(length + 1) else {
                return Vec::new();
            };
            (1..=last_start)
                .filter_map(|start| {
                    let end = start.checked_add(length)?;
                    (canonical_base_at(&bases, start - 1) && canonical_base_at(&bases, end))
                        .then_some(Opportunity::Deletion { start, end })
                })
                .collect()
        }
    }
}

fn canonical_base_at(bases: &[u8], offset: u64) -> bool {
    usize::try_from(offset)
        .ok()
        .and_then(|offset| bases.get(offset))
        .is_some_and(|base| matches!(*base, b'A' | b'C' | b'G' | b'T'))
}

fn random_insertion_bases(
    rng: &mut LocalRng,
    length: u64,
) -> Result<Vec<u8>, IndelApplicationError> {
    let capacity =
        usize::try_from(length).map_err(|_| IndelApplicationError::ArithmeticOverflow)?;
    let mut bases = Vec::with_capacity(capacity);
    for _ in 0..length {
        let base = match uniform_u64_below(rng, 4).map_err(IndelApplicationError::Random)? {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        bases.push(base);
    }
    Ok(bases)
}
