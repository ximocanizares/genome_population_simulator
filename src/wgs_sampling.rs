//! Muestreo WGS mínimo desde la vista haplotípica faseada en memoria.
//!
//! Este módulo no escribe FASTQ ni conserva procedencia durable. Cada llamada
//! al iterador produce una única lectura verdadera y la siguiente descarta su
//! búfer, por lo que el estado residente no crece con el número de lecturas.

use crate::application::{EffectiveTechnicalEventTarget, ILLUMINA_SE_MINIMAL_READ_LENGTH};
use crate::diagnostics::Diagnostic;
use crate::haplotype::{HaplotypeCopy, HaplotypeCopyId, Individual};
use crate::identity::{EventId, IndividualId, LibraryId, SequenceId};
use crate::materialization::{MaterializationError, materialize_copy_interval};
use crate::paired_end::{PairedEndFragment, PairedEndPair, emit_paired_end};
use crate::randomness::{
    F4TechnicalDecision, RandomSamplingError, derive_f4_technical_subseed, f4_truncated_normal,
    f4_uniform_index, poisson_count, uniform_u64_below,
};
use crate::reference::{Reference, ReferenceSequence};
use crate::technical_identity::{
    F4TechnicalIdentityDeriver, TechnicalIdentityPlan, TechnicalTemplateIdentity, TemplateOrdinal,
};
use crate::technical_profile::F4TechnicalProfile;

pub const WGS_READ_LENGTH: u64 = ILLUMINA_SE_MINIMAL_READ_LENGTH as u64;
pub const WGS_COUNT_ALGORITHM: &str = "wgs-coverage-poisson-v1";

/// Orientación de la secuencia verdadera de una lectura single-end.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReadOrientation {
    Forward,
    ReverseComplement,
}

/// Intervalo semiabierto en las coordenadas vigentes de una copia haplotípica.
///
/// Después de un indel estas coordenadas no se fuerzan a ser un intervalo de
/// referencia: conservan exactamente el espacio desde el que se extrajo la
/// secuencia verdadera.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HaplotypeReadInterval {
    sequence_id: SequenceId,
    start: u64,
    end: u64,
    sequence_length: u64,
}

impl HaplotypeReadInterval {
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
    pub const fn sequence_length(&self) -> u64 {
        self.sequence_length
    }
}

/// Origen transitorio de una plantilla WGS antes de F3.1.4.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SampledWgsRead {
    template: TechnicalTemplateIdentity,
    individual_id: IndividualId,
    haplotype_copy_id: HaplotypeCopyId,
    interval: HaplotypeReadInterval,
    orientation: ReadOrientation,
    true_sequence: Vec<u8>,
}

impl SampledWgsRead {
    #[must_use]
    pub const fn template(&self) -> &TechnicalTemplateIdentity {
        &self.template
    }

    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }

    #[must_use]
    pub const fn haplotype_copy_id(&self) -> &HaplotypeCopyId {
        &self.haplotype_copy_id
    }

    #[must_use]
    pub const fn interval(&self) -> &HaplotypeReadInterval {
        &self.interval
    }

    #[must_use]
    pub const fn orientation(&self) -> ReadOrientation {
        self.orientation
    }

    /// Secuencia verdadera ya orientada como se entregará al siguiente paso.
    #[must_use]
    pub fn true_sequence(&self) -> &[u8] {
        &self.true_sequence
    }
}

/// Resultado estable de resolver el objetivo de una biblioteca.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WgsReadPlan {
    read_count: u64,
    haploid_eligible_length: u64,
}

impl WgsReadPlan {
    #[must_use]
    pub const fn read_count(self) -> u64 {
        self.read_count
    }

    /// Denominador WGS: cada secuencia de la referencia se cuenta una vez,
    /// sin multiplicar por la ploidía.
    #[must_use]
    pub const fn haploid_eligible_length(self) -> u64 {
        self.haploid_eligible_length
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum WgsSamplingError {
    TechnicalPlanReferenceMismatch,
    IndividualReferenceMismatch,
    MissingTechnicalLibrary,
    NoEligibleSequence { read_length: u64 },
    LengthOverflow,
    InvalidCoverageMean,
    Random(RandomSamplingError),
    TechnicalIdentity(Diagnostic),
    Materialization(MaterializationError),
}

struct EligibleCopy<'a> {
    copy: &'a HaplotypeCopy,
    current_length: u64,
    start_count: u64,
}

/// Flujo por biblioteca de orígenes WGS. Su memoria es O(copias + segmentos),
/// independiente de `read_count`.
pub struct WgsLibrarySampler<'a> {
    reference: &'a Reference,
    individual: &'a Individual,
    identities: &'a TechnicalIdentityPlan,
    eligible_copies: Vec<EligibleCopy<'a>>,
    total_start_count: u64,
    plan: WgsReadPlan,
    next_ordinal: u64,
    remaining: u64,
    failed: bool,
}

impl<'a> WgsLibrarySampler<'a> {
    pub fn new(
        reference: &'a Reference,
        individual: &'a Individual,
        identities: &'a TechnicalIdentityPlan,
        target: EffectiveTechnicalEventTarget,
    ) -> Result<Self, WgsSamplingError> {
        if identities.reference_id() != reference.id() {
            return Err(WgsSamplingError::TechnicalPlanReferenceMismatch);
        }
        if individual.reference_id() != reference.id() {
            return Err(WgsSamplingError::IndividualReferenceMismatch);
        }
        if identities.library(individual.id()).is_none() {
            return Err(WgsSamplingError::MissingTechnicalLibrary);
        }
        let haploid_eligible_length = eligible_reference_length(reference)?;
        let (eligible_copies, total_start_count) = eligible_copies(reference, individual)?;
        if total_start_count == 0 {
            return Err(WgsSamplingError::NoEligibleSequence {
                read_length: WGS_READ_LENGTH,
            });
        }
        let read_count =
            resolve_read_count(identities, individual.id(), target, haploid_eligible_length)?;
        Ok(Self {
            reference,
            individual,
            identities,
            eligible_copies,
            total_start_count,
            plan: WgsReadPlan {
                read_count,
                haploid_eligible_length,
            },
            next_ordinal: 0,
            remaining: read_count,
            failed: false,
        })
    }

    #[must_use]
    pub const fn plan(&self) -> WgsReadPlan {
        self.plan
    }

    /// Identidad del individuo cuya única biblioteca recorre este flujo.
    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        self.individual.id()
    }

    /// Evento técnico que causó las identidades emitidas por este flujo.
    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        self.identities.event_id()
    }

    /// Biblioteca técnica de este flujo.
    ///
    /// La construcción ya comprobó que existe una biblioteca para el
    /// individuo; `None` solo protege esta frontera frente a una regresión de
    /// esa invariante compartida.
    #[must_use]
    pub fn library_id(&self) -> Option<&LibraryId> {
        self.identities
            .library(self.individual.id())
            .map(|library| library.library_id())
    }

    #[must_use]
    pub const fn resident_candidate_count(&self) -> usize {
        self.eligible_copies.len()
    }

    fn sample_next(&mut self) -> Result<SampledWgsRead, WgsSamplingError> {
        let ordinal = TemplateOrdinal::new(self.next_ordinal);
        let template = self
            .identities
            .derive_template(self.individual.id(), ordinal)
            .map_err(WgsSamplingError::TechnicalIdentity)?;
        let mut rng = template.subseed().local_rng();
        let candidate = self.pick_copy(&mut rng)?;
        let start =
            uniform_u64_below(&mut rng, candidate.start_count).map_err(WgsSamplingError::Random)?;
        let end = start
            .checked_add(WGS_READ_LENGTH)
            .ok_or(WgsSamplingError::LengthOverflow)?;
        let mut true_sequence =
            materialize_copy_interval(candidate.copy, self.reference, start, end)
                .map_err(WgsSamplingError::Materialization)?;
        let orientation = match uniform_u64_below(&mut rng, 2).map_err(WgsSamplingError::Random)? {
            0 => ReadOrientation::Forward,
            1 => {
                reverse_complement_in_place(&mut true_sequence);
                ReadOrientation::ReverseComplement
            }
            _ => unreachable!("a two-element uniform sample is bounded"),
        };
        Ok(SampledWgsRead {
            template,
            individual_id: self.individual.id().clone(),
            haplotype_copy_id: candidate.copy.id().clone(),
            interval: HaplotypeReadInterval {
                sequence_id: candidate.copy.id().sequence_id().clone(),
                start,
                end,
                sequence_length: candidate.current_length,
            },
            orientation,
            true_sequence,
        })
    }

    fn pick_copy(
        &self,
        rng: &mut crate::randomness::LocalRng,
    ) -> Result<&EligibleCopy<'a>, WgsSamplingError> {
        let selected =
            uniform_u64_below(rng, self.total_start_count).map_err(WgsSamplingError::Random)?;
        let mut cursor = 0_u64;
        for candidate in &self.eligible_copies {
            cursor = cursor
                .checked_add(candidate.start_count)
                .ok_or(WgsSamplingError::LengthOverflow)?;
            if selected < cursor {
                return Ok(candidate);
            }
        }
        Err(WgsSamplingError::LengthOverflow)
    }
}

impl Iterator for WgsLibrarySampler<'_> {
    type Item = Result<SampledWgsRead, WgsSamplingError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 || self.failed {
            return None;
        }
        let result = self.sample_next();
        match result {
            Ok(read) => {
                self.remaining -= 1;
                self.next_ordinal = self.next_ordinal.saturating_add(1);
                Some(Ok(read))
            }
            Err(error) => {
                self.failed = true;
                Some(Err(error))
            }
        }
    }
}

fn eligible_reference_length(reference: &Reference) -> Result<u64, WgsSamplingError> {
    reference
        .sequences()
        .iter()
        .filter(|sequence| sequence.length() >= WGS_READ_LENGTH)
        .try_fold(0_u64, |total, sequence| {
            total
                .checked_add(sequence.length())
                .ok_or(WgsSamplingError::LengthOverflow)
        })
}

fn eligible_copies<'a>(
    reference: &'a Reference,
    individual: &'a Individual,
) -> Result<(Vec<EligibleCopy<'a>>, u64), WgsSamplingError> {
    let mut copies = Vec::new();
    let mut total_start_count = 0_u64;
    for genotype in individual.sequences() {
        let reference_sequence = reference_sequence(reference, genotype.sequence_id())
            .ok_or(WgsSamplingError::IndividualReferenceMismatch)?;
        if reference_sequence.length() < WGS_READ_LENGTH {
            continue;
        }
        for copy in genotype.copies() {
            let current_length = copy
                .segments()
                .current_length()
                .map_err(|_| WgsSamplingError::LengthOverflow)?;
            let Some(start_count) = current_length
                .checked_sub(WGS_READ_LENGTH)
                .and_then(|value| value.checked_add(1))
            else {
                continue;
            };
            total_start_count = total_start_count
                .checked_add(start_count)
                .ok_or(WgsSamplingError::LengthOverflow)?;
            copies.push(EligibleCopy {
                copy,
                current_length,
                start_count,
            });
        }
    }
    Ok((copies, total_start_count))
}

fn resolve_read_count(
    identities: &TechnicalIdentityPlan,
    individual_id: &IndividualId,
    target: EffectiveTechnicalEventTarget,
    haploid_eligible_length: u64,
) -> Result<u64, WgsSamplingError> {
    match target {
        EffectiveTechnicalEventTarget::ReadCount(read_count) => Ok(read_count),
        EffectiveTechnicalEventTarget::CoverageWgs(coverage) => {
            let lambda = coverage * (haploid_eligible_length as f64) / (WGS_READ_LENGTH as f64);
            if !lambda.is_finite() {
                return Err(WgsSamplingError::InvalidCoverageMean);
            }
            // El ámbito exclusivo de conteo es específico de la biblioteca y
            // no comparte ni estado ni subsemilla con sus plantillas.
            let count_seed = identities
                .derive_library_count_subseed(individual_id)
                .map_err(WgsSamplingError::TechnicalIdentity)?;
            poisson_count(&mut count_seed.local_rng(), lambda).map_err(WgsSamplingError::Random)
        }
    }
}

fn reference_sequence<'a>(
    reference: &'a Reference,
    sequence_id: &SequenceId,
) -> Option<&'a ReferenceSequence> {
    reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id() == sequence_id)
}

fn reverse_complement_in_place(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence {
        *base = match *base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            _ => unreachable!("the reference and inserted-base domains are validated"),
        };
    }
}

/// Muestreador F4 de fragmentos PE antes de cobertura, errores y artefactos.
///
/// El número de parejas no forma parte de la configuración F4.1. Este tipo lo
/// recibe como dato interno y verificable para que las pruebas y F4.2 puedan
/// producir fragmentos sin abrir un objetivo alternativo de usuario.
pub struct WgsPairedEndSampler<'a> {
    reference: &'a Reference,
    individual: &'a Individual,
    identities: &'a TechnicalIdentityPlan,
    profile: F4TechnicalProfile,
    read_length: u16,
    remaining_pairs: u64,
    next_ordinal: u64,
    identity_deriver: F4TechnicalIdentityDeriver,
    failed: bool,
}

impl<'a> WgsPairedEndSampler<'a> {
    pub fn new_for_pair_count(
        reference: &'a Reference,
        individual: &'a Individual,
        identities: &'a TechnicalIdentityPlan,
        profile: F4TechnicalProfile,
        read_length: u16,
        pair_count: u64,
    ) -> Result<Self, WgsSamplingError> {
        if identities.reference_id() != reference.id() {
            return Err(WgsSamplingError::TechnicalPlanReferenceMismatch);
        }
        if individual.reference_id() != reference.id() {
            return Err(WgsSamplingError::IndividualReferenceMismatch);
        }
        if identities.library(individual.id()).is_none() {
            return Err(WgsSamplingError::MissingTechnicalLibrary);
        }
        if !profile.common_read_length().contains(read_length)
            || read_length > profile.insert_length().minimum_bases()
        {
            return Err(WgsSamplingError::NoEligibleSequence {
                read_length: u64::from(read_length),
            });
        }
        if pair_count > 0
            && !has_copy_for_insert(individual, profile.insert_length().minimum_bases())?
        {
            return Err(WgsSamplingError::NoEligibleSequence {
                read_length: u64::from(profile.insert_length().minimum_bases()),
            });
        }
        Ok(Self {
            reference,
            individual,
            identities,
            profile,
            read_length,
            remaining_pairs: pair_count,
            next_ordinal: 0,
            identity_deriver: F4TechnicalIdentityDeriver::new(
                identities.root_seed(),
                identities.event_id().clone(),
            ),
            failed: false,
        })
    }

    fn sample_next_pair(&mut self, ordinal: u64) -> Result<PairedEndPair, WgsSamplingError> {
        let template = self
            .identities
            .derive_template(self.individual.id(), TemplateOrdinal::new(ordinal))
            .map_err(WgsSamplingError::TechnicalIdentity)?;
        let library_id = template.library_id();
        let event_id = self.identities.event_id();
        let distribution = self.profile.insert_length();
        let insert_seed = derive_f4_technical_subseed(
            self.identities.root_seed(),
            F4TechnicalDecision::fragment_insert_length(event_id, library_id, ordinal),
        )
        .map_err(WgsSamplingError::TechnicalIdentity)?;
        let insert = f4_truncated_normal(
            &mut insert_seed.local_rng(),
            f64::from(distribution.mean_bases()),
            f64::from(distribution.standard_deviation_bases()),
            i64::from(distribution.minimum_bases()),
            i64::from(distribution.maximum_bases()),
        )
        .map_err(WgsSamplingError::Random)?;
        let insert = u64::try_from(insert).map_err(|_| WgsSamplingError::LengthOverflow)?;
        let candidates = fragment_candidates(self.reference, self.individual, insert)?;
        let total = candidates.iter().try_fold(0_u64, |sum, candidate| {
            sum.checked_add(candidate.start_count)
                .ok_or(WgsSamplingError::LengthOverflow)
        })?;
        let copy_seed = derive_f4_technical_subseed(
            self.identities.root_seed(),
            F4TechnicalDecision::fragment_copy(event_id, library_id, ordinal),
        )
        .map_err(WgsSamplingError::TechnicalIdentity)?;
        let selected = f4_uniform_index(&mut copy_seed.local_rng(), total)
            .map_err(WgsSamplingError::Random)?;
        let candidate = select_fragment_candidate(&candidates, selected)?;
        let start_seed = derive_f4_technical_subseed(
            self.identities.root_seed(),
            F4TechnicalDecision::fragment_start(event_id, library_id, ordinal),
        )
        .map_err(WgsSamplingError::TechnicalIdentity)?;
        let start = f4_uniform_index(&mut start_seed.local_rng(), candidate.start_count)
            .map_err(WgsSamplingError::Random)?;
        let end = start
            .checked_add(insert)
            .ok_or(WgsSamplingError::LengthOverflow)?;
        let sequence = materialize_copy_interval(candidate.copy, self.reference, start, end)
            .map_err(WgsSamplingError::Materialization)?;
        let fragment =
            PairedEndFragment::new(template, candidate.copy.id().clone(), start, sequence)
                .map_err(WgsSamplingError::TechnicalIdentity)?;
        let pre_hop_identity = self
            .identity_deriver
            .derive_own_pre_hop_pair(fragment.template().clone())
            .map_err(WgsSamplingError::TechnicalIdentity)?;
        emit_paired_end(fragment, pre_hop_identity, self.read_length)
            .map_err(WgsSamplingError::TechnicalIdentity)
    }
}

impl Iterator for WgsPairedEndSampler<'_> {
    type Item = Result<PairedEndPair, WgsSamplingError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed || self.remaining_pairs == 0 {
            return None;
        }
        let result = self.sample_next_pair(self.next_ordinal);
        match result {
            Ok(pair) => {
                self.next_ordinal = self.next_ordinal.saturating_add(1);
                self.remaining_pairs -= 1;
                Some(Ok(pair))
            }
            Err(error) => {
                self.failed = true;
                Some(Err(error))
            }
        }
    }
}

fn has_copy_for_insert(individual: &Individual, insert: u16) -> Result<bool, WgsSamplingError> {
    let insert = u64::from(insert);
    for genotype in individual.sequences() {
        for copy in genotype.copies() {
            if copy
                .segments()
                .current_length()
                .map_err(|_| WgsSamplingError::LengthOverflow)?
                >= insert
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn fragment_candidates<'a>(
    reference: &'a Reference,
    individual: &'a Individual,
    insert: u64,
) -> Result<Vec<EligibleCopy<'a>>, WgsSamplingError> {
    let mut candidates = Vec::new();
    for genotype in individual.sequences() {
        if reference_sequence(reference, genotype.sequence_id()).is_none() {
            return Err(WgsSamplingError::IndividualReferenceMismatch);
        }
        for copy in genotype.copies() {
            let current_length = copy
                .segments()
                .current_length()
                .map_err(|_| WgsSamplingError::LengthOverflow)?;
            let Some(start_count) = current_length
                .checked_sub(insert)
                .and_then(|value| value.checked_add(1))
            else {
                continue;
            };
            candidates.push(EligibleCopy {
                copy,
                current_length,
                start_count,
            });
        }
    }
    if candidates.is_empty() {
        return Err(WgsSamplingError::NoEligibleSequence {
            read_length: insert,
        });
    }
    Ok(candidates)
}

fn select_fragment_candidate<'a>(
    candidates: &'a [EligibleCopy<'a>],
    selected: u64,
) -> Result<&'a EligibleCopy<'a>, WgsSamplingError> {
    let mut cursor = 0_u64;
    for candidate in candidates {
        cursor = cursor
            .checked_add(candidate.start_count)
            .ok_or(WgsSamplingError::LengthOverflow)?;
        if selected < cursor {
            return Ok(candidate);
        }
    }
    Err(WgsSamplingError::LengthOverflow)
}
