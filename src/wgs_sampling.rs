//! Muestreo WGS mínimo desde la vista haplotípica faseada en memoria.
//!
//! Este módulo no escribe FASTQ ni conserva procedencia durable. Cada llamada
//! al iterador produce una única lectura verdadera y la siguiente descarta su
//! búfer, por lo que el estado residente no crece con el número de lecturas.

use crate::application::{
    EffectiveTechnicalEvent, EffectiveTechnicalEventTarget, ILLUMINA_SE_MINIMAL_READ_LENGTH,
};
use crate::chimera::{ChimeraAccumulator, ChimeraReport, ChimeraSettings, transform_paired_end};
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

/// Plan cuantitativo de una biblioteca Illumina PE antes de formar fragmentos.
///
/// El denominador conserva la referencia estructural haploide completa: las
/// secuencias demasiado cortas para el inserto mínimo siguen contando aquí,
/// pero se declaran aparte con cobertura básica cero. Así una deleción o una
/// referencia de mapeo no puede cambiar el significado de la cobertura
/// solicitada.
#[derive(Clone, Debug, PartialEq)]
pub struct WgsPairedEndPlan {
    target_coverage: f64,
    pair_count: u64,
    haploid_reference_length: u64,
    emitted_bases_per_pair: u64,
    equivalent_coverage: f64,
    zero_coverage_sequences: Vec<SequenceId>,
}

impl WgsPairedEndPlan {
    /// Cobertura común solicitada para el evento, en veces la referencia
    /// estructural haploide.
    #[must_use]
    pub const fn target_coverage(&self) -> f64 {
        self.target_coverage
    }

    /// Número entero de parejas derivado del objetivo, nunca una entrada de
    /// configuración independiente.
    #[must_use]
    pub const fn pair_count(&self) -> u64 {
        self.pair_count
    }

    /// Denominador estructural: todas las secuencias de referencia una vez,
    /// sin ploidía y sin filtrar las que no admiten un inserto PE.
    #[must_use]
    pub const fn haploid_reference_length(&self) -> u64 {
        self.haploid_reference_length
    }

    /// Bases emitidas por pareja. Incluye R1 y R2 aunque sus intervalos se
    /// solapen sobre el mismo inserto.
    #[must_use]
    pub const fn emitted_bases_per_pair(&self) -> u64 {
        self.emitted_bases_per_pair
    }

    /// Cobertura que equivale al número entero de parejas derivado.
    #[must_use]
    pub const fn equivalent_coverage(&self) -> f64 {
        self.equivalent_coverage
    }

    /// Secuencias estructurales excluidas del muestreo PE por no alojar el
    /// inserto mínimo. Su cobertura WGS básica es exactamente cero.
    #[must_use]
    pub fn zero_coverage_sequences(&self) -> &[SequenceId] {
        &self.zero_coverage_sequences
    }
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
    InvalidPairedEndCoverage,
    PairCountOverflow,
    Random(RandomSamplingError),
    TechnicalIdentity(Diagnostic),
    Materialization(MaterializationError),
    Chimera(crate::chimera::ChimeraError),
}

/// Resuelve la cobertura F4 a parejas Illumina PE completas.
///
/// La conversión es determinista: `target_x * haploid_length / (2 * read_length)`
/// se redondea a la pareja más cercana y un empate de media pareja asciende.
/// No usa el contador Poisson single-end de F3: la estocasticidad básica PE
/// pertenece a los fragmentos, copias y posiciones seleccionados después.
pub fn resolve_wgs_paired_end_plan(
    reference: &Reference,
    profile: F4TechnicalProfile,
    read_length: u16,
    target_coverage: f64,
) -> Result<WgsPairedEndPlan, WgsSamplingError> {
    if !target_coverage.is_finite() || target_coverage < 0.0 {
        return Err(WgsSamplingError::InvalidPairedEndCoverage);
    }
    if !profile.common_read_length().contains(read_length)
        || read_length > profile.insert_length().minimum_bases()
    {
        return Err(WgsSamplingError::NoEligibleSequence {
            read_length: u64::from(read_length),
        });
    }

    let haploid_reference_length =
        reference
            .sequences()
            .iter()
            .try_fold(0_u64, |total, sequence| {
                total
                    .checked_add(sequence.length())
                    .ok_or(WgsSamplingError::LengthOverflow)
            })?;
    let emitted_bases_per_pair = u64::from(read_length)
        .checked_mul(2)
        .ok_or(WgsSamplingError::LengthOverflow)?;
    let raw_pair_count =
        target_coverage * (haploid_reference_length as f64) / (emitted_bases_per_pair as f64);
    if !raw_pair_count.is_finite() || raw_pair_count > u64::MAX as f64 {
        return Err(WgsSamplingError::PairCountOverflow);
    }
    // `round` se aleja de cero en el empate; el dominio ya garantiza que el
    // valor es no negativo, por lo que implementa exactamente empate hacia
    // arriba sin una tolerancia numérica ad hoc.
    let pair_count = raw_pair_count.round() as u64;
    let equivalent_coverage = if haploid_reference_length == 0 {
        0.0
    } else {
        (pair_count as f64) * (emitted_bases_per_pair as f64) / (haploid_reference_length as f64)
    };
    let minimum_insert = u64::from(profile.insert_length().minimum_bases());
    let zero_coverage_sequences = reference
        .sequences()
        .iter()
        .filter(|sequence| sequence.length() < minimum_insert)
        .map(|sequence| sequence.id().clone())
        .collect();

    Ok(WgsPairedEndPlan {
        target_coverage,
        pair_count,
        haploid_reference_length,
        emitted_bases_per_pair,
        equivalent_coverage,
        zero_coverage_sequences,
    })
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
    effective_insert_maximum_bases: u16,
    coverage_plan: Option<WgsPairedEndPlan>,
    remaining_pairs: u64,
    next_ordinal: u64,
    identity_deriver: F4TechnicalIdentityDeriver,
    chimera_settings: ChimeraSettings,
    chimera_accumulator: ChimeraAccumulator,
    failed: bool,
}

struct PairCountConstruction {
    pair_count: u64,
    coverage_plan: Option<WgsPairedEndPlan>,
    allow_empty_eligible_space: bool,
}

impl<'a> WgsPairedEndSampler<'a> {
    /// Construye el único muestreador PE admisible en F4 a partir de la
    /// cobertura común del evento. No expone un objetivo por número de parejas.
    pub fn new_for_coverage(
        reference: &'a Reference,
        individual: &'a Individual,
        identities: &'a TechnicalIdentityPlan,
        profile: F4TechnicalProfile,
        read_length: u16,
        target_coverage: f64,
    ) -> Result<Self, WgsSamplingError> {
        let plan = resolve_wgs_paired_end_plan(reference, profile, read_length, target_coverage)?;
        Self::new_with_pair_count(
            reference,
            individual,
            identities,
            profile,
            read_length,
            PairCountConstruction {
                pair_count: plan.pair_count(),
                coverage_plan: Some(plan),
                allow_empty_eligible_space: true,
            },
        )
    }

    /// Alias explícito para consumidores que nombran el campo de configuración
    /// `target_coverage_wgs`.
    pub fn new_for_target_coverage(
        reference: &'a Reference,
        individual: &'a Individual,
        identities: &'a TechnicalIdentityPlan,
        profile: F4TechnicalProfile,
        read_length: u16,
        target_coverage: f64,
    ) -> Result<Self, WgsSamplingError> {
        Self::new_for_coverage(
            reference,
            individual,
            identities,
            profile,
            read_length,
            target_coverage,
        )
    }

    /// Plan de cobertura que originó esta biblioteca F4.
    ///
    /// El constructor interno `new_for_pair_count` de F4.1 no tiene uno: su
    /// entero solo permite probar la formación de pares antes de F4.2.
    #[must_use]
    pub fn coverage_plan(&self) -> Option<&WgsPairedEndPlan> {
        self.coverage_plan.as_ref()
    }

    /// Límite superior efectivo del inserto sorteado para esta biblioteca.
    ///
    /// Si alguna copia elegible es más corta que el máximo del perfil, el
    /// sorteo se condiciona al intervalo que puede alojar al menos una copia.
    /// Este valor permite que el informe distinga ese condicionamiento de los
    /// límites configurados por el perfil. La elección posterior de copia y
    /// comienzo sigue filtrando las candidatas que alojan cada inserto concreto.
    #[must_use]
    pub const fn effective_insert_maximum_bases(&self) -> u16 {
        self.effective_insert_maximum_bases
    }

    pub fn new_for_pair_count(
        reference: &'a Reference,
        individual: &'a Individual,
        identities: &'a TechnicalIdentityPlan,
        profile: F4TechnicalProfile,
        read_length: u16,
        pair_count: u64,
    ) -> Result<Self, WgsSamplingError> {
        Self::new_with_pair_count(
            reference,
            individual,
            identities,
            profile,
            read_length,
            PairCountConstruction {
                pair_count,
                coverage_plan: None,
                allow_empty_eligible_space: false,
            },
        )
    }

    fn new_with_pair_count(
        reference: &'a Reference,
        individual: &'a Individual,
        identities: &'a TechnicalIdentityPlan,
        profile: F4TechnicalProfile,
        read_length: u16,
        construction: PairCountConstruction,
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
        let effective_insert_maximum_bases = effective_insert_maximum_bases(individual, profile)?;
        let has_eligible_copy = effective_insert_maximum_bases.is_some();
        if construction.pair_count > 0
            && !has_eligible_copy
            && !construction.allow_empty_eligible_space
        {
            return Err(WgsSamplingError::NoEligibleSequence {
                read_length: u64::from(profile.insert_length().minimum_bases()),
            });
        }
        let chimera_settings =
            ChimeraSettings::new(profile.chimera_rate()).map_err(WgsSamplingError::Chimera)?;
        Ok(Self {
            reference,
            individual,
            identities,
            profile,
            read_length,
            // Si hay al menos una copia de tamaño mínimo, el máximo está
            // acotado tanto por el perfil como por la copia más larga. Así no
            // puede salir un inserto sin candidatas (p. ej. 301--450 sobre
            // una biblioteca cuya copia más larga mide 300 pb).
            effective_insert_maximum_bases: effective_insert_maximum_bases
                .unwrap_or_else(|| profile.insert_length().minimum_bases()),
            coverage_plan: construction.coverage_plan,
            // Una referencia puede estar formada solo por secuencias más
            // cortas que el inserto mínimo. F4 las registra a cero y no aborta
            // el evento; no se inventan extremos, inserciones ni reintentos.
            remaining_pairs: if has_eligible_copy {
                construction.pair_count
            } else {
                0
            },
            next_ordinal: 0,
            identity_deriver: F4TechnicalIdentityDeriver::new(
                identities.root_seed(),
                identities.event_id().clone(),
            ),
            chimera_settings,
            chimera_accumulator: ChimeraAccumulator::default(),
            failed: false,
        })
    }

    /// Configura la transformación de quimeras para este flujo. El constructor
    /// conserva el valor de perfil inicial (`0`) hasta que una configuración F4
    /// posterior autorice una tasa explícita.
    #[must_use]
    pub fn with_chimera_settings(mut self, settings: ChimeraSettings) -> Self {
        self.chimera_settings = settings;
        self
    }

    /// Conecta la tasa efectiva del evento al productor PE mediante el tipo
    /// cerrado `ChimeraSettings`.
    #[must_use]
    pub fn with_effective_technical_event(mut self, event: &EffectiveTechnicalEvent) -> Self {
        self.chimera_settings = event.chimera_settings();
        self
    }

    /// Recuento agregado de quimeras ya emitidas por este flujo.
    #[must_use]
    pub const fn chimera_report(&self) -> ChimeraReport {
        self.chimera_accumulator.report()
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
            i64::from(self.effective_insert_maximum_bases),
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
        let pair = emit_paired_end(fragment, pre_hop_identity, self.read_length)
            .map_err(WgsSamplingError::TechnicalIdentity)?;
        let transformation = transform_paired_end(
            self.identities.root_seed(),
            event_id,
            self.reference,
            self.individual,
            pair,
            self.chimera_settings,
        )
        .map_err(WgsSamplingError::Chimera)?;
        self.chimera_accumulator
            .record(&transformation)
            .map_err(WgsSamplingError::Chimera)?;
        Ok(transformation.into_pair())
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

/// Devuelve el máximo de la distribución que puede alojar alguna copia actual.
///
/// Se conserva el mínimo del perfil como frontera de elegibilidad; por debajo
/// de él no hay espacio PE. Sobre ese mínimo, condicionar la normal truncada
/// evita que una muestra válida en el perfil sea imposible para todas las
/// copias de esta biblioteca.
fn effective_insert_maximum_bases(
    individual: &Individual,
    profile: F4TechnicalProfile,
) -> Result<Option<u16>, WgsSamplingError> {
    let minimum = u64::from(profile.insert_length().minimum_bases());
    let maximum = u64::from(profile.insert_length().maximum_bases());
    let mut largest = None;
    for genotype in individual.sequences() {
        for copy in genotype.copies() {
            let current_length = copy
                .segments()
                .current_length()
                .map_err(|_| WgsSamplingError::LengthOverflow)?;
            if current_length >= minimum {
                largest = Some(largest.unwrap_or(0_u64).max(current_length));
            }
        }
    }
    let Some(largest) = largest else {
        return Ok(None);
    };
    u16::try_from(largest.min(maximum))
        .map(Some)
        .map_err(|_| WgsSamplingError::LengthOverflow)
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
