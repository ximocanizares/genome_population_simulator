//! Muestreo WGS mínimo desde la vista haplotípica faseada en memoria.
//!
//! Este módulo no escribe FASTQ ni conserva procedencia durable. Cada llamada
//! al iterador produce una única lectura verdadera y la siguiente descarta su
//! búfer, por lo que el estado residente no crece con el número de lecturas.

use crate::application::{
    EffectiveTechnicalEvent, EffectiveTechnicalEventTarget, ILLUMINA_SE_MINIMAL_READ_LENGTH,
};
use crate::chimera::{ChimeraAccumulator, ChimeraReport, ChimeraSettings, transform_paired_end};
use crate::clonal_stutter::LineageNetSteps;
use crate::diagnostics::Diagnostic;
use crate::haplotype::{HaplotypeCopy, HaplotypeCopyId, Individual};
use crate::identity::{EventId, IndividualId, LibraryId, SequenceId};
use crate::materialization::{MaterializationError, materialize_copy_interval};
use crate::paired_end::{PairedEndFragment, PairedEndPair, emit_paired_end};
use crate::pcr_amplification::{PcrAccumulator, PcrPlan, PcrReport};
use crate::randomness::{
    F4TechnicalDecision, FragmentRetryDecision, RandomSamplingError, categorical_f64_with_source,
    derive_f4_technical_subseed, f4_bernoulli, f4_truncated_normal, f4_uniform_index,
    poisson_count, uniform_u64_below,
};
use crate::reference::{Reference, ReferenceSequence};
use crate::stutter::{
    StutterAccumulator, StutterError, StutterPlan, StutterReport, StutterTractIndex,
};
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

// `Eq` decae a `PartialEq` porque el rechazo del tartamudeo nombra el valor
// numérico que lo causó, y un `f64` no admite igualdad total.
#[derive(Clone, Debug, PartialEq)]
pub enum WgsSamplingError {
    TechnicalPlanReferenceMismatch,
    IndividualReferenceMismatch,
    MissingTechnicalLibrary,
    NoEligibleSequence {
        read_length: u64,
    },
    LengthOverflow,
    InvalidCoverageMean,
    InvalidPairedEndCoverage,
    PairCountOverflow,
    /// DT-006: un sorteo del inserto de este individuo agotó la cota de
    /// intentos de `f4_truncated_normal` sin acertar en el intervalo
    /// admitido. Se traduce desde
    /// [`RandomSamplingError::TruncatedNormalAttemptsExhausted`] en vez de
    /// dejarlo como un fallo aleatorio genérico, para nombrar al individuo,
    /// su copia más larga y el mínimo del perfil.
    ImpracticalInsertSampling {
        individual: IndividualId,
        longest_copy_bases: u64,
        minimum_insert_bases: u16,
    },
    Random(RandomSamplingError),
    TechnicalIdentity(Diagnostic),
    Materialization(MaterializationError),
    Chimera(crate::chimera::ChimeraError),
    Stutter(crate::stutter::StutterError),
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
    /// Plan de tartamudeo del evento. El plan por omisión está apagado, de
    /// modo que un muestreador construido sin él se comporta exactamente como
    /// antes de que el tartamudeo existiera.
    stutter_plan: StutterPlan,
    /// Tractos del plan situados en cada copia de la individua, construidos
    /// una sola vez al conectar el plan (F6.4.a.7). Guarda el error, si lo
    /// hubo, para devolverlo en la primera pareja que lo necesite, que es
    /// cuando se notaba antes.
    stutter_tracts: Result<StutterTractIndex, StutterError>,
    stutter_accumulator: StutterAccumulator,
    /// Plan de amplificación del evento. Apagado por omisión, y con él apagado
    /// **no se abre ni una sola decisión nueva**: es lo que conserva byte a
    /// byte los paquetes anteriores a F4.9.
    pcr_plan: PcrPlan,
    pcr_accumulator: PcrAccumulator,
    /// Molécula original de la familia en curso y parejas que le quedan por
    /// emitir. Es todo el estado que la amplificación añade entre llamadas, y
    /// no crece con el número de parejas ni con el genoma.
    family: Option<FamilyState>,
    failed: bool,
}

/// Unidades de operación que cuesta la genealogía de una familia.
///
/// Es el término `K × (2·c·T + c + 1)` de la cota publicada: por cada copia,
/// dos decisiones por ciclo y tracto —el resbalón y, como mucho, su paso—, más
/// el recorrido de los `c + 1` pesos de la divergencia y la elección del
/// progenitor. Se cuenta **en el peor caso**, no en el realizado: la cota que
/// se publica tiene que seguir siendo una cota cuando todos los ciclos
/// resbalen.
fn lineage_units(cycles: u32, family_size: u64, tracts: u64) -> u64 {
    if tracts == 0 {
        return 0;
    }
    let cycles = u64::from(cycles);
    family_size.saturating_mul(
        2_u64
            .saturating_mul(cycles)
            .saturating_mul(tracts)
            .saturating_add(cycles)
            .saturating_add(1),
    )
}

/// Molécula original de una familia, su genealogía y lo que queda por emitir.
#[derive(Clone, Debug)]
struct FamilyState {
    molecule: SampledMolecule,
    remaining: u64,
    /// Índice de la **copia siguiente** dentro de la familia. Es el que entra
    /// en las claves del linaje, de modo que la primera pareja de cada familia
    /// es siempre la copia `0` y no depende de cuántas parejas emitieran las
    /// familias anteriores.
    next_member: usize,
    /// Pasos netos de tartamudeo clonal de cada copia. Vacía cuando el
    /// escenario no declara loci de repetición, que es lo que hace que una
    /// biblioteca amplificada sin microsatélites no abra ni una decisión
    /// nueva.
    lineage: crate::clonal_stutter::FamilyLineage,
}

/// Lo que una familia recién abierta deja anotado para el acumulador.
#[derive(Clone, Copy, Debug)]
struct FamilyAccounting {
    attempts: u64,
    attempts_exhausted: bool,
    emitted_pairs: u32,
    truncated: bool,
    operation_units: u64,
}

/// Una molécula original ya materializada: es lo que hoy produce el tramo de
/// `fragment` del muestreador, extraído a un tipo propio para que una familia
/// entera pueda emitirse desde **la misma** molécula sin volver a sortearla.
#[derive(Clone, Debug)]
struct SampledMolecule {
    copy_id: HaplotypeCopyId,
    start: u64,
    sequence: Vec<u8>,
    /// Copias candidatas enumeradas para elegir esta molécula, y bases del
    /// inserto. Solo alimentan el contador de unidades de la cota publicada.
    candidate_count: u64,
    insert_bases: u64,
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
            stutter_plan: StutterPlan::default(),
            stutter_tracts: Ok(StutterTractIndex::default()),
            stutter_accumulator: StutterAccumulator::default(),
            pcr_plan: PcrPlan::default(),
            pcr_accumulator: PcrAccumulator::default(),
            family: None,
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

    /// Conecta el plan de tartamudeo del evento al productor PE.
    #[must_use]
    pub fn with_stutter_plan(mut self, plan: StutterPlan) -> Self {
        // Con el plan apagado no se sitúa nada: ni el tartamudeo de base ni el
        // linaje consultan el índice.
        self.stutter_tracts = if plan.is_active() {
            StutterTractIndex::build(self.individual, &plan)
        } else {
            Ok(StutterTractIndex::default())
        };
        self.stutter_plan = plan;
        self
    }

    /// El índice de tractos de la individua, o el error que dio construirlo.
    fn stutter_tracts(&self) -> Result<&StutterTractIndex, WgsSamplingError> {
        self.stutter_tracts
            .as_ref()
            .map_err(|error| WgsSamplingError::Stutter(error.clone()))
    }

    /// Conecta el plan de amplificación del evento al productor PE.
    ///
    /// Las unidades que costó construir las tablas de los primeros ciclos se
    /// anotan **aquí y una sola vez**, que es donde se pagan: el plan se
    /// resuelve una vez por evento y las tablas viajan con él. Contarlas por
    /// familia era justo la discrepancia entre el código y la cota que ninguna
    /// prueba podía ver.
    #[must_use]
    pub fn with_pcr_plan(mut self, plan: PcrPlan) -> Self {
        if let Some(settings) = plan.settings() {
            self.pcr_accumulator
                .record_operation_units(settings.early_copy_build_units());
        }
        self.pcr_plan = plan;
        self
    }

    /// Recuentos agregados de la amplificación de este flujo.
    #[must_use]
    pub fn pcr_report(&self) -> PcrReport {
        self.pcr_accumulator.clone().report()
    }

    /// Recuento agregado de quimeras ya emitidas por este flujo.
    #[must_use]
    pub const fn chimera_report(&self) -> ChimeraReport {
        self.chimera_accumulator.report()
    }

    /// Recuentos agregados del tartamudeo de este flujo.
    #[must_use]
    pub fn stutter_report(&self) -> StutterReport {
        self.stutter_accumulator.clone().report()
    }

    /// Traduce el fallo del sorteo del inserto de DT-006 a un diagnóstico que
    /// nombra el individuo, su copia más larga y el mínimo del perfil; deja
    /// pasar sin tocar cualquier otro fallo aleatorio, que sigue su vía
    /// genérica de siempre.
    fn insert_sampling_error(&self, error: RandomSamplingError) -> WgsSamplingError {
        if matches!(error, RandomSamplingError::TruncatedNormalAttemptsExhausted) {
            WgsSamplingError::ImpracticalInsertSampling {
                individual: self.individual.id().clone(),
                longest_copy_bases: u64::from(self.effective_insert_maximum_bases),
                minimum_insert_bases: self.profile.insert_length().minimum_bases(),
            }
        } else {
            WgsSamplingError::Random(error)
        }
    }

    /// Sortea una molécula original con las claves de `fragment`.
    ///
    /// El intento `0` usa **literalmente** las claves de hoy; a partir del `1`,
    /// que solo existe cuando la aceptación por composición rechazó, las claves
    /// llevan además el componente `attempt` y nacen en el catálogo v5. Es esa
    /// separación la que hace que una biblioteca amplificada con curva neutra
    /// consuma exactamente lo mismo que una sin amplificar.
    fn sample_molecule(
        &self,
        ordinal: u64,
        attempt: u64,
    ) -> Result<SampledMolecule, WgsSamplingError> {
        let library_id = self
            .identities
            .library(self.individual.id())
            .ok_or(WgsSamplingError::MissingTechnicalLibrary)?
            .library_id();
        let event_id = self.identities.event_id();
        let distribution = self.profile.insert_length();
        let retry = std::num::NonZeroU64::new(attempt);
        let insert_decision = match retry {
            None => F4TechnicalDecision::fragment_insert_length(event_id, library_id, ordinal),
            Some(attempt) => F4TechnicalDecision::fragment_retry(
                event_id,
                library_id,
                ordinal,
                FragmentRetryDecision::InsertLength,
                attempt,
            ),
        };
        let insert_seed = derive_f4_technical_subseed(self.identities.root_seed(), insert_decision)
            .map_err(WgsSamplingError::TechnicalIdentity)?;
        let insert = f4_truncated_normal(
            &mut insert_seed.local_rng(),
            f64::from(distribution.mean_bases()),
            f64::from(distribution.standard_deviation_bases()),
            i64::from(distribution.minimum_bases()),
            i64::from(self.effective_insert_maximum_bases),
        )
        .map_err(|error| self.insert_sampling_error(error))?;
        let insert = u64::try_from(insert).map_err(|_| WgsSamplingError::LengthOverflow)?;
        let candidates = fragment_candidates(self.reference, self.individual, insert)?;
        let total = candidates.iter().try_fold(0_u64, |sum, candidate| {
            sum.checked_add(candidate.start_count)
                .ok_or(WgsSamplingError::LengthOverflow)
        })?;
        let copy_decision = match retry {
            None => F4TechnicalDecision::fragment_copy(event_id, library_id, ordinal),
            Some(attempt) => F4TechnicalDecision::fragment_retry(
                event_id,
                library_id,
                ordinal,
                FragmentRetryDecision::Copy,
                attempt,
            ),
        };
        let copy_seed = derive_f4_technical_subseed(self.identities.root_seed(), copy_decision)
            .map_err(WgsSamplingError::TechnicalIdentity)?;
        let selected = f4_uniform_index(&mut copy_seed.local_rng(), total)
            .map_err(WgsSamplingError::Random)?;
        let candidate = select_fragment_candidate(&candidates, selected)?;
        let start_decision = match retry {
            None => F4TechnicalDecision::fragment_start(event_id, library_id, ordinal),
            Some(attempt) => F4TechnicalDecision::fragment_retry(
                event_id,
                library_id,
                ordinal,
                FragmentRetryDecision::Start,
                attempt,
            ),
        };
        let start_seed = derive_f4_technical_subseed(self.identities.root_seed(), start_decision)
            .map_err(WgsSamplingError::TechnicalIdentity)?;
        let start = f4_uniform_index(&mut start_seed.local_rng(), candidate.start_count)
            .map_err(WgsSamplingError::Random)?;
        let end = start
            .checked_add(insert)
            .ok_or(WgsSamplingError::LengthOverflow)?;
        let sequence = materialize_copy_interval(candidate.copy, self.reference, start, end)
            .map_err(WgsSamplingError::Materialization)?;
        Ok(SampledMolecule {
            copy_id: candidate.copy.id().clone(),
            start,
            sequence,
            candidate_count: candidates.len() as u64,
            insert_bases: insert,
        })
    }

    /// Molécula de la que sale la pareja de este ordinal.
    ///
    /// Sin amplificación es la de hoy, sorteada una por pareja. Con
    /// amplificación, la de la familia en curso; si la familia se agotó, abre
    /// una nueva en este mismo ordinal, que pasa a ser su **ordinal de
    /// familia**.
    fn molecule_for(
        &mut self,
        ordinal: u64,
    ) -> Result<(SampledMolecule, Option<LineageNetSteps>), WgsSamplingError> {
        if !self.pcr_plan.is_active() {
            return Ok((self.sample_molecule(ordinal, 0)?, None));
        }
        if self
            .family
            .as_ref()
            .is_none_or(|family| family.remaining == 0)
        {
            let (family, accounting) = self.open_family(ordinal)?;
            self.pcr_accumulator
                .record_gc_attempts(accounting.attempts, accounting.attempts_exhausted);
            self.pcr_accumulator
                .record_family(accounting.emitted_pairs, accounting.truncated);
            self.pcr_accumulator
                .record_operation_units(accounting.operation_units);
            self.pcr_accumulator.record_lineage(&family.lineage);
            self.family = Some(family);
        }
        let exhausts = self
            .family
            .as_ref()
            .is_some_and(|family| family.remaining == 1);
        if exhausts {
            // La última pareja de la familia se lleva la molécula en vez de
            // copiarla. Con familias de tamaño 1 —el caso normal con pocos
            // duplicados— eso significa que la amplificación no añade ni una
            // copia de secuencia por pareja.
            let family = self.family.take().ok_or(WgsSamplingError::LengthOverflow)?;
            let steps = family.lineage.member(family.next_member);
            return Ok((family.molecule, steps));
        }
        let family = self
            .family
            .as_mut()
            .ok_or(WgsSamplingError::LengthOverflow)?;
        family.remaining -= 1;
        let steps = family.lineage.member(family.next_member);
        family.next_member += 1;
        Ok((family.molecule.clone(), steps))
    }

    /// Abre una familia: clase de eficiencia, copias tras los primeros ciclos,
    /// molécula aceptada por composición y tamaño.
    ///
    /// **La aceptación por composición va al elegir la molécula**, no al
    /// repartir sus copias. Con pocas parejas por molécula —el caso normal—
    /// casi todas las familias tienen tamaño 1 tanto en regiones ricas como
    /// pobres en GC, de modo que meter el sesgo en el tamaño de familia no
    /// produciría infrarrepresentación ninguna. Lo que la produce de verdad son
    /// las moléculas que **no llegan a observarse**, y eso solo puede
    /// expresarse aquí.
    ///
    /// No toca el acumulador: devuelve lo que hay que anotar para que el
    /// llamador lo haga y esta función pueda trabajar con `&self`, sin clonar
    /// los ajustes una vez por pareja.
    fn open_family(
        &self,
        ordinal: u64,
    ) -> Result<(FamilyState, FamilyAccounting), WgsSamplingError> {
        let settings = self
            .pcr_plan
            .settings()
            .ok_or(WgsSamplingError::LengthOverflow)?;
        let library_id = self
            .identities
            .library(self.individual.id())
            .ok_or(WgsSamplingError::MissingTechnicalLibrary)?
            .library_id();
        let event_id = self.identities.event_id();
        let root = self.identities.root_seed();
        let mut units = 0_u64;

        // 1. Clase de eficiencia. Con una sola clase no se consume nada.
        let classes = settings.efficiency_classes();
        units = units.saturating_add(classes.len() as u64);
        let class_index = if classes.len() == 1 {
            0
        } else {
            let weights: Vec<f64> = classes.iter().map(|class| class.weight()).collect();
            let seed = derive_f4_technical_subseed(
                root,
                F4TechnicalDecision::amplification_efficiency_class(event_id, library_id, ordinal),
            )
            .map_err(WgsSamplingError::TechnicalIdentity)?;
            categorical_f64_with_source(&mut seed.local_rng(), &weights)
                .map_err(WgsSamplingError::Random)?
        };
        let efficiency = classes[class_index].efficiency();
        // La tabla de esta clase está **precalculada**: el triángulo de Pascal
        // y la propagación ciclo a ciclo se pagaron una vez al construir el
        // plan, no una vez por familia. Aquí solo se recorre.
        let early_table = settings
            .early_copy_table(class_index)
            .ok_or(WgsSamplingError::LengthOverflow)?;

        // 2. Copias tras los primeros ciclos, de la distribución **sesgada por
        //    tamaño**: una familia observada es más probable cuantas más copias
        //    tenga su molécula.
        let mean_pairs = settings.mean_pairs_per_molecule();
        let early_copies = if settings.early_cycles() == 0 || mean_pairs == 0.0 {
            1_u64
        } else {
            let weights = early_table.size_biased_weights();
            units = units.saturating_add(weights.len() as u64);
            let seed = derive_f4_technical_subseed(
                root,
                F4TechnicalDecision::amplification_early_copies(event_id, library_id, ordinal),
            )
            .map_err(WgsSamplingError::TechnicalIdentity)?;
            let index = categorical_f64_with_source(&mut seed.local_rng(), weights)
                .map_err(WgsSamplingError::Random)?;
            (index as u64) + 1
        };

        // 3. Molécula original, con aceptación por composición.
        let mut attempt = 0_u64;
        let mut attempts_exhausted = false;
        let (molecule, composition) = loop {
            let molecule = self.sample_molecule(ordinal, attempt)?;
            if attempt > 0 {
                // Solo el trabajo **añadido** entra en el contador: el primer
                // intento es la materialización que F4 ya hacía, y está en el
                // coste por pareja que esta cota deja fuera a propósito.
                units = units.saturating_add(
                    molecule
                        .candidate_count
                        .saturating_add(3_u64.saturating_mul(molecule.insert_bases)),
                );
            }
            let relative = settings
                .gc_curve()
                .evaluate(crate::pcr_amplification::gc_fraction(&molecule.sequence));
            let factor = crate::pcr_amplification::composition_factor(
                efficiency,
                relative,
                settings.cycles(),
            );
            attempt += 1;
            if factor == 1.0 {
                // Curva neutra: el ensayo de Bernoulli no se abre siquiera, de
                // modo que no se deriva ni una subsemilla.
                break (molecule, factor);
            }
            let seed = derive_f4_technical_subseed(
                root,
                F4TechnicalDecision::amplification_gc_accept(
                    event_id,
                    library_id,
                    ordinal,
                    attempt - 1,
                ),
            )
            .map_err(WgsSamplingError::TechnicalIdentity)?;
            if f4_bernoulli(&mut seed.local_rng(), factor).map_err(WgsSamplingError::Random)? {
                break (molecule, factor);
            }
            if attempt >= u64::from(settings.gc_max_attempts()) {
                // Agotado el tope **se acepta el último intento y se cuenta**.
                // El recuento es lo que hace honesto el atajo: el usuario ve
                // cuántas veces ocurrió y juzga si el sesgo le afecta.
                attempts_exhausted = true;
                break (molecule, factor);
            }
        };

        // 4. Tamaño de familia.
        let mean_early = crate::pcr_amplification::EarlyCopyTable::mean_copies(early_table);
        let lambda = mean_pairs * ((early_copies as f64) / mean_early) * composition;
        let drawn = if lambda == 0.0 {
            1_u64
        } else {
            let weights =
                crate::pcr_amplification::family_size_weights(lambda, settings.max_family_size());
            units = units.saturating_add(2_u64.saturating_mul(weights.len() as u64));
            let seed = derive_f4_technical_subseed(
                root,
                F4TechnicalDecision::amplification_family_size(event_id, library_id, ordinal),
            )
            .map_err(WgsSamplingError::TechnicalIdentity)?;
            let index = categorical_f64_with_source(&mut seed.local_rng(), &weights)
                .map_err(WgsSamplingError::Random)?;
            (index as u64) + 1
        };
        // `remaining_pairs` cuenta esta pareja, de modo que aquí nunca vale
        // cero: el iterador no habría llegado hasta aquí.
        let emitted = drawn.min(self.remaining_pairs.max(1));
        let emitted_pairs = u32::try_from(emitted).map_err(|_| WgsSamplingError::LengthOverflow)?;

        // 5. Genealogía de la familia y pasos netos de tartamudeo clonal.
        //
        //    Se resuelve **aquí y para la familia entera** porque la herencia
        //    de una copia necesita el linaje de su progenitor, que es otra
        //    copia de esta misma familia. Resolverlo pareja a pareja obligaría
        //    a rehacer el árbol en cada una.
        let lineage = self.resolve_lineage(ordinal, efficiency, emitted, &molecule)?;
        units = units.saturating_add(lineage_units(
            settings.cycles(),
            emitted,
            lineage.loci_count() as u64,
        ));
        Ok((
            FamilyState {
                molecule,
                remaining: emitted,
                next_member: 0,
                lineage,
            },
            FamilyAccounting {
                attempts: attempt,
                attempts_exhausted,
                emitted_pairs,
                truncated: emitted < drawn,
                operation_units: units,
            },
        ))
    }

    /// Resuelve la genealogía de una familia sobre su molécula original.
    ///
    /// # Qué molécula se mira, y por qué esa
    ///
    /// La **original de la familia**, la misma de la que salen todas sus
    /// parejas. Es lo que hace que las guardas de viabilidad sean coherentes
    /// entre una copia y su progenitor: los dos descienden de la misma
    /// molécula y tienen la misma tirada de partida. Si una pareja resulta
    /// **quimérica**, el paso neto se vuelve a comprobar al aplicarlo, y ahí
    /// puede no caber; eso se decide en `crate::stutter`, no aquí.
    ///
    /// Sin loci declarados, sin tartamudeo de base activo o sin tasas por
    /// ciclo, devuelve la genealogía vacía **sin consumir una sola palabra**.
    fn resolve_lineage(
        &self,
        ordinal: u64,
        efficiency: f64,
        family_size: u64,
        molecule: &SampledMolecule,
    ) -> Result<crate::clonal_stutter::FamilyLineage, WgsSamplingError> {
        let Some(settings) = self.pcr_plan.settings() else {
            return Ok(crate::clonal_stutter::FamilyLineage::default());
        };
        let Some(clonal) = settings.clonal_stutter() else {
            return Ok(crate::clonal_stutter::FamilyLineage::default());
        };
        if !self.stutter_plan.is_active() {
            return Ok(crate::clonal_stutter::FamilyLineage::default());
        }
        let Some(weights) = self.stutter_plan.settings().step_weights() else {
            return Ok(crate::clonal_stutter::FamilyLineage::default());
        };
        let length =
            u64::try_from(molecule.sequence.len()).map_err(|_| WgsSamplingError::LengthOverflow)?;
        let placements = crate::stutter::place_tracts(
            self.stutter_tracts()?,
            &self.stutter_plan,
            &molecule.copy_id,
            molecule.start,
            length,
        )
        .map_err(WgsSamplingError::Stutter)?;
        let mut tracts = Vec::with_capacity(placements.len());
        for placed in &placements {
            let (start, end) = placed.window();
            let start = usize::try_from(start).map_err(|_| WgsSamplingError::LengthOverflow)?;
            let end = usize::try_from(end).map_err(|_| WgsSamplingError::LengthOverflow)?;
            if end > molecule.sequence.len() {
                return Err(WgsSamplingError::LengthOverflow);
            }
            let observed = crate::tandem_repeat::observe_tract(
                &molecule.sequence[start..end],
                placed.locus().motifs(),
            );
            if observed.longest_run_units() == 0 {
                // Sin ninguna copia completa dentro de la molécula no hay
                // tirada donde resbalar, y este tracto no abre ni una
                // decisión.
                continue;
            }
            let motif_length = observed.longest_run_motif_length();
            let Some(probability) =
                clonal.probability_for(motif_length, observed.longest_run_units())
            else {
                // Longitud de motivo fuera del tramo calibrado. No se
                // extrapola; la admisión ya lo rechaza antes de llegar aquí.
                continue;
            };
            tracts.push(crate::clonal_stutter::LineageTract::new(
                placed.locus().id(),
                observed.longest_run_units(),
                motif_length,
                probability,
            ));
        }
        if tracts.is_empty() {
            return Ok(crate::clonal_stutter::FamilyLineage::default());
        }
        let library_id = self
            .identities
            .library(self.individual.id())
            .ok_or(WgsSamplingError::MissingTechnicalLibrary)?
            .library_id();
        let divergence = crate::clonal_stutter::divergence_weights(clonal.cycles(), efficiency);
        crate::clonal_stutter::resolve_family_lineage(
            crate::clonal_stutter::LineageKey {
                root: self.identities.root_seed(),
                event_id: self.identities.event_id(),
                library_id,
                family_ordinal: ordinal,
            },
            family_size,
            &divergence,
            clonal.cycles(),
            &tracts,
            weights,
            crate::clonal_stutter::LineageMolecule {
                length_bases: length,
                read_length_bases: u64::from(self.read_length),
            },
        )
        .map_err(WgsSamplingError::Stutter)
    }

    fn sample_next_pair(&mut self, ordinal: u64) -> Result<PairedEndPair, WgsSamplingError> {
        let template = self
            .identities
            .derive_template(self.individual.id(), TemplateOrdinal::new(ordinal))
            .map_err(WgsSamplingError::TechnicalIdentity)?;
        let event_id = self.identities.event_id();
        let (molecule, lineage) = self.molecule_for(ordinal)?;
        let fragment = PairedEndFragment::new(
            template,
            molecule.copy_id.clone(),
            molecule.start,
            molecule.sequence,
        )
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
        // **Después de la quimera y antes de calidad y errores.** La quimera
        // compara byte a byte el prefijo del fragmento contra la copia
        // haplotípica: un fragmento ya tartamudeado rompería esa comprobación.
        let (pair, chimera_origin) = transformation.into_parts();
        let stuttered = crate::stutter::transform_paired_end(
            self.identities.root_seed(),
            event_id,
            self.stutter_tracts()?,
            &self.stutter_plan,
            lineage.as_ref(),
            chimera_origin.as_ref(),
            pair,
            self.read_length,
        )
        .map_err(WgsSamplingError::Stutter)?;
        self.stutter_accumulator
            .record(&stuttered)
            .map_err(WgsSamplingError::Stutter)?;
        // Los pasos de linaje que **de verdad se escribieron**. Se cuentan
        // aquí, después de aplicarlos, y no al sortearlos: entre una cosa y la
        // otra un paso puede quedarse fuera porque no quepa en la molécula que
        // se va a emitir.
        self.pcr_accumulator
            .record_applied_lineage(stuttered.lineage_applied());
        Ok(stuttered.into_pair())
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
