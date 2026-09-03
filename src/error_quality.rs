//! Transformación transitoria de calidad y error para el perfil F4.
//!
//! La pareja verdadera nunca se modifica. Cada invocación construye una vista
//! observada, la entrega al escritor y la descarta; el único estado retenible
//! son contadores agregados por mate. Calidad es una función pura de ciclo y
//! ajustes, mientras que el Bernoulli y la sustitución usan subsemillas F4
//! distintas por base.

use std::fmt;

use crate::identity::{EventId, LibraryId};
use crate::paired_end::PairedEndPair;
use crate::randomness::{
    F4Mate, F4TechnicalDecision, RandomSamplingError, RootSeed, derive_f4_technical_subseed,
    f4_bernoulli, f4_uniform_index,
};
use crate::technical_identity::TechnicalMate;
use crate::technical_profile::{F4TechnicalProfile, ProfileMate, QualityCurve};

/// Semántica efectiva de las sustituciones F4. `Disabled` es un control
/// experimental: conserva el perfil de calidades, pero esas calidades ya no
/// expresan una probabilidad de error realizada.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SequencingErrorsMode {
    Enabled,
    Disabled,
}

impl SequencingErrorsMode {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Enabled => "enabled",
            Self::Disabled => "disabled",
        }
    }

    #[must_use]
    pub const fn qualities_represent_realized_error_probability(self) -> bool {
        matches!(self, Self::Enabled)
    }
}

/// Ajustes efectivos globales del evento. Sus neutros se materializan siempre
/// aunque la futura configuración los haya omitido.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ErrorQualitySettings {
    sequencing_errors: SequencingErrorsMode,
    quality_shift_r1: i16,
    quality_shift_r2: i16,
    error_scale_r1: f64,
    error_scale_r2: f64,
}

impl Default for ErrorQualitySettings {
    fn default() -> Self {
        Self {
            sequencing_errors: SequencingErrorsMode::Enabled,
            quality_shift_r1: 0,
            quality_shift_r2: 0,
            error_scale_r1: 1.0,
            error_scale_r2: 1.0,
        }
    }
}

impl ErrorQualitySettings {
    pub fn new(
        quality_shift_r1: i16,
        quality_shift_r2: i16,
        error_scale_r1: f64,
        error_scale_r2: f64,
    ) -> Result<Self, ErrorQualityError> {
        if !error_scale_r1.is_finite() || error_scale_r1 < 0.0 {
            return Err(ErrorQualityError::InvalidErrorScale {
                mate: TechnicalMate::R1,
                value: error_scale_r1,
            });
        }
        if !error_scale_r2.is_finite() || error_scale_r2 < 0.0 {
            return Err(ErrorQualityError::InvalidErrorScale {
                mate: TechnicalMate::R2,
                value: error_scale_r2,
            });
        }
        Ok(Self {
            sequencing_errors: SequencingErrorsMode::Enabled,
            quality_shift_r1,
            quality_shift_r2,
            error_scale_r1,
            error_scale_r2,
        })
    }

    #[must_use]
    pub const fn with_sequencing_errors(mut self, mode: SequencingErrorsMode) -> Self {
        self.sequencing_errors = mode;
        self
    }

    #[must_use]
    pub const fn sequencing_errors(self) -> SequencingErrorsMode {
        self.sequencing_errors
    }

    #[must_use]
    pub const fn quality_shift(self, mate: TechnicalMate) -> i16 {
        match mate {
            TechnicalMate::R1 => self.quality_shift_r1,
            TechnicalMate::R2 => self.quality_shift_r2,
        }
    }

    #[must_use]
    pub const fn error_scale(self, mate: TechnicalMate) -> f64 {
        match mate {
            TechnicalMate::R1 => self.error_scale_r1,
            TechnicalMate::R2 => self.error_scale_r2,
        }
    }
}

/// Una lectura observada efímera, apta para codificación FASTQ. No conserva su
/// diferencia con la lectura verdadera.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ObservedMate {
    mate: TechnicalMate,
    sequence: Vec<u8>,
    phred: Vec<u8>,
}

impl ObservedMate {
    #[must_use]
    pub const fn mate(&self) -> TechnicalMate {
        self.mate
    }
    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }
    #[must_use]
    pub fn phred(&self) -> &[u8] {
        &self.phred
    }
}

/// Pareja observada temporal que conserva la cardinalidad PE completa.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ObservedPair {
    r1: ObservedMate,
    r2: ObservedMate,
}

impl ObservedPair {
    #[must_use]
    pub const fn r1(&self) -> &ObservedMate {
        &self.r1
    }
    #[must_use]
    pub const fn r2(&self) -> &ObservedMate {
        &self.r2
    }
}

/// Contadores de una transformación para un único mate.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct MateErrorCounts {
    opportunities: u64,
    errors: u64,
}

impl MateErrorCounts {
    #[must_use]
    pub const fn opportunities(self) -> u64 {
        self.opportunities
    }
    #[must_use]
    pub const fn errors(self) -> u64 {
        self.errors
    }
    #[must_use]
    pub fn realized_rate(self) -> f64 {
        if self.opportunities == 0 {
            0.0
        } else {
            self.errors as f64 / self.opportunities as f64
        }
    }

    pub fn add_assign(&mut self, other: Self) -> Result<(), ErrorQualityError> {
        self.opportunities = self
            .opportunities
            .checked_add(other.opportunities)
            .ok_or(ErrorQualityError::CounterOverflow)?;
        self.errors = self
            .errors
            .checked_add(other.errors)
            .ok_or(ErrorQualityError::CounterOverflow)?;
        Ok(())
    }
}

/// Resultado agregado por pareja, sin origen ni diferencias por base.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ErrorQualityReport {
    r1: MateErrorCounts,
    r2: MateErrorCounts,
}

impl ErrorQualityReport {
    #[must_use]
    pub const fn r1(self) -> MateErrorCounts {
        self.r1
    }
    #[must_use]
    pub const fn r2(self) -> MateErrorCounts {
        self.r2
    }
}

/// Acumulador del evento, deliberadamente agregado por mate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ErrorQualityAccumulator(ErrorQualityReport);

impl ErrorQualityAccumulator {
    pub fn record(&mut self, report: ErrorQualityReport) -> Result<(), ErrorQualityError> {
        self.0.r1.add_assign(report.r1)?;
        self.0.r2.add_assign(report.r2)
    }
    #[must_use]
    pub const fn report(self) -> ErrorQualityReport {
        self.0
    }
}

/// Resultado de aplicar el perfil a una pareja: vista temporal y métricas
/// agregables. Ningún campo representa el diff con su verdad.
#[derive(Clone, Debug, PartialEq)]
pub struct ObservedPairTransformation {
    observed: ObservedPair,
    report: ErrorQualityReport,
}

impl ObservedPairTransformation {
    #[must_use]
    pub const fn observed(&self) -> &ObservedPair {
        &self.observed
    }
    #[must_use]
    pub const fn report(&self) -> ErrorQualityReport {
        self.report
    }

    /// Transfiere la vista observada al siguiente consumidor una vez que sus
    /// métricas ya se han registrado. Evita duplicar los buffers R1/R2.
    #[must_use]
    pub fn into_observed(self) -> ObservedPair {
        self.observed
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum ErrorQualityError {
    InvalidErrorScale {
        mate: TechnicalMate,
        value: f64,
    },
    InvalidEffectivePhred {
        mate: TechnicalMate,
        cycle: u16,
        value: i32,
    },
    InvalidFinalProbability {
        mate: TechnicalMate,
        cycle: u16,
        value: f64,
    },
    InvalidFinalPhred {
        mate: TechnicalMate,
        cycle: u16,
        value: f64,
    },
    InvalidBase {
        mate: TechnicalMate,
        cycle: u16,
        base: u8,
    },
    CounterOverflow,
    Random(RandomSamplingError),
    Subseed(crate::diagnostics::Diagnostic),
}

impl fmt::Display for ErrorQualityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidErrorScale { .. } => f.write_str("escala de error F4 inválida"),
            Self::InvalidEffectivePhred { .. } => f.write_str("Phred efectivo F4 fuera de 0..=93"),
            Self::InvalidFinalProbability { .. } => {
                f.write_str("probabilidad final de error F4 inválida")
            }
            Self::InvalidFinalPhred { .. } => f.write_str("Phred final emitido F4 fuera de 0..=93"),
            Self::InvalidBase { .. } => f.write_str("base verdadera F4 inválida"),
            Self::CounterOverflow => f.write_str("contador agregado F4 desbordado"),
            Self::Random(error) => write!(f, "muestreo F4 inválido: {error:?}"),
            Self::Subseed(error) => write!(f, "subsemilla F4 inválida: {error:?}"),
        }
    }
}

impl std::error::Error for ErrorQualityError {}

/// Calidad efectiva por ciclo, sin RNG. La interpolación usa aritmética entera
/// y redondeo half-up del recurso; el desplazamiento se aplica después.
pub fn effective_phred_for_cycle(
    profile: F4TechnicalProfile,
    mate: TechnicalMate,
    cycle: u16,
    settings: ErrorQualitySettings,
) -> Result<u8, ErrorQualityError> {
    let curve = profile.quality_curve(profile_mate(mate));
    let base =
        interpolated_phred(curve, cycle).ok_or(ErrorQualityError::InvalidEffectivePhred {
            mate,
            cycle,
            value: -1,
        })?;
    let effective = i32::from(base) + i32::from(settings.quality_shift(mate));
    u8::try_from(effective)
        .ok()
        .filter(|value| *value <= 93)
        .ok_or(ErrorQualityError::InvalidEffectivePhred {
            mate,
            cycle,
            value: effective,
        })
}

/// Transforma una pareja verdadera bajo sus claves estables. R1/R2 y cada
/// ciclo reinician sus decisiones, por lo que solapamiento, orden de emisión y
/// una sustitución previa no cambian ningún otro resultado.
pub fn transform_paired_end(
    root_seed: RootSeed,
    event_id: &EventId,
    pair: &PairedEndPair,
    profile: F4TechnicalProfile,
    settings: ErrorQualitySettings,
) -> Result<ObservedPairTransformation, ErrorQualityError> {
    let template = pair.pre_hop_identity().template();
    let library_id = template.library_id();
    let ordinal = template.template_ordinal().get();
    let (r1, r1_counts) = transform_mate(
        root_seed,
        event_id,
        library_id,
        ordinal,
        pair.r1().mate(),
        pair.r1().sequence(),
        profile,
        settings,
    )?;
    let (r2, r2_counts) = transform_mate(
        root_seed,
        event_id,
        library_id,
        ordinal,
        pair.r2().mate(),
        pair.r2().sequence(),
        profile,
        settings,
    )?;
    Ok(ObservedPairTransformation {
        observed: ObservedPair { r1, r2 },
        report: ErrorQualityReport {
            r1: r1_counts,
            r2: r2_counts,
        },
    })
}

#[allow(clippy::too_many_arguments)]
fn transform_mate(
    root_seed: RootSeed,
    event_id: &EventId,
    library_id: &LibraryId,
    template_ordinal: u64,
    mate: TechnicalMate,
    true_sequence: &[u8],
    profile: F4TechnicalProfile,
    settings: ErrorQualitySettings,
) -> Result<(ObservedMate, MateErrorCounts), ErrorQualityError> {
    let mut sequence = Vec::with_capacity(true_sequence.len());
    let mut phred = Vec::with_capacity(true_sequence.len());
    let mut counts = MateErrorCounts::default();
    for (index, base) in true_sequence.iter().copied().enumerate() {
        let cycle = u16::try_from(index + 1).map_err(|_| ErrorQualityError::CounterOverflow)?;
        if base == b'N' {
            sequence.push(b'N');
            phred.push(0);
            continue;
        }
        if !matches!(base, b'A' | b'C' | b'G' | b'T') {
            return Err(ErrorQualityError::InvalidBase { mate, cycle, base });
        }
        let adjusted_phred = effective_phred_for_cycle(profile, mate, cycle, settings)?;
        counts.opportunities = counts
            .opportunities
            .checked_add(1)
            .ok_or(ErrorQualityError::CounterOverflow)?;
        let (observed, emitted_phred) = if settings.sequencing_errors()
            == SequencingErrorsMode::Disabled
        {
            // No se usa escala cero: se conservan literalmente bases y Phred
            // del perfil y no se deriva ninguna decisión aleatoria de error.
            (base, adjusted_phred)
        } else {
            let probability = final_error_probability(adjusted_phred, settings.error_scale(mate));
            let emitted_phred = emitted_phred_from_probability(probability, mate, cycle)?;
            let f4_mate = rng_mate(mate);
            let decision = F4TechnicalDecision::error_bernoulli(
                event_id,
                library_id,
                template_ordinal,
                f4_mate,
                u64::from(cycle),
            );
            let mut error_rng = derive_f4_technical_subseed(root_seed, decision)
                .map_err(ErrorQualityError::Subseed)?
                .local_rng();
            let errored =
                f4_bernoulli(&mut error_rng, probability).map_err(ErrorQualityError::Random)?;
            let observed = if errored {
                counts.errors = counts
                    .errors
                    .checked_add(1)
                    .ok_or(ErrorQualityError::CounterOverflow)?;
                let decision = F4TechnicalDecision::error_substitution(
                    event_id,
                    library_id,
                    template_ordinal,
                    f4_mate,
                    u64::from(cycle),
                );
                let mut substitution_rng = derive_f4_technical_subseed(root_seed, decision)
                    .map_err(ErrorQualityError::Subseed)?
                    .local_rng();
                alternate_base(
                    base,
                    f4_uniform_index(&mut substitution_rng, 3)
                        .map_err(ErrorQualityError::Random)?,
                )
            } else {
                base
            };
            (observed, emitted_phred)
        };
        sequence.push(observed);
        phred.push(emitted_phred);
    }
    Ok((
        ObservedMate {
            mate,
            sequence,
            phred,
        },
        counts,
    ))
}

fn interpolated_phred(curve: QualityCurve, cycle: u16) -> Option<u8> {
    let segment = curve
        .segments()
        .iter()
        .copied()
        .find(|segment| (segment.first_cycle()..=segment.last_cycle()).contains(&cycle))?;
    let span = i64::from(segment.last_cycle() - segment.first_cycle());
    if span == 0 {
        return Some(segment.first_phred());
    }
    let offset = i64::from(cycle - segment.first_cycle());
    let numerator = i64::from(segment.first_phred()) * span
        + (i64::from(segment.last_phred()) - i64::from(segment.first_phred())) * offset;
    let rounded = (2 * numerator + span) / (2 * span);
    u8::try_from(rounded).ok()
}

fn final_error_probability(phred: u8, scale: f64) -> f64 {
    libm::pow(10.0, -f64::from(phred) / 10.0) * scale
}

fn emitted_phred_from_probability(
    probability: f64,
    mate: TechnicalMate,
    cycle: u16,
) -> Result<u8, ErrorQualityError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) || probability == 0.0 {
        return Err(ErrorQualityError::InvalidFinalProbability {
            mate,
            cycle,
            value: probability,
        });
    }

    let phred = -10.0 * libm::log10(probability);
    let rounded = phred.round();
    if !rounded.is_finite() || !(0.0..=93.0).contains(&rounded) {
        return Err(ErrorQualityError::InvalidFinalPhred {
            mate,
            cycle,
            value: rounded,
        });
    }

    u8::try_from(rounded as i32).map_err(|_| ErrorQualityError::InvalidFinalPhred {
        mate,
        cycle,
        value: rounded,
    })
}

fn profile_mate(mate: TechnicalMate) -> ProfileMate {
    match mate {
        TechnicalMate::R1 => ProfileMate::R1,
        TechnicalMate::R2 => ProfileMate::R2,
    }
}

fn rng_mate(mate: TechnicalMate) -> F4Mate {
    match mate {
        TechnicalMate::R1 => F4Mate::R1,
        TechnicalMate::R2 => F4Mate::R2,
    }
}

fn alternate_base(base: u8, choice: u64) -> u8 {
    let alternatives = match base {
        b'A' => b"CGT",
        b'C' => b"AGT",
        b'G' => b"ACT",
        b'T' => b"ACG",
        _ => unreachable!("validated A/C/G/T base"),
    };
    alternatives[usize::try_from(choice).expect("uniform choice below three")]
}
