//! Admisión y configuración efectiva mínima de F1.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::configuration::{
    DeclaredGeneralConfig, DeclaredPopulationSnvs, DeclaredReductionTarget,
    DeclaredReferenceSource, ExpectedConfigKind, SequenceClass,
};
use crate::diagnostics::Diagnostic;
use crate::exact_decimal::ExactDecimal;
use crate::haplotype::Ploidy;
use crate::identity::{ReferenceId, SequenceId};
use crate::population_snvs::{
    FrequencyAssignmentMode, FrequencySampling, FrequencySpectrumComponent,
    GeneratedFrequencySpectrum, GeneratedLocusCount, GeneratedPopulationFrequency,
    GeneratedSpatialMode, GeneratedSubstitutionProfile, GeneratedWeightedRegion,
    PopulationSnvsGenerationConfig, validate_generation_config,
};
use crate::randomness::{
    EntropySource, ReproducibilityMetadata, RuntimeSubseedScope, Subseed, SystemEntropy,
    derive_subseed, resolve_root_seed,
};
use crate::reference_generation::normalize_composition_probabilities;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ValidatedGeneralSummary {
    pub total_genome_size: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveSequence {
    id: String,
    class: SequenceClass,
    length: u64,
}

impl EffectiveSequence {
    #[must_use]
    pub fn id(&self) -> &str {
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

#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveComposition {
    values: [f64; 4],
}

/// Parámetro efectivo compartido por las normalizaciones DTC-049 consumidas en F2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffectiveNormalization {
    tolerance: f64,
}

impl EffectiveNormalization {
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }
}

impl EffectiveComposition {
    #[must_use]
    pub const fn values(&self) -> &[f64; 4] {
        &self.values
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EffectiveReferenceSource {
    Simulated,
    Supplied { original_path: PathBuf },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EffectiveReductionTarget {
    Bases { declared: u64 },
    ReductionPercent { literal: String },
}

impl EffectiveReductionTarget {
    #[must_use]
    pub const fn mode(&self) -> &'static str {
        match self {
            Self::Bases { .. } => "target_bases",
            Self::ReductionPercent { .. } => "reduction_percent",
        }
    }

    #[must_use]
    pub fn declared_value(&self) -> String {
        match self {
            Self::Bases { declared } => declared.to_string(),
            Self::ReductionPercent { literal } => literal.clone(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveReducedReference {
    id: String,
    source_reference_id: String,
    region_size_bases: u64,
    target: EffectiveReductionTarget,
    target_bases: u64,
    selection_subseed: Subseed,
}

impl EffectiveReducedReference {
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn source_reference_id(&self) -> &str {
        &self.source_reference_id
    }

    #[must_use]
    pub const fn region_size_bases(&self) -> u64 {
        self.region_size_bases
    }

    #[must_use]
    pub const fn target(&self) -> &EffectiveReductionTarget {
        &self.target
    }

    #[must_use]
    pub const fn target_bases(&self) -> u64 {
        self.target_bases
    }

    #[must_use]
    pub const fn selection_subseed(&self) -> Subseed {
        self.selection_subseed
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveConfigurationInput {
    original_path: String,
    bytes: Arc<[u8]>,
}

impl EffectiveConfigurationInput {
    #[must_use]
    pub fn original_path(&self) -> &str {
        &self.original_path
    }

    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveConfig {
    source_path: String,
    reference_id: String,
    ploidy: Ploidy,
    normalization: EffectiveNormalization,
    reference_source: EffectiveReferenceSource,
    composition: Option<EffectiveComposition>,
    reduced_references: Vec<EffectiveReducedReference>,
    total_genome_size: u64,
    sequences: Vec<EffectiveSequence>,
    configuration_input: EffectiveConfigurationInput,
    reproducibility: ReproducibilityMetadata,
    population_snvs_generation: Option<PopulationSnvsGenerationConfig>,
}

impl EffectiveConfig {
    #[must_use]
    pub fn source_path(&self) -> &str {
        &self.source_path
    }

    #[must_use]
    pub fn reference_id(&self) -> &str {
        &self.reference_id
    }

    #[must_use]
    pub const fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    #[must_use]
    pub const fn normalization(&self) -> EffectiveNormalization {
        self.normalization
    }

    #[must_use]
    pub const fn reference_source(&self) -> &EffectiveReferenceSource {
        &self.reference_source
    }

    #[must_use]
    pub const fn composition(&self) -> Option<&EffectiveComposition> {
        self.composition.as_ref()
    }

    #[must_use]
    pub fn reduced_references(&self) -> &[EffectiveReducedReference] {
        &self.reduced_references
    }

    #[must_use]
    pub const fn total_genome_size(&self) -> u64 {
        self.total_genome_size
    }

    #[must_use]
    pub fn sequences(&self) -> &[EffectiveSequence] {
        &self.sequences
    }

    #[must_use]
    pub const fn configuration_input(&self) -> &EffectiveConfigurationInput {
        &self.configuration_input
    }

    #[must_use]
    pub const fn reproducibility(&self) -> &ReproducibilityMetadata {
        &self.reproducibility
    }
    #[must_use]
    pub fn population_snvs_generation(&self) -> Option<&PopulationSnvsGenerationConfig> {
        self.population_snvs_generation.as_ref()
    }
}

pub fn load_and_validate_config(path: &Path) -> Result<EffectiveConfig, Diagnostic> {
    load_and_validate_config_with_entropy(path, &mut SystemEntropy)
}

pub fn load_and_validate_config_with_entropy(
    path: &Path,
    entropy: &mut impl EntropySource,
) -> Result<EffectiveConfig, Diagnostic> {
    let bytes = std::fs::read(path).map_err(|error| input_error(path, error))?;
    let admitted = admit_toml_header(path, bytes, ExpectedConfigKind::General)?;
    let general = parse_declared_config(admitted)?;
    build_effective_config_with_entropy(&general, entropy)
}

pub fn validate_general_config(
    general: &DeclaredGeneralConfig,
) -> Result<ValidatedGeneralSummary, Diagnostic> {
    effective_ploidy(general)?;
    ReferenceId::parse(general.reference().id.value(), "configuration")?;
    if general.reference().sequences.is_empty() {
        return Err(config_error(
            general,
            "reference.sequences",
            "se exige al menos una secuencia",
        ));
    }
    let mut ids = BTreeSet::new();
    let mut total = 0_u64;
    for sequence in &general.reference().sequences {
        SequenceId::parse(sequence.id.value(), general.reference().id.value())?;
        if !ids.insert(sequence.id.value()) {
            return Err(config_error(
                general,
                "reference.sequences.id",
                format!("identificador repetido: {}", sequence.id.value()),
            ));
        }
        if *sequence.length.value() == 0 {
            return Err(config_error(
                general,
                "reference.sequences.length",
                format!("{} tiene longitud cero", sequence.id.value()),
            ));
        }
        total = total.checked_add(*sequence.length.value()).ok_or_else(|| {
            config_error(
                general,
                "reference.sequences.length",
                "el tamaño total del genoma desborda u64",
            )
        })?;
    }
    if let DeclaredReferenceSource::Simulated { composition } = &general.reference().source {
        normalize_composition(general, composition)?;
    }
    if !general.reduced_references().is_empty()
        && !matches!(
            &general.reference().source,
            DeclaredReferenceSource::Supplied { .. }
        )
    {
        return Err(config_error(
            general,
            "reduced_references",
            "F1 solo permite crear referencias reducidas desde una referencia suministrada",
        ));
    }
    let mut reference_ids = BTreeSet::from([general.reference().id.value().as_str()]);
    for reduced in general.reduced_references() {
        ReferenceId::parse(reduced.id.value(), "configuration")?;
        if !reference_ids.insert(reduced.id.value()) {
            return Err(config_error(
                general,
                "reduced_references.id",
                format!("identificador repetido: {}", reduced.id.value()),
            ));
        }
        if reduced.source_reference_id.value() != general.reference().id.value() {
            return Err(config_error(
                general,
                "reduced_references.source_reference_id",
                "F1 solo permite reducir la referencia inicial declarada",
            ));
        }
        if *reduced.region_size_bases.value() == 0 {
            return Err(config_error(
                general,
                "reduced_references.region_size_bases",
                "el tamaño de región debe ser positivo",
            ));
        }
    }
    Ok(ValidatedGeneralSummary {
        total_genome_size: total,
    })
}

pub fn build_effective_config(
    general: &DeclaredGeneralConfig,
) -> Result<EffectiveConfig, Diagnostic> {
    build_effective_config_with_entropy(general, &mut SystemEntropy)
}

pub fn build_effective_config_with_entropy(
    general: &DeclaredGeneralConfig,
    entropy: &mut impl EntropySource,
) -> Result<EffectiveConfig, Diagnostic> {
    let summary = validate_general_config(general)?;
    let ploidy = effective_ploidy(general)?;
    let normalization = EffectiveNormalization {
        tolerance: effective_normalization_tolerance(general)?,
    };
    let (root_seed, seed_origin) = resolve_root_seed(
        general.seed().map(|declared| declared.value().as_str()),
        entropy,
    )
    .map_err(|diagnostic| {
        let mut diagnostic = diagnostic.with("path", general.source().path().display().to_string());
        if let Some(seed) = general.seed() {
            diagnostic = diagnostic.with("location", seed.location().to_string());
        }
        diagnostic
    })?;
    let composition = match &general.reference().source {
        DeclaredReferenceSource::Simulated { composition } => Some(EffectiveComposition {
            values: normalize_composition(general, composition)?,
        }),
        DeclaredReferenceSource::Supplied { .. } => None,
    };
    let sequences: Vec<EffectiveSequence> = general
        .reference()
        .sequences
        .iter()
        .map(|sequence| EffectiveSequence {
            id: sequence.id.value().clone(),
            class: *sequence.class.value(),
            length: *sequence.length.value(),
        })
        .collect();
    let mut reduced_references = Vec::with_capacity(general.reduced_references().len());
    for reduction in general.reduced_references() {
        let (target, target_bases) = match &reduction.target {
            DeclaredReductionTarget::Bases(value) => (
                EffectiveReductionTarget::Bases {
                    declared: *value.value(),
                },
                *value.value(),
            ),
            DeclaredReductionTarget::ReductionPercent(value) => {
                let exact = ExactDecimal::parse("reduction_percent", value.literal()).map_err(
                    |diagnostic| {
                        diagnostic
                            .with("path", general.source().path().display().to_string())
                            .with("location", value.location().to_string())
                    },
                )?;
                (
                    EffectiveReductionTarget::ReductionPercent {
                        literal: value.literal().to_owned(),
                    },
                    exact
                        .retained_after_percent(summary.total_genome_size)
                        .map_err(|diagnostic| {
                            diagnostic
                                .with("path", general.source().path().display().to_string())
                                .with("location", value.location().to_string())
                        })?,
                )
            }
        };
        if target_bases == 0 || target_bases >= summary.total_genome_size {
            return Err(Diagnostic::configuration(
                "GPS-CONF-021",
                "El objetivo de la referencia reducida no es menor que el espacio elegible.",
            )
            .with("path", general.source().path().display().to_string())
            .with("reference", reduction.id.value())
            .with("target_bases", target_bases.to_string())
            .with("eligible_bases", summary.total_genome_size.to_string()));
        }
        reduced_references.push(EffectiveReducedReference {
            id: reduction.id.value().clone(),
            source_reference_id: reduction.source_reference_id.value().clone(),
            region_size_bases: *reduction.region_size_bases.value(),
            target,
            target_bases,
            selection_subseed: derive_subseed(
                root_seed,
                RuntimeSubseedScope::ReferenceReduction {
                    source_reference_id: reduction.source_reference_id.value(),
                    derived_reference_id: reduction.id.value(),
                },
            )?,
        });
    }
    let reference_source = match &general.reference().source {
        DeclaredReferenceSource::Simulated { .. } => EffectiveReferenceSource::Simulated,
        DeclaredReferenceSource::Supplied { fasta } => EffectiveReferenceSource::Supplied {
            original_path: resolve_link(general.source().path(), fasta.value()),
        },
    };
    let reproducibility = ReproducibilityMetadata::new(root_seed, seed_origin);
    let population_snvs_generation = effective_population_snvs_generation(general)?;
    if let Some(generation) = &population_snvs_generation {
        validate_generation_config(generation).map_err(|error| {
            config_error(general, "population_snvs.generate", error.to_string())
        })?;
    }
    if let Some(PopulationSnvsGenerationConfig {
        spatial: GeneratedSpatialMode::WeightedRegions { regions, .. },
        ..
    }) = &population_snvs_generation
    {
        for region in regions {
            let valid = sequences
                .iter()
                .find(|sequence| sequence.id() == region.contig)
                .is_some_and(|sequence| {
                    region.start_1based > 0 && region.end_1based <= sequence.length()
                });
            if !valid {
                return Err(config_error(
                    general,
                    "population_snvs.generate.regions",
                    format!(
                        "región fuera de la referencia efectiva: {}:{}..{}",
                        region.contig, region.start_1based, region.end_1based
                    ),
                ));
            }
        }
    }
    Ok(EffectiveConfig {
        source_path: general.source().path().display().to_string(),
        reference_id: general.reference().id.value().clone(),
        ploidy,
        normalization,
        reference_source,
        composition,
        reduced_references,
        total_genome_size: summary.total_genome_size,
        sequences,
        configuration_input: EffectiveConfigurationInput {
            original_path: general.source().path().display().to_string(),
            bytes: Arc::from(general.source().bytes()),
        },
        reproducibility,
        population_snvs_generation,
    })
}

fn effective_population_snvs_generation(
    general: &DeclaredGeneralConfig,
) -> Result<Option<PopulationSnvsGenerationConfig>, Diagnostic> {
    let Some(DeclaredPopulationSnvs::Generate { parameters }) = general.population_snvs() else {
        return Ok(None);
    };
    let decimal = |field: &str, value: &crate::configuration::DeclaredNumber| {
        value
            .literal()
            .parse::<f64>()
            .map_err(|_| config_error(general, field, "debe ser un f64 válido"))
    };
    let count = if let Some(value) = &parameters.locus_count {
        GeneratedLocusCount::Exact(*value.value())
    } else if let Some(value) = &parameters.mean_loci {
        GeneratedLocusCount::PoissonMean(decimal("population_snvs.generate.mean_loci", value)?)
    } else {
        GeneratedLocusCount::PoissonDensityPerBase(decimal(
            "population_snvs.generate.mean_density_per_base",
            parameters.mean_density_per_base.as_ref().ok_or_else(|| {
                config_error(general, "population_snvs.generate", "recuento ausente")
            })?,
        )?)
    };
    let spatial = if parameters.spatial_mode.value() == "uniform" {
        GeneratedSpatialMode::Uniform
    } else {
        GeneratedSpatialMode::WeightedRegions {
            outside_weight: decimal(
                "population_snvs.generate.outside_weight",
                parameters.outside_weight.as_ref().ok_or_else(|| {
                    config_error(
                        general,
                        "population_snvs.generate.outside_weight",
                        "falta outside_weight",
                    )
                })?,
            )?,
            regions: parameters
                .regions
                .iter()
                .map(|region| {
                    Ok(GeneratedWeightedRegion {
                        contig: region.contig.value().clone(),
                        start_1based: *region.start.value(),
                        end_1based: *region.end.value(),
                        weight: decimal("population_snvs.generate.regions.weight", &region.weight)?,
                    })
                })
                .collect::<Result<Vec<_>, Diagnostic>>()?,
        }
    };
    Ok(Some(PopulationSnvsGenerationConfig {
        count,
        spatial,
        substitution_profile: match parameters.substitution_profile.value().as_str() {
            "uniform_snv12_v1" => GeneratedSubstitutionProfile::UniformSnv12V1,
            "human_harris2017_snv12_approx_v1" => {
                GeneratedSubstitutionProfile::HumanHarris2017Snv12ApproxV1
            }
            _ => return Err(config_error(general, "population_snvs.generate.substitution_profile", "perfil de sustitución desconocido")),
        },
        populations: parameters
            .populations
            .iter()
            .map(|population| {
                let spectrum = match population.frequency_model.value().as_str() {
                    "uniform" => None,
                    "human_1000g_global_sfs_approx_v1" => Some(human_1000g_spectrum(
                        population.frequency_mode.as_ref().map_or("stochastic", |mode| mode.value()),
                    )?),
                    "spectrum" => Some(GeneratedFrequencySpectrum {
                        mode: frequency_assignment_mode(
                            population.frequency_mode.as_ref().map_or("stochastic", |mode| mode.value()),
                        )?,
                        components: population.spectrum.iter().map(|component| {
                            let weight = decimal("population_snvs.generate.populations.spectrum.weight", &component.weight)?;
                            match component.kind.value().as_str() {
                                "fixed" => Ok(FrequencySpectrumComponent::Fixed {
                                    value: decimal("population_snvs.generate.populations.spectrum.fixed_alt_af", component.fixed_alt_af.as_ref().ok_or_else(|| config_error(general, "population_snvs.generate.populations.spectrum.fixed_alt_af", "fixed exige fixed_alt_af"))?)?,
                                    weight,
                                }),
                                "interval" => Ok(FrequencySpectrumComponent::Interval {
                                    min: decimal("population_snvs.generate.populations.spectrum.min_alt_af", component.min_alt_af.as_ref().ok_or_else(|| config_error(general, "population_snvs.generate.populations.spectrum.min_alt_af", "interval exige min_alt_af"))?)?,
                                    max: decimal("population_snvs.generate.populations.spectrum.max_alt_af", component.max_alt_af.as_ref().ok_or_else(|| config_error(general, "population_snvs.generate.populations.spectrum.max_alt_af", "interval exige max_alt_af"))?)?,
                                    weight,
                                    sampling: match component.sampling.as_ref().map_or("", |sampling| sampling.value().as_str()) {
                                        "linear_uniform" => FrequencySampling::LinearUniform,
                                        "log_uniform" => FrequencySampling::LogUniform,
                                        _ => return Err(config_error(general, "population_snvs.generate.populations.spectrum.sampling", "interval exige linear_uniform o log_uniform")),
                                    },
                                }),
                                _ => Err(config_error(general, "population_snvs.generate.populations.spectrum.kind", "debe ser fixed o interval")),
                            }
                        }).collect::<Result<Vec<_>, Diagnostic>>()?,
                    }),
                    _ => return Err(config_error(general, "population_snvs.generate.populations.frequency_model", "modelo de frecuencia desconocido")),
                };
                Ok(GeneratedPopulationFrequency {
                    population_id: population.population_id.value().clone(),
                    min_alt_af: population.min_alt_af.as_ref().map(|value| decimal(
                        "population_snvs.generate.populations.min_alt_af",
                        value,
                    )).transpose()?.unwrap_or(0.0),
                    max_alt_af: population.max_alt_af.as_ref().map(|value| decimal(
                        "population_snvs.generate.populations.max_alt_af",
                        value,
                    )).transpose()?.unwrap_or(1.0),
                    spectrum,
                })
            })
            .collect::<Result<Vec<_>, Diagnostic>>()?,
    }))
}

fn frequency_assignment_mode(value: &str) -> Result<FrequencyAssignmentMode, Diagnostic> {
    match value {
        "stochastic" => Ok(FrequencyAssignmentMode::Stochastic),
        "exact" => Ok(FrequencyAssignmentMode::Exact),
        _ => Err(Diagnostic::configuration(
            "GPS-CONF-008",
            "El modo de espectro de frecuencias no es válido.",
        )),
    }
}

fn human_1000g_spectrum(mode: &str) -> Result<GeneratedFrequencySpectrum, Diagnostic> {
    Ok(GeneratedFrequencySpectrum {
        mode: frequency_assignment_mode(mode)?,
        components: vec![
            FrequencySpectrumComponent::Interval {
                min: 1.0 / 5008.0,
                max: 0.005,
                weight: 64.0,
                sampling: FrequencySampling::LogUniform,
            },
            FrequencySpectrumComponent::Interval {
                min: 0.005,
                max: 0.05,
                weight: 12.0,
                sampling: FrequencySampling::LogUniform,
            },
            FrequencySpectrumComponent::Interval {
                min: 0.05,
                max: 1.0,
                weight: 8.0,
                sampling: FrequencySampling::LinearUniform,
            },
        ],
    })
}

fn effective_normalization_tolerance(general: &DeclaredGeneralConfig) -> Result<f64, Diagnostic> {
    const DEFAULT: f64 = 1e-9;
    const MAXIMUM: f64 = 1e-6;
    let Some(declared) = general.normalization_tolerance() else {
        return Ok(DEFAULT);
    };
    let value = declared
        .literal()
        .parse::<f64>()
        .map_err(|_| config_error(general, "normalization_tolerance", "debe ser un f64 válido"))?;
    if !value.is_finite() || !(0.0..=MAXIMUM).contains(&value) {
        return Err(config_error(
            general,
            "normalization_tolerance",
            "debe ser finita y pertenecer a [0,1e-6]",
        ));
    }
    Ok(value)
}

fn effective_ploidy(general: &DeclaredGeneralConfig) -> Result<Ploidy, Diagnostic> {
    let Some(declared) = general.ploidy() else {
        return Ok(Ploidy::DIPLOID);
    };
    let value = *declared.value();
    let value = u8::try_from(value).map_err(|_| {
        config_error(
            general,
            "ploidy",
            format!("debe pertenecer a {}..={}", Ploidy::MIN, Ploidy::MAX),
        )
    })?;
    Ploidy::new(value).map_err(|_| {
        config_error(
            general,
            "ploidy",
            format!("debe pertenecer a {}..={}", Ploidy::MIN, Ploidy::MAX),
        )
    })
}

fn normalize_composition(
    general: &DeclaredGeneralConfig,
    declared: &[crate::configuration::DeclaredCompositionValue],
) -> Result<[f64; 4], Diagnostic> {
    if declared.len() != 4 {
        return Err(config_error(
            general,
            "reference.composition",
            "se exigen A, C, G y T exactamente una vez",
        ));
    }
    let mut parsed = [0.0_f64; 4];
    let mut seen = [false; 4];
    for item in declared {
        let index = match item.base.value().as_str() {
            "A" => 0,
            "C" => 1,
            "G" => 2,
            "T" => 3,
            _ => {
                return Err(config_error(
                    general,
                    "reference.composition",
                    format!("base no admitida: {}", item.base.value()),
                ));
            }
        };
        if std::mem::replace(&mut seen[index], true) {
            return Err(config_error(
                general,
                "reference.composition",
                format!("base repetida: {}", item.base.value()),
            ));
        }
        let value = item.value.literal().parse::<f64>().map_err(|_| {
            config_error(
                general,
                "reference.composition",
                format!("valor decimal inválido para {}", item.base.value()),
            )
        })?;
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(config_error(
                general,
                "reference.composition",
                format!("{} debe ser finito y pertenecer a [0,1]", item.base.value()),
            ));
        }
        parsed[index] = value;
    }
    normalize_composition_probabilities(parsed)
        .map_err(|detail| config_error(general, "reference.composition", detail))
}

fn resolve_link(general_path: &Path, declared: &str) -> PathBuf {
    let declared = Path::new(declared);
    if declared.is_absolute() {
        declared.to_path_buf()
    } else {
        general_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(declared)
    }
}

fn input_error(path: &Path, error: std::io::Error) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-001",
        "No se puede abrir o leer una entrada declarada.",
    )
    .with("path", path.display().to_string())
    .with("detail", error.to_string())
}

fn config_error(
    general: &DeclaredGeneralConfig,
    fields: &str,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", general.source().path().display().to_string())
    .with("fields", fields)
    .with("detail", detail.into())
}

#[cfg(test)]
mod profile_resource_tests {
    use super::{
        FrequencyAssignmentMode, FrequencySampling, FrequencySpectrumComponent,
        human_1000g_spectrum,
    };

    #[test]
    fn human_1000g_spectrum_matches_its_versioned_resource_literals() {
        let resource = include_str!("../resources/f2_2_6/human_1000g_global_sfs_approx_v1.tsv");
        let rows = resource
            .lines()
            .skip(1)
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 3);
        let spectrum = human_1000g_spectrum("exact").unwrap();
        assert_eq!(spectrum.mode, FrequencyAssignmentMode::Exact);
        for (component, fields) in spectrum.components.iter().zip(rows) {
            assert_eq!(fields[0], "interval");
            let (min, max, weight, sampling) = match component {
                FrequencySpectrumComponent::Interval {
                    min,
                    max,
                    weight,
                    sampling,
                } => (min, max, weight, sampling),
                FrequencySpectrumComponent::Fixed { .. } => {
                    panic!("the SFS resource has only intervals")
                }
            };
            assert_eq!(min.to_bits(), fields[1].parse::<f64>().unwrap().to_bits());
            assert_eq!(max.to_bits(), fields[2].parse::<f64>().unwrap().to_bits());
            assert_eq!(
                weight.to_bits(),
                fields[3].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(
                *sampling,
                match fields[4] {
                    "log_uniform" => FrequencySampling::LogUniform,
                    "linear_uniform" => FrequencySampling::LinearUniform,
                    _ => panic!("unexpected versioned SFS sampling literal"),
                }
            );
        }
    }
}
