//! Admisión y configuración efectiva mínima de F1.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::configuration::{
    DeclaredGeneralConfig, DeclaredPopulationSnvs, DeclaredReductionTarget,
    DeclaredReferenceSource, DeclaredTechnicalEventTarget, ExpectedConfigKind, SequenceClass,
};
use crate::diagnostics::Diagnostic;
use crate::exact_decimal::ExactDecimal;
use crate::haplotype::{HaplotypePhase, Ploidy};
use crate::identity::{
    AutomaticIdRequest, EventId, FounderKey, IdClass, IdentityRegistry, IndividualId, PopulationId,
    ReferenceId, SequenceId,
};
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
use crate::vcf_contract::VcfOutputMode;

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
    unknown_base_fraction: f64,
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

    /// Fracción independiente de posiciones ambiguas (`N`) de la referencia generada.
    #[must_use]
    pub const fn unknown_base_fraction(&self) -> f64 {
        self.unknown_base_fraction
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
    population_snvs: Option<EffectivePopulationSnvs>,
    initial_populations: Vec<EffectiveInitialPopulation>,
    direct_events: Vec<EffectiveDirectEvent>,
    technical_event: Option<EffectiveTechnicalEvent>,
    vcf_output: Option<VcfOutputMode>,
}

/// Perfil sintético único que F3 genera antes de que F4 incorpore perfiles
/// cuantitativos o tecnologías alternativas.
pub const ILLUMINA_SE_MINIMAL_V1: &str = "illumina-se-minimal-v1";
pub const ILLUMINA_SE_MINIMAL_READ_LENGTH: u16 = 150;
pub const ILLUMINA_SE_MINIMAL_ACGT_PHRED: u8 = 40;
pub const ILLUMINA_SE_MINIMAL_N_PHRED: u8 = 0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EffectiveTechnicalEventTarget {
    CoverageWgs(f64),
    ReadCount(u64),
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveTechnicalEvent {
    event_id: EventId,
    individual_ids: Vec<IndividualId>,
    target: EffectiveTechnicalEventTarget,
}

impl EffectiveTechnicalEvent {
    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        &self.event_id
    }

    #[must_use]
    pub fn individual_ids(&self) -> &[IndividualId] {
        &self.individual_ids
    }

    #[must_use]
    pub const fn target(&self) -> EffectiveTechnicalEventTarget {
        self.target
    }

    #[must_use]
    pub const fn profile(&self) -> &'static str {
        ILLUMINA_SE_MINIMAL_V1
    }

    #[must_use]
    pub const fn read_length(&self) -> u16 {
        ILLUMINA_SE_MINIMAL_READ_LENGTH
    }

    #[must_use]
    pub const fn acgt_phred(&self) -> u8 {
        ILLUMINA_SE_MINIMAL_ACGT_PHRED
    }

    #[must_use]
    pub const fn n_phred(&self) -> u8 {
        ILLUMINA_SE_MINIMAL_N_PHRED
    }
}

/// Origen efectivo del recurso canónico F2.2.7.
#[derive(Clone, Debug, PartialEq)]
pub enum EffectivePopulationSnvs {
    Supplied {
        input_path: PathBuf,
    },
    Reuse {
        metadata_path: PathBuf,
    },
    Generate {
        parameters: PopulationSnvsGenerationConfig,
    },
}

/// Modelo de genotipo solicitado para una cohorte directa.
///
/// F2.3.1 solo conserva esta elección validada. El sorteo de composición y la
/// fase se incorporan en F2.3.2 y F2.3.3, respectivamente.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EffectiveInitialGenotypeModel {
    HardyWeinberg,
    FDiploid { inbreeding_coefficient: f64 },
}

/// Cambio directo de F2.3.4 ya ligado a una identidad individual y a una
/// fase existente. Las coordenadas de indels mantienen el convenio actual
/// del compositor, pues se evalúan en orden declarado tras los SNP de matriz.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EffectiveDirectEvent {
    Snp {
        event_id: EventId,
        individual_id: IndividualId,
        phase: HaplotypePhase,
        sequence_id: SequenceId,
        position_1based: u64,
        reference: u8,
        alternate: u8,
    },
    Insertion {
        event_id: EventId,
        individual_id: IndividualId,
        phase: HaplotypePhase,
        sequence_id: SequenceId,
        point: u64,
        bases: Arc<[u8]>,
    },
    Deletion {
        event_id: EventId,
        individual_id: IndividualId,
        phase: HaplotypePhase,
        sequence_id: SequenceId,
        start: u64,
        end: u64,
    },
}

impl EffectiveDirectEvent {
    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        match self {
            Self::Snp { event_id, .. }
            | Self::Insertion { event_id, .. }
            | Self::Deletion { event_id, .. } => event_id,
        }
    }

    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        match self {
            Self::Snp { individual_id, .. }
            | Self::Insertion { individual_id, .. }
            | Self::Deletion { individual_id, .. } => individual_id,
        }
    }

    #[must_use]
    pub const fn phase(&self) -> HaplotypePhase {
        match self {
            Self::Snp { phase, .. }
            | Self::Insertion { phase, .. }
            | Self::Deletion { phase, .. } => *phase,
        }
    }

    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        match self {
            Self::Snp { sequence_id, .. }
            | Self::Insertion { sequence_id, .. }
            | Self::Deletion { sequence_id, .. } => sequence_id,
        }
    }
}

/// Identidades y parámetros ya resueltos de una cohorte directa.
///
/// La población de origen es un perfil de frecuencias de la matriz, no una
/// población evolutiva ni una estructura reproductiva.
#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveInitialPopulation {
    population_id: PopulationId,
    individual_ids: Vec<IndividualId>,
    founder_keys: Option<Vec<FounderKey>>,
    ploidy: Ploidy,
    genotype_model: EffectiveInitialGenotypeModel,
}

impl EffectiveInitialPopulation {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }

    #[must_use]
    pub fn individual_ids(&self) -> &[IndividualId] {
        &self.individual_ids
    }

    #[must_use]
    pub fn founder_keys(&self) -> Option<&[FounderKey]> {
        self.founder_keys.as_deref()
    }

    #[must_use]
    pub const fn individual_count(&self) -> usize {
        self.individual_ids.len()
    }

    #[must_use]
    pub const fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    #[must_use]
    pub const fn genotype_model(&self) -> EffectiveInitialGenotypeModel {
        self.genotype_model
    }
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
        match self.population_snvs.as_ref() {
            Some(EffectivePopulationSnvs::Generate { parameters }) => Some(parameters),
            _ => None,
        }
    }

    #[must_use]
    pub fn population_snvs(&self) -> Option<&EffectivePopulationSnvs> {
        self.population_snvs.as_ref()
    }

    /// Cohortes directas listas para que F2.3.2 consuma sus identidades.
    #[must_use]
    pub fn initial_populations(&self) -> &[EffectiveInitialPopulation] {
        &self.initial_populations
    }

    /// Eventos directos tipados y en su orden de declaración global.
    #[must_use]
    pub fn direct_events(&self) -> &[EffectiveDirectEvent] {
        &self.direct_events
    }

    #[must_use]
    pub const fn technical_event(&self) -> Option<&EffectiveTechnicalEvent> {
        self.technical_event.as_ref()
    }

    #[must_use]
    pub const fn vcf_output(&self) -> Option<VcfOutputMode> {
        self.vcf_output
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
    if let DeclaredReferenceSource::Simulated {
        composition,
        unknown_base_fraction,
    } = &general.reference().source
    {
        normalize_composition(general, composition)?;
        parse_unknown_base_fraction(general, unknown_base_fraction)?;
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
        DeclaredReferenceSource::Simulated {
            composition,
            unknown_base_fraction,
        } => Some(EffectiveComposition {
            values: normalize_composition(general, composition)?,
            unknown_base_fraction: parse_unknown_base_fraction(general, unknown_base_fraction)?,
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
    let initial_populations = effective_initial_populations(general, root_seed)?;
    let population_snvs = effective_population_snvs(general)?;
    let direct_events = effective_direct_events(general, &initial_populations, &sequences)?;
    let technical_event = effective_technical_event(general, &initial_populations, &direct_events)?;
    let vcf_output = general
        .vcf_output()
        .map(|output| match output.mode.value().as_str() {
            "text" => VcfOutputMode::Text,
            "bgzf" => VcfOutputMode::Bgzf,
            _ => unreachable!("the TOML adapter validates vcf_output.mode"),
        });
    if !direct_events.is_empty() && population_snvs.is_none() {
        return Err(config_error(
            general,
            "direct_events,population_snvs",
            "direct_events exige una matriz SNV poblacional para ligar genotipo y fase",
        ));
    }
    if let Some(EffectivePopulationSnvs::Generate {
        parameters: generation,
    }) = &population_snvs
    {
        validate_generation_config(generation).map_err(|error| {
            config_error(general, "population_snvs.generate", error.to_string())
        })?;
    }
    if let Some(PopulationSnvsGenerationConfig {
        spatial: GeneratedSpatialMode::WeightedRegions { regions, .. },
        ..
    }) = population_snvs_generation_from_effective(population_snvs.as_ref())
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
        population_snvs,
        initial_populations,
        direct_events,
        technical_event,
        vcf_output,
    })
}

fn effective_technical_event(
    general: &DeclaredGeneralConfig,
    populations: &[EffectiveInitialPopulation],
    direct_events: &[EffectiveDirectEvent],
) -> Result<Option<EffectiveTechnicalEvent>, Diagnostic> {
    let Some(declared) = general.technical_event() else {
        return Ok(None);
    };
    let event_id = EventId::parse(declared.event_id.value().clone(), "technical_event")
        .map_err(|error| with_initial_population_path(general, error))?;
    if direct_events
        .iter()
        .any(|event| event.event_id() == &event_id)
    {
        return Err(config_error(
            general,
            "technical_event.event_id",
            format!("identificador de evento repetido en el ámbito de ejecución: {event_id}"),
        ));
    }
    let declared_individuals = declared.individual_ids.as_slice();
    let mut known_individuals = BTreeSet::new();
    for population in populations {
        known_individuals.extend(population.individual_ids().iter().cloned());
    }
    let mut individual_ids = Vec::with_capacity(declared_individuals.len());
    let mut selected = BTreeSet::new();
    for declared_id in declared_individuals {
        let individual_id = IndividualId::parse(declared_id.value().clone(), "technical_event")
            .map_err(|error| with_initial_population_path(general, error))?;
        if !known_individuals.contains(&individual_id) {
            return Err(config_error(
                general,
                "technical_event.individual_ids",
                format!("individuo no declarado en initial_populations: {individual_id}"),
            ));
        }
        if !selected.insert(individual_id.clone()) {
            return Err(config_error(
                general,
                "technical_event.individual_ids",
                format!("individuo repetido: {individual_id}"),
            ));
        }
        individual_ids.push(individual_id);
    }
    let target = match &declared.target {
        DeclaredTechnicalEventTarget::CoverageWgs(value) => {
            let coverage = value.literal().parse::<f64>().map_err(|_| {
                config_error(
                    general,
                    "technical_event.target_coverage_wgs",
                    "debe ser un f64 válido",
                )
            })?;
            if !coverage.is_finite() || coverage < 0.0 {
                return Err(config_error(
                    general,
                    "technical_event.target_coverage_wgs",
                    "debe ser finita y mayor o igual que cero",
                ));
            }
            EffectiveTechnicalEventTarget::CoverageWgs(coverage)
        }
        DeclaredTechnicalEventTarget::ReadCount(value) => {
            EffectiveTechnicalEventTarget::ReadCount(*value.value())
        }
    };
    Ok(Some(EffectiveTechnicalEvent {
        event_id,
        individual_ids,
        target,
    }))
}

fn effective_direct_events(
    general: &DeclaredGeneralConfig,
    populations: &[EffectiveInitialPopulation],
    sequences: &[EffectiveSequence],
) -> Result<Vec<EffectiveDirectEvent>, Diagnostic> {
    let mut individuals = std::collections::BTreeMap::new();
    for population in populations {
        for individual in population.individual_ids() {
            individuals.insert(individual.clone(), population.ploidy());
        }
    }
    let mut registry = IdentityRegistry::default();
    general
        .direct_events()
        .iter()
        .map(|declared| {
            let event_id = EventId::parse(declared.event_id.value().clone(), "direct_events")
                .map_err(|error| with_initial_population_path(general, error))?;
            let event_id = registry
                .register_event("execution", event_id)
                .map_err(|error| with_initial_population_path(general, error))?;
            let individual_id =
                IndividualId::parse(declared.individual_id.value().clone(), "direct_events")
                    .map_err(|error| with_initial_population_path(general, error))?;
            let ploidy = individuals.get(&individual_id).ok_or_else(|| {
                config_error(
                    general,
                    "direct_events.individual_id",
                    format!("individuo no declarado en initial_populations: {individual_id}"),
                )
            })?;
            let phase_index = u8::try_from(*declared.phase.value()).map_err(|_| {
                config_error(general, "direct_events.phase", "no se puede representar")
            })?;
            if phase_index >= ploidy.copy_count() {
                return Err(config_error(
                    general,
                    "direct_events.phase",
                    format!(
                        "la fase {phase_index} no existe para ploidía {}",
                        ploidy.copy_count()
                    ),
                ));
            }
            let phase = HaplotypePhase::from_copy_index(phase_index)
                .map_err(|_| config_error(general, "direct_events.phase", "fase inválida"))?;
            let sequence_id = SequenceId::parse(declared.contig.value().clone(), "direct_events")
                .map_err(|error| with_initial_population_path(general, error))?;
            if !sequences
                .iter()
                .any(|sequence| sequence.id() == sequence_id.as_str())
            {
                return Err(config_error(
                    general,
                    "direct_events.contig",
                    format!("contig ausente de la referencia efectiva: {sequence_id}"),
                ));
            }
            match declared.kind.value().as_str() {
                "snp" => {
                    let position_1based = *declared
                        .position_1based
                        .as_ref()
                        .expect("parser checked")
                        .value();
                    let reference = single_canonical_base(
                        general,
                        "direct_events.ref",
                        declared.reference.as_ref().expect("parser checked").value(),
                    )?;
                    let alternate = single_canonical_base(
                        general,
                        "direct_events.alt",
                        declared.alternate.as_ref().expect("parser checked").value(),
                    )?;
                    if reference == alternate {
                        return Err(config_error(
                            general,
                            "direct_events.ref,direct_events.alt",
                            "deben ser bases canónicas distintas",
                        ));
                    }
                    Ok(EffectiveDirectEvent::Snp {
                        event_id,
                        individual_id,
                        phase,
                        sequence_id,
                        position_1based,
                        reference,
                        alternate,
                    })
                }
                "insertion" => {
                    let bases = declared
                        .bases
                        .as_ref()
                        .expect("parser checked")
                        .value()
                        .as_bytes();
                    if bases.is_empty()
                        || bases.iter().any(|base| {
                            !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
                        })
                    {
                        return Err(config_error(
                            general,
                            "direct_events.bases",
                            "debe contener una o más bases A/C/G/T",
                        ));
                    }
                    Ok(EffectiveDirectEvent::Insertion {
                        event_id,
                        individual_id,
                        phase,
                        sequence_id,
                        point: *declared.point.as_ref().expect("parser checked").value(),
                        bases: Arc::from(
                            bases.iter().map(u8::to_ascii_uppercase).collect::<Vec<_>>(),
                        ),
                    })
                }
                "deletion" => Ok(EffectiveDirectEvent::Deletion {
                    event_id,
                    individual_id,
                    phase,
                    sequence_id,
                    start: *declared.start.as_ref().expect("parser checked").value(),
                    end: *declared.end.as_ref().expect("parser checked").value(),
                }),
                _ => unreachable!("parser validates direct event kind"),
            }
        })
        .collect()
}

fn single_canonical_base(
    general: &DeclaredGeneralConfig,
    field: &str,
    value: &str,
) -> Result<u8, Diagnostic> {
    let bytes = value.as_bytes();
    if bytes.len() != 1 || !matches!(bytes[0].to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') {
        return Err(config_error(
            general,
            field,
            "debe ser exactamente una base A/C/G/T",
        ));
    }
    Ok(bytes[0].to_ascii_uppercase())
}

fn effective_initial_populations(
    general: &DeclaredGeneralConfig,
    root_seed: crate::randomness::RootSeed,
) -> Result<Vec<EffectiveInitialPopulation>, Diagnostic> {
    let mut population_ids = BTreeSet::new();
    let mut derivation_registry = IdentityRegistry::default();
    let mut execution_registry = IdentityRegistry::default();
    let mut effective = Vec::with_capacity(general.initial_populations().len());

    for declared in general.initial_populations() {
        let population_id = PopulationId::parse(
            declared.population_id.value().clone(),
            "initial_populations",
        )
        .map_err(|error| with_initial_population_path(general, error))?;
        if !population_ids.insert(population_id.clone()) {
            return Err(config_error(
                general,
                "initial_populations.population_id",
                format!("identificador repetido: {population_id}"),
            ));
        }

        let individual_count = *declared.individual_count.value();
        if individual_count == 0 {
            return Err(config_error(
                general,
                "initial_populations.individual_count",
                "debe ser mayor que cero",
            ));
        }
        let individual_count = usize::try_from(individual_count).map_err(|_| {
            config_error(
                general,
                "initial_populations.individual_count",
                "no se puede representar en esta plataforma",
            )
        })?;
        let ploidy = u8::try_from(*declared.ploidy.value())
            .ok()
            .and_then(|value| Ploidy::new(value).ok())
            .ok_or_else(|| {
                config_error(
                    general,
                    "initial_populations.ploidy",
                    format!("debe pertenecer a {}..={}", Ploidy::MIN, Ploidy::MAX),
                )
            })?;
        let genotype_model = match declared.genotype_model.value().as_str() {
            "hwe" => EffectiveInitialGenotypeModel::HardyWeinberg,
            "f_diploid" => {
                if ploidy != Ploidy::DIPLOID {
                    return Err(config_error(
                        general,
                        "initial_populations.genotype_model,initial_populations.ploidy",
                        "f_diploid solo se admite con ploidía 2",
                    ));
                }
                let coefficient = declared
                    .inbreeding_coefficient
                    .as_ref()
                    .ok_or_else(|| {
                        config_error(
                            general,
                            "initial_populations.inbreeding_coefficient",
                            "f_diploid exige inbreeding_coefficient",
                        )
                    })?
                    .literal()
                    .parse::<f64>()
                    .map_err(|_| {
                        config_error(
                            general,
                            "initial_populations.inbreeding_coefficient",
                            "debe ser un f64 válido",
                        )
                    })?;
                if !coefficient.is_finite() || !(0.0..=1.0).contains(&coefficient) {
                    return Err(config_error(
                        general,
                        "initial_populations.inbreeding_coefficient",
                        "debe ser finito y pertenecer a [0,1]",
                    ));
                }
                EffectiveInitialGenotypeModel::FDiploid {
                    inbreeding_coefficient: coefficient,
                }
            }
            _ => unreachable!("parser validates initial genotype model"),
        };

        let (individual_ids, founder_keys) = match (
            &declared.individual_ids,
            &declared.founder_keys,
        ) {
            (Some(provided), None) => {
                if provided.len() != individual_count {
                    return Err(config_error(
                        general,
                        "initial_populations.individual_ids,initial_populations.individual_count",
                        "individual_ids debe contener exactamente individual_count identificadores",
                    ));
                }
                let individual_ids = provided
                    .iter()
                    .map(|id| {
                        let id = IndividualId::parse(id.value().clone(), "initial_populations")
                            .map_err(|error| with_initial_population_path(general, error))?;
                        execution_registry
                            .register_individual("execution", id)
                            .map_err(|error| with_initial_population_path(general, error))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                (individual_ids, None)
            }
            (None, Some(provided)) => {
                if provided.len() != individual_count {
                    return Err(config_error(
                        general,
                        "initial_populations.founder_keys,initial_populations.individual_count",
                        "founder_keys debe contener exactamente individual_count claves",
                    ));
                }
                let mut founder_key_registry = BTreeSet::new();
                let founder_keys = provided
                    .iter()
                    .map(|key| {
                        let key = FounderKey::parse(key.value().clone(), "initial_populations")
                            .map_err(|error| with_initial_population_path(general, error))?;
                        if !founder_key_registry.insert(key.clone()) {
                            return Err(config_error(
                                general,
                                "initial_populations.founder_keys",
                                format!("founder_key repetida: {key}"),
                            ));
                        }
                        Ok(key)
                    })
                    .collect::<Result<Vec<_>, Diagnostic>>()?;
                let individual_ids = founder_keys
                    .iter()
                    .map(|founder_key| {
                        let request = AutomaticIdRequest::new(
                            IdClass::Individual,
                            population_id.as_str(),
                            &[("founder_key", founder_key.as_str())],
                        )
                        .map_err(|error| with_initial_population_path(general, error))?;
                        let id = derivation_registry
                            .derive_individual(root_seed, &request)
                            .map_err(|error| with_initial_population_path(general, error))?;
                        execution_registry
                            .register_individual("execution", id)
                            .map_err(|error| with_initial_population_path(general, error))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                (individual_ids, Some(founder_keys))
            }
            (Some(_), Some(_)) => {
                return Err(config_error(
                    general,
                    "initial_populations.individual_ids,initial_populations.founder_keys",
                    "individual_ids y founder_keys son mutuamente excluyentes",
                ));
            }
            (None, None) => {
                return Err(config_error(
                    general,
                    "initial_populations.founder_keys",
                    "founder_keys es obligatorio cuando se omite individual_ids",
                ));
            }
        };
        effective.push(EffectiveInitialPopulation {
            population_id,
            individual_ids,
            founder_keys,
            ploidy,
            genotype_model,
        });
    }
    Ok(effective)
}

fn with_initial_population_path(
    general: &DeclaredGeneralConfig,
    diagnostic: Diagnostic,
) -> Diagnostic {
    diagnostic.with("path", general.source().path().display().to_string())
}

fn population_snvs_generation_from_effective(
    value: Option<&EffectivePopulationSnvs>,
) -> Option<&PopulationSnvsGenerationConfig> {
    match value {
        Some(EffectivePopulationSnvs::Generate { parameters }) => Some(parameters),
        _ => None,
    }
}

fn effective_population_snvs(
    general: &DeclaredGeneralConfig,
) -> Result<Option<EffectivePopulationSnvs>, Diagnostic> {
    let Some(declared) = general.population_snvs() else {
        return Ok(None);
    };
    match declared {
        DeclaredPopulationSnvs::Supplied { input_path } => {
            return Ok(Some(EffectivePopulationSnvs::Supplied {
                input_path: resolve_link(general.source().path(), input_path.value()),
            }));
        }
        DeclaredPopulationSnvs::Reuse { metadata_path } => {
            return Ok(Some(EffectivePopulationSnvs::Reuse {
                metadata_path: resolve_link(general.source().path(), metadata_path.value()),
            }));
        }
        DeclaredPopulationSnvs::Generate { parameters: _ } => {}
    }
    let DeclaredPopulationSnvs::Generate { parameters } = declared else {
        unreachable!("handled non-generation population_snvs modes")
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
    Ok(Some(EffectivePopulationSnvs::Generate { parameters: PopulationSnvsGenerationConfig {
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
    }}))
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

fn parse_unknown_base_fraction(
    general: &DeclaredGeneralConfig,
    declared: &Option<crate::configuration::DeclaredNumber>,
) -> Result<f64, Diagnostic> {
    let Some(declared) = declared else {
        return Ok(0.0);
    };
    let value = declared.literal().parse::<f64>().map_err(|_| {
        config_error(
            general,
            "reference.composition.N",
            "valor decimal inválido para N",
        )
    })?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(config_error(
            general,
            "reference.composition.N",
            "N debe ser finito y pertenecer a [0,1]",
        ));
    }
    Ok(value)
}

fn resolve_link(general_path: &Path, declared: &str) -> PathBuf {
    // Los campos de ruta que introduce Marimo son texto. Recortamos únicamente
    // espacio accidental en los extremos antes de decidir si es absoluta; los
    // espacios internos se conservan como parte legítima del nombre.
    let declared = Path::new(declared.trim());
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
