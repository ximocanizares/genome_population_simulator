//! Semillas, flujos aleatorios locales y contratos de reproducibilidad.

use std::fmt;
use std::num::NonZeroU64;

use rand_chacha::ChaCha12Rng;
use rand_chacha::rand_core::{Rng, SeedableRng};

use crate::diagnostics::Diagnostic;
use crate::identity::{EventId, LibraryId, RepeatLocusId};
use crate::identity::{F5BlockId, F5SiteId, IndividualId, PopulationId, ScenarioId};

pub const DERIVATION_VERSION: u32 = 1;
pub const F2_DERIVATION_VERSION: u32 = 2;
pub const SUBSEED_DERIVATION_VERSIONS: &[u32] = &[DERIVATION_VERSION, F2_DERIVATION_VERSION];
pub const SUBSEED_CATALOG_VERSION: &str = "9";
pub const DERIVATION_CONTEXT: &str = "genome_population_simulator/v1/subseed";
pub const RANDOM_ALGORITHM: &str = "ChaCha12";
pub const REFERENCE_REDUCTION_SELECTOR: &str = "direct_proportional_selection";

/// Semilla raíz canónica de 256 bits.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct RootSeed([u8; 32]);

impl RootSeed {
    pub fn parse(field: &str, found: &str) -> Result<Self, Diagnostic> {
        if found.len() != 64
            || !found
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return Err(Diagnostic::configuration(
                "GPS-CONF-017",
                "La semilla raíz no contiene exactamente 256 bits en hexadecimal minúsculo.",
            )
            .with("field", field)
            .with("found", found)
            .with("expected", "64 hexadecimal lowercase characters"));
        }
        let mut bytes = [0_u8; 32];
        for (index, pair) in found.as_bytes().chunks_exact(2).enumerate() {
            bytes[index] = decode_nibble(pair[0]) * 16 + decode_nibble(pair[1]);
        }
        Ok(Self(bytes))
    }

    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    #[must_use]
    pub fn to_hex(self) -> String {
        encode_hex(&self.0)
    }
}

impl fmt::Display for RootSeed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&encode_hex(&self.0))
    }
}

fn decode_nibble(byte: u8) -> u8 {
    match byte {
        b'0'..=b'9' => byte - b'0',
        b'a'..=b'f' => byte - b'a' + 10,
        _ => unreachable!("RootSeed::parse validated every hexadecimal digit"),
    }
}

fn encode_hex(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        encoded.push(char::from(HEX[usize::from(byte >> 4)]));
        encoded.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    encoded
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SeedOrigin {
    Declared,
    SystemCsprng,
}

impl SeedOrigin {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Declared => "declared",
            Self::SystemCsprng => "system_csprng",
        }
    }
}

pub trait EntropySource {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String>;
}

#[derive(Debug, Default)]
pub struct SystemEntropy;

impl EntropySource for SystemEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        getrandom::fill(destination).map_err(|error| error.to_string())
    }
}

pub fn resolve_root_seed(
    declared: Option<&str>,
    entropy: &mut impl EntropySource,
) -> Result<(RootSeed, SeedOrigin), Diagnostic> {
    if let Some(declared) = declared {
        return RootSeed::parse("seed", declared).map(|seed| (seed, SeedOrigin::Declared));
    }
    let mut bytes = [0_u8; 32];
    entropy.fill_root_seed(&mut bytes).map_err(|detail| {
        Diagnostic::resource(
            "GPS-RES-001",
            "El CSPRNG del sistema operativo no pudo producir la semilla raíz.",
        )
        .with("operation", "generate_root_seed")
        .with("detail", detail)
    })?;
    Ok((RootSeed::from_bytes(bytes), SeedOrigin::SystemCsprng))
}

/// Ámbitos que tienen un consumidor real durante F1.
///
/// Los catálogos publicados v1–v3 contienen también filas reservadas para fases
/// posteriores. Esas filas se verifican como datos históricos, pero no forman
/// parte de esta frontera de ejecución.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RuntimeSubseedScope<'a> {
    ReferenceGeneration {
        reference_id: &'a str,
    },
    ReferenceReduction {
        source_reference_id: &'a str,
        derived_reference_id: &'a str,
    },
    /// Ámbito del único evento técnico F3, ya publicado como fila histórica
    /// del catálogo v1 y activado ahora por su primer consumidor real.
    TechnicalEvent {
        event_id: &'a str,
        reference_id: &'a str,
    },
    MatrixCandidateCount {
        scenario_id: &'a str,
        reference_id: &'a str,
    },
    MatrixLocusSelection {
        scenario_id: &'a str,
        reference_id: &'a str,
    },
    MatrixAlternateAllele {
        scenario_id: &'a str,
        reference_id: &'a str,
        sequence_id: &'a str,
        position_1based: u64,
        reference_base: &'a str,
    },
    MatrixExactAssignment {
        scenario_id: &'a str,
        reference_id: &'a str,
        population_id: &'a str,
    },
    MatrixFrequencyComponent {
        scenario_id: &'a str,
        reference_id: &'a str,
        population_id: &'a str,
        variant_id: &'a str,
    },
    MatrixFrequencyValue {
        scenario_id: &'a str,
        reference_id: &'a str,
        population_id: &'a str,
        variant_id: &'a str,
    },
    IndelCount {
        scenario_id: &'a str,
        reference_id: &'a str,
        individual_id: &'a str,
        sequence_id: &'a str,
        variant_class: &'a str,
        spatial_scope_id: &'a str,
    },
    IndelPosition {
        scenario_id: &'a str,
        reference_id: &'a str,
        individual_id: &'a str,
        sequence_id: &'a str,
        variant_class: &'a str,
        spatial_scope_id: &'a str,
    },
    InsertionContent {
        scenario_id: &'a str,
        reference_id: &'a str,
        individual_id: &'a str,
        sequence_id: &'a str,
        spatial_scope_id: &'a str,
        event_index: u64,
    },
    InitialGenotype {
        scenario_id: &'a str,
        reference_id: &'a str,
        population_id: &'a str,
        individual_id: &'a str,
        sequence_id: &'a str,
        position_1based: u64,
        reference_base: &'a str,
        alternate_alleles: &'a str,
    },
    /// Fila histórica v4, reservada antes de cerrar la fase multialélica.
    InitialPhase {
        scenario_id: &'a str,
        reference_id: &'a str,
        population_id: &'a str,
        individual_id: &'a str,
        variant_id: &'a str,
    },
    InitialPhaseLocus {
        scenario_id: &'a str,
        reference_id: &'a str,
        population_id: &'a str,
        individual_id: &'a str,
        sequence_id: &'a str,
        position_1based: u64,
        reference_base: &'a str,
        alternate_alleles: &'a str,
    },
    /// Ámbito propio de una plantilla F3. El ordinal es lógico dentro de una
    /// biblioteca, no una posición accidental de iteración o planificación.
    TechnicalTemplate {
        event_id: &'a str,
        library_id: &'a str,
        template_ordinal: u64,
    },
    /// Recuento estocástico WGS por biblioteca. Es distinto del sorteo de
    /// cada plantilla para conservar independencia e invariancia por muestra.
    TechnicalLibraryCount {
        event_id: &'a str,
        library_id: &'a str,
    },
    /// Sorteo del alelo de un locus de repetición en una copia haplotípica.
    ///
    /// Sus componentes son **identidades estables**: individuo, índice de copia
    /// e identificador del locus. Ni la posición del locus en la referencia ni
    /// el orden en que se declaró entran, de modo que insertar un individuo al
    /// principio de la lista, reordenarla o añadir un locus que ordene antes no
    /// cambia ningún genotipo ya sorteado.
    InitialRepeatAllele {
        scenario_id: &'a str,
        reference_id: &'a str,
        population_id: &'a str,
        individual_id: &'a str,
        copy_index: u8,
        locus_id: &'a str,
    },
}

pub const F1_ACTIVE_SUBSEED_SCOPES: [&str; 2] = ["reference_generation", "reference_reduction"];

fn canonical_material(root: RootSeed, components: &[(&str, &str)]) -> Result<Vec<u8>, Diagnostic> {
    let component_count =
        u32::try_from(components.len()).map_err(|_| material_length_error("component_count"))?;
    let mut material = Vec::new();
    material.extend_from_slice(root.as_bytes());
    material.extend_from_slice(&DERIVATION_VERSION.to_be_bytes());
    material.extend_from_slice(&component_count.to_be_bytes());
    for (name, value) in components {
        append_length_prefixed(&mut material, name.as_bytes(), "component_name")?;
        append_length_prefixed(&mut material, value.as_bytes(), "component_value")?;
    }
    Ok(material)
}

fn canonical_material_v2(
    root: RootSeed,
    scope_name: &str,
    components: &[(&str, &str)],
) -> Result<Vec<u8>, Diagnostic> {
    let component_count =
        u32::try_from(components.len()).map_err(|_| material_length_error("component_count"))?;
    let mut material = Vec::new();
    material.extend_from_slice(root.as_bytes());
    material.extend_from_slice(&F2_DERIVATION_VERSION.to_be_bytes());
    append_length_prefixed(&mut material, scope_name.as_bytes(), "scope_name")?;
    material.extend_from_slice(&component_count.to_be_bytes());
    for (name, value) in components {
        append_length_prefixed(&mut material, name.as_bytes(), "component_name")?;
        append_length_prefixed(&mut material, value.as_bytes(), "component_value")?;
    }
    Ok(material)
}

fn append_length_prefixed(
    destination: &mut Vec<u8>,
    value: &[u8],
    field: &str,
) -> Result<(), Diagnostic> {
    let length = u32::try_from(value.len()).map_err(|_| material_length_error(field))?;
    destination.extend_from_slice(&length.to_be_bytes());
    destination.extend_from_slice(value);
    Ok(())
}

fn material_length_error(field: &str) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "subseed_derivation")
    .with("operation", "encode_canonical_material")
    .with("detail", format!("{field} exceeds u32"))
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Subseed([u8; 32]);

impl Subseed {
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    #[must_use]
    pub fn to_hex(self) -> String {
        encode_hex(&self.0)
    }

    #[must_use]
    pub fn local_rng(self) -> LocalRng {
        LocalRng(ChaCha12Rng::from_seed(self.0))
    }
}

pub fn derive_subseed(
    root: RootSeed,
    scope: RuntimeSubseedScope<'_>,
) -> Result<Subseed, Diagnostic> {
    let material = match scope {
        RuntimeSubseedScope::ReferenceGeneration { reference_id } => {
            canonical_material(root, &[("reference_id", reference_id)])?
        }
        RuntimeSubseedScope::ReferenceReduction {
            source_reference_id,
            derived_reference_id,
        } => canonical_material(
            root,
            &[
                ("source_reference_id", source_reference_id),
                ("derived_reference_id", derived_reference_id),
                ("selector", REFERENCE_REDUCTION_SELECTOR),
            ],
        )?,
        RuntimeSubseedScope::TechnicalEvent {
            event_id,
            reference_id,
        } => canonical_material(
            root,
            &[("event_id", event_id), ("reference_id", reference_id)],
        )?,
        RuntimeSubseedScope::MatrixCandidateCount {
            scenario_id,
            reference_id,
        } => canonical_material_v2(
            root,
            "matrix_candidate_count",
            &[("scenario_id", scenario_id), ("reference_id", reference_id)],
        )?,
        RuntimeSubseedScope::MatrixLocusSelection {
            scenario_id,
            reference_id,
        } => canonical_material_v2(
            root,
            "matrix_locus_selection",
            &[("scenario_id", scenario_id), ("reference_id", reference_id)],
        )?,
        RuntimeSubseedScope::MatrixAlternateAllele {
            scenario_id,
            reference_id,
            sequence_id,
            position_1based,
            reference_base,
        } => canonical_material_v2(
            root,
            "matrix_alternate_allele",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("sequence_id", sequence_id),
                ("position_1based", &position_1based.to_string()),
                ("ref", reference_base),
            ],
        )?,
        RuntimeSubseedScope::MatrixExactAssignment {
            scenario_id,
            reference_id,
            population_id,
        } => canonical_material_v2(
            root,
            "matrix_exact_assignment",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("population_id", population_id),
            ],
        )?,
        RuntimeSubseedScope::MatrixFrequencyComponent {
            scenario_id,
            reference_id,
            population_id,
            variant_id,
        } => canonical_material_v2(
            root,
            "matrix_frequency_component",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("population_id", population_id),
                ("variant_id", variant_id),
            ],
        )?,
        RuntimeSubseedScope::MatrixFrequencyValue {
            scenario_id,
            reference_id,
            population_id,
            variant_id,
        } => canonical_material_v2(
            root,
            "matrix_frequency_value",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("population_id", population_id),
                ("variant_id", variant_id),
            ],
        )?,
        RuntimeSubseedScope::IndelCount {
            scenario_id,
            reference_id,
            individual_id,
            sequence_id,
            variant_class,
            spatial_scope_id,
        } => canonical_material_v2(
            root,
            "indel_count",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("individual_id", individual_id),
                ("sequence_id", sequence_id),
                ("variant_class", variant_class),
                ("spatial_scope_id", spatial_scope_id),
            ],
        )?,
        RuntimeSubseedScope::IndelPosition {
            scenario_id,
            reference_id,
            individual_id,
            sequence_id,
            variant_class,
            spatial_scope_id,
        } => canonical_material_v2(
            root,
            "indel_position",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("individual_id", individual_id),
                ("sequence_id", sequence_id),
                ("variant_class", variant_class),
                ("spatial_scope_id", spatial_scope_id),
            ],
        )?,
        RuntimeSubseedScope::InsertionContent {
            scenario_id,
            reference_id,
            individual_id,
            sequence_id,
            spatial_scope_id,
            event_index,
        } => canonical_material_v2(
            root,
            "insertion_content",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("individual_id", individual_id),
                ("sequence_id", sequence_id),
                ("spatial_scope_id", spatial_scope_id),
                ("event_index", &event_index.to_string()),
            ],
        )?,
        RuntimeSubseedScope::InitialGenotype {
            scenario_id,
            reference_id,
            population_id,
            individual_id,
            sequence_id,
            position_1based,
            reference_base,
            alternate_alleles,
        } => canonical_material_v2(
            root,
            "initial_genotype",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("population_id", population_id),
                ("individual_id", individual_id),
                ("sequence_id", sequence_id),
                ("position_1based", &position_1based.to_string()),
                ("ref", reference_base),
                ("alternate_alleles", alternate_alleles),
            ],
        )?,
        RuntimeSubseedScope::InitialPhase {
            scenario_id,
            reference_id,
            population_id,
            individual_id,
            variant_id,
        } => canonical_material_v2(
            root,
            "initial_phase",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("population_id", population_id),
                ("individual_id", individual_id),
                ("variant_id", variant_id),
            ],
        )?,
        RuntimeSubseedScope::InitialPhaseLocus {
            scenario_id,
            reference_id,
            population_id,
            individual_id,
            sequence_id,
            position_1based,
            reference_base,
            alternate_alleles,
        } => canonical_material_v2(
            root,
            "initial_phase_locus",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("population_id", population_id),
                ("individual_id", individual_id),
                ("sequence_id", sequence_id),
                ("position_1based", &position_1based.to_string()),
                ("ref", reference_base),
                ("alternate_alleles", alternate_alleles),
            ],
        )?,
        RuntimeSubseedScope::TechnicalTemplate {
            event_id,
            library_id,
            template_ordinal,
        } => canonical_material_v2(
            root,
            "technical_template",
            &[
                ("event_id", event_id),
                ("library_id", library_id),
                ("template_ordinal", &template_ordinal.to_string()),
            ],
        )?,
        RuntimeSubseedScope::InitialRepeatAllele {
            scenario_id,
            reference_id,
            population_id,
            individual_id,
            copy_index,
            locus_id,
        } => canonical_material_v2(
            root,
            "initial_repeat_allele",
            &[
                ("scenario_id", scenario_id),
                ("reference_id", reference_id),
                ("population_id", population_id),
                ("individual_id", individual_id),
                ("copy_index", &copy_index.to_string()),
                ("locus_id", locus_id),
            ],
        )?,
        RuntimeSubseedScope::TechnicalLibraryCount {
            event_id,
            library_id,
        } => canonical_material_v2(
            root,
            "technical_library_count",
            &[("event_id", event_id), ("library_id", library_id)],
        )?,
    };
    Ok(Subseed(blake3::derive_key(DERIVATION_CONTEXT, &material)))
}

/// Generador local no clonable cuya propiedad queda ligada a un único ámbito.
#[derive(Debug)]
pub struct LocalRng(ChaCha12Rng);

impl LocalRng {
    pub fn next_u64(&mut self) -> u64 {
        self.0.next_u64()
    }

    pub fn fill_bytes(&mut self, destination: &mut [u8]) {
        self.0.fill_bytes(destination);
    }
}

/// Fuente instrumentable de palabras RNG para consumidores normativos.
///
/// Las implementaciones de producción usan [`LocalRng`]. Esta frontera permite
/// verificar vectores de palabras y consumo sin acoplar la especificación a
/// ChaCha ni a una semilla concreta.
pub trait RandomWordSource {
    fn next_word(&mut self) -> u64;
}

impl RandomWordSource for LocalRng {
    fn next_word(&mut self) -> u64 {
        self.next_u64()
    }
}

/// Decisión atómica del catálogo `f5-rng-catalog-v1`.
///
/// Los constructores fijan tanto el orden como la presencia de los componentes
/// del catálogo. `stage` no es una entrada libre: es el literal de decisión.
#[derive(Clone, Copy, Debug)]
pub struct F5RandomDecision<'a>(F5RandomDecisionKind<'a>);

#[derive(Clone, Copy, Debug)]
enum F5RandomDecisionKind<'a> {
    LayoutLength {
        scenario_id: &'a ScenarioId,
        attempt_index: u64,
        block_index: u64,
    },
    LayoutPosition {
        scenario_id: &'a ScenarioId,
        attempt_index: u64,
        block_index: u64,
    },
    LayoutGap {
        scenario_id: &'a ScenarioId,
        attempt_index: u64,
        block_index: u64,
    },
    SiteAlt {
        scenario_id: &'a ScenarioId,
        block_id: &'a F5BlockId,
        site_id: &'a F5SiteId,
    },
    MarginalSpectrum {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        site_id: &'a F5SiteId,
    },
    InitialColumn {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        attempt_index: u64,
        site_index: u64,
    },
    Polarity {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        attempt_index: u64,
        site_index: u64,
    },
    SwapProposal {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        attempt_index: u64,
    },
    FounderSample {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        individual_id: &'a IndividualId,
        copy_index: u32,
    },
    FounderF {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        individual_id: &'a IndividualId,
    },
    ControlledQuotaTie {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
    },
    FounderControlledPermutation {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
    },
    BlockHaplotypeCount {
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
    },
}

impl<'a> F5RandomDecision<'a> {
    #[must_use]
    pub const fn layout_length(
        scenario_id: &'a ScenarioId,
        attempt_index: u64,
        block_index: u64,
    ) -> Self {
        Self(F5RandomDecisionKind::LayoutLength {
            scenario_id,
            attempt_index,
            block_index,
        })
    }
    #[must_use]
    pub const fn layout_position(
        scenario_id: &'a ScenarioId,
        attempt_index: u64,
        block_index: u64,
    ) -> Self {
        Self(F5RandomDecisionKind::LayoutPosition {
            scenario_id,
            attempt_index,
            block_index,
        })
    }
    #[must_use]
    pub const fn layout_gap(
        scenario_id: &'a ScenarioId,
        attempt_index: u64,
        block_index: u64,
    ) -> Self {
        Self(F5RandomDecisionKind::LayoutGap {
            scenario_id,
            attempt_index,
            block_index,
        })
    }
    #[must_use]
    pub const fn site_alt(
        scenario_id: &'a ScenarioId,
        block_id: &'a F5BlockId,
        site_id: &'a F5SiteId,
    ) -> Self {
        Self(F5RandomDecisionKind::SiteAlt {
            scenario_id,
            block_id,
            site_id,
        })
    }
    #[must_use]
    pub const fn marginal_spectrum(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        site_id: &'a F5SiteId,
    ) -> Self {
        Self(F5RandomDecisionKind::MarginalSpectrum {
            scenario_id,
            population_id,
            block_id,
            site_id,
        })
    }
    #[must_use]
    pub const fn initial_column(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        attempt_index: u64,
        site_index: u64,
    ) -> Self {
        Self(F5RandomDecisionKind::InitialColumn {
            scenario_id,
            population_id,
            block_id,
            attempt_index,
            site_index,
        })
    }
    #[must_use]
    pub const fn polarity(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        attempt_index: u64,
        site_index: u64,
    ) -> Self {
        Self(F5RandomDecisionKind::Polarity {
            scenario_id,
            population_id,
            block_id,
            attempt_index,
            site_index,
        })
    }
    #[must_use]
    pub const fn swap_proposal(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        attempt_index: u64,
    ) -> Self {
        Self(F5RandomDecisionKind::SwapProposal {
            scenario_id,
            population_id,
            block_id,
            attempt_index,
        })
    }
    #[must_use]
    pub const fn founder_sample(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        individual_id: &'a IndividualId,
        copy_index: u32,
    ) -> Self {
        Self(F5RandomDecisionKind::FounderSample {
            scenario_id,
            population_id,
            block_id,
            individual_id,
            copy_index,
        })
    }
    #[must_use]
    pub const fn founder_f(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
        individual_id: &'a IndividualId,
    ) -> Self {
        Self(F5RandomDecisionKind::FounderF {
            scenario_id,
            population_id,
            block_id,
            individual_id,
        })
    }
    #[must_use]
    pub const fn controlled_quota_tie(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
    ) -> Self {
        Self(F5RandomDecisionKind::ControlledQuotaTie {
            scenario_id,
            population_id,
            block_id,
        })
    }
    #[must_use]
    pub const fn founder_controlled_permutation(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
    ) -> Self {
        Self(F5RandomDecisionKind::FounderControlledPermutation {
            scenario_id,
            population_id,
            block_id,
        })
    }

    /// Ámbito del recuento haplotípico de un bloque×subpoblación.
    ///
    /// No lleva `attempt_index`: H es una dimensión de la unidad, no del
    /// intento. Si dependiera del intento, dos intentos de la misma unidad
    /// producirían matrices de alturas distintas y ni los fundadores ni los
    /// artefactos podrían indexar un rango haplotípico estable.
    #[must_use]
    pub const fn block_haplotype_count(
        scenario_id: &'a ScenarioId,
        population_id: &'a PopulationId,
        block_id: &'a F5BlockId,
    ) -> Self {
        Self(F5RandomDecisionKind::BlockHaplotypeCount {
            scenario_id,
            population_id,
            block_id,
        })
    }

    fn canonical_components(self) -> Vec<(&'static str, String)> {
        let mut components = Vec::with_capacity(7);
        let (stage, salt) = match self.0 {
            F5RandomDecisionKind::LayoutLength {
                scenario_id,
                attempt_index,
                block_index,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("stage", "layout_length".to_owned()),
                    ("attempt_index_u64", attempt_index.to_string()),
                    ("block_index_u64", block_index.to_string()),
                ]);
                ("layout_length", None)
            }
            F5RandomDecisionKind::LayoutPosition {
                scenario_id,
                attempt_index,
                block_index,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("stage", "layout_position".to_owned()),
                    ("attempt_index_u64", attempt_index.to_string()),
                    ("block_index_u64", block_index.to_string()),
                ]);
                ("layout_position", None)
            }
            F5RandomDecisionKind::LayoutGap {
                scenario_id,
                attempt_index,
                block_index,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("stage", "layout_gap".to_owned()),
                    ("attempt_index_u64", attempt_index.to_string()),
                    ("block_index_u64", block_index.to_string()),
                ]);
                ("layout_gap", None)
            }
            F5RandomDecisionKind::SiteAlt {
                scenario_id,
                block_id,
                site_id,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("site_id", site_id.as_str().to_owned()),
                    ("stage", "site_alt".to_owned()),
                ]);
                ("site_alt", None)
            }
            F5RandomDecisionKind::MarginalSpectrum {
                scenario_id,
                population_id,
                block_id,
                site_id,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("site_id", site_id.as_str().to_owned()),
                    ("stage", "marginal_spectrum".to_owned()),
                ]);
                ("marginal_spectrum", None)
            }
            F5RandomDecisionKind::InitialColumn {
                scenario_id,
                population_id,
                block_id,
                attempt_index,
                site_index,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("stage", "initial_column".to_owned()),
                    ("attempt_index_u64", attempt_index.to_string()),
                ]);
                ("initial_column", Some(site_index.to_string()))
            }
            F5RandomDecisionKind::Polarity {
                scenario_id,
                population_id,
                block_id,
                attempt_index,
                site_index,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("stage", "polarity".to_owned()),
                    ("attempt_index_u64", attempt_index.to_string()),
                ]);
                ("polarity", Some(site_index.to_string()))
            }
            F5RandomDecisionKind::SwapProposal {
                scenario_id,
                population_id,
                block_id,
                attempt_index,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("stage", "swap_proposal".to_owned()),
                    ("attempt_index_u64", attempt_index.to_string()),
                ]);
                ("swap_proposal", None)
            }
            F5RandomDecisionKind::FounderSample {
                scenario_id,
                population_id,
                block_id,
                individual_id,
                copy_index,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("individual_id", individual_id.as_str().to_owned()),
                    ("copy_index_u32", copy_index.to_string()),
                    ("stage", "founder_sample".to_owned()),
                ]);
                ("founder_sample", None)
            }
            F5RandomDecisionKind::FounderF {
                scenario_id,
                population_id,
                block_id,
                individual_id,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("individual_id", individual_id.as_str().to_owned()),
                    ("stage", "founder_f".to_owned()),
                ]);
                ("founder_f", None)
            }
            F5RandomDecisionKind::ControlledQuotaTie {
                scenario_id,
                population_id,
                block_id,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("stage", "controlled_quota_tie".to_owned()),
                ]);
                ("controlled_quota_tie", None)
            }
            F5RandomDecisionKind::FounderControlledPermutation {
                scenario_id,
                population_id,
                block_id,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("stage", "founder_controlled_permutation".to_owned()),
                ]);
                ("founder_controlled_permutation", None)
            }
            F5RandomDecisionKind::BlockHaplotypeCount {
                scenario_id,
                population_id,
                block_id,
            } => {
                components.extend([
                    ("scenario_id", scenario_id.as_str().to_owned()),
                    ("population_id", population_id.as_str().to_owned()),
                    ("block_id", block_id.as_str().to_owned()),
                    ("stage", "block_haplotype_count".to_owned()),
                ]);
                ("block_haplotype_count", None)
            }
        };
        debug_assert_eq!(
            components
                .iter()
                .find(|(name, _)| *name == "stage")
                .map(|(_, value)| value.as_str()),
            Some(stage)
        );
        if let Some(site_index) = salt {
            components.push(("salt.site_index_u64", site_index));
        }
        components
    }
}

/// Deriva un ChaCha12 local de una decisión F5, sin flujo global compartido.
pub fn derive_f5_subseed(
    root: RootSeed,
    decision: F5RandomDecision<'_>,
) -> Result<Subseed, Diagnostic> {
    let components = decision.canonical_components();
    let borrowed = components
        .iter()
        .map(|(name, value)| (*name, value.as_str()))
        .collect::<Vec<_>>();
    Ok(Subseed(blake3::derive_key(
        DERIVATION_CONTEXT,
        &canonical_material(root, &borrowed)?,
    )))
}

/// Codificación de desempate/layout de DTC-089: cuatro `u64` big-endian.
#[must_use]
pub fn f5_placed_layout_bytes(rows: &[(u64, u64, u64, u64)]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(rows.len().saturating_mul(32));
    for (sequence_ordinal, start_1based, end_1based, length_bp) in rows {
        for value in [sequence_ordinal, start_1based, end_1based, length_bp] {
            bytes.extend_from_slice(&value.to_be_bytes());
        }
    }
    bytes
}

/// Codificación de matriz binaria F5: H/S big-endian, filas por rango y bits
/// MSB primero con relleno bajo a cero.
pub fn f5_matrix_bytes(rows: &[Vec<bool>]) -> Result<Vec<u8>, RandomSamplingError> {
    let haplotype_count =
        u64::try_from(rows.len()).map_err(|_| RandomSamplingError::NumericFailure)?;
    let site_count = rows.first().map_or(0, Vec::len);
    if rows.iter().any(|row| row.len() != site_count) {
        return Err(RandomSamplingError::InvalidMatrix);
    }
    let site_count_u64 =
        u64::try_from(site_count).map_err(|_| RandomSamplingError::NumericFailure)?;
    let bytes_per_row = site_count.div_ceil(8);
    let payload = rows
        .len()
        .checked_mul(bytes_per_row)
        .ok_or(RandomSamplingError::NumericFailure)?;
    let mut bytes = Vec::with_capacity(
        16usize
            .checked_add(payload)
            .ok_or(RandomSamplingError::NumericFailure)?,
    );
    bytes.extend_from_slice(&haplotype_count.to_be_bytes());
    bytes.extend_from_slice(&site_count_u64.to_be_bytes());
    for row in rows {
        for chunk in row.chunks(8) {
            let mut byte = 0_u8;
            for (offset, allele) in chunk.iter().enumerate() {
                if *allele {
                    byte |= 1_u8 << (7 - offset);
                }
            }
            bytes.push(byte);
        }
    }
    Ok(bytes)
}

pub fn f5_uniform_index(rng: &mut LocalRng, n: u64) -> Result<u64, RandomSamplingError> {
    f5_uniform_index_with_source(rng, n)
}

pub(crate) fn f5_uniform_index_with_source(
    rng: &mut impl RandomWordSource,
    n: u64,
) -> Result<u64, RandomSamplingError> {
    if n == 0 {
        return Err(RandomSamplingError::InvalidBound);
    }
    if n == 1 {
        return Ok(0);
    }
    let remainder = (u64::MAX % n + 1) % n;
    let limit = (remainder != 0).then(|| u64::MAX - (remainder - 1));
    loop {
        let value = rng.next_word();
        if limit.is_none_or(|bound| value < bound) {
            return Ok(value % n);
        }
    }
}

/// Los 53 bits altos de una palabra F5: `bits53 = next_u64() >> 11`.
///
/// Se expone aparte del uniforme en binary64 porque Bernoulli compara enteros
/// y no debe perder ningún bit convirtiendo antes a `f64`.
pub fn f5_bits_53(rng: &mut impl RandomWordSource) -> u64 {
    rng.next_word() >> 11
}

/// U53 normativo de DTC-089: `U = bits53 / 2^53`, cerrado en 0 y abierto en 1.
///
/// A diferencia de [`f4_uniform_open_53`] no rechaza ninguna palabra ni
/// desplaza el valor media unidad: consume exactamente una palabra, devuelve 0
/// exacto y nunca devuelve 1. Es el único uniforme continuo de los ámbitos F5.
pub fn f5_uniform_53(rng: &mut impl RandomWordSource) -> f64 {
    f5_bits_53(rng) as f64 / F5_FIFTY_THREE_BIT_BOUND as f64
}

/// Bernoulli DTC-089. Los extremos no abren ni avanzan el RNG.
pub fn f5_bernoulli(rng: &mut LocalRng, probability: f64) -> Result<bool, RandomSamplingError> {
    f5_bernoulli_with_source(rng, probability)
}

pub(crate) fn f5_bernoulli_with_source(
    rng: &mut impl RandomWordSource,
    probability: f64,
) -> Result<bool, RandomSamplingError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(RandomSamplingError::InvalidProbability);
    }
    if probability == 0.0 {
        return Ok(false);
    }
    if probability == 1.0 {
        return Ok(true);
    }
    let threshold = (probability * F5_FIFTY_THREE_BIT_BOUND as f64).floor() as u64;
    Ok(f5_bits_53(rng) < threshold)
}

pub fn f5_fisher_yates<T>(rng: &mut LocalRng, values: &mut [T]) -> Result<(), RandomSamplingError> {
    f5_fisher_yates_with_source(rng, values)
}

pub(crate) fn f5_fisher_yates_with_source<T>(
    rng: &mut impl RandomWordSource,
    values: &mut [T],
) -> Result<(), RandomSamplingError> {
    for index in (1..values.len()).rev() {
        let other = usize::try_from(f5_uniform_index_with_source(rng, (index + 1) as u64)?)
            .map_err(|_| RandomSamplingError::NumericFailure)?;
        values.swap(index, other);
    }
    Ok(())
}

pub const POISSON_COUNT_ALGORITHM: &str = "poisson-count-v1";
/// Mayor `lambda` que admite `poisson-count-v1`: `2^52`.
///
/// DTC-067 la fija y exige que un valor mayor **se rechace antes de consumir
/// RNG**, con valor, límite y ámbito. La constante vive aquí para que el
/// muestreador y quien valide una lambda declarada consuman la misma cifra: si
/// solo el muestreador la conociera, un validador podría admitir una lambda
/// que después muere como error interno a mitad de ejecución.
pub const POISSON_COUNT_MAX_LAMBDA: f64 = 4_503_599_627_370_496.0;
pub const HWE_COMPOSITION_ALGORITHM: &str = "hwe-composition-v1";
pub const F_DIPLOID_GENOTYPE_ALGORITHM: &str = "f-diploid-genotype-v1";
pub const PHASE_PERMUTATION_ALGORITHM: &str = "phase-permutation-v1";
/// Elección ponderada de una categoría entre varias (DTC-072).
pub const CATEGORICAL_F64_ALGORITHM: &str = "categorical-f64-v1";
pub const F4_UNIFORM_INDEX_ALGORITHM: &str = "uniform-index-v1";
pub const F4_UNIFORM_OPEN_53_ALGORITHM: &str = "uniform-open-53-v1";
pub const F4_UNIFORM_CLOSED_53_ALGORITHM: &str = "uniform-closed-53-bit-v1";
pub const F4_BERNOULLI_ALGORITHM: &str = "bernoulli-v1";
pub const F4_PARTIAL_FISHER_YATES_ALGORITHM: &str = "partial-fisher-yates-v1";
pub const F4_TRUNCATED_NORMAL_ALGORITHM: &str = "normal-truncated-box-muller-v1";
/// Catálogo normativo de claves técnicas F4 de DTC-086.
///
/// **Es además el ámbito de derivación** de cada una de sus trece decisiones:
/// el nombre entra literalmente en el material que `canonical_material_v2`
/// amasa. Por eso no puede reescribirse al subir de versión sin reescribir
/// todas las subsemillas ya publicadas.
pub const F4_TECHNICAL_RNG_CATALOG_V3: &str = "f4-technical-v3";
/// Catálogo que **extiende** el anterior con la etapa `stutter` (DTC-090).
///
/// # Por qué la extensión no desplaza nada
///
/// Un catálogo nuevo podría entenderse de dos maneras: como un ámbito de
/// derivación nuevo para **todas** las decisiones, o como el vigente **más**
/// una etapa. Aquí es lo segundo, y no por comodidad: DTC-081 exige que cada
/// sorteo arranque de su propia clave y que «no existe consumo acumulado entre
/// decisiones», y la decisión S-11 saca de ahí la consecuencia explícita de que
/// **añadir el tartamudeo no puede desplazar ninguna decisión anterior**.
///
/// En consecuencia, las trece decisiones heredadas siguen derivándose bajo
/// [`F4_TECHNICAL_RNG_CATALOG_V3`] y solo las dos decisiones de la etapa
/// `stutter` usan este nombre como ámbito. Lo que cambia al subir de versión es
/// **qué claves existen**, no qué produce ninguna de las que ya existían.
///
/// Una ejecución publica el catálogo del que salen **sus** decisiones, como ya
/// hace F5 con `f5-rng-catalog-v1`/`v2`: un evento técnico sin tartamudeo
/// activo no abre la etapa `stutter` y publica el v3, que es lo que mantiene
/// byte a byte los paquetes anteriores a esta ampliación.
pub const F4_TECHNICAL_RNG_CATALOG_V4: &str = "f4-technical-v4";
/// Catálogo que **extiende** los anteriores con la etapa `amplification`
/// (F4.9).
///
/// # Por qué una versión nueva y no más claves dentro del v4
///
/// Por la misma razón que el v4 no se metió dentro del v3: una ejecución
/// **publica el catálogo del que salen sus decisiones**, y ese nombre es
/// además el ámbito de derivación. Si la amplificación naciera en el v4, una
/// biblioteca amplificada sin tartamudeo publicaría un catálogo que el
/// contrato ata al tartamudeo, y un lector no podría distinguir qué etapa abrió
/// de verdad. El propio diseño con PCR lo dejó escrito: «`f4-technical-v4` si
/// se hace junto con `stutter`, **la siguiente si después**».
///
/// Como en la ampliación anterior, **ninguna decisión ya publicada se mueve**:
/// las trece heredadas siguen en el v3, las dos del tartamudeo en el v4, y solo
/// las claves que nacen con la amplificación usan este nombre como ámbito.
/// Entre ellas están las tres claves de `fragment` **con componente `attempt`**,
/// que solo existen cuando la aceptación por composición rechaza un intento: el
/// intento `0` conserva literalmente las claves de hoy, y eso es lo que permite
/// que una biblioteca amplificada con parámetros neutros produzca exactamente
/// los mismos bytes que una sin amplificar.
pub const F4_TECHNICAL_RNG_CATALOG_V5: &str = "f4-technical-v5";
/// Catálogo vigente para una ejecución que no abre la etapa `stutter`.
pub const F4_TECHNICAL_RNG_CATALOG: &str = F4_TECHNICAL_RNG_CATALOG_V3;
/// Catálogo inmutable de ámbitos aleatorios de F5 (DTC-089).
///
/// Este catálogo no amplía `subseed-scopes-v8`: F5 publica sus decisiones en
/// un artefacto independiente para no reescribir los vectores heredados.
///
/// v2 añade `block_haplotype_count`, el ámbito del recuento haplotípico por
/// bloque×subpoblación. Las doce filas de v1 se conservan **literales**: el
/// material, la subsemilla y las cuatro primeras palabras de cada decisión
/// heredada son las mismas, de modo que subir de catálogo no reescribe
/// ninguna decisión ya congelada.
///
/// Una ejecución publica el catálogo del que salen **sus** decisiones, no el
/// más alto que el binario conoce: un escenario cuyo perfil no declara la
/// distribución del recuento haplotípico no puede abrir
/// `block_haplotype_count` y publica [`F5_RNG_CATALOG_V1`]. Eso es lo que
/// mantiene byte a byte los paquetes anteriores a esta ampliación; publicar
/// v2 en todos habría cambiado el `manifest.json` de cada uno de ellos sin
/// que ninguna decisión cambiase.
pub const F5_RNG_CATALOG_V1: &str = "f5-rng-catalog-v1";
/// Catálogo de una ejecución que sí declara la distribución del recuento.
pub const F5_RNG_CATALOG_V2: &str = "f5-rng-catalog-v2";
/// Catálogo más alto que este binario implementa.
pub const F5_RNG_CATALOG: &str = F5_RNG_CATALOG_V2;
pub const F5_UNIFORM_INDEX_ALGORITHM: &str = "uniform-index-v1";
pub const F5_BERNOULLI_ALGORITHM: &str = "bernoulli-v1";
pub const F5_FISHER_YATES_ALGORITHM: &str = "fisher-yates-v1";

const F5_FIFTY_THREE_BIT_BOUND: u64 = 1_u64 << 53;
const F4_FIFTY_THREE_BIT_BOUND: u64 = 1_u64 << 53;
const F4_FIFTY_THREE_BIT_MAX: u64 = F4_FIFTY_THREE_BIT_BOUND - 1;
const F4_PI: f64 = f64::from_bits(0x4009_21fb_5444_2d18);

/// Mate de una pareja PE en las decisiones de error F4.
///
/// No se admite un literal arbitrario para impedir que una clave de error
/// represente un mate ajeno al contrato PE de F4.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F4Mate {
    R1,
    R2,
}

impl F4Mate {
    const fn canonical_value(self) -> &'static str {
        match self {
            Self::R1 => "R1",
            Self::R2 => "R2",
        }
    }
}

/// Identificador canónico de una instancia previa al index hopping de F4.
///
/// La validación permanece junto a la frontera del catálogo RNG: la identidad
/// de dominio que lo produzca se introducirá en una fase posterior. Mientras
/// tanto, solo se admite la categoría `PH_` y la sintaxis ASCII canónica de
/// identificadores aportados.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F4IndexHopInstance<'a>(&'a str);

impl<'a> F4IndexHopInstance<'a> {
    pub fn parse(value: &'a str) -> Result<Self, Diagnostic> {
        let suffix = value.strip_prefix("PH_");
        let is_valid = suffix.is_some_and(|suffix| {
            !suffix.is_empty()
                && value.len() <= 64
                && suffix
                    .as_bytes()
                    .first()
                    .is_some_and(u8::is_ascii_alphanumeric)
                && suffix
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
        });
        if is_valid {
            return Ok(Self(value));
        }
        Err(Diagnostic::domain(
            "GPS-DOMAIN-001",
            "Un identificador aportado incumple la sintaxis o el ámbito de su categoría.",
        )
        .with("identifier", value)
        .with("category", "f4_index_hop_instance")
        .with("scope", "f4_technical_rng"))
    }

    #[must_use]
    pub const fn as_str(self) -> &'a str {
        self.0
    }
}

/// Decisión atómica del catálogo `f4-technical-v3`.
///
/// Cada variante fija el ámbito, el discriminador y los componentes que DTC-086
/// exige. Los campos permanecen privados para que no pueda construirse una
/// clave incompleta, ni combinar una decisión con componentes de otra etapa.
#[derive(Clone, Copy, Debug)]
pub struct F4TechnicalDecision<'a>(F4TechnicalDecisionKind<'a>);

#[derive(Clone, Copy, Debug)]
enum F4TechnicalDecisionKind<'a> {
    FragmentInsertLength {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    },
    FragmentCopy {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    },
    FragmentStart {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    },
    ChimeraBernoulli {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    },
    ChimeraBreakpoint {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    },
    ChimeraSecondSource {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        breakpoint: NonZeroU64,
    },
    ErrorBernoulli {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        mate: F4Mate,
        cycle: u64,
    },
    ErrorSubstitution {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        mate: F4Mate,
        cycle: u64,
    },
    ContaminationPercent {
        event_id: &'a EventId,
        recipient_library_id: &'a LibraryId,
        recipient_id: &'a IndividualId,
    },
    ContaminationDonors {
        event_id: &'a EventId,
        recipient_library_id: &'a LibraryId,
        recipient_id: &'a IndividualId,
    },
    ContaminationRecipients {
        event_id: &'a EventId,
        recipient_library_id: &'a LibraryId,
        recipient_id: &'a IndividualId,
    },
    IndexHopBernoulli {
        event_id: &'a EventId,
        pre_hop_library_id: &'a LibraryId,
        pre_hop_instance_id: F4IndexHopInstance<'a>,
    },
    IndexHopDestination {
        event_id: &'a EventId,
        pre_hop_library_id: &'a LibraryId,
        pre_hop_instance_id: F4IndexHopInstance<'a>,
    },
    StutterBernoulli {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        locus_id: &'a RepeatLocusId,
    },
    StutterStep {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        locus_id: &'a RepeatLocusId,
    },
    /// Sorteo de `fragment` de un intento de aceptación `a >= 1`. El intento
    /// `0` usa las tres variantes heredadas, sin componente `attempt`.
    FragmentRetry {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        decision: FragmentRetryDecision,
        attempt: NonZeroU64,
    },
    AmplificationEfficiencyClass {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        family_ordinal: u64,
    },
    AmplificationEarlyCopies {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        family_ordinal: u64,
    },
    AmplificationGcAccept {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        family_ordinal: u64,
        attempt: u64,
    },
    AmplificationFamilySize {
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        family_ordinal: u64,
    },
}

/// Cuál de los tres sorteos de `fragment` se repite en un reintento.
///
/// Es una enumeración y no una cadena para que no pueda inventarse un nombre
/// de decisión que el catálogo no contiene.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FragmentRetryDecision {
    InsertLength,
    Copy,
    Start,
}

impl FragmentRetryDecision {
    const fn canonical_value(self) -> &'static str {
        match self {
            Self::InsertLength => "insert_length",
            Self::Copy => "copy",
            Self::Start => "start",
        }
    }
}

impl<'a> F4TechnicalDecision<'a> {
    #[must_use]
    pub const fn fragment_insert_length(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::FragmentInsertLength {
            event_id,
            library_id,
            template_ordinal,
        })
    }

    #[must_use]
    pub const fn fragment_copy(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::FragmentCopy {
            event_id,
            library_id,
            template_ordinal,
        })
    }

    #[must_use]
    pub const fn fragment_start(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::FragmentStart {
            event_id,
            library_id,
            template_ordinal,
        })
    }

    #[must_use]
    pub const fn chimera_bernoulli(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::ChimeraBernoulli {
            event_id,
            library_id,
            template_ordinal,
        })
    }

    #[must_use]
    pub const fn chimera_breakpoint(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::ChimeraBreakpoint {
            event_id,
            library_id,
            template_ordinal,
        })
    }

    #[must_use]
    pub const fn chimera_second_source(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        breakpoint: NonZeroU64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::ChimeraSecondSource {
            event_id,
            library_id,
            template_ordinal,
            breakpoint,
        })
    }

    #[must_use]
    pub const fn error_bernoulli(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        mate: F4Mate,
        cycle: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::ErrorBernoulli {
            event_id,
            library_id,
            template_ordinal,
            mate,
            cycle,
        })
    }

    #[must_use]
    pub const fn error_substitution(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        mate: F4Mate,
        cycle: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::ErrorSubstitution {
            event_id,
            library_id,
            template_ordinal,
            mate,
            cycle,
        })
    }

    #[must_use]
    pub const fn contamination_percent(
        event_id: &'a EventId,
        recipient_library_id: &'a LibraryId,
        recipient_id: &'a IndividualId,
    ) -> Self {
        Self(F4TechnicalDecisionKind::ContaminationPercent {
            event_id,
            recipient_library_id,
            recipient_id,
        })
    }

    #[must_use]
    pub const fn contamination_donors(
        event_id: &'a EventId,
        recipient_library_id: &'a LibraryId,
        recipient_id: &'a IndividualId,
    ) -> Self {
        Self(F4TechnicalDecisionKind::ContaminationDonors {
            event_id,
            recipient_library_id,
            recipient_id,
        })
    }

    #[must_use]
    pub const fn contamination_recipients(
        event_id: &'a EventId,
        recipient_library_id: &'a LibraryId,
        recipient_id: &'a IndividualId,
    ) -> Self {
        Self(F4TechnicalDecisionKind::ContaminationRecipients {
            event_id,
            recipient_library_id,
            recipient_id,
        })
    }

    #[must_use]
    pub const fn index_hop_bernoulli(
        event_id: &'a EventId,
        pre_hop_library_id: &'a LibraryId,
        pre_hop_instance_id: F4IndexHopInstance<'a>,
    ) -> Self {
        Self(F4TechnicalDecisionKind::IndexHopBernoulli {
            event_id,
            pre_hop_library_id,
            pre_hop_instance_id,
        })
    }

    #[must_use]
    pub const fn index_hop_destination(
        event_id: &'a EventId,
        pre_hop_library_id: &'a LibraryId,
        pre_hop_instance_id: F4IndexHopInstance<'a>,
    ) -> Self {
        Self(F4TechnicalDecisionKind::IndexHopDestination {
            event_id,
            pre_hop_library_id,
            pre_hop_instance_id,
        })
    }

    /// Ensayo binario «¿resbala la polimerasa sobre este tracto?».
    ///
    /// Se indexa por **identificador de locus**, nunca por la posición del
    /// tracto dentro de la molécula: así, declarar un locus nuevo en otro
    /// cromosoma no desplaza las decisiones de los que ya había (S-11).
    #[must_use]
    pub const fn stutter_bernoulli(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        locus_id: &'a RepeatLocusId,
    ) -> Self {
        Self(F4TechnicalDecisionKind::StutterBernoulli {
            event_id,
            library_id,
            template_ordinal,
            locus_id,
        })
    }

    /// Elección del paso en unidades enteras entre los candidatos **ya
    /// filtrados**. No es un recuento, de modo que DTC-067 no le aplica.
    #[must_use]
    pub const fn stutter_step(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        locus_id: &'a RepeatLocusId,
    ) -> Self {
        Self(F4TechnicalDecisionKind::StutterStep {
            event_id,
            library_id,
            template_ordinal,
            locus_id,
        })
    }

    /// Sorteo de `fragment` repetido en el intento `attempt >= 1` de la
    /// aceptación por composición.
    #[must_use]
    pub const fn fragment_retry(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        template_ordinal: u64,
        decision: FragmentRetryDecision,
        attempt: NonZeroU64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::FragmentRetry {
            event_id,
            library_id,
            template_ordinal,
            decision,
            attempt,
        })
    }

    /// Clase de eficiencia de la molécula que da origen a una familia.
    #[must_use]
    pub const fn amplification_efficiency_class(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        family_ordinal: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::AmplificationEfficiencyClass {
            event_id,
            library_id,
            family_ordinal,
        })
    }

    /// Copias tras los primeros ciclos, sorteadas de la distribución sesgada
    /// por tamaño.
    #[must_use]
    pub const fn amplification_early_copies(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        family_ordinal: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::AmplificationEarlyCopies {
            event_id,
            library_id,
            family_ordinal,
        })
    }

    /// Aceptación por composición de un intento concreto.
    #[must_use]
    pub const fn amplification_gc_accept(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        family_ordinal: u64,
        attempt: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::AmplificationGcAccept {
            event_id,
            library_id,
            family_ordinal,
            attempt,
        })
    }

    /// Tamaño de la familia.
    #[must_use]
    pub const fn amplification_family_size(
        event_id: &'a EventId,
        library_id: &'a LibraryId,
        family_ordinal: u64,
    ) -> Self {
        Self(F4TechnicalDecisionKind::AmplificationFamilySize {
            event_id,
            library_id,
            family_ordinal,
        })
    }

    /// Catálogo del que sale **esta** decisión, que es su ámbito de derivación.
    ///
    /// Las trece decisiones de DTC-086 conservan el v3 para que una etapa nueva
    /// no reescriba ninguna subsemilla ya publicada; solo las dos del
    /// tartamudeo nacen en el v4.
    /// **Las trece heredadas se enumeran una a una, sin comodín.** Con un
    /// comodín al final, una decisión **nueva de la familia del tartamudeo**
    /// caería en el v3 en silencio, que es exactamente el error que toda esta
    /// pieza existe para evitar. Enumerándolas, añadir una variante no compila
    /// hasta que alguien decide a qué catálogo pertenece.
    #[must_use]
    pub const fn catalog(self) -> &'static str {
        match self.0 {
            F4TechnicalDecisionKind::StutterBernoulli { .. }
            | F4TechnicalDecisionKind::StutterStep { .. } => F4_TECHNICAL_RNG_CATALOG_V4,
            F4TechnicalDecisionKind::FragmentRetry { .. }
            | F4TechnicalDecisionKind::AmplificationEfficiencyClass { .. }
            | F4TechnicalDecisionKind::AmplificationEarlyCopies { .. }
            | F4TechnicalDecisionKind::AmplificationGcAccept { .. }
            | F4TechnicalDecisionKind::AmplificationFamilySize { .. } => {
                F4_TECHNICAL_RNG_CATALOG_V5
            }
            F4TechnicalDecisionKind::FragmentInsertLength { .. }
            | F4TechnicalDecisionKind::FragmentCopy { .. }
            | F4TechnicalDecisionKind::FragmentStart { .. }
            | F4TechnicalDecisionKind::ChimeraBernoulli { .. }
            | F4TechnicalDecisionKind::ChimeraBreakpoint { .. }
            | F4TechnicalDecisionKind::ChimeraSecondSource { .. }
            | F4TechnicalDecisionKind::ErrorBernoulli { .. }
            | F4TechnicalDecisionKind::ErrorSubstitution { .. }
            | F4TechnicalDecisionKind::ContaminationPercent { .. }
            | F4TechnicalDecisionKind::ContaminationDonors { .. }
            | F4TechnicalDecisionKind::ContaminationRecipients { .. }
            | F4TechnicalDecisionKind::IndexHopBernoulli { .. }
            | F4TechnicalDecisionKind::IndexHopDestination { .. } => F4_TECHNICAL_RNG_CATALOG_V3,
        }
    }

    fn canonical_components(self) -> Vec<(&'static str, String)> {
        match self.0 {
            F4TechnicalDecisionKind::FragmentInsertLength {
                event_id,
                library_id,
                template_ordinal,
            } => f4_template_components(
                event_id,
                library_id,
                "fragment",
                template_ordinal,
                "insert_length",
                None,
            ),
            F4TechnicalDecisionKind::FragmentCopy {
                event_id,
                library_id,
                template_ordinal,
            } => f4_template_components(
                event_id,
                library_id,
                "fragment",
                template_ordinal,
                "copy",
                None,
            ),
            F4TechnicalDecisionKind::FragmentStart {
                event_id,
                library_id,
                template_ordinal,
            } => f4_template_components(
                event_id,
                library_id,
                "fragment",
                template_ordinal,
                "start",
                None,
            ),
            F4TechnicalDecisionKind::ChimeraBernoulli {
                event_id,
                library_id,
                template_ordinal,
            } => f4_template_components(
                event_id,
                library_id,
                "chimera",
                template_ordinal,
                "bernoulli",
                None,
            ),
            F4TechnicalDecisionKind::ChimeraBreakpoint {
                event_id,
                library_id,
                template_ordinal,
            } => f4_template_components(
                event_id,
                library_id,
                "chimera",
                template_ordinal,
                "breakpoint",
                None,
            ),
            F4TechnicalDecisionKind::ChimeraSecondSource {
                event_id,
                library_id,
                template_ordinal,
                breakpoint,
            } => f4_template_components(
                event_id,
                library_id,
                "chimera",
                template_ordinal,
                "second_source",
                Some(breakpoint.get()),
            ),
            F4TechnicalDecisionKind::ErrorBernoulli {
                event_id,
                library_id,
                template_ordinal,
                mate,
                cycle,
            } => f4_error_components(
                event_id,
                library_id,
                template_ordinal,
                mate,
                cycle,
                "bernoulli",
            ),
            F4TechnicalDecisionKind::ErrorSubstitution {
                event_id,
                library_id,
                template_ordinal,
                mate,
                cycle,
            } => f4_error_components(
                event_id,
                library_id,
                template_ordinal,
                mate,
                cycle,
                "substitution",
            ),
            F4TechnicalDecisionKind::ContaminationPercent {
                event_id,
                recipient_library_id,
                recipient_id,
            } => f4_contamination_components(
                event_id,
                recipient_library_id,
                "contamination_percent",
                recipient_id,
                "percent",
            ),
            F4TechnicalDecisionKind::ContaminationDonors {
                event_id,
                recipient_library_id,
                recipient_id,
            } => f4_contamination_components(
                event_id,
                recipient_library_id,
                "contamination_donors",
                recipient_id,
                "select_donors",
            ),
            F4TechnicalDecisionKind::ContaminationRecipients {
                event_id,
                recipient_library_id,
                recipient_id,
            } => f4_contamination_components(
                event_id,
                recipient_library_id,
                "contamination_remove",
                recipient_id,
                "select_recipients",
            ),
            F4TechnicalDecisionKind::IndexHopBernoulli {
                event_id,
                pre_hop_library_id,
                pre_hop_instance_id,
            } => f4_index_hop_components(
                event_id,
                pre_hop_library_id,
                pre_hop_instance_id,
                "bernoulli",
            ),
            F4TechnicalDecisionKind::IndexHopDestination {
                event_id,
                pre_hop_library_id,
                pre_hop_instance_id,
            } => f4_index_hop_components(
                event_id,
                pre_hop_library_id,
                pre_hop_instance_id,
                "destination",
            ),
            F4TechnicalDecisionKind::StutterBernoulli {
                event_id,
                library_id,
                template_ordinal,
                locus_id,
            } => f4_stutter_components(
                event_id,
                library_id,
                template_ordinal,
                locus_id,
                "bernoulli",
            ),
            F4TechnicalDecisionKind::StutterStep {
                event_id,
                library_id,
                template_ordinal,
                locus_id,
            } => f4_stutter_components(event_id, library_id, template_ordinal, locus_id, "step"),
            F4TechnicalDecisionKind::FragmentRetry {
                event_id,
                library_id,
                template_ordinal,
                decision,
                attempt,
            } => {
                let mut components = f4_template_components(
                    event_id,
                    library_id,
                    "fragment",
                    template_ordinal,
                    decision.canonical_value(),
                    None,
                );
                components.push(("attempt", attempt.get().to_string()));
                components
            }
            F4TechnicalDecisionKind::AmplificationEfficiencyClass {
                event_id,
                library_id,
                family_ordinal,
            } => f4_amplification_components(
                event_id,
                library_id,
                family_ordinal,
                "efficiency_class",
                None,
            ),
            F4TechnicalDecisionKind::AmplificationEarlyCopies {
                event_id,
                library_id,
                family_ordinal,
            } => f4_amplification_components(
                event_id,
                library_id,
                family_ordinal,
                "early_copies",
                None,
            ),
            F4TechnicalDecisionKind::AmplificationGcAccept {
                event_id,
                library_id,
                family_ordinal,
                attempt,
            } => f4_amplification_components(
                event_id,
                library_id,
                family_ordinal,
                "gc_accept",
                Some(attempt),
            ),
            F4TechnicalDecisionKind::AmplificationFamilySize {
                event_id,
                library_id,
                family_ordinal,
            } => f4_amplification_components(
                event_id,
                library_id,
                family_ordinal,
                "family_size",
                None,
            ),
        }
    }
}

/// Componentes de la etapa `amplification`, en el orden literal de DTC-086:
/// `[event_id, library_id, amplification, family_ordinal, decision]`, con el
/// `attempt` detrás cuando la decisión lo distingue.
///
/// El **ordinal de familia** es el ordinal de plantilla de su primera pareja.
/// Eso hace que la clave de una familia no dependa de cuántas parejas emitió
/// ninguna familia anterior más allá de lo que ya fija el recorrido, y que
/// añadir una biblioteca no mueva nada de las demás.
fn f4_amplification_components(
    event_id: &EventId,
    library_id: &LibraryId,
    family_ordinal: u64,
    decision: &'static str,
    attempt: Option<u64>,
) -> Vec<(&'static str, String)> {
    let mut components = vec![
        ("event_id", event_id.as_str().to_owned()),
        ("library_id", library_id.as_str().to_owned()),
        ("stage", "amplification".to_owned()),
        ("family_ordinal", family_ordinal.to_string()),
        ("decision", decision.to_owned()),
    ];
    if let Some(attempt) = attempt {
        components.push(("attempt", attempt.to_string()));
    }
    components
}

fn f4_template_components(
    event_id: &EventId,
    library_id: &LibraryId,
    stage: &'static str,
    template_ordinal: u64,
    decision: &'static str,
    breakpoint: Option<u64>,
) -> Vec<(&'static str, String)> {
    let mut components = vec![
        ("event_id", event_id.as_str().to_owned()),
        ("library_id", library_id.as_str().to_owned()),
        ("stage", stage.to_owned()),
        ("template_ordinal", template_ordinal.to_string()),
        ("decision", decision.to_owned()),
    ];
    if let Some(breakpoint) = breakpoint {
        components.push(("breakpoint", breakpoint.to_string()));
    }
    components
}

fn f4_error_components(
    event_id: &EventId,
    library_id: &LibraryId,
    template_ordinal: u64,
    mate: F4Mate,
    cycle: u64,
    decision: &'static str,
) -> Vec<(&'static str, String)> {
    vec![
        ("event_id", event_id.as_str().to_owned()),
        ("library_id", library_id.as_str().to_owned()),
        ("stage", "error".to_owned()),
        ("template_ordinal", template_ordinal.to_string()),
        ("mate", mate.canonical_value().to_owned()),
        ("cycle", cycle.to_string()),
        ("decision", decision.to_owned()),
    ]
}

/// Componentes de la etapa `stutter`, en el orden literal de DTC-090:
/// `[event_id, library_id, stutter, template_ordinal, locus_id, decision]`.
fn f4_stutter_components(
    event_id: &EventId,
    library_id: &LibraryId,
    template_ordinal: u64,
    locus_id: &RepeatLocusId,
    decision: &'static str,
) -> Vec<(&'static str, String)> {
    vec![
        ("event_id", event_id.as_str().to_owned()),
        ("library_id", library_id.as_str().to_owned()),
        ("stage", "stutter".to_owned()),
        ("template_ordinal", template_ordinal.to_string()),
        ("locus_id", locus_id.as_str().to_owned()),
        ("decision", decision.to_owned()),
    ]
}

fn f4_contamination_components(
    event_id: &EventId,
    recipient_library_id: &LibraryId,
    stage: &'static str,
    recipient_id: &IndividualId,
    decision: &'static str,
) -> Vec<(&'static str, String)> {
    vec![
        ("event_id", event_id.as_str().to_owned()),
        (
            "recipient_library_id",
            recipient_library_id.as_str().to_owned(),
        ),
        ("stage", stage.to_owned()),
        ("recipient_id", recipient_id.as_str().to_owned()),
        ("decision", decision.to_owned()),
    ]
}

fn f4_index_hop_components(
    event_id: &EventId,
    pre_hop_library_id: &LibraryId,
    pre_hop_instance_id: F4IndexHopInstance<'_>,
    decision: &'static str,
) -> Vec<(&'static str, String)> {
    vec![
        ("event_id", event_id.as_str().to_owned()),
        ("pre_hop_library_id", pre_hop_library_id.as_str().to_owned()),
        ("stage", "index_hop".to_owned()),
        (
            "pre_hop_instance_id",
            pre_hop_instance_id.as_str().to_owned(),
        ),
        ("decision", decision.to_owned()),
    ]
}

/// Deriva la subsemilla exclusiva de una decisión F4.
///
/// El nombre del catálogo constituye el ámbito de derivación v2; los
/// componentes posteriores siguen el orden literal de DTC-086. Cada llamada
/// devuelve una nueva semilla ChaCha12 y nunca avanza un flujo compartido.
///
/// El ámbito lo elige **la decisión**, no el binario: las heredadas conservan
/// el v3 y las del tartamudeo nacen en el v4. Ver
/// [`F4_TECHNICAL_RNG_CATALOG_V4`].
pub fn derive_f4_technical_subseed(
    root: RootSeed,
    decision: F4TechnicalDecision<'_>,
) -> Result<Subseed, Diagnostic> {
    let catalog = decision.catalog();
    let components = decision.canonical_components();
    let components = components
        .iter()
        .map(|(name, value)| (*name, value.as_str()))
        .collect::<Vec<_>>();
    Ok(Subseed(blake3::derive_key(
        DERIVATION_CONTEXT,
        &canonical_material_v2(root, catalog, &components)?,
    )))
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum RandomSamplingError {
    InvalidBound,
    InvalidProbability,
    InvalidWeights,
    InvalidInterval,
    InvalidLambda,
    InvalidPloidy,
    InvalidDosage,
    InvalidInbreedingCoefficient,
    NumericFailure,
    InsufficientCandidates,
    InvalidMatrix,
    /// DT-006: un sorteo de la normal truncada F4 agotó
    /// [`F4_TRUNCATED_NORMAL_MAX_ATTEMPTS`] intentos sin acertar en el
    /// intervalo admitido.
    TruncatedNormalAttemptsExhausted,
}

/// Uniforme abierto normativo, compartido por los consumidores F2.
#[must_use]
pub fn open_unit_f64(rng: &mut LocalRng) -> f64 {
    open_unit_f64_with_source(rng)
}

fn open_unit_f64_with_source(rng: &mut impl RandomWordSource) -> f64 {
    (((rng.next_word() >> 12) as f64) + 0.5) * (1.0 / 4_503_599_627_370_496.0)
}

pub fn uniform_u64_below(rng: &mut LocalRng, n: u64) -> Result<u64, RandomSamplingError> {
    uniform_u64_below_with_source(rng, n)
}

/// Uniforme discreto F4 definido por DTC-088.
///
/// A diferencia del consumidor histórico [`uniform_u64_below`], esta primitiva
/// rechaza la cola superior de un espacio completo de 64 bits. Su contrato y
/// sus vectores no se comparten con fases anteriores.
pub fn f4_uniform_index(rng: &mut LocalRng, n: u64) -> Result<u64, RandomSamplingError> {
    f4_uniform_index_with_source(rng, n)
}

fn f4_uniform_index_with_source(
    rng: &mut impl RandomWordSource,
    n: u64,
) -> Result<u64, RandomSamplingError> {
    if n == 0 {
        return Err(RandomSamplingError::InvalidBound);
    }
    if n == 1 {
        return Ok(0);
    }
    // L = 2^64 - (2^64 mod n). Si el resto es cero, L no cabe en `u64` y
    // todo el espacio es aceptable. En otro caso, L sí es representable.
    let remainder = (u64::MAX % n + 1) % n;
    let limit = (remainder != 0).then(|| u64::MAX - (remainder - 1));
    loop {
        let value = rng.next_word();
        if limit.is_none_or(|bound| value < bound) {
            return Ok(value % n);
        }
    }
}

/// Uniforme abierto F4 con 53 bits efectivos, `(0, 1)`.
#[must_use]
pub fn f4_uniform_open_53(rng: &mut LocalRng) -> f64 {
    f4_uniform_open_53_with_source(rng)
}

fn f4_uniform_open_53_with_source(rng: &mut impl RandomWordSource) -> f64 {
    let value = f4_uniform_index_with_source(rng, F4_FIFTY_THREE_BIT_BOUND)
        .expect("2^53 is a valid fixed F4 uniform bound");
    (value as f64 + 0.5) / F4_FIFTY_THREE_BIT_BOUND as f64
}

/// Bernoulli F4. Las probabilidades exactas 0 y 1 no consumen aleatoriedad.
pub fn f4_bernoulli(rng: &mut LocalRng, probability: f64) -> Result<bool, RandomSamplingError> {
    f4_bernoulli_with_source(rng, probability)
}

fn f4_bernoulli_with_source(
    rng: &mut impl RandomWordSource,
    probability: f64,
) -> Result<bool, RandomSamplingError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(RandomSamplingError::InvalidProbability);
    }
    if probability == 0.0 {
        return Ok(false);
    }
    if probability == 1.0 {
        return Ok(true);
    }
    let numerator = f4_uniform_index_with_source(rng, F4_FIFTY_THREE_BIT_BOUND)?;
    Ok(numerator as f64 / (F4_FIFTY_THREE_BIT_BOUND as f64) < probability)
}

/// Uniforme cerrado F4 para un intervalo finito `[minimum, maximum]`.
pub fn f4_uniform_closed(
    rng: &mut LocalRng,
    minimum: f64,
    maximum: f64,
) -> Result<f64, RandomSamplingError> {
    f4_uniform_closed_with_source(rng, minimum, maximum)
}

fn f4_uniform_closed_with_source(
    rng: &mut impl RandomWordSource,
    minimum: f64,
    maximum: f64,
) -> Result<f64, RandomSamplingError> {
    if !minimum.is_finite() || !maximum.is_finite() || minimum > maximum {
        return Err(RandomSamplingError::InvalidInterval);
    }
    if minimum == maximum {
        return Ok(minimum);
    }
    let numerator = f4_uniform_index_with_source(rng, F4_FIFTY_THREE_BIT_BOUND)?;
    if numerator == 0 {
        return Ok(minimum);
    }
    if numerator == F4_FIFTY_THREE_BIT_MAX {
        return Ok(maximum);
    }
    let fraction = numerator as f64 / F4_FIFTY_THREE_BIT_MAX as f64;
    Ok((1.0 - fraction) * minimum + fraction * maximum)
}

/// Permuta sin reemplazo solo el prefijo solicitado de una lista canónica F4.
pub fn f4_partial_fisher_yates<T>(
    rng: &mut LocalRng,
    values: &mut [T],
    selected_count: usize,
) -> Result<(), RandomSamplingError> {
    f4_partial_fisher_yates_with_source(rng, values, selected_count)
}

fn f4_partial_fisher_yates_with_source<T>(
    rng: &mut impl RandomWordSource,
    values: &mut [T],
    selected_count: usize,
) -> Result<(), RandomSamplingError> {
    if selected_count > values.len() {
        return Err(RandomSamplingError::InsufficientCandidates);
    }
    for index in 0..selected_count {
        let remaining = values.len() - index;
        let offset = f4_uniform_index_with_source(
            rng,
            u64::try_from(remaining).map_err(|_| RandomSamplingError::NumericFailure)?,
        )?;
        let other = index
            .checked_add(usize::try_from(offset).map_err(|_| RandomSamplingError::NumericFailure)?)
            .ok_or(RandomSamplingError::NumericFailure)?;
        values.swap(index, other);
    }
    Ok(())
}

/// Cota de intentos de un único sorteo de la normal truncada F4 (DT-006).
///
/// `f4_truncated_normal` sortea por rechazo: genera un candidato y lo
/// descarta hasta que cae en el intervalo admitido. Cuando ese intervalo
/// queda hundido en la cola de la distribución, el número de intentos se
/// dispara. Medido sobre el perfil vigente (`gps-hiseqx-pcrfree-pe150-v1`,
/// media 350 pb, desviación 20 pb): un techo efectivo de 271 pb necesita del
/// orden de **23 000** intentos por sorteo, el peor caso legítimo medido; un
/// techo de 250 pb, el propio mínimo del perfil, necesita del orden de
/// **13 000 000** —dos parejas tardaron 111 segundos antes de que existiera
/// esta cota—. Un millón deja unas **cuarenta veces** de margen sobre el
/// peor caso legítimo conocido y sigue cazando el desbocamiento de trece
/// millones muy por debajo de su coste real.
const F4_TRUNCATED_NORMAL_MAX_ATTEMPTS: u32 = 1_000_000;

/// Normal truncada F4 de DTC-088, redondeada al entero más próximo.
///
/// Un único sorteo que agote [`F4_TRUNCATED_NORMAL_MAX_ATTEMPTS`] intentos
/// sin acertar en el intervalo admitido se aborta con
/// [`RandomSamplingError::TruncatedNormalAttemptsExhausted`] en vez de seguir
/// girando (DT-006). El corte no consume ningún número de la cinta de azar
/// que un intento no hubiera consumido ya: cuenta, no sortea.
pub fn f4_truncated_normal(
    rng: &mut LocalRng,
    mean: f64,
    standard_deviation: f64,
    minimum: i64,
    maximum: i64,
) -> Result<i64, RandomSamplingError> {
    f4_truncated_normal_with_source(rng, mean, standard_deviation, minimum, maximum)
}

fn f4_truncated_normal_with_source(
    rng: &mut impl RandomWordSource,
    mean: f64,
    standard_deviation: f64,
    minimum: i64,
    maximum: i64,
) -> Result<i64, RandomSamplingError> {
    if !mean.is_finite()
        || !standard_deviation.is_finite()
        || standard_deviation <= 0.0
        || minimum > maximum
    {
        return Err(RandomSamplingError::InvalidInterval);
    }
    let mut attempts: u32 = 0;
    loop {
        attempts += 1;
        if attempts > F4_TRUNCATED_NORMAL_MAX_ATTEMPTS {
            return Err(RandomSamplingError::TruncatedNormalAttemptsExhausted);
        }
        let first = f4_uniform_open_53_with_source(rng);
        let second = f4_uniform_open_53_with_source(rng);
        let z = libm::sqrt(-2.0 * libm::log(first)) * libm::cos(2.0 * F4_PI * second);
        let rounded = (mean + standard_deviation * z).round();
        if !rounded.is_finite() || rounded < i64::MIN as f64 || rounded > i64::MAX as f64 {
            return Err(RandomSamplingError::NumericFailure);
        }
        let candidate = rounded as i64;
        if (minimum..=maximum).contains(&candidate) {
            return Ok(candidate);
        }
    }
}

fn uniform_u64_below_with_source(
    rng: &mut impl RandomWordSource,
    n: u64,
) -> Result<u64, RandomSamplingError> {
    if n == 0 {
        return Err(RandomSamplingError::InvalidBound);
    }
    if n == 1 {
        return Ok(0);
    }
    let threshold = n.wrapping_neg() % n;
    loop {
        let value = rng.next_word();
        if value >= threshold {
            return Ok(value % n);
        }
    }
}

pub fn categorical_f64(rng: &mut LocalRng, weights: &[f64]) -> Result<usize, RandomSamplingError> {
    categorical_f64_with_source(rng, weights)
}

pub(crate) fn categorical_f64_with_source(
    rng: &mut impl RandomWordSource,
    weights: &[f64],
) -> Result<usize, RandomSamplingError> {
    if weights
        .iter()
        .any(|weight| !weight.is_finite() || *weight < 0.0)
    {
        return Err(RandomSamplingError::InvalidWeights);
    }
    let total = neumaier_sum(weights.iter().copied());
    if !total.is_finite() || total <= 0.0 {
        return Err(RandomSamplingError::InvalidWeights);
    }
    let positives = weights.iter().filter(|weight| **weight > 0.0).count();
    if positives == 1 {
        return Ok(weights
            .iter()
            .position(|weight| *weight > 0.0)
            .expect("positive count checked"));
    }
    let target = open_unit_f64_with_source(rng) * total;
    let mut cumulative = 0.0;
    let mut last_positive = 0;
    for (index, weight) in weights.iter().enumerate() {
        if *weight > 0.0 {
            last_positive = index;
        }
        cumulative += weight;
        if cumulative > target {
            return Ok(index);
        }
    }
    Ok(last_positive)
}

pub fn fisher_yates<T>(rng: &mut LocalRng, values: &mut [T]) -> Result<(), RandomSamplingError> {
    fisher_yates_with_source(rng, values)
}

fn fisher_yates_with_source<T>(
    rng: &mut impl RandomWordSource,
    values: &mut [T],
) -> Result<(), RandomSamplingError> {
    for index in (1..values.len()).rev() {
        let other = uniform_u64_below_with_source(rng, (index + 1) as u64)? as usize;
        values.swap(index, other);
    }
    Ok(())
}

/// Asigna una composición alélica completa a fases de copia mediante una
/// permutación normativa.
///
/// `allele_dosages` está en el orden canónico del locus: REF primero y ALT
/// ASCII ascendente después. La composición ya decide esas dosis; esta operación
/// crea la lista repetida de índices de alelo, la permuta una sola vez y asigna
/// su resultado secuencialmente a las fases `0..P`. Un locus monomórfico no
/// tiene elección de fase y no consume palabras RNG.
pub fn phase_permutation(
    rng: &mut LocalRng,
    allele_dosages: &[u8],
) -> Result<Vec<u8>, RandomSamplingError> {
    phase_permutation_with_source(rng, allele_dosages)
}

/// Variante instrumentable de [`phase_permutation`] para vectores normativos.
///
/// El preflight y los casos degenerados no consumen RNG. Un locus polimórfico
/// consume exactamente lo que consume `fisher-yates-v1` para `P` alelos,
/// incluidas sus repeticiones de rechazo si las hubiera.
pub fn phase_permutation_with_source(
    rng: &mut impl RandomWordSource,
    allele_dosages: &[u8],
) -> Result<Vec<u8>, RandomSamplingError> {
    if allele_dosages.is_empty() {
        return Err(RandomSamplingError::InvalidDosage);
    }
    let copy_count = allele_dosages.iter().try_fold(0_u8, |total, dosage| {
        total
            .checked_add(*dosage)
            .ok_or(RandomSamplingError::InvalidDosage)
    })?;
    if !(1..=10).contains(&copy_count) {
        return Err(RandomSamplingError::InvalidPloidy);
    }
    let mut alleles = Vec::with_capacity(usize::from(copy_count));
    for (allele_index, dosage) in allele_dosages.iter().copied().enumerate() {
        let allele_index =
            u8::try_from(allele_index).map_err(|_| RandomSamplingError::InvalidDosage)?;
        alleles.extend(std::iter::repeat_n(allele_index, usize::from(dosage)));
    }
    if alleles.iter().all(|allele| *allele == alleles[0]) {
        return Ok(alleles);
    }
    fisher_yates_with_source(rng, &mut alleles)?;
    Ok(alleles)
}

pub fn weighted_without_replacement(
    rng: &mut LocalRng,
    weights: &[f64],
    selected_count: usize,
) -> Result<Vec<usize>, RandomSamplingError> {
    weighted_without_replacement_with_source(rng, weights, selected_count)
}

fn weighted_without_replacement_with_source(
    rng: &mut impl RandomWordSource,
    weights: &[f64],
    selected_count: usize,
) -> Result<Vec<usize>, RandomSamplingError> {
    if selected_count == 0 {
        return Ok(Vec::new());
    }
    if weights.iter().filter(|weight| **weight > 0.0).count() < selected_count {
        return Err(RandomSamplingError::InsufficientCandidates);
    }
    let mut remaining = weights.to_vec();
    let mut selected = Vec::with_capacity(selected_count);
    for _ in 0..selected_count {
        let index = categorical_f64_with_source(rng, &remaining)?;
        selected.push(index);
        remaining[index] = 0.0;
    }
    Ok(selected)
}

pub fn linear_uniform(rng: &mut LocalRng, min: f64, max: f64) -> Result<f64, RandomSamplingError> {
    linear_uniform_with_source(rng, min, max)
}

fn linear_uniform_with_source(
    rng: &mut impl RandomWordSource,
    min: f64,
    max: f64,
) -> Result<f64, RandomSamplingError> {
    if !min.is_finite()
        || !max.is_finite()
        || min < 0.0
        || min >= max
        || max > 1.0
        || max <= f64::from_bits(1)
    {
        return Err(RandomSamplingError::InvalidInterval);
    }
    let mut value = min + open_unit_f64_with_source(rng) * (max - min);
    if value >= max {
        value = f64::from_bits(max.to_bits() - 1);
    }
    if value == 0.0 {
        value = f64::from_bits(1);
    }
    Ok(value)
}

pub fn log_uniform(rng: &mut LocalRng, min: f64, max: f64) -> Result<f64, RandomSamplingError> {
    log_uniform_with_source(rng, min, max)
}

fn log_uniform_with_source(
    rng: &mut impl RandomWordSource,
    min: f64,
    max: f64,
) -> Result<f64, RandomSamplingError> {
    if min <= 0.0 {
        return Err(RandomSamplingError::InvalidInterval);
    }
    linear_uniform_preflight(min, max)?;
    let mut value = libm::exp(
        libm::log(min) + open_unit_f64_with_source(rng) * (libm::log(max) - libm::log(min)),
    );
    if value >= max {
        value = f64::from_bits(max.to_bits() - 1);
    }
    if value < min {
        value = min;
    }
    Ok(value)
}

fn linear_uniform_preflight(min: f64, max: f64) -> Result<(), RandomSamplingError> {
    if !min.is_finite()
        || !max.is_finite()
        || min < 0.0
        || min >= max
        || max > 1.0
        || max <= f64::from_bits(1)
    {
        Err(RandomSamplingError::InvalidInterval)
    } else {
        Ok(())
    }
}

/// Suma compensada de Neumaier, el algoritmo `neumaier_sum_v1` que DTC-049
/// exige para normalizar pesos.
///
/// Se expone al resto de la biblioteca a propósito: es **la misma** que consume
/// el sorteo categórico normativo, de modo que quien valide un reparto de
/// frecuencias y quien después lo sortee usan exactamente la misma aritmética y
/// no pueden discrepar por el último bit.
pub(crate) fn neumaier_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        if sum.abs() >= value.abs() {
            correction += (sum - next) + value;
        } else {
            correction += (value - next) + sum;
        }
        sum = next;
    }
    sum + correction
}

/// Muestrea las dosis de cada alelo bajo HWE sin asignar fase a las copias.
///
/// Las frecuencias deben estar normalizadas y en el orden canónico del locus:
/// primero REF y después ALT. Cada copia es una categoría independiente; un
/// locus monomórfico no consume RNG.
pub fn hwe_composition(
    rng: &mut LocalRng,
    allele_frequencies: &[f64],
    copy_count: u8,
) -> Result<Vec<u8>, RandomSamplingError> {
    hwe_composition_with_source(rng, allele_frequencies, copy_count)
}

/// Variante instrumentable de [`hwe_composition`] para vectores normativos.
///
/// Consume una palabra por copia solo si el locus tiene más de una categoría
/// con peso positivo; las comprobaciones de entrada no consumen palabras.
pub fn hwe_composition_with_source(
    rng: &mut impl RandomWordSource,
    allele_frequencies: &[f64],
    copy_count: u8,
) -> Result<Vec<u8>, RandomSamplingError> {
    if !(1..=10).contains(&copy_count) {
        return Err(RandomSamplingError::InvalidPloidy);
    }
    validate_allele_frequencies(allele_frequencies)?;
    let mut dosage = vec![0_u8; allele_frequencies.len()];
    for _ in 0..copy_count {
        let allele = categorical_f64_with_source(rng, allele_frequencies)?;
        dosage[allele] = dosage[allele]
            .checked_add(1)
            .ok_or(RandomSamplingError::NumericFailure)?;
    }
    Ok(dosage)
}

/// Muestrea un genotipo diploide no ordenado con el coeficiente F normativo.
///
/// El resultado son dos índices de alelo sin fase. La enumeración es estable:
/// homocigotos `i,i` y después heterocigotos `i,j` con `i < j`.
pub fn f_diploid_genotype(
    rng: &mut LocalRng,
    allele_frequencies: &[f64],
    inbreeding_coefficient: f64,
) -> Result<[u8; 2], RandomSamplingError> {
    f_diploid_genotype_with_source(rng, allele_frequencies, inbreeding_coefficient)
}

/// Variante instrumentable de [`f_diploid_genotype`] para vectores normativos.
///
/// La enumeración, fórmulas y consumo son idénticos a la API de producción.
pub fn f_diploid_genotype_with_source(
    rng: &mut impl RandomWordSource,
    allele_frequencies: &[f64],
    inbreeding_coefficient: f64,
) -> Result<[u8; 2], RandomSamplingError> {
    if !inbreeding_coefficient.is_finite() || !(0.0..=1.0).contains(&inbreeding_coefficient) {
        return Err(RandomSamplingError::InvalidInbreedingCoefficient);
    }
    validate_allele_frequencies(allele_frequencies)?;
    let mut genotypes =
        Vec::with_capacity(allele_frequencies.len() * (allele_frequencies.len() + 1) / 2);
    let mut weights = Vec::with_capacity(genotypes.capacity());
    for (index, frequency) in allele_frequencies.iter().copied().enumerate() {
        let allele = u8::try_from(index).map_err(|_| RandomSamplingError::NumericFailure)?;
        genotypes.push([allele, allele]);
        weights
            .push(frequency * frequency + inbreeding_coefficient * frequency * (1.0 - frequency));
    }
    for (left, left_frequency) in allele_frequencies.iter().copied().enumerate() {
        for (right, right_frequency) in allele_frequencies
            .iter()
            .copied()
            .enumerate()
            .skip(left + 1)
        {
            genotypes.push([
                u8::try_from(left).map_err(|_| RandomSamplingError::NumericFailure)?,
                u8::try_from(right).map_err(|_| RandomSamplingError::NumericFailure)?,
            ]);
            weights.push(2.0 * (1.0 - inbreeding_coefficient) * left_frequency * right_frequency);
        }
    }
    Ok(genotypes[categorical_f64_with_source(rng, &weights)?])
}

fn validate_allele_frequencies(allele_frequencies: &[f64]) -> Result<(), RandomSamplingError> {
    if allele_frequencies.is_empty()
        || allele_frequencies
            .iter()
            .any(|frequency| !frequency.is_finite() || *frequency < 0.0)
    {
        return Err(RandomSamplingError::InvalidWeights);
    }
    let total = neumaier_sum(allele_frequencies.iter().copied());
    if !total.is_finite() || (total - 1.0).abs() > 1e-12 {
        return Err(RandomSamplingError::InvalidWeights);
    }
    Ok(())
}

pub fn poisson_count(rng: &mut LocalRng, lambda: f64) -> Result<u64, RandomSamplingError> {
    poisson_count_with_source(rng, lambda, None)
}

fn poisson_count_with_source(
    rng: &mut impl RandomWordSource,
    lambda: f64,
    candidate_override: Option<f64>,
) -> Result<u64, RandomSamplingError> {
    if !lambda.is_finite() || !(0.0..=POISSON_COUNT_MAX_LAMBDA).contains(&lambda) {
        return Err(RandomSamplingError::InvalidLambda);
    }
    if lambda == 0.0 {
        return Ok(0);
    }
    if lambda < 10.0 {
        let threshold = libm::exp(-lambda);
        let mut product = 1.0;
        let mut count = 0_u64;
        loop {
            product *= open_unit_f64_with_source(rng);
            if product <= threshold {
                return Ok(count);
            }
            count = count
                .checked_add(1)
                .ok_or(RandomSamplingError::NumericFailure)?;
        }
    }
    let sqrt_lambda = libm::sqrt(lambda);
    let b = 0.931 + 2.53 * sqrt_lambda;
    let a = -0.059 + 0.02483 * b;
    let inverse_alpha = 1.1239 + 1.1328 / (b - 3.4);
    let vr = 0.9277 - 3.6224 / (b - 2.0);
    let log_lambda = libm::log(lambda);
    loop {
        let u = open_unit_f64_with_source(rng) - 0.5;
        let v = open_unit_f64_with_source(rng);
        let us = 0.5 - u.abs();
        let candidate = candidate_override
            .unwrap_or_else(|| libm::floor((2.0 * a / us + b) * u + lambda + 0.43));
        if !(0.0..=9_007_199_254_740_991.0).contains(&candidate) {
            return Err(RandomSamplingError::NumericFailure);
        }
        if us >= 0.07 && v <= vr {
            return Ok(candidate as u64);
        }
        if us < 0.013 && v > us {
            continue;
        }
        let left = libm::log(v * inverse_alpha / (a / (us * us) + b));
        let right = -lambda + candidate * log_lambda - libm::lgamma(candidate + 1.0);
        if left <= right {
            return Ok(candidate as u64);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        F4_FIFTY_THREE_BIT_MAX, F4_TRUNCATED_NORMAL_MAX_ATTEMPTS, RandomSamplingError,
        RandomWordSource, categorical_f64_with_source, f4_bernoulli_with_source,
        f4_partial_fisher_yates_with_source, f4_truncated_normal_with_source,
        f4_uniform_closed_with_source, f4_uniform_index_with_source,
        f4_uniform_open_53_with_source, f5_bernoulli_with_source, f5_fisher_yates_with_source,
        f5_matrix_bytes, f5_placed_layout_bytes, f5_uniform_index_with_source,
        fisher_yates_with_source, linear_uniform_with_source, log_uniform_with_source,
        neumaier_sum, open_unit_f64_with_source, poisson_count_with_source,
        uniform_u64_below_with_source, weighted_without_replacement_with_source,
    };

    struct Scripted {
        words: Vec<u64>,
        consumed: usize,
    }

    impl Scripted {
        fn new(words: &[u64]) -> Self {
            Self {
                words: words.to_vec(),
                consumed: 0,
            }
        }
    }

    impl RandomWordSource for Scripted {
        fn next_word(&mut self) -> u64 {
            let value = self.words[self.consumed];
            self.consumed += 1;
            value
        }
    }

    /// Fuente sin fin para el corte de intentos de DT-006: cada intento
    /// produce el par `(0, 0)`, que con `mean=0.0, standard_deviation=1.0,
    /// minimum=1, maximum=1` cae siempre fuera del intervalo admitido, salvo
    /// en el intento `succeed_at_attempt` (si lo hay), donde produce el par
    /// `(half, 0)` que sí acierta. No asigna un vector gigante: calcula cada
    /// palabra a partir de cuántas lleva consumidas.
    struct ExhaustingSource {
        consumed: usize,
        succeed_at_attempt: Option<u32>,
    }

    impl ExhaustingSource {
        const HALF: u64 = 1_u64 << 52;

        fn never_succeeds() -> Self {
            Self {
                consumed: 0,
                succeed_at_attempt: None,
            }
        }

        fn succeeding_on_attempt(attempt: u32) -> Self {
            Self {
                consumed: 0,
                succeed_at_attempt: Some(attempt),
            }
        }
    }

    impl RandomWordSource for ExhaustingSource {
        fn next_word(&mut self) -> u64 {
            let word_index = self.consumed;
            self.consumed += 1;
            let attempt_is_one_based = (word_index / 2) as u32 + 1;
            let is_first_word_of_attempt = word_index.is_multiple_of(2);
            if self.succeed_at_attempt == Some(attempt_is_one_based) {
                if is_first_word_of_attempt {
                    Self::HALF
                } else {
                    0
                }
            } else {
                0
            }
        }
    }

    #[test]
    fn frozen_rng_utility_vectors_cover_values_rejection_and_consumption() {
        let mut open = Scripted::new(&[0]);
        assert_eq!(
            open_unit_f64_with_source(&mut open).to_bits(),
            0x3ca0_0000_0000_0000
        );
        assert_eq!(open.consumed, 1);

        let mut bounded = Scripted::new(&[0, 5]);
        assert_eq!(uniform_u64_below_with_source(&mut bounded, 3), Ok(2));
        assert_eq!(bounded.consumed, 2);
        assert_eq!(
            uniform_u64_below_with_source(&mut Scripted::new(&[]), 0),
            Err(RandomSamplingError::InvalidBound)
        );

        let mut categorical = Scripted::new(&[1 << 63]);
        assert_eq!(
            categorical_f64_with_source(&mut categorical, &[1.0, 2.0, 3.0]),
            Ok(2)
        );
        assert_eq!(categorical.consumed, 1);

        let mut shuffle = Scripted::new(&[1, 3]);
        let mut values = [0, 1, 2];
        assert_eq!(fisher_yates_with_source(&mut shuffle, &mut values), Ok(()));
        assert_eq!(values, [0, 2, 1]);
        assert_eq!(shuffle.consumed, 2);

        let mut weighted = Scripted::new(&[1 << 63, 0]);
        assert_eq!(
            weighted_without_replacement_with_source(&mut weighted, &[1.0, 2.0, 3.0], 2),
            Ok(vec![2, 0])
        );
        assert_eq!(weighted.consumed, 2);

        let mut linear = Scripted::new(&[0]);
        assert_eq!(
            linear_uniform_with_source(&mut linear, 0.0, 1.0).map(f64::to_bits),
            Ok(0x3ca0_0000_0000_0000)
        );
        assert_eq!(linear.consumed, 1);
        let mut logarithmic = Scripted::new(&[1 << 63]);
        assert_eq!(
            log_uniform_with_source(&mut logarithmic, 0.25, 1.0).map(f64::to_bits),
            Ok(0x3fe0_0000_0000_0001)
        );
        assert_eq!(logarithmic.consumed, 1);

        let mut knuth = Scripted::new(&[0]);
        assert_eq!(poisson_count_with_source(&mut knuth, 1.0, None), Ok(0));
        assert_eq!(knuth.consumed, 1);
        let mut ptrs = Scripted::new(&[1 << 63, 0]);
        assert_eq!(poisson_count_with_source(&mut ptrs, 10.0, None), Ok(10));
        assert_eq!(ptrs.consumed, 2);
    }

    #[test]
    fn f4_primitives_follow_dtc088_values_bounds_and_consumption() {
        let mut index = Scripted::new(&[u64::MAX, 5]);
        assert_eq!(f4_uniform_index_with_source(&mut index, 3), Ok(2));
        assert_eq!(index.consumed, 2);
        assert_eq!(
            f4_uniform_index_with_source(&mut Scripted::new(&[]), 0),
            Err(RandomSamplingError::InvalidBound)
        );
        let mut singleton = Scripted::new(&[]);
        assert_eq!(f4_uniform_index_with_source(&mut singleton, 1), Ok(0));
        assert_eq!(singleton.consumed, 0);

        let mut open = Scripted::new(&[0]);
        assert_eq!(
            f4_uniform_open_53_with_source(&mut open).to_bits(),
            (0.5_f64 / (1_u64 << 53) as f64).to_bits()
        );
        assert_eq!(open.consumed, 1);

        let mut bernoulli = Scripted::new(&[0]);
        assert_eq!(f4_bernoulli_with_source(&mut bernoulli, 0.5), Ok(true));
        assert_eq!(bernoulli.consumed, 1);
        let mut endpoints = Scripted::new(&[]);
        assert_eq!(f4_bernoulli_with_source(&mut endpoints, 0.0), Ok(false));
        assert_eq!(f4_bernoulli_with_source(&mut endpoints, 1.0), Ok(true));
        assert_eq!(endpoints.consumed, 0);
        assert_eq!(
            f4_bernoulli_with_source(&mut endpoints, f64::NAN),
            Err(RandomSamplingError::InvalidProbability)
        );

        let mut closed_lower = Scripted::new(&[0]);
        assert_eq!(
            f4_uniform_closed_with_source(&mut closed_lower, 1.0, 3.0),
            Ok(1.0)
        );
        assert_eq!(closed_lower.consumed, 1);
        let mut closed = Scripted::new(&[F4_FIFTY_THREE_BIT_MAX]);
        assert_eq!(
            f4_uniform_closed_with_source(&mut closed, 1.0, 3.0),
            Ok(3.0)
        );
        assert_eq!(closed.consumed, 1);
        let mut closed_singleton = Scripted::new(&[]);
        assert_eq!(
            f4_uniform_closed_with_source(&mut closed_singleton, 2.0, 2.0),
            Ok(2.0)
        );
        assert_eq!(closed_singleton.consumed, 0);

        let mut permutation = Scripted::new(&[3, 1]);
        let mut values = [0, 1, 2, 3];
        assert_eq!(
            f4_partial_fisher_yates_with_source(&mut permutation, &mut values, 2),
            Ok(())
        );
        assert_eq!(values, [3, 2, 1, 0]);
        assert_eq!(
            f4_partial_fisher_yates_with_source(&mut Scripted::new(&[]), &mut values, 5),
            Err(RandomSamplingError::InsufficientCandidates)
        );

        let half = 1_u64 << 52;
        let mut normal = Scripted::new(&[half, half, half, 0]);
        assert_eq!(0.5_f64.round(), 1.0, "ties round away from zero");
        assert_eq!((-0.5_f64).round(), -1.0, "ties round away from zero");
        assert_eq!(
            f4_truncated_normal_with_source(&mut normal, 0.0, 1.0, 1, 1),
            Ok(1)
        );
        assert_eq!(normal.consumed, 4);
        let mut invalid_normal = Scripted::new(&[]);
        assert_eq!(
            f4_truncated_normal_with_source(&mut invalid_normal, 0.0, 0.0, 1, 1),
            Err(RandomSamplingError::InvalidInterval)
        );
        assert_eq!(invalid_normal.consumed, 0);
    }

    /// DT-006: un sorteo que nunca acertaría se aborta exactamente en la cota
    /// de intentos, sin consumir ni un número de más.
    ///
    /// Regresión: sin el corte, esta llamada no terminaría nunca (el arnés
    /// que descubrió DT-006 midió 111 s para dos parejas con el mismo
    /// intervalo degenerado). La aserción sobre `consumed` es la que
    /// realmente ejercita el corte: sin él, la prueba no terminaría en vez de
    /// fallar, así que cualquier prueba de integración que solo comprobara
    /// `is_err()` bajo un límite de tiempo generoso no lo demostraría.
    #[test]
    fn f4_truncated_normal_aborts_at_the_attempt_cap_without_one_extra_draw() {
        let mut source = ExhaustingSource::never_succeeds();
        assert_eq!(
            f4_truncated_normal_with_source(&mut source, 0.0, 1.0, 1, 1),
            Err(RandomSamplingError::TruncatedNormalAttemptsExhausted)
        );
        // Exactamente dos palabras por intento, ni una palabra del intento
        // que habría seguido al que agotó la cota.
        assert_eq!(
            source.consumed,
            2 * F4_TRUNCATED_NORMAL_MAX_ATTEMPTS as usize
        );
    }

    /// El intento número `F4_TRUNCATED_NORMAL_MAX_ATTEMPTS`, el último
    /// permitido, todavía puede acertar: la cota no le roba su oportunidad.
    #[test]
    fn f4_truncated_normal_still_admits_a_success_on_the_last_allowed_attempt() {
        let mut source = ExhaustingSource::succeeding_on_attempt(F4_TRUNCATED_NORMAL_MAX_ATTEMPTS);
        assert_eq!(
            f4_truncated_normal_with_source(&mut source, 0.0, 1.0, 1, 1),
            Ok(1)
        );
        assert_eq!(
            source.consumed,
            2 * F4_TRUNCATED_NORMAL_MAX_ATTEMPTS as usize
        );
    }

    #[test]
    fn f5_primitives_follow_dtc089_rejection_consumption_and_canonical_bytes() {
        let mut index = Scripted::new(&[u64::MAX, 5]);
        assert_eq!(f5_uniform_index_with_source(&mut index, 3), Ok(2));
        assert_eq!(index.consumed, 2);
        let mut invalid_index = Scripted::new(&[]);
        assert_eq!(
            f5_uniform_index_with_source(&mut invalid_index, 0),
            Err(RandomSamplingError::InvalidBound)
        );
        assert_eq!(invalid_index.consumed, 0);
        let mut singleton = Scripted::new(&[]);
        assert_eq!(f5_uniform_index_with_source(&mut singleton, 1), Ok(0));
        assert_eq!(singleton.consumed, 0);

        let mut bernoulli = Scripted::new(&[0, u64::MAX]);
        assert_eq!(f5_bernoulli_with_source(&mut bernoulli, 0.5), Ok(true));
        assert_eq!(f5_bernoulli_with_source(&mut bernoulli, 0.5), Ok(false));
        assert_eq!(bernoulli.consumed, 2);
        let mut endpoints = Scripted::new(&[]);
        assert_eq!(f5_bernoulli_with_source(&mut endpoints, 0.0), Ok(false));
        assert_eq!(f5_bernoulli_with_source(&mut endpoints, 1.0), Ok(true));
        assert_eq!(
            f5_bernoulli_with_source(&mut endpoints, f64::NAN),
            Err(RandomSamplingError::InvalidProbability)
        );
        assert_eq!(endpoints.consumed, 0);

        let mut permutation = Scripted::new(&[1, 0, 1]);
        let mut values = [0, 1, 2, 3];
        assert_eq!(
            f5_fisher_yates_with_source(&mut permutation, &mut values),
            Ok(())
        );
        assert_eq!(values, [2, 3, 0, 1]);
        assert_eq!(permutation.consumed, 3);
        let mut trivial = Scripted::new(&[]);
        let mut one = [42];
        assert_eq!(f5_fisher_yates_with_source(&mut trivial, &mut one), Ok(()));
        assert_eq!(trivial.consumed, 0);

        assert_eq!(
            f5_placed_layout_bytes(&[(1, 2, 5, 4), (255, 256, 511, 256)]),
            [
                0, 0, 0, 0, 0, 0, 0, 1, //
                0, 0, 0, 0, 0, 0, 0, 2, //
                0, 0, 0, 0, 0, 0, 0, 5, //
                0, 0, 0, 0, 0, 0, 0, 4, //
                0, 0, 0, 0, 0, 0, 0, 255, //
                0, 0, 0, 0, 0, 0, 1, 0, //
                0, 0, 0, 0, 0, 0, 1, 255, //
                0, 0, 0, 0, 0, 0, 1, 0,
            ]
        );
        assert_eq!(
            f5_matrix_bytes(&[
                vec![true, false, true, false, false, false, false, true, true],
                vec![false, true, false, true, false, true, false, true, false],
            ]),
            Ok(vec![
                0, 0, 0, 0, 0, 0, 0, 2, //
                0, 0, 0, 0, 0, 0, 0, 9, //
                0xa1, 0x80, 0x55, 0x00,
            ])
        );
        assert_eq!(
            f5_matrix_bytes(&[vec![true], vec![false, true]]),
            Err(RandomSamplingError::InvalidMatrix)
        );
    }

    #[test]
    fn categorical_uses_neumaier_total_for_adversarial_canonical_weights() {
        let weights = [1e16, 1.0, 1.0];
        assert_eq!(weights.iter().sum::<f64>(), 1e16);
        assert_eq!(neumaier_sum(weights), 10_000_000_000_000_002.0);
        let mut rng = Scripted::new(&[u64::MAX]);
        assert_eq!(categorical_f64_with_source(&mut rng, &weights), Ok(2));
        assert_eq!(rng.consumed, 1);
    }

    #[test]
    fn poisson_frozen_boundary_vectors_have_exact_output_and_consumption() {
        for lambda in [f64::from_bits(1), 10.0_f64.next_down(), 1.0, 5.0] {
            let mut rng = Scripted::new(&[0]);
            assert_eq!(poisson_count_with_source(&mut rng, lambda, None), Ok(0));
            assert_eq!(rng.consumed, 1, "lambda={lambda}");
        }
        for lambda in [
            10.0,
            100.0,
            4_503_599_627_370_496.0_f64.next_down(),
            4_503_599_627_370_496.0,
        ] {
            let mut rng = Scripted::new(&[1 << 63, 0]);
            let expected = if lambda == 4_503_599_627_370_496.0_f64.next_down() {
                4_503_599_627_370_496
            } else {
                lambda as u64
            };
            assert_eq!(
                poisson_count_with_source(&mut rng, lambda, None),
                Ok(expected)
            );
            assert_eq!(rng.consumed, 2, "lambda={lambda}");
        }
    }

    #[test]
    fn ptrs_ordinary_rejection_then_acceptance_consumes_four_words() {
        let mut rng = Scripted::new(&[1 << 63, u64::MAX, 1 << 63, 0]);
        assert_eq!(poisson_count_with_source(&mut rng, 10.0, None), Ok(10));
        assert_eq!(rng.consumed, 4);
    }

    #[test]
    fn ptrs_out_of_domain_candidate_fails_after_its_two_uniforms() {
        let mut observed = Scripted::new(&[1 << 63, 1 << 63]);
        assert_eq!(
            poisson_count_with_source(&mut observed, 10.0, Some(f64::NAN)),
            Err(RandomSamplingError::NumericFailure)
        );
        assert_eq!(observed.consumed, 2);
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReproducibilityMetadata {
    root_seed: RootSeed,
    seed_origin: SeedOrigin,
    random_algorithm: &'static str,
    derivation_versions: &'static [u32],
    catalog_version: &'static str,
}

impl ReproducibilityMetadata {
    #[must_use]
    pub const fn new(root_seed: RootSeed, seed_origin: SeedOrigin) -> Self {
        Self {
            root_seed,
            seed_origin,
            random_algorithm: RANDOM_ALGORITHM,
            derivation_versions: SUBSEED_DERIVATION_VERSIONS,
            catalog_version: SUBSEED_CATALOG_VERSION,
        }
    }

    #[must_use]
    pub const fn root_seed(&self) -> RootSeed {
        self.root_seed
    }

    #[must_use]
    pub const fn seed_origin(&self) -> SeedOrigin {
        self.seed_origin
    }

    #[must_use]
    pub const fn random_algorithm(&self) -> &'static str {
        self.random_algorithm
    }

    #[must_use]
    pub const fn derivation_versions(&self) -> &'static [u32] {
        self.derivation_versions
    }

    #[must_use]
    pub const fn catalog_version(&self) -> &'static str {
        self.catalog_version
    }
}
