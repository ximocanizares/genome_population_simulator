//! Semillas, flujos aleatorios locales y contratos de reproducibilidad.

use std::fmt;

use rand_chacha::ChaCha12Rng;
use rand_chacha::rand_core::{Rng, SeedableRng};

use crate::diagnostics::Diagnostic;

pub const DERIVATION_VERSION: u32 = 1;
pub const SUBSEED_CATALOG_VERSION: &str = "1";
pub const DERIVATION_CONTEXT: &str = "genome_population_simulator/v1/subseed";
pub const RANDOM_ALGORITHM: &str = "ChaCha12";

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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SubseedRequest {
    scope: &'static str,
    components: Vec<(&'static str, String)>,
}

impl SubseedRequest {
    pub fn from_scope(scope: &str, components: &[(&str, &str)]) -> Result<Self, Diagnostic> {
        let (canonical_scope, names): (&'static str, &[&'static str]) = match scope {
            "initial_population" => ("initial_population", &["scenario_id", "population_id"]),
            "biological_process" => (
                "biological_process",
                &["generation", "population_id", "process"],
            ),
            "sampling" => ("sampling", &["generation", "sampling_id"]),
            "technical_event" => ("technical_event", &["event_id", "reference_id"]),
            _ => return Err(unknown_scope(scope)),
        };
        if components.len() != names.len()
            || components
                .iter()
                .zip(names)
                .any(|((found, _), expected)| found != expected)
        {
            return Err(invalid_scope_components(scope, components, names));
        }
        Ok(Self {
            scope: canonical_scope,
            components: names
                .iter()
                .zip(components)
                .map(|(name, (_, value))| (*name, (*value).to_owned()))
                .collect(),
        })
    }

    #[must_use]
    pub const fn scope(&self) -> &'static str {
        self.scope
    }

    pub fn components(&self) -> impl Iterator<Item = (&'static str, &str)> {
        self.components
            .iter()
            .map(|(name, value)| (*name, value.as_str()))
    }
}

fn unknown_scope(scope: &str) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-018",
        "El ámbito solicitado no existe en el catálogo normativo de subsemillas.",
    )
    .with("scope", scope)
    .with("catalog_version", SUBSEED_CATALOG_VERSION)
}

fn invalid_scope_components(
    scope: &str,
    components: &[(&str, &str)],
    expected: &[&str],
) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "subseed_derivation")
    .with("operation", "construct_scope_request")
    .with(
        "detail",
        format!(
            "scope {scope} expected components {expected:?}, found {:?}",
            components.iter().map(|(name, _)| *name).collect::<Vec<_>>()
        ),
    )
}

pub fn canonical_material(root: RootSeed, request: &SubseedRequest) -> Result<Vec<u8>, Diagnostic> {
    let component_count = u32::try_from(request.components.len())
        .map_err(|_| material_length_error("component_count"))?;
    let mut material = Vec::new();
    material.extend_from_slice(root.as_bytes());
    material.extend_from_slice(&DERIVATION_VERSION.to_be_bytes());
    material.extend_from_slice(&component_count.to_be_bytes());
    for (name, value) in &request.components {
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

pub fn derive_subseed(root: RootSeed, request: &SubseedRequest) -> Result<Subseed, Diagnostic> {
    let material = canonical_material(root, request)?;
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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EventSubseed {
    event_id: String,
    subseed: Subseed,
}

impl EventSubseed {
    #[must_use]
    pub fn event_id(&self) -> &str {
        &self.event_id
    }

    #[must_use]
    pub const fn subseed(&self) -> Subseed {
        self.subseed
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReproducibilityMetadata {
    root_seed: RootSeed,
    seed_origin: SeedOrigin,
    random_algorithm: &'static str,
    derivation_version: u32,
    catalog_version: &'static str,
    event_subseeds: Vec<EventSubseed>,
}

impl ReproducibilityMetadata {
    pub fn new(
        root_seed: RootSeed,
        seed_origin: SeedOrigin,
        reference_id: &str,
        event_ids: &[String],
    ) -> Result<Self, Diagnostic> {
        let mut event_subseeds = event_ids
            .iter()
            .map(|event_id| {
                let request = SubseedRequest::from_scope(
                    "technical_event",
                    &[("event_id", event_id), ("reference_id", reference_id)],
                )?;
                Ok(EventSubseed {
                    event_id: event_id.clone(),
                    subseed: derive_subseed(root_seed, &request)?,
                })
            })
            .collect::<Result<Vec<_>, Diagnostic>>()?;
        event_subseeds.sort_by(|first, second| first.event_id.cmp(&second.event_id));
        Ok(Self {
            root_seed,
            seed_origin,
            random_algorithm: RANDOM_ALGORITHM,
            derivation_version: DERIVATION_VERSION,
            catalog_version: SUBSEED_CATALOG_VERSION,
            event_subseeds,
        })
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
    pub const fn derivation_version(&self) -> u32 {
        self.derivation_version
    }

    #[must_use]
    pub const fn catalog_version(&self) -> &'static str {
        self.catalog_version
    }

    #[must_use]
    pub fn event_subseeds(&self) -> &[EventSubseed] {
        &self.event_subseeds
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GuaranteeStatus {
    Guaranteed,
    NotGuaranteed,
    NotApplicable,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GuaranteeAssessment {
    status: GuaranteeStatus,
    conditions: &'static [&'static str],
}

impl GuaranteeAssessment {
    #[must_use]
    pub const fn status(&self) -> GuaranteeStatus {
        self.status
    }

    #[must_use]
    pub const fn conditions(&self) -> &'static [&'static str] {
        self.conditions
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReproducibilityContext {
    pub simulator_version: String,
    pub platform: String,
    pub root_seed: RootSeed,
    pub effective_config: String,
    pub input_identity: String,
    pub profile: String,
    pub random_algorithm: String,
    pub derivation_version: u32,
    pub catalog_version: String,
    /// Dependencias que intervienen en la derivación o en los flujos RNG.
    pub dependency_versions: Vec<String>,
    /// Dependencias adicionales que intervienen únicamente en resultados técnicos.
    pub technical_dependency_versions: Vec<String>,
    pub external_tools: Vec<String>,
}

/// Evidencia previa y concreta que permite comparar bytes de una clase de
/// artefacto. Su presencia no declara por sí sola que una salida real cumpla
/// REQ-CONF-115; esa declaración corresponde al manifiesto y al validador del
/// artefacto consumidor.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BinaryIdentityEvidence {
    artifact_class: String,
    verification_policy: String,
    evidence_id: String,
}

impl BinaryIdentityEvidence {
    #[must_use]
    pub fn new(
        artifact_class: impl Into<String>,
        verification_policy: impl Into<String>,
        evidence_id: impl Into<String>,
    ) -> Option<Self> {
        let evidence = Self {
            artifact_class: artifact_class.into(),
            verification_policy: verification_policy.into(),
            evidence_id: evidence_id.into(),
        };
        (!evidence.artifact_class.trim().is_empty()
            && !evidence.verification_policy.trim().is_empty()
            && !evidence.evidence_id.trim().is_empty())
        .then_some(evidence)
    }

    #[must_use]
    pub fn artifact_class(&self) -> &str {
        &self.artifact_class
    }

    #[must_use]
    pub fn verification_policy(&self) -> &str {
        &self.verification_policy
    }

    #[must_use]
    pub fn evidence_id(&self) -> &str {
        &self.evidence_id
    }
}

/// Entrada preparatoria para comparar hechos ya validados por consumidores.
/// No sustituye un resultado, un manifiesto ni evidencia sobre un artefacto real.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReproducibilityObservation {
    pub context: ReproducibilityContext,
    pub biological_decisions: Option<String>,
    pub technical_content: Option<String>,
    pub semantic_content: Option<String>,
    pub binary_content: Option<Vec<u8>>,
    pub event_subseeds: Vec<EventSubseed>,
    pub binary_identity_evidence: Option<BinaryIdentityEvidence>,
}

fn base_conditions_match(first: &ReproducibilityContext, second: &ReproducibilityContext) -> bool {
    first.simulator_version == second.simulator_version
        && first.platform == second.platform
        && first.root_seed == second.root_seed
        && first.effective_config == second.effective_config
        && first.input_identity == second.input_identity
        && first.profile == second.profile
        && first.random_algorithm == second.random_algorithm
        && first.derivation_version == second.derivation_version
        && first.catalog_version == second.catalog_version
        && same_unordered(&first.dependency_versions, &second.dependency_versions)
}

fn same_unordered<T: PartialEq>(first: &[T], second: &[T]) -> bool {
    if first.len() != second.len() {
        return false;
    }
    let mut matched = vec![false; second.len()];
    first.iter().all(|value| {
        second
            .iter()
            .enumerate()
            .find(|(index, candidate)| !matched[*index] && value == *candidate)
            .is_some_and(|(index, _)| {
                matched[index] = true;
                true
            })
    })
}

pub fn compare_biological(
    first: &ReproducibilityObservation,
    second: &ReproducibilityObservation,
) -> GuaranteeAssessment {
    let status = match (&first.biological_decisions, &second.biological_decisions) {
        (Some(first_value), Some(second_value))
            if base_conditions_match(&first.context, &second.context)
                && first_value == second_value =>
        {
            GuaranteeStatus::Guaranteed
        }
        (Some(_), Some(_)) | (Some(_), None) | (None, Some(_)) => GuaranteeStatus::NotGuaranteed,
        (None, None) => GuaranteeStatus::NotApplicable,
    };
    GuaranteeAssessment {
        status,
        conditions: &[
            "simulator_version",
            "platform",
            "root_seed",
            "effective_config",
            "input_identity",
            "profile",
            "random_algorithm",
            "derivation_version",
            "catalog_version",
            "dependency_versions",
            "biological_decisions",
        ],
    }
}

pub fn compare_technical(
    first: &ReproducibilityObservation,
    second: &ReproducibilityObservation,
) -> GuaranteeAssessment {
    let status = match (&first.technical_content, &second.technical_content) {
        (Some(first_value), Some(second_value))
            if base_conditions_match(&first.context, &second.context)
                && same_unordered(
                    &first.context.technical_dependency_versions,
                    &second.context.technical_dependency_versions,
                )
                && same_unordered(&first.event_subseeds, &second.event_subseeds)
                && first_value == second_value =>
        {
            GuaranteeStatus::Guaranteed
        }
        (Some(_), Some(_)) | (Some(_), None) | (None, Some(_)) => GuaranteeStatus::NotGuaranteed,
        (None, None) => GuaranteeStatus::NotApplicable,
    };
    GuaranteeAssessment {
        status,
        conditions: &[
            "biological_base_conditions",
            "event_subseeds",
            "technical_dependency_versions",
            "technical_content",
        ],
    }
}

pub fn compare_semantic(
    first: &ReproducibilityObservation,
    second: &ReproducibilityObservation,
) -> GuaranteeAssessment {
    let status = match (&first.semantic_content, &second.semantic_content) {
        (Some(first_value), Some(second_value))
            if base_conditions_match(&first.context, &second.context)
                && same_unordered(
                    &first.context.external_tools,
                    &second.context.external_tools,
                )
                && first_value == second_value =>
        {
            GuaranteeStatus::Guaranteed
        }
        (Some(_), Some(_)) | (Some(_), None) | (None, Some(_)) => GuaranteeStatus::NotGuaranteed,
        (None, None) => GuaranteeStatus::NotApplicable,
    };
    GuaranteeAssessment {
        status,
        conditions: &[
            "biological_base_conditions",
            "external_tools",
            "semantic_content",
        ],
    }
}

pub fn compare_binary(
    first: &ReproducibilityObservation,
    second: &ReproducibilityObservation,
) -> GuaranteeAssessment {
    let status = match (&first.binary_content, &second.binary_content) {
        (Some(first_value), Some(second_value))
            if first.binary_identity_evidence.is_some()
                && first.binary_identity_evidence == second.binary_identity_evidence
                && base_conditions_match(&first.context, &second.context)
                && same_unordered(
                    &first.context.technical_dependency_versions,
                    &second.context.technical_dependency_versions,
                )
                && same_unordered(
                    &first.context.external_tools,
                    &second.context.external_tools,
                )
                && same_unordered(&first.event_subseeds, &second.event_subseeds)
                && first_value == second_value =>
        {
            GuaranteeStatus::Guaranteed
        }
        (Some(_), Some(_)) | (Some(_), None) | (None, Some(_)) => GuaranteeStatus::NotGuaranteed,
        (None, None) => GuaranteeStatus::NotApplicable,
    };
    GuaranteeAssessment {
        status,
        conditions: &[
            "binary_artifact_class",
            "binary_verification_policy",
            "binary_verification_evidence",
            "biological_base_conditions",
            "event_subseeds",
            "technical_dependency_versions",
            "external_tools",
            "binary_content",
        ],
    }
}
