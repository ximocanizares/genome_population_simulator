//! Identidades tipadas, etiquetas humanas y derivación automática versionada.

use std::collections::BTreeSet;
use std::fmt;

use crate::diagnostics::Diagnostic;
use crate::randomness::{DERIVATION_VERSION, RootSeed};

pub const ID_CONTEXT: &str = "genome_population_simulator/v1/id";
pub const ID_CATALOG_VERSION: &str = "1";

macro_rules! typed_id {
    ($name:ident, $category:literal) => {
        #[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
        pub struct $name(String);

        impl $name {
            pub fn parse(value: impl Into<String>, scope: &str) -> Result<Self, Diagnostic> {
                let value = value.into();
                validate_provided_id(&value, $category, scope)?;
                Ok(Self(value))
            }

            fn from_automatic(value: String) -> Self {
                Self(value)
            }

            #[must_use]
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(&self.0)
            }
        }
    };
}

typed_id!(ExecutionId, "execution");
typed_id!(ReferenceId, "reference");
typed_id!(IndividualId, "individual");
typed_id!(EventId, "event");

macro_rules! provided_typed_id {
    ($name:ident, $category:literal) => {
        #[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
        pub struct $name(String);

        impl $name {
            pub fn parse(value: impl Into<String>, scope: &str) -> Result<Self, Diagnostic> {
                let value = value.into();
                validate_provided_id(&value, $category, scope)?;
                Ok(Self(value))
            }

            #[must_use]
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(&self.0)
            }
        }
    };
}

provided_typed_id!(SequenceId, "sequence");
provided_typed_id!(ScenarioId, "scenario");
provided_typed_id!(PopulationId, "population");
provided_typed_id!(SubpopulationId, "subpopulation");
provided_typed_id!(PhaseId, "phase");
provided_typed_id!(GroupId, "group");
provided_typed_id!(BlockId, "block");
provided_typed_id!(CrossId, "cross");

/// Etiqueta Unicode opcional y sin autoridad de identidad.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HumanLabel(String);

impl HumanLabel {
    #[must_use]
    pub fn new(value: impl Into<String>) -> Self {
        Self(value.into())
    }

    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

fn validate_provided_id(value: &str, category: &str, scope: &str) -> Result<(), Diagnostic> {
    let is_valid = (1..=64).contains(&value.len())
        && value.is_ascii()
        && value
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_alphanumeric)
        && value.as_bytes()[1..]
            .iter()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(*byte, b'.' | b'_' | b'-'));
    if is_valid {
        return Ok(());
    }
    Err(Diagnostic::domain(
        "GPS-DOMAIN-001",
        "Un identificador aportado incumple la sintaxis o el ámbito de su categoría.",
    )
    .with("identifier", value)
    .with("category", category)
    .with("scope", scope))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IdClass {
    Execution,
    Reference,
    Individual,
    Event,
}

impl IdClass {
    const fn name(self) -> &'static str {
        match self {
            Self::Execution => "execution",
            Self::Reference => "reference",
            Self::Individual => "individual",
            Self::Event => "event",
        }
    }

    const fn prefix(self) -> &'static str {
        match self {
            Self::Execution => "RUN_",
            Self::Reference => "REF_",
            Self::Individual => "IND_",
            Self::Event => "EVT_",
        }
    }

    const fn cause_names(self) -> &'static [&'static str] {
        match self {
            Self::Execution => &["config_sha256"],
            Self::Reference => &["source_kind", "source_digest"],
            Self::Individual => &["founder_key"],
            Self::Event => &["event_key"],
        }
    }
}

/// Petición que solo puede construirse con el orden literal del catálogo A6.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AutomaticIdRequest {
    class: IdClass,
    scope: String,
    causes: Vec<(&'static str, String)>,
}

impl AutomaticIdRequest {
    pub fn new(
        class: IdClass,
        scope: impl Into<String>,
        causes: &[(&str, &str)],
    ) -> Result<Self, Diagnostic> {
        let expected = class.cause_names();
        if causes.len() != expected.len()
            || causes
                .iter()
                .zip(expected)
                .any(|((found, _), expected)| found != expected)
        {
            return Err(Diagnostic::internal(
                "GPS-INTERNAL-001",
                "Una invariante interna considerada imposible se ha roto.",
            )
            .with("component", "identity")
            .with("operation", "construct_automatic_id_request")
            .with("detail", format!("expected causes {expected:?}")));
        }
        Ok(Self {
            class,
            scope: scope.into(),
            causes: expected
                .iter()
                .zip(causes)
                .map(|(name, (_, value))| (*name, (*value).to_owned()))
                .collect(),
        })
    }

    pub fn canonical_material(&self, root: RootSeed) -> Result<Vec<u8>, Diagnostic> {
        let component_count = u32::try_from(self.causes.len() + 2)
            .map_err(|_| identity_material_error("component_count"))?;
        let mut material = Vec::new();
        material.extend_from_slice(root.as_bytes());
        material.extend_from_slice(&DERIVATION_VERSION.to_be_bytes());
        material.extend_from_slice(&component_count.to_be_bytes());
        append_component(&mut material, "class", self.class.name())?;
        append_component(&mut material, "scope", &self.scope)?;
        for (name, value) in &self.causes {
            append_component(&mut material, name, value)?;
        }
        Ok(material)
    }

    fn derive_text(&self, root: RootSeed) -> Result<String, Diagnostic> {
        let digest = blake3::derive_key(ID_CONTEXT, &self.canonical_material(root)?);
        Ok(format!(
            "{}{}",
            self.class.prefix(),
            encode_base32(&digest[..16])
        ))
    }
}

fn append_component(destination: &mut Vec<u8>, name: &str, value: &str) -> Result<(), Diagnostic> {
    append_length_prefixed(destination, name.as_bytes(), "component_name")?;
    append_length_prefixed(destination, value.as_bytes(), "component_value")
}

fn append_length_prefixed(
    destination: &mut Vec<u8>,
    value: &[u8],
    field: &str,
) -> Result<(), Diagnostic> {
    let length = u32::try_from(value.len()).map_err(|_| identity_material_error(field))?;
    destination.extend_from_slice(&length.to_be_bytes());
    destination.extend_from_slice(value);
    Ok(())
}

fn identity_material_error(field: &str) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "identity")
    .with("operation", "encode_canonical_material")
    .with("detail", format!("{field} exceeds u32"))
}

fn encode_base32(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 32] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
    let mut encoded = String::with_capacity(bytes.len().div_ceil(5) * 8);
    let mut buffer = 0_u32;
    let mut bits = 0_u8;
    for byte in bytes {
        buffer = (buffer << 8) | u32::from(*byte);
        bits += 8;
        while bits >= 5 {
            bits -= 5;
            encoded.push(char::from(ALPHABET[((buffer >> bits) & 0x1f) as usize]));
        }
    }
    if bits > 0 {
        encoded.push(char::from(
            ALPHABET[((buffer << (5 - bits)) & 0x1f) as usize],
        ));
    }
    encoded
}

/// Registro de unicidad por categoría y ámbito. Las categorías no se mezclan.
#[derive(Debug, Default)]
pub struct IdentityRegistry {
    registered: BTreeSet<(IdClassKey, String, String)>,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum IdClassKey {
    Execution,
    Reference,
    Individual,
    Event,
}

impl From<IdClass> for IdClassKey {
    fn from(value: IdClass) -> Self {
        match value {
            IdClass::Execution => Self::Execution,
            IdClass::Reference => Self::Reference,
            IdClass::Individual => Self::Individual,
            IdClass::Event => Self::Event,
        }
    }
}

impl IdentityRegistry {
    fn register_text(
        &mut self,
        class: IdClass,
        scope: &str,
        identifier: &str,
        causes: &str,
    ) -> Result<(), Diagnostic> {
        if self.registered.insert((
            IdClassKey::from(class),
            scope.to_owned(),
            identifier.to_owned(),
        )) {
            return Ok(());
        }
        Err(Diagnostic::domain(
            "GPS-DOMAIN-002",
            "La derivación de una identidad produjo una colisión.",
        )
        .with("identifier", identifier)
        .with("category", class.name())
        .with("causes", causes))
    }

    pub fn register_individual(
        &mut self,
        scope: &str,
        id: IndividualId,
    ) -> Result<IndividualId, Diagnostic> {
        self.register_text(IdClass::Individual, scope, id.as_str(), "provided")?;
        Ok(id)
    }

    pub fn register_event(&mut self, scope: &str, id: EventId) -> Result<EventId, Diagnostic> {
        self.register_text(IdClass::Event, scope, id.as_str(), "provided")?;
        Ok(id)
    }

    pub fn derive_execution(
        &mut self,
        root: RootSeed,
        request: &AutomaticIdRequest,
    ) -> Result<ExecutionId, Diagnostic> {
        self.require_class(request, IdClass::Execution)?;
        let text = request.derive_text(root)?;
        self.register_text(request.class, &request.scope, &text, "automatic")?;
        Ok(ExecutionId::from_automatic(text))
    }

    pub fn derive_reference(
        &mut self,
        root: RootSeed,
        request: &AutomaticIdRequest,
    ) -> Result<ReferenceId, Diagnostic> {
        self.require_class(request, IdClass::Reference)?;
        let text = request.derive_text(root)?;
        self.register_text(request.class, &request.scope, &text, "automatic")?;
        Ok(ReferenceId::from_automatic(text))
    }

    pub fn derive_individual(
        &mut self,
        root: RootSeed,
        request: &AutomaticIdRequest,
    ) -> Result<IndividualId, Diagnostic> {
        self.require_class(request, IdClass::Individual)?;
        let text = request.derive_text(root)?;
        self.register_text(request.class, &request.scope, &text, "automatic")?;
        Ok(IndividualId::from_automatic(text))
    }

    pub fn derive_event(
        &mut self,
        root: RootSeed,
        request: &AutomaticIdRequest,
    ) -> Result<EventId, Diagnostic> {
        self.require_class(request, IdClass::Event)?;
        let text = request.derive_text(root)?;
        self.register_text(request.class, &request.scope, &text, "automatic")?;
        Ok(EventId::from_automatic(text))
    }

    fn require_class(
        &self,
        request: &AutomaticIdRequest,
        expected: IdClass,
    ) -> Result<(), Diagnostic> {
        if request.class == expected {
            return Ok(());
        }
        Err(Diagnostic::internal(
            "GPS-INTERNAL-001",
            "Una invariante interna considerada imposible se ha roto.",
        )
        .with("component", "identity")
        .with("operation", "derive_typed_id")
        .with(
            "detail",
            "automatic request category does not match return type",
        ))
    }
}

/// Deriva la identidad de ejecución con el único material canónico de F1.
pub fn derive_execution_id(
    root: RootSeed,
    configuration_sha256: &str,
) -> Result<ExecutionId, Diagnostic> {
    let request = AutomaticIdRequest::new(
        IdClass::Execution,
        "execution",
        &[("config_sha256", configuration_sha256)],
    )?;
    IdentityRegistry::default().derive_execution(root, &request)
}

/// Las categorías no son intercambiables aunque el texto sea igual.
///
/// ```compile_fail
/// use genome_population_simulator::identity::{EventId, IndividualId};
/// fn needs_individual(_: IndividualId) {}
/// let event = EventId::parse("same", "execution").unwrap();
/// needs_individual(event);
/// ```
pub struct TypedIdentityInvariant;
