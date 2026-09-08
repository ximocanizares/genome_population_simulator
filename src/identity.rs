//! Identidades tipadas, etiquetas humanas y derivación automática versionada.

use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::diagnostics::Diagnostic;
use crate::randomness::{DERIVATION_VERSION, RootSeed};

pub const ID_CONTEXT: &str = "genome_population_simulator/v1/id";
pub const ID_CATALOG_VERSION: &str = "3";
pub const F5_ID_CATALOG: &str = "id-catalog-v4";
pub const TECHNICAL_ID_DERIVATION_VERSION: u32 = 2;

/// El texto de un identificador vive en una asignación compartida e inmutable.
///
/// Un identificador no se modifica nunca después de construirse, y F4 clona el
/// mismo valor millones de veces: la biblioteca de una pareja, el individuo de
/// cada clave limpia o la copia que el registro de colisiones conserva. Con
/// `Arc<str>` cada clon es un contador, no una asignación nueva, y el registro
/// comparte el texto que ya sostiene la clave en lugar de duplicarlo. El orden,
/// la igualdad, el hash y el `Debug` siguen siendo los del texto, de modo que
/// ningún sorteo ni artefacto cambia.
macro_rules! typed_id {
    ($name:ident, $category:literal) => {
        #[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
        pub struct $name(Arc<str>);

        impl $name {
            pub fn parse(value: impl Into<String>, scope: &str) -> Result<Self, Diagnostic> {
                let value = value.into();
                validate_provided_id(&value, $category, scope)?;
                Ok(Self(Arc::from(value)))
            }

            /// Acepta el texto recién derivado o el que ya vive compartido,
            /// de modo que el registro de colisiones y el identificador
            /// puedan sostener una única asignación.
            fn from_automatic(value: impl Into<Arc<str>>) -> Self {
                Self(value.into())
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
typed_id!(VariantId, "variant");
typed_id!(BatchId, "batch");
typed_id!(LibraryId, "library");
typed_id!(ReadGroupId, "read_group");
typed_id!(TemplateId, "template");
typed_id!(PreHopInstanceId, "pre_hop_instance");
typed_id!(TechnicalPairId, "technical_pair");
// F5 produces these identities. They deliberately are different domain types
// from supplied `BlockId` and from each other, even though all serialize as a
// compact ASCII identifier.
typed_id!(F5BlockId, "f5_block");
typed_id!(F5SiteId, "f5_site");
typed_id!(F5HaplotypeId, "f5_haplotype");
typed_id!(FounderCopyId, "f5_founder_copy");
typed_id!(
    FounderHaplotypeAssignmentId,
    "f5_block_haplotype_assignment"
);

/// Construye el ID de plantilla de F3 a partir de una biblioteca ya única en
/// el evento y de su ordinal lógico canónico en decimal.
///
/// Esta es la excepción acotada de DTC-039: no usa BLAKE3 ni requiere registrar
/// IDs previos, porque la pareja `(library_id, template_ordinal)` es única por
/// construcción. Los componentes ya cumplen la sintaxis de identificador y la
/// longitud máxima queda acotada a 55 caracteres (`TPL_` + 30 + `_` + 20).
#[must_use]
pub(crate) fn compose_template_id(library_id: &LibraryId, template_ordinal: u64) -> TemplateId {
    TemplateId::from_automatic(format!("TPL_{}_{template_ordinal}", library_id.as_str()))
}

macro_rules! provided_typed_id {
    ($name:ident, $category:literal) => {
        #[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
        pub struct $name(Arc<str>);

        impl $name {
            pub fn parse(value: impl Into<String>, scope: &str) -> Result<Self, Diagnostic> {
                let value = value.into();
                validate_provided_id(&value, $category, scope)?;
                Ok(Self(Arc::from(value)))
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
provided_typed_id!(FounderKey, "founder_key");
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
    Batch,
    Library,
    ReadGroup,
}

impl IdClass {
    const fn name(self) -> &'static str {
        match self {
            Self::Execution => "execution",
            Self::Reference => "reference",
            Self::Individual => "individual",
            Self::Event => "event",
            Self::Batch => "batch",
            Self::Library => "library",
            Self::ReadGroup => "read_group",
        }
    }

    const fn prefix(self) -> &'static str {
        match self {
            Self::Execution => "RUN_",
            Self::Reference => "REF_",
            Self::Individual => "IND_",
            Self::Event => "EVT_",
            Self::Batch => "BAT_",
            Self::Library => "LIB_",
            Self::ReadGroup => "RG_",
        }
    }

    const fn cause_names(self) -> &'static [&'static str] {
        match self {
            Self::Execution => &["config_sha256"],
            Self::Reference => &["source_kind", "source_digest"],
            Self::Individual => &["founder_key"],
            Self::Event => &["event_key"],
            Self::Batch => &[],
            Self::Library => &["individual_id"],
            Self::ReadGroup => &["batch_id", "library_id"],
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

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum TechnicalIdClass {
    PreHopInstance,
    TechnicalPair,
}

impl TechnicalIdClass {
    const fn name(self) -> &'static str {
        match self {
            Self::PreHopInstance => "PreHopInstanceId",
            Self::TechnicalPair => "TechnicalPairId",
        }
    }

    const fn prefix(self) -> &'static str {
        match self {
            Self::PreHopInstance => "PH_",
            Self::TechnicalPair => "TP_",
        }
    }
}

/// Registro de colisiones de las identidades técnicas transitorias de F4.
///
/// Las dos clases usan la derivación v2 de `id-catalog-v3`; no comparten el
/// material ni el registro de las clases heredadas v1.
///
/// El registro no conserva el texto de un identificador derivado, sino sus
/// 128 bits.
///
/// Un ID automático es el prefijo de su clase seguido del Base32 de los
/// primeros 128 bits del BLAKE3 (DTC-039), de modo que esos 16 bytes y su
/// texto se determinan mutuamente: retener los bytes detecta exactamente las
/// mismas colisiones y ahorra el texto de cada pareja del evento, que en F4.5
/// nadie más sostiene.
///
/// Un ID aportado no tiene por qué seguir esa forma, y DTC-039 exige detectar
/// también sus colisiones, así que se conserva su texto. La clasificación es
/// canónica y no ambigua: todo texto con la forma derivada va a `derived`,
/// venga de una derivación o de una configuración, y solo lo que no puede
/// proceder de una derivación va a `provided`. Por eso un aportado nunca puede
/// esconder la colisión con un derivado.
#[derive(Debug, Default)]
pub struct TechnicalIdentityRegistry {
    /// `(clase, 128 bits)` ocupa 17 bytes sin relleno, frente a los 16 de un
    /// puntero compartido más su asignación de texto.
    derived: BTreeSet<(TechnicalIdClass, [u8; 16])>,
    provided: BTreeSet<(TechnicalIdClass, Arc<str>)>,
}

impl TechnicalIdentityRegistry {
    /// Incorpora una instancia ya emitida por otra frontera del mismo evento.
    /// La causa original no hace falta para detectar una colisión: el ID es el
    /// único valor que puede aparecer después en FASTQ.
    pub fn register_existing_pre_hop_instance(
        &mut self,
        pre_hop_instance_id: &PreHopInstanceId,
    ) -> Result<(), Diagnostic> {
        self.register_text(
            TechnicalIdClass::PreHopInstance,
            &pre_hop_instance_id.0,
            "existing",
        )
    }

    /// Deriva una instancia propia antes de *index hopping*.
    pub fn derive_own_pre_hop_instance(
        &mut self,
        root: RootSeed,
        event_id: &EventId,
        library_id: &LibraryId,
        template_id: &TemplateId,
    ) -> Result<PreHopInstanceId, Diagnostic> {
        let (text, digest) = derive_technical_id(
            root,
            TechnicalIdClass::PreHopInstance,
            &[
                ("class", TechnicalIdClass::PreHopInstance.name()),
                ("event_id", event_id.as_str()),
                ("library_id", library_id.as_str()),
                ("origin", "own"),
                ("template_id", template_id.as_str()),
            ],
        )?;
        self.register_derived(TechnicalIdClass::PreHopInstance, digest, "own")?;
        Ok(PreHopInstanceId::from_automatic(text))
    }

    /// Deriva la copia de una instancia fuente en una biblioteca receptora.
    pub fn derive_contamination_pre_hop_instance(
        &mut self,
        root: RootSeed,
        event_id: &EventId,
        recipient_library_id: &LibraryId,
        source_pre_hop_instance_id: &PreHopInstanceId,
    ) -> Result<PreHopInstanceId, Diagnostic> {
        let (text, digest) = derive_technical_id(
            root,
            TechnicalIdClass::PreHopInstance,
            &[
                ("class", TechnicalIdClass::PreHopInstance.name()),
                ("event_id", event_id.as_str()),
                ("library_id", recipient_library_id.as_str()),
                ("origin", "contamination"),
                (
                    "source_pre_hop_instance_id",
                    source_pre_hop_instance_id.as_str(),
                ),
            ],
        )?;
        self.register_derived(TechnicalIdClass::PreHopInstance, digest, "contamination")?;
        Ok(PreHopInstanceId::from_automatic(text))
    }

    /// Deriva la identidad de una pareja después del demultiplexado.
    pub fn derive_technical_pair(
        &mut self,
        root: RootSeed,
        event_id: &EventId,
        final_library_id: &LibraryId,
        pre_hop_instance_id: &PreHopInstanceId,
    ) -> Result<TechnicalPairId, Diagnostic> {
        let (text, digest) = derive_technical_id(
            root,
            TechnicalIdClass::TechnicalPair,
            &[
                ("class", TechnicalIdClass::TechnicalPair.name()),
                ("event_id", event_id.as_str()),
                ("final_library_id", final_library_id.as_str()),
                ("pre_hop_instance_id", pre_hop_instance_id.as_str()),
            ],
        )?;
        self.register_derived(TechnicalIdClass::TechnicalPair, digest, "final")?;
        Ok(TechnicalPairId::from_automatic(text))
    }

    /// Anota un identificador recién derivado por sus 128 bits.
    fn register_derived(
        &mut self,
        class: TechnicalIdClass,
        digest: [u8; 16],
        causes: &str,
    ) -> Result<(), Diagnostic> {
        if self.derived.insert((class, digest)) {
            return Ok(());
        }
        Err(collision(class, &technical_id_text(class, &digest), causes))
    }

    /// Anota un identificador del que solo se conoce su texto. Si tiene la
    /// forma que produce una derivación, se anota entre los derivados —de otro
    /// modo la misma identidad viviría en dos conjuntos y su colisión pasaría
    /// inadvertida—; si no, se conserva el texto, que es lo único exacto.
    fn register_text(
        &mut self,
        class: TechnicalIdClass,
        identifier: &Arc<str>,
        causes: &str,
    ) -> Result<(), Diagnostic> {
        if let Some(digest) = decode_derived_id(class, identifier) {
            return self.register_derived(class, digest, causes);
        }
        if self.provided.insert((class, Arc::clone(identifier))) {
            return Ok(());
        }
        Err(collision(class, identifier, causes))
    }
}

fn collision(class: TechnicalIdClass, identifier: &str, causes: &str) -> Diagnostic {
    Diagnostic::domain(
        "GPS-DOMAIN-002",
        "La derivación de una identidad produjo una colisión.",
    )
    .with("identifier", identifier)
    .with("category", class.name())
    .with("causes", causes)
}

#[cfg(test)]
mod shared_text_tests {
    use super::{
        EventId, LibraryId, PreHopInstanceId, TechnicalIdClass, TechnicalIdentityRegistry,
        compose_template_id, decode_derived_id, technical_id_text,
    };
    use crate::randomness::RootSeed;
    use std::sync::Arc;

    /// Clonar una identidad no vuelve a asignar su texto. F4 depende de esto:
    /// cada clave limpia clona el individuo de su biblioteca.
    #[test]
    fn a_cloned_identity_shares_its_text() {
        let original = PreHopInstanceId::parse("PH_SHARED", "test").expect("valid identifier");
        let copy = original.clone();

        assert!(Arc::ptr_eq(&original.0, &copy.0));
        assert_eq!(original, copy);
        assert_eq!(copy.as_str(), "PH_SHARED");
    }

    /// Un identificador aportado que no puede proceder de una derivación
    /// conserva su texto, porque es lo único que lo identifica con exactitud.
    #[test]
    fn a_provided_identifier_keeps_its_text() {
        let existing = PreHopInstanceId::parse("PH_EXISTING", "test").expect("valid identifier");
        let mut registry = TechnicalIdentityRegistry::default();

        registry
            .register_existing_pre_hop_instance(&existing)
            .expect("a first registration is accepted");

        assert!(registry.derived.is_empty());
        let stored = registry
            .provided
            .iter()
            .map(|(_, text)| text)
            .find(|text| text.as_ref() == "PH_EXISTING")
            .expect("the instance is registered");
        assert!(Arc::ptr_eq(stored, &existing.0));

        let error = registry
            .register_existing_pre_hop_instance(&existing)
            .expect_err("a repeated identifier collides");
        assert_eq!(error.code(), "GPS-DOMAIN-002");
    }

    /// Una instancia derivada se anota por sus 128 bits, sin conservar su
    /// texto, y sigue detectando una segunda anotación del mismo identificador.
    #[test]
    fn a_derived_instance_is_registered_by_its_digest_and_still_detects_a_collision() {
        let root = RootSeed::from_bytes([0x11; 32]);
        let event_id = EventId::parse("EVT_shared", "test").expect("valid identifier");
        let library_id = LibraryId::parse("LIB_shared", "test").expect("valid identifier");
        let template_id = compose_template_id(&library_id, 0);
        let mut registry = TechnicalIdentityRegistry::default();

        let derived = registry
            .derive_own_pre_hop_instance(root, &event_id, &library_id, &template_id)
            .expect("the first derivation is registered");

        // El texto no se retiene: lo que queda son sus bits, y reconstruirlos
        // devuelve exactamente el mismo identificador.
        assert!(registry.provided.is_empty());
        let (_, digest) = registry
            .derived
            .iter()
            .find(|(class, _)| *class == TechnicalIdClass::PreHopInstance)
            .expect("the derived instance is registered");
        assert_eq!(
            technical_id_text(TechnicalIdClass::PreHopInstance, digest),
            derived.as_str()
        );

        let error = registry
            .register_existing_pre_hop_instance(&derived)
            .expect_err("a repeated identifier collides");
        assert_eq!(error.code(), "GPS-DOMAIN-002");
        assert_eq!(
            error.context().get("identifier").map(String::as_str),
            Some(derived.as_str())
        );
    }

    /// El riesgo de separar los dos conjuntos sería que la misma identidad
    /// viviera en ambos y su colisión pasara inadvertida. No ocurre: un texto
    /// con la forma derivada se anota entre los derivados aunque llegue como
    /// identificador aportado.
    #[test]
    fn a_provided_identifier_with_the_derived_shape_collides_with_the_derivation() {
        let root = RootSeed::from_bytes([0x22; 32]);
        let event_id = EventId::parse("EVT_shape", "test").expect("valid identifier");
        let library_id = LibraryId::parse("LIB_shape", "test").expect("valid identifier");
        let template_id = compose_template_id(&library_id, 0);

        let mut deriving = TechnicalIdentityRegistry::default();
        let derived = deriving
            .derive_own_pre_hop_instance(root, &event_id, &library_id, &template_id)
            .expect("the derivation is registered");

        // El mismo texto, esta vez recorriendo la vía de un ID aportado.
        let reparsed = PreHopInstanceId::parse(derived.as_str(), "test").expect("valid identifier");
        let mut registry = TechnicalIdentityRegistry::default();
        registry
            .register_existing_pre_hop_instance(&reparsed)
            .expect("a first registration is accepted");
        assert!(registry.provided.is_empty(), "it is stored as derived");

        let error = registry
            .derive_own_pre_hop_instance(root, &event_id, &library_id, &template_id)
            .expect_err("deriving the same identity collides");
        assert_eq!(error.code(), "GPS-DOMAIN-002");
    }

    /// Solo el texto que una derivación pudo producir se reconoce como tal.
    #[test]
    fn only_a_canonical_text_decodes_to_its_digest() {
        let digest = [0x5a_u8; 16];
        let text = technical_id_text(TechnicalIdClass::TechnicalPair, &digest);
        assert_eq!(
            decode_derived_id(TechnicalIdClass::TechnicalPair, &text),
            Some(digest)
        );

        // La clase forma parte de la identidad: el mismo texto no pertenece a
        // la otra clase.
        assert_eq!(
            decode_derived_id(TechnicalIdClass::PreHopInstance, &text),
            None
        );
        for rejected in [
            "TP_SHORT",
            "TP_AAAAAAAAAAAAAAAAAAAAAAAAAA1",
            "TP_AAAAAAAAAAAAAAAAAAAAAAAA1A",
            "PH_EXISTING",
            "",
        ] {
            assert_eq!(
                decode_derived_id(TechnicalIdClass::TechnicalPair, rejected),
                None,
                "{rejected} must not pass as a derived identifier"
            );
        }

        // Los dos bits de relleno a uno no los produce ninguna derivación: el
        // último símbolo canónico de estos bits es `A`, no `B`.
        let canonical = technical_id_text(TechnicalIdClass::TechnicalPair, &[0_u8; 16]);
        assert!(canonical.ends_with('A'));
        let padded = format!("{}B", &canonical[..canonical.len() - 1]);
        assert_eq!(
            decode_derived_id(TechnicalIdClass::TechnicalPair, &padded),
            None
        );
    }
}

/// Deriva un identificador técnico y entrega su texto junto a los 128 bits de
/// los que procede, para que el registro de colisiones pueda retener solo estos.
fn derive_technical_id(
    root: RootSeed,
    class: TechnicalIdClass,
    components: &[(&str, &str)],
) -> Result<(String, [u8; 16]), Diagnostic> {
    let component_count =
        u32::try_from(components.len()).map_err(|_| identity_material_error("component_count"))?;
    let mut material = Vec::new();
    material.extend_from_slice(root.as_bytes());
    material.extend_from_slice(&TECHNICAL_ID_DERIVATION_VERSION.to_be_bytes());
    material.extend_from_slice(&component_count.to_be_bytes());
    for (name, value) in components {
        append_component(&mut material, name, value)?;
    }
    let digest = blake3::derive_key(ID_CONTEXT, &material);
    let mut truncated = [0_u8; 16];
    truncated.copy_from_slice(&digest[..16]);
    Ok((technical_id_text(class, &truncated), truncated))
}

/// Compone el texto de un identificador técnico a partir de sus 128 bits.
fn technical_id_text(class: TechnicalIdClass, digest: &[u8; 16]) -> String {
    format!("{}{}", class.prefix(), encode_base32(digest))
}

/// Número de símbolos Base32 que codifican los 128 bits de un ID automático.
const DERIVED_SYMBOL_COUNT: usize = 26;

/// Instancia `PH_` retenida por sus símbolos, sin asignación propia.
///
/// El inventario de parejas limpias conserva una instancia por pareja de todo
/// el evento, y cada texto costaba una asignación. Como el prefijo es constante
/// para la clase, guardar solo los 26 símbolos basta para reconstruir el
/// identificador, y compararlos equivale a comparar el texto: el orden del
/// inventario —que gobierna qué parejas sustituye la contaminación— no se
/// mueve.
///
/// Solo admite la forma que produce una derivación. En una ejecución no hay
/// otra: una instancia `PH_` nunca se declara en configuración, siempre se
/// deriva.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct CompactPreHopInstanceId([u8; DERIVED_SYMBOL_COUNT]);

impl CompactPreHopInstanceId {
    /// Compacta una instancia derivada; rechaza cualquier texto que una
    /// derivación no haya podido producir.
    pub fn new(instance: &PreHopInstanceId) -> Result<Self, Diagnostic> {
        let class = TechnicalIdClass::PreHopInstance;
        let text = instance.as_str();
        if decode_derived_id(class, text).is_none() {
            return Err(Diagnostic::domain(
                "GPS-DOMAIN-001",
                "Un identificador técnico no tiene la forma que produce su derivación.",
            )
            .with("identifier", text)
            .with("category", class.name()));
        }
        let mut symbols = [0_u8; DERIVED_SYMBOL_COUNT];
        symbols.copy_from_slice(&text.as_bytes()[class.prefix().len()..]);
        Ok(Self(symbols))
    }

    /// Reconstruye el identificador completo.
    #[must_use]
    pub fn to_instance_id(self) -> PreHopInstanceId {
        let class = TechnicalIdClass::PreHopInstance;
        let mut text = String::with_capacity(class.prefix().len() + DERIVED_SYMBOL_COUNT);
        text.push_str(class.prefix());
        // Los símbolos provienen de un texto validado, así que son ASCII.
        text.push_str(&String::from_utf8_lossy(&self.0));
        PreHopInstanceId::from_automatic(text)
    }
}

/// Devuelve los 128 bits de un texto que solo una derivación pudo producir.
///
/// Exige el prefijo de la clase y que recodificar los bytes reproduzca el texto
/// exacto. Esa comprobación descarta cualquier variante no canónica —longitud
/// distinta, un carácter fuera del alfabeto o los dos bits de relleno a uno—,
/// de modo que un texto aceptado aquí procede necesariamente de los bytes que
/// se devuelven, y ninguno rechazado puede coincidir con un identificador
/// derivado.
fn decode_derived_id(class: TechnicalIdClass, text: &str) -> Option<[u8; 16]> {
    const ALPHABET: &[u8; 32] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
    let encoded = text.strip_prefix(class.prefix())?.as_bytes();
    if encoded.len() != DERIVED_SYMBOL_COUNT {
        return None;
    }
    let mut digest = [0_u8; 16];
    let mut buffer = 0_u32;
    let mut bits = 0_u8;
    let mut written = 0_usize;
    for symbol in encoded {
        let value = ALPHABET.iter().position(|candidate| candidate == symbol)?;
        buffer = (buffer << 5) | u32::try_from(value).ok()?;
        bits += 5;
        if bits >= 8 {
            bits -= 8;
            if written == digest.len() {
                return None;
            }
            digest[written] = u8::try_from((buffer >> bits) & 0xff).ok()?;
            written += 1;
        }
    }
    if written != digest.len() || encode_base32(&digest) != text[class.prefix().len()..] {
        return None;
    }
    Some(digest)
}

/// Registro de las cinco clases añadidas por `id-catalog-v4`.
///
/// Es independiente de las clases heredadas: publicar F5 no modifica sus
/// materiales, vectores ni ámbitos de colisión.
#[derive(Debug, Default)]
pub struct F5IdentityRegistry {
    registered: BTreeSet<(F5IdClass, String)>,
}

#[derive(Clone, Copy, Debug)]
pub struct F5IdentityContext<'a> {
    root: RootSeed,
    execution_id: &'a ExecutionId,
    scenario_id: &'a ScenarioId,
    reference_id: &'a ReferenceId,
}

impl<'a> F5IdentityContext<'a> {
    #[must_use]
    pub const fn new(
        root: RootSeed,
        execution_id: &'a ExecutionId,
        scenario_id: &'a ScenarioId,
        reference_id: &'a ReferenceId,
    ) -> Self {
        Self {
            root,
            execution_id,
            scenario_id,
            reference_id,
        }
    }

    #[must_use]
    pub const fn root_seed(self) -> RootSeed {
        self.root
    }

    #[must_use]
    pub const fn scenario_id(self) -> &'a ScenarioId {
        self.scenario_id
    }
}

#[derive(Clone, Copy, Debug)]
pub struct F5PopulationIdentityContext<'a> {
    base: F5IdentityContext<'a>,
    population_id: &'a PopulationId,
}

impl<'a> F5PopulationIdentityContext<'a> {
    #[must_use]
    pub const fn new(base: F5IdentityContext<'a>, population_id: &'a PopulationId) -> Self {
        Self {
            base,
            population_id,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct F5FounderCopyScope<'a> {
    individual_id: &'a IndividualId,
    sequence_id: &'a SequenceId,
    phase_id: &'a PhaseId,
}

impl<'a> F5FounderCopyScope<'a> {
    #[must_use]
    pub const fn new(
        individual_id: &'a IndividualId,
        sequence_id: &'a SequenceId,
        phase_id: &'a PhaseId,
    ) -> Self {
        Self {
            individual_id,
            sequence_id,
            phase_id,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct F5BlockHaplotypeAssignmentScope<'a> {
    block_index: u64,
    founder_copy_id: &'a FounderCopyId,
    haplotype_id: &'a F5HaplotypeId,
}

impl<'a> F5BlockHaplotypeAssignmentScope<'a> {
    #[must_use]
    pub const fn new(
        block_index: u64,
        founder_copy_id: &'a FounderCopyId,
        haplotype_id: &'a F5HaplotypeId,
    ) -> Self {
        Self {
            block_index,
            founder_copy_id,
            haplotype_id,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum F5IdClass {
    Block,
    Site,
    Haplotype,
    FounderCopy,
    BlockHaplotypeAssignment,
}

impl F5IdClass {
    const fn name(self) -> &'static str {
        match self {
            Self::Block => "f5_block",
            Self::Site => "f5_site",
            Self::Haplotype => "f5_haplotype",
            Self::FounderCopy => "f5_founder_copy",
            Self::BlockHaplotypeAssignment => "f5_block_haplotype_assignment",
        }
    }
    const fn prefix(self) -> &'static str {
        match self {
            Self::Block => "BLK_",
            Self::Site => "SITE_",
            Self::Haplotype => "HAP_",
            Self::FounderCopy => "FCP_",
            Self::BlockHaplotypeAssignment => "FHA_",
        }
    }
}

impl F5IdentityRegistry {
    pub fn derive_block(
        &mut self,
        context: F5IdentityContext<'_>,
        block_index: u64,
    ) -> Result<F5BlockId, Diagnostic> {
        let block_index = block_index.to_string();
        self.derive(
            context.root,
            F5IdClass::Block,
            &[
                ("execution_id", context.execution_id.as_str()),
                ("scenario_id", context.scenario_id.as_str()),
                ("reference_id", context.reference_id.as_str()),
                ("block_index_u64", &block_index),
            ],
        )
        .map(F5BlockId::from_automatic)
    }

    pub fn derive_site(
        &mut self,
        context: F5IdentityContext<'_>,
        block_index: u64,
        site_index: u64,
    ) -> Result<F5SiteId, Diagnostic> {
        let block_index = block_index.to_string();
        let site_index = site_index.to_string();
        self.derive(
            context.root,
            F5IdClass::Site,
            &[
                ("execution_id", context.execution_id.as_str()),
                ("scenario_id", context.scenario_id.as_str()),
                ("reference_id", context.reference_id.as_str()),
                ("block_index_u64", &block_index),
                ("site_index_u64", &site_index),
            ],
        )
        .map(F5SiteId::from_automatic)
    }

    pub fn derive_haplotype(
        &mut self,
        context: F5PopulationIdentityContext<'_>,
        block_index: u64,
        haplotype_rank: u64,
        alleles: &str,
    ) -> Result<F5HaplotypeId, Diagnostic> {
        let block_index = block_index.to_string();
        let haplotype_rank = haplotype_rank.to_string();
        self.derive(
            context.base.root,
            F5IdClass::Haplotype,
            &[
                ("execution_id", context.base.execution_id.as_str()),
                ("scenario_id", context.base.scenario_id.as_str()),
                ("reference_id", context.base.reference_id.as_str()),
                ("population_id", context.population_id.as_str()),
                ("block_index_u64", &block_index),
                ("haplotype_rank_u64", &haplotype_rank),
                ("alleles", alleles),
            ],
        )
        .map(F5HaplotypeId::from_automatic)
    }

    pub fn derive_founder_copy(
        &mut self,
        context: F5PopulationIdentityContext<'_>,
        scope: F5FounderCopyScope<'_>,
    ) -> Result<FounderCopyId, Diagnostic> {
        self.derive(
            context.base.root,
            F5IdClass::FounderCopy,
            &[
                ("execution_id", context.base.execution_id.as_str()),
                ("scenario_id", context.base.scenario_id.as_str()),
                ("reference_id", context.base.reference_id.as_str()),
                ("population_id", context.population_id.as_str()),
                ("individual_id", scope.individual_id.as_str()),
                ("sequence_id", scope.sequence_id.as_str()),
                ("phase_id", scope.phase_id.as_str()),
            ],
        )
        .map(FounderCopyId::from_automatic)
    }

    pub fn derive_block_haplotype_assignment(
        &mut self,
        context: F5PopulationIdentityContext<'_>,
        scope: F5BlockHaplotypeAssignmentScope<'_>,
    ) -> Result<FounderHaplotypeAssignmentId, Diagnostic> {
        let block_index = scope.block_index.to_string();
        self.derive(
            context.base.root,
            F5IdClass::BlockHaplotypeAssignment,
            &[
                ("execution_id", context.base.execution_id.as_str()),
                ("scenario_id", context.base.scenario_id.as_str()),
                ("reference_id", context.base.reference_id.as_str()),
                ("population_id", context.population_id.as_str()),
                ("block_index_u64", &block_index),
                ("founder_copy_id", scope.founder_copy_id.as_str()),
                ("haplotype_id", scope.haplotype_id.as_str()),
            ],
        )
        .map(FounderHaplotypeAssignmentId::from_automatic)
    }

    fn derive(
        &mut self,
        root: RootSeed,
        class: F5IdClass,
        causes: &[(&str, &str)],
    ) -> Result<String, Diagnostic> {
        self.derive_with_digest(root, class, causes, |material| {
            blake3::derive_key(ID_CONTEXT, material)
        })
    }

    fn derive_with_digest(
        &mut self,
        root: RootSeed,
        class: F5IdClass,
        causes: &[(&str, &str)],
        digest: impl FnOnce(&[u8]) -> [u8; 32],
    ) -> Result<String, Diagnostic> {
        let mut material = Vec::new();
        material.extend_from_slice(root.as_bytes());
        material.extend_from_slice(&DERIVATION_VERSION.to_be_bytes());
        material.extend_from_slice(
            &u32::try_from(causes.len() + 1)
                .map_err(|_| identity_material_error("component_count"))?
                .to_be_bytes(),
        );
        append_component(&mut material, "class", class.name())?;
        for (name, value) in causes {
            append_component(&mut material, name, value)?;
        }
        let digest = digest(&material);
        let identifier = format!("{}{}", class.prefix(), encode_base32(&digest[..16]));
        if self.registered.insert((class, identifier.clone())) {
            return Ok(identifier);
        }
        Err(Diagnostic::domain(
            "GPS-DOMAIN-002",
            "La derivación de una identidad produjo una colisión.",
        )
        .with("identifier", identifier)
        .with("category", class.name())
        .with("causes", "automatic"))
    }
}

#[cfg(test)]
mod f5_collision_tests {
    use super::{F5IdClass, F5IdentityRegistry};
    use crate::randomness::RootSeed;

    #[test]
    fn distinct_f5_causes_reject_a_forced_digest_collision() {
        let root = RootSeed::from_bytes([0x5a; 32]);
        let forced_digest = [0xa5; 32];
        let mut registry = F5IdentityRegistry::default();
        let first_causes = [
            ("execution_id", "RUN_collision"),
            ("scenario_id", "scenario-f5"),
            ("reference_id", "reference-main"),
            ("population_id", "population-A"),
            ("individual_id", "individual-A"),
            ("sequence_id", "sequence-1"),
            ("phase_id", "phase-0"),
        ];
        let second_causes = [
            ("execution_id", "RUN_collision"),
            ("scenario_id", "scenario-f5"),
            ("reference_id", "reference-main"),
            ("population_id", "population-A"),
            ("individual_id", "individual-A"),
            ("sequence_id", "sequence-1"),
            ("phase_id", "phase-1"),
        ];

        let first = registry
            .derive_with_digest(root, F5IdClass::FounderCopy, &first_causes, |_| {
                forced_digest
            })
            .expect("first forced digest is registered");
        let error = registry
            .derive_with_digest(root, F5IdClass::FounderCopy, &second_causes, |_| {
                forced_digest
            })
            .expect_err("different phase with the same forced digest must collide");

        assert_eq!(error.code(), "GPS-DOMAIN-002");
        assert_eq!(error.context().get("identifier"), Some(&first));
        assert_eq!(
            error.context().get("category").map(String::as_str),
            Some("f5_founder_copy")
        );
        assert_ne!(first_causes, second_causes);
    }
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
    Batch,
    Library,
    ReadGroup,
}

impl From<IdClass> for IdClassKey {
    fn from(value: IdClass) -> Self {
        match value {
            IdClass::Execution => Self::Execution,
            IdClass::Reference => Self::Reference,
            IdClass::Individual => Self::Individual,
            IdClass::Event => Self::Event,
            IdClass::Batch => Self::Batch,
            IdClass::Library => Self::Library,
            IdClass::ReadGroup => Self::ReadGroup,
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

    pub fn derive_batch(
        &mut self,
        root: RootSeed,
        request: &AutomaticIdRequest,
    ) -> Result<BatchId, Diagnostic> {
        self.require_class(request, IdClass::Batch)?;
        let text = request.derive_text(root)?;
        self.register_text(request.class, &request.scope, &text, "automatic")?;
        Ok(BatchId::from_automatic(text))
    }

    pub fn derive_library(
        &mut self,
        root: RootSeed,
        request: &AutomaticIdRequest,
    ) -> Result<LibraryId, Diagnostic> {
        self.require_class(request, IdClass::Library)?;
        let text = request.derive_text(root)?;
        self.register_text(request.class, &request.scope, &text, "automatic")?;
        Ok(LibraryId::from_automatic(text))
    }

    pub fn derive_read_group(
        &mut self,
        root: RootSeed,
        request: &AutomaticIdRequest,
    ) -> Result<ReadGroupId, Diagnostic> {
        self.require_class(request, IdClass::ReadGroup)?;
        let text = request.derive_text(root)?;
        self.register_text(request.class, &request.scope, &text, "automatic")?;
        Ok(ReadGroupId::from_automatic(text))
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

/// Deriva el ID de una variante de matriz a partir de su locus canónico.
///
/// A diferencia de las entidades sorteadas, este ID no incorpora semilla: debe
/// ser idéntico al reutilizar una misma matriz bajo otra ejecución. Las causas
/// están completas y ordenadas para que no dependa de población, fila, ruta ni
/// orden físico.
pub fn derive_matrix_variant_id(
    reference_id: &ReferenceId,
    sequence_id: &SequenceId,
    position_1based: u64,
    reference: u8,
    alternate: u8,
) -> Result<VariantId, Diagnostic> {
    let position = position_1based.to_string();
    let reference = char::from(reference).to_string();
    let alternate = char::from(alternate).to_string();
    let components = [
        ("class", "variant"),
        ("reference_id", reference_id.as_str()),
        ("sequence_id", sequence_id.as_str()),
        ("position_1based", position.as_str()),
        ("ref", reference.as_str()),
        ("alt", alternate.as_str()),
    ];
    let mut material = Vec::new();
    material.extend_from_slice(
        &u32::try_from(components.len())
            .map_err(|_| identity_material_error("component_count"))?
            .to_be_bytes(),
    );
    for (name, value) in components {
        append_component(&mut material, name, value)?;
    }
    let digest = blake3::derive_key("genome_population_simulator/v1/variant-id", &material);
    Ok(VariantId::from_automatic(format!(
        "VAR_{}",
        encode_base32(&digest[..16])
    )))
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
