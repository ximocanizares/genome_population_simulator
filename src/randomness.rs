//! Semillas, flujos aleatorios locales y contratos de reproducibilidad.

use std::fmt;

use rand_chacha::ChaCha12Rng;
use rand_chacha::rand_core::{Rng, SeedableRng};

use crate::diagnostics::Diagnostic;

pub const DERIVATION_VERSION: u32 = 1;
pub const SUBSEED_CATALOG_VERSION: &str = "3";
pub const DERIVATION_CONTEXT: &str = "genome_population_simulator/v1/subseed";
pub const RANDOM_ALGORITHM: &str = "ChaCha12";
pub const REFERENCE_REDUCTION_ALGORITHM: &str = "reduced_reference_selection_v1";

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
                ("algorithm", REFERENCE_REDUCTION_ALGORITHM),
            ],
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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReproducibilityMetadata {
    root_seed: RootSeed,
    seed_origin: SeedOrigin,
    random_algorithm: &'static str,
    derivation_version: u32,
    catalog_version: &'static str,
}

impl ReproducibilityMetadata {
    #[must_use]
    pub const fn new(root_seed: RootSeed, seed_origin: SeedOrigin) -> Self {
        Self {
            root_seed,
            seed_origin,
            random_algorithm: RANDOM_ALGORITHM,
            derivation_version: DERIVATION_VERSION,
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
    pub const fn derivation_version(&self) -> u32 {
        self.derivation_version
    }

    #[must_use]
    pub const fn catalog_version(&self) -> &'static str {
        self.catalog_version
    }
}
