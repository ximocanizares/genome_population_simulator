//! Perfil técnico incorporado de F4.
//!
//! Este módulo es deliberadamente una frontera de datos: valida el único
//! recurso versionado aprobado por DTPD-106 y expone los hechos que necesitarán
//! la construcción de pares. No selecciona el perfil desde TOML ni materializa
//! calidades, errores, contaminación o FASTQ; esas responsabilidades pertenecen
//! a hitos posteriores.

use crate::guard_trace::rejected;
use std::fmt;

use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::pcr_amplification::{
    PCR_AMPLIFICATION_ALGORITHM, PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS,
    PCR_DUPLICATE_FRACTION_INVERSION_METHOD,
};

/// Identificador del único perfil incorporado inicialmente en F4.
pub const GPS_HISEQX_PCRFREE_PE150_V1: &str = "gps-hiseqx-pcrfree-pe150-v1";
/// Versión semántica del perfil, separada de la versión de su recurso JSON.
pub const GPS_HISEQX_PCRFREE_PE150_PROFILE_VERSION: u8 = 1;
/// Ruta reproducible del recurso dentro del árbol fuente.
pub const GPS_HISEQX_PCRFREE_PE150_RESOURCE: &str = "resources/f4/gps-hiseqx-pcrfree-pe150-v1.json";
/// SHA-256 normativo fijado por DTPD-106.
pub const GPS_HISEQX_PCRFREE_PE150_SHA256: &str =
    "3d8878ab920237dae03705276a26dbb0823941830d531374d35aa1c540f6deea";

/// Identificador del perfil **hermano** con preparación amplificada (F4.9.1).
///
/// # Por qué es un recurso aparte y no una preparación más dentro del vigente
///
/// El nombre `pcrfree` del recurso vigente es **literal** en tres esquemas
/// publicados y su SHA-256 está congelado en el código y en los documentos ya
/// emitidos. Meter las dos preparaciones en un solo recurso cambiaría esos
/// bytes y rompería los valores fijos de los esquemas. La decisión S-25 lo
/// resuelve con un hermano: el vigente no se toca ni en un byte.
pub const GPS_HISEQX_PCR_PE150_V1: &str = "gps-hiseqx-pcr-pe150-v1";
/// Versión semántica del perfil amplificado.
pub const GPS_HISEQX_PCR_PE150_PROFILE_VERSION: u8 = 1;
/// Ruta reproducible del recurso amplificado dentro del árbol fuente.
pub const GPS_HISEQX_PCR_PE150_RESOURCE: &str = "resources/f4/gps-hiseqx-pcr-pe150-v1.json";
/// SHA-256 normativo del recurso amplificado.
pub const GPS_HISEQX_PCR_PE150_SHA256: &str =
    "83b00fb8d084df8f6aef0464fd2b0260cfefe3f89fc523b135da6bbf6595a03c";

// ---------------------------------------------------------------------------
// Constantes de contrato de la amplificación (F4.9.2).
//
// **Una sola constante de Rust por cada una.** El recurso hermano, los dos
// esquemas y la interfaz repiten estos mismos números en forma de `maximum`, y
// un juez —`tests/pcr_contract.rs`— compara cada copia con su
// constante. La copia decimal silenciosa de un máximo en un solo esquema es lo
// que sobrevivió al cambio de contrato la vez anterior.
// ---------------------------------------------------------------------------

/// `C_MAX`: ciclos de PCR admitidos como máximo.
pub const PCR_MAXIMUM_CYCLES: u32 = 30;
/// `E_MAX`: ciclos explícitos admitidos como máximo.
pub const PCR_MAXIMUM_EARLY_CYCLES: u32 = 6;
/// `K_MAX`: tamaño máximo de familia admitido como máximo.
pub const PCR_MAXIMUM_FAMILY_SIZE: u32 = 256;
/// `A_MAX`: intentos de aceptación por composición admitidos como máximo.
pub const PCR_MAXIMUM_GC_ATTEMPTS: u32 = 16;
/// `Q_MAX`: clases de eficiencia admitidas como máximo.
pub const PCR_MAXIMUM_EFFICIENCY_CLASSES: u32 = 4;
/// `M_MAX`: parejas por molécula que la inversión numérica puede alcanzar.
pub const PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE: f64 = 64.0;

const RESOURCE_BYTES: &[u8] = include_bytes!("../resources/f4/gps-hiseqx-pcrfree-pe150-v1.json");
const PCR_RESOURCE_BYTES: &[u8] = include_bytes!("../resources/f4/gps-hiseqx-pcr-pe150-v1.json");

/// Preparación de biblioteca de un perfil incorporado.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LibraryPreparation {
    PcrFree,
    Pcr,
}

impl LibraryPreparation {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::PcrFree => "pcr_free",
            Self::Pcr => "pcr",
        }
    }

    /// Lee la preparación declarada por el autor del escenario.
    #[must_use]
    pub fn parse(value: &str) -> Option<Self> {
        match value {
            "pcr_free" => Some(Self::PcrFree),
            "pcr" => Some(Self::Pcr),
            _ => None,
        }
    }

    /// Identificador del perfil incorporado que corresponde a la preparación.
    #[must_use]
    pub const fn profile_id(self) -> &'static str {
        match self {
            Self::PcrFree => GPS_HISEQX_PCRFREE_PE150_V1,
            Self::Pcr => GPS_HISEQX_PCR_PE150_V1,
        }
    }
}

/// Disposición de lectura determinada por el recurso, no por la configuración.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReadLayout {
    PairedEnd,
}

/// Orientación del par respecto al inserto.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairOrientation {
    Fr,
}

/// Límites de una longitud común para ambos mates.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CommonReadLengthBounds {
    minimum_bases: u16,
    maximum_bases: u16,
}

impl CommonReadLengthBounds {
    #[must_use]
    pub const fn minimum_bases(self) -> u16 {
        self.minimum_bases
    }

    #[must_use]
    pub const fn maximum_bases(self) -> u16 {
        self.maximum_bases
    }

    #[must_use]
    pub const fn contains(self, length_bases: u16) -> bool {
        length_bases >= self.minimum_bases && length_bases <= self.maximum_bases
    }
}

/// Hechos normativos de la distribución de inserto que consumirá F4.1.3.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct InsertLengthDistribution {
    algorithm: &'static str,
    mean_bases: u16,
    standard_deviation_bases: u16,
    minimum_bases: u16,
    maximum_bases: u16,
}

/// Un tramo inclusivo de la curva de calidad por ciclo del perfil.
///
/// Los extremos son enteros del recurso DTPD-106. La interpolación y el
/// redondeo pertenecen al consumidor de calidad: este tipo solo impide que el
/// consumidor replique o invente una tabla distinta.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct QualitySegment {
    first_cycle: u16,
    last_cycle: u16,
    first_phred: u8,
    last_phred: u8,
}

impl QualitySegment {
    #[must_use]
    pub const fn first_cycle(self) -> u16 {
        self.first_cycle
    }

    #[must_use]
    pub const fn last_cycle(self) -> u16 {
        self.last_cycle
    }

    #[must_use]
    pub const fn first_phred(self) -> u8 {
        self.first_phred
    }

    #[must_use]
    pub const fn last_phred(self) -> u8 {
        self.last_phred
    }
}

/// Curva versionada de calidad de uno de los mates PE.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct QualityCurve {
    segments: &'static [QualitySegment],
}

impl QualityCurve {
    #[must_use]
    pub const fn segments(self) -> &'static [QualitySegment] {
        self.segments
    }
}

/// Mate seleccionado por el perfil técnico. Se mantiene separado del rol de
/// identidad para que el recurso no dependa de las identidades transitorias.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProfileMate {
    R1,
    R2,
}

impl InsertLengthDistribution {
    #[must_use]
    pub const fn algorithm(self) -> &'static str {
        self.algorithm
    }

    #[must_use]
    pub const fn mean_bases(self) -> u16 {
        self.mean_bases
    }

    #[must_use]
    pub const fn standard_deviation_bases(self) -> u16 {
        self.standard_deviation_bases
    }

    #[must_use]
    pub const fn minimum_bases(self) -> u16 {
        self.minimum_bases
    }

    #[must_use]
    pub const fn maximum_bases(self) -> u16 {
        self.maximum_bases
    }
}

/// Perfil F4 ya validado, sin estado de una ejecución concreta.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F4TechnicalProfile {
    id: &'static str,
    version: u8,
    resource_path: &'static str,
    resource_sha256: &'static str,
    preparation: LibraryPreparation,
    layout: ReadLayout,
    orientation: PairOrientation,
    common_read_length: CommonReadLengthBounds,
    insert_length: InsertLengthDistribution,
    quality_r1: QualityCurve,
    quality_r2: QualityCurve,
}

impl F4TechnicalProfile {
    #[must_use]
    pub const fn id(self) -> &'static str {
        self.id
    }

    #[must_use]
    pub const fn version(self) -> u8 {
        self.version
    }

    #[must_use]
    pub const fn resource_path(self) -> &'static str {
        self.resource_path
    }

    #[must_use]
    pub const fn resource_sha256(self) -> &'static str {
        self.resource_sha256
    }

    #[must_use]
    pub const fn preparation(self) -> LibraryPreparation {
        self.preparation
    }

    #[must_use]
    pub const fn layout(self) -> ReadLayout {
        self.layout
    }

    #[must_use]
    pub const fn orientation(self) -> PairOrientation {
        self.orientation
    }

    #[must_use]
    pub const fn common_read_length(self) -> CommonReadLengthBounds {
        self.common_read_length
    }

    #[must_use]
    pub const fn insert_length(self) -> InsertLengthDistribution {
        self.insert_length
    }

    /// Frecuencia predeterminada de quimeras del perfil incorporado. DTPD-106
    /// fija el neutro: una tasa distinta requiere configuración F4 explícita.
    #[must_use]
    pub const fn chimera_rate(self) -> f64 {
        0.0
    }

    #[must_use]
    pub const fn quality_curve(self, mate: ProfileMate) -> QualityCurve {
        match mate {
            ProfileMate::R1 => self.quality_r1,
            ProfileMate::R2 => self.quality_r2,
        }
    }
}

const QUALITY_R1_SEGMENTS: [QualitySegment; 3] = [
    QualitySegment {
        first_cycle: 1,
        last_cycle: 10,
        first_phred: 36,
        last_phred: 36,
    },
    QualitySegment {
        first_cycle: 11,
        last_cycle: 100,
        first_phred: 36,
        last_phred: 32,
    },
    QualitySegment {
        first_cycle: 101,
        last_cycle: 150,
        first_phred: 32,
        last_phred: 28,
    },
];

const QUALITY_R2_SEGMENTS: [QualitySegment; 3] = [
    QualitySegment {
        first_cycle: 1,
        last_cycle: 10,
        first_phred: 35,
        last_phred: 35,
    },
    QualitySegment {
        first_cycle: 11,
        last_cycle: 100,
        first_phred: 35,
        last_phred: 30,
    },
    QualitySegment {
        first_cycle: 101,
        last_cycle: 150,
        first_phred: 30,
        last_phred: 25,
    },
];

const INITIAL_PROFILE: F4TechnicalProfile = F4TechnicalProfile {
    id: GPS_HISEQX_PCRFREE_PE150_V1,
    version: GPS_HISEQX_PCRFREE_PE150_PROFILE_VERSION,
    resource_path: GPS_HISEQX_PCRFREE_PE150_RESOURCE,
    resource_sha256: GPS_HISEQX_PCRFREE_PE150_SHA256,
    preparation: LibraryPreparation::PcrFree,
    layout: ReadLayout::PairedEnd,
    orientation: PairOrientation::Fr,
    common_read_length: CommonReadLengthBounds {
        minimum_bases: 1,
        maximum_bases: 150,
    },
    insert_length: InsertLengthDistribution {
        algorithm: "normal_truncated_box_muller_v1",
        mean_bases: 350,
        standard_deviation_bases: 20,
        minimum_bases: 250,
        maximum_bases: 450,
    },
    quality_r1: QualityCurve {
        segments: &QUALITY_R1_SEGMENTS,
    },
    quality_r2: QualityCurve {
        segments: &QUALITY_R2_SEGMENTS,
    },
};

/// Perfil hermano con preparación amplificada.
///
/// **Comparte literalmente** las curvas de calidad y la distribución de
/// inserto del vigente: son las mismas constantes de Rust, no copias, de modo
/// que no pueden separarse por descuido.
const PCR_PROFILE: F4TechnicalProfile = F4TechnicalProfile {
    id: GPS_HISEQX_PCR_PE150_V1,
    version: GPS_HISEQX_PCR_PE150_PROFILE_VERSION,
    resource_path: GPS_HISEQX_PCR_PE150_RESOURCE,
    resource_sha256: GPS_HISEQX_PCR_PE150_SHA256,
    preparation: LibraryPreparation::Pcr,
    layout: ReadLayout::PairedEnd,
    orientation: PairOrientation::Fr,
    common_read_length: INITIAL_PROFILE.common_read_length,
    insert_length: INITIAL_PROFILE.insert_length,
    quality_r1: INITIAL_PROFILE.quality_r1,
    quality_r2: INITIAL_PROFILE.quality_r2,
};

/// Hechos publicables de un perfil incorporado, **sin recalcular su huella**.
///
/// # Por qué existe aparte de [`resolve_f4_technical_profile`]
///
/// Resolver un perfil vuelve a leer el recurso incorporado y a calcular su
/// SHA-256 entero. Eso está bien cuando hay que **comprobarlo**, que es al
/// admitir la configuración, y está mal cuando solo hay que **publicar** su
/// identificador y su ruta, que es lo que hacen el documento efectivo y el
/// manifiesto cada vez que se escriben. Esta función devuelve los literales
/// desde las constantes, sin tocar el recurso ni poder fallar.
#[must_use]
pub const fn f4_profile_facts(preparation: LibraryPreparation) -> F4ProfileFacts {
    match preparation {
        LibraryPreparation::PcrFree => F4ProfileFacts {
            id: GPS_HISEQX_PCRFREE_PE150_V1,
            resource_path: GPS_HISEQX_PCRFREE_PE150_RESOURCE,
            sha256: GPS_HISEQX_PCRFREE_PE150_SHA256,
            preparation: LibraryPreparation::PcrFree,
            insert_length: INITIAL_PROFILE.insert_length,
        },
        LibraryPreparation::Pcr => F4ProfileFacts {
            id: GPS_HISEQX_PCR_PE150_V1,
            resource_path: GPS_HISEQX_PCR_PE150_RESOURCE,
            sha256: GPS_HISEQX_PCR_PE150_SHA256,
            preparation: LibraryPreparation::Pcr,
            insert_length: PCR_PROFILE.insert_length,
        },
    }
}

/// Los mismos hechos, a partir del identificador publicado. `None` cuando el
/// identificador no es el de ningún perfil incorporado.
#[must_use]
pub fn f4_profile_facts_by_id(profile_id: &str) -> Option<F4ProfileFacts> {
    match profile_id {
        GPS_HISEQX_PCRFREE_PE150_V1 => Some(f4_profile_facts(LibraryPreparation::PcrFree)),
        GPS_HISEQX_PCR_PE150_V1 => Some(f4_profile_facts(LibraryPreparation::Pcr)),
        _ => None,
    }
}

/// La terna que identifica un perfil incorporado, siempre junta.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F4ProfileFacts {
    id: &'static str,
    resource_path: &'static str,
    sha256: &'static str,
    preparation: LibraryPreparation,
    insert_length: InsertLengthDistribution,
}

impl F4ProfileFacts {
    #[must_use]
    pub const fn id(self) -> &'static str {
        self.id
    }

    #[must_use]
    pub const fn resource_path(self) -> &'static str {
        self.resource_path
    }

    #[must_use]
    pub const fn sha256(self) -> &'static str {
        self.sha256
    }

    #[must_use]
    pub const fn preparation(self) -> LibraryPreparation {
        self.preparation
    }

    /// La distribución de inserto, tomada del **mismo** perfil que el
    /// validador del recurso comprueba: no es una copia con sus propios
    /// literales.
    #[must_use]
    pub const fn insert_length(self) -> InsertLengthDistribution {
        self.insert_length
    }
}

/// Fallo estable al resolver o comprobar un perfil incorporado.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum F4TechnicalProfileError {
    UnknownProfile { requested: String },
    ResourceDigestMismatch { observed: String },
    InvalidResource { detail: &'static str },
}

impl fmt::Display for F4TechnicalProfileError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownProfile { requested } => write!(
                formatter,
                "el perfil técnico F4 no está incorporado: {requested}"
            ),
            Self::ResourceDigestMismatch { observed } => write!(
                formatter,
                "el recurso técnico F4 no coincide con el SHA-256 normativo: {observed}"
            ),
            Self::InvalidResource { detail } => {
                write!(
                    formatter,
                    "el recurso técnico F4 no cumple DTPD-106: {detail}"
                )
            }
        }
    }
}

impl std::error::Error for F4TechnicalProfileError {}

/// Resuelve un perfil incorporado de F4 después de comprobar su recurso.
///
/// La frontera no interpreta aún un campo de configuración: permite que F4.1.2
/// y F4.1.3 reciban un perfil tipado sin hacer del identificador una cadena
/// implícita en sus algoritmos. Desde F4.9.1 hay **dos** perfiles hermanos, y
/// el identificador elige entre ellos.
pub fn resolve_f4_technical_profile(
    profile_id: &str,
) -> Result<F4TechnicalProfile, F4TechnicalProfileError> {
    match profile_id {
        GPS_HISEQX_PCRFREE_PE150_V1 => {
            verify_initial_profile_resource()?;
            Ok(INITIAL_PROFILE)
        }
        GPS_HISEQX_PCR_PE150_V1 => {
            verify_pcr_profile_resource()?;
            Ok(PCR_PROFILE)
        }
        _ => Err(rejected(F4TechnicalProfileError::UnknownProfile {
            requested: profile_id.to_owned(),
        })),
    }
}

/// Resuelve el perfil incorporado que corresponde a una preparación.
pub fn resolve_f4_technical_profile_for(
    preparation: LibraryPreparation,
) -> Result<F4TechnicalProfile, F4TechnicalProfileError> {
    resolve_f4_technical_profile(preparation.profile_id())
}

/// Comprueba el recurso literal y su SHA-256 antes de que un consumidor F4 lo
/// use. Se mantiene separada de la configuración hasta la transición atómica
/// de contratos de F4.1.3.
pub fn verify_initial_profile_resource() -> Result<(), F4TechnicalProfileError> {
    // Inalcanzable el rechazo de huella y el de JSON ilegible que siguen: se
    // comparan los bytes incrustados en el binario con una constante y el JSON se
    // lee de esos mismos bytes. La prueba `the_embedded_resource_digests_and_json_checks_cannot_fail`
    // comprueba las dos premisas. Si alguien cambia el recurso sin cambiar la
    // constante, esta guarda salta en toda ejecución, no solo ante una entrada.
    let observed = format!("{:x}", Sha256::digest(RESOURCE_BYTES));
    if observed != GPS_HISEQX_PCRFREE_PE150_SHA256 {
        return Err(rejected(F4TechnicalProfileError::ResourceDigestMismatch {
            observed,
        }));
    }
    let value: Value = serde_json::from_slice(RESOURCE_BYTES).map_err(|_| {
        rejected(F4TechnicalProfileError::InvalidResource {
            detail: "JSON no legible",
        })
    })?;
    validate_resource(&value)
}

/// Comprueba el recurso hermano de la preparación amplificada.
pub fn verify_pcr_profile_resource() -> Result<(), F4TechnicalProfileError> {
    // Inalcanzable, igual que en `verify_initial_profile_resource`: bytes
    // incrustados contra una constante (`the_embedded_resource_digests_and_json_checks_cannot_fail`).
    let observed = format!("{:x}", Sha256::digest(PCR_RESOURCE_BYTES));
    if observed != GPS_HISEQX_PCR_PE150_SHA256 {
        return Err(rejected(F4TechnicalProfileError::ResourceDigestMismatch {
            observed,
        }));
    }
    let value: Value = serde_json::from_slice(PCR_RESOURCE_BYTES).map_err(|_| {
        rejected(F4TechnicalProfileError::InvalidResource {
            detail: "JSON no legible",
        })
    })?;
    validate_pcr_resource(&value)
}

fn validate_resource(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let root = object(value, "raíz")?;
    exact_fields(
        root,
        &[
            "schema_version",
            "profile_id",
            "classification",
            "description",
            "external_reference",
            "read",
            "insert",
            "quality",
            "error",
            "reproducibility",
        ],
        "raíz",
    )?;
    required_text(root, "schema_version", "1")?;
    required_text(root, "profile_id", GPS_HISEQX_PCRFREE_PE150_V1)?;
    required_text(
        root,
        "classification",
        "approximate_synthetic_not_empirical_calibration",
    )?;
    validate_external_reference(required(root, "external_reference")?)?;
    validate_read(required(root, "read")?)?;
    validate_insert(required(root, "insert")?)?;
    validate_quality(required(root, "quality")?)?;
    validate_error(required(root, "error")?)?;
    validate_reproducibility(required(root, "reproducibility")?, "f4-technical-v3")
}

/// Valida el recurso hermano.
///
/// Las cuatro secciones que el plan exige **idénticas** —lectura, inserto,
/// calidad y error— se comprueban con las **mismas** funciones que el recurso
/// vigente, no con copias: si alguien tocara un valor de una de ellas en el
/// hermano, este validador lo rechazaría aquí, que es exactamente lo que
/// significa «idéntico» en el contrato.
fn validate_pcr_resource(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let root = object(value, "raíz")?;
    exact_fields(
        root,
        &[
            "schema_version",
            "profile_id",
            "classification",
            "description",
            "external_reference",
            "preparation",
            "read",
            "insert",
            "quality",
            "error",
            "amplification",
            "reproducibility",
        ],
        "raíz",
    )?;
    required_text(root, "schema_version", "1")?;
    required_text(root, "profile_id", GPS_HISEQX_PCR_PE150_V1)?;
    required_text(
        root,
        "classification",
        "approximate_synthetic_not_empirical_calibration",
    )?;
    required_text(root, "preparation", LibraryPreparation::Pcr.as_str())?;
    validate_external_reference(required(root, "external_reference")?)?;
    validate_read(required(root, "read")?)?;
    validate_insert(required(root, "insert")?)?;
    validate_quality(required(root, "quality")?)?;
    validate_error(required(root, "error")?)?;
    validate_amplification(required(root, "amplification")?)?;
    validate_reproducibility(required(root, "reproducibility")?, "f4-technical-v5")
}

/// Los dominios de la amplificación, comparados **contra las constantes de
/// contrato**, nunca contra literales repetidos aquí.
fn validate_amplification(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let object = object(value, "amplification")?;
    exact_fields(
        object,
        &[
            "algorithm",
            "duplicate_fraction_inversion",
            "duplicate_fraction_inversion_iterations",
            "maximum_cycles",
            "maximum_early_cycles",
            "maximum_family_size",
            "maximum_gc_attempts",
            "maximum_efficiency_classes",
            "maximum_mean_pairs_per_molecule",
        ],
        "amplification",
    )?;
    required_text(object, "algorithm", PCR_AMPLIFICATION_ALGORITHM)?;
    required_text(
        object,
        "duplicate_fraction_inversion",
        PCR_DUPLICATE_FRACTION_INVERSION_METHOD,
    )?;
    required_u64(
        object,
        "duplicate_fraction_inversion_iterations",
        u64::from(PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS),
    )?;
    required_u64(object, "maximum_cycles", u64::from(PCR_MAXIMUM_CYCLES))?;
    required_u64(
        object,
        "maximum_early_cycles",
        u64::from(PCR_MAXIMUM_EARLY_CYCLES),
    )?;
    required_u64(
        object,
        "maximum_family_size",
        u64::from(PCR_MAXIMUM_FAMILY_SIZE),
    )?;
    required_u64(
        object,
        "maximum_gc_attempts",
        u64::from(PCR_MAXIMUM_GC_ATTEMPTS),
    )?;
    required_u64(
        object,
        "maximum_efficiency_classes",
        u64::from(PCR_MAXIMUM_EFFICIENCY_CLASSES),
    )?;
    // El máximo de parejas por molécula es el único dominio con parte
    // fraccionaria posible, de modo que se compara como número, no como entero.
    let declared = object
        .get("maximum_mean_pairs_per_molecule")
        .and_then(Value::as_f64);
    if declared != Some(PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE) {
        return invalid("amplification.maximum_mean_pairs_per_molecule no coincide con M_MAX");
    }
    Ok(())
}

fn validate_external_reference(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let object = object(value, "external_reference")?;
    exact_fields(object, &["name", "use"], "external_reference")?;
    required_text(object, "name", "ART HSXn")?;
    required_text(
        object,
        "use",
        "comparative_reference_only_not_incorporated_or_redistributed",
    )
}

fn validate_read(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let object = object(value, "read")?;
    exact_fields(
        object,
        &[
            "layout",
            "orientation",
            "common_length_bases",
            "minimum_supported_length_bases",
            "maximum_supported_length_bases",
        ],
        "read",
    )?;
    required_text(object, "layout", "paired_end")?;
    required_text(object, "orientation", "FR")?;
    required_u64(object, "common_length_bases", 150)?;
    required_u64(object, "minimum_supported_length_bases", 1)?;
    required_u64(object, "maximum_supported_length_bases", 150)
}

fn validate_insert(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let object = object(value, "insert")?;
    exact_fields(
        object,
        &[
            "distribution",
            "mean_bases",
            "standard_deviation_bases",
            "minimum_bases",
            "maximum_bases",
            "integer_rounding",
            "uniform",
            "no_cached_normal_variate",
        ],
        "insert",
    )?;
    required_text(object, "distribution", "normal_truncated_box_muller_v1")?;
    required_u64(object, "mean_bases", 350)?;
    required_u64(object, "standard_deviation_bases", 20)?;
    required_u64(object, "minimum_bases", 250)?;
    required_u64(object, "maximum_bases", 450)?;
    required_text(object, "integer_rounding", "nearest_ties_away_from_zero")?;
    required_text(object, "uniform", "open_53_bit_v1")?;
    if object.get("no_cached_normal_variate") != Some(&Value::Bool(true)) {
        return invalid("insert.no_cached_normal_variate debe ser true");
    }
    Ok(())
}

fn validate_quality(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let object = object(value, "quality")?;
    exact_fields(object, &["curve_kind", "r1", "r2"], "quality")?;
    required_text(
        object,
        "curve_kind",
        "piecewise_linear_phred_integer_round_half_up_v1",
    )?;
    validate_quality_segments(
        required(object, "r1")?,
        &[(1, 10, 36, 36), (11, 100, 36, 32), (101, 150, 32, 28)],
        "quality.r1",
    )?;
    validate_quality_segments(
        required(object, "r2")?,
        &[(1, 10, 35, 35), (11, 100, 35, 30), (101, 150, 30, 25)],
        "quality.r2",
    )
}

fn validate_quality_segments(
    value: &Value,
    expected: &[(u64, u64, u64, u64)],
    name: &'static str,
) -> Result<(), F4TechnicalProfileError> {
    let segments = value.as_array().ok_or_else(|| {
        rejected(F4TechnicalProfileError::InvalidResource {
            detail: "los segmentos de calidad deben ser una lista",
        })
    })?;
    if segments.len() != expected.len() {
        return invalid("el número de segmentos de calidad no coincide");
    }
    for (segment, (first_cycle, last_cycle, first_phred, last_phred)) in
        segments.iter().zip(expected)
    {
        let object = object(segment, name)?;
        exact_fields(
            object,
            &["first_cycle", "last_cycle", "first_phred", "last_phred"],
            name,
        )?;
        if object.get("first_cycle").and_then(Value::as_u64) != Some(*first_cycle)
            || object.get("last_cycle").and_then(Value::as_u64) != Some(*last_cycle)
            || object.get("first_phred").and_then(Value::as_u64) != Some(*first_phred)
            || object.get("last_phred").and_then(Value::as_u64) != Some(*last_phred)
        {
            return invalid("un segmento de calidad no coincide con DTPD-106");
        }
    }
    Ok(())
}

fn validate_error(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let object = object(value, "error")?;
    exact_fields(
        object,
        &[
            "probability_from_phred",
            "alternate_base",
            "substitution_matrix",
        ],
        "error",
    )?;
    required_text(object, "probability_from_phred", "10^(-Q/10)")?;
    required_text(object, "alternate_base", "uniform_among_other_acgt")?;
    required_text(object, "substitution_matrix", "deferred_dtf_053")
}

fn validate_reproducibility(
    value: &Value,
    catalog: &'static str,
) -> Result<(), F4TechnicalProfileError> {
    let object = object(value, "reproducibility")?;
    exact_fields(
        object,
        &["rng_catalog", "profile_algorithm_version"],
        "reproducibility",
    )?;
    required_text(object, "rng_catalog", catalog)?;
    required_text(object, "profile_algorithm_version", "1")
}

fn object<'a>(
    value: &'a Value,
    name: &'static str,
) -> Result<&'a serde_json::Map<String, Value>, F4TechnicalProfileError> {
    value.as_object().ok_or_else(|| {
        rejected(F4TechnicalProfileError::InvalidResource {
            detail: match name {
                "raíz" => "la raíz debe ser un objeto",
                _ => "un objeto requerido no es un objeto JSON",
            },
        })
    })
}

fn required<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &'static str,
) -> Result<&'a Value, F4TechnicalProfileError> {
    // Inalcanzable: todo campo que se pide con `required` figura en la lista de
    // `exact_fields` del mismo objeto, y esa comprobación va siempre antes. Lo
    // ejecuta `a_missing_required_field_is_always_caught_by_the_exact_fields_check_first`.
    object.get(field).ok_or_else(|| {
        rejected(F4TechnicalProfileError::InvalidResource {
            detail: "falta un campo requerido",
        })
    })
}

fn exact_fields(
    object: &serde_json::Map<String, Value>,
    expected: &[&str],
    name: &'static str,
) -> Result<(), F4TechnicalProfileError> {
    if object.len() != expected.len() || expected.iter().any(|field| !object.contains_key(*field)) {
        return Err(rejected(F4TechnicalProfileError::InvalidResource {
            detail: match name {
                "raíz" => "los campos de raíz no coinciden con el recurso normativo",
                _ => "los campos de un objeto no coinciden con el recurso normativo",
            },
        }));
    }
    Ok(())
}

fn required_text(
    object: &serde_json::Map<String, Value>,
    field: &'static str,
    expected: &'static str,
) -> Result<(), F4TechnicalProfileError> {
    if object.get(field).and_then(Value::as_str) != Some(expected) {
        return invalid("un literal de texto no coincide con DTPD-106");
    }
    Ok(())
}

fn required_u64(
    object: &serde_json::Map<String, Value>,
    field: &'static str,
    expected: u64,
) -> Result<(), F4TechnicalProfileError> {
    if object.get(field).and_then(Value::as_u64) != Some(expected) {
        return invalid("un literal numérico no coincide con DTPD-106");
    }
    Ok(())
}

#[track_caller]
fn invalid(detail: &'static str) -> Result<(), F4TechnicalProfileError> {
    crate::guard_trace::record();
    Err(rejected(F4TechnicalProfileError::InvalidResource {
        detail,
    }))
}

#[cfg(test)]
mod tests {
    //! Pruebas de las guardas del validador de los recursos incorporados.
    //!
    //! `verify_*_profile_resource` solo ve los bytes incrustados, así que desde
    //! fuera no hay forma de presentarle un recurso alterado. Estas pruebas
    //! llaman a los validadores privados con el recurso real **alterado en un
    //! solo punto** y exigen el motivo exacto del rechazo, para que una guarda
    //! rota no la tape otra anterior (DT-001, DT-032).

    use super::*;
    use serde_json::json;

    const EXACT_FIELDS_ROOT: &str = "los campos de raíz no coinciden con el recurso normativo";
    const EXACT_FIELDS_OBJECT: &str =
        "los campos de un objeto no coinciden con el recurso normativo";
    const TEXT_LITERAL: &str = "un literal de texto no coincide con DTPD-106";
    const NUMBER_LITERAL: &str = "un literal numérico no coincide con DTPD-106";

    #[derive(Clone, Copy, Debug)]
    enum Which {
        PcrFree,
        Pcr,
    }

    fn resource(which: Which) -> Value {
        let bytes = match which {
            Which::PcrFree => RESOURCE_BYTES,
            Which::Pcr => PCR_RESOURCE_BYTES,
        };
        serde_json::from_slice(bytes).expect("the embedded resource is JSON")
    }

    fn validate(which: Which, value: &Value) -> Result<(), F4TechnicalProfileError> {
        match which {
            Which::PcrFree => validate_resource(value),
            Which::Pcr => validate_pcr_resource(value),
        }
    }

    fn invalid_detail(detail: &'static str) -> Result<(), F4TechnicalProfileError> {
        Err(F4TechnicalProfileError::InvalidResource { detail })
    }

    #[test]
    fn the_unaltered_resources_pass_their_own_validators() {
        for which in [Which::PcrFree, Which::Pcr] {
            assert_eq!(validate(which, &resource(which)), Ok(()), "{which:?}");
        }
    }

    #[test]
    fn every_single_point_alteration_of_a_resource_is_rejected_for_its_own_reason() {
        // Cada caso: recurso, puntero JSON, valor nuevo (`None` = quitar la
        // clave) y motivo exacto. El puntero apunta a **un** sitio, de modo
        // que el único rechazo posible es el de la guarda que lo vigila.
        let cases: Vec<(&str, Which, &str, Option<Value>, &str)> = vec![
            // Línea de `exact_fields`, con la raíz y con un objeto anidado.
            (
                "extra root field",
                Which::PcrFree,
                "/unexpected",
                Some(json!(1)),
                EXACT_FIELDS_ROOT,
            ),
            (
                "missing root field",
                Which::Pcr,
                "/preparation",
                None,
                EXACT_FIELDS_ROOT,
            ),
            (
                "extra nested field",
                Which::PcrFree,
                "/read/unexpected",
                Some(json!(1)),
                EXACT_FIELDS_OBJECT,
            ),
            (
                "missing nested field",
                Which::Pcr,
                "/amplification/algorithm",
                None,
                EXACT_FIELDS_OBJECT,
            ),
            // Mismo número de campos, uno con otro nombre: solo lo caza la
            // mitad de `exact_fields` que busca cada campo esperado.
            (
                "renamed nested field",
                Which::PcrFree,
                "/read",
                Some(json!({
                    "layoutx": "paired_end",
                    "orientation": "FR",
                    "common_length_bases": 150,
                    "minimum_supported_length_bases": 1,
                    "maximum_supported_length_bases": 150
                })),
                EXACT_FIELDS_OBJECT,
            ),
            // Línea de `object`: un objeto requerido que no lo es.
            (
                "section that is not an object",
                Which::PcrFree,
                "/read",
                Some(json!(5)),
                "un objeto requerido no es un objeto JSON",
            ),
            // Línea de `required_text`.
            (
                "wrong schema version",
                Which::PcrFree,
                "/schema_version",
                Some(json!("2")),
                TEXT_LITERAL,
            ),
            // Línea de `required_u64`.
            (
                "wrong common read length",
                Which::PcrFree,
                "/read/common_length_bases",
                Some(json!(149)),
                NUMBER_LITERAL,
            ),
            // Un campo cada vez de `validate_quality_segments`: `last_phred`
            // lo cubre el caso de más abajo; aquí los otros tres.
            (
                "quality segment first cycle altered",
                Which::PcrFree,
                "/quality/r1/1/first_cycle",
                Some(json!(12)),
                "un segmento de calidad no coincide con DTPD-106",
            ),
            (
                "quality segment last cycle altered",
                Which::PcrFree,
                "/quality/r1/1/last_cycle",
                Some(json!(99)),
                "un segmento de calidad no coincide con DTPD-106",
            ),
            (
                "quality segment first phred altered",
                Which::PcrFree,
                "/quality/r1/1/first_phred",
                Some(json!(37)),
                "un segmento de calidad no coincide con DTPD-106",
            ),
            // Dos literales concretos que pasan por los ayudantes compartidos.
            (
                "wrong read orientation",
                Which::PcrFree,
                "/read/orientation",
                Some(json!("RF")),
                TEXT_LITERAL,
            ),
            (
                "wrong maximum PCR cycles",
                Which::Pcr,
                "/amplification/maximum_cycles",
                Some(json!(29)),
                NUMBER_LITERAL,
            ),
            // Línea de `no_cached_normal_variate`.
            (
                "cached normal variate allowed",
                Which::PcrFree,
                "/insert/no_cached_normal_variate",
                Some(json!(false)),
                "insert.no_cached_normal_variate debe ser true",
            ),
            // Tres líneas de `validate_quality_segments`.
            (
                "quality segments not a list",
                Which::PcrFree,
                "/quality/r1",
                Some(json!({})),
                "los segmentos de calidad deben ser una lista",
            ),
            (
                "quality segment missing",
                Which::Pcr,
                "/quality/r2/2",
                None,
                "el número de segmentos de calidad no coincide",
            ),
            (
                "quality segment value altered",
                Which::PcrFree,
                "/quality/r1/1/last_phred",
                Some(json!(31)),
                "un segmento de calidad no coincide con DTPD-106",
            ),
            // Línea de M_MAX: es el único dominio que se compara como número.
            (
                "wrong maximum mean pairs per molecule",
                Which::Pcr,
                "/amplification/maximum_mean_pairs_per_molecule",
                Some(json!(63.5)),
                "amplification.maximum_mean_pairs_per_molecule no coincide con M_MAX",
            ),
        ];
        for (name, which, pointer, replacement, detail) in cases {
            let mut value = resource(which);
            let (parent, key) = pointer.rsplit_once('/').expect("pointer has a parent");
            let parent = value.pointer_mut(parent).expect("parent exists");
            match (&replacement, parent) {
                (Some(new), Value::Object(map)) => {
                    map.insert(key.to_owned(), new.clone());
                }
                (None, Value::Object(map)) => {
                    map.remove(key).expect("key exists");
                }
                (Some(new), Value::Array(items)) => {
                    let index: usize = key.parse().expect("array index");
                    items[index] = new.clone();
                }
                (None, Value::Array(items)) => {
                    items.remove(key.parse::<usize>().expect("array index"));
                }
                _ => panic!("{name}: the pointer parent is a scalar"),
            }
            assert_eq!(validate(which, &value), invalid_detail(detail), "{name}");
        }
    }

    #[test]
    fn a_resource_that_is_not_an_object_is_rejected_at_the_root() {
        for which in [Which::PcrFree, Which::Pcr] {
            assert_eq!(
                validate(which, &json!([])),
                invalid_detail("la raíz debe ser un objeto"),
                "{which:?}"
            );
        }
    }

    /// `required` no puede disparar nunca su rechazo: todo campo que se pide
    /// con `required(..)` está en la lista de `exact_fields` del mismo objeto
    /// y esa comprobación va antes. Esta prueba lo **ejecuta** quitando cada
    /// campo de la raíz y de `quality`: el rechazo es siempre el de
    /// `exact_fields`, nunca «falta un campo requerido». No discrimina la
    /// guarda de `required` (no se puede alcanzar); documenta por qué.
    #[test]
    fn a_missing_required_field_is_always_caught_by_the_exact_fields_check_first() {
        for which in [Which::PcrFree, Which::Pcr] {
            let original = resource(which);
            for field in original.as_object().expect("root").keys() {
                let mut value = original.clone();
                value.as_object_mut().expect("root").remove(field);
                assert_eq!(
                    validate(which, &value),
                    invalid_detail(EXACT_FIELDS_ROOT),
                    "{which:?} without {field}"
                );
            }
            for field in ["r1", "r2", "curve_kind"] {
                let mut value = original.clone();
                value["quality"]
                    .as_object_mut()
                    .expect("quality")
                    .remove(field);
                assert_eq!(
                    validate(which, &value),
                    invalid_detail(EXACT_FIELDS_OBJECT),
                    "{which:?} quality without {field}"
                );
            }
        }
    }

    /// Los rechazos por huella y por JSON ilegible de
    /// `verify_*_profile_resource` son inalcanzables desde cualquier entrada:
    /// comparan los bytes **incrustados en el binario** con una constante, y
    /// el JSON se lee de esos mismos bytes. Esta prueba comprueba las dos
    /// premisas (la huella coincide y el JSON se lee), que son las que hacen
    /// imposible disparar esas cuatro guardas sin cambiar el recurso.
    #[test]
    fn the_embedded_resource_digests_and_json_checks_cannot_fail() {
        assert_eq!(
            format!("{:x}", Sha256::digest(RESOURCE_BYTES)),
            GPS_HISEQX_PCRFREE_PE150_SHA256
        );
        assert_eq!(
            format!("{:x}", Sha256::digest(PCR_RESOURCE_BYTES)),
            GPS_HISEQX_PCR_PE150_SHA256
        );
        assert!(serde_json::from_slice::<Value>(RESOURCE_BYTES).is_ok());
        assert!(serde_json::from_slice::<Value>(PCR_RESOURCE_BYTES).is_ok());
        assert_eq!(verify_initial_profile_resource(), Ok(()));
        assert_eq!(verify_pcr_profile_resource(), Ok(()));
    }
}
