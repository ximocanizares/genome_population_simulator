//! Perfil técnico incorporado de F4.
//!
//! Este módulo es deliberadamente una frontera de datos: valida el único
//! recurso versionado aprobado por DTPD-106 y expone los hechos que necesitarán
//! la construcción de pares. No selecciona el perfil desde TOML ni materializa
//! calidades, errores, contaminación o FASTQ; esas responsabilidades pertenecen
//! a hitos posteriores.

use std::fmt;

use serde_json::Value;
use sha2::{Digest, Sha256};

/// Identificador del único perfil incorporado inicialmente en F4.
pub const GPS_HISEQX_PCRFREE_PE150_V1: &str = "gps-hiseqx-pcrfree-pe150-v1";
/// Versión semántica del perfil, separada de la versión de su recurso JSON.
pub const GPS_HISEQX_PCRFREE_PE150_PROFILE_VERSION: u8 = 1;
/// Ruta reproducible del recurso dentro del árbol fuente.
pub const GPS_HISEQX_PCRFREE_PE150_RESOURCE: &str = "resources/f4/gps-hiseqx-pcrfree-pe150-v1.json";
/// SHA-256 normativo fijado por DTPD-106.
pub const GPS_HISEQX_PCRFREE_PE150_SHA256: &str =
    "3d8878ab920237dae03705276a26dbb0823941830d531374d35aa1c540f6deea";

const RESOURCE_BYTES: &[u8] = include_bytes!("../resources/f4/gps-hiseqx-pcrfree-pe150-v1.json");

/// Preparación de biblioteca fija para el perfil inicial.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LibraryPreparation {
    PcrFree,
}

impl LibraryPreparation {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::PcrFree => "pcr_free",
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

/// Resuelve el único perfil de F4 después de comprobar el recurso incorporado.
///
/// La frontera no interpreta aún un campo de configuración: permite que F4.1.2
/// y F4.1.3 reciban un perfil tipado sin hacer del identificador una cadena
/// implícita en sus algoritmos.
pub fn resolve_f4_technical_profile(
    profile_id: &str,
) -> Result<F4TechnicalProfile, F4TechnicalProfileError> {
    if profile_id != GPS_HISEQX_PCRFREE_PE150_V1 {
        return Err(F4TechnicalProfileError::UnknownProfile {
            requested: profile_id.to_owned(),
        });
    }
    verify_initial_profile_resource()?;
    Ok(INITIAL_PROFILE)
}

/// Comprueba el recurso literal y su SHA-256 antes de que un consumidor F4 lo
/// use. Se mantiene separada de la configuración hasta la transición atómica
/// de contratos de F4.1.3.
pub fn verify_initial_profile_resource() -> Result<(), F4TechnicalProfileError> {
    let observed = format!("{:x}", Sha256::digest(RESOURCE_BYTES));
    if observed != GPS_HISEQX_PCRFREE_PE150_SHA256 {
        return Err(F4TechnicalProfileError::ResourceDigestMismatch { observed });
    }
    let value: Value = serde_json::from_slice(RESOURCE_BYTES).map_err(|_| {
        F4TechnicalProfileError::InvalidResource {
            detail: "JSON no legible",
        }
    })?;
    validate_resource(&value)
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
    validate_reproducibility(required(root, "reproducibility")?)
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
    let segments = value
        .as_array()
        .ok_or(F4TechnicalProfileError::InvalidResource {
            detail: "los segmentos de calidad deben ser una lista",
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

fn validate_reproducibility(value: &Value) -> Result<(), F4TechnicalProfileError> {
    let object = object(value, "reproducibility")?;
    exact_fields(
        object,
        &["rng_catalog", "profile_algorithm_version"],
        "reproducibility",
    )?;
    required_text(object, "rng_catalog", "f4-technical-v3")?;
    required_text(object, "profile_algorithm_version", "1")
}

fn object<'a>(
    value: &'a Value,
    name: &'static str,
) -> Result<&'a serde_json::Map<String, Value>, F4TechnicalProfileError> {
    value
        .as_object()
        .ok_or(F4TechnicalProfileError::InvalidResource {
            detail: match name {
                "raíz" => "la raíz debe ser un objeto",
                _ => "un objeto requerido no es un objeto JSON",
            },
        })
}

fn required<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &'static str,
) -> Result<&'a Value, F4TechnicalProfileError> {
    object
        .get(field)
        .ok_or(F4TechnicalProfileError::InvalidResource {
            detail: "falta un campo requerido",
        })
}

fn exact_fields(
    object: &serde_json::Map<String, Value>,
    expected: &[&str],
    name: &'static str,
) -> Result<(), F4TechnicalProfileError> {
    if object.len() != expected.len() || expected.iter().any(|field| !object.contains_key(*field)) {
        return Err(F4TechnicalProfileError::InvalidResource {
            detail: match name {
                "raíz" => "los campos de raíz no coinciden con el recurso normativo",
                _ => "los campos de un objeto no coinciden con el recurso normativo",
            },
        });
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

fn invalid(detail: &'static str) -> Result<(), F4TechnicalProfileError> {
    Err(F4TechnicalProfileError::InvalidResource { detail })
}
