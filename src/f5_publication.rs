//! F5.4.3. Configuración efectiva v5, manifiesto general v5, `f5-failure.json`
//! y publicación atómica.
//!
//! Este módulo **no decide nada biológico**: instancia los dos documentos que
//! DTC-089 cierra término a término y los publica, o no publica nada. Lo que
//! escribe procede de tres fuentes y de ninguna más: la admisión sellada
//! (F5.0), el estado que F5.1–F5.3 cerraron y los bytes que F5.4.1 acaba de
//! publicar, releídos por el oráculo independiente de F5.4.2.
//!
//! ## Cuatro interpretaciones del contrato, declaradas
//!
//! - **`effective-config.json` es un documento raíz, no un artefacto.** El
//!   objeto `artifacts` del manifiesto v5 es una lista cerrada de exactamente
//!   los siete TSV —`prefixItems` con `items:false`— y no tiene ningún campo
//!   donde enlazar la configuración efectiva. Publicarla dentro del catálogo
//!   la dejaría sin declarar en el manifiesto, que es justo lo que v4 evitaba
//!   con `effective_config_artifact_id`. Por eso viaja como documento raíz
//!   junto a `manifest.json`, escrita con la misma disciplina de temporal y
//!   renombrado.
//! - **Las rutas publicadas son relativas a la raíz de autoría**, que es el
//!   directorio del general TOML. `$defs/path` prohíbe la ruta absoluta, de
//!   modo que la ruta física del sistema de ficheros no es publicable: vive,
//!   completa y con su SHA-256 observado, en `resource_preflight.input_identities`,
//!   cuyo `normalized_path` no tiene patrón. `source_path` es por tanto el
//!   nombre del general, y escenario, perfil, máscara y referencia se componen
//!   con los enlaces relativos que los propios documentos declaran. Una ruta
//!   de autoría con caracteres fuera de `[A-Za-z0-9._-]` **no es
//!   representable** y se rechaza con `GPS-FMT-006` antes de abrir destino: es
//!   una consecuencia del esquema, no una restricción añadida aquí.
//! - **Los resúmenes MAN-70 son una entrada, no un cálculo.** `recounted_*`
//!   solo puede nacer de releer `outputs/truth.jsonl`, y REQ-MAN-70 exige que
//!   el recuento no reutilice acumuladores del productor. F5.4.2 ya tiene ese
//!   lector; este módulo le entrega la declaración y el árbol publicado y
//!   copia lo que devuelve. Nunca construye un `F5InitializationSummary`: sus
//!   campos son privados del oráculo, de modo que el tipo mismo impide que el
//!   productor se autocertifique.
//! - **Un fallo controlado no publica ningún artefacto biológico.** El
//!   manifiesto fallido enlaza exactamente `outputs/f5-failure.json` y la
//!   puerta física comprueba, sobre el árbol que se va a renombrar, que **no
//!   existe** ninguno de los siete TSV ni `outputs/truth.jsonl`, existieran o
//!   no como temporales. No es una promesa: es una comprobación.
//!
//! ## Terreno compartido con F5.4.1 y F5.4.2
//!
//! `F5_TSV_CONTRACTS`, `F5_TSV_SCHEMA_VERSION` y los cinco `F5_TRUTH_*` los
//! importan **los tres módulos**: el escritor los usa para producir, el
//! oráculo para juzgar y este para declarar. Son una sola definición en
//! `f5_artifacts.rs` precisamente para que no puedan divergir: si el contrato
//! cambiara un identificador, una clase o una ruta, los tres se romperían a la
//! vez y ninguno podría aprobar lo que otro escribió. Los enumerados de código
//! y etapa de un fallo, en cambio, viven aquí y el oráculo tiene los suyos:
//! son el contrato del **esquema**, no del escritor, y duplicarlos es lo que
//! permite que el oráculo no herede la lectura del productor.
//!
//! ## El oráculo queda fuera de este flujo, y por qué se dice
//!
//! `f5-resource-accounting-v9` no reserva ni un byte para el oráculo, cuyo
//! pico medido en F5.4.2 es del orden de 4,2 veces `working_memory_bytes`.
//! F5.4.3 **no lo ejecuta**: [`publish_f5_success`] recibe el recuento MAN-70
//! como una función que el llamador aporta, junto con la declaración que este
//! módulo construye. Quien la ejecute dentro de una ejecución productiva
//! —hoy, nadie: la vía `run` de F5 la abre F5.6— deberá antes darle un término
//! propio en la fórmula o sustituirla por un lector por flujo del registro
//! `haplotype_copy`. La condición no se resuelve aquí porque hacerlo exigiría
//! reabrir la contabilidad de recursos, cuyo literal es un `const` de
//! `effective-config-v5` al que el manifiesto se refiere por `$ref`. Sigue
//! declarada como deuda con destino F9.7: la corrección de v9 fue de la cota de
//! **salida** y no tocó los dos términos de memoria que faltan.

use std::collections::BTreeSet;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use serde::ser::{SerializeMap, SerializeSeq};
use serde::{Serialize, Serializer};
use serde_json::{Map, Value, json};
use sha2::{Digest, Sha256};

use crate::application::EffectiveF5Config;
use crate::configuration::DeclaredF5Target;
use crate::diagnostics::Diagnostic;
use crate::f5_artifacts::{
    F5_TRUTH_ARTIFACT_ID, F5_TRUTH_ARTIFACT_KIND, F5_TRUTH_PATH, F5_TRUTH_PROJECTION_VERSION,
    F5_TRUTH_SCHEMA_VERSION, F5_TSV_CONTRACTS, F5_TSV_SCHEMA_VERSION, F5ReferenceSequences,
    F5SuccessArtifacts, F5SuccessRequest, write_f5_success_artifacts,
};
use crate::f5_founders::{
    F5FounderAssignment, F5FounderIdentities, F5RealizedFStatus, F5SamplingMode,
};
use crate::f5_layout::{F5BestLayout, F5LayoutFailure, F5PhysicalLayout};
use crate::f5_oracle::{
    F5_MAN70_RECOUNT_SOURCE, F5DeclaredAttemptSummary, F5DeclaredBounds, F5DeclaredHaplotypeCount,
    F5DeclaredHaplotypeCountDistribution, F5DeclaredLdBin, F5DeclaredLink, F5DeclaredMode,
    F5DeclaredRealizedF, F5DeclaredSampling, F5DeclaredSpectrumBand, F5DeclaredSubpopulation,
    F5DeclaredTruthLink, F5InitializationSummary, F5PackageDeclaration,
};
use crate::f5_preflight::{
    F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS, F5_DEFAULT_MAX_HAPLOTYPE_CELLS,
    F5_DEFAULT_MAX_SWAP_PROPOSALS, F5ResolvedTarget, SealedF5Admission,
};
use crate::f5_search::{F5AttemptRecord, F5HaplotypeSearch, F5SearchFailure};
use crate::identity::{ExecutionId, F5_ID_CATALOG};
use crate::integrity::identify_file;
use crate::publication::ExecutionWorkspace;
use crate::randomness::{F5_RNG_CATALOG_V1, F5_RNG_CATALOG_V2, POISSON_COUNT_ALGORITHM};
use crate::results::{ArtifactRecord, json_document_bytes};

// ---------------------------------------------------------------------------
// Contratos literales
// ---------------------------------------------------------------------------

/// Documento raíz de la configuración efectiva publicada.
pub const F5_EFFECTIVE_CONFIG_PATH: &str = "effective-config.json";
/// Documento raíz del manifiesto general publicado.
pub const F5_MANIFEST_PATH: &str = "manifest.json";
/// `schema_version` de la raíz de ambos documentos v5.
pub const F5_DOCUMENT_SCHEMA_VERSION: &str = "5";
/// `schema_version` del objeto `f5` de ambos documentos.
pub const F5_OBJECT_SCHEMA_VERSION: &str = "1";
/// Diagnóstico publicado de un fallo controlado.
pub const F5_FAILURE_PATH: &str = "outputs/f5-failure.json";
pub const F5_FAILURE_ARTIFACT_ID: &str = "f5-00-failure-v1";
pub const F5_FAILURE_ARTIFACT_KIND: &str = "f5_failure";
pub const F5_FAILURE_SCHEMA_VERSION: &str = "1";
pub const F5_FAILURE_DIAGNOSTIC_KIND: &str = "f5_failure";
/// Rutas publicadas de los catálogos normativos F5.
pub const F5_RNG_CATALOG_V1_PATH: &str = "catalogs/f5-rng-catalog-v1.tsv";
pub const F5_RNG_CATALOG_V2_PATH: &str = "catalogs/f5-rng-catalog-v2.tsv";
pub const F5_ID_CATALOG_PATH: &str = "catalogs/id-catalog-v4.tsv";

const RNG_CATALOG_V1_BYTES: &[u8] = include_bytes!("../catalogs/f5-rng-catalog-v1.tsv");
const RNG_CATALOG_V2_BYTES: &[u8] = include_bytes!("../catalogs/f5-rng-catalog-v2.tsv");
const ID_CATALOG_BYTES: &[u8] = include_bytes!("../catalogs/id-catalog-v4.tsv");

/// Catálogo RNG **de la ejecución**: v2 solo cuando el perfil declara la
/// distribución del recuento haplotípico y, por tanto, el ámbito
/// `block_haplotype_count` puede abrirse. Un escenario que no la declara
/// publica el mismo catálogo, con la misma huella, que antes de que este
/// ámbito existiera.
#[must_use]
pub fn f5_rng_catalog_for(
    config: &EffectiveF5Config,
) -> (&'static str, &'static str, &'static [u8]) {
    if config.profile().haplotype_count_distribution().is_some() {
        (
            F5_RNG_CATALOG_V2,
            F5_RNG_CATALOG_V2_PATH,
            RNG_CATALOG_V2_BYTES,
        )
    } else {
        (
            F5_RNG_CATALOG_V1,
            F5_RNG_CATALOG_V1_PATH,
            RNG_CATALOG_V1_BYTES,
        )
    }
}

/// Los siete algoritmos exactos del objeto `algorithms`, en el orden de
/// DTC-089. Ninguno es configurable: F5 solo publica estos.
const F5_ALGORITHMS: [(&str, &str); 7] = [
    ("layout", "f5-layout-v1"),
    ("marginals", "f5-spectrum-v1"),
    ("haplotypes", "weighted-prefix-swap-v1"),
    ("ld", "phased-r2-neumaier-v1"),
    ("controlled", "largest-remainder-v1"),
    ("inbreeding", "diploid-f-v1"),
    ("truth_projection", F5_TRUTH_PROJECTION_VERSION),
];

/// Las nueve claves del objeto `f5` que el manifiesto **reutiliza exactamente**
/// del efectivo, según DTC-089.
const SHARED_OBJECTS: [&str; 9] = [
    "profile",
    "mask",
    "target",
    "algorithms",
    "rng_catalog",
    "id_catalog",
    "limits",
    "resource_preflight",
    "normalizations",
];

/// `recount_source` obligatorio de REQ-MAN-70, el mismo literal que el
/// oráculo escribe en cada resumen. No se redeclara: se importa, para que un
/// cambio de versión de proyección no pueda dejar dos literales distintos.
const RECOUNT_SOURCE: &str = F5_MAN70_RECOUNT_SOURCE;

/// Los cuatro conjuntos normalizados, en el orden cerrado del contrato.
const NORMALIZATION_ORDER: [&str; 4] = [
    "block_lengths",
    "snp_gaps",
    "snp_spectrum",
    "haplotype_weights",
];

#[track_caller]
fn document_error(path: &str, detail: impl Into<String>) -> Diagnostic {
    crate::guard_trace::record();
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path)
    .with("schema", "effective-config-v5.schema.json")
    .with("detail", detail.into())
}

#[track_caller]
fn manifest_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    crate::guard_trace::record();
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path.display().to_string())
    .with("schema", "manifest-general-v5.schema.json")
    .with("detail", detail.into())
}

/// Rechazo de un término que el diagnóstico de fallo y el manifiesto publican
/// distinto. Es una función y no una clausura para que lleve `#[track_caller]`:
/// cada uso es su propio punto de rechazo para la puerta de cobertura.
#[track_caller]
fn failure_mismatch(manifest_path: &Path, term: &str) -> Diagnostic {
    crate::guard_trace::record();
    manifest_error(
        manifest_path,
        format!("{F5_FAILURE_PATH}/{term}: el diagnóstico y el manifiesto discrepan"),
    )
}

// ---------------------------------------------------------------------------
// Rutas publicables
// ---------------------------------------------------------------------------

/// Convierte una ruta absoluta admitida en la ruta relativa a la raíz de
/// autoría que `$defs/path` admite.
///
/// La raíz de autoría es el directorio del general TOML. Una entrada fuera de
/// ese árbol, o con un componente fuera de `[A-Za-z0-9._-]`, no es
/// representable en el contrato publicado y se rechaza antes de abrir destino.
fn publishable_path(root: &Path, absolute: &Path, field: &str) -> Result<String, Diagnostic> {
    let relative = absolute.strip_prefix(root).map_err(|_| {
        document_error(
            F5_EFFECTIVE_CONFIG_PATH,
            format!(
                "{field}: {} no cuelga de la raíz de autoría {} y $defs/path prohíbe la ruta absoluta",
                absolute.display(),
                root.display()
            ),
        )
    })?;
    publishable_relative(relative, field)
}

fn publishable_relative(relative: &Path, field: &str) -> Result<String, Diagnostic> {
    let text = relative
        .to_str()
        .ok_or_else(|| document_error(F5_EFFECTIVE_CONFIG_PATH, format!("{field}: no es UTF-8")))?;
    let representable = !text.is_empty()
        && text.split('/').all(|component| {
            !component.is_empty()
                && component != "."
                && component != ".."
                && component
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
        });
    if !representable {
        return Err(document_error(
            F5_EFFECTIVE_CONFIG_PATH,
            format!("{field}: {text:?} no es representable en $defs/path"),
        ));
    }
    Ok(text.to_owned())
}

// ---------------------------------------------------------------------------
// Enlaces compartidos por los dos documentos
// ---------------------------------------------------------------------------

/// Los objetos que DTC-089 obliga a que el manifiesto reutilice **exactamente**
/// del efectivo, construidos una sola vez.
struct SharedContract {
    profile: Value,
    mask: Value,
    target: Value,
    algorithms: Value,
    rng_catalog: Value,
    id_catalog: Value,
    limits: Value,
    resource_preflight: Value,
    normalizations: Value,
}

impl SharedContract {
    fn build(config: &EffectiveF5Config, sealed: &SealedF5Admission) -> Result<Self, Diagnostic> {
        let root = config
            .general_path()
            .parent()
            .ok_or_else(|| {
                document_error(
                    F5_EFFECTIVE_CONFIG_PATH,
                    "el general admitido no tiene directorio de autoría",
                )
            })?
            .to_path_buf();
        let receipt = sealed.receipt();
        let metadata = receipt
            .identity(crate::f5_admission::InputRole::ProfileMetadata)
            .ok_or_else(|| {
                document_error(
                    F5_EFFECTIVE_CONFIG_PATH,
                    "la admisión sellada no conserva la identidad del metadata de perfil",
                )
            })?;
        let profile = json!({
            "profile_id": sealed.profile_id(),
            "path": publishable_path(&root, metadata.normalized_path(), "profile.path")?,
            "size_bytes": metadata.observed_size_bytes(),
            "sha256": metadata.observed_sha256(),
        });
        let mask = match receipt.identity(crate::f5_admission::InputRole::Mask) {
            None => Value::Null,
            Some(identity) => json!({
                "path": publishable_path(&root, identity.normalized_path(), "mask.path")?,
                "size_bytes": identity.observed_size_bytes(),
                "sha256": identity.observed_sha256(),
            }),
        };
        Ok(Self {
            profile,
            mask,
            target: target_json(config, sealed)?,
            algorithms: algorithms_json(),
            rng_catalog: {
                let (catalog_id, path, bytes) = f5_rng_catalog_for(config);
                catalog_json(catalog_id, path, bytes)
            },
            id_catalog: catalog_json(F5_ID_CATALOG, F5_ID_CATALOG_PATH, ID_CATALOG_BYTES),
            limits: limits_json(config, sealed)?,
            resource_preflight: sealed.to_contract_json(),
            normalizations: normalizations_json(config)?,
        })
    }

    /// Los nueve valores en el orden exacto de [`SHARED_OBJECTS`], que es la
    /// única lista de esas claves en el módulo: emparejarlos aquí impide que
    /// una clave nueva se añada en un sitio y no en el otro.
    const fn values(&self) -> [&Value; SHARED_OBJECTS.len()] {
        [
            &self.profile,
            &self.mask,
            &self.target,
            &self.algorithms,
            &self.rng_catalog,
            &self.id_catalog,
            &self.limits,
            &self.resource_preflight,
            &self.normalizations,
        ]
    }

    fn insert_into(&self, object: &mut Map<String, Value>) {
        for (key, value) in SHARED_OBJECTS.into_iter().zip(self.values()) {
            object.insert(key.to_owned(), value.clone());
        }
    }
}

fn algorithms_json() -> Value {
    let mut object = Map::new();
    for (key, value) in F5_ALGORITHMS {
        object.insert(key.to_owned(), json!(value));
    }
    Value::Object(object)
}

fn catalog_json(catalog_id: &str, path: &str, bytes: &[u8]) -> Value {
    json!({
        "catalog_id": catalog_id,
        "path": path,
        "size_bytes": bytes.len() as u64,
        "sha256": format!("{:x}", Sha256::digest(bytes)),
    })
}

fn target_json(
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
) -> Result<Value, Diagnostic> {
    let declared_value = match &config.declared_scenario().target {
        DeclaredF5Target::BlockCount(value) => json!(*value.value()),
        DeclaredF5Target::CoveredBases(value) => json!(*value.value()),
        DeclaredF5Target::CoverageFraction(value) => json!(value.literal()),
    };
    let (kind, numerator, denominator) = match sealed.resolved_target() {
        F5ResolvedTarget::BlockCount { .. } => ("block_count", Value::Null, Value::Null),
        F5ResolvedTarget::CoveredBases { .. } => ("covered_bases", Value::Null, Value::Null),
        F5ResolvedTarget::CoverageFraction {
            numerator,
            denominator,
        } => ("coverage_fraction", json!(numerator), json!(denominator)),
    };
    Ok(json!({
        "kind": kind,
        "declared_value": declared_value,
        "numerator": numerator,
        "denominator": denominator,
        "target_bases": match sealed.target_bases() {
            Some(bases) => json!(bases),
            None => Value::Null,
        },
    }))
}

/// Los siete límites del contrato, derivados de la autoría admitida.
///
/// Es `pub(crate)` porque el auditor del paquete de fallo los necesita **sin
/// leerlos del manifiesto**: si los tomara de ahí compararía la serialización
/// del productor consigo misma. La regla vive aquí una sola vez.
pub(crate) fn limits_json(
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
) -> Result<Value, Diagnostic> {
    let scenario = config.declared_scenario();
    let declared = |value: Option<&crate::configuration::Located<u64>>, fallback: u64| {
        value.map_or(fallback, |located| *located.value())
    };
    let max_layout_attempts = sealed.max_layout_attempts();
    let declared_layout = declared(
        scenario.max_layout_attempts.as_ref(),
        crate::f5_preflight::F5_DEFAULT_MAX_LAYOUT_ATTEMPTS,
    );
    if declared_layout != max_layout_attempts {
        return Err(document_error(
            F5_EFFECTIVE_CONFIG_PATH,
            format!(
                "limits.max_layout_attempts: el escenario declara {declared_layout} y la admisión selló {max_layout_attempts}"
            ),
        ));
    }
    Ok(json!({
        "max_layout_attempts": max_layout_attempts,
        "max_haplotype_attempts": declared(
            scenario.max_haplotype_attempts.as_ref(),
            F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS,
        ),
        "max_swap_proposals": declared(
            scenario.max_swap_proposals.as_ref(),
            F5_DEFAULT_MAX_SWAP_PROPOSALS,
        ),
        "max_haplotype_cells": declared(
            scenario.max_haplotype_cells.as_ref(),
            F5_DEFAULT_MAX_HAPLOTYPE_CELLS,
        ),
        "max_operation_units": sealed.operation_limit(),
        "max_working_memory_bytes": sealed.working_memory_limit_bytes(),
        "max_output_bytes": sealed.output_limit_bytes(),
    }))
}

fn normalizations_json(config: &EffectiveF5Config) -> Result<Value, Diagnostic> {
    let audits = config.profile().normalizations();
    if audits.len() != NORMALIZATION_ORDER.len() {
        return Err(document_error(
            F5_EFFECTIVE_CONFIG_PATH,
            format!(
                "normalizations: el perfil admitido conserva {} auditorías y el contrato exige {}",
                audits.len(),
                NORMALIZATION_ORDER.len()
            ),
        ));
    }
    let mut rows = Vec::with_capacity(audits.len());
    for (audit, expected) in audits.iter().zip(NORMALIZATION_ORDER) {
        if audit.set().as_str() != expected {
            return Err(document_error(
                F5_EFFECTIVE_CONFIG_PATH,
                format!(
                    "normalizations: el orden cerrado exige {expected} y el perfil trae {}",
                    audit.set().as_str()
                ),
            ));
        }
        rows.push(json!({
            "set_kind": audit.set().as_str(),
            "canonical_order": audit.set().canonical_order(),
            "algorithm": audit.algorithm(),
            "tolerance": audit.tolerance(),
            "input_values": named_numbers(audit.input_values()),
            "input_sum": audit.input_sum(),
            "input_delta": audit.input_delta(),
            "effective_values": named_numbers(audit.effective_values()),
            "output_sum": audit.output_sum(),
            "output_delta": audit.output_delta(),
        }));
    }
    Ok(Value::Array(rows))
}

fn named_numbers(values: &[crate::f5_profile::F5NamedProbability]) -> Value {
    Value::Array(
        values
            .iter()
            .map(|value| json!({"key": value.key(), "value": value.value()}))
            .collect(),
    )
}

// ---------------------------------------------------------------------------
// Listas comunes: subpoblaciones y muestreo
// ---------------------------------------------------------------------------

fn optional_number(value: Option<f64>) -> Value {
    value.map_or(Value::Null, |value| json!(value))
}

/// `sampling` del efectivo y del manifiesto: la misma lista, una fila por
/// bloque×subpoblación, ordenada por `population_id` y `block_id`.
fn sampling_json(
    declaration: &F5PackageDeclaration,
    assignment: &F5FounderAssignment,
) -> Result<Value, Diagnostic> {
    let mut rows = Vec::with_capacity(assignment.sampling().len());
    for (declared, row) in declaration.sampling.iter().zip(assignment.sampling()) {
        let block_index = u32::try_from(row.block_index()).map_err(|_| {
            document_error(
                F5_EFFECTIVE_CONFIG_PATH,
                format!(
                    "sampling.block_index: {} excede el dominio del contrato",
                    row.block_index()
                ),
            )
        })?;
        rows.push(json!({
            "population_id": declared.population_id,
            "block_id": declared.block_id,
            "block_index": block_index,
            "requested_mode": declared.requested_mode.as_str(),
            "effective_mode": declared.effective_mode.as_str(),
            "source": row.source().as_str(),
            "controlled_max_abs_error_declared": optional_number(
                declared.controlled_max_abs_error_declared,
            ),
            "controlled_max_abs_error_effective": optional_number(
                declared.controlled_max_abs_error_effective,
            ),
        }));
    }
    Ok(Value::Array(rows))
}

/// `subpopulations` del efectivo: censo completo con claves fundadoras.
fn effective_subpopulations_json(assignment: &F5FounderAssignment) -> Value {
    Value::Array(
        assignment
            .rosters()
            .iter()
            .map(|roster| {
                json!({
                    "population_id": roster.population_id().as_str(),
                    "individual_count": roster.individual_count(),
                    "inbreeding_f": optional_number(roster.inbreeding_f()),
                    "inbreeding_f_source": roster.inbreeding_f_source().as_str(),
                    "sampling": roster.default_sampling().as_str(),
                    "founders": Value::Array(
                        roster
                            .individuals()
                            .iter()
                            .map(|individual| {
                                json!({
                                    "individual_id": individual.individual_id().as_str(),
                                    "founder_key": individual
                                        .founder_key()
                                        .map_or(Value::Null, |key| json!(key.as_str())),
                                })
                            })
                            .collect(),
                    ),
                })
            })
            .collect(),
    )
}

// ---------------------------------------------------------------------------
// Documentos v5
// ---------------------------------------------------------------------------

fn effective_document(
    request: &F5PublicationRequest<'_>,
    declaration: &F5PackageDeclaration,
    shared: &SharedContract,
) -> Result<Value, Diagnostic> {
    let (config, sealed, assignment) = (request.config, request.sealed, request.assignment);
    let mut f5 = Map::new();
    f5.insert("schema_version".to_owned(), json!(F5_OBJECT_SCHEMA_VERSION));
    f5.insert("scenario_id".to_owned(), json!(config.scenario_id()));
    f5.insert("reference_id".to_owned(), json!(config.reference_id()));
    f5.insert(
        "reference_sha256".to_owned(),
        json!(sealed.reference_sha256()),
    );
    f5.insert("seed".to_owned(), json!(config.seed()));
    f5.insert("ploidy".to_owned(), json!(config.ploidy()));
    f5.insert(
        "founder_count".to_owned(),
        json!(*config.declared_scenario().founder_count.value()),
    );
    f5.insert(
        "subpopulations".to_owned(),
        effective_subpopulations_json(assignment),
    );
    f5.insert(
        "sampling".to_owned(),
        sampling_json(declaration, assignment)?,
    );
    shared.insert_into(&mut f5);
    let source = config.general_path().file_name().ok_or_else(|| {
        document_error(
            F5_EFFECTIVE_CONFIG_PATH,
            "source_path: el general admitido no tiene nombre de fichero",
        )
    })?;
    Ok(json!({
        "effective_config_kind": "gps_effective_config",
        "schema_version": F5_DOCUMENT_SCHEMA_VERSION,
        "source_path": publishable_relative(Path::new(source), "source_path")?,
        "seed": config.seed(),
        "f5": Value::Object(f5),
    }))
}

/// Serializa los resúmenes MAN-70 **por flujo**, sin construir ningún árbol.
///
/// Es el término dominante del manifiesto: crece con `B*N_pop*S` y, con la
/// característica `arbitrary_precision` activa, cada número de un
/// `serde_json::Value` es a su vez un objeto de un miembro. Materializarlos
/// costaba del orden de catorce veces el texto publicado. Aquí cada fila se
/// emite directamente al serializador, en el mismo orden alfabético de claves
/// que produce `serde_json::Map` —que sin `preserve_order` es un `BTreeMap`—,
/// de modo que los bytes publicados no cambian.
struct Man70Rows<'a>(&'a [F5InitializationSummary]);

impl Serialize for Man70Rows<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut rows = serializer.serialize_seq(Some(self.0.len()))?;
        for summary in self.0 {
            rows.serialize_element(&Man70Row(summary))?;
        }
        rows.end()
    }
}

struct Man70Row<'a>(&'a F5InitializationSummary);

/// Un vector `{dosage,value}` emitido por flujo.
struct Dosages<'a, T>(&'a [T]);

impl<T: Serialize> Serialize for Dosages<'_, T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut rows = serializer.serialize_seq(Some(self.0.len()))?;
        for (dosage, value) in self.0.iter().enumerate() {
            rows.serialize_element(&Dosage {
                dosage: dosage as u64,
                value,
            })?;
        }
        rows.end()
    }
}

#[derive(Serialize)]
struct Dosage<'a, T> {
    dosage: u64,
    value: &'a T,
}

impl Serialize for Man70Row<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let summary = self.0;
        let mut row = serializer.serialize_map(Some(18))?;
        row.serialize_entry("block_id", summary.block_id())?;
        row.serialize_entry("effective_mode", summary.effective_mode().as_str())?;
        row.serialize_entry(
            "expected_genotype_counts",
            &Dosages(summary.expected_genotype_counts()),
        )?;
        row.serialize_entry(
            "expected_genotype_frequencies",
            &Dosages(summary.expected_genotype_frequencies()),
        )?;
        row.serialize_entry(
            "expected_genotype_probabilities",
            &Dosages(summary.expected_genotype_probabilities()),
        )?;
        row.serialize_entry("input_alt_af", &summary.input_alt_af())?;
        row.serialize_entry("model", summary.model().as_str())?;
        row.serialize_entry("ploidy", &summary.ploidy())?;
        row.serialize_entry("population_id", summary.population_id())?;
        row.serialize_entry("realized_alt_af", &summary.realized_alt_af())?;
        row.serialize_entry(
            "realized_genotype_counts",
            &Dosages(summary.realized_genotype_counts()),
        )?;
        row.serialize_entry(
            "realized_genotype_frequencies",
            &Dosages(summary.realized_genotype_frequencies()),
        )?;
        row.serialize_entry("recount_matches", &summary.recount_matches())?;
        row.serialize_entry("recount_source", summary.recount_source())?;
        row.serialize_entry(
            "recounted_genotype_counts",
            &Dosages(summary.recounted_genotype_counts()),
        )?;
        row.serialize_entry(
            "recounted_genotype_frequencies",
            &Dosages(summary.recounted_genotype_frequencies()),
        )?;
        row.serialize_entry("requested_mode", summary.requested_mode().as_str())?;
        row.serialize_entry("site_id", summary.site_id())?;
        row.end()
    }
}

/// El manifiesto de un éxito, con su objeto `f5` ya construido **sin** los
/// resúmenes y estos emitidos por flujo en su posición alfabética.
struct SuccessManifest<'a> {
    execution_id: &'a str,
    f5: &'a Map<String, Value>,
    summaries: &'a [F5InitializationSummary],
}

impl Serialize for SuccessManifest<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut root = serializer.serialize_map(Some(5))?;
        root.serialize_entry("execution_id", self.execution_id)?;
        root.serialize_entry(
            "f5",
            &SuccessF5Object {
                map: self.f5,
                summaries: self.summaries,
            },
        )?;
        root.serialize_entry("manifest_kind", "gps_execution_manifest")?;
        root.serialize_entry("schema_version", F5_DOCUMENT_SCHEMA_VERSION)?;
        root.serialize_entry("state", "complete")?;
        root.end()
    }
}

struct SuccessF5Object<'a> {
    map: &'a Map<String, Value>,
    summaries: &'a [F5InitializationSummary],
}

impl Serialize for SuccessF5Object<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        const KEY: &str = "initialization_summaries";
        let mut object = serializer.serialize_map(Some(self.map.len() + 1))?;
        let mut emitted = false;
        for (key, value) in self.map {
            if !emitted && key.as_str() > KEY {
                object.serialize_entry(KEY, &Man70Rows(self.summaries))?;
                emitted = true;
            }
            object.serialize_entry(key, value)?;
        }
        if !emitted {
            object.serialize_entry(KEY, &Man70Rows(self.summaries))?;
        }
        object.end()
    }
}

/// Guarda previa: ningún resumen puede declarar que su recuento
/// independiente no coincide con el realizado.
///
/// El tipo hace hoy esta rama inalcanzable —el oráculo rechaza antes de
/// devolver un resumen desigual—, pero el recuento lo aporta una función
/// ajena y el contrato no puede descansar en que ninguna futura la relaje.
fn ensure_recounts_match(summaries: &[F5InitializationSummary]) -> Result<(), Diagnostic> {
    for summary in summaries {
        if !summary.recount_matches() {
            return Err(manifest_error(
                Path::new(F5_MANIFEST_PATH),
                format!(
                    "initialization_summaries: el recuento independiente de {}×{} no coincide con el realizado",
                    summary.population_id(),
                    summary.site_id()
                ),
            ));
        }
    }
    Ok(())
}

fn quota_ties_json(assignment: &F5FounderAssignment) -> Value {
    Value::Array(
        assignment
            .units()
            .iter()
            .filter_map(|unit| {
                unit.quota()
                    .and_then(|quota| quota.tie())
                    .map(|tie| (unit, tie))
            })
            .map(|(unit, tie)| {
                json!({
                    "population_id": unit.population_id().as_str(),
                    "block_id": unit.block_id().as_str(),
                    "remainder": tie.remainder(),
                    "candidate_ranks": tie.candidate_ranks(),
                    "select_count": tie.select_count(),
                    "selected_ranks": tie.selected_ranks(),
                })
            })
            .collect(),
    )
}

/// El objeto `f5` de un éxito, **sin** `initialization_summaries`: esos se
/// emiten por flujo al serializar, en su posición alfabética.
fn success_manifest_f5(
    request: &F5PublicationRequest<'_>,
    declaration: &F5PackageDeclaration,
    summaries: &[F5InitializationSummary],
    shared: &SharedContract,
) -> Result<Map<String, Value>, Diagnostic> {
    let assignment = request.assignment;
    let mut f5 = Map::new();
    // Todo lo que la declaración ya contiene se **deriva de ella**, no se
    // vuelve a leer del productor: el oráculo juzga esa declaración, y si el
    // manifiesto se construyera en paralelo desde las mismas estructuras el
    // oráculo estaría aprobando un documento distinto del publicado. Hoy
    // coincidirían, pero nada estructural lo garantizaría.
    f5.insert("schema_version".to_owned(), json!(F5_OBJECT_SCHEMA_VERSION));
    f5.insert("scenario_id".to_owned(), json!(declaration.scenario_id));
    f5.insert("seed".to_owned(), json!(declaration.seed));
    f5.insert("reference_id".to_owned(), json!(declaration.reference_id));
    f5.insert(
        "reference_sha256".to_owned(),
        json!(declaration.reference_sha256),
    );
    f5.insert(
        "eligible_bases".to_owned(),
        json!(declaration.bounds.eligible_bases),
    );
    f5.insert("ploidy".to_owned(), json!(declaration.ploidy));
    f5.insert(
        "founder_count".to_owned(),
        json!(*request.config.declared_scenario().founder_count.value()),
    );
    // `inbreeding_f_source` y el `block_index` y `source` de cada fila de
    // muestreo son los tres únicos términos que la declaración de F5.4.2 no
    // transporta —su contrato es cerrado y este hito no lo reabre—, de modo
    // que siguen viniendo del censo y de la resolución de muestreo. El resto
    // de cada objeto sale de la declaración.
    f5.insert(
        "subpopulations".to_owned(),
        Value::Array(
            declaration
                .subpopulations
                .iter()
                .zip(assignment.rosters())
                .map(|(declared, roster)| {
                    json!({
                        "population_id": declared.population_id,
                        "individual_count": declared.individual_count,
                        "inbreeding_f": optional_number(declared.inbreeding_f),
                        "inbreeding_f_source": roster.inbreeding_f_source().as_str(),
                    })
                })
                .collect(),
        ),
    );
    // Las dos claves del recuento haplotípico solo existen cuando el perfil
    // declara la distribución. Un escenario heredado publica el mismo objeto
    // `f5`, byte a byte, que antes de que esta decisión existiera.
    if let Some(distribution) = declaration.haplotype_count_distribution {
        f5.insert(
            "haplotype_count_distribution".to_owned(),
            json!({
                "algorithm": POISSON_COUNT_ALGORITHM,
                "lambda": distribution.lambda,
                "h_min": distribution.minimum,
                "h_max": distribution.maximum,
            }),
        );
        f5.insert(
            "haplotype_counts".to_owned(),
            Value::Array(
                declaration
                    .haplotype_counts
                    .iter()
                    .map(|row| {
                        json!({
                            "population_id": row.population_id,
                            "block_id": row.block_id,
                            "haplotype_count": row.haplotype_count,
                        })
                    })
                    .collect(),
            ),
        );
    }
    f5.insert(
        "attempt_summary".to_owned(),
        json!({
            "layout_attempts_used": declaration.attempt_summary.layout_attempts_used,
            "haplotype_attempts_total": declaration.attempt_summary.haplotype_attempts_total,
            "max_attempt_index": declaration.attempt_summary.max_attempt_index,
        }),
    );
    f5.insert(
        "sampling".to_owned(),
        sampling_json(declaration, assignment)?,
    );
    f5.insert("quota_ties".to_owned(), quota_ties_json(assignment));
    f5.insert(
        "f_realized".to_owned(),
        Value::Array(
            declaration
                .f_realized
                .iter()
                .map(|row| {
                    json!({
                        "population_id": row.population_id,
                        "block_id": row.block_id,
                        "target_f": optional_number(row.target_f),
                        "status": row.status,
                        "value": optional_number(row.value),
                    })
                })
                .collect(),
        ),
    );
    ensure_recounts_match(summaries)?;
    f5.insert(
        "truth".to_owned(),
        json!({
            "artifact_id": declaration.truth.artifact_id,
            "artifact_kind": declaration.truth.artifact_kind,
            "schema_version": declaration.truth.schema_version,
            "projection_version": declaration.truth.projection_version,
            "path": declaration.truth.path,
            "size_bytes": declaration.truth.size_bytes,
            "sha256": declaration.truth.sha256,
        }),
    );
    f5.insert(
        "artifacts".to_owned(),
        Value::Array(
            declaration
                .artifacts
                .iter()
                .map(|link| {
                    json!({
                        "artifact_id": link.artifact_id,
                        "artifact_kind": link.artifact_kind,
                        "schema_version": link.schema_version,
                        "path": link.path,
                        "size_bytes": link.size_bytes,
                        "sha256": link.sha256,
                    })
                })
                .collect(),
        ),
    );
    shared.insert_into(&mut f5);
    Ok(f5)
}

// ---------------------------------------------------------------------------
// Fallo controlado
// ---------------------------------------------------------------------------

/// Etapa en la que un fallo controlado F5 se produjo.
///
/// `Preflight` existe porque el esquema la enumera, pero **no tiene productor
/// en esta vía y no puede tenerlo**: DTC-089 exige que un rechazo de preflight
/// devuelva el diagnóstico DTC-044 y no abra ni publique destino, de modo que
/// nunca hay un manifiesto fallido que la declare. Se conserva para que un
/// lector del enumerado reconozca el literal del esquema, no para usarla.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F5FailureStage {
    Preflight,
    Layout,
    Haplotypes,
    Founders,
    Publication,
}

impl F5FailureStage {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Preflight => "preflight",
            Self::Layout => "layout",
            Self::Haplotypes => "haplotypes",
            Self::Founders => "founders",
            Self::Publication => "publication",
        }
    }

    /// Las cinco etapas exactas del contrato; cualquier otro literal es un
    /// documento que su propio esquema rechazaría.
    #[must_use]
    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "preflight" => Some(Self::Preflight),
            "layout" => Some(Self::Layout),
            "haplotypes" => Some(Self::Haplotypes),
            "founders" => Some(Self::Founders),
            "publication" => Some(Self::Publication),
            _ => None,
        }
    }
}

/// Catálogo cerrado de `f5-failure-v1.schema.json#/$defs/code`.
///
/// Existe porque el esquema es cerrado y el catálogo público de diagnósticos
/// no lo es: `GPS-INTERNAL-001`, `GPS-IO-001` o `GPS-FMT-006` son causas
/// perfectamente posibles de un fallo controlado —y las naturales de la etapa
/// `publication`— y ninguna pertenece a este enumerado. Publicar una de ellas
/// produciría, atómicamente y con destino final creado, un manifiesto que su
/// propio esquema rechaza. El tipo impide construirlo.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F5FailureCode {
    Configuration,
    Format,
    Domain101,
    Domain102,
    Domain103,
    Domain104,
    Domain105,
    Domain106,
    Resource,
    Internal,
}

impl F5FailureCode {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Configuration => "GPS-CONF-101",
            Self::Format => "GPS-FMT-101",
            Self::Domain101 => "GPS-DOMAIN-101",
            Self::Domain102 => "GPS-DOMAIN-102",
            Self::Domain103 => "GPS-DOMAIN-103",
            Self::Domain104 => "GPS-DOMAIN-104",
            Self::Domain105 => "GPS-DOMAIN-105",
            Self::Domain106 => "GPS-DOMAIN-106",
            Self::Resource => "GPS-RES-101",
            Self::Internal => "GPS-INTERNAL-101",
        }
    }

    #[must_use]
    pub fn parse(text: &str) -> Option<Self> {
        [
            Self::Configuration,
            Self::Format,
            Self::Domain101,
            Self::Domain102,
            Self::Domain103,
            Self::Domain104,
            Self::Domain105,
            Self::Domain106,
            Self::Resource,
            Self::Internal,
        ]
        .into_iter()
        .find(|code| code.as_str() == text)
    }
}

/// Reduce cualquier causa al enumerado cerrado del esquema de fallo.
///
/// El catálogo público de diagnósticos es más ancho que
/// `f5-failure-v1.schema.json#/$defs/code`. Una causa que no pertenezca se
/// publica como `GPS-INTERNAL-101` —que sí pertenece— y el diagnóstico
/// conserva **el código real**, anotando bajo `published_failure_code` cuál se
/// publicó en su lugar. Las tres vías de construcción usan esta función: sin
/// ella, dos de ellas degradaban el código en silencio.
fn publishable_code(diagnostic: Diagnostic) -> (F5FailureCode, Diagnostic) {
    let observed = diagnostic.code();
    match F5FailureCode::parse(observed) {
        Some(code) => (code, diagnostic),
        None => (
            F5FailureCode::Internal,
            diagnostic.with("published_failure_code", F5FailureCode::Internal.as_str()),
        ),
    }
}

/// Contadores exactos de `f5-failure.json`.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct F5FailureCounters {
    layout_attempts_used: u64,
    haplotype_attempts_used: u64,
    swap_proposals: u64,
    swap_accepted: u64,
}

/// Un fallo controlado ya reducido a los términos que el contrato publica.
///
/// No se construye desde texto diagnóstico: cada constructor recibe el testigo
/// tipado que su etapa produjo —`F5LayoutFailure` conserva el mejor layout y
/// `F5SearchFailure` el historial completo de intentos— y de ahí saca código,
/// etapa, identificadores y contadores.
#[derive(Clone, Debug, PartialEq)]
pub struct F5FailureReport {
    code: F5FailureCode,
    stage: F5FailureStage,
    scenario_id: String,
    population_id: Option<String>,
    block_id: Option<String>,
    counters: F5FailureCounters,
    best_layout: Option<F5BestLayout>,
    best_attempt: Option<F5AttemptRecord>,
    diagnostic: Diagnostic,
}

impl F5FailureReport {
    /// Fallo de layout. `layout_attempts_used` lo dice el propio fallo: no se
    /// deduce de la presencia de un mejor candidato, que existe desde el
    /// primer intento porque `consider_best` corre antes de materializar.
    #[must_use]
    pub fn from_layout(config: &EffectiveF5Config, failure: &F5LayoutFailure) -> Self {
        let (code, diagnostic) = publishable_code(failure.diagnostic().clone());
        Self {
            code,
            stage: F5FailureStage::Layout,
            scenario_id: config.scenario_id().to_owned(),
            population_id: failure.context().get("population_id").cloned(),
            block_id: failure.context().get("block_id").cloned(),
            counters: F5FailureCounters {
                layout_attempts_used: failure.attempts_used(),
                ..F5FailureCounters::default()
            },
            best_layout: failure.best_layout().cloned(),
            best_attempt: None,
            diagnostic,
        }
    }

    /// Fallo de la búsqueda haplotípica. Los contadores describen la unidad
    /// que el diagnóstico nombra: los tres agregados de otra unidad no
    /// explicarían por qué esta se agotó.
    #[must_use]
    pub fn from_search(
        config: &EffectiveF5Config,
        layout: &F5PhysicalLayout,
        failure: &F5SearchFailure,
    ) -> Self {
        let population_id = failure.context().get("population_id").cloned();
        let block_id = failure.context().get("block_id").cloned();
        let unit = |record: &F5AttemptRecord| {
            population_id
                .as_deref()
                .is_some_and(|population| record.population_id().as_str() == population)
                && block_id
                    .as_deref()
                    .is_some_and(|block| record.block_id().as_str() == block)
        };
        let mut counters = F5FailureCounters {
            layout_attempts_used: layout.successful_attempt_index().saturating_add(1),
            ..F5FailureCounters::default()
        };
        for record in failure.records().iter().filter(|record| unit(record)) {
            counters.haplotype_attempts_used = counters.haplotype_attempts_used.saturating_add(1);
            counters.swap_proposals = counters.swap_proposals.saturating_add(record.proposals());
            counters.swap_accepted = counters.swap_accepted.saturating_add(record.accepted());
        }
        let best_attempt = failure
            .records()
            .iter()
            .find(|record| unit(record) && record.selected_best())
            .cloned();
        let (code, diagnostic) = publishable_code(failure.diagnostic().clone());
        Self {
            code,
            stage: F5FailureStage::Haplotypes,
            scenario_id: config.scenario_id().to_owned(),
            population_id,
            block_id,
            counters,
            best_layout: None,
            best_attempt,
            diagnostic,
        }
    }

    /// Cualquier otro fallo controlado posterior a la admisión: cuota
    /// inviable, identidad o publicación. No hay mejor intento que declarar.
    #[must_use]
    pub fn from_diagnostic(
        config: &EffectiveF5Config,
        layout: Option<&F5PhysicalLayout>,
        stage: F5FailureStage,
        diagnostic: Diagnostic,
    ) -> Self {
        let (code, diagnostic) = publishable_code(diagnostic);
        Self {
            code,
            stage,
            scenario_id: config.scenario_id().to_owned(),
            population_id: diagnostic.context().get("population_id").cloned(),
            block_id: diagnostic.context().get("block_id").cloned(),
            counters: F5FailureCounters {
                layout_attempts_used: layout.map_or(0, |layout| {
                    layout.successful_attempt_index().saturating_add(1)
                }),
                ..F5FailureCounters::default()
            },
            best_layout: None,
            best_attempt: None,
            diagnostic,
        }
    }

    /// Diagnóstico público estable del fallo, para el canal DTC-044.
    #[must_use]
    pub const fn diagnostic(&self) -> &Diagnostic {
        &self.diagnostic
    }

    /// Código **publicable**: siempre uno del enumerado cerrado del esquema.
    /// Cuando la causa real no pertenece a él, [`Self::diagnostic`] la
    /// conserva íntegra.
    #[must_use]
    pub const fn code(&self) -> F5FailureCode {
        self.code
    }

    #[must_use]
    pub const fn stage(&self) -> F5FailureStage {
        self.stage
    }

    fn identifier(value: Option<&String>) -> Value {
        value.map_or(Value::Null, |value| json!(value))
    }

    fn counters_json(&self) -> Value {
        json!({
            "layout_attempts_used": self.counters.layout_attempts_used,
            "haplotype_attempts_used": self.counters.haplotype_attempts_used,
            "swap_proposals": self.counters.swap_proposals,
            "swap_accepted": self.counters.swap_accepted,
        })
    }

    fn best_layout_json(&self) -> Value {
        self.best_layout.as_ref().map_or(Value::Null, |best| {
            json!({
                "attempt_index": best.attempt_index(),
                "placed_blocks": best.placed_blocks(),
                "unplaced_blocks": best.unplaced_blocks(),
                "placed_bases": best.placed_bases(),
                "remaining_target_bases": best.remaining_target_bases(),
                "layout_bytes_sha256": best.layout_bytes_sha256(),
            })
        })
    }

    fn best_attempt_json(&self) -> Value {
        self.best_attempt.as_ref().map_or(Value::Null, |best| {
            json!({
                "attempt_index": best.attempt_index(),
                "duplicate_rows": best.duplicate_rows(),
                "bins_outside": best.bins_outside(),
                "signs_failed": best.signs_failed(),
                "max_r2_error": best.max_r2_error(),
                "sum_r2_error": best.sum_r2_error(),
                "max_marginal_error": best.max_marginal_error(),
                "matrix_bytes_sha256": best.matrix_bytes_sha256(),
            })
        })
    }

    /// Los nueve términos que `f5-failure.json` y el objeto `f5_failure` del
    /// manifiesto comparten literalmente, salvo su lista `artifacts`.
    fn common_json(&self, limits: &Value) -> Map<String, Value> {
        let mut object = Map::new();
        object.insert("code".to_owned(), json!(self.code.as_str()));
        object.insert("stage".to_owned(), json!(self.stage.as_str()));
        object.insert("scenario_id".to_owned(), json!(self.scenario_id));
        object.insert(
            "population_id".to_owned(),
            Self::identifier(self.population_id.as_ref()),
        );
        object.insert(
            "block_id".to_owned(),
            Self::identifier(self.block_id.as_ref()),
        );
        object.insert("limits".to_owned(), limits.clone());
        object.insert("counters".to_owned(), self.counters_json());
        object.insert("best_layout".to_owned(), self.best_layout_json());
        object.insert("best_attempt".to_owned(), self.best_attempt_json());
        object
    }

    fn diagnostic_document(&self, limits: &Value) -> Value {
        let mut object = Map::new();
        object.insert(
            "diagnostic_kind".to_owned(),
            json!(F5_FAILURE_DIAGNOSTIC_KIND),
        );
        object.insert(
            "schema_version".to_owned(),
            json!(F5_FAILURE_SCHEMA_VERSION),
        );
        object.extend(self.common_json(limits));
        object.insert("artifacts".to_owned(), Value::Array(Vec::new()));
        Value::Object(object)
    }

    fn manifest_document(&self, execution_id: &ExecutionId, limits: &Value, link: Value) -> Value {
        let mut failure = self.common_json(limits);
        failure.insert("artifacts".to_owned(), Value::Array(vec![link]));
        json!({
            "manifest_kind": "gps_execution_manifest",
            "schema_version": F5_DOCUMENT_SCHEMA_VERSION,
            "execution_id": execution_id.as_str(),
            "state": "failed",
            "f5": {
                "schema_version": F5_OBJECT_SCHEMA_VERSION,
                "f5_failure": Value::Object(failure),
            },
        })
    }
}

// ---------------------------------------------------------------------------
// Frontera de publicación
// ---------------------------------------------------------------------------

/// Todo lo que un éxito F5 necesita para publicarse, enlazado entre sí.
///
/// Las longitudes de secuencia ya no son una entrada del llamador: se cablean
/// desde la admisión sellada, que es el único punto donde el escaneo de
/// referencia las midió.
#[derive(Debug)]
pub struct F5PublicationRequest<'a> {
    config: &'a EffectiveF5Config,
    sealed: &'a SealedF5Admission,
    execution_id: &'a ExecutionId,
    layout: &'a F5PhysicalLayout,
    search: &'a F5HaplotypeSearch,
    assignment: &'a F5FounderAssignment,
    sequences: F5ReferenceSequences,
}

impl<'a> F5PublicationRequest<'a> {
    /// # Errors
    ///
    /// `GPS-CONF-101` si la admisión sellada no conserva longitudes de
    /// secuencia utilizables.
    pub fn new(
        config: &'a EffectiveF5Config,
        sealed: &'a SealedF5Admission,
        execution_id: &'a ExecutionId,
        layout: &'a F5PhysicalLayout,
        search: &'a F5HaplotypeSearch,
        assignment: &'a F5FounderAssignment,
    ) -> Result<Self, Diagnostic> {
        Ok(Self {
            config,
            sealed,
            execution_id,
            layout,
            search,
            assignment,
            sequences: F5ReferenceSequences::from_sealed(sealed)?,
        })
    }

    /// La configuración efectiva que esta publicación usa.
    ///
    /// Es la que el paquete publica y la que la admisión sellada admite. Se
    /// expone para que la clausura de F5.5.1 pueda congelar en su testigo la
    /// clase fundadora de cada subpoblación **en el momento de publicar**, en
    /// vez de leerla de un objeto que le entreguen después.
    #[must_use]
    pub const fn config(&self) -> &'a EffectiveF5Config {
        self.config
    }

    /// La admisión sellada que esta publicación consume.
    #[must_use]
    pub const fn sealed(&self) -> &'a SealedF5Admission {
        self.sealed
    }

    /// La declaración que el oráculo independiente contrasta contra los bytes
    /// y que el manifiesto copiará.
    fn declaration(&self, artifacts: &F5SuccessArtifacts) -> F5PackageDeclaration {
        let profile = self.config.profile();
        let scenario = self.config.declared_scenario();
        let subpopulations = self
            .assignment
            .rosters()
            .iter()
            .map(|roster| F5DeclaredSubpopulation {
                population_id: roster.population_id().as_str().to_owned(),
                individual_count: roster.individual_count(),
                inbreeding_f: roster.inbreeding_f(),
            })
            .collect::<Vec<_>>();
        let ld_bins = subpopulations
            .iter()
            .flat_map(|subpopulation| {
                profile
                    .ld_bins(&subpopulation.population_id)
                    .unwrap_or_default()
                    .iter()
                    .map(|bin| F5DeclaredLdBin {
                        population_id: subpopulation.population_id.clone(),
                        distance_start_bp: bin.distance_start_bp(),
                        distance_end_bp: bin.distance_end_bp(),
                        target_mean_r2: bin.target_mean_r2(),
                        max_abs_error: bin.max_abs_error(),
                    })
            })
            .collect();
        let truth = artifacts.truth();
        F5PackageDeclaration {
            scenario_id: self.config.scenario_id().to_owned(),
            execution_id: self.execution_id.as_str().to_owned(),
            reference_id: self.config.reference_id().to_owned(),
            reference_sha256: self.sealed.reference_sha256().to_owned(),
            seed: self.config.seed().to_owned(),
            ploidy: self.config.ploidy(),
            profile_id: self.sealed.profile_id().to_owned(),
            profile_metadata_sha256: self.sealed.profile_metadata_sha256().to_owned(),
            haplotype_weights: profile.haplotype_weights().to_vec(),
            haplotype_count_distribution: profile.haplotype_count_distribution().map(
                |distribution| F5DeclaredHaplotypeCountDistribution {
                    lambda: distribution.lambda(),
                    minimum: distribution.minimum(),
                    maximum: distribution.maximum(),
                },
            ),
            // El H de cada unidad se publica **explícito** aunque sea
            // recontable desde `f5-haplotypes.tsv`: así el oráculo compara su
            // propio sorteo contra una cifra declarada, en vez de contra un
            // recuento de filas que él mismo tendría que reconstruir. Sin
            // distribución declarada la lista va vacía y el documento no gana
            // ninguna clave.
            haplotype_counts: if profile.haplotype_count_distribution().is_some() {
                declared_haplotype_counts(self.search)
            } else {
                Vec::new()
            },
            snp_spectrum: profile
                .snp_spectrum()
                .iter()
                .map(|band| F5DeclaredSpectrumBand {
                    alt_af_min: band.minimum(),
                    alt_af_max: band.maximum(),
                    weight: band.weight(),
                })
                .collect(),
            subpopulations,
            sampling: self
                .assignment
                .sampling()
                .iter()
                .map(|row| F5DeclaredSampling {
                    population_id: row.population_id().as_str().to_owned(),
                    block_id: row.block_id().as_str().to_owned(),
                    requested_mode: declared_mode(row.requested_mode()),
                    effective_mode: declared_mode(row.effective_mode()),
                    controlled_max_abs_error_declared: row.controlled_max_abs_error_declared(),
                    controlled_max_abs_error_effective: row.controlled_max_abs_error_effective(),
                })
                .collect(),
            f_realized: self
                .assignment
                .f_realized()
                .iter()
                .map(|row| F5DeclaredRealizedF {
                    population_id: row.population_id().as_str().to_owned(),
                    block_id: row.block_id().as_str().to_owned(),
                    target_f: row.target_f(),
                    status: match row.status() {
                        F5RealizedFStatus::Evaluated => "evaluated".to_owned(),
                        F5RealizedFStatus::NotApplicable => "not_applicable".to_owned(),
                    },
                    value: row.value(),
                })
                .collect(),
            ld_bins,
            normalizations: profile.normalizations().to_vec(),
            artifacts: artifacts
                .artifacts()
                .iter()
                .map(|link| F5DeclaredLink {
                    artifact_id: link.artifact_id().to_owned(),
                    artifact_kind: link.artifact_kind().to_owned(),
                    schema_version: link.schema_version().to_owned(),
                    path: link.path().to_owned(),
                    size_bytes: link.size_bytes(),
                    sha256: link.sha256().to_owned(),
                })
                .collect(),
            truth: F5DeclaredTruthLink {
                artifact_id: truth.artifact_id().to_owned(),
                artifact_kind: truth.artifact_kind().to_owned(),
                schema_version: truth.schema_version().to_owned(),
                projection_version: truth.projection_version().to_owned(),
                path: truth.path().to_owned(),
                size_bytes: truth.size_bytes(),
                sha256: truth.sha256().to_owned(),
            },
            bounds: F5DeclaredBounds {
                block_count: self.sealed.block_count(),
                eligible_bases: self.sealed.eligible_bases(),
                site_bound_per_block: self.sealed.site_bound_per_block(),
                bin_count_total: self.sealed.bin_count_total(),
                max_haplotype_attempts: scenario
                    .max_haplotype_attempts
                    .as_ref()
                    .map_or(F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS, |value| *value.value()),
                max_swap_proposals: scenario
                    .max_swap_proposals
                    .as_ref()
                    .map_or(F5_DEFAULT_MAX_SWAP_PROPOSALS, |value| *value.value()),
                tsv_rows_bound: self.sealed.tsv_rows_bound(),
                json_rows_bound: self.sealed.json_rows_bound(),
                output_bytes: self.sealed.output_bytes(),
                output_limit_bytes: self.sealed.output_limit_bytes(),
            },
            id_catalog_id: F5_ID_CATALOG.to_owned(),
            rng_catalog_id: f5_rng_catalog_for(self.config).0.to_owned(),
            attempt_summary: F5DeclaredAttemptSummary {
                layout_attempts_used: self.layout.successful_attempt_index().saturating_add(1),
                haplotype_attempts_total: self.search.haplotype_attempts_total(),
                max_attempt_index: self.search.max_attempt_index(),
            },
        }
    }
}

/// H publicado de cada bloque×subpoblación, en el orden canónico
/// `(population_id, block_id)` que usa la lista `sampling`.
fn declared_haplotype_counts(search: &F5HaplotypeSearch) -> Vec<F5DeclaredHaplotypeCount> {
    let mut rows = search
        .blocks()
        .iter()
        .map(|block| F5DeclaredHaplotypeCount {
            population_id: block.population_id().as_str().to_owned(),
            block_id: block.block_id().as_str().to_owned(),
            haplotype_count: block.matrix().haplotype_count() as u64,
        })
        .collect::<Vec<_>>();
    rows.sort_by(|left, right| {
        (&left.population_id, &left.block_id).cmp(&(&right.population_id, &right.block_id))
    });
    rows
}

const fn declared_mode(mode: F5SamplingMode) -> F5DeclaredMode {
    match mode {
        F5SamplingMode::Probabilistic => F5DeclaredMode::Probabilistic,
        F5SamplingMode::Controlled => F5DeclaredMode::Controlled,
    }
}

/// Publica un éxito F5 completo: los siete TSV, `outputs/truth.jsonl`,
/// `effective-config.json`, `manifest.json` y el renombrado atómico.
///
/// El recuento MAN-70 **no se ejecuta aquí**: `recount` lo aporta el llamador y
/// recibe la raíz de carga útil ya escrita y la declaración que este módulo
/// construyó. Es la frontera que deja el pico del oráculo fuera de la fórmula
/// v9, que no lo presupuesta.
///
/// Nada se publica a medias: el destino final solo aparece cuando los dos
/// documentos raíz y los ocho ficheros han pasado la puerta física. Si
/// cualquier paso falla, el `ExecutionWorkspace` se destruye con su provisional
/// marcado `INCOMPLETE` y el destino no llega a existir.
///
/// # Errors
///
/// Cualquier diagnóstico de F5.4.1, el que devuelva `recount`, `GPS-FMT-006`
/// si un documento producido no satisface su contrato o no corresponde al
/// paquete físico, y `GPS-IO-003` si el destino ya existe.
pub fn publish_f5_success<R>(
    mut workspace: ExecutionWorkspace,
    request: &F5PublicationRequest<'_>,
    identities: &mut F5FounderIdentities,
    recount: R,
) -> Result<PathBuf, Box<F5PublicationFailure>>
where
    R: FnOnce(&Path, &F5PackageDeclaration) -> Result<Vec<F5InitializationSummary>, Diagnostic>,
{
    let attempt = |workspace: &mut ExecutionWorkspace| -> Result<(), Diagnostic> {
        publish_f5_success_into(workspace, request, identities, recount)
    };
    match attempt(&mut workspace) {
        Ok(()) => workspace.publish().map_err(|diagnostic| {
            Box::new(F5PublicationFailure {
                workspace: None,
                diagnostic,
            })
        }),
        Err(diagnostic) => Err(Box::new(F5PublicationFailure {
            workspace: Some(workspace),
            diagnostic,
        })),
    }
}

/// Fallo de [`publish_f5_success`] que **devuelve el provisional** al llamador.
///
/// DTC-045 distingue tres resultados y solo uno deja el provisional
/// `INCOMPLETE`: la **interrupción no controlable**. Un fallo controlado
/// posterior al comienzo —un recuento MAN-70 que rechaza, una identidad que no
/// cuadra, una puerta física que tumba el paquete— tiene que «publicar
/// atómicamente el directorio final como fallido». Con la firma anterior eso
/// era imposible: el `ExecutionWorkspace` entraba por valor y al devolver
/// `Err` ya estaba destruido, de modo que el llamador se quedaba sin con qué
/// publicar y el único rastro era un provisional marcado como si el proceso
/// hubiera muerto. Se devuelve el provisional **vivo** para que el llamador
/// pueda cerrarlo con [`publish_f5_failure`].
///
/// **Qué no puede publicar un fallo, y queda escrito.** `workspace` es `None`
/// en cuanto la publicación entra en [`ExecutionWorkspace::publish`], que toma
/// `self` por valor: no solo cuando el renombrado ya empezó —destino existente
/// o sistema de ficheros sin garantía no-clobber—, sino **también** cuando
/// falla cualquiera de sus comprobaciones previas: el recuento del catálogo,
/// el marcador `INCOMPLETE` y la **segunda pasada** de la puerta física. Ese
/// último caso es justamente un `GPS-FMT-006` de etapa `publication`, y no
/// poder publicarlo es defendible: significa que el provisional se alteró
/// entre dos revalidaciones, de modo que ya no se sabe qué se estaría
/// declarando. A esos casos se suma que publicar el fallo vuelva a fallar: un
/// `ENOSPC` real no deja escribir ni el diagnóstico. Todos **son** la
/// interrupción de DTC-045 y dejan el provisional `INCOMPLETE`, que es
/// exactamente lo que esa decisión reserva para ellos.
#[derive(Debug)]
pub struct F5PublicationFailure {
    workspace: Option<ExecutionWorkspace>,
    diagnostic: Diagnostic,
}

impl F5PublicationFailure {
    /// El provisional vivo, si el fallo admite publicar un paquete `failed`.
    #[must_use]
    pub fn into_workspace(self) -> Option<ExecutionWorkspace> {
        self.workspace
    }

    #[must_use]
    pub const fn diagnostic(&self) -> &Diagnostic {
        &self.diagnostic
    }

    #[must_use]
    pub fn into_diagnostic(self) -> Diagnostic {
        self.diagnostic
    }

    /// Descompone el fallo en sus dos mitades, que es lo que necesita un
    /// llamador que va a publicar el paquete `failed`.
    #[must_use]
    pub fn into_parts(self) -> (Option<ExecutionWorkspace>, Diagnostic) {
        (self.workspace, self.diagnostic)
    }
}

fn publish_f5_success_into<R>(
    workspace: &mut ExecutionWorkspace,
    request: &F5PublicationRequest<'_>,
    identities: &mut F5FounderIdentities,
    recount: R,
) -> Result<(), Diagnostic>
where
    R: FnOnce(&Path, &F5PackageDeclaration) -> Result<Vec<F5InitializationSummary>, Diagnostic>,
{
    let artifacts = {
        let success = F5SuccessRequest::new(
            request.config,
            request.sealed,
            request.execution_id,
            request.layout,
            request.search,
            request.assignment,
            &request.sequences,
        )?;
        write_f5_success_artifacts(workspace, &success, identities)?
    };
    let declaration = request.declaration(&artifacts);
    let summaries = recount(workspace.payload_root(), &declaration)?;
    let shared = SharedContract::build(request.config, request.sealed)?;
    let effective = effective_document(request, &declaration, &shared)?;
    let f5 = success_manifest_f5(request, &declaration, &summaries, &shared)?;
    // La puerta de F5.0 compara cada escalar de `resource_preflight` contra la
    // proyección sellada, no una copia contra la otra: dos copias iguales
    // entre sí pero distintas de la admisión también se rechazan. El
    // manifiesto entra como documento mínimo con **el mismo valor** que
    // `insert_into` acaba de copiar en él, porque el objeto completo ya no
    // existe como árbol: su término dominante se emite por flujo.
    request.sealed.ensure_contract_copies_match(
        Path::new(F5_EFFECTIVE_CONFIG_PATH),
        &effective,
        Path::new(F5_MANIFEST_PATH),
        &json!({"f5": {"resource_preflight": shared.resource_preflight.clone()}}),
    )?;
    // Cada árbol se suelta en cuanto existen sus bytes: la puerta física relee
    // el documento publicado y no necesita —ni debe usar— la copia en memoria
    // que lo produjo.
    let effective_bytes = json_document_bytes(&effective)?;
    drop(effective);
    workspace.write_root_document(F5_EFFECTIVE_CONFIG_PATH, &effective_bytes)?;
    drop(effective_bytes);
    workspace.write_f5_manifest(&SuccessManifest {
        execution_id: request.execution_id.as_str(),
        f5: &f5,
        summaries: &summaries,
    })?;
    // El documento se suelta **antes** de juzgarlo: la puerta relee el
    // fichero publicado y mantenerlo vivo sumaría su tamaño al pico.
    drop(f5);
    drop(summaries);
    workspace.close_f5_manifest()
}

/// Publica un fallo controlado: `outputs/f5-failure.json`, el manifiesto
/// `failed` que lo enlaza y **ningún artefacto biológico**.
///
/// El `ExecutionWorkspace` puede llegar con artefactos ya escritos: se retiran
/// antes de formar el manifiesto y la puerta física comprueba, sobre el árbol
/// que se va a renombrar, que no queda ninguno de los siete TSV ni
/// `outputs/truth.jsonl`.
///
/// # Errors
///
/// `GPS-FMT-006` si el diagnóstico producido no satisface su contrato, y los
/// diagnósticos de E/S y publicación de `ExecutionWorkspace`.
pub fn publish_f5_failure(
    mut workspace: ExecutionWorkspace,
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
    execution_id: &ExecutionId,
    report: &F5FailureReport,
) -> Result<PathBuf, Diagnostic> {
    workspace.discard_all_artifacts()?;
    let limits = limits_json(config, sealed)?;
    let diagnostic = report.diagnostic_document(&limits);
    let record = workspace.write_artifact(
        F5_FAILURE_ARTIFACT_ID,
        F5_FAILURE_ARTIFACT_KIND,
        Path::new(F5_FAILURE_PATH),
        &json_document_bytes(&diagnostic)?,
    )?;
    let link = json!({
        "artifact_id": F5_FAILURE_ARTIFACT_ID,
        "artifact_kind": F5_FAILURE_ARTIFACT_KIND,
        "schema_version": F5_FAILURE_SCHEMA_VERSION,
        "path": F5_FAILURE_PATH,
        "size_bytes": record.size_bytes,
        "sha256": record.digest,
    });
    let manifest = report.manifest_document(execution_id, &limits, link);
    workspace.write_f5_manifest(&manifest)?;
    drop(manifest);
    workspace.close_f5_manifest()?;
    workspace.publish()
}
// ---------------------------------------------------------------------------
// Puerta física del paquete v5
// ---------------------------------------------------------------------------

/// Documentos de la puerta: tipados, no árboles `serde_json::Value`.
///
/// La primera versión parseaba los dos documentos **enteros** a `Value` y
/// picaba 9,57 veces el pico sellado con `B=32` y 12,59 con `B=128`: el árbol
/// de `initialization_summaries` crece con `B*N_pop*S` y cuesta del orden de
/// siete veces su texto. La puerta no necesita ese árbol —solo mira dos campos
/// por resumen—, de modo que los tipos lo declaran y `serde` lo salta sin
/// asignar nada. Lo que sí se materializa es exactamente lo que el contrato
/// obliga a comparar entre efectivo y manifiesto, que crece con
/// `B*N_pop`.
///
/// `deny_unknown_fields` en cada objeto cerrado es además la puerta que el
/// esquema no tenía dentro del productor: una clave de más en `/f5` o en un
/// enlace se rechaza aquí y no solo en el validador externo.
mod gate {
    use crate::guard_trace::rejected;
    use serde::de::value::MapAccessDeserializer;
    use serde::de::{Error as _, IgnoredAny, MapAccess, Unexpected, Visitor};
    use serde::{Deserialize, Deserializer};
    use serde_json::Value;
    use std::ops::Deref;

    /// Un objeto JSON que se lee como tal y **solo** como tal.
    ///
    /// `derive(Deserialize)` acepta también una lista con los valores en el
    /// orden de los campos (`["gps_execution_manifest", ...]`), que el
    /// esquema `manifest-general-v5` rechaza (DT-033). Este envoltorio pide
    /// un mapa y delega en `T` con el mismo flujo, sin materializar nada:
    /// `deny_unknown_fields` y los mensajes de `T` siguen intactos. Usa
    /// `deserialize_any` y no `deserialize_map` porque este último rechaza la
    /// lista por su cuenta, sin pasar por `rejected` y, por tanto, sin que la
    /// puerta de cobertura vea el rechazo.
    pub(super) struct Strict<T>(pub(super) T);

    impl<T> Deref for Strict<T> {
        type Target = T;
        fn deref(&self) -> &T {
            &self.0
        }
    }

    impl<'de, T: Deserialize<'de>> Deserialize<'de> for Strict<T> {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            struct Object<T>(std::marker::PhantomData<T>);
            impl<'de, T: Deserialize<'de>> Visitor<'de> for Object<T> {
                type Value = Strict<T>;
                fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    formatter.write_str("a JSON object")
                }
                fn visit_map<A: MapAccess<'de>>(self, map: A) -> Result<Strict<T>, A::Error> {
                    T::deserialize(MapAccessDeserializer::new(map)).map(Strict)
                }
                fn visit_seq<A: serde::de::SeqAccess<'de>>(
                    self,
                    _: A,
                ) -> Result<Strict<T>, A::Error> {
                    Err(rejected(A::Error::invalid_type(
                        Unexpected::Seq,
                        &"an object, not a list of its values",
                    )))
                }
            }
            deserializer.deserialize_any(Object(std::marker::PhantomData))
        }
    }

    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct Root {
        pub(super) manifest_kind: String,
        pub(super) schema_version: String,
        pub(super) execution_id: String,
        pub(super) state: String,
        #[allow(dead_code)]
        f5: IgnoredAny,
    }

    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct Link {
        pub(super) artifact_id: String,
        pub(super) artifact_kind: String,
        pub(super) schema_version: String,
        pub(super) path: String,
        pub(super) size_bytes: u64,
        pub(super) sha256: String,
    }

    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct TruthLink {
        pub(super) artifact_id: String,
        pub(super) artifact_kind: String,
        pub(super) schema_version: String,
        pub(super) projection_version: String,
        pub(super) path: String,
        pub(super) size_bytes: u64,
        pub(super) sha256: String,
    }

    /// Los dos únicos campos que REQ-MAN-70 obliga a comprobar en cada
    /// resumen. El resto del objeto se salta sin materializarse.
    #[derive(Deserialize)]
    pub(super) struct Summary {
        pub(super) recount_source: String,
        pub(super) recount_matches: bool,
    }

    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct SuccessDocument {
        #[allow(dead_code)]
        manifest_kind: IgnoredAny,
        #[allow(dead_code)]
        schema_version: IgnoredAny,
        #[allow(dead_code)]
        execution_id: IgnoredAny,
        #[allow(dead_code)]
        state: IgnoredAny,
        pub(super) f5: Strict<SuccessF5>,
    }

    /// Las veinticinco claves exactas del objeto `f5` de un éxito, más las
    /// dos que solo existen cuando el perfil declara la distribución del
    /// recuento haplotípico.
    ///
    /// Las dos nuevas son `Option` y **no** pueden aparecer por separado: la
    /// puerta lo comprueba tras deserializar, porque `serde` no sabe expresar
    /// «ambas o ninguna».
    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct SuccessF5 {
        pub(super) schema_version: Value,
        pub(super) scenario_id: Value,
        pub(super) seed: String,
        pub(super) reference_id: Value,
        pub(super) reference_sha256: Value,
        pub(super) profile: Value,
        pub(super) mask: Value,
        pub(super) eligible_bases: u64,
        pub(super) target: Value,
        pub(super) ploidy: Value,
        pub(super) founder_count: Value,
        pub(super) subpopulations: Vec<Value>,
        pub(super) algorithms: Value,
        pub(super) rng_catalog: Value,
        pub(super) id_catalog: Value,
        pub(super) limits: Value,
        pub(super) resource_preflight: Value,
        pub(super) normalizations: Value,
        #[serde(default)]
        pub(super) haplotype_count_distribution: Option<Value>,
        #[serde(default)]
        pub(super) haplotype_counts: Option<Vec<Value>>,
        #[allow(dead_code)]
        attempt_summary: IgnoredAny,
        pub(super) sampling: Value,
        /// Declarado **no observable** desde los bytes publicados, igual que
        /// en F5.4.2: ningún TSV enumera los rangos empatados de una cuota
        /// controlada, de modo que ni la puerta ni el oráculo pueden
        /// contrastarlos. Queda dicho aquí para que no se confunda con un
        /// olvido.
        #[allow(dead_code)]
        quota_ties: IgnoredAny,
        #[allow(dead_code)]
        f_realized: IgnoredAny,
        pub(super) initialization_summaries: Vec<Strict<Summary>>,
        pub(super) truth: Strict<TruthLink>,
        pub(super) artifacts: Vec<Strict<Link>>,
    }

    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct FailedDocument {
        #[allow(dead_code)]
        manifest_kind: IgnoredAny,
        #[allow(dead_code)]
        schema_version: IgnoredAny,
        #[allow(dead_code)]
        execution_id: IgnoredAny,
        #[allow(dead_code)]
        state: IgnoredAny,
        pub(super) f5: Strict<FailedF5>,
    }

    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct FailedF5 {
        pub(super) schema_version: String,
        pub(super) f5_failure: Strict<Failure>,
    }

    /// Los nueve términos comunes que `f5-failure.json` y el objeto
    /// `f5_failure` del manifiesto tienen que publicar **iguales**.
    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct Failure {
        pub(super) code: String,
        pub(super) stage: String,
        pub(super) scenario_id: Value,
        pub(super) population_id: Value,
        pub(super) block_id: Value,
        pub(super) limits: Value,
        pub(super) counters: Value,
        pub(super) best_layout: Value,
        pub(super) best_attempt: Value,
        pub(super) artifacts: Vec<Strict<Link>>,
    }

    /// `outputs/f5-failure.json`, con sus doce claves exactas.
    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct FailureDiagnostic {
        pub(super) diagnostic_kind: String,
        pub(super) schema_version: String,
        pub(super) code: String,
        pub(super) stage: String,
        pub(super) scenario_id: Value,
        pub(super) population_id: Value,
        pub(super) block_id: Value,
        pub(super) limits: Value,
        pub(super) counters: Value,
        pub(super) best_layout: Value,
        pub(super) best_attempt: Value,
        pub(super) artifacts: Vec<IgnoredAny>,
    }

    /// Las dieciocho claves exactas del objeto `f5` del efectivo.
    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct EffectiveF5 {
        pub(super) schema_version: Value,
        pub(super) scenario_id: Value,
        pub(super) reference_id: Value,
        pub(super) reference_sha256: Value,
        pub(super) seed: String,
        pub(super) ploidy: Value,
        pub(super) founder_count: Value,
        pub(super) subpopulations: Vec<Value>,
        pub(super) profile: Value,
        pub(super) mask: Value,
        pub(super) target: Value,
        pub(super) sampling: Value,
        pub(super) normalizations: Value,
        pub(super) limits: Value,
        pub(super) resource_preflight: Value,
        pub(super) algorithms: Value,
        pub(super) rng_catalog: Value,
        pub(super) id_catalog: Value,
    }

    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    pub(super) struct EffectiveDocument {
        pub(super) effective_config_kind: String,
        pub(super) schema_version: String,
        #[allow(dead_code)]
        source_path: String,
        pub(super) seed: String,
        pub(super) f5: Strict<EffectiveF5>,
    }
}

/// Analiza el documento **desde su fichero**, sin materializar su texto.
///
/// El manifiesto v5 crece con `B*N_pop*S` y su texto, solo él, excede el pico
/// que la admisión sella. Leerlo por flujo es lo que permite que la puerta
/// quepa en ese presupuesto.
fn parse<T: serde::de::DeserializeOwned>(
    document_path: &Path,
    path: &Path,
    pointer: &str,
) -> Result<T, Diagnostic> {
    let file = File::open(document_path).map_err(|error| {
        manifest_error(
            path,
            format!(
                "{pointer}: no se pudo abrir {}: {error}",
                document_path.display()
            ),
        )
    })?;
    serde_json::from_reader(BufReader::new(file))
        .map_err(|error| manifest_error(path, format!("{pointer}: {error}")))
}

/// Reabre el fichero enlazado y recomputa su tamaño y SHA-256.
fn verify_link(
    payload_root: &Path,
    manifest_path: &Path,
    link: &gate::Link,
    pointer: &str,
    catalog: &mut Vec<(String, String)>,
) -> Result<(), Diagnostic> {
    if link.schema_version != F5_TSV_SCHEMA_VERSION {
        return Err(manifest_error(
            manifest_path,
            format!(
                "{pointer}/schema_version: el contrato fija {F5_TSV_SCHEMA_VERSION} y el enlace publica {}",
                link.schema_version
            ),
        ));
    }
    let identity = identify_file(&payload_root.join(&link.path), link.path.as_str())?;
    if identity.size_bytes() != link.size_bytes || identity.digest() != link.sha256 {
        return Err(manifest_error(
            manifest_path,
            format!(
                "{pointer}: el enlace declara {} bytes y {} y el fichero publicado mide {} con {}",
                link.size_bytes,
                link.sha256,
                identity.size_bytes(),
                identity.digest()
            ),
        ));
    }
    catalog.push((link.artifact_id.clone(), link.path.clone()));
    Ok(())
}

/// Enumera la carga útil y exige que sea **exactamente** el conjunto
/// declarado.
///
/// La puerta anterior comprobaba una lista literal de rutas prohibidas con
/// `exists()`, de modo que cualquier fichero que ningún manifiesto declarase
/// —un documento raíz escrito por el llamador, un residuo— llegaba al destino
/// final sin ruta, tamaño ni digest, en los **dos** estados. No es
/// hipotético: DTC-035 exige conservar la configuración original como copia
/// TOML exacta y en v5 esa copia solo puede ser un documento raíz, porque
/// `artifacts` es `prefixItems`+`items:false`. El primer consumidor del hueco
/// será F5.6.
fn verify_payload_inventory(
    payload_root: &Path,
    manifest_path: &Path,
    declared: &BTreeSet<String>,
) -> Result<(), Diagnostic> {
    let mut found = BTreeSet::new();
    collect_payload_files(payload_root, payload_root, manifest_path, &mut found)?;
    if &found == declared {
        return Ok(());
    }
    let undeclared = found.difference(declared).cloned().collect::<Vec<_>>();
    let missing = declared.difference(&found).cloned().collect::<Vec<_>>();
    Err(manifest_error(
        manifest_path,
        format!(
            "la carga útil no coincide con lo declarado: sobran {undeclared:?} y faltan {missing:?}"
        ),
    ))
}

fn collect_payload_files(
    payload_root: &Path,
    directory: &Path,
    manifest_path: &Path,
    found: &mut BTreeSet<String>,
) -> Result<(), Diagnostic> {
    let entries = std::fs::read_dir(directory).map_err(|error| {
        manifest_error(
            manifest_path,
            format!("no se pudo recorrer {}: {error}", directory.display()),
        )
    })?;
    for entry in entries {
        let entry = entry.map_err(|error| {
            manifest_error(
                manifest_path,
                format!(
                    "no se pudo leer una entrada de {}: {error}",
                    directory.display()
                ),
            )
        })?;
        let path = entry.path();
        // `file_type` no sigue enlaces, y `is_dir` sí. Un paquete que declara
        // cada fichero por su huella no puede admitir un enlace en su árbol:
        // el destino quedaría fuera de la carga útil y su contenido no sería
        // el que se hasheó. Hoy es inalcanzable —el provisional es 0o700, lo
        // crea el propio proceso y ninguna pieza F5 crea enlaces—, pero la
        // recursión tampoco tiene tope y un ciclo la colgaría.
        let kind = entry.file_type().map_err(|error| {
            manifest_error(
                manifest_path,
                format!("no se pudo tipar {}: {error}", path.display()),
            )
        })?;
        if kind.is_dir() {
            collect_payload_files(payload_root, &path, manifest_path, found)?;
            continue;
        }
        // La guarda no cambia el veredicto —el inventario ya vería el enlace
        // como fichero sin declarar— sino la **razón**: dice que es un enlace
        // en vez de limitarse a llamarlo no declarado. Lo que sí cambia el
        // veredicto es tipar con `file_type` en lugar de `is_dir`, que seguía
        // el enlace y recorría lo que hubiera detrás.
        if !kind.is_file() {
            return Err(manifest_error(
                manifest_path,
                format!(
                    "{} no es un fichero regular ni un directorio: un paquete declarado por huella no admite enlaces",
                    path.display()
                ),
            ));
        }
        let relative = path.strip_prefix(payload_root).map_err(|_| {
            manifest_error(
                manifest_path,
                format!("{} no cuelga de la carga útil", path.display()),
            )
        })?;
        let text = relative
            .to_str()
            .ok_or_else(|| {
                manifest_error(
                    manifest_path,
                    format!("{} no es una ruta UTF-8", relative.display()),
                )
            })?
            .replace('\\', "/");
        found.insert(text);
    }
    Ok(())
}

/// Juzga el paquete v5 completo desde los bytes que se van a publicar.
///
/// Se ejecuta dos veces: al cerrar el manifiesto y otra vez inmediatamente
/// antes del renombrado, porque entre ambos momentos el árbol provisional
/// todavía puede alterarse. Es pública porque un paquete **ya publicado** se
/// juzga con la misma puerta: no existe una segunda lectura más laxa para el
/// destino final.
///
/// # Errors
///
/// `GPS-FMT-006` en cuanto un término del manifiesto no corresponde al
/// contrato v5, a los bytes del árbol o a la copia efectiva; `GPS-IO-001` si
/// un fichero enlazado no puede reabrirse.
pub fn verify_published_f5_package(
    payload_root: &Path,
    manifest_path: &Path,
    artifacts: &[ArtifactRecord],
) -> Result<(), Diagnostic> {
    let root: gate::Strict<gate::Root> =
        parse(manifest_path, manifest_path, "la raíz del manifiesto")?;
    if root.manifest_kind != "gps_execution_manifest"
        || root.schema_version != F5_DOCUMENT_SCHEMA_VERSION
    {
        return Err(manifest_error(
            manifest_path,
            "clase o versión incompatible con manifest-general-v5",
        ));
    }
    if ExecutionId::parse(&root.execution_id, "manifest").is_err() {
        return Err(manifest_error(
            manifest_path,
            format!(
                "execution_id: {} no es una identidad tipada",
                root.execution_id
            ),
        ));
    }
    let mut catalog = Vec::new();
    match root.state.as_str() {
        "complete" => verify_success(payload_root, manifest_path, &mut catalog)?,
        "failed" => verify_failure(payload_root, manifest_path, &mut catalog)?,
        other => {
            return Err(manifest_error(
                manifest_path,
                format!("state: {other:?} no es complete ni failed"),
            ));
        }
    }
    let mut observed = artifacts
        .iter()
        .map(|record| (record.id.clone(), record.path.clone()))
        .collect::<Vec<_>>();
    observed.sort_unstable();
    catalog.sort_unstable();
    if observed != catalog {
        return Err(manifest_error(
            manifest_path,
            format!(
                "el catálogo del provisional contiene {observed:?} y el manifiesto enlaza {catalog:?}"
            ),
        ));
    }
    Ok(())
}

/// Los siete TSV, la verdad, el efectivo y su igualdad exacta con el
/// manifiesto.
fn verify_success(
    payload_root: &Path,
    manifest_path: &Path,
    catalog: &mut Vec<(String, String)>,
) -> Result<(), Diagnostic> {
    // `Strict` aquí es defensa redundante, inalcanzable desde las pruebas: la raíz,
    // que lee el mismo fichero envuelta, rechaza antes la lista.
    let document: gate::Strict<gate::SuccessDocument> = parse(manifest_path, manifest_path, "/f5")?;
    let f5 = document.0.f5.0;
    if f5.artifacts.len() != F5_TSV_CONTRACTS.len() {
        return Err(manifest_error(
            manifest_path,
            format!(
                "/f5/artifacts: el contrato exige {} enlaces y el manifiesto publica {}",
                F5_TSV_CONTRACTS.len(),
                f5.artifacts.len()
            ),
        ));
    }
    for (index, (link, contract)) in f5.artifacts.iter().zip(&F5_TSV_CONTRACTS).enumerate() {
        let pointer = format!("/f5/artifacts/{index}");
        if link.artifact_id != contract.artifact_id()
            || link.artifact_kind != contract.artifact_kind()
            || link.path != contract.path()
        {
            return Err(manifest_error(
                manifest_path,
                format!(
                    "{pointer}: el contrato fija {}/{}/{} en esta posición",
                    contract.artifact_id(),
                    contract.artifact_kind(),
                    contract.path()
                ),
            ));
        }
        verify_link(payload_root, manifest_path, link, &pointer, catalog)?;
    }
    if f5.truth.artifact_id != F5_TRUTH_ARTIFACT_ID
        || f5.truth.artifact_kind != F5_TRUTH_ARTIFACT_KIND
        || f5.truth.schema_version != F5_TRUTH_SCHEMA_VERSION
        || f5.truth.projection_version != F5_TRUTH_PROJECTION_VERSION
        || f5.truth.path != F5_TRUTH_PATH
    {
        return Err(manifest_error(
            manifest_path,
            "/f5/truth: identidad, clase, versión, versión de proyección o ruta fuera del contrato",
        ));
    }
    let identity = identify_file(&payload_root.join(F5_TRUTH_PATH), F5_TRUTH_PATH)?;
    if identity.size_bytes() != f5.truth.size_bytes || identity.digest() != f5.truth.sha256 {
        return Err(manifest_error(
            manifest_path,
            format!(
                "/f5/truth: el enlace declara {} bytes y {} y el fichero publicado mide {} con {}",
                f5.truth.size_bytes,
                f5.truth.sha256,
                identity.size_bytes(),
                identity.digest()
            ),
        ));
    }
    catalog.push((f5.truth.artifact_id.clone(), f5.truth.path.clone()));

    // REQ-MAN-70 obliga a enlazar la fuente del recuento y a exigir la
    // igualdad. Un resumen que no cite `truth.jsonl` v2, o que confiese no
    // coincidir, invalida la publicación entera.
    for (index, summary) in f5.initialization_summaries.iter().enumerate() {
        if summary.recount_source != RECOUNT_SOURCE || !summary.recount_matches {
            return Err(manifest_error(
                manifest_path,
                format!(
                    "/f5/initialization_summaries/{index}: recount_source debe ser {RECOUNT_SOURCE} y recount_matches true"
                ),
            ));
        }
    }

    // Las dos claves del recuento haplotípico van juntas o no van. Una
    // distribución sin recuentos dejaría al oráculo sin nada que comparar, y
    // unos recuentos sin distribución no serían rederivables desde ningún
    // ámbito. `serde` no sabe expresarlo, así que lo dice esta puerta.
    if f5.haplotype_count_distribution.is_some() != f5.haplotype_counts.is_some() {
        return Err(manifest_error(
            manifest_path,
            "/f5: haplotype_count_distribution y haplotype_counts se publican juntos o no se publican",
        ));
    }
    if let Some(counts) = &f5.haplotype_counts {
        let units = f5.sampling.as_array().map_or(0, Vec::len);
        if counts.len() != units {
            return Err(manifest_error(
                manifest_path,
                format!(
                    "/f5/haplotype_counts: {} filas frente a las {units} unidades de /f5/sampling",
                    counts.len()
                ),
            ));
        }
    }

    // `eligible_bases` no viaja al efectivo y el oráculo no puede recomputarlo
    // —no abre la referencia—, pero el propio manifiesto publica la cota
    // superior en `resource_preflight.reference_bases_processed`: las bases
    // elegibles son un subconjunto de las recorridas. La cota inferior
    // —las bases que el layout colocó— la comprueba el oráculo.
    let processed = f5
        .resource_preflight
        .get("reference_bases_processed")
        .and_then(Value::as_u64)
        .ok_or_else(|| {
            manifest_error(
                manifest_path,
                "/f5/resource_preflight/reference_bases_processed: falta o no es u64",
            )
        })?;
    if f5.eligible_bases > processed {
        return Err(manifest_error(
            manifest_path,
            format!(
                "/f5/eligible_bases: {} excede las {processed} bases que el escaneo recorrió",
                f5.eligible_bases
            ),
        ));
    }

    // El árbol publicado es exactamente el conjunto declarado: los siete TSV,
    // la verdad y los dos documentos raíz. Ni uno más.
    let mut declared = F5_TSV_CONTRACTS
        .iter()
        .map(|contract| contract.path().to_owned())
        .collect::<BTreeSet<_>>();
    declared.insert(F5_TRUTH_PATH.to_owned());
    declared.insert(F5_MANIFEST_PATH.to_owned());
    declared.insert(F5_EFFECTIVE_CONFIG_PATH.to_owned());
    verify_payload_inventory(payload_root, manifest_path, &declared)?;

    let effective_path = payload_root.join(F5_EFFECTIVE_CONFIG_PATH);
    let effective: gate::Strict<gate::EffectiveDocument> =
        parse(&effective_path, manifest_path, F5_EFFECTIVE_CONFIG_PATH)?;
    verify_effective_equality(manifest_path, &f5, &effective)
}

/// La igualdad exacta que DTC-089 exige entre efectivo y manifiesto, releída
/// de los **dos ficheros publicados** y no de la memoria que los produjo.
fn verify_effective_equality(
    manifest_path: &Path,
    manifest: &gate::SuccessF5,
    effective: &gate::EffectiveDocument,
) -> Result<(), Diagnostic> {
    if effective.effective_config_kind != "gps_effective_config"
        || effective.schema_version != F5_DOCUMENT_SCHEMA_VERSION
    {
        return Err(manifest_error(
            manifest_path,
            format!("{F5_EFFECTIVE_CONFIG_PATH}: clase o versión incompatible con v5"),
        ));
    }
    // Las tres copias de la semilla: la raíz del efectivo, su objeto `f5` y el
    // objeto `f5` del manifiesto.
    if effective.seed != effective.f5.seed || effective.seed != manifest.seed {
        return Err(manifest_error(
            manifest_path,
            format!(
                "seed: la raíz efectiva publica {}, su objeto f5 {} y el manifiesto {}",
                effective.seed, effective.f5.seed, manifest.seed
            ),
        ));
    }
    let effective_f5 = &effective.f5;
    for (key, left, right) in [
        (
            "schema_version",
            &effective_f5.schema_version,
            &manifest.schema_version,
        ),
        (
            "scenario_id",
            &effective_f5.scenario_id,
            &manifest.scenario_id,
        ),
        (
            "reference_id",
            &effective_f5.reference_id,
            &manifest.reference_id,
        ),
        (
            "reference_sha256",
            &effective_f5.reference_sha256,
            &manifest.reference_sha256,
        ),
        ("ploidy", &effective_f5.ploidy, &manifest.ploidy),
        (
            "founder_count",
            &effective_f5.founder_count,
            &manifest.founder_count,
        ),
        ("sampling", &effective_f5.sampling, &manifest.sampling),
        ("profile", &effective_f5.profile, &manifest.profile),
        ("mask", &effective_f5.mask, &manifest.mask),
        ("target", &effective_f5.target, &manifest.target),
        ("algorithms", &effective_f5.algorithms, &manifest.algorithms),
        (
            "rng_catalog",
            &effective_f5.rng_catalog,
            &manifest.rng_catalog,
        ),
        ("id_catalog", &effective_f5.id_catalog, &manifest.id_catalog),
        ("limits", &effective_f5.limits, &manifest.limits),
        (
            "resource_preflight",
            &effective_f5.resource_preflight,
            &manifest.resource_preflight,
        ),
        (
            "normalizations",
            &effective_f5.normalizations,
            &manifest.normalizations,
        ),
    ] {
        if left != right {
            return Err(manifest_error(
                manifest_path,
                format!("/f5/{key}: la copia del manifiesto difiere de la efectiva"),
            ));
        }
    }
    verify_subpopulation_equality(
        manifest_path,
        &effective_f5.subpopulations,
        &manifest.subpopulations,
    )
}

/// El manifiesto copia cuatro de los seis campos de cada subpoblación
/// efectiva, en el mismo orden y con los mismos valores.
///
/// Las dos listas se prestan, no se clonan.
fn verify_subpopulation_equality(
    manifest_path: &Path,
    effective: &[Value],
    manifest: &[Value],
) -> Result<(), Diagnostic> {
    const COPIED: [&str; 4] = [
        "population_id",
        "individual_count",
        "inbreeding_f",
        "inbreeding_f_source",
    ];
    if effective.len() != manifest.len() {
        return Err(manifest_error(
            manifest_path,
            format!(
                "/f5/subpopulations: el efectivo publica {} y el manifiesto {}",
                effective.len(),
                manifest.len()
            ),
        ));
    }
    for (index, (left, right)) in effective.iter().zip(manifest).enumerate() {
        let right = right.as_object().ok_or_else(|| {
            manifest_error(
                manifest_path,
                format!("/f5/subpopulations/{index}: no es un objeto"),
            )
        })?;
        if right.len() != COPIED.len() || COPIED.iter().any(|key| !right.contains_key(*key)) {
            let observed = right.keys().cloned().collect::<Vec<_>>().join(",");
            return Err(manifest_error(
                manifest_path,
                format!(
                    "/f5/subpopulations/{index}: el objeto cerrado exige {COPIED:?} y trae [{observed}]"
                ),
            ));
        }
        for key in COPIED {
            if left.get(key) != right.get(key) {
                return Err(manifest_error(
                    manifest_path,
                    format!(
                        "/f5/subpopulations/{index}/{key}: la copia del manifiesto difiere de la efectiva"
                    ),
                ));
            }
        }
    }
    Ok(())
}

/// El fallo publica el diagnóstico y **cero** artefactos biológicos.
fn verify_failure(
    payload_root: &Path,
    manifest_path: &Path,
    catalog: &mut Vec<(String, String)>,
) -> Result<(), Diagnostic> {
    // Igual que en el éxito: defensa redundante, inalcanzable desde las pruebas.
    let document: gate::Strict<gate::FailedDocument> =
        parse(manifest_path, manifest_path, "/f5/f5_failure")?;
    let f5 = document.0.f5.0;
    if f5.schema_version != F5_OBJECT_SCHEMA_VERSION {
        return Err(manifest_error(
            manifest_path,
            format!(
                "/f5/schema_version: el contrato fija {F5_OBJECT_SCHEMA_VERSION} y el manifiesto publica {}",
                f5.schema_version
            ),
        ));
    }
    let failure = f5.f5_failure.0;
    // El catálogo público del esquema es cerrado: un código ajeno produciría
    // un manifiesto que su propio esquema rechaza.
    if F5FailureCode::parse(&failure.code).is_none() {
        return Err(manifest_error(
            manifest_path,
            format!(
                "/f5/f5_failure/code: {} no pertenece al catálogo cerrado de f5-failure-v1",
                failure.code
            ),
        ));
    }
    if F5FailureStage::parse(&failure.stage).is_none() {
        return Err(manifest_error(
            manifest_path,
            format!(
                "/f5/f5_failure/stage: {} no pertenece a las cinco etapas del contrato",
                failure.stage
            ),
        ));
    }
    let [link] = failure.artifacts.as_slice() else {
        return Err(manifest_error(
            manifest_path,
            format!(
                "/f5/f5_failure/artifacts: el contrato exige exactamente un enlace y el manifiesto publica {}",
                failure.artifacts.len()
            ),
        ));
    };
    if link.artifact_id != F5_FAILURE_ARTIFACT_ID
        || link.artifact_kind != F5_FAILURE_ARTIFACT_KIND
        || link.path != F5_FAILURE_PATH
    {
        return Err(manifest_error(
            manifest_path,
            "/f5/f5_failure/artifacts/0: el contrato fija el enlace no biológico a f5-failure.json",
        ));
    }
    verify_link(
        payload_root,
        manifest_path,
        link,
        "/f5/f5_failure/artifacts/0",
        catalog,
    )?;
    // El árbol publicado es exactamente el manifiesto y el diagnóstico:
    // ningún artefacto biológico sobrevive, existiera o no como temporal, y
    // ningún fichero sin declarar se cuela. Enumerarlo es lo que distingue
    // esta comprobación de la lista literal de rutas prohibidas que tenía
    // antes, que no veía nada fuera de ella.
    let declared = BTreeSet::from([F5_MANIFEST_PATH.to_owned(), F5_FAILURE_PATH.to_owned()]);
    verify_payload_inventory(payload_root, manifest_path, &declared)?;
    verify_failure_diagnostic(payload_root, manifest_path, &failure)
}

/// El diagnóstico publicado y su copia en el manifiesto dicen **lo mismo**.
///
/// La puerta anterior leía `code`, `stage` y el enlace del manifiesto y no
/// abría nunca `outputs/f5-failure.json`: un paquete cuyo diagnóstico
/// contradijera a su propio manifiesto en causa, etapa, contadores y
/// artefactos pasaba, incluso rehaciendo el enlace con el tamaño y el SHA-256
/// nuevos. El bucle de los nueve términos estaba escrito en la batería; aquí
/// se promueve al productor, igual que F5.4.2 promovió el lector de F5.4.1.
fn verify_failure_diagnostic(
    payload_root: &Path,
    manifest_path: &Path,
    manifest: &gate::Failure,
) -> Result<(), Diagnostic> {
    let document_path = payload_root.join(F5_FAILURE_PATH);
    let document: gate::Strict<gate::FailureDiagnostic> =
        parse(&document_path, manifest_path, F5_FAILURE_PATH)?;
    if document.diagnostic_kind != F5_FAILURE_DIAGNOSTIC_KIND
        || document.schema_version != F5_FAILURE_SCHEMA_VERSION
    {
        return Err(manifest_error(
            manifest_path,
            format!("{F5_FAILURE_PATH}: clase o versión incompatible con f5-failure-v1"),
        ));
    }
    if !document.artifacts.is_empty() {
        return Err(manifest_error(
            manifest_path,
            format!(
                "{F5_FAILURE_PATH}/artifacts: el contrato exige la lista vacía y publica {}",
                document.artifacts.len()
            ),
        ));
    }
    if document.code != manifest.code {
        return Err(failure_mismatch(manifest_path, "code"));
    }
    if document.stage != manifest.stage {
        return Err(failure_mismatch(manifest_path, "stage"));
    }
    for (term, left, right) in [
        ("scenario_id", &document.scenario_id, &manifest.scenario_id),
        (
            "population_id",
            &document.population_id,
            &manifest.population_id,
        ),
        ("block_id", &document.block_id, &manifest.block_id),
        ("limits", &document.limits, &manifest.limits),
        ("counters", &document.counters, &manifest.counters),
        ("best_layout", &document.best_layout, &manifest.best_layout),
        (
            "best_attempt",
            &document.best_attempt,
            &manifest.best_attempt,
        ),
    ] {
        if left != right {
            return Err(failure_mismatch(manifest_path, term));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;

    static NONCE: AtomicU64 = AtomicU64::new(0);

    fn temporary(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gps-f5-publication-unit-{}-{}-{label}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn digest(bytes: &[u8]) -> String {
        format!("{:x}", Sha256::digest(bytes))
    }

    /// Objetos compartidos sintéticos: la puerta física no interpreta su
    /// contenido, solo exige que las dos copias sean el mismo valor. Eso
    /// permite ejercitar el renombrado sin montar una cadena F5 entera.
    fn synthetic_shared() -> Map<String, Value> {
        let mut shared = Map::new();
        for key in SHARED_OBJECTS {
            shared.insert(key.to_owned(), json!({"synthetic": key}));
        }
        // `resource_preflight` es el único cuyo contenido la puerta sí mira:
        // de él sale la cota superior de `eligible_bases`.
        shared.insert(
            "resource_preflight".to_owned(),
            json!({"reference_bases_processed": 1024_u64}),
        );
        shared
    }

    fn synthetic_documents(links: &[(String, u64, String)]) -> (Value, Value) {
        let seed = "a".repeat(64);
        let mut effective_f5 = synthetic_shared();
        let mut manifest_f5 = synthetic_shared();
        for f5 in [&mut effective_f5, &mut manifest_f5] {
            f5.insert("schema_version".to_owned(), json!(F5_OBJECT_SCHEMA_VERSION));
            f5.insert("scenario_id".to_owned(), json!("scenario"));
            f5.insert("reference_id".to_owned(), json!("reference"));
            f5.insert("reference_sha256".to_owned(), json!("b".repeat(64)));
            f5.insert("seed".to_owned(), json!(seed));
            f5.insert("ploidy".to_owned(), json!(2_u64));
            f5.insert("founder_count".to_owned(), json!(1_u64));
            f5.insert("sampling".to_owned(), json!([]));
        }
        // Las cinco claves que solo el manifiesto tiene. Sin ellas la puerta
        // rechaza el objeto por incompleto, que es exactamente lo que debe
        // hacer: `deny_unknown_fields` y la obligatoriedad de cada campo
        // cierran el objeto `/f5` dentro del productor y no solo en el
        // validador externo.
        manifest_f5.insert("eligible_bases".to_owned(), json!(0_u64));
        manifest_f5.insert(
            "attempt_summary".to_owned(),
            json!({
                "layout_attempts_used": 1,
                "haplotype_attempts_total": 1,
                "max_attempt_index": 0,
            }),
        );
        manifest_f5.insert("quota_ties".to_owned(), json!([]));
        manifest_f5.insert("f_realized".to_owned(), json!([]));
        effective_f5.insert(
            "subpopulations".to_owned(),
            json!([{
                "population_id": "P1",
                "individual_count": 1,
                "inbreeding_f": 0.0,
                "inbreeding_f_source": "population",
                "sampling": "probabilistic",
                "founders": [{"individual_id": "I1", "founder_key": null}],
            }]),
        );
        manifest_f5.insert(
            "subpopulations".to_owned(),
            json!([{
                "population_id": "P1",
                "individual_count": 1,
                "inbreeding_f": 0.0,
                "inbreeding_f_source": "population",
            }]),
        );
        manifest_f5.insert("initialization_summaries".to_owned(), json!([]));
        let mut artifacts = Vec::new();
        for (contract, (_, size, sha)) in F5_TSV_CONTRACTS.iter().zip(links) {
            artifacts.push(json!({
                "artifact_id": contract.artifact_id(),
                "artifact_kind": contract.artifact_kind(),
                "schema_version": "1",
                "path": contract.path(),
                "size_bytes": size,
                "sha256": sha,
            }));
        }
        manifest_f5.insert("artifacts".to_owned(), Value::Array(artifacts));
        let (_, truth_size, truth_sha) = links.last().expect("enlace de verdad");
        manifest_f5.insert(
            "truth".to_owned(),
            json!({
                "artifact_id": "truth-jsonl",
                "artifact_kind": "truth_jsonl",
                "schema_version": "2",
                "projection_version": F5_TRUTH_PROJECTION_VERSION,
                "path": F5_TRUTH_PATH,
                "size_bytes": truth_size,
                "sha256": truth_sha,
            }),
        );
        (
            json!({
                "effective_config_kind": "gps_effective_config",
                "schema_version": F5_DOCUMENT_SCHEMA_VERSION,
                "source_path": "general.toml",
                "seed": seed,
                "f5": Value::Object(effective_f5),
            }),
            json!({
                "manifest_kind": "gps_execution_manifest",
                "schema_version": F5_DOCUMENT_SCHEMA_VERSION,
                "execution_id": "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4",
                "state": "complete",
                "f5": Value::Object(manifest_f5),
            }),
        )
    }

    /// Escribe un paquete v5 sintético en un provisional y lo deja con su
    /// manifiesto ya cerrado, listo para publicar.
    fn staged_workspace(label: &str) -> (PathBuf, ExecutionWorkspace) {
        staged_workspace_with(label, false).expect("paquete sintético")
    }

    fn staged_workspace_with(
        label: &str,
        extra_artifact: bool,
    ) -> Result<(PathBuf, ExecutionWorkspace), Diagnostic> {
        let root = temporary(label);
        std::fs::create_dir_all(&root).expect("raíz");
        let destination = root.join("published");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("provisional");
        let mut links = Vec::new();
        for contract in &F5_TSV_CONTRACTS {
            let bytes = format!("{}\n", contract.header()).into_bytes();
            workspace
                .write_artifact(
                    contract.artifact_id(),
                    contract.artifact_kind(),
                    Path::new(contract.path()),
                    &bytes,
                )
                .expect("artefacto");
            links.push((
                contract.path().to_owned(),
                bytes.len() as u64,
                digest(&bytes),
            ));
        }
        let truth = b"{\"record_kind\":\"header\"}\n".to_vec();
        workspace
            .write_artifact(
                "truth-jsonl",
                "truth_jsonl",
                Path::new(F5_TRUTH_PATH),
                &truth,
            )
            .expect("verdad");
        links.push((F5_TRUTH_PATH.to_owned(), truth.len() as u64, digest(&truth)));
        if extra_artifact {
            workspace
                .write_artifact(
                    "f5-08-extra-v1",
                    "f5_extra",
                    Path::new("outputs/f5-extra.tsv"),
                    b"extra\n",
                )
                .expect("artefacto ajeno");
        }
        let (effective, manifest) = synthetic_documents(&links);
        workspace
            .write_root_document(
                F5_EFFECTIVE_CONFIG_PATH,
                &json_document_bytes(&effective).expect("efectivo"),
            )
            .expect("documento raíz");
        workspace.write_f5_manifest(&manifest)?;
        workspace.close_f5_manifest()?;
        Ok((root, workspace))
    }

    /// Un artefacto del provisional que el manifiesto no enlaza impide cerrar
    /// el catálogo: el paquete publicado no puede contener nada sin declarar.
    ///
    /// Quien lo ve primero es el inventario del árbol, que cita la **ruta**;
    /// la comparación con el catálogo, que cita el identificador, es la
    /// segunda red. Se comprueba la ruta porque es lo que la puerta dice.
    #[test]
    fn an_artifact_the_manifest_does_not_link_blocks_the_closing_gate() {
        let rejection = staged_workspace_with("extra-artifact", true)
            .expect_err("un artefacto ajeno debe tumbar el cierre");
        assert_eq!(rejection.code(), "GPS-FMT-006");
        assert!(
            rejection
                .context()
                .get("detail")
                .is_some_and(|detail| detail.contains("outputs/f5-extra.tsv")),
            "{:?}",
            rejection.context()
        );
    }

    /// La publicación vuelve a juzgar el paquete **inmediatamente antes** del
    /// renombrado: un artefacto alterado entre el cierre del manifiesto y el
    /// renombrado no llega al destino.
    ///
    /// La corrupción actualiza también el catálogo del provisional, de modo
    /// que la única puerta que puede verla es la que recomputa el SHA-256 de
    /// los bytes contra el enlace **declarado en el manifiesto**.
    #[test]
    fn publication_revalidates_the_f5_package_immediately_before_the_rename() {
        let (root, mut workspace) = staged_workspace("revalidate");
        let destination = workspace.destination().to_path_buf();
        workspace
            .replace_artifact_bytes_for_validation_test("f5-02-sites-v1", b"corrupto\n")
            .expect("corromper el artefacto ya catalogado");
        let rejection = workspace
            .publish()
            .expect_err("un artefacto alterado tras el cierre no puede publicarse");
        assert_eq!(rejection.code(), "GPS-FMT-006");
        assert!(
            rejection
                .context()
                .get("detail")
                .is_some_and(|detail| detail.contains("/f5/artifacts/1")),
            "{:?}",
            rejection.context()
        );
        assert!(!destination.exists());
        let _ = std::fs::remove_dir_all(&root);
    }

    /// Ninguna causa ajena al enumerado cerrado del esquema puede publicarse.
    ///
    /// El catálogo público de diagnósticos es más ancho que
    /// `f5-failure-v1.schema.json#/$defs/code`, y las causas naturales de la
    /// etapa `publication` —`GPS-INTERNAL-001`, `GPS-IO-001`,
    /// `GPS-FMT-006`— no pertenecen a él. Antes de la corrección se publicaba
    /// atómicamente, con destino final creado, un manifiesto que su propio
    /// esquema rechaza.
    #[test]
    fn a_cause_outside_the_closed_catalogue_is_published_as_the_internal_code() {
        assert_eq!(F5FailureCode::parse("GPS-INTERNAL-001"), None);
        assert_eq!(F5FailureCode::parse("GPS-IO-001"), None);
        assert_eq!(F5FailureCode::parse("GPS-FMT-006"), None);
        assert_eq!(
            F5FailureCode::parse("GPS-DOMAIN-104"),
            Some(F5FailureCode::Domain104)
        );
        // Los diez literales del esquema, y solo ellos.
        for code in [
            "GPS-CONF-101",
            "GPS-FMT-101",
            "GPS-DOMAIN-101",
            "GPS-DOMAIN-102",
            "GPS-DOMAIN-103",
            "GPS-DOMAIN-104",
            "GPS-DOMAIN-105",
            "GPS-DOMAIN-106",
            "GPS-RES-101",
            "GPS-INTERNAL-101",
        ] {
            assert_eq!(
                F5FailureCode::parse(code).map(F5FailureCode::as_str),
                Some(code)
            );
        }
        // Y las cinco etapas exactas.
        for stage in [
            "preflight",
            "layout",
            "haplotypes",
            "founders",
            "publication",
        ] {
            assert_eq!(
                F5FailureStage::parse(stage).map(F5FailureStage::as_str),
                Some(stage)
            );
        }
        assert_eq!(F5FailureStage::parse("publicación"), None);
    }

    /// El mismo paquete, intacto, sí se publica: la prueba anterior mide la
    /// corrupción, no una puerta que rechace siempre.
    #[test]
    fn an_intact_synthetic_package_reaches_its_destination() {
        let (root, workspace) = staged_workspace("intact");
        let destination = workspace.destination().to_path_buf();
        let published = workspace.publish().expect("paquete íntegro publicado");
        assert_eq!(published, destination);
        assert!(destination.join(F5_MANIFEST_PATH).is_file());
        assert!(destination.join(F5_EFFECTIVE_CONFIG_PATH).is_file());
        assert!(!destination.join("INCOMPLETE").exists());
        let _ = std::fs::remove_dir_all(&root);
    }
}
