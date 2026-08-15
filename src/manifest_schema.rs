//! Puerta estructural runtime del manifiesto general.
//!
//! El productor valida exactamente el esquema JSON Schema Draft 2020-12 que
//! publica. La dependencia queda encapsulada aquí: ninguna estructura o error
//! del crate forma parte del contrato público del simulador.

use std::sync::OnceLock;

use jsonschema::{Registry, Validator};
use serde_json::Value;

const COMMON_SCHEMA: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/schemas/results-common-v1.schema.json"
));
const GENERAL_MANIFEST_SCHEMA: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/schemas/manifest-general-v1.schema.json"
));
const EFFECTIVE_CONFIG_SCHEMA: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/schemas/effective-config-v1.schema.json"
));
const EVENT_MANIFEST_SCHEMA: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/schemas/manifest-event-v1.schema.json"
));
const RESOLUTION_EVIDENCE_SCHEMA: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/schemas/resolution-evidence-v1.schema.json"
));
const SCHEMA_BASE: &str =
    "https://genome-population-simulator.invalid/schemas/manifest-general-v1.schema.json";
const COMMON_SCHEMA_URI: &str =
    "https://genome-population-simulator.invalid/schemas/results-common-v1.schema.json";

static GENERAL_MANIFEST_VALIDATOR: OnceLock<Result<Validator, String>> = OnceLock::new();
static EFFECTIVE_CONFIG_VALIDATOR: OnceLock<Result<Validator, String>> = OnceLock::new();
static EVENT_MANIFEST_VALIDATOR: OnceLock<Result<Validator, String>> = OnceLock::new();
static RESOLUTION_EVIDENCE_VALIDATOR: OnceLock<Result<Validator, String>> = OnceLock::new();

pub(crate) fn validate_general_manifest_structure(document: &Value) -> Result<(), String> {
    let validator = general_manifest_validator()?;
    validator
        .validate(document)
        .map_err(|error| format!("{}: {error}", error.instance_path()))
}

pub(crate) fn validate_effective_config_structure(document: &Value) -> Result<(), String> {
    validate_with_embedded_schema(
        &EFFECTIVE_CONFIG_VALIDATOR,
        EFFECTIVE_CONFIG_SCHEMA,
        "effective-config-v1.schema.json",
        document,
    )
}

pub(crate) fn validate_event_manifest_structure(document: &Value) -> Result<(), String> {
    validate_with_embedded_schema(
        &EVENT_MANIFEST_VALIDATOR,
        EVENT_MANIFEST_SCHEMA,
        "manifest-event-v1.schema.json",
        document,
    )
}

pub(crate) fn validate_resolution_evidence_structure(document: &Value) -> Result<(), String> {
    validate_with_embedded_schema(
        &RESOLUTION_EVIDENCE_VALIDATOR,
        RESOLUTION_EVIDENCE_SCHEMA,
        "resolution-evidence-v1.schema.json",
        document,
    )
}

fn validate_with_embedded_schema(
    cell: &'static OnceLock<Result<Validator, String>>,
    schema_text: &str,
    schema_name: &str,
    document: &Value,
) -> Result<(), String> {
    let validator = cell
        .get_or_init(|| {
            let schema: Value = serde_json::from_str(schema_text)
                .map_err(|error| format!("{schema_name} no es JSON: {error}"))?;
            let common: Value = serde_json::from_str(COMMON_SCHEMA)
                .map_err(|error| format!("el esquema común publicado no es JSON: {error}"))?;
            jsonschema::draft202012::meta::validate(&schema)
                .map_err(|error| format!("{schema_name} no es Draft 2020-12: {error}"))?;
            let registry = Registry::new()
                .add(COMMON_SCHEMA_URI, &common)
                .map_err(|error| format!("no se pudo registrar el esquema común: {error}"))?
                .prepare()
                .map_err(|error| format!("no se pudo preparar el esquema común: {error}"))?;
            jsonschema::draft202012::options()
                .with_registry(&registry)
                .with_base_uri(format!(
                    "https://genome-population-simulator.invalid/schemas/{schema_name}"
                ))
                .build(&schema)
                .map_err(|error| format!("no se pudo compilar {schema_name}: {error}"))
        })
        .as_ref()
        .map_err(Clone::clone)?;
    validator
        .validate(document)
        .map_err(|error| format!("{}: {error}", error.instance_path()))
}

fn general_manifest_validator() -> Result<&'static Validator, String> {
    GENERAL_MANIFEST_VALIDATOR
        .get_or_init(|| {
            let schema: Value = serde_json::from_str(GENERAL_MANIFEST_SCHEMA)
                .map_err(|error| format!("el esquema publicado no es JSON: {error}"))?;
            let common: Value = serde_json::from_str(COMMON_SCHEMA)
                .map_err(|error| format!("el esquema común publicado no es JSON: {error}"))?;
            jsonschema::draft202012::meta::validate(&schema)
                .map_err(|error| format!("el esquema publicado no es Draft 2020-12: {error}"))?;
            jsonschema::draft202012::meta::validate(&common).map_err(|error| {
                format!("el esquema común publicado no es Draft 2020-12: {error}")
            })?;
            let registry = Registry::new()
                .add(COMMON_SCHEMA_URI, &common)
                .map_err(|error| format!("no se pudo registrar el esquema común: {error}"))?
                .prepare()
                .map_err(|error| format!("no se pudo preparar el esquema común: {error}"))?;
            jsonschema::draft202012::options()
                .with_registry(&registry)
                .with_base_uri(SCHEMA_BASE)
                .build(&schema)
                .map_err(|error| format!("no se pudo compilar el esquema publicado: {error}"))
        })
        .as_ref()
        .map_err(Clone::clone)
}

#[cfg(test)]
mod tests {
    use super::{
        COMMON_SCHEMA, GENERAL_MANIFEST_SCHEMA, RESOLUTION_EVIDENCE_SCHEMA,
        general_manifest_validator,
    };

    #[test]
    fn published_schema_is_valid_draft_2020_12_and_compiles() {
        let schema: serde_json::Value = serde_json::from_str(GENERAL_MANIFEST_SCHEMA).unwrap();
        let common: serde_json::Value = serde_json::from_str(COMMON_SCHEMA).unwrap();
        let resolution: serde_json::Value =
            serde_json::from_str(RESOLUTION_EVIDENCE_SCHEMA).unwrap();
        jsonschema::draft202012::meta::validate(&schema).unwrap();
        jsonschema::draft202012::meta::validate(&common).unwrap();
        jsonschema::draft202012::meta::validate(&resolution).unwrap();
        general_manifest_validator().unwrap();
    }
}
