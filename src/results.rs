//! Documentos JSON mínimos producidos por F1.

use std::collections::BTreeSet;
use std::fs::File;
use std::io::Write;
use std::path::{Component, Path};

use serde_json::{Value, json};

use crate::application::{EffectiveConfig, EffectiveReducedReference};
use crate::diagnostics::Diagnostic;
use crate::integrity::{ArtifactIdentity, identify_file};
use crate::reference::{Reference, SequenceClass};
use crate::reference_reduction::MaterializedPlan;

pub const RESULT_SCHEMA_VERSION: &str = "1";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum UnitState {
    Incomplete,
    Complete,
    Failed,
}

impl UnitState {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Incomplete => "incomplete",
            Self::Complete => "complete",
            Self::Failed => "failed",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArtifactRecord {
    pub id: String,
    pub kind: String,
    pub path: String,
    pub size_bytes: u64,
    pub digest: String,
    pub state: UnitState,
    pub validated: bool,
}

impl ArtifactRecord {
    pub fn from_identity(
        id: impl Into<String>,
        kind: impl Into<String>,
        identity: &ArtifactIdentity,
    ) -> Result<Self, Diagnostic> {
        let path = identity
            .relative_path()
            .to_string_lossy()
            .replace('\\', "/");
        validate_relative_link(&path)?;
        Ok(Self {
            id: id.into(),
            kind: kind.into(),
            path,
            size_bytes: identity.size_bytes(),
            digest: identity.digest().to_owned(),
            state: UnitState::Complete,
            validated: true,
        })
    }

    #[must_use]
    pub fn to_json(&self) -> Value {
        json!({
            "id": self.id,
            "kind": self.kind,
            "path": self.path,
            "size_bytes": self.size_bytes,
            "checksum": {"algorithm": "sha256", "digest": self.digest},
            "state": self.state.as_str(),
            "validated": self.validated,
        })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConfigurationInputRecord {
    pub kind: String,
    pub original_path: String,
    pub copy: ArtifactRecord,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UnitRecord {
    pub id: String,
    pub kind: String,
    pub state: UnitState,
    pub exportable: bool,
    pub diagnostic_code: Option<String>,
}

#[derive(Clone, Copy)]
pub struct ReferenceOutput<'a> {
    pub reference: &'a Reference,
    pub fasta: &'a ArtifactRecord,
}

#[derive(Clone, Copy)]
pub struct ReducedReferenceOutput<'a> {
    pub effective: &'a EffectiveReducedReference,
    pub reference: &'a Reference,
    pub fasta: &'a ArtifactRecord,
    pub map: &'a ArtifactRecord,
    pub plan: &'a MaterializedPlan,
    pub source_fasta: &'a ArtifactRecord,
    pub state: UnitState,
}

pub struct GeneralManifest<'a> {
    pub execution_id: &'a str,
    pub state: UnitState,
    pub effective: &'a EffectiveConfig,
    pub configuration_input: &'a ConfigurationInputRecord,
    pub effective_config: Option<&'a ArtifactRecord>,
    pub reference: Option<ReferenceOutput<'a>>,
    pub reduced_references: &'a [ReducedReferenceOutput<'a>],
    pub artifacts: &'a [ArtifactRecord],
    pub units: &'a [UnitRecord],
    pub diagnostics: &'a [Diagnostic],
}

pub fn effective_config_document_with_reference(
    config: &EffectiveConfig,
    output: ReferenceOutput<'_>,
) -> Result<Value, Diagnostic> {
    effective_config_document_with_reductions(config, output, &[])
}

pub fn effective_config_document_with_reductions(
    config: &EffectiveConfig,
    output: ReferenceOutput<'_>,
    reductions: &[ReducedReferenceOutput<'_>],
) -> Result<Value, Diagnostic> {
    validate_reference_output(output)?;
    let composition = config.composition().map(|composition| {
        let values = composition.values();
        json!({"A": values[0], "C": values[1], "G": values[2], "T": values[3]})
    });
    Ok(json!({
        "effective_config_kind": "gps_effective_config",
        "schema_version": RESULT_SCHEMA_VERSION,
        "source_path": config.source_path(),
        "seed": {
            "root": config.reproducibility().root_seed().to_hex(),
            "origin": config.reproducibility().seed_origin().as_str(),
        },
        "reference": reference_json(config, output, composition),
        "reduced_references": reductions.iter().map(reduction_json).collect::<Vec<_>>(),
    }))
}

#[must_use]
pub fn general_manifest_document(manifest: &GeneralManifest<'_>) -> Value {
    json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": RESULT_SCHEMA_VERSION,
        "execution_id": manifest.execution_id,
        "state": manifest.state.as_str(),
        "configuration_input": {
            "kind": manifest.configuration_input.kind,
            "original_path": manifest.configuration_input.original_path,
            "copy": manifest.configuration_input.copy.to_json(),
        },
        "effective_config": manifest.effective_config.map(ArtifactRecord::to_json),
        "references": manifest.reference.map(|reference| reference_link_json(reference)).into_iter()
            .chain(manifest.reduced_references.iter().filter(|reduction| reduction.state == UnitState::Complete).map(|reduction| reference_link_json(ReferenceOutput { reference: reduction.reference, fasta: reduction.fasta })))
            .collect::<Vec<_>>(),
        "reference_reductions": manifest.reduced_references.iter().map(reduction_json).collect::<Vec<_>>(),
        "artifacts": manifest.artifacts.iter().map(ArtifactRecord::to_json).collect::<Vec<_>>(),
        "units": manifest.units.iter().map(unit_json).collect::<Vec<_>>(),
        "diagnostics": manifest.diagnostics.iter().map(|diagnostic| json!({
            "code": diagnostic.code(),
            "message": diagnostic.message(),
        })).collect::<Vec<_>>(),
    })
}

pub fn write_json_document(
    path: &Path,
    relative_path: &Path,
    document: &Value,
) -> Result<ArtifactIdentity, Diagnostic> {
    let bytes = json_document_bytes(document)?;
    let mut file = File::create(path).map_err(|error| io_error(path, "create", error))?;
    file.write_all(&bytes)
        .map_err(|error| io_error(path, "write", error))?;
    file.flush()
        .map_err(|error| io_error(path, "flush", error))?;
    drop(file);
    identify_file(path, relative_path)
}

pub fn json_document_bytes(document: &Value) -> Result<Vec<u8>, Diagnostic> {
    let mut bytes = serde_json::to_vec_pretty(document).map_err(json_error)?;
    bytes.push(b'\n');
    Ok(bytes)
}

pub struct ValidatedGeneralManifest {
    bytes: Vec<u8>,
    artifacts: Vec<ArtifactRecord>,
    workspace_token: u64,
    workspace_snapshot: Vec<ArtifactRecord>,
}

impl std::fmt::Debug for ValidatedGeneralManifest {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ValidatedGeneralManifest")
            .field("size_bytes", &self.bytes.len())
            .field("artifact_count", &self.artifacts.len())
            .finish_non_exhaustive()
    }
}

impl ValidatedGeneralManifest {
    pub(crate) fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub(crate) fn artifacts(&self) -> &[ArtifactRecord] {
        &self.artifacts
    }

    pub(crate) const fn workspace_token(&self) -> u64 {
        self.workspace_token
    }

    pub(crate) fn workspace_snapshot(&self) -> &[ArtifactRecord] {
        &self.workspace_snapshot
    }
}

pub fn validate_general_manifest_bytes(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    parse_general_manifest_bytes(path, bytes).map(|_| ())
}

pub(crate) fn bind_general_manifest_to_workspace(
    path: &Path,
    bytes: &[u8],
    workspace_token: u64,
    workspace_snapshot: Vec<ArtifactRecord>,
) -> Result<ValidatedGeneralManifest, Diagnostic> {
    let (bytes, artifacts) = parse_general_manifest_bytes(path, bytes)?;
    Ok(ValidatedGeneralManifest {
        bytes,
        artifacts,
        workspace_token,
        workspace_snapshot,
    })
}

pub(crate) fn parse_general_manifest_bytes(
    path: &Path,
    bytes: &[u8],
) -> Result<(Vec<u8>, Vec<ArtifactRecord>), Diagnostic> {
    let document: Value = serde_json::from_slice(bytes)
        .map_err(|error| manifest_error(path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| manifest_error(path, "la raíz no es un objeto"))?;
    const FIELDS: [&str; 11] = [
        "manifest_kind",
        "schema_version",
        "execution_id",
        "state",
        "configuration_input",
        "effective_config",
        "references",
        "reference_reductions",
        "artifacts",
        "units",
        "diagnostics",
    ];
    validate_closed_object(object, &FIELDS, "manifest.json", path)?;
    if object.get("manifest_kind").and_then(Value::as_str) != Some("gps_execution_manifest")
        || object.get("schema_version").and_then(Value::as_str) != Some("1")
    {
        return Err(manifest_error(path, "clase o versión incompatible"));
    }
    required_non_empty_text(object, "execution_id", path)?;
    match object.get("state").and_then(Value::as_str) {
        Some("complete" | "failed") => {}
        _ => return Err(manifest_error(path, "estado de ejecución inválido")),
    }
    validate_configuration_input(required_object(object, "configuration_input", path)?, path)?;
    match object
        .get("effective_config")
        .ok_or_else(|| manifest_error(path, "falta effective_config"))?
    {
        Value::Null => {}
        value => {
            validate_artifact(value, path)?;
        }
    }
    validate_references(required_array(object, "references", path)?, path)?;
    validate_reference_reductions(required_array(object, "reference_reductions", path)?, path)?;
    let values = required_array(object, "artifacts", path)?;
    let mut ids = BTreeSet::new();
    let mut paths = BTreeSet::new();
    let mut records = Vec::with_capacity(values.len());
    for value in values {
        let record = validate_artifact(value, path)?;
        if !ids.insert(record.id.clone()) || !paths.insert(record.path.clone()) {
            return Err(manifest_error(path, "ID o ruta de artefacto repetidos"));
        }
        records.push(record);
    }
    validate_units(required_array(object, "units", path)?, path)?;
    validate_diagnostics(required_array(object, "diagnostics", path)?, path)?;
    Ok((bytes.to_vec(), records))
}

fn validate_configuration_input(
    object: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 3] = ["kind", "original_path", "copy"];
    validate_closed_object(object, &FIELDS, "configuration_input", path)?;
    if object.get("kind").and_then(Value::as_str) != Some("general") {
        return Err(manifest_error(path, "configuration_input.kind inválido"));
    }
    required_non_empty_text(object, "original_path", path)?;
    validate_artifact(
        object
            .get("copy")
            .ok_or_else(|| manifest_error(path, "falta copy"))?,
        path,
    )?;
    Ok(())
}

fn validate_references(values: &[Value], path: &Path) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 2] = ["id", "fasta"];
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "una referencia no es objeto"))?;
        validate_closed_object(object, &FIELDS, "references[]", path)?;
        required_non_empty_text(object, "id", path)?;
        validate_artifact(
            object
                .get("fasta")
                .ok_or_else(|| manifest_error(path, "falta fasta"))?,
            path,
        )?;
    }
    Ok(())
}

fn validate_reference_reductions(values: &[Value], path: &Path) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 10] = [
        "id",
        "source_reference_id",
        "region_size_bases",
        "target",
        "selection_subseed",
        "selection_order",
        "state",
        "fasta",
        "map",
        "source_fasta_id",
    ];
    const TARGET_FIELDS: [&str; 3] = ["mode", "declared", "bases"];
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "una reducción no es objeto"))?;
        validate_closed_object(object, &FIELDS, "reference_reductions[]", path)?;
        for field in ["id", "source_reference_id", "source_fasta_id"] {
            required_non_empty_text(object, field, path)?;
        }
        required_positive_u64(object, "region_size_bases", path)?;
        let target = required_object(object, "target", path)?;
        validate_closed_object(
            target,
            &TARGET_FIELDS,
            "reference_reductions[].target",
            path,
        )?;
        match target.get("mode").and_then(Value::as_str) {
            Some("target_bases" | "reduction_percent") => {}
            _ => return Err(manifest_error(path, "modo de reducción inválido")),
        }
        required_non_empty_text(target, "declared", path)?;
        required_positive_u64(target, "bases", path)?;
        let subseed = required_non_empty_text(object, "selection_subseed", path)?;
        validate_lower_hex_64(subseed, "selection_subseed", path)?;
        for order in required_array(object, "selection_order", path)? {
            if order.as_u64().is_none_or(|value| value == 0) {
                return Err(manifest_error(path, "selection_order inválido"));
            }
        }
        match object.get("state").and_then(Value::as_str) {
            Some("incomplete" | "complete" | "failed") => {}
            _ => return Err(manifest_error(path, "estado de reducción inválido")),
        }
        for field in ["fasta", "map"] {
            validate_artifact(
                object
                    .get(field)
                    .ok_or_else(|| manifest_error(path, format!("falta {field}")))?,
                path,
            )?;
        }
    }
    Ok(())
}

fn validate_artifact(value: &Value, path: &Path) -> Result<ArtifactRecord, Diagnostic> {
    const FIELDS: [&str; 7] = [
        "id",
        "kind",
        "path",
        "size_bytes",
        "checksum",
        "state",
        "validated",
    ];
    const CHECKSUM_FIELDS: [&str; 2] = ["algorithm", "digest"];
    let artifact = value
        .as_object()
        .ok_or_else(|| manifest_error(path, "un artefacto no es objeto"))?;
    validate_closed_object(artifact, &FIELDS, "artifact", path)?;
    let id = required_non_empty_text(artifact, "id", path)?;
    let kind = required_non_empty_text(artifact, "kind", path)?;
    let relative = required_non_empty_text(artifact, "path", path)?;
    validate_relative_link(relative)?;
    let checksum = required_object(artifact, "checksum", path)?;
    validate_closed_object(checksum, &CHECKSUM_FIELDS, "artifact.checksum", path)?;
    if checksum.get("algorithm").and_then(Value::as_str) != Some("sha256") {
        return Err(manifest_error(path, "checksum no usa sha256"));
    }
    let digest = required_non_empty_text(checksum, "digest", path)?;
    validate_lower_hex_64(digest, "digest", path)?;
    let state = match artifact.get("state").and_then(Value::as_str) {
        Some("complete") => UnitState::Complete,
        Some("failed") => UnitState::Failed,
        _ => return Err(manifest_error(path, "estado de artefacto inválido")),
    };
    let validated = artifact
        .get("validated")
        .and_then(Value::as_bool)
        .ok_or_else(|| manifest_error(path, "validated inválido"))?;
    if validated != (state == UnitState::Complete) {
        return Err(manifest_error(
            path,
            "state y validated del artefacto son incompatibles",
        ));
    }
    Ok(ArtifactRecord {
        id: id.to_owned(),
        kind: kind.to_owned(),
        path: relative.to_owned(),
        size_bytes: artifact
            .get("size_bytes")
            .and_then(Value::as_u64)
            .ok_or_else(|| manifest_error(path, "size_bytes inválido"))?,
        digest: digest.to_owned(),
        state,
        validated,
    })
}

fn validate_units(values: &[Value], path: &Path) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 5] = ["id", "kind", "state", "exportable", "diagnostic_code"];
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "una unidad no es objeto"))?;
        validate_closed_object(object, &FIELDS, "units[]", path)?;
        required_non_empty_text(object, "id", path)?;
        required_non_empty_text(object, "kind", path)?;
        match object.get("state").and_then(Value::as_str) {
            Some("incomplete" | "complete" | "failed") => {}
            _ => return Err(manifest_error(path, "estado de unidad inválido")),
        }
        required_bool(object, "exportable", path)?;
        match object.get("diagnostic_code") {
            Some(Value::Null) => {}
            Some(Value::String(code)) => validate_diagnostic_code(code, path)?,
            _ => return Err(manifest_error(path, "diagnostic_code inválido")),
        }
    }
    Ok(())
}

fn validate_diagnostics(values: &[Value], path: &Path) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 2] = ["code", "message"];
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "un diagnóstico no es objeto"))?;
        validate_closed_object(object, &FIELDS, "diagnostics[]", path)?;
        validate_diagnostic_code(required_non_empty_text(object, "code", path)?, path)?;
        required_non_empty_text(object, "message", path)?;
    }
    Ok(())
}

fn validate_closed_object(
    object: &serde_json::Map<String, Value>,
    fields: &[&str],
    context: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if object.keys().any(|key| !fields.contains(&key.as_str())) {
        return Err(manifest_error(
            path,
            format!("campo desconocido en {context}"),
        ));
    }
    Ok(())
}

fn required_object<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a serde_json::Map<String, Value>, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_object)
        .ok_or_else(|| manifest_error(path, format!("{field} no es objeto")))
}

fn required_array<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a [Value], Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .ok_or_else(|| manifest_error(path, format!("{field} no es lista")))
}

fn required_non_empty_text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    let value = required_text(object, field, path)?;
    if value.is_empty() {
        return Err(manifest_error(path, format!("{field} está vacío")));
    }
    Ok(value)
}

fn required_positive_u64(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<u64, Diagnostic> {
    match object.get(field).and_then(Value::as_u64) {
        Some(value) if value > 0 => Ok(value),
        _ => Err(manifest_error(path, format!("{field} inválido"))),
    }
}

fn required_bool(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<bool, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_bool)
        .ok_or_else(|| manifest_error(path, format!("{field} inválido")))
}

fn validate_lower_hex_64(value: &str, field: &str, path: &Path) -> Result<(), Diagnostic> {
    if value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Ok(());
    }
    Err(manifest_error(path, format!("{field} SHA-256 inválido")))
}

fn validate_diagnostic_code(code: &str, path: &Path) -> Result<(), Diagnostic> {
    let mut parts = code.split('-');
    let valid_category = matches!(parts.next(), Some("GPS"))
        && matches!(
            parts.next(),
            Some("CONF" | "INPUT" | "FMT" | "DOMAIN" | "RES" | "TOOL" | "IO" | "INTERNAL")
        );
    let valid_number = parts.next().is_some_and(|number| {
        number.len() == 3 && number.bytes().all(|byte| byte.is_ascii_digit())
    });
    if valid_category && valid_number && parts.next().is_none() {
        return Ok(());
    }
    Err(manifest_error(path, "código de diagnóstico inválido"))
}

fn required_text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| manifest_error(path, format!("falta {field}")))
}

fn reference_json(
    config: &EffectiveConfig,
    output: ReferenceOutput<'_>,
    composition: Option<Value>,
) -> Value {
    json!({
        "id": config.reference_id(),
        "source": match config.reference_source() {
            crate::application::EffectiveReferenceSource::Simulated => "simulated",
            crate::application::EffectiveReferenceSource::Supplied { .. } => "supplied",
        },
        "total_genome_size": config.total_genome_size(),
        "composition": composition,
        "sequences": output.reference.sequences().iter().map(|sequence| json!({
            "id": sequence.id().as_str(),
            "class": sequence_class(sequence.class()),
            "length": sequence.length(),
            "unknown_base_count": sequence.unknown_base_count(),
        })).collect::<Vec<_>>(),
        "fasta": output.fasta.to_json(),
    })
}

fn reference_link_json(output: ReferenceOutput<'_>) -> Value {
    json!({
        "id": output.reference.id().as_str(),
        "fasta": output.fasta.to_json(),
    })
}

fn reduction_json(output: &ReducedReferenceOutput<'_>) -> Value {
    json!({
        "id": output.effective.id(),
        "source_reference_id": output.effective.source_reference_id(),
        "region_size_bases": output.effective.region_size_bases(),
        "target": {
            "mode": output.effective.target().mode(),
            "declared": output.effective.target().declared_value(),
            "bases": output.effective.target_bases(),
        },
        "selection_subseed": output.effective.selection_subseed().to_hex(),
        "selection_order": output.plan.segments().iter().map(|segment| segment.selection_order()).collect::<Vec<_>>(),
        "state": output.state.as_str(),
        "fasta": output.fasta.to_json(),
        "map": output.map.to_json(),
        "source_fasta_id": output.source_fasta.id,
    })
}

fn validate_reference_output(output: ReferenceOutput<'_>) -> Result<(), Diagnostic> {
    if output.reference.authoritative_resource().is_none()
        || output
            .reference
            .authoritative_resource()
            .and_then(Path::to_str)
            != Some(&output.fasta.path)
    {
        return Err(Diagnostic::internal_invariant(
            "results",
            "reference artifact is not the authoritative resource",
        ));
    }
    Ok(())
}

fn sequence_class(class: SequenceClass) -> &'static str {
    match class {
        SequenceClass::Chromosome => "chromosome",
        SequenceClass::Scaffold => "scaffold",
    }
}

fn unit_json(unit: &UnitRecord) -> Value {
    json!({
        "id": unit.id,
        "kind": unit.kind,
        "state": unit.state.as_str(),
        "exportable": unit.exportable,
        "diagnostic_code": unit.diagnostic_code,
    })
}

fn validate_relative_link(path: &str) -> Result<(), Diagnostic> {
    let candidate = Path::new(path);
    if path.is_empty()
        || candidate.is_absolute()
        || path.contains("//")
        || path.contains('\\')
        || path.ends_with('/')
        || path.bytes().any(|byte| byte <= b'\x1f' || byte == b'\x7f')
        || path
            .split('/')
            .any(|component| matches!(component, "." | ".."))
        || candidate
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(Diagnostic::format(
            "GPS-FMT-005",
            "Un enlace de manifiesto no es una ruta relativa canónica.",
        )
        .with("field", "path")
        .with("value", path));
    }
    Ok(())
}

fn manifest_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-006",
        "El manifiesto general no cumple su contrato mínimo.",
    )
    .with("path", path.display().to_string())
    .with("detail", detail.into())
}

fn json_error(error: serde_json::Error) -> Diagnostic {
    Diagnostic::internal_invariant("results", format!("JSON serialization failed: {error}"))
}

fn io_error(path: &Path, operation: &str, error: std::io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}
