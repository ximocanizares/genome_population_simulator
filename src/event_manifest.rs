//! Contrato y escritura del primer manifiesto de un evento técnico F3.
//!
//! La unidad no ensambla ni publica el árbol de una ejecución: recibe hechos
//! cerrados por sus productores. F3.5.3 será responsable de calcular estados,
//! integrarlo en el manifiesto general y publicar ambos documentos.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use serde_json::{Value, json};

use crate::application::EffectiveTechnicalEvent;
use crate::diagnostics::{Diagnostic, diagnostic_matches_catalog};
use crate::identity::{EventId, IndividualId, LibraryId, ReadGroupId, ReferenceId};
use crate::results::{ArtifactRecord, UnitState, json_document_bytes, validate_relative_link};
use crate::technical_identity::TechnicalIdentityPlan;

pub const EVENT_MANIFEST_SCHEMA_VERSION: &str = "1";
pub const EVENT_MANIFEST_KIND: &str = "gps_event_manifest";
pub const TECHNICAL_SEQUENCING_EVENT_TYPE: &str = "technical_sequencing";
/// Tecnología del único perfil técnico que F3 admite. El perfil conserva sus
/// parámetros efectivos; este identificador distingue la tecnología del evento.
pub const ILLUMINA_SHORT_READ_TECHNOLOGY_V1: &str = "illumina-short-read-v1";
const REFERENCE_FASTA_ARTIFACT_KIND: &str = "reference_fasta";
const FASTQ_ARTIFACT_KIND: &str = "fastq";

/// Referencia común usada por el evento. Se enlaza, no se copia dentro del
/// directorio de evento.
#[derive(Clone, Copy, Debug)]
pub struct EventManifestReference<'a> {
    pub reference_id: &'a ReferenceId,
    pub fasta: &'a ArtifactRecord,
}

/// Estado de una salida propia del evento. Las salidas no completas no tienen
/// ruta, tamaño ni digest porque nunca se catalogan como artefactos válidos.
#[derive(Clone, Copy, Debug)]
pub struct EventManifestArtifact<'a> {
    pub id: &'a str,
    pub kind: &'a str,
    pub state: UnitState,
    pub individual_id: &'a IndividualId,
    pub library_id: &'a LibraryId,
    pub artifact: Option<&'a ArtifactRecord>,
    pub diagnostic: Option<&'a Diagnostic>,
}

/// Hechos del primer evento técnico que el escritor registra sin inferirlos de
/// rutas, nombres de fichero o mensajes.
pub struct EventManifest<'a> {
    pub event: &'a EffectiveTechnicalEvent,
    pub state: UnitState,
    pub identities: &'a TechnicalIdentityPlan,
    pub reference: EventManifestReference<'a>,
    pub artifacts: &'a [EventManifestArtifact<'a>],
    pub diagnostics: &'a [Diagnostic],
}

/// Ruta canónica del manifiesto dentro de la raíz de una ejecución.
#[must_use]
pub fn event_manifest_relative_path(event_id: &EventId) -> PathBuf {
    PathBuf::from("events")
        .join(event_id.as_str())
        .join("manifest.json")
}

/// Forma el documento JSON y verifica que los hechos recibidos sean coherentes
/// con la identidad técnica ya derivada.
pub fn event_manifest_document(manifest: &EventManifest<'_>) -> Result<Value, Diagnostic> {
    validate_writer_facts(manifest)?;
    let selected_individuals = manifest
        .identities
        .libraries()
        .map(|library| {
            json!({
                "individual_id": library.individual_id().as_str(),
                "library_id": library.library_id().as_str(),
                "read_group_id": library.read_group_id().as_str(),
            })
        })
        .collect::<Vec<_>>();
    let mut artifacts = manifest.artifacts.to_vec();
    artifacts.sort_by_key(|artifact| {
        (
            artifact.individual_id.as_str(),
            artifact.library_id.as_str(),
            artifact.id,
        )
    });
    let artifacts = artifacts
        .into_iter()
        .map(|artifact| {
            json!({
                "id": artifact.id,
                "kind": artifact.kind,
                "state": artifact.state.as_str(),
                "individual_id": artifact.individual_id.as_str(),
                "library_id": artifact.library_id.as_str(),
                "artifact": artifact.artifact.map(ArtifactRecord::to_json),
                "diagnostic_code": artifact.diagnostic.map(Diagnostic::code),
            })
        })
        .collect::<Vec<_>>();
    let document = json!({
        "manifest_kind": EVENT_MANIFEST_KIND,
        "schema_version": EVENT_MANIFEST_SCHEMA_VERSION,
        "event_id": manifest.event.event_id().as_str(),
        "event_type": TECHNICAL_SEQUENCING_EVENT_TYPE,
        "technology": ILLUMINA_SHORT_READ_TECHNOLOGY_V1,
        "state": manifest.state.as_str(),
        "profile": {
            "id": manifest.event.profile(),
            "read_length": manifest.event.read_length(),
            "quality": {
                "acgt_phred": manifest.event.acgt_phred(),
                "n_phred": manifest.event.n_phred(),
            },
        },
        "technical_batch_id": manifest.identities.batch_id().as_str(),
        "selected_individuals": selected_individuals,
        "event_subseed": manifest.identities.event_subseed().to_hex(),
        "reference": {
            "reference_id": manifest.reference.reference_id.as_str(),
            "fasta": manifest.reference.fasta.to_json(),
        },
        "artifacts": artifacts,
        "diagnostics": canonical_diagnostics(manifest.diagnostics),
    });
    validate_event_manifest_bytes(
        &event_manifest_relative_path(manifest.event.event_id()),
        &json_document_bytes(&document)?,
    )?;
    Ok(document)
}

/// Serializa el manifiesto ya validado. F3.5.3 lo escribirá como artefacto
/// temporal del `ExecutionWorkspace` y lo finalizará con la política DTC-045.
pub fn event_manifest_bytes(manifest: &EventManifest<'_>) -> Result<Vec<u8>, Diagnostic> {
    json_document_bytes(&event_manifest_document(manifest)?)
}

/// Valida el esquema propio y las relaciones internas que JSON Schema no puede
/// expresar (identidades, orden estable y enlaces entre bibliotecas y FASTQ).
pub fn validate_event_manifest_bytes(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    let document: Value = serde_json::from_slice(bytes)
        .map_err(|error| event_error(path, format!("JSON no legible: {error}")))?;
    let root = required_object(&document, "raíz", path)?;
    validate_closed(root, ROOT_FIELDS, "raíz", path)?;
    require_text(root, "manifest_kind", path, EVENT_MANIFEST_KIND)?;
    require_text(root, "schema_version", path, EVENT_MANIFEST_SCHEMA_VERSION)?;
    let event_id = parse_id::<EventId>(required_string(root, "event_id", path)?, "event", path)?;
    require_text(root, "event_type", path, TECHNICAL_SEQUENCING_EVENT_TYPE)?;
    require_text(root, "technology", path, ILLUMINA_SHORT_READ_TECHNOLOGY_V1)?;
    let state = parse_state(required_string(root, "state", path)?, path)?;
    validate_profile(required_value(root, "profile", path)?, path)?;
    let _batch = parse_id::<crate::identity::BatchId>(
        required_string(root, "technical_batch_id", path)?,
        "batch",
        path,
    )?;
    let selected =
        validate_selected_individuals(required_value(root, "selected_individuals", path)?, path)?;
    validate_hex(
        required_string(root, "event_subseed", path)?,
        "event_subseed",
        path,
    )?;
    let reference = validate_reference(required_value(root, "reference", path)?, path)?;
    let diagnostics = validate_diagnostics(required_value(root, "diagnostics", path)?, path)?;
    let artifact_states = validate_artifacts(
        required_value(root, "artifacts", path)?,
        &selected,
        &reference,
        &diagnostics,
        path,
    )?;
    validate_event_state(state, &artifact_states, &diagnostics, path)?;
    let _ = event_id;
    Ok(())
}

const ROOT_FIELDS: &[&str] = &[
    "manifest_kind",
    "schema_version",
    "event_id",
    "event_type",
    "technology",
    "state",
    "profile",
    "technical_batch_id",
    "selected_individuals",
    "event_subseed",
    "reference",
    "artifacts",
    "diagnostics",
];

fn validate_writer_facts(manifest: &EventManifest<'_>) -> Result<(), Diagnostic> {
    let path = event_manifest_relative_path(manifest.event.event_id());
    if manifest.identities.event_id() != manifest.event.event_id()
        || manifest.identities.reference_id() != manifest.reference.reference_id
    {
        return Err(event_error(
            &path,
            "el plan técnico no coincide con el evento o referencia",
        ));
    }
    if manifest.reference.fasta.kind != REFERENCE_FASTA_ARTIFACT_KIND {
        return Err(event_error(
            &path,
            "el enlace de referencia no apunta a un FASTA",
        ));
    }
    let configured = manifest
        .event
        .individual_ids()
        .iter()
        .collect::<BTreeSet<_>>();
    let planned = manifest
        .identities
        .libraries()
        .map(|library| library.individual_id())
        .collect::<BTreeSet<_>>();
    if configured.len() != manifest.event.individual_ids().len() || configured != planned {
        return Err(event_error(
            &path,
            "individuos técnicos y bibliotecas no coinciden",
        ));
    }
    let diagnostics = manifest
        .diagnostics
        .iter()
        .map(Diagnostic::code)
        .collect::<BTreeSet<_>>();
    if diagnostics.len() != manifest.diagnostics.len() {
        return Err(event_error(&path, "diagnósticos repetidos"));
    }
    let mut artifact_libraries = BTreeSet::new();
    for artifact in manifest.artifacts {
        let library = manifest
            .identities
            .library(artifact.individual_id)
            .ok_or_else(|| {
                event_error(&path, "un artefacto refiere un individuo no seleccionado")
            })?;
        if library.library_id() != artifact.library_id {
            return Err(event_error(
                &path,
                "un artefacto refiere una biblioteca ajena",
            ));
        }
        if !artifact_libraries.insert(artifact.individual_id.as_str()) {
            return Err(event_error(
                &path,
                "una biblioteca tiene más de un resultado FASTQ",
            ));
        }
        match artifact.state {
            UnitState::Complete => {
                let record = artifact.artifact.ok_or_else(|| {
                    event_error(&path, "un artefacto completo no tiene identidad física")
                })?;
                if artifact.diagnostic.is_some()
                    || record.id != artifact.id
                    || record.kind != artifact.kind
                {
                    return Err(event_error(&path, "artefacto completo incoherente"));
                }
            }
            UnitState::Incomplete | UnitState::Failed => {
                let diagnostic = artifact
                    .diagnostic
                    .ok_or_else(|| event_error(&path, "artefacto no completo sin diagnóstico"))?;
                if artifact.artifact.is_some() || !diagnostics.contains(diagnostic.code()) {
                    return Err(event_error(&path, "artefacto no completo incoherente"));
                }
            }
        }
    }
    if artifact_libraries.len() != planned.len()
        || !planned
            .iter()
            .all(|individual| artifact_libraries.contains(individual.as_str()))
    {
        return Err(event_error(
            &path,
            "falta el resultado FASTQ de una biblioteca seleccionada",
        ));
    }
    Ok(())
}

fn validate_profile(value: &Value, path: &Path) -> Result<(), Diagnostic> {
    let profile = required_object(value, "profile", path)?;
    validate_closed(profile, &["id", "read_length", "quality"], "profile", path)?;
    require_text(
        profile,
        "id",
        path,
        crate::application::ILLUMINA_SE_MINIMAL_V1,
    )?;
    if profile.get("read_length").and_then(Value::as_u64)
        != Some(u64::from(
            crate::application::ILLUMINA_SE_MINIMAL_READ_LENGTH,
        ))
    {
        return Err(event_error(path, "longitud de lectura de perfil inválida"));
    }
    let quality = required_object(
        required_value(profile, "quality", path)?,
        "profile.quality",
        path,
    )?;
    validate_closed(quality, &["acgt_phred", "n_phred"], "profile.quality", path)?;
    if quality.get("acgt_phred").and_then(Value::as_u64)
        != Some(u64::from(
            crate::application::ILLUMINA_SE_MINIMAL_ACGT_PHRED,
        ))
        || quality.get("n_phred").and_then(Value::as_u64)
            != Some(u64::from(crate::application::ILLUMINA_SE_MINIMAL_N_PHRED))
    {
        return Err(event_error(path, "calidad de perfil inválida"));
    }
    Ok(())
}

fn validate_selected_individuals(
    value: &Value,
    path: &Path,
) -> Result<BTreeMap<String, String>, Diagnostic> {
    let values = required_array(value, "selected_individuals", path)?;
    if values.is_empty() {
        return Err(event_error(path, "el evento no selecciona individuos"));
    }
    let mut selected = BTreeMap::new();
    let mut previous = None::<String>;
    let mut read_groups = BTreeSet::new();
    for value in values {
        let item = required_object(value, "selected_individuals[]", path)?;
        validate_closed(
            item,
            &["individual_id", "library_id", "read_group_id"],
            "selected_individuals[]",
            path,
        )?;
        let individual = parse_id::<IndividualId>(
            required_string(item, "individual_id", path)?,
            "individual",
            path,
        )?;
        let library =
            parse_id::<LibraryId>(required_string(item, "library_id", path)?, "library", path)?;
        let read_group = parse_id::<ReadGroupId>(
            required_string(item, "read_group_id", path)?,
            "read_group",
            path,
        )?;
        if previous
            .as_deref()
            .is_some_and(|prior| prior >= individual.as_str())
            || selected
                .insert(individual.as_str().to_owned(), library.as_str().to_owned())
                .is_some()
            || !read_groups.insert(read_group.as_str().to_owned())
        {
            return Err(event_error(
                path,
                "individuos, bibliotecas o grupos no son únicos y ordenados",
            ));
        }
        previous = Some(individual.as_str().to_owned());
    }
    Ok(selected)
}

fn validate_reference(value: &Value, path: &Path) -> Result<ArtifactRecord, Diagnostic> {
    let reference = required_object(value, "reference", path)?;
    validate_closed(reference, &["reference_id", "fasta"], "reference", path)?;
    let _ = parse_id::<ReferenceId>(
        required_string(reference, "reference_id", path)?,
        "reference",
        path,
    )?;
    let fasta = validate_artifact(required_value(reference, "fasta", path)?, path)?;
    if fasta.kind != REFERENCE_FASTA_ARTIFACT_KIND {
        return Err(event_error(path, "la referencia no enlaza un FASTA"));
    }
    Ok(fasta)
}

fn validate_diagnostics(value: &Value, path: &Path) -> Result<BTreeSet<String>, Diagnostic> {
    let values = required_array(value, "diagnostics", path)?;
    let mut codes = BTreeSet::new();
    let mut previous = None::<String>;
    for value in values {
        let diagnostic = required_object(value, "diagnostics[]", path)?;
        validate_closed(diagnostic, &["code", "message"], "diagnostics[]", path)?;
        let code = required_string(diagnostic, "code", path)?;
        let message = required_string(diagnostic, "message", path)?;
        if previous.as_deref().is_some_and(|prior| prior >= code)
            || !diagnostic_matches_catalog(code, message)
            || !codes.insert(code.to_owned())
        {
            return Err(event_error(
                path,
                "diagnóstico desconocido, alterado o repetido",
            ));
        }
        previous = Some(code.to_owned());
    }
    Ok(codes)
}

fn validate_artifacts(
    value: &Value,
    selected: &BTreeMap<String, String>,
    reference: &ArtifactRecord,
    diagnostics: &BTreeSet<String>,
    path: &Path,
) -> Result<Vec<UnitState>, Diagnostic> {
    let values = required_array(value, "artifacts", path)?;
    if values.is_empty() {
        return Err(event_error(path, "el evento no tiene resultados FASTQ"));
    }
    let mut ids = BTreeSet::new();
    let mut paths = BTreeSet::new();
    let mut artifact_libraries = BTreeSet::new();
    let mut states = Vec::with_capacity(values.len());
    let mut previous = None::<(String, String, String)>;
    for value in values {
        let artifact = required_object(value, "artifacts[]", path)?;
        validate_closed(
            artifact,
            &[
                "id",
                "kind",
                "state",
                "individual_id",
                "library_id",
                "artifact",
                "diagnostic_code",
            ],
            "artifacts[]",
            path,
        )?;
        let id = required_string(artifact, "id", path)?;
        let kind = required_string(artifact, "kind", path)?;
        if kind != FASTQ_ARTIFACT_KIND {
            return Err(event_error(path, "F3 solo admite artefactos FASTQ"));
        }
        let _ = parse_id::<crate::identity::TemplateId>(id, "artifact", path)?;
        let state = parse_state(required_string(artifact, "state", path)?, path)?;
        let individual = required_string(artifact, "individual_id", path)?;
        let library = required_string(artifact, "library_id", path)?;
        let ordering = (individual.to_owned(), library.to_owned(), id.to_owned());
        if selected.get(individual).map(String::as_str) != Some(library)
            || !ids.insert(id.to_owned())
            || !artifact_libraries.insert(individual.to_owned())
            || previous.as_ref().is_some_and(|prior| prior >= &ordering)
        {
            return Err(event_error(
                path,
                "artefactos FASTQ no son únicos, completos por biblioteca u ordenados",
            ));
        }
        previous = Some(ordering);
        match state {
            UnitState::Complete => {
                let record =
                    validate_artifact(required_non_null(artifact, "artifact", path)?, path)?;
                if record.id != id
                    || record.kind != kind
                    || !paths.insert(record.path.clone())
                    || record.id == reference.id
                    || record.path == reference.path
                    || artifact.get("diagnostic_code") != Some(&Value::Null)
                {
                    return Err(event_error(path, "artefacto FASTQ completo incoherente"));
                }
            }
            UnitState::Incomplete | UnitState::Failed => {
                if artifact.get("artifact") != Some(&Value::Null) {
                    return Err(event_error(
                        path,
                        "artefacto no completo no puede tener identidad física",
                    ));
                }
                let code = required_string(artifact, "diagnostic_code", path)?;
                if !diagnostics.contains(code) {
                    return Err(event_error(path, "diagnóstico de artefacto ausente"));
                }
            }
        }
        states.push(state);
    }
    if artifact_libraries.len() != selected.len()
        || !selected
            .keys()
            .all(|individual| artifact_libraries.contains(individual))
    {
        return Err(event_error(
            path,
            "falta el resultado FASTQ de una biblioteca seleccionada",
        ));
    }
    Ok(states)
}

fn validate_event_state(
    state: UnitState,
    artifact_states: &[UnitState],
    diagnostics: &BTreeSet<String>,
    path: &Path,
) -> Result<(), Diagnostic> {
    match state {
        UnitState::Complete
            if !diagnostics.is_empty()
                || artifact_states
                    .iter()
                    .any(|artifact_state| *artifact_state != UnitState::Complete) =>
        {
            Err(event_error(
                path,
                "un evento completo exige FASTQ completos y ningún diagnóstico",
            ))
        }
        UnitState::Incomplete | UnitState::Failed
            if diagnostics.is_empty()
                || artifact_states
                    .iter()
                    .all(|artifact_state| *artifact_state == UnitState::Complete) =>
        {
            Err(event_error(
                path,
                "un evento no completo exige diagnóstico y un FASTQ no completo",
            ))
        }
        _ => Ok(()),
    }
}

fn canonical_diagnostics(diagnostics: &[Diagnostic]) -> Vec<Value> {
    let mut diagnostics = diagnostics.iter().collect::<Vec<_>>();
    diagnostics.sort_by_key(|diagnostic| diagnostic.code());
    diagnostics
        .into_iter()
        .map(|diagnostic| {
            json!({
                "code": diagnostic.code(), "message": diagnostic.message(),
            })
        })
        .collect()
}

fn validate_artifact(value: &Value, path: &Path) -> Result<ArtifactRecord, Diagnostic> {
    let artifact = required_object(value, "artifact", path)?;
    validate_closed(
        artifact,
        &["id", "kind", "path", "size_bytes", "checksum"],
        "artifact",
        path,
    )?;
    let id = required_string(artifact, "id", path)?;
    let kind = required_string(artifact, "kind", path)?;
    let artifact_path = required_string(artifact, "path", path)?;
    validate_relative_link(artifact_path)
        .map_err(|_| event_error(path, "ruta de artefacto no canónica"))?;
    let size_bytes = artifact
        .get("size_bytes")
        .and_then(Value::as_u64)
        .ok_or_else(|| event_error(path, "tamaño de artefacto inválido"))?;
    let checksum = required_object(
        required_value(artifact, "checksum", path)?,
        "checksum",
        path,
    )?;
    validate_closed(checksum, &["algorithm", "digest"], "checksum", path)?;
    require_text(checksum, "algorithm", path, "sha256")?;
    let digest = required_string(checksum, "digest", path)?;
    validate_hex(digest, "checksum.digest", path)?;
    Ok(ArtifactRecord {
        id: id.to_owned(),
        kind: kind.to_owned(),
        path: artifact_path.to_owned(),
        size_bytes,
        digest: digest.to_owned(),
    })
}

fn parse_state(value: &str, path: &Path) -> Result<UnitState, Diagnostic> {
    match value {
        "incomplete" => Ok(UnitState::Incomplete),
        "complete" => Ok(UnitState::Complete),
        "failed" => Ok(UnitState::Failed),
        _ => Err(event_error(path, "estado inválido")),
    }
}

fn parse_id<T>(value: &str, category: &str, path: &Path) -> Result<T, Diagnostic>
where
    T: ParseTypedId,
{
    T::parse_id(value).map_err(|_| event_error(path, format!("ID {category} inválido")))
}

trait ParseTypedId: Sized {
    fn parse_id(value: &str) -> Result<Self, Diagnostic>;
}
macro_rules! typed_id_parser { ($($type:ty),+ $(,)?) => { $(impl ParseTypedId for $type { fn parse_id(value: &str) -> Result<Self, Diagnostic> { Self::parse(value, "event_manifest") } })+ }; }
typed_id_parser!(
    EventId,
    IndividualId,
    LibraryId,
    ReadGroupId,
    ReferenceId,
    crate::identity::BatchId,
    crate::identity::TemplateId
);

fn validate_hex(value: &str, field: &str, path: &Path) -> Result<(), Diagnostic> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'a'..=b'f'))
    {
        return Err(event_error(
            path,
            format!("{field} debe ser SHA-256 hexadecimal minúsculo"),
        ));
    }
    Ok(())
}

fn required_object<'a>(
    value: &'a Value,
    context: &str,
    path: &Path,
) -> Result<&'a serde_json::Map<String, Value>, Diagnostic> {
    value
        .as_object()
        .ok_or_else(|| event_error(path, format!("{context} no es objeto")))
}
fn required_array<'a>(
    value: &'a Value,
    context: &str,
    path: &Path,
) -> Result<&'a [Value], Diagnostic> {
    value
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| event_error(path, format!("{context} no es lista")))
}
fn required_value<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a Value, Diagnostic> {
    object
        .get(field)
        .ok_or_else(|| event_error(path, format!("falta {field}")))
}
fn required_non_null<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a Value, Diagnostic> {
    let value = required_value(object, field, path)?;
    if value.is_object() {
        Ok(value)
    } else {
        Err(event_error(path, format!("{field} no puede ser nulo")))
    }
}
fn required_string<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    required_value(object, field, path)?
        .as_str()
        .filter(|value| !value.is_empty())
        .ok_or_else(|| event_error(path, format!("{field} debe ser texto no vacío")))
}
fn require_text(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
    expected: &str,
) -> Result<(), Diagnostic> {
    if required_string(object, field, path)? == expected {
        Ok(())
    } else {
        Err(event_error(path, format!("{field} incompatible")))
    }
}
fn validate_closed(
    object: &serde_json::Map<String, Value>,
    fields: &[&str],
    context: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if object.len() == fields.len() && object.keys().all(|field| fields.contains(&field.as_str())) {
        Ok(())
    } else {
        Err(event_error(path, format!("campos inválidos en {context}")))
    }
}
fn event_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path.display().to_string())
    .with("schema", "manifest-event-v1.schema.json")
    .with("detail", detail.into())
}
