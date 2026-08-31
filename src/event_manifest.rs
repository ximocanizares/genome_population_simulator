//! Contrato y escritura del primer manifiesto de un evento técnico F3.
//!
//! La unidad no ensambla ni publica el árbol de una ejecución: recibe hechos
//! cerrados por sus productores. F3.5.3 será responsable de calcular estados,
//! integrarlo en el manifiesto general y publicar ambos documentos.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use serde_json::{Value, json};

use crate::application::EffectiveTechnicalEvent;
use crate::biological_contamination::{
    BiologicalContaminationManifestRecipient, contamination_target_pair_count,
};
use crate::diagnostics::{Diagnostic, diagnostic_matches_catalog};
use crate::fastq_output::paired_fastq_relative_paths;
use crate::fastq_validation::F4CoverageReport;
use crate::identity::{EventId, IndividualId, LibraryId, ReadGroupId, ReferenceId};
use crate::results::{ArtifactRecord, UnitState, json_document_bytes, validate_relative_link};
use crate::technical_identity::TechnicalIdentityPlan;
use crate::technical_identity::TechnicalMate;

/// DTC-082 activa esta versión junto con effective-config v4 y manifest v4.
pub const EVENT_MANIFEST_SCHEMA_VERSION: &str = "2";
pub const EVENT_MANIFEST_KIND: &str = "gps_event_manifest";
pub const TECHNICAL_SEQUENCING_EVENT_TYPE: &str = "technical_sequencing";
/// Tecnología del único perfil técnico que F3 admite. El perfil conserva sus
/// parámetros efectivos; este identificador distingue la tecnología del evento.
pub const ILLUMINA_SHORT_READ_TECHNOLOGY_V1: &str = "illumina-short-read-pe-v1";
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
    pub mate: TechnicalMate,
    pub artifact: Option<&'a ArtifactRecord>,
    pub diagnostic: Option<&'a Diagnostic>,
}

/// Resumen de cobertura final de una biblioteca, sin procedencia por lectura.
#[derive(Clone, Copy, Debug)]
pub struct EventManifestCoverage<'a> {
    pub individual_id: &'a IndividualId,
    pub library_id: &'a LibraryId,
    pub report: F4CoverageReport,
}

/// Hechos del primer evento técnico que el escritor registra sin inferirlos de
/// rutas, nombres de fichero o mensajes.
pub struct EventManifest<'a> {
    pub event: &'a EffectiveTechnicalEvent,
    pub state: UnitState,
    pub identities: &'a TechnicalIdentityPlan,
    pub reference: EventManifestReference<'a>,
    pub artifacts: &'a [EventManifestArtifact<'a>],
    pub coverage: &'a [EventManifestCoverage<'a>],
    /// Hecho agregado de F4.3. Nunca contiene origen, ID ni coordenadas de
    /// una pareja; esas piezas solo existen durante la transformación.
    pub chimera_rate: f64,
    pub realized_chimera_pairs: u64,
    /// Vista agregada del plan y la materialización F4.4. No contiene claves
    /// de pareja, IDs de instancia ni datos materializables individuales.
    pub biological_contamination: &'a [BiologicalContaminationManifestRecipient],
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
            let configured = manifest
                .event
                .libraries()
                .iter()
                .find(|configured| configured.individual_id() == library.individual_id())
                .expect("effective event invariants require one UDI per library");
            json!({
                "individual_id": library.individual_id().as_str(),
                "library_id": library.library_id().as_str(),
                "read_group_id": library.read_group_id().as_str(),
                "i7": configured.i7(),
                "i5": configured.i5(),
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
                "mate": artifact.mate.number(),
                "artifact": artifact.artifact.map(ArtifactRecord::to_json),
                "diagnostic_code": artifact.diagnostic.map(Diagnostic::code),
            })
        })
        .collect::<Vec<_>>();
    let mut coverage = manifest.coverage.to_vec();
    coverage.sort_by_key(|entry| (entry.individual_id.as_str(), entry.library_id.as_str()));
    let coverage = coverage
        .into_iter()
        .map(|entry| {
            json!({
                "individual_id": entry.individual_id.as_str(),
                "library_id": entry.library_id.as_str(),
                "target_x": entry.report.target_x(),
                "basic_x": entry.report.basic_wgs_x(),
                "final_x": {
                    "own_x": entry.report.final_fastq_own_x(),
                    "contaminant_x": entry.report.final_fastq_contaminant_x(),
                },
                "denominator_bases": entry.report.denominator_bases(),
                "n_bases": entry.report.n_bases(),
            })
        })
        .collect::<Vec<_>>();
    let document = json!({
        "manifest_kind": EVENT_MANIFEST_KIND,
        "schema_version": EVENT_MANIFEST_SCHEMA_VERSION,
        "event_id": manifest.event.event_id().as_str(),
        "event_type": TECHNICAL_SEQUENCING_EVENT_TYPE,
        "technology": ILLUMINA_SHORT_READ_TECHNOLOGY_V1,
        "rng_catalog": crate::randomness::F4_TECHNICAL_RNG_CATALOG,
        "id_catalog": f4_id_catalog(),
        "state": manifest.state.as_str(),
        "profile": {
            "id": manifest.event.profile(),
            "resource_path": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_RESOURCE,
            "sha256": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256,
            "preparation": "pcr_free",
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
        "coverage": coverage,
        "artifacts": {
            "chimera": {
                "rate": manifest.chimera_rate,
                "realized_pairs": manifest.realized_chimera_pairs,
            },
            "biological_contamination": manifest.biological_contamination.iter().map(|entry| json!({
                "recipient_id": entry.recipient_id().as_str(),
                "recipient_library_id": entry.recipient_library_id().as_str(),
                "percent_range": {
                    "minimum": entry.minimum_percent(),
                    "maximum": entry.maximum_percent(),
                },
                "sampled_percent": entry.sampled_percent(),
                "clean_pair_count": entry.clean_pair_count(),
                "target_pair_count": entry.target_pair_count(),
                "removed_own_pair_count": entry.removed_own_pair_count(),
                "inserted_contaminant_pair_count": entry.inserted_contaminant_pair_count(),
                "retained_own_pair_count": entry.retained_own_pair_count(),
                "final_pair_count": entry.final_pair_count(),
                "subseeds": {
                    "percent": entry.percent_subseed().to_hex(),
                    "donors": entry.donor_subseed().to_hex(),
                    "removal": entry.removal_subseed().to_hex(),
                },
                "donor_contributions": entry.donor_contributions().iter().map(|(donor, count)| json!({
                    "individual_id": donor.as_str(),
                    "pairs": count,
                })).collect::<Vec<_>>(),
            })).collect::<Vec<_>>(),
            "pcr": {"applicable": manifest.event.pcr_duplicates_applicable()},
            "optical_duplicates": {"realized_pairs": manifest.event.realized_optical_duplicate_pairs()},
            "gc": {"mode": "neutral"},
            "adapters": {"enabled": false},
        },
        "outputs": artifacts,
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
    require_text(
        root,
        "rng_catalog",
        path,
        crate::randomness::F4_TECHNICAL_RNG_CATALOG,
    )?;
    validate_f4_id_catalog(required_value(root, "id_catalog", path)?, path)?;
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
    let coverage = validate_coverage(required_value(root, "coverage", path)?, &selected, path)?;
    let diagnostics = validate_diagnostics(required_value(root, "diagnostics", path)?, path)?;
    validate_technical_artifacts(required_value(root, "artifacts", path)?, path)?;
    let artifact_states = validate_outputs(
        required_value(root, "outputs", path)?,
        &event_id,
        &selected,
        &reference,
        &diagnostics,
        path,
    )?;
    validate_event_state(state, &artifact_states, &diagnostics, path)?;
    if state == UnitState::Complete && coverage.len() != selected.len() {
        return Err(event_error(
            path,
            "falta cobertura de una biblioteca completa",
        ));
    }
    if state != UnitState::Complete && !coverage.is_empty() {
        return Err(event_error(
            path,
            "un evento no completo no publica cobertura F4",
        ));
    }
    Ok(())
}

const ROOT_FIELDS: &[&str] = &[
    "manifest_kind",
    "schema_version",
    "event_id",
    "event_type",
    "technology",
    "rng_catalog",
    "id_catalog",
    "state",
    "profile",
    "technical_batch_id",
    "selected_individuals",
    "event_subseed",
    "reference",
    "coverage",
    "artifacts",
    "outputs",
    "diagnostics",
];

fn validate_writer_facts(manifest: &EventManifest<'_>) -> Result<(), Diagnostic> {
    let path = event_manifest_relative_path(manifest.event.event_id());
    if manifest.chimera_rate != manifest.event.chimera_rate()
        || !manifest.chimera_rate.is_finite()
        || !(0.0..=1.0).contains(&manifest.chimera_rate)
        || (manifest.chimera_rate == 0.0 && manifest.realized_chimera_pairs != 0)
    {
        return Err(event_error(
            &path,
            "la tasa de quimera no coincide con la configuración efectiva",
        ));
    }
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
    let mut artifact_mates = BTreeSet::new();
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
        if !artifact_mates.insert((artifact.individual_id.as_str(), artifact.mate)) {
            return Err(event_error(&path, "una biblioteca repite un mate FASTQ"));
        }
        match artifact.state {
            UnitState::Complete => {
                let record = artifact.artifact.ok_or_else(|| {
                    event_error(&path, "un artefacto completo no tiene identidad física")
                })?;
                if artifact.diagnostic.is_some()
                    || record.id != artifact.id
                    || record.kind != artifact.kind
                    || record.path
                        != paired_path(
                            manifest.event.event_id().as_str(),
                            artifact.individual_id.as_str(),
                            artifact.mate,
                        )
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
    if artifact_mates.len() != planned.len() * 2
        || !planned.iter().all(|individual| {
            TechnicalMate::all()
                .into_iter()
                .all(|mate| artifact_mates.contains(&(individual.as_str(), mate)))
        })
    {
        return Err(event_error(
            &path,
            "falta el resultado FASTQ de una biblioteca seleccionada",
        ));
    }
    if manifest.state == UnitState::Complete {
        let coverage_libraries = manifest
            .coverage
            .iter()
            .map(|entry| (entry.individual_id.as_str(), entry.library_id.as_str()))
            .collect::<BTreeSet<_>>();
        let planned_libraries = manifest
            .identities
            .libraries()
            .map(|library| {
                (
                    library.individual_id().as_str(),
                    library.library_id().as_str(),
                )
            })
            .collect::<BTreeSet<_>>();
        if coverage_libraries.len() != manifest.coverage.len()
            || coverage_libraries != planned_libraries
        {
            return Err(event_error(
                &path,
                "cobertura por biblioteca incompleta o repetida",
            ));
        }
    } else if !manifest.coverage.is_empty() {
        return Err(event_error(
            &path,
            "un evento no completo no publica cobertura F4",
        ));
    }
    Ok(())
}

fn f4_id_catalog() -> Value {
    json!({
        "version": crate::identity::ID_CATALOG_VERSION,
        "derivations": {
            "v1": ["execution", "reference", "individual", "event", "batch", "library", "read_group", "template"],
            "v2": ["PreHopInstanceId", "TechnicalPairId"],
        }
    })
}

fn validate_f4_id_catalog(value: &Value, path: &Path) -> Result<(), Diagnostic> {
    let catalog = required_object(value, "id_catalog", path)?;
    validate_closed(catalog, &["version", "derivations"], "id_catalog", path)?;
    require_text(
        catalog,
        "version",
        path,
        crate::identity::ID_CATALOG_VERSION,
    )?;
    let derivations = required_object(
        required_value(catalog, "derivations", path)?,
        "id_catalog.derivations",
        path,
    )?;
    validate_closed(derivations, &["v1", "v2"], "id_catalog.derivations", path)?;
    let expected_v1 = [
        "execution",
        "reference",
        "individual",
        "event",
        "batch",
        "library",
        "read_group",
        "template",
    ];
    let expected_v2 = ["PreHopInstanceId", "TechnicalPairId"];
    for (field, expected) in [
        ("v1", expected_v1.as_slice()),
        ("v2", expected_v2.as_slice()),
    ] {
        let values = required_array(required_value(derivations, field, path)?, field, path)?;
        if values.iter().filter_map(Value::as_str).collect::<Vec<_>>() != expected {
            return Err(event_error(path, "id_catalog incompatible"));
        }
    }
    Ok(())
}

fn validate_profile(value: &Value, path: &Path) -> Result<(), Diagnostic> {
    let profile = required_object(value, "profile", path)?;
    validate_closed(
        profile,
        &[
            "id",
            "resource_path",
            "sha256",
            "preparation",
            "read_length",
            "quality",
        ],
        "profile",
        path,
    )?;
    require_text(
        profile,
        "id",
        path,
        crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_V1,
    )?;
    require_text(
        profile,
        "resource_path",
        path,
        crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_RESOURCE,
    )?;
    require_text(
        profile,
        "sha256",
        path,
        crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256,
    )?;
    require_text(profile, "preparation", path, "pcr_free")?;
    if !matches!(
        profile.get("read_length").and_then(Value::as_u64),
        Some(1..=150)
    ) {
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
            &["individual_id", "library_id", "read_group_id", "i7", "i5"],
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
        for field in ["i7", "i5"] {
            let udi = required_string(item, field, path)?;
            if udi.len() != 8
                || !udi
                    .bytes()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            {
                return Err(event_error(path, "UDI de biblioteca inválido"));
            }
        }
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

fn validate_coverage(
    value: &Value,
    selected: &BTreeMap<String, String>,
    path: &Path,
) -> Result<BTreeSet<(String, String)>, Diagnostic> {
    let values = required_array(value, "coverage", path)?;
    let mut coverage_libraries = BTreeSet::new();
    let mut previous = None::<(String, String)>;
    for value in values {
        let entry = required_object(value, "coverage[]", path)?;
        validate_closed(
            entry,
            &[
                "individual_id",
                "library_id",
                "target_x",
                "basic_x",
                "final_x",
                "denominator_bases",
                "n_bases",
            ],
            "coverage[]",
            path,
        )?;
        let individual = required_string(entry, "individual_id", path)?;
        let library = required_string(entry, "library_id", path)?;
        let ordering = (individual.to_owned(), library.to_owned());
        if selected.get(individual).map(String::as_str) != Some(library)
            || previous.as_ref().is_some_and(|prior| prior >= &ordering)
            || !coverage_libraries.insert(ordering.clone())
        {
            return Err(event_error(
                path,
                "coberturas no son únicas, ordenadas o compatibles con las bibliotecas",
            ));
        }
        previous = Some(ordering);
        for field in ["target_x", "basic_x"] {
            validate_nonnegative_finite(entry, field, path)?;
        }
        let final_x = required_object(
            required_value(entry, "final_x", path)?,
            "coverage[].final_x",
            path,
        )?;
        validate_closed(
            final_x,
            &["own_x", "contaminant_x"],
            "coverage[].final_x",
            path,
        )?;
        for field in ["own_x", "contaminant_x"] {
            validate_nonnegative_finite(final_x, field, path)?;
        }
        let _denominator = entry
            .get("denominator_bases")
            .and_then(Value::as_u64)
            .filter(|value| *value > 0)
            .ok_or_else(|| event_error(path, "denominador de cobertura inválido"))?;
        let _n_bases = entry
            .get("n_bases")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento N de cobertura inválido"))?;
    }
    Ok(coverage_libraries)
}

fn validate_nonnegative_finite(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if !matches!(object.get(field).and_then(Value::as_f64), Some(value) if value.is_finite() && value >= 0.0)
    {
        return Err(event_error(path, format!("{field} de cobertura inválida")));
    }
    Ok(())
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

/// Valida los únicos artefactos técnicos agregados activos en F4.3. El
/// contrato es deliberadamente cerrado para que no se filtren listas de
/// parejas, identificadores o coordenadas al manifiesto durable.
fn validate_technical_artifacts(value: &Value, path: &Path) -> Result<(), Diagnostic> {
    let artifacts = required_object(value, "artifacts", path)?;
    validate_closed(
        artifacts,
        &[
            "chimera",
            "biological_contamination",
            "pcr",
            "optical_duplicates",
            "gc",
            "adapters",
        ],
        "artifacts",
        path,
    )?;
    let chimera = required_object(
        required_value(artifacts, "chimera", path)?,
        "artifacts.chimera",
        path,
    )?;
    validate_biological_contamination_manifest(
        required_value(artifacts, "biological_contamination", path)?,
        path,
    )?;
    validate_closed(
        chimera,
        &["rate", "realized_pairs"],
        "artifacts.chimera",
        path,
    )?;
    let rate = chimera.get("rate").and_then(Value::as_f64);
    let realized_pairs = chimera.get("realized_pairs").and_then(Value::as_u64);
    if !matches!(rate, Some(rate) if rate.is_finite() && (0.0..=1.0).contains(&rate))
        || realized_pairs.is_none()
        || (rate == Some(0.0) && realized_pairs != Some(0))
    {
        return Err(event_error(path, "artefacto de quimera inválido"));
    }
    let pcr = required_object(
        required_value(artifacts, "pcr", path)?,
        "artifacts.pcr",
        path,
    )?;
    let optical = required_object(
        required_value(artifacts, "optical_duplicates", path)?,
        "artifacts.optical_duplicates",
        path,
    )?;
    let gc = required_object(required_value(artifacts, "gc", path)?, "artifacts.gc", path)?;
    let adapters = required_object(
        required_value(artifacts, "adapters", path)?,
        "artifacts.adapters",
        path,
    )?;
    validate_closed(pcr, &["applicable"], "artifacts.pcr", path)?;
    validate_closed(
        optical,
        &["realized_pairs"],
        "artifacts.optical_duplicates",
        path,
    )?;
    validate_closed(gc, &["mode"], "artifacts.gc", path)?;
    validate_closed(adapters, &["enabled"], "artifacts.adapters", path)?;
    if pcr.get("applicable") != Some(&Value::Bool(false))
        || optical.get("realized_pairs").and_then(Value::as_u64) != Some(0)
        || gc.get("mode").and_then(Value::as_str) != Some("neutral")
        || adapters.get("enabled") != Some(&Value::Bool(false))
    {
        return Err(event_error(path, "artefactos técnicos fijos inválidos"));
    }
    Ok(())
}

fn validate_biological_contamination_manifest(
    value: &Value,
    path: &Path,
) -> Result<(), Diagnostic> {
    let records = required_array(value, "artifacts.biological_contamination", path)?;
    let mut recipients = BTreeSet::new();
    for value in records {
        let record = required_object(value, "artifacts.biological_contamination[]", path)?;
        validate_closed(
            record,
            &[
                "recipient_id",
                "recipient_library_id",
                "percent_range",
                "sampled_percent",
                "clean_pair_count",
                "target_pair_count",
                "removed_own_pair_count",
                "inserted_contaminant_pair_count",
                "retained_own_pair_count",
                "final_pair_count",
                "subseeds",
                "donor_contributions",
            ],
            "artifacts.biological_contamination[]",
            path,
        )?;
        let recipient = required_string(record, "recipient_id", path)?;
        let _ = parse_id::<IndividualId>(recipient, "individual", path)?;
        let _ = parse_id::<LibraryId>(
            required_string(record, "recipient_library_id", path)?,
            "library",
            path,
        )?;
        if !recipients.insert(recipient.to_owned()) {
            return Err(event_error(
                path,
                "receptora de contaminación biológica repetida",
            ));
        }
        let range = required_object(
            required_value(record, "percent_range", path)?,
            "percent_range",
            path,
        )?;
        validate_closed(range, &["minimum", "maximum"], "percent_range", path)?;
        let minimum = range.get("minimum").and_then(Value::as_f64);
        let maximum = range.get("maximum").and_then(Value::as_f64);
        let sampled = record.get("sampled_percent").and_then(Value::as_f64);
        if !matches!((minimum, maximum, sampled), (Some(minimum), Some(maximum), Some(sampled)) if minimum.is_finite() && maximum.is_finite() && sampled.is_finite() && (0.0..=100.0).contains(&minimum) && (0.0..=100.0).contains(&maximum) && minimum <= sampled && sampled <= maximum)
        {
            return Err(event_error(
                path,
                "porcentaje de contaminación biológica inválido",
            ));
        }
        let sampled = sampled
            .ok_or_else(|| event_error(path, "porcentaje de contaminación biológica ausente"))?;
        let clean = record
            .get("clean_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento limpio de contaminación inválido"))?;
        let target = record
            .get("target_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "objetivo de contaminación inválido"))?;
        let removed = record
            .get("removed_own_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento retirado de contaminación inválido"))?;
        let inserted = record
            .get("inserted_contaminant_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento insertado de contaminación inválido"))?;
        let retained = record
            .get("retained_own_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento retenido de contaminación inválido"))?;
        let final_pairs = record
            .get("final_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento final de contaminación inválido"))?;
        if target > clean {
            return Err(event_error(
                path,
                "objetivo de contaminación supera parejas limpias",
            ));
        }
        let expected_target = contamination_target_pair_count(clean, sampled)
            .map_err(|_| event_error(path, "objetivo de contaminación no es recomputable"))?;
        if target != expected_target {
            return Err(event_error(
                path,
                "objetivo de contaminación no aplica el redondeo normativo",
            ));
        }
        if removed != target
            || inserted != target
            || retained.checked_add(removed) != Some(clean)
            || retained.checked_add(inserted) != Some(final_pairs)
            || final_pairs != clean
        {
            return Err(event_error(
                path,
                "los contadores realizados de contaminación no conservan las parejas",
            ));
        }
        let subseeds =
            required_object(required_value(record, "subseeds", path)?, "subseeds", path)?;
        validate_closed(
            subseeds,
            &["percent", "donors", "removal"],
            "subseeds",
            path,
        )?;
        for value in ["percent", "donors", "removal"] {
            validate_subseed_hex(required_string(subseeds, value, path)?, value, path)?;
        }
        let contributions = required_array(
            required_value(record, "donor_contributions", path)?,
            "donor_contributions",
            path,
        )?;
        let mut donors = BTreeSet::new();
        let mut total = 0_u64;
        for value in contributions {
            let contribution = required_object(value, "donor_contributions[]", path)?;
            validate_closed(
                contribution,
                &["individual_id", "pairs"],
                "donor_contributions[]",
                path,
            )?;
            let donor = required_string(contribution, "individual_id", path)?;
            let _ = parse_id::<IndividualId>(donor, "individual", path)?;
            let pairs = contribution
                .get("pairs")
                .and_then(Value::as_u64)
                .ok_or_else(|| event_error(path, "aporte de donante inválido"))?;
            if donor == recipient || pairs == 0 || !donors.insert(donor.to_owned()) {
                return Err(event_error(path, "aporte de donante inválido"));
            }
            total = total
                .checked_add(pairs)
                .ok_or_else(|| event_error(path, "desbordamiento de aportes"))?;
        }
        if total != target {
            return Err(event_error(
                path,
                "los aportes no suman el objetivo de contaminación",
            ));
        }
    }
    Ok(())
}

fn validate_subseed_hex(value: &str, field: &str, path: &Path) -> Result<(), Diagnostic> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        return Err(event_error(
            path,
            format!("{field} no es una subsemilla hexadecimal de 64 caracteres"),
        ));
    }
    Ok(())
}

fn validate_outputs(
    value: &Value,
    event_id: &EventId,
    selected: &BTreeMap<String, String>,
    reference: &ArtifactRecord,
    diagnostics: &BTreeSet<String>,
    path: &Path,
) -> Result<Vec<UnitState>, Diagnostic> {
    let values = required_array(value, "outputs", path)?;
    if values.is_empty() {
        return Err(event_error(path, "el evento no tiene resultados FASTQ"));
    }
    let mut ids = BTreeSet::new();
    let mut paths = BTreeSet::new();
    let mut artifact_mates = BTreeSet::new();
    let mut states = Vec::with_capacity(values.len());
    let mut previous = None::<(String, String, u8, String)>;
    for value in values {
        let artifact = required_object(value, "outputs[]", path)?;
        validate_closed(
            artifact,
            &[
                "id",
                "kind",
                "state",
                "individual_id",
                "library_id",
                "mate",
                "artifact",
                "diagnostic_code",
            ],
            "outputs[]",
            path,
        )?;
        let id = required_string(artifact, "id", path)?;
        let kind = required_string(artifact, "kind", path)?;
        if kind != FASTQ_ARTIFACT_KIND {
            return Err(event_error(path, "F3 solo admite artefactos FASTQ"));
        }
        let state = parse_state(required_string(artifact, "state", path)?, path)?;
        let individual = required_string(artifact, "individual_id", path)?;
        let library = required_string(artifact, "library_id", path)?;
        let mate = parse_mate(artifact, path)?;
        let ordering = (
            individual.to_owned(),
            library.to_owned(),
            mate.number(),
            id.to_owned(),
        );
        if selected.get(individual).map(String::as_str) != Some(library)
            || !ids.insert(id.to_owned())
            || !artifact_mates.insert((individual.to_owned(), mate))
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
                    || record.path != paired_path(event_id.as_str(), individual, mate)
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
    if artifact_mates.len() != selected.len() * 2
        || !selected.keys().all(|individual| {
            TechnicalMate::all()
                .into_iter()
                .all(|mate| artifact_mates.contains(&(individual.clone(), mate)))
        })
    {
        return Err(event_error(
            path,
            "falta el resultado FASTQ de una biblioteca seleccionada",
        ));
    }
    Ok(states)
}

fn paired_path(event_id: &str, individual_id: &str, mate: TechnicalMate) -> String {
    let (r1, r2) = paired_fastq_relative_paths(event_id, individual_id);
    match mate {
        TechnicalMate::R1 => r1,
        TechnicalMate::R2 => r2,
    }
    .to_string_lossy()
    .into_owned()
}

fn parse_mate(
    artifact: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<TechnicalMate, Diagnostic> {
    match artifact.get("mate").and_then(Value::as_u64) {
        Some(1) => Ok(TechnicalMate::R1),
        Some(2) => Ok(TechnicalMate::R2),
        _ => Err(event_error(path, "mate FASTQ inválido")),
    }
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
    .with("schema", "manifest-event-v2.schema.json")
    .with("detail", detail.into())
}
