//! Preparación, ensamblado y cierre de una ejecución F1.4 validada.

use std::path::{Path, PathBuf};

use crate::application::{EffectiveConfig, EffectiveConfigurationInput};
use crate::diagnostics::Diagnostic;
use crate::identity::{AutomaticIdRequest, IdClass, IdentityRegistry};
use crate::integrity::Sha256Writer;
use crate::publication::{
    ArtifactStage, AtomicPublisher, ExecutionWorkspace, LinuxAtomicPublisher, PublishableUnit,
};
use crate::results::{
    ArtifactRecord, ConfigurationInputRecord, EventManifest, EventManifestLink, GeneralManifest,
    UnitRecord, UnitState, effective_config_document, event_manifest_document,
    general_manifest_document, json_document_bytes,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ExecutionFault {
    None,
    AfterConfigurationCopy,
    AfterEffectiveConfig,
    InterruptedAfterWorkspace,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ExecutionOutcome {
    Published {
        destination: PathBuf,
        state: UnitState,
    },
    Interrupted {
        provisional: PathBuf,
    },
}

pub fn execute_validated_config(
    effective: &EffectiveConfig,
    destination: &Path,
) -> Result<ExecutionOutcome, Diagnostic> {
    execute_validated_config_with(
        effective,
        destination,
        ExecutionFault::None,
        &LinuxAtomicPublisher,
    )
}

pub fn execute_validated_config_with(
    effective: &EffectiveConfig,
    destination: &Path,
    fault: ExecutionFault,
    publisher: &impl AtomicPublisher,
) -> Result<ExecutionOutcome, Diagnostic> {
    let effective_document = effective_config_document(effective);
    let effective_bytes = json_document_bytes(&effective_document)?;
    let configuration_digest = configuration_digest(effective.configuration_inputs())?;
    let execution_request = AutomaticIdRequest::new(
        IdClass::Execution,
        "execution",
        &[("config_sha256", configuration_digest.as_str())],
    )?;
    let mut identities = IdentityRegistry::default();
    let execution_id =
        identities.derive_execution(effective.reproducibility().root_seed(), &execution_request)?;

    let mut workspace = ExecutionWorkspace::create(destination)?;
    if fault == ExecutionFault::InterruptedAfterWorkspace {
        return Ok(ExecutionOutcome::Interrupted {
            provisional: workspace.root().to_path_buf(),
        });
    }

    let mut configuration_inputs = copy_configuration_inputs(effective, &mut workspace)?;
    configuration_inputs.sort_by(|first, second| {
        (&first.kind, &first.copy.id).cmp(&(&second.kind, &second.copy.id))
    });
    if fault == ExecutionFault::AfterConfigurationCopy {
        return close_failed_execution(
            effective,
            workspace,
            publisher,
            execution_id.as_str(),
            &configuration_inputs,
            None,
            &[],
            injected_failure("after_configuration_copy"),
        );
    }

    let effective_record = workspace.write_artifact(
        "effective-config",
        "effective_config",
        Path::new("effective-config.json"),
        &effective_bytes,
    )?;
    if fault == ExecutionFault::AfterEffectiveConfig {
        return close_failed_execution(
            effective,
            workspace,
            publisher,
            execution_id.as_str(),
            &configuration_inputs,
            Some(&effective_record),
            &[],
            injected_failure("after_effective_config"),
        );
    }

    let event_links = write_event_manifests(effective, &mut workspace)?;
    let mut execution_unit = PublishableUnit::new(
        execution_id.as_str(),
        "execution",
        std::iter::once("effective-config".to_owned())
            .chain(
                configuration_inputs
                    .iter()
                    .map(|input| input.copy.id.clone()),
            )
            .chain(event_links.iter().map(|event| event.manifest.id.clone())),
    );
    execution_unit.record_artifact("effective-config", ArtifactStage::Cataloged)?;
    for input in &configuration_inputs {
        execution_unit.record_artifact(&input.copy.id, ArtifactStage::Cataloged)?;
    }
    for event in &event_links {
        execution_unit.record_artifact(&event.manifest.id, ArtifactStage::Cataloged)?;
    }
    execution_unit.complete()?;
    let units = std::iter::once(execution_unit.record())
        .chain(event_links.iter().map(|event| UnitRecord {
            id: event.event_id.clone(),
            kind: "event".to_owned(),
            state: UnitState::Complete,
            exportable: true,
            diagnostic_code: None,
        }))
        .collect::<Vec<_>>();
    let artifacts = workspace.artifacts().to_vec();
    let general = GeneralManifest {
        execution_id: execution_id.as_str(),
        state: UnitState::Complete,
        effective,
        configuration_inputs: &configuration_inputs,
        effective_config: Some(&effective_record),
        events: &event_links,
        artifacts: &artifacts,
        units: &units,
        diagnostics: &[],
    };
    let bytes = json_document_bytes(&general_manifest_document(&general))?;
    workspace.write_artifact(
        "manifest-general",
        "manifest_general",
        Path::new("manifest.json"),
        &bytes,
    )?;
    let destination = workspace.publish(publisher)?;
    Ok(ExecutionOutcome::Published {
        destination,
        state: UnitState::Complete,
    })
}

fn copy_configuration_inputs(
    effective: &EffectiveConfig,
    workspace: &mut ExecutionWorkspace,
) -> Result<Vec<ConfigurationInputRecord>, Diagnostic> {
    let mut inputs = Vec::new();
    for input in effective.configuration_inputs() {
        let (artifact_id, artifact_kind, relative) = if input.kind() == "general" {
            (
                "config-general".to_owned(),
                "config_general_toml",
                PathBuf::from("inputs/config/general.toml"),
            )
        } else {
            (
                format!("config-scenario-{}", input.id()),
                "config_scenario_toml",
                PathBuf::from(format!("inputs/config/scenarios/{}.toml", input.id())),
            )
        };
        let copy =
            workspace.write_artifact(&artifact_id, artifact_kind, &relative, input.bytes())?;
        inputs.push(ConfigurationInputRecord {
            kind: input.kind().to_owned(),
            original_path: input.original_path().to_owned(),
            copy,
        });
    }
    Ok(inputs)
}

fn configuration_digest(inputs: &[EffectiveConfigurationInput]) -> Result<String, Diagnostic> {
    let mut hasher = Sha256Writer::default();
    hasher.update(b"genome_population_simulator/v1/configuration-digest\0")?;
    hasher.update(
        &u64::try_from(inputs.len())
            .map_err(|_| digest_material_error("input_count"))?
            .to_be_bytes(),
    )?;
    for input in inputs {
        update_length_prefixed(&mut hasher, input.kind().as_bytes(), "input_kind")?;
        update_length_prefixed(&mut hasher, input.id().as_bytes(), "input_id")?;
        update_length_prefixed(&mut hasher, input.bytes(), "input_bytes")?;
    }
    Ok(hasher.finish("configuration-inputs").digest().to_owned())
}

fn update_length_prefixed(
    hasher: &mut Sha256Writer,
    value: &[u8],
    field: &str,
) -> Result<(), Diagnostic> {
    let length = u64::try_from(value.len()).map_err(|_| digest_material_error(field))?;
    hasher.update(&length.to_be_bytes())?;
    hasher.update(value)
}

fn digest_material_error(field: &str) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "execution")
    .with("operation", "configuration_digest")
    .with("detail", format!("{field} exceeds u64"))
}

fn write_event_manifests(
    effective: &EffectiveConfig,
    workspace: &mut ExecutionWorkspace,
) -> Result<Vec<EventManifestLink>, Diagnostic> {
    let mut links = Vec::with_capacity(effective.reproducibility().event_subseeds().len());
    for event in effective.reproducibility().event_subseeds() {
        let subseed = event.subseed().to_hex();
        let event_unit = UnitRecord {
            id: event.event_id().to_owned(),
            kind: "event".to_owned(),
            state: UnitState::Complete,
            exportable: true,
            diagnostic_code: None,
        };
        let manifest = EventManifest {
            event_id: event.event_id(),
            state: UnitState::Complete,
            subseed: &subseed,
            artifacts: &[],
            units: &[event_unit],
            diagnostics: &[],
        };
        let bytes = json_document_bytes(&event_manifest_document(&manifest))?;
        let relative = PathBuf::from(format!("events/{}/manifest.json", event.event_id()));
        let record = workspace.write_artifact(
            &format!("manifest-event-{}", event.event_id()),
            "manifest_event",
            &relative,
            &bytes,
        )?;
        links.push(EventManifestLink {
            event_id: event.event_id().to_owned(),
            subseed,
            state: UnitState::Complete,
            manifest: record,
        });
    }
    links.sort_by(|first, second| first.event_id.cmp(&second.event_id));
    Ok(links)
}

#[allow(clippy::too_many_arguments)]
fn close_failed_execution(
    effective: &EffectiveConfig,
    mut workspace: ExecutionWorkspace,
    publisher: &impl AtomicPublisher,
    execution_id: &str,
    configuration_inputs: &[ConfigurationInputRecord],
    effective_config: Option<&ArtifactRecord>,
    events: &[EventManifestLink],
    diagnostic: Diagnostic,
) -> Result<ExecutionOutcome, Diagnostic> {
    let units = vec![UnitRecord {
        id: execution_id.to_owned(),
        kind: "execution".to_owned(),
        state: UnitState::Failed,
        exportable: false,
        diagnostic_code: Some(diagnostic.code().to_owned()),
    }];
    let artifacts = workspace.artifacts().to_vec();
    let diagnostics = [diagnostic];
    let manifest = GeneralManifest {
        execution_id,
        state: UnitState::Failed,
        effective,
        configuration_inputs,
        effective_config,
        events,
        artifacts: &artifacts,
        units: &units,
        diagnostics: &diagnostics,
    };
    let bytes = json_document_bytes(&general_manifest_document(&manifest))?;
    workspace.write_artifact(
        "manifest-general",
        "manifest_general",
        Path::new("manifest.json"),
        &bytes,
    )?;
    let destination = workspace.publish(publisher)?;
    Ok(ExecutionOutcome::Published {
        destination,
        state: UnitState::Failed,
    })
}

fn injected_failure(stage: &str) -> Diagnostic {
    Diagnostic::domain(
        "GPS-DOMAIN-003",
        "Una invariante del estado biológico o técnico válido no se cumple.",
    )
    .with("entity", "execution")
    .with("invariant", "controlled_failure_injection")
    .with("detail", stage)
}
