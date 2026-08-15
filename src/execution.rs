//! Preparación, ensamblado y cierre de una ejecución F1.4 validada.

use std::fs;
use std::path::{Path, PathBuf};

use crate::application::{EffectiveConfig, EffectiveConfigurationInput, EffectiveReferenceSource};
use crate::diagnostics::Diagnostic;
use crate::identity::{AutomaticIdRequest, IdClass, IdentityRegistry, ReferenceId};
use crate::integrity::Sha256Writer;
use crate::publication::{
    ArtifactStage, AtomicPublisher, ExecutionWorkspace, LinuxAtomicPublisher, PublishableUnit,
};
use crate::reference::Reference;
use crate::reference_fasta::{
    SequenceExpectation, ValidatedSuppliedReference, validate_supplied_fasta,
    validate_supplied_fasta_with_origin, write_generated_fasta,
};
use crate::reference_generation::{BaseFrequency, GeneratedSequenceSpec, generate_reference};
use crate::reference_reduction::{
    EligibleSequence, MaterializedPlan, ReductionError, materialize_selection,
    preflight_reduced_reference, select_reduced_reference,
};
use crate::reference_reduction_io::{
    ReductionMapParameters, build_reduced_reference, reduction_map_document,
    validate_reduced_reference, validate_reduction_map_document,
};
use crate::results::{
    ArtifactRecord, ConfigurationInputRecord, EventManifest, EventManifestLink, GeneralManifest,
    ReducedReferenceOutput, ReferenceOutput, UnitRecord, UnitState,
    effective_config_document_with_reductions, event_manifest_document, general_manifest_document,
    json_document_bytes, resolution_evidence_document,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) enum ReductionFaultStage {
    Fasta,
    Map,
    Metadata,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) enum ExecutionFault {
    None,
    AfterConfigurationCopy,
    AfterEffectiveConfig,
    AtReduction {
        reduction_index: usize,
        stage: ReductionFaultStage,
    },
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

pub(crate) fn execute_validated_config_with(
    effective: &EffectiveConfig,
    destination: &Path,
    fault: ExecutionFault,
    publisher: &impl AtomicPublisher,
) -> Result<ExecutionOutcome, Diagnostic> {
    let mut prepared_reference = prepare_reference(effective)?;
    let resolution_evidence_bytes = json_document_bytes(&resolution_evidence_document(effective))?;
    let configuration_digest =
        configuration_digest(effective.configuration_inputs(), &resolution_evidence_bytes)?;
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
    let resolution_evidence = workspace.write_artifact(
        "resolution-evidence",
        "resolution_evidence",
        Path::new("inputs/resolution/evidence.json"),
        &resolution_evidence_bytes,
    )?;
    if fault == ExecutionFault::AfterConfigurationCopy {
        return close_failed_execution(
            effective,
            workspace,
            publisher,
            execution_id.as_str(),
            &configuration_inputs,
            &resolution_evidence,
            None,
            None,
            &[],
            &[],
            None,
            injected_failure("after_configuration_copy"),
        );
    }

    let reference_artifact = write_reference(effective, &prepared_reference, &mut workspace)?;
    prepared_reference = load_authoritative_reference(
        effective,
        &workspace,
        &reference_artifact,
        prepared_reference,
    )?;
    let reference_output = ReferenceOutput {
        reference: prepared_reference.reference(),
        fasta: &reference_artifact,
    };
    let prepared_reductions = prepare_reductions(effective, reference_output.reference)?;
    let mut published_reductions = Vec::with_capacity(prepared_reductions.len());
    for (effective_index, mut prepared) in prepared_reductions.into_iter().enumerate() {
        let effective_reduction = &effective.reduced_references()[effective_index];
        let fasta_relative = PathBuf::from(format!(
            "resources/references/{}.fa",
            effective_reduction.id()
        ));
        let fasta = workspace.write_streamed_artifact(
            &format!("reference-fasta-{}", effective_reduction.id()),
            "reference_fasta",
            &fasta_relative,
            |file| write_generated_fasta(&prepared.reference, file),
        )?;
        prepared
            .reference
            .bind_authoritative_resource(PathBuf::from(&fasta.path))
            .map_err(|error| reduction_internal("bind_reduced_authoritative_resource", error))?;
        if reduction_fault_matches(fault, effective_index, ReductionFaultStage::Fasta) {
            let completed_reductions =
                reduction_outputs(effective, &published_reductions, &reference_artifact, None);
            return close_failed_execution(
                effective,
                workspace,
                publisher,
                execution_id.as_str(),
                &configuration_inputs,
                &resolution_evidence,
                None,
                Some(reference_output),
                &[],
                &completed_reductions,
                Some(effective_index),
                injected_failure("after_reduced_fasta"),
            );
        }
        let map_document = reduction_map_document(&prepared.plan);
        let map_bytes = json_document_bytes(&map_document)?;
        let map_relative = PathBuf::from(format!(
            "resources/references/{}.map.json",
            effective_reduction.id()
        ));
        let map = workspace.write_artifact(
            &format!("reference-map-{}", effective_reduction.id()),
            "reference_reduction_map",
            &map_relative,
            &map_bytes,
        )?;
        let serialized_map = fs::read(workspace.root().join(&map.path)).map_err(|error| {
            Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
                .with("path", map.path.clone())
                .with("operation", "read_serialized_reduction_map")
                .with("detail", error.to_string())
        })?;
        let reparsed_map: serde_json::Value = serde_json::from_slice(&serialized_map).map_err(
            |error| {
                Diagnostic::format(
                    "GPS-FMT-007",
                    "El mapa de reducción no permite reconstruir exactamente la referencia derivada.",
                )
                .with("path", map.path.clone())
                .with("reference", effective_reduction.id())
                .with("invariant", "serialized_json")
                .with("detail", error.to_string())
            },
        )?;
        validate_reduction_map_document(
            reference_output.reference,
            &prepared.reference,
            &reparsed_map,
            ReductionMapParameters {
                region_size_bases: effective_reduction.region_size_bases(),
                target_bases: effective_reduction.target_bases(),
            },
        )?;
        if reduction_fault_matches(fault, effective_index, ReductionFaultStage::Map) {
            let completed_reductions =
                reduction_outputs(effective, &published_reductions, &reference_artifact, None);
            return close_failed_execution(
                effective,
                workspace,
                publisher,
                execution_id.as_str(),
                &configuration_inputs,
                &resolution_evidence,
                None,
                Some(reference_output),
                &[],
                &completed_reductions,
                Some(effective_index),
                injected_failure("after_reduced_map"),
            );
        }
        published_reductions.push(PublishedReduction {
            effective_index,
            reference: prepared.reference,
            plan: prepared.plan,
            fasta,
            map,
        });
    }
    let metadata_failure_index = match fault {
        ExecutionFault::AtReduction {
            reduction_index,
            stage: ReductionFaultStage::Metadata,
        } if reduction_index < published_reductions.len() => Some(reduction_index),
        _ => None,
    };
    let reduction_outputs = reduction_outputs(
        effective,
        &published_reductions,
        &reference_artifact,
        metadata_failure_index,
    );
    let effective_document =
        effective_config_document_with_reductions(effective, reference_output, &reduction_outputs)?;
    let effective_bytes = json_document_bytes(&effective_document)?;

    let effective_record = workspace.write_artifact(
        "effective-config",
        "effective_config",
        Path::new("effective-config.json"),
        &effective_bytes,
    )?;
    if fault == ExecutionFault::AfterEffectiveConfig || metadata_failure_index.is_some() {
        return close_failed_execution(
            effective,
            workspace,
            publisher,
            execution_id.as_str(),
            &configuration_inputs,
            &resolution_evidence,
            Some(&effective_record),
            Some(reference_output),
            &[],
            &reduction_outputs,
            metadata_failure_index,
            injected_failure(if metadata_failure_index.is_some() {
                "after_reduced_metadata"
            } else {
                "after_effective_config"
            }),
        );
    }

    let event_links = write_event_manifests(effective, reference_output, &mut workspace)?;
    let mut execution_unit = PublishableUnit::new(
        execution_id.as_str(),
        "execution",
        std::iter::once("effective-config".to_owned())
            .chain(std::iter::once("resolution-evidence".to_owned()))
            .chain(
                configuration_inputs
                    .iter()
                    .map(|input| input.copy.id.clone()),
            )
            .chain(event_links.iter().map(|event| event.manifest.id.clone()))
            .chain(std::iter::once(reference_artifact.id.clone()))
            .chain(
                published_reductions
                    .iter()
                    .flat_map(|reduction| [reduction.fasta.id.clone(), reduction.map.id.clone()]),
            ),
    );
    execution_unit.record_artifact("effective-config", ArtifactStage::Cataloged)?;
    execution_unit.record_artifact("resolution-evidence", ArtifactStage::Cataloged)?;
    for input in &configuration_inputs {
        execution_unit.record_artifact(&input.copy.id, ArtifactStage::Cataloged)?;
    }
    for event in &event_links {
        execution_unit.record_artifact(&event.manifest.id, ArtifactStage::Cataloged)?;
    }
    execution_unit.record_artifact(&reference_artifact.id, ArtifactStage::Cataloged)?;
    for reduction in &published_reductions {
        execution_unit.record_artifact(&reduction.fasta.id, ArtifactStage::Cataloged)?;
        execution_unit.record_artifact(&reduction.map.id, ArtifactStage::Cataloged)?;
    }
    execution_unit.complete()?;
    let units = std::iter::once(execution_unit.record())
        .chain(
            effective
                .reduced_references()
                .iter()
                .map(|reduction| UnitRecord {
                    id: reduction.id().to_owned(),
                    kind: "reference_reduction".to_owned(),
                    state: UnitState::Complete,
                    exportable: true,
                    diagnostic_code: None,
                }),
        )
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
        resolution_evidence: &resolution_evidence,
        effective_config: Some(&effective_record),
        reference: Some(reference_output),
        reduced_references: &reduction_outputs,
        events: &event_links,
        artifacts: &artifacts,
        units: &units,
        diagnostics: &[],
    };
    write_validated_general_manifest(&mut workspace, &general_manifest_document(&general))?;
    let destination = workspace.publish_with(publisher)?;
    Ok(ExecutionOutcome::Published {
        destination,
        state: UnitState::Complete,
    })
}

enum PreparedReference {
    Supplied(ValidatedSuppliedReference),
    Generated(Reference),
}

struct PreparedReduction {
    reference: Reference,
    plan: MaterializedPlan,
}

struct PublishedReduction {
    effective_index: usize,
    reference: Reference,
    plan: MaterializedPlan,
    fasta: ArtifactRecord,
    map: ArtifactRecord,
}

fn reduction_fault_matches(
    fault: ExecutionFault,
    reduction_index: usize,
    stage: ReductionFaultStage,
) -> bool {
    matches!(
        fault,
        ExecutionFault::AtReduction {
            reduction_index: selected,
            stage: selected_stage,
        } if selected == reduction_index && selected_stage == stage
    )
}

fn reduction_outputs<'a>(
    effective: &'a EffectiveConfig,
    published: &'a [PublishedReduction],
    source_fasta: &'a ArtifactRecord,
    failure_start: Option<usize>,
) -> Vec<ReducedReferenceOutput<'a>> {
    published
        .iter()
        .map(|reduction| ReducedReferenceOutput {
            effective: &effective.reduced_references()[reduction.effective_index],
            reference: &reduction.reference,
            fasta: &reduction.fasta,
            map: &reduction.map,
            plan: &reduction.plan,
            source_fasta,
            state: if failure_start.is_some_and(|start| reduction.effective_index >= start) {
                UnitState::Failed
            } else {
                UnitState::Complete
            },
        })
        .collect()
}

impl PreparedReference {
    fn reference(&self) -> &Reference {
        match self {
            Self::Supplied(supplied) => supplied.reference(),
            Self::Generated(reference) => reference,
        }
    }

    fn reference_mut(&mut self) -> &mut Reference {
        match self {
            Self::Supplied(supplied) => supplied.reference_mut(),
            Self::Generated(reference) => reference,
        }
    }
}

fn prepare_reference(effective: &EffectiveConfig) -> Result<PreparedReference, Diagnostic> {
    match effective.reference_source() {
        EffectiveReferenceSource::Supplied { original_path } => {
            let expected = sequence_expectations(effective);
            validate_supplied_fasta(original_path, effective.reference_id(), &expected)
                .map(PreparedReference::Supplied)
        }
        EffectiveReferenceSource::Simulated => {
            let specs = effective
                .sequences()
                .iter()
                .map(|sequence| GeneratedSequenceSpec {
                    id: sequence.id(),
                    class: sequence.class(),
                    length: sequence.length(),
                })
                .collect::<Vec<_>>();
            let normalization = effective
                .normalized_collections()
                .iter()
                .find(|collection| collection.parameter == "reference.composition")
                .ok_or_else(|| {
                    Diagnostic::internal(
                        "GPS-INTERNAL-001",
                        "Una invariante interna considerada imposible se ha roto.",
                    )
                    .with("component", "execution")
                    .with("operation", "prepare_reference")
                    .with("detail", "simulated reference has no effective composition")
                })?;
            let composition = normalization
                .values
                .iter()
                .map(|value| {
                    let base = match value.id.as_str() {
                        "A" => b'A',
                        "C" => b'C',
                        "G" => b'G',
                        "T" => b'T',
                        other => {
                            return Err(Diagnostic::internal(
                                "GPS-INTERNAL-001",
                                "Una invariante interna considerada imposible se ha roto.",
                            )
                            .with("component", "execution")
                            .with("operation", "prepare_reference")
                            .with("detail", format!("unexpected admitted base ID {other}")));
                        }
                    };
                    Ok(BaseFrequency {
                        base,
                        probability: value.after,
                    })
                })
                .collect::<Result<Vec<_>, Diagnostic>>()?;
            generate_reference(
                effective.reference_id(),
                effective.reproducibility().root_seed(),
                &specs,
                &composition,
            )
            .map(PreparedReference::Generated)
        }
    }
}

fn prepare_reductions(
    effective: &EffectiveConfig,
    source: &Reference,
) -> Result<Vec<PreparedReduction>, Diagnostic> {
    if effective.reduced_references().is_empty() {
        return Ok(Vec::new());
    }
    let source_id = source.id().clone();
    let eligible = source
        .sequences()
        .iter()
        .map(|sequence| {
            EligibleSequence::new(sequence.id().clone(), sequence.class(), sequence.length())
                .map_err(|error| reduction_error("source", effective, error))
        })
        .collect::<Result<Vec<_>, _>>()?;
    for reduction in effective.reduced_references() {
        preflight_reduced_reference(
            &eligible,
            reduction.target_bases(),
            reduction.region_size_bases(),
        )
        .map_err(|error| reduction_error(reduction.id(), effective, error))?;
    }
    let mut prepared = Vec::with_capacity(effective.reduced_references().len());
    for reduction in effective.reduced_references() {
        let derived_id = ReferenceId::parse(reduction.id(), "execution")?;
        let mut rng = reduction.selection_subseed().local_rng();
        let selection = select_reduced_reference(
            &source_id,
            &eligible,
            reduction.target_bases(),
            reduction.region_size_bases(),
            &mut rng,
        )
        .map_err(|error| reduction_error(reduction.id(), effective, error))?;
        let plan = materialize_selection(&selection, &derived_id, &eligible)
            .map_err(|error| reduction_error(reduction.id(), effective, error))?;
        let reference = build_reduced_reference(source, &plan)?;
        validate_reduced_reference(source, &reference, &plan)?;
        prepared.push(PreparedReduction { reference, plan });
    }
    Ok(prepared)
}

fn reduction_error(
    reference: &str,
    effective: &EffectiveConfig,
    error: ReductionError,
) -> Diagnostic {
    match error {
        ReductionError::AllocationFailed => {
            return Diagnostic::resource(
                "GPS-RES-002",
                "Los recursos disponibles no permiten completar el plan sin degradar resultados obligatorios.",
            )
            .with("resource", "memory")
            .with("required", "complete immutable reduction plan")
            .with("available", "allocation rejected");
        }
        ReductionError::ArithmeticOverflow
        | ReductionError::TooManySegments
        | ReductionError::Coordinate(_)
        | ReductionError::SelectionDoesNotMatchCatalog => {
            return reduction_internal("plan_reduced_reference", error);
        }
        ReductionError::EmptyEligibleSpace
        | ReductionError::ZeroLengthSequence(_)
        | ReductionError::DuplicateSequenceId(_)
        | ReductionError::ZeroRegionSize
        | ReductionError::ZeroTarget
        | ReductionError::TargetNotSmallerThanEligible { .. }
        | ReductionError::InfeasibleLayout => {}
    }
    Diagnostic::configuration(
        "GPS-CONF-022",
        "No existe una disposición completa de intervalos para la reducción solicitada.",
    )
    .with("reference", reference)
    .with(
        "target_bases",
        effective
            .reduced_references()
            .iter()
            .find(|candidate| candidate.id() == reference)
            .map_or_else(
                || "unknown".to_owned(),
                |candidate| candidate.target_bases().to_string(),
            ),
    )
    .with(
        "region_size_bases",
        effective
            .reduced_references()
            .iter()
            .find(|candidate| candidate.id() == reference)
            .map_or_else(
                || "unknown".to_owned(),
                |candidate| candidate.region_size_bases().to_string(),
            ),
    )
    .with("eligible_bases", effective.total_genome_size().to_string())
    .with("detail", format!("{error:?}"))
}

fn reduction_internal(operation: &str, error: impl std::fmt::Debug) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "execution")
    .with("operation", operation)
    .with("detail", format!("{error:?}"))
}

fn sequence_expectations(effective: &EffectiveConfig) -> Vec<SequenceExpectation<'_>> {
    effective
        .sequences()
        .iter()
        .map(|sequence| SequenceExpectation {
            id: sequence.id(),
            class: sequence.class(),
            length: sequence.length(),
        })
        .collect()
}

fn load_authoritative_reference(
    effective: &EffectiveConfig,
    workspace: &ExecutionWorkspace,
    artifact: &ArtifactRecord,
    prepared: PreparedReference,
) -> Result<PreparedReference, Diagnostic> {
    let mut authoritative = match (effective.reference_source(), prepared) {
        (EffectiveReferenceSource::Supplied { original_path }, PreparedReference::Supplied(_)) => {
            let read_path = workspace.root().join(&artifact.path);
            let expected = sequence_expectations(effective);
            let supplied = validate_supplied_fasta_with_origin(
                &read_path,
                original_path,
                effective.reference_id(),
                &expected,
            )?;
            if supplied.source_identity().size_bytes() != artifact.size_bytes
                || supplied.source_identity().digest() != artifact.digest
            {
                return Err(Diagnostic::internal(
                    "GPS-INTERNAL-001",
                    "Una invariante interna considerada imposible se ha roto.",
                )
                .with("component", "execution")
                .with("operation", "load_authoritative_reference")
                .with(
                    "detail",
                    "the cataloged reference changed before authoritative loading",
                ));
            }
            PreparedReference::Supplied(supplied)
        }
        (EffectiveReferenceSource::Simulated, PreparedReference::Generated(reference)) => {
            PreparedReference::Generated(reference)
        }
        _ => {
            return Err(Diagnostic::internal(
                "GPS-INTERNAL-001",
                "Una invariante interna considerada imposible se ha roto.",
            )
            .with("component", "execution")
            .with("operation", "load_authoritative_reference")
            .with(
                "detail",
                "prepared reference origin changed before authoritative loading",
            ));
        }
    };
    authoritative
        .reference_mut()
        .bind_authoritative_resource(PathBuf::from(&artifact.path))
        .map_err(|error| {
            Diagnostic::internal(
                "GPS-INTERNAL-001",
                "Una invariante interna considerada imposible se ha roto.",
            )
            .with("component", "execution")
            .with("operation", "bind_authoritative_reference")
            .with("detail", format!("{error:?}"))
        })?;
    Ok(authoritative)
}

fn write_reference(
    effective: &EffectiveConfig,
    prepared: &PreparedReference,
    workspace: &mut ExecutionWorkspace,
) -> Result<ArtifactRecord, Diagnostic> {
    let relative = PathBuf::from(format!(
        "resources/references/{}.fa",
        effective.reference_id()
    ));
    match (effective.reference_source(), prepared) {
        (
            EffectiveReferenceSource::Supplied { original_path },
            PreparedReference::Supplied(supplied),
        ) => workspace.copy_verified_artifact(
            "reference-fasta",
            "reference_fasta",
            original_path,
            &relative,
            supplied.source_identity(),
        ),
        (EffectiveReferenceSource::Simulated, PreparedReference::Generated(reference)) => workspace
            .write_streamed_artifact("reference-fasta", "reference_fasta", &relative, |file| {
                write_generated_fasta(reference, file)
            }),
        _ => Err(Diagnostic::internal(
            "GPS-INTERNAL-001",
            "Una invariante interna considerada imposible se ha roto.",
        )
        .with("component", "execution")
        .with("operation", "write_reference")
        .with(
            "detail",
            "prepared reference origin changed before publication",
        )),
    }
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

fn configuration_digest(
    inputs: &[EffectiveConfigurationInput],
    resolution_evidence: &[u8],
) -> Result<String, Diagnostic> {
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
    update_length_prefixed(
        &mut hasher,
        resolution_evidence,
        "resolution_evidence_bytes",
    )?;
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
    reference: ReferenceOutput<'_>,
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
            reference,
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
    resolution_evidence: &ArtifactRecord,
    effective_config: Option<&ArtifactRecord>,
    reference: Option<ReferenceOutput<'_>>,
    events: &[EventManifestLink],
    reduced_references: &[ReducedReferenceOutput<'_>],
    reduction_failure_start: Option<usize>,
    diagnostic: Diagnostic,
) -> Result<ExecutionOutcome, Diagnostic> {
    let mut units = vec![UnitRecord {
        id: execution_id.to_owned(),
        kind: "execution".to_owned(),
        state: UnitState::Failed,
        exportable: false,
        diagnostic_code: Some(diagnostic.code().to_owned()),
    }];
    if let Some(failure_start) = reduction_failure_start {
        units.extend(effective.reduced_references().iter().enumerate().map(
            |(index, reduction)| {
                let state = if index < failure_start {
                    UnitState::Complete
                } else {
                    UnitState::Failed
                };
                UnitRecord {
                    id: reduction.id().to_owned(),
                    kind: "reference_reduction".to_owned(),
                    state,
                    exportable: state == UnitState::Complete,
                    diagnostic_code: (state == UnitState::Failed)
                        .then(|| diagnostic.code().to_owned()),
                }
            },
        ));
    } else {
        units.extend(reduced_references.iter().map(|reduction| UnitRecord {
            id: reduction.effective.id().to_owned(),
            kind: "reference_reduction".to_owned(),
            state: reduction.state,
            exportable: reduction.state == UnitState::Complete,
            diagnostic_code:
                (reduction.state == UnitState::Failed).then(|| diagnostic.code().to_owned()),
        }));
    }
    let reduction_artifact_ids = effective
        .reduced_references()
        .iter()
        .enumerate()
        .filter(|(index, _)| reduction_failure_start.is_some_and(|start| *index >= start))
        .flat_map(|(_, reduction)| {
            [
                format!("reference-fasta-{}", reduction.id()),
                format!("reference-map-{}", reduction.id()),
            ]
        })
        .collect::<std::collections::BTreeSet<_>>();
    if reduction_failure_start.is_some() {
        workspace.mark_existing_artifacts_failed(&reduction_artifact_ids);
    }
    let artifacts = workspace.artifacts().to_vec();
    let diagnostics = [diagnostic];
    let manifest = GeneralManifest {
        execution_id,
        state: UnitState::Failed,
        effective,
        configuration_inputs,
        resolution_evidence,
        effective_config,
        reference,
        reduced_references,
        events,
        artifacts: &artifacts,
        units: &units,
        diagnostics: &diagnostics,
    };
    write_validated_general_manifest(&mut workspace, &general_manifest_document(&manifest))?;
    let destination = workspace.publish_with(publisher)?;
    Ok(ExecutionOutcome::Published {
        destination,
        state: UnitState::Failed,
    })
}

fn write_validated_general_manifest(
    workspace: &mut ExecutionWorkspace,
    document: &serde_json::Value,
) -> Result<(), Diagnostic> {
    let bytes = json_document_bytes(document)?;
    let manifest = workspace.validate_general_manifest(Path::new("manifest.json"), &bytes)?;
    workspace.write_general_manifest(&manifest)?;
    Ok(())
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

#[cfg(test)]
mod tests {
    use std::cell::{Cell, RefCell};
    use std::collections::{BTreeMap, BTreeSet};
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::process::Command;
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::{
        ExecutionFault, ExecutionOutcome, ReductionFaultStage, execute_validated_config_with,
    };
    use crate::adapters::{admit_toml_header, parse_declared_config};
    use crate::application::{
        PhaseGenerationInterval, ResolutionDomain, ResolutionRequest, build_effective_config,
        build_effective_config_with_domain, load_and_validate_config,
    };
    use crate::configuration::{DeclaredConfig, ExpectedConfigKind};
    use crate::diagnostics::Diagnostic;
    use crate::publication::{AtomicPublisher, LinuxAtomicPublisher};
    use crate::results::{UnitState, validate_general_manifest_bytes};
    use crate::rules::{ParameterKey, ResolutionContext};

    static NONCE: AtomicU64 = AtomicU64::new(0);

    struct ResolutionEvidenceMutatingPublisher {
        provisional: RefCell<Option<PathBuf>>,
        rename_called: Cell<bool>,
    }

    impl AtomicPublisher for ResolutionEvidenceMutatingPublisher {
        fn before_final_validation(&self, provisional: &Path) -> Result<(), Diagnostic> {
            self.provisional.replace(Some(provisional.to_path_buf()));
            fs::write(
                provisional.join("inputs/resolution/evidence.json"),
                b"{\"forged\":true}\n",
            )
            .expect("mutación controlada de evidencia");
            Ok(())
        }

        fn publish_no_replace(
            &self,
            _provisional: &Path,
            _destination: &Path,
        ) -> Result<(), Diagnostic> {
            self.rename_called.set(true);
            panic!("la mutación de evidencia no debe alcanzar rename")
        }
    }

    #[test]
    fn resolution_evidence_toctou_is_rejected_before_publication() {
        let effective =
            load_and_validate_config(Path::new("tests/fixtures/config/general-seeded.toml"))
                .expect("configuración válida");
        let destination = temporary("resolution-evidence-toctou");
        let publisher = ResolutionEvidenceMutatingPublisher {
            provisional: RefCell::new(None),
            rename_called: Cell::new(false),
        };
        let error = execute_validated_config_with(
            &effective,
            &destination,
            ExecutionFault::None,
            &publisher,
        )
        .expect_err("la revalidación final debe detectar la mutación");
        assert_eq!(error.code(), "GPS-IO-005");
        assert!(!publisher.rename_called.get());
        assert!(!destination.exists());
        let provisional = publisher
            .provisional
            .borrow()
            .clone()
            .expect("ruta provisional observada");
        assert!(provisional.join("INCOMPLETE").is_file());
        fs::remove_dir_all(provisional).expect("limpieza TOCTOU");
    }

    const GLOBAL_LINK_CORPUS_CASES: &[&str] = &[
        "positive-zero-complete",
        "positive-zero-failed",
        "positive-one-complete",
        "positive-one-failed",
        "positive-two-complete",
        "positive-two-failed",
        "zero-complete-reference-digest",
        "zero-complete-authoritative-resource",
        "zero-complete-effective-config-digest",
        "one-complete-reference-digest",
        "one-complete-authoritative-resource",
        "one-complete-effective-config-digest",
        "one-failed-reference-digest",
        "one-failed-authoritative-resource",
        "one-failed-effective-config-digest",
        "two-complete-reference-digest",
        "two-complete-authoritative-resource",
        "two-complete-effective-config-digest",
        "two-failed-reference-digest",
        "two-failed-authoritative-resource",
        "two-failed-effective-config-digest",
        "zero-failed-effective-config-digest",
        "zero-reference-id",
        "zero-reference-kind",
        "zero-reference-path",
        "zero-reference-size",
        "zero-effective-config-id",
        "zero-effective-config-kind",
        "zero-effective-config-path",
        "zero-effective-config-size",
        "zero-event-manifest-id",
        "zero-event-manifest-kind",
        "zero-event-manifest-path",
        "zero-event-manifest-size",
        "zero-complete-event-digest",
    ];

    const MANIFEST_GRAPH_CORPUS_CASES: &[&str] = &[
        "graph-positive-zero-complete",
        "graph-positive-zero-failed",
        "graph-positive-one-complete",
        "graph-positive-one-failed",
        "graph-positive-two-complete",
        "graph-positive-two-failed",
        "graph-positive-multiple-complete",
        "graph-positive-multiple-failed",
        "graph-configuration-copy-digest",
        "graph-configuration-copy-id",
        "graph-configuration-copy-class",
        "graph-configuration-copy-path",
        "graph-configuration-copy-size",
        "graph-configuration-copy-reused",
        "graph-configuration-input-duplicate",
        "graph-configuration-kind-incompatible",
        "graph-general-configuration-missing",
        "graph-resolution-evidence-digest",
        "graph-resolution-evidence-id",
        "graph-resolution-evidence-kind",
        "graph-resolution-evidence-path",
        "graph-resolution-evidence-size",
        "graph-event-unit-crossed",
        "graph-event-unit-orphan",
        "graph-event-unit-missing",
        "graph-execution-unit-crossed",
        "graph-output-artifact-missing",
        "graph-output-wrong-category",
        "graph-artifact-path-duplicate",
        "graph-orphan-event-artifact",
        "graph-failed-configuration-copy-digest",
        "graph-failed-execution-unit-crossed",
        "graph-one-complete-output-missing",
        "graph-one-failed-output-missing",
        "graph-two-complete-output-missing",
        "graph-two-failed-output-missing",
        "graph-two-complete-artifact-path-duplicate",
    ];

    const MANIFEST_DYNAMIC_CORPUS_CASES: &[&str] = &[
        "dynamic-positive-zero-complete",
        "dynamic-positive-zero-failed",
        "dynamic-positive-one-complete",
        "dynamic-positive-one-failed",
        "dynamic-positive-two-complete",
        "dynamic-positive-two-failed",
        "dynamic-positive-percent-complete",
        "dynamic-diagnostic-execution-absent",
        "dynamic-diagnostic-reduction-absent",
        "dynamic-diagnostic-invented",
        "dynamic-diagnostic-category",
        "dynamic-diagnostic-meaning",
        "dynamic-diagnostic-required-context",
        "dynamic-diagnostic-duplicate",
        "dynamic-diagnostic-orphan",
        "dynamic-reference-sequence-duplicate",
        "dynamic-reference-total",
        "dynamic-reference-unknown-count",
        "dynamic-reduction-partition-overlap",
        "dynamic-reduction-partition-missing",
        "dynamic-reduction-selection-order",
        "dynamic-reduction-sequence-missing",
        "dynamic-reduction-out-of-range",
        "dynamic-reduction-overlap",
        "dynamic-reduction-global-sum",
        "dynamic-reduction-per-sequence-sum",
        "dynamic-reduction-eligible",
        "dynamic-reduction-declared-target",
        "dynamic-reduction-derived-target",
        "dynamic-reduction-region",
        "dynamic-reduction-residual",
        "dynamic-reduction-derived-class",
        "dynamic-reduction-derived-ploidy",
        "dynamic-reduction-percent-target",
    ];

    const MANIFEST_EXHAUSTIVE_CORPUS_CASES: &[&str] = &[
        "exhaustive-positive-percent-above-integer-boundary",
        "exhaustive-positive-percent-exponent-notation",
        "exhaustive-positive-percent-below-integer-boundary",
        "exhaustive-positive-non-utf8-header",
        "exhaustive-positive-multiple-events-normalizations",
        "exhaustive-positive-two-reductions",
        "exhaustive-reduction-root-derived",
        "exhaustive-reduction-subseed-derived",
        "exhaustive-derived-reference-order",
        "exhaustive-header-token",
        "exhaustive-header-description",
        "exhaustive-header-description-non-utf8",
        "exhaustive-normalization-sum-before",
        "exhaustive-normalization-delta-before",
        "exhaustive-normalization-before",
        "exhaustive-normalization-after",
        "exhaustive-normalization-sum-after",
        "exhaustive-normalization-delta-after",
        "exhaustive-normalization-decision",
        "exhaustive-normalization-order",
        "exhaustive-normalization-duplicate-value",
        "exhaustive-normalization-duplicate-parameter",
        "exhaustive-event-general-subseed",
        "exhaustive-event-effective-subseed",
        "exhaustive-event-manifest-subseed",
        "exhaustive-event-correlated-derived-subseed",
    ];

    const MANIFEST_R17_CORPUS_CASES: &[&str] = &[
        "r17-positive-percent-separators",
        "r17-positive-percent-plus-exponent",
        "r17-positive-percent-high-precision",
        "r17-positive-fasta-non-utf8",
        "r17-positive-effective-multiple-scopes",
        "r17-positive-two-reductions",
        "r17-percent-leading-separator",
        "r17-percent-trailing-separator",
        "r17-percent-double-separator",
        "r17-percent-separator-before-dot",
        "r17-percent-separator-after-dot",
        "r17-percent-missing-integer",
        "r17-percent-missing-fraction",
        "r17-percent-leading-zero-integer",
        "r17-percent-leading-zero-float",
        "r17-percent-signed-leading-zero",
        "r17-percent-empty-exponent",
        "r17-percent-empty-signed-exponent",
        "r17-percent-leading-exponent-separator",
        "r17-percent-trailing-exponent-separator",
        "r17-fasta-correlated-id-header",
        "r17-fasta-correlated-description",
        "r17-fasta-correlated-length",
        "r17-fasta-correlated-unknown-count",
        "r17-fasta-correlated-order",
        "r17-fasta-intermediate-symlink",
        "r17-effective-duplicate-identical",
        "r17-effective-duplicate-contradictory",
        "r17-effective-reordered",
        "r17-effective-parameter-cross",
        "r17-effective-value-cross",
        "r17-effective-winner-cross",
        "r17-effective-selectors-cross",
        "r17-effective-provenance-path",
        "r17-effective-provenance-location",
        "r17-effective-automatic-missing",
        "r17-effective-correlated-ploidy",
        "r17-effective-source-path",
        "r17-effective-controlled-duplicate",
        "r17-effective-controlled-count",
    ];

    const MANIFEST_R18_CORPUS_CASES: &[&str] = &[
        "r18-effective-dominated-base-copy",
        "r18-effective-foreign-object",
    ];

    const MANIFEST_R19_R21_CORPUS_CASES: &[&str] = &[
        "r19-biological-parameter-technical-context",
        "r21-orphan-position-with-individual",
    ];

    const MANIFEST_R22_CORPUS_CASES: &[&str] = &[
        "r22-coverage-locus",
        "r22-recombination-region",
        "r22-ploidy-locus",
        "r22-reference-object",
    ];

    const MANIFEST_R23_R24_CORPUS_CASES: &[&str] = &[
        "r23-complete-missing-phase-domain-proof",
        "r23-failed-missing-phase-domain-proof",
        "r24-complete-scenario-copy-id",
        "r24-failed-scenario-copy-id",
    ];

    const MANIFEST_R25_CORPUS_CASES: &[&str] = &[
        "r25-complete-correlated-domain-provider-interval",
        "r25-failed-correlated-domain-provider-interval",
        "r25-correlated-global-output-to-event",
    ];

    const MANIFEST_R26_CORPUS_CASES: &[&str] = &[
        "r26-complete-duplicate-phase-across-scenarios",
        "r26-failed-duplicate-phase-across-scenarios",
    ];

    const MANIFEST_R27_CORPUS_CASES: &[&str] =
        &["r27-complete-parameter-order", "r27-failed-parameter-order"];

    fn temporary(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gps-execution-unit-{}-{}-{label}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/reference_reduction")
            .join(name)
    }

    fn reduction_fixture(count: usize, label: &str) -> (PathBuf, PathBuf) {
        assert!((1..=2).contains(&count));
        let root = temporary(label);
        if root.exists() {
            fs::remove_dir_all(&root).expect("stale reduction corpus cleanup");
        }
        fs::create_dir_all(&root).expect("reduction fixture root");
        fs::copy(fixture("source.fa"), root.join("source.fa")).expect("source fixture copy");
        let mut config = fs::read_to_string(fixture("general.toml")).expect("fixture TOML");
        if count == 2 {
            config.push_str(
                r#"

[[reduced_references]]
id = "reference-tiny"
source_reference_id = "reference-main"
region_size_bases = 4
target_bases = 7
"#,
            );
        }
        let path = root.join("general.toml");
        fs::write(&path, config).expect("write reduction fixture");
        (root, path)
    }

    fn python_schema_accepts(path: &Path, schema: &str) -> bool {
        let schema_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("schemas");
        let base_uri = format!("file://{}/", schema_root.display());
        Command::new("python3")
            .args(["-m", "jsonschema", "--base-uri", &base_uri, "-i"])
            .arg(path)
            .arg(schema)
            .current_dir(schema_root)
            .output()
            .expect("python-jsonschema")
            .status
            .success()
    }

    fn reduction_unit<'a>(manifest: &'a serde_json::Value, id: &str) -> &'a serde_json::Value {
        manifest["units"]
            .as_array()
            .expect("units")
            .iter()
            .find(|unit| unit["kind"] == "reference_reduction" && unit["id"] == id)
            .expect("reduction unit")
    }

    fn artifact<'a>(manifest: &'a serde_json::Value, id: &str) -> Option<&'a serde_json::Value> {
        manifest["artifacts"]
            .as_array()
            .expect("artifacts")
            .iter()
            .find(|artifact| artifact["id"] == id)
    }

    #[test]
    fn internal_execution_faults_cover_failed_closure_and_interruption() {
        let seeded = load_and_validate_config(
            &Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/config/general-seeded.toml"),
        )
        .expect("seeded fixture");
        for (label, fault) in [
            ("configuration-copy", ExecutionFault::AfterConfigurationCopy),
            ("effective-config", ExecutionFault::AfterEffectiveConfig),
        ] {
            let destination = temporary(label);
            let outcome =
                execute_validated_config_with(&seeded, &destination, fault, &LinuxAtomicPublisher)
                    .expect("controlled failed closure");
            assert!(matches!(
                outcome,
                ExecutionOutcome::Published {
                    state: UnitState::Failed,
                    ..
                }
            ));
            let manifest: serde_json::Value = serde_json::from_slice(
                &fs::read(destination.join("manifest.json")).expect("failed manifest"),
            )
            .expect("failed manifest JSON");
            assert_eq!(manifest["state"], "failed");
            assert!(!destination.join("INCOMPLETE").exists());
            fs::remove_dir_all(destination).expect("failed result cleanup");
        }

        let destination = temporary("interrupted");
        let outcome = execute_validated_config_with(
            &seeded,
            &destination,
            ExecutionFault::InterruptedAfterWorkspace,
            &LinuxAtomicPublisher,
        )
        .expect("controlled interruption");
        let ExecutionOutcome::Interrupted { provisional } = outcome else {
            panic!("interruption must keep a provisional workspace")
        };
        assert!(provisional.join("INCOMPLETE").is_file());
        assert!(!destination.exists());
        fs::remove_dir_all(provisional).expect("interrupted cleanup");
    }

    #[test]
    fn reduction_fault_matrix_preserves_prior_complete_and_fails_target_and_later_units() {
        let stages = [
            ("after-fasta", ReductionFaultStage::Fasta),
            ("after-map", ReductionFaultStage::Map),
            ("after-metadata", ReductionFaultStage::Metadata),
        ];
        for count in 1..=2 {
            for target_index in 0..count {
                for (stage_label, stage) in stages {
                    let label = format!("matrix-{count}-{target_index}-{stage_label}");
                    let (root, config) = reduction_fixture(count, &label);
                    let effective =
                        load_and_validate_config(&config).expect("reduction fixture configuration");
                    let destination = root.join("output");
                    let outcome = execute_validated_config_with(
                        &effective,
                        &destination,
                        ExecutionFault::AtReduction {
                            reduction_index: target_index,
                            stage,
                        },
                        &LinuxAtomicPublisher,
                    )
                    .expect("controlled reduction failure");
                    assert!(matches!(
                        outcome,
                        ExecutionOutcome::Published {
                            state: UnitState::Failed,
                            ..
                        }
                    ));

                    let manifest_path = destination.join("manifest.json");
                    let manifest_bytes = fs::read(&manifest_path).expect("failed manifest");
                    let manifest: serde_json::Value =
                        serde_json::from_slice(&manifest_bytes).expect("failed manifest JSON");
                    assert_eq!(manifest["state"], "failed");
                    validate_general_manifest_bytes(&manifest_path, &manifest_bytes)
                        .expect("Rust validates failed manifest after reloading bytes");
                    assert!(
                        python_schema_accepts(&manifest_path, "manifest-general-v1.schema.json"),
                        "failed manifest schema: {label}"
                    );

                    let reduction_ids = ["reference-small", "reference-tiny"];
                    for (index, id) in reduction_ids[..count].iter().enumerate() {
                        let expected_state = if index < target_index {
                            "complete"
                        } else {
                            "failed"
                        };
                        let unit = reduction_unit(&manifest, id);
                        assert_eq!(unit["state"], expected_state, "unit {id}: {label}");
                        assert_eq!(
                            unit["exportable"],
                            expected_state == "complete",
                            "exportability {id}: {label}"
                        );
                        assert_eq!(
                            manifest["references"]
                                .as_array()
                                .expect("references")
                                .iter()
                                .any(|reference| reference["id"] == *id),
                            expected_state == "complete",
                            "derived reference presence {id}: {label}"
                        );
                        for artifact_id in [
                            format!("reference-fasta-{id}"),
                            format!("reference-map-{id}"),
                        ] {
                            if let Some(record) = artifact(&manifest, &artifact_id) {
                                assert_eq!(
                                    record["state"], expected_state,
                                    "artifact state {artifact_id}: {label}"
                                );
                                assert_eq!(
                                    record["validated"],
                                    expected_state == "complete",
                                    "artifact validation {artifact_id}: {label}"
                                );
                            } else {
                                assert!(
                                    index >= target_index,
                                    "a prior complete reduction cannot lose {artifact_id}: {label}"
                                );
                            }
                        }
                    }

                    if stage == ReductionFaultStage::Metadata {
                        let effective_path = destination.join("effective-config.json");
                        let effective_document: serde_json::Value = serde_json::from_slice(
                            &fs::read(&effective_path).expect("effective configuration"),
                        )
                        .expect("effective configuration JSON");
                        assert!(
                            python_schema_accepts(
                                &effective_path,
                                "effective-config-v1.schema.json"
                            ),
                            "effective configuration schema: {label}"
                        );
                        for (index, id) in reduction_ids[..count].iter().enumerate() {
                            let expected_state = if index < target_index {
                                "complete"
                            } else {
                                "failed"
                            };
                            let reduction = effective_document["reduced_references"]
                                .as_array()
                                .expect("effective reductions")
                                .iter()
                                .find(|reduction| reduction["derived_reference_id"] == *id)
                                .expect("effective reduction");
                            assert_eq!(
                                reduction["derived_fasta"]["state"], expected_state,
                                "effective FASTA state for {id}: {label}"
                            );
                            assert_eq!(
                                reduction["map"]["state"], expected_state,
                                "effective map state for {id}: {label}"
                            );
                            assert_eq!(
                                reduction_unit(&manifest, id)["state"],
                                expected_state,
                                "effective config and manifest unit disagree for {id}: {label}"
                            );
                        }
                    } else {
                        assert!(
                            !destination.join("effective-config.json").exists(),
                            "metadata must not be claimed before its boundary: {label}"
                        );
                    }
                    assert!(!destination.join("INCOMPLETE").exists());
                    fs::remove_dir_all(root).expect("reduction matrix cleanup");
                }
            }
        }
    }

    fn write_manifest_case(root: &Path, name: &str, document: &serde_json::Value) -> PathBuf {
        let state_directory = match document["state"].as_str() {
            Some("complete") => "complete",
            Some("failed") => "failed",
            _ => "complete",
        };
        let path = root.join(state_directory).join(format!("{name}.json"));
        fs::write(
            &path,
            serde_json::to_vec_pretty(document).expect("manifest serialization"),
        )
        .expect("manifest corpus write");
        path
    }

    fn independent_text<'a>(value: &'a serde_json::Value, field: &str) -> Result<&'a str, String> {
        value[field]
            .as_str()
            .ok_or_else(|| format!("{field} is not text"))
    }

    fn independent_string_set<'a>(
        value: &'a serde_json::Value,
        context: &str,
    ) -> Result<BTreeSet<&'a str>, String> {
        let values = value
            .as_array()
            .ok_or_else(|| format!("{context} is not an array"))?;
        let mut result = BTreeSet::new();
        for item in values {
            let item = item
                .as_str()
                .ok_or_else(|| format!("{context} contains a non-string"))?;
            if !result.insert(item) {
                return Err(format!("{context} contains duplicate {item}"));
            }
        }
        Ok(result)
    }

    type IndependentDynamicInvariant = fn(&serde_json::Value) -> Result<(), String>;

    const INDEPENDENT_DYNAMIC_INVARIANTS: &[(&str, IndependentDynamicInvariant)] = &[
        (
            "diagnostics-v2-and-failed-unit-edges",
            independent_diagnostic_graph,
        ),
        (
            "reference-catalog-and-reduction-equations",
            independent_reference_reduction_graph,
        ),
        (
            "normalization-neumaier-and-canonical-order",
            independent_normalization_graph,
        ),
    ];

    fn independent_diagnostic_graph(document: &serde_json::Value) -> Result<(), String> {
        let catalog = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/catalogs/diagnostics-v2.tsv"
        ));
        let mut contracts = BTreeMap::<&str, (&str, &str, BTreeSet<&str>)>::new();
        for line in catalog
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
        {
            let fields = line.split('|').collect::<Vec<_>>();
            if fields.len() != 6 {
                return Err("independent diagnostic catalog field count".to_owned());
            }
            let required = fields[3].split(',').collect::<BTreeSet<_>>();
            if contracts
                .insert(fields[0], (fields[1], fields[2], required))
                .is_some()
            {
                return Err(format!("independent duplicate catalog code {}", fields[0]));
            }
        }
        let mut embedded = BTreeMap::<&str, &serde_json::Value>::new();
        for diagnostic in document["diagnostics"]
            .as_array()
            .ok_or_else(|| "independent diagnostics is not an array".to_owned())?
        {
            let code = independent_text(diagnostic, "code")?;
            let (category, meaning, required) = contracts
                .get(code)
                .ok_or_else(|| format!("independent unknown diagnostic {code}"))?;
            if independent_text(diagnostic, "category")? != *category
                || independent_text(diagnostic, "message")? != *meaning
            {
                return Err(format!("independent diagnostic contract mismatch {code}"));
            }
            let context = diagnostic["context"]
                .as_object()
                .ok_or_else(|| format!("independent diagnostic context missing {code}"))?;
            if required.iter().any(|field| {
                context
                    .get(*field)
                    .and_then(serde_json::Value::as_str)
                    .is_none_or(str::is_empty)
            }) {
                return Err(format!("independent diagnostic context incomplete {code}"));
            }
            if embedded.insert(code, diagnostic).is_some() {
                return Err(format!("independent duplicate embedded diagnostic {code}"));
            }
        }
        let mut referenced = BTreeSet::new();
        for unit in document["units"]
            .as_array()
            .ok_or_else(|| "independent units is not an array".to_owned())?
        {
            match independent_text(unit, "state")? {
                "complete" if unit["diagnostic_code"].is_null() => {}
                "failed" => {
                    let code = independent_text(unit, "diagnostic_code")?;
                    if !embedded.contains_key(code) {
                        return Err(format!("independent unresolved unit diagnostic {code}"));
                    }
                    referenced.insert(code);
                }
                _ => return Err("independent unit diagnostic state mismatch".to_owned()),
            }
        }
        if embedded.keys().copied().collect::<BTreeSet<_>>() != referenced {
            return Err("independent orphan diagnostic or missing unit edge".to_owned());
        }
        Ok(())
    }

    #[derive(Clone)]
    struct IndependentSequence {
        class: String,
        length: u64,
        ploidy: u64,
    }

    #[derive(Clone)]
    struct IndependentReference {
        total: u64,
        order: Vec<String>,
        sequences: BTreeMap<String, IndependentSequence>,
    }

    fn independent_reference_reduction_graph(document: &serde_json::Value) -> Result<(), String> {
        let mut references = BTreeMap::<String, IndependentReference>::new();
        for reference in document["references"]
            .as_array()
            .ok_or_else(|| "independent references is not an array".to_owned())?
        {
            let reference_id = independent_text(reference, "id")?.to_owned();
            let mut sequences = BTreeMap::new();
            let mut order = Vec::new();
            let mut total = 0_u64;
            for sequence in reference["sequences"]
                .as_array()
                .ok_or_else(|| format!("independent sequences missing {reference_id}"))?
            {
                let id = independent_text(sequence, "id")?.to_owned();
                independent_sequence_header(sequence, &id)?;
                let length = sequence["length"]
                    .as_u64()
                    .ok_or_else(|| format!("independent sequence length missing {id}"))?;
                if sequence["unknown_base_count"]
                    .as_u64()
                    .is_none_or(|count| count > length)
                {
                    return Err(format!("independent unknown bases exceed length {id}"));
                }
                total = total
                    .checked_add(length)
                    .ok_or_else(|| "independent reference total overflow".to_owned())?;
                let summary = IndependentSequence {
                    class: independent_text(sequence, "class")?.to_owned(),
                    length,
                    ploidy: sequence["ploidy"]
                        .as_u64()
                        .ok_or_else(|| format!("independent ploidy missing {id}"))?,
                };
                if sequences.insert(id.clone(), summary).is_some() {
                    return Err(format!(
                        "independent duplicate sequence {reference_id}/{id}"
                    ));
                }
                order.push(id);
            }
            if reference["total_genome_size"].as_u64() != Some(total) {
                return Err(format!(
                    "independent reference total mismatch {reference_id}"
                ));
            }
            if references
                .insert(
                    reference_id.clone(),
                    IndependentReference {
                        total,
                        order,
                        sequences,
                    },
                )
                .is_some()
            {
                return Err(format!("independent duplicate reference {reference_id}"));
            }
        }

        let reduction_states = document["units"]
            .as_array()
            .ok_or_else(|| "independent units missing".to_owned())?
            .iter()
            .filter(|unit| unit["kind"] == "reference_reduction")
            .map(|unit| {
                Ok((
                    independent_text(unit, "id")?.to_owned(),
                    independent_text(unit, "state")?.to_owned(),
                ))
            })
            .collect::<Result<BTreeMap<_, _>, String>>()?;
        for reduction in document["reference_reductions"]
            .as_array()
            .ok_or_else(|| "independent reductions is not an array".to_owned())?
        {
            let source_id = independent_text(reduction, "source_reference_id")?;
            let derived_id = independent_text(reduction, "derived_reference_id")?;
            let source = references
                .get(source_id)
                .ok_or_else(|| format!("independent source missing {source_id}"))?;
            let eligible = reduction["eligible_bases"]
                .as_u64()
                .ok_or_else(|| format!("independent eligible missing {derived_id}"))?;
            let region = reduction["region_size_bases"]
                .as_u64()
                .ok_or_else(|| format!("independent region missing {derived_id}"))?;
            let target = reduction["target"]["derived_bases"]
                .as_u64()
                .ok_or_else(|| format!("independent target missing {derived_id}"))?;
            if eligible != source.total || target == 0 || target >= eligible || region == 0 {
                return Err(format!("independent eligible/target mismatch {derived_id}"));
            }
            let declared = independent_text(&reduction["target"], "declared_value")?;
            let expected_target = match independent_text(&reduction["target"], "mode")? {
                "target_bases" => declared
                    .replace('_', "")
                    .parse::<u64>()
                    .map_err(|_| "independent target literal".to_owned())?,
                "reduction_percent" => independent_retained_after_percent(declared, eligible)?,
                other => return Err(format!("independent target mode {other}")),
            };
            if expected_target != target {
                return Err(format!("independent declared target mismatch {derived_id}"));
            }
            let included = independent_string_set(
                &reduction["included_sequence_ids"],
                "independent included",
            )?
            .into_iter()
            .map(str::to_owned)
            .collect::<BTreeSet<_>>();
            let omitted =
                independent_string_set(&reduction["omitted_sequence_ids"], "independent omitted")?
                    .into_iter()
                    .map(str::to_owned)
                    .collect::<BTreeSet<_>>();
            let source_ids = source.sequences.keys().cloned().collect::<BTreeSet<_>>();
            if !included.is_disjoint(&omitted)
                || included.union(&omitted).cloned().collect::<BTreeSet<_>>() != source_ids
            {
                return Err(format!(
                    "independent source partition mismatch {derived_id}"
                ));
            }

            let selection = reduction["selection_order"]
                .as_array()
                .ok_or_else(|| format!("independent selection missing {derived_id}"))?;
            let full_count = target / region;
            let residual = target % region;
            let expected_count = full_count + u64::from(residual > 0);
            if u64::try_from(selection.len()).ok() != Some(expected_count)
                || reduction["has_residual"] != serde_json::json!(residual > 0)
            {
                return Err(format!("independent residual/count mismatch {derived_id}"));
            }
            let mut selected_ids = BTreeSet::new();
            let mut selected_lengths = BTreeMap::<String, u64>::new();
            let mut intervals = BTreeMap::<String, Vec<(u64, u64)>>::new();
            let mut total = 0_u64;
            for (index, segment) in selection.iter().enumerate() {
                let order = segment["selection_order"]
                    .as_u64()
                    .ok_or_else(|| format!("independent selection order {derived_id}"))?;
                let expected_order = u64::try_from(index).map_err(|_| "index overflow")? + 1;
                if order != expected_order {
                    return Err(format!("independent nonconsecutive selection {derived_id}"));
                }
                let sequence_id = independent_text(segment, "sequence_id")?;
                let sequence = source.sequences.get(sequence_id).ok_or_else(|| {
                    format!("independent selected sequence missing {sequence_id}")
                })?;
                let start = segment["start"]
                    .as_u64()
                    .ok_or_else(|| "independent selection start".to_owned())?;
                let end = segment["end"]
                    .as_u64()
                    .ok_or_else(|| "independent selection end".to_owned())?;
                if start == 0 || start > end || end > sequence.length {
                    return Err(format!("independent interval out of range {sequence_id}"));
                }
                let length = end - start + 1;
                let expected_length = if order <= full_count {
                    region
                } else if residual > 0 && order == expected_count {
                    residual
                } else {
                    return Err(format!("independent residual position {derived_id}"));
                };
                if length != expected_length {
                    return Err(format!("independent segment length {derived_id}"));
                }
                total = total
                    .checked_add(length)
                    .ok_or_else(|| "independent selected total overflow".to_owned())?;
                *selected_lengths.entry(sequence_id.to_owned()).or_default() += length;
                intervals
                    .entry(sequence_id.to_owned())
                    .or_default()
                    .push((start, end));
                selected_ids.insert(sequence_id.to_owned());
            }
            for values in intervals.values_mut() {
                values.sort_unstable();
                if values.windows(2).any(|pair| pair[1].0 <= pair[0].1) {
                    return Err(format!("independent overlapping selection {derived_id}"));
                }
            }
            if total != target || selected_ids != included {
                return Err(format!("independent selected target mismatch {derived_id}"));
            }
            match reduction_states.get(derived_id).map(String::as_str) {
                Some("complete") => {
                    let derived = references
                        .get(derived_id)
                        .ok_or_else(|| format!("independent derived missing {derived_id}"))?;
                    if derived.total != target
                        || derived.sequences.keys().cloned().collect::<BTreeSet<_>>() != included
                        || derived.order
                            != source
                                .order
                                .iter()
                                .filter(|id| included.contains(*id))
                                .cloned()
                                .collect::<Vec<_>>()
                    {
                        return Err(format!("independent derived catalog mismatch {derived_id}"));
                    }
                    for (id, sequence) in &derived.sequences {
                        let source_sequence = &source.sequences[id];
                        if selected_lengths.get(id) != Some(&sequence.length)
                            || sequence.class != source_sequence.class
                            || sequence.ploidy != source_sequence.ploidy
                        {
                            return Err(format!("independent derived sequence mismatch {id}"));
                        }
                    }
                }
                Some("failed") if references.contains_key(derived_id) => {
                    return Err(format!("independent failed derived published {derived_id}"));
                }
                Some("failed") => {}
                _ => return Err(format!("independent reduction unit missing {derived_id}")),
            }
        }
        Ok(())
    }

    fn independent_sequence_header(
        sequence: &serde_json::Value,
        sequence_id: &str,
    ) -> Result<(), String> {
        let raw = independent_reversible_bytes(&sequence["header"])?;
        let first = raw
            .iter()
            .position(|byte| !byte.is_ascii_whitespace())
            .ok_or_else(|| "independent empty FASTA header".to_owned())?;
        let end = raw[first..]
            .iter()
            .position(u8::is_ascii_whitespace)
            .map_or(raw.len(), |offset| first + offset);
        let token = std::str::from_utf8(&raw[first..end])
            .map_err(|_| "independent non-UTF8 FASTA ID".to_owned())?;
        if token != sequence_id {
            return Err(format!("independent header token mismatch {sequence_id}"));
        }
        let observed_description = raw.get(end.saturating_add(1)..).unwrap_or_default();
        if observed_description.is_empty() {
            if !sequence["description"].is_null() {
                return Err(format!("independent unexpected description {sequence_id}"));
            }
        } else if sequence["description"].is_null()
            || independent_reversible_bytes(&sequence["description"])? != observed_description
        {
            return Err(format!("independent description mismatch {sequence_id}"));
        }
        Ok(())
    }

    fn independent_reversible_bytes(value: &serde_json::Value) -> Result<Vec<u8>, String> {
        match independent_text(value, "encoding")? {
            "utf-8" => Ok(independent_text(value, "value")?.as_bytes().to_vec()),
            "bytes-hex" => {
                let encoded = independent_text(value, "value")?.as_bytes();
                if encoded.len() % 2 != 0 {
                    return Err("independent odd reversible hex".to_owned());
                }
                encoded
                    .chunks_exact(2)
                    .map(|pair| {
                        let decode = |byte| match byte {
                            b'0'..=b'9' => Ok(byte - b'0'),
                            b'a'..=b'f' => Ok(byte - b'a' + 10),
                            _ => Err("independent invalid reversible hex".to_owned()),
                        };
                        Ok(decode(pair[0])? * 16 + decode(pair[1])?)
                    })
                    .collect()
            }
            other => Err(format!("independent reversible encoding {other}")),
        }
    }

    fn independent_normalization_graph(document: &serde_json::Value) -> Result<(), String> {
        let collections = document["normalizations"]
            .as_array()
            .ok_or_else(|| "independent normalizations missing".to_owned())?;
        let mut parameters = BTreeSet::new();
        for collection in collections {
            let parameter = independent_text(collection, "parameter")?;
            if !parameters.insert(parameter) {
                return Err(format!("independent duplicate normalization {parameter}"));
            }
            let tolerance = collection["tolerance"]
                .as_f64()
                .ok_or_else(|| "independent normalization tolerance".to_owned())?;
            if !(0.0..=1e-6).contains(&tolerance) {
                return Err("independent normalization tolerance range".to_owned());
            }
            let values = collection["values"]
                .as_array()
                .ok_or_else(|| "independent normalization values".to_owned())?;
            let mut previous: Option<(&str, &str)> = None;
            let mut before = Vec::new();
            let mut after = Vec::new();
            for value in values {
                let identity = (
                    independent_text(value, "id")?,
                    independent_text(value, "category")?,
                );
                if previous.is_some_and(|prior| prior >= identity) {
                    return Err("independent normalization order/duplicate".to_owned());
                }
                previous = Some(identity);
                before.push(
                    value["before"]
                        .as_f64()
                        .ok_or_else(|| "independent normalization before".to_owned())?,
                );
                after.push(
                    value["after"]
                        .as_f64()
                        .ok_or_else(|| "independent normalization after".to_owned())?,
                );
            }
            let sum_before = independent_neumaier(&before);
            let delta_before = (sum_before - 1.0).abs();
            if collection["sum_before"].as_f64().map(f64::to_bits) != Some(sum_before.to_bits())
                || collection["delta_before"].as_f64().map(f64::to_bits)
                    != Some(delta_before.to_bits())
                || delta_before > tolerance
                || sum_before <= 0.0
            {
                return Err("independent normalization pre-sum".to_owned());
            }
            let expected_after = before
                .iter()
                .map(|value| *value / sum_before)
                .collect::<Vec<_>>();
            if after
                .iter()
                .zip(&expected_after)
                .any(|(observed, expected)| observed.to_bits() != expected.to_bits())
            {
                return Err("independent normalization values-after".to_owned());
            }
            let sum_after = independent_neumaier(&expected_after);
            let delta_after = (sum_after - 1.0).abs();
            if collection["sum_after"].as_f64().map(f64::to_bits) != Some(sum_after.to_bits())
                || collection["delta_after"].as_f64().map(f64::to_bits)
                    != Some(delta_after.to_bits())
                || delta_after > tolerance
            {
                return Err("independent normalization post-sum".to_owned());
            }
            let changed = before
                .iter()
                .zip(&expected_after)
                .any(|(raw, normalized)| raw.to_bits() != normalized.to_bits());
            let decision = if changed {
                "normalized"
            } else {
                "accepted_without_change"
            };
            if independent_text(collection, "decision")? != decision {
                return Err("independent normalization decision".to_owned());
            }
        }
        Ok(())
    }

    fn independent_neumaier(values: &[f64]) -> f64 {
        let (sum, correction) =
            values
                .iter()
                .fold((0.0_f64, 0.0_f64), |(sum, correction), value| {
                    let next = sum + value;
                    let adjustment = if sum.abs() >= value.abs() {
                        (sum - next) + value
                    } else {
                        (value - next) + sum
                    };
                    (next, correction + adjustment)
                });
        sum + correction
    }

    fn independent_retained_after_percent(literal: &str, eligible: u64) -> Result<u64, String> {
        let (negative, integer, fraction, exponent) = independent_toml_decimal_parts(literal)?;
        if negative {
            return Err("independent negative percentage".to_owned());
        }
        let coefficient = format!("{integer}{fraction}");
        let coefficient = coefficient.trim_start_matches('0');
        let coefficient = if coefficient.is_empty() {
            "0"
        } else {
            coefficient
        };
        let fraction_digits = i64::try_from(fraction.len())
            .map_err(|_| "independent percentage scale overflow".to_owned())?;
        let power = exponent
            .checked_sub(fraction_digits)
            .ok_or_else(|| "independent percentage scale overflow".to_owned())?;
        if coefficient == "0" || independent_decimal_at_least_hundred(coefficient, power) {
            return Err("independent percentage outside (0,100)".to_owned());
        }

        // Independientemente del tipo de producción, se calcula primero
        // ceil(S*p/100) mediante aritmética larga decimal y luego S-removed.
        let mut product = independent_multiply_decimal(coefficient, eligible);
        let divisor_places = 2_i64
            .checked_sub(power)
            .ok_or_else(|| "independent percentage divisor overflow".to_owned())?;
        let removed = if divisor_places <= 0 {
            let zeros = usize::try_from(divisor_places.unsigned_abs())
                .map_err(|_| "independent percentage expansion overflow".to_owned())?;
            product.extend(std::iter::repeat_n('0', zeros));
            product
                .parse::<u64>()
                .map_err(|_| "independent removed bases overflow".to_owned())?
        } else {
            let places = usize::try_from(divisor_places)
                .map_err(|_| "independent percentage divisor overflow".to_owned())?;
            if places >= product.len() {
                u64::from(product.bytes().any(|byte| byte != b'0'))
            } else {
                let split = product.len() - places;
                let quotient = product[..split]
                    .parse::<u64>()
                    .map_err(|_| "independent quotient overflow".to_owned())?;
                quotient
                    .checked_add(u64::from(product[split..].bytes().any(|byte| byte != b'0')))
                    .ok_or_else(|| "independent quotient overflow".to_owned())?
            }
        };
        eligible
            .checked_sub(removed)
            .ok_or_else(|| "independent percentage removes more than eligible".to_owned())
    }

    fn independent_toml_decimal_parts(
        literal: &str,
    ) -> Result<(bool, String, String, i64), String> {
        let bytes = literal.as_bytes();
        let mut position = 0_usize;
        let negative = match bytes.first() {
            Some(b'+') => {
                position = 1;
                false
            }
            Some(b'-') => {
                position = 1;
                true
            }
            _ => false,
        };
        let integer = independent_decimal_digits(bytes, &mut position)?;
        if integer.len() > 1 && integer.starts_with('0') {
            return Err("independent leading zero in percentage".to_owned());
        }
        let fraction = if bytes.get(position) == Some(&b'.') {
            position += 1;
            independent_decimal_digits(bytes, &mut position)?
        } else {
            String::new()
        };
        let exponent = if matches!(bytes.get(position), Some(b'e' | b'E')) {
            position += 1;
            let negative_exponent = match bytes.get(position) {
                Some(b'+') => {
                    position += 1;
                    false
                }
                Some(b'-') => {
                    position += 1;
                    true
                }
                _ => false,
            };
            let digits = independent_decimal_digits(bytes, &mut position)?;
            let magnitude = digits
                .parse::<i64>()
                .map_err(|_| "independent percentage exponent overflow".to_owned())?;
            if negative_exponent {
                magnitude
                    .checked_neg()
                    .ok_or_else(|| "independent percentage exponent overflow".to_owned())?
            } else {
                magnitude
            }
        } else {
            0
        };
        if position != bytes.len() {
            return Err("independent trailing percentage syntax".to_owned());
        }
        Ok((negative, integer, fraction, exponent))
    }

    fn independent_decimal_digits(bytes: &[u8], position: &mut usize) -> Result<String, String> {
        let mut digits = String::new();
        let mut expecting_digit = true;
        while let Some(byte) = bytes.get(*position).copied() {
            match byte {
                b'0'..=b'9' => {
                    digits.push(char::from(byte));
                    expecting_digit = false;
                    *position += 1;
                }
                b'_' if !expecting_digit
                    && bytes.get(*position + 1).is_some_and(u8::is_ascii_digit) =>
                {
                    expecting_digit = true;
                    *position += 1;
                }
                _ => break,
            }
        }
        if digits.is_empty() || expecting_digit {
            Err("independent malformed decimal digit group".to_owned())
        } else {
            Ok(digits)
        }
    }

    fn independent_decimal_at_least_hundred(coefficient: &str, power: i64) -> bool {
        let Ok(digits) = i64::try_from(coefficient.len()) else {
            return true;
        };
        let Some(integer_digits) = digits.checked_add(power) else {
            return power.is_positive();
        };
        match integer_digits.cmp(&3) {
            std::cmp::Ordering::Less => false,
            std::cmp::Ordering::Greater => true,
            std::cmp::Ordering::Equal => {
                let mut prefix = coefficient.bytes().take(3).collect::<Vec<_>>();
                prefix.resize(3, b'0');
                prefix.as_slice() >= b"100".as_slice()
            }
        }
    }

    fn independent_multiply_decimal(digits: &str, multiplier: u64) -> String {
        let mut output = Vec::with_capacity(digits.len() + 20);
        let mut carry = 0_u128;
        for byte in digits.bytes().rev() {
            let value = u128::from(byte - b'0') * u128::from(multiplier) + carry;
            output.push(b'0' + u8::try_from(value % 10).expect("decimal digit"));
            carry = value / 10;
        }
        while carry != 0 {
            output.push(b'0' + u8::try_from(carry % 10).expect("decimal digit"));
            carry /= 10;
        }
        output.reverse();
        String::from_utf8(output).expect("ASCII decimal product")
    }

    fn independent_semantic_manifest_oracle(path: &Path, bytes: &[u8]) -> Result<(), String> {
        let document: serde_json::Value =
            serde_json::from_slice(bytes).map_err(|error| error.to_string())?;
        for (name, invariant) in INDEPENDENT_DYNAMIC_INVARIANTS {
            invariant(&document).map_err(|detail| format!("{name}: {detail}"))?;
        }
        independent_linked_document_graph(path, &document)?;
        let array = |field: &str| {
            document[field]
                .as_array()
                .ok_or_else(|| format!("{field} is not an array"))
        };
        let mut artifacts = BTreeMap::<&str, &serde_json::Value>::new();
        let mut artifact_paths = BTreeMap::<&str, &str>::new();
        for artifact in array("artifacts")? {
            let id = independent_text(artifact, "id")?;
            if artifacts.insert(id, artifact).is_some() {
                return Err(format!("duplicate artifact {id}"));
            }
            let path = independent_text(artifact, "path")?;
            if let Some(previous) = artifact_paths.insert(path, id) {
                return Err(format!(
                    "artifact path {path} shared by {previous} and {id}"
                ));
            }
        }
        let mut reachable_artifacts = BTreeSet::new();
        let mut configuration_copies = BTreeSet::new();
        let mut scenario_ids = BTreeSet::new();
        let mut general_count = 0_usize;
        for input in array("configuration_inputs")? {
            let input_kind = independent_text(input, "kind")?;
            let copy = &input["copy"];
            let copy_id = independent_text(copy, "id")?;
            let catalog = artifacts
                .get(copy_id)
                .ok_or_else(|| format!("configuration copy absent from catalog: {copy_id}"))?;
            if **catalog != *copy {
                return Err(format!("configuration copy contradicts catalog: {copy_id}"));
            }
            if !configuration_copies.insert(copy_id) {
                return Err(format!("configuration copy reused: {copy_id}"));
            }
            reachable_artifacts.insert(copy_id);
            match input_kind {
                "general" => {
                    general_count += 1;
                    if copy_id != "config-general"
                        || independent_text(copy, "kind")? != "config_general_toml"
                        || independent_text(copy, "path")? != "inputs/config/general.toml"
                    {
                        return Err("non-canonical general configuration copy".to_owned());
                    }
                }
                "scenario" => {
                    let scenario_id = copy_id
                        .strip_prefix("config-scenario-")
                        .filter(|id| !id.is_empty())
                        .ok_or_else(|| "scenario copy has no scenario identity".to_owned())?;
                    if independent_text(copy, "kind")? != "config_scenario_toml"
                        || independent_text(copy, "path")?
                            != format!("inputs/config/scenarios/{scenario_id}.toml")
                        || !scenario_ids.insert(scenario_id)
                    {
                        return Err(format!("non-canonical scenario copy {copy_id}"));
                    }
                }
                other => return Err(format!("unknown configuration kind {other}")),
            }
        }
        if general_count != 1 {
            return Err(format!(
                "expected one general configuration, got {general_count}"
            ));
        }
        let resolution_link = &document["resolution_evidence"];
        let resolution_id = independent_text(resolution_link, "id")?;
        if resolution_id != "resolution-evidence"
            || independent_text(resolution_link, "kind")? != "resolution_evidence"
            || independent_text(resolution_link, "path")? != "inputs/resolution/evidence.json"
        {
            return Err("non-canonical resolution evidence link".to_owned());
        }
        let resolution_catalog = artifacts
            .get(resolution_id)
            .ok_or_else(|| "resolution evidence absent from catalog".to_owned())?;
        if **resolution_catalog != *resolution_link {
            return Err("resolution evidence contradicts catalog".to_owned());
        }
        reachable_artifacts.insert(resolution_id);
        let mut references = BTreeMap::<&str, &serde_json::Value>::new();
        let mut fasta_owners = BTreeMap::<&str, &str>::new();
        for reference in array("references")? {
            let id = independent_text(reference, "id")?;
            if references.insert(id, reference).is_some() {
                return Err(format!("duplicate reference {id}"));
            }
            let fasta_id = independent_text(&reference["fasta"], "id")?;
            if let Some(previous) = fasta_owners.insert(fasta_id, id) {
                return Err(format!("FASTA {fasta_id} reused by {previous} and {id}"));
            }
            let catalog = artifacts
                .get(fasta_id)
                .ok_or_else(|| format!("reference FASTA absent from catalog: {fasta_id}"))?;
            if **catalog != reference["fasta"] {
                return Err(format!("reference {id} contradicts artifact {fasta_id}"));
            }
            reachable_artifacts.insert(fasta_id);
            if independent_text(reference, "authoritative_resource")?
                != independent_text(&reference["fasta"], "path")?
            {
                return Err(format!(
                    "reference {id} authoritative resource contradicts its FASTA"
                ));
            }
        }
        if !document["effective_config"].is_null() {
            let effective_id = independent_text(&document["effective_config"], "id")?;
            if effective_id != "effective-config"
                || independent_text(&document["effective_config"], "kind")? != "effective_config"
            {
                return Err("effective configuration has a foreign identity".to_owned());
            }
            let catalog = artifacts
                .get(effective_id)
                .ok_or_else(|| "effective configuration is absent from catalog".to_owned())?;
            if **catalog != document["effective_config"] {
                return Err("effective configuration contradicts catalog".to_owned());
            }
            reachable_artifacts.insert(effective_id);
        }
        let mut event_ids = BTreeSet::new();
        let mut event_manifests = BTreeSet::new();
        for event in array("events")? {
            let event_id = independent_text(event, "event_id")?;
            if !event_ids.insert(event_id) {
                return Err(format!("duplicate event {event_id}"));
            }
            let manifest_id = independent_text(&event["manifest"], "id")?;
            if manifest_id != format!("manifest-event-{event_id}")
                || independent_text(&event["manifest"], "kind")? != "manifest_event"
            {
                return Err(format!("event {event_id} has a foreign manifest identity"));
            }
            if !event_manifests.insert(manifest_id) {
                return Err(format!("event manifest reused: {manifest_id}"));
            }
            let catalog = artifacts
                .get(manifest_id)
                .ok_or_else(|| format!("event manifest absent from catalog: {manifest_id}"))?;
            if **catalog != event["manifest"] {
                return Err(format!("event {event_id} contradicts catalog"));
            }
            reachable_artifacts.insert(manifest_id);
        }
        let mut reductions = BTreeMap::<&str, &serde_json::Value>::new();
        for reduction in array("reference_reductions")? {
            let id = independent_text(reduction, "derived_reference_id")?;
            if reductions.insert(id, reduction).is_some() {
                return Err(format!("duplicate reduction {id}"));
            }
        }
        let mut units = BTreeMap::<&str, &serde_json::Value>::new();
        let mut execution_units = Vec::new();
        let mut event_units = BTreeSet::new();
        for unit in array("units")? {
            let kind = independent_text(unit, "kind")?;
            let id = independent_text(unit, "id")?;
            match kind {
                "reference_reduction" => {
                    if units.insert(id, unit).is_some() {
                        return Err(format!("duplicate reduction unit {id}"));
                    }
                }
                "event" => {
                    if !event_units.insert(id) {
                        return Err(format!("duplicate event unit {id}"));
                    }
                }
                "execution" => execution_units.push(unit),
                other => return Err(format!("unknown unit kind {other}")),
            }
        }
        if execution_units.len() != 1
            || independent_text(execution_units[0], "id")?
                != independent_text(&document, "execution_id")?
            || independent_text(execution_units[0], "state")?
                != independent_text(&document, "state")?
        {
            return Err("execution unit does not identify the manifest execution".to_owned());
        }
        if event_units != event_ids {
            return Err(format!(
                "event units do not equal event links: {event_units:?}/{event_ids:?}"
            ));
        }

        let mut claimed_references = BTreeSet::new();
        for (derived_id, unit) in &units {
            let state = independent_text(unit, "state")?;
            let Some(reduction) = reductions.get(derived_id) else {
                if state == "failed" {
                    for expected_id in [
                        format!("reference-fasta-{derived_id}"),
                        format!("reference-map-{derived_id}"),
                    ] {
                        if let Some(artifact) = artifacts.get(expected_id.as_str()) {
                            if independent_text(artifact, "state")? != "failed"
                                || artifact["validated"] != serde_json::json!(false)
                            {
                                return Err(format!(
                                    "partial reduction artifact is not failed: {expected_id}"
                                ));
                            }
                            reachable_artifacts.insert(independent_text(artifact, "id")?);
                        }
                    }
                    continue;
                }
                return Err(format!("complete unit without reduction {derived_id}"));
            };
            let source_id = independent_text(reduction, "source_reference_id")?;
            let source = references
                .get(source_id)
                .ok_or_else(|| format!("missing source reference {source_id}"))?;
            if source["origin"]["kind"] == "reduced" {
                return Err(format!("source reference {source_id} is itself reduced"));
            }
            if source["fasta"] != reduction["source_fasta"] {
                return Err(format!("source FASTA identity crossed for {derived_id}"));
            }

            let derived_fasta_id = independent_text(&reduction["derived_fasta"], "id")?;
            let map_id = independent_text(&reduction["map"], "id")?;
            if derived_fasta_id != format!("reference-fasta-{derived_id}")
                || map_id != format!("reference-map-{derived_id}")
            {
                return Err(format!(
                    "derived artifact identity crossed for {derived_id}"
                ));
            }
            for (id, link) in [
                (derived_fasta_id, &reduction["derived_fasta"]),
                (map_id, &reduction["map"]),
            ] {
                let catalog = artifacts
                    .get(id)
                    .ok_or_else(|| format!("reduction artifact absent from catalog: {id}"))?;
                if **catalog != *link {
                    return Err(format!("reduction link contradicts catalog artifact {id}"));
                }
                reachable_artifacts.insert(id);
            }

            match state {
                "complete" => {
                    let derived = references
                        .get(derived_id)
                        .ok_or_else(|| format!("missing derived reference {derived_id}"))?;
                    if derived["fasta"] != reduction["derived_fasta"]
                        || derived["origin"]["kind"] != "reduced"
                        || derived["origin"]["source_reference_id"] != *source_id
                    {
                        return Err(format!(
                            "derived reference identity crossed for {derived_id}"
                        ));
                    }
                    claimed_references.insert(*derived_id);
                }
                "failed" if references.contains_key(derived_id) => {
                    return Err(format!(
                        "failed reduction owns complete reference {derived_id}"
                    ));
                }
                "failed" => {}
                other => return Err(format!("unsupported reduction state {other}")),
            }
        }
        for derived_id in reductions.keys() {
            if !units.contains_key(derived_id) {
                return Err(format!("reduction without unit {derived_id}"));
            }
        }
        for (reference_id, reference) in references {
            if reference["origin"]["kind"] == "reduced"
                && !claimed_references.contains(reference_id)
            {
                return Err(format!("orphan reduced reference {reference_id}"));
            }
        }
        let observed_minimum =
            independent_string_set(&document["outputs"]["minimum"], "outputs.minimum")?;
        let observed_explicit = independent_string_set(
            &document["outputs"]["explicitly_requested"],
            "outputs.explicitly_requested",
        )?;
        let mut expected_minimum = BTreeSet::from(["manifest-general"]);
        let mut expected_explicit = BTreeSet::new();
        for (id, artifact) in &artifacts {
            let kind = independent_text(artifact, "kind")?;
            if kind == "reference_reduction_map"
                || (kind == "reference_fasta" && id.starts_with("reference-fasta-"))
            {
                expected_explicit.insert(*id);
            } else {
                expected_minimum.insert(*id);
            }
            if !matches!(
                kind,
                "config_general_toml"
                    | "config_scenario_toml"
                    | "resolution_evidence"
                    | "effective_config"
                    | "reference_fasta"
                    | "reference_reduction_map"
                    | "manifest_event"
            ) {
                return Err(format!("unsupported artifact class {kind}/{id}"));
            }
            if !reachable_artifacts.contains(id) {
                return Err(format!("artifact is unreachable from graph roots: {id}"));
            }
        }
        if observed_minimum != expected_minimum || observed_explicit != expected_explicit {
            return Err(format!(
                "output categories do not partition artifacts: {observed_minimum:?}/{expected_minimum:?}; {observed_explicit:?}/{expected_explicit:?}"
            ));
        }
        Ok(())
    }

    fn independent_linked_document_graph(
        manifest_path: &Path,
        document: &serde_json::Value,
    ) -> Result<(), String> {
        if document["effective_config"].is_null() {
            return Ok(());
        }
        let root = manifest_path.parent().unwrap_or_else(|| Path::new("."));
        let resolution_path = independent_text(&document["resolution_evidence"], "path")?;
        let resolution_evidence: serde_json::Value = serde_json::from_slice(
            &fs::read(root.join(resolution_path))
                .map_err(|error| format!("independent resolution evidence read: {error}"))?,
        )
        .map_err(|error| format!("independent resolution evidence JSON: {error}"))?;
        let effective_path = independent_text(&document["effective_config"], "path")?;
        let effective: serde_json::Value = serde_json::from_slice(
            &fs::read(root.join(effective_path))
                .map_err(|error| format!("independent effective read: {error}"))?,
        )
        .map_err(|error| format!("independent effective JSON: {error}"))?;
        independent_normalization_graph(&effective)?;
        independent_authoritative_fasta_graph(root, document)?;
        independent_effective_exclusive_graph(root, document, &effective, &resolution_evidence)?;
        if document["route"] != effective["route"]
            || document["normalizations"] != effective["normalizations"]
            || document["reference_reductions"] != effective["reduced_references"]
        {
            return Err("independent effective/general cross-document mismatch".to_owned());
        }
        let reference_id = independent_text(&effective["reference"], "id")?;
        let general_reference = document["references"]
            .as_array()
            .and_then(|references| {
                references
                    .iter()
                    .find(|reference| reference["id"].as_str() == Some(reference_id))
            })
            .ok_or_else(|| "independent effective reference absent".to_owned())?;
        if general_reference != &effective["reference"] {
            return Err("independent effective reference mismatch".to_owned());
        }
        let root_hex = independent_text(&effective["randomness"], "root_seed")?;
        for reduction in document["reference_reductions"]
            .as_array()
            .ok_or_else(|| "independent reductions missing".to_owned())?
        {
            if independent_text(&reduction["randomness"], "root_seed")? != root_hex {
                return Err("independent reduction root mismatch".to_owned());
            }
            let components = [
                (
                    "source_reference_id",
                    independent_text(reduction, "source_reference_id")?,
                ),
                (
                    "derived_reference_id",
                    independent_text(reduction, "derived_reference_id")?,
                ),
                ("algorithm", independent_text(reduction, "algorithm")?),
            ];
            let expected = independent_subseed(root_hex, &components)?;
            if independent_text(&reduction["randomness"], "subseed")? != expected {
                return Err("independent reduction subseed mismatch".to_owned());
            }
        }

        let effective_event_ids = effective["event_ids"]
            .as_array()
            .ok_or_else(|| "independent effective event_ids".to_owned())?;
        let effective_subseeds = effective["randomness"]["event_subseeds"]
            .as_array()
            .ok_or_else(|| "independent effective event_subseeds".to_owned())?;
        if effective_event_ids.len() != effective_subseeds.len() {
            return Err("independent effective event cardinality".to_owned());
        }
        let mut subseeds = BTreeMap::new();
        let mut previous: Option<&str> = None;
        for (event_id, record) in effective_event_ids.iter().zip(effective_subseeds) {
            let event_id = event_id
                .as_str()
                .ok_or_else(|| "independent event ID text".to_owned())?;
            if previous.is_some_and(|value| value >= event_id)
                || independent_text(record, "event_id")? != event_id
            {
                return Err("independent event order/identity".to_owned());
            }
            previous = Some(event_id);
            let expected = independent_subseed(
                root_hex,
                &[("event_id", event_id), ("reference_id", reference_id)],
            )?;
            if independent_text(record, "subseed")? != expected {
                return Err("independent effective event subseed".to_owned());
            }
            subseeds.insert(event_id, expected);
        }
        let general_events = document["events"]
            .as_array()
            .ok_or_else(|| "independent general events".to_owned())?;
        if document["state"] == "complete" && general_events.len() != subseeds.len() {
            return Err("independent complete event cardinality".to_owned());
        }
        for event in general_events {
            let event_id = independent_text(event, "event_id")?;
            let expected = subseeds
                .get(event_id)
                .ok_or_else(|| "independent general event unknown".to_owned())?;
            if independent_text(event, "subseed")? != expected {
                return Err("independent general event subseed".to_owned());
            }
            let event_path = independent_text(&event["manifest"], "path")?;
            let linked: serde_json::Value = serde_json::from_slice(
                &fs::read(root.join(event_path))
                    .map_err(|error| format!("independent event read: {error}"))?,
            )
            .map_err(|error| format!("independent event JSON: {error}"))?;
            if independent_text(&linked, "event_id")? != event_id
                || independent_text(&linked, "subseed")? != expected
                || linked["units"][0]["id"].as_str() != Some(event_id)
            {
                return Err("independent linked event identity/subseed/unit".to_owned());
            }
            let expected_reference = serde_json::json!({
                "id": effective["reference"]["id"],
                "fasta": effective["reference"]["fasta"],
                "indexes": effective["reference"]["indexes"],
            });
            if linked["references"][0] != expected_reference {
                return Err("independent linked event reference".to_owned());
            }
        }
        Ok(())
    }

    fn independent_authoritative_fasta_graph(
        root: &Path,
        document: &serde_json::Value,
    ) -> Result<(), String> {
        for reference in document["references"]
            .as_array()
            .ok_or_else(|| "independent references missing".to_owned())?
        {
            let path = independent_text(reference, "authoritative_resource")?;
            let bytes = fs::read(root.join(path))
                .map_err(|error| format!("independent FASTA read: {error}"))?;
            let mut observed = Vec::<(Vec<u8>, String, u64, u64)>::new();
            let mut current: Option<(Vec<u8>, String, u64, u64)> = None;
            for source_line in bytes.split_inclusive(|byte| *byte == b'\n') {
                let mut line = source_line.strip_suffix(b"\n").unwrap_or(source_line);
                line = line.strip_suffix(b"\r").unwrap_or(line);
                if let Some(header) = line.strip_prefix(b">") {
                    if let Some(record) = current.take() {
                        observed.push(record);
                    }
                    let first = header
                        .iter()
                        .position(|byte| !byte.is_ascii_whitespace())
                        .ok_or_else(|| "independent FASTA empty header".to_owned())?;
                    let end = header[first..]
                        .iter()
                        .position(u8::is_ascii_whitespace)
                        .map_or(header.len(), |offset| first + offset);
                    let id = std::str::from_utf8(&header[first..end])
                        .map_err(|_| "independent FASTA non-UTF8 ID".to_owned())?
                        .to_owned();
                    current = Some((header.to_vec(), id, 0, 0));
                } else {
                    let record = current
                        .as_mut()
                        .ok_or_else(|| "independent FASTA bases before header".to_owned())?;
                    for byte in line {
                        let base = byte.to_ascii_uppercase();
                        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                            return Err("independent FASTA alphabet".to_owned());
                        }
                        record.2 = record
                            .2
                            .checked_add(1)
                            .ok_or_else(|| "independent FASTA length overflow".to_owned())?;
                        record.3 = record
                            .3
                            .checked_add(u64::from(base == b'N'))
                            .ok_or_else(|| "independent FASTA N overflow".to_owned())?;
                    }
                }
            }
            if let Some(record) = current {
                observed.push(record);
            }
            let declared = reference["sequences"]
                .as_array()
                .ok_or_else(|| "independent FASTA catalog missing".to_owned())?;
            if observed.len() != declared.len() {
                return Err("independent FASTA/catalog cardinality".to_owned());
            }
            for ((header, id, length, unknown), sequence) in observed.iter().zip(declared) {
                if independent_text(sequence, "id")? != id
                    || independent_reversible_bytes(&sequence["header"])? != *header
                    || sequence["length"].as_u64() != Some(*length)
                    || sequence["unknown_base_count"].as_u64() != Some(*unknown)
                {
                    return Err("independent FASTA/catalog metadata mismatch".to_owned());
                }
            }
        }
        Ok(())
    }

    #[derive(Debug)]
    struct IndependentRuleRecord {
        parameter: String,
        value: serde_json::Value,
        selectors: serde_json::Value,
        location: String,
        automatic_object: Option<String>,
    }

    #[derive(Debug, Default)]
    struct IndependentResolutionIds {
        reference: Option<String>,
        sequences: BTreeMap<String, u64>,
        individuals: BTreeSet<String>,
        populations: BTreeSet<String>,
        subpopulations: BTreeSet<String>,
        phases: BTreeSet<String>,
        groups: BTreeSet<String>,
        blocks: BTreeSet<String>,
        events: BTreeSet<String>,
        scenarios: BTreeSet<String>,
    }

    fn independent_effective_exclusive_graph(
        root: &Path,
        document: &serde_json::Value,
        effective: &serde_json::Value,
        resolution_evidence: &serde_json::Value,
    ) -> Result<(), String> {
        let inputs = document["configuration_inputs"]
            .as_array()
            .ok_or_else(|| "independent configuration inputs".to_owned())?;
        let general = inputs
            .iter()
            .find(|input| input["kind"] == "general")
            .ok_or_else(|| "independent general input".to_owned())?;
        let source_path = independent_text(general, "original_path")?;
        if independent_text(effective, "source_path")? != source_path {
            return Err("independent effective source path".to_owned());
        }
        let config_path = independent_text(&general["copy"], "path")?;
        let config = fs::read(root.join(config_path))
            .map_err(|error| format!("independent general copy read: {error}"))?;
        let rules = independent_rule_catalog(&config)?;
        let mut resolution_ids = independent_resolution_ids(&config)?;
        independent_add_scenario_resolution_ids(root, inputs, &mut resolution_ids)?;
        independent_declared_exclusive_collections(&config, effective)?;
        if resolution_evidence["resolution_evidence_kind"].as_str()
            != Some("gps_resolution_evidence")
            || resolution_evidence["schema_version"].as_str() != Some("1")
            || resolution_evidence["parameter_order_version"].as_str() != Some("1")
            || effective["resolution_evidence_id"].as_str() != Some("resolution-evidence")
            || effective["parameter_order_version"].as_str() != Some("1")
        {
            return Err("independent resolution evidence header/link".to_owned());
        }
        let resolution_domain = resolution_evidence
            .get("resolution_domain")
            .ok_or_else(|| "independent resolution domain missing".to_owned())?;
        if effective.get("resolution_domain") != Some(resolution_domain) {
            return Err("independent effective/domain evidence mismatch".to_owned());
        }
        independent_resolution_domain_shape(&resolution_ids, resolution_domain)?;

        for field in ["scenario_ids", "event_ids", "individual_ids"] {
            let values = effective[field]
                .as_array()
                .ok_or_else(|| format!("independent effective {field}"))?;
            let mut previous: Option<&str> = None;
            for value in values {
                let value = value
                    .as_str()
                    .ok_or_else(|| format!("independent {field} text"))?;
                if previous.is_some_and(|prior| prior >= value) {
                    return Err(format!("independent {field} order/duplicate"));
                }
                previous = Some(value);
            }
        }
        if effective["scenario_ids"]
            != serde_json::json!(resolution_ids.scenarios.iter().collect::<Vec<_>>())
            || effective["event_ids"]
                != serde_json::json!(resolution_ids.events.iter().collect::<Vec<_>>())
        {
            return Err("independent scenario/event catalog mismatch".to_owned());
        }

        let controlled = effective["controlled_frequency_sets"]
            .as_array()
            .ok_or_else(|| "independent controlled sets".to_owned())?;
        let mut previous_set: Option<&str> = None;
        for set in controlled {
            let id = independent_text(set, "id")?;
            if previous_set.is_some_and(|prior| prior >= id) {
                return Err("independent controlled set order/duplicate".to_owned());
            }
            previous_set = Some(id);
            let mut previous_count: Option<(&str, &str)> = None;
            for count in set["counts"]
                .as_array()
                .ok_or_else(|| "independent controlled counts".to_owned())?
            {
                let identity = (
                    independent_text(count, "id")?,
                    independent_text(count, "category")?,
                );
                if previous_count.is_some_and(|prior| prior >= identity) {
                    return Err("independent controlled count order/duplicate".to_owned());
                }
                previous_count = Some(identity);
            }
        }

        let evidence_requests = resolution_evidence["resolution_requests"]
            .as_array()
            .ok_or_else(|| "independent resolution requests".to_owned())?;
        let mut previous_request: Option<(&str, usize)> = None;
        for request in evidence_requests {
            let object_id = independent_text(request, "object_id")?;
            let parameter = independent_text(request, "parameter")?;
            let identity = (object_id, independent_parameter_rank(parameter)?);
            if previous_request.is_some_and(|prior| prior >= identity)
                || request.get("resolution_context").is_none()
            {
                return Err("independent resolution request order/duplicate/context".to_owned());
            }
            previous_request = Some(identity);
        }

        let values = effective["effective_values"]
            .as_array()
            .ok_or_else(|| "independent effective values".to_owned())?;
        if effective["parameter_order_version"].as_str() != Some("1") {
            return Err("independent parameter order version".to_owned());
        }
        let mut previous: Option<(&str, usize)> = None;
        for value in values {
            let object_id = independent_text(value, "object_id")?;
            let parameter = independent_text(value, "parameter")?;
            let context = value
                .get("resolution_context")
                .ok_or_else(|| "independent resolution context missing".to_owned())?;
            let identity = (object_id, independent_parameter_rank(parameter)?);
            if previous.is_some_and(|prior| prior >= identity) {
                return Err("independent effective value order/duplicate".to_owned());
            }
            previous = Some(identity);
            let winner = independent_text(value, "winning_rule_id")?;
            independent_validate_resolution_object(
                &resolution_ids,
                object_id,
                parameter,
                context,
                resolution_domain,
            )?;
            let matching = rules
                .iter()
                .filter(|(_, rule)| {
                    rule.parameter == parameter
                        && independent_rule_matches(&rule.selectors, context)
                })
                .collect::<Vec<_>>();
            let maximal = matching
                .iter()
                .copied()
                .filter(|(_, candidate)| {
                    !matching.iter().any(|(_, other)| {
                        independent_selector_dominates(&other.selectors, &candidate.selectors)
                    })
                })
                .collect::<Vec<_>>();
            let [(resolved_id, rule)] = maximal.as_slice() else {
                return Err("independent effective resolution is absent or ambiguous".to_owned());
            };
            if winner != resolved_id.as_str() {
                return Err("independent effective winner is not the dominant rule".to_owned());
            }
            if parameter != rule.parameter
                || value["value"] != rule.value
                || value["selectors"] != rule.selectors
                || value["provenance"]
                    != serde_json::json!({"path": source_path, "location": rule.location})
            {
                return Err("independent effective winner/value/selectors/provenance".to_owned());
            }
            if parameter == "ploidy" {
                let sequence = effective["reference"]["sequences"]
                    .as_array()
                    .and_then(|sequences| {
                        sequences
                            .iter()
                            .find(|sequence| sequence["id"].as_str() == Some(object_id))
                    })
                    .ok_or_else(|| "independent ploidy sequence".to_owned())?;
                if sequence["ploidy"] != value["value"] {
                    return Err("independent effective/reference ploidy".to_owned());
                }
            }
            let is_authoritative_request = evidence_requests.iter().any(|request| {
                request["object_id"].as_str() == Some(object_id)
                    && request["parameter"].as_str() == Some(parameter)
                    && request.get("resolution_context") == Some(context)
            });
            let parameter_has_authoritative_request = evidence_requests
                .iter()
                .any(|request| request["parameter"].as_str() == Some(parameter));
            let is_expected_automatic = rules.values().any(|rule| {
                rule.parameter == parameter
                    && rule.automatic_object.as_deref() == Some(object_id)
                    && if parameter == "ploidy" {
                        context["sequence"].as_str() == Some(object_id)
                    } else {
                        !parameter_has_authoritative_request
                            && independent_context_is_default(context)
                    }
            });
            if !is_authoritative_request && !is_expected_automatic {
                return Err(
                    "independent effective value lacks authoritative request/evidence".to_owned(),
                );
            }
        }
        for request in evidence_requests {
            if !values.iter().any(|value| {
                value["object_id"] == request["object_id"]
                    && value["parameter"] == request["parameter"]
                    && value["resolution_context"] == request["resolution_context"]
            }) {
                return Err("independent authoritative resolution request missing".to_owned());
            }
        }
        for (rule_id, rule) in &rules {
            if let Some(object_id) = &rule.automatic_object {
                let present = values.iter().any(|value| {
                    value["object_id"].as_str() == Some(object_id)
                        && value["parameter"].as_str() == Some(rule.parameter.as_str())
                        && value["winning_rule_id"].as_str() == Some(rule_id.as_str())
                        && value.get("resolution_context").is_some_and(|context| {
                            if rule.parameter == "ploidy" {
                                context["sequence"].as_str() == Some(object_id)
                            } else {
                                independent_context_is_default(context)
                            }
                        })
                });
                let replaced_by_explicit_request = rule.parameter != "ploidy"
                    && evidence_requests.iter().any(|request| {
                        request["parameter"].as_str() == Some(rule.parameter.as_str())
                    });
                if (replaced_by_explicit_request && present)
                    || (!replaced_by_explicit_request && !present)
                {
                    return Err("independent automatic effective value missing".to_owned());
                }
            }
        }
        Ok(())
    }

    fn independent_parameter_rank(parameter: &str) -> Result<usize, String> {
        match parameter {
            "ploidy" => Ok(0),
            "variant_rate" => Ok(1),
            "selfing_rate" => Ok(2),
            "recombination_rate" => Ok(3),
            "coverage" => Ok(4),
            "error_rate" => Ok(5),
            "contamination_rate" => Ok(6),
            "output_enabled" => Ok(7),
            _ => Err(format!("independent unknown parameter: {parameter}")),
        }
    }

    fn independent_context_is_default(context: &serde_json::Value) -> bool {
        context.as_object().is_some_and(|context| {
            context.len() == 12 && context.values().all(serde_json::Value::is_null)
        })
    }

    fn independent_rule_matches(
        selectors: &serde_json::Value,
        context: &serde_json::Value,
    ) -> bool {
        let text_matches =
            |field: &str| selectors[field].is_null() || selectors[field] == context[field];
        if !text_matches("individual")
            || !text_matches("group")
            || !text_matches("block")
            || !text_matches("technical_event")
        {
            return false;
        }
        if let Some(population) = selectors["population"].as_object() {
            let field = match population.get("kind").and_then(serde_json::Value::as_str) {
                Some("population") => "population",
                Some("subpopulation") => "subpopulation",
                _ => return false,
            };
            if population.get("id") != context.get(field) {
                return false;
            }
        }
        if let Some(temporal) = selectors["temporal"].as_object() {
            let matches = match temporal.get("kind").and_then(serde_json::Value::as_str) {
                Some("phase") => temporal.get("id") == context.get("phase"),
                Some("generation") => temporal.get("generation") == context.get("generation"),
                _ => false,
            };
            if !matches {
                return false;
            }
        }
        if let Some(spatial) = selectors["spatial"].as_object() {
            if spatial.get("reference") != context.get("reference")
                || spatial.get("sequence") != context.get("sequence")
            {
                return false;
            }
            let matches = match spatial.get("kind").and_then(serde_json::Value::as_str) {
                Some("sequence") => true,
                Some("region") => context["position"].as_u64().is_some_and(|position| {
                    spatial["start"]
                        .as_u64()
                        .is_some_and(|start| start <= position)
                        && spatial["end"].as_u64().is_some_and(|end| position <= end)
                }),
                Some("locus") => spatial.get("locus") == context.get("locus"),
                _ => false,
            };
            if !matches {
                return false;
            }
        }
        true
    }

    fn independent_selector_dominates(
        candidate: &serde_json::Value,
        other: &serde_json::Value,
    ) -> bool {
        let ranks = |selectors: &serde_json::Value| {
            [
                u8::from(!selectors["individual"].is_null()),
                match selectors["population"]["kind"].as_str() {
                    Some("population") => 1,
                    Some("subpopulation") => 2,
                    _ => 0,
                },
                match selectors["temporal"]["kind"].as_str() {
                    Some("phase") => 1,
                    Some("generation") => 2,
                    _ => 0,
                },
                match selectors["spatial"]["kind"].as_str() {
                    Some("sequence") => 1,
                    Some("region") => 2,
                    Some("locus") => 3,
                    _ => 0,
                },
                u8::from(!selectors["group"].is_null()),
                u8::from(!selectors["block"].is_null()),
                u8::from(!selectors["technical_event"].is_null()),
            ]
        };
        let candidate = ranks(candidate);
        let other = ranks(other);
        candidate
            .iter()
            .zip(other)
            .all(|(left, right)| *left >= right)
            && candidate
                .iter()
                .zip(other)
                .any(|(left, right)| *left > right)
    }

    fn independent_validate_resolution_object(
        ids: &IndependentResolutionIds,
        object_id: &str,
        parameter: &str,
        context: &serde_json::Value,
        domain: &serde_json::Value,
    ) -> Result<(), String> {
        let object = context
            .as_object()
            .ok_or_else(|| "independent context object".to_owned())?;
        if object.len() != 12 {
            return Err("independent context field set".to_owned());
        }
        if independent_context_is_default(context) {
            return if object_id == "global" && parameter != "ploidy" {
                Ok(())
            } else {
                Err("independent global object/class".to_owned())
            };
        }
        let allowed_axes = independent_allowed_resolution_axes(parameter)?;
        for (field, present) in [
            ("individual", !context["individual"].is_null()),
            (
                "population",
                !context["population"].is_null() || !context["subpopulation"].is_null(),
            ),
            (
                "temporal",
                !context["phase"].is_null() || !context["generation"].is_null(),
            ),
            (
                "spatial",
                !context["reference"].is_null()
                    || !context["sequence"].is_null()
                    || !context["position"].is_null()
                    || !context["locus"].is_null(),
            ),
            ("group", !context["group"].is_null()),
            ("block", !context["block"].is_null()),
            ("technical_event", !context["technical_event"].is_null()),
        ] {
            if present && !allowed_axes.contains(field) {
                return Err(format!("independent forbidden resolution axis: {field}"));
            }
        }
        let spatial_level = if !context["locus"].is_null() {
            Some("locus")
        } else if !context["position"].is_null() {
            Some("region")
        } else if !context["reference"].is_null() || !context["sequence"].is_null() {
            Some("sequence")
        } else {
            None
        };
        if let Some(level) = spatial_level
            && !independent_allowed_spatial_levels(parameter)?.contains(level)
        {
            return Err(format!(
                "independent forbidden resolution spatial level: {level}"
            ));
        }
        if !context["subpopulation"].is_null() && context["population"].is_null() {
            return Err("independent subpopulation without population".to_owned());
        }
        if !context["generation"].is_null() && context["phase"].is_null() {
            return Err("independent generation without phase".to_owned());
        }
        let has_position = !context["position"].is_null();
        let has_locus = !context["locus"].is_null();
        if (has_position || has_locus)
            && (context["reference"].is_null() || context["sequence"].is_null())
        {
            return Err("independent orphan spatial coordinate".to_owned());
        }
        if has_locus && !has_position {
            return Err("independent locus without position".to_owned());
        }
        if let (Some(population), Some(subpopulation)) = (
            context["population"].as_str(),
            context["subpopulation"].as_str(),
        ) && !domain["population_memberships"]
            .as_array()
            .is_some_and(|memberships| {
                memberships.iter().any(|membership| {
                    membership["population"].as_str() == Some(population)
                        && membership["subpopulation"].as_str() == Some(subpopulation)
                })
            })
        {
            return Err("independent unproven population hierarchy".to_owned());
        }
        if let (Some(phase), Some(generation)) =
            (context["phase"].as_str(), context["generation"].as_u64())
        {
            let proven = domain["phase_generation_intervals"]
                .as_array()
                .is_some_and(|intervals| {
                    intervals.iter().any(|interval| {
                        interval["phase"].as_str() == Some(phase)
                            && interval["start"]
                                .as_u64()
                                .is_some_and(|start| start <= generation)
                            && interval["end"]
                                .as_u64()
                                .is_some_and(|end| generation <= end)
                    })
                });
            if !proven {
                return Err("independent unproven temporal hierarchy".to_owned());
            }
        }
        if let (Some(locus), Some(reference), Some(sequence), Some(position)) = (
            context["locus"].as_str(),
            context["reference"].as_str(),
            context["sequence"].as_str(),
            context["position"].as_u64(),
        ) {
            let proven = domain["loci"].as_array().is_some_and(|loci| {
                loci.iter().any(|entry| {
                    entry["id"].as_str() == Some(locus)
                        && entry["reference"].as_str() == Some(reference)
                        && entry["sequence"].as_str() == Some(sequence)
                        && entry["position"].as_u64() == Some(position)
                })
            });
            if !proven {
                return Err("independent unproven locus hierarchy".to_owned());
            }
        }
        let mut candidates = BTreeSet::new();
        for (field, declared) in [
            ("individual", &ids.individuals),
            ("population", &ids.populations),
            ("subpopulation", &ids.subpopulations),
            ("phase", &ids.phases),
            ("group", &ids.groups),
            ("block", &ids.blocks),
            ("technical_event", &ids.events),
        ] {
            if let Some(id) = context[field].as_str() {
                if !declared.contains(id) {
                    return Err(format!("independent undeclared {field}"));
                }
                candidates.insert(id);
            }
        }
        match (context["reference"].as_str(), context["sequence"].as_str()) {
            (Some(reference), Some(sequence))
                if ids.reference.as_deref() == Some(reference)
                    && ids.sequences.contains_key(sequence) =>
            {
                if let Some(position) = context["position"].as_u64()
                    && (position == 0
                        || position
                            > *ids
                                .sequences
                                .get(sequence)
                                .ok_or_else(|| "independent sequence length".to_owned())?)
                {
                    return Err("independent position outside sequence".to_owned());
                }
                if let Some(locus) = context["locus"].as_str() {
                    candidates.insert(locus);
                } else {
                    candidates.insert(sequence);
                }
            }
            (None, None) => {}
            _ => return Err("independent undeclared reference/sequence".to_owned()),
        }
        if !candidates.contains(object_id) {
            return Err("independent foreign object/class".to_owned());
        }
        if parameter == "ploidy" && context["sequence"].as_str() != Some(object_id) {
            return Err("independent ploidy object/class".to_owned());
        }
        Ok(())
    }

    fn independent_allowed_resolution_axes(
        parameter: &str,
    ) -> Result<BTreeSet<&'static str>, String> {
        let axes: &[&str] = match parameter {
            "ploidy" => &["spatial"],
            "variant_rate" => &[
                "individual",
                "population",
                "temporal",
                "spatial",
                "group",
                "block",
            ],
            "selfing_rate" => &["population", "temporal"],
            "recombination_rate" => &["spatial", "group"],
            "coverage" | "error_rate" => &["technical_event", "individual", "spatial"],
            "contamination_rate" | "output_enabled" => &["technical_event", "individual"],
            _ => return Err("independent unknown parameter axes".to_owned()),
        };
        Ok(axes.iter().copied().collect())
    }

    fn independent_allowed_spatial_levels(
        parameter: &str,
    ) -> Result<BTreeSet<&'static str>, String> {
        let levels: &[&str] = match parameter {
            "ploidy" | "recombination_rate" => &["sequence"],
            "variant_rate" => &["sequence", "region", "locus"],
            "coverage" | "error_rate" => &["sequence", "region"],
            "selfing_rate" | "contamination_rate" | "output_enabled" => &[],
            _ => return Err("independent unknown parameter spatial levels".to_owned()),
        };
        Ok(levels.iter().copied().collect())
    }

    fn independent_resolution_ids(bytes: &[u8]) -> Result<IndependentResolutionIds, String> {
        use toml::de::DeTable;

        let text = std::str::from_utf8(bytes)
            .map_err(|_| "independent resolution TOML UTF-8".to_owned())?;
        let parsed = DeTable::parse(text)
            .map_err(|error| format!("independent resolution TOML parse: {error}"))?;
        let table = parsed.get_ref();
        let mut result = IndependentResolutionIds::default();
        let reference = table
            .get("reference")
            .and_then(|value| value.get_ref().as_table())
            .ok_or_else(|| "independent reference table".to_owned())?;
        result.reference = reference
            .get("id")
            .and_then(|value| value.get_ref().as_str())
            .map(str::to_owned);
        if let Some(sequences) = reference.get("sequences") {
            for sequence in sequences
                .get_ref()
                .as_array()
                .ok_or_else(|| "independent reference sequences".to_owned())?
            {
                let id = sequence
                    .get_ref()
                    .as_table()
                    .and_then(|sequence| sequence.get("id"))
                    .and_then(|value| value.get_ref().as_str())
                    .ok_or_else(|| "independent sequence id".to_owned())?;
                let length = sequence
                    .get_ref()
                    .as_table()
                    .and_then(|sequence| sequence.get("length"))
                    .and_then(|value| value.get_ref().as_integer())
                    .and_then(|value| value.as_str().replace('_', "").parse::<u64>().ok())
                    .ok_or_else(|| "independent sequence length".to_owned())?;
                result.sequences.insert(id.to_owned(), length);
            }
        }
        let collect = |field: &str| -> Result<BTreeSet<String>, String> {
            let mut ids = BTreeSet::new();
            if let Some(values) = table.get(field) {
                for value in values
                    .get_ref()
                    .as_array()
                    .ok_or_else(|| format!("independent {field} array"))?
                {
                    let id = value
                        .get_ref()
                        .as_table()
                        .and_then(|value| value.get("id"))
                        .and_then(|value| value.get_ref().as_str())
                        .ok_or_else(|| format!("independent {field} id"))?;
                    ids.insert(id.to_owned());
                }
            }
            Ok(ids)
        };
        result.individuals = collect("individuals")?;
        result.populations = collect("populations")?;
        result.subpopulations = collect("subpopulations")?;
        result.phases = collect("phases")?;
        result.groups = collect("groups")?;
        result.blocks = collect("blocks")?;
        result.events = collect("events")?;
        result.scenarios = collect("scenarios")?;
        Ok(result)
    }

    fn independent_add_scenario_resolution_ids(
        root: &Path,
        inputs: &[serde_json::Value],
        ids: &mut IndependentResolutionIds,
    ) -> Result<(), String> {
        use toml::de::DeTable;

        let mut observed = BTreeSet::new();
        for input in inputs
            .iter()
            .filter(|input| input["kind"].as_str() == Some("scenario"))
        {
            let copy = &input["copy"];
            let input_id = independent_text(copy, "id")?
                .strip_prefix("config-scenario-")
                .filter(|id| !id.is_empty())
                .ok_or_else(|| "independent scenario copy identity".to_owned())?;
            let path = independent_text(copy, "path")?;
            if path != format!("inputs/config/scenarios/{input_id}.toml")
                || !ids.scenarios.contains(input_id)
                || !observed.insert(input_id.to_owned())
            {
                return Err("independent scenario input identity/path mismatch".to_owned());
            }
            let bytes = fs::read(root.join(path))
                .map_err(|error| format!("independent scenario copy read: {error}"))?;
            let text = std::str::from_utf8(&bytes)
                .map_err(|_| "independent scenario TOML UTF-8".to_owned())?;
            let parsed = DeTable::parse(text)
                .map_err(|error| format!("independent scenario TOML parse: {error}"))?;
            let table = parsed.get_ref();
            if table
                .get("config_kind")
                .and_then(|value| value.get_ref().as_str())
                != Some("scenario")
                || table.get("id").and_then(|value| value.get_ref().as_str()) != Some(input_id)
            {
                return Err("independent scenario copy class/id mismatch".to_owned());
            }
            for (field, target) in [("phases", &mut ids.phases), ("events", &mut ids.events)] {
                if let Some(values) = table.get(field) {
                    for value in values
                        .get_ref()
                        .as_array()
                        .ok_or_else(|| format!("independent scenario {field} array"))?
                    {
                        let id = value
                            .get_ref()
                            .as_table()
                            .and_then(|entry| entry.get("id"))
                            .and_then(|id| id.get_ref().as_str())
                            .ok_or_else(|| format!("independent scenario {field} id"))?;
                        if !target.insert(id.to_owned()) {
                            return Err(format!("independent duplicate scenario {field} id"));
                        }
                    }
                }
            }
        }
        if observed != ids.scenarios {
            return Err("independent scenario inputs do not cover general links".to_owned());
        }
        Ok(())
    }

    fn independent_resolution_domain_shape(
        ids: &IndependentResolutionIds,
        domain: &serde_json::Value,
    ) -> Result<(), String> {
        if domain["class"].as_str() != Some("resolution_domain")
            || domain["version"].as_str() != Some("1")
            || domain["provider"].as_str().is_none_or(str::is_empty)
        {
            return Err("independent resolution domain header".to_owned());
        }
        for (field, key) in [
            ("population_memberships", "subpopulation"),
            ("phase_generation_intervals", "phase"),
            ("loci", "id"),
        ] {
            let values = domain[field]
                .as_array()
                .ok_or_else(|| format!("independent resolution domain {field}"))?;
            let mut previous: Option<&str> = None;
            for value in values {
                let identity = value[key]
                    .as_str()
                    .ok_or_else(|| format!("independent resolution domain {field} identity"))?;
                if previous.is_some_and(|prior| prior >= identity) {
                    return Err(format!(
                        "independent resolution domain {field} order/duplicate"
                    ));
                }
                previous = Some(identity);
                if field == "phase_generation_intervals"
                    && !matches!(
                        (value["start"].as_u64(), value["end"].as_u64()),
                        (Some(start), Some(end)) if start <= end && end <= u64::from(u32::MAX)
                    )
                {
                    return Err("independent resolution phase interval".to_owned());
                }
                match field {
                    "population_memberships"
                        if !value["population"]
                            .as_str()
                            .is_some_and(|id| ids.populations.contains(id))
                            || !value["subpopulation"]
                                .as_str()
                                .is_some_and(|id| ids.subpopulations.contains(id)) =>
                    {
                        return Err("independent resolution population identity".to_owned());
                    }
                    "phase_generation_intervals"
                        if !value["phase"]
                            .as_str()
                            .is_some_and(|id| ids.phases.contains(id)) =>
                    {
                        return Err("independent resolution phase identity".to_owned());
                    }
                    "loci" => {
                        let position = value["position"]
                            .as_u64()
                            .ok_or_else(|| "independent resolution locus position".to_owned())?;
                        let sequence = value["sequence"]
                            .as_str()
                            .ok_or_else(|| "independent resolution locus sequence".to_owned())?;
                        if position == 0
                            || ids.reference.as_deref() != value["reference"].as_str()
                            || !ids
                                .sequences
                                .get(sequence)
                                .is_some_and(|length| position <= *length)
                        {
                            return Err("independent resolution locus coordinate".to_owned());
                        }
                    }
                    _ => {}
                }
            }
        }
        Ok(())
    }

    fn independent_declared_exclusive_collections(
        bytes: &[u8],
        effective: &serde_json::Value,
    ) -> Result<(), String> {
        use toml::de::DeTable;

        let text = std::str::from_utf8(bytes)
            .map_err(|_| "independent exclusive TOML UTF-8".to_owned())?;
        let parsed = DeTable::parse(text)
            .map_err(|error| format!("independent exclusive TOML parse: {error}"))?;
        let table = parsed.get_ref();
        let mut individual_ids = table
            .get("individuals")
            .map(|value| {
                value
                    .get_ref()
                    .as_array()
                    .ok_or_else(|| "independent individuals array".to_owned())?
                    .iter()
                    .map(|entry| {
                        entry
                            .get_ref()
                            .as_table()
                            .and_then(|entry| entry.get("id"))
                            .and_then(|id| id.get_ref().as_str())
                            .map(str::to_owned)
                            .ok_or_else(|| "independent individual id".to_owned())
                    })
                    .collect::<Result<Vec<_>, String>>()
            })
            .transpose()?
            .unwrap_or_default();
        individual_ids.sort();
        if effective["individual_ids"] != serde_json::json!(individual_ids) {
            return Err("independent declared/effective individual IDs".to_owned());
        }

        let mut controlled = Vec::new();
        if let Some(sets) = table.get("frequency_sets") {
            for set in sets
                .get_ref()
                .as_array()
                .ok_or_else(|| "independent frequency sets array".to_owned())?
            {
                let set = set
                    .get_ref()
                    .as_table()
                    .ok_or_else(|| "independent frequency set table".to_owned())?;
                let mode = set
                    .get("mode")
                    .and_then(|value| value.get_ref().as_str())
                    .ok_or_else(|| "independent frequency mode".to_owned())?;
                if mode != "controlled" {
                    continue;
                }
                let id = set
                    .get("id")
                    .and_then(|value| value.get_ref().as_str())
                    .ok_or_else(|| "independent controlled id".to_owned())?;
                let mut counts = set
                    .get("counts")
                    .and_then(|value| value.get_ref().as_array())
                    .ok_or_else(|| "independent controlled counts array".to_owned())?
                    .iter()
                    .map(|count| {
                        let count = count
                            .get_ref()
                            .as_table()
                            .ok_or_else(|| "independent controlled count table".to_owned())?;
                        let id = count
                            .get("id")
                            .and_then(|value| value.get_ref().as_str())
                            .ok_or_else(|| "independent controlled count id".to_owned())?;
                        let category = count
                            .get("category")
                            .and_then(|value| value.get_ref().as_str())
                            .ok_or_else(|| "independent controlled count category".to_owned())?;
                        let value = count
                            .get("count")
                            .and_then(|value| value.get_ref().as_integer())
                            .and_then(|value| value.as_str().replace('_', "").parse::<u64>().ok())
                            .ok_or_else(|| "independent controlled count value".to_owned())?;
                        Ok(serde_json::json!({"id": id, "category": category, "count": value}))
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                counts.sort_by(|first, second| {
                    (first["id"].as_str(), first["category"].as_str())
                        .cmp(&(second["id"].as_str(), second["category"].as_str()))
                });
                controlled.push(serde_json::json!({"id": id, "counts": counts}));
            }
        }
        controlled.sort_by(|first, second| first["id"].as_str().cmp(&second["id"].as_str()));
        if effective["controlled_frequency_sets"] != serde_json::Value::Array(controlled) {
            return Err("independent declared/effective controlled sets".to_owned());
        }
        Ok(())
    }

    fn independent_rule_catalog(
        bytes: &[u8],
    ) -> Result<BTreeMap<String, IndependentRuleRecord>, String> {
        let text = std::str::from_utf8(bytes)
            .map_err(|_| "independent general TOML not UTF-8".to_owned())?;
        let mut catalog = BTreeMap::new();
        let mut fields = BTreeMap::<String, String>::new();
        let mut location: Option<String> = None;
        let mut offset = 0_usize;
        let finish = |fields: &mut BTreeMap<String, String>,
                      location: &mut Option<String>,
                      catalog: &mut BTreeMap<String, IndependentRuleRecord>|
         -> Result<(), String> {
            let Some(rule_location) = location.take() else {
                return Ok(());
            };
            let id = independent_toml_string(
                fields
                    .remove("id")
                    .ok_or_else(|| "independent rule id".to_owned())?
                    .as_str(),
            )?;
            let parameter = independent_toml_string(
                fields
                    .remove("parameter")
                    .ok_or_else(|| "independent rule parameter".to_owned())?
                    .as_str(),
            )?;
            let value = independent_toml_scalar(
                fields
                    .remove("value")
                    .ok_or_else(|| "independent rule value".to_owned())?
                    .as_str(),
            )?;
            let string = |key: &str| -> Result<Option<String>, String> {
                fields
                    .get(key)
                    .map(|value| independent_toml_string(value))
                    .transpose()
            };
            let integer = |key: &str| -> Result<Option<u64>, String> {
                fields
                    .get(key)
                    .map(|value| {
                        value
                            .replace('_', "")
                            .parse::<u64>()
                            .map_err(|_| format!("independent selector integer {key}"))
                    })
                    .transpose()
            };
            let reference = string("reference")?;
            let sequence = string("sequence")?;
            let spatial = match (
                reference.as_deref(),
                sequence.as_deref(),
                integer("start")?,
                integer("end")?,
                integer("locus")?,
            ) {
                (Some(reference), Some(sequence), _, _, Some(locus)) => Some(serde_json::json!({
                    "kind": "locus", "reference": reference, "sequence": sequence, "locus": locus
                })),
                (Some(reference), Some(sequence), Some(start), Some(end), None) => {
                    Some(serde_json::json!({
                        "kind": "region", "reference": reference, "sequence": sequence, "start": start, "end": end
                    }))
                }
                (Some(reference), Some(sequence), None, None, None) => Some(serde_json::json!({
                    "kind": "sequence", "reference": reference, "sequence": sequence
                })),
                (None, None, None, None, None) => None,
                _ => return Err("independent spatial selector shape".to_owned()),
            };
            let population = if let Some(id) = string("population")? {
                Some(serde_json::json!({"kind": "population", "id": id}))
            } else {
                string("subpopulation")?
                    .map(|id| serde_json::json!({"kind": "subpopulation", "id": id}))
            };
            let temporal = if let Some(id) = string("phase")? {
                Some(serde_json::json!({"kind": "phase", "id": id}))
            } else {
                integer("generation")?
                    .map(|generation| serde_json::json!({"kind": "generation", "generation": generation}))
            };
            let selectors = serde_json::json!({
                "individual": string("individual")?,
                "population": population,
                "temporal": temporal,
                "spatial": spatial,
                "group": string("group")?,
                "block": string("block")?,
                "technical_event": string("technical_event")?,
            });
            let automatic_object = if parameter != "ploidy"
                && selectors
                    .as_object()
                    .is_some_and(|selectors| selectors.values().all(serde_json::Value::is_null))
            {
                Some("global".to_owned())
            } else if parameter == "ploidy" {
                selectors["spatial"]
                    .as_object()
                    .filter(|spatial| {
                        spatial.get("kind").and_then(serde_json::Value::as_str) == Some("sequence")
                    })
                    .and_then(|spatial| spatial.get("sequence"))
                    .and_then(serde_json::Value::as_str)
                    .map(str::to_owned)
            } else {
                None
            };
            if catalog
                .insert(
                    id,
                    IndependentRuleRecord {
                        parameter,
                        value,
                        selectors,
                        location: rule_location,
                        automatic_object,
                    },
                )
                .is_some()
            {
                return Err("independent duplicate rule id".to_owned());
            }
            fields.clear();
            Ok(())
        };
        for (line_index, source_line) in text.split_inclusive('\n').enumerate() {
            let line = source_line.strip_suffix('\n').unwrap_or(source_line);
            let trimmed = line.trim();
            if trimmed == "[[rules]]" {
                finish(&mut fields, &mut location, &mut catalog)?;
                let column = line.find("[[rules]]").unwrap_or(0) + 1;
                let start = offset + column - 1;
                location = Some(format!(
                    "line {}, column {column}, bytes {start}..{}",
                    line_index + 1,
                    start + "[[rules]]".len()
                ));
            } else if location.is_some() && trimmed.starts_with('[') {
                finish(&mut fields, &mut location, &mut catalog)?;
            } else if location.is_some()
                && !trimmed.is_empty()
                && !trimmed.starts_with('#')
                && let Some((key, value)) = trimmed.split_once('=')
            {
                fields.insert(key.trim().to_owned(), value.trim().to_owned());
            }
            offset += source_line.len();
        }
        finish(&mut fields, &mut location, &mut catalog)?;
        Ok(catalog)
    }

    fn independent_toml_string(value: &str) -> Result<String, String> {
        value
            .strip_prefix('"')
            .and_then(|value| value.strip_suffix('"'))
            .map(str::to_owned)
            .ok_or_else(|| "independent TOML basic string".to_owned())
    }

    fn independent_toml_scalar(value: &str) -> Result<serde_json::Value, String> {
        if value == "true" || value == "false" {
            return Ok(serde_json::json!(value == "true"));
        }
        if value.starts_with('"') {
            return independent_toml_string(value).map(serde_json::Value::String);
        }
        let compact = value.replace('_', "");
        if compact.contains(['.', 'e', 'E']) {
            let number = compact
                .parse::<f64>()
                .map_err(|_| "independent TOML float".to_owned())?;
            serde_json::Number::from_f64(number)
                .map(serde_json::Value::Number)
                .ok_or_else(|| "independent non-finite TOML float".to_owned())
        } else {
            compact
                .parse::<i64>()
                .map(serde_json::Number::from)
                .map(serde_json::Value::Number)
                .map_err(|_| "independent TOML integer".to_owned())
        }
    }

    fn independent_subseed(root_hex: &str, components: &[(&str, &str)]) -> Result<String, String> {
        if root_hex.len() != 64 {
            return Err("independent root seed length".to_owned());
        }
        let mut root = [0_u8; 32];
        for (destination, pair) in root.iter_mut().zip(root_hex.as_bytes().chunks_exact(2)) {
            let decode = |byte| match byte {
                b'0'..=b'9' => Ok(byte - b'0'),
                b'a'..=b'f' => Ok(byte - b'a' + 10),
                _ => Err("independent root seed hex".to_owned()),
            };
            *destination = decode(pair[0])? * 16 + decode(pair[1])?;
        }
        let mut material = root.to_vec();
        material.extend_from_slice(&1_u32.to_be_bytes());
        material.extend_from_slice(
            &u32::try_from(components.len())
                .map_err(|_| "independent component count".to_owned())?
                .to_be_bytes(),
        );
        for (name, value) in components {
            for bytes in [name.as_bytes(), value.as_bytes()] {
                material.extend_from_slice(
                    &u32::try_from(bytes.len())
                        .map_err(|_| "independent component length".to_owned())?
                        .to_be_bytes(),
                );
                material.extend_from_slice(bytes);
            }
        }
        let digest = blake3::derive_key("genome_population_simulator/v1/subseed", &material);
        const HEX: &[u8; 16] = b"0123456789abcdef";
        let mut encoded = String::with_capacity(64);
        for byte in digest {
            encoded.push(char::from(HEX[usize::from(byte >> 4)]));
            encoded.push(char::from(HEX[usize::from(byte & 15)]));
        }
        Ok(encoded)
    }

    fn complete_and_failed_manifest_corpus(
        count: usize,
        label: &str,
    ) -> (PathBuf, serde_json::Value, serde_json::Value) {
        let (root, config) = reduction_fixture(count, label);
        let effective = load_and_validate_config(&config).expect("corpus configuration");
        let complete_destination = root.join("complete");
        execute_validated_config_with(
            &effective,
            &complete_destination,
            ExecutionFault::None,
            &LinuxAtomicPublisher,
        )
        .expect("complete corpus execution");
        let failed_destination = root.join("failed");
        execute_validated_config_with(
            &effective,
            &failed_destination,
            ExecutionFault::AtReduction {
                reduction_index: count - 1,
                stage: ReductionFaultStage::Metadata,
            },
            &LinuxAtomicPublisher,
        )
        .expect("failed corpus execution");
        let complete = serde_json::from_slice(
            &fs::read(complete_destination.join("manifest.json")).expect("complete manifest"),
        )
        .expect("complete manifest JSON");
        let failed = serde_json::from_slice(
            &fs::read(failed_destination.join("manifest.json")).expect("failed manifest"),
        )
        .expect("failed manifest JSON");
        (root, complete, failed)
    }

    fn percentage_manifest_corpus(literal: &str, label: &str) -> (PathBuf, serde_json::Value) {
        let (root, config_path) = reduction_fixture(1, label);
        let config = fs::read_to_string(&config_path)
            .expect("percentage configuration")
            .replace(
                "target_bases = 13",
                &format!("reduction_percent = {literal}"),
            );
        fs::write(&config_path, config).expect("percentage configuration write");
        let effective =
            load_and_validate_config(&config_path).expect("percentage corpus configuration");
        let destination = root.join("complete");
        execute_validated_config_with(
            &effective,
            &destination,
            ExecutionFault::None,
            &LinuxAtomicPublisher,
        )
        .expect("percentage corpus execution");
        let manifest = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("percentage manifest"),
        )
        .expect("percentage manifest JSON");
        (root, manifest)
    }

    fn non_utf8_header_manifest_corpus(label: &str) -> (PathBuf, serde_json::Value) {
        let root = temporary(label);
        if root.exists() {
            fs::remove_dir_all(&root).expect("stale non-UTF8 corpus cleanup");
        }
        fs::create_dir_all(&root).expect("non-UTF8 corpus root");
        fs::write(root.join("source.fa"), b">chr1 description \xff\nACGTN\n")
            .expect("non-UTF8 FASTA");
        fs::write(
            root.join("general.toml"),
            format!(
                r#"config_kind = "general"
schema_version = "1"
route = "direct"
seed = "{}"

[reference]
kind = "supplied"
id = "reference-main"
fasta = "source.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 5
ploidy = 2
"#,
                "00".repeat(32)
            ),
        )
        .expect("non-UTF8 config");
        let effective =
            load_and_validate_config(&root.join("general.toml")).expect("non-UTF8 effective");
        let destination = root.join("complete");
        execute_validated_config_with(
            &effective,
            &destination,
            ExecutionFault::None,
            &LinuxAtomicPublisher,
        )
        .expect("non-UTF8 corpus execution");
        let manifest = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("non-UTF8 manifest"),
        )
        .expect("non-UTF8 manifest JSON");
        (root, manifest)
    }

    fn supplied_two_sequence_manifest_corpus(label: &str) -> (PathBuf, serde_json::Value) {
        let root = temporary(label);
        if root.exists() {
            fs::remove_dir_all(&root).expect("stale supplied corpus cleanup");
        }
        fs::create_dir_all(&root).expect("supplied corpus root");
        fs::write(
            root.join("source.fa"),
            b">chr1 first\nACGTN\n>chr2 second \xff\nNNAC\n",
        )
        .expect("supplied FASTA");
        fs::write(
            root.join("general.toml"),
            r#"config_kind = "general"
schema_version = "1"
route = "direct"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "source.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 5
ploidy = 2

[[reference.sequences]]
id = "chr2"
class = "scaffold"
length = 4
ploidy = 2
"#,
        )
        .expect("supplied configuration");
        let effective = load_and_validate_config(&root.join("general.toml"))
            .expect("supplied configuration validation");
        let destination = root.join("complete");
        execute_validated_config_with(
            &effective,
            &destination,
            ExecutionFault::None,
            &LinuxAtomicPublisher,
        )
        .expect("supplied corpus execution");
        let manifest = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("supplied manifest"),
        )
        .expect("supplied manifest JSON");
        (root, manifest)
    }

    fn scoped_effective_manifest_corpus(
        label: &str,
    ) -> (PathBuf, serde_json::Value, serde_json::Value) {
        let root = temporary(label);
        if root.exists() {
            fs::remove_dir_all(&root).expect("stale effective corpus cleanup");
        }
        fs::create_dir_all(&root).expect("effective corpus root");
        let config_path = root.join("general.toml");
        let mut config = fs::read_to_string("tests/fixtures/config/general-effective.toml")
            .expect("effective fixture");
        config.push_str(
            r#"
[[events]]
id = "event-1"

[[rules]]
id = "output-base"
parameter = "output_enabled"
value = true

[[rules]]
id = "coverage-base"
parameter = "coverage"
value = 30.0

[[rules]]
id = "recombination-base"
parameter = "recombination_rate"
value = 0.01
"#,
        );
        fs::write(&config_path, config).expect("effective scoped configuration");
        let admitted = admit_toml_header(
            &config_path,
            fs::read(&config_path).expect("effective fixture"),
            ExpectedConfigKind::General,
        )
        .expect("effective header");
        let DeclaredConfig::General(general) =
            parse_declared_config(admitted).expect("effective declared config")
        else {
            panic!("expected general effective fixture");
        };
        let requests = [
            ResolutionRequest {
                object_id: "ind-1".to_owned(),
                parameter: ParameterKey::VariantRate,
                context: ResolutionContext {
                    individual: Some("ind-1".to_owned()),
                    ..ResolutionContext::default()
                },
            },
            ResolutionRequest {
                object_id: "event-1".to_owned(),
                parameter: ParameterKey::Coverage,
                context: ResolutionContext {
                    technical_event: Some("event-1".to_owned()),
                    ..ResolutionContext::default()
                },
            },
            ResolutionRequest {
                object_id: "chr1".to_owned(),
                parameter: ParameterKey::RecombinationRate,
                context: ResolutionContext {
                    reference: Some("ref-main".to_owned()),
                    sequence: Some("chr1".to_owned()),
                    ..ResolutionContext::default()
                },
            },
        ];
        let effective = build_effective_config(&general, &[], &requests)
            .expect("effective multiple-scope build");
        let destination = root.join("complete");
        execute_validated_config_with(
            &effective,
            &destination,
            ExecutionFault::None,
            &LinuxAtomicPublisher,
        )
        .expect("effective multiple-scope execution");
        let manifest = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("effective manifest"),
        )
        .expect("effective manifest JSON");
        let effective_document = serde_json::from_slice(
            &fs::read(destination.join("effective-config.json")).expect("effective document"),
        )
        .expect("effective document JSON");
        (root, manifest, effective_document)
    }

    fn zero_reduction_manifest_corpus(
        label: &str,
    ) -> (PathBuf, serde_json::Value, serde_json::Value) {
        let root = temporary(label);
        if root.exists() {
            fs::remove_dir_all(&root).expect("stale zero-reduction corpus cleanup");
        }
        fs::create_dir_all(&root).expect("zero-reduction corpus root");
        let effective = load_and_validate_config(
            &Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/config/general-seeded.toml"),
        )
        .expect("seeded zero-reduction configuration");
        let destination = root.join("complete");
        execute_validated_config_with(
            &effective,
            &destination,
            ExecutionFault::None,
            &LinuxAtomicPublisher,
        )
        .expect("zero-reduction corpus execution");
        let complete = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("zero-reduction manifest"),
        )
        .expect("zero-reduction manifest JSON");
        let failed_destination = root.join("failed");
        execute_validated_config_with(
            &effective,
            &failed_destination,
            ExecutionFault::AfterEffectiveConfig,
            &LinuxAtomicPublisher,
        )
        .expect("zero-reduction failed corpus execution");
        let failed = serde_json::from_slice(
            &fs::read(failed_destination.join("manifest.json"))
                .expect("zero-reduction failed manifest"),
        )
        .expect("zero-reduction failed manifest JSON");
        (root, complete, failed)
    }

    fn parameter_order_manifest_corpus(
        label: &str,
    ) -> (PathBuf, serde_json::Value, serde_json::Value) {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("parameter order corpus root");
        let config_path = root.join("general.toml");
        let mut config = fs::read_to_string("tests/fixtures/config/general-effective.toml")
            .expect("effective fixture");
        config.push_str(
            r#"

[[rules]]
id = "output-base"
parameter = "output_enabled"
value = false

[[rules]]
id = "output-ind-1"
parameter = "output_enabled"
value = true
individual = "ind-1"
"#,
        );
        fs::write(&config_path, config).expect("parameter order configuration");
        let admitted = admit_toml_header(
            &config_path,
            fs::read(&config_path).expect("parameter order fixture"),
            ExpectedConfigKind::General,
        )
        .expect("parameter order header");
        let DeclaredConfig::General(general) =
            parse_declared_config(admitted).expect("parameter order config")
        else {
            panic!("expected parameter order general");
        };
        let requests = [
            ResolutionRequest {
                object_id: "ind-1".to_owned(),
                parameter: ParameterKey::OutputEnabled,
                context: ResolutionContext {
                    individual: Some("ind-1".to_owned()),
                    ..ResolutionContext::default()
                },
            },
            ResolutionRequest {
                object_id: "ind-1".to_owned(),
                parameter: ParameterKey::VariantRate,
                context: ResolutionContext {
                    individual: Some("ind-1".to_owned()),
                    ..ResolutionContext::default()
                },
            },
        ];
        let effective = build_effective_config(&general, &[], &requests)
            .expect("same-object parameters resolve");
        for (state, fault) in [
            ("complete", ExecutionFault::None),
            ("failed", ExecutionFault::AfterEffectiveConfig),
        ] {
            execute_validated_config_with(
                &effective,
                &root.join(state),
                fault,
                &LinuxAtomicPublisher,
            )
            .expect("parameter order execution");
        }
        let read = |state: &str| {
            serde_json::from_slice(
                &fs::read(root.join(state).join("manifest.json"))
                    .expect("parameter order manifest"),
            )
            .expect("parameter order manifest JSON")
        };
        (root.clone(), read("complete"), read("failed"))
    }

    fn multiple_configuration_and_event_manifest_corpus(
        label: &str,
    ) -> (PathBuf, serde_json::Value, serde_json::Value) {
        let root = temporary(label);
        if root.exists() {
            fs::remove_dir_all(&root).expect("stale multiple-configuration corpus cleanup");
        }
        fs::create_dir_all(&root).expect("multiple-configuration corpus root");
        let fixture_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/config");
        let mut general =
            fs::read_to_string(fixture_root.join("general-evolution.toml")).expect("general");
        general.push_str(
            r#"

[[scenarios]]
id = "evolution-other"
path = "scenario-evolution-other.toml"
route = "evolution"

[[frequency_sets]]
id = "frequency-alpha"
mode = "probabilistic"
values = [
  { id = "A", category = "allele", value = 0.5 },
  { id = "B", category = "allele", value = 0.5 },
]

[[frequency_sets]]
id = "frequency-beta"
mode = "probabilistic"
values = [
  { id = "A", category = "allele", value = 0.25 },
  { id = "B", category = "allele", value = 0.75 },
]

[[rules]]
id = "selfing-generation-5"
parameter = "selfing_rate"
value = 0.5
generation = 5
"#,
        );
        fs::write(root.join("general.toml"), general).expect("write multiple scenario general");
        fs::copy(
            fixture_root.join("scenario-evolution.toml"),
            root.join("scenario-evolution.toml"),
        )
        .expect("copy first scenario");
        let second = fs::read_to_string(fixture_root.join("scenario-evolution.toml"))
            .expect("second scenario template")
            .replace("evolution-small", "evolution-other")
            .replace("phase-1", "phase-2")
            .replace("sample-generation-10", "sample-generation-20");
        fs::write(root.join("scenario-evolution-other.toml"), second)
            .expect("write second scenario");

        let general_path = root.join("general.toml");
        let admitted = admit_toml_header(
            &general_path,
            fs::read(&general_path).expect("multiple general bytes"),
            ExpectedConfigKind::General,
        )
        .expect("multiple general header");
        let DeclaredConfig::General(general) =
            parse_declared_config(admitted).expect("multiple general declaration")
        else {
            panic!("expected multiple general configuration");
        };
        let mut scenarios = Vec::new();
        for (id, path) in [
            ("evolution-small", "scenario-evolution.toml"),
            ("evolution-other", "scenario-evolution-other.toml"),
        ] {
            let scenario_path = root.join(path);
            let admitted = admit_toml_header(
                &scenario_path,
                fs::read(&scenario_path).expect("multiple scenario bytes"),
                ExpectedConfigKind::Scenario,
            )
            .expect("multiple scenario header");
            let DeclaredConfig::Scenario(scenario) =
                parse_declared_config(admitted).expect("multiple scenario declaration")
            else {
                panic!("expected scenario {id}");
            };
            assert_eq!(scenario.id().value(), id);
            scenarios.push(scenario);
        }
        let scenario_refs = scenarios.iter().map(Box::as_ref).collect::<Vec<_>>();
        let domain = ResolutionDomain::validated(
            "synthetic-scenario-phase-test",
            [],
            [PhaseGenerationInterval {
                phase: "phase-1".to_owned(),
                start: 1,
                end: 10,
            }],
            [],
        )
        .expect("synthetic scenario phase domain");
        let effective = build_effective_config_with_domain(
            &general,
            &scenario_refs,
            &[ResolutionRequest {
                object_id: "phase-1".to_owned(),
                parameter: ParameterKey::SelfingRate,
                context: ResolutionContext {
                    phase: Some("phase-1".to_owned()),
                    generation: Some(5),
                    ..ResolutionContext::default()
                },
            }],
            &domain,
        )
        .expect("scenario-only phase resolution");
        let complete_destination = root.join("complete");
        execute_validated_config_with(
            &effective,
            &complete_destination,
            ExecutionFault::None,
            &LinuxAtomicPublisher,
        )
        .expect("multiple-configuration complete execution");
        let failed_destination = root.join("failed");
        execute_validated_config_with(
            &effective,
            &failed_destination,
            ExecutionFault::AfterEffectiveConfig,
            &LinuxAtomicPublisher,
        )
        .expect("multiple-configuration failed execution");
        let read_manifest = |destination: &Path| {
            serde_json::from_slice(
                &fs::read(destination.join("manifest.json")).expect("graph manifest"),
            )
            .expect("graph manifest JSON")
        };
        (
            root,
            read_manifest(&complete_destination),
            read_manifest(&failed_destination),
        )
    }

    fn unit_index(document: &serde_json::Value, kind: &str) -> usize {
        document["units"]
            .as_array()
            .expect("units")
            .iter()
            .position(|unit| unit["kind"] == kind)
            .unwrap_or_else(|| panic!("missing {kind} unit"))
    }

    fn artifact_index(document: &serde_json::Value, id: &str) -> usize {
        document["artifacts"]
            .as_array()
            .expect("artifacts")
            .iter()
            .position(|artifact| artifact["id"] == id)
            .unwrap_or_else(|| panic!("missing artifact {id}"))
    }

    #[test]
    fn manifest_graph_rejects_every_dynamic_edge_contradiction() {
        let (root, complete, failed) =
            multiple_configuration_and_event_manifest_corpus("graph-regression");
        let (zero_root, zero_complete, zero_failed) = zero_reduction_manifest_corpus("graph-zero");
        let (one_root, one_complete, one_failed) =
            complete_and_failed_manifest_corpus(1, "graph-one");
        let (two_root, two_complete, two_failed) =
            complete_and_failed_manifest_corpus(2, "graph-two");
        let mut executed_cases = BTreeSet::new();
        for (case_root, name, document) in [
            (&zero_root, "graph-positive-zero-complete", &zero_complete),
            (&zero_root, "graph-positive-zero-failed", &zero_failed),
            (&one_root, "graph-positive-one-complete", &one_complete),
            (&one_root, "graph-positive-one-failed", &one_failed),
            (&two_root, "graph-positive-two-complete", &two_complete),
            (&two_root, "graph-positive-two-failed", &two_failed),
            (&root, "graph-positive-multiple-complete", &complete),
            (&root, "graph-positive-multiple-failed", &failed),
        ] {
            assert_global_manifest_positive(case_root, name, document);
            assert!(executed_cases.insert(name.to_owned()));
        }

        let scenario_copy = complete["configuration_inputs"][1]["copy"].clone();
        let general_copy = complete["configuration_inputs"][0]["copy"].clone();
        let event_unit = complete["units"][unit_index(&complete, "event")].clone();

        let mut mutations = Vec::new();
        for (name, pointer, value) in [
            (
                "graph-configuration-copy-digest",
                "/configuration_inputs/1/copy/checksum/digest",
                serde_json::json!("0".repeat(64)),
            ),
            (
                "graph-configuration-copy-id",
                "/configuration_inputs/1/copy/id",
                serde_json::json!("config-scenario-other"),
            ),
            (
                "graph-configuration-copy-class",
                "/configuration_inputs/1/copy/kind",
                serde_json::json!("config_general_toml"),
            ),
            (
                "graph-configuration-copy-path",
                "/configuration_inputs/1/copy/path",
                serde_json::json!("inputs/config/scenarios/other.toml"),
            ),
            (
                "graph-configuration-copy-size",
                "/configuration_inputs/1/copy/size_bytes",
                serde_json::json!(999_999_999_u64),
            ),
        ] {
            let mut document = complete.clone();
            *document
                .pointer_mut(pointer)
                .expect("graph mutation pointer") = value;
            mutations.push((&root, name, document));
        }

        let mut reused_general = complete.clone();
        reused_general["configuration_inputs"][1]["copy"] = general_copy.clone();
        mutations.push((&root, "graph-configuration-copy-reused", reused_general));
        let mut duplicate_configuration = complete.clone();
        duplicate_configuration["configuration_inputs"]
            .as_array_mut()
            .expect("configuration inputs")
            .push(complete["configuration_inputs"][1].clone());
        mutations.push((
            &root,
            "graph-configuration-input-duplicate",
            duplicate_configuration,
        ));
        let mut incompatible_kind = complete.clone();
        incompatible_kind["configuration_inputs"][0]["kind"] = serde_json::json!("scenario");
        mutations.push((
            &root,
            "graph-configuration-kind-incompatible",
            incompatible_kind,
        ));
        let mut missing_general = complete.clone();
        missing_general["configuration_inputs"]
            .as_array_mut()
            .expect("configuration inputs")
            .remove(0);
        mutations.push((
            &root,
            "graph-general-configuration-missing",
            missing_general,
        ));

        for (name, pointer, value) in [
            (
                "graph-resolution-evidence-digest",
                "/resolution_evidence/checksum/digest",
                serde_json::json!("0".repeat(64)),
            ),
            (
                "graph-resolution-evidence-id",
                "/resolution_evidence/id",
                serde_json::json!("resolution-evidence-other"),
            ),
            (
                "graph-resolution-evidence-kind",
                "/resolution_evidence/kind",
                serde_json::json!("effective_config_json"),
            ),
            (
                "graph-resolution-evidence-path",
                "/resolution_evidence/path",
                serde_json::json!("inputs/resolution/other.json"),
            ),
            (
                "graph-resolution-evidence-size",
                "/resolution_evidence/size_bytes",
                serde_json::json!(999_999_999_u64),
            ),
        ] {
            let mut document = complete.clone();
            *document
                .pointer_mut(pointer)
                .expect("resolution evidence mutation pointer") = value;
            mutations.push((&root, name, document));
        }

        let mut crossed_event_unit = complete.clone();
        let event_index = unit_index(&crossed_event_unit, "event");
        crossed_event_unit["units"][event_index]["id"] = serde_json::json!("event-other");
        mutations.push((&root, "graph-event-unit-crossed", crossed_event_unit));
        let mut orphan_event_unit = complete.clone();
        let mut orphan = event_unit;
        orphan["id"] = serde_json::json!("event-orphan");
        orphan_event_unit["units"]
            .as_array_mut()
            .expect("units")
            .push(orphan);
        mutations.push((&root, "graph-event-unit-orphan", orphan_event_unit));
        let mut missing_event_unit = complete.clone();
        let event_index = unit_index(&missing_event_unit, "event");
        missing_event_unit["units"]
            .as_array_mut()
            .expect("units")
            .remove(event_index);
        mutations.push((&root, "graph-event-unit-missing", missing_event_unit));
        let mut crossed_execution_unit = complete.clone();
        let execution_index = unit_index(&crossed_execution_unit, "execution");
        let execution_id = crossed_execution_unit["execution_id"]
            .as_str()
            .expect("execution ID");
        let replacement = if execution_id.ends_with('A') {
            'B'
        } else {
            'A'
        };
        let mut other_execution_id = execution_id[..execution_id.len() - 1].to_owned();
        other_execution_id.push(replacement);
        crossed_execution_unit["units"][execution_index]["id"] =
            serde_json::json!(other_execution_id);
        mutations.push((
            &root,
            "graph-execution-unit-crossed",
            crossed_execution_unit,
        ));

        let mut missing_output = complete.clone();
        missing_output["outputs"]["minimum"]
            .as_array_mut()
            .expect("minimum outputs")
            .push(serde_json::json!("artifact-missing"));
        mutations.push((&root, "graph-output-artifact-missing", missing_output));
        let mut wrong_category = complete.clone();
        let minimum = wrong_category["outputs"]["minimum"]
            .as_array_mut()
            .expect("minimum outputs");
        let config_position = minimum
            .iter()
            .position(|id| id == "config-general")
            .expect("general config output");
        let config_id = minimum.remove(config_position);
        wrong_category["outputs"]["explicitly_requested"]
            .as_array_mut()
            .expect("explicit outputs")
            .push(config_id);
        mutations.push((&root, "graph-output-wrong-category", wrong_category));
        let mut duplicate_path = complete.clone();
        let scenario_id = scenario_copy["id"].as_str().expect("scenario copy ID");
        let scenario_artifact_index = artifact_index(&duplicate_path, scenario_id);
        duplicate_path["artifacts"][scenario_artifact_index]["path"] = general_copy["path"].clone();
        duplicate_path["configuration_inputs"][1]["copy"] =
            duplicate_path["artifacts"][scenario_artifact_index].clone();
        mutations.push((&root, "graph-artifact-path-duplicate", duplicate_path));
        let mut orphan_event_artifact = complete.clone();
        let event_artifact_index = orphan_event_artifact["artifacts"]
            .as_array()
            .expect("artifacts")
            .iter()
            .position(|artifact| artifact["kind"] == "manifest_event")
            .expect("event artifact");
        let mut orphan_artifact = orphan_event_artifact["artifacts"][event_artifact_index].clone();
        orphan_artifact["id"] = serde_json::json!("manifest-event-orphan");
        orphan_artifact["path"] = serde_json::json!("events/orphan/manifest.json");
        orphan_event_artifact["artifacts"]
            .as_array_mut()
            .expect("artifacts")
            .push(orphan_artifact);
        orphan_event_artifact["outputs"]["minimum"]
            .as_array_mut()
            .expect("minimum outputs")
            .push(serde_json::json!("manifest-event-orphan"));
        mutations.push((&root, "graph-orphan-event-artifact", orphan_event_artifact));

        let mut failed_configuration_digest = failed.clone();
        failed_configuration_digest["configuration_inputs"][1]["copy"]["checksum"]["digest"] =
            serde_json::json!("0".repeat(64));
        mutations.push((
            &root,
            "graph-failed-configuration-copy-digest",
            failed_configuration_digest,
        ));
        let mut failed_execution = failed.clone();
        let failed_execution_index = unit_index(&failed_execution, "execution");
        let failed_id = failed_execution["execution_id"]
            .as_str()
            .expect("failed ID");
        let mut other_failed_id = failed_id[..failed_id.len() - 1].to_owned();
        other_failed_id.push(if failed_id.ends_with('A') { 'B' } else { 'A' });
        failed_execution["units"][failed_execution_index]["id"] =
            serde_json::json!(other_failed_id);
        mutations.push((
            &root,
            "graph-failed-execution-unit-crossed",
            failed_execution,
        ));

        for (case_root, name, source) in [
            (
                &one_root,
                "graph-one-complete-output-missing",
                &one_complete,
            ),
            (&one_root, "graph-one-failed-output-missing", &one_failed),
            (
                &two_root,
                "graph-two-complete-output-missing",
                &two_complete,
            ),
            (&two_root, "graph-two-failed-output-missing", &two_failed),
        ] {
            let mut document = source.clone();
            document["outputs"]["minimum"]
                .as_array_mut()
                .expect("minimum outputs")
                .push(serde_json::json!("artifact-missing"));
            mutations.push((case_root, name, document));
        }
        let mut two_duplicate_path = two_complete.clone();
        let config_artifact_index = artifact_index(&two_duplicate_path, "config-general");
        let reference_path = two_duplicate_path["references"][0]["fasta"]["path"].clone();
        two_duplicate_path["artifacts"][config_artifact_index]["path"] = reference_path;
        two_duplicate_path["configuration_inputs"][0]["copy"] =
            two_duplicate_path["artifacts"][config_artifact_index].clone();
        mutations.push((
            &two_root,
            "graph-two-complete-artifact-path-duplicate",
            two_duplicate_path,
        ));

        for (case_root, name, document) in mutations {
            assert_global_manifest_mutation(case_root, name, &document);
            assert!(executed_cases.insert(name.to_owned()));
        }
        assert_eq!(
            executed_cases,
            MANIFEST_GRAPH_CORPUS_CASES
                .iter()
                .map(|case| (*case).to_owned())
                .collect(),
            "the graph corpus catalog must equal the cases actually executed"
        );
        for case_root in [root, zero_root, one_root, two_root] {
            fs::remove_dir_all(case_root).expect("graph corpus cleanup");
        }
    }

    #[test]
    fn manifest_dynamic_graph_rejects_diagnostic_and_reference_inconsistencies() {
        let (zero_root, zero_complete, zero_failed) =
            zero_reduction_manifest_corpus("dynamic-zero");
        let (one_root, one_complete, one_failed) =
            complete_and_failed_manifest_corpus(1, "dynamic-one");
        let (two_root, two_complete, two_failed) =
            complete_and_failed_manifest_corpus(2, "dynamic-two");
        let (percent_root, percent_complete) = percentage_manifest_corpus("35", "dynamic-percent");

        let mut executed_cases = BTreeSet::new();
        for (root, name, document) in [
            (&zero_root, "dynamic-positive-zero-complete", &zero_complete),
            (&zero_root, "dynamic-positive-zero-failed", &zero_failed),
            (&one_root, "dynamic-positive-one-complete", &one_complete),
            (&one_root, "dynamic-positive-one-failed", &one_failed),
            (&two_root, "dynamic-positive-two-complete", &two_complete),
            (&two_root, "dynamic-positive-two-failed", &two_failed),
            (
                &percent_root,
                "dynamic-positive-percent-complete",
                &percent_complete,
            ),
        ] {
            assert_global_manifest_positive(root, name, document);
            assert!(executed_cases.insert(name.to_owned()));
        }

        let mut mutations = Vec::<(&Path, &str, serde_json::Value)>::new();
        let mut execution_absent = zero_failed.clone();
        let execution_index = unit_index(&execution_absent, "execution");
        execution_absent["units"][execution_index]["diagnostic_code"] =
            serde_json::json!("GPS-CONF-001");
        mutations.push((
            &zero_root,
            "dynamic-diagnostic-execution-absent",
            execution_absent,
        ));

        let mut reduction_absent = one_failed.clone();
        let failed_reduction_index = reduction_absent["units"]
            .as_array()
            .expect("units")
            .iter()
            .position(|unit| unit["kind"] == "reference_reduction" && unit["state"] == "failed")
            .expect("failed reduction unit");
        reduction_absent["units"][failed_reduction_index]["diagnostic_code"] =
            serde_json::json!("GPS-CONF-001");
        mutations.push((
            &one_root,
            "dynamic-diagnostic-reduction-absent",
            reduction_absent,
        ));

        let mut invented = two_failed.clone();
        invented["diagnostics"][0]["code"] = serde_json::json!("GPS-DOMAIN-999");
        for unit in invented["units"].as_array_mut().expect("units") {
            if unit["state"] == "failed" {
                unit["diagnostic_code"] = serde_json::json!("GPS-DOMAIN-999");
            }
        }
        mutations.push((&two_root, "dynamic-diagnostic-invented", invented));

        let mut category = zero_failed.clone();
        category["diagnostics"][0]["category"] = serde_json::json!("CONF");
        mutations.push((&zero_root, "dynamic-diagnostic-category", category));
        let mut meaning = zero_failed.clone();
        meaning["diagnostics"][0]["message"] = serde_json::json!("mensaje inventado");
        mutations.push((&zero_root, "dynamic-diagnostic-meaning", meaning));
        let mut required_context = zero_failed.clone();
        required_context["diagnostics"][0]["context"]
            .as_object_mut()
            .expect("diagnostic context")
            .remove("entity");
        mutations.push((
            &zero_root,
            "dynamic-diagnostic-required-context",
            required_context,
        ));
        let mut duplicate_diagnostic = zero_failed.clone();
        duplicate_diagnostic["diagnostics"]
            .as_array_mut()
            .expect("diagnostics")
            .push(zero_failed["diagnostics"][0].clone());
        mutations.push((
            &zero_root,
            "dynamic-diagnostic-duplicate",
            duplicate_diagnostic,
        ));
        let mut orphan_diagnostic = zero_failed.clone();
        orphan_diagnostic["diagnostics"]
            .as_array_mut()
            .expect("diagnostics")
            .push(serde_json::json!({
                "code": "GPS-CONF-001",
                "category": "CONF",
                "message": "El documento de configuración no tiene sintaxis TOML válida.",
                "context": {"path": "x", "location": "x:1", "detail": "x"}
            }));
        mutations.push((&zero_root, "dynamic-diagnostic-orphan", orphan_diagnostic));

        let mut duplicate_sequence = one_complete.clone();
        duplicate_sequence["references"][0]["sequences"][1]["id"] =
            duplicate_sequence["references"][0]["sequences"][0]["id"].clone();
        mutations.push((
            &one_root,
            "dynamic-reference-sequence-duplicate",
            duplicate_sequence,
        ));
        let mut wrong_reference_total = zero_complete.clone();
        wrong_reference_total["references"][0]["total_genome_size"] =
            serde_json::json!(999_999_999_u64);
        mutations.push((&zero_root, "dynamic-reference-total", wrong_reference_total));
        let mut unknown_count = one_complete.clone();
        let first_length = unknown_count["references"][0]["sequences"][0]["length"]
            .as_u64()
            .expect("sequence length");
        unknown_count["references"][0]["sequences"][0]["unknown_base_count"] =
            serde_json::json!(first_length + 1);
        mutations.push((&one_root, "dynamic-reference-unknown-count", unknown_count));

        let mut partition_overlap = one_complete.clone();
        let first_included =
            partition_overlap["reference_reductions"][0]["included_sequence_ids"][0].clone();
        partition_overlap["reference_reductions"][0]["omitted_sequence_ids"]
            .as_array_mut()
            .expect("omitted")
            .push(first_included);
        mutations.push((
            &one_root,
            "dynamic-reduction-partition-overlap",
            partition_overlap,
        ));
        let mut partition_missing = one_complete.clone();
        partition_missing["reference_reductions"][0]["included_sequence_ids"]
            .as_array_mut()
            .expect("included")
            .pop();
        mutations.push((
            &one_root,
            "dynamic-reduction-partition-missing",
            partition_missing,
        ));
        let mut selection_order = one_complete.clone();
        selection_order["reference_reductions"][0]["selection_order"][1]["selection_order"] =
            serde_json::json!(1);
        mutations.push((
            &one_root,
            "dynamic-reduction-selection-order",
            selection_order,
        ));
        let mut missing_sequence = one_complete.clone();
        missing_sequence["reference_reductions"][0]["selection_order"][0]["sequence_id"] =
            serde_json::json!("missing-sequence");
        mutations.push((
            &one_root,
            "dynamic-reduction-sequence-missing",
            missing_sequence,
        ));
        let mut out_of_range = one_complete.clone();
        out_of_range["reference_reductions"][0]["selection_order"][0]["end"] =
            serde_json::json!(999_999_u64);
        mutations.push((&one_root, "dynamic-reduction-out-of-range", out_of_range));
        let mut overlap = one_complete.clone();
        let first_start = overlap["reference_reductions"][0]["selection_order"][0]["start"].clone();
        let first_end = overlap["reference_reductions"][0]["selection_order"][0]["end"].clone();
        overlap["reference_reductions"][0]["selection_order"][1]["start"] = first_start;
        overlap["reference_reductions"][0]["selection_order"][1]["end"] = first_end;
        mutations.push((&one_root, "dynamic-reduction-overlap", overlap));
        let mut global_sum = one_complete.clone();
        let last_selection = global_sum["reference_reductions"][0]["selection_order"]
            .as_array()
            .expect("selection")
            .len()
            - 1;
        let last_end =
            global_sum["reference_reductions"][0]["selection_order"][last_selection]["end"]
                .as_u64()
                .expect("selection end");
        global_sum["reference_reductions"][0]["selection_order"][last_selection]["end"] =
            serde_json::json!(last_end - 1);
        mutations.push((&one_root, "dynamic-reduction-global-sum", global_sum));

        let mut per_sequence_sum = one_complete.clone();
        let first_derived_length = per_sequence_sum["references"][1]["sequences"][0]["length"]
            .as_u64()
            .expect("derived length");
        let second_derived_length = per_sequence_sum["references"][1]["sequences"][1]["length"]
            .as_u64()
            .expect("derived length");
        per_sequence_sum["references"][1]["sequences"][0]["length"] =
            serde_json::json!(first_derived_length - 1);
        per_sequence_sum["references"][1]["sequences"][1]["length"] =
            serde_json::json!(second_derived_length + 1);
        mutations.push((
            &one_root,
            "dynamic-reduction-per-sequence-sum",
            per_sequence_sum,
        ));
        let mut eligible = one_complete.clone();
        eligible["reference_reductions"][0]["eligible_bases"] = serde_json::json!(21);
        mutations.push((&one_root, "dynamic-reduction-eligible", eligible));
        let mut declared_target = one_complete.clone();
        declared_target["reference_reductions"][0]["target"]["declared_value"] =
            serde_json::json!("12");
        mutations.push((
            &one_root,
            "dynamic-reduction-declared-target",
            declared_target,
        ));
        let mut derived_target = one_complete.clone();
        derived_target["reference_reductions"][0]["target"]["declared_value"] =
            serde_json::json!("12");
        derived_target["reference_reductions"][0]["target"]["derived_bases"] =
            serde_json::json!(12);
        mutations.push((
            &one_root,
            "dynamic-reduction-derived-target",
            derived_target,
        ));
        let mut region = one_complete.clone();
        region["reference_reductions"][0]["region_size_bases"] = serde_json::json!(4);
        mutations.push((&one_root, "dynamic-reduction-region", region));
        let mut residual = one_complete.clone();
        residual["reference_reductions"][0]["has_residual"] = serde_json::json!(false);
        mutations.push((&one_root, "dynamic-reduction-residual", residual));
        let mut derived_class = one_complete.clone();
        derived_class["references"][1]["sequences"][0]["class"] = serde_json::json!("scaffold");
        mutations.push((&one_root, "dynamic-reduction-derived-class", derived_class));
        let mut derived_ploidy = one_complete.clone();
        derived_ploidy["references"][1]["sequences"][0]["ploidy"] = serde_json::json!(3);
        mutations.push((
            &one_root,
            "dynamic-reduction-derived-ploidy",
            derived_ploidy,
        ));
        let mut percent_target = percent_complete;
        percent_target["reference_reductions"][0]["target"]["declared_value"] =
            serde_json::json!("30");
        mutations.push((
            &percent_root,
            "dynamic-reduction-percent-target",
            percent_target,
        ));

        for (root, name, document) in mutations {
            assert_global_manifest_mutation(root, name, &document);
            assert!(executed_cases.insert(name.to_owned()));
        }
        assert_eq!(
            executed_cases,
            MANIFEST_DYNAMIC_CORPUS_CASES
                .iter()
                .map(|name| (*name).to_owned())
                .collect(),
            "the dynamic invariant catalog must equal the cases actually executed"
        );
        for root in [zero_root, one_root, two_root, percent_root] {
            fs::remove_dir_all(root).expect("dynamic graph corpus cleanup");
        }
    }

    #[test]
    fn manifest_bundle_closes_exact_percent_seed_header_normalization_and_event_edges() {
        for (literal, eligible, expected) in [
            ("30.0000000000000000001", 20, 13),
            ("3e1", 20, 14),
            ("29.9999999999999999999", 20, 14),
            ("3.00000000000000000001e1", 20, 13),
        ] {
            assert_eq!(
                independent_retained_after_percent(literal, eligible).unwrap(),
                expected,
                "independent frozen exact-decimal vector {literal}"
            );
        }

        let (above_root, above) =
            percentage_manifest_corpus("30.0000000000000000001", "exhaustive-percent-above");
        let (exponent_root, exponent) =
            percentage_manifest_corpus("3e1", "exhaustive-percent-exponent");
        let (below_root, below) =
            percentage_manifest_corpus("29.9999999999999999999", "exhaustive-percent-below");
        let (non_utf8_root, non_utf8) = non_utf8_header_manifest_corpus("exhaustive-non-utf8");
        let (multiple_root, multiple, multiple_failed) =
            multiple_configuration_and_event_manifest_corpus("exhaustive-multiple");
        let (reduction_root, reduction, _) =
            complete_and_failed_manifest_corpus(1, "exhaustive-reduction");
        let (two_reduction_root, two_reductions, _) =
            complete_and_failed_manifest_corpus(2, "exhaustive-two-reductions");

        let mut executed = BTreeSet::new();
        for (root, name, document) in [
            (
                &above_root,
                "exhaustive-positive-percent-above-integer-boundary",
                &above,
            ),
            (
                &exponent_root,
                "exhaustive-positive-percent-exponent-notation",
                &exponent,
            ),
            (
                &below_root,
                "exhaustive-positive-percent-below-integer-boundary",
                &below,
            ),
            (
                &non_utf8_root,
                "exhaustive-positive-non-utf8-header",
                &non_utf8,
            ),
            (
                &multiple_root,
                "exhaustive-positive-multiple-events-normalizations",
                &multiple,
            ),
            (
                &two_reduction_root,
                "exhaustive-positive-two-reductions",
                &two_reductions,
            ),
        ] {
            assert_global_manifest_positive(root, name, document);
            executed.insert(name.to_owned());
        }

        let reduction_effective_path = reduction_root.join("complete/effective-config.json");
        let reduction_effective: serde_json::Value = serde_json::from_slice(
            &fs::read(&reduction_effective_path).expect("reduction effective"),
        )
        .expect("reduction effective JSON");

        let mut root_manifest = reduction.clone();
        let mut root_effective = reduction_effective.clone();
        root_manifest["reference_reductions"][0]["randomness"]["root_seed"] =
            serde_json::json!("11".repeat(32));
        root_effective["reduced_references"][0]["randomness"]["root_seed"] =
            serde_json::json!("11".repeat(32));
        assert_correlated_bundle_mutation(
            &reduction_root,
            "exhaustive-reduction-root-derived",
            &root_manifest,
            vec![(
                PathBuf::from("effective-config.json"),
                root_effective,
                "effective-config-v1.schema.json",
            )],
        );
        executed.insert("exhaustive-reduction-root-derived".to_owned());

        let mut subseed_manifest = reduction.clone();
        let mut subseed_effective = reduction_effective.clone();
        subseed_manifest["reference_reductions"][0]["randomness"]["subseed"] =
            serde_json::json!("00".repeat(32));
        subseed_effective["reduced_references"][0]["randomness"]["subseed"] =
            serde_json::json!("00".repeat(32));
        assert_correlated_bundle_mutation(
            &reduction_root,
            "exhaustive-reduction-subseed-derived",
            &subseed_manifest,
            vec![(
                PathBuf::from("effective-config.json"),
                subseed_effective,
                "effective-config-v1.schema.json",
            )],
        );
        executed.insert("exhaustive-reduction-subseed-derived".to_owned());

        let mut order = reduction.clone();
        order["references"][1]["sequences"]
            .as_array_mut()
            .expect("derived sequences")
            .swap(0, 1);
        assert_global_manifest_mutation(
            &reduction_root,
            "exhaustive-derived-reference-order",
            &order,
        );
        executed.insert("exhaustive-derived-reference-order".to_owned());

        let non_utf8_effective_path = non_utf8_root.join("complete/effective-config.json");
        let non_utf8_effective: serde_json::Value = serde_json::from_slice(
            &fs::read(&non_utf8_effective_path).expect("non-UTF8 effective"),
        )
        .expect("non-UTF8 effective JSON");
        let mut header_manifest = non_utf8.clone();
        let mut header_effective = non_utf8_effective.clone();
        for reference in [
            &mut header_manifest["references"][0],
            &mut header_effective["reference"],
        ] {
            reference["sequences"][0]["header"] =
                serde_json::json!({"encoding": "utf-8", "value": "other description"});
            reference["sequences"][0]["description"] =
                serde_json::json!({"encoding": "utf-8", "value": "description"});
        }
        assert_correlated_bundle_mutation(
            &non_utf8_root,
            "exhaustive-header-token",
            &header_manifest,
            vec![(
                PathBuf::from("effective-config.json"),
                header_effective,
                "effective-config-v1.schema.json",
            )],
        );
        executed.insert("exhaustive-header-token".to_owned());

        let mut description_manifest = non_utf8.clone();
        let mut description_effective = non_utf8_effective.clone();
        for reference in [
            &mut description_manifest["references"][0],
            &mut description_effective["reference"],
        ] {
            reference["sequences"][0]["description"] =
                serde_json::json!({"encoding": "utf-8", "value": "other"});
        }
        assert_correlated_bundle_mutation(
            &non_utf8_root,
            "exhaustive-header-description",
            &description_manifest,
            vec![(
                PathBuf::from("effective-config.json"),
                description_effective,
                "effective-config-v1.schema.json",
            )],
        );
        executed.insert("exhaustive-header-description".to_owned());

        let mut bytes_description_manifest = non_utf8.clone();
        let mut bytes_description_effective = non_utf8_effective;
        for reference in [
            &mut bytes_description_manifest["references"][0],
            &mut bytes_description_effective["reference"],
        ] {
            reference["sequences"][0]["description"] =
                serde_json::json!({"encoding": "bytes-hex", "value": "00"});
        }
        assert_correlated_bundle_mutation(
            &non_utf8_root,
            "exhaustive-header-description-non-utf8",
            &bytes_description_manifest,
            vec![(
                PathBuf::from("effective-config.json"),
                bytes_description_effective,
                "effective-config-v1.schema.json",
            )],
        );
        executed.insert("exhaustive-header-description-non-utf8".to_owned());

        let multiple_effective_path = multiple_root.join("complete/effective-config.json");
        let multiple_effective: serde_json::Value = serde_json::from_slice(
            &fs::read(&multiple_effective_path).expect("multiple effective"),
        )
        .expect("multiple effective JSON");
        assert!(!multiple["normalizations"].as_array().unwrap().is_empty());
        assert!(multiple["normalizations"].as_array().unwrap().len() >= 3);
        assert!(multiple["events"].as_array().unwrap().len() >= 2);
        assert!(
            multiple_effective["effective_values"]
                .as_array()
                .is_some_and(|values| values.iter().any(|value| {
                    value["object_id"] == "phase-1"
                        && value["winning_rule_id"] == "selfing-generation-5"
                        && value["resolution_context"]["generation"] == 5
                }))
        );

        for (name, manifest) in [
            ("r23-complete-missing-phase-domain-proof", &multiple),
            ("r23-failed-missing-phase-domain-proof", &multiple_failed),
        ] {
            let state = manifest["state"].as_str().expect("multiple state");
            let mut effective: serde_json::Value = serde_json::from_slice(
                &fs::read(multiple_root.join(state).join("effective-config.json"))
                    .expect("multiple effective by state"),
            )
            .expect("multiple effective JSON by state");
            effective["resolution_domain"]["phase_generation_intervals"] = serde_json::json!([]);
            assert_correlated_bundle_mutation(
                &multiple_root,
                name,
                manifest,
                vec![(
                    PathBuf::from("effective-config.json"),
                    effective,
                    "effective-config-v1.schema.json",
                )],
            );
            executed.insert(name.to_owned());
        }

        for (name, manifest) in [
            (
                "r25-complete-correlated-domain-provider-interval",
                &multiple,
            ),
            (
                "r25-failed-correlated-domain-provider-interval",
                &multiple_failed,
            ),
        ] {
            let state = manifest["state"].as_str().expect("multiple state");
            let mut effective: serde_json::Value = serde_json::from_slice(
                &fs::read(multiple_root.join(state).join("effective-config.json"))
                    .expect("multiple effective by state"),
            )
            .expect("multiple effective JSON by state");
            effective["resolution_domain"]["provider"] = serde_json::json!("forged-provider");
            effective["resolution_domain"]["phase_generation_intervals"][0]["start"] =
                serde_json::json!(5);
            assert_correlated_bundle_mutation(
                &multiple_root,
                name,
                manifest,
                vec![(
                    PathBuf::from("effective-config.json"),
                    effective,
                    "effective-config-v1.schema.json",
                )],
            );
            executed.insert(name.to_owned());
        }

        for (name, manifest) in [
            ("r24-complete-scenario-copy-id", &multiple),
            ("r24-failed-scenario-copy-id", &multiple_failed),
        ] {
            let state = manifest["state"].as_str().expect("multiple state");
            let relative = PathBuf::from("inputs/config/scenarios/evolution-small.toml");
            let original = fs::read_to_string(multiple_root.join(state).join(&relative))
                .expect("scenario copy for R24");
            let mutated =
                original.replacen("id = \"evolution-small\"", "id = \"evolution-foreign\"", 1);
            assert_ne!(mutated, original);
            assert_correlated_raw_bundle_mutation(
                &multiple_root,
                name,
                manifest,
                &relative,
                mutated.as_bytes(),
            );
            executed.insert(name.to_owned());
        }

        for (name, manifest) in [
            ("r26-complete-duplicate-phase-across-scenarios", &multiple),
            (
                "r26-failed-duplicate-phase-across-scenarios",
                &multiple_failed,
            ),
        ] {
            let state = manifest["state"].as_str().expect("multiple state");
            let relative = PathBuf::from("inputs/config/scenarios/evolution-other.toml");
            let original = fs::read_to_string(multiple_root.join(state).join(&relative))
                .expect("second scenario copy for R26");
            let mutated = original.replacen("phase-2", "phase-1", 1);
            assert_ne!(mutated, original);
            assert_correlated_raw_bundle_mutation(
                &multiple_root,
                name,
                manifest,
                &relative,
                mutated.as_bytes(),
            );
            executed.insert(name.to_owned());
        }

        let mut normalization_cases = Vec::new();
        for name in [
            "exhaustive-normalization-sum-before",
            "exhaustive-normalization-delta-before",
            "exhaustive-normalization-before",
            "exhaustive-normalization-after",
            "exhaustive-normalization-sum-after",
            "exhaustive-normalization-delta-after",
            "exhaustive-normalization-decision",
            "exhaustive-normalization-order",
            "exhaustive-normalization-duplicate-value",
            "exhaustive-normalization-duplicate-parameter",
        ] {
            let mut manifest = multiple.clone();
            let mut effective = multiple_effective.clone();
            match name {
                "exhaustive-normalization-sum-before" => {
                    manifest["normalizations"][0]["sum_before"] = serde_json::json!(0.5);
                    effective["normalizations"][0]["sum_before"] = serde_json::json!(0.5);
                }
                "exhaustive-normalization-delta-before" => {
                    manifest["normalizations"][0]["delta_before"] = serde_json::json!(0.5);
                    effective["normalizations"][0]["delta_before"] = serde_json::json!(0.5);
                }
                "exhaustive-normalization-before" => {
                    manifest["normalizations"][0]["values"][0]["before"] = serde_json::json!(0.4);
                    effective["normalizations"][0]["values"][0]["before"] = serde_json::json!(0.4);
                }
                "exhaustive-normalization-after" => {
                    manifest["normalizations"][0]["values"][0]["after"] = serde_json::json!(0.4);
                    effective["normalizations"][0]["values"][0]["after"] = serde_json::json!(0.4);
                }
                "exhaustive-normalization-sum-after" => {
                    manifest["normalizations"][0]["sum_after"] = serde_json::json!(0.5);
                    effective["normalizations"][0]["sum_after"] = serde_json::json!(0.5);
                }
                "exhaustive-normalization-delta-after" => {
                    manifest["normalizations"][0]["delta_after"] = serde_json::json!(0.5);
                    effective["normalizations"][0]["delta_after"] = serde_json::json!(0.5);
                }
                "exhaustive-normalization-decision" => {
                    let replacement = if manifest["normalizations"][0]["decision"] == "normalized" {
                        "accepted_without_change"
                    } else {
                        "normalized"
                    };
                    manifest["normalizations"][0]["decision"] = serde_json::json!(replacement);
                    effective["normalizations"][0]["decision"] = serde_json::json!(replacement);
                }
                "exhaustive-normalization-order" => {
                    manifest["normalizations"][0]["values"]
                        .as_array_mut()
                        .unwrap()
                        .swap(0, 1);
                    effective["normalizations"][0]["values"]
                        .as_array_mut()
                        .unwrap()
                        .swap(0, 1);
                }
                "exhaustive-normalization-duplicate-value" => {
                    let duplicate = manifest["normalizations"][0]["values"][0].clone();
                    manifest["normalizations"][0]["values"]
                        .as_array_mut()
                        .unwrap()
                        .push(duplicate.clone());
                    effective["normalizations"][0]["values"]
                        .as_array_mut()
                        .unwrap()
                        .push(duplicate);
                }
                "exhaustive-normalization-duplicate-parameter" => {
                    let duplicate = manifest["normalizations"][0].clone();
                    manifest["normalizations"]
                        .as_array_mut()
                        .unwrap()
                        .push(duplicate.clone());
                    effective["normalizations"]
                        .as_array_mut()
                        .unwrap()
                        .push(duplicate);
                }
                _ => unreachable!(),
            }
            normalization_cases.push((name, manifest, effective));
        }
        for (name, manifest, effective) in normalization_cases {
            assert_correlated_bundle_mutation(
                &multiple_root,
                name,
                &manifest,
                vec![(
                    PathBuf::from("effective-config.json"),
                    effective,
                    "effective-config-v1.schema.json",
                )],
            );
            executed.insert(name.to_owned());
        }

        let first_event_id = multiple["events"][0]["event_id"]
            .as_str()
            .expect("first event ID")
            .to_owned();
        let first_event_path = PathBuf::from(
            multiple["events"][0]["manifest"]["path"]
                .as_str()
                .expect("first event path"),
        );
        let event_document: serde_json::Value = serde_json::from_slice(
            &fs::read(multiple_root.join("complete").join(&first_event_path))
                .expect("event document"),
        )
        .expect("event JSON");

        let mut general_event = multiple.clone();
        general_event["events"][0]["subseed"] = serde_json::json!("00".repeat(32));
        assert_global_manifest_mutation(
            &multiple_root,
            "exhaustive-event-general-subseed",
            &general_event,
        );
        executed.insert("exhaustive-event-general-subseed".to_owned());

        let mut effective_event = multiple_effective.clone();
        effective_event["randomness"]["event_subseeds"][0]["subseed"] =
            serde_json::json!("00".repeat(32));
        assert_correlated_bundle_mutation(
            &multiple_root,
            "exhaustive-event-effective-subseed",
            &multiple,
            vec![(
                PathBuf::from("effective-config.json"),
                effective_event,
                "effective-config-v1.schema.json",
            )],
        );
        executed.insert("exhaustive-event-effective-subseed".to_owned());

        let mut event_manifest = event_document.clone();
        event_manifest["subseed"] = serde_json::json!("00".repeat(32));
        assert_correlated_bundle_mutation(
            &multiple_root,
            "exhaustive-event-manifest-subseed",
            &multiple,
            vec![(
                first_event_path.clone(),
                event_manifest,
                "manifest-event-v1.schema.json",
            )],
        );
        executed.insert("exhaustive-event-manifest-subseed".to_owned());

        let mut correlated_manifest = multiple.clone();
        correlated_manifest["events"][0]["subseed"] = serde_json::json!("00".repeat(32));
        let mut correlated_effective = multiple_effective;
        let effective_event_index = correlated_effective["randomness"]["event_subseeds"]
            .as_array()
            .unwrap()
            .iter()
            .position(|event| event["event_id"] == first_event_id)
            .expect("matching effective event");
        correlated_effective["randomness"]["event_subseeds"][effective_event_index]["subseed"] =
            serde_json::json!("00".repeat(32));
        let mut correlated_event = event_document;
        correlated_event["subseed"] = serde_json::json!("00".repeat(32));
        assert_correlated_bundle_mutation(
            &multiple_root,
            "exhaustive-event-correlated-derived-subseed",
            &correlated_manifest,
            vec![
                (
                    PathBuf::from("effective-config.json"),
                    correlated_effective,
                    "effective-config-v1.schema.json",
                ),
                (
                    first_event_path,
                    correlated_event,
                    "manifest-event-v1.schema.json",
                ),
            ],
        );
        executed.insert("exhaustive-event-correlated-derived-subseed".to_owned());

        assert_eq!(
            executed,
            MANIFEST_EXHAUSTIVE_CORPUS_CASES
                .iter()
                .chain(MANIFEST_R23_R24_CORPUS_CASES)
                .chain(MANIFEST_R26_CORPUS_CASES)
                .chain(
                    MANIFEST_R25_CORPUS_CASES
                        .iter()
                        .filter(|name| { **name != "r25-correlated-global-output-to-event" })
                )
                .map(|name| (*name).to_owned())
                .collect(),
            "the exhaustive dynamic catalog must equal the cases actually executed"
        );
        for root in [
            above_root,
            exponent_root,
            below_root,
            non_utf8_root,
            multiple_root,
            reduction_root,
            two_reduction_root,
        ] {
            fs::remove_dir_all(root).expect("exhaustive corpus cleanup");
        }
    }

    #[test]
    fn manifest_bundle_closes_toml_fasta_and_effective_value_edges() {
        let (separator_root, separator) =
            percentage_manifest_corpus("1_0.0_0", "r17-percent-separators");
        let (plus_root, plus) = percentage_manifest_corpus("+1.25e1", "r17-percent-plus-exponent");
        let (precision_root, precision) =
            percentage_manifest_corpus("30.000_000_000_000_000_000_1", "r17-percent-precision");
        let (non_utf8_root, non_utf8) = non_utf8_header_manifest_corpus("r17-fasta-non-utf8");
        let (supplied_root, supplied) =
            supplied_two_sequence_manifest_corpus("r17-fasta-two-sequences");
        let (effective_root, effective_manifest, effective) =
            scoped_effective_manifest_corpus("r17-effective-scopes");
        let (two_root, two, _) = complete_and_failed_manifest_corpus(2, "r17-two-reductions");
        let (parameter_order_root, parameter_order_complete, parameter_order_failed) =
            parameter_order_manifest_corpus("r27-parameter-order");
        let mut executed = BTreeSet::new();

        for (root, name, document) in [
            (
                &separator_root,
                "r17-positive-percent-separators",
                &separator,
            ),
            (&plus_root, "r17-positive-percent-plus-exponent", &plus),
            (
                &precision_root,
                "r17-positive-percent-high-precision",
                &precision,
            ),
            (&non_utf8_root, "r17-positive-fasta-non-utf8", &non_utf8),
            (
                &effective_root,
                "r17-positive-effective-multiple-scopes",
                &effective_manifest,
            ),
            (&two_root, "r17-positive-two-reductions", &two),
        ] {
            assert_global_manifest_positive(root, name, document);
            executed.insert(name.to_owned());
        }

        let separator_effective: serde_json::Value = serde_json::from_slice(
            &fs::read(separator_root.join("complete/effective-config.json"))
                .expect("separator effective"),
        )
        .expect("separator effective JSON");
        for (name, literal) in [
            ("r17-percent-leading-separator", "_10"),
            ("r17-percent-trailing-separator", "10_"),
            ("r17-percent-double-separator", "1__0"),
            ("r17-percent-separator-before-dot", "1_.0"),
            ("r17-percent-separator-after-dot", "1._0"),
            ("r17-percent-missing-integer", ".5"),
            ("r17-percent-missing-fraction", "5."),
            ("r17-percent-leading-zero-integer", "00"),
            ("r17-percent-leading-zero-float", "01.0"),
            ("r17-percent-signed-leading-zero", "+01"),
            ("r17-percent-empty-exponent", "1e"),
            ("r17-percent-empty-signed-exponent", "1e+"),
            ("r17-percent-leading-exponent-separator", "1e_2"),
            ("r17-percent-trailing-exponent-separator", "1e2_"),
        ] {
            let mut manifest = separator.clone();
            let mut linked = separator_effective.clone();
            manifest["reference_reductions"][0]["target"]["declared_value"] =
                serde_json::json!(literal);
            linked["reduced_references"][0]["target"]["declared_value"] =
                serde_json::json!(literal);
            assert_correlated_bundle_mutation(
                &separator_root,
                name,
                &manifest,
                vec![(
                    PathBuf::from("effective-config.json"),
                    linked,
                    "effective-config-v1.schema.json",
                )],
            );
            executed.insert(name.to_owned());
        }

        let non_utf8_effective: serde_json::Value = serde_json::from_slice(
            &fs::read(non_utf8_root.join("complete/effective-config.json"))
                .expect("non-UTF8 effective"),
        )
        .expect("non-UTF8 effective JSON");
        for name in [
            "r17-fasta-correlated-id-header",
            "r17-fasta-correlated-description",
            "r17-fasta-correlated-length",
            "r17-fasta-correlated-unknown-count",
        ] {
            let mut manifest = non_utf8.clone();
            let mut linked = non_utf8_effective.clone();
            match name {
                "r17-fasta-correlated-id-header" => {
                    for reference in [&mut manifest["references"][0], &mut linked["reference"]] {
                        reference["sequences"][0]["id"] = serde_json::json!("other");
                        reference["sequences"][0]["header"] = serde_json::json!({
                            "encoding": "utf-8", "value": "other altered"
                        });
                        reference["sequences"][0]["description"] = serde_json::json!({
                            "encoding": "utf-8", "value": "altered"
                        });
                    }
                }
                "r17-fasta-correlated-description" => {
                    for reference in [&mut manifest["references"][0], &mut linked["reference"]] {
                        reference["sequences"][0]["header"] = serde_json::json!({
                            "encoding": "utf-8", "value": "chr1 altered"
                        });
                        reference["sequences"][0]["description"] = serde_json::json!({
                            "encoding": "utf-8", "value": "altered"
                        });
                    }
                }
                "r17-fasta-correlated-length" => {
                    for reference in [&mut manifest["references"][0], &mut linked["reference"]] {
                        reference["sequences"][0]["length"] = serde_json::json!(6);
                        reference["total_genome_size"] = serde_json::json!(6);
                    }
                }
                "r17-fasta-correlated-unknown-count" => {
                    for reference in [&mut manifest["references"][0], &mut linked["reference"]] {
                        reference["sequences"][0]["unknown_base_count"] = serde_json::json!(0);
                    }
                }
                _ => unreachable!(),
            }
            assert_correlated_bundle_mutation(
                &non_utf8_root,
                name,
                &manifest,
                vec![(
                    PathBuf::from("effective-config.json"),
                    linked,
                    "effective-config-v1.schema.json",
                )],
            );
            executed.insert(name.to_owned());
        }

        let supplied_effective: serde_json::Value = serde_json::from_slice(
            &fs::read(supplied_root.join("complete/effective-config.json"))
                .expect("supplied effective"),
        )
        .expect("supplied effective JSON");
        let mut order_manifest = supplied.clone();
        let mut order_effective = supplied_effective;
        order_manifest["references"][0]["sequences"]
            .as_array_mut()
            .expect("manifest sequences")
            .swap(0, 1);
        order_effective["reference"]["sequences"]
            .as_array_mut()
            .expect("effective sequences")
            .swap(0, 1);
        assert_correlated_bundle_mutation(
            &supplied_root,
            "r17-fasta-correlated-order",
            &order_manifest,
            vec![(
                PathBuf::from("effective-config.json"),
                order_effective,
                "effective-config-v1.schema.json",
            )],
        );
        executed.insert("r17-fasta-correlated-order".to_owned());

        let references = supplied_root.join("complete/resources/references");
        let real_references = supplied_root.join("complete/resources/references-real");
        fs::rename(&references, &real_references).expect("move real reference directory");
        std::os::unix::fs::symlink("references-real", &references)
            .expect("intermediate reference symlink");
        let symlink_path =
            write_manifest_case(&supplied_root, "r17-fasta-intermediate-symlink", &supplied);
        let symlink_bytes = fs::read(&symlink_path).expect("symlink manifest");
        assert!(
            fs::read(supplied_root.join("complete/resources/references/ref-main.fa")).is_ok(),
            "ordinary path traversal must still reach the symlink target"
        );
        assert!(
            validate_general_manifest_bytes(&symlink_path, &symlink_bytes).is_err(),
            "production must reject an intermediate FASTA symlink"
        );
        fs::remove_file(&references).expect("remove test symlink");
        fs::rename(&real_references, &references).expect("restore reference directory");
        executed.insert("r17-fasta-intermediate-symlink".to_owned());

        let values = effective["effective_values"]
            .as_array()
            .expect("effective values");
        assert!(values.len() >= 3, "multiple effective scopes are required");
        let value_index = |object_id: &str| {
            values
                .iter()
                .position(|value| value["object_id"] == object_id)
                .unwrap_or_else(|| panic!("missing effective value for {object_id}"))
        };
        let global_index = value_index("global");
        let individual_index = value_index("ind-1");
        assert!(
            values.iter().any(|value| {
                value["object_id"] == "event-1"
                    && value["parameter"] == "coverage"
                    && value["winning_rule_id"] == "coverage-base"
                    && value["resolution_context"]["technical_event"] == "event-1"
            }),
            "the executable corpus needs a positive technical-axis witness"
        );
        let ploidy_index = values
            .iter()
            .position(|value| value["parameter"] == "ploidy")
            .expect("ploidy effective value");
        let coverage_index = values
            .iter()
            .position(|value| value["parameter"] == "coverage")
            .expect("coverage effective value");
        let recombination_index = values
            .iter()
            .position(|value| value["parameter"] == "recombination_rate")
            .expect("recombination effective value");
        let output_index = values
            .iter()
            .position(|value| {
                value["object_id"] == "global" && value["parameter"] == "output_enabled"
            })
            .expect("global output effective value");
        let mut effective_cases = Vec::new();

        let mut forged_output_request = effective.clone();
        forged_output_request["effective_values"][output_index]["object_id"] =
            serde_json::json!("event-1");
        forged_output_request["effective_values"][output_index]["resolution_context"] = serde_json::json!({
            "individual": null,
            "population": null,
            "subpopulation": null,
            "phase": null,
            "generation": null,
            "reference": null,
            "sequence": null,
            "position": null,
            "locus": null,
            "group": null,
            "block": null,
            "technical_event": "event-1",
        });
        effective_cases.push((
            "r25-correlated-global-output-to-event",
            effective_manifest.clone(),
            forged_output_request,
        ));

        let mut duplicate = effective.clone();
        duplicate["effective_values"]
            .as_array_mut()
            .unwrap()
            .push(effective["effective_values"][individual_index].clone());
        effective_cases.push((
            "r17-effective-duplicate-identical",
            effective_manifest.clone(),
            duplicate,
        ));

        let mut contradictory = effective.clone();
        let mut contradictory_value = effective["effective_values"][individual_index].clone();
        contradictory_value["value"] = serde_json::json!(0.03);
        contradictory["effective_values"]
            .as_array_mut()
            .unwrap()
            .push(contradictory_value);
        effective_cases.push((
            "r17-effective-duplicate-contradictory",
            effective_manifest.clone(),
            contradictory,
        ));

        let mut reordered = effective.clone();
        reordered["effective_values"]
            .as_array_mut()
            .unwrap()
            .swap(0, 1);
        effective_cases.push((
            "r17-effective-reordered",
            effective_manifest.clone(),
            reordered,
        ));

        let mut parameter = effective.clone();
        parameter["effective_values"][individual_index]["parameter"] =
            serde_json::json!("coverage");
        effective_cases.push((
            "r17-effective-parameter-cross",
            effective_manifest.clone(),
            parameter,
        ));

        let mut crossed_value = effective.clone();
        crossed_value["effective_values"][individual_index]["value"] = serde_json::json!(0.03);
        effective_cases.push((
            "r17-effective-value-cross",
            effective_manifest.clone(),
            crossed_value,
        ));

        let mut winner = effective.clone();
        winner["effective_values"][individual_index]["winning_rule_id"] =
            effective["effective_values"][global_index]["winning_rule_id"].clone();
        effective_cases.push((
            "r17-effective-winner-cross",
            effective_manifest.clone(),
            winner,
        ));

        let mut selectors = effective.clone();
        selectors["effective_values"][individual_index]["selectors"] =
            effective["effective_values"][global_index]["selectors"].clone();
        effective_cases.push((
            "r17-effective-selectors-cross",
            effective_manifest.clone(),
            selectors,
        ));

        let mut provenance_path = effective.clone();
        provenance_path["effective_values"][individual_index]["provenance"]["path"] =
            serde_json::json!("other.toml");
        effective_cases.push((
            "r17-effective-provenance-path",
            effective_manifest.clone(),
            provenance_path,
        ));

        let mut provenance_location = effective.clone();
        provenance_location["effective_values"][individual_index]["provenance"]["location"] =
            serde_json::json!("line 1, column 1, bytes 0..1");
        effective_cases.push((
            "r17-effective-provenance-location",
            effective_manifest.clone(),
            provenance_location,
        ));

        let mut missing = effective.clone();
        missing["effective_values"]
            .as_array_mut()
            .unwrap()
            .remove(global_index);
        effective_cases.push((
            "r17-effective-automatic-missing",
            effective_manifest.clone(),
            missing,
        ));

        let mut ploidy_manifest = effective_manifest.clone();
        let mut ploidy = effective.clone();
        ploidy["effective_values"][ploidy_index]["value"] = serde_json::json!(3);
        ploidy["reference"]["sequences"][0]["ploidy"] = serde_json::json!(3);
        ploidy_manifest["references"][0]["sequences"][0]["ploidy"] = serde_json::json!(3);
        effective_cases.push(("r17-effective-correlated-ploidy", ploidy_manifest, ploidy));

        let mut source_path = effective.clone();
        source_path["source_path"] = serde_json::json!("other.toml");
        effective_cases.push((
            "r17-effective-source-path",
            effective_manifest.clone(),
            source_path,
        ));

        let mut controlled_duplicate = effective.clone();
        let controlled = controlled_duplicate["controlled_frequency_sets"][0].clone();
        controlled_duplicate["controlled_frequency_sets"]
            .as_array_mut()
            .unwrap()
            .push(controlled);
        effective_cases.push((
            "r17-effective-controlled-duplicate",
            effective_manifest.clone(),
            controlled_duplicate,
        ));

        let mut controlled_count = effective.clone();
        controlled_count["controlled_frequency_sets"][0]["counts"][0]["count"] =
            serde_json::json!(99);
        effective_cases.push((
            "r17-effective-controlled-count",
            effective_manifest.clone(),
            controlled_count,
        ));

        let general_copy = effective_manifest["configuration_inputs"]
            .as_array()
            .expect("configuration inputs")
            .iter()
            .find(|input| input["kind"] == "general")
            .expect("general input");
        let general_copy_path = general_copy["copy"]["path"]
            .as_str()
            .expect("general copy path");
        let declared_rules = independent_rule_catalog(
            &fs::read(effective_root.join("complete").join(general_copy_path))
                .expect("general copy for R18"),
        )
        .expect("R18 independent rule catalog");
        let base_rule = declared_rules
            .get("variant-base")
            .expect("variant base rule");
        let base_effective_value = |object_id: &str| {
            serde_json::json!({
                "object_id": object_id,
                "parameter": base_rule.parameter,
                "resolution_context": effective["effective_values"][individual_index]["resolution_context"],
                "value": base_rule.value,
                "winning_rule_id": "variant-base",
                "selectors": base_rule.selectors,
                "provenance": {
                    "path": effective["source_path"],
                    "location": base_rule.location,
                },
            })
        };
        let mut dominated_base = effective.clone();
        let copied_base = base_effective_value("ind-1");
        dominated_base["effective_values"][individual_index] = copied_base;
        effective_cases.push((
            "r18-effective-dominated-base-copy",
            effective_manifest.clone(),
            dominated_base,
        ));

        let mut foreign_object = effective.clone();
        let copied_foreign = base_effective_value("zzz-foreign");
        foreign_object["effective_values"][individual_index] = copied_foreign;
        effective_cases.push((
            "r18-effective-foreign-object",
            effective_manifest.clone(),
            foreign_object,
        ));

        let mut forbidden_axis = effective.clone();
        let mut forbidden_value = base_effective_value("event-1");
        forbidden_value["resolution_context"] = serde_json::json!({
            "individual": null,
            "population": null,
            "subpopulation": null,
            "phase": null,
            "generation": null,
            "reference": null,
            "sequence": null,
            "position": null,
            "locus": null,
            "group": null,
            "block": null,
            "technical_event": "event-1",
        });
        forbidden_axis["effective_values"][individual_index] = forbidden_value;
        effective_cases.push((
            "r19-biological-parameter-technical-context",
            effective_manifest.clone(),
            forbidden_axis,
        ));

        let mut orphan_position = effective.clone();
        orphan_position["effective_values"][individual_index]["resolution_context"]["position"] =
            serde_json::json!(15);
        effective_cases.push((
            "r21-orphan-position-with-individual",
            effective_manifest.clone(),
            orphan_position,
        ));

        let mut coverage_locus = effective.clone();
        coverage_locus["effective_values"][coverage_index]["resolution_context"] = serde_json::json!({
            "individual": null,
            "population": null,
            "subpopulation": null,
            "phase": null,
            "generation": null,
            "reference": "ref-main",
            "sequence": "chr1",
            "position": 15,
            "locus": "locus-15",
            "group": null,
            "block": null,
            "technical_event": null,
        });
        effective_cases.push((
            "r22-coverage-locus",
            effective_manifest.clone(),
            coverage_locus,
        ));

        let mut recombination_region = effective.clone();
        recombination_region["effective_values"][recombination_index]["resolution_context"]["position"] =
            serde_json::json!(15);
        effective_cases.push((
            "r22-recombination-region",
            effective_manifest.clone(),
            recombination_region,
        ));

        let mut ploidy_locus = effective.clone();
        ploidy_locus["effective_values"][ploidy_index]["resolution_context"]["position"] =
            serde_json::json!(15);
        ploidy_locus["effective_values"][ploidy_index]["resolution_context"]["locus"] =
            serde_json::json!("locus-15");
        effective_cases.push(("r22-ploidy-locus", effective_manifest.clone(), ploidy_locus));

        let mut reference_object = effective.clone();
        reference_object["effective_values"][individual_index]["object_id"] =
            serde_json::json!("ref-main");
        reference_object["effective_values"][individual_index]["resolution_context"] = serde_json::json!({
            "individual": null,
            "population": null,
            "subpopulation": null,
            "phase": null,
            "generation": null,
            "reference": "ref-main",
            "sequence": "chr1",
            "position": null,
            "locus": null,
            "group": null,
            "block": null,
            "technical_event": null,
        });
        effective_cases.push((
            "r22-reference-object",
            effective_manifest.clone(),
            reference_object,
        ));

        for (name, manifest, linked) in effective_cases {
            assert_correlated_bundle_mutation(
                &effective_root,
                name,
                &manifest,
                vec![(
                    PathBuf::from("effective-config.json"),
                    linked,
                    "effective-config-v1.schema.json",
                )],
            );
            executed.insert(name.to_owned());
        }

        for (name, manifest) in [
            ("r27-complete-parameter-order", &parameter_order_complete),
            ("r27-failed-parameter-order", &parameter_order_failed),
        ] {
            let state = manifest["state"].as_str().expect("parameter-order state");
            let mut linked: serde_json::Value = serde_json::from_slice(
                &fs::read(
                    parameter_order_root
                        .join(state)
                        .join("effective-config.json"),
                )
                .expect("parameter-order effective config"),
            )
            .expect("parameter-order effective JSON");
            assert_eq!(linked["parameter_order_version"], serde_json::json!("1"));
            let values = linked["effective_values"]
                .as_array_mut()
                .expect("parameter-order values");
            let variant_index = values
                .iter()
                .position(|value| {
                    value["object_id"] == "ind-1" && value["parameter"] == "variant_rate"
                })
                .expect("individual variant value");
            let output_index = values
                .iter()
                .position(|value| {
                    value["object_id"] == "ind-1" && value["parameter"] == "output_enabled"
                })
                .expect("individual output value");
            assert!(
                variant_index < output_index,
                "canonical literal order witness"
            );
            values.swap(variant_index, output_index);
            assert_correlated_bundle_mutation(
                &parameter_order_root,
                name,
                manifest,
                vec![(
                    PathBuf::from("effective-config.json"),
                    linked,
                    "effective-config-v1.schema.json",
                )],
            );
            executed.insert(name.to_owned());
        }

        assert_eq!(
            executed,
            MANIFEST_R17_CORPUS_CASES
                .iter()
                .chain(MANIFEST_R18_CORPUS_CASES)
                .chain(MANIFEST_R19_R21_CORPUS_CASES)
                .chain(MANIFEST_R22_CORPUS_CASES)
                .chain(
                    MANIFEST_R25_CORPUS_CASES
                        .iter()
                        .filter(|name| { **name == "r25-correlated-global-output-to-event" })
                )
                .chain(MANIFEST_R27_CORPUS_CASES)
                .map(|name| (*name).to_owned())
                .collect(),
            "R17-R21 catalogs must equal the executable corpus"
        );
        for root in [
            separator_root,
            plus_root,
            precision_root,
            non_utf8_root,
            supplied_root,
            effective_root,
            two_root,
            parameter_order_root,
        ] {
            fs::remove_dir_all(root).expect("R17 corpus cleanup");
        }
    }

    #[test]
    fn global_manifest_links_must_match_catalog_exactly_without_reductions() {
        let (root, complete, _) = zero_reduction_manifest_corpus("global-links-regression");
        let mut reference_digest = complete.clone();
        reference_digest["references"][0]["fasta"]["checksum"]["digest"] =
            serde_json::json!("0".repeat(64));
        let mut authoritative_resource = complete.clone();
        authoritative_resource["references"][0]["authoritative_resource"] =
            serde_json::json!("resources/references/other.fa");
        let mut event_digest = complete.clone();
        event_digest["events"][0]["manifest"]["checksum"]["digest"] =
            serde_json::json!("0".repeat(64));
        let mut effective_config_digest = complete;
        effective_config_digest["effective_config"]["checksum"]["digest"] =
            serde_json::json!("0".repeat(64));

        for (name, document) in [
            ("reference-digest", reference_digest),
            ("authoritative-resource", authoritative_resource),
            ("event-digest", event_digest),
            ("effective-config-digest", effective_config_digest),
        ] {
            let path = write_manifest_case(&root, name, &document);
            let bytes = fs::read(&path).expect("global-link mutation bytes");
            assert!(
                python_schema_accepts(&path, "manifest-general-v1.schema.json"),
                "the individual cross-collection mutation remains structurally valid: {name}"
            );
            let error = validate_general_manifest_bytes(&path, &bytes)
                .expect_err("the semantic gate must reject a non-catalog global link");
            assert_eq!(error.code(), "GPS-FMT-006", "case={name}");
        }
        fs::remove_dir_all(root).expect("global-link regression cleanup");
    }

    fn assert_global_manifest_positive(root: &Path, name: &str, document: &serde_json::Value) {
        let path = write_manifest_case(root, name, document);
        let bytes = fs::read(&path).expect("positive global-link bytes");
        assert!(
            python_schema_accepts(&path, "manifest-general-v1.schema.json"),
            "Python rejected positive {name}"
        );
        independent_semantic_manifest_oracle(&path, &bytes)
            .unwrap_or_else(|error| panic!("independent oracle rejected {name}: {error}"));
        validate_general_manifest_bytes(&path, &bytes)
            .unwrap_or_else(|error| panic!("Rust rejected {name}: {error:?}"));
    }

    fn assert_global_manifest_mutation(root: &Path, name: &str, document: &serde_json::Value) {
        let path = write_manifest_case(root, name, document);
        let bytes = fs::read(&path).expect("mutated global-link bytes");
        assert!(
            python_schema_accepts(&path, "manifest-general-v1.schema.json"),
            "individual semantic mutation must remain structurally valid: {name}"
        );
        assert!(
            independent_semantic_manifest_oracle(&path, &bytes).is_err(),
            "independent semantic oracle accepted {name}"
        );
        let error = validate_general_manifest_bytes(&path, &bytes)
            .expect_err("Rust semantic gate must reject the individual mutation");
        assert_eq!(error.code(), "GPS-FMT-006", "case={name}");
    }

    fn assert_correlated_bundle_mutation(
        root: &Path,
        name: &str,
        document: &serde_json::Value,
        linked_documents: Vec<(PathBuf, serde_json::Value, &'static str)>,
    ) {
        let state_directory = document["state"].as_str().expect("manifest state");
        let base = root.join(state_directory);
        let mut document = document.clone();
        let mut originals = Vec::new();
        for (relative, linked, schema) in &linked_documents {
            let path = base.join(relative);
            originals.push((path.clone(), fs::read(&path).expect("linked original")));
            let linked_bytes =
                serde_json::to_vec_pretty(linked).expect("linked mutation serialization");
            fs::write(&path, &linked_bytes).expect("linked mutation write");
            let mut identity = crate::integrity::Sha256Writer::default();
            identity
                .update(&linked_bytes)
                .expect("linked mutation digest");
            let identity = identity.finish(relative);
            rebind_artifact_identity(
                &mut document,
                relative.to_str().expect("linked relative UTF-8"),
                identity.size_bytes(),
                identity.digest(),
            );
            assert!(
                python_schema_accepts(&path, schema),
                "linked mutation must remain structurally valid: {name}/{schema}"
            );
        }
        let path = write_manifest_case(root, name, &document);
        let bytes = fs::read(&path).expect("correlated manifest bytes");
        assert!(
            python_schema_accepts(&path, "manifest-general-v1.schema.json"),
            "correlated mutation must preserve manifest structure: {name}"
        );
        assert!(
            independent_semantic_manifest_oracle(&path, &bytes).is_err(),
            "independent oracle accepted correlated mutation {name}"
        );
        let error = validate_general_manifest_bytes(&path, &bytes)
            .expect_err("production accepted correlated bundle mutation");
        assert_eq!(error.code(), "GPS-FMT-006", "case={name}");
        for (path, bytes) in originals {
            fs::write(path, bytes).expect("linked document restore");
        }
    }

    fn assert_correlated_raw_bundle_mutation(
        root: &Path,
        name: &str,
        document: &serde_json::Value,
        relative: &Path,
        linked_bytes: &[u8],
    ) {
        let state_directory = document["state"].as_str().expect("manifest state");
        let path = root.join(state_directory).join(relative);
        let original = fs::read(&path).expect("raw linked original");
        fs::write(&path, linked_bytes).expect("raw linked mutation write");
        let mut identity = crate::integrity::Sha256Writer::default();
        identity.update(linked_bytes).expect("raw mutation digest");
        let identity = identity.finish(relative);
        let mut document = document.clone();
        rebind_artifact_identity(
            &mut document,
            relative.to_str().expect("raw relative UTF-8"),
            identity.size_bytes(),
            identity.digest(),
        );
        let manifest_path = write_manifest_case(root, name, &document);
        let bytes = fs::read(&manifest_path).expect("raw correlated manifest bytes");
        assert!(
            python_schema_accepts(&manifest_path, "manifest-general-v1.schema.json"),
            "raw correlated mutation must preserve manifest structure: {name}"
        );
        assert!(
            independent_semantic_manifest_oracle(&manifest_path, &bytes).is_err(),
            "independent oracle accepted raw correlated mutation {name}"
        );
        let error = validate_general_manifest_bytes(&manifest_path, &bytes)
            .expect_err("production accepted raw correlated bundle mutation");
        assert_eq!(error.code(), "GPS-FMT-006", "case={name}");
        fs::write(path, original).expect("raw linked document restore");
    }

    fn rebind_artifact_identity(
        value: &mut serde_json::Value,
        path: &str,
        size_bytes: u64,
        digest: &str,
    ) {
        match value {
            serde_json::Value::Array(values) => {
                for value in values {
                    rebind_artifact_identity(value, path, size_bytes, digest);
                }
            }
            serde_json::Value::Object(object) => {
                if object.get("path").and_then(serde_json::Value::as_str) == Some(path)
                    && object.contains_key("checksum")
                {
                    object.insert("size_bytes".to_owned(), serde_json::json!(size_bytes));
                    object.insert(
                        "checksum".to_owned(),
                        serde_json::json!({"algorithm": "sha256", "digest": digest}),
                    );
                }
                for value in object.values_mut() {
                    rebind_artifact_identity(value, path, size_bytes, digest);
                }
            }
            _ => {}
        }
    }

    #[test]
    fn global_manifest_links_match_catalog_for_complete_failed_zero_one_and_two_reductions() {
        let (zero_root, zero_complete, zero_failed) =
            zero_reduction_manifest_corpus("global-links-zero");
        let (one_root, one_complete, one_failed) =
            complete_and_failed_manifest_corpus(1, "global-links-one");
        let (two_root, two_complete, two_failed) =
            complete_and_failed_manifest_corpus(2, "global-links-two");

        let mut executed_cases = BTreeSet::new();
        for (root, name, document) in [
            (&zero_root, "positive-zero-complete", &zero_complete),
            (&zero_root, "positive-zero-failed", &zero_failed),
            (&one_root, "positive-one-complete", &one_complete),
            (&one_root, "positive-one-failed", &one_failed),
            (&two_root, "positive-two-complete", &two_complete),
            (&two_root, "positive-two-failed", &two_failed),
        ] {
            assert_global_manifest_positive(root, name, document);
            assert!(executed_cases.insert(name.to_owned()));
        }

        let mut mutations = Vec::new();
        for (prefix, root, document) in [
            ("zero-complete", &zero_root, &zero_complete),
            ("one-complete", &one_root, &one_complete),
            ("one-failed", &one_root, &one_failed),
            ("two-complete", &two_root, &two_complete),
            ("two-failed", &two_root, &two_failed),
        ] {
            let mut reference_digest = document.clone();
            reference_digest["references"][0]["fasta"]["checksum"]["digest"] =
                serde_json::json!("0".repeat(64));
            mutations.push((root, format!("{prefix}-reference-digest"), reference_digest));
            let mut authoritative_resource = document.clone();
            authoritative_resource["references"][0]["authoritative_resource"] =
                serde_json::json!("resources/references/other.fa");
            mutations.push((
                root,
                format!("{prefix}-authoritative-resource"),
                authoritative_resource,
            ));
            let mut effective_config_digest = document.clone();
            effective_config_digest["effective_config"]["checksum"]["digest"] =
                serde_json::json!("0".repeat(64));
            mutations.push((
                root,
                format!("{prefix}-effective-config-digest"),
                effective_config_digest,
            ));
        }
        let mut zero_failed_effective = zero_failed.clone();
        zero_failed_effective["effective_config"]["checksum"]["digest"] =
            serde_json::json!("0".repeat(64));
        mutations.push((
            &zero_root,
            "zero-failed-effective-config-digest".to_owned(),
            zero_failed_effective,
        ));
        for (name, pointer, value) in [
            (
                "zero-reference-id",
                "/references/0/fasta/id",
                serde_json::json!("reference-fasta-other"),
            ),
            (
                "zero-reference-kind",
                "/references/0/fasta/kind",
                serde_json::json!("other"),
            ),
            (
                "zero-reference-path",
                "/references/0/fasta/path",
                serde_json::json!("resources/references/other.fa"),
            ),
            (
                "zero-reference-size",
                "/references/0/fasta/size_bytes",
                serde_json::json!(999_999_999_u64),
            ),
            (
                "zero-effective-config-id",
                "/effective_config/id",
                serde_json::json!("effective-config-other"),
            ),
            (
                "zero-effective-config-kind",
                "/effective_config/kind",
                serde_json::json!("other"),
            ),
            (
                "zero-effective-config-path",
                "/effective_config/path",
                serde_json::json!("other-effective-config.json"),
            ),
            (
                "zero-effective-config-size",
                "/effective_config/size_bytes",
                serde_json::json!(999_999_999_u64),
            ),
            (
                "zero-event-manifest-id",
                "/events/0/manifest/id",
                serde_json::json!("manifest-event-other"),
            ),
            (
                "zero-event-manifest-kind",
                "/events/0/manifest/kind",
                serde_json::json!("other"),
            ),
            (
                "zero-event-manifest-path",
                "/events/0/manifest/path",
                serde_json::json!("events/other/manifest.json"),
            ),
            (
                "zero-event-manifest-size",
                "/events/0/manifest/size_bytes",
                serde_json::json!(999_999_999_u64),
            ),
        ] {
            let mut document = zero_complete.clone();
            *document
                .pointer_mut(pointer)
                .unwrap_or_else(|| panic!("missing mutation pointer {pointer}")) = value;
            mutations.push((&zero_root, name.to_owned(), document));
        }
        let mut event_digest = zero_complete;
        event_digest["events"][0]["manifest"]["checksum"]["digest"] =
            serde_json::json!("0".repeat(64));
        mutations.push((
            &zero_root,
            "zero-complete-event-digest".to_owned(),
            event_digest,
        ));

        for (root, name, document) in mutations {
            assert_global_manifest_mutation(root, &name, &document);
            assert!(executed_cases.insert(name));
        }
        assert_eq!(
            executed_cases,
            GLOBAL_LINK_CORPUS_CASES
                .iter()
                .map(|case| (*case).to_owned())
                .collect(),
            "the recorded global-link corpus must equal the cases actually executed"
        );
        for root in [zero_root, one_root, two_root] {
            fs::remove_dir_all(root).expect("global-link corpus cleanup");
        }
    }

    #[test]
    fn draft202012_runtime_matches_python_for_complete_failed_and_contextual_corpus() {
        let (root, complete, failed) = complete_and_failed_manifest_corpus(1, "draft-corpus");
        let mut complete_failed_fasta = complete.clone();
        complete_failed_fasta["reference_reductions"][0]["derived_fasta"]["state"] =
            serde_json::json!("failed");
        complete_failed_fasta["reference_reductions"][0]["derived_fasta"]["validated"] =
            serde_json::json!(false);
        let mut complete_failed_map = complete.clone();
        complete_failed_map["reference_reductions"][0]["map"]["state"] =
            serde_json::json!("failed");
        complete_failed_map["reference_reductions"][0]["map"]["validated"] =
            serde_json::json!(false);
        let mut complete_failed_unit = complete.clone();
        complete_failed_unit["units"][1]["state"] = serde_json::json!("failed");
        complete_failed_unit["units"][1]["exportable"] = serde_json::json!(false);
        complete_failed_unit["units"][1]["diagnostic_code"] = serde_json::json!("GPS-DOMAIN-003");
        let mut failed_without_diagnostic = failed.clone();
        failed_without_diagnostic["diagnostics"] = serde_json::json!([]);
        let mut failed_complete_execution = failed.clone();
        failed_complete_execution["units"][0]["state"] = serde_json::json!("complete");
        failed_complete_execution["units"][0]["exportable"] = serde_json::json!(true);
        failed_complete_execution["units"][0]["diagnostic_code"] = serde_json::Value::Null;
        let mut failed_exportable_reduction = failed.clone();
        failed_exportable_reduction["units"][1]["exportable"] = serde_json::json!(true);

        for (name, document, expected) in [
            ("positive-complete", complete, true),
            ("positive-failed", failed, true),
            ("complete-failed-fasta", complete_failed_fasta, false),
            ("complete-failed-map", complete_failed_map, false),
            ("complete-failed-unit", complete_failed_unit, false),
            (
                "failed-without-diagnostic",
                failed_without_diagnostic,
                false,
            ),
            (
                "failed-complete-execution",
                failed_complete_execution,
                false,
            ),
            (
                "failed-exportable-reduction",
                failed_exportable_reduction,
                false,
            ),
        ] {
            let path = write_manifest_case(&root, name, &document);
            let bytes = fs::read(&path).expect("corpus case bytes");
            let python = python_schema_accepts(&path, "manifest-general-v1.schema.json");
            let rust = validate_general_manifest_bytes(&path, &bytes).is_ok();
            assert_eq!(python, expected, "Python result for {name}");
            assert_eq!(rust, expected, "Rust result for {name}");
            assert_eq!(rust, python, "Draft 2020-12 disagreement for {name}");
        }
        fs::remove_dir_all(root).expect("Draft corpus cleanup");
    }

    #[test]
    fn semantic_manifest_gate_rejects_cross_collection_contradictions_from_reloaded_bytes() {
        let (root, complete, failed) = complete_and_failed_manifest_corpus(1, "semantic-corpus");
        for (name, positive) in [
            ("positive-complete", &complete),
            ("positive-failed", &failed),
        ] {
            let path = write_manifest_case(&root, name, positive);
            let bytes = fs::read(&path).expect("positive semantic corpus bytes");
            independent_semantic_manifest_oracle(&path, &bytes)
                .unwrap_or_else(|error| panic!("independent oracle rejected {name}: {error}"));
            validate_general_manifest_bytes(&path, &bytes)
                .unwrap_or_else(|error| panic!("Rust rejected {name}: {error:?}"));
        }
        let mut wrong_fasta_digest = complete.clone();
        wrong_fasta_digest["reference_reductions"][0]["derived_fasta"]["checksum"]["digest"] =
            serde_json::json!("0".repeat(64));
        let mut wrong_map_path = complete.clone();
        wrong_map_path["reference_reductions"][0]["map"]["path"] =
            serde_json::json!("resources/references/other.map.json");
        let mut wrong_source = complete.clone();
        wrong_source["reference_reductions"][0]["source_reference_id"] =
            serde_json::json!("reference-missing");
        let mut wrong_origin = complete.clone();
        wrong_origin["references"][1]["origin"]["source_reference_id"] =
            serde_json::json!("reference-small");
        let mut orphan_reduced_reference = complete.clone();
        orphan_reduced_reference["references"][1]["id"] = serde_json::json!("reference-orphan");
        let mut reused_fasta = complete.clone();
        reused_fasta["references"][1]["fasta"] = reused_fasta["references"][0]["fasta"].clone();
        let mut incompatible_origin_kind = complete.clone();
        incompatible_origin_kind["references"][1]["origin"] =
            serde_json::json!({"kind": "generated"});
        let mut crossed_derived_fasta_id = complete.clone();
        crossed_derived_fasta_id["reference_reductions"][0]["derived_fasta"]["id"] =
            serde_json::json!("reference-fasta-reference-main");
        let mut crossed_map_id = complete.clone();
        crossed_map_id["reference_reductions"][0]["map"]["id"] =
            serde_json::json!("reference-map-reference-other");
        let mut crossed_derived_reference_id = complete;
        crossed_derived_reference_id["reference_reductions"][0]["derived_reference_id"] =
            serde_json::json!("reference-main");

        for (name, document) in [
            ("wrong-fasta-digest", wrong_fasta_digest),
            ("wrong-map-path", wrong_map_path),
            ("wrong-source", wrong_source),
            ("wrong-derived-origin", wrong_origin),
            ("orphan-reduced-reference", orphan_reduced_reference),
            ("reused-reference-fasta", reused_fasta),
            ("incompatible-origin-kind", incompatible_origin_kind),
            ("crossed-derived-fasta-id", crossed_derived_fasta_id),
            ("crossed-map-id", crossed_map_id),
            ("crossed-derived-reference-id", crossed_derived_reference_id),
        ] {
            let path = write_manifest_case(&root, name, &document);
            let bytes = fs::read(&path).expect("semantic corpus bytes");
            assert!(
                python_schema_accepts(&path, "manifest-general-v1.schema.json"),
                "JSON Schema should accept the structurally valid contradiction {name}"
            );
            assert!(
                independent_semantic_manifest_oracle(&path, &bytes).is_err(),
                "independent semantic oracle accepted contradiction {name}"
            );
            let error = validate_general_manifest_bytes(&path, &bytes)
                .expect_err("semantic Rust gate must reject cross-collection contradiction");
            assert_eq!(error.code(), "GPS-FMT-006", "case={name}");
        }
        fs::remove_dir_all(root).expect("semantic corpus cleanup");
    }
}
