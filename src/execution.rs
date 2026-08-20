//! Recorrido de ejecución F1 sin escenarios ni unidades poblacionales sintéticas.

use std::fs;
use std::path::{Path, PathBuf};

use crate::application::{EffectiveConfig, EffectiveReferenceSource};
use crate::diagnostics::Diagnostic;
use crate::identity::{AutomaticIdRequest, IdClass, IdentityRegistry, ReferenceId};
use crate::integrity::Sha256Writer;
use crate::publication::{
    ArtifactStage, ExecutionWorkspace, LinuxAtomicPublisher, PublishableUnit,
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
    ArtifactRecord, ConfigurationInputRecord, GeneralManifest, ReducedReferenceOutput,
    ReferenceOutput, UnitRecord, UnitState, effective_config_document_with_reductions,
    general_manifest_document, json_document_bytes,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ExecutionOutcome {
    Published {
        destination: PathBuf,
        state: UnitState,
    },
}

pub fn execute_validated_config(
    effective: &EffectiveConfig,
    destination: &Path,
) -> Result<ExecutionOutcome, Diagnostic> {
    let mut prepared_reference = prepare_reference(effective)?;
    let configuration_digest = configuration_digest(effective.configuration_input().bytes())?;
    let request = AutomaticIdRequest::new(
        IdClass::Execution,
        "execution",
        &[("config_sha256", configuration_digest.as_str())],
    )?;
    let mut identities = IdentityRegistry::default();
    let execution_id =
        identities.derive_execution(effective.reproducibility().root_seed(), &request)?;

    let mut workspace = ExecutionWorkspace::create(destination)?;
    let configuration_input = copy_configuration_input(effective, &mut workspace)?;
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
    for (index, mut prepared) in prepared_reductions.into_iter().enumerate() {
        let reduction = &effective.reduced_references()[index];
        let fasta_path = PathBuf::from(format!("resources/references/{}.fa", reduction.id()));
        let fasta = workspace.write_streamed_artifact(
            &format!("reference-fasta-{}", reduction.id()),
            "reference_fasta",
            &fasta_path,
            |file| write_generated_fasta(&prepared.reference, file),
        )?;
        prepared
            .reference
            .bind_authoritative_resource(PathBuf::from(&fasta.path))
            .map_err(|error| internal("bind_reduced_authoritative_resource", error))?;
        let map_document = reduction_map_document(&prepared.plan);
        let map_bytes = json_document_bytes(&map_document)?;
        let map_path = PathBuf::from(format!("resources/references/{}.map.json", reduction.id()));
        let map = workspace.write_artifact(
            &format!("reference-map-{}", reduction.id()),
            "reference_reduction_map",
            &map_path,
            &map_bytes,
        )?;
        let reparsed: serde_json::Value = serde_json::from_slice(
            &fs::read(workspace.root().join(&map.path)).map_err(|error| {
                io_error(
                    &workspace.root().join(&map.path),
                    "read_reduction_map",
                    error,
                )
            })?,
        )
        .map_err(|error| {
            Diagnostic::format(
                "GPS-FMT-007",
                "El mapa de reducción no permite reconstruir exactamente la referencia derivada.",
            )
            .with("path", &map.path)
            .with("detail", error.to_string())
        })?;
        validate_reduction_map_document(
            reference_output.reference,
            &prepared.reference,
            &reparsed,
            ReductionMapParameters {
                region_size_bases: reduction.region_size_bases(),
                target_bases: reduction.target_bases(),
            },
        )?;
        published_reductions.push(PublishedReduction {
            effective_index: index,
            reference: prepared.reference,
            plan: prepared.plan,
            fasta,
            map,
        });
    }
    let reduction_outputs =
        reduction_outputs(effective, &published_reductions, &reference_artifact);
    let effective_document =
        effective_config_document_with_reductions(effective, reference_output, &reduction_outputs)?;
    let effective_record = workspace.write_artifact(
        "effective-config",
        "effective_config",
        Path::new("effective-config.json"),
        &json_document_bytes(&effective_document)?,
    )?;

    let artifact_ids = std::iter::once(configuration_input.copy.id.clone())
        .chain(std::iter::once(reference_artifact.id.clone()))
        .chain(std::iter::once(effective_record.id.clone()))
        .chain(
            published_reductions
                .iter()
                .flat_map(|reduction| [reduction.fasta.id.clone(), reduction.map.id.clone()]),
        )
        .collect::<Vec<_>>();
    let mut execution_unit =
        PublishableUnit::new(execution_id.as_str(), "execution", artifact_ids.clone());
    for id in &artifact_ids {
        execution_unit.record_artifact(id, ArtifactStage::Cataloged)?;
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
        .collect::<Vec<_>>();
    let artifacts = workspace.artifacts().to_vec();
    let manifest = GeneralManifest {
        execution_id: execution_id.as_str(),
        state: UnitState::Complete,
        effective,
        configuration_input: &configuration_input,
        effective_config: Some(&effective_record),
        reference: Some(reference_output),
        reduced_references: &reduction_outputs,
        artifacts: &artifacts,
        units: &units,
        diagnostics: &[],
    };
    let manifest_bytes = json_document_bytes(&general_manifest_document(&manifest))?;
    let validated =
        workspace.validate_general_manifest(Path::new("manifest.json"), &manifest_bytes)?;
    workspace.write_general_manifest(&validated)?;
    let destination = workspace.publish_with(&LinuxAtomicPublisher)?;
    Ok(ExecutionOutcome::Published {
        destination,
        state: UnitState::Complete,
    })
}

enum PreparedReference {
    Supplied(ValidatedSuppliedReference),
    Generated(Reference),
}

impl PreparedReference {
    fn reference(&self) -> &Reference {
        match self {
            Self::Supplied(value) => value.reference(),
            Self::Generated(value) => value,
        }
    }

    fn reference_mut(&mut self) -> &mut Reference {
        match self {
            Self::Supplied(value) => value.reference_mut(),
            Self::Generated(value) => value,
        }
    }
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

fn prepare_reference(effective: &EffectiveConfig) -> Result<PreparedReference, Diagnostic> {
    match effective.reference_source() {
        EffectiveReferenceSource::Supplied { original_path } => validate_supplied_fasta(
            original_path,
            effective.reference_id(),
            &sequence_expectations(effective),
        )
        .map(PreparedReference::Supplied),
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
            let values = effective
                .composition()
                .ok_or_else(|| internal("prepare_reference", "missing composition"))?
                .values();
            let composition = (*b"ACGT")
                .into_iter()
                .zip(values)
                .map(|(base, probability)| BaseFrequency {
                    base,
                    probability: *probability,
                })
                .collect::<Vec<_>>();
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
    let eligible = source
        .sequences()
        .iter()
        .map(|sequence| {
            EligibleSequence::new(sequence.id().clone(), sequence.class(), sequence.length())
                .map_err(|error| reduction_error("source", effective, error))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut prepared = Vec::with_capacity(effective.reduced_references().len());
    for reduction in effective.reduced_references() {
        preflight_reduced_reference(
            &eligible,
            reduction.target_bases(),
            reduction.region_size_bases(),
        )
        .map_err(|error| reduction_error(reduction.id(), effective, error))?;
        let mut rng = reduction.selection_subseed().local_rng();
        let selection = select_reduced_reference(
            source.id(),
            &eligible,
            reduction.target_bases(),
            reduction.region_size_bases(),
            &mut rng,
        )
        .map_err(|error| reduction_error(reduction.id(), effective, error))?;
        let derived_id = ReferenceId::parse(reduction.id(), "execution")?;
        let plan = materialize_selection(&selection, &derived_id, &eligible)
            .map_err(|error| reduction_error(reduction.id(), effective, error))?;
        let reference = build_reduced_reference(source, &plan)?;
        validate_reduced_reference(source, &reference, &plan)?;
        prepared.push(PreparedReduction { reference, plan });
    }
    Ok(prepared)
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
            let supplied = validate_supplied_fasta_with_origin(
                &workspace.root().join(&artifact.path),
                original_path,
                effective.reference_id(),
                &sequence_expectations(effective),
            )?;
            if supplied.source_identity().size_bytes() != artifact.size_bytes
                || supplied.source_identity().digest() != artifact.digest
            {
                return Err(internal(
                    "load_authoritative_reference",
                    "cataloged reference changed",
                ));
            }
            PreparedReference::Supplied(supplied)
        }
        (EffectiveReferenceSource::Simulated, PreparedReference::Generated(reference)) => {
            PreparedReference::Generated(reference)
        }
        _ => return Err(internal("load_authoritative_reference", "origin changed")),
    };
    authoritative
        .reference_mut()
        .bind_authoritative_resource(PathBuf::from(&artifact.path))
        .map_err(|error| internal("bind_authoritative_reference", error))?;
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
        _ => Err(internal("write_reference", "origin changed")),
    }
}

fn copy_configuration_input(
    effective: &EffectiveConfig,
    workspace: &mut ExecutionWorkspace,
) -> Result<ConfigurationInputRecord, Diagnostic> {
    let input = effective.configuration_input();
    let copy = workspace.write_artifact(
        "config-general",
        "config_general_toml",
        Path::new("inputs/config/general.toml"),
        input.bytes(),
    )?;
    Ok(ConfigurationInputRecord {
        kind: "general".to_owned(),
        original_path: input.original_path().to_owned(),
        copy,
    })
}

fn reduction_outputs<'a>(
    effective: &'a EffectiveConfig,
    reductions: &'a [PublishedReduction],
    source_fasta: &'a ArtifactRecord,
) -> Vec<ReducedReferenceOutput<'a>> {
    reductions
        .iter()
        .map(|reduction| ReducedReferenceOutput {
            effective: &effective.reduced_references()[reduction.effective_index],
            reference: &reduction.reference,
            fasta: &reduction.fasta,
            map: &reduction.map,
            plan: &reduction.plan,
            source_fasta,
            state: UnitState::Complete,
        })
        .collect()
}

fn configuration_digest(bytes: &[u8]) -> Result<String, Diagnostic> {
    let mut hasher = Sha256Writer::default();
    hasher.update(b"genome_population_simulator/v1/configuration-digest\0")?;
    hasher.update(
        &u64::try_from(bytes.len())
            .map_err(|_| internal("configuration_digest", "input exceeds u64"))?
            .to_be_bytes(),
    )?;
    hasher.update(bytes)?;
    Ok(hasher.finish("configuration-input").digest().to_owned())
}

fn reduction_error(
    reference: &str,
    effective: &EffectiveConfig,
    error: ReductionError,
) -> Diagnostic {
    match error {
        ReductionError::AllocationFailed => Diagnostic::resource(
            "GPS-RES-002",
            "Los recursos disponibles no permiten completar el plan sin degradar resultados obligatorios.",
        )
        .with("resource", "memory")
        .with("required", "complete immutable reduction plan")
        .with("available", "allocation rejected"),
        ReductionError::ArithmeticOverflow
        | ReductionError::TooManySegments
        | ReductionError::Coordinate(_)
        | ReductionError::SelectionDoesNotMatchCatalog => internal("plan_reduced_reference", error),
        _ => Diagnostic::configuration(
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
                .map_or_else(|| "unknown".to_owned(), |value| value.target_bases().to_string()),
        )
        .with("eligible_bases", effective.total_genome_size().to_string())
        .with("detail", format!("{error:?}")),
    }
}

fn internal(operation: &str, detail: impl std::fmt::Debug) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "execution")
    .with("operation", operation)
    .with("detail", format!("{detail:?}"))
}

fn io_error(path: &Path, operation: &str, error: std::io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}
