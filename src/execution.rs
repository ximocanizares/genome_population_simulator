//! Recorrido de ejecución F1 sin escenarios ni unidades poblacionales sintéticas.

use std::fs;
use std::path::{Path, PathBuf};

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::application::{EffectiveConfig, EffectiveReferenceSource};
use crate::configuration::ExpectedConfigKind;
use crate::diagnostics::Diagnostic;
use crate::identity::{ReferenceId, derive_execution_id};
use crate::integrity::configuration_digest;
use crate::population_snvs::generate_population_snvs_from_reference;
use crate::publication::ExecutionWorkspace;
use crate::reference::Reference;
use crate::reference_fasta::{
    SequenceExpectation, ValidatedSuppliedReference, validate_supplied_fasta,
    validate_supplied_fasta_with_origin, write_generated_fasta,
};
use crate::reference_generation::{BaseFrequency, GeneratedSequenceSpec, generate_reference};
use crate::reference_reduction::{
    EligibleSequence, MaterializedPlan, REDUCED_REFERENCE_SELECTION_ALGORITHM, ReductionError,
    materialize_selection, preflight_reduced_reference, select_reduced_reference,
};
use crate::reference_reduction_io::{
    REDUCTION_MAP_DOCUMENT_KIND, REDUCTION_MAP_SCHEMA_VERSION, ReductionMapParameters,
    build_reduced_reference, reduction_map_document, validate_reduced_reference,
    validate_reduction_map_document,
};
use crate::results::{
    ArtifactRecord, ConfigurationInputRecord, GeneralManifest, ManifestReductionLink,
    ManifestReferenceLink, ReducedReferenceOutput, ReferenceOutput, UnitState,
    effective_config_document_with_reductions, general_manifest_document, json_document_bytes,
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
    execute_validated_config_with_effective_document(
        effective,
        destination,
        |_| Ok(()),
        std::convert::identity,
        |_| Ok(()),
    )
}

fn execute_validated_config_with_effective_document(
    effective: &EffectiveConfig,
    destination: &Path,
    after_prepare_reference: impl FnOnce(&mut PreparedReference) -> Result<(), Diagnostic>,
    transform_effective: impl FnOnce(serde_json::Value) -> serde_json::Value,
    before_manifest: impl FnOnce(&mut ExecutionWorkspace) -> Result<(), Diagnostic>,
) -> Result<ExecutionOutcome, Diagnostic> {
    let mut prepared_reference = prepare_reference(effective)?;
    after_prepare_reference(&mut prepared_reference)?;
    // F2.2.5 materializa únicamente la vista en memoria; F2.2.7 será la dueña
    // de escribir/publicar el par canónico.
    if let Some(generation) = effective.population_snvs_generation() {
        let _matrix = generate_population_snvs_from_reference(
            generation,
            prepared_reference.reference(),
            effective.reproducibility().root_seed(),
        )
        .map_err(|error| {
            Diagnostic::configuration(
                "GPS-CONF-008",
                "La generación de la matriz SNV no cumple su configuración efectiva.",
            )
            .with("component", "population_snvs.generate")
            .with("detail", error.to_string())
        })?;
    }
    // Todo error de dominio o preflight ocurre antes de crear el provisional y
    // antes de que exista cualquier RNG local de selección.
    preflight_reductions(effective, prepared_reference.reference())?;
    let configuration_digest = configuration_digest(effective.configuration_input().bytes())?;
    let execution_id = derive_execution_id(
        effective.reproducibility().root_seed(),
        &configuration_digest,
    )?;

    let mut workspace = ExecutionWorkspace::create(destination)?;
    let mut configuration_input = None;
    let mut reference_artifact = None;
    let mut manifest_references = Vec::new();
    let mut manifest_reductions = Vec::new();
    let mut published_reductions = Vec::with_capacity(effective.reduced_references().len());
    let started = (|| -> Result<(), Diagnostic> {
        configuration_input = Some(copy_configuration_input(effective, &mut workspace)?);
        reference_artifact = Some(write_reference(
            effective,
            &prepared_reference,
            &mut workspace,
        )?);
        let source_artifact = reference_artifact
            .as_ref()
            .ok_or_else(|| internal("execute", "missing source artifact"))?;
        manifest_references.push(ManifestReferenceLink {
            id: effective.reference_id().to_owned(),
            fasta_artifact_id: source_artifact.id.clone(),
        });
        prepared_reference = load_authoritative_reference(
            effective,
            &workspace,
            source_artifact,
            prepared_reference,
        )?;
        // A partir de aquí la copia recargada es la única autoridad para
        // seleccionar intervalos y materializar bases derivadas.
        let prepared_reductions = prepare_reductions(effective, prepared_reference.reference())?;
        manifest_reductions = manifest_reduction_links(effective, &prepared_reductions);
        for reduction in &mut manifest_reductions {
            reduction.source_fasta_artifact_id = Some(source_artifact.id.clone());
        }
        let reference_output = ReferenceOutput {
            reference: prepared_reference.reference(),
            fasta: source_artifact,
        };
        for (index, prepared) in prepared_reductions.into_iter().enumerate() {
            let reduction = &effective.reduced_references()[index];
            let attempt = (|| -> Result<PublishedReduction, Diagnostic> {
                let fasta_path =
                    PathBuf::from(format!("resources/references/{}.fa", reduction.id()));
                let fasta = workspace.write_streamed_artifact(
                    &format!("reference-fasta-{}", reduction.id()),
                    "reference_fasta",
                    &fasta_path,
                    |file| write_generated_fasta(&prepared.reference, file),
                )?;
                manifest_reductions[index].fasta_artifact_id = Some(fasta.id.clone());
                manifest_references.push(ManifestReferenceLink {
                    id: reduction.id().to_owned(),
                    fasta_artifact_id: fasta.id.clone(),
                });
                let authoritative_derived = load_authoritative_reduced_reference(
                    &workspace,
                    &fasta,
                    &prepared.reference,
                    reference_output.reference.id(),
                )?;
                let map_document = reduction_map_document(&prepared.plan);
                let map_bytes = json_document_bytes(&map_document)?;
                let map_path =
                    PathBuf::from(format!("resources/references/{}.map.json", reduction.id()));
                let map = workspace.write_artifact(
                    &format!("reference-map-{}", reduction.id()),
                    "reference_reduction_map",
                    &map_path,
                    &map_bytes,
                )?;
                manifest_reductions[index].map_artifact_id = Some(map.id.clone());
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
                    &authoritative_derived,
                    &reparsed,
                    ReductionMapParameters {
                        region_size_bases: reduction.region_size_bases(),
                        target_bases: reduction.target_bases(),
                    },
                )?;
                Ok(PublishedReduction {
                    effective_index: index,
                    reference: authoritative_derived,
                    plan: prepared.plan,
                    fasta,
                    map,
                })
            })();
            match attempt {
                Ok(published) => {
                    manifest_reductions[index].state = UnitState::Complete;
                    published_reductions.push(published);
                }
                Err(error) => {
                    manifest_reductions[index].state = UnitState::Failed;
                    return Err(error);
                }
            }
        }
        let reduction_outputs =
            reduction_outputs(effective, &published_reductions, source_artifact);
        let effective_document = transform_effective(effective_config_document_with_reductions(
            effective,
            reference_output,
            &reduction_outputs,
        )?);
        let effective_record = workspace.write_artifact(
            "effective-config",
            "effective_config",
            Path::new("effective-config.json"),
            &json_document_bytes(&effective_document)?,
        )?;
        before_manifest(&mut workspace)?;

        let artifacts = workspace.artifacts().to_vec();
        let manifest = GeneralManifest {
            execution_id: execution_id.as_str(),
            state: UnitState::Complete,
            configuration_input: configuration_input.as_ref(),
            effective_config: Some(&effective_record),
            references: &manifest_references,
            reduced_references: &manifest_reductions,
            artifacts: &artifacts,
            diagnostics: &[],
        };
        let manifest_bytes = json_document_bytes(&general_manifest_document(&manifest))?;
        if let Err(error) = workspace.write_general_manifest(&manifest_bytes) {
            workspace.discard_last_artifact(&effective_record)?;
            return Err(error);
        }
        Ok(())
    })();

    match started {
        Ok(()) => {
            let destination = workspace.publish()?;
            Ok(ExecutionOutcome::Published {
                destination,
                state: UnitState::Complete,
            })
        }
        Err(cause) => {
            for reduction in &mut manifest_reductions {
                reduction.state = UnitState::Failed;
                reduction.fasta_artifact_id = None;
                reduction.map_artifact_id = None;
                reduction.source_fasta_artifact_id = None;
            }
            publish_controlled_failure(
                &mut workspace,
                execution_id.as_str(),
                effective,
                &manifest_reductions,
                &cause,
            )?;
            workspace.publish()?;
            Err(cause)
        }
    }
}

fn publish_controlled_failure(
    workspace: &mut ExecutionWorkspace,
    execution_id: &str,
    effective: &EffectiveConfig,
    reductions: &[ManifestReductionLink],
    cause: &Diagnostic,
) -> Result<(), Diagnostic> {
    let configuration_artifact = workspace
        .artifacts()
        .iter()
        .find(|artifact| artifact.id == "config-general")
        .cloned()
        .filter(|artifact| valid_failure_configuration_copy(workspace, effective, artifact));
    workspace.retain_only_artifact(
        configuration_artifact
            .as_ref()
            .map(|artifact| artifact.id.as_str()),
    )?;
    let artifacts = workspace.artifacts().to_vec();
    let configuration_input = configuration_artifact.map(|copy| ConfigurationInputRecord {
        kind: "general".to_owned(),
        original_path: effective.configuration_input().original_path().to_owned(),
        copy,
    });
    let manifest = GeneralManifest {
        execution_id,
        state: UnitState::Failed,
        configuration_input: configuration_input.as_ref(),
        effective_config: None,
        references: &[],
        reduced_references: reductions,
        artifacts: &artifacts,
        diagnostics: std::slice::from_ref(cause),
    };
    workspace.write_general_manifest(&json_document_bytes(&general_manifest_document(&manifest))?)
}

fn valid_failure_configuration_copy(
    workspace: &ExecutionWorkspace,
    effective: &EffectiveConfig,
    artifact: &ArtifactRecord,
) -> bool {
    let Ok(bytes) = fs::read(workspace.root().join(&artifact.path)) else {
        return false;
    };
    if bytes != effective.configuration_input().bytes() {
        return false;
    }
    admit_toml_header(
        Path::new(effective.configuration_input().original_path()),
        bytes,
        ExpectedConfigKind::General,
    )
    .and_then(parse_declared_config)
    .is_ok()
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

fn manifest_reduction_links(
    effective: &EffectiveConfig,
    prepared: &[PreparedReduction],
) -> Vec<ManifestReductionLink> {
    effective
        .reduced_references()
        .iter()
        .zip(prepared)
        .map(|(reduction, prepared)| ManifestReductionLink {
            id: reduction.id().to_owned(),
            source_reference_id: reduction.source_reference_id().to_owned(),
            eligible_bases: effective.total_genome_size(),
            region_size_bases: reduction.region_size_bases(),
            target_mode: reduction.target().mode().to_owned(),
            target_declared: reduction.target().declared_value().to_owned(),
            target_bases: reduction.target_bases(),
            root_seed: effective.reproducibility().root_seed().to_hex(),
            selection_algorithm: REDUCED_REFERENCE_SELECTION_ALGORITHM.to_owned(),
            selection_subseed: reduction.selection_subseed().to_hex(),
            selection_order: prepared
                .plan
                .segments()
                .iter()
                .map(|segment| segment.selection_order())
                .collect(),
            state: UnitState::Incomplete,
            fasta_artifact_id: None,
            map_artifact_id: None,
            source_fasta_artifact_id: None,
            map_document_kind: REDUCTION_MAP_DOCUMENT_KIND.to_owned(),
            map_schema_version: REDUCTION_MAP_SCHEMA_VERSION.to_owned(),
        })
        .collect()
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
    let eligible = eligible_sequences(source, effective)?;
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

fn preflight_reductions(effective: &EffectiveConfig, source: &Reference) -> Result<(), Diagnostic> {
    let eligible = eligible_sequences(source, effective)?;
    for reduction in effective.reduced_references() {
        preflight_reduced_reference(
            &eligible,
            reduction.target_bases(),
            reduction.region_size_bases(),
        )
        .map_err(|error| reduction_error(reduction.id(), effective, error))?;
    }
    Ok(())
}

fn eligible_sequences(
    source: &Reference,
    effective: &EffectiveConfig,
) -> Result<Vec<EligibleSequence>, Diagnostic> {
    source
        .sequences()
        .iter()
        .map(|sequence| {
            EligibleSequence::new(sequence.id().clone(), sequence.class(), sequence.length())
                .map_err(|error| reduction_error("source", effective, error))
        })
        .collect()
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
        (EffectiveReferenceSource::Simulated, PreparedReference::Generated(_)) => {
            let physical_path = workspace.root().join(&artifact.path);
            let parsed = validate_supplied_fasta(
                &physical_path,
                effective.reference_id(),
                &sequence_expectations(effective),
            )?;
            if parsed.source_identity().size_bytes() != artifact.size_bytes
                || parsed.source_identity().digest() != artifact.digest
            {
                return Err(internal(
                    "load_authoritative_reference",
                    "cataloged generated reference changed",
                ));
            }
            let reference = Reference::new(
                ReferenceId::parse(effective.reference_id(), "execution")?,
                crate::reference::ReferenceOrigin::Generated,
                parsed.reference().sequences().to_vec(),
            )
            .map_err(|error| internal("reload_generated_reference", error))?;
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

fn load_authoritative_reduced_reference(
    workspace: &ExecutionWorkspace,
    artifact: &ArtifactRecord,
    expected: &Reference,
    source_reference_id: &ReferenceId,
) -> Result<Reference, Diagnostic> {
    let expectations = expected
        .sequences()
        .iter()
        .map(|sequence| SequenceExpectation {
            id: sequence.id().as_str(),
            class: sequence.class(),
            length: sequence.length(),
        })
        .collect::<Vec<_>>();
    let physical_path = workspace.root().join(&artifact.path);
    let parsed = validate_supplied_fasta(&physical_path, expected.id().as_str(), &expectations)?;
    if parsed.source_identity().size_bytes() != artifact.size_bytes
        || parsed.source_identity().digest() != artifact.digest
    {
        return Err(internal(
            "load_authoritative_reduced_reference",
            "cataloged derived reference changed",
        ));
    }
    let mut authoritative = Reference::new(
        expected.id().clone(),
        crate::reference::ReferenceOrigin::Reduced {
            source_reference_id: source_reference_id.clone(),
        },
        parsed.reference().sequences().to_vec(),
    )
    .map_err(|error| internal("reload_reduced_reference", error))?;
    authoritative
        .bind_authoritative_resource(PathBuf::from(&artifact.path))
        .map_err(|error| internal("bind_reduced_authoritative_resource", error))?;
    Ok(authoritative)
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

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use serde_json::json;

    use super::*;
    use crate::application::load_and_validate_config;

    static COUNTER: AtomicU64 = AtomicU64::new(0);

    #[test]
    fn invalid_effective_document_is_removed_before_failed_fallback_publication() {
        let root = std::env::temp_dir().join(format!(
            "gps-f14-invalid-effective-fallback-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("test root");
        let destination = root.join("result");
        let config =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/config/general-seeded.toml");
        let effective = load_and_validate_config(&config).expect("valid input config");

        let error = execute_validated_config_with_effective_document(
            &effective,
            &destination,
            |_| Ok(()),
            |mut document| {
                document["effective_config_kind"] = json!("unknown");
                document
            },
            |_| Ok(()),
        )
        .expect_err("invalid produced effective config");

        assert_eq!(error.code(), "GPS-FMT-006");
        assert_eq!(error.category(), "FMT");
        assert_eq!(error.exit_code(), 1);
        let manifest: serde_json::Value = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("failed manifest"),
        )
        .expect("manifest JSON");
        assert_eq!(manifest["state"], "failed");
        assert_eq!(
            manifest["effective_config_artifact_id"],
            serde_json::Value::Null
        );
        assert_eq!(manifest["diagnostics"][0]["code"], "GPS-FMT-006");
        assert!(!destination.join("effective-config.json").exists());
        assert!(
            !destination
                .join("resources/references/ref-main.fa")
                .exists()
        );
        assert_eq!(manifest["references"], json!([]));
        assert_eq!(manifest["reference_reductions"], json!([]));
        assert_eq!(
            manifest["output_categories"]["minimum"],
            json!(["config-general"])
        );
        assert_eq!(manifest["output_categories"]["requested"], json!([]));
        assert!(
            manifest["artifacts"]
                .as_array()
                .expect("artifact catalog")
                .iter()
                .all(|artifact| artifact["id"] == "config-general")
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn physical_resource_failure_marks_reduction_failed_without_partial_artifacts() {
        let root = std::env::temp_dir().join(format!(
            "gps-f14-physical-authority-fallback-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("test root");
        let destination = root.join("result");
        let config = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/reference_reduction/general.toml");
        let effective = load_and_validate_config(&config).expect("valid input config");

        let error = execute_validated_config_with_effective_document(
            &effective,
            &destination,
            |_| Ok(()),
            std::convert::identity,
            |workspace| {
                let artifact = workspace
                    .artifacts()
                    .iter()
                    .find(|artifact| artifact.id == "reference-fasta-reference-small")
                    .expect("derived reference artifact");
                let mut bytes =
                    fs::read(workspace.root().join(&artifact.path)).expect("derived FASTA");
                let sequence_start = bytes
                    .iter()
                    .position(|byte| *byte == b'\n')
                    .expect("FASTA header")
                    + 1;
                let base = bytes[sequence_start];
                bytes[sequence_start] = if base == b'A' { b'C' } else { b'A' };
                workspace.replace_artifact_bytes_for_validation_test(
                    "reference-fasta-reference-small",
                    &bytes,
                )
            },
        )
        .expect_err("physical bytes disagree with generated authority");

        assert_eq!(error.code(), "GPS-FMT-006");
        let manifest: serde_json::Value = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("failed manifest"),
        )
        .expect("manifest JSON");
        assert_eq!(manifest["state"], "failed");
        assert_eq!(
            manifest["effective_config_artifact_id"],
            serde_json::Value::Null
        );
        assert_eq!(manifest["diagnostics"][0]["code"], "GPS-FMT-006");
        assert!(!destination.join("effective-config.json").exists());
        for relative in [
            "resources/references/reference-main.fa",
            "resources/references/reference-small.fa",
            "resources/references/reference-small.map.json",
        ] {
            assert!(!destination.join(relative).exists(), "retained {relative}");
        }
        assert_eq!(manifest["references"], json!([]));
        assert_eq!(
            manifest["reference_reductions"].as_array().unwrap().len(),
            1
        );
        let reduction = &manifest["reference_reductions"][0];
        assert_eq!(reduction["id"], "reference-small");
        assert_eq!(reduction["state"], "failed");
        assert_eq!(reduction["fasta_artifact_id"], serde_json::Value::Null);
        assert_eq!(reduction["map_artifact_id"], serde_json::Value::Null);
        assert_eq!(
            reduction["source_fasta_artifact_id"],
            serde_json::Value::Null
        );
        assert_eq!(
            manifest["output_categories"],
            json!({
                "minimum": ["config-general"],
                "defaulted": [],
                "requested": [],
                "disabled": []
            })
        );
        assert!(
            manifest["artifacts"]
                .as_array()
                .expect("artifact catalog")
                .iter()
                .all(|artifact| artifact["id"] == "config-general")
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn reduction_selection_uses_reloaded_copy_not_the_stale_origin_object() {
        use crate::reference::{DnaSequence, ReferenceSequence};

        let root = std::env::temp_dir().join(format!(
            "gps-f16-authoritative-copy-selection-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("test root");
        let destination = root.join("result");
        let config = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/reference_reduction/general.toml");
        let effective = load_and_validate_config(&config).expect("valid input config");

        execute_validated_config_with_effective_document(
            &effective,
            &destination,
            |prepared| {
                let PreparedReference::Supplied(supplied) = prepared else {
                    return Err(internal("test_stale_origin", "expected supplied reference"));
                };
                let reference = supplied.reference();
                let mut sequences = reference.sequences().to_vec();
                let first = &sequences[0];
                let mut bases = first.sequence().as_bytes().to_vec();
                for base in &mut bases {
                    *base = if *base == b'N' { b'N' } else { b'A' };
                }
                let dna = DnaSequence::try_from_bytes(&bases)
                    .map_err(|error| internal("test_stale_origin_dna", error))?;
                sequences[0] = ReferenceSequence::new(
                    first.id().clone(),
                    first.header().clone(),
                    first.class(),
                    &dna,
                )
                .map_err(|error| internal("test_stale_origin_sequence", error))?;
                let stale = Reference::new(
                    reference.id().clone(),
                    reference.origin().clone(),
                    sequences,
                )
                .map_err(|error| internal("test_stale_origin_reference", error))?;
                *supplied.reference_mut() = stale;
                Ok(())
            },
            std::convert::identity,
            |_| Ok(()),
        )
        .expect("selection must ignore stale origin object after the authoritative copy exists");

        let copied = fs::read(destination.join("resources/references/reference-main.fa"))
            .expect("authoritative source copy");
        assert_eq!(
            copied,
            fs::read(config.parent().unwrap().join("source.fa")).unwrap()
        );
        let effective_document: serde_json::Value = serde_json::from_slice(
            &fs::read(destination.join("effective-config.json")).expect("effective config"),
        )
        .expect("effective JSON");
        assert_eq!(
            effective_document["reduced_references"][0]["source_fasta_artifact_id"],
            "reference-fasta"
        );
        let baseline = root.join("baseline");
        execute_validated_config(&effective, &baseline).expect("unpoisoned baseline");
        for relative in [
            "resources/references/reference-small.fa",
            "resources/references/reference-small.map.json",
        ] {
            assert_eq!(
                fs::read(destination.join(relative)).expect("stale-origin run artifact"),
                fs::read(baseline.join(relative)).expect("baseline artifact"),
                "{relative} must depend only on the authoritative copy"
            );
        }
        fs::remove_dir_all(root).expect("cleanup");
    }
}
