//! Recorrido de ejecución F1 sin escenarios ni unidades poblacionales sintéticas.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::application::{
    EffectiveConfig, EffectiveDirectEvent, EffectiveReferenceSource, EffectiveTechnicalEventTarget,
};
use crate::biological_contamination::{
    BiologicalContaminationMaterializer, CleanLibraryPairs, CleanPairKey,
    plan_biological_contamination,
};
use crate::configuration::ExpectedConfigKind;
use crate::coordinates::BasePosition;
use crate::demultiplexing::{Demultiplexer, DemultiplexingLibrary, PostContaminationPair, Udi};
use crate::diagnostics::Diagnostic;
use crate::error_quality::{
    ErrorQualityAccumulator, ErrorQualitySettings, ObservedPairTransformation, transform_paired_end,
};
use crate::event_manifest::{
    EventManifest, EventManifestArtifact, EventManifestCoverage, EventManifestDemultiplexing,
    EventManifestErrorQuality, EventManifestFinalFastqSummary, EventManifestReference,
    event_manifest_bytes, event_manifest_relative_path, validate_event_manifest_bytes,
};
use crate::fastq_output::write_final_paired_fastq;
use crate::fastq_validation::F4CoverageReport;
use crate::identity::{ReferenceId, derive_execution_id};
use crate::initial_individual_integration::IncompatibleGenotypeAuthorization;
use crate::initial_individual_materialization::{
    InitialIndividualsMaterialization, InitialIndividualsMaterializationError,
    materialize_initial_individuals,
};
use crate::integrity::{Sha256Writer, configuration_digest};
use crate::materialization::SnpLocus;
use crate::population_snvs::{
    PopulationSnvsProvenance, canonical_metadata_document, canonical_tsv_bytes,
    generate_population_snvs_from_reference_with_tolerance, normalization_audit_document,
    parse_supplied_tsv, prepare_population_snvs_from_reference,
    validate_reusable_population_snvs_with_matrix_sha256,
};
use crate::publication::ExecutionWorkspace;
use crate::reference::Reference;
use crate::reference_fasta::{
    SequenceExpectation, ValidatedSuppliedReference, validate_supplied_fasta,
    validate_supplied_fasta_with_origin, write_generated_fasta,
};
use crate::reference_generation::{
    BaseFrequency, GeneratedSequenceSpec, generate_reference_with_unknown_fraction,
};
use crate::reference_reduction::{
    EligibleSequence, MaterializedPlan, ReductionError, materialize_selection,
    preflight_reduced_reference, select_reduced_reference,
};
use crate::reference_reduction_io::{
    REDUCTION_MAP_DOCUMENT_KIND, REDUCTION_MAP_SCHEMA_VERSION, ReductionMapParameters,
    build_reduced_reference, reduction_map_document, validate_reduction_map_document,
};
use crate::results::{
    ArtifactRecord, ConfigurationInputRecord, GeneralManifest, ManifestEventLink,
    ManifestReductionLink, ManifestReferenceLink, PopulationSnvsManifestLink,
    ReducedReferenceOutput, ReferenceOutput, UnitState, VcfManifestLink,
    effective_config_document_with_population_snvs, general_manifest_document, json_document_bytes,
};
use crate::segments::{CanonicalBase, ModificationProvenance};
use crate::small_event_composition::{
    DeclaredSmallEvent, TargetedCurrentDeletion, TargetedCurrentInsertion,
};
use crate::snp_application::{ApplicableSnp, TargetedSnp};
use crate::technical_identity::derive_technical_identity_plan;
use crate::technical_profile::resolve_f4_technical_profile;
use crate::truth_projection::{TruthProjectionInput, project_truth, write_truth_artifact};
use crate::truth_schema::{
    NamedDigest, TRUTH_COORDINATE_CONVENTION, TRUTH_DOCUMENT_KIND, TRUTH_SCHEMA_VERSION,
    TruthDocument, TruthHeader, TruthSourceKind, TruthSourceRef,
};
use crate::truth_validation::reference_content_sha256;
use crate::vcf_publication::publish_truth_vcf;
use crate::wgs_sampling::WgsPairedEndSampler;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ExecutionOutcome {
    Published {
        destination: PathBuf,
        state: UnitState,
    },
}

/// Confirmaciones explícitas de la invocación, separadas del TOML autorado.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ExecutionAuthorizations {
    pub allow_empty_snv_matrix: bool,
    pub allow_incompatible_genotypes: bool,
}

pub fn execute_validated_config(
    effective: &EffectiveConfig,
    destination: &Path,
) -> Result<ExecutionOutcome, Diagnostic> {
    execute_validated_config_with_authorizations(
        effective,
        destination,
        ExecutionAuthorizations::default(),
    )
}

pub fn execute_validated_config_with_authorizations(
    effective: &EffectiveConfig,
    destination: &Path,
    authorizations: ExecutionAuthorizations,
) -> Result<ExecutionOutcome, Diagnostic> {
    execute_validated_config_with_effective_document(
        effective,
        destination,
        authorizations,
        |_| Ok(()),
        std::convert::identity,
        |_| Ok(()),
    )
}

/// Ejecuta todo el preflight de recursos y creación directa sin crear un
/// workspace ni publicar artefactos. `validate` y `run` comparten esta frontera
/// para que las autorizaciones y los rechazos F2.3.4 no diverjan.
pub fn preflight_validated_config_with_authorizations(
    effective: &EffectiveConfig,
    authorizations: ExecutionAuthorizations,
) -> Result<(), Diagnostic> {
    let prepared_reference = prepare_reference(effective)?;
    preflight_prepared_config(effective, authorizations, &prepared_reference).map(|_| ())
}

fn execute_validated_config_with_effective_document(
    effective: &EffectiveConfig,
    destination: &Path,
    authorizations: ExecutionAuthorizations,
    after_prepare_reference: impl FnOnce(&mut PreparedReference) -> Result<(), Diagnostic>,
    transform_effective: impl FnOnce(serde_json::Value) -> serde_json::Value,
    before_manifest: impl FnOnce(&mut ExecutionWorkspace) -> Result<(), Diagnostic>,
) -> Result<ExecutionOutcome, Diagnostic> {
    let mut prepared_reference = prepare_reference(effective)?;
    after_prepare_reference(&mut prepared_reference)?;
    // Todo error de dominio o preflight ocurre antes de crear el provisional y
    // antes de que exista cualquier RNG local de selección.
    let prepared_execution =
        preflight_prepared_config(effective, authorizations, &prepared_reference)?;
    let configuration_digest = configuration_digest(effective.configuration_input().bytes())?;
    let execution_id = derive_execution_id(
        effective.reproducibility().root_seed(),
        &configuration_digest,
    )?;

    // Una ejecución carga el genoma tres veces en momentos distintos: al
    // prepararlo, al releer la copia publicada como autoridad, y al validar
    // físicamente esa publicación. Cada carga retiene todas las bases, de modo
    // que dos que se solapen doblan el pico de memoria.
    //
    // Aquí se envuelve la preparada en `Option` para poder soltarla antes de
    // escribir el manifiesto, que es lo que dispara la validación física. La
    // otra liberación está en `load_authoritative_reference`. Hacen falta las
    // dos: si solo se elimina uno de los dos solapamientos, el otro sostiene el
    // mismo pico y la mejora no se observa.
    let mut prepared_reference = Some(prepared_reference);
    let mut workspace = ExecutionWorkspace::create(destination)?;
    let mut configuration_input = None;
    let mut reference_artifact = None;
    let mut manifest_references = Vec::new();
    let mut manifest_reductions = Vec::new();
    let mut published_reductions = Vec::with_capacity(effective.reduced_references().len());
    let mut population_snvs = None;
    let mut vcf_outputs = Vec::new();
    let mut event_links: Option<Vec<ManifestEventLink>> = None;
    let started = (|| -> Result<(), Diagnostic> {
        configuration_input = Some(copy_configuration_input(effective, &mut workspace)?);
        reference_artifact = Some(write_reference(
            effective,
            prepared_reference
                .as_ref()
                .ok_or_else(|| internal("execute", "reference already released"))?,
            &mut workspace,
        )?);
        let source_artifact = reference_artifact
            .as_ref()
            .ok_or_else(|| internal("execute", "missing source artifact"))?;
        manifest_references.push(ManifestReferenceLink {
            id: effective.reference_id().to_owned(),
            fasta_artifact_id: source_artifact.id.clone(),
        });
        prepared_reference = Some(load_authoritative_reference(
            effective,
            &workspace,
            source_artifact,
            prepared_reference
                .take()
                .ok_or_else(|| internal("execute", "reference already released"))?,
        )?);
        let reference = prepared_reference
            .as_ref()
            .ok_or_else(|| internal("execute", "reference already released"))?
            .reference();
        population_snvs = publish_prepared_population_snvs(
            prepared_execution.population_snvs.as_ref(),
            source_artifact,
            &mut workspace,
        )?;
        if let (Some(link), Some(initial_individuals)) = (
            &mut population_snvs,
            &prepared_execution.initial_individuals,
        ) {
            link.initial_population_statistics = Some(initial_individuals.statistics().clone());
        }
        // A partir de aquí la copia recargada es la única autoridad para
        // seleccionar intervalos y materializar bases derivadas.
        let prepared_reductions = prepare_reductions(effective, reference)?;
        manifest_reductions = manifest_reduction_links(effective, &prepared_reductions);
        for reduction in &mut manifest_reductions {
            reduction.source_fasta_artifact_id = Some(source_artifact.id.clone());
        }
        let reference_output = ReferenceOutput {
            reference,
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
        let truth = write_execution_truth(
            &mut workspace,
            execution_id.as_str(),
            reference,
            prepared_execution.population_snvs.as_ref(),
            prepared_execution.initial_individuals.as_ref(),
            source_artifact,
            configuration_input
                .as_ref()
                .ok_or_else(|| internal("execute", "missing configuration artifact"))?,
            effective,
        )?;
        if let Some(mode) = effective.vcf_output() {
            let published = publish_truth_vcf(&mut workspace, reference, &truth, mode)?;
            vcf_outputs.push(VcfManifestLink {
                state: UnitState::Complete,
                reference_id: published.reference_id,
                samples: published.samples,
                mode: published.mode,
                vcf_artifact_id: published.vcf.id,
                index_kind: published.index_kind,
                index_artifact_id: published.index.map(|artifact| artifact.id),
            });
        }
        if let Some(event) = effective.technical_event() {
            match publish_technical_event(
                &mut workspace,
                effective,
                event,
                reference,
                source_artifact,
                prepared_execution.initial_individuals.as_ref(),
            ) {
                Ok(link) => event_links = Some(vec![link]),
                Err(failure) => {
                    event_links = failure.link.map(|link| vec![link]);
                    return Err(failure.cause);
                }
            }
        }
        let effective_document =
            transform_effective(effective_config_document_with_population_snvs(
                effective,
                authorizations,
                reference_output,
                &reduction_outputs,
                population_snvs.as_ref(),
            )?);
        // Ningún artefacto posterior necesita ya las bases de la referencia. Se
        // libera antes de escribir el manifiesto, porque esa escritura dispara
        // la validación física, que relee el FASTA publicado y construye otra
        // referencia completa. Sin esto ambas conviven y el pico de memoria
        // dobla el tamaño del genoma de partida.
        prepared_reference = None;
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
            population_snvs: population_snvs.as_ref(),
            vcf_outputs: &vcf_outputs,
            events: event_links.as_deref(),
            authorizations,
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
            let failure_references = if event_links.is_some() {
                manifest_references
                    .iter()
                    .filter(|reference| {
                        !manifest_reductions
                            .iter()
                            .any(|reduction| reduction.id == reference.id)
                    })
                    .cloned()
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            };
            publish_controlled_failure(
                &mut workspace,
                execution_id.as_str(),
                effective,
                authorizations,
                &manifest_reductions,
                event_links.as_deref(),
                &failure_references,
                &cause,
            )?;
            workspace.publish()?;
            Err(cause)
        }
    }
}

/// Publica siempre la proyección obligatoria de verdad. Las ejecuciones F1 que
/// aún no materializan población conservan un documento vacío y reconstruible;
/// cuando F2.3 preparó matriz e individuos, reutiliza exclusivamente esa
/// verdad faseada para escribir todos los registros.
#[allow(clippy::too_many_arguments)]
fn write_execution_truth(
    workspace: &mut ExecutionWorkspace,
    execution_id: &str,
    reference: &Reference,
    population_snvs: Option<&PreparedPopulationSnvsForExecution>,
    initial_individuals: Option<&InitialIndividualsMaterialization>,
    reference_artifact: &ArtifactRecord,
    configuration_input: &ConfigurationInputRecord,
    effective: &EffectiveConfig,
) -> Result<TruthDocument, Diagnostic> {
    let document = match (population_snvs, initial_individuals) {
        (Some(prepared), Some(individuals)) => {
            let matrix_digest = prepared_matrix_digest(prepared, workspace)
                .ok_or_else(|| internal("truth", "missing matrix identity"))?;
            let mut source_digests = vec![named_digest(
                TruthSourceKind::Reference,
                reference.id().as_str(),
                &reference_artifact.digest,
            )];
            source_digests.push(named_digest(
                TruthSourceKind::Matrix,
                &prepared.matrix.canonical_fingerprint(),
                &matrix_digest,
            ));
            source_digests.extend(prepared.matrix.rows().iter().map(|row| {
                named_digest(
                    TruthSourceKind::Matrix,
                    matrix_source_id(row.source_id()),
                    &matrix_digest,
                )
            }));
            source_digests.extend(effective.direct_events().iter().map(|event| {
                named_digest(
                    TruthSourceKind::Explicit,
                    event.event_id().as_str(),
                    &configuration_input.copy.digest,
                )
            }));
            deduplicate_digests(&mut source_digests);
            let group_digests = prepared
                .matrix
                .population_ids()
                .iter()
                .map(|population_id| {
                    named_digest(
                        TruthSourceKind::Matrix,
                        population_id.as_str(),
                        &matrix_digest,
                    )
                })
                .collect::<Vec<_>>();
            project_truth(TruthProjectionInput {
                execution_id,
                reference,
                individuals,
                matrix: &prepared.matrix,
                source_digests: &source_digests,
                group_digests: &group_digests,
            })
            .map_err(|error| {
                Diagnostic::format("GPS-FMT-006", "No se pudo proyectar la verdad genómica.")
                    .with("detail", error.to_string())
            })?
        }
        (Some(prepared), None) => empty_truth_document(
            execution_id,
            reference,
            reference_artifact,
            Some((prepared, workspace)),
        )?,
        (None, None) => empty_truth_document(execution_id, reference, reference_artifact, None)?,
        _ => return Err(internal("truth", "matrix and initial individuals diverged")),
    };
    write_truth_artifact(workspace, &document)?;
    Ok(document)
}

fn empty_truth_document(
    execution_id: &str,
    reference: &Reference,
    reference_artifact: &ArtifactRecord,
    matrix: Option<(&PreparedPopulationSnvsForExecution, &ExecutionWorkspace)>,
) -> Result<TruthDocument, Diagnostic> {
    let mut source_digests = vec![named_digest(
        TruthSourceKind::Reference,
        reference.id().as_str(),
        &reference_artifact.digest,
    )];
    let mut group_digests = Vec::new();
    if let Some((prepared, workspace)) = matrix {
        let matrix_digest = prepared_matrix_digest(prepared, workspace)
            .ok_or_else(|| internal("truth", "missing matrix identity"))?;
        source_digests.push(named_digest(
            TruthSourceKind::Matrix,
            &prepared.matrix.canonical_fingerprint(),
            &matrix_digest,
        ));
        source_digests.extend(prepared.matrix.rows().iter().map(|row| {
            named_digest(
                TruthSourceKind::Matrix,
                matrix_source_id(row.source_id()),
                &matrix_digest,
            )
        }));
        group_digests.extend(
            prepared
                .matrix
                .population_ids()
                .iter()
                .map(|population_id| {
                    named_digest(
                        TruthSourceKind::Matrix,
                        population_id.as_str(),
                        &matrix_digest,
                    )
                }),
        );
        deduplicate_digests(&mut source_digests);
    }
    Ok(TruthDocument {
        header: TruthHeader {
            document_kind: TRUTH_DOCUMENT_KIND.to_owned(),
            schema_version: TRUTH_SCHEMA_VERSION.to_owned(),
            execution_id: execution_id.to_owned(),
            reference_id: reference.id().as_str().to_owned(),
            reference_content_sha256: reference_content_sha256(reference),
            coordinate_convention: TRUTH_COORDINATE_CONVENTION.to_owned(),
            copy_count: 0,
            genotype_count: 0,
            source_digests,
            group_digests,
        },
        records: vec![],
    })
}

fn prepared_matrix_digest(
    prepared: &PreparedPopulationSnvsForExecution,
    workspace: &ExecutionWorkspace,
) -> Option<String> {
    workspace
        .artifacts()
        .iter()
        .find(|artifact| artifact.id == "population-snvs-tsv")
        .map(|artifact| artifact.digest.clone())
        .or_else(|| match &prepared.publication {
            PopulationSnvsPublication::Reuse { matrix_sha256, .. } => Some(matrix_sha256.clone()),
            PopulationSnvsPublication::Publish { .. } => None,
        })
}

struct TechnicalEventAssemblyFailure {
    cause: Diagnostic,
    link: Option<ManifestEventLink>,
}

struct TechnicalEventArtifact {
    individual_id: crate::identity::IndividualId,
    library_id: crate::identity::LibraryId,
    mate: crate::technical_identity::TechnicalMate,
    id: String,
    state: UnitState,
    artifact: Option<ArtifactRecord>,
    diagnostic: Option<Diagnostic>,
}

impl TechnicalEventArtifact {
    fn complete(
        individual_id: crate::identity::IndividualId,
        library_id: crate::identity::LibraryId,
        mate: crate::technical_identity::TechnicalMate,
        artifact: ArtifactRecord,
    ) -> Self {
        Self {
            individual_id,
            library_id,
            mate,
            id: artifact.id.clone(),
            state: UnitState::Complete,
            artifact: Some(artifact),
            diagnostic: None,
        }
    }
}

#[cfg(test)]
thread_local! {
    static FORCE_FASTQ_REVALIDATION_FAILURE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Ensambla el único evento técnico real de F3 exclusivamente desde FASTQ ya
/// cerrados y validados por el `ExecutionWorkspace`. El manifiesto de evento
/// se vuelve un artefacto más del catálogo antes de que la raíz lo enlace.
fn publish_technical_event(
    workspace: &mut ExecutionWorkspace,
    effective: &EffectiveConfig,
    event: &crate::application::EffectiveTechnicalEvent,
    reference: &Reference,
    reference_fasta: &ArtifactRecord,
    initial_individuals: Option<&InitialIndividualsMaterialization>,
) -> Result<ManifestEventLink, Box<TechnicalEventAssemblyFailure>> {
    let identities = derive_technical_identity_plan(
        effective.reproducibility().root_seed(),
        reference.id(),
        event.event_id().clone(),
        event.individual_ids(),
    )
    .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
    let individuals = initial_individuals
        .ok_or_else(|| {
            technical_event_error("initial_individuals", "faltan individuos materializados")
        })
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?
        .individuals()
        .iter()
        .map(|value| {
            (
                value.individual().id().clone(),
                (value.individual(), value.population_id()),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let libraries = identities.libraries().collect::<Vec<_>>();
    let profile = resolve_f4_technical_profile(event.profile()).map_err(|error| {
        Box::new(TechnicalEventAssemblyFailure {
            cause: technical_event_error("technical_profile", error.to_string()),
            link: None,
        })
    })?;
    let mut udis = BTreeMap::new();
    let demultiplexing_libraries = event
        .libraries()
        .iter()
        .map(|declared| {
            let library = identities
                .library(declared.individual_id())
                .ok_or_else(|| {
                    technical_event_error("technical_event.libraries", "biblioteca técnica ausente")
                })?;
            let udi = Udi::new(declared.i7(), declared.i5()).map_err(|error| {
                technical_event_error("technical_event.udi", format!("{error:?}"))
            })?;
            udis.insert(library.library_id().clone(), udi.clone());
            Ok(DemultiplexingLibrary::new(
                library.library_id().clone(),
                udi,
            ))
        })
        .collect::<Result<Vec<_>, Diagnostic>>()
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;

    // La primera pasada conserva exclusivamente las claves PH y los hechos de
    // quimera; las secuencias se regeneran bajo demanda en la materialización.
    let mut clean_libraries = Vec::with_capacity(libraries.len());
    let mut realized_chimera_pairs = 0_u64;
    for library in &libraries {
        let (individual, population_id) =
            individuals.get(library.individual_id()).ok_or_else(|| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error(
                        "technical_event.individual_ids",
                        "un individuo seleccionado no tiene verdad materializada",
                    ),
                    link: None,
                })
            })?;
        let mut sampler = paired_sampler(reference, individual, &identities, profile, event)
            .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        let mut keys = Vec::new();
        for pair in &mut sampler {
            let pair = pair.map_err(|error| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error("wgs_paired_sampler", format!("{error:?}")),
                    link: None,
                })
            })?;
            keys.push(CleanPairKey::new(
                individual.id().clone(),
                pair.pre_hop_identity().pre_hop_instance_id().clone(),
            ));
        }
        realized_chimera_pairs = realized_chimera_pairs
            .checked_add(sampler.chimera_report().realized_pairs())
            .ok_or_else(|| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error("chimera", "contador de quimeras desbordado"),
                    link: None,
                })
            })?;
        clean_libraries.push(
            CleanLibraryPairs::new(
                individual.id().clone(),
                (*population_id).clone(),
                library.library_id().clone(),
                keys,
            )
            .map_err(|error| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error("biological_contamination", format!("{error:?}")),
                    link: None,
                })
            })?,
        );
    }
    let contamination_plan = plan_biological_contamination(
        effective.reproducibility().root_seed(),
        event.event_id(),
        event.biological_contamination(),
        &clean_libraries,
    )
    .map_err(|error| {
        Box::new(TechnicalEventAssemblyFailure {
            cause: technical_event_error("biological_contamination", format!("{error:?}")),
            link: None,
        })
    })?;
    let mut materializer = BiologicalContaminationMaterializer::new(
        effective.reproducibility().root_seed(),
        event.event_id().clone(),
        &clean_libraries,
    )
    .map_err(|error| {
        Box::new(TechnicalEventAssemblyFailure {
            cause: technical_event_error("biological_contamination", format!("{error:?}")),
            link: None,
        })
    })?;
    let mut error_quality = ErrorQualityAccumulator::default();
    let mut post_contamination = Vec::new();
    let mut biological_manifest = Vec::new();
    for recipient_plan in contamination_plan.recipients() {
        let clean = clean_libraries
            .iter()
            .find(|library| library.individual_id() == recipient_plan.recipient_id())
            .ok_or_else(|| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error("biological_contamination", "receptora ausente"),
                    link: None,
                })
            })?;
        let mut regenerated = RegeneratedPairs::build(
            clean,
            reference,
            &individuals,
            &identities,
            profile,
            event,
            effective.reproducibility().root_seed(),
        )
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        let report = materializer
            .emit_materialized_recipient(
                recipient_plan,
                clean,
                |key| {
                    // Las parejas propias salen de la pasada única; las de una
                    // donante pertenecen a otra serie y se reproducen aparte.
                    regenerated.take(key).map_or_else(
                        || {
                            regenerate_observed_pair(
                                key,
                                reference,
                                &individuals,
                                &identities,
                                profile,
                                event,
                                effective.reproducibility().root_seed(),
                            )
                        },
                        Ok,
                    )
                },
                |materialized| {
                    let report = materialized.payload().report();
                    error_quality.record(report).map_err(|error| {
                        technical_event_error("error_quality", format!("{error:?}"))
                    })?;
                    let recipient_library_id = materialized.recipient_library_id().clone();
                    let udi = udis.get(&recipient_library_id).cloned().ok_or_else(|| {
                        technical_event_error("demultiplexing", "UDI de receptora ausente")
                    })?;
                    let pre_hop_instance_id = materialized.pre_hop_instance_id().clone();
                    let is_contaminant = materialized.is_contaminant();
                    let observed = materialized.into_payload().into_observed();
                    post_contamination.push(PostContaminationPair::new(
                        recipient_library_id,
                        pre_hop_instance_id,
                        udi,
                        is_contaminant,
                        observed,
                    ));
                    Ok(())
                },
            )
            .map_err(|error| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error("biological_contamination", format!("{error:?}")),
                    link: None,
                })
            })?;
        biological_manifest.push(
            recipient_plan
                .realized_manifest_recipient(&report)
                .ok_or_else(|| {
                    Box::new(TechnicalEventAssemblyFailure {
                        cause: technical_event_error(
                            "biological_contamination",
                            "informe realizado incoherente",
                        ),
                        link: None,
                    })
                })?,
        );
    }
    // Las bibliotecas fuera de un subconjunto explícito de receptoras no se
    // pierden: pasan intactas a demultiplexado, regeneradas bajo demanda igual
    // que las propias retenidas de una receptora.
    for clean in &clean_libraries {
        if contamination_plan
            .recipients()
            .iter()
            .any(|plan| plan.recipient_id() == clean.individual_id())
        {
            continue;
        }
        let udi = udis.get(clean.library_id()).cloned().ok_or_else(|| {
            Box::new(TechnicalEventAssemblyFailure {
                cause: technical_event_error("demultiplexing", "UDI de biblioteca ausente"),
                link: None,
            })
        })?;
        let mut regenerated = RegeneratedPairs::build(
            clean,
            reference,
            &individuals,
            &identities,
            profile,
            event,
            effective.reproducibility().root_seed(),
        )
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        for key in clean.own_pairs() {
            let transformed = regenerated.take(key).ok_or_else(|| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error(
                        "biological_contamination",
                        "la clave limpia no se puede regenerar",
                    ),
                    link: None,
                })
            })?;
            error_quality
                .record(transformed.report())
                .map_err(|error| {
                    Box::new(TechnicalEventAssemblyFailure {
                        cause: technical_event_error("error_quality", format!("{error:?}")),
                        link: None,
                    })
                })?;
            post_contamination.push(PostContaminationPair::new(
                clean.library_id().clone(),
                key.pre_hop_instance_id().clone(),
                udi.clone(),
                false,
                transformed.into_observed(),
            ));
        }
    }
    let mut demultiplexer = Demultiplexer::new(
        effective.reproducibility().root_seed(),
        event.event_id().clone(),
        event.index_hopping_rate(),
        demultiplexing_libraries,
    )
    .map_err(|error| {
        Box::new(TechnicalEventAssemblyFailure {
            cause: technical_event_error("demultiplexing", format!("{error:?}")),
            link: None,
        })
    })?;
    let demultiplexed = demultiplexer
        .demultiplex(post_contamination)
        .map_err(|error| {
            Box::new(TechnicalEventAssemblyFailure {
                cause: technical_event_error("demultiplexing", format!("{error:?}")),
                link: None,
            })
        })?;
    let demultiplexing = EventManifestDemultiplexing::from_report(demultiplexed.report().clone());
    let mut final_libraries = demultiplexed.into_final_libraries();
    let mut outputs = Vec::with_capacity(libraries.len() * 2);
    let mut coverage = Vec::with_capacity(libraries.len());
    let mut summaries = Vec::with_capacity(libraries.len());
    for library in &libraries {
        let individual = individuals
            .get(library.individual_id())
            .expect("selected above")
            .0;
        let pairs = final_libraries
            .remove(library.library_id())
            .ok_or_else(|| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error("demultiplexing", "biblioteca final ausente"),
                    link: None,
                })
            })?;
        let publication = write_final_paired_fastq(
            workspace,
            event.event_id().as_str(),
            individual.id().as_str(),
            library.library_id(),
            pairs,
        )
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        let basic_bases = clean_libraries
            .iter()
            .find(|clean| clean.library_id() == library.library_id())
            .and_then(|clean| u64::try_from(clean.own_pairs().len()).ok())
            .and_then(|pairs| pairs.checked_mul(2))
            .and_then(|reads| reads.checked_mul(u64::from(event.read_length())))
            .ok_or_else(|| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error("coverage", "bases WGS básicas desbordan"),
                    link: None,
                })
            })?;
        let target_x = match event.target() {
            EffectiveTechnicalEventTarget::CoverageWgs(value) => value,
            EffectiveTechnicalEventTarget::ReadCount(_) => {
                basic_bases as f64 / reference.total_length() as f64
            }
        };
        let report = F4CoverageReport::from_base_counts(
            target_x,
            basic_bases,
            publication.own_base_count(),
            publication.contaminant_base_count(),
            reference.total_length(),
            publication.n_base_count(),
        )
        .map_err(|error| {
            Box::new(TechnicalEventAssemblyFailure {
                cause: technical_event_error("coverage", error.to_string()),
                link: None,
            })
        })?;
        outputs.extend([
            TechnicalEventArtifact::complete(
                individual.id().clone(),
                library.library_id().clone(),
                crate::technical_identity::TechnicalMate::R1,
                publication.r1().clone(),
            ),
            TechnicalEventArtifact::complete(
                individual.id().clone(),
                library.library_id().clone(),
                crate::technical_identity::TechnicalMate::R2,
                publication.r2().clone(),
            ),
        ]);
        coverage.push(EventManifestCoverage {
            individual_id: individual.id(),
            library_id: library.library_id(),
            report,
        });
        summaries.push(EventManifestFinalFastqSummary {
            individual_id: individual.id(),
            library_id: library.library_id(),
            pair_count: publication.pair_count(),
            read_count: publication.read_count(),
            base_count: publication.base_count(),
            own_pair_count: publication.own_pair_count(),
            contaminant_pair_count: publication.contaminant_pair_count(),
            own_base_count: publication.own_base_count(),
            contaminant_base_count: publication.contaminant_base_count(),
            n_base_count: publication.n_base_count(),
        });
    }
    write_event_manifest(
        workspace,
        event,
        &identities,
        reference_fasta,
        &outputs,
        &coverage,
        &summaries,
        EventManifestErrorQuality::from_report(
            ErrorQualitySettings::default().with_sequencing_errors(event.sequencing_errors()),
            error_quality.report(),
        ),
        demultiplexing,
        realized_chimera_pairs,
        &biological_manifest,
    )
    .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))
}

fn paired_sampler<'a>(
    reference: &'a Reference,
    individual: &'a crate::haplotype::Individual,
    identities: &'a crate::technical_identity::TechnicalIdentityPlan,
    profile: crate::technical_profile::F4TechnicalProfile,
    event: &crate::application::EffectiveTechnicalEvent,
) -> Result<WgsPairedEndSampler<'a>, Diagnostic> {
    let sampler = match event.target() {
        EffectiveTechnicalEventTarget::CoverageWgs(target) => {
            WgsPairedEndSampler::new_for_coverage(
                reference,
                individual,
                identities,
                profile,
                event.read_length(),
                target,
            )
        }
        // La configuración histórica de F3 se mantiene solo como cantidad
        // interna de parejas; la salida sigue siendo PE y el manifiesto deja
        // explícita la cobertura equivalente realizada.
        EffectiveTechnicalEventTarget::ReadCount(pair_count) => {
            WgsPairedEndSampler::new_for_pair_count(
                reference,
                individual,
                identities,
                profile,
                event.read_length(),
                pair_count,
            )
        }
    };
    sampler
        .map(|sampler| sampler.with_effective_technical_event(event))
        .map_err(|error| technical_event_error("wgs_paired_sampler", format!("{error:?}")))
}

/// Reproduce una pareja concreta recorriendo la serie de su individuo desde el
/// principio.
///
/// Solo se usa para las parejas de una donante: son una fracción pequeña y su
/// individuo lo indica la propia clave. Las parejas propias de la receptora, que
/// son la mayoría, salen de [`RegeneratedPairs`] en una sola pasada.
fn regenerate_observed_pair(
    key: &CleanPairKey,
    reference: &Reference,
    individuals: &BTreeMap<
        crate::identity::IndividualId,
        (
            &crate::haplotype::Individual,
            &crate::identity::PopulationId,
        ),
    >,
    identities: &crate::technical_identity::TechnicalIdentityPlan,
    profile: crate::technical_profile::F4TechnicalProfile,
    event: &crate::application::EffectiveTechnicalEvent,
    root_seed: crate::randomness::RootSeed,
) -> Result<ObservedPairTransformation, Diagnostic> {
    let individual = individuals
        .get(key.source_individual_id())
        .map(|(individual, _)| *individual)
        .ok_or_else(|| technical_event_error("biological_contamination", "donante ausente"))?;
    let sampler = paired_sampler(reference, individual, identities, profile, event)?;
    for pair in sampler {
        let pair = pair
            .map_err(|error| technical_event_error("wgs_paired_sampler", format!("{error:?}")))?;
        if pair.pre_hop_identity().pre_hop_instance_id() == key.pre_hop_instance_id() {
            return transform_paired_end(
                root_seed,
                event.event_id(),
                &pair,
                profile,
                ErrorQualitySettings::default().with_sequencing_errors(event.sequencing_errors()),
            )
            .map_err(|error| technical_event_error("error_quality", format!("{error:?}")));
        }
    }
    Err(technical_event_error(
        "biological_contamination",
        "la clave limpia no se puede regenerar",
    ))
}

/// Reproduce en una sola pasada las parejas de una biblioteca, indexadas por
/// identidad.
///
/// `CleanLibraryPairs` ordena sus claves para poder buscarlas por dicotomía, de
/// modo que llegan en orden de identificador y no en el de generación: no sirve
/// avanzar el muestreador clave a clave. Reiniciarlo por cada una, en cambio,
/// recorre la serie entera cada vez y convierte el coste en cuadrático sobre el
/// número de parejas.
///
/// Retener las parejas transformadas no empeora el pico de memoria: el llamador
/// las acumula igualmente en `post_contamination` antes de demultiplexar.
struct RegeneratedPairs {
    pairs: BTreeMap<crate::identity::PreHopInstanceId, ObservedPairTransformation>,
}

impl RegeneratedPairs {
    fn build(
        clean: &CleanLibraryPairs,
        reference: &Reference,
        individuals: &BTreeMap<
            crate::identity::IndividualId,
            (
                &crate::haplotype::Individual,
                &crate::identity::PopulationId,
            ),
        >,
        identities: &crate::technical_identity::TechnicalIdentityPlan,
        profile: crate::technical_profile::F4TechnicalProfile,
        event: &crate::application::EffectiveTechnicalEvent,
        root_seed: crate::randomness::RootSeed,
    ) -> Result<Self, Diagnostic> {
        let individual = individuals
            .get(clean.individual_id())
            .map(|(individual, _)| *individual)
            .ok_or_else(|| {
                technical_event_error("biological_contamination", "individuo ausente")
            })?;
        let wanted = clean
            .own_pairs()
            .iter()
            .map(|key| key.pre_hop_instance_id().clone())
            .collect::<BTreeSet<_>>();
        let sampler = paired_sampler(reference, individual, identities, profile, event)?;
        let mut pairs = BTreeMap::new();
        for pair in sampler {
            let pair = pair.map_err(|error| {
                technical_event_error("wgs_paired_sampler", format!("{error:?}"))
            })?;
            let instance = pair.pre_hop_identity().pre_hop_instance_id();
            if !wanted.contains(instance) {
                continue;
            }
            let instance = instance.clone();
            let transformed = transform_paired_end(
                root_seed,
                event.event_id(),
                &pair,
                profile,
                ErrorQualitySettings::default().with_sequencing_errors(event.sequencing_errors()),
            )
            .map_err(|error| technical_event_error("error_quality", format!("{error:?}")))?;
            pairs.insert(instance, transformed);
        }
        Ok(Self { pairs })
    }

    /// Entrega la pareja de una clave propia desde la pasada única.
    ///
    /// Una clave de otra individua pertenece a una donante: no está en esta
    /// serie y se reproduce aparte. Devolver un error aquí ocultaría ese caso
    /// legítimo de la contaminación.
    fn take(&mut self, key: &CleanPairKey) -> Option<ObservedPairTransformation> {
        self.pairs.remove(key.pre_hop_instance_id())
    }
}

#[allow(clippy::too_many_arguments)] // hechos ya cerrados de las fases F4; agruparlos escondería contratos distintos.
fn write_event_manifest(
    workspace: &mut ExecutionWorkspace,
    event: &crate::application::EffectiveTechnicalEvent,
    identities: &crate::technical_identity::TechnicalIdentityPlan,
    reference_fasta: &ArtifactRecord,
    outputs: &[TechnicalEventArtifact],
    coverage: &[EventManifestCoverage<'_>],
    final_fastq_summary: &[EventManifestFinalFastqSummary<'_>],
    error_quality: crate::event_manifest::EventManifestErrorQuality,
    demultiplexing: EventManifestDemultiplexing,
    realized_chimera_pairs: u64,
    biological_contamination: &[crate::biological_contamination::BiologicalContaminationManifestRecipient],
) -> Result<ManifestEventLink, Diagnostic> {
    let artifacts = outputs
        .iter()
        .map(|output| EventManifestArtifact {
            id: &output.id,
            kind: "fastq",
            state: output.state,
            individual_id: &output.individual_id,
            library_id: &output.library_id,
            mate: output.mate,
            artifact: output.artifact.as_ref(),
            diagnostic: output.diagnostic.as_ref(),
        })
        .collect::<Vec<_>>();
    let document = EventManifest {
        event,
        state: UnitState::Complete,
        identities,
        reference: EventManifestReference {
            reference_id: identities.reference_id(),
            fasta: reference_fasta,
        },
        artifacts: &artifacts,
        coverage,
        final_fastq_summary,
        error_quality,
        demultiplexing,
        chimera_rate: event.chimera_rate(),
        realized_chimera_pairs,
        biological_contamination,
        diagnostics: &[],
    };
    let bytes = event_manifest_bytes(&document)?;
    let relative_path = event_manifest_relative_path(event.event_id());
    let artifact = workspace.write_streamed_artifact_with_validation(
        &format!("event-manifest-{}", event.event_id().as_str()),
        "event_manifest",
        &relative_path,
        |file| {
            use std::io::Write;
            file.write_all(&bytes)
                .map_err(|error| io_error(&relative_path, "write_event_manifest", error))
        },
        |temporary| {
            let bytes = fs::read(temporary)
                .map_err(|error| io_error(temporary, "read_event_manifest", error))?;
            validate_event_manifest_bytes(&relative_path, &bytes)
        },
    )?;
    Ok(ManifestEventLink {
        event_id: event.event_id().as_str().to_owned(),
        state: UnitState::Complete,
        manifest_artifact_id: artifact.id,
        chimera_rate: event.chimera_rate(),
        realized_chimera_pairs,
    })
}

fn technical_event_error(component: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::internal_invariant("technical_event", component).with("detail", detail.into())
}

fn matrix_source_id(source_id: &str) -> &str {
    if source_id.is_empty() {
        "MATRIX_DEFAULT"
    } else {
        source_id
    }
}

fn deduplicate_digests(digests: &mut Vec<NamedDigest>) {
    digests.sort_by(|left, right| left.source.cmp(&right.source));
    digests.dedup_by(|left, right| left.source == right.source);
}

fn named_digest(kind: TruthSourceKind, source_id: &str, sha256: &str) -> NamedDigest {
    NamedDigest {
        source: TruthSourceRef {
            kind,
            source_id: source_id.to_owned(),
        },
        sha256: sha256.to_owned(),
    }
}

fn preflight_prepared_config(
    effective: &EffectiveConfig,
    authorizations: ExecutionAuthorizations,
    prepared_reference: &PreparedReference,
) -> Result<PreparedExecution, Diagnostic> {
    preflight_vcf_samples(effective)?;
    preflight_reductions(effective, prepared_reference.reference())?;
    let prepared_population_snvs =
        prepare_population_snvs_for_execution(effective, authorizations, prepared_reference)?;
    let initial_individuals = preflight_initial_individuals(
        effective,
        prepared_reference.reference(),
        prepared_population_snvs
            .as_ref()
            .map(|prepared| &prepared.matrix),
        authorizations,
    )?;
    Ok(PreparedExecution {
        population_snvs: prepared_population_snvs,
        initial_individuals,
    })
}

/// La proyección F2 solo representa genotipos de individuos iniciales. Un VCF
/// solicitado sin ninguno produciría una cabecera sin muestras, que no es una
/// vista semántica del conjunto de verdad. Se rechaza antes de crear el
/// provisional de publicación, tanto para `validate` como para `run`.
fn preflight_vcf_samples(effective: &EffectiveConfig) -> Result<(), Diagnostic> {
    if effective.vcf_output().is_some()
        && !effective
            .initial_populations()
            .iter()
            .any(|population| !population.individual_ids().is_empty())
    {
        return Err(Diagnostic::configuration(
            "GPS-CONF-008",
            "La salida VCF exige al menos una muestra de una población inicial.",
        )
        .with("component", "vcf_output,initial_populations")
        .with(
            "detail",
            "truth.jsonl no contiene copias ni genotipos cuando initial_populations está vacío",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)] // El cierre DTC-045 necesita todas las relaciones ya validadas.
fn publish_controlled_failure(
    workspace: &mut ExecutionWorkspace,
    execution_id: &str,
    effective: &EffectiveConfig,
    authorizations: ExecutionAuthorizations,
    reductions: &[ManifestReductionLink],
    events: Option<&[ManifestEventLink]>,
    references: &[ManifestReferenceLink],
    cause: &Diagnostic,
) -> Result<(), Diagnostic> {
    let configuration_artifact = workspace
        .artifacts()
        .iter()
        .find(|artifact| artifact.id == "config-general")
        .cloned()
        .filter(|artifact| valid_failure_configuration_copy(workspace, effective, artifact));
    let mut retained = BTreeSet::new();
    if let Some(artifact) = &configuration_artifact {
        retained.insert(artifact.id.clone());
    }
    if events.is_some() {
        retained.extend(
            references
                .iter()
                .map(|reference| reference.fasta_artifact_id.clone()),
        );
        retained.extend(
            workspace
                .artifacts()
                .iter()
                .filter(|artifact| artifact.kind == "fastq")
                .map(|artifact| artifact.id.clone()),
        );
        retained.extend(
            events
                .into_iter()
                .flatten()
                .map(|event| event.manifest_artifact_id.clone()),
        );
    }
    workspace.retain_only_artifacts(&retained)?;
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
        references,
        reduced_references: reductions,
        population_snvs: None,
        vcf_outputs: &[],
        events,
        authorizations,
        artifacts: &artifacts,
        diagnostics: std::slice::from_ref(cause),
    };
    workspace.write_general_manifest(&json_document_bytes(&general_manifest_document(&manifest))?)
}

struct PreparedPopulationSnvsForExecution {
    matrix: crate::population_snvs::PreparedPopulationSnvs,
    publication: PopulationSnvsPublication,
}

struct PreparedExecution {
    population_snvs: Option<PreparedPopulationSnvsForExecution>,
    initial_individuals: Option<InitialIndividualsMaterialization>,
}

enum PopulationSnvsPublication {
    Publish {
        provenance: PopulationSnvsProvenance,
    },
    Reuse {
        metadata_sha256: String,
        matrix_sha256: String,
        normalization_audit: serde_json::Value,
    },
}

fn prepare_population_snvs_for_execution(
    effective: &EffectiveConfig,
    authorizations: ExecutionAuthorizations,
    prepared_reference: &PreparedReference,
) -> Result<Option<PreparedPopulationSnvsForExecution>, Diagnostic> {
    let Some(source) = effective.population_snvs() else {
        return Ok(None);
    };
    let reference = prepared_reference.reference();
    let prepared = match source {
        crate::application::EffectivePopulationSnvs::Generate { parameters } => {
            let matrix = generate_population_snvs_from_reference_with_tolerance(
                parameters,
                reference,
                effective.reproducibility().root_seed(),
                effective.normalization().tolerance(),
            )
            .map_err(|error| population_snvs_error("population_snvs.generate", error.to_string()))?
            .authorize_empty_matrix(authorizations.allow_empty_snv_matrix)
            .map_err(|error| {
                population_snvs_error("population_snvs.generate", error.to_string())
            })?;
            PreparedPopulationSnvsForExecution {
                matrix,
                publication: PopulationSnvsPublication::Publish {
                    provenance: PopulationSnvsProvenance::Generated,
                },
            }
        }
        crate::application::EffectivePopulationSnvs::Supplied { input_path } => {
            let bytes = fs::read(input_path)
                .map_err(|error| io_error(input_path, "read_population_snvs_input", error))?;
            let table = parse_supplied_tsv(&bytes).map_err(|error| {
                population_snvs_error("population_snvs.input_path", error.to_string())
            })?;
            let matrix = prepare_population_snvs_from_reference(
                &table,
                reference,
                effective.normalization().tolerance(),
                authorizations.allow_empty_snv_matrix,
            )
            .map_err(|error| {
                population_snvs_error("population_snvs.input_path", error.to_string())
            })?;
            PreparedPopulationSnvsForExecution {
                matrix,
                publication: PopulationSnvsPublication::Publish {
                    provenance: PopulationSnvsProvenance::Supplied,
                },
            }
        }
        crate::application::EffectivePopulationSnvs::Reuse { metadata_path } => {
            let (matrix, _, metadata_sha256, normalization_audit, matrix_sha256) =
                validate_reusable_population_snvs_with_matrix_sha256(
                    metadata_path,
                    reference,
                    &prepared_reference_fasta_digest(prepared_reference)?,
                    effective.normalization().tolerance(),
                    authorizations.allow_empty_snv_matrix,
                )
                .map_err(|error| {
                    population_snvs_error("population_snvs.metadata_path", error.to_string())
                })?;
            PreparedPopulationSnvsForExecution {
                matrix,
                publication: PopulationSnvsPublication::Reuse {
                    metadata_sha256,
                    matrix_sha256,
                    normalization_audit,
                },
            }
        }
    };
    Ok(Some(prepared))
}

fn publish_prepared_population_snvs(
    prepared: Option<&PreparedPopulationSnvsForExecution>,
    reference_fasta: &ArtifactRecord,
    workspace: &mut ExecutionWorkspace,
) -> Result<Option<PopulationSnvsManifestLink>, Diagnostic> {
    let Some(prepared) = prepared else {
        return Ok(None);
    };
    match &prepared.publication {
        PopulationSnvsPublication::Publish { provenance } => {
            publish_population_snvs(workspace, &prepared.matrix, reference_fasta, *provenance)
                .map(Some)
        }
        PopulationSnvsPublication::Reuse {
            metadata_sha256,
            matrix_sha256: _,
            normalization_audit,
        } => Ok(Some(PopulationSnvsManifestLink {
            mode: "reuse".to_owned(),
            metadata_sha256: metadata_sha256.clone(),
            metadata_artifact_id: None,
            normalization_audit: normalization_audit.clone(),
            initial_population_statistics: None,
        })),
    }
}

fn prepared_reference_fasta_digest(
    prepared_reference: &PreparedReference,
) -> Result<String, Diagnostic> {
    match prepared_reference {
        PreparedReference::Supplied(value) => Ok(value.source_identity().digest().to_owned()),
        PreparedReference::Generated(reference) => {
            let mut bytes = Vec::new();
            write_generated_fasta(reference, &mut bytes)?;
            let mut digest = Sha256Writer::default();
            digest.update(&bytes)?;
            Ok(digest.finish("reference.fa").digest().to_owned())
        }
    }
}

/// Recorre F2.3.2–F2.3.5 antes de abrir el workspace: integra cada fundador,
/// materializa sus copias bajo demanda y deriva el resumen sin alterar matriz ni
/// genotipos. `validate` descarta el resultado; `run` solo publica su resumen.
fn preflight_initial_individuals(
    effective: &EffectiveConfig,
    reference: &Reference,
    matrix: Option<&crate::population_snvs::PreparedPopulationSnvs>,
    authorizations: ExecutionAuthorizations,
) -> Result<Option<InitialIndividualsMaterialization>, Diagnostic> {
    if effective.initial_populations().is_empty() {
        return Ok(None);
    }
    let matrix = matrix.ok_or_else(|| {
        initial_individual_error(
            "direct_events,population_snvs",
            "la matriz SNV falta después de la validación efectiva",
        )
    })?;
    let authorization = if authorizations.allow_incompatible_genotypes {
        IncompatibleGenotypeAuthorization::allowed()
    } else {
        IncompatibleGenotypeAuthorization::denied()
    };
    materialize_initial_individuals(
        effective,
        reference,
        matrix,
        authorization,
        |individual_id| {
            effective
                .direct_events()
                .iter()
                .filter(|event| event.individual_id() == individual_id)
                .map(|event| {
                    declared_small_event(event, reference).map_err(|error| {
                        InitialIndividualsMaterializationError::Sampling(error.render_human())
                    })
                })
                .collect()
        },
    )
    .map(Some)
    .map_err(|error| initial_individual_error("initial_individuals", error.to_string()))
}

pub(crate) fn declared_small_event(
    event: &EffectiveDirectEvent,
    reference: &Reference,
) -> Result<DeclaredSmallEvent, Diagnostic> {
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id() == event.sequence_id())
        .ok_or_else(|| initial_individual_error("direct_events.contig", "contig ausente"))?;
    let provenance = ModificationProvenance::new(event.event_id().clone());
    match event {
        EffectiveDirectEvent::Snp {
            phase,
            sequence_id,
            position_1based,
            reference: expected,
            alternate,
            ..
        } => {
            let position = BasePosition::from_base_1(
                reference.id().clone(),
                sequence_id.clone(),
                *position_1based,
                sequence.length(),
            )
            .map_err(|error| {
                initial_individual_error("direct_events.position_1based", format!("{error:?}"))
            })?;
            let locus = SnpLocus::new(
                position,
                CanonicalBase::from_byte(*expected).map_err(|error| {
                    initial_individual_error("direct_events.ref", format!("{error:?}"))
                })?,
                CanonicalBase::from_byte(*alternate).map_err(|error| {
                    initial_individual_error("direct_events.alt", format!("{error:?}"))
                })?,
            )
            .map_err(|error| initial_individual_error("direct_events", format!("{error:?}")))?;
            Ok(DeclaredSmallEvent::Snp(TargetedSnp::new(
                *phase,
                ApplicableSnp::explicit(locus, event.event_id().clone()),
            )))
        }
        EffectiveDirectEvent::Insertion {
            phase,
            sequence_id,
            point,
            bases,
            ..
        } => Ok(DeclaredSmallEvent::Insertion(
            TargetedCurrentInsertion::new(
                sequence_id.clone(),
                *phase,
                *point,
                bases.clone(),
                provenance,
            ),
        )),
        EffectiveDirectEvent::Deletion {
            phase,
            sequence_id,
            start,
            end,
            ..
        } => Ok(DeclaredSmallEvent::Deletion(TargetedCurrentDeletion::new(
            sequence_id.clone(),
            *phase,
            *start,
            *end,
            provenance,
        ))),
    }
}

fn initial_individual_error(component: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "Los cambios directos iniciales no cumplen el contrato F2.3.4.",
    )
    .with("component", component)
    .with("detail", detail.into())
}

fn publish_population_snvs(
    workspace: &mut ExecutionWorkspace,
    prepared: &crate::population_snvs::PreparedPopulationSnvs,
    reference_fasta: &ArtifactRecord,
    provenance: PopulationSnvsProvenance,
) -> Result<PopulationSnvsManifestLink, Diagnostic> {
    let matrix = workspace.write_artifact(
        "population-snvs-tsv",
        "population_snvs_tsv",
        Path::new("population_snvs.tsv"),
        &canonical_tsv_bytes(prepared),
    )?;
    let metadata = canonical_metadata_document(
        prepared,
        &reference_fasta.digest,
        &matrix.digest,
        provenance,
    );
    let metadata = workspace.write_artifact(
        "population-snvs-metadata",
        "population_snvs_metadata",
        Path::new("population_snvs.metadata.json"),
        &json_document_bytes(&metadata)?,
    )?;
    Ok(PopulationSnvsManifestLink {
        mode: provenance.as_str().to_owned(),
        metadata_sha256: metadata.digest.clone(),
        metadata_artifact_id: Some(metadata.id),
        normalization_audit: normalization_audit_document(prepared),
        initial_population_statistics: None,
    })
}

fn population_snvs_error(component: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La matriz SNV poblacional no cumple el contrato canónico.",
    )
    .with("component", component)
    .with("detail", detail.into())
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
            let composition = effective
                .composition()
                .ok_or_else(|| internal("prepare_reference", "missing composition"))?;
            let frequencies = (*b"ACGT")
                .into_iter()
                .zip(composition.values())
                .map(|(base, probability)| BaseFrequency {
                    base,
                    probability: *probability,
                })
                .collect::<Vec<_>>();
            generate_reference_with_unknown_fraction(
                effective.reference_id(),
                effective.reproducibility().root_seed(),
                &specs,
                &frequencies,
                composition.unknown_base_fraction(),
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
    // Ningún brazo necesita la referencia ya preparada, solo su origen. Se
    // libera antes de releer la copia publicada porque, de lo contrario, la
    // preparada y la autoritativa conviven y el pico dobla el genoma.
    let prepared_was_supplied = matches!(prepared, PreparedReference::Supplied(_));
    drop(prepared);
    let mut authoritative = match (effective.reference_source(), prepared_was_supplied) {
        (EffectiveReferenceSource::Supplied { original_path }, true) => {
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
        (EffectiveReferenceSource::Simulated, false) => {
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
            ExecutionAuthorizations::default(),
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
            ExecutionAuthorizations::default(),
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
    #[ignore = "F4.5.2 sustituye la ruta SE histórica; la vía PE se cubre en f3_event_assembly"]
    #[allow(unreachable_code)]
    fn controlled_failure_after_event_discards_derived_references_with_failed_reductions() {
        let root = std::env::temp_dir().join(format!(
            "gps-f3-event-reduction-failure-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("test root");
        fs::write(
            root.join("reference.fa"),
            format!(">chr1\n{}\n", "A".repeat(200)),
        )
        .expect("reference");
        let config = root.join("general.toml");
        fs::write(
            &config,
            r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[reference]
kind = "supplied"
id = "reference-main"
fasta = "reference.fa"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 200
[[reduced_references]]
id = "reference-small"
source_reference_id = "reference-main"
region_size_bases = 10
target_bases = 100
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["sample-A"]
ploidy = 2
genotype_model = "hwe"
[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 1.0
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
        )
        .expect("config");
        let effective = load_and_validate_config(&config).expect("effective config");
        let destination = root.join("result");

        let error = execute_validated_config_with_effective_document(
            &effective,
            &destination,
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            |_| {
                Err(Diagnostic::internal_invariant(
                    "test",
                    "before-manifest-failure",
                ))
            },
        )
        .expect_err("injected failure");
        assert_eq!(error.code(), "GPS-CONF-008");
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("cleanup");
        return;

        let manifest_path = destination.join("manifest.json");
        let manifest_bytes = fs::read(&manifest_path).expect("failed manifest");
        crate::results::validate_general_manifest_bytes(&manifest_path, &manifest_bytes)
            .expect("controlled failure remains publishable");
        let manifest: serde_json::Value = serde_json::from_slice(&manifest_bytes).expect("JSON");
        assert_eq!(manifest["state"], "failed");
        assert_eq!(
            manifest["references"],
            json!([{
                "id": "reference-main",
                "fasta_artifact_id": "reference-fasta"
            }])
        );
        assert_eq!(manifest["reference_reductions"][0]["state"], "failed");
        assert_eq!(
            manifest["reference_reductions"][0]["fasta_artifact_id"],
            serde_json::Value::Null
        );
        assert!(
            manifest["artifacts"]
                .as_array()
                .expect("catalog")
                .iter()
                .all(
                    |artifact| artifact["id"] != "reference-fasta-reference-small"
                        && artifact["id"] != "reference-map-reference-small"
                )
        );
        assert!(
            destination
                .join("events/reads-main/manifest.json")
                .is_file()
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    #[ignore = "F4.5.2 elimina la validación SE histórica; la publicación PE valida ambos mates"]
    #[allow(unreachable_code)]
    fn failed_independent_fastq_validation_discards_fastq_before_failed_event_manifest() {
        let root = std::env::temp_dir().join(format!(
            "gps-f3-fastq-revalidation-failure-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        if root.exists() {
            fs::remove_dir_all(&root).expect("remove stale test root");
        }
        fs::create_dir(&root).expect("test root");
        let config = root.join("general.toml");
        fs::write(
            &config,
            r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[reference]
kind = "simulated"
id = "reference-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 200
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["sample-A"]
ploidy = 2
genotype_model = "hwe"
[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 1.0
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
        )
        .expect("config");
        let effective = load_and_validate_config(&config).expect("effective config");
        let destination = root.join("result");
        FORCE_FASTQ_REVALIDATION_FAILURE.set(true);
        let error = execute_validated_config(&effective, &destination)
            .expect_err("injected FASTQ validation failure");
        assert_eq!(error.code(), "GPS-CONF-008");
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("cleanup");
        return;
        assert!(!FORCE_FASTQ_REVALIDATION_FAILURE.get());

        let manifest: serde_json::Value = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("failed manifest"),
        )
        .expect("manifest JSON");
        assert_eq!(manifest["state"], "failed");
        assert_eq!(manifest["events"][0]["state"], "failed");
        assert!(
            destination
                .join("events/reads-main/manifest.json")
                .is_file()
        );
        assert!(
            !destination
                .join("outputs/fastq/reads-main/sample-A.fastq.gz")
                .exists(),
            "a FASTQ rejected by independent replay must not be retained"
        );
        assert!(
            manifest["artifacts"]
                .as_array()
                .expect("catalog")
                .iter()
                .all(|artifact| artifact["kind"] != "fastq")
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
            ExecutionAuthorizations::default(),
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

    #[test]
    fn vcf_index_io_failure_publishes_only_the_controlled_failed_result() {
        let root = std::env::temp_dir().join(format!(
            "gps-f2-5-5-index-failure-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("test root");
        let destination = root.join("result");
        let config = root.join("config.toml");
        let mut bytes = fs::read(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/config/general-seeded.toml"),
        )
        .expect("fixture");
        bytes.extend_from_slice(
            b"\n[population_snvs]\nmode = \"generate\"\n[population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 1\nspatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n[[population_snvs.generate.populations]]\npopulation_id = \"P1\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nindividual_ids = [\"IND_01\"]\nploidy = 2\ngenotype_model = \"hwe\"\n[vcf_output]\nmode = \"bgzf\"\n",
        );
        fs::write(&config, bytes).expect("VCF config");
        let effective = load_and_validate_config(&config).expect("valid VCF config");
        let injected = Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
            .with("path", "outputs/truth.vcf.gz.tbi")
            .with("operation", "write_index")
            .with("detail", "injected index writer failure");

        let error = execute_validated_config_with_effective_document(
            &effective,
            &destination,
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            |_| Err(injected.clone()),
        )
        .expect_err("the VCF index I/O failure must fail the publication");

        assert_eq!(error.code(), "GPS-IO-001");
        assert_eq!(error.exit_code(), 5);
        let manifest: serde_json::Value = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("failed manifest"),
        )
        .expect("manifest JSON");
        assert_eq!(manifest["state"], "failed");
        assert_eq!(manifest["diagnostics"][0]["code"], "GPS-IO-001");
        assert!(
            manifest["vcf_outputs"]
                .as_array()
                .expect("VCF outputs")
                .is_empty()
        );
        assert!(!destination.join("outputs/truth.vcf.gz").exists());
        assert!(!destination.join("outputs/truth.vcf.gz.tbi").exists());
        fs::remove_dir_all(root).expect("cleanup");
    }
}
