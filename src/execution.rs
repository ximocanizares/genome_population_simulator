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
use crate::fastq_output::FinalPairedFastqStream;
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
    SequenceExpectation, ValidatedSuppliedReference, reference_from_summary,
    validate_supplied_fasta, validate_supplied_fasta_summary, validate_supplied_fasta_with_origin,
    write_generated_fasta,
};
use crate::reference_generation::{
    BaseFrequency, DeclaredRegion, GeneratedSequenceSpec, generate_reference_with_regions,
};
use crate::reference_reduction::{
    EligibleSequence, MaterializedPlan, ReductionError, materialize_selection,
    preflight_reduced_reference, select_reduced_reference,
};
use crate::reference_reduction_io::{
    REDUCTION_MAP_DOCUMENT_KIND, REDUCTION_MAP_SCHEMA_VERSION, ReductionMapParameters,
    build_reduced_reference, reduction_map_document, source_for, validate_reduction_map_document,
};
use crate::reference_tract_insertion::{
    COORDINATE_MAP_DOCUMENT_KIND, CoordinateMap, DeclaredInsertion, coordinate_map_document,
    insert_tracts,
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
use crate::tandem_repeat::{TandemRepeatTract, observe_tract};
use crate::technical_identity::derive_technical_identity_plan;
use crate::technical_profile::resolve_f4_technical_profile;
use crate::truth_projection::{
    TruthProjectionError, TruthProjectionInput, project_truth, write_truth_artifact,
};
use crate::truth_schema::{
    NamedDigest, TRUTH_COORDINATE_CONVENTION, TRUTH_DOCUMENT_KIND, TRUTH_SCHEMA_VERSION,
    TruthDocument, TruthHeader, TruthSourceKind, TruthSourceRef,
};
use crate::truth_validation::reference_content_sha256;
use crate::vcf_publication::{TractOutputCap, publish_truth_vcf_within};
use crate::wgs_sampling::{WgsPairedEndSampler, WgsSamplingError};

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
    // `validate` tiene que ver **lo mismo** que `run`, y eso incluye el abanico
    // que el deslizamiento biológico genera: sin generarlo aquí, un escenario
    // con loci generados pasaba la validación con los loci vacíos y solo
    // fallaba al ejecutar, que es tarde.
    let generated = effective.with_generated_repeat_alleles(prepared_reference.reference())?;
    let effective: &EffectiveConfig = generated.as_ref().unwrap_or(effective);
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
    // **El generador de alelos por deslizamiento corre aquí, y no antes.** El
    // alelo del que parte cada linaje son las bases que la referencia lleva en
    // el locus, y la referencia no existe hasta ahora. De aquí en adelante
    // manda la configuración con los alelos dentro; si ningún locus declara
    // generador, `generated` es `None` y se sigue usando la original, que es lo
    // que garantiza que un escenario que no lo usa publique los mismos bytes.
    //
    // Sigue siendo **antes** de crear el provisional, de modo que un rechazo
    // del generador no deja nada escrito.
    let generated = effective.with_generated_repeat_alleles(prepared_reference.reference())?;
    let effective: &EffectiveConfig = generated.as_ref().unwrap_or(effective);
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
    if let Some(crate::application::EffectivePopulationSnvs::Reuse { metadata_path }) =
        effective.population_snvs()
    {
        workspace.set_reused_snv_metadata(metadata_path);
    }
    let mut configuration_input = None;
    let mut reference_artifact = None;
    let mut manifest_references = Vec::new();
    let mut manifest_reductions = Vec::new();
    let mut published_reductions = Vec::with_capacity(effective.reduced_references().len());
    let mut population_snvs = None;
    let mut vcf_outputs = Vec::new();
    let mut event_links: Option<Vec<ManifestEventLink>> = None;
    let mut realized_tract_output = crate::results::RealizedTractOutput::default();
    let mut inserted_publication: Option<PublishedInsertion> = None;
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
        // **Con inserción, la referencia final se publica junto a la copia
        // suministrada** (DTC-098), y es la que reciben todos los consumidores
        // de aquí en adelante.
        inserted_publication = write_inserted_reference(
            prepared_reference
                .as_ref()
                .ok_or_else(|| internal("execute", "reference already released"))?,
            &mut workspace,
        )?;
        if let Some(inserted) = &inserted_publication {
            manifest_references.push(ManifestReferenceLink {
                id: effective.simulation_reference_id().to_owned(),
                fasta_artifact_id: inserted.fasta.id.clone(),
            });
        }
        prepared_reference = Some(load_authoritative_reference(
            effective,
            &workspace,
            source_artifact,
            inserted_publication
                .as_ref()
                .map(|inserted| &inserted.fasta),
            prepared_reference
                .take()
                .ok_or_else(|| internal("execute", "reference already released"))?,
        )?);
        let prepared = prepared_reference
            .as_ref()
            .ok_or_else(|| internal("execute", "reference already released"))?;
        let reference = prepared.reference();
        let declared_reference = prepared.declared_reference();
        let simulation_artifact = inserted_publication
            .as_ref()
            .map_or(source_artifact, |inserted| &inserted.fasta);
        population_snvs = publish_prepared_population_snvs(
            prepared_execution.population_snvs.as_ref(),
            simulation_artifact,
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
        let source_fasta = workspace.root().join(&source_artifact.path);
        let prepared_reductions = prepare_reductions(effective, declared_reference, &source_fasta)?;
        manifest_reductions = manifest_reduction_links(effective, &prepared_reductions);
        for reduction in &mut manifest_reductions {
            reduction.source_fasta_artifact_id = Some(source_artifact.id.clone());
        }
        let reference_output = ReferenceOutput {
            reference: declared_reference,
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
                    source_for(reference_output.reference, Some(&source_fasta))?.as_ref(),
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
            simulation_artifact,
            configuration_input
                .as_ref()
                .ok_or_else(|| internal("execute", "missing configuration artifact"))?,
            effective,
        )?;
        // Las bases de tracto realizadas se cuentan **donde se producen**: la
        // de la verdad sobre el documento recién construido, la del VCF sobre
        // los registros recién escritos. Ninguna de las dos pasa por la cota.
        realized_tract_output.truth_bases = crate::results::realized_truth_tract_bases(&truth);
        if let Some(mode) = effective.vcf_output() {
            let cap = effective
                .repeat_output_budget()
                .declared_cap()
                .map(|cap| TractOutputCap { cap });
            let published =
                publish_truth_vcf_within(&mut workspace, reference, &truth, mode, cap.as_ref())?;
            realized_tract_output.vcf_bases = published.tract_bases;
            realized_tract_output.vcf_by_locus = published.tract_bases_by_locus;
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
                simulation_artifact,
                prepared_execution.initial_individuals.as_ref(),
            ) {
                Ok(link) => event_links = Some(vec![link]),
                Err(failure) => {
                    event_links = failure.link.map(|link| vec![link]);
                    return Err(failure.cause);
                }
            }
        }
        let mut effective_document = effective_config_document_with_population_snvs(
            effective,
            authorizations,
            reference_output,
            &reduction_outputs,
            population_snvs.as_ref(),
        )?;
        crate::results::publish_measured_vcf_tract_bases(
            &mut effective_document,
            &realized_tract_output,
        );
        if let Some(inserted) = &inserted_publication {
            effective_document["inserted_reference"] = crate::results::inserted_reference_json(
                effective,
                reference,
                &inserted.fasta,
                &inserted.map,
                source_artifact,
            )?;
        }
        let effective_document = transform_effective(effective_document);
        // Ningún artefacto posterior necesita ya las bases de la referencia. Se
        // libera antes de escribir el manifiesto, porque esa escritura dispara
        // la validación física, que relee el FASTA publicado y construye otra
        // referencia completa. Sin esto ambas conviven y el pico de memoria
        // dobla el tamaño del genoma de partida.
        prepared_reference = None;
        // Tampoco hace falta ya la verdad en memoria: el juez de los recursos
        // publicados relee la del paquete. Viva, convivía con la releída y con
        // muchos casos subía el pico de memoria un 41 % (DT-020).
        drop(truth);
        let effective_record = workspace.write_artifact(
            "effective-config",
            "effective_config",
            Path::new("effective-config.json"),
            &json_document_bytes(&effective_document)?,
        )?;
        before_manifest(&mut workspace)?;

        let artifacts = workspace.artifacts().to_vec();
        // El sorteo de los alelos de repetición ya ocurrió en el preflight; sus
        // resultados son lo único que permite publicar la frecuencia
        // **realizada** en vez de repetir la declarada.
        //
        // **La sección se construye desde la configuración, no desde el
        // sorteo.** Un escenario puede declarar loci y no tener ninguna
        // población inicial, y entonces no hay una sola copia sorteada; el
        // manifiesto tiene que registrar igualmente los loci y sus parámetros,
        // con cero copias y sin frecuencia realizada. Atarlo a que existieran
        // individuos materializados dejaba ese escenario sin registro y
        // emparejaba una configuración efectiva con loci contra un manifiesto
        // sin ellos.
        let repeat_loci = crate::results::repeat_loci_manifest_json(
            effective,
            prepared_execution
                .initial_individuals
                .as_ref()
                .map_or(&[][..], InitialIndividualsMaterialization::repeat_draws),
            &realized_tract_output,
        )?;
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
            repeat_loci: repeat_loci.as_ref(),
        };
        let mut manifest_document = general_manifest_document(&manifest);
        if let Some(inserted) = &inserted_publication {
            manifest_document["reference_insertion"] = crate::results::manifest_insertion_json(
                effective,
                UnitState::Complete,
                &inserted.fasta,
                &inserted.map,
                source_artifact,
            )?;
        }
        let manifest_bytes = json_document_bytes(&manifest_document)?;
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
            // **Con lecturas ya publicadas, la referencia final se conserva con
            // su mapa** (DTC-098): las lecturas salieron de ella, y el paquete
            // fallido tiene que poder relacionarla con la suministrada
            // (REQ-MAN-07). Sin evento no se conserva nada de ella.
            let failure_insertion = event_links
                .as_ref()
                .and(inserted_publication.as_ref())
                .zip(reference_artifact.as_ref());
            publish_controlled_failure(
                &mut workspace,
                execution_id.as_str(),
                effective,
                authorizations,
                &manifest_reductions,
                event_links.as_deref(),
                &failure_references,
                failure_insertion,
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
            // El alelo de un locus de repetición lo **sortea** el simulador de
            // las frecuencias declaradas, de modo que su procedencia es
            // aleatoria y su huella es la del fichero que declaró esas
            // frecuencias. Un escenario sin loci no añade ninguna fila, y por
            // eso su cabecera sigue siendo byte a byte la de antes.
            let applied_repeat_loci = individuals.applied_repeat_locus_ids();
            source_digests.extend(applied_repeat_loci.iter().map(|locus_id| {
                named_digest(
                    TruthSourceKind::Random,
                    locus_id,
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
            .map_err(|error| match error {
                // Que la verdad no describa la copia simulada es un defecto
                // del programa, no un fallo de formato de la entrada.
                TruthProjectionError::IndependentReconstruction { .. }
                | TruthProjectionError::ReconstructionMismatch { .. }
                | TruthProjectionError::CopySetMismatch
                | TruthProjectionError::GenotypeJudgement(_) => {
                    internal("truth_reconstruction", error)
                }
                // Un documento que el propio programa acaba de proyectar y
                // que incumple el contrato de la verdad —un portador de otro
                // individuo, un identificador fuera de léxico, dos casos con
                // el mismo identificador, un alelo fuera de su catálogo— es un
                // defecto del productor, no un error de la entrada (DT-021).
                // La lectura de una verdad ya escrita no pasa por aquí.
                //
                // Solo `Contract` lo ejercita una prueba. `DuplicateCaseId` y
                // `RepeatAlleleOutsideItsCase` **no se alcanzan**: la
                // proyección arma cada catálogo y cada identificador con las
                // mismas copias que después consulta, y ninguna entrada puede
                // hacerlos divergir. Están aquí para que, si algún día se
                // alcanzan, salgan como lo que serían, un defecto interno.
                // Devolverlos a `GPS-FMT-006` no lo detecta ninguna prueba.
                TruthProjectionError::Contract(_)
                | TruthProjectionError::DuplicateCaseId(_)
                | TruthProjectionError::RepeatAlleleOutsideItsCase(_) => {
                    internal("truth_projection", error)
                }
                TruthProjectionError::BasesNotMaterialized
                | TruthProjectionError::MissingMatrix
                | TruthProjectionError::MatrixReferenceMismatch
                | TruthProjectionError::MatrixSourceMismatch
                | TruthProjectionError::MissingDigest(_)
                | TruthProjectionError::InvalidUtf8Insertion { .. }
                | TruthProjectionError::Genotype(_)
                | TruthProjectionError::Materialization(_)
                | TruthProjectionError::MissingRepeatSequence(_) => {
                    Diagnostic::format("GPS-FMT-006", "No se pudo proyectar la verdad genómica.")
                        .with("detail", error.to_string())
                }
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

/// Huella del contenido de la referencia.
///
/// Cuando se validó sin conservar las bases, la huella se calculó durante esa
/// misma pasada y viaja con ella: recalcularla exigiría volver a leer el FASTA.
fn reference_content_digest(reference: &Reference) -> Result<String, Diagnostic> {
    reference.content_sha256().map_or_else(
        || {
            reference_content_sha256(reference)
                .map_err(|error| internal("truth_reference_digest", format!("{error:?}")))
        },
        |digest| Ok(digest.to_owned()),
    )
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
            reference_content_sha256: reference_content_digest(reference)?,
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

    // El padrón fija el puesto de cada individua antes de recorrer nada, de
    // modo que el orden del inventario no dependa de en qué orden se generen
    // las parejas. Es lo que permite que una clave guarde un puesto en lugar de
    // repetir el identificador de su individua un millón de veces.
    let roster = crate::identity::IndividualRoster::new(individuals.keys().cloned())
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
    // El plan de tartamudeo se resuelve **una vez** y se presta a cada
    // muestreador. Las parejas se regeneran varias veces —por la
    // contaminación—, y un plan reconstruido en cada pasada podría no ser el
    // mismo objeto pero sí tendría que dar la misma decisión; resolverlo aquí
    // hace que esa igualdad no dependa de nadie.
    let stutter_plan = effective.stutter_plan();
    // El plan de amplificación se resuelve **una vez**, por la misma razón que
    // el de tartamudeo: las parejas se regeneran varias veces y las dos
    // pasadas tienen que abrir exactamente las mismas familias.
    let pcr_plan = effective.pcr_plan();
    // La primera pasada conserva exclusivamente las claves PH y los hechos de
    // quimera; las secuencias se regeneran bajo demanda en la materialización.
    let mut clean_libraries = Vec::with_capacity(libraries.len());
    let mut realized_chimera_pairs = 0_u64;
    let mut stutter_accumulator = crate::stutter::StutterAccumulator::default();
    let mut pcr_accumulator = crate::pcr_amplification::PcrAccumulator::default();
    // La cota se acumula **antes de iterar** cada muestreador, desde la
    // configuración y la geometría ya conocidas: es un valor calculado sin
    // ejecutar, no una medida.
    let mut pcr_bound = crate::pcr_amplification::PcrResourceBound::empty();
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
        let seat = roster.seat(individual.id()).ok_or_else(|| {
            Box::new(TechnicalEventAssemblyFailure {
                cause: technical_event_error(
                    "technical_event.individual_ids",
                    "una individua seleccionada no figura en el padrón",
                ),
                link: None,
            })
        })?;
        let mut sampler = paired_sampler(
            reference,
            individual,
            &identities,
            profile,
            event,
            &stutter_plan,
            &pcr_plan,
        )
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        if let Some(settings) = pcr_plan.settings() {
            let haplotype_copies = individual
                .sequences()
                .iter()
                .map(|genotype| genotype.copies().len() as u64)
                .sum::<u64>();
            pcr_bound = pcr_bound.combine(crate::pcr_amplification::resource_bound(
                settings,
                sampler
                    .coverage_plan()
                    .map_or(0, crate::wgs_sampling::WgsPairedEndPlan::pair_count),
                haplotype_copies,
                u64::from(profile.insert_length().maximum_bases()),
                // Los tractos de repetición **no cuentan todavía**: la
                // preparación amplificada rechaza por ahora los escenarios que
                // declaran loci, así que este término es cero hasta que el
                // tartamudeo clonal exista.
                0,
            ));
        }
        let mut keys = Vec::new();
        for pair in &mut sampler {
            let pair = pair.map_err(|error| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: wgs_paired_sampler_error(error),
                    link: None,
                })
            })?;
            keys.push(
                CleanPairKey::new(seat, pair.pre_hop_identity().pre_hop_instance_id()).map_err(
                    |cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }),
                )?,
            );
        }
        // El informe del tartamudeo se recoge **de la primera pasada**, igual
        // que el de la quimera: las regeneraciones posteriores repiten las
        // mismas decisiones y contarlas otra vez las duplicaría.
        stutter_accumulator
            .merge(&sampler.stutter_report())
            .map_err(|error| {
                Box::new(TechnicalEventAssemblyFailure {
                    cause: technical_event_error("stutter", format!("{error:?}")),
                    link: None,
                })
            })?;
        // Los recuentos de la amplificación se recogen **de la primera
        // pasada**, igual que los de quimera y tartamudeo: las regeneraciones
        // posteriores repiten las mismas familias y contarlas otra vez las
        // duplicaría.
        pcr_accumulator.merge(&sampler.pcr_report());
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
                seat,
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
    let mut biological_manifest = Vec::new();
    // El demultiplexado y la escritura final se preparan antes de generar: la
    // biblioteca final de una pareja depende solo de su biblioteca e instancia
    // anteriores al salto, de modo que cada pareja puede clasificarse y
    // escribirse en cuanto existe. Así el pico no crece con la cobertura.
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
    let mut streams: BTreeMap<crate::identity::LibraryId, FinalPairedFastqStream> = BTreeMap::new();
    for library in &libraries {
        let individual = individuals
            .get(library.individual_id())
            .expect("selected above")
            .0;
        let udi = udis.get(library.library_id()).ok_or_else(|| {
            Box::new(TechnicalEventAssemblyFailure {
                cause: technical_event_error("technical_event.udi", "biblioteca sin UDI declarado"),
                link: None,
            })
        })?;
        let stream = FinalPairedFastqStream::open(
            workspace,
            event.event_id().as_str(),
            individual.id().as_str(),
            library.library_id(),
            udi,
            event.read_length(),
        )
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        streams.insert(library.library_id().clone(), stream);
    }
    // Las parejas que aporta una donante se conocen antes de materializar: el
    // plan ya dice cuál sustituye a cuál. Se agrupan por individua para recorrer
    // cada serie una sola vez, en lugar de recorrerla entera por cada pareja
    // contaminada, que volvería el coste cuadrático igual que en las propias.
    let mut donor_requests: BTreeMap<
        crate::identity::IndividualSeat,
        BTreeSet<crate::identity::CompactPreHopInstanceId>,
    > = BTreeMap::new();
    for recipient_plan in contamination_plan.recipients() {
        for assignment in recipient_plan.assignments() {
            let source = assignment.source_pair();
            donor_requests
                .entry(source.source_individual())
                .or_default()
                .insert(source.compact_pre_hop_instance_id());
        }
    }
    let mut donor_pairs = BTreeMap::new();
    for (seat, wanted) in &donor_requests {
        let individual_id = roster.individual_id(*seat).ok_or_else(|| {
            Box::new(TechnicalEventAssemblyFailure {
                cause: technical_event_error(
                    "biological_contamination",
                    "una donante no figura en el padrón",
                ),
                link: None,
            })
        })?;
        let collected = RegeneratedPairs::collect(
            individual_id,
            wanted,
            reference,
            &individuals,
            &identities,
            profile,
            event,
            effective.reproducibility().root_seed(),
            &stutter_plan,
            &pcr_plan,
        )
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        donor_pairs.insert(*seat, collected);
    }
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
            &stutter_plan,
            &pcr_plan,
        )
        .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        let report = materializer
            .emit_materialized_recipient(
                recipient_plan,
                clean,
                |key| {
                    // Las parejas propias salen de la pasada única de esta
                    // biblioteca; las de una donante, de la serie que ya se
                    // recorrió para su individua.
                    regenerated
                        .take(key)
                        .or_else(|| {
                            donor_pairs
                                .get(&key.source_individual())
                                .and_then(|pairs| pairs.get(&key.compact_pre_hop_instance_id()))
                                .cloned()
                        })
                        .ok_or_else(|| {
                            technical_event_error(
                                "biological_contamination",
                                "la clave limpia no se puede regenerar",
                            )
                        })
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
                    write_demultiplexed_pair(
                        &mut demultiplexer,
                        &mut streams,
                        PostContaminationPair::new(
                            recipient_library_id,
                            pre_hop_instance_id,
                            udi,
                            is_contaminant,
                            observed,
                        ),
                    )
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
            &stutter_plan,
            &pcr_plan,
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
            write_demultiplexed_pair(
                &mut demultiplexer,
                &mut streams,
                PostContaminationPair::new(
                    clean.library_id().clone(),
                    key.pre_hop_instance_id().clone(),
                    udi.clone(),
                    false,
                    transformed.into_observed(),
                ),
            )
            .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))?;
        }
    }
    let demultiplexing = EventManifestDemultiplexing::from_report(demultiplexer.take_report());
    let mut outputs = Vec::with_capacity(libraries.len() * 2);
    let mut coverage = Vec::with_capacity(libraries.len());
    let mut summaries = Vec::with_capacity(libraries.len());
    for library in &libraries {
        let individual = individuals
            .get(library.individual_id())
            .expect("selected above")
            .0;
        let stream = streams.remove(library.library_id()).ok_or_else(|| {
            Box::new(TechnicalEventAssemblyFailure {
                cause: technical_event_error("demultiplexing", "biblioteca final ausente"),
                link: None,
            })
        })?;
        let publication = stream
            .finish(workspace)
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
        // La sección solo existe cuando el evento abre de verdad la etapa
        // `stutter`: con el tartamudeo apagado, o sin loci declarados, el
        // manifiesto sale idéntico al de antes de esta ampliación.
        stutter_plan
            .is_active()
            .then(|| crate::event_manifest::EventManifestStutter {
                settings: stutter_plan.settings().clone(),
                report: stutter_accumulator.report(),
                loci: stutter_plan
                    .loci()
                    .iter()
                    .map(|locus| {
                        (
                            locus.id().as_str().to_owned(),
                            locus
                                .motifs()
                                .iter()
                                .map(crate::tandem_repeat::RepeatMotif::length)
                                .collect(),
                        )
                    })
                    .collect(),
            }),
        // La sección solo existe con preparación amplificada: con la
        // preparación sin PCR, `artifacts.pcr` conserva exactamente los bytes
        // de antes de esta ampliación.
        event.pcr().map(|settings| {
            let report = pcr_accumulator.clone().report();
            crate::event_manifest::EventManifestPcr {
                settings: settings.clone(),
                report,
                bound: pcr_bound,
                basic_wgs_coverage_x: coverage
                    .iter()
                    .map(|entry| entry.report.basic_wgs_x())
                    .sum(),
            }
        }),
    )
    .map_err(|cause| Box::new(TechnicalEventAssemblyFailure { cause, link: None }))
}

/// Clasifica una pareja y la escribe en el FASTQ de su biblioteca final.
///
/// Una pareja descartada por UDI ajeno a la tanda no se escribe: solo queda
/// contada en los agregados del demultiplexado, igual que en la vía anterior.
fn write_demultiplexed_pair(
    demultiplexer: &mut Demultiplexer,
    streams: &mut BTreeMap<crate::identity::LibraryId, FinalPairedFastqStream>,
    pair: PostContaminationPair<crate::error_quality::ObservedPair>,
) -> Result<(), Diagnostic> {
    let Some(demultiplexed) = demultiplexer
        .accept(pair)
        .map_err(|error| technical_event_error("demultiplexing", format!("{error:?}")))?
    else {
        return Ok(());
    };
    streams
        .get_mut(demultiplexed.final_library_id())
        .ok_or_else(|| {
            technical_event_error("demultiplexing", "biblioteca final sin flujo abierto")
        })?
        .write_pair(demultiplexed)
}

fn paired_sampler<'a>(
    reference: &'a Reference,
    individual: &'a crate::haplotype::Individual,
    identities: &'a crate::technical_identity::TechnicalIdentityPlan,
    profile: crate::technical_profile::F4TechnicalProfile,
    event: &crate::application::EffectiveTechnicalEvent,
    stutter: &crate::stutter::StutterPlan,
    pcr: &crate::pcr_amplification::PcrPlan,
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
        .map(|sampler| {
            sampler
                .with_effective_technical_event(event)
                .with_stutter_plan(stutter.clone())
                .with_pcr_plan(pcr.clone())
        })
        // El rechazo de admisión de DT-006 (`ImpracticalInsertSampling`) nunca
        // llega aquí: se dispara dentro de `sample_next_pair`, al iterar, no
        // al construir el muestreador. Los dos puntos donde sí se puede
        // alcanzar usan `wgs_paired_sampler_error`.
        .map_err(|error| technical_event_error("wgs_paired_sampler", format!("{error:?}")))
}

/// Traduce el fallo de una iteración del muestreador PE a un diagnóstico
/// público; es el único punto de conversión para las dos pasadas que
/// consumen el iterador (la de claves y la de regeneración final), de modo
/// que las dos ven el mismo tratamiento de DT-006.
///
/// El rechazo de admisión de DT-006 —`ImpracticalInsertSampling`, un sorteo
/// del inserto que agotó su cota de intentos— necesita un diagnóstico claro
/// que nombre al individuo, su copia más larga, el mínimo del perfil y la
/// salida disponible, en vez de volcar la variante de Rust que la produjo.
/// El resto de variantes de `WgsSamplingError` no debería alcanzarse nunca en
/// una ejecución con configuración ya validada, así que conservan el
/// tratamiento uniforme de invariante interna.
///
/// Reutiliza `GPS-DOMAIN-003` en vez de abrir un código nuevo: el manifiesto
/// general valida sus diagnósticos contra el catálogo histórico congelado en
/// `diagnostics-v3.tsv` (DTC-036 fija su esquema en la versión `"2"`), que no
/// conoce ningún código añadido después de F5. Un código nuevo se vería
/// rechazado al publicar precisamente el fallo controlado que debía informar
/// —la limitación en sí está registrada aparte, en DT-007—.
fn wgs_paired_sampler_error(error: WgsSamplingError) -> Diagnostic {
    match error {
        WgsSamplingError::ImpracticalInsertSampling {
            individual,
            longest_copy_bases,
            minimum_insert_bases,
        } => Diagnostic::domain(
            "GPS-DOMAIN-003",
            "Una invariante del estado biológico o técnico válido no se cumple.",
        )
        .with("entity", individual.as_str())
        .with("invariant", "inserto_practicable")
        .with(
            "detail",
            format!(
                "el individuo {individual} tiene su copia haplotípica más larga en {longest_copy_bases} pb, \
y el perfil técnico exige un inserto mínimo de {minimum_insert_bases} pb: ese margen hace impracticable el \
sorteo del inserto. Usa una referencia más larga para este individuo."
            ),
        ),
        other => technical_event_error("wgs_paired_sampler", format!("{other:?}")),
    }
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
/// Se indexan por la forma compacta de la instancia, que es la que sostiene el
/// inventario: así el bucle no reconstruye un texto por pareja para buscarla.
struct RegeneratedPairs {
    pairs: BTreeMap<crate::identity::CompactPreHopInstanceId, ObservedPairTransformation>,
}

impl RegeneratedPairs {
    #[allow(clippy::too_many_arguments)] // los mismos hechos cerrados de F4 que `collect`.
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
        stutter: &crate::stutter::StutterPlan,
        pcr: &crate::pcr_amplification::PcrPlan,
    ) -> Result<Self, Diagnostic> {
        let wanted = clean
            .own_pairs()
            .iter()
            .map(CleanPairKey::compact_pre_hop_instance_id)
            .collect::<BTreeSet<_>>();
        Ok(Self {
            pairs: Self::collect(
                clean.individual_id(),
                &wanted,
                reference,
                individuals,
                identities,
                profile,
                event,
                root_seed,
                stutter,
                pcr,
            )?,
        })
    }

    /// Recorre la serie de una individua una sola vez y transforma únicamente
    /// las parejas pedidas.
    #[allow(clippy::too_many_arguments)] // los mismos hechos cerrados de F4 que el resto de la fase.
    fn collect(
        individual_id: &crate::identity::IndividualId,
        wanted: &BTreeSet<crate::identity::CompactPreHopInstanceId>,
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
        stutter: &crate::stutter::StutterPlan,
        pcr: &crate::pcr_amplification::PcrPlan,
    ) -> Result<
        BTreeMap<crate::identity::CompactPreHopInstanceId, ObservedPairTransformation>,
        Diagnostic,
    > {
        let individual = individuals
            .get(individual_id)
            .map(|(individual, _)| *individual)
            .ok_or_else(|| {
                technical_event_error("biological_contamination", "individuo ausente")
            })?;
        let sampler = paired_sampler(
            reference, individual, identities, profile, event, stutter, pcr,
        )?;
        let mut pairs = BTreeMap::new();
        for pair in sampler {
            let pair = pair.map_err(wgs_paired_sampler_error)?;
            let instance = crate::identity::CompactPreHopInstanceId::new(
                pair.pre_hop_identity().pre_hop_instance_id(),
            )?;
            if !wanted.contains(&instance) {
                continue;
            }
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
        Ok(pairs)
    }

    /// Entrega la pareja de una clave propia desde la pasada única.
    ///
    /// Una clave de otra individua pertenece a una donante: no está en esta
    /// serie y se reproduce aparte. Devolver un error aquí ocultaría ese caso
    /// legítimo de la contaminación.
    fn take(&mut self, key: &CleanPairKey) -> Option<ObservedPairTransformation> {
        self.pairs.remove(&key.compact_pre_hop_instance_id())
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
    stutter: Option<crate::event_manifest::EventManifestStutter>,
    pcr: Option<crate::event_manifest::EventManifestPcr>,
) -> Result<ManifestEventLink, Diagnostic> {
    let stutter_present = stutter.is_some();
    let pcr_present = pcr.is_some();
    let profile_id = event.profile();
    let profile_sha256 = match event.preparation() {
        crate::technical_profile::LibraryPreparation::PcrFree => {
            crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256
        }
        crate::technical_profile::LibraryPreparation::Pcr => {
            crate::technical_profile::GPS_HISEQX_PCR_PE150_SHA256
        }
    };
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
        stutter,
        pcr,
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
        profile_id,
        profile_sha256,
        rng_catalog: if pcr_present {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V5
        } else if stutter_present {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
        } else {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
        },
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
    admit_repeat_loci_against_small_variants(
        effective,
        prepared_population_snvs
            .as_ref()
            .map(|prepared| &prepared.matrix),
    )?;
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
    insertion: Option<(&PublishedInsertion, &ArtifactRecord)>,
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
        if let Some((inserted, _)) = insertion {
            retained.insert(inserted.map.id.clone());
        }
    }
    workspace.retain_only_artifacts(&retained)?;
    let artifacts = workspace.artifacts().to_vec();
    let configuration_input = configuration_artifact.map(|copy| ConfigurationInputRecord {
        kind: "general".to_owned(),
        original_path: effective.configuration_input().original_path().to_owned(),
        copy,
    });
    // Una ejecución fallida no ha publicado tracto alguno: sus cifras
    // realizadas son cero y la cota sigue siendo la que la configuración fijó.
    let repeat_loci = crate::results::repeat_loci_manifest_json(
        effective,
        &[],
        &crate::results::RealizedTractOutput::default(),
    )?;
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
        // Una ejecución fallida no ha sorteado nada, pero **sí ha declarado**
        // los loci, y el manifiesto los registra igual que en un éxito: con
        // cero copias sorteadas y sin frecuencia realizada. Lo que se registra
        // es la configuración, que existe tanto si la ejecución terminó como
        // si no.
        repeat_loci: repeat_loci.as_ref(),
    };
    let mut document = general_manifest_document(&manifest);
    if let Some((inserted, source_fasta)) = insertion {
        document["reference_insertion"] = crate::results::manifest_insertion_json(
            effective,
            UnitState::Failed,
            &inserted.fasta,
            &inserted.map,
            source_fasta,
        )?;
    }
    workspace.write_general_manifest(&json_document_bytes(&document)?)
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
        PreparedReference::Generated(reference)
        | PreparedReference::Inserted {
            derived: reference, ..
        } => {
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
    .map_err(|error| match error {
        // El rechazo de un locus de repetición viaja con su diagnóstico
        // entero y se publica tal cual: envolverlo dentro del detalle de otro
        // era lo que producía el mensaje ilegible.
        InitialIndividualsMaterializationError::Repeat(diagnostic) => diagnostic,
        other => initial_individual_error("initial_individuals", other.to_string()),
    })
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
    /// Referencia suministrada con tractos insertados (DTC-098). La
    /// **simulación** usa la final; la suministrada solo se conserva por su
    /// catálogo e identidad, sin bases, para publicar su copia y describirla.
    Inserted {
        supplied: ValidatedSuppliedReference,
        derived: Reference,
        map: CoordinateMap,
    },
}

impl PreparedReference {
    /// La referencia **de la simulación**: la final cuando hay inserción.
    fn reference(&self) -> &Reference {
        match self {
            Self::Supplied(value) => value.reference(),
            Self::Generated(value) => value,
            Self::Inserted { derived, .. } => derived,
        }
    }

    fn reference_mut(&mut self) -> &mut Reference {
        match self {
            Self::Supplied(value) => value.reference_mut(),
            Self::Generated(value) => value,
            Self::Inserted { derived, .. } => derived,
        }
    }

    /// La referencia **declarada**: la suministrada aunque haya inserción. Es
    /// la que describe la sección `reference` de la configuración efectiva.
    fn declared_reference(&self) -> &Reference {
        match self {
            Self::Inserted { supplied, .. } => supplied.reference(),
            _ => self.reference(),
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

/// Indica si la configuración necesita las bases de la referencia principal.
///
/// Una reducción no las necesita para planificar: le bastan identificador,
/// clase y longitud, y pide después los tramos elegidos. Cualquier otra salida
/// —matriz SNV, individuos, eventos directos, evento técnico o VCF— trabaja
/// sobre la secuencia y obliga a materializarla.
///
/// **Los loci de repetición también obligan.** La guarda de admisión tiene que
/// leer el tramo declarado de la referencia para comprobar que allí hay de
/// verdad un tracto de los motivos declarados, y sin bases no hay nada que
/// mirar. Mientras la guarda no existió, esta línea se dejó fuera a propósito:
/// materializar el genoma entero sin que nadie lo leyera habría sido pagar la
/// memoria por nada.
fn requires_materialized_bases(effective: &EffectiveConfig) -> bool {
    effective.population_snvs().is_some()
        || !effective.initial_populations().is_empty()
        || !effective.direct_events().is_empty()
        || effective.technical_event().is_some()
        || effective.vcf_output().is_some()
        || !effective.repeat_loci().is_empty()
}

fn prepare_reference(effective: &EffectiveConfig) -> Result<PreparedReference, Diagnostic> {
    let prepared = prepare_reference_without_repeat_admission(effective)?;
    admit_declared_repeat_loci(effective, prepared.reference())?;
    Ok(prepared)
}

/// **Guarda de admisión de los loci sobre la referencia.**
///
/// Al declarar un locus sobre una referencia **suministrada o reducida** hay que
/// comprobar que en esa posición hay realmente un tracto de los motivos
/// declarados. Sin ella se anotaría como microsatélite un trozo cualquiera del
/// genoma, y todo lo que viniera después —los alelos, la verdad, el VCF— hablaría
/// de una repetición que no está ahí.
///
/// El criterio es el mismo par de mínimos que ya se aplica al alelo declarado:
/// **pureza mínima** y **unidades mínimas en la tirada más larga**, los del locus
/// si los declara y si no los generales. Los dos hacen falta y ninguno sobra: un
/// tramo puede ser bastante puro y no tener ni una tirada aprovechable, y al
/// revés.
///
/// **Con referencia simulada también se aplica, desde el 01/10/2026**, sobre el
/// tracto que la preparación acaba de colocar. Hasta ese día se saltaba con el
/// argumento de la decisión S-06, «el tracto no existe hasta que se crea»; pero
/// desde que la preparación lo coloca (S-03) sí existe antes de simular, y
/// saltarse la guarda publicaba tractos que el programa rechaza en cualquier
/// otro sitio: `ACG` sobre ocho bases coloca `ACGACGAC`, de pureza 0,75, que la
/// misma referencia suministrada no pasaría. Un tracto de un motivo sin copia
/// parcial, o uno compuesto de bloques enteros, tiene pureza 1 y la tirada que
/// declara: solo se rechaza lo que de verdad no cumple.
fn admit_declared_repeat_loci(
    effective: &EffectiveConfig,
    reference: &Reference,
) -> Result<(), Diagnostic> {
    if effective.repeat_loci().is_empty() {
        return Ok(());
    }
    for locus in effective.repeat_loci() {
        let interval = locus.declaration().interval();
        let sequence = reference
            .sequences()
            .iter()
            .find(|sequence| sequence.id() == interval.sequence_id())
            .ok_or_else(|| {
                repeat_admission_error(
                    locus.id(),
                    format!(
                        "la referencia no contiene la secuencia {}",
                        interval.sequence_id()
                    ),
                )
            })?;
        let bases = sequence.sequence().as_bytes().map_err(|error| {
            internal(
                "admit_declared_repeat_loci",
                format!("las bases de la referencia no están materializadas: {error:?}"),
            )
        })?;
        let (start, end) = (interval.start() as usize, interval.end() as usize);
        if end > bases.len() {
            return Err(repeat_admission_error(
                locus.id(),
                format!(
                    "el intervalo declarado rebasa las {} bases de {}",
                    bases.len(),
                    interval.sequence_id()
                ),
            ));
        }
        let observed = observe_tract(&bases[start..end], locus.motifs());
        let admission = locus.admission();
        if observed.purity_as_f64() < admission.min_purity() {
            return Err(repeat_admission_error(
                locus.id(),
                format!(
                    "el tramo declarado tiene pureza {} y el mínimo es {}",
                    observed.purity_as_f64(),
                    admission.min_purity()
                ),
            ));
        }
        if observed.longest_run_units() < admission.min_longest_run_units() {
            return Err(repeat_admission_error(
                locus.id(),
                format!(
                    "el tramo declarado tiene una tirada más larga de {} unidades y el mínimo es {}",
                    observed.longest_run_units(),
                    admission.min_longest_run_units()
                ),
            ));
        }
    }
    Ok(())
}

/// **Un cambio pequeño dentro del intervalo de un locus rechaza el escenario**
/// (decisión S-32), y se rechaza **antes de simular**.
///
/// # Por qué tiene que ser aquí y no en el estado
///
/// El alelo de un microsatélite describe el tracto **tal y como está en la
/// referencia**. Si otra modificación ya tocó ese tramo, la anotación publicada
/// describiría un tracto que ya no existe.
///
/// La cadena de segmentos tiene su propia guarda, pero **no basta por sí
/// sola**: solo se consulta cuando hay un alelo que aplicar, y el alelo se
/// sortea. Con el alelo de la referencia entre los declarados, que el escenario
/// se aceptara o no acababa dependiendo de la semilla, y cuando se aceptaba la
/// verdad publicaba que esa copia «conserva el tracto de la referencia» cuando
/// no lo conserva. Por eso la decisión vive aquí, donde no hay sorteo ninguno,
/// y la del estado queda como red interna.
///
/// # Qué se mira aquí, y por qué **no** las inserciones ni las deleciones
///
/// Se miran los **cambios de base** declarados y **cada fila de la matriz de
/// variantes**, venga suministrada o generada. Las dos cosas declaran su
/// posición en **coordenadas de la referencia**, que es el mismo espacio en el
/// que vive el intervalo del locus, de modo que la comparación es exacta: ni
/// deja pasar una invasión ni inventa una que no existe. Las generadas eligen
/// su posición con un sorteo propio, así que sin esta comprobación declarar un
/// microsatélite y pedir variantes generadas salía bien o mal según la suerte.
///
/// **Una inserción o una deleción no se juzgan aquí**, y es deliberado: su
/// punto se declara en **coordenadas vigentes** —las de la secuencia después de
/// los cambios anteriores de esa misma copia—, no en las de la referencia. Con
/// un solo indel los dos espacios coinciden; encadenando varios se separan, y
/// entonces compararlos contra el intervalo **se equivoca en las dos
/// direcciones**: deja escapar una invasión real y rechaza un escenario
/// legítimo diciendo que algo cae dentro del tracto cuando cae fuera.
///
/// Las dos direcciones están comprobadas con pruebas propias, y por eso este
/// código no las mira. Quien sí las juzga con exactitud es la **cadena de
/// segmentos**: allí el intervalo se busca en coordenadas de referencia sobre
/// el estado ya compuesto, de modo que no hay deriva posible. Esa comprobación
/// se consulta **para toda copia y todo locus**, se aplique o no un alelo, y
/// por tanto tampoco depende de la semilla.
///
/// La alternativa era traducir aquí el punto declarado a coordenadas de la
/// referencia descontando los indeles anteriores. Sería un **segundo modelo de
/// coordenadas** junto al que la cadena de segmentos ya sostiene, y este
/// proyecto ya ha pagado caro tener dos implementaciones del mismo concepto.
fn admit_repeat_loci_against_small_variants(
    effective: &EffectiveConfig,
    matrix: Option<&crate::population_snvs::PreparedPopulationSnvs>,
) -> Result<(), Diagnostic> {
    if effective.repeat_loci().is_empty() {
        return Ok(());
    }
    for locus in effective.repeat_loci() {
        let interval = locus.declaration().interval();
        let sequence_id = interval.sequence_id().as_str();
        let (start, end) = (interval.start(), interval.end());
        let invades = |other_sequence: &str, other_start: u64, other_end: u64| {
            other_sequence == sequence_id && other_start < end && start < other_end
        };
        for event in effective.direct_events() {
            // **El identificador del evento que coloca un alelo es el del
            // propio locus**, de modo que un cambio directo llamado igual
            // haría ambiguo el origen de la modificación en la verdad. Se
            // rechaza aquí, sobre lo que el usuario escribió, en vez de morir
            // al final con un error de formato sobre un documento producido.
            if event.event_id().as_str() == locus.id() {
                return Err(repeat_admission_error(
                    locus.id(),
                    format!(
                        "el cambio directo {} lleva el mismo identificador que el locus, y el locus lo usa como identificador del evento que coloca su alelo",
                        event.event_id()
                    ),
                ));
            }
            // Solo el cambio de base declara su posición en coordenadas de la
            // referencia. La inserción y la deleción la declaran en las
            // vigentes, y compararlas aquí se equivoca en las dos direcciones:
            // las juzga la cadena de segmentos, que no tiene deriva posible.
            let EffectiveDirectEvent::Snp {
                sequence_id: event_sequence,
                position_1based,
                ..
            } = event
            else {
                continue;
            };
            if invades(
                event_sequence.as_str(),
                position_1based.saturating_sub(1),
                *position_1based,
            ) {
                return Err(repeat_invasion_error(
                    locus,
                    &format!(
                        "el cambio de base {} en {event_sequence}:{position_1based}",
                        event.event_id()
                    ),
                ));
            }
        }
        for row in matrix.map(|matrix| matrix.rows()).unwrap_or_default() {
            if invades(
                row.contig().as_str(),
                row.position_1based().saturating_sub(1),
                row.position_1based(),
            ) {
                return Err(repeat_invasion_error(
                    locus,
                    &format!(
                        "la variante {} de la matriz en {}:{}",
                        row.variant_id(),
                        row.contig(),
                        row.position_1based()
                    ),
                ));
            }
        }
    }
    Ok(())
}

/// El diagnóstico que el usuario va a leer cuando invade un microsatélite.
///
/// Dice **qué** lo invade y **qué hacer en su lugar**, porque la decisión S-32
/// deja al usuario sin la vía obvia: una variante puntual dentro de un
/// microsatélite hay que expresarla como parte del alelo del locus.
fn repeat_invasion_error(
    locus: &crate::application::EffectiveRepeatLocus,
    invader: &str,
) -> Diagnostic {
    let interval = locus.declaration().interval();
    repeat_admission_error(
        locus.id(),
        format!(
            "{invader} cae dentro del tracto declarado ({}:{} a {}). Un alelo de repetición describe el tracto que hay en la referencia, de modo que ninguna otra modificación puede tocarlo: exprésalo como parte de un alelo del locus, con su secuencia o con sus interrupciones",
            interval.sequence_id(),
            interval.start() + 1,
            interval.end()
        ),
    )
}

fn repeat_admission_error(locus_id: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", "effective-config")
    .with("fields", "reference.repeat_loci")
    .with("detail", format!("locus {locus_id}: {}", detail.into()))
}

fn prepare_reference_without_repeat_admission(
    effective: &EffectiveConfig,
) -> Result<PreparedReference, Diagnostic> {
    match effective.reference_source() {
        // Una ejecución que solo deriva referencias reducidas no necesita las
        // bases del genoma de partida: planifica con el catálogo y pide después
        // los tramos del plan. Validarlo sin retenerlas mantiene la memoria
        // independiente del tamaño de ese genoma.
        EffectiveReferenceSource::Supplied { original_path }
            if !requires_materialized_bases(effective) =>
        {
            let summary = validate_supplied_fasta_summary(
                original_path,
                effective.reference_id(),
                &sequence_expectations(effective),
            )?;
            let reference = reference_from_summary(
                &summary,
                effective.reference_id(),
                crate::reference::ReferenceOrigin::Supplied {
                    original_path: original_path.clone(),
                },
            )?;
            Ok(PreparedReference::Supplied(
                ValidatedSuppliedReference::from_parts(
                    reference,
                    summary.source_identity().clone(),
                ),
            ))
        }
        EffectiveReferenceSource::Supplied { original_path } => {
            let supplied = validate_supplied_fasta(
                original_path,
                effective.reference_id(),
                &sequence_expectations(effective),
            )?;
            if effective.tract_insertion().is_some() {
                insert_declared_tracts(effective, supplied)
            } else {
                Ok(PreparedReference::Supplied(supplied))
            }
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
            let tracts = simulated_repeat_tracts(effective)?;
            let regions = tracts
                .iter()
                .map(|tract| DeclaredRegion {
                    sequence_id: tract.sequence_id,
                    start: tract.start,
                    bases: &tract.bases,
                })
                .collect::<Vec<_>>();
            generate_reference_with_regions(
                effective.reference_id(),
                effective.reproducibility().root_seed(),
                &specs,
                &frequencies,
                composition.unknown_base_fraction(),
                &regions,
            )
            .map(PreparedReference::Generated)
        }
    }
}

/// **Construye la referencia final insertando los tractos declarados**
/// (REQ-REF-37, DTC-098).
///
/// Cada punto `insert_after` está en base uno del genoma de partida y dice
/// tras qué base entra el tracto; el módulo de inserción recibe el punto
/// interbase en base cero, que es el mismo número: «antes de la base cero `P`»
/// es «después de la base uno `P`». Todas las posiciones se leen en el genoma
/// de partida, también la segunda y las siguientes.
///
/// Una lista vacía se rechaza aquí: el módulo la aceptaría y produciría una
/// «derivada» idéntica a la suministrada con otro nombre. La resolución de la
/// configuración no crea un plan sin puntos, así que esta rama no se alcanza
/// desde una configuración; se conserva porque es la frontera del módulo.
///
/// Después de construirla se sueltan las bases de la suministrada: lo que
/// queda de ella es su catálogo y su identidad, que es lo único que hace falta
/// para publicar su copia y describirla.
fn insert_declared_tracts(
    effective: &EffectiveConfig,
    supplied: ValidatedSuppliedReference,
) -> Result<PreparedReference, Diagnostic> {
    let plan = effective
        .tract_insertion()
        .ok_or_else(|| internal("insert_declared_tracts", "falta el plan de inserción"))?;
    let insertions = plan
        .points()
        .iter()
        .map(|point| DeclaredInsertion {
            sequence_id: point.sequence_id(),
            source_offset: point.insert_after(),
            bases: point.bases(),
        })
        .collect::<Vec<_>>();
    if insertions.is_empty() {
        return Err(internal(
            "insert_declared_tracts",
            "el plan de inserción no tiene ningún punto",
        ));
    }
    let (derived, map) = insert_tracts(
        supplied.reference(),
        plan.derived_reference_id(),
        &insertions,
    )?
    .into_parts();
    let catalog = reference_without_bases(supplied.reference())?;
    Ok(PreparedReference::Inserted {
        supplied: ValidatedSuppliedReference::from_parts(
            catalog,
            supplied.source_identity().clone(),
        ),
        derived,
        map,
    })
}

/// La misma referencia, con su catálogo y sin sus bases.
fn reference_without_bases(reference: &Reference) -> Result<Reference, Diagnostic> {
    let sequences = reference
        .sequences()
        .iter()
        .map(|sequence| {
            let summary = crate::reference::DnaSequence::summary(
                sequence.length(),
                sequence.unknown_base_count(),
            )
            .map_err(|error| internal("reference_without_bases", error))?;
            crate::reference::ReferenceSequence::from_owned(
                sequence.id().clone(),
                sequence.header().clone(),
                sequence.class(),
                summary,
            )
            .map_err(|error| internal("reference_without_bases", error))
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;
    Reference::new(
        reference.id().clone(),
        reference.origin().clone(),
        sequences,
    )
    .map_err(|error| internal("reference_without_bases", error))
}

/// FASTA y mapa de la referencia final, ya publicados.
struct PublishedInsertion {
    fasta: ArtifactRecord,
    map: ArtifactRecord,
}

/// Publica el FASTA canónico de la referencia final y su mapa de coordenadas
/// (DTC-098), con la misma forma de ruta que una reducción.
fn write_inserted_reference(
    prepared: &PreparedReference,
    workspace: &mut ExecutionWorkspace,
) -> Result<Option<PublishedInsertion>, Diagnostic> {
    let PreparedReference::Inserted { derived, map, .. } = prepared else {
        return Ok(None);
    };
    let id = derived.id().as_str();
    let fasta = workspace.write_streamed_artifact(
        &format!("reference-fasta-{id}"),
        "reference_fasta",
        &PathBuf::from(format!("resources/references/{id}.fa")),
        |file| write_generated_fasta(derived, file),
    )?;
    let map_bytes = json_document_bytes(&coordinate_map_document(map)?)?;
    let map = workspace.write_artifact(
        &format!("reference-map-{id}"),
        COORDINATE_MAP_DOCUMENT_KIND,
        &PathBuf::from(format!("resources/references/{id}.map.json")),
        &map_bytes,
    )?;
    Ok(Some(PublishedInsertion { fasta, map }))
}

/// El tracto que un locus coloca en una referencia simulada, en base cero.
struct SimulatedRepeatTract<'a> {
    sequence_id: &'a str,
    start: u64,
    bases: Vec<u8>,
}

/// **Los tractos que la referencia simulada lleva antes del relleno al azar**
/// (decisión S-03, con la forma que fijó el responsable el 01/10/2026).
///
/// En un locus de **un solo motivo**, el tramo declarado es el tracto mismo
/// —de `start` a `end` inclusive, sin flancos: es el tramo que la guarda de
/// admisión mide sobre una referencia real y del que el generador de alelos
/// lee el ancestral— y se rellena con el motivo repetido desde su primera
/// base. Si el tramo no es múltiplo del motivo, la última copia queda parcial.
/// El tipo del locus garantiza que cabe al menos una copia entera.
///
/// En un locus **compuesto**, cada bloque se coloca con su motivo repetido
/// tantas veces como diga `reference_units`, en el orden de los motivos
/// (DT-016). La resolución de la configuración ya comprobó que el campo existe
/// y que su suma es exactamente la longitud del tramo.
///
/// Hasta el 01/10/2026 esta función no existía y la preparación llamaba al
/// generador **sin regiones**: el FASTA publicado llevaba bases al azar en el
/// tramo de cada locus y ningún alelo declarado coincidía con la referencia
/// (DT-013). Sin loci la lista sale vacía, y con la lista vacía el generador
/// consume y escribe exactamente lo mismo que antes.
fn simulated_repeat_tracts(
    effective: &EffectiveConfig,
) -> Result<Vec<SimulatedRepeatTract<'_>>, Diagnostic> {
    let mut tracts = Vec::with_capacity(effective.repeat_loci().len());
    for locus in effective.repeat_loci() {
        let interval = locus.declaration().interval();
        let bases = match (locus.single_motif(), locus.reference_units()) {
            (Some(motif), _) => {
                let length = usize::try_from(interval.len()).map_err(|_| {
                    internal(
                        "simulated_repeat_tracts",
                        format!("el tramo del locus {} no cabe en memoria", locus.id()),
                    )
                })?;
                motif.bases().iter().copied().cycle().take(length).collect()
            }
            (None, Some(units)) => {
                let blocks = locus
                    .motifs()
                    .iter()
                    .cloned()
                    .zip(units.iter().copied())
                    .collect::<Vec<_>>();
                TandemRepeatTract::pure_blocks(&blocks)
                    .map_err(|error| {
                        internal(
                            "simulated_repeat_tracts",
                            format!("el tracto del locus {}: {}", locus.id(), error.detail()),
                        )
                    })?
                    .bases()
                    .to_vec()
            }
            // La resolución de la configuración exige `reference_units` a todo
            // locus compuesto sobre referencia simulada, de modo que aquí no se
            // llega desde una configuración. Si se llegara, colocar nada
            // dejaría el tramo al azar en silencio, que es justo DT-016.
            (None, None) => {
                return Err(internal(
                    "simulated_repeat_tracts",
                    format!(
                        "el locus compuesto {} llega sin unidades por bloque",
                        locus.id()
                    ),
                ));
            }
        };
        tracts.push(SimulatedRepeatTract {
            sequence_id: interval.sequence_id().as_str(),
            start: interval.start(),
            bases,
        });
    }
    Ok(tracts)
}

fn prepare_reductions(
    effective: &EffectiveConfig,
    source: &Reference,
    source_fasta: &Path,
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
        let source_bases = source_for(source, Some(source_fasta))?;
        let reference = build_reduced_reference(source_bases.as_ref(), &plan)?;
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
    inserted_artifact: Option<&ArtifactRecord>,
    prepared: PreparedReference,
) -> Result<PreparedReference, Diagnostic> {
    // Ningún brazo necesita la referencia ya preparada, solo su origen. Se
    // libera antes de releer la copia publicada porque, de lo contrario, la
    // preparada y la autoritativa conviven y el pico dobla el genoma.
    let prepared_was_supplied = matches!(prepared, PreparedReference::Supplied(_));
    // Consumir la preparada aquí la suelta antes de releer: con inserción solo
    // sobrevive su mapa, que no lleva bases.
    let inserted_map = match prepared {
        PreparedReference::Inserted { map, .. } => Some(map),
        _ => None,
    };
    if let Some(map) = inserted_map {
        return load_authoritative_inserted_reference(
            effective,
            workspace,
            artifact,
            inserted_artifact.ok_or_else(|| {
                internal(
                    "load_authoritative_reference",
                    "la referencia final no está publicada",
                )
            })?,
            map,
        );
    }
    let mut authoritative = match (effective.reference_source(), prepared_was_supplied) {
        (EffectiveReferenceSource::Supplied { original_path }, true) => {
            let physical_path = workspace.root().join(&artifact.path);
            // La copia publicada se relee para comprobar que coincide con lo
            // catalogado. Cuando la configuración no necesita las bases basta el
            // catálogo: releerlas solo para descartarlas volvería a atar la
            // memoria al tamaño del genoma.
            let supplied = if requires_materialized_bases(effective) {
                validate_supplied_fasta_with_origin(
                    &physical_path,
                    original_path,
                    effective.reference_id(),
                    &sequence_expectations(effective),
                )?
            } else {
                let summary = validate_supplied_fasta_summary(
                    &physical_path,
                    effective.reference_id(),
                    &sequence_expectations(effective),
                )?;
                let reference = reference_from_summary(
                    &summary,
                    effective.reference_id(),
                    crate::reference::ReferenceOrigin::Supplied {
                        original_path: original_path.clone(),
                    },
                )?;
                ValidatedSuppliedReference::from_parts(reference, summary.source_identity().clone())
            };
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
            // Se consumen las secuencias en lugar de clonarlas: clonarlas
            // mantendría vivas dos copias del genoma a la vez.
            let reference = Reference::new(
                ReferenceId::parse(effective.reference_id(), "execution")?,
                crate::reference::ReferenceOrigin::Generated,
                parsed.into_reference().into_sequences(),
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

/// Relee las dos copias publicadas de una ejecución con inserción (DTC-098).
///
/// La suministrada se relee **sin bases**: basta para comprobar que la copia
/// es la catalogada y para describirla. La final se relee entera desde su FASTA
/// publicado, que desde aquí es la única autoridad de la simulación; su huella
/// se compara con la del catálogo, igual que la de cualquier copia.
fn load_authoritative_inserted_reference(
    effective: &EffectiveConfig,
    workspace: &ExecutionWorkspace,
    supplied_artifact: &ArtifactRecord,
    inserted_artifact: &ArtifactRecord,
    map: CoordinateMap,
) -> Result<PreparedReference, Diagnostic> {
    let EffectiveReferenceSource::Supplied { original_path } = effective.reference_source() else {
        return Err(internal(
            "load_authoritative_inserted_reference",
            "origin changed",
        ));
    };
    let plan = effective.tract_insertion().ok_or_else(|| {
        internal(
            "load_authoritative_inserted_reference",
            "falta el plan de inserción",
        )
    })?;
    let supplied_path = workspace.root().join(&supplied_artifact.path);
    let summary = validate_supplied_fasta_summary(
        &supplied_path,
        effective.reference_id(),
        &sequence_expectations(effective),
    )?;
    if summary.source_identity().size_bytes() != supplied_artifact.size_bytes
        || summary.source_identity().digest() != supplied_artifact.digest
    {
        return Err(internal(
            "load_authoritative_inserted_reference",
            "cataloged reference changed",
        ));
    }
    let mut catalog = reference_from_summary(
        &summary,
        effective.reference_id(),
        crate::reference::ReferenceOrigin::Supplied {
            original_path: original_path.clone(),
        },
    )?;
    catalog
        .bind_authoritative_resource(PathBuf::from(&supplied_artifact.path))
        .map_err(|error| internal("bind_authoritative_reference", error))?;
    let derived_path = workspace.root().join(&inserted_artifact.path);
    let expectations = plan
        .sequences()
        .iter()
        .map(|sequence| SequenceExpectation {
            id: sequence.id(),
            class: sequence.class(),
            length: sequence.length(),
        })
        .collect::<Vec<_>>();
    let parsed =
        validate_supplied_fasta(&derived_path, plan.derived_reference_id(), &expectations)?;
    if parsed.source_identity().size_bytes() != inserted_artifact.size_bytes
        || parsed.source_identity().digest() != inserted_artifact.digest
    {
        return Err(internal(
            "load_authoritative_inserted_reference",
            "cataloged inserted reference changed",
        ));
    }
    let mut derived = Reference::new(
        ReferenceId::parse(plan.derived_reference_id(), "execution")?,
        crate::reference::ReferenceOrigin::TractInsertion {
            source_reference_id: ReferenceId::parse(effective.reference_id(), "execution")?,
        },
        parsed.into_reference().into_sequences(),
    )
    .map_err(|error| internal("reload_inserted_reference", error))?;
    derived
        .bind_authoritative_resource(PathBuf::from(&inserted_artifact.path))
        .map_err(|error| internal("bind_authoritative_reference", error))?;
    Ok(PreparedReference::Inserted {
        supplied: ValidatedSuppliedReference::from_parts(
            catalog,
            summary.source_identity().clone(),
        ),
        derived,
        map,
    })
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
            PreparedReference::Supplied(supplied) | PreparedReference::Inserted { supplied, .. },
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

    /// **Un paquete fallido con lecturas conserva la referencia final con su
    /// mapa y su sección** (DTC-098, REQ-MAN-07). El fallo se provoca justo
    /// antes del manifiesto, cuando las lecturas ya salieron de la final: el
    /// paquete fallido las conserva, y con ellas la final; si conservara el
    /// FASTA sin el mapa ni `reference_insertion`, el manifiesto nombraría una
    /// referencia derivada sin decir de dónde deriva.
    #[test]
    fn a_failed_run_with_reads_keeps_the_final_reference_with_its_map() {
        let root = std::env::temp_dir().join(format!(
            "gps-insertion-failure-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("test root");
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let bases = (0..2000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                char::from(b"ACGT"[(state % 4) as usize])
            })
            .collect::<String>();
        fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
        let config = root.join("general.toml");
        fs::write(
            &config,
            r#"config_kind = "general"
schema_version = "4"
seed = "00000000000000000000000000000000000000000000000000000000000000aa"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
outside_weight = 0.0
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.regions]]
contig = "chr1"
start = 1500
end = 1600
weight = 1.0

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"
inserted_reference_id = "ref-main-str"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 2000

[[reference.repeat_loci]]
id = "STR-INS"
sequence = "chr1"
insert_after = 1000
motif = "CAG"
reference_units = [10]
[[reference.repeat_loci.alleles]]
units = 10
frequency = 1.0

[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["I1"]
ploidy = 2
genotype_model = "hwe"

[technical_event]
event_id = "reads-main"
individual_ids = ["I1"]
target_coverage_wgs = 2
sequencing_errors = "disabled"
[[technical_event.libraries]]
individual_id = "I1"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
        )
        .expect("configuración");
        let effective = load_and_validate_config(&config).expect("configuración válida");
        let destination = root.join("result");
        let error = execute_validated_config_with_effective_document(
            &effective,
            &destination,
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            |_| {
                Err(internal(
                    "test_failure",
                    "fallo provocado antes del manifiesto",
                ))
            },
        )
        .expect_err("la ejecución falla");
        assert!(format!("{error:?}").contains("fallo provocado"));
        let manifest: serde_json::Value = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("manifiesto fallido"),
        )
        .expect("JSON");
        assert_eq!(manifest["state"], "failed");
        let insertion = &manifest["reference_insertion"];
        assert_eq!(insertion["state"], "failed");
        assert_eq!(insertion["id"], "ref-main-str");
        let kinds = manifest["artifacts"]
            .as_array()
            .expect("artefactos")
            .iter()
            .map(|artifact| artifact["id"].as_str().expect("id").to_owned())
            .collect::<Vec<_>>();
        assert!(
            kinds.contains(&"reference-fasta-ref-main-str".to_owned()),
            "{kinds:?}"
        );
        assert!(
            kinds.contains(&"reference-map-ref-main-str".to_owned()),
            "{kinds:?}"
        );
        assert!(
            destination
                .join("resources/references/ref-main-str.map.json")
                .exists()
        );
        fs::remove_dir_all(root).ok();
    }

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
                // La referencia de un escenario que solo deriva reducciones se
                // valida sin conservar las bases, así que el objeto obsoleto se
                // fabrica alterando la ficha: si la selección lo usara en lugar
                // de la copia recargada, el recuento falseado llegaría a la
                // configuración efectiva y su validación física lo rechazaría.
                let dna = first
                    .sequence()
                    .as_bytes()
                    .map_or_else(
                        |_| DnaSequence::summary(first.length(), first.unknown_base_count() + 1),
                        |bases| {
                            let stale = bases
                                .iter()
                                .map(|base| if *base == b'N' { b'N' } else { b'A' })
                                .collect::<Vec<_>>();
                            DnaSequence::try_from_bytes(&stale)
                        },
                    )
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

    const TRUTH_JUDGE_CONFIG: &[u8] = b"\n[population_snvs]\nmode = \"generate\"\n[population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 1\nspatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n[[population_snvs.generate.populations]]\npopulation_id = \"P1\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nindividual_ids = [\"IND_01\"]\nploidy = 2\ngenotype_model = \"hwe\"\n";

    /// Ejecuta de extremo a extremo, con el productor de la verdad alterado,
    /// una configuración con un caso genotípico, con o sin `[vcf_output]`.
    /// Devuelve el diagnóstico con que acaba.
    fn run_with_tampered_producer(
        with_vcf: bool,
        tamper: impl Fn(&mut crate::truth_schema::TruthDocument) + 'static,
    ) -> Diagnostic {
        let root = std::env::temp_dir().join(format!(
            "gps-dt019-{}-{}",
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
        bytes.extend_from_slice(TRUTH_JUDGE_CONFIG);
        if with_vcf {
            bytes.extend_from_slice(b"[vcf_output]\nmode = \"bgzf\"\n");
        }
        fs::write(&config, bytes).expect("config");
        let effective = load_and_validate_config(&config).expect("valid config");
        let _guard = crate::truth_projection::test_tamper::install(tamper);
        let error = execute_validated_config_with_effective_document(
            &effective,
            &destination,
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            |_| Ok(()),
        )
        .expect_err("a tampered producer must not publish");
        fs::remove_dir_all(root).expect("cleanup");
        error
    }

    fn shift_first_genotype_locus(document: &mut crate::truth_schema::TruthDocument) {
        for record in &mut document.records {
            if let crate::truth_schema::TruthRecord::Genotype(genotype) = record {
                genotype.locus.start += 1;
                genotype.locus.end += 1;
                return;
            }
        }
        panic!("the scenario must publish a genotype case");
    }

    /// Intercambia, de forma coherente, qué alelo lleva cada copia y la dosis
    /// de cada alelo: la estructura del documento sigue siendo válida y sus
    /// totales no cambian, de modo que solo el juez semántico puede verlo.
    fn change_first_genotype_dosage(document: &mut crate::truth_schema::TruthDocument) {
        for record in &mut document.records {
            if let crate::truth_schema::TruthRecord::Genotype(genotype) = record {
                assert_eq!(genotype.allele_dosages.len(), 2, "REF and one ALT");
                genotype.allele_dosages.swap(0, 1);
                for carrier in &mut genotype.carriers {
                    carrier.allele_index = 1 - carrier.allele_index;
                }
                return;
            }
        }
        panic!("the scenario must publish a genotype case");
    }

    /// DT-019. Sin `[vcf_output]` —la configuración por omisión— nadie más que
    /// el juez de la publicación contrasta los casos genotípicos: un locus
    /// desplazado tiene que acabar rechazado como defecto interno.
    #[test]
    fn shifted_case_locus_without_vcf_is_rejected_as_internal_defect() {
        let error = run_with_tampered_producer(false, shift_first_genotype_locus);
        assert_eq!(error.code(), "GPS-INTERNAL-001");
        assert!(
            format!("{error:?}").contains("GenotypeJudgement"),
            "the judge, not another check, must reject: {error:?}"
        );
    }

    /// DT-019, con la dosis y los portadores. Si el escenario fijado dejara
    /// de producir un caso con REF y un ALT, el `assert_eq!` de la ayuda lo
    /// diría en vez de fingir que discrimina. Intercambiar los índices de
    /// alelo de los portadores siempre cambia quién lleva qué, también con
    /// dosis simétricas.
    #[test]
    fn changed_case_dosage_without_vcf_is_rejected_as_internal_defect() {
        let error = run_with_tampered_producer(false, change_first_genotype_dosage);
        assert_eq!(error.code(), "GPS-INTERNAL-001");
        assert!(
            format!("{error:?}").contains("GenotypeJudgement"),
            "the judge, not another check, must reject: {error:?}"
        );
    }

    /// Con VCF el productor alterado también se rechaza, ya en la proyección y
    /// como defecto interno, antes de escribir ningún VCF.
    #[test]
    fn shifted_case_locus_with_vcf_is_rejected_before_any_vcf() {
        let error = run_with_tampered_producer(true, shift_first_genotype_locus);
        assert_eq!(error.code(), "GPS-INTERNAL-001");
        assert!(
            format!("{error:?}").contains("GenotypeJudgement"),
            "{error:?}"
        );
    }

    // ---- DT-019: un escenario con varios casos de cada clase, sin VCF ----

    use crate::truth_schema::{
        GenotypeCaseKind, GenotypeRecord, TruthAllele, TruthDocument, TruthRecord,
        TruthStructuralEvent,
    };

    /// Ejecuta sin VCF el escenario de efecto neto nulo con repeticiones,
    /// ploidía 4 y dos individuos, al que se añaden una inserción y una
    /// deleción más, con el gancho de alteración instalado.
    fn run_direct_cancel_without_vcf(
        tamper: impl Fn(&mut TruthDocument) + 'static,
    ) -> Result<(), Diagnostic> {
        let root = std::env::temp_dir().join(format!(
            "gps-dt019-multi-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("root");
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/vcf-null-effect/direct-cancel");
        for name in ["reference.fa", "snvs.tsv"] {
            fs::copy(fixture.join(name), root.join(name)).expect("fixture file");
        }
        let mut config = fs::read_to_string(fixture.join("general.toml")).expect("config");
        config = config.replace("[vcf_output]\nmode = \"text\"\n", "");
        assert!(
            !config.contains("vcf_output"),
            "the scenario must have no VCF"
        );
        config.push_str(
            "\n[[direct_events]]\nevent_id = \"extra-ins\"\ncontig = \"chr1\"\nindividual_id = \"ind-1\"\nphase = 0\nkind = \"insertion\"\npoint = 60\nbases = \"GGA\"\n\n[[direct_events]]\nevent_id = \"extra-del\"\ncontig = \"chr1\"\nindividual_id = \"ind-1\"\nphase = 2\nkind = \"deletion\"\nstart = 250\nend = 253\n",
        );
        let path = root.join("general.toml");
        fs::write(&path, config).expect("write config");
        let effective = load_and_validate_config(&path).expect("valid config");
        let _guard = crate::truth_projection::test_tamper::install(tamper);
        let result = execute_validated_config_with_effective_document(
            &effective,
            &root.join("out"),
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            |_| Ok(()),
        )
        .map(|_| ());
        fs::remove_dir_all(root).expect("cleanup");
        result
    }

    fn assert_judge_rejects(result: Result<(), Diagnostic>) {
        let error = result.expect_err("the tampered producer must not publish");
        assert_eq!(error.code(), "GPS-INTERNAL-001", "{error:?}");
        assert!(
            format!("{error:?}").contains("GenotypeJudgement"),
            "the judge, not another check, must reject: {error:?}"
        );
    }

    fn genotypes_mut(document: &mut TruthDocument) -> Vec<&mut GenotypeRecord> {
        document
            .records
            .iter_mut()
            .filter_map(|record| match record {
                TruthRecord::Genotype(genotype) => Some(genotype),
                TruthRecord::HaplotypeCopy(_) => None,
            })
            .collect()
    }

    /// Reordena los casos en el orden canónico, para que un locus alterado no
    /// lo cace antes la comprobación de orden del contrato.
    fn sort_genotypes(document: &mut TruthDocument) {
        let copies = document
            .records
            .iter()
            .filter(|record| matches!(record, TruthRecord::HaplotypeCopy(_)))
            .count();
        document.records[copies..].sort_by(|left, right| match (left, right) {
            (TruthRecord::Genotype(a), TruthRecord::Genotype(b)) => (
                &a.reference_id,
                &a.individual_id,
                &a.locus,
                &a.case_id,
            )
                .cmp(&(&b.reference_id, &b.individual_id, &b.locus, &b.case_id)),
            _ => unreachable!("genotypes follow copies"),
        });
    }

    fn swap_two_differing_carriers(genotype: &mut GenotypeRecord) -> bool {
        for i in 0..genotype.carriers.len() {
            for j in i + 1..genotype.carriers.len() {
                if genotype.carriers[i].allele_index != genotype.carriers[j].allele_index {
                    let first = genotype.carriers[i].allele_index;
                    genotype.carriers[i].allele_index = genotype.carriers[j].allele_index;
                    genotype.carriers[j].allele_index = first;
                    return true;
                }
            }
        }
        false
    }

    /// Un locus de indel desplazado. Es lo que DT-019 dejó pasar de verdad.
    /// Lo cazan solo el juez de los casos indel: ignorarlos lo deja pasar.
    #[test]
    fn shifted_indel_locus_without_vcf_is_rejected() {
        assert_judge_rejects(run_direct_cancel_without_vcf(|document| {
            let genotype = genotypes_mut(document)
                .into_iter()
                .find(|genotype| genotype.case_kind == GenotypeCaseKind::Indel)
                .expect("an indel case");
            genotype.locus.start += 1;
            genotype.locus.end += 1;
            sort_genotypes(document);
        }));
    }

    /// Las bases de una inserción cambiadas **junto con** las de su evento
    /// fuente: el documento es coherente consigo mismo y solo contrastarlo con
    /// la copia reconstruida lo delata.
    #[test]
    fn changed_insertion_bases_with_their_source_event_without_vcf_are_rejected() {
        assert_judge_rejects(run_direct_cancel_without_vcf(|document| {
            let mut changed = false;
            'outer: for genotype in genotypes_mut(document) {
                if genotype.case_kind != GenotypeCaseKind::Indel {
                    continue;
                }
                for allele in &mut genotype.alleles {
                    if let TruthAllele::Insertion {
                        bases,
                        source_event: TruthStructuralEvent::Insertion { bases: source, .. },
                        ..
                    } = allele
                    {
                        let first = if bases.starts_with('T') { 'A' } else { 'T' };
                        bases.replace_range(..1, &first.to_string());
                        source.clone_from(bases);
                        changed = true;
                        break 'outer;
                    }
                }
            }
            assert!(changed, "the scenario must carry an insertion case");
        }));
    }

    /// Los portadores de un caso de repetición intercambiados: quién lleva qué
    /// alelo de la repetición. Ignorar los casos de repetición lo deja pasar.
    #[test]
    fn swapped_repeat_carriers_without_vcf_are_rejected() {
        assert_judge_rejects(run_direct_cancel_without_vcf(|document| {
            let swapped = genotypes_mut(document)
                .into_iter()
                .filter(|genotype| genotype.case_kind == GenotypeCaseKind::Repeat)
                .any(swap_two_differing_carriers);
            assert!(swapped, "a heterozygous repeat case must exist");
        }));
    }

    /// Una alteración en el **último** caso del documento, que no es el
    /// primero: juzgar solo el primer caso, o parar en él, lo deja pasar.
    #[test]
    fn tampering_only_the_last_case_without_vcf_is_rejected() {
        assert_judge_rejects(run_direct_cancel_without_vcf(|document| {
            let mut genotypes = genotypes_mut(document);
            assert!(genotypes.len() > 2, "the scenario must carry several cases");
            let last = genotypes.pop().expect("last case");
            last.locus.start += 1;
            last.locus.end += 1;
            let altered = last.case_id.clone();
            sort_genotypes(document);
            let position = genotypes_mut(document)
                .iter()
                .position(|genotype| genotype.case_id == altered)
                .expect("altered case");
            assert!(position > 0, "the altered case must not be the first");
        }));
    }

    /// DT-019, referencia. Una huella de referencia falsa en la cabecera no la
    /// ve ninguna otra comprobación de la publicación: el contrato estructural
    /// no conoce las bases de la referencia.
    #[test]
    fn wrong_reference_digest_in_the_header_without_vcf_is_rejected() {
        assert_judge_rejects(run_direct_cancel_without_vcf(|document| {
            let digest = &mut document.header.reference_content_sha256;
            let first = if digest.starts_with('0') { "1" } else { "0" };
            digest.replace_range(..1, first);
        }));
    }

    // ---- DT-020 y DT-021: el conjunto de casos de la verdad ----

    use crate::truth_schema::{TruthHistoryEvent, TruthModification};

    /// Posición, en base 1, de la fila de matriz que este escenario añade con
    /// frecuencia alternativa casi nula: cada individuo de la población tiene que
    /// publicar ahí un caso todo-REF, aunque ninguna copia lleve el ALT.
    const ALL_REFERENCE_LOCUS_1BASED: u64 = 310;

    /// Posición, en base 1, del cambio de base directo del escenario, que no
    /// es fila de la matriz.
    const DIRECT_SNP_1BASED: u64 = 40;

    /// El escenario de efecto neto nulo con repeticiones, ploidía 4 y dos
    /// individuos, con tres eventos más y una fila de matriz de frecuencia
    /// casi nula, con o sin VCF y con el gancho de alteración instalado.
    ///
    /// - `extra-ins`: inserción `GGA` en el punto 60 de la copia 0 de `ind-1`;
    /// - `extra-del`: deleción 250..253 de la copia 2 de `ind-1`, y
    ///   `extra-del-2`, 280..282 en la misma copia, que sigue en el historial
    ///   cuando la primera se traslada a los segmentos;
    /// - `seg-ins`: inserción `TTG` en el punto 300 de la copia 0 de `ind-1`,
    ///   aguas abajo de todo lo demás, para poder trasladarla a los segmentos;
    /// - `extra-snp` y `extra-snp-2`: dos cambios de base directos en la base
    ///   40 de las copias 0 y 1 de `ind-2`, con alternativos distintos, que no
    ///   es fila de la matriz. La copia 0 lleva el mayor en orden de base: el
    ///   catálogo del caso los ordena por base, no por copia;
    /// - la fila `chr1:310` con `alt_af = 0.01` (una fila con frecuencia 0 la
    ///   retira la preparación de la matriz: su ALT no existe en la población).
    fn run_case_set_scenario(
        with_vcf: bool,
        tamper: impl Fn(&mut TruthDocument) + 'static,
    ) -> Result<(), Diagnostic> {
        run_fixture_with_tamper(
            "direct-cancel",
            with_vcf,
            |config, snvs, reference| {
                config.push_str(
                    "\n[[direct_events]]\nevent_id = \"extra-ins\"\ncontig = \"chr1\"\nindividual_id = \"ind-1\"\nphase = 0\nkind = \"insertion\"\npoint = 60\nbases = \"GGA\"\n\n[[direct_events]]\nevent_id = \"extra-del\"\ncontig = \"chr1\"\nindividual_id = \"ind-1\"\nphase = 2\nkind = \"deletion\"\nstart = 250\nend = 253\n\n[[direct_events]]\nevent_id = \"extra-del-2\"\ncontig = \"chr1\"\nindividual_id = \"ind-1\"\nphase = 2\nkind = \"deletion\"\nstart = 280\nend = 282\n\n[[direct_events]]\nevent_id = \"seg-ins\"\ncontig = \"chr1\"\nindividual_id = \"ind-1\"\nphase = 0\nkind = \"insertion\"\npoint = 300\nbases = \"TTG\"\n",
                );
                let offset = usize::try_from(ALL_REFERENCE_LOCUS_1BASED - 1).expect("offset");
                let base = char::from(reference[offset]);
                let alternate = if base == 'A' { 'C' } else { 'A' };
                snvs.push_str(&format!(
                    "chr1\t{ALL_REFERENCE_LOCUS_1BASED}\t{base}\t{alternate}\tP1\t0.01\n"
                ));
                let offset = usize::try_from(DIRECT_SNP_1BASED - 1).expect("offset");
                let base = char::from(reference[offset]);
                let alternates = ['A', 'C', 'G', 'T']
                    .into_iter()
                    .filter(|candidate| *candidate != base)
                    .collect::<Vec<_>>();
                // La copia 0 lleva el mayor y la 1 el menor: en el orden de
                // las copias salen al revés que en el de base.
                for (event, phase, alternate) in [
                    ("extra-snp", 0, alternates[2]),
                    ("extra-snp-2", 1, alternates[0]),
                ] {
                    config.push_str(&format!(
                        "\n[[direct_events]]\nevent_id = \"{event}\"\ncontig = \"chr1\"\nindividual_id = \"ind-2\"\nphase = {phase}\nkind = \"snp\"\nposition_1based = {DIRECT_SNP_1BASED}\nref = \"{base}\"\nalt = \"{alternate}\"\n"
                    ));
                }
            },
            tamper,
        )
    }

    /// Copia un escenario de `tests/fixtures/vcf-null-effect/`, deja que
    /// `edit` cambie su configuración y su matriz (recibe además las bases de
    /// la referencia), quita o deja el VCF y lo ejecuta con `tamper`.
    fn run_fixture_with_tamper(
        name: &str,
        with_vcf: bool,
        edit: impl FnOnce(&mut String, &mut String, &[u8]),
        tamper: impl Fn(&mut TruthDocument) + 'static,
    ) -> Result<(), Diagnostic> {
        let root = std::env::temp_dir().join(format!(
            "gps-dt020-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("root");
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/vcf-null-effect")
            .join(name);
        fs::copy(fixture.join("reference.fa"), root.join("reference.fa")).expect("FASTA");
        let reference = fs::read_to_string(fixture.join("reference.fa"))
            .expect("FASTA")
            .lines()
            .filter(|line| !line.starts_with('>'))
            .collect::<String>()
            .into_bytes();
        let mut config = fs::read_to_string(fixture.join("general.toml")).expect("config");
        let mut snvs = fs::read_to_string(fixture.join("snvs.tsv")).expect("matrix");
        edit(&mut config, &mut snvs, &reference);
        if !with_vcf {
            config = config.replace("[vcf_output]\nmode = \"text\"\n", "");
        }
        assert_eq!(
            config.contains("[vcf_output]"),
            with_vcf,
            "the scenario must have VCF exactly when asked"
        );
        fs::write(root.join("snvs.tsv"), snvs).expect("write matrix");
        let path = root.join("general.toml");
        fs::write(&path, config).expect("write config");
        let effective = load_and_validate_config(&path).expect("valid config");
        let _guard = crate::truth_projection::test_tamper::install(tamper);
        let result = execute_validated_config_with_effective_document(
            &effective,
            &root.join("out"),
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            |_| Ok(()),
        )
        .map(|_| ());
        fs::remove_dir_all(root).expect("cleanup");
        result
    }

    /// El rechazo tiene que ser el del juez del **conjunto** de casos, no el
    /// de otra comprobación: lo exige `reason`, que es el nombre de su error
    /// (`MissingCase`, que solo da el conjunto) o el texto exacto de su motivo
    /// (los `GenotypeMismatch` del conjunto, cuyo nombre comparte con el
    /// contraste caso a caso).
    fn assert_case_set_rejects(result: Result<(), Diagnostic>, reason: &str) {
        let error = result.expect_err("the tampered producer must not publish");
        assert_eq!(error.code(), "GPS-INTERNAL-001", "{error:?}");
        let debug = format!("{error:?}");
        assert!(
            debug.contains("GenotypeJudgement") && debug.contains(reason),
            "the case-set judge must reject with {reason}: {error:?}"
        );
    }

    const UNANCHORED: &str = "no corresponde a ningún cambio";
    const DUPLICATED: &str = "publica ya ese locus o ese evento";

    /// Quita el primer caso que cumpla `select` y baja el recuento de la
    /// cabecera, para que el contrato estructural no lo vea.
    fn remove_case(
        document: &mut TruthDocument,
        select: impl Fn(&GenotypeRecord) -> bool,
    ) -> GenotypeRecord {
        let index = document
            .records
            .iter()
            .position(
                |record| matches!(record, TruthRecord::Genotype(genotype) if select(genotype)),
            )
            .expect("the scenario must publish the selected case");
        let TruthRecord::Genotype(removed) = document.records.remove(index) else {
            unreachable!("selected a genotype");
        };
        document.header.genotype_count -= 1;
        removed
    }

    fn alternate_event_id(genotype: &GenotypeRecord) -> Option<&str> {
        genotype.alleles.iter().find_map(|allele| match allele {
            TruthAllele::Insertion { event_id, .. } | TruthAllele::Deletion { event_id, .. } => {
                Some(event_id.as_str())
            }
            TruthAllele::Reference | TruthAllele::Snp { .. } | TruthAllele::Repeat { .. } => None,
        })
    }

    fn is_event_case(event_id: &'static str) -> impl Fn(&GenotypeRecord) -> bool {
        move |genotype| alternate_event_id(genotype) == Some(event_id)
    }

    /// El caso todo-REF de `ind-1` en la fila de frecuencia nula.
    fn is_all_reference_case(genotype: &GenotypeRecord) -> bool {
        genotype.case_kind == GenotypeCaseKind::Snp
            && genotype.individual_id == "ind-1"
            && genotype.locus.start == ALL_REFERENCE_LOCUS_1BASED - 1
    }

    fn all_reference_case(document: &TruthDocument) -> GenotypeRecord {
        let genotype = document
            .records
            .iter()
            .find_map(|record| match record {
                TruthRecord::Genotype(genotype) if is_all_reference_case(genotype) => {
                    Some(genotype.clone())
                }
                _ => None,
            })
            .expect("the zero-frequency matrix row must publish a case");
        assert!(
            genotype
                .allele_dosages
                .iter()
                .skip(1)
                .all(|dosage| *dosage == 0)
                && genotype.allele_dosages[0] > 0,
            "the case must be all-REF: {:?}",
            genotype.allele_dosages
        );
        genotype
    }

    /// DT-020. Sin el caso de una inserción que sigue en la copia, la verdad
    /// dejaría de explicar esas bases. Sin VCF se publicaba; con VCF lo
    /// paraba el agrupador, no el juez de la verdad.
    #[test]
    fn omitted_insertion_case_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    remove_case(document, is_event_case("extra-ins"));
                }),
                "MissingCase",
            );
        }
    }

    /// DT-020, con una deleción: las bases que faltan en la copia no las
    /// explica ningún caso.
    #[test]
    fn omitted_deletion_case_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    remove_case(document, is_event_case("extra-del"));
                }),
                "MissingCase",
            );
        }
    }

    /// DT-020, con una repetición: el tracto realizado de una copia sin el
    /// caso que lo publica.
    #[test]
    fn omitted_repeat_case_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    remove_case(document, |genotype| {
                        genotype.case_kind == GenotypeCaseKind::Repeat
                    });
                }),
                "MissingCase",
            );
        }
    }

    /// DT-021, segunda limitación. Un caso SNP cuyas cuatro copias llevan el
    /// alternativo se quitaba sin que nadie lo viera, **también con VCF**,
    /// porque el caso de otro individuo en el mismo locus le daba al VCF el
    /// alternativo que buscaba. La ayuda comprueba que el escenario lo
    /// reproduce de verdad: dosis `[0, 4]` y otro individuo con el ALT ahí.
    #[test]
    fn omitted_all_alternate_snp_case_covered_by_another_individual_is_rejected() {
        fn tamper(document: &mut TruthDocument) {
            let covered = |genotype: &GenotypeRecord, document: &TruthDocument| {
                document.records.iter().any(|record| {
                    matches!(record, TruthRecord::Genotype(other)
                        if other.case_kind == GenotypeCaseKind::Snp
                            && other.individual_id != genotype.individual_id
                            && other.locus == genotype.locus
                            && other.allele_dosages.iter().skip(1).any(|dosage| *dosage > 0))
                })
            };
            let snapshot = document.clone();
            let removed = remove_case(document, |genotype| {
                genotype.case_kind == GenotypeCaseKind::Snp
                    && genotype.allele_dosages[0] == 0
                    && genotype.allele_dosages.iter().sum::<u8>() == 4
                    && covered(genotype, &snapshot)
            });
            assert_eq!(removed.allele_dosages.iter().sum::<u8>(), 4);
        }
        for with_vcf in [false, true] {
            assert_case_set_rejects(run_case_set_scenario(with_vcf, tamper), "MissingCase");
        }
    }

    /// DT-021, primera limitación. Un caso todo-REF desplazado a un sitio en
    /// que todas las copias son referencia pasaba el contraste de portadores:
    /// nada lo ataba a su fila de matriz. Lo caza primero que la fila se quede
    /// sin caso; que el sitio nuevo no sea fila ni cambio de ninguna copia lo
    /// ejercita la prueba del caso fabricado, que no deja ninguna fila sin
    /// caso.
    #[test]
    fn shifted_all_reference_case_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    all_reference_case(document);
                    for genotype in genotypes_mut(document) {
                        if is_all_reference_case(genotype) {
                            genotype.locus.start = 3;
                            genotype.locus.end = 4;
                        }
                    }
                    sort_genotypes(document);
                }),
                "MissingCase",
            );
        }
    }

    /// DT-021, primera limitación: un caso todo-REF inventado en un sitio que
    /// no es fila de matriz ni cambio de ninguna copia.
    #[test]
    fn fabricated_all_reference_case_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    let mut fabricated = all_reference_case(document);
                    fabricated.case_id = "CASE_fabricated".to_owned();
                    fabricated.locus.start = 3;
                    fabricated.locus.end = 4;
                    document.records.push(TruthRecord::Genotype(fabricated));
                    document.header.genotype_count += 1;
                    sort_genotypes(document);
                }),
                UNANCHORED,
            );
        }
    }

    /// DT-021, primera limitación: el mismo caso todo-REF publicado dos veces
    /// con otro identificador. Cada copia se contrasta igual en los dos, así
    /// que solo el conjunto lo ve.
    #[test]
    fn duplicated_case_under_another_case_id_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    let mut duplicate = all_reference_case(document);
                    let last = duplicate.case_id.pop().expect("case id");
                    duplicate.case_id.push(if last == 'f' { 'e' } else { 'f' });
                    document.records.push(TruthRecord::Genotype(duplicate));
                    document.header.genotype_count += 1;
                    sort_genotypes(document);
                }),
                DUPLICATED,
            );
        }
    }

    /// DT-021, tercera limitación. Un portador que es copia de **otro**
    /// individuo es una verdad que el propio programa ha escrito mal: un
    /// defecto interno, no un error de formato (`GPS-FMT-006`) que señale a la
    /// entrada del usuario.
    #[test]
    fn carrier_of_another_individual_is_an_internal_defect_with_and_without_vcf() {
        for with_vcf in [false, true] {
            let error = run_case_set_scenario(with_vcf, |document| {
                let genotype = genotypes_mut(document)
                    .into_iter()
                    .find(|genotype| genotype.case_kind == GenotypeCaseKind::Snp)
                    .expect("a SNP case");
                let other = if genotype.individual_id == "ind-1" {
                    "ind-2"
                } else {
                    "ind-1"
                };
                genotype.carriers[0].copy.individual_id = other.to_owned();
            })
            .expect_err("a structurally corrupt truth must not publish");
            assert_eq!(error.code(), "GPS-INTERNAL-001", "{error:?}");
            assert!(format!("{error:?}").contains("Contract"), "{error:?}");
        }
    }

    /// El caso de la deleción que borra entera una inserción anterior no
    /// cambia la copia, pero el evento ocurrió y su caso tiene que estar
    /// (`schemas/truth-v4.md`, «Un caso sin registro VCF»). Ninguna
    /// diferencia de secuencia lo delata: solo el historial.
    #[test]
    fn omitted_null_effect_deletion_case_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_fixture_with_tamper(
                    "inserted-then-deleted",
                    with_vcf,
                    |_, _, _| {},
                    |document| {
                        let removed = remove_case(document, is_event_case("b"));
                        assert_eq!(removed.locus.start, removed.locus.end, "null deletion");
                    },
                ),
                "MissingCase",
            );
        }
    }

    /// Control: la inserción borrada entera **no** tiene caso, a propósito, y
    /// el escenario publica con VCF y sin él.
    #[test]
    fn inserted_then_deleted_publishes_with_and_without_vcf() {
        for with_vcf in [false, true] {
            run_fixture_with_tamper("inserted-then-deleted", with_vcf, |_, _, _| {}, |document| {
                assert!(
                    !document.records.iter().any(|record| matches!(record,
                        TruthRecord::Genotype(genotype) if alternate_event_id(genotype) == Some("a"))),
                    "the deleted insertion has no case"
                );
            })
            .expect("a legitimate truth publishes");
        }
    }

    /// Lleva un evento del historial de una copia a sus segmentos, como un
    /// cambio heredado, y quita su caso. La copia reconstruida es la misma
    /// —coinciden las huellas— y el historial ya no tiene el evento, de modo
    /// que ni el contraste de copias ni la exigencia de un caso por evento lo
    /// ven: **solo** que una diferencia de la copia con la referencia quede
    /// sin caso. Es el escenario saturado de esa comprobación.
    fn move_event_to_segments(document: &mut TruthDocument, event: &'static str) {
        let case = remove_case(document, is_event_case(event));
        let mut moved = 0;
        for record in &mut document.records {
            let TruthRecord::HaplotypeCopy(copy) = record else {
                continue;
            };
            let Some(index) = copy.history.iter().position(|historical| {
                matches!(historical,
                    TruthHistoryEvent::Insertion { event_id, .. }
                    | TruthHistoryEvent::Deletion { event_id, .. } if event_id == event)
            }) else {
                continue;
            };
            let modification = match copy.history.remove(index) {
                TruthHistoryEvent::Insertion {
                    event_id,
                    bases,
                    provenance,
                    ..
                } => TruthModification::Insertion {
                    event_id,
                    point: case.locus.start,
                    bases,
                    provenance,
                },
                TruthHistoryEvent::Deletion {
                    event_id,
                    provenance,
                    ..
                } => TruthModification::Deletion {
                    event_id,
                    start: case.locus.start,
                    end: case.locus.end,
                    provenance,
                },
                TruthHistoryEvent::Snp { .. } | TruthHistoryEvent::RepeatAllele { .. } => {
                    unreachable!("selected a structural event")
                }
            };
            assert_eq!(copy.segments.len(), 1, "one reference segment");
            let modifications = &mut copy.segments[0].modifications;
            let start = |modification: &TruthModification| match modification {
                TruthModification::Snp { position, .. } => *position,
                TruthModification::Insertion { point, .. } => *point,
                TruthModification::Deletion { start, .. } => *start,
            };
            modifications.push(modification);
            modifications.sort_by_key(start);
            moved += 1;
        }
        assert_eq!(moved, 1, "the event lives in exactly one copy");
    }

    /// Bases insertadas que ningún evento del historial ni ningún caso
    /// explican. Sin la comprobación de que todo fragmento insertado tenga su
    /// evento, el caso quitado no lo echaría en falta nadie.
    #[test]
    fn insertion_moved_to_segments_without_its_case_is_rejected() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    move_event_to_segments(document, "seg-ins");
                }),
                "MissingCase",
            );
        }
    }

    /// Bases de referencia que faltan en el estado heredado de la copia, sin
    /// deleción del historial ni caso que lo explique. Solo la comprobación del
    /// estado heredado lo ve.
    #[test]
    fn deletion_moved_to_segments_without_its_case_is_rejected() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    move_event_to_segments(document, "extra-del");
                }),
                "MissingCase",
            );
        }
    }

    /// Control de los escenarios anteriores: sin alteración, publican.
    #[test]
    fn case_set_scenario_publishes_untampered_with_and_without_vcf() {
        for with_vcf in [false, true] {
            run_case_set_scenario(with_vcf, |document| {
                all_reference_case(document);
            })
            .expect("the untampered scenario publishes");
        }
    }

    // ---- Segunda revisión de DT-020 y DT-021 ----

    fn snp_case_mut<'a>(
        document: &'a mut TruthDocument,
        individual: &str,
        position_1based: u64,
    ) -> &'a mut GenotypeRecord {
        genotypes_mut(document)
            .into_iter()
            .find(|genotype| {
                genotype.case_kind == GenotypeCaseKind::Snp
                    && genotype.individual_id == individual
                    && genotype.locus.start == position_1based - 1
            })
            .unwrap_or_else(|| panic!("no SNP case of {individual} at {position_1based}"))
    }

    const WRONG_ALLELES: &str = "los alelos del caso SNP no son";

    /// Un caso todo-REF atado a su posición pero con otro alternativo que el
    /// de su fila: nadie lo lleva, así que el contraste de portadores no lo
    /// ve, y la posición sigue siendo la de la fila.
    #[test]
    fn all_reference_case_with_another_alternate_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    all_reference_case(document);
                    let genotype = snp_case_mut(document, "ind-1", ALL_REFERENCE_LOCUS_1BASED);
                    let TruthAllele::Snp { base } = &mut genotype.alleles[1] else {
                        panic!("the row alternate");
                    };
                    *base = if base == "G" { "T" } else { "G" }.to_owned();
                }),
                WRONG_ALLELES,
            );
        }
    }

    /// Un alternativo de más, que no lleva ninguna copia, en el caso de un
    /// cambio de base directo y en el de una fila de la matriz.
    #[test]
    fn snp_case_with_an_extra_alternate_is_rejected_with_and_without_vcf() {
        for (individual, position) in [
            ("ind-2", DIRECT_SNP_1BASED),
            ("ind-1", ALL_REFERENCE_LOCUS_1BASED),
        ] {
            for with_vcf in [false, true] {
                assert_case_set_rejects(
                    run_case_set_scenario(with_vcf, move |document| {
                        let genotype = snp_case_mut(document, individual, position);
                        let present = genotype
                            .alleles
                            .iter()
                            .filter_map(|allele| match allele {
                                TruthAllele::Snp { base } => Some(base.clone()),
                                _ => None,
                            })
                            .collect::<Vec<_>>();
                        let extra = ["A", "C", "G", "T"]
                            .into_iter()
                            .find(|base| !present.iter().any(|known| known == base))
                            .expect("a free base");
                        genotype.alleles.push(TruthAllele::Snp {
                            base: extra.to_owned(),
                        });
                        genotype.allele_dosages.push(0);
                    }),
                    WRONG_ALLELES,
                );
            }
        }
    }

    /// Una deleción heredada de un segmento escondida dentro de otra del
    /// historial: la del historial se acorta una base y la del segmento quita
    /// la que falta, de modo que la copia y el locus son los mismos. El locus
    /// de la del historial abarca la base escondida; solo mirar el estado
    /// heredado lo ve.
    #[test]
    fn inherited_deletion_hidden_inside_a_history_deletion_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    let locus = genotypes_mut(document)
                        .into_iter()
                        .find(|genotype| alternate_event_id(genotype) == Some("extra-del"))
                        .expect("extra-del case")
                        .locus
                        .clone();
                    for genotype in genotypes_mut(document) {
                        for allele in &mut genotype.alleles {
                            if let TruthAllele::Deletion {
                                event_id,
                                end,
                                source_event:
                                    TruthStructuralEvent::Deletion {
                                        end: source_end, ..
                                    },
                                ..
                            } = allele
                                && event_id == "extra-del"
                            {
                                *end -= 1;
                                *source_end -= 1;
                            }
                        }
                    }
                    let mut hidden = 0;
                    for record in &mut document.records {
                        let TruthRecord::HaplotypeCopy(copy) = record else {
                            continue;
                        };
                        let mut provenance = None;
                        for event in &mut copy.history {
                            if let TruthHistoryEvent::Deletion {
                                event_id,
                                end,
                                provenance: source,
                                ..
                            } = event
                                && event_id == "extra-del"
                            {
                                *end -= 1;
                                provenance = Some(source.clone());
                            }
                        }
                        let Some(provenance) = provenance else {
                            continue;
                        };
                        let modifications = &mut copy.segments[0].modifications;
                        modifications.push(TruthModification::Deletion {
                            event_id: "hidden".to_owned(),
                            start: locus.start + 1,
                            end: locus.start + 2,
                            provenance,
                        });
                        modifications.sort_by_key(|modification| match modification {
                            TruthModification::Snp { position, .. } => *position,
                            TruthModification::Insertion { point, .. } => *point,
                            TruthModification::Deletion { start, .. } => *start,
                        });
                        hidden += 1;
                    }
                    assert_eq!(hidden, 1, "one copy carries extra-del");
                }),
                "estado heredado",
            );
        }
    }

    /// Las copias de `ind-1` declaran otra población que la suya: todas, o
    /// solo la última, para que mirar solo la primera no baste.
    #[test]
    fn copies_declaring_another_population_are_rejected_with_and_without_vcf() {
        for only_last in [false, true] {
            for with_vcf in [false, true] {
                assert_case_set_rejects(
                    run_case_set_scenario(with_vcf, move |document| {
                        for record in &mut document.records {
                            if let TruthRecord::HaplotypeCopy(copy) = record
                                && copy.identity.individual_id == "ind-1"
                                && (!only_last || copy.identity.copy_index == 3)
                            {
                                copy.provenance.origin = "initial_population:P2".to_owned();
                            }
                        }
                    }),
                    "CopyOriginMismatch",
                );
            }
        }
    }

    /// Los casos todo-REF de dos individuos en la misma fila intercambiados:
    /// alelos, portadores y dosis cuadran en los dos, y solo los delata que el
    /// identificador de cada caso es el del otro individuo.
    #[test]
    fn all_reference_cases_swapped_between_individuals_are_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    let first = snp_case_mut(document, "ind-1", ALL_REFERENCE_LOCUS_1BASED).clone();
                    let second =
                        snp_case_mut(document, "ind-2", ALL_REFERENCE_LOCUS_1BASED).clone();
                    assert_ne!(first.case_id, second.case_id);
                    assert!(
                        [&first, &second].iter().all(|genotype| genotype
                            .allele_dosages
                            .iter()
                            .skip(1)
                            .all(|dosage| *dosage == 0)),
                        "both cases are all-REF"
                    );
                    snp_case_mut(document, "ind-1", ALL_REFERENCE_LOCUS_1BASED)
                        .case_id
                        .clone_from(&second.case_id);
                    snp_case_mut(document, "ind-2", ALL_REFERENCE_LOCUS_1BASED)
                        .case_id
                        .clone_from(&first.case_id);
                }),
                "el identificador del caso no es",
            );
        }
    }

    /// Un caso de repetición inventado en un intervalo en que ninguna copia
    /// lleva un tracto: todas sus copias son referencia ahí, así que el
    /// contraste caso a caso lo acepta, y solo el conjunto lo rechaza.
    #[test]
    fn fabricated_repeat_case_is_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    let mut fabricated = genotypes_mut(document)
                        .into_iter()
                        .find(|genotype| genotype.case_kind == GenotypeCaseKind::Repeat)
                        .expect("a repeat case")
                        .clone();
                    fabricated.case_id = "CASE_fabricated_repeat".to_owned();
                    fabricated.locus.start += 1;
                    fabricated.locus.end += 1;
                    for carrier in &mut fabricated.carriers {
                        carrier.allele_index = 0;
                    }
                    let copies = u8::try_from(fabricated.carriers.len()).expect("ploidy");
                    fabricated.allele_dosages.fill(0);
                    fabricated.allele_dosages[0] = copies;
                    document.records.push(TruthRecord::Genotype(fabricated));
                    document.header.genotype_count += 1;
                    sort_genotypes(document);
                }),
                UNANCHORED,
            );
        }
    }

    /// Los alelos de un caso de repetición permutados, con sus portadores y su
    /// identificador rehechos: todo cuadra salvo el orden, que deja de ser el
    /// del productor y hace que el identificador ya no sea canónico.
    #[test]
    fn permuted_repeat_alleles_are_rejected_with_and_without_vcf() {
        for with_vcf in [false, true] {
            assert_case_set_rejects(
                run_case_set_scenario(with_vcf, |document| {
                    let genotype = genotypes_mut(document)
                        .into_iter()
                        .find(|genotype| {
                            genotype.case_kind == GenotypeCaseKind::Repeat
                                && genotype.alleles.len() >= 3
                        })
                        .expect("a repeat case with two alternates");
                    genotype.alleles.swap(1, 2);
                    genotype.allele_dosages.swap(1, 2);
                    for carrier in &mut genotype.carriers {
                        carrier.allele_index = match carrier.allele_index {
                            1 => 2,
                            2 => 1,
                            other => other,
                        };
                    }
                    genotype.case_id = crate::truth_projection::case_id(
                        "P1",
                        &genotype.individual_id,
                        &genotype.locus.sequence_id,
                        genotype.locus.start + 1,
                        &genotype.alleles,
                    );
                    sort_genotypes(document);
                }),
                "los alelos del caso de repetición",
            );
        }
    }
}
