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
use crate::identity::{ReferenceId, SequenceId, derive_execution_id};
use crate::initial_individual_integration::IncompatibleGenotypeAuthorization;
use crate::initial_individual_materialization::{
    InitialIndividualsMaterialization, InitialIndividualsMaterializationError,
    materialize_initial_individuals,
};
use crate::integrity::{Sha256Writer, configuration_digest};
use crate::materialization::SnpLocus;
use crate::population_snvs::{
    PopulationSnvsProvenance, canonical_metadata_document, canonical_tsv_bytes,
    generate_population_snvs_from_reference_excluding, normalization_audit_document,
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
    // `validate` tiene que ver **lo mismo** que `run`, y eso incluye el
    // catálogo de microsatélites y el abanico que el deslizamiento biológico
    // genera: sin hacerlo aquí, un escenario con loci generados pasaba la
    // validación con los loci vacíos y solo fallaba al ejecutar, que es tarde.
    let (resolved, population_snvs) =
        resolve_repeat_catalog(effective, authorizations, &prepared_reference)?;
    let effective: &EffectiveConfig = resolved.as_ref().unwrap_or(effective);
    preflight_prepared_config(
        effective,
        authorizations,
        &prepared_reference,
        population_snvs,
    )
    .map(|_| ())
}

/// **Las tres etapas de una simulación, en su orden** (REQ-REF-38): con la
/// referencia final ya preparada, detectar sus microsatélites, preparar la
/// matriz de SNV —la generada excluye los loci del catálogo—, retirar los
/// detectados que una variante suministrada o escrita a mano invade
/// (REQ-CONF-143) y generar los abanicos de alelos por deslizamiento. Al
/// final, sobre el catálogo que queda, la regla de la PCR y la desconexión del
/// tartamudeo clonal si no queda ningún locus (`with_final_repeat_catalog`).
///
/// Devuelve la configuración resultante, o `None` si ninguna etapa la cambió,
/// y la matriz ya preparada. Sin la sección de microsatélites, detectar y
/// retirar no hacen nada, la exclusión está vacía y todo es como antes.
///
/// Corre **antes** de crear nada publicado, en `validate` y en `run`.
fn resolve_repeat_catalog(
    effective: &EffectiveConfig,
    authorizations: ExecutionAuthorizations,
    prepared_reference: &PreparedReference,
) -> Result<
    (
        Option<EffectiveConfig>,
        Option<PreparedPopulationSnvsForExecution>,
    ),
    Diagnostic,
> {
    preflight_vcf_samples(effective)?;
    preflight_reductions(effective, prepared_reference.reference())?;
    let detected = effective.with_detected_repeat_loci(prepared_reference.reference())?;
    let with_catalog: &EffectiveConfig = detected.as_ref().unwrap_or(effective);
    // La guarda de admisión, sobre el catálogo entero: un detectado la pasa
    // por construcción; un declarado, con sus mínimos propios, sobre su
    // tramo; y un declarado de 7 pb o más, que la detección no ve, solo por
    // ella.
    admit_declared_repeat_loci(with_catalog, prepared_reference.reference())?;
    let population_snvs =
        prepare_population_snvs_for_execution(with_catalog, authorizations, prepared_reference)?;
    let rows: Vec<(&str, u64, &str)> = population_snvs
        .as_ref()
        .map(|prepared| {
            prepared
                .matrix
                .rows()
                .iter()
                .map(|row| {
                    (
                        row.contig().as_str(),
                        row.position_1based(),
                        row.variant_id().as_str(),
                    )
                })
                .collect()
        })
        .unwrap_or_default();
    let retired = with_catalog.with_retired_repeat_loci(&rows)?;
    drop(rows);
    let without_retired: &EffectiveConfig = retired.as_ref().unwrap_or(with_catalog);
    let generated =
        without_retired.with_generated_repeat_alleles(prepared_reference.reference())?;
    let with_alleles: &EffectiveConfig = generated.as_ref().unwrap_or(without_retired);
    // Sobre el catálogo final: la PCR y el tartamudeo clonal (ver
    // `with_final_repeat_catalog`).
    let finished = with_alleles.with_final_repeat_catalog()?;
    Ok((
        finished.or(generated).or(retired).or(detected),
        population_snvs,
    ))
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
    let (resolved, population_snvs) =
        resolve_repeat_catalog(effective, authorizations, &prepared_reference)?;
    let effective: &EffectiveConfig = resolved.as_ref().unwrap_or(effective);
    // Todo error de dominio o preflight ocurre antes de crear el provisional y
    // antes de que exista cualquier RNG local de selección.
    let prepared_execution = preflight_prepared_config(
        effective,
        authorizations,
        &prepared_reference,
        population_snvs,
    )?;
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
    // La unidad genotípica (DT-038) con los registros de cada artefacto **tal
    // como se validaron al escribirlos**: la verdad por un lado y, por otro, la
    // matriz de SNV y el VCF con su índice. Las referencias de las que depende
    // ya tienen sus propios registros (`reference_artifact` y
    // `inserted_publication`).
    let mut truth_record: Option<ArtifactRecord> = None;
    let mut genotype_records: Vec<ArtifactRecord> = Vec::new();
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
        genotype_records.extend(
            workspace
                .artifacts()
                .iter()
                .filter(|artifact| {
                    matches!(
                        artifact.id.as_str(),
                        "population-snvs-tsv" | "population-snvs-metadata"
                    )
                })
                .cloned(),
        );
        // **El catálogo de microsatélites aparte** (REQ-MAN-98), solo con la
        // sección: sin ella no se escribe nada y los bytes son los de antes.
        let repeat_catalog_link =
            publish_repeat_catalog(effective, simulation_artifact, &mut workspace)?;
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
        let (truth, written_truth) = write_execution_truth(
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
        truth_record = Some(written_truth);
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
            genotype_records.push(published.vcf.clone());
            genotype_records.extend(published.index.clone());
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
        if let Some(link) = repeat_catalog_link {
            effective_document["repeat_detection"]["catalog"] = link;
        }
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
            // REQ-SAL-25, DTC-045 y REQ-SAL-27, **una sola regla**: una reducción
            // completa cuyos bytes en disco siguen siendo los que se validaron al
            // escribirlos se conserva entera, haya evento o no y falle lo que
            // falle; la que no pasa esa comprobación, o la que causó el fallo,
            // queda fallida y sin artefactos. Solo se comprueba que los bytes en
            // disco son los que se validaron al escribirlos; no se repite el juicio
            // final (reconstruir la derivada desde subsemilla y selección y
            // compararla base a base), de modo que un defecto del propio productor
            // en la selección no se detectaría aquí.
            let kept = published_reductions
                .iter()
                .filter(|published| {
                    manifest_reductions[published.effective_index].state == UnitState::Complete
                        && workspace.disk_still_matches(&published.fasta)
                        && workspace.disk_still_matches(&published.map)
                })
                .map(|published| published.effective_index)
                .collect::<BTreeSet<_>>();
            for (index, reduction) in manifest_reductions.iter_mut().enumerate() {
                if kept.contains(&index) {
                    continue;
                }
                reduction.state = UnitState::Failed;
                reduction.fasta_artifact_id = None;
                reduction.map_artifact_id = None;
                reduction.source_fasta_artifact_id = None;
            }
            // DT-038, la **unidad genotípica**: la verdad, la matriz de SNV que
            // cita, el VCF con su índice si se pidió y la referencia de la
            // simulación (la suministrada o, con inserción, la final con su
            // mapa). Se conserva **entera o nada**, con la misma regla que una
            // reducción: cada artefacto en disco tiene que ser byte a byte el
            // que se validó al escribirlo. Si el VCF se pidió y no llegó a
            // publicarse, la unidad entera queda fuera: una verdad sin su VCF
            // no es lo que se pidió, y además podría superar el tope de bases
            // de tracto que el VCF hizo saltar. Como en las reducciones, no se
            // repite el juicio final: una verdad que el juez final rechaza por
            // un defecto del productor, con los bytes intactos, se conserva.
            //
            // La verdad cita por su huella la copia de la configuración (las
            // fuentes `explicit` de los eventos directos y `random` de los
            // loci de repetición), así que la unidad exige también que esa
            // copia se conserve, con **la misma** comprobación con que el
            // paquete fallido decide conservarla. Es lo único que la unidad
            // cita fuera de sí misma: la matriz y sus metadatos citan la
            // referencia de la simulación, y el VCF, la verdad y la referencia.
            let genotype_unit_kept = truth_record.is_some()
                && (effective.vcf_output().is_none() || !vcf_outputs.is_empty())
                && configuration_input.as_ref().is_some_and(|input| {
                    valid_failure_configuration_copy(&workspace, effective, &input.copy)
                })
                && truth_record
                    .iter()
                    .chain(&genotype_records)
                    .chain(reference_artifact.iter())
                    .chain(
                        inserted_publication
                            .iter()
                            .flat_map(|inserted| [&inserted.fasta, &inserted.map]),
                    )
                    .all(|record| workspace.disk_still_matches(record));
            let failure_references = manifest_references
                .iter()
                .filter(|reference| {
                    match manifest_reductions
                        .iter()
                        .find(|reduction| reduction.id == reference.id)
                    {
                        Some(reduction) => reduction.state == UnitState::Complete,
                        // Sin evento ni unidad genotípica solo sobrevive la
                        // suministrada, y solo si alguna reducción se conserva
                        // y la nombra de origen.
                        None => {
                            event_links.is_some()
                                || genotype_unit_kept
                                || (!kept.is_empty() && reference.id == effective.reference_id())
                        }
                    }
                })
                .cloned()
                .collect::<Vec<_>>();
            // **Con lecturas o verdad ya publicadas, la referencia final se
            // conserva con su mapa** (DTC-098): salieron de ella, y el paquete
            // fallido tiene que poder relacionarla con la suministrada
            // (REQ-MAN-07). Sin ninguna de las dos no se conserva nada de ella.
            let failure_insertion = (event_links.is_some() || genotype_unit_kept)
                .then_some(())
                .and(inserted_publication.as_ref())
                .zip(reference_artifact.as_ref());
            let genotype_unit = genotype_unit_kept.then(|| FailedGenotypeUnit {
                artifact_ids: truth_record
                    .iter()
                    .chain(&genotype_records)
                    .map(|record| record.id.clone())
                    .collect(),
                population_snvs: population_snvs.as_ref(),
                vcf_outputs: &vcf_outputs,
                repeat_draws: prepared_execution
                    .initial_individuals
                    .as_ref()
                    .map_or(&[][..], InitialIndividualsMaterialization::repeat_draws),
                realized_tract_output: &realized_tract_output,
            });
            let first = publish_controlled_failure(
                &mut workspace,
                execution_id.as_str(),
                effective,
                authorizations,
                &manifest_reductions,
                event_links.as_deref(),
                &failure_references,
                failure_insertion,
                genotype_unit.as_ref(),
                std::slice::from_ref(&cause),
            );
            // DT-039, **el último recurso**. Si el primer paquete fallido no
            // pasa sus jueces —porque alguien alteró en disco un artefacto que
            // conservaba, como la referencia suministrada—, se publica uno
            // mínimo: la copia de la configuración si sigue intacta, las
            // reducciones marcadas `failed`, ningún evento ni unidad
            // genotípica, y en los diagnósticos la causa y detrás el fallo
            // secundario. Se devuelve **siempre la causa original**; si ni
            // siquiera el mínimo o la publicación salen, con el código del
            // fallo secundario como contexto.
            if let Err(secondary) = first {
                for reduction in &mut manifest_reductions {
                    reduction.state = UnitState::Failed;
                    reduction.fasta_artifact_id = None;
                    reduction.map_artifact_id = None;
                    reduction.source_fasta_artifact_id = None;
                }
                if publish_controlled_failure(
                    &mut workspace,
                    execution_id.as_str(),
                    effective,
                    authorizations,
                    &manifest_reductions,
                    None,
                    &[],
                    None,
                    None,
                    &[cause.clone(), secondary.clone()],
                )
                .is_err()
                {
                    return Err(cause.with("secondary_failure_code", secondary.code()));
                }
            }
            if let Err(secondary) = workspace.publish() {
                return Err(cause.with("secondary_failure_code", secondary.code()));
            }
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
) -> Result<(TruthDocument, ArtifactRecord), Diagnostic> {
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
    let record = write_truth_artifact(workspace, &document)?;
    Ok((document, record))
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

/// Un fallo del evento técnico y, si pudo escribirse, el enlace con estado
/// `failed` a su manifiesto de evento fallido (DT-036).
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

/// Publica el único evento técnico real de F3 y devuelve su enlace.
///
/// Si el ensamblado falla, el evento **queda enlazado igualmente**, con estado
/// `failed` y un manifiesto de evento fallido que registra la causa
/// (REQ-SAL-21, DT-036). Solo si ni siquiera ese manifiesto puede escribirse
/// —o si falla la derivación de identidades, sin la cual no hay manifiesto
/// posible— el fallo sale sin enlace y el paquete fallido no nombra el evento.
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
    assemble_technical_event(
        workspace,
        effective,
        event,
        &identities,
        reference,
        reference_fasta,
        initial_individuals,
    )
    .map_err(|cause| {
        let link = write_failed_event_manifest(
            workspace,
            effective,
            event,
            &identities,
            reference_fasta,
            &cause,
        )
        .ok();
        Box::new(TechnicalEventAssemblyFailure { cause, link })
    })
}

/// Ensambla el evento exclusivamente desde FASTQ ya cerrados y validados por
/// el `ExecutionWorkspace`. El manifiesto de evento se vuelve un artefacto más
/// del catálogo antes de que la raíz lo enlace.
fn assemble_technical_event(
    workspace: &mut ExecutionWorkspace,
    effective: &EffectiveConfig,
    event: &crate::application::EffectiveTechnicalEvent,
    identities: &crate::technical_identity::TechnicalIdentityPlan,
    reference: &Reference,
    reference_fasta: &ArtifactRecord,
    initial_individuals: Option<&InitialIndividualsMaterialization>,
) -> Result<ManifestEventLink, Diagnostic> {
    let individuals = initial_individuals
        .ok_or_else(|| {
            technical_event_error("initial_individuals", "faltan individuos materializados")
        })?
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
    let profile = resolve_f4_technical_profile(event.profile())
        .map_err(|error| technical_event_error("technical_profile", error.to_string()))?;
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
        .collect::<Result<Vec<_>, Diagnostic>>()?;

    // El padrón fija el puesto de cada individua antes de recorrer nada, de
    // modo que el orden del inventario no dependa de en qué orden se generen
    // las parejas. Es lo que permite que una clave guarde un puesto en lugar de
    // repetir el identificador de su individua un millón de veces.
    let roster = crate::identity::IndividualRoster::new(individuals.keys().cloned())?;
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
                technical_event_error(
                    "technical_event.individual_ids",
                    "un individuo seleccionado no tiene verdad materializada",
                )
            })?;
        let seat = roster.seat(individual.id()).ok_or_else(|| {
            technical_event_error(
                "technical_event.individual_ids",
                "una individua seleccionada no figura en el padrón",
            )
        })?;
        let mut sampler = paired_sampler(
            reference,
            individual,
            identities,
            profile,
            event,
            &stutter_plan,
            &pcr_plan,
        )?;
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
            let pair = pair.map_err(|error| wgs_paired_sampler_error(error, individual.id()))?;
            keys.push(CleanPairKey::new(
                seat,
                pair.pre_hop_identity().pre_hop_instance_id(),
            )?);
        }
        // El informe del tartamudeo se recoge **de la primera pasada**, igual
        // que el de la quimera: las regeneraciones posteriores repiten las
        // mismas decisiones y contarlas otra vez las duplicaría.
        stutter_accumulator
            .merge(&sampler.stutter_report())
            .map_err(|error| technical_event_error("stutter", format!("{error:?}")))?;
        // Los recuentos de la amplificación se recogen **de la primera
        // pasada**, igual que los de quimera y tartamudeo: las regeneraciones
        // posteriores repiten las mismas familias y contarlas otra vez las
        // duplicaría.
        pcr_accumulator.merge(&sampler.pcr_report());
        realized_chimera_pairs = realized_chimera_pairs
            .checked_add(sampler.chimera_report().realized_pairs())
            .ok_or_else(|| technical_event_error("chimera", "contador de quimeras desbordado"))?;
        clean_libraries.push(
            CleanLibraryPairs::new(
                individual.id().clone(),
                seat,
                (*population_id).clone(),
                library.library_id().clone(),
                keys,
            )
            .map_err(|error| {
                technical_event_error("biological_contamination", format!("{error:?}"))
            })?,
        );
    }
    let contamination_plan = plan_biological_contamination(
        effective.reproducibility().root_seed(),
        event.event_id(),
        event.biological_contamination(),
        &clean_libraries,
    )
    .map_err(|error| technical_event_error("biological_contamination", format!("{error:?}")))?;
    let mut materializer = BiologicalContaminationMaterializer::new(
        effective.reproducibility().root_seed(),
        event.event_id().clone(),
        &clean_libraries,
    )
    .map_err(|error| technical_event_error("biological_contamination", format!("{error:?}")))?;
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
    .map_err(|error| technical_event_error("demultiplexing", format!("{error:?}")))?;
    let mut streams: BTreeMap<crate::identity::LibraryId, FinalPairedFastqStream> = BTreeMap::new();
    for library in &libraries {
        let individual = individuals
            .get(library.individual_id())
            .expect("selected above")
            .0;
        let udi = udis.get(library.library_id()).ok_or_else(|| {
            technical_event_error("technical_event.udi", "biblioteca sin UDI declarado")
        })?;
        let stream = FinalPairedFastqStream::open(
            workspace,
            event.event_id().as_str(),
            individual.id().as_str(),
            library.library_id(),
            udi,
            event.read_length(),
        )?;
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
            technical_event_error(
                "biological_contamination",
                "una donante no figura en el padrón",
            )
        })?;
        let collected = RegeneratedPairs::collect(
            individual_id,
            wanted,
            reference,
            &individuals,
            identities,
            profile,
            event,
            effective.reproducibility().root_seed(),
            &stutter_plan,
            &pcr_plan,
        )?;
        donor_pairs.insert(*seat, collected);
    }
    for recipient_plan in contamination_plan.recipients() {
        let clean = clean_libraries
            .iter()
            .find(|library| library.individual_id() == recipient_plan.recipient_id())
            .ok_or_else(|| {
                technical_event_error("biological_contamination", "receptora ausente")
            })?;
        let mut regenerated = RegeneratedPairs::build(
            clean,
            reference,
            &individuals,
            identities,
            profile,
            event,
            effective.reproducibility().root_seed(),
            &stutter_plan,
            &pcr_plan,
        )?;
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
                technical_event_error("biological_contamination", format!("{error:?}"))
            })?;
        biological_manifest.push(
            recipient_plan
                .realized_manifest_recipient(&report)
                .ok_or_else(|| {
                    technical_event_error(
                        "biological_contamination",
                        "informe realizado incoherente",
                    )
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
        let udi = udis
            .get(clean.library_id())
            .cloned()
            .ok_or_else(|| technical_event_error("demultiplexing", "UDI de biblioteca ausente"))?;
        let mut regenerated = RegeneratedPairs::build(
            clean,
            reference,
            &individuals,
            identities,
            profile,
            event,
            effective.reproducibility().root_seed(),
            &stutter_plan,
            &pcr_plan,
        )?;
        for key in clean.own_pairs() {
            let transformed = regenerated.take(key).ok_or_else(|| {
                technical_event_error(
                    "biological_contamination",
                    "la clave limpia no se puede regenerar",
                )
            })?;
            error_quality
                .record(transformed.report())
                .map_err(|error| technical_event_error("error_quality", format!("{error:?}")))?;
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
            )?;
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
        let stream = streams
            .remove(library.library_id())
            .ok_or_else(|| technical_event_error("demultiplexing", "biblioteca final ausente"))?;
        let publication = stream.finish(workspace)?;
        let basic_bases = clean_libraries
            .iter()
            .find(|clean| clean.library_id() == library.library_id())
            .and_then(|clean| u64::try_from(clean.own_pairs().len()).ok())
            .and_then(|pairs| pairs.checked_mul(2))
            .and_then(|reads| reads.checked_mul(u64::from(event.read_length())))
            .ok_or_else(|| technical_event_error("coverage", "bases WGS básicas desbordan"))?;
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
        .map_err(|error| technical_event_error("coverage", error.to_string()))?;
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
    let basic_wgs_coverage_x = coverage
        .iter()
        .map(|entry| entry.report.basic_wgs_x())
        .sum();
    write_event_manifest(
        workspace,
        event,
        identities,
        reference_fasta,
        UnitState::Complete,
        &[],
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
        event_manifest_stutter(&stutter_plan, stutter_accumulator.report()),
        event_manifest_pcr(
            event,
            pcr_accumulator.report(),
            pcr_bound,
            basic_wgs_coverage_x,
        ),
    )
}

/// Sección del tartamudeo del manifiesto de evento. Solo existe cuando el
/// evento abre de verdad la etapa `stutter`: con el tartamudeo apagado, o sin
/// loci declarados, el manifiesto sale idéntico al de antes de esta ampliación.
fn event_manifest_stutter(
    stutter_plan: &crate::stutter::StutterPlan,
    report: crate::stutter::StutterReport,
) -> Option<crate::event_manifest::EventManifestStutter> {
    stutter_plan
        .is_active()
        .then(|| crate::event_manifest::EventManifestStutter {
            settings: stutter_plan.settings().clone(),
            report,
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
        })
}

/// Sección de la amplificación del manifiesto de evento. Solo existe con
/// preparación amplificada: sin PCR, `artifacts.pcr` conserva exactamente los
/// bytes de antes de esta ampliación.
fn event_manifest_pcr(
    event: &crate::application::EffectiveTechnicalEvent,
    report: crate::pcr_amplification::PcrReport,
    bound: crate::pcr_amplification::PcrResourceBound,
    basic_wgs_coverage_x: f64,
) -> Option<crate::event_manifest::EventManifestPcr> {
    event
        .pcr()
        .map(|settings| crate::event_manifest::EventManifestPcr {
            settings: settings.clone(),
            report,
            bound,
            basic_wgs_coverage_x,
        })
}

/// Manifiesto de un evento cuyo ensamblado falló (REQ-SAL-21, DT-036).
///
/// Lleva estado `failed`, la causa en `diagnostics` y las 2n salidas FASTQ
/// marcadas `failed`, sin artefacto y con el código de la causa. Cobertura y
/// resúmenes van vacíos y todas las cifras realizadas a cero; el modo de
/// errores, la tasa de quimera y la de salto de índice son las declaradas.
///
/// **No conserva ningún FASTQ, ni siquiera los ya cerrados**: cada uno depende
/// del evento entero —el demultiplexado y la contaminación reparten parejas
/// entre bibliotecas—, así que un FASTQ suelto de un evento fallido no
/// describiría nada. Se retiran del catálogo antes de escribir el manifiesto,
/// y es lo que ya se hacía antes de enlazar el evento fallido (DT-002).
fn write_failed_event_manifest(
    workspace: &mut ExecutionWorkspace,
    effective: &EffectiveConfig,
    event: &crate::application::EffectiveTechnicalEvent,
    identities: &crate::technical_identity::TechnicalIdentityPlan,
    reference_fasta: &ArtifactRecord,
    cause: &Diagnostic,
) -> Result<ManifestEventLink, Diagnostic> {
    let without_fastq = workspace
        .artifacts()
        .iter()
        .filter(|artifact| artifact.kind != "fastq")
        .map(|artifact| artifact.id.clone())
        .collect::<BTreeSet<_>>();
    workspace.retain_only_artifacts(&without_fastq)?;
    let outputs = identities
        .libraries()
        .flat_map(|library| {
            crate::technical_identity::TechnicalMate::all().map(|mate| TechnicalEventArtifact {
                individual_id: library.individual_id().clone(),
                library_id: library.library_id().clone(),
                mate,
                id: format!("fastq-{}-r{}", library.library_id().as_str(), mate.number()),
                state: UnitState::Failed,
                artifact: None,
                diagnostic: Some(cause.clone()),
            })
        })
        .collect::<Vec<_>>();
    write_event_manifest(
        workspace,
        event,
        identities,
        reference_fasta,
        UnitState::Failed,
        std::slice::from_ref(cause),
        &outputs,
        &[],
        &[],
        EventManifestErrorQuality::from_report(
            ErrorQualitySettings::default().with_sequencing_errors(event.sequencing_errors()),
            ErrorQualityAccumulator::default().report(),
        ),
        EventManifestDemultiplexing::without_pairs(event.index_hopping_rate()),
        0,
        &[],
        event_manifest_stutter(
            &effective.stutter_plan(),
            crate::stutter::StutterAccumulator::default().report(),
        ),
        event_manifest_pcr(
            event,
            crate::pcr_amplification::PcrAccumulator::default().report(),
            crate::pcr_amplification::PcrResourceBound::empty(),
            0.0,
        ),
    )
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
fn wgs_paired_sampler_error(
    error: WgsSamplingError,
    individual: &crate::identity::IndividualId,
) -> Diagnostic {
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
        // DT-031: con una tasa de quimeras alta, una referencia corta puede no
        // ofrecer a alguna pareja seleccionada un segundo origen que no toque
        // su prefijo (DTC-084). Es un fallo previsto que se declara antes de
        // publicar, sin reintento ni reducción silenciosa de la tasa, y no un
        // defecto interno. Reutiliza `GPS-DOMAIN-003`, como DT-006, y no una
        // salida de configuración: depende de la semilla (la misma
        // configuración publica con la mayoría de las semillas y `validate` la
        // acepta), y DTC-044 reserva la salida 2 a la configuración inválida.
        WgsSamplingError::Chimera(crate::chimera::ChimeraError::NoSecondSource) => {
            Diagnostic::domain(
                "GPS-DOMAIN-003",
                "Una invariante del estado biológico o técnico válido no se cumple.",
            )
            .with("entity", individual.as_str())
            .with("invariant", "segundo_origen_de_quimera")
            .with(
                "detail",
                format!(
                    "con la tasa de quimeras pedida, a alguna pareja seleccionada del individuo {individual} \
no le queda en su referencia un segundo origen que no toque su prefijo. Depende de la semilla y de la \
longitud de las copias del individuo: otra semilla puede publicar. Usa una referencia más larga o una \
tasa de quimeras menor; el simulador no la reduce por su cuenta."
                ),
            )
        }
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
            let pair = pair.map_err(|error| wgs_paired_sampler_error(error, individual_id))?;
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
    state: UnitState,
    diagnostics: &[Diagnostic],
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
    // El catálogo del enlace en el manifiesto general es el de la amplificación
    // cuando la hay: v6 con tartamudeo clonal, v5 sin él. Hasta el 12/10/2026
    // aquí se escribía siempre el v5 (DT-040).
    let pcr_catalog = pcr.as_ref().map(|pcr| pcr.settings.rng_catalog());
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
        state,
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
        diagnostics,
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
        state,
        manifest_artifact_id: artifact.id,
        chimera_rate: event.chimera_rate(),
        realized_chimera_pairs,
        profile_id,
        profile_sha256,
        rng_catalog: if let Some(catalog) = pcr_catalog {
            catalog
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
    prepared_population_snvs: Option<PreparedPopulationSnvsForExecution>,
) -> Result<PreparedExecution, Diagnostic> {
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

/// Lo que un paquete fallido conserva de la unidad genotípica (DT-038).
struct FailedGenotypeUnit<'a> {
    /// La verdad, la matriz de SNV y el VCF con su índice, ya comprobados en
    /// disco. Las referencias viajan aparte, con las demás.
    artifact_ids: Vec<String>,
    population_snvs: Option<&'a PopulationSnvsManifestLink>,
    vcf_outputs: &'a [VcfManifestLink],
    repeat_draws: &'a [crate::repeat_genotypes::RepeatAlleleDraw],
    realized_tract_output: &'a crate::results::RealizedTractOutput,
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
    genotype_unit: Option<&FailedGenotypeUnit<'_>>,
    diagnostics: &[Diagnostic],
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
    retained.extend(
        references
            .iter()
            .map(|reference| reference.fasta_artifact_id.clone()),
    );
    retained.extend(
        reductions
            .iter()
            .filter_map(|reduction| reduction.map_artifact_id.clone()),
    );
    if events.is_some() {
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
    if let Some((inserted, _)) = insertion {
        retained.insert(inserted.map.id.clone());
    }
    if let Some(unit) = genotype_unit {
        retained.extend(unit.artifact_ids.iter().cloned());
    }
    workspace.retain_only_artifacts(&retained)?;
    let artifacts = workspace.artifacts().to_vec();
    let configuration_input = configuration_artifact.map(|copy| ConfigurationInputRecord {
        kind: "general".to_owned(),
        original_path: effective.configuration_input().original_path().to_owned(),
        copy,
    });
    // Sin unidad genotípica el paquete fallido no publica tracto alguno: sus
    // cifras realizadas son cero y la cota sigue siendo la que la
    // configuración fijó. Con ella, son los sorteos y las bases de tracto que
    // la verdad y el VCF conservados llevan de verdad.
    let no_tract_output = crate::results::RealizedTractOutput::default();
    let repeat_loci = crate::results::repeat_loci_manifest_json(
        effective,
        genotype_unit.map_or(&[][..], |unit| unit.repeat_draws),
        genotype_unit.map_or(&no_tract_output, |unit| unit.realized_tract_output),
    )?;
    let manifest = GeneralManifest {
        execution_id,
        state: UnitState::Failed,
        configuration_input: configuration_input.as_ref(),
        effective_config: None,
        references,
        reduced_references: reductions,
        population_snvs: genotype_unit.and_then(|unit| unit.population_snvs),
        vcf_outputs: genotype_unit.map_or(&[][..], |unit| unit.vcf_outputs),
        events,
        authorizations,
        artifacts: &artifacts,
        diagnostics,
        // Sin unidad genotípica no se publica ningún sorteo, pero los loci
        // **sí están declarados**, y el manifiesto los registra igual que en
        // un éxito: con cero copias sorteadas y sin frecuencia realizada. Lo
        // que se registra es la configuración, que existe tanto si la
        // ejecución terminó como si no.
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
            let matrix = generate_population_snvs_from_reference_excluding(
                parameters,
                reference,
                effective.reproducibility().root_seed(),
                effective.normalization().tolerance(),
                &repeat_catalog_exclusion(effective),
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

/// Los loci del catálogo, por secuencia, como tramos que el generador de la
/// matriz no puede elegir (REQ-CONF-143). Sin loci, vacío: el generador hace
/// exactamente lo de siempre.
fn repeat_catalog_exclusion(
    effective: &EffectiveConfig,
) -> crate::population_snvs::ExcludedIntervals {
    let mut excluded = crate::population_snvs::ExcludedIntervals::new();
    for locus in effective.repeat_loci() {
        let interval = locus.declaration().interval();
        excluded
            .entry(interval.sequence_id().as_str().to_owned())
            .or_default()
            .push((interval.start(), interval.end()));
    }
    for intervals in excluded.values_mut() {
        intervals.sort_unstable();
    }
    excluded
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

/// Escribe la tabla y el metadato del catálogo de microsatélites, con la
/// huella del FASTA final y la de la tabla dentro del metadato, y devuelve el
/// enlace con huella que la efectiva publica en `repeat_detection.catalog`.
/// Sin la sección de microsatélites no escribe nada.
fn publish_repeat_catalog(
    effective: &EffectiveConfig,
    reference_fasta: &ArtifactRecord,
    workspace: &mut ExecutionWorkspace,
) -> Result<Option<serde_json::Value>, Diagnostic> {
    if effective.repeat_detection().is_none() {
        return Ok(None);
    }
    let table = workspace.write_artifact(
        "repeat-catalog-tsv",
        "repeat_catalog_tsv",
        Path::new(crate::results::REPEAT_CATALOG_TSV_PATH),
        &crate::results::repeat_catalog_tsv_bytes(effective)?,
    )?;
    let metadata = crate::results::repeat_catalog_metadata_document(
        effective,
        &reference_fasta.digest,
        &table.digest,
    )?;
    let metadata = workspace.write_artifact(
        "repeat-catalog-metadata",
        "repeat_catalog_metadata",
        Path::new(crate::results::REPEAT_CATALOG_METADATA_PATH),
        &json_document_bytes(&metadata)?,
    )?;
    Ok(Some(crate::results::repeat_catalog_link_json(
        &metadata.digest,
    )))
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
        || effective.repeat_detection().is_some()
}

/// Prepara la referencia final. La guarda de admisión de los loci declarados
/// ya no corre aquí sino después de detectar (`resolve_repeat_catalog`): con
/// la detección, un declarado que no coincide con ningún tramo tiene que
/// rechazarse con el diagnóstico de la coincidencia (REQ-REF-40), que nombra
/// lo detectado más cercano y la etapa que no supera.
fn prepare_reference(effective: &EffectiveConfig) -> Result<PreparedReference, Diagnostic> {
    prepare_reference_without_repeat_admission(effective)
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
    // **Índices construidos una vez**, en vez de comparar cada locus con cada
    // cambio directo y con cada fila de la matriz: con miles de loci y una
    // matriz grande, eso crecía con el producto de los dos (paso F6.4.a.8h
    // del plan de detección de microsatélites). Cada índice guarda, por
    // posición, **el primero en su orden** —de declaración o de la matriz—,
    // y de los que caen dentro de un locus se nombra el primero en ese orden:
    // el mismo que encontraba el recorrido de antes.
    //
    // Una posición en base 1 igual a cero no cae dentro de ningún intervalo
    // (su fin, en base 0, es el cero); por eso no entra en los índices.
    let direct_events = effective.direct_events();
    let mut events_by_id = BTreeMap::<&str, usize>::new();
    let mut snp_events_by_position = BTreeMap::<(&str, u64), (usize, &SequenceId, u64)>::new();
    for (index, event) in direct_events.iter().enumerate() {
        events_by_id
            .entry(event.event_id().as_str())
            .or_insert(index);
        // Solo el cambio de base declara su posición en coordenadas de la
        // referencia. La inserción y la deleción la declaran en las
        // vigentes, y compararlas aquí se equivoca en las dos direcciones:
        // las juzga la cadena de segmentos, que no tiene deriva posible.
        if let EffectiveDirectEvent::Snp {
            sequence_id,
            position_1based,
            ..
        } = event
            && *position_1based > 0
        {
            snp_events_by_position
                .entry((sequence_id.as_str(), position_1based - 1))
                .or_insert((index, sequence_id, *position_1based));
        }
    }
    let rows = matrix.map(|matrix| matrix.rows()).unwrap_or_default();
    let mut rows_by_position = BTreeMap::<(&str, u64), usize>::new();
    for (index, row) in rows.iter().enumerate() {
        if row.position_1based() > 0 {
            rows_by_position
                .entry((row.contig().as_str(), row.position_1based() - 1))
                .or_insert(index);
        }
    }
    for locus in effective.repeat_loci() {
        let interval = locus.declaration().interval();
        let sequence_id = interval.sequence_id().as_str();
        let (start, end) = (interval.start(), interval.end());
        // **El identificador del evento que coloca un alelo es el del propio
        // locus**, de modo que un cambio directo llamado igual haría ambiguo
        // el origen de la modificación en la verdad. Se rechaza aquí, sobre
        // lo que el usuario escribió, en vez de morir al final con un error
        // de formato sobre un documento producido.
        //
        // Los dos rechazos por un cambio directo se comprueban juntos y gana
        // el del cambio declarado antes, como cuando se recorrían en orden.
        let named_like_the_locus = events_by_id.get(locus.id()).copied();
        let invading = first_inside(&snp_events_by_position, sequence_id, start, end);
        match (named_like_the_locus, invading) {
            (Some(named), invading) if invading.is_none_or(|(invading, ..)| named <= invading) => {
                return Err(repeat_admission_error(
                    locus.id(),
                    format!(
                        "el cambio directo {} lleva el mismo identificador que el locus, y el locus lo usa como identificador del evento que coloca su alelo",
                        direct_events[named].event_id()
                    ),
                ));
            }
            (_, Some((invading, event_sequence, position_1based))) => {
                return Err(repeat_invasion_error(
                    locus,
                    &format!(
                        "el cambio de base {} en {event_sequence}:{position_1based}",
                        direct_events[invading].event_id()
                    ),
                ));
            }
            (_, None) => {}
        }
        if let Some(invading) = first_inside(&rows_by_position, sequence_id, start, end) {
            let row = &rows[invading];
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
    Ok(())
}

/// El primero, en su orden, de los que un índice por posición tiene dentro de
/// `[start, end)` de la secuencia.
fn first_inside<'index, V: Copy + Ord>(
    index: &BTreeMap<(&'index str, u64), V>,
    sequence_id: &'index str,
    start: u64,
    end: u64,
) -> Option<V> {
    index
        .range((sequence_id, start)..(sequence_id, end))
        .map(|(_, value)| *value)
        .min()
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
            crate::test_detection_section::with_detection_section(
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
            ),
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

    /// La configuración efectiva inválida se retira, pero la verdad y la
    /// referencia que cita, ya validadas y sin tocar, se conservan como unidad
    /// genotípica (DT-038): este escenario no pide VCF ni matriz de SNV.
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
            destination
                .join("resources/references/ref-main.fa")
                .is_file()
        );
        assert!(destination.join("outputs/truth.jsonl").is_file());
        assert_eq!(
            manifest["references"],
            json!([{"id": "ref-main", "fasta_artifact_id": "reference-fasta"}])
        );
        assert_eq!(manifest["reference_reductions"], json!([]));
        assert_eq!(
            manifest["output_categories"]["minimum"],
            json!(["config-general", "reference-fasta", "truth-jsonl"])
        );
        assert_eq!(manifest["output_categories"]["requested"], json!([]));
        assert_eq!(
            catalog_ids(&manifest),
            ["config-general", "reference-fasta", "truth-jsonl"]
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// La reducción alterada queda `failed` y sin artefactos; la verdad y la
    /// referencia suministrada, intactas, se conservan como unidad genotípica
    /// (DT-038).
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
            "resources/references/reference-small.fa",
            "resources/references/reference-small.map.json",
        ] {
            assert!(!destination.join(relative).exists(), "retained {relative}");
        }
        assert!(
            destination
                .join("resources/references/reference-main.fa")
                .is_file()
        );
        assert!(destination.join("outputs/truth.jsonl").is_file());
        assert_eq!(
            manifest["references"],
            json!([{"id": "reference-main", "fasta_artifact_id": "reference-fasta"}])
        );
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
                "minimum": ["config-general", "reference-fasta", "truth-jsonl"],
                "defaulted": [],
                "requested": [],
                "disabled": []
            })
        );
        assert_eq!(
            catalog_ids(&manifest),
            ["config-general", "reference-fasta", "truth-jsonl"]
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// Configuración sobre una referencia suministrada de 2000 pb con **dos**
    /// referencias reducidas derivadas (`ref-small` y `ref-small-b`) y, si
    /// `individuals` no está vacío, un evento técnico de la vía de parejas.
    fn reductions_config(root: &Path, individuals: &[(&str, &str, &str)]) -> PathBuf {
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
        let (population_ids, event) = if individuals.is_empty() {
            ("\"I1\"".to_owned(), String::new())
        } else {
            let ids = individuals
                .iter()
                .map(|(id, _, _)| format!("\"{id}\""))
                .collect::<Vec<_>>()
                .join(", ");
            let libraries = individuals
                .iter()
                .map(|(id, i7, i5)| {
                    format!(
                        "[[technical_event.libraries]]\nindividual_id = \"{id}\"\ni7 = \"{i7}\"\ni5 = \"{i5}\"\n"
                    )
                })
                .collect::<String>();
            (
                ids.clone(),
                format!(
                    "\n[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [{ids}]\ntarget_coverage_wgs = 2\nsequencing_errors = \"disabled\"\n{libraries}"
                ),
            )
        };
        let count = individuals.len().max(1);
        let config = root.join("general.toml");
        fs::write(
            &config,
            crate::test_detection_section::with_detection_section(&format!(
                r#"config_kind = "general"
schema_version = "4"
seed = "00000000000000000000000000000000000000000000000000000000000000bb"

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

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 2000

[[reduced_references]]
id = "ref-small"
source_reference_id = "ref-main"
region_size_bases = 100
target_bases = 1000

[[reduced_references]]
id = "ref-small-b"
source_reference_id = "ref-main"
region_size_bases = 100
target_bases = 800

[[initial_populations]]
population_id = "P1"
individual_count = {count}
individual_ids = [{population_ids}]
ploidy = 2
genotype_model = "hwe"
{event}"#
            )),
        )
        .expect("config");
        config
    }

    #[derive(Clone, Copy)]
    enum FailureScenario {
        /// El gancho previo al manifiesto falla con un error de E/S propio.
        HookFails,
        /// Una base del FASTA de `ref-small` cambia **después** de validarse.
        TamperedFasta,
        /// El mapa de `ref-small` gana un campo ajeno después de validarse.
        TamperedMap,
    }

    /// Ejecuta `reductions_config` con o sin evento hasta un fallo controlado y
    /// devuelve el error, el manifiesto fallido y el destino.
    fn run_failed_with_reductions(
        with_event: bool,
        scenario: FailureScenario,
    ) -> (Diagnostic, serde_json::Value, PathBuf, PathBuf) {
        let root = std::env::temp_dir().join(format!(
            "gps-failure-reductions-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("test root");
        let individuals: &[(&str, &str, &str)] = if with_event {
            &[("I1", "ACGTACGT", "TGCATGCA")]
        } else {
            &[]
        };
        let config = reductions_config(&root, individuals);
        let effective = load_and_validate_config(&config).expect("valid config");
        let destination = root.join("result");
        let error = execute_validated_config_with_effective_document(
            &effective,
            &destination,
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            |workspace| match scenario {
                FailureScenario::HookFails => Err(Diagnostic::io(
                    "GPS-IO-001",
                    "Una operación de lectura o escritura falló.",
                )
                .with("detail", "fallo provocado tras las reducciones")),
                FailureScenario::TamperedFasta => {
                    let artifact = workspace
                        .artifacts()
                        .iter()
                        .find(|artifact| artifact.id == "reference-fasta-ref-small")
                        .expect("derived reference artifact")
                        .clone();
                    let mut bytes =
                        fs::read(workspace.root().join(&artifact.path)).expect("derived FASTA");
                    let start = bytes
                        .iter()
                        .position(|byte| *byte == b'\n')
                        .expect("header")
                        + 1;
                    bytes[start] = if bytes[start] == b'A' { b'C' } else { b'A' };
                    workspace.replace_artifact_bytes_for_validation_test(
                        "reference-fasta-ref-small",
                        &bytes,
                    )
                }
                FailureScenario::TamperedMap => {
                    let artifact = workspace
                        .artifacts()
                        .iter()
                        .find(|artifact| artifact.id == "reference-map-ref-small")
                        .expect("derived map artifact")
                        .clone();
                    let mut document: serde_json::Value = serde_json::from_slice(
                        &fs::read(workspace.root().join(&artifact.path)).expect("derived map"),
                    )
                    .expect("map JSON");
                    document["foreign_field"] = json!(1);
                    workspace.replace_artifact_bytes_for_validation_test(
                        "reference-map-ref-small",
                        &serde_json::to_vec(&document).expect("map bytes"),
                    )
                }
            },
        )
        .expect_err("la ejecución falla");
        let manifest: serde_json::Value = serde_json::from_slice(
            &fs::read(destination.join("manifest.json")).expect("manifiesto fallido"),
        )
        .expect("JSON");
        (error, manifest, destination, root)
    }

    fn assert_reduction_kept(manifest: &serde_json::Value, destination: &Path, id: &str) {
        let reduction = manifest["reference_reductions"]
            .as_array()
            .expect("reductions")
            .iter()
            .find(|reduction| reduction["id"] == id)
            .expect("reduction");
        assert_eq!(reduction["state"], "complete", "{id}");
        assert_eq!(
            reduction["fasta_artifact_id"],
            format!("reference-fasta-{id}")
        );
        assert_eq!(reduction["map_artifact_id"], format!("reference-map-{id}"));
        assert_eq!(reduction["source_fasta_artifact_id"], "reference-fasta");
        assert!(
            manifest["references"]
                .as_array()
                .expect("references")
                .contains(&json!({"id": id, "fasta_artifact_id": format!("reference-fasta-{id}")})),
            "{id} in {}",
            manifest["references"]
        );
        for relative in [
            format!("resources/references/{id}.fa"),
            format!("resources/references/{id}.map.json"),
        ] {
            assert!(destination.join(&relative).is_file(), "missing {relative}");
        }
    }

    fn assert_reduction_discarded(manifest: &serde_json::Value, destination: &Path, id: &str) {
        let reduction = manifest["reference_reductions"]
            .as_array()
            .expect("reductions")
            .iter()
            .find(|reduction| reduction["id"] == id)
            .expect("reduction");
        assert_eq!(reduction["state"], "failed", "{id}");
        assert_eq!(reduction["fasta_artifact_id"], serde_json::Value::Null);
        assert_eq!(reduction["map_artifact_id"], serde_json::Value::Null);
        assert_eq!(
            reduction["source_fasta_artifact_id"],
            serde_json::Value::Null
        );
        assert!(
            manifest["references"]
                .as_array()
                .expect("references")
                .iter()
                .all(|reference| reference["id"] != id),
            "{id} in {}",
            manifest["references"]
        );
        for relative in [
            format!("resources/references/{id}.fa"),
            format!("resources/references/{id}.map.json"),
        ] {
            assert!(!destination.join(&relative).exists(), "retained {relative}");
        }
    }

    /// REQ-SAL-25 y DTC-045. Un fallo **después** de publicar el evento no
    /// invalida lo ya publicado y validado que sea independiente: el paquete
    /// fallido conserva la referencia suministrada, el manifiesto del evento y
    /// **las dos** referencias reducidas completas (FASTA y mapa en disco,
    /// enlaces y estado `complete`), y registra la ejecución como fallida con su
    /// causa. Las aserciones sobre la lista de artefactos y sobre la presencia
    /// en disco **no son independientes del juez del catálogo**
    /// (`require_artifact_kind`, `verify_catalog`): sacar la reducción de
    /// `references` o no retener su mapa cae antes por `GPS-FMT-006` en la
    /// primera aserción. Las que sí cazan por sí solas son las del estado, la
    /// lista de referencias y la de que se conservan las dos.
    #[test]
    fn a_failure_after_the_event_keeps_the_completed_reduced_references() {
        let (error, manifest, destination, root) =
            run_failed_with_reductions(true, FailureScenario::HookFails);
        assert_eq!(error.code(), "GPS-IO-001");
        assert!(format!("{error:?}").contains("fallo provocado"));
        assert_eq!(manifest["state"], "failed");
        // Al manifiesto solo llegan el código y el mensaje de catálogo, así que
        // una causa distinta se detecta por un código distinto del genérico.
        assert_eq!(manifest["diagnostics"][0]["code"], "GPS-IO-001");
        assert_eq!(manifest["diagnostics"].as_array().expect("array").len(), 1);
        assert_eq!(manifest["events"][0]["event_id"], "reads-main");
        assert_eq!(
            manifest["references"],
            json!([
                {"id": "ref-main", "fasta_artifact_id": "reference-fasta"},
                {"id": "ref-small", "fasta_artifact_id": "reference-fasta-ref-small"},
                {"id": "ref-small-b", "fasta_artifact_id": "reference-fasta-ref-small-b"}
            ])
        );
        assert_reduction_kept(&manifest, &destination, "ref-small");
        assert_reduction_kept(&manifest, &destination, "ref-small-b");
        assert!(
            destination
                .join("events/reads-main/manifest.json")
                .is_file()
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// REQ-SAL-25 y DTC-045 sin evento: la misma regla. Una reducción completa
    /// se conserva aunque el fallo llegue antes de cualquier lectura, junto con
    /// la referencia suministrada de la que deriva.
    #[test]
    fn a_failure_without_event_keeps_the_completed_reduced_references_and_their_source() {
        let (error, manifest, destination, root) =
            run_failed_with_reductions(false, FailureScenario::HookFails);
        assert_eq!(error.code(), "GPS-IO-001");
        assert_eq!(manifest["state"], "failed");
        assert_eq!(manifest["events"], serde_json::Value::Null);
        assert_eq!(
            manifest["references"],
            json!([
                {"id": "ref-main", "fasta_artifact_id": "reference-fasta"},
                {"id": "ref-small", "fasta_artifact_id": "reference-fasta-ref-small"},
                {"id": "ref-small-b", "fasta_artifact_id": "reference-fasta-ref-small-b"}
            ])
        );
        assert_reduction_kept(&manifest, &destination, "ref-small");
        assert_reduction_kept(&manifest, &destination, "ref-small-b");
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// REQ-SAL-27. Una reducción cuyos bytes cambian **después** de validarse
    /// —el juez final la rechaza— no se publica como completa en el paquete
    /// fallido, con evento o sin él, ni cambiando el FASTA ni el mapa: queda
    /// `failed`, sin enlaces y sin ficheros, y la otra reducción, intacta, se
    /// conserva.
    #[test]
    fn a_reduction_altered_after_validation_is_discarded_and_the_other_is_kept() {
        for with_event in [true, false] {
            for scenario in [FailureScenario::TamperedFasta, FailureScenario::TamperedMap] {
                let (error, manifest, destination, root) =
                    run_failed_with_reductions(with_event, scenario);
                assert_eq!(error.code(), "GPS-FMT-006");
                assert_eq!(manifest["state"], "failed");
                assert_reduction_discarded(&manifest, &destination, "ref-small");
                assert_reduction_kept(&manifest, &destination, "ref-small-b");
                fs::remove_dir_all(root).expect("cleanup");
            }
        }
    }

    /// REQ-SAL-21, REQ-SAL-22, REQ-SAL-23 y REQ-SAL-25 (DT-036). Si el
    /// ensamblado del evento falla en la **segunda** biblioteca —la validación
    /// conjunta de sus mates—, la primera ya cerró sus FASTQ y aun así el
    /// paquete fallido no conserva ninguno. El evento **queda enlazado** con
    /// estado `failed` y un manifiesto de evento fallido: la causa en sus
    /// diagnósticos y las cuatro salidas FASTQ `failed`, sin artefacto y con el
    /// código de la causa. Se conservan además las reducciones completas, su
    /// referencia de origen y la unidad genotípica, porque este escenario no
    /// pide VCF. El manifiesto general solo guarda el código y el mensaje de
    /// catálogo de la causa; la que distingue causas entre sí es la prueba del
    /// fallo tras el evento, con un código propio.
    #[test]
    fn a_failed_event_assembly_links_a_failed_event_manifest_and_keeps_no_fastq() {
        use crate::fastq_output::FORCE_FINAL_PAIR_VALIDATION_FAILURE;
        let root = test_root("event-assembly");
        let config = reductions_config(
            &root,
            &[
                ("I1", "ACGTACGT", "TGCATGCA"),
                ("I2", "TGCATGCA", "ACGTACGT"),
            ],
        );
        let destination = root.join("result");
        FORCE_FINAL_PAIR_VALIDATION_FAILURE.set(2);
        let (error, manifest) = run_failed(&config, &destination, |_| Ok(()));
        assert_eq!(
            FORCE_FINAL_PAIR_VALIDATION_FAILURE.get(),
            0,
            "the injected failure must have fired on the second library"
        );
        assert!(
            format!("{error:?}").contains("injected final-pair validation failure"),
            "{error:?}"
        );
        assert_eq!(manifest["state"], "failed");
        assert_eq!(diagnostic_codes(&manifest), [error.code()]);
        assert_eq!(
            manifest["events"],
            json!([{
                "event_id": "reads-main",
                "state": "failed",
                "manifest_artifact_id": "event-manifest-reads-main"
            }])
        );
        assert_eq!(
            manifest["technical_events"][0]["event_manifest"]["id"],
            "event-manifest-reads-main"
        );
        assert_eq!(
            manifest["references"],
            json!([
                {"id": "ref-main", "fasta_artifact_id": "reference-fasta"},
                {"id": "ref-small", "fasta_artifact_id": "reference-fasta-ref-small"},
                {"id": "ref-small-b", "fasta_artifact_id": "reference-fasta-ref-small-b"}
            ])
        );
        assert_reduction_kept(&manifest, &destination, "ref-small");
        assert_reduction_kept(&manifest, &destination, "ref-small-b");
        assert_eq!(
            catalog_ids(&manifest),
            [
                "config-general",
                "event-manifest-reads-main",
                "population-snvs-metadata",
                "population-snvs-tsv",
                "reference-fasta",
                "reference-fasta-ref-small",
                "reference-fasta-ref-small-b",
                "reference-map-ref-small",
                "reference-map-ref-small-b",
                "truth-jsonl"
            ]
        );
        assert_failed_event_manifest(
            &read_json(&destination.join("events/reads-main/manifest.json")),
            error.code(),
            2,
        );
        assert_no_fastq(&destination);
        assert_external_schema(
            &destination.join("manifest.json"),
            "manifest-general-v4.schema.json",
        );
        assert_external_schema(
            &destination.join("events/reads-main/manifest.json"),
            "manifest-event-v2.schema.json",
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// Un manifiesto de evento fallido según DT-036: estado `failed`, la causa
    /// como único diagnóstico, las 2n salidas `failed` sin artefacto y con su
    /// código, cobertura y resúmenes vacíos y todas las cifras realizadas a
    /// cero.
    fn assert_failed_event_manifest(event: &serde_json::Value, code: &str, libraries: usize) {
        assert_eq!(event["state"], "failed");
        assert_eq!(diagnostic_codes(event), [code]);
        let outputs = event["outputs"].as_array().expect("outputs");
        assert_eq!(outputs.len(), 2 * libraries, "{outputs:?}");
        for output in outputs {
            assert_eq!(output["state"], "failed", "{output}");
            assert_eq!(output["artifact"], serde_json::Value::Null, "{output}");
            assert_eq!(output["diagnostic_code"], code, "{output}");
        }
        assert_eq!(event["coverage"], json!([]));
        assert_eq!(event["final_fastq_summary"], json!([]));
        assert_eq!(event["artifacts"]["chimera"]["realized_pairs"], 0);
        assert_eq!(event["artifacts"]["biological_contamination"], json!([]));
        for field in ["reassigned_pairs", "discarded_pairs", "unassigned_pairs"] {
            assert_eq!(event["demultiplexing"][field], 0, "{field}");
        }
    }

    fn assert_no_fastq(destination: &Path) {
        for relative in files_under(destination) {
            assert!(!relative.contains("fastq"), "retained {relative}");
        }
    }

    /// Los cuatro escenarios de evento de `tests/event_manifest_guards.rs`
    /// —contaminación con salto de índice, amplificación, amplificación con
    /// tartamudeo clonal y tartamudeo sin amplificar—, con el fallo inyectado
    /// en la última biblioteca. Devuelve la configuración, el número de
    /// bibliotecas y los esquemas del manifiesto general y del de evento.
    fn failed_event_scenario(
        root: &Path,
        scenario: &str,
    ) -> (PathBuf, usize, &'static str, &'static str) {
        let short_and_long = |short: usize, long: usize| {
            let short = "CA".repeat(short);
            let long = "CA".repeat(long);
            let head = format!("{}{short}GGGG{long}", "G".repeat(60));
            (
                format!("{head}{}", "G".repeat(300 - head.len())),
                format!(
                    "[[reference.repeat_loci]]\nid = \"STR-SHORT\"\nsequence = \"chr1\"\nstart = 61\nend = {}\nmotif = \"CA\"\n[[reference.repeat_loci.alleles]]\nsequence = \"{short}\"\nfrequency = 1.0\n[[reference.repeat_loci]]\nid = \"STR-LONG\"\nsequence = \"chr1\"\nstart = {}\nend = {}\nmotif = \"CA\"\n[[reference.repeat_loci.alleles]]\nsequence = \"{long}\"\nfrequency = 1.0\n",
                    60 + short.len(),
                    65 + short.len(),
                    64 + short.len() + long.len()
                ),
            )
        };
        let stutter = |rate: f64| {
            format!(
                "[technical_event.stutter]\nenabled = true\n[[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = {rate:?}\n[[technical_event.stutter.steps]]\nunits = -3\nweight = 0.05\n[[technical_event.stutter.steps]]\nunits = -2\nweight = 0.15\n[[technical_event.stutter.steps]]\nunits = -1\nweight = 0.6\n[[technical_event.stutter.steps]]\nunits = 1\nweight = 0.2\n"
            )
        };
        let pcr = "library_preparation = \"pcr\"\n[technical_event.pcr]\nduplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n";
        let library_a = "[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n";
        let (bases, loci, individuals, event, libraries, schemas) = match scenario {
            "contamination" => (
                "A".repeat(1000),
                String::new(),
                "individual_count = 2\nploidy = 2\nindividual_ids = [\"sample-A\", \"sample-B\"]",
                format!(
                    "individual_ids = [\"sample-A\", \"sample-B\"]\ntarget_coverage_wgs = 1\ncontamination_min_percent = 50\ncontamination_max_percent = 50\ncontamination_recipient_ids = [\"sample-A\"]\nindex_hopping_rate = 1\n{library_a}[[technical_event.libraries]]\nindividual_id = \"sample-B\"\ni7 = \"CGTACGTA\"\ni5 = \"GCATGCAT\"\n"
                ),
                2,
                (
                    "manifest-general-v4.schema.json",
                    "manifest-event-v2.schema.json",
                ),
            ),
            "pcr" => (
                "ACGTTGCAAGCT".repeat(84)[..1000].to_owned(),
                String::new(),
                "individual_count = 1\nploidy = 1\nindividual_ids = [\"sample-A\"]",
                format!(
                    "individual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 4\nsequencing_errors = \"disabled\"\n{pcr}{library_a}"
                ),
                1,
                (
                    "manifest-general-v4.schema.json",
                    "manifest-event-v3.schema.json",
                ),
            ),
            "clonal" => {
                let (bases, loci) = short_and_long(10, 17);
                (
                    bases,
                    loci,
                    "individual_count = 1\nploidy = 1\nindividual_ids = [\"sample-A\"]",
                    format!(
                        "individual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 20\nsequencing_errors = \"disabled\"\n{pcr}{}{library_a}",
                        stutter(0.3)
                    ),
                    1,
                    (
                        "manifest-general-v6.schema.json",
                        "manifest-event-v3.schema.json",
                    ),
                )
            }
            "stutter" => {
                let (bases, loci) = short_and_long(5, 17);
                (
                    bases,
                    loci,
                    "individual_count = 1\nploidy = 1\nindividual_ids = [\"sample-A\"]",
                    format!(
                        "individual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 20\nsequencing_errors = \"disabled\"\n{}{library_a}",
                        stutter(1.0)
                    ),
                    1,
                    (
                        "manifest-general-v6.schema.json",
                        "manifest-event-v3.schema.json",
                    ),
                )
            }
            other => panic!("unknown scenario {other}"),
        };
        fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
        fs::write(
            root.join("population-snvs.tsv"),
            format!(
                "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\t{}\tT\tP1\t0.5\n",
                &bases[19..20]
            ),
        )
        .expect("SNV matrix");
        let config = root.join("general.toml");
        fs::write(
            &config,
            crate::test_detection_section::with_detection_section(&format!(
                r#"config_kind = "general"
schema_version = "4"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[population_snvs]
mode = "supplied"
input_path = "population-snvs.tsv"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {}
{loci}
[[initial_populations]]
population_id = "P1"
genotype_model = "hwe"
{individuals}

[technical_event]
event_id = "reads-main"
{event}"#,
                bases.len()
            )),
        )
        .expect("config");
        (config, libraries, schemas.0, schemas.1)
    }

    /// DT-036 con cada etapa que cambia la forma del manifiesto de evento: la
    /// contaminación con salto de índice (v2), la amplificación (v3, catálogo
    /// v5), la amplificación con tartamudeo clonal (v3, catálogo v6, manifiesto
    /// general v6) y el tartamudeo sin amplificar (v3, catálogo v4). El
    /// manifiesto de evento fallido, con sus secciones de etapa en cero, y el
    /// manifiesto general que lo enlaza pasan los **esquemas publicados**
    /// leídos por el validador externo, no solo los jueces de Rust.
    #[test]
    fn failed_event_manifests_of_every_stage_pass_the_published_schemas() {
        use crate::fastq_output::FORCE_FINAL_PAIR_VALIDATION_FAILURE;
        for scenario in ["contamination", "pcr", "clonal", "stutter"] {
            let root = test_root(&format!("event-stage-{scenario}"));
            let (config, libraries, general_schema, event_schema) =
                failed_event_scenario(&root, scenario);
            let destination = root.join("result");
            FORCE_FINAL_PAIR_VALIDATION_FAILURE.set(u32::try_from(libraries).expect("count"));
            let (error, manifest) = run_failed(&config, &destination, |_| Ok(()));
            assert_eq!(FORCE_FINAL_PAIR_VALIDATION_FAILURE.get(), 0, "{scenario}");
            assert_eq!(manifest["events"][0]["state"], "failed", "{scenario}");
            let event_path = destination.join("events/reads-main/manifest.json");
            let event = read_json(&event_path);
            assert_failed_event_manifest(&event, error.code(), libraries);
            match scenario {
                "pcr" | "clonal" => {
                    assert_eq!(event["artifacts"]["pcr"]["pairs"], 0, "{scenario}");
                    assert_eq!(event["artifacts"]["pcr"]["families"], 0, "{scenario}");
                }
                "contamination" => {
                    assert_eq!(event["demultiplexing"]["index_hopping_rate"], 1.0);
                }
                _ => {}
            }
            if matches!(scenario, "clonal" | "stutter") {
                assert_eq!(
                    event["artifacts"]["stutter"]["evaluated_tracts"], 0,
                    "{scenario}"
                );
            }
            assert_no_fastq(&destination);
            assert_external_schema(&destination.join("manifest.json"), general_schema);
            assert_external_schema(&event_path, event_schema);
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    /// DT-036 con un fallo **real**, sin inyección: un inserto impracticable
    /// para una copia de 255 pb (DT-006) agota la cota de intentos del
    /// muestreador y da `GPS-DOMAIN-003`. El evento queda enlazado y fallido
    /// con esa causa en su manifiesto y en cada salida.
    #[test]
    fn an_impractical_insert_leaves_the_event_linked_as_failed_with_its_cause() {
        let root = test_root("event-impractical-insert");
        let config = root.join("general.toml");
        fs::write(
            &config,
            crate::test_detection_section::with_detection_section(
                r#"config_kind = "general"
schema_version = "2"
seed = "0f0e0d0c0b0a09080706050403020100f0e0d0c0b0a090807060504030201000"
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
length = 255
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
ploidy = 1
genotype_model = "hwe"
[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 5
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
            ),
        )
        .expect("config");
        let destination = root.join("result");
        let (error, manifest) = run_failed(&config, &destination, |_| Ok(()));
        assert_eq!(error.code(), "GPS-DOMAIN-003");
        assert_eq!(diagnostic_codes(&manifest), ["GPS-DOMAIN-003"]);
        assert_eq!(manifest["events"][0]["state"], "failed");
        let event_path = destination.join("events/reads-main/manifest.json");
        assert_failed_event_manifest(&read_json(&event_path), "GPS-DOMAIN-003", 1);
        assert_no_fastq(&destination);
        assert_external_schema(
            &destination.join("manifest.json"),
            "manifest-general-v4.schema.json",
        );
        assert_external_schema(&event_path, "manifest-event-v2.schema.json");
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

    fn test_root(label: &str) -> PathBuf {
        let root = std::env::temp_dir().join(format!(
            "gps-failed-package-{label}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("test root");
        root
    }

    fn read_json(path: &Path) -> serde_json::Value {
        serde_json::from_slice(&fs::read(path).unwrap_or_else(|error| panic!("{path:?}: {error}")))
            .expect("JSON")
    }

    /// Identificadores del catálogo de un manifiesto general, ordenados.
    fn catalog_ids(manifest: &serde_json::Value) -> Vec<String> {
        let mut ids = manifest["artifacts"]
            .as_array()
            .expect("catalog")
            .iter()
            .map(|artifact| artifact["id"].as_str().expect("id").to_owned())
            .collect::<Vec<_>>();
        ids.sort();
        ids
    }

    fn diagnostic_codes(document: &serde_json::Value) -> Vec<String> {
        document["diagnostics"]
            .as_array()
            .expect("diagnostics")
            .iter()
            .map(|diagnostic| diagnostic["code"].as_str().expect("code").to_owned())
            .collect()
    }

    /// Rutas relativas de todos los ficheros bajo `root`, ordenadas.
    fn files_under(root: &Path) -> Vec<String> {
        let mut files = Vec::new();
        let mut pending = vec![root.to_path_buf()];
        while let Some(directory) = pending.pop() {
            for entry in fs::read_dir(&directory).expect("read directory") {
                let path = entry.expect("entry").path();
                if path.is_dir() {
                    pending.push(path);
                } else {
                    files.push(
                        path.strip_prefix(root)
                            .expect("relative")
                            .to_string_lossy()
                            .into_owned(),
                    );
                }
            }
        }
        files.sort();
        files
    }

    /// Valida un documento publicado con el validador externo de JSON Schema
    /// (`jsonschema` de Python), un proceso ajeno al productor de Rust.
    fn assert_external_schema(path: &Path, schema: &str) {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let output = std::process::Command::new("python3")
            .arg(root.join("tests/external_jsonschema_validator.py"))
            .arg(root.join("schemas").join(schema))
            .arg(path)
            .output()
            .expect("external JSON Schema validator");
        assert!(
            output.status.success(),
            "external validator for {path:?} against {schema}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    fn provoked_failure() -> Diagnostic {
        Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
            .with("detail", "fallo provocado antes del manifiesto")
    }

    /// Ejecuta `config` hasta un fallo controlado, con `before_manifest` como
    /// gancho previo al manifiesto, y devuelve el error y el manifiesto
    /// fallido publicado en `destination`.
    fn run_failed(
        config: &Path,
        destination: &Path,
        before_manifest: impl FnOnce(&mut ExecutionWorkspace) -> Result<(), Diagnostic>,
    ) -> (Diagnostic, serde_json::Value) {
        let effective = load_and_validate_config(config).expect("valid config");
        let error = execute_validated_config_with_effective_document(
            &effective,
            destination,
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            before_manifest,
        )
        .expect_err("the run fails");
        (error, read_json(&destination.join("manifest.json")))
    }

    /// Ejecuta `config` sin fallo y devuelve el manifiesto completo.
    fn run_complete(config: &Path, destination: &Path) -> serde_json::Value {
        let effective = load_and_validate_config(config).expect("valid config");
        execute_validated_config(&effective, destination).expect("complete run");
        read_json(&destination.join("manifest.json"))
    }

    fn assert_same_bytes(left: &Path, right: &Path, relative: &str) {
        assert_eq!(
            fs::read(left.join(relative)).unwrap_or_else(|error| panic!("{relative}: {error}")),
            fs::read(right.join(relative)).unwrap_or_else(|error| panic!("{relative}: {error}")),
            "{relative}"
        );
    }

    /// Referencia suministrada de 2000 pb, una matriz de SNV generada y dos
    /// individuos diploides; con un evento técnico sobre los dos si se pide y
    /// con VCF en el modo pedido.
    fn genotype_unit_config(root: &Path, with_event: bool, vcf: Option<&str>) -> PathBuf {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let bases = (0..2000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                char::from(b"ACGT"[(state % 4) as usize])
            })
            .collect::<String>();
        fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
        let event = if with_event {
            "\n[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [\"I1\", \"I2\"]\ntarget_coverage_wgs = 2\nsequencing_errors = \"disabled\"\n[[technical_event.libraries]]\nindividual_id = \"I1\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n[[technical_event.libraries]]\nindividual_id = \"I2\"\ni7 = \"TGCATGCA\"\ni5 = \"ACGTACGT\"\n"
        } else {
            ""
        };
        let vcf = vcf.map_or_else(String::new, |mode| {
            format!("\n[vcf_output]\nmode = \"{mode}\"\n")
        });
        let config = root.join("general.toml");
        fs::write(
            &config,
            crate::test_detection_section::with_detection_section(&format!(
                r#"config_kind = "general"
schema_version = "4"
seed = "00000000000000000000000000000000000000000000000000000000000000cc"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 3
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 2000

[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"
{event}{vcf}"#
            )),
        )
        .expect("config");
        config
    }

    /// Lo que la unidad genotípica de `genotype_unit_config` con VCF BGZF
    /// publica, aparte de la referencia.
    const GENOTYPE_UNIT_FILES: [&str; 5] = [
        "outputs/truth.jsonl",
        "outputs/truth.vcf.gz",
        "outputs/truth.vcf.gz.tbi",
        "population_snvs.tsv",
        "population_snvs.metadata.json",
    ];

    /// DT-038, REQ-SAL-25 y DTC-045. Un fallo **después** de publicar el VCF
    /// conserva la unidad genotípica entera —la verdad, la matriz de SNV, el
    /// VCF con su índice y la referencia que citan—, con evento y sin él, y
    /// byte a byte la misma que publica una ejecución completa de la misma
    /// configuración. Es la mitad de la antigua
    /// `vcf_index_io_failure_publishes_only_the_controlled_failed_result`, que
    /// inyectaba este mismo fallo y exigía lo contrario, descartar el VCF.
    #[test]
    fn a_failure_after_the_vcf_keeps_the_whole_genotype_unit() {
        for with_event in [false, true] {
            let root = test_root("unit-kept");
            let config = genotype_unit_config(&root, with_event, Some("bgzf"));
            let complete = root.join("complete");
            run_complete(&config, &complete);
            let destination = root.join("result");
            let (error, manifest) = run_failed(&config, &destination, |_| Err(provoked_failure()));
            assert_eq!(error.code(), "GPS-IO-001");
            assert!(format!("{error:?}").contains("fallo provocado"));
            assert_eq!(manifest["state"], "failed");
            assert_eq!(diagnostic_codes(&manifest), ["GPS-IO-001"]);
            let mut expected = vec![
                "config-general",
                "population-snvs-metadata",
                "population-snvs-tsv",
                "reference-fasta",
                "truth-jsonl",
                "truth-vcf",
                "truth-vcf-index",
            ];
            if with_event {
                expected.push("event-manifest-reads-main");
                expected.sort_unstable();
            }
            let ids = catalog_ids(&manifest);
            assert_eq!(
                ids.iter()
                    .filter(|id| !id.starts_with("fastq-"))
                    .collect::<Vec<_>>(),
                expected,
                "with_event={with_event}"
            );
            assert_eq!(
                ids.iter().filter(|id| id.starts_with("fastq-")).count(),
                if with_event { 4 } else { 0 }
            );
            assert_eq!(
                manifest["references"],
                json!([{"id": "ref-main", "fasta_artifact_id": "reference-fasta"}])
            );
            assert_eq!(manifest["vcf_outputs"][0]["state"], "complete");
            assert_eq!(manifest["vcf_outputs"][0]["reference_id"], "ref-main");
            assert_eq!(manifest["vcf_outputs"][0]["vcf_artifact_id"], "truth-vcf");
            assert_eq!(
                manifest["vcf_outputs"][0]["index_artifact_id"],
                "truth-vcf-index"
            );
            assert_eq!(
                manifest["population_snvs"]["metadata_artifact_id"],
                "population-snvs-metadata"
            );
            assert_eq!(
                manifest["output_categories"]["requested"],
                json!(["truth-vcf", "truth-vcf-index"])
            );
            for relative in GENOTYPE_UNIT_FILES
                .into_iter()
                .chain(["resources/references/ref-main.fa"])
            {
                assert_same_bytes(&destination, &complete, relative);
            }
            assert_external_schema(
                &destination.join("manifest.json"),
                "manifest-general-v4.schema.json",
            );
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    /// DT-038. Un fallo **real** al escribir el índice, con el VCF ya escrito
    /// y validado, deja el VCF pedido sin publicar, y con él cae la unidad
    /// genotípica entera: ni VCF, ni índice, ni verdad, ni matriz. Es la otra
    /// mitad de la antigua
    /// `vcf_index_io_failure_publishes_only_the_controlled_failed_result`.
    ///
    /// **No puede fallar contra el código anterior a DT-038**, que no
    /// conservaba nunca nada de la unidad. Lo que discrimina es una unidad que
    /// se conservara sin mirar si el VCF pedido llegó a publicarse.
    #[test]
    fn a_real_vcf_index_failure_drops_the_whole_genotype_unit() {
        use crate::vcf_publication::FORCE_INDEX_WRITE_FAILURE;
        let root = test_root("unit-index");
        let config = genotype_unit_config(&root, false, Some("bgzf"));
        let destination = root.join("result");
        FORCE_INDEX_WRITE_FAILURE.set(true);
        let (error, manifest) = run_failed(&config, &destination, |_| Ok(()));
        assert!(
            !FORCE_INDEX_WRITE_FAILURE.get(),
            "the injected index failure must have fired"
        );
        assert_eq!(error.code(), "GPS-IO-001");
        assert_eq!(error.exit_code(), 5);
        assert!(format!("{error:?}").contains("injected index writer failure"));
        assert_eq!(manifest["state"], "failed");
        assert_eq!(diagnostic_codes(&manifest), ["GPS-IO-001"]);
        assert_eq!(manifest["vcf_outputs"], json!([]));
        assert_eq!(manifest["population_snvs"], serde_json::Value::Null);
        assert_eq!(manifest["references"], json!([]));
        assert_eq!(catalog_ids(&manifest), ["config-general"]);
        assert_eq!(
            files_under(&destination),
            ["inputs/config/general.toml", "manifest.json"]
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// Cambia el último byte de un artefacto ya escrito y **rehace su entrada
    /// del catálogo**, de modo que el catálogo del momento del fallo coincide
    /// con el disco y solo el registro tomado al validar el artefacto delata
    /// el cambio. Es el escenario saturado que hace decidir a la regla de
    /// conservación: con una alteración bruta, sin rehacer el catálogo, el
    /// primer paquete fallido no pasaría sus jueces y el último recurso de
    /// DT-039 quitaría la unidad de todos modos, decidiera lo que decidiera la
    /// regla.
    fn alter_after_validation(workspace: &mut ExecutionWorkspace, id: &str) {
        let artifact = workspace
            .artifacts()
            .iter()
            .find(|artifact| artifact.id == id)
            .unwrap_or_else(|| panic!("{id} absent"))
            .clone();
        let mut bytes = fs::read(workspace.root().join(&artifact.path)).expect("artifact bytes");
        let last = bytes.len() - 1;
        bytes[last] ^= 0x01;
        workspace
            .replace_artifact_bytes_for_validation_test(id, &bytes)
            .expect("replace artifact");
    }

    /// DT-038, la regla de conservación pieza a pieza, con su pareja. Cada
    /// artefacto de la unidad genotípica —la verdad, la matriz y sus
    /// metadatos, el VCF, su índice y la referencia suministrada—, cambiado
    /// **después** de validarse y con el catálogo rehecho, deja fuera la
    /// unidad entera. La pareja, sin alterar nada, la conserva; sin ella, un
    /// paquete fallido que no conservara nunca la unidad pasaría la prueba.
    #[test]
    fn each_genotype_unit_piece_altered_after_validation_drops_the_whole_unit() {
        for altered in [
            None,
            Some("truth-jsonl"),
            Some("population-snvs-tsv"),
            Some("population-snvs-metadata"),
            Some("truth-vcf"),
            Some("truth-vcf-index"),
            Some("reference-fasta"),
        ] {
            let root = test_root("unit-piece");
            let config = genotype_unit_config(&root, false, Some("bgzf"));
            let destination = root.join("result");
            let (error, manifest) = run_failed(&config, &destination, move |workspace| {
                if let Some(id) = altered {
                    alter_after_validation(workspace, id);
                }
                Err(provoked_failure())
            });
            assert_eq!(error.code(), "GPS-IO-001", "{altered:?}");
            assert_eq!(manifest["state"], "failed");
            let unit_files = GENOTYPE_UNIT_FILES
                .into_iter()
                .chain(["resources/references/ref-main.fa"]);
            if altered.is_some() {
                assert_eq!(catalog_ids(&manifest), ["config-general"], "{altered:?}");
                assert_eq!(manifest["vcf_outputs"], json!([]), "{altered:?}");
                assert_eq!(manifest["population_snvs"], serde_json::Value::Null);
                assert_eq!(manifest["references"], json!([]), "{altered:?}");
                for relative in unit_files {
                    assert!(
                        !destination.join(relative).exists(),
                        "{altered:?} retained {relative}"
                    );
                }
            } else {
                for relative in unit_files {
                    assert!(destination.join(relative).is_file(), "missing {relative}");
                }
            }
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    /// DT-038 con inserción, la misma regla y su pareja. La referencia final,
    /// su mapa, la suministrada de la que deriva o la verdad, cambiadas
    /// **después** de validarse y con el catálogo rehecho, dejan fuera la
    /// unidad entera, y con ella la sección `reference_insertion`.
    #[test]
    fn each_inserted_reference_piece_altered_after_validation_drops_the_whole_unit() {
        for altered in [
            None,
            Some("reference-fasta-ref-main-str"),
            Some("reference-map-ref-main-str"),
            Some("reference-fasta"),
            Some("truth-jsonl"),
        ] {
            let root = test_root("unit-insertion-piece");
            let config = insertion_unit_config(&root);
            let destination = root.join("result");
            let (error, manifest) = run_failed(&config, &destination, move |workspace| {
                if let Some(id) = altered {
                    alter_after_validation(workspace, id);
                }
                Err(provoked_failure())
            });
            assert_eq!(error.code(), "GPS-IO-001", "{altered:?}");
            if altered.is_some() {
                assert_eq!(catalog_ids(&manifest), ["config-general"], "{altered:?}");
                assert!(manifest.get("reference_insertion").is_none(), "{altered:?}");
                assert_eq!(
                    files_under(&destination),
                    ["inputs/config/general.toml", "manifest.json"],
                    "{altered:?}"
                );
            } else {
                assert_eq!(manifest["reference_insertion"]["id"], "ref-main-str");
                assert!(destination.join("outputs/truth.jsonl").is_file());
                assert!(
                    destination
                        .join("resources/references/ref-main-str.map.json")
                        .is_file()
                );
            }
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    /// DT-038 con su pareja. Con loci de repetición, la cabecera de la verdad
    /// cita **la copia de la configuración** por su huella (fuente `random`).
    /// Si esa copia cambia en disco antes del cierre, el paquete fallido no la
    /// conserva, y entonces tampoco la verdad que la cita ni el resto de la
    /// unidad: entera o nada. La pareja, con la copia intacta, conserva las
    /// dos.
    #[test]
    fn an_altered_configuration_copy_drops_the_genotype_unit_that_cites_it() {
        for tamper in [false, true] {
            let root = test_root("unit-configuration");
            let config = repeat_unit_config(&root, SINGLE_REPEAT_ALLELE, None);
            let destination = root.join("result");
            let (error, manifest) = run_failed(&config, &destination, move |workspace| {
                if tamper {
                    let path = workspace.root().join("inputs/config/general.toml");
                    let mut bytes = fs::read(&path).expect("configuration copy");
                    bytes.extend_from_slice(b"\n# tampered\n");
                    fs::write(&path, bytes).expect("tamper");
                }
                Err(provoked_failure())
            });
            assert_eq!(error.code(), "GPS-IO-001");
            assert_eq!(manifest["state"], "failed");
            if tamper {
                assert_eq!(manifest["configuration_input"], serde_json::Value::Null);
                assert!(catalog_ids(&manifest).is_empty(), "{manifest}");
                assert_eq!(files_under(&destination), ["manifest.json"]);
            } else {
                assert_eq!(
                    manifest["configuration_input"]["artifact_id"],
                    "config-general"
                );
                let truth = fs::read_to_string(destination.join("outputs/truth.jsonl"))
                    .expect("kept truth");
                let header = truth.lines().next().expect("header");
                let copy = manifest["artifacts"]
                    .as_array()
                    .expect("catalog")
                    .iter()
                    .find(|artifact| artifact["id"] == "config-general")
                    .expect("configuration copy");
                let digest = copy["checksum"]["digest"].as_str().expect("copy digest");
                assert!(
                    header.contains(digest),
                    "control: the kept truth cites the kept configuration copy"
                );
            }
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    /// El juez del manifiesto general exige la clase y la ruta de la verdad
    /// también en un paquete **fallido** que la conserva (DT-038), no solo en
    /// uno completo. Cada cambio se prueba solo, sobre un manifiesto fallido
    /// real que el juez acepta sin él.
    #[test]
    fn a_failed_manifest_with_a_misclassified_truth_is_rejected() {
        let root = test_root("judge-truth-class");
        let config = genotype_unit_config(&root, false, None);
        let destination = root.join("result");
        let (_, manifest) = run_failed(&config, &destination, |_| Err(provoked_failure()));
        assert_eq!(manifest["state"], "failed");
        let path = destination.join("manifest.json");
        let judge = |document: &serde_json::Value| {
            crate::results::parse_general_manifest_bytes(
                &path,
                &serde_json::to_vec_pretty(document).expect("bytes"),
            )
            .map(|_| ())
        };
        judge(&manifest).expect("control: the untouched failed manifest is accepted");
        for (field, value) in [
            ("kind", "population_snvs_tsv"),
            ("path", "outputs/other-truth.jsonl"),
        ] {
            let mut altered = manifest.clone();
            let truth = altered["artifacts"]
                .as_array_mut()
                .expect("catalog")
                .iter_mut()
                .find(|artifact| artifact["id"] == "truth-jsonl")
                .expect("kept truth");
            truth[field] = json!(value);
            let error = judge(&altered).expect_err(field);
            assert!(
                format!("{error:?}").contains(
                    "truth-jsonl debe tener clase truth_jsonl y ruta outputs/truth.jsonl"
                ),
                "{field}: {error:?}"
            );
        }
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// Referencia de 300 pb con un microsatélite `(CA)50` en 101-200, una
    /// matriz generada, dos individuos diploides y VCF de texto; los loci y el
    /// tope de bases de tracto, a elección. Es el escenario de
    /// `tests/repeat_tract_budget.rs`.
    fn repeat_unit_config(root: &Path, loci: &str, cap: Option<u64>) -> PathBuf {
        let bases = format!("{}{}{}", "G".repeat(100), "CA".repeat(50), "G".repeat(100));
        fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
        let cap = cap.map_or_else(String::new, |value| {
            format!("max_tract_output_bases = {value}\n")
        });
        let config = root.join("general.toml");
        fs::write(
            &config,
            crate::test_detection_section::with_detection_section(&format!(
                r#"config_kind = "general"
schema_version = "4"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

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

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"
{cap}
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {}
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"

[vcf_output]
mode = "text"
"#,
                bases.len()
            )),
        )
        .expect("config");
        config
    }

    const SINGLE_REPEAT_ALLELE: &str = "\n[[reference.repeat_loci]]\nid = \"STR-1\"\nsequence = \"chr1\"\nstart = 101\nend = 200\nmotif = \"CA\"\n[[reference.repeat_loci.alleles]]\nunits = 45\nfrequency = 1.0\n";

    const THREE_REPEAT_ALLELES: &str = "\n[[reference.repeat_loci]]\nid = \"STR-1\"\nsequence = \"chr1\"\nstart = 101\nend = 200\nmotif = \"CA\"\n[[reference.repeat_loci.alleles]]\nunits = 45\nfrequency = 0.4\n[[reference.repeat_loci.alleles]]\nunits = 40\nfrequency = 0.35\n[[reference.repeat_loci.alleles]]\nunits = 35\nfrequency = 0.25\n";

    /// DT-038 con su pareja. Cuando el tope de bases de tracto detiene el VCF,
    /// la verdad —que cabía sola— no se conserva: se pidió con su VCF, y sola
    /// no es la unidad pedida. La pareja, con el tope justo en lo que la
    /// ejecución completa publica y un fallo provocado después, sí la conserva.
    #[test]
    fn the_tract_cap_that_stops_the_vcf_also_drops_the_truth() {
        let root = test_root("unit-cap");
        let complete = run_complete(
            &repeat_unit_config(&root, SINGLE_REPEAT_ALLELE, None),
            &root.join("complete"),
        );
        let total = complete["repeat_loci"]["tract_output"]["realized_total_tract_bases"]
            .as_u64()
            .expect("realized total");
        assert!(
            complete["repeat_loci"]["tract_output"]["realized_vcf_tract_bases"]
                .as_u64()
                .expect("VCF")
                > 0,
            "control: the VCF carries tract bases"
        );

        let exact_root = test_root("unit-cap-exact");
        let destination = exact_root.join("result");
        let (error, manifest) = run_failed(
            &repeat_unit_config(&exact_root, SINGLE_REPEAT_ALLELE, Some(total)),
            &destination,
            |_| Err(provoked_failure()),
        );
        assert_eq!(error.code(), "GPS-IO-001");
        assert!(catalog_ids(&manifest).contains(&"truth-jsonl".to_owned()));
        assert!(destination.join("outputs/truth.vcf").is_file());
        assert_external_schema(
            &destination.join("manifest.json"),
            "manifest-general-v6.schema.json",
        );

        let tight_root = test_root("unit-cap-tight");
        let destination = tight_root.join("result");
        let (error, manifest) = run_failed(
            &repeat_unit_config(&tight_root, SINGLE_REPEAT_ALLELE, Some(total - 1)),
            &destination,
            |_| Ok(()),
        );
        assert_eq!(error.code(), "GPS-CONF-008");
        assert_eq!(catalog_ids(&manifest), ["config-general"]);
        assert_eq!(manifest["vcf_outputs"], json!([]));
        assert!(!destination.join("outputs/truth.jsonl").exists());
        assert_eq!(
            manifest["repeat_loci"]["tract_output"]["realized_total_tract_bases"],
            0
        );
        assert_external_schema(
            &destination.join("manifest.json"),
            "manifest-general-v6.schema.json",
        );
        for root in [root, exact_root, tight_root] {
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    /// Bases de tracto de una verdad y de un VCF de texto publicados, contadas
    /// por un proceso ajeno al productor (Python) con el mismo criterio que
    /// `tests/repeat_tract_budget.rs`.
    fn tract_bases_counted_with_python(destination: &Path) -> (u64, u64) {
        let script = r##"
import json, sys
truth = 0
for line in open(sys.argv[1]):
    r = json.loads(line)
    for piece in r.get("current_pieces", []) + r.get("history", []):
        if piece.get("kind") == "repeat_allele":
            truth += len(piece["bases"])
    for allele in r.get("alleles", []):
        if allele.get("kind") == "repeat":
            truth += len(allele["bases"])
vcf = 0
for line in open(sys.argv[2]):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if any(k.startswith("RN=") for k in f[7].split(";")):
        vcf += len(f[3]) + sum(len(a) for a in f[4].split(","))
print(truth, vcf)
"##;
        let output = std::process::Command::new("python3")
            .arg("-c")
            .arg(script)
            .arg(destination.join("outputs/truth.jsonl"))
            .arg(destination.join("outputs/truth.vcf"))
            .output()
            .expect("independent counter");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let text = String::from_utf8(output.stdout).expect("UTF-8");
        let mut values = text
            .split_whitespace()
            .map(|value| value.parse::<u64>().expect("integer"));
        (values.next().expect("truth"), values.next().expect("VCF"))
    }

    /// DT-038 con loci de repetición y el manifiesto general v6. El paquete
    /// fallido que conserva la unidad publica en `repeat_loci` **los sorteos y
    /// las bases de tracto reales**, no ceros: la sección entera es la de la
    /// ejecución completa de la misma configuración, y las bases realizadas
    /// son las que un contador ajeno encuentra en la verdad y el VCF
    /// conservados.
    #[test]
    fn a_failed_package_with_repeat_loci_publishes_the_realized_tract_bases_it_keeps() {
        let root = test_root("unit-repeats");
        let config = repeat_unit_config(&root, THREE_REPEAT_ALLELES, None);
        let complete = run_complete(&config, &root.join("complete"));
        let destination = root.join("result");
        let (error, manifest) = run_failed(&config, &destination, |_| Err(provoked_failure()));
        assert_eq!(error.code(), "GPS-IO-001");
        assert_eq!(manifest["schema_version"], "6");
        assert_eq!(manifest["repeat_loci"], complete["repeat_loci"]);
        let accounting = &manifest["repeat_loci"]["tract_output"];
        let (truth, vcf) = tract_bases_counted_with_python(&destination);
        assert!(truth > 0 && vcf > 0, "control: both documents carry tract");
        assert_eq!(accounting["realized_truth_tract_bases"], truth);
        assert_eq!(accounting["realized_vcf_tract_bases"], vcf);
        assert_eq!(accounting["realized_total_tract_bases"], truth + vcf);
        assert_eq!(manifest["repeat_loci"]["loci"][0]["drawn_copies"], 4);
        assert_external_schema(
            &destination.join("manifest.json"),
            "manifest-general-v6.schema.json",
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// Referencia suministrada de 2000 pb con un microsatélite **insertado**
    /// (DTC-098), una matriz generada, un individuo diploide y VCF de texto,
    /// sin evento: la unidad genotípica cita la referencia final.
    fn insertion_unit_config(root: &Path) -> PathBuf {
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
            crate::test_detection_section::with_detection_section(
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

[vcf_output]
mode = "text"
"#,
            ),
        )
        .expect("config");
        config
    }

    /// DT-038 con inserción de tractos (DTC-098). Sin evento, el paquete
    /// fallido que conserva la unidad conserva también la referencia final
    /// —la del VCF y la verdad— con su mapa y la sección
    /// `reference_insertion`, junto a la suministrada de la que deriva.
    #[test]
    fn a_failed_package_with_an_inserted_reference_keeps_the_final_reference_of_its_vcf() {
        let root = test_root("unit-insertion");
        let config = insertion_unit_config(&root);
        let complete = root.join("complete");
        run_complete(&config, &complete);
        let destination = root.join("result");
        let (error, manifest) = run_failed(&config, &destination, |_| Err(provoked_failure()));
        assert_eq!(error.code(), "GPS-IO-001");
        assert_eq!(manifest["events"], serde_json::Value::Null);
        assert_eq!(
            manifest["references"],
            json!([
                {"id": "ref-main", "fasta_artifact_id": "reference-fasta"},
                {"id": "ref-main-str", "fasta_artifact_id": "reference-fasta-ref-main-str"}
            ])
        );
        assert_eq!(manifest["reference_insertion"]["id"], "ref-main-str");
        assert_eq!(manifest["vcf_outputs"][0]["reference_id"], "ref-main-str");
        assert!(catalog_ids(&manifest).contains(&"reference-map-ref-main-str".to_owned()));
        for relative in [
            "outputs/truth.jsonl",
            "outputs/truth.vcf",
            "resources/references/ref-main-str.fa",
            "resources/references/ref-main-str.map.json",
            "resources/references/ref-main.fa",
        ] {
            assert_same_bytes(&destination, &complete, relative);
        }
        assert_external_schema(
            &destination.join("manifest.json"),
            "manifest-general-v6.schema.json",
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// Ejecuta `reductions_config` alterando en disco, sin tocar el catálogo,
    /// la copia de la referencia suministrada justo antes del manifiesto: la
    /// alteración real de DT-039. El gancho da además su propio error si se
    /// pide.
    fn run_with_altered_supplied_reference(
        with_event: bool,
        hook_fails: bool,
    ) -> (Diagnostic, serde_json::Value, PathBuf, PathBuf) {
        let root = test_root("altered-reference");
        let individuals: &[(&str, &str, &str)] = if with_event {
            &[("I1", "ACGTACGT", "TGCATGCA")]
        } else {
            &[]
        };
        let config = reductions_config(&root, individuals);
        let destination = root.join("result");
        let (error, manifest) = run_failed(&config, &destination, move |workspace| {
            let artifact = workspace
                .artifacts()
                .iter()
                .find(|artifact| artifact.id == "reference-fasta")
                .expect("supplied reference copy")
                .clone();
            let path = workspace.root().join(&artifact.path);
            let mut bytes = fs::read(&path).expect("FASTA");
            let start = bytes
                .iter()
                .position(|byte| *byte == b'\n')
                .expect("header")
                + 1;
            bytes[start] = if bytes[start] == b'A' { b'C' } else { b'A' };
            fs::write(&path, bytes).expect("alteration");
            if hook_fails {
                Err(provoked_failure())
            } else {
                Ok(())
            }
        });
        (error, manifest, destination, root)
    }

    /// El paquete fallido mínimo de DT-039: solo la copia de la configuración,
    /// las reducciones `failed` y sin enlaces, ni evento ni unidad genotípica,
    /// y en los diagnósticos la causa y detrás el fallo secundario.
    fn assert_minimal_failed_package(
        manifest: &serde_json::Value,
        destination: &Path,
        diagnostics: &[&str],
    ) {
        assert_eq!(manifest["state"], "failed");
        assert_eq!(diagnostic_codes(manifest), diagnostics);
        assert_eq!(catalog_ids(manifest), ["config-general"]);
        assert_eq!(
            files_under(destination),
            ["inputs/config/general.toml", "manifest.json"]
        );
        assert_eq!(manifest["events"], serde_json::Value::Null);
        assert_eq!(manifest["technical_events"], json!([]));
        assert_eq!(manifest["references"], json!([]));
        assert_eq!(manifest["vcf_outputs"], json!([]));
        assert_eq!(manifest["population_snvs"], serde_json::Value::Null);
        for id in ["ref-small", "ref-small-b"] {
            assert_reduction_discarded(manifest, destination, id);
        }
        assert_external_schema(
            &destination.join("manifest.json"),
            "manifest-general-v4.schema.json",
        );
    }

    /// DT-039. Con la copia de la referencia suministrada alterada en disco y
    /// el gancho dando su propio error, el primer paquete fallido no pasa sus
    /// jueces (`GPS-IO-002`), con evento y, sin él, por las reducciones que
    /// conservaría. Se publica entonces el paquete mínimo, y se devuelve **la
    /// causa original**, no la del fallo secundario.
    #[test]
    fn an_altered_supplied_reference_falls_back_to_a_minimal_package_with_the_original_cause() {
        for with_event in [true, false] {
            let (error, manifest, destination, root) =
                run_with_altered_supplied_reference(with_event, true);
            assert_eq!(error.code(), "GPS-IO-001", "with_event={with_event}");
            assert!(
                format!("{error:?}").contains("fallo provocado"),
                "{error:?}"
            );
            assert!(!error.context().contains_key("secondary_failure_code"));
            assert_minimal_failed_package(&manifest, &destination, &["GPS-IO-001", "GPS-IO-002"]);
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    /// DT-039 con la alteración sola: la causa es la propia alteración, que el
    /// juez del manifiesto completo detecta (`GPS-IO-002`), y el paquete
    /// fallido mínimo se publica igual en vez de no dejar ningún destino.
    #[test]
    fn an_altered_supplied_reference_alone_still_publishes_a_minimal_failed_package() {
        let (error, manifest, destination, root) = run_with_altered_supplied_reference(true, false);
        assert_eq!(error.code(), "GPS-IO-002");
        assert_minimal_failed_package(&manifest, &destination, &["GPS-IO-002", "GPS-IO-002"]);
        fs::remove_dir_all(root).expect("cleanup");
    }

    /// DT-039. Si el paquete fallido no puede publicarse —aquí, porque el
    /// destino aparece mientras tanto—, se devuelve igualmente la causa
    /// original, con el código del fallo secundario como contexto.
    #[test]
    fn a_failed_package_that_cannot_be_published_still_returns_the_original_cause() {
        let root = test_root("unpublishable");
        let config =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/config/general-seeded.toml");
        let effective = load_and_validate_config(&config).expect("valid config");
        let destination = root.join("result");
        let squatter = destination.clone();
        let error = execute_validated_config_with_effective_document(
            &effective,
            &destination,
            ExecutionAuthorizations::default(),
            |_| Ok(()),
            std::convert::identity,
            move |_| {
                fs::create_dir(&squatter).expect("squatter destination");
                Err(provoked_failure())
            },
        )
        .expect_err("the run fails");
        assert_eq!(error.code(), "GPS-IO-001");
        assert!(
            format!("{error:?}").contains("fallo provocado"),
            "{error:?}"
        );
        assert_eq!(
            error
                .context()
                .get("secondary_failure_code")
                .map(String::as_str),
            Some("GPS-IO-003")
        );
        assert!(!destination.join("manifest.json").exists());
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
        fs::write(
            &path,
            crate::test_detection_section::with_detection_section(&config),
        )
        .expect("write config");
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
        fs::write(
            &path,
            crate::test_detection_section::with_detection_section(&config),
        )
        .expect("write config");
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
