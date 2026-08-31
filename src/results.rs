//! Documentos JSON mínimos producidos por F1.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Component, Path, PathBuf};

use serde_json::{Value, json};

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::application::{
    EffectiveConfig, EffectiveDirectEvent, EffectiveInitialGenotypeModel,
    EffectiveReducedReference, EffectiveReferenceSource, EffectiveTechnicalEvent,
    EffectiveTechnicalEventTarget, build_effective_config_with_entropy,
};
use crate::biological_contamination::contamination_target_pair_count;
use crate::configuration::ExpectedConfigKind;
use crate::diagnostics::{Diagnostic, diagnostic_matches_catalog};
use crate::event_manifest::{event_manifest_relative_path, validate_event_manifest_bytes};
use crate::execution::declared_small_event;
use crate::identity::{
    EventId, FounderKey, IndividualId, LibraryId, PopulationId, ReferenceId, SequenceId,
    derive_execution_id,
};
use crate::initial_individual_integration::IncompatibleGenotypeAuthorization;
use crate::initial_individual_materialization::{
    InitialIndividualsMaterializationError, materialize_initial_individuals,
};
use crate::integrity::{ArtifactIdentity, configuration_digest, identify_file};
use crate::population_snvs::{
    PreparedPopulationSnvs, validate_reusable_population_snvs,
    validate_reusable_population_snvs_bytes,
};
use crate::randomness::{
    EntropySource, F4TechnicalDecision, RootSeed, RuntimeSubseedScope, derive_f4_technical_subseed,
    derive_subseed, f4_uniform_closed,
};
use crate::reference::{
    AlphabetValidation, Reference, ReferenceOrigin, ReferenceSequence, ReversibleText,
    SequenceClass,
};
use crate::reference_fasta::{
    CANONICAL_FASTA_VERSION, CANONICAL_LINE_BASE_COUNT, SequenceExpectation,
    validate_supplied_fasta,
};
use crate::reference_generation::COMPOSITION_TOLERANCE;
use crate::reference_generation::{BaseFrequency, GeneratedSequenceSpec};
use crate::reference_reduction::{
    EligibleSequence, MaterializedPlan, materialize_selection, select_reduced_reference,
};
use crate::reference_reduction_io::{
    REDUCTION_MAP_DOCUMENT_KIND, REDUCTION_MAP_SCHEMA_VERSION, ReductionMapParameters,
    build_reduced_reference, reduction_map_document, validate_reduction_map_document,
};
use crate::technical_identity::derive_technical_identity_plan;
use crate::vcf_contract::{TBI_MAXIMUM_POSITION_EXCLUSIVE, VcfIndexKind, VcfOutputMode};
use crate::vcf_publication::{
    VCF_ARTIFACT_ID, VCF_ARTIFACT_KIND, VCF_BGZF_ARTIFACT_PATH, VCF_CSI_ARTIFACT_PATH,
    VCF_INDEX_ARTIFACT_ID, VCF_TBI_ARTIFACT_PATH, VCF_TEXT_ARTIFACT_PATH,
};

/// Versión vigente del manifiesto general. DTC-082 la activa junto con el
/// primer productor PE y no mantiene lectores F3.
pub const GENERAL_MANIFEST_SCHEMA_VERSION: &str = "4";

/// Versión vigente de la configuración efectiva para el productor PE F4.
pub const EFFECTIVE_CONFIG_SCHEMA_VERSION: &str = "4";

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

/// Única descripción física de un fichero publicado por una ejecución.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArtifactRecord {
    pub id: String,
    pub kind: String,
    pub path: String,
    pub size_bytes: u64,
    pub digest: String,
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
        })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConfigurationInputRecord {
    pub kind: String,
    pub original_path: String,
    pub copy: ArtifactRecord,
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

/// Enlace semántico de una referencia al catálogo físico único.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ManifestReferenceLink {
    pub id: String,
    pub fasta_artifact_id: String,
}

/// Referencia F1 resuelta desde un manifiesto general ya validado.
///
/// Esta frontera conserva juntos el ID solicitado, el enlace semántico del
/// manifiesto y la identidad física comprobada del FASTA. Los consumidores F2
/// no deben aceptar un ID ni un FASTA alternativos fuera de esta resolución.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ResolvedManifestReference {
    reference_id: ReferenceId,
    fasta_path: PathBuf,
    fasta: ArtifactRecord,
}

impl ResolvedManifestReference {
    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub fn fasta_path(&self) -> &Path {
        &self.fasta_path
    }

    #[must_use]
    pub const fn fasta(&self) -> &ArtifactRecord {
        &self.fasta
    }
}

/// Estado alcanzado por una reducción solicitada, incluso si la ejecución falla.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ManifestReductionLink {
    pub id: String,
    pub source_reference_id: String,
    pub eligible_bases: u64,
    pub region_size_bases: u64,
    pub target_mode: String,
    pub target_declared: String,
    pub target_bases: u64,
    pub root_seed: String,
    pub selection_subseed: String,
    pub selection_order: Vec<u64>,
    pub state: UnitState,
    pub fasta_artifact_id: Option<String>,
    pub map_artifact_id: Option<String>,
    pub source_fasta_artifact_id: Option<String>,
    pub map_document_kind: String,
    pub map_schema_version: String,
}

/// Enlace verificable al metadato de la matriz SNV canónica. En `reuse` el
/// recurso permanece externo e inmutable, por eso no hay identificador de
/// artefacto publicado.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PopulationSnvsManifestLink {
    pub mode: String,
    pub metadata_sha256: String,
    pub metadata_artifact_id: Option<String>,
    pub normalization_audit: Value,
    pub initial_population_statistics: Option<Value>,
}

/// Enlace de un VCF opcional a su única referencia, muestras e inventario
/// físico. Las identidades SHA-256 y tamaños viven en el catálogo compartido.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfManifestLink {
    pub state: UnitState,
    pub reference_id: String,
    pub samples: Vec<String>,
    pub mode: VcfOutputMode,
    pub vcf_artifact_id: String,
    pub index_kind: Option<VcfIndexKind>,
    pub index_artifact_id: Option<String>,
}

/// Enlace del manifiesto general al único manifiesto de un evento real. El
/// detalle técnico continúa en el documento de evento; la raíz solo enlaza su
/// identidad, estado y artefacto físico validado.
#[derive(Clone, Debug, PartialEq)]
pub struct ManifestEventLink {
    pub event_id: String,
    pub state: UnitState,
    pub manifest_artifact_id: String,
    /// Resumen agregado duplicado en la raíz para inspección sin abrir el
    /// documento de evento. No transporta procedencia por pareja.
    pub chimera_rate: f64,
    pub realized_chimera_pairs: u64,
}

pub struct GeneralManifest<'a> {
    pub execution_id: &'a str,
    pub state: UnitState,
    pub configuration_input: Option<&'a ConfigurationInputRecord>,
    pub effective_config: Option<&'a ArtifactRecord>,
    pub references: &'a [ManifestReferenceLink],
    pub reduced_references: &'a [ManifestReductionLink],
    pub population_snvs: Option<&'a PopulationSnvsManifestLink>,
    pub vcf_outputs: &'a [VcfManifestLink],
    pub events: Option<&'a [ManifestEventLink]>,
    pub authorizations: crate::execution::ExecutionAuthorizations,
    pub artifacts: &'a [ArtifactRecord],
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
    validate_reference_output(config, output)?;
    let composition = config.composition().map(|composition| {
        let values = composition.values();
        json!({"A": values[0], "C": values[1], "G": values[2], "T": values[3], "N": composition.unknown_base_fraction()})
    });
    let mut document = json!({
        "effective_config_kind": "gps_effective_config",
        "schema_version": EFFECTIVE_CONFIG_SCHEMA_VERSION,
        "source_path": config.source_path(),
        "ploidy": config.ploidy().copy_count(),
        "normalization": {
            "tolerance": config.normalization().tolerance(),
            "canonical_order": "population_id,contig,position,ref,alt,variant_id,source_id",
            "sum_algorithm": "neumaier_sum_v1",
        },
        "seed": {
            "root": config.reproducibility().root_seed().to_hex(),
            "origin": config.reproducibility().seed_origin().as_str(),
        },
        "reference": reference_json(config, output, composition),
        "reduced_references": reductions.iter().map(reduction_json).collect::<Vec<_>>(),
        "initial_populations": effective_initial_populations_json(config),
        "direct_events": config.direct_events().iter().map(direct_event_json).collect::<Vec<_>>(),
        "vcf_output": config.vcf_output().map(|mode| match mode {
            VcfOutputMode::Text => json!({"mode": "text"}),
            VcfOutputMode::Bgzf => json!({"mode": "bgzf"}),
        }),
    });
    if let Some(event) = config.technical_event() {
        document["technical_event"] = technical_event_json(config, event);
    }
    Ok(document)
}

fn effective_initial_populations_json(config: &EffectiveConfig) -> Value {
    Value::Array(
        config
            .initial_populations()
            .iter()
            .map(|population| {
                let genotype_model = match population.genotype_model() {
                    EffectiveInitialGenotypeModel::HardyWeinberg => json!({"kind": "hwe"}),
                    EffectiveInitialGenotypeModel::FDiploid {
                        inbreeding_coefficient,
                    } => json!({"kind": "f_diploid", "inbreeding_coefficient": inbreeding_coefficient}),
                };
                let mut result = json!({
                    "population_id": population.population_id().as_str(),
                    "individual_count": population.individual_count(),
                    "individual_ids": population.individual_ids().iter().map(|id| id.as_str()).collect::<Vec<_>>(),
                    "ploidy": population.ploidy().copy_count(),
                    "genotype_model": genotype_model,
                });
                if let Some(founder_keys) = population.founder_keys() {
                    result["founder_keys"] =
                        json!(founder_keys.iter().map(|key| key.as_str()).collect::<Vec<_>>());
                }
                result
            })
            .collect(),
    )
}

fn direct_event_json(event: &EffectiveDirectEvent) -> Value {
    let common = |kind: &str| {
        json!({
            "event_id": event.event_id().as_str(),
            "individual_id": event.individual_id().as_str(),
            "phase": event.phase().copy_index(),
            "kind": kind,
            "contig": event.sequence_id().as_str(),
        })
    };
    match event {
        EffectiveDirectEvent::Snp {
            position_1based,
            reference,
            alternate,
            ..
        } => {
            let mut value = common("snp");
            let object = value.as_object_mut().expect("JSON object");
            object.insert("position_1based".to_owned(), json!(position_1based));
            object.insert("ref".to_owned(), json!(char::from(*reference)));
            object.insert("alt".to_owned(), json!(char::from(*alternate)));
            value
        }
        EffectiveDirectEvent::Insertion { point, bases, .. } => {
            let mut value = common("insertion");
            let object = value.as_object_mut().expect("JSON object");
            object.insert("point".to_owned(), json!(point));
            object.insert(
                "bases".to_owned(),
                json!(std::str::from_utf8(bases).expect("effective bases are canonical ASCII")),
            );
            value
        }
        EffectiveDirectEvent::Deletion { start, end, .. } => {
            let mut value = common("deletion");
            let object = value.as_object_mut().expect("JSON object");
            object.insert("start".to_owned(), json!(start));
            object.insert("end".to_owned(), json!(end));
            value
        }
    }
}

fn technical_event_json(config: &EffectiveConfig, event: &EffectiveTechnicalEvent) -> Value {
    let target = match event.target() {
        EffectiveTechnicalEventTarget::CoverageWgs(value) => json!({"coverage_wgs": value}),
        EffectiveTechnicalEventTarget::ReadCount(value) => json!({"read_count": value}),
    };
    json!({
        "event_id": event.event_id().as_str(),
        "individual_ids": event.individual_ids().iter().map(|id| id.as_str()).collect::<Vec<_>>(),
        "profile": {
            "id": event.profile(),
            "resource_path": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_RESOURCE,
            "sha256": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256,
            "preparation": "pcr_free",
            "insert_length": {
                "minimum_bases": 250,
                "maximum_bases": 450,
            },
        },
        "read_length": event.read_length(),
        "quality": {"acgt_phred": event.acgt_phred(), "n_phred": event.n_phred()},
        "artifacts": {
            "chimera": {
                "rate": event.chimera_rate(),
            },
            "pcr": {"applicable": event.pcr_duplicates_applicable()},
            "optical_duplicates": {
                "realized_pairs": event.realized_optical_duplicate_pairs(),
            },
            "gc": {"mode": "neutral"},
            "adapters": {"enabled": false},
        },
        "biological_contamination": {
            "percent_range": {
                "minimum": event.biological_contamination().minimum_percent(),
                "maximum": event.biological_contamination().maximum_percent(),
            },
            "recipient_ids": resolved_contamination_recipients(event),
        },
        "rng_catalog": crate::randomness::F4_TECHNICAL_RNG_CATALOG,
        "id_catalog": f4_id_catalog_value(),
        "libraries": resolved_technical_libraries(config, event),
        "target": target,
    })
}

fn resolved_technical_libraries(
    config: &EffectiveConfig,
    event: &EffectiveTechnicalEvent,
) -> Vec<Value> {
    let reference_id = ReferenceId::parse(config.reference_id(), "effective_config")
        .expect("effective configuration preserves a valid reference ID");
    let identities = derive_technical_identity_plan(
        config.reproducibility().root_seed(),
        &reference_id,
        event.event_id().clone(),
        event.individual_ids(),
    )
    .expect("effective technical event preserves derivable identities");
    let recipients = resolved_contamination_recipients(event)
        .into_iter()
        .filter_map(|value| value.as_str().map(str::to_owned))
        .collect::<std::collections::BTreeSet<_>>();
    identities
        .libraries()
        .map(|identity| {
            let configured = event
                .libraries()
                .iter()
                .find(|library| library.individual_id() == identity.individual_id())
                .expect("effective event has one configured library per individual");
            json!({
                "individual_id": identity.individual_id().as_str(),
                "library_id": identity.library_id().as_str(),
                "i7": configured.i7(),
                "i5": configured.i5(),
                "is_contamination_recipient": recipients.contains(identity.individual_id().as_str()),
            })
        })
        .collect()
}

fn resolved_contamination_recipients(event: &EffectiveTechnicalEvent) -> Vec<Value> {
    event
        .biological_contamination()
        .configured_recipients()
        .map(|recipients| recipients.to_vec())
        .unwrap_or_else(|| event.individual_ids().to_vec())
        .into_iter()
        .map(|id| json!(id.as_str()))
        .collect()
}

fn f4_id_catalog_value() -> Value {
    json!({
        "version": crate::identity::ID_CATALOG_VERSION,
        "derivations": {
            "v1": ["execution", "reference", "individual", "event", "batch", "library", "read_group", "template"],
            "v2": ["PreHopInstanceId", "TechnicalPairId"],
        }
    })
}

fn f4_technical_profile_object() -> Value {
    json!({
        "id": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_V1,
        "resource_path": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_RESOURCE,
        "sha256": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256,
        "preparation": "pcr_free",
        "insert_length": {
            "minimum_bases": 250,
            "maximum_bases": 450,
        },
    })
}

pub fn effective_config_document_with_population_snvs(
    config: &EffectiveConfig,
    authorizations: crate::execution::ExecutionAuthorizations,
    output: ReferenceOutput<'_>,
    reductions: &[ReducedReferenceOutput<'_>],
    population_snvs: Option<&PopulationSnvsManifestLink>,
) -> Result<Value, Diagnostic> {
    let mut document = effective_config_document_with_reductions(config, output, reductions)?;
    let object = document.as_object_mut().ok_or_else(|| {
        Diagnostic::internal_invariant(
            "effective-config",
            "effective_config_document_is_not_object",
        )
    })?;
    object.insert(
        "authorizations".to_owned(),
        json!({
            "allow_empty_snv_matrix": authorizations.allow_empty_snv_matrix,
            "allow_incompatible_genotypes": authorizations.allow_incompatible_genotypes,
        }),
    );
    object.insert(
        "population_snvs".to_owned(),
        population_snvs.map_or(Value::Null, population_snvs_link_json),
    );
    Ok(document)
}

#[must_use]
pub fn general_manifest_document(manifest: &GeneralManifest<'_>) -> Value {
    let mut document = json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": GENERAL_MANIFEST_SCHEMA_VERSION,
        "execution_id": manifest.execution_id,
        "state": manifest.state.as_str(),
        "configuration_input": manifest.configuration_input.map(|input| json!({
            "kind": input.kind,
            "original_path": input.original_path,
            "artifact_id": input.copy.id,
        })),
        "effective_config_artifact_id": manifest.effective_config.map(|artifact| artifact.id.as_str()),
        "references": manifest.references.iter().map(manifest_reference_json).collect::<Vec<_>>(),
        "reference_reductions": manifest.reduced_references.iter().map(manifest_reduction_json).collect::<Vec<_>>(),
        "population_snvs": manifest.population_snvs.map(population_snvs_link_json),
        "vcf_outputs": manifest.vcf_outputs.iter().map(vcf_output_json).collect::<Vec<_>>(),
        "authorizations": {
            "allow_empty_snv_matrix": manifest.authorizations.allow_empty_snv_matrix,
            "allow_incompatible_genotypes": manifest.authorizations.allow_incompatible_genotypes,
        },
        "output_categories": output_categories(manifest),
        "artifacts": manifest.artifacts.iter().map(ArtifactRecord::to_json).collect::<Vec<_>>(),
        "diagnostics": manifest.diagnostics.iter().map(|diagnostic| json!({
            "code": diagnostic.code(),
            "message": diagnostic.message(),
        })).collect::<Vec<_>>(),
        "technical_events": manifest.events.map_or_else(Vec::new, |events| {
            events.iter().map(|link| technical_event_manifest_json(link, manifest)).collect()
        }),
    });
    if let Some(events) = manifest.events {
        document["events"] = Value::Array(events.iter().map(manifest_event_json).collect());
    }
    document
}

fn technical_event_manifest_json(
    link: &ManifestEventLink,
    manifest: &GeneralManifest<'_>,
) -> Value {
    let event_manifest = manifest
        .artifacts
        .iter()
        .find(|artifact| artifact.id == link.manifest_artifact_id);
    json!({
        "event_id": link.event_id,
        "effective_config": manifest.effective_config.map(ArtifactRecord::to_json),
        "profile": {
            "id": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_V1,
            "sha256": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256,
        },
        "rng_catalog": crate::randomness::F4_TECHNICAL_RNG_CATALOG,
        "id_catalog": f4_id_catalog_value(),
        "artifacts": {
            "chimera": {
                "rate": link.chimera_rate,
                "realized_pairs": link.realized_chimera_pairs,
            }
        },
        "event_manifest": event_manifest.map(ArtifactRecord::to_json),
    })
}

fn manifest_event_json(link: &ManifestEventLink) -> Value {
    json!({
        "event_id": link.event_id,
        "state": link.state.as_str(),
        "manifest_artifact_id": link.manifest_artifact_id,
    })
}

fn vcf_output_json(link: &VcfManifestLink) -> Value {
    json!({
        "state": link.state.as_str(),
        "reference_id": link.reference_id,
        "samples": link.samples,
        "mode": match link.mode { VcfOutputMode::Text => "text", VcfOutputMode::Bgzf => "bgzf" },
        "vcf_artifact_id": link.vcf_artifact_id,
        "index_kind": link.index_kind.map(|kind| match kind { VcfIndexKind::Tbi => "tbi", VcfIndexKind::Csi => "csi" }),
        "index_artifact_id": link.index_artifact_id,
    })
}

fn population_snvs_link_json(link: &PopulationSnvsManifestLink) -> Value {
    let mut object = serde_json::Map::from_iter([
        ("mode".to_owned(), json!(link.mode)),
        ("metadata_sha256".to_owned(), json!(link.metadata_sha256)),
        (
            "metadata_artifact_id".to_owned(),
            json!(link.metadata_artifact_id),
        ),
        (
            "normalization_audit".to_owned(),
            link.normalization_audit.clone(),
        ),
    ]);
    if let Some(statistics) = &link.initial_population_statistics {
        object.insert(
            "initial_population_statistics".to_owned(),
            statistics.clone(),
        );
    }
    Value::Object(object)
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

/// Analiza y valida el contrato propio y la resolución de todos los enlaces.
pub fn validate_general_manifest_bytes(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    parse_general_manifest_bytes(path, bytes).map(|_| ())
}

/// Relee cada manifiesto de evento enlazado por la raíz y comprueba que sus
/// hechos físicos siguen siendo exactamente los del catálogo general. Esta
/// puerta es deliberadamente independiente de la configuración efectiva para
/// que también proteja un fallo controlado que ya no la publica.
pub fn validate_published_event_manifests(
    root: &Path,
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    let records = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let document: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "la raíz no es un objeto"))?;
    let catalog = records
        .iter()
        .map(|record| (record.id.as_str(), record))
        .collect::<BTreeMap<_, _>>();
    let Some(events) = object.get("events") else {
        return Ok(());
    };
    let mut linked_event_artifacts = BTreeSet::new();
    let mut linked_fastq_artifacts = BTreeSet::new();
    for event in events
        .as_array()
        .ok_or_else(|| manifest_error(manifest_path, "events no es una lista"))?
    {
        let link = event
            .as_object()
            .ok_or_else(|| manifest_error(manifest_path, "un enlace de evento no es objeto"))?;
        let event_id = required_non_empty_text(link, "event_id", manifest_path)?;
        let state = required_non_empty_text(link, "state", manifest_path)?;
        let manifest_id = required_non_empty_text(link, "manifest_artifact_id", manifest_path)?;
        linked_event_artifacts.insert(manifest_id.to_owned());
        let artifact = catalog.get(manifest_id).ok_or_else(|| {
            manifest_error(manifest_path, "el artefacto del manifiesto de evento falta")
        })?;
        let event_path = root.join(&artifact.path);
        let observed = identify_file(&event_path, Path::new(&artifact.path))?;
        if observed.size_bytes() != artifact.size_bytes || observed.digest() != artifact.digest {
            return Err(manifest_error(
                manifest_path,
                "el manifiesto de evento no coincide con su identidad publicada",
            ));
        }
        let bytes = fs::read(&event_path)
            .map_err(|error| io_error(&event_path, "read_event_manifest", error))?;
        validate_event_manifest_bytes(Path::new(&artifact.path), &bytes)?;
        let event_document: Value = serde_json::from_slice(&bytes).map_err(|error| {
            manifest_error(
                manifest_path,
                format!("el manifiesto de evento no es JSON legible: {error}"),
            )
        })?;
        if event_document.get("event_id").and_then(Value::as_str) != Some(event_id)
            || event_document.get("state").and_then(Value::as_str) != Some(state)
        {
            return Err(manifest_error(
                manifest_path,
                "el enlace general no coincide con el manifiesto de evento",
            ));
        }
        let reference = event_document
            .get("reference")
            .and_then(Value::as_object)
            .and_then(|reference| reference.get("fasta"))
            .ok_or_else(|| manifest_error(manifest_path, "falta FASTA de evento"))?;
        validate_event_artifact_against_catalog(reference, &catalog, manifest_path)?;
        for artifact in event_document
            .get("outputs")
            .and_then(Value::as_array)
            .ok_or_else(|| manifest_error(manifest_path, "faltan salidas de evento"))?
        {
            if artifact.get("state").and_then(Value::as_str) == Some("complete") {
                let physical = artifact.get("artifact").ok_or_else(|| {
                    manifest_error(manifest_path, "artefacto completo de evento ausente")
                })?;
                validate_event_artifact_against_catalog(physical, &catalog, manifest_path)?;
                let id = physical.get("id").and_then(Value::as_str).ok_or_else(|| {
                    manifest_error(manifest_path, "artefacto completo de evento sin ID")
                })?;
                linked_fastq_artifacts.insert(id.to_owned());
            }
        }
    }
    let cataloged_event_artifacts = records
        .iter()
        .filter(|artifact| artifact.kind == "event_manifest")
        .map(|artifact| artifact.id.clone())
        .collect::<BTreeSet<_>>();
    let cataloged_fastq_artifacts = records
        .iter()
        .filter(|artifact| artifact.kind == "fastq")
        .map(|artifact| artifact.id.clone())
        .collect::<BTreeSet<_>>();
    if cataloged_event_artifacts != linked_event_artifacts
        || cataloged_fastq_artifacts != linked_fastq_artifacts
    {
        return Err(manifest_error(
            manifest_path,
            "el catálogo contiene o omite artefactos de evento respecto de sus manifiestos enlazados",
        ));
    }
    Ok(())
}

fn validate_event_artifact_against_catalog(
    value: &Value,
    catalog: &BTreeMap<&str, &ArtifactRecord>,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let object = value
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "artefacto de evento inválido"))?;
    let id = object
        .get("id")
        .and_then(Value::as_str)
        .ok_or_else(|| manifest_error(manifest_path, "artefacto de evento sin ID"))?;
    let record = catalog
        .get(id)
        .ok_or_else(|| manifest_error(manifest_path, "artefacto de evento no catalogado"))?;
    if value != &record.to_json() {
        return Err(manifest_error(
            manifest_path,
            "artefacto de evento no coincide con el catálogo físico general",
        ));
    }
    Ok(())
}

/// Resuelve una referencia exclusivamente desde su manifiesto general F1 y
/// comprueba que el FASTA enlazado sigue teniendo los bytes declarados.
///
/// Es la frontera reutilizable de DTC-060: el llamador aporta conjuntamente el
/// manifiesto y el `reference_id`; ni una ruta FASTA libre ni un ID aislado
/// pueden llegar al consumidor posterior.
pub fn resolve_manifest_reference(
    manifest_path: &Path,
    reference_id: &ReferenceId,
) -> Result<ResolvedManifestReference, Diagnostic> {
    let bytes = fs::read(manifest_path).map_err(|error| {
        Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
            .with("path", manifest_path.display().to_string())
            .with("operation", "read_reference_manifest")
            .with("detail", error.to_string())
    })?;
    let artifacts = parse_general_manifest_bytes(manifest_path, &bytes)?;
    let parent = manifest_path.parent().ok_or_else(|| {
        manifest_error(
            manifest_path,
            "el manifiesto no tiene un directorio padre para resolver recursos",
        )
    })?;
    let document: Value = serde_json::from_slice(&bytes).map_err(|error| {
        manifest_error(
            manifest_path,
            format!("JSON no legible tras validación: {error}"),
        )
    })?;
    let references = document
        .get("references")
        .and_then(Value::as_array)
        .ok_or_else(|| manifest_error(manifest_path, "falta references tras validación"))?;
    let fasta_artifact_id = references
        .iter()
        .filter_map(Value::as_object)
        .find(|reference| {
            reference.get("id").and_then(Value::as_str) == Some(reference_id.as_str())
        })
        .and_then(|reference| reference.get("fasta_artifact_id").and_then(Value::as_str))
        .ok_or_else(|| {
            manifest_error(
                manifest_path,
                format!(
                    "la referencia solicitada no existe: {}",
                    reference_id.as_str()
                ),
            )
        })?;
    let fasta = artifacts
        .iter()
        .find(|artifact| artifact.id == fasta_artifact_id && artifact.kind == "reference_fasta")
        .cloned()
        .ok_or_else(|| {
            manifest_error(
                manifest_path,
                format!("el FASTA de referencia no existe: {fasta_artifact_id}"),
            )
        })?;
    let fasta_path = parent.join(&fasta.path);
    let observed = identify_file(&fasta_path, &fasta.path)?;
    if observed.size_bytes() != fasta.size_bytes || observed.digest() != fasta.digest {
        return Err(Diagnostic::io(
            "GPS-IO-002",
            "El tamaño o SHA-256 de un artefacto no coincide con el valor esperado.",
        )
        .with("path", fasta_path.display().to_string())
        .with("expected_size", fasta.size_bytes.to_string())
        .with("observed_size", observed.size_bytes().to_string())
        .with("expected_digest", &fasta.digest)
        .with("observed_digest", observed.digest()));
    }

    // El digest del FASTA seleccionado no autentica por sí solo una referencia
    // reducida: su mapa puede haber sido sustituido conservando intactos el
    // FASTA y sus metadatos. Reutilizamos la puerta física F1 para reconstruir
    // todas las relaciones de la publicación antes de entregar cualquier
    // enlace al consumidor F2.
    let effective_id = document
        .get("effective_config_artifact_id")
        .and_then(Value::as_str)
        .ok_or_else(|| manifest_error(manifest_path, "falta effective_config_artifact_id"))?;
    let effective = artifacts
        .iter()
        .find(|artifact| artifact.id == effective_id && artifact.kind == "effective_config")
        .ok_or_else(|| {
            manifest_error(
                manifest_path,
                format!("la configuración efectiva no existe: {effective_id}"),
            )
        })?;
    let effective_path = parent.join(&effective.path);
    let effective_bytes = fs::read(&effective_path)
        .map_err(|error| io_error(&effective_path, "read_effective_config", error))?;
    validate_published_resources(
        parent,
        &effective_path,
        &effective_bytes,
        manifest_path,
        &bytes,
    )?;
    Ok(ResolvedManifestReference {
        reference_id: reference_id.clone(),
        fasta_path,
        fasta,
    })
}

/// Valida la configuración efectiva releyéndola y cruza todos sus enlaces con
/// el catálogo y los enlaces semánticos del manifiesto ya validado.
pub fn validate_effective_config_against_manifest_bytes(
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    match effective_schema_version(effective_path, effective_bytes)? {
        EFFECTIVE_CONFIG_SCHEMA_VERSION => validate_effective_config_v3_against_manifest_bytes(
            effective_path,
            effective_bytes,
            manifest_path,
            manifest_bytes,
        ),
        _ => unreachable!("effective_schema_version only returns supported versions"),
    }
}

fn validate_effective_config_v3_against_manifest_bytes(
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    let artifacts = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let manifest = manifest
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "la raíz no es un objeto"))?;
    if manifest.get("schema_version").and_then(Value::as_str)
        != Some(GENERAL_MANIFEST_SCHEMA_VERSION)
    {
        return Err(effective_error(
            effective_path,
            "el manifiesto no usa la versión general vigente",
        ));
    }
    let catalog = artifacts
        .iter()
        .map(|record| (record.id.as_str(), record.kind.as_str()))
        .collect::<BTreeMap<_, _>>();
    let effective_id = optional_text(manifest, "effective_config_artifact_id", manifest_path)?
        .ok_or_else(|| {
            effective_error(
                effective_path,
                "el manifiesto no enlaza la configuración efectiva",
            )
        })?;
    require_artifact_kind(&catalog, effective_id, "effective_config", manifest_path)?;
    let manifest_references = validate_references(
        required_array(manifest, "references", manifest_path)?,
        &catalog,
        manifest_path,
    )?;

    let document: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "la raíz no es un objeto"))?;
    const LEGACY_FIELDS: [&str; 11] = [
        "effective_config_kind",
        "schema_version",
        "source_path",
        "ploidy",
        "normalization",
        "seed",
        "reference",
        "population_snvs",
        "authorizations",
        "reduced_references",
        "initial_populations",
    ];
    let mut fields = LEGACY_FIELDS.to_vec();
    for field in ["direct_events", "vcf_output", "technical_event"] {
        if object.contains_key(field) {
            fields.push(field);
        }
    }
    validate_effective_closed_object(object, &fields, "effective-config.json", effective_path)?;
    if object.get("effective_config_kind").and_then(Value::as_str) != Some("gps_effective_config")
        || object.get("schema_version").and_then(Value::as_str)
            != Some(EFFECTIVE_CONFIG_SCHEMA_VERSION)
    {
        return Err(effective_error(
            effective_path,
            "clase o versión incompatible",
        ));
    }
    let source_path = effective_non_empty_text(object, "source_path", effective_path)?;
    let effective_ploidy = effective_positive_u64(object, "ploidy", effective_path)?;
    if effective_ploidy > 10 {
        return Err(effective_error(effective_path, "ploidy fuera de rango"));
    }
    validate_effective_normalization(
        effective_object(object, "normalization", effective_path)?,
        effective_path,
    )?;
    let configuration_input = required_object(manifest, "configuration_input", manifest_path)?;
    if required_non_empty_text(configuration_input, "original_path", manifest_path)? != source_path
    {
        return Err(effective_error(
            effective_path,
            "source_path no coincide con la procedencia de la copia TOML",
        ));
    }
    let seed = effective_object(object, "seed", effective_path)?;
    validate_effective_seed(seed, effective_path)?;
    let effective_root_seed = effective_non_empty_text(seed, "root", effective_path)?;
    let effective_reference = effective_object(object, "reference", effective_path)?;
    let reference = validate_effective_reference(effective_reference, &catalog, effective_path)?;
    let eligible_bases =
        effective_positive_u64(effective_reference, "total_genome_size", effective_path)?;
    if manifest_references.get(&reference.0).map(String::as_str) != Some(reference.1.as_str()) {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no coincide con el enlace del manifiesto",
        ));
    }
    let mut effective_references = BTreeMap::from([(reference.0.clone(), reference.1.clone())]);
    validate_effective_population_snvs(object, manifest, &catalog, effective_path, manifest_path)?;
    validate_effective_authorizations(object, manifest, effective_path, manifest_path)?;
    let initial_individual_ids = validate_effective_initial_populations(object, effective_path)?;
    validate_effective_vcf_against_manifest(
        object,
        manifest,
        &catalog,
        effective_reference,
        effective_path,
        manifest_path,
    )?;
    let direct_event_ids = if object.contains_key("direct_events") {
        validate_effective_direct_events(object, effective_path)?
    } else {
        BTreeSet::new()
    };
    if object.contains_key("technical_event") {
        validate_effective_technical_event(
            object,
            &initial_individual_ids,
            &direct_event_ids,
            effective_path,
        )?;
    }
    validate_effective_technical_event_link(object, manifest, effective_path, manifest_path)?;

    let manifest_reductions = required_array(manifest, "reference_reductions", manifest_path)?
        .iter()
        .map(|value| {
            let reduction = value
                .as_object()
                .ok_or_else(|| manifest_error(manifest_path, "una reducción no es objeto"))?;
            let id = required_non_empty_text(reduction, "id", manifest_path)?;
            Ok((id.to_owned(), reduction))
        })
        .collect::<Result<BTreeMap<_, _>, Diagnostic>>()?;
    let effective_reductions = effective_array(object, "reduced_references", effective_path)?;
    if effective_reductions.len() != manifest_reductions.len() {
        return Err(effective_error(
            effective_path,
            "las reducciones efectivas no coinciden con el manifiesto",
        ));
    }
    let mut effective_reduction_ids = BTreeSet::new();
    for value in effective_reductions {
        let reduction = validate_effective_reduction(value, &catalog, effective_path)?;
        if !effective_reduction_ids.insert(reduction.id.clone()) {
            return Err(effective_error(
                effective_path,
                "ID de reducción efectiva repetido",
            ));
        }
        let manifest_reduction = manifest_reductions.get(&reduction.id).ok_or_else(|| {
            effective_error(
                effective_path,
                format!("la reducción {} no existe en el manifiesto", reduction.id),
            )
        })?;
        let effective_reduction = value
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        if manifest_reduction.get("state").and_then(Value::as_str) != Some("complete")
            || manifest_reduction
                .get("eligible_bases")
                .and_then(Value::as_u64)
                != Some(eligible_bases)
            || manifest_reduction.get("root_seed").and_then(Value::as_str)
                != Some(effective_root_seed)
            || manifest_reduction
                .get("map_document_kind")
                .and_then(Value::as_str)
                != Some(REDUCTION_MAP_DOCUMENT_KIND)
            || manifest_reduction
                .get("map_schema_version")
                .and_then(Value::as_str)
                != Some(REDUCTION_MAP_SCHEMA_VERSION)
            || manifest_reduction
                .get("source_reference_id")
                .and_then(Value::as_str)
                != Some(reduction.source_reference_id.as_str())
            || manifest_reduction
                .get("fasta_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.fasta_artifact_id.as_str())
            || manifest_reduction
                .get("map_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.map_artifact_id.as_str())
            || manifest_reduction
                .get("source_fasta_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.source_fasta_artifact_id.as_str())
            || [
                "region_size_bases",
                "target",
                "selection_subseed",
                "selection_order",
            ]
            .iter()
            .any(|field| manifest_reduction.get(*field) != effective_reduction.get(*field))
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "la reducción efectiva {} no coincide con el manifiesto",
                    reduction.id
                ),
            ));
        }
        if effective_references
            .insert(reduction.id.clone(), reduction.fasta_artifact_id.clone())
            .is_some()
        {
            return Err(effective_error(
                effective_path,
                "una reducción reutiliza un ID de referencia",
            ));
        }
    }
    if effective_references != manifest_references {
        return Err(effective_error(
            effective_path,
            "el conjunto de referencias del manifiesto no coincide exactamente con la configuración efectiva",
        ));
    }
    Ok(())
}

fn validate_effective_technical_event_link(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    if !effective.contains_key("technical_event") {
        return match manifest.get("events") {
            None => Ok(()),
            Some(Value::Array(events)) if events.is_empty() => Ok(()),
            _ => Err(effective_error(
                effective_path,
                "un manifiesto con eventos exige technical_event en la configuración efectiva",
            )),
        };
    }
    let event = effective_object(effective, "technical_event", effective_path)?;
    let event_id = effective_non_empty_text(event, "event_id", effective_path)?;
    let events = manifest
        .get("events")
        .and_then(Value::as_array)
        .ok_or_else(|| {
            effective_error(
                effective_path,
                "technical_event no tiene enlace en el manifiesto general",
            )
        })?;
    if events.len() != 1 {
        return Err(effective_error(
            effective_path,
            "F3 exige exactamente un enlace de evento técnico",
        ));
    }
    let link = events[0]
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "enlace de evento inválido"))?;
    if link.get("event_id").and_then(Value::as_str) != Some(event_id)
        || link.get("state").and_then(Value::as_str) != Some("complete")
    {
        return Err(effective_error(
            effective_path,
            "el enlace del evento no coincide con technical_event efectivo",
        ));
    }
    let manifest_id = link
        .get("manifest_artifact_id")
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(effective_path, "enlace de evento sin manifiesto"))?;
    let artifacts = required_array(manifest, "artifacts", manifest_path)?;
    if !artifacts.iter().any(|artifact| {
        artifact.get("id").and_then(Value::as_str) == Some(manifest_id)
            && artifact.get("kind").and_then(Value::as_str) == Some("event_manifest")
    }) {
        return Err(effective_error(
            effective_path,
            "el manifiesto de evento enlazado no está catalogado",
        ));
    }
    Ok(())
}

/// Relee los artefactos físicos enlazados y reconstruye sus relaciones sin
/// aceptar effective/manifest como prueba de sus propias declaraciones.
pub fn validate_published_resources(
    root: &Path,
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    match effective_schema_version(effective_path, effective_bytes)? {
        EFFECTIVE_CONFIG_SCHEMA_VERSION => validate_published_resources_v3(
            root,
            effective_path,
            effective_bytes,
            manifest_path,
            manifest_bytes,
        ),
        _ => unreachable!("effective_schema_version only returns supported versions"),
    }
}

fn validate_published_resources_v3(
    root: &Path,
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    validate_published_event_manifests(root, manifest_path, manifest_bytes)?;
    validate_effective_config_against_manifest_bytes(
        effective_path,
        effective_bytes,
        manifest_path,
        manifest_bytes,
    )?;
    let artifacts = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let artifact_paths = artifacts
        .iter()
        .map(|artifact| (artifact.id.as_str(), artifact.path.as_str()))
        .collect::<BTreeMap<_, _>>();
    let effective: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let effective = effective
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "la raíz no es objeto"))?;
    let physical_config = validate_effective_targets_against_configuration(
        root,
        effective,
        &artifact_paths,
        manifest_path,
        manifest_bytes,
        effective_path,
    )?;
    let root_seed = RootSeed::parse(
        "effective-config.seed.root",
        effective["seed"]["root"]
            .as_str()
            .ok_or_else(|| effective_error(effective_path, "seed.root inválida"))?,
    )
    .map_err(|error| physical_error(effective_path, "root_seed", &error))?;
    let reference_value = effective_object(effective, "reference", effective_path)?;
    let source_reference =
        read_physical_reference(root, reference_value, &artifact_paths, effective_path)?;
    validate_physical_technical_event(
        root,
        effective,
        manifest_path,
        manifest_bytes,
        &artifacts,
        &root_seed,
        &source_reference,
        effective_path,
    )?;
    validate_physical_vcf_output(
        root,
        &physical_config,
        manifest_bytes,
        &artifacts,
        &source_reference,
        effective_path,
        manifest_path,
    )?;
    let allow_empty_snv_matrix = effective
        .get("authorizations")
        .and_then(|value| value.get("allow_empty_snv_matrix"))
        .and_then(Value::as_bool)
        .ok_or_else(|| effective_error(effective_path, "authorizations inválido"))?;
    let matrix = validate_physical_population_snvs(
        root,
        effective,
        &physical_config,
        &artifacts,
        &source_reference,
        effective_path,
        allow_empty_snv_matrix,
    )?;
    validate_physical_initial_population_statistics(
        effective,
        &physical_config,
        &source_reference,
        matrix.as_ref(),
        effective_path,
    )?;
    if reference_value.get("source").and_then(Value::as_str) == Some("simulated") {
        validate_generated_reference(
            reference_value,
            root_seed,
            &source_reference,
            effective_path,
        )?;
    }
    let mut references =
        BTreeMap::from([(source_reference.id().as_str().to_owned(), source_reference)]);

    for reduction_value in effective_array(effective, "reduced_references", effective_path)? {
        let reduction = reduction_value
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        let id = effective_non_empty_text(reduction, "id", effective_path)?;
        let source_id = effective_non_empty_text(reduction, "source_reference_id", effective_path)?;
        let source = references.get(source_id).ok_or_else(|| {
            effective_error(
                effective_path,
                format!("la referencia fuente física {source_id} no existe"),
            )
        })?;
        let derived_subseed = derive_subseed(
            root_seed,
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: source_id,
                derived_reference_id: id,
            },
        )
        .map_err(|error| physical_error(effective_path, "derive_reduction_subseed", &error))?;
        if effective_non_empty_text(reduction, "selection_subseed", effective_path)?
            != derived_subseed.to_hex()
        {
            return Err(effective_error(
                effective_path,
                format!("selection_subseed de {id} no reproduce el ámbito canónico"),
            ));
        }
        let eligible = source
            .sequences()
            .iter()
            .map(|sequence| {
                EligibleSequence::new(sequence.id().clone(), sequence.class(), sequence.length())
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| {
                effective_error(
                    effective_path,
                    format!("catálogo elegible inválido: {error:?}"),
                )
            })?;
        let target_bases = effective_positive_u64(
            effective_object(reduction, "target", effective_path)?,
            "bases",
            effective_path,
        )?;
        let region_size_bases =
            effective_positive_u64(reduction, "region_size_bases", effective_path)?;
        let mut rng = derived_subseed.local_rng();
        let selection = select_reduced_reference(
            source.id(),
            &eligible,
            target_bases,
            region_size_bases,
            &mut rng,
        )
        .map_err(|error| {
            effective_error(
                effective_path,
                format!("selección física no reproducible: {error:?}"),
            )
        })?;
        let derived_id = ReferenceId::parse(id, "published-resource-validation")
            .map_err(|error| physical_error(effective_path, "derived_reference_id", &error))?;
        let plan = materialize_selection(&selection, &derived_id, &eligible).map_err(|error| {
            effective_error(
                effective_path,
                format!("plan físico no materializable: {error:?}"),
            )
        })?;
        let declared_order = effective_array(reduction, "selection_order", effective_path)?
            .iter()
            .map(Value::as_u64)
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| effective_error(effective_path, "selection_order inválido"))?;
        let expected_order = plan
            .segments()
            .iter()
            .map(|segment| segment.selection_order())
            .collect::<Vec<_>>();
        if declared_order != expected_order {
            return Err(effective_error(
                effective_path,
                format!("selection_order de {id} no coincide con la subsemilla"),
            ));
        }
        let expected_derived = build_reduced_reference(source, &plan)
            .map_err(|error| physical_error(effective_path, "build_reduced_reference", &error))?;
        let fasta_id = effective_non_empty_text(reduction, "fasta_artifact_id", effective_path)?;
        let derived = read_reference_with_expected(
            root,
            id,
            fasta_id,
            effective_array(reduction, "sequences", effective_path)?,
            expected_derived.sequences(),
            &artifact_paths,
            effective_path,
        )?;
        compare_reference_bases(&expected_derived, &derived, effective_path, id)?;
        let map_id = effective_non_empty_text(reduction, "map_artifact_id", effective_path)?;
        let map_relative = artifact_paths.get(map_id).ok_or_else(|| {
            effective_error(effective_path, format!("mapa físico ausente: {map_id}"))
        })?;
        let map_path = root.join(map_relative);
        let map_bytes = std::fs::read(&map_path)
            .map_err(|error| io_error(&map_path, "read_reduction_map", error))?;
        let map: Value = serde_json::from_slice(&map_bytes).map_err(|error| {
            effective_error(
                effective_path,
                format!("mapa {map_id} no es JSON legible: {error}"),
            )
        })?;
        if map != reduction_map_document(&plan) {
            return Err(effective_error(
                effective_path,
                format!("mapa físico {map_id} no coincide con subsemilla y selección"),
            ));
        }
        validate_reduction_map_document(
            source,
            &derived,
            &map,
            ReductionMapParameters {
                region_size_bases,
                target_bases,
            },
        )
        .map_err(|error| physical_error(effective_path, "validate_reduction_map", &error))?;
        references.insert(id.to_owned(), derived);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn validate_physical_technical_event(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    manifest_path: &Path,
    manifest_bytes: &[u8],
    artifacts: &[ArtifactRecord],
    root_seed: &RootSeed,
    reference: &Reference,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let Some(technical_event) = effective.get("technical_event") else {
        return Ok(());
    };
    let technical_event = technical_event
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "technical_event inválido"))?;
    let event_id = EventId::parse(
        effective_non_empty_text(technical_event, "event_id", effective_path)?,
        "effective_config",
    )
    .map_err(|_| effective_error(effective_path, "event_id técnico inválido"))?;
    let individuals = effective_array(technical_event, "individual_ids", effective_path)?
        .iter()
        .map(|value| {
            let value = value
                .as_str()
                .ok_or_else(|| effective_error(effective_path, "individual_id técnico inválido"))?;
            IndividualId::parse(value, "effective_config")
                .map_err(|_| effective_error(effective_path, "individual_id técnico inválido"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let identities =
        derive_technical_identity_plan(*root_seed, reference.id(), event_id.clone(), &individuals)
            .map_err(|error| {
                physical_error(effective_path, "technical_event_identities", &error)
            })?;
    let general: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let event_link = general
        .get("events")
        .and_then(Value::as_array)
        .and_then(|events| events.first())
        .and_then(Value::as_object)
        .ok_or_else(|| effective_error(effective_path, "falta enlace físico de technical_event"))?;
    let manifest_id = required_non_empty_text(event_link, "manifest_artifact_id", manifest_path)?;
    let artifact = artifacts
        .iter()
        .find(|artifact| artifact.id == manifest_id)
        .ok_or_else(|| effective_error(effective_path, "manifiesto técnico no catalogado"))?;
    let event_path = root.join(&artifact.path);
    let event: Value = serde_json::from_slice(
        &fs::read(&event_path)
            .map_err(|error| io_error(&event_path, "read_event_manifest", error))?,
    )
    .map_err(|error| {
        effective_error(
            effective_path,
            format!("manifiesto técnico ilegible: {error}"),
        )
    })?;
    let configured_libraries = effective_array(technical_event, "libraries", effective_path)?;
    let selected = identities
        .libraries()
        .map(|library| {
            let configured = configured_libraries
                .iter()
                .filter_map(Value::as_object)
                .find(|configured| {
                    configured.get("individual_id").and_then(Value::as_str)
                        == Some(library.individual_id().as_str())
                })
                .ok_or_else(|| {
                    effective_error(
                        effective_path,
                        "falta la UDI de una biblioteca técnica en la configuración efectiva",
                    )
                })?;
            Ok(json!({
                "individual_id": library.individual_id().as_str(),
                "library_id": library.library_id().as_str(),
                "read_group_id": library.read_group_id().as_str(),
                "i7": effective_non_empty_text(configured, "i7", effective_path)?,
                "i5": effective_non_empty_text(configured, "i5", effective_path)?,
            }))
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;
    let expected_event_subseed = identities.event_subseed().to_hex();
    let profile_matches = event
        .get("profile")
        .and_then(Value::as_object)
        .is_some_and(|profile| {
            technical_event
                .get("profile")
                .and_then(Value::as_object)
                .is_some_and(|technical_profile| {
                    profile.get("id") == technical_profile.get("id")
                        && profile.get("read_length") == technical_event.get("read_length")
                        && profile.get("quality") == technical_event.get("quality")
                })
        });
    let effective_chimera_rate = technical_event
        .get("artifacts")
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("chimera"))
        .and_then(Value::as_object)
        .and_then(|chimera| chimera.get("rate"));
    let event_chimera = event
        .get("artifacts")
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("chimera"));
    let general_chimera = general
        .get("technical_events")
        .and_then(Value::as_array)
        .and_then(|events| events.first())
        .and_then(Value::as_object)
        .and_then(|technical| technical.get("artifacts"))
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("chimera"));
    validate_f4_coverage_against_effective(&event, technical_event, reference, effective_path)?;
    validate_physical_biological_contamination_summary(
        &event,
        technical_event,
        effective,
        root_seed,
        &event_id,
        effective_path,
    )?;
    if event.get("event_id").and_then(Value::as_str) != Some(event_id.as_str())
        || event.get("event_subseed").and_then(Value::as_str)
            != Some(expected_event_subseed.as_str())
        || event.get("technical_batch_id").and_then(Value::as_str)
            != Some(identities.batch_id().as_str())
        || event
            .get("reference")
            .and_then(Value::as_object)
            .and_then(|reference| reference.get("reference_id"))
            .and_then(Value::as_str)
            != Some(reference.id().as_str())
        || event.get("selected_individuals").and_then(Value::as_array) != Some(&selected)
        || !profile_matches
        || event_chimera
            .and_then(Value::as_object)
            .and_then(|chimera| chimera.get("rate"))
            != effective_chimera_rate
        || event_chimera != general_chimera
    {
        return Err(effective_error(
            effective_path,
            "el manifiesto técnico físico no coincide con la configuración efectiva, semilla o referencia",
        ));
    }
    Ok(())
}

/// Contrasta el resumen agregado F4.4 con la configuración efectiva sin
/// retener ni reconstruir parejas. La materialización solo entrega los hechos
/// auditables del plan; aquí no se reintroduce procedencia individual.
fn validate_physical_biological_contamination_summary(
    event: &Value,
    technical_event: &serde_json::Map<String, Value>,
    effective: &serde_json::Map<String, Value>,
    root_seed: &RootSeed,
    event_id: &EventId,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let contamination =
        effective_object(technical_event, "biological_contamination", effective_path)?;
    let range = effective_object(contamination, "percent_range", effective_path)?;
    let minimum = range
        .get("minimum")
        .and_then(Value::as_f64)
        .ok_or_else(|| {
            effective_error(effective_path, "mínimo efectivo de contaminación inválido")
        })?;
    let maximum = range
        .get("maximum")
        .and_then(Value::as_f64)
        .ok_or_else(|| {
            effective_error(effective_path, "máximo efectivo de contaminación inválido")
        })?;
    let configured_recipients = effective_array(contamination, "recipient_ids", effective_path)?
        .iter()
        .map(|value| {
            value.as_str().map(str::to_owned).ok_or_else(|| {
                effective_error(
                    effective_path,
                    "receptora efectiva de contaminación inválida",
                )
            })
        })
        .collect::<Result<BTreeSet<_>, _>>()?;
    let selected_libraries = effective_array(technical_event, "libraries", effective_path)?
        .iter()
        .map(|value| {
            let library = value.as_object().ok_or_else(|| {
                effective_error(effective_path, "biblioteca técnica efectiva inválida")
            })?;
            Ok((
                effective_non_empty_text(library, "individual_id", effective_path)?.to_owned(),
                effective_non_empty_text(library, "library_id", effective_path)?.to_owned(),
            ))
        })
        .collect::<Result<BTreeMap<_, _>, Diagnostic>>()?;
    let populations = effective_array(effective, "initial_populations", effective_path)?
        .iter()
        .map(|value| {
            let population = value.as_object().ok_or_else(|| {
                effective_error(effective_path, "población inicial efectiva inválida")
            })?;
            let population_id =
                effective_non_empty_text(population, "population_id", effective_path)?;
            let individuals = effective_array(population, "individual_ids", effective_path)?;
            individuals
                .iter()
                .map(|value| {
                    value
                        .as_str()
                        .map(|individual| (individual.to_owned(), population_id.to_owned()))
                        .ok_or_else(|| {
                            effective_error(effective_path, "individual inicial efectivo inválido")
                        })
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?
        .into_iter()
        .flatten()
        .collect::<BTreeMap<_, _>>();
    let records = event
        .get("artifacts")
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("biological_contamination"))
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(effective_path, "falta resumen físico de contaminación"))?;
    // Un resumen neutro puede omitirse. Cualquier intervalo positivo exige el
    // plan entero antes de que la publicación PE (F4.5) pueda abrir una ruta.
    if records.is_empty() {
        if maximum == 0.0 {
            return Ok(());
        }
        return Err(effective_error(
            effective_path,
            "un rango positivo de contaminación exige resumen por receptora",
        ));
    }
    let mut observed_recipients = BTreeSet::new();
    for value in records {
        let record = value.as_object().ok_or_else(|| {
            effective_error(effective_path, "resumen físico de contaminación inválido")
        })?;
        let recipient = effective_non_empty_text(record, "recipient_id", effective_path)?;
        let recipient_library =
            effective_non_empty_text(record, "recipient_library_id", effective_path)?;
        if !configured_recipients.contains(recipient)
            || selected_libraries.get(recipient).map(String::as_str) != Some(recipient_library)
            || !observed_recipients.insert(recipient.to_owned())
        {
            return Err(effective_error(
                effective_path,
                "receptora o biblioteca del resumen de contaminación no coincide con la configuración efectiva",
            ));
        }
        let record_range = effective_object(record, "percent_range", effective_path)?;
        if record_range.get("minimum").and_then(Value::as_f64) != Some(minimum)
            || record_range.get("maximum").and_then(Value::as_f64) != Some(maximum)
        {
            return Err(effective_error(
                effective_path,
                "rango del resumen de contaminación no coincide con la configuración efectiva",
            ));
        }
        let sampled = record
            .get("sampled_percent")
            .and_then(Value::as_f64)
            .ok_or_else(|| effective_error(effective_path, "porcentaje contaminante inválido"))?;
        let clean = record
            .get("clean_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento limpio contaminante inválido")
            })?;
        let target = record
            .get("target_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| effective_error(effective_path, "objetivo contaminante inválido"))?;
        let removed = record
            .get("removed_own_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento retirado contaminante inválido")
            })?;
        let inserted = record
            .get("inserted_contaminant_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento insertado contaminante inválido")
            })?;
        let retained = record
            .get("retained_own_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento retenido contaminante inválido")
            })?;
        let final_pairs = record
            .get("final_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento final contaminante inválido")
            })?;
        let recipient_id = IndividualId::parse(recipient, "effective_config")
            .map_err(|_| effective_error(effective_path, "receptora contaminante inválida"))?;
        let library_id = LibraryId::parse(recipient_library, "effective_config").map_err(|_| {
            effective_error(effective_path, "biblioteca receptora contaminante inválida")
        })?;
        let percent_subseed = derive_f4_technical_subseed(
            *root_seed,
            F4TechnicalDecision::contamination_percent(event_id, &library_id, &recipient_id),
        )
        .map_err(|_| {
            effective_error(
                effective_path,
                "no se pudo rederivar una subsemilla contaminante",
            )
        })?;
        let expected_subseeds = [
            ("percent", Ok(percent_subseed)),
            (
                "donors",
                derive_f4_technical_subseed(
                    *root_seed,
                    F4TechnicalDecision::contamination_donors(event_id, &library_id, &recipient_id),
                ),
            ),
            (
                "removal",
                derive_f4_technical_subseed(
                    *root_seed,
                    F4TechnicalDecision::contamination_recipients(
                        event_id,
                        &library_id,
                        &recipient_id,
                    ),
                ),
            ),
        ];
        let subseeds = effective_object(record, "subseeds", effective_path)?;
        for (name, expected) in expected_subseeds {
            let expected = expected.map_err(|_| {
                effective_error(
                    effective_path,
                    "no se pudo rederivar una subsemilla contaminante",
                )
            })?;
            if subseeds.get(name).and_then(Value::as_str) != Some(expected.to_hex().as_str()) {
                return Err(effective_error(
                    effective_path,
                    "subsemilla de contaminación no coincide con el catálogo F4",
                ));
            }
        }
        let expected_sampled =
            f4_uniform_closed(&mut percent_subseed.local_rng(), minimum, maximum).map_err(
                |_| {
                    effective_error(
                        effective_path,
                        "no se pudo remuestrear el porcentaje contaminante",
                    )
                },
            )?;
        if sampled.to_bits() != expected_sampled.to_bits() {
            return Err(effective_error(
                effective_path,
                "porcentaje contaminante no coincide con la subsemilla F4",
            ));
        }
        if contamination_target_pair_count(clean, sampled).ok() != Some(target) {
            return Err(effective_error(
                effective_path,
                "objetivo contaminante no coincide con el redondeo normativo",
            ));
        }
        if removed != target
            || inserted != target
            || retained.checked_add(removed) != Some(clean)
            || retained.checked_add(inserted) != Some(final_pairs)
            || final_pairs != clean
        {
            return Err(effective_error(
                effective_path,
                "los contadores físicos de contaminación no conservan las parejas",
            ));
        }
        let recipient_population = populations.get(recipient).ok_or_else(|| {
            effective_error(
                effective_path,
                "receptora contaminante fuera de población seleccionada",
            )
        })?;
        let contributions = effective_array(record, "donor_contributions", effective_path)?;
        let contribution_total = contributions.iter().try_fold(0_u64, |total, value| {
            let contribution = value
                .as_object()
                .ok_or_else(|| effective_error(effective_path, "aporte contaminante inválido"))?;
            let donor = effective_non_empty_text(contribution, "individual_id", effective_path)?;
            let pairs = contribution
                .get("pairs")
                .and_then(Value::as_u64)
                .ok_or_else(|| {
                    effective_error(effective_path, "pares de aporte contaminante inválidos")
                })?;
            if donor == recipient
                || !selected_libraries.contains_key(donor)
                || populations.get(donor) != Some(recipient_population)
            {
                return Err(effective_error(
                    effective_path,
                    "donante contaminante fuera de la población seleccionada",
                ));
            }
            total.checked_add(pairs).ok_or_else(|| {
                effective_error(effective_path, "desbordamiento de aportes contaminantes")
            })
        })?;
        if contribution_total != target {
            return Err(effective_error(
                effective_path,
                "los aportes físicos de contaminación no suman el objetivo",
            ));
        }
    }
    if observed_recipients != configured_recipients {
        return Err(effective_error(
            effective_path,
            "el resumen de contaminación no cubre exactamente las receptoras efectivas",
        ));
    }
    Ok(())
}

/// Comprueba los hechos que el manifiesto F4 puede contrastar sin convertir la
/// observación FASTQ final en una nueva configuración u objetivo.
fn validate_f4_coverage_against_effective(
    event: &Value,
    technical_event: &serde_json::Map<String, Value>,
    reference: &Reference,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let target = effective_object(technical_event, "target", effective_path)?;
    let target_x = target
        .get("coverage_wgs")
        .and_then(Value::as_f64)
        .filter(|value| value.is_finite() && *value >= 0.0)
        .ok_or_else(|| effective_error(effective_path, "objetivo WGS técnico inválido"))?;
    let expected_denominator = reference.total_length();
    let coverage = event
        .get("coverage")
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(effective_path, "falta cobertura F4 de evento"))?;
    for entry in coverage {
        let entry = entry
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "cobertura F4 de evento inválida"))?;
        if entry.get("target_x").and_then(Value::as_f64) != Some(target_x)
            || entry.get("denominator_bases").and_then(Value::as_u64) != Some(expected_denominator)
        {
            return Err(effective_error(
                effective_path,
                "la cobertura F4 no coincide con el objetivo o denominador estructural efectivo",
            ));
        }
    }
    Ok(())
}

#[allow(dead_code)]
fn validate_published_resources_v1(
    root: &Path,
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    validate_effective_config_against_manifest_bytes(
        effective_path,
        effective_bytes,
        manifest_path,
        manifest_bytes,
    )?;
    let artifacts = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let artifact_paths = artifacts
        .iter()
        .map(|artifact| (artifact.id.as_str(), artifact.path.as_str()))
        .collect::<BTreeMap<_, _>>();
    let effective: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let effective = effective
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "la raíz no es objeto"))?;
    validate_effective_targets_against_configuration_v1(
        root,
        effective,
        &artifact_paths,
        manifest_path,
        manifest_bytes,
        effective_path,
    )?;
    let root_seed = RootSeed::parse(
        "effective-config.seed.root",
        effective["seed"]["root"]
            .as_str()
            .ok_or_else(|| effective_error(effective_path, "seed.root inválida"))?,
    )
    .map_err(|error| physical_error(effective_path, "root_seed", &error))?;
    let reference_value = effective_object(effective, "reference", effective_path)?;
    let source_reference =
        read_physical_reference(root, reference_value, &artifact_paths, effective_path)?;
    if reference_value.get("source").and_then(Value::as_str) == Some("simulated") {
        validate_generated_reference(
            reference_value,
            root_seed,
            &source_reference,
            effective_path,
        )?;
    }
    let mut references =
        BTreeMap::from([(source_reference.id().as_str().to_owned(), source_reference)]);

    for reduction_value in effective_array(effective, "reduced_references", effective_path)? {
        let reduction = reduction_value
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        let id = effective_non_empty_text(reduction, "id", effective_path)?;
        let source_id = effective_non_empty_text(reduction, "source_reference_id", effective_path)?;
        let source = references.get(source_id).ok_or_else(|| {
            effective_error(
                effective_path,
                format!("la referencia fuente física {source_id} no existe"),
            )
        })?;
        let derived_subseed = derive_subseed(
            root_seed,
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: source_id,
                derived_reference_id: id,
            },
        )
        .map_err(|error| physical_error(effective_path, "derive_reduction_subseed", &error))?;
        if effective_non_empty_text(reduction, "selection_subseed", effective_path)?
            != derived_subseed.to_hex()
        {
            return Err(effective_error(
                effective_path,
                format!("selection_subseed de {id} no reproduce el ámbito canónico"),
            ));
        }
        let eligible = source
            .sequences()
            .iter()
            .map(|sequence| {
                EligibleSequence::new(sequence.id().clone(), sequence.class(), sequence.length())
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| {
                effective_error(
                    effective_path,
                    format!("catálogo elegible inválido: {error:?}"),
                )
            })?;
        let target_bases = effective_positive_u64(
            effective_object(reduction, "target", effective_path)?,
            "bases",
            effective_path,
        )?;
        let region_size_bases =
            effective_positive_u64(reduction, "region_size_bases", effective_path)?;
        let mut rng = derived_subseed.local_rng();
        let selection = select_reduced_reference(
            source.id(),
            &eligible,
            target_bases,
            region_size_bases,
            &mut rng,
        )
        .map_err(|error| {
            effective_error(
                effective_path,
                format!("selección física no reproducible: {error:?}"),
            )
        })?;
        let derived_id = ReferenceId::parse(id, "published-resource-validation")
            .map_err(|error| physical_error(effective_path, "derived_reference_id", &error))?;
        let plan = materialize_selection(&selection, &derived_id, &eligible).map_err(|error| {
            effective_error(
                effective_path,
                format!("plan físico no materializable: {error:?}"),
            )
        })?;
        let declared_order = effective_array(reduction, "selection_order", effective_path)?
            .iter()
            .map(Value::as_u64)
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| effective_error(effective_path, "selection_order inválido"))?;
        let expected_order = plan
            .segments()
            .iter()
            .map(|segment| segment.selection_order())
            .collect::<Vec<_>>();
        if declared_order != expected_order {
            return Err(effective_error(
                effective_path,
                format!("selection_order de {id} no coincide con la subsemilla"),
            ));
        }
        let expected_derived = build_reduced_reference(source, &plan)
            .map_err(|error| physical_error(effective_path, "build_reduced_reference", &error))?;
        let fasta_id = effective_non_empty_text(reduction, "fasta_artifact_id", effective_path)?;
        let derived = read_reference_with_expected(
            root,
            id,
            fasta_id,
            effective_array(reduction, "sequences", effective_path)?,
            expected_derived.sequences(),
            &artifact_paths,
            effective_path,
        )?;
        compare_reference_bases(&expected_derived, &derived, effective_path, id)?;
        let map_id = effective_non_empty_text(reduction, "map_artifact_id", effective_path)?;
        let map_relative = artifact_paths.get(map_id).ok_or_else(|| {
            effective_error(effective_path, format!("mapa físico ausente: {map_id}"))
        })?;
        let map_path = root.join(map_relative);
        let map_bytes = std::fs::read(&map_path)
            .map_err(|error| io_error(&map_path, "read_reduction_map", error))?;
        let map: Value = serde_json::from_slice(&map_bytes).map_err(|error| {
            effective_error(
                effective_path,
                format!("mapa {map_id} no es JSON legible: {error}"),
            )
        })?;
        if map != reduction_map_document(&plan) {
            return Err(effective_error(
                effective_path,
                format!("mapa físico {map_id} no coincide con subsemilla y selección"),
            ));
        }
        validate_reduction_map_document(
            source,
            &derived,
            &map,
            ReductionMapParameters {
                region_size_bases,
                target_bases,
            },
        )
        .map_err(|error| physical_error(effective_path, "validate_reduction_map", &error))?;
        references.insert(id.to_owned(), derived);
    }
    Ok(())
}

#[allow(dead_code)]
fn validate_effective_targets_against_configuration_v1(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    artifact_paths: &BTreeMap<&str, &str>,
    manifest_path: &Path,
    manifest_bytes: &[u8],
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let configuration_id = manifest["configuration_input"]["artifact_id"]
        .as_str()
        .ok_or_else(|| effective_error(effective_path, "configuration_input no enlaza TOML"))?;
    let relative = artifact_paths.get(configuration_id).ok_or_else(|| {
        effective_error(
            effective_path,
            format!("TOML físico ausente: {configuration_id}"),
        )
    })?;
    let configuration_path = root.join(relative);
    let bytes = std::fs::read(&configuration_path)
        .map_err(|error| io_error(&configuration_path, "read_configuration", error))?;
    let digest = configuration_digest(&bytes)
        .map_err(|error| physical_error(effective_path, "configuration_digest", &error))?;
    let original_configuration_path = manifest["configuration_input"]["original_path"]
        .as_str()
        .ok_or_else(|| effective_error(effective_path, "configuration_input sin procedencia"))?;
    let admitted = admit_toml_header(
        Path::new(original_configuration_path),
        bytes,
        ExpectedConfigKind::General,
    )
    .map_err(|error| physical_error(effective_path, "admit_configuration", &error))?;
    let declared = parse_declared_config(admitted)
        .map_err(|error| physical_error(effective_path, "parse_configuration", &error))?;
    let seed = effective_object(effective, "seed", effective_path)?;
    let root_seed = RootSeed::parse(
        "effective-config.seed.root",
        effective_non_empty_text(seed, "root", effective_path)?,
    )
    .map_err(|error| physical_error(effective_path, "root_seed", &error))?;
    let mut entropy = ReplayEntropy(root_seed);
    let expected = build_effective_config_with_entropy(&declared, &mut entropy)
        .map_err(|error| physical_error(effective_path, "rebuild_effective_config", &error))?;
    let expected_execution = derive_execution_id(root_seed, &digest)
        .map_err(|error| physical_error(effective_path, "derive_execution_id", &error))?;
    if manifest.get("execution_id").and_then(Value::as_str) != Some(expected_execution.as_str()) {
        return Err(effective_error(
            effective_path,
            "execution_id no coincide con la semilla y los bytes TOML",
        ));
    }

    if seed.get("root").and_then(Value::as_str)
        != Some(expected.reproducibility().root_seed().to_hex().as_str())
        || seed.get("origin").and_then(Value::as_str)
            != Some(expected.reproducibility().seed_origin().as_str())
    {
        return Err(effective_error(
            effective_path,
            "la semilla efectiva no coincide con la declaración TOML",
        ));
    }

    let reference = effective_object(effective, "reference", effective_path)?;
    let expected_source = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => "simulated",
        EffectiveReferenceSource::Supplied { .. } => "supplied",
    };
    if effective_non_empty_text(reference, "id", effective_path)? != expected.reference_id()
        || reference.get("source").and_then(Value::as_str) != Some(expected_source)
        || reference.get("total_genome_size").and_then(Value::as_u64)
            != Some(expected.total_genome_size())
    {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no coincide con la declaración TOML",
        ));
    }
    let expected_original_path = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => Value::Null,
        EffectiveReferenceSource::Supplied { original_path } => {
            Value::String(original_path.to_string_lossy().into_owned())
        }
    };
    if reference.get("original_fasta_path") != Some(&expected_original_path) {
        return Err(effective_error(
            effective_path,
            "la procedencia FASTA efectiva no coincide con la declaración TOML",
        ));
    }
    let expected_canonical_format = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => canonical_fasta_json(),
        EffectiveReferenceSource::Supplied { .. } => Value::Null,
    };
    if reference.get("canonical_format") != Some(&expected_canonical_format) {
        return Err(effective_error(
            effective_path,
            "el formato FASTA efectivo no corresponde al origen declarado",
        ));
    }

    let expected_composition = expected.composition().map(|composition| {
        let values = composition.values();
        json!({"A": values[0], "C": values[1], "G": values[2], "T": values[3]})
    });
    if reference.get("composition") != Some(&expected_composition.unwrap_or(Value::Null)) {
        return Err(effective_error(
            effective_path,
            "la composición efectiva no coincide con la declaración TOML",
        ));
    }

    let effective_sequences = effective_array(reference, "sequences", effective_path)?;
    if effective_sequences.len() != expected.sequences().len() {
        return Err(effective_error(
            effective_path,
            "el número de secuencias efectivo no coincide con el TOML",
        ));
    }
    for (observed, expected) in effective_sequences.iter().zip(expected.sequences()) {
        let observed = observed
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        if effective_non_empty_text(observed, "id", effective_path)? != expected.id()
            || observed.get("class").and_then(Value::as_str)
                != Some(sequence_class(expected.class()))
            || observed.get("length").and_then(Value::as_u64) != Some(expected.length())
        {
            return Err(effective_error(
                effective_path,
                "orden, ID, clase o longitud de secuencia no coincide con el TOML",
            ));
        }
    }

    let effective_reductions = effective_array(effective, "reduced_references", effective_path)?;
    if effective_reductions.len() != expected.reduced_references().len() {
        return Err(effective_error(
            effective_path,
            "las reducciones efectivas no coinciden con las declaradas en TOML",
        ));
    }
    for (expected, effective) in expected
        .reduced_references()
        .iter()
        .zip(effective_reductions)
    {
        let effective = effective
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        let target = effective_object(effective, "target", effective_path)?;
        let expected_declared = expected.target().declared_value();
        if effective_non_empty_text(effective, "id", effective_path)? != expected.id()
            || effective_non_empty_text(effective, "source_reference_id", effective_path)?
                != expected.source_reference_id()
            || effective.get("region_size_bases").and_then(Value::as_u64)
                != Some(expected.region_size_bases())
            || target.get("mode").and_then(Value::as_str) != Some(expected.target().mode())
            || target.get("declared").and_then(Value::as_str) != Some(expected_declared.as_str())
            || target.get("bases").and_then(Value::as_u64) != Some(expected.target_bases())
            || effective.get("selection_subseed").and_then(Value::as_str)
                != Some(expected.selection_subseed().to_hex().as_str())
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "la reducción efectiva {} no coincide completamente con el TOML",
                    expected.id()
                ),
            ));
        }
    }
    Ok(())
}

fn validate_physical_population_snvs(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    physical_config: &EffectiveConfig,
    artifacts: &[ArtifactRecord],
    reference: &Reference,
    effective_path: &Path,
    allow_empty_snv_matrix: bool,
) -> Result<Option<PreparedPopulationSnvs>, Diagnostic> {
    let Some(link) = effective.get("population_snvs").and_then(Value::as_object) else {
        return Ok(None);
    };
    if link.get("mode").and_then(Value::as_str) == Some("reuse") {
        let crate::application::EffectivePopulationSnvs::Reuse { metadata_path } =
            physical_config.population_snvs().ok_or_else(|| {
                effective_error(effective_path, "reutilización SNV sin fuente física")
            })?
        else {
            return Err(effective_error(
                effective_path,
                "modo SNV efectivo no coincide con reutilización",
            ));
        };
        let reference_id = effective_object(effective, "reference", effective_path)?
            .get("fasta_artifact_id")
            .and_then(Value::as_str)
            .ok_or_else(|| effective_error(effective_path, "referencia sin artefacto FASTA"))?;
        let reference_artifact = artifacts
            .iter()
            .find(|artifact| artifact.id == reference_id)
            .ok_or_else(|| effective_error(effective_path, "FASTA SNV ausente del catálogo"))?;
        let (prepared, _, metadata_sha256, normalization_audit) =
            validate_reusable_population_snvs(
                metadata_path,
                reference,
                &reference_artifact.digest,
                physical_config.normalization().tolerance(),
                allow_empty_snv_matrix,
            )
            .map_err(|error| {
                effective_error(
                    effective_path,
                    format!("reuso SNV físico inválido: {error}"),
                )
            })?;
        if link.get("metadata_sha256").and_then(Value::as_str) != Some(metadata_sha256.as_str())
            || link.get("normalization_audit") != Some(&normalization_audit)
        {
            return Err(effective_error(
                effective_path,
                "enlace SNV reutilizado no coincide con el recurso físico",
            ));
        }
        return Ok(Some(prepared));
    }
    let metadata_id = link
        .get("metadata_artifact_id")
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(effective_path, "population_snvs publicado sin metadata"))?;
    let metadata = artifacts
        .iter()
        .find(|artifact| artifact.id == metadata_id)
        .ok_or_else(|| effective_error(effective_path, "metadata SNV ausente del catálogo"))?;
    let matrix = artifacts
        .iter()
        .find(|artifact| artifact.id == "population-snvs-tsv")
        .ok_or_else(|| effective_error(effective_path, "TSV SNV ausente del catálogo"))?;
    let reference_id = effective_object(effective, "reference", effective_path)?
        .get("fasta_artifact_id")
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(effective_path, "referencia sin artefacto FASTA"))?;
    let reference_artifact = artifacts
        .iter()
        .find(|artifact| artifact.id == reference_id)
        .ok_or_else(|| effective_error(effective_path, "FASTA SNV ausente del catálogo"))?;
    let metadata_path = root.join(&metadata.path);
    let matrix_path = root.join(&matrix.path);
    let metadata_bytes = read_verified_artifact_bytes(&metadata_path, metadata, effective_path)?;
    let matrix_bytes = read_verified_artifact_bytes(&matrix_path, matrix, effective_path)?;
    let (prepared, provenance, metadata_sha256, normalization_audit) =
        validate_reusable_population_snvs_bytes(
            &metadata_bytes,
            &matrix_bytes,
            reference,
            &reference_artifact.digest,
            effective
                .get("normalization")
                .and_then(|value| value.get("tolerance"))
                .and_then(Value::as_f64)
                .ok_or_else(|| {
                    effective_error(effective_path, "normalization.tolerance inválida")
                })?,
            allow_empty_snv_matrix,
        )
        .map_err(|error| {
            effective_error(effective_path, format!("par SNV físico inválido: {error}"))
        })?;
    if link.get("metadata_sha256").and_then(Value::as_str) != Some(metadata_sha256.as_str())
        || link.get("mode").and_then(Value::as_str) != Some(provenance.as_str())
        || link.get("normalization_audit") != Some(&normalization_audit)
    {
        return Err(effective_error(
            effective_path,
            "enlace population_snvs no coincide con el par físico",
        ));
    }
    Ok(Some(prepared))
}

fn validate_physical_initial_population_statistics(
    effective: &serde_json::Map<String, Value>,
    physical_config: &EffectiveConfig,
    reference: &Reference,
    matrix: Option<&PreparedPopulationSnvs>,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let declared = effective
        .get("population_snvs")
        .and_then(Value::as_object)
        .and_then(|link| link.get("initial_population_statistics"));
    if physical_config.initial_populations().is_empty() {
        if declared.is_some() {
            return Err(effective_error(
                effective_path,
                "hay estadísticas iniciales sin poblaciones iniciales",
            ));
        }
        return Ok(());
    }
    let matrix = matrix.ok_or_else(|| {
        effective_error(
            effective_path,
            "faltan matriz y estadísticas para poblaciones iniciales",
        )
    })?;
    let authorizations = effective
        .get("authorizations")
        .and_then(Value::as_object)
        .ok_or_else(|| effective_error(effective_path, "authorizations inválido"))?;
    let authorization = if authorizations
        .get("allow_incompatible_genotypes")
        .and_then(Value::as_bool)
        .ok_or_else(|| effective_error(effective_path, "allow_incompatible_genotypes inválido"))?
    {
        IncompatibleGenotypeAuthorization::allowed()
    } else {
        IncompatibleGenotypeAuthorization::denied()
    };
    let rebuilt = materialize_initial_individuals(
        physical_config,
        reference,
        matrix,
        authorization,
        |individual_id| {
            physical_config
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
    .map_err(|error| {
        effective_error(
            effective_path,
            format!("estadísticas iniciales no reconstruibles: {error}"),
        )
    })?;
    if declared != Some(rebuilt.statistics()) {
        return Err(effective_error(
            effective_path,
            "estadísticas iniciales no coinciden con matriz, semilla y eventos físicos",
        ));
    }
    Ok(())
}

fn read_verified_artifact_bytes(
    path: &Path,
    artifact: &ArtifactRecord,
    effective_path: &Path,
) -> Result<Vec<u8>, Diagnostic> {
    let bytes = fs::read(path).map_err(|error| io_error(path, "read_population_snvs", error))?;
    let observed = identify_file(path, Path::new(&artifact.path))?;
    if observed.size_bytes() != artifact.size_bytes || observed.digest() != artifact.digest {
        return Err(effective_error(
            effective_path,
            format!(
                "artefacto SNV físico {} no coincide con su SHA-256",
                artifact.id
            ),
        ));
    }
    Ok(bytes)
}

fn validate_physical_vcf_output(
    root: &Path,
    effective: &EffectiveConfig,
    manifest_bytes: &[u8],
    artifacts: &[ArtifactRecord],
    reference: &Reference,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let manifest = manifest
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "la raíz no es un objeto"))?;
    if manifest.get("schema_version").and_then(Value::as_str)
        != Some(GENERAL_MANIFEST_SCHEMA_VERSION)
    {
        return Err(effective_error(
            effective_path,
            "el manifiesto no usa la versión general vigente",
        ));
    }
    let outputs = required_array(manifest, "vcf_outputs", manifest_path)?;
    let expected = match effective.vcf_output() {
        None => Vec::new(),
        Some(VcfOutputMode::Text) => {
            vec![(VCF_ARTIFACT_ID, VCF_ARTIFACT_KIND, VCF_TEXT_ARTIFACT_PATH)]
        }
        Some(VcfOutputMode::Bgzf) => {
            let index = if reference
                .sequences()
                .iter()
                .all(|sequence| sequence.length() < TBI_MAXIMUM_POSITION_EXCLUSIVE)
            {
                ("vcf_tbi_index", VCF_TBI_ARTIFACT_PATH)
            } else {
                ("vcf_csi_index", VCF_CSI_ARTIFACT_PATH)
            };
            vec![
                (VCF_ARTIFACT_ID, VCF_ARTIFACT_KIND, VCF_BGZF_ARTIFACT_PATH),
                (VCF_INDEX_ARTIFACT_ID, index.0, index.1),
            ]
        }
    };
    if (effective.vcf_output().is_none() && !outputs.is_empty())
        || (effective.vcf_output().is_some() && outputs.len() != 1)
    {
        return Err(effective_error(
            effective_path,
            "vcf_outputs no coincide con la activación efectiva",
        ));
    }
    let actual = artifacts
        .iter()
        .filter(|artifact| {
            matches!(
                artifact.kind.as_str(),
                "truth_vcf" | "vcf_tbi_index" | "vcf_csi_index"
            )
        })
        .collect::<Vec<_>>();
    if actual.len() != expected.len() {
        return Err(effective_error(
            effective_path,
            "el catálogo físico VCF no contiene exactamente los artefactos esperados",
        ));
    }
    for (id, kind, relative_path) in expected {
        let artifact = artifacts
            .iter()
            .find(|artifact| artifact.id == id)
            .ok_or_else(|| {
                effective_error(
                    effective_path,
                    format!("falta el artefacto VCF esperado {id}"),
                )
            })?;
        if artifact.kind != kind || artifact.path != relative_path {
            return Err(effective_error(
                effective_path,
                format!("el artefacto {id} no tiene clase o ruta VCF canónica"),
            ));
        }
        let physical_path = root.join(relative_path);
        let identity = identify_file(&physical_path, Path::new(relative_path))?;
        if identity.size_bytes() != artifact.size_bytes || identity.digest() != artifact.digest {
            return Err(Diagnostic::io(
                "GPS-IO-002",
                "El tamaño o SHA-256 de un artefacto no coincide con el valor esperado.",
            )
            .with("path", physical_path.display().to_string())
            .with("artifact_id", id));
        }
    }
    Ok(())
}

fn validate_effective_targets_against_configuration(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    artifact_paths: &BTreeMap<&str, &str>,
    manifest_path: &Path,
    manifest_bytes: &[u8],
    effective_path: &Path,
) -> Result<EffectiveConfig, Diagnostic> {
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let configuration_id = manifest["configuration_input"]["artifact_id"]
        .as_str()
        .ok_or_else(|| effective_error(effective_path, "configuration_input no enlaza TOML"))?;
    let relative = artifact_paths.get(configuration_id).ok_or_else(|| {
        effective_error(
            effective_path,
            format!("TOML físico ausente: {configuration_id}"),
        )
    })?;
    let configuration_path = root.join(relative);
    let bytes = std::fs::read(&configuration_path)
        .map_err(|error| io_error(&configuration_path, "read_configuration", error))?;
    let digest = configuration_digest(&bytes)
        .map_err(|error| physical_error(effective_path, "configuration_digest", &error))?;
    let original_configuration_path = manifest["configuration_input"]["original_path"]
        .as_str()
        .ok_or_else(|| effective_error(effective_path, "configuration_input sin procedencia"))?;
    let admitted = admit_toml_header(
        Path::new(original_configuration_path),
        bytes,
        ExpectedConfigKind::General,
    )
    .map_err(|error| physical_error(effective_path, "admit_configuration", &error))?;
    let declared = parse_declared_config(admitted)
        .map_err(|error| physical_error(effective_path, "parse_configuration", &error))?;
    let seed = effective_object(effective, "seed", effective_path)?;
    let root_seed = RootSeed::parse(
        "effective-config.seed.root",
        effective_non_empty_text(seed, "root", effective_path)?,
    )
    .map_err(|error| physical_error(effective_path, "root_seed", &error))?;
    let mut entropy = ReplayEntropy(root_seed);
    let expected = build_effective_config_with_entropy(&declared, &mut entropy)
        .map_err(|error| physical_error(effective_path, "rebuild_effective_config", &error))?;
    let expected_execution = derive_execution_id(root_seed, &digest)
        .map_err(|error| physical_error(effective_path, "derive_execution_id", &error))?;
    if manifest.get("execution_id").and_then(Value::as_str) != Some(expected_execution.as_str()) {
        return Err(effective_error(
            effective_path,
            "execution_id no coincide con la semilla y los bytes TOML",
        ));
    }

    if seed.get("root").and_then(Value::as_str)
        != Some(expected.reproducibility().root_seed().to_hex().as_str())
        || seed.get("origin").and_then(Value::as_str)
            != Some(expected.reproducibility().seed_origin().as_str())
    {
        return Err(effective_error(
            effective_path,
            "la semilla efectiva no coincide con la declaración TOML",
        ));
    }

    if effective.get("ploidy").and_then(Value::as_u64)
        != Some(u64::from(expected.ploidy().copy_count()))
    {
        return Err(effective_error(
            effective_path,
            "la ploidía efectiva no coincide con la declaración TOML",
        ));
    }
    if effective.get("initial_populations") != Some(&effective_initial_populations_json(&expected))
    {
        return Err(effective_error(
            effective_path,
            "las cohortes iniciales efectivas no coinciden exactamente con el TOML",
        ));
    }
    let normalization = effective_object(effective, "normalization", effective_path)?;
    if normalization
        .get("tolerance")
        .and_then(Value::as_f64)
        .map(f64::to_bits)
        != Some(expected.normalization().tolerance().to_bits())
    {
        return Err(effective_error(
            effective_path,
            "la tolerancia efectiva no coincide con la declaración TOML",
        ));
    }

    let reference = effective_object(effective, "reference", effective_path)?;
    let expected_source = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => "simulated",
        EffectiveReferenceSource::Supplied { .. } => "supplied",
    };
    if effective_non_empty_text(reference, "id", effective_path)? != expected.reference_id()
        || reference.get("source").and_then(Value::as_str) != Some(expected_source)
        || reference.get("total_genome_size").and_then(Value::as_u64)
            != Some(expected.total_genome_size())
    {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no coincide con la declaración TOML",
        ));
    }
    let expected_original_path = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => Value::Null,
        EffectiveReferenceSource::Supplied { original_path } => {
            Value::String(original_path.to_string_lossy().into_owned())
        }
    };
    if reference.get("original_fasta_path") != Some(&expected_original_path) {
        return Err(effective_error(
            effective_path,
            "la procedencia FASTA efectiva no coincide con la declaración TOML",
        ));
    }
    let expected_canonical_format = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => canonical_fasta_json(),
        EffectiveReferenceSource::Supplied { .. } => Value::Null,
    };
    if reference.get("canonical_format") != Some(&expected_canonical_format) {
        return Err(effective_error(
            effective_path,
            "el formato FASTA efectivo no corresponde al origen declarado",
        ));
    }

    let expected_composition = expected.composition().map(|composition| {
        let values = composition.values();
        json!({"A": values[0], "C": values[1], "G": values[2], "T": values[3], "N": composition.unknown_base_fraction()})
    });
    if reference.get("composition") != Some(&expected_composition.unwrap_or(Value::Null)) {
        return Err(effective_error(
            effective_path,
            "la composición efectiva no coincide con la declaración TOML",
        ));
    }

    let effective_sequences = effective_array(reference, "sequences", effective_path)?;
    if effective_sequences.len() != expected.sequences().len() {
        return Err(effective_error(
            effective_path,
            "el número de secuencias efectivo no coincide con el TOML",
        ));
    }
    for (observed, expected) in effective_sequences.iter().zip(expected.sequences()) {
        let observed = observed
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        if effective_non_empty_text(observed, "id", effective_path)? != expected.id()
            || observed.get("class").and_then(Value::as_str)
                != Some(sequence_class(expected.class()))
            || observed.get("length").and_then(Value::as_u64) != Some(expected.length())
        {
            return Err(effective_error(
                effective_path,
                "orden, ID, clase o longitud de secuencia no coincide con el TOML",
            ));
        }
    }

    let effective_reductions = effective_array(effective, "reduced_references", effective_path)?;
    if effective_reductions.len() != expected.reduced_references().len() {
        return Err(effective_error(
            effective_path,
            "las reducciones efectivas no coinciden con las declaradas en TOML",
        ));
    }
    for (expected, effective) in expected
        .reduced_references()
        .iter()
        .zip(effective_reductions)
    {
        let effective = effective
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        let target = effective_object(effective, "target", effective_path)?;
        let expected_declared = expected.target().declared_value();
        if effective_non_empty_text(effective, "id", effective_path)? != expected.id()
            || effective_non_empty_text(effective, "source_reference_id", effective_path)?
                != expected.source_reference_id()
            || effective.get("region_size_bases").and_then(Value::as_u64)
                != Some(expected.region_size_bases())
            || target.get("mode").and_then(Value::as_str) != Some(expected.target().mode())
            || target.get("declared").and_then(Value::as_str) != Some(expected_declared.as_str())
            || target.get("bases").and_then(Value::as_u64) != Some(expected.target_bases())
            || effective.get("selection_subseed").and_then(Value::as_str)
                != Some(expected.selection_subseed().to_hex().as_str())
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "la reducción efectiva {} no coincide completamente con el TOML",
                    expected.id()
                ),
            ));
        }
    }
    let expected_direct_events = Value::Array(
        expected
            .direct_events()
            .iter()
            .map(direct_event_json)
            .collect(),
    );
    if effective.get("direct_events") != Some(&expected_direct_events) {
        return Err(effective_error(
            effective_path,
            "direct_events efectivo no coincide exactamente con el TOML resuelto",
        ));
    }
    let expected_technical_event = expected
        .technical_event()
        .map(|event| technical_event_json(&expected, event));
    if effective.get("technical_event") != expected_technical_event.as_ref() {
        return Err(effective_error(
            effective_path,
            "technical_event efectivo no coincide exactamente con el TOML resuelto",
        ));
    }
    let expected_vcf = expected.vcf_output().map(|mode| match mode {
        VcfOutputMode::Text => json!({"mode": "text"}),
        VcfOutputMode::Bgzf => json!({"mode": "bgzf"}),
    });
    if effective.get("vcf_output") != Some(&expected_vcf.unwrap_or(Value::Null)) {
        return Err(effective_error(
            effective_path,
            "vcf_output efectivo no coincide exactamente con el TOML resuelto",
        ));
    }
    Ok(expected)
}

struct ReplayEntropy(RootSeed);

impl EntropySource for ReplayEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        destination.copy_from_slice(self.0.as_bytes());
        Ok(())
    }
}

struct EffectiveReductionLinks {
    id: String,
    source_reference_id: String,
    fasta_artifact_id: String,
    map_artifact_id: String,
    source_fasta_artifact_id: String,
}

fn read_physical_reference(
    root: &Path,
    reference: &serde_json::Map<String, Value>,
    artifact_paths: &BTreeMap<&str, &str>,
    effective_path: &Path,
) -> Result<Reference, Diagnostic> {
    let id = effective_non_empty_text(reference, "id", effective_path)?;
    let sequences = effective_array(reference, "sequences", effective_path)?;
    let artifact_id = effective_non_empty_text(reference, "fasta_artifact_id", effective_path)?;
    read_physical_sequence_catalog(
        root,
        id,
        artifact_id,
        sequences,
        "unknown_base_count",
        artifact_paths,
        effective_path,
    )
}

fn read_physical_sequence_catalog(
    root: &Path,
    reference_id: &str,
    artifact_id: &str,
    sequences: &[Value],
    n_count_field: &str,
    artifact_paths: &BTreeMap<&str, &str>,
    effective_path: &Path,
) -> Result<Reference, Diagnostic> {
    let mut expectations = Vec::with_capacity(sequences.len());
    for sequence in sequences {
        let sequence = sequence
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        expectations.push(SequenceExpectation {
            id: effective_non_empty_text(sequence, "id", effective_path)?,
            class: effective_sequence_class(sequence, effective_path)?,
            length: effective_positive_u64(sequence, "length", effective_path)?,
        });
    }
    let relative = artifact_paths.get(artifact_id).ok_or_else(|| {
        effective_error(
            effective_path,
            format!("FASTA físico ausente: {artifact_id}"),
        )
    })?;
    let path = root.join(relative);
    let observed = validate_supplied_fasta(&path, reference_id, &expectations)
        .map_err(|error| physical_error(effective_path, "validate_reference_fasta", &error))?
        .into_reference();
    for (declared, physical) in sequences.iter().zip(observed.sequences()) {
        let declared = declared
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        let unknown = declared
            .get(n_count_field)
            .and_then(Value::as_u64)
            .ok_or_else(|| effective_error(effective_path, format!("{n_count_field} inválido")))?;
        if unknown != physical.unknown_base_count() {
            return Err(effective_error(
                effective_path,
                format!(
                    "{n_count_field} físico de {} no coincide",
                    physical.id().as_str()
                ),
            ));
        }
        if declared.get("header") != Some(&reversible_bytes_json(physical.header().raw_bytes()))
            || declared.get("description")
                != Some(&reversible_description_json(
                    physical.header().description(),
                ))
            || declared.get("alphabet_validation").and_then(Value::as_str)
                != Some(alphabet_validation(physical.alphabet_validation()))
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "cabecera, descripción o validación de alfabeto física de {} no coincide",
                    physical.id().as_str()
                ),
            ));
        }
    }
    Ok(observed)
}

fn read_reference_with_expected(
    root: &Path,
    reference_id: &str,
    artifact_id: &str,
    declared_sequences: &[Value],
    expected: &[ReferenceSequence],
    artifact_paths: &BTreeMap<&str, &str>,
    effective_path: &Path,
) -> Result<Reference, Diagnostic> {
    let observed = read_physical_sequence_catalog(
        root,
        reference_id,
        artifact_id,
        declared_sequences,
        "n_count",
        artifact_paths,
        effective_path,
    )?;
    if observed.sequences().len() != expected.len() {
        return Err(effective_error(
            effective_path,
            "el catálogo físico derivado no coincide con el plan",
        ));
    }
    Ok(observed)
}

fn validate_generated_reference(
    reference: &serde_json::Map<String, Value>,
    root_seed: RootSeed,
    physical: &Reference,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let sequences = effective_array(reference, "sequences", effective_path)?;
    let mut specs = Vec::with_capacity(sequences.len());
    for sequence in sequences {
        let sequence = sequence
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        specs.push(GeneratedSequenceSpec {
            id: effective_non_empty_text(sequence, "id", effective_path)?,
            class: effective_sequence_class(sequence, effective_path)?,
            length: effective_positive_u64(sequence, "length", effective_path)?,
        });
    }
    let composition = effective_object(reference, "composition", effective_path)?;
    let frequencies = (*b"ACGT")
        .into_iter()
        .zip(["A", "C", "G", "T"])
        .map(|(base, label)| {
            composition
                .get(label)
                .and_then(Value::as_f64)
                .map(|probability| BaseFrequency { base, probability })
                .ok_or_else(|| effective_error(effective_path, "composición efectiva inválida"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let unknown_base_fraction = composition.get("N").map_or(Ok(0.0), |value| {
        value
            .as_f64()
            .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
            .ok_or_else(|| effective_error(effective_path, "fracción N efectiva inválida"))
    })?;
    let expected = crate::reference_generation::generate_reference_with_unknown_fraction(
        effective_non_empty_text(reference, "id", effective_path)?,
        root_seed,
        &specs,
        &frequencies,
        unknown_base_fraction,
    )
    .map_err(|error| physical_error(effective_path, "regenerate_reference", &error))?;
    compare_reference_bases(&expected, physical, effective_path, expected.id().as_str())
}

fn compare_reference_bases(
    expected: &Reference,
    physical: &Reference,
    effective_path: &Path,
    reference_id: &str,
) -> Result<(), Diagnostic> {
    if expected.id() != physical.id()
        || expected.sequences().len() != physical.sequences().len()
        || expected
            .sequences()
            .iter()
            .zip(physical.sequences())
            .any(|(expected, physical)| {
                expected.id() != physical.id()
                    || expected.class() != physical.class()
                    || expected.length() != physical.length()
                    || expected.unknown_base_count() != physical.unknown_base_count()
                    || expected.sequence().as_bytes() != physical.sequence().as_bytes()
            })
    {
        return Err(effective_error(
            effective_path,
            format!("FASTA físico de {reference_id} no coincide base a base"),
        ));
    }
    Ok(())
}

fn effective_sequence_class(
    sequence: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<SequenceClass, Diagnostic> {
    match sequence.get("class").and_then(Value::as_str) {
        Some("chromosome") => Ok(SequenceClass::Chromosome),
        Some("scaffold") => Ok(SequenceClass::Scaffold),
        _ => Err(effective_error(path, "clase de secuencia inválida")),
    }
}

fn physical_error(path: &Path, operation: &str, error: &Diagnostic) -> Diagnostic {
    effective_error(
        path,
        format!(
            "autoridad física rechazada en {operation}: {} {}",
            error.code(),
            error.message()
        ),
    )
}

pub(crate) fn parse_general_manifest_bytes(
    path: &Path,
    bytes: &[u8],
) -> Result<Vec<ArtifactRecord>, Diagnostic> {
    let schema_version = serde_json::from_slice::<Value>(bytes)
        .ok()
        .and_then(|document| {
            document
                .get("schema_version")
                .and_then(Value::as_str)
                .map(str::to_owned)
        });
    parse_general_manifest_bytes_unretagged(path, bytes).map_err(|error| {
        if schema_version.as_deref() == Some(GENERAL_MANIFEST_SCHEMA_VERSION) {
            retag_format_schema(error, "manifest-general-v4.schema.json")
        } else {
            error
        }
    })
}

fn parse_general_manifest_bytes_unretagged(
    path: &Path,
    bytes: &[u8],
) -> Result<Vec<ArtifactRecord>, Diagnostic> {
    let document: Value = serde_json::from_slice(bytes)
        .map_err(|error| manifest_error(path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| manifest_error(path, "la raíz no es un objeto"))?;
    let schema_version = object.get("schema_version").and_then(Value::as_str);
    match schema_version {
        Some(GENERAL_MANIFEST_SCHEMA_VERSION) => {}
        observed => {
            return Err(manifest_error_for_schema(
                path,
                "manifest-general-v4.schema.json",
                format!(
                    "versión de manifiesto no admitida: {:?}; se admite exclusivamente \"{}\"",
                    observed, GENERAL_MANIFEST_SCHEMA_VERSION
                ),
            ));
        }
    }
    const FIELDS: [&str; 16] = [
        "manifest_kind",
        "schema_version",
        "execution_id",
        "state",
        "configuration_input",
        "effective_config_artifact_id",
        "references",
        "reference_reductions",
        "population_snvs",
        "vcf_outputs",
        "authorizations",
        "output_categories",
        "artifacts",
        "diagnostics",
        "events",
        "technical_events",
    ];
    validate_closed_object(object, &FIELDS, "manifest.json", path)?;
    if object.get("manifest_kind").and_then(Value::as_str) != Some("gps_execution_manifest")
        || schema_version != Some(GENERAL_MANIFEST_SCHEMA_VERSION)
    {
        return Err(manifest_error(path, "clase o versión incompatible"));
    }
    required_non_empty_text(object, "execution_id", path)?;
    let state = match object.get("state").and_then(Value::as_str) {
        Some("complete") => UnitState::Complete,
        Some("failed") => UnitState::Failed,
        _ => return Err(manifest_error(path, "estado de ejecución inválido")),
    };

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
    let catalog = records
        .iter()
        .map(|record| (record.id.as_str(), record.kind.as_str()))
        .collect::<BTreeMap<_, _>>();
    if state == UnitState::Complete {
        require_truth_artifact(&records, path)?;
    }

    let configuration_input = object
        .get("configuration_input")
        .ok_or_else(|| manifest_error(path, "falta configuration_input"))?;
    let configuration_artifact_id = match configuration_input {
        Value::Null => None,
        Value::Object(value) => Some(validate_configuration_input(value, &catalog, path)?),
        _ => return Err(manifest_error(path, "configuration_input inválido")),
    };
    let effective_id = optional_text(object, "effective_config_artifact_id", path)?;
    if let Some(id) = effective_id {
        require_artifact_kind(&catalog, id, "effective_config", path)?;
    }
    let references =
        validate_references(required_array(object, "references", path)?, &catalog, path)?;
    let vcf_requested = validate_vcf_outputs(
        required_array(object, "vcf_outputs", path)?,
        &catalog,
        &references,
        path,
    )?;
    let reductions = validate_reference_reductions(
        required_array(object, "reference_reductions", path)?,
        &catalog,
        &references,
        path,
    )?;
    let population_snvs =
        validate_manifest_population_snvs(object.get("population_snvs"), &catalog, path)?;
    let events = validate_manifest_events(object.get("events"), &records, &catalog, path)?;
    validate_manifest_technical_events(
        required_array(object, "technical_events", path)?,
        &events,
        effective_id,
        &records,
        path,
    )?;
    validate_manifest_authorizations(object.get("authorizations"), path)?;
    if let Some(id) = population_snvs.as_deref() {
        let declared = object
            .get("population_snvs")
            .and_then(Value::as_object)
            .and_then(|link| link.get("metadata_sha256"))
            .and_then(Value::as_str);
        let observed = records
            .iter()
            .find(|artifact| artifact.id == id)
            .map(|artifact| artifact.digest.as_str());
        if declared != observed {
            return Err(manifest_error(
                path,
                "population_snvs.metadata_sha256 no coincide con el artefacto",
            ));
        }
    }
    let mut expected_minimum = BTreeSet::new();
    if let Some(id) = configuration_artifact_id.as_deref() {
        expected_minimum.insert(id);
    }
    if let Some(id) = effective_id {
        expected_minimum.insert(id);
    }
    if let Some(id) = population_snvs.as_deref() {
        expected_minimum.insert(id);
        require_artifact_kind(&catalog, "population-snvs-tsv", "population_snvs_tsv", path)?;
        expected_minimum.insert("population-snvs-tsv");
    }
    if state == UnitState::Complete {
        expected_minimum.insert("truth-jsonl");
    }
    if !events.is_empty() {
        expected_minimum.extend(
            records
                .iter()
                .filter(|artifact| artifact.kind == "fastq" || artifact.kind == "event_manifest")
                .map(|artifact| artifact.id.as_str()),
        );
    }
    for (reference_id, artifact_id) in &references {
        if !reductions.derived_reference_ids.contains(reference_id) {
            expected_minimum.insert(artifact_id.as_str());
        }
    }
    validate_output_categories(
        required_object(object, "output_categories", path)?,
        &catalog,
        &expected_minimum,
        &reductions
            .requested_artifact_ids
            .union(&vcf_requested)
            .cloned()
            .collect(),
        path,
    )?;
    let diagnostics = required_array(object, "diagnostics", path)?;
    validate_diagnostics(diagnostics, path)?;

    match state {
        UnitState::Complete => {
            if configuration_input.is_null()
                || effective_id.is_none()
                || required_array(object, "references", path)?.is_empty()
                || reductions
                    .states
                    .iter()
                    .any(|state| *state != UnitState::Complete)
                || !diagnostics.is_empty()
                || events
                    .iter()
                    .any(|(_, state, _)| *state != UnitState::Complete)
            {
                return Err(manifest_error(
                    path,
                    "una ejecución complete exige entrada, configuración efectiva, referencia y ningún diagnóstico",
                ));
            }
        }
        UnitState::Failed => {
            if diagnostics.is_empty() {
                return Err(manifest_error(
                    path,
                    "una ejecución failed exige al menos un diagnóstico",
                ));
            }
        }
        UnitState::Incomplete => unreachable!("manifest parser never constructs incomplete"),
    }
    Ok(records)
}

/// Reads the frozen F1 manifest contract.  It deliberately does not route a
/// v1 document through F2 requirements such as `truth.jsonl`, population-SNV
/// metadata, authorizations, or VCF outputs.
#[allow(dead_code)]
fn parse_general_manifest_v1_object(
    path: &Path,
    object: &serde_json::Map<String, Value>,
) -> Result<Vec<ArtifactRecord>, Diagnostic> {
    const FIELDS: [&str; 11] = [
        "manifest_kind",
        "schema_version",
        "execution_id",
        "state",
        "configuration_input",
        "effective_config_artifact_id",
        "references",
        "reference_reductions",
        "output_categories",
        "artifacts",
        "diagnostics",
    ];
    validate_closed_object(object, &FIELDS, "manifest.json", path)?;
    if object.get("manifest_kind").and_then(Value::as_str) != Some("gps_execution_manifest")
        || object.get("schema_version").and_then(Value::as_str) != Some("1")
    {
        return Err(manifest_error(path, "clase o versión incompatible"));
    }
    required_non_empty_text(object, "execution_id", path)?;
    let state = match object.get("state").and_then(Value::as_str) {
        Some("complete") => UnitState::Complete,
        Some("failed") => UnitState::Failed,
        _ => return Err(manifest_error(path, "estado de ejecución inválido")),
    };

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
    let catalog = records
        .iter()
        .map(|record| (record.id.as_str(), record.kind.as_str()))
        .collect::<BTreeMap<_, _>>();
    let configuration_input = object
        .get("configuration_input")
        .ok_or_else(|| manifest_error(path, "falta configuration_input"))?;
    let configuration_artifact_id = match configuration_input {
        Value::Null => None,
        Value::Object(value) => Some(validate_configuration_input(value, &catalog, path)?),
        _ => return Err(manifest_error(path, "configuration_input inválido")),
    };
    let effective_id = optional_text(object, "effective_config_artifact_id", path)?;
    if let Some(id) = effective_id {
        require_artifact_kind(&catalog, id, "effective_config", path)?;
    }
    let references =
        validate_references(required_array(object, "references", path)?, &catalog, path)?;
    let reductions = validate_reference_reductions(
        required_array(object, "reference_reductions", path)?,
        &catalog,
        &references,
        path,
    )?;
    let mut expected_minimum = BTreeSet::new();
    if let Some(id) = configuration_artifact_id.as_deref() {
        expected_minimum.insert(id);
    }
    if let Some(id) = effective_id {
        expected_minimum.insert(id);
    }
    for (reference_id, artifact_id) in &references {
        if !reductions.derived_reference_ids.contains(reference_id) {
            expected_minimum.insert(artifact_id.as_str());
        }
    }
    validate_output_categories(
        required_object(object, "output_categories", path)?,
        &catalog,
        &expected_minimum,
        &reductions.requested_artifact_ids,
        path,
    )?;
    let diagnostics = required_array(object, "diagnostics", path)?;
    validate_diagnostics(diagnostics, path)?;
    match state {
        UnitState::Complete
            if configuration_input.is_null()
                || effective_id.is_none()
                || required_array(object, "references", path)?.is_empty()
                || reductions
                    .states
                    .iter()
                    .any(|state| *state != UnitState::Complete)
                || !diagnostics.is_empty() =>
        {
            Err(manifest_error(
                path,
                "una ejecución complete exige entrada, configuración efectiva, referencia y ningún diagnóstico",
            ))
        }
        UnitState::Failed if diagnostics.is_empty() => Err(manifest_error(
            path,
            "una ejecución failed exige al menos un diagnóstico",
        )),
        UnitState::Incomplete => unreachable!("manifest parser never constructs incomplete"),
        _ => Ok(records),
    }
}

#[allow(dead_code)]
fn validate_effective_config_v1_against_manifest_bytes(
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    let artifacts = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let manifest = manifest
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "la raíz no es un objeto"))?;
    if manifest.get("schema_version").and_then(Value::as_str) != Some("1") {
        return Err(effective_error(
            effective_path,
            "el manifiesto no usa la misma versión de contrato que la configuración efectiva",
        ));
    }
    let catalog = artifacts
        .iter()
        .map(|record| (record.id.as_str(), record.kind.as_str()))
        .collect::<BTreeMap<_, _>>();
    let effective_id = optional_text(manifest, "effective_config_artifact_id", manifest_path)?
        .ok_or_else(|| {
            effective_error(
                effective_path,
                "el manifiesto no enlaza la configuración efectiva",
            )
        })?;
    require_artifact_kind(&catalog, effective_id, "effective_config", manifest_path)?;
    let manifest_references = validate_references(
        required_array(manifest, "references", manifest_path)?,
        &catalog,
        manifest_path,
    )?;
    let document: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "la raíz no es un objeto"))?;
    const FIELDS: [&str; 6] = [
        "effective_config_kind",
        "schema_version",
        "source_path",
        "seed",
        "reference",
        "reduced_references",
    ];
    validate_effective_closed_object(object, &FIELDS, "effective-config.json", effective_path)?;
    if object.get("effective_config_kind").and_then(Value::as_str) != Some("gps_effective_config")
        || object.get("schema_version").and_then(Value::as_str) != Some("1")
    {
        return Err(effective_error(
            effective_path,
            "clase o versión incompatible",
        ));
    }
    let source_path = effective_non_empty_text(object, "source_path", effective_path)?;
    let configuration_input = required_object(manifest, "configuration_input", manifest_path)?;
    if required_non_empty_text(configuration_input, "original_path", manifest_path)? != source_path
    {
        return Err(effective_error(
            effective_path,
            "source_path no coincide con la procedencia de la copia TOML",
        ));
    }
    let seed = effective_object(object, "seed", effective_path)?;
    validate_effective_seed(seed, effective_path)?;
    let effective_root_seed = effective_non_empty_text(seed, "root", effective_path)?;
    let effective_reference = effective_object(object, "reference", effective_path)?;
    let reference = validate_effective_reference(effective_reference, &catalog, effective_path)?;
    let eligible_bases =
        effective_positive_u64(effective_reference, "total_genome_size", effective_path)?;
    if manifest_references.get(&reference.0).map(String::as_str) != Some(reference.1.as_str()) {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no coincide con el enlace del manifiesto",
        ));
    }
    let mut effective_references = BTreeMap::from([(reference.0.clone(), reference.1.clone())]);
    let manifest_reductions = required_array(manifest, "reference_reductions", manifest_path)?
        .iter()
        .map(|value| {
            let reduction = value
                .as_object()
                .ok_or_else(|| manifest_error(manifest_path, "una reducción no es objeto"))?;
            let id = required_non_empty_text(reduction, "id", manifest_path)?;
            Ok((id.to_owned(), reduction))
        })
        .collect::<Result<BTreeMap<_, _>, Diagnostic>>()?;
    let effective_reductions = effective_array(object, "reduced_references", effective_path)?;
    if effective_reductions.len() != manifest_reductions.len() {
        return Err(effective_error(
            effective_path,
            "las reducciones efectivas no coinciden con el manifiesto",
        ));
    }
    for value in effective_reductions {
        let reduction = validate_effective_reduction(value, &catalog, effective_path)?;
        let manifest_reduction = manifest_reductions.get(&reduction.id).ok_or_else(|| {
            effective_error(
                effective_path,
                format!("la reducción {} no existe en el manifiesto", reduction.id),
            )
        })?;
        let effective_reduction = value
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        if manifest_reduction.get("state").and_then(Value::as_str) != Some("complete")
            || manifest_reduction
                .get("eligible_bases")
                .and_then(Value::as_u64)
                != Some(eligible_bases)
            || manifest_reduction.get("root_seed").and_then(Value::as_str)
                != Some(effective_root_seed)
            || manifest_reduction
                .get("source_reference_id")
                .and_then(Value::as_str)
                != Some(reduction.source_reference_id.as_str())
            || manifest_reduction
                .get("fasta_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.fasta_artifact_id.as_str())
            || manifest_reduction
                .get("map_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.map_artifact_id.as_str())
            || manifest_reduction
                .get("source_fasta_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.source_fasta_artifact_id.as_str())
            || [
                "region_size_bases",
                "target",
                "selection_subseed",
                "selection_order",
            ]
            .iter()
            .any(|field| manifest_reduction.get(*field) != effective_reduction.get(*field))
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "la reducción efectiva {} no coincide con el manifiesto",
                    reduction.id
                ),
            ));
        }
        if effective_references
            .insert(reduction.id.clone(), reduction.fasta_artifact_id.clone())
            .is_some()
        {
            return Err(effective_error(
                effective_path,
                "una reducción reutiliza un ID de referencia",
            ));
        }
    }
    if effective_references != manifest_references {
        return Err(effective_error(
            effective_path,
            "el conjunto de referencias del manifiesto no coincide exactamente con la configuración efectiva",
        ));
    }
    Ok(())
}

fn require_truth_artifact(records: &[ArtifactRecord], path: &Path) -> Result<(), Diagnostic> {
    let Some(truth) = records.iter().find(|artifact| artifact.id == "truth-jsonl") else {
        return Err(manifest_error(
            path,
            "falta el artefacto obligatorio truth-jsonl",
        ));
    };
    if truth.kind != "truth_jsonl" || truth.path != "outputs/truth.jsonl" {
        return Err(manifest_error(
            path,
            "truth-jsonl debe tener clase truth_jsonl y ruta outputs/truth.jsonl",
        ));
    }
    Ok(())
}

fn validate_output_categories(
    object: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    expected_minimum: &BTreeSet<&str>,
    expected_requested: &BTreeSet<String>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 4] = ["minimum", "defaulted", "requested", "disabled"];
    validate_closed_object(object, &FIELDS, "output_categories", path)?;
    let mut categories = BTreeMap::new();
    let mut globally_classified = BTreeSet::new();
    for field in FIELDS {
        let mut ids = BTreeSet::new();
        for value in required_array(object, field, path)? {
            let id = value
                .as_str()
                .filter(|id| !id.is_empty())
                .ok_or_else(|| manifest_error(path, format!("{field} contiene un ID inválido")))?;
            if !catalog.contains_key(id) {
                return Err(manifest_error(
                    path,
                    format!("la categoría {field} enlaza un artefacto ausente: {id}"),
                ));
            }
            if !ids.insert(id) {
                return Err(manifest_error(
                    path,
                    format!("el artefacto {id} se repite en {field}"),
                ));
            }
            if !globally_classified.insert(id) {
                return Err(manifest_error(
                    path,
                    format!("el artefacto {id} aparece en más de una categoría"),
                ));
            }
        }
        categories.insert(field, ids);
    }
    if !categories["defaulted"].is_empty() || !categories["disabled"].is_empty() {
        return Err(manifest_error(
            path,
            "defaulted y disabled deben estar vacías en F1",
        ));
    }
    let expected_requested = expected_requested
        .iter()
        .map(String::as_str)
        .collect::<BTreeSet<_>>();
    if !expected_minimum.is_disjoint(&expected_requested) {
        return Err(manifest_error(
            path,
            "un mismo artefacto no puede tener simultáneamente rol minimum y requested",
        ));
    }
    if categories["minimum"] != *expected_minimum || categories["requested"] != expected_requested {
        return Err(manifest_error(
            path,
            "minimum/requested no coinciden con los roles normativos de F1",
        ));
    }
    let expected_catalog = expected_minimum
        .iter()
        .copied()
        .chain(expected_requested.iter().copied())
        .collect::<BTreeSet<_>>();
    if expected_catalog != catalog.keys().copied().collect::<BTreeSet<_>>() {
        return Err(manifest_error(
            path,
            "el catálogo contiene un artefacto sin rol de salida permitido en F1",
        ));
    }
    Ok(())
}

fn validate_manifest_events(
    value: Option<&Value>,
    records: &[ArtifactRecord],
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<Vec<(String, UnitState, String)>, Diagnostic> {
    let Some(value) = value else {
        return Ok(Vec::new());
    };
    let values = value
        .as_array()
        .ok_or_else(|| manifest_error(path, "events no es una lista"))?;
    let mut previous = None;
    let mut links = Vec::with_capacity(values.len());
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "un enlace de evento no es objeto"))?;
        const FIELDS: [&str; 3] = ["event_id", "state", "manifest_artifact_id"];
        validate_closed_object(object, &FIELDS, "events[]", path)?;
        let event_id = required_non_empty_text(object, "event_id", path)?;
        EventId::parse(event_id.to_owned(), "manifest.events")
            .map_err(|error| manifest_error(path, error.render_human()))?;
        if previous.is_some_and(|previous: &str| previous >= event_id) {
            return Err(manifest_error(path, "events no está ordenado por event_id"));
        }
        previous = Some(event_id);
        let state = match required_non_empty_text(object, "state", path)? {
            "complete" => UnitState::Complete,
            "failed" => UnitState::Failed,
            _ => return Err(manifest_error(path, "estado de evento inválido")),
        };
        let artifact_id = required_non_empty_text(object, "manifest_artifact_id", path)?;
        require_artifact_kind(catalog, artifact_id, "event_manifest", path)?;
        let artifact = records
            .iter()
            .find(|artifact| artifact.id == artifact_id)
            .ok_or_else(|| manifest_error(path, "manifiesto de evento ausente"))?;
        let expected_path = event_manifest_relative_path(
            &EventId::parse(event_id.to_owned(), "manifest.events")
                .map_err(|error| manifest_error(path, error.render_human()))?,
        );
        if artifact.path != expected_path.to_string_lossy() {
            return Err(manifest_error(
                path,
                "la ruta del manifiesto de evento no es canónica",
            ));
        }
        links.push((event_id.to_owned(), state, artifact_id.to_owned()));
    }
    Ok(links)
}

fn validate_manifest_technical_events(
    values: &[Value],
    events: &[(String, UnitState, String)],
    effective_id: Option<&str>,
    records: &[ArtifactRecord],
    path: &Path,
) -> Result<(), Diagnostic> {
    if values.len() != events.len() {
        return Err(manifest_error(
            path,
            "technical_events no coincide con events",
        ));
    }
    for (value, (event_id, _, manifest_id)) in values.iter().zip(events) {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "technical_event inválido"))?;
        const FIELDS: [&str; 7] = [
            "event_id",
            "effective_config",
            "profile",
            "rng_catalog",
            "id_catalog",
            "artifacts",
            "event_manifest",
        ];
        validate_closed_object(object, &FIELDS, "technical_events[]", path)?;
        if required_non_empty_text(object, "event_id", path)? != event_id
            || object.get("profile")
                != Some(
                    &json!({"id": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_V1, "sha256": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256}),
                )
            || object.get("rng_catalog").and_then(Value::as_str)
                != Some(crate::randomness::F4_TECHNICAL_RNG_CATALOG)
            || object.get("id_catalog") != Some(&f4_id_catalog_value())
        {
            return Err(manifest_error(path, "inventario técnico F4 incompatible"));
        }
        validate_embedded_artifact(object.get("effective_config"), effective_id, records, path)?;
        validate_embedded_artifact(
            object.get("event_manifest"),
            Some(manifest_id),
            records,
            path,
        )?;
        let artifacts = object
            .get("artifacts")
            .and_then(Value::as_object)
            .ok_or_else(|| manifest_error(path, "artefactos técnicos agregados inválidos"))?;
        validate_closed_object(
            artifacts,
            &["chimera"],
            "technical_events[].artifacts",
            path,
        )?;
        let chimera = artifacts
            .get("chimera")
            .and_then(Value::as_object)
            .ok_or_else(|| manifest_error(path, "resumen de quimera inválido"))?;
        validate_closed_object(
            chimera,
            &["rate", "realized_pairs"],
            "technical_events[].artifacts.chimera",
            path,
        )?;
        let rate = chimera.get("rate").and_then(Value::as_f64);
        let realized_pairs = chimera.get("realized_pairs").and_then(Value::as_u64);
        if !matches!(rate, Some(rate) if rate.is_finite() && (0.0..=1.0).contains(&rate))
            || realized_pairs.is_none()
            || (rate == Some(0.0) && realized_pairs != Some(0))
        {
            return Err(manifest_error(path, "resumen de quimera inválido"));
        }
    }
    Ok(())
}

fn validate_embedded_artifact(
    value: Option<&Value>,
    expected_id: Option<&str>,
    records: &[ArtifactRecord],
    path: &Path,
) -> Result<(), Diagnostic> {
    match (value, expected_id) {
        (Some(Value::Null), None) => Ok(()),
        (Some(value), Some(expected_id))
            if value
                == &records
                    .iter()
                    .find(|record| record.id == expected_id)
                    .map_or(Value::Null, ArtifactRecord::to_json) =>
        {
            Ok(())
        }
        _ => Err(manifest_error(
            path,
            "artefacto técnico embebido no coincide con el catálogo",
        )),
    }
}

fn validate_manifest_population_snvs(
    value: Option<&Value>,
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<Option<String>, Diagnostic> {
    let Some(value) = value else { return Ok(None) };
    if value.is_null() {
        return Ok(None);
    }
    let object = value
        .as_object()
        .ok_or_else(|| manifest_error(path, "population_snvs inválido"))?;
    const FIELDS: [&str; 5] = [
        "mode",
        "metadata_sha256",
        "metadata_artifact_id",
        "normalization_audit",
        "initial_population_statistics",
    ];
    let fields = if object.contains_key("initial_population_statistics") {
        &FIELDS[..]
    } else {
        &FIELDS[..4]
    };
    validate_closed_object(object, fields, "population_snvs", path)?;
    let mode = match object.get("mode").and_then(Value::as_str) {
        Some(mode @ ("supplied" | "generated" | "reuse")) => mode,
        _ => return Err(manifest_error(path, "population_snvs.mode inválido")),
    };
    let digest = required_non_empty_text(object, "metadata_sha256", path)?;
    if digest.len() != 64
        || !digest
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        return Err(manifest_error(
            path,
            "population_snvs.metadata_sha256 inválido",
        ));
    }
    validate_serialized_normalization_audit(
        object.get("normalization_audit"),
        "population_snvs.normalization_audit inválida",
    )
    .map_err(|detail| manifest_error(path, detail))?;
    if let Some(statistics) = object.get("initial_population_statistics") {
        validate_initial_population_statistics(statistics)
            .map_err(|detail| manifest_error(path, detail))?;
    }
    match (mode, object.get("metadata_artifact_id")) {
        ("reuse", None | Some(Value::Null)) => Ok(None),
        ("supplied" | "generated", Some(Value::String(id))) if !id.is_empty() => {
            require_artifact_kind(catalog, id, "population_snvs_metadata", path)?;
            Ok(Some(id.clone()))
        }
        _ => Err(manifest_error(
            path,
            "population_snvs.metadata_artifact_id inválido",
        )),
    }
}

fn validate_initial_population_statistics(value: &Value) -> Result<(), &'static str> {
    let object = value
        .as_object()
        .ok_or("estadísticas iniciales inválidas")?;
    const FIELDS: [&str; 4] = ["kind", "version", "provenance", "populations"];
    if object.len() != FIELDS.len() || object.keys().any(|field| !FIELDS.contains(&field.as_str()))
    {
        return Err("estadísticas iniciales inválidas");
    }
    if object.get("kind").and_then(Value::as_str) != Some("initial_population_statistics")
        || object.get("version").and_then(Value::as_str) != Some("1")
    {
        return Err("estadísticas iniciales inválidas");
    }
    let provenance = object
        .get("provenance")
        .and_then(Value::as_object)
        .ok_or("procedencia de estadísticas iniciales inválida")?;
    validate_initial_population_statistics_provenance(provenance)?;
    let populations = object
        .get("populations")
        .and_then(Value::as_array)
        .ok_or("estadísticas iniciales inválidas")?;
    if populations.is_empty() {
        return Err("estadísticas iniciales sin poblaciones");
    }
    let mut population_ids = BTreeSet::new();
    for population in populations {
        validate_initial_population_statistics_population(population)?;
        let population_id = population
            .get("population_id")
            .and_then(Value::as_str)
            .ok_or("población estadística inválida")?;
        if !population_ids.insert(population_id) {
            return Err("población estadística repetida");
        }
    }
    Ok(())
}

fn validate_initial_population_statistics_provenance(
    provenance: &serde_json::Map<String, Value>,
) -> Result<(), &'static str> {
    const FIELDS: [&str; 6] = [
        "matrix_reference_id",
        "matrix_reference_fingerprint",
        "matrix_canonical_tsv_sha256",
        "genotype_algorithms",
        "missing_locus_policy",
        "explicit_snp_constraints",
    ];
    if provenance.len() != FIELDS.len()
        || provenance
            .keys()
            .any(|field| !FIELDS.contains(&field.as_str()))
    {
        return Err("procedencia de estadísticas iniciales inválida");
    }
    let reference_id = provenance
        .get("matrix_reference_id")
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty());
    if reference_id.is_none()
        || !valid_lower_sha256(provenance.get("matrix_reference_fingerprint"))
        || !valid_lower_sha256(provenance.get("matrix_canonical_tsv_sha256"))
        || provenance
            .get("missing_locus_policy")
            .and_then(Value::as_str)
            != Some("missing_by_deletion_excluded_from_allele_frequency_denominator_v1")
    {
        return Err("procedencia de estadísticas iniciales inválida");
    }
    let algorithms = provenance
        .get("genotype_algorithms")
        .and_then(Value::as_array)
        .ok_or("algoritmos de genotipo inválidos")?;
    if algorithms.len() != 3
        || algorithms[0].as_str() != Some("hwe-composition-v1")
        || algorithms[1].as_str() != Some("f-diploid-genotype-v1")
        || algorithms[2].as_str() != Some("phase-permutation-v1")
    {
        return Err("algoritmos de genotipo inválidos");
    }
    let constraints = provenance
        .get("explicit_snp_constraints")
        .and_then(Value::as_array)
        .ok_or("restricciones SNV explícitas inválidas")?;
    let mut constraint_ids = BTreeSet::new();
    for constraint in constraints {
        let constraint = constraint
            .as_object()
            .ok_or("restricción SNV explícita inválida")?;
        const CONSTRAINT_FIELDS: [&str; 8] = [
            "population_id",
            "individual_id",
            "event_id",
            "phase",
            "matrix_variant_id",
            "sampled_alternate",
            "explicit_alternate",
            "incompatible",
        ];
        if constraint.len() != CONSTRAINT_FIELDS.len()
            || constraint
                .keys()
                .any(|field| !CONSTRAINT_FIELDS.contains(&field.as_str()))
            || ![
                "population_id",
                "individual_id",
                "event_id",
                "matrix_variant_id",
            ]
            .iter()
            .all(|field| {
                constraint
                    .get(*field)
                    .and_then(Value::as_str)
                    .is_some_and(|value| !value.is_empty())
            })
            || !matches!(constraint.get("phase").and_then(Value::as_u64), Some(0..=9))
            || !constraint.get("sampled_alternate").is_some_and(|value| {
                value.is_null() || value.as_str().is_some_and(is_canonical_base_text)
            })
            || !constraint
                .get("explicit_alternate")
                .and_then(Value::as_str)
                .is_some_and(is_canonical_base_text)
            || !constraint
                .get("incompatible")
                .is_some_and(Value::is_boolean)
        {
            return Err("restricción SNV explícita inválida");
        }
        let key = (
            constraint["population_id"].as_str().unwrap_or_default(),
            constraint["individual_id"].as_str().unwrap_or_default(),
            constraint["event_id"].as_str().unwrap_or_default(),
        );
        if !constraint_ids.insert(key) {
            return Err("restricción SNV explícita repetida");
        }
    }
    Ok(())
}

fn validate_initial_population_statistics_population(value: &Value) -> Result<(), &'static str> {
    let population = value.as_object().ok_or("población estadística inválida")?;
    const FIELDS: [&str; 5] = [
        "population_id",
        "genotype_model",
        "individual_count",
        "ploidy",
        "loci",
    ];
    if population.len() != FIELDS.len()
        || population.keys().any(|key| !FIELDS.contains(&key.as_str()))
    {
        return Err("población estadística inválida");
    }
    let population_id = population
        .get("population_id")
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty());
    let individual_count = population.get("individual_count").and_then(Value::as_u64);
    let ploidy = population.get("ploidy").and_then(Value::as_u64);
    if population_id.is_none()
        || individual_count == Some(0)
        || individual_count.is_none()
        || !matches!(ploidy, Some(1..=10))
    {
        return Err("población estadística inválida");
    }
    validate_initial_population_statistics_model(
        population
            .get("genotype_model")
            .and_then(Value::as_object)
            .ok_or("modelo genotípico estadístico inválido")?,
        ploidy.unwrap_or_default(),
    )?;
    let loci = population
        .get("loci")
        .and_then(Value::as_array)
        .ok_or("loci estadísticos inválidos")?;
    let mut locus_keys = BTreeSet::new();
    for locus in loci {
        validate_initial_population_statistics_locus(
            locus,
            individual_count.unwrap_or_default(),
            ploidy.unwrap_or_default(),
        )?;
        let key = (
            locus
                .get("contig")
                .and_then(Value::as_str)
                .unwrap_or_default(),
            locus
                .get("position_1based")
                .and_then(Value::as_u64)
                .unwrap_or_default(),
            locus.get("ref").and_then(Value::as_str).unwrap_or_default(),
        );
        if !locus_keys.insert(key) {
            return Err("locus estadístico repetido");
        }
    }
    Ok(())
}

fn validate_initial_population_statistics_model(
    model: &serde_json::Map<String, Value>,
    ploidy: u64,
) -> Result<(), &'static str> {
    let kind = model.get("kind").and_then(Value::as_str);
    match kind {
        Some("hwe") if model.len() == 1 => Ok(()),
        Some("f_diploid")
            if ploidy == 2
                && model.len() == 2
                && model
                    .get("inbreeding_coefficient")
                    .and_then(Value::as_f64)
                    .is_some_and(|value| value.is_finite() && (0.0..=1.0).contains(&value)) =>
        {
            Ok(())
        }
        _ => Err("modelo genotípico estadístico inválido"),
    }
}

fn validate_initial_population_statistics_locus(
    value: &Value,
    individual_count: u64,
    ploidy: u64,
) -> Result<(), &'static str> {
    let locus = value.as_object().ok_or("locus estadístico inválido")?;
    const FIELDS: [&str; 9] = [
        "contig",
        "position_1based",
        "ref",
        "total_copies",
        "observed_copies",
        "missing_by_deletion",
        "alleles",
        "complete_genotype_observations",
        "genotypes",
    ];
    if locus.len() != FIELDS.len() || locus.keys().any(|key| !FIELDS.contains(&key.as_str())) {
        return Err("locus estadístico inválido");
    }
    if !locus
        .get("contig")
        .and_then(Value::as_str)
        .is_some_and(|value| !value.is_empty())
        || !matches!(
            locus.get("position_1based").and_then(Value::as_u64),
            Some(1..)
        )
        || !locus
            .get("ref")
            .and_then(Value::as_str)
            .is_some_and(is_canonical_base_text)
    {
        return Err("locus estadístico inválido");
    }
    let total = locus
        .get("total_copies")
        .and_then(Value::as_u64)
        .ok_or("total de copias inválido")?;
    let observed = locus
        .get("observed_copies")
        .and_then(Value::as_u64)
        .ok_or("copias observadas inválidas")?;
    let missing = locus
        .get("missing_by_deletion")
        .and_then(Value::as_u64)
        .ok_or("copias ausentes inválidas")?;
    if total
        != individual_count
            .checked_mul(ploidy)
            .ok_or("total de copias inválido")?
        || observed.checked_add(missing) != Some(total)
    {
        return Err("las copias observadas y ausentes no suman el total");
    }
    let alleles = locus
        .get("alleles")
        .and_then(Value::as_array)
        .ok_or("alelos estadísticos inválidos")?;
    if alleles.len() < 2 {
        return Err("alelos estadísticos inválidos");
    }
    let reference = locus.get("ref").and_then(Value::as_str).unwrap_or_default();
    let mut allele_texts = BTreeSet::new();
    let mut count_sum = 0_u64;
    let mut requested_sum = 0.0;
    for (index, allele) in alleles.iter().enumerate() {
        let allele = allele.as_object().ok_or("alelo estadístico inválido")?;
        const REF_FIELDS: [&str; 5] = [
            "allele",
            "kind",
            "requested_frequency",
            "realized_count",
            "realized_frequency",
        ];
        const ALT_FIELDS: [&str; 7] = [
            "allele",
            "kind",
            "variant_id",
            "source_id",
            "requested_frequency",
            "realized_count",
            "realized_frequency",
        ];
        let kind = allele.get("kind").and_then(Value::as_str);
        let fields = match kind {
            Some("ref") if index == 0 => &REF_FIELDS[..],
            Some("alt") if index > 0 => &ALT_FIELDS[..],
            _ => return Err("clase u orden de alelo estadístico inválido"),
        };
        if allele.len() != fields.len()
            || allele.keys().any(|field| !fields.contains(&field.as_str()))
        {
            return Err("alelo estadístico inválido");
        }
        let allele_text = allele
            .get("allele")
            .and_then(Value::as_str)
            .filter(|value| is_canonical_base_text(value))
            .ok_or("alelo estadístico inválido")?;
        if !allele_texts.insert(allele_text) || (index == 0 && allele_text != reference) {
            return Err("alelo estadístico duplicado o referencia incoherente");
        }
        if index > 0
            && (!allele
                .get("variant_id")
                .and_then(Value::as_str)
                .is_some_and(|value| !value.is_empty())
                || allele.get("source_id").and_then(Value::as_str).is_none())
        {
            return Err("procedencia de alelo alternativo inválida");
        }
        let count = allele
            .get("realized_count")
            .and_then(Value::as_u64)
            .ok_or("recuento alélico inválido")?;
        let requested = allele
            .get("requested_frequency")
            .and_then(Value::as_f64)
            .ok_or("frecuencia objetivo inválida")?;
        let realized = allele
            .get("realized_frequency")
            .and_then(Value::as_f64)
            .ok_or("frecuencia realizada inválida")?;
        if !requested.is_finite()
            || !realized.is_finite()
            || !(0.0..=1.0).contains(&requested)
            || !(0.0..=1.0).contains(&realized)
        {
            return Err("frecuencia estadística no finita o fuera de rango");
        }
        let expected_realized = if observed == 0 {
            0.0
        } else {
            count as f64 / observed as f64
        };
        if realized.to_bits() != expected_realized.to_bits() {
            return Err("frecuencia realizada no coincide con recuento y denominador");
        }
        requested_sum += requested;
        count_sum = count_sum
            .checked_add(count)
            .ok_or("recuento alélico desbordado")?;
    }
    if count_sum != observed {
        return Err("los recuentos alélicos no suman las copias observadas");
    }
    if !requested_sum.is_finite() || (requested_sum - 1.0).abs() > 1e-12 {
        return Err("las frecuencias objetivo no suman uno");
    }
    let complete = locus
        .get("complete_genotype_observations")
        .and_then(Value::as_u64)
        .ok_or("observaciones genotípicas inválidas")?;
    if complete > individual_count {
        return Err("observaciones genotípicas fuera de rango");
    }
    let genotypes = locus
        .get("genotypes")
        .and_then(Value::as_array)
        .ok_or("genotipos estadísticos inválidos")?;
    if genotypes.is_empty() {
        return Err("genotipos estadísticos inválidos");
    }
    let mut genotype_keys = BTreeSet::new();
    let mut expected_probability_sum = 0.0;
    let mut realized_genotype_sum = 0_u64;
    for genotype in genotypes {
        let genotype = genotype
            .as_object()
            .ok_or("genotipo estadístico inválido")?;
        const GENOTYPE_FIELDS: [&str; 4] = [
            "allele_dosages",
            "expected_probability",
            "expected_count",
            "realized_count",
        ];
        if genotype.len() != GENOTYPE_FIELDS.len()
            || genotype
                .keys()
                .any(|field| !GENOTYPE_FIELDS.contains(&field.as_str()))
        {
            return Err("genotipo estadístico inválido");
        }
        let dosages = genotype
            .get("allele_dosages")
            .and_then(Value::as_array)
            .ok_or("dosis genotípicas inválidas")?;
        let dosages = dosages
            .iter()
            .map(Value::as_u64)
            .collect::<Option<Vec<_>>>()
            .ok_or("dosis genotípicas inválidas")?;
        let dosage_sum = dosages.iter().try_fold(0_u64, |sum, dosage| {
            sum.checked_add(*dosage)
                .ok_or("dosis genotípicas desbordadas")
        })?;
        if dosages.len() != alleles.len() || dosage_sum != ploidy || !genotype_keys.insert(dosages)
        {
            return Err("dosis genotípicas inválidas o repetidas");
        }
        let probability = genotype
            .get("expected_probability")
            .and_then(Value::as_f64)
            .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
            .ok_or("probabilidad genotípica inválida")?;
        genotype
            .get("expected_count")
            .and_then(Value::as_f64)
            .filter(|value| value.is_finite() && *value >= 0.0)
            .ok_or("recuento genotípico esperado inválido")?;
        // El recuento esperado se conserva como una magnitud de salida. Su
        // correspondencia exacta con la probabilidad y la población se vuelve
        // a calcular desde TOML, semilla y matriz en la puerta física; aquí se
        // valida sólo su representación cerrada y finita para no confundir el
        // contrato de esquema con ese oráculo independiente.
        let realized = genotype
            .get("realized_count")
            .and_then(Value::as_u64)
            .ok_or("recuento genotípico realizado inválido")?;
        expected_probability_sum += probability;
        realized_genotype_sum = realized_genotype_sum
            .checked_add(realized)
            .ok_or("recuento genotípico realizado desbordado")?;
    }
    if !expected_probability_sum.is_finite() || (expected_probability_sum - 1.0).abs() > 1e-12 {
        return Err("probabilidades genotípicas no suman uno");
    }
    if realized_genotype_sum != complete {
        return Err("recuentos genotípicos realizados incoherentes");
    }
    Ok(())
}

fn valid_lower_sha256(value: Option<&Value>) -> bool {
    value.and_then(Value::as_str).is_some_and(|digest| {
        digest.len() == 64
            && digest
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    })
}

fn is_canonical_base_text(value: &str) -> bool {
    matches!(value, "A" | "C" | "G" | "T")
}

fn validate_serialized_normalization_audit(
    value: Option<&Value>,
    detail: &'static str,
) -> Result<(), &'static str> {
    let groups = value.and_then(Value::as_array).ok_or(detail)?;
    let mut previous_group = None;
    for group in groups {
        let group = group.as_object().ok_or(detail)?;
        const GROUP_FIELDS: [&str; 14] = [
            "population_id",
            "contig",
            "position",
            "ref",
            "tolerance",
            "canonical_order",
            "sum_algorithm",
            "original_sum",
            "original_delta",
            "effective_sum",
            "effective_delta",
            "reference_frequency",
            "decision",
            "values",
        ];
        if group.len() != GROUP_FIELDS.len()
            || group
                .keys()
                .any(|field| !GROUP_FIELDS.contains(&field.as_str()))
        {
            return Err(detail);
        }
        let population = group
            .get("population_id")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .ok_or(detail)?;
        let contig = group
            .get("contig")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .ok_or(detail)?;
        let position = group
            .get("position")
            .and_then(Value::as_str)
            .and_then(|value| value.parse::<u64>().ok())
            .filter(|value| *value > 0)
            .ok_or(detail)?;
        let reference = group
            .get("ref")
            .and_then(Value::as_str)
            .and_then(single_canonical_base)
            .ok_or(detail)?;
        if group
            .get("tolerance")
            .and_then(Value::as_f64)
            .filter(|value| value.is_finite() && (0.0..=1e-6).contains(value))
            .is_none()
            || group.get("canonical_order").and_then(Value::as_str)
                != Some("population_id,contig,position,ref,alt,variant_id,source_id")
            || group.get("sum_algorithm").and_then(Value::as_str) != Some("neumaier_sum_v1")
            || [
                "original_sum",
                "original_delta",
                "effective_sum",
                "effective_delta",
            ]
            .iter()
            .any(|field| {
                group
                    .get(*field)
                    .and_then(Value::as_f64)
                    .filter(|value| value.is_finite() && *value >= 0.0)
                    .is_none()
            })
            || group
                .get("reference_frequency")
                .and_then(Value::as_f64)
                .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
                .is_none()
            || !matches!(
                group.get("decision").and_then(Value::as_str),
                Some("reference_remainder" | "normalized")
            )
        {
            return Err(detail);
        }
        let group_key = (population, contig, position, reference);
        if previous_group
            .as_ref()
            .is_some_and(|last| last >= &group_key)
        {
            return Err(detail);
        }
        previous_group = Some(group_key);
        let values = group
            .get("values")
            .and_then(Value::as_array)
            .filter(|values| !values.is_empty())
            .ok_or(detail)?;
        let mut previous_value = None;
        for value in values {
            let value = value.as_object().ok_or(detail)?;
            const VALUE_FIELDS: [&str; 5] =
                ["alt", "variant_id", "source_id", "original", "effective"];
            if value.len() != VALUE_FIELDS.len()
                || value
                    .keys()
                    .any(|field| !VALUE_FIELDS.contains(&field.as_str()))
            {
                return Err(detail);
            }
            let alternate = value
                .get("alt")
                .and_then(Value::as_str)
                .and_then(single_canonical_base)
                .filter(|alternate| *alternate != reference)
                .ok_or(detail)?;
            let variant_id = value
                .get("variant_id")
                .and_then(Value::as_str)
                .filter(|value| !value.is_empty())
                .ok_or(detail)?;
            let source_id = value
                .get("source_id")
                .and_then(Value::as_str)
                .ok_or(detail)?;
            if ["original", "effective"].iter().any(|field| {
                value
                    .get(*field)
                    .and_then(Value::as_f64)
                    .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
                    .is_none()
            }) {
                return Err(detail);
            }
            let value_key = (alternate, variant_id, source_id);
            if previous_value
                .as_ref()
                .is_some_and(|last| last >= &value_key)
            {
                return Err(detail);
            }
            previous_value = Some(value_key);
        }
    }
    Ok(())
}

fn single_canonical_base(value: &str) -> Option<u8> {
    (value.len() == 1)
        .then(|| value.as_bytes()[0])
        .filter(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

fn validate_effective_population_snvs(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let effective_value = effective.get("population_snvs").unwrap_or(&Value::Null);
    let manifest_value = manifest.get("population_snvs").unwrap_or(&Value::Null);
    if effective_value != manifest_value {
        return Err(effective_error(
            effective_path,
            "population_snvs no coincide con el manifiesto",
        ));
    }
    if effective_value.is_null() {
        return Ok(());
    }
    let object = effective_value
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "population_snvs inválido"))?;
    let id = validate_manifest_population_snvs(Some(effective_value), catalog, manifest_path)?;
    if object.get("mode").and_then(Value::as_str) != Some("reuse") && id.is_none() {
        return Err(effective_error(
            effective_path,
            "population_snvs publicado requiere metadata_artifact_id",
        ));
    }
    Ok(())
}

fn validate_effective_authorizations(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let effective_value = effective
        .get("authorizations")
        .ok_or_else(|| effective_error(effective_path, "falta authorizations"))?;
    let manifest_value = manifest
        .get("authorizations")
        .ok_or_else(|| manifest_error(manifest_path, "falta authorizations"))?;
    if effective_value != manifest_value {
        return Err(effective_error(
            effective_path,
            "authorizations no coincide con el manifiesto",
        ));
    }
    validate_authorizations_value(effective_value)
        .map_err(|detail| effective_error(effective_path, detail))
}

fn validate_effective_vcf_output(
    effective: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<Option<VcfOutputMode>, Diagnostic> {
    let Some(value) = effective.get("vcf_output") else {
        // Los documentos F1/F2 previos a F2.5.5 no incluyen la clave. La
        // reconstrucción física contra el TOML actual sí exige `null` explícito.
        return Ok(None);
    };
    if value.is_null() {
        return Ok(None);
    }
    let output = value
        .as_object()
        .ok_or_else(|| effective_error(path, "vcf_output inválido"))?;
    const FIELDS: [&str; 1] = ["mode"];
    validate_effective_closed_object(output, &FIELDS, "vcf_output", path)?;
    match output.get("mode").and_then(Value::as_str) {
        Some("text") => Ok(Some(VcfOutputMode::Text)),
        Some("bgzf") => Ok(Some(VcfOutputMode::Bgzf)),
        _ => Err(effective_error(path, "vcf_output.mode inválido")),
    }
}

fn validate_effective_vcf_against_manifest(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    reference: &serde_json::Map<String, Value>,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let mode = validate_effective_vcf_output(effective, effective_path)?;
    let outputs = required_array(manifest, "vcf_outputs", manifest_path)?;
    let Some(mode) = mode else {
        if !outputs.is_empty() {
            return Err(effective_error(
                effective_path,
                "el manifiesto declara VCF aunque la configuración efectiva lo desactiva",
            ));
        }
        return Ok(());
    };
    if outputs.len() != 1 {
        return Err(effective_error(
            effective_path,
            "la configuración efectiva VCF exige exactamente un enlace en el manifiesto",
        ));
    }
    let output = outputs[0]
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "enlace VCF no es objeto"))?;
    let effective_reference_id = effective_non_empty_text(reference, "id", effective_path)?;
    let expected_mode = match mode {
        VcfOutputMode::Text => "text",
        VcfOutputMode::Bgzf => "bgzf",
    };
    if output.get("reference_id").and_then(Value::as_str) != Some(effective_reference_id)
        || output.get("mode").and_then(Value::as_str) != Some(expected_mode)
        || output.get("vcf_artifact_id").and_then(Value::as_str) != Some(VCF_ARTIFACT_ID)
        || catalog.get(VCF_ARTIFACT_ID) != Some(&VCF_ARTIFACT_KIND)
    {
        return Err(effective_error(
            effective_path,
            "referencia, modo o artefacto VCF del manifiesto no coincide con la configuración efectiva",
        ));
    }
    match mode {
        VcfOutputMode::Text => {
            if !output.get("index_kind").is_some_and(Value::is_null)
                || !output.get("index_artifact_id").is_some_and(Value::is_null)
            {
                return Err(effective_error(
                    effective_path,
                    "VCF de texto efectivo no admite índice en el manifiesto",
                ));
            }
        }
        VcfOutputMode::Bgzf => {
            let expected_index = expected_vcf_index_kind(reference, effective_path)?;
            let (index_kind, artifact_kind) = match expected_index {
                VcfIndexKind::Tbi => ("tbi", "vcf_tbi_index"),
                VcfIndexKind::Csi => ("csi", "vcf_csi_index"),
            };
            if output.get("index_kind").and_then(Value::as_str) != Some(index_kind)
                || output.get("index_artifact_id").and_then(Value::as_str)
                    != Some(VCF_INDEX_ARTIFACT_ID)
                || catalog.get(VCF_INDEX_ARTIFACT_ID) != Some(&artifact_kind)
            {
                return Err(effective_error(
                    effective_path,
                    "índice VCF del manifiesto no coincide con referencia y configuración efectiva",
                ));
            }
        }
    }
    Ok(())
}

fn expected_vcf_index_kind(
    reference: &serde_json::Map<String, Value>,
    effective_path: &Path,
) -> Result<VcfIndexKind, Diagnostic> {
    let sequences = effective_array(reference, "sequences", effective_path)?;
    if sequences.is_empty() {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no tiene secuencias para seleccionar índice VCF",
        ));
    }
    if sequences.iter().any(|sequence| {
        sequence
            .as_object()
            .and_then(|sequence| sequence.get("length"))
            .and_then(Value::as_u64)
            .is_none_or(|length| length >= TBI_MAXIMUM_POSITION_EXCLUSIVE)
    }) {
        Ok(VcfIndexKind::Csi)
    } else {
        Ok(VcfIndexKind::Tbi)
    }
}

fn validate_effective_initial_populations(
    effective: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<BTreeSet<IndividualId>, Diagnostic> {
    let populations = effective_array(effective, "initial_populations", path)?;
    let mut population_ids = BTreeSet::new();
    let mut individual_ids = BTreeSet::new();
    for value in populations {
        let population = value
            .as_object()
            .ok_or_else(|| effective_error(path, "initial_populations[] inválida"))?;
        const FIELDS: [&str; 5] = [
            "population_id",
            "individual_count",
            "individual_ids",
            "ploidy",
            "genotype_model",
        ];
        const DERIVED_FIELDS: [&str; 6] = [
            "population_id",
            "individual_count",
            "individual_ids",
            "ploidy",
            "genotype_model",
            "founder_keys",
        ];
        let fields = if population.contains_key("founder_keys") {
            &DERIVED_FIELDS[..]
        } else {
            &FIELDS[..]
        };
        validate_effective_closed_object(population, fields, "initial_populations[]", path)?;
        let population_id = effective_non_empty_text(population, "population_id", path)?;
        PopulationId::parse(population_id, "effective_config")
            .map_err(|_| effective_error(path, "population_id inicial inválido"))?;
        if !population_ids.insert(population_id) {
            return Err(effective_error(path, "population_id inicial repetido"));
        }
        let count = effective_positive_u64(population, "individual_count", path)?;
        let ploidy = effective_positive_u64(population, "ploidy", path)?;
        if ploidy > 10 {
            return Err(effective_error(path, "ploidy inicial fuera de rango"));
        }
        let ids = effective_array(population, "individual_ids", path)?;
        if usize::try_from(count).ok() != Some(ids.len()) {
            return Err(effective_error(
                path,
                "individual_count inicial no coincide con individual_ids",
            ));
        }
        for id in ids {
            let id = id
                .as_str()
                .filter(|id| !id.is_empty())
                .ok_or_else(|| effective_error(path, "individual_id inicial inválido"))?;
            let individual_id = IndividualId::parse(id, "effective_config")
                .map_err(|_| effective_error(path, "individual_id inicial inválido"))?;
            if !individual_ids.insert(individual_id) {
                return Err(effective_error(path, "individual_id inicial repetido"));
            }
        }
        if let Some(founder_keys_value) = population.get("founder_keys") {
            let mut founder_keys = BTreeSet::new();
            let founder_key_values = founder_keys_value
                .as_array()
                .ok_or_else(|| effective_error(path, "founder_keys inicial inválida"))?;
            if founder_key_values.len() != ids.len() {
                return Err(effective_error(
                    path,
                    "founder_keys inicial no coincide con individual_ids",
                ));
            }
            for key in founder_key_values {
                let key = key
                    .as_str()
                    .ok_or_else(|| effective_error(path, "founder_key inicial inválida"))?;
                FounderKey::parse(key, "effective_config")
                    .map_err(|_| effective_error(path, "founder_key inicial inválida"))?;
                if !founder_keys.insert(key) {
                    return Err(effective_error(path, "founder_key inicial repetida"));
                }
            }
        }
        let model = effective_object(population, "genotype_model", path)?;
        match model.get("kind").and_then(Value::as_str) {
            Some("hwe") if model.len() == 1 => {}
            Some("f_diploid") if model.len() == 2 && ploidy == 2 => {
                let coefficient = model
                    .get("inbreeding_coefficient")
                    .and_then(Value::as_f64)
                    .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
                    .ok_or_else(|| {
                        effective_error(path, "inbreeding_coefficient inicial inválido")
                    })?;
                let _ = coefficient;
            }
            _ => return Err(effective_error(path, "genotype_model inicial inválido")),
        }
    }
    if !populations.is_empty()
        && !effective
            .get("population_snvs")
            .and_then(Value::as_object)
            .is_some_and(|link| link.contains_key("initial_population_statistics"))
    {
        return Err(effective_error(
            path,
            "faltan estadísticas iniciales para cohortes iniciales",
        ));
    }
    Ok(individual_ids)
}

fn validate_effective_direct_events(
    effective: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<BTreeSet<EventId>, Diagnostic> {
    let events = effective_array(effective, "direct_events", path)?;
    let mut event_ids = BTreeSet::new();
    for value in events {
        let event = value
            .as_object()
            .ok_or_else(|| effective_error(path, "direct_events[] inválido"))?;
        let kind = effective_non_empty_text(event, "kind", path)?;
        let fields: &[&str] = match kind {
            "snp" => &[
                "event_id",
                "individual_id",
                "phase",
                "kind",
                "contig",
                "position_1based",
                "ref",
                "alt",
            ],
            "insertion" => &[
                "event_id",
                "individual_id",
                "phase",
                "kind",
                "contig",
                "point",
                "bases",
            ],
            "deletion" => &[
                "event_id",
                "individual_id",
                "phase",
                "kind",
                "contig",
                "start",
                "end",
            ],
            _ => return Err(effective_error(path, "direct_events[].kind inválido")),
        };
        validate_effective_closed_object(event, fields, "direct_events[]", path)?;
        let event_id = effective_non_empty_text(event, "event_id", path)?;
        let event_id = EventId::parse(event_id, "effective_config")
            .map_err(|_| effective_error(path, "event_id directo inválido"))?;
        if !event_ids.insert(event_id) {
            return Err(effective_error(path, "event_id directo repetido"));
        }
        IndividualId::parse(
            effective_non_empty_text(event, "individual_id", path)?,
            "effective_config",
        )
        .map_err(|_| effective_error(path, "individual_id directo inválido"))?;
        SequenceId::parse(
            effective_non_empty_text(event, "contig", path)?,
            "effective_config",
        )
        .map_err(|_| effective_error(path, "contig directo inválido"))?;
        if effective_nonnegative_u64(event, "phase", path)? > 9 {
            return Err(effective_error(path, "fase directa fuera de rango"));
        }
        match kind {
            "snp" => {
                let _ = effective_positive_u64(event, "position_1based", path)?;
                for field in ["ref", "alt"] {
                    let base = effective_non_empty_text(event, field, path)?;
                    if base.len() != 1 || !matches!(base.as_bytes()[0], b'A' | b'C' | b'G' | b'T') {
                        return Err(effective_error(path, "alelo directo inválido"));
                    }
                }
            }
            "insertion" => {
                let _ = effective_nonnegative_u64(event, "point", path)?;
                if effective_non_empty_text(event, "bases", path)?
                    .bytes()
                    .any(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
                {
                    return Err(effective_error(path, "bases de inserción inválidas"));
                }
            }
            "deletion" => {
                let _ = effective_nonnegative_u64(event, "start", path)?;
                let _ = effective_nonnegative_u64(event, "end", path)?;
            }
            _ => unreachable!("kind was checked"),
        }
    }
    Ok(event_ids)
}

fn validate_effective_technical_event(
    effective: &serde_json::Map<String, Value>,
    initial_individual_ids: &BTreeSet<IndividualId>,
    direct_event_ids: &BTreeSet<EventId>,
    path: &Path,
) -> Result<(), Diagnostic> {
    let event = effective_object(effective, "technical_event", path)?;
    const FIELDS: [&str; 11] = [
        "event_id",
        "individual_ids",
        "profile",
        "read_length",
        "quality",
        "artifacts",
        "biological_contamination",
        "libraries",
        "rng_catalog",
        "id_catalog",
        "target",
    ];
    validate_effective_closed_object(event, &FIELDS, "technical_event", path)?;
    let event_id = EventId::parse(
        effective_non_empty_text(event, "event_id", path)?,
        "effective_config",
    )
    .map_err(|_| effective_error(path, "event_id técnico inválido"))?;
    if direct_event_ids.contains(&event_id) {
        return Err(effective_error(
            path,
            "event_id técnico colisiona con un evento directo",
        ));
    }
    let read_length = effective_nonnegative_u64(event, "read_length", path)?;
    let profile = effective_object(event, "profile", path)?;
    if Value::Object(profile.clone()) != f4_technical_profile_object()
        || !(1..=150).contains(&read_length)
    {
        return Err(effective_error(path, "perfil técnico mínimo inválido"));
    }
    if effective_non_empty_text(event, "rng_catalog", path)?
        != crate::randomness::F4_TECHNICAL_RNG_CATALOG
        || event.get("id_catalog") != Some(&f4_id_catalog_value())
    {
        return Err(effective_error(path, "catálogo técnico F4 inválido"));
    }
    let contamination = effective_object(event, "biological_contamination", path)?;
    validate_effective_closed_object(
        contamination,
        &["percent_range", "recipient_ids"],
        "technical_event.biological_contamination",
        path,
    )?;
    let percent_range = effective_object(contamination, "percent_range", path)?;
    validate_effective_closed_object(
        percent_range,
        &["minimum", "maximum"],
        "technical_event.biological_contamination.percent_range",
        path,
    )?;
    let minimum_percent = percent_range.get("minimum").and_then(Value::as_f64);
    let maximum_percent = percent_range.get("maximum").and_then(Value::as_f64);
    if !matches!((minimum_percent, maximum_percent), (Some(minimum), Some(maximum)) if minimum.is_finite() && maximum.is_finite() && (0.0..=100.0).contains(&minimum) && (0.0..=100.0).contains(&maximum) && minimum <= maximum)
    {
        return Err(effective_error(
            path,
            "rango de contaminación biológica inválido",
        ));
    }
    let quality = effective_object(event, "quality", path)?;
    validate_effective_closed_object(
        quality,
        &["acgt_phred", "n_phred"],
        "technical_event.quality",
        path,
    )?;
    if effective_nonnegative_u64(quality, "acgt_phred", path)? != 40
        || effective_nonnegative_u64(quality, "n_phred", path)? != 0
    {
        return Err(effective_error(path, "calidad técnica mínima inválida"));
    }
    let artifacts = effective_object(event, "artifacts", path)?;
    validate_effective_closed_object(
        artifacts,
        &["chimera", "pcr", "optical_duplicates", "gc", "adapters"],
        "technical_event.artifacts",
        path,
    )?;
    let chimera = effective_object(artifacts, "chimera", path)?;
    validate_effective_closed_object(
        chimera,
        &["rate"],
        "technical_event.artifacts.chimera",
        path,
    )?;
    if !matches!(chimera.get("rate").and_then(Value::as_f64), Some(value) if value.is_finite() && (0.0..=1.0).contains(&value))
    {
        return Err(effective_error(
            path,
            "artefacto de quimera efectivo inválido",
        ));
    }
    let pcr = effective_object(artifacts, "pcr", path)?;
    validate_effective_closed_object(pcr, &["applicable"], "technical_event.artifacts.pcr", path)?;
    let optical_duplicates = effective_object(artifacts, "optical_duplicates", path)?;
    validate_effective_closed_object(
        optical_duplicates,
        &["realized_pairs"],
        "technical_event.artifacts.optical_duplicates",
        path,
    )?;
    if pcr.get("applicable") != Some(&Value::Bool(false))
        || effective_nonnegative_u64(optical_duplicates, "realized_pairs", path)? != 0
    {
        return Err(effective_error(
            path,
            "artefactos de duplicación del perfil PCR-free inválidos",
        ));
    }
    let gc = effective_object(artifacts, "gc", path)?;
    validate_effective_closed_object(gc, &["mode"], "technical_event.artifacts.gc", path)?;
    let adapters = effective_object(artifacts, "adapters", path)?;
    validate_effective_closed_object(
        adapters,
        &["enabled"],
        "technical_event.artifacts.adapters",
        path,
    )?;
    if gc.get("mode").and_then(Value::as_str) != Some("neutral")
        || adapters.get("enabled") != Some(&Value::Bool(false))
    {
        return Err(effective_error(
            path,
            "artefactos GC o adaptadores del perfil inicial inválidos",
        ));
    }
    let mut individuals = BTreeSet::new();
    let individual_ids = effective_array(event, "individual_ids", path)?;
    if individual_ids.is_empty() {
        return Err(effective_error(path, "technical_event sin individuos"));
    }
    for individual_id in individual_ids {
        let individual_id = individual_id
            .as_str()
            .ok_or_else(|| effective_error(path, "individual_id técnico inválido"))?;
        let individual_id = IndividualId::parse(individual_id, "effective_config")
            .map_err(|_| effective_error(path, "individual_id técnico inválido"))?;
        if !individuals.insert(individual_id.clone()) {
            return Err(effective_error(path, "individual_id técnico repetido"));
        }
        if !initial_individual_ids.contains(&individual_id) {
            return Err(effective_error(
                path,
                "individual_id técnico fuera de las poblaciones iniciales",
            ));
        }
    }
    let recipient_ids = effective_array(contamination, "recipient_ids", path)?;
    if recipient_ids.is_empty() {
        return Err(effective_error(
            path,
            "contaminación biológica sin receptoras",
        ));
    }
    let mut recipients = BTreeSet::new();
    for recipient_id in recipient_ids {
        let recipient_id = recipient_id.as_str().ok_or_else(|| {
            effective_error(path, "receptora de contaminación biológica inválida")
        })?;
        let recipient_id = IndividualId::parse(recipient_id, "effective_config")
            .map_err(|_| effective_error(path, "receptora de contaminación biológica inválida"))?;
        if !individuals.contains(&recipient_id) || !recipients.insert(recipient_id) {
            return Err(effective_error(
                path,
                "receptoras de contaminación inválidas o repetidas",
            ));
        }
    }
    let libraries = effective_array(event, "libraries", path)?;
    if libraries.len() != individuals.len() {
        return Err(effective_error(
            path,
            "libraries no coincide con individual_ids",
        ));
    }
    let mut library_individuals = BTreeSet::new();
    let mut i7 = BTreeSet::new();
    let mut i5 = BTreeSet::new();
    for library in libraries {
        let library = library
            .as_object()
            .ok_or_else(|| effective_error(path, "library inválida"))?;
        validate_effective_closed_object(
            library,
            &[
                "individual_id",
                "library_id",
                "i7",
                "i5",
                "is_contamination_recipient",
            ],
            "technical_event.libraries",
            path,
        )?;
        let individual_id = IndividualId::parse(
            effective_non_empty_text(library, "individual_id", path)?,
            "effective_config",
        )
        .map_err(|_| effective_error(path, "library individual_id inválido"))?;
        let expected_library = derive_technical_identity_plan(
            RootSeed::parse(
                "root",
                effective_non_empty_text(effective_object(effective, "seed", path)?, "root", path)?,
            )?,
            &ReferenceId::parse(
                effective_non_empty_text(
                    effective_object(effective, "reference", path)?,
                    "id",
                    path,
                )?,
                "effective_config",
            )?,
            event_id.clone(),
            &individuals.iter().cloned().collect::<Vec<_>>(),
        )
        .map_err(|_| effective_error(path, "biblioteca técnica no derivable"))?
        .library(&individual_id)
        .ok_or_else(|| effective_error(path, "library técnica no derivable"))?
        .library_id()
        .as_str()
        .to_owned();
        let valid_udi = |field: &str| {
            effective_non_empty_text(library, field, path).map(|value| {
                value.len() == 8
                    && value
                        .bytes()
                        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            })
        };
        if !individuals.contains(&individual_id)
            || !library_individuals.insert(individual_id.clone())
            || effective_non_empty_text(library, "library_id", path)? != expected_library
            || library.get("is_contamination_recipient")
                != Some(&Value::Bool(recipients.contains(&individual_id)))
            || !valid_udi("i7")?
            || !valid_udi("i5")?
            || !i7.insert(effective_non_empty_text(library, "i7", path)?.to_owned())
            || !i5.insert(effective_non_empty_text(library, "i5", path)?.to_owned())
        {
            return Err(effective_error(path, "biblioteca UDI inválida o repetida"));
        }
    }
    let target = effective_object(event, "target", path)?;
    validate_effective_closed_object(target, &["coverage_wgs"], "technical_event.target", path)?;
    if !matches!(target.get("coverage_wgs").and_then(Value::as_f64), Some(value) if value.is_finite() && value >= 0.0)
    {
        return Err(effective_error(path, "cobertura WGS técnica inválida"));
    }
    Ok(())
}

fn validate_manifest_authorizations(value: Option<&Value>, path: &Path) -> Result<(), Diagnostic> {
    let value = value.ok_or_else(|| manifest_error(path, "falta authorizations"))?;
    validate_authorizations_value(value).map_err(|detail| manifest_error(path, detail))
}

fn validate_authorizations_value(value: &Value) -> Result<(), &'static str> {
    let object = value.as_object().ok_or("authorizations inválido")?;
    if object.len() != 2
        || !matches!(object.get("allow_empty_snv_matrix"), Some(Value::Bool(_)))
        || !matches!(
            object.get("allow_incompatible_genotypes"),
            Some(Value::Bool(_))
        )
    {
        return Err("authorizations inválido");
    }
    Ok(())
}

fn validate_configuration_input(
    object: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<String, Diagnostic> {
    const FIELDS: [&str; 3] = ["kind", "original_path", "artifact_id"];
    validate_closed_object(object, &FIELDS, "configuration_input", path)?;
    if object.get("kind").and_then(Value::as_str) != Some("general") {
        return Err(manifest_error(path, "configuration_input.kind inválido"));
    }
    required_non_empty_text(object, "original_path", path)?;
    let id = required_non_empty_text(object, "artifact_id", path)?;
    require_artifact_kind(catalog, id, "config_general_toml", path)?;
    Ok(id.to_owned())
}

fn validate_references(
    values: &[Value],
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<BTreeMap<String, String>, Diagnostic> {
    const FIELDS: [&str; 2] = ["id", "fasta_artifact_id"];
    let mut references = BTreeMap::new();
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "una referencia no es objeto"))?;
        validate_closed_object(object, &FIELDS, "references[]", path)?;
        let id = required_non_empty_text(object, "id", path)?;
        let fasta_id = required_non_empty_text(object, "fasta_artifact_id", path)?;
        if references
            .insert(id.to_owned(), fasta_id.to_owned())
            .is_some()
        {
            return Err(manifest_error(path, "ID de referencia repetido"));
        }
        require_artifact_kind(catalog, fasta_id, "reference_fasta", path)?;
    }
    Ok(references)
}

struct ValidatedManifestReductions {
    states: Vec<UnitState>,
    derived_reference_ids: BTreeSet<String>,
    requested_artifact_ids: BTreeSet<String>,
}

fn validate_reference_reductions(
    values: &[Value],
    catalog: &BTreeMap<&str, &str>,
    references: &BTreeMap<String, String>,
    path: &Path,
) -> Result<ValidatedManifestReductions, Diagnostic> {
    const FIELDS: [&str; 14] = [
        "id",
        "source_reference_id",
        "eligible_bases",
        "region_size_bases",
        "target",
        "root_seed",
        "selection_subseed",
        "selection_order",
        "state",
        "fasta_artifact_id",
        "map_artifact_id",
        "source_fasta_artifact_id",
        "map_document_kind",
        "map_schema_version",
    ];
    const TARGET_FIELDS: [&str; 3] = ["mode", "declared", "bases"];
    let mut reduction_ids = BTreeSet::new();
    let mut states = Vec::with_capacity(values.len());
    let mut derived_reference_ids = BTreeSet::new();
    let mut requested_artifact_ids = BTreeSet::new();
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "una reducción no es objeto"))?;
        validate_closed_object(object, &FIELDS, "reference_reductions[]", path)?;
        let id = required_non_empty_text(object, "id", path)?;
        let source_id = required_non_empty_text(object, "source_reference_id", path)?;
        if !reduction_ids.insert(id) {
            return Err(manifest_error(path, "ID de reducción repetido"));
        }
        derived_reference_ids.insert(id.to_owned());
        required_positive_u64(object, "eligible_bases", path)?;
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
        validate_lower_hex_64(
            required_non_empty_text(object, "root_seed", path)?,
            "root_seed",
            path,
        )?;
        if required_non_empty_text(object, "map_document_kind", path)?
            != REDUCTION_MAP_DOCUMENT_KIND
            || required_non_empty_text(object, "map_schema_version", path)?
                != REDUCTION_MAP_SCHEMA_VERSION
        {
            return Err(manifest_error(
                path,
                "algoritmo o contrato de mapa de reducción incompatible",
            ));
        }
        validate_lower_hex_64(
            required_non_empty_text(object, "selection_subseed", path)?,
            "selection_subseed",
            path,
        )?;
        let selection_order = required_array(object, "selection_order", path)?;
        if selection_order.is_empty() {
            return Err(manifest_error(path, "selection_order está vacío"));
        }
        for order in selection_order {
            if order.as_u64().is_none_or(|value| value == 0) {
                return Err(manifest_error(path, "selection_order inválido"));
            }
        }
        let state = match object.get("state").and_then(Value::as_str) {
            Some("incomplete") => UnitState::Incomplete,
            Some("complete") => UnitState::Complete,
            Some("failed") => UnitState::Failed,
            _ => return Err(manifest_error(path, "estado de reducción inválido")),
        };
        states.push(state);
        let mut linked = BTreeMap::new();
        for (field, kind) in [
            ("fasta_artifact_id", "reference_fasta"),
            ("map_artifact_id", "reference_reduction_map"),
            ("source_fasta_artifact_id", "reference_fasta"),
        ] {
            let id = optional_text(object, field, path)?;
            if let Some(id) = id {
                require_artifact_kind(catalog, id, kind, path)?;
                if field != "source_fasta_artifact_id" {
                    requested_artifact_ids.insert(id.to_owned());
                }
            }
            linked.insert(field, id);
        }
        if state == UnitState::Complete && linked.values().any(Option::is_none) {
            return Err(manifest_error(
                path,
                "una reducción complete exige sus tres enlaces de artefacto",
            ));
        }
        if linked["fasta_artifact_id"].is_some()
            && linked["fasta_artifact_id"] == linked["source_fasta_artifact_id"]
        {
            return Err(manifest_error(
                path,
                "una reducción no puede reutilizar el mismo FASTA como fuente y derivado",
            ));
        }
        if let Some(fasta_id) = linked["fasta_artifact_id"]
            && references.get(id).map(String::as_str) != Some(fasta_id)
        {
            return Err(manifest_error(
                path,
                "la reducción no coincide con el enlace de su referencia derivada",
            ));
        }
        if let Some(source_fasta_id) = linked["source_fasta_artifact_id"]
            && references.get(source_id).map(String::as_str) != Some(source_fasta_id)
        {
            return Err(manifest_error(
                path,
                "la reducción no coincide con el enlace de su referencia fuente",
            ));
        }
    }
    Ok(ValidatedManifestReductions {
        states,
        derived_reference_ids,
        requested_artifact_ids,
    })
}

fn validate_vcf_outputs(
    values: &[Value],
    catalog: &BTreeMap<&str, &str>,
    references: &BTreeMap<String, String>,
    path: &Path,
) -> Result<BTreeSet<String>, Diagnostic> {
    if values.len() > 1 {
        return Err(manifest_error(
            path,
            "solo se admite un VCF final por ejecución",
        ));
    }
    let mut requested = BTreeSet::new();
    let mut reference_ids = BTreeSet::new();
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "vcf_outputs contiene un valor no objeto"))?;
        const FIELDS: [&str; 7] = [
            "state",
            "reference_id",
            "samples",
            "mode",
            "vcf_artifact_id",
            "index_kind",
            "index_artifact_id",
        ];
        validate_closed_object(object, &FIELDS, "vcf_outputs[]", path)?;
        if object.get("state").and_then(Value::as_str) != Some("complete") {
            return Err(manifest_error(path, "un VCF publicado debe estar complete"));
        }
        let reference_id = required_non_empty_text(object, "reference_id", path)?;
        if !references.contains_key(reference_id) || !reference_ids.insert(reference_id.to_owned())
        {
            return Err(manifest_error(
                path,
                "VCF con referencia ausente o duplicada",
            ));
        }
        let samples = required_array(object, "samples", path)?;
        let mut sample_ids = BTreeSet::new();
        for sample in samples {
            let sample = sample
                .as_str()
                .filter(|sample| !sample.is_empty())
                .ok_or_else(|| manifest_error(path, "muestra VCF inválida"))?;
            if !sample_ids.insert(sample) {
                return Err(manifest_error(path, "muestra VCF duplicada"));
            }
        }
        let vcf_id = required_non_empty_text(object, "vcf_artifact_id", path)?;
        require_artifact_kind(catalog, vcf_id, "truth_vcf", path)?;
        match object.get("mode").and_then(Value::as_str) {
            Some("text") => {
                if !object.get("index_kind").is_some_and(Value::is_null)
                    || !object.get("index_artifact_id").is_some_and(Value::is_null)
                {
                    return Err(manifest_error(path, "VCF de texto no admite índice"));
                }
            }
            Some("bgzf") => {
                let expected_kind = match object.get("index_kind").and_then(Value::as_str) {
                    Some("tbi") => "vcf_tbi_index",
                    Some("csi") => "vcf_csi_index",
                    _ => return Err(manifest_error(path, "VCF BGZF exige índice TBI o CSI")),
                };
                let index_id = optional_text(object, "index_artifact_id", path)?
                    .ok_or_else(|| manifest_error(path, "VCF BGZF exige artefacto de índice"))?;
                require_artifact_kind(catalog, index_id, expected_kind, path)?;
                requested.insert(index_id.to_owned());
            }
            _ => return Err(manifest_error(path, "modo VCF inválido")),
        }
        requested.insert(vcf_id.to_owned());
    }
    Ok(requested)
}

fn validate_artifact(value: &Value, path: &Path) -> Result<ArtifactRecord, Diagnostic> {
    const FIELDS: [&str; 5] = ["id", "kind", "path", "size_bytes", "checksum"];
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
    Ok(ArtifactRecord {
        id: id.to_owned(),
        kind: kind.to_owned(),
        path: relative.to_owned(),
        size_bytes: artifact
            .get("size_bytes")
            .and_then(Value::as_u64)
            .ok_or_else(|| manifest_error(path, "size_bytes inválido"))?,
        digest: digest.to_owned(),
    })
}

fn validate_diagnostics(values: &[Value], path: &Path) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 2] = ["code", "message"];
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "un diagnóstico no es objeto"))?;
        validate_closed_object(object, &FIELDS, "diagnostics[]", path)?;
        let code = required_non_empty_text(object, "code", path)?;
        validate_diagnostic_code(code, path)?;
        let message = required_non_empty_text(object, "message", path)?;
        if !diagnostic_matches_catalog(code, message) {
            return Err(manifest_error(
                path,
                format!("el diagnóstico {code} no coincide con diagnostics-v3.tsv"),
            ));
        }
    }
    Ok(())
}

fn require_artifact_kind(
    catalog: &BTreeMap<&str, &str>,
    id: &str,
    kind: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    match catalog.get(id) {
        Some(observed) if *observed == kind => Ok(()),
        Some(observed) => Err(manifest_error(
            path,
            format!("el enlace {id} exige {kind}, pero el catálogo declara {observed}"),
        )),
        None => Err(manifest_error(path, format!("el enlace {id} no existe"))),
    }
}

fn require_effective_artifact_kind(
    catalog: &BTreeMap<&str, &str>,
    id: &str,
    kind: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    match catalog.get(id) {
        Some(observed) if *observed == kind => Ok(()),
        Some(observed) => Err(effective_error(
            path,
            format!("el enlace {id} exige {kind}, pero el catálogo declara {observed}"),
        )),
        None => Err(effective_error(path, format!("el enlace {id} no existe"))),
    }
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
    let value = object
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| manifest_error(path, format!("falta {field}")))?;
    if value.is_empty() {
        return Err(manifest_error(path, format!("{field} está vacío")));
    }
    Ok(value)
}

fn optional_text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<Option<&'a str>, Diagnostic> {
    match object.get(field) {
        Some(Value::Null) => Ok(None),
        Some(Value::String(value)) if !value.is_empty() => Ok(Some(value)),
        _ => Err(manifest_error(path, format!("{field} inválido"))),
    }
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

fn validate_effective_seed(
    object: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 2] = ["root", "origin"];
    validate_effective_closed_object(object, &FIELDS, "seed", path)?;
    validate_lower_hex_64(
        effective_non_empty_text(object, "root", path)?,
        "seed.root",
        path,
    )?;
    match object.get("origin").and_then(Value::as_str) {
        Some("declared" | "system_csprng") => Ok(()),
        _ => Err(effective_error(path, "seed.origin inválido")),
    }
}

fn validate_effective_normalization(
    object: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 3] = ["tolerance", "canonical_order", "sum_algorithm"];
    validate_effective_closed_object(object, &FIELDS, "normalization", path)?;
    let tolerance = object
        .get("tolerance")
        .and_then(Value::as_f64)
        .filter(|value| value.is_finite() && (0.0..=1e-6).contains(value))
        .ok_or_else(|| effective_error(path, "normalization.tolerance inválida"))?;
    let _ = tolerance;
    if object.get("canonical_order").and_then(Value::as_str)
        != Some("population_id,contig,position,ref,alt,variant_id,source_id")
        || object.get("sum_algorithm").and_then(Value::as_str) != Some("neumaier_sum_v1")
    {
        return Err(effective_error(
            path,
            "metadatos de normalización incompatibles",
        ));
    }
    Ok(())
}

fn validate_effective_reference(
    object: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<(String, String), Diagnostic> {
    const FIELDS: [&str; 8] = [
        "id",
        "source",
        "original_fasta_path",
        "canonical_format",
        "total_genome_size",
        "composition",
        "sequences",
        "fasta_artifact_id",
    ];
    validate_effective_closed_object(object, &FIELDS, "reference", path)?;
    let id = effective_non_empty_text(object, "id", path)?.to_owned();
    let source = match object.get("source").and_then(Value::as_str) {
        Some(source @ ("simulated" | "supplied")) => source,
        _ => return Err(effective_error(path, "reference.source inválido")),
    };
    let declared_total = effective_positive_u64(object, "total_genome_size", path)?;
    let composition_is_object = match object.get("composition") {
        Some(Value::Null) => false,
        Some(Value::Object(composition)) => {
            const BASES: [&str; 4] = ["A", "C", "G", "T"];
            const BASES_WITH_N: [&str; 5] = ["A", "C", "G", "T", "N"];
            if composition
                .keys()
                .any(|key| !BASES_WITH_N.contains(&key.as_str()))
            {
                return Err(effective_error(
                    path,
                    "campos desconocidos en reference.composition",
                ));
            }
            for base in BASES {
                let frequency = composition
                    .get(base)
                    .and_then(Value::as_f64)
                    .ok_or_else(|| {
                        effective_error(path, format!("reference.composition.{base} inválida"))
                    })?;
                if !(0.0..=1.0).contains(&frequency) {
                    return Err(effective_error(
                        path,
                        format!("reference.composition.{base} fuera de rango"),
                    ));
                }
            }
            let sum = BASES
                .iter()
                .map(|base| composition[*base].as_f64().expect("validated frequency"))
                .sum::<f64>();
            if !sum.is_finite() || (sum - 1.0).abs() > COMPOSITION_TOLERANCE {
                return Err(effective_error(
                    path,
                    "reference.composition no suma 1 dentro de la tolerancia",
                ));
            }
            if let Some(value) = composition.get("N") {
                let value = value
                    .as_f64()
                    .ok_or_else(|| effective_error(path, "reference.composition.N inválida"))?;
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    return Err(effective_error(
                        path,
                        "reference.composition.N fuera de rango",
                    ));
                }
            }
            true
        }
        _ => return Err(effective_error(path, "reference.composition inválida")),
    };
    if (source == "simulated") != composition_is_object {
        return Err(effective_error(
            path,
            "reference.composition no corresponde al origen declarado",
        ));
    }
    match (source, object.get("original_fasta_path")) {
        ("simulated", Some(Value::Null)) => {}
        ("supplied", Some(Value::String(value))) if !value.is_empty() => {}
        _ => {
            return Err(effective_error(
                path,
                "reference.original_fasta_path no corresponde al origen declarado",
            ));
        }
    }
    match (source, object.get("canonical_format")) {
        ("simulated", Some(Value::Object(format))) => {
            validate_canonical_fasta(format, path)?;
        }
        ("supplied", Some(Value::Null)) => {}
        _ => {
            return Err(effective_error(
                path,
                "reference.canonical_format no corresponde al origen declarado",
            ));
        }
    }
    let sequences = effective_array(object, "sequences", path)?;
    if sequences.is_empty() {
        return Err(effective_error(path, "reference.sequences está vacío"));
    }
    let mut sequence_ids = BTreeSet::new();
    let mut summed_length = 0_u64;
    for value in sequences {
        let sequence = value
            .as_object()
            .ok_or_else(|| effective_error(path, "una secuencia no es objeto"))?;
        const SEQUENCE_FIELDS: [&str; 7] = [
            "id",
            "class",
            "length",
            "header",
            "description",
            "alphabet_validation",
            "unknown_base_count",
        ];
        validate_effective_closed_object(
            sequence,
            &SEQUENCE_FIELDS,
            "reference.sequences[]",
            path,
        )?;
        let sequence_id = effective_non_empty_text(sequence, "id", path)?;
        if !sequence_ids.insert(sequence_id) {
            return Err(effective_error(path, "ID de secuencia repetido"));
        }
        match sequence.get("class").and_then(Value::as_str) {
            Some("chromosome" | "scaffold") => {}
            _ => return Err(effective_error(path, "clase de secuencia inválida")),
        }
        let length = effective_positive_u64(sequence, "length", path)?;
        summed_length = summed_length
            .checked_add(length)
            .ok_or_else(|| effective_error(path, "la suma de longitudes desborda u64"))?;
        let unknown = sequence
            .get("unknown_base_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| effective_error(path, "unknown_base_count inválido"))?;
        if unknown > length {
            return Err(effective_error(path, "unknown_base_count excede length"));
        }
        validate_reversible_text(
            sequence
                .get("header")
                .ok_or_else(|| effective_error(path, "header ausente"))?,
            "reference.sequences[].header",
            path,
        )?;
        match sequence.get("description") {
            Some(Value::Null) => {}
            Some(value) => {
                validate_reversible_text(value, "reference.sequences[].description", path)?
            }
            None => return Err(effective_error(path, "description ausente")),
        }
        if sequence.get("alphabet_validation").and_then(Value::as_str)
            != Some("ACGTN-case-insensitive")
        {
            return Err(effective_error(path, "alphabet_validation inválida"));
        }
    }
    if summed_length != declared_total {
        return Err(effective_error(
            path,
            "reference.total_genome_size no coincide con la suma de secuencias",
        ));
    }
    let fasta_id = effective_non_empty_text(object, "fasta_artifact_id", path)?.to_owned();
    require_effective_artifact_kind(catalog, &fasta_id, "reference_fasta", path)?;
    Ok((id, fasta_id))
}

fn validate_effective_reduction(
    value: &Value,
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<EffectiveReductionLinks, Diagnostic> {
    let object = value
        .as_object()
        .ok_or_else(|| effective_error(path, "una reducción no es objeto"))?;
    const FIELDS: [&str; 11] = [
        "id",
        "source_reference_id",
        "region_size_bases",
        "target",
        "selection_subseed",
        "selection_order",
        "sequences",
        "state",
        "fasta_artifact_id",
        "map_artifact_id",
        "source_fasta_artifact_id",
    ];
    validate_effective_closed_object(object, &FIELDS, "reduced_references[]", path)?;
    let id = effective_non_empty_text(object, "id", path)?.to_owned();
    let source_reference_id =
        effective_non_empty_text(object, "source_reference_id", path)?.to_owned();
    effective_positive_u64(object, "region_size_bases", path)?;
    let target = effective_object(object, "target", path)?;
    const TARGET_FIELDS: [&str; 3] = ["mode", "declared", "bases"];
    validate_effective_closed_object(target, &TARGET_FIELDS, "reduced_references[].target", path)?;
    match target.get("mode").and_then(Value::as_str) {
        Some("target_bases" | "reduction_percent") => {}
        _ => return Err(effective_error(path, "target.mode inválido")),
    }
    effective_non_empty_text(target, "declared", path)?;
    let target_bases = effective_positive_u64(target, "bases", path)?;
    validate_lower_hex_64(
        effective_non_empty_text(object, "selection_subseed", path)?,
        "selection_subseed",
        path,
    )?;
    let selection_order = effective_array(object, "selection_order", path)?;
    if selection_order.is_empty()
        || selection_order
            .iter()
            .any(|value| value.as_u64().is_none_or(|value| value == 0))
    {
        return Err(effective_error(path, "selection_order inválido"));
    }
    validate_effective_derived_sequences(object, target_bases, path)?;
    if object.get("state").and_then(Value::as_str) != Some("complete") {
        return Err(effective_error(
            path,
            "estado de reducción efectiva inválido",
        ));
    }
    let fasta_artifact_id = effective_non_empty_text(object, "fasta_artifact_id", path)?.to_owned();
    let map_artifact_id = effective_non_empty_text(object, "map_artifact_id", path)?.to_owned();
    let source_fasta_artifact_id =
        effective_non_empty_text(object, "source_fasta_artifact_id", path)?.to_owned();
    require_effective_artifact_kind(catalog, &fasta_artifact_id, "reference_fasta", path)?;
    require_effective_artifact_kind(catalog, &map_artifact_id, "reference_reduction_map", path)?;
    require_effective_artifact_kind(catalog, &source_fasta_artifact_id, "reference_fasta", path)?;
    Ok(EffectiveReductionLinks {
        id,
        source_reference_id,
        fasta_artifact_id,
        map_artifact_id,
        source_fasta_artifact_id,
    })
}

fn validate_effective_derived_sequences(
    reduction: &serde_json::Map<String, Value>,
    target_bases: u64,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 7] = [
        "id",
        "class",
        "length",
        "header",
        "description",
        "alphabet_validation",
        "n_count",
    ];
    let sequences = effective_array(reduction, "sequences", path)?;
    if sequences.is_empty() {
        return Err(effective_error(path, "catálogo derivado vacío"));
    }
    let mut ids = BTreeSet::new();
    let mut total = 0_u64;
    for value in sequences {
        let sequence = value
            .as_object()
            .ok_or_else(|| effective_error(path, "una secuencia derivada no es objeto"))?;
        validate_effective_closed_object(
            sequence,
            &FIELDS,
            "reduced_references[].sequences[]",
            path,
        )?;
        let id = effective_non_empty_text(sequence, "id", path)?;
        if !ids.insert(id) {
            return Err(effective_error(path, "ID de secuencia derivada repetido"));
        }
        match sequence.get("class").and_then(Value::as_str) {
            Some("chromosome" | "scaffold") => {}
            _ => {
                return Err(effective_error(
                    path,
                    "clase de secuencia derivada inválida",
                ));
            }
        }
        let length = effective_positive_u64(sequence, "length", path)?;
        total = total
            .checked_add(length)
            .ok_or_else(|| effective_error(path, "la longitud derivada desborda u64"))?;
        let n_count = sequence
            .get("n_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| effective_error(path, "n_count derivado inválido"))?;
        if n_count > length {
            return Err(effective_error(path, "n_count derivado excede length"));
        }
        validate_reversible_text(
            sequence
                .get("header")
                .ok_or_else(|| effective_error(path, "header derivado ausente"))?,
            "reduced_references[].sequences[].header",
            path,
        )?;
        match sequence.get("description") {
            Some(Value::Null) => {}
            Some(value) => validate_reversible_text(
                value,
                "reduced_references[].sequences[].description",
                path,
            )?,
            None => return Err(effective_error(path, "description derivada ausente")),
        }
        if sequence.get("alphabet_validation").and_then(Value::as_str)
            != Some("ACGTN-case-insensitive")
        {
            return Err(effective_error(
                path,
                "alphabet_validation derivada inválida",
            ));
        }
    }
    if total != target_bases {
        return Err(effective_error(
            path,
            "la suma del catálogo derivado no coincide con target.bases",
        ));
    }
    Ok(())
}

fn validate_effective_closed_object(
    object: &serde_json::Map<String, Value>,
    fields: &[&str],
    context: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if object.len() != fields.len() || object.keys().any(|key| !fields.contains(&key.as_str())) {
        return Err(effective_error(
            path,
            format!("campos ausentes o desconocidos en {context}"),
        ));
    }
    Ok(())
}

fn effective_object<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a serde_json::Map<String, Value>, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_object)
        .ok_or_else(|| effective_error(path, format!("{field} no es objeto")))
}

fn effective_array<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a [Value], Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .ok_or_else(|| effective_error(path, format!("{field} no es lista")))
}

fn effective_non_empty_text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty())
        .ok_or_else(|| effective_error(path, format!("{field} falta o está vacío")))
}

fn effective_positive_u64(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<u64, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_u64)
        .filter(|value| *value > 0)
        .ok_or_else(|| effective_error(path, format!("{field} inválido")))
}

fn effective_nonnegative_u64(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<u64, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_u64)
        .ok_or_else(|| effective_error(path, format!("{field} inválido")))
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
        "original_fasta_path": match config.reference_source() {
            crate::application::EffectiveReferenceSource::Simulated => Value::Null,
            crate::application::EffectiveReferenceSource::Supplied { original_path } => {
                Value::String(original_path.to_string_lossy().into_owned())
            },
        },
        "canonical_format": match config.reference_source() {
            crate::application::EffectiveReferenceSource::Simulated => canonical_fasta_json(),
            crate::application::EffectiveReferenceSource::Supplied { .. } => Value::Null,
        },
        "total_genome_size": config.total_genome_size(),
        "composition": composition,
        "sequences": output.reference.sequences().iter().map(|sequence| json!({
            "id": sequence.id().as_str(),
            "class": sequence_class(sequence.class()),
            "length": sequence.length(),
            "header": reversible_bytes_json(sequence.header().raw_bytes()),
            "description": reversible_description_json(sequence.header().description()),
            "alphabet_validation": alphabet_validation(sequence.alphabet_validation()),
            "unknown_base_count": sequence.unknown_base_count(),
        })).collect::<Vec<_>>(),
        "fasta_artifact_id": output.fasta.id,
    })
}

fn canonical_fasta_json() -> Value {
    json!({
        "version": CANONICAL_FASTA_VERSION,
        "line_base_count": CANONICAL_LINE_BASE_COUNT,
        "line_ending": "LF",
        "base_case": "uppercase",
    })
}

fn validate_canonical_fasta(
    object: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 4] = ["version", "line_base_count", "line_ending", "base_case"];
    validate_effective_closed_object(object, &FIELDS, "reference.canonical_format", path)?;
    if object.get("version").and_then(Value::as_str) != Some(CANONICAL_FASTA_VERSION)
        || object.get("line_base_count").and_then(Value::as_u64)
            != Some(CANONICAL_LINE_BASE_COUNT as u64)
        || object.get("line_ending").and_then(Value::as_str) != Some("LF")
        || object.get("base_case").and_then(Value::as_str) != Some("uppercase")
    {
        return Err(effective_error(path, "reference.canonical_format inválido"));
    }
    Ok(())
}

fn reversible_bytes_json(bytes: &[u8]) -> Value {
    match std::str::from_utf8(bytes) {
        Ok(text) => json!({"encoding": "utf-8", "value": text}),
        Err(_) => json!({"encoding": "bytes-hex", "value": encode_hex(bytes)}),
    }
}

fn reversible_description_json(description: Option<&ReversibleText>) -> Value {
    description.map_or(
        Value::Null,
        |text| json!({"encoding": text.encoding(), "value": text.encoded()}),
    )
}

fn validate_reversible_text(value: &Value, field: &str, path: &Path) -> Result<(), Diagnostic> {
    let object = value
        .as_object()
        .ok_or_else(|| effective_error(path, format!("{field} no es objeto")))?;
    const FIELDS: [&str; 2] = ["encoding", "value"];
    validate_effective_closed_object(object, &FIELDS, field, path)?;
    let encoded = object
        .get("value")
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(path, format!("{field}.value inválido")))?;
    if encoded.is_empty() {
        return Err(effective_error(path, format!("{field}.value está vacío")));
    }
    match object.get("encoding").and_then(Value::as_str) {
        Some("utf-8") => Ok(()),
        Some("bytes-hex")
            if encoded.len() % 2 == 0
                && encoded
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)) =>
        {
            Ok(())
        }
        _ => Err(effective_error(path, format!("{field}.encoding inválido"))),
    }
}

fn alphabet_validation(validation: AlphabetValidation) -> &'static str {
    match validation {
        AlphabetValidation::AcgtnCaseInsensitive => "ACGTN-case-insensitive",
    }
}

fn encode_hex(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(char::from(HEX[usize::from(byte >> 4)]));
        output.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    output
}

fn manifest_reference_json(output: &ManifestReferenceLink) -> Value {
    json!({
        "id": output.id,
        "fasta_artifact_id": output.fasta_artifact_id,
    })
}

fn manifest_reduction_json(output: &ManifestReductionLink) -> Value {
    json!({
        "id": output.id,
        "source_reference_id": output.source_reference_id,
        "eligible_bases": output.eligible_bases,
        "region_size_bases": output.region_size_bases,
        "target": {
            "mode": output.target_mode,
            "declared": output.target_declared,
            "bases": output.target_bases,
        },
        "root_seed": output.root_seed,
        "selection_subseed": output.selection_subseed,
        "selection_order": output.selection_order,
        "state": output.state.as_str(),
        "fasta_artifact_id": output.fasta_artifact_id,
        "map_artifact_id": output.map_artifact_id,
        "source_fasta_artifact_id": output.source_fasta_artifact_id,
        "map_document_kind": output.map_document_kind,
        "map_schema_version": output.map_schema_version,
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
        "sequences": output.reference.sequences().iter().map(|sequence| json!({
            "id": sequence.id().as_str(),
            "class": sequence_class(sequence.class()),
            "length": sequence.length(),
            "header": reversible_bytes_json(sequence.header().raw_bytes()),
            "description": reversible_description_json(sequence.header().description()),
            "alphabet_validation": alphabet_validation(sequence.alphabet_validation()),
            "n_count": sequence.unknown_base_count(),
        })).collect::<Vec<_>>(),
        "state": output.state.as_str(),
        "fasta_artifact_id": output.fasta.id,
        "map_artifact_id": output.map.id,
        "source_fasta_artifact_id": output.source_fasta.id,
    })
}

fn output_categories(manifest: &GeneralManifest<'_>) -> Value {
    let requested = manifest
        .reduced_references
        .iter()
        .flat_map(|reduction| {
            [
                reduction.fasta_artifact_id.as_deref(),
                reduction.map_artifact_id.as_deref(),
            ]
        })
        .flatten()
        .chain(manifest.vcf_outputs.iter().flat_map(|output| {
            std::iter::once(output.vcf_artifact_id.as_str())
                .chain(output.index_artifact_id.as_deref())
        }))
        .collect::<BTreeSet<_>>();
    let minimum = manifest
        .artifacts
        .iter()
        .filter(|artifact| !requested.contains(artifact.id.as_str()))
        .map(|artifact| artifact.id.as_str())
        .collect::<Vec<_>>();
    json!({
        "minimum": minimum,
        "defaulted": [],
        "requested": requested,
        "disabled": [],
    })
}

fn validate_reference_output(
    config: &EffectiveConfig,
    output: ReferenceOutput<'_>,
) -> Result<(), Diagnostic> {
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
    let origin_matches = match (config.reference_source(), output.reference.origin()) {
        (EffectiveReferenceSource::Simulated, ReferenceOrigin::Generated) => true,
        (
            EffectiveReferenceSource::Supplied { original_path },
            ReferenceOrigin::Supplied {
                original_path: observed,
            },
        ) => original_path == observed,
        _ => false,
    };
    if !origin_matches
        || output.reference.id().as_str() != config.reference_id()
        || output.reference.total_length() != config.total_genome_size()
        || output.reference.sequences().len() != config.sequences().len()
        || output
            .reference
            .sequences()
            .iter()
            .zip(config.sequences())
            .any(|(observed, expected)| {
                observed.id().as_str() != expected.id()
                    || observed.class() != expected.class()
                    || observed.length() != expected.length()
            })
    {
        return Err(Diagnostic::internal_invariant(
            "results",
            "reference domain does not match effective configuration",
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

pub(crate) fn validate_relative_link(path: &str) -> Result<(), Diagnostic> {
    let candidate = Path::new(path);
    if path.is_empty()
        || path == "INCOMPLETE"
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
        return Err(manifest_error(
            Path::new("manifest.json"),
            format!("la ruta de artefacto no es canónica: {path:?}"),
        ));
    }
    Ok(())
}

fn effective_schema_version(path: &Path, bytes: &[u8]) -> Result<&'static str, Diagnostic> {
    let document: Value = serde_json::from_slice(bytes)
        .map_err(|error| effective_error(path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| effective_error(path, "la raíz no es un objeto"))?;
    match object.get("schema_version").and_then(Value::as_str) {
        Some(EFFECTIVE_CONFIG_SCHEMA_VERSION) => Ok(EFFECTIVE_CONFIG_SCHEMA_VERSION),
        observed => Err(effective_error_for_schema(
            path,
            "effective-config-v4.schema.json",
            format!(
                "versión de configuración efectiva no admitida: {:?}; se admite exclusivamente \"{}\"",
                observed, EFFECTIVE_CONFIG_SCHEMA_VERSION
            ),
        )),
    }
}

fn retag_format_schema(error: Diagnostic, schema: &'static str) -> Diagnostic {
    if error.code() != "GPS-FMT-006" {
        return error;
    }
    let context = error.context();
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with(
        "path",
        context
            .get("path")
            .cloned()
            .unwrap_or_else(|| "documento JSON".to_owned()),
    )
    .with("schema", schema)
    .with(
        "detail",
        context
            .get("detail")
            .cloned()
            .unwrap_or_else(|| "documento inválido".to_owned()),
    )
}

#[allow(dead_code)]
fn retag_v1_pair_error(
    error: Diagnostic,
    effective_path: &Path,
    manifest_path: &Path,
) -> Diagnostic {
    if error.code() != "GPS-FMT-006" {
        return error;
    }
    let context = error.context();
    let error_path = context.get("path").map(String::as_str);
    let error_schema = context.get("schema").map(String::as_str);
    if error_path == Some(&manifest_path.display().to_string())
        || error_schema.is_some_and(|schema| schema.starts_with("manifest-general-"))
    {
        retag_format_schema(error, "manifest-general-v1.schema.json")
    } else if error_path == Some(&effective_path.display().to_string())
        || error_schema.is_some_and(|schema| schema.starts_with("effective-config-"))
    {
        retag_format_schema(error, "effective-config-v1.schema.json")
    } else {
        error
    }
}

fn manifest_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    manifest_error_for_schema(path, "manifest-general-v4.schema.json", detail)
}

fn manifest_error_for_schema(
    path: &Path,
    schema: &'static str,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path.display().to_string())
    .with("schema", schema)
    .with("detail", detail.into())
}

fn effective_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    effective_error_for_schema(path, "effective-config-v4.schema.json", detail)
}

fn effective_error_for_schema(
    path: &Path,
    schema: &'static str,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path.display().to_string())
    .with("schema", schema)
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
