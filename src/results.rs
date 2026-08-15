//! Documentos JSON versionados de configuración efectiva y manifiestos F1.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Component, Path};

use rustix::fs::{CWD, FileType, Mode, OFlags, fstat, openat};
use serde_json::{Value, json};

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::application::{
    DomainLocus, EffectiveConfig, EffectiveReducedReference, EffectiveValue,
    PhaseGenerationInterval, PopulationMembership, ResolutionDomain, ResolutionRequest,
    build_effective_config_with_entropy_and_domain, validate_resolution_request,
};
use crate::configuration::{DeclaredConfig, ExpectedConfigKind, SequenceClass};
use crate::diagnostics::{Diagnostic, diagnostic_contracts};
use crate::exact_decimal::ExactDecimal;
use crate::identity::SequenceId;
use crate::integrity::{ArtifactIdentity, Sha256Writer, identify_file};
use crate::numeric::{MAX_NORMALIZATION_TOLERANCE, NormalizedCollection, neumaier};
use crate::randomness::{EntropySource, EventSubseed, RootSeed, SubseedRequest, derive_subseed};
use crate::reference::{
    AlphabetValidation, FastaHeader, Reference, ReferenceOrigin,
    SequenceClass as ReferenceSequenceClass,
};
use crate::reference_fasta::{CANONICAL_FASTA_VERSION, CANONICAL_LINE_BASE_COUNT};
use crate::reference_reduction::{
    EligibleSequence, MaterializedPlan, REDUCED_REFERENCE_SELECTION_ALGORITHM,
    ReductionSegmentClass, preflight_reduced_reference,
};
use crate::rules::{
    ParameterKey, PopulationSelector, ResolutionContext, RuleSet, RuleValue, SelectorSet,
    SpatialSelector, TemporalSelector,
};

pub const RESULT_SCHEMA_VERSION: &str = "1";
pub const GENERAL_MANIFEST_RUNTIME_VALIDATION_VERSION: &str =
    "manifest-general-runtime-validation-v1";

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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArtifactRecord {
    pub id: String,
    pub kind: String,
    pub path: String,
    pub size_bytes: u64,
    pub digest: String,
    pub state: UnitState,
    pub validated: bool,
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
            state: UnitState::Complete,
            validated: true,
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
            "state": self.state.as_str(),
            "validated": self.validated,
        })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConfigurationInputRecord {
    pub kind: String,
    pub original_path: String,
    pub copy: ArtifactRecord,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EventManifestLink {
    pub event_id: String,
    pub subseed: String,
    pub state: UnitState,
    pub manifest: ArtifactRecord,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UnitRecord {
    pub id: String,
    pub kind: String,
    pub state: UnitState,
    pub exportable: bool,
    pub diagnostic_code: Option<String>,
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

pub struct GeneralManifest<'a> {
    pub execution_id: &'a str,
    pub state: UnitState,
    pub effective: &'a EffectiveConfig,
    pub configuration_inputs: &'a [ConfigurationInputRecord],
    pub resolution_evidence: &'a ArtifactRecord,
    pub effective_config: Option<&'a ArtifactRecord>,
    pub reference: Option<ReferenceOutput<'a>>,
    pub reduced_references: &'a [ReducedReferenceOutput<'a>],
    pub events: &'a [EventManifestLink],
    pub artifacts: &'a [ArtifactRecord],
    pub units: &'a [UnitRecord],
    pub diagnostics: &'a [Diagnostic],
}

pub struct EventManifest<'a> {
    pub event_id: &'a str,
    pub state: UnitState,
    pub subseed: &'a str,
    pub reference: ReferenceOutput<'a>,
    pub artifacts: &'a [ArtifactRecord],
    pub units: &'a [UnitRecord],
    pub diagnostics: &'a [Diagnostic],
}

fn effective_config_base(config: &EffectiveConfig) -> Value {
    let sequences = config
        .sequences()
        .iter()
        .map(|sequence| {
            json!({
                "id": sequence.id(),
                "class": match sequence.class() {
                    SequenceClass::Chromosome => "chromosome",
                    SequenceClass::Scaffold => "scaffold",
                },
                "length": sequence.length(),
                "ploidy": sequence.ploidy(),
            })
        })
        .collect::<Vec<_>>();
    let event_subseeds = config
        .reproducibility()
        .event_subseeds()
        .iter()
        .map(event_subseed_json)
        .collect::<Vec<_>>();
    let controlled = config
        .controlled_sets()
        .iter()
        .map(|set| {
            json!({
                "id": set.id(),
                "counts": set.counts().iter().map(|count| {
                    let (id, category) = count.identity();
                    json!({"id": id, "category": category, "count": count.count()})
                }).collect::<Vec<_>>()
            })
        })
        .collect::<Vec<_>>();

    json!({
        "effective_config_kind": "gps_effective_config",
        "schema_version": RESULT_SCHEMA_VERSION,
        "source_path": config.source_path(),
        "route": config.route().as_str(),
        "reference": {
            "id": config.reference_id(),
            "total_genome_size": config.total_genome_size(),
            "sequences": sequences,
        },
        "scenario_ids": config.scenario_ids(),
        "event_ids": config.explicit_event_ids(),
        "individual_ids": config.individual_ids(),
        "parameter_order_version": crate::rules::PARAMETER_ORDER_VERSION,
        "resolution_evidence_id": "resolution-evidence",
        "resolution_domain": resolution_domain_json(config.resolution_domain()),
        "randomness": {
            "root_seed": config.reproducibility().root_seed().to_hex(),
            "seed_origin": config.reproducibility().seed_origin().as_str(),
            "algorithm": config.reproducibility().random_algorithm(),
            "derivation_version": config.reproducibility().derivation_version(),
            "catalog_version": config.reproducibility().catalog_version(),
            "event_subseeds": event_subseeds,
        },
        "normalizations": config.normalized_collections().iter().map(normalization_json).collect::<Vec<_>>(),
        "controlled_frequency_sets": controlled,
        "effective_values": config.values().iter().map(effective_value_json).collect::<Vec<_>>(),
    })
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
    validate_reference_output(output)?;
    for reduction in reductions {
        validate_reference_output(ReferenceOutput {
            reference: reduction.reference,
            fasta: reduction.fasta,
        })?;
    }
    let mut document = effective_config_base(config);
    document["reference"] = reference_json(config, output);
    document["reduced_references"] = Value::Array(
        reductions
            .iter()
            .map(|reduction| reduction_json(config, *reduction))
            .collect(),
    );
    Ok(document)
}

#[must_use]
pub fn general_manifest_document(manifest: &GeneralManifest<'_>) -> Value {
    let dependencies = [
        ("blake3", "1.8.6"),
        ("getrandom", "0.4.3"),
        ("jsonschema", "0.49.9"),
        ("noodles-fasta", "0.66.0"),
        ("rand_chacha", "0.10.0"),
        ("rustix", "1.1.4"),
        ("serde_json", "1.0.149"),
        ("sha2", "0.10.9"),
        ("toml", "1.1.4"),
    ]
    .into_iter()
    .map(|(name, version)| json!({"name": name, "version": version}))
    .collect::<Vec<_>>();
    let configuration_inputs = manifest
        .configuration_inputs
        .iter()
        .map(|input| {
            json!({
                "kind": input.kind,
                "original_path": input.original_path,
                "copy": input.copy.to_json(),
            })
        })
        .collect::<Vec<_>>();
    let events = manifest
        .events
        .iter()
        .map(|event| {
            json!({
                "event_id": event.event_id,
                "subseed": event.subseed,
                "state": event.state.as_str(),
                "manifest": event.manifest.to_json(),
            })
        })
        .collect::<Vec<_>>();

    json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": RESULT_SCHEMA_VERSION,
        "execution_id": manifest.execution_id,
        "state": manifest.state.as_str(),
        "route": manifest.effective.route().as_str(),
        "program": {"name": env!("CARGO_PKG_NAME"), "version": env!("CARGO_PKG_VERSION")},
        "platform": {"os": std::env::consts::OS, "architecture": std::env::consts::ARCH, "official_target": "x86_64-unknown-linux-gnu"},
        "toolchain": {"rustc": "1.97.1", "cargo": "1.97.1", "edition": "2024", "lock": "Cargo.lock"},
        "dependencies": dependencies,
        "tools": [],
        "profiles": [],
        "checks": [
            {"name": "configuration_validation", "status": "passed"},
            {"name": "artifact_size_and_sha256", "status": "passed"},
            {"name": "result_json_internal_contract", "status": "passed"}
        ],
        "configuration_inputs": configuration_inputs,
        "resolution_evidence": manifest.resolution_evidence.to_json(),
        "effective_config": manifest.effective_config.map(ArtifactRecord::to_json),
        "references": manifest.reference.map(|reference| reference_json(manifest.effective, reference)).into_iter().chain(manifest.reduced_references.iter().filter(|reduction| reduction.state == UnitState::Complete).map(|reduction| reference_json(manifest.effective, ReferenceOutput { reference: reduction.reference, fasta: reduction.fasta }))).collect::<Vec<_>>(),
        "reference_reductions": manifest.reduced_references.iter().map(|reduction| reduction_json(manifest.effective, *reduction)).collect::<Vec<_>>(),
        "events": events,
        "artifacts": manifest.artifacts.iter().map(ArtifactRecord::to_json).collect::<Vec<_>>(),
        "outputs": general_output_categories(manifest),
        "normalizations": manifest.effective.normalized_collections().iter().map(normalization_json).collect::<Vec<_>>(),
        "reproducibility": reproducibility_levels(),
        "units": manifest.units.iter().map(unit_json).collect::<Vec<_>>(),
        "diagnostics": manifest.diagnostics.iter().map(diagnostic_json).collect::<Vec<_>>(),
    })
}

fn general_output_categories(manifest: &GeneralManifest<'_>) -> Value {
    let reduction_artifact_ids = manifest
        .effective
        .reduced_references()
        .iter()
        .flat_map(|reduction| {
            [
                format!("reference-fasta-{}", reduction.id()),
                format!("reference-map-{}", reduction.id()),
            ]
        })
        .collect::<BTreeSet<_>>();
    let mut minimum = vec!["manifest-general".to_owned()];
    let mut explicitly_requested = Vec::new();
    for artifact in manifest.artifacts {
        if reduction_artifact_ids.contains(&artifact.id) {
            explicitly_requested.push(artifact.id.clone());
        } else {
            minimum.push(artifact.id.clone());
        }
    }
    minimum.sort();
    explicitly_requested.sort();
    json!({
        "minimum": minimum,
        "default_enabled": [],
        "explicitly_requested": explicitly_requested,
        "disabled": ["truth", "vcf", "fastq", "bam", "cram"],
    })
}

#[must_use]
pub fn event_manifest_document(manifest: &EventManifest<'_>) -> Value {
    json!({
        "manifest_kind": "gps_event_manifest",
        "schema_version": RESULT_SCHEMA_VERSION,
        "event_id": manifest.event_id,
        "state": manifest.state.as_str(),
        "subseed": manifest.subseed,
        "references": [reference_link_json(manifest.reference)],
        "artifacts": manifest.artifacts.iter().map(ArtifactRecord::to_json).collect::<Vec<_>>(),
        "outputs": {
            "minimum": ["manifest_event"],
            "default_enabled": [],
            "explicitly_requested": [],
            "disabled": ["truth", "vcf", "fastq", "bam", "cram"]
        },
        "units": manifest.units.iter().map(unit_json).collect::<Vec<_>>(),
        "diagnostics": manifest.diagnostics.iter().map(diagnostic_json).collect::<Vec<_>>(),
    })
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

/// Capacidad opaca que acredita las puertas estructural y semántica del
/// manifiesto general v1 sobre una secuencia concreta de bytes.
pub struct ValidatedGeneralManifest {
    bytes: Vec<u8>,
    artifacts: Vec<ArtifactRecord>,
    workspace_token: u64,
    workspace_snapshot: Vec<ArtifactRecord>,
}

impl std::fmt::Debug for ValidatedGeneralManifest {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ValidatedGeneralManifest")
            .field("size_bytes", &self.bytes.len())
            .field("artifact_count", &self.artifacts.len())
            .field("workspace_bound", &true)
            .finish_non_exhaustive()
    }
}

impl ValidatedGeneralManifest {
    pub(crate) fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub(crate) fn artifacts(&self) -> &[ArtifactRecord] {
        &self.artifacts
    }

    pub(crate) const fn workspace_token(&self) -> u64 {
        self.workspace_token
    }

    pub(crate) fn workspace_snapshot(&self) -> &[ArtifactRecord] {
        &self.workspace_snapshot
    }
}

/// Ejecuta la puerta estructural v1, las igualdades dinámicas y la relectura
/// semántica de la configuración efectiva y los manifiestos de evento enlazados.
/// `path` debe identificar la ubicación real del manifiesto dentro del paquete.
///
/// Esta puerta pública no concede capacidad de escritura o publicación. Solo un
/// [`crate::publication::ExecutionWorkspace`] puede ligar unos bytes validados a
/// su instancia y catálogo privados.
pub fn validate_general_manifest_bytes(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    parse_general_manifest_bytes(path, bytes).map(|_| ())
}

pub(crate) fn bind_general_manifest_to_workspace(
    path: &Path,
    bytes: &[u8],
    workspace_token: u64,
    workspace_snapshot: Vec<ArtifactRecord>,
) -> Result<ValidatedGeneralManifest, Diagnostic> {
    let (bytes, artifacts) = parse_general_manifest_bytes(path, bytes)?;
    Ok(ValidatedGeneralManifest {
        bytes,
        artifacts,
        workspace_token,
        workspace_snapshot,
    })
}

pub(crate) fn parse_general_manifest_bytes(
    path: &Path,
    bytes: &[u8],
) -> Result<(Vec<u8>, Vec<ArtifactRecord>), Diagnostic> {
    let document: Value = serde_json::from_slice(bytes)
        .map_err(|error| semantic_manifest_error(path, format!("JSON no legible: {error}")))?;
    crate::manifest_schema::validate_general_manifest_structure(&document)
        .map_err(|detail| semantic_manifest_error(path, format!("estructura v1: {detail}")))?;
    validate_general_manifest_graph(&document)
        .map_err(|detail| semantic_manifest_error(path, format!("semántica v1: {detail}")))?;
    validate_linked_manifest_documents(path, &document).map_err(|detail| {
        semantic_manifest_error(path, format!("documentos enlazados v1: {detail}"))
    })?;
    let artifacts = manifest_artifact_records(&document)
        .map_err(|detail| semantic_manifest_error(path, format!("catálogo v1: {detail}")))?;
    Ok((bytes.to_vec(), artifacts))
}

fn manifest_artifact_records(document: &Value) -> Result<Vec<ArtifactRecord>, String> {
    let manifest = object(document, "manifest")?;
    array(manifest, "artifacts", "manifest")?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let context = format!("artifacts[{index}]");
            let artifact = object(value, &context)?;
            let checksum = artifact
                .get("checksum")
                .ok_or_else(|| format!("falta {context}.checksum"))?;
            let checksum = object(checksum, &format!("{context}.checksum"))?;
            if text(checksum, "algorithm", &format!("{context}.checksum"))? != "sha256" {
                return Err(format!("{context} no usa SHA-256"));
            }
            let state = match text(artifact, "state", &context)? {
                "complete" => UnitState::Complete,
                "failed" => UnitState::Failed,
                other => return Err(format!("estado de artefacto no admitido: {other}")),
            };
            let size_bytes = artifact
                .get("size_bytes")
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("falta o no es u64 {context}.size_bytes"))?;
            Ok(ArtifactRecord {
                id: text(artifact, "id", &context)?.to_owned(),
                kind: text(artifact, "kind", &context)?.to_owned(),
                path: text(artifact, "path", &context)?.to_owned(),
                size_bytes,
                digest: text(checksum, "digest", &format!("{context}.checksum"))?.to_owned(),
                state,
                validated: boolean(artifact, "validated", &context)?,
            })
        })
        .collect()
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ManifestSequenceCatalog {
    class: String,
    length: u64,
    ploidy: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ManifestReferenceCatalog {
    total_genome_size: u64,
    order: Vec<String>,
    sequences: BTreeMap<String, ManifestSequenceCatalog>,
}

fn manifest_reference_catalog(
    reference: &Value,
    context: &str,
) -> Result<ManifestReferenceCatalog, String> {
    let reference = object(reference, context)?;
    let sequence_values = array(reference, "sequences", context)?;
    let mut sequences = BTreeMap::new();
    let mut order = Vec::with_capacity(sequence_values.len());
    let mut calculated_total = 0_u64;
    for (index, sequence) in sequence_values.iter().enumerate() {
        let sequence_context = format!("{context}.sequences[{index}]");
        let sequence = object(sequence, &sequence_context)?;
        let id = text(sequence, "id", &sequence_context)?;
        validate_manifest_sequence_header(sequence, id, &sequence_context)?;
        let length = unsigned(sequence, "length", &sequence_context)?;
        let unknown_base_count = unsigned(sequence, "unknown_base_count", &sequence_context)?;
        if unknown_base_count > length {
            return Err(format!(
                "{sequence_context}.unknown_base_count excede la longitud"
            ));
        }
        calculated_total = calculated_total
            .checked_add(length)
            .ok_or_else(|| format!("{context}.total_genome_size desborda u64"))?;
        let summary = ManifestSequenceCatalog {
            class: text(sequence, "class", &sequence_context)?.to_owned(),
            length,
            ploidy: unsigned(sequence, "ploidy", &sequence_context)?,
        };
        if sequences.insert(id.to_owned(), summary).is_some() {
            return Err(format!("ID de secuencia duplicado en {context}: {id}"));
        }
        order.push(id.to_owned());
    }
    let declared_total = unsigned(reference, "total_genome_size", context)?;
    if declared_total != calculated_total {
        return Err(format!(
            "{context}.total_genome_size no suma sus secuencias: {declared_total}/{calculated_total}"
        ));
    }
    Ok(ManifestReferenceCatalog {
        total_genome_size: declared_total,
        order,
        sequences,
    })
}

fn validate_manifest_sequence_header(
    sequence: &serde_json::Map<String, Value>,
    sequence_id: &str,
    context: &str,
) -> Result<(), String> {
    let header_value = sequence
        .get("header")
        .ok_or_else(|| format!("falta {context}.header"))?;
    let header_bytes = reversible_text_bytes(header_value, &format!("{context}.header"))?;
    let header = FastaHeader::parse(&header_bytes)
        .map_err(|error| format!("{context}.header no es FASTA reversible: {error:?}"))?;
    if header.id_text() != sequence_id {
        return Err(format!(
            "{context}.header no comienza por el ID de secuencia"
        ));
    }
    let declared_description = sequence
        .get("description")
        .ok_or_else(|| format!("falta {context}.description"))?;
    match (header.description(), declared_description.is_null()) {
        (None, true) => Ok(()),
        (Some(description), false) => {
            let object = object(declared_description, &format!("{context}.description"))?;
            let encoding = text(object, "encoding", &format!("{context}.description"))?;
            let value = text(object, "value", &format!("{context}.description"))?;
            if encoding == description.encoding() && value == description.encoded() {
                Ok(())
            } else {
                Err(format!(
                    "{context}.description no reconstruye la descripción de header"
                ))
            }
        }
        _ => Err(format!(
            "{context}.description no coincide con la presencia de descripción en header"
        )),
    }
}

fn reversible_text_bytes(value: &Value, context: &str) -> Result<Vec<u8>, String> {
    let value = object(value, context)?;
    let encoding = text(value, "encoding", context)?;
    let encoded = text(value, "value", context)?;
    match encoding {
        "utf-8" => Ok(encoded.as_bytes().to_vec()),
        "bytes-hex" => {
            let bytes = encoded.as_bytes();
            if bytes.len() % 2 != 0 {
                return Err(format!("{context}.value no contiene pares hexadecimales"));
            }
            bytes
                .chunks_exact(2)
                .map(|pair| {
                    let high = hex_nibble(pair[0])?;
                    let low = hex_nibble(pair[1])?;
                    Ok(high * 16 + low)
                })
                .collect()
        }
        other => Err(format!("{context}.encoding no admitido: {other}")),
    }
}

fn hex_nibble(byte: u8) -> Result<u8, String> {
    match byte {
        b'0'..=b'9' => Ok(byte - b'0'),
        b'a'..=b'f' => Ok(byte - b'a' + 10),
        _ => Err("dígito hexadecimal reversible no canónico".to_owned()),
    }
}

fn validate_manifest_diagnostic_graph(
    manifest: &serde_json::Map<String, Value>,
    units: &[Value],
) -> Result<(), String> {
    let contracts = diagnostic_contracts()
        .map_err(|detail| format!("catálogo diagnostics-v2 inválido: {detail}"))?;
    let mut embedded = BTreeSet::<String>::new();
    for (index, diagnostic) in array(manifest, "diagnostics", "manifest")?
        .iter()
        .enumerate()
    {
        let context = format!("diagnostics[{index}]");
        let diagnostic = object(diagnostic, &context)?;
        let code = text(diagnostic, "code", &context)?;
        let contract = contracts
            .get(code)
            .ok_or_else(|| format!("{context}.code no existe en diagnostics-v2: {code}"))?;
        if text(diagnostic, "category", &context)? != contract.category {
            return Err(format!("{context}.category contradice diagnostics-v2"));
        }
        if text(diagnostic, "message", &context)? != contract.meaning {
            return Err(format!("{context}.message contradice diagnostics-v2"));
        }
        let diagnostic_context = diagnostic
            .get("context")
            .ok_or_else(|| format!("falta {context}.context"))?;
        let diagnostic_context = object(diagnostic_context, &format!("{context}.context"))?;
        for required in &contract.required_context {
            if diagnostic_context
                .get(*required)
                .and_then(Value::as_str)
                .is_none_or(str::is_empty)
            {
                return Err(format!(
                    "{context}.context no aporta el campo requerido {required}"
                ));
            }
        }
        if !embedded.insert(code.to_owned()) {
            return Err(format!("código diagnóstico duplicado: {code}"));
        }
    }

    let mut referenced = BTreeSet::new();
    for (index, unit) in units.iter().enumerate() {
        let context = format!("units[{index}]");
        let unit = object(unit, &context)?;
        let state = text(unit, "state", &context)?;
        let diagnostic_code = unit
            .get("diagnostic_code")
            .ok_or_else(|| format!("falta {context}.diagnostic_code"))?;
        match state {
            "complete" if diagnostic_code.is_null() => {}
            "failed" => {
                let code = diagnostic_code
                    .as_str()
                    .ok_or_else(|| format!("{context} fallida sin diagnostic_code"))?;
                if !embedded.contains(code) {
                    return Err(format!(
                        "{context}.diagnostic_code no resuelve un diagnóstico embebido: {code}"
                    ));
                }
                referenced.insert(code.to_owned());
            }
            _ => {
                return Err(format!(
                    "{context} tiene estado y diagnostic_code incompatibles"
                ));
            }
        }
    }
    for code in &embedded {
        if !referenced.contains(code) {
            return Err(format!("diagnóstico embebido huérfano: {code}"));
        }
    }
    Ok(())
}

fn validate_reduction_metadata(
    reduction: &serde_json::Map<String, Value>,
    source: &ManifestReferenceCatalog,
    derived: Option<&ManifestReferenceCatalog>,
    context: &str,
) -> Result<(), String> {
    let eligible_bases = unsigned(reduction, "eligible_bases", context)?;
    if eligible_bases != source.total_genome_size {
        return Err(format!(
            "{context}.eligible_bases contradice el catálogo fuente"
        ));
    }
    let region_size = unsigned(reduction, "region_size_bases", context)?;
    let target = reduction
        .get("target")
        .ok_or_else(|| format!("falta {context}.target"))?;
    let target = object(target, &format!("{context}.target"))?;
    let derived_bases = unsigned(target, "derived_bases", &format!("{context}.target"))?;
    if derived_bases == 0 || derived_bases >= eligible_bases {
        return Err(format!("{context}.target queda fuera del espacio elegible"));
    }
    let declared_value = text(target, "declared_value", &format!("{context}.target"))?;
    match text(target, "mode", &format!("{context}.target"))? {
        "target_bases" => {
            let declared = declared_value
                .replace('_', "")
                .parse::<u64>()
                .map_err(|_| format!("{context}.target.declared_value no es un u64"))?;
            if declared != derived_bases {
                return Err(format!(
                    "{context}.target.declared_value contradice derived_bases"
                ));
            }
        }
        "reduction_percent" => {
            let retained = ExactDecimal::parse("reduction_percent", declared_value)
                .and_then(|value| value.retained_after_percent(eligible_bases))
                .map_err(|diagnostic| diagnostic.render_human())?;
            if retained != derived_bases {
                return Err(format!(
                    "{context}.target.declared_value contradice el cálculo decimal exacto"
                ));
            }
        }
        other => return Err(format!("{context}.target.mode no admitido: {other}")),
    }
    let eligible_sequences = source
        .sequences
        .iter()
        .map(|(id, sequence)| {
            let id =
                SequenceId::parse(id, context).map_err(|diagnostic| diagnostic.render_human())?;
            let class = match sequence.class.as_str() {
                "chromosome" => ReferenceSequenceClass::Chromosome,
                "scaffold" => ReferenceSequenceClass::Scaffold,
                other => return Err(format!("{context} contiene clase de secuencia {other}")),
            };
            EligibleSequence::new(id, class, sequence.length)
                .map_err(|error| format!("{context} catálogo no elegible: {error:?}"))
        })
        .collect::<Result<Vec<_>, String>>()?;
    preflight_reduced_reference(&eligible_sequences, derived_bases, region_size)
        .map_err(|error| format!("{context} no supera preflight canónico: {error:?}"))?;

    let included = manifest_string_set(
        reduction
            .get("included_sequence_ids")
            .ok_or_else(|| format!("falta {context}.included_sequence_ids"))?,
        &format!("{context}.included_sequence_ids"),
    )?;
    let omitted = manifest_string_set(
        reduction
            .get("omitted_sequence_ids")
            .ok_or_else(|| format!("falta {context}.omitted_sequence_ids"))?,
        &format!("{context}.omitted_sequence_ids"),
    )?;
    let source_ids = source.sequences.keys().cloned().collect::<BTreeSet<_>>();
    if !included.is_disjoint(&omitted)
        || included.union(&omitted).cloned().collect::<BTreeSet<_>>() != source_ids
    {
        return Err(format!(
            "{context}.included_sequence_ids/omitted_sequence_ids no forman la partición fuente"
        ));
    }

    let selection = array(reduction, "selection_order", context)?;
    let full_count = derived_bases / region_size;
    let residual = derived_bases % region_size;
    let expected_count = full_count
        .checked_add(u64::from(residual > 0))
        .ok_or_else(|| format!("{context}.selection_order desborda u64"))?;
    if u64::try_from(selection.len()).ok() != Some(expected_count)
        || boolean(reduction, "has_residual", context)? != (residual > 0)
    {
        return Err(format!(
            "{context}.selection_order/has_residual contradice región y objetivo"
        ));
    }
    let mut selected_ids = BTreeSet::new();
    let mut intervals = BTreeMap::<String, Vec<(u64, u64)>>::new();
    let mut selected_lengths = BTreeMap::<String, u64>::new();
    let mut selected_total = 0_u64;
    for (index, entry) in selection.iter().enumerate() {
        let entry_context = format!("{context}.selection_order[{index}]");
        let entry = object(entry, &entry_context)?;
        let order = unsigned(entry, "selection_order", &entry_context)?;
        let expected_order = u64::try_from(index)
            .ok()
            .and_then(|value| value.checked_add(1))
            .ok_or_else(|| format!("{entry_context}.selection_order desborda u64"))?;
        if order != expected_order {
            return Err(format!("{entry_context}.selection_order no es consecutivo"));
        }
        let sequence_id = text(entry, "sequence_id", &entry_context)?;
        let source_sequence = source.sequences.get(sequence_id).ok_or_else(|| {
            format!("{entry_context}.sequence_id no existe en la referencia fuente")
        })?;
        let start = unsigned(entry, "start", &entry_context)?;
        let end = unsigned(entry, "end", &entry_context)?;
        if start == 0 || start > end || end > source_sequence.length {
            return Err(format!(
                "{entry_context} queda fuera de la secuencia fuente"
            ));
        }
        let length = end
            .checked_sub(start)
            .and_then(|value| value.checked_add(1))
            .ok_or_else(|| format!("{entry_context} tiene longitud inválida"))?;
        let expected_length = if order <= full_count {
            region_size
        } else if residual > 0 && order == expected_count {
            residual
        } else {
            return Err(format!("{entry_context} no corresponde a región/residual"));
        };
        if length != expected_length {
            return Err(format!("{entry_context} contradice región/residual"));
        }
        selected_total = selected_total
            .checked_add(length)
            .ok_or_else(|| format!("{context}.selection_order desborda el total"))?;
        let selected_length = selected_lengths.entry(sequence_id.to_owned()).or_default();
        *selected_length = selected_length
            .checked_add(length)
            .ok_or_else(|| format!("{entry_context} desborda la suma por secuencia"))?;
        intervals
            .entry(sequence_id.to_owned())
            .or_default()
            .push((start, end));
        selected_ids.insert(sequence_id.to_owned());
    }
    if selected_total != derived_bases || selected_ids != included {
        return Err(format!(
            "{context}.selection_order no suma el objetivo o no coincide con included"
        ));
    }
    for (sequence_id, values) in &mut intervals {
        values.sort_unstable();
        if values.windows(2).any(|pair| pair[1].0 <= pair[0].1) {
            return Err(format!(
                "{context}.selection_order solapa intervalos de {sequence_id}"
            ));
        }
    }

    if let Some(derived) = derived {
        let expected_derived_order = source
            .order
            .iter()
            .filter(|id| included.contains(*id))
            .cloned()
            .collect::<Vec<_>>();
        if derived.total_genome_size != derived_bases
            || derived.sequences.keys().cloned().collect::<BTreeSet<_>>() != included
            || derived.order != expected_derived_order
        {
            return Err(format!(
                "{context} no coincide con el catálogo u orden de referencia derivada"
            ));
        }
        for (id, derived_sequence) in &derived.sequences {
            let source_sequence = &source.sequences[id];
            if selected_lengths.get(id) != Some(&derived_sequence.length)
                || derived_sequence.class != source_sequence.class
                || derived_sequence.ploidy != source_sequence.ploidy
            {
                return Err(format!(
                    "{context} contradice longitud/clase/ploidía derivada de {id}"
                ));
            }
        }
    }
    Ok(())
}

fn manifest_string_set(value: &Value, context: &str) -> Result<BTreeSet<String>, String> {
    let values = value
        .as_array()
        .ok_or_else(|| format!("{context} no es un array"))?;
    let mut result = BTreeSet::new();
    for (index, value) in values.iter().enumerate() {
        let value = value
            .as_str()
            .ok_or_else(|| format!("{context}[{index}] no es texto"))?;
        if !result.insert(value.to_owned()) {
            return Err(format!("{context} repite {value}"));
        }
    }
    Ok(result)
}

fn validate_general_manifest_graph(document: &Value) -> Result<(), String> {
    let manifest = object(document, "manifest")?;
    if text(manifest, "manifest_kind", "manifest")? != "gps_execution_manifest" {
        return Err("manifest_kind no identifica un manifiesto general".to_owned());
    }
    let manifest_state = text(manifest, "state", "manifest")?;
    if !matches!(manifest_state, "complete" | "failed") {
        return Err(format!(
            "estado de manifiesto no admitido: {manifest_state}"
        ));
    }
    validate_normalizations(
        manifest
            .get("normalizations")
            .ok_or_else(|| "falta manifest.normalizations".to_owned())?,
        "manifest.normalizations",
    )?;

    let artifact_values = array(manifest, "artifacts", "manifest")?;
    let mut artifacts = BTreeMap::<String, &Value>::new();
    let mut artifact_paths = BTreeMap::<String, String>::new();
    for (index, artifact) in artifact_values.iter().enumerate() {
        let context = format!("artifacts[{index}]");
        let (id, _, _, _) = artifact_fields(artifact, &context)?;
        if artifacts.insert(id.to_owned(), artifact).is_some() {
            return Err(format!("ID de artefacto duplicado en catálogo: {id}"));
        }
        let path = text(object(artifact, &context)?, "path", &context)?;
        if let Some(previous) = artifact_paths.insert(path.to_owned(), id.to_owned()) {
            return Err(format!(
                "ruta de artefacto reutilizada por {previous} y {id}: {path}"
            ));
        }
    }
    if artifacts.contains_key("manifest-general") {
        return Err("manifest-general es una salida raíz reservada, no un artefacto".to_owned());
    }

    let configuration_values = array(manifest, "configuration_inputs", "manifest")?;
    let mut general_configuration_count = 0_usize;
    let mut scenario_ids = BTreeSet::new();
    let mut configuration_copy_owners = BTreeMap::<String, String>::new();
    for (index, configuration) in configuration_values.iter().enumerate() {
        let context = format!("configuration_inputs[{index}]");
        let configuration = object(configuration, &context)?;
        let input_kind = text(configuration, "kind", &context)?;
        let copy = configuration
            .get("copy")
            .ok_or_else(|| format!("falta {context}.copy"))?;
        let copy_context = format!("{context}.copy");
        let (copy_id, copy_kind, copy_state, copy_validated) =
            artifact_fields(copy, &copy_context)?;
        if copy_state != "complete" || !copy_validated {
            return Err(format!("{copy_context} no está completo y validado"));
        }
        require_catalog_identity(copy_id, copy, &artifacts, &copy_context)?;
        let copy_path = text(object(copy, &copy_context)?, "path", &copy_context)?;
        match input_kind {
            "general" => {
                general_configuration_count += 1;
                if copy_id != "config-general"
                    || copy_kind != "config_general_toml"
                    || copy_path != "inputs/config/general.toml"
                {
                    return Err(format!("{copy_context} no es la copia general canónica"));
                }
            }
            "scenario" => {
                let scenario_id = copy_id.strip_prefix("config-scenario-").ok_or_else(|| {
                    format!("{copy_context} no identifica una configuración de escenario")
                })?;
                if scenario_id.is_empty()
                    || copy_kind != "config_scenario_toml"
                    || copy_path != format!("inputs/config/scenarios/{scenario_id}.toml")
                {
                    return Err(format!(
                        "{copy_context} no es una copia de escenario canónica"
                    ));
                }
                if !scenario_ids.insert(scenario_id.to_owned()) {
                    return Err(format!(
                        "escenario de configuración duplicado: {scenario_id}"
                    ));
                }
            }
            other => return Err(format!("clase de configuración no admitida: {other}")),
        }
        if let Some(previous) = configuration_copy_owners.insert(copy_id.to_owned(), context) {
            return Err(format!(
                "copia de configuración {copy_id} reutilizada por {previous}"
            ));
        }
    }
    if general_configuration_count != 1 {
        return Err("debe existir exactamente una configuración general".to_owned());
    }
    for (id, artifact) in &artifacts {
        let (_, kind, _, _) = artifact_fields(artifact, &format!("artifact {id}"))?;
        if matches!(kind, "config_general_toml" | "config_scenario_toml")
            && !configuration_copy_owners.contains_key(id)
        {
            return Err(format!("copia de configuración huérfana: {id}"));
        }
    }

    let resolution_evidence = manifest
        .get("resolution_evidence")
        .ok_or_else(|| "falta resolution_evidence".to_owned())?;
    let (resolution_id, resolution_kind, resolution_state, resolution_validated) =
        artifact_fields(resolution_evidence, "resolution_evidence")?;
    let resolution_path = text(
        object(resolution_evidence, "resolution_evidence")?,
        "path",
        "resolution_evidence",
    )?;
    if resolution_id != "resolution-evidence"
        || resolution_kind != "resolution_evidence"
        || resolution_path != "inputs/resolution/evidence.json"
        || resolution_state != "complete"
        || !resolution_validated
    {
        return Err("resolution_evidence no tiene identidad canónica completa".to_owned());
    }
    require_catalog_identity(
        resolution_id,
        resolution_evidence,
        &artifacts,
        "resolution_evidence",
    )?;

    validate_general_output_categories(manifest, &artifacts)?;

    let reference_values = array(manifest, "references", "manifest")?;
    let mut references = BTreeMap::<String, &Value>::new();
    let mut reference_catalogs = BTreeMap::<String, ManifestReferenceCatalog>::new();
    let mut reference_fasta_owners = BTreeMap::<String, String>::new();
    for (index, reference) in reference_values.iter().enumerate() {
        let reference_object = object(reference, &format!("references[{index}]"))?;
        let id = text(reference_object, "id", &format!("references[{index}]"))?;
        if references.insert(id.to_owned(), reference).is_some() {
            return Err(format!("ID de referencia duplicado: {id}"));
        }
        let catalog = manifest_reference_catalog(reference, &format!("references[{index}]"))?;
        reference_catalogs.insert(id.to_owned(), catalog);
        let fasta = reference_object
            .get("fasta")
            .ok_or_else(|| format!("falta references[{index}].fasta"))?;
        let fasta_context = format!("references[{index}].fasta");
        let (fasta_id, fasta_kind, fasta_state, fasta_validated) =
            artifact_fields(fasta, &fasta_context)?;
        if fasta_kind != "reference_fasta" || fasta_state != "complete" || !fasta_validated {
            return Err(format!(
                "{fasta_context} no es un FASTA completo y validado"
            ));
        }
        require_catalog_identity(fasta_id, fasta, &artifacts, &fasta_context)?;
        let authoritative_resource = text(
            reference_object,
            "authoritative_resource",
            &format!("references[{index}]"),
        )?;
        let fasta_object = object(fasta, &fasta_context)?;
        if authoritative_resource != text(fasta_object, "path", &fasta_context)? {
            return Err(format!(
                "references[{index}].authoritative_resource discrepa de su FASTA"
            ));
        }
        if let Some(previous) = reference_fasta_owners.insert(fasta_id.to_owned(), id.to_owned()) {
            return Err(format!(
                "FASTA {fasta_id} reutilizado por referencias {previous} y {id}"
            ));
        }
    }

    let effective_config_id = if let Some(effective_config) = manifest
        .get("effective_config")
        .filter(|value| !value.is_null())
    {
        let context = "effective_config";
        let (id, kind, _, _) = artifact_fields(effective_config, context)?;
        if id != "effective-config" || kind != "effective_config" {
            return Err("effective_config tiene un ID o clase ajenos".to_owned());
        }
        require_catalog_identity(id, effective_config, &artifacts, context)?;
        Some(id.to_owned())
    } else {
        None
    };

    let event_values = array(manifest, "events", "manifest")?;
    let mut event_ids = BTreeSet::new();
    let mut event_manifest_owners = BTreeMap::<String, String>::new();
    for (index, event) in event_values.iter().enumerate() {
        let event_context = format!("events[{index}]");
        let event = object(event, &event_context)?;
        let event_id = text(event, "event_id", &event_context)?;
        if !event_ids.insert(event_id.to_owned()) {
            return Err(format!("ID de evento duplicado: {event_id}"));
        }
        let manifest_link = event
            .get("manifest")
            .ok_or_else(|| format!("falta {event_context}.manifest"))?;
        let link_context = format!("{event_context}.manifest");
        let (manifest_id, manifest_kind, manifest_state, manifest_validated) =
            artifact_fields(manifest_link, &link_context)?;
        let expected_manifest_id = format!("manifest-event-{event_id}");
        if manifest_id != expected_manifest_id
            || manifest_kind != "manifest_event"
            || manifest_state != "complete"
            || !manifest_validated
        {
            return Err(format!("{link_context} tiene identidad o estado ajenos"));
        }
        require_catalog_identity(manifest_id, manifest_link, &artifacts, &link_context)?;
        if let Some(previous) =
            event_manifest_owners.insert(manifest_id.to_owned(), event_id.to_owned())
        {
            return Err(format!(
                "manifiesto {manifest_id} reutilizado por eventos {previous} y {event_id}"
            ));
        }
    }

    let unit_values = array(manifest, "units", "manifest")?;
    validate_manifest_diagnostic_graph(manifest, unit_values)?;
    let mut reduction_units = BTreeMap::<String, &Value>::new();
    let mut event_unit_ids = BTreeSet::new();
    let mut unit_identities = BTreeSet::new();
    let mut execution_units = Vec::new();
    for (index, unit) in unit_values.iter().enumerate() {
        let unit_object = object(unit, &format!("units[{index}]"))?;
        let id = text(unit_object, "id", &format!("units[{index}]"))?;
        let kind = text(unit_object, "kind", &format!("units[{index}]"))?;
        if !unit_identities.insert((kind.to_owned(), id.to_owned())) {
            return Err(format!("unidad duplicada: {kind}/{id}"));
        }
        match kind {
            "reference_reduction" => {
                reduction_units.insert(id.to_owned(), unit);
            }
            "event" => {
                event_unit_ids.insert(id.to_owned());
            }
            "execution" => execution_units.push(unit),
            _ => {}
        }
    }
    if execution_units.len() != 1 {
        return Err("debe existir una única unidad de ejecución".to_owned());
    }
    let execution = object(execution_units[0], "execution unit")?;
    if text(execution, "id", "execution unit")? != text(manifest, "execution_id", "manifest")? {
        return Err("la unidad de ejecución contradice execution_id".to_owned());
    }
    if text(execution, "state", "execution unit")? != manifest_state {
        return Err("el estado de la unidad de ejecución contradice el manifiesto".to_owned());
    }
    if event_unit_ids != event_ids {
        return Err(format!(
            "las unidades de evento no coinciden 1:1 con los eventos: units={event_unit_ids:?}, events={event_ids:?}"
        ));
    }

    let reduction_values = array(manifest, "reference_reductions", "manifest")?;
    let mut reductions = BTreeMap::<String, &Value>::new();
    for (index, reduction) in reduction_values.iter().enumerate() {
        let reduction_object = object(reduction, &format!("reference_reductions[{index}]"))?;
        let id = text(
            reduction_object,
            "derived_reference_id",
            &format!("reference_reductions[{index}]"),
        )?;
        if reductions.insert(id.to_owned(), reduction).is_some() {
            return Err(format!("registro de reducción duplicado: {id}"));
        }
    }

    let mut claimed_artifacts = BTreeSet::new();
    let mut claimed_references = BTreeSet::new();
    for (derived_id, unit_value) in &reduction_units {
        let context = format!("reference_reduction {derived_id}");
        let unit = object(unit_value, &context)?;
        let unit_state = text(unit, "state", &context)?;
        let exportable = boolean(unit, "exportable", &context)?;
        let diagnostic = unit
            .get("diagnostic_code")
            .ok_or_else(|| format!("falta {context}.diagnostic_code"))?;
        validate_reduction_unit_state(
            derived_id,
            unit_state,
            exportable,
            diagnostic,
            manifest_state,
        )?;

        let expected_fasta_id = format!("reference-fasta-{derived_id}");
        let expected_map_id = format!("reference-map-{derived_id}");
        match reductions.get(derived_id) {
            Some(reduction_value) => validate_reduction_record(
                derived_id,
                reduction_value,
                unit_state,
                &expected_fasta_id,
                &expected_map_id,
                &artifacts,
                &references,
                &reference_catalogs,
                &mut claimed_artifacts,
                &mut claimed_references,
            )?,
            None => validate_partial_failed_reduction(
                derived_id,
                unit_state,
                &expected_fasta_id,
                &expected_map_id,
                &artifacts,
                &references,
                &mut claimed_artifacts,
            )?,
        }
    }

    for derived_id in reductions.keys() {
        if !reduction_units.contains_key(derived_id) {
            return Err(format!(
                "registro de reducción sin unidad correspondiente: {derived_id}"
            ));
        }
    }

    for (reference_id, reference) in &references {
        let reference = object(reference, &format!("reference {reference_id}"))?;
        let origin = reference
            .get("origin")
            .ok_or_else(|| format!("falta reference {reference_id}.origin"))?;
        let origin = object(origin, &format!("reference {reference_id}.origin"))?;
        if text(origin, "kind", &format!("reference {reference_id}.origin"))? == "reduced"
            && !claimed_references.contains(reference_id)
        {
            return Err(format!(
                "referencia reducida huérfana, sin reducción propietaria: {reference_id}"
            ));
        }
    }

    for (id, artifact) in &artifacts {
        let (_, kind, state, _) = artifact_fields(artifact, &format!("artifact {id}"))?;
        let reduction_shaped = reduction_shaped_artifact(id, kind);
        if (state == "failed" || reduction_shaped) && !claimed_artifacts.contains(id) {
            return Err(format!("artefacto de reducción huérfano: {id}"));
        }
        match kind {
            "config_general_toml" | "config_scenario_toml"
                if configuration_copy_owners.contains_key(id) => {}
            "resolution_evidence" if id == resolution_id => {}
            "effective_config" if effective_config_id.as_deref() == Some(id) => {}
            "reference_fasta"
                if reference_fasta_owners.contains_key(id) || claimed_artifacts.contains(id) => {}
            "reference_reduction_map" if claimed_artifacts.contains(id) => {}
            "manifest_event" if event_manifest_owners.contains_key(id) => {}
            "config_general_toml"
            | "config_scenario_toml"
            | "resolution_evidence"
            | "effective_config"
            | "reference_fasta"
            | "reference_reduction_map"
            | "manifest_event" => return Err(format!("artefacto sin propietario semántico: {id}")),
            other => {
                return Err(format!(
                    "clase de artefacto no admitida en F1: {other}/{id}"
                ));
            }
        }
    }
    Ok(())
}

fn validate_linked_manifest_documents(path: &Path, document: &Value) -> Result<(), String> {
    let manifest = object(document, "manifest")?;
    let resolution_link = manifest
        .get("resolution_evidence")
        .ok_or_else(|| "falta resolution_evidence enlazada".to_owned())?;
    let resolution_path = text(
        object(resolution_link, "resolution_evidence")?,
        "path",
        "resolution_evidence",
    )?;
    let root = path.parent().unwrap_or_else(|| Path::new("."));
    let resolution_bytes = read_package_file_nofollow(root, Path::new(resolution_path))
        .map_err(|error| format!("no se pudo releer {resolution_path}: {error}"))?;
    validate_reloaded_artifact_identity(resolution_link, &resolution_bytes, "resolution_evidence")?;
    let resolution_evidence: Value = serde_json::from_slice(&resolution_bytes)
        .map_err(|error| format!("{resolution_path} no es JSON: {error}"))?;
    crate::manifest_schema::validate_resolution_evidence_structure(&resolution_evidence)
        .map_err(|detail| format!("{resolution_path} no cumple el esquema: {detail}"))?;
    let (authoritative_domain, authoritative_requests) =
        parse_resolution_evidence(&resolution_evidence)?;
    let Some(effective_link) = manifest
        .get("effective_config")
        .filter(|value| !value.is_null())
    else {
        return Ok(());
    };
    let effective_path = text(
        object(effective_link, "effective_config")?,
        "path",
        "effective_config",
    )?;
    let effective_bytes = read_package_file_nofollow(root, Path::new(effective_path))
        .map_err(|error| format!("no se pudo releer {effective_path}: {error}"))?;
    let effective: Value = serde_json::from_slice(&effective_bytes)
        .map_err(|error| format!("{effective_path} no es JSON: {error}"))?;
    crate::manifest_schema::validate_effective_config_structure(&effective)
        .map_err(|detail| format!("{effective_path} no cumple el esquema: {detail}"))?;
    validate_effective_config_semantics(&effective)?;

    if manifest.get("route") != effective.get("route") {
        return Err("route discrepa entre manifiesto y configuración efectiva".to_owned());
    }
    if manifest.get("normalizations") != effective.get("normalizations") {
        return Err("normalizations discrepa entre manifiesto y configuración efectiva".to_owned());
    }
    if manifest.get("reference_reductions") != effective.get("reduced_references") {
        return Err("reference_reductions discrepa de reduced_references efectivos".to_owned());
    }

    let effective_reference = effective
        .get("reference")
        .ok_or_else(|| "falta effective.reference".to_owned())?;
    let effective_reference_id =
        independent_value_text(effective_reference, "id", "effective.reference")?;
    let general_reference = array(manifest, "references", "manifest")?
        .iter()
        .find(|reference| {
            reference
                .get("id")
                .and_then(Value::as_str)
                .is_some_and(|id| id == effective_reference_id)
        })
        .ok_or_else(|| "la referencia efectiva no existe en el manifiesto".to_owned())?;
    if general_reference != effective_reference {
        return Err("la referencia fuente discrepa de la configuración efectiva".to_owned());
    }
    validate_authoritative_fastas(root, manifest)?;

    let effective_scenarios = string_array(&effective, "scenario_ids", "effective")?;
    let observed_scenarios = array(manifest, "configuration_inputs", "manifest")?
        .iter()
        .filter(|input| input.get("kind").and_then(Value::as_str) == Some("scenario"))
        .map(|input| {
            input["copy"]["id"]
                .as_str()
                .and_then(|id| id.strip_prefix("config-scenario-"))
                .map(str::to_owned)
        })
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| "entrada de escenario sin identidad canónica".to_owned())?;
    if observed_scenarios != effective_scenarios {
        return Err("scenario_ids discrepa de configuration_inputs".to_owned());
    }

    let randomness = object(
        effective
            .get("randomness")
            .ok_or_else(|| "falta effective.randomness".to_owned())?,
        "effective.randomness",
    )?;
    let root_hex = text(randomness, "root_seed", "effective.randomness")?;
    let root_seed =
        RootSeed::parse("root_seed", root_hex).map_err(|diagnostic| diagnostic.render_human())?;
    validate_effective_against_configuration_inputs(
        root,
        manifest,
        &effective,
        root_seed,
        &authoritative_domain,
        &authoritative_requests,
    )?;
    validate_reduction_randomness(manifest, root_seed, root_hex)?;
    validate_linked_events(
        root,
        manifest,
        &effective,
        effective_reference,
        effective_reference_id,
        root_seed,
    )
}

fn validate_reloaded_artifact_identity(
    link: &Value,
    bytes: &[u8],
    context: &str,
) -> Result<(), String> {
    let link = object(link, context)?;
    let expected_size = unsigned(link, "size_bytes", context)?;
    let checksum = object(
        link.get("checksum")
            .ok_or_else(|| format!("falta {context}.checksum"))?,
        &format!("{context}.checksum"),
    )?;
    let expected_digest = text(checksum, "digest", &format!("{context}.checksum"))?;
    let mut hasher = Sha256Writer::default();
    hasher
        .update(bytes)
        .map_err(|diagnostic| diagnostic.render_human())?;
    let identity = hasher.finish(context);
    if identity.size_bytes() != expected_size || identity.digest() != expected_digest {
        return Err(format!(
            "{context} releído contradice tamaño o SHA-256 enlazados"
        ));
    }
    Ok(())
}

#[derive(Debug)]
struct ManifestEntropy([u8; 32]);

impl EntropySource for ManifestEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        *destination = self.0;
        Ok(())
    }
}

fn validate_effective_against_configuration_inputs(
    root: &Path,
    manifest: &serde_json::Map<String, Value>,
    effective: &Value,
    root_seed: RootSeed,
    resolution_domain: &ResolutionDomain,
    resolution_requests: &[ResolutionRequest],
) -> Result<(), String> {
    let inputs = array(manifest, "configuration_inputs", "manifest")?;
    let general_input = inputs
        .iter()
        .find(|input| input.get("kind").and_then(Value::as_str) == Some("general"))
        .ok_or_else(|| "falta configuración general enlazada".to_owned())?;
    let general_path = independent_value_text(general_input, "original_path", "general input")?;
    if independent_value_text(effective, "source_path", "effective")? != general_path {
        return Err("effective.source_path discrepa de la configuración general".to_owned());
    }
    let general_copy = general_input
        .get("copy")
        .ok_or_else(|| "falta copia de configuración general".to_owned())?;
    let general_relative = independent_value_text(general_copy, "path", "general copy")?;
    let general_bytes = read_package_file_nofollow(root, Path::new(general_relative))?;
    let admitted = admit_toml_header(
        Path::new(general_path),
        general_bytes,
        ExpectedConfigKind::General,
    )
    .map_err(|diagnostic| diagnostic.render_human())?;
    let DeclaredConfig::General(general) =
        parse_declared_config(admitted).map_err(|diagnostic| diagnostic.render_human())?
    else {
        return Err("la copia general no contiene una configuración general".to_owned());
    };

    let mut scenarios = Vec::new();
    for input in inputs
        .iter()
        .filter(|input| input.get("kind").and_then(Value::as_str) == Some("scenario"))
    {
        let original = independent_value_text(input, "original_path", "scenario input")?;
        let copy = input
            .get("copy")
            .ok_or_else(|| "falta copia de escenario".to_owned())?;
        let relative = independent_value_text(copy, "path", "scenario copy")?;
        let bytes = read_package_file_nofollow(root, Path::new(relative))?;
        let admitted = admit_toml_header(Path::new(original), bytes, ExpectedConfigKind::Scenario)
            .map_err(|diagnostic| diagnostic.render_human())?;
        let DeclaredConfig::Scenario(scenario) =
            parse_declared_config(admitted).map_err(|diagnostic| diagnostic.render_human())?
        else {
            return Err("una copia de escenario no contiene un escenario".to_owned());
        };
        scenarios.push(scenario);
    }
    let scenario_refs = scenarios.iter().map(Box::as_ref).collect::<Vec<_>>();
    if effective.get("resolution_domain") != Some(&resolution_domain_json(resolution_domain)) {
        return Err("effective.resolution_domain discrepa de la evidencia autoritativa".to_owned());
    }
    if effective["resolution_evidence_id"].as_str() != Some("resolution-evidence") {
        return Err("effective no enlaza la evidencia de resolución canónica".to_owned());
    }
    let mut entropy = ManifestEntropy(*root_seed.as_bytes());
    let reconstructed = build_effective_config_with_entropy_and_domain(
        &general,
        &scenario_refs,
        resolution_requests,
        &mut entropy,
        resolution_domain,
    )
    .map_err(|diagnostic| diagnostic.render_human())?;
    let expected = effective_config_base(&reconstructed);
    for field in [
        "source_path",
        "route",
        "scenario_ids",
        "event_ids",
        "individual_ids",
        "parameter_order_version",
        "resolution_domain",
        "randomness",
        "normalizations",
        "controlled_frequency_sets",
    ] {
        if effective.get(field) != expected.get(field) {
            return Err(format!(
                "effective.{field} discrepa de las copias TOML autoritativas"
            ));
        }
    }
    validate_effective_values(
        effective,
        &general,
        &scenario_refs,
        resolution_domain,
        &expected,
    )
}

fn validate_effective_values(
    effective: &Value,
    general: &crate::configuration::DeclaredGeneralConfig,
    scenarios: &[&crate::configuration::DeclaredScenarioConfig],
    domain: &ResolutionDomain,
    reconstructed: &Value,
) -> Result<(), String> {
    let values = effective["effective_values"]
        .as_array()
        .ok_or_else(|| "effective.effective_values no es array".to_owned())?;
    let mut grouped = BTreeMap::<ParameterKey, Vec<crate::rules::Rule>>::new();
    for rule in general.rules() {
        grouped
            .entry(rule.parameter())
            .or_default()
            .push(rule.clone());
    }
    let rule_sets = grouped
        .into_iter()
        .map(|(parameter, rules)| {
            RuleSet::new_with_fallback(parameter, rules, parameter == ParameterKey::Ploidy)
                .map(|set| (parameter, set))
                .map_err(|diagnostic| diagnostic.render_human())
        })
        .collect::<Result<BTreeMap<_, _>, _>>()?;
    let mut observed = BTreeMap::<(String, ParameterKey), (&Value, ResolutionContext)>::new();
    let mut previous: Option<(&str, usize)> = None;
    for (index, value) in values.iter().enumerate() {
        let context = format!("effective.effective_values[{index}]");
        let value_object = object(value, &context)?;
        let object_id = text(value_object, "object_id", &context)?;
        let parameter_text = text(value_object, "parameter", &context)?;
        let parameter = ParameterKey::parse(parameter_text)
            .ok_or_else(|| format!("{context}.parameter no pertenece al catálogo"))?;
        let identity = (object_id, parameter.canonical_rank());
        if previous.is_some_and(|prior| prior >= identity) {
            return Err(
                "effective.effective_values no está ordenado o repite identidad".to_owned(),
            );
        }
        previous = Some(identity);
        let resolution_context = parse_resolution_context(
            value_object
                .get("resolution_context")
                .ok_or_else(|| format!("falta {context}.resolution_context"))?,
            &context,
        )?;
        validate_resolution_request(
            general,
            scenarios,
            &ResolutionRequest {
                object_id: object_id.to_owned(),
                parameter,
                context: resolution_context.clone(),
            },
            domain,
        )
        .map_err(|diagnostic| format!("{context}: {}", diagnostic.render_human()))?;
        let rule = rule_sets
            .get(&parameter)
            .ok_or_else(|| format!("{context} usa un parámetro sin reglas"))?
            .resolve(&resolution_context)
            .map_err(|diagnostic| diagnostic.render_human())?;
        let expected = json!({
            "object_id": object_id,
            "parameter": rule.parameter().as_str(),
            "resolution_context": resolution_context_json(&resolution_context),
            "value": rule_value_json(rule.value()),
            "winning_rule_id": rule.id(),
            "selectors": selectors_json(rule.selectors()),
            "provenance": {
                "path": rule.provenance().path,
                "location": rule.provenance().location,
            },
        });
        if value != &expected {
            return Err(format!(
                "{context} contradice valor, selector o procedencia de la regla ganadora"
            ));
        }
        if parameter == ParameterKey::Ploidy {
            let sequence = effective["reference"]["sequences"]
                .as_array()
                .and_then(|sequences| {
                    sequences
                        .iter()
                        .find(|sequence| sequence["id"].as_str() == Some(object_id))
                })
                .ok_or_else(|| format!("{context} aplica ploidía a una secuencia inexistente"))?;
            if sequence.get("ploidy") != value.get("value") {
                return Err(format!("{context} discrepa de reference.sequences.ploidy"));
            }
        }
        observed.insert(
            (object_id.to_owned(), parameter),
            (value, resolution_context),
        );
    }
    let reconstructed_values = reconstructed["effective_values"]
        .as_array()
        .ok_or_else(|| "reconstrucción efectiva sin valores".to_owned())?;
    if observed.len() != reconstructed_values.len() {
        return Err(
            "effective.effective_values no coincide con el conjunto autoritativo de solicitudes"
                .to_owned(),
        );
    }
    for expected in reconstructed_values {
        let object_id = independent_value_text(expected, "object_id", "expected value")?;
        let parameter = ParameterKey::parse(independent_value_text(
            expected,
            "parameter",
            "expected value",
        )?)
        .ok_or_else(|| "parámetro reconstruido ajeno al catálogo".to_owned())?;
        let exact = observed
            .get(&(object_id.to_owned(), parameter))
            .is_some_and(|(observed, _)| *observed == expected);
        if !exact {
            return Err(format!(
                "effective.effective_values omite o altera el valor automático {object_id}:{}",
                parameter.as_str()
            ));
        }
    }
    Ok(())
}

fn parse_resolution_context(value: &Value, label: &str) -> Result<ResolutionContext, String> {
    let context = value
        .as_object()
        .ok_or_else(|| format!("{label}.resolution_context no es objeto"))?;
    let optional_text = |field: &str| -> Result<Option<String>, String> {
        match context.get(field) {
            Some(Value::Null) => Ok(None),
            Some(Value::String(value)) if !value.is_empty() => Ok(Some(value.clone())),
            _ => Err(format!(
                "{label}.resolution_context.{field} no es texto/null"
            )),
        }
    };
    let optional_u64 = |field: &str| -> Result<Option<u64>, String> {
        match context.get(field) {
            Some(Value::Null) => Ok(None),
            Some(value) => value
                .as_u64()
                .map(Some)
                .ok_or_else(|| format!("{label}.resolution_context.{field} no es entero/null")),
            None => Err(format!("falta {label}.resolution_context.{field}")),
        }
    };
    Ok(ResolutionContext {
        individual: optional_text("individual")?,
        population: optional_text("population")?,
        subpopulation: optional_text("subpopulation")?,
        phase: optional_text("phase")?,
        generation: optional_u64("generation")?
            .map(|value| {
                u32::try_from(value)
                    .map_err(|_| format!("{label}.resolution_context.generation no cabe en u32"))
            })
            .transpose()?,
        reference: optional_text("reference")?,
        sequence: optional_text("sequence")?,
        position: optional_u64("position")?,
        locus: optional_text("locus")?,
        group: optional_text("group")?,
        block: optional_text("block")?,
        technical_event: optional_text("technical_event")?,
    })
}

fn independent_value_text<'a>(
    value: &'a Value,
    field: &str,
    context: &str,
) -> Result<&'a str, String> {
    value
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("falta o no es texto {context}.{field}"))
}

fn string_array(value: &Value, field: &str, context: &str) -> Result<Vec<String>, String> {
    value
        .get(field)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("falta o no es array {context}.{field}"))?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            value
                .as_str()
                .map(str::to_owned)
                .ok_or_else(|| format!("{context}.{field}[{index}] no es texto"))
        })
        .collect()
}

fn validate_effective_config_semantics(effective: &Value) -> Result<(), String> {
    validate_normalizations(
        effective
            .get("normalizations")
            .ok_or_else(|| "falta effective.normalizations".to_owned())?,
        "effective.normalizations",
    )?;
    let reference = effective
        .get("reference")
        .ok_or_else(|| "falta effective.reference".to_owned())?;
    manifest_reference_catalog(reference, "effective.reference")?;
    Ok(())
}

fn validate_reduction_randomness(
    manifest: &serde_json::Map<String, Value>,
    root_seed: RootSeed,
    root_hex: &str,
) -> Result<(), String> {
    for (index, reduction) in array(manifest, "reference_reductions", "manifest")?
        .iter()
        .enumerate()
    {
        let context = format!("reference_reductions[{index}]");
        let randomness = object(
            reduction
                .get("randomness")
                .ok_or_else(|| format!("falta {context}.randomness"))?,
            &format!("{context}.randomness"),
        )?;
        if text(randomness, "root_seed", &format!("{context}.randomness"))? != root_hex {
            return Err(format!(
                "{context}.randomness.root_seed discrepa de la raíz efectiva"
            ));
        }
        let source = independent_value_text(reduction, "source_reference_id", &context)?;
        let derived = independent_value_text(reduction, "derived_reference_id", &context)?;
        let algorithm = independent_value_text(reduction, "algorithm", &context)?;
        let request = SubseedRequest::from_scope(
            "reference_reduction",
            &[
                ("source_reference_id", source),
                ("derived_reference_id", derived),
                ("algorithm", algorithm),
            ],
        )
        .map_err(|diagnostic| diagnostic.render_human())?;
        let expected = derive_subseed(root_seed, &request)
            .map_err(|diagnostic| diagnostic.render_human())?
            .to_hex();
        if text(randomness, "subseed", &format!("{context}.randomness"))? != expected {
            return Err(format!(
                "{context}.randomness.subseed no se deriva del catálogo v3"
            ));
        }
    }
    Ok(())
}

fn validate_linked_events(
    root: &Path,
    manifest: &serde_json::Map<String, Value>,
    effective: &Value,
    effective_reference: &Value,
    reference_id: &str,
    root_seed: RootSeed,
) -> Result<(), String> {
    let effective_event_ids = string_array(effective, "event_ids", "effective")?;
    let event_subseeds = effective["randomness"]["event_subseeds"]
        .as_array()
        .ok_or_else(|| "effective.randomness.event_subseeds no es array".to_owned())?;
    if event_subseeds.len() != effective_event_ids.len() {
        return Err("event_ids y event_subseeds no tienen la misma cardinalidad".to_owned());
    }
    let general_events = array(manifest, "events", "manifest")?;
    if text(manifest, "state", "manifest")? == "complete"
        && general_events.len() != effective_event_ids.len()
    {
        return Err("events no coincide en cardinalidad con event_ids".to_owned());
    }
    let mut previous: Option<&str> = None;
    let mut expected_by_id = BTreeMap::new();
    for (index, event_id) in effective_event_ids.iter().enumerate() {
        if previous.is_some_and(|value| value >= event_id.as_str()) {
            return Err("event_ids no está en orden canónico o contiene duplicados".to_owned());
        }
        previous = Some(event_id);
        let effective_subseed = &event_subseeds[index];
        if independent_value_text(effective_subseed, "event_id", "event_subseed")? != event_id {
            return Err("event_subseeds no sigue event_ids".to_owned());
        }
        let request = SubseedRequest::from_scope(
            "technical_event",
            &[("event_id", event_id), ("reference_id", reference_id)],
        )
        .map_err(|diagnostic| diagnostic.render_human())?;
        let expected_subseed = derive_subseed(root_seed, &request)
            .map_err(|diagnostic| diagnostic.render_human())?
            .to_hex();
        if independent_value_text(effective_subseed, "subseed", "event_subseed")?
            != expected_subseed
        {
            return Err(format!(
                "subsemilla efectiva de {event_id} no deriva del catálogo v3"
            ));
        }
        expected_by_id.insert(event_id.as_str(), expected_subseed);
    }
    let mut observed_event_ids = BTreeSet::new();
    for general_event in general_events {
        let event_id = independent_value_text(general_event, "event_id", "event")?;
        if !observed_event_ids.insert(event_id) {
            return Err(format!("evento general duplicado: {event_id}"));
        }
        let expected_subseed = expected_by_id
            .get(event_id)
            .ok_or_else(|| format!("evento general ajeno a event_ids: {event_id}"))?;
        if independent_value_text(general_event, "event_id", "event")? != event_id
            || independent_value_text(general_event, "subseed", "event")?
                != expected_subseed.as_str()
        {
            return Err(format!(
                "evento general {event_id} discrepa de la configuración efectiva"
            ));
        }
        let manifest_link = general_event
            .get("manifest")
            .ok_or_else(|| format!("falta manifiesto de evento {event_id}"))?;
        let event_path = independent_value_text(manifest_link, "path", "event.manifest")?;
        let bytes = read_package_file_nofollow(root, Path::new(event_path))
            .map_err(|error| format!("no se pudo releer {event_path}: {error}"))?;
        let event: Value = serde_json::from_slice(&bytes)
            .map_err(|error| format!("{event_path} no es JSON: {error}"))?;
        crate::manifest_schema::validate_event_manifest_structure(&event)
            .map_err(|detail| format!("{event_path} no cumple el esquema: {detail}"))?;
        if independent_value_text(&event, "event_id", "event manifest")? != event_id
            || independent_value_text(&event, "subseed", "event manifest")?
                != expected_subseed.as_str()
            || event["units"][0]["id"].as_str() != Some(event_id)
        {
            return Err(format!(
                "manifiesto de evento {event_id} contradice ID/subsemilla/unidad"
            ));
        }
        let expected_reference_link = json!({
            "id": effective_reference["id"],
            "fasta": effective_reference["fasta"],
            "indexes": effective_reference["indexes"],
        });
        if event["references"][0] != expected_reference_link {
            return Err(format!(
                "manifiesto de evento {event_id} enlaza otra referencia"
            ));
        }
    }
    Ok(())
}

fn read_package_file_nofollow(root: &Path, relative: &Path) -> Result<Vec<u8>, String> {
    if relative.is_absolute()
        || relative
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(format!("ruta no canónica: {}", relative.display()));
    }
    let mut directory = openat(
        CWD,
        root,
        OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
        Mode::empty(),
    )
    .map_err(|error| format!("no se puede abrir raíz no-follow: {error}"))?;
    let mut components = relative.components().peekable();
    let mut file = None;
    while let Some(Component::Normal(name)) = components.next() {
        if components.peek().is_some() {
            directory = openat(
                &directory,
                name,
                OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
                Mode::empty(),
            )
            .map_err(|error| {
                format!("no se puede abrir directorio intermedio no-follow: {error}")
            })?;
        } else {
            file = Some(
                openat(
                    &directory,
                    name,
                    OFlags::RDONLY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
                    Mode::empty(),
                )
                .map_err(|error| format!("no se puede abrir fichero no-follow: {error}"))?,
            );
        }
    }
    let file = file.ok_or_else(|| "ruta de fichero vacía".to_owned())?;
    let metadata = fstat(&file).map_err(|error| format!("no se puede inspeccionar fd: {error}"))?;
    if FileType::from_raw_mode(metadata.st_mode) != FileType::RegularFile {
        return Err("la ruta enlazada no es un fichero regular".to_owned());
    }
    let mut bytes = Vec::new();
    File::from(file)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("no se puede leer fd: {error}"))?;
    Ok(bytes)
}

#[derive(Debug, Eq, PartialEq)]
struct ObservedFastaSequence {
    header: Vec<u8>,
    id: String,
    length: u64,
    unknown_base_count: u64,
}

fn validate_authoritative_fastas(
    root: &Path,
    manifest: &serde_json::Map<String, Value>,
) -> Result<(), String> {
    for (reference_index, reference) in array(manifest, "references", "manifest")?
        .iter()
        .enumerate()
    {
        let context = format!("references[{reference_index}]");
        let reference = object(reference, &context)?;
        let relative = text(reference, "authoritative_resource", &context)?;
        let bytes = read_package_file_nofollow(root, Path::new(relative))
            .map_err(|error| format!("{context} no puede releer su FASTA: {error}"))?;
        let observed = parse_authoritative_fasta(&bytes, &context)?;
        let declared = array(reference, "sequences", &context)?;
        if observed.len() != declared.len() {
            return Err(format!(
                "{context}.sequences no coincide en cardinalidad con el FASTA"
            ));
        }
        for (index, (observed, declared)) in observed.iter().zip(declared).enumerate() {
            let sequence_context = format!("{context}.sequences[{index}]");
            let declared = object(declared, &sequence_context)?;
            let declared_header = reversible_text_bytes(
                declared
                    .get("header")
                    .ok_or_else(|| format!("falta {sequence_context}.header"))?,
                &format!("{sequence_context}.header"),
            )?;
            if text(declared, "id", &sequence_context)? != observed.id
                || declared_header != observed.header
                || unsigned(declared, "length", &sequence_context)? != observed.length
                || unsigned(declared, "unknown_base_count", &sequence_context)?
                    != observed.unknown_base_count
            {
                return Err(format!(
                    "{sequence_context} discrepa de los bytes del FASTA autoritativo"
                ));
            }
        }
    }
    Ok(())
}

fn parse_authoritative_fasta(
    bytes: &[u8],
    context: &str,
) -> Result<Vec<ObservedFastaSequence>, String> {
    if bytes.starts_with(&[0x1f, 0x8b]) {
        return Err(format!("{context} enlaza un FASTA comprimido"));
    }
    let mut records = Vec::new();
    let mut current: Option<ObservedFastaSequence> = None;
    for raw_line in bytes.split_inclusive(|byte| *byte == b'\n') {
        let line = raw_line
            .strip_suffix(b"\n")
            .unwrap_or(raw_line)
            .strip_suffix(b"\r")
            .unwrap_or_else(|| raw_line.strip_suffix(b"\n").unwrap_or(raw_line));
        if let Some(raw_header) = line.strip_prefix(b">") {
            if let Some(record) = current.take() {
                if record.length == 0 {
                    return Err(format!("{context} contiene una secuencia vacía"));
                }
                records.push(record);
            }
            let header = FastaHeader::parse(raw_header)
                .map_err(|error| format!("{context} contiene cabecera inválida: {error:?}"))?;
            current = Some(ObservedFastaSequence {
                header: raw_header.to_vec(),
                id: header.id_text().to_owned(),
                length: 0,
                unknown_base_count: 0,
            });
        } else {
            let record = current
                .as_mut()
                .ok_or_else(|| format!("{context} contiene bases antes de la cabecera"))?;
            for byte in line {
                let base = byte.to_ascii_uppercase();
                if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                    return Err(format!("{context} contiene una base FASTA no admitida"));
                }
                record.length = record
                    .length
                    .checked_add(1)
                    .ok_or_else(|| format!("{context} desborda la longitud FASTA"))?;
                record.unknown_base_count = record
                    .unknown_base_count
                    .checked_add(u64::from(base == b'N'))
                    .ok_or_else(|| format!("{context} desborda el recuento de N"))?;
            }
        }
    }
    let record = current.ok_or_else(|| format!("{context} enlaza un FASTA vacío"))?;
    if record.length == 0 {
        return Err(format!("{context} contiene una secuencia vacía"));
    }
    records.push(record);
    Ok(records)
}

fn validate_normalizations(value: &Value, context: &str) -> Result<(), String> {
    let collections = value
        .as_array()
        .ok_or_else(|| format!("{context} no es un array"))?;
    let mut parameters = BTreeSet::new();
    for (index, collection) in collections.iter().enumerate() {
        let item_context = format!("{context}[{index}]");
        let collection = object(collection, &item_context)?;
        let parameter = text(collection, "parameter", &item_context)?;
        if !parameters.insert(parameter) {
            return Err(format!("{context} repite el parámetro {parameter}"));
        }
        let tolerance = finite_number(collection, "tolerance", &item_context)?;
        if !(0.0..=MAX_NORMALIZATION_TOLERANCE).contains(&tolerance) {
            return Err(format!("{item_context}.tolerance queda fuera de DTC-049"));
        }
        let values = array(collection, "values", &item_context)?;
        let mut previous: Option<(&str, &str)> = None;
        let mut before = Vec::with_capacity(values.len());
        let mut declared_after = Vec::with_capacity(values.len());
        for (value_index, value) in values.iter().enumerate() {
            let value_context = format!("{item_context}.values[{value_index}]");
            let value = object(value, &value_context)?;
            let identity = (
                text(value, "id", &value_context)?,
                text(value, "category", &value_context)?,
            );
            if previous.is_some_and(|prior| prior >= identity) {
                return Err(format!(
                    "{item_context}.values no está ordenado o repite identidad"
                ));
            }
            previous = Some(identity);
            let raw = finite_number(value, "before", &value_context)?;
            let normalized = finite_number(value, "after", &value_context)?;
            if !(0.0..=1.0).contains(&raw) || !(0.0..=1.0).contains(&normalized) {
                return Err(format!("{value_context} queda fuera de [0,1]"));
            }
            before.push(raw);
            declared_after.push(normalized);
        }
        let sum_before = neumaier(before.iter().copied());
        let delta_before = (sum_before - 1.0).abs();
        if !sum_before.is_finite()
            || sum_before <= 0.0
            || delta_before > tolerance
            || finite_number(collection, "sum_before", &item_context)?.to_bits()
                != sum_before.to_bits()
            || finite_number(collection, "delta_before", &item_context)?.to_bits()
                != delta_before.to_bits()
        {
            return Err(format!("{item_context} contradice la suma Neumaier previa"));
        }
        let expected_after = before
            .iter()
            .map(|value| *value / sum_before)
            .collect::<Vec<_>>();
        if declared_after
            .iter()
            .zip(&expected_after)
            .any(|(declared, expected)| declared.to_bits() != expected.to_bits())
        {
            return Err(format!(
                "{item_context}.values.after no aplica la normalización"
            ));
        }
        let sum_after = neumaier(expected_after.iter().copied());
        let delta_after = (sum_after - 1.0).abs();
        if delta_after > tolerance
            || finite_number(collection, "sum_after", &item_context)?.to_bits()
                != sum_after.to_bits()
            || finite_number(collection, "delta_after", &item_context)?.to_bits()
                != delta_after.to_bits()
        {
            return Err(format!(
                "{item_context} contradice la suma Neumaier posterior"
            ));
        }
        let changed = before
            .iter()
            .zip(&expected_after)
            .any(|(raw, normalized)| raw.to_bits() != normalized.to_bits());
        let expected_decision = if changed {
            "normalized"
        } else {
            "accepted_without_change"
        };
        if text(collection, "decision", &item_context)? != expected_decision {
            return Err(format!("{item_context}.decision contradice los valores"));
        }
    }
    Ok(())
}

fn finite_number(
    object: &serde_json::Map<String, Value>,
    field: &str,
    context: &str,
) -> Result<f64, String> {
    let value = object
        .get(field)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("falta o no es número {context}.{field}"))?;
    value
        .is_finite()
        .then_some(value)
        .ok_or_else(|| format!("{context}.{field} no es finito"))
}

fn validate_general_output_categories(
    manifest: &serde_json::Map<String, Value>,
    artifacts: &BTreeMap<String, &Value>,
) -> Result<(), String> {
    let outputs = manifest
        .get("outputs")
        .ok_or_else(|| "falta manifest.outputs".to_owned())?;
    let outputs = object(outputs, "outputs")?;
    let observed_minimum = output_id_set(outputs, "minimum")?;
    let observed_explicit = output_id_set(outputs, "explicitly_requested")?;
    let observed_default = output_id_set(outputs, "default_enabled")?;
    if !observed_default.is_empty() {
        return Err("outputs.default_enabled debe estar vacío en F1".to_owned());
    }

    let mut expected_minimum = BTreeSet::from(["manifest-general".to_owned()]);
    let mut expected_explicit = BTreeSet::new();
    for (id, artifact) in artifacts {
        let (_, kind, _, _) = artifact_fields(artifact, &format!("artifact {id}"))?;
        if reduction_shaped_artifact(id, kind) {
            expected_explicit.insert(id.to_owned());
        } else {
            expected_minimum.insert(id.to_owned());
        }
    }
    if observed_minimum != expected_minimum || observed_explicit != expected_explicit {
        return Err(format!(
            "las categorías de salida no resuelven exactamente el catálogo: minimum={observed_minimum:?}/{expected_minimum:?}, explicitly_requested={observed_explicit:?}/{expected_explicit:?}"
        ));
    }
    Ok(())
}

fn output_id_set(
    outputs: &serde_json::Map<String, Value>,
    field: &str,
) -> Result<BTreeSet<String>, String> {
    let values = array(outputs, field, "outputs")?;
    let mut ids = BTreeSet::new();
    for (index, value) in values.iter().enumerate() {
        let id = value
            .as_str()
            .ok_or_else(|| format!("outputs.{field}[{index}] no es texto"))?;
        if !ids.insert(id.to_owned()) {
            return Err(format!("salida duplicada en outputs.{field}: {id}"));
        }
    }
    Ok(ids)
}

fn reduction_shaped_artifact(id: &str, kind: &str) -> bool {
    kind == "reference_reduction_map"
        || (kind == "reference_fasta" && id.starts_with("reference-fasta-"))
}

#[allow(clippy::too_many_arguments)]
fn validate_reduction_record(
    derived_id: &str,
    reduction_value: &Value,
    unit_state: &str,
    expected_fasta_id: &str,
    expected_map_id: &str,
    artifacts: &BTreeMap<String, &Value>,
    references: &BTreeMap<String, &Value>,
    reference_catalogs: &BTreeMap<String, ManifestReferenceCatalog>,
    claimed_artifacts: &mut BTreeSet<String>,
    claimed_references: &mut BTreeSet<String>,
) -> Result<(), String> {
    let context = format!("reference_reduction {derived_id}");
    let reduction = object(reduction_value, &context)?;
    let source_reference_id = text(reduction, "source_reference_id", &context)?;
    let source_fasta = reduction
        .get("source_fasta")
        .ok_or_else(|| format!("falta {context}.source_fasta"))?;
    let derived_fasta = reduction
        .get("derived_fasta")
        .ok_or_else(|| format!("falta {context}.derived_fasta"))?;
    let map = reduction
        .get("map")
        .ok_or_else(|| format!("falta {context}.map"))?;

    let (source_id, source_kind, source_state, source_validated) =
        artifact_fields(source_fasta, &format!("{context}.source_fasta"))?;
    if source_kind != "reference_fasta" || source_state != "complete" || !source_validated {
        return Err(format!("{context}.source_fasta no es completo y validado"));
    }
    require_catalog_identity(source_id, source_fasta, artifacts, &context)?;
    let source_reference = references.get(source_reference_id).ok_or_else(|| {
        format!("{context} enlaza una referencia origen ausente: {source_reference_id}")
    })?;
    require_reference_fasta(source_reference, source_fasta, &context)?;
    require_non_reduced_source(source_reference, source_reference_id, &context)?;
    let source_catalog = reference_catalogs.get(source_reference_id).ok_or_else(|| {
        format!("{context} no tiene catálogo de secuencias origen: {source_reference_id}")
    })?;

    let (derived_fasta_id, derived_kind, derived_state, derived_validated) =
        artifact_fields(derived_fasta, &format!("{context}.derived_fasta"))?;
    let (map_id, map_kind, map_state, map_validated) =
        artifact_fields(map, &format!("{context}.map"))?;
    if derived_fasta_id != expected_fasta_id || derived_kind != "reference_fasta" {
        return Err(format!(
            "{context}.derived_fasta tiene un ID o clase ajenos"
        ));
    }
    if map_id != expected_map_id || map_kind != "reference_reduction_map" {
        return Err(format!("{context}.map tiene un ID o clase ajenos"));
    }
    require_artifact_state(&context, derived_state, derived_validated, unit_state)?;
    require_artifact_state(&context, map_state, map_validated, unit_state)?;
    require_catalog_identity(derived_fasta_id, derived_fasta, artifacts, &context)?;
    require_catalog_identity(map_id, map, artifacts, &context)?;
    claim(derived_fasta_id, claimed_artifacts)?;
    claim(map_id, claimed_artifacts)?;

    match unit_state {
        "complete" => {
            let derived_reference = references.get(derived_id).ok_or_else(|| {
                format!("reducción completa sin referencia derivada: {derived_id}")
            })?;
            let derived_catalog = reference_catalogs
                .get(derived_id)
                .ok_or_else(|| format!("reducción completa sin catálogo derivado: {derived_id}"))?;
            require_reference_fasta(derived_reference, derived_fasta, &context)?;
            require_reduced_origin(derived_reference, derived_id, source_reference_id, &context)?;
            validate_reduction_metadata(
                reduction,
                source_catalog,
                Some(derived_catalog),
                &context,
            )?;
            if !claimed_references.insert(derived_id.to_owned()) {
                return Err(format!(
                    "referencia derivada perteneciente a más de una reducción: {derived_id}"
                ));
            }
        }
        "failed" if references.contains_key(derived_id) => {
            return Err(format!(
                "reducción fallida presentada como referencia completa: {derived_id}"
            ));
        }
        "failed" => validate_reduction_metadata(reduction, source_catalog, None, &context)?,
        _ => return Err(format!("estado de reducción no admitido: {unit_state}")),
    }
    Ok(())
}

fn require_reduced_origin(
    reference: &Value,
    derived_id: &str,
    source_reference_id: &str,
    context: &str,
) -> Result<(), String> {
    let reference = object(reference, context)?;
    if text(reference, "id", context)? != derived_id {
        return Err(format!("{context} contradice derived_reference_id"));
    }
    let origin = reference
        .get("origin")
        .ok_or_else(|| format!("falta {context}.origin"))?;
    let origin = object(origin, &format!("{context}.origin"))?;
    if text(origin, "kind", &format!("{context}.origin"))? != "reduced"
        || text(origin, "source_reference_id", &format!("{context}.origin"))? != source_reference_id
    {
        return Err(format!(
            "{context}.origin no enlaza exactamente la referencia origen"
        ));
    }
    Ok(())
}

fn require_non_reduced_source(
    reference: &Value,
    source_reference_id: &str,
    context: &str,
) -> Result<(), String> {
    let reference = object(reference, context)?;
    if text(reference, "id", context)? != source_reference_id {
        return Err(format!("{context} contradice source_reference_id"));
    }
    let origin = reference
        .get("origin")
        .ok_or_else(|| format!("falta {context}.origin"))?;
    let origin = object(origin, &format!("{context}.origin"))?;
    if text(origin, "kind", &format!("{context}.origin"))? == "reduced" {
        return Err(format!(
            "{context} intenta encadenar una reducción no admitida en F1"
        ));
    }
    Ok(())
}

fn validate_partial_failed_reduction(
    derived_id: &str,
    unit_state: &str,
    expected_fasta_id: &str,
    expected_map_id: &str,
    artifacts: &BTreeMap<String, &Value>,
    references: &BTreeMap<String, &Value>,
    claimed_artifacts: &mut BTreeSet<String>,
) -> Result<(), String> {
    if unit_state != "failed" {
        return Err(format!(
            "reducción completa sin registro de metadatos: {derived_id}"
        ));
    }
    if references.contains_key(derived_id) {
        return Err(format!(
            "reducción parcial fallida presentada como referencia: {derived_id}"
        ));
    }
    let fasta = artifacts.get(expected_fasta_id);
    let map = artifacts.get(expected_map_id);
    if map.is_some() && fasta.is_none() {
        return Err(format!(
            "reducción parcial {derived_id} cataloga mapa sin FASTA"
        ));
    }
    for (expected_id, expected_kind, artifact) in [
        (expected_fasta_id, "reference_fasta", fasta),
        (expected_map_id, "reference_reduction_map", map),
    ] {
        if let Some(artifact) = artifact {
            let (id, kind, state, validated) =
                artifact_fields(artifact, &format!("partial reduction {derived_id}"))?;
            if id != expected_id || kind != expected_kind || state != "failed" || validated {
                return Err(format!(
                    "artefacto parcial incoherente para reducción {derived_id}: {id}"
                ));
            }
            claim(id, claimed_artifacts)?;
        }
    }
    Ok(())
}

fn validate_reduction_unit_state(
    id: &str,
    state: &str,
    exportable: bool,
    diagnostic: &Value,
    manifest_state: &str,
) -> Result<(), String> {
    match state {
        "complete" if exportable && diagnostic.is_null() => {
            if manifest_state == "complete" || manifest_state == "failed" {
                Ok(())
            } else {
                Err(format!("estado de manifiesto inválido para reducción {id}"))
            }
        }
        "failed" if !exportable => {
            diagnostic
                .as_str()
                .ok_or_else(|| format!("reducción fallida {id} sin diagnóstico"))?;
            Ok(())
        }
        _ => Err(format!(
            "estado/exportable/diagnóstico incoherentes para reducción {id}"
        )),
    }
}

fn require_artifact_state(
    context: &str,
    state: &str,
    validated: bool,
    unit_state: &str,
) -> Result<(), String> {
    match unit_state {
        "complete" if state == "complete" && validated => Ok(()),
        "failed" if state == "failed" && !validated => Ok(()),
        _ => Err(format!(
            "{context} contradice el estado {unit_state} de su unidad"
        )),
    }
}

fn require_catalog_identity(
    id: &str,
    embedded: &Value,
    artifacts: &BTreeMap<String, &Value>,
    context: &str,
) -> Result<(), String> {
    let cataloged = artifacts
        .get(id)
        .ok_or_else(|| format!("{context} enlaza un artefacto no catalogado: {id}"))?;
    if *cataloged != embedded {
        return Err(format!(
            "{context} y catálogo discrepan para el artefacto {id}"
        ));
    }
    Ok(())
}

fn require_reference_fasta(
    reference: &Value,
    expected_fasta: &Value,
    context: &str,
) -> Result<(), String> {
    let reference = object(reference, context)?;
    let fasta = reference
        .get("fasta")
        .ok_or_else(|| format!("falta {context}.fasta"))?;
    if fasta != expected_fasta {
        return Err(format!("{context} y referencia discrepan en FASTA"));
    }
    Ok(())
}

fn claim(id: &str, claimed: &mut BTreeSet<String>) -> Result<(), String> {
    if claimed.insert(id.to_owned()) {
        Ok(())
    } else {
        Err(format!(
            "artefacto de reducción perteneciente a más de una unidad: {id}"
        ))
    }
}

fn artifact_fields<'a>(
    value: &'a Value,
    context: &str,
) -> Result<(&'a str, &'a str, &'a str, bool), String> {
    let artifact = object(value, context)?;
    Ok((
        text(artifact, "id", context)?,
        text(artifact, "kind", context)?,
        text(artifact, "state", context)?,
        boolean(artifact, "validated", context)?,
    ))
}

fn object<'a>(
    value: &'a Value,
    context: &str,
) -> Result<&'a serde_json::Map<String, Value>, String> {
    value
        .as_object()
        .ok_or_else(|| format!("{context} no es un objeto"))
}

fn array<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    context: &str,
) -> Result<&'a [Value], String> {
    object
        .get(field)
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .ok_or_else(|| format!("falta o no es lista {context}.{field}"))
}

fn text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    context: &str,
) -> Result<&'a str, String> {
    object
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("falta o no es texto {context}.{field}"))
}

fn unsigned(
    object: &serde_json::Map<String, Value>,
    field: &str,
    context: &str,
) -> Result<u64, String> {
    object
        .get(field)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("falta o no es entero sin signo {context}.{field}"))
}

fn boolean(
    object: &serde_json::Map<String, Value>,
    field: &str,
    context: &str,
) -> Result<bool, String> {
    object
        .get(field)
        .and_then(Value::as_bool)
        .ok_or_else(|| format!("falta o no es booleano {context}.{field}"))
}

fn semantic_manifest_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path.display().to_string())
    .with(
        "schema",
        format!("manifest-general-v1.schema.json + {GENERAL_MANIFEST_RUNTIME_VALIDATION_VERSION}"),
    )
    .with("detail", detail.into())
}

fn effective_value_json(value: &EffectiveValue) -> Value {
    json!({
        "object_id": value.object_id(),
        "parameter": value.parameter().as_str(),
        "resolution_context": resolution_context_json(value.resolution_context()),
        "value": rule_value_json(value.value()),
        "winning_rule_id": value.winning_rule_id(),
        "selectors": selectors_json(value.selectors()),
        "provenance": {"path": value.provenance().path, "location": value.provenance().location},
    })
}

fn resolution_context_json(context: &ResolutionContext) -> Value {
    json!({
        "individual": context.individual,
        "population": context.population,
        "subpopulation": context.subpopulation,
        "phase": context.phase,
        "generation": context.generation,
        "reference": context.reference,
        "sequence": context.sequence,
        "position": context.position,
        "locus": context.locus,
        "group": context.group,
        "block": context.block,
        "technical_event": context.technical_event,
    })
}

/// Documento autoritativo previo a la configuración efectiva. Conserva el
/// dominio validado y el conjunto exacto de solicitudes que el productor
/// resolvió; no contiene valores ganadores derivados.
#[must_use]
pub fn resolution_evidence_document(config: &EffectiveConfig) -> Value {
    json!({
        "resolution_evidence_kind": "gps_resolution_evidence",
        "schema_version": "1",
        "parameter_order_version": crate::rules::PARAMETER_ORDER_VERSION,
        "resolution_domain": resolution_domain_json(config.resolution_domain()),
        "resolution_requests": config.resolution_requests().iter().map(|request| {
            json!({
                "object_id": request.object_id,
                "parameter": request.parameter.as_str(),
                "resolution_context": resolution_context_json(&request.context),
            })
        }).collect::<Vec<_>>(),
    })
}

fn resolution_domain_json(domain: &ResolutionDomain) -> Value {
    let populations = domain
        .population_memberships()
        .iter()
        .map(|(subpopulation, population)| {
            json!({"population": population, "subpopulation": subpopulation})
        })
        .collect::<Vec<_>>();
    let phases = domain
        .phase_generation_intervals()
        .iter()
        .map(|(phase, (start, end))| json!({"phase": phase, "start": start, "end": end}))
        .collect::<Vec<_>>();
    let loci = domain
        .loci()
        .values()
        .map(|locus| {
            json!({
                "id": locus.id,
                "reference": locus.reference,
                "sequence": locus.sequence,
                "position": locus.position,
            })
        })
        .collect::<Vec<_>>();
    json!({
        "class": "resolution_domain",
        "version": "1",
        "provider": domain.provider(),
        "population_memberships": populations,
        "phase_generation_intervals": phases,
        "loci": loci,
    })
}

fn parse_resolution_domain(value: &Value) -> Result<ResolutionDomain, String> {
    let domain = object(value, "effective.resolution_domain")?;
    if text(domain, "class", "effective.resolution_domain")? != "resolution_domain"
        || text(domain, "version", "effective.resolution_domain")? != "1"
    {
        return Err("effective.resolution_domain usa clase o versión incompatible".to_owned());
    }
    let provider = text(domain, "provider", "effective.resolution_domain")?;
    let populations = array(
        domain,
        "population_memberships",
        "effective.resolution_domain",
    )?
    .iter()
    .enumerate()
    .map(|(index, value)| {
        let context = format!("effective.resolution_domain.population_memberships[{index}]");
        let entry = object(value, &context)?;
        Ok(PopulationMembership {
            population: text(entry, "population", &context)?.to_owned(),
            subpopulation: text(entry, "subpopulation", &context)?.to_owned(),
        })
    })
    .collect::<Result<Vec<_>, String>>()?;
    let phases = array(
        domain,
        "phase_generation_intervals",
        "effective.resolution_domain",
    )?
    .iter()
    .enumerate()
    .map(|(index, value)| {
        let context = format!("effective.resolution_domain.phase_generation_intervals[{index}]");
        let entry = object(value, &context)?;
        Ok(PhaseGenerationInterval {
            phase: text(entry, "phase", &context)?.to_owned(),
            start: u32::try_from(unsigned(entry, "start", &context)?)
                .map_err(|_| format!("{context}.start no cabe en u32"))?,
            end: u32::try_from(unsigned(entry, "end", &context)?)
                .map_err(|_| format!("{context}.end no cabe en u32"))?,
        })
    })
    .collect::<Result<Vec<_>, String>>()?;
    let loci = array(domain, "loci", "effective.resolution_domain")?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let context = format!("effective.resolution_domain.loci[{index}]");
            let entry = object(value, &context)?;
            Ok(DomainLocus {
                id: text(entry, "id", &context)?.to_owned(),
                reference: text(entry, "reference", &context)?.to_owned(),
                sequence: text(entry, "sequence", &context)?.to_owned(),
                position: unsigned(entry, "position", &context)?,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    ResolutionDomain::validated(provider, populations, phases, loci)
}

fn parse_resolution_evidence(
    value: &Value,
) -> Result<(ResolutionDomain, Vec<ResolutionRequest>), String> {
    let evidence = object(value, "resolution_evidence")?;
    if text(evidence, "resolution_evidence_kind", "resolution_evidence")?
        != "gps_resolution_evidence"
        || text(evidence, "schema_version", "resolution_evidence")? != "1"
        || text(evidence, "parameter_order_version", "resolution_evidence")? != "1"
    {
        return Err("resolution_evidence usa clase o versión incompatible".to_owned());
    }
    let domain = parse_resolution_domain(
        evidence
            .get("resolution_domain")
            .ok_or_else(|| "falta resolution_evidence.resolution_domain".to_owned())?,
    )?;
    let mut requests = Vec::new();
    let mut previous: Option<(&str, usize)> = None;
    for (index, request) in array(evidence, "resolution_requests", "resolution_evidence")?
        .iter()
        .enumerate()
    {
        let context = format!("resolution_evidence.resolution_requests[{index}]");
        let request = object(request, &context)?;
        let object_id = text(request, "object_id", &context)?;
        let parameter = ParameterKey::parse(text(request, "parameter", &context)?)
            .ok_or_else(|| format!("{context}.parameter no pertenece al catálogo"))?;
        let identity = (object_id, parameter.canonical_rank());
        if previous.is_some_and(|prior| prior >= identity) {
            return Err(
                "resolution_evidence.resolution_requests no está ordenado o repite identidad"
                    .to_owned(),
            );
        }
        previous = Some(identity);
        requests.push(ResolutionRequest {
            object_id: object_id.to_owned(),
            parameter,
            context: parse_resolution_context(
                request
                    .get("resolution_context")
                    .ok_or_else(|| format!("falta {context}.resolution_context"))?,
                &context,
            )?,
        });
    }
    Ok((domain, requests))
}

fn rule_value_json(value: &RuleValue) -> Value {
    match value {
        RuleValue::Float(value) => json!(value),
        RuleValue::Integer(value) => json!(value),
        RuleValue::Boolean(value) => json!(value),
        RuleValue::Text(value) => json!(value),
    }
}

fn selectors_json(selectors: &SelectorSet) -> Value {
    json!({
        "individual": selectors.individual,
        "population": selectors.population.as_ref().map(|selector| match selector {
            PopulationSelector::Population(id) => json!({"kind": "population", "id": id}),
            PopulationSelector::Subpopulation(id) => json!({"kind": "subpopulation", "id": id}),
        }),
        "temporal": selectors.temporal.as_ref().map(|selector| match selector {
            TemporalSelector::Phase(id) => json!({"kind": "phase", "id": id}),
            TemporalSelector::Generation(generation) => json!({"kind": "generation", "generation": generation}),
        }),
        "spatial": selectors.spatial.as_ref().map(|selector| match selector {
            SpatialSelector::Sequence { reference, sequence } => json!({"kind": "sequence", "reference": reference, "sequence": sequence}),
            SpatialSelector::Region(region) => json!({"kind": "region", "reference": region.reference, "sequence": region.sequence, "start": region.start, "end": region.end}),
            SpatialSelector::Locus { reference, sequence, locus } => json!({"kind": "locus", "reference": reference, "sequence": sequence, "locus": locus}),
        }),
        "group": selectors.group,
        "block": selectors.block,
        "technical_event": selectors.technical_event,
    })
}

fn event_subseed_json(event: &EventSubseed) -> Value {
    json!({"event_id": event.event_id(), "subseed": event.subseed().to_hex()})
}

fn normalization_json(normalization: &NormalizedCollection) -> Value {
    let changed = normalization
        .values
        .iter()
        .any(|value| value.before.to_bits() != value.after.to_bits());
    json!({
        "parameter": normalization.parameter,
        "tolerance": normalization.tolerance,
        "canonical_order": normalization.canonical_order,
        "algorithm": normalization.algorithm,
        "algorithm_version": normalization.algorithm_version,
        "sum_before": normalization.sum_before,
        "delta_before": normalization.delta_before,
        "decision": if changed { "normalized" } else { "accepted_without_change" },
        "values": normalization.values.iter().map(|value| json!({"id": value.id, "category": value.category, "before": value.before, "after": value.after})).collect::<Vec<_>>(),
        "sum_after": normalization.sum_after,
        "delta_after": normalization.delta_after,
        "specific_tolerances": [],
    })
}

fn reproducibility_levels() -> Vec<Value> {
    vec![
        json!({"level": "biological", "status": "guaranteed", "conditions": ["simulator_version", "platform", "root_seed", "effective_configuration", "input_contents", "profiles", "random_algorithm", "subseed_catalog_v3", "reference_fasta_sha256"]}),
        json!({"level": "technical", "status": "not_applicable", "conditions": []}),
        json!({"level": "semantic", "status": "not_guaranteed", "conditions": ["comparison_required"]}),
        json!({"level": "binary", "status": "not_guaranteed", "conditions": ["format_specific_evidence_required"]}),
    ]
}

fn reference_json(config: &EffectiveConfig, output: ReferenceOutput<'_>) -> Value {
    let origin = match output.reference.origin() {
        ReferenceOrigin::Generated => json!({"kind": "generated"}),
        ReferenceOrigin::Supplied { original_path } => json!({
            "kind": "supplied",
            "original_path": original_path.to_string_lossy(),
        }),
        ReferenceOrigin::Reduced {
            source_reference_id,
        } => json!({
            "kind": "reduced",
            "source_reference_id": source_reference_id.as_str(),
        }),
    };
    let canonical_format = matches!(
        output.reference.origin(),
        ReferenceOrigin::Generated | ReferenceOrigin::Reduced { .. }
    )
    .then(|| {
        json!({
            "version": CANONICAL_FASTA_VERSION,
            "line_base_count": CANONICAL_LINE_BASE_COUNT,
            "line_ending": "LF",
            "base_case": "uppercase",
        })
    });
    let sequences = output
        .reference
        .sequences()
        .iter()
        .map(|sequence| {
            let ploidy = config
                .sequences()
                .iter()
                .find(|declared| declared.id() == sequence.id().as_str())
                .map_or(1, |declared| declared.ploidy());
            let header = reversible_bytes_json(sequence.header().raw_bytes());
            let description = sequence.header().description().map(|text| {
                json!({
                    "encoding": text.encoding(),
                    "value": text.encoded(),
                })
            });
            json!({
                "id": sequence.id().as_str(),
                "class": match sequence.class() {
                    ReferenceSequenceClass::Chromosome => "chromosome",
                    ReferenceSequenceClass::Scaffold => "scaffold",
                },
                "length": sequence.length(),
                "ploidy": ploidy,
                "header": header,
                "description": description,
                "alphabet_validation": match sequence.alphabet_validation() {
                    AlphabetValidation::AcgtnCaseInsensitive => "ACGTN-case-insensitive",
                },
                "unknown_base_count": sequence.unknown_base_count(),
            })
        })
        .collect::<Vec<_>>();
    json!({
        "id": output.reference.id().as_str(),
        "total_genome_size": output.reference.total_length(),
        "origin": origin,
        "authoritative_resource": output.reference.authoritative_resource().map(|path| path.to_string_lossy()),
        "fasta": output.fasta.to_json(),
        "canonical_format": canonical_format,
        "indexes": [],
        "sequences": sequences,
    })
}

fn reduction_json(config: &EffectiveConfig, output: ReducedReferenceOutput<'_>) -> Value {
    let mut included = output
        .reference
        .sequences()
        .iter()
        .map(|sequence| sequence.id().as_str())
        .collect::<Vec<_>>();
    included.sort_unstable();
    let mut omitted = config
        .sequences()
        .iter()
        .filter(|sequence| !included.contains(&sequence.id()))
        .map(|sequence| sequence.id())
        .collect::<Vec<_>>();
    omitted.sort_unstable();
    let mut selected = output.plan.segments().iter().collect::<Vec<_>>();
    selected.sort_by_key(|segment| segment.selection_order());
    let selection_order = selected
        .into_iter()
        .map(|segment| {
            json!({
                "selection_order": segment.selection_order(),
                "sequence_id": segment.source_interval().sequence_id().as_str(),
                "start": segment.source_interval().start() + 1,
                "end": segment.source_interval().end(),
            })
        })
        .collect::<Vec<_>>();
    json!({
        "source_reference_id": output.effective.source_reference_id(),
        "derived_reference_id": output.effective.id(),
        "source_fasta": output.source_fasta.to_json(),
        "derived_fasta": artifact_json_with_state(output.fasta, output.state),
        "eligible_bases": config.total_genome_size(),
        "region_size_bases": output.effective.region_size_bases(),
        "target": {
            "mode": output.effective.target().mode(),
            "declared_value": output.effective.target().declared_value(),
            "derived_bases": output.effective.target_bases(),
        },
        "randomness": {
            "root_seed": config.reproducibility().root_seed().to_hex(),
            "subseed": output.effective.selection_subseed().to_hex(),
        },
        "algorithm": REDUCED_REFERENCE_SELECTION_ALGORITHM,
        "map": artifact_json_with_state(output.map, output.state),
        "included_sequence_ids": included,
        "omitted_sequence_ids": omitted,
        "has_residual": output.plan.segments().iter().any(|segment| segment.class() == ReductionSegmentClass::Residual),
        "selection_order": selection_order,
    })
}

fn artifact_json_with_state(artifact: &ArtifactRecord, state: UnitState) -> Value {
    let mut document = artifact.to_json();
    if state == UnitState::Failed {
        document["state"] = json!("failed");
        document["validated"] = json!(false);
    }
    document
}

fn validate_reference_output(output: ReferenceOutput<'_>) -> Result<(), Diagnostic> {
    let authoritative = output.reference.authoritative_resource().ok_or_else(|| {
        Diagnostic::internal(
            "GPS-INTERNAL-001",
            "Una invariante interna considerada imposible se ha roto.",
        )
        .with("component", "results_json")
        .with("operation", "effective_config_document_with_reference")
        .with("detail", "reference has no authoritative resource")
    })?;
    if authoritative != Path::new(&output.fasta.path) {
        return Err(Diagnostic::internal(
            "GPS-INTERNAL-001",
            "Una invariante interna considerada imposible se ha roto.",
        )
        .with("component", "results_json")
        .with("operation", "effective_config_document_with_reference")
        .with(
            "detail",
            "reference authoritative resource differs from FASTA artifact",
        ));
    }
    Ok(())
}

fn reference_link_json(output: ReferenceOutput<'_>) -> Value {
    json!({
        "id": output.reference.id().as_str(),
        "fasta": output.fasta.to_json(),
        "indexes": [],
    })
}

fn reversible_bytes_json(bytes: &[u8]) -> Value {
    match std::str::from_utf8(bytes) {
        Ok(text) => json!({"encoding": "utf-8", "value": text}),
        Err(_) => json!({"encoding": "bytes-hex", "value": encode_hex(bytes)}),
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

fn unit_json(unit: &UnitRecord) -> Value {
    json!({
        "id": unit.id,
        "kind": unit.kind,
        "state": unit.state.as_str(),
        "exportable": unit.exportable,
        "diagnostic_code": unit.diagnostic_code,
    })
}

fn diagnostic_json(diagnostic: &Diagnostic) -> Value {
    json!({
        "code": diagnostic.code(),
        "category": diagnostic.category(),
        "message": diagnostic.message(),
        "context": diagnostic.context(),
    })
}

fn validate_relative_link(path: &str) -> Result<(), Diagnostic> {
    let path = Path::new(path);
    if path.is_absolute()
        || path
            .components()
            .any(|component| matches!(component, std::path::Component::ParentDir))
    {
        return Err(Diagnostic::io(
            "GPS-IO-005",
            "Una unidad incompleta o inválida no puede publicarse como completa.",
        )
        .with("unit", path.display().to_string())
        .with("state", "invalid_link")
        .with(
            "detail",
            "artifact links must remain under the execution root",
        ));
    }
    Ok(())
}

fn json_error(error: serde_json::Error) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "results_json")
    .with("operation", "serialize")
    .with("detail", error.to_string())
}

fn io_error(path: &Path, operation: &str, error: std::io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}
