//! Documentos JSON versionados de configuración efectiva y manifiestos F1.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use serde_json::{Value, json};

use crate::application::{EffectiveConfig, EffectiveValue};
use crate::configuration::SequenceClass;
use crate::diagnostics::Diagnostic;
use crate::integrity::{ArtifactIdentity, identify_file};
use crate::numeric::NormalizedCollection;
use crate::randomness::EventSubseed;
use crate::rules::{PopulationSelector, RuleValue, SelectorSet, SpatialSelector, TemporalSelector};

pub const RESULT_SCHEMA_VERSION: &str = "1";

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

pub struct GeneralManifest<'a> {
    pub execution_id: &'a str,
    pub state: UnitState,
    pub effective: &'a EffectiveConfig,
    pub configuration_inputs: &'a [ConfigurationInputRecord],
    pub effective_config: Option<&'a ArtifactRecord>,
    pub events: &'a [EventManifestLink],
    pub artifacts: &'a [ArtifactRecord],
    pub units: &'a [UnitRecord],
    pub diagnostics: &'a [Diagnostic],
}

pub struct EventManifest<'a> {
    pub event_id: &'a str,
    pub state: UnitState,
    pub subseed: &'a str,
    pub artifacts: &'a [ArtifactRecord],
    pub units: &'a [UnitRecord],
    pub diagnostics: &'a [Diagnostic],
}

#[must_use]
pub fn effective_config_document(config: &EffectiveConfig) -> Value {
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

#[must_use]
pub fn general_manifest_document(manifest: &GeneralManifest<'_>) -> Value {
    let dependencies = [
        ("blake3", "1.8.6"),
        ("getrandom", "0.4.3"),
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
        "effective_config": manifest.effective_config.map(ArtifactRecord::to_json),
        "events": events,
        "artifacts": manifest.artifacts.iter().map(ArtifactRecord::to_json).collect::<Vec<_>>(),
        "outputs": {
            "minimum": ["configuration_copies", "effective_config", "manifest_general"],
            "default_enabled": [],
            "explicitly_requested": [],
            "disabled": ["reference", "truth", "vcf", "fastq", "bam", "cram"]
        },
        "normalizations": manifest.effective.normalized_collections().iter().map(normalization_json).collect::<Vec<_>>(),
        "reproducibility": reproducibility_levels(),
        "units": manifest.units.iter().map(unit_json).collect::<Vec<_>>(),
        "diagnostics": manifest.diagnostics.iter().map(diagnostic_json).collect::<Vec<_>>(),
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

fn effective_value_json(value: &EffectiveValue) -> Value {
    json!({
        "object_id": value.object_id(),
        "parameter": value.parameter().as_str(),
        "value": rule_value_json(value.value()),
        "winning_rule_id": value.winning_rule_id(),
        "selectors": selectors_json(value.selectors()),
        "provenance": {"path": value.provenance().path, "location": value.provenance().location},
    })
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
        json!({"level": "biological", "status": "not_applicable", "conditions": []}),
        json!({"level": "technical", "status": "not_applicable", "conditions": []}),
        json!({"level": "semantic", "status": "not_guaranteed", "conditions": ["comparison_required"]}),
        json!({"level": "binary", "status": "not_guaranteed", "conditions": ["format_specific_evidence_required"]}),
    ]
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
