//! Contratos F4.3.1–F4.3.2: el único perfil es PCR-free, GC neutro y sin
//! adaptadores. Los rechazos se prueban en la frontera TOML, antes de resolver
//! o generar cualquier pareja.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::results::{
    ArtifactRecord, ReferenceOutput, effective_config_document_with_reference,
    validate_general_manifest_bytes,
};
use serde_json::json;

use genome_population_simulator::results::validate_effective_config_against_manifest_bytes;

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn initial_profile_has_only_neutral_non_configurable_artifacts() {
    let effective = build_effective_config(&declared(valid_event())).expect("effective config");
    let event = effective.technical_event().expect("technical event");

    assert!(!event.pcr_duplicates_applicable());
    assert_eq!(event.realized_optical_duplicate_pairs(), 0);
    assert!(event.read_length() < 250);

    let mut reference = Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("sequence ID"),
                FastaHeader::parse(b"chr1").expect("FASTA header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(&[b'A'; 100]).expect("DNA"),
            )
            .expect("reference sequence"),
        ],
    )
    .expect("reference");
    let fasta = ArtifactRecord {
        id: "reference-main-fasta".to_owned(),
        kind: "reference_fasta".to_owned(),
        path: "references/reference-main.fasta".to_owned(),
        size_bytes: 100,
        digest: "a".repeat(64),
    };
    reference
        .bind_authoritative_resource(fasta.path.clone().into())
        .expect("authoritative resource");
    let document = effective_config_document_with_reference(
        &effective,
        ReferenceOutput {
            reference: &reference,
            fasta: &fasta,
        },
    )
    .expect("effective document");
    assert_eq!(
        document["technical_event"]["profile"]["insert_length"],
        json!({"minimum_bases": 250, "maximum_bases": 450})
    );
    assert_eq!(
        document["technical_event"]["artifacts"],
        json!({
            "chimera": {"rate": 0.0},
            "index_hopping": {"rate": 0.0},
            "pcr": {"applicable": false},
            "optical_duplicates": {"realized_pairs": 0},
            "gc": {"mode": "neutral"},
            "adapters": {"enabled": false},
        })
    );
}

#[test]
fn unsupported_profile_artifact_configuration_is_rejected_before_effective_resolution() {
    for (field, declaration) in [
        ("profile", "profile = \"amplified-wgs-v1\""),
        ("pcr_duplicate_rate", "pcr_duplicate_rate = 0.1"),
        ("optical_duplicate_rate", "optical_duplicate_rate = 0.1"),
        ("gc_bias_intensity", "gc_bias_intensity = 0.1"),
        ("gc_bias_curve", "gc_bias_curve = \"calibrated-v1\""),
        ("adapters_enabled", "adapters_enabled = true"),
    ] {
        let error = parse_declared_config(admitted(&format!("{}\n{declaration}\n", valid_event())))
            .expect_err("F4.3 must reject unsupported profile artifact configuration");
        assert_eq!(error.code(), "GPS-CONF-006", "{field}");
        assert_eq!(
            error.context().get("field"),
            Some(&field.to_owned()),
            "{field}"
        );
    }
}

#[test]
fn observable_pe_categories_and_direct_alignment_are_rejected_before_publication() {
    for (field, declaration, deferred_to) in [
        (
            "orientation_anomaly",
            "orientation_anomaly = true",
            "DTF-055",
        ),
        (
            "orientation_anomaly_rate",
            "orientation_anomaly_rate = 0.1",
            "DTF-055",
        ),
        (
            "observed_insert_size",
            "observed_insert_size = 300",
            "DTF-055",
        ),
        (
            "observed_insert_length",
            "observed_insert_length = 300",
            "DTF-055",
        ),
        ("unaligned_mate", "unaligned_mate = true", "DTF-055"),
        (
            "unaligned_mate_rate",
            "unaligned_mate_rate = 0.1",
            "DTF-055",
        ),
        (
            "paired_end_observables",
            "[technical_event.paired_end_observables]",
            "DTF-055",
        ),
        (
            "orientation_anomalies",
            "[technical_event.orientation_anomalies]",
            "DTF-055",
        ),
        (
            "observed_insert_sizes",
            "[technical_event.observed_insert_sizes]",
            "DTF-055",
        ),
        (
            "observed_insert_lengths",
            "[technical_event.observed_insert_lengths]",
            "DTF-055",
        ),
        (
            "unaligned_mate",
            "[technical_event.unaligned_mate]",
            "DTF-055",
        ),
        (
            "unaligned_mates",
            "[technical_event.unaligned_mates]",
            "DTF-055",
        ),
        ("direct_alignment", "direct_alignment = true", "DTF-051"),
        (
            "direct_alignment",
            "[technical_event.direct_alignment]",
            "DTF-051",
        ),
        (
            "alignment_category",
            "alignment_category = \"concordant\"",
            "DTF-051",
        ),
        (
            "alignment_categories",
            "[technical_event.alignment_categories]",
            "DTF-051",
        ),
        (
            "direct_alignment_categories",
            "[technical_event.direct_alignment_categories]",
            "DTF-051",
        ),
    ] {
        let source = format!("{}\n{declaration}\n", valid_event());
        let error = parse_declared_config(admitted(&source))
            .expect_err("F4 must reject deferred PE and direct-alignment requests");
        assert_eq!(error.code(), "GPS-CONF-023", "{field}");
        assert_eq!(error.context().get("field"), Some(&field.to_owned()));
        assert_eq!(
            error.context().get("deferred_to"),
            Some(&deferred_to.to_owned())
        );

        let root = temporary("deferred-alignment-request");
        fs::write(&root, source).expect("invalid configuration");
        let output = root.with_extension("output");
        let run = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
            .args(["run"])
            .arg(&root)
            .arg("--output")
            .arg(&output)
            .output()
            .expect("rejected run");
        assert!(
            !run.status.success(),
            "{field} must fail before publication"
        );
        assert!(
            !output.exists(),
            "{field} must not create an output directory before rejection"
        );
        fs::remove_file(&root).expect("invalid configuration cleanup");
    }
}

#[test]
fn published_f4_effective_config_rejects_mutated_duplicate_artifacts() {
    let root = temporary("published-f4-duplicate-artifacts");
    if root.exists() {
        fs::remove_dir_all(&root).expect("stale test directory");
    }
    fs::create_dir_all(&root).expect("test directory");
    let config_path = root.join("general.toml");
    fs::write(
        &config_path,
        valid_published_base().replacen(
            "schema_version = \"2\"\n\n",
            "schema_version = \"2\"\n\n[population_snvs]\nmode = \"generate\"\n[population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 1\nspatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n[[population_snvs.generate.populations]]\npopulation_id = \"P1\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n\n",
            1,
        ),
    )
    .expect("configuration");
    let output = root.join("output");
    let run = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run"])
        .arg(&config_path)
        .arg("--output")
        .arg(&output)
        .output()
        .expect("published run");
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );

    let effective_path = output.join("effective-config.json");
    let manifest_path = output.join("manifest.json");
    let effective_bytes = fs::read(&effective_path).expect("effective configuration");
    let manifest_bytes = fs::read(&manifest_path).expect("manifest");
    let mut effective: serde_json::Value =
        serde_json::from_slice(&effective_bytes).expect("effective JSON");
    effective["technical_event"] = f4_technical_event_document(
        effective["seed"]["root"]
            .as_str()
            .expect("published root seed"),
    );
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&manifest_bytes).expect("manifest JSON");
    add_f4_technical_event_link(&mut manifest);
    let manifest_bytes = serde_json::to_vec(&manifest).expect("manifest JSON");
    let effective_bytes = serde_json::to_vec(&effective).expect("effective JSON");
    validate_effective_config_against_manifest_bytes(
        &effective_path,
        &effective_bytes,
        &manifest_path,
        &manifest_bytes,
    )
    .expect("valid published F4 effective configuration");
    validate_general_manifest_bytes(&manifest_path, &manifest_bytes)
        .expect("valid general F4 chimera summary");
    let mut zero_rate_with_realized_pairs = manifest.clone();
    zero_rate_with_realized_pairs["technical_events"][0]["artifacts"]["chimera"]["realized_pairs"] =
        json!(1);
    assert!(
        validate_general_manifest_bytes(
            &manifest_path,
            &serde_json::to_vec(&zero_rate_with_realized_pairs).expect("manifest JSON"),
        )
        .is_err(),
        "zero chimera rate must prohibit realized pairs in the general manifest"
    );
    for (name, mutation) in [
        ("pcr-applicable", DuplicateArtifactMutation::PcrApplicable),
        (
            "optical-realized-pairs",
            DuplicateArtifactMutation::OpticalRealizedPairs,
        ),
        (
            "artifacts-removed",
            DuplicateArtifactMutation::ArtifactsRemoved,
        ),
        ("extra-rate", DuplicateArtifactMutation::ExtraRate),
        ("gc-not-neutral", DuplicateArtifactMutation::GcNotNeutral),
        (
            "adapters-enabled",
            DuplicateArtifactMutation::AdaptersEnabled,
        ),
        (
            "insert-not-greater-than-read",
            DuplicateArtifactMutation::InsertNotGreaterThanRead,
        ),
    ] {
        let mut mutated = effective.clone();
        mutation.apply(&mut mutated);
        assert!(
            validate_effective_config_against_manifest_bytes(
                &effective_path,
                &serde_json::to_vec(&mutated).expect("effective JSON"),
                &manifest_path,
                &manifest_bytes,
            )
            .is_err(),
            "{name} must be rejected",
        );
    }

    fs::remove_dir_all(root).expect("cleanup");
}

fn add_f4_technical_event_link(manifest: &mut serde_json::Value) {
    let event_manifest = json!({
        "id": "event-reads-main-manifest",
        "kind": "event_manifest",
        "path": "events/reads-main/manifest.json",
        "size_bytes": 1,
        "checksum": {"algorithm": "sha256", "digest": "0".repeat(64)},
    });
    manifest["artifacts"]
        .as_array_mut()
        .expect("artifact catalog")
        .push(event_manifest.clone());
    manifest["output_categories"]["minimum"]
        .as_array_mut()
        .expect("minimum output category")
        .push(json!("event-reads-main-manifest"));
    let effective_config = manifest["artifacts"]
        .as_array()
        .expect("artifact catalog")
        .iter()
        .find(|artifact| artifact["id"] == manifest["effective_config_artifact_id"])
        .expect("effective config artifact")
        .clone();
    manifest["events"] = json!([{
        "event_id": "reads-main",
        "state": "complete",
        "manifest_artifact_id": "event-reads-main-manifest",
    }]);
    manifest["technical_events"] = json!([{
        "event_id": "reads-main",
        "effective_config": effective_config,
        "profile": {
            "id": "gps-hiseqx-pcrfree-pe150-v1",
            "sha256": "3d8878ab920237dae03705276a26dbb0823941830d531374d35aa1c540f6deea",
        },
        "rng_catalog": "f4-technical-v3",
        "id_catalog": {"version":"3", "derivations":{"v1":["execution","reference","individual","event","batch","library","read_group","template"],"v2":["PreHopInstanceId","TechnicalPairId"]}},
        "artifacts": {"chimera": {"rate": 0.0, "realized_pairs": 0}},
        "event_manifest": event_manifest,
    }]);
}

fn f4_technical_event_document(root_seed: &str) -> serde_json::Value {
    let source = valid_event().replacen(
        "schema_version = \"2\"\n",
        &format!("schema_version = \"2\"\nseed = \"{root_seed}\"\n"),
        1,
    );
    let effective = build_effective_config(&declared(&source)).expect("effective config");
    let mut reference = Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("sequence ID"),
                FastaHeader::parse(b"chr1").expect("FASTA header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(&[b'A'; 100]).expect("DNA"),
            )
            .expect("reference sequence"),
        ],
    )
    .expect("reference");
    let fasta = ArtifactRecord {
        id: "reference-main-fasta".to_owned(),
        kind: "reference_fasta".to_owned(),
        path: "references/reference-main.fasta".to_owned(),
        size_bytes: 100,
        digest: "a".repeat(64),
    };
    reference
        .bind_authoritative_resource(fasta.path.clone().into())
        .expect("authoritative resource");
    effective_config_document_with_reference(
        &effective,
        ReferenceOutput {
            reference: &reference,
            fasta: &fasta,
        },
    )
    .expect("effective document")["technical_event"]
        .clone()
}

#[derive(Clone, Copy)]
enum DuplicateArtifactMutation {
    PcrApplicable,
    OpticalRealizedPairs,
    ArtifactsRemoved,
    ExtraRate,
    GcNotNeutral,
    AdaptersEnabled,
    InsertNotGreaterThanRead,
}

impl DuplicateArtifactMutation {
    fn apply(self, document: &mut serde_json::Value) {
        match self {
            Self::PcrApplicable => {
                document["technical_event"]["artifacts"]["pcr"]["applicable"] = json!(true);
            }
            Self::OpticalRealizedPairs => {
                document["technical_event"]["artifacts"]["optical_duplicates"]["realized_pairs"] =
                    json!(1);
            }
            Self::ArtifactsRemoved => {
                document["technical_event"]
                    .as_object_mut()
                    .expect("technical event")
                    .remove("artifacts");
            }
            Self::ExtraRate => {
                document["technical_event"]["artifacts"]["pcr"]["rate"] = json!(0.1);
            }
            Self::GcNotNeutral => {
                document["technical_event"]["artifacts"]["gc"]["mode"] = json!("biased");
            }
            Self::AdaptersEnabled => {
                document["technical_event"]["artifacts"]["adapters"]["enabled"] = json!(true);
            }
            Self::InsertNotGreaterThanRead => {
                document["technical_event"]["profile"]["insert_length"]["minimum_bases"] =
                    json!(150);
            }
        }
    }
}

fn declared(source: &str) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    parse_declared_config(admitted(source)).expect("declared config")
}

fn admitted(source: &str) -> genome_population_simulator::configuration::AdmittedToml {
    admit_toml_header(
        Path::new("f4-profile-rejections.toml"),
        source.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("TOML header")
}

fn valid_event() -> &'static str {
    r#"config_kind = "general"
schema_version = "2"

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
length = 100

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 1

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
}

fn valid_published_base() -> &'static str {
    r#"config_kind = "general"
schema_version = "2"

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
length = 100

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A"]
"#
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-profile-rejections-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
