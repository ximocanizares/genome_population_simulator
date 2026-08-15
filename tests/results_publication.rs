use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;
use std::time::{Duration, Instant};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::publication::{
    ArtifactLifecycle, ArtifactStage, ExecutionWorkspace, PublishableUnit,
};
use genome_population_simulator::reference::SequenceClass as ReferenceSequenceClass;
use genome_population_simulator::reference_generation::{
    BaseFrequency, GeneratedSequenceSpec, generate_reference,
};
use genome_population_simulator::results::{
    ArtifactRecord, ReferenceOutput, UnitState, effective_config_document_with_reference,
    validate_general_manifest_bytes,
};

const BINARY: &str = env!("CARGO_BIN_EXE_genome-population-simulator");
static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f14-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn structurally_valid_manifest_bytes() -> &'static [u8] {
    &manifest_fixture().manifest
}

struct ManifestFixture {
    manifest: Vec<u8>,
    artifacts: Vec<(String, String, PathBuf, Vec<u8>)>,
}

fn manifest_fixture() -> &'static ManifestFixture {
    static FIXTURE: OnceLock<ManifestFixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let destination = temporary("structural-manifest-fixture");
        let effective = load_and_validate_config(&fixture("general-seeded.toml"))
            .expect("fixture configuration");
        execute_validated_config(&effective, &destination).expect("fixture execution");
        let manifest = fs::read(destination.join("manifest.json")).expect("fixture manifest");
        let document: serde_json::Value = serde_json::from_slice(&manifest).unwrap();
        let artifacts = document["artifacts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|artifact| {
                let id = artifact["id"].as_str().unwrap().to_owned();
                let kind = artifact["kind"].as_str().unwrap().to_owned();
                let path = PathBuf::from(artifact["path"].as_str().unwrap());
                let bytes = fs::read(destination.join(&path)).unwrap();
                (id, kind, path, bytes)
            })
            .collect();
        fs::remove_dir_all(destination).expect("fixture cleanup");
        ManifestFixture {
            manifest,
            artifacts,
        }
    })
}

fn matching_workspace(
    destination: &Path,
) -> (
    ExecutionWorkspace,
    genome_population_simulator::results::ValidatedGeneralManifest,
) {
    let mut workspace = ExecutionWorkspace::create(destination).expect("workspace");
    for (id, kind, path, bytes) in &manifest_fixture().artifacts {
        workspace
            .write_artifact(id, kind, path, bytes)
            .expect("fixture payload artifact");
    }
    let manifest = workspace
        .validate_general_manifest(
            Path::new("manifest.json"),
            structurally_valid_manifest_bytes(),
        )
        .expect("structurally, semantically and physically valid manifest fixture");
    (workspace, manifest)
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/config")
        .join(name)
}

fn provisional_siblings(destination: &Path) -> Vec<PathBuf> {
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    let destination_name = destination
        .file_name()
        .and_then(|name| name.to_str())
        .expect("destination name");
    let prefix = format!(".{destination_name}.gps-incomplete-");
    let mut matches = fs::read_dir(parent)
        .expect("provisional parent")
        .filter_map(Result::ok)
        .filter(|entry| entry.file_name().to_string_lossy().starts_with(&prefix))
        .map(|entry| entry.path())
        .collect::<Vec<_>>();
    matches.sort();
    matches
}

fn provisional_workspace(destination: &Path) -> PathBuf {
    let siblings = provisional_siblings(destination);
    assert_eq!(siblings.len(), 1, "exactly one provisional workspace");
    siblings[0].clone()
}

#[test]
fn artifact_and_unit_state_machines_reject_every_premature_completion() {
    let mut artifact = ArtifactLifecycle::new("a");
    assert_eq!(artifact.stage(), ArtifactStage::Open);
    assert_eq!(
        artifact.validate().expect_err("open is not valid").code(),
        "GPS-IO-005"
    );
    assert_eq!(
        artifact
            .catalog()
            .expect_err("open is not catalogable")
            .code(),
        "GPS-IO-005"
    );
    artifact.close().expect("close");
    assert_eq!(
        artifact.close().expect_err("double close").code(),
        "GPS-IO-005"
    );
    artifact.validate().expect("validate");
    artifact.catalog().expect("catalog");
    assert_eq!(
        artifact.catalog().expect_err("double catalog").code(),
        "GPS-IO-005"
    );
    assert_eq!(
        artifact
            .fail()
            .expect_err("cataloged artifact remains valid")
            .code(),
        "GPS-IO-005"
    );

    let mut failed = ArtifactLifecycle::new("failed");
    failed.fail().expect("open may fail");
    assert_eq!(
        failed.close().expect_err("failed cannot close").code(),
        "GPS-IO-005"
    );

    let mut unit = PublishableUnit::new("run", "execution", ["a".to_owned(), "b".to_owned()]);
    unit.record_artifact("a", ArtifactStage::Cataloged)
        .expect("known");
    assert_eq!(
        unit.complete().expect_err("b remains open").code(),
        "GPS-IO-005"
    );
    unit.record_artifact("b", ArtifactStage::Validated)
        .expect("known");
    assert_eq!(
        unit.record_artifact("b", ArtifactStage::Closed)
            .expect_err("validated artifact cannot regress")
            .code(),
        "GPS-IO-005"
    );
    assert_eq!(
        unit.complete().expect_err("b not cataloged").code(),
        "GPS-IO-005"
    );
    unit.record_artifact("b", ArtifactStage::Cataloged)
        .expect("known");
    unit.complete().expect("all cataloged");
    assert_eq!(unit.record().state, UnitState::Complete);
    assert_eq!(
        unit.record_artifact("a", ArtifactStage::Open)
            .expect_err("complete unit cannot regress")
            .code(),
        "GPS-IO-005"
    );
}

#[test]
fn admitted_configuration_bytes_are_frozen_and_execution_identity_ignores_source_paths() {
    let original = fs::read(fixture("general-seeded.toml")).expect("fixture");
    let first_source = temporary("first-source.toml");
    let second_source = temporary("second-source.toml");
    fs::write(&first_source, &original).expect("first source");
    fs::write(&second_source, &original).expect("second source");
    let first = load_and_validate_config(&first_source).expect("first admission");
    let second = load_and_validate_config(&second_source).expect("second admission");
    fs::write(&first_source, b"invalid after admission").expect("mutate source");

    let first_destination = temporary("first-output");
    let second_destination = temporary("second-output");
    execute_validated_config(&first, &first_destination).expect("first execution");
    execute_validated_config(&second, &second_destination).expect("second execution");
    assert_eq!(
        fs::read(first_destination.join("inputs/config/general.toml")).expect("frozen copy"),
        original
    );
    let first_manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(first_destination.join("manifest.json")).expect("first manifest"),
    )
    .expect("first JSON");
    let second_manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(second_destination.join("manifest.json")).expect("second manifest"),
    )
    .expect("second JSON");
    assert_eq!(
        first_manifest["execution_id"],
        second_manifest["execution_id"]
    );

    fs::remove_file(first_source).expect("first source cleanup");
    fs::remove_file(second_source).expect("second source cleanup");
    fs::remove_dir_all(first_destination).expect("first destination cleanup");
    fs::remove_dir_all(second_destination).expect("second destination cleanup");
}

#[test]
fn scenario_inputs_are_copied_exactly_and_invalid_scenario_ids_create_no_workspace() {
    let destination = temporary("scenario-output");
    let general_source = fixture("general-evolution.toml");
    let scenario_source = fixture("scenario-evolution.toml");
    let output = Command::new(BINARY)
        .args(["run"])
        .arg(&general_source)
        .args(["--output"])
        .arg(&destination)
        .output()
        .expect("scenario run");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read(destination.join("inputs/config/general.toml")).expect("general copy"),
        fs::read(&general_source).expect("general source")
    );
    assert_eq!(
        fs::read(destination.join("inputs/config/scenarios/evolution-small.toml"))
            .expect("scenario copy"),
        fs::read(&scenario_source).expect("scenario source")
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
            .expect("JSON");
    assert_eq!(
        manifest["configuration_inputs"]
            .as_array()
            .expect("configuration inputs")
            .len(),
        2
    );
    assert_eq!(manifest["events"].as_array().expect("events").len(), 1);
    fs::remove_dir_all(destination).expect("destination cleanup");

    let invalid_root = temporary("invalid-scenario-root");
    fs::create_dir(&invalid_root).expect("invalid root");
    let invalid_general = invalid_root.join("general-evolution.toml");
    let invalid_scenario = invalid_root.join("scenario-evolution.toml");
    fs::write(
        &invalid_general,
        include_str!("fixtures/config/general-evolution.toml")
            .replace("evolution-small", "../evil"),
    )
    .expect("invalid general");
    fs::write(
        &invalid_scenario,
        include_str!("fixtures/config/scenario-evolution.toml")
            .replace("evolution-small", "../evil"),
    )
    .expect("invalid scenario");
    let rejected_destination = invalid_root.join("invalid-scenario-output");
    let rejected = Command::new(BINARY)
        .args(["run"])
        .arg(&invalid_general)
        .args(["--output"])
        .arg(&rejected_destination)
        .output()
        .expect("invalid scenario run");
    assert_eq!(rejected.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("GPS-DOMAIN-001"));
    assert!(!rejected_destination.exists());
    assert_eq!(
        fs::read_dir(&invalid_root)
            .expect("parent entries")
            .filter_map(Result::ok)
            .filter(|entry| {
                entry
                    .file_name()
                    .to_string_lossy()
                    .contains("invalid-scenario-output.gps-incomplete")
            })
            .count(),
        0
    );
    fs::remove_dir_all(invalid_root).expect("invalid root cleanup");
}

#[test]
fn complete_zero_one_and_two_event_results_have_resolvable_hashed_links_and_valid_json() {
    for (name, expected_events) in [
        ("general-minimal.toml", 0_usize),
        ("general-seeded.toml", 1),
        ("general-events.toml", 2),
    ] {
        let destination = temporary(name);
        let output = Command::new(BINARY)
            .args(["run"])
            .arg(fixture(name))
            .args(["--output"])
            .arg(&destination)
            .output()
            .expect("run CLI");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(!destination.join("INCOMPLETE").exists());
        validate_json(
            &destination.join("effective-config.json"),
            "effective-config-v1.schema.json",
            true,
        );
        validate_json(
            &destination.join("manifest.json"),
            "manifest-general-v1.schema.json",
            true,
        );
        validate_json(
            &destination.join("inputs/resolution/evidence.json"),
            "resolution-evidence-v1.schema.json",
            true,
        );
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
                .expect("JSON");
        assert_eq!(manifest["state"], "complete");
        assert_eq!(
            manifest["events"].as_array().expect("events").len(),
            expected_events
        );
        for artifact in manifest["artifacts"].as_array().expect("artifacts") {
            let path = destination.join(artifact["path"].as_str().expect("path"));
            assert!(path.starts_with(&destination));
            assert!(path.is_file(), "{}", path.display());
            let external = Command::new("sha256sum")
                .arg(&path)
                .output()
                .expect("sha256sum");
            assert!(external.status.success());
            let digest = String::from_utf8(external.stdout).expect("ASCII");
            assert_eq!(
                artifact["checksum"]["digest"].as_str().expect("digest"),
                digest.split_whitespace().next().expect("external digest")
            );
            assert_eq!(
                artifact["size_bytes"].as_u64().expect("size"),
                fs::metadata(path).expect("metadata").len()
            );
        }
        for event in manifest["events"].as_array().expect("events") {
            let path = destination.join(event["manifest"]["path"].as_str().expect("event path"));
            validate_json(&path, "manifest-event-v1.schema.json", true);
        }
        fs::remove_dir_all(&destination).expect("cleanup");
    }
}

#[test]
fn effective_json_is_semantically_stable_under_toml_reordering() {
    let path = temporary("reordered.toml");
    let seed_line = format!("route = \"direct\"\nseed = \"{}\"", "00".repeat(32));
    let first = include_str!("fixtures/config/general-minimal.toml").replacen(
        "route = \"direct\"",
        &seed_line,
        1,
    );
    fs::write(&path, first).expect("first config");
    let first_config = load_and_validate_config(&path).expect("first");
    let first_destination = temporary("first-effective");
    execute_validated_config(&first_config, &first_destination).expect("first execution");
    let first: serde_json::Value = serde_json::from_slice(
        &fs::read(first_destination.join("effective-config.json")).expect("first effective"),
    )
    .expect("first JSON");
    let second = include_str!("fixtures/config/general-minimal-reordered.toml").replacen(
        "route = \"direct\"",
        &seed_line,
        1,
    );
    fs::write(&path, second).expect("second config");
    let second_config = load_and_validate_config(&path).expect("second");
    let second_destination = temporary("second-effective");
    execute_validated_config(&second_config, &second_destination).expect("second execution");
    let second: serde_json::Value = serde_json::from_slice(
        &fs::read(second_destination.join("effective-config.json")).expect("second effective"),
    )
    .expect("second JSON");
    assert_eq!(first, second);
    fs::remove_file(path).expect("cleanup");
    fs::remove_dir_all(first_destination).expect("first cleanup");
    fs::remove_dir_all(second_destination).expect("second cleanup");
}

#[test]
fn manifest_preserves_every_dtc049_value_and_exact_f64_roundtrip() {
    let source = fixture("numeric-vectors.toml");
    let effective = load_and_validate_config(&source).expect("numeric vectors");
    let expected = effective.normalized_collections().to_vec();
    let destination = temporary("normalization");
    let output = Command::new(BINARY)
        .args(["run"])
        .arg(&source)
        .args(["--output"])
        .arg(&destination)
        .output()
        .expect("run");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
            .expect("JSON");
    let observed = manifest["normalizations"]
        .as_array()
        .expect("normalizations");
    assert_eq!(observed.len(), expected.len());
    for (document, expected) in observed.iter().zip(&expected) {
        assert_eq!(document["algorithm"], expected.algorithm);
        assert_eq!(document["algorithm_version"], expected.algorithm_version);
        assert_eq!(document["canonical_order"], expected.canonical_order);
        for (field, value) in [
            ("tolerance", expected.tolerance),
            ("sum_before", expected.sum_before),
            ("delta_before", expected.delta_before),
            ("sum_after", expected.sum_after),
            ("delta_after", expected.delta_after),
        ] {
            assert_eq!(
                document[field].as_f64().expect("f64").to_bits(),
                value.to_bits(),
                "{}:{field}",
                expected.parameter
            );
        }
        for (document_value, expected_value) in document["values"]
            .as_array()
            .expect("values")
            .iter()
            .zip(&expected.values)
        {
            assert_eq!(
                document_value["before"].as_f64().expect("before").to_bits(),
                expected_value.before.to_bits()
            );
            assert_eq!(
                document_value["after"].as_f64().expect("after").to_bits(),
                expected_value.after.to_bits()
            );
        }
        assert_eq!(document["specific_tolerances"], serde_json::json!([]));
    }
    fs::remove_dir_all(destination).expect("cleanup");
}

#[test]
fn schemas_reject_truncation_unknown_fields_and_wrong_class_or_version() {
    let root = temporary("schema-negative");
    fs::create_dir_all(&root).expect("root");
    let execution = root.join("execution");
    execute_validated_config(
        &load_and_validate_config(&fixture("general-seeded.toml")).expect("effective"),
        &execution,
    )
    .expect("execution");
    let effective: serde_json::Value = serde_json::from_slice(
        &fs::read(execution.join("effective-config.json")).expect("effective output"),
    )
    .expect("effective JSON");
    let valid = root.join("valid.json");
    fs::write(
        &valid,
        serde_json::to_vec_pretty(&effective).expect("serialize"),
    )
    .expect("valid");
    validate_json(&valid, "effective-config-v1.schema.json", true);
    fs::write(root.join("truncated.json"), b"{\"schema_version\":").expect("truncated");
    validate_json(
        &root.join("truncated.json"),
        "effective-config-v1.schema.json",
        false,
    );
    for (name, field, value) in [
        ("unknown", "unknown_field", serde_json::json!(true)),
        ("class", "effective_config_kind", serde_json::json!("wrong")),
        ("version", "schema_version", serde_json::json!("2")),
    ] {
        let mut corrupt = effective.clone();
        corrupt
            .as_object_mut()
            .expect("object")
            .insert(field.to_owned(), value);
        let path = root.join(format!("{name}.json"));
        fs::write(
            &path,
            serde_json::to_vec_pretty(&corrupt).expect("serialize"),
        )
        .expect("corrupt");
        validate_json(&path, "effective-config-v1.schema.json", false);
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn public_effective_document_builder_is_schema_valid_before_one_nested_mutation() {
    let config = load_and_validate_config(&fixture("general-seeded.toml")).expect("effective");
    let specs = config
        .sequences()
        .iter()
        .map(|sequence| GeneratedSequenceSpec {
            id: sequence.id(),
            class: match sequence.class() {
                genome_population_simulator::configuration::SequenceClass::Chromosome => {
                    ReferenceSequenceClass::Chromosome
                }
                genome_population_simulator::configuration::SequenceClass::Scaffold => {
                    ReferenceSequenceClass::Scaffold
                }
            },
            length: sequence.length(),
        })
        .collect::<Vec<_>>();
    let composition = config.normalized_collections()[0]
        .values
        .iter()
        .map(|value| BaseFrequency {
            base: match value.id.as_str() {
                "A" => b'A',
                "C" => b'C',
                "G" => b'G',
                "T" => b'T',
                other => panic!("unexpected admitted base ID {other}"),
            },
            probability: value.after,
        })
        .collect::<Vec<_>>();
    let mut reference = generate_reference(
        config.reference_id(),
        config.reproducibility().root_seed(),
        &specs,
        &composition,
    )
    .expect("reference");
    reference
        .bind_authoritative_resource(PathBuf::from("resources/references/ref-main.fa"))
        .expect("resource");
    let artifact = ArtifactRecord {
        id: "reference-fasta".to_owned(),
        kind: "reference_fasta".to_owned(),
        path: "resources/references/ref-main.fa".to_owned(),
        size_bytes: 1019,
        digest: "0".repeat(64),
        state: UnitState::Complete,
        validated: true,
    };
    let valid = effective_config_document_with_reference(
        &config,
        ReferenceOutput {
            reference: &reference,
            fasta: &artifact,
        },
    )
    .expect("public builder contract");
    let root = temporary("public-effective-builder");
    fs::create_dir(&root).expect("root");
    let valid_path = root.join("valid.json");
    fs::write(
        &valid_path,
        serde_json::to_vec_pretty(&valid).expect("serialize"),
    )
    .expect("valid document");
    validate_json(&valid_path, "effective-config-v1.schema.json", true);

    let mut invalid = valid;
    invalid["reference"]["sequences"][0]["unexpected"] = serde_json::json!(true);
    let invalid_path = root.join("invalid.json");
    fs::write(
        &invalid_path,
        serde_json::to_vec_pretty(&invalid).expect("serialize invalid"),
    )
    .expect("invalid document");
    validate_json(&invalid_path, "effective-config-v1.schema.json", false);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn run_rejection_creates_no_destination_and_existing_destination_is_not_clobbered() {
    let rejected = temporary("rejected");
    let invalid = Command::new(BINARY)
        .args(["run"])
        .arg(fixture("invalid-unknown-field.toml"))
        .args(["--output"])
        .arg(&rejected)
        .output()
        .expect("run invalid");
    assert_eq!(invalid.status.code(), Some(2));
    assert!(!rejected.exists());

    let destination = temporary("existing");
    fs::create_dir_all(&destination).expect("existing destination");
    fs::write(destination.join("sentinel"), b"keep").expect("sentinel");
    let output = Command::new(BINARY)
        .args(["run"])
        .arg(fixture("general-seeded.toml"))
        .args(["--output"])
        .arg(&destination)
        .output()
        .expect("run existing");
    assert_eq!(output.status.code(), Some(5));
    assert_eq!(
        fs::read(destination.join("sentinel")).expect("sentinel"),
        b"keep"
    );
    let provisional = provisional_siblings(&destination);
    assert_eq!(provisional.len(), 1);
    assert!(provisional[0].join("INCOMPLETE").is_file());
    fs::remove_dir_all(&provisional[0]).expect("provisional cleanup");
    fs::remove_dir_all(destination).expect("cleanup");
}

#[test]
fn two_processes_racing_for_one_destination_publish_exactly_one_result() {
    let destination = temporary("race");
    let mut first = Command::new(BINARY)
        .args(["run"])
        .arg(fixture("general-seeded.toml"))
        .args(["--output"])
        .arg(&destination)
        .spawn()
        .expect("first");
    let mut second = Command::new(BINARY)
        .args(["run"])
        .arg(fixture("general-seeded.toml"))
        .args(["--output"])
        .arg(&destination)
        .spawn()
        .expect("second");
    let statuses = [
        first.wait().expect("first status"),
        second.wait().expect("second status"),
    ];
    assert_eq!(statuses.iter().filter(|status| status.success()).count(), 1);
    assert_eq!(
        statuses
            .iter()
            .filter(|status| status.code() == Some(5))
            .count(),
        1
    );
    validate_json(
        &destination.join("manifest.json"),
        "manifest-general-v1.schema.json",
        true,
    );
    let provisional = provisional_siblings(&destination);
    assert_eq!(provisional.len(), 1);
    assert!(provisional[0].join("INCOMPLETE").is_file());
    fs::remove_dir_all(&provisional[0]).expect("losing provisional cleanup");
    fs::remove_dir_all(destination).expect("cleanup");
}

#[test]
fn runtime_structure_matches_jsonschema_and_generic_manifest_paths_are_closed() {
    let root = temporary("runtime-structure-corpus");
    fs::create_dir_all(&root).unwrap();
    for (_, _, relative, bytes) in &manifest_fixture().artifacts {
        let path = root.join(relative);
        fs::create_dir_all(path.parent().expect("fixture artifact parent")).unwrap();
        fs::write(path, bytes).unwrap();
    }
    let valid_path = root.join("valid.json");
    fs::write(&valid_path, structurally_valid_manifest_bytes()).unwrap();
    validate_json(&valid_path, "manifest-general-v1.schema.json", true);
    validate_general_manifest_bytes(&valid_path, structurally_valid_manifest_bytes())
        .expect("runtime accepts the structural positive");

    let valid: serde_json::Value =
        serde_json::from_slice(structurally_valid_manifest_bytes()).unwrap();
    let incomplete = serde_json::json!({
        "manifest_kind": "gps_execution_manifest",
        "state": "failed",
        "references": [],
        "reference_reductions": [],
        "artifacts": [],
        "units": [],
        "diagnostics": []
    });
    let mut missing_required = valid.clone();
    missing_required.as_object_mut().unwrap().remove("program");
    let mut unknown_root = valid.clone();
    unknown_root["unknown_runtime_bypass"] = serde_json::json!(true);
    let mut wrong_type = valid.clone();
    wrong_type["effective_config"]["validated"] = serde_json::json!("true");
    let mut invalid_nested_enum = valid.clone();
    invalid_nested_enum["references"][0]["origin"]["kind"] = serde_json::json!("invented");
    let mut wrong_const = valid.clone();
    wrong_const["manifest_kind"] = serde_json::json!("other_manifest");
    let mut wrong_pattern = valid.clone();
    wrong_pattern["execution_id"] = serde_json::json!("run_invalid");
    let mut zero_minimum = valid.clone();
    zero_minimum["references"][0]["total_genome_size"] = serde_json::json!(0);
    let mut empty_min_length = valid.clone();
    empty_min_length["program"]["name"] = serde_json::json!("");
    let mut invalid_version_pattern = valid.clone();
    invalid_version_pattern["program"]["version"] = serde_json::json!("not-a-version");
    let mut duplicate_unique_item = valid.clone();
    duplicate_unique_item["outputs"]["minimum"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!("manifest-general"));
    let mut wrong_prefix_items = valid.clone();
    wrong_prefix_items["outputs"]["disabled"]
        .as_array_mut()
        .unwrap()
        .swap(0, 1);
    let mut missing_dependency_item = valid.clone();
    missing_dependency_item["dependencies"]
        .as_array_mut()
        .unwrap()
        .pop();
    let mut inconsistent_if_else = valid.clone();
    inconsistent_if_else["references"][0]["origin"] = serde_json::json!({
        "kind": "supplied",
        "original_path": "input.fa"
    });
    for (name, mutation) in [
        ("incomplete", incomplete),
        ("missing-required", missing_required),
        ("unknown-root", unknown_root),
        ("wrong-type", wrong_type),
        ("invalid-nested-enum", invalid_nested_enum),
        ("wrong-const", wrong_const),
        ("wrong-pattern", wrong_pattern),
        ("zero-minimum", zero_minimum),
        ("empty-min-length", empty_min_length),
        ("invalid-version-pattern", invalid_version_pattern),
        ("duplicate-unique-item", duplicate_unique_item),
        ("wrong-prefix-items", wrong_prefix_items),
        ("missing-dependency-item", missing_dependency_item),
        ("inconsistent-if-else", inconsistent_if_else),
    ] {
        let path = root.join(format!("{name}.json"));
        let bytes = serde_json::to_vec_pretty(&mutation).unwrap();
        fs::write(&path, &bytes).unwrap();
        validate_json(&path, "manifest-general-v1.schema.json", false);
        let error = match validate_general_manifest_bytes(&path, &bytes) {
            Err(error) => error,
            Ok(_) => panic!("runtime structural gate accepted Python rejection: {name}"),
        };
        assert_eq!(error.code(), "GPS-FMT-006");
    }

    let destination = root.join("generic-bypass");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let direct = workspace
        .write_artifact(
            "manifest-general",
            "manifest_general",
            Path::new("manifest.json"),
            structurally_valid_manifest_bytes(),
        )
        .expect_err("raw manifest bytes cannot use generic writer");
    assert_eq!(direct.code(), "GPS-IO-005");
    let streamed = workspace
        .write_streamed_artifact(
            "manifest-general",
            "manifest_general",
            Path::new("manifest.json"),
            |_| Ok(()),
        )
        .expect_err("streaming cannot bypass the capability");
    assert_eq!(streamed.code(), "GPS-IO-005");
    workspace
        .write_artifact(
            "manifest-disguised",
            "other",
            Path::new("manifest.json"),
            structurally_valid_manifest_bytes(),
        )
        .expect("non-manifest artifact for publication rejection");
    let provisional = provisional_workspace(&destination);
    let publish = workspace
        .publish()
        .expect_err("publication requires the capability-backed artifact class");
    assert_eq!(publish.code(), "GPS-IO-005");
    fs::remove_dir_all(provisional).expect("provisional cleanup");
    fs::remove_dir_all(root).expect("corpus cleanup");
}

#[test]
fn manifest_capability_is_bound_to_the_exact_workspace_catalog_and_files() {
    let root = temporary("workspace-binding");
    fs::create_dir_all(&root).unwrap();
    let source_destination = root.join("source");
    let (source_workspace, manifest) = matching_workspace(&source_destination);

    let identical_destination = root.join("identical");
    let (mut identical, identical_manifest) = matching_workspace(&identical_destination);
    assert_eq!(
        identical
            .write_general_manifest(&manifest)
            .expect_err("identical bytes and catalog do not transfer workspace authority")
            .code(),
        "GPS-IO-005"
    );
    identical
        .write_general_manifest(&identical_manifest)
        .expect("the capability created by this workspace remains valid");
    fs::remove_dir_all(provisional_workspace(&identical_destination)).unwrap();

    let empty_destination = root.join("empty");
    let mut empty = ExecutionWorkspace::create(&empty_destination).unwrap();
    assert_eq!(
        empty
            .write_general_manifest(&manifest)
            .expect_err("a foreign capability cannot authorize an empty workspace")
            .code(),
        "GPS-IO-005"
    );
    fs::remove_dir_all(provisional_workspace(&empty_destination)).unwrap();

    for (missing_index, (_, missing_kind, _, _)) in manifest_fixture().artifacts.iter().enumerate()
    {
        let destination = root.join(format!("missing-{missing_kind}"));
        let mut workspace = ExecutionWorkspace::create(&destination).unwrap();
        for (index, (id, kind, path, bytes)) in manifest_fixture().artifacts.iter().enumerate() {
            if index != missing_index {
                workspace.write_artifact(id, kind, path, bytes).unwrap();
            }
        }
        assert_eq!(
            workspace
                .write_general_manifest(&manifest)
                .expect_err("every artifact class is required")
                .code(),
            "GPS-IO-005"
        );
        fs::remove_dir_all(provisional_workspace(&destination)).unwrap();
    }

    let extra_destination = root.join("extra");
    let (mut extra, extra_manifest) = matching_workspace(&extra_destination);
    extra
        .write_artifact("extra", "extra_kind", Path::new("extra.bin"), b"extra")
        .unwrap();
    assert_eq!(
        extra
            .write_general_manifest(&extra_manifest)
            .expect_err("workspace extras are forbidden")
            .code(),
        "GPS-IO-005"
    );
    fs::remove_dir_all(provisional_workspace(&extra_destination)).unwrap();

    for mutation in ["class", "path"] {
        let destination = root.join(mutation);
        let mut workspace = ExecutionWorkspace::create(&destination).unwrap();
        for (index, (id, kind, path, bytes)) in manifest_fixture().artifacts.iter().enumerate() {
            let changed_kind = if index == 0 && mutation == "class" {
                "different_kind"
            } else {
                kind
            };
            let changed_path = if index == 0 && mutation == "path" {
                Path::new("different/path.bin")
            } else {
                path
            };
            workspace
                .write_artifact(id, changed_kind, changed_path, bytes)
                .unwrap();
        }
        assert_eq!(
            workspace
                .write_general_manifest(&manifest)
                .expect_err("catalog metadata must match exactly")
                .code(),
            "GPS-IO-005"
        );
        fs::remove_dir_all(provisional_workspace(&destination)).unwrap();
    }

    fs::remove_dir_all(provisional_workspace(&source_destination)).unwrap();
    drop(source_workspace);
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn truncated_temporary_is_never_cataloged_or_published() {
    let destination = temporary("truncated-destination");
    let workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let provisional = provisional_workspace(&destination);
    let truncated = workspace
        .create_truncated_temporary(Path::new("artifacts/truncated.json"), b"{\"x\":")
        .expect("truncated fixture");
    assert!(truncated.is_file());
    assert!(
        truncated
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(".gps-tmp"))
    );
    assert!(workspace.artifacts().is_empty());
    let error = workspace
        .publish()
        .expect_err("missing validated manifest cannot publish");
    assert_eq!(error.code(), "GPS-IO-005");
    assert!(provisional.join("INCOMPLETE").is_file());
    assert!(!destination.exists());
    fs::remove_dir_all(provisional).expect("cleanup");
}

#[test]
fn uncontrolled_sigkill_leaves_a_provisional_incomplete_marker() {
    if std::env::var_os("GPS_F14_SIGKILL_CHILD").is_some() {
        let destination =
            PathBuf::from(std::env::var_os("GPS_F14_DESTINATION").expect("destination"));
        let ready = PathBuf::from(std::env::var_os("GPS_F14_READY").expect("ready"));
        let _workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        fs::write(
            &ready,
            provisional_workspace(&destination).display().to_string(),
        )
        .expect("ready");
        loop {
            thread::sleep(Duration::from_secs(1));
        }
    }

    let destination = temporary("sigkill-destination");
    let ready = temporary("sigkill-ready");
    let mut child = Command::new(std::env::current_exe().expect("current test binary"))
        .args([
            "--exact",
            "uncontrolled_sigkill_leaves_a_provisional_incomplete_marker",
            "--nocapture",
        ])
        .env("GPS_F14_SIGKILL_CHILD", "1")
        .env("GPS_F14_DESTINATION", &destination)
        .env("GPS_F14_READY", &ready)
        .spawn()
        .expect("child");
    let deadline = Instant::now() + Duration::from_secs(10);
    while !ready.exists() && Instant::now() < deadline {
        thread::sleep(Duration::from_millis(10));
    }
    assert!(ready.exists(), "child did not create provisional");
    let provisional = PathBuf::from(fs::read_to_string(&ready).expect("ready path"));
    let kill = Command::new("kill")
        .args(["-9", &child.id().to_string()])
        .status()
        .expect("kill");
    assert!(kill.success());
    let status = child.wait().expect("killed status");
    assert!(!status.success());
    assert!(!destination.exists());
    assert!(provisional.join("INCOMPLETE").is_file());
    fs::remove_file(ready).expect("ready cleanup");
    fs::remove_dir_all(provisional).expect("provisional cleanup");
}

fn validate_json(instance: &Path, schema: &str, expected_success: bool) {
    let output = python_jsonschema_output(instance, schema);
    assert_eq!(
        output.status.success(),
        expected_success,
        "schema={schema} instance={} stderr={}",
        instance.display(),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn python_jsonschema_output(instance: &Path, schema: &str) -> std::process::Output {
    let schema_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("schemas");
    let base_uri = format!("file://{}/", schema_root.display());
    Command::new("python3")
        .args(["-m", "jsonschema", "--base-uri", &base_uri, "-i"])
        .arg(instance)
        .arg(schema)
        .current_dir(schema_root)
        .output()
        .expect("independent JSON Schema validator from A7")
}
