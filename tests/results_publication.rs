use std::cell::Cell;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;
use std::time::{Duration, Instant};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::{
    ExecutionFault, ExecutionOutcome, execute_validated_config, execute_validated_config_with,
};
use genome_population_simulator::publication::{
    ArtifactLifecycle, ArtifactStage, AtomicPublisher, ExecutionWorkspace, LinuxAtomicPublisher,
    PublishableUnit, UnsupportedAtomicPublisher,
};
use genome_population_simulator::reference::SequenceClass as ReferenceSequenceClass;
use genome_population_simulator::reference_generation::{
    BaseFrequency, GeneratedSequenceSpec, generate_reference,
};
use genome_population_simulator::results::{
    ArtifactRecord, ReferenceOutput, UnitState, effective_config_document_with_reference,
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
    execute_validated_config_with(
        &first,
        &first_destination,
        ExecutionFault::None,
        &LinuxAtomicPublisher,
    )
    .expect("first execution");
    execute_validated_config_with(
        &second,
        &second_destination,
        ExecutionFault::None,
        &LinuxAtomicPublisher,
    )
    .expect("second execution");
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
fn schemas_reject_unknown_fields_in_every_nested_f1_4_contract() {
    let root = temporary("schema-nested-negative");
    fs::create_dir(&root).expect("root");

    let complete_destination = root.join("complete");
    let output = Command::new(BINARY)
        .args(["run"])
        .arg(fixture("general-seeded.toml"))
        .args(["--output"])
        .arg(&complete_destination)
        .output()
        .expect("complete run");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let general: serde_json::Value = serde_json::from_slice(
        &fs::read(complete_destination.join("manifest.json")).expect("general manifest"),
    )
    .expect("general JSON");
    let event_path = complete_destination.join(
        general["events"][0]["manifest"]["path"]
            .as_str()
            .expect("event manifest path"),
    );
    let event: serde_json::Value =
        serde_json::from_slice(&fs::read(event_path).expect("event manifest")).expect("event JSON");

    for (name, pointer) in [
        ("general-program", "/program"),
        ("general-platform", "/platform"),
        ("general-toolchain", "/toolchain"),
        ("general-dependency", "/dependencies/0"),
        ("general-check", "/checks/0"),
        ("general-input", "/configuration_inputs/0"),
        ("general-event", "/events/0"),
        ("general-outputs", "/outputs"),
        ("general-unit", "/units/0"),
    ] {
        reject_nested_unknown(
            &root,
            name,
            general.clone(),
            pointer,
            "manifest-general-v1.schema.json",
        );
    }
    reject_missing_nested_required(
        &root,
        "general-program-missing-name",
        general.clone(),
        "/program",
        "name",
        "manifest-general-v1.schema.json",
    );
    reject_nested_value(
        &root,
        "general-input-wrong-path-type",
        general.clone(),
        "/configuration_inputs/0/original_path",
        serde_json::json!(false),
        "manifest-general-v1.schema.json",
    );
    for (name, pointer) in [("event-outputs", "/outputs"), ("event-unit", "/units/0")] {
        reject_nested_unknown(
            &root,
            name,
            event.clone(),
            pointer,
            "manifest-event-v1.schema.json",
        );
    }
    reject_missing_nested_required(
        &root,
        "event-outputs-missing-disabled",
        event.clone(),
        "/outputs",
        "disabled",
        "manifest-event-v1.schema.json",
    );
    reject_nested_value(
        &root,
        "event-unit-wrong-exportable-type",
        event.clone(),
        "/units/0/exportable",
        serde_json::json!("true"),
        "manifest-event-v1.schema.json",
    );

    let rich_destination = root.join("rich-effective-result");
    execute_validated_config(
        &load_and_validate_config(&fixture("general-effective.toml")).expect("rich effective"),
        &rich_destination,
    )
    .expect("rich execution");
    let rich_effective: serde_json::Value = serde_json::from_slice(
        &fs::read(rich_destination.join("effective-config.json")).expect("rich output"),
    )
    .expect("rich JSON");
    for (name, pointer) in [
        ("effective-reference", "/reference"),
        ("effective-randomness", "/randomness"),
        ("effective-controlled-set", "/controlled_frequency_sets/0"),
        (
            "effective-controlled-count",
            "/controlled_frequency_sets/0/counts/0",
        ),
        ("effective-value", "/effective_values/0"),
        ("effective-selectors", "/effective_values/0/selectors"),
        ("effective-provenance", "/effective_values/0/provenance"),
    ] {
        reject_nested_unknown(
            &root,
            name,
            rich_effective.clone(),
            pointer,
            "effective-config-v1.schema.json",
        );
    }
    reject_missing_nested_required(
        &root,
        "effective-provenance-missing-location",
        rich_effective.clone(),
        "/effective_values/0/provenance",
        "location",
        "effective-config-v1.schema.json",
    );
    reject_nested_value(
        &root,
        "effective-controlled-negative-count",
        rich_effective,
        "/controlled_frequency_sets/0/counts/0/count",
        serde_json::json!(-1),
        "effective-config-v1.schema.json",
    );

    let failed_destination = root.join("failed");
    let effective = load_and_validate_config(&fixture("general-seeded.toml")).expect("effective");
    execute_validated_config_with(
        &effective,
        &failed_destination,
        ExecutionFault::AfterEffectiveConfig,
        &LinuxAtomicPublisher,
    )
    .expect("failed result");
    let failed: serde_json::Value = serde_json::from_slice(
        &fs::read(failed_destination.join("manifest.json")).expect("failed manifest"),
    )
    .expect("failed JSON");
    reject_nested_unknown(
        &root,
        "general-diagnostic",
        failed,
        "/diagnostics/0",
        "manifest-general-v1.schema.json",
    );

    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn schemas_reject_semantically_contradictory_f1_4_documents() {
    let root = temporary("schema-semantic-negative");
    fs::create_dir(&root).expect("root");

    let complete_destination = root.join("complete");
    let output = Command::new(BINARY)
        .args(["run"])
        .arg(fixture("general-seeded.toml"))
        .args(["--output"])
        .arg(&complete_destination)
        .output()
        .expect("complete run");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let complete: serde_json::Value = serde_json::from_slice(
        &fs::read(complete_destination.join("manifest.json")).expect("general manifest"),
    )
    .expect("general JSON");
    let event_path = complete_destination.join(
        complete["events"][0]["manifest"]["path"]
            .as_str()
            .expect("event manifest path"),
    );
    let event: serde_json::Value =
        serde_json::from_slice(&fs::read(event_path).expect("event manifest")).expect("event JSON");

    reject_nested_value(
        &root,
        "complete-declared-failed-without-diagnostic",
        complete.clone(),
        "/state",
        serde_json::json!("failed"),
        "manifest-general-v1.schema.json",
    );

    let mut duplicate_level = complete.clone();
    *duplicate_level
        .pointer_mut("/reproducibility/0/level")
        .expect("biological level") = serde_json::json!("technical");
    *duplicate_level
        .pointer_mut("/reproducibility/0/status")
        .expect("biological status") = serde_json::json!("guaranteed");
    write_corrupt_and_reject(
        &root,
        "reproducibility-duplicates-technical-and-omits-biological",
        &duplicate_level,
        "manifest-general-v1.schema.json",
    );

    reject_nested_value(
        &root,
        "complete-without-effective-config",
        complete.clone(),
        "/effective_config",
        serde_json::Value::Null,
        "manifest-general-v1.schema.json",
    );

    let mut complete_with_failed_unit = complete.clone();
    *complete_with_failed_unit
        .pointer_mut("/units/0/state")
        .expect("execution state") = serde_json::json!("failed");
    *complete_with_failed_unit
        .pointer_mut("/units/0/exportable")
        .expect("execution exportability") = serde_json::json!(false);
    *complete_with_failed_unit
        .pointer_mut("/units/0/diagnostic_code")
        .expect("execution diagnostic") = serde_json::json!("GPS-DOMAIN-003");
    write_corrupt_and_reject(
        &root,
        "complete-with-failed-execution-unit",
        &complete_with_failed_unit,
        "manifest-general-v1.schema.json",
    );

    let failed_destination = root.join("failed");
    let effective = load_and_validate_config(&fixture("general-seeded.toml")).expect("effective");
    execute_validated_config_with(
        &effective,
        &failed_destination,
        ExecutionFault::AfterEffectiveConfig,
        &LinuxAtomicPublisher,
    )
    .expect("failed result");
    let failed: serde_json::Value = serde_json::from_slice(
        &fs::read(failed_destination.join("manifest.json")).expect("failed manifest"),
    )
    .expect("failed JSON");

    reject_nested_value(
        &root,
        "complete-with-diagnostic",
        complete.clone(),
        "/diagnostics",
        failed["diagnostics"].clone(),
        "manifest-general-v1.schema.json",
    );

    reject_nested_value(
        &root,
        "failed-without-diagnostics",
        failed.clone(),
        "/diagnostics",
        serde_json::json!([]),
        "manifest-general-v1.schema.json",
    );
    reject_nested_value(
        &root,
        "failed-execution-unit-without-code",
        failed.clone(),
        "/units/0/diagnostic_code",
        serde_json::Value::Null,
        "manifest-general-v1.schema.json",
    );
    reject_nested_value(
        &root,
        "failed-execution-unit-exportable",
        failed.clone(),
        "/units/0/exportable",
        serde_json::json!(true),
        "manifest-general-v1.schema.json",
    );

    let mut failed_with_complete_unit = failed;
    *failed_with_complete_unit
        .pointer_mut("/units/0/state")
        .expect("execution state") = serde_json::json!("complete");
    *failed_with_complete_unit
        .pointer_mut("/units/0/exportable")
        .expect("execution exportability") = serde_json::json!(true);
    *failed_with_complete_unit
        .pointer_mut("/units/0/diagnostic_code")
        .expect("execution diagnostic") = serde_json::Value::Null;
    write_corrupt_and_reject(
        &root,
        "failed-with-complete-execution-unit",
        &failed_with_complete_unit,
        "manifest-general-v1.schema.json",
    );

    reject_nested_value(
        &root,
        "event-unit-not-exportable",
        event.clone(),
        "/units/0/exportable",
        serde_json::json!(false),
        "manifest-event-v1.schema.json",
    );
    reject_nested_value(
        &root,
        "event-unit-wrong-kind",
        event.clone(),
        "/units/0/kind",
        serde_json::json!("execution"),
        "manifest-event-v1.schema.json",
    );
    reject_nested_value(
        &root,
        "event-unit-wrong-state",
        event,
        "/units/0/state",
        serde_json::json!("failed"),
        "manifest-event-v1.schema.json",
    );

    let mut duplicate_dependency = complete.clone();
    *duplicate_dependency
        .pointer_mut("/dependencies/1/name")
        .expect("second dependency name") = duplicate_dependency["dependencies"][0]["name"].clone();
    *duplicate_dependency
        .pointer_mut("/dependencies/1/version")
        .expect("second dependency version") = serde_json::json!("9.9.9");
    write_corrupt_and_reject(
        &root,
        "dependencies-duplicate-one-and-omit-one",
        &duplicate_dependency,
        "manifest-general-v1.schema.json",
    );
    reject_nested_value(
        &root,
        "dependency-version-outside-f1-4-lock",
        complete,
        "/dependencies/0/version",
        serde_json::json!("9.9.9"),
        "manifest-general-v1.schema.json",
    );

    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn controlled_failure_publishes_failed_manifest_and_interruption_keeps_marker() {
    let effective = load_and_validate_config(&fixture("general-seeded.toml")).expect("effective");
    for fault in [
        ExecutionFault::AfterConfigurationCopy,
        ExecutionFault::AfterEffectiveConfig,
    ] {
        let destination = temporary("controlled-failure");
        let outcome =
            execute_validated_config_with(&effective, &destination, fault, &LinuxAtomicPublisher)
                .expect("controlled failure closes coherently");
        assert_eq!(
            outcome,
            ExecutionOutcome::Published {
                destination: destination.clone(),
                state: UnitState::Failed
            }
        );
        validate_json(
            &destination.join("manifest.json"),
            "manifest-general-v1.schema.json",
            true,
        );
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
                .expect("JSON");
        assert_eq!(manifest["state"], "failed");
        assert_eq!(manifest["diagnostics"][0]["code"], "GPS-DOMAIN-003");
        assert!(!destination.join("INCOMPLETE").exists());
        fs::remove_dir_all(destination).expect("cleanup");
    }

    let destination = temporary("interrupted");
    let outcome = execute_validated_config_with(
        &effective,
        &destination,
        ExecutionFault::InterruptedAfterWorkspace,
        &LinuxAtomicPublisher,
    )
    .expect("interruption fixture");
    let ExecutionOutcome::Interrupted { provisional } = outcome else {
        panic!("expected provisional");
    };
    assert!(!destination.exists());
    assert!(provisional.join("INCOMPLETE").is_file());
    fs::remove_dir_all(provisional).expect("cleanup");
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
fn unsupported_publisher_never_falls_back_or_removes_incomplete_marker() {
    let destination = temporary("unsupported");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    workspace
        .write_artifact("a", "test", Path::new("a"), b"a")
        .expect("artifact");
    workspace
        .write_artifact(
            "manifest-general",
            "manifest_general",
            Path::new("manifest.json"),
            b"{}\n",
        )
        .expect("cataloged manifest fixture");
    let provisional = workspace.root().to_path_buf();
    let error = workspace
        .publish(&UnsupportedAtomicPublisher)
        .expect_err("unsupported");
    assert_eq!(error.code(), "GPS-IO-004");
    assert!(provisional.join("INCOMPLETE").is_file());
    assert!(!destination.exists());
    fs::remove_dir_all(provisional).expect("cleanup");
}

#[test]
fn incomplete_marker_is_present_during_atomic_publication_and_survives_publisher_failure() {
    for fail in [false, true] {
        let destination = temporary(if fail {
            "marker-publisher-failure"
        } else {
            "marker-publisher-success"
        });
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        workspace
            .write_artifact(
                "manifest-general",
                "manifest_general",
                Path::new("manifest.json"),
                b"{}\n",
            )
            .expect("cataloged manifest fixture");
        let provisional = workspace.root().to_path_buf();
        let publisher = MarkerObservingPublisher::new(fail);
        let outcome = workspace.publish(&publisher);
        assert!(
            publisher.marker_present_at_entry.get(),
            "INCOMPLETE must cover the entire pre-publication interval"
        );
        if fail {
            assert_eq!(outcome.expect_err("simulated failure").code(), "GPS-IO-004");
            assert!(provisional.join("INCOMPLETE").is_file());
            assert!(!destination.exists());
            fs::remove_dir_all(provisional).expect("failed provisional cleanup");
        } else {
            assert_eq!(outcome.expect("publication"), destination);
            assert!(!provisional.exists());
            assert!(!destination.join("INCOMPLETE").exists());
            fs::remove_dir_all(destination).expect("published cleanup");
        }
    }
}

#[test]
fn marker_restoration_and_removal_failures_preserve_the_incomplete_invariant() {
    let failed_destination = temporary("marker-restoration");
    let mut workspace = ExecutionWorkspace::create(&failed_destination).expect("workspace");
    workspace
        .write_artifact(
            "manifest-general",
            "manifest_general",
            Path::new("manifest.json"),
            b"{}\n",
        )
        .expect("manifest");
    let provisional = workspace.root().to_path_buf();
    let error = workspace
        .publish(&MarkerRemovingFailurePublisher)
        .expect_err("publisher failure");
    assert_eq!(error.code(), "GPS-IO-004");
    assert!(provisional.join("INCOMPLETE").is_file());
    assert!(!failed_destination.exists());
    fs::remove_dir_all(provisional).expect("restored provisional cleanup");

    let removal_destination = temporary("marker-removal-failure");
    let mut workspace = ExecutionWorkspace::create(&removal_destination).expect("workspace");
    workspace
        .write_artifact(
            "manifest-general",
            "manifest_general",
            Path::new("manifest.json"),
            b"{}\n",
        )
        .expect("manifest");
    let error = workspace
        .publish(&ReadOnlyDestinationPublisher)
        .expect_err("marker removal must be reported");
    assert_eq!(error.code(), "GPS-IO-001");
    assert!(removal_destination.join("INCOMPLETE").is_file());
    let mut permissions = fs::metadata(&removal_destination)
        .expect("destination metadata")
        .permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(&removal_destination, permissions).expect("restore permissions");
    fs::remove_dir_all(removal_destination).expect("removal failure cleanup");
}

#[test]
fn truncated_temporary_is_never_cataloged_or_published() {
    let destination = temporary("truncated-destination");
    let workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let provisional = workspace.root().to_path_buf();
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
        .publish(&LinuxAtomicPublisher)
        .expect_err("missing validated manifest cannot publish");
    assert_eq!(error.code(), "GPS-IO-005");
    assert!(provisional.join("INCOMPLETE").is_file());
    assert!(!destination.exists());
    fs::remove_dir_all(provisional).expect("cleanup");
}

#[test]
fn controlled_failure_process_exits_one_and_publishes_failed_manifest() {
    if std::env::var_os("GPS_F14_CONTROLLED_CHILD").is_some() {
        let destination =
            PathBuf::from(std::env::var_os("GPS_F14_DESTINATION").expect("destination"));
        let effective =
            load_and_validate_config(&fixture("general-seeded.toml")).expect("effective");
        let outcome = execute_validated_config_with(
            &effective,
            &destination,
            ExecutionFault::AfterEffectiveConfig,
            &LinuxAtomicPublisher,
        )
        .expect("controlled closure");
        assert!(matches!(
            outcome,
            ExecutionOutcome::Published {
                state: UnitState::Failed,
                ..
            }
        ));
        std::process::exit(1);
    }

    let destination = temporary("controlled-process");
    let status = Command::new(std::env::current_exe().expect("current test binary"))
        .args([
            "--exact",
            "controlled_failure_process_exits_one_and_publishes_failed_manifest",
            "--nocapture",
        ])
        .env("GPS_F14_CONTROLLED_CHILD", "1")
        .env("GPS_F14_DESTINATION", &destination)
        .status()
        .expect("controlled child");
    assert_eq!(status.code(), Some(1));
    validate_json(
        &destination.join("manifest.json"),
        "manifest-general-v1.schema.json",
        true,
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
            .expect("JSON");
    assert_eq!(manifest["state"], "failed");
    assert!(!destination.join("INCOMPLETE").exists());
    fs::remove_dir_all(destination).expect("cleanup");
}

#[test]
fn uncontrolled_sigkill_leaves_a_provisional_incomplete_marker() {
    if std::env::var_os("GPS_F14_SIGKILL_CHILD").is_some() {
        let destination =
            PathBuf::from(std::env::var_os("GPS_F14_DESTINATION").expect("destination"));
        let ready = PathBuf::from(std::env::var_os("GPS_F14_READY").expect("ready"));
        let workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        fs::write(&ready, workspace.root().display().to_string()).expect("ready");
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
    let schema_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("schemas");
    let base_uri = format!("file://{}/", schema_root.display());
    let output = Command::new("python3")
        .args(["-m", "jsonschema", "--base-uri", &base_uri, "-i"])
        .arg(instance)
        .arg(schema)
        .current_dir(schema_root)
        .output()
        .expect("independent JSON Schema validator from A7");
    assert_eq!(
        output.status.success(),
        expected_success,
        "schema={schema} instance={} stderr={}",
        instance.display(),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn reject_nested_unknown(
    root: &Path,
    name: &str,
    mut document: serde_json::Value,
    pointer: &str,
    schema: &str,
) {
    write_base_and_accept(root, name, &document, schema);
    document
        .pointer_mut(pointer)
        .unwrap_or_else(|| panic!("missing JSON pointer {pointer}"))
        .as_object_mut()
        .unwrap_or_else(|| panic!("JSON pointer {pointer} is not an object"))
        .insert(
            "unexpected_nested_field".to_owned(),
            serde_json::json!(true),
        );
    let path = root.join(format!("{name}.json"));
    fs::write(
        &path,
        serde_json::to_vec_pretty(&document).expect("serialize corruption"),
    )
    .expect("write corruption");
    validate_json(&path, schema, false);
}

fn reject_missing_nested_required(
    root: &Path,
    name: &str,
    mut document: serde_json::Value,
    pointer: &str,
    field: &str,
    schema: &str,
) {
    write_base_and_accept(root, name, &document, schema);
    document
        .pointer_mut(pointer)
        .unwrap_or_else(|| panic!("missing JSON pointer {pointer}"))
        .as_object_mut()
        .unwrap_or_else(|| panic!("JSON pointer {pointer} is not an object"))
        .remove(field)
        .unwrap_or_else(|| panic!("missing required field {field} at {pointer}"));
    write_corrupt_and_reject(root, name, &document, schema);
}

fn reject_nested_value(
    root: &Path,
    name: &str,
    mut document: serde_json::Value,
    pointer: &str,
    value: serde_json::Value,
    schema: &str,
) {
    write_base_and_accept(root, name, &document, schema);
    *document
        .pointer_mut(pointer)
        .unwrap_or_else(|| panic!("missing JSON pointer {pointer}")) = value;
    write_corrupt_and_reject(root, name, &document, schema);
}

fn write_base_and_accept(root: &Path, name: &str, document: &serde_json::Value, schema: &str) {
    let path = root.join(format!("{name}-base-valid.json"));
    fs::write(
        &path,
        serde_json::to_vec_pretty(document).expect("serialize valid base"),
    )
    .expect("write valid base");
    validate_json(&path, schema, true);
}

fn write_corrupt_and_reject(root: &Path, name: &str, document: &serde_json::Value, schema: &str) {
    let path = root.join(format!("{name}.json"));
    fs::write(
        &path,
        serde_json::to_vec_pretty(document).expect("serialize corruption"),
    )
    .expect("write corruption");
    validate_json(&path, schema, false);
}

struct MarkerObservingPublisher {
    marker_present_at_entry: Cell<bool>,
    fail: bool,
}

impl MarkerObservingPublisher {
    fn new(fail: bool) -> Self {
        Self {
            marker_present_at_entry: Cell::new(false),
            fail,
        }
    }
}

impl AtomicPublisher for MarkerObservingPublisher {
    fn publish_no_replace(
        &self,
        provisional: &Path,
        destination: &Path,
    ) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
        self.marker_present_at_entry
            .set(provisional.join("INCOMPLETE").is_file());
        if self.fail {
            UnsupportedAtomicPublisher.publish_no_replace(provisional, destination)
        } else {
            LinuxAtomicPublisher.publish_no_replace(provisional, destination)
        }
    }
}

struct MarkerRemovingFailurePublisher;

impl AtomicPublisher for MarkerRemovingFailurePublisher {
    fn publish_no_replace(
        &self,
        provisional: &Path,
        destination: &Path,
    ) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
        fs::remove_file(provisional.join("INCOMPLETE")).expect("remove marker in failure double");
        UnsupportedAtomicPublisher.publish_no_replace(provisional, destination)
    }
}

struct ReadOnlyDestinationPublisher;

impl AtomicPublisher for ReadOnlyDestinationPublisher {
    fn publish_no_replace(
        &self,
        provisional: &Path,
        destination: &Path,
    ) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
        LinuxAtomicPublisher.publish_no_replace(provisional, destination)?;
        let mut permissions = fs::metadata(destination)
            .expect("destination metadata")
            .permissions();
        permissions.set_mode(0o555);
        fs::set_permissions(destination, permissions).expect("read-only destination");
        Ok(())
    }
}
