use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::results::validate_effective_config_against_manifest_bytes;
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn runtime_rejects_truncation_unknown_class_absent_id_and_wrong_artifact_kind() {
    let destination = temporary("effective-cross-links");
    let config = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/reference_reduction/general.toml");
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run"])
        .arg(config)
        .arg("--output")
        .arg(&destination)
        .output()
        .expect("real execution");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let effective_bytes = fs::read(&effective_path).expect("effective config");
    let manifest_bytes = fs::read(&manifest_path).expect("manifest");
    validate_effective_config_against_manifest_bytes(
        &effective_path,
        &effective_bytes,
        &manifest_path,
        &manifest_bytes,
    )
    .expect("valid producer documents");

    let mut truncated = effective_bytes.clone();
    truncated.truncate(truncated.len() / 2);
    assert_rejected(&effective_path, &truncated, &manifest_path, &manifest_bytes);

    let mut unknown_class = parse(&effective_bytes);
    unknown_class["effective_config_kind"] = json!("unknown");
    assert_rejected_json(
        &effective_path,
        &unknown_class,
        &manifest_path,
        &manifest_bytes,
    );

    let mut absent_id = parse(&effective_bytes);
    absent_id["reference"]["fasta_artifact_id"] = json!("missing-artifact");
    assert_rejected_json(&effective_path, &absent_id, &manifest_path, &manifest_bytes);

    let manifest = parse(&manifest_bytes);
    let config_id = manifest["configuration_input"]["artifact_id"]
        .as_str()
        .expect("config artifact ID");
    let mut wrong_kind = parse(&effective_bytes);
    wrong_kind["reference"]["fasta_artifact_id"] = json!(config_id);
    assert_rejected_json(
        &effective_path,
        &wrong_kind,
        &manifest_path,
        &manifest_bytes,
    );

    fs::remove_dir_all(destination).expect("cleanup");
}

#[test]
fn every_reduction_field_and_exact_reference_set_are_cross_checked() {
    let destination = run_fixture(
        "reduction-representation",
        "tests/fixtures/reference_reduction/general.toml",
    );
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let effective_bytes = fs::read(&effective_path).expect("effective config");
    let manifest_bytes = fs::read(&manifest_path).expect("manifest");

    for (name, pointer, replacement) in [
        (
            "region-size",
            "/reduced_references/0/region_size_bases",
            json!(2),
        ),
        (
            "target-mode",
            "/reduced_references/0/target/mode",
            json!("reduction_percent"),
        ),
        (
            "target-declared",
            "/reduced_references/0/target/declared",
            json!("99"),
        ),
        (
            "target-bases",
            "/reduced_references/0/target/bases",
            json!(12),
        ),
        (
            "selection-subseed",
            "/reduced_references/0/selection_subseed",
            json!("2".repeat(64)),
        ),
        (
            "selection-order",
            "/reduced_references/0/selection_order",
            json!([99]),
        ),
    ] {
        let mut changed = parse(&effective_bytes);
        *changed.pointer_mut(pointer).expect(name) = replacement;
        assert_rejected_json(&effective_path, &changed, &manifest_path, &manifest_bytes);
    }

    for (name, pointer, replacement, expected_schema) in [
        (
            "eligible-bases",
            "/reference_reductions/0/eligible_bases",
            json!(19),
            "effective-config-v3.schema.json",
        ),
        (
            "root-seed",
            "/reference_reductions/0/root_seed",
            json!("2".repeat(64)),
            "effective-config-v3.schema.json",
        ),
        (
            "map-document-kind",
            "/reference_reductions/0/map_document_kind",
            json!("future_map"),
            "manifest-general-v3.schema.json",
        ),
        (
            "map-schema-version",
            "/reference_reductions/0/map_schema_version",
            json!("2"),
            "manifest-general-v3.schema.json",
        ),
    ] {
        let mut changed = parse(&manifest_bytes);
        *changed.pointer_mut(pointer).expect(name) = replacement;
        assert_rejected_manifest_json(
            &effective_path,
            &parse(&effective_bytes),
            &manifest_path,
            &changed,
            expected_schema,
        );
    }

    let effective = parse(&effective_bytes);
    let mut orphan_manifest = parse(&manifest_bytes);
    orphan_manifest["references"]
        .as_array_mut()
        .expect("references")
        .push(json!({
            "id": "orphan-reference",
            "fasta_artifact_id": effective["reference"]["fasta_artifact_id"]
        }));
    assert_rejected_json(
        &effective_path,
        &effective,
        &manifest_path,
        &serde_json::to_vec(&orphan_manifest).expect("manifest JSON"),
    );
    fs::remove_dir_all(destination).expect("cleanup");
}

#[test]
fn reference_composition_and_total_size_are_recomputed() {
    let destination = run_fixture(
        "reference-invariants",
        "tests/fixtures/config/general-seeded.toml",
    );
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let effective_bytes = fs::read(&effective_path).expect("effective config");
    let manifest_bytes = fs::read(&manifest_path).expect("manifest");

    let mut wrong_sum = parse(&effective_bytes);
    wrong_sum["reference"]["composition"]["A"] = json!(0.9);
    assert_rejected_json(&effective_path, &wrong_sum, &manifest_path, &manifest_bytes);

    let mut wrong_total = parse(&effective_bytes);
    let total = wrong_total["reference"]["total_genome_size"]
        .as_u64()
        .expect("total");
    wrong_total["reference"]["total_genome_size"] = json!(total + 1);
    assert_rejected_json(
        &effective_path,
        &wrong_total,
        &manifest_path,
        &manifest_bytes,
    );
    fs::remove_dir_all(destination).expect("cleanup");
}

#[test]
fn runtime_rejects_effective_initial_population_ids_that_violate_the_typed_contract() {
    let destination = run_fixture(
        "effective-initial-population-ids",
        "tests/fixtures/config/general-seeded.toml",
    );
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let effective = parse(&fs::read(&effective_path).expect("effective config"));
    let manifest = fs::read(&manifest_path).expect("manifest");

    for (population_id, individual_id) in [
        ("invalid population", "IND_valid"),
        ("P1", "invalid individual"),
    ] {
        let mut changed = effective.clone();
        changed["initial_populations"] = json!([{
            "population_id": population_id,
            "individual_count": 1,
            "individual_ids": [individual_id],
            "ploidy": 2,
            "genotype_model": {"kind": "hwe"}
        }]);
        assert_rejected_json(&effective_path, &changed, &manifest_path, &manifest);
    }
    fs::remove_dir_all(destination).expect("cleanup");
}

#[test]
fn runtime_rejects_corrupt_effective_initial_population_cardinality_duplicates_and_f_model() {
    let destination = run_fixture(
        "effective-initial-population-invariants",
        "tests/fixtures/config/general-seeded.toml",
    );
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let effective = parse(&fs::read(&effective_path).expect("effective config"));
    let manifest = fs::read(&manifest_path).expect("manifest");

    for population in [
        json!({
            "population_id": "P1",
            "individual_count": 2,
            "individual_ids": ["IND_1"],
            "ploidy": 2,
            "genotype_model": {"kind": "hwe"}
        }),
        json!({
            "population_id": "P1",
            "individual_count": 2,
            "individual_ids": ["IND_1", "IND_1"],
            "ploidy": 2,
            "genotype_model": {"kind": "hwe"}
        }),
        json!({
            "population_id": "P1",
            "individual_count": 2,
            "individual_ids": ["IND_1", "IND_2"],
            "ploidy": 3,
            "genotype_model": {"kind": "f_diploid", "inbreeding_coefficient": 0.25}
        }),
    ] {
        let mut changed = effective.clone();
        changed["initial_populations"] = json!([population]);
        assert_rejected_json(&effective_path, &changed, &manifest_path, &manifest);
    }
    fs::remove_dir_all(destination).expect("cleanup");
}

fn assert_rejected_json(
    effective_path: &Path,
    document: &Value,
    manifest_path: &Path,
    manifest_bytes: &[u8],
) {
    assert_rejected(
        effective_path,
        &serde_json::to_vec(document).expect("JSON"),
        manifest_path,
        manifest_bytes,
    );
}

fn assert_rejected_manifest_json(
    effective_path: &Path,
    effective: &Value,
    manifest_path: &Path,
    manifest: &Value,
    expected_schema: &str,
) {
    let error = validate_effective_config_against_manifest_bytes(
        effective_path,
        &serde_json::to_vec(effective).expect("effective JSON"),
        manifest_path,
        &serde_json::to_vec(manifest).expect("manifest JSON"),
    )
    .expect_err("invalid effective/manifest pair");
    assert_eq!(error.code(), "GPS-FMT-006");
    assert_eq!(error.category(), "FMT");
    assert_eq!(error.exit_code(), 1);
    assert_eq!(
        error.context().get("schema").map(String::as_str),
        Some(expected_schema)
    );
}

fn assert_rejected(
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) {
    let error = validate_effective_config_against_manifest_bytes(
        effective_path,
        effective_bytes,
        manifest_path,
        manifest_bytes,
    )
    .expect_err("invalid effective/manifest pair");
    assert_eq!(error.code(), "GPS-FMT-006");
    assert_eq!(error.category(), "FMT");
    assert_eq!(error.exit_code(), 1);
    assert_eq!(
        error.context().get("schema").map(String::as_str),
        Some("effective-config-v3.schema.json")
    );
}

fn parse(bytes: &[u8]) -> Value {
    serde_json::from_slice(bytes).expect("JSON fixture")
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f14-runtime-effective-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn run_fixture(label: &str, fixture: &str) -> PathBuf {
    let destination = temporary(label);
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run", fixture, "--output"])
        .arg(&destination)
        .output()
        .expect("real execution");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    destination
}
