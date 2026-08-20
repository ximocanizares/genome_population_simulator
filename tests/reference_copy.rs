use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::{ExecutionOutcome, execute_validated_config};
use genome_population_simulator::publication::ExecutionWorkspace;
use genome_population_simulator::reference::{ReferenceOrigin, SequenceClass};
use genome_population_simulator::reference_fasta::{
    SequenceExpectation, validate_supplied_fasta, validate_supplied_fasta_with_origin,
};

static COUNTER: AtomicU64 = AtomicU64::new(0);

#[test]
fn supplied_run_copies_the_validated_fasta_byte_for_byte() {
    let directory = TestDirectory::new("exact-copy");
    let fasta = directory.path().join("source.fa");
    let original = b">chr1 description\r\naCgTn\r\n";
    fs::write(&fasta, original).expect("write FASTA");
    let config = directory.path().join("config.toml");
    write_supplied_config(&config, "source.fa", 5);
    let effective = load_and_validate_config(&config).expect("valid configuration");
    let destination = directory.path().join("result");

    assert!(matches!(
        execute_validated_config(&effective, &destination).expect("published run"),
        ExecutionOutcome::Published { .. }
    ));
    let copy = destination.join("resources/references/ref-supplied.fa");
    assert_eq!(fs::read(&copy).expect("read copied FASTA"), original);
    assert_eq!(fs::read(&fasta).expect("read original FASTA"), original);
    assert!(destination.join("manifest.json").is_file());

    let effective_json: serde_json::Value = serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("effective config"),
    )
    .expect("effective JSON");
    let reference = &effective_json["reference"];
    assert_eq!(reference["source"], "supplied");
    assert_eq!(reference["composition"], serde_json::Value::Null);
    assert_eq!(reference["sequences"][0]["unknown_base_count"], 1);
    assert_eq!(reference["sequences"][0]["length"], 5);
    assert_eq!(reference["sequences"][0]["class"], "chromosome");

    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
            .expect("manifest JSON");
    assert_eq!(manifest["references"][0]["id"], reference["id"]);
    assert_eq!(manifest["references"][0]["fasta"], reference["fasta"]);
    let linked = destination.join(
        reference["fasta"]["path"]
            .as_str()
            .expect("relative FASTA link"),
    );
    assert_eq!(linked, copy);
    fs::write(&fasta, b">chr1 changed\nAAAAA\n").expect("mutate original after publication");
    assert_eq!(
        fs::read(&copy).expect("authoritative copy after source mutation"),
        original,
        "post-copy source changes cannot alter the authoritative resource"
    );
    let external = Command::new("sha256sum")
        .arg(&linked)
        .output()
        .expect("independent SHA-256 tool");
    assert!(external.status.success());
    let digest = String::from_utf8(external.stdout)
        .expect("sha256sum UTF-8")
        .split_whitespace()
        .next()
        .expect("digest")
        .to_owned();
    assert_eq!(reference["fasta"]["checksum"]["digest"], digest);
}

#[test]
fn invalid_supplied_fasta_is_rejected_before_creating_a_destination() {
    let directory = TestDirectory::new("invalid-preflight");
    fs::write(directory.path().join("source.fa"), b">chr1\nACR\n").expect("write FASTA");
    let config = directory.path().join("config.toml");
    write_supplied_config(&config, "source.fa", 3);
    let effective = load_and_validate_config(&config).expect("configuration syntax is valid");
    let destination = directory.path().join("result");

    let error = execute_validated_config(&effective, &destination)
        .expect_err("invalid FASTA must reject execution admission");
    assert_eq!(error.code(), "GPS-FMT-005");
    assert_eq!(error.exit_code(), 2);
    assert!(!destination.exists());
    assert_eq!(
        fs::read_dir(directory.path())
            .expect("read test directory")
            .filter_map(Result::ok)
            .filter(|entry| entry
                .file_name()
                .to_string_lossy()
                .contains("gps-incomplete"))
            .count(),
        0
    );
}

#[test]
fn a_source_change_between_validation_and_copy_cannot_be_cataloged() {
    let directory = TestDirectory::new("changed-source");
    let source = directory.path().join("source.fa");
    fs::write(&source, b">chr1\nACGT\n").expect("write FASTA");
    let expectation = [SequenceExpectation {
        id: "chr1",
        class: SequenceClass::Chromosome,
        length: 4,
    }];
    let validated =
        validate_supplied_fasta(&source, "ref-main", &expectation).expect("initial FASTA is valid");
    fs::write(&source, b">chr1\nACGA\n").expect("mutate source");
    let destination = directory.path().join("result");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");

    let error = workspace
        .copy_verified_artifact(
            "reference-fasta",
            "reference_fasta",
            &source,
            Path::new("resources/references/ref-main.fa"),
            validated.source_identity(),
        )
        .expect_err("changed source must not be cataloged");
    assert_eq!(error.code(), "GPS-INPUT-002");
    assert!(workspace.artifacts().is_empty());
    assert!(
        provisional_workspace(&destination)
            .join("INCOMPLETE")
            .is_file()
    );
    assert!(!destination.exists());
}

#[test]
fn authoritative_representation_is_loaded_from_the_copy_and_keeps_original_provenance() {
    let directory = TestDirectory::new("authoritative-copy");
    let source = directory.path().join("source.fa");
    let original = b">chr1 copied header\r\naCgTn\r\n";
    fs::write(&source, original).expect("write FASTA");
    let expectation = [SequenceExpectation {
        id: "chr1",
        class: SequenceClass::Chromosome,
        length: 5,
    }];
    let validated = validate_supplied_fasta(&source, "ref-main", &expectation)
        .expect("original FASTA passes preflight validation");
    let destination = directory.path().join("result");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let artifact = workspace
        .copy_verified_artifact(
            "reference-fasta",
            "reference_fasta",
            &source,
            Path::new("resources/references/ref-main.fa"),
            validated.source_identity(),
        )
        .expect("verified copy");

    fs::write(&source, b">chr1 mutated header\nAAAAA\n")
        .expect("mutate original before authoritative loading");
    let copy = provisional_workspace(&destination).join(&artifact.path);
    let authoritative =
        validate_supplied_fasta_with_origin(&copy, &source, "ref-main", &expectation)
            .expect("authoritative representation is loaded from copied bytes");

    assert_eq!(authoritative.source_identity().digest(), artifact.digest);
    assert_eq!(
        authoritative.source_identity().size_bytes(),
        artifact.size_bytes
    );
    assert_eq!(
        authoritative.reference().origin(),
        &ReferenceOrigin::Supplied {
            original_path: source.clone()
        },
        "provenance must remain the original input path"
    );
    let sequence = &authoritative.reference().sequences()[0];
    assert_eq!(sequence.header().raw_bytes(), b"chr1 copied header");
    assert_eq!(sequence.sequence().as_bytes(), b"ACGTN");
    assert_eq!(
        fs::read(&source).expect("mutated original"),
        b">chr1 mutated header\nAAAAA\n",
        "the original differs before the authoritative representation is consumed"
    );
    assert_eq!(fs::read(copy).expect("copied FASTA"), original);
}

fn provisional_workspace(destination: &Path) -> PathBuf {
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    let name = destination
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap();
    let prefix = format!(".{name}.gps-incomplete-");
    let mut matches = fs::read_dir(parent)
        .unwrap()
        .filter_map(Result::ok)
        .filter(|entry| entry.file_name().to_string_lossy().starts_with(&prefix))
        .map(|entry| entry.path())
        .collect::<Vec<_>>();
    matches.sort();
    assert_eq!(matches.len(), 1);
    matches.remove(0)
}

fn write_supplied_config(path: &Path, fasta: &str, length: u64) {
    fs::write(
        path,
        format!(
            r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[reference]
kind = "supplied"
id = "ref-supplied"
fasta = "{fasta}"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {length}
"#
        ),
    )
    .expect("write configuration");
}

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new(label: &str) -> Self {
        let counter = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "gps-reference-copy-{label}-{}-{counter}",
            std::process::id()
        ));
        fs::create_dir(&path).expect("create isolated test directory");
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
