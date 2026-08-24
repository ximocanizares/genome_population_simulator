use std::fs;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_fasta::{SequenceExpectation, validate_supplied_fasta};
use genome_population_simulator::results::validate_general_manifest_bytes;
use genome_population_simulator::truth_validation::validate_and_reconstruct_truth_jsonl;
use serde_json::json;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f2-4-4-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

#[test]
fn publication_catalogs_truth_with_its_complete_identity_and_a_broken_link_is_rejected() {
    let root = temporary();
    fs::create_dir_all(&root).expect("root");
    let config = root.join("config.toml");
    fs::write(
        &config,
        r#"config_kind = "general"
schema_version = "1"
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
kind = "simulated"
id = "ref-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 16
[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["IND_01"]
ploidy = 2
genotype_model = "hwe"
"#,
    )
    .expect("configuration");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("effective configuration"),
        &destination,
    )
    .expect("publication");

    let manifest_path = destination.join("manifest.json");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("manifest bytes"))
            .expect("manifest JSON");
    let truth = manifest["artifacts"]
        .as_array()
        .expect("catalog")
        .iter()
        .find(|artifact| artifact["id"] == "truth-jsonl")
        .expect("truth catalog entry");
    assert_eq!(truth["kind"], "truth_jsonl");
    assert_eq!(truth["path"], "outputs/truth.jsonl");
    let observed = identify_file(
        &destination.join("outputs/truth.jsonl"),
        PathBuf::from("outputs/truth.jsonl"),
    )
    .expect("truth identity");
    assert_eq!(truth["size_bytes"], json!(observed.size_bytes()));
    assert_eq!(truth["checksum"]["algorithm"], "sha256");
    assert_eq!(truth["checksum"]["digest"], observed.digest());

    let reference = validate_supplied_fasta(
        &destination.join("resources/references/ref-main.fa"),
        "ref-main",
        &[SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 16,
        }],
    )
    .expect("published authoritative reference");
    let truth_file = fs::File::open(destination.join("outputs/truth.jsonl")).expect("truth file");
    validate_and_reconstruct_truth_jsonl(&mut BufReader::new(truth_file), reference.reference())
        .expect("truth remains independently reconstructible");

    manifest["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .retain(|artifact| artifact["id"] != "truth-jsonl");
    assert!(
        validate_general_manifest_bytes(
            &manifest_path,
            &serde_json::to_vec(&manifest).expect("broken manifest")
        )
        .is_err()
    );

    fs::remove_dir_all(root).expect("cleanup");
}
