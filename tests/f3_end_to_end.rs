//! Cierre PE actual: matriz válida, contaminación y hopping sin procedencia por pareja.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::fastq_validation::validate_final_paired_fastq;
use genome_population_simulator::results::{ArtifactRecord, validate_published_resources};
use genome_population_simulator::technical_identity::TechnicalMate;

const ROOT_SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn valid_matrix_contamination_and_hopping_publish_only_final_paired_fastq() {
    let root = temporary();
    fs::create_dir_all(&root).expect("temporary root");
    let fasta = root.join("reference.fa");
    let matrix = root.join("population-snvs.tsv");
    let config = root.join("paired-end.toml");
    fs::write(&fasta, format!(">chr-known\n{}\n", "A".repeat(300))).expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr-known\t75\tA\tT\tP1\t0.5\n",
    )
    .expect("matrix");
    fs::write(&config, configuration(&fasta, &matrix)).expect("configuration");

    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("effective configuration"),
        &destination,
    )
    .expect("current F4 PE publication");

    let manifest_path = destination.join("manifest.json");
    let effective_path = destination.join("effective-config.json");
    let manifest_bytes = fs::read(&manifest_path).expect("general manifest");
    let effective_bytes = fs::read(&effective_path).expect("effective configuration");
    validate_published_resources(
        &destination,
        &effective_path,
        &effective_bytes,
        &manifest_path,
        &manifest_bytes,
    )
    .expect("published PE resources are mutually verifiable");
    assert_external_manifest_schema(&manifest_path, "manifest-general-v4.schema.json");

    let event_path = destination.join("events/reads-known/manifest.json");
    assert_external_manifest_schema(&event_path, "manifest-event-v2.schema.json");
    let event: serde_json::Value =
        serde_json::from_slice(&fs::read(&event_path).expect("event manifest")).expect("JSON");
    assert_eq!(event["state"], "complete");
    assert_eq!(
        event["final_fastq_summary"].as_array().map(Vec::len),
        Some(2)
    );
    assert_eq!(event["demultiplexing"]["index_hopping_rate"], 1.0);
    let reassigned = event["demultiplexing"]["reassigned_pairs"]
        .as_u64()
        .expect("reassigned total");
    let published_pairs = event["final_fastq_summary"]
        .as_array()
        .expect("summary")
        .iter()
        .map(|entry| entry["pair_count"].as_u64().expect("pair count"))
        .sum::<u64>();
    assert_eq!(
        reassigned, published_pairs,
        "hopping at rate one loses no pairs"
    );
    assert_eq!(
        event["demultiplexing"]["reassigned_by_sample"]
            .as_array()
            .expect("aggregates")
            .iter()
            .map(|entry| entry["pair_count"].as_u64().expect("aggregate count"))
            .sum::<u64>(),
        reassigned
    );
    assert!(
        event["final_fastq_summary"]
            .as_array()
            .expect("summary")
            .iter()
            .any(|entry| entry["contaminant_pair_count"]
                .as_u64()
                .expect("contaminant count")
                > 0)
    );

    for summary in event["final_fastq_summary"].as_array().expect("summary") {
        let individual = summary["individual_id"].as_str().expect("individual");
        let outputs = event["outputs"].as_array().expect("outputs");
        let r1 = output_artifact(outputs, individual, TechnicalMate::R1);
        let r2 = output_artifact(outputs, individual, TechnicalMate::R2);
        validate_final_paired_fastq(
            &destination.join(&r1.path),
            &r1,
            &destination.join(&r2.path),
            &r2,
        )
        .expect("R1/R2 final paired FASTQ");
    }
    let serialized = serde_json::to_string(&event).expect("serialize manifest");
    for forbidden in ["PH_", "TP_", "technical_pair_id", "pair_ids"] {
        assert!(
            !serialized.contains(forbidden),
            "durable provenance: {forbidden}"
        );
    }
    fs::remove_dir_all(root).expect("cleanup");
}

fn configuration(fasta: &Path, matrix: &Path) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "2"
seed = "{ROOT_SEED}"
[population_snvs]
mode = "supplied"
input_path = "{}"
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["sample-A", "sample-B"]
ploidy = 2
genotype_model = "hwe"
[reference]
kind = "supplied"
id = "reference-known"
fasta = "{}"
[[reference.sequences]]
id = "chr-known"
class = "chromosome"
length = 300
[technical_event]
event_id = "reads-known"
individual_ids = ["sample-A", "sample-B"]
target_coverage_wgs = 1.0
contamination_min_percent = 50
contamination_max_percent = 50
contamination_recipient_ids = ["sample-A"]
index_hopping_rate = 1.0

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"

[[technical_event.libraries]]
individual_id = "sample-B"
i7 = "CGTACGTA"
i5 = "GCATGCAT"
"#,
        matrix.display(),
        fasta.display()
    )
}

fn output_artifact(
    outputs: &[serde_json::Value],
    individual_id: &str,
    mate: TechnicalMate,
) -> ArtifactRecord {
    let mate_number = match mate {
        TechnicalMate::R1 => 1,
        TechnicalMate::R2 => 2,
    };
    let output = outputs
        .iter()
        .find(|entry| {
            entry["individual_id"] == individual_id
                && entry["mate"] == mate_number
                && entry["state"] == "complete"
        })
        .expect("completed mate output");
    let artifact = &output["artifact"];
    ArtifactRecord {
        id: artifact["id"].as_str().expect("artifact ID").to_owned(),
        kind: artifact["kind"].as_str().expect("artifact kind").to_owned(),
        path: artifact["path"].as_str().expect("artifact path").to_owned(),
        size_bytes: artifact["size_bytes"].as_u64().expect("artifact size"),
        digest: artifact["checksum"]["digest"]
            .as_str()
            .expect("artifact digest")
            .to_owned(),
    }
}

fn assert_external_manifest_schema(path: &Path, schema: &str) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(path)
        .output()
        .expect("external JSON Schema validator");
    assert!(
        output.status.success(),
        "external validator for {path:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f3-end-to-end-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}
