//! F2.5.5: publicación local VCF, catálogo y modos físicos.

use std::fs;
use std::io::{self, Read, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::{ExecutionId, ReferenceId};
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::results::{
    validate_effective_config_against_manifest_bytes, validate_published_resources,
};
use genome_population_simulator::vcf_contract::{
    TBI_MAXIMUM_POSITION_EXCLUSIVE, VcfContig, VcfHeaderContract, VcfHeaderIdentity,
    write_vcf_index, write_vcf_records,
};
use genome_population_simulator::vcf_grouping::{GroupedVcfRecord, PhasedGenotype};
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn header_identity() -> VcfHeaderIdentity {
    VcfHeaderIdentity::new(
        ReferenceId::parse("ref-main", "VCF publication test").expect("reference ID"),
        ExecutionId::parse("run-main", "VCF publication test").expect("execution ID"),
    )
}

struct FailingWriter;

impl Write for FailingWriter {
    fn write(&mut self, _: &[u8]) -> io::Result<usize> {
        Err(io::Error::other("injected index failure"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f2-5-5-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn configuration(mode: Option<&str>) -> String {
    let vcf = mode.map_or(String::new(), |mode| {
        format!("\n[vcf_output]\nmode = \"{mode}\"\n")
    });
    format!(
        "config_kind = \"general\"\nschema_version = \"1\"\nseed = \"000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f\"\n[reference]\nkind = \"simulated\"\nid = \"ref-main\"\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 16\n{vcf}"
    )
}

fn publish(mode: Option<&str>) -> (PathBuf, Value) {
    let root = temporary("publication");
    fs::create_dir_all(&root).expect("root");
    let config = root.join("config.toml");
    fs::write(&config, configuration(mode)).expect("config");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("valid config"),
        &destination,
    )
    .expect("publication");
    let manifest =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
            .expect("manifest JSON");
    (root, manifest)
}

#[test]
fn vcf_is_disabled_by_default_and_not_cataloged() {
    let (root, manifest) = publish(None);
    assert!(
        manifest["vcf_outputs"]
            .as_array()
            .expect("outputs")
            .is_empty()
    );
    assert!(!root.join("published/outputs/truth.vcf").exists());
    assert!(!root.join("published/outputs/truth.vcf.gz").exists());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn existing_destination_prevents_requested_vcf_publication() {
    let root = temporary("existing-destination");
    fs::create_dir_all(&root).expect("root");
    let config = root.join("config.toml");
    fs::write(&config, configuration(Some("bgzf"))).expect("config");
    let destination = root.join("published");
    fs::create_dir_all(&destination).expect("existing destination");
    fs::write(destination.join("sentinel"), b"keep").expect("sentinel");

    let error = execute_validated_config(
        &load_and_validate_config(&config).expect("valid config"),
        &destination,
    )
    .expect_err("an existing destination must never be overwritten");

    assert_eq!(error.code(), "GPS-IO-003");
    assert_eq!(
        fs::read(destination.join("sentinel")).expect("sentinel"),
        b"keep"
    );
    assert!(!destination.join("outputs/truth.vcf.gz").exists());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn explicit_text_vcf_is_complete_and_has_no_index() {
    let (root, manifest) = publish(Some("text"));
    let destination = root.join("published");
    let vcf = destination.join("outputs/truth.vcf");
    let text = fs::read_to_string(&vcf).expect("text VCF");
    assert!(text.starts_with("##fileformat=VCFv4.5\n"));
    assert!(text.contains("##INFO=<ID=GPS_CASE,Number=.,Type=String,"));
    assert!(!destination.join("outputs/truth.vcf.gz.tbi").exists());
    assert!(!destination.join("outputs/truth.vcf.gz.csi").exists());
    let output = &manifest["vcf_outputs"][0];
    assert_eq!(output["mode"], "text");
    assert!(text.contains(&format!(
        "##reference={}\n",
        output["reference_id"].as_str().expect("reference ID")
    )));
    assert!(text.contains(&format!(
        "##gps_execution_id={}\n",
        manifest["execution_id"].as_str().expect("execution ID")
    )));
    assert!(text.contains("##gps_qual=.\n"));
    assert!(text.contains("##FILTER=<ID=PASS,Description=\"All filters passed\">\n"));
    assert!(output["index_kind"].is_null());
    assert!(output["index_artifact_id"].is_null());
    let artifact = manifest["artifacts"]
        .as_array()
        .expect("catalog")
        .iter()
        .find(|artifact| artifact["id"] == "truth-vcf")
        .expect("VCF catalog entry");
    let identity = identify_file(&vcf, PathBuf::from("outputs/truth.vcf")).expect("SHA");
    assert_eq!(artifact["size_bytes"], identity.size_bytes());
    assert_eq!(artifact["checksum"]["digest"], identity.digest());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn effective_vcf_request_manifest_and_catalog_are_bound_bidirectionally() {
    let (root, manifest) = publish(Some("text"));
    let destination = root.join("published");
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let effective = fs::read(&effective_path).expect("effective");

    let mut missing_catalog = manifest.clone();
    missing_catalog["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .retain(|artifact| artifact["id"] != "truth-vcf");
    assert!(
        validate_effective_config_against_manifest_bytes(
            &effective_path,
            &effective,
            &manifest_path,
            &serde_json::to_vec(&missing_catalog).expect("manifest JSON"),
        )
        .is_err()
    );

    let mut missing_link = manifest.clone();
    missing_link["vcf_outputs"] = serde_json::json!([]);
    assert!(
        validate_effective_config_against_manifest_bytes(
            &effective_path,
            &effective,
            &manifest_path,
            &serde_json::to_vec(&missing_link).expect("manifest JSON"),
        )
        .is_err()
    );

    let mut other_valid_reference = manifest.clone();
    let source_fasta = other_valid_reference["references"][0]["fasta_artifact_id"].clone();
    other_valid_reference["references"]
        .as_array_mut()
        .expect("references")
        .push(serde_json::json!({
            "id": "ref-alternate",
            "fasta_artifact_id": source_fasta,
        }));
    other_valid_reference["vcf_outputs"][0]["reference_id"] =
        Value::String("ref-alternate".to_owned());
    assert!(
        validate_effective_config_against_manifest_bytes(
            &effective_path,
            &effective,
            &manifest_path,
            &serde_json::to_vec(&other_valid_reference).expect("manifest JSON"),
        )
        .is_err()
    );

    let mut wrong_route = manifest;
    wrong_route["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .iter_mut()
        .find(|artifact| artifact["id"] == "truth-vcf")
        .expect("VCF artifact")["path"] = Value::String("outputs/not-truth.vcf".to_owned());
    assert!(
        validate_published_resources(
            &destination,
            &effective_path,
            &effective,
            &manifest_path,
            &serde_json::to_vec(&wrong_route).expect("manifest JSON"),
        )
        .is_err()
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn toml_without_vcf_cannot_claim_an_effective_vcf_request() {
    let (root, manifest) = publish(None);
    let destination = root.join("published");
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let mut effective: Value =
        serde_json::from_slice(&fs::read(&effective_path).expect("effective")).expect("JSON");
    effective["vcf_output"] = serde_json::json!({"mode": "text"});
    let error = validate_published_resources(
        &destination,
        &effective_path,
        &serde_json::to_vec(&effective).expect("effective JSON"),
        &manifest_path,
        &serde_json::to_vec(&manifest).expect("manifest JSON"),
    )
    .expect_err("the copied TOML keeps VCF disabled");
    assert_eq!(error.code(), "GPS-FMT-006");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn explicit_bgzf_vcf_has_one_readable_tbi_and_cataloged_sha() {
    let (root, manifest) = publish(Some("bgzf"));
    let destination = root.join("published");
    let vcf = destination.join("outputs/truth.vcf.gz");
    let index = destination.join("outputs/truth.vcf.gz.tbi");
    assert!(vcf.is_file());
    assert!(index.is_file());
    assert!(!destination.join("outputs/truth.vcf.gz.csi").exists());
    let output = &manifest["vcf_outputs"][0];
    let mut reader = noodles_bgzf::io::Reader::new(fs::File::open(&vcf).expect("BGZF"));
    let mut text = String::new();
    reader.read_to_string(&mut text).expect("read header");
    assert!(text.starts_with("##fileformat=VCFv4.5\n"));
    assert!(text.contains(&format!(
        "##reference={}\n",
        output["reference_id"].as_str().expect("reference ID")
    )));
    assert!(text.contains(&format!(
        "##gps_execution_id={}\n",
        manifest["execution_id"].as_str().expect("execution ID")
    )));
    assert!(text.contains("##gps_qual=.\n"));
    assert!(text.contains("##FILTER=<ID=PASS,Description=\"All filters passed\">\n"));
    noodles_tabix::io::Reader::new(fs::File::open(&index).expect("TBI"))
        .read_index()
        .expect("readable TBI");
    assert_eq!(output["mode"], "bgzf");
    assert_eq!(output["index_kind"], "tbi");
    let vcf_identity = identify_file(&vcf, PathBuf::from("outputs/truth.vcf.gz")).expect("VCF SHA");
    let index_identity =
        identify_file(&index, PathBuf::from("outputs/truth.vcf.gz.tbi")).expect("TBI SHA");
    for (id, identity) in [
        ("truth-vcf", vcf_identity),
        ("truth-vcf-index", index_identity),
    ] {
        let artifact = manifest["artifacts"]
            .as_array()
            .expect("catalog")
            .iter()
            .find(|artifact| artifact["id"] == id)
            .expect("catalog entry");
        assert_eq!(artifact["size_bytes"], identity.size_bytes());
        assert_eq!(artifact["checksum"]["digest"], identity.digest());
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn csi_indexes_real_records_and_is_queryable_by_region() {
    let root = temporary("csi");
    fs::create_dir_all(&root).expect("root");
    let header = VcfHeaderContract::new(
        header_identity(),
        vec![VcfContig::new("chr1", TBI_MAXIMUM_POSITION_EXCLUSIVE).expect("contig")],
        vec!["sample-1".to_owned()],
    )
    .expect("header");
    let records = [GroupedVcfRecord {
        sequence_id: "chr1".to_owned(),
        position: 4,
        reference: "A".to_owned(),
        alternates: vec!["C".to_owned()],
        samples: vec![(
            "sample-1".to_owned(),
            PhasedGenotype {
                alleles: vec![1, 0],
                allele_dosages: vec![1, 1],
            },
        )],
        gps_cases: vec!["CASE_1".to_owned()],
        alternate_cases: vec![vec!["CASE_1".to_owned()]],
    }];
    let vcf = root.join("truth.vcf.gz");
    let writer = noodles_bgzf::io::Writer::new(fs::File::create(&vcf).expect("VCF"));
    write_vcf_records(writer, &header, &records)
        .expect("write VCF")
        .finish()
        .expect("finish BGZF");
    let index = root.join("truth.vcf.gz.csi");
    let kind =
        write_vcf_index(&vcf, &header, fs::File::create(&index).expect("CSI")).expect("CSI index");
    assert_eq!(
        kind,
        genome_population_simulator::vcf_contract::VcfIndexKind::Csi
    );
    let index = noodles_csi::io::Reader::new(fs::File::open(&index).expect("CSI"))
        .read_index()
        .expect("read CSI");
    use noodles_csi::BinningIndex as _;
    let position = noodles_core::Position::try_from(5).expect("position");
    let chunks = index
        .query(0, (position..=position).into())
        .expect("regional CSI query");
    assert!(!chunks.is_empty());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn index_write_failure_is_reported_before_any_publication_step() {
    let root = temporary("index-failure");
    fs::create_dir_all(&root).expect("root");
    let header = VcfHeaderContract::new(
        header_identity(),
        vec![VcfContig::new("chr1", 12).expect("contig")],
        vec![],
    )
    .expect("header");
    let vcf = root.join("truth.vcf.gz");
    let writer = noodles_bgzf::io::Writer::new(fs::File::create(&vcf).expect("VCF"));
    write_vcf_records(writer, &header, &[])
        .expect("write VCF")
        .finish()
        .expect("finish BGZF");
    let error = write_vcf_index(&vcf, &header, FailingWriter)
        .expect_err("the failing index destination must be propagated");
    assert!(error.to_string().contains("injected index failure"));
    fs::remove_dir_all(root).expect("cleanup");
}
