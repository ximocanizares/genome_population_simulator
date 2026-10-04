//! F2.5.5: publicación local VCF, catálogo y modos físicos.

use std::fs;
use std::io::{self, Read, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::{
    ExecutionAuthorizations, execute_validated_config, execute_validated_config_with_authorizations,
};
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
use genome_population_simulator::vcf_semantic_validation::validate_published_vcf_semantics;
use serde_json::Value;
use sha2::{Digest, Sha256};

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
        "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f\"\n[reference]\nkind = \"simulated\"\nid = \"ref-main\"\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 16\n[population_snvs]\nmode = \"generate\"\n[population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 1\nspatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n[[population_snvs.generate.populations]]\npopulation_id = \"P1\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nindividual_ids = [\"IND_01\"]\nploidy = 2\ngenotype_model = \"hwe\"\n{vcf}"
    )
}

fn configuration_without_populations(mode: &str) -> String {
    format!(
        "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f\"\n[reference]\nkind = \"simulated\"\nid = \"ref-main\"\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 16\n[vcf_output]\nmode = \"{mode}\"\n"
    )
}

fn publish(mode: Option<&str>) -> (PathBuf, Value) {
    publish_configuration(configuration(mode))
}

fn publish_zero_variant_vcf(mode: &str) -> (PathBuf, Value) {
    let root = temporary("zero-variant-publication");
    fs::create_dir_all(&root).expect("root");
    let config = root.join("config.toml");
    fs::write(
        &config,
        configuration(Some(mode)).replacen("locus_count = 1", "locus_count = 0", 1),
    )
    .expect("config");
    let destination = root.join("published");
    execute_validated_config_with_authorizations(
        &load_and_validate_config(&config).expect("valid config"),
        &destination,
        ExecutionAuthorizations {
            allow_empty_snv_matrix: true,
            allow_incompatible_genotypes: false,
        },
    )
    .expect("publication");
    let manifest =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
            .expect("manifest JSON");
    (root, manifest)
}

fn publish_configuration(configuration: String) -> (PathBuf, Value) {
    let root = temporary("publication");
    fs::create_dir_all(&root).expect("root");
    let config = root.join("config.toml");
    fs::write(&config, configuration).expect("config");
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
fn requested_text_vcf_without_populations_rejects_before_destination_or_manifest() {
    let root = temporary("vcf-without-samples");
    fs::create_dir_all(&root).expect("root");
    let config = root.join("config.toml");
    fs::write(&config, configuration_without_populations("text")).expect("config");
    let destination = root.join("published");

    let error = execute_validated_config(
        &load_and_validate_config(&config).expect("valid declaration"),
        &destination,
    )
    .expect_err("a VCF without truth samples must fail during preflight");

    assert_eq!(error.code(), "GPS-CONF-008");
    assert_eq!(
        error.context().get("component").map(String::as_str),
        Some("vcf_output,initial_populations")
    );
    assert!(!destination.exists());
    assert_eq!(
        fs::read_dir(&root)
            .expect("read root")
            .filter_map(Result::ok)
            .filter(|entry| entry
                .file_name()
                .to_string_lossy()
                .contains("gps-incomplete"))
            .count(),
        0
    );
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
    validate_published_vcf_semantics(
        &destination.join("resources/references/ref-main.fa"),
        &destination.join("outputs/truth.jsonl"),
        &vcf,
    )
    .expect("the published text VCF is semantically compatible with truth");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn fixed_fixture_keeps_the_text_vcf_byte_and_semantic_contract() {
    // SHA-256 congelado de la salida de HEAD anterior a PERF-VCF.1, obtenida
    // con esta configuración y revisada antes de compactar la ruta VCF.
    const PRE_COMPACTION_VCF_SHA256: &str =
        "5ce67541a3fe9a3e4edc1dc734debbd9bd1a920d3150f5000a32669b0c464253";
    let (root, _) = publish(Some("text"));
    let destination = root.join("published");
    let vcf = destination.join("outputs/truth.vcf");

    assert_eq!(
        format!("{:x}", Sha256::digest(fs::read(&vcf).expect("VCF"))),
        PRE_COMPACTION_VCF_SHA256,
        "la ruta VCF compacta conserva los bytes de la salida anterior"
    );
    validate_published_vcf_semantics(
        &destination.join("resources/references/ref-main.fa"),
        &destination.join("outputs/truth.jsonl"),
        &vcf,
    )
    .expect("VCF compatible con truth");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn frozen_pre_compaction_vcf_fixture_covers_snp_insertion_deletion_samples_phase_gt_and_gps_case() {
    // SHA-256 congelado de la salida anterior, no calculado por producción.
    const PRE_COMPACTION_VCF_SHA256: &str =
        "e2303de9754ac9cae96be181524fea2faba382690507c0c9cf9db7aa026cc4cf";
    let header = VcfHeaderContract::new(
        header_identity(),
        vec![VcfContig::new("chr1", 16).expect("contig")],
        vec!["IND_A".to_owned(), "IND_B".to_owned()],
    )
    .expect("header");
    let records = vec![
        golden_grouped_record(1, "A", "G", "CASE_SNP", vec![1, 0], vec![0, 1]),
        golden_grouped_record(3, "T", "TAA", "CASE_INS", vec![0, 1], vec![1, 0]),
        golden_grouped_record(5, "CG", "C", "CASE_DEL", vec![1, 0], vec![0, 1]),
    ];
    let bytes = write_vcf_records(Vec::new(), &header, &records).expect("VCF bytes");

    assert_eq!(
        format!("{:x}", Sha256::digest(&bytes)),
        PRE_COMPACTION_VCF_SHA256
    );
    let text = String::from_utf8(bytes).expect("text VCF");
    assert!(text.contains("\tA\tG\t.\tPASS\tGPS_CASE=CASE_SNP\tGT\t1|0\t0|1\n"));
    assert!(text.contains("\tT\tTAA\t.\tPASS\tGPS_CASE=CASE_INS\tGT\t0|1\t1|0\n"));
    assert!(text.contains("\tCG\tC\t.\tPASS\tGPS_CASE=CASE_DEL\tGT\t1|0\t0|1\n"));
}

#[test]
fn zero_variant_publication_is_complete_in_text_and_bgzf_modes() {
    for mode in ["text", "bgzf"] {
        let (root, manifest) = publish_zero_variant_vcf(mode);
        let destination = root.join("published");
        let relative_vcf = if mode == "text" {
            "outputs/truth.vcf"
        } else {
            "outputs/truth.vcf.gz"
        };
        let vcf = destination.join(relative_vcf);
        let text = if mode == "text" {
            fs::read_to_string(&vcf).expect("text VCF")
        } else {
            let mut reader = noodles_bgzf::io::Reader::new(fs::File::open(&vcf).expect("BGZF VCF"));
            let mut text = String::new();
            reader.read_to_string(&mut text).expect("read BGZF VCF");
            text
        };
        assert!(text.starts_with("##fileformat=VCFv4.5\n"), "{mode}");
        assert_eq!(
            text.lines().filter(|line| !line.starts_with('#')).count(),
            0,
            "{mode} must not emit variant records"
        );
        assert_eq!(manifest["vcf_outputs"][0]["mode"], mode, "{mode}");
        if mode == "bgzf" {
            assert!(
                destination.join("outputs/truth.vcf.gz.tbi").exists(),
                "{mode}"
            );
        }
        fs::remove_dir_all(root).expect("cleanup");
    }
}

fn golden_grouped_record(
    position: u64,
    reference: &str,
    alternate: &str,
    case: &str,
    first_sample: Vec<u32>,
    second_sample: Vec<u32>,
) -> GroupedVcfRecord {
    GroupedVcfRecord {
        sequence_id: "chr1".to_owned(),
        position,
        reference: reference.to_owned(),
        alternates: vec![alternate.to_owned()],
        samples: vec![
            (
                "IND_A".to_owned(),
                PhasedGenotype {
                    alleles: first_sample,
                    allele_dosages: vec![1, 1],
                },
            ),
            (
                "IND_B".to_owned(),
                PhasedGenotype {
                    alleles: second_sample,
                    allele_dosages: vec![1, 1],
                },
            ),
        ],
        gps_cases: vec![case.to_owned()],
        alternate_cases: vec![vec![case.to_owned()]],
        alternate_repeats: vec![Vec::new()],
    }
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
            None,
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
        None,
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
    validate_published_vcf_semantics(
        &destination.join("resources/references/ref-main.fa"),
        &destination.join("outputs/truth.jsonl"),
        &vcf,
    )
    .expect("the published BGZF VCF is semantically compatible with truth");
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
        alternate_repeats: vec![Vec::new()],
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
