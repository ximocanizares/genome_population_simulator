use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::{
    ExecutionAuthorizations, execute_validated_config, execute_validated_config_with_authorizations,
};
use genome_population_simulator::identity::{SequenceId, derive_matrix_variant_id};
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::population_snvs::{
    TSV_HEADER, parse_canonical_tsv, parse_metadata_json, validate_reusable_population_snvs_bytes,
    validate_reusable_population_snvs_with_hooks,
};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_generation::{
    BaseFrequency, GeneratedSequenceSpec, generate_reference,
};
use genome_population_simulator::results::validate_published_resources;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f2-2-7-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn general(population_snvs: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
{population_snvs}
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
length = 20
"#
    )
}

#[test]
fn generation_publishes_the_closed_canonical_pair_and_reuse_only_consumes_it() {
    let root = temporary();
    fs::create_dir_all(&root).expect("root");
    let generated_config = root.join("generated.toml");
    fs::write(
        &generated_config,
        general(
            r#"[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 2
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P2"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9"#,
        ),
    )
    .expect("config");
    let generated = root.join("generated");
    execute_validated_config(
        &load_and_validate_config(&generated_config).expect("effective"),
        &generated,
    )
    .expect("generated publication");

    let matrix_path = generated.join("population_snvs.tsv");
    let metadata_path = generated.join("population_snvs.metadata.json");
    let matrix = fs::read(&matrix_path).expect("matrix");
    let table = parse_canonical_tsv(&matrix).expect("canonical TSV");
    assert_eq!(table.rows().len(), 4);
    assert_eq!(
        matrix.split(|byte| *byte == b'\n').next().expect("header"),
        TSV_HEADER.join("\t").as_bytes()
    );
    let metadata =
        parse_metadata_json(&fs::read(&metadata_path).expect("metadata")).expect("closed metadata");
    let matrix_digest =
        identify_file(&matrix_path, PathBuf::from("population_snvs.tsv")).expect("digest");
    assert_eq!(
        metadata.document()["matrix"]["sha256"],
        matrix_digest.digest()
    );
    assert_eq!(metadata.document()["populations"][0]["population_id"], "P1");
    assert_eq!(metadata.document()["provenance"], "generated");

    let repeated = root.join("generated-repeat");
    execute_validated_config(
        &load_and_validate_config(&generated_config).expect("same effective"),
        &repeated,
    )
    .expect("repeated generated publication");
    assert_eq!(
        fs::read(repeated.join("population_snvs.tsv")).expect("repeated TSV"),
        matrix
    );
    assert_eq!(
        fs::read(repeated.join("population_snvs.metadata.json")).expect("repeated metadata"),
        fs::read(&metadata_path).expect("metadata")
    );

    let reference = generate_reference(
        "ref-main",
        RootSeed::from_bytes([
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31,
        ]),
        &[GeneratedSequenceSpec {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 20,
        }],
        &[
            BaseFrequency {
                base: b'A',
                probability: 0.25,
            },
            BaseFrequency {
                base: b'C',
                probability: 0.25,
            },
            BaseFrequency {
                base: b'G',
                probability: 0.25,
            },
            BaseFrequency {
                base: b'T',
                probability: 0.25,
            },
        ],
    )
    .expect("same reference");
    let metadata_original = fs::read(&metadata_path).expect("metadata snapshot");
    let matrix_original = fs::read(&matrix_path).expect("matrix snapshot");
    let expected_metadata_sha = identify_file(
        &metadata_path,
        PathBuf::from("population_snvs.metadata.json"),
    )
    .expect("original metadata SHA")
    .digest()
    .to_owned();
    let reference_sha = metadata.document()["reference_fasta"]["sha256"]
        .as_str()
        .expect("sha")
        .to_owned();
    let (_, _, snapshot_metadata_sha, _) = validate_reusable_population_snvs_with_hooks(
        &metadata_path,
        &reference,
        &reference_sha,
        1e-9,
        false,
        || {
            fs::write(&metadata_path, b"corrupted after metadata snapshot")
                .expect("mutate metadata")
        },
        || fs::write(&matrix_path, b"corrupted after matrix snapshot").expect("mutate matrix"),
    )
    .expect("validated snapshot must not reopen either source");
    assert_eq!(snapshot_metadata_sha, expected_metadata_sha);
    fs::write(&metadata_path, &metadata_original).expect("restore metadata");
    fs::write(&matrix_path, &matrix_original).expect("restore matrix");
    let reject = |metadata: &[u8], matrix: &[u8]| {
        validate_reusable_population_snvs_bytes(
            metadata,
            matrix,
            &reference,
            &reference_sha,
            1e-9,
            false,
        )
        .is_err()
    };
    let original_document: serde_json::Value =
        serde_json::from_slice(&metadata_original).expect("JSON");
    let mut changed_reference = original_document.clone();
    changed_reference["reference"]["canonical_fingerprint"] = serde_json::json!("0".repeat(64));
    assert!(reject(
        &serde_json::to_vec(&changed_reference).expect("JSON"),
        &matrix_original
    ));
    let mut changed_reference_fasta = original_document.clone();
    changed_reference_fasta["reference_fasta"]["sha256"] = serde_json::json!("0".repeat(64));
    assert!(reject(
        &serde_json::to_vec(&changed_reference_fasta).expect("JSON"),
        &matrix_original
    ));
    let mut changed_provenance = original_document.clone();
    changed_provenance["provenance"] = serde_json::json!("unknown");
    assert!(reject(
        &serde_json::to_vec(&changed_provenance).expect("JSON"),
        &matrix_original
    ));
    let mut changed_audit = original_document.clone();
    changed_audit["normalization_audit"][0]["original_sum"] = serde_json::json!(0.0);
    assert!(reject(
        &serde_json::to_vec(&changed_audit).expect("JSON"),
        &matrix_original
    ));
    let mut bogus_group = original_document.clone();
    let mut group = bogus_group["normalization_audit"][0].clone();
    group["population_id"] = serde_json::json!("P999");
    group["original_sum"] = serde_json::json!(0.0);
    group["original_delta"] = serde_json::json!(1.0);
    group["effective_sum"] = serde_json::json!(0.0);
    group["effective_delta"] = serde_json::json!(1.0);
    group["reference_frequency"] = serde_json::json!(1.0);
    for value in group["values"].as_array_mut().expect("values") {
        value["original"] = serde_json::json!(0.0);
        value["effective"] = serde_json::json!(0.0);
    }
    bogus_group["normalization_audit"]
        .as_array_mut()
        .expect("audit")
        .push(group);
    assert!(reject(
        &serde_json::to_vec(&bogus_group).expect("JSON"),
        &matrix_original
    ));
    let mut bogus_zero = original_document.clone();
    let group = &mut bogus_zero["normalization_audit"][0];
    let reference_base = group["ref"].as_str().expect("ref").as_bytes()[0];
    let used_alternates = group["values"]
        .as_array()
        .expect("values")
        .iter()
        .filter_map(|value| value["alt"].as_str())
        .collect::<Vec<_>>();
    let alternate = *b"ACGT"
        .iter()
        .find(|base| {
            **base != reference_base
                && !used_alternates.contains(&std::str::from_utf8(&[**base]).expect("base"))
        })
        .expect("unused alternate");
    let variant_id = derive_matrix_variant_id(
        reference.id(),
        &SequenceId::parse(group["contig"].as_str().expect("contig"), "test").expect("contig"),
        group["position"]
            .as_str()
            .expect("position")
            .parse()
            .expect("position"),
        reference_base,
        alternate,
    )
    .expect("canonical VariantId");
    group["values"]
        .as_array_mut()
        .expect("values")
        .push(serde_json::json!({
            "alt": char::from(alternate).to_string(),
            "variant_id": variant_id.as_str(),
            "source_id": "bogus-zero",
            "original": 0.0,
            "effective": 0.0,
        }));
    group["values"]
        .as_array_mut()
        .expect("values")
        .sort_by(|left, right| {
            (
                left["alt"].as_str(),
                left["variant_id"].as_str(),
                left["source_id"].as_str(),
            )
                .cmp(&(
                    right["alt"].as_str(),
                    right["variant_id"].as_str(),
                    right["source_id"].as_str(),
                ))
        });
    assert!(reject(
        &serde_json::to_vec(&bogus_zero).expect("JSON"),
        &matrix_original
    ));
    let mut changed_populations = original_document.clone();
    changed_populations["populations"]
        .as_array_mut()
        .expect("populations")
        .reverse();
    assert!(reject(
        &serde_json::to_vec(&changed_populations).expect("JSON"),
        &matrix_original
    ));
    let mut changed_af = matrix_original.clone();
    let first_row = changed_af
        .iter()
        .position(|byte| *byte == b'\n')
        .expect("header")
        + 1;
    let af = changed_af[first_row..]
        .iter()
        .enumerate()
        .filter(|(_, byte)| **byte == b'\t')
        .nth(5)
        .map(|(offset, _)| first_row + offset + 1)
        .expect("AF");
    changed_af[af] = if changed_af[af] == b'9' { b'8' } else { b'9' };
    assert!(reject(&metadata_original, &changed_af));
    let mut changed_id = matrix_original.clone();
    let variant_id = changed_id
        .iter()
        .position(|byte| *byte == b'\n')
        .expect("header")
        + 1;
    changed_id[variant_id] = b'x';
    assert!(reject(&metadata_original, &changed_id));
    let mut lines = matrix_original
        .split(|byte| *byte == b'\n')
        .collect::<Vec<_>>();
    lines.swap(1, 2);
    let mut changed_order = lines[..lines.len() - 1].join(&b'\n');
    changed_order.push(b'\n');
    assert!(reject(&metadata_original, &changed_order));

    let reuse_config = root.join("reuse.toml");
    fs::write(
        &reuse_config,
        general(&format!(
            "[population_snvs]\nmode = \"reuse\"\nmetadata_path = \"{}\"",
            metadata_path.display()
        )),
    )
    .expect("reuse config");
    let reused = root.join("reused");
    execute_validated_config(
        &load_and_validate_config(&reuse_config).expect("reuse effective"),
        &reused,
    )
    .expect("reuse publication");
    assert!(!reused.join("population_snvs.tsv").exists());
    assert!(!reused.join("population_snvs.metadata.json").exists());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(reused.join("manifest.json")).expect("manifest"))
            .expect("JSON");
    assert_eq!(manifest["population_snvs"]["mode"], "reuse");
    assert!(manifest["population_snvs"]["metadata_artifact_id"].is_null());
    let source_matrix_sha256 = identify_file(
        &generated.join("population_snvs.tsv"),
        PathBuf::from("population_snvs.tsv"),
    )
    .expect("source matrix identity")
    .digest()
    .to_owned();
    let truth_header: serde_json::Value = serde_json::from_slice(
        fs::read(reused.join("outputs/truth.jsonl"))
            .expect("truth")
            .split(|byte| *byte == b'\n')
            .next()
            .expect("truth header"),
    )
    .expect("truth header JSON");
    for digest in truth_header["source_digests"]
        .as_array()
        .expect("source digests")
        .iter()
        .chain(
            truth_header["group_digests"]
                .as_array()
                .expect("group digests"),
        )
        .filter(|digest| digest["source"]["kind"] == "matrix")
    {
        assert_eq!(digest["sha256"], source_matrix_sha256);
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn discarded_variant_requires_every_declared_population_in_reuse_and_physical_validation() {
    let root = temporary();
    fs::create_dir_all(&root).expect("root");
    let fasta = root.join("reference.fa");
    fs::write(&fasta, b">chr1\nAC\n").expect("FASTA");
    let input = root.join("all-zero.tsv");
    fs::write(
        &input,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\tsource_id\nchr1\t1\tA\tC\tP1\t0\tzero-p1\nchr1\t1\tA\tC\tP2\t0\tzero-p2\n",
    )
    .expect("all-zero matrix");
    let config = root.join("supplied.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "2"
[population_snvs]
mode = "supplied"
input_path = "{}"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 2
"#,
            input.display(),
            fasta.display()
        ),
    )
    .expect("config");
    let authorizations = ExecutionAuthorizations {
        allow_empty_snv_matrix: true,
        allow_incompatible_genotypes: false,
    };
    let published = root.join("published");
    execute_validated_config_with_authorizations(
        &load_and_validate_config(&config).expect("effective"),
        &published,
        authorizations,
    )
    .expect("authorized all-zero publication");
    let metadata_path = published.join("population_snvs.metadata.json");
    let mut metadata: serde_json::Value =
        serde_json::from_slice(&fs::read(&metadata_path).expect("metadata")).expect("JSON");
    assert_eq!(
        metadata["populations"]
            .as_array()
            .expect("populations")
            .len(),
        2
    );
    assert_eq!(
        metadata["discarded_variant_pairs"]
            .as_array()
            .expect("pairs")
            .len(),
        2
    );
    metadata["normalization_audit"]
        .as_array_mut()
        .expect("audit")
        .retain(|group| group["population_id"] != "P2");
    fs::write(
        &metadata_path,
        serde_json::to_vec_pretty(&metadata).expect("bad metadata"),
    )
    .expect("write bad metadata");
    let reuse = root.join("reuse.toml");
    fs::write(
        &reuse,
        format!(
            r#"config_kind = "general"
schema_version = "2"
[population_snvs]
mode = "reuse"
metadata_path = "{}"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 2
"#,
            metadata_path.display(),
            fasta.display()
        ),
    )
    .expect("reuse config");
    assert!(
        execute_validated_config_with_authorizations(
            &load_and_validate_config(&reuse).expect("reuse effective"),
            &root.join("reused"),
            authorizations,
        )
        .is_err(),
        "falta el par P2×VariantId descartado"
    );
    let effective_path = published.join("effective-config.json");
    let manifest_path = published.join("manifest.json");
    let metadata_identity = identify_file(
        &metadata_path,
        PathBuf::from("population_snvs.metadata.json"),
    )
    .expect("mutated metadata identity");
    let mut effective: serde_json::Value =
        serde_json::from_slice(&fs::read(&effective_path).expect("effective")).expect("JSON");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("manifest")).expect("JSON");
    let audit = metadata["normalization_audit"].clone();
    for document in [&mut effective, &mut manifest] {
        document["population_snvs"]["metadata_sha256"] =
            serde_json::json!(metadata_identity.digest());
        document["population_snvs"]["normalization_audit"] = audit.clone();
    }
    for artifact in manifest["artifacts"].as_array_mut().expect("catalog") {
        if artifact["id"] == "population-snvs-metadata" {
            artifact["size_bytes"] = serde_json::json!(metadata_identity.size_bytes());
            artifact["checksum"]["digest"] = serde_json::json!(metadata_identity.digest());
        }
    }
    let physical_error = validate_published_resources(
        &published,
        &effective_path,
        &serde_json::to_vec(&effective).expect("forged effective"),
        &manifest_path,
        &serde_json::to_vec(&manifest).expect("forged manifest"),
    )
    .expect_err("el validador físico debe rechazar metadata sin el grupo P2");
    assert!(
        physical_error.context()["detail"].contains("pares población"),
        "el rechazo físico debe proceder de la auditoría de descartes: {}",
        physical_error.render_human()
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn supplied_matrix_is_prepared_into_the_same_closed_resource_shape() {
    let root = temporary();
    fs::create_dir_all(&root).expect("root");
    let fasta = root.join("reference.fa");
    fs::write(&fasta, b">chr1\nACGTACGTACGTACGTACGT\n").expect("reference");
    let input = root.join("input.tsv");
    fs::write(
        &input,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\tsource_id\nchr1\t1\tA\tC\tP2\t0.2\tstudy-b\nchr1\t1\tA\tC\tP1\t0.1\tstudy-a\n",
    )
    .expect("input");
    let config = root.join("supplied.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[population_snvs]
mode = "supplied"
input_path = "{}"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 20
"#,
            input.display(),
            fasta.display()
        ),
    )
    .expect("config");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("effective"),
        &destination,
    )
    .expect("supplied publication");
    let matrix = fs::read(destination.join("population_snvs.tsv")).expect("matrix");
    let table = parse_canonical_tsv(&matrix).expect("canonical TSV");
    assert_eq!(table.rows()[0].field(5), "P1");
    assert!(!table.rows()[0].field(0).is_empty());
    let metadata = parse_metadata_json(
        &fs::read(destination.join("population_snvs.metadata.json")).expect("metadata"),
    )
    .expect("metadata");
    assert_eq!(metadata.document()["provenance"], "supplied");
    let source_matrix = matrix.clone();
    let source_metadata =
        fs::read(destination.join("population_snvs.metadata.json")).expect("metadata");
    let reuse_config = root.join("supplied-reuse.toml");
    fs::write(
        &reuse_config,
        format!(
            r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[population_snvs]
mode = "reuse"
metadata_path = "{}"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 20
"#,
            destination.join("population_snvs.metadata.json").display(),
            fasta.display()
        ),
    )
    .expect("reuse config");
    let reused = root.join("supplied-reused");
    execute_validated_config(
        &load_and_validate_config(&reuse_config).expect("reuse effective"),
        &reused,
    )
    .expect("reuse");
    assert_eq!(
        fs::read(destination.join("population_snvs.tsv")).expect("source unchanged"),
        source_matrix
    );
    assert_eq!(
        fs::read(destination.join("population_snvs.metadata.json")).expect("source unchanged"),
        source_metadata
    );
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let mut effective: serde_json::Value =
        serde_json::from_slice(&fs::read(&effective_path).expect("effective")).expect("JSON");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("manifest")).expect("JSON");
    let mut correlated_metadata: serde_json::Value = serde_json::from_slice(
        &fs::read(destination.join("population_snvs.metadata.json")).expect("metadata"),
    )
    .expect("JSON");
    let mut mutated_matrix = matrix.clone();
    let variant_id = mutated_matrix
        .iter()
        .position(|byte| *byte == b'\n')
        .expect("header")
        + 1;
    mutated_matrix[variant_id] = b'x';
    fs::write(destination.join("population_snvs.tsv"), mutated_matrix).expect("mutate matrix");
    correlated_metadata["matrix"]["sha256"] = serde_json::Value::String(
        identify_file(
            &destination.join("population_snvs.tsv"),
            PathBuf::from("population_snvs.tsv"),
        )
        .expect("digest")
        .digest()
        .to_owned(),
    );
    fs::write(
        destination.join("population_snvs.metadata.json"),
        serde_json::to_vec_pretty(&correlated_metadata).expect("JSON"),
    )
    .expect("mutate metadata");
    let matrix_identity = identify_file(
        &destination.join("population_snvs.tsv"),
        PathBuf::from("population_snvs.tsv"),
    )
    .expect("matrix identity");
    let metadata_identity = identify_file(
        &destination.join("population_snvs.metadata.json"),
        PathBuf::from("population_snvs.metadata.json"),
    )
    .expect("metadata identity");
    effective["population_snvs"]["metadata_sha256"] = serde_json::json!(metadata_identity.digest());
    manifest["population_snvs"]["metadata_sha256"] = serde_json::json!(metadata_identity.digest());
    for artifact in manifest["artifacts"].as_array_mut().expect("catalog") {
        let identity = match artifact["id"].as_str() {
            Some("population-snvs-tsv") => Some(&matrix_identity),
            Some("population-snvs-metadata") => Some(&metadata_identity),
            _ => None,
        };
        if let Some(identity) = identity {
            artifact["size_bytes"] = serde_json::json!(identity.size_bytes());
            artifact["checksum"]["digest"] = serde_json::json!(identity.digest());
        }
    }
    let effective = serde_json::to_vec_pretty(&effective).expect("effective bytes");
    let manifest = serde_json::to_vec_pretty(&manifest).expect("manifest bytes");
    let error = validate_published_resources(
        &destination,
        &effective_path,
        &effective,
        &manifest_path,
        &manifest,
    )
    .expect_err("correlated canonical mutation must fail before publication");
    assert!(
        error.context()["detail"].contains("canónico"),
        "must reject the semantic/canonical mismatch, not a stale checksum: {}",
        error.render_human()
    );
    fs::remove_dir_all(root).expect("cleanup");
}
