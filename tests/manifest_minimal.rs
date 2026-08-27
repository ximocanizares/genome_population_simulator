use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::results::{
    ArtifactRecord, GeneralManifest, ManifestReductionLink, ManifestReferenceLink, UnitState,
    general_manifest_document, validate_general_manifest_bytes,
};
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn complete_and_failed_minimal_manifests_pass_rust_and_external_validators() {
    assert_accepted("complete", &complete_manifest());
    assert_accepted("failed", &failed_manifest());
}

#[test]
fn truncation_unknown_field_and_absent_link_fail_both_validators() {
    let complete = complete_manifest();
    let mut bytes = serde_json::to_vec(&complete).expect("manifest JSON");
    bytes.truncate(bytes.len() / 2);
    assert!(validate_general_manifest_bytes(Path::new("manifest.json"), &bytes).is_err());
    assert_external_bytes("truncated", &bytes, false);

    let mut unknown = complete.clone();
    unknown
        .as_object_mut()
        .expect("root object")
        .insert("sets".to_owned(), json!([]));
    assert_rejected("unknown-root-field", &unknown);

    let mut absent = complete;
    absent["references"][0]["fasta_artifact_id"] = json!("artifact-absent");
    assert_rejected("absent-artifact-link", &absent);
}

#[test]
fn complete_manifest_requires_the_canonical_truth_artifact_and_minimum_category() {
    let complete = complete_manifest();
    let mut missing = complete.clone();
    missing["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .retain(|artifact| artifact["id"] != "truth-jsonl");
    missing["output_categories"]["minimum"] =
        json!(["config-general", "effective-config", "reference-fasta"]);
    assert_rejected("complete-manifest-missing-truth", &missing);

    let mut wrong_path = complete;
    wrong_path["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .iter_mut()
        .find(|artifact| artifact["id"] == "truth-jsonl")
        .expect("truth artifact")["path"] = json!("truth.jsonl");
    assert_rejected("complete-manifest-mutated-truth-path", &wrong_path);
}

#[test]
fn v3_complete_manifest_requires_complete_event_links_in_rust_and_external_schema() {
    let mut document = complete_manifest();
    document["schema_version"] = json!("3");
    document["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .push(artifact(
            "event-manifest-reads-main",
            "event_manifest",
            "events/reads-main/manifest.json",
        ));
    document["output_categories"]["minimum"]
        .as_array_mut()
        .expect("minimum")
        .push(json!("event-manifest-reads-main"));
    document["events"] = json!([{
        "event_id": "reads-main",
        "state": "failed",
        "manifest_artifact_id": "event-manifest-reads-main"
    }]);

    assert_rejected("v3-complete-manifest-failed-event", &document);
}

#[test]
fn v3_manifest_validation_reports_the_v3_schema_label() {
    let mut document = complete_manifest();
    document["schema_version"] = json!("3");
    document["state"] = json!("unexpected");
    let error = validate_general_manifest_bytes(
        Path::new("manifest-v3.json"),
        &serde_json::to_vec(&document).expect("manifest JSON"),
    )
    .expect_err("invalid v3 manifest");
    assert_eq!(
        error.context().get("schema").map(String::as_str),
        Some("manifest-general-v3.schema.json")
    );
}

#[test]
fn population_snvs_link_is_discriminated_by_mode_in_rust_and_schema() {
    let mut generated_without_id = complete_manifest();
    generated_without_id["population_snvs"] = json!({
        "mode": "generated", "metadata_sha256": "a".repeat(64), "normalization_audit": []
    });
    assert_rejected(
        "generated-population-snvs-requires-an-artifact",
        &generated_without_id,
    );

    let mut generated_with_catalog = complete_manifest();
    generated_with_catalog["artifacts"].as_array_mut().expect("catalog").extend([
        json!({"id":"population-snvs-tsv", "kind":"population_snvs_tsv", "path":"population_snvs.tsv", "size_bytes":1, "checksum":{"algorithm":"sha256", "digest":"b".repeat(64)}}),
        json!({"id":"population-snvs-metadata", "kind":"population_snvs_metadata", "path":"population_snvs.metadata.json", "size_bytes":1, "checksum":{"algorithm":"sha256", "digest":"a".repeat(64)}}),
    ]);
    generated_with_catalog["output_categories"]["minimum"] = json!([
        "config-general",
        "effective-config",
        "reference-fasta",
        "truth-jsonl",
        "population-snvs-tsv",
        "population-snvs-metadata"
    ]);
    generated_with_catalog["population_snvs"] = json!({
        "mode": "generated", "metadata_sha256": "a".repeat(64), "metadata_artifact_id": "population-snvs-metadata", "normalization_audit": []
    });
    assert_accepted(
        "generated-population-snvs-accepts-a-nonempty-metadata-artifact",
        &generated_with_catalog,
    );

    let mut malformed_audit = generated_with_catalog.clone();
    malformed_audit["population_snvs"]["normalization_audit"] = json!([{
        "population_id": "P1", "contig": "chr1", "position": "0", "ref": "N",
        "tolerance": 1e-9, "canonical_order": "population_id,contig,position,ref,alt,variant_id,source_id",
        "sum_algorithm": "neumaier_sum_v1", "original_sum": 0.0, "original_delta": 1.0,
        "effective_sum": 0.0, "effective_delta": 1.0, "reference_frequency": 1.0,
        "decision": "reference_remainder", "values": [{
            "alt": "N", "variant_id": "VAR_BOGUS", "source_id": "", "original": 0.0, "effective": 0.0
        }]
    }]);
    assert_rejected(
        "population-snvs-audit-rejects-invalid-position-and-alleles",
        &malformed_audit,
    );

    let mut reuse_with_id = generated_with_catalog.clone();
    reuse_with_id["population_snvs"]["mode"] = json!("reuse");
    assert_rejected(
        "reuse-population-snvs-must-not-name-an-artifact",
        &reuse_with_id,
    );
}

#[test]
fn initial_population_statistics_are_closed_in_rust_and_schema() {
    let mut valid = complete_manifest();
    valid["artifacts"].as_array_mut().expect("catalog").extend([
        json!({"id":"population-snvs-tsv", "kind":"population_snvs_tsv", "path":"population_snvs.tsv", "size_bytes":1, "checksum":{"algorithm":"sha256", "digest":"b".repeat(64)}}),
        json!({"id":"population-snvs-metadata", "kind":"population_snvs_metadata", "path":"population_snvs.metadata.json", "size_bytes":1, "checksum":{"algorithm":"sha256", "digest":"a".repeat(64)}}),
    ]);
    valid["output_categories"]["minimum"] = json!([
        "config-general",
        "effective-config",
        "reference-fasta",
        "truth-jsonl",
        "population-snvs-tsv",
        "population-snvs-metadata"
    ]);
    valid["population_snvs"] = json!({
        "mode": "generated", "metadata_sha256": "a".repeat(64), "metadata_artifact_id": "population-snvs-metadata", "normalization_audit": [],
        "initial_population_statistics": initial_population_statistics_fixture(),
    });
    assert_accepted("initial-population-statistics", &valid);

    let mut too_few_algorithms = valid.clone();
    too_few_algorithms["population_snvs"]["initial_population_statistics"]["provenance"]["genotype_algorithms"] =
        json!(["hwe-composition-v1", "f-diploid-genotype-v1"]);
    assert_rejected("statistics-require-three-algorithms", &too_few_algorithms);

    let mut duplicate_population = valid.clone();
    let population =
        duplicate_population["population_snvs"]["initial_population_statistics"]["populations"][0]
            .clone();
    duplicate_population["population_snvs"]["initial_population_statistics"]["populations"]
        .as_array_mut()
        .expect("populations")
        .push(population);
    assert_rejected("statistics-duplicate-population", &duplicate_population);

    let mut frequency_mismatch = valid.clone();
    frequency_mismatch["population_snvs"]["initial_population_statistics"]["populations"][0]["loci"]
        [0]["alleles"][0]["realized_frequency"] = json!(0.25);
    let bytes = serde_json::to_vec(&frequency_mismatch).expect("manifest JSON");
    assert!(
        validate_general_manifest_bytes(Path::new("manifest.json"), &bytes).is_err(),
        "Rust accepted an inconsistent frequency/count pair"
    );
    assert_external_bytes("statistics-frequency-count-mismatch", &bytes, true);

    let mut invalid_f_model = valid;
    invalid_f_model["population_snvs"]["initial_population_statistics"]["populations"][0]["genotype_model"] =
        json!({"kind": "f_diploid", "inbreeding_coefficient": 1.1});
    assert_rejected("statistics-f-model-out-of-range", &invalid_f_model);
}

#[test]
fn duplicate_ids_paths_and_output_categories_are_rejected() {
    let mut duplicate_id = complete_manifest();
    duplicate_id["artifacts"][1]["id"] = duplicate_id["artifacts"][0]["id"].clone();
    assert_rejected("duplicate-artifact-id", &duplicate_id);

    let mut duplicate_path = complete_manifest();
    duplicate_path["artifacts"][1]["path"] = duplicate_path["artifacts"][0]["path"].clone();
    assert_rejected("duplicate-artifact-path", &duplicate_path);

    let mut repeated_category = complete_manifest();
    repeated_category["output_categories"]["requested"] = json!(["reference-fasta"]);
    assert_rejected("artifact-in-two-categories", &repeated_category);

    let mut unclassified = complete_manifest();
    unclassified["output_categories"]["minimum"] = json!(["config-general"]);
    assert_rejected("unclassified-artifact", &unclassified);
}

#[test]
fn failed_diagnostics_must_exist_in_catalog_v2_with_canonical_meaning() {
    let mut unknown = failed_manifest();
    unknown["diagnostics"][0]["code"] = json!("GPS-FMT-999");
    assert_rejected("unknown-diagnostic-code", &unknown);

    let mut altered = failed_manifest();
    altered["diagnostics"][0]["message"] = json!("mensaje autofirmado");
    assert_rejected("altered-diagnostic-meaning", &altered);

    let mut full_catalog_positive = failed_manifest();
    full_catalog_positive["diagnostics"][0] = json!({
        "code": "GPS-DOMAIN-003",
        "message": "Una invariante del estado biológico o técnico válido no se cumple."
    });
    assert_accepted("full-catalog-domain-003", &full_catalog_positive);

    assert_accepted("cataloged-failed-diagnostic", &failed_manifest());
}

#[test]
fn output_categories_cannot_be_swapped_while_preserving_a_partition() {
    let mut document = complete_manifest();
    document["output_categories"]["minimum"] = json!(["config-general", "effective-config"]);
    document["output_categories"]["requested"] = json!(["reference-fasta"]);
    assert_rejected("source-reference-swapped-to-requested", &document);
}

#[test]
fn source_and_derived_fasta_cannot_alias_across_output_categories() {
    let artifacts = vec![
        record("shared", "reference_fasta", "resources/shared.fa"),
        record(
            "map",
            "reference_reduction_map",
            "resources/derived.map.json",
        ),
    ];
    let references = vec![
        ManifestReferenceLink {
            id: "main".to_owned(),
            fasta_artifact_id: "shared".to_owned(),
        },
        ManifestReferenceLink {
            id: "small".to_owned(),
            fasta_artifact_id: "shared".to_owned(),
        },
    ];
    let reductions = vec![ManifestReductionLink {
        id: "small".to_owned(),
        source_reference_id: "main".to_owned(),
        eligible_bases: 20,
        region_size_bases: 10,
        target_mode: "target_bases".to_owned(),
        target_declared: "10".to_owned(),
        target_bases: 10,
        root_seed: "0".repeat(64),
        selection_subseed: "1".repeat(64),
        selection_order: vec![1],
        state: UnitState::Complete,
        fasta_artifact_id: Some("shared".to_owned()),
        map_artifact_id: Some("map".to_owned()),
        source_fasta_artifact_id: Some("shared".to_owned()),
        map_document_kind: "reference_reduction_map".to_owned(),
        map_schema_version: "1".to_owned(),
    }];
    let diagnostic = validate_general_manifest_bytes(Path::new("truncated.json"), b"{")
        .expect_err("controlled format diagnostic");
    let manifest = GeneralManifest {
        execution_id: "execution-alias",
        state: UnitState::Failed,
        configuration_input: None,
        effective_config: None,
        references: &references,
        reduced_references: &reductions,
        population_snvs: None,
        vcf_outputs: &[],
        events: None,
        authorizations: genome_population_simulator::execution::ExecutionAuthorizations::default(),
        artifacts: &artifacts,
        diagnostics: std::slice::from_ref(&diagnostic),
    };
    assert_rejected(
        "source-derived-fasta-alias",
        &general_manifest_document(&manifest),
    );
}

#[test]
fn artifact_paths_are_canonical_and_digests_are_lowercase_sha256() {
    for value in [
        json!("../outside"),
        json!("a//b"),
        json!("a\\b"),
        json!("trailing/"),
        json!("INCOMPLETE"),
    ] {
        let mut document = complete_manifest();
        document["artifacts"][0]["path"] = value.clone();
        assert_rejected(&format!("invalid-path-{value}"), &document);
    }
    let mut digest = complete_manifest();
    digest["artifacts"][0]["checksum"]["digest"] = json!("A".repeat(64));
    assert_rejected("uppercase-digest", &digest);
}

#[test]
fn failed_reduction_keeps_reached_links_and_requested_classification() {
    let artifacts = vec![
        record("source", "reference_fasta", "resources/source.fa"),
        record("derived", "reference_fasta", "resources/derived.fa"),
        record(
            "map",
            "reference_reduction_map",
            "resources/derived.map.json",
        ),
    ];
    let references = vec![
        ManifestReferenceLink {
            id: "main".to_owned(),
            fasta_artifact_id: "source".to_owned(),
        },
        ManifestReferenceLink {
            id: "small".to_owned(),
            fasta_artifact_id: "derived".to_owned(),
        },
    ];
    let reductions = vec![ManifestReductionLink {
        id: "small".to_owned(),
        source_reference_id: "main".to_owned(),
        eligible_bases: 20,
        region_size_bases: 10,
        target_mode: "target_bases".to_owned(),
        target_declared: "10".to_owned(),
        target_bases: 10,
        root_seed: "0".repeat(64),
        selection_subseed: "1".repeat(64),
        selection_order: vec![1],
        state: UnitState::Failed,
        fasta_artifact_id: Some("derived".to_owned()),
        map_artifact_id: Some("map".to_owned()),
        source_fasta_artifact_id: Some("source".to_owned()),
        map_document_kind: "reference_reduction_map".to_owned(),
        map_schema_version: "1".to_owned(),
    }];
    let diagnostic = validate_general_manifest_bytes(Path::new("truncated.json"), b"{")
        .expect_err("controlled format diagnostic");
    let manifest = GeneralManifest {
        execution_id: "execution-partial",
        state: UnitState::Failed,
        configuration_input: None,
        effective_config: None,
        references: &references,
        reduced_references: &reductions,
        population_snvs: None,
        vcf_outputs: &[],
        events: None,
        authorizations: genome_population_simulator::execution::ExecutionAuthorizations::default(),
        artifacts: &artifacts,
        diagnostics: std::slice::from_ref(&diagnostic),
    };
    let document = general_manifest_document(&manifest);
    assert_eq!(
        document["output_categories"]["requested"],
        json!(["derived", "map"])
    );
    assert_eq!(document["output_categories"]["minimum"], json!(["source"]));
    assert_accepted("failed-partial", &document);
}

fn record(id: &str, kind: &str, path: &str) -> ArtifactRecord {
    ArtifactRecord {
        id: id.to_owned(),
        kind: kind.to_owned(),
        path: path.to_owned(),
        size_bytes: 1,
        digest: "0".repeat(64),
    }
}

fn complete_manifest() -> Value {
    json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": "2",
        "execution_id": "execution-001",
        "state": "complete",
        "configuration_input": {
            "kind": "general",
            "original_path": "config/general.toml",
            "artifact_id": "config-general"
        },
        "effective_config_artifact_id": "effective-config",
        "references": [{"id": "ref", "fasta_artifact_id": "reference-fasta"}],
        "reference_reductions": [],
        "vcf_outputs": [],
        "authorizations": {"allow_empty_snv_matrix": false, "allow_incompatible_genotypes": false},
        "output_categories": {
            "minimum": ["config-general", "effective-config", "reference-fasta", "truth-jsonl"],
            "defaulted": [],
            "requested": [],
            "disabled": []
        },
        "artifacts": [
            artifact("config-general", "config_general_toml", "inputs/config/general.toml"),
            artifact("reference-fasta", "reference_fasta", "resources/references/ref.fa"),
            artifact("effective-config", "effective_config", "effective-config.json"),
            artifact("truth-jsonl", "truth_jsonl", "outputs/truth.jsonl")
        ],
        "diagnostics": []
    })
}

fn failed_manifest() -> Value {
    json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": "2",
        "execution_id": "execution-failed-001",
        "state": "failed",
        "configuration_input": {
            "kind": "general",
            "original_path": "config/general.toml",
            "artifact_id": "config-general"
        },
        "effective_config_artifact_id": null,
        "references": [],
        "reference_reductions": [],
        "vcf_outputs": [],
        "authorizations": {"allow_empty_snv_matrix": false, "allow_incompatible_genotypes": false},
        "output_categories": {
            "minimum": ["config-general"],
            "defaulted": [],
            "requested": [],
            "disabled": []
        },
        "artifacts": [
            artifact("config-general", "config_general_toml", "inputs/config/general.toml")
        ],
        "diagnostics": [{
            "code": "GPS-IO-001",
            "message": "Una operación de lectura o escritura falló."
        }]
    })
}

fn artifact(id: &str, kind: &str, path: &str) -> Value {
    json!({
        "id": id,
        "kind": kind,
        "path": path,
        "size_bytes": 1,
        "checksum": {"algorithm": "sha256", "digest": "0".repeat(64)}
    })
}

fn initial_population_statistics_fixture() -> Value {
    json!({
        "kind": "initial_population_statistics",
        "version": "1",
        "provenance": {
            "matrix_reference_id": "ref",
            "matrix_reference_fingerprint": "b".repeat(64),
            "matrix_canonical_tsv_sha256": "c".repeat(64),
            "genotype_algorithms": ["hwe-composition-v1", "f-diploid-genotype-v1", "phase-permutation-v1"],
            "missing_locus_policy": "missing_by_deletion_excluded_from_allele_frequency_denominator_v1",
            "explicit_snp_constraints": []
        },
        "populations": [{
            "population_id": "P1", "genotype_model": {"kind": "hwe"}, "individual_count": 1, "ploidy": 2,
            "loci": [{
                "contig": "chr1", "position_1based": 1, "ref": "A", "total_copies": 2, "observed_copies": 2, "missing_by_deletion": 0,
                "alleles": [
                    {"allele": "A", "kind": "ref", "requested_frequency": 0.5, "realized_count": 1, "realized_frequency": 0.5},
                    {"allele": "C", "kind": "alt", "variant_id": "v1", "source_id": "test", "requested_frequency": 0.5, "realized_count": 1, "realized_frequency": 0.5}
                ],
                "complete_genotype_observations": 1,
                "genotypes": [
                    {"allele_dosages": [2, 0], "expected_probability": 0.25, "expected_count": 0.25, "realized_count": 0},
                    {"allele_dosages": [1, 1], "expected_probability": 0.5, "expected_count": 0.5, "realized_count": 1},
                    {"allele_dosages": [0, 2], "expected_probability": 0.25, "expected_count": 0.25, "realized_count": 0}
                ]
            }]
        }]
    })
}

fn assert_accepted(name: &str, document: &Value) {
    let bytes = serde_json::to_vec(document).expect("manifest JSON");
    validate_general_manifest_bytes(Path::new("manifest.json"), &bytes)
        .unwrap_or_else(|error| panic!("Rust rejected positive {name}: {error:?}"));
    assert_external_bytes(name, &bytes, true);
}

fn assert_rejected(name: &str, document: &Value) {
    let bytes = serde_json::to_vec(document).expect("manifest JSON");
    assert!(
        validate_general_manifest_bytes(Path::new("manifest.json"), &bytes).is_err(),
        "Rust accepted negative {name}"
    );
    assert_external_bytes(name, &bytes, false);
}

fn assert_external_bytes(name: &str, bytes: &[u8], expected_success: bool) {
    let instance = temporary(name);
    fs::write(&instance, bytes).expect("write external instance");
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let schema = match serde_json::from_slice::<Value>(bytes)
        .ok()
        .and_then(|document| {
            document
                .get("schema_version")
                .and_then(Value::as_str)
                .map(str::to_owned)
        })
        .as_deref()
    {
        Some("3") => "manifest-general-v3.schema.json",
        _ => "manifest-general-v2.schema.json",
    };
    let output = Command::new("python3")
        .arg(root.join("tests/external_manifest_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(&instance)
        .arg(root.join("catalogs/diagnostics-v2.tsv"))
        .output()
        .expect("independent Python validator");
    let _ = fs::remove_file(&instance);
    assert_eq!(
        output.status.success(),
        expected_success,
        "external validator mismatch for {name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn temporary(_label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f14-manifest-{}-{}.json",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
