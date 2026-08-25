use std::fs;
use std::path::Path;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn real_effective_shapes_keep_their_two_valid_source_modes() {
    assert_valid(
        "simulated-without-reductions",
        &effective_config("simulated", false),
    );
    assert_valid(
        "supplied-with-reduction",
        &effective_config("supplied", true),
    );
}

#[test]
fn simulated_reference_cannot_claim_a_reduced_reference() {
    assert_invalid(
        "simulated-with-reduction",
        &effective_config("simulated", true),
    );
}

#[test]
fn population_snvs_effective_link_is_discriminated_by_mode() {
    let mut reuse = effective_config("simulated", false);
    reuse["population_snvs"] =
        json!({"mode":"reuse", "metadata_sha256":"a".repeat(64), "normalization_audit": []});
    assert_valid("reuse-population-snvs-may-omit-artifact-id", &reuse);

    let mut reuse_with_id = reuse.clone();
    reuse_with_id["population_snvs"]["metadata_artifact_id"] = json!("metadata");
    assert_invalid("reuse-population-snvs-rejects-artifact-id", &reuse_with_id);

    let mut generated = effective_config("simulated", false);
    generated["population_snvs"] =
        json!({"mode":"generated", "metadata_sha256":"a".repeat(64), "normalization_audit": []});
    assert_invalid("generated-population-snvs-requires-artifact-id", &generated);

    let mut malformed_audit = reuse.clone();
    malformed_audit["population_snvs"]["normalization_audit"] = json!([{
        "population_id": "P1", "contig": "chr1", "position": "0", "ref": "N",
        "tolerance": 1e-9, "canonical_order": "population_id,contig,position,ref,alt,variant_id,source_id",
        "sum_algorithm": "neumaier_sum_v1", "original_sum": 0.0, "original_delta": 1.0,
        "effective_sum": 0.0, "effective_delta": 1.0, "reference_frequency": 1.0,
        "decision": "reference_remainder", "values": [{
            "alt": "N", "variant_id": "VAR_BOGUS", "source_id": "", "original": 0.0, "effective": 0.0
        }]
    }]);
    assert_invalid(
        "population-snvs-audit-rejects-invalid-position-and-alleles",
        &malformed_audit,
    );
}

#[test]
fn initial_population_statistics_schema_is_closed_and_typed() {
    let mut valid = effective_config("simulated", false);
    valid["initial_populations"] = json!([{
        "population_id": "P1", "individual_count": 1, "individual_ids": ["I1"], "ploidy": 2,
        "genotype_model": {"kind": "hwe"}
    }]);
    valid["population_snvs"] = json!({
        "mode": "reuse",
        "metadata_sha256": "a".repeat(64),
        "normalization_audit": [],
        "initial_population_statistics": initial_population_statistics_fixture(),
    });
    assert_valid("initial-population-statistics", &valid);

    let mut missing_statistics = valid.clone();
    missing_statistics["population_snvs"]
        .as_object_mut()
        .expect("population SNVs")
        .remove("initial_population_statistics");
    assert_invalid(
        "initial-populations-require-initial-statistics",
        &missing_statistics,
    );

    let cases = [
        (
            "statistics-unknown-provenance",
            "/population_snvs/initial_population_statistics/provenance/future",
            json!(true),
        ),
        (
            "statistics-invalid-fingerprint",
            "/population_snvs/initial_population_statistics/provenance/matrix_reference_fingerprint",
            json!("A".repeat(64)),
        ),
        (
            "statistics-invalid-missing-policy",
            "/population_snvs/initial_population_statistics/provenance/missing_locus_policy",
            json!("future"),
        ),
        (
            "statistics-alt-without-source",
            "/population_snvs/initial_population_statistics/populations/0/loci/0/alleles/1/source_id",
            Value::Null,
        ),
        (
            "statistics-invalid-dosage",
            "/population_snvs/initial_population_statistics/populations/0/loci/0/genotypes/0/allele_dosages/0",
            json!(-1),
        ),
        (
            "statistics-too-few-algorithms",
            "/population_snvs/initial_population_statistics/provenance/genotype_algorithms",
            json!(["hwe-composition-v1", "f-diploid-genotype-v1"]),
        ),
    ];
    for (name, pointer, replacement) in cases {
        let mut invalid = valid.clone();
        if name == "statistics-unknown-provenance" {
            invalid["population_snvs"]["initial_population_statistics"]["provenance"]["future"] =
                replacement;
        } else if replacement.is_null() {
            let (parent, field) = pointer.rsplit_once('/').expect("field pointer");
            invalid
                .pointer_mut(parent)
                .and_then(Value::as_object_mut)
                .expect("parent object")
                .remove(field);
        } else {
            *invalid.pointer_mut(pointer).expect("fixture pointer") = replacement;
        }
        assert_invalid(name, &invalid);
    }
}

#[test]
fn effective_initial_population_schema_enforces_local_bounds_and_f_diploidy() {
    let schema: Value =
        serde_json::from_str(include_str!("../schemas/effective-config-v2.schema.json"))
            .expect("effective schema JSON");
    let cohort_schema = &schema["$defs"]["initial_population"];
    assert_eq!(
        cohort_schema["properties"]["individual_count"]["maximum"],
        json!(u64::MAX)
    );
    assert_eq!(cohort_schema["properties"]["individual_ids"]["minItems"], 1);

    let mut valid = effective_config("simulated", false);
    valid["initial_populations"] = json!([{
        "population_id": "P1",
        "individual_count": 2,
        "individual_ids": ["IND_1", "IND_2"],
        "ploidy": 2,
        "genotype_model": {"kind": "f_diploid", "inbreeding_coefficient": 0.25}
    }]);
    add_initial_population_statistics(&mut valid);
    assert_valid("effective-initial-population-valid", &valid);

    let mut derived = valid.clone();
    derived["initial_populations"][0]["founder_keys"] = json!(["founder-A", "founder-B"]);
    assert_valid("effective-initial-population-derived-ids", &derived);

    let mut same_key_in_another_population = derived.clone();
    same_key_in_another_population["initial_populations"] = json!([
        derived["initial_populations"][0].clone(),
        {
            "population_id": "P2",
            "individual_count": 2,
            "individual_ids": ["P2_IND_1", "P2_IND_2"],
            "founder_keys": ["founder-A", "founder-B"],
            "ploidy": 2,
            "genotype_model": {"kind": "hwe"}
        }
    ]);
    assert_valid(
        "effective-initial-population-founder-keys-are-local-to-population",
        &same_key_in_another_population,
    );

    let mut duplicate_founder_keys = derived;
    duplicate_founder_keys["initial_populations"][0]["founder_keys"] =
        json!(["founder-A", "founder-A"]);
    assert_invalid(
        "effective-initial-population-duplicate-founder-keys",
        &duplicate_founder_keys,
    );

    let mut empty_ids = valid.clone();
    empty_ids["initial_populations"][0]["individual_ids"] = json!([]);
    assert_invalid("effective-initial-population-empty-ids", &empty_ids);

    let mut duplicated_ids = valid.clone();
    duplicated_ids["initial_populations"][0]["individual_ids"] = json!(["IND_1", "IND_1"]);
    assert_invalid(
        "effective-initial-population-duplicate-ids",
        &duplicated_ids,
    );

    let mut f_outside_diploidy = valid;
    f_outside_diploidy["initial_populations"][0]["ploidy"] = json!(3);
    assert_invalid(
        "effective-initial-population-f-outside-diploidy",
        &f_outside_diploidy,
    );
}

#[test]
fn effective_selection_order_is_strictly_one_based() {
    let mut document = effective_config("supplied", true);
    *document
        .pointer_mut("/reduced_references/0/selection_order/0")
        .expect("selection order") = json!(0);
    assert_invalid("zero-selection-order", &document);
}

#[test]
fn effective_reduction_selection_algorithm_is_literal_and_required() {
    let mut wrong = effective_config("supplied", true);
    wrong["reduced_references"][0]["selection_algorithm"] = json!("future_algorithm");
    assert_invalid("wrong-reduction-selection-algorithm", &wrong);

    let mut absent = effective_config("supplied", true);
    absent["reduced_references"][0]
        .as_object_mut()
        .expect("reduction object")
        .remove("selection_algorithm");
    assert_invalid("absent-reduction-selection-algorithm", &absent);
}

#[test]
fn derived_sequence_catalog_is_required_closed_and_reversible() {
    let document = effective_config("supplied", true);
    assert_valid("derived-sequence-catalog", &document);

    for (name, pointer, replacement) in [
        (
            "derived-id",
            "/reduced_references/0/sequences/0/id",
            json!(""),
        ),
        (
            "derived-class",
            "/reduced_references/0/sequences/0/class",
            json!("contig"),
        ),
        (
            "derived-length",
            "/reduced_references/0/sequences/0/length",
            json!(0),
        ),
        (
            "derived-header",
            "/reduced_references/0/sequences/0/header/value",
            json!(""),
        ),
        (
            "derived-description",
            "/reduced_references/0/sequences/0/description/encoding",
            json!("latin-1"),
        ),
        (
            "derived-alphabet",
            "/reduced_references/0/sequences/0/alphabet_validation",
            json!("ACGT"),
        ),
        (
            "derived-n-count",
            "/reduced_references/0/sequences/0/n_count",
            json!(-1),
        ),
    ] {
        let mut changed = document.clone();
        *changed.pointer_mut(pointer).expect(name) = replacement;
        assert_invalid(name, &changed);
    }

    let mut unknown = document.clone();
    unknown["reduced_references"][0]["sequences"][0]["future"] = json!(true);
    assert_invalid("derived-sequence-unknown-field", &unknown);
}

#[test]
fn every_effective_artifact_link_is_a_non_empty_catalog_id() {
    for id in ["reference-main", "a", "referencia-ñ"] {
        let mut document = effective_config("supplied", true);
        *document
            .pointer_mut("/reference/fasta_artifact_id")
            .expect("reference artifact ID") = json!(id);
        assert_valid(&format!("artifact-id-{id:?}"), &document);
    }

    for value in [json!(""), json!(null), json!(1)] {
        let mut document = effective_config("supplied", true);
        *document
            .pointer_mut("/reference/fasta_artifact_id")
            .expect("reference artifact ID") = value.clone();
        assert_invalid(&format!("invalid-artifact-id-{value:?}"), &document);
    }
}

#[test]
fn semantic_ids_source_path_and_selection_order_cannot_be_empty() {
    for (name, pointer) in [
        ("source-path", "/source_path"),
        ("reference-id", "/reference/id"),
        ("sequence-id", "/reference/sequences/0/id"),
        ("reduction-id", "/reduced_references/0/id"),
        (
            "source-reference-id",
            "/reduced_references/0/source_reference_id",
        ),
        ("target-declared", "/reduced_references/0/target/declared"),
    ] {
        let mut document = effective_config("supplied", true);
        *document.pointer_mut(pointer).expect("fixture pointer") = json!("");
        assert_invalid(name, &document);
    }

    let mut document = effective_config("supplied", true);
    *document
        .pointer_mut("/reduced_references/0/selection_order")
        .expect("selection order") = json!([]);
    assert_invalid("empty-selection-order", &document);
}

#[test]
fn effective_ploidy_is_required_and_limited_to_the_f2_domain() {
    let mut missing = effective_config("simulated", false);
    missing
        .as_object_mut()
        .expect("effective config object")
        .remove("ploidy");
    assert_invalid("missing-ploidy", &missing);

    for value in [json!(0), json!(11), json!(1.5), json!("2")] {
        let mut document = effective_config("simulated", false);
        document["ploidy"] = value.clone();
        assert_invalid(&format!("invalid-ploidy-{value}"), &document);
    }
}

#[test]
fn initial_populations_publish_only_resolved_ids_and_validated_hwe_or_f_models() {
    let mut document = effective_config("simulated", false);
    document["initial_populations"] = json!([
        {
            "population_id": "P1",
            "individual_count": 2,
            "individual_ids": ["IND_A", "IND_B"],
            "ploidy": 1,
            "genotype_model": {"kind": "hwe"}
        },
        {
            "population_id": "P2",
            "individual_count": 2,
            "individual_ids": ["IND_C", "IND_D"],
            "ploidy": 2,
            "genotype_model": {"kind": "f_diploid", "inbreeding_coefficient": 0.5}
        }
    ]);
    add_initial_population_statistics(&mut document);
    assert_valid("initial-populations", &document);

    let mut wrong_f_ploidy = document.clone();
    wrong_f_ploidy["initial_populations"][1]["ploidy"] = json!(3);
    assert_invalid("initial-populations-f-outside-diploidy", &wrong_f_ploidy);
}

#[test]
fn reference_origin_and_reversible_metadata_are_closed_and_source_specific() {
    let mut reversible_bytes = effective_config("supplied", false);
    reversible_bytes["reference"]["sequences"][0]["header"] =
        json!({"encoding": "bytes-hex", "value": "63687231ff"});
    reversible_bytes["reference"]["sequences"][0]["description"] =
        json!({"encoding": "bytes-hex", "value": "ff"});
    assert_valid("reversible-non-utf8-header", &reversible_bytes);

    let mut simulated_path = effective_config("simulated", false);
    simulated_path["reference"]["original_fasta_path"] = json!("source.fa");
    assert_invalid("simulated-original-path", &simulated_path);

    let mut supplied_format = effective_config("supplied", false);
    supplied_format["reference"]["canonical_format"] = json!({
        "version": "1",
        "line_base_count": 80,
        "line_ending": "LF",
        "base_case": "uppercase"
    });
    assert_invalid("supplied-canonical-format", &supplied_format);

    let mut wrong_width = effective_config("simulated", false);
    wrong_width["reference"]["canonical_format"]["line_base_count"] = json!(79);
    assert_invalid("wrong-canonical-width", &wrong_width);

    let mut invalid_hex = effective_config("supplied", false);
    invalid_hex["reference"]["sequences"][0]["header"] =
        json!({"encoding": "bytes-hex", "value": "0G"});
    assert_invalid("invalid-header-hex", &invalid_hex);

    let mut wrong_alphabet = effective_config("supplied", false);
    wrong_alphabet["reference"]["sequences"][0]["alphabet_validation"] = json!("ACGT");
    assert_invalid("wrong-alphabet-contract", &wrong_alphabet);
}

fn effective_config(source: &str, with_reduction: bool) -> Value {
    let composition = match source {
        "simulated" => json!({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}),
        "supplied" => Value::Null,
        _ => panic!("unsupported fixture source"),
    };
    let reductions = if with_reduction {
        vec![json!({
            "id": "reduced",
            "source_reference_id": "main",
            "region_size_bases": 100,
            "target": {"mode": "target_bases", "declared": "200", "bases": 200},
            "selection_algorithm": "reduced_reference_selection_v1",
            "selection_subseed": "1".repeat(64),
            "selection_order": [1, 2],
            "sequences": [{
                "id": "chr1",
                "class": "chromosome",
                "length": 200,
                "header": {"encoding": "utf-8", "value": "chr1"},
                "description": {"encoding": "bytes-hex", "value": "ff"},
                "alphabet_validation": "ACGTN-case-insensitive",
                "n_count": 1
            }],
            "state": "complete",
            "fasta_artifact_id": "reference-reduced",
            "map_artifact_id": "reference-map",
            "source_fasta_artifact_id": "reference-main"
        })]
    } else {
        Vec::new()
    };
    json!({
        "effective_config_kind": "gps_effective_config",
        "schema_version": "2",
        "source_path": "config/general.toml",
        "ploidy": 2,
        "normalization": {"tolerance": 1e-9, "canonical_order": "population_id,contig,position,ref,alt,variant_id,source_id", "sum_algorithm": "neumaier_sum_v1"},
        "seed": {
            "root": "0".repeat(64),
            "origin": "declared"
        },
        "authorizations": {"allow_empty_snv_matrix": false, "allow_incompatible_genotypes": false},
        "reference": {
            "id": "main",
            "source": source,
            "original_fasta_path": match source {
                "simulated" => Value::Null,
                "supplied" => json!("inputs/reference.fa"),
                _ => unreachable!(),
            },
            "canonical_format": match source {
                "simulated" => json!({
                    "version": "1",
                    "line_base_count": 80,
                    "line_ending": "LF",
                    "base_case": "uppercase"
                }),
                "supplied" => Value::Null,
                _ => unreachable!(),
            },
            "total_genome_size": 1000,
            "composition": composition,
            "sequences": [{
                "id": "chr1",
                "class": "chromosome",
                "length": 1000,
                "header": {"encoding": "utf-8", "value": "chr1"},
                "description": null,
                "alphabet_validation": "ACGTN-case-insensitive",
                "unknown_base_count": 0
            }],
            "fasta_artifact_id": "reference-main"
        },
        "reduced_references": reductions,
        "initial_populations": []
    })
}

fn initial_population_statistics_fixture() -> Value {
    json!({
        "kind": "initial_population_statistics",
        "version": "1",
        "provenance": {
            "matrix_reference_id": "main",
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

fn add_initial_population_statistics(document: &mut Value) {
    document["population_snvs"] = json!({
        "mode": "reuse",
        "metadata_sha256": "a".repeat(64),
        "normalization_audit": [],
        "initial_population_statistics": initial_population_statistics_fixture(),
    });
}

fn assert_valid(name: &str, document: &Value) {
    assert_external(name, document, true);
}

fn assert_invalid(name: &str, document: &Value) {
    assert_external(name, document, false);
}

fn assert_external(name: &str, document: &Value, expected_success: bool) {
    let instance = std::env::temp_dir().join(format!(
        "gps-effective-schema-{}-{}-{name}.json",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::write(&instance, serde_json::to_vec(document).expect("JSON")).expect("instance");
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas/effective-config-v2.schema.json"))
        .arg(&instance)
        .output()
        .expect("independent Python validator");
    let _ = fs::remove_file(instance);
    assert_eq!(
        output.status.success(),
        expected_success,
        "effective-config schema mismatch for {name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}
