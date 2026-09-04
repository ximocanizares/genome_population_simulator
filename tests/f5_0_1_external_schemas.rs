#![recursion_limit = "256"]

use std::fs;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

use serde_json::{Value, json};

static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

fn assert_external_cases(schema: &str, dependencies: &[&str], cases: Vec<(Value, bool)>) {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let serial = NEXT_TEMP.fetch_add(1, Ordering::Relaxed);
    let bundle =
        std::env::temp_dir().join(format!("gps-f5-schema-{}-{serial}", std::process::id()));
    fs::create_dir(&bundle).expect("isolated schema bundle");

    for name in std::iter::once(schema).chain(dependencies.iter().copied()) {
        assert!(
            !name.contains('/'),
            "schema bundle accepts only sibling names"
        );
        fs::copy(root.join("schemas").join(name), bundle.join(name)).expect("copy schema");
    }

    let case_path = bundle.join("cases.json");
    let case_document = cases
        .into_iter()
        .map(|(instance, expected)| json!({"instance": instance, "expected": expected}))
        .collect::<Vec<_>>();
    fs::write(
        &case_path,
        serde_json::to_vec(&case_document).expect("JSON cases"),
    )
    .expect("write cases");

    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(bundle.join(schema))
        .arg(&case_path)
        .arg("--cases")
        .stderr(Stdio::piped())
        .output()
        .expect("external validator");
    fs::remove_dir_all(&bundle).expect("remove isolated schema bundle");
    assert!(
        output.status.success(),
        "{schema}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn limits() -> Value {
    json!({
        "max_layout_attempts": 1,
        "max_haplotype_attempts": 1,
        "max_swap_proposals": 1,
        "max_haplotype_cells": 1,
        "max_operation_units": 1,
        "max_working_memory_bytes": 1,
        "max_output_bytes": 1
    })
}

fn target() -> Value {
    json!({
        "kind": "block_count",
        "declared_value": 1,
        "numerator": null,
        "denominator": null,
        "target_bases": null
    })
}

fn algorithms() -> Value {
    json!({
        "layout": "f5-layout-v1",
        "marginals": "f5-spectrum-v1",
        "haplotypes": "weighted-prefix-swap-v1",
        "ld": "phased-r2-neumaier-v1",
        "controlled": "largest-remainder-v1",
        "inbreeding": "diploid-f-v1",
        "truth_projection": "f5-founder-truth-projection-v2"
    })
}

fn normalization(set_kind: &str, canonical_order: &str) -> Value {
    json!({
        "set_kind": set_kind,
        "canonical_order": canonical_order,
        "algorithm": "neumaier-f64-v1",
        "tolerance": 0.0,
        "input_values": [{"key": "row-0", "value": 1.0}],
        "input_sum": 1.0,
        "input_delta": 0.0,
        "effective_values": [{"key": "row-0", "value": 1.0}],
        "output_sum": 1.0,
        "output_delta": 0.0
    })
}

fn normalizations() -> Value {
    json!([
        normalization("block_lengths", "length_min_bp"),
        normalization("snp_gaps", "gap_min_bp"),
        normalization("snp_spectrum", "alt_af_min"),
        normalization("haplotype_weights", "haplotype_rank")
    ])
}

fn preflight() -> Value {
    json!({
        "algorithm": "f5-resource-accounting-v6",
        "reference_bases_processed": 1,
        "fasta_bytes_processed": 1,
        "mask_bytes_processed": 0,
        "general_toml_bytes": 1,
        "scenario_toml_bytes": 1,
        "profile_metadata_bytes": 1,
        "profile_tsv_bytes": 5,
        "input_parser_peak_bytes": 16,
        "tract_count": 1,
        "retained_sequence_ids": 1,
        "site_bound_per_block": 1,
        "bin_count_total": 1,
        "profile_rows": 1,
        "normalization_rows": 4,
        "tsv_rows_bound": 1,
        "json_rows_bound": 1,
        "other_operation_units": 0,
        "other_working_memory_bytes": 0,
        "other_output_bytes": 0,
        "operation_units": 1,
        "working_memory_bytes": 1,
        "output_bytes": 1,
        "operation_limit": 1,
        "working_memory_limit_bytes": 1,
        "output_limit_bytes": 1
    })
}

fn profile_link() -> Value {
    json!({
        "profile_id": "profile",
        "path": "profiles/profile.metadata.json",
        "size_bytes": 1,
        "sha256": "a".repeat(64)
    })
}

fn catalog_link(catalog_id: &str, path: &str, digest: char) -> Value {
    json!({
        "catalog_id": catalog_id,
        "path": path,
        "size_bytes": 1,
        "sha256": digest.to_string().repeat(64)
    })
}

fn effective() -> Value {
    json!({
        "effective_config_kind": "gps_effective_config",
        "schema_version": "5",
        "source_path": "authoring/general.toml",
        "seed": "b".repeat(64),
        "f5": {
            "schema_version": "1",
            "scenario_id": "scenario",
            "reference_id": "reference",
            "reference_sha256": "c".repeat(64),
            "seed": "b".repeat(64),
            "ploidy": 2,
            "founder_count": 1,
            "subpopulations": [{
                "population_id": "P1",
                "individual_count": 1,
                "inbreeding_f": 0.0,
                "inbreeding_f_source": "population",
                "sampling": "probabilistic",
                "founders": [{"individual_id": "I1", "founder_key": null}]
            }],
            "profile": profile_link(),
            "mask": null,
            "target": target(),
            "sampling": [],
            "normalizations": normalizations(),
            "limits": limits(),
            "resource_preflight": preflight(),
            "algorithms": algorithms(),
            "rng_catalog": catalog_link("f5-rng-catalog-v1", "catalogs/f5-rng-catalog-v1.tsv", 'd'),
            "id_catalog": catalog_link("id-catalog-v4", "catalogs/id-catalog-v4.tsv", 'e')
        }
    })
}

fn artifact_link(id: &str, kind: &str, path: &str, digest: char) -> Value {
    json!({
        "artifact_id": id,
        "artifact_kind": kind,
        "schema_version": "1",
        "path": path,
        "size_bytes": 1,
        "sha256": digest.to_string().repeat(64)
    })
}

fn artifacts() -> Value {
    json!([
        artifact_link(
            "f5-01-block-layout-v1",
            "f5_block_layout",
            "outputs/f5-block-layout.tsv",
            '1'
        ),
        artifact_link("f5-02-sites-v1", "f5_sites", "outputs/f5-sites.tsv", '2'),
        artifact_link(
            "f5-03-haplotypes-v1",
            "f5_haplotypes",
            "outputs/f5-haplotypes.tsv",
            '3'
        ),
        artifact_link(
            "f5-04-frequencies-v1",
            "f5_frequencies",
            "outputs/f5-frequencies.tsv",
            '4'
        ),
        artifact_link(
            "f5-05-ld-evaluation-v1",
            "f5_ld_evaluation",
            "outputs/f5-ld-evaluation.tsv",
            '5'
        ),
        artifact_link(
            "f5-06-founders-v1",
            "f5_founders",
            "outputs/f5-founders.tsv",
            '6'
        ),
        artifact_link(
            "f5-07-attempts-v1",
            "f5_attempts",
            "outputs/f5-attempts.tsv",
            '7'
        )
    ])
}

fn failure_body(code: &str, limits: Value, artifact: bool) -> Value {
    json!({
        "code": code,
        "stage": "layout",
        "scenario_id": "scenario",
        "population_id": null,
        "block_id": null,
        "limits": limits,
        "counters": {
            "layout_attempts_used": 0,
            "haplotype_attempts_used": 0,
            "swap_proposals": 0,
            "swap_accepted": 0
        },
        "best_layout": null,
        "best_attempt": null,
        "artifacts": if artifact {
            json!([artifact_link(
                "f5-00-failure-v1",
                "f5_failure",
                "outputs/f5-failure.json",
                'f'
            )])
        } else {
            json!([])
        }
    })
}

fn failure_document(code: &str, limits: Value) -> Value {
    let mut value = failure_body(code, limits, false);
    value
        .as_object_mut()
        .expect("failure object")
        .insert("diagnostic_kind".to_owned(), json!("f5_failure"));
    value
        .as_object_mut()
        .expect("failure object")
        .insert("schema_version".to_owned(), json!("1"));
    value
}

fn success_manifest() -> Value {
    let dosage_numbers = json!([
        {"dosage": 0, "value": 0.25},
        {"dosage": 1, "value": 0.5},
        {"dosage": 2, "value": 0.25}
    ]);
    let dosage_counts = json!([
        {"dosage": 0, "value": 0},
        {"dosage": 1, "value": 1},
        {"dosage": 2, "value": 0}
    ]);
    json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": "5",
        "execution_id": "execution",
        "state": "complete",
        "f5": {
            "schema_version": "1",
            "scenario_id": "scenario",
            "seed": "b".repeat(64),
            "reference_id": "reference",
            "reference_sha256": "c".repeat(64),
            "profile": profile_link(),
            "mask": null,
            "eligible_bases": 1,
            "target": target(),
            "ploidy": 2,
            "founder_count": 1,
            "subpopulations": [{
                "population_id": "P1",
                "individual_count": 1,
                "inbreeding_f": 0.0,
                "inbreeding_f_source": "population"
            }],
            "algorithms": algorithms(),
            "rng_catalog": catalog_link("f5-rng-catalog-v1", "catalogs/f5-rng-catalog-v1.tsv", 'd'),
            "id_catalog": catalog_link("id-catalog-v4", "catalogs/id-catalog-v4.tsv", 'e'),
            "limits": limits(),
            "resource_preflight": preflight(),
            "normalizations": normalizations(),
            "attempt_summary": {
                "layout_attempts_used": 1,
                "haplotype_attempts_total": 1,
                "max_attempt_index": 0
            },
            "sampling": [],
            "quota_ties": [{
                "population_id": "P1",
                "block_id": "BLK_1",
                "remainder": 0.5,
                "candidate_ranks": [0, 1],
                "select_count": 1,
                "selected_ranks": [1]
            }],
            "f_realized": [{
                "population_id": "P1",
                "block_id": "BLK_1",
                "target_f": 0.0,
                "status": "evaluated",
                "value": 0.0
            }],
            "initialization_summaries": [{
                "population_id": "P1",
                "block_id": "BLK_1",
                "site_id": "SITE_1",
                "ploidy": 2,
                "model": "hwe",
                "requested_mode": "probabilistic",
                "effective_mode": "probabilistic",
                "input_alt_af": 0.5,
                "realized_alt_af": 0.5,
                "expected_genotype_probabilities": dosage_numbers,
                "expected_genotype_frequencies": dosage_numbers,
                "expected_genotype_counts": dosage_numbers,
                "realized_genotype_frequencies": dosage_numbers,
                "realized_genotype_counts": dosage_counts,
                "recounted_genotype_frequencies": dosage_numbers,
                "recounted_genotype_counts": dosage_counts,
                "recount_source": "outputs/truth.jsonl#f5-founder-truth-projection-v2",
                "recount_matches": true
            }],
            "truth": {
                "artifact_id": "truth-jsonl",
                "artifact_kind": "truth_jsonl",
                "schema_version": "2",
                "projection_version": "f5-founder-truth-projection-v2",
                "path": "outputs/truth.jsonl",
                "size_bytes": 1,
                "sha256": "8".repeat(64)
            },
            "artifacts": artifacts()
        }
    })
}

fn failed_manifest(code: &str, limits: Value) -> Value {
    json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": "5",
        "execution_id": "execution",
        "state": "failed",
        "f5": {
            "schema_version": "1",
            "f5_failure": failure_body(code, limits, true)
        }
    })
}

#[test]
fn externally_validates_general_scenario_and_profile_in_isolated_bundles() {
    let general = json!({
        "config_kind": "general",
        "schema_version": "3",
        "seed": "a".repeat(64),
        "ploidy": 2,
        "reference_id": "reference",
        "reference_path": "references/reference.fa",
        "reference_sha256": "b".repeat(64),
        "f5_scenario": {
            "scenario_id": "scenario",
            "scenario_path": "scenarios/scenario.toml"
        }
    });
    let mut invalid_general = general.clone();
    invalid_general["unknown"] = json!(true);
    assert_external_cases(
        "config-general-v3.schema.json",
        &[],
        vec![(general, true), (invalid_general, false)],
    );

    let scenario = json!({
        "config_kind": "f5_scenario",
        "schema_version": "1",
        "scenario_id": "scenario",
        "founder_count": 1,
        "f5": {
            "profile_metadata_path": "profiles/profile.metadata.json",
            "profile_metadata_sha256": "c".repeat(64),
            "block_count": 1
        },
        "subpopulations": [{
            "population_id": "P1",
            "individual_count": 1,
            "founder_ids": ["I1"]
        }]
    });
    let mut invalid_scenario = scenario.clone();
    invalid_scenario["f5"]["block_count"] = json!(0);
    assert_external_cases(
        "config-f5-scenario-v1.schema.json",
        &[],
        vec![(scenario, true), (invalid_scenario, false)],
    );

    let profile = json!({
        "profile_kind": "f5_generation_profile",
        "schema_version": "1",
        "profile_id": "profile",
        "provenance_kind": "built_in",
        "source_citation": "citation",
        "source_license": null,
        "source_sha256": null,
        "files": [
            {"artifact_kind":"f5_block_length_profile","schema_version":"1","path":"block_lengths.tsv","size_bytes":1,"sha256":"a".repeat(64)},
            {"artifact_kind":"f5_snp_gap_profile","schema_version":"1","path":"snp_gaps.tsv","size_bytes":1,"sha256":"b".repeat(64)},
            {"artifact_kind":"f5_snp_spectrum","schema_version":"1","path":"snp_spectrum.tsv","size_bytes":1,"sha256":"c".repeat(64)},
            {"artifact_kind":"f5_haplotype_weights","schema_version":"1","path":"haplotype_weights.tsv","size_bytes":1,"sha256":"d".repeat(64)},
            {"artifact_kind":"f5_ld_profile","schema_version":"1","path":"ld_profile.tsv","size_bytes":1,"sha256":"e".repeat(64)}
        ]
    });
    let mut invalid_profile = profile.clone();
    invalid_profile["files"][0]["path"] = json!("snp_gaps.tsv");
    assert_external_cases(
        "profile-metadata-f5-v1.schema.json",
        &[],
        vec![(profile, true), (invalid_profile, false)],
    );
}

#[test]
fn authoring_schemas_match_u64_and_normalized_path_boundaries() {
    let general = json!({
        "config_kind": "general",
        "schema_version": "3",
        "ploidy": 2,
        "reference_id": "reference",
        "reference_path": "references/reference.fa",
        "reference_sha256": "b".repeat(64),
        "f5_scenario": {
            "scenario_id": "scenario",
            "scenario_path": "scenarios/scenario.toml"
        }
    });
    let scenario = json!({
        "config_kind": "f5_scenario",
        "schema_version": "1",
        "scenario_id": "scenario",
        "founder_count": 1,
        "f5": {
            "profile_metadata_path": "profiles/profile.metadata.json",
            "profile_metadata_sha256": "c".repeat(64),
            "excluded_regions_path": "masks/excluded.bed",
            "excluded_regions_sha256": "d".repeat(64),
            "covered_bases": 1
        },
        "subpopulations": [{
            "population_id": "P1",
            "individual_count": 1,
            "founder_ids": ["I1"]
        }]
    });
    let invalid_paths = [
        "scenarios//scenario.toml",
        "scenarios/./scenario.toml",
        "scenarios/.",
        "scenarios/",
        "scenarios/../scenario.toml",
        "/scenarios/scenario.toml",
        "scenarios\\scenario.toml",
        "scenarios/scenario file.toml",
    ];

    let mut general_cases = vec![(general.clone(), true)];
    for invalid in invalid_paths {
        let mut candidate = general.clone();
        candidate["f5_scenario"]["scenario_path"] = json!(invalid);
        general_cases.push((candidate, false));
    }
    assert_external_cases("config-general-v3.schema.json", &[], general_cases);

    let mut scenario_cases = vec![(scenario.clone(), true)];
    for invalid in invalid_paths {
        for field in ["profile_metadata_path", "excluded_regions_path"] {
            let mut candidate = scenario.clone();
            candidate["f5"][field] = json!(invalid);
            scenario_cases.push((candidate, false));
        }
    }
    let max = serde_json::from_str::<Value>("18446744073709551615").expect("u64::MAX JSON");
    let max_plus_one =
        serde_json::from_str::<Value>("18446744073709551616").expect("u64::MAX + 1 JSON");
    for pointer in [
        "/founder_count",
        "/f5/covered_bases",
        "/subpopulations/0/individual_count",
    ] {
        let mut accepted = scenario.clone();
        *accepted.pointer_mut(pointer).expect("u64 field") = max.clone();
        scenario_cases.push((accepted, true));

        let mut rejected = scenario.clone();
        *rejected.pointer_mut(pointer).expect("u64 field") = max_plus_one.clone();
        scenario_cases.push((rejected, false));
    }
    assert_external_cases("config-f5-scenario-v1.schema.json", &[], scenario_cases);
}

#[test]
fn effective_profile_and_catalog_links_are_closed_and_instantiable() {
    let valid = effective();
    let mut cases = vec![(valid.clone(), true)];
    for field in [
        "fasta_bytes_processed",
        "mask_bytes_processed",
        "general_toml_bytes",
        "scenario_toml_bytes",
        "profile_metadata_bytes",
        "profile_tsv_bytes",
        "input_parser_peak_bytes",
    ] {
        let mut missing = valid.clone();
        missing["f5"]["resource_preflight"]
            .as_object_mut()
            .expect("preflight object")
            .remove(field);
        cases.push((missing, false));
    }
    let mut stale_algorithm = valid.clone();
    stale_algorithm["f5"]["resource_preflight"]["algorithm"] = json!("f5-resource-accounting-v5");
    cases.push((stale_algorithm, false));
    for (object, id_field) in [
        ("profile", "profile_id"),
        ("rng_catalog", "catalog_id"),
        ("id_catalog", "catalog_id"),
    ] {
        for (field, invalid) in [
            (id_field, json!("")),
            ("path", json!("../escape")),
            ("size_bytes", json!(-1)),
            ("sha256", json!("not-a-sha")),
        ] {
            let mut mutation = valid.clone();
            mutation["f5"][object][field] = invalid;
            cases.push((mutation, false));
        }
    }
    assert_external_cases("effective-config-v5.schema.json", &[], cases);
}

#[test]
fn failure_limits_and_codes_use_exact_f5_contract() {
    let bounds = [
        ("max_layout_attempts", 100_000_u64),
        ("max_haplotype_attempts", 10_000),
        ("max_swap_proposals", 10_000_000),
        ("max_haplotype_cells", 67_108_864),
        ("max_operation_units", 268_435_456),
        ("max_working_memory_bytes", 1_073_741_824),
        ("max_output_bytes", 2_147_483_648),
    ];
    let mut cases = vec![(failure_document("GPS-CONF-101", limits()), true)];
    for (field, maximum) in bounds {
        for (value, expected) in [(0, false), (maximum, true), (maximum + 1, false)] {
            let mut mutation = failure_document("GPS-CONF-101", limits());
            mutation["limits"][field] = json!(value);
            cases.push((mutation, expected));
        }
    }
    for code in [
        "GPS-CONF-101",
        "GPS-FMT-101",
        "GPS-DOMAIN-101",
        "GPS-DOMAIN-102",
        "GPS-DOMAIN-103",
        "GPS-DOMAIN-104",
        "GPS-DOMAIN-105",
        "GPS-DOMAIN-106",
        "GPS-RES-101",
        "GPS-INTERNAL-101",
    ] {
        cases.push((failure_document(code, limits()), true));
    }
    cases.push((failure_document("GPS-CONF-999", limits()), false));
    assert_external_cases("f5-failure-v1.schema.json", &[], cases);
}

#[test]
fn manifest_success_and_failure_branches_are_closed_and_literal() {
    let success = success_manifest();
    let failure = failed_manifest("GPS-DOMAIN-104", limits());
    let mut cases = vec![(success.clone(), true), (failure.clone(), true)];

    let mut wrong_truth_version = success.clone();
    wrong_truth_version["f5"]["truth"]["schema_version"] = json!("1");
    cases.push((wrong_truth_version, false));
    let mut open_truth = success.clone();
    open_truth["f5"]["truth"]["unknown"] = json!(true);
    cases.push((open_truth, false));

    for collection in ["quota_ties", "f_realized", "initialization_summaries"] {
        let mut mutation = success.clone();
        mutation["f5"][collection][0]["unknown"] = json!(true);
        cases.push((mutation, false));
    }
    for index in 0..7 {
        let mut mutation = success.clone();
        mutation["f5"]["artifacts"][index]["artifact_id"] = json!("wrong");
        cases.push((mutation, false));
    }

    let bounds = [
        ("max_layout_attempts", 100_000_u64),
        ("max_haplotype_attempts", 10_000),
        ("max_swap_proposals", 10_000_000),
        ("max_haplotype_cells", 67_108_864),
        ("max_operation_units", 268_435_456),
        ("max_working_memory_bytes", 1_073_741_824),
        ("max_output_bytes", 2_147_483_648),
    ];
    for (field, maximum) in bounds {
        for (value, expected) in [(0, false), (maximum, true), (maximum + 1, false)] {
            let mut mutation = failure.clone();
            mutation["f5"]["f5_failure"]["limits"][field] = json!(value);
            cases.push((mutation, expected));
        }
    }
    let mut unknown_code = failure;
    unknown_code["f5"]["f5_failure"]["code"] = json!("GPS-DOMAIN-999");
    cases.push((unknown_code, false));

    assert_external_cases(
        "manifest-general-v5.schema.json",
        &[
            "effective-config-v5.schema.json",
            "f5-failure-v1.schema.json",
        ],
        cases,
    );
}

#[test]
fn truth_v2_validates_standalone_with_historical_closed_definitions() {
    let schema = include_str!("../schemas/truth-v2.schema.json");
    assert!(
        !schema.contains("truth-v1.schema.json"),
        "truth v2 must not resolve its historical definitions from a sibling"
    );

    let lines = include_str!("fixtures/truth-v1/valid.jsonl")
        .lines()
        .collect::<Vec<_>>();
    let mut header: Value = serde_json::from_str(lines[0]).expect("header");
    header["schema_version"] = json!("2");
    header["projection"] = json!({
        "kind": "f5_founder",
        "projection_version": "f5-founder-truth-projection-v2",
        "scenario_id": "scenario"
    });
    header["populations"] = json!([{"population_id":"P1","individual_count":1,"ploidy":2}]);

    let mut copy: Value = serde_json::from_str(lines[1]).expect("copy");
    copy["founder_membership"] = json!({"population_id":"P1","founder_copy_id":"FCP_01"});

    let mut genotype: Value = serde_json::from_str(lines[3]).expect("genotype");
    genotype["f5_locus"] = json!({
        "population_id": "P1",
        "block_id": "BLK_01",
        "site_id": "SITE_01",
        "ref": "A",
        "alt": "G"
    });

    let mut invalid_header = header.clone();
    invalid_header["projection"]["kind"] = json!("other");
    let mut invalid_copy = copy.clone();
    invalid_copy["founder_membership"]["unknown"] = json!(true);
    let mut invalid_historical_copy = copy.clone();
    invalid_historical_copy["identity"]
        .as_object_mut()
        .expect("identity")
        .remove("reference_id");
    let mut invalid_genotype = genotype.clone();
    invalid_genotype["f5_locus"]["alt"] = json!("N");

    assert_external_cases(
        "truth-v2.schema.json",
        &[],
        vec![
            (header, true),
            (copy, true),
            (genotype, true),
            (invalid_header, false),
            (invalid_copy, false),
            (invalid_historical_copy, false),
            (invalid_genotype, false),
        ],
    );
}
