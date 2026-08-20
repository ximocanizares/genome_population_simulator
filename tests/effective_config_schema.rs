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
fn effective_selection_order_is_strictly_one_based() {
    let mut document = effective_config("supplied", true);
    *document
        .pointer_mut("/reduced_references/0/selection_order/0")
        .expect("selection order") = json!(0);
    assert_invalid("zero-selection-order", &document);
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
            "selection_subseed": "1".repeat(64),
            "selection_order": [1, 2],
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
        "schema_version": "1",
        "source_path": "config/general.toml",
        "seed": {
            "root": "0".repeat(64),
            "origin": "declared"
        },
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
        "reduced_references": reductions
    })
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
        .arg(root.join("schemas/effective-config-v1.schema.json"))
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
