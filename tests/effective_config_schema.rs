use jsonschema::{Draft, Validator};
use serde_json::{Value, json};

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
fn every_effective_artifact_uses_the_canonical_relative_path_pattern() {
    for path in [
        "references/main/reference.fasta",
        "a",
        "a/.hidden",
        "a/name..txt",
        "referencias/ñ/reference.fasta",
    ] {
        let mut document = effective_config("supplied", true);
        *document
            .pointer_mut("/reference/fasta/path")
            .expect("reference path") = json!(path);
        assert_valid(&format!("canonical-path-{path:?}"), &document);
    }

    for path in [
        "",
        "/absolute",
        "../outside",
        "a/../outside",
        "./relative",
        "a/./relative",
        "a//duplicate",
        "a\\windows",
        "a\nb",
        "a\u{7f}b",
        "trailing/",
    ] {
        let mut document = effective_config("supplied", true);
        *document
            .pointer_mut("/reference/fasta/path")
            .expect("reference path") = json!(path);
        assert_invalid(&format!("noncanonical-path-{path:?}"), &document);
    }
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
            "fasta": artifact("reference-reduced", "references/reduced/reference.fasta"),
            "map": artifact("reference-map", "references/reduced/map.tsv"),
            "source_fasta_id": "reference-main"
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
            "total_genome_size": 1000,
            "composition": composition,
            "sequences": [{
                "id": "chr1",
                "class": "chromosome",
                "length": 1000,
                "unknown_base_count": 0
            }],
            "fasta": artifact("reference-main", "references/main/reference.fasta")
        },
        "reduced_references": reductions
    })
}

fn artifact(id: &str, path: &str) -> Value {
    json!({
        "id": id,
        "kind": "reference_fasta",
        "path": path,
        "size_bytes": 100,
        "checksum": {"algorithm": "sha256", "digest": "0".repeat(64)},
        "state": "complete",
        "validated": true
    })
}

fn schema_validator() -> Validator {
    let schema: Value =
        serde_json::from_str(include_str!("../schemas/effective-config-v1.schema.json"))
            .expect("effective-config schema JSON");
    jsonschema::options()
        .with_draft(Draft::Draft202012)
        .build(&schema)
        .expect("Draft 2020-12 schema")
}

fn assert_valid(name: &str, document: &Value) {
    assert!(
        schema_validator().is_valid(document),
        "effective-config schema rejected positive {name}: {document}"
    );
}

fn assert_invalid(name: &str, document: &Value) {
    assert!(
        !schema_validator().is_valid(document),
        "effective-config schema accepted negative {name}: {document}"
    );
}
