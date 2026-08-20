use std::path::Path;

use genome_population_simulator::results::validate_general_manifest_bytes;
use jsonschema::{Draft, Validator};
use serde_json::{Value, json};

#[test]
fn complete_and_failed_minimal_manifests_pass_both_validators() {
    for (name, document) in [
        ("complete", complete_manifest()),
        ("failed", failed_manifest()),
    ] {
        assert_accepted(name, &document);
    }
}

#[test]
fn exact_reviewer_bypass_is_rejected_by_both_validators() {
    let bypass = json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": "1",
        "execution_id": "",
        "state": "complete",
        "configuration_input": {},
        "effective_config": {},
        "references": [42],
        "reference_reductions": [false],
        "artifacts": [],
        "units": ["not-unit"],
        "diagnostics": [17]
    });
    assert_rejected("exact-reviewer-bypass", &bypass);
}

#[test]
fn every_root_field_and_collection_has_a_negative_witness() {
    let base = complete_manifest();
    for field in [
        "manifest_kind",
        "schema_version",
        "execution_id",
        "state",
        "configuration_input",
        "effective_config",
        "references",
        "reference_reductions",
        "artifacts",
        "units",
        "diagnostics",
    ] {
        let mut document = base.clone();
        document.as_object_mut().expect("root").remove(field);
        assert_rejected(&format!("missing-root-{field}"), &document);
    }

    for (name, pointer, replacement) in [
        ("manifest-kind", "/manifest_kind", json!("other")),
        ("schema-version", "/schema_version", json!(2)),
        ("execution-id", "/execution_id", json!("")),
        ("state", "/state", json!("incomplete")),
        ("configuration-input", "/configuration_input", json!([])),
        ("effective-config", "/effective_config", json!({})),
        ("references", "/references", json!({})),
        ("reference-reductions", "/reference_reductions", json!({})),
        ("artifacts", "/artifacts", json!({})),
        ("units", "/units", json!({})),
        ("diagnostics", "/diagnostics", json!({})),
    ] {
        let mut document = base.clone();
        *document.pointer_mut(pointer).expect("root pointer") = replacement;
        assert_rejected(name, &document);
    }

    for (name, pointer, replacement) in [
        ("reference-item", "/references/0", json!(42)),
        ("reduction-item", "/reference_reductions/0", json!(false)),
        ("artifact-item", "/artifacts/0", json!("not-artifact")),
        ("unit-item", "/units/0", json!("not-unit")),
        ("diagnostic-item", "/diagnostics/0", json!(17)),
        (
            "selection-order-item",
            "/reference_reductions/0/selection_order/0",
            json!(0),
        ),
    ] {
        let mut document = base.clone();
        *document.pointer_mut(pointer).expect("collection pointer") = replacement;
        assert_rejected(name, &document);
    }
}

#[test]
fn every_nested_required_field_has_a_negative_witness() {
    let base = complete_manifest();
    for (object_pointer, fields) in [
        (
            "/configuration_input",
            &["kind", "original_path", "copy"][..],
        ),
        (
            "/configuration_input/copy",
            &[
                "id",
                "kind",
                "path",
                "size_bytes",
                "checksum",
                "state",
                "validated",
            ][..],
        ),
        (
            "/configuration_input/copy/checksum",
            &["algorithm", "digest"][..],
        ),
        ("/references/0", &["id", "fasta"][..]),
        (
            "/reference_reductions/0",
            &[
                "id",
                "source_reference_id",
                "region_size_bases",
                "target",
                "selection_subseed",
                "selection_order",
                "state",
                "fasta",
                "map",
                "source_fasta_id",
            ][..],
        ),
        (
            "/reference_reductions/0/target",
            &["mode", "declared", "bases"][..],
        ),
        (
            "/units/0",
            &["id", "kind", "state", "exportable", "diagnostic_code"][..],
        ),
        ("/diagnostics/0", &["code", "message"][..]),
    ] {
        for field in fields {
            let mut document = base.clone();
            document
                .pointer_mut(object_pointer)
                .and_then(Value::as_object_mut)
                .expect("object pointer")
                .remove(*field);
            assert_rejected(
                &format!("missing-{}-{field}", object_pointer.replace('/', "-")),
                &document,
            );
        }
    }
}

#[test]
fn every_object_kind_is_closed_and_field_types_are_concrete() {
    let base = complete_manifest();
    for pointer in [
        "",
        "/configuration_input",
        "/configuration_input/copy",
        "/configuration_input/copy/checksum",
        "/references/0",
        "/reference_reductions/0",
        "/reference_reductions/0/target",
        "/units/0",
        "/diagnostics/0",
    ] {
        let mut document = base.clone();
        let object = if pointer.is_empty() {
            document.as_object_mut().expect("root")
        } else {
            document
                .pointer_mut(pointer)
                .and_then(Value::as_object_mut)
                .expect("object pointer")
        };
        object.insert("unexpected".to_owned(), json!(true));
        assert_rejected(
            &format!("unknown-field-{}", pointer.replace('/', "-")),
            &document,
        );
    }

    for (name, pointer, replacement) in [
        (
            "configuration-kind",
            "/configuration_input/kind",
            json!("scenario"),
        ),
        (
            "configuration-path",
            "/configuration_input/original_path",
            json!(7),
        ),
        ("artifact-id", "/artifacts/0/id", json!("")),
        ("artifact-kind", "/artifacts/0/kind", json!(false)),
        ("artifact-path", "/artifacts/0/path", json!("")),
        ("artifact-size", "/artifacts/0/size_bytes", json!(-1)),
        (
            "checksum-algorithm",
            "/artifacts/0/checksum/algorithm",
            json!("md5"),
        ),
        (
            "checksum-digest",
            "/artifacts/0/checksum/digest",
            json!("A".repeat(64)),
        ),
        ("artifact-state", "/artifacts/0/state", json!("incomplete")),
        ("artifact-validated", "/artifacts/0/validated", json!("yes")),
        ("reference-id", "/references/0/id", json!(0)),
        (
            "reduction-source",
            "/reference_reductions/0/source_reference_id",
            json!(null),
        ),
        (
            "reduction-region",
            "/reference_reductions/0/region_size_bases",
            json!(0),
        ),
        (
            "target-mode",
            "/reference_reductions/0/target/mode",
            json!("other"),
        ),
        (
            "target-declared",
            "/reference_reductions/0/target/declared",
            json!(1),
        ),
        (
            "target-bases",
            "/reference_reductions/0/target/bases",
            json!(0),
        ),
        (
            "subseed",
            "/reference_reductions/0/selection_subseed",
            json!("A".repeat(64)),
        ),
        (
            "reduction-state",
            "/reference_reductions/0/state",
            json!(true),
        ),
        (
            "source-fasta-id",
            "/reference_reductions/0/source_fasta_id",
            json!(""),
        ),
        ("unit-id", "/units/0/id", json!(null)),
        ("unit-kind", "/units/0/kind", json!("")),
        ("unit-state", "/units/0/state", json!("other")),
        ("unit-exportable", "/units/0/exportable", json!(1)),
        ("unit-diagnostic", "/units/0/diagnostic_code", json!("BAD")),
        ("diagnostic-code", "/diagnostics/0/code", json!("BAD")),
        ("diagnostic-message", "/diagnostics/0/message", json!("")),
    ] {
        let mut document = base.clone();
        *document.pointer_mut(pointer).expect("field pointer") = replacement;
        assert_rejected(name, &document);
    }

    let mut incompatible_artifact = base;
    *incompatible_artifact
        .pointer_mut("/artifacts/0/validated")
        .expect("validated") = json!(false);
    assert_rejected("artifact-state-validation-mismatch", &incompatible_artifact);
}

#[test]
fn artifact_paths_have_the_same_canonical_contract_in_schema_and_rust() {
    for path in [
        "inputs/config/general.toml",
        "a",
        "a/.hidden",
        "a/name..txt",
        "referencias/ñ/reference.fasta",
    ] {
        let mut document = complete_manifest();
        *document
            .pointer_mut("/artifacts/0/path")
            .expect("artifact path") = json!(path);
        assert_accepted(&format!("canonical-path-{path:?}"), &document);
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
        let mut document = complete_manifest();
        *document
            .pointer_mut("/artifacts/0/path")
            .expect("artifact path") = json!(path);
        assert_rejected(&format!("noncanonical-path-{path:?}"), &document);
    }
}

fn complete_manifest() -> Value {
    let config = artifact(
        "config-general",
        "config_general_toml",
        "inputs/config/general.toml",
    );
    let reference = artifact(
        "reference-fasta-ref",
        "reference_fasta",
        "references/ref/reference.fasta",
    );
    let reduced = artifact(
        "reference-fasta-reduced",
        "reference_fasta",
        "references/reduced/reference.fasta",
    );
    let map = artifact(
        "reference-map-reduced",
        "reference_reduction_map",
        "references/reduced/map.tsv",
    );
    let effective = artifact(
        "effective-config",
        "effective_config",
        "effective-config.json",
    );
    json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": "1",
        "execution_id": "execution-001",
        "state": "complete",
        "configuration_input": {
            "kind": "general",
            "original_path": "config/general.toml",
            "copy": config
        },
        "effective_config": effective,
        "references": [
            {"id": "ref", "fasta": reference},
            {"id": "reduced", "fasta": reduced}
        ],
        "reference_reductions": [{
            "id": "reduced",
            "source_reference_id": "ref",
            "region_size_bases": 100,
            "target": {"mode": "target_bases", "declared": "200", "bases": 200},
            "selection_subseed": "1".repeat(64),
            "selection_order": [1, 2],
            "state": "complete",
            "fasta": reduced,
            "map": map,
            "source_fasta_id": "reference-fasta-ref"
        }],
        "artifacts": [config, reference, reduced, map, effective],
        "units": [{
            "id": "execution-001",
            "kind": "execution",
            "state": "complete",
            "exportable": true,
            "diagnostic_code": null
        }],
        "diagnostics": [{
            "code": "GPS-FMT-006",
            "message": "Testigo estructurado de diagnóstico."
        }]
    })
}

fn failed_manifest() -> Value {
    let config = artifact(
        "config-general",
        "config_general_toml",
        "inputs/config/general.toml",
    );
    json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": "1",
        "execution_id": "execution-failed-001",
        "state": "failed",
        "configuration_input": {
            "kind": "general",
            "original_path": "config/general.toml",
            "copy": config
        },
        "effective_config": null,
        "references": [],
        "reference_reductions": [],
        "artifacts": [config],
        "units": [{
            "id": "execution-failed-001",
            "kind": "execution",
            "state": "failed",
            "exportable": false,
            "diagnostic_code": "GPS-IO-001"
        }],
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
        "checksum": {"algorithm": "sha256", "digest": "0".repeat(64)},
        "state": "complete",
        "validated": true
    })
}

fn schema_validator() -> Validator {
    let schema: Value =
        serde_json::from_str(include_str!("../schemas/manifest-general-v1.schema.json"))
            .expect("manifest schema JSON");
    jsonschema::options()
        .with_draft(Draft::Draft202012)
        .build(&schema)
        .expect("Draft 2020-12 schema")
}

fn assert_accepted(name: &str, document: &Value) {
    assert!(
        schema_validator().is_valid(document),
        "JSON Schema rejected positive {name}"
    );
    let bytes = serde_json::to_vec(document).expect("manifest JSON");
    validate_general_manifest_bytes(Path::new("manifest.json"), &bytes)
        .unwrap_or_else(|error| panic!("Rust validator rejected positive {name}: {error:?}"));
}

fn assert_rejected(name: &str, document: &Value) {
    assert!(
        !schema_validator().is_valid(document),
        "JSON Schema accepted negative {name}: {document}"
    );
    let bytes = serde_json::to_vec(document).expect("manifest JSON");
    assert!(
        validate_general_manifest_bytes(Path::new("manifest.json"), &bytes).is_err(),
        "Rust validator accepted negative {name}: {document}"
    );
}
