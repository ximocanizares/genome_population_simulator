use std::collections::BTreeSet;

use serde_json::Value;

const AUTHORING_SCHEMA: &str = include_str!("../schemas/config-general-v1.schema.json");
const DIAGNOSTIC_SCHEMA: &str = include_str!("../schemas/diagnostics-v1.schema.json");
const DIAGNOSTIC_CATALOG: &str = include_str!("../catalogs/diagnostics-v2.tsv");

fn keys_at(document: &Value, pointer: &str) -> BTreeSet<String> {
    document
        .pointer(pointer)
        .and_then(Value::as_object)
        .unwrap_or_else(|| panic!("missing object at {pointer}"))
        .keys()
        .cloned()
        .collect()
}

fn strings_at(document: &Value, pointer: &str) -> BTreeSet<String> {
    document
        .pointer(pointer)
        .and_then(Value::as_array)
        .unwrap_or_else(|| panic!("missing array at {pointer}"))
        .iter()
        .map(|value| {
            value
                .as_str()
                .unwrap_or_else(|| panic!("non-string value at {pointer}"))
                .to_owned()
        })
        .collect()
}

#[test]
fn f1_authoring_schema_exposes_exactly_the_fields_accepted_by_rust() {
    let schema: Value = serde_json::from_str(AUTHORING_SCHEMA).expect("valid JSON schema");

    assert_eq!(
        keys_at(&schema, "/properties"),
        BTreeSet::from([
            "config_kind".to_owned(),
            "reduced_references".to_owned(),
            "reference".to_owned(),
            "schema_version".to_owned(),
            "seed".to_owned(),
        ])
    );
    assert_eq!(
        strings_at(&schema, "/required"),
        BTreeSet::from([
            "config_kind".to_owned(),
            "reference".to_owned(),
            "schema_version".to_owned(),
        ])
    );
    assert_eq!(
        keys_at(&schema, "/$defs/simulated_reference/properties"),
        BTreeSet::from([
            "composition".to_owned(),
            "id".to_owned(),
            "kind".to_owned(),
            "sequences".to_owned(),
        ])
    );
    assert_eq!(
        keys_at(&schema, "/$defs/supplied_reference/properties"),
        BTreeSet::from([
            "fasta".to_owned(),
            "id".to_owned(),
            "kind".to_owned(),
            "sequences".to_owned(),
        ])
    );
    assert_eq!(
        keys_at(&schema, "/$defs/sequence/properties"),
        BTreeSet::from(["class".to_owned(), "id".to_owned(), "length".to_owned()])
    );
    assert_eq!(
        keys_at(&schema, "/$defs/composition/properties"),
        BTreeSet::from([
            "A".to_owned(),
            "C".to_owned(),
            "G".to_owned(),
            "T".to_owned(),
        ])
    );
    assert_eq!(
        keys_at(&schema, "/$defs/reduced_reference/properties"),
        BTreeSet::from([
            "id".to_owned(),
            "reduction_percent".to_owned(),
            "region_size_bases".to_owned(),
            "source_reference_id".to_owned(),
            "target_bases".to_owned(),
        ])
    );
    assert_eq!(
        strings_at(&schema, "/$defs/simulated_reference/required"),
        BTreeSet::from([
            "composition".to_owned(),
            "id".to_owned(),
            "kind".to_owned(),
            "sequences".to_owned(),
        ])
    );
    assert_eq!(
        strings_at(&schema, "/$defs/supplied_reference/required"),
        BTreeSet::from([
            "fasta".to_owned(),
            "id".to_owned(),
            "kind".to_owned(),
            "sequences".to_owned(),
        ])
    );
    assert_eq!(
        strings_at(&schema, "/$defs/sequence/required"),
        BTreeSet::from(["class".to_owned(), "id".to_owned(), "length".to_owned()])
    );
    assert_eq!(
        strings_at(&schema, "/$defs/composition/required"),
        BTreeSet::from([
            "A".to_owned(),
            "C".to_owned(),
            "G".to_owned(),
            "T".to_owned(),
        ])
    );
    assert_eq!(
        strings_at(&schema, "/$defs/reduced_reference/required"),
        BTreeSet::from([
            "id".to_owned(),
            "region_size_bases".to_owned(),
            "source_reference_id".to_owned(),
        ])
    );
    let target_modes = schema["$defs"]["reduced_reference"]["oneOf"]
        .as_array()
        .expect("reduction target union");
    assert_eq!(target_modes.len(), 2);
    assert_eq!(
        strings_at(&schema, "/$defs/reduced_reference/oneOf/0/required"),
        BTreeSet::from(["target_bases".to_owned()])
    );
    assert_eq!(
        strings_at(&schema, "/$defs/reduced_reference/oneOf/1/required"),
        BTreeSet::from(["reduction_percent".to_owned()])
    );

    for future_field in [
        "scenario_id",
        "scenarios",
        "route",
        "ploidy",
        "individuals",
        "populations",
    ] {
        assert!(
            !AUTHORING_SCHEMA.contains(&format!("\"{future_field}\"")),
            "future F2 field leaked into F2.0.1: {future_field}"
        );
    }
}

#[test]
fn authoring_resources_are_versioned_and_keep_rust_as_semantic_authority() {
    let schema: Value = serde_json::from_str(AUTHORING_SCHEMA).expect("valid authoring schema");
    let diagnostics: Value =
        serde_json::from_str(DIAGNOSTIC_SCHEMA).expect("valid diagnostic schema");

    assert_eq!(schema["properties"]["config_kind"]["const"], "general");
    assert_eq!(schema["properties"]["schema_version"]["const"], "1");
    assert_eq!(schema["additionalProperties"], false);
    assert!(schema["description"].as_str().is_some_and(|description| {
        description.contains("binario Rust conserva la validación semántica")
    }));
    assert_eq!(diagnostics["properties"]["schema_version"]["const"], "1");
    assert!(DIAGNOSTIC_CATALOG.contains("# diagnostic_catalog_version=2"));

    for path in [
        "/properties/seed/x-gps-control",
        "/properties/reference/x-gps-tab",
        "/properties/reduced_references/x-gps-control",
        "/$defs/sequence/properties/class/x-gps-control",
        "/$defs/sequence/properties/length/x-gps-unit",
    ] {
        assert!(
            schema.pointer(path).is_some(),
            "missing UI annotation at {path}"
        );
    }
}
