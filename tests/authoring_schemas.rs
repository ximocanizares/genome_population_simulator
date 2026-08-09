const GENERAL_SCHEMA: &str = include_str!("../schemas/config-general-v1.md");
const SCENARIO_SCHEMA: &str = include_str!("../schemas/config-scenario-v1.md");

#[test]
fn authoring_schemas_fix_class_version_and_explicit_list_identity() {
    for schema in [GENERAL_SCHEMA, SCENARIO_SCHEMA] {
        assert!(schema.contains("config_kind"));
        assert!(schema.contains("schema_version = \"1\""));
        assert!(schema.contains("`id` explícito"));
        assert!(schema.contains("campos desconocidos") || schema.contains("claves desconocidas"));
    }
}

#[test]
fn valid_examples_cover_zero_one_many_events_and_both_scenario_routes() {
    let zero = include_str!("fixtures/config/general-minimal.toml");
    let many = include_str!("fixtures/config/general-events.toml");
    let evolution = include_str!("fixtures/config/scenario-evolution.toml");
    let families = include_str!("fixtures/config/scenario-families.toml");
    let supplied = include_str!("fixtures/config/general-supplied.toml");
    let evolution_wrapper = include_str!("fixtures/config/general-evolution.toml");
    let families_wrapper = include_str!("fixtures/config/general-families.toml");

    assert!(!zero.contains("[[events]]"));
    assert_eq!(many.matches("[[events]]").count(), 2);
    assert!(evolution.contains("route = \"evolution\""));
    assert!(families.contains("route = \"families\""));
    assert!(supplied.contains("kind = \"supplied\""));
    assert!(evolution_wrapper.contains("path = \"scenario-evolution.toml\""));
    assert!(families_wrapper.contains("path = \"scenario-families.toml\""));
}

#[test]
fn invalid_examples_cover_header_boundary_without_partial_semantics() {
    let fixtures = [
        include_str!("fixtures/config/invalid-missing-kind.toml"),
        include_str!("fixtures/config/invalid-missing-version.toml"),
        include_str!("fixtures/config/invalid-unknown-kind.toml"),
        include_str!("fixtures/config/invalid-crossed-kind.toml"),
        include_str!("fixtures/config/invalid-version-zero.toml"),
        include_str!("fixtures/config/invalid-version-two.toml"),
        include_str!("fixtures/config/invalid-version-integer.toml"),
        include_str!("fixtures/config/invalid-malformed.toml"),
    ];
    assert_eq!(fixtures.len(), 8);
}
