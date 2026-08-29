const GENERAL_SCHEMA: &str = include_str!("../schemas/config-general-v2.md");
const SCENARIO_SCHEMA: &str = include_str!("../schemas/config-scenario-v1.md");

#[test]
fn active_authoring_schema_is_v2_and_scenario_schema_is_documental_only() {
    for literal in [
        "config_kind = \"general\"",
        "schema_version = \"2\"",
        "La raíz es cerrada",
        "`reference` es obligatorio",
        "`reduced_references`",
        "Evento técnico F4",
        "target_coverage_wgs",
        "read_count` no forma parte de v2 y se rechaza",
    ] {
        assert!(GENERAL_SCHEMA.contains(literal), "missing {literal}");
    }
    assert!(GENERAL_SCHEMA.contains("library_id` se deriva"));
    assert!(SCENARIO_SCHEMA.contains("config_kind = \"scenario\""));
    assert!(SCENARIO_SCHEMA.contains("No forma parte del contrato ejecutable de F1"));
    assert!(!SCENARIO_SCHEMA.contains("F1.2"));
    assert!(!SCENARIO_SCHEMA.contains("admitidas en F1"));
}
