const GENERAL_SCHEMA: &str = include_str!("../schemas/config-general-v1.md");
const SCENARIO_SCHEMA: &str = include_str!("../schemas/config-scenario-v1.md");

#[test]
fn f2_2_schema_is_closed_and_scenario_schema_is_documental_only() {
    for literal in [
        "config_kind = \"general\"",
        "schema_version = \"1\"",
        "`seed`, `ploidy` y\n`normalization_tolerance` opcionales",
        "`reference` obligatorio",
        "`reduced_references` opcional",
        "A, C, G y T exactamente una vez",
        "1e-9",
        "F2.2.2 valida sus sumas ALT multialélicas",
    ] {
        assert!(GENERAL_SCHEMA.contains(literal), "missing {literal}");
    }
    assert!(GENERAL_SCHEMA.contains("binario F1 no lo carga ni lo interpreta"));
    assert!(SCENARIO_SCHEMA.contains("config_kind = \"scenario\""));
    assert!(SCENARIO_SCHEMA.contains("No forma parte del contrato ejecutable de F1"));
    assert!(!SCENARIO_SCHEMA.contains("F1.2"));
    assert!(!SCENARIO_SCHEMA.contains("admitidas en F1"));
}
