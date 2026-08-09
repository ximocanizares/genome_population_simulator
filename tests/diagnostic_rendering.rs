use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::configuration::ExpectedConfigKind;

#[test]
fn human_and_json_are_rendered_from_the_same_structured_diagnostic() {
    let admitted = admit_toml_header(
        Path::new("fixtures/unknown.toml"),
        include_bytes!("fixtures/config/invalid-unknown-field.toml").to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("la cabecera es válida");
    let diagnostic = parse_declared_config(admitted).expect_err("el campo es desconocido");
    let human = diagnostic.render_human();
    let json = diagnostic.render_json_document();

    for literal in [
        diagnostic.code(),
        diagnostic.category(),
        diagnostic.message(),
    ] {
        assert!(human.contains(literal));
        assert!(json.contains(literal));
    }
    for (key, value) in diagnostic.context() {
        assert!(human.contains(key));
        assert!(human.contains(value));
        assert!(json.contains(key));
        assert!(json.contains(value));
    }
    assert_eq!(json.matches("\"diagnostic_kind\"").count(), 1);
    assert_eq!(json.matches("\"diagnostics\"").count(), 1);
}

#[test]
fn json_rendering_escapes_path_control_characters_without_emitting_jsonl() {
    let bytes = b"config_kind = \"general\"\nschema_version = \"2\"\n".to_vec();
    let diagnostic = admit_toml_header(
        Path::new("fixture-\"quoted\"\nline.toml"),
        bytes,
        ExpectedConfigKind::General,
    )
    .expect_err("versión no admitida");
    let json = diagnostic.render_json_document();

    assert!(json.contains("fixture-\\\"quoted\\\"\\nline.toml"));
    assert!(!json.contains('\n'));
    assert!(json.starts_with('{') && json.ends_with('}'));
}
