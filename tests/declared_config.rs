use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::DeclaredPopulationSnvs;
use genome_population_simulator::configuration::ExpectedConfigKind;

fn parse(bytes: &[u8]) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let admitted = admit_toml_header(
        Path::new("general.toml"),
        bytes.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    parse_declared_config(admitted).expect("semantics")
}

#[test]
fn reordered_toml_has_the_same_declared_semantics_but_keeps_its_own_bytes() {
    let first = include_bytes!("fixtures/config/general-minimal.toml");
    let second = include_bytes!("fixtures/config/general-minimal-reordered.toml");
    let first_config = parse(first);
    let second_config = parse(second);
    assert_eq!(
        first_config.reference().id.value(),
        second_config.reference().id.value()
    );
    assert_eq!(
        first_config.reference().sequences[0].id.value(),
        second_config.reference().sequences[0].id.value(),
    );
    assert_eq!(
        first_config.reference().sequences[0].length.value(),
        second_config.reference().sequences[0].length.value(),
    );
    assert_ne!(
        first_config.source().bytes(),
        second_config.source().bytes()
    );
    assert_eq!(first_config.source().bytes(), first);
    assert_eq!(second_config.source().bytes(), second);
}

#[test]
fn unknown_fields_are_rejected_with_location() {
    let fixture = include_bytes!("fixtures/config/invalid-unknown-field.toml");
    let admitted = admit_toml_header(
        Path::new("unknown.toml"),
        fixture.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    let error = parse_declared_config(admitted).expect_err("unknown field");
    assert_eq!(error.code(), "GPS-CONF-006");
    assert!(error.context().contains_key("location"));
}

#[test]
fn population_snv_modes_are_exclusive_at_the_toml_boundary() {
    let supplied = parse(
        br#"config_kind = "general"
schema_version = "2"
[population_snvs]
mode = "supplied"
input_path = "input.tsv"
[reference]
kind = "simulated"
id = "ref-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1
"#,
    );
    assert!(matches!(
        supplied.population_snvs(),
        Some(DeclaredPopulationSnvs::Supplied { input_path }) if input_path.value() == "input.tsv"
    ));

    let invalid = br#"config_kind = "general"
schema_version = "2"
[population_snvs]
mode = "reuse"
input_path = "input.tsv"
metadata_path = "matrix.json"
[reference]
kind = "simulated"
id = "ref-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1
"#;
    let admitted = admit_toml_header(
        Path::new("invalid.toml"),
        invalid.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    assert!(parse_declared_config(admitted).is_err());
}

#[test]
fn dtc049_tolerance_defaults_and_is_limited_in_the_general_toml() {
    let base = r#"config_kind = "general"
schema_version = "2"
[reference]
kind = "simulated"
id = "ref-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1
"#;
    assert_eq!(
        build_effective_config(&parse(base.as_bytes()))
            .unwrap()
            .normalization()
            .tolerance()
            .to_bits(),
        1e-9_f64.to_bits()
    );
    let configured = base.replacen(
        "[reference]",
        "normalization_tolerance = 0.000001\n[reference]",
        1,
    );
    assert_eq!(
        build_effective_config(&parse(configured.as_bytes()))
            .unwrap()
            .normalization()
            .tolerance()
            .to_bits(),
        1e-6_f64.to_bits()
    );
    let invalid = base.replacen(
        "[reference]",
        "normalization_tolerance = 0.0000011\n[reference]",
        1,
    );
    assert!(build_effective_config(&parse(invalid.as_bytes())).is_err());
}
