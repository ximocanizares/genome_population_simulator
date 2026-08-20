use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
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
