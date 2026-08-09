use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::configuration::{DeclaredConfig, ExpectedConfigKind};

fn parse_general(
    bytes: &[u8],
) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let admitted = admit_toml_header(
        Path::new("general.toml"),
        bytes.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("header must be admitted");
    match parse_declared_config(admitted).expect("semantics must parse") {
        DeclaredConfig::General(config) => *config,
        DeclaredConfig::Scenario(_) => panic!("expected general configuration"),
    }
}

#[test]
fn reordered_toml_has_the_same_declared_semantics_but_keeps_its_own_bytes() {
    let first_bytes = include_bytes!("fixtures/config/general-minimal.toml");
    let second_bytes = include_bytes!("fixtures/config/general-minimal-reordered.toml");
    let first = parse_general(first_bytes);
    let second = parse_general(second_bytes);

    assert_eq!(first.route().value(), second.route().value());
    assert_eq!(first.reference().id.value(), second.reference().id.value());
    assert_eq!(
        first.reference().sequences[0].id.value(),
        second.reference().sequences[0].id.value()
    );
    assert_eq!(
        first.reference().sequences[0].class.value(),
        second.reference().sequences[0].class.value()
    );
    assert_eq!(
        first.reference().sequences[0].length.value(),
        second.reference().sequences[0].length.value()
    );
    assert_eq!(first.events(), second.events());
    assert_ne!(first.source().bytes(), second.source().bytes());
    assert_eq!(first.source().bytes(), first_bytes);
    assert_eq!(second.source().bytes(), second_bytes);
}

#[test]
fn lists_require_explicit_unique_ids() {
    for fixture in [
        include_bytes!("fixtures/config/invalid-event-missing-id.toml").as_slice(),
        include_bytes!("fixtures/config/invalid-event-duplicate-id.toml").as_slice(),
    ] {
        let admitted = admit_toml_header(
            Path::new("invalid.toml"),
            fixture.to_vec(),
            ExpectedConfigKind::General,
        )
        .expect("header remains valid");
        let error = parse_declared_config(admitted).expect_err("ID defect must be rejected");
        assert_eq!(error.code(), "GPS-CONF-008");
        assert_eq!(error.exit_code(), 2);
    }
}

#[test]
fn unknown_fields_are_rejected_with_location_after_header_admission() {
    let fixture = include_bytes!("fixtures/config/invalid-unknown-field.toml");
    let admitted = admit_toml_header(
        Path::new("unknown.toml"),
        fixture.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("header remains valid");
    let error = parse_declared_config(admitted).expect_err("unknown field must be rejected");

    assert_eq!(error.code(), "GPS-CONF-006");
    assert_eq!(
        error.context().get("field").map(String::as_str),
        Some("mystery")
    );
    assert!(error.context().contains_key("location"));
}

#[test]
fn scenario_lists_preserve_normative_order_and_explicit_identity() {
    let bytes = include_bytes!("fixtures/config/scenario-evolution.toml").to_vec();
    let admitted = admit_toml_header(
        Path::new("scenario.toml"),
        bytes,
        ExpectedConfigKind::Scenario,
    )
    .expect("valid scenario header");
    let DeclaredConfig::Scenario(scenario) =
        parse_declared_config(admitted).expect("valid declared scenario")
    else {
        panic!("expected scenario");
    };

    assert_eq!(scenario.id().value(), "evolution-small");
    assert_eq!(scenario.phases()[0].id().value(), "phase-1");
    assert_eq!(scenario.events()[0].id().value(), "sample-generation-10");
}
