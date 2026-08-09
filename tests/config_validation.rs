use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::validate_general_config;
use genome_population_simulator::configuration::{
    DeclaredConfig, DeclaredGeneralConfig, DeclaredReferenceSource, DeclaredScenarioConfig,
    ExpectedConfigKind, SequenceClass,
};

fn general(bytes: &[u8]) -> DeclaredGeneralConfig {
    let admitted = admit_toml_header(
        Path::new("general.toml"),
        bytes.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("valid header");
    let DeclaredConfig::General(config) = parse_declared_config(admitted).expect("declared config")
    else {
        panic!("expected general")
    };
    *config
}

fn scenario(bytes: &[u8]) -> DeclaredScenarioConfig {
    let admitted = admit_toml_header(
        Path::new("scenario.toml"),
        bytes.to_vec(),
        ExpectedConfigKind::Scenario,
    )
    .expect("valid header");
    let DeclaredConfig::Scenario(config) =
        parse_declared_config(admitted).expect("declared config")
    else {
        panic!("expected scenario")
    };
    *config
}

#[test]
fn known_sum_classes_and_zero_explicit_events_are_preserved() {
    let config = general(include_bytes!("fixtures/config/general-minimal.toml"));
    let summary = validate_general_config(&config, &[]).expect("valid direct config");

    assert_eq!(summary.total_genome_size, 1000);
    assert_eq!(
        summary.sequence_classes.get("chr1"),
        Some(&SequenceClass::Chromosome)
    );
    assert!(summary.explicit_event_ids.is_empty());
}

#[test]
fn simulated_reference_preserves_sequence_count_lengths_classes_and_composition() {
    let source = br#"config_kind = "general"
schema_version = "1"
route = "direct"
[reference]
kind = "simulated"
id = "ref-two"
default_ploidy = 2
composition = [
  { id = "A", category = "base", value = 0.6 },
  { id = "C", category = "base", value = 0.4 },
]
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 7
[[reference.sequences]]
id = "scaffold-1"
class = "scaffold"
length = 5
"#;
    let config = general(source);
    assert_eq!(config.reference().sequences.len(), 2);
    assert_eq!(*config.reference().sequences[0].length.value(), 7);
    assert_eq!(*config.reference().sequences[1].length.value(), 5);
    assert_eq!(
        *config.reference().sequences[0].class.value(),
        SequenceClass::Chromosome
    );
    assert_eq!(
        *config.reference().sequences[1].class.value(),
        SequenceClass::Scaffold
    );
    let DeclaredReferenceSource::Simulated { composition, .. } = &config.reference().source else {
        panic!("la referencia debe ser simulada")
    };
    assert_eq!(composition.len(), 2);
    assert_eq!(composition[0].id.value(), "A");
    assert_eq!(composition[0].value.literal(), "0.6");
    assert_eq!(composition[1].id.value(), "C");
    assert_eq!(composition[1].value.literal(), "0.4");

    let summary = validate_general_config(&config, &[]).expect("referencia simulada válida");
    assert_eq!(summary.total_genome_size, 12);
}

#[test]
fn total_size_overflow_is_rejected_without_wrapping() {
    let config = general(include_bytes!("fixtures/config/invalid-size-overflow.toml"));
    let error = validate_general_config(&config, &[]).expect_err("sum must overflow");
    assert_eq!(error.code(), "GPS-CONF-008");
    assert!(error.context()["detail"].contains("desborda"));
}

#[test]
fn supplied_sequence_without_effective_ploidy_is_rejected() {
    let config = general(include_bytes!(
        "fixtures/config/invalid-supplied-ploidy-missing.toml"
    ));
    let error = validate_general_config(&config, &[]).expect_err("ploidy is mandatory");
    assert_eq!(error.code(), "GPS-CONF-014");
}

#[test]
fn incompatible_frequency_forms_fail_during_declared_admission() {
    let admitted = admit_toml_header(
        Path::new("frequency.toml"),
        include_bytes!("fixtures/config/invalid-frequency-forms.toml").to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("valid header");
    let error = parse_declared_config(admitted).expect_err("mixed forms must fail");
    assert_eq!(error.code(), "GPS-CONF-008");
}

#[test]
fn route_and_scenario_must_agree() {
    let evolution_general = b"config_kind = \"general\"\nschema_version = \"1\"\nroute = \"evolution\"\n[[scenarios]]\nid = \"families-small\"\npath = \"scenario.toml\"\nroute = \"evolution\"\n[reference]\nkind = \"simulated\"\nid = \"ref\"\ndefault_ploidy = 2\ncomposition = [{ id = \"A\", category = \"base\", value = 1.0 }]\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 10\n";
    let general = general(evolution_general);
    let families = scenario(include_bytes!("fixtures/config/scenario-families.toml"));
    let error = validate_general_config(&general, &[&families])
        .expect_err("scenario route mismatch must fail");
    assert_eq!(error.code(), "GPS-CONF-008");
}

#[test]
fn validation_does_not_create_a_final_directory() {
    let config = general(include_bytes!("fixtures/config/general-minimal.toml"));
    let destination =
        std::env::temp_dir().join(format!("gps-f1.2-no-publication-{}", std::process::id()));
    assert!(!destination.exists());
    validate_general_config(&config, &[]).expect("valid config");
    assert!(!destination.exists());
}
