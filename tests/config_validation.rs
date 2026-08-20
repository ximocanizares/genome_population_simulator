use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{build_effective_config, validate_general_config};
use genome_population_simulator::configuration::{
    DeclaredReferenceSource, ExpectedConfigKind, SequenceClass,
};

fn general(bytes: &[u8]) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let admitted = admit_toml_header(
        Path::new("general.toml"),
        bytes.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("valid header");
    parse_declared_config(admitted).expect("declared config")
}

#[test]
fn simulated_reference_preserves_sequence_count_lengths_classes_and_composition() {
    let config = general(include_bytes!("fixtures/config/general-minimal.toml"));
    let summary = validate_general_config(&config).expect("valid config");
    assert_eq!(summary.total_genome_size, 1000);
    assert_eq!(config.reference().sequences.len(), 1);
    assert_eq!(
        *config.reference().sequences[0].class.value(),
        SequenceClass::Chromosome
    );
    assert_eq!(*config.reference().sequences[0].length.value(), 1000);
    let DeclaredReferenceSource::Simulated { composition } = &config.reference().source else {
        panic!("expected simulated reference")
    };
    assert_eq!(
        composition
            .iter()
            .map(|value| value.base.value().as_str())
            .collect::<Vec<_>>(),
        ["A", "C", "G", "T"]
    );
}

#[test]
fn multiple_sequences_preserve_order_classes_lengths_and_total_size() {
    let source = br#"config_kind = "general"
schema_version = "1"
[reference]
kind = "simulated"
id = "ref-multiple"
[reference.composition]
A = 0.1
C = 0.2
G = 0.3
T = 0.4
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
    let summary = validate_general_config(&config).expect("valid multi-sequence config");

    assert_eq!(summary.total_genome_size, 12);
    assert_eq!(config.reference().sequences.len(), 2);
    assert_eq!(config.reference().sequences[0].id.value(), "chr1");
    assert_eq!(
        *config.reference().sequences[0].class.value(),
        SequenceClass::Chromosome
    );
    assert_eq!(*config.reference().sequences[0].length.value(), 7);
    assert_eq!(config.reference().sequences[1].id.value(), "scaffold-1");
    assert_eq!(
        *config.reference().sequences[1].class.value(),
        SequenceClass::Scaffold
    );
    assert_eq!(*config.reference().sequences[1].length.value(), 5);
    let DeclaredReferenceSource::Simulated { composition } = &config.reference().source else {
        panic!("expected simulated reference")
    };
    assert_eq!(
        composition
            .iter()
            .map(|value| value.value.literal())
            .collect::<Vec<_>>(),
        ["0.1", "0.2", "0.3", "0.4"]
    );
}

#[test]
fn reduced_reference_requires_a_supplied_source_reference() {
    let source = br#"config_kind = "general"
schema_version = "1"
seed = "0000000000000000000000000000000000000000000000000000000000000000"
[reference]
kind = "simulated"
id = "ref-simulated"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 100
[[reduced_references]]
id = "ref-reduced"
source_reference_id = "ref-simulated"
region_size_bases = 10
target_bases = 50
"#;

    let error = validate_general_config(&general(source))
        .expect_err("a simulated reference cannot be reduced in F1");
    assert_eq!(error.code(), "GPS-CONF-008");
    assert_eq!(error.context()["fields"], "reduced_references");
    assert!(error.context()["detail"].contains("suministrada"));
}

#[test]
fn composition_uses_fixed_tolerance_and_effective_values_are_normalized() {
    let source = br#"config_kind = "general"
schema_version = "1"
[reference]
kind = "simulated"
id = "ref"
[reference.composition]
A = 0.2500000001
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 10
"#;
    let config = general(source);
    let effective = build_effective_config(&config).expect("within fixed tolerance");
    let sum: f64 = effective
        .composition()
        .expect("composition")
        .values()
        .iter()
        .sum();
    assert!((sum - 1.0).abs() <= 1e-9);
}

#[test]
fn total_size_overflow_is_rejected_without_wrapping() {
    let source = br#"config_kind = "general"
schema_version = "1"
[reference]
kind = "supplied"
id = "ref"
fasta = "source.fa"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 18446744073709551615
[[reference.sequences]]
id = "chr2"
class = "chromosome"
length = 1
"#;
    let error = validate_general_config(&general(source)).expect_err("sum must overflow");
    assert_eq!(error.code(), "GPS-CONF-008");
    assert!(error.context()["detail"].contains("desborda"));
}

#[test]
fn future_root_and_nested_fields_are_rejected() {
    let base = r#"config_kind = "general"
schema_version = "1"
[reference]
kind = "supplied"
id = "ref"
fasta = "source.fa"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1
"#;
    for future_field in [
        "route = \"direct\"",
        "default_ploidy = 2",
        "scenarios = []",
        "events = []",
        "individuals = []",
        "populations = []",
        "phases = []",
        "groups = []",
        "blocks = []",
        "rules = []",
        "frequency_sets = []",
        "normalization_tolerance = 1e-9",
    ] {
        let source = base.replacen("[reference]", &format!("{future_field}\n[reference]"), 1);
        let admitted = admit_toml_header(
            Path::new("future.toml"),
            source.into_bytes(),
            ExpectedConfigKind::General,
        )
        .unwrap();
        let error = parse_declared_config(admitted).expect_err("future field must fail");
        assert_eq!(error.code(), "GPS-CONF-006");
    }

    let nested = br#"config_kind = "general"
schema_version = "1"
[reference]
kind = "supplied"
id = "ref"
fasta = "source.fa"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1
ploidy = 2
"#;
    let admitted = admit_toml_header(
        Path::new("future.toml"),
        nested.to_vec(),
        ExpectedConfigKind::General,
    )
    .unwrap();
    let error = parse_declared_config(admitted).expect_err("ploidy must fail in F1");
    assert_eq!(error.code(), "GPS-CONF-006");
}

#[test]
fn reference_source_fields_are_mutually_exclusive() {
    let simulated_with_fasta = br#"config_kind = "general"
schema_version = "1"
[reference]
kind = "simulated"
id = "ref"
fasta = "source.fa"
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
    let supplied_with_composition = br#"config_kind = "general"
schema_version = "1"
[reference]
kind = "supplied"
id = "ref"
fasta = "source.fa"
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
    for source in [
        simulated_with_fasta.as_slice(),
        supplied_with_composition.as_slice(),
    ] {
        let admitted = admit_toml_header(
            Path::new("exclusive.toml"),
            source.to_vec(),
            ExpectedConfigKind::General,
        )
        .unwrap();
        let error = parse_declared_config(admitted).expect_err("exclusive fields must fail");
        assert_eq!(error.code(), "GPS-CONF-008");
    }
}
