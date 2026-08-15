use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::{
    DeclaredConfig, DeclaredReductionTarget, ExpectedConfigKind,
};

const BASE: &str = r#"config_kind = "general"
schema_version = "1"
route = "direct"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[reference]
kind = "supplied"
id = "reference-main"
fasta = "source.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 10
ploidy = 2

[[reference.sequences]]
id = "scaffold-1"
class = "scaffold"
length = 7
ploidy = 2
"#;

fn declared(
    reduction: &str,
) -> Box<genome_population_simulator::configuration::DeclaredGeneralConfig> {
    let bytes = format!("{BASE}\n{reduction}\n").into_bytes();
    let admitted = admit_toml_header(
        Path::new("config/general.toml"),
        bytes,
        ExpectedConfigKind::General,
    )
    .unwrap();
    let DeclaredConfig::General(config) = parse_declared_config(admitted).unwrap() else {
        panic!("expected general config")
    };
    config
}

#[test]
fn declared_reduction_preserves_target_variant_literal_and_location() {
    let config = declared(
        r#"[[reduced_references]]
id = "reference-small"
source_reference_id = "reference-main"
region_size_bases = 4
reduction_percent = 2.5e1"#,
    );
    let reduction = &config.reduced_references()[0];
    assert_eq!(reduction.id.value(), "reference-small");
    assert_eq!(*reduction.region_size_bases.value(), 4);
    let DeclaredReductionTarget::ReductionPercent(percent) = &reduction.target else {
        panic!("expected percentage target")
    };
    assert_eq!(percent.literal(), "2.5e1");
    assert!(percent.location().line > 1);
}

#[test]
fn exactly_one_target_and_known_fields_are_required() {
    for fields in [
        "target_bases = 8\nreduction_percent = 25",
        "",
        "target_bases = 8\nunknown = true",
    ] {
        let text = format!(
            "[[reduced_references]]\nid = \"reference-small\"\nsource_reference_id = \"reference-main\"\nregion_size_bases = 4\n{fields}"
        );
        let bytes = format!("{BASE}\n{text}\n").into_bytes();
        let admitted = admit_toml_header(
            Path::new("config/general.toml"),
            bytes,
            ExpectedConfigKind::General,
        )
        .unwrap();
        let error = parse_declared_config(admitted).expect_err("invalid reduction must fail");
        assert!(matches!(error.code(), "GPS-CONF-006" | "GPS-CONF-008"));
        assert_eq!(error.exit_code(), 2);
    }
}

#[test]
fn effective_reduction_uses_exact_decimal_and_own_subseed() {
    let config = declared(
        r#"[[reduced_references]]
id = "reference-small"
source_reference_id = "reference-main"
region_size_bases = 4
reduction_percent = 33.333333333333333333"#,
    );
    let effective = build_effective_config(&config, &[], &[]).unwrap();
    let reduction = &effective.reduced_references()[0];
    assert_eq!(reduction.target_bases(), 11);
    assert_eq!(reduction.target().declared_value(), "33.333333333333333333");
    assert_eq!(
        reduction.selection_subseed().to_hex(),
        "7fa724d52b2db97a043920dbc3062a34c6377eab76f1d9d6e4ff1d513747decd"
    );
}

#[test]
fn zero_endpoints_source_mismatch_and_non_reduction_are_rejected() {
    for reduction in [
        "[[reduced_references]]\nid = \"small\"\nsource_reference_id = \"reference-main\"\nregion_size_bases = 0\ntarget_bases = 1",
        "[[reduced_references]]\nid = \"small\"\nsource_reference_id = \"other\"\nregion_size_bases = 4\ntarget_bases = 1",
        "[[reduced_references]]\nid = \"small\"\nsource_reference_id = \"reference-main\"\nregion_size_bases = 4\nreduction_percent = 0",
        "[[reduced_references]]\nid = \"small\"\nsource_reference_id = \"reference-main\"\nregion_size_bases = 4\nreduction_percent = 100",
        "[[reduced_references]]\nid = \"small\"\nsource_reference_id = \"reference-main\"\nregion_size_bases = 4\nreduction_percent = nan",
        "[[reduced_references]]\nid = \"small\"\nsource_reference_id = \"reference-main\"\nregion_size_bases = 4\nreduction_percent = -1",
        "[[reduced_references]]\nid = \"small\"\nsource_reference_id = \"reference-main\"\nregion_size_bases = 4\ntarget_bases = 17",
    ] {
        let config = declared(reduction);
        let error = build_effective_config(&config, &[], &[]).expect_err("preflight must reject");
        assert_eq!(error.exit_code(), 2);
        if matches!(error.code(), "GPS-CONF-019" | "GPS-CONF-020") {
            assert_eq!(
                error.context().get("path").map(String::as_str),
                Some("config/general.toml")
            );
            assert!(error.context().contains_key("location"));
        }
    }
}
