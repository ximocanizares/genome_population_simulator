use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;

#[test]
fn changing_a_toml_parameter_changes_the_effective_value_without_code_changes() {
    let base = include_str!("fixtures/config/general-minimal.toml");
    let changed = base.replace("length = 1000", "length = 1001");
    let build = |source: &[u8]| {
        let admitted = admit_toml_header(
            Path::new("general.toml"),
            source.to_vec(),
            ExpectedConfigKind::General,
        )
        .expect("header");
        let declared = parse_declared_config(admitted).expect("declared");
        build_effective_config(&declared).expect("effective")
    };
    assert_eq!(build(base.as_bytes()).total_genome_size(), 1000);
    assert_eq!(build(changed.as_bytes()).total_genome_size(), 1001);
}
