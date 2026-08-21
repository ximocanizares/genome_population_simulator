use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::configuration::ExpectedConfigKind;

fn parse(text: &str) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
    let admitted = admit_toml_header(
        Path::new("population-snvs.toml"),
        text.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )?;
    parse_declared_config(admitted).map(|_| ())
}

fn config(mode: &str, path_key: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "1"
[population_snvs]
mode = "{mode}"
{path_key} = ""
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
"#
    )
}

#[test]
fn supplied_and_reuse_paths_must_not_be_empty() {
    for (mode, field) in [("supplied", "input_path"), ("reuse", "metadata_path")] {
        let error = parse(&config(mode, field)).expect_err("empty path must fail");
        assert_eq!(error.code(), "GPS-CONF-008");
        assert_eq!(
            error.context()["fields"],
            format!("population_snvs.{field}")
        );
        assert!(error.context()["detail"].contains("vacía"));
        assert!(error.context().contains_key("location"));
    }
}
