use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;

fn parse(text: &str) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
    let admitted = admit_toml_header(
        Path::new("population-snvs.toml"),
        text.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )?;
    parse_declared_config(admitted).map(|_| ())
}

fn build(text: &str) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
    let admitted = admit_toml_header(
        Path::new("population-snvs.toml"),
        text.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )?;
    let declared = parse_declared_config(admitted)?;
    build_effective_config(&declared).map(|_| ())
}

fn generated_config(parameters: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "1"
[population_snvs]
mode = "generate"
[population_snvs.generate]
{parameters}
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
length = 10
"#
    )
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

#[test]
fn generate_requires_the_closed_f2_2_5_contract() {
    let text = r#"config_kind = "general"
schema_version = "1"
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 2
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9
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
length = 2
"#;
    let admitted = admit_toml_header(
        Path::new("population-snvs.toml"),
        text.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("admitted");
    let effective = build_effective_config(&parse_declared_config(admitted).expect("declared"))
        .expect("effective");
    assert!(effective.population_snvs_generation().is_some());
}

#[test]
fn generate_rejects_incompatible_count_and_frequency_contracts() {
    let invalid = config("generate", "generate");
    assert!(parse(&invalid).is_err());
    let invalid_frequency = r#"config_kind = "general"
schema_version = "1"
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P"
frequency_model = "log"
min_alt_af = 0.1
max_alt_af = 0.9
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
    assert!(parse(invalid_frequency).is_err());
}

#[test]
fn effective_config_rejects_generation_domain_errors_before_generation() {
    let cases = [
        (
            "intervalo AF invertido",
            r#"count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.8
max_alt_af = 0.1"#,
            "InvalidFrequencyInterval",
        ),
        (
            "intervalo AF igualado",
            r#"count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.8
max_alt_af = 0.8"#,
            "InvalidFrequencyInterval",
        ),
        (
            "región con inicio posterior al fin",
            r#"count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
substitution_profile = "uniform_snv12_v1"
outside_weight = 1.0
[[population_snvs.generate.regions]]
contig = "chr1"
start = 4
end = 3
weight = 1.0
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9"#,
            "InvalidRegion",
        ),
        (
            "regiones solapadas",
            r#"count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
substitution_profile = "uniform_snv12_v1"
outside_weight = 1.0
[[population_snvs.generate.regions]]
contig = "chr1"
start = 1
end = 4
weight = 1.0
[[population_snvs.generate.regions]]
contig = "chr1"
start = 4
end = 5
weight = 1.0
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9"#,
            "OverlappingRegions",
        ),
        (
            "población duplicada",
            r#"count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8"#,
            "DuplicatePopulation",
        ),
        (
            "peso negativo",
            r#"count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
substitution_profile = "uniform_snv12_v1"
outside_weight = -1.0
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9"#,
            "InvalidWeight",
        ),
        (
            "peso no finito",
            r#"count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
substitution_profile = "uniform_snv12_v1"
outside_weight = nan
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9"#,
            "InvalidWeight",
        ),
        (
            "media negativa",
            r#"count_mode = "poisson"
mean_loci = -1.0
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9"#,
            "InvalidCount",
        ),
        (
            "media no finita",
            r#"count_mode = "poisson"
mean_loci = nan
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9"#,
            "InvalidCount",
        ),
    ];

    for (name, parameters, expected) in cases {
        let error = build(&generated_config(parameters)).expect_err(name);
        assert_eq!(error.code(), "GPS-CONF-008", "{name}");
        assert_eq!(
            error.context()["fields"],
            "population_snvs.generate",
            "{name}"
        );
        assert!(
            error.context()["detail"].contains(expected),
            "{name}: {error:?}"
        );
    }
}
