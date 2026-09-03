use std::path::Path;

use genome_population_simulator::adapters::{
    admit_toml_header, parse_declared_f5_general, parse_declared_f5_scenario,
    validate_declared_f5_pair,
};
use genome_population_simulator::configuration::ExpectedConfigKind;
use sha2::Digest;

const HASH: &str = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

fn general(seed: &str) -> String {
    format!(
        "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{seed}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{HASH}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n"
    )
}

fn scenario(extra: &str) -> String {
    format!(
        "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = 2\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{HASH}\"\nblock_count = 1\n{extra}\n[[subpopulations]]\npopulation_id = \"P1\"\nindividual_count = 2\nfounder_ids = [\"I1\", \"I2\"]\n"
    )
}

fn parse_pair(general_source: &str, scenario_source: &str) -> Result<(), String> {
    let general = parse_declared_f5_general(
        admit_toml_header(
            Path::new("authoring/general.toml"),
            general_source.as_bytes().to_vec(),
            ExpectedConfigKind::F5General,
        )
        .map_err(|error| error.render_human())?,
    )
    .map_err(|error| error.render_human())?;
    let scenario = parse_declared_f5_scenario(
        admit_toml_header(
            Path::new("authoring/scenarios/scenario.toml"),
            scenario_source.as_bytes().to_vec(),
            ExpectedConfigKind::F5Scenario,
        )
        .map_err(|error| error.render_human())?,
    )
    .map_err(|error| error.render_human())?;
    validate_declared_f5_pair(&general, &scenario).map_err(|error| error.render_human())
}

fn parse_scenario(scenario_source: &str) -> Result<(), String> {
    parse_declared_f5_scenario(
        admit_toml_header(
            Path::new("authoring/scenarios/scenario.toml"),
            scenario_source.as_bytes().to_vec(),
            ExpectedConfigKind::F5Scenario,
        )
        .map_err(|error| error.render_human())?,
    )
    .map(|_| ())
    .map_err(|error| error.render_human())
}

#[test]
fn admits_minimal_closed_f5_pair() {
    assert!(parse_pair(&general(&"a".repeat(64)), &scenario("")).is_ok());
}

#[test]
fn rejects_legacy_general_and_noncanonical_paths() {
    let legacy =
        general(&"a".repeat(64)).replace("schema_version = \"3\"", "schema_version = \"2\"");
    assert!(
        admit_toml_header(
            Path::new("general.toml"),
            legacy.into_bytes(),
            ExpectedConfigKind::F5General
        )
        .is_err()
    );
    let escaped = general(&"a".repeat(64)).replace("scenarios/scenario.toml", "../scenario.toml");
    assert!(parse_pair(&escaped, &scenario("")).is_err());
}

#[test]
fn all_f5_links_reject_noncanonical_path_spellings_in_toml() {
    let invalid_paths = [
        "scenarios//scenario.toml",
        "scenarios/./scenario.toml",
        "scenarios/.",
        "scenarios/",
        "scenarios/../scenario.toml",
        "/scenarios/scenario.toml",
        "scenarios\\scenario.toml",
        "scenarios/scenario file.toml",
    ];

    for invalid in invalid_paths {
        let invalid_general = general(&"a".repeat(64)).replace("scenarios/scenario.toml", invalid);
        assert!(
            parse_pair(&invalid_general, &scenario("")).is_err(),
            "general→scenario must reject {invalid:?}"
        );

        let invalid_reference = general(&"a".repeat(64)).replace("references/ref.fa", invalid);
        assert!(
            parse_pair(&invalid_reference, &scenario("")).is_err(),
            "general→reference must reject {invalid:?}"
        );

        let invalid_profile = scenario("").replace("profiles/profile.metadata.json", invalid);
        assert!(
            parse_pair(&general(&"a".repeat(64)), &invalid_profile).is_err(),
            "scenario→profile must reject {invalid:?}"
        );

        let invalid_mask = scenario(&format!(
            "excluded_regions_path = \"{invalid}\"\nexcluded_regions_sha256 = \"{HASH}\"\n"
        ));
        assert!(
            parse_pair(&general(&"a".repeat(64)), &invalid_mask).is_err(),
            "scenario→optional mask must reject {invalid:?}"
        );
    }
}

#[test]
fn scenario_u64_fields_accept_max_and_reject_max_plus_one_in_toml() {
    const MAX: &str = "18446744073709551615";
    const MAX_PLUS_ONE: &str = "18446744073709551616";

    for (field, original) in [
        ("founder_count", "founder_count = 2"),
        ("covered_bases", "block_count = 1"),
        ("individual_count", "individual_count = 2"),
    ] {
        let max_assignment = format!("{field} = {MAX}");
        let admitted = scenario("").replace(original, &max_assignment);
        assert!(
            parse_scenario(&admitted).is_ok(),
            "{field} must accept u64::MAX"
        );

        let overflow_assignment = format!("{field} = {MAX_PLUS_ONE}");
        let rejected = scenario("").replace(original, &overflow_assignment);
        assert!(
            parse_scenario(&rejected).is_err(),
            "{field} must reject u64::MAX + 1"
        );
    }
}

#[test]
fn rejects_target_mask_founder_and_f_inconsistencies() {
    assert!(parse_pair(&general(&"a".repeat(64)), &scenario("covered_bases = 2\n")).is_err());
    assert!(
        parse_pair(
            &general(&"a".repeat(64)),
            &scenario("excluded_regions_path = \"mask.bed\"\n")
        )
        .is_err()
    );
    assert!(parse_pair(&general(&"a".repeat(64)), &scenario("founder_count = 3\n")).is_err());
    let haploid = general(&"a".repeat(64)).replace("ploidy = 2", "ploidy = 1");
    assert!(parse_pair(&haploid, &scenario("inbreeding_f = 0.2\n")).is_err());
}

#[test]
fn coverage_fraction_uses_the_exact_decimal_literal_and_override_is_u32() {
    for literal in ["0.1", "0.000000000000000001", "1.0", "1.000000000000000000"] {
        let candidate = scenario(&format!("block_count = 1\ncoverage_fraction = {literal}\n"))
            .replace("block_count = 1\nblock_count = 1\n", "");
        assert!(
            parse_pair(&general(&"a".repeat(64)), &candidate).is_ok(),
            "must admit {literal}"
        );
    }

    for literal in ["0.0", "1.000000000000000001", "1.1", "+0.5", "-0.5", "1e-1"] {
        let candidate = scenario(&format!("block_count = 1\ncoverage_fraction = {literal}\n"))
            .replace("block_count = 1\nblock_count = 1\n", "");
        assert!(
            parse_pair(&general(&"a".repeat(64)), &candidate).is_err(),
            "must reject {literal}"
        );
    }

    let oversized = scenario("").replace(
        "founder_ids = [\"I1\", \"I2\"]\n",
        "founder_ids = [\"I1\", \"I2\"]\n\n[[sampling_overrides]]\npopulation_id = \"P1\"\nblock_index = 4294967296\nmode = \"probabilistic\"\n",
    );
    assert!(parse_pair(&general(&"a".repeat(64)), &oversized).is_err());
}

#[test]
fn rejects_unknown_fields_duplicate_founders_and_noncanonical_seed() {
    assert!(parse_pair(&general(&"A".repeat(64)), &scenario("")).is_err());
    assert!(parse_pair(&general(&"a".repeat(64)), &scenario("mystery = true\n")).is_err());
    let duplicated = scenario("").replace("[\"I1\", \"I2\"]", "[\"I1\", \"I1\"]");
    assert!(parse_pair(&general(&"a".repeat(64)), &duplicated).is_err());
}

#[test]
fn legacy_schemas_remain_byte_for_byte_historical() {
    let schemas = [
        (
            &include_bytes!("../schemas/config-general-v2.schema.json")[..],
            "76ec2d384d8260ed728edfa0627ae5b9cf02fb0914594f0db14277ceec506e63",
        ),
        (
            &include_bytes!("../schemas/effective-config-v4.schema.json")[..],
            "a1ec75343e88098e196af275f0346b1617bf8da2902346af924f3b0be438c700",
        ),
        (
            &include_bytes!("../schemas/manifest-general-v4.schema.json")[..],
            "ada06880f5766c835a8f29609e3af1a969b40ff4423e0723f4829efb148f7b41",
        ),
        (
            &include_bytes!("../schemas/truth-v1.schema.json")[..],
            "e3d49abd11447370d34e4acf29b755a288fb310641474218086ed60099233fd7",
        ),
        (
            &include_bytes!("../catalogs/diagnostics-v3.tsv")[..],
            "6e8a8fa3c8c1306ed360bad6de402649020d109a10fc0c5b4072d47df66f64f9",
        ),
    ];
    for (bytes, expected) in schemas {
        assert_eq!(format!("{:x}", sha2::Sha256::digest(bytes)), expected);
    }
}
