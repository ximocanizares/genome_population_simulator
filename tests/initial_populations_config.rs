use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{
    EffectiveInitialGenotypeModel, build_effective_config,
};
use genome_population_simulator::configuration::ExpectedConfigKind;
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);

fn build(
    source: &str,
) -> Result<
    genome_population_simulator::application::EffectiveConfig,
    genome_population_simulator::diagnostics::Diagnostic,
> {
    let admitted = admit_toml_header(
        Path::new("initial-populations.toml"),
        source.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )?;
    build_effective_config(&parse_declared_config(admitted)?)
}

fn config(populations: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
{populations}
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

fn cohort(population_id: &str, ploidy: u8, model: &str, extra: &str) -> String {
    let identities = if extra.contains("individual_ids") || extra.contains("founder_keys") {
        extra.to_owned()
    } else {
        format!(
            "founder_keys = [\"{population_id}-founder-1\", \"{population_id}-founder-2\"]\n{extra}"
        )
    };
    format!(
        r#"[[initial_populations]]
population_id = "{population_id}"
individual_count = 2
ploidy = {ploidy}
genotype_model = "{model}"
{identities}
"#
    )
}

#[test]
fn hwe_accepts_every_effective_ploidy_and_derives_stable_unique_ids() {
    let populations = (1..=10)
        .map(|ploidy| cohort(&format!("P{ploidy}"), ploidy, "hwe", ""))
        .collect::<String>();
    let effective = build(&config(&populations)).expect("all F2 ploidies are valid under HWE");
    assert_eq!(effective.initial_populations().len(), 10);

    let ids = effective
        .initial_populations()
        .iter()
        .flat_map(|population| population.individual_ids())
        .map(|id| id.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(ids.len(), 20);
    assert!(
        effective
            .initial_populations()
            .iter()
            .enumerate()
            .all(
                |(index, population)| population.ploidy().copy_count() == (index + 1) as u8
                    && matches!(
                        population.genotype_model(),
                        EffectiveInitialGenotypeModel::HardyWeinberg
                    )
            )
    );

    let reordered = (1..=10)
        .rev()
        .map(|ploidy| cohort(&format!("P{ploidy}"), ploidy, "hwe", ""))
        .collect::<String>();
    let reordered = build(&config(&reordered)).expect("reordered cohorts");
    for population in effective.initial_populations() {
        let counterpart = reordered
            .initial_populations()
            .iter()
            .find(|other| other.population_id() == population.population_id())
            .expect("same population");
        assert_eq!(population.individual_ids(), counterpart.individual_ids());
    }
}

#[test]
fn provided_ids_replace_derivation_only_with_exact_cardinality_and_global_uniqueness() {
    let supplied = cohort(
        "P1",
        2,
        "hwe",
        "individual_ids = [\"sample-A\", \"sample-B\"]",
    );
    let effective = build(&config(&supplied)).expect("provided IDs");
    assert_eq!(
        effective.initial_populations()[0]
            .individual_ids()
            .iter()
            .map(|id| id.as_str())
            .collect::<Vec<_>>(),
        ["sample-A", "sample-B"]
    );

    for (extra, expected_code, expected_field) in [
        (
            "individual_ids = [\"only-one\"]",
            "GPS-CONF-008",
            "individual_ids",
        ),
        (
            "individual_ids = [\"same\", \"same\"]",
            "GPS-DOMAIN-002",
            "",
        ),
    ] {
        let error = build(&config(&cohort("P1", 2, "hwe", extra))).expect_err("invalid IDs");
        assert_eq!(error.code(), expected_code);
        if !expected_field.is_empty() {
            assert!(error.context()["fields"].contains(expected_field));
        }
    }

    let duplicated_across_cohorts = format!(
        "{}{}",
        cohort("P1", 2, "hwe", "individual_ids = [\"shared\", \"one\"]"),
        cohort("P2", 2, "hwe", "individual_ids = [\"shared\", \"two\"]"),
    );
    let error = build(&config(&duplicated_across_cohorts)).expect_err("global duplicate");
    assert_eq!(error.code(), "GPS-DOMAIN-002");
}

#[test]
fn omitted_ids_require_unique_founder_keys_and_derive_from_them() {
    let source = cohort(
        "P1",
        2,
        "hwe",
        "founder_keys = [\"founder-alpha\", \"founder-beta\"]",
    );
    let first = build(&config(&source)).expect("founder keys derive IDs");
    let reordered = build(&config(&cohort(
        "P1",
        2,
        "hwe",
        "founder_keys = [\"founder-beta\", \"founder-alpha\"]",
    )))
    .expect("reordered founder keys");
    let first_by_key = first.initial_populations()[0]
        .founder_keys()
        .expect("derived keys")
        .iter()
        .zip(first.initial_populations()[0].individual_ids())
        .map(|(key, id)| (key.as_str(), id.as_str()))
        .collect::<std::collections::BTreeMap<_, _>>();
    let reordered_by_key = reordered.initial_populations()[0]
        .founder_keys()
        .expect("derived keys")
        .iter()
        .zip(reordered.initial_populations()[0].individual_ids())
        .map(|(key, id)| (key.as_str(), id.as_str()))
        .collect::<std::collections::BTreeMap<_, _>>();
    assert_eq!(first_by_key, reordered_by_key);

    let missing_founder_keys =
        cohort("P1", 2, "hwe", "founder_keys = []").replace("founder_keys = []\n", "");
    let missing = build(&config(&missing_founder_keys))
        .expect_err("founder keys are required when IDs are omitted");
    assert_eq!(missing.code(), "GPS-CONF-008");
    assert!(missing.context()["fields"].contains("founder_keys"));

    for (extra, field) in [
        ("founder_keys = []", "founder_keys"),
        ("founder_keys = [\"only-one\"]", "founder_keys"),
        ("founder_keys = [\"same\", \"same\"]", "founder_keys"),
        (
            "individual_ids = [\"I1\", \"I2\"]\nfounder_keys = [\"a\", \"b\"]",
            "individual_ids,initial_populations.founder_keys",
        ),
    ] {
        let error = build(&config(&cohort("P1", 2, "hwe", extra)))
            .expect_err("invalid founder-key contract");
        assert_eq!(error.code(), "GPS-CONF-008");
        assert!(error.context()["fields"].contains(field));
    }

    let repeated_across_cohorts = format!(
        "{}{}",
        cohort("P1", 2, "hwe", "founder_keys = [\"shared\", \"P1-only\"]"),
        cohort("P2", 2, "hwe", "founder_keys = [\"shared\", \"P2-only\"]"),
    );
    let repeated = build(&config(&repeated_across_cohorts))
        .expect("the same founder key is valid in distinct populations");
    assert_ne!(
        repeated.initial_populations()[0].individual_ids()[0],
        repeated.initial_populations()[1].individual_ids()[0],
        "the population remains part of the deterministic derivation scope"
    );
}

#[test]
fn f_diploid_is_exclusive_to_diploids_and_validates_its_coefficient() {
    for coefficient in ["0", "1", "0.25"] {
        let population = cohort(
            "P2",
            2,
            "f_diploid",
            &format!("inbreeding_coefficient = {coefficient}"),
        );
        let effective = build(&config(&population)).expect("valid F coefficient");
        assert!(matches!(
            effective.initial_populations()[0].genotype_model(),
            EffectiveInitialGenotypeModel::FDiploid { .. }
        ));
    }
    for ploidy in [1, 3, 10] {
        let population = cohort("P", ploidy, "f_diploid", "inbreeding_coefficient = 0.2");
        let error = build(&config(&population)).expect_err("F outside diploidy");
        assert_eq!(error.code(), "GPS-CONF-008");
        assert!(error.context()["fields"].contains("ploidy"));
    }
    for extra in [
        "",
        "inbreeding_coefficient = -0.1",
        "inbreeding_coefficient = 1.1",
    ] {
        let population = cohort("P", 2, "f_diploid", extra);
        assert!(build(&config(&population)).is_err());
    }
    let hwe_with_f = cohort("P", 2, "hwe", "inbreeding_coefficient = 0.2");
    assert!(build(&config(&hwe_with_f)).is_err());
}

#[test]
fn authoring_schema_and_rust_accept_the_same_initial_cohort_shape() {
    let schema: serde_json::Value =
        serde_json::from_str(include_str!("../schemas/config-general-v2.schema.json"))
            .expect("authoring schema JSON");
    let cohort_schema = &schema["$defs"]["initial_population"];
    assert_eq!(
        cohort_schema["required"]
            .as_array()
            .expect("required")
            .len(),
        4
    );
    assert_eq!(cohort_schema["properties"]["ploidy"]["minimum"], 1);
    assert_eq!(cohort_schema["properties"]["ploidy"]["maximum"], 10);
    assert_eq!(
        cohort_schema["properties"]["individual_count"]["maximum"],
        json!(u64::MAX)
    );
    assert_eq!(cohort_schema["properties"]["individual_ids"]["minItems"], 1);

    assert_authoring_and_rust(
        "valid-f-diploid",
        json!({
            "population_id": "P", "individual_count": 2, "ploidy": 2,
            "genotype_model": "f_diploid", "inbreeding_coefficient": 0.5,
            "individual_ids": ["sample-A", "sample-B"]
        }),
        &cohort(
            "P",
            2,
            "f_diploid",
            "inbreeding_coefficient = 0.5\nindividual_ids = [\"sample-A\", \"sample-B\"]",
        ),
        true,
        true,
    );
    assert_authoring_and_rust(
        "valid-derived-ids-from-founder-keys",
        json!({
            "population_id": "P", "individual_count": 2, "ploidy": 2,
            "genotype_model": "hwe", "founder_keys": ["founder-A", "founder-B"]
        }),
        &cohort(
            "P",
            2,
            "hwe",
            "founder_keys = [\"founder-A\", \"founder-B\"]",
        ),
        true,
        true,
    );
    assert_authoring_and_rust(
        "missing-ids-and-founder-keys",
        json!({
            "population_id": "P", "individual_count": 2, "ploidy": 2,
            "genotype_model": "hwe"
        }),
        &cohort("P", 2, "hwe", "founder_keys = []").replace("founder_keys = []\n", ""),
        false,
        false,
    );
    assert_authoring_and_rust(
        "provided-ids-and-founder-keys-are-exclusive",
        json!({
            "population_id": "P", "individual_count": 2, "ploidy": 2,
            "genotype_model": "hwe", "individual_ids": ["I1", "I2"],
            "founder_keys": ["founder-A", "founder-B"]
        }),
        &cohort(
            "P",
            2,
            "hwe",
            "individual_ids = [\"I1\", \"I2\"]\nfounder_keys = [\"founder-A\", \"founder-B\"]",
        ),
        false,
        false,
    );
    assert_authoring_and_rust(
        "empty-provided-ids",
        json!({
            "population_id": "P", "individual_count": 2, "ploidy": 2,
            "genotype_model": "hwe", "individual_ids": []
        }),
        &cohort("P", 2, "hwe", "individual_ids = []"),
        false,
        false,
    );
    assert_authoring_and_rust(
        "duplicate-provided-ids",
        json!({
            "population_id": "P", "individual_count": 2, "ploidy": 2,
            "genotype_model": "hwe", "individual_ids": ["same", "same"]
        }),
        &cohort("P", 2, "hwe", "individual_ids = [\"same\", \"same\"]"),
        false,
        false,
    );
    assert_authoring_and_rust(
        "f-outside-diploidy",
        json!({
            "population_id": "P", "individual_count": 2, "ploidy": 3,
            "genotype_model": "f_diploid", "inbreeding_coefficient": 0.5
        }),
        &cohort("P", 3, "f_diploid", "inbreeding_coefficient = 0.5"),
        false,
        false,
    );

    // JSON Schema no puede expresar esta igualdad entre dos campos; Rust sigue
    // siendo la autoridad semántica de cardinalidad.
    assert_authoring_and_rust(
        "provided-id-cardinality-is-rust-semantic",
        json!({
            "population_id": "P", "individual_count": 2, "ploidy": 2,
            "genotype_model": "hwe", "individual_ids": ["only-one"]
        }),
        &cohort("P", 2, "hwe", "individual_ids = [\"only-one\"]"),
        true,
        false,
    );
}

#[test]
fn invalid_local_initial_cohort_stops_run_before_creating_destination() {
    let root = temporary("invalid-initial-cohort-run");
    fs::create_dir(&root).expect("test root");
    let input = root.join("invalid.toml");
    fs::write(
        &input,
        config(&cohort("P", 2, "hwe", "individual_ids = [\"only-one\"]")),
    )
    .expect("local configuration");
    let destination = root.join("result");

    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run"])
        .arg(&input)
        .args(["--output"])
        .arg(&destination)
        .output()
        .expect("run invalid local cohort");

    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("GPS-CONF-008"));
    assert!(!destination.exists());
    fs::remove_dir_all(root).expect("cleanup");
}

fn assert_authoring_and_rust(
    name: &str,
    initial_population: Value,
    toml_population: &str,
    schema_valid: bool,
    rust_valid: bool,
) {
    assert_authoring_schema(name, initial_population, schema_valid);
    assert_eq!(
        build(&config(toml_population)).is_ok(),
        rust_valid,
        "Rust authoring validation mismatch for {name}"
    );
}

fn assert_authoring_schema(name: &str, initial_population: Value, expected_success: bool) {
    let instance = temporary(&format!("authoring-schema-{name}")).with_extension("json");
    fs::write(
        &instance,
        serde_json::to_vec(&json!({
            "config_kind": "general",
            "schema_version": "2",
            "reference": {
                "kind": "simulated",
                "id": "ref-main",
                "composition": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                "sequences": [{"id": "chr1", "class": "chromosome", "length": 10}]
            },
            "initial_populations": [initial_population]
        }))
        .expect("JSON"),
    )
    .expect("instance");
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas/config-general-v2.schema.json"))
        .arg(&instance)
        .output()
        .expect("independent Python validator");
    let _ = fs::remove_file(&instance);
    assert_eq!(
        output.status.success(),
        expected_success,
        "authoring schema mismatch for {name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f2-3-1-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
