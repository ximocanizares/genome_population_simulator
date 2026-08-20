use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use genome_population_simulator::application::load_and_validate_config_with_entropy;
use genome_population_simulator::randomness::{
    DERIVATION_VERSION, EntropySource, ReproducibilityMetadata, RootSeed, RuntimeSubseedScope,
    SUBSEED_CATALOG_VERSION, SeedOrigin, derive_subseed,
};

#[derive(Debug)]
struct FixedEntropy([u8; 32]);

impl EntropySource for FixedEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        *destination = self.0;
        Ok(())
    }
}

fn root() -> RootSeed {
    RootSeed::from_bytes(core::array::from_fn(|index| index as u8))
}

#[test]
fn f1_3_matrix_closes_one_requirement_and_marks_three_results_as_partial() {
    let matrix = include_str!("../evidence/f1-reimplementation/f1.3-test-req-matrix.tsv");
    let source = include_str!("f1_3_traceability.rs");
    let traceability = include_str!("../design/trazabilidad-requisitos-0.1.md");
    let rows = matrix
        .lines()
        .skip(1)
        .map(|row| row.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();

    assert_eq!(rows.len(), 4);
    assert_eq!(
        rows.iter().map(|row| row[1]).collect::<BTreeSet<_>>(),
        BTreeSet::from(["REQ-CONF-26", "REQ-CONF-27", "REQ-REP-22", "REQ-REP-23"])
    );
    for row in rows {
        assert_eq!(row.len(), 6);
        let symbol = row[4]
            .split_once("::")
            .expect("Rust evidence with a symbol")
            .1;
        assert!(source.contains(&format!("fn {symbol}(")));
        let traceability_row = traceability
            .lines()
            .find(|line| line.starts_with(&format!("| {} |", row[1])))
            .expect("cada requisito de la matriz aparece una vez en trazabilidad");
        let traceability_columns = traceability_row
            .split('|')
            .map(str::trim)
            .collect::<Vec<_>>();
        assert_eq!(
            traceability_columns[3], row[3],
            "hito divergente para {}",
            row[1]
        );
        assert!(
            traceability_columns[4].contains(symbol),
            "símbolo de evidencia divergente para {}",
            row[1]
        );
        if row[1] == "REQ-REP-23" {
            assert_eq!((row[2], row[3]), ("pass", "F1.3"));
            assert!(!row[5].starts_with("evidencia parcial:"));
        } else {
            let expected_milestone = match row[1] {
                "REQ-CONF-26" => "F1.5",
                "REQ-CONF-27" | "REQ-REP-22" => "F9.7",
                unexpected => panic!("requisito parcial inesperado: {unexpected}"),
            };
            assert_eq!((row[2], row[3]), ("reassigned", expected_milestone));
            assert!(row[5].starts_with("evidencia parcial:"));
            assert!(symbol.starts_with("partial_test_req_"));
        }
    }
}

#[test]
fn partial_test_req_conf_26_root_is_effective_before_any_f1_rng_scope() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/config/general-minimal.toml");
    let expected = [0xa5; 32];
    let effective = load_and_validate_config_with_entropy(&path, &mut FixedEntropy(expected))
        .expect("effective configuration with generated root");
    let metadata = effective.reproducibility();
    assert_eq!(metadata.root_seed(), RootSeed::from_bytes(expected));
    assert_eq!(metadata.seed_origin(), SeedOrigin::SystemCsprng);

    let subseed = derive_subseed(
        metadata.root_seed(),
        RuntimeSubseedScope::ReferenceGeneration {
            reference_id: effective.reference_id(),
        },
    )
    .expect("the active scope is constructed from the effective root");
    assert_ne!(subseed.as_bytes(), &[0; 32]);
}

#[test]
fn partial_test_req_conf_27_active_scope_derivation_is_logically_repeatable() {
    let generation = RuntimeSubseedScope::ReferenceGeneration {
        reference_id: "reference-main",
    };
    let reduction = RuntimeSubseedScope::ReferenceReduction {
        source_reference_id: "reference-main",
        derived_reference_id: "reference-small",
    };
    for scope in [generation, reduction] {
        let first = derive_subseed(root(), scope).expect("first derivation");
        let second = derive_subseed(root(), scope).expect("repeated derivation");
        assert_eq!(first, second);

        let mut first_rng = first.local_rng();
        let mut second_rng = second.local_rng();
        assert_eq!(first_rng.next_u64(), second_rng.next_u64());
        assert_eq!(first_rng.next_u64(), second_rng.next_u64());
    }
}

#[test]
fn partial_test_req_rep_22_effective_metadata_identifies_the_rng_contract() {
    let metadata = ReproducibilityMetadata::new(root(), SeedOrigin::Declared);
    assert_eq!(metadata.root_seed(), root());
    assert_eq!(metadata.seed_origin(), SeedOrigin::Declared);
    assert_eq!(metadata.random_algorithm(), "ChaCha12");
    assert_eq!(metadata.derivation_version(), DERIVATION_VERSION);
    assert_eq!(metadata.catalog_version(), SUBSEED_CATALOG_VERSION);
}

#[derive(Clone, Copy)]
struct RngContract<'a> {
    simulator_version: &'a str,
    random_algorithm: &'a str,
    derivation_version: u32,
    catalog_version: &'a str,
}

fn same_reproduction_contract(first: RngContract<'_>, second: RngContract<'_>) -> bool {
    first.simulator_version == second.simulator_version
        && first.random_algorithm == second.random_algorithm
        && first.derivation_version == second.derivation_version
        && first.catalog_version == second.catalog_version
}

#[test]
fn test_req_rep_23_requires_exact_simulator_and_rng_contract_versions() {
    let current = RngContract {
        simulator_version: env!("CARGO_PKG_VERSION"),
        random_algorithm: "ChaCha12",
        derivation_version: DERIVATION_VERSION,
        catalog_version: SUBSEED_CATALOG_VERSION,
    };
    assert!(same_reproduction_contract(current, current));

    for changed in [
        RngContract {
            simulator_version: "another-version",
            ..current
        },
        RngContract {
            random_algorithm: "another-algorithm",
            ..current
        },
        RngContract {
            derivation_version: DERIVATION_VERSION + 1,
            ..current
        },
        RngContract {
            catalog_version: "another-catalog",
            ..current
        },
    ] {
        assert!(!same_reproduction_contract(current, changed));
    }
}

#[test]
fn f1_runtime_has_no_preparatory_rng_scope_or_event_comparator() {
    fn append_rust_sources(path: &Path, destination: &mut String) {
        let mut entries = fs::read_dir(path)
            .expect("directorio fuente legible")
            .collect::<Result<Vec<_>, _>>()
            .expect("entradas fuente legibles");
        entries.sort_by_key(|entry| entry.path());
        for entry in entries {
            let path = entry.path();
            if path.is_dir() {
                append_rust_sources(&path, destination);
            } else if path.extension().and_then(|value| value.to_str()) == Some("rs") {
                destination.push_str(&fs::read_to_string(path).expect("fuente Rust UTF-8"));
            }
        }
    }

    let mut runtime = String::new();
    append_rust_sources(
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("src"),
        &mut runtime,
    );
    for forbidden in [
        "\"initial_population\"",
        "\"biological_process\"",
        "\"sampling\"",
        "\"technical_event\"",
        "EventSubseed",
        "event_subseeds",
        "compare_biological",
        "compare_technical",
        "compare_semantic",
        "compare_binary",
    ] {
        assert!(
            !runtime.contains(forbidden),
            "forbidden F1 runtime surface: {forbidden}"
        );
    }
}
