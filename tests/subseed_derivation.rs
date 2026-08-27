use genome_population_simulator::randomness::{
    DERIVATION_CONTEXT, DERIVATION_VERSION, F1_ACTIVE_SUBSEED_SCOPES, F2_DERIVATION_VERSION,
    RootSeed, RuntimeSubseedScope, SUBSEED_CATALOG_VERSION, derive_subseed,
};
use sha2::{Digest, Sha256};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn decode_hex(value: &str) -> Vec<u8> {
    assert_eq!(value.len() % 2, 0, "hex literal with complete bytes");
    value
        .as_bytes()
        .chunks_exact(2)
        .map(|pair| {
            let text = std::str::from_utf8(pair).expect("hexadecimal ASCII");
            u8::from_str_radix(text, 16).expect("hexadecimal literal")
        })
        .collect()
}

fn encode_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn independent_material(
    root: &[u8],
    derivation_version: u32,
    scope_name: Option<&str>,
    components: &[(&str, &str)],
) -> Vec<u8> {
    let mut result = root.to_vec();
    result.extend_from_slice(&derivation_version.to_be_bytes());
    if let Some(scope_name) = scope_name {
        result.extend_from_slice(
            &u32::try_from(scope_name.len())
                .expect("small literal scope name")
                .to_be_bytes(),
        );
        result.extend_from_slice(scope_name.as_bytes());
    }
    result.extend_from_slice(
        &u32::try_from(components.len())
            .expect("small literal component count")
            .to_be_bytes(),
    );
    for (name, value) in components {
        for bytes in [name.as_bytes(), value.as_bytes()] {
            result.extend_from_slice(
                &u32::try_from(bytes.len())
                    .expect("small literal component")
                    .to_be_bytes(),
            );
            result.extend_from_slice(bytes);
        }
    }
    result
}

fn catalog_scope_derivation_version(catalog: &str, scope_name: &str) -> u32 {
    let mapping = catalog
        .lines()
        .find_map(|line| line.strip_prefix("# scope_derivation_versions="))
        .expect("v4 declares each row derivation version");
    mapping
        .split(';')
        .find_map(|entry| {
            let (scope, version) = entry.split_once(':').expect("scope:version");
            (scope == scope_name).then(|| version.parse::<u32>().expect("u32 version"))
        })
        .expect("every catalog row has a declared derivation version")
}

fn verify_catalog_with_independent_oracle(catalog: &str, derivation_version: u32) {
    for row in catalog
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let columns = row.split('|').collect::<Vec<_>>();
        let [
            scope_name,
            component_text,
            root_hex,
            material_hex,
            subseed_hex,
        ] = columns.as_slice()
        else {
            panic!("malformed frozen catalog row: {row}");
        };
        let components = component_text
            .split(';')
            .map(|component| component.split_once('=').expect("name=value"))
            .collect::<Vec<_>>();
        let root = decode_hex(root_hex);
        let row_derivation_version = if derivation_version == F2_DERIVATION_VERSION {
            catalog_scope_derivation_version(catalog, scope_name)
        } else {
            derivation_version
        };
        assert!(
            [DERIVATION_VERSION, F2_DERIVATION_VERSION].contains(&row_derivation_version),
            "unsupported row derivation version"
        );
        let is_v2 = row_derivation_version == F2_DERIVATION_VERSION;
        let material = independent_material(
            &root,
            row_derivation_version,
            is_v2.then_some(*scope_name),
            &components,
        );
        assert_eq!(encode_hex(&material), *material_hex);
        assert_eq!(
            encode_hex(&blake3::derive_key(DERIVATION_CONTEXT, &material)),
            *subseed_hex
        );
        for (mutated_index, (name, value)) in components.iter().enumerate() {
            let mut mutated = components.clone();
            let mutated_value = format!("{value}-mutated");
            mutated[mutated_index] = (name, &mutated_value);
            let mutated_material = independent_material(
                &root,
                row_derivation_version,
                is_v2.then_some(*scope_name),
                &mutated,
            );
            assert_ne!(
                encode_hex(&blake3::derive_key(DERIVATION_CONTEXT, &mutated_material)),
                *subseed_hex,
                "mutating {scope_name}.{name} must isolate its subseed"
            );
        }
    }
}

#[test]
fn published_catalogs_are_immutable_and_v8_has_an_independent_oracle() {
    let catalogs = [
        (
            include_bytes!("../catalogs/subseed-scopes-v1.tsv").as_slice(),
            "b6c56848b1ae48d5646dcdc515709259ea65fc5c8b6a407bf403e4da8477675d",
        ),
        (
            include_bytes!("../catalogs/subseed-scopes-v2.tsv").as_slice(),
            "c40f812d1645bb11a2e2641cc76a5b6db583f0505e246cb41a24bc0aa554ceff",
        ),
        (
            include_bytes!("../catalogs/subseed-scopes-v3.tsv").as_slice(),
            "415ead5f7083220953f3795f51da1cc464c41783174d11359d3e0e5759460066",
        ),
    ];
    for (bytes, expected_sha256) in catalogs {
        assert_eq!(encode_hex(&Sha256::digest(bytes)), expected_sha256);
        verify_catalog_with_independent_oracle(
            std::str::from_utf8(bytes).expect("published catalogs are UTF-8"),
            DERIVATION_VERSION,
        );
    }
    let v4 = include_str!("../catalogs/subseed-scopes-v4.tsv");
    assert_eq!(
        encode_hex(&Sha256::digest(v4.as_bytes())),
        "c14c582af00d309580846264184baeafc0194891c3795e3af05bc62a925fb74e"
    );
    assert!(v4.starts_with(
        "# derivation_versions=1,2\n# catalog_version=4\n# context=genome_population_simulator/v1/subseed\n"
    ));
    verify_catalog_with_independent_oracle(v4, F2_DERIVATION_VERSION);
    let v5 = include_str!("../catalogs/subseed-scopes-v5.tsv");
    assert!(v5.starts_with(
        "# derivation_versions=1,2\n# catalog_version=5\n# context=genome_population_simulator/v1/subseed\n"
    ));
    let v4_rows = v4
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect::<Vec<_>>();
    let v5_rows = v5
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(&v5_rows[..v4_rows.len()], v4_rows);
    assert_eq!(v5_rows.len(), v4_rows.len() + 1);
    verify_catalog_with_independent_oracle(v5, F2_DERIVATION_VERSION);

    let v6 = include_str!("../catalogs/subseed-scopes-v6.tsv");
    assert!(v6.starts_with(
        "# derivation_versions=1,2\n# catalog_version=6\n# context=genome_population_simulator/v1/subseed\n"
    ));
    let v6_rows = v6
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(v6_rows.len(), v5_rows.len());
    for (v5_row, v6_row) in v5_rows.iter().zip(&v6_rows) {
        if v5_row.starts_with("reference_reduction|") {
            assert!(v6_row.starts_with(
                "reference_reduction|source_reference_id=reference-main;derived_reference_id=reference-small;selector=direct_proportional_selection|"
            ));
        } else {
            assert_eq!(v6_row, v5_row);
        }
    }
    verify_catalog_with_independent_oracle(v6, F2_DERIVATION_VERSION);

    let v7 = include_str!("../catalogs/subseed-scopes-v7.tsv");
    assert_eq!(
        encode_hex(&Sha256::digest(v7.as_bytes())),
        "bddd45fdb968e312676028585333b39a05c331cdcddae6d8bffda6edad05dbf8"
    );
    assert!(v7.starts_with(
        "# derivation_versions=1,2\n# catalog_version=7\n# context=genome_population_simulator/v1/subseed\n"
    ));
    let v7_rows = v7
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(&v7_rows[..v6_rows.len()], v6_rows);
    assert_eq!(v7_rows.len(), v6_rows.len() + 1);
    assert!(
        v7_rows
            .last()
            .is_some_and(|row| row.starts_with("technical_template|"))
    );
    verify_catalog_with_independent_oracle(v7, F2_DERIVATION_VERSION);

    let v8 = include_str!("../catalogs/subseed-scopes-v8.tsv");
    assert_eq!(
        encode_hex(&Sha256::digest(v8.as_bytes())),
        "cdf431dd9cc63858788d9377cc914506013f134052b45171cba827c69d2ddbbd"
    );
    assert!(v8.starts_with(
        "# derivation_versions=1,2\n# catalog_version=8\n# context=genome_population_simulator/v1/subseed\n"
    ));
    let v8_rows = v8
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(&v8_rows[..v7_rows.len()], v7_rows);
    assert_eq!(v8_rows.len(), v7_rows.len() + 1);
    assert!(
        v8_rows
            .last()
            .is_some_and(|row| row.starts_with("technical_library_count|"))
    );
    verify_catalog_with_independent_oracle(v8, F2_DERIVATION_VERSION);
}

#[test]
fn f1_runtime_scopes_remain_on_derivation_v1_and_match_active_catalog_v8() {
    let root = RootSeed::parse("root", ROOT).expect("frozen 256-bit root");
    assert_eq!(DERIVATION_VERSION, 1);
    assert_eq!(SUBSEED_CATALOG_VERSION, "8");
    assert_eq!(
        derive_subseed(
            root,
            RuntimeSubseedScope::ReferenceGeneration {
                reference_id: "reference-main",
            },
        )
        .expect("generation subseed")
        .to_hex(),
        "8a933c0e70340b518294b23a825ae36c75b3160ca824fa739ca19adcb47fd5a4"
    );
    assert_eq!(
        derive_subseed(
            root,
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: "reference-main",
                derived_reference_id: "reference-small",
            },
        )
        .expect("reduction subseed")
        .to_hex(),
        "20c4a19ada8e8cdfce862db5b3423fd8250be2a61cf83321b4a6fb9210963169"
    );
    let active_row = include_str!("../catalogs/subseed-scopes-v8.tsv")
        .lines()
        .find(|row| row.starts_with("reference_reduction|"))
        .expect("active reduction catalog row");
    assert!(active_row.contains("selector=direct_proportional_selection"));
    assert!(
        active_row.ends_with("20c4a19ada8e8cdfce862db5b3423fd8250be2a61cf83321b4a6fb9210963169")
    );
}

fn v4_catalog_subseed(scope_name: &str) -> &str {
    include_str!("../catalogs/subseed-scopes-v4.tsv")
        .lines()
        .find_map(|row| {
            let fields = row.split('|').collect::<Vec<_>>();
            if fields.len() == 5 && fields[0] == scope_name {
                Some(fields[4])
            } else {
                None
            }
        })
        .expect("frozen v4 vector")
}

fn v5_catalog_subseed(scope_name: &str) -> &str {
    include_str!("../catalogs/subseed-scopes-v5.tsv")
        .lines()
        .find_map(|row| {
            let fields = row.split('|').collect::<Vec<_>>();
            if fields.len() == 5 && fields[0] == scope_name {
                Some(fields[4])
            } else {
                None
            }
        })
        .expect("frozen v5 vector")
}

#[test]
fn every_f2_runtime_scope_matches_its_frozen_v4_or_v5_vector() {
    let root = RootSeed::parse("root", ROOT).expect("frozen 256-bit root");
    let vectors = [
        (
            "matrix_candidate_count",
            RuntimeSubseedScope::MatrixCandidateCount {
                scenario_id: "direct",
                reference_id: "reference-main",
            },
        ),
        (
            "matrix_locus_selection",
            RuntimeSubseedScope::MatrixLocusSelection {
                scenario_id: "direct",
                reference_id: "reference-main",
            },
        ),
        (
            "matrix_alternate_allele",
            RuntimeSubseedScope::MatrixAlternateAllele {
                scenario_id: "direct",
                reference_id: "reference-main",
                sequence_id: "sequence-main",
                position_1based: 42,
                reference_base: "A",
            },
        ),
        (
            "matrix_exact_assignment",
            RuntimeSubseedScope::MatrixExactAssignment {
                scenario_id: "direct",
                reference_id: "reference-main",
                population_id: "population-main",
            },
        ),
        (
            "matrix_frequency_component",
            RuntimeSubseedScope::MatrixFrequencyComponent {
                scenario_id: "direct",
                reference_id: "reference-main",
                population_id: "population-main",
                variant_id: "variant-main",
            },
        ),
        (
            "matrix_frequency_value",
            RuntimeSubseedScope::MatrixFrequencyValue {
                scenario_id: "direct",
                reference_id: "reference-main",
                population_id: "population-main",
                variant_id: "variant-main",
            },
        ),
        (
            "indel_count",
            RuntimeSubseedScope::IndelCount {
                scenario_id: "direct",
                reference_id: "reference-main",
                individual_id: "individual-main",
                sequence_id: "sequence-main",
                variant_class: "insertion",
                spatial_scope_id: "global",
            },
        ),
        (
            "indel_position",
            RuntimeSubseedScope::IndelPosition {
                scenario_id: "direct",
                reference_id: "reference-main",
                individual_id: "individual-main",
                sequence_id: "sequence-main",
                variant_class: "insertion",
                spatial_scope_id: "global",
            },
        ),
        (
            "insertion_content",
            RuntimeSubseedScope::InsertionContent {
                scenario_id: "direct",
                reference_id: "reference-main",
                individual_id: "individual-main",
                sequence_id: "sequence-main",
                spatial_scope_id: "global",
                event_index: 0,
            },
        ),
        (
            "initial_genotype",
            RuntimeSubseedScope::InitialGenotype {
                scenario_id: "direct",
                reference_id: "reference-main",
                population_id: "population-main",
                individual_id: "individual-main",
                sequence_id: "sequence-main",
                position_1based: 42,
                reference_base: "A",
                alternate_alleles: "C,G",
            },
        ),
        (
            "initial_phase",
            RuntimeSubseedScope::InitialPhase {
                scenario_id: "direct",
                reference_id: "reference-main",
                population_id: "population-main",
                individual_id: "individual-main",
                variant_id: "variant-main",
            },
        ),
    ];
    for (scope_name, scope) in vectors {
        assert_eq!(
            derive_subseed(root, scope).expect("F2 subseed").to_hex(),
            v4_catalog_subseed(scope_name),
            "{scope_name}"
        );
    }
    assert_eq!(
        derive_subseed(
            root,
            RuntimeSubseedScope::InitialPhaseLocus {
                scenario_id: "direct",
                reference_id: "reference-main",
                population_id: "population-main",
                individual_id: "individual-main",
                sequence_id: "sequence-main",
                position_1based: 42,
                reference_base: "C",
                alternate_alleles: "A,G",
            },
        )
        .expect("F2.3 phase subseed")
        .to_hex(),
        v5_catalog_subseed("initial_phase_locus")
    );
}

#[test]
fn derivation_is_repeatable_and_generation_and_reduction_are_isolated() {
    let root = RootSeed::from_bytes([11; 32]);
    let generation = |reference_id| {
        derive_subseed(
            root,
            RuntimeSubseedScope::ReferenceGeneration { reference_id },
        )
        .expect("generation subseed")
    };
    let reduction = |derived_reference_id| {
        derive_subseed(
            root,
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: "reference-main",
                derived_reference_id,
            },
        )
        .expect("reduction subseed")
    };

    assert_eq!(generation("reference-main"), generation("reference-main"));
    assert_eq!(reduction("reference-small"), reduction("reference-small"));

    let unchanged_reduction = reduction("reference-small");
    assert_ne!(generation("reference-main"), generation("reference-other"));
    assert_eq!(reduction("reference-small"), unchanged_reduction);

    let unchanged_generation = generation("reference-main");
    assert_ne!(reduction("reference-small"), reduction("reference-other"));
    assert_eq!(generation("reference-main"), unchanged_generation);
}

#[test]
fn every_active_scope_owns_a_fresh_local_chacha12_stream() {
    let root = RootSeed::from_bytes([7; 32]);
    let generation = derive_subseed(
        root,
        RuntimeSubseedScope::ReferenceGeneration {
            reference_id: "reference-main",
        },
    )
    .expect("generation subseed");
    let reduction = derive_subseed(
        root,
        RuntimeSubseedScope::ReferenceReduction {
            source_reference_id: "reference-main",
            derived_reference_id: "reference-small",
        },
    )
    .expect("reduction subseed");

    let mut first = generation.local_rng();
    let mut repeated = generation.local_rng();
    let generation_values = [first.next_u64(), first.next_u64(), first.next_u64()];
    assert_eq!(
        generation_values,
        [
            repeated.next_u64(),
            repeated.next_u64(),
            repeated.next_u64(),
        ]
    );
    let mut reduction_rng = reduction.local_rng();
    assert_ne!(generation_values[0], reduction_rng.next_u64());
}

#[test]
fn the_f1_runtime_catalog_is_closed_to_exactly_two_consumers() {
    assert_eq!(
        F1_ACTIVE_SUBSEED_SCOPES,
        ["reference_generation", "reference_reduction"]
    );
}
