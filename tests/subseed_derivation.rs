use genome_population_simulator::randomness::{
    DERIVATION_CONTEXT, DERIVATION_VERSION, F1_ACTIVE_SUBSEED_SCOPES,
    REFERENCE_REDUCTION_ALGORITHM, RootSeed, RuntimeSubseedScope, SUBSEED_CATALOG_VERSION,
    derive_subseed,
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

fn independent_material(root: &[u8], components: &[(&str, &str)]) -> Vec<u8> {
    let mut result = root.to_vec();
    result.extend_from_slice(&1_u32.to_be_bytes());
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

fn verify_catalog_with_independent_oracle(catalog: &str) {
    for row in catalog
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let columns = row.split('|').collect::<Vec<_>>();
        let [_, component_text, root_hex, material_hex, subseed_hex] = columns.as_slice() else {
            panic!("malformed frozen catalog row: {row}");
        };
        let components = component_text
            .split(';')
            .map(|component| component.split_once('=').expect("name=value"))
            .collect::<Vec<_>>();
        let root = decode_hex(root_hex);
        let material = independent_material(&root, &components);
        assert_eq!(encode_hex(&material), *material_hex);
        assert_eq!(
            encode_hex(&blake3::derive_key(DERIVATION_CONTEXT, &material)),
            *subseed_hex
        );
    }
}

#[test]
fn published_v1_v3_catalogs_remain_byte_exact_and_all_vectors_have_an_independent_oracle() {
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
        );
    }
}

#[test]
fn both_active_runtime_scopes_match_their_frozen_v3_vectors() {
    let root = RootSeed::parse("root", ROOT).expect("frozen 256-bit root");
    assert_eq!(DERIVATION_VERSION, 1);
    assert_eq!(SUBSEED_CATALOG_VERSION, "3");
    assert_eq!(
        REFERENCE_REDUCTION_ALGORITHM,
        "reduced_reference_selection_v1"
    );
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
        "7fa724d52b2db97a043920dbc3062a34c6377eab76f1d9d6e4ff1d513747decd"
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
