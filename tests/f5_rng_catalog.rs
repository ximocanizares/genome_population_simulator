use std::collections::BTreeMap;

use genome_population_simulator::identity::{
    F5BlockId, F5SiteId, IndividualId, PopulationId, ScenarioId,
};
use genome_population_simulator::randomness::{
    F5_RNG_CATALOG, F5_UNIFORM_INDEX_ALGORITHM, F5RandomDecision, RootSeed, derive_f5_subseed,
    f5_bernoulli, f5_fisher_yates, f5_matrix_bytes, f5_placed_layout_bytes, f5_uniform_index,
};

const F5_RNG_CATALOG_TSV: &str = include_str!("../catalogs/f5-rng-catalog-v1.tsv");

fn root() -> RootSeed {
    RootSeed::from_bytes([
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e,
        0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d,
        0x1e, 0x1f,
    ])
}

fn expected_rows() -> BTreeMap<&'static str, (&'static str, [u64; 4])> {
    let mut rows = BTreeMap::new();
    for line in F5_RNG_CATALOG_TSV
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let fields = line.split('|').collect::<Vec<_>>();
        assert_eq!(fields.len(), 9);
        let words = std::array::from_fn(|index| {
            let encoded = fields[index + 5];
            assert_eq!(encoded.len(), 16);
            assert!(
                encoded
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
            );
            u64::from_str_radix(encoded, 16).expect("u64 vector must be lowercase hex")
        });
        rows.insert(fields[0], (fields[4], words));
    }
    rows
}

#[test]
fn f5_rng_catalog_constants_name_the_frozen_contract() {
    assert_eq!(F5_RNG_CATALOG, "f5-rng-catalog-v1");
    assert_eq!(F5_UNIFORM_INDEX_ALGORITHM, "uniform-index-v1");
    assert!(F5_RNG_CATALOG_TSV.starts_with("# catalog_id=f5-rng-catalog-v1\n"));
}

#[test]
fn f5_rng_catalog_v1_has_literal_subseed_and_chacha12_vectors() {
    let scenario_id = ScenarioId::parse("scenario-f5", "test").expect("scenario");
    let population_id = PopulationId::parse("population-A", "test").expect("population");
    let block_id = F5BlockId::parse("BLK_blockA", "test").expect("block");
    let site_id = F5SiteId::parse("SITE_siteA", "test").expect("site");
    let individual_id = IndividualId::parse("individual-A", "test").expect("individual");
    let decisions = [
        (
            "layout_length",
            F5RandomDecision::layout_length(&scenario_id, 2, 7),
        ),
        (
            "layout_position",
            F5RandomDecision::layout_position(&scenario_id, 2, 7),
        ),
        (
            "layout_gap",
            F5RandomDecision::layout_gap(&scenario_id, 2, 7),
        ),
        (
            "site_alt",
            F5RandomDecision::site_alt(&scenario_id, &block_id, &site_id),
        ),
        (
            "marginal_spectrum",
            F5RandomDecision::marginal_spectrum(&scenario_id, &population_id, &block_id, &site_id),
        ),
        (
            "initial_column",
            F5RandomDecision::initial_column(&scenario_id, &population_id, &block_id, 3, 5),
        ),
        (
            "polarity",
            F5RandomDecision::polarity(&scenario_id, &population_id, &block_id, 3, 5),
        ),
        (
            "swap_proposal",
            F5RandomDecision::swap_proposal(&scenario_id, &population_id, &block_id, 3),
        ),
        (
            "founder_sample",
            F5RandomDecision::founder_sample(
                &scenario_id,
                &population_id,
                &block_id,
                &individual_id,
                1,
            ),
        ),
        (
            "founder_f",
            F5RandomDecision::founder_f(&scenario_id, &population_id, &block_id, &individual_id),
        ),
        (
            "controlled_quota_tie",
            F5RandomDecision::controlled_quota_tie(&scenario_id, &population_id, &block_id),
        ),
        (
            "founder_controlled_permutation",
            F5RandomDecision::founder_controlled_permutation(
                &scenario_id,
                &population_id,
                &block_id,
            ),
        ),
    ];
    let expected = expected_rows();
    assert_eq!(expected.len(), decisions.len());

    for (decision_name, decision) in decisions {
        let (expected_subseed, expected_words) =
            expected.get(decision_name).expect("decision row exists");
        let subseed = derive_f5_subseed(root(), decision).expect("subseed");
        assert_eq!(subseed.to_hex(), *expected_subseed, "{decision_name}");
        let mut rng = subseed.local_rng();
        let observed_words = std::array::from_fn(|_| rng.next_u64());
        assert_eq!(observed_words, *expected_words, "{decision_name}");
    }
}

#[test]
fn f5_runtime_primitives_preserve_extreme_consumption_contracts() {
    let scenario_id = ScenarioId::parse("scenario-f5", "test").expect("scenario");
    let subseed = derive_f5_subseed(root(), F5RandomDecision::layout_length(&scenario_id, 2, 7))
        .expect("subseed");
    let mut singleton = subseed.local_rng();
    assert_eq!(f5_uniform_index(&mut singleton, 1), Ok(0));
    assert_eq!(
        singleton.next_u64(),
        subseed.local_rng().next_u64(),
        "n=1 must not consume the F5 stream"
    );

    let mut endpoints = subseed.local_rng();
    assert_eq!(f5_bernoulli(&mut endpoints, 0.0), Ok(false));
    assert_eq!(f5_bernoulli(&mut endpoints, 1.0), Ok(true));
    assert_eq!(
        endpoints.next_u64(),
        subseed.local_rng().next_u64(),
        "Bernoulli endpoints must not consume the F5 stream"
    );

    let mut shuffle = subseed.local_rng();
    let mut values = [0, 1, 2, 3];
    f5_fisher_yates(&mut shuffle, &mut values).expect("shuffle");
    assert_eq!(values, [0, 2, 3, 1]);
}

#[test]
fn f5_layout_and_matrix_bytes_are_canonical_independent_of_later_writers() {
    assert_eq!(
        f5_placed_layout_bytes(&[(1, 2, 5, 4)]),
        [1_u64, 2, 5, 4]
            .into_iter()
            .flat_map(u64::to_be_bytes)
            .collect::<Vec<u8>>()
    );
    assert_eq!(
        f5_matrix_bytes(&[
            vec![true, false, true, false, false, false, false, true, true],
            vec![false, true, false, true, false, true, false, true, false],
        ]),
        Ok(vec![
            0, 0, 0, 0, 0, 0, 0, 2, //
            0, 0, 0, 0, 0, 0, 0, 9, //
            0xa1, 0x80, 0x55, 0x00,
        ])
    );
}
