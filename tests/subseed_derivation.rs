use std::collections::BTreeMap;

use genome_population_simulator::randomness::{
    DERIVATION_VERSION, RootSeed, SUBSEED_CATALOG_VERSION, SubseedRequest, canonical_material,
    derive_subseed,
};

fn decode_hex(value: &str) -> Vec<u8> {
    value
        .as_bytes()
        .chunks_exact(2)
        .map(|pair| {
            let text = std::str::from_utf8(pair).expect("hexadecimal ASCII");
            u8::from_str_radix(text, 16).expect("hexadecimal válido")
        })
        .collect()
}

#[test]
fn every_a5_vector_matches_literal_material_and_blake3_subseed() {
    for row in include_str!("../catalogs/subseed-scopes-v1.tsv")
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let columns = row.split('|').collect::<Vec<_>>();
        let [scope, component_text, root_hex, material_hex, subseed_hex] = columns.as_slice()
        else {
            panic!("fila A5 mal formada: {row}");
        };
        let owned = component_text
            .split(';')
            .map(|component| component.split_once('=').expect("nombre=valor"))
            .collect::<Vec<_>>();
        let request = SubseedRequest::from_scope(scope, &owned).expect("ámbito A5 conocido");
        let root = RootSeed::parse("root", root_hex).expect("raíz A5");
        assert_eq!(
            canonical_material(root, &request).unwrap(),
            decode_hex(material_hex)
        );
        assert_eq!(
            derive_subseed(root, &request).unwrap().to_hex(),
            *subseed_hex
        );
    }
}

#[test]
fn canonical_encoding_keeps_empty_unicode_lengths_and_normative_order() {
    let root = RootSeed::from_bytes([0; 32]);
    let request = SubseedRequest::from_scope(
        "sampling",
        &[("generation", ""), ("sampling_id", "muestra-ñ")],
    )
    .unwrap();
    let material = canonical_material(root, &request).unwrap();
    let mut expected = vec![0; 32];
    expected.extend_from_slice(&DERIVATION_VERSION.to_be_bytes());
    expected.extend_from_slice(&2_u32.to_be_bytes());
    for (name, value) in [("generation", ""), ("sampling_id", "muestra-ñ")] {
        expected.extend_from_slice(&u32::try_from(name.len()).unwrap().to_be_bytes());
        expected.extend_from_slice(name.as_bytes());
        expected.extend_from_slice(&u32::try_from(value.len()).unwrap().to_be_bytes());
        expected.extend_from_slice(value.as_bytes());
    }
    assert_eq!(material, expected);

    let reversed = SubseedRequest::from_scope(
        "sampling",
        &[("sampling_id", "muestra-ñ"), ("generation", "")],
    )
    .expect_err("el orden del catálogo no es permutable");
    assert_eq!(reversed.code(), "GPS-INTERNAL-001");
    assert_eq!(reversed.exit_code(), 70);
}

#[test]
fn unknown_scope_is_rejected_with_catalog_version() {
    let error = SubseedRequest::from_scope("thread-7", &[]).expect_err("ámbito no normativo");
    assert_eq!(error.code(), "GPS-CONF-018");
    assert_eq!(
        error.context().get("catalog_version").map(String::as_str),
        Some(SUBSEED_CATALOG_VERSION)
    );
}

fn observed_streams(order: &[&str]) -> BTreeMap<String, [u64; 3]> {
    let root = RootSeed::from_bytes([7; 32]);
    order
        .iter()
        .map(|population| {
            let request = SubseedRequest::from_scope(
                "initial_population",
                &[("scenario_id", "scenario-a"), ("population_id", population)],
            )
            .unwrap();
            let mut rng = derive_subseed(root, &request).unwrap().local_rng();
            (
                (*population).to_owned(),
                [rng.next_u64(), rng.next_u64(), rng.next_u64()],
            )
        })
        .collect()
}

#[test]
fn local_chacha12_stream_is_independent_of_added_removed_or_reordered_consumers() {
    let alone = observed_streams(&["observed"]);
    let added_before = observed_streams(&["unrelated-a", "observed", "unrelated-b"]);
    let reordered = observed_streams(&["unrelated-b", "unrelated-a", "observed"]);
    assert_eq!(alone["observed"], added_before["observed"]);
    assert_eq!(alone["observed"], reordered["observed"]);
    assert_ne!(added_before["unrelated-a"], added_before["observed"]);
}
