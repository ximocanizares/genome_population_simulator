use genome_population_simulator::identity::{EventId, IndividualId, LibraryId};
use genome_population_simulator::randomness::{
    F4IndexHopInstance, F4Mate, F4TechnicalDecision, RootSeed, derive_f4_technical_subseed,
    f4_uniform_closed,
};
use std::num::NonZeroU64;

fn root() -> RootSeed {
    RootSeed::from_bytes([
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e,
        0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d,
        0x1e, 0x1f,
    ])
}

fn inputs() -> (EventId, LibraryId, IndividualId) {
    (
        EventId::parse("event-4", "test").expect("event"),
        LibraryId::parse("library-4", "test").expect("library"),
        IndividualId::parse("sample-4", "test").expect("individual"),
    )
}

fn pre_hop_instance(value: &str) -> F4IndexHopInstance<'_> {
    F4IndexHopInstance::parse(value).expect("canonical pre-hop instance")
}

#[test]
fn f4_technical_v3_has_literal_subseed_and_chacha12_vectors() {
    let (event_id, library_id, recipient_id) = inputs();
    let decisions = [
        F4TechnicalDecision::fragment_insert_length(&event_id, &library_id, 17),
        F4TechnicalDecision::fragment_copy(&event_id, &library_id, 17),
        F4TechnicalDecision::fragment_start(&event_id, &library_id, 17),
        F4TechnicalDecision::chimera_bernoulli(&event_id, &library_id, 17),
        F4TechnicalDecision::chimera_breakpoint(&event_id, &library_id, 17),
        F4TechnicalDecision::chimera_second_source(
            &event_id,
            &library_id,
            17,
            NonZeroU64::new(43).expect("nonzero breakpoint"),
        ),
        F4TechnicalDecision::error_bernoulli(&event_id, &library_id, 17, F4Mate::R1, 149),
        F4TechnicalDecision::error_substitution(&event_id, &library_id, 17, F4Mate::R2, 149),
        F4TechnicalDecision::contamination_percent(&event_id, &library_id, &recipient_id),
        F4TechnicalDecision::contamination_donors(&event_id, &library_id, &recipient_id),
        F4TechnicalDecision::contamination_recipients(&event_id, &library_id, &recipient_id),
        F4TechnicalDecision::index_hop_bernoulli(
            &event_id,
            &library_id,
            pre_hop_instance("PH_pre-hop-4"),
        ),
        F4TechnicalDecision::index_hop_destination(
            &event_id,
            &library_id,
            pre_hop_instance("PH_pre-hop-4"),
        ),
    ];

    let expected = [
        (
            "43067d5742e5d7fb984f922f24fca19dec22f0b6b0af4e521be686e847e563ab",
            0x42ee04d0125b20c5,
        ),
        (
            "969f08953f5ee5b287fc835c6c6365bd1d61b2067e759a998a5d8ae8a47bebbf",
            0x8c8772b47d9d1998,
        ),
        (
            "bcad565f4329918c146c4a9195f7d7054594dfb89449de50b69ca1b79dca64c5",
            0x5bb645578eb00fcf,
        ),
        (
            "851c204d515c63b94e4505499e477be2ed2dbbac3d825c53dee02ae9a545d729",
            0x15609b5fd11adaff,
        ),
        (
            "b71b71ce2fd14bcce10be376e23d5f7fbfe4deadd2ea050b8d26436d3fda6296",
            0xe58135cdcfedee09,
        ),
        (
            "51873105f524b76797c314d538349674e22be1a539bfdc892b9de4fe65de364d",
            0xd60fc0eba9f8b1f3,
        ),
        (
            "86bada22301869d99674d7b4b4e0496c97efc3eacd922791a375b860eed602e9",
            0x751a660bc3a08f37,
        ),
        (
            "8870de6e4e0846702c87ba822e921edfd32afea8a97d0af78e0c86e8027e866b",
            0xc7a52771a101fd7b,
        ),
        (
            "7d3e6d6a82e99fc67a666f2582206a790e2f981230e46364e9f7c1caf05ff23b",
            0x21a7c242e38e2db8,
        ),
        (
            "43c9344c239f7c9a47f03820be07b369cadacb2cbf1ddb448e0931a4e7c599bf",
            0xa2d66c98b07449c7,
        ),
        (
            "934893e2592ef921a1800c5b3a4e4c1b2b6537f4583094e08e672d5847169e87",
            0xdeb97793dcf1b926,
        ),
        (
            "b306deae558fbb4ad281ada8993640fca9593d466d1a0842e4555c3019cbbbec",
            0x15e54ee714e8e7a5,
        ),
        (
            "2d0ee64e51a2582c9972a583d0ec57b2d2ad0d598c836f1a492e8b97f8202970",
            0xfc6fdac77aa4aaac,
        ),
    ];
    for (decision, (expected_subseed, expected_word)) in decisions.into_iter().zip(expected) {
        let subseed = derive_f4_technical_subseed(root(), decision).expect("subseed");
        assert_eq!(subseed.to_hex(), expected_subseed);
        assert_eq!(subseed.local_rng().next_u64(), expected_word);
    }
}

#[test]
fn decisions_regenerate_and_are_isolated_from_component_changes() {
    let (event_id, library_id, recipient_id) = inputs();
    let baseline =
        F4TechnicalDecision::error_bernoulli(&event_id, &library_id, 17, F4Mate::R1, 149);
    let repeated = derive_f4_technical_subseed(root(), baseline).expect("repeat");
    assert_eq!(
        repeated,
        derive_f4_technical_subseed(root(), baseline).expect("same cause")
    );

    let changed_cycle =
        F4TechnicalDecision::error_bernoulli(&event_id, &library_id, 17, F4Mate::R1, 148);
    assert_ne!(
        repeated,
        derive_f4_technical_subseed(root(), changed_cycle).expect("cycle")
    );
    assert_ne!(
        derive_f4_technical_subseed(
            root(),
            F4TechnicalDecision::contamination_donors(&event_id, &library_id, &recipient_id),
        )
        .expect("donors"),
        derive_f4_technical_subseed(
            root(),
            F4TechnicalDecision::contamination_recipients(&event_id, &library_id, &recipient_id),
        )
        .expect("recipients")
    );
    assert_ne!(
        derive_f4_technical_subseed(
            root(),
            F4TechnicalDecision::index_hop_bernoulli(
                &event_id,
                &library_id,
                pre_hop_instance("PH_a"),
            ),
        )
        .expect("hop a"),
        derive_f4_technical_subseed(
            root(),
            F4TechnicalDecision::index_hop_bernoulli(
                &event_id,
                &library_id,
                pre_hop_instance("PH_b"),
            ),
        )
        .expect("hop b")
    );
}

#[test]
fn reordering_other_decisions_does_not_change_a_key() {
    let (event_id, library_id, recipient_id) = inputs();
    let target = F4TechnicalDecision::chimera_second_source(
        &event_id,
        &library_id,
        17,
        NonZeroU64::new(43).expect("nonzero breakpoint"),
    );
    let expected = derive_f4_technical_subseed(root(), target).expect("target");
    let decisions = [
        F4TechnicalDecision::fragment_start(&event_id, &library_id, 17),
        F4TechnicalDecision::contamination_percent(&event_id, &library_id, &recipient_id),
        target,
        F4TechnicalDecision::index_hop_destination(
            &event_id,
            &library_id,
            pre_hop_instance("PH_pre-hop-4"),
        ),
    ];
    for decision in decisions.into_iter().rev() {
        let _ = derive_f4_technical_subseed(root(), decision).expect("independent key");
    }
    assert_eq!(
        derive_f4_technical_subseed(root(), target).expect("target after reordering"),
        expected
    );
}

#[test]
fn f4_uniform_closed_extreme_interval_with_interior_numerator_stays_finite() {
    let (event_id, library_id, _) = inputs();
    let decision = F4TechnicalDecision::fragment_insert_length(&event_id, &library_id, 17);
    let subseed = derive_f4_technical_subseed(root(), decision).expect("subseed");
    let mut numerator_rng = subseed.local_rng();
    let numerator = numerator_rng.next_u64() % (1_u64 << 53);
    assert_eq!(numerator, 0x000e_04d0_125b_20c5);
    assert_ne!(numerator, 0);
    assert_ne!(numerator, (1_u64 << 53) - 1);

    let mut rng = subseed.local_rng();
    let minimum = -f64::MAX;
    let maximum = f64::MAX;
    let sampled = f4_uniform_closed(&mut rng, minimum, maximum).expect("valid interval");
    assert!(sampled.is_finite());
    assert!(sampled > minimum);
    assert!(sampled < maximum);
}

#[test]
fn index_hop_instance_requires_the_canonical_ph_category() {
    let accepted = F4IndexHopInstance::parse("PH_pre-hop-4").expect("canonical PH literal");
    assert_eq!(accepted.as_str(), "PH_pre-hop-4");

    for invalid in ["", "PH_", "PH_pre hop", "PH_pre-hóp", "XX_pre-hop-4"] {
        let error = F4IndexHopInstance::parse(invalid).expect_err("must reject invalid instance");
        assert_eq!(error.code(), "GPS-DOMAIN-001");
        assert_eq!(error.context().get("identifier"), Some(&invalid.to_owned()));
        assert_eq!(
            error.context().get("category"),
            Some(&"f4_index_hop_instance".to_owned())
        );
    }
}
