use genome_population_simulator::identity::{
    EventId, LibraryId, PreHopInstanceId, TechnicalPairId, TemplateId,
};
use genome_population_simulator::identity::{IndividualId, ReferenceId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::technical_identity::{
    F4TechnicalIdentityDeriver, TechnicalMate, TechnicalPairIdentity, TemplateOrdinal,
    derive_technical_identity_plan,
};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn root() -> RootSeed {
    RootSeed::parse("root", ROOT).expect("literal root")
}

fn event() -> EventId {
    EventId::parse("event-illumina", "test").expect("event")
}

fn library(value: &str) -> LibraryId {
    LibraryId::parse(value, "test").expect("library")
}

fn parse_template(value: &str) -> TemplateId {
    TemplateId::parse(value, "test").expect("template")
}

fn individual(value: &str) -> IndividualId {
    IndividualId::parse(value, "test").expect("individual")
}

fn reference() -> ReferenceId {
    ReferenceId::parse("reference-main", "test").expect("reference")
}

#[test]
fn f4_technical_identity_matches_literal_blake3_vectors() {
    let source_library = library("LIB_SOURCE");
    let recipient_library = library("LIB_RECIPIENT");
    let template = parse_template("TPL_LIB_SOURCE_7");
    let mut deriver = F4TechnicalIdentityDeriver::new(root(), event());

    let own = deriver
        .derive_own_pre_hop_instance(&source_library, &template)
        .expect("own instance");
    let contamination = deriver
        .derive_contamination_pre_hop_instance(&recipient_library, &own)
        .expect("contamination copy");
    let pair = deriver
        .derive_technical_pair(&recipient_library, &contamination)
        .expect("final pair");

    assert_eq!(own.as_str(), "PH_2JEXAKCSWMBIFBQUNT6I2I4LQQ");
    assert_eq!(contamination.as_str(), "PH_FA54UQIJX35IRARTKLJWWGC2RE");
    assert_eq!(pair.as_str(), "TP_PZIA5RTAPW2KPG3VYQW7OWLTLE");
}

#[test]
fn technical_ids_have_fixed_valid_encoding_and_distinct_recipient_copies() {
    let source = library(&format!("L{}", "s".repeat(63)));
    let recipient_a = library(&format!("R{}", "a".repeat(63)));
    let recipient_b = library(&format!("R{}", "b".repeat(63)));
    let template = parse_template(&format!("T{}", "t".repeat(63)));
    let mut deriver = F4TechnicalIdentityDeriver::new(root(), event());
    let own = deriver
        .derive_own_pre_hop_instance(&source, &template)
        .expect("own instance");
    let copy_a = deriver
        .derive_contamination_pre_hop_instance(&recipient_a, &own)
        .expect("first recipient copy");
    let copy_b = deriver
        .derive_contamination_pre_hop_instance(&recipient_b, &own)
        .expect("second recipient copy");
    let pair = deriver
        .derive_technical_pair(&recipient_a, &copy_a)
        .expect("pair");

    for (id, prefix) in [
        (own.as_str(), "PH_"),
        (copy_a.as_str(), "PH_"),
        (copy_b.as_str(), "PH_"),
        (pair.as_str(), "TP_"),
    ] {
        assert!(id.starts_with(prefix));
        assert_eq!(id.len(), prefix.len() + 26);
        assert!(
            id[prefix.len()..]
                .bytes()
                .all(|byte| byte.is_ascii_uppercase() || (b'2'..=b'7').contains(&byte))
        );
    }
    assert_ne!(copy_a, copy_b);
}

#[test]
fn technical_id_derivation_is_repeatable_and_independent_of_call_order() {
    let source = library("LIB_SOURCE");
    let recipient = library("LIB_RECIPIENT");
    let template = parse_template("TPL_LIB_SOURCE_7");
    let mut first = F4TechnicalIdentityDeriver::new(root(), event());
    let first_own = first
        .derive_own_pre_hop_instance(&source, &template)
        .expect("first own");
    let first_copy = first
        .derive_contamination_pre_hop_instance(&recipient, &first_own)
        .expect("first copy");
    let first_pair = first
        .derive_technical_pair(&recipient, &first_copy)
        .expect("first pair");

    let unrelated_library = library("LIB_UNRELATED");
    let unrelated_template = parse_template("TPL_LIB_UNRELATED_9");
    let mut reordered = F4TechnicalIdentityDeriver::new(root(), event());
    let _ = reordered
        .derive_own_pre_hop_instance(&unrelated_library, &unrelated_template)
        .expect("unrelated own");
    let reordered_own = reordered
        .derive_own_pre_hop_instance(&source, &template)
        .expect("reordered own");
    let reordered_copy = reordered
        .derive_contamination_pre_hop_instance(&recipient, &reordered_own)
        .expect("reordered copy");
    let reordered_pair = reordered
        .derive_technical_pair(&recipient, &reordered_copy)
        .expect("reordered pair");

    assert_eq!(first_own, reordered_own);
    assert_eq!(first_copy, reordered_copy);
    assert_eq!(first_pair, reordered_pair);
}

#[test]
fn technical_identity_collisions_fail_for_both_new_classes() {
    let source = library("LIB_SOURCE");
    let recipient = library("LIB_RECIPIENT");
    let template = parse_template("TPL_LIB_SOURCE_7");
    let mut deriver = F4TechnicalIdentityDeriver::new(root(), event());
    let own = deriver
        .derive_own_pre_hop_instance(&source, &template)
        .expect("own instance");
    assert_eq!(
        deriver
            .derive_own_pre_hop_instance(&source, &template)
            .expect_err("same own instance collides")
            .code(),
        "GPS-DOMAIN-002"
    );

    let pair = deriver
        .derive_technical_pair(&recipient, &own)
        .expect("pair");
    assert_eq!(
        deriver
            .derive_technical_pair(&recipient, &own)
            .expect_err("same pair collides")
            .code(),
        "GPS-DOMAIN-002"
    );

    let _: PreHopInstanceId = own;
    let _: TechnicalPairId = pair;
}

#[test]
fn paired_pre_hop_identity_reuses_the_f3_template_without_cross_library_mismatch() {
    let sample_a = individual("sample-A");
    let only_a = derive_technical_identity_plan(
        root(),
        &reference(),
        event(),
        std::slice::from_ref(&sample_a),
    )
    .expect("single-library technical plan");
    let reordered = derive_technical_identity_plan(
        root(),
        &reference(),
        event(),
        &[individual("sample-B"), sample_a.clone()],
    )
    .expect("reordered technical plan");

    let template = only_a
        .derive_template(&sample_a, TemplateOrdinal::new(7))
        .expect("F3 template");
    assert_eq!(
        template.template_id().as_str(),
        "TPL_LIB_EB6VHCIWUACHDO7B5IVJ5RXTCQ_7"
    );
    assert_eq!(
        template.library_id(),
        only_a
            .library(&sample_a)
            .expect("source library")
            .library_id()
    );

    let mut first = F4TechnicalIdentityDeriver::new(root(), event());
    let own = first
        .derive_own_pre_hop_pair(template)
        .expect("own pre-hop pair");
    let reordered_template = reordered
        .derive_template(&sample_a, TemplateOrdinal::new(7))
        .expect("same F3 template after reordering");
    let mut after_reordering = F4TechnicalIdentityDeriver::new(root(), event());
    let reordered_own = after_reordering
        .derive_own_pre_hop_pair(reordered_template)
        .expect("same pre-hop pair after reordering");

    assert_eq!(own, reordered_own);
    assert_eq!(own.mate_roles(), [TechnicalMate::R1, TechnicalMate::R2]);
    assert_eq!(own.mate_roles().map(TechnicalMate::number), [1, 2]);
    assert!(own.pre_hop_instance_id().as_str().starts_with("PH_"));
    assert_eq!(
        own.pre_hop_instance_id().as_str(),
        "PH_DSRV252THT4ES2M33JP2QBTFJQ"
    );
}

#[test]
fn final_pair_identity_always_exposes_both_mates_with_the_same_pair_id() {
    let source = library("LIB_SOURCE");
    let template = parse_template("TPL_LIB_SOURCE_7");
    let mut deriver = F4TechnicalIdentityDeriver::new(root(), event());
    let pre_hop = deriver
        .derive_own_pre_hop_instance(&source, &template)
        .expect("pre-hop instance from the frozen F4.0 vector");
    let pair_id = deriver
        .derive_technical_pair(&source, &pre_hop)
        .expect("F4.0 derivation vector, not FASTQ emission");
    let pair = TechnicalPairIdentity::new(pair_id.clone());
    let mates = pair.mates();

    assert_eq!(pair.technical_pair_id(), &pair_id);
    assert_eq!(
        mates.each_ref().map(|mate| mate.technical_pair_id()),
        [&pair_id, &pair_id]
    );
    assert_eq!(
        mates.each_ref().map(|mate| mate.mate()),
        [TechnicalMate::R1, TechnicalMate::R2]
    );
    assert_ne!(mates[0].mate(), mates[1].mate());
}
