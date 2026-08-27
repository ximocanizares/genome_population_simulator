use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId};
use genome_population_simulator::randomness::{RootSeed, RuntimeSubseedScope, derive_subseed};
use genome_population_simulator::technical_identity::{
    TemplateOrdinal, derive_technical_identity_plan,
};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn root() -> RootSeed {
    RootSeed::parse("root", ROOT).expect("literal root")
}

fn event() -> EventId {
    EventId::parse("event-illumina", "test").expect("event")
}

fn reference() -> ReferenceId {
    ReferenceId::parse("reference-main", "test").expect("reference")
}

fn individual(id: &str) -> IndividualId {
    IndividualId::parse(id, "test").expect("individual")
}

fn plan(
    individuals: &[&str],
) -> genome_population_simulator::technical_identity::TechnicalIdentityPlan {
    let individual_ids = individuals
        .iter()
        .map(|id| individual(id))
        .collect::<Vec<_>>();
    derive_technical_identity_plan(root(), &reference(), event(), &individual_ids)
        .expect("technical plan")
}

#[test]
fn technical_identities_and_template_subseed_are_repeatable() {
    let sample_a = individual("sample-A");
    let first = plan(&["sample-A"]);
    let repeated = plan(&["sample-A"]);

    assert_eq!(first.event_id().as_str(), repeated.event_id().as_str());
    assert_eq!(first.batch_id().as_str(), repeated.batch_id().as_str());
    assert_eq!(first.event_subseed(), repeated.event_subseed());
    assert_eq!(
        first.library(&sample_a).expect("library"),
        repeated.library(&sample_a).expect("library")
    );
    let first_template = first
        .derive_template(&sample_a, TemplateOrdinal::new(7))
        .expect("template");
    assert_eq!(first.batch_id().as_str(), "BAT_LVFBOBZHZAZORHNNGRATZ55AFE");
    assert_eq!(
        first
            .library(&sample_a)
            .expect("library")
            .library_id()
            .as_str(),
        "LIB_EB6VHCIWUACHDO7B5IVJ5RXTCQ"
    );
    assert_eq!(
        first
            .library(&sample_a)
            .expect("library")
            .read_group_id()
            .as_str(),
        "RG_FGR44MPG223B52YE6ROSTNX3B4"
    );
    assert_eq!(
        first_template.template_id().as_str(),
        "TPL_LIB_EB6VHCIWUACHDO7B5IVJ5RXTCQ_7"
    );
    assert_eq!(
        first_template.subseed().to_hex(),
        "dcbd175161ab253763525c07e556de419087647a29ff6d8fac39077cb12d1ccc"
    );
    let repeated_template = repeated
        .derive_template(&sample_a, TemplateOrdinal::new(7))
        .expect("template");
    assert_eq!(first_template, repeated_template);
}

#[test]
fn existing_template_is_invariant_to_individual_addition_and_reordering() {
    let sample_a = individual("sample-A");
    let only_a = plan(&["sample-A"]);
    let reordered_and_extended = plan(&["sample-B", "sample-A"]);

    assert_eq!(
        only_a.event_subseed(),
        reordered_and_extended.event_subseed()
    );
    assert_eq!(only_a.batch_id(), reordered_and_extended.batch_id());
    assert_eq!(
        only_a.library(&sample_a).expect("library"),
        reordered_and_extended.library(&sample_a).expect("library")
    );
    assert_eq!(
        only_a
            .derive_template(&sample_a, TemplateOrdinal::new(3))
            .expect("template"),
        reordered_and_extended
            .derive_template(&sample_a, TemplateOrdinal::new(3))
            .expect("template")
    );
}

#[test]
fn technical_identifiers_are_valid_and_template_ids_are_unique_by_construction() {
    let sample_a = individual("sample-A");
    let plan = plan(&["sample-A"]);
    for id in [
        plan.batch_id().as_str(),
        plan.library(&sample_a)
            .expect("library")
            .library_id()
            .as_str(),
        plan.library(&sample_a)
            .expect("library")
            .read_group_id()
            .as_str(),
    ] {
        assert!(id.len() <= 64);
        assert!(
            id.bytes()
                .all(|byte| { byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-') })
        );
    }
    let template = plan
        .derive_template(&sample_a, TemplateOrdinal::new(0))
        .expect("template");
    assert_eq!(
        template.template_id().as_str(),
        "TPL_LIB_EB6VHCIWUACHDO7B5IVJ5RXTCQ_0"
    );
    assert!(template.template_id().as_str().len() <= 64);
    assert!(
        template
            .template_id()
            .as_str()
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
    );
    assert_eq!(
        plan.derive_template(&sample_a, TemplateOrdinal::new(0))
            .expect("same canonical causes repeat deterministically"),
        template
    );
    assert_ne!(
        template.template_id(),
        plan.derive_template(&sample_a, TemplateOrdinal::new(1))
            .expect("next ordinal")
            .template_id()
    );
    let last_ordinal = plan
        .derive_template(&sample_a, TemplateOrdinal::new(u64::MAX))
        .expect("maximum u64 ordinal");
    assert_eq!(
        last_ordinal.template_id().as_str(),
        "TPL_LIB_EB6VHCIWUACHDO7B5IVJ5RXTCQ_18446744073709551615"
    );
    assert_eq!(last_ordinal.template_id().as_str().len(), 55);
}

#[test]
fn template_derivation_keeps_no_per_template_registry() {
    let sample_a = individual("sample-A");
    let plan = plan(&["sample-A"]);

    // `derive_template` only borrows the plan. Iterar muchos ordinales prueba
    // el contrato de streaming: no existe estado mutable donde retener IDs.
    for ordinal in 0..100_000 {
        let template = plan
            .derive_template(&sample_a, TemplateOrdinal::new(ordinal))
            .expect("known library");
        assert!(template.template_id().as_str().starts_with("TPL_LIB_"));
    }
    assert_eq!(
        plan.derive_template(&sample_a, TemplateOrdinal::new(99_999))
            .expect("repeated ordinal")
            .template_id()
            .as_str(),
        "TPL_LIB_EB6VHCIWUACHDO7B5IVJ5RXTCQ_99999"
    );
}

#[test]
fn technical_event_subseed_matches_historical_vector_and_isolates_causes() {
    let plan = plan(&["sample-A"]);
    assert_eq!(
        plan.event_subseed().to_hex(),
        "b6a891fdb16514ffa614eb9731a702c6de5e5c205a704cddf890f4ff94a7f6f3"
    );
    assert_eq!(
        plan.event_subseed(),
        derive_subseed(
            root(),
            RuntimeSubseedScope::TechnicalEvent {
                event_id: "event-illumina",
                reference_id: "reference-main",
            },
        )
        .expect("historical technical_event vector")
    );
    let changed_event = derive_technical_identity_plan(
        root(),
        &reference(),
        EventId::parse("event-other", "test").expect("event"),
        &[individual("sample-A")],
    )
    .expect("changed event plan");
    let changed_reference = derive_technical_identity_plan(
        root(),
        &ReferenceId::parse("reference-other", "test").expect("reference"),
        event(),
        &[individual("sample-A")],
    )
    .expect("changed reference plan");
    assert_ne!(plan.event_subseed(), changed_event.event_subseed());
    assert_ne!(plan.event_subseed(), changed_reference.event_subseed());
}

#[test]
fn invalid_reference_text_is_rejected_before_the_typed_frontier() {
    let error = ReferenceId::parse(
        "references/reference-main.fa",
        "technical_event.reference_id",
    )
    .expect_err("paths are not reference IDs");
    assert_eq!(error.code(), "GPS-DOMAIN-001");
}
