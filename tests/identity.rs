use genome_population_simulator::identity::{
    AutomaticIdRequest, EventId, HumanLabel, IdClass, IdentityRegistry, IndividualId, ReferenceId,
    SequenceId, derive_matrix_variant_id,
};
use genome_population_simulator::randomness::RootSeed;

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn root() -> RootSeed {
    RootSeed::parse("root", ROOT).expect("literal root")
}

#[test]
fn provided_ids_accept_the_exact_ascii_boundary_and_labels_are_unicode() {
    let one = IndividualId::parse("A", "execution").expect("one character");
    let sixty_four =
        IndividualId::parse(format!("A{}", "z".repeat(63)), "execution").expect("64 characters");
    assert_eq!(one.as_str(), "A");
    assert_eq!(sixty_four.as_str().len(), 64);
    assert!(IndividualId::parse(format!("A{}", "z".repeat(64)), "execution").is_err());
    assert!(IndividualId::parse("A B", "execution").is_err());
    assert!(IndividualId::parse("../A", "execution").is_err());
    assert!(IndividualId::parse("Á", "execution").is_err());
    assert!(IndividualId::parse("Sample.UPPER_1-ok", "execution").is_ok());
    assert_eq!(HumanLabel::new("Muestra α").as_str(), "Muestra α");
}

#[test]
fn registry_separates_categories_but_rejects_reuse_in_the_same_scope() {
    let mut registry = IdentityRegistry::default();
    registry
        .register_individual(
            "execution",
            IndividualId::parse("same", "execution").expect("individual"),
        )
        .expect("first registration");
    registry
        .register_event(
            "execution",
            EventId::parse("same", "execution").expect("event"),
        )
        .expect("different category");
    let duplicate = registry.register_individual(
        "execution",
        IndividualId::parse("same", "execution").expect("individual"),
    );
    assert_eq!(duplicate.expect_err("duplicate").code(), "GPS-DOMAIN-002");
}

#[test]
fn all_a6_vectors_match_material_and_identifier() {
    let cases = [
        (
            IdClass::Execution,
            "execution",
            vec![(
                "config_sha256",
                "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
            )],
            "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4",
        ),
        (
            IdClass::Reference,
            "execution",
            vec![
                ("source_kind", "simulated"),
                (
                    "source_digest",
                    "89abcdef0123456789abcdef0123456789abcdef0123456789abcdef01234567",
                ),
            ],
            "REF_J4RT4NORBLKFQ4QGQF2JFO6YVY",
        ),
        (
            IdClass::Individual,
            "population-main",
            vec![("founder_key", "founder-alpha")],
            "IND_FQEFG3NP2UP2GJLSNGJLFOXGUE",
        ),
        (
            IdClass::Event,
            "execution",
            vec![("event_key", "event-illumina")],
            "EVT_2V4A26ZPZMCRUEKKJQIZJ7EEVU",
        ),
    ];

    for (class, scope, causes, expected) in cases {
        let request = AutomaticIdRequest::new(class, scope, &causes).expect("catalog request");
        let mut registry = IdentityRegistry::default();
        let observed = match class {
            IdClass::Execution => registry
                .derive_execution(root(), &request)
                .map(|id| id.to_string()),
            IdClass::Reference => registry
                .derive_reference(root(), &request)
                .map(|id| id.to_string()),
            IdClass::Individual => registry
                .derive_individual(root(), &request)
                .map(|id| id.to_string()),
            IdClass::Event => registry
                .derive_event(root(), &request)
                .map(|id| id.to_string()),
        }
        .expect("derive vector");
        assert_eq!(observed, expected);
        assert_eq!(observed.len(), class_prefix_len(class) + 26);
        assert!(
            observed[class_prefix_len(class)..]
                .bytes()
                .all(|byte| byte.is_ascii_uppercase() || (b'2'..=b'7').contains(&byte))
        );
    }
}

fn class_prefix_len(class: IdClass) -> usize {
    match class {
        IdClass::Execution | IdClass::Reference => 4,
        IdClass::Individual | IdClass::Event => 4,
    }
}

#[test]
fn matrix_variant_id_is_stable_and_excludes_population_and_row_order() {
    let reference = ReferenceId::parse("ref-main", "test").expect("reference");
    let sequence = SequenceId::parse("chr1", "test").expect("sequence");
    let first = derive_matrix_variant_id(&reference, &sequence, 7, b'A', b'C').expect("variant");
    let repeated = derive_matrix_variant_id(&reference, &sequence, 7, b'A', b'C').expect("variant");
    let changed_alt =
        derive_matrix_variant_id(&reference, &sequence, 7, b'A', b'G').expect("variant");
    assert_eq!(first.as_str(), "VAR_ETWKMO4NA6HZ4YXCBOLCN34GTY");
    assert!(include_str!("../catalogs/variant-id-v1.tsv").contains(first.as_str()));
    assert_eq!(first, repeated);
    assert_ne!(first, changed_alt);
    assert!(first.as_str().starts_with("VAR_"));
}

#[test]
fn automatic_collision_is_rejected_and_types_remain_distinct() {
    let request = AutomaticIdRequest::new(
        IdClass::Event,
        "execution",
        &[("event_key", "event-illumina")],
    )
    .expect("request");
    let mut registry = IdentityRegistry::default();
    registry.derive_event(root(), &request).expect("first");
    assert_eq!(
        registry
            .derive_event(root(), &request)
            .expect_err("collision")
            .code(),
        "GPS-DOMAIN-002"
    );

    let reference = ReferenceId::parse("same", "execution").expect("reference");
    let individual = IndividualId::parse("same", "execution").expect("individual");
    assert_eq!(reference.as_str(), individual.as_str());
}
