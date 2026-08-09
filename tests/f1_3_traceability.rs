use std::collections::BTreeSet;

use genome_population_simulator::randomness::{
    BinaryIdentityEvidence, EntropySource, GuaranteeStatus, ReproducibilityContext,
    ReproducibilityMetadata, ReproducibilityObservation, RootSeed, SeedOrigin, SubseedRequest,
    compare_binary, compare_biological, compare_semantic, compare_technical, derive_subseed,
    resolve_root_seed,
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

fn event_subseeds(
    root_seed: RootSeed,
    event: &str,
) -> Vec<genome_population_simulator::randomness::EventSubseed> {
    ReproducibilityMetadata::new(
        root_seed,
        SeedOrigin::Declared,
        "reference-main",
        &[event.to_owned()],
    )
    .unwrap()
    .event_subseeds()
    .to_vec()
}

fn observation(root_seed: RootSeed) -> ReproducibilityObservation {
    ReproducibilityObservation {
        context: ReproducibilityContext {
            simulator_version: "0.1.0".to_owned(),
            platform: "x86_64-unknown-linux-gnu".to_owned(),
            root_seed,
            effective_config: "effective-config-sha256".to_owned(),
            input_identity: "input-sha256".to_owned(),
            profile: "portable".to_owned(),
            random_algorithm: "ChaCha12+subseed-v1".to_owned(),
            derivation_version: 1,
            catalog_version: "1".to_owned(),
            dependency_versions: vec!["blake3=1.8.6".to_owned(), "rand_chacha=0.10.0".to_owned()],
            technical_dependency_versions: vec!["read-simulator=1".to_owned()],
            external_tools: vec!["samtools=1.22".to_owned()],
        },
        biological_decisions: Some("parents=2,7;mutations=11,43".to_owned()),
        technical_content: Some("reads-and-qualities-logical-content".to_owned()),
        semantic_content: Some("records:A,C,G,T".to_owned()),
        binary_content: Some(vec![1, 2, 3]),
        event_subseeds: event_subseeds(root_seed, "event-illumina"),
        binary_identity_evidence: None,
    }
}

fn biological_decisions(root_seed: RootSeed) -> [u64; 4] {
    let request = SubseedRequest::from_scope(
        "biological_process",
        &[
            ("generation", "1"),
            ("population_id", "population-main"),
            ("process", "reproduction"),
        ],
    )
    .unwrap();
    let mut rng = derive_subseed(root_seed, &request).unwrap().local_rng();
    [
        rng.next_u64(),
        rng.next_u64(),
        rng.next_u64(),
        rng.next_u64(),
    ]
}

#[test]
fn f1_3_matrix_audits_ten_requirements_with_one_pass_and_nine_partial_reassignments() {
    let matrix = include_str!("../evidence/f1.3/test-req-matrix.tsv");
    let source = include_str!("f1_3_traceability.rs");
    let roadmap = include_str!("../design/hoja-de-ruta-0.1.md");
    let traceability = include_str!("../design/trazabilidad-requisitos-0.1.md");
    let rows = matrix
        .lines()
        .skip(1)
        .map(|row| row.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 10);
    assert_eq!(
        rows.iter().map(|row| row[1]).collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "REQ-CONF-26",
            "REQ-CONF-27",
            "REQ-CONF-111",
            "REQ-CONF-112",
            "REQ-CONF-113",
            "REQ-CONF-114",
            "REQ-CONF-115",
            "REQ-REP-22",
            "REQ-REP-23",
            "REQ-VAL-44",
        ])
    );

    for row in rows {
        assert_eq!(row.len(), 6);
        let symbol = row[4]
            .split_once("::")
            .expect("evidencia Rust con símbolo")
            .1;
        assert!(source.contains(&format!("fn {symbol}(")));
        assert!(traceability.contains(&format!("| {} |", row[1])));
        if row[1] == "REQ-REP-23" {
            assert_eq!((row[2], row[3]), ("pass", "F1.3"));
            assert!(!row[5].starts_with("evidencia parcial:"));
        } else {
            assert_eq!(row[2], "reassigned");
            assert!(matches!(row[3], "F1.5" | "F4.6" | "F9.7"));
            assert!(roadmap.contains(&format!("**{}.", row[3])));
            assert!(row[5].starts_with("evidencia parcial:"));
            assert!(symbol.starts_with("partial_test_req_"));
        }
    }
}

#[test]
fn rng_dependency_derivation_and_catalog_versions_are_order_independent_biological_conditions() {
    let baseline = observation(root());
    let mut reordered = baseline.clone();
    reordered.context.dependency_versions.reverse();
    assert_eq!(
        compare_biological(&baseline, &reordered).status(),
        GuaranteeStatus::Guaranteed
    );

    let mut changed_rand_chacha = baseline.clone();
    changed_rand_chacha.context.dependency_versions =
        vec!["blake3=1.8.6".to_owned(), "rand_chacha=changed".to_owned()];
    let mut changed_blake3 = baseline.clone();
    changed_blake3.context.dependency_versions =
        vec!["blake3=changed".to_owned(), "rand_chacha=0.10.0".to_owned()];
    let mut changed_derivation = baseline.clone();
    changed_derivation.context.derivation_version = 2;
    let mut changed_catalog = baseline.clone();
    changed_catalog.context.catalog_version = "2".to_owned();

    for changed in [
        changed_rand_chacha,
        changed_blake3,
        changed_derivation,
        changed_catalog,
    ] {
        assert_eq!(
            compare_biological(&baseline, &changed).status(),
            GuaranteeStatus::NotGuaranteed
        );
    }
    for condition in [
        "dependency_versions",
        "derivation_version",
        "catalog_version",
    ] {
        assert!(
            compare_biological(&baseline, &baseline)
                .conditions()
                .contains(&condition)
        );
    }
}

#[test]
fn partial_test_req_conf_26_root_is_declared_or_generated_before_test_decisions() {
    let declared = resolve_root_seed(
        Some("000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"),
        &mut FixedEntropy([9; 32]),
    )
    .unwrap();
    let generated = resolve_root_seed(None, &mut FixedEntropy([9; 32])).unwrap();
    assert_eq!(declared.1, SeedOrigin::Declared);
    assert_eq!(
        generated,
        (RootSeed::from_bytes([9; 32]), SeedOrigin::SystemCsprng)
    );
    assert_ne!(
        biological_decisions(declared.0),
        biological_decisions(generated.0)
    );
}

#[test]
fn partial_test_req_conf_27_repeats_an_isolated_rng_decision_vector() {
    let first = biological_decisions(root());
    let second = biological_decisions(root());
    assert_eq!(first, second);
    assert_eq!(
        first,
        [
            5_753_905_640_567_314_474,
            14_586_710_614_320_764_305,
            8_589_965_884_717_974_334,
            3_243_491_804_722_205_367,
        ]
    );
    assert_eq!(
        compare_biological(&observation(root()), &observation(root())).status(),
        GuaranteeStatus::Guaranteed
    );

    let baseline = observation(root());
    let mut changed = baseline.clone();
    changed.context.root_seed = RootSeed::from_bytes([99; 32]);
    assert_eq!(
        compare_biological(&baseline, &changed).status(),
        GuaranteeStatus::NotGuaranteed
    );
    changed = baseline.clone();
    changed.context.effective_config = "changed-config".to_owned();
    assert_eq!(
        compare_biological(&baseline, &changed).status(),
        GuaranteeStatus::NotGuaranteed
    );
    changed = baseline.clone();
    changed.context.input_identity = "changed-input".to_owned();
    assert_eq!(
        compare_biological(&baseline, &changed).status(),
        GuaranteeStatus::NotGuaranteed
    );
    changed = baseline.clone();
    changed.context.profile = "changed-profile".to_owned();
    assert_eq!(
        compare_biological(&baseline, &changed).status(),
        GuaranteeStatus::NotGuaranteed
    );
    changed = baseline.clone();
    changed.context.platform = "unverified-platform".to_owned();
    assert_eq!(
        compare_biological(&baseline, &changed).status(),
        GuaranteeStatus::NotGuaranteed
    );
    changed = baseline.clone();
    changed.context.random_algorithm = "changed-algorithm".to_owned();
    assert_eq!(
        compare_biological(&baseline, &changed).status(),
        GuaranteeStatus::NotGuaranteed
    );
}

#[test]
fn partial_test_req_conf_111_preparatory_comparators_do_not_imply_each_other() {
    let first = observation(root());
    let mut second = first.clone();
    second.binary_content = Some(vec![9, 9, 9]);
    assert_eq!(
        compare_biological(&first, &second).status(),
        GuaranteeStatus::Guaranteed
    );
    assert_eq!(
        compare_technical(&first, &second).status(),
        GuaranteeStatus::Guaranteed
    );
    assert_eq!(
        compare_semantic(&first, &second).status(),
        GuaranteeStatus::Guaranteed
    );
    assert_eq!(
        compare_binary(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
}

#[test]
fn partial_test_req_conf_112_rng_stream_isolated_from_technical_consumption() {
    let expected = biological_decisions(root());
    let mut technical = derive_subseed(
        root(),
        &SubseedRequest::from_scope(
            "technical_event",
            &[
                ("event_id", "new-observation"),
                ("reference_id", "reference-main"),
            ],
        )
        .unwrap(),
    )
    .unwrap()
    .local_rng();
    for _ in 0..100 {
        let _ = technical.next_u64();
    }
    assert_eq!(biological_decisions(root()), expected);
    assert_eq!(
        compare_biological(&observation(root()), &observation(root())).status(),
        GuaranteeStatus::Guaranteed
    );
}

#[test]
fn partial_test_req_conf_113_comparator_requires_subseed_and_technical_dependencies() {
    let first = observation(root());
    let mut second = first.clone();
    assert_eq!(
        compare_technical(&first, &second).status(),
        GuaranteeStatus::Guaranteed
    );
    second.context.technical_dependency_versions =
        vec!["dependency-b=2".to_owned(), "read-simulator=1".to_owned()];
    let mut two_dependencies = first.clone();
    two_dependencies.context.technical_dependency_versions =
        vec!["read-simulator=1".to_owned(), "dependency-b=2".to_owned()];
    assert_eq!(
        compare_technical(&two_dependencies, &second).status(),
        GuaranteeStatus::Guaranteed
    );
    let events = ReproducibilityMetadata::new(
        root(),
        SeedOrigin::Declared,
        "reference-main",
        &["event-a".to_owned(), "event-b".to_owned()],
    )
    .unwrap()
    .event_subseeds()
    .to_vec();
    let mut events_first = first.clone();
    events_first.event_subseeds = events.clone();
    let mut events_reordered = first.clone();
    events_reordered.event_subseeds = events.into_iter().rev().collect();
    assert_eq!(
        compare_technical(&events_first, &events_reordered).status(),
        GuaranteeStatus::Guaranteed
    );
    second.event_subseeds = event_subseeds(root(), "event-nanopore");
    assert_eq!(
        compare_technical(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
    second = first.clone();
    second.context.technical_dependency_versions = vec!["read-simulator=changed".to_owned()];
    assert_eq!(
        compare_technical(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
}

#[test]
fn partial_test_req_conf_114_comparator_separates_logical_content_from_bytes() {
    let first = observation(root());
    let mut second = first.clone();
    second.binary_content = Some(vec![7, 8, 9]);
    assert_eq!(
        compare_semantic(&first, &second).status(),
        GuaranteeStatus::Guaranteed
    );
    let mut tools_first = first.clone();
    tools_first.context.external_tools =
        vec!["samtools=1.22".to_owned(), "bcftools=1.22".to_owned()];
    let mut tools_reordered = first.clone();
    tools_reordered.context.external_tools =
        vec!["bcftools=1.22".to_owned(), "samtools=1.22".to_owned()];
    assert_eq!(
        compare_semantic(&tools_first, &tools_reordered).status(),
        GuaranteeStatus::Guaranteed
    );
    second.context.external_tools = vec!["samtools=changed".to_owned()];
    assert_eq!(
        compare_semantic(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
    assert!(
        compare_semantic(&first, &first)
            .conditions()
            .contains(&"external_tools")
    );
}

#[test]
fn partial_test_req_conf_115_preparatory_binary_comparison_requires_artifact_evidence() {
    let first = observation(root());
    let mut second = first.clone();
    assert!(BinaryIdentityEvidence::new("", "policy", "evidence").is_none());
    assert_eq!(
        compare_binary(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
    second.binary_identity_evidence = BinaryIdentityEvidence::new(
        "effective-config-json-v1",
        "canonical-json-bytes-v1",
        "TEST-EVIDENCE-001",
    );
    let mut verified = first.clone();
    verified.binary_identity_evidence = second.binary_identity_evidence.clone();
    let evidence = verified
        .binary_identity_evidence
        .as_ref()
        .expect("evidencia binaria concreta");
    assert_eq!(evidence.artifact_class(), "effective-config-json-v1");
    assert_eq!(evidence.verification_policy(), "canonical-json-bytes-v1");
    assert_eq!(evidence.evidence_id(), "TEST-EVIDENCE-001");
    assert_eq!(
        compare_binary(&verified, &second).status(),
        GuaranteeStatus::Guaranteed
    );
    second.binary_identity_evidence = BinaryIdentityEvidence::new(
        "manifest-json-v1",
        "canonical-json-bytes-v1",
        "TEST-EVIDENCE-001",
    );
    assert_eq!(
        compare_binary(&verified, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
    second.binary_identity_evidence = verified.binary_identity_evidence.clone();
    second.binary_content = Some(vec![3, 2, 1]);
    assert_eq!(
        compare_binary(&verified, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
}

#[test]
fn partial_test_req_rep_22_rng_metadata_reconstructs_an_event_subseed() {
    let metadata = ReproducibilityMetadata::new(
        root(),
        SeedOrigin::Declared,
        "reference-main",
        &["event-illumina".to_owned()],
    )
    .unwrap();
    let request = SubseedRequest::from_scope(
        "technical_event",
        &[
            ("event_id", "event-illumina"),
            ("reference_id", "reference-main"),
        ],
    )
    .unwrap();
    assert_eq!(
        derive_subseed(metadata.root_seed(), &request).unwrap(),
        metadata.event_subseeds()[0].subseed()
    );
    assert_eq!(metadata.random_algorithm(), "ChaCha12");
    assert_eq!(metadata.derivation_version(), 1);
    assert_eq!(metadata.catalog_version(), "1");

    let first = observation(root());
    let mut changed_input = first.clone();
    changed_input.context.input_identity = "other-input".to_owned();
    assert_eq!(
        compare_biological(&first, &changed_input).status(),
        GuaranteeStatus::NotGuaranteed
    );
    let mut changed_rng_dependency = first.clone();
    changed_rng_dependency.context.dependency_versions = vec!["rand_chacha=changed".to_owned()];
    assert_eq!(
        compare_biological(&first, &changed_rng_dependency).status(),
        GuaranteeStatus::NotGuaranteed
    );
}

#[test]
fn test_req_rep_23_version_change_removes_every_reproduction_guarantee() {
    let first = observation(root());
    let mut second = first.clone();
    second.context.simulator_version = "0.2.0".to_owned();
    assert_eq!(
        compare_biological(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
    assert_eq!(
        compare_technical(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
    assert_eq!(
        compare_semantic(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
    assert_eq!(
        compare_binary(&first, &second).status(),
        GuaranteeStatus::NotGuaranteed
    );
}

#[test]
fn partial_test_req_val_44_preparatory_assessments_expose_status_and_conditions() {
    let first = observation(root());
    let mut missing = first.clone();
    missing.technical_content = None;
    assert_eq!(
        compare_technical(&first, &missing).status(),
        GuaranteeStatus::NotGuaranteed
    );
    let mut absent = first.clone();
    absent.technical_content = None;
    assert_eq!(
        compare_technical(&absent, &absent).status(),
        GuaranteeStatus::NotApplicable
    );
    assert!(
        compare_biological(&first, &first)
            .conditions()
            .contains(&"biological_decisions")
    );
    assert!(
        compare_technical(&first, &first)
            .conditions()
            .contains(&"event_subseeds")
    );
    assert!(
        compare_semantic(&first, &first)
            .conditions()
            .contains(&"semantic_content")
    );
    assert!(
        compare_binary(&first, &first)
            .conditions()
            .contains(&"binary_verification_evidence")
    );
    assert_eq!(
        compare_binary(&first, &first).status(),
        GuaranteeStatus::NotGuaranteed
    );
}
