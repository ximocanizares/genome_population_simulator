use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::demultiplexing::{
    Demultiplexer, DemultiplexingCause, DemultiplexingError, DemultiplexingLibrary,
    PostContaminationPair, Udi,
};
use genome_population_simulator::identity::{EventId, LibraryId, PreHopInstanceId};
use genome_population_simulator::randomness::{
    F4IndexHopInstance, F4TechnicalDecision, RootSeed, derive_f4_technical_subseed,
    f4_uniform_index,
};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn root() -> RootSeed {
    RootSeed::parse("root", ROOT).expect("literal root")
}

fn event() -> EventId {
    EventId::parse("demultiplex-main", "test").expect("event")
}

fn library(value: &str) -> LibraryId {
    LibraryId::parse(value, "test").expect("library")
}

fn instance(value: &str) -> PreHopInstanceId {
    PreHopInstanceId::parse(value, "test").expect("pre-hop instance")
}

fn udi(i7: &str, i5: &str) -> Udi {
    Udi::new(i7, i5).expect("valid UDI")
}

fn libraries(order: &[&str]) -> Vec<DemultiplexingLibrary> {
    order
        .iter()
        .map(|name| match *name {
            "A" => DemultiplexingLibrary::new(library("LIB_A"), udi("ACGTACGT", "TGCATGCA")),
            "B" => DemultiplexingLibrary::new(library("LIB_B"), udi("CGTACGTA", "GCATGCAT")),
            "C" => DemultiplexingLibrary::new(library("LIB_C"), udi("GTACGTAC", "CATGCATG")),
            _ => unreachable!("test library"),
        })
        .collect()
}

fn pair(
    name: &str,
    contaminant: bool,
    payload: &'static str,
) -> PostContaminationPair<&'static str> {
    match name {
        "A" => PostContaminationPair::new(
            library("LIB_A"),
            instance("PH_A"),
            udi("ACGTACGT", "TGCATGCA"),
            contaminant,
            payload,
        ),
        "B" => PostContaminationPair::new(
            library("LIB_B"),
            instance("PH_B"),
            udi("CGTACGTA", "GCATGCAT"),
            contaminant,
            payload,
        ),
        "C" => PostContaminationPair::new(
            library("LIB_C"),
            instance("PH_C"),
            udi("GTACGTAC", "CATGCATG"),
            contaminant,
            payload,
        ),
        _ => unreachable!("test pair"),
    }
}

#[test]
fn udi_requires_two_literal_eight_base_acgt_components_and_preserves_i5_orientation() {
    let literal = udi("ACGTACGT", "TGCATGCA");
    assert_eq!(literal.i7(), "ACGTACGT");
    assert_eq!(literal.i5(), "TGCATGCA");
    assert_eq!(literal.literal(), "ACGTACGT+TGCATGCA");

    for (i7, i5) in [
        ("0", "TGCATGCA"),
        ("ACGTACG", "TGCATGCA"),
        ("ACGTACGTA", "TGCATGCA"),
        ("ACGTNCGT", "TGCATGCA"),
        ("acgtacgt", "TGCATGCA"),
        ("ACGTACGT", "0"),
    ] {
        assert!(matches!(
            Udi::new(i7, i5),
            Err(DemultiplexingError::InvalidUdi)
        ));
    }
}

#[test]
fn configuration_defaults_and_validates_index_hopping_rate_and_udi_components() {
    let defaulted = effective_event(technical_toml(None, libraries_toml()));
    assert_eq!(defaulted.index_hopping_rate(), 0.0);
    assert_eq!(defaulted.libraries()[0].i5(), "TGCATGCA");

    let configured = effective_event(technical_toml(Some("1"), libraries_toml()));
    assert_eq!(configured.index_hopping_rate(), 1.0);

    for rate in ["-0.1", "1.1", "nan"] {
        let declared = declared(technical_toml(Some(rate), libraries_toml()));
        assert_eq!(
            build_effective_config(&declared)
                .expect_err("invalid index hopping rate")
                .code(),
            "GPS-CONF-008"
        );
    }
    for invalid_libraries in [
        libraries_toml().replace("i7 = \"CGTACGTA\"", "i7 = \"ACGTACGT\""),
        libraries_toml().replace("i5 = \"GCATGCAT\"", "i5 = \"TGCATGCA\""),
        libraries_toml().replace("i7 = \"ACGTACGT\"", "i7 = \"0\""),
        libraries_toml().replace("i5 = \"TGCATGCA\"", "i5 = \"ACGTN\""),
    ] {
        assert_eq!(
            build_effective_config(&declared(technical_toml(None, invalid_libraries)))
                .expect_err("invalid UDI declaration")
                .code(),
            "GPS-CONF-008"
        );
    }
}

#[test]
fn positive_index_hopping_rate_requires_two_libraries_during_effective_resolution() {
    let source = technical_toml(
        Some("0.1"),
        r#"[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        .to_owned(),
    )
    .replace("individual_count = 2", "individual_count = 1")
    .replace(
        "individual_ids = [\"sample-A\", \"sample-B\"]",
        "individual_ids = [\"sample-A\"]",
    );
    let error = build_effective_config(&declared(source))
        .expect_err("a positive hopping rate needs another declared library");
    assert_eq!(error.code(), "GPS-CONF-008");
    assert_eq!(
        error.context().get("fields"),
        Some(&"technical_event.index_hopping_rate".to_owned())
    );
}

#[test]
fn rate_zero_retains_each_pair_in_its_declared_udi_without_copying_content() {
    let mut demultiplexer =
        Demultiplexer::new(root(), event(), 0.0, libraries(&["A", "B"])).expect("demultiplexer");
    let outcome = demultiplexer
        .demultiplex(vec![
            pair("A", false, "own-A"),
            pair("B", true, "contaminant-B"),
        ])
        .expect("demultiplexing");

    assert_eq!(outcome.report().reassigned_pairs(), 0);
    assert_eq!(outcome.report().discarded_pairs(), 0);
    let a = &outcome.final_libraries()[&library("LIB_A")];
    let b = &outcome.final_libraries()[&library("LIB_B")];
    assert_eq!(a.len(), 1);
    assert_eq!(b.len(), 1);
    assert_eq!(a[0].payload(), &"own-A");
    assert_eq!(b[0].payload(), &"contaminant-B");
    assert!(!a[0].is_contaminant());
    assert!(b[0].is_contaminant());
    assert_eq!(a[0].final_udi().literal(), "ACGTACGT+TGCATGCA");
    assert_eq!(
        a[0].technical_identity().technical_pair_id().as_str(),
        "TP_J5UDQCBKI55ML4J2J6S74GQFO4"
    );
}

#[test]
fn prevalidates_pre_hop_identity_and_declared_observed_udi_at_rate_zero() {
    let invalid_identity = PostContaminationPair::new(
        library("LIB_A"),
        instance("BAD"),
        udi("ACGTACGT", "TGCATGCA"),
        false,
        "invalid-identity",
    );
    let mismatched_declared_udi = PostContaminationPair::new(
        library("LIB_A"),
        instance("PH_MISMATCH"),
        udi("CGTACGTA", "GCATGCAT"),
        false,
        "mismatched-udi",
    );

    let mut demultiplexer =
        Demultiplexer::new(root(), event(), 0.0, libraries(&["A", "B"])).expect("demultiplexer");
    assert!(matches!(
        demultiplexer.demultiplex(vec![invalid_identity]),
        Err(DemultiplexingError::TechnicalIdentity(ref diagnostic))
            if diagnostic.code() == "GPS-DOMAIN-001"
    ));

    let mut demultiplexer =
        Demultiplexer::new(root(), event(), 0.0, libraries(&["A", "B"])).expect("demultiplexer");
    assert!(matches!(
        demultiplexer.demultiplex(vec![mismatched_declared_udi]),
        Err(DemultiplexingError::InvalidUdi)
    ));
}

#[test]
fn rate_one_transfers_each_post_contamination_pair_once_to_another_library() {
    let mut demultiplexer = Demultiplexer::new(root(), event(), 1.0, libraries(&["A", "B", "C"]))
        .expect("demultiplexer");
    let outcome = demultiplexer
        .demultiplex(vec![
            pair("A", false, "own-A"),
            pair("B", true, "contaminant-B"),
            pair("C", false, "own-C"),
        ])
        .expect("demultiplexing");

    assert_eq!(outcome.report().reassigned_pairs(), 3);
    assert_eq!(outcome.report().discarded_pairs(), 0);
    let reassigned = outcome.report().reassigned_by_sample();
    assert_eq!(
        reassigned
            .iter()
            .map(|entry| entry.pair_count())
            .sum::<u64>(),
        3
    );
    assert!(reassigned.iter().all(|entry| {
        entry.cause() == DemultiplexingCause::IndexHopping
            && entry
                .destination_library_id()
                .is_some_and(|destination| destination != entry.source_library_id())
    }));
    let all = outcome
        .final_libraries()
        .values()
        .flat_map(|pairs| pairs.iter())
        .collect::<Vec<_>>();
    assert_eq!(all.len(), 3, "a transfer does not duplicate a pair");
    assert_eq!(
        all.iter()
            .map(|pair| *pair.payload())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from(["contaminant-B", "own-A", "own-C"])
    );
    assert!(all.iter().any(|pair| pair.is_contaminant()));
    assert!(all.iter().all(|pair| {
        (pair.payload() != &"own-A" || pair.final_library_id() != &library("LIB_A"))
            && (pair.payload() != &"contaminant-B" || pair.final_library_id() != &library("LIB_B"))
            && (pair.payload() != &"own-C" || pair.final_library_id() != &library("LIB_C"))
    }));
    assert!(all.iter().all(|pair| {
        pair.final_udi().literal()
            == match pair.final_library_id().as_str() {
                "LIB_A" => "ACGTACGT+TGCATGCA".to_owned(),
                "LIB_B" => "CGTACGTA+GCATGCAT".to_owned(),
                "LIB_C" => "GTACGTAC+CATGCATG".to_owned(),
                _ => unreachable!("declared final library"),
            }
    }));
}

#[test]
fn hopping_destination_matches_the_canonical_candidate_list() {
    let mut demultiplexer = Demultiplexer::new(root(), event(), 1.0, libraries(&["A", "B", "C"]))
        .expect("demultiplexer");
    let outcome = demultiplexer
        .demultiplex(vec![
            pair("A", false, "own-A"),
            pair("B", false, "own-B"),
            pair("C", false, "own-C"),
        ])
        .expect("demultiplexing");

    for (source, pre_hop_instance, payload) in [
        ("LIB_A", "PH_A", "own-A"),
        ("LIB_B", "PH_B", "own-B"),
        ("LIB_C", "PH_C", "own-C"),
    ] {
        let actual = outcome
            .final_libraries()
            .values()
            .flat_map(|pairs| pairs.iter())
            .find(|pair| pair.payload() == &payload)
            .expect("demultiplexed pair");
        assert_eq!(
            actual.final_library_id(),
            &canonical_destination(source, pre_hop_instance),
            "{payload}"
        );
    }
}

#[test]
fn unexpected_udi_is_defensively_discarded_and_only_aggregate_counts_remain() {
    let unexpected = PostContaminationPair::new(
        library("LIB_A"),
        instance("PH_UNEXPECTED"),
        udi("AAAAAAAA", "CCCCCCCC"),
        false,
        "discarded",
    );
    let mut demultiplexer =
        Demultiplexer::new(root(), event(), 0.0, libraries(&["A", "B"])).expect("demultiplexer");
    let outcome = demultiplexer
        .demultiplex(vec![unexpected])
        .expect("defensive discard");
    assert_eq!(outcome.report().discarded_pairs(), 1);
    assert_eq!(outcome.report().unassigned_pairs(), 1);
    assert_eq!(outcome.report().reassigned_pairs(), 0);
    for aggregates in [
        outcome.report().discarded_by_sample(),
        outcome.report().unassigned_by_sample(),
    ] {
        assert_eq!(aggregates.len(), 1);
        let aggregate = &aggregates[0];
        assert_eq!(aggregate.source_library_id(), &library("LIB_A"));
        assert_eq!(aggregate.destination_library_id(), None);
        assert_eq!(aggregate.cause(), DemultiplexingCause::UnrecognizedUdi);
        assert_eq!(aggregate.pair_count(), 1);
    }
    assert!(outcome.final_libraries().values().all(Vec::is_empty));
}

#[test]
fn nonzero_rate_needs_another_library_and_duplicate_components_are_rejected() {
    assert!(matches!(
        Demultiplexer::new(root(), event(), 0.1, libraries(&["A"])),
        Err(DemultiplexingError::MissingAlternativeLibrary)
    ));
    let duplicate_i7 = vec![
        DemultiplexingLibrary::new(library("LIB_A"), udi("ACGTACGT", "TGCATGCA")),
        DemultiplexingLibrary::new(library("LIB_B"), udi("ACGTACGT", "GCATGCAT")),
    ];
    assert!(matches!(
        Demultiplexer::new(root(), event(), 0.0, duplicate_i7),
        Err(DemultiplexingError::DuplicateI7)
    ));
    let duplicate_i5 = vec![
        DemultiplexingLibrary::new(library("LIB_A"), udi("ACGTACGT", "TGCATGCA")),
        DemultiplexingLibrary::new(library("LIB_B"), udi("CGTACGTA", "TGCATGCA")),
    ];
    assert!(matches!(
        Demultiplexer::new(root(), event(), 0.0, duplicate_i5),
        Err(DemultiplexingError::DuplicateI5)
    ));
}

#[test]
fn repeated_final_cause_is_a_technical_id_collision_not_a_silent_rename() {
    let mut demultiplexer =
        Demultiplexer::new(root(), event(), 0.0, libraries(&["A", "B"])).expect("demultiplexer");
    let error = demultiplexer
        .demultiplex(vec![pair("A", false, "first"), pair("A", false, "second")])
        .expect_err("the same final library and PH identity collide");
    assert!(matches!(
        error,
        DemultiplexingError::TechnicalIdentity(ref diagnostic) if diagnostic.code() == "GPS-DOMAIN-002"
    ));
}

#[test]
fn decisions_and_final_ids_are_invariant_to_library_and_input_permutations() {
    let forward = observed_assignments(&["A", "B", "C"], &["A", "B", "C"]);
    let reversed = observed_assignments(&["C", "B", "A"], &["C", "A", "B"]);
    assert_eq!(forward, reversed);
}

#[test]
#[ignore = "manual release performance regression check"]
fn release_scale_hopping_keeps_each_pair_once_with_canonical_destinations() {
    const LIBRARY_COUNT: usize = 64;
    const PAIR_COUNT: usize = 200_000;

    let libraries = (0..LIBRARY_COUNT)
        .map(|index| {
            DemultiplexingLibrary::new(
                library(&format!("LIB_{index:02}")),
                udi(&udi_component(index), &udi_component(index + LIBRARY_COUNT)),
            )
        })
        .collect();
    let pairs = (0..PAIR_COUNT)
        .map(|ordinal| {
            let source = ordinal % LIBRARY_COUNT;
            PostContaminationPair::new(
                library(&format!("LIB_{source:02}")),
                instance(&format!("PH_SCALE_{ordinal}")),
                udi(
                    &udi_component(source),
                    &udi_component(source + LIBRARY_COUNT),
                ),
                false,
                ordinal,
            )
        })
        .collect();
    let mut demultiplexer =
        Demultiplexer::new(root(), event(), 1.0, libraries).expect("demultiplexer");
    let outcome = demultiplexer.demultiplex(pairs).expect("demultiplexing");

    assert_eq!(outcome.report().reassigned_pairs(), PAIR_COUNT as u64);
    assert_eq!(
        outcome
            .final_libraries()
            .values()
            .map(Vec::len)
            .sum::<usize>(),
        PAIR_COUNT
    );
}

fn observed_assignments(
    library_order: &[&str],
    pair_order: &[&str],
) -> BTreeMap<String, BTreeSet<String>> {
    let pairs = pair_order
        .iter()
        .map(|name| {
            pair(
                name,
                *name == "B",
                match *name {
                    "A" => "own-A",
                    "B" => "contaminant-B",
                    "C" => "own-C",
                    _ => unreachable!("test pair"),
                },
            )
        })
        .collect();
    let mut demultiplexer =
        Demultiplexer::new(root(), event(), 1.0, libraries(library_order)).expect("demultiplexer");
    demultiplexer
        .demultiplex(pairs)
        .expect("demultiplexing")
        .final_libraries()
        .iter()
        .map(|(library_id, pairs)| {
            (
                library_id.as_str().to_owned(),
                pairs
                    .iter()
                    .map(|pair| {
                        pair.technical_identity()
                            .technical_pair_id()
                            .as_str()
                            .to_owned()
                    })
                    .collect(),
            )
        })
        .collect()
}

fn udi_component(mut value: usize) -> String {
    let mut bases = [b'A'; 8];
    for base in bases.iter_mut().rev() {
        *base = match value % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        value /= 4;
    }
    String::from_utf8(bases.to_vec()).expect("ASCII bases")
}

fn canonical_destination(source: &str, pre_hop_instance: &str) -> LibraryId {
    let source = library(source);
    let instance = F4IndexHopInstance::parse(pre_hop_instance).expect("pre-hop identity");
    let event = event();
    let decision = F4TechnicalDecision::index_hop_destination(&event, &source, instance);
    let mut rng = derive_f4_technical_subseed(root(), decision)
        .expect("subseed")
        .local_rng();
    let index = f4_uniform_index(&mut rng, 2).expect("candidate index") as usize;
    [library("LIB_A"), library("LIB_B"), library("LIB_C")]
        .into_iter()
        .filter(|candidate| candidate != &source)
        .nth(index)
        .expect("canonical candidate")
}

fn libraries_toml() -> String {
    r#"[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"

[[technical_event.libraries]]
individual_id = "sample-B"
i7 = "CGTACGTA"
i5 = "GCATGCAT"
"#
    .to_owned()
}

fn technical_toml(index_hopping_rate: Option<&str>, libraries: String) -> String {
    let rate = index_hopping_rate
        .map(|rate| format!("index_hopping_rate = {rate}\n"))
        .unwrap_or_default();
    format!(
        r#"config_kind = "general"
schema_version = "2"

[reference]
kind = "simulated"
id = "reference-main"

[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1000

[[initial_populations]]
population_id = "P1"
individual_count = 2
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A", "sample-B"]

[technical_event]
event_id = "demultiplex-main"
individual_ids = ["sample-A", "sample-B"]
target_coverage_wgs = 1
{rate}
{libraries}"#
    )
}

fn declared(text: String) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    parse_declared_config(
        admit_toml_header(
            Path::new("f4-demultiplexing.toml"),
            text.into_bytes(),
            ExpectedConfigKind::General,
        )
        .expect("header"),
    )
    .expect("declared configuration")
}

fn effective_event(
    text: String,
) -> genome_population_simulator::application::EffectiveTechnicalEvent {
    build_effective_config(&declared(text))
        .expect("effective configuration")
        .technical_event()
        .expect("technical event")
        .clone()
}
