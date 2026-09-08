//! Oráculos de F4.4.1: el plan conserva claves y recuentos, nunca lecturas.

use std::collections::BTreeSet;
use std::fs;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::biological_contamination::{
    BiologicalContaminationError, BiologicalContaminationMaterializer,
    BiologicalContaminationSettings, CleanLibraryPairs, CleanPairKey,
    plan_biological_contamination,
};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::{
    EventId, IndividualId, LibraryId, PopulationId, PreHopInstanceId,
};
use genome_population_simulator::randomness::{
    F4TechnicalDecision, LocalRng, RootSeed, derive_f4_technical_subseed, f4_uniform_closed,
};

const CLOSED_53_MAX: u64 = (1_u64 << 53) - 1;

// Oráculo deliberadamente literal e independiente del planificador: no llama
// a ninguna primitiva F4 de producción para el índice, porcentaje o permuta.
fn literal_uniform_index(rng: &mut LocalRng, n: u64) -> u64 {
    assert!(n > 0);
    if n == 1 {
        return 0;
    }
    let remainder = (u64::MAX % n + 1) % n;
    let limit = (remainder != 0).then(|| u64::MAX - (remainder - 1));
    loop {
        let word = rng.next_u64();
        if limit.is_none_or(|bound| word < bound) {
            return word % n;
        }
    }
}

fn literal_closed_percent(minimum: f64, maximum: f64, numerator: u64) -> f64 {
    assert!((0.0..=100.0).contains(&minimum));
    assert!((0.0..=100.0).contains(&maximum));
    assert!(minimum <= maximum);
    assert!(numerator <= CLOSED_53_MAX);
    if minimum == maximum || numerator == 0 {
        return minimum;
    }
    if numerator == CLOSED_53_MAX {
        return maximum;
    }
    let fraction = numerator as f64 / CLOSED_53_MAX as f64;
    (1.0 - fraction) * minimum + fraction * maximum
}

fn literal_partial_permutation<T>(rng: &mut LocalRng, values: &mut [T], count: usize) {
    assert!(count <= values.len());
    for index in 0..count {
        let offset = literal_uniform_index(rng, (values.len() - index) as u64) as usize;
        values.swap(index, index + offset);
    }
}

fn id(value: &str) -> IndividualId {
    IndividualId::parse(value, "test").expect("individual")
}
fn population(value: &str) -> PopulationId {
    PopulationId::parse(value, "test").expect("population")
}
fn library(value: &str) -> LibraryId {
    LibraryId::parse(value, "test").expect("library")
}
fn instance(value: &str) -> PreHopInstanceId {
    PreHopInstanceId::parse(value, "test").expect("instance")
}
fn root(last: u8) -> RootSeed {
    let mut bytes = [0_u8; 32];
    bytes[31] = last;
    RootSeed::from_bytes(bytes)
}

/// Instancia con la forma que produce una derivación: el prefijo y 26 símbolos
/// Base32, con los dos bits de relleno a cero.
///
/// El inventario solo admite esa forma, porque en una ejecución una instancia
/// `PH_` siempre se deriva. Persona y ordinal se codifican en `A`-`Z`, cuyo
/// orden ASCII coincide con el numérico, de modo que las claves conservan el
/// mismo orden relativo que cuando estas pruebas usaban nombres legibles.
fn derived_instance(person: &str, ordinal: u8) -> PreHopInstanceId {
    let mut symbols = [b'A'; 26];
    symbols[0] = person.as_bytes()[0];
    symbols[1] = b'A' + ordinal;
    let symbols = std::str::from_utf8(&symbols).expect("ASCII symbols");
    instance(&format!("PH_{symbols}"))
}

fn clean(person: &str, pop: &str, count: u8) -> CleanLibraryPairs {
    let individual = id(person);
    let pairs = (0..count)
        .map(|ordinal| {
            CleanPairKey::new(individual.clone(), &derived_instance(person, ordinal))
                .expect("a derived instance is accepted")
        })
        .collect();
    CleanLibraryPairs::new(
        individual,
        population(pop),
        library(&format!("LIB_{person}")),
        pairs,
    )
    .expect("clean inventory")
}

#[test]
fn percentage_is_closed_reproducible_and_rounds_clean_pairs_nearest_ties_up() {
    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![clean("A", "P1", 5), clean("B", "P1", 5)];
    let settings = BiologicalContaminationSettings::new(50.0, 50.0, None).expect("settings");
    let first =
        plan_biological_contamination(root(7), &event, &settings, &libraries).expect("plan");
    let second = plan_biological_contamination(root(7), &event, &settings, &libraries)
        .expect("repeated plan");
    assert_eq!(first, second);
    assert!(first.recipients().iter().all(|plan| {
        plan.sampled_percent() == 50.0
            && plan.clean_pair_count() == 5
            && plan.target_pair_count() == 3
            && plan.assignments().len() == 3
    }));
    // Un intervalo degenerado no consume el RNG: la semilla queda registrada,
    // pero el resultado se obtiene sin extraer ninguna palabra.
    let recipient = &first.recipients()[0];
    let mut manual_rng = recipient.percent_subseed().local_rng();
    assert_eq!(
        f4_uniform_closed(&mut manual_rng, 50.0, 50.0).expect("closed fixed"),
        recipient.sampled_percent()
    );
    assert_eq!(
        manual_rng.next_u64(),
        recipient.percent_subseed().local_rng().next_u64()
    );

    let variable = BiologicalContaminationSettings::new(12.0, 87.0, Some(vec![id("A")]))
        .expect("variable settings");
    let plan = plan_biological_contamination(root(9), &event, &variable, &libraries).expect("plan");
    let recipient = &plan.recipients()[0];
    let expected_seed = derive_f4_technical_subseed(
        root(9),
        F4TechnicalDecision::contamination_percent(&event, &library("LIB_A"), &id("A")),
    )
    .expect("percent seed");
    let expected =
        f4_uniform_closed(&mut expected_seed.local_rng(), 12.0, 87.0).expect("manual percent");
    assert_eq!(recipient.sampled_percent().to_bits(), expected.to_bits());
    assert!((12.0..=87.0).contains(&expected));
}

#[test]
fn literal_oracle_covers_closed_extremes_and_the_assignment_permutation() {
    assert_eq!(literal_closed_percent(12.0, 87.0, 0), 12.0);
    assert_eq!(literal_closed_percent(12.0, 87.0, CLOSED_53_MAX), 87.0);
    assert_eq!(literal_closed_percent(33.0, 33.0, 0), 33.0);

    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![
        clean("A", "P1", 4),
        clean("B", "P1", 4),
        clean("C", "P1", 4),
    ];
    let settings =
        BiologicalContaminationSettings::new(50.0, 50.0, Some(vec![id("A")])).expect("settings");
    let plan =
        plan_biological_contamination(root(31), &event, &settings, &libraries).expect("plan");
    let recipient = plan.recipients().first().expect("one recipient");

    let mut donors = libraries[1]
        .own_pairs()
        .iter()
        .chain(libraries[2].own_pairs())
        .cloned()
        .collect::<Vec<_>>();
    donors.sort();
    let mut removed = libraries[0].own_pairs().to_vec();
    let mut donor_rng = recipient.donor_subseed().local_rng();
    let mut removal_rng = recipient.removal_subseed().local_rng();
    literal_partial_permutation(&mut donor_rng, &mut donors, 2);
    literal_partial_permutation(&mut removal_rng, &mut removed, 2);
    let expected = donors
        .into_iter()
        .take(2)
        .zip(removed.into_iter().take(2))
        .collect::<Vec<_>>();
    let actual = recipient
        .assignments()
        .iter()
        .map(|assignment| {
            (
                assignment.source_pair().clone(),
                assignment.removed_recipient_pair().clone(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);
}

#[test]
fn plan_uses_only_same_population_clean_pairs_without_replacement_or_chains() {
    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![
        clean("A", "P1", 4),
        clean("B", "P1", 4),
        clean("C", "P1", 4),
        clean("D", "P2", 8),
    ];
    let settings =
        BiologicalContaminationSettings::new(50.0, 50.0, Some(vec![id("A"), id("B"), id("C")]))
            .expect("settings");
    let plan =
        plan_biological_contamination(root(21), &event, &settings, &libraries).expect("plan");
    assert_eq!(
        plan.recipients().len(),
        3,
        "the explicit subset excludes P2"
    );
    for recipient in plan.recipients() {
        let sources = recipient
            .assignments()
            .iter()
            .map(|item| item.source_pair())
            .collect::<BTreeSet<_>>();
        let removed = recipient
            .assignments()
            .iter()
            .map(|item| item.removed_recipient_pair())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            sources.len(),
            recipient.assignments().len(),
            "one source cannot repeat inside a recipient"
        );
        assert_eq!(
            removed.len(),
            recipient.assignments().len(),
            "one own pair is removed once"
        );
        for assignment in recipient.assignments() {
            assert_ne!(
                assignment.source_pair().source_individual_id(),
                recipient.recipient_id()
            );
            assert_eq!(
                assignment.removed_recipient_pair().source_individual_id(),
                recipient.recipient_id()
            );
            assert_ne!(
                assignment.source_pair().source_individual_id().as_str(),
                "D"
            );
        }
    }
}

#[test]
fn preflight_rejects_missing_donor_and_insufficient_capacity_before_any_materialization() {
    let event = EventId::parse("contamination-main", "test").expect("event");
    let alone = vec![clean("A", "P1", 1)];
    let positive = BiologicalContaminationSettings::new(100.0, 100.0, None).expect("settings");
    assert_eq!(
        plan_biological_contamination(root(1), &event, &positive, &alone),
        Err(BiologicalContaminationError::MissingDonor { recipient: id("A") }),
    );
    let insufficient = vec![clean("A", "P1", 3), clean("B", "P1", 2)];
    assert_eq!(
        plan_biological_contamination(root(1), &event, &positive, &insufficient),
        Err(BiologicalContaminationError::InsufficientDonorCapacity {
            recipient: id("A"),
            requested: 3,
            available: 2,
        }),
    );
    let zero = BiologicalContaminationSettings::disabled();
    let plan =
        plan_biological_contamination(root(1), &event, &zero, &alone).expect("0% needs no donor");
    assert_eq!(plan.recipients()[0].target_pair_count(), 0);
}

#[test]
fn one_clean_source_can_be_reused_by_different_recipients_but_never_from_contaminated_output() {
    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![
        clean("A", "P1", 5),
        clean("B", "P1", 5),
        clean("C", "P1", 5),
    ];
    let settings = BiologicalContaminationSettings::new(60.0, 60.0, Some(vec![id("A"), id("B")]))
        .expect("subset");
    let reused = (0..=u8::MAX)
        .find_map(|last| {
            let plan =
                plan_biological_contamination(root(last), &event, &settings, &libraries).ok()?;
            let a = plan.recipients()[0]
                .assignments()
                .iter()
                .map(|item| item.source_pair())
                .collect::<BTreeSet<_>>();
            let b = plan.recipients()[1]
                .assignments()
                .iter()
                .map(|item| item.source_pair())
                .collect::<BTreeSet<_>>();
            let shared = a.intersection(&b).next().map(|pair| (*pair).clone());
            shared.map(|pair| (plan, pair))
        })
        .expect("at least one independent seed reuses an original clean source");
    let (plan, reused) = reused;
    assert!(
        plan.recipients()
            .iter()
            .all(
                |recipient| recipient.assignments().iter().all(|assignment| {
                    assignment.source_pair().source_individual_id() != recipient.recipient_id()
                })
            )
    );
    assert!(
        plan.recipients()
            .iter()
            .filter(|recipient| {
                recipient
                    .assignments()
                    .iter()
                    .any(|assignment| assignment.source_pair() == &reused)
            })
            .count()
            >= 2
    );
}

#[test]
fn manifest_record_requires_successful_materialization_and_records_realized_aggregates() {
    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![clean("A", "P1", 2), clean("B", "P1", 2)];
    let settings =
        BiologicalContaminationSettings::new(100.0, 100.0, Some(vec![id("A")])).expect("settings");
    let plan = plan_biological_contamination(root(2), &event, &settings, &libraries).expect("plan");
    let recipient = plan.recipients().first().expect("recipient");
    let mut materializer =
        BiologicalContaminationMaterializer::new(root(2), event, &libraries).expect("materializer");
    let report = materializer
        .emit_materialized_recipient(
            recipient,
            &libraries[0],
            |_| {
                Ok::<_, ()>((
                    b"R1".to_vec(),
                    b"II".to_vec(),
                    b"R2".to_vec(),
                    b"JJ".to_vec(),
                ))
            },
            |_| Ok::<_, ()>(()),
        )
        .expect("materialized");
    let record = recipient
        .realized_manifest_recipient(&report)
        .expect("successful materialization creates the manifest record");
    assert_eq!(record.recipient_id(), &id("A"));
    assert_eq!(record.clean_pair_count(), 2);
    assert_eq!(record.target_pair_count(), 2);
    assert_eq!(record.removed_own_pair_count(), 2);
    assert_eq!(record.inserted_contaminant_pair_count(), 2);
    assert_eq!(record.retained_own_pair_count(), 0);
    assert_eq!(record.final_pair_count(), 2);
    assert_eq!(record.donor_contributions().get(&id("B")), Some(&2));
    assert_ne!(record.percent_subseed(), record.donor_subseed());
    assert_ne!(record.donor_subseed(), record.removal_subseed());
    assert!(!format!("{record:?}").contains("PH_"));
    // El registro realizado no ofrece claves ni asignaciones temporales.
}

#[test]
fn materialization_replaces_recipient_pairs_with_new_recipient_instances_without_retaining_links() {
    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![clean("A", "P1", 2), clean("B", "P1", 2)];
    let settings =
        BiologicalContaminationSettings::new(50.0, 50.0, Some(vec![id("A")])).expect("settings");
    let plan =
        plan_biological_contamination(root(37), &event, &settings, &libraries).expect("plan");
    let recipient = plan.recipients().first().expect("recipient");
    let mut materializer = BiologicalContaminationMaterializer::new(root(37), event, &libraries)
        .expect("materializer");
    let mut emitted = Vec::new();
    let mut regenerated_sources = 0;
    let report = materializer
        .emit_materialized_recipient(
            recipient,
            &libraries[0],
            |source| {
                let donor = source.source_individual_id() == &id("B");
                if donor {
                    regenerated_sources += 1;
                }
                Ok::<_, ()>(if donor {
                    format!("observed:{}", source.pre_hop_instance_id())
                } else {
                    format!("own:{}", source.pre_hop_instance_id())
                })
            },
            |pair| {
                emitted.push((
                    pair.recipient_library_id().clone(),
                    pair.pre_hop_instance_id().clone(),
                    pair.is_contaminant(),
                    pair.into_payload(),
                ));
                Ok::<_, ()>(())
            },
        )
        .expect("materialize recipient");

    assert_eq!(report.removed_own_pairs(), 1);
    assert_eq!(report.inserted_contaminant_pairs(), 1);
    assert_eq!(report.retained_own_pairs(), 1);
    assert_eq!(
        regenerated_sources, 1,
        "each source is regenerated once on demand"
    );
    assert_eq!(
        emitted.len(),
        2,
        "replacement preserves PE-pair cardinality"
    );
    assert_eq!(
        emitted
            .iter()
            .filter(|(_, _, contaminant, _)| *contaminant)
            .count(),
        1
    );
    let contaminant = emitted
        .iter()
        .find(|(_, _, contaminant, _)| *contaminant)
        .expect("one new contaminant");
    assert_eq!(contaminant.0, library("LIB_A"));
    assert!(contaminant.1.as_str().starts_with("PH_"));
    assert!(
        libraries[1].own_pairs().iter().any(|source| {
            contaminant.3 == format!("observed:{}", source.pre_hop_instance_id().as_str())
        }),
        "the contaminant payload comes from a pair of the donor inventory"
    );
    assert!(
        !libraries[1]
            .own_pairs()
            .iter()
            .any(|source| source.pre_hop_instance_id() == contaminant.1),
        "the recipient gets a new PH instance and the donor inventory is untouched"
    );
    let aggregate = recipient
        .realized_manifest_recipient(&report)
        .expect("only the successful report makes a durable aggregate");
    assert_eq!(aggregate.donor_contributions().get(&id("B")), Some(&1));
    assert_eq!(aggregate.final_pair_count(), 2);
    // The aggregate manifest view deliberately exposes neither source keys nor
    // the fresh PH instance: individual provenance dies with the emission.
}

#[test]
fn materialization_preserves_paired_donor_bytes_at_zero_and_hundred_percent() {
    #[derive(Clone, Debug, Eq, PartialEq)]
    struct PairedPayload {
        r1_sequence: Vec<u8>,
        r1_quality: Vec<u8>,
        r2_sequence: Vec<u8>,
        r2_quality: Vec<u8>,
    }

    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![clean("A", "P1", 2), clean("B", "P1", 2)];
    let donor_payload = PairedPayload {
        r1_sequence: b"ACTG".to_vec(),
        r1_quality: b"!I#J".to_vec(),
        r2_sequence: b"CAGT".to_vec(),
        r2_quality: b"J#!I".to_vec(),
    };
    for (percent, expected_contaminants) in [(0.0, 0_u64), (100.0, 2_u64)] {
        let settings = BiologicalContaminationSettings::new(percent, percent, Some(vec![id("A")]))
            .expect("settings");
        let plan =
            plan_biological_contamination(root(44), &event, &settings, &libraries).expect("plan");
        let recipient = plan.recipients().first().expect("recipient");
        let mut materializer =
            BiologicalContaminationMaterializer::new(root(44), event.clone(), &libraries)
                .expect("materializer");
        let mut emitted = Vec::new();
        let report = materializer
            .emit_materialized_recipient(
                recipient,
                &libraries[0],
                |key| {
                    Ok::<_, ()>(if key.source_individual_id() == &id("B") {
                        donor_payload.clone()
                    } else {
                        PairedPayload {
                            r1_sequence: b"AAAA".to_vec(),
                            r1_quality: b"IIII".to_vec(),
                            r2_sequence: b"TTTT".to_vec(),
                            r2_quality: b"JJJJ".to_vec(),
                        }
                    })
                },
                |pair| {
                    emitted.push((pair.is_contaminant(), pair.into_payload()));
                    Ok::<_, ()>(())
                },
            )
            .expect("materialized");
        assert_eq!(report.inserted_contaminant_pairs(), expected_contaminants);
        assert_eq!(report.final_pair_count(), 2);
        for (_, payload) in emitted.iter().filter(|(contaminant, _)| *contaminant) {
            assert_eq!(
                payload, &donor_payload,
                "paired sequence and quality copy exactly"
            );
        }
    }
}

#[test]
fn materialization_keeps_only_one_payload_live_at_a_time_with_retained_and_contaminant_pairs() {
    struct LivePayload(Arc<AtomicUsize>);
    impl Drop for LivePayload {
        fn drop(&mut self) {
            self.0.fetch_sub(1, Ordering::SeqCst);
        }
    }
    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![clean("A", "P1", 8), clean("B", "P1", 8)];
    let settings =
        BiologicalContaminationSettings::new(50.0, 50.0, Some(vec![id("A")])).expect("settings");
    let plan =
        plan_biological_contamination(root(55), &event, &settings, &libraries).expect("plan");
    let live = Arc::new(AtomicUsize::new(0));
    let peak = Arc::new(AtomicUsize::new(0));
    let mut materializer = BiologicalContaminationMaterializer::new(root(55), event, &libraries)
        .expect("materializer");
    materializer
        .emit_materialized_recipient(
            plan.recipients().first().expect("recipient"),
            &libraries[0],
            |_| {
                let current = live.fetch_add(1, Ordering::SeqCst) + 1;
                peak.fetch_max(current, Ordering::SeqCst);
                Ok::<_, ()>(LivePayload(Arc::clone(&live)))
            },
            |_| Ok::<_, ()>(()),
        )
        .expect("streamed materialization");
    assert_eq!(
        live.load(Ordering::SeqCst),
        0,
        "the emitter retained no payload"
    );
    assert_eq!(peak.load(Ordering::SeqCst), 1, "payload memory is bounded");
}

#[test]
fn materialization_report_cannot_certify_another_recipient_with_equal_counts() {
    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![
        clean("A", "P1", 2),
        clean("B", "P1", 2),
        clean("C", "P1", 2),
    ];
    let settings = BiologicalContaminationSettings::new(100.0, 100.0, Some(vec![id("A"), id("B")]))
        .expect("settings");
    let plan =
        plan_biological_contamination(root(67), &event, &settings, &libraries).expect("plan");
    let mut materializer = BiologicalContaminationMaterializer::new(root(67), event, &libraries)
        .expect("materializer");
    let report_a = materializer
        .emit_materialized_recipient(
            &plan.recipients()[0],
            &libraries[0],
            |_| Ok::<_, ()>(()),
            |_| Ok::<_, ()>(()),
        )
        .expect("A materialized");

    assert_eq!(
        report_a.final_pair_count(),
        plan.recipients()[1].clean_pair_count(),
        "the regression uses equal aggregate counts"
    );
    assert!(
        plan.recipients()[0]
            .realized_manifest_recipient(&report_a)
            .is_some()
    );
    assert!(
        plan.recipients()[1]
            .realized_manifest_recipient(&report_a)
            .is_none()
    );
}

#[test]
fn reused_clean_source_materializes_as_distinct_recipient_instances() {
    let event = EventId::parse("contamination-main", "test").expect("event");
    let libraries = vec![
        clean("A", "P1", 1),
        clean("B", "P1", 1),
        clean("C", "P1", 1),
    ];
    let settings = BiologicalContaminationSettings::new(100.0, 100.0, Some(vec![id("A"), id("B")]))
        .expect("settings");
    let (seed, plan) = (0..=u8::MAX)
        .find_map(|last| {
            let plan =
                plan_biological_contamination(root(last), &event, &settings, &libraries).ok()?;
            (plan.recipients()[0].assignments()[0].source_pair()
                == plan.recipients()[1].assignments()[0].source_pair())
            .then_some((last, plan))
        })
        .expect("one seed reuses the sole common clean source");
    let source_a = plan.recipients()[0].assignments()[0].source_pair();
    let source_b = plan.recipients()[1].assignments()[0].source_pair();
    assert_eq!(
        source_a, source_b,
        "both recipients reuse the same clean source"
    );

    let mut materializer = BiologicalContaminationMaterializer::new(root(seed), event, &libraries)
        .expect("materializer");
    let mut instances = Vec::new();
    for (recipient, inventory) in plan.recipients().iter().zip(&libraries[..2]) {
        materializer
            .emit_materialized_recipient(
                recipient,
                inventory,
                |_| Ok::<_, ()>(()),
                |pair| {
                    instances.push((
                        pair.recipient_library_id().clone(),
                        pair.pre_hop_instance_id().clone(),
                    ));
                    Ok::<_, ()>(())
                },
            )
            .expect("materialized recipient");
    }
    assert_eq!(instances.len(), 2);
    assert_ne!(instances[0].1, instances[1].1);
    assert_ne!(instances[0].0, instances[1].0);
}

#[test]
fn configuration_resolves_default_or_explicit_recipients_and_rejects_external_matrix_and_weights() {
    let base = "config_kind = \"general\"\nschema_version = \"2\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1000\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\", \"sample-B\"]\n\n[technical_event]\nevent_id = \"contamination-main\"\nindividual_ids = [\"sample-A\", \"sample-B\"]\ntarget_coverage_wgs = 1\n\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n\n[[technical_event.libraries]]\nindividual_id = \"sample-B\"\ni7 = \"CGTACGTA\"\ni5 = \"GCATGCAT\"\n";
    let explicit = base.replacen(
        "\n[[technical_event.libraries]]",
        "\ncontamination_min_percent = 10\ncontamination_max_percent = 20\ncontamination_recipient_ids = [\"sample-B\"]\n\n[[technical_event.libraries]]",
        1,
    );
    let declared = parse_declared_config(
        admit_toml_header(
            std::path::Path::new("contamination.toml"),
            explicit.into_bytes(),
            ExpectedConfigKind::General,
        )
        .expect("header"),
    )
    .expect("declared");
    let effective = build_effective_config(&declared).expect("effective");
    let settings = effective
        .technical_event()
        .expect("event")
        .biological_contamination();
    assert_eq!(
        (settings.minimum_percent(), settings.maximum_percent()),
        (10.0, 20.0)
    );
    assert_eq!(
        settings.configured_recipients(),
        Some(&[id("sample-B")][..])
    );

    let defaulted = build_effective_config(
        &parse_declared_config(
            admit_toml_header(
                std::path::Path::new("default-contamination.toml"),
                base.as_bytes().to_vec(),
                ExpectedConfigKind::General,
            )
            .expect("header"),
        )
        .expect("declared"),
    )
    .expect("effective");
    let defaulted = defaulted
        .technical_event()
        .expect("event")
        .biological_contamination();
    assert_eq!(
        (defaulted.minimum_percent(), defaulted.maximum_percent()),
        (0.0, 0.0)
    );
    assert_eq!(defaulted.configured_recipients(), None);

    for forbidden in [
        "external_fastq = \"source.fq.gz\"",
        "matrix = \"directed.tsv\"",
        "donor_weights = [1, 2]",
    ] {
        let invalid = base.replacen(
            "\n[[technical_event.libraries]]",
            &format!("\ncontamination_min_percent = 0\ncontamination_max_percent = 0\n{forbidden}\n\n[[technical_event.libraries]]"),
            1,
        );
        assert!(
            parse_declared_config(
                admit_toml_header(
                    std::path::Path::new("invalid-contamination.toml"),
                    invalid.into_bytes(),
                    ExpectedConfigKind::General,
                )
                .expect("header")
            )
            .is_err(),
            "{forbidden}"
        );
    }
    let nested = base.replacen(
        "\n[[technical_event.libraries]]",
        "\n[technical_event.biological_contamination]\nminimum_percent = 10\nmaximum_percent = 20\n\n[[technical_event.libraries]]",
        1,
    );
    assert!(
        parse_declared_config(
            admit_toml_header(
                std::path::Path::new("nested-contamination.toml"),
                nested.into_bytes(),
                ExpectedConfigKind::General,
            )
            .expect("header"),
        )
        .is_err(),
        "la forma anidada histórica no debe volver a entrar por el parser"
    );
}

#[test]
fn positive_contamination_and_hopping_publish_reconciled_aggregate_manifest() {
    let input_root = std::env::temp_dir().join(format!(
        "gps-positive-contamination-input-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ));
    fs::create_dir_all(&input_root).expect("input root");
    let fasta = input_root.join("reference.fa");
    let matrix = input_root.join("population-snvs.tsv");
    fs::write(&fasta, format!(">chr1\n{}\n", "A".repeat(1000))).expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t500\tA\tT\tP1\t0.5\n",
    )
    .expect("matrix");
    let text = "config_kind = \"general\"\nschema_version = \"2\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1000\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\", \"sample-B\"]\n\n[technical_event]\nevent_id = \"contamination-main\"\nindividual_ids = [\"sample-A\", \"sample-B\"]\ntarget_coverage_wgs = 1\nsequencing_errors = \"disabled\"\ncontamination_min_percent = 50\ncontamination_max_percent = 50\ncontamination_recipient_ids = [\"sample-A\"]\nindex_hopping_rate = 1\n\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n\n[[technical_event.libraries]]\nindividual_id = \"sample-B\"\ni7 = \"CGTACGTA\"\ni5 = \"GCATGCAT\"\n";
    let text = text.replacen(
        "[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n",
        &format!(
            "[population_snvs]\nmode = \"supplied\"\ninput_path = \"{}\"\n\n[reference]\nkind = \"supplied\"\nid = \"reference-main\"\nfasta = \"{}\"\n",
            matrix.display(),
            fasta.display()
        ),
        1,
    );
    let declared = parse_declared_config(
        admit_toml_header(
            std::path::Path::new("positive-contamination.toml"),
            text.as_bytes().to_vec(),
            ExpectedConfigKind::General,
        )
        .expect("header"),
    )
    .expect("declared");
    let effective = build_effective_config(&declared).expect("effective");
    let destination = std::env::temp_dir().join(format!(
        "gps-positive-contamination-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ));
    execute_validated_config(&effective, &destination).expect("F4 final PE publication");
    let manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(destination.join("events/contamination-main/manifest.json")).expect("manifest"),
    )
    .expect("JSON");
    assert_eq!(manifest["error_quality"]["effective_mode"], "disabled");
    assert_eq!(
        manifest["error_quality"]["qualities_represent_realized_error_probability"],
        false
    );
    for mate in ["r1", "r2"] {
        assert_eq!(manifest["error_quality"][mate]["errors"], 0);
    }
    let demultiplexing = &manifest["demultiplexing"];
    let reassigned = demultiplexing["reassigned_pairs"].as_u64().expect("total");
    let reassigned_entries = demultiplexing["reassigned_by_sample"]
        .as_array()
        .expect("entries");
    assert!(!reassigned_entries.is_empty());
    assert_eq!(
        reassigned_entries
            .iter()
            .map(|entry| entry["pair_count"].as_u64().expect("count"))
            .sum::<u64>(),
        reassigned
    );
    assert!(
        reassigned_entries
            .iter()
            .all(|entry| entry["cause"] == "index_hopping"
                && entry["source_library_id"] != entry["destination_library_id"])
    );
    let fastq_pairs = manifest["final_fastq_summary"]
        .as_array()
        .expect("summary")
        .iter()
        .map(|entry| entry["pair_count"].as_u64().expect("pair count"))
        .sum::<u64>();
    assert_eq!(
        reassigned, fastq_pairs,
        "rate-one hopping preserves all final pairs"
    );
    assert!(
        manifest["final_fastq_summary"]
            .as_array()
            .expect("summary")
            .iter()
            .any(|entry| entry["contaminant_pair_count"]
                .as_u64()
                .expect("contaminants")
                > 0)
    );
    let outputs = manifest["outputs"].as_array().expect("outputs");
    assert_eq!(outputs.len(), 4, "two samples publish R1 and R2 only");
    for mate in [1, 2] {
        assert_eq!(
            outputs
                .iter()
                .filter(|output| output["mate"] == mate && output["state"] == "complete")
                .count(),
            2,
            "each sample has a complete R{mate}"
        );
    }
    let serialized = serde_json::to_string(&manifest).expect("serialize");
    for forbidden in ["PH_", "TP_", "technical_pair_id", "pair_ids"] {
        assert!(
            !serialized.contains(forbidden),
            "durable provenance: {forbidden}"
        );
    }
    fs::remove_dir_all(destination).expect("cleanup");
    fs::remove_dir_all(input_root).expect("input cleanup");
}
