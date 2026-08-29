//! F4.0.4: aislamiento entre las decisiones técnicas ya publicadas.
//!
//! Esta prueba permanece deliberadamente en la frontera de claves y primitivas:
//! F4.0 no tiene todavía productores de cobertura, fragmentos o parejas. El
//! vector representa sus causas canónicas, no lecturas hipotéticas.

use std::collections::BTreeMap;
use std::thread;

use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId};
use genome_population_simulator::randomness::{
    F4IndexHopInstance, F4Mate, F4TechnicalDecision, RootSeed, derive_f4_technical_subseed,
};
use genome_population_simulator::technical_identity::{
    F4TechnicalIdentityDeriver, TemplateOrdinal, derive_technical_identity_plan,
};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
const TEMPLATE_ORDINAL: u64 = 17;

#[derive(Clone, Debug, Eq, PartialEq)]
struct DecisionVector {
    coverage: (String, u64),
    fragment: (String, u64),
    chimera: (String, u64),
    error: (String, u64),
    contamination: (String, u64),
    index_hop: (String, u64),
}

fn root() -> RootSeed {
    RootSeed::parse("root", ROOT).expect("literal root")
}

fn event() -> EventId {
    EventId::parse("event-isolation", "test").expect("event")
}

fn reference() -> ReferenceId {
    ReferenceId::parse("reference-isolation", "test").expect("reference")
}

fn individual(value: &str) -> IndividualId {
    IndividualId::parse(value, "test").expect("individual")
}

fn observed_vector(individual_order: &[&str]) -> DecisionVector {
    let individual_ids = individual_order
        .iter()
        .map(|value| individual(value))
        .collect::<Vec<_>>();
    let plan = derive_technical_identity_plan(root(), &reference(), event(), &individual_ids)
        .expect("technical identities");
    let source_individual = individual("sample-A");
    let recipient_individual = individual("sample-B");
    let source_library = plan.library(&source_individual).expect("source library");
    let recipient_library = plan
        .library(&recipient_individual)
        .expect("recipient library");
    let template = plan
        .derive_template(&source_individual, TemplateOrdinal::new(TEMPLATE_ORDINAL))
        .expect("template");
    let pre_hop_id = F4TechnicalIdentityDeriver::new(root(), event())
        .derive_own_pre_hop_instance(source_library.library_id(), template.template_id())
        .expect("pre-hop identity");
    let pre_hop = F4IndexHopInstance::parse(pre_hop_id.as_str()).expect("PH identity");

    let observe = |decision| {
        let subseed = derive_f4_technical_subseed(root(), decision).expect("F4 decision seed");
        let word = subseed.local_rng().next_u64();
        (subseed.to_hex(), word)
    };

    DecisionVector {
        coverage: {
            let subseed = plan
                .derive_library_count_subseed(&source_individual)
                .expect("coverage seed");
            let word = subseed.local_rng().next_u64();
            (subseed.to_hex(), word)
        },
        fragment: observe(F4TechnicalDecision::fragment_insert_length(
            plan.event_id(),
            source_library.library_id(),
            TEMPLATE_ORDINAL,
        )),
        chimera: observe(F4TechnicalDecision::chimera_bernoulli(
            plan.event_id(),
            source_library.library_id(),
            TEMPLATE_ORDINAL,
        )),
        error: observe(F4TechnicalDecision::error_bernoulli(
            plan.event_id(),
            source_library.library_id(),
            TEMPLATE_ORDINAL,
            F4Mate::R2,
            11,
        )),
        contamination: observe(F4TechnicalDecision::contamination_percent(
            plan.event_id(),
            recipient_library.library_id(),
            &recipient_individual,
        )),
        index_hop: observe(F4TechnicalDecision::index_hop_destination(
            plan.event_id(),
            source_library.library_id(),
            pre_hop,
        )),
    }
}

#[test]
fn f4_decisions_regenerate_and_remain_isolated_from_library_order_and_threads() {
    let expected = DecisionVector {
        coverage: (
            "143fc606a46110ea9529475287cb7212f23d8c3266d83320fb2f509614d4090f".to_owned(),
            16_274_260_773_103_446_090,
        ),
        fragment: (
            "a151890f2eb70d0b455a1483bf61d33151a21684bc8d6f9cbf44e7d50922ec37".to_owned(),
            18_134_805_061_096_423_665,
        ),
        chimera: (
            "ddff1bcc48edd5c8ef16caa90abf85df6eda4865f201fa9a4f834b6178ffe131".to_owned(),
            18_118_388_711_401_500_271,
        ),
        error: (
            "60f8baea120021612a1d5cd4d943d403e99200a6ae60107ff4bad57c449d8b5f".to_owned(),
            16_277_303_102_483_717_294,
        ),
        contamination: (
            "fe3ff2943f5e3bf0bd8b91d1a5ca50c5905846344ab1a2f581476c5665479003".to_owned(),
            5_141_541_243_848_133_472,
        ),
        index_hop: (
            "6de8c8493eae769e94676b4d71a817aa91b6c44cc9e7f1f0290dafb5553d5de7".to_owned(),
            15_919_244_004_831_785_186,
        ),
    };
    let baseline = observed_vector(&["sample-A", "sample-B", "sample-C"]);

    // Regenerar las mismas causas no depende de una ejecución previa.
    assert_eq!(
        baseline,
        observed_vector(&["sample-A", "sample-B", "sample-C"])
    );
    // Las identidades de biblioteca son canónicas, no posiciones de entrada.
    assert_eq!(
        baseline,
        observed_vector(&["sample-C", "sample-B", "sample-A"])
    );

    // Dos tareas independientes representan una partición de trabajo distinta.
    let first = thread::spawn(|| observed_vector(&["sample-A", "sample-B", "sample-C"]));
    let second = thread::spawn(|| observed_vector(&["sample-C", "sample-A", "sample-B"]));
    assert_eq!(first.join().expect("first task"), baseline);
    assert_eq!(second.join().expect("second task"), baseline);

    assert_eq!(baseline, expected, "manual F4.0.4 vector");
}

#[derive(Clone, Copy)]
struct SnapshotInstance<'a> {
    pre_hop_library: &'a str,
    pre_hop_instance: &'a str,
    destination: &'a str,
}

fn apply_simultaneous_transfers<'a>(
    snapshot: &[SnapshotInstance<'a>],
    evaluation_order: &[usize],
) -> BTreeMap<&'a str, Vec<&'a str>> {
    let mut final_libraries = BTreeMap::<&'a str, Vec<&'a str>>::new();
    for instance in snapshot {
        final_libraries.entry(instance.pre_hop_library).or_default();
        final_libraries.entry(instance.destination).or_default();
    }
    // Cada decisión se evalúa contra la instantánea pre-salto. En particular,
    // `pre_hop_library` no se actualiza cuando otra transferencia llega antes.
    let transfers = evaluation_order
        .iter()
        .map(|index| snapshot[*index])
        .collect::<Vec<_>>();
    for transfer in transfers {
        final_libraries
            .get_mut(transfer.destination)
            .expect("destination from snapshot")
            .push(transfer.pre_hop_instance);
    }
    final_libraries
}

#[test]
fn index_hop_fixture_applies_a_snapshot_simultaneously_without_chains() {
    let snapshot = [
        SnapshotInstance {
            pre_hop_library: "LIB_A",
            pre_hop_instance: "PH_A",
            destination: "LIB_B",
        },
        SnapshotInstance {
            pre_hop_library: "LIB_B",
            pre_hop_instance: "PH_B",
            destination: "LIB_C",
        },
    ];

    let forward = apply_simultaneous_transfers(&snapshot, &[0, 1]);
    let reverse = apply_simultaneous_transfers(&snapshot, &[1, 0]);
    assert_eq!(forward, reverse, "transfer order has no authority");
    assert_eq!(forward["LIB_A"], Vec::<&str>::new());
    assert_eq!(forward["LIB_B"], vec!["PH_A"]);
    assert_eq!(forward["LIB_C"], vec!["PH_B"]);
}
