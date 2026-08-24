use genome_population_simulator::coordinates::BasePosition;
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{IndividualId, ReferenceId, ScenarioId, SequenceId};
use genome_population_simulator::indel_application::{
    IndelApplicationError, IndelClass, RandomIndelRequest, RandomSmallIndelRequest,
    SmallIndelLengths, apply_random_indels, apply_random_small_indels, apply_sampled_indels,
    apply_sampled_small_indels,
};
use genome_population_simulator::materialization::{
    GenotypeQueryError, SnpLocus, derive_snp_genotype, materialize_copy,
};
use genome_population_simulator::randomness::{RootSeed, RuntimeSubseedScope, derive_subseed};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{CanonicalBase, CurrentIndelEvent, ModificationOrigin};
use std::sync::OnceLock;

fn reference(bases: &[u8]) -> Reference {
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("sequence"),
                FastaHeader::parse(b"chr1").expect("header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(bases).expect("DNA"),
            )
            .expect("sequence"),
        ],
    )
    .expect("reference")
}

fn individual(reference: &Reference) -> Individual {
    Individual::from_reference(
        IndividualId::parse("individual-1", "test").expect("individual"),
        reference,
        Ploidy::new(2).expect("ploidy"),
    )
}

fn request<'a>(reference: &'a Reference, class: IndelClass, rate: f64) -> RandomIndelRequest<'a> {
    RandomIndelRequest {
        scenario_id: scenario(),
        sequence_id: reference.sequences()[0].id(),
        phase: HaplotypePhase::from_copy_index(0).expect("phase"),
        spatial_scope_id: "global",
        class,
        rate_per_opportunity: rate,
    }
}

fn scenario() -> &'static ScenarioId {
    static SCENARIO: OnceLock<ScenarioId> = OnceLock::new();
    SCENARIO.get_or_init(|| ScenarioId::parse("direct", "test").expect("scenario"))
}

fn root() -> RootSeed {
    RootSeed::from_bytes([42; 32])
}

fn first_copy_bases(individual: &Individual, reference: &Reference) -> Vec<u8> {
    materialize_copy(
        &individual
            .copies_for(reference.sequences()[0].id())
            .expect("copies")[0],
        reference,
    )
    .expect("materialization")
}

#[test]
fn exact_indel_lengths_are_independent_and_stay_below_the_effective_threshold() {
    assert_eq!(SmallIndelLengths::default().insertion_length(), 1);
    assert_eq!(SmallIndelLengths::default().deletion_length(), 1);
    let lengths = SmallIndelLengths::new(1, 2, 5).expect("valid independent lengths");
    assert_eq!(lengths.insertion_length(), 1);
    assert_eq!(lengths.deletion_length(), 2);
    let threshold_minus_one = SmallIndelLengths::new(4, 1, 5).expect("threshold - 1");
    assert_eq!(threshold_minus_one.insertion_length(), 4);
    assert!(matches!(
        SmallIndelLengths::new(0, 1, 5),
        Err(IndelApplicationError::InvalidLength {
            class: IndelClass::Insertion,
            ..
        })
    ));
    assert!(matches!(
        SmallIndelLengths::new(1, 5, 5),
        Err(IndelApplicationError::InvalidLength {
            class: IndelClass::Deletion,
            ..
        })
    ));
}

#[test]
fn threshold_minus_one_is_applied_as_an_exact_insertion_length() {
    let reference = reference(b"ACG");
    let original = individual(&reference);
    let request = request(&reference, IndelClass::Insertion, 1.0);
    let result = apply_sampled_indels(
        &original,
        &reference,
        RootSeed::from_bytes([0; 32]),
        SmallIndelLengths::new(4, 1, 5).expect("valid length"),
        &request,
        1,
    )
    .expect("one exact insertion");
    assert_eq!(first_copy_bases(result.individual(), &reference).len(), 7);
    let CurrentIndelEvent::Insertion {
        bases, provenance, ..
    } = &result
        .individual()
        .copies_for(reference.sequences()[0].id())
        .expect("copies")[0]
        .segments()
        .current_indel_events()[0]
    else {
        panic!("expected insertion");
    };
    assert_eq!(bases.len(), 4);
    assert_eq!(provenance.origin(), ModificationOrigin::Random);
}

#[test]
fn insertion_reenumerates_the_current_state_and_can_use_an_anchor_it_created() {
    let reference = reference(b"AC");
    let original = individual(&reference);
    let request = request(&reference, IndelClass::Insertion, 1.0);
    let result = apply_sampled_indels(
        &original,
        &reference,
        RootSeed::from_bytes([1; 32]),
        SmallIndelLengths::default(),
        &request,
        2,
    )
    .expect("two dynamic insertions");
    assert_eq!(result.accepted_count(), 2);
    assert_eq!(first_copy_bases(result.individual(), &reference).len(), 4);
    let events = result
        .individual()
        .copies_for(reference.sequences()[0].id())
        .expect("copies")[0]
        .segments()
        .current_indel_events();
    assert_eq!(events.len(), 2);
    assert!(events.iter().all(
        |event| matches!(event, CurrentIndelEvent::Insertion { bases, .. } if bases.len() == 1)
    ));
    let CurrentIndelEvent::Insertion { point, .. } = &events[1] else {
        panic!("expected second insertion");
    };
    assert_eq!(
        *point, 2,
        "the second event anchors after the inserted base"
    );
    assert_eq!(original, individual(&reference));
}

#[test]
fn deletion_reenumeration_eliminates_opportunities_and_rolls_back_the_full_class() {
    let reference = reference(b"ACGT");
    let original = individual(&reference);
    let request = request(&reference, IndelClass::Deletion, 1.0);
    let error = apply_sampled_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &request,
        3,
    )
    .expect_err("the third one cannot be accepted");
    assert_eq!(
        error,
        IndelApplicationError::ExhaustedOpportunities {
            sampled_count: 3,
            accepted_before_failure: 2,
            active_before_rollback: 2,
            opportunities_at_failure: 0,
        }
    );
    assert_eq!(first_copy_bases(&original, &reference), b"ACGT");
}

#[test]
fn terminal_deletions_and_n_anchors_are_not_opportunities_but_insertion_content_is_acgt() {
    let short_reference = reference(b"AC");
    let short_individual = individual(&short_reference);
    let deletion = request(&short_reference, IndelClass::Deletion, 1.0);
    assert!(matches!(
        apply_sampled_indels(
            &short_individual,
            &short_reference,
            root(),
            SmallIndelLengths::default(),
            &deletion,
            1,
        ),
        Err(IndelApplicationError::ExhaustedOpportunities {
            accepted_before_failure: 0,
            ..
        })
    ));

    let n_reference = reference(b"NAC");
    let n_individual = individual(&n_reference);
    let insertion = request(&n_reference, IndelClass::Insertion, 1.0);
    let result = apply_sampled_indels(
        &n_individual,
        &n_reference,
        root(),
        SmallIndelLengths::default(),
        &insertion,
        1,
    )
    .expect("only the canonical anchor after A is usable");
    let events = result
        .individual()
        .copies_for(n_reference.sequences()[0].id())
        .expect("copies")[0]
        .segments()
        .current_indel_events();
    let CurrentIndelEvent::Insertion { point, bases, .. } = &events[0] else {
        panic!("expected insertion");
    };
    assert_eq!(*point, 2);
    assert!(
        bases
            .iter()
            .all(|base| matches!(*base, b'A' | b'C' | b'G' | b'T'))
    );
}

#[test]
fn zero_rate_leaves_the_state_unchanged_without_a_count_or_position_consumer() {
    let reference = reference(b"ACGT");
    let original = individual(&reference);
    let request = request(&reference, IndelClass::Insertion, 0.0);
    let result = apply_random_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &request,
    )
    .expect("zero is valid");
    assert_eq!(result.individual(), &original);
    assert_eq!(result.sampled_count(), 0);
    assert_eq!(result.accepted_count(), 0);
}

#[test]
fn joint_entrypoint_always_applies_insertions_before_deletions() {
    let reference = reference(b"ACG");
    let original = individual(&reference);
    let scenario = ScenarioId::parse("direct", "test").expect("scenario");
    let request = RandomSmallIndelRequest {
        scenario_id: &scenario,
        sequence_id: reference.sequences()[0].id(),
        phase: HaplotypePhase::from_copy_index(0).expect("phase"),
        spatial_scope_id: "global",
        insertion_rate_per_opportunity: 0.0,
        deletion_rate_per_opportunity: 0.0,
    };
    let result = apply_random_small_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &request,
    )
    .expect("joint zero-rate call");
    assert_eq!(result.insertion().sampled_count(), 0);
    assert_eq!(result.deletion().sampled_count(), 0);
}

#[test]
fn joint_sampled_path_uses_the_insertion_state_to_create_a_deletion_opportunity() {
    let reference = reference(b"AC");
    let original = individual(&reference);
    let request = RandomSmallIndelRequest {
        scenario_id: scenario(),
        sequence_id: reference.sequences()[0].id(),
        phase: HaplotypePhase::from_copy_index(0).expect("phase"),
        spatial_scope_id: "global",
        insertion_rate_per_opportunity: 1.0,
        deletion_rate_per_opportunity: 1.0,
    };
    let result = apply_sampled_small_indels(
        &original,
        &reference,
        RootSeed::from_bytes([1; 32]),
        SmallIndelLengths::default(),
        &request,
        1,
        1,
    )
    .expect("the insertion creates the only internal deletion opportunity");
    assert_eq!(result.insertion().initial_opportunities(), 1);
    assert_eq!(result.insertion().sampled_count(), 1);
    assert_eq!(result.deletion().initial_opportunities(), 1);
    assert_eq!(result.deletion().sampled_count(), 1);
    assert_eq!(result.deletion().accepted_count(), 1);
    assert_eq!(first_copy_bases(result.individual(), &reference), b"AC");
}

#[test]
fn a_two_base_deletion_has_exact_length_and_leaves_both_terminal_bases() {
    let reference = reference(b"ACGTA");
    let original = individual(&reference);
    let request = request(&reference, IndelClass::Deletion, 1.0);
    let result = apply_sampled_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::new(1, 2, 5).expect("valid lengths"),
        &request,
        1,
    )
    .expect("one deletion of two bases");
    assert_eq!(first_copy_bases(result.individual(), &reference).len(), 3);
    let CurrentIndelEvent::Deletion { start, end, .. } = &result
        .individual()
        .copies_for(reference.sequences()[0].id())
        .expect("copies")[0]
        .segments()
        .current_indel_events()[0]
    else {
        panic!("expected deletion");
    };
    assert_eq!(end - start, 2);
    let bases = first_copy_bases(result.individual(), &reference);
    assert_eq!(bases.first(), Some(&b'A'));
    assert_eq!(bases.last(), Some(&b'A'));
}

#[test]
fn random_indel_classes_are_reproducible_with_their_own_count_and_position_scopes() {
    let reference = reference(b"ACGTACGTAC");
    let original = individual(&reference);
    let insertion = request(&reference, IndelClass::Insertion, 0.4);
    let deletion = request(&reference, IndelClass::Deletion, 0.4);
    let first_insertion = apply_random_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &insertion,
    )
    .expect("random insertions");
    let repeated_insertion = apply_random_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &insertion,
    )
    .expect("repeated random insertions");
    let first_deletion = apply_random_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &deletion,
    )
    .expect("random deletions");
    let repeated_deletion = apply_random_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &deletion,
    )
    .expect("repeated random deletions");
    assert_eq!(first_insertion, repeated_insertion);
    assert_eq!(first_deletion, repeated_deletion);
    assert_eq!(first_insertion.initial_opportunities(), 9);
    assert_eq!(first_insertion.sampled_count(), 5);
    assert_eq!(
        first_copy_bases(first_insertion.individual(), &reference),
        b"AGCGTACGTCTCAGC",
    );
    assert_eq!(first_deletion.initial_opportunities(), 8);
    assert_eq!(first_deletion.sampled_count(), 5);
    assert_eq!(
        first_copy_bases(first_deletion.individual(), &reference),
        b"AAGAC",
    );
    let count_scope = RuntimeSubseedScope::IndelCount {
        scenario_id: "direct",
        reference_id: "ref-main",
        individual_id: "individual-1",
        sequence_id: "chr1",
        variant_class: "insertion",
        spatial_scope_id: "global",
    };
    let position_scope = RuntimeSubseedScope::IndelPosition {
        scenario_id: "direct",
        reference_id: "ref-main",
        individual_id: "individual-1",
        sequence_id: "chr1",
        variant_class: "insertion",
        spatial_scope_id: "global",
    };
    assert_ne!(
        derive_subseed(root(), count_scope).expect("count scope"),
        derive_subseed(root(), position_scope).expect("position scope"),
    );
}

#[test]
fn later_calls_keep_inserted_provenance_distinct_and_a_deleted_source_locus_is_absent() {
    let reference = reference(b"ACGT");
    let original = individual(&reference);
    let insertion = request(&reference, IndelClass::Insertion, 1.0);
    let once = apply_sampled_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &insertion,
        1,
    )
    .expect("first insertion");
    let twice = apply_sampled_indels(
        once.individual(),
        &reference,
        root(),
        SmallIndelLengths::default(),
        &insertion,
        1,
    )
    .expect("second insertion");
    let insertion_events = twice
        .individual()
        .copies_for(reference.sequences()[0].id())
        .expect("copies")[0]
        .segments()
        .current_indel_events();
    assert_ne!(
        insertion_events[0].event_id(),
        insertion_events[1].event_id()
    );

    let deletion = request(&reference, IndelClass::Deletion, 1.0);
    let deleted = apply_sampled_indels(
        &original,
        &reference,
        root(),
        SmallIndelLengths::default(),
        &deletion,
        1,
    )
    .expect("one deletion");
    let locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            reference.sequences()[0].id().clone(),
            1,
            4,
        )
        .expect("position"),
        CanonicalBase::C,
        CanonicalBase::T,
    )
    .expect("locus");
    assert!(matches!(
        derive_snp_genotype(deleted.individual(), &reference, &locus),
        Err(GenotypeQueryError::LocusDeletedFromCopy { .. })
    ));
}
