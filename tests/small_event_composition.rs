use genome_population_simulator::coordinates::BasePosition;
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{
    EventId, IndividualId, ReferenceId, ScenarioId, SequenceId,
};
use genome_population_simulator::indel_application::SmallIndelLengths;
use genome_population_simulator::materialization::{
    GenotypeQueryError, SnpLocus, derive_snp_genotype, materialize_copy,
};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{CanonicalBase, ModificationProvenance, SegmentError};
use genome_population_simulator::small_event_composition::{
    AppliedSmallEvent, DeclaredSmallEvent, DeclaredSmallEventKind, DeclaredSmallIndelRequest,
    RandomSmallEventBlock, SampledSmallIndelCounts, SmallEventCompositionError,
    SmallEventSourceError, TargetedCurrentDeletion, TargetedCurrentInsertion,
    apply_declared_small_events, preflight_declared_small_events,
};
use genome_population_simulator::snp_application::{ApplicableSnp, TargetedSnp};

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

fn snp(
    reference: &Reference,
    offset: u64,
    alternate: CanonicalBase,
    event_id: &str,
) -> TargetedSnp {
    let sequence = &reference.sequences()[0];
    let reference_base = CanonicalBase::from_byte(
        sequence.sequence().as_bytes()[usize::try_from(offset).expect("offset")],
    )
    .expect("canonical base");
    TargetedSnp::new(
        HaplotypePhase::from_copy_index(0).expect("phase"),
        ApplicableSnp::explicit(
            SnpLocus::new(
                BasePosition::new(
                    reference.id().clone(),
                    sequence.id().clone(),
                    offset,
                    sequence.sequence().len(),
                )
                .expect("position"),
                reference_base,
                alternate,
            )
            .expect("locus"),
            EventId::parse(event_id, "test").expect("event"),
        ),
    )
}

fn random_indels(
    reference: &Reference,
    insertion_count: u64,
    deletion_count: u64,
) -> RandomSmallEventBlock {
    RandomSmallEventBlock::SampledIndels {
        lengths: SmallIndelLengths::default(),
        request: DeclaredSmallIndelRequest::new(
            ScenarioId::parse("direct", "test").expect("scenario"),
            reference.sequences()[0].id().clone(),
            HaplotypePhase::from_copy_index(0).expect("phase"),
            "global".to_owned(),
            1.0,
            1.0,
        ),
        counts: SampledSmallIndelCounts::new(insertion_count, deletion_count),
    }
}

fn insertion(
    reference: &Reference,
    point: u64,
    bases: &[u8],
    event_id: &str,
) -> DeclaredSmallEvent {
    DeclaredSmallEvent::Insertion(TargetedCurrentInsertion::new(
        reference.sequences()[0].id().clone(),
        HaplotypePhase::from_copy_index(0).expect("phase"),
        point,
        std::sync::Arc::<[u8]>::from(bases),
        ModificationProvenance::new(EventId::parse(event_id, "test").expect("event")),
    ))
}

fn deletion(reference: &Reference, start: u64, end: u64, event_id: &str) -> DeclaredSmallEvent {
    DeclaredSmallEvent::Deletion(TargetedCurrentDeletion::new(
        reference.sequences()[0].id().clone(),
        HaplotypePhase::from_copy_index(0).expect("phase"),
        start,
        end,
        ModificationProvenance::new(EventId::parse(event_id, "test").expect("event")),
    ))
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
fn declared_snp_then_deletion_removes_the_snp_without_publishing_a_prefix() {
    let reference = reference(b"ACG");
    let original = individual(&reference);
    let events = [
        DeclaredSmallEvent::Snp(snp(&reference, 1, CanonicalBase::T, "snp-first")),
        deletion(&reference, 1, 2, "delete-second-base"),
    ];

    preflight_declared_small_events(
        &original,
        &reference,
        RootSeed::from_bytes([42; 32]),
        &events,
        &[],
    )
    .expect("preflight");
    assert_eq!(original, individual(&reference));
    let result = apply_declared_small_events(
        &original,
        &reference,
        RootSeed::from_bytes([42; 32]),
        &events,
        &[],
    )
    .expect("composition");

    assert_eq!(first_copy_bases(result.individual(), &reference), b"AG");
    let locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            reference.sequences()[0].id().clone(),
            1,
            3,
        )
        .expect("position"),
        CanonicalBase::C,
        CanonicalBase::T,
    )
    .expect("locus");
    assert!(matches!(
        derive_snp_genotype(result.individual(), &reference, &locus),
        Err(GenotypeQueryError::LocusDeletedFromCopy { .. })
    ));
    assert!(matches!(
        result.applied_events(),
        [AppliedSmallEvent::Snp, AppliedSmallEvent::Deletion]
    ));
    assert!(result
        .individual()
        .copies_for(reference.sequences()[0].id())
        .expect("copies")[0]
        .segments()
        .current_modifications()
        .iter()
        .all(|piece| matches!(piece, genome_population_simulator::segments::CurrentIndelPiece::Source(source) if source.snps().is_empty())));
    assert_eq!(
        result
            .individual()
            .copies_for(reference.sequences()[0].id())
            .expect("copies")[0]
            .segments()
            .small_event_history()
            .len(),
        2
    );
}

#[test]
fn declared_insertions_are_composed_before_a_later_snp_on_the_same_copy() {
    let reference = reference(b"AC");
    let original = individual(&reference);
    let events = [
        insertion(&reference, 1, b"G", "insert-first"),
        insertion(&reference, 2, b"T", "insert-neighbor"),
        DeclaredSmallEvent::Snp(snp(&reference, 1, CanonicalBase::T, "snp-after-indels")),
    ];

    let result = apply_declared_small_events(
        &original,
        &reference,
        RootSeed::from_bytes([1; 32]),
        &events,
        &[],
    )
    .expect("composition");
    assert_eq!(first_copy_bases(result.individual(), &reference).len(), 4);
    let locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            reference.sequences()[0].id().clone(),
            1,
            2,
        )
        .expect("position"),
        CanonicalBase::C,
        CanonicalBase::T,
    )
    .expect("locus");
    assert_eq!(
        derive_snp_genotype(result.individual(), &reference, &locus)
            .expect("genotype")
            .alternate_dosage(),
        1
    );
    let current = result
        .individual()
        .copies_for(reference.sequences()[0].id())
        .expect("copies")[0]
        .segments()
        .current_modifications();
    assert_eq!(
        current
            .iter()
            .filter(|piece| matches!(
                piece,
                genome_population_simulator::segments::CurrentIndelPiece::Inserted(_)
            ))
            .count(),
        2
    );
    assert!(current.iter().any(
        |piece| matches!(piece, genome_population_simulator::segments::CurrentIndelPiece::Source(source) if source.snps().len() == 1)
    ));
}

#[test]
fn deletion_then_snp_is_an_indexed_explicit_conflict() {
    let reference = reference(b"ACG");
    let original = individual(&reference);
    let events = [
        deletion(&reference, 1, 2, "delete-first"),
        DeclaredSmallEvent::Snp(snp(&reference, 1, CanonicalBase::T, "snp-deleted")),
    ];
    assert_eq!(
        apply_declared_small_events(
            &original,
            &reference,
            RootSeed::from_bytes([0; 32]),
            &events,
            &[],
        ),
        Err(SmallEventCompositionError::Event {
            declared_index: 1,
            kind: DeclaredSmallEventKind::Snp,
            source: SmallEventSourceError::Snp(
                genome_population_simulator::snp_application::SnpApplicationError::Segment(
                    SegmentError::SnpPositionCoveredByDeletion,
                ),
            ),
        })
    );
    assert_eq!(original, individual(&reference));
}

#[test]
fn random_blocks_run_only_after_the_full_explicit_prefix() {
    let reference = reference(b"ACG");
    let original = individual(&reference);
    let explicit = [DeclaredSmallEvent::Snp(snp(
        &reference,
        1,
        CanonicalBase::T,
        "explicit-before-random",
    ))];
    let result = apply_declared_small_events(
        &original,
        &reference,
        RootSeed::from_bytes([42; 32]),
        &explicit,
        &[random_indels(&reference, 0, 1)],
    )
    .expect("the random deletion observes the explicit SNP state");
    assert_eq!(first_copy_bases(result.individual(), &reference), b"AG");
    assert!(matches!(
        result.applied_events(),
        [
            AppliedSmallEvent::Snp,
            AppliedSmallEvent::Indels {
                deletion_accepted: 1,
                ..
            }
        ]
    ));
}

#[test]
fn conflicting_or_incompatible_events_name_the_declared_element_and_leave_input_intact() {
    let reference = reference(b"ACG");
    let original = individual(&reference);
    let conflict = [
        DeclaredSmallEvent::Snp(snp(&reference, 1, CanonicalBase::T, "first")),
        DeclaredSmallEvent::Snp(snp(&reference, 1, CanonicalBase::A, "second")),
    ];
    assert_eq!(
        apply_declared_small_events(
            &original,
            &reference,
            RootSeed::from_bytes([0; 32]),
            &conflict,
            &[],
        ),
        Err(SmallEventCompositionError::Event {
            declared_index: 1,
            kind: DeclaredSmallEventKind::Snp,
            source: SmallEventSourceError::Snp(
                genome_population_simulator::snp_application::SnpApplicationError::Segment(
                    SegmentError::SnpPositionAlreadyModified,
                ),
            ),
        })
    );
    assert_eq!(original, individual(&reference));

    let foreign_reference = reference_with_id("ref-other", b"ACG");
    let incompatible = [DeclaredSmallEvent::Snp(snp(
        &foreign_reference,
        1,
        CanonicalBase::T,
        "foreign",
    ))];
    assert!(matches!(
        preflight_declared_small_events(
            &original,
            &reference,
            RootSeed::from_bytes([0; 32]),
            &incompatible,
            &[],
        ),
        Err(SmallEventCompositionError::Event {
            declared_index: 0,
            kind: DeclaredSmallEventKind::Snp,
            source: SmallEventSourceError::Snp(_),
        })
    ));
    assert_eq!(original, individual(&reference));
}

#[test]
fn short_chains_reject_an_explicit_deletion_during_preflight() {
    let reference = reference(b"AC");
    let original = individual(&reference);
    assert!(matches!(
        preflight_declared_small_events(
            &original,
            &reference,
            RootSeed::from_bytes([0; 32]),
            &[deletion(&reference, 1, 2, "too-short")],
            &[],
        ),
        Err(SmallEventCompositionError::Event {
            declared_index: 0,
            kind: DeclaredSmallEventKind::Deletion,
            source: SmallEventSourceError::Segment(_),
        })
    ));
    assert_eq!(original, individual(&reference));
}

fn reference_with_id(id: &str, bases: &[u8]) -> Reference {
    Reference::new(
        ReferenceId::parse(id, "test").expect("reference"),
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
