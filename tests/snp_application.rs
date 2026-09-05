use genome_population_simulator::coordinates::{BasePosition, GenomicInterval, InternalDeletion};
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::materialization::{SnpLocus, derive_snp_genotype};
use genome_population_simulator::population_snvs::{
    TSV_HEADER, parse_canonical_tsv, prepare_population_snvs_from_reference,
};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{
    CanonicalBase, DeletionModification, ModificationProvenance, Segment, SegmentChain,
    SegmentError, SmallModification,
};
use genome_population_simulator::small_event_composition::{
    DeclaredSmallEvent, DeclaredSmallEventKind, SmallEventCompositionError, SmallEventSourceError,
    TargetedCurrentInsertion, apply_declared_small_events,
};
use genome_population_simulator::snp_application::{
    ApplicableSnp, SnpApplicationError, SnpOrigin, TargetedSnp, apply_snps, preflight_snps,
};

fn reference(bases: &[u8]) -> Reference {
    reference_with_id("ref-main", bases)
}

fn reference_with_id(id: &str, bases: &[u8]) -> Reference {
    let id = ReferenceId::parse(id, "test").expect("reference");
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence");
    let dna = DnaSequence::try_from_bytes(bases).expect("DNA");
    Reference::new(
        id,
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                sequence_id,
                FastaHeader::parse(b"chr1").expect("header"),
                SequenceClass::Chromosome,
                &dna,
            )
            .expect("sequence"),
        ],
    )
    .expect("reference")
}

fn two_sequence_reference() -> Reference {
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("first sequence ID"),
                FastaHeader::parse(b"chr1").expect("first header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(b"ACGT").expect("first DNA"),
            )
            .expect("first sequence"),
            ReferenceSequence::new(
                SequenceId::parse("chr2", "test").expect("second sequence ID"),
                FastaHeader::parse(b"chr2").expect("second header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(b"TGCA").expect("second DNA"),
            )
            .expect("second sequence"),
        ],
    )
    .expect("reference")
}

fn explicit_on_sequence(
    reference: &Reference,
    sequence_id: &SequenceId,
    offset: u64,
    alternate: CanonicalBase,
    event: &str,
) -> ApplicableSnp {
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id() == sequence_id)
        .expect("known sequence");
    let reference_base = CanonicalBase::from_byte(
        sequence
            .sequence()
            .as_bytes()
            .expect("bases materializadas")[usize::try_from(offset).expect("offset")],
    )
    .expect("canonical reference");
    ApplicableSnp::explicit(
        SnpLocus::new(
            BasePosition::new(
                reference.id().clone(),
                sequence_id.clone(),
                offset,
                sequence.sequence().len(),
            )
            .expect("position"),
            reference_base,
            alternate,
        )
        .expect("locus"),
        EventId::parse(event, "test").expect("event"),
    )
}

fn individual(reference: &Reference, id: &str) -> Individual {
    Individual::from_reference(
        IndividualId::parse(id, "test").expect("individual"),
        reference,
        Ploidy::new(2).expect("ploidy"),
    )
}

fn explicit(
    reference: &Reference,
    offset: u64,
    alternate: CanonicalBase,
    event: &str,
) -> ApplicableSnp {
    let sequence = &reference.sequences()[0];
    let reference_base = CanonicalBase::from_byte(
        sequence
            .sequence()
            .as_bytes()
            .expect("bases materializadas")[usize::try_from(offset).expect("offset")],
    )
    .expect("canonical reference");
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
        EventId::parse(event, "test").expect("event"),
    )
}

fn full_interval(reference: &Reference) -> GenomicInterval {
    let sequence = &reference.sequences()[0];
    GenomicInterval::new(
        reference.id().clone(),
        sequence.id().clone(),
        0,
        sequence.sequence().len(),
        sequence.sequence().len(),
    )
    .expect("full interval")
}

fn with_phase_chain(
    individual: &Individual,
    reference: &Reference,
    phase: HaplotypePhase,
    chain: SegmentChain,
) -> Individual {
    individual
        .with_copy_segments(reference.sequences()[0].id(), phase, chain)
        .expect("updated copy")
}

#[test]
fn explicit_snp_preflights_then_updates_only_the_target_phase() {
    let reference = reference(b"ACGT");
    let original = individual(&reference, "individual-1");
    let application = TargetedSnp::new(
        HaplotypePhase::from_copy_index(1).expect("phase"),
        explicit(&reference, 1, CanonicalBase::T, "event-explicit"),
    );

    preflight_snps(&original, &reference, std::slice::from_ref(&application)).expect("preflight");
    let updated = apply_snps(&original, &reference, &[application]).expect("application");
    assert_eq!(original, individual(&reference, "individual-1"));

    let locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            SequenceId::parse("chr1", "test").expect("sequence"),
            1,
            4,
        )
        .expect("position"),
        CanonicalBase::C,
        CanonicalBase::T,
    )
    .expect("locus");
    let genotype = derive_snp_genotype(&updated, &reference, &locus).expect("genotype");
    assert_eq!(genotype.alternate_dosage(), 1);
    assert_eq!(genotype.alleles()[1].phase().copy_index(), 1);
}

#[test]
fn matrix_snp_preserves_canonical_variant_and_population_provenance() {
    let reference = reference(b"ACGT");
    let header = TSV_HEADER.join("\t");
    let table = parse_canonical_tsv(
        format!("{header}\nignored\tchr1\t2\tC\tT\tP1\t0.25\tpanel-row\n").as_bytes(),
    )
    .expect("matrix TSV");
    let prepared = prepare_population_snvs_from_reference(&table, &reference, 1e-9, false)
        .expect("prepared matrix");
    let row = &prepared.rows()[0];
    let application = ApplicableSnp::from_population_matrix(
        &reference,
        row,
        EventId::parse("event-matrix", "test").expect("event"),
    )
    .expect("matrix application");
    assert!(matches!(
        application.origin(),
        SnpOrigin::PopulationMatrix { variant_id, population_id, source_id }
            if variant_id == row.variant_id() && population_id == row.population_id() && source_id == row.source_id()
    ));

    let updated = apply_snps(
        &individual(&reference, "individual-1"),
        &reference,
        &[TargetedSnp::new(
            HaplotypePhase::from_copy_index(0).expect("phase"),
            application,
        )],
    )
    .expect("application");
    let modification = updated.copies_for(row.contig()).expect("copies")[0]
        .segments()
        .segments()[0]
        .modifications()
        .first()
        .expect("SNP");
    let SmallModification::Snp(snp) = modification else {
        panic!("expected SNP");
    };
    assert_eq!(snp.reference(), CanonicalBase::C);
    assert_eq!(snp.alternate(), CanonicalBase::T);
    assert_eq!(snp.matrix_variant_id(), Some(row.variant_id()));
    assert_eq!(snp.matrix_population_id(), Some(row.population_id()));
    assert_eq!(snp.matrix_source_id(), Some("panel-row"));
    assert_eq!(snp.provenance().event_id().as_str(), "event-matrix");
}

#[test]
fn preflight_rejects_conflicts_without_partially_updating_a_copy() {
    let reference = reference(b"ACGT");
    let original = individual(&reference, "individual-1");
    let phase = HaplotypePhase::from_copy_index(0).expect("phase");
    let first = TargetedSnp::new(phase, explicit(&reference, 1, CanonicalBase::T, "event-1"));
    let duplicate = TargetedSnp::new(phase, explicit(&reference, 1, CanonicalBase::A, "event-2"));

    assert_eq!(
        preflight_snps(&original, &reference, &[first.clone(), duplicate.clone()]),
        Err(SnpApplicationError::DuplicateTarget { phase, position: 1 })
    );
    assert_eq!(
        apply_snps(&original, &reference, &[first, duplicate]),
        Err(SnpApplicationError::DuplicateTarget { phase, position: 1 })
    );
    assert_eq!(original, individual(&reference, "individual-1"));
}

#[test]
fn preflight_rejects_an_explicit_snp_on_unknown_reference_base() {
    let reference = reference(b"ANGT");
    let original = individual(&reference, "individual-1");
    let locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            SequenceId::parse("chr1", "test").expect("sequence"),
            1,
            4,
        )
        .expect("position"),
        CanonicalBase::A,
        CanonicalBase::T,
    )
    .expect("locus");
    let application = TargetedSnp::new(
        HaplotypePhase::from_copy_index(0).expect("phase"),
        ApplicableSnp::explicit(locus, EventId::parse("event-n", "test").expect("event")),
    );

    assert_eq!(
        preflight_snps(&original, &reference, &[application]),
        Err(SnpApplicationError::ReferenceBaseUnknown { position: 1 })
    );
}

#[test]
fn explicit_snp_rejects_a_declared_ref_that_disagrees_with_the_authoritative_reference() {
    let reference = reference(b"ACGT");
    let original = individual(&reference, "individual-1");
    let locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            SequenceId::parse("chr1", "test").expect("sequence"),
            1,
            4,
        )
        .expect("position"),
        CanonicalBase::A,
        CanonicalBase::T,
    )
    .expect("locus");
    let application = TargetedSnp::new(
        HaplotypePhase::from_copy_index(0).expect("phase"),
        ApplicableSnp::explicit(
            locus,
            EventId::parse("event-explicit-ref-mismatch", "test").expect("event"),
        ),
    );

    assert_eq!(
        apply_snps(&original, &reference, &[application]),
        Err(SnpApplicationError::ReferenceBaseMismatch {
            position: 1,
            expected: CanonicalBase::A,
            found: b'C',
        })
    );
    assert_eq!(original, individual(&reference, "individual-1"));
}

#[test]
fn matrix_rejects_a_declared_ref_that_disagrees_with_the_authoritative_fasta() {
    let reference = reference(b"ACGT");
    let original = individual(&reference, "individual-1");
    let header = TSV_HEADER.join("\t");
    let table = parse_canonical_tsv(
        format!("{header}\nignored\tchr1\t2\tA\tT\tP1\t0.25\tsource\n").as_bytes(),
    )
    .expect("matrix TSV");

    assert!(prepare_population_snvs_from_reference(&table, &reference, 1e-9, false).is_err());
    assert_eq!(original, individual(&reference, "individual-1"));
}

#[test]
fn preflight_accepts_first_and_last_base_without_mutating_the_input() {
    let reference = reference(b"ACGT");
    let original = individual(&reference, "individual-1");
    let phase = HaplotypePhase::from_copy_index(0).expect("phase");
    let updated = apply_snps(
        &original,
        &reference,
        &[
            TargetedSnp::new(
                phase,
                explicit(&reference, 0, CanonicalBase::C, "event-first"),
            ),
            TargetedSnp::new(
                phase,
                explicit(&reference, 3, CanonicalBase::A, "event-last"),
            ),
        ],
    )
    .expect("terminal SNPs are valid");

    assert_eq!(original, individual(&reference, "individual-1"));
    assert_eq!(
        updated
            .copies_for(reference.sequences()[0].id())
            .expect("copies")[0]
            .segments()
            .segments()[0]
            .modifications()
            .len(),
        2
    );
}

#[test]
fn preflight_rejects_missing_copy_and_preserves_the_input() {
    let reference = reference(b"ACGT");
    let original = individual(&reference, "individual-1");
    let missing_phase = HaplotypePhase::from_copy_index(2).expect("phase domain");
    let application = TargetedSnp::new(
        missing_phase,
        explicit(&reference, 1, CanonicalBase::T, "event-missing-copy"),
    );

    assert_eq!(
        apply_snps(&original, &reference, &[application]),
        Err(SnpApplicationError::MissingTargetCopy {
            phase: missing_phase
        })
    );
    assert_eq!(original, individual(&reference, "individual-1"));
}

#[test]
fn preflight_rejects_incompatible_individual_and_locus_references_transactionally() {
    let reference = reference(b"ACGT");
    let other_reference = reference_with_id("ref-other", b"ACGT");
    let original = individual(&reference, "individual-1");
    let phase = HaplotypePhase::from_copy_index(0).expect("phase");
    let application = TargetedSnp::new(
        phase,
        explicit(&reference, 1, CanonicalBase::T, "event-reference"),
    );
    assert_eq!(
        apply_snps(
            &original,
            &other_reference,
            std::slice::from_ref(&application)
        ),
        Err(SnpApplicationError::IndividualReferenceMismatch)
    );
    assert_eq!(original, individual(&reference, "individual-1"));

    let foreign_locus = SnpLocus::new(
        BasePosition::new(
            ReferenceId::parse("ref-other", "test").expect("reference"),
            SequenceId::parse("chr1", "test").expect("sequence"),
            1,
            4,
        )
        .expect("position"),
        CanonicalBase::C,
        CanonicalBase::T,
    )
    .expect("locus");
    let foreign_application = TargetedSnp::new(
        phase,
        ApplicableSnp::explicit(
            foreign_locus,
            EventId::parse("event-foreign-locus", "test").expect("event"),
        ),
    );
    assert_eq!(
        apply_snps(&original, &reference, &[foreign_application]),
        Err(SnpApplicationError::LocusReferenceMismatch)
    );
    assert_eq!(original, individual(&reference, "individual-1"));
}

#[test]
fn preflight_rejects_missing_ambiguous_and_conflicting_loci_transactionally() {
    let reference = reference(b"ACGT");
    let phase = HaplotypePhase::from_copy_index(0).expect("phase");
    let application = TargetedSnp::new(phase, explicit(&reference, 1, CanonicalBase::T, "event"));
    let original = individual(&reference, "individual-1");

    let missing = with_phase_chain(
        &original,
        &reference,
        phase,
        SegmentChain::new(vec![
            Segment::new(
                GenomicInterval::new(
                    reference.id().clone(),
                    reference.sequences()[0].id().clone(),
                    0,
                    1,
                    4,
                )
                .expect("first segment"),
                vec![],
            )
            .expect("segment"),
            Segment::new(
                GenomicInterval::new(
                    reference.id().clone(),
                    reference.sequences()[0].id().clone(),
                    2,
                    4,
                    4,
                )
                .expect("second segment"),
                vec![],
            )
            .expect("segment"),
        ])
        .expect("chain"),
    );
    assert_eq!(
        apply_snps(&missing, &reference, std::slice::from_ref(&application)),
        Err(SnpApplicationError::Segment(
            SegmentError::SnpPositionNotCovered
        ))
    );
    assert_eq!(
        missing,
        with_phase_chain(
            &original,
            &reference,
            phase,
            SegmentChain::new(vec![
                Segment::new(
                    GenomicInterval::new(
                        reference.id().clone(),
                        reference.sequences()[0].id().clone(),
                        0,
                        1,
                        4
                    )
                    .expect("first"),
                    vec![]
                )
                .expect("segment"),
                Segment::new(
                    GenomicInterval::new(
                        reference.id().clone(),
                        reference.sequences()[0].id().clone(),
                        2,
                        4,
                        4
                    )
                    .expect("second"),
                    vec![]
                )
                .expect("segment"),
            ])
            .expect("chain"),
        )
    );

    let ambiguous = with_phase_chain(
        &original,
        &reference,
        phase,
        SegmentChain::new(vec![
            Segment::new(
                GenomicInterval::new(
                    reference.id().clone(),
                    reference.sequences()[0].id().clone(),
                    0,
                    3,
                    4,
                )
                .expect("first"),
                vec![],
            )
            .expect("segment"),
            Segment::new(
                GenomicInterval::new(
                    reference.id().clone(),
                    reference.sequences()[0].id().clone(),
                    1,
                    4,
                    4,
                )
                .expect("second"),
                vec![],
            )
            .expect("segment"),
        ])
        .expect("chain"),
    );
    assert_eq!(
        apply_snps(&ambiguous, &reference, std::slice::from_ref(&application)),
        Err(SnpApplicationError::Segment(
            SegmentError::AmbiguousSnpPosition
        ))
    );
    assert_eq!(
        ambiguous,
        with_phase_chain(
            &original,
            &reference,
            phase,
            SegmentChain::new(vec![
                Segment::new(
                    GenomicInterval::new(
                        reference.id().clone(),
                        reference.sequences()[0].id().clone(),
                        0,
                        3,
                        4
                    )
                    .expect("first"),
                    vec![]
                )
                .expect("segment"),
                Segment::new(
                    GenomicInterval::new(
                        reference.id().clone(),
                        reference.sequences()[0].id().clone(),
                        1,
                        4,
                        4
                    )
                    .expect("second"),
                    vec![]
                )
                .expect("segment"),
            ])
            .expect("chain"),
        )
    );

    let prior_snp =
        apply_snps(&original, &reference, std::slice::from_ref(&application)).expect("first SNP");
    assert_eq!(
        apply_snps(&prior_snp, &reference, std::slice::from_ref(&application)),
        Err(SnpApplicationError::Segment(
            SegmentError::SnpPositionAlreadyModified
        ))
    );
    assert_eq!(
        prior_snp,
        apply_snps(&original, &reference, std::slice::from_ref(&application)).expect("first SNP")
    );

    let deletion = DeletionModification::new(
        InternalDeletion::new(
            GenomicInterval::new(
                reference.id().clone(),
                reference.sequences()[0].id().clone(),
                1,
                2,
                4,
            )
            .expect("deletion interval"),
        )
        .expect("internal deletion"),
        ModificationProvenance::new(EventId::parse("deletion", "test").expect("event")),
    );
    let deleted = with_phase_chain(
        &original,
        &reference,
        phase,
        SegmentChain::new(vec![
            Segment::new(
                full_interval(&reference),
                vec![SmallModification::Deletion(deletion)],
            )
            .expect("deleted segment"),
        ])
        .expect("chain"),
    );
    assert_eq!(
        apply_snps(&deleted, &reference, &[application]),
        Err(SnpApplicationError::Segment(
            SegmentError::SnpPositionCoveredByDeletion
        ))
    );
    assert_eq!(
        deleted,
        with_phase_chain(
            &original,
            &reference,
            phase,
            SegmentChain::new(vec![
                Segment::new(
                    full_interval(&reference),
                    vec![SmallModification::Deletion(DeletionModification::new(
                        InternalDeletion::new(
                            GenomicInterval::new(
                                reference.id().clone(),
                                reference.sequences()[0].id().clone(),
                                1,
                                2,
                                4
                            )
                            .expect("interval")
                        )
                        .expect("deletion"),
                        ModificationProvenance::new(
                            EventId::parse("deletion", "test").expect("event")
                        ),
                    ))],
                )
                .expect("segment"),
            ])
            .expect("chain"),
        )
    );
}

#[test]
fn one_canonical_matrix_row_can_target_fixed_individuals_and_phases_without_rng() {
    let reference = reference(b"ACGT");
    let header = TSV_HEADER.join("\t");
    let table = parse_canonical_tsv(
        format!("{header}\nignored\tchr1\t2\tC\tT\tP1\t0.25\tpanel-row\n").as_bytes(),
    )
    .expect("matrix TSV");
    let prepared = prepare_population_snvs_from_reference(&table, &reference, 1e-9, false)
        .expect("prepared matrix");
    let row = &prepared.rows()[0];
    let source = ApplicableSnp::from_population_matrix(
        &reference,
        row,
        EventId::parse("matrix-event", "test").expect("event"),
    )
    .expect("matrix application");
    let first = apply_snps(
        &individual(&reference, "individual-1"),
        &reference,
        &[TargetedSnp::new(
            HaplotypePhase::from_copy_index(0).expect("phase"),
            source.clone(),
        )],
    )
    .expect("first application");
    let second = apply_snps(
        &individual(&reference, "individual-2"),
        &reference,
        &[TargetedSnp::new(
            HaplotypePhase::from_copy_index(1).expect("phase"),
            source,
        )],
    )
    .expect("second application");
    for (individual, phase) in [(&first, 0_usize), (&second, 1_usize)] {
        let modification = individual.copies_for(row.contig()).expect("copies")[phase]
            .segments()
            .segments()[0]
            .modifications()[0]
            .clone();
        let SmallModification::Snp(snp) = modification else {
            panic!("SNP")
        };
        assert_eq!(snp.matrix_variant_id(), Some(row.variant_id()));
        assert_eq!(snp.matrix_population_id(), Some(row.population_id()));
        assert_eq!(snp.matrix_source_id(), Some(row.source_id()));
    }
}

#[test]
fn matrix_snp_batch_matches_sequential_application_across_sequences_and_phases() {
    let reference = two_sequence_reference();
    let original = individual(&reference, "individual-1");
    let header = TSV_HEADER.join("\t");
    let table = parse_canonical_tsv(
        format!(
            "{header}\nignored\tchr1\t1\tA\tC\tP1\t0.5\tmatrix-1\nignored\tchr1\t4\tT\tA\tP1\t0.5\tmatrix-2\nignored\tchr2\t2\tG\tA\tP1\t0.5\tmatrix-3\n"
        )
        .as_bytes(),
    )
    .expect("matrix TSV");
    let prepared = prepare_population_snvs_from_reference(&table, &reference, 1e-9, false)
        .expect("prepared matrix");
    let phase_zero = HaplotypePhase::from_copy_index(0).expect("phase zero");
    let phase_one = HaplotypePhase::from_copy_index(1).expect("phase one");
    let applications = prepared
        .rows()
        .iter()
        .zip([phase_zero, phase_one, phase_zero])
        .enumerate()
        .map(|(index, (row, phase))| {
            TargetedSnp::new(
                phase,
                ApplicableSnp::from_population_matrix(
                    &reference,
                    row,
                    EventId::parse(format!("matrix-event-{}", index + 1), "test")
                        .expect("matrix event"),
                )
                .expect("matrix SNP"),
            )
        })
        .collect::<Vec<_>>();

    let batched = apply_snps(&original, &reference, &applications).expect("matrix batch");
    let repeated = apply_snps(&original, &reference, &applications).expect("repeated matrix batch");
    let sequential = applications
        .iter()
        .try_fold(original.clone(), |current, application| {
            apply_snps(&current, &reference, std::slice::from_ref(application))
        })
        .expect("sequential matrix SNPs");
    assert_eq!(batched, sequential);
    assert_eq!(batched, repeated);
    assert_eq!(original, individual(&reference, "individual-1"));

    let duplicate_matrix_target = [
        DeclaredSmallEvent::Snp(applications[0].clone()),
        DeclaredSmallEvent::Snp(applications[0].clone()),
    ];
    assert_eq!(
        apply_declared_small_events(
            &original,
            &reference,
            genome_population_simulator::randomness::RootSeed::from_bytes([7; 32]),
            &duplicate_matrix_target,
            &[],
        ),
        Err(SmallEventCompositionError::Event {
            declared_index: 1,
            kind: DeclaredSmallEventKind::Snp,
            source: SmallEventSourceError::Snp(SnpApplicationError::Segment(
                SegmentError::SnpPositionAlreadyModified,
            )),
        })
    );
    assert_eq!(original, individual(&reference, "individual-1"));

    let chr1 = SequenceId::parse("chr1", "test").expect("chr1 ID");
    let later_events = [
        DeclaredSmallEvent::Snp(TargetedSnp::new(
            phase_zero,
            explicit_on_sequence(
                &reference,
                &chr1,
                2,
                CanonicalBase::T,
                "explicit-after-matrix",
            ),
        )),
        DeclaredSmallEvent::Insertion(TargetedCurrentInsertion::new(
            chr1.clone(),
            phase_zero,
            2,
            b"GG".as_slice(),
            ModificationProvenance::new(
                EventId::parse("indel-after-matrix", "test").expect("indel event"),
            ),
        )),
    ];
    let root_seed = genome_population_simulator::randomness::RootSeed::from_bytes([7; 32]);
    let after_batched =
        apply_declared_small_events(&batched, &reference, root_seed, &later_events, &[])
            .expect("later explicit events after batch");
    let after_sequential =
        apply_declared_small_events(&sequential, &reference, root_seed, &later_events, &[])
            .expect("later explicit events after sequential SNPs");
    let after_repeated =
        apply_declared_small_events(&batched, &reference, root_seed, &later_events, &[])
            .expect("repeated later explicit events");

    assert_eq!(after_batched, after_sequential);
    assert_eq!(after_batched, after_repeated);
    for sequence in reference.sequences() {
        for phase in [phase_zero, phase_one] {
            let batched_copy = after_batched
                .individual()
                .copies_for(sequence.id())
                .expect("batched copies")
                .iter()
                .find(|copy| copy.phase() == phase)
                .expect("batched phase");
            let sequential_copy = after_sequential
                .individual()
                .copies_for(sequence.id())
                .expect("sequential copies")
                .iter()
                .find(|copy| copy.phase() == phase)
                .expect("sequential phase");
            assert_eq!(
                genome_population_simulator::materialization::materialize_copy(
                    batched_copy,
                    &reference
                )
                .expect("batched genome"),
                genome_population_simulator::materialization::materialize_copy(
                    sequential_copy,
                    &reference
                )
                .expect("sequential genome")
            );
            assert_eq!(
                batched_copy.segments().current_modifications(),
                sequential_copy.segments().current_modifications()
            );
            assert_eq!(
                batched_copy.segments().small_event_history(),
                sequential_copy.segments().small_event_history()
            );
        }
    }
}
