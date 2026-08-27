use std::sync::Arc;

use genome_population_simulator::coordinates::{
    BasePosition, GenomicInterval, InsertionAnchor, InternalDeletion,
};
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::materialization::{
    GenotypeQueryError, MaterializationError, SnpLocus, SnpLocusError, derive_snp_genotype,
    materialize_copy, materialize_copy_interval,
};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{
    CanonicalBase, DeletionModification, InsertionModification, ModificationProvenance, Segment,
    SegmentChain, SmallModification, SnpModification,
};

fn reference_with_bases(id: &str, bases: &[u8]) -> Reference {
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let sequence = ReferenceSequence::new(
        sequence_id,
        FastaHeader::parse(b"chr1").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(bases).expect("DNA"),
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse(id, "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn reference() -> Reference {
    reference_with_bases("reference-main", b"NACGTN")
}

fn provenance() -> ModificationProvenance {
    ModificationProvenance::new(EventId::parse("event-1", "test").expect("event ID"))
}

fn full_interval(reference: &Reference) -> GenomicInterval {
    GenomicInterval::new(
        reference.id().clone(),
        SequenceId::parse("chr1", "test").expect("sequence ID"),
        0,
        6,
        6,
    )
    .expect("full interval")
}

fn edited_chain(reference: &Reference) -> SegmentChain {
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let insertion = SmallModification::Insertion(
        InsertionModification::new(
            InsertionAnchor::from_base_1(reference.id().clone(), sequence_id.clone(), 2, 6)
                .expect("anchor"),
            Arc::<[u8]>::from(&b"NN"[..]),
            provenance(),
        )
        .expect("insertion"),
    );
    let snp = SmallModification::Snp(
        SnpModification::new(
            BasePosition::new(reference.id().clone(), sequence_id.clone(), 2, 6)
                .expect("SNP position"),
            CanonicalBase::C,
            CanonicalBase::T,
            provenance(),
        )
        .expect("SNP"),
    );
    let deletion = SmallModification::Deletion(DeletionModification::new(
        InternalDeletion::new(
            GenomicInterval::new(reference.id().clone(), sequence_id, 3, 4, 6)
                .expect("deletion interval"),
        )
        .expect("internal deletion"),
        provenance(),
    ));
    SegmentChain::new(vec![
        Segment::new(full_interval(reference), vec![insertion, snp, deletion])
            .expect("ordered edits"),
    ])
    .expect("chain")
}

fn individual(reference: &Reference, ploidy: u8) -> Individual {
    Individual::from_reference(
        IndividualId::parse("individual-1", "test").expect("individual ID"),
        reference,
        Ploidy::new(ploidy).expect("ploidy"),
    )
}

fn locus(reference: &Reference) -> SnpLocus {
    SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            SequenceId::parse("chr1", "test").expect("sequence ID"),
            2,
            6,
        )
        .expect("position"),
        CanonicalBase::C,
        CanonicalBase::T,
    )
    .expect("valid SNP locus")
}

#[test]
fn materializes_a_haploid_copy_from_manual_snp_indel_truth_and_preserves_n() {
    let reference = reference();
    let individual = individual(&reference, 1)
        .with_copy_segments(
            &SequenceId::parse("chr1", "test").expect("sequence ID"),
            HaplotypePhase::from_copy_index(0).expect("phase"),
            edited_chain(&reference),
        )
        .expect("updated individual");

    let sequence = materialize_copy(&individual.sequences()[0].copies()[0], &reference)
        .expect("materialized copy");
    assert_eq!(sequence, b"NANNTTN");
}

#[test]
fn interval_materialization_matches_the_manual_full_truth_across_edits_and_edges() {
    let reference = reference();
    let individual = individual(&reference, 1)
        .with_copy_segments(
            &SequenceId::parse("chr1", "test").expect("sequence ID"),
            HaplotypePhase::from_copy_index(0).expect("phase"),
            edited_chain(&reference),
        )
        .expect("updated individual");
    let copy = &individual.sequences()[0].copies()[0];
    // Verdad construida manualmente: NACGTN, inserción NN tras A, C→T y
    // deleción de G = NANNTTN. No deriva las ventanas del materializador.
    let truth = b"NANNTTN";
    let full = materialize_copy(copy, &reference).expect("full copy");
    assert_eq!(full, truth);

    for (start, end, label) in [
        (0, 2, "extremo inicial"),
        (1, 5, "inserción y SNP"),
        (3, 6, "límite de deleción"),
        (5, 7, "extremo final"),
    ] {
        assert_eq!(
            materialize_copy_interval(copy, &reference, start, end).expect(label),
            full[usize::try_from(start).expect("start")..usize::try_from(end).expect("end")],
            "{label}",
        );
    }
    for start in 0..u64::try_from(truth.len()).expect("length") {
        for end in start + 1..=u64::try_from(truth.len()).expect("length") {
            assert_eq!(
                materialize_copy_interval(copy, &reference, start, end).expect("window"),
                full[usize::try_from(start).expect("start")..usize::try_from(end).expect("end")],
            );
        }
    }
}

#[test]
fn derives_diploid_phase_in_copy_order_and_dose_from_canonical_edits() {
    let reference = reference();
    let individual = individual(&reference, 2)
        .with_copy_segments(
            &SequenceId::parse("chr1", "test").expect("sequence ID"),
            HaplotypePhase::from_copy_index(1).expect("phase"),
            edited_chain(&reference),
        )
        .expect("updated individual");

    let genotype =
        derive_snp_genotype(&individual, &reference, &locus(&reference)).expect("phased genotype");
    assert_eq!(genotype.alternate_dosage(), 1);
    assert_eq!(
        genotype
            .alleles()
            .iter()
            .map(|allele| (allele.phase().copy_index(), allele.allele()))
            .collect::<Vec<_>>(),
        vec![(0, CanonicalBase::C), (1, CanonicalBase::T)]
    );
}

#[test]
fn derives_polyploid_dose_as_the_number_of_carrier_copies() {
    let reference = reference();
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let individual = individual(&reference, 3)
        .with_copy_segments(
            &sequence_id,
            HaplotypePhase::from_copy_index(0).expect("phase"),
            edited_chain(&reference),
        )
        .expect("first edit")
        .with_copy_segments(
            &sequence_id,
            HaplotypePhase::from_copy_index(2).expect("phase"),
            edited_chain(&reference),
        )
        .expect("second edit");

    let genotype =
        derive_snp_genotype(&individual, &reference, &locus(&reference)).expect("phased genotype");
    assert_eq!(genotype.alternate_dosage(), 2);
    assert_eq!(
        genotype
            .alleles()
            .iter()
            .map(|allele| allele.phase().copy_index())
            .collect::<Vec<_>>(),
        vec![0, 1, 2]
    );
}

#[test]
fn rejects_reference_and_sequence_incompatibilities_instead_of_reinterpreting_coordinates() {
    let reference = reference();
    let individual = individual(&reference, 1);
    assert_eq!(
        materialize_copy(
            &individual.sequences()[0].copies()[0],
            &reference_with_bases("other-reference", b"NACGTN"),
        ),
        Err(MaterializationError::CopyReferenceMismatch)
    );
    assert_eq!(
        materialize_copy(
            &individual.sequences()[0].copies()[0],
            &reference_with_bases("reference-main", b"NACGTNN"),
        ),
        Err(MaterializationError::SourceLengthMismatch {
            expected: 6,
            found: 7,
        })
    );
}

#[test]
fn rejects_a_locus_deleted_from_one_copy() {
    let reference = reference();
    let individual = individual(&reference, 1)
        .with_copy_segments(
            &SequenceId::parse("chr1", "test").expect("sequence ID"),
            HaplotypePhase::from_copy_index(0).expect("phase"),
            edited_chain(&reference),
        )
        .expect("updated individual");
    let deleted_locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            SequenceId::parse("chr1", "test").expect("sequence ID"),
            3,
            6,
        )
        .expect("position"),
        CanonicalBase::G,
        CanonicalBase::A,
    )
    .expect("valid SNP locus");
    assert_eq!(
        derive_snp_genotype(&individual, &reference, &deleted_locus),
        Err(GenotypeQueryError::LocusDeletedFromCopy {
            phase: HaplotypePhase::from_copy_index(0).expect("phase"),
        })
    );
}

#[test]
fn rejects_equal_snp_alleles_before_any_genotype_or_dose_is_derived() {
    let reference = reference();
    let position = BasePosition::new(
        reference.id().clone(),
        SequenceId::parse("chr1", "test").expect("sequence ID"),
        2,
        6,
    )
    .expect("position");
    assert_eq!(
        SnpLocus::new(position, CanonicalBase::C, CanonicalBase::C),
        Err(SnpLocusError::AllelesEqual)
    );
}

#[test]
fn rejects_a_chain_with_the_right_ids_but_the_wrong_source_length() {
    let reference = reference();
    let wrong_length_chain = SegmentChain::new(vec![
        Segment::new(
            GenomicInterval::new(
                reference.id().clone(),
                SequenceId::parse("chr1", "test").expect("sequence ID"),
                0,
                6,
                7,
            )
            .expect("interval in incompatible source space"),
            vec![],
        )
        .expect("segment"),
    ])
    .expect("chain");
    assert_eq!(
        individual(&reference, 1).with_copy_segments(
            &SequenceId::parse("chr1", "test").expect("sequence ID"),
            HaplotypePhase::from_copy_index(0).expect("phase"),
            wrong_length_chain,
        ),
        Err(genome_population_simulator::haplotype::HaplotypeStateError::SegmentSourceLengthMismatch {
            expected: 6,
            found: 7,
        })
    );
}
