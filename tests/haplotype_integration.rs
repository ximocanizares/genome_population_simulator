//! Cierre integrador de F2.1 sobre las APIs públicas del estado haplotípico.
//!
//! La verdad esperada se escribe de forma literal para no reutilizar el
//! materializador como oráculo. F5.3 ampliará la ploidía efectiva por secuencia;
//! esta prueba cubre únicamente la ploidía uniforme provisional de F2.1.

use std::sync::Arc;

use genome_population_simulator::coordinates::{
    BasePosition, GenomicInterval, InsertionAnchor, InternalDeletion,
};
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::materialization::{
    SnpLocus, derive_snp_genotype, materialize_copy,
};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{
    CanonicalBase, DeletionModification, InsertionModification, ModificationProvenance, Segment,
    SegmentChain, SmallModification, SnpModification,
};

fn reference() -> Reference {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "integration test").expect("sequence ID"),
        FastaHeader::parse(b"chr1").expect("FASTA header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(b"NACGTN").expect("valid DNA with N"),
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("reference-integration", "integration test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn edited_chain(reference: &Reference) -> SegmentChain {
    let sequence_id = SequenceId::parse("chr1", "integration test").expect("sequence ID");
    let provenance = ModificationProvenance::new(
        EventId::parse("event-integration", "integration test").expect("event ID"),
    );
    let insertion = SmallModification::Insertion(
        InsertionModification::new(
            InsertionAnchor::from_base_1(reference.id().clone(), sequence_id.clone(), 2, 6)
                .expect("base-1 insertion anchor"),
            Arc::<[u8]>::from(&b"NN"[..]),
            provenance.clone(),
        )
        .expect("insertion"),
    );
    let snp = SmallModification::Snp(
        SnpModification::new(
            BasePosition::new(reference.id().clone(), sequence_id.clone(), 2, 6)
                .expect("base-0 SNP position"),
            CanonicalBase::C,
            CanonicalBase::T,
            provenance.clone(),
        )
        .expect("SNP"),
    );
    let deletion = SmallModification::Deletion(DeletionModification::new(
        InternalDeletion::new(
            GenomicInterval::new(reference.id().clone(), sequence_id.clone(), 3, 4, 6)
                .expect("base-0 half-open deletion"),
        )
        .expect("internal deletion"),
        provenance,
    ));
    SegmentChain::new(vec![
        Segment::new(
            GenomicInterval::new(reference.id().clone(), sequence_id, 0, 6, 6)
                .expect("base-0 half-open source"),
            vec![insertion, snp, deletion],
        )
        .expect("ordered non-overlapping edits"),
    ])
    .expect("single-source chain")
}

fn individual(reference: &Reference, ploidy: u8) -> Individual {
    Individual::from_reference(
        IndividualId::parse("individual-integration", "integration test").expect("individual ID"),
        reference,
        Ploidy::new(ploidy).expect("contract ploidy"),
    )
}

fn snp_locus(reference: &Reference) -> SnpLocus {
    SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            SequenceId::parse("chr1", "integration test").expect("sequence ID"),
            2,
            6,
        )
        .expect("base-0 locus"),
        CanonicalBase::C,
        CanonicalBase::T,
    )
    .expect("valid SNP locus")
}

#[test]
fn haplotype_state_reconstructs_known_truth_preserves_n_and_derives_phase_and_dose() {
    let reference = reference();
    let sequence_id = SequenceId::parse("chr1", "integration test").expect("sequence ID");
    let edited = edited_chain(&reference);

    let haploid = individual(&reference, 1)
        .with_copy_segments(
            &sequence_id,
            HaplotypePhase::from_copy_index(0).expect("phase"),
            edited.clone(),
        )
        .expect("edited haploid");
    assert_eq!(
        materialize_copy(
            &haploid.copies_for(&sequence_id).expect("copies")[0],
            &reference
        )
        .expect("haploid reconstruction"),
        b"NANNTTN"
    );

    let diploid = individual(&reference, 2)
        .with_copy_segments(
            &sequence_id,
            HaplotypePhase::from_copy_index(1).expect("phase"),
            edited.clone(),
        )
        .expect("edited diploid");
    let diploid_genotype = derive_snp_genotype(&diploid, &reference, &snp_locus(&reference))
        .expect("diploid genotype from canonical state");
    assert_eq!(diploid_genotype.alternate_dosage(), 1);
    assert_eq!(
        diploid_genotype
            .alleles()
            .iter()
            .map(|allele| (allele.phase().copy_index(), allele.allele()))
            .collect::<Vec<_>>(),
        vec![(0, CanonicalBase::C), (1, CanonicalBase::T)]
    );

    let polyploid = individual(&reference, 4)
        .with_copy_segments(
            &sequence_id,
            HaplotypePhase::from_copy_index(0).expect("phase"),
            edited.clone(),
        )
        .expect("first edited copy")
        .with_copy_segments(
            &sequence_id,
            HaplotypePhase::from_copy_index(3).expect("phase"),
            edited,
        )
        .expect("second edited copy");
    let polyploid_genotype = derive_snp_genotype(&polyploid, &reference, &snp_locus(&reference))
        .expect("polyploid genotype from canonical state");
    assert_eq!(polyploid_genotype.alternate_dosage(), 2);
    assert_eq!(
        polyploid_genotype
            .alleles()
            .iter()
            .map(|allele| allele.phase().copy_index())
            .collect::<Vec<_>>(),
        vec![0, 1, 2, 3]
    );

    let intact = individual(&reference, 4);
    let intact_copies = intact.copies_for(&sequence_id).expect("copies");
    assert!(
        intact_copies
            .windows(2)
            .all(|pair| pair[0].segments().is_shared_with(pair[1].segments()))
    );
    assert!(intact_copies.iter().all(|copy| {
        copy.segments().segments()[0].modifications().is_empty()
            && copy.segments().segments()[0].source().len() == 6
    }));
    assert_eq!(
        materialize_copy(&intact_copies[0], &reference).expect("intact reconstruction"),
        b"NACGTN"
    );
}
