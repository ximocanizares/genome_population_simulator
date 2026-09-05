use std::collections::BTreeSet;

use genome_population_simulator::haplotype::{Individual, Ploidy, PloidyError};
use genome_population_simulator::identity::{IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};

fn reference() -> Reference {
    let sequences = [
        ("chr1", b"NACGTN".as_slice(), SequenceClass::Chromosome),
        ("scaffold1", b"N".as_slice(), SequenceClass::Scaffold),
    ]
    .into_iter()
    .map(|(id, bases, class)| {
        ReferenceSequence::new(
            SequenceId::parse(id, "reference").expect("sequence id"),
            FastaHeader::parse(id.as_bytes()).expect("header"),
            class,
            &DnaSequence::try_from_bytes(bases).expect("DNA"),
        )
        .expect("reference sequence")
    })
    .collect();
    Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference id"),
        ReferenceOrigin::Generated,
        sequences,
    )
    .expect("reference")
}

#[test]
fn ploidy_accepts_every_contract_boundary_and_rejects_values_outside_it() {
    for value in [1, 2, 10] {
        assert_eq!(
            Ploidy::new(value).expect("valid ploidy").copy_count(),
            value
        );
    }
    assert_eq!(Ploidy::new(0), Err(PloidyError::OutOfRange { value: 0 }));
    assert_eq!(Ploidy::new(11), Err(PloidyError::OutOfRange { value: 11 }));
}

#[test]
fn individual_has_every_sequence_and_one_distinct_phased_copy_per_ploidy_index() {
    let reference = reference();
    let individual = Individual::from_reference(
        IndividualId::parse("individual-1", "test").expect("individual id"),
        &reference,
        Ploidy::new(10).expect("ploidy"),
    );

    assert_eq!(individual.reference_id(), reference.id());
    assert_eq!(individual.sequences().len(), 2);
    let copy_ids: BTreeSet<_> = individual
        .sequences()
        .iter()
        .flat_map(|sequence| sequence.copies())
        .map(|copy| copy.id().clone())
        .collect();
    assert_eq!(copy_ids.len(), 20);
    for sequence in individual.sequences() {
        assert_eq!(sequence.copies().len(), 10);
        for (expected, copy) in sequence.copies().iter().enumerate() {
            assert_eq!(copy.id().reference_id(), reference.id());
            assert_eq!(copy.id().individual_id(), individual.id());
            assert_eq!(copy.id().sequence_id(), sequence.sequence_id());
            assert_eq!(copy.id().copy_index(), expected as u8);
            assert_eq!(copy.phase().copy_index(), expected as u8);
        }
    }
}

#[test]
fn building_copies_keeps_unknown_bases_only_in_the_shared_reference() {
    let reference = reference();
    let before = reference.sequences()[0]
        .sequence()
        .as_bytes()
        .expect("bases materializadas")
        .to_vec();
    let individual = Individual::from_reference(
        IndividualId::parse("individual-n", "test").expect("individual id"),
        &reference,
        Ploidy::new(2).expect("ploidy"),
    );

    assert_eq!(
        reference.sequences()[0]
            .sequence()
            .as_bytes()
            .expect("bases materializadas"),
        before
    );
    assert_eq!(reference.sequences()[0].unknown_base_count(), 2);
    assert_eq!(
        individual
            .copies_for(reference.sequences()[0].id())
            .unwrap()
            .len(),
        2
    );
    let copies = individual
        .copies_for(reference.sequences()[0].id())
        .expect("copies for chr1");
    assert_eq!(copies[0].segments().segments().len(), 1);
    assert_eq!(copies[0].segments().segments()[0].source().len(), 6);
    assert!(copies[0].segments().is_shared_with(copies[1].segments()));
}
