use std::path::PathBuf;

use genome_population_simulator::configuration::SequenceClass;
use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    AlphabetValidation, DnaSequence, FastaHeader, Reference, ReferenceError, ReferenceIndex,
    ReferenceIndexKind, ReferenceOrigin, ReferenceSequence, ReversibleText,
};

#[test]
fn dna_accepts_both_cases_and_normalizes_only_the_internal_value() {
    let sequence = DnaSequence::try_from_bytes(b"aCgTnN").expect("valid DNA");

    assert_eq!(
        sequence.as_bytes().expect("bases materializadas"),
        b"ACGTNN"
    );
    assert_eq!(sequence.unknown_base_count(), 2);
    assert_eq!(
        sequence.validation(),
        AlphabetValidation::AcgtnCaseInsensitive
    );
}

#[test]
fn dna_reports_the_exact_one_based_position_and_raw_symbol() {
    assert_eq!(
        DnaSequence::try_from_bytes(b"ACGRN"),
        Err(ReferenceError::InvalidDnaSymbol {
            position: 4,
            symbol: b'R',
        })
    );
}

#[test]
fn fasta_header_keeps_raw_bytes_and_a_utf8_description_separately() {
    let header = FastaHeader::parse(b"chr1  chromosome one").expect("valid header");

    assert_eq!(header.id_text(), "chr1");
    assert_eq!(header.raw_bytes(), b"chr1  chromosome one");
    assert_eq!(
        header.description(),
        Some(&ReversibleText::Utf8(" chromosome one".to_owned()))
    );
}

#[test]
fn fasta_header_preserves_non_utf8_description_reversibly() {
    let header = FastaHeader::parse(b"scaf-1 description \xff").expect("valid header");

    assert_eq!(header.raw_bytes(), b"scaf-1 description \xff");
    assert_eq!(
        header.description(),
        Some(&ReversibleText::Bytes(vec![
            b'd', b'e', b's', b'c', b'r', b'i', b'p', b't', b'i', b'o', b'n', b' ', 0xff,
        ]))
    );
    assert_eq!(header.description_encoding(), Some("bytes-hex"));
    assert_eq!(
        header.encoded_description(),
        Some("6465736372697074696f6e20ff".to_owned())
    );
}

#[test]
fn fasta_header_preserves_description_boundary_whitespace_reversibly() {
    let header = FastaHeader::parse(b"chr1  padded description \t").expect("valid header");

    assert_eq!(header.raw_bytes(), b"chr1  padded description \t");
    assert_eq!(
        header.description(),
        Some(&ReversibleText::Utf8(" padded description \t".to_owned()))
    );
}

#[test]
fn reference_sequence_retains_the_validated_normalized_bases() {
    let sequence = sequence("chr1", SequenceClass::Chromosome, b"aCgTn");

    assert_eq!(
        sequence
            .sequence()
            .as_bytes()
            .expect("bases materializadas"),
        b"ACGTN"
    );
    assert_eq!(sequence.length(), 5);
    assert_eq!(sequence.unknown_base_count(), 1);
}

#[test]
fn reference_preserves_sequence_order_origin_and_zero_or_one_f1_index() {
    let reference_id = ReferenceId::parse("ref-main", "execution").expect("reference ID");
    let first = sequence("chr1", SequenceClass::Chromosome, b"ACGTN");
    let second = sequence("scaffold.2", SequenceClass::Scaffold, b"NNAA");
    let mut reference = Reference::new(
        reference_id,
        ReferenceOrigin::Supplied {
            original_path: PathBuf::from("inputs/reference.fa"),
        },
        vec![first, second],
    )
    .expect("valid reference");

    assert_eq!(
        reference
            .sequences()
            .iter()
            .map(|sequence| sequence.id().as_str())
            .collect::<Vec<_>>(),
        ["chr1", "scaffold.2"]
    );
    assert_eq!(reference.total_length(), 9);
    assert!(reference.indexes().is_empty());

    reference
        .add_index(ReferenceIndex::new(
            ReferenceIndexKind::Fai,
            PathBuf::from("resources/ref-main.fa.fai"),
        ))
        .expect("first FAI");
    assert_eq!(reference.indexes().len(), 1);
    assert_eq!(
        reference.add_index(ReferenceIndex::new(
            ReferenceIndexKind::Fai,
            PathBuf::from("resources/duplicate.fai"),
        )),
        Err(ReferenceError::DuplicateIndexKind(ReferenceIndexKind::Fai))
    );
}

#[test]
fn reference_rejects_empty_and_duplicate_sequence_catalogs() {
    let reference_id = ReferenceId::parse("ref-main", "execution").expect("reference ID");
    assert!(matches!(
        Reference::new(reference_id, ReferenceOrigin::Generated, Vec::new()),
        Err(ReferenceError::EmptyReference)
    ));

    let reference_id = ReferenceId::parse("ref-main", "execution").expect("reference ID");
    assert!(matches!(
        Reference::new(
            reference_id,
            ReferenceOrigin::Generated,
            vec![
                sequence("chr1", SequenceClass::Chromosome, b"A"),
                sequence("chr1", SequenceClass::Scaffold, b"C"),
            ],
        ),
        Err(ReferenceError::DuplicateSequenceId(id)) if id == "chr1"
    ));
}

fn sequence(id: &str, class: SequenceClass, bases: &[u8]) -> ReferenceSequence {
    let id = SequenceId::parse(id, "ref-main").expect("sequence ID");
    let header = FastaHeader::parse(id.as_str().as_bytes()).expect("header");
    let dna = DnaSequence::try_from_bytes(bases).expect("DNA");
    ReferenceSequence::new(id, header, class, &dna).expect("sequence metadata")
}
