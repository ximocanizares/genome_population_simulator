use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_fasta::{
    CANONICAL_FASTA_VERSION, CANONICAL_LINE_BASE_COUNT, write_generated_fasta,
};
use genome_population_simulator::reference_generation::{
    BaseFrequency, GeneratedSequenceSpec, generate_reference,
};

const ROOT: RootSeed = RootSeed::from_bytes([
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31,
]);

#[test]
fn generated_reference_respects_order_lengths_classes_and_degenerate_composition() {
    let specs = [
        GeneratedSequenceSpec {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 5,
        },
        GeneratedSequenceSpec {
            id: "scaf-2",
            class: SequenceClass::Scaffold,
            length: 3,
        },
    ];
    let reference = generate_reference(
        "reference-main",
        ROOT,
        &specs,
        &composition(1.0, 0.0, 0.0, 0.0),
    )
    .expect("generated reference");

    assert_eq!(reference.total_length(), 8);
    assert_eq!(reference.sequences()[0].sequence().as_bytes(), b"AAAAA");
    assert_eq!(reference.sequences()[1].sequence().as_bytes(), b"AAA");
    assert_eq!(reference.sequences()[1].class(), SequenceClass::Scaffold);
    assert_eq!(reference.sequences()[0].unknown_base_count(), 0);
}

#[test]
fn generated_reference_is_repeated_by_seed_and_changes_with_another_seed() {
    let specs = [GeneratedSequenceSpec {
        id: "chr1",
        class: SequenceClass::Chromosome,
        length: 32,
    }];
    let frequencies = composition(0.25, 0.25, 0.25, 0.25);
    let first = generate_reference("reference-main", ROOT, &specs, &frequencies).unwrap();
    let repeated = generate_reference("reference-main", ROOT, &specs, &frequencies).unwrap();
    let changed = generate_reference(
        "reference-main",
        RootSeed::from_bytes([9; 32]),
        &specs,
        &frequencies,
    )
    .unwrap();

    let first = first.sequences()[0].sequence().as_bytes();
    assert_eq!(first, repeated.sequences()[0].sequence().as_bytes());
    assert_ne!(first, changed.sequences()[0].sequence().as_bytes());
    assert_eq!(first, b"CCCATGGGGAGGTCTAATCCAAGCGACTTACA");
    assert!(
        first
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
    );
}

#[test]
fn generated_reference_rejects_a_non_acgt_catalog() {
    let specs = [GeneratedSequenceSpec {
        id: "chr1",
        class: SequenceClass::Chromosome,
        length: 1,
    }];
    let error = generate_reference(
        "reference-main",
        ROOT,
        &specs,
        &[
            BaseFrequency {
                base: b'A',
                probability: 0.25,
            },
            BaseFrequency {
                base: b'C',
                probability: 0.25,
            },
            BaseFrequency {
                base: b'G',
                probability: 0.25,
            },
            BaseFrequency {
                base: b'N',
                probability: 0.25,
            },
        ],
    )
    .expect_err("generated N is not an admitted composition outcome");
    assert_eq!(error.code(), "GPS-CONF-008");
}

#[test]
fn generated_fasta_v1_has_literal_headers_lf_order_and_eighty_base_lines() {
    assert_eq!(CANONICAL_FASTA_VERSION, "1");
    assert_eq!(CANONICAL_LINE_BASE_COUNT, 80);
    let specs = [
        GeneratedSequenceSpec {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 81,
        },
        GeneratedSequenceSpec {
            id: "scaf-2",
            class: SequenceClass::Scaffold,
            length: 2,
        },
    ];
    let reference = generate_reference(
        "reference-main",
        ROOT,
        &specs,
        &composition(1.0, 0.0, 0.0, 0.0),
    )
    .unwrap();
    let mut bytes = Vec::new();
    write_generated_fasta(&reference, &mut bytes).expect("canonical FASTA");
    let expected = [b">chr1\n".as_slice(), &[b'A'; 80], b"\nA\n>scaf-2\nAA\n"].concat();
    assert_eq!(bytes, expected);
    assert!(!bytes.contains(&b'\r'));
}

fn composition(a: f64, c: f64, g: f64, t: f64) -> [BaseFrequency; 4] {
    [
        BaseFrequency {
            base: b'A',
            probability: a,
        },
        BaseFrequency {
            base: b'C',
            probability: c,
        },
        BaseFrequency {
            base: b'G',
            probability: g,
        },
        BaseFrequency {
            base: b'T',
            probability: t,
        },
    ]
}
