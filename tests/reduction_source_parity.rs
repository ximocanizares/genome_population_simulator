//! Las dos fuentes de bases de una reducción deben ser indistinguibles.
//!
//! `MaterializedSource` sirve los tramos desde la referencia en memoria e
//! `IndexedSource` desde el FASTA en disco mediante su FAI. Si divergieran, una
//! reducción reconstruida por una no coincidiría con la validada por la otra, y
//! la garantía de que la derivada se corresponde con su origen dejaría de valer.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_fasta::{SequenceExpectation, validate_supplied_fasta};
use genome_population_simulator::reference_reduction_io::{
    IndexedSource, MaterializedSource, ReductionSource,
};

static COUNTER: AtomicU64 = AtomicU64::new(0);

struct TestDirectory {
    path: PathBuf,
}

impl TestDirectory {
    fn new(label: &str) -> Self {
        let unique = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!("gps-reduction-source-{label}-{unique}"));
        fs::create_dir_all(&path).expect("create test directory");
        Self { path }
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

fn fixture() -> (TestDirectory, PathBuf) {
    let directory = TestDirectory::new("parity");
    let path = directory.path.join("reference.fa");
    // Dos secuencias, líneas de distinto ancho y minúsculas que la validación
    // normaliza: el índice debe resolver las mismas bases que la memoria.
    let bytes = b">chr1 first\nACGTACGTAC\nGTNNacgtac\n>chr2 second\nTTTTGGGGCC\nAAAA\n";
    fs::write(&path, bytes).expect("write fixture");
    (directory, path)
}

fn expectations() -> Vec<SequenceExpectation<'static>> {
    vec![
        SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 20,
        },
        SequenceExpectation {
            id: "chr2",
            class: SequenceClass::Scaffold,
            length: 14,
        },
    ]
}

#[test]
fn indexed_and_materialized_sources_return_the_same_bases() {
    let (_directory, path) = fixture();
    let expected = expectations();
    let validated = validate_supplied_fasta(&path, "ref-main", &expected).expect("validate");
    let materialized = MaterializedSource::new(validated.reference());
    let catalog = materialized.catalog().to_vec();
    let indexed = IndexedSource::new(
        ReferenceId::parse("ref-main", "test").expect("reference id"),
        catalog,
        &path,
    )
    .expect("build index");

    assert_eq!(indexed.reference_id(), materialized.reference_id());
    assert_eq!(indexed.catalog(), materialized.catalog());

    for (sequence, length) in [("chr1", 20_usize), ("chr2", 14)] {
        let id = SequenceId::parse(sequence, "ref-main").expect("sequence id");
        for start in 0..length {
            for end in (start + 1)..=length {
                let from_memory = materialized
                    .interval(&id, start, end)
                    .expect("materialized interval");
                let from_index = indexed.interval(&id, start, end).expect("indexed interval");
                assert_eq!(
                    from_index, from_memory,
                    "{sequence}:{start}-{end} difiere entre memoria e índice"
                );
            }
        }
    }
}

#[test]
fn both_sources_reject_an_interval_outside_the_sequence() {
    let (_directory, path) = fixture();
    let expected = expectations();
    let validated = validate_supplied_fasta(&path, "ref-main", &expected).expect("validate");
    let materialized = MaterializedSource::new(validated.reference());
    let indexed = IndexedSource::new(
        ReferenceId::parse("ref-main", "test").expect("reference id"),
        materialized.catalog().to_vec(),
        &path,
    )
    .expect("build index");
    let id = SequenceId::parse("chr1", "ref-main").expect("sequence id");

    assert!(materialized.interval(&id, 0, 21).is_err());
    assert!(indexed.interval(&id, 0, 21).is_err());
    // Un intervalo vacío no es un error en la fuente en memoria: devuelve un
    // fragmento vacío. La indexada conserva ese comportamiento para que ambas
    // sean indistinguibles, aunque el índice no sepa representar esa región.
    assert_eq!(
        materialized.interval(&id, 5, 5).expect("empty interval"),
        Vec::<u8>::new()
    );
    assert_eq!(
        indexed.interval(&id, 5, 5).expect("empty interval"),
        Vec::<u8>::new()
    );
}

#[test]
fn both_sources_reject_an_unknown_sequence() {
    let (_directory, path) = fixture();
    let expected = expectations();
    let validated = validate_supplied_fasta(&path, "ref-main", &expected).expect("validate");
    let materialized = MaterializedSource::new(validated.reference());
    let indexed = IndexedSource::new(
        ReferenceId::parse("ref-main", "test").expect("reference id"),
        materialized.catalog().to_vec(),
        &path,
    )
    .expect("build index");
    let id = SequenceId::parse("ausente", "ref-main").expect("sequence id");

    assert!(materialized.interval(&id, 0, 4).is_err());
    assert!(indexed.interval(&id, 0, 4).is_err());
}
