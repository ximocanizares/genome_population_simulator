use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::reference::{ReversibleText, SequenceClass};
use genome_population_simulator::reference_fasta::{SequenceExpectation, validate_supplied_fasta};

static COUNTER: AtomicU64 = AtomicU64::new(0);

#[test]
fn supplied_fasta_preserves_n_case_crlf_and_non_utf8_header_bytes() {
    let directory = TestDirectory::new("valid");
    let path = directory.path().join("reference.fa");
    let bytes = b">chr1 description \xff\r\naCgTn\r\n>scaf-2\tboundary \r\nNNaa\r\n";
    fs::write(&path, bytes).expect("write fixture");
    let expected = [
        SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 5,
        },
        SequenceExpectation {
            id: "scaf-2",
            class: SequenceClass::Scaffold,
            length: 4,
        },
    ];

    let validated =
        validate_supplied_fasta(&path, "ref-main", &expected).expect("valid supplied FASTA");
    let sequences = validated.reference().sequences();
    assert_eq!(sequences[0].sequence().as_bytes(), b"ACGTN");
    assert_eq!(sequences[0].header().raw_bytes(), b"chr1 description \xff");
    assert!(matches!(
        sequences[0].header().description(),
        Some(ReversibleText::Bytes(bytes)) if bytes == b"description \xff"
    ));
    assert_eq!(sequences[1].header().raw_bytes(), b"scaf-2\tboundary ");
    assert_eq!(validated.source_identity().size_bytes(), bytes.len() as u64);
}

#[test]
fn supplied_fasta_rejects_empty_file_and_empty_sequence() {
    assert_error(b"", &[], "GPS-FMT-002");
    assert_error(
        b">chr1\n>chr2\nAC\n",
        &[
            expected("chr1", 1),
            SequenceExpectation {
                id: "chr2",
                class: SequenceClass::Scaffold,
                length: 2,
            },
        ],
        "GPS-FMT-002",
    );
}

#[test]
fn supplied_fasta_rejects_empty_and_duplicate_ids() {
    assert_error(b">   \nAC\n", &[expected("chr1", 2)], "GPS-FMT-003");
    assert_error(
        b">chr1\nAC\n>chr1 duplicate\nGT\n",
        &[expected("chr1", 2), expected("chr1", 2)],
        "GPS-FMT-004",
    );
}

#[test]
fn supplied_fasta_rejects_gzip_gap_iupac_and_other_symbol() {
    assert_error(&[0x1f, 0x8b, 0x08, 0x00], &[], "GPS-FMT-001");
    for symbol in *b"-R. " {
        let bytes = [b">chr1\nAC".as_slice(), &[symbol], b"T\n"].concat();
        let directory = TestDirectory::new("symbol");
        let path = directory.path().join("reference.fa");
        fs::write(&path, bytes).expect("write fixture");
        let error = validate_supplied_fasta(&path, "ref-main", &[expected("chr1", 4)])
            .expect_err("symbol must be rejected");
        assert_eq!(error.code(), "GPS-FMT-005");
        assert_eq!(
            error.context().get("position").map(String::as_str),
            Some("3")
        );
        assert_eq!(
            error.context().get("found").map(String::as_str),
            Some(format!("0x{symbol:02x}").as_str())
        );
    }
}

#[test]
fn supplied_fasta_rejects_sequence_catalog_or_length_mismatch() {
    assert_error(b">other\nAC\n", &[expected("chr1", 2)], "GPS-CONF-008");
    assert_error(b">chr1\nAC\n", &[expected("chr1", 3)], "GPS-CONF-008");
}

fn expected(id: &str, length: u64) -> SequenceExpectation<'_> {
    SequenceExpectation {
        id,
        class: SequenceClass::Chromosome,
        length,
    }
}

fn assert_error(bytes: &[u8], expected: &[SequenceExpectation<'_>], code: &str) {
    let directory = TestDirectory::new("invalid");
    let path = directory.path().join("reference.fa");
    fs::write(&path, bytes).expect("write fixture");
    let error = validate_supplied_fasta(&path, "ref-main", expected)
        .expect_err("invalid FASTA must be rejected");
    assert_eq!(error.code(), code);
    assert_eq!(error.exit_code(), 2);
}

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new(label: &str) -> Self {
        let counter = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "gps-fasta-validation-{label}-{}-{counter}",
            std::process::id()
        ));
        fs::create_dir(&path).expect("create isolated test directory");
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
