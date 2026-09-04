use std::fs;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use genome_population_simulator::f5_admission::{
    AdmissionPhase, BoundedInput, F5_RESOURCE_ACCOUNTING_ALGORITHM, F5ResourceLedger, InputRole,
    LedgerCharge, LedgerLimits, OperationCost, OutputReservation, WorkingMemoryPeak,
};
use sha2::{Digest, Sha256};

fn temporary(label: &str) -> PathBuf {
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock after epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("gps-f5-bounded-input-{label}-{unique}"))
}

fn independent_file_identity(path: &Path) -> (u64, String) {
    let mut file = fs::File::open(path).expect("literal file for independent oracle");
    let mut digest = Sha256::new();
    let mut size = 0_u64;
    let mut byte = [0_u8; 1];
    loop {
        let read = file.read(&mut byte).expect("oracle read");
        if read == 0 {
            break;
        }
        size = size.checked_add(read as u64).expect("oracle size fits u64");
        digest.update(&byte[..read]);
    }
    (size, format!("{:x}", digest.finalize()))
}

struct FragmentedReader {
    bytes: Vec<u8>,
    offset: usize,
    observed: Arc<AtomicUsize>,
}

impl Read for FragmentedReader {
    fn read(&mut self, destination: &mut [u8]) -> io::Result<usize> {
        if self.offset == self.bytes.len() {
            return Ok(0);
        }
        destination[0] = self.bytes[self.offset];
        self.offset += 1;
        self.observed.fetch_add(1, Ordering::Relaxed);
        Ok(1)
    }
}

struct FailingReader {
    bytes_before_error: usize,
    offset: usize,
}

struct SpyReader {
    reads: Arc<AtomicUsize>,
}

impl Read for SpyReader {
    fn read(&mut self, _destination: &mut [u8]) -> io::Result<usize> {
        self.reads.fetch_add(1, Ordering::Relaxed);
        Ok(0)
    }
}

impl Read for FailingReader {
    fn read(&mut self, destination: &mut [u8]) -> io::Result<usize> {
        if self.offset == self.bytes_before_error {
            return Err(io::Error::new(
                io::ErrorKind::Interrupted,
                "directed read failure",
            ));
        }
        destination[0] = b'x';
        self.offset += 1;
        Ok(1)
    }
}

fn bounded<R: Read>(reader: R, limit: u64) -> BoundedInput<R> {
    BoundedInput::new(
        InputRole::ProfileMetadata,
        AdmissionPhase::ProfileMetadata,
        PathBuf::from("/authoring/profiles/profile.metadata.json"),
        limit,
        reader,
    )
}

fn bootstrap_bounded<R: Read>(reader: R, limit: u64) -> BoundedInput<R> {
    BoundedInput::new(
        InputRole::GeneralToml,
        AdmissionPhase::Bootstrap,
        PathBuf::from("/authoring/general.toml"),
        limit,
        reader,
    )
}

fn profile_ledger(working_memory_limit: u64) -> F5ResourceLedger {
    let mut ledger = F5ResourceLedger::new(LedgerLimits::new(
        OperationCost::new(u64::MAX),
        WorkingMemoryPeak::new(working_memory_limit),
        OutputReservation::new(u64::MAX),
    ));
    ledger
        .transition_to(AdmissionPhase::ProfileMetadata)
        .expect("profile metadata follows bootstrap");
    ledger
}

#[test]
fn fragmented_reader_reaches_exact_eof_and_hashes_the_same_literal_bytes() {
    let bytes = b"one-byte fragments".to_vec();
    let observed = Arc::new(AtomicUsize::new(0));
    let mut ledger = profile_ledger(bytes.len() as u64);
    let read = bounded(
        FragmentedReader {
            bytes: bytes.clone(),
            offset: 0,
            observed: Arc::clone(&observed),
        },
        bytes.len() as u64,
    )
    .read_to_end_with_ledger(Some(&mut ledger))
    .expect("exact EOF is inclusive");

    let mut oracle = Sha256::new();
    oracle.update(&bytes);
    assert_eq!(read.bytes(), bytes);
    assert_eq!(observed.load(Ordering::Relaxed), bytes.len());
    assert_eq!(read.identity().role(), InputRole::ProfileMetadata);
    assert_eq!(
        read.identity().normalized_path(),
        Path::new("/authoring/profiles/profile.metadata.json")
    );
    assert_eq!(read.identity().observed_size_bytes(), bytes.len() as u64);
    assert_eq!(
        read.identity().observed_sha256(),
        format!("{:x}", oracle.finalize())
    );
}

#[test]
fn excess_byte_is_observed_only_to_reject_and_never_enters_the_identity() {
    let observed = Arc::new(AtomicUsize::new(0));
    let mut ledger = profile_ledger(5);
    let error = bounded(
        FragmentedReader {
            bytes: b"abcd!".to_vec(),
            offset: 0,
            observed: Arc::clone(&observed),
        },
        4,
    )
    .read_to_end_with_ledger(Some(&mut ledger))
    .expect_err("the fifth byte exceeds the inclusive limit");

    assert_eq!(observed.load(Ordering::Relaxed), 5);
    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(
        error.context().get("algorithm"),
        Some(&F5_RESOURCE_ACCOUNTING_ALGORITHM.to_owned())
    );
    assert_eq!(
        error.context().get("phase"),
        Some(&"profile_metadata".to_owned())
    );
    assert_eq!(
        error.context().get("role"),
        Some(&"profile_metadata".to_owned())
    );
    assert_eq!(
        error.context().get("dimension"),
        Some(&"physical_input_bytes".to_owned())
    );
    assert_eq!(error.context().get("accumulated"), Some(&"5".to_owned()));
    assert_eq!(error.context().get("limit"), Some(&"4".to_owned()));
    assert_eq!(ledger.phase(), AdmissionPhase::Rejected);
}

#[test]
fn intermediate_reader_error_is_preserved() {
    let mut ledger = profile_ledger(4);
    let error = bounded(
        FailingReader {
            bytes_before_error: 2,
            offset: 0,
        },
        4,
    )
    .read_to_end_with_ledger(Some(&mut ledger))
    .expect_err("reader failure must not become EOF or an over-limit success");

    assert_eq!(error.code(), "GPS-INPUT-001");
    assert!(error.render_human().contains("directed read failure"));
}

#[test]
fn non_bootstrap_read_without_ledger_does_not_touch_the_reader() {
    let reads = Arc::new(AtomicUsize::new(0));
    let error = bounded(
        SpyReader {
            reads: Arc::clone(&reads),
        },
        1,
    )
    .read_to_end()
    .expect_err("profile input requires a ledger before the first read");

    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(reads.load(Ordering::Relaxed), 0);
}

#[test]
fn ledger_role_or_phase_mismatch_does_not_touch_the_reader() {
    let reads = Arc::new(AtomicUsize::new(0));
    let mut ledger = F5ResourceLedger::new(LedgerLimits::new(
        OperationCost::new(10),
        WorkingMemoryPeak::new(10),
        OutputReservation::new(10),
    ));
    let error = bounded(
        SpyReader {
            reads: Arc::clone(&reads),
        },
        1,
    )
    .read_to_end_with_ledger(Some(&mut ledger))
    .expect_err("profile metadata cannot be read during bootstrap");

    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(reads.load(Ordering::Relaxed), 0);
    assert_eq!(ledger.phase(), AdmissionPhase::Bootstrap);

    let role_reads = Arc::new(AtomicUsize::new(0));
    let mut profile_ledger = profile_ledger(10);
    let error = BoundedInput::new(
        InputRole::ScenarioToml,
        AdmissionPhase::ProfileMetadata,
        PathBuf::from("/authoring/scenario.toml"),
        1,
        SpyReader {
            reads: Arc::clone(&role_reads),
        },
    )
    .read_to_end_with_ledger(Some(&mut profile_ledger))
    .expect_err("a bootstrap role cannot be read during profile metadata");

    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(role_reads.load(Ordering::Relaxed), 0);
    assert_eq!(profile_ledger.phase(), AdmissionPhase::ProfileMetadata);
}

#[test]
fn directed_overflow_is_rejected_without_changing_accumulated_ledger() {
    let mut ledger = F5ResourceLedger::new(LedgerLimits::new(
        OperationCost::new(u64::MAX),
        WorkingMemoryPeak::new(10),
        OutputReservation::new(10),
    ));
    ledger
        .charge(
            InputRole::GeneralToml,
            LedgerCharge::Operation(OperationCost::new(u64::MAX)),
        )
        .expect("maximum exact operation cost is inclusive");
    let before = ledger.totals();
    let error = ledger
        .charge(
            InputRole::GeneralToml,
            LedgerCharge::Operation(OperationCost::new(1)),
        )
        .expect_err("checked addition must reject overflow");

    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(
        error.context().get("dimension"),
        Some(&"operation_units".to_owned())
    );
    assert_eq!(ledger.totals(), before);
    assert_eq!(ledger.phase(), AdmissionPhase::Rejected);
    assert!(
        ledger
            .charge(
                InputRole::GeneralToml,
                LedgerCharge::Operation(OperationCost::new(0)),
            )
            .is_err()
    );
    assert!(
        ledger
            .transition_to(AdmissionPhase::ProfileMetadata)
            .is_err()
    );
    let reads = Arc::new(AtomicUsize::new(0));
    assert!(
        bootstrap_bounded(
            SpyReader {
                reads: Arc::clone(&reads),
            },
            1,
        )
        .read_to_end_with_ledger(Some(&mut ledger))
        .is_err()
    );
    assert_eq!(reads.load(Ordering::Relaxed), 0);
    assert_eq!(ledger.totals(), before);
}

#[test]
fn metadata_size_is_not_an_admission_authority_and_file_oracle_recomputes_identity() {
    let root = temporary("metadata-does-not-authorize");
    fs::create_dir_all(&root).expect("fixture directory");
    let path = root.join("profile.metadata.json");
    fs::write(&path, b"old").expect("first literal");
    let stale_metadata_size = fs::metadata(&path).expect("metadata only for test").len();
    fs::write(&path, b"newer literal bytes").expect("mutated literal");

    let (oracle_size, oracle_sha) = independent_file_identity(&path);
    let mut ledger = profile_ledger(oracle_size);
    let read = BoundedInput::new(
        InputRole::ProfileMetadata,
        AdmissionPhase::ProfileMetadata,
        path.clone(),
        oracle_size,
        fs::File::open(&path).expect("same handle used for BoundedInput"),
    )
    .read_to_end_with_ledger(Some(&mut ledger))
    .expect("admission follows observed bytes, not stale metadata");

    assert_ne!(stale_metadata_size, oracle_size);
    assert_eq!(read.identity().observed_size_bytes(), oracle_size);
    assert_eq!(read.identity().observed_sha256(), oracle_sha);
    fs::remove_dir_all(root).expect("fixture cleanup");
}

#[test]
fn ledger_is_cumulative_and_rejection_does_not_reserve_above_the_limit() {
    let mut ledger = F5ResourceLedger::new(LedgerLimits::new(
        OperationCost::new(5),
        WorkingMemoryPeak::new(3),
        OutputReservation::new(5),
    ));
    ledger
        .charge(
            InputRole::GeneralToml,
            LedgerCharge::Operation(OperationCost::new(2)),
        )
        .expect("first cumulative charge");
    ledger
        .charge(
            InputRole::ScenarioToml,
            LedgerCharge::Operation(OperationCost::new(3)),
        )
        .expect("inclusive total limit");
    let before = ledger.totals();
    let error = ledger
        .charge(
            InputRole::ProfileMetadata,
            LedgerCharge::Output(OutputReservation::new(6)),
        )
        .expect_err("prospective rejection");

    assert_eq!(
        error.context().get("dimension"),
        Some(&"output_bytes".to_owned())
    );
    assert_eq!(error.context().get("accumulated"), Some(&"6".to_owned()));
    assert_eq!(ledger.totals(), before);
    assert_eq!(ledger.phase(), AdmissionPhase::Rejected);
    assert!(
        ledger
            .charge(
                InputRole::ProfileMetadata,
                LedgerCharge::Output(OutputReservation::new(0)),
            )
            .is_err()
    );
    assert!(
        ledger
            .transition_to(AdmissionPhase::ProfileMetadata)
            .is_err()
    );
    let reads = Arc::new(AtomicUsize::new(0));
    assert!(
        bootstrap_bounded(
            SpyReader {
                reads: Arc::clone(&reads),
            },
            1,
        )
        .read_to_end_with_ledger(Some(&mut ledger))
        .is_err()
    );
    assert_eq!(reads.load(Ordering::Relaxed), 0);

    let mut memory_ledger = F5ResourceLedger::new(LedgerLimits::new(
        OperationCost::new(5),
        WorkingMemoryPeak::new(3),
        OutputReservation::new(5),
    ));
    let before_memory = memory_ledger.totals();
    let error = bootstrap_bounded(io::Cursor::new(b"four".to_vec()), 4)
        .read_to_end_with_ledger(Some(&mut memory_ledger))
        .expect_err("the prospective buffer peak exceeds three bytes");
    assert_eq!(
        error.context().get("dimension"),
        Some(&"working_memory_bytes".to_owned())
    );
    assert_eq!(
        memory_ledger.totals(),
        before_memory,
        "failed growth did not reserve memory"
    );
}
