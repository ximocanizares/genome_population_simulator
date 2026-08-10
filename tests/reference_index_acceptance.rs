use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::reference::ReferenceIndexKind;
use genome_population_simulator::reference_fasta::{
    create_fai_for_acceptance, query_fai_for_acceptance,
};

static COUNTER: AtomicU64 = AtomicU64::new(0);

#[test]
fn normal_f1_run_has_zero_indexes_and_acceptance_query_creates_fai_on_demand() {
    let directory = TestDirectory::new();
    let effective =
        load_and_validate_config(Path::new("tests/fixtures/config/general-seeded.toml"))
            .expect("effective configuration");
    let result = directory.path().join("result");
    execute_validated_config(&effective, &result).expect("normal F1 run");
    assert_eq!(count_fai(&result), 0, "normal F1 must publish zero indexes");

    let internal_fasta = directory.path().join("internal.fa");
    fs::write(&internal_fasta, b">chr1\nACGTACGT\n").expect("known FASTA");
    let internal_fai = directory.path().join("internal.fa.fai");
    let index = create_fai_for_acceptance(&internal_fasta, &internal_fai).expect("on-demand FAI");
    assert_eq!(index.kind(), ReferenceIndexKind::Fai);
    assert!(internal_fai.is_file());

    let internal = query_fai_for_acceptance(&internal_fasta, &internal_fai, "chr1:2-5")
        .expect("indexed acceptance query");
    assert_eq!(internal, b"CGTA", "literal biological oracle");

    let external_fasta = directory.path().join("external.fa");
    fs::copy(&internal_fasta, &external_fasta).expect("independent FASTA copy");
    let external_fai = directory.path().join("external.fa.fai");
    assert!(!external_fai.exists());
    let external = Command::new("samtools")
        .args(["faidx"])
        .arg(&external_fasta)
        .arg("chr1:2-5")
        .output()
        .expect("approved samtools validator");
    assert!(
        external.status.success(),
        "{}",
        String::from_utf8_lossy(&external.stderr)
    );
    let external_bases = external
        .stdout
        .split(|byte| *byte == b'\n')
        .skip(1)
        .flatten()
        .copied()
        .collect::<Vec<_>>();
    assert_eq!(internal, external_bases);
    assert!(external_fai.is_file(), "samtools must build its own FAI");
    assert_ne!(internal_fai, external_fai);
}

fn count_fai(root: &Path) -> usize {
    fs::read_dir(root)
        .expect("read result tree")
        .filter_map(Result::ok)
        .map(|entry| {
            if entry.file_type().expect("file type").is_dir() {
                count_fai(&entry.path())
            } else {
                usize::from(
                    entry
                        .path()
                        .extension()
                        .is_some_and(|extension| extension == "fai"),
                )
            }
        })
        .sum()
}

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new() -> Self {
        let counter = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "gps-reference-index-{}-{counter}",
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
