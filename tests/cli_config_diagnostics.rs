use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

const BINARY: &str = env!("CARGO_BIN_EXE_genome-population-simulator");
static NONCE: AtomicU64 = AtomicU64::new(0);

fn run(command: &str, path: &Path, mode: &str, output: &Path) -> Output {
    Command::new(BINARY)
        .arg(command)
        .arg(path)
        .args(["--diagnostics", mode, "--output"])
        .arg(output)
        .output()
        .expect("CLI")
}

fn destination(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f12-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

#[test]
fn invalid_header_and_future_fields_have_equivalent_diagnostics_without_publication() {
    let cases = [
        ("invalid-missing-kind.toml", "GPS-CONF-002"),
        ("invalid-missing-version.toml", "GPS-CONF-004"),
        ("invalid-unknown-kind.toml", "GPS-CONF-003"),
        ("invalid-crossed-kind.toml", "GPS-CONF-003"),
        ("invalid-version-one.toml", "GPS-CONF-005"),
        ("invalid-malformed.toml", "GPS-CONF-001"),
        ("invalid-duplicate-key.toml", "GPS-CONF-007"),
        ("invalid-unknown-field.toml", "GPS-CONF-006"),
    ];
    for (fixture, code) in cases {
        let input = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/config")
            .join(fixture);
        let output = destination(fixture);
        let human = run("validate", &input, "human", &output);
        let json = run("validate", &input, "json", &output);
        assert_eq!(human.status.code(), Some(2), "{fixture}");
        assert_eq!(json.status.code(), Some(2), "{fixture}");
        assert!(String::from_utf8_lossy(&human.stderr).contains(code));
        assert!(String::from_utf8_lossy(&json.stdout).contains(code));
        assert!(!output.exists());
    }
}

#[test]
fn validation_is_read_only() {
    for fixture in ["general-minimal.toml", "general-seeded.toml"] {
        let input = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/config")
            .join(fixture);
        let output = destination(fixture);
        let result = run("validate", &input, "human", &output);
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(!output.exists());
    }
    let supplied_root = destination("supplied-read-only");
    fs::create_dir(&supplied_root).expect("supplied fixture root");
    let input = supplied_root.join("general.toml");
    fs::write(
        supplied_root.join("reference.fa"),
        format!(">chr1\n{}\n", "A".repeat(1000)),
    )
    .expect("supplied FASTA");
    fs::write(
        &input,
        include_str!("fixtures/config/general-supplied.toml"),
    )
    .expect("supplied config");
    let output = supplied_root.join("result");
    let result = run("validate", &input, "human", &output);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!output.exists());
    fs::remove_dir_all(supplied_root).expect("supplied fixture cleanup");
}

#[test]
fn root_scenario_field_is_rejected_by_run_before_opening_links() {
    let input = destination("future.toml");
    let missing_scenario = input.with_extension("missing-scenario.toml");
    assert!(!missing_scenario.exists());
    let general = include_str!("fixtures/config/general-minimal.toml");
    let root_scenario = format!(
        "scenarios = [{{ id = \"s\", path = \"{}\", route = \"evolution\" }}]\n[reference]",
        missing_scenario.display()
    );
    let source = general.replacen("[reference]", &root_scenario, 1);
    fs::write(&input, source).expect("fixture");
    let output = destination("future-output");
    let result = run("run", &input, "human", &output);
    assert_eq!(result.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("GPS-CONF-006"));
    assert!(!stderr.contains("GPS-INPUT-001"));
    assert!(!missing_scenario.exists());
    assert!(!output.exists());
    fs::remove_file(input).expect("cleanup");
}

#[test]
fn simulated_reduction_is_rejected_by_run_without_publication() {
    let input = destination("simulated-reduction.toml");
    fs::write(
        &input,
        r#"config_kind = "general"
schema_version = "2"
seed = "0000000000000000000000000000000000000000000000000000000000000000"
[reference]
kind = "simulated"
id = "ref-simulated"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 100
[[reduced_references]]
id = "ref-reduced"
source_reference_id = "ref-simulated"
region_size_bases = 10
target_bases = 50
"#,
    )
    .expect("fixture");
    let output = destination("simulated-reduction-output");

    let result = run("run", &input, "human", &output);

    assert_eq!(result.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&result.stderr).contains("GPS-CONF-008"));
    assert!(!output.exists());
    fs::remove_file(input).expect("cleanup");
}
