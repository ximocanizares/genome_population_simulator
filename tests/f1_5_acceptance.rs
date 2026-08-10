use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

const BINARY: &str = env!("CARGO_BIN_EXE_genome-population-simulator");
static COUNTER: AtomicU64 = AtomicU64::new(0);

#[test]
fn simulated_cli_runs_repeat_logically_and_pass_independent_oracles() {
    let root = temporary("simulated");
    fs::create_dir(&root).expect("root");
    let config = Path::new("tests/fixtures/config/general-seeded.toml");
    let first = root.join("first");
    let second = root.join("second");
    run_cli(config, &first, true);
    run_cli(config, &second, true);

    let first_fasta = first.join("resources/references/ref-main.fa");
    let second_fasta = second.join("resources/references/ref-main.fa");
    assert_eq!(
        fs::read(&first_fasta).expect("first FASTA"),
        fs::read(&second_fasta).expect("second FASTA")
    );
    assert_eq!(
        read_json(&first.join("effective-config.json")),
        read_json(&second.join("effective-config.json"))
    );
    assert_eq!(
        read_json(&first.join("manifest.json")),
        read_json(&second.join("manifest.json"))
    );
    assert!(
        fs::read(&first_fasta)
            .expect("FASTA")
            .split(|base| *base == b'\n')
            .skip(1)
            .flatten()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
    );

    validate_json(
        &first.join("effective-config.json"),
        "effective-config-v1.schema.json",
    );
    validate_json(
        &first.join("manifest.json"),
        "manifest-general-v1.schema.json",
    );
    assert!(!first_fasta.with_extension("fa.fai").exists());
    let acceptance_fasta = root.join("acceptance.fa");
    fs::copy(&first_fasta, &acceptance_fasta).expect("acceptance FASTA");
    let samtools = Command::new("samtools")
        .args(["faidx"])
        .arg(&acceptance_fasta)
        .status()
        .expect("samtools");
    assert!(samtools.success());
    assert!(acceptance_fasta.with_extension("fa.fai").is_file());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn supplied_cli_preserves_n_and_exact_input_while_invalid_input_has_no_destination() {
    let root = temporary("supplied");
    fs::create_dir(&root).expect("root");
    let source = root.join("reference.fa");
    let original = b">chrN exact description\r\naCnTN\r\n";
    fs::write(&source, original).expect("source");
    let config = root.join("general.toml");
    write_supplied_config(&config, "reference.fa", 5);
    let destination = root.join("valid-result");
    run_cli(&config, &destination, true);
    let copy = destination.join("resources/references/ref-supplied.fa");
    assert_eq!(fs::read(&copy).expect("copy"), original);
    assert_eq!(
        read_json(&destination.join("effective-config.json"))["reference"]["sequences"][0]["unknown_base_count"],
        2
    );
    assert!(
        Command::new("cmp")
            .args(["-s"])
            .arg(&source)
            .arg(&copy)
            .status()
            .expect("cmp")
            .success()
    );

    fs::write(&source, b">chrN\nACRTA\n").expect("invalid source");
    let rejected = root.join("rejected");
    let output = run_cli(&config, &rejected, false);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("GPS-FMT-005"));
    assert!(!rejected.exists());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn partial_ac_composition_is_completed_with_zero_g_t_before_cli_execution() {
    let root = temporary("partial-composition");
    fs::create_dir(&root).expect("root");
    let config = root.join("general.toml");
    fs::write(
        &config,
        r#"config_kind = "general"
schema_version = "1"
route = "direct"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[reference]
kind = "simulated"
id = "ref-partial"
default_ploidy = 2
composition = [
  { id = "A", category = "base", value = 0.75 },
  { id = "C", category = "base", value = 0.25 },
]
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 64
"#,
    )
    .expect("config");
    let destination = root.join("result");
    run_cli(&config, &destination, true);
    let effective = read_json(&destination.join("effective-config.json"));
    let values = effective["normalizations"][0]["values"]
        .as_array()
        .expect("composition");
    assert_eq!(values.len(), 4);
    assert_eq!(
        values[2],
        serde_json::json!({"id": "G", "category": "base", "before": 0.0, "after": 0.0})
    );
    assert_eq!(
        values[3],
        serde_json::json!({"id": "T", "category": "base", "before": 0.0, "after": 0.0})
    );
    let fasta = fs::read(destination.join("resources/references/ref-partial.fa")).expect("FASTA");
    assert!(fasta.iter().all(|byte| !matches!(byte, b'G' | b'T')));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn generated_composition_rejects_nonexact_ids_categories_and_duplicates_at_admission() {
    for (label, composition) in [
        (
            "lowercase",
            r#"[{ id = "a", category = "base", value = 1.0 }]"#,
        ),
        (
            "category",
            r#"[{ id = "A", category = "nucleotide", value = 1.0 }]"#,
        ),
        (
            "duplicate",
            r#"[{ id = "A", category = "base", value = 0.5 }, { id = "A", category = "base", value = 0.5 }]"#,
        ),
    ] {
        let root = temporary(label);
        fs::create_dir(&root).expect("root");
        let config = root.join("general.toml");
        fs::write(
            &config,
            format!(
                r#"config_kind = "general"
schema_version = "1"
route = "direct"
[reference]
kind = "simulated"
id = "ref-invalid"
default_ploidy = 2
composition = {composition}
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 8
"#
            ),
        )
        .expect("config");
        let destination = root.join("result");
        let output = run_cli(&config, &destination, false);
        assert_eq!(output.status.code(), Some(2), "{label}");
        assert!(String::from_utf8_lossy(&output.stderr).contains("GPS-CONF-008"));
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("cleanup");
    }
}

fn run_cli(config: &Path, destination: &Path, success: bool) -> std::process::Output {
    let output = Command::new(BINARY)
        .arg("run")
        .arg(config)
        .arg("--output")
        .arg(destination)
        .output()
        .expect("CLI");
    assert_eq!(
        output.status.success(),
        success,
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn read_json(path: &Path) -> serde_json::Value {
    serde_json::from_slice(&fs::read(path).expect("JSON document")).expect("valid JSON")
}

fn validate_json(instance: &Path, schema: &str) {
    let schema_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("schemas");
    let base_uri = format!("file://{}/", schema_root.display());
    let output = Command::new("python3")
        .args(["-m", "jsonschema", "--base-uri", &base_uri, "-i"])
        .arg(instance)
        .arg(schema)
        .current_dir(schema_root)
        .output()
        .expect("python-jsonschema");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn write_supplied_config(path: &Path, fasta: &str, length: u64) {
    fs::write(
        path,
        format!(
            r#"config_kind = "general"
schema_version = "1"
route = "direct"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[reference]
kind = "supplied"
id = "ref-supplied"
fasta = "{fasta}"

[[reference.sequences]]
id = "chrN"
class = "chromosome"
length = {length}
ploidy = 2
"#
        ),
    )
    .expect("config");
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f15-{label}-{}-{}",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ))
}
