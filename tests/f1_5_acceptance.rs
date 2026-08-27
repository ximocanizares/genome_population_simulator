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
        "effective-config-v3.schema.json",
    );
    validate_json(
        &first.join("manifest.json"),
        "manifest-general-v3.schema.json",
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
fn declared_and_system_csprng_roots_are_effective_before_reference_generation() {
    let root = temporary("declared-and-csprng");
    fs::create_dir(&root).expect("root");

    let declared_destination = root.join("declared");
    run_cli(
        Path::new("tests/fixtures/config/general-seeded.toml"),
        &declared_destination,
        true,
    );
    let declared = read_json(&declared_destination.join("effective-config.json"));
    assert_eq!(declared["seed"]["origin"], "declared");
    assert_eq!(
        declared["seed"]["root"],
        "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
    );

    let seedless_config = Path::new("tests/fixtures/config/general-minimal.toml");
    let generated_destination = root.join("generated");
    run_cli(seedless_config, &generated_destination, true);
    let generated = read_json(&generated_destination.join("effective-config.json"));
    assert_eq!(generated["seed"]["origin"], "system_csprng");
    let generated_root = generated["seed"]["root"]
        .as_str()
        .expect("generated root seed");
    assert_eq!(generated_root.len(), 64);
    assert!(
        generated_root
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    );

    let replay_config = root.join("replay.toml");
    let seedless_source = fs::read_to_string(seedless_config).expect("seedless config");
    fs::write(
        &replay_config,
        format!("seed = \"{generated_root}\"\n{seedless_source}"),
    )
    .expect("replay config");
    let replay_destination = root.join("replay");
    run_cli(&replay_config, &replay_destination, true);
    let replay = read_json(&replay_destination.join("effective-config.json"));
    assert_eq!(replay["seed"]["origin"], "declared");
    assert_eq!(replay["seed"]["root"], generated_root);
    assert_eq!(
        fs::read(generated_destination.join("resources/references/ref-main.fa"))
            .expect("generated FASTA"),
        fs::read(replay_destination.join("resources/references/ref-main.fa"))
            .expect("replayed FASTA")
    );

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
    let acceptance_fasta = root.join("supplied-acceptance.fa");
    fs::copy(&copy, &acceptance_fasta).expect("independent supplied FASTA copy");
    let samtools = Command::new("samtools")
        .args(["faidx"])
        .arg(&acceptance_fasta)
        .status()
        .expect("samtools");
    assert!(samtools.success());
    assert!(acceptance_fasta.with_extension("fa.fai").is_file());
    validate_json(
        &destination.join("effective-config.json"),
        "effective-config-v3.schema.json",
    );
    validate_json(
        &destination.join("manifest.json"),
        "manifest-general-v3.schema.json",
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
fn partial_ac_composition_is_rejected_before_cli_execution() {
    let root = temporary("partial-composition");
    fs::create_dir(&root).expect("root");
    let config = root.join("general.toml");
    fs::write(
        &config,
        r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[reference]
kind = "simulated"
id = "ref-partial"
[reference.composition]
A = 0.75
C = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 64
"#,
    )
    .expect("config");
    let destination = root.join("result");
    let output = run_cli(&config, &destination, false);
    assert_eq!(output.status.code(), Some(2));
    assert!(!destination.exists());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn generated_composition_rejects_unknown_and_duplicate_keys_at_admission() {
    for (label, composition, expected_code) in [
        (
            "lowercase",
            "a = 0.25\nC = 0.25\nG = 0.25\nT = 0.25",
            "GPS-CONF-006",
        ),
        (
            "legacy-category-field",
            "A = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\ncategory = \"base\"",
            "GPS-CONF-006",
        ),
        (
            "duplicate-a",
            "A = 0.5\nA = 0.5\nC = 0.0\nG = 0.0\nT = 0.0",
            "GPS-CONF-007",
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
[reference]
kind = "simulated"
id = "ref-invalid"
[reference.composition]
{composition}
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
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(expected_code),
            "{label}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("cleanup");
    }
}

#[test]
fn generated_composition_accepts_an_independent_n_fraction() {
    let root = temporary("generated-n");
    fs::create_dir(&root).expect("root");
    let config = root.join("general.toml");
    fs::write(
        &config,
        r#"config_kind = "general"
schema_version = "1"
[reference]
kind = "simulated"
id = "ref-n"
[reference.composition]
A = 1.0
C = 0.0
G = 0.0
T = 0.0
N = 1.0
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 8
"#,
    )
    .expect("config");
    let destination = root.join("result");
    run_cli(&config, &destination, true);
    assert_eq!(
        fs::read(destination.join("resources/references/ref-n.fa")).expect("generated FASTA"),
        b">chr1\nNNNNNNNN\n"
    );
    fs::remove_dir_all(root).expect("cleanup");
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
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[reference]
kind = "supplied"
id = "ref-supplied"
fasta = "{fasta}"

[[reference.sequences]]
id = "chrN"
class = "chromosome"
length = {length}
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
