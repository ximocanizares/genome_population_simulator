use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use sha2::{Digest, Sha256};

use genome_population_simulator::application::load_and_validate_f5_config;

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

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn write_f5_profile(profile_dir: &Path) -> String {
    fs::create_dir_all(profile_dir).expect("profile directory");
    let files = [
        (
            "block_lengths.tsv",
            "f5_block_length_profile",
            "length_min_bp\tlength_max_bp\tweight\n2\t10\t1\n",
        ),
        (
            "snp_gaps.tsv",
            "f5_snp_gap_profile",
            "gap_min_bp\tgap_max_bp\tweight\n1\t4\t1\n",
        ),
        (
            "snp_spectrum.tsv",
            "f5_snp_spectrum",
            "alt_af_min\talt_af_max\tweight\n0\t1\t1\n",
        ),
        (
            "haplotype_weights.tsv",
            "f5_haplotype_weights",
            "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        ),
        (
            "ld_profile.tsv",
            "f5_ld_profile",
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t9\t1\t0.05\n",
        ),
    ];
    let links = files
        .iter()
        .map(|(name, kind, content)| {
            fs::write(profile_dir.join(name), content).expect("profile TSV");
            format!(
                r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#,
                content.len(),
                sha256(content.as_bytes())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let metadata = format!(
        r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"test_profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
    );
    fs::write(profile_dir.join("profile.metadata.json"), &metadata).expect("profile metadata");
    sha256(metadata.as_bytes())
}

fn write_f5_cli_fixture(root: &Path) -> PathBuf {
    let scenario_dir = root.join("authoring/scenarios");
    let profile_hash = write_f5_profile(&scenario_dir.join("profiles"));
    fs::create_dir_all(&scenario_dir).expect("scenario directory");
    let reference_dir = root.join("authoring/references");
    fs::create_dir_all(&reference_dir).expect("reference directory");
    let reference = b">chr1\nAAAAAAAAAA\n";
    fs::write(reference_dir.join("ref.fa"), reference).expect("reference FASTA");
    let general_path = root.join("authoring/general.toml");
    fs::write(
        &general_path,
        format!(
            "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
            "a".repeat(64), sha256(reference)
        ),
    )
    .expect("general");
    fs::write(
        scenario_dir.join("scenario.toml"),
        format!(
            "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = 1\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{profile_hash}\"\ncoverage_fraction = 1.0\n\n[[subpopulations]]\npopulation_id = \"P1\"\nindividual_count = 1\nfounder_ids = [\"I1\"]\n"
        ),
    )
    .expect("scenario");
    general_path
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

#[test]
fn cli_dispatches_general_v3_from_an_unrelated_three_level_working_directory() {
    let root = destination("f5-cli-relative");
    let foreign_cwd = root.join("elsewhere/a/b");
    fs::create_dir_all(&foreign_cwd).expect("foreign cwd");
    let general = write_f5_cli_fixture(&root);
    let effective = load_and_validate_f5_config(&general).expect("effective F5 config");
    assert_eq!(effective.profile().profile_id(), "test_profile");
    assert_eq!(effective.profile().normalizations().len(), 4);
    let relative = Path::new("../../../authoring/general.toml");
    let validate = Command::new(BINARY)
        .current_dir(&foreign_cwd)
        .args(["validate", relative.to_str().expect("relative path")])
        .output()
        .expect("F5 validate");
    assert!(
        validate.status.success(),
        "{}",
        String::from_utf8_lossy(&validate.stderr)
    );
    let output = root.join("out");
    let run = Command::new(BINARY)
        .current_dir(&foreign_cwd)
        .args(["run", relative.to_str().expect("relative path"), "--output"])
        .arg(&output)
        .output()
        .expect("F5 run");
    // F5.6.1 encadena `run`: el mismo par, invocado por una ruta relativa
    // desde un directorio ajeno, publica su paquete en el destino pedido.
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    for relative_artifact in [
        "manifest.json",
        "effective-config.json",
        "outputs/f5-block-layout.tsv",
        "outputs/f5-sites.tsv",
        "outputs/f5-haplotypes.tsv",
        "outputs/f5-frequencies.tsv",
        "outputs/f5-ld-evaluation.tsv",
        "outputs/f5-founders.tsv",
        "outputs/f5-attempts.tsv",
        "outputs/truth.jsonl",
    ] {
        assert!(
            output.join(relative_artifact).is_file(),
            "falta {relative_artifact} en el paquete publicado"
        );
    }
    // El nombre del general se publica: la ruta de invocación es relativa y
    // ajena, y aun así `source_path` es exactamente el nombre del documento.
    let effective_document =
        fs::read_to_string(output.join("effective-config.json")).expect("efectivo publicado");
    assert!(
        effective_document.contains("\"source_path\": \"general.toml\""),
        "{effective_document}"
    );
    // DTC-045: el destino no se sobrescribe ni se mezcla.
    let again = Command::new(BINARY)
        .current_dir(&foreign_cwd)
        .args(["run", relative.to_str().expect("relative path"), "--output"])
        .arg(&output)
        .output()
        .expect("F5 run repetido");
    assert!(!again.status.success());
    assert!(
        String::from_utf8_lossy(&again.stderr).contains("GPS-IO-003"),
        "{}",
        String::from_utf8_lossy(&again.stderr)
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn f5_validate_rejects_missing_corrupt_or_mismatched_profile_before_reference() {
    for case in ["absent", "corrupt", "hash"] {
        let root = destination(case);
        let general = write_f5_cli_fixture(&root);
        let profile = root.join("authoring/scenarios/profiles/profile.metadata.json");
        let scenario = root.join("authoring/scenarios/scenario.toml");
        match case {
            "absent" => fs::remove_file(&profile).expect("remove profile"),
            "corrupt" => {
                let corrupt = b"not-json";
                fs::write(&profile, corrupt).expect("corrupt profile");
                let source = fs::read_to_string(&scenario).expect("scenario source");
                let old_hash = source
                    .lines()
                    .find_map(|line| {
                        line.strip_prefix("profile_metadata_sha256 = \"")?
                            .strip_suffix('"')
                    })
                    .expect("declared hash");
                fs::write(&scenario, source.replace(old_hash, &sha256(corrupt)))
                    .expect("updated scenario");
            }
            "hash" => {
                let source = fs::read_to_string(&scenario).expect("scenario source");
                let old_hash = source
                    .lines()
                    .find_map(|line| {
                        line.strip_prefix("profile_metadata_sha256 = \"")?
                            .strip_suffix('"')
                    })
                    .expect("declared hash");
                fs::write(&scenario, source.replace(old_hash, &"0".repeat(64)))
                    .expect("wrong scenario hash");
            }
            _ => unreachable!(),
        }
        let output = root.join("out");
        let result = run("validate", &general, "human", &output);
        assert_eq!(
            result.status.code(),
            Some(2),
            "{case}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(!output.exists());
        fs::remove_dir_all(root).expect("cleanup");
    }
}

#[test]
fn f5_validate_executes_physical_admission_without_destination() {
    let root = destination("f5-cli-preflight");
    let general = write_f5_cli_fixture(&root);
    fs::write(
        root.join("authoring/references/ref.fa"),
        b">chr1\nAAAAAAAATA\n",
    )
    .expect("mutated reference");
    let output = root.join("must-not-exist");
    let result = run("validate", &general, "human", &output);
    assert_eq!(result.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&result.stderr).contains("GPS-INPUT-002"));
    assert!(!output.exists());

    let general = write_f5_cli_fixture(&root);
    let scenario = root.join("authoring/scenarios/scenario.toml");
    let source = fs::read_to_string(&scenario).expect("scenario source");
    fs::write(
        &scenario,
        source.replace(
            "coverage_fraction = 1.0",
            "coverage_fraction = 1.0\nmax_operation_units = 1",
        ),
    )
    .expect("resource limit");
    let result = run("validate", &general, "json", &output);
    assert_eq!(result.status.code(), Some(2));
    assert!(!output.exists());
    fs::remove_dir_all(root).expect("cleanup");
}
