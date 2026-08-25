use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::publication::ExecutionWorkspace;
use genome_population_simulator::results::{
    ConfigurationInputRecord, GeneralManifest, UnitState, general_manifest_document,
    json_document_bytes, validate_general_manifest_bytes,
};

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn configuration_copy_is_byte_exact_linked_and_matches_external_sha256() {
    let root = temporary("configuration-copy");
    fs::create_dir(&root).unwrap();
    let destination = root.join("result");
    let source = Path::new("tests/fixtures/config/general-seeded.toml");
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run"])
        .arg(source)
        .arg("--output")
        .arg(&destination)
        .output()
        .expect("complete process");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let copied = destination.join("inputs/config/general.toml");
    assert_eq!(fs::read(source).unwrap(), fs::read(&copied).unwrap());

    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["configuration_input"]["kind"], "general");
    assert_eq!(
        manifest["configuration_input"]["original_path"],
        source.to_string_lossy().as_ref()
    );
    let artifact_id = manifest["configuration_input"]["artifact_id"]
        .as_str()
        .unwrap();
    let artifact = manifest["artifacts"]
        .as_array()
        .unwrap()
        .iter()
        .find(|artifact| artifact["id"] == artifact_id)
        .expect("configuration_input resolves exactly once");
    assert_eq!(artifact["kind"], "config_general_toml");
    assert_eq!(artifact["path"], "inputs/config/general.toml");
    assert_eq!(
        artifact["size_bytes"].as_u64().unwrap(),
        fs::metadata(&copied).unwrap().len()
    );
    assert_eq!(artifact["checksum"]["algorithm"], "sha256");
    let sha256sum = Command::new("sha256sum").arg(&copied).output().unwrap();
    assert!(sha256sum.status.success());
    let digest = String::from_utf8(sha256sum.stdout)
        .unwrap()
        .split_whitespace()
        .next()
        .unwrap()
        .to_owned();
    assert_eq!(artifact["checksum"]["digest"], digest);

    fs::remove_dir_all(root).unwrap();
}

#[test]
fn catalog_files_are_unique_resolvable_and_match_external_sha256() {
    let root = temporary("catalog");
    fs::create_dir(&root).unwrap();
    let destination = root.join("result");
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args([
            "run",
            "tests/fixtures/config/general-seeded.toml",
            "--output",
        ])
        .arg(&destination)
        .output()
        .expect("complete process");
    assert!(output.status.success());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).unwrap()).unwrap();
    let artifacts = manifest["artifacts"].as_array().unwrap();
    let mut ids = std::collections::BTreeSet::new();
    let mut paths = std::collections::BTreeSet::new();
    for artifact in artifacts {
        assert!(ids.insert(artifact["id"].as_str().unwrap()));
        assert!(paths.insert(artifact["path"].as_str().unwrap()));
        let path = destination.join(artifact["path"].as_str().unwrap());
        assert_eq!(
            fs::metadata(&path).unwrap().len(),
            artifact["size_bytes"].as_u64().unwrap()
        );
        let output = Command::new("sha256sum").arg(&path).output().unwrap();
        assert!(output.status.success());
        let digest = String::from_utf8(output.stdout)
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .to_owned();
        assert_eq!(artifact["checksum"]["algorithm"], "sha256");
        assert_eq!(artifact["checksum"]["digest"], digest);
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn complete_failed_and_incomplete_are_real_minimal_trees() {
    let root = temporary("three-trees");
    fs::create_dir(&root).unwrap();

    let complete = root.join("complete");
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args([
            "run",
            "tests/fixtures/config/general-seeded.toml",
            "--output",
        ])
        .arg(&complete)
        .output()
        .expect("complete process");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_final_tree(&complete, "complete");

    let failed = root.join("failed");
    let (failed_config, failed_reference) = write_slow_supplied_config(&root);
    let failed_child = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run"])
        .arg(&failed_config)
        .arg("--output")
        .arg(&failed)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("real failed execution");
    wait_for_provisionals(&failed, 1);
    fs::rename(&failed_reference, root.join("reference-removed.fa"))
        .expect("remove input after preflight");
    let failed_output = failed_child.wait_with_output().expect("failed child");
    assert_eq!(failed_output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&failed_output.stderr).contains("GPS-INPUT-001"));
    assert_final_tree(&failed, "failed");
    let failed_manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(failed.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(failed_manifest["diagnostics"][0]["code"], "GPS-INPUT-001");
    let failed_reductions = failed_manifest["reference_reductions"]
        .as_array()
        .expect("failed reductions");
    assert!(
        failed_reductions.is_empty(),
        "the source copy failed before any reduction selection was reached"
    );
    assert!(failed_manifest["references"].as_array().unwrap().is_empty());
    assert_eq!(failed_manifest["artifacts"].as_array().unwrap().len(), 1);
    assert_eq!(failed_manifest["artifacts"][0]["id"], "config-general");
    assert_eq!(
        failed_manifest["output_categories"],
        serde_json::json!({
            "minimum": ["config-general"],
            "defaulted": [],
            "requested": [],
            "disabled": []
        })
    );

    let interrupted = root.join("interrupted");
    let mut child = child("wait", &interrupted, None)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .expect("interruptible child");
    let provisional = wait_for_provisionals(&interrupted, 1).remove(0);
    let kill = Command::new("kill")
        .args(["-KILL", &child.id().to_string()])
        .status()
        .expect("kill child");
    assert!(kill.success());
    let status = child.wait().expect("reap killed child");
    assert!(!status.success());
    assert!(!interrupted.exists());
    assert!(provisional.join("INCOMPLETE").is_file());
    assert_eq!(
        fs::metadata(&provisional).unwrap().permissions().mode() & 0o777,
        0o700
    );

    fs::remove_dir_all(root).unwrap();
}

#[test]
fn destination_is_never_clobbered_and_two_processes_have_one_winner() {
    let root = temporary("race");
    fs::create_dir(&root).unwrap();

    let existing = root.join("existing");
    fs::create_dir(&existing).unwrap();
    fs::write(existing.join("sentinel"), b"preserve").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args([
            "run",
            "tests/fixtures/config/general-seeded.toml",
            "--output",
        ])
        .arg(&existing)
        .output()
        .expect("existing destination process");
    assert_eq!(output.status.code(), Some(5));
    assert_eq!(fs::read(existing.join("sentinel")).unwrap(), b"preserve");

    let destination = root.join("winner");
    let barrier = root.join("go");
    let first = child("race", &destination, Some(&barrier))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("first racer");
    let second = child("race", &destination, Some(&barrier))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("second racer");
    let provisionals = wait_for_provisionals(&destination, 2);
    fs::write(&barrier, b"go\n").unwrap();
    let outputs = [
        first.wait_with_output().unwrap(),
        second.wait_with_output().unwrap(),
    ];
    assert_eq!(
        outputs
            .iter()
            .filter(|output| output.status.success())
            .count(),
        1
    );
    assert_eq!(
        outputs
            .iter()
            .filter(|output| output.status.code() == Some(5))
            .count(),
        1
    );
    let loser = outputs
        .iter()
        .find(|output| !output.status.success())
        .unwrap();
    assert!(String::from_utf8_lossy(&loser.stderr).contains("GPS-IO-003"));
    assert_final_tree(&destination, "failed");
    assert_eq!(
        provisionals.iter().filter(|path| path.exists()).count(),
        1,
        "the losing process keeps its marked provisional"
    );
    assert!(
        provisionals
            .iter()
            .filter(|path| path.exists())
            .all(|path| path.join("INCOMPLETE").is_file())
    );

    fs::remove_dir_all(root).unwrap();
}

#[test]
fn invalid_configuration_exits_two_without_destination_or_provisional() {
    let root = temporary("exit-two");
    fs::create_dir(&root).unwrap();
    let destination = root.join("rejected");
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args([
            "run",
            "tests/fixtures/config/invalid-version-two.toml",
            "--output",
        ])
        .arg(&destination)
        .output()
        .expect("invalid configuration process");
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("GPS-CONF-005"));
    assert!(!destination.exists());
    assert!(wait_for_no_provisionals(&destination));
    fs::remove_dir_all(root).unwrap();
}

#[test]
#[ignore = "subprocess helper"]
fn process_child() {
    let mode = std::env::var("GPS_F14_CHILD_MODE").expect("child mode");
    let destination = PathBuf::from(std::env::var_os("GPS_F14_DESTINATION").expect("destination"));
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    match mode.as_str() {
        "race" => {
            let diagnostic = controlled_diagnostic();
            write_failed_manifest(&mut workspace, None, &diagnostic);
            let barrier = PathBuf::from(std::env::var_os("GPS_F14_BARRIER").expect("barrier"));
            while !barrier.exists() {
                std::thread::sleep(Duration::from_millis(5));
            }
            match workspace.publish() {
                Ok(_) => {}
                Err(diagnostic) => {
                    eprintln!("{}", diagnostic.render_human());
                    std::process::exit(diagnostic.exit_code().into());
                }
            }
        }
        "wait" => loop {
            std::thread::park();
        },
        other => panic!("unknown child mode {other}"),
    }
}

fn write_slow_supplied_config(root: &Path) -> (PathBuf, PathBuf) {
    let reference = root.join("reference.fa");
    let mut fasta = String::from(">chr1\n");
    for _ in 0..20 {
        fasta.push_str(&"A".repeat(50));
        fasta.push('\n');
    }
    fs::write(&reference, fasta).expect("reference");
    let config = root.join("slow-general.toml");
    let mut text = format!(
        "config_kind = \"general\"\nschema_version = \"1\"\nseed = \"{}\"\n\
         [reference]\nkind = \"supplied\"\nid = \"reference-main\"\nfasta = \"reference.fa\"\n\
         [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1000\n\
         [[reduced_references]]\nid = \"reference-small\"\nsource_reference_id = \"reference-main\"\nregion_size_bases = 10\ntarget_bases = 100\n",
        "1".repeat(64)
    );
    // Tras crear el provisional, el orquestador copia estos bytes antes de
    // abrir la referencia; da una ventana determinista al proceso padre.
    text.push_str(&"# padding para fallo controlado\n".repeat(750_000));
    fs::write(&config, text).expect("large config");
    (config, reference)
}

fn controlled_diagnostic() -> Diagnostic {
    validate_general_manifest_bytes(Path::new("controlled-failure.json"), b"{")
        .expect_err("truncated JSON produces a controlled diagnostic")
}

fn write_failed_manifest(
    workspace: &mut ExecutionWorkspace,
    configuration: Option<&ConfigurationInputRecord>,
    diagnostic: &Diagnostic,
) {
    let artifacts = workspace.artifacts().to_vec();
    let manifest = GeneralManifest {
        execution_id: "execution-process-probe",
        state: UnitState::Failed,
        configuration_input: configuration,
        effective_config: None,
        references: &[],
        reduced_references: &[],
        population_snvs: None,
        vcf_outputs: &[],
        authorizations: genome_population_simulator::execution::ExecutionAuthorizations::default(),
        artifacts: &artifacts,
        diagnostics: std::slice::from_ref(diagnostic),
    };
    let bytes = json_document_bytes(&general_manifest_document(&manifest)).expect("manifest bytes");
    workspace
        .write_general_manifest(&bytes)
        .expect("write failed manifest");
}

fn child(mode: &str, destination: &Path, barrier: Option<&Path>) -> Command {
    let mut command = Command::new(std::env::current_exe().expect("current test executable"));
    command
        .args(["--ignored", "--exact", "process_child", "--nocapture"])
        .env("GPS_F14_CHILD_MODE", mode)
        .env("GPS_F14_DESTINATION", destination);
    if let Some(barrier) = barrier {
        command.env("GPS_F14_BARRIER", barrier);
    }
    command
}

fn wait_for_provisionals(destination: &Path, count: usize) -> Vec<PathBuf> {
    let deadline = Instant::now() + Duration::from_secs(10);
    loop {
        let parent = destination.parent().unwrap();
        let prefix = format!(
            ".{}.gps-incomplete-",
            destination.file_name().unwrap().to_string_lossy()
        );
        let mut paths = fs::read_dir(parent)
            .unwrap()
            .filter_map(Result::ok)
            .filter(|entry| entry.file_name().to_string_lossy().starts_with(&prefix))
            .map(|entry| entry.path())
            .filter(|path| path.join("INCOMPLETE").is_file())
            .collect::<Vec<_>>();
        paths.sort();
        if paths.len() >= count {
            return paths;
        }
        assert!(
            Instant::now() < deadline,
            "timeout waiting for provisionals"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
}

fn wait_for_no_provisionals(destination: &Path) -> bool {
    let parent = destination.parent().unwrap();
    let prefix = format!(
        ".{}.gps-incomplete-",
        destination.file_name().unwrap().to_string_lossy()
    );
    !fs::read_dir(parent)
        .unwrap()
        .filter_map(Result::ok)
        .any(|entry| entry.file_name().to_string_lossy().starts_with(&prefix))
}

fn assert_final_tree(destination: &Path, expected_state: &str) {
    assert!(destination.is_dir());
    assert!(!destination.join("INCOMPLETE").exists());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["state"], expected_state);
    assert!(manifest.get("units").is_none());
    assert!(manifest.get("events").is_none());
    assert!(manifest.get("sets").is_none());
    for artifact in manifest["artifacts"].as_array().unwrap() {
        assert!(
            destination
                .join(artifact["path"].as_str().unwrap())
                .is_file()
        );
    }
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f14-process-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
