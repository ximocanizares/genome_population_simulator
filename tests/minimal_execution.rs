use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn minimal_cli_never_publishes_resolution_evidence_or_event_units() {
    let root = temporary();
    fs::create_dir(&root).expect("root");
    let destination = root.join("result");
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args([
            "run",
            "tests/fixtures/config/general-seeded.toml",
            "--output",
        ])
        .arg(&destination)
        .output()
        .expect("CLI");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!destination.join("inputs/resolution/evidence.json").exists());
    let paths = collect_paths(&destination);
    assert!(paths.iter().all(|path| {
        let text = path.to_string_lossy();
        !text.contains("resolution") && !text.contains("event")
    }));
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifest"))
            .expect("JSON");
    assert!(manifest.get("resolution_evidence").is_none());
    assert!(manifest.get("events").is_none());
    assert!(
        manifest["units"]
            .as_array()
            .expect("units")
            .iter()
            .all(|unit| unit["kind"] != "event")
    );
    fs::remove_dir_all(root).expect("cleanup");
}

fn temporary() -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "gps-minimal-execution-{}-{}-{nanos}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn collect_paths(root: &Path) -> Vec<PathBuf> {
    let mut pending = vec![root.to_path_buf()];
    let mut result = Vec::new();
    while let Some(directory) = pending.pop() {
        for entry in fs::read_dir(directory).expect("read tree") {
            let path = entry.expect("entry").path();
            if path.is_dir() {
                pending.push(path.clone());
            }
            result.push(path);
        }
    }
    result
}
