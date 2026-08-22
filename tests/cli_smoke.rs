use std::fs;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn help_and_version_are_successful_cli_queries() {
    for argument in ["--help", "-h", "--version"] {
        let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
            .arg(argument)
            .output()
            .expect("el binario de integración debe arrancar");

        assert!(
            output.status.success(),
            "{argument} terminó con {:?}: {}",
            output.status.code(),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(!output.stdout.is_empty(), "{argument} no produjo salida");
    }
}

#[test]
fn empty_population_matrix_requires_and_records_its_explicit_authorization() {
    let root = std::env::temp_dir().join(format!(
        "gps-empty-matrix-cli-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    if root.exists() {
        fs::remove_dir_all(&root).expect("stale test root cleanup");
    }
    fs::create_dir_all(&root).expect("temporary root");
    let matrix = root.join("empty.tsv");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\n",
    )
    .expect("empty canonical matrix");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "1"
[population_snvs]
mode = "supplied"
input_path = "{}"
[reference]
kind = "simulated"
id = "ref-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1
"#,
            matrix.display()
        ),
    )
    .expect("configuration");
    let binary = env!("CARGO_BIN_EXE_genome-population-simulator");
    let rejected = Command::new(binary)
        .args(["run", config.to_str().expect("path"), "--output"])
        .arg(root.join("rejected"))
        .output()
        .expect("run without authorization");
    assert!(!rejected.status.success());

    let output = root.join("authorized");
    let authorized = Command::new(binary)
        .args([
            "run",
            config.to_str().expect("path"),
            "--allow-empty-snv-matrix",
            "--output",
        ])
        .arg(&output)
        .output()
        .expect("run with authorization");
    assert!(
        authorized.status.success(),
        "{}",
        String::from_utf8_lossy(&authorized.stderr)
    );
    for name in ["effective-config.json", "manifest.json"] {
        let document: serde_json::Value =
            serde_json::from_slice(&fs::read(output.join(name)).expect("recorded output"))
                .expect("JSON");
        assert_eq!(document["authorizations"]["allow_empty_snv_matrix"], true);
    }
    fs::remove_dir_all(root).expect("cleanup");
}
