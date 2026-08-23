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
            "--allow-incompatible-genotypes",
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
        assert_eq!(
            document["authorizations"]["allow_incompatible_genotypes"],
            true
        );
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn incompatible_genotype_authorization_is_named_and_force_is_not_a_cli_option() {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/config/general-minimal.toml");
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["validate", fixture.to_str().expect("UTF-8 path"), "--force"])
        .output()
        .expect("CLI invocation");
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("argumento desconocido"));
}

#[test]
fn direct_events_are_preflighted_before_publication_and_the_narrow_flag_is_traced() {
    let root = std::env::temp_dir().join(format!(
        "gps-direct-events-cli-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("temporary root");
    let fasta = root.join("reference.fa");
    let matrix = root.join("matrix.tsv");
    fs::write(&fasta, b">chr1\nACGT\n").expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t1\tA\tG\tP1\t0.5\nchr1\t1\tA\tT\tP1\t0.5\n",
    )
    .expect("matrix");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[population_snvs]
mode = "supplied"
input_path = "{}"
[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["I1"]
ploidy = 2
genotype_model = "f_diploid"
inbreeding_coefficient = 1
[[direct_events]]
event_id = "force-g"
individual_id = "I1"
phase = 0
kind = "snp"
contig = "chr1"
position_1based = 1
ref = "A"
alt = "G"
[[direct_events]]
event_id = "force-t"
individual_id = "I1"
phase = 1
kind = "snp"
contig = "chr1"
position_1based = 1
ref = "A"
alt = "T"
[[direct_events]]
event_id = "insert-c"
individual_id = "I1"
phase = 0
kind = "insertion"
contig = "chr1"
point = 1
bases = "C"
[[direct_events]]
event_id = "delete-current"
individual_id = "I1"
phase = 0
kind = "deletion"
contig = "chr1"
start = 1
end = 2
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 4
"#,
            matrix.display(),
            fasta.display(),
        ),
    )
    .expect("configuration");
    let binary = env!("CARGO_BIN_EXE_genome-population-simulator");
    let validation = Command::new(binary)
        .args(["validate", config.to_str().expect("path")])
        .output()
        .expect("rejected validation");
    assert!(!validation.status.success());
    let validation = Command::new(binary)
        .args([
            "validate",
            config.to_str().expect("path"),
            "--allow-incompatible-genotypes",
        ])
        .output()
        .expect("authorized validation");
    assert!(
        validation.status.success(),
        "{}",
        String::from_utf8_lossy(&validation.stderr)
    );
    let rejected = root.join("rejected");
    let output = Command::new(binary)
        .args(["run", config.to_str().expect("path"), "--output"])
        .arg(&rejected)
        .output()
        .expect("rejected run");
    assert!(!output.status.success());
    assert!(
        !rejected.exists(),
        "preflight must precede workspace creation"
    );

    let authorized = root.join("authorized");
    let output = Command::new(binary)
        .args([
            "run",
            config.to_str().expect("path"),
            "--allow-incompatible-genotypes",
            "--output",
        ])
        .arg(&authorized)
        .output()
        .expect("authorized run");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    for name in ["effective-config.json", "manifest.json"] {
        let document: serde_json::Value =
            serde_json::from_slice(&fs::read(authorized.join(name)).expect("published JSON"))
                .expect("JSON");
        assert_eq!(
            document["authorizations"]["allow_incompatible_genotypes"],
            true
        );
    }
    let effective: serde_json::Value = serde_json::from_slice(
        &fs::read(authorized.join("effective-config.json")).expect("effective configuration"),
    )
    .expect("JSON");
    assert_eq!(effective["direct_events"][0]["event_id"], "force-g");
    fs::remove_dir_all(root).expect("cleanup");
}
