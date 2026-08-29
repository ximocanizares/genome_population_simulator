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

    let version = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .arg("--version")
        .output()
        .expect("el binario de integración debe informar su versión");
    assert_eq!(
        String::from_utf8(version.stdout)
            .expect("versión UTF-8")
            .trim(),
        "genome-population-simulator 0.2.0"
    );
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
schema_version = "2"
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
fn cli_publishes_founder_keys_when_it_derives_initial_individual_ids() {
    let root = std::env::temp_dir().join(format!(
        "gps-founder-keys-cli-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("temporary root");
    let fasta = root.join("reference.fa");
    let matrix = root.join("matrix.tsv");
    fs::write(&fasta, b">chr1\nACGT\n").expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t1\tA\tG\tP1\t0.5\n",
    )
    .expect("matrix");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "2"
seed = "{}"
[population_snvs]
mode = "supplied"
input_path = "{}"
[[initial_populations]]
population_id = "P1"
individual_count = 1
founder_keys = ["founder-I1"]
ploidy = 2
genotype_model = "hwe"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 4
"#,
            "0".repeat(64),
            matrix.display(),
            fasta.display(),
        ),
    )
    .expect("configuration");
    let output = root.join("result");
    let result = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args([
            "run",
            config.to_str().expect("path"),
            "--output",
            output.to_str().expect("path"),
        ])
        .output()
        .expect("CLI run");
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let effective: serde_json::Value = serde_json::from_slice(
        &fs::read(output.join("effective-config.json")).expect("effective configuration"),
    )
    .expect("JSON");
    assert_eq!(
        effective["initial_populations"][0]["founder_keys"],
        serde_json::json!(["founder-I1"])
    );
    assert_ne!(
        effective["initial_populations"][0]["individual_ids"][0], "founder-I1",
        "the published ID is derived from population and founder key"
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn normalized_multiallelic_matrix_keeps_zero_ref_through_validate_and_run() {
    let root = std::env::temp_dir().join(format!(
        "gps-normalized-multiallelic-cli-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    if root.exists() {
        fs::remove_dir_all(&root).expect("stale temporary root cleanup");
    }
    fs::create_dir_all(&root).expect("temporary root");
    let fasta = root.join("reference.fa");
    let matrix = root.join("matrix.tsv");
    fs::write(&fasta, b">chr1\nA\n").expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\n\
         chr1\t1\tA\tC\tP1\t0.000001\n\
         chr1\t1\tA\tG\tP1\t0.002992\n\
         chr1\t1\tA\tT\tP1\t0.997007001\n",
    )
    .expect("matrix");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "2"
seed = "{}"
normalization_tolerance = 0.000001
[population_snvs]
mode = "supplied"
input_path = "{}"
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1
"#,
            "1".repeat(64),
            matrix.display(),
            fasta.display(),
        ),
    )
    .expect("configuration");
    let binary = env!("CARGO_BIN_EXE_genome-population-simulator");
    let validation = Command::new(binary)
        .args(["validate", config.to_str().expect("path")])
        .output()
        .expect("CLI validate");
    assert!(
        validation.status.success(),
        "{}",
        String::from_utf8_lossy(&validation.stderr)
    );
    let output = root.join("result");
    let run = Command::new(binary)
        .args([
            "run",
            config.to_str().expect("path"),
            "--output",
            output.to_str().expect("path"),
        ])
        .output()
        .expect("CLI run");
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let effective: serde_json::Value = serde_json::from_slice(
        &fs::read(output.join("effective-config.json")).expect("effective configuration"),
    )
    .expect("JSON");
    assert_eq!(
        effective["population_snvs"]["normalization_audit"][0]["reference_frequency"],
        0.0
    );
    let locus =
        &effective["population_snvs"]["initial_population_statistics"]["populations"][0]["loci"][0];
    assert_eq!(locus["alleles"][0]["requested_frequency"], 0.0);
    assert_eq!(locus["alleles"][0]["realized_count"], 0);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn globally_absent_alt_is_omitted_from_active_projection_for_supplied_and_reused_cli() {
    let root = std::env::temp_dir().join(format!(
        "gps-zero-alt-projection-cli-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    if root.exists() {
        fs::remove_dir_all(&root).expect("stale temporary root cleanup");
    }
    fs::create_dir_all(&root).expect("temporary root");
    let fasta = root.join("reference.fa");
    let matrix = root.join("matrix.tsv");
    fs::write(&fasta, b">chr1\nA\n").expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\n\
         chr1\t1\tA\tC\tP1\t0\n\
         chr1\t1\tA\tG\tP1\t0.5\n",
    )
    .expect("matrix");
    let supplied = root.join("supplied.toml");
    let common = format!(
        r#"config_kind = "general"
schema_version = "2"
seed = "{}"
[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["I1"]
ploidy = 2
genotype_model = "hwe"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1
"#,
        "2".repeat(64),
        fasta.display(),
    );
    fs::write(
        &supplied,
        format!(
            "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"{}\"\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"{}\"\n{}",
            "2".repeat(64),
            matrix.display(),
            common
                .splitn(4, '\n')
                .skip(3)
                .collect::<Vec<_>>()
                .join("\n"),
        ),
    )
    .expect("supplied config");
    let binary = env!("CARGO_BIN_EXE_genome-population-simulator");
    let supplied_output = root.join("supplied-result");
    for arguments in [
        vec!["validate".to_owned(), supplied.display().to_string()],
        vec![
            "run".to_owned(),
            supplied.display().to_string(),
            "--output".to_owned(),
            supplied_output.display().to_string(),
        ],
    ] {
        let result = Command::new(binary).args(arguments).output().expect("CLI");
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    let effective: serde_json::Value = serde_json::from_slice(
        &fs::read(supplied_output.join("effective-config.json")).expect("effective"),
    )
    .expect("JSON");
    let locus =
        &effective["population_snvs"]["initial_population_statistics"]["populations"][0]["loci"][0];
    assert_eq!(locus["alleles"].as_array().expect("alleles").len(), 2);
    assert_eq!(locus["alleles"][1]["allele"], "G");

    let reused = root.join("reused.toml");
    fs::write(
        &reused,
        format!(
            "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"{}\"\n\
             [population_snvs]\nmode = \"reuse\"\nmetadata_path = \"{}\"\n{}",
            "2".repeat(64),
            supplied_output
                .join("population_snvs.metadata.json")
                .display(),
            common
                .splitn(4, '\n')
                .skip(3)
                .collect::<Vec<_>>()
                .join("\n"),
        ),
    )
    .expect("reuse config");
    let reused_output = root.join("reuse-result");
    for arguments in [
        vec!["validate".to_owned(), reused.display().to_string()],
        vec![
            "run".to_owned(),
            reused.display().to_string(),
            "--output".to_owned(),
            reused_output.display().to_string(),
        ],
    ] {
        let result = Command::new(binary).args(arguments).output().expect("CLI");
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    let metadata_path = supplied_output.join("population_snvs.metadata.json");
    let mut hostile: serde_json::Value =
        serde_json::from_slice(&fs::read(&metadata_path).expect("published metadata"))
            .expect("JSON");
    hostile["normalization_audit"][0]["values"][0]["effective"] = serde_json::json!(0.1);
    fs::write(
        &metadata_path,
        serde_json::to_vec(&hostile).expect("hostile JSON"),
    )
    .expect("hostile metadata");
    let rejected = Command::new(binary)
        .args(["validate", reused.to_str().expect("path")])
        .output()
        .expect("hostile reuse validation");
    assert!(
        !rejected.status.success(),
        "an unmatched nonzero discarded allele must be rejected"
    );
    fs::remove_dir_all(root).expect("cleanup");
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
schema_version = "2"
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
    let statistics = &effective["population_snvs"]["initial_population_statistics"];
    assert_eq!(statistics["kind"], "initial_population_statistics");
    assert!(
        statistics["populations"][0]["loci"][0]["observed_copies"]
            .as_u64()
            .is_some()
    );
    fs::remove_dir_all(root).expect("cleanup");
}
