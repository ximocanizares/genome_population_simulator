use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{
    CliValidatedConfig, load_and_validate_cli_config,
    load_and_validate_cli_config_with_bootstrap_open_observer, load_and_validate_config,
    load_and_validate_config_with_entropy, load_and_validate_f5_config,
    load_and_validate_f5_config_with_entropy,
};
use genome_population_simulator::randomness::EntropySource;
use sha2::{Digest, Sha256};

const TOML_LIMIT: usize = 65_536;
static NONCE: AtomicU64 = AtomicU64::new(0);

#[derive(Debug)]
struct PanicEntropy;

impl EntropySource for PanicEntropy {
    fn fill_root_seed(&mut self, _: &mut [u8; 32]) -> Result<(), String> {
        panic!("el bootstrap rechazado no puede abrir RNG")
    }
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-bootstrap-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn write_profile(directory: &Path) -> String {
    fs::create_dir_all(directory).expect("profile directory");
    let files = [
        (
            "block_lengths.tsv",
            "f5_block_length_profile",
            "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
        ),
        (
            "snp_gaps.tsv",
            "f5_snp_gap_profile",
            "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n",
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
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
        ),
    ];
    let links = files
        .iter()
        .map(|(name, kind, contents)| {
            fs::write(directory.join(name), contents).expect("profile table");
            format!(
                r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#,
                contents.len(),
                sha256(contents.as_bytes())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let metadata = format!(
        r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
    );
    fs::write(directory.join("profile.metadata.json"), &metadata).expect("profile metadata");
    sha256(metadata.as_bytes())
}

fn general() -> String {
    let fasta = b">chr1\nAA\n";
    format!(
        "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
        "a".repeat(64),
        sha256(fasta)
    )
}

fn scenario(profile_hash: &str, memory_limit: &str) -> String {
    format!(
        "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = 1\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{profile_hash}\"\nmax_working_memory_bytes = {memory_limit}\nblock_count = 1\n\n[[subpopulations]]\npopulation_id = \"P1\"\nindividual_count = 1\nfounder_ids = [\"I1\"]\n"
    )
}

fn pad_to_toml_limit(mut source: String) -> String {
    assert!(source.len() <= TOML_LIMIT);
    source.push_str(&"#".repeat(TOML_LIMIT - source.len()));
    source
}

fn write_fixture(root: &Path, create_profile: bool, scenario_source: String) -> PathBuf {
    let authoring = root.join("authoring");
    let scenarios = authoring.join("scenarios");
    let _profile_hash = if create_profile {
        write_profile(&scenarios.join("profiles"))
    } else {
        "b".repeat(64)
    };
    fs::create_dir_all(authoring.join("references")).expect("reference directory");
    fs::write(authoring.join("references/ref.fa"), b">chr1\nAA\n").expect("reference");
    fs::write(authoring.join("general.toml"), general()).expect("general");
    fs::create_dir_all(&scenarios).expect("scenario directory");
    fs::write(scenarios.join("scenario.toml"), scenario_source).expect("scenario");
    authoring.join("general.toml")
}

#[test]
fn dispatcher_classifies_discriminators_after_4096_bytes_and_preserves_non_f5() {
    let root = temporary("late-discriminators");
    let profile_hash = write_profile(&root.join("authoring/scenarios/profiles"));
    fs::create_dir_all(root.join("authoring/references")).expect("reference directory");
    fs::write(root.join("authoring/references/ref.fa"), b">chr1\nAA\n").expect("reference");
    let late_general = format!("{}\n{}", "#".repeat(4_097), general());
    let general_path = root.join("authoring/general.toml");
    fs::write(&general_path, late_general).expect("late general");
    fs::write(
        root.join("authoring/scenarios/scenario.toml"),
        scenario(&profile_hash, "999999"),
    )
    .expect("scenario");

    let dispatched = load_and_validate_cli_config(&general_path).expect("F5 after 4096 bytes");
    assert!(matches!(dispatched, CliValidatedConfig::F5(_)));

    let legacy = root.join("legacy.toml");
    fs::write(
        &legacy,
        include_str!("fixtures/config/general-minimal.toml"),
    )
    .expect("legacy general");
    assert!(matches!(
        load_and_validate_cli_config(&legacy).expect("legacy remains legacy"),
        CliValidatedConfig::Legacy(_)
    ));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn bootstrap_accepts_exact_general_and_scenario_limits_and_rejects_byte_65537() {
    for (role, pad_general) in [("general", true), ("scenario", false)] {
        let root = temporary(role);
        let profile_hash = write_profile(&root.join("authoring/scenarios/profiles"));
        fs::create_dir_all(root.join("authoring/references")).expect("reference directory");
        fs::write(root.join("authoring/references/ref.fa"), b">chr1\nAA\n").expect("reference");
        let general_source = if pad_general {
            pad_to_toml_limit(general())
        } else {
            general()
        };
        let scenario_source = if pad_general {
            scenario(&profile_hash, "999999")
        } else {
            pad_to_toml_limit(scenario(&profile_hash, "999999"))
        };
        let general_path = root.join("authoring/general.toml");
        fs::write(&general_path, &general_source).expect("general");
        fs::write(
            root.join("authoring/scenarios/scenario.toml"),
            &scenario_source,
        )
        .expect("scenario");

        let exact = load_and_validate_cli_config(&general_path).expect("exact limit");
        assert!(matches!(exact, CliValidatedConfig::F5(_)));

        let oversized = if pad_general {
            &general_path
        } else {
            &root.join("authoring/scenarios/scenario.toml")
        };
        fs::OpenOptions::new()
            .append(true)
            .open(oversized)
            .expect("oversized document")
            .write_all(b"#")
            .expect("byte 65537");
        let error = load_and_validate_cli_config(&general_path).expect_err("65537 must reject");
        assert_eq!(error.code(), "GPS-INPUT-001");
        assert!(!error.context().contains_key("scenario_id"));
        assert!(error.render_human().contains("byte 65537"));
        fs::remove_dir_all(root).expect("cleanup");
    }
}

#[test]
fn dispatcher_opens_each_bootstrap_toml_once_and_records_identity_from_written_bytes() {
    let root = temporary("physical-identity-and-opens");
    // El dispatcher público completa ahora el perfil; el observador sigue
    // cubriendo exclusivamente la frontera TOML y debe registrar una apertura
    // por cada documento de autoría, no las aperturas del perfil.
    let profile_hash = write_profile(&root.join("authoring/scenarios/profiles"));
    let general_bytes = pad_to_toml_limit(general()).into_bytes();
    let scenario_bytes = pad_to_toml_limit(scenario(&profile_hash, "999999")).into_bytes();
    let general_path = root.join("authoring/general.toml");
    let scenario_path = root.join("authoring/scenarios/scenario.toml");
    fs::create_dir_all(scenario_path.parent().expect("scenario parent")).expect("authoring");
    fs::create_dir_all(root.join("authoring/references")).expect("reference directory");
    fs::write(root.join("authoring/references/ref.fa"), b">chr1\nAA\n").expect("reference");
    fs::write(&general_path, &general_bytes).expect("general");
    fs::write(&scenario_path, &scenario_bytes).expect("scenario");

    let mut opened = Vec::new();
    let dispatched =
        load_and_validate_cli_config_with_bootstrap_open_observer(&general_path, &mut |path| {
            opened.push(path.to_path_buf())
        })
        .expect("bootstrap");
    let CliValidatedConfig::F5(bootstrap) = dispatched else {
        panic!("general v3 must dispatch to F5");
    };

    let documents = bootstrap.bootstrap_documents();
    assert_eq!(documents.general().size_bytes(), general_bytes.len() as u64);
    assert_eq!(documents.general().sha256(), sha256(&general_bytes));
    assert_eq!(
        documents.scenario().size_bytes(),
        scenario_bytes.len() as u64
    );
    assert_eq!(documents.scenario().sha256(), sha256(&scenario_bytes));
    assert_eq!(
        opened,
        vec![
            fs::canonicalize(&general_path).expect("canonical general"),
            fs::canonicalize(&scenario_path).expect("canonical scenario"),
        ],
        "el observador cuenta solo los dos TOML; la admisión .4 abre después perfil y referencia sin reabrirlos"
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn dispatcher_rejects_scenario_limit_plus_one_and_very_large_after_only_two_opens() {
    let root = temporary("scenario-limit-plus-one");
    let profile_hash = "b".repeat(64);
    let general_path = root.join("authoring/general.toml");
    let scenario_path = root.join("authoring/scenarios/scenario.toml");
    fs::create_dir_all(scenario_path.parent().expect("scenario parent")).expect("authoring");
    fs::write(&general_path, general()).expect("general");

    for (label, suffix) in [
        ("65537", "#".to_owned()),
        ("very-large", "#".repeat(TOML_LIMIT * 16)),
    ] {
        let mut scenario_bytes = pad_to_toml_limit(scenario(&profile_hash, "999999")).into_bytes();
        scenario_bytes.extend_from_slice(suffix.as_bytes());
        fs::write(&scenario_path, scenario_bytes).expect("oversized scenario");
        let mut opened = Vec::new();
        let error =
            load_and_validate_cli_config_with_bootstrap_open_observer(&general_path, &mut |path| {
                opened.push(path.to_path_buf())
            })
            .expect_err(label);
        assert_eq!(error.code(), "GPS-INPUT-001", "{label}");
        assert!(error.render_human().contains("byte 65537"), "{label}");
        assert_eq!(
            opened,
            vec![
                fs::canonicalize(&general_path).expect("canonical general"),
                fs::canonicalize(&scenario_path).expect("canonical scenario"),
            ],
            "{label}: no hay una segunda apertura ni recursos posteriores",
        );
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn every_public_general_loader_applies_the_v3_bound_before_classification_or_materialization() {
    let root = temporary("all-public-loaders");
    let profile_hash = write_profile(&root.join("authoring/scenarios/profiles"));
    fs::create_dir_all(root.join("authoring/references")).expect("reference directory");
    fs::write(root.join("authoring/references/ref.fa"), b">chr1\nAA\n").expect("reference");
    let general_path = root.join("authoring/general.toml");
    let exact_v3 = pad_to_toml_limit(general());
    fs::write(&general_path, &exact_v3).expect("exact v3 general");
    fs::write(
        root.join("authoring/scenarios/scenario.toml"),
        scenario(&profile_hash, "999999"),
    )
    .expect("scenario");

    // El cargador histórico no adopta v3, pero solo lo decide después de la
    // misma lectura física acotada. Los demás conservan sus contratos v3.
    assert_eq!(
        load_and_validate_config(&general_path)
            .expect_err("the legacy API still rejects v3")
            .code(),
        "GPS-CONF-005"
    );
    assert_eq!(
        load_and_validate_config_with_entropy(&general_path, &mut PanicEntropy)
            .expect_err("v3 must reject before entropy")
            .code(),
        "GPS-CONF-005"
    );
    assert!(matches!(
        load_and_validate_cli_config(&general_path).expect("CLI exact v3"),
        CliValidatedConfig::F5(_)
    ));
    assert_eq!(
        load_and_validate_f5_config(&general_path)
            .expect("F5 exact v3")
            .bootstrap_documents()
            .general()
            .size_bytes(),
        TOML_LIMIT as u64
    );
    assert_eq!(
        load_and_validate_f5_config_with_entropy(&general_path, &mut PanicEntropy)
            .expect("F5 exact v3 with injected entropy")
            .bootstrap_documents()
            .general()
            .size_bytes(),
        TOML_LIMIT as u64
    );

    for (label, source) in [
        ("byte-65537", format!("{exact_v3}[")),
        (
            "very-large",
            format!("{exact_v3}{}", "#".repeat(TOML_LIMIT * 16)),
        ),
    ] {
        fs::write(&general_path, source).expect("oversized v3 general");
        let errors = [
            load_and_validate_config(&general_path).map(|_| ()),
            load_and_validate_config_with_entropy(&general_path, &mut PanicEntropy).map(|_| ()),
            load_and_validate_cli_config(&general_path).map(|_| ()),
            load_and_validate_f5_config(&general_path).map(|_| ()),
            load_and_validate_f5_config_with_entropy(&general_path, &mut PanicEntropy).map(|_| ()),
        ];
        for result in errors {
            let error = result.expect_err(label);
            assert_eq!(error.code(), "GPS-INPUT-001", "{label}");
            assert!(!error.context().contains_key("scenario_id"), "{label}");
            assert!(error.render_human().contains("byte 65537"), "{label}");
            assert!(!error.render_human().contains("F5"), "{label}");
        }
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn oversized_documents_are_neutral_before_v2_v3_or_invalid_classification() {
    let root = temporary("neutral-limit-diagnostic");
    for (label, prefix) in [
        (
            "v2",
            include_str!("fixtures/config/general-minimal.toml").to_owned(),
        ),
        ("v3", general()),
        (
            "unclassifiable",
            "# no configuration discriminators\n".to_owned(),
        ),
    ] {
        let path = root.join(format!("{label}.toml"));
        fs::create_dir_all(&root).expect("root");
        fs::write(&path, format!("{prefix}{}", "#".repeat(TOML_LIMIT * 16)))
            .expect("large document");
        let error = load_and_validate_cli_config(&path).expect_err("limit before dispatch");
        assert_eq!(error.code(), "GPS-INPUT-001", "{label}");
        assert!(!error.context().contains_key("scenario_id"), "{label}");
        assert!(error.render_human().contains("byte 65537"), "{label}");
        assert!(!error.render_human().contains("F5"), "{label}");
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn dispatcher_preserves_relative_v2_authoring_path_while_opening_the_normalized_path() {
    let requested = Path::new("tests/fixtures/config/general-minimal.toml");
    let CliValidatedConfig::Legacy(effective) =
        load_and_validate_cli_config(requested).expect("relative v2 through dispatcher")
    else {
        panic!("general v2 must remain legacy");
    };
    assert_eq!(
        effective.source_path(),
        requested.to_str().expect("UTF-8 path")
    );
    assert_eq!(
        effective.configuration_input().original_path(),
        requested.to_str().expect("UTF-8 path")
    );
}

#[test]
fn bootstrap_rejects_invalid_final_toml_before_profile_reference_mask_destination_or_rng() {
    let root = temporary("invalid-final");
    let source = format!("{}invalid = [", scenario(&"b".repeat(64), "999999"));
    let general = write_fixture(&root, false, source);
    let error = load_and_validate_f5_config_with_entropy(&general, &mut PanicEntropy)
        .expect_err("invalid final TOML must stop bootstrap");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert!(
        !root
            .join("authoring/scenarios/profiles/profile.metadata.json")
            .exists(),
        "the test spy leaves the profile path absent: its open would replace the TOML diagnostic"
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn fixed_width_bootstrap_memory_limit_is_inclusive_and_minus_one_stops_before_profile_and_rng() {
    let root = temporary("fixed-point");
    let profile_hash = write_profile(&root.join("authoring/scenarios/profiles"));
    fs::create_dir_all(root.join("authoring/references")).expect("reference directory");
    fs::write(root.join("authoring/references/ref.fa"), b">chr1\nAA\n").expect("reference");
    let general_source = general();
    let scenario_source = pad_to_toml_limit(scenario(&profile_hash, "999999"));
    let exact_limit = 2_u64
        * u64::try_from(general_source.len() + scenario_source.len()).expect("bootstrap size");
    assert_eq!(exact_limit.to_string().len(), 6, "fixed decimal width");
    let exact_scenario = scenario_source.replace(
        "max_working_memory_bytes = 999999",
        &format!("max_working_memory_bytes = {exact_limit}"),
    );
    assert_eq!(exact_scenario.len(), scenario_source.len());
    let general_path = root.join("authoring/general.toml");
    fs::write(&general_path, &general_source).expect("general");
    fs::write(
        root.join("authoring/scenarios/scenario.toml"),
        &exact_scenario,
    )
    .expect("scenario");
    load_and_validate_cli_config(&general_path).expect("exact bootstrap bound is inclusive");

    let minus_one = exact_scenario.replace(
        &format!("max_working_memory_bytes = {exact_limit}"),
        &format!("max_working_memory_bytes = {}", exact_limit - 1),
    );
    assert_eq!(minus_one.len(), exact_scenario.len(), "fixed-width −1");
    fs::write(root.join("authoring/scenarios/scenario.toml"), minus_one).expect("minus one");
    fs::remove_file(root.join("authoring/scenarios/profiles/profile.metadata.json"))
        .expect("remove profile spy");
    let error = load_and_validate_f5_config_with_entropy(&general_path, &mut PanicEntropy)
        .expect_err("minus one must stop before profile/RNG");
    assert_eq!(error.code(), "GPS-RES-101");
    assert!(error.render_human().contains("phase=bootstrap"));
    fs::remove_dir_all(root).expect("cleanup");
}
