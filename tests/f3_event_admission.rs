use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use serde_json::json;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{
    EffectiveTechnicalEventTarget, ILLUMINA_SE_MINIMAL_ACGT_PHRED, ILLUMINA_SE_MINIMAL_N_PHRED,
    ILLUMINA_SE_MINIMAL_READ_LENGTH, ILLUMINA_SE_MINIMAL_V1, build_effective_config,
};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::results::validate_effective_config_against_manifest_bytes;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn config(event: &str) -> String {
    format!(
        "config_kind = \"general\"\nschema_version = \"1\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1000\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\", \"sample-B\"]\n\n{event}\n"
    )
}

fn declared(event: &str) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let admitted = admit_toml_header(
        Path::new("f3-event.toml"),
        config(event).into_bytes(),
        ExpectedConfigKind::General,
    )
    .expect("TOML header");
    parse_declared_config(admitted).expect("declared configuration")
}

#[test]
fn f3_minimal_technical_event_resolves_fixed_profile_and_exactly_one_target() {
    let effective = build_effective_config(&declared(
        "[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [\"sample-A\", \"sample-B\"]\ntarget_coverage_wgs = 2.5",
    ))
    .expect("effective configuration");
    let event = effective.technical_event().expect("technical event");

    assert_eq!(event.event_id().as_str(), "reads-main");
    assert_eq!(
        event
            .individual_ids()
            .iter()
            .map(|id| id.as_str())
            .collect::<Vec<_>>(),
        ["sample-A", "sample-B"]
    );
    assert!(matches!(
        event.target(),
        EffectiveTechnicalEventTarget::CoverageWgs(value) if value.to_bits() == 2.5_f64.to_bits()
    ));
    assert_eq!(event.profile(), ILLUMINA_SE_MINIMAL_V1);
    assert_eq!(event.read_length(), ILLUMINA_SE_MINIMAL_READ_LENGTH);
    assert_eq!(event.acgt_phred(), ILLUMINA_SE_MINIMAL_ACGT_PHRED);
    assert_eq!(event.n_phred(), ILLUMINA_SE_MINIMAL_N_PHRED);

    let count = build_effective_config(&declared(
        "[technical_event]\nevent_id = \"reads-count\"\nindividual_ids = [\"sample-A\"]\nread_count = 0",
    ))
    .expect("count target");
    assert!(matches!(
        count.technical_event().expect("event").target(),
        EffectiveTechnicalEventTarget::ReadCount(0)
    ));
}

#[test]
fn f3_technical_event_rejects_invalid_target_selection_and_individuals() {
    for (name, event) in [
        (
            "both-targets",
            "[technical_event]\nevent_id = \"reads\"\nindividual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 1\nread_count = 2",
        ),
        (
            "no-target",
            "[technical_event]\nevent_id = \"reads\"\nindividual_ids = [\"sample-A\"]",
        ),
        (
            "empty-individuals",
            "[technical_event]\nevent_id = \"reads\"\nindividual_ids = []\nread_count = 1",
        ),
    ] {
        let admitted = admit_toml_header(
            Path::new(name),
            config(event).into_bytes(),
            ExpectedConfigKind::General,
        )
        .expect("header");
        assert_eq!(
            parse_declared_config(admitted).expect_err(name).code(),
            "GPS-CONF-008"
        );
    }

    for (name, event) in [
        (
            "unknown-individual",
            "[technical_event]\nevent_id = \"reads\"\nindividual_ids = [\"absent\"]\nread_count = 1",
        ),
        (
            "duplicate-individual",
            "[technical_event]\nevent_id = \"reads\"\nindividual_ids = [\"sample-A\", \"sample-A\"]\nread_count = 1",
        ),
        (
            "negative-coverage",
            "[technical_event]\nevent_id = \"reads\"\nindividual_ids = [\"sample-A\"]\ntarget_coverage_wgs = -1",
        ),
    ] {
        assert_eq!(
            build_effective_config(&declared(event))
                .expect_err(name)
                .code(),
            "GPS-CONF-008"
        );
    }
}

#[test]
fn f3_technical_event_rejects_f4_fields_and_shared_event_id_collision() {
    for field in [
        "profile = \"future\"",
        "paired_end = true",
        "contamination = 0.1",
        "mapping_reference = \"reference-main\"",
    ] {
        let event = format!(
            "[technical_event]\nevent_id = \"reads\"\nindividual_ids = [\"sample-A\"]\nread_count = 1\n{field}"
        );
        let admitted = admit_toml_header(
            Path::new("f4-field.toml"),
            config(&event).into_bytes(),
            ExpectedConfigKind::General,
        )
        .expect("header");
        assert_eq!(
            parse_declared_config(admitted).expect_err(field).code(),
            "GPS-CONF-006"
        );
    }

    let collision = "[[direct_events]]\nevent_id = \"same-event\"\nindividual_id = \"sample-A\"\nphase = 0\nkind = \"snp\"\ncontig = \"chr1\"\nposition_1based = 1\nref = \"A\"\nalt = \"C\"\n\n[technical_event]\nevent_id = \"same-event\"\nindividual_ids = [\"sample-A\"]\nread_count = 1";
    assert_eq!(
        build_effective_config(&declared(collision))
            .expect_err("event collision")
            .code(),
        "GPS-CONF-008"
    );
}

#[test]
fn f3_technical_event_is_written_to_and_revalidated_from_effective_v3() {
    let root = temporary("published-effective");
    fs::create_dir_all(&root).expect("test directory");
    let config_path = root.join("general.toml");
    fs::write(
        &config_path,
        config("[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [\"sample-A\"]\nread_count = 7")
            .replacen(
                "schema_version = \"1\"\n\n",
                "schema_version = \"1\"\n\n[population_snvs]\nmode = \"generate\"\n[population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 1\nspatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n[[population_snvs.generate.populations]]\npopulation_id = \"P1\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n\n",
                1,
            ),
    )
    .expect("configuration");
    let output = root.join("output");
    let run = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run"])
        .arg(&config_path)
        .arg("--output")
        .arg(&output)
        .output()
        .expect("published run");
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );

    let effective_path = output.join("effective-config.json");
    let manifest_path = output.join("manifest.json");
    let effective_bytes = fs::read(&effective_path).expect("effective configuration");
    let manifest_bytes = fs::read(&manifest_path).expect("manifest");
    validate_effective_config_against_manifest_bytes(
        &effective_path,
        &effective_bytes,
        &manifest_path,
        &manifest_bytes,
    )
    .expect("valid published v3 documents");
    let effective: serde_json::Value =
        serde_json::from_slice(&effective_bytes).expect("effective JSON");
    assert_eq!(effective["schema_version"], "3");
    assert_eq!(
        effective["technical_event"],
        json!({
            "event_id": "reads-main",
            "individual_ids": ["sample-A"],
            "profile": "illumina-se-minimal-v1",
            "read_length": 150,
            "quality": {"acgt_phred": 40, "n_phred": 0},
            "target": {"read_count": 7},
        })
    );

    let mut unknown_individual = effective.clone();
    unknown_individual["technical_event"]["individual_ids"][0] = json!("absent");
    assert!(
        validate_effective_config_against_manifest_bytes(
            &effective_path,
            &serde_json::to_vec(&unknown_individual).expect("effective JSON"),
            &manifest_path,
            &manifest_bytes,
        )
        .is_err()
    );

    let mut colliding_event = effective.clone();
    colliding_event["direct_events"] = json!([{
        "event_id": "reads-main",
        "individual_id": "sample-A",
        "phase": 0,
        "kind": "snp",
        "contig": "chr1",
        "position_1based": 1,
        "ref": "A",
        "alt": "C",
    }]);
    assert!(
        validate_effective_config_against_manifest_bytes(
            &effective_path,
            &serde_json::to_vec(&colliding_event).expect("effective JSON"),
            &manifest_path,
            &manifest_bytes,
        )
        .is_err()
    );
    fs::remove_dir_all(root).expect("cleanup");
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f3-event-admission-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
