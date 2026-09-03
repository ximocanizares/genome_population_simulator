//! Regresiones del control PE neutro VE-D02.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::error_quality::SequencingErrorsMode;
use genome_population_simulator::execution::execute_validated_config;
use serde_json::Value;

static TEMPORARY_COUNTER: AtomicU64 = AtomicU64::new(0);

#[test]
fn default_mode_remains_enabled() {
    let effective = build_effective_config(&declared(None)).expect("effective default");
    let explicit = build_effective_config(&declared(Some("enabled"))).expect("effective enabled");
    assert_eq!(
        effective
            .technical_event()
            .expect("technical event")
            .sequencing_errors(),
        SequencingErrorsMode::Enabled
    );
    assert_eq!(effective.technical_event(), explicit.technical_event());
}

#[test]
fn disabled_mode_preserves_truth_and_marks_effective_and_event_manifest() {
    let effective =
        build_effective_config(&declared(Some("disabled"))).expect("effective disabled");
    assert_eq!(
        effective
            .technical_event()
            .expect("technical event")
            .sequencing_errors(),
        SequencingErrorsMode::Disabled
    );
    let destination = temporary("disabled");
    execute_validated_config(&effective, &destination).expect("disabled PE publication");

    let effective_document = read_json(&destination.join("effective-config.json"));
    let effective_errors = &effective_document["technical_event"]["sequencing_errors"];
    assert_eq!(effective_errors["requested_mode"], "disabled");
    assert_eq!(effective_errors["effective_mode"], "disabled");
    assert_eq!(
        effective_errors["qualities_represent_realized_error_probability"],
        false
    );

    let event = read_json(&destination.join("events/neutral-main/manifest.json"));
    let errors = &event["error_quality"];
    assert_eq!(errors["requested_mode"], "disabled");
    assert_eq!(errors["effective_mode"], "disabled");
    assert_eq!(
        errors["qualities_represent_realized_error_probability"],
        false
    );
    for mate in ["r1", "r2"] {
        let mate_errors = errors[mate].as_object().expect("mate error fields");
        assert_eq!(mate_errors.len(), 4);
        assert_eq!(errors[mate]["quality_shift"], 0);
        assert_eq!(errors[mate]["error_scale"], 1.0);
        assert!(
            errors[mate]["opportunities"]
                .as_u64()
                .expect("opportunities")
                > 0
        );
        assert_eq!(errors[mate]["errors"], 0);
    }
    fs::remove_dir_all(destination).expect("temporary output cleanup");
}

#[test]
fn invalid_or_ambiguous_sequencing_error_modes_are_rejected() {
    assert!(build_effective_config(&declared(Some("off"))).is_err());

    let duplicate = source(Some("disabled")).replace(
        "sequencing_errors = \"disabled\"",
        "sequencing_errors = \"disabled\"\nsequencing_errors = \"enabled\"",
    );
    assert!(
        admit_toml_header(
            Path::new("duplicate-sequencing-errors.toml"),
            duplicate.into_bytes(),
            ExpectedConfigKind::General,
        )
        .is_err()
    );
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("JSON")).expect("valid JSON")
}

fn declared(
    mode: Option<&str>,
) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let text = source(mode);
    parse_declared_config(
        admit_toml_header(
            Path::new("sequencing-errors.toml"),
            text.into_bytes(),
            ExpectedConfigKind::General,
        )
        .expect("header"),
    )
    .expect("declared")
}

fn source(mode: Option<&str>) -> String {
    let sequencing_errors = mode.map_or_else(String::new, |mode| {
        format!("sequencing_errors = \"{mode}\"\n")
    });
    format!(
        "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 1\nC = 0\nG = 0\nT = 0\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 500\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 1\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\"]\n\n[population_snvs]\nmode = \"generate\"\n\n[population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 1\nspatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n\n[[population_snvs.generate.populations]]\npopulation_id = \"P1\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n\n[technical_event]\nevent_id = \"neutral-main\"\nindividual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 1\nread_length_bases = 1\n{sequencing_errors}\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n"
    )
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-sequencing-errors-{label}-{}-{}",
        std::process::id(),
        TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed)
    ))
}
