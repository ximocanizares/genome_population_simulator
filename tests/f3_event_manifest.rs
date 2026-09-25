use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::event_manifest::{
    EventManifest, EventManifestArtifact, EventManifestCoverage, EventManifestDemultiplexing,
    EventManifestFinalFastqSummary, EventManifestReference, event_manifest_bytes,
    event_manifest_document, event_manifest_relative_path, validate_event_manifest_bytes,
};
use genome_population_simulator::fastq_validation::F4CoverageReport;
use genome_population_simulator::identity::ReferenceId;
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::results::{ArtifactRecord, UnitState};
use genome_population_simulator::technical_identity::TechnicalMate;
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use serde_json::json;

static NONCE: AtomicU64 = AtomicU64::new(0);

const ROOT_SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
const DIGEST_A: &str = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
const DIGEST_B: &str = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";

#[test]
fn complete_technical_event_manifest_is_written_and_accepted_by_both_validators() {
    let root = temporary("complete");
    let (manifest, event_id) = complete_manifest();
    let bytes = event_manifest_bytes(&manifest).expect("event JSON");
    let document: serde_json::Value = serde_json::from_slice(&bytes).expect("event document");
    assert_eq!(document["technology"], "illumina-short-read-pe-v1");
    validate_event_manifest_bytes(Path::new("events/reads-main/manifest.json"), &bytes)
        .expect("Rust validator");

    let manifest_path = root.join("events").join(event_id).join("manifest.json");
    fs::create_dir_all(manifest_path.parent().expect("manifest parent")).expect("manifest parent");
    fs::write(&manifest_path, &bytes).expect("test manifest");
    assert_eq!(
        event_manifest_relative_path(manifest.event.event_id()).to_string_lossy(),
        "events/reads-main/manifest.json"
    );
    assert_external(&manifest_path, true);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn invalid_schema_and_invalid_library_link_are_rejected() {
    let root = temporary("invalid");
    fs::create_dir_all(&root).expect("test directory");
    let (manifest, _) = complete_manifest();
    let mut malformed = event_manifest_document(&manifest).expect("event document");
    malformed["outputs"][0]["artifact"]["checksum"]["digest"] = json!("A".repeat(64));
    let malformed_bytes = serde_json::to_vec(&malformed).expect("malformed JSON");
    assert!(validate_event_manifest_bytes(Path::new("manifest.json"), &malformed_bytes).is_err());
    let malformed_path = root.join("malformed.json");
    fs::write(&malformed_path, malformed_bytes).expect("write malformed");
    assert_external(&malformed_path, false);

    let mut wrong_library = event_manifest_document(&manifest).expect("event document");
    wrong_library["outputs"][0]["library_id"] = json!("LIB_not_selected");
    let wrong_library_bytes = serde_json::to_vec(&wrong_library).expect("wrong library JSON");
    assert!(
        validate_event_manifest_bytes(Path::new("manifest.json"), &wrong_library_bytes).is_err()
    );

    let mut single_end = event_manifest_document(&manifest).expect("event document");
    single_end["outputs"] = json!([single_end["outputs"][0].clone()]);
    assert_rejected_by_both(&root, "single-end-output", &single_end);

    let mut single_end_profile = event_manifest_document(&manifest).expect("event document");
    single_end_profile["profile"]["id"] = json!("illumina-se-minimal-v1");
    assert_rejected_by_both(&root, "single-end-profile", &single_end_profile);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn chimera_summary_is_closed_and_rejects_invalid_rate_counter_or_provenance() {
    let root = temporary("chimera-summary");
    fs::create_dir_all(&root).expect("test directory");
    let (manifest, _) = complete_manifest();
    let document = event_manifest_document(&manifest).expect("event document");

    let mut invalid_rate = document.clone();
    invalid_rate["artifacts"]["chimera"]["rate"] = json!(1.1);
    assert_rejected_by_both(&root, "chimera-invalid-rate", &invalid_rate);

    let mut invalid_counter = document.clone();
    invalid_counter["artifacts"]["chimera"]["realized_pairs"] = json!("four");
    assert_rejected_by_both(&root, "chimera-invalid-counter", &invalid_counter);

    let mut zero_rate_with_realized_pairs = document.clone();
    zero_rate_with_realized_pairs["artifacts"]["chimera"]["realized_pairs"] = json!(1);
    assert_rejected_by_both(
        &root,
        "chimera-zero-rate-with-realized-pairs",
        &zero_rate_with_realized_pairs,
    );

    let mut leaked_provenance = document;
    leaked_provenance["artifacts"]["chimera"]["pair_ids"] = json!(["TP_hidden"]);
    assert_rejected_by_both(&root, "chimera-provenance", &leaked_provenance);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn error_quality_adjustments_are_fixed_for_each_mate_in_rust_and_schema() {
    let root = temporary("fixed-error-quality-adjustments");
    fs::create_dir_all(&root).expect("test directory");
    let (manifest, _) = complete_manifest();
    let document = event_manifest_document(&manifest).expect("event document");

    for (label, mate, field, value) in [
        ("r1-quality-shift", "r1", "quality_shift", json!(1)),
        ("r1-error-scale", "r1", "error_scale", json!(2.0)),
        ("r2-quality-shift", "r2", "quality_shift", json!(-1)),
        ("r2-error-scale", "r2", "error_scale", json!(0.5)),
    ] {
        let mut mutated = document.clone();
        mutated["error_quality"][mate][field] = value;
        assert_rejected_by_both(&root, label, &mutated);
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn each_selected_library_has_one_canonical_fastq_result() {
    let (manifest, _) = complete_manifest();
    let document = event_manifest_document(&manifest).expect("event document");

    let mut missing = document.clone();
    missing["outputs"] = json!([]);
    assert!(validate_document(&missing).is_err());

    let mut duplicate = document.clone();
    let mut repeated = duplicate["outputs"][0].clone();
    repeated["id"] = json!("fastq-duplicate");
    repeated["artifact"]["id"] = json!("fastq-duplicate");
    repeated["artifact"]["path"] = json!("outputs/fastq/reads-main/sample-A_R1.fastq.gz");
    duplicate["outputs"]
        .as_array_mut()
        .expect("outputs")
        .push(repeated);
    assert!(validate_document(&duplicate).is_err());

    let mut two_libraries = document.clone();
    two_libraries["selected_individuals"]
        .as_array_mut()
        .expect("individuals")
        .push(json!({
            "individual_id": "sample-B", "library_id": "LIB_sample-B", "read_group_id": "RG_sample-B",
            "i7": "GGGGGGGG", "i5": "CCCCCCCC",
        }));
    let mut second = two_libraries["outputs"][0].clone();
    second["id"] = json!("fastq-sample-B");
    second["individual_id"] = json!("sample-B");
    second["library_id"] = json!("LIB_sample-B");
    second["artifact"]["id"] = json!("fastq-sample-B");
    second["artifact"]["path"] = json!("outputs/fastq/reads-main/sample-B_R1.fastq.gz");
    two_libraries["outputs"]
        .as_array_mut()
        .expect("outputs")
        .push(second);
    let mut second_r2 = two_libraries["outputs"][1].clone();
    second_r2["id"] = json!("fastq-sample-B-r2");
    second_r2["individual_id"] = json!("sample-B");
    second_r2["library_id"] = json!("LIB_sample-B");
    second_r2["mate"] = json!(2);
    second_r2["artifact"]["id"] = json!("fastq-sample-B-r2");
    second_r2["artifact"]["path"] = json!("outputs/fastq/reads-main/sample-B_R2.fastq.gz");
    two_libraries["outputs"]
        .as_array_mut()
        .expect("outputs")
        .push(second_r2);
    let mut second_coverage = two_libraries["coverage"][0].clone();
    second_coverage["individual_id"] = json!("sample-B");
    second_coverage["library_id"] = json!("LIB_sample-B");
    two_libraries["coverage"]
        .as_array_mut()
        .expect("coverage")
        .push(second_coverage);
    let mut second_summary = two_libraries["final_fastq_summary"][0].clone();
    second_summary["individual_id"] = json!("sample-B");
    second_summary["library_id"] = json!("LIB_sample-B");
    two_libraries["final_fastq_summary"]
        .as_array_mut()
        .expect("summary")
        .push(second_summary);
    assert!(validate_document(&two_libraries).is_ok());
    two_libraries["outputs"]
        .as_array_mut()
        .expect("outputs")
        .reverse();
    assert!(validate_document(&two_libraries).is_err());
}

#[test]
fn complete_fastq_path_has_the_canonical_shape_and_correlates_to_its_event_and_individual() {
    let root = temporary("fastq-path");
    fs::create_dir_all(&root).expect("test directory");
    let (manifest, _) = complete_manifest();
    let document = event_manifest_document(&manifest).expect("event document");

    let mut malformed_shape = document.clone();
    malformed_shape["outputs"][0]["artifact"]["path"] =
        json!("outputs/fastq/reads-main/nested/sample-A_R1.fastq.gz");
    assert_rejected_by_both(&root, "malformed-fastq-path", &malformed_shape);

    let mut wrong_event = document;
    wrong_event["outputs"][0]["artifact"]["path"] =
        json!("outputs/fastq/other-event/sample-A_R1.fastq.gz");
    assert!(validate_document(&wrong_event).is_err());
    let wrong_event_path = root.join("wrong-event-fastq-path.json");
    fs::write(
        &wrong_event_path,
        serde_json::to_vec(&wrong_event).expect("JSON"),
    )
    .expect("test document");
    assert_external(&wrong_event_path, true);

    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn state_diagnostics_and_artifact_coherence_is_checked_externally() {
    let root = temporary("states");
    fs::create_dir_all(&root).expect("test directory");
    let (manifest, _) = complete_manifest();
    let complete = event_manifest_document(&manifest).expect("event document");

    let mut complete_without_coverage = complete.clone();
    complete_without_coverage["coverage"] = json!([]);
    assert_rejected_by_both(
        &root,
        "complete-without-coverage",
        &complete_without_coverage,
    );

    let mut complete_with_diagnostic = complete.clone();
    complete_with_diagnostic["diagnostics"] = json!([format_diagnostic()]);
    assert_rejected_by_both(&root, "complete-with-diagnostic", &complete_with_diagnostic);

    let mut failed_without_diagnostic = complete.clone();
    failed_without_diagnostic["state"] = json!("failed");
    assert_rejected_by_both(
        &root,
        "failed-without-diagnostic",
        &failed_without_diagnostic,
    );

    let failed = failed_document(&complete);
    assert_accepted_by_both(&root, "failed", &failed);

    let mut failed_with_coverage = failed.clone();
    failed_with_coverage["coverage"] = complete["coverage"].clone();
    assert_rejected_by_both(&root, "failed-with-coverage", &failed_with_coverage);

    let mut ordered_diagnostics = failed.clone();
    ordered_diagnostics["diagnostics"]
        .as_array_mut()
        .expect("diagnostics")
        .push(json!({
            "code": "GPS-IO-001",
            "message": "Una operación de lectura o escritura falló.",
        }));
    assert!(validate_document(&ordered_diagnostics).is_ok());
    ordered_diagnostics["diagnostics"]
        .as_array_mut()
        .expect("diagnostics")
        .reverse();
    assert!(validate_document(&ordered_diagnostics).is_err());

    let mut invalid_code = failed;
    invalid_code["diagnostics"][0]["code"] = json!("GPS-BOGUS-001");
    invalid_code["outputs"][0]["diagnostic_code"] = json!("GPS-BOGUS-001");
    assert_rejected_by_both(&root, "invalid-diagnostic-code", &invalid_code);

    let mut f5_only_diagnostic = failed_document(&complete);
    f5_only_diagnostic["diagnostics"][0] = json!({
        "code": "GPS-CONF-101",
        "message": "El escenario F5, su objetivo, F/P o sus límites no satisfacen el contrato de admisión.",
    });
    f5_only_diagnostic["outputs"][0]["diagnostic_code"] = json!("GPS-CONF-101");
    assert_rejected_by_both(&root, "f4-event-f5-diagnostic", &f5_only_diagnostic);

    fs::remove_dir_all(root).expect("cleanup");
}

fn complete_manifest() -> (EventManifest<'static>, &'static str) {
    let effective = Box::leak(Box::new(
        build_effective_config(&declared()).expect("effective"),
    ));
    let event = effective.technical_event().expect("technical event");
    let reference = Box::leak(Box::new(
        ReferenceId::parse("reference-main", "test").expect("reference ID"),
    ));
    let identities = Box::leak(Box::new(
        derive_technical_identity_plan(
            RootSeed::parse("root", ROOT_SEED).expect("seed"),
            reference,
            event.event_id().clone(),
            event.individual_ids(),
        )
        .expect("technical identities"),
    ));
    let fasta = Box::leak(Box::new(record(
        "reference-main-fasta",
        "reference_fasta",
        "resources/references/reference-main.fa",
        DIGEST_A,
    )));
    let individual = event.individual_ids().first().expect("individual");
    let library = identities.library(individual).expect("library");
    let artifact_id =
        Box::leak(format!("fastq-{}", library.library_id().as_str()).into_boxed_str());
    let fastq_r1 = Box::leak(Box::new(record(
        artifact_id,
        "fastq",
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz",
        DIGEST_B,
    )));
    let artifact_r2 = Box::leak(format!("{artifact_id}-r2").into_boxed_str());
    let fastq_r2 = Box::leak(Box::new(record(
        artifact_r2,
        "fastq",
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz",
        DIGEST_A,
    )));
    let artifacts = Box::leak(Box::new([
        EventManifestArtifact {
            id: artifact_id,
            kind: "fastq",
            state: UnitState::Complete,
            individual_id: individual,
            library_id: library.library_id(),
            mate: TechnicalMate::R1,
            artifact: Some(fastq_r1),
            diagnostic: None,
        },
        EventManifestArtifact {
            id: artifact_r2,
            kind: "fastq",
            state: UnitState::Complete,
            individual_id: individual,
            library_id: library.library_id(),
            mate: TechnicalMate::R2,
            artifact: Some(fastq_r2),
            diagnostic: None,
        },
    ]));
    let coverage = Box::leak(Box::new([EventManifestCoverage {
        individual_id: individual,
        library_id: library.library_id(),
        report: F4CoverageReport::from_base_counts(30.0, 30_000, 30_000, 0, 1_000, 0)
            .expect("coverage"),
    }]));
    let summary = Box::leak(Box::new([EventManifestFinalFastqSummary {
        individual_id: individual,
        library_id: library.library_id(),
        pair_count: 100,
        read_count: 200,
        base_count: 30_000,
        own_pair_count: 100,
        contaminant_pair_count: 0,
        own_base_count: 30_000,
        contaminant_base_count: 0,
        n_base_count: 0,
    }]));
    (
        EventManifest {
            event,
            state: UnitState::Complete,
            identities,
            reference: EventManifestReference {
                reference_id: reference,
                fasta,
            },
            artifacts,
            coverage,
            final_fastq_summary: summary,
            error_quality:
                genome_population_simulator::event_manifest::EventManifestErrorQuality::neutral(),
            demultiplexing: EventManifestDemultiplexing::disabled(),
            chimera_rate: event.chimera_rate(),
            realized_chimera_pairs: 0,
            biological_contamination: &[],
            stutter: None,
            pcr: None,
            diagnostics: &[],
        },
        "reads-main",
    )
}

fn declared() -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let text = b"config_kind = \"general\"\nschema_version = \"2\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1000\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\"]\n\n[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 30\n\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n";
    let admitted = admit_toml_header(
        Path::new("f3-event-manifest.toml"),
        text.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    parse_declared_config(admitted).expect("declared configuration")
}

fn record(id: &str, kind: &str, path: &str, digest: &str) -> ArtifactRecord {
    ArtifactRecord {
        id: id.to_owned(),
        kind: kind.to_owned(),
        path: path.to_owned(),
        size_bytes: 1,
        digest: digest.to_owned(),
    }
}

fn assert_external(path: &Path, accepted: bool) {
    let output = Command::new("python3")
        .arg("tests/external_jsonschema_validator.py")
        .arg("schemas/manifest-event-v2.schema.json")
        .arg(path)
        .output()
        .expect("external JSON Schema validator");
    assert_eq!(
        output.status.success(),
        accepted,
        "external validator for {path:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn validate_document(
    document: &serde_json::Value,
) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
    validate_event_manifest_bytes(
        Path::new("manifest.json"),
        &serde_json::to_vec(document).expect("JSON"),
    )
}

fn failed_document(complete: &serde_json::Value) -> serde_json::Value {
    let mut failed = complete.clone();
    failed["state"] = json!("failed");
    failed["outputs"][0]["state"] = json!("failed");
    failed["outputs"][0]["artifact"] = serde_json::Value::Null;
    failed["outputs"][0]["diagnostic_code"] = json!("GPS-FMT-006");
    failed["coverage"] = json!([]);
    failed["final_fastq_summary"] = json!([]);
    failed["diagnostics"] = json!([format_diagnostic()]);
    failed
}

fn format_diagnostic() -> serde_json::Value {
    json!({
        "code": "GPS-FMT-006",
        "message": "Un documento JSON producido no satisface su esquema publicado.",
    })
}

fn assert_accepted_by_both(root: &Path, label: &str, document: &serde_json::Value) {
    assert!(
        validate_document(document).is_ok(),
        "Rust validator: {label}"
    );
    let path = root.join(format!("{label}.json"));
    fs::write(&path, serde_json::to_vec(document).expect("JSON")).expect("test document");
    assert_external(&path, true);
}

fn assert_rejected_by_both(root: &Path, label: &str, document: &serde_json::Value) {
    assert!(
        validate_document(document).is_err(),
        "Rust validator: {label}"
    );
    let path = root.join(format!("{label}.json"));
    fs::write(&path, serde_json::to_vec(document).expect("JSON")).expect("test document");
    assert_external(&path, false);
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f3-event-manifest-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
