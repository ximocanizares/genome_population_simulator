use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{build_effective_config, load_and_validate_config};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::event_manifest::{
    EventManifest, EventManifestArtifact, EventManifestDemultiplexing, EventManifestReference,
    event_manifest_bytes, event_manifest_relative_path,
};
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::ReferenceId;
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::publication::ExecutionWorkspace;
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::results::{
    GeneralManifest, ManifestEventLink, ManifestReferenceLink, UnitState,
    general_manifest_document, json_document_bytes, validate_published_event_manifests,
    validate_published_resources,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::{
    adapters::{admit_toml_header, parse_declared_config},
    execution::ExecutionAuthorizations,
};

static NONCE: AtomicU64 = AtomicU64::new(0);
const ROOT_SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

#[test]
fn complete_execution_assembles_and_links_the_validated_event_manifest() {
    let root = temporary("complete");
    fs::create_dir_all(&root).expect("root");
    let config = root.join("config.toml");
    fs::write(&config, configuration()).expect("config");
    let destination = root.join("published");

    execute_validated_config(
        &load_and_validate_config(&config).expect("effective config"),
        &destination,
    )
    .expect("complete publication");

    let manifest_path = destination.join("manifest.json");
    let manifest_bytes = fs::read(&manifest_path).expect("general manifest");
    let manifest: serde_json::Value = serde_json::from_slice(&manifest_bytes).expect("JSON");
    assert_eq!(manifest["state"], "complete");
    let events = manifest["events"].as_array().expect("real event link");
    assert_eq!(events.len(), 1);
    assert_eq!(events[0]["event_id"], "reads-main");
    assert_eq!(events[0]["state"], "complete");
    let event_path = destination.join("events/reads-main/manifest.json");
    assert!(event_path.is_file());
    let event: serde_json::Value =
        serde_json::from_slice(&fs::read(&event_path).expect("event bytes")).expect("event JSON");
    assert_eq!(event["state"], "complete");
    let fastq = event["outputs"].as_array().expect("FASTQ rows");
    assert_eq!(fastq.len(), 4);
    for artifact in fastq {
        assert_eq!(artifact["state"], "complete");
        let path = artifact["artifact"]["path"].as_str().expect("FASTQ path");
        assert!(destination.join(path).is_file());
    }
    validate_published_event_manifests(&destination, &manifest_path, &manifest_bytes)
        .expect("physical event link");
    let stray_path = destination.join("outputs/fastq/reads-main/stray.fastq.gz");
    fs::write(&stray_path, b"@stray\nA\n+\nI\n").expect("stray FASTQ");
    let stray_identity = identify_file(
        &stray_path,
        Path::new("outputs/fastq/reads-main/stray.fastq.gz"),
    )
    .expect("stray identity");
    let mut with_stray = manifest.clone();
    with_stray["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .push(serde_json::json!({
            "id": "fastq-stray",
            "kind": "fastq",
            "path": "outputs/fastq/reads-main/stray.fastq.gz",
            "size_bytes": stray_identity.size_bytes(),
            "checksum": {"algorithm": "sha256", "digest": stray_identity.digest()},
        }));
    with_stray["output_categories"]["minimum"]
        .as_array_mut()
        .expect("minimum")
        .push(serde_json::json!("fastq-stray"));
    assert!(
        validate_published_event_manifests(
            &destination,
            &manifest_path,
            &serde_json::to_vec(&with_stray).expect("stray manifest"),
        )
        .is_err(),
        "a physical FASTQ must be linked by exactly one complete event artifact"
    );
    let mut unlinked = manifest.clone();
    unlinked
        .as_object_mut()
        .expect("manifest object")
        .remove("events");
    let removed_ids = unlinked["artifacts"]
        .as_array()
        .expect("catalog")
        .iter()
        .filter(|artifact| matches!(artifact["kind"].as_str(), Some("fastq" | "event_manifest")))
        .filter_map(|artifact| artifact["id"].as_str().map(str::to_owned))
        .collect::<Vec<_>>();
    unlinked["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .retain(|artifact| !removed_ids.iter().any(|id| artifact["id"] == *id));
    unlinked["output_categories"]["minimum"]
        .as_array_mut()
        .expect("minimum")
        .retain(|id| !removed_ids.iter().any(|removed| *id == *removed));
    assert!(
        validate_published_resources(
            &destination,
            &destination.join("effective-config.json"),
            &fs::read(destination.join("effective-config.json")).expect("effective bytes"),
            &manifest_path,
            &serde_json::to_vec(&unlinked).expect("unlinked manifest"),
            None,
        )
        .is_err()
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
#[ignore = "F4.1 retira el ensamblado/publicación FASTQ SE de F3"]
fn controlled_failed_event_keeps_only_closed_fastq_and_its_failed_manifest() {
    let root = temporary("failed");
    fs::create_dir_all(&root).expect("root");
    let destination = root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let reference = workspace
        .write_artifact(
            "reference-main-fasta",
            "reference_fasta",
            Path::new("resources/references/reference-main.fa"),
            b">chr1\nACGT\n",
        )
        .expect("reference");
    let effective = effective_event();
    let event = effective.technical_event().expect("technical event");
    let reference_id = ReferenceId::parse("reference-main", "test").expect("reference ID");
    let identities = derive_technical_identity_plan(
        RootSeed::parse("root", ROOT_SEED).expect("seed"),
        &reference_id,
        event.event_id().clone(),
        event.individual_ids(),
    )
    .expect("identities");
    let libraries = identities.libraries().collect::<Vec<_>>();
    let library = libraries.first().expect("first library");
    let missing_library = libraries.get(1).expect("second library");
    let fastq = workspace
        .write_artifact(
            &format!("fastq-{}", library.library_id().as_str()),
            "fastq",
            Path::new("outputs/fastq/reads-main/sample-A.fastq.gz"),
            b"@read\nA\n+\nI\n",
        )
        .expect("closed independent FASTQ");
    let failure = Diagnostic::internal_invariant("test", "controlled-event-failure");
    let failed_artifacts = [
        EventManifestArtifact {
            id: &fastq.id,
            kind: "fastq",
            state: UnitState::Complete,
            individual_id: library.individual_id(),
            library_id: library.library_id(),
            mate: genome_population_simulator::technical_identity::TechnicalMate::R1,
            artifact: Some(&fastq),
            diagnostic: None,
        },
        EventManifestArtifact {
            id: "fastq-unclosed-library",
            kind: "fastq",
            state: UnitState::Failed,
            individual_id: missing_library.individual_id(),
            library_id: missing_library.library_id(),
            mate: genome_population_simulator::technical_identity::TechnicalMate::R1,
            artifact: None,
            diagnostic: Some(&failure),
        },
    ];
    let event_document = EventManifest {
        event,
        state: UnitState::Failed,
        identities: &identities,
        reference: EventManifestReference {
            reference_id: &reference_id,
            fasta: &reference,
        },
        artifacts: &failed_artifacts,
        coverage: &[],
        final_fastq_summary: &[],
        error_quality:
            genome_population_simulator::event_manifest::EventManifestErrorQuality::neutral(),
        demultiplexing: EventManifestDemultiplexing::disabled(),
        chimera_rate: event.chimera_rate(),
        realized_chimera_pairs: 0,
        biological_contamination: &[],
        stutter: None,
        pcr: None,
        diagnostics: std::slice::from_ref(&failure),
    };
    let bytes = event_manifest_bytes(&event_document).expect("failed event bytes");
    let event_record = workspace
        .write_streamed_artifact_with_validation(
            "event-manifest-reads-main",
            "event_manifest",
            &event_manifest_relative_path(event.event_id()),
            |file| {
                use std::io::Write;
                file.write_all(&bytes).map_err(|error| {
                    Diagnostic::internal_invariant("test", "write-failed-event-manifest")
                        .with("detail", error.to_string())
                })
            },
            |path| {
                let bytes = fs::read(path).map_err(|error| {
                    Diagnostic::internal_invariant("test", "read-failed-event-manifest")
                        .with("detail", error.to_string())
                })?;
                genome_population_simulator::event_manifest::validate_event_manifest_bytes(
                    Path::new("events/reads-main/manifest.json"),
                    &bytes,
                )
            },
        )
        .expect("closed failed event manifest");
    let events = [ManifestEventLink {
        event_id: "reads-main".to_owned(),
        state: UnitState::Failed,
        manifest_artifact_id: event_record.id,
        chimera_rate: event.chimera_rate(),
        realized_chimera_pairs: 0,
        profile_id: genome_population_simulator::technical_profile::GPS_HISEQX_PCRFREE_PE150_V1,
        profile_sha256:
            genome_population_simulator::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256,
        rng_catalog: genome_population_simulator::randomness::F4_TECHNICAL_RNG_CATALOG_V3,
    }];
    let references = [ManifestReferenceLink {
        id: "reference-main".to_owned(),
        fasta_artifact_id: reference.id,
    }];
    let catalog = workspace.artifacts().to_vec();
    let general = GeneralManifest {
        execution_id: "failed-event-execution",
        state: UnitState::Failed,
        configuration_input: None,
        effective_config: None,
        references: &references,
        reduced_references: &[],
        population_snvs: None,
        vcf_outputs: &[],
        events: Some(&events),
        authorizations: ExecutionAuthorizations::default(),
        artifacts: &catalog,
        diagnostics: std::slice::from_ref(&failure),
        repeat_loci: None,
    };
    workspace
        .write_general_manifest(
            &json_document_bytes(&general_manifest_document(&general)).expect("general bytes"),
        )
        .expect("failed general manifest");
    workspace.publish().expect("atomic failed publication");
    let general_path = destination.join("manifest.json");
    let general_bytes = fs::read(&general_path).expect("general bytes");
    validate_published_event_manifests(&destination, &general_path, &general_bytes)
        .expect("failed event remains coherent");
    assert!(destination.join(&fastq.path).is_file());
    fs::remove_dir_all(root).expect("cleanup");
}

fn effective_event() -> genome_population_simulator::application::EffectiveConfig {
    let admitted = admit_toml_header(
        Path::new("f3-event-assembly.toml"),
        configuration().into_bytes(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    build_effective_config(&parse_declared_config(admitted).expect("declared")).expect("effective")
}

fn configuration() -> String {
    format!(
        r#"config_kind = "general"
schema_version = "2"
seed = "{ROOT_SEED}"
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
[reference]
kind = "simulated"
id = "reference-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 500
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["sample-A", "sample-B"]
ploidy = 2
genotype_model = "hwe"
[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A", "sample-B"]
target_coverage_wgs = 1
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
[[technical_event.libraries]]
individual_id = "sample-B"
i7 = "GGGGGGGG"
i5 = "CCCCCCCC"
"#
    )
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f3-event-assembly-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
    ))
}
