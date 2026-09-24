//! Caso controlado de F4.2.3: ninguna de las tres coberturas sustituye a otra.

use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::event_manifest::{
    EventManifest, EventManifestArtifact, EventManifestCoverage, EventManifestDemultiplexing,
    EventManifestReference, event_manifest_document, validate_event_manifest_bytes,
};
use genome_population_simulator::fastq_validation::F4CoverageReport;
use genome_population_simulator::identity::ReferenceId;
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::results::{ArtifactRecord, UnitState};
use genome_population_simulator::technical_identity::{
    TechnicalMate, derive_technical_identity_plan,
};
use serde_json::json;

const ROOT_SEED: &str = "000102030405060708090a0b1c1d1e1f101112131415161718191a1b1c1d1e1f";
const DIGEST_A: &str = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
const DIGEST_B: &str = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";

#[test]
fn target_basic_and_final_own_contaminant_coverages_are_published_separately() {
    let effective = build_effective_config(&declared()).expect("effective configuration");
    let event = effective.technical_event().expect("technical event");
    let reference_id = ReferenceId::parse("reference-main", "test").expect("reference ID");
    let identities = derive_technical_identity_plan(
        RootSeed::parse("root", ROOT_SEED).expect("root seed"),
        &reference_id,
        event.event_id().clone(),
        event.individual_ids(),
    )
    .expect("identities");
    let individual = event.individual_ids().first().expect("individual");
    let library = identities.library(individual).expect("library");
    let reference = record(
        "reference-main-fasta",
        "reference_fasta",
        "resources/references/reference-main.fa",
        DIGEST_A,
    );
    let r1 = record(
        "fastq-LIB_sample-A",
        "fastq",
        "outputs/fastq/coverage-main/sample-A_R1.fastq.gz",
        DIGEST_A,
    );
    let r2 = record(
        "fastq-LIB_sample-A-r2",
        "fastq",
        "outputs/fastq/coverage-main/sample-A_R2.fastq.gz",
        DIGEST_B,
    );
    let artifacts = [
        EventManifestArtifact {
            id: &r1.id,
            kind: "fastq",
            state: UnitState::Complete,
            individual_id: individual,
            library_id: library.library_id(),
            mate: TechnicalMate::R1,
            artifact: Some(&r1),
            diagnostic: None,
        },
        EventManifestArtifact {
            id: &r2.id,
            kind: "fastq",
            state: UnitState::Complete,
            individual_id: individual,
            library_id: library.library_id(),
            mate: TechnicalMate::R2,
            artifact: Some(&r2),
            diagnostic: None,
        },
    ];
    // Oráculo controlado: 10× solicitado, 13× tras WGS básica y 9× propias
    // más 6× contaminantes en el FASTQ final. Cada pareja mantiene 2×150 pb.
    let report =
        F4CoverageReport::from_base_counts(10.0, 1_300, 900, 600, 100, 3).expect("coverage report");
    let coverage = [EventManifestCoverage {
        individual_id: individual,
        library_id: library.library_id(),
        report,
    }];
    let summary = [
        genome_population_simulator::event_manifest::EventManifestFinalFastqSummary {
            individual_id: individual,
            library_id: library.library_id(),
            pair_count: 5,
            read_count: 10,
            base_count: 1_500,
            own_pair_count: 3,
            contaminant_pair_count: 2,
            own_base_count: 900,
            contaminant_base_count: 600,
            n_base_count: 3,
        },
    ];
    let manifest = EventManifest {
        event,
        state: UnitState::Complete,
        identities: &identities,
        reference: EventManifestReference {
            reference_id: &reference_id,
            fasta: &reference,
        },
        artifacts: &artifacts,
        coverage: &coverage,
        final_fastq_summary: &summary,
        error_quality:
            genome_population_simulator::event_manifest::EventManifestErrorQuality::neutral(),
        demultiplexing: EventManifestDemultiplexing::disabled(),
        chimera_rate: event.chimera_rate(),
        realized_chimera_pairs: 4,
        biological_contamination: &[],
        stutter: None,
        diagnostics: &[],
    };

    let document = event_manifest_document(&manifest).expect("event manifest");
    let coverage = &document["coverage"][0];
    assert_eq!(coverage["target_x"], json!(10.0));
    assert_eq!(coverage["basic_x"], json!(13.0));
    assert_eq!(coverage["final_x"]["own_x"], json!(9.0));
    assert_eq!(coverage["final_x"]["contaminant_x"], json!(6.0));
    assert_eq!(coverage["denominator_bases"], json!(100));
    assert_eq!(coverage["n_bases"], json!(3));
    assert_eq!(document["artifacts"]["chimera"]["rate"], json!(1.0));
    assert_eq!(document["artifacts"]["chimera"]["realized_pairs"], json!(4));
    assert!(document["artifacts"]["chimera"].get("pair_ids").is_none());
    assert!(
        document["artifacts"]["chimera"]
            .get("coordinates")
            .is_none()
    );

    // `N` pertenece a las bases FASTQ observadas: a 2× puede superar el
    // denominador de referencia sin que deje de ser una cobertura válida.
    assert_eq!(
        F4CoverageReport::from_base_counts(10.0, 200, 200, 0, 100, 150)
            .expect("N coverage above structural denominator")
            .n_bases(),
        150
    );

    let bytes = serde_json::to_vec(&document).expect("JSON");
    validate_event_manifest_bytes(Path::new("events/coverage-main/manifest.json"), &bytes)
        .expect("manifest validator");

    let mut collapsed_final = document;
    collapsed_final["coverage"][0]["final_x"] = json!(12.0);
    assert!(
        validate_event_manifest_bytes(
            Path::new("events/coverage-main/manifest.json"),
            &serde_json::to_vec(&collapsed_final).expect("JSON"),
        )
        .is_err(),
        "the final FASTQ coverage has mandatory own and contaminant components"
    );
}

fn declared() -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let text = b"config_kind = \"general\"\nschema_version = \"2\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 100\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\"]\n\n[technical_event]\nevent_id = \"coverage-main\"\nindividual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 10\nchimera_rate = 1\n\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n";
    let admitted = admit_toml_header(
        Path::new("f4-coverage-reporting.toml"),
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
