use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::demultiplexing::{
    Demultiplexer, DemultiplexingLibrary, PostContaminationPair, Udi,
};
use genome_population_simulator::error_quality::{ErrorQualitySettings, transform_paired_end};
use genome_population_simulator::event_manifest::{
    EventManifest, EventManifestArtifact, EventManifestCoverage, EventManifestDemultiplexing,
    EventManifestFinalFastqSummary, EventManifestReference, event_manifest_document,
    validate_event_manifest_bytes,
};
use genome_population_simulator::fastq_output::write_final_paired_fastq;
use genome_population_simulator::fastq_validation::F4CoverageReport;
use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::publication::ExecutionWorkspace;
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::results::UnitState;
use genome_population_simulator::technical_identity::{
    TechnicalMate, derive_technical_identity_plan,
};
use genome_population_simulator::technical_profile::resolve_f4_technical_profile;
use genome_population_simulator::wgs_sampling::WgsPairedEndSampler;
use sha2::{Digest, Sha256};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static TEMPORARY_COUNTER: AtomicU64 = AtomicU64::new(0);

#[test]
fn final_paired_fastq_is_gzip_hashed_paired_and_manifested_without_individual_provenance() {
    let root = temporary("final");
    fs::create_dir_all(&root).expect("root");
    let effective = build_effective_config(&declared()).expect("effective configuration");
    let event = effective.technical_event().expect("technical event");
    let reference = reference();
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual"),
        &reference,
        Ploidy::new(1).expect("ploidy"),
    );
    let seed = RootSeed::parse("root", ROOT).expect("root seed");
    let identities = derive_technical_identity_plan(
        seed,
        reference.id(),
        event.event_id().clone(),
        &[individual.id().clone()],
    )
    .expect("identities");
    let library = identities.library(individual.id()).expect("library");
    let profile = resolve_f4_technical_profile(event.profile()).expect("profile");
    let pair = WgsPairedEndSampler::new_for_pair_count(
        &reference,
        &individual,
        &identities,
        profile,
        event.read_length(),
        1,
    )
    .expect("sampler")
    .next()
    .expect("one pair")
    .expect("sampled pair");
    let observed = transform_paired_end(
        seed,
        event.event_id(),
        &pair,
        profile,
        ErrorQualitySettings::default(),
    )
    .expect("observed pair")
    .observed()
    .clone();
    let udi = Udi::new("ACGTACGT", "TGCATGCA").expect("UDI");
    let mut demultiplexer = Demultiplexer::new(
        seed,
        event.event_id().clone(),
        0.0,
        vec![DemultiplexingLibrary::new(
            library.library_id().clone(),
            udi.clone(),
        )],
    )
    .expect("demultiplexer");
    let outcome = demultiplexer
        .demultiplex(vec![PostContaminationPair::new(
            library.library_id().clone(),
            pair.pre_hop_identity().pre_hop_instance_id().clone(),
            udi,
            false,
            observed,
        )])
        .expect("demultiplex");
    let demultiplexing = EventManifestDemultiplexing::from_report(outcome.report().clone());
    let pairs = outcome
        .into_final_libraries()
        .remove(library.library_id())
        .expect("final library");

    let destination = root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let fasta = workspace
        .write_artifact(
            "reference-fasta",
            "reference_fasta",
            Path::new("resources/reference.fa"),
            b">chr1\nA\n",
        )
        .expect("FASTA");
    let publication = write_final_paired_fastq(
        &mut workspace,
        event.event_id().as_str(),
        individual.id().as_str(),
        library.library_id(),
        pairs,
    )
    .expect("final FASTQ");

    // Oráculo independiente: gzip del sistema, gramática FASTQ, pares y
    // SHA-256 de los bytes comprimidos, sin reutilizar el validador del núcleo.
    let r1 = artifact_path(&root, &destination, &publication.r1().path);
    let r2 = artifact_path(&root, &destination, &publication.r2().path);
    let r1_records = independently_read_fastq(&r1, TechnicalMate::R1);
    let r2_records = independently_read_fastq(&r2, TechnicalMate::R2);
    assert_eq!(r1_records, r2_records, "no final paired-end orphans");
    assert_eq!(r1_records.len(), 1);
    assert_eq!(sha256_file(&r1), publication.r1().digest);
    assert_eq!(sha256_file(&r2), publication.r2().digest);
    assert_eq!(
        fs::metadata(&r1).expect("R1 metadata").len(),
        publication.r1().size_bytes
    );
    assert_eq!(
        fs::metadata(&r2).expect("R2 metadata").len(),
        publication.r2().size_bytes
    );
    assert_eq!(
        workspace.artifacts().len(),
        3,
        "only FASTQ R1/R2 and the referenced FASTA are durable"
    );
    assert!(
        workspace
            .artifacts()
            .iter()
            .all(|artifact| !artifact.path.contains("clean")
                && !artifact.path.contains("provenance"))
    );

    let artifacts = [
        EventManifestArtifact {
            id: &publication.r1().id,
            kind: "fastq",
            state: UnitState::Complete,
            individual_id: individual.id(),
            library_id: library.library_id(),
            mate: TechnicalMate::R1,
            artifact: Some(publication.r1()),
            diagnostic: None,
        },
        EventManifestArtifact {
            id: &publication.r2().id,
            kind: "fastq",
            state: UnitState::Complete,
            individual_id: individual.id(),
            library_id: library.library_id(),
            mate: TechnicalMate::R2,
            artifact: Some(publication.r2()),
            diagnostic: None,
        },
    ];
    let summary = [EventManifestFinalFastqSummary {
        individual_id: individual.id(),
        library_id: library.library_id(),
        pair_count: publication.pair_count(),
        read_count: publication.read_count(),
        base_count: publication.base_count(),
        own_pair_count: publication.own_pair_count(),
        contaminant_pair_count: publication.contaminant_pair_count(),
        own_base_count: publication.own_base_count(),
        contaminant_base_count: publication.contaminant_base_count(),
        n_base_count: publication.n_base_count(),
    }];
    let coverage = [EventManifestCoverage {
        individual_id: individual.id(),
        library_id: library.library_id(),
        report: F4CoverageReport::from_base_counts(
            1.0,
            publication.own_base_count(),
            publication.own_base_count(),
            publication.contaminant_base_count(),
            reference.total_length(),
            publication.n_base_count(),
        )
        .expect("coverage"),
    }];
    let manifest = EventManifest {
        event,
        state: UnitState::Complete,
        identities: &identities,
        reference: EventManifestReference {
            reference_id: reference.id(),
            fasta: &fasta,
        },
        artifacts: &artifacts,
        coverage: &coverage,
        final_fastq_summary: &summary,
        error_quality:
            genome_population_simulator::event_manifest::EventManifestErrorQuality::neutral(),
        demultiplexing,
        chimera_rate: event.chimera_rate(),
        realized_chimera_pairs: 0,
        biological_contamination: &[],
        stutter: None,
        diagnostics: &[],
    };
    let document = event_manifest_document(&manifest).expect("event manifest");
    assert_eq!(document["final_fastq_summary"][0]["pair_count"], 1);
    assert_eq!(document["demultiplexing"]["clean_fastq_published"], false);
    assert_eq!(
        document["demultiplexing"]["individual_provenance_published"],
        false
    );
    let serialized = serde_json::to_string(&document).expect("manifest JSON");
    assert!(!serialized.contains("PH_"));
    assert!(!serialized.contains("technical_pair_id"));
    assert!(!serialized.contains("pair_ids"));
    let mut inconsistent_summary = document.clone();
    inconsistent_summary["final_fastq_summary"][0]["own_base_count"] =
        serde_json::json!(publication.own_base_count() + 1);
    assert_manifest_rejected(inconsistent_summary, "FASTQ summary must remain coherent");

    let mut source_outside_event = document.clone();
    source_outside_event["demultiplexing"] = serde_json::json!({
        "index_hopping_rate": 1.0,
        "reassigned_pairs": 1,
        "discarded_pairs": 0,
        "unassigned_pairs": 0,
        "reassigned_by_sample": [{
            "source_library_id": "LIB_outside_event",
            "destination_library_id": library.library_id().as_str(),
            "cause": "index_hopping",
            "pair_count": 1,
        }],
        "discarded_by_sample": [],
        "unassigned_by_sample": [],
        "clean_fastq_published": false,
        "individual_provenance_published": false,
    });
    assert_manifest_rejected(source_outside_event, "source library outside the event");

    let mut destination_outside_event = document.clone();
    destination_outside_event["demultiplexing"] = serde_json::json!({
        "index_hopping_rate": 1.0,
        "reassigned_pairs": 1,
        "discarded_pairs": 0,
        "unassigned_pairs": 0,
        "reassigned_by_sample": [{
            "source_library_id": library.library_id().as_str(),
            "destination_library_id": "LIB_outside_event",
            "cause": "index_hopping",
            "pair_count": 1,
        }],
        "discarded_by_sample": [],
        "unassigned_by_sample": [],
        "clean_fastq_published": false,
        "individual_provenance_published": false,
    });
    assert_manifest_rejected(
        destination_outside_event,
        "destination library outside the event",
    );

    let mut reassigned_over_final = document.clone();
    reassigned_over_final["demultiplexing"]["index_hopping_rate"] = serde_json::json!(1.0);
    reassigned_over_final["demultiplexing"]["reassigned_pairs"] = serde_json::json!(2);
    assert_manifest_rejected(reassigned_over_final, "reassignments exceed final pairs");

    let mut zero_rate_reassignment = document.clone();
    zero_rate_reassignment["demultiplexing"]["reassigned_pairs"] = serde_json::json!(1);
    assert_manifest_rejected(zero_rate_reassignment, "zero hopping rate reassigns a pair");

    let mut unreconciled_discards = document;
    unreconciled_discards["demultiplexing"]["discarded_pairs"] = serde_json::json!(2);
    unreconciled_discards["demultiplexing"]["unassigned_pairs"] = serde_json::json!(1);
    unreconciled_discards["demultiplexing"]["discarded_by_sample"] = serde_json::json!([{
        "source_library_id": library.library_id().as_str(),
        "cause": "unrecognized_udi",
        "pair_count": 2,
    }]);
    unreconciled_discards["demultiplexing"]["unassigned_by_sample"] = serde_json::json!([{
        "source_library_id": library.library_id().as_str(),
        "cause": "unrecognized_udi",
        "pair_count": 1,
    }]);
    assert_manifest_rejected(unreconciled_discards, "discard and unassigned maps differ");
    fs::remove_dir_all(root).expect("cleanup");
}

fn assert_manifest_rejected(document: serde_json::Value, detail: &str) {
    assert!(
        validate_event_manifest_bytes(
            Path::new("events/final-main/manifest.json"),
            &serde_json::to_vec(&document).expect("mutated JSON"),
        )
        .is_err(),
        "{detail} must be rejected"
    );
}

fn independently_read_fastq(path: &Path, mate: TechnicalMate) -> Vec<String> {
    let output = std::process::Command::new("gzip")
        .args(["-cd", path.to_str().expect("UTF-8 path")])
        .output()
        .expect("gzip installed");
    assert!(
        output.status.success(),
        "gzip decoder failed: {:?}",
        output.stderr
    );
    let lines = output
        .stdout
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
        .collect::<Vec<_>>();
    assert_eq!(lines.len() % 4, 0);
    lines
        .chunks_exact(4)
        .map(|record| {
            let header = std::str::from_utf8(record[0]).expect("header UTF-8");
            assert!(header.starts_with("@GPS:final-main:"));
            assert!(header.contains(&format!(" {}:N:0:ACGTACGT+TGCATGCA", mate.number())));
            assert_eq!(record[2], b"+");
            assert_eq!(record[1].len(), record[3].len());
            assert!(
                record[1]
                    .iter()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
            );
            assert!(
                record[3]
                    .iter()
                    .all(|quality| (b'!'..=b'~').contains(quality))
            );
            header
                .split(' ')
                .next()
                .and_then(|identity| identity.rsplit(':').next())
                .expect("pair ID")
                .to_owned()
        })
        .collect()
}

fn sha256_file(path: &Path) -> String {
    let mut file = fs::File::open(path).expect("artifact");
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 4096];
    loop {
        let read = file.read(&mut buffer).expect("read");
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    format!("{:x}", hasher.finalize())
}

fn artifact_path(root: &Path, destination: &Path, relative_path: &str) -> PathBuf {
    let destination_name = destination
        .file_name()
        .expect("destination name")
        .to_string_lossy();
    let provisional = fs::read_dir(root)
        .expect("workspace root")
        .map(|entry| entry.expect("entry").path())
        .find(|path| {
            path.file_name().is_some_and(|name| {
                name.to_string_lossy()
                    .starts_with(&format!(".{destination_name}.gps-incomplete-"))
            })
        })
        .expect("workspace provisional");
    provisional.join("payload").join(relative_path)
}

fn reference() -> Reference {
    let bases = (0..500).map(|index| b"ACGT"[index % 4]).collect::<Vec<_>>();
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "test").expect("sequence"),
        FastaHeader::parse(b"chr1").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(&bases).expect("DNA"),
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn declared() -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let bytes = b"config_kind = \"general\"\nschema_version = \"2\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 500\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 1\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\"]\n\n[technical_event]\nevent_id = \"final-main\"\nindividual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 1\n\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n";
    parse_declared_config(
        admit_toml_header(
            Path::new("final.toml"),
            bytes.to_vec(),
            ExpectedConfigKind::General,
        )
        .expect("header"),
    )
    .expect("declared")
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-final-{label}-{}-{}",
        std::process::id(),
        TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed)
    ))
}
