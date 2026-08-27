use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::EffectiveTechnicalEventTarget;
use genome_population_simulator::fastq_output::{FASTQ_ARTIFACT_KIND, write_wgs_fastq};
use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::publication::ExecutionWorkspace;
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::wgs_sampling::WgsLibrarySampler;

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static TEMPORARY_COUNTER: AtomicU64 = AtomicU64::new(0);

fn temporary(name: &str) -> PathBuf {
    let counter = TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!(
        "gps-f3-fastq-{name}-{}-{counter}",
        std::process::id()
    ))
}

fn reference(bases: &[u8]) -> Reference {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "test").expect("sequence ID"),
        FastaHeader::parse(b"chr1").expect("FASTA header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(bases).expect("DNA"),
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn individual(reference: &Reference) -> Individual {
    Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual ID"),
        reference,
        Ploidy::new(1).expect("ploidy"),
    )
}

fn plan(
    reference: &Reference,
    individual: &Individual,
) -> genome_population_simulator::technical_identity::TechnicalIdentityPlan {
    derive_technical_identity_plan(
        RootSeed::parse("root", ROOT).expect("root seed"),
        reference.id(),
        EventId::parse("reads-main", "test").expect("event ID"),
        &[individual.id().clone()],
    )
    .expect("technical identities")
}

fn artifact_bytes(root: &Path, destination: &Path, relative_path: &str) -> Vec<u8> {
    let destination_name = destination
        .file_name()
        .expect("destination name")
        .to_string_lossy();
    let provisional = fs::read_dir(root)
        .expect("temporary root")
        .map(|entry| entry.expect("directory entry").path())
        .find(|path| {
            path.file_name().is_some_and(|name| {
                name.to_string_lossy()
                    .starts_with(&format!(".{destination_name}.gps-incomplete-"))
            })
        })
        .expect("workspace provisional");
    fs::read(provisional.join("payload").join(relative_path)).expect("FASTQ artifact")
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            _ => panic!("test fixture only contains DNA"),
        })
        .collect()
}

#[test]
fn fastq_stream_preserves_literal_unknown_bases_technical_ids_and_phred_values() {
    let root = temporary("unknown-bases");
    fs::create_dir_all(&root).expect("root");
    let mut bases = vec![b'A'; 150];
    bases[0] = b'N';
    bases[76] = b'N';
    bases[149] = b'N';
    let reference = reference(&bases);
    let individual = individual(&reference);
    let plan = plan(&reference, &individual);
    let mut sampler = WgsLibrarySampler::new(
        &reference,
        &individual,
        &plan,
        EffectiveTechnicalEventTarget::ReadCount(1),
    )
    .expect("sampler");
    let candidate_count = sampler.resident_candidate_count();
    let destination = root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");

    let publication = write_wgs_fastq(&mut workspace, &mut sampler).expect("FASTQ");
    let bytes = artifact_bytes(&root, &destination, &publication.artifact().path);
    let library = plan.library(individual.id()).expect("library");
    let expected_header = format!("@TPL_{}_0\n", library.library_id().as_str());
    let expected_forward = bases;
    let expected_reverse = reverse_complement(&expected_forward);
    let lines = bytes
        .split_inclusive(|byte| *byte == b'\n')
        .collect::<Vec<_>>();

    assert_eq!(lines.len(), 4);
    assert_eq!(lines[0], expected_header.as_bytes());
    let observed_sequence = &lines[1][..lines[1].len() - 1];
    assert!(observed_sequence == expected_forward || observed_sequence == expected_reverse);
    assert_eq!(lines[2], b"+\n");
    assert_eq!(lines[3].len(), 151);
    for (base, quality) in observed_sequence.iter().zip(&lines[3][..150]) {
        assert_eq!(*quality, if *base == b'N' { b'!' } else { b'I' });
    }
    assert_eq!(publication.read_count(), 1);
    assert_eq!(publication.base_count(), 150);
    assert_eq!(publication.artifact().kind, FASTQ_ARTIFACT_KIND);
    assert_eq!(sampler.resident_candidate_count(), candidate_count);
    assert_eq!(workspace.artifacts().len(), 1);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn fastq_bytes_are_reproducible_and_the_stream_keeps_only_aggregate_output_state() {
    let mut bases = vec![b'C'; 150];
    bases[32] = b'N';
    let reference = reference(&bases);
    let individual = individual(&reference);
    let plan = plan(&reference, &individual);
    let mut outputs = Vec::new();

    for repeat in 0..2 {
        let root = temporary(&format!("repeated-{repeat}"));
        fs::create_dir_all(&root).expect("root");
        let destination = root.join("published");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let mut sampler = WgsLibrarySampler::new(
            &reference,
            &individual,
            &plan,
            EffectiveTechnicalEventTarget::ReadCount(7),
        )
        .expect("sampler");
        let candidates = sampler.resident_candidate_count();
        let publication = write_wgs_fastq(&mut workspace, &mut sampler).expect("FASTQ");
        assert_eq!(publication.read_count(), 7);
        assert_eq!(publication.base_count(), 7 * 150);
        assert_eq!(sampler.resident_candidate_count(), candidates);
        let bytes = artifact_bytes(&root, &destination, &publication.artifact().path);
        let library = plan.library(individual.id()).expect("library");
        let headers = bytes
            .split_inclusive(|byte| *byte == b'\n')
            .step_by(4)
            .collect::<Vec<_>>();
        assert_eq!(headers.len(), 7);
        for (ordinal, header) in headers.iter().enumerate() {
            assert_eq!(
                *header,
                format!("@TPL_{}_{}\n", library.library_id().as_str(), ordinal).as_bytes()
            );
        }
        outputs.push((bytes, publication.artifact().digest.clone(), root));
    }

    assert_eq!(outputs[0].0, outputs[1].0);
    assert_eq!(outputs[0].1, outputs[1].1);
    for (_, _, root) in outputs {
        fs::remove_dir_all(root).expect("cleanup");
    }
}

#[test]
fn failed_temporary_validation_removes_the_unfinalized_artifact() {
    let root = temporary("failed-temporary");
    fs::create_dir_all(&root).expect("root");
    let destination = root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let diagnostic = workspace
        .write_streamed_artifact_with_validation(
            "synthetic-failure",
            "test",
            Path::new("outputs/failure.fastq"),
            |file| {
                use std::io::Write;
                file.write_all(b"partial\n").map_err(|error| {
                    genome_population_simulator::diagnostics::Diagnostic::internal_invariant(
                        "test", "write",
                    )
                    .with("detail", error.to_string())
                })
            },
            |_| {
                Err(
                    genome_population_simulator::diagnostics::Diagnostic::internal_invariant(
                        "test", "validate",
                    ),
                )
            },
        )
        .expect_err("validation must fail");

    assert_eq!(
        diagnostic.context().get("operation"),
        Some(&"validate".to_owned())
    );
    assert!(workspace.artifacts().is_empty());
    let payload = fs::read_dir(&root)
        .expect("root entries")
        .map(|entry| entry.expect("entry").path())
        .find(|path| path.is_dir())
        .expect("provisional")
        .join("payload/outputs");
    assert!(
        !payload.exists()
            || fs::read_dir(payload)
                .expect("output entries")
                .next()
                .is_none()
    );
    fs::remove_dir_all(root).expect("cleanup");
}
