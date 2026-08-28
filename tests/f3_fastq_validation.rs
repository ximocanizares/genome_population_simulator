use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::EffectiveTechnicalEventTarget;
use genome_population_simulator::coordinates::BasePosition;
use genome_population_simulator::fastq_output::write_wgs_fastq;
use genome_population_simulator::fastq_output::{FASTQ_ARTIFACT_KIND, fastq_relative_path};
use genome_population_simulator::fastq_validation::{FastqValidationError, validate_wgs_fastq};
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::materialization::SnpLocus;
use genome_population_simulator::publication::ExecutionWorkspace;
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::CanonicalBase;
use genome_population_simulator::snp_application::{ApplicableSnp, TargetedSnp, apply_snps};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::wgs_sampling::WgsLibrarySampler;
use noodles_bgzf as bgzf;

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static TEMPORARY_COUNTER: AtomicU64 = AtomicU64::new(0);

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f3-fastq-validation-{}-{}-{label}",
        std::process::id(),
        TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed)
    ))
}

fn reference() -> Reference {
    // Verdad sintética escrita a mano: los extremos y el interior contienen N
    // y la secuencia no es palindrómica para que cambiar orientación cambie la
    // lectura observada.
    // Referencia congelada, mayor que una lectura y no periódica: cada
    // ventana posible es distinguible. Conserva N en extremos e interior.
    let bases = b"NGCTGTTTCCTATGCTACTCACTTGGGGGCTAGACACGACTGGTACAAGACGTTTGAACTTCCTCGTCGCACCNCATGAGAGCGTAAAGGCGTACTTAATAGTATGCGTCCCCTCAAAGCCCTGAACCGACTATCAATCCCGCCAAGTANAAGCCATATACGGACGTTTGCGCGTATTCTTTTGGTGGCCGAATTCCGATGTGAATGTCATTCGAACGAAACTCTAGTACTTGCTTATAAAATCGATGCCGTTTGGN".to_vec();
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr-synthetic", "test").expect("sequence ID"),
        FastaHeader::parse(b"chr-synthetic").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(&bases).expect("synthetic DNA"),
    )
    .expect("sequence");
    Reference::new(
        ReferenceId::parse("reference-synthetic", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn truth_state(reference: &Reference) -> Individual {
    Individual::from_reference(
        IndividualId::parse("sample-synthetic", "test").expect("individual ID"),
        reference,
        Ploidy::new(1).expect("ploidy"),
    )
}

fn differentiated_truth_state(reference: &Reference) -> Individual {
    let baseline = Individual::from_reference(
        IndividualId::parse("sample-synthetic", "test").expect("individual ID"),
        reference,
        Ploidy::DIPLOID,
    );
    let sequence = &reference.sequences()[0];
    let alternate = |base| match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        b'T' => b'A',
        b'N' => b'A',
        _ => panic!("synthetic truth contains only DNA"),
    };
    let applications = [10_u64, 35, 60, 85, 110, 135]
        .into_iter()
        .map(|offset| {
            let base = sequence.sequence().as_bytes()[usize::try_from(offset).expect("offset")];
            let locus = SnpLocus::new(
                BasePosition::new(
                    reference.id().clone(),
                    sequence.id().clone(),
                    offset,
                    sequence.length(),
                )
                .expect("position"),
                CanonicalBase::from_byte(base).expect("canonical reference base"),
                CanonicalBase::from_byte(alternate(base)).expect("canonical alternate base"),
            )
            .expect("locus");
            TargetedSnp::new(
                HaplotypePhase::from_copy_index(1).expect("phase"),
                ApplicableSnp::explicit(
                    locus,
                    EventId::parse("fixture-snp", "test").expect("event ID"),
                ),
            )
        })
        .collect::<Vec<_>>();
    apply_snps(&baseline, reference, &applications).expect("phased SNP fixture")
}

fn all_unknown_truth_reference() -> Reference {
    let bases = vec![b'N'; 150];
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr-unknown", "test").expect("sequence ID"),
        FastaHeader::parse(b"chr-unknown").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(&bases).expect("synthetic DNA"),
    )
    .expect("sequence");
    Reference::new(
        ReferenceId::parse("reference-unknown", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn identities(
    reference: &Reference,
    individual: &Individual,
) -> genome_population_simulator::technical_identity::TechnicalIdentityPlan {
    derive_technical_identity_plan(
        RootSeed::parse("root", ROOT).expect("root seed"),
        reference.id(),
        EventId::parse("event-synthetic", "test").expect("event ID"),
        &[individual.id().clone()],
    )
    .expect("identities")
}

fn write_fastq(
    root: &Path,
    reference: &Reference,
    individual: &Individual,
    identities: &genome_population_simulator::technical_identity::TechnicalIdentityPlan,
    target: EffectiveTechnicalEventTarget,
) -> (
    PathBuf,
    genome_population_simulator::results::ArtifactRecord,
) {
    fs::create_dir_all(root).expect("root");
    let destination = root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let mut sampler =
        WgsLibrarySampler::new(reference, individual, identities, target).expect("sampler");
    let publication = write_wgs_fastq(&mut workspace, &mut sampler).expect("FASTQ");
    let destination_name = destination
        .file_name()
        .expect("destination name")
        .to_string_lossy();
    let provisional = fs::read_dir(root)
        .expect("temporary root")
        .map(|entry| entry.expect("entry").path())
        .find(|path| {
            path.file_name().is_some_and(|name| {
                name.to_string_lossy()
                    .starts_with(&format!(".{destination_name}.gps-incomplete-"))
            })
        })
        .expect("workspace provisional");
    (
        provisional
            .join("payload")
            .join(&publication.artifact().path),
        publication.artifact().clone(),
    )
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
            _ => panic!("synthetic truth contains only ACGTN"),
        })
        .collect()
}

fn write_independent_unknown_fastq(
    path: &Path,
    individual: &Individual,
    identities: &genome_population_simulator::technical_identity::TechnicalIdentityPlan,
    count: u64,
) {
    let library = identities.library(individual.id()).expect("library");
    let mut bytes = Vec::new();
    for ordinal in 0..count {
        bytes.extend_from_slice(b"@");
        bytes.extend_from_slice(
            format!("TPL_{}_{}", library.library_id().as_str(), ordinal).as_bytes(),
        );
        bytes.extend_from_slice(b"\n");
        bytes.extend(std::iter::repeat_n(b'N', 150));
        bytes.extend_from_slice(b"\n+\n");
        bytes.extend(std::iter::repeat_n(b'!', 150));
        bytes.extend_from_slice(b"\n");
    }
    write_gzip_fastq(path, &bytes);
}

fn write_gzip_fastq(path: &Path, bytes: &[u8]) {
    let file = fs::File::create(path).expect("FASTQ gzip fixture");
    let mut writer = bgzf::io::Writer::new(file);
    writer.write_all(bytes).expect("compress FASTQ fixture");
    writer.finish().expect("finish FASTQ gzip fixture");
}

fn read_gzip_fastq(path: &Path) -> Vec<u8> {
    let file = fs::File::open(path).expect("FASTQ gzip fixture");
    let mut reader = bgzf::io::Reader::new(file);
    let mut bytes = Vec::new();
    reader
        .read_to_end(&mut bytes)
        .expect("decompress FASTQ fixture");
    bytes
}

fn expected_artifact(
    path: &Path,
    identities: &genome_population_simulator::technical_identity::TechnicalIdentityPlan,
    individual: &Individual,
) -> genome_population_simulator::results::ArtifactRecord {
    let identity = identify_file(
        path,
        fastq_relative_path(identities.event_id().as_str(), individual.id().as_str()),
    )
    .expect("file identity");
    genome_population_simulator::results::ArtifactRecord::from_identity(
        format!(
            "fastq-{}",
            identities
                .library(individual.id())
                .expect("library")
                .library_id()
                .as_str()
        ),
        FASTQ_ARTIFACT_KIND,
        &identity,
    )
    .expect("artifact record")
}

#[test]
fn validation_replays_synthetic_truth_and_reports_coverage_without_a_tolerance() {
    let root = temporary("report");
    fs::create_dir_all(&root).expect("root");
    let reference = all_unknown_truth_reference();
    let individual = truth_state(&reference);
    let identities = identities(&reference, &individual);
    let target = EffectiveTechnicalEventTarget::CoverageWgs(2.0);
    let count = WgsLibrarySampler::new(&reference, &individual, &identities, target)
        .expect("sampler")
        .plan()
        .read_count();
    let fastq_path = root.join("independent.fastq.gz");
    write_independent_unknown_fastq(&fastq_path, &individual, &identities, count);
    let artifact = expected_artifact(&fastq_path, &identities, &individual);

    let report = validate_wgs_fastq(
        &fastq_path,
        &artifact,
        &reference,
        &individual,
        &identities,
        target,
    )
    .expect("valid FASTQ");
    assert_eq!(report.fastq_sha256(), artifact.digest);
    assert_eq!(report.verified_origin_count(), report.read_count());
    assert_eq!(report.base_count(), report.read_count() * 150);
    assert_eq!(report.unknown_base_count(), report.read_count() * 150);
    assert_eq!(report.coverage().target_coverage(), Some(2.0));
    assert!(report.coverage().absolute_deviation().is_some());
    assert!(report.coverage().relative_deviation().is_some());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn structural_mutations_of_id_origin_sequence_and_quality_fail_without_repairing_fastq() {
    let root = temporary("mutations");
    let reference = reference();
    let individual = truth_state(&reference);
    let identities = identities(&reference, &individual);
    let (fastq_path, artifact) = write_fastq(
        &root,
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(1),
    );
    let original = read_gzip_fastq(&fastq_path);
    let header_end = original
        .iter()
        .position(|base| *base == b'\n')
        .expect("header end");
    let sequence_start = header_end + 1;
    let sequence_end = sequence_start + 150;
    let quality_start = sequence_end + 3;

    let mut changed_id = original.clone();
    changed_id[header_end - 1] = b'9';
    let mut changed_sequence = original.clone();
    changed_sequence[sequence_start + 7] = match changed_sequence[sequence_start + 7] {
        b'A' => b'C',
        _ => b'A',
    };
    changed_sequence[quality_start + 7] = b'I';
    let mut changed_origin = original.clone();
    let expected_sequence = &original[sequence_start..sequence_end];
    let altered_origin_sequence = reverse_complement(expected_sequence);
    changed_origin[sequence_start..sequence_end].copy_from_slice(&altered_origin_sequence);
    for (base, quality) in altered_origin_sequence
        .iter()
        .zip(&mut changed_origin[quality_start..quality_start + 150])
    {
        *quality = if *base == b'N' { b'!' } else { b'I' };
    }
    let mut changed_quality = original.clone();
    changed_quality[quality_start] = match changed_quality[quality_start] {
        b'I' => b'!',
        b'!' => b'I',
        quality => panic!("unexpected literal quality byte: {quality:?}"),
    };
    let mut oversized_header = original.clone();
    oversized_header.splice(header_end..header_end, std::iter::repeat_n(b'X', 64 * 1024));

    for (name, bytes, expected_detail) in [
        (
            "id",
            changed_id,
            "FASTQ template ID does not match the reproduced origin",
        ),
        (
            "origin",
            changed_origin,
            "FASTQ sequence does not match the reproduced haplotype origin",
        ),
        (
            "sequence",
            changed_sequence,
            "FASTQ sequence does not match the reproduced haplotype origin",
        ),
        (
            "quality",
            changed_quality,
            "FASTQ Phred quality does not match its literal base",
        ),
        (
            "oversized-header",
            oversized_header,
            "FASTQ line exceeds the minimal profile bound",
        ),
    ] {
        write_gzip_fastq(&fastq_path, &bytes);
        let error = validate_wgs_fastq(
            &fastq_path,
            &artifact,
            &reference,
            &individual,
            &identities,
            EffectiveTechnicalEventTarget::ReadCount(1),
        )
        .expect_err(name);
        match error {
            FastqValidationError::Structural { detail, .. } => assert_eq!(detail, expected_detail),
            other => panic!("{name}: expected structural error, found {other:?}"),
        }
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn independently_frozen_haplotypic_origins_detect_copy_start_and_orientation_regressions() {
    // Fixture written by hand from the frozen v1 contract. Its non-periodic
    // reference makes the exact windows and their orientations observable.
    const FROZEN_FASTQ: &[u8] = b"@TPL_LIB_K3VL4L2HVRF7HOXF4TOCA6KR7E_0\n\
ATAGTATGCGTCCCCTCAAAGCCCTGAACCGACTATCAATCCCGCCAAGTANAAGCCATATACGGACGTTTGCGCGTATTCTTTTGGTGGCCGAATTCCGATGTGAATGTCATTCGAACGAAACTCTAGTACTTGCTTATAAAATCGATG\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@TPL_LIB_K3VL4L2HVRF7HOXF4TOCA6KR7E_1\n\
ACTACCTCGTCGCACCNCATGAGAGCGTCAAGGCGTACTTAATAGTATGCGTCGCCTCAAAGCCCTGAACCGACTATCCATCCCGCCAAGTANAAGCCATATACGGACGTTTGCGCGTATTCTTTTGGTGGCCGAATTCCGATGTGAATG\n\
+\n\
IIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@TPL_LIB_K3VL4L2HVRF7HOXF4TOCA6KR7E_2\n\
GTATGCGTCGCCTCAAAGCCCTGAACCGACTATCCATCCCGCCAAGTANAAGCCATATACGGACGTTTGCGCGTATTCTTTTGGTGGCCGAATTCCGATGTGAATGTCATTCGAACGAAACTCTAGTACTTGCTTATAAAATCGATGCCG\n\
+\n\
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let root = temporary("frozen");
    fs::create_dir_all(&root).expect("root");
    let reference = reference();
    let individual = differentiated_truth_state(&reference);
    let identities = identities(&reference, &individual);
    let origins = WgsLibrarySampler::new(
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(3),
    )
    .expect("sampler")
    .map(|read| read.expect("sampled read").interval().start())
    .collect::<Vec<_>>();
    assert_eq!(origins, vec![98, 57, 101]);
    assert!(origins.iter().any(|start| *start != 0));
    let fastq_path = root.join("frozen.fastq.gz");
    write_gzip_fastq(&fastq_path, FROZEN_FASTQ);
    let artifact = expected_artifact(&fastq_path, &identities, &individual);

    let report = validate_wgs_fastq(
        &fastq_path,
        &artifact,
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(3),
    )
    .expect("frozen origins validate");
    assert_eq!(report.read_count(), 3);
    assert_eq!(report.verified_origin_count(), 3);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn validation_rejects_missing_extra_truncated_and_malformed_records() {
    let root = temporary("record-contract");
    let reference = reference();
    let individual = truth_state(&reference);
    let identities = identities(&reference, &individual);
    let (fastq_path, artifact) = write_fastq(
        &root,
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(1),
    );
    let valid = read_gzip_fastq(&fastq_path);
    let header_end = valid
        .iter()
        .position(|base| *base == b'\n')
        .expect("header");
    let sequence_start = header_end + 1;
    let plus_start = sequence_start + 151;
    let quality_start = plus_start + 2;

    let mut malformed_separator = valid.clone();
    malformed_separator[plus_start] = b'-';
    let mut short_sequence = valid.clone();
    short_sequence.remove(sequence_start);
    let mut short_quality = valid.clone();
    short_quality.remove(quality_start);
    let mut extra_record = valid.clone();
    extra_record.extend_from_slice(&valid);

    for (name, bytes, ordinal, detail) in [
        ("missing", Vec::new(), 0, "FASTQ header line is truncated"),
        (
            "extra",
            extra_record,
            1,
            "FASTQ contains records beyond the reproduced read plan",
        ),
        (
            "truncated",
            valid[..valid.len() - 1].to_vec(),
            0,
            "FASTQ quality line is truncated",
        ),
        (
            "separator",
            malformed_separator,
            0,
            "FASTQ separator line is invalid",
        ),
        (
            "sequence-length",
            short_sequence,
            0,
            "FASTQ sequence length does not match the minimal profile",
        ),
        (
            "quality-length",
            short_quality,
            0,
            "FASTQ quality length does not match the minimal profile",
        ),
    ] {
        write_gzip_fastq(&fastq_path, &bytes);
        let error = validate_wgs_fastq(
            &fastq_path,
            &artifact,
            &reference,
            &individual,
            &identities,
            EffectiveTechnicalEventTarget::ReadCount(1),
        )
        .expect_err(name);
        match error {
            FastqValidationError::Structural {
                record_ordinal,
                detail: observed,
                ..
            } => {
                assert_eq!(record_ordinal, ordinal, "{name}");
                assert_eq!(observed, detail, "{name}");
            }
            other => panic!("{name}: expected structural error, found {other:?}"),
        }
    }
    write_gzip_fastq(&fastq_path, &valid);
    let error = validate_wgs_fastq(
        &fastq_path,
        &artifact,
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(2),
    )
    .expect_err("one record cannot satisfy a two-record plan");
    assert!(matches!(
        error,
        FastqValidationError::Structural {
            record_ordinal: 1,
            detail: "FASTQ header line is truncated",
            ..
        }
    ));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn validation_rejects_a_physically_truncated_bgzf_eof_marker() {
    let root = temporary("truncated-bgzf-eof");
    let reference = reference();
    let individual = truth_state(&reference);
    let identities = identities(&reference, &individual);
    let (fastq_path, artifact) = write_fastq(
        &root,
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(1),
    );
    let mut bytes = fs::read(&fastq_path).expect("compressed FASTQ");
    bytes.truncate(bytes.len() - 15);
    fs::write(&fastq_path, bytes).expect("truncated compressed FASTQ");

    let error = validate_wgs_fastq(
        &fastq_path,
        &artifact,
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(1),
    )
    .expect_err("a partial BGZF EOF marker must not be accepted as end of input");
    assert!(matches!(
        error,
        FastqValidationError::Structural {
            detail: "FASTQ BGZF stream is missing its canonical EOF marker",
            ..
        }
    ));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn validation_binds_exact_artifact_record_and_opening_errors() {
    let root = temporary("artifact-contract");
    let reference = reference();
    let individual = truth_state(&reference);
    let identities = identities(&reference, &individual);
    let (fastq_path, artifact) = write_fastq(
        &root,
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(1),
    );
    for (field, mutate) in [
        (
            "id",
            Box::new(
                |record: &mut genome_population_simulator::results::ArtifactRecord| {
                    record.id.push('x')
                },
            ) as Box<dyn Fn(&mut genome_population_simulator::results::ArtifactRecord)>,
        ),
        (
            "kind",
            Box::new(
                |record: &mut genome_population_simulator::results::ArtifactRecord| {
                    record.kind.push('x')
                },
            ),
        ),
        (
            "path",
            Box::new(
                |record: &mut genome_population_simulator::results::ArtifactRecord| {
                    record.path.push('x')
                },
            ),
        ),
        (
            "size_bytes",
            Box::new(
                |record: &mut genome_population_simulator::results::ArtifactRecord| {
                    record.size_bytes += 1
                },
            ),
        ),
        (
            "digest",
            Box::new(
                |record: &mut genome_population_simulator::results::ArtifactRecord| {
                    record.digest.replace_range(..1, "0")
                },
            ),
        ),
    ] {
        let mut invalid = artifact.clone();
        mutate(&mut invalid);
        let error = validate_wgs_fastq(
            &fastq_path,
            &invalid,
            &reference,
            &individual,
            &identities,
            EffectiveTechnicalEventTarget::ReadCount(1),
        )
        .expect_err(field);
        match error {
            FastqValidationError::Artifact {
                field: observed, ..
            } => assert_eq!(observed, field),
            other => panic!("{field}: expected artifact error, found {other:?}"),
        }
    }
    let error = validate_wgs_fastq(
        &root.join("missing.fastq.gz"),
        &artifact,
        &reference,
        &individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(1),
    )
    .expect_err("opening a missing FASTQ must fail");
    assert!(matches!(
        error,
        FastqValidationError::Io {
            operation: "open",
            ..
        }
    ));
    fs::remove_dir_all(root).expect("cleanup");
}
