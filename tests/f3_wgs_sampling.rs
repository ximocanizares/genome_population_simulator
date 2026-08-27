use std::time::Instant;

use genome_population_simulator::application::EffectiveTechnicalEventTarget;
use genome_population_simulator::coordinates::{BasePosition, GenomicInterval};
use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{
    CanonicalBase, ModificationProvenance, Segment, SegmentChain, SnpModification,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::wgs_sampling::{
    ReadOrientation, WGS_COUNT_ALGORITHM, WGS_READ_LENGTH, WgsLibrarySampler, WgsSamplingError,
};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn reference(sequences: &[(&str, &[u8])]) -> Reference {
    reference_with_id("reference-main", sequences)
}

fn reference_with_id(id: &str, sequences: &[(&str, &[u8])]) -> Reference {
    let sequences = sequences
        .iter()
        .map(|(id, bases)| {
            ReferenceSequence::new(
                SequenceId::parse(*id, "test").expect("sequence ID"),
                FastaHeader::parse(id.as_bytes()).expect("FASTA header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(bases).expect("DNA"),
            )
            .expect("reference sequence")
        })
        .collect();
    Reference::new(
        ReferenceId::parse(id, "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        sequences,
    )
    .expect("reference")
}

fn individual(reference: &Reference, id: &str, ploidy: u8) -> Individual {
    Individual::from_reference(
        IndividualId::parse(id, "test").expect("individual ID"),
        reference,
        Ploidy::new(ploidy).expect("ploidy"),
    )
}

fn identities(
    reference: &Reference,
    individuals: &[&Individual],
) -> genome_population_simulator::technical_identity::TechnicalIdentityPlan {
    let ids = individuals
        .iter()
        .map(|individual| individual.id().clone())
        .collect::<Vec<_>>();
    derive_technical_identity_plan(
        RootSeed::parse("root", ROOT).expect("root seed"),
        reference.id(),
        EventId::parse("reads-main", "test").expect("event ID"),
        &ids,
    )
    .expect("identities")
}

fn reads(
    reference: &Reference,
    individual: &Individual,
    identities: &genome_population_simulator::technical_identity::TechnicalIdentityPlan,
    target: EffectiveTechnicalEventTarget,
) -> Vec<genome_population_simulator::wgs_sampling::SampledWgsRead> {
    WgsLibrarySampler::new(reference, individual, identities, target)
        .expect("sampler")
        .collect::<Result<Vec<_>, _>>()
        .expect("sampled reads")
}

#[test]
fn sampling_uses_haploid_eligible_denominator_and_exact_read_count() {
    let long = vec![b'A'; 300];
    let short = vec![b'C'; 149];
    let reference = reference(&[("chr-long", &long), ("chr-short", &short)]);
    let sample = individual(&reference, "sample-A", 4);
    let identities = identities(&reference, &[&sample]);
    let sampler = WgsLibrarySampler::new(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(7),
    )
    .expect("sampler");

    assert_eq!(sampler.plan().haploid_eligible_length(), 300);
    assert_eq!(sampler.plan().read_count(), 7);
    let reads = sampler.collect::<Result<Vec<_>, _>>().expect("reads");
    assert_eq!(reads.len(), 7);
    assert!(reads.iter().all(|read| {
        read.haplotype_copy_id().sequence_id().as_str() == "chr-long"
            && read.interval().end() - read.interval().start() == WGS_READ_LENGTH
            && read.interval().end() <= 300
    }));
    assert_eq!(
        reads
            .iter()
            .map(|read| read.template().template_ordinal().get())
            .collect::<Vec<_>>(),
        (0..7).collect::<Vec<_>>()
    );
}

#[test]
fn sampling_reads_the_current_phased_haplotype_without_truth_materialization() {
    let bases = vec![b'C'; 150];
    let reference = reference(&[("chr1", &bases)]);
    let initial = individual(&reference, "sample-A", 1);
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let snp = SnpModification::new(
        BasePosition::new(reference.id().clone(), sequence_id.clone(), 37, 150).expect("SNP"),
        CanonicalBase::C,
        CanonicalBase::A,
        ModificationProvenance::new(EventId::parse("edit-1", "test").expect("event ID")),
    )
    .expect("SNP modification");
    let chain = SegmentChain::new(vec![
        Segment::new(
            GenomicInterval::new(reference.id().clone(), sequence_id.clone(), 0, 150, 150)
                .expect("full interval"),
            vec![genome_population_simulator::segments::SmallModification::Snp(snp)],
        )
        .expect("segment"),
    ])
    .expect("chain");
    let sample = initial
        .with_copy_segments(
            &sequence_id,
            genome_population_simulator::haplotype::HaplotypePhase::from_copy_index(0)
                .expect("phase"),
            chain,
        )
        .expect("edited copy");
    let identities = identities(&reference, &[&sample]);
    let read = reads(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(1),
    )
    .pop()
    .expect("one read");

    match read.orientation() {
        ReadOrientation::Forward => {
            assert_eq!(read.true_sequence()[37], b'A');
            assert!(
                read.true_sequence()
                    .iter()
                    .enumerate()
                    .all(|(index, base)| index == 37 || *base == b'C')
            );
        }
        ReadOrientation::ReverseComplement => {
            assert_eq!(read.true_sequence()[112], b'T');
            assert!(
                read.true_sequence()
                    .iter()
                    .enumerate()
                    .all(|(index, base)| index == 112 || *base == b'G')
            );
        }
    }
}

#[test]
fn coverage_is_poisson_per_library_with_a_haploid_denominator() {
    let chr1 = vec![b'A'; 150];
    let chr2 = vec![b'C'; 300];
    let reference = reference(&[("chr1", &chr1), ("chr2", &chr2)]);
    let sample = individual(&reference, "sample-A", 3);
    let identities = identities(&reference, &[&sample]);
    let first = WgsLibrarySampler::new(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::CoverageWgs(4.0),
    )
    .expect("sampler");
    let repeated = WgsLibrarySampler::new(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::CoverageWgs(4.0),
    )
    .expect("same sampler");

    // 150 + 300, una vez cada cromosoma: no 3 × (150 + 300).
    assert_eq!(first.plan().haploid_eligible_length(), 450);
    assert_eq!(first.plan().read_count(), repeated.plan().read_count());
    assert_eq!(WGS_COUNT_ALGORITHM, "wgs-coverage-poisson-v1");
}

#[test]
fn copies_contribute_by_dose_and_eligible_starts() {
    let bases = vec![b'A'; 151];
    let reference = reference(&[("chr1", &bases)]);
    let sample = individual(&reference, "sample-A", 4);
    let identities = identities(&reference, &[&sample]);
    let reads = reads(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(4_000),
    );
    let mut by_copy = [0_usize; 4];
    for read in &reads {
        by_copy[usize::from(read.haplotype_copy_id().copy_index())] += 1;
        // Una secuencia de 151 bases tiene exactamente los inicios 0 y 1.
        assert!(read.interval().start() <= 1);
    }
    assert!(by_copy.iter().all(|count| (700..=1_300).contains(count)));
}

#[test]
fn unknown_bases_are_eligible_and_preserved_through_both_orientations() {
    let bases = vec![b'N'; 150];
    let reference = reference(&[("chr1", &bases)]);
    let sample = individual(&reference, "sample-A", 1);
    let identities = identities(&reference, &[&sample]);
    let reads = reads(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(32),
    );

    assert_eq!(reads.len(), 32);
    assert!(reads.iter().all(|read| read.true_sequence() == bases));
    assert!(
        reads
            .iter()
            .any(|read| read.orientation() == ReadOrientation::Forward)
    );
    assert!(
        reads
            .iter()
            .any(|read| read.orientation() == ReadOrientation::ReverseComplement)
    );
}

#[test]
fn short_references_have_no_eligible_interval() {
    let bases = vec![b'A'; 149];
    let reference = reference(&[("short", &bases)]);
    let sample = individual(&reference, "sample-A", 2);
    let identities = identities(&reference, &[&sample]);
    assert!(matches!(
        WgsLibrarySampler::new(
            &reference,
            &sample,
            &identities,
            EffectiveTechnicalEventTarget::ReadCount(1),
        ),
        Err(WgsSamplingError::NoEligibleSequence {
            read_length: WGS_READ_LENGTH,
        })
    ));
}

#[test]
fn zero_target_yields_an_empty_stream() {
    let bases = vec![b'A'; 150];
    let reference = reference(&[("chr1", &bases)]);
    let sample = individual(&reference, "sample-A", 1);
    let identities = identities(&reference, &[&sample]);
    let sampler = WgsLibrarySampler::new(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::CoverageWgs(0.0),
    )
    .expect("zero target is valid");
    assert_eq!(sampler.plan().read_count(), 0);
    assert_eq!(sampler.count(), 0);
}

#[test]
fn sampling_is_reproducible_and_invariant_to_individual_order_or_addition() {
    let bases = (0..151)
        .map(|index| match index % 5 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        })
        .collect::<Vec<_>>();
    let reference = reference(&[("chr1", &bases)]);
    let sample_a = individual(&reference, "sample-A", 2);
    let sample_b = individual(&reference, "sample-B", 2);
    let only_a = identities(&reference, &[&sample_a]);
    let reordered = identities(&reference, &[&sample_b, &sample_a]);
    let target = EffectiveTechnicalEventTarget::ReadCount(64);

    let baseline = reads(&reference, &sample_a, &only_a, target);
    assert_eq!(baseline, reads(&reference, &sample_a, &only_a, target));
    assert_eq!(baseline, reads(&reference, &sample_a, &reordered, target));

    let coverage_target = EffectiveTechnicalEventTarget::CoverageWgs(8.0);
    let baseline_coverage = reads(&reference, &sample_a, &only_a, coverage_target);
    assert_eq!(
        baseline_coverage,
        reads(&reference, &sample_a, &reordered, coverage_target)
    );
}

#[test]
fn sampler_rejects_a_plan_for_a_different_reference_before_sampling() {
    let bases = vec![b'A'; 150];
    let state_reference = reference_with_id("reference-state", &[("chr1", &bases)]);
    let plan_reference = reference_with_id("reference-plan", &[("chr1", &bases)]);
    let sample = individual(&state_reference, "sample-A", 1);
    let identities = identities(&plan_reference, &[&sample]);

    assert!(matches!(
        WgsLibrarySampler::new(
            &state_reference,
            &sample,
            &identities,
            EffectiveTechnicalEventTarget::ReadCount(1),
        ),
        Err(WgsSamplingError::TechnicalPlanReferenceMismatch)
    ));
}

#[test]
fn sampled_origin_vector_is_frozen_for_the_template_rng_order() {
    let bases = vec![b'A'; 151];
    let reference = reference(&[("chr1", &bases)]);
    let sample = individual(&reference, "sample-A", 2);
    let identities = identities(&reference, &[&sample]);
    let origins = reads(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(12),
    )
    .into_iter()
    .map(|read| {
        (
            read.haplotype_copy_id().copy_index(),
            read.interval().start(),
            read.orientation(),
        )
    })
    .collect::<Vec<_>>();

    // Vector de orígenes resuelto desde la verdad manual: dos copias intactas
    // de 151 bases, por lo que solo existen los inicios 0 y 1.
    assert_eq!(
        origins,
        vec![
            (1, 1, ReadOrientation::Forward),
            (0, 0, ReadOrientation::Forward),
            (0, 0, ReadOrientation::Forward),
            (1, 1, ReadOrientation::Forward),
            (1, 1, ReadOrientation::ReverseComplement),
            (0, 1, ReadOrientation::ReverseComplement),
            (1, 1, ReadOrientation::ReverseComplement),
            (0, 0, ReadOrientation::ReverseComplement),
            (1, 0, ReadOrientation::Forward),
            (0, 1, ReadOrientation::Forward),
            (1, 1, ReadOrientation::Forward),
            (0, 0, ReadOrientation::ReverseComplement),
        ]
    );
}

#[test]
fn sampling_many_ordinals_does_not_append_candidates() {
    let bases = vec![b'A'; 150];
    let reference = reference(&[("chr1", &bases)]);
    let sample = individual(&reference, "sample-A", 2);
    let identities = identities(&reference, &[&sample]);
    let mut sampler = WgsLibrarySampler::new(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(100_000),
    )
    .expect("sampler");
    let candidate_count = sampler.resident_candidate_count();
    for _ in 0..100_000 {
        let read = sampler.next().expect("read ordinal").expect("sampled read");
        assert_eq!(read.true_sequence().len(), 150);
    }
    assert_eq!(sampler.resident_candidate_count(), candidate_count);
    assert!(sampler.next().is_none());
}

#[test]
#[ignore = "experimento de rendimiento: ejecutar en release y con --nocapture"]
fn release_materialization_and_sampling_experiment_with_representative_snv_density() {
    const REFERENCE_LENGTH: usize = 1_000_000;
    const SNV_SPACING: usize = 1_000;
    const READ_COUNT: u64 = 10_000;

    let bases = vec![b'A'; REFERENCE_LENGTH];
    let reference = reference(&[("chr1", &bases)]);
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let modifications = (SNV_SPACING..REFERENCE_LENGTH)
        .step_by(SNV_SPACING)
        .map(|position| {
            genome_population_simulator::segments::SmallModification::Snp(
                SnpModification::new(
                    BasePosition::new(
                        reference.id().clone(),
                        sequence_id.clone(),
                        u64::try_from(position).expect("position fits"),
                        u64::try_from(REFERENCE_LENGTH).expect("length fits"),
                    )
                    .expect("SNP position"),
                    CanonicalBase::A,
                    CanonicalBase::C,
                    ModificationProvenance::new(
                        EventId::parse("benchmark-snp", "test").expect("event ID"),
                    ),
                )
                .expect("SNP"),
            )
        })
        .collect::<Vec<_>>();
    let chain = SegmentChain::new(vec![
        Segment::new(
            GenomicInterval::new(
                reference.id().clone(),
                sequence_id.clone(),
                0,
                u64::try_from(REFERENCE_LENGTH).expect("length fits"),
                u64::try_from(REFERENCE_LENGTH).expect("length fits"),
            )
            .expect("full interval"),
            modifications,
        )
        .expect("segment"),
    ])
    .expect("chain");
    let sample = individual(&reference, "sample-A", 2)
        .with_copy_segments(
            &sequence_id,
            genome_population_simulator::haplotype::HaplotypePhase::from_copy_index(0)
                .expect("phase"),
            chain,
        )
        .expect("edited copy");
    let identities = identities(&reference, &[&sample]);
    let copy = &sample.sequences()[0].copies()[0];

    let materialization_started = Instant::now();
    let materialized =
        genome_population_simulator::materialization::materialize_copy(copy, &reference)
            .expect("full materialization");
    let materialization_elapsed = materialization_started.elapsed();
    let sampling_started = Instant::now();
    let mut checksum = 0_u64;
    for sampled in WgsLibrarySampler::new(
        &reference,
        &sample,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(READ_COUNT),
    )
    .expect("sampler")
    {
        let read = sampled.expect("read");
        checksum = checksum.wrapping_add(u64::from(read.true_sequence()[0]));
    }
    let sampling_elapsed = sampling_started.elapsed();
    assert_eq!(materialized.len(), REFERENCE_LENGTH);
    assert_ne!(checksum, 0);
    println!(
        "scenario=reference_bases:{REFERENCE_LENGTH},snv_spacing:{SNV_SPACING},ploidy:2,reads:{READ_COUNT}; materialize_full_ms={:.3}; sampling_ms={:.3}; checksum={checksum}",
        materialization_elapsed.as_secs_f64() * 1_000.0,
        sampling_elapsed.as_secs_f64() * 1_000.0,
    );
}
