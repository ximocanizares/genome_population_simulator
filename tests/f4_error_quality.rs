//! Oráculos externos de F4.2.2: los valores esperados se calculan con los
//! literales DTPD-106 y no llaman al interpolador ni al productor.

use genome_population_simulator::error_quality::{
    ErrorQualitySettings, effective_phred_for_cycle, transform_paired_end,
};
use genome_population_simulator::fastq_output::encode_fastq_sanger_qualities;
use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::technical_identity::{
    TechnicalMate, derive_technical_identity_plan,
};
use genome_population_simulator::technical_profile::resolve_f4_technical_profile;
use genome_population_simulator::wgs_sampling::WgsPairedEndSampler;

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn pair_from(
    bases: &[u8],
) -> (
    RootSeed,
    EventId,
    genome_population_simulator::paired_end::PairedEndPair,
    genome_population_simulator::technical_profile::F4TechnicalProfile,
) {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "test").expect("sequence"),
        FastaHeader::parse(b"chr1").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(bases).expect("DNA"),
    )
    .expect("reference sequence");
    let reference = Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference");
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual"),
        &reference,
        Ploidy::new(1).expect("ploidy"),
    );
    let root = RootSeed::parse("root", ROOT).expect("root");
    let event = EventId::parse("error-quality", "test").expect("event");
    let identities = derive_technical_identity_plan(
        root,
        reference.id(),
        event.clone(),
        &[individual.id().clone()],
    )
    .expect("identities");
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let pair = WgsPairedEndSampler::new_for_pair_count(
        &reference,
        &individual,
        &identities,
        profile,
        1,
        1,
    )
    .expect("sampler")
    .next()
    .expect("one pair")
    .expect("sampled pair");
    (root, event, pair, profile)
}

#[test]
fn quality_is_pure_piecewise_half_up_and_needs_no_rng_source() {
    let (_, _, _, profile) = pair_from(&vec![b'A'; 250]);
    let neutral = ErrorQualitySettings::default();

    // Oráculo independiente: los tramos son 36→36, 36→32 y 32→28 en R1;
    // y 35→35, 35→30, 30→25 en R2. 55 en R1 da
    // round_half_up(36 - 4 * 44 / 89) = 34.
    assert_eq!(
        effective_phred_for_cycle(profile, TechnicalMate::R1, 10, neutral),
        Ok(36)
    );
    assert_eq!(
        effective_phred_for_cycle(profile, TechnicalMate::R1, 55, neutral),
        Ok(34)
    );
    assert_eq!(
        effective_phred_for_cycle(profile, TechnicalMate::R1, 150, neutral),
        Ok(28)
    );
    assert_eq!(
        effective_phred_for_cycle(profile, TechnicalMate::R2, 55, neutral),
        Ok(33)
    );
    assert_eq!(
        effective_phred_for_cycle(profile, TechnicalMate::R2, 150, neutral),
        Ok(25)
    );

    // La API de calidad recibe solo perfil, mate, ciclo y ajustes: no acepta
    // RNG. Por ello repetirla no puede avanzar ningún flujo de error.
    let shifted = ErrorQualitySettings::new(-36, -35, 1.0, 1.0).expect("settings");
    assert_eq!(
        effective_phred_for_cycle(profile, TechnicalMate::R1, 1, shifted),
        Ok(0)
    );
    assert_eq!(
        effective_phred_for_cycle(profile, TechnicalMate::R2, 1, shifted),
        Ok(0)
    );
}

#[test]
fn forced_errors_substitute_only_other_acgt_and_publish_mate_aggregates() {
    let (root, event, pair, profile) = pair_from(&vec![b'A'; 250]);
    // Q=0 implica p=10^0×1=1, de modo que el oráculo no depende de suerte.
    let settings = ErrorQualitySettings::new(-36, -35, 1.0, 1.0).expect("settings");
    let transformed =
        transform_paired_end(root, &event, &pair, profile, settings).expect("transform");
    let observed = transformed.observed();

    assert_eq!(transformed.report().r1().opportunities(), 1);
    assert_eq!(transformed.report().r1().errors(), 1);
    assert_eq!(transformed.report().r1().realized_rate(), 1.0);
    assert_eq!(transformed.report().r2().opportunities(), 1);
    assert_eq!(transformed.report().r2().errors(), 1);
    assert_eq!(transformed.report().r2().realized_rate(), 1.0);
    assert_ne!(observed.r1().sequence(), pair.r1().sequence());
    assert_ne!(observed.r2().sequence(), pair.r2().sequence());
    assert!(matches!(observed.r1().sequence()[0], b'C' | b'G' | b'T'));
    assert!(matches!(observed.r2().sequence()[0], b'A' | b'C' | b'G'));
    assert_eq!(observed.r1().phred(), &[0]);
    assert_eq!(observed.r2().phred(), &[0]);
    assert_eq!(
        encode_fastq_sanger_qualities(observed.r1().phred()).expect("Sanger"),
        b"!"
    );
}

#[test]
fn n_is_preserved_at_phred_zero_without_an_error_opportunity() {
    let (root, event, pair, profile) = pair_from(&vec![b'N'; 250]);
    let settings = ErrorQualitySettings::new(-36, -35, 1.0, 1.0).expect("settings");
    let transformed =
        transform_paired_end(root, &event, &pair, profile, settings).expect("transform");

    assert_eq!(transformed.observed().r1().sequence(), b"N");
    assert_eq!(transformed.observed().r2().sequence(), b"N");
    assert_eq!(transformed.observed().r1().phred(), &[0]);
    assert_eq!(transformed.observed().r2().phred(), &[0]);
    assert_eq!(transformed.report().r1().opportunities(), 0);
    assert_eq!(transformed.report().r1().errors(), 0);
    assert_eq!(transformed.report().r2().opportunities(), 0);
    assert_eq!(transformed.report().r2().errors(), 0);
}

#[test]
fn mate_error_streams_are_isolated_from_r1_adjustments() {
    let (root, event, pair, profile) = pair_from(&vec![b'A'; 250]);
    let r1_disabled = ErrorQualitySettings::new(-36, -35, 0.0, 1.0).expect("settings");
    let both_enabled = ErrorQualitySettings::new(-36, -35, 1.0, 1.0).expect("settings");
    let left = transform_paired_end(root, &event, &pair, profile, r1_disabled).expect("left");
    let right = transform_paired_end(root, &event, &pair, profile, both_enabled).expect("right");

    assert_eq!(left.observed().r2(), right.observed().r2());
    assert_eq!(left.report().r2(), right.report().r2());
    assert_eq!(left.report().r1().errors(), 0);
    assert_eq!(right.report().r1().errors(), 1);
}

#[test]
fn invalid_probability_and_fastq_phred_are_rejected_without_saturation() {
    let (root, event, pair, profile) = pair_from(&vec![b'A'; 250]);
    let invalid_probability = ErrorQualitySettings::new(-36, -35, 2.0, 1.0).expect("settings");
    assert!(transform_paired_end(root, &event, &pair, profile, invalid_probability).is_err());
    assert!(encode_fastq_sanger_qualities(&[94]).is_err());
}
