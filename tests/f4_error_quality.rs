//! Oráculos externos de F4.2.2: los valores esperados se calculan con los
//! literales DTPD-106 y no llaman al interpolador ni al productor.

use genome_population_simulator::error_quality::{
    ErrorQualityError, ErrorQualitySettings, SequencingErrorsMode, effective_phred_for_cycle,
    transform_paired_end,
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
fn disabled_mode_preserves_truth_and_profile_quality_with_zero_errors_per_mate() {
    let (root, event, pair, profile) = pair_from(&vec![b'A'; 250]);
    let disabled =
        ErrorQualitySettings::default().with_sequencing_errors(SequencingErrorsMode::Disabled);
    let transformed =
        transform_paired_end(root, &event, &pair, profile, disabled).expect("disabled transform");

    // Oráculo: la pareja fuente es la verdad antes de error/calidad. En el
    // modo neutro ambas secuencias deben conservarse literalmente, incluidos
    // los mates de orientación inversa, y las calidades son las del perfil.
    assert_eq!(transformed.observed().r1().sequence(), pair.r1().sequence());
    assert_eq!(transformed.observed().r2().sequence(), pair.r2().sequence());
    assert_eq!(transformed.observed().r1().phred(), &[36]);
    assert_eq!(transformed.observed().r2().phred(), &[35]);
    assert_eq!(
        format!(
            "@r1\n{}\n+\n{}\n",
            std::str::from_utf8(transformed.observed().r1().sequence()).expect("DNA"),
            std::str::from_utf8(
                &encode_fastq_sanger_qualities(transformed.observed().r1().phred())
                    .expect("Sanger"),
            )
            .expect("quality"),
        ),
        "@r1\nA\n+\nE\n"
    );
    assert_eq!(transformed.report().r1().opportunities(), 1);
    assert_eq!(transformed.report().r2().opportunities(), 1);
    assert_eq!(transformed.report().r1().errors(), 0);
    assert_eq!(transformed.report().r2().errors(), 0);
    assert!(
        !disabled
            .sequencing_errors()
            .qualities_represent_realized_error_probability()
    );

    let other_root = RootSeed::parse(
        "other root",
        "f00102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
    )
    .expect("other root");
    let with_other_seed = transform_paired_end(other_root, &event, &pair, profile, disabled)
        .expect("disabled transform");
    assert_eq!(with_other_seed, transformed);
}

#[test]
fn mate_error_streams_are_isolated_from_r1_adjustments() {
    let (root, event, pair, profile) = pair_from(&vec![b'A'; 250]);
    let lower_r1_scale = ErrorQualitySettings::new(-35, -35, 0.5, 1.0).expect("settings");
    let neutral_r1_scale = ErrorQualitySettings::new(-35, -35, 1.0, 1.0).expect("settings");
    let left = transform_paired_end(root, &event, &pair, profile, lower_r1_scale).expect("left");
    let right =
        transform_paired_end(root, &event, &pair, profile, neutral_r1_scale).expect("right");

    assert_eq!(left.observed().r2(), right.observed().r2());
    assert_eq!(left.report().r2(), right.report().r2());
}

#[test]
fn error_scale_recalculates_the_emitted_fastq_phred_from_final_probability() {
    let (root, event, pair, profile) = pair_from(&vec![b'A'; 250]);
    let doubled = ErrorQualitySettings::new(0, 0, 2.0, 2.0).expect("settings");
    let transformed =
        transform_paired_end(root, &event, &pair, profile, doubled).expect("transform");

    // Oráculo independiente: Q36 × 2 da
    // -10 log10(2 × 10^-3.6) = 32.9897..., que redondea a Q33.
    // Para R2, Q35 × 2 da 31.9897..., que redondea a Q32.
    assert_eq!(transformed.observed().r1().phred(), &[33]);
    assert_eq!(transformed.observed().r2().phred(), &[32]);
    assert_eq!(
        encode_fastq_sanger_qualities(transformed.observed().r1().phred()).expect("Sanger"),
        b"B"
    );
}

#[test]
fn final_probability_and_emitted_phred_boundaries_are_rejected_without_saturation() {
    let (root, event, pair, profile) = pair_from(&vec![b'A'; 250]);

    // Q36 con escala 10^-5.7 da exactamente Q93 tras la conversión inversa.
    let maximum_representable =
        ErrorQualitySettings::new(0, 0, 1.995_262_314_968_879_5e-6, 1.0).expect("settings");
    let transformed =
        transform_paired_end(root, &event, &pair, profile, maximum_representable).expect("Q93");
    assert_eq!(transformed.observed().r1().phred(), &[93]);
    assert_eq!(
        encode_fastq_sanger_qualities(transformed.observed().r1().phred()).expect("Sanger"),
        b"~"
    );

    let zero_probability = ErrorQualitySettings::new(0, 0, 0.0, 1.0).expect("settings");
    assert!(matches!(
        transform_paired_end(root, &event, &pair, profile, zero_probability),
        Err(ErrorQualityError::InvalidFinalProbability {
            mate: TechnicalMate::R1,
            cycle: 1,
            value,
        }) if value == 0.0
    ));

    let above_one_probability = ErrorQualitySettings::new(-36, -35, 2.0, 1.0).expect("settings");
    assert!(matches!(
        transform_paired_end(root, &event, &pair, profile, above_one_probability),
        Err(ErrorQualityError::InvalidFinalProbability {
            mate: TechnicalMate::R1,
            cycle: 1,
            value,
        }) if value == 2.0
    ));

    // Q36 con escala 10^-5.8 da una calidad final Q94, fuera del rango
    // FASTQ permitido; no se admite saturarla a Q93.
    let unrepresentable_phred =
        ErrorQualitySettings::new(0, 0, 1.584_893_192_461_114e-6, 1.0).expect("settings");
    assert!(matches!(
        transform_paired_end(root, &event, &pair, profile, unrepresentable_phred),
        Err(ErrorQualityError::InvalidFinalPhred {
            mate: TechnicalMate::R1,
            cycle: 1,
            value,
        }) if value == 94.0
    ));

    assert!(matches!(
        ErrorQualitySettings::new(0, 0, f64::INFINITY, 1.0),
        Err(ErrorQualityError::InvalidErrorScale {
            mate: TechnicalMate::R1,
            value,
        }) if value.is_infinite()
    ));
    assert!(encode_fastq_sanger_qualities(&[94]).is_err());
}
