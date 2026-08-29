//! Oráculos de F4.2.1 separados del productor: las fórmulas se calculan con
//! los literales del caso y no se reutilizan desde `WgsPairedEndPlan`.

use std::collections::BTreeMap;

use genome_population_simulator::coordinates::{GenomicInterval, InternalDeletion};
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{
    DeletionModification, ModificationProvenance, Segment, SegmentChain, SmallModification,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::technical_profile::resolve_f4_technical_profile;
use genome_population_simulator::wgs_sampling::{WgsPairedEndSampler, resolve_wgs_paired_end_plan};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn reference(sequences: &[(&str, &[u8])]) -> Reference {
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
        ReferenceId::parse("reference-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        sequences,
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
        EventId::parse("wgs-pe", "test").expect("event ID"),
        &[individual.id().clone()],
    )
    .expect("identities")
}

fn shorten_first_copy_to(
    individual: &Individual,
    reference: &Reference,
    sequence_id: &SequenceId,
    length: u64,
) -> Individual {
    let source_length = reference.sequences()[0].length();
    let removed = source_length
        .checked_sub(length)
        .expect("the fixture only shortens copies");
    let deletion = DeletionModification::new(
        InternalDeletion::new(
            GenomicInterval::new(
                reference.id().clone(),
                sequence_id.clone(),
                1,
                1 + removed,
                source_length,
            )
            .expect("deletion interval"),
        )
        .expect("internal deletion"),
        ModificationProvenance::new(EventId::parse("shorten-copy", "test").expect("event ID")),
    );
    let chain = SegmentChain::new(vec![
        Segment::new(
            GenomicInterval::new(
                reference.id().clone(),
                sequence_id.clone(),
                0,
                source_length,
                source_length,
            )
            .expect("full interval"),
            vec![SmallModification::Deletion(deletion)],
        )
        .expect("segment"),
    ])
    .expect("chain");
    individual
        .with_copy_segments(
            sequence_id,
            HaplotypePhase::from_copy_index(0).expect("first phase"),
            chain,
        )
        .expect("shortened first copy")
}

#[test]
fn paired_end_coverage_plan_uses_full_haploid_denominator_and_rounds_half_up() {
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let structural_reference = reference(&[
        ("chr-long", &vec![b'N'; 500]),
        ("chr-short", &vec![b'A'; 249]),
    ]);

    // Oráculo: el denominador estructural suma 500 + 249 una sola vez. El
    // numerador de una PE100 es 100 + 100, aun si los mates solapan.
    let denominator = 500_u64 + 249;
    let emitted_per_pair = 100_u64 + 100;
    let target = 2.0;
    let expected_pairs = (target * denominator as f64 / emitted_per_pair as f64).round() as u64;
    let plan =
        resolve_wgs_paired_end_plan(&structural_reference, profile, 100, target).expect("plan");
    assert_eq!(plan.haploid_reference_length(), denominator);
    assert_eq!(plan.emitted_bases_per_pair(), emitted_per_pair);
    assert_eq!(plan.pair_count(), expected_pairs);
    assert_eq!(plan.zero_coverage_sequences().len(), 1);
    assert_eq!(plan.zero_coverage_sequences()[0].as_str(), "chr-short");

    let tie_reference = reference(&[("chr1", &vec![b'A'; 300])]);
    let tie = resolve_wgs_paired_end_plan(&tie_reference, profile, 150, 0.5).expect("tie plan");
    // 0.5 * 300 / (150 + 150) = 0.5: el empate debe producir una pareja.
    assert_eq!(tie.pair_count(), 1);
    assert_eq!(tie.equivalent_coverage(), 1.0);
}

#[test]
fn coverage_sampler_preserves_copy_proportions_and_keeps_n_eligible() {
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let reference = reference(&[("chr-n", &vec![b'N'; 500]), ("chr-short", &vec![b'A'; 249])]);
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual ID"),
        &reference,
        Ploidy::new(2).expect("ploidy"),
    );
    let identities = identities(&reference, &individual);

    // Oráculo independiente: 600 * (500 + 249) / (150 + 150), redondeado a
    // parejas. La ploidía no entra en el denominador.
    let expected_pairs = (600.0_f64 * 749.0 / 300.0).round() as usize;
    let sampler = WgsPairedEndSampler::new_for_coverage(
        &reference,
        &individual,
        &identities,
        profile,
        150,
        600.0,
    )
    .expect("sampler");
    let plan = sampler.coverage_plan().expect("coverage plan");
    assert_eq!(plan.haploid_reference_length(), 749);
    assert_eq!(plan.pair_count(), expected_pairs as u64);

    let pairs = sampler.collect::<Result<Vec<_>, _>>().expect("pairs");
    assert_eq!(pairs.len(), expected_pairs);
    let mut by_copy = BTreeMap::new();
    for pair in &pairs {
        assert_eq!(
            pair.fragment().haplotype_copy_id().sequence_id().as_str(),
            "chr-n"
        );
        assert!(pair.r1().sequence().iter().all(|base| *base == b'N'));
        assert!(pair.r2().sequence().iter().all(|base| *base == b'N'));
        *by_copy
            .entry(pair.fragment().haplotype_copy_id().copy_index())
            .or_insert(0_usize) += 1;
    }
    let copy_0 = by_copy[&0];
    let copy_1 = by_copy[&1];
    assert!(
        copy_0 > 0 && copy_1 > 0,
        "las dos copias elegibles aportan pares"
    );
    assert!(
        copy_0.abs_diff(copy_1) < expected_pairs / 3,
        "el reparto por copia sigue el espacio uniforme de inicios"
    );
}

#[test]
fn sequences_shorter_than_the_minimum_insert_are_reported_zero_without_aborting() {
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let reference = reference(&[("chr-short", &vec![b'N'; 249])]);
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual ID"),
        &reference,
        Ploidy::new(1).expect("ploidy"),
    );
    let identities = identities(&reference, &individual);

    let sampler = WgsPairedEndSampler::new_for_coverage(
        &reference,
        &individual,
        &identities,
        profile,
        150,
        1.0,
    )
    .expect("the inaccessible structural sequence is not a fatal configuration error");
    let plan = sampler.coverage_plan().expect("coverage plan");
    assert_eq!(plan.zero_coverage_sequences()[0].as_str(), "chr-short");
    assert!(
        sampler
            .collect::<Result<Vec<_>, _>>()
            .expect("empty stream")
            .is_empty()
    );
}

#[test]
fn intermediate_length_conditions_insert_sampling_and_realizes_every_planned_pair() {
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let reference = reference(&[("chr-300", &vec![b'A'; 300])]);
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual ID"),
        &reference,
        Ploidy::new(1).expect("ploidy"),
    );
    let identities = identities(&reference, &individual);

    // Oráculo: con PE150, 100x sobre 300 pb exige exactamente 100 parejas.
    // El perfil permite hasta 450 pb, pero 301--450 no caben en esta única
    // copia: el sorteo debe condicionarse a 250..=300 y no abortar a mitad.
    let sampler = WgsPairedEndSampler::new_for_coverage(
        &reference,
        &individual,
        &identities,
        profile,
        150,
        100.0,
    )
    .expect("sampler");
    assert_eq!(sampler.coverage_plan().expect("plan").pair_count(), 100);
    assert_eq!(sampler.effective_insert_maximum_bases(), 300);
    let pairs = sampler.collect::<Result<Vec<_>, _>>().expect("all pairs");
    assert_eq!(pairs.len(), 100);
    assert!(
        pairs
            .iter()
            .all(|pair| (250..=300).contains(&(pair.fragment().sequence().len() as u64)))
    );
}

#[test]
fn mixed_length_copies_filter_candidates_after_each_conditioned_insert() {
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let reference = reference(&[("chr-mixed", &vec![b'A'; 450])]);
    let sequence_id = reference.sequences()[0].id().clone();
    let baseline = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual ID"),
        &reference,
        Ploidy::new(2).expect("ploidy"),
    );
    let individual = shorten_first_copy_to(&baseline, &reference, &sequence_id, 390);
    let identities = identities(&reference, &individual);

    let sampler = WgsPairedEndSampler::new_for_coverage(
        &reference,
        &individual,
        &identities,
        profile,
        150,
        100.0,
    )
    .expect("sampler");
    assert_eq!(sampler.coverage_plan().expect("plan").pair_count(), 150);
    assert_eq!(sampler.effective_insert_maximum_bases(), 450);
    let pairs = sampler.collect::<Result<Vec<_>, _>>().expect("all pairs");
    assert_eq!(pairs.len(), 150);

    let mut copies = BTreeMap::new();
    let mut long_insert_on_long_copy = false;
    for pair in pairs {
        let copy = pair.fragment().haplotype_copy_id().copy_index();
        let insert = pair.fragment().sequence().len() as u64;
        *copies.entry(copy).or_insert(0_usize) += 1;
        if copy == 0 {
            assert!(
                insert <= 390,
                "the shortened copy cannot host a longer insert"
            );
        }
        if insert > 390 {
            assert_eq!(copy, 1, "only the intact copy can host this insert");
            long_insert_on_long_copy = true;
        }
    }
    assert!(
        copies.contains_key(&0),
        "the shortened copy remains eligible"
    );
    assert!(copies.contains_key(&1), "the intact copy remains eligible");
    assert!(
        long_insert_on_long_copy,
        "the fixture exercises an incompatible insert"
    );
}

#[test]
fn f4_configuration_rejects_coverage_substitutions_before_sampling() {
    let toml = br#"
config_kind = "general"
schema_version = "2"

[reference]
id = "reference-main"
source = "simulated"
total_length = 500

[technical_event]
event_id = "wgs-pe"
individual_ids = ["sample-A"]
target_coverage_wgs = 1
coverage_by_individual = { sample-A = 2 }
"#;
    let path = std::env::temp_dir().join("gps-f4-coverage-substitution.toml");
    let admitted = genome_population_simulator::adapters::admit_toml_header(
        &path,
        toml.to_vec(),
        genome_population_simulator::configuration::ExpectedConfigKind::General,
    )
    .expect("F4 TOML header");
    assert!(genome_population_simulator::adapters::parse_declared_config(admitted).is_err());
}
