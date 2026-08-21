use genome_population_simulator::randomness::{
    RandomSamplingError, RootSeed, RuntimeSubseedScope, categorical_f64, derive_subseed,
    fisher_yates, linear_uniform, log_uniform, open_unit_f64, poisson_count, uniform_u64_below,
    weighted_without_replacement,
};

fn root() -> RootSeed {
    RootSeed::from_bytes([7; 32])
}

#[test]
fn v4_scope_name_separates_each_previously_colliding_pair() {
    let cases = [
        (
            RuntimeSubseedScope::MatrixCandidateCount {
                scenario_id: "direct",
                reference_id: "r",
            },
            RuntimeSubseedScope::MatrixLocusSelection {
                scenario_id: "direct",
                reference_id: "r",
            },
        ),
        (
            RuntimeSubseedScope::MatrixFrequencyComponent {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                variant_id: "v",
            },
            RuntimeSubseedScope::MatrixFrequencyValue {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                variant_id: "v",
            },
        ),
        (
            RuntimeSubseedScope::IndelCount {
                scenario_id: "direct",
                reference_id: "r",
                individual_id: "i",
                sequence_id: "s",
                variant_class: "insertion",
                spatial_scope_id: "global",
            },
            RuntimeSubseedScope::IndelPosition {
                scenario_id: "direct",
                reference_id: "r",
                individual_id: "i",
                sequence_id: "s",
                variant_class: "insertion",
                spatial_scope_id: "global",
            },
        ),
        (
            RuntimeSubseedScope::InitialGenotype {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                individual_id: "i",
                variant_id: "v",
            },
            RuntimeSubseedScope::InitialPhase {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                individual_id: "i",
                variant_id: "v",
            },
        ),
    ];
    for (left, right) in cases {
        assert_ne!(
            derive_subseed(root(), left).unwrap(),
            derive_subseed(root(), right).unwrap()
        );
    }
}

#[test]
fn every_f2_runtime_scope_changes_when_its_production_components_change() {
    let cases = [
        (
            RuntimeSubseedScope::MatrixCandidateCount {
                scenario_id: "direct",
                reference_id: "r",
            },
            RuntimeSubseedScope::MatrixCandidateCount {
                scenario_id: "changed",
                reference_id: "r",
            },
        ),
        (
            RuntimeSubseedScope::MatrixLocusSelection {
                scenario_id: "direct",
                reference_id: "r",
            },
            RuntimeSubseedScope::MatrixLocusSelection {
                scenario_id: "changed",
                reference_id: "r",
            },
        ),
        (
            RuntimeSubseedScope::MatrixAlternateAllele {
                scenario_id: "direct",
                reference_id: "r",
                sequence_id: "s",
                position_1based: 1,
                reference_base: "A",
            },
            RuntimeSubseedScope::MatrixAlternateAllele {
                scenario_id: "direct",
                reference_id: "r",
                sequence_id: "s",
                position_1based: 2,
                reference_base: "A",
            },
        ),
        (
            RuntimeSubseedScope::MatrixExactAssignment {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
            },
            RuntimeSubseedScope::MatrixExactAssignment {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "changed",
            },
        ),
        (
            RuntimeSubseedScope::MatrixFrequencyComponent {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                variant_id: "v",
            },
            RuntimeSubseedScope::MatrixFrequencyComponent {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                variant_id: "changed",
            },
        ),
        (
            RuntimeSubseedScope::MatrixFrequencyValue {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                variant_id: "v",
            },
            RuntimeSubseedScope::MatrixFrequencyValue {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                variant_id: "changed",
            },
        ),
        (
            RuntimeSubseedScope::IndelCount {
                scenario_id: "direct",
                reference_id: "r",
                individual_id: "i",
                sequence_id: "s",
                variant_class: "insertion",
                spatial_scope_id: "global",
            },
            RuntimeSubseedScope::IndelCount {
                scenario_id: "direct",
                reference_id: "r",
                individual_id: "changed",
                sequence_id: "s",
                variant_class: "insertion",
                spatial_scope_id: "global",
            },
        ),
        (
            RuntimeSubseedScope::IndelPosition {
                scenario_id: "direct",
                reference_id: "r",
                individual_id: "i",
                sequence_id: "s",
                variant_class: "insertion",
                spatial_scope_id: "global",
            },
            RuntimeSubseedScope::IndelPosition {
                scenario_id: "direct",
                reference_id: "r",
                individual_id: "changed",
                sequence_id: "s",
                variant_class: "insertion",
                spatial_scope_id: "global",
            },
        ),
        (
            RuntimeSubseedScope::InsertionContent {
                scenario_id: "direct",
                reference_id: "r",
                individual_id: "i",
                sequence_id: "s",
                spatial_scope_id: "global",
                event_index: 0,
            },
            RuntimeSubseedScope::InsertionContent {
                scenario_id: "direct",
                reference_id: "r",
                individual_id: "i",
                sequence_id: "s",
                spatial_scope_id: "global",
                event_index: 1,
            },
        ),
        (
            RuntimeSubseedScope::InitialGenotype {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                individual_id: "i",
                variant_id: "v",
            },
            RuntimeSubseedScope::InitialGenotype {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                individual_id: "i",
                variant_id: "changed",
            },
        ),
        (
            RuntimeSubseedScope::InitialPhase {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                individual_id: "i",
                variant_id: "v",
            },
            RuntimeSubseedScope::InitialPhase {
                scenario_id: "direct",
                reference_id: "r",
                population_id: "p",
                individual_id: "i",
                variant_id: "changed",
            },
        ),
    ];
    for (base, mutated) in cases {
        assert_ne!(
            derive_subseed(root(), base).unwrap(),
            derive_subseed(root(), mutated).unwrap()
        );
    }
}

fn mutated_scope(scope: &str, field: &str) -> RuntimeSubseedScope<'static> {
    let changed = |name| if field == name { "changed" } else { name };
    match scope {
        "matrix_candidate_count" => RuntimeSubseedScope::MatrixCandidateCount {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
        },
        "matrix_locus_selection" => RuntimeSubseedScope::MatrixLocusSelection {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
        },
        "matrix_alternate_allele" => RuntimeSubseedScope::MatrixAlternateAllele {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            sequence_id: changed("sequence_id"),
            position_1based: if field == "position_1based" { 2 } else { 1 },
            reference_base: changed("ref"),
        },
        "matrix_exact_assignment" => RuntimeSubseedScope::MatrixExactAssignment {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            population_id: changed("population_id"),
        },
        "matrix_frequency_component" => RuntimeSubseedScope::MatrixFrequencyComponent {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            population_id: changed("population_id"),
            variant_id: changed("variant_id"),
        },
        "matrix_frequency_value" => RuntimeSubseedScope::MatrixFrequencyValue {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            population_id: changed("population_id"),
            variant_id: changed("variant_id"),
        },
        "indel_count" => RuntimeSubseedScope::IndelCount {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            individual_id: changed("individual_id"),
            sequence_id: changed("sequence_id"),
            variant_class: changed("variant_class"),
            spatial_scope_id: changed("spatial_scope_id"),
        },
        "indel_position" => RuntimeSubseedScope::IndelPosition {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            individual_id: changed("individual_id"),
            sequence_id: changed("sequence_id"),
            variant_class: changed("variant_class"),
            spatial_scope_id: changed("spatial_scope_id"),
        },
        "insertion_content" => RuntimeSubseedScope::InsertionContent {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            individual_id: changed("individual_id"),
            sequence_id: changed("sequence_id"),
            spatial_scope_id: changed("spatial_scope_id"),
            event_index: if field == "event_index" { 1 } else { 0 },
        },
        "initial_genotype" => RuntimeSubseedScope::InitialGenotype {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            population_id: changed("population_id"),
            individual_id: changed("individual_id"),
            variant_id: changed("variant_id"),
        },
        "initial_phase" => RuntimeSubseedScope::InitialPhase {
            scenario_id: changed("scenario_id"),
            reference_id: changed("reference_id"),
            population_id: changed("population_id"),
            individual_id: changed("individual_id"),
            variant_id: changed("variant_id"),
        },
        _ => unreachable!("closed F2 scope fixture"),
    }
}

#[test]
fn every_runtime_scope_component_mutation_reaches_production_derivation() {
    let fixtures: &[(&str, &[&str])] = &[
        ("matrix_candidate_count", &["scenario_id", "reference_id"]),
        ("matrix_locus_selection", &["scenario_id", "reference_id"]),
        (
            "matrix_alternate_allele",
            &[
                "scenario_id",
                "reference_id",
                "sequence_id",
                "position_1based",
                "ref",
            ],
        ),
        (
            "matrix_exact_assignment",
            &["scenario_id", "reference_id", "population_id"],
        ),
        (
            "matrix_frequency_component",
            &["scenario_id", "reference_id", "population_id", "variant_id"],
        ),
        (
            "matrix_frequency_value",
            &["scenario_id", "reference_id", "population_id", "variant_id"],
        ),
        (
            "indel_count",
            &[
                "scenario_id",
                "reference_id",
                "individual_id",
                "sequence_id",
                "variant_class",
                "spatial_scope_id",
            ],
        ),
        (
            "indel_position",
            &[
                "scenario_id",
                "reference_id",
                "individual_id",
                "sequence_id",
                "variant_class",
                "spatial_scope_id",
            ],
        ),
        (
            "insertion_content",
            &[
                "scenario_id",
                "reference_id",
                "individual_id",
                "sequence_id",
                "spatial_scope_id",
                "event_index",
            ],
        ),
        (
            "initial_genotype",
            &[
                "scenario_id",
                "reference_id",
                "population_id",
                "individual_id",
                "variant_id",
            ],
        ),
        (
            "initial_phase",
            &[
                "scenario_id",
                "reference_id",
                "population_id",
                "individual_id",
                "variant_id",
            ],
        ),
    ];
    for (scope, fields) in fixtures {
        for field in *fields {
            let base = mutated_scope(scope, "");
            let mutated = mutated_scope(scope, field);
            assert_ne!(
                derive_subseed(root(), base).unwrap(),
                derive_subseed(root(), mutated).unwrap(),
                "{scope}.{field}"
            );
        }
    }
}

#[test]
fn poisson_empirical_moments_and_small_pmf_match_an_independent_oracle() {
    let lambda = 1.0;
    let mut counts = [0_u64; 8];
    let mut values = Vec::with_capacity(4_096);
    for sample in 0_u64..4_096 {
        let mut bytes = [0_u8; 32];
        bytes[..8].copy_from_slice(&sample.to_be_bytes());
        let mut rng = derive_subseed(
            RootSeed::from_bytes(bytes),
            RuntimeSubseedScope::MatrixCandidateCount {
                scenario_id: "direct",
                reference_id: "poisson-oracle",
            },
        )
        .unwrap()
        .local_rng();
        let value = poisson_count(&mut rng, lambda).unwrap();
        if let Some(slot) = counts.get_mut(value as usize) {
            *slot += 1;
        }
        values.push(value as f64);
    }
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / n;
    assert!((mean - lambda).abs() < 0.08, "mean={mean}");
    assert!((variance - lambda).abs() < 0.12, "variance={variance}");
    let mut pmf = libm::exp(-lambda);
    for (k, observed) in counts.iter().take(6).enumerate() {
        assert!(((*observed as f64 / n) - pmf).abs() < 0.035, "k={k}");
        pmf *= lambda / (k + 1) as f64;
    }

    let lambda = 100.0;
    let mut values = Vec::with_capacity(4_096);
    let mut at_mean = 0_u64;
    let mut numeric_failures = 0_u64;
    for sample in 0_u64..4_096 {
        let mut bytes = [0_u8; 32];
        bytes[..8].copy_from_slice(&(sample + 4_096).to_be_bytes());
        let mut rng = derive_subseed(
            RootSeed::from_bytes(bytes),
            RuntimeSubseedScope::MatrixCandidateCount {
                scenario_id: "direct",
                reference_id: "poisson-ptrs-oracle",
            },
        )
        .unwrap()
        .local_rng();
        match poisson_count(&mut rng, lambda) {
            Ok(value) => {
                at_mean += u64::from(value == 100);
                values.push(value as f64);
            }
            Err(RandomSamplingError::NumericFailure) => numeric_failures += 1,
            Err(error) => panic!("unexpected PTRS error: {error:?}"),
        }
    }
    assert_eq!(numeric_failures, 27, "frozen 4096-seed corpus");
    let n = values.len() as f64;
    assert!(n > 4_000.0, "PTRS numeric failures must remain exceptional");
    let mean = values.iter().sum::<f64>() / n;
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / n;
    assert!((mean - lambda).abs() < 1.5, "PTRS mean={mean}");
    assert!((variance - lambda).abs() < 4.0, "PTRS variance={variance}");
    let mut pmf_at_mean = libm::exp(-lambda);
    for k in 1..=100 {
        pmf_at_mean *= lambda / k as f64;
    }
    assert!(
        (at_mean as f64 / n - pmf_at_mean).abs() < 0.012,
        "PTRS p100"
    );
}

#[test]
fn zero_consumers_preserve_the_local_stream() {
    let seed = derive_subseed(
        root(),
        RuntimeSubseedScope::MatrixCandidateCount {
            scenario_id: "direct",
            reference_id: "r",
        },
    )
    .unwrap();
    let mut observed = seed.local_rng();
    let mut untouched = seed.local_rng();
    assert_eq!(poisson_count(&mut observed, 0.0), Ok(0));
    assert_eq!(
        weighted_without_replacement(&mut observed, &[1.0], 0),
        Ok(Vec::new())
    );
    assert_eq!(observed.next_u64(), untouched.next_u64());
}

fn local_rng() -> genome_population_simulator::randomness::LocalRng {
    derive_subseed(
        root(),
        RuntimeSubseedScope::MatrixCandidateCount {
            scenario_id: "direct",
            reference_id: "rng-fixture",
        },
    )
    .expect("fixture subseed")
    .local_rng()
}

fn assert_unchanged_after(
    action: impl FnOnce(&mut genome_population_simulator::randomness::LocalRng),
) {
    let mut observed = local_rng();
    let mut untouched = local_rng();
    action(&mut observed);
    assert_eq!(observed.next_u64(), untouched.next_u64());
}

#[test]
fn primitive_rng_utilities_obey_preflight_ranges_and_stream_consumption() {
    let mut rng = local_rng();
    let open = open_unit_f64(&mut rng);
    assert!((0.0..1.0).contains(&open));
    assert_unchanged_after(|rng| {
        assert_eq!(
            uniform_u64_below(rng, 0),
            Err(RandomSamplingError::InvalidBound)
        )
    });
    assert_unchanged_after(|rng| assert_eq!(uniform_u64_below(rng, 1), Ok(0)));
    let mut rng = local_rng();
    assert!(uniform_u64_below(&mut rng, 7).expect("bounded") < 7);

    assert_unchanged_after(|rng| {
        assert_eq!(
            categorical_f64(rng, &[f64::NAN]),
            Err(RandomSamplingError::InvalidWeights)
        )
    });
    assert_unchanged_after(|rng| assert_eq!(categorical_f64(rng, &[0.0, 3.0, 0.0]), Ok(1)));
    let mut rng = local_rng();
    assert!(categorical_f64(&mut rng, &[1.0, 2.0, 3.0]).expect("categorical") < 3);

    assert_unchanged_after(|rng| {
        assert_eq!(
            weighted_without_replacement(rng, &[1.0], 2),
            Err(RandomSamplingError::InsufficientCandidates)
        )
    });
    let mut rng = local_rng();
    let selected = weighted_without_replacement(&mut rng, &[1.0, 2.0, 3.0], 2).expect("selection");
    assert_eq!(selected.len(), 2);
    assert_ne!(selected[0], selected[1]);

    assert_unchanged_after(|rng| {
        assert_eq!(
            linear_uniform(rng, 1.0, 1.0),
            Err(RandomSamplingError::InvalidInterval)
        )
    });
    assert_unchanged_after(|rng| {
        assert_eq!(
            log_uniform(rng, 0.0, 1.0),
            Err(RandomSamplingError::InvalidInterval)
        )
    });
    let mut rng = local_rng();
    assert!((0.0..1.0).contains(&linear_uniform(&mut rng, 0.0, 1.0).expect("linear")));
    assert!((0.1..0.9).contains(&log_uniform(&mut rng, 0.1, 0.9).expect("log")));
}

#[test]
fn shuffle_and_poisson_cover_contract_boundaries_without_hidden_consumption() {
    assert_unchanged_after(|rng| {
        let mut values = [1_u8];
        assert_eq!(fisher_yates(rng, &mut values), Ok(()));
    });
    let mut rng = local_rng();
    let mut values = [0_u8, 1, 2, 3];
    fisher_yates(&mut rng, &mut values).expect("shuffle");
    values.sort_unstable();
    assert_eq!(values, [0, 1, 2, 3]);

    for lambda in [
        -1.0,
        f64::NAN,
        f64::INFINITY,
        4_503_599_627_370_496.0_f64.next_up(),
    ] {
        assert_unchanged_after(|rng| {
            assert_eq!(
                poisson_count(rng, lambda),
                Err(RandomSamplingError::InvalidLambda)
            )
        });
    }
    assert_unchanged_after(|rng| assert_eq!(poisson_count(rng, 0.0), Ok(0)));
    for lambda in [f64::from_bits(1), 9.999_999, 10.0, 4_503_599_627_370_496.0] {
        let mut rng = local_rng();
        assert!(poisson_count(&mut rng, lambda).is_ok(), "lambda={lambda}");
    }
}
