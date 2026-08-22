use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::ReferenceId;
use genome_population_simulator::population_snvs::{
    FrequencyAssignmentMode, FrequencySampling, FrequencySpectrumComponent,
    GeneratedFrequencySpectrum, GeneratedLocusCount, GeneratedPopulationFrequency,
    GeneratedSpatialMode, GeneratedSubstitutionProfile, GeneratedWeightedRegion,
    PopulationSnvsGenerationConfig, PopulationSnvsGenerationError, generate_population_snvs,
};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_fasta::{SequenceExpectation, validate_supplied_fasta};
use genome_population_simulator::results::resolve_manifest_reference;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn root() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f2-2-5-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn reference(root: &Path) -> genome_population_simulator::results::ResolvedManifestReference {
    fs::create_dir_all(root).unwrap();
    fs::write(root.join("source.fa"), b">chr1\nANCGTN\n").unwrap();
    fs::write(
        root.join("config.toml"),
        r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "source.fa"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 6
"#,
    )
    .unwrap();
    let effective = load_and_validate_config(&root.join("config.toml")).unwrap();
    execute_validated_config(&effective, &root.join("published")).unwrap();
    resolve_manifest_reference(
        &root.join("published/manifest.json"),
        &ReferenceId::parse("ref-main", "test").unwrap(),
    )
    .unwrap()
}

fn config(count: GeneratedLocusCount) -> PopulationSnvsGenerationConfig {
    PopulationSnvsGenerationConfig {
        count,
        spatial: GeneratedSpatialMode::Uniform,
        substitution_profile: GeneratedSubstitutionProfile::UniformSnv12V1,
        populations: vec![
            GeneratedPopulationFrequency {
                population_id: "P1".into(),
                min_alt_af: 0.1,
                max_alt_af: 0.9,
                spectrum: None,
            },
            GeneratedPopulationFrequency {
                population_id: "P2".into(),
                min_alt_af: 0.2,
                max_alt_af: 0.8,
                spectrum: None,
            },
        ],
    }
}

#[test]
fn generation_excludes_n_is_repeatable_and_populations_are_isolated() {
    let directory = root();
    let reference = reference(&directory);
    let seed = RootSeed::from_bytes([7; 32]);
    let first =
        generate_population_snvs(&config(GeneratedLocusCount::Exact(4)), &reference, seed).unwrap();
    let second =
        generate_population_snvs(&config(GeneratedLocusCount::Exact(4)), &reference, seed).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.rows().len(), 8);
    assert!(first.rows().iter().all(|row| row.reference() != b'N'));
    assert!(
        first
            .rows()
            .iter()
            .all(|row| row.alternate() != row.reference())
    );
    assert!(
        first
            .rows()
            .iter()
            .all(|row| row.alt_af() > 0.0 && row.alt_af() < 1.0)
    );
    assert!(first.rows().iter().any(|row| {
        row.population_id().as_str() == "P1"
            && first.rows().iter().any(|other| {
                other.variant_id() == row.variant_id()
                    && other.population_id().as_str() == "P2"
                    && other.alt_af() != row.alt_af()
            })
    }));
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn generation_rejects_empty_or_insufficient_space_and_poisson_zero_is_empty() {
    let directory = root();
    let reference = reference(&directory);
    let seed = RootSeed::from_bytes([8; 32]);
    assert!(matches!(
        generate_population_snvs(&config(GeneratedLocusCount::Exact(5)), &reference, seed),
        Err(PopulationSnvsGenerationError::InsufficientEligibleSpace { .. })
    ));
    let zero = generate_population_snvs(
        &config(GeneratedLocusCount::PoissonMean(0.0)),
        &reference,
        seed,
    )
    .unwrap();
    assert!(zero.rows().is_empty());
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn poisson_positive_density_and_semiopen_frequency_ranges_are_observable() {
    let directory = root();
    let reference = reference(&directory);
    let seed = RootSeed::from_bytes([12; 32]);
    let positive = generate_population_snvs(
        &config(GeneratedLocusCount::PoissonMean(2.0)),
        &reference,
        seed,
    )
    .unwrap();
    assert!(
        positive
            .rows()
            .iter()
            .all(|row| row.alt_af() > 0.0 && row.alt_af() < 1.0)
    );
    let density = generate_population_snvs(
        &config(GeneratedLocusCount::PoissonDensityPerBase(0.5)),
        &reference,
        seed,
    )
    .unwrap();
    assert!(
        density
            .rows()
            .iter()
            .all(|row| row.alt_af() > 0.0 && row.alt_af() < 1.0)
    );
    let mut ranges = config(GeneratedLocusCount::Exact(2));
    ranges.populations[0].min_alt_af = 0.2;
    ranges.populations[0].max_alt_af = 0.3;
    let matrix = generate_population_snvs(&ranges, &reference, seed).unwrap();
    assert!(
        matrix
            .rows()
            .iter()
            .filter(|row| row.population_id().as_str() == "P1")
            .all(|row| row.alt_af() >= 0.2 && row.alt_af() < 0.3)
    );
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn population_frequency_streams_are_invariant_to_population_order_and_addition() {
    let directory = root();
    let reference = reference(&directory);
    let seed = RootSeed::from_bytes([13; 32]);
    let baseline =
        generate_population_snvs(&config(GeneratedLocusCount::Exact(3)), &reference, seed).unwrap();
    let mut altered = config(GeneratedLocusCount::Exact(3));
    altered.populations.reverse();
    altered.populations.push(GeneratedPopulationFrequency {
        population_id: "P3".into(),
        min_alt_af: 0.4,
        max_alt_af: 0.5,
        spectrum: None,
    });
    let result = generate_population_snvs(&altered, &reference, seed).unwrap();
    for row in baseline.rows() {
        let observed = result
            .rows()
            .iter()
            .find(|candidate| {
                candidate.variant_id() == row.variant_id()
                    && candidate.population_id() == row.population_id()
            })
            .unwrap();
        assert_eq!(observed.alt_af().to_bits(), row.alt_af().to_bits());
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn spectrum_exact_uses_largest_remainders_and_fixed_extremes() {
    let directory = root();
    let reference = reference(&directory);
    let mut parameters = config(GeneratedLocusCount::Exact(4));
    parameters.populations.truncate(1);
    parameters.populations[0].spectrum = Some(GeneratedFrequencySpectrum {
        mode: FrequencyAssignmentMode::Exact,
        components: vec![
            FrequencySpectrumComponent::Fixed {
                value: 0.0,
                weight: 1.0,
            },
            FrequencySpectrumComponent::Fixed {
                value: 1.0,
                weight: 1.0,
            },
            FrequencySpectrumComponent::Interval {
                min: 0.1,
                max: 0.2,
                weight: 2.0,
                sampling: FrequencySampling::LogUniform,
            },
        ],
    });
    let matrix =
        generate_population_snvs(&parameters, &reference, RootSeed::from_bytes([14; 32])).unwrap();
    let values = matrix
        .rows()
        .iter()
        .map(|row| row.alt_af())
        .collect::<Vec<_>>();
    assert_eq!(values.iter().filter(|value| **value == 0.0).count(), 1);
    assert_eq!(values.iter().filter(|value| **value == 1.0).count(), 1);
    assert_eq!(
        values
            .iter()
            .filter(|value| **value >= 0.1 && **value < 0.2)
            .count(),
        2
    );
    let mut continuous = config(GeneratedLocusCount::Exact(4));
    continuous.populations.truncate(1);
    continuous.populations[0].spectrum = Some(GeneratedFrequencySpectrum {
        mode: FrequencyAssignmentMode::Exact,
        components: vec![FrequencySpectrumComponent::Interval {
            min: 0.0,
            max: 1.0,
            weight: 1.0,
            sampling: FrequencySampling::LinearUniform,
        }],
    });
    let continuous =
        generate_population_snvs(&continuous, &reference, RootSeed::from_bytes([16; 32])).unwrap();
    assert!(
        continuous
            .rows()
            .iter()
            .all(|row| row.alt_af() > 0.0 && row.alt_af() < 1.0)
    );
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn rejects_overlapping_or_zero_weight_spectra() {
    let directory = root();
    let reference = reference(&directory);
    let mut parameters = config(GeneratedLocusCount::Exact(1));
    parameters.populations[0].spectrum = Some(GeneratedFrequencySpectrum {
        mode: FrequencyAssignmentMode::Stochastic,
        components: vec![
            FrequencySpectrumComponent::Interval {
                min: 0.1,
                max: 0.5,
                weight: 0.0,
                sampling: FrequencySampling::LinearUniform,
            },
            FrequencySpectrumComponent::Interval {
                min: 0.4,
                max: 0.6,
                weight: 0.0,
                sampling: FrequencySampling::LinearUniform,
            },
        ],
    });
    assert!(matches!(
        generate_population_snvs(&parameters, &reference, RootSeed::from_bytes([15; 32])),
        Err(PopulationSnvsGenerationError::InvalidFrequencySpectrum { .. })
    ));
    fs::remove_dir_all(directory).unwrap();
}

fn spectrum_validation_error(components: Vec<FrequencySpectrumComponent>) {
    let directory = root();
    let reference = reference(&directory);
    let mut parameters = config(GeneratedLocusCount::Exact(1));
    parameters.populations[0].spectrum = Some(GeneratedFrequencySpectrum {
        mode: FrequencyAssignmentMode::Stochastic,
        components,
    });
    assert!(matches!(
        generate_population_snvs(&parameters, &reference, RootSeed::from_bytes([18; 32])),
        Err(PopulationSnvsGenerationError::InvalidFrequencySpectrum { .. })
    ));
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn spectrum_rejects_overlapping_intervals_in_isolation() {
    spectrum_validation_error(vec![
        FrequencySpectrumComponent::Interval {
            min: 0.1,
            max: 0.6,
            weight: 1.0,
            sampling: FrequencySampling::LinearUniform,
        },
        FrequencySpectrumComponent::Interval {
            min: 0.5,
            max: 0.9,
            weight: 1.0,
            sampling: FrequencySampling::LinearUniform,
        },
    ]);
}

#[test]
fn spectrum_rejects_all_zero_negative_and_nonfinite_weights_in_isolation() {
    spectrum_validation_error(vec![FrequencySpectrumComponent::Fixed {
        value: 0.1,
        weight: 0.0,
    }]);
    spectrum_validation_error(vec![FrequencySpectrumComponent::Fixed {
        value: 0.1,
        weight: -1.0,
    }]);
    spectrum_validation_error(vec![FrequencySpectrumComponent::Fixed {
        value: 0.1,
        weight: f64::NAN,
    }]);
    spectrum_validation_error(vec![FrequencySpectrumComponent::Fixed {
        value: 0.1,
        weight: f64::INFINITY,
    }]);
}

#[test]
fn spectrum_rejects_invalid_log_and_fixed_components_in_isolation() {
    spectrum_validation_error(vec![FrequencySpectrumComponent::Interval {
        min: 0.0,
        max: 0.5,
        weight: 1.0,
        sampling: FrequencySampling::LogUniform,
    }]);
    spectrum_validation_error(vec![FrequencySpectrumComponent::Fixed {
        value: 1.1,
        weight: 1.0,
    }]);
    spectrum_validation_error(vec![
        FrequencySpectrumComponent::Fixed {
            value: 0.1,
            weight: 1.0,
        },
        FrequencySpectrumComponent::Fixed {
            value: 0.1,
            weight: 1.0,
        },
    ]);
    spectrum_validation_error(vec![
        FrequencySpectrumComponent::Fixed {
            value: 0.2,
            weight: 1.0,
        },
        FrequencySpectrumComponent::Interval {
            min: 0.1,
            max: 0.3,
            weight: 1.0,
            sampling: FrequencySampling::LinearUniform,
        },
    ]);
}

#[test]
fn spectrum_allows_adjacent_semiopen_intervals() {
    let directory = root();
    let reference = reference(&directory);
    let mut parameters = config(GeneratedLocusCount::Exact(4));
    parameters.populations.truncate(1);
    parameters.populations[0].spectrum = Some(GeneratedFrequencySpectrum {
        mode: FrequencyAssignmentMode::Exact,
        components: vec![
            FrequencySpectrumComponent::Interval {
                min: 0.0,
                max: 0.5,
                weight: 1.0,
                sampling: FrequencySampling::LinearUniform,
            },
            FrequencySpectrumComponent::Interval {
                min: 0.5,
                max: 1.0,
                weight: 1.0,
                sampling: FrequencySampling::LinearUniform,
            },
        ],
    });
    assert!(
        generate_population_snvs(&parameters, &reference, RootSeed::from_bytes([19; 32])).is_ok()
    );
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn exact_spectrum_uses_fractional_largest_remainders_and_stable_tie_order() {
    let directory = root();
    let reference = reference(&directory);
    let mut parameters = config(GeneratedLocusCount::Exact(4));
    parameters.populations.truncate(1);
    parameters.populations[0].spectrum = Some(GeneratedFrequencySpectrum {
        mode: FrequencyAssignmentMode::Exact,
        components: vec![
            FrequencySpectrumComponent::Fixed {
                value: 0.11,
                weight: 1.0,
            },
            FrequencySpectrumComponent::Fixed {
                value: 0.22,
                weight: 1.0,
            },
            FrequencySpectrumComponent::Fixed {
                value: 0.33,
                weight: 1.0,
            },
        ],
    });
    let matrix =
        generate_population_snvs(&parameters, &reference, RootSeed::from_bytes([17; 32])).unwrap();
    let values = matrix
        .rows()
        .iter()
        .map(|row| row.alt_af())
        .collect::<Vec<_>>();
    assert_eq!(values.iter().filter(|value| **value == 0.11).count(), 2);
    assert_eq!(values.iter().filter(|value| **value == 0.22).count(), 1);
    assert_eq!(values.iter().filter(|value| **value == 0.33).count(), 1);
    assert_eq!(values, vec![0.11, 0.33, 0.22, 0.11]);
    fs::remove_dir_all(directory).unwrap();
}

fn exact_population(
    id: &str,
    fixed: f64,
    min: f64,
    max: f64,
    include_zero_quota_component: bool,
) -> GeneratedPopulationFrequency {
    let mut components = vec![
        FrequencySpectrumComponent::Fixed {
            value: fixed,
            weight: 1.0,
        },
        FrequencySpectrumComponent::Interval {
            min,
            max,
            weight: 2.0,
            sampling: FrequencySampling::LinearUniform,
        },
    ];
    if include_zero_quota_component {
        // Para N=4 su cuota es cero: no debe alterar la asignación ni
        // disparar un sorteo de valor para esta componente fija.
        components.push(FrequencySpectrumComponent::Fixed {
            value: 0.99,
            weight: 0.0,
        });
    }
    GeneratedPopulationFrequency {
        population_id: id.to_owned(),
        min_alt_af: 0.0,
        max_alt_af: 1.0,
        spectrum: Some(GeneratedFrequencySpectrum {
            mode: FrequencyAssignmentMode::Exact,
            components,
        }),
    }
}

#[test]
fn exact_frequency_assignment_is_variant_id_stable_across_population_order_and_addition() {
    let directory = root();
    let reference = reference(&directory);
    let mut baseline = config(GeneratedLocusCount::Exact(4));
    baseline.populations = vec![
        exact_population("P1", 0.11, 0.40, 0.50, true),
        exact_population("P2", 0.22, 0.60, 0.70, true),
    ];
    let seed = RootSeed::from_bytes([21; 32]);
    let first = generate_population_snvs(&baseline, &reference, seed).unwrap();
    let mut reordered = baseline.clone();
    reordered.populations.reverse();
    reordered
        .populations
        .push(exact_population("P3", 0.33, 0.70, 0.80, true));
    let second = generate_population_snvs(&reordered, &reference, seed).unwrap();
    for row in first.rows() {
        let other = second
            .rows()
            .iter()
            .find(|candidate| {
                candidate.population_id() == row.population_id()
                    && candidate.variant_id() == row.variant_id()
            })
            .unwrap();
        assert_eq!(other.alt_af().to_bits(), row.alt_af().to_bits());
    }

    // Una componente de cuota cero no participa ni desplaza las corrientes de
    // asignación/valor de las componentes activas.
    let mut without_zero_quota = config(GeneratedLocusCount::Exact(4));
    without_zero_quota.populations = vec![
        exact_population("P1", 0.11, 0.40, 0.50, false),
        exact_population("P2", 0.22, 0.60, 0.70, false),
    ];
    let without_zero_quota =
        generate_population_snvs(&without_zero_quota, &reference, seed).unwrap();
    for row in first.rows() {
        let other = without_zero_quota
            .rows()
            .iter()
            .find(|candidate| {
                candidate.population_id() == row.population_id()
                    && candidate.variant_id() == row.variant_id()
            })
            .unwrap();
        assert_eq!(other.alt_af().to_bits(), row.alt_af().to_bits());
    }

    // Cambiar exclusivamente el valor fijo conserva las tres extracciones de
    // intervalo por VariantId: fixed no consume la corriente de valores.
    let mut changed_fixed = baseline.clone();
    changed_fixed.populations[0] = exact_population("P1", 0.12, 0.40, 0.50, true);
    let changed_fixed = generate_population_snvs(&changed_fixed, &reference, seed).unwrap();
    for row in first
        .rows()
        .iter()
        .filter(|row| row.population_id().as_str() == "P1")
    {
        let other = changed_fixed
            .rows()
            .iter()
            .find(|candidate| {
                candidate.population_id() == row.population_id()
                    && candidate.variant_id() == row.variant_id()
            })
            .unwrap();
        let expected = if row.alt_af() == 0.11 {
            0.12
        } else {
            row.alt_af()
        };
        assert_eq!(other.alt_af().to_bits(), expected.to_bits());
    }

    let p1 = first
        .rows()
        .iter()
        .filter(|row| row.population_id().as_str() == "P1")
        .map(|row| (row.variant_id().as_str(), row.alt_af()))
        .collect::<Vec<_>>();
    let p2 = first
        .rows()
        .iter()
        .filter(|row| row.population_id().as_str() == "P2")
        .map(|row| (row.variant_id().as_str(), row.alt_af()))
        .collect::<Vec<_>>();
    assert_eq!(
        p1,
        vec![
            ("VAR_USYUVYBQYFQOIUM2Z6LW4K7HKA", 0.11),
            ("VAR_SEF5VD7KSUS6KWUEZIK2ZXAIWI", 0.479_584_908_530_667_76),
            ("VAR_QJXGBAVGTGRCFFBQIJ2H5M5HCQ", 0.496_784_172_816_939_54),
            ("VAR_JE6MDR4M7NEHP6GJ4WCUEEDK7E", 0.472_674_930_901_828_16),
        ]
    );
    assert_eq!(
        p2,
        vec![
            ("VAR_USYUVYBQYFQOIUM2Z6LW4K7HKA", 0.682_954_222_217_499_1),
            ("VAR_SEF5VD7KSUS6KWUEZIK2ZXAIWI", 0.22),
            ("VAR_QJXGBAVGTGRCFFBQIJ2H5M5HCQ", 0.642_199_956_246_785_5),
            ("VAR_JE6MDR4M7NEHP6GJ4WCUEEDK7E", 0.619_627_744_784_811_5),
        ]
    );
    // Con las mismas coordenadas, las cuotas fijas recaen en VariantId
    // distintos entre P1 y P2; por tanto no siguen el orden canónico.
    assert_ne!(p1[0].0, p2[1].0);
    for (population, fixed, min, max) in [("P1", 0.11, 0.40, 0.50), ("P2", 0.22, 0.60, 0.70)] {
        let values = first
            .rows()
            .iter()
            .filter(|row| row.population_id().as_str() == population)
            .map(|row| row.alt_af())
            .collect::<Vec<_>>();
        assert_eq!(values.iter().filter(|value| **value == fixed).count(), 1);
        assert_eq!(
            values
                .iter()
                .filter(|value| **value >= min && **value < max)
                .count(),
            3
        );
        assert!(!values.contains(&0.99));
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn weighted_regions_select_only_positive_bases() {
    let directory = root();
    let reference = reference(&directory);
    let mut weighted = config(GeneratedLocusCount::Exact(1));
    weighted.spatial = GeneratedSpatialMode::WeightedRegions {
        outside_weight: 0.0,
        regions: vec![GeneratedWeightedRegion {
            contig: "chr1".into(),
            start_1based: 3,
            end_1based: 3,
            weight: 1.0,
        }],
    };
    let matrix =
        generate_population_snvs(&weighted, &reference, RootSeed::from_bytes([9; 32])).unwrap();
    assert!(matrix.rows().iter().all(|row| row.position_1based() == 3));
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn weighted_regions_reject_unknown_contigs_and_out_of_reference_coordinates() {
    let directory = root();
    let reference = reference(&directory);
    for (contig, end) in [("unknown", 1), ("chr1", 7)] {
        let mut invalid = config(GeneratedLocusCount::Exact(1));
        invalid.spatial = GeneratedSpatialMode::WeightedRegions {
            outside_weight: 0.0,
            regions: vec![GeneratedWeightedRegion {
                contig: contig.into(),
                start_1based: 1,
                end_1based: end,
                weight: 1.0,
            }],
        };
        assert!(matches!(
            generate_population_snvs(&invalid, &reference, RootSeed::from_bytes([11; 32])),
            Err(PopulationSnvsGenerationError::InvalidRegion { .. })
        ));
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn resolved_and_direct_validated_reference_share_the_same_generation_contract() {
    let directory = root();
    let resolved = reference(&directory);
    let direct = validate_supplied_fasta(
        &directory.join("source.fa"),
        "ref-main",
        &[SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 6,
        }],
    )
    .unwrap()
    .into_reference();
    let parameters = config(GeneratedLocusCount::Exact(3));
    let seed = RootSeed::from_bytes([10; 32]);
    assert_eq!(
        generate_population_snvs(&parameters, &resolved, seed).unwrap(),
        genome_population_simulator::population_snvs::generate_population_snvs_from_reference(
            &parameters,
            &direct,
            seed
        )
        .unwrap(),
    );
    fs::remove_dir_all(directory).unwrap();
}
