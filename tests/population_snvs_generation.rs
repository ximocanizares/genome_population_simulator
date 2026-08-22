use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::ReferenceId;
use genome_population_simulator::population_snvs::{
    GeneratedLocusCount, GeneratedPopulationFrequency, GeneratedSpatialMode,
    GeneratedWeightedRegion, PopulationSnvsGenerationConfig, PopulationSnvsGenerationError,
    generate_population_snvs,
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
        populations: vec![
            GeneratedPopulationFrequency {
                population_id: "P1".into(),
                min_alt_af: 0.1,
                max_alt_af: 0.9,
            },
            GeneratedPopulationFrequency {
                population_id: "P2".into(),
                min_alt_af: 0.2,
                max_alt_af: 0.8,
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
