use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{EffectiveConfig, build_effective_config};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::initial_genotypes::{
    InitialGenotypeSamplingError, sample_initial_genotypes,
};
use genome_population_simulator::population_snvs::{
    PreparedPopulationSnvs, PreparedReferenceIdentity, TSV_HEADER, parse_canonical_tsv,
    prepare_population_snvs_from_reference,
};
use genome_population_simulator::randomness::{
    RandomSamplingError, RandomWordSource, RootSeed, RuntimeSubseedScope, derive_subseed,
    f_diploid_genotype, f_diploid_genotype_with_source, hwe_composition,
    hwe_composition_with_source,
};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};

fn build(source: &str) -> EffectiveConfig {
    let admitted = admit_toml_header(
        Path::new("initial-genotypes.toml"),
        source.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("valid TOML header");
    build_effective_config(&parse_declared_config(admitted).expect("declared config"))
        .expect("effective config")
}

fn config(populations: &str) -> EffectiveConfig {
    config_with_reference("ref-main", populations)
}

fn config_with_reference(reference_id: &str, populations: &str) -> EffectiveConfig {
    build(&format!(
        r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
{populations}
[reference]
kind = "simulated"
id = "{reference_id}"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 8
"#
    ))
}

struct ScriptedWords {
    words: Vec<u64>,
    consumed: usize,
}

impl ScriptedWords {
    fn new(words: &[u64]) -> Self {
        Self {
            words: words.to_vec(),
            consumed: 0,
        }
    }
}

impl RandomWordSource for ScriptedWords {
    fn next_word(&mut self) -> u64 {
        let word = self.words[self.consumed];
        self.consumed += 1;
        word
    }
}

fn cohort(id: &str, ploidy: u8, model: &str, extra: &str) -> String {
    format!(
        r#"[[initial_populations]]
population_id = "{id}"
individual_count = 1
individual_ids = ["I-{id}"]
ploidy = {ploidy}
genotype_model = "{model}"
{extra}
"#
    )
}

fn reference() -> Reference {
    let id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let dna = DnaSequence::try_from_bytes(b"ACGTACGT").expect("DNA");
    let sequence = ReferenceSequence::new(
        id,
        FastaHeader::parse(b"chr1").expect("header"),
        SequenceClass::Chromosome,
        &dna,
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn reference_with_distinct_tail() -> Reference {
    let id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let dna = DnaSequence::try_from_bytes(b"ACGTACGA").expect("DNA");
    let sequence = ReferenceSequence::new(
        id,
        FastaHeader::parse(b"chr1").expect("header"),
        SequenceClass::Chromosome,
        &dna,
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn reference_with_nonlexicographic_contigs() -> Reference {
    let sequences = ["chr2", "chr1"]
        .into_iter()
        .map(|name| {
            let id = SequenceId::parse(name, "test").expect("sequence ID");
            let dna = DnaSequence::try_from_bytes(b"ACGT").expect("DNA");
            ReferenceSequence::new(
                id,
                FastaHeader::parse(name.as_bytes()).expect("header"),
                SequenceClass::Chromosome,
                &dna,
            )
            .expect("reference sequence")
        })
        .collect();
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        sequences,
    )
    .expect("reference")
}

fn matrix(reference: &Reference, rows: &[&str]) -> PreparedPopulationSnvs {
    let mut text = TSV_HEADER.join("\t");
    text.push('\n');
    text.push_str(&rows.join("\n"));
    text.push('\n');
    let table = parse_canonical_tsv(text.as_bytes()).expect("canonical matrix");
    prepare_population_snvs_from_reference(&table, reference, 1e-9, false)
        .expect("prepared population SNVs")
}

fn alternate_dosage(
    value: &genome_population_simulator::initial_genotypes::InitialLocusComposition,
    allele: u8,
) -> u8 {
    value
        .alternate_dosages()
        .iter()
        .find(|dose| dose.alternate() == allele)
        .expect("ALT present")
        .dosage()
}

#[test]
fn sampling_requires_the_reference_bound_to_the_prepared_matrix() {
    let reference = reference();
    let prepared = matrix(&reference, &["VAR_IGNORED\tchr1\t1\tA\tC\tP1\t0.25\ttest"]);
    assert_eq!(prepared.reference_id().as_str(), "ref-main");
    assert!(
        !prepared.reference_fingerprint().is_empty(),
        "preparation retains the validated canonical reference fingerprint"
    );

    let populations = cohort("P1", 2, "hwe", "");
    assert!(
        sample_initial_genotypes(
            &config_with_reference("ref-main", &populations),
            &reference,
            &prepared,
        )
        .is_ok(),
        "the configuration may consume the matrix bound to ref-main"
    );

    let error = sample_initial_genotypes(
        &config_with_reference("ref-other", &populations),
        &reference,
        &prepared,
    )
    .expect_err("a differently identified reference must fail before genotype RNG");
    assert_eq!(
        error,
        InitialGenotypeSamplingError::ReferenceMismatch {
            configured_reference_id: "ref-other".to_owned(),
            sampled_reference: PreparedReferenceIdentity::from_reference(&reference),
            prepared_reference: prepared.prepared_reference().clone(),
        }
    );
}

#[test]
fn sampling_rejects_same_id_with_a_different_reference_fingerprint_before_rng() {
    let prepared_reference = reference();
    let prepared = matrix(
        &prepared_reference,
        &["VAR_IGNORED\tchr1\t1\tA\tC\tP1\t0.25\ttest"],
    );
    let sampled_reference = reference_with_distinct_tail();
    assert_eq!(prepared_reference.id(), sampled_reference.id());
    assert_ne!(
        PreparedReferenceIdentity::from_reference(&prepared_reference),
        PreparedReferenceIdentity::from_reference(&sampled_reference),
        "a base outside the sampled locus still changes the canonical reference fingerprint"
    );

    let error = sample_initial_genotypes(
        &config_with_reference("ref-main", &cohort("P1", 2, "hwe", "")),
        &sampled_reference,
        &prepared,
    )
    .expect_err("the prepared matrix must reject a distinct reference before genotype RNG");
    assert_eq!(
        error,
        InitialGenotypeSamplingError::ReferenceMismatch {
            configured_reference_id: "ref-main".to_owned(),
            sampled_reference: PreparedReferenceIdentity::from_reference(&sampled_reference),
            prepared_reference: prepared.prepared_reference().clone(),
        }
    );
}

#[test]
fn hwe_multiallelic_normative_word_vectors_cover_ploidies_and_consumption() {
    const SENTINEL: u64 = 0xdeaf_beef_cafe_babe;
    const WORDS: [u64; 10] = [
        0x0000_0000_0000_0000,
        0x4000_0000_0000_0000,
        0x8000_0000_0000_0000,
        0xc000_0000_0000_0000,
        0x0000_0000_0000_0000,
        0x4000_0000_0000_0000,
        0x8000_0000_0000_0000,
        0xc000_0000_0000_0000,
        0x0000_0000_0000_0000,
        0x4000_0000_0000_0000,
    ];
    const EXPECTED: [[u8; 3]; 10] = [
        [1, 0, 0],
        [1, 1, 0],
        [1, 1, 1],
        [1, 1, 2],
        [2, 1, 2],
        [2, 2, 2],
        [2, 2, 3],
        [2, 2, 4],
        [3, 2, 4],
        [3, 3, 4],
    ];
    let frequencies = [0.2, 0.3, 0.5];

    for (ploidy, expected) in (1_u8..=10).zip(EXPECTED) {
        let mut input = WORDS[..usize::from(ploidy)].to_vec();
        input.push(SENTINEL);
        let mut source = ScriptedWords::new(&input);
        let observed = hwe_composition_with_source(&mut source, &frequencies, ploidy)
            .expect("literal HWE vector");
        assert_eq!(observed, expected, "P={ploidy}");
        assert_eq!(
            observed,
            independent_multinomial(&frequencies, &WORDS[..usize::from(ploidy)]),
            "independent multinomial oracle for P={ploidy}"
        );
        assert_eq!(source.next_word(), SENTINEL, "P={ploidy} consumption");
    }
}

#[test]
fn f_multiallelic_normative_word_vectors_cover_f_categories_and_consumption() {
    const SENTINEL: u64 = 0xdeaf_beef_cafe_babe;
    let frequencies = [0.2, 0.3, 0.5];
    let vectors = [
        (
            0.0,
            [
                (0x0000_0000_0000_0000, [0, 0]),
                (0x2000_0000_0000_0000, [1, 1]),
                (0x4000_0000_0000_0000, [2, 2]),
                (0x7000_0000_0000_0000, [0, 1]),
                (0x9000_0000_0000_0000, [0, 2]),
                (0xd000_0000_0000_0000, [1, 2]),
            ],
        ),
        (
            0.4,
            [
                (0x0000_0000_0000_0000, [0, 0]),
                (0x2000_0000_0000_0000, [1, 1]),
                (0x5000_0000_0000_0000, [2, 2]),
                (0xb000_0000_0000_0000, [0, 1]),
                (0xc000_0000_0000_0000, [0, 2]),
                (0xe000_0000_0000_0000, [1, 2]),
            ],
        ),
        (
            1.0,
            [
                (0x0000_0000_0000_0000, [0, 0]),
                (0x4000_0000_0000_0000, [1, 1]),
                (0xc000_0000_0000_0000, [2, 2]),
                (0xc000_0000_0000_0000, [2, 2]),
                (0xc000_0000_0000_0000, [2, 2]),
                (0xc000_0000_0000_0000, [2, 2]),
            ],
        ),
    ];

    for (inbreeding, cases) in vectors {
        for (word, expected) in cases {
            let mut source = ScriptedWords::new(&[word, SENTINEL]);
            let observed = f_diploid_genotype_with_source(&mut source, &frequencies, inbreeding)
                .expect("literal F vector");
            assert_eq!(observed, expected, "F={inbreeding}, word={word:#x}");
            assert_eq!(
                observed,
                independent_f_genotype(&frequencies, inbreeding, open_from_word(word)),
                "independent F oracle for F={inbreeding}, word={word:#x}"
            );
            assert_eq!(source.next_word(), SENTINEL, "F={inbreeding} consumption");
        }
    }
}

#[test]
fn f_degenerate_and_invalid_vectors_preserve_the_rng_sentinel() {
    const SENTINEL: u64 = 0xdeaf_beef_cafe_babe;
    let mut fixed = ScriptedWords::new(&[SENTINEL]);
    assert_eq!(
        f_diploid_genotype_with_source(&mut fixed, &[0.0, 1.0, 0.0], 0.4),
        Ok([1, 1])
    );
    assert_eq!(
        fixed.next_word(),
        SENTINEL,
        "single positive F consumes zero words"
    );

    let mut invalid = ScriptedWords::new(&[SENTINEL]);
    assert_eq!(
        f_diploid_genotype_with_source(&mut invalid, &[0.2, 0.3, 0.5], -0.1),
        Err(RandomSamplingError::InvalidInbreedingCoefficient)
    );
    assert_eq!(
        invalid.next_word(),
        SENTINEL,
        "invalid F consumes zero words"
    );
}

#[test]
fn hwe_and_f_preflight_degenerate_loci_and_formulas_without_poisson() {
    let root = RootSeed::from_bytes([17; 32]);
    let seed = derive_subseed(
        root,
        RuntimeSubseedScope::InitialGenotype {
            scenario_id: "direct",
            reference_id: "ref-main",
            population_id: "P",
            individual_id: "I",
            sequence_id: "chr1",
            position_1based: 1,
            reference_base: "A",
            alternate_alleles: "C",
        },
    )
    .expect("subseed");
    let mut fixed = seed.local_rng();
    let mut untouched = seed.local_rng();
    assert_eq!(
        hwe_composition(&mut fixed, &[1.0, 0.0], 10),
        Ok(vec![10, 0])
    );
    assert_eq!(
        fixed.next_u64(),
        untouched.next_u64(),
        "fixed HWE consumes zero words"
    );

    let mut invalid = seed.local_rng();
    let mut invalid_untouched = seed.local_rng();
    assert!(hwe_composition(&mut invalid, &[0.9, 0.2], 2).is_err());
    assert_eq!(
        invalid.next_u64(),
        invalid_untouched.next_u64(),
        "preflight consumes zero words"
    );

    let frequencies = [0.5, 0.2, 0.3];
    let inbreeding = 0.4;
    let mut producer = seed.local_rng();
    let mut oracle_rng = seed.local_rng();
    let open = (((oracle_rng.next_u64() >> 12) as f64) + 0.5) / 4_503_599_627_370_496.0;
    let expected = independent_f_genotype(&frequencies, inbreeding, open);
    assert_eq!(
        f_diploid_genotype(&mut producer, &frequencies, inbreeding).expect("F genotype"),
        expected,
        "independent F formula and category oracle"
    );
}

#[test]
fn compositions_cover_hwe_ploidies_one_two_and_ten_with_frozen_doses() {
    let reference = reference();
    let populations = format!(
        "{}{}{}",
        cohort("P1", 1, "hwe", ""),
        cohort("P2", 2, "hwe", ""),
        cohort("P10", 10, "hwe", ""),
    );
    let effective = config(&populations);
    let prepared = matrix(
        &reference,
        &[
            "VAR_IGNORED\tchr1\t1\tA\tC\tP1\t0.25\ttest",
            "VAR_IGNORED\tchr1\t1\tA\tC\tP2\t0.25\ttest",
            "VAR_IGNORED\tchr1\t1\tA\tC\tP10\t0.25\ttest",
        ],
    );
    let sampled =
        sample_initial_genotypes(&effective, &reference, &prepared).expect("compositions");
    assert_eq!(sampled.loci().len(), 3);
    let observed = sampled
        .loci()
        .iter()
        .map(|locus| {
            (
                locus.population_id().as_str(),
                locus.reference_dosage(),
                alternate_dosage(locus, b'C'),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(observed, [("P1", 0, 1), ("P10", 5, 5), ("P2", 1, 1)]);
    for locus in sampled.loci() {
        let total = locus.reference_dosage()
            + locus
                .alternate_dosages()
                .iter()
                .map(|dose| dose.dosage())
                .sum::<u8>();
        let ploidy = effective
            .initial_populations()
            .iter()
            .find(|population| population.population_id() == locus.population_id())
            .expect("cohort")
            .ploidy()
            .copy_count();
        assert_eq!(total, ploidy);
    }
}

#[test]
fn multiallelic_f_uses_one_canonical_locus_key_and_preserves_all_doses() {
    let reference = reference();
    let effective = config(&cohort(
        "P2",
        2,
        "f_diploid",
        "inbreeding_coefficient = 0.4",
    ));
    let prepared = matrix(
        &reference,
        &[
            "VAR_IGNORED\tchr1\t2\tC\tG\tP2\t0.3\ttest",
            "VAR_IGNORED\tchr1\t2\tC\tA\tP2\t0.2\ttest",
        ],
    );
    let sampled =
        sample_initial_genotypes(&effective, &reference, &prepared).expect("multiallelic F");
    let locus = sampled.loci().first().expect("one locus");
    assert_eq!(locus.reference(), b'C');
    assert_eq!(
        locus
            .alternate_dosages()
            .iter()
            .map(|dose| dose.alternate())
            .collect::<Vec<_>>(),
        vec![b'A', b'G'],
        "ALT list is canonical rather than input order"
    );
    assert_eq!(
        locus.reference_dosage()
            + locus
                .alternate_dosages()
                .iter()
                .map(|dose| dose.dosage())
                .sum::<u8>(),
        2
    );
    let root = effective.reproducibility().root_seed();
    let individual = effective.initial_populations()[0].individual_ids()[0].as_str();
    let seed = derive_subseed(
        root,
        RuntimeSubseedScope::InitialGenotype {
            scenario_id: "direct",
            reference_id: "ref-main",
            population_id: "P2",
            individual_id: individual,
            sequence_id: "chr1",
            position_1based: 2,
            reference_base: "C",
            alternate_alleles: "A,G",
        },
    )
    .expect("canonical locus subseed");
    let mut oracle_rng = seed.local_rng();
    let open = (((oracle_rng.next_u64() >> 12) as f64) + 0.5) / 4_503_599_627_370_496.0;
    let expected = independent_f_genotype(&[0.5, 0.2, 0.3], 0.4, open);
    let observed = [
        locus.reference_dosage(),
        alternate_dosage(locus, b'A'),
        alternate_dosage(locus, b'G'),
    ];
    let mut expected_doses = [0_u8; 3];
    expected_doses[usize::from(expected[0])] += 1;
    expected_doses[usize::from(expected[1])] += 1;
    assert_eq!(observed, expected_doses);
}

#[test]
fn composition_is_reproducible_and_isolated_from_population_and_row_changes() {
    let reference = reference();
    let base_config = config(&(cohort("P1", 2, "hwe", "") + &cohort("P2", 2, "hwe", "")));
    let rows = [
        "VAR_IGNORED\tchr1\t1\tA\tC\tP1\t0.25\ttest",
        "VAR_IGNORED\tchr1\t1\tA\tC\tP2\t0.50\ttest",
    ];
    let base = sample_initial_genotypes(&base_config, &reference, &matrix(&reference, &rows))
        .expect("base");
    let reordered = sample_initial_genotypes(
        &base_config,
        &reference,
        &matrix(&reference, &[rows[1], rows[0]]),
    )
    .expect("reordered");
    assert_eq!(base, reordered, "canonical preparation removes row order");

    let expanded_config = config(&format!(
        "{}{}{}",
        cohort("P1", 2, "hwe", ""),
        cohort("P2", 2, "hwe", ""),
        cohort("P3", 2, "hwe", ""),
    ));
    let expanded = sample_initial_genotypes(
        &expanded_config,
        &reference,
        &matrix(
            &reference,
            &[
                rows[0],
                rows[1],
                "VAR_IGNORED\tchr1\t1\tA\tC\tP3\t0.75\ttest",
            ],
        ),
    )
    .expect("expanded");
    let base_p1 = base
        .loci()
        .iter()
        .filter(|locus| locus.population_id().as_str() == "P1")
        .collect::<Vec<_>>();
    let expanded_p1 = expanded
        .loci()
        .iter()
        .filter(|locus| locus.population_id().as_str() == "P1")
        .collect::<Vec<_>>();
    assert_eq!(
        base_p1, expanded_p1,
        "P1 has a dedicated population/individual/locus stream"
    );
}

#[test]
fn composition_keeps_the_prepared_reference_contig_order() {
    let reference = reference_with_nonlexicographic_contigs();
    let effective = config(&cohort("P1", 2, "hwe", ""));
    let prepared = matrix(
        &reference,
        &[
            "VAR_IGNORED\tchr1\t1\tA\tC\tP1\t0.25\ttest",
            "VAR_IGNORED\tchr2\t1\tA\tC\tP1\t0.25\ttest",
        ],
    );
    let sampled =
        sample_initial_genotypes(&effective, &reference, &prepared).expect("compositions");
    let contigs = sampled
        .loci()
        .iter()
        .map(|locus| locus.contig().as_str())
        .collect::<Vec<_>>();
    assert_eq!(contigs, ["chr2", "chr1"]);
}

#[test]
fn hwe_empirical_dosage_statistics_are_frozen_and_match_binomial_expectation() {
    let mut histogram = [0_u64; 11];
    for value in 0_u64..4_096 {
        let mut root = [0_u8; 32];
        root[..8].copy_from_slice(&value.to_be_bytes());
        let mut rng = derive_subseed(
            RootSeed::from_bytes(root),
            RuntimeSubseedScope::InitialGenotype {
                scenario_id: "direct",
                reference_id: "ref-main",
                population_id: "P10",
                individual_id: "I-P10",
                sequence_id: "chr1",
                position_1based: 1,
                reference_base: "A",
                alternate_alleles: "C",
            },
        )
        .expect("subseed")
        .local_rng();
        let dosage = hwe_composition(&mut rng, &[0.75, 0.25], 10).expect("HWE")[1];
        histogram[usize::from(dosage)] += 1;
    }
    assert_eq!(
        histogram,
        [245, 767, 1153, 1051, 578, 230, 54, 17, 1, 0, 0],
        "frozen 4096-root corpus"
    );
    let observed_mean = histogram
        .iter()
        .enumerate()
        .map(|(dose, count)| dose as f64 * *count as f64)
        .sum::<f64>()
        / 4096.0;
    assert!((observed_mean - 2.5).abs() < 0.08, "mean={observed_mean}");
}

fn independent_f_genotype(frequencies: &[f64], inbreeding: f64, open: f64) -> [u8; 2] {
    let mut states = Vec::new();
    let mut weights = Vec::new();
    for (index, frequency) in frequencies.iter().copied().enumerate() {
        states.push([index as u8, index as u8]);
        weights.push(frequency * frequency + inbreeding * frequency * (1.0 - frequency));
    }
    for (left, left_frequency) in frequencies.iter().copied().enumerate() {
        for (right, right_frequency) in frequencies.iter().copied().enumerate().skip(left + 1) {
            states.push([left as u8, right as u8]);
            weights.push(2.0 * (1.0 - inbreeding) * left_frequency * right_frequency);
        }
    }
    let total = weights.iter().sum::<f64>();
    let target = open * total;
    let mut cumulative = 0.0;
    for (state, weight) in states.iter().copied().zip(weights) {
        cumulative += weight;
        if cumulative > target {
            return state;
        }
    }
    *states.last().expect("one allele creates one homozygote")
}

fn independent_multinomial(frequencies: &[f64], words: &[u64]) -> Vec<u8> {
    let mut dosage = vec![0_u8; frequencies.len()];
    for word in words {
        let target = open_from_word(*word);
        let mut cumulative = 0.0;
        for (index, frequency) in frequencies.iter().copied().enumerate() {
            cumulative += frequency;
            if cumulative > target {
                dosage[index] += 1;
                break;
            }
        }
    }
    dosage
}

fn open_from_word(word: u64) -> f64 {
    (((word >> 12) as f64) + 0.5) / 4_503_599_627_370_496.0
}
