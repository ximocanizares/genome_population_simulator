//! Oráculos independientes para F2.3.5: bytes finales y frecuencias finitas.

use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{EffectiveConfig, build_effective_config};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::haplotype::HaplotypePhase;
use genome_population_simulator::identity::{EventId, ReferenceId, SequenceId};
use genome_population_simulator::initial_individual_integration::IncompatibleGenotypeAuthorization;
use genome_population_simulator::initial_individual_materialization::materialize_initial_individuals;
use genome_population_simulator::population_snvs::{
    PreparedPopulationSnvs, TSV_HEADER, parse_canonical_tsv, prepare_population_snvs_from_reference,
};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::ModificationProvenance;
use genome_population_simulator::small_event_composition::{
    DeclaredSmallEvent, TargetedCurrentDeletion,
};

fn reference() -> Reference {
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("sequence ID"),
                FastaHeader::parse(b"chr1").expect("header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(b"ACGT").expect("DNA"),
            )
            .expect("reference sequence"),
        ],
    )
    .expect("reference")
}

fn config(populations: &str) -> EffectiveConfig {
    let source = format!(
        r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
{populations}
[reference]
kind = "simulated"
id = "ref-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 4
"#,
    );
    let admitted = admit_toml_header(
        Path::new("initial-individual-materialization.toml"),
        source.into_bytes(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    build_effective_config(&parse_declared_config(admitted).expect("declared")).expect("effective")
}

fn matrix(reference: &Reference, rows: &[&str]) -> PreparedPopulationSnvs {
    let bytes = format!("{}\n{}\n", TSV_HEADER.join("\t"), rows.join("\n"));
    prepare_population_snvs_from_reference(
        &parse_canonical_tsv(bytes.as_bytes()).expect("TSV"),
        reference,
        1e-9,
        false,
    )
    .expect("prepared matrix")
}

#[test]
fn deletion_is_explicitly_missing_and_uses_only_observed_copies_for_frequencies() {
    let reference = reference();
    let config = config(
        r#"[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe""#,
    );
    let matrix = matrix(&reference, &["VAR_IGNORED\tchr1\t2\tC\tT\tP1\t1\tfixed"]);
    let output = materialize_initial_individuals(
        &config,
        &reference,
        &matrix,
        IncompatibleGenotypeAuthorization::denied(),
        |individual_id| {
            if individual_id.as_str() == "I1" {
                Ok(vec![DeclaredSmallEvent::Deletion(
                    TargetedCurrentDeletion::new(
                        SequenceId::parse("chr1", "test").expect("sequence ID"),
                        HaplotypePhase::from_copy_index(0).expect("phase"),
                        1,
                        2,
                        ModificationProvenance::new(
                            EventId::parse("delete-c", "test").expect("event"),
                        ),
                    ),
                )])
            } else {
                Ok(Vec::new())
            }
        },
    )
    .expect("materialized founders");

    let i1 = output
        .individuals()
        .iter()
        .find(|individual| individual.individual().id().as_str() == "I1")
        .expect("I1");
    // Oráculo manual: todas las copias parten de ACGT y C→T; la fase 0 de I1
    // pierde exactamente la segunda base, sin reajustar la otra copia.
    assert_eq!(i1.copies()[0].bases(), b"AGT");
    assert_eq!(i1.copies()[1].bases(), b"ATGT");

    let locus = &output.statistics()["populations"][0]["loci"][0];
    assert_eq!(locus["total_copies"], 4);
    assert_eq!(locus["observed_copies"], 3);
    assert_eq!(locus["missing_by_deletion"], 1);
    assert_eq!(locus["alleles"][1]["realized_count"], 3);
    assert_eq!(locus["alleles"][1]["realized_frequency"], 1.0);
    assert_eq!(locus["complete_genotype_observations"], 1);
}

#[test]
fn adding_or_reordering_entities_does_not_change_an_existing_individual() {
    let reference = reference();
    let base = config(
        r#"[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["I1"]
ploidy = 2
genotype_model = "hwe""#,
    );
    let extended = config(
        r#"[[initial_populations]]
population_id = "P2"
individual_count = 1
individual_ids = ["J1"]
ploidy = 2
genotype_model = "hwe"
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I2", "I1"]
ploidy = 2
genotype_model = "hwe""#,
    );
    let matrix = matrix(
        &reference,
        &[
            "VAR_IGNORED\tchr1\t2\tC\tT\tP1\t0.5\tfixed-p1",
            "VAR_IGNORED\tchr1\t2\tC\tT\tP2\t0.5\tfixed-p2",
        ],
    );
    let founders = |config: &EffectiveConfig| {
        materialize_initial_individuals(
            config,
            &reference,
            &matrix,
            IncompatibleGenotypeAuthorization::denied(),
            |_| Ok(Vec::new()),
        )
        .expect("founders")
    };
    let base_i1 = founders(&base)
        .individuals()
        .iter()
        .find(|individual| individual.individual().id().as_str() == "I1")
        .expect("base I1")
        .copies()
        .iter()
        .map(|copy| copy.bases().to_vec())
        .collect::<Vec<_>>();
    let extended_i1 = founders(&extended)
        .individuals()
        .iter()
        .find(|individual| individual.individual().id().as_str() == "I1")
        .expect("extended I1")
        .copies()
        .iter()
        .map(|copy| copy.bases().to_vec())
        .collect::<Vec<_>>();
    assert_eq!(base_i1, extended_i1);
}

#[test]
fn independent_hwe_and_f_oracles_cover_models_ploidies_and_zero_observations() {
    for ploidy in 1..=10 {
        let config = config(&format!(
            r#"[[initial_populations]]
population_id = "P1"
individual_count = 3
individual_ids = ["I1", "I2", "I3"]
ploidy = {ploidy}
genotype_model = "hwe""#
        ));
        let reference = reference();
        let matrix = matrix(&reference, &["VAR_IGNORED\tchr1\t2\tC\tT\tP1\t0.25\tfixed"]);
        let output = materialize_initial_individuals(
            &config,
            &reference,
            &matrix,
            IncompatibleGenotypeAuthorization::denied(),
            |_| Ok(Vec::new()),
        )
        .expect("HWE founders");
        let genotypes = output.statistics()["populations"][0]["loci"][0]["genotypes"]
            .as_array()
            .expect("HWE genotypes");
        let sum = genotypes
            .iter()
            .map(|row| row["expected_probability"].as_f64().expect("probability"))
            .sum::<f64>();
        assert!((sum - 1.0).abs() < 1e-12, "P={ploidy}");
        let one_alt = genotypes
            .iter()
            .find(|row| row["allele_dosages"] == serde_json::json!([ploidy - 1, 1]))
            .expect("one ALT dosage");
        let expected = f64::from(ploidy) * 0.25 * 0.75_f64.powi(ploidy - 1);
        assert!(
            (one_alt["expected_probability"]
                .as_f64()
                .expect("probability")
                - expected)
                .abs()
                < 1e-12
        );
    }

    let reference = reference();
    let f_config = config(
        r#"[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "f_diploid"
inbreeding_coefficient = 0.4"#,
    );
    let f_matrix = matrix(
        &reference,
        &[
            "VAR_IGNORED\tchr1\t2\tC\tG\tP1\t0.2\tfixed-g",
            "VAR_IGNORED\tchr1\t2\tC\tT\tP1\t0.3\tfixed-t",
        ],
    );
    let f_output = materialize_initial_individuals(
        &f_config,
        &reference,
        &f_matrix,
        IncompatibleGenotypeAuthorization::denied(),
        |_| Ok(Vec::new()),
    )
    .expect("F founders");
    let genotypes = f_output.statistics()["populations"][0]["loci"][0]["genotypes"]
        .as_array()
        .expect("F genotypes");
    let heterozygote = genotypes
        .iter()
        .find(|row| row["allele_dosages"] == serde_json::json!([1, 1, 0]))
        .expect("REF/G heterozygote");
    assert!(
        (heterozygote["expected_probability"]
            .as_f64()
            .expect("probability")
            - 0.12)
            .abs()
            < 1e-12
    );

    let deleted = materialize_initial_individuals(
        &config(
            r#"[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["I1"]
ploidy = 2
genotype_model = "hwe""#,
        ),
        &reference,
        &matrix(&reference, &["VAR_IGNORED\tchr1\t2\tC\tT\tP1\t1\tfixed"]),
        IncompatibleGenotypeAuthorization::denied(),
        |_| {
            Ok((0..2)
                .map(|phase| {
                    DeclaredSmallEvent::Deletion(TargetedCurrentDeletion::new(
                        SequenceId::parse("chr1", "test").expect("sequence ID"),
                        HaplotypePhase::from_copy_index(phase).expect("phase"),
                        1,
                        2,
                        ModificationProvenance::new(
                            EventId::parse(format!("delete-{phase}"), "test").expect("event"),
                        ),
                    ))
                })
                .collect())
        },
    )
    .expect("all deleted founders");
    let locus = &deleted.statistics()["populations"][0]["loci"][0];
    assert_eq!(locus["observed_copies"], 0);
    assert_eq!(locus["missing_by_deletion"], 2);
    assert!(
        locus["alleles"]
            .as_array()
            .expect("alleles")
            .iter()
            .all(|allele| allele["realized_frequency"] == 0.0)
    );
}
