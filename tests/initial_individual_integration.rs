//! Integración F2.3.4 con oráculos de secuencia escritos a mano.

use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{EffectiveConfig, build_effective_config};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::coordinates::BasePosition;
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::initial_genotypes::{
    assign_initial_phases, sample_initial_genotypes,
};
use genome_population_simulator::initial_individual_integration::{
    IncompatibleGenotypeAuthorization, InitialIndividualIntegrationError,
    InitialIndividualIntegrationRequest, integrate_initial_individual,
};
use genome_population_simulator::materialization::{SnpLocus, materialize_copy};
use genome_population_simulator::population_snvs::{
    PreparedPopulationSnvs, TSV_HEADER, parse_canonical_tsv, prepare_population_snvs_from_reference,
};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{CanonicalBase, ModificationProvenance};
use genome_population_simulator::small_event_composition::{
    DeclaredSmallEvent, TargetedCurrentDeletion, TargetedCurrentInsertion,
};
use genome_population_simulator::snp_application::{ApplicableSnp, TargetedSnp};

fn reference() -> Reference {
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("sequence ID"),
                FastaHeader::parse(b"chr1").expect("header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(b"ACGTAC").expect("DNA"),
            )
            .expect("reference sequence"),
        ],
    )
    .expect("reference")
}

fn config() -> EffectiveConfig {
    let source = r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["I-P1"]
ploidy = 2
genotype_model = "hwe"
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
length = 6
"#;
    let admitted = admit_toml_header(
        Path::new("initial-individual-integration.toml"),
        source.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    build_effective_config(&parse_declared_config(admitted).expect("declared config"))
        .expect("effective config")
}

fn multi_individual_config() -> EffectiveConfig {
    let source = r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I-P1-A", "I-P1-B"]
ploidy = 2
genotype_model = "hwe"
[[initial_populations]]
population_id = "P2"
individual_count = 1
individual_ids = ["I-P2"]
ploidy = 2
genotype_model = "hwe"
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
length = 6
"#;
    let admitted = admit_toml_header(
        Path::new("multi-initial-individual-integration.toml"),
        source.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    build_effective_config(&parse_declared_config(admitted).expect("declared config"))
        .expect("effective config")
}

fn matrix(reference: &Reference) -> PreparedPopulationSnvs {
    matrix_from_rows(reference, &["VAR_IGNORED\tchr1\t2\tC\tT\tP1\t1\tfixed-alt"])
}

fn matrix_from_rows(reference: &Reference, rows: &[&str]) -> PreparedPopulationSnvs {
    let mut text = TSV_HEADER.join("\t");
    text.push('\n');
    text.push_str(&rows.join("\n"));
    text.push('\n');
    prepare_population_snvs_from_reference(
        &parse_canonical_tsv(text.as_bytes()).expect("canonical TSV"),
        reference,
        1e-9,
        false,
    )
    .expect("prepared matrix")
}

fn phased_genotypes(
    config: &EffectiveConfig,
    reference: &Reference,
    matrix: &PreparedPopulationSnvs,
) -> genome_population_simulator::initial_genotypes::InitialPhasedGenotypes {
    assign_initial_phases(
        config.reproducibility().root_seed(),
        &sample_initial_genotypes(config, reference, matrix).expect("sampled genotypes"),
    )
    .expect("phased genotypes")
}

fn individual(reference: &Reference) -> Individual {
    Individual::from_reference(
        IndividualId::parse("I-P1", "test").expect("individual ID"),
        reference,
        Ploidy::DIPLOID,
    )
}

fn explicit_snp(alternate: u8) -> DeclaredSmallEvent {
    explicit_snp_on(0, alternate)
}

fn explicit_snp_on(phase: u8, alternate: u8) -> DeclaredSmallEvent {
    let reference_id = ReferenceId::parse("ref-main", "test").expect("reference ID");
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let position = BasePosition::new(reference_id, sequence_id, 1, 6).expect("position");
    let locus = SnpLocus::new(
        position,
        CanonicalBase::from_byte(b'C').expect("reference base"),
        CanonicalBase::from_byte(alternate).expect("alternate base"),
    )
    .expect("SNP locus");
    DeclaredSmallEvent::Snp(TargetedSnp::new(
        HaplotypePhase::from_copy_index(phase).expect("phase"),
        ApplicableSnp::explicit(
            locus,
            EventId::parse("explicit-snp", "test").expect("event ID"),
        ),
    ))
}

#[test]
fn f_equal_one_rejects_a_forced_heterozygote_before_applying_any_copy_change() {
    let reference = reference();
    let config = config();
    let matrix = matrix_from_rows(
        &reference,
        &[
            "VAR_IGNORED\tchr1\t2\tC\tG\tP1\t0.5\tfixed-g",
            "VAR_IGNORED\tchr1\t2\tC\tT\tP1\t0.5\tfixed-t",
        ],
    );
    let phased = phased_genotypes(&config, &reference, &matrix);
    let original = individual(&reference);
    let events = [explicit_snp_on(0, b'G'), explicit_snp_on(1, b'T')];
    let request = InitialIndividualIntegrationRequest::new(
        &matrix,
        &phased,
        config.reproducibility().root_seed(),
        &events,
        &[],
        IncompatibleGenotypeAuthorization::denied(),
        genome_population_simulator::application::EffectiveInitialGenotypeModel::FDiploid {
            inbreeding_coefficient: 1.0,
        },
    );
    assert!(matches!(
        integrate_initial_individual(&original, &reference, request),
        Err(InitialIndividualIntegrationError::IncompatibleExplicitSnp { .. })
    ));
    assert_eq!(copy_bases(&original, &reference, 0), b"ACGTAC");
}

fn copy_bases(individual: &Individual, reference: &Reference, phase: u8) -> Vec<u8> {
    materialize_copy(
        individual
            .copies_for(&SequenceId::parse("chr1", "test").expect("sequence ID"))
            .expect("copies")
            .iter()
            .find(|copy| copy.phase().copy_index() == phase)
            .expect("phase"),
        reference,
    )
    .expect("materialized test oracle")
}

fn request<'a>(
    config: &EffectiveConfig,
    matrix: &'a PreparedPopulationSnvs,
    phased: &'a genome_population_simulator::initial_genotypes::InitialPhasedGenotypes,
    explicit_events: &'a [DeclaredSmallEvent],
    authorization: IncompatibleGenotypeAuthorization,
) -> InitialIndividualIntegrationRequest<'a> {
    InitialIndividualIntegrationRequest::new(
        matrix,
        phased,
        config.reproducibility().root_seed(),
        explicit_events,
        &[],
        authorization,
        genome_population_simulator::application::EffectiveInitialGenotypeModel::HardyWeinberg,
    )
}

#[test]
fn phase_views_are_strictly_scoped_per_individual_across_cohorts() {
    let reference = reference();
    let config = multi_individual_config();
    let matrix = matrix_from_rows(
        &reference,
        &[
            "VAR_IGNORED\tchr1\t2\tC\tT\tP1\t1\tfixed-p1",
            "VAR_IGNORED\tchr1\t2\tC\tT\tP2\t1\tfixed-p2",
        ],
    );
    let phased = phased_genotypes(&config, &reference, &matrix);
    assert_eq!(phased.loci().len(), 3);
    for population in config.initial_populations() {
        for individual_id in population.individual_ids() {
            let view = phased.for_individual(individual_id);
            assert_eq!(view.loci().len(), 1);
            assert!(
                view.loci()
                    .iter()
                    .all(|locus| locus.composition().individual_id() == individual_id)
            );
            let individual =
                Individual::from_reference(individual_id.clone(), &reference, population.ploidy());
            integrate_initial_individual(
                &individual,
                &reference,
                InitialIndividualIntegrationRequest::new(
                    &matrix,
                    &view,
                    config.reproducibility().root_seed(),
                    &[],
                    &[],
                    IncompatibleGenotypeAuthorization::denied(),
                    population.genotype_model(),
                ),
            )
            .expect("each individual consumes only its own phased loci");
        }
    }
}

#[test]
fn phased_loci_reject_a_matrix_with_different_frequency_fingerprint() {
    let reference = reference();
    let config = config();
    let sampled_matrix =
        matrix_from_rows(&reference, &["VAR_IGNORED\tchr1\t2\tC\tT\tP1\t0.5\tfixed"]);
    let mismatched_matrix =
        matrix_from_rows(&reference, &["VAR_IGNORED\tchr1\t2\tC\tT\tP1\t0.75\tfixed"]);
    let phased = phased_genotypes(&config, &reference, &sampled_matrix);
    assert!(matches!(
        integrate_initial_individual(
            &individual(&reference),
            &reference,
            request(
                &config,
                &mismatched_matrix,
                &phased,
                &[],
                IncompatibleGenotypeAuthorization::denied(),
            ),
        ),
        Err(InitialIndividualIntegrationError::PhasedMatrixMismatch { .. })
    ));
}

#[test]
fn projected_matrix_snps_and_declared_indel_change_only_the_targeted_phased_copies() {
    let reference = reference();
    let config = config();
    let matrix = matrix(&reference);
    let phased = phased_genotypes(&config, &reference, &matrix);
    let original = individual(&reference);
    let insertion = DeclaredSmallEvent::Insertion(TargetedCurrentInsertion::new(
        SequenceId::parse("chr1", "test").expect("sequence ID"),
        HaplotypePhase::from_copy_index(0).expect("phase"),
        1,
        b"GG".as_slice(),
        ModificationProvenance::new(
            EventId::parse("explicit-insertion", "test").expect("event ID"),
        ),
    ));

    let integrated = integrate_initial_individual(
        &original,
        &reference,
        request(
            &config,
            &matrix,
            &phased,
            &[insertion],
            IncompatibleGenotypeAuthorization::denied(),
        ),
    )
    .expect("matrix SNP plus explicit indel");

    assert_eq!(
        copy_bases(integrated.individual(), &reference, 0),
        b"AGGTGTAC"
    );
    assert_eq!(
        copy_bases(integrated.individual(), &reference, 1),
        b"ATGTAC"
    );
    assert!(integrated.explicit_snp_constraints().is_empty());
    assert_eq!(copy_bases(&original, &reference, 0), b"ACGTAC");
    assert_eq!(copy_bases(&original, &reference, 1), b"ACGTAC");
}

#[test]
fn positive_probability_multiallelic_hwe_constraint_needs_no_authorization() {
    let reference = reference();
    let config = config();
    let matrix = matrix_from_rows(
        &reference,
        &[
            "VAR_IGNORED\tchr1\t2\tC\tG\tP1\t0.5\tfixed-g",
            "VAR_IGNORED\tchr1\t2\tC\tT\tP1\t0.5\tfixed-t",
        ],
    );
    let phased = phased_genotypes(&config, &reference, &matrix);
    let original = individual(&reference);
    let locus = &phased.loci()[0];
    let sampled = locus
        .alleles_by_phase()
        .iter()
        .find(|allele| allele.phase().copy_index() == 0)
        .and_then(|allele| {
            (allele.allele_index() > 0).then(|| {
                locus.composition().alternate_dosages()[usize::from(allele.allele_index() - 1)]
                    .alternate()
            })
        })
        .expect("all copies sample one of the two ALT alleles");
    let forced = if sampled == b'G' { b'T' } else { b'G' };
    let explicit = explicit_snp(forced);

    let integrated = integrate_initial_individual(
        &original,
        &reference,
        request(
            &config,
            &matrix,
            &phased,
            &[explicit],
            IncompatibleGenotypeAuthorization::denied(),
        ),
    )
    .expect("a positive-probability HWE restriction needs no authorization");
    assert_eq!(
        copy_bases(integrated.individual(), &reference, 0)[1],
        forced
    );
    assert_eq!(integrated.explicit_snp_constraints().len(), 1);
    assert!(!integrated.explicit_snp_constraints()[0].is_incompatible());
    assert_eq!(
        integrated.explicit_snp_constraints()[0].sampled_alternate(),
        Some(sampled)
    );
}

#[test]
fn authorization_never_relaxes_the_matrix_allele_catalog() {
    let reference = reference();
    let config = config();
    let matrix = matrix(&reference);
    let phased = phased_genotypes(&config, &reference, &matrix);
    let original = individual(&reference);

    assert!(matches!(
        integrate_initial_individual(
            &original,
            &reference,
            request(
                &config,
                &matrix,
                &phased,
                &[explicit_snp(b'G')],
                IncompatibleGenotypeAuthorization::allowed(),
            ),
        ),
        Err(
            InitialIndividualIntegrationError::ExplicitSnpAlleleOutsideMatrix {
                explicit_alternate: b'G',
                ..
            }
        )
    ));
    assert_eq!(copy_bases(&original, &reference, 0), b"ACGTAC");
}

#[test]
fn direct_event_list_cannot_relabel_a_matrix_snp_as_explicit() {
    let reference = reference();
    let config = config();
    let matrix = matrix(&reference);
    let phased = phased_genotypes(&config, &reference, &matrix);
    let original = individual(&reference);
    let matrix_snp = ApplicableSnp::from_population_matrix(
        &reference,
        &matrix.rows()[0],
        EventId::parse("mislabelled-matrix", "test").expect("event ID"),
    )
    .expect("matrix SNP");

    assert!(matches!(
        integrate_initial_individual(
            &original,
            &reference,
            request(
                &config,
                &matrix,
                &phased,
                &[DeclaredSmallEvent::Snp(TargetedSnp::new(
                    HaplotypePhase::from_copy_index(0).expect("phase"),
                    matrix_snp,
                ))],
                IncompatibleGenotypeAuthorization::denied(),
            ),
        ),
        Err(InitialIndividualIntegrationError::NonExplicitDirectSnp { .. })
    ));
}

#[test]
fn later_invalid_indel_rejects_the_complete_batch_without_publishing_matrix_snps() {
    let reference = reference();
    let config = config();
    let matrix = matrix(&reference);
    let phased = phased_genotypes(&config, &reference, &matrix);
    let original = individual(&reference);
    let invalid_deletion = DeclaredSmallEvent::Deletion(TargetedCurrentDeletion::new(
        SequenceId::parse("chr1", "test").expect("sequence ID"),
        HaplotypePhase::from_copy_index(0).expect("phase"),
        1,
        100,
        ModificationProvenance::new(EventId::parse("invalid-deletion", "test").expect("event ID")),
    ));

    assert!(matches!(
        integrate_initial_individual(
            &original,
            &reference,
            request(
                &config,
                &matrix,
                &phased,
                &[invalid_deletion],
                IncompatibleGenotypeAuthorization::denied(),
            ),
        ),
        Err(InitialIndividualIntegrationError::EventComposition(_))
    ));
    assert_eq!(copy_bases(&original, &reference, 0), b"ACGTAC");
    assert_eq!(copy_bases(&original, &reference, 1), b"ACGTAC");
}
