use genome_population_simulator::fastq_output::paired_fastq_relative_paths;
use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::technical_profile::resolve_f4_technical_profile;
use genome_population_simulator::wgs_sampling::WgsPairedEndSampler;

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn fixture(bases: &[u8]) -> (Reference, Individual) {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "test").expect("sequence id"),
        FastaHeader::parse(b"chr1").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(bases).expect("DNA"),
    )
    .expect("reference sequence");
    let reference = Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference id"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference");
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual id"),
        &reference,
        Ploidy::new(1).expect("ploidy"),
    );
    (reference, individual)
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            _ => panic!("fixture DNA"),
        })
        .collect()
}

#[test]
fn paired_end_fragments_have_manual_fr_orientation_shared_template_and_no_orphans() {
    let bases = (0..500)
        .map(|index| match index % 5 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        })
        .collect::<Vec<_>>();
    let (reference, individual) = fixture(&bases);
    let identities = derive_technical_identity_plan(
        RootSeed::parse("root", ROOT).expect("root"),
        reference.id(),
        EventId::parse("event-pe", "test").expect("event"),
        &[individual.id().clone()],
    )
    .expect("identities");
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let pairs = WgsPairedEndSampler::new_for_pair_count(
        &reference,
        &individual,
        &identities,
        profile,
        150,
        3,
    )
    .expect("sampler")
    .collect::<Result<Vec<_>, _>>()
    .expect("pairs");

    assert_eq!(pairs.len(), 3);
    for pair in pairs {
        assert_eq!(
            pair.pre_hop_identity().template(),
            pair.fragment().template(),
            "la identidad PH se deriva una vez desde la misma plantilla",
        );
        assert!(
            pair.pre_hop_identity()
                .pre_hop_instance_id()
                .as_str()
                .starts_with("PH_"),
        );
        let fragment = pair.fragment().sequence();
        assert_eq!(pair.r1().sequence(), &fragment[..150]);
        assert_eq!(
            pair.r2().sequence(),
            reverse_complement(&fragment[fragment.len() - 150..])
        );
        assert_eq!(pair.r1().mate().number(), 1);
        assert_eq!(pair.r2().mate().number(), 2);
        assert_eq!(pair.r1().sequence().len(), 150);
        assert_eq!(pair.r2().sequence().len(), 150);
        assert_eq!(
            pair.fragment().end() - pair.fragment().start(),
            fragment.len() as u64
        );
        assert!(
            pair.r1()
                .phred()
                .iter()
                .zip(pair.r1().sequence())
                .all(|(q, b)| *q == if *b == b'N' { 0 } else { 40 })
        );
        assert!(
            pair.r2()
                .phred()
                .iter()
                .zip(pair.r2().sequence())
                .all(|(q, b)| *q == if *b == b'N' { 0 } else { 40 })
        );
    }
}

#[test]
fn paired_end_rejects_an_internal_pair_count_as_a_configuration_concept() {
    // La única frontera pública de autoría no admite `read_count`; el entero
    // aquí es exclusivamente un parámetro interno del productor para F4.2 y
    // para este oráculo de fragmentos.
    let toml = br#"
config_kind = "general"
schema_version = "2"

[reference]
id = "reference-main"
source = "simulated"
total_length = 500

[technical_event]
event_id = "event-pe"
individual_ids = ["sample-A"]
read_count = 1
"#;
    let path = std::env::temp_dir().join("gps-f4-read-count-rejected.toml");
    let admitted = genome_population_simulator::adapters::admit_toml_header(
        &path,
        toml.to_vec(),
        genome_population_simulator::configuration::ExpectedConfigKind::General,
    )
    .expect("F4 TOML header");
    assert!(genome_population_simulator::adapters::parse_declared_config(admitted).is_err());
}

#[test]
fn paired_end_uses_two_distinct_final_fastq_paths_without_a_single_end_alias() {
    let (r1, r2) = paired_fastq_relative_paths("event-pe", "sample-A");
    assert_eq!(
        r1.to_string_lossy(),
        "outputs/fastq/event-pe/sample-A_R1.fastq.gz"
    );
    assert_eq!(
        r2.to_string_lossy(),
        "outputs/fastq/event-pe/sample-A_R2.fastq.gz"
    );
    assert_ne!(r1, r2);
}
