//! Contrato F4.3.1: el único perfil es PCR-free y no expone un modelo de
//! duplicados. Los rechazos se prueban en la frontera TOML, antes de resolver
//! o generar cualquier pareja.

use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::results::{
    ArtifactRecord, ReferenceOutput, effective_config_document_with_reference,
};
use serde_json::json;

#[test]
fn only_pcr_free_profile_has_non_configurable_zero_duplicate_artifacts() {
    let effective = build_effective_config(&declared(valid_event())).expect("effective config");
    let event = effective.technical_event().expect("technical event");

    assert!(!event.pcr_duplicates_applicable());
    assert_eq!(event.realized_optical_duplicate_pairs(), 0);

    let mut reference = Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("sequence ID"),
                FastaHeader::parse(b"chr1").expect("FASTA header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(&[b'A'; 100]).expect("DNA"),
            )
            .expect("reference sequence"),
        ],
    )
    .expect("reference");
    let fasta = ArtifactRecord {
        id: "reference-main-fasta".to_owned(),
        kind: "reference_fasta".to_owned(),
        path: "references/reference-main.fasta".to_owned(),
        size_bytes: 100,
        digest: "a".repeat(64),
    };
    reference
        .bind_authoritative_resource(fasta.path.clone().into())
        .expect("authoritative resource");
    let document = effective_config_document_with_reference(
        &effective,
        ReferenceOutput {
            reference: &reference,
            fasta: &fasta,
        },
    )
    .expect("effective document");
    assert_eq!(
        document["technical_event"]["artifacts"],
        json!({"pcr": {"applicable": false}, "optical_duplicates": {"realized_pairs": 0}})
    );
}

#[test]
fn profiles_and_duplicate_rates_are_rejected_before_effective_resolution() {
    for (field, declaration) in [
        ("profile", "profile = \"amplified-wgs-v1\""),
        ("pcr_duplicate_rate", "pcr_duplicate_rate = 0.1"),
        ("optical_duplicate_rate", "optical_duplicate_rate = 0.1"),
    ] {
        let error = parse_declared_config(admitted(&format!("{}\n{declaration}\n", valid_event())))
            .expect_err("F4.3.1 must reject unsupported duplicate configuration");
        assert_eq!(error.code(), "GPS-CONF-006", "{field}");
        assert_eq!(
            error.context().get("field"),
            Some(&field.to_owned()),
            "{field}"
        );
    }
}

fn declared(source: &str) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    parse_declared_config(admitted(source)).expect("declared config")
}

fn admitted(source: &str) -> genome_population_simulator::configuration::AdmittedToml {
    admit_toml_header(
        Path::new("f4-profile-rejections.toml"),
        source.as_bytes().to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("TOML header")
}

fn valid_event() -> &'static str {
    r#"config_kind = "general"
schema_version = "2"

[reference]
kind = "simulated"
id = "reference-main"

[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 100

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 1

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
}
