use genome_population_simulator::identity::{
    ExecutionId, F5_ID_CATALOG, F5BlockHaplotypeAssignmentScope, F5FounderCopyScope,
    F5IdentityContext, F5IdentityRegistry, F5PopulationIdentityContext, IndividualId, PopulationId,
    ReferenceId, ScenarioId, SequenceId,
};
use genome_population_simulator::randomness::RootSeed;

fn root() -> RootSeed {
    RootSeed::from_bytes([
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e,
        0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d,
        0x1e, 0x1f,
    ])
}

struct Inputs {
    execution_id: ExecutionId,
    scenario_id: ScenarioId,
    reference_id: ReferenceId,
    population_id: PopulationId,
    individual_id: IndividualId,
    sequence_id: SequenceId,
}

fn inputs() -> Inputs {
    Inputs {
        execution_id: ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test")
            .expect("execution"),
        scenario_id: ScenarioId::parse("scenario-f5", "test").expect("scenario"),
        reference_id: ReferenceId::parse("reference-main", "test").expect("reference"),
        population_id: PopulationId::parse("population-A", "test").expect("population"),
        individual_id: IndividualId::parse("individual-A", "test").expect("individual"),
        sequence_id: SequenceId::parse("sequence-1", "test").expect("sequence"),
    }
}

#[test]
fn f5_id_catalog_v4_derives_distinct_typed_prefixes_for_founder_outputs() {
    assert_eq!(F5_ID_CATALOG, "id-catalog-v4");
    let ids = inputs();
    let mut registry = F5IdentityRegistry::default();
    let base = F5IdentityContext::new(
        root(),
        &ids.execution_id,
        &ids.scenario_id,
        &ids.reference_id,
    );
    let population = F5PopulationIdentityContext::new(base, &ids.population_id);
    let block = registry.derive_block(base, 7).expect("block");
    assert_eq!(block.as_str(), "BLK_7AFGJODKAHE5N5Q35VXIEFNE2U");

    let site = registry.derive_site(base, 7, 5).expect("site");
    assert_eq!(site.as_str(), "SITE_5UJKXIUFLAKFCUYRMUXC7ZIXHU");

    let haplotype = registry
        .derive_haplotype(population, 7, 3, "101001")
        .expect("haplotype");
    assert_eq!(haplotype.as_str(), "HAP_LWRZFIHLKXNC3GHSPCIH6SJV7A");

    let founder_copy = registry
        .derive_founder_copy(
            population,
            F5FounderCopyScope::new(&ids.individual_id, &ids.sequence_id, 0),
        )
        .expect("founder copy");
    assert_eq!(founder_copy.as_str(), "FCP_3MYOBVN22OO7JPHBBXMUFWFJQQ");

    let assignment = registry
        .derive_block_haplotype_assignment(
            population,
            F5BlockHaplotypeAssignmentScope::new(7, &founder_copy, &haplotype),
        )
        .expect("assignment");
    assert_eq!(assignment.as_str(), "FHA_3H6PQOV7XJRT7HA3RGRDRDONFY");
}

#[test]
fn f5_founder_ids_are_sensitive_to_copy_index_and_assignment_links() {
    let ids = inputs();
    let mut registry = F5IdentityRegistry::default();
    let base = F5IdentityContext::new(
        root(),
        &ids.execution_id,
        &ids.scenario_id,
        &ids.reference_id,
    );
    let population = F5PopulationIdentityContext::new(base, &ids.population_id);
    let first = registry
        .derive_founder_copy(
            population,
            F5FounderCopyScope::new(&ids.individual_id, &ids.sequence_id, 0),
        )
        .expect("copy-zero founder copy");
    let second = registry
        .derive_founder_copy(
            population,
            F5FounderCopyScope::new(&ids.individual_id, &ids.sequence_id, 1),
        )
        .expect("copy-one founder copy");
    assert_ne!(first, second);

    let first_haplotype = registry
        .derive_haplotype(population, 7, 3, "101001")
        .expect("first haplotype");
    let second_haplotype = registry
        .derive_haplotype(population, 7, 4, "010110")
        .expect("second haplotype");
    let first_assignment = registry
        .derive_block_haplotype_assignment(
            population,
            F5BlockHaplotypeAssignmentScope::new(7, &first, &first_haplotype),
        )
        .expect("first assignment");
    let second_assignment = registry
        .derive_block_haplotype_assignment(
            population,
            F5BlockHaplotypeAssignmentScope::new(7, &first, &second_haplotype),
        )
        .expect("second assignment");
    let other_copy_assignment = registry
        .derive_block_haplotype_assignment(
            population,
            F5BlockHaplotypeAssignmentScope::new(7, &second, &first_haplotype),
        )
        .expect("other-copy assignment");
    assert_ne!(first_assignment, second_assignment);
    assert_ne!(first_assignment, other_copy_assignment);
}
