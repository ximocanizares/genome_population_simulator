//! Cierre integrador de F2.2 sin introducir una segunda vía de ejecución.
//!
//! El recorrido ejercita los tres orígenes de la matriz común publicados por
//! F2.2.7 y el núcleo de cambios pequeños de F2.2.8--F2.2.10. Los oráculos no
//! reconstruyen los datos mediante los productores: inspeccionan los bytes
//! publicados, la secuencia FASTA y una verdad haplotípica escrita de forma
//! directa. La selección de individuos/copia en TOML queda deliberadamente
//! fuera de esta prueba: es la frontera observable de F2.3.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{
    EventId, IndividualId, ReferenceId, ScenarioId, SequenceId,
};
use genome_population_simulator::indel_application::SmallIndelLengths;
use genome_population_simulator::materialization::materialize_copy;
use genome_population_simulator::population_snvs::{
    TSV_HEADER, parse_canonical_tsv, parse_metadata_json, prepare_population_snvs_from_reference,
};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::results::validate_published_resources;
use genome_population_simulator::small_event_composition::{
    AppliedSmallEvent, DeclaredSmallIndelRequest, RandomSmallEventBlock, SampledSmallIndelCounts,
    apply_declared_small_events,
};
use genome_population_simulator::snp_application::{ApplicableSnp, SnpOrigin, TargetedSnp};

static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f2-2-12-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn generated_config() -> String {
    r#"config_kind = "general"
schema_version = "1"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
normalization_tolerance = 0.000000000001
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 2
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P2"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9
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
length = 12
"#
    .to_owned()
}

fn execute(config_path: &Path, destination: &Path) {
    execute_validated_config(
        &load_and_validate_config(config_path).expect("valid integrated configuration"),
        destination,
    )
    .expect("publication must succeed");
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    validate_published_resources(
        destination,
        &effective_path,
        &fs::read(&effective_path).expect("effective bytes"),
        &manifest_path,
        &fs::read(&manifest_path).expect("manifest bytes"),
    )
    .expect("the published physical resources validate together");
}

fn fasta_bases(path: &Path) -> Vec<u8> {
    fs::read_to_string(path)
        .expect("FASTA text")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(str::bytes)
        .collect()
}

#[test]
fn matrix_sources_are_reproducible_canonical_and_reusable_without_republication() {
    let root = temporary();
    fs::create_dir_all(&root).expect("root");
    let generated_config_path = root.join("generated.toml");
    fs::write(&generated_config_path, generated_config()).expect("generation TOML");

    let generated = root.join("generated");
    let generated_repeat = root.join("generated-repeat");
    execute(&generated_config_path, &generated);
    execute(&generated_config_path, &generated_repeat);

    let generated_matrix = fs::read(generated.join("population_snvs.tsv")).expect("matrix");
    assert_eq!(
        generated_matrix,
        fs::read(generated_repeat.join("population_snvs.tsv")).expect("repeated matrix"),
        "la misma semilla y configuración deben conservar las mismas decisiones F2.2"
    );
    let table = parse_canonical_tsv(&generated_matrix).expect("canonical generated TSV");
    let fasta = fasta_bases(&generated.join("resources/references/ref-main.fa"));
    assert_eq!(table.rows().len(), 4, "dos loci por dos poblaciones");
    for row in table.rows() {
        let position = row.field(2).parse::<usize>().expect("small test position") - 1;
        let reference = row.field(3).as_bytes()[0];
        let alternate = row.field(4).as_bytes()[0];
        let alt_af = row.field(6).parse::<f64>().expect("finite frequency");
        assert_eq!(reference, fasta[position], "oráculo FASTA independiente");
        assert!(matches!(reference, b'A' | b'C' | b'G' | b'T'));
        assert_ne!(reference, alternate);
        assert!((0.1..0.9).contains(&alt_af));
    }
    let metadata = parse_metadata_json(
        &fs::read(generated.join("population_snvs.metadata.json")).expect("metadata"),
    )
    .expect("closed metadata");
    assert_eq!(metadata.document()["provenance"], "generated");
    assert_eq!(metadata.document()["populations"][0]["population_id"], "P1");
    let generated_audit = metadata.document()["normalization_audit"].clone();
    assert_eq!(
        generated_audit.as_array().expect("generated audit").len(),
        4
    );
    assert_eq!(
        generated_audit[0]["tolerance"]
            .as_f64()
            .expect("audit tolerance")
            .to_bits(),
        1e-12_f64.to_bits(),
        "la generación debe registrar la tolerancia efectiva, no asumir 1e-9"
    );
    let generated_effective: serde_json::Value = serde_json::from_slice(
        &fs::read(generated.join("effective-config.json")).expect("generated effective"),
    )
    .expect("generated effective JSON");
    let generated_manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(generated.join("manifest.json")).expect("generated manifest"),
    )
    .expect("generated manifest JSON");
    assert_eq!(
        generated_effective["population_snvs"]["normalization_audit"],
        generated_audit
    );
    assert_eq!(
        generated_manifest["population_snvs"]["normalization_audit"],
        generated_audit
    );

    let reuse_config = root.join("reuse.toml");
    fs::write(
        &reuse_config,
        format!(
            "{}[population_snvs]\nmode = \"reuse\"\nmetadata_path = \"{}\"\n",
            generated_config()
                .replace("[population_snvs]\nmode = \"generate\"\n[population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 2\nspatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n[[population_snvs.generate.populations]]\npopulation_id = \"P2\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n[[population_snvs.generate.populations]]\npopulation_id = \"P1\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n", ""),
            generated.join("population_snvs.metadata.json").display()
        ),
    )
    .expect("reuse TOML");
    let source_snapshot = generated_matrix.clone();
    let reused = root.join("reused");
    execute(&reuse_config, &reused);
    assert_eq!(
        fs::read(generated.join("population_snvs.tsv")).expect("source matrix"),
        source_snapshot
    );
    assert!(!reused.join("population_snvs.tsv").exists());
    assert!(!reused.join("population_snvs.metadata.json").exists());
    let reuse_manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(reused.join("manifest.json")).expect("reuse manifest"))
            .expect("reuse manifest JSON");
    assert_eq!(reuse_manifest["population_snvs"]["mode"], "reuse");
    assert_eq!(
        reuse_manifest["population_snvs"]["normalization_audit"],
        generated_audit
    );

    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn supplied_matrix_provenance_reaches_an_explicit_snp_before_reproducible_random_indels() {
    let root = temporary();
    fs::create_dir_all(&root).expect("root");
    let fasta = root.join("reference.fa");
    fs::write(&fasta, b">chr1\nACGTNAC\n").expect("supplied FASTA");
    let supplied = root.join("supplied.tsv");
    let source_tsv = b"contig\tposition\tref\talt\tpopulation_id\talt_af\tsource_id\nchr1\t1\tA\tC\tP2\t0.4\tstudy-b-c\nchr1\t1\tA\tT\tP2\t0.6\tstudy-b-t\nchr1\t1\tA\tC\tP1\t0.50000000025\tstudy-a-c\nchr1\t1\tA\tT\tP1\t0.50000000025\tstudy-a-t\n";
    fs::write(&supplied, source_tsv).expect("source matrix");
    let config = root.join("supplied.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "1"
seed = "1111111111111111111111111111111111111111111111111111111111111111"
[population_snvs]
mode = "supplied"
input_path = "{}"
[reference]
kind = "supplied"
id = "ref-main"
fasta = "{}"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 7
"#,
            supplied.display(),
            fasta.display()
        ),
    )
    .expect("supplied TOML");
    let published = root.join("published");
    execute(&config, &published);
    assert_eq!(fs::read(&supplied).expect("input preserved"), source_tsv);
    let matrix = fs::read(published.join("population_snvs.tsv")).expect("published matrix");
    assert_eq!(
        matrix.split(|byte| *byte == b'\n').next(),
        Some(TSV_HEADER.join("\t").as_bytes())
    );
    let table = parse_canonical_tsv(&matrix).expect("canonical supplied matrix");
    assert_eq!(table.rows()[0].field(5), "P1");
    assert_eq!(table.rows()[0].field(7), "study-a-c");
    let metadata = parse_metadata_json(
        &fs::read(published.join("population_snvs.metadata.json")).expect("metadata"),
    )
    .expect("metadata");
    assert_eq!(metadata.document()["provenance"], "supplied");
    let audit = &metadata.document()["normalization_audit"];
    let p1_audit = audit
        .as_array()
        .expect("normalization audit list")
        .iter()
        .find(|group| group["population_id"] == "P1")
        .expect("P1 audit group");
    assert_eq!(p1_audit["decision"], "normalized");
    assert!(p1_audit["original_sum"].as_f64().expect("original sum") > 1.0);
    assert_eq!(p1_audit["reference_frequency"], 0.0);
    let mut forged_effective: serde_json::Value = serde_json::from_slice(
        &fs::read(published.join("effective-config.json")).expect("effective bytes"),
    )
    .expect("effective JSON");
    let mut forged_manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(published.join("manifest.json")).expect("manifest bytes"))
            .expect("manifest JSON");
    forged_effective["population_snvs"]["normalization_audit"][0]["original_sum"] =
        serde_json::json!(0.0);
    forged_manifest["population_snvs"]["normalization_audit"][0]["original_sum"] =
        serde_json::json!(0.0);
    assert!(
        validate_published_resources(
            &published,
            &published.join("effective-config.json"),
            &serde_json::to_vec(&forged_effective).expect("forged effective bytes"),
            &published.join("manifest.json"),
            &serde_json::to_vec(&forged_manifest).expect("forged manifest bytes"),
        )
        .is_err(),
        "la auditoría ligada debe coincidir con la metadata canónica"
    );

    let reference = in_memory_reference();
    let original = Individual::from_reference(
        IndividualId::parse("individual-1", "integration").expect("individual"),
        &reference,
        Ploidy::new(2).expect("ploidy"),
    );
    let prepared = prepare_population_snvs_from_reference(&table, &reference, 1e-9, false)
        .expect("published canonical matrix is semantically valid");
    let row = &prepared.rows()[0];
    assert_eq!(row.variant_id().as_str(), "VAR_QH5HOMDDZN3RBNXHNDST2X37DY");
    assert_eq!(row.population_id().as_str(), "P1");
    assert_eq!(row.source_id(), "study-a-c");
    let matrix_snp = ApplicableSnp::from_population_matrix(
        &reference,
        row,
        EventId::parse("matrix-snp", "integration").expect("event"),
    )
    .expect("matrix row applies to the authoritative reference");
    assert!(matches!(
        matrix_snp.origin(),
        SnpOrigin::PopulationMatrix {
            variant_id,
            population_id,
            source_id,
        } if variant_id == row.variant_id()
            && population_id == row.population_id()
            && source_id == row.source_id()
    ));
    let explicit = [
        genome_population_simulator::small_event_composition::DeclaredSmallEvent::Snp(
            TargetedSnp::new(
                HaplotypePhase::from_copy_index(0).expect("phase"),
                matrix_snp,
            ),
        ),
    ];
    let random = [RandomSmallEventBlock::SampledIndels {
        lengths: SmallIndelLengths::default(),
        request: DeclaredSmallIndelRequest::new(
            ScenarioId::parse("direct", "integration").expect("scenario"),
            SequenceId::parse("chr1", "integration").expect("sequence"),
            HaplotypePhase::from_copy_index(0).expect("phase"),
            "global".to_owned(),
            1.0,
            1.0,
        ),
        counts: SampledSmallIndelCounts::new(1, 1),
    }];
    let seed = RootSeed::from_bytes([23; 32]);
    let first = apply_declared_small_events(&original, &reference, seed, &explicit, &random)
        .expect("explicit prefix and sampled random suffix");
    let repeated = apply_declared_small_events(&original, &reference, seed, &explicit, &random)
        .expect("same root seed");
    assert_eq!(
        first, repeated,
        "el bloque aleatorio debe ser reproducible con su semilla"
    );
    assert_eq!(
        original,
        Individual::from_reference(
            IndividualId::parse("individual-1", "integration").expect("individual"),
            &reference,
            Ploidy::new(2).expect("ploidy"),
        )
    );
    assert!(matches!(
        first.applied_events(),
        [
            AppliedSmallEvent::Snp,
            AppliedSmallEvent::Indels {
                insertion_accepted: 1,
                deletion_accepted: 1,
                ..
            }
        ]
    ));
    let copy = &first
        .individual()
        .copies_for(&SequenceId::parse("chr1", "integration").expect("sequence"))
        .expect("copies")[0];
    let materialized = materialize_copy(copy, &reference).expect("materialized truth");
    // Oráculo independiente: para la semilla, eventos y referencia fijados,
    // la secuencia materializada esperada es la verdad concreta del flujo.
    assert_eq!(materialized, b"CGTANAC");

    fs::remove_dir_all(root).expect("cleanup");
}

fn in_memory_reference() -> Reference {
    let sequence_id = SequenceId::parse("chr1", "integration").expect("sequence");
    Reference::new(
        ReferenceId::parse("ref-main", "integration").expect("reference"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                sequence_id,
                FastaHeader::parse(b"chr1").expect("header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(b"ACGTNAC").expect("DNA"),
            )
            .expect("sequence"),
        ],
    )
    .expect("reference")
}
