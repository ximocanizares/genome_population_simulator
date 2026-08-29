//! Cierre F3.5.4: una verdad biológica pequeña, conocida e independiente de
//! la evidencia recorre la configuración, referencia, individuo, `truth.jsonl`
//! y FASTQ final publicado.

use std::fs;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{
    EffectiveTechnicalEventTarget, load_and_validate_config,
};
use genome_population_simulator::coordinates::BasePosition;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::fastq_validation::validate_wgs_fastq;
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, SequenceId};
use genome_population_simulator::materialization::SnpLocus;
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_fasta::{SequenceExpectation, validate_supplied_fasta};
use genome_population_simulator::results::{ArtifactRecord, validate_published_resources};
use genome_population_simulator::segments::CanonicalBase;
use genome_population_simulator::snp_application::{ApplicableSnp, TargetedSnp, apply_snps};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::truth_validation::validate_and_reconstruct_truth_jsonl;
use noodles_bgzf as bgzf;

const ROOT_SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
const READ_LENGTH: usize = 150;
const UNKNOWN_POSITION_1BASED: usize = 10;
const MATRIX_POSITION_1BASED: usize = 25;
const ALTERED_POSITION_1BASED: usize = 75;

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn known_truth_flows_from_configuration_to_a_verifiable_published_fastq() {
    let root = temporary();
    fs::create_dir_all(&root).expect("temporary root");
    let mut reference_bases = vec![b'A'; READ_LENGTH];
    reference_bases[UNKNOWN_POSITION_1BASED - 1] = b'N';
    let mut matrix_haplotype = reference_bases.clone();
    matrix_haplotype[MATRIX_POSITION_1BASED - 1] = b'T';
    let mut expected_haplotype = matrix_haplotype.clone();
    expected_haplotype[ALTERED_POSITION_1BASED - 1] = b'G';
    let fasta = root.join("reference.fa");
    let matrix = root.join("population-snvs.tsv");
    let config = root.join("known-truth.toml");
    fs::write(
        &fasta,
        format!(
            ">chr-known\n{}\n",
            std::str::from_utf8(&reference_bases).expect("canonical reference")
        ),
    )
    .expect("known reference");
    fs::write(
        &matrix,
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr-known\t{MATRIX_POSITION_1BASED}\tA\tT\tP1\t1\n"
        ),
    )
    .expect("fixed matrix at a distinct known locus");
    fs::write(&config, configuration(&fasta, &matrix)).expect("configuration");

    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("effective configuration"),
        &destination,
    )
    .expect("complete publication");

    let manifest_path = destination.join("manifest.json");
    let effective_path = destination.join("effective-config.json");
    let manifest_bytes = fs::read(&manifest_path).expect("general manifest");
    let effective_bytes = fs::read(&effective_path).expect("effective configuration");
    validate_published_resources(
        &destination,
        &effective_path,
        &effective_bytes,
        &manifest_path,
        &manifest_bytes,
    )
    .expect("published reference, truth, event and FASTQ are mutually verifiable");
    assert_external_manifest_schema(&manifest_path, "manifest-general-v3.schema.json");

    let published_reference = validate_supplied_fasta(
        &destination.join("resources/references/reference-known.fa"),
        "reference-known",
        &[SequenceExpectation {
            id: "chr-known",
            class: SequenceClass::Chromosome,
            length: READ_LENGTH as u64,
        }],
    )
    .expect("published authoritative reference");
    let truth_path = destination.join("outputs/truth.jsonl");
    let reconstructed = validate_and_reconstruct_truth_jsonl(
        &mut BufReader::new(fs::File::open(&truth_path).expect("truth JSONL")),
        published_reference.reference(),
    )
    .expect("truth reconstructs independently from the published reference");
    assert_eq!(reconstructed.document.header.copy_count, 2);
    assert_eq!(reconstructed.document.header.genotype_count, 2);
    assert_eq!(reconstructed.copies.len(), 2);
    assert_eq!(
        reconstructed
            .copies
            .values()
            .filter(|copy| *copy == &expected_haplotype)
            .count(),
        1,
        "the direct phase-0 SNP must alter exactly one copy"
    );
    assert_eq!(
        reconstructed
            .copies
            .values()
            .filter(|copy| *copy == &matrix_haplotype)
            .count(),
        1,
        "the fixed matrix must leave the other copy without the direct SNP"
    );

    let event_path = destination.join("events/reads-known/manifest.json");
    assert_external_manifest_schema(&event_path, "manifest-event-v1.schema.json");
    let event: serde_json::Value =
        serde_json::from_slice(&fs::read(&event_path).expect("event manifest")).expect("JSON");
    assert_eq!(event["state"], "complete");
    assert_eq!(event["event_id"], "reads-known");
    assert_eq!(event["event_type"], "technical_sequencing");
    assert_eq!(event["reference"]["reference_id"], "reference-known");
    assert_eq!(
        event["selected_individuals"][0]["individual_id"],
        "individual-known"
    );
    let artifact = event["artifacts"]
        .as_array()
        .expect("event artifacts")
        .first()
        .expect("one selected library");
    assert_eq!(artifact["state"], "complete");
    assert_eq!(artifact["individual_id"], "individual-known");
    let fastq_artifact = artifact_record(&artifact["artifact"]);
    let fastq_path = destination.join(fastq_artifact.path.as_str());
    assert!(fastq_artifact.path.ends_with(".fastq.gz"));
    let mut gzip_reader =
        bgzf::io::Reader::new(fs::File::open(&fastq_path).expect("published FASTQ"));
    let mut fastq = Vec::new();
    gzip_reader
        .read_to_end(&mut fastq)
        .expect("decompress published FASTQ");
    let lines = fastq
        .split_inclusive(|byte| *byte == b'\n')
        .collect::<Vec<_>>();
    assert_eq!(lines.len(), 4);
    assert!(
        lines[0].starts_with(b"@TPL_") && lines[0].ends_with(b"_0\n"),
        "the published record preserves the canonical first template ordinal"
    );
    let observed = &lines[1][..READ_LENGTH];
    assert!(
        observed == expected_haplotype || observed == reverse_complement(&expected_haplotype),
        "the one eligible 150-base interval must come from the independently known haplotype"
    );
    assert_eq!(lines[2], b"+\n");
    let expected_qualities = observed
        .iter()
        .map(|base| if *base == b'N' { b'!' } else { b'I' })
        .chain(std::iter::once(b'\n'))
        .collect::<Vec<_>>();
    assert_eq!(lines[3], expected_qualities);
    assert_eq!(
        observed.iter().filter(|base| **base == b'N').count(),
        1,
        "the unknown base from the authoritative reference remains literal in FASTQ"
    );

    let known_individual = known_individual(published_reference.reference());
    let identities = derive_technical_identity_plan(
        RootSeed::parse("test root seed", ROOT_SEED).expect("root seed"),
        published_reference.reference().id(),
        EventId::parse("reads-known", "test").expect("event ID"),
        &[IndividualId::parse("individual-known", "test").expect("individual ID")],
    )
    .expect("known technical identities");
    let report = validate_wgs_fastq(
        &fastq_path,
        &fastq_artifact,
        published_reference.reference(),
        &known_individual,
        &identities,
        EffectiveTechnicalEventTarget::ReadCount(1),
    )
    .expect("independent FASTQ replay against known truth");
    assert_eq!(report.read_count(), 1);
    assert_eq!(report.verified_origin_count(), 1);
    assert_eq!(report.unknown_base_count(), 1);

    // La WGS básica se calcula antes de escribir el FASTQ: un único fragmento
    // de 150 bases sobre una referencia haploide de 150 bases es 1x. La
    // cobertura final se vuelve a medir desde el FASTQ publicado. En F3 ambas
    // coinciden porque el contrato mínimo no descarta ni añade lecturas; no se
    // confunden ni se interpreta esta igualdad como cobertura alineada.
    let haploid_length = published_reference.reference().total_length();
    let basic_wgs_coverage = READ_LENGTH as f64 / haploid_length as f64;
    let final_fastq_coverage = report.base_count() as f64 / haploid_length as f64;
    assert_eq!(basic_wgs_coverage, 1.0);
    assert_eq!(final_fastq_coverage, 1.0);
    assert_eq!(report.coverage().realized_coverage(), final_fastq_coverage);

    // La verdad se valida antes y después de dañar la evidencia: no se deriva
    // del FASTQ ni el FASTQ puede repararla. Su validación independiente sí
    // rechaza inmediatamente la evidencia dañada.
    fs::write(&fastq_path, b"@broken\nA\n+\nI\n").expect("corrupt FASTQ");
    validate_and_reconstruct_truth_jsonl(
        &mut BufReader::new(fs::File::open(&truth_path).expect("truth JSONL after corruption")),
        published_reference.reference(),
    )
    .expect("truth remains independently valid after FASTQ corruption");
    assert!(
        validate_wgs_fastq(
            &fastq_path,
            &fastq_artifact,
            published_reference.reference(),
            &known_individual,
            &identities,
            EffectiveTechnicalEventTarget::ReadCount(1),
        )
        .is_err(),
        "the independent FASTQ replay must reject the corrupted evidence"
    );

    fs::remove_dir_all(root).expect("cleanup");
}

fn configuration(fasta: &Path, matrix: &Path) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "1"
seed = "{ROOT_SEED}"
[population_snvs]
mode = "supplied"
input_path = "{}"
[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["individual-known"]
ploidy = 2
genotype_model = "hwe"
[[direct_events]]
event_id = "known-phase-0-snp"
individual_id = "individual-known"
phase = 0
kind = "snp"
contig = "chr-known"
position_1based = {ALTERED_POSITION_1BASED}
ref = "A"
alt = "G"
[reference]
kind = "supplied"
id = "reference-known"
fasta = "{}"
[[reference.sequences]]
id = "chr-known"
class = "chromosome"
length = {READ_LENGTH}
[technical_event]
event_id = "reads-known"
individual_ids = ["individual-known"]
read_count = 1
"#,
        matrix.display(),
        fasta.display(),
    )
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
            _ => panic!("known truth contains only DNA"),
        })
        .collect()
}

fn known_individual(reference: &genome_population_simulator::reference::Reference) -> Individual {
    let individual_id = IndividualId::parse("individual-known", "test").expect("individual ID");
    let baseline = Individual::from_reference(individual_id, reference, Ploidy::DIPLOID);
    let sequence_id = SequenceId::parse("chr-known", "test").expect("sequence ID");
    let matrix_locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            sequence_id.clone(),
            (MATRIX_POSITION_1BASED - 1) as u64,
            READ_LENGTH as u64,
        )
        .expect("known matrix position"),
        CanonicalBase::from_byte(b'A').expect("known matrix reference base"),
        CanonicalBase::from_byte(b'T').expect("known matrix alternate base"),
    )
    .expect("known matrix SNP locus");
    let matrix_applications = [0_u8, 1]
        .into_iter()
        .map(|phase| {
            TargetedSnp::new(
                HaplotypePhase::from_copy_index(phase).expect("known matrix phase"),
                ApplicableSnp::explicit(
                    matrix_locus.clone(),
                    EventId::parse("known-matrix-oracle", "test").expect("matrix oracle event ID"),
                ),
            )
        })
        .collect::<Vec<_>>();
    let matrix_individual =
        apply_snps(&baseline, reference, &matrix_applications).expect("independent matrix truth");
    let locus = SnpLocus::new(
        BasePosition::new(
            reference.id().clone(),
            sequence_id,
            (ALTERED_POSITION_1BASED - 1) as u64,
            READ_LENGTH as u64,
        )
        .expect("known position"),
        CanonicalBase::from_byte(b'A').expect("known reference base"),
        CanonicalBase::from_byte(b'G').expect("known alternate base"),
    )
    .expect("known SNP locus");
    let applications = [0_u8]
        .into_iter()
        .map(|phase| {
            TargetedSnp::new(
                HaplotypePhase::from_copy_index(phase).expect("known phase"),
                ApplicableSnp::explicit(
                    locus.clone(),
                    EventId::parse("known-truth-oracle", "test").expect("oracle event ID"),
                ),
            )
        })
        .collect::<Vec<_>>();
    apply_snps(&matrix_individual, reference, &applications).expect("independent known truth")
}

fn artifact_record(value: &serde_json::Value) -> ArtifactRecord {
    ArtifactRecord {
        id: value["id"].as_str().expect("artifact ID").to_owned(),
        kind: value["kind"].as_str().expect("artifact kind").to_owned(),
        path: value["path"].as_str().expect("artifact path").to_owned(),
        size_bytes: value["size_bytes"].as_u64().expect("artifact size"),
        digest: value["checksum"]["digest"]
            .as_str()
            .expect("artifact digest")
            .to_owned(),
    }
}

fn assert_external_manifest_schema(path: &Path, schema: &str) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(path)
        .output()
        .expect("external JSON Schema validator");
    assert!(
        output.status.success(),
        "external validator for {path:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f3-end-to-end-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
    ))
}
