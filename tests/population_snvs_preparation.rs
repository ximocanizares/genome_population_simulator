use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::ReferenceId;
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::population_snvs::{
    PopulationSnvsPreparationError, TSV_HEADER, parse_canonical_tsv, parse_supplied_tsv,
    prepare_supplied_population_snvs,
};
use genome_population_simulator::results::resolve_manifest_reference;
use serde_json::json;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f2-2-4-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

#[test]
fn supplied_parser_accepts_only_six_or_seven_columns() {
    let six = parse_supplied_tsv(
        b"contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t1\tA\tC\tP\t0\n",
    )
    .expect("six");
    assert_eq!(six.rows()[0].field(0), "");
    let seven = parse_supplied_tsv(
        b"contig\tposition\tref\talt\tpopulation_id\talt_af\tsource_id\nchr1\t1\tA\tC\tP\t0\ts\n",
    )
    .expect("seven");
    assert_eq!(seven.rows()[0].field(7), "s");
    assert!(parse_supplied_tsv(b"contig\tposition\tref\talt\tpopulation_id\n").is_err());
}

fn reduction_fixture_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/reference_reduction")
        .join(name)
}

/// Publica un bundle F1 completo, con fuente, configuración efectiva, mapa y
/// FASTA derivado. Así F2 ejerce la misma frontera física que reutilizará.
fn publish_complete_reduction_fixture(root: &Path) -> PathBuf {
    fs::create_dir_all(root).expect("fixture root");
    fs::copy(reduction_fixture_path("source.fa"), root.join("source.fa")).expect("source FASTA");
    fs::copy(
        reduction_fixture_path("general.toml"),
        root.join("general.toml"),
    )
    .expect("general configuration");
    let destination = root.join("published");
    let effective = load_and_validate_config(&root.join("general.toml"))
        .expect("valid F1 reduction configuration");
    execute_validated_config(&effective, &destination).expect("complete F1 publication");
    destination
}

fn publish_f1_bundle(root: &Path, source_fasta: &[u8], configuration: &str) -> PathBuf {
    fs::create_dir_all(root).expect("fixture root");
    fs::write(root.join("source.fa"), source_fasta).expect("source FASTA");
    fs::write(root.join("general.toml"), configuration).expect("general configuration");
    let destination = root.join("published");
    let effective = load_and_validate_config(&root.join("general.toml"))
        .expect("valid F1 reduction configuration");
    execute_validated_config(&effective, &destination).expect("complete F1 publication");
    destination
}

fn resolved(
    root: &Path,
    reference_id: &str,
) -> genome_population_simulator::results::ResolvedManifestReference {
    let destination = publish_complete_reduction_fixture(root);
    resolve_manifest_reference(
        &destination.join("manifest.json"),
        &ReferenceId::parse(reference_id, "test").expect("ID"),
        None,
    )
    .expect("resolved manifest")
}

fn table(rows: &[&str]) -> genome_population_simulator::population_snvs::PopulationSnvsTsv {
    let mut bytes = TSV_HEADER.join("\t");
    bytes.push('\n');
    bytes.push_str(&rows.join("\n"));
    bytes.push('\n');
    parse_canonical_tsv(bytes.as_bytes()).expect("shape")
}

#[test]
fn manifest_resolution_rejects_changed_fasta_bytes() {
    let root = temporary();
    let resolved = resolved(&root, "reference-main");
    fs::write(resolved.fasta_path(), b">chr1\nAATG\n").expect("mutate");
    let error = resolve_manifest_reference(
        &root.join("published/manifest.json"),
        resolved.reference_id(),
        None,
    )
    .expect_err("digest mismatch");
    assert_eq!(error.code(), "GPS-IO-002");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn manifest_resolution_rejects_corrupted_reduction_map_with_intact_fasta_and_hashes() {
    let root = temporary();
    let destination = publish_complete_reduction_fixture(&root);
    let manifest_path = destination.join("manifest.json");
    let original = resolve_manifest_reference(
        &manifest_path,
        &ReferenceId::parse("reference-small", "test").expect("ID"),
        None,
    )
    .expect("complete derived reference resolves before corruption");
    let original_fasta = identify_file(original.fasta_path(), Path::new("derived.fa"))
        .expect("derived FASTA identity");

    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("manifest bytes"))
            .expect("manifest JSON");
    let map_artifact = manifest["reference_reductions"][0]["map_artifact_id"]
        .as_str()
        .expect("map artifact ID")
        .to_owned();
    let map_relative = manifest["artifacts"]
        .as_array()
        .expect("artifact catalog")
        .iter()
        .find(|artifact| artifact["id"] == map_artifact)
        .and_then(|artifact| artifact["path"].as_str())
        .expect("map path")
        .to_owned();
    let map_path = destination.join(&map_relative);
    let mut map: serde_json::Value =
        serde_json::from_slice(&fs::read(&map_path).expect("map bytes")).expect("map JSON");
    map["algorithm"] = json!("corrupted-map-v1");
    fs::write(
        &map_path,
        serde_json::to_vec(&map).expect("corrupted map JSON"),
    )
    .expect("write corrupted map");
    let identity = identify_file(&map_path, Path::new(&map_relative)).expect("map identity");
    let artifact = manifest["artifacts"]
        .as_array_mut()
        .expect("artifact catalog")
        .iter_mut()
        .find(|artifact| artifact["id"] == map_artifact)
        .expect("map artifact record");
    artifact["size_bytes"] = json!(identity.size_bytes());
    artifact["checksum"]["digest"] = json!(identity.digest());
    fs::write(
        &manifest_path,
        serde_json::to_vec(&manifest).expect("updated manifest JSON"),
    )
    .expect("updated manifest");

    let unchanged_fasta = identify_file(original.fasta_path(), Path::new("derived.fa"))
        .expect("derived FASTA identity after map mutation");
    assert_eq!(
        unchanged_fasta, original_fasta,
        "map mutation must not touch FASTA"
    );
    let error = resolve_manifest_reference(
        &manifest_path,
        &ReferenceId::parse("reference-small", "test").expect("ID"),
        None,
    )
    .expect_err("physical F1 validation must reject the corrupted map");
    assert_eq!(error.code(), "GPS-FMT-006");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn preparation_uses_independent_fasta_truth_and_canonical_order() {
    let root = temporary();
    let reference = resolved(&root, "reference-main");
    let prepared = prepare_supplied_population_snvs(
        &table(&[
            "x\tchr1\t2\tC\tA\tP2\t0.2\t",
            "x\tchr1\t2\tC\tA\tP1\t0.2\t",
            "x\tscaffold-1\t1\tT\tG\tP2\t1\t",
            "x\tscaffold-1\t1\tT\tG\tP1\t1\t",
        ]),
        &reference,
        1e-9,
        false,
    )
    .expect("prepared");
    assert_eq!(prepared.reference_id(), reference.reference_id());
    assert_eq!(prepared.reference_fingerprint().len(), 64);
    assert_eq!(prepared.rows().len(), 4);
    // Oráculo independiente: la fila reordenada físicamente conserva su orden canónico.
    assert_eq!(prepared.rows()[0].contig().as_str(), "chr1");
    assert_eq!(prepared.rows()[0].population_id().as_str(), "P1");
    let wrong = prepare_supplied_population_snvs(
        &table(&["x\tchr1\t2\tA\tC\tP1\t1\t"]),
        &reference,
        1e-9,
        false,
    )
    .expect_err("REF wrong");
    assert!(matches!(
        wrong,
        PopulationSnvsPreparationError::ReferenceMismatch { .. }
    ));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn preparation_preserves_nonlexical_fasta_order_and_nonutf8_description() {
    let root = temporary();
    let configuration = r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[reference]
kind = "supplied"
id = "reference-main"
fasta = "source.fa"

[[reference.sequences]]
id = "zeta"
class = "chromosome"
length = 4

[[reference.sequences]]
id = "alpha"
class = "scaffold"
length = 4

[[reduced_references]]
id = "reference-small"
source_reference_id = "reference-main"
region_size_bases = 4
target_bases = 4
"#;
    let destination = publish_f1_bundle(
        &root,
        b">zeta description \xff\nACGT\n>alpha\nTTAA\n",
        configuration,
    );
    let reference = resolve_manifest_reference(
        &destination.join("manifest.json"),
        &ReferenceId::parse("reference-main", "test").expect("ID"),
        None,
    )
    .expect("nonlexical F1 reference with reversible non-UTF-8 description");
    let prepared = prepare_supplied_population_snvs(
        &table(&["a\talpha\t1\tT\tC\tP\t1\t", "z\tzeta\t1\tA\tC\tP\t1\t"]),
        &reference,
        1e-9,
        false,
    )
    .expect("prepared in FASTA catalog order");
    assert_eq!(prepared.rows()[0].contig().as_str(), "zeta");
    assert_eq!(prepared.rows()[1].contig().as_str(), "alpha");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn preparation_rejects_duplicate_and_incomplete_multiallelic_coverage() {
    let root = temporary();
    let reference = resolved(&root, "reference-main");
    let duplicate = prepare_supplied_population_snvs(
        &table(&["x\tchr1\t1\tA\tC\tP\t0.1\t", "y\tchr1\t1\tA\tC\tP\t0.1\t"]),
        &reference,
        1e-9,
        false,
    )
    .expect_err("duplicate");
    assert!(matches!(
        duplicate,
        PopulationSnvsPreparationError::DuplicateAltPopulation { .. }
    ));
    let incomplete = prepare_supplied_population_snvs(
        &table(&[
            "x\tchr1\t1\tA\tC\tP1\t0.1\t",
            "x\tchr1\t1\tA\tC\tP2\t0.1\t",
            "y\tchr1\t1\tA\tG\tP1\t0.1\t",
        ]),
        &reference,
        1e-9,
        false,
    )
    .expect_err("coverage");
    assert!(matches!(
        incomplete,
        PopulationSnvsPreparationError::IncompletePopulationCoverage { .. }
    ));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn preparation_normalizes_and_filters_zero_without_dropping_fixed_or_authorization() {
    let root = temporary();
    let reference = resolved(&root, "reference-main");
    let prepared = prepare_supplied_population_snvs(
        &table(&[
            "x\tchr1\t1\tA\tC\tP1\t0\t",
            "x\tchr1\t1\tA\tC\tP2\t0\t",
            "y\tchr1\t1\tA\tG\tP1\t0.5000005\t",
            "y\tchr1\t1\tA\tG\tP2\t1\t",
            "z\tchr1\t1\tA\tT\tP1\t0.5\t",
            "z\tchr1\t1\tA\tT\tP2\t0\t",
        ]),
        &reference,
        1e-6,
        false,
    )
    .expect("normalizes");
    assert_eq!(prepared.discarded_globally_absent().len(), 1);
    assert_eq!(prepared.rows().len(), 4); // G fixed in P2 remains; multiallelic C disappeared only.
    let zero = table(&["x\tchr1\t1\tA\tC\tP\t0\t"]);
    assert!(matches!(
        prepare_supplied_population_snvs(&zero, &reference, 1e-9, false),
        Err(PopulationSnvsPreparationError::EmptyActiveMatrixRequiresAuthorization)
    ));
    let authorized = prepare_supplied_population_snvs(&zero, &reference, 1e-9, true)
        .expect("specific authorization");
    assert!(authorized.rows().is_empty() && authorized.empty_authorized());
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn preparation_treats_every_negative_zero_alt_frequency_as_globally_absent() {
    let root = temporary();
    let reference = resolved(&root, "reference-main");
    let all_negative_zero = table(&["x\tchr1\t1\tA\tC\tP1\t-0\t", "x\tchr1\t1\tA\tC\tP2\t-0\t"]);
    assert!(matches!(
        prepare_supplied_population_snvs(&all_negative_zero, &reference, 1e-9, false),
        Err(PopulationSnvsPreparationError::EmptyActiveMatrixRequiresAuthorization)
    ));
    let authorized = prepare_supplied_population_snvs(&all_negative_zero, &reference, 1e-9, true)
        .expect("all negative zero frequencies are absent");
    assert!(authorized.rows().is_empty());
    assert_eq!(authorized.discarded_globally_absent().len(), 1);
    fs::remove_dir_all(root).expect("cleanup");
}
