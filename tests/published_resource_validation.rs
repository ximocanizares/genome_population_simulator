use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::results::{
    json_document_bytes, validate_general_manifest_bytes, validate_published_resources,
};
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);

#[test]
fn correlated_documents_cannot_replace_physical_authority() {
    for (label, old, new) in [
        (
            "toml-seed",
            "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
            "1111111111111111111111111111111111111111111111111111111111111111",
        ),
        (
            "toml-region-size",
            "region_size_bases = 5",
            "region_size_bases = 4",
        ),
        (
            "toml-reference-id",
            "id = \"reference-main\"",
            "id = \"reference-other\"",
        ),
        (
            "toml-sequence-class",
            "class = \"chromosome\"",
            "class = \"scaffold\"",
        ),
        (
            "toml-byte-digest",
            "config_kind = \"general\"",
            "config_kind = \"general\" # same semantics, different bytes",
        ),
    ] {
        assert_physical_rejection(label, supplied_fixture(), |root, _effective, manifest| {
            replace_configuration_text(root, manifest, old, new);
        });
    }

    assert_physical_rejection(
        "toml-composition",
        simulated_fixture(),
        |root, _effective, manifest| {
            replace_configuration_text(root, manifest, "A = 0.25", "A = 0.1");
            replace_configuration_text(root, manifest, "C = 0.25", "C = 0.4");
        },
    );

    assert_physical_rejection(
        "region-size",
        supplied_fixture(),
        |_root, effective, manifest| {
            effective["reduced_references"][0]["region_size_bases"] = json!(4);
            manifest["reference_reductions"][0]["region_size_bases"] = json!(4);
        },
    );

    assert_physical_rejection(
        "target-bases-vs-toml",
        supplied_fixture(),
        |_root, effective, manifest| {
            effective["reduced_references"][0]["target"]["declared"] = json!("12");
            effective["reduced_references"][0]["target"]["bases"] = json!(12);
            manifest["reference_reductions"][0]["target"]["declared"] = json!("12");
            manifest["reference_reductions"][0]["target"]["bases"] = json!(12);
        },
    );

    let percent_root = temporary("percent-fixture");
    let percent_config = write_percentage_fixture(&percent_root);
    assert_physical_rejection(
        "percent-99-9-vs-toml",
        percent_config,
        |_root, effective, manifest| {
            effective["reduced_references"][0]["target"]["declared"] = json!("99.9");
            effective["reduced_references"][0]["target"]["bases"] = json!(1);
            manifest["reference_reductions"][0]["target"]["declared"] = json!("99.9");
            manifest["reference_reductions"][0]["target"]["bases"] = json!(1);
        },
    );
    fs::remove_dir_all(percent_root).expect("percent fixture cleanup");

    assert_physical_rejection(
        "reference-ids",
        supplied_fixture(),
        |root, effective, manifest| {
            effective["reference"]["id"] = json!("reference-other");
            effective["reduced_references"][0]["source_reference_id"] = json!("reference-other");
            manifest["references"][0]["id"] = json!("reference-other");
            manifest["reference_reductions"][0]["source_reference_id"] = json!("reference-other");
            let map_path = root.join(artifact_path(manifest, "reference-map-reference-small"));
            let mut map = read_json(&map_path);
            map["source_reference_id"] = json!("reference-other");
            write_json(&map_path, &map);
            refresh_artifact(root, manifest, "reference-map-reference-small");
        },
    );

    assert_physical_rejection(
        "fasta-length",
        supplied_fixture(),
        |root, effective, manifest| {
            let fasta_path = root.join(artifact_path(manifest, "reference-fasta"));
            let fasta = fs::read_to_string(&fasta_path)
                .expect("source FASTA")
                .replace("ACGTNNACGTAC", "ACGTNNACGTA");
            fs::write(&fasta_path, fasta).expect("mutated source FASTA");
            effective["reference"]["sequences"][0]["length"] = json!(11);
            effective["reference"]["total_genome_size"] = json!(19);
            refresh_artifact(root, manifest, "reference-fasta");
        },
    );

    assert_physical_rejection(
        "composition",
        simulated_fixture(),
        |_root, effective, _manifest| {
            effective["reference"]["composition"] =
                json!({"A": 0.4, "C": 0.1, "G": 0.25, "T": 0.25});
        },
    );

    for (label, pointer, replacement) in [
        (
            "reference-original-path",
            "/reference/original_fasta_path",
            json!("other.fa"),
        ),
        (
            "reference-header",
            "/reference/sequences/0/header/value",
            json!("changed-header"),
        ),
        (
            "reference-description",
            "/reference/sequences/0/description",
            json!({"encoding": "utf-8", "value": "invented"}),
        ),
        (
            "reference-alphabet",
            "/reference/sequences/0/alphabet_validation",
            json!("ACGTN-case-insensitive-correlated"),
        ),
    ] {
        assert_physical_rejection(
            label,
            supplied_fixture(),
            move |_root, effective, _manifest| {
                *effective
                    .pointer_mut(pointer)
                    .expect("reference metadata pointer") = replacement;
            },
        );
    }

    assert_physical_rejection(
        "canonical-format",
        simulated_fixture(),
        |_root, effective, _manifest| {
            effective["reference"]["canonical_format"]["line_base_count"] = json!(79);
        },
    );

    assert_physical_rejection(
        "reduction-map",
        supplied_fixture(),
        |root, _effective, manifest| {
            let map_path = root.join(artifact_path(manifest, "reference-map-reference-small"));
            let mut map = read_json(&map_path);
            map["algorithm"] = json!("correlated-but-invalid");
            write_json(&map_path, &map);
            refresh_artifact(root, manifest, "reference-map-reference-small");
        },
    );
}

#[test]
fn every_derived_sequence_catalog_field_is_checked_independently() {
    for (label, pointer, replacement) in [
        (
            "derived-sequence-id",
            "/reduced_references/0/sequences/0/id",
            json!("changed-sequence"),
        ),
        (
            "derived-sequence-length",
            "/reduced_references/0/sequences/0/length",
            json!(1),
        ),
        (
            "derived-sequence-header",
            "/reduced_references/0/sequences/0/header/value",
            json!("changed-header"),
        ),
        (
            "derived-sequence-description",
            "/reduced_references/0/sequences/0/description",
            json!({"encoding": "utf-8", "value": "invented"}),
        ),
        (
            "derived-sequence-alphabet",
            "/reduced_references/0/sequences/0/alphabet_validation",
            json!("ACGT"),
        ),
    ] {
        assert_physical_rejection(
            label,
            supplied_fixture(),
            move |_root, effective, _manifest| {
                *effective.pointer_mut(pointer).expect(label) = replacement;
            },
        );
    }

    assert_physical_rejection(
        "derived-sequence-class",
        supplied_fixture(),
        |_root, effective, _manifest| {
            let class = &mut effective["reduced_references"][0]["sequences"][0]["class"];
            *class = if class == "chromosome" {
                json!("scaffold")
            } else {
                json!("chromosome")
            };
        },
    );
    assert_physical_rejection(
        "derived-sequence-n-count",
        supplied_fixture(),
        |_root, effective, _manifest| {
            let n_count = &mut effective["reduced_references"][0]["sequences"][0]["n_count"];
            let observed = n_count.as_u64().expect("n_count");
            *n_count = json!(if observed == 0 { 1 } else { observed - 1 });
        },
    );
}

#[test]
fn execution_id_is_bound_to_csprng_seed_and_exact_toml_digest() {
    let parent = temporary("execution-id-csprng");
    fs::create_dir(&parent).expect("temporary root");
    let config =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/config/general-minimal.toml");
    let first = parent.join("first");
    let second = parent.join("second");
    execute_validated_config(
        &load_and_validate_config(&config).expect("first seedless config"),
        &first,
    )
    .expect("first execution");
    execute_validated_config(
        &load_and_validate_config(&config).expect("second seedless config"),
        &second,
    )
    .expect("second execution");

    let mut first_manifest = read_json(&first.join("manifest.json"));
    let second_manifest = read_json(&second.join("manifest.json"));
    assert_ne!(
        first_manifest["execution_id"], second_manifest["execution_id"],
        "independent CSPRNG roots must derive distinct execution IDs"
    );
    for kind in ["effective_config", "reference_fasta"] {
        let first_id = artifact_id_by_kind(&first_manifest, kind).to_owned();
        let second_id = artifact_id_by_kind(&second_manifest, kind);
        assert_eq!(first_id, second_id, "fixture artifact IDs remain stable");
        let from = second.join(artifact_path(&second_manifest, second_id));
        let to = first.join(artifact_path(&first_manifest, &first_id));
        fs::copy(from, to).expect("replace correlated artifact");
        refresh_artifact(&first, &mut first_manifest, &first_id);
    }
    write_json(&first.join("manifest.json"), &first_manifest);

    let effective_id = artifact_id_by_kind(&first_manifest, "effective_config");
    let effective_path = first.join(artifact_path(&first_manifest, effective_id));
    let manifest_path = first.join("manifest.json");
    let error = validate_published_resources(
        &first,
        &effective_path,
        &fs::read(&effective_path).expect("effective bytes"),
        &manifest_path,
        &fs::read(&manifest_path).expect("manifest bytes"),
    )
    .expect_err("execution ID from the first root cannot authenticate the second root");
    assert_eq!(error.code(), "GPS-FMT-006");
    fs::remove_dir_all(parent).expect("cleanup");
}

#[test]
fn direct_events_must_match_the_typed_toml_resolution_exactly() {
    let fixture_root = temporary("direct-events-fixture");
    fs::create_dir(&fixture_root).expect("fixture root");
    let config = write_direct_events_fixture(&fixture_root);
    for (label, pointer, replacement) in [
        (
            "direct-individual",
            "/direct_events/0/individual_id",
            json!("I2"),
        ),
        ("direct-phase", "/direct_events/0/phase", json!(1)),
        ("direct-contig", "/direct_events/0/contig", json!("chr2")),
        (
            "direct-position",
            "/direct_events/0/position_1based",
            json!(2),
        ),
        ("direct-reference", "/direct_events/0/ref", json!("C")),
        ("direct-alternate", "/direct_events/0/alt", json!("T")),
    ] {
        assert_physical_rejection(label, config.clone(), move |_root, effective, _manifest| {
            *effective.pointer_mut(pointer).expect(label) = replacement;
        });
    }
    assert_physical_rejection(
        "initial-statistics-realized-count",
        config.clone(),
        |_root, effective, manifest| {
            let count = effective["population_snvs"]["initial_population_statistics"]
                ["populations"][0]["loci"][0]["genotypes"][0]["expected_count"]
                .as_f64()
                .expect("expected count");
            effective["population_snvs"]["initial_population_statistics"]["populations"][0]["loci"]
                [0]["genotypes"][0]["expected_count"] = json!(count + 1.0);
            manifest["population_snvs"]["initial_population_statistics"] =
                effective["population_snvs"]["initial_population_statistics"].clone();
        },
    );
    assert_physical_rejection(
        "initial-populations-must-match-toml-exactly",
        config.clone(),
        |_root, effective, _manifest| {
            effective["initial_populations"][0]["individual_ids"][0] = json!("I2");
        },
    );
    assert_physical_rejection(
        "initial-populations-require-statistics",
        config.clone(),
        |_root, effective, manifest| {
            effective["population_snvs"]
                .as_object_mut()
                .expect("effective population SNVs")
                .remove("initial_population_statistics");
            manifest["population_snvs"]
                .as_object_mut()
                .expect("manifest population SNVs")
                .remove("initial_population_statistics");
        },
    );
    assert_physical_rejection("direct-class", config, |_root, effective, _manifest| {
        effective["direct_events"][0] = json!({
            "event_id": "event-1",
            "individual_id": "I1",
            "phase": 0,
            "kind": "deletion",
            "contig": "chr1",
            "start": 1,
            "end": 2,
        });
    });
    fs::remove_dir_all(fixture_root).expect("fixture cleanup");
}

#[test]
fn published_founder_keys_are_revalidated_against_the_toml_contract() {
    let fixture_root = temporary("founder-keys-fixture");
    fs::create_dir(&fixture_root).expect("fixture root");
    let config = write_founder_keys_fixture(&fixture_root);
    let destination = temporary("founder-keys-result");
    let effective_config = load_and_validate_config(&config).expect("valid fixture");
    execute_validated_config(&effective_config, &destination).expect("complete baseline");

    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let mut effective = read_json(&effective_path);
    let mut manifest = read_json(&manifest_path);
    assert_eq!(
        effective["initial_populations"][0]["founder_keys"],
        json!(["founder-I1"])
    );
    effective["initial_populations"][0]["founder_keys"][0] = json!("other-founder");
    write_json(&effective_path, &effective);
    refresh_artifact(&destination, &mut manifest, "effective-config");
    write_json(&manifest_path, &manifest);

    let error = validate_published_resources(
        &destination,
        &effective_path,
        &fs::read(&effective_path).expect("effective bytes"),
        &manifest_path,
        &fs::read(&manifest_path).expect("manifest bytes"),
    )
    .expect_err("revalidation must bind founder keys to the typed TOML resolution");
    assert_eq!(error.code(), "GPS-FMT-006");

    fs::remove_dir_all(destination).expect("result cleanup");
    fs::remove_dir_all(fixture_root).expect("fixture cleanup");
}

#[test]
fn hostile_genotype_dosage_overflow_is_a_controlled_revalidation_error() {
    let fixture_root = temporary("hostile-dosage-fixture");
    fs::create_dir(&fixture_root).expect("fixture root");
    let config = write_direct_events_fixture(&fixture_root);
    let destination = temporary("hostile-dosage-result");
    let effective_config = load_and_validate_config(&config).expect("valid fixture");
    execute_validated_config(&effective_config, &destination).expect("complete baseline");

    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let mut effective = read_json(&effective_path);
    let mut manifest = read_json(&manifest_path);
    effective["population_snvs"]["initial_population_statistics"]["populations"][0]["loci"][0]["genotypes"]
        [0]["allele_dosages"] = json!([u64::MAX, 3]);
    manifest["population_snvs"]["initial_population_statistics"] =
        effective["population_snvs"]["initial_population_statistics"].clone();
    write_json(&effective_path, &effective);
    refresh_artifact(&destination, &mut manifest, "effective-config");
    write_json(&manifest_path, &manifest);

    let effective_bytes = fs::read(&effective_path).expect("effective bytes");
    let manifest_bytes = fs::read(&manifest_path).expect("manifest bytes");
    assert!(
        validate_published_resources(
            &destination,
            &effective_path,
            &effective_bytes,
            &manifest_path,
            &manifest_bytes,
        )
        .is_err(),
        "la revalidación no debe aceptar ni desbordar con dosis hostiles"
    );

    fs::remove_dir_all(destination).expect("result cleanup");
    fs::remove_dir_all(fixture_root).expect("fixture cleanup");
}

fn assert_physical_rejection(
    label: &str,
    config: PathBuf,
    mutate: impl FnOnce(&Path, &mut Value, &mut Value),
) {
    let parent = temporary(label);
    fs::create_dir(&parent).expect("temporary root");
    let destination = parent.join("result");
    let effective_config = load_and_validate_config(&config).expect("valid fixture");
    execute_validated_config(&effective_config, &destination).expect("complete baseline");

    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let mut effective = read_json(&effective_path);
    let mut manifest = read_json(&manifest_path);
    mutate(&destination, &mut effective, &mut manifest);
    write_json(&effective_path, &effective);
    refresh_artifact(&destination, &mut manifest, "effective-config");
    write_json(&manifest_path, &manifest);

    let effective_bytes = fs::read(&effective_path).expect("effective bytes");
    let manifest_bytes = fs::read(&manifest_path).expect("manifest bytes");
    validate_general_manifest_bytes(&manifest_path, &manifest_bytes)
        .expect("correlated catalog remains structurally valid");
    let validation = validate_published_resources(
        &destination,
        &effective_path,
        &effective_bytes,
        &manifest_path,
        &manifest_bytes,
    );
    let error = match validation {
        Ok(()) => panic!("case {label}: physical authority accepted correlated declarations"),
        Err(error) => error,
    };
    assert_eq!(error.code(), "GPS-FMT-006", "case {label}: {error:?}");

    fs::remove_dir_all(parent).expect("cleanup");
}

fn refresh_artifact(root: &Path, manifest: &mut Value, id: &str) {
    let artifacts = manifest["artifacts"]
        .as_array_mut()
        .expect("artifact catalog");
    let artifact = artifacts
        .iter_mut()
        .find(|artifact| artifact["id"] == id)
        .unwrap_or_else(|| panic!("artifact {id}"));
    let relative = artifact["path"].as_str().expect("artifact path");
    let identity = identify_file(&root.join(relative), Path::new(relative)).expect("identity");
    artifact["size_bytes"] = json!(identity.size_bytes());
    artifact["checksum"]["digest"] = json!(identity.digest());
}

fn replace_configuration_text(root: &Path, manifest: &mut Value, old: &str, new: &str) {
    let id = manifest["configuration_input"]["artifact_id"]
        .as_str()
        .expect("configuration artifact ID")
        .to_owned();
    let path = root.join(artifact_path(manifest, &id));
    let original = fs::read_to_string(&path).expect("configuration copy");
    let changed = original.replace(old, new);
    assert_ne!(
        changed, original,
        "configuration mutation must be effective"
    );
    fs::write(&path, changed).expect("mutated configuration copy");
    refresh_artifact(root, manifest, &id);
}

fn artifact_path<'a>(manifest: &'a Value, id: &str) -> &'a str {
    manifest["artifacts"]
        .as_array()
        .expect("artifact catalog")
        .iter()
        .find(|artifact| artifact["id"] == id)
        .unwrap_or_else(|| panic!("artifact {id}"))["path"]
        .as_str()
        .expect("artifact path")
}

fn artifact_id_by_kind<'a>(manifest: &'a Value, kind: &str) -> &'a str {
    manifest["artifacts"]
        .as_array()
        .expect("artifact catalog")
        .iter()
        .find(|artifact| artifact["kind"] == kind)
        .unwrap_or_else(|| panic!("artifact kind {kind}"))["id"]
        .as_str()
        .expect("artifact ID")
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("JSON bytes")).expect("JSON document")
}

fn write_json(path: &Path, value: &Value) {
    fs::write(path, json_document_bytes(value).expect("canonical JSON")).expect("write JSON");
}

fn supplied_fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/reference_reduction/general.toml")
}

fn simulated_fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/config/general-seeded.toml")
}

fn write_percentage_fixture(root: &Path) -> PathBuf {
    fs::create_dir(root).expect("percent fixture root");
    let mut fasta = String::from(">chr1\n");
    fasta.push_str(&"A".repeat(1000));
    fasta.push('\n');
    fs::write(root.join("source.fa"), fasta).expect("percent fixture FASTA");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"{}\"\n\
             [reference]\nkind = \"supplied\"\nid = \"reference-main\"\nfasta = \"source.fa\"\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1000\n\
             [[reduced_references]]\nid = \"reference-small\"\nsource_reference_id = \"reference-main\"\n\
             region_size_bases = 10\nreduction_percent = 25\n",
            "1".repeat(64)
        ),
    )
    .expect("percent fixture config");
    config
}

fn write_direct_events_fixture(root: &Path) -> PathBuf {
    let fasta = root.join("reference.fa");
    let matrix = root.join("matrix.tsv");
    fs::write(&fasta, b">chr1\nACGT\n").expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t1\tA\tG\tP1\t0.5\n",
    )
    .expect("matrix");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"{}\"\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"{}\"\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nindividual_ids = [\"I1\"]\nploidy = 2\ngenotype_model = \"hwe\"\n\
             [[direct_events]]\nevent_id = \"event-1\"\nindividual_id = \"I1\"\nphase = 0\nkind = \"snp\"\ncontig = \"chr1\"\nposition_1based = 1\nref = \"A\"\nalt = \"G\"\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"{}\"\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 4\n",
            "0".repeat(64),
            matrix.display(),
            fasta.display(),
        ),
    )
    .expect("direct events config");
    config
}

fn write_founder_keys_fixture(root: &Path) -> PathBuf {
    let fasta = root.join("reference.fa");
    let matrix = root.join("matrix.tsv");
    fs::write(&fasta, b">chr1\nACGT\n").expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t1\tA\tG\tP1\t0.5\n",
    )
    .expect("matrix");
    let config = root.join("founder-keys.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"{}\"\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"{}\"\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nfounder_keys = [\"founder-I1\"]\nploidy = 2\ngenotype_model = \"hwe\"\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"{}\"\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 4\n",
            "0".repeat(64),
            matrix.display(),
            fasta.display(),
        ),
    )
    .expect("founder keys config");
    config
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f14-physical-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
