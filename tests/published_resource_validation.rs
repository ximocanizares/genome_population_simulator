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
        None,
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
        None,
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
            None,
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
        None,
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

/// DT-020 y DT-021 sobre los **bytes publicados**: el juez de los recursos
/// publicados lee la matriz publicada y la población de cada individuo de la
/// configuración efectiva y pasa el juez completo de la verdad sobre el
/// `truth.jsonl` del paquete, y exige además exactamente las copias de la
/// simulación. Hasta el 3 de octubre de 2026 solo miraba su cabecera, y
/// aceptaba un caso todo-REF inventado, desplazado u omitido, un portador
/// cambiado con su dosis, un caso *indel* quitado, un individuo o una copia
/// quitados y un individuo renombrado.
#[test]
fn published_truth_is_judged_with_the_published_matrix() {
    let fixture_root = temporary("truth-judge-fixture");
    fs::create_dir(&fixture_root).expect("fixture root");
    let config = write_truth_judge_fixture(&fixture_root);
    type Tamper = fn(&mut Vec<Value>);
    let cases: [(&str, Tamper, &str); 9] = [
        (
            "fabricated-all-ref",
            |records| {
                let mut fabricated = records[all_reference_index(records)].clone();
                fabricated["case_id"] = json!(format!("CASE_{}", "f".repeat(59)));
                fabricated["locus"]["start"] = json!(10);
                fabricated["locus"]["end"] = json!(11);
                records.push(fabricated);
            },
            "no corresponde",
        ),
        (
            "shifted-all-ref",
            |records| {
                let index = all_reference_index(records);
                records[index]["locus"]["start"] = json!(10);
                records[index]["locus"]["end"] = json!(11);
            },
            "MissingCase",
        ),
        (
            "omitted-all-ref",
            |records| {
                let index = all_reference_index(records);
                records.remove(index);
            },
            "MissingCase",
        ),
        (
            "flipped-carrier-and-dosage",
            |records| {
                let case = records
                    .iter_mut()
                    .find(|record| {
                        record["case_kind"] == "snp"
                            && record["allele_dosages"][1].as_u64().unwrap_or(0) > 0
                    })
                    .expect("a SNP case with the alternate");
                let index = case["carriers"][0]["allele_index"].as_u64().expect("index");
                case["carriers"][0]["allele_index"] = json!(1 - index);
                recount(case);
            },
            "las copias portadoras no coinciden",
        ),
        (
            "omitted-indel",
            |records| {
                let index = records
                    .iter()
                    .position(|record| record["case_kind"] == "indel")
                    .expect("an indel case");
                records.remove(index);
            },
            "MissingCase",
        ),
        (
            "dropped-individual",
            |records| {
                records.retain(|record| individual_of(record) != "I2");
            },
            "le falta una copia",
        ),
        (
            "dropped-copy-with-rebuilt-carriers",
            |records| {
                records.retain(|record| {
                    !(record["record_type"] == "haplotype_copy"
                        && record["identity"]["individual_id"] == "I1"
                        && record["identity"]["sequence_id"] == "chr1"
                        && record["identity"]["copy_index"] == 1)
                });
                for record in records.iter_mut() {
                    if record["record_type"] == "genotype"
                        && record["individual_id"] == "I1"
                        && record["locus"]["sequence_id"] == "chr1"
                    {
                        record["carriers"]
                            .as_array_mut()
                            .expect("carriers")
                            .retain(|carrier| carrier["copy"]["copy_index"] != 1);
                        recount(record);
                    }
                }
            },
            "le falta una copia",
        ),
        (
            "extra-copy",
            |records| {
                let mut extra = records
                    .iter()
                    .find(|record| {
                        record["record_type"] == "haplotype_copy"
                            && record["identity"]["individual_id"] == "I1"
                            && record["identity"]["sequence_id"] == "chr2"
                            && record["identity"]["copy_index"] == 1
                    })
                    .expect("copy 1 of I1 on chr2")
                    .clone();
                extra["identity"]["copy_index"] = json!(2);
                extra["phase"] = json!(2);
                records.push(extra);
            },
            "trae una copia que la simulación no tiene",
        ),
        (
            "unknown-individual",
            |records| {
                for record in records.iter_mut() {
                    if record["record_type"] == "haplotype_copy"
                        && record["identity"]["individual_id"] == "I2"
                    {
                        record["identity"]["individual_id"] = json!("I9");
                    }
                    if record["record_type"] == "genotype" && record["individual_id"] == "I2" {
                        record["individual_id"] = json!("I9");
                        for carrier in record["carriers"].as_array_mut().expect("carriers") {
                            carrier["copy"]["individual_id"] = json!("I9");
                        }
                    }
                }
            },
            "el individuo no tiene población conocida",
        ),
    ];
    for (label, tamper, reason) in cases {
        let parent = temporary(label);
        fs::create_dir(&parent).expect("temporary root");
        let destination = parent.join("result");
        let effective_config = load_and_validate_config(&config).expect("valid fixture");
        execute_validated_config(&effective_config, &destination).expect("complete baseline");
        let effective_path = destination.join("effective-config.json");
        let manifest_path = destination.join("manifest.json");
        let mut manifest = read_json(&manifest_path);
        let truth_path = destination.join(artifact_path(&manifest, "truth-jsonl"));
        let mut lines = fs::read_to_string(&truth_path)
            .expect("truth")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("truth line"))
            .collect::<Vec<_>>();
        let mut records = lines.split_off(1);
        tamper(&mut records);
        let (mut copies, mut genotypes): (Vec<_>, Vec<_>) = records
            .into_iter()
            .partition(|record| record["record_type"] == "haplotype_copy");
        let identity = |copy: &Value| {
            (
                copy["identity"]["individual_id"]
                    .as_str()
                    .map(str::to_owned),
                copy["identity"]["sequence_id"].as_str().map(str::to_owned),
                copy["identity"]["copy_index"].as_u64(),
            )
        };
        copies.sort_by_key(identity);
        genotypes.sort_by_key(|genotype| {
            (
                genotype["individual_id"].as_str().map(str::to_owned),
                genotype["locus"]["sequence_id"].as_str().map(str::to_owned),
                genotype["locus"]["start"].as_u64(),
                genotype["locus"]["end"].as_u64(),
                genotype["case_id"].as_str().map(str::to_owned),
            )
        });
        lines[0]["counts"]["copies"] = json!(copies.len());
        lines[0]["counts"]["genotypes"] = json!(genotypes.len());
        lines.extend(copies);
        lines.extend(genotypes);
        let text = lines
            .iter()
            .map(|line| serde_json::to_string(line).expect("line") + "\n")
            .collect::<String>();
        fs::write(&truth_path, text).expect("tampered truth");
        refresh_artifact(&destination, &mut manifest, "truth-jsonl");
        write_json(&manifest_path, &manifest);
        let effective_bytes = fs::read(&effective_path).expect("effective bytes");
        let manifest_bytes = fs::read(&manifest_path).expect("manifest bytes");
        let error = validate_published_resources(
            &destination,
            &effective_path,
            &effective_bytes,
            &manifest_path,
            &manifest_bytes,
            None,
        )
        .expect_err(label);
        let debug = format!("{error:?}");
        assert!(debug.contains(reason), "case {label}: {error:?}");
        fs::remove_dir_all(parent).expect("cleanup");
    }
    fs::remove_dir_all(fixture_root).expect("fixture cleanup");
}

fn individual_of(record: &Value) -> &str {
    if record["record_type"] == "haplotype_copy" {
        record["identity"]["individual_id"]
            .as_str()
            .expect("individual")
    } else {
        record["individual_id"].as_str().expect("individual")
    }
}

fn recount(case: &mut Value) {
    let mut dosages = vec![0_u64; case["alleles"].as_array().expect("alleles").len()];
    for carrier in case["carriers"].as_array().expect("carriers") {
        dosages
            [usize::try_from(carrier["allele_index"].as_u64().expect("index")).expect("index")] +=
            1;
    }
    case["allele_dosages"] = json!(dosages);
}

/// El caso todo-REF de `I1` en la fila de frecuencia casi nula.
fn all_reference_index(records: &[Value]) -> usize {
    records
        .iter()
        .position(|record| {
            record["case_kind"] == "snp"
                && record["individual_id"] == "I1"
                && record["locus"]["sequence_id"] == "chr1"
                && record["locus"]["start"] == 7
        })
        .inspect(|index| {
            assert!(
                records[*index]["allele_dosages"][1] == 0,
                "the low-frequency row must stay all-REF: {}",
                records[*index]
            );
        })
        .expect("the all-REF case")
}

/// Dos individuos diploides, dos secuencias, una fila de matriz frecuente y
/// otra de frecuencia 0,01 (caso todo-REF), una deleción directa y, en `I2`,
/// un cambio de base directo en el mismo desplazamiento de cada secuencia, con
/// alternativos distintos: un juez que mezclara las copias de las dos
/// secuencias al catalogar los alelos de un caso rechazaría el paquete bueno.
fn write_truth_judge_fixture(root: &Path) -> PathBuf {
    let fasta = root.join("reference.fa");
    let matrix = root.join("matrix.tsv");
    fs::write(&fasta, b">chr1\nACGTACGTACGT\n>chr2\nACGTACGTACGT\n").expect("reference");
    fs::write(
        &matrix,
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t2\tC\tG\tP1\t0.5\nchr1\t8\tT\tA\tP1\t0.01\n",
    )
    .expect("matrix");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"{}\"\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"{}\"\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nindividual_ids = [\"I1\", \"I2\"]\nploidy = 2\ngenotype_model = \"hwe\"\n\
             [[direct_events]]\nevent_id = \"event-del\"\nindividual_id = \"I1\"\nphase = 0\nkind = \"deletion\"\ncontig = \"chr1\"\nstart = 4\nend = 6\n\
             [[direct_events]]\nevent_id = \"snp-chr1\"\nindividual_id = \"I2\"\nphase = 0\nkind = \"snp\"\ncontig = \"chr1\"\nposition_1based = 11\nref = \"G\"\nalt = \"T\"\n\
             [[direct_events]]\nevent_id = \"snp-chr2\"\nindividual_id = \"I2\"\nphase = 1\nkind = \"snp\"\ncontig = \"chr2\"\nposition_1based = 11\nref = \"G\"\nalt = \"A\"\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"{}\"\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 12\n\
             [[reference.sequences]]\nid = \"chr2\"\nclass = \"chromosome\"\nlength = 12\n",
            "0".repeat(64),
            matrix.display(),
            fasta.display(),
        ),
    )
    .expect("truth judge config");
    config
}

/// Dos loci de repetición en las mismas coordenadas de dos secuencias, con
/// tractos distintos: cada caso cataloga solo los tractos de su secuencia. Un
/// juez que mezclara las copias de las dos al ordenar los alelos de un caso
/// de repetición rechazaría este escenario, que es válido.
#[test]
fn repeat_loci_at_the_same_coordinates_on_two_sequences_publish() {
    let root = temporary("two-sequence-repeats");
    fs::create_dir(&root).expect("root");
    let flank = "GTTGATGTTGATGTTGATGTTG";
    let sequence = format!("GTTGATGTTG{}{flank}", "CA".repeat(4));
    fs::write(
        root.join("reference.fa"),
        format!(">chr1\n{sequence}\n>chr2\n{sequence}\n"),
    )
    .expect("reference");
    fs::write(
        root.join("matrix.tsv"),
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t2\tT\tC\tP1\t0.5\n",
    )
    .expect("matrix");
    let length = sequence.len();
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{}\"\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"matrix.tsv\"\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = {length}\n\
             [[reference.sequences]]\nid = \"chr2\"\nclass = \"chromosome\"\nlength = {length}\n\
             [[reference.repeat_loci]]\nid = \"STR-1\"\nsequence = \"chr1\"\nstart = 11\nend = 18\nmotif = \"CA\"\n\
             [[reference.repeat_loci.alleles]]\nunits = 5\nfrequency = 1.0\n\
             [[reference.repeat_loci]]\nid = \"STR-2\"\nsequence = \"chr2\"\nstart = 11\nend = 18\nmotif = \"CA\"\n\
             [[reference.repeat_loci.alleles]]\nunits = 6\nfrequency = 1.0\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nindividual_ids = [\"I1\"]\nploidy = 2\ngenotype_model = \"hwe\"\n",
            "0".repeat(64),
        ),
    )
    .expect("config");
    let destination = root.join("result");
    let effective = load_and_validate_config(&config).expect("valid config");
    execute_validated_config(&effective, &destination).expect("a valid scenario publishes");
    let truth = fs::read_to_string(destination.join("outputs/truth.jsonl")).expect("truth");
    let repeat_cases = truth
        .lines()
        .filter(|line| line.contains("\"case_kind\":\"repeat\""))
        .count();
    assert_eq!(repeat_cases, 2, "one repeat case per sequence");
    fs::remove_dir_all(root).expect("cleanup");
}

/// DT-027: un paquete **sin poblaciones** también pasa la cabecera de la verdad
/// por el contrato. Hasta el 6 de octubre de 2026 el juez de los recursos
/// publicados no leía en ese caso el documento: solo comparaba nombre, huella
/// de contenido y huellas de fuente, y aceptaba cualquier otra alteración de
/// la cabecera que el lector de la verdad rechaza.
///
/// Se repite con las dos referencias que existen sin poblaciones: la
/// **simulada** y la **suministrada**, que el juez lee solo con su ficha, sin
/// bases. La rama sin bases es la que protege el caso de la huella del
/// contenido (`content-digest-flip` y `content-digest-zeros`), que antes se
/// omitía en silencio.
///
/// Los casos `control-*` (`reference_id`, huella del contenido a ceros y huella
/// de la referencia en mayúsculas) ya los cazaba la comparación de nombres con
/// bases: son el control de que el arnés rechaza por la verdad y no por otra
/// cosa. Lo que esta prueba no discrimina: con bases y sin poblaciones, el
/// juez semántico completo no aporta nada medible sobre el lector.
#[test]
fn truth_header_of_a_package_without_populations_is_judged() {
    type Tamper = fn(&mut Vec<Value>);
    let cases: [(&str, Tamper); 13] = [
        ("control-reference-id", |lines| {
            lines[0]["reference_id"] = json!("ref-other");
        }),
        ("content-digest-zeros", |lines| {
            lines[0]["reference_content_sha256"] = json!("0".repeat(64));
        }),
        ("content-digest-flip", |lines| {
            let digest = lines[0]["reference_content_sha256"]
                .as_str()
                .expect("digest")
                .to_owned();
            let first = if digest.starts_with('0') { "1" } else { "0" };
            lines[0]["reference_content_sha256"] = json!(format!("{first}{}", &digest[1..]));
        }),
        ("document-kind", |lines| {
            lines[0]["document_kind"] = json!("gps_other");
        }),
        ("schema-version", |lines| {
            lines[0]["schema_version"] = json!("3");
        }),
        ("promoted-schema-version", |lines| {
            lines[0]["schema_version"] = json!("5");
        }),
        ("coordinate-convention", |lines| {
            lines[0]["coordinate_convention"] = json!("one_based_closed");
        }),
        ("unknown-header-key", |lines| {
            lines[0]["zz_extra"] = json!(1);
        }),
        ("copy-count", |lines| {
            lines[0]["counts"]["copies"] = json!(1);
        }),
        ("control-uppercase-reference-digest", |lines| {
            let digest = lines[0]["source_digests"][0]["sha256"]
                .as_str()
                .expect("digest")
                .to_uppercase();
            lines[0]["source_digests"][0]["sha256"] = json!(digest);
        }),
        ("malformed-group-digest", |lines| {
            lines[0]["group_digests"] = json!([{
                "sha256": "not-a-digest",
                "source": {"kind": "founder", "source_id": "founder-1"}
            }]);
        }),
        ("unknown-record-type", |lines| {
            lines.push(json!({"record_type": "zz_extra"}));
        }),
        ("second-header", |lines| {
            let header = lines[0].clone();
            lines.push(header);
        }),
    ];
    for (fixture_name, fixture) in [
        ("simulated", simulated_fixture()),
        ("supplied", supplied_fixture()),
    ] {
        for (label, tamper) in cases {
            let label = &format!("{fixture_name}-{label}");
            let parent = temporary(label);
            fs::create_dir(&parent).expect("temporary root");
            let destination = parent.join("result");
            let effective_config = load_and_validate_config(&fixture).expect("fixture");
            execute_validated_config(&effective_config, &destination).expect("complete baseline");
            let effective_path = destination.join("effective-config.json");
            let manifest_path = destination.join("manifest.json");
            let mut manifest = read_json(&manifest_path);
            let truth_path = destination.join(artifact_path(&manifest, "truth-jsonl"));
            let text = fs::read_to_string(&truth_path).expect("truth");
            let mut lines = text
                .lines()
                .map(|line| serde_json::from_str::<Value>(line).expect("truth line"))
                .collect::<Vec<_>>();
            assert_eq!(
                lines.len(),
                1,
                "a package without populations has no records"
            );
            tamper(&mut lines);
            let tampered = lines
                .iter()
                .map(|line| serde_json::to_string(line).expect("line") + "\n")
                .collect::<String>();
            assert_ne!(
                tampered, text,
                "case {label}: the tamper must change the truth"
            );
            fs::write(&truth_path, tampered).expect("tampered truth");
            refresh_artifact(&destination, &mut manifest, "truth-jsonl");
            write_json(&manifest_path, &manifest);
            let effective_bytes = fs::read(&effective_path).expect("effective bytes");
            let manifest_bytes = fs::read(&manifest_path).expect("manifest bytes");
            let error = validate_published_resources(
                &destination,
                &effective_path,
                &effective_bytes,
                &manifest_path,
                &manifest_bytes,
                None,
            )
            .expect_err(label);
            assert_eq!(error.code(), "GPS-FMT-006", "case {label}: {error:?}");
            fs::remove_dir_all(parent).expect("cleanup");
        }
    }
}
