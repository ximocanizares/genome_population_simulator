use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::{ExecutionOutcome, execute_validated_config};
use genome_population_simulator::reference::{Reference, ReferenceSequence, SequenceClass};
use genome_population_simulator::reference_fasta::{SequenceExpectation, validate_supplied_fasta};
use genome_population_simulator::reference_reduction_io::{
    ReductionMapParameters, validate_reduction_map_document,
};
use genome_population_simulator::results::UnitState;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/reference_reduction")
        .join(name)
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f1-6-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn write_config_with_reductions(
    root: &Path,
    file_name: &str,
    root_fields: &str,
    reductions: &str,
) -> PathBuf {
    fs::create_dir_all(root).unwrap();
    fs::copy(fixture("source.fa"), root.join("source.fa")).unwrap();
    let original = fs::read_to_string(fixture("general.toml")).unwrap();
    let base = original
        .split("[[reduced_references]]")
        .next()
        .expect("fixture base");
    let base = base.replace("\n[reference]", &format!("\n{root_fields}\n[reference]"));
    let path = root.join(file_name);
    fs::write(&path, format!("{base}\n{reductions}\n")).unwrap();
    path
}

fn read_fasta_independently(path: &Path) -> BTreeMap<String, Vec<u8>> {
    let text = fs::read_to_string(path).expect("FASTA text");
    let mut sequences = BTreeMap::new();
    let mut current = None::<String>;
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let id = header
                .split_whitespace()
                .next()
                .expect("FASTA ID")
                .to_owned();
            sequences.insert(id.clone(), Vec::new());
            current = Some(id);
        } else {
            sequences
                .get_mut(current.as_deref().expect("sequence before bases"))
                .expect("known sequence")
                .extend_from_slice(line.as_bytes());
        }
    }
    sequences
}

fn read_fasta_order_independently(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .expect("FASTA text")
        .lines()
        .filter_map(|line| {
            line.strip_prefix('>').map(|header| {
                header
                    .split_whitespace()
                    .next()
                    .expect("FASTA ID")
                    .to_owned()
            })
        })
        .collect()
}

fn reconstruct_canonical_fasta_independently(
    source_path: &Path,
    map: &serde_json::Value,
    sequence_classes: &[(String, &'static str)],
    derived_catalog: &[(String, &'static str)],
    region_size_bases: u64,
    target_bases: u64,
) -> Result<Vec<u8>, String> {
    let source = read_fasta_independently(source_path);
    let segments = map
        .get("segments")
        .and_then(serde_json::Value::as_array)
        .ok_or("segments is not an array")?;
    let full_count = target_bases / region_size_bases;
    let residual = target_bases % region_size_bases;
    let expected_count = full_count + u64::from(residual > 0);
    if u64::try_from(segments.len()).ok() != Some(expected_count) {
        return Err("segment count differs from effective parameters".to_owned());
    }
    let mut by_sequence = BTreeMap::<String, Vec<u8>>::new();
    let mut previous_source_end = BTreeMap::<String, u64>::new();
    let mut previous_derived_end = BTreeMap::<String, u64>::new();
    let mut selection_orders = Vec::new();
    let mut materialized_catalog = Vec::<(String, &'static str)>::new();
    let mut previous_sequence_index = None;
    for (index, segment) in segments.iter().enumerate() {
        if segment["materialization_order"].as_u64()
            != u64::try_from(index)
                .ok()
                .and_then(|value| value.checked_add(1))
        {
            return Err("materialization order is not contiguous".to_owned());
        }
        if segment["orientation"] != "forward" {
            return Err("orientation is not forward".to_owned());
        }
        let selection_order = segment["selection_order"]
            .as_u64()
            .ok_or("selection order is not u64")?;
        selection_orders.push(selection_order);
        let sequence_id = segment["sequence_id"]
            .as_str()
            .ok_or("sequence ID is not text")?;
        let effective_class = sequence_classes
            .iter()
            .find(|(id, _)| id == sequence_id)
            .map(|(_, class)| *class);
        if segment["sequence_class"].as_str() != effective_class {
            return Err("sequence class differs from the effective catalog".to_owned());
        }
        let sequence_index = sequence_classes
            .iter()
            .position(|(id, _)| id == sequence_id)
            .ok_or("sequence absent from source catalog")?;
        if previous_sequence_index.is_some_and(|previous| sequence_index < previous) {
            return Err("materialized sequences invert source catalog order".to_owned());
        }
        previous_sequence_index = Some(sequence_index);
        if materialized_catalog
            .last()
            .is_none_or(|(id, _)| id != sequence_id)
        {
            materialized_catalog.push((
                sequence_id.to_owned(),
                effective_class.ok_or("missing effective sequence class")?,
            ));
        }
        let source_start = segment["source_interval"]["start"]
            .as_u64()
            .ok_or("source start is not u64")?;
        let source_end = segment["source_interval"]["end"]
            .as_u64()
            .ok_or("source end is not u64")?;
        let derived_start = segment["derived_interval"]["start"]
            .as_u64()
            .ok_or("derived start is not u64")?;
        let derived_end = segment["derived_interval"]["end"]
            .as_u64()
            .ok_or("derived end is not u64")?;
        if source_start == 0 || source_start > source_end || derived_start == 0 {
            return Err("invalid closed interval".to_owned());
        }
        if previous_source_end
            .get(sequence_id)
            .is_some_and(|end| source_start <= *end)
        {
            return Err("source overlap or order inversion".to_owned());
        }
        let expected_derived_start =
            previous_derived_end.get(sequence_id).copied().unwrap_or(0) + 1;
        if derived_start != expected_derived_start
            || derived_end - derived_start != source_end - source_start
        {
            return Err("derived gap, overlap or length change".to_owned());
        }
        let (expected_class, expected_length) = if selection_order <= full_count {
            ("full", region_size_bases)
        } else if residual > 0 && selection_order == expected_count {
            ("residual", residual)
        } else {
            return Err("selection order contradicts effective parameters".to_owned());
        };
        if segment["segment_class"] != expected_class
            || source_end - source_start + 1 != expected_length
        {
            return Err("full/residual label contradicts effective parameters".to_owned());
        }
        let source_bases = source.get(sequence_id).ok_or("unknown source sequence")?;
        let start = usize::try_from(source_start - 1).map_err(|_| "start does not fit")?;
        let end = usize::try_from(source_end).map_err(|_| "end does not fit")?;
        let fragment = source_bases
            .get(start..end)
            .ok_or("source interval out of bounds")?;
        by_sequence
            .entry(sequence_id.to_owned())
            .or_default()
            .extend(fragment.iter().map(u8::to_ascii_uppercase));
        previous_source_end.insert(sequence_id.to_owned(), source_end);
        previous_derived_end.insert(sequence_id.to_owned(), derived_end);
    }
    selection_orders.sort_unstable();
    if selection_orders != (1..=expected_count).collect::<Vec<_>>() {
        return Err("selection order is incomplete".to_owned());
    }
    if by_sequence.values().map(Vec::len).sum::<usize>()
        != usize::try_from(target_bases).map_err(|_| "target does not fit")?
    {
        return Err("rebuilt target differs".to_owned());
    }
    if materialized_catalog != derived_catalog {
        return Err("derived catalog order or class differs from the map".to_owned());
    }
    let mut fasta = Vec::new();
    for (sequence_id, _) in sequence_classes {
        let Some(bases) = by_sequence.get(sequence_id) else {
            continue;
        };
        fasta.extend_from_slice(b">");
        fasta.extend_from_slice(sequence_id.as_bytes());
        fasta.push(b'\n');
        for line in bases.chunks(80) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }
    }
    Ok(fasta)
}

fn validate_json(instance: &Path, schema: &str, expected_success: bool) {
    let schema_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("schemas");
    let base_uri = format!("file://{}/", schema_root.display());
    let output = Command::new("python3")
        .args(["-m", "jsonschema", "--base-uri", &base_uri, "-i"])
        .arg(instance)
        .arg(schema)
        .current_dir(schema_root)
        .output()
        .expect("independent JSON Schema validator");
    assert_eq!(
        output.status.success(),
        expected_success,
        "schema={schema} stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_artifacts_are_classified_once(
    manifest: &serde_json::Value,
    expected_reduction_units: usize,
    expected_explicit_artifacts: usize,
) {
    let outputs = manifest["outputs"].as_object().expect("outputs object");
    let mut classified = BTreeMap::<String, usize>::new();
    for category in ["minimum", "default_enabled", "explicitly_requested"] {
        for artifact_id in outputs[category].as_array().expect("output category") {
            *classified
                .entry(artifact_id.as_str().expect("artifact ID").to_owned())
                .or_default() += 1;
        }
    }
    assert_eq!(classified.get("manifest-general"), Some(&1));
    for artifact in manifest["artifacts"].as_array().expect("artifacts") {
        let id = artifact["id"].as_str().expect("artifact ID");
        assert_eq!(classified.get(id), Some(&1), "classification for {id}");
    }
    assert_eq!(
        outputs["explicitly_requested"]
            .as_array()
            .expect("explicit outputs")
            .len(),
        expected_explicit_artifacts
    );
    assert_eq!(
        manifest["units"]
            .as_array()
            .expect("units")
            .iter()
            .filter(|unit| unit["kind"] == "reference_reduction")
            .count(),
        expected_reduction_units
    );
}

fn load_map_references(
    root: &Path,
) -> (
    genome_population_simulator::reference::Reference,
    genome_population_simulator::reference::Reference,
) {
    let source_path = root.join("resources/references/reference-main.fa");
    let source_expectations = [
        SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 12,
        },
        SequenceExpectation {
            id: "scaffold-1",
            class: SequenceClass::Scaffold,
            length: 8,
        },
    ];
    let source = validate_supplied_fasta(&source_path, "reference-main", &source_expectations)
        .expect("independent source load")
        .into_reference();
    let derived_path = root.join("resources/references/reference-small.fa");
    let derived_bytes = read_fasta_independently(&derived_path);
    let derived_expectations = source
        .sequences()
        .iter()
        .filter_map(|sequence| {
            derived_bytes
                .get(sequence.id().as_str())
                .map(|bases| SequenceExpectation {
                    id: sequence.id().as_str(),
                    class: sequence.class(),
                    length: u64::try_from(bases.len()).expect("small fixture"),
                })
        })
        .collect::<Vec<_>>();
    let derived = validate_supplied_fasta(&derived_path, "reference-small", &derived_expectations)
        .expect("independent derived load")
        .into_reference();
    (source, derived)
}

#[test]
fn target_reduction_publishes_reconstructible_fasta_map_and_manifest() {
    let effective = load_and_validate_config(&fixture("general.toml")).expect("valid config");
    let first = temporary("first");
    let second = temporary("second");
    for destination in [&first, &second] {
        assert_eq!(
            execute_validated_config(&effective, destination).expect("complete run"),
            ExecutionOutcome::Published {
                destination: destination.clone(),
                state: UnitState::Complete,
            }
        );
    }

    let source = read_fasta_independently(&first.join("resources/references/reference-main.fa"));
    let derived = read_fasta_independently(&first.join("resources/references/reference-small.fa"));
    let map: serde_json::Value = serde_json::from_slice(
        &fs::read(first.join("resources/references/reference-small.map.json")).expect("map"),
    )
    .expect("map JSON");
    assert_eq!(map["document_kind"], "reference_reduction_map");
    assert_eq!(map["schema_version"], "1");
    assert_eq!(map["algorithm"], "reduced_reference_selection_v1");
    let sequence_classes = effective
        .sequences()
        .iter()
        .map(|sequence| {
            (
                sequence.id().to_owned(),
                match sequence.class() {
                    SequenceClass::Chromosome => "chromosome",
                    SequenceClass::Scaffold => "scaffold",
                },
            )
        })
        .collect::<Vec<_>>();
    let derived_catalog =
        read_fasta_order_independently(&first.join("resources/references/reference-small.fa"))
            .into_iter()
            .map(|id| {
                let class = sequence_classes
                    .iter()
                    .find(|(source_id, _)| source_id == &id)
                    .map(|(_, class)| *class)
                    .expect("derived ID exists in effective catalog");
                (id, class)
            })
            .collect::<Vec<_>>();
    let reconstructed_fasta = reconstruct_canonical_fasta_independently(
        &first.join("resources/references/reference-main.fa"),
        &map,
        &sequence_classes,
        &derived_catalog,
        effective.reduced_references()[0].region_size_bases(),
        effective.reduced_references()[0].target_bases(),
    )
    .expect("independent canonical reconstruction");
    assert_eq!(
        reconstructed_fasta,
        fs::read(first.join("resources/references/reference-small.fa")).unwrap()
    );

    let mut rebuilt = BTreeMap::<String, Vec<u8>>::new();
    let mut previous = None;
    let mut residuals = 0;
    for segment in map["segments"].as_array().expect("segments") {
        let order = segment["materialization_order"].as_u64().expect("order");
        assert_eq!(order, previous.map_or(1, |value| value + 1));
        previous = Some(order);
        residuals += usize::from(segment["segment_class"] == "residual");
        assert_eq!(segment["orientation"], "forward");
        let id = segment["sequence_id"].as_str().expect("sequence ID");
        let start = usize::try_from(segment["source_interval"]["start"].as_u64().unwrap()).unwrap();
        let end = usize::try_from(segment["source_interval"]["end"].as_u64().unwrap()).unwrap();
        rebuilt
            .entry(id.to_owned())
            .or_default()
            .extend_from_slice(&source[id][start - 1..end]);
    }
    assert!(residuals <= 1);
    assert_eq!(rebuilt, derived);
    assert_eq!(derived.values().map(Vec::len).sum::<usize>(), 13);
    assert!(derived.values().flatten().any(|base| *base == b'N'));

    let effective_json: serde_json::Value = serde_json::from_slice(
        &fs::read(first.join("effective-config.json")).expect("effective config"),
    )
    .unwrap();
    let reduction = &effective_json["reduced_references"][0];
    assert_eq!(reduction["target"]["declared_value"], "13");
    assert_eq!(reduction["target"]["derived_bases"], 13);
    assert_eq!(reduction["eligible_bases"], 20);
    assert_eq!(
        reduction["randomness"]["subseed"],
        effective.reduced_references()[0]
            .selection_subseed()
            .to_hex()
    );

    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(first.join("manifest.json")).expect("manifest")).unwrap();
    assert_eq!(manifest["references"].as_array().unwrap().len(), 2);
    assert_eq!(
        manifest["reference_reductions"].as_array().unwrap().len(),
        1
    );
    assert_eq!(
        manifest["reference_reductions"][0]["map"]["validated"],
        true
    );
    validate_json(
        &first.join("effective-config.json"),
        "effective-config-v1.schema.json",
        true,
    );
    validate_json(
        &first.join("manifest.json"),
        "manifest-general-v1.schema.json",
        true,
    );
    validate_json(
        &first.join("resources/references/reference-small.map.json"),
        "reference-reduction-map-v1.schema.json",
        true,
    );
    let derived_path = first.join("resources/references/reference-small.fa");
    let sha = Command::new("sha256sum")
        .arg(&derived_path)
        .output()
        .expect("sha256sum");
    assert!(sha.status.success());
    let observed_digest = String::from_utf8(sha.stdout)
        .unwrap()
        .split_whitespace()
        .next()
        .unwrap()
        .to_owned();
    assert_eq!(
        observed_digest,
        manifest["reference_reductions"][0]["derived_fasta"]["checksum"]["digest"]
            .as_str()
            .unwrap()
    );
    let reconstructed_path = temporary("reconstructed.fa");
    fs::write(&reconstructed_path, &reconstructed_fasta).unwrap();
    let reconstructed_sha = Command::new("sha256sum")
        .arg(&reconstructed_path)
        .output()
        .expect("sha256sum reconstructed FASTA");
    assert!(reconstructed_sha.status.success());
    assert_eq!(
        String::from_utf8(reconstructed_sha.stdout)
            .unwrap()
            .split_whitespace()
            .next(),
        Some(observed_digest.as_str())
    );
    let mut altered_reconstruction = reconstructed_fasta.clone();
    let base = altered_reconstruction
        .iter_mut()
        .find(|byte| matches!(**byte, b'A' | b'C' | b'G' | b'T' | b'N'))
        .expect("reconstructed base");
    *base = if *base == b'A' { b'C' } else { b'A' };
    let altered_path = temporary("reconstructed-altered.fa");
    fs::write(&altered_path, altered_reconstruction).unwrap();
    let altered_sha = Command::new("sha256sum")
        .arg(&altered_path)
        .output()
        .expect("sha256sum altered reconstruction");
    assert!(altered_sha.status.success());
    assert_ne!(
        String::from_utf8(altered_sha.stdout)
            .unwrap()
            .split_whitespace()
            .next(),
        Some(observed_digest.as_str())
    );
    assert_eq!(
        fs::read(fixture("source.fa")).unwrap(),
        fs::read(first.join("resources/references/reference-main.fa")).unwrap()
    );
    let faidx_copy = temporary("faidx.fa");
    fs::copy(&derived_path, &faidx_copy).unwrap();
    let faidx = Command::new("samtools")
        .arg("faidx")
        .arg(&faidx_copy)
        .status()
        .expect("samtools faidx");
    assert!(faidx.success());
    assert!(PathBuf::from(format!("{}.fai", faidx_copy.display())).is_file());
    assert_eq!(
        fs::read(first.join("resources/references/reference-small.fa")).unwrap(),
        fs::read(second.join("resources/references/reference-small.fa")).unwrap()
    );
    assert_eq!(
        fs::read(first.join("resources/references/reference-small.map.json")).unwrap(),
        fs::read(second.join("resources/references/reference-small.map.json")).unwrap()
    );

    fs::remove_dir_all(first).unwrap();
    fs::remove_dir_all(second).unwrap();
    fs::remove_file(&faidx_copy).unwrap();
    fs::remove_file(format!("{}.fai", faidx_copy.display())).unwrap();
    fs::remove_file(reconstructed_path).unwrap();
    fs::remove_file(altered_path).unwrap();
}

#[test]
fn serialized_map_validator_rejects_truncation_overlap_and_out_of_range() {
    let effective = load_and_validate_config(&fixture("general.toml")).expect("valid config");
    let destination = temporary("map-corruption");
    execute_validated_config(&effective, &destination).expect("complete run");
    let map_path = destination.join("resources/references/reference-small.map.json");
    let bytes = fs::read(&map_path).unwrap();
    let map: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
    let (source, derived) = load_map_references(&destination);
    let parameters = ReductionMapParameters {
        region_size_bases: 5,
        target_bases: 13,
    };
    validate_reduction_map_document(&source, &derived, &map, parameters)
        .expect("valid serialized map");

    let mut reordered_sequences = derived.sequences().to_vec();
    reordered_sequences.swap(0, 1);
    let reordered_derived = Reference::new(
        derived.id().clone(),
        derived.origin().clone(),
        reordered_sequences,
    )
    .expect("reordered derived fixture remains a valid reference");
    assert!(
        validate_reduction_map_document(&source, &reordered_derived, &map, parameters).is_err()
    );

    let first = &derived.sequences()[0];
    let changed_class = match first.class() {
        SequenceClass::Chromosome => SequenceClass::Scaffold,
        SequenceClass::Scaffold => SequenceClass::Chromosome,
    };
    let replacement = ReferenceSequence::new(
        first.id().clone(),
        first.header().clone(),
        changed_class,
        first.sequence(),
    )
    .expect("class-mutated derived sequence");
    let mut changed_class_sequences = derived.sequences().to_vec();
    changed_class_sequences[0] = replacement;
    let changed_class_derived = Reference::new(
        derived.id().clone(),
        derived.origin().clone(),
        changed_class_sequences,
    )
    .expect("class-mutated derived fixture remains structurally valid");
    assert!(
        validate_reduction_map_document(&source, &changed_class_derived, &map, parameters).is_err()
    );

    assert!(serde_json::from_slice::<serde_json::Value>(&bytes[..bytes.len() / 2]).is_err());

    let mut out_of_range = map.clone();
    out_of_range["segments"][0]["source_interval"]["end"] = serde_json::json!(999);
    assert!(validate_reduction_map_document(&source, &derived, &out_of_range, parameters).is_err());

    let mut unknown_root = map.clone();
    unknown_root
        .as_object_mut()
        .expect("root object")
        .insert("unknown_root".to_owned(), serde_json::json!(true));
    assert!(validate_reduction_map_document(&source, &derived, &unknown_root, parameters).is_err());
    let mut unknown_segment = map.clone();
    unknown_segment["segments"][0]
        .as_object_mut()
        .unwrap()
        .insert("unknown_segment".to_owned(), serde_json::json!(true));
    assert!(
        validate_reduction_map_document(&source, &derived, &unknown_segment, parameters).is_err()
    );
    let mut unknown_interval = map.clone();
    unknown_interval["segments"][0]["source_interval"]
        .as_object_mut()
        .unwrap()
        .insert("unknown_interval".to_owned(), serde_json::json!(true));
    assert!(
        validate_reduction_map_document(&source, &derived, &unknown_interval, parameters).is_err()
    );
    for (name, document) in [
        ("unknown-root", &unknown_root),
        ("unknown-segment", &unknown_segment),
        ("unknown-interval", &unknown_interval),
    ] {
        let path = temporary(&format!("{name}.map.json"));
        fs::write(&path, serde_json::to_vec(document).unwrap()).unwrap();
        validate_json(&path, "reference-reduction-map-v1.schema.json", false);
        fs::remove_file(path).unwrap();
    }

    let full_index = map["segments"]
        .as_array()
        .unwrap()
        .iter()
        .position(|segment| segment["segment_class"] == "full")
        .unwrap();
    let residual_index = map["segments"]
        .as_array()
        .unwrap()
        .iter()
        .position(|segment| segment["segment_class"] == "residual")
        .unwrap();
    let mut swapped_labels = map.clone();
    swapped_labels["segments"][full_index]["segment_class"] = serde_json::json!("residual");
    swapped_labels["segments"][residual_index]["segment_class"] = serde_json::json!("full");
    assert!(
        validate_reduction_map_document(&source, &derived, &swapped_labels, parameters).is_err()
    );

    let mut wrong_order = map.clone();
    wrong_order["segments"][0]["materialization_order"] = serde_json::json!(2);
    assert!(validate_reduction_map_document(&source, &derived, &wrong_order, parameters).is_err());

    let mut derived_gap = map.clone();
    let start = derived_gap["segments"][0]["derived_interval"]["start"]
        .as_u64()
        .unwrap();
    derived_gap["segments"][0]["derived_interval"]["start"] = serde_json::json!(start + 1);
    assert!(validate_reduction_map_document(&source, &derived, &derived_gap, parameters).is_err());

    let mut derived_overlap = map.clone();
    let end = derived_overlap["segments"][0]["derived_interval"]["end"]
        .as_u64()
        .unwrap();
    derived_overlap["segments"][0]["derived_interval"]["end"] = serde_json::json!(end + 1);
    assert!(
        validate_reduction_map_document(&source, &derived, &derived_overlap, parameters).is_err()
    );

    let segments = map["segments"].as_array().unwrap();
    let mut pair = None;
    for first in 0..segments.len() {
        for second in first + 1..segments.len() {
            if segments[first]["sequence_id"] == segments[second]["sequence_id"] {
                pair = Some((first, second));
            }
        }
    }
    let (first, second) = pair.expect("three segments across two sequences contain a pair");
    let mut overlap = map.clone();
    overlap["segments"][second]["source_interval"] =
        overlap["segments"][first]["source_interval"].clone();
    assert!(validate_reduction_map_document(&source, &derived, &overlap, parameters).is_err());

    let mut derived_overlap_between_segments = map.clone();
    let previous_end =
        derived_overlap_between_segments["segments"][first]["derived_interval"]["end"]
            .as_u64()
            .unwrap();
    let later_start =
        derived_overlap_between_segments["segments"][second]["derived_interval"]["start"]
            .as_u64()
            .unwrap();
    let later_end = derived_overlap_between_segments["segments"][second]["derived_interval"]["end"]
        .as_u64()
        .unwrap();
    let later_length = later_end - later_start + 1;
    derived_overlap_between_segments["segments"][second]["derived_interval"]["start"] =
        serde_json::json!(previous_end);
    derived_overlap_between_segments["segments"][second]["derived_interval"]["end"] =
        serde_json::json!(previous_end + later_length - 1);
    assert!(
        validate_reduction_map_document(
            &source,
            &derived,
            &derived_overlap_between_segments,
            parameters
        )
        .is_err()
    );

    fs::remove_dir_all(destination).unwrap();
}

#[test]
fn percentage_mode_is_exact_and_infeasible_layout_never_creates_a_derived_artifact() {
    let root = temporary("percentage-root");
    fs::create_dir_all(&root).unwrap();
    fs::copy(fixture("source.fa"), root.join("source.fa")).unwrap();
    let original = fs::read_to_string(fixture("general.toml")).unwrap();
    let percentage = original.replace("target_bases = 13", "reduction_percent = 3.5e1");
    let percentage_config = root.join("percentage.toml");
    fs::write(&percentage_config, percentage).unwrap();
    let effective = load_and_validate_config(&percentage_config).expect("percentage config");
    assert_eq!(effective.reduced_references()[0].target_bases(), 13);
    let percentage_output = root.join("percentage-output");
    execute_validated_config(&effective, &percentage_output).expect("percentage run");
    assert_eq!(
        read_fasta_independently(
            &percentage_output.join("resources/references/reference-small.fa")
        )
        .values()
        .map(Vec::len)
        .sum::<usize>(),
        13
    );
    let alternate_seed = percentage_config.with_file_name("alternate-seed.toml");
    fs::write(
        &alternate_seed,
        fs::read_to_string(&percentage_config).unwrap().replace(
            "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
        ),
    )
    .unwrap();
    let alternate = load_and_validate_config(&alternate_seed).expect("alternate seed config");
    let alternate_output = root.join("alternate-output");
    execute_validated_config(&alternate, &alternate_output).expect("alternate seed run");
    assert_ne!(
        fs::read(percentage_output.join("resources/references/reference-small.map.json")).unwrap(),
        fs::read(alternate_output.join("resources/references/reference-small.map.json")).unwrap()
    );

    let infeasible = original
        .replace("region_size_bases = 5", "region_size_bases = 10")
        .replace("target_bases = 13", "target_bases = 19");
    let infeasible_config = root.join("infeasible.toml");
    fs::write(&infeasible_config, infeasible).unwrap();
    let effective = load_and_validate_config(&infeasible_config).expect("valid declaration");
    let destination = root.join("infeasible-output");
    let diagnostic = execute_validated_config(&effective, &destination).unwrap_err();
    assert_eq!(diagnostic.code(), "GPS-CONF-022");
    assert!(!destination.exists());
    for entry in fs::read_dir(&root).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() && path != percentage_output && path != alternate_output {
            assert!(
                !path
                    .join("resources/references/reference-small.fa")
                    .exists()
            );
            assert!(
                !path
                    .join("resources/references/reference-small.map.json")
                    .exists()
            );
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn output_categories_and_reduction_units_cover_zero_one_and_multiple_reductions() {
    let first = r#"[[reduced_references]]
id = "reference-small"
source_reference_id = "reference-main"
region_size_bases = 5
target_bases = 13"#;
    let second = r#"[[reduced_references]]
id = "reference-tiny"
source_reference_id = "reference-main"
region_size_bases = 4
target_bases = 7"#;
    for (count, reductions) in [
        (0, "".to_owned()),
        (1, first.to_owned()),
        (2, format!("{first}\n\n{second}")),
    ] {
        let root = temporary(&format!("classification-{count}"));
        let config = write_config_with_reductions(&root, "general.toml", "", &reductions);
        let effective = load_and_validate_config(&config).expect("classification config");
        let destination = root.join("output");
        execute_validated_config(&effective, &destination).expect("complete classified run");
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(destination.join("manifest.json")).unwrap()).unwrap();
        assert_artifacts_are_classified_once(&manifest, count, count * 2);
        assert!(
            manifest["units"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|unit| unit["kind"] == "reference_reduction")
                .all(|unit| unit["state"] == "complete" && unit["exportable"] == true)
        );
        validate_json(
            &destination.join("manifest.json"),
            "manifest-general-v1.schema.json",
            true,
        );
        fs::remove_dir_all(root).unwrap();
    }
}

#[test]
fn omitted_sequence_and_unrelated_consumers_do_not_change_the_reduction_stream() {
    let reduction = r#"[[reduced_references]]
id = "reference-one-block"
source_reference_id = "reference-main"
region_size_bases = 5
target_bases = 5"#;
    let baseline_root = temporary("isolation-baseline");
    let baseline_config =
        write_config_with_reductions(&baseline_root, "general.toml", "", reduction);
    let baseline = load_and_validate_config(&baseline_config).expect("baseline config");
    let baseline_output = baseline_root.join("output");
    execute_validated_config(&baseline, &baseline_output).expect("baseline run");

    let consumers_root = temporary("isolation-consumers");
    let consumers_config = write_config_with_reductions(
        &consumers_root,
        "general.toml",
        "individuals = [{ id = \"individual-extra\" }]\nevents = [{ id = \"event-extra\" }]\n\n[[rules]]\nid = \"output-extra\"\nparameter = \"output_enabled\"\nvalue = true",
        reduction,
    );
    let consumers = load_and_validate_config(&consumers_config).expect("consumer config");
    assert!(
        consumers
            .values()
            .iter()
            .any(|value| value.parameter().as_str() == "output_enabled")
    );
    let consumers_output = consumers_root.join("output");
    execute_validated_config(&consumers, &consumers_output).expect("consumer run");

    assert_eq!(
        baseline.reduced_references()[0].selection_subseed(),
        consumers.reduced_references()[0].selection_subseed()
    );
    assert_eq!(
        fs::read(baseline_output.join("resources/references/reference-one-block.map.json"))
            .unwrap(),
        fs::read(consumers_output.join("resources/references/reference-one-block.map.json"))
            .unwrap()
    );
    assert_eq!(
        fs::read(baseline_output.join("resources/references/reference-one-block.fa")).unwrap(),
        fs::read(consumers_output.join("resources/references/reference-one-block.fa")).unwrap()
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(baseline_output.join("manifest.json")).unwrap()).unwrap();
    let reduction_record = &manifest["reference_reductions"][0];
    assert_eq!(
        reduction_record["included_sequence_ids"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        reduction_record["omitted_sequence_ids"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_ne!(
        manifest["reference_reductions"][0]["included_sequence_ids"][0],
        manifest["reference_reductions"][0]["omitted_sequence_ids"][0]
    );

    fs::remove_dir_all(baseline_root).unwrap();
    fs::remove_dir_all(consumers_root).unwrap();
}
