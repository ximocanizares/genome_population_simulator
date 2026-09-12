use std::collections::BTreeSet;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{
    EffectiveF5Config, load_and_validate_f5_config, load_and_validate_f5_config_with_entropy,
    load_and_validate_f5_config_with_entropy_and_resources,
};
use genome_population_simulator::f5_layout::{
    F5BlockLengthPlan, F5LayoutRequest, build_f5_physical_layout,
};
use genome_population_simulator::f5_preflight::SealedF5Admission;
use genome_population_simulator::f5_preflight::{
    F5_HARD_OPERATION_LIMIT, F5_HARD_OUTPUT_LIMIT, F5_HARD_WORKING_MEMORY_LIMIT, F5OtherResources,
    recover_sealed_f5_admission,
};
use genome_population_simulator::identity::{
    ExecutionId, F5IdentityContext, F5IdentityRegistry, ReferenceId, ScenarioId,
};
use genome_population_simulator::randomness::{
    EntropySource, F5RandomDecision, RootSeed, derive_f5_subseed,
};
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

/// Recupera la admisión sellada con el nombre que estas pruebas usaban antes
/// del sellado de F5.1.2.5. Ningún consumidor productivo la conserva: layout y
/// contrato exigen ya `SealedF5Admission`.
fn preflight_f5_resources(
    config: &EffectiveF5Config,
    other: F5OtherResources,
) -> Result<SealedF5Admission, genome_population_simulator::diagnostics::Diagnostic> {
    recover_sealed_f5_admission(config, other)
}

#[derive(Debug)]
struct FixedEntropy([u8; 32]);

impl EntropySource for FixedEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        *destination = self.0;
        Ok(())
    }
}

#[derive(Debug, Default)]
struct CountingEntropy {
    calls: usize,
}

impl EntropySource for CountingEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        self.calls += 1;
        *destination = [0xbb; 32];
        Ok(())
    }
}

fn load_with_resources(
    path: &Path,
    other: F5OtherResources,
) -> Result<EffectiveF5Config, genome_population_simulator::diagnostics::Diagnostic> {
    load_and_validate_f5_config_with_entropy_and_resources(
        path,
        &mut FixedEntropy([0xaa; 32]),
        other,
    )
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-preflight-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn assert_preflight_projection_validates_through_effective_and_manifest(
    root: &Path,
    projection: &serde_json::Value,
) {
    let bundle = root.join("schema-projection");
    fs::create_dir_all(&bundle).expect("schema projection directory");
    let repository = Path::new(env!("CARGO_MANIFEST_DIR"));
    for schema in [
        "effective-config-v5.schema.json",
        "manifest-general-v5.schema.json",
        "f5-failure-v1.schema.json",
    ] {
        fs::copy(repository.join("schemas").join(schema), bundle.join(schema))
            .expect("schema dependency");
    }
    let instance = bundle.join("preflight.json");
    fs::write(
        &instance,
        serde_json::to_vec(projection).expect("serialized public preflight projection"),
    )
    .expect("preflight instance");
    for (wrapper, contract_ref) in [
        (
            "effective-preflight.schema.json",
            "effective-config-v5.schema.json#/$defs/preflight",
        ),
        (
            "manifest-preflight.schema.json",
            "manifest-general-v5.schema.json#/$defs/success/properties/resource_preflight",
        ),
    ] {
        let wrapper_path = bundle.join(wrapper);
        fs::write(
            &wrapper_path,
            format!(
                r#"{{"$schema":"https://json-schema.org/draft/2020-12/schema","$ref":"{contract_ref}"}}"#
            ),
        )
        .expect("wrapper schema");
        let output = Command::new("python3")
            .arg(repository.join("tests/external_jsonschema_validator.py"))
            .arg(&wrapper_path)
            .arg(&instance)
            .output()
            .expect("external JSON Schema validator");
        assert!(
            output.status.success(),
            "{wrapper}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

fn assert_f5_1_2_5_v8_contract_key_set(projection: &serde_json::Value) {
    // Oráculo literal independiente de `to_contract_json` y de los esquemas.
    // F5.1.2.5 cierra el objeto con `input_identities` y `phases` junto a los
    // escalares ya materializados en F5.1.2.3.
    let expected = [
        "algorithm",
        "input_identities",
        "phases",
        "reference_bases_processed",
        "fasta_bytes_processed",
        "mask_bytes_processed",
        "general_toml_bytes",
        "scenario_toml_bytes",
        "profile_metadata_bytes",
        "profile_tsv_bytes",
        "tract_count",
        "retained_sequence_ids",
        "gap_support_cardinality",
        "gap_selector_operation_units",
        "site_bound_per_block",
        "bin_count_total",
        "profile_rows",
        "normalization_rows",
        "tsv_rows_bound",
        "json_rows_bound",
        "other_operation_units",
        "other_working_memory_bytes",
        "other_output_bytes",
        "operation_units",
        "working_memory_bytes",
        "output_bytes",
        "operation_limit",
        "working_memory_limit_bytes",
        "output_limit_bytes",
    ]
    .into_iter()
    .collect::<BTreeSet<_>>();
    let observed = projection
        .as_object()
        .expect("public v8 preflight object")
        .keys()
        .map(String::as_str)
        .collect::<BTreeSet<_>>();
    assert_eq!(observed, expected);
    for identity in projection["input_identities"]
        .as_array()
        .expect("sealed identity list")
    {
        let keys = identity
            .as_object()
            .expect("identity object")
            .keys()
            .map(String::as_str)
            .collect::<BTreeSet<_>>();
        assert_eq!(
            keys,
            [
                "role",
                "normalized_path",
                "observed_size_bytes",
                "observed_sha256",
                "declared_size_bytes",
                "declared_sha256",
            ]
            .into_iter()
            .collect::<BTreeSet<_>>()
        );
    }
    for record in projection["phases"].as_array().expect("sealed phase list") {
        let keys = record
            .as_object()
            .expect("phase object")
            .keys()
            .map(String::as_str)
            .collect::<BTreeSet<_>>();
        assert_eq!(
            keys,
            [
                "phase",
                "operation_units_cumulative",
                "retained_before_bytes",
                "transient_bytes",
                "retained_after_bytes",
                "phase_peak_bytes",
                "output_bytes_reserved",
            ]
            .into_iter()
            .collect::<BTreeSet<_>>()
        );
    }
}

fn write_profile(path: &Path) -> String {
    write_custom_profile(
        path,
        "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
    )
}

fn write_custom_profile(
    path: &Path,
    block_lengths: &str,
    haplotype_weights: &str,
    ld_profile: &str,
) -> String {
    write_custom_profile_with_gap(
        path,
        block_lengths,
        "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n",
        haplotype_weights,
        ld_profile,
    )
}

fn write_custom_profile_with_gap(
    path: &Path,
    block_lengths: &str,
    snp_gaps: &str,
    haplotype_weights: &str,
    ld_profile: &str,
) -> String {
    fs::create_dir_all(path).expect("profile directory");
    let files = [
        (
            "block_lengths.tsv",
            "f5_block_length_profile",
            block_lengths,
        ),
        ("snp_gaps.tsv", "f5_snp_gap_profile", snp_gaps),
        (
            "snp_spectrum.tsv",
            "f5_snp_spectrum",
            "alt_af_min\talt_af_max\tweight\n0\t1\t1\n",
        ),
        (
            "haplotype_weights.tsv",
            "f5_haplotype_weights",
            haplotype_weights,
        ),
        ("ld_profile.tsv", "f5_ld_profile", ld_profile),
    ];
    let links = files.iter().map(|(name, kind, contents)| {
        fs::write(path.join(name), contents).expect("profile file");
        format!(r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#, contents.len(), sha256(contents.as_bytes()))
    }).collect::<Vec<_>>().join(",");
    let metadata = format!(
        r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
    );
    fs::write(path.join("profile.metadata.json"), &metadata).expect("metadata");
    sha256(metadata.as_bytes())
}

fn configuration(root: &Path, mask: Option<&str>, limits: &str, bases: &str) -> PathBuf {
    let scenario_directory = root.join("authoring/scenarios");
    let profile_hash = write_profile(&scenario_directory.join("profiles"));
    write_configuration(root, mask, limits, bases, &profile_hash)
}

fn custom_configuration(
    root: &Path,
    limits: &str,
    bases: &str,
    block_lengths: &str,
    haplotype_weights: &str,
    ld_profile: &str,
) -> PathBuf {
    let scenario_directory = root.join("authoring/scenarios");
    let profile_hash = write_custom_profile(
        &scenario_directory.join("profiles"),
        block_lengths,
        haplotype_weights,
        ld_profile,
    );
    write_configuration(root, None, limits, bases, &profile_hash)
}

fn custom_configuration_with_gap(
    root: &Path,
    limits: &str,
    bases: &str,
    block_lengths: &str,
    snp_gaps: &str,
    haplotype_weights: &str,
    ld_profile: &str,
) -> PathBuf {
    let scenario_directory = root.join("authoring/scenarios");
    let profile_hash = write_custom_profile_with_gap(
        &scenario_directory.join("profiles"),
        block_lengths,
        snp_gaps,
        haplotype_weights,
        ld_profile,
    );
    write_configuration(root, None, limits, bases, &profile_hash)
}

fn write_configuration(
    root: &Path,
    mask: Option<&str>,
    limits: &str,
    bases: &str,
    profile_hash: &str,
) -> PathBuf {
    let scenario_directory = root.join("authoring/scenarios");
    let mask_fields = mask.map_or_else(String::new, |contents| {
        fs::write(scenario_directory.join("excluded_regions.tsv"), contents).expect("mask");
        format!(
            "excluded_regions_path = \"excluded_regions.tsv\"\nexcluded_regions_sha256 = \"{}\"\n",
            sha256(contents.as_bytes())
        )
    });
    let general = root.join("authoring/general.toml");
    fs::create_dir_all(root.join("authoring/references")).expect("reference directory");
    let fasta = format!(">chr1\n{bases}\n");
    fs::write(root.join("authoring/references/ref.fa"), &fasta).expect("reference");
    fs::write(&general, format!("config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n", "a".repeat(64), sha256(fasta.as_bytes()))).expect("general");
    fs::write(scenario_directory.join("scenario.toml"), format!("config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = 1\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{profile_hash}\"\n{mask_fields}max_layout_attempts = 1\nmax_haplotype_attempts = 1\nmax_swap_proposals = 1\nmax_haplotype_cells = 16\n{limits}block_count = 1\n\n[[subpopulations]]\npopulation_id = \"P1\"\nindividual_count = 1\nfounder_ids = [\"I1\"]\n")).expect("scenario");
    general
}

fn replace_scenario(root: &Path, from: &str, to: &str) {
    let path = root.join("authoring/scenarios/scenario.toml");
    let source = fs::read_to_string(&path).expect("scenario source");
    assert!(source.contains(from), "missing scenario fragment: {from}");
    fs::write(path, source.replace(from, to)).expect("scenario rewrite");
}

fn set_scenario_operation_limit(root: &Path, operation_limit: u64) {
    let path = root.join("authoring/scenarios/scenario.toml");
    let source = fs::read_to_string(&path).expect("scenario source");
    let mut replaced = false;
    let rewritten = source
        .lines()
        .map(|line| {
            if line.starts_with("max_operation_units = ") {
                replaced = true;
                format!("max_operation_units = {operation_limit}")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert!(replaced, "scenario operation limit");
    fs::write(path, format!("{rewritten}\n")).expect("scenario operation-limit rewrite");
}

fn literal_fixture_operation_oracle(root: &Path, mask: Option<&[u8]>) -> u64 {
    let authoring = root.join("authoring");
    let scenario = authoring.join("scenarios");
    let profile = scenario.join("profiles");
    let fasta = fs::read(authoring.join("references/ref.fa")).expect("literal FASTA bytes");
    assert_eq!(fasta, b">chr1\nAAA\n");
    match mask {
        Some(expected) => assert_eq!(
            fs::read(scenario.join("excluded_regions.tsv")).expect("literal mask bytes"),
            expected
        ),
        None => assert!(!scenario.join("excluded_regions.tsv").exists()),
    }

    let general_bytes = fs::read(authoring.join("general.toml"))
        .expect("literal general TOML bytes")
        .len() as u64;
    let scenario_bytes = fs::read(scenario.join("scenario.toml"))
        .expect("literal scenario TOML bytes")
        .len() as u64;
    let metadata_bytes = fs::read(profile.join("profile.metadata.json"))
        .expect("literal profile metadata bytes")
        .len() as u64;
    let profile_table_bytes = [
        "block_lengths.tsv",
        "snp_gaps.tsv",
        "snp_spectrum.tsv",
        "haplotype_weights.tsv",
        "ld_profile.tsv",
    ]
    .into_iter()
    .map(|name| {
        fs::read(profile.join(name))
            .expect("literal profile TSV bytes")
            .len() as u64
    })
    .sum::<u64>();
    let mask_bytes = mask.map_or(0, |bytes| bytes.len() as u64);

    // Oráculo deliberadamente desplegado para este fixture literal. No usa
    // `F5ResourcePreflight`, helpers aritméticos ni la fórmula productiva.
    [
        0,                   // other_operation_units
        general_bytes,       // lectura/hash del TOML general
        general_bytes,       // parser del TOML general
        scenario_bytes,      // lectura/hash del escenario
        scenario_bytes,      // parser del escenario
        metadata_bytes,      // lectura/hash de metadata
        metadata_bytes,      // parser de metadata
        profile_table_bytes, // lectura/hash de los cinco TSV
        profile_table_bytes, // parser de los cinco TSV
        fasta.len() as u64,  // primer recorrido FASTA
        fasta.len() as u64,  // recuperación REF posterior
        3,                   // tres bases de referencia literales
        mask_bytes,          // único recorrido físico de la máscara
        1,                   // ID FASTA retenido para unicidad
        1,                   // tramo emitido
        1,                   // visita de capacidad del tramo
        1,                   // visita de validación del tramo
        1,                   // primera consulta directa de layout
        1,                   // segunda consulta directa de layout
        1,                   // copia del tramo en el intento de layout
        3,                   // desplazamientos remove/insert del layout
        6,                   // cuatro TSV con una fila y dos haplotipos
        5,                   // 1 + 1 + 1 + 2 filas normalizadas
        1,                   // pasada Neumaier del soporte de longitudes
        1,                   // pasada de prefijo del soporte de longitudes
        1,
        1,
        1,
        1, // dos sorteos de gap, dos pasadas literales cada uno
        1,
        1,
        1,  // B + S visitas de request para B=1, S=2
        8,  // evaluación inicial: dos haplotipos, cuatro pares
        2,  // propuesta: un swap por dos sitios
        2,  // dos copias fundadoras del bloque
        4,  // dos copias por una secuencia y dos sitios
        2,  // dos genotipos SNP máximos
        2,  // dos resúmenes de sitio
        11, // filas TSV del fixture
        30, // elementos JSON repetidos del fixture
    ]
    .into_iter()
    .try_fold(0_u64, u64::checked_add)
    .expect("literal fixture operation oracle fits u64")
}

fn remove_declared_seed(root: &Path) {
    let path = root.join("authoring/general.toml");
    let source = fs::read_to_string(&path).expect("general source");
    fs::write(
        path,
        source.replace(&format!("seed = \"{}\"\n", "a".repeat(64)), ""),
    )
    .expect("seedless general");
}

fn replace_mask(root: &Path, contents: &[u8]) {
    fs::write(
        root.join("authoring/scenarios/excluded_regions.tsv"),
        contents,
    )
    .expect("mask rewrite");
    let scenario = root.join("authoring/scenarios/scenario.toml");
    let source = fs::read_to_string(&scenario).expect("scenario source");
    let rewritten = source
        .lines()
        .map(|line| {
            if line.starts_with("excluded_regions_sha256 = ") {
                format!("excluded_regions_sha256 = \"{}\"", sha256(contents))
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(scenario, format!("{rewritten}\n")).expect("scenario mask digest rewrite");
}

fn rename_scenario(root: &Path, current: &str, renamed: &str) {
    for path in [
        root.join("authoring/general.toml"),
        root.join("authoring/scenarios/scenario.toml"),
    ] {
        let source = fs::read_to_string(&path).expect("scenario-bearing source");
        fs::write(
            path,
            source.replace(
                &format!("scenario_id = \"{current}\""),
                &format!("scenario_id = \"{renamed}\""),
            ),
        )
        .expect("scenario rename");
    }
}

fn replace_reference(root: &Path, fasta: &[u8]) {
    fs::write(root.join("authoring/references/ref.fa"), fasta).expect("reference rewrite");
    let path = root.join("authoring/general.toml");
    let source = fs::read_to_string(&path).expect("general source");
    let rewritten = source
        .lines()
        .map(|line| {
            if line.starts_with("reference_sha256 = ") {
                format!("reference_sha256 = \"{}\"", sha256(fasta))
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(path, format!("{rewritten}\n")).expect("general rewrite");
}

fn replace_profile_source_citation(root: &Path, citation: &str) {
    let metadata_path = root.join("authoring/scenarios/profiles/profile.metadata.json");
    let metadata = fs::read_to_string(&metadata_path).expect("profile metadata");
    let replacement = format!(r#""source_citation":"{citation}""#);
    let metadata = metadata.replace(r#""source_citation":"test""#, &replacement);
    fs::write(&metadata_path, &metadata).expect("rewritten profile metadata");

    let scenario_path = root.join("authoring/scenarios/scenario.toml");
    let scenario = fs::read_to_string(&scenario_path).expect("scenario source");
    let scenario = scenario
        .lines()
        .map(|line| {
            if line.starts_with("profile_metadata_sha256 = ") {
                format!(
                    "profile_metadata_sha256 = \"{}\"",
                    sha256(metadata.as_bytes())
                )
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(scenario_path, format!("{scenario}\n")).expect("rewritten scenario");
}

fn replace_reference_with_repeated_base(root: &Path, base_count: usize) {
    let fasta_path = root.join("authoring/references/ref.fa");
    let mut fasta = fs::File::create(&fasta_path).expect("large reference");
    let mut digest = Sha256::new();
    let header = b">chr1\n";
    fasta.write_all(header).expect("large reference header");
    digest.update(header);
    let chunk = [b'A'; 65_536];
    let mut remaining = base_count;
    while remaining != 0 {
        let count = remaining.min(chunk.len());
        fasta
            .write_all(&chunk[..count])
            .expect("large reference bases");
        digest.update(&chunk[..count]);
        remaining -= count;
    }
    fasta.write_all(b"\n").expect("large reference terminator");
    digest.update(b"\n");
    drop(fasta);

    let path = root.join("authoring/general.toml");
    let source = fs::read_to_string(&path).expect("general source");
    let observed = format!("{:x}", digest.finalize());
    let rewritten = source
        .lines()
        .map(|line| {
            if line.starts_with("reference_sha256 = ") {
                format!("reference_sha256 = \"{observed}\"")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(path, format!("{rewritten}\n")).expect("general rewrite");
}

#[test]
fn seedless_entropy_is_resolved_only_after_complete_physical_admission() {
    let success_root = temporary("entropy-success");
    fs::create_dir_all(&success_root).expect("success root");
    let success = configuration(&success_root, None, "", "AAAA");
    remove_declared_seed(&success_root);
    let mut entropy = CountingEntropy::default();
    let admitted = load_and_validate_f5_config_with_entropy(&success, &mut entropy)
        .expect("successful seedless admission");
    assert_eq!(entropy.calls, 1);
    assert_eq!(admitted.seed(), "bb".repeat(32));

    // Una referencia vacía o truncada la rechaza el escáner de este hito, no el
    // lector de F1: son lectores distintos y la cobertura de uno no acredita al
    // otro.
    for (label, prepare) in [
        ("fasta-missing", 0_u8),
        ("fasta-invalid", 1),
        ("fasta-sha", 2),
        ("fasta-empty", 3),
        ("fasta-header-only", 4),
        ("fasta-truncated-header", 5),
    ] {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("failure root");
        let general = configuration(&root, None, "", "AAAA");
        remove_declared_seed(&root);
        let fasta = root.join("authoring/references/ref.fa");
        match prepare {
            0 => fs::remove_file(fasta).expect("remove FASTA"),
            1 => replace_reference(&root, b">chr1\nA?AA\n"),
            2 => fs::write(fasta, b">chr1\nAAAT\n").expect("mutate FASTA without digest"),
            3 => replace_reference(&root, b""),
            4 => replace_reference(&root, b">chr1\n"),
            5 => replace_reference(&root, b">ch"),
            _ => unreachable!(),
        }
        let mut entropy = CountingEntropy::default();
        load_and_validate_f5_config_with_entropy(&general, &mut entropy)
            .expect_err("invalid FASTA admission");
        assert_eq!(entropy.calls, 0, "{label}");
        fs::remove_dir_all(root).expect("failure cleanup");
    }

    let valid_mask = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\n";
    for (label, prepare) in [
        ("mask-missing", 0_u8),
        ("mask-invalid", 1),
        ("mask-sha", 2),
        ("mask-empty", 3),
        ("mask-truncated-row", 4),
    ] {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("failure root");
        let general = configuration(&root, Some(valid_mask), "", "AAAA");
        remove_declared_seed(&root);
        let mask = root.join("authoring/scenarios/excluded_regions.tsv");
        match prepare {
            0 => fs::remove_file(mask).expect("remove mask"),
            1 => replace_mask(&root, b"bad-header\n"),
            2 => fs::write(mask, valid_mask.replace("2\t2", "3\t3"))
                .expect("mutate mask without digest"),
            3 => replace_mask(&root, b""),
            4 => replace_mask(
                &root,
                b"region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2",
            ),
            _ => unreachable!(),
        }
        let mut entropy = CountingEntropy::default();
        let error = load_and_validate_f5_config_with_entropy(&general, &mut entropy)
            .expect_err("invalid mask admission");
        assert_eq!(
            error.code(),
            match prepare {
                0 => "GPS-INPUT-001",
                // Una máscara vacía o con la última fila sin LF es un fallo de
                // formato, no de identidad: se clasifica antes de comparar el
                // resumen declarado.
                1 | 3 | 4 => "GPS-FMT-101",
                2 => "GPS-INPUT-002",
                _ => unreachable!(),
            },
            "{label}"
        );
        assert_eq!(entropy.calls, 0, "{label}");
        fs::remove_dir_all(root).expect("failure cleanup");
    }

    for (label, mask) in [("fasta-budget", None), ("mask-budget", Some(valid_mask))] {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("budget root");
        let general = configuration(&root, mask, "max_operation_units = 9999\n", "AAAA");
        remove_declared_seed(&root);
        let mut baseline_entropy = CountingEntropy::default();
        let baseline = load_and_validate_f5_config_with_entropy(&general, &mut baseline_entropy)
            .expect("baseline admission");
        let required = preflight_f5_resources(&baseline, F5OtherResources::default())
            .expect("baseline preflight")
            .operation_units();
        replace_scenario(
            &root,
            "max_operation_units = 9999\n",
            &format!("max_operation_units = {}\n", required - 1),
        );
        let mut entropy = CountingEntropy::default();
        load_and_validate_f5_config_with_entropy(&general, &mut entropy)
            .expect_err("operation limit minus one");
        assert_eq!(entropy.calls, 0, "{label}");
        fs::remove_dir_all(root).expect("budget cleanup");
    }
    fs::remove_dir_all(success_root).expect("success cleanup");
}

#[test]
fn known_external_exhaustion_rejects_before_reference_or_mask_and_entropy() {
    let root = temporary("known-other-before-input");
    fs::create_dir_all(&root).expect("root");
    let mask = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\n";
    let general = configuration(&root, Some(mask), "", "AAAA");
    remove_declared_seed(&root);
    fs::remove_file(root.join("authoring/references/ref.fa")).expect("remove reference sentinel");
    fs::remove_file(root.join("authoring/scenarios/excluded_regions.tsv"))
        .expect("remove mask sentinel");
    let mut entropy = CountingEntropy::default();
    let error = load_and_validate_f5_config_with_entropy_and_resources(
        &general,
        &mut entropy,
        F5OtherResources::new(F5_HARD_OPERATION_LIMIT + 1, 0, 0),
    )
    .expect_err("known external work exceeds the ledger");
    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(error.context()["dimension"], "operation_units");
    assert_eq!(entropy.calls, 0);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn public_preflight_reuses_the_initial_other_binding_without_late_io() {
    let root = temporary("public-other-binding");
    fs::create_dir_all(&root).expect("root");
    let mask = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\n";
    let general = configuration(&root, Some(mask), "", "AAAA");
    remove_declared_seed(&root);
    let other = F5OtherResources::new(7, 11, 13);
    let mut entropy = CountingEntropy::default();
    let config =
        load_and_validate_f5_config_with_entropy_and_resources(&general, &mut entropy, other)
            .expect("initial staged admission");
    assert_eq!(entropy.calls, 1);
    let first = preflight_f5_resources(&config, other).expect("matching public recovery");

    fs::remove_file(root.join("authoring/references/ref.fa")).expect("reference sentinel");
    fs::remove_file(root.join("authoring/scenarios/excluded_regions.tsv")).expect("mask sentinel");
    let second = preflight_f5_resources(&config, other)
        .expect("matching recovery must neither reopen nor recalculate");
    assert_eq!(second, first);
    assert_eq!(
        second.eligible_tracts().as_ptr(),
        first.eligible_tracts().as_ptr()
    );
    assert_eq!(entropy.calls, 1);

    let error = preflight_f5_resources(&config, F5OtherResources::new(8, 11, 13))
        .expect_err("a late contribution must not authorize another admission");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(error.context()["fields"], "preflight.other_resources");
    assert_eq!(entropy.calls, 1);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn sealed_receipt_binds_every_input_in_contractual_role_order_with_mask() {
    let root = temporary("scan-identities");
    fs::create_dir_all(&root).expect("root");
    let mask = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\n";
    let general = configuration(&root, Some(mask), "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("admitted configuration");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("sealed admission");
    let identities = sealed.receipt().input_identities();
    // Oráculo literal del orden contractual, independiente del enum productivo.
    assert_eq!(
        identities
            .iter()
            .map(|identity| identity.role().as_str())
            .collect::<Vec<_>>(),
        vec![
            "general_toml",
            "scenario_toml",
            "profile_metadata",
            "profile_block_lengths",
            "profile_snp_gaps",
            "profile_snp_spectrum",
            "profile_haplotype_weights",
            "profile_ld",
            "reference_fasta",
            "mask",
        ]
    );
    let fasta = b">chr1\nAAAA\n";
    let reference = &identities[8];
    assert_eq!(reference.observed_size_bytes(), fasta.len() as u64);
    assert_eq!(reference.observed_sha256(), sha256(fasta));
    assert_eq!(reference.declared_size_bytes(), None);
    assert_eq!(reference.declared_sha256(), Some(sha256(fasta).as_str()));
    assert_eq!(
        reference.normalized_path(),
        root.join("authoring/references/ref.fa")
    );
    let sealed_mask = &identities[9];
    assert_eq!(sealed_mask.observed_size_bytes(), mask.len() as u64);
    assert_eq!(sealed_mask.observed_sha256(), sha256(mask.as_bytes()));
    assert_eq!(sealed_mask.declared_size_bytes(), None);
    assert_eq!(
        sealed_mask.declared_sha256(),
        Some(sha256(mask.as_bytes()).as_str())
    );
    assert_eq!(
        sealed_mask.normalized_path(),
        root.join("authoring/scenarios/excluded_regions.tsv")
    );
    // Los dos TOML de autoría no tienen documento propietario que declare su
    // identidad; los cinco TSV sí la tienen completa en el metadata.
    for identity in &identities[..2] {
        assert_eq!(identity.declared_size_bytes(), None);
        assert_eq!(identity.declared_sha256(), None);
    }
    assert_eq!(identities[2].declared_size_bytes(), None);
    assert!(identities[2].declared_sha256().is_some());
    for identity in &identities[3..8] {
        assert_eq!(
            identity.declared_size_bytes(),
            Some(identity.observed_size_bytes())
        );
        assert_eq!(identity.declared_sha256(), Some(identity.observed_sha256()));
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn sealed_receipt_without_a_declared_mask_stops_at_the_reference() {
    let root = temporary("scan-identities-no-mask");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("admitted configuration");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("sealed admission");
    let identities = sealed.receipt().input_identities();
    assert_eq!(identities.len(), 9);
    assert_eq!(identities[8].role().as_str(), "reference_fasta");
    assert!(
        identities
            .iter()
            .all(|identity| identity.role().as_str() != "mask")
    );
    fs::remove_dir_all(root).expect("cleanup");
}

fn replace_reference_with_long_description(root: &Path, description_bytes: usize) {
    let fasta_path = root.join("authoring/references/ref.fa");
    let mut fasta = fs::File::create(&fasta_path).expect("long-header reference");
    let mut digest = Sha256::new();
    let prefix = b">chr1 ";
    fasta.write_all(prefix).expect("long-header prefix");
    digest.update(prefix);
    let chunk = [b'd'; 65_536];
    let mut remaining = description_bytes;
    while remaining != 0 {
        let count = remaining.min(chunk.len());
        fasta
            .write_all(&chunk[..count])
            .expect("long-header description");
        digest.update(&chunk[..count]);
        remaining -= count;
    }
    let suffix = b"\nAA\n";
    fasta.write_all(suffix).expect("long-header sequence");
    digest.update(suffix);
    drop(fasta);
    replace_reference_digest(root, format!("{:x}", digest.finalize()));
}

fn replace_reference_with_alternating_bases(root: &Path, base_count: usize) {
    let fasta_path = root.join("authoring/references/ref.fa");
    let mut fasta = fs::File::create(&fasta_path).expect("alternating reference");
    let mut digest = Sha256::new();
    let header = b">chr1\n";
    fasta.write_all(header).expect("alternating header");
    digest.update(header);
    let chunk = b"AN".repeat(32_768);
    let mut remaining = base_count;
    while remaining != 0 {
        let count = remaining.min(chunk.len());
        fasta.write_all(&chunk[..count]).expect("alternating bases");
        digest.update(&chunk[..count]);
        remaining -= count;
    }
    fasta.write_all(b"\n").expect("alternating terminator");
    digest.update(b"\n");
    drop(fasta);
    replace_reference_digest(root, format!("{:x}", digest.finalize()));
}

fn replace_reference_with_many_maximum_ids(root: &Path, sequence_count: usize, base: char) {
    assert!(matches!(base, 'A' | 'N'));
    let fasta_path = root.join("authoring/references/ref.fa");
    let mut fasta = fs::File::create(&fasta_path).expect("many-sequence reference");
    let mut digest = Sha256::new();
    for ordinal in 0..sequence_count {
        let id = format!("s{ordinal:015x}{}", "x".repeat(48));
        assert_eq!(id.len(), 64);
        let record = format!(">{id}\n{base}\n");
        fasta.write_all(record.as_bytes()).expect("sequence record");
        digest.update(record.as_bytes());
    }
    drop(fasta);
    replace_reference_digest(root, format!("{:x}", digest.finalize()));
}

fn replace_reference_with_n_only_ids_and_eligible_tail(root: &Path, n_only_count: usize) {
    let fasta_path = root.join("authoring/references/ref.fa");
    let mut fasta = fs::File::create(&fasta_path).expect("many-sequence reference");
    let mut digest = Sha256::new();
    for ordinal in 0..n_only_count {
        let record = format!(">n{ordinal:015x}\nN\n");
        fasta
            .write_all(record.as_bytes())
            .expect("N-only sequence record");
        digest.update(record.as_bytes());
    }
    let tail = b">eligible\nAA\n";
    fasta.write_all(tail).expect("eligible sequence record");
    digest.update(tail);
    drop(fasta);
    replace_reference_digest(root, format!("{:x}", digest.finalize()));
}

fn replace_reference_digest(root: &Path, observed: String) {
    let path = root.join("authoring/general.toml");
    let source = fs::read_to_string(&path).expect("general source");
    let rewritten = source
        .lines()
        .map(|line| {
            if line.starts_with("reference_sha256 = ") {
                format!("reference_sha256 = \"{observed}\"")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(path, format!("{rewritten}\n")).expect("general rewrite");
}

struct FinalFormulaMemoryInputs {
    block_count: u64,
    tract_count: u64,
    retained_sequence_ids: u64,
    profile_rows: u64,
    normalization_rows: u64,
    population_count: u64,
    haplotype_count: u64,
    site_bound: u64,
    bin_count_total: u64,
    founder_copy_count: u64,
    /// `sequence_count`: el segundo eje de la identidad de copia. El
    /// repositorio conserva una copia por individuo **y** secuencia, de modo
    /// que `X = C_total * sequence_count` y no `C_total`.
    sequence_count: u64,
}

fn final_formula_memory_oracle(inputs: FinalFormulaMemoryInputs) -> u64 {
    let pair_bound = inputs.site_bound * (inputs.site_bound - 1) / 2;
    let hs_bytes = (inputs.haplotype_count * inputs.site_bound).div_ceil(8);
    let copy_sequences = inputs.founder_copy_count * inputs.sequence_count;
    let seats = inputs.block_count * inputs.founder_copy_count;
    let registry_entries = copy_sequences
        + inputs.block_count * inputs.population_count * inputs.haplotype_count
        + seats;
    256 * (4 * inputs.block_count
        + 2 * inputs.tract_count
        + inputs.retained_sequence_ids
        + inputs.profile_rows
        + inputs.normalization_rows)
        + inputs.block_count
            * inputs.population_count
            * (hs_bytes
                + (160 + inputs.site_bound) * inputs.haplotype_count
                + 256 * inputs.site_bound
                + 192 * pair_bound)
        + 256 * inputs.block_count * inputs.bin_count_total
        + 192 * copy_sequences
        + 64 * seats
        + 128 * registry_entries
        + 65_536
}

#[test]
fn fragmented_mask_changes_only_tract_terms_and_is_bounded_inclusively() {
    let root = temporary("fragmented");
    fs::create_dir_all(&root).expect("root");
    let length = 4_096 * 2 + 4_095;
    let bases = "A".repeat(length);
    let compact_mask = format!(
        "region_id\tsequence_id\tstart_1based\tend_1based\ncompact\tchr1\t8193\t{length}\n"
    );
    let compact_path = configuration(&root.join("compact"), Some(&compact_mask), "", &bases);
    let compact_config = load_and_validate_f5_config(&compact_path).expect("compact configuration");
    let compact = preflight_f5_resources(&compact_config, F5OtherResources::default())
        .expect("compact preflight");
    let mask = format!(
        "region_id\tsequence_id\tstart_1based\tend_1based\n{}",
        (0..4_095)
            .map(|index| format!("m{index}\tchr1\t{}\t{}\n", index * 3 + 3, index * 3 + 3))
            .collect::<String>()
    );
    let fragmented_path = configuration(&root.join("fragmented"), Some(&mask), "", &bases);
    let fragmented_config =
        load_and_validate_f5_config(&fragmented_path).expect("fragmented configuration");
    let fragmented = preflight_f5_resources(&fragmented_config, F5OtherResources::default())
        .expect("fragmented preflight");
    assert_eq!(compact.tract_count(), 1);
    assert_eq!(compact.retained_sequence_ids(), 1);
    assert_eq!(fragmented.tract_count(), 4_096);
    let recovered = preflight_f5_resources(&fragmented_config, F5OtherResources::default())
        .expect("public recovery of the admitted fragmented geometry");
    assert_eq!(
        fragmented.eligible_tracts().as_ptr(),
        recovered.eligible_tracts().as_ptr(),
        "public recovery must share the single admitted tract allocation"
    );
    assert_eq!(fragmented.reference_bases_processed(), length as u64);
    assert_eq!(
        compact.fasta_bytes_processed(),
        fs::metadata(root.join("compact/authoring/references/ref.fa"))
            .expect("compact FASTA metadata")
            .len()
    );
    assert_eq!(
        fragmented.fasta_bytes_processed(),
        fs::metadata(root.join("fragmented/authoring/references/ref.fa"))
            .expect("fragmented FASTA metadata")
            .len()
    );
    assert_eq!(compact.eligible_bases(), fragmented.eligible_bases());
    assert_eq!(compact.profile_rows(), 6);
    assert_eq!(compact.normalization_rows(), 5);
    assert_eq!(compact.bin_count_total(), 1);
    assert_eq!(compact.tsv_rows_bound(), 11);
    assert_eq!(compact.json_rows_bound(), 30);
    let operation_without_gap_selector = 41_580;
    assert_eq!(compact.gap_selector_operation_units(), 4);
    assert_eq!(
        compact.operation_units(),
        operation_without_gap_selector + compact.gap_selector_operation_units()
    );
    let compact_final_formula = final_formula_memory_oracle(FinalFormulaMemoryInputs {
        block_count: 1,
        tract_count: 1,
        retained_sequence_ids: 1,
        profile_rows: 6,
        normalization_rows: 5,
        population_count: 1,
        haplotype_count: 2,
        site_bound: 2,
        bin_count_total: 1,
        founder_copy_count: 2,
        sequence_count: 1,
    });
    let independent_scan_peak = 65_536 + 65_536 + 256 * (6 + 5 + 1 + 1);
    assert_eq!(
        compact.working_memory_bytes(),
        compact_final_formula.max(independent_scan_peak)
    );
    assert_eq!(compact.output_bytes(), 1_306_624);
    assert_eq!(
        fragmented.operation_units() - compact.operation_units(),
        4_095 * (1 + 2 + 2 + 1 + 3) + fragmented.mask_bytes_processed()
            - compact.mask_bytes_processed()
    );
    let fragmented_final_formula = final_formula_memory_oracle(FinalFormulaMemoryInputs {
        block_count: 1,
        tract_count: 4_096,
        retained_sequence_ids: 1,
        profile_rows: 6,
        normalization_rows: 5,
        population_count: 1,
        haplotype_count: 2,
        site_bound: 2,
        bin_count_total: 1,
        founder_copy_count: 2,
        sequence_count: 1,
    });
    assert_eq!(fragmented.working_memory_bytes(), fragmented_final_formula);
    assert_eq!(
        fragmented.working_memory_bytes() - compact.working_memory_bytes(),
        fragmented_final_formula - independent_scan_peak
    );
    assert_eq!(fragmented.output_bytes(), compact.output_bytes());

    for (label, field, placeholder) in [
        ("operations", "max_operation_units", 999_999_u64),
        ("memory", "max_working_memory_bytes", 9_999_999_u64),
    ] {
        let placeholder_line = format!("{field} = {placeholder}\n");
        let path = configuration(&root.join(label), Some(&mask), &placeholder_line, &bases);
        let provisional = load_and_validate_f5_config(&path).expect("provisional configuration");
        let provisional = preflight_f5_resources(&provisional, F5OtherResources::default())
            .expect("provisional bound");
        let exact = if field == "max_operation_units" {
            provisional.operation_units()
        } else {
            provisional.working_memory_bytes()
        };
        assert_eq!(exact.to_string().len(), placeholder.to_string().len());
        replace_scenario(
            &root.join(label),
            &placeholder_line,
            &format!("{field} = {exact}\n"),
        );
        let exact_config = load_and_validate_f5_config(&path).expect("exact configuration");
        preflight_f5_resources(&exact_config, F5OtherResources::default())
            .expect("inclusive fragmented bound");
        replace_scenario(
            &root.join(label),
            &format!("{field} = {exact}\n"),
            &format!("{field} = {}\n", exact - 1),
        );
        assert_eq!(
            load_and_validate_f5_config(&path)
                .expect_err("bound minus one must reject")
                .code(),
            "GPS-RES-101"
        );
    }

    let exact_path = configuration(
        &root.join("operations-exact-layout"),
        Some(&mask),
        "max_operation_units = 999999\n",
        &bases,
    );
    let provisional = load_and_validate_f5_config(&exact_path).expect("layout provisional config");
    let exact_operation = preflight_f5_resources(&provisional, F5OtherResources::default())
        .expect("layout provisional bound")
        .operation_units();
    replace_scenario(
        &root.join("operations-exact-layout"),
        "max_operation_units = 999999\n",
        &format!("max_operation_units = {exact_operation}\n"),
    );
    let exact_config = load_and_validate_f5_config(&exact_path).expect("exact layout config");
    let exact_preflight = preflight_f5_resources(&exact_config, F5OtherResources::default())
        .expect("the exact three-pass operation budget is admitted");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    build_f5_physical_layout(
        F5LayoutRequest::new(&exact_config, &exact_preflight, &execution)
            .expect("exact preflight binding"),
    )
    .expect("the admitted layout must complete its counted REF pass");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn the_sealed_reference_identity_is_reused_without_reopening_the_fasta() {
    let root = temporary("reference-sha");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("configuration");
    fs::write(root.join("authoring/references/ref.fa"), b">chr1\nAAAT\n")
        .expect("reference mutation");
    let preflight = preflight_f5_resources(&config, F5OtherResources::default())
        .expect("recovering the sealed admission never reopens the reference");
    assert_eq!(
        preflight.receipt().input_identities()[8].observed_sha256(),
        sha256(b">chr1\nAAAA\n"),
        "the receipt keeps the bytes the single physical pass observed"
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn exact_fraction_avoids_intermediate_overflow_and_floors_once() {
    let root = temporary("exact-fraction");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", &"A".repeat(100));
    replace_scenario(
        &root,
        "block_count = 1",
        "coverage_fraction = 0.500000000000000000",
    );
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let result =
        preflight_f5_resources(&config, F5OtherResources::default()).expect("exact fraction");
    assert_eq!(result.target_bases(), Some(50));
    assert_eq!(result.block_count(), 25);

    replace_scenario(
        &root,
        "coverage_fraction = 0.500000000000000000",
        "coverage_fraction = 0.999999999999999999",
    );
    replace_reference(&root, format!(">chr1\n{}\n", "A".repeat(101)).as_bytes());
    let config = load_and_validate_f5_config(&general).expect("extreme configuration");
    let result = preflight_f5_resources(&config, F5OtherResources::default())
        .expect("extreme exact fraction");
    assert_eq!(result.target_bases(), Some(100));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn impossible_targets_and_override_outside_resolved_blocks_are_rejected() {
    for (label, target) in [
        ("target-over-eligible", "covered_bases = 5"),
        ("discontinuous-support", "covered_bases = 3"),
    ] {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("root");
        let general = configuration(&root, None, "", "AAAA");
        replace_scenario(&root, "block_count = 1", target);
        assert_eq!(
            load_and_validate_f5_config(&general)
                .expect_err("impossible target is rejected during admission")
                .code(),
            "GPS-DOMAIN-101"
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    let root = temporary("override-block");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    replace_scenario(
        &root,
        "founder_ids = [\"I1\"]",
        "founder_ids = [\"I1\"]\n\n[[sampling_overrides]]\npopulation_id = \"P1\"\nblock_index = 1\nmode = \"probabilistic\"",
    );
    assert_eq!(
        load_and_validate_f5_config(&general)
            .expect_err("override outside B is rejected during admission")
            .code(),
        "GPS-CONF-101"
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn covered_bases_reserves_its_known_formula_before_opening_reference() {
    let root = temporary("covered-known-before-reference");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration(
        &root,
        "max_operation_units = 999999\n",
        "AAAAAAAAAA",
        "length_min_bp\tlength_max_bp\tweight\n2\t4\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t3\t0.5\t0.05\n",
    );
    replace_scenario(&root, "block_count = 1", "covered_bases = 10");
    let provisional = load_and_validate_f5_config(&general).expect("covered-bases admission");
    let preflight = preflight_f5_resources(&provisional, F5OtherResources::default())
        .expect("covered-bases preflight");
    assert_eq!(preflight.block_count(), 3);

    // Oráculo independiente de los términos que solo se aprenden al
    // recorrer esta referencia de un tramo y una secuencia. El resto debe
    // estar ya cargado al ledger antes de File::open.
    let block_count = preflight.block_count();
    let tract_count = preflight.tract_count();
    let layout_attempts = preflight.max_layout_attempts();
    let site_bound = preflight.site_bound_per_block();
    let scan_learned = 2 * preflight.fasta_bytes_processed()
        + preflight.reference_bases_processed()
        + preflight.mask_bytes_processed()
        + preflight.retained_sequence_ids()
        + tract_count
        + 2 * tract_count
        + 2 * layout_attempts * block_count * tract_count
        + layout_attempts * (tract_count + 3 * block_count * tract_count)
        + 2 * preflight.retained_sequence_ids() * block_count * site_bound;
    let known_at_original_width = preflight.operation_units() - scan_learned;

    let scenario_path = root.join("authoring/scenarios/scenario.toml");
    let original = fs::read_to_string(&scenario_path).expect("scenario source");
    let original_len = original.len() as i128;
    let mut exact = known_at_original_width;
    loop {
        let rewritten = original.replace(
            "max_operation_units = 999999",
            &format!("max_operation_units = {exact}"),
        );
        let adjusted = (i128::from(known_at_original_width)
            + 2 * (rewritten.len() as i128 - original_len)) as u64;
        if adjusted == exact {
            fs::write(&scenario_path, rewritten).expect("exact known-cost scenario");
            break;
        }
        exact = adjusted;
    }
    fs::remove_file(root.join("authoring/references/ref.fa")).expect("reference sentinel");
    let exact_error = load_and_validate_f5_config(&general)
        .expect_err("the exact known bound reaches the reference open");
    assert_eq!(exact_error.code(), "GPS-INPUT-001");

    replace_scenario(
        &root,
        &format!("max_operation_units = {exact}"),
        &format!("max_operation_units = {}", exact - 1),
    );
    let minus_one = load_and_validate_f5_config(&general)
        .expect_err("one less than the known formula rejects before reference");
    assert_eq!(minus_one.code(), "GPS-RES-101");
    assert_eq!(minus_one.context()["dimension"], "operation_units");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn invalid_covered_bases_and_site_bounds_win_before_absent_reference() {
    let invalid_target_root = temporary("invalid-covered-before-reference");
    fs::create_dir_all(&invalid_target_root).expect("target root");
    let general = custom_configuration(
        &invalid_target_root,
        "",
        "AAAAAA",
        "length_min_bp\tlength_max_bp\tweight\n3\t4\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t3\t0.5\t0.05\n",
    );
    replace_scenario(&invalid_target_root, "block_count = 1", "covered_bases = 5");
    fs::remove_file(invalid_target_root.join("authoring/references/ref.fa"))
        .expect("target reference sentinel");
    let error = load_and_validate_f5_config(&general)
        .expect_err("invalid covered target must precede reference open");
    assert_eq!(error.code(), "GPS-DOMAIN-101");

    let site_root = temporary("invalid-site-bound-before-reference");
    fs::create_dir_all(&site_root).expect("site root");
    let general = custom_configuration(
        &site_root,
        "",
        &"A".repeat(4_097),
        "length_min_bp\tlength_max_bp\tweight\n4097\t4097\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t4096\t0.5\t0.05\n",
    );
    replace_scenario(
        &site_root,
        "max_haplotype_cells = 16",
        "max_haplotype_cells = 8192",
    );
    fs::remove_file(site_root.join("authoring/references/ref.fa"))
        .expect("site reference sentinel");
    let error =
        load_and_validate_f5_config(&general).expect_err("S_bound must precede reference open");
    assert_eq!(error.code(), "GPS-RES-101");
    fs::remove_dir_all(invalid_target_root).expect("target cleanup");
    fs::remove_dir_all(site_root).expect("site cleanup");
}

#[test]
fn mask_scanner_is_ordered_across_sequences_and_has_a_bounded_row() {
    let root = temporary("multi-sequence-mask");
    fs::create_dir_all(&root).expect("root");
    let mask = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\nm2\tchr2\t3\t3\n";
    let general = configuration(&root, Some(mask), "", "AAAA");
    replace_reference(&root, b">chr1\nAAAA\n>chr2\nAAAA\n");
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let result =
        preflight_f5_resources(&config, F5OtherResources::default()).expect("multi-sequence mask");
    assert_eq!(result.reference_bases_processed(), 8);
    assert_eq!(result.eligible_bases(), 6);
    assert_eq!(result.tract_count(), 4);
    fs::remove_dir_all(root).expect("cleanup");

    let root = temporary("bounded-mask-row");
    fs::create_dir_all(&root).expect("root");
    let oversized = format!(
        "region_id\tsequence_id\tstart_1based\tend_1based\n{}\tchr1\t1\t1\n",
        "r".repeat(65_536)
    );
    let general = configuration(&root, Some(&oversized), "", "AAAA");
    assert_eq!(
        load_and_validate_f5_config(&general)
            .expect_err("oversized row is rejected during admission")
            .code(),
        "GPS-FMT-101"
    );
    fs::remove_dir_all(root).expect("cleanup");

    // La frontera de la fila es inclusiva y exacta en 65536 bytes contando el
    // LF, que es la capacidad reservada para el búfer: una fila de ese tamaño
    // se admite y la siguiente se rechaza, de modo que nunca hay realloc.
    let header = "region_id\tsequence_id\tstart_1based\tend_1based\n";
    let suffix = "\tchr1\t1\t1\n";
    for (label, row_bytes, accepted) in [("exact", 65_536_usize, true), ("excess", 65_537, false)] {
        let root = temporary(&format!("mask-row-{label}"));
        fs::create_dir_all(&root).expect("root");
        let identifier = "r".repeat(row_bytes - suffix.len());
        let mask = format!("{header}{identifier}{suffix}");
        assert_eq!(mask.len() - header.len(), row_bytes, "{label}");
        let general = configuration(&root, Some(&mask), "", "AAAA");
        let outcome = load_and_validate_f5_config(&general);
        if accepted {
            outcome.expect("a row of exactly the bounded size is admitted");
        } else {
            assert_eq!(
                outcome
                    .expect_err("one byte past the bound is rejected")
                    .code(),
                "GPS-FMT-101"
            );
        }
        fs::remove_dir_all(root).expect("cleanup");
    }
}

#[test]
fn site_and_cell_bounds_bins_and_checked_formula_overflow_are_accounted() {
    let root = temporary("site-bound");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration(
        &root,
        "",
        &"A".repeat(4_097),
        "length_min_bp\tlength_max_bp\tweight\n4097\t4097\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t4096\t0.5\t0.05\n",
    );
    assert_eq!(
        load_and_validate_f5_config(&general)
            .expect_err("S bound is rejected before reference I/O")
            .code(),
        "GPS-RES-101"
    );
    fs::remove_dir_all(root).expect("cleanup");

    let root = temporary("cell-bound");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration(
        &root,
        "",
        "AAAAAAAAAA",
        "length_min_bp\tlength_max_bp\tweight\n10\t10\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t9\t0.5\t0.05\n",
    );
    assert_eq!(
        load_and_validate_f5_config(&general)
            .expect_err("H*S bound is rejected during admission")
            .code(),
        "GPS-RES-101"
    );
    fs::remove_dir_all(root).expect("cleanup");

    let root = temporary("bins");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration(
        &root,
        "",
        "AAAAAA",
        "length_min_bp\tlength_max_bp\tweight\n3\t3\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t2\t0.5\t0.05\nP2\t1\t1\t0.4\t0.05\nP2\t2\t2\t0.3\t0.05\n",
    );
    replace_scenario(&root, "founder_count = 1", "founder_count = 2");
    replace_scenario(
        &root,
        "founder_ids = [\"I1\"]",
        "founder_ids = [\"I1\"]\n\n[[subpopulations]]\npopulation_id = \"P2\"\nindividual_count = 1\nfounder_ids = [\"I2\"]",
    );
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let result = preflight_f5_resources(&config, F5OtherResources::default()).expect("bins");
    assert_eq!(result.bin_count_total(), 3);
    assert_eq!(result.profile_rows(), 8);
    fs::remove_dir_all(root).expect("cleanup");

    let root = temporary("checked-overflow");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    for other in [
        F5OtherResources::new(u64::MAX, 0, 0),
        F5OtherResources::new(0, u64::MAX, 0),
        F5OtherResources::new(0, 0, u64::MAX),
    ] {
        assert_eq!(
            load_with_resources(&general, other)
                .expect_err("checked formula overflow")
                .code(),
            "GPS-RES-101"
        );
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn site_bound_covers_one_and_the_inclusive_4096_boundary() {
    let root = temporary("site-bound-one");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration_with_gap(
        &root,
        "",
        "AA",
        "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
        "gap_min_bp\tgap_max_bp\tweight\n2\t2\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
    );
    let config = load_and_validate_f5_config(&general).expect("S_bound=1 configuration");
    let admitted =
        preflight_f5_resources(&config, F5OtherResources::default()).expect("S_bound=1 preflight");
    assert_eq!(admitted.site_bound_per_block(), 1);
    fs::remove_dir_all(root).expect("cleanup");

    let root = temporary("site-bound-4096");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration(
        &root,
        "",
        &"A".repeat(4_096),
        "length_min_bp\tlength_max_bp\tweight\n4096\t4096\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t4095\t0.5\t0.05\n",
    );
    replace_scenario(
        &root,
        "max_haplotype_cells = 16",
        "max_haplotype_cells = 8192",
    );
    let error = load_and_validate_f5_config(&general)
        .expect_err("the independent global memory bound remains applicable");
    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(error.context()["dimension"], "working_memory_bytes");
    assert!(
        !error
            .render_human()
            .contains("site_bound_per_block supera 4096")
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[cfg(target_os = "linux")]
#[test]
fn large_fasta_with_small_budget_is_rejected_without_materializing_all_bases() {
    const CHILD_ROOT: &str = "GPS_F5_STREAMING_MEMORY_CHILD_ROOT";
    if let Some(root) = std::env::var_os(CHILD_ROOT) {
        let general = PathBuf::from(root).join("authoring/general.toml");
        let error = load_and_validate_f5_config(&general)
            .expect_err("declared memory budget must reject during admission");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(error.context()["dimension"], "working_memory_bytes");
        return;
    }

    let root = temporary("large-streaming-memory");
    fs::create_dir_all(&root).expect("root");
    configuration(&root, None, "max_working_memory_bytes = 65536\n", "AA");
    // Una sola línea fuerza al antiguo lector materializante a retener a la
    // vez la línea y el vector completo de bases. El hijo dispone de 128 MiB.
    replace_reference_with_repeated_base(&root, 80 * 1_024 * 1_024);

    let executable = std::env::current_exe().expect("test executable");
    let output = Command::new("sh")
        .args([
            "-c",
            "ulimit -v 131072; exec \"$1\" --exact large_fasta_with_small_budget_is_rejected_without_materializing_all_bases --nocapture",
            "f5-streaming-child",
        ])
        .arg(executable)
        .env(CHILD_ROOT, &root)
        .output()
        .expect("bounded child process");
    assert!(
        output.status.success(),
        "bounded child failed:\nstdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[cfg(target_os = "linux")]
#[test]
fn long_fasta_description_is_discarded_under_a_bounded_address_space() {
    const CHILD_ROOT: &str = "GPS_F5_LONG_HEADER_MEMORY_CHILD_ROOT";
    if let Some(root) = std::env::var_os(CHILD_ROOT) {
        let general = PathBuf::from(root).join("authoring/general.toml");
        let error = load_and_validate_f5_config(&general)
            .expect_err("declared memory budget must reject during admission");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(error.context()["dimension"], "working_memory_bytes");
        return;
    }

    let root = temporary("long-header-memory");
    fs::create_dir_all(&root).expect("root");
    configuration(&root, None, "max_working_memory_bytes = 65536\n", "AA");
    replace_reference_with_long_description(&root, 80 * 1_024 * 1_024);

    let executable = std::env::current_exe().expect("test executable");
    let output = Command::new("sh")
        .args([
            "-c",
            "ulimit -v 131072; exec \"$1\" --exact long_fasta_description_is_discarded_under_a_bounded_address_space --nocapture",
            "f5-long-header-child",
        ])
        .arg(executable)
        .env(CHILD_ROOT, &root)
        .output()
        .expect("bounded child process");
    assert!(
        output.status.success(),
        "bounded child failed:\nstdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn fasta_descriptions_raise_v8_operation_budget_and_preserve_productive_rng_path() {
    const DESCRIPTION_BYTES: usize = 80 * 1_024 * 1_024;
    let short_root = temporary("fasta-bytes-short");
    let long_root = temporary("fasta-bytes-long");
    fs::create_dir_all(&short_root).expect("short root");
    fs::create_dir_all(&long_root).expect("long root");

    let short_general = configuration(&short_root, None, "max_operation_units = 199999999\n", "AA");
    let long_general = configuration(&long_root, None, "max_operation_units = 199999999\n", "AA");
    replace_reference_with_long_description(&long_root, DESCRIPTION_BYTES);

    let short_config = load_and_validate_f5_config(&short_general).expect("short configuration");
    let short = preflight_f5_resources(&short_config, F5OtherResources::default())
        .expect("short preflight");
    let long_config = load_and_validate_f5_config(&long_general).expect("long configuration");
    let long =
        preflight_f5_resources(&long_config, F5OtherResources::default()).expect("long preflight");

    assert_eq!(
        short.reference_bases_processed(),
        long.reference_bases_processed()
    );
    assert_eq!(short.eligible_bases(), long.eligible_bases());
    assert_eq!(short.tract_count(), long.tract_count());
    assert_eq!(
        short.fasta_bytes_processed(),
        fs::metadata(short_root.join("authoring/references/ref.fa"))
            .expect("short FASTA metadata")
            .len()
    );
    assert_eq!(
        long.fasta_bytes_processed(),
        fs::metadata(long_root.join("authoring/references/ref.fa"))
            .expect("long FASTA metadata")
            .len()
    );
    assert_eq!(
        long.fasta_bytes_processed() - short.fasta_bytes_processed(),
        u64::try_from(DESCRIPTION_BYTES + 1).expect("description length fits u64")
    );
    assert_eq!(
        long.operation_units() - short.operation_units(),
        2 * (long.fasta_bytes_processed() - short.fasta_bytes_processed()),
        "the raw header is traversed by the preflight and productive REF recovery passes"
    );

    // El literal sustituido debe conservar la anchura del marcador: los bytes del
    // escenario entran en la cota como `2*scenario_toml_bytes`, de modo que una
    // anchura distinta desplazaría la frontera que esta prueba acredita.
    assert_eq!(long.operation_units().to_string().len(), 9);
    replace_scenario(
        &long_root,
        "max_haplotype_cells = 16\nmax_operation_units = 199999999\n",
        &format!(
            "max_haplotype_cells = 16\nmax_operation_units = {}\n",
            long.operation_units()
        ),
    );
    let exact_config = load_and_validate_f5_config(&long_general).expect("exact configuration");
    let exact = preflight_f5_resources(&exact_config, F5OtherResources::default())
        .expect("exact v8 operation limit is inclusive");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let short_layout = build_f5_physical_layout(
        F5LayoutRequest::new(&short_config, &short, &execution).expect("short request"),
    )
    .expect("short productive layout");
    let long_layout = build_f5_physical_layout(
        F5LayoutRequest::new(&exact_config, &exact, &execution).expect("long request"),
    )
    .expect("long productive layout");
    assert_eq!(
        long_layout, short_layout,
        "FASTA descriptions do not affect layout or its named RNG streams"
    );

    replace_scenario(
        &long_root,
        &format!(
            "max_haplotype_cells = 16\nmax_operation_units = {}\n",
            long.operation_units()
        ),
        &format!(
            "max_haplotype_cells = 16\nmax_operation_units = {}\n",
            long.operation_units() - 1
        ),
    );
    assert_eq!(
        load_and_validate_f5_config(&long_general)
            .expect_err("one operation below the long-header bound must reject")
            .code(),
        "GPS-RES-101"
    );

    fs::remove_dir_all(short_root).expect("short cleanup");
    fs::remove_dir_all(long_root).expect("long cleanup");
}

#[test]
fn mask_bytes_raise_v8_operation_budget_with_equal_geometry_and_exact_limit() {
    let short_root = temporary("mask-bytes-short");
    let long_root = temporary("mask-bytes-long");
    fs::create_dir_all(&short_root).expect("short root");
    fs::create_dir_all(&long_root).expect("long root");
    let short_mask = "region_id\tsequence_id\tstart_1based\tend_1based\na\tchr1\t3\t3\n";
    let long_mask = format!(
        "region_id\tsequence_id\tstart_1based\tend_1based\n{}\tchr1\t3\t3\n",
        "m".repeat(60_000)
    );
    let short_general = configuration(
        &short_root,
        Some(short_mask),
        "max_operation_units = 99999\n",
        "AAA",
    );
    let long_general = configuration(
        &long_root,
        Some(&long_mask),
        "max_operation_units = 99999\n",
        "AAA",
    );
    let short_config = load_and_validate_f5_config(&short_general).expect("short configuration");
    let short = preflight_f5_resources(&short_config, F5OtherResources::default())
        .expect("short preflight");
    let long_config = load_and_validate_f5_config(&long_general).expect("long configuration");
    let long =
        preflight_f5_resources(&long_config, F5OtherResources::default()).expect("long preflight");

    assert_eq!(
        short.reference_bases_processed(),
        long.reference_bases_processed()
    );
    assert_eq!(short.eligible_bases(), long.eligible_bases());
    assert_eq!(short.tract_count(), long.tract_count());
    assert_eq!(
        short.mask_bytes_processed(),
        fs::metadata(short_root.join("authoring/scenarios/excluded_regions.tsv"))
            .expect("short mask metadata")
            .len()
    );
    assert_eq!(
        long.mask_bytes_processed(),
        fs::metadata(long_root.join("authoring/scenarios/excluded_regions.tsv"))
            .expect("long mask metadata")
            .len()
    );
    assert_eq!(
        long.operation_units() - short.operation_units(),
        long.mask_bytes_processed() - short.mask_bytes_processed(),
        "row and region-ID bytes are scanned even when mask geometry is identical"
    );

    assert_eq!(long.operation_units().to_string().len(), 5);
    replace_scenario(
        &long_root,
        "max_haplotype_cells = 16\nmax_operation_units = 99999\n",
        &format!(
            "max_haplotype_cells = 16\nmax_operation_units = {}\n",
            long.operation_units()
        ),
    );
    let exact_config = load_and_validate_f5_config(&long_general).expect("exact configuration");
    preflight_f5_resources(&exact_config, F5OtherResources::default())
        .expect("exact mask-byte limit is inclusive");
    replace_scenario(
        &long_root,
        &format!(
            "max_haplotype_cells = 16\nmax_operation_units = {}\n",
            long.operation_units()
        ),
        &format!(
            "max_haplotype_cells = 16\nmax_operation_units = {}\n",
            long.operation_units() - 1
        ),
    );
    assert_eq!(
        load_and_validate_f5_config(&long_general)
            .expect_err("one operation below the mask-byte bound must reject")
            .code(),
        "GPS-RES-101"
    );

    fs::remove_dir_all(short_root).expect("short cleanup");
    fs::remove_dir_all(long_root).expect("long cleanup");
}

#[test]
fn literal_byte_oracle_accepts_exact_operation_formula_and_rejects_minus_one() {
    let mask = b"region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t3\t3\n";
    for (label, mask) in [("without-mask", None), ("with-mask", Some(mask.as_slice()))] {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("root");
        let general = configuration(
            &root,
            mask.map(|bytes| std::str::from_utf8(bytes).expect("UTF-8 mask")),
            "max_operation_units = 999999\n",
            "AAA",
        );

        let baseline_config =
            load_and_validate_f5_config(&general).expect("high-limit literal fixture");
        let baseline = preflight_f5_resources(&baseline_config, F5OtherResources::default())
            .expect("high-limit public preflight");
        assert_eq!(
            baseline.operation_units(),
            literal_fixture_operation_oracle(&root, mask),
            "term-by-term literal oracle differs before the exact-bound rewrite: {label}"
        );

        let mut candidate = 999_999_u64;
        let exact = loop {
            set_scenario_operation_limit(&root, candidate);
            let observed = literal_fixture_operation_oracle(&root, mask);
            if observed == candidate {
                break observed;
            }
            candidate = observed;
        };
        let exact_config = load_and_validate_f5_config(&general).expect("exact literal fixture");
        let exact_preflight = preflight_f5_resources(&exact_config, F5OtherResources::default())
            .expect("exact public preflight");
        assert_eq!(exact_preflight.operation_units(), exact, "{label}");
        assert_eq!(
            exact_preflight.fasta_bytes_processed(),
            b">chr1\nAAA\n".len() as u64,
            "{label}"
        );
        assert_eq!(
            exact_preflight.mask_bytes_processed(),
            mask.map_or(0, |bytes| bytes.len() as u64),
            "{label}"
        );

        assert_eq!(
            exact.to_string().len(),
            (exact - 1).to_string().len(),
            "the minus-one fixture must preserve the scenario byte count"
        );
        set_scenario_operation_limit(&root, exact - 1);
        let error = load_and_validate_f5_config(&general)
            .expect_err("one operation below the literal formula must reject");
        assert_eq!(error.code(), "GPS-RES-101", "{label}");
        assert_eq!(error.context()["dimension"], "operation_units", "{label}");
        fs::remove_dir_all(root).expect("cleanup");
    }
}

#[test]
fn valid_mask_identity_changes_are_input_mismatches_before_entropy_or_destination() {
    let declared = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\n";
    let changed = [
        (
            "shorter",
            b"region_id\tsequence_id\tstart_1based\tend_1based\n".as_slice(),
        ),
        (
            "longer",
            b"region_id\tsequence_id\tstart_1based\tend_1based\nmask-longer\tchr1\t2\t2\n"
                .as_slice(),
        ),
    ];
    for (label, changed) in changed {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("root");
        let general = configuration(&root, Some(declared), "", "AAAA");
        remove_declared_seed(&root);
        fs::write(
            root.join("authoring/scenarios/excluded_regions.tsv"),
            changed,
        )
        .expect("mask changed without updating its declared digest");

        let destination = root.join("must-not-exist");
        let mut entropy = CountingEntropy::default();
        let error = load_and_validate_f5_config_with_entropy(&general, &mut entropy)
            .expect_err("stale mask identity");
        assert_eq!(error.code(), "GPS-INPUT-002", "{label}");
        assert_eq!(
            error.context()["expected_digest"],
            sha256(declared.as_bytes())
        );
        assert_eq!(error.context()["observed_digest"], sha256(changed));
        assert_eq!(entropy.calls, 0, "{label}");
        let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
            .args([
                "validate",
                general.to_str().expect("UTF-8 path"),
                "--output",
                destination.to_str().expect("UTF-8 destination"),
            ])
            .output()
            .expect("CLI validation of stale mask identity");
        assert_eq!(output.status.code(), Some(2), "{label}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("GPS-INPUT-002"),
            "{label}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(!destination.exists(), "{label}");
        fs::remove_dir_all(root).expect("cleanup");
    }
}

#[cfg(target_os = "linux")]
#[test]
fn fragmented_tracts_are_budgeted_before_vector_growth() {
    const CHILD_ROOT: &str = "GPS_F5_FRAGMENTED_TRACTS_CHILD_ROOT";
    if let Some(root) = std::env::var_os(CHILD_ROOT) {
        let general = PathBuf::from(root).join("authoring/general.toml");
        let error = load_and_validate_f5_config(&general)
            .expect_err("scan-time memory lower bound must reject");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(error.context()["dimension"], "working_memory_bytes");
        return;
    }

    let root = temporary("fragmented-tract-memory");
    fs::create_dir_all(&root).expect("root");
    configuration(&root, None, "max_working_memory_bytes = 131072\n", "AA");
    replace_reference_with_alternating_bases(&root, 16 * 1_024 * 1_024);

    let executable = std::env::current_exe().expect("test executable");
    let output = Command::new("sh")
        .args([
            "-c",
            "ulimit -v 131072; exec \"$1\" --exact fragmented_tracts_are_budgeted_before_vector_growth --nocapture",
            "f5-fragmented-tract-child",
        ])
        .arg(executable)
        .env(CHILD_ROOT, &root)
        .output()
        .expect("bounded child process");
    assert!(
        output.status.success(),
        "bounded child failed:\nstdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[cfg(target_os = "linux")]
#[test]
fn maximum_length_sequence_ids_are_budgeted_before_metadata_growth() {
    const CHILD_ROOT: &str = "GPS_F5_MANY_SEQUENCE_IDS_CHILD_ROOT";
    if let Some(root) = std::env::var_os(CHILD_ROOT) {
        let general = PathBuf::from(root).join("authoring/general.toml");
        let error = load_and_validate_f5_config(&general)
            .expect_err("scan-time ID memory lower bound must reject");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(error.context()["dimension"], "working_memory_bytes");
        assert_eq!(error.context()["accumulated"], "131328");
        assert_eq!(error.context()["limit"], "131072");
        return;
    }

    let root = temporary("many-sequence-id-memory");
    fs::create_dir_all(&root).expect("root");
    configuration(&root, None, "max_working_memory_bytes = 131072\n", "AA");
    replace_reference_with_many_maximum_ids(&root, 750_000, 'N');

    let executable = std::env::current_exe().expect("test executable");
    let output = Command::new("sh")
        .args([
            "-c",
            "ulimit -v 131072; exec \"$1\" --exact maximum_length_sequence_ids_are_budgeted_before_metadata_growth --nocapture",
            "f5-many-ids-child",
        ])
        .arg(executable)
        .env(CHILD_ROOT, &root)
        .output()
        .expect("bounded child process");
    assert!(
        output.status.success(),
        "bounded child failed:\nstdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn retained_sequence_id_peak_is_persisted_in_the_published_memory_bound() {
    let root = temporary("retained-sequence-id-result");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AA");
    replace_reference_with_n_only_ids_and_eligible_tail(&root, 128);
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let result = preflight_f5_resources(&config, F5OtherResources::default())
        .expect("N-only IDs plus one eligible sequence");

    assert_eq!(result.retained_sequence_ids(), 129);
    assert_eq!(result.tract_count(), 1);
    assert_eq!(
        result.working_memory_bytes(),
        final_formula_memory_oracle(FinalFormulaMemoryInputs {
            block_count: 1,
            tract_count: 1,
            retained_sequence_ids: 129,
            profile_rows: 6,
            normalization_rows: 5,
            population_count: 1,
            haplotype_count: 2,
            site_bound: 2,
            bin_count_total: 1,
            founder_copy_count: 2,
            sequence_count: 129,
        })
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn resource_limit_includes_other_outputs_and_rejects_before_any_destination() {
    let root = temporary("limit");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let other = F5OtherResources::new(7, 11, 13);
    let config = load_with_resources(&general, other).expect("configuration");
    let baseline = preflight_f5_resources(&config, other).expect("baseline");
    assert_eq!(baseline.other_resources(), other);
    let rejected_general = configuration(
        &root.join("rejected"),
        None,
        &format!("max_output_bytes = {}\n", baseline.output_bytes() - 1),
        "AAAA",
    );
    let error = load_with_resources(&rejected_general, other)
        .expect_err("output budget must reject before returning a configuration");
    assert_eq!(error.code(), "GPS-RES-101");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn hard_resource_bounds_are_inclusive_and_the_next_byte_is_rejected() {
    let root = temporary("hard-bounds");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let base =
        preflight_f5_resources(&config, F5OtherResources::default()).expect("base preflight");
    let exact = F5OtherResources::new(
        F5_HARD_OPERATION_LIMIT - base.operation_units(),
        F5_HARD_WORKING_MEMORY_LIMIT - base.working_memory_bytes(),
        F5_HARD_OUTPUT_LIMIT - base.output_bytes(),
    );
    let exact_config = load_with_resources(&general, exact).expect("inclusive bounds");
    let admitted = preflight_f5_resources(&exact_config, exact).expect("bound preflight");
    assert_eq!(admitted.operation_units(), F5_HARD_OPERATION_LIMIT);
    assert_eq!(
        admitted.working_memory_bytes(),
        F5_HARD_WORKING_MEMORY_LIMIT
    );
    assert_eq!(admitted.output_bytes(), F5_HARD_OUTPUT_LIMIT);
    let error = preflight_f5_resources(
        &exact_config,
        F5OtherResources::new(
            exact.operation_units(),
            exact.working_memory_bytes(),
            exact.output_bytes() + 1,
        ),
    )
    .expect_err("one byte beyond the hard limit must fail");
    assert_eq!(error.code(), "GPS-CONF-101");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn zero_fractional_target_is_a_domain_failure_without_rng() {
    let root = temporary("zero-target");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let scenario = root.join("authoring/scenarios/scenario.toml");
    let source = fs::read_to_string(&scenario)
        .expect("scenario text")
        .replace(
            "block_count = 1",
            "coverage_fraction = 0.000000000000000001",
        );
    fs::write(&scenario, source).expect("scenario rewrite");
    let error = load_and_validate_f5_config(&general).expect_err("T=0 must fail during admission");
    assert_eq!(error.code(), "GPS-DOMAIN-101");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn block_length_plan_rejects_a_preflight_crossed_with_another_profile() {
    let root = temporary("crossed-length-profile");
    fs::create_dir_all(&root).expect("root");
    let first_ld_profile = "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n";
    let second_ld_profile = "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t2\t0.5\t0.05\n";
    let haplotypes = "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n";
    let first_path = custom_configuration(
        &root.join("first"),
        "",
        "AAA",
        "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
        haplotypes,
        first_ld_profile,
    );
    let second_path = custom_configuration(
        &root.join("second"),
        "",
        "AAA",
        "length_min_bp\tlength_max_bp\tweight\n3\t3\t1\n",
        haplotypes,
        second_ld_profile,
    );
    let first = load_and_validate_f5_config(&first_path).expect("first configuration");
    let second = load_and_validate_f5_config(&second_path).expect("second configuration");
    let first_preflight =
        preflight_f5_resources(&first, F5OtherResources::default()).expect("first preflight");

    F5BlockLengthPlan::from_preflight(&first, &first_preflight)
        .expect("the bound configuration and preflight must be accepted");
    let error = F5BlockLengthPlan::from_preflight(&second, &first_preflight)
        .expect_err("a preflight from another profile must be rejected");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(error.context()["fields"], "preflight.binding");
    assert!(error.context()["detail"].contains("no pertenece"));
    assert_ne!(
        first_preflight.profile_metadata_sha256(),
        second.declared_scenario().profile_metadata_sha256.value()
    );
    assert_eq!(first_preflight.block_length_minimum(), 2);
    assert_eq!(first_preflight.block_length_maximum(), 2);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn block_length_plan_binds_same_profile_preflights_to_the_declared_objective() {
    let root = temporary("crossed-length-target");
    fs::create_dir_all(&root).expect("root");
    let block_path = configuration(&root.join("block"), None, "", "AAAA");
    let covered_path = configuration(&root.join("covered"), None, "", "AAAA");
    let fraction_path = configuration(&root.join("fraction"), None, "", "AAAA");
    let equivalent_fraction_path =
        configuration(&root.join("equivalent-fraction"), None, "", "AAAA");
    replace_scenario(
        &root.join("covered"),
        "block_count = 1",
        "covered_bases = 2",
    );
    replace_scenario(
        &root.join("fraction"),
        "block_count = 1",
        "coverage_fraction = 0.5",
    );
    replace_scenario(
        &root.join("equivalent-fraction"),
        "block_count = 1",
        "coverage_fraction = 0.50",
    );

    let block = load_and_validate_f5_config(&block_path).expect("block configuration");
    let covered = load_and_validate_f5_config(&covered_path).expect("covered configuration");
    let fraction = load_and_validate_f5_config(&fraction_path).expect("fraction configuration");
    let equivalent_fraction = load_and_validate_f5_config(&equivalent_fraction_path)
        .expect("equivalent fraction configuration");
    assert_eq!(
        block.declared_scenario().profile_metadata_sha256.value(),
        covered.declared_scenario().profile_metadata_sha256.value(),
        "the test must cross only objectives, not profiles"
    );

    let block_preflight =
        preflight_f5_resources(&block, F5OtherResources::default()).expect("block preflight");
    let fraction_preflight =
        preflight_f5_resources(&fraction, F5OtherResources::default()).expect("fraction preflight");
    for crossed in [&covered, &fraction] {
        let error = F5BlockLengthPlan::from_preflight(crossed, &block_preflight)
            .expect_err("a same-profile preflight must not cross target kinds");
        assert_eq!(error.code(), "GPS-CONF-101");
        assert_eq!(error.context()["fields"], "preflight.binding");
    }
    let error = F5BlockLengthPlan::from_preflight(&covered, &fraction_preflight)
        .expect_err("covered_bases and coverage_fraction remain distinct objectives");
    assert_eq!(error.code(), "GPS-CONF-101");
    F5BlockLengthPlan::from_preflight(&equivalent_fraction, &fraction_preflight)
        .expect("equivalent decimal spellings bind to the same exact rational");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn wide_literal_length_support_is_rejected_by_preflight_before_enumeration() {
    let root = temporary("wide-length-support");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration_with_gap(
        &root,
        "",
        "AAA",
        "length_min_bp\tlength_max_bp\tweight\n2\t134217729\t1\n",
        "gap_min_bp\tgap_max_bp\tweight\n134217729\t134217729\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t134217728\t0.5\t0.05\n",
    );
    let error = load_and_validate_f5_config(&general)
        .expect_err("the literal two-pass support cannot exceed the operation budget");
    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(error.context()["dimension"], "operation_units");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn public_layout_preserves_configured_scenario_in_plan_subseed_and_ids() {
    let root = temporary("configured-scenario-identity");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    rename_scenario(&root, "scenario", "configured-scenario");
    let configured = load_and_validate_f5_config(&general).expect("renamed configuration");
    let preflight = preflight_f5_resources(&configured, F5OtherResources::default())
        .expect("renamed preflight");
    let plan = F5BlockLengthPlan::from_preflight(&configured, &preflight).expect("bound plan");
    assert_eq!(plan.scenario_id().as_str(), "configured-scenario");

    let root_seed = RootSeed::parse("seed", configured.seed()).expect("root seed");
    let scenario = ScenarioId::parse(configured.scenario_id(), "test").expect("scenario");
    let reference = ReferenceId::parse(configured.reference_id(), "test").expect("reference");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let plan_subseed = derive_f5_subseed(
        root_seed,
        F5RandomDecision::layout_length(plan.scenario_id(), 0, 0),
    )
    .expect("plan subseed");
    let configured_subseed =
        derive_f5_subseed(root_seed, F5RandomDecision::layout_length(&scenario, 0, 0))
            .expect("configured subseed");
    assert_eq!(plan_subseed, configured_subseed);

    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&configured, &preflight, &execution).expect("request"),
    )
    .expect("layout");
    let expected_block = F5IdentityRegistry::default()
        .derive_block(
            F5IdentityContext::new(root_seed, &execution, &scenario, &reference),
            0,
        )
        .expect("independent block ID");
    assert_eq!(layout.blocks()[0].block_id(), &expected_block);

    rename_scenario(&root, "configured-scenario", "different-scenario");
    let different = load_and_validate_f5_config(&general).expect("second renamed configuration");
    let different_preflight =
        preflight_f5_resources(&different, F5OtherResources::default()).expect("second preflight");
    let different_plan =
        F5BlockLengthPlan::from_preflight(&different, &different_preflight).expect("second plan");
    let different_subseed = derive_f5_subseed(
        root_seed,
        F5RandomDecision::layout_length(different_plan.scenario_id(), 0, 0),
    )
    .expect("different subseed");
    let different_layout = build_f5_physical_layout(
        F5LayoutRequest::new(&different, &different_preflight, &execution).expect("second request"),
    )
    .expect("second layout");
    assert_ne!(plan_subseed, different_subseed);
    assert_ne!(
        layout.blocks()[0].block_id(),
        different_layout.blocks()[0].block_id(),
        "scenario_id participates in the public block identity"
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn gap_support_cardinality_is_discontinuous_bounded_and_inclusive_before_rng() {
    let root = temporary("gap-support-budget");
    fs::create_dir_all(&root).expect("root");
    let other = F5OtherResources::new(100_000, 0, 0);
    let general = custom_configuration_with_gap(
        &root,
        "max_operation_units = 999999\n",
        "AAAA",
        "length_min_bp\tlength_max_bp\tweight\n4\t4\t1\n",
        "gap_min_bp\tgap_max_bp\tweight\n1\t2\t0.5\n10\t12\t0.5\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t3\t0.5\t0.05\n",
    );
    let unrestricted = load_with_resources(&general, other).expect("discontinuous support");
    let baseline = preflight_f5_resources(&unrestricted, other).expect("bounded preflight");
    assert_eq!(baseline.gap_support_cardinality(), 5);
    assert_eq!(baseline.site_bound_per_block(), 4);
    let maximum_gap_draws =
        baseline.max_layout_attempts() * baseline.block_count() * baseline.site_bound_per_block();
    assert_eq!(
        baseline.gap_selector_operation_units(),
        2 * maximum_gap_draws * baseline.gap_support_cardinality()
    );
    assert_eq!(baseline.operation_units().to_string().len(), 6);

    replace_scenario(
        &root,
        "max_operation_units = 999999\n",
        &format!("max_operation_units = {}\n", baseline.operation_units()),
    );
    let exact = load_with_resources(&general, other).expect("exact configuration");
    assert_eq!(
        preflight_f5_resources(&exact, other)
            .expect("inclusive operation bound")
            .operation_units(),
        baseline.operation_units()
    );
    replace_scenario(
        &root,
        &format!("max_operation_units = {}\n", baseline.operation_units()),
        &format!("max_operation_units = {}\n", baseline.operation_units() - 1),
    );
    let error = load_with_resources(&general, other)
        .expect_err("gap traversal bound rejects before any layout RNG");
    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(error.context()["dimension"], "operation_units");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn public_v8_preflight_serializes_gap_terms_and_validates_both_contracts() {
    let root = temporary("public-v8-projection");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration_with_gap(
        &root,
        "max_operation_units = 999999\n",
        "AAAA",
        "length_min_bp\tlength_max_bp\tweight\n4\t4\t1\n",
        "gap_min_bp\tgap_max_bp\tweight\n1\t2\t0.5\n10\t12\t0.5\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t3\t0.5\t0.05\n",
    );
    let config = load_and_validate_f5_config(&general).expect("public configuration");
    let preflight = preflight_f5_resources(&config, F5OtherResources::default())
        .expect("public resource calculation");
    assert_eq!(
        preflight.algorithm(),
        genome_population_simulator::f5_admission::F5_RESOURCE_ACCOUNTING_ALGORITHM
    );
    assert_eq!(preflight.algorithm(), "f5-resource-accounting-v8");

    let projection = preflight.to_contract_json();
    assert_f5_1_2_5_v8_contract_key_set(&projection);
    assert_eq!(projection["algorithm"], preflight.algorithm());
    assert_eq!(
        projection["gap_support_cardinality"],
        preflight.gap_support_cardinality()
    );
    assert_eq!(
        projection["gap_selector_operation_units"],
        preflight.gap_selector_operation_units()
    );
    assert_preflight_projection_validates_through_effective_and_manifest(&root, &projection);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn trillion_wide_gap_support_is_rejected_without_enumeration() {
    let root = temporary("trillion-gap-support");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration_with_gap(
        &root,
        "",
        "AAAA",
        "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
        "gap_min_bp\tgap_max_bp\tweight\n1\t1000000000000\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
    );
    let error = load_and_validate_f5_config(&general)
        .expect_err("the literal selector cost exceeds the hard operation bound");
    assert_eq!(error.code(), "GPS-RES-101");
    assert_eq!(error.context()["dimension"], "operation_units");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn maximum_singleton_gap_exhausts_layout_normally_without_overflow() {
    let root = temporary("maximum-singleton-gap");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration_with_gap(
        &root,
        "",
        "AA",
        "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
        "gap_min_bp\tgap_max_bp\tweight\n18446744073709551615\t18446744073709551615\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
    );
    let config = load_and_validate_f5_config(&general).expect("maximum singleton gap");
    let preflight = preflight_f5_resources(&config, F5OtherResources::default())
        .expect("singleton support stays bounded");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let error = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &preflight, &execution).expect("request"),
    )
    .expect_err("one-site attempt is insufficient");
    assert_eq!(error.code(), "GPS-DOMAIN-106");
    assert!(!error.diagnostic().render_human().contains("desborda"));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn disjoint_gap_profile_flows_from_preflight_to_reproducible_layout() {
    let root = temporary("disjoint-gaps-layout");
    fs::create_dir_all(&root).expect("root");
    let general = custom_configuration_with_gap(
        &root,
        "",
        "AAAA",
        "length_min_bp\tlength_max_bp\tweight\n4\t4\t1\n",
        "gap_min_bp\tgap_max_bp\tweight\n1\t1\t0.5\n3\t3\t0.5\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t3\t0.5\t0.05\n",
    );
    let config = load_and_validate_f5_config(&general).expect("disjoint gaps profile");
    let preflight = preflight_f5_resources(&config, F5OtherResources::default())
        .expect("preflight accepts disjoint gaps without GPS-DOMAIN-101");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let first = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &preflight, &execution).expect("bound request"),
    )
    .expect("layout");
    let second = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &preflight, &execution).expect("bound request"),
    )
    .expect("same layout");
    assert_eq!(first, second, "the exact root seed scopes the gap RNG");
    assert_eq!(first.blocks().len(), 1);
    assert!(first.blocks()[0].sites().len() >= 2);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn layout_request_rejects_preflight_from_other_reference_before_layout_rng() {
    let first_root = temporary("layout-preflight-reference-a");
    let second_root = temporary("layout-preflight-reference-b");
    fs::create_dir_all(&first_root).expect("first root");
    fs::create_dir_all(&second_root).expect("second root");
    let first_general = configuration(&first_root, None, "", "AAAA");
    let second_general = configuration(&second_root, None, "", "CCCC");
    let first = load_and_validate_f5_config(&first_general).expect("first config");
    let preflight =
        preflight_f5_resources(&first, F5OtherResources::default()).expect("first preflight");
    let second = load_and_validate_f5_config(&second_general).expect("second config");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let error = F5LayoutRequest::new(&second, &preflight, &execution)
        .expect_err("a preflight cannot cross a reference SHA boundary");
    assert_eq!(error.code(), "GPS-CONF-101");
    replace_scenario(
        &first_root,
        "max_layout_attempts = 1",
        "max_layout_attempts = 2",
    );
    let changed_attempt_limit = load_and_validate_f5_config(&first_general)
        .expect("configuration with a different attempt limit");
    let error = F5LayoutRequest::new(&changed_attempt_limit, &preflight, &execution)
        .expect_err("a preflight cannot be reused with more layout attempts");
    assert_eq!(error.code(), "GPS-CONF-101");
    fs::remove_dir_all(first_root).expect("first cleanup");
    fs::remove_dir_all(second_root).expect("second cleanup");
}

#[test]
fn preflight_binding_rejects_equal_capacity_masks_and_changed_resource_limits() {
    let first_root = temporary("binding-mask-first");
    let second_root = temporary("binding-mask-second");
    let limits_root = temporary("binding-limits");
    fs::create_dir_all(&first_root).expect("first root");
    fs::create_dir_all(&second_root).expect("second root");
    fs::create_dir_all(&limits_root).expect("limits root");
    let first_mask = "region_id\tsequence_id\tstart_1based\tend_1based\na\tchr1\t1\t2\n";
    let second_mask = "region_id\tsequence_id\tstart_1based\tend_1based\nb\tchr1\t5\t6\n";
    let first_general = configuration(&first_root, Some(first_mask), "", "AAAAAA");
    let second_general = configuration(&second_root, Some(second_mask), "", "AAAAAA");
    let first = load_and_validate_f5_config(&first_general).expect("first masked configuration");
    let second = load_and_validate_f5_config(&second_general).expect("second masked configuration");
    let first_preflight =
        preflight_f5_resources(&first, F5OtherResources::default()).expect("first preflight");
    let second_preflight =
        preflight_f5_resources(&second, F5OtherResources::default()).expect("second preflight");
    assert_eq!(
        first_preflight.eligible_bases(),
        second_preflight.eligible_bases()
    );
    assert_eq!(
        first_preflight.tract_count(),
        second_preflight.tract_count()
    );
    assert_ne!(
        first_preflight.eligible_tracts(),
        second_preflight.eligible_tracts(),
        "the regression needs different masks with the same accounting totals"
    );
    let error = F5BlockLengthPlan::from_preflight(&second, &first_preflight)
        .expect_err("a preflight must retain the exact mask identity");
    assert_eq!(error.code(), "GPS-CONF-101");

    let limits_general = configuration(&limits_root, None, "max_operation_units = 9999\n", "AAAA");
    let unrestricted = load_and_validate_f5_config(&limits_general).expect("unrestricted config");
    let unrestricted_preflight = preflight_f5_resources(&unrestricted, F5OtherResources::default())
        .expect("unrestricted preflight");
    replace_scenario(
        &limits_root,
        "max_haplotype_cells = 16\nmax_operation_units = 9999\n",
        &format!(
            "max_haplotype_cells = 16\nmax_operation_units = {}\n",
            unrestricted_preflight.operation_units()
        ),
    );
    let limited = load_and_validate_f5_config(&limits_general).expect("same geometry, own limit");
    preflight_f5_resources(&limited, F5OtherResources::default())
        .expect("the exact operation bound remains admitted");
    let error = F5BlockLengthPlan::from_preflight(&limited, &unrestricted_preflight)
        .expect_err("a preflight must retain its operation limit");
    assert_eq!(error.code(), "GPS-CONF-101");

    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let error = F5LayoutRequest::new(&second, &first_preflight, &execution)
        .expect_err("request must fail before opening layout RNG for a crossed mask");
    assert_eq!(error.code(), "GPS-CONF-101");
    let error = F5LayoutRequest::new(&limited, &unrestricted_preflight, &execution)
        .expect_err("request must fail before opening layout RNG for a crossed limit");
    assert_eq!(error.code(), "GPS-CONF-101");
    fs::remove_dir_all(first_root).expect("first cleanup");
    fs::remove_dir_all(second_root).expect("second cleanup");
    fs::remove_dir_all(limits_root).expect("limits cleanup");
}

#[test]
fn preflight_identity_binds_seed_and_seed_origin_before_layout_rng() {
    let root = temporary("binding-seed-and-origin");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "max_working_memory_bytes = 99999\n", "AAAA");
    let first = load_and_validate_f5_config(&general).expect("declared seed configuration");
    let preflight =
        preflight_f5_resources(&first, F5OtherResources::default()).expect("first preflight");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");

    let general_source = fs::read_to_string(&general).expect("general source");
    fs::write(
        &general,
        general_source.replace(&"a".repeat(64), &"b".repeat(64)),
    )
    .expect("changed seed");
    let changed_seed = load_and_validate_f5_config(&general).expect("changed seed configuration");
    let error = F5LayoutRequest::new(&changed_seed, &preflight, &execution)
        .expect_err("a preflight must not cross an effective seed");
    assert_eq!(error.code(), "GPS-CONF-101");

    fs::write(
        &general,
        general_source.replace(&format!("seed = \"{}\"\n", "a".repeat(64)), ""),
    )
    .expect("generated seed configuration");
    let generated =
        load_and_validate_f5_config_with_entropy(&general, &mut FixedEntropy([0xaa; 32]))
            .expect("same seed bytes with generated origin");
    assert_eq!(generated.seed(), first.seed());
    assert_ne!(generated.seed_origin(), first.seed_origin());
    let error = F5LayoutRequest::new(&generated, &preflight, &execution)
        .expect_err("a preflight must not cross seed provenance");
    assert_eq!(error.code(), "GPS-CONF-101");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn preflight_identity_binds_semantic_subpopulations_and_overrides_before_layout_rng() {
    let root = temporary("binding-subpopulations-and-overrides");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let first = load_and_validate_f5_config(&general).expect("base configuration");
    let preflight =
        preflight_f5_resources(&first, F5OtherResources::default()).expect("base preflight");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let scenario = root.join("authoring/scenarios/scenario.toml");

    let source = fs::read_to_string(&scenario).expect("scenario source");
    fs::write(
        &scenario,
        source.replace(
            "population_id = \"P1\"\nindividual_count = 1\nfounder_ids = [\"I1\"]\n",
            "population_id = \"P1\"\nindividual_count = 1\nsampling = \"probabilistic\"\ninbreeding_f = 0.25\nfounder_ids = [\"I2\"]\n",
        ),
    )
    .expect("semantic subpopulation change");
    let changed_subpopulation =
        load_and_validate_f5_config(&general).expect("changed subpopulation configuration");
    let error = F5LayoutRequest::new(&changed_subpopulation, &preflight, &execution)
        .expect_err("a preflight must bind sampling, F and founders");
    assert_eq!(error.code(), "GPS-CONF-101");

    fs::write(
        &scenario,
        format!(
            "{source}\n[[sampling_overrides]]\npopulation_id = \"P1\"\nblock_index = 0\nmode = \"controlled\"\n"
        ),
    )
    .expect("sampling override");
    let changed_override = load_and_validate_f5_config(&general).expect("override configuration");
    let error = F5LayoutRequest::new(&changed_override, &preflight, &execution)
        .expect_err("a preflight must bind override population, block, mode and tolerance");
    assert_eq!(error.code(), "GPS-CONF-101");

    fs::write(
        &scenario,
        source.replace(
            "founder_ids = [\"I1\"]\n",
            "sampling = \"probabilistic\"\ninbreeding_f = 0.0\nfounder_ids = [\"I1\"]\n",
        ),
    )
    .expect("equivalent effective defaults");
    let equivalent = load_and_validate_f5_config(&general).expect("equivalent configuration");
    // La mitad semántica del binding sigue siendo canónica: el plan de
    // longitudes acepta la reescritura con los valores por omisión explícitos.
    F5BlockLengthPlan::from_preflight(&equivalent, &preflight)
        .expect("canonical effective identity accepts equivalent defaults");
    // La mitad física, en cambio, no: el escenario reescrito tiene otro SHA y
    // F5.1.2.5 exige las dos mitades para abrir layout.
    let error = F5LayoutRequest::new(&equivalent, &preflight, &execution)
        .expect_err("the sealed receipt binds the exact bytes of the scenario it admitted");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(error.context()["fields"], "preflight.binding");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn layout_request_binds_other_resource_contributions() {
    let root = temporary("binding-other-resources");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let other = F5OtherResources::new(7, 11, 13);
    let config = load_with_resources(&general, other).expect("configuration");
    let preflight = preflight_f5_resources(&config, other).expect("preflight with other costs");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let error = F5LayoutRequest::new(&config, &preflight, &execution)
        .expect_err("the default request must not erase non-F5 contributions");
    assert_eq!(error.code(), "GPS-CONF-101");
    let error = F5LayoutRequest::new_with_other_resources(
        &config,
        &preflight,
        F5OtherResources::new(8, 11, 13),
        &execution,
    )
    .expect_err("a distinct external operation contribution must be rejected");
    assert_eq!(error.code(), "GPS-CONF-101");
    F5LayoutRequest::new_with_other_resources(&config, &preflight, other, &execution)
        .expect("the exact external contributions remain bound");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn layout_streams_large_fasta_and_keeps_the_admitted_working_bound() {
    let root = temporary("stream-large-fasta");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "max_working_memory_bytes = 99999\n", "AAAA");
    let reference_bases = 2 * 1024 * 1024;
    replace_reference_with_repeated_base(&root, reference_bases);
    let unrestricted = load_and_validate_f5_config(&general).expect("large configuration");
    let baseline = preflight_f5_resources(&unrestricted, F5OtherResources::default())
        .expect("large-reference preflight");
    assert!(
        usize::try_from(baseline.working_memory_bytes()).expect("bound usize")
            < reference_bases / 16,
        "the admitted working bound is intentionally far below a materialized reference"
    );
    replace_scenario(
        &root,
        "max_haplotype_cells = 16\nmax_working_memory_bytes = 99999\n",
        &format!(
            "max_haplotype_cells = 16\nmax_working_memory_bytes = {}\n",
            baseline.working_memory_bytes()
        ),
    );
    let config = load_and_validate_f5_config(&general).expect("working-bound configuration");
    let preflight = preflight_f5_resources(&config, F5OtherResources::default())
        .expect("the exact low working bound is admitted");
    assert_eq!(
        preflight.working_memory_bytes(),
        baseline.working_memory_bytes()
    );
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &preflight, &execution).expect("bound request"),
    )
    .expect("streaming REF recovery must not materialize the large sequence");
    assert_eq!(layout.blocks().len(), 1);
    assert!(
        layout.blocks()[0]
            .sites()
            .iter()
            .all(|site| site.reference_base() == b'A')
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn profile_padding_changes_its_phase_without_summing_noncoexistent_peaks() {
    let short_root = temporary("profile-input-short");
    let long_root = temporary("profile-input-long");
    fs::create_dir_all(&short_root).expect("short root");
    fs::create_dir_all(&long_root).expect("long root");
    let short_general = configuration(
        &short_root,
        None,
        "max_working_memory_bytes = 99999\n",
        "AA",
    );
    let long_general = configuration(&long_root, None, "max_working_memory_bytes = 99999\n", "AA");
    replace_profile_source_citation(&long_root, &"c".repeat(4_096));

    let short_config = load_and_validate_f5_config(&short_general).expect("short configuration");
    let short = preflight_f5_resources(&short_config, F5OtherResources::default())
        .expect("short preflight");
    let long_config = load_and_validate_f5_config(&long_general).expect("long configuration");
    let long =
        preflight_f5_resources(&long_config, F5OtherResources::default()).expect("long preflight");
    assert_eq!(short.eligible_bases(), long.eligible_bases());
    assert_eq!(short.profile_rows(), long.profile_rows());
    assert_eq!(
        long.profile_metadata_bytes() - short.profile_metadata_bytes(),
        4_092,
        "the citation replaces the four-byte fixture literal"
    );
    assert_eq!(
        long.input_parser_peak_bytes() - short.input_parser_peak_bytes(),
        8_184,
        "the metadata phase reserves raw bytes plus its parser representation"
    );
    assert_eq!(
        long.working_memory_bytes(),
        short.working_memory_bytes(),
        "the unchanged final-formula phase remains the global maximum"
    );

    replace_scenario(
        &long_root,
        "max_working_memory_bytes = 99999\n",
        &format!(
            "max_working_memory_bytes = {}\n",
            long.working_memory_bytes()
        ),
    );
    let exact_config = load_and_validate_f5_config(&long_general).expect("exact config");
    let exact = preflight_f5_resources(&exact_config, F5OtherResources::default())
        .expect("the exact v8 memory limit is inclusive");
    assert_eq!(exact.working_memory_bytes(), long.working_memory_bytes());
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    assert_eq!(
        build_f5_physical_layout(
            F5LayoutRequest::new(&short_config, &short, &execution).expect("short request"),
        )
        .expect("short layout"),
        build_f5_physical_layout(
            F5LayoutRequest::new(&exact_config, &exact, &execution).expect("long request"),
        )
        .expect("long layout"),
        "admitted source citations do not affect named RNG streams or layout"
    );

    replace_scenario(
        &long_root,
        &format!(
            "max_working_memory_bytes = {}\n",
            long.working_memory_bytes()
        ),
        &format!(
            "max_working_memory_bytes = {}\n",
            long.working_memory_bytes() - 1
        ),
    );
    assert_eq!(
        load_and_validate_f5_config(&long_general)
            .expect_err("admission owns the profile and scan memory budget")
            .code(),
        "GPS-RES-101"
    );

    fs::remove_dir_all(short_root).expect("short cleanup");
    fs::remove_dir_all(long_root).expect("long cleanup");
}

#[test]
fn cli_enforces_the_physical_toml_limit_before_materializing_a_comment() {
    const TOML_LIMIT: usize = 65_536;
    let root = temporary("toml-input-bound");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AA");
    let original = fs::read(&general).expect("general TOML");
    let mut exact = original.clone();
    exact.extend(std::iter::repeat_n(b'#', TOML_LIMIT - exact.len()));
    fs::write(&general, &exact).expect("exact TOML");

    let exact_output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["validate", general.to_str().expect("UTF-8 path")])
        .output()
        .expect("exact CLI validation");
    assert!(
        exact_output.status.success(),
        "exact physical TOML limit failed: {}",
        String::from_utf8_lossy(&exact_output.stderr)
    );
    let exact_config = load_and_validate_f5_config(&general).expect("exact configuration");
    let exact_preflight = preflight_f5_resources(&exact_config, F5OtherResources::default())
        .expect("exact preflight");

    let mut oversized = exact;
    oversized.push(b'#');
    fs::write(&general, oversized).expect("oversized TOML");
    let rejected = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["validate", general.to_str().expect("UTF-8 path")])
        .output()
        .expect("oversized CLI validation");
    assert_eq!(rejected.status.code(), Some(2));
    let rejected_stderr = String::from_utf8_lossy(&rejected.stderr);
    assert!(rejected_stderr.contains("GPS-INPUT-001"));
    assert!(!rejected_stderr.contains("scenario_id="));
    assert!(!rejected_stderr.contains("F5"));

    fs::write(&general, original).expect("restore minimal TOML");
    let minimal_config = load_and_validate_f5_config(&general).expect("minimal configuration");
    let minimal_preflight = preflight_f5_resources(&minimal_config, F5OtherResources::default())
        .expect("minimal preflight");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    assert_eq!(
        build_f5_physical_layout(
            F5LayoutRequest::new(&exact_config, &exact_preflight, &execution)
                .expect("exact request"),
        )
        .expect("exact layout"),
        build_f5_physical_layout(
            F5LayoutRequest::new(&minimal_config, &minimal_preflight, &execution)
                .expect("minimal request"),
        )
        .expect("minimal layout"),
        "accepted TOML comments do not alter F5 RNG or layout"
    );
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn cli_validate_executes_point_four_without_creating_a_destination() {
    let cases = [
        "declared",
        "missing-fasta",
        "malformed-fasta",
        "invalid-mask",
    ];
    for case in cases {
        let root = temporary(case);
        fs::create_dir_all(&root).expect("root");
        let mask = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\n";
        let general = configuration(&root, (case == "invalid-mask").then_some(mask), "", "AAAA");
        match case {
            "missing-fasta" => {
                fs::remove_file(root.join("authoring/references/ref.fa")).expect("remove FASTA");
            }
            "malformed-fasta" => replace_reference(&root, b">chr1\nAXAA\n"),
            "invalid-mask" => replace_mask(&root, b"invalid\n"),
            "declared" => {}
            _ => unreachable!(),
        }
        // El destino se pide de verdad: comprobar que no existe uno que nunca se
        // solicitó no acredita nada.
        let destination = root.join("must-not-exist");
        let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
            .args([
                "validate",
                general.to_str().expect("UTF-8 path"),
                "--output",
                destination.to_str().expect("UTF-8 path"),
            ])
            .output()
            .expect("CLI validate");
        assert_eq!(output.status.success(), case == "declared", "{case}");
        assert!(
            !destination.exists(),
            "validate never creates a destination"
        );
        if case != "declared" {
            assert!(
                !output.stderr.is_empty(),
                "{case} must report its rejection"
            );
        }
        fs::remove_dir_all(root).expect("cleanup");
    }
}

/// F5.1.2.5: la mitad semántica del binding no distingue dos TOML iguales que
/// difieran en un comentario o en su relleno; la mitad física sí, porque el
/// receipt sella el SHA observado de cada documento.
#[test]
fn sealed_receipt_rejects_semantically_equal_toml_with_different_comments() {
    let root = temporary("sealed-toml-crossing");
    fs::create_dir_all(&root).expect("root");
    let general_path = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general_path).expect("configuration");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("sealed admission");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("own binding"),
    )
    .expect("the sealed admission opens layout for its own bytes");

    for (label, path) in [
        ("scenario", root.join("authoring/scenarios/scenario.toml")),
        ("general", general_path.clone()),
    ] {
        let original = fs::read_to_string(&path).expect("document source");
        let padded = format!("{original}# {label} padding sin efecto semántico\n");
        fs::write(&path, &padded).expect("padded document");
        let reloaded = load_and_validate_f5_config(&general_path).expect("padded configuration");
        // La identidad semántica no cambia: el plan de longitudes sigue ligado.
        F5BlockLengthPlan::from_preflight(&reloaded, &sealed)
            .expect("a comment does not change the semantic identity");
        let error = F5LayoutRequest::new(&reloaded, &sealed, &execution)
            .expect_err("the sealed receipt binds the exact bytes it admitted");
        assert_eq!(error.code(), "GPS-CONF-101");
        assert_eq!(error.context()["fields"], "preflight.binding");
        // El cruce recíproco tampoco se admite.
        let reloaded_seal = recover_sealed_f5_admission(&reloaded, F5OtherResources::default())
            .expect("padded admission");
        let error = F5LayoutRequest::new(&config, &reloaded_seal, &execution)
            .expect_err("the reciprocal crossing is equally inadmissible");
        assert_eq!(error.code(), "GPS-CONF-101");
        fs::write(&path, original).expect("restored document");
    }
    fs::remove_dir_all(root).expect("cleanup");
}

/// El metadata declara tamaño y SHA de cada TSV, de modo que intercambiar dos
/// tablas se rechaza al cargarlas; y una admisión sellada sobre el paquete
/// original nunca puede abrir layout sobre el intercambiado.
#[test]
fn swapped_profile_tables_are_rejected_on_load_and_never_cross_a_sealed_receipt() {
    let root = temporary("swapped-profile-tables");
    fs::create_dir_all(&root).expect("root");
    let general_path = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general_path).expect("configuration");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("sealed admission");
    let identities = sealed.receipt().input_identities();
    let gaps_digest = identities[4].observed_sha256().to_owned();
    let spectrum_digest = identities[5].observed_sha256().to_owned();
    assert_ne!(gaps_digest, spectrum_digest);

    let profiles = root.join("authoring/scenarios/profiles");
    let gaps = profiles.join("snp_gaps.tsv");
    let spectrum = profiles.join("snp_spectrum.tsv");
    let gaps_source = fs::read_to_string(&gaps).expect("gap table");
    let spectrum_source = fs::read_to_string(&spectrum).expect("spectrum table");
    fs::write(&gaps, &spectrum_source).expect("swapped gap table");
    fs::write(&spectrum, &gaps_source).expect("swapped spectrum table");

    let error = load_and_validate_f5_config(&general_path)
        .expect_err("a swapped table contradicts the identity its metadata declares");
    assert_eq!(error.code(), "GPS-INPUT-002");
    assert!(
        error.context()["path"].ends_with("snp_gaps.tsv"),
        "the diagnostic must name the swapped table: {:?}",
        error.context()["path"]
    );

    // El paquete intercambiado se reconstruye con su propio metadata: sus
    // identidades físicas son legítimas, pero no las que el receipt selló.
    let swapped_general = {
        let swapped_root = root.join("swapped");
        let path = custom_configuration_with_gap(
            &swapped_root,
            "",
            "AAAA",
            "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
            "gap_min_bp\tgap_max_bp\tweight\n1\t2\t1\n",
            "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
        );
        load_and_validate_f5_config(&path).expect("independent profile package")
    };
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let error = F5LayoutRequest::new(&swapped_general, &sealed, &execution)
        .expect_err("a sealed receipt does not admit another profile package");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(error.context()["fields"], "preflight.binding");
    fs::remove_dir_all(root).expect("cleanup");
}

/// La recuperación REF exige tamaño y SHA idénticos a los sellados y nunca lee
/// más de `tamaño + 1` bytes, de modo que un crecimiento se prueba con un solo
/// byte en lugar de recorrerse.
#[test]
fn second_fasta_pass_rejects_modified_truncated_and_grown_references() {
    let sealed_fasta = b">chr1\nAAAA\n";
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let cases: [(&str, Vec<u8>, u64); 3] = [
        (
            "modified",
            b">chr1\nAAAT\n".to_vec(),
            sealed_fasta.len() as u64,
        ),
        (
            "truncated",
            b">chr1\nAAA\n".to_vec(),
            (sealed_fasta.len() - 1) as u64,
        ),
        (
            "grown",
            format!(">chr1\n{}\n", "A".repeat(20_000)).into_bytes(),
            sealed_fasta.len() as u64 + 1,
        ),
    ];
    for (label, replacement, expected_observed_size) in cases {
        let root = temporary(&format!("second-pass-{label}"));
        fs::create_dir_all(&root).expect("root");
        let general_path = configuration(&root, None, "", "AAAA");
        let config = load_and_validate_f5_config(&general_path).expect("configuration");
        let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
            .expect("sealed admission");
        let reference = sealed
            .receipt()
            .identity(genome_population_simulator::f5_admission::InputRole::ReferenceFasta)
            .expect("sealed reference identity");
        assert_eq!(reference.observed_size_bytes(), sealed_fasta.len() as u64);
        assert_eq!(reference.observed_sha256(), sha256(sealed_fasta));

        fs::write(root.join("authoring/references/ref.fa"), &replacement)
            .expect("reference mutated between passes");
        let error = build_f5_physical_layout(
            F5LayoutRequest::new(&config, &sealed, &execution).expect("sealed binding"),
        )
        .expect_err("the second pass must not accept REF from another reference");
        assert_eq!(error.code(), "GPS-INPUT-002", "case {label}");
        assert_eq!(
            error.context()["expected_size_bytes"],
            sealed_fasta.len().to_string(),
            "case {label}"
        );
        // El crecimiento se acredita con el byte sobrante, no recorriendo los
        // 20.007 bytes del fichero sustituido.
        assert_eq!(
            error.context()["observed_size_bytes"],
            expected_observed_size.to_string(),
            "case {label}"
        );
        assert_eq!(
            error.context()["expected_digest"],
            sha256(sealed_fasta),
            "case {label}"
        );
        fs::remove_dir_all(root).expect("cleanup");
    }
}

/// Una referencia intacta entre pasadas sigue admitiéndose, de modo que la
/// cota `tamaño + 1` no rechaza el recorrido legítimo por su último byte.
#[test]
fn second_fasta_pass_accepts_the_exact_sealed_size() {
    let root = temporary("second-pass-exact");
    fs::create_dir_all(&root).expect("root");
    let general_path = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general_path).expect("configuration");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("sealed admission");
    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("sealed binding"),
    )
    .expect("an untouched reference recovers REF under the sealed bound");
    assert_eq!(layout.blocks().len(), 1);
    fs::remove_dir_all(root).expect("cleanup");
}

/// Las fases selladas siguen el orden cerrado de la máquina de admisión, son
/// continuas en retención y cierran exactamente los totales publicados.
#[test]
fn sealed_phases_follow_the_closed_admission_order_and_close_the_totals() {
    let root = temporary("sealed-phases");
    fs::create_dir_all(&root).expect("root");
    let mask = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\n";
    let general = configuration(&root, Some(mask), "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("sealed admission");
    let phases = sealed.receipt().phases();
    // Oráculo literal del orden cerrado, independiente del enum productivo.
    assert_eq!(
        phases
            .iter()
            .map(|record| record.phase().as_str())
            .collect::<Vec<_>>(),
        vec![
            "bootstrap",
            "profile_metadata",
            "profile_table_block_lengths",
            "profile_table_snp_gaps",
            "profile_table_snp_spectrum",
            "profile_table_haplotype_weights",
            "profile_table_ld",
            "reference_and_mask_scan",
            "final_formula",
        ]
    );
    assert_eq!(phases[0].retained_before_bytes(), 0);
    for pair in phases.windows(2) {
        assert!(
            pair[0].operation_units_cumulative() <= pair[1].operation_units_cumulative(),
            "cumulative operations must not decrease"
        );
        assert!(pair[0].output_bytes_reserved() <= pair[1].output_bytes_reserved());
        assert_eq!(
            pair[0].retained_after_bytes(),
            pair[1].retained_before_bytes(),
            "retention must be continuous across phases"
        );
    }
    for record in phases {
        assert!(
            record.phase_peak_bytes()
                >= record
                    .retained_before_bytes()
                    .max(record.retained_after_bytes()),
            "a phase peak never falls below its live retention"
        );
        assert!(record.phase_peak_bytes() <= sealed.working_memory_bytes());
    }
    let last = phases.last().expect("final formula record");
    assert_eq!(last.operation_units_cumulative(), sealed.operation_units());
    assert_eq!(last.output_bytes_reserved(), sealed.output_bytes());
    assert_eq!(
        phases
            .iter()
            .map(genome_population_simulator::f5_admission::AdmissionPhaseRecord::phase_peak_bytes)
            .max(),
        Some(sealed.working_memory_bytes())
    );
    fs::remove_dir_all(root).expect("cleanup");
}

/// El receipt sellado publica el algoritmo v8 y su proyección reproduce, una a
/// una, las identidades y fases que la admisión observó.
#[test]
fn sealed_projection_publishes_identities_and_phases_in_receipt_order() {
    let root = temporary("sealed-projection");
    fs::create_dir_all(&root).expect("root");
    let mask = "region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t2\t2\n";
    let general = configuration(&root, Some(mask), "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("sealed admission");
    assert_eq!(sealed.receipt().algorithm(), "f5-resource-accounting-v8");

    let projection = sealed.to_contract_json();
    assert_f5_1_2_5_v8_contract_key_set(&projection);
    let identities = projection["input_identities"]
        .as_array()
        .expect("identity list");
    assert_eq!(identities.len(), 10);
    for (published, sealed_identity) in identities.iter().zip(sealed.receipt().input_identities()) {
        assert_eq!(published["role"], sealed_identity.role().as_str());
        assert_eq!(
            published["normalized_path"],
            sealed_identity.normalized_path().display().to_string()
        );
        assert_eq!(
            published["observed_size_bytes"],
            sealed_identity.observed_size_bytes()
        );
        assert_eq!(
            published["observed_sha256"],
            sealed_identity.observed_sha256()
        );
        match sealed_identity.declared_size_bytes() {
            Some(size) => assert_eq!(published["declared_size_bytes"], size),
            None => assert!(published["declared_size_bytes"].is_null()),
        }
        match sealed_identity.declared_sha256() {
            Some(digest) => assert_eq!(published["declared_sha256"], digest),
            None => assert!(published["declared_sha256"].is_null()),
        }
    }
    let published_phases = projection["phases"].as_array().expect("phase list");
    assert_eq!(published_phases.len(), 9);
    for (published, record) in published_phases.iter().zip(sealed.receipt().phases()) {
        assert_eq!(published["phase"], record.phase().as_str());
        assert_eq!(
            published["operation_units_cumulative"],
            record.operation_units_cumulative()
        );
        assert_eq!(
            published["retained_before_bytes"],
            record.retained_before_bytes()
        );
        assert_eq!(published["transient_bytes"], record.transient_bytes());
        assert_eq!(
            published["retained_after_bytes"],
            record.retained_after_bytes()
        );
        assert_eq!(published["phase_peak_bytes"], record.phase_peak_bytes());
        assert_eq!(
            published["output_bytes_reserved"],
            record.output_bytes_reserved()
        );
    }
    assert_preflight_projection_validates_through_effective_and_manifest(&root, &projection);
    fs::remove_dir_all(root).expect("cleanup");
}

// ---------------------------------------------------------------------------
// F5.1.2.6 — segundo oráculo de `resource_preflight` v8.
//
// Reconstruye el objeto contractual completo —algoritmo, identidades físicas,
// las nueve fases y los veintisiete escalares— desde los bytes de las entradas
// y los hechos literales del fixture. No usa `SealedF5Admission`,
// `to_contract_json`, el ledger ni ninguna otra función productiva: cada
// término se despliega aquí con su justificación.
//
// `literal_fixture_operation_oracle` cubre una sola magnitud escalar; este
// segundo oráculo cubre el objeto entero, incluida la contabilidad por fase
// que F5.1.2.5 publicó.
// ---------------------------------------------------------------------------

/// Hechos literales del fixture `configuration(root, mask, limits, "AAA")`.
/// Ninguno procede de la admisión: todos se deducen de los documentos que el
/// propio fixture escribe.
struct SealedProjectionOracle {
    /// `>chr1\nAAA\n`: diez bytes físicos y tres bases.
    fasta_bytes: u64,
    reference_bases: u64,
    /// Un único ID FASTA retenido para comprobar unicidad.
    retained_sequence_ids: u64,
    /// Con o sin máscara el fixture deja un solo tramo elegible: `1..=3` sin
    /// máscara y `1..=2` cuando la fila `m1 chr1 3 3` excluye la última base.
    tract_count: u64,
    /// `block_count = 1`, una subpoblación, un fundador y ploidía 2.
    block_count: u64,
    population_count: u64,
    founder_count: u64,
    ploidy: u64,
    /// Perfil: longitudes `[2,2]`, gaps `[1,1]`, dos rangos haplotípicos, una
    /// fila de espectro y una de LD.
    haplotype_count: u64,
    gap_support_cardinality: u64,
    /// Un único bin de LD declarado.
    bin_count_total: u64,
    /// Filas retenidas y filas normalizadas del paquete de perfiles.
    profile_rows: u64,
    normalization_rows: u64,
    /// Cotas de filas TSV y de elementos JSON del fixture.
    tsv_rows_bound: u64,
    json_rows_bound: u64,
    /// Dos pasadas literales por cada uno de los dos sorteos de gap.
    gap_selector_operation_units: u64,
    /// Límites declarados o duros vigentes.
    operation_limit: u64,
    working_memory_limit_bytes: u64,
    output_limit_bytes: u64,
}

impl SealedProjectionOracle {
    fn literal_fixture(operation_limit: u64) -> Self {
        Self {
            fasta_bytes: 10,
            reference_bases: 3,
            retained_sequence_ids: 1,
            tract_count: 1,
            block_count: 1,
            population_count: 1,
            founder_count: 1,
            ploidy: 2,
            haplotype_count: 2,
            gap_support_cardinality: 1,
            bin_count_total: 1,
            profile_rows: 6,
            normalization_rows: 5,
            tsv_rows_bound: 11,
            json_rows_bound: 30,
            gap_selector_operation_units: 4,
            operation_limit,
            working_memory_limit_bytes: 1 << 30,
            output_limit_bytes: 1 << 31,
        }
    }

    /// `S_bound = 1 + (L_max - 1) / gap_min` con `L_max = 2` y `gap_min = 1`.
    const fn site_bound(&self) -> u64 {
        2
    }

    /// `C_total = founder_count * ploidy`.
    const fn founder_copies(&self) -> u64 {
        self.founder_count * self.ploidy
    }

    /// `K_bound = B * S_bound`.
    const fn site_bound_total(&self) -> u64 {
        self.block_count * self.site_bound()
    }

    /// `X = C_total * sequence_count`.
    const fn copy_sequences(&self) -> u64 {
        self.founder_copies() * self.retained_sequence_ids
    }

    /// `G_bound = founder_count * K_bound`.
    const fn founder_site_bound(&self) -> u64 {
        self.founder_count * self.site_bound_total()
    }

    /// Estructuras tipadas derivables del escenario: matriz haplotípica, bins
    /// de LD, copias fundadoras, asientos de asignación y el registro de
    /// identidades que detecta toda colisión.
    ///
    /// `copies` es `X` en la fórmula final y `C_total` en la reserva anterior
    /// al escaneo, donde `sequence_count` todavía no se ha observado y
    /// `C_total` es su cota inferior.
    fn derived_structures(&self, copies: u64) -> u64 {
        let site_bound = self.site_bound();
        let pairs = site_bound * (site_bound - 1) / 2;
        let packed = (self.haplotype_count * site_bound).div_ceil(8);
        let seats = self.block_count * self.founder_copies();
        let registry =
            copies + self.block_count * self.population_count * self.haplotype_count + seats;
        self.block_count
            * self.population_count
            * (packed + (160 + site_bound) * self.haplotype_count + 256 * site_bound + 192 * pairs)
            + 256 * self.block_count * self.bin_count_total
            + 192 * copies
            + 64 * seats
            + 128 * registry
    }

    /// Transitorio conocido de la fórmula final, reservado antes de abrir el
    /// FASTA: no incluye todavía tramos ni IDs.
    fn known_final_transient(&self) -> u64 {
        1_024 * self.block_count + self.derived_structures(self.founder_copies()) + 65_536
    }

    /// Pico de la fórmula final: vectores de bloques, tramos, IDs y filas del
    /// perfil junto a las estructuras derivadas y el búfer fijo del escáner.
    fn final_formula_peak(&self) -> u64 {
        256 * (4 * self.block_count
            + 2 * self.tract_count
            + self.retained_sequence_ids
            + self.profile_rows
            + self.normalization_rows)
            + self.derived_structures(self.copy_sequences())
            + 65_536
    }

    /// Retención tipada del perfil que sobrevive al paquete de tablas.
    fn profile_retained_bytes(&self) -> u64 {
        256 * (self.profile_rows + self.normalization_rows)
    }

    /// Retención al cerrar el escaneo: perfil más tramos e IDs.
    fn scan_retained_bytes(&self) -> u64 {
        self.profile_retained_bytes() + 256 * (self.tract_count + self.retained_sequence_ids)
    }

    /// Envolvente de salida ya conocida antes de cerrar la fórmula final.
    fn known_output_bytes(&self) -> u64 {
        1_048_576
            + 4_096 * (self.tsv_rows_bound + self.json_rows_bound)
            + 65_536
            + 4_096 * self.founder_site_bound()
    }

    /// Envolvente final de salida, con la cota de verdad completa.
    fn output_bytes(&self) -> u64 {
        let truth = 65_536
            + self.copy_sequences() * (4_096 + 2_048 * self.site_bound_total())
            + 4_096 * self.founder_site_bound();
        1_048_576 + 4_096 * (self.tsv_rows_bound + self.json_rows_bound) + truth
    }

    /// Reconstruye el objeto `resource_preflight` v8 completo.
    fn project(&self, root: &Path, mask: Option<&[u8]>) -> serde_json::Value {
        let authoring = root.join("authoring");
        let scenario_directory = authoring.join("scenarios");
        let profile_directory = scenario_directory.join("profiles");
        let read = |path: &Path| fs::read(path).expect("byte del fixture");

        let general_path = authoring.join("general.toml");
        let scenario_path = scenario_directory.join("scenario.toml");
        let metadata_path = profile_directory.join("profile.metadata.json");
        let fasta_path = authoring.join("references/ref.fa");
        let mask_path = scenario_directory.join("excluded_regions.tsv");

        let general = read(&general_path);
        let scenario = read(&scenario_path);
        let metadata = read(&metadata_path);
        let fasta = read(&fasta_path);
        assert_eq!(fasta.len() as u64, self.fasta_bytes);

        // Orden contractual de las cinco tablas y filas que cada una retiene:
        // sus filas propias más sus pesos normalizados. El LD no normaliza.
        let tables = [
            (
                "block_lengths.tsv",
                "profile_block_lengths",
                "profile_table_block_lengths",
                1 + 1,
            ),
            (
                "snp_gaps.tsv",
                "profile_snp_gaps",
                "profile_table_snp_gaps",
                1 + 1,
            ),
            (
                "snp_spectrum.tsv",
                "profile_snp_spectrum",
                "profile_table_snp_spectrum",
                1 + 1,
            ),
            (
                "haplotype_weights.tsv",
                "profile_haplotype_weights",
                "profile_table_haplotype_weights",
                2 + 2,
            ),
            // El LD no normaliza pesos: retiene solo su fila.
            ("ld_profile.tsv", "profile_ld", "profile_table_ld", 1),
        ];
        let table_bytes = tables
            .iter()
            .map(|(name, _, _, _)| read(&profile_directory.join(name)))
            .collect::<Vec<_>>();
        let profile_tsv_bytes = table_bytes
            .iter()
            .map(|bytes| bytes.len() as u64)
            .sum::<u64>();
        let largest_table_bytes = table_bytes
            .iter()
            .map(|bytes| bytes.len() as u64)
            .max()
            .expect("cinco tablas");
        assert_eq!(
            tables.iter().map(|(_, _, _, rows)| *rows).sum::<u64>(),
            self.profile_rows + self.normalization_rows,
            "las filas retenidas por tabla suman las filas del perfil"
        );

        let general_bytes = general.len() as u64;
        let scenario_bytes = scenario.len() as u64;
        let metadata_bytes = metadata.len() as u64;
        let mask_bytes = mask.map_or(0, |bytes| bytes.len() as u64);
        if let Some(expected) = mask {
            assert_eq!(read(&mask_path), expected, "la máscara literal del fixture");
        } else {
            assert!(!mask_path.exists());
        }

        // Identidades físicas en orden de rol contractual. Solo el metadata,
        // los cinco TSV, la referencia y la máscara tienen SHA declarado; solo
        // los TSV declaran además su tamaño.
        let mut identities = vec![
            serde_json::json!({
                "role": "general_toml",
                "normalized_path": general_path.display().to_string(),
                "observed_size_bytes": general_bytes,
                "observed_sha256": sha256(&general),
                "declared_size_bytes": serde_json::Value::Null,
                "declared_sha256": serde_json::Value::Null,
            }),
            serde_json::json!({
                "role": "scenario_toml",
                "normalized_path": scenario_path.display().to_string(),
                "observed_size_bytes": scenario_bytes,
                "observed_sha256": sha256(&scenario),
                "declared_size_bytes": serde_json::Value::Null,
                "declared_sha256": serde_json::Value::Null,
            }),
            serde_json::json!({
                "role": "profile_metadata",
                "normalized_path": metadata_path.display().to_string(),
                "observed_size_bytes": metadata_bytes,
                "observed_sha256": sha256(&metadata),
                "declared_size_bytes": serde_json::Value::Null,
                "declared_sha256": sha256(&metadata),
            }),
        ];
        for ((name, role, _, _), bytes) in tables.iter().zip(&table_bytes) {
            identities.push(serde_json::json!({
                "role": role,
                "normalized_path": profile_directory.join(name).display().to_string(),
                "observed_size_bytes": bytes.len() as u64,
                "observed_sha256": sha256(bytes),
                "declared_size_bytes": bytes.len() as u64,
                "declared_sha256": sha256(bytes),
            }));
        }
        identities.push(serde_json::json!({
            "role": "reference_fasta",
            "normalized_path": fasta_path.display().to_string(),
            "observed_size_bytes": self.fasta_bytes,
            "observed_sha256": sha256(&fasta),
            "declared_size_bytes": serde_json::Value::Null,
            "declared_sha256": sha256(&fasta),
        }));
        if let Some(bytes) = mask {
            identities.push(serde_json::json!({
                "role": "mask",
                "normalized_path": mask_path.display().to_string(),
                "observed_size_bytes": mask_bytes,
                "observed_sha256": sha256(bytes),
                "declared_size_bytes": serde_json::Value::Null,
                "declared_sha256": sha256(bytes),
            }));
        }

        let operation_units = literal_fixture_operation_oracle(root, mask);

        // Picos por fase de la máquina de admisión, cada uno con las únicas
        // estructuras vivas de su fase.
        let bootstrap_peak = 2 * (general_bytes + scenario_bytes);
        let metadata_peak = 2 * metadata_bytes;
        let table_peak = 2 * largest_table_bytes + self.profile_retained_bytes();
        let scan_transient = 65_536 + if mask.is_some() { 65_536 } else { 0 };
        let scan_peak = scan_transient + self.scan_retained_bytes();
        let working_memory_bytes = [
            bootstrap_peak,
            metadata_peak,
            table_peak,
            scan_peak,
            self.final_formula_peak(),
        ]
        .into_iter()
        .max()
        .expect("cinco fases de memoria");

        let mut phases = vec![
            serde_json::json!({
                "phase": "bootstrap",
                "operation_units_cumulative": bootstrap_peak,
                "retained_before_bytes": 0,
                "transient_bytes": bootstrap_peak,
                "retained_after_bytes": 0,
                "phase_peak_bytes": bootstrap_peak,
                "output_bytes_reserved": 0,
            }),
            serde_json::json!({
                "phase": "profile_metadata",
                "operation_units_cumulative":
                    bootstrap_peak + 2 * metadata_bytes + 2 * profile_tsv_bytes,
                "retained_before_bytes": 0,
                // Antes de abrir la primera tabla se reserva el mayor búfer
                // de las cinco: sus lecturas no coexisten.
                "transient_bytes": 2 * largest_table_bytes,
                "retained_after_bytes": 0,
                "phase_peak_bytes": metadata_peak.max(2 * largest_table_bytes),
                "output_bytes_reserved": 0,
            }),
        ];
        let profile_cumulative = bootstrap_peak + 2 * metadata_bytes + 2 * profile_tsv_bytes;
        let mut retained_rows = 0_u64;
        for ((_, _, phase, rows), bytes) in tables.iter().zip(&table_bytes) {
            let before = 256 * retained_rows;
            retained_rows += rows;
            let after = 256 * retained_rows;
            let transient = 2 * bytes.len() as u64;
            phases.push(serde_json::json!({
                "phase": phase,
                "operation_units_cumulative": profile_cumulative,
                "retained_before_bytes": before,
                "transient_bytes": transient,
                "retained_after_bytes": after,
                "phase_peak_bytes": after + transient,
                "output_bytes_reserved": 0,
            }));
        }
        let known_final_transient = self.known_final_transient();
        phases.push(serde_json::json!({
            "phase": "reference_and_mask_scan",
            "operation_units_cumulative": operation_units,
            "retained_before_bytes": self.profile_retained_bytes(),
            "transient_bytes": known_final_transient.max(scan_transient),
            "retained_after_bytes": self.scan_retained_bytes(),
            // La reserva prospectiva de la fórmula final coexiste con la
            // retención del perfil; el búfer del escáner, con tramos e IDs.
            "phase_peak_bytes": (self.profile_retained_bytes() + known_final_transient)
                .max(self.scan_retained_bytes() + scan_transient),
            "output_bytes_reserved": self.known_output_bytes(),
        }));
        phases.push(serde_json::json!({
            "phase": "final_formula",
            "operation_units_cumulative": operation_units,
            "retained_before_bytes": self.scan_retained_bytes(),
            "transient_bytes": 0,
            "retained_after_bytes": self.scan_retained_bytes(),
            "phase_peak_bytes": working_memory_bytes.max(self.scan_retained_bytes()),
            "output_bytes_reserved": self.output_bytes(),
        }));

        serde_json::json!({
            "algorithm": "f5-resource-accounting-v8",
            "input_identities": identities,
            "phases": phases,
            "reference_bases_processed": self.reference_bases,
            "fasta_bytes_processed": self.fasta_bytes,
            "mask_bytes_processed": mask_bytes,
            "general_toml_bytes": general_bytes,
            "scenario_toml_bytes": scenario_bytes,
            "profile_metadata_bytes": metadata_bytes,
            "profile_tsv_bytes": profile_tsv_bytes,
            "tract_count": self.tract_count,
            "retained_sequence_ids": self.retained_sequence_ids,
            "gap_support_cardinality": self.gap_support_cardinality,
            "gap_selector_operation_units": self.gap_selector_operation_units,
            "site_bound_per_block": self.site_bound(),
            "bin_count_total": self.bin_count_total,
            "profile_rows": self.profile_rows,
            "normalization_rows": self.normalization_rows,
            "tsv_rows_bound": self.tsv_rows_bound,
            "json_rows_bound": self.json_rows_bound,
            "other_operation_units": 0,
            "other_working_memory_bytes": 0,
            "other_output_bytes": 0,
            "operation_units": operation_units,
            "working_memory_bytes": working_memory_bytes,
            "output_bytes": self.output_bytes(),
            "operation_limit": self.operation_limit,
            "working_memory_limit_bytes": self.working_memory_limit_bytes,
            "output_limit_bytes": self.output_limit_bytes,
        })
    }
}

/// F5.1.2.6: la proyección v8 completa se reconstruye desde las entradas y los
/// hechos literales del fixture, sin tipos ni funciones productivas.
#[test]
fn an_independent_oracle_reconstructs_the_whole_v8_projection_from_inputs() {
    let mask = b"region_id\tsequence_id\tstart_1based\tend_1based\nm1\tchr1\t3\t3\n";
    for (label, mask) in [
        ("oracle-without-mask", None),
        ("oracle-with-mask", Some(mask.as_slice())),
    ] {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("root");
        let general = configuration(
            &root,
            mask.map(|bytes| std::str::from_utf8(bytes).expect("UTF-8 mask")),
            "max_operation_units = 999999\n",
            "AAA",
        );
        let config = load_and_validate_f5_config(&general).expect("configuración admitida");
        let sealed =
            preflight_f5_resources(&config, F5OtherResources::default()).expect("admisión sellada");
        let projection = sealed.to_contract_json();
        let oracle = SealedProjectionOracle::literal_fixture(999_999);
        assert_eq!(
            projection,
            oracle.project(&root, mask),
            "el oráculo independiente debe reconstruir la proyección v8: {label}"
        );
        // La reconstrucción también debe superar los dos esquemas publicados.
        assert_preflight_projection_validates_through_effective_and_manifest(&root, &projection);
        fs::remove_dir_all(root).expect("cleanup");
    }
}

/// F5.1.2.6: el comparador productivo exige que las dos copias publicadas de
/// `resource_preflight` sean la proyección sellada. Los esquemas comparten el
/// mismo `$def` y no pueden observar una desigualdad dirigida; esta puerta sí.
#[test]
fn effective_and_manifest_are_rejected_when_one_published_copy_diverges() {
    let root = temporary("contract-copies");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAA");
    let config = load_and_validate_f5_config(&general).expect("configuración admitida");
    let sealed =
        preflight_f5_resources(&config, F5OtherResources::default()).expect("admisión sellada");
    let projection = sealed.to_contract_json();

    let effective_path = root.join("effective-config.json");
    let manifest_path = root.join("manifest.json");
    let document = |projection: &serde_json::Value| serde_json::json!({ "f5": { "resource_preflight": projection.clone() } });
    let effective = document(&projection);
    let manifest = document(&projection);
    sealed
        .ensure_contract_copies_match(&effective_path, &effective, &manifest_path, &manifest)
        .expect("dos copias fieles de la admisión sellada se aceptan");

    // Desigualdad dirigida en cada clase de término del contrato v8.
    let directed = [
        ("/operation_units", serde_json::json!(1_u64)),
        ("/other_output_bytes", serde_json::json!(1_u64)),
        (
            "/input_identities/8/observed_sha256",
            serde_json::json!("a".repeat(64)),
        ),
        (
            "/input_identities/3/declared_size_bytes",
            serde_json::json!(2_u64),
        ),
        ("/phases/7/phase_peak_bytes", serde_json::json!(2_u64)),
        ("/phases/0/phase", serde_json::json!("profile_metadata")),
    ];
    for (pointer, replacement) in directed {
        for divergent_copy in ["effective", "manifest"] {
            let mut effective = document(&projection);
            let mut manifest = document(&projection);
            let target = if divergent_copy == "effective" {
                &mut effective
            } else {
                &mut manifest
            };
            let slot = target
                .pointer_mut(&format!("/f5/resource_preflight{pointer}"))
                .expect("puntero dirigido");
            assert_ne!(*slot, replacement, "la mutación debe cambiar el valor");
            *slot = replacement.clone();
            let error = sealed
                .ensure_contract_copies_match(
                    &effective_path,
                    &effective,
                    &manifest_path,
                    &manifest,
                )
                .expect_err("una copia discrepante debe rechazarse");
            assert_eq!(error.code(), "GPS-FMT-006");
            assert!(
                error.context()["detail"].contains(pointer),
                "el detalle debe citar {pointer}: {:?}",
                error.context()["detail"]
            );
            let expected_path = if divergent_copy == "effective" {
                &effective_path
            } else {
                &manifest_path
            };
            assert_eq!(
                error.context()["path"],
                expected_path.display().to_string(),
                "el rechazo debe citar la copia divergente"
            );
        }
    }

    // Una divergencia de forma se nombra por su puntero y su tipo, sin volcar
    // el subárbol: ese detalle acaba publicado en el manifiesto de fallo.
    for (pointer, replacement) in [
        ("", serde_json::Value::Null),
        ("/phases", serde_json::json!("no es una lista")),
        ("/input_identities/0", serde_json::json!(7_u64)),
    ] {
        let mut effective = document(&projection);
        let manifest = document(&projection);
        *effective
            .pointer_mut(&format!("/f5/resource_preflight{pointer}"))
            .expect("puntero de forma") = replacement;
        let error = sealed
            .ensure_contract_copies_match(&effective_path, &effective, &manifest_path, &manifest)
            .expect_err("una copia con otra forma debe rechazarse");
        assert_eq!(error.code(), "GPS-FMT-006");
        let detail = error.context()["detail"].clone();
        assert!(
            detail.len() < 200,
            "el detalle no puede volcar la proyección entera: {} bytes",
            detail.len()
        );
        let cited = if pointer.is_empty() {
            "/f5/resource_preflight"
        } else {
            pointer
        };
        assert!(
            detail.contains(cited),
            "el detalle debe citar {cited}: {detail}"
        );
    }

    // Dos copias iguales entre sí pero distintas de la admisión también se
    // rechazan: el comparador no compara copia contra copia.
    let mut both = projection.clone();
    both["operation_units"] = serde_json::json!(
        projection["operation_units"]
            .as_u64()
            .expect("operation_units u64")
            + 1
    );
    let error = sealed
        .ensure_contract_copies_match(
            &effective_path,
            &document(&both),
            &manifest_path,
            &document(&both),
        )
        .expect_err("dos copias iguales pero ajenas a la admisión se rechazan");
    assert_eq!(error.code(), "GPS-FMT-006");

    // Un término no contratado y una copia ausente también rechazan.
    let mut extra = document(&projection);
    extra["f5"]["resource_preflight"]["uncontracted_term"] = serde_json::json!(1_u64);
    let error = sealed
        .ensure_contract_copies_match(&effective_path, &extra, &manifest_path, &manifest)
        .expect_err("un término no contratado se rechaza");
    assert_eq!(error.code(), "GPS-FMT-006");
    assert!(error.context()["detail"].contains("uncontracted_term"));

    let error = sealed
        .ensure_contract_copies_match(
            &effective_path,
            &serde_json::json!({ "f5": {} }),
            &manifest_path,
            &manifest,
        )
        .expect_err("una copia ausente se rechaza");
    assert_eq!(error.code(), "GPS-FMT-006");
    assert!(error.context()["detail"].contains("/f5/resource_preflight"));

    fs::remove_dir_all(root).expect("cleanup");
}

/// F5.1.2.6: dos contigs, N interior y adyacencia recorridos por la ruta
/// productiva completa —carga, admisión, binding, segunda pasada FASTA y
/// layout—, con una verdad independiente leída del propio FASTA.
#[test]
fn two_contigs_with_interior_n_reach_layout_through_the_production_path() {
    let root = temporary("two-contigs");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    // `chr1` queda partido por una N interior y `chr2` es un tramo íntegro
    // adyacente a su cabecera; ninguno de los tres tramos se solapa.
    let fasta = b">chr1\nAANAA\n>chr2\nCCCC\n";
    replace_reference(&root, fasta);
    let config = load_and_validate_f5_config(&general).expect("configuración admitida");
    let sealed =
        preflight_f5_resources(&config, F5OtherResources::default()).expect("admisión sellada");

    // Verdad independiente: los tramos maximales sin N de los dos contigs.
    let tracts = sealed
        .eligible_tracts()
        .iter()
        .map(|tract| {
            (
                tract.sequence_ordinal(),
                tract.start_1based(),
                tract.end_1based(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(tracts, vec![(0, 1, 2), (0, 4, 5), (1, 1, 4)]);
    assert_eq!(sealed.reference_bases_processed(), 9);
    assert_eq!(sealed.retained_sequence_ids(), 2);
    assert_eq!(sealed.tract_count(), 3);
    assert_eq!(sealed.fasta_bytes_processed(), fasta.len() as u64);

    let execution =
        ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("binding sellado"),
    )
    .expect("el layout admitido se completa sobre dos contigs");

    // Cada bloque cae íntegro dentro de un tramo elegible y sus bases de
    // referencia coinciden con el FASTA leído aparte, no con el escáner.
    let sequences = [(0_u64, "chr1", b"AANAA".as_slice()), (1, "chr2", b"CCCC")];
    for block in layout.blocks() {
        let (_, expected_id, bases) = sequences
            .iter()
            .find(|(ordinal, _, _)| *ordinal == block.sequence_ordinal())
            .expect("ordinal conocido");
        assert_eq!(block.sequence_id().as_str(), *expected_id);
        assert!(
            tracts.iter().any(|(ordinal, start, end)| {
                *ordinal == block.sequence_ordinal()
                    && *start <= block.start_1based()
                    && block.end_1based() <= *end
            }),
            "todo bloque vive dentro de un tramo elegible: {block:?}"
        );
        for site in block.sites() {
            let index = usize::try_from(site.position_1based() - 1).expect("posición usize");
            assert_eq!(site.reference_base(), bases[index]);
            assert_ne!(site.reference_base(), b'N');
            assert_ne!(site.reference_base(), site.alternate_base());
        }
    }
    fs::remove_dir_all(root).expect("cleanup");
}
