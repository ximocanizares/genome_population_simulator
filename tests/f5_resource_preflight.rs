use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_f5_config;
use genome_population_simulator::f5_preflight::{
    F5_HARD_OPERATION_LIMIT, F5_HARD_OUTPUT_LIMIT, F5_HARD_WORKING_MEMORY_LIMIT, F5OtherResources,
    preflight_f5_resources,
};
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

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
    assert_eq!(fragmented.reference_bases_processed(), length as u64);
    assert_eq!(compact.eligible_bases(), fragmented.eligible_bases());
    assert_eq!(compact.profile_rows(), 6);
    assert_eq!(compact.normalization_rows(), 5);
    assert_eq!(compact.bin_count_total(), 1);
    assert_eq!(compact.tsv_rows_bound(), 11);
    assert_eq!(compact.json_rows_bound(), 30);
    assert_eq!(compact.operation_units(), 24_651);
    assert_eq!(compact.working_memory_bytes(), 70_465);
    assert_eq!(compact.output_bytes(), 1_306_624);
    assert_eq!(
        fragmented.operation_units() - compact.operation_units(),
        4_095 * (1 + 1)
    );
    assert_eq!(
        fragmented.working_memory_bytes() - compact.working_memory_bytes(),
        256 * 4_095
    );
    assert_eq!(fragmented.output_bytes(), compact.output_bytes());

    for (label, limit, expected_ok) in [
        (
            "operations-exact",
            format!("max_operation_units = {}\n", fragmented.operation_units()),
            true,
        ),
        (
            "operations-minus-one",
            format!(
                "max_operation_units = {}\n",
                fragmented.operation_units() - 1
            ),
            false,
        ),
        (
            "memory-exact",
            format!(
                "max_working_memory_bytes = {}\n",
                fragmented.working_memory_bytes()
            ),
            true,
        ),
        (
            "memory-minus-one",
            format!(
                "max_working_memory_bytes = {}\n",
                fragmented.working_memory_bytes() - 1
            ),
            false,
        ),
    ] {
        let path = configuration(&root.join(label), Some(&mask), &limit, &bases);
        let config = load_and_validate_f5_config(&path).expect("bounded configuration");
        let result = preflight_f5_resources(&config, F5OtherResources::default());
        if expected_ok {
            result.expect("inclusive fragmented bound");
        } else {
            assert_eq!(
                result.expect_err("bound minus one must reject").code(),
                "GPS-RES-101"
            );
        }
    }
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn reference_path_and_exact_sha_are_both_required_before_accounting() {
    let root = temporary("reference-sha");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("configuration");
    fs::write(root.join("authoring/references/ref.fa"), b">chr1\nAAAT\n")
        .expect("reference mutation");
    let error = preflight_f5_resources(&config, F5OtherResources::default())
        .expect_err("same ID with changed bytes must fail");
    assert_eq!(error.code(), "GPS-INPUT-002");
    assert!(error.render_human().contains("expected_digest"));
    assert!(error.render_human().contains("observed_digest"));
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
        let config = load_and_validate_f5_config(&general).expect("configuration");
        assert_eq!(
            preflight_f5_resources(&config, F5OtherResources::default())
                .expect_err("impossible target")
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
    let config = load_and_validate_f5_config(&general).expect("configuration");
    assert_eq!(
        preflight_f5_resources(&config, F5OtherResources::default())
            .expect_err("override outside B")
            .code(),
        "GPS-CONF-101"
    );
    fs::remove_dir_all(root).expect("cleanup");
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
    let config = load_and_validate_f5_config(&general).expect("configuration");
    assert_eq!(
        preflight_f5_resources(&config, F5OtherResources::default())
            .expect_err("oversized row")
            .code(),
        "GPS-FMT-101"
    );
    fs::remove_dir_all(root).expect("cleanup");
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
    let config = load_and_validate_f5_config(&general).expect("configuration");
    assert_eq!(
        preflight_f5_resources(&config, F5OtherResources::default())
            .expect_err("S bound")
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
    let config = load_and_validate_f5_config(&general).expect("configuration");
    assert_eq!(
        preflight_f5_resources(&config, F5OtherResources::default())
            .expect_err("H*S bound")
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
    let config = load_and_validate_f5_config(&general).expect("configuration");
    for other in [
        F5OtherResources::new(u64::MAX, 0, 0),
        F5OtherResources::new(0, u64::MAX, 0),
        F5OtherResources::new(0, 0, u64::MAX),
    ] {
        assert_eq!(
            preflight_f5_resources(&config, other)
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
    let config = load_and_validate_f5_config(&general).expect("S_bound=4096 configuration");
    let error = preflight_f5_resources(&config, F5OtherResources::default())
        .expect_err("the independent global memory bound remains applicable");
    assert_eq!(error.code(), "GPS-RES-101");
    assert!(error.render_human().contains("working_memory_bytes="));
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
        let config = load_and_validate_f5_config(&general).expect("child configuration");
        let error = preflight_f5_resources(&config, F5OtherResources::default())
            .expect_err("declared memory budget must reject");
        assert_eq!(error.code(), "GPS-RES-101");
        assert!(error.render_human().contains("working_memory_bytes="));
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
        let config = load_and_validate_f5_config(&general).expect("child configuration");
        let error = preflight_f5_resources(&config, F5OtherResources::default())
            .expect_err("declared memory budget must reject");
        assert_eq!(error.code(), "GPS-RES-101");
        assert!(error.render_human().contains("working_memory_bytes="));
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

#[cfg(target_os = "linux")]
#[test]
fn fragmented_tracts_are_budgeted_before_vector_growth() {
    const CHILD_ROOT: &str = "GPS_F5_FRAGMENTED_TRACTS_CHILD_ROOT";
    if let Some(root) = std::env::var_os(CHILD_ROOT) {
        let general = PathBuf::from(root).join("authoring/general.toml");
        let config = load_and_validate_f5_config(&general).expect("child configuration");
        let error = preflight_f5_resources(&config, F5OtherResources::default())
            .expect_err("scan-time memory lower bound must reject");
        assert_eq!(error.code(), "GPS-RES-101");
        assert!(error.render_human().contains("working_memory_bytes="));
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
        let config = load_and_validate_f5_config(&general).expect("child configuration");
        let error = preflight_f5_resources(&config, F5OtherResources::default())
            .expect_err("scan-time ID memory lower bound must reject");
        assert_eq!(error.code(), "GPS-RES-101");
        let rendered = error.render_human();
        assert!(rendered.contains("working_memory_bytes=131328/131072"));
        assert!(rendered.contains("tract_count=0"));
        assert!(rendered.contains("retained_sequence_ids=257"));
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
    assert_eq!(result.working_memory_bytes(), 103_233);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn resource_limit_includes_other_outputs_and_rejects_before_any_destination() {
    let root = temporary("limit");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let baseline =
        preflight_f5_resources(&config, F5OtherResources::new(7, 11, 13)).expect("baseline");
    assert_eq!(baseline.other_resources(), F5OtherResources::new(7, 11, 13));
    let rejected_general = configuration(
        &root.join("rejected"),
        None,
        &format!("max_output_bytes = {}\n", baseline.output_bytes() - 1),
        "AAAA",
    );
    let rejected = load_and_validate_f5_config(&rejected_general).expect("limited configuration");
    let error = preflight_f5_resources(&rejected, F5OtherResources::new(7, 11, 13))
        .expect_err("output budget must reject");
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
    let admitted = preflight_f5_resources(&config, exact).expect("inclusive bounds");
    assert_eq!(admitted.operation_units(), F5_HARD_OPERATION_LIMIT);
    assert_eq!(
        admitted.working_memory_bytes(),
        F5_HARD_WORKING_MEMORY_LIMIT
    );
    assert_eq!(admitted.output_bytes(), F5_HARD_OUTPUT_LIMIT);
    let error = preflight_f5_resources(
        &config,
        F5OtherResources::new(
            exact.operation_units(),
            exact.working_memory_bytes(),
            exact.output_bytes() + 1,
        ),
    )
    .expect_err("one byte beyond the hard limit must fail");
    assert_eq!(error.code(), "GPS-RES-101");
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
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let error =
        preflight_f5_resources(&config, F5OtherResources::default()).expect_err("T=0 must fail");
    assert_eq!(error.code(), "GPS-DOMAIN-101");
    fs::remove_dir_all(root).expect("cleanup");
}
