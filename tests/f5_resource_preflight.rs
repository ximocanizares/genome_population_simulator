use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{
    load_and_validate_f5_config, load_and_validate_f5_config_with_entropy,
};
use genome_population_simulator::f5_layout::{
    F5BlockLengthPlan, F5LayoutRequest, build_f5_physical_layout,
};
use genome_population_simulator::f5_preflight::{
    F5_HARD_OPERATION_LIMIT, F5_HARD_OUTPUT_LIMIT, F5_HARD_WORKING_MEMORY_LIMIT, F5OtherResources,
    preflight_f5_resources,
};
use genome_population_simulator::identity::ExecutionId;
use genome_population_simulator::randomness::EntropySource;
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

#[derive(Debug)]
struct FixedEntropy([u8; 32]);

impl EntropySource for FixedEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        *destination = self.0;
        Ok(())
    }
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
    assert_eq!(compact.operation_units(), 39_305);
    assert_eq!(compact.working_memory_bytes(), 76_037);
    assert_eq!(compact.output_bytes(), 1_306_624);
    assert_eq!(
        fragmented.operation_units() - compact.operation_units(),
        4_095 * (1 + 2 + 2 + 1 + 3) + fragmented.mask_bytes_processed()
            - compact.mask_bytes_processed()
    );
    assert_eq!(
        fragmented.working_memory_bytes() - compact.working_memory_bytes(),
        512 * 4_095
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
        let rejected = load_and_validate_f5_config(&path).expect("minus-one configuration");
        assert_eq!(
            preflight_f5_resources(&rejected, F5OtherResources::default())
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

#[test]
fn fasta_descriptions_raise_v5_operation_budget_and_preserve_productive_rng_path() {
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
        .expect("exact v5 operation limit is inclusive");
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
    let rejected_config =
        load_and_validate_f5_config(&long_general).expect("minus-one configuration");
    assert_eq!(
        preflight_f5_resources(&rejected_config, F5OtherResources::default())
            .expect_err("one operation below the long-header bound must reject")
            .code(),
        "GPS-RES-101"
    );

    fs::remove_dir_all(short_root).expect("short cleanup");
    fs::remove_dir_all(long_root).expect("long cleanup");
}

#[test]
fn mask_bytes_raise_v5_operation_budget_with_equal_geometry_and_exact_limit() {
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
    let rejected_config =
        load_and_validate_f5_config(&long_general).expect("minus-one configuration");
    assert_eq!(
        preflight_f5_resources(&rejected_config, F5OtherResources::default())
            .expect_err("one operation below the mask-byte bound must reject")
            .code(),
        "GPS-RES-101"
    );

    fs::remove_dir_all(short_root).expect("short cleanup");
    fs::remove_dir_all(long_root).expect("long cleanup");
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
    assert_eq!(result.working_memory_bytes(), 108_525);
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
    let config = load_and_validate_f5_config(&general).expect("wide support configuration");
    let error = preflight_f5_resources(&config, F5OtherResources::default())
        .expect_err("the literal two-pass support cannot exceed the operation budget");
    assert_eq!(error.code(), "GPS-RES-101");
    assert!(error.context()["detail"].contains("length_support_cardinality"));
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
    F5LayoutRequest::new(&equivalent, &preflight, &execution)
        .expect("canonical effective identity accepts equivalent defaults");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn layout_request_binds_other_resource_contributions() {
    let root = temporary("binding-other-resources");
    fs::create_dir_all(&root).expect("root");
    let general = configuration(&root, None, "", "AAAA");
    let config = load_and_validate_f5_config(&general).expect("configuration");
    let other = F5OtherResources::new(7, 11, 13);
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
fn cli_admits_or_rejects_profile_document_bytes_before_materializing_them() {
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
        "v6 reserves raw bytes plus the parser-retention copy"
    );
    assert_eq!(
        long.working_memory_bytes() - short.working_memory_bytes(),
        8_184
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
        .expect("the exact v6 memory limit is inclusive");
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
    let rejected = load_and_validate_f5_config(&long_general).expect("minus-one configuration");
    assert_eq!(
        preflight_f5_resources(&rejected, F5OtherResources::default())
            .expect_err("the preflight owns the profile budget before F5.1.2.3")
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
