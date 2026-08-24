//! F2.5.1: contrato físico, sin proyección de variantes ni lectura de verdad.

use std::fs;
use std::io::BufRead;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::vcf_contract::{
    TBI_MAXIMUM_POSITION_EXCLUSIVE, VcfContig, VcfContractError, VcfHeaderContract, VcfIndexKind,
    VcfOutputMode, write_vcf_contract,
};

static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary(label: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!(
        "gps-vcf-contract-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
    ))
}

fn contract(length: u64) -> VcfHeaderContract {
    VcfHeaderContract::new(
        vec![VcfContig::new("chr1", length).expect("contig")],
        vec!["sample-a".to_owned(), "sample-b".to_owned()],
    )
    .expect("contract")
}

#[test]
fn text_mode_has_exact_vcf_45_header_and_no_index() {
    let root = temporary("text");
    fs::create_dir_all(&root).expect("directory");
    let output = write_vcf_contract(&root.join("truth"), VcfOutputMode::Text, &contract(12))
        .expect("text VCF");
    let text = fs::read_to_string(&output.vcf).expect("VCF text");
    assert!(text.starts_with("##fileformat=VCFv4.5\n"));
    assert!(text.contains("##contig=<ID=chr1,length=12>\n"));
    assert!(
        text.ends_with(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-a\tsample-b\n"
        )
    );
    assert_eq!(output.index, None);
    assert_eq!(output.index_kind, None);
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn bgzf_tbi_is_selected_below_normative_limit_and_is_readable() {
    let root = temporary("tbi");
    fs::create_dir_all(&root).expect("directory");
    let output = write_vcf_contract(&root.join("truth"), VcfOutputMode::Bgzf, &contract(12))
        .expect("BGZF VCF");
    assert_eq!(output.index_kind, Some(VcfIndexKind::Tbi));
    let index = output.index.expect("TBI");
    assert_eq!(
        index.extension().and_then(|value| value.to_str()),
        Some("tbi")
    );
    let mut index_reader = noodles_tabix::io::Reader::new(fs::File::open(&index).expect("TBI"));
    index_reader.read_index().expect("readable TBI");
    let mut reader = noodles_bgzf::io::Reader::new(fs::File::open(&output.vcf).expect("BGZF"));
    let mut first_line = String::new();
    reader.read_line(&mut first_line).expect("header line");
    assert_eq!(first_line, "##fileformat=VCFv4.5\n");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn bgzf_csi_is_selected_at_normative_limit_and_is_compressed() {
    let root = temporary("csi");
    fs::create_dir_all(&root).expect("directory");
    let output = write_vcf_contract(
        &root.join("truth"),
        VcfOutputMode::Bgzf,
        &contract(TBI_MAXIMUM_POSITION_EXCLUSIVE),
    )
    .expect("BGZF CSI VCF");
    assert_eq!(output.index_kind, Some(VcfIndexKind::Csi));
    let index = output.index.expect("CSI");
    assert_eq!(
        index.extension().and_then(|value| value.to_str()),
        Some("csi")
    );
    let mut index_reader = noodles_csi::io::Reader::new(fs::File::open(&index).expect("CSI"));
    index_reader.read_index().expect("readable CSI");
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn contract_rejects_ambiguous_or_invalid_headers() {
    assert_eq!(
        VcfHeaderContract::new(vec![], vec![]).expect_err("contigs required"),
        VcfContractError::NoContigs
    );
    for invalid in [
        "chr 1", "chr,1", "chr\"1", "chr\\1", "chr[1]", "chr<1>", "chr*1", "chr=1",
    ] {
        assert!(
            matches!(
                VcfContig::new(invalid, 1),
                Err(VcfContractError::InvalidContigIdentifier { .. })
            ),
            "{invalid} must not enter a VCF contig header"
        );
    }
    let contig = VcfContig::new("chr1", 1).expect("contig");
    assert!(matches!(
        VcfHeaderContract::new(vec![contig.clone(), contig], vec![]),
        Err(VcfContractError::DuplicateContig(_))
    ));
    assert!(matches!(
        VcfHeaderContract::new(
            vec![VcfContig::new("chr1", 1).expect("contig")],
            vec!["s".into(), "s".into()]
        ),
        Err(VcfContractError::DuplicateSample(_))
    ));
}

#[test]
fn sequence_id_backed_contigs_admit_the_project_identifier_domain() {
    let contract = VcfHeaderContract::new(
        vec![
            VcfContig::new("chr1", 12).expect("contig"),
            VcfContig::new("scaffold_02", 13).expect("contig"),
            VcfContig::new("A-1.2", TBI_MAXIMUM_POSITION_EXCLUSIVE).expect("contig"),
        ],
        vec![],
    )
    .expect("header");
    assert_eq!(contract.selected_index(), VcfIndexKind::Csi);
}

#[test]
fn index_selection_considers_every_contig_not_only_the_first() {
    let under_limit = VcfHeaderContract::new(
        vec![
            VcfContig::new("chr1", 12).expect("contig"),
            VcfContig::new("chr2", TBI_MAXIMUM_POSITION_EXCLUSIVE - 1).expect("contig"),
        ],
        vec![],
    )
    .expect("header");
    assert_eq!(under_limit.selected_index(), VcfIndexKind::Tbi);

    let at_limit = VcfHeaderContract::new(
        vec![
            VcfContig::new("chr1", 12).expect("contig"),
            VcfContig::new("chr2", TBI_MAXIMUM_POSITION_EXCLUSIVE).expect("contig"),
        ],
        vec![],
    )
    .expect("header");
    assert_eq!(at_limit.selected_index(), VcfIndexKind::Csi);
}

#[test]
fn external_validator_probe_is_informative_and_never_an_acceptance_gate() {
    let bcftools = version_line("bcftools", &["--version"]);
    let tabix = version_line("tabix", &["--version"]);
    eprintln!(
        "DTC-066 sonda informativa: bcftools={bcftools:?}, tabix={tabix:?}; la aceptación exige scripts/validate_vcf_external.sh."
    );
}

fn version_line(program: &str, args: &[&str]) -> Option<String> {
    let output = Command::new(program).args(args).output().ok()?;
    let text = String::from_utf8_lossy(&output.stdout);
    text.lines().next().map(str::to_owned)
}
