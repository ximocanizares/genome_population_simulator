//! F2.5.6: corpus independiente de la proyección y agrupación productoras.

use std::fs;
use std::io::{Cursor, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::truth_schema::{
    TruthAllele, TruthCurrentPiece, TruthModification, TruthRecord,
};
use genome_population_simulator::truth_validation::{
    read_truth_jsonl, reference_content_sha256, validate_and_reconstruct_truth,
};
use genome_population_simulator::vcf_semantic_validation::{
    VcfSemanticValidationError, validate_published_vcf_semantics, validate_vcf_against_truth,
};

static NONCE: AtomicU64 = AtomicU64::new(0);

fn reference() -> Reference {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "semantic VCF test").expect("sequence ID"),
        FastaHeader::parse(b"chr1").expect("FASTA header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(b"AAAAAAAA").expect("DNA"),
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("REF_main", "semantic VCF test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn truth() -> String {
    include_str!("fixtures/truth-v1/valid.jsonl").replacen(
        "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        &reference_content_sha256(&reference()),
        1,
    )
}

fn valid_vcf() -> String {
    concat!(
        "##fileformat=VCFv4.5\n",
        "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_01\n",
        "chr1\t3\t.\tA\tG\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n",
    )
    .to_owned()
}

fn parsed_truth() -> genome_population_simulator::truth_schema::TruthDocument {
    read_truth_jsonl(&mut Cursor::new(truth())).expect("valid independent truth fixture")
}

fn two_alt_two_case_truth() -> genome_population_simulator::truth_schema::TruthDocument {
    let mut document = parsed_truth();
    let mut changed_copy = match &document.records[0] {
        TruthRecord::HaplotypeCopy(copy) => copy.clone(),
        TruthRecord::Genotype(_) => panic!("fixture begins with copy"),
    };
    changed_copy.identity.individual_id = "IND_02".to_owned();
    changed_copy.identity.copy_index = 0;
    for segment in &mut changed_copy.segments {
        for modification in &mut segment.modifications {
            if let TruthModification::Snp {
                event_id,
                alternate,
                ..
            } = modification
            {
                *event_id = "EVT_02".to_owned();
                *alternate = "T".to_owned();
            }
        }
    }
    for piece in &mut changed_copy.current_pieces {
        if let TruthCurrentPiece::Source { snps, .. } = piece {
            for snp in snps {
                snp.event_id = "EVT_02".to_owned();
                snp.alternate = "T".to_owned();
            }
        }
    }
    for event in &mut changed_copy.history {
        if let genome_population_simulator::truth_schema::TruthHistoryEvent::Snp {
            event_id, ..
        } = event
        {
            *event_id = "EVT_02".to_owned();
        }
    }

    let mut reference_copy = match &document.records[1] {
        TruthRecord::HaplotypeCopy(copy) => copy.clone(),
        TruthRecord::Genotype(_) => panic!("fixture has a second copy"),
    };
    reference_copy.identity.individual_id = "IND_02".to_owned();
    reference_copy.identity.copy_index = 1;

    let mut genotype = match &document.records[2] {
        TruthRecord::Genotype(genotype) => genotype.clone(),
        TruthRecord::HaplotypeCopy(_) => panic!("fixture ends with genotype"),
    };
    genotype.case_id = "CASE_02".to_owned();
    genotype.individual_id = "IND_02".to_owned();
    genotype.alleles[1] = TruthAllele::Snp {
        base: "T".to_owned(),
    };
    for carrier in &mut genotype.carriers {
        carrier.copy.individual_id = "IND_02".to_owned();
    }
    document.records = vec![
        TruthRecord::HaplotypeCopy(match &document.records[0] {
            TruthRecord::HaplotypeCopy(copy) => copy.clone(),
            TruthRecord::Genotype(_) => unreachable!(),
        }),
        TruthRecord::HaplotypeCopy(match &document.records[1] {
            TruthRecord::HaplotypeCopy(copy) => copy.clone(),
            TruthRecord::Genotype(_) => unreachable!(),
        }),
        TruthRecord::HaplotypeCopy(changed_copy),
        TruthRecord::HaplotypeCopy(reference_copy),
        TruthRecord::Genotype(match &document.records[2] {
            TruthRecord::Genotype(genotype) => genotype.clone(),
            TruthRecord::HaplotypeCopy(_) => unreachable!(),
        }),
        TruthRecord::Genotype(genotype),
    ];
    document.header.copy_count = 4;
    document.header.genotype_count = 2;
    validate_and_reconstruct_truth(document.clone(), &reference()).expect("two ALT truth");
    document
}

fn composed_truth() -> (
    Reference,
    genome_population_simulator::truth_schema::TruthDocument,
) {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "composed oracle").expect("sequence"),
        FastaHeader::parse(b"chr1").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(b"CAAAAAAA").expect("DNA"),
    )
    .expect("reference sequence");
    let reference = Reference::new(
        ReferenceId::parse("REF_main", "composed oracle").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference");
    let mut lines = include_str!("fixtures/truth-v1/valid.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).expect("JSON"))
        .collect::<Vec<_>>();
    lines[0]["reference_content_sha256"] = serde_json::json!(reference_content_sha256(&reference));
    lines[0]["counts"]["genotypes"] = serde_json::json!(2);
    let source = lines[1]["segments"][0]["source"].clone();
    lines[1]["segments"][0]["source"]["end"] = serde_json::json!(8);
    lines[1]["segments"][0]["source"]["sequence_length"] = serde_json::json!(8);
    let provenance = lines[1]["segments"][0]["modifications"][0]["provenance"].clone();
    lines[1]["segments"][0]["modifications"] = serde_json::json!([
        {"kind":"snp","event_id":"SNP_03","position":3,"reference":"A","alternate":"T","provenance":provenance}
    ]);
    let first_source = serde_json::json!({
        "kind":"source","source":{"source":source["source"].clone(),"reference_id":"REF_main","sequence_id":"chr1","start":0,"end":6,"sequence_length":8},
        "snps":[{"event_id":"SNP_03","position":3,"reference":"A","alternate":"T","provenance":lines[1]["segments"][0]["modifications"][0]["provenance"].clone()}]
    });
    let second_source = serde_json::json!({
        "kind":"source","source":{"source":source["source"].clone(),"reference_id":"REF_main","sequence_id":"chr1","start":7,"end":8,"sequence_length":8},"snps":[]
    });
    lines[1]["current_pieces"] = serde_json::json!([first_source, second_source]);
    lines[1]["history"] = serde_json::json!([
        {"kind":"snp","event_id":"SNP_03","position":3,"provenance":lines[1]["segments"][0]["modifications"][0]["provenance"].clone()},
        {"kind":"deletion","event_id":"DEL_67","start":6,"end":7,"provenance":lines[1]["segments"][0]["modifications"][0]["provenance"].clone()}
    ]);
    lines[2]["segments"][0]["source"]["end"] = serde_json::json!(8);
    lines[2]["segments"][0]["source"]["sequence_length"] = serde_json::json!(8);
    lines[2]["current_pieces"][0]["source"]["end"] = serde_json::json!(8);
    lines[2]["current_pieces"][0]["source"]["sequence_length"] = serde_json::json!(8);
    lines[3]["case_id"] = serde_json::json!("CASE_SNP");
    lines[3]["locus"] = serde_json::json!({"sequence_id":"chr1","start":3,"end":4});
    lines[3]["alleles"] = serde_json::json!([{"kind":"reference"},{"kind":"snp","base":"T"}]);
    let deletion = serde_json::json!({
        "record_type":"genotype","case_id":"CASE_DEL","reference_id":"REF_main","individual_id":"IND_01","case_kind":"indel","locus":{"sequence_id":"chr1","start":6,"end":7},
        "alleles":[{"kind":"reference"},{"kind":"deletion","event_id":"DEL_67","start":6,"end":7,"final_haplotype_bases":"CAATAAA","source_event":{"kind":"deletion","event_id":"DEL_67","start":6,"end":7},"piece_links":[]}],
        "allele_dosages":[1,1],"carriers":[{"allele_index":1,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":0}},{"allele_index":0,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":1}}]
    });
    lines.push(deletion);
    let jsonl = lines
        .into_iter()
        .map(|line| serde_json::to_string(&line).expect("serialize"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    (
        reference,
        read_truth_jsonl(&mut Cursor::new(jsonl)).expect("composed truth"),
    )
}

fn one_case_several_records_truth() -> genome_population_simulator::truth_schema::TruthDocument {
    let mut document = parsed_truth();
    let mut copy = match &document.records[0] {
        TruthRecord::HaplotypeCopy(copy) => copy.clone(),
        TruthRecord::Genotype(_) => unreachable!(),
    };
    let provenance = match &copy.segments[0].modifications[0] {
        TruthModification::Snp { provenance, .. } => provenance.clone(),
        _ => unreachable!(),
    };
    copy.segments[0].modifications.push(TruthModification::Snp {
        event_id: "EVT_05".to_owned(),
        position: 5,
        reference: "A".to_owned(),
        alternate: "T".to_owned(),
        provenance: provenance.clone(),
    });
    if let TruthCurrentPiece::Source { snps, .. } = &mut copy.current_pieces[0] {
        snps.push(genome_population_simulator::truth_schema::TruthCurrentSnp {
            event_id: "EVT_05".to_owned(),
            position: 5,
            reference: "A".to_owned(),
            alternate: "T".to_owned(),
            provenance: provenance.clone(),
        });
    }
    copy.history.push(
        genome_population_simulator::truth_schema::TruthHistoryEvent::Snp {
            event_id: "EVT_05".to_owned(),
            position: 5,
            provenance,
        },
    );
    let mut second_case = match &document.records[2] {
        TruthRecord::Genotype(genotype) => genotype.clone(),
        TruthRecord::HaplotypeCopy(_) => unreachable!(),
    };
    second_case.locus.start = 5;
    second_case.locus.end = 6;
    second_case.alleles[1] = TruthAllele::Snp {
        base: "T".to_owned(),
    };
    document.records[0] = TruthRecord::HaplotypeCopy(copy);
    document.records.push(TruthRecord::Genotype(second_case));
    document.header.genotype_count = 2;
    validate_and_reconstruct_truth(document.clone(), &reference()).expect("one case two loci");
    document
}

fn insertion_bridge_truth() -> (
    Reference,
    genome_population_simulator::truth_schema::TruthDocument,
) {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "insertion bridge").expect("sequence"),
        FastaHeader::parse(b"chr1").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(b"ACG").expect("DNA"),
    )
    .expect("reference sequence");
    let reference = Reference::new(
        ReferenceId::parse("REF_main", "insertion bridge").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference");
    let mut lines = include_str!("fixtures/truth-v1/valid.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).expect("JSON"))
        .collect::<Vec<_>>();
    lines[0]["reference_content_sha256"] = serde_json::json!(reference_content_sha256(&reference));
    lines[0]["counts"] = serde_json::json!({"copies":1,"genotypes":3});
    let source = serde_json::json!({"source":{"kind":"reference","source_id":"REF_main"},"reference_id":"REF_main","sequence_id":"chr1","start":0,"end":3,"sequence_length":3});
    let provenance = lines[1]["segments"][0]["modifications"][0]["provenance"].clone();
    lines[1]["segments"] = serde_json::json!([{"source":source,"modifications":[
        {"kind":"snp","event_id":"SNP_01","position":1,"reference":"C","alternate":"T","provenance":provenance},
        {"kind":"snp","event_id":"SNP_02","position":2,"reference":"G","alternate":"A","provenance":lines[1]["segments"][0]["modifications"][0]["provenance"].clone()}
    ]}]);
    let source_02 = serde_json::json!({"source":{"kind":"reference","source_id":"REF_main"},"reference_id":"REF_main","sequence_id":"chr1","start":0,"end":2,"sequence_length":3});
    let source_23 = serde_json::json!({"source":{"kind":"reference","source_id":"REF_main"},"reference_id":"REF_main","sequence_id":"chr1","start":2,"end":3,"sequence_length":3});
    lines[1]["current_pieces"] = serde_json::json!([
        {"kind":"source","source":source_02,"snps":[{"event_id":"SNP_01","position":1,"reference":"C","alternate":"T","provenance":lines[1]["segments"][0]["modifications"][0]["provenance"].clone()}]},
        {"kind":"inserted","bases":"A","provenance":lines[1]["segments"][0]["modifications"][0]["provenance"].clone()},
        {"kind":"source","source":source_23,"snps":[{"event_id":"SNP_02","position":2,"reference":"G","alternate":"A","provenance":lines[1]["segments"][0]["modifications"][1]["provenance"].clone()}]}
    ]);
    lines[1]["history"] = serde_json::json!([
        {"kind":"snp","event_id":"SNP_01","position":1,"provenance":lines[1]["segments"][0]["modifications"][0]["provenance"].clone()},
        {"kind":"snp","event_id":"SNP_02","position":2,"provenance":lines[1]["segments"][0]["modifications"][1]["provenance"].clone()},
        {"kind":"insertion","event_id":"INS_02","point":2,"bases":"A","provenance":lines[1]["segments"][0]["modifications"][0]["provenance"].clone()}
    ]);
    let copy = lines[1].clone();
    let carrier = serde_json::json!({"allele_index":1,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":0}});
    let genotype = |case_id: &str, locus: serde_json::Value, allele: serde_json::Value| {
        serde_json::json!({
            "record_type":"genotype","case_id":case_id,"reference_id":"REF_main","individual_id":"IND_01","case_kind": if case_id == "CASE_INS" { "indel" } else { "snp" },"locus":locus,
            "alleles":[{"kind":"reference"},allele],"allele_dosages":[0,1],"carriers":[carrier.clone()]
        })
    };
    lines = vec![
        copy,
        genotype(
            "CASE_SNP1",
            serde_json::json!({"sequence_id":"chr1","start":1,"end":2}),
            serde_json::json!({"kind":"snp","base":"T"}),
        ),
        genotype(
            "CASE_INS",
            serde_json::json!({"sequence_id":"chr1","start":2,"end":2}),
            serde_json::json!({"kind":"insertion","event_id":"INS_02","bases":"A","final_haplotype_bases":"ATAA","source_event":{"kind":"insertion","event_id":"INS_02","point":2,"bases":"A"},"piece_links":[{"piece_index":1,"start":2,"end":3}]}),
        ),
        genotype(
            "CASE_SNP2",
            serde_json::json!({"sequence_id":"chr1","start":2,"end":3}),
            serde_json::json!({"kind":"snp","base":"A"}),
        ),
    ];
    lines.insert(0, serde_json::json!({"record_type":"header","document_kind":"gps_truth","schema_version":"1","execution_id":"RUN_01","reference_id":"REF_main","reference_content_sha256":reference_content_sha256(&reference),"coordinate_convention":"zero_based_half_open_v1","counts":{"copies":1,"genotypes":3},"source_digests":[{"source":{"kind":"reference","source_id":"REF_main"},"sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"},{"source":{"kind":"matrix","source_id":"MATRIX_01"},"sha256":"cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"}],"group_digests":[{"source":{"kind":"matrix","source_id":"POP_01"},"sha256":"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"}]}));
    let jsonl = lines
        .into_iter()
        .map(|line| serde_json::to_string(&line).expect("serialize"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    (
        reference,
        read_truth_jsonl(&mut Cursor::new(jsonl)).expect("bridge truth"),
    )
}

fn truth_with_ploidy(ploidy: u8) -> genome_population_simulator::truth_schema::TruthDocument {
    let source = parsed_truth();
    let altered = match &source.records[0] {
        TruthRecord::HaplotypeCopy(copy) => copy.clone(),
        TruthRecord::Genotype(_) => unreachable!(),
    };
    let mut genotype = match &source.records[2] {
        TruthRecord::Genotype(genotype) => genotype.clone(),
        TruthRecord::HaplotypeCopy(_) => unreachable!(),
    };
    let copies = (0..ploidy)
        .map(|index| {
            let mut copy = altered.clone();
            copy.identity.copy_index = index;
            copy.phase = index;
            TruthRecord::HaplotypeCopy(copy)
        })
        .collect::<Vec<_>>();
    genotype.allele_dosages = vec![0, ploidy];
    genotype.carriers = (0..ploidy)
        .map(|index| {
            let mut carrier = genotype.carriers[0].clone();
            carrier.allele_index = 1;
            carrier.copy.copy_index = index;
            carrier
        })
        .collect();
    let mut document = genome_population_simulator::truth_schema::TruthDocument {
        header: source.header,
        records: copies,
    };
    document.records.push(TruthRecord::Genotype(genotype));
    document.header.copy_count = u64::from(ploidy);
    validate_and_reconstruct_truth(document.clone(), &reference()).expect("ploidy truth");
    document
}

#[test]
fn independent_oracle_accepts_ref_alt_gt_phase_dosage_samples_and_gps_case() {
    validate_vcf_against_truth(&reference(), &parsed_truth(), &mut Cursor::new(valid_vcf()))
        .expect("literal VCF agrees with truth");
}

#[test]
fn independent_oracle_checks_two_alt_two_case_mapping_and_sample_dosages() {
    let vcf = concat!(
        "##fileformat=VCFv4.5\n",
        "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_01\tIND_02\n",
        "chr1\t3\t.\tA\tG,T\t.\tPASS\tGPS_CASE=CASE_01,CASE_02\tGT\t1|0\t2|0\n",
    );
    validate_vcf_against_truth(
        &reference(),
        &two_alt_two_case_truth(),
        &mut Cursor::new(vcf),
    )
    .expect("multi-ALT cases and phased doses");
}

#[test]
fn independent_oracle_composes_final_haplotype_components_before_case_mapping() {
    let (reference, truth) = composed_truth();
    let vcf = concat!(
        "##fileformat=VCFv4.5\n",
        "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_01\n",
        "chr1\t4\t.\tAA\tT\t.\tPASS\tGPS_CASE=CASE_DEL,CASE_SNP\tGT\t1|0\n",
    );
    validate_vcf_against_truth(&reference, &truth, &mut Cursor::new(vcf))
        .expect("SNP plus crossed left-aligned deletion is one final ALT");
}

#[test]
fn independent_oracle_keeps_an_insertion_bridge_until_all_three_changes_compose() {
    let (reference, truth) = insertion_bridge_truth();
    let vcf = concat!(
        "##fileformat=VCFv4.5\n",
        "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_01\n",
        "chr1\t2\t.\tCG\tTAA\t.\tPASS\tGPS_CASE=CASE_INS,CASE_SNP1,CASE_SNP2\tGT\t1\n",
    );
    validate_vcf_against_truth(&reference, &truth, &mut Cursor::new(vcf))
        .expect("right zero-width insertion bridge joins both SNP components");
}

#[test]
fn independent_oracle_accepts_fixed_alt_for_haploid_diploid_and_decaploid_truth() {
    for ploidy in [1_u8, 2, 10] {
        let gt = std::iter::repeat_n("1", usize::from(ploidy))
            .collect::<Vec<_>>()
            .join("|");
        let vcf = format!(
            "##fileformat=VCFv4.5\n##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_01\nchr1\t3\t.\tA\tG\t.\tPASS\tGPS_CASE=CASE_01\tGT\t{gt}\n"
        );
        validate_vcf_against_truth(
            &reference(),
            &truth_with_ploidy(ploidy),
            &mut Cursor::new(vcf),
        )
        .expect("fixed ALT across supported ploidies");
    }
}

#[test]
fn independent_oracle_allows_one_case_to_repeat_on_several_records() {
    let vcf = concat!(
        "##fileformat=VCFv4.5\n",
        "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_01\n",
        "chr1\t3\t.\tA\tG\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n",
        "chr1\t6\t.\tA\tT\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n",
    );
    validate_vcf_against_truth(
        &reference(),
        &one_case_several_records_truth(),
        &mut Cursor::new(vcf),
    )
    .expect("one source case repeats for its two records");
}

#[test]
fn isolated_mutations_reject_each_f2_semantic_invariant() {
    for (label, mutated) in [
        ("REF", valid_vcf().replacen("\tA\tG\t", "\tC\tG\t", 1)),
        ("ALT", valid_vcf().replacen("\tA\tG\t", "\tA\tT\t", 1)),
        ("GT", valid_vcf().replacen("\t1|0\n", "\t0|1\n", 1)),
        ("dosage", valid_vcf().replacen("\t1|0\n", "\t1|1\n", 1)),
        ("phase", valid_vcf().replacen("\t1|0\n", "\t1/0\n", 1)),
        (
            "GPS_CASE",
            valid_vcf().replacen("CASE_01", "CASE_UNKNOWN", 1),
        ),
        (
            "prohibited complement",
            valid_vcf().replacen("GPS_CASE=CASE_01", "GPS_CASE=CASE_01;GPS_COMPLEMENT=x", 1),
        ),
        (
            "prohibited loss",
            valid_vcf().replacen("GPS_CASE=CASE_01", "GPS_CASE=CASE_01;GPS_LOSS=x", 1),
        ),
        (
            "duplicate position",
            format!(
                "{}{}",
                valid_vcf(),
                "chr1\t3\t.\tA\tG\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n"
            ),
        ),
        (
            "missing record",
            valid_vcf().lines().take(3).collect::<Vec<_>>().join("\n") + "\n",
        ),
        ("ID", valid_vcf().replacen("\t.\tA\t", "\tID_1\tA\t", 1)),
        (
            "QUAL",
            valid_vcf().replacen("\tG\t.\tPASS", "\tG\t1\tPASS", 1),
        ),
        (
            "FILTER",
            valid_vcf().replacen("\t.\tPASS\tGPS", "\t.\tq10\tGPS", 1),
        ),
        (
            "INFO extra",
            valid_vcf().replacen("GPS_CASE=CASE_01", "GPS_CASE=CASE_01;X=1", 1),
        ),
        (
            "INFO duplicate",
            valid_vcf().replacen("GPS_CASE=CASE_01", "GPS_CASE=CASE_01,CASE_01", 1),
        ),
        ("VCF version", valid_vcf().replacen("VCFv4.5", "VCFv4.4", 1)),
        (
            "GPS_CASE header Number/Type",
            valid_vcf().replacen("Number=.", "Number=1", 1).replacen(
                "Type=String",
                "Type=Integer",
                1,
            ),
        ),
        (
            "GPS_CASE header description",
            valid_vcf().replacen(
                "Canonical truth case identifiers contributing an ALT allele",
                "misleading description",
                1,
            ),
        ),
        (
            "duplicate GPS_CASE header",
            valid_vcf().replacen(
                "#CHROM",
                "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n#CHROM",
                1,
            ),
        ),
    ] {
        assert!(
            matches!(
                validate_vcf_against_truth(
                    &reference(),
                    &parsed_truth(),
                    &mut Cursor::new(mutated)
                ),
                Err(VcfSemanticValidationError::Vcf { .. })
            ),
            "{label}"
        );
    }
}

#[test]
fn independent_oracle_rejects_sample_order_record_order_and_disordered_gps_case() {
    let samples = concat!(
        "##fileformat=VCFv4.5\n",
        "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_02\tIND_01\n",
        "chr1\t3\t.\tA\tG,T\t.\tPASS\tGPS_CASE=CASE_01,CASE_02\tGT\t2|0\t1|0\n",
    );
    assert!(
        validate_vcf_against_truth(
            &reference(),
            &two_alt_two_case_truth(),
            &mut Cursor::new(samples),
        )
        .is_err()
    );
    let reversed = concat!(
        "##fileformat=VCFv4.5\n",
        "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_01\n",
        "chr1\t6\t.\tA\tT\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n",
        "chr1\t3\t.\tA\tG\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n",
    );
    assert!(
        validate_vcf_against_truth(
            &reference(),
            &one_case_several_records_truth(),
            &mut Cursor::new(reversed),
        )
        .is_err()
    );
    let (reference, truth) = composed_truth();
    let gps = concat!(
        "##fileformat=VCFv4.5\n",
        "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND_01\n",
        "chr1\t4\t.\tAA\tT\t.\tPASS\tGPS_CASE=CASE_SNP,CASE_DEL\tGT\t1|0\n",
    );
    assert!(validate_vcf_against_truth(&reference, &truth, &mut Cursor::new(gps)).is_err());
}

#[test]
fn independent_phase_check_rejects_noncontiguous_copy_indices_before_truth_reconstruction() {
    let mut truth = parsed_truth();
    for (index, record) in truth.records.iter_mut().enumerate() {
        if let TruthRecord::HaplotypeCopy(copy) = record {
            copy.phase = if index == 0 { 0 } else { 2 };
        }
    }
    assert!(matches!(
        validate_vcf_against_truth(&reference(), &truth, &mut Cursor::new(valid_vcf())),
        Err(VcfSemanticValidationError::Vcf { line: 0, .. })
    ));

    let mut haploid = parsed_truth();
    haploid.records.retain(|record| {
        !matches!(record, TruthRecord::HaplotypeCopy(copy) if copy.identity.copy_index == 1)
    });
    if let TruthRecord::HaplotypeCopy(copy) = &mut haploid.records[0] {
        copy.phase = 1;
    }
    if let TruthRecord::Genotype(genotype) = &mut haploid.records[1] {
        genotype.allele_dosages = vec![0, 1];
        genotype
            .carriers
            .retain(|carrier| carrier.allele_index == 1);
    }
    haploid.header.copy_count = 1;
    assert!(matches!(
        validate_vcf_against_truth(&reference(), &haploid, &mut Cursor::new(valid_vcf())),
        Err(VcfSemanticValidationError::Vcf { line: 0, .. })
    ));
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f2-5-6-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

#[test]
fn path_oracle_reads_truth_before_opening_vcf() {
    let root = temporary("truth-first");
    fs::create_dir_all(&root).expect("directory");
    let fasta = root.join("reference.fa");
    let bad_truth = root.join("truth.jsonl");
    fs::write(&fasta, ">chr1\nAAAAAAAA\n").expect("FASTA");
    fs::write(&bad_truth, b"{not json}\n").expect("bad truth");

    let error = validate_published_vcf_semantics(&fasta, &bad_truth, &root.join("absent.vcf"))
        .expect_err("a VCF must never supply or repair truth");
    assert!(matches!(error, VcfSemanticValidationError::Truth(_)));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn path_oracle_reconstructs_truth_against_fasta_before_opening_vcf() {
    let root = temporary("reference-before-vcf");
    fs::create_dir_all(&root).expect("directory");
    let fasta = root.join("reference.fa");
    let truth_path = root.join("truth.jsonl");
    fs::write(&fasta, ">chr1\nCAAAAAAA\n").expect("different FASTA");
    fs::write(&truth_path, truth()).expect("valid truth for A FASTA");

    assert!(matches!(
        validate_published_vcf_semantics(&fasta, &truth_path, &root.join("absent.vcf")),
        Err(VcfSemanticValidationError::Truth(
            genome_population_simulator::truth_validation::TruthValidationError::ReferenceDigestMismatch
        ))
    ));
    fs::remove_dir_all(root).expect("cleanup");
}

#[test]
fn parser_rejects_a_gps_case_declaration_after_a_record() {
    let vcf = format!(
        "{}##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
        valid_vcf()
    );
    assert!(matches!(
        validate_vcf_against_truth(&reference(), &parsed_truth(), &mut Cursor::new(vcf)),
        Err(VcfSemanticValidationError::Vcf { .. })
    ));
}

#[test]
fn path_oracle_accepts_the_published_text_and_bgzf_forms() {
    let root = temporary("published-forms");
    fs::create_dir_all(&root).expect("directory");
    let fasta = root.join("reference.fa");
    let truth_path = root.join("truth.jsonl");
    let text = root.join("truth.vcf");
    let bgzf = root.join("truth.vcf.gz");
    fs::write(&fasta, ">chr1\nAAAAAAAA\n").expect("FASTA");
    fs::write(&truth_path, truth()).expect("truth");
    fs::write(&text, valid_vcf()).expect("text VCF");
    let mut writer = noodles_bgzf::io::Writer::new(fs::File::create(&bgzf).expect("BGZF"));
    writer
        .write_all(valid_vcf().as_bytes())
        .expect("write BGZF");
    writer.finish().expect("finish BGZF");

    validate_published_vcf_semantics(&fasta, &truth_path, &text).expect("text VCF");
    validate_published_vcf_semantics(&fasta, &truth_path, &bgzf).expect("BGZF VCF");
    fs::remove_dir_all(root).expect("cleanup");
}
