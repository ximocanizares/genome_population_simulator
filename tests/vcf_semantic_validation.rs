//! F2.5.6: corpus independiente de la proyección y agrupación productoras.

use std::fs;
use std::io::{Cursor, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::identity::{ExecutionId, ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::truth_schema::{
    TruthAllele, TruthCurrentPiece, TruthModification, TruthRecord,
};
use genome_population_simulator::truth_validation::{
    read_truth_jsonl, reference_content_sha256, validate_and_reconstruct_truth,
};
use genome_population_simulator::vcf_contract::{
    VcfContig, VcfHeaderContract, VcfHeaderIdentity, write_vcf_records,
};
use genome_population_simulator::vcf_grouping::group_normalized_alleles;
use genome_population_simulator::vcf_projection::project_normalized_alleles;
use genome_population_simulator::vcf_representation::prepare_conditional_vcf_records;
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
    include_str!("fixtures/truth-v4/valid.jsonl").replacen(
        "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        &reference_content_sha256(&reference()).expect("huella"),
        1,
    )
}

fn valid_vcf() -> String {
    canonical_vcf(
        "IND_01",
        "chr1\t3\t.\tA\tG\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n",
    )
}

fn canonical_vcf(samples: &str, records: &str) -> String {
    format!(
        concat!(
            "##fileformat=VCFv4.5\n",
            "##reference=REF_main\n",
            "##gps_execution_id=RUN_01\n",
            "##gps_qual=.\n",
            "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
            "##INFO=<ID=GPS_CASE,Number=.,Type=String,Description=\"Canonical truth case identifiers contributing an ALT allele\">\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n{records}"
        ),
        samples = samples,
        records = records,
    )
}

fn parsed_truth() -> genome_population_simulator::truth_schema::TruthDocument {
    read_truth_jsonl(&mut Cursor::new(truth())).expect("valid independent truth fixture")
}

fn truth_with_fully_deleted_snp() -> genome_population_simulator::truth_schema::TruthDocument {
    let mut lines = truth()
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).expect("truth JSON"))
        .collect::<Vec<_>>();
    let source = lines[1]["current_pieces"][0]["source"].clone();
    let provenance = lines[1]["history"][0]["provenance"].clone();
    lines[1]["current_pieces"] = serde_json::json!([
        {"kind":"source", "source": {
            "source": source["source"].clone(), "reference_id":"REF_main", "sequence_id":"chr1",
            "start":0, "end":2, "sequence_length":8
        }, "snps": []},
        {"kind":"source", "source": {
            "source": source["source"].clone(), "reference_id":"REF_main", "sequence_id":"chr1",
            "start":3, "end":8, "sequence_length":8
        }, "snps": []}
    ]);
    lines[1]["history"] = serde_json::json!([
        {"kind":"snp", "event_id":"EVT_01", "position":2, "provenance":provenance.clone()},
        {"kind":"deletion", "event_id":"DEL_01", "start":2, "end":3, "provenance":provenance}
    ]);
    lines[3] = serde_json::json!({
        "record_type":"genotype", "case_id":"CASE_DEL", "reference_id":"REF_main",
        "individual_id":"IND_01", "case_kind":"indel",
        "locus":{"sequence_id":"chr1", "start":2, "end":3},
        "alleles":[
            {"kind":"reference"},
            {"kind":"deletion", "event_id":"DEL_01", "start":2, "end":3,
             "source_event":{"kind":"deletion", "event_id":"DEL_01", "start":2, "end":3},
             "piece_links":[]}
        ],
        "allele_dosages":[1,1],
        "carriers":[
            {"allele_index":1, "copy":{"reference_id":"REF_main", "individual_id":"IND_01", "sequence_id":"chr1", "copy_index":0}},
            {"allele_index":0, "copy":{"reference_id":"REF_main", "individual_id":"IND_01", "sequence_id":"chr1", "copy_index":1}}
        ]
    });
    let jsonl = lines
        .into_iter()
        .map(|line| serde_json::to_string(&line).expect("serialize truth"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    read_truth_jsonl(&mut Cursor::new(jsonl)).expect("truth with a fully deleted SNP")
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
    let mut lines = include_str!("fixtures/truth-v4/valid.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).expect("JSON"))
        .collect::<Vec<_>>();
    lines[0]["reference_content_sha256"] =
        serde_json::json!(reference_content_sha256(&reference).expect("huella"));
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
        "alleles":[{"kind":"reference"},{"kind":"deletion","event_id":"DEL_67","start":6,"end":7,"source_event":{"kind":"deletion","event_id":"DEL_67","start":6,"end":7},"piece_links":[]}],
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
    let mut lines = include_str!("fixtures/truth-v4/valid.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).expect("JSON"))
        .collect::<Vec<_>>();
    lines[0]["reference_content_sha256"] =
        serde_json::json!(reference_content_sha256(&reference).expect("huella"));
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
            serde_json::json!({"kind":"insertion","event_id":"INS_02","bases":"A","source_event":{"kind":"insertion","event_id":"INS_02","point":2,"bases":"A"},"piece_links":[{"piece_index":1,"start":2,"end":3}]}),
        ),
        genotype(
            "CASE_SNP2",
            serde_json::json!({"sequence_id":"chr1","start":2,"end":3}),
            serde_json::json!({"kind":"snp","base":"A"}),
        ),
    ];
    lines.insert(0, serde_json::json!({"record_type":"header","document_kind":"gps_truth","schema_version":"4","execution_id":"RUN_01","reference_id":"REF_main","reference_content_sha256":reference_content_sha256(&reference).expect("huella"),"coordinate_convention":"zero_based_half_open_v1","counts":{"copies":1,"genotypes":3},"source_digests":[{"source":{"kind":"reference","source_id":"REF_main"},"sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"},{"source":{"kind":"matrix","source_id":"MATRIX_01"},"sha256":"cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"}],"group_digests":[{"source":{"kind":"matrix","source_id":"POP_01"},"sha256":"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"}]}));
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
fn fully_deleted_snp_is_absent_from_truth_and_vcf() {
    let truth = truth_with_fully_deleted_snp();
    assert!(
        truth
            .records
            .iter()
            .filter_map(|record| match record {
                TruthRecord::Genotype(genotype) => Some(genotype),
                TruthRecord::HaplotypeCopy(_) => None,
            })
            .all(|genotype| genotype.case_id == "CASE_DEL")
    );
    let normalized = project_normalized_alleles(&truth, &reference()).expect("project deletion");
    let grouped =
        group_normalized_alleles(&truth, &reference(), &normalized).expect("group deletion");
    let records = prepare_conditional_vcf_records(grouped)
        .expect("represent deletion")
        .records;
    let header = VcfHeaderContract::new(
        VcfHeaderIdentity::new(
            reference().id().clone(),
            ExecutionId::parse("RUN_01", "deleted SNP test").expect("execution ID"),
        ),
        vec![VcfContig::new("chr1", 8).expect("contig")],
        vec!["IND_01".to_owned()],
    )
    .expect("header");
    let vcf = String::from_utf8(
        write_vcf_records(Vec::new(), &header, &records).expect("write deletion VCF"),
    )
    .expect("UTF-8 VCF");
    assert!(!vcf.contains("\t3\t"));
    validate_vcf_against_truth(&reference(), &truth, &mut Cursor::new(vcf))
        .expect("a SNP fully deleted from current truth has no VCF record");
}

#[test]
fn independent_oracle_checks_two_alt_two_case_mapping_and_sample_dosages() {
    let vcf = canonical_vcf(
        "IND_01\tIND_02",
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
    let vcf = canonical_vcf(
        "IND_01",
        "chr1\t4\t.\tAA\tT\t.\tPASS\tGPS_CASE=CASE_DEL,CASE_SNP\tGT\t1|0\n",
    );
    validate_vcf_against_truth(&reference, &truth, &mut Cursor::new(vcf))
        .expect("SNP plus crossed left-aligned deletion is one final ALT");
}

#[test]
fn independent_oracle_keeps_an_insertion_bridge_until_all_three_changes_compose() {
    let (reference, truth) = insertion_bridge_truth();
    let vcf = canonical_vcf(
        "IND_01",
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
        let vcf = canonical_vcf(
            "IND_01",
            &format!("chr1\t3\t.\tA\tG\t.\tPASS\tGPS_CASE=CASE_01\tGT\t{gt}\n"),
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
    let vcf = canonical_vcf(
        "IND_01",
        concat!(
            "chr1\t3\t.\tA\tG\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n",
            "chr1\t6\t.\tA\tT\t.\tPASS\tGPS_CASE=CASE_01\tGT\t1|0\n",
        ),
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
            "FORMAT GT:DP",
            valid_vcf().replacen("\tGT\t1|0\n", "\tGT:DP\t1|0:10\n", 1),
        ),
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
            "missing reference identity",
            valid_vcf().replacen("##reference=REF_main\n", "", 1),
        ),
        (
            "swapped reference identity",
            valid_vcf().replacen("##reference=REF_main", "##reference=REF_other", 1),
        ),
        (
            "missing execution identity",
            valid_vcf().replacen("##gps_execution_id=RUN_01\n", "", 1),
        ),
        (
            "swapped execution identity",
            valid_vcf().replacen("##gps_execution_id=RUN_01", "##gps_execution_id=RUN_02", 1),
        ),
        (
            "missing QUAL declaration",
            valid_vcf().replacen("##gps_qual=.\n", "", 1),
        ),
        (
            "swapped QUAL declaration",
            valid_vcf().replacen("##gps_qual=.", "##gps_qual=1", 1),
        ),
        (
            "swapped FILTER declaration",
            valid_vcf().replacen("All filters passed", "Quality below 10", 1),
        ),
        (
            "missing FILTER declaration",
            valid_vcf().replacen(
                "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
                "",
                1,
            ),
        ),
        (
            "missing GT declaration",
            valid_vcf().replacen(
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "",
                1,
            ),
        ),
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

// ---------------------------------------------------------------------------
// Corpus de corrupciones **con motivo exacto**
//
// El corpus de `isolated_mutations_reject_each_f2_semantic_invariant` exige
// solo que el oráculo rechace, no **por qué**. Eso basta para lo que aquel
// cubre —cada corrupción cae por una guarda distinta y ninguna otra la caza
// antes—, pero no sirve para las guardas que este corpus añade: varias de
// ellas conviven con otras que rechazarían el mismo documento por otro motivo,
// y una prueba que solo mirase `is_err()` no distinguiría cuál trabajó.
//
// Por eso este corpus imita al de `tests/repeat_vcf_projection.rs`: afirma que
// la corrupción **cambió el documento** y exige el **motivo exacto** del
// rechazo.
// ---------------------------------------------------------------------------

/// Exige que cada corrupción cambie el documento y se rechace **por su motivo**.
fn reject_each_vcf(original: &str, corrupted: &[(&str, &str, String)]) {
    for (reason, expected_detail, content) in corrupted {
        assert_ne!(
            content, original,
            "la corrupción «{reason}» tiene que cambiar el documento"
        );
        let error = validate_vcf_against_truth(
            &reference(),
            &parsed_truth(),
            &mut Cursor::new(content.clone()),
        )
        .expect_err("el oráculo tiene que rechazar la corrupción");
        let VcfSemanticValidationError::Vcf { detail, .. } = &error else {
            panic!("«{reason}» se rechazó fuera del contrato VCF: {error:?}");
        };
        assert_eq!(
            detail, expected_detail,
            "«{reason}» se rechazó por otra razón"
        );
    }
}

/// **El contrato de cabecera y de columnas del VCF, guarda por guarda.**
///
/// *Qué hace que pueda fallar:* el documento de partida se acepta —la primera
/// aserción lo fija, para que un oráculo que rechazara todo no pasara—, y cada
/// fila exige el **motivo exacto**. Si una de estas guardas se borrase,
/// invirtiese o mirase otro campo, la fila correspondiente dejaría de citar su
/// motivo y la prueba moriría.
#[test]
fn the_header_and_column_contract_rejects_each_corruption_by_its_own_reason() {
    let original = valid_vcf();
    validate_vcf_against_truth(
        &reference(),
        &parsed_truth(),
        &mut Cursor::new(original.clone()),
    )
    .expect("el documento de partida se acepta");

    reject_each_vcf(
        &original,
        &[
            (
                "declaración FORMAT con otra descripción",
                "declaración FORMAT distinta de GT fijo en F2",
                original.replacen("Description=\"Genotype\"", "Description=\"Genotipo\"", 1),
            ),
            (
                "GPS_LOSS declarado en la cabecera",
                "F2 prohíbe declarar GPS_LOSS y GPS_COMPLEMENT",
                original.replacen(
                    "#CHROM",
                    "##INFO=<ID=GPS_LOSS,Number=.,Type=String,Description=\"x\">\n#CHROM",
                    1,
                ),
            ),
            (
                "línea de una sola almohadilla tras #CHROM",
                "una línea # no puede aparecer tras #CHROM",
                format!("{original}#cola\n"),
            ),
            (
                "cabecera de columnas con FORMAT mal escrito",
                "cabecera de columnas VCF distinta del contrato F2",
                original.replacen("\tINFO\tFORMAT\t", "\tINFO\tFORMATO\t", 1),
            ),
            (
                "línea de una sola almohadilla antes de #CHROM",
                "línea de cabecera VCF no reconocida",
                original.replacen("#CHROM", "#comentario\n#CHROM", 1),
            ),
            (
                "registro con una columna de muestra de más",
                "número de columnas VCF incompatible con las muestras",
                original.replacen("\tGT\t1|0\n", "\tGT\t1|0\t1|0\n", 1),
            ),
            (
                "registro en una posición que la verdad no declara",
                "registro CHROM+POS ausente de la verdad",
                original.replacen("chr1\t3\t", "chr1\t5\t", 1),
            ),
            (
                "INFO que no empieza por GPS_CASE",
                "INFO debe declarar GPS_CASE en primer lugar",
                original.replacen("GPS_CASE=CASE_01", "GPS_CASEX=CASE_01", 1),
            ),
            (
                "GPS_CASE sin ningún caso",
                "INFO debe contener GPS_CASE canónico sin valores vacíos",
                original.replacen("GPS_CASE=CASE_01", "GPS_CASE=", 1),
            ),
        ],
    );

    // Un VCF con la cabecera entera y **ningún registro** no cae por ninguna
    // de las guardas anteriores: su rechazo es el recuento final.
    reject_each_vcf(
        &original,
        &[(
            "cabecera completa sin ningún registro",
            "faltan registros VCF que exige la verdad",
            canonical_vcf("IND_01", ""),
        )],
    );
}

/// El cuerpo de una función, desde su firma hasta la llave de cierre de
/// columna cero.
fn function_body<'a>(source: &'a str, signature: &str) -> &'a str {
    let start = source
        .find(signature)
        .unwrap_or_else(|| panic!("{signature} tiene que existir"));
    let rest = &source[start..];
    let end = rest
        .find("\n}\n")
        .unwrap_or_else(|| panic!("{signature} tiene que terminar"));
    &rest[..end]
}

/// **Todo rechazo que este verificador decide nace en una de sus constructoras.**
///
/// *Por qué existe.* Las constructoras son el punto por el que pasa el barrido
/// de guardas, de modo que un rechazo construido a mano queda **fuera de la
/// puerta**: no se cuenta en el denominador y no se anota cuando se dispara,
/// así que podría romperse sin que nadie lo notase. Trece rechazos de este
/// módulo estaban así, uno de ellos dentro de la propia función que ata la
/// columna de referencia al genoma real.
///
/// *Guarda léxica, no semántica*: lee el texto del módulo, igual que las que la
/// auditoría de identidad usa para lo suyo. Su valor entero está en los
/// controles positivos, que impiden que pase por estar mal escrita.
///
/// *Lo que esta prueba NO cubre, dicho sin adornos.* La variante que envuelve
/// el veredicto del lector independiente de la verdad queda fuera a propósito:
/// no es un rechazo que decida este verificador, sino uno ajeno que se
/// reenvía, y sus guardas viven en el otro módulo. Quien quisiera colar un
/// rechazo propio por ahí tendría que fabricar antes un error del otro tipo.
#[test]
fn the_vcf_oracle_builds_every_rejection_through_its_constructors() {
    let source = include_str!("../src/vcf_semantic_validation.rs");
    let module = source
        .split("\n#[cfg(test)]\n")
        .next()
        .expect("el módulo tiene una parte que no es de pruebas");
    for (variant, constructor) in [
        ("VcfSemanticValidationError::Vcf {", "fn vcf_reject("),
        ("VcfSemanticValidationError::Fasta {", "fn fasta_reject("),
        ("VcfSemanticValidationError::Io {", "fn io_reject("),
        (
            "VcfSemanticValidationError::BasesNotMaterialized",
            "fn bases_not_materialized(",
        ),
    ] {
        assert_eq!(
            module.matches(variant).count(),
            1,
            "«{variant}» solo puede construirse en {constructor}"
        );
        assert!(
            function_body(module, constructor).contains(variant),
            "la única construcción de «{variant}» tiene que estar dentro de {constructor}"
        );
        // Y esa constructora tiene que participar en el barrido: sin la
        // anotación, la posición anotada sería la suya y no la de la guarda.
        let declaration = module
            .split(constructor)
            .next()
            .expect("la constructora existe");
        assert!(
            declaration.ends_with("#[track_caller]\n"),
            "{constructor} tiene que llevar #[track_caller] justo encima"
        );
    }
    // Control positivo 1: el instrumento encuentra el cuerpo correcto y
    // distingue una constructora de otra.
    assert!(function_body(module, "fn vcf_reject(").contains("line,"));
    assert!(!function_body(module, "fn vcf_reject(").contains("Fasta"));
    // Control positivo 2: la misma búsqueda sí ve una construcción a mano
    // donde la hay.
    assert_eq!(
        "let x = VcfSemanticValidationError::Vcf { line: 1, detail: d };"
            .matches("VcfSemanticValidationError::Vcf {")
            .count(),
        1
    );
    // Control positivo 3: las cinco constructoras existen y llaman al barrido.
    for constructor in [
        "fn vcf_reject(",
        "fn fasta_reject(",
        "fn io_reject(",
        "fn bases_not_materialized(",
    ] {
        assert!(
            function_body(module, constructor).contains("crate::guard_trace::record();"),
            "{constructor} tiene que anotar su disparo"
        );
    }
    // Y las tres que delegan lo hacen sin construir nada por su cuenta.
    for (delegating, target) in [
        ("fn oracle_error(", "vcf_reject(0, detail)"),
        ("fn vcf_error<T>(", "Err(vcf_reject(line, detail))"),
        ("fn io_error(", "io_reject("),
    ] {
        assert!(
            function_body(module, delegating).contains(target),
            "{delegating} tiene que delegar en su constructora"
        );
    }
}

/// **Los rechazos que hasta ahora se construían a mano, y por eso quedaban
/// fuera del barrido: ficheros ausentes, ilegibles y FASTA roto.**
///
/// *Por qué no existían.* No es que nadie hubiera querido escribirlos: el
/// enumerador de la puerta solo contaba los rechazos que pasan por una
/// constructora, y estos se construían en el sitio. No se contaban en el
/// denominador y no se anotaban al dispararse, de modo que romperlos no habría
/// movido ninguna cifra.
///
/// *Qué hace que pueda fallar:* la primera aserción fija que el trío íntegro se
/// acepta, y cada fila exige la clase de rechazo y, cuando el rechazo la fija,
/// su motivo exacto.
#[test]
fn the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta() {
    let root = temporary("path-rejections");
    fs::create_dir_all(&root).expect("directorio");
    let fasta = root.join("reference.fa");
    let truth_path = root.join("truth.jsonl");
    let vcf = root.join("truth.vcf");
    let bgzf = root.join("truth.vcf.gz");
    fs::write(&fasta, ">chr1\nAAAAAAAA\n").expect("FASTA");
    fs::write(&truth_path, truth()).expect("verdad");
    fs::write(&vcf, valid_vcf()).expect("VCF");
    let mut writer = noodles_bgzf::io::Writer::new(fs::File::create(&bgzf).expect("BGZF"));
    writer
        .write_all(valid_vcf().as_bytes())
        .expect("escribir BGZF");
    writer.finish().expect("cerrar BGZF");
    validate_published_vcf_semantics(&fasta, &truth_path, &vcf).expect("el trío íntegro se acepta");

    let missing = root.join("no-existe");
    for (reason, arguments) in [
        ("la verdad no está", (&fasta, &missing, &vcf)),
        ("el FASTA no está", (&missing, &truth_path, &vcf)),
        ("el VCF de texto no está", (&fasta, &truth_path, &missing)),
    ] {
        let error = validate_published_vcf_semantics(arguments.0, arguments.1, arguments.2)
            .expect_err(reason);
        assert!(
            matches!(error, VcfSemanticValidationError::Io { .. }),
            "«{reason}» tiene que ser un rechazo de entrada y salida: {error:?}"
        );
    }
    // El VCF comprimido lo abre otra rama, con su propio rechazo.
    let missing_bgzf = root.join("no-existe.gz");
    let error = validate_published_vcf_semantics(&fasta, &truth_path, &missing_bgzf)
        .expect_err("el VCF comprimido no está");
    assert!(matches!(error, VcfSemanticValidationError::Io { .. }));

    // Ficheros que existen pero no son texto.
    let broken_fasta = root.join("no-texto.fa");
    fs::write(&broken_fasta, b">chr1\n\xFFAAAAAAA\n".as_slice()).expect("FASTA no texto");
    let error = validate_published_vcf_semantics(&broken_fasta, &truth_path, &vcf)
        .expect_err("el FASTA no es texto");
    assert!(
        matches!(error, VcfSemanticValidationError::Io { .. }),
        "{error:?}"
    );

    let broken_vcf = root.join("no-texto.vcf");
    let mut bytes = valid_vcf().into_bytes();
    bytes.insert(1, 0xFF);
    fs::write(&broken_vcf, &bytes).expect("VCF no texto");
    let error = validate_published_vcf_semantics(&fasta, &truth_path, &broken_vcf)
        .expect_err("el VCF no es texto");
    let VcfSemanticValidationError::Io { path, .. } = &error else {
        panic!("tenía que ser un rechazo de entrada y salida: {error:?}");
    };
    assert_eq!(
        path, "VCF",
        "y el que falla es el lector del VCF, no un fichero"
    );

    // FASTA sintácticamente roto, de cuatro formas distintas.
    for (reason, contents) in [
        (
            "bases antes de la primera cabecera",
            "AAAA\n>chr1\nAAAAAAAA\n",
        ),
        ("cabecera FASTA vacía", ">\nAAAAAAAA\n"),
        ("bases que no son ADN", ">chr1\nAAAAZAAA\n"),
        ("una secuencia sin ninguna base", ">chr1\n"),
    ] {
        let path = root.join("roto.fa");
        fs::write(&path, contents).expect("FASTA roto");
        let error = validate_published_vcf_semantics(&path, &truth_path, &vcf).expect_err(reason);
        assert!(
            matches!(error, VcfSemanticValidationError::Fasta { .. }),
            "«{reason}» tiene que ser un rechazo del FASTA: {error:?}"
        );
    }

    // Una forma más de FASTA roto que no es de sintaxis de línea, sino del
    // conjunto: dos secuencias con el mismo nombre.
    let duplicated = root.join("duplicada.fa");
    fs::write(&duplicated, ">chr1\nAAAAAAAA\n>chr1\nAAAAAAAA\n").expect("FASTA duplicado");
    let error = validate_published_vcf_semantics(&duplicated, &truth_path, &vcf)
        .expect_err("dos secuencias con el mismo nombre");
    assert!(
        matches!(error, VcfSemanticValidationError::Fasta { .. }),
        "{error:?}"
    );

    // Y dos rechazos del contrato VCF que tampoco pasaban por la constructora.
    for (reason, expected_detail, contents) in [
        (
            "un registro antes de #CHROM",
            "registro VCF antes de #CHROM",
            valid_vcf().replacen("#CHROM", "chr1\t3\t.\tA\tG\t.\tPASS\tX\tGT\t1|0\n#CHROM", 1),
        ),
        (
            "POS que no es un entero positivo",
            "POS debe ser un entero positivo",
            valid_vcf().replacen("chr1\t3\t", "chr1\t0\t", 1),
        ),
    ] {
        let path = root.join("roto.vcf");
        fs::write(&path, &contents).expect("VCF roto");
        let error = validate_published_vcf_semantics(&fasta, &truth_path, &path).expect_err(reason);
        let VcfSemanticValidationError::Vcf { detail, .. } = &error else {
            panic!("«{reason}» se rechazó fuera del contrato VCF: {error:?}");
        };
        assert_eq!(
            detail, expected_detail,
            "«{reason}» se rechazó por otra razón"
        );
    }
    fs::remove_dir_all(root).expect("limpieza");
}
