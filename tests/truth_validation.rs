//! F2.4.3: el oráculo consume solo JSONL y una referencia, nunca VCF.

use std::io::Cursor;

use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::truth_validation::{
    TruthValidationError, read_truth_jsonl, reference_content_sha256,
    validate_and_reconstruct_truth_jsonl, validate_truth_jsonl_streaming,
};
use serde_json::{Value, json};

fn reference() -> Reference {
    reference_with_bases(b"AAAAAAAA")
}

fn reference_with_bases(bases: &[u8]) -> Reference {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "truth validation test").expect("sequence ID"),
        FastaHeader::parse(b"chr1").expect("FASTA header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(bases).expect("DNA"),
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("REF_main", "truth validation test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn fixture() -> String {
    include_str!("fixtures/truth-v1/valid.jsonl").replacen(
        "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        &reference_content_sha256(&reference()),
        1,
    )
}

fn insertion_fixture() -> String {
    let mut lines = fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    lines[0]["source_digests"]
        .as_array_mut()
        .expect("digests")
        .push(json!({
            "source":{"kind":"explicit","source_id":"INS_01"},
            "sha256":"d".repeat(64),
        }));
    let provenance =
        json!({"kind":"explicit","version":"1","source":{"kind":"explicit","source_id":"INS_01"}});
    lines[1]["segments"][0]["modifications"] = json!([]);
    let source = lines[1]["current_pieces"][0]["source"].clone();
    let first_source = json!({
        "kind":"source",
        "source": {"source":source["source"].clone(),"reference_id":"REF_main","sequence_id":"chr1","start":0,"end":2,"sequence_length":8},
        "snps":[]
    });
    let second_source = json!({
        "kind":"source",
        "source": {"source":source["source"].clone(),"reference_id":"REF_main","sequence_id":"chr1","start":2,"end":8,"sequence_length":8},
        "snps":[]
    });
    lines[1]["current_pieces"] = json!([
        first_source,
        {"kind":"inserted","bases":"CT","provenance":provenance.clone()},
        second_source,
    ]);
    lines[1]["history"] = json!([{
        "kind":"insertion","event_id":"INS_01","point":2,"bases":"CT","provenance":provenance,
    }]);
    lines[3] = json!({
        "record_type":"genotype","case_id":"CASE_INS","reference_id":"REF_main","individual_id":"IND_01",
        "case_kind":"indel","locus":{"sequence_id":"chr1","start":2,"end":2},
        "alleles":[
            {"kind":"reference"},
            {"kind":"insertion","event_id":"INS_01","bases":"CT","final_haplotype_bases":"AACTAAAAAA",
             "source_event":{"kind":"insertion","event_id":"INS_01","point":2,"bases":"CT"},
             "piece_links":[{"piece_index":1,"start":2,"end":4}]}
        ],
        "allele_dosages":[1,1],
        "carriers":[
            {"allele_index":1,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":0}},
            {"allele_index":0,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":1}}
        ]
    });
    lines
        .into_iter()
        .map(|line| serde_json::to_string(&line).expect("serialize fixture"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

#[test]
fn reconstructs_known_sequences_and_genotype_with_vcf_disabled() {
    let mut input = Cursor::new(fixture());
    let reconstructed = validate_and_reconstruct_truth_jsonl(&mut input, &reference())
        .expect("truth JSONL alone is sufficient");

    assert_eq!(
        reconstructed
            .copies
            .values()
            .map(|bases| std::str::from_utf8(bases).expect("DNA"))
            .collect::<Vec<_>>(),
        ["AAGAAAAA", "AAAAAAAA"]
    );
    assert_eq!(reconstructed.document.header.genotype_count, 1);
}

#[test]
fn streaming_validator_releases_each_copy_without_a_contract_invalid_line_cap() {
    let mut reconstructed = Vec::new();
    let summary = validate_truth_jsonl_streaming(
        &mut Cursor::new(fixture()),
        &reference(),
        |identity, bases| reconstructed.push((identity.clone(), bases.to_vec())),
    )
    .expect("streaming truth validation");
    assert_eq!(summary.copy_count, 2);
    assert_eq!(reconstructed.len(), 2);

    let oversized = vec![b'x'; 16 * 1024 * 1024 + 1];
    assert!(matches!(
        read_truth_jsonl(&mut Cursor::new(oversized)),
        Err(TruthValidationError::MissingFinalNewline)
    ));
}

#[test]
fn streaming_rejects_structural_genotype_corruption_before_semantics() {
    let bad_dosage = fixture().replacen("\"allele_dosages\":[1,1]", "\"allele_dosages\":[2,0]", 1);
    assert!(matches!(
        validate_truth_jsonl_streaming(&mut Cursor::new(bad_dosage), &reference(), |_, _| {}),
        Err(TruthValidationError::Contract(
            genome_population_simulator::truth_schema::TruthContractError::DosageMismatch { .. }
        ))
    ));
    let duplicate_carrier = fixture().replacen(
        "\"allele_index\":0,\"copy\":{\"reference_id\":\"REF_main\",\"individual_id\":\"IND_01\",\"sequence_id\":\"chr1\",\"copy_index\":1}",
        "\"allele_index\":0,\"copy\":{\"reference_id\":\"REF_main\",\"individual_id\":\"IND_01\",\"sequence_id\":\"chr1\",\"copy_index\":0}",
        1,
    );
    assert!(matches!(
        validate_truth_jsonl_streaming(
            &mut Cursor::new(duplicate_carrier),
            &reference(),
            |_, _| {}
        ),
        Err(TruthValidationError::Contract(
            genome_population_simulator::truth_schema::TruthContractError::DuplicateCarrier { .. }
        ))
    ));
}

#[test]
fn reconstructs_known_insertion_and_validates_piece_links_without_vcf() {
    let reconstructed =
        validate_and_reconstruct_truth_jsonl(&mut Cursor::new(insertion_fixture()), &reference())
            .expect("typed insertion truth is self-contained");
    assert_eq!(
        reconstructed
            .copies
            .values()
            .map(|bases| std::str::from_utf8(bases).expect("DNA"))
            .collect::<Vec<_>>(),
        ["AACTAAAAAA", "AAAAAAAA"]
    );
}

#[test]
fn schema_reader_rejects_isolated_json_schema_corruption() {
    let corrupt = fixture().replacen(
        "\"record_type\":\"header\"",
        "\"record_type\":\"header\",\"unknown\":true",
        1,
    );
    let error = read_truth_jsonl(&mut Cursor::new(corrupt)).expect_err("unknown key must fail");
    assert!(matches!(
        error,
        TruthValidationError::Schema { line: 1, .. }
    ));

    let duplicate = fixture().replacen("\"copy_index\":0", "\"copy_index\":0,\"copy_index\":1", 1);
    assert!(matches!(
        read_truth_jsonl(&mut Cursor::new(duplicate)),
        Err(TruthValidationError::Schema { line: 2, .. })
    ));
}

#[test]
fn numeric_grammar_preserves_exact_json_lexemes_and_u64_bounds() {
    let exact_float = fixture().replacen("\"copies\":2", "\"copies\":2.0", 1);
    read_truth_jsonl(&mut Cursor::new(exact_float)).expect("2.0 is a JSON Schema integer");
    let exponent = fixture().replacen("\"copies\":2", "\"copies\":2e0", 1);
    read_truth_jsonl(&mut Cursor::new(exponent)).expect("2e0 is an exact JSON integer");
    let fractional = fixture().replacen("\"copies\":2", "\"copies\":2.5", 1);
    assert!(matches!(
        read_truth_jsonl(&mut Cursor::new(fractional)),
        Err(TruthValidationError::Schema { .. })
    ));
    let fractional_exponent = fixture().replacen("\"copies\":2", "\"copies\":1e-1", 1);
    assert!(matches!(
        read_truth_jsonl(&mut Cursor::new(fractional_exponent)),
        Err(TruthValidationError::Schema { .. })
    ));
    for literal in ["-0.0", "0e999999999999999999999999999999999999999999"] {
        let zero = fixture().replacen("\"copies\":2", &format!("\"copies\":{literal}"), 1);
        assert!(matches!(
            read_truth_jsonl(&mut Cursor::new(zero)),
            Err(TruthValidationError::Contract(
                genome_population_simulator::truth_schema::TruthContractError::CopyCountMismatch {
                    declared: 0,
                    ..
                }
            ))
        ));
    }
    let above_u64 = fixture().replacen("\"copies\":2", "\"copies\":18446744073709551616", 1);
    assert!(matches!(
        read_truth_jsonl(&mut Cursor::new(above_u64)),
        Err(TruthValidationError::Schema { .. })
    ));
    let wide_exact = fixture().replacen("\"copies\":2", "\"copies\":18446744073709551615.0", 1);
    assert!(matches!(
        read_truth_jsonl(&mut Cursor::new(wide_exact)),
        Err(TruthValidationError::Contract(
            genome_population_simulator::truth_schema::TruthContractError::CopyCountMismatch {
                declared: u64::MAX,
                ..
            }
        ))
    ));
    let maximum = fixture().replacen("\"position\":2", "\"position\":18446744073709551615", 1);
    assert!(matches!(
        read_truth_jsonl(&mut Cursor::new(maximum)),
        Err(TruthValidationError::Contract(_))
    ));
}

#[test]
fn semantic_reconstructor_rejects_each_independent_relation_corruption() {
    let cases = [
        (
            "current SNP disagrees with the reference",
            fixture().replacen("\"reference\":\"A\",\"alternate\":\"G\"", "\"reference\":\"C\",\"alternate\":\"G\"", 1),
        ),
        (
            "history event is detached from historical modification",
            fixture().replacen("\"kind\":\"snp\",\"event_id\":\"EVT_01\",\"position\":2,\"provenance\"", "\"kind\":\"snp\",\"event_id\":\"EVT_02\",\"position\":2,\"provenance\"", 1),
        ),
        (
            "genotype carriers no longer match reconstructed haplotypes",
            fixture()
                .replacen("\"allele_index\":1,\"copy\":{\"reference_id\":\"REF_main\",\"individual_id\":\"IND_01\",\"sequence_id\":\"chr1\",\"copy_index\":0}", "\"allele_index\":0,\"copy\":{\"reference_id\":\"REF_main\",\"individual_id\":\"IND_01\",\"sequence_id\":\"chr1\",\"copy_index\":0}", 1)
                .replacen("\"allele_index\":0,\"copy\":{\"reference_id\":\"REF_main\",\"individual_id\":\"IND_01\",\"sequence_id\":\"chr1\",\"copy_index\":1}", "\"allele_index\":1,\"copy\":{\"reference_id\":\"REF_main\",\"individual_id\":\"IND_01\",\"sequence_id\":\"chr1\",\"copy_index\":1}", 1),
        ),
    ];
    for (name, corrupt) in cases {
        assert!(
            validate_and_reconstruct_truth_jsonl(&mut Cursor::new(corrupt), &reference(),).is_err(),
            "{name}"
        );
    }

    let bad_link = insertion_fixture().replacen("\"piece_index\":1", "\"piece_index\":0", 1);
    assert!(
        validate_and_reconstruct_truth_jsonl(&mut Cursor::new(bad_link), &reference(),).is_err(),
        "insertion links must identify the surviving inserted piece"
    );
}

#[test]
fn replay_accepts_inherited_snp_without_history_and_a_later_deleted_snp() {
    let mut lines = fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    lines[1]["history"] = json!([]);
    let source = lines[1]["current_pieces"][0]["source"].clone();
    let inherited_snp = lines[1]["current_pieces"][0]["snps"][0].clone();
    lines[1]["current_pieces"] = json!([
        source_piece(&source, 0, 2),
        {"kind":"source","source":{"source":source["source"].clone(),"reference_id":"REF_main","sequence_id":"chr1","start":2,"end":3,"sequence_length":8},"snps":[inherited_snp]},
        source_piece(&source, 3, 8),
    ]);
    let inherited_snp = lines_to_jsonl(&lines);
    validate_and_reconstruct_truth_jsonl(&mut Cursor::new(inherited_snp), &reference())
        .expect("a segment SNP is part of the inherited initial state");

    let mut lines = fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    let source = lines[1]["current_pieces"][0]["source"].clone();
    lines[1]["current_pieces"] = json!([source_piece(&source, 0, 2), source_piece(&source, 3, 8),]);
    let provenance = lines[1]["history"][0]["provenance"].clone();
    lines[1]["history"] = json!([
        {"kind":"snp","event_id":"EVT_01","position":2,"provenance":provenance.clone()},
        {"kind":"deletion","event_id":"DEL_01","start":2,"end":3,"provenance" : provenance}
    ]);
    lines[3]["allele_dosages"] = json!([1, 0]);
    lines[3]["carriers"] = json!([
        {"allele_index":0,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":1}}
    ]);
    validate_and_reconstruct_truth_jsonl(&mut Cursor::new(lines_to_jsonl(&lines)), &reference())
        .expect("the SNP is replayed, then removed by its later deletion");
}

#[test]
fn replay_rejects_permuted_or_out_of_range_history_and_changed_reference_digest() {
    let mut permuted = insertion_fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    permuted[1]["history"] = json!([{
        "kind":"deletion","event_id":"DEL_01","start":3,"end":5,
        "provenance":{"kind":"explicit","version":"1","source":{"kind":"explicit","source_id":"INS_01"}}
    }, {
        "kind":"insertion","event_id":"INS_01","point":2,"bases":"CT",
        "provenance":{"kind":"explicit","version":"1","source":{"kind":"explicit","source_id":"INS_01"}}
    }]);
    assert!(
        validate_and_reconstruct_truth_jsonl(
            &mut Cursor::new(lines_to_jsonl(&permuted)),
            &reference(),
        )
        .is_err()
    );

    let out_of_range = insertion_fixture().replacen("\"point\":2", "\"point\":99", 1);
    assert!(
        validate_and_reconstruct_truth_jsonl(&mut Cursor::new(out_of_range), &reference(),)
            .is_err()
    );

    let mut omitted_piece = insertion_fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    omitted_piece[1]["current_pieces"]
        .as_array_mut()
        .expect("pieces")
        .remove(1);
    assert!(
        validate_and_reconstruct_truth_jsonl(
            &mut Cursor::new(lines_to_jsonl(&omitted_piece)),
            &reference(),
        )
        .is_err()
    );

    let mut reordered_piece = insertion_fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    reordered_piece[1]["current_pieces"]
        .as_array_mut()
        .expect("pieces")
        .swap(0, 1);
    assert!(
        validate_and_reconstruct_truth_jsonl(
            &mut Cursor::new(lines_to_jsonl(&reordered_piece)),
            &reference(),
        )
        .is_err()
    );

    let mut invented_indel = insertion_fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    invented_indel[0]["source_digests"]
        .as_array_mut()
        .expect("digests")
        .push(json!({"source":{"kind":"explicit","source_id":"INS_02"},"sha256":"e".repeat(64)}));
    invented_indel[1]["history"]
        .as_array_mut()
        .expect("history")
        .push(json!({
            "kind":"insertion","event_id":"INS_02","point":3,"bases":"G",
            "provenance":{"kind":"explicit","version":"1","source":{"kind":"explicit","source_id":"INS_02"}}
        }));
    assert!(
        validate_and_reconstruct_truth_jsonl(
            &mut Cursor::new(lines_to_jsonl(&invented_indel)),
            &reference(),
        )
        .is_err()
    );

    assert!(matches!(
        validate_and_reconstruct_truth_jsonl(
            &mut Cursor::new(fixture()),
            &reference_with_bases(b"CAAAAAAA"),
        ),
        Err(TruthValidationError::ReferenceDigestMismatch)
    ));
}

#[test]
fn replay_accepts_nested_insertion_with_partial_survival_after_deletion() {
    let mut lines = insertion_fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    lines[0]["source_digests"]
        .as_array_mut()
        .expect("digests")
        .push(json!({
            "source":{"kind":"explicit","source_id":"INS_02"},
            "sha256":"e".repeat(64),
        }));
    let source = lines[1]["current_pieces"][0]["source"].clone();
    let first_provenance =
        json!({"kind":"explicit","version":"1","source":{"kind":"explicit","source_id":"INS_01"}});
    let nested_provenance =
        json!({"kind":"explicit","version":"1","source":{"kind":"explicit","source_id":"INS_02"}});
    lines[1]["current_pieces"] = json!([
        source_piece(&source, 0, 2),
        {"kind":"inserted","bases":"C","provenance":first_provenance.clone()},
        source_piece(&source, 2, 8),
    ]);
    lines[1]["history"] = json!([
        {"kind":"insertion","event_id":"INS_01","point":2,"bases":"CT","provenance":first_provenance},
        {"kind":"insertion","event_id":"INS_02","point":3,"bases":"G","provenance":nested_provenance.clone()},
        {"kind":"deletion","event_id":"DEL_01","start":3,"end":5,"provenance":nested_provenance},
    ]);
    lines[3]["alleles"][1]["final_haplotype_bases"] = json!("AACAAAAAA");
    lines[3]["alleles"][1]["piece_links"] = json!([{"piece_index":1,"start":2,"end":3}]);
    validate_and_reconstruct_truth_jsonl(&mut Cursor::new(lines_to_jsonl(&lines)), &reference())
        .expect("nested insertion and partial deletion replay exactly");
}

fn source_piece(source: &Value, start: u64, end: u64) -> Value {
    json!({
        "kind":"source",
        "source":{"source":source["source"].clone(),"reference_id":"REF_main","sequence_id":"chr1","start":start,"end":end,"sequence_length":8},
        "snps":[]
    })
}

fn lines_to_jsonl(lines: &[Value]) -> String {
    lines
        .iter()
        .map(|line| serde_json::to_string(line).expect("serialize fixture"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}
