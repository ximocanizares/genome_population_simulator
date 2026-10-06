//! F2.4.3: el oráculo consume solo JSONL y una referencia, nunca VCF.

use std::collections::BTreeSet;
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
    include_str!("fixtures/truth-v4/valid.jsonl").replacen(
        "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        &reference_content_sha256(&reference()).expect("huella"),
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
            {"kind":"insertion","event_id":"INS_01","bases":"CT",
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

/// DT-028: el esquema de la verdad cierra cada alelo con
/// `additionalProperties: false`, el de referencia incluido. Un alelo
/// `{"kind":"reference"}` con una clave de más se aceptaba, porque la
/// deserialización de una variante sin campos de un enum etiquetado por dentro
/// no aplica `deny_unknown_fields`. La mutación de verificación es volver a la
/// variante sin campos, y debe matar el primer caso.
#[test]
fn schema_reader_rejects_unknown_keys_in_every_allele_kind() {
    let cases = [
        (
            "reference",
            fixture(),
            "{\"kind\":\"reference\"}",
            "{\"kind\":\"reference\",\"zz_extra\":1}",
        ),
        (
            "snp",
            fixture(),
            "\"kind\":\"snp\",\"base\"",
            "\"kind\":\"snp\",\"zz_extra\":1,\"base\"",
        ),
        (
            "insertion",
            insertion_fixture(),
            // `serde_json` ordena las claves: `kind` va tras `event_id`.
            "\"event_id\":\"INS_01\",\"kind\":\"insertion\",\"piece_links\"",
            "\"event_id\":\"INS_01\",\"kind\":\"insertion\",\"zz_extra\":1,\"piece_links\"",
        ),
    ];
    for (label, text, old, new) in cases {
        read_truth_jsonl(&mut Cursor::new(text.clone())).unwrap_or_else(|error| {
            panic!("control {label}: the untouched fixture reads: {error:?}")
        });
        assert!(
            text.contains(old),
            "case {label}: the fixture carries the allele"
        );
        let corrupt = text.replacen(old, new, 1);
        let error = read_truth_jsonl(&mut Cursor::new(corrupt))
            .expect_err(&format!("case {label}: an unknown allele key must fail"));
        assert!(
            matches!(error, TruthValidationError::Schema { line: 4, .. }),
            "case {label}: {error:?}"
        );
    }
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
    // Antes lo cazaba también la estructura, contra la longitud de
    // `final_haplotype_bases`; sin ese campo, solo la copia reconstruida
    // —que mide 10 bases— puede decir que el enlace se sale de ella.
    let beyond_copy = insertion_fixture().replacen(
        "\"piece_links\":[{\"end\":4,\"piece_index\":1,\"start\":2}]",
        "\"piece_links\":[{\"end\":11,\"piece_index\":1,\"start\":2}]",
        1,
    );
    assert_ne!(
        beyond_copy,
        insertion_fixture(),
        "the corruption must apply"
    );
    read_truth_jsonl(&mut Cursor::new(beyond_copy.clone()))
        .expect("the structure alone cannot see past the end of the copy");
    assert!(
        matches!(
            validate_and_reconstruct_truth_jsonl(&mut Cursor::new(beyond_copy), &reference()),
            Err(TruthValidationError::GenotypeMismatch { .. })
        ),
        "insertion link beyond the reconstructed copy"
    );
    // Una deleción no deja fragmentos propios, así que no puede enlazar
    // ninguno. La estructura no lo prohíbe; solo la reconstrucción lo sabe.
    let mut lines = deletion_fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    assert_eq!(lines[3]["alleles"][1]["kind"], "deletion");
    lines[3]["alleles"][1]["piece_links"] = json!([{"piece_index":0,"start":0,"end":1}]);
    let linked_deletion = lines_to_jsonl(&lines);
    read_truth_jsonl(&mut Cursor::new(linked_deletion.clone()))
        .expect("the structure alone admits a deletion link");
    assert!(
        matches!(
            validate_and_reconstruct_truth_jsonl(&mut Cursor::new(linked_deletion), &reference()),
            Err(TruthValidationError::GenotypeMismatch { .. })
        ),
        "a deletion that links a piece of the reconstructed copy"
    );
}

/// **El lector no admite la copia entera del haplotipo en un alelo.**
///
/// El mismo documento sin el campo se reconstruye entero (la prueba de la
/// inserción conocida), de modo que el rechazo solo puede deberse al campo.
/// Se prueba en la inserción y en la deleción, cada una en su documento.
#[test]
fn reader_rejects_a_current_document_carrying_the_copied_final_haplotype() {
    for (name, document, copied) in [
        ("insertion", insertion_fixture(), "AACTAAAAAA"),
        ("deletion", deletion_fixture(), "AAAAAAA"),
    ] {
        validate_and_reconstruct_truth_jsonl(&mut Cursor::new(document.clone()), &reference())
            .unwrap_or_else(|error| panic!("without the field the {name} is valid: {error:?}"));
        let mut lines = document
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
            .collect::<Vec<_>>();
        assert_eq!(lines[3]["alleles"][1]["kind"], name);
        lines[3]["alleles"][1]["final_haplotype_bases"] = json!(copied);
        assert!(
            matches!(
                read_truth_jsonl(&mut Cursor::new(lines_to_jsonl(&lines))),
                Err(TruthValidationError::Schema { line: 4, .. })
            ),
            "{name}"
        );
    }
}

/// Una deleción de una base en la copia 0, con su caso genotípico.
fn deletion_fixture() -> String {
    let mut lines = fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    lines[0]["source_digests"]
        .as_array_mut()
        .expect("digests")
        .push(json!({
            "source":{"kind":"explicit","source_id":"DEL_01"},
            "sha256":"d".repeat(64),
        }));
    let provenance =
        json!({"kind":"explicit","version":"1","source":{"kind":"explicit","source_id":"DEL_01"}});
    lines[1]["segments"][0]["modifications"] = json!([]);
    let source = lines[1]["current_pieces"][0]["source"].clone();
    lines[1]["current_pieces"] = json!([source_piece(&source, 0, 3), source_piece(&source, 4, 8)]);
    lines[1]["history"] = json!([{
        "kind":"deletion","event_id":"DEL_01","start":3,"end":4,"provenance":provenance,
    }]);
    lines[3] = json!({
        "record_type":"genotype","case_id":"CASE_DEL","reference_id":"REF_main","individual_id":"IND_01",
        "case_kind":"indel","locus":{"sequence_id":"chr1","start":3,"end":4},
        "alleles":[
            {"kind":"reference"},
            {"kind":"deletion","event_id":"DEL_01","start":3,"end":4,
             "source_event":{"kind":"deletion","event_id":"DEL_01","start":3,"end":4},
             "piece_links":[]}
        ],
        "allele_dosages":[1,1],
        "carriers":[
            {"allele_index":1,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":0}},
            {"allele_index":0,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":1}}
        ]
    });
    lines_to_jsonl(&lines)
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
    // Toda deleción del historial publica su caso (DT-020).
    push_deletion_case(&mut lines, "DEL_01", (2, 3), (2, 3));
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
    lines[3]["alleles"][1]["piece_links"] = json!([{"piece_index":1,"start":2,"end":3}]);
    // La deleción solo elimina bases insertadas: su locus es el punto vacío
    // de su ancla, y su caso se publica igual (DT-020). La inserción anidada,
    // borrada entera, no tiene caso.
    push_deletion_case(&mut lines, "DEL_01", (2, 2), (3, 5));
    validate_and_reconstruct_truth_jsonl(&mut Cursor::new(lines_to_jsonl(&lines)), &reference())
        .expect("nested insertion and partial deletion replay exactly");
}

/// Añade el caso de una deleción de la copia 0 de `IND_01`, con su locus de
/// referencia y sus coordenadas vigentes, y deja los casos en orden canónico.
fn push_deletion_case(
    lines: &mut Vec<Value>,
    event_id: &str,
    (locus_start, locus_end): (u64, u64),
    (start, end): (u64, u64),
) {
    lines.push(json!({
        "record_type":"genotype","case_id":format!("CASE_{event_id}"),"reference_id":"REF_main","individual_id":"IND_01",
        "case_kind":"indel","locus":{"sequence_id":"chr1","start":locus_start,"end":locus_end},
        "alleles":[
            {"kind":"reference"},
            {"kind":"deletion","event_id":event_id,"start":start,"end":end,
             "source_event":{"kind":"deletion","event_id":event_id,"start":start,"end":end},
             "piece_links":[]}
        ],
        "allele_dosages":[1,1],
        "carriers":[
            {"allele_index":1,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":0}},
            {"allele_index":0,"copy":{"reference_id":"REF_main","individual_id":"IND_01","sequence_id":"chr1","copy_index":1}}
        ]
    }));
    let genotypes = lines[0]["counts"]["genotypes"].as_u64().expect("count") + 1;
    lines[0]["counts"]["genotypes"] = json!(genotypes);
    lines[3..].sort_by_key(|line| {
        (
            line["locus"]["start"].as_u64(),
            line["locus"]["end"].as_u64(),
            line["case_id"].as_str().map(str::to_owned),
        )
    });
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

// ---- DT-020 y DT-021: el conjunto de casos, en los dos lectores ----

/// Las dos formas de leer una verdad, completa y por flujo, sobre los mismos
/// bytes.
fn both_readers(jsonl: &str) -> [Result<(), TruthValidationError>; 2] {
    [
        validate_and_reconstruct_truth_jsonl(&mut Cursor::new(jsonl.to_owned()), &reference())
            .map(|_| ()),
        validate_truth_jsonl_streaming(&mut Cursor::new(jsonl.to_owned()), &reference(), |_, _| {})
            .map(|_| ()),
    ]
}

/// La verdad de la deleción con un segundo individuo, `IND_02`, igual que el
/// primero salvo el identificador de su deleción, `DEL_02`: así, juzgar a un
/// individuo con las copias del otro le exigiría un caso que no tiene.
fn two_individual_deletion_lines() -> Vec<Value> {
    let lines = deletion_fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    let rename = |line: &Value| {
        let text = serde_json::to_string(line)
            .expect("JSON")
            .replace("IND_01", "IND_02")
            .replace("CASE_DEL", "CASE_DEL_2")
            .replace("DEL_01", "DEL_02");
        serde_json::from_str::<Value>(&text).expect("JSON")
    };
    let mut two = vec![lines[0].clone(), lines[1].clone(), lines[2].clone()];
    two[0]["source_digests"]
        .as_array_mut()
        .expect("digests")
        .push(json!({"source":{"kind":"explicit","source_id":"DEL_02"},"sha256":"d".repeat(64)}));
    two.push(rename(&lines[1]));
    two.push(rename(&lines[2]));
    two.push(lines[3].clone());
    two.push(rename(&lines[3]));
    two[0]["counts"] = json!({"copies": 4, "genotypes": 2});
    two
}

/// Control de las dos pruebas siguientes: con dos individuos completos, los
/// dos lectores aceptan. El lector por flujo juzga a cada individuo al pasar
/// al siguiente; si lo juzgara sin sus casos, rechazaría aquí.
#[test]
fn both_readers_accept_two_complete_individuals() {
    for result in both_readers(&lines_to_jsonl(&two_individual_deletion_lines())) {
        result.expect("both individuals publish every case");
    }
}

/// DT-020. Una deleción del historial sin su caso: la copia pierde una base
/// de referencia que nada explica. Sin casos en el documento, el lector por
/// flujo solo la ve al final, al juzgar a los individuos que no aparecieron.
#[test]
fn both_readers_reject_a_deletion_without_its_case() {
    let mut lines = deletion_fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    lines.truncate(3);
    lines[0]["counts"]["genotypes"] = json!(0);
    for result in both_readers(&lines_to_jsonl(&lines)) {
        assert!(
            matches!(result, Err(TruthValidationError::MissingCase { ref individual_id, start: 3, end: 4, .. }) if individual_id == "IND_01"),
            "{result:?}"
        );
    }
}

/// DT-021. El caso de un individuo publicado dos veces con otro
/// identificador, primero en `IND_01` y después en `IND_02`. El lector por
/// flujo tiene que juzgar a `IND_01` al pasar a `IND_02`, y a `IND_02`, el
/// último, al acabar el flujo: los dos momentos, no solo uno.
#[test]
fn both_readers_reject_a_case_published_twice_in_any_individual() {
    for duplicated in [5, 6] {
        let mut lines = two_individual_deletion_lines();
        let mut duplicate = lines[duplicated].clone();
        duplicate["case_id"] = json!("CASE_DEL_DUPLICATE");
        lines.insert(duplicated + 1, duplicate);
        lines[0]["counts"]["genotypes"] = json!(3);
        for result in both_readers(&lines_to_jsonl(&lines)) {
            assert!(
                matches!(result, Err(TruthValidationError::GenotypeMismatch { ref case_id, detail })
                    if case_id.starts_with("CASE_DEL") && detail.contains("ya ese locus")),
                "duplicated line {duplicated}: {result:?}"
            );
        }
    }
}

/// Un caso *indel* cuyo evento no está en el historial de ninguna copia: el
/// contraste de portadores lo da por todo-REF —ninguna copia lleva el
/// evento— y solo el conjunto lo rechaza, porque no corresponde a nada.
#[test]
fn both_readers_reject_an_indel_case_without_its_event() {
    let mut lines = two_individual_deletion_lines();
    let mut invented = lines[5].clone();
    invented["case_id"] = json!("CASE_INVENTED");
    invented["locus"] = json!({"sequence_id":"chr1","start":5,"end":6});
    invented["alleles"][1] = json!({"kind":"deletion","event_id":"DEL_99","start":5,"end":6,
        "source_event":{"kind":"deletion","event_id":"DEL_99","start":5,"end":6},"piece_links":[]});
    invented["allele_dosages"] = json!([2, 0]);
    invented["carriers"][0]["allele_index"] = json!(0);
    lines.insert(6, invented);
    lines[0]["counts"]["genotypes"] = json!(3);
    for result in both_readers(&lines_to_jsonl(&lines)) {
        assert!(
            matches!(result, Err(TruthValidationError::GenotypeMismatch { ref case_id, detail })
                if case_id == "CASE_INVENTED" && detail.contains("no corresponde")),
            "{result:?}"
        );
    }
}

/// DT-021, segunda limitación, en los lectores: sin la matriz, el caso de un
/// cambio de base vigente solo lo exige que la copia lo lleve. Quitarlo deja
/// la base cambiada sin caso.
#[test]
fn both_readers_reject_a_base_change_without_its_case() {
    let mut lines = fixture()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    assert_eq!(lines[1]["current_pieces"][0]["snps"][0]["position"], 2);
    lines.truncate(3);
    lines[0]["counts"]["genotypes"] = json!(0);
    for result in both_readers(&lines_to_jsonl(&lines)) {
        assert!(
            matches!(result, Err(TruthValidationError::MissingCase { start: 2, end: 3, detail, .. })
                if detail.contains("cambio de base")),
            "{result:?}"
        );
    }
}

/// Un segmento recortado al principio o al final: la copia pierde bases de la
/// referencia en su estado heredado, sin ninguna deleción del historial ni
/// caso que lo explique. La copia 1 de la verdad de la deleción no lleva
/// ningún cambio, de modo que nada más la toca.
#[test]
fn both_readers_reject_a_segment_trimmed_at_either_end() {
    for (trim, expected) in [("head", (0, 2)), ("tail", (6, 8))] {
        let mut lines = deletion_fixture()
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
            .collect::<Vec<_>>();
        assert_eq!(lines[2]["identity"]["copy_index"], 1);
        let field = if trim == "head" { "start" } else { "end" };
        let value = if trim == "head" { 2 } else { 6 };
        lines[2]["segments"][0]["source"][field] = json!(value);
        lines[2]["current_pieces"][0]["source"][field] = json!(value);
        for result in both_readers(&lines_to_jsonl(&lines)) {
            assert!(
                matches!(result, Err(TruthValidationError::MissingCase { start, end, detail, .. })
                    if (start, end) == expected && detail.contains("estado heredado")),
                "{trim}: {result:?}"
            );
        }
    }
}

// ---- DT-032: ninguna estructura puede venir escrita como lista ----

/// Texto del rechazo de la comprobación nueva (`ObjectShape`). El barrido lo
/// exige en cada caso: sin eso, un caso que el propio tipo ya rechazaba por
/// otro motivo contaría como cubierto aunque la comprobación no actuase.
const LIST_REJECTION: &str = "not a list of its values";

/// Cada posición de objeto del documento: la ruta de claves e índices hasta él.
fn object_paths(value: &Value, path: &mut Vec<String>, found: &mut Vec<Vec<String>>) {
    match value {
        Value::Object(map) => {
            found.push(path.clone());
            for (key, child) in map {
                path.push(key.clone());
                object_paths(child, path, found);
                path.pop();
            }
        }
        Value::Array(items) => {
            for (index, child) in items.iter().enumerate() {
                path.push(index.to_string());
                object_paths(child, path, found);
                path.pop();
            }
        }
        _ => {}
    }
}

/// Sustituye el objeto de `path` por la lista de sus valores **en orden
/// alfabético de claves** (el de `serde_json::Map` sin `preserve_order`), que
/// no es el orden de los campos de Rust. No hace falta que `serde` aceptara esa
/// lista: lo que se exige es el rechazo de `ObjectShape`, que va antes del tipo.
fn object_as_list(line: &Value, path: &[String]) -> Value {
    let mut altered = line.clone();
    let mut node = &mut altered;
    for step in path {
        node = match node {
            Value::Object(map) => map.get_mut(step).expect("object step"),
            Value::Array(items) => &mut items[step.parse::<usize>().expect("index")],
            _ => unreachable!("a path only crosses containers"),
        };
    }
    let Value::Object(map) = node.take() else {
        unreachable!("the path ends at an object");
    };
    *node = Value::Array(map.into_iter().map(|(_, value)| value).collect());
    altered
}

/// Una verdad con su referencia, si el juez completo puede leerla.
struct SweepCorpus {
    label: &'static str,
    jsonl: String,
    with_reference: bool,
}

fn sweep_corpus() -> Vec<SweepCorpus> {
    vec![
        SweepCorpus {
            label: "snp-matrix",
            jsonl: fixture(),
            with_reference: true,
        },
        SweepCorpus {
            label: "insertion",
            jsonl: insertion_fixture(),
            with_reference: true,
        },
        SweepCorpus {
            label: "deletion",
            jsonl: deletion_fixture(),
            with_reference: true,
        },
        // Una verdad publicada de verdad: cuatro individuos con matriz, un
        // locus compuesto de tres bloques, alelos de repetición. Su
        // referencia no es el fixture de ocho bases, así que solo la lee el
        // lector de esquema.
        SweepCorpus {
            label: "repeat-population",
            jsonl: include_str!("fixtures/truth-v5/repeat-population.jsonl").to_owned(),
            with_reference: false,
        },
    ]
}

/// DT-032: `serde` acepta una estructura escrita como la lista de sus valores
/// en el orden de sus campos (`"counts":[6,8]`), y `deny_unknown_fields` no lo
/// impide. El esquema JSON la rechaza, así que el lector tiene que rechazarla.
///
/// Recorre **cada objeto de cada línea** de cuatro verdades, lo sustituye por
/// la lista de sus valores y exige el rechazo, por el lector de esquema y,
/// cuando hay referencia, también por el juez completo y el de flujo, y con el
/// texto del rechazo de `ObjectShape` (`LIST_REJECTION`), no con cualquier
/// error. Los controles sin alterar se aceptan. Lo que no ejercita: un campo
/// de lista del esquema escrito como objeto (lo vigila la prueba de nombres de
/// `src/truth_validation.rs`), ni un escalar escrito como lista (`units` como
/// lista: lo rechaza su tipo y el barrido solo sustituye objetos). Como la lista de valores sale en orden
/// alfabético y no en el de los campos de Rust, el propio tipo no la aceptaría
/// en muchos casos; da igual, porque `ObjectShape` va antes y es lo que se
/// exige.
#[test]
fn every_struct_written_as_a_list_is_rejected_by_the_truth_readers() {
    let mut variants = 0_usize;
    let mut signatures = BTreeSet::new();
    for corpus in sweep_corpus() {
        let label = corpus.label;
        read_truth_jsonl(&mut Cursor::new(corpus.jsonl.clone())).unwrap_or_else(|error| {
            panic!("control {label}: the schema reader reads it: {error:?}")
        });
        if corpus.with_reference {
            for (index, result) in both_readers(&corpus.jsonl).into_iter().enumerate() {
                result.unwrap_or_else(|error| panic!("control {label}/{index}: {error:?}"));
            }
        }
        let lines = corpus
            .jsonl
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("line JSON"))
            .collect::<Vec<_>>();
        for (line_index, line) in lines.iter().enumerate() {
            let mut paths = Vec::new();
            object_paths(line, &mut Vec::new(), &mut paths);
            for path in paths {
                let mut altered_lines = lines.clone();
                altered_lines[line_index] = object_as_list(line, &path);
                let altered = lines_to_jsonl(&altered_lines);
                let place = format!("{label} line {} at {path:?}", line_index + 1);
                let mut node = line;
                for step in &path {
                    node = match node {
                        Value::Object(map) => &map[step],
                        Value::Array(items) => &items[step.parse::<usize>().expect("index")],
                        _ => unreachable!(),
                    };
                }
                signatures.insert(
                    node.as_object()
                        .expect("object")
                        .keys()
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(","),
                );
                let error = read_truth_jsonl(&mut Cursor::new(altered.clone()))
                    .expect_err(&format!("the schema reader accepted a list at {place}"));
                assert!(
                    matches!(&error, TruthValidationError::Schema { line, detail }
                        if *line == line_index as u64 + 1 && detail.contains(LIST_REJECTION)),
                    "{place}: {error:?}"
                );
                if corpus.with_reference {
                    for (index, result) in both_readers(&altered).into_iter().enumerate() {
                        let error = result
                            .expect_err(&format!("reader {index} accepted a list at {place}"));
                        assert!(
                            matches!(&error, TruthValidationError::Schema { detail, .. }
                                if detail.contains(LIST_REJECTION)),
                            "reader {index} at {place}: {error:?}"
                        );
                    }
                }
                variants += 1;
            }
        }
    }
    // Un suelo contra un barrido vacío, y la lista de los tipos cubiertos —cada
    // uno por el conjunto ordenado de sus claves—: si un fixture deja de
    // llevar alguno, el barrido deja de probarlo.
    assert!(variants >= 340, "variants: {variants}");
    let expected: &[&str] = &[
        "allele_dosages,alleles,carriers,case_id,case_kind,individual_id,locus,record_type,reference_id",
        "allele_index,copy",
        "alternate,event_id,kind,position,provenance,reference",
        "alternate,event_id,position,provenance,reference",
        "base,kind",
        "bases,blocks,end,event_id,kind,locus_id,provenance,start",
        "bases,blocks,kind,locus_id",
        "bases,blocks,kind,locus_id,provenance,source",
        "bases,event_id,kind,piece_links,source_event",
        "bases,event_id,kind,point",
        "bases,event_id,kind,point,provenance",
        "bases,kind,provenance",
        "bases,motif,motif_length,units",
        "coordinate_convention,counts,document_kind,execution_id,group_digests,record_type,reference_content_sha256,reference_id,schema_version,source_digests",
        "copies,genotypes",
        "copy_index,individual_id,reference_id,sequence_id",
        "current_pieces,history,identity,phase,provenance,record_type,segments",
        "end,event_id,kind,piece_links,source_event,start",
        "end,event_id,kind,provenance,start",
        "end,event_id,kind,start",
        "end,piece_index,start",
        "end,reference_id,sequence_id,sequence_length,source,start",
        "end,sequence_id,start",
        "event_id,kind,position,provenance",
        "kind",
        "kind,matrix,source,version",
        "kind,snps,source",
        "kind,source,version",
        "kind,source_id",
        "matrix_source_id,population_id,variant_id",
        "modifications,source",
        "origin,sources",
        "sha256,source",
    ];
    for expected in expected {
        assert!(
            signatures.contains(*expected),
            "the corpus never turned `{expected}` into a list"
        );
    }
}
