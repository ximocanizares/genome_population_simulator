use std::fs;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::truth_schema::*;
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);

fn header() -> TruthHeader {
    TruthHeader {
        document_kind: TRUTH_DOCUMENT_KIND.into(),
        schema_version: TRUTH_SCHEMA_VERSION.into(),
        execution_id: "RUN_01".into(),
        reference_id: "REF_main".into(),
        coordinate_convention: TRUTH_COORDINATE_CONVENTION.into(),
        copy_count: 2,
        genotype_count: 1,
        source_digests: vec![NamedDigest {
            source: TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "REF_main".into(),
            },
            sha256: "a".repeat(64),
        }],
        group_digests: vec![NamedDigest {
            source: TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "group-main".into(),
            },
            sha256: "b".repeat(64),
        }],
    }
}

fn copy(index: u8) -> HaplotypeCopyRecord {
    HaplotypeCopyRecord {
        identity: CopyIdentity {
            reference_id: "REF_main".into(),
            individual_id: "IND_01".into(),
            sequence_id: "chr1".into(),
            copy_index: index,
        },
        phase: index,
        provenance: CopyProvenance {
            origin: "initial_population".into(),
            sources: vec![TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "REF_main".into(),
            }],
        },
        segments: vec![TruthSegment {
            source: SegmentSource {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Reference,
                    source_id: "REF_main".into(),
                },
                reference_id: "REF_main".into(),
                sequence_id: "chr1".into(),
                start: 0,
                end: 8,
                sequence_length: 8,
            },
            modifications: vec![],
        }],
    }
}

fn genotype() -> GenotypeRecord {
    GenotypeRecord {
        case_id: "CASE_01".into(),
        reference_id: "REF_main".into(),
        individual_id: "IND_01".into(),
        locus: TruthLocus {
            sequence_id: "chr1".into(),
            start: 2,
            end: 3,
        },
        alleles: vec!["A".into(), "G".into()],
        allele_dosages: vec![1, 1],
        carriers: vec![
            AlleleCarrier {
                allele_index: 0,
                copy: copy(0).identity,
            },
            AlleleCarrier {
                allele_index: 1,
                copy: copy(1).identity,
            },
        ],
    }
}

fn valid_document() -> TruthDocument {
    TruthDocument {
        header: header(),
        records: vec![
            TruthRecord::HaplotypeCopy(copy(0)),
            TruthRecord::HaplotypeCopy(copy(1)),
            TruthRecord::Genotype(genotype()),
        ],
    }
}

#[test]
fn manual_fixture_lines_are_valid_against_the_external_schema() {
    let schema: Value =
        serde_json::from_str(include_str!("../schemas/truth-v1.schema.json")).expect("schema JSON");
    assert_eq!(
        schema["$schema"],
        "https://json-schema.org/draft/2020-12/schema"
    );
    let lines = include_str!("fixtures/truth-v1/valid.jsonl")
        .lines()
        .collect::<Vec<_>>();
    assert_eq!(lines.len(), 4);
    assert!(include_str!("fixtures/truth-v1/valid.jsonl").ends_with('\n'));
    for (index, line) in lines.iter().copied().enumerate() {
        let value: Value = serde_json::from_str(line).expect("manual fixture JSON");
        assert!(value["record_type"].is_string());
        assert_schema_result(&format!("fixture-line-{index}"), &value, true);
    }
    assert_eq!(
        serde_json::from_str::<Value>(lines[0]).expect("header")["record_type"],
        "header"
    );
    assert_eq!(
        serde_json::from_str::<Value>(lines[1]).expect("copy")["record_type"],
        "haplotype_copy"
    );
    assert_eq!(
        serde_json::from_str::<Value>(lines[2]).expect("copy")["record_type"],
        "haplotype_copy"
    );
    assert_eq!(
        serde_json::from_str::<Value>(lines[3]).expect("genotype")["record_type"],
        "genotype"
    );
}

#[test]
fn external_schema_rejects_missing_fields_and_unknown_record_keys() {
    let mut header: Value = serde_json::from_str(
        include_str!("fixtures/truth-v1/valid.jsonl")
            .lines()
            .next()
            .expect("header fixture"),
    )
    .expect("header JSON");
    assert_schema_result("header-valid", &header, true);
    header
        .as_object_mut()
        .expect("header object")
        .remove("execution_id");
    assert_schema_result("header-missing-execution", &header, false);

    let mut copy: Value = serde_json::from_str(
        include_str!("fixtures/truth-v1/valid.jsonl")
            .lines()
            .nth(1)
            .expect("copy fixture"),
    )
    .expect("copy JSON");
    copy["unexpected"] = json!(true);
    assert_schema_result("copy-unknown-field", &copy, false);
}

#[test]
fn external_parser_rejects_duplicate_keys_in_each_record_discriminator() {
    let lines = include_str!("fixtures/truth-v1/valid.jsonl")
        .lines()
        .collect::<Vec<_>>();
    for (name, line, needle, duplicate) in [
        (
            "header",
            lines[0],
            "\"execution_id\":\"RUN_01\"",
            "\"execution_id\":\"RUN_01\",\"execution_id\":\"RUN_02\"",
        ),
        (
            "copy",
            lines[1],
            "\"copy_index\":0",
            "\"copy_index\":0,\"copy_index\":1",
        ),
        (
            "genotype",
            lines[3],
            "\"start\":2",
            "\"start\":2,\"start\":1",
        ),
    ] {
        assert_raw_schema_result(name, &line.replacen(needle, duplicate, 1), false);
    }
}

#[test]
fn valid_fixture_shape_has_all_required_external_fields() {
    let records: Vec<Value> = include_str!("fixtures/truth-v1/valid.jsonl")
        .lines()
        .map(|line| serde_json::from_str(line).expect("JSON"))
        .collect();
    assert_eq!(records[0]["schema_version"], "1");
    assert!(records[1]["identity"]["copy_index"].is_u64());
    assert!(records[1]["provenance"]["sources"].is_array());
    assert_eq!(
        records[1]["segments"][0]["modifications"][0]["provenance"]["kind"],
        "matrix"
    );
    assert!(records[3]["allele_dosages"].is_array());
    assert!(records[3]["carriers"].is_array());
    assert!(valid_document().validate_structure().is_ok());
}

#[test]
fn rejects_missing_and_unknown_versions() {
    let mut missing = valid_document();
    missing.header.schema_version.clear();
    assert_eq!(
        missing.validate_structure(),
        Err(TruthContractError::UnsupportedVersion)
    );
    let mut future = valid_document();
    future.header.schema_version = "2".into();
    assert_eq!(
        future.validate_structure(),
        Err(TruthContractError::UnsupportedVersion)
    );
}

#[test]
fn rejects_duplicate_copy_and_out_of_canonical_order() {
    let mut duplicate = valid_document();
    duplicate
        .records
        .insert(2, TruthRecord::HaplotypeCopy(copy(1)));
    assert!(matches!(
        duplicate.validate_structure(),
        Err(TruthContractError::DuplicateCopy { .. })
    ));
    let mut unordered = valid_document();
    unordered.records.swap(0, 1);
    assert_eq!(
        unordered.validate_structure(),
        Err(TruthContractError::RecordsOutOfOrder)
    );
}

#[test]
fn rejects_broken_copy_genotype_relations_and_dosage() {
    let mut missing_copy = valid_document();
    if let TruthRecord::Genotype(genotype) = &mut missing_copy.records[2] {
        genotype.carriers[1].copy.copy_index = 7;
    }
    assert!(matches!(
        missing_copy.validate_structure(),
        Err(TruthContractError::UnknownCarrierCopy { .. })
    ));
    let mut wrong_scope = valid_document();
    if let TruthRecord::Genotype(genotype) = &mut wrong_scope.records[2] {
        genotype.carriers[1].copy.individual_id = "IND_02".into();
    }
    assert!(matches!(
        wrong_scope.validate_structure(),
        Err(TruthContractError::CarrierOutsideGenotype { .. })
    ));
    let mut wrong_dosage = valid_document();
    if let TruthRecord::Genotype(genotype) = &mut wrong_dosage.records[2] {
        genotype.allele_dosages[0] = 2;
    }
    assert!(matches!(
        wrong_dosage.validate_structure(),
        Err(TruthContractError::DosageMismatch { .. })
    ));
}

#[test]
fn rejects_declared_counts_and_duplicate_digests() {
    let mut counts = valid_document();
    counts.header.copy_count = 3;
    assert!(matches!(
        counts.validate_structure(),
        Err(TruthContractError::CopyCountMismatch { .. })
    ));
    let mut digests = valid_document();
    digests
        .header
        .source_digests
        .push(digests.header.source_digests[0].clone());
    assert!(matches!(
        digests.validate_structure(),
        Err(TruthContractError::DuplicateDigest { .. })
    ));
}

#[test]
fn validates_local_segment_invariants_and_typed_modification_provenance() {
    let mut document = valid_document();
    document.header.source_digests.extend([
        NamedDigest {
            source: TruthSourceRef {
                kind: TruthSourceKind::Explicit,
                source_id: "EVT_source".into(),
            },
            sha256: "c".repeat(64),
        },
        NamedDigest {
            source: TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: "MATRIX_01".into(),
            },
            sha256: "d".repeat(64),
        },
    ]);
    document.header.group_digests.push(NamedDigest {
        source: TruthSourceRef {
            kind: TruthSourceKind::Matrix,
            source_id: "POP_01".into(),
        },
        sha256: "e".repeat(64),
    });
    let explicit = ModificationProvenance {
        kind: ModificationProvenanceKind::Explicit,
        version: "1".into(),
        source: TruthSourceRef {
            kind: TruthSourceKind::Explicit,
            source_id: "EVT_source".into(),
        },
        matrix: None,
    };
    let matrix = ModificationProvenance {
        kind: ModificationProvenanceKind::Matrix,
        version: "1".into(),
        source: TruthSourceRef {
            kind: TruthSourceKind::Matrix,
            source_id: "MATRIX_01".into(),
        },
        matrix: Some(MatrixProvenance {
            matrix_source_id: "MATRIX_01".into(),
            population_id: "POP_01".into(),
            variant_id: "VAR_01".into(),
        }),
    };
    if let TruthRecord::HaplotypeCopy(copy) = &mut document.records[0] {
        copy.segments[0].modifications = vec![
            TruthModification::Snp {
                event_id: "EVT_01".into(),
                position: 1,
                reference: "A".into(),
                alternate: "G".into(),
                provenance: explicit,
            },
            TruthModification::Snp {
                event_id: "EVT_02".into(),
                position: 3,
                reference: "C".into(),
                alternate: "T".into(),
                provenance: matrix,
            },
        ];
    }
    assert!(document.validate_structure().is_ok());
    let mut bad_interval = document.clone();
    if let TruthRecord::HaplotypeCopy(copy) = &mut bad_interval.records[0] {
        copy.segments[0].source.end = 9;
    }
    assert!(matches!(
        bad_interval.validate_structure(),
        Err(TruthContractError::InvalidSegmentInterval { .. })
    ));
    let mut same_allele = document.clone();
    if let TruthRecord::HaplotypeCopy(copy) = &mut same_allele.records[0]
        && let TruthModification::Snp { alternate, .. } = &mut copy.segments[0].modifications[0]
    {
        *alternate = "A".into();
    }
    assert!(matches!(
        same_allele.validate_structure(),
        Err(TruthContractError::SnpDoesNotChangeBase { .. })
    ));
    let mut missing_digest = document;
    if let TruthRecord::HaplotypeCopy(copy) = &mut missing_digest.records[0]
        && let TruthModification::Snp { provenance, .. } = &mut copy.segments[0].modifications[1]
    {
        provenance.source.source_id = "MISSING".into();
    }
    assert!(matches!(
        missing_digest.validate_structure(),
        Err(TruthContractError::InvalidModificationProvenance { .. })
            | Err(TruthContractError::MissingDigest { .. })
    ));
}

#[test]
fn focal_contract_regressions_cover_boundaries_segments_and_groups() {
    let mut document = valid_document();
    document.header.group_digests.clear();
    assert!(
        document.validate_structure().is_ok(),
        "groups are optional without matrix provenance"
    );
    if let TruthRecord::HaplotypeCopy(copy) = &mut document.records[0] {
        copy.segments.clear();
    }
    assert!(matches!(
        document.validate_structure(),
        Err(TruthContractError::EmptySegment { .. })
    ));

    let mut terminal = valid_document();
    terminal.header.source_digests.push(NamedDigest {
        source: TruthSourceRef {
            kind: TruthSourceKind::Explicit,
            source_id: "EXPLICIT_01".into(),
        },
        sha256: "f".repeat(64),
    });
    let provenance = ModificationProvenance {
        kind: ModificationProvenanceKind::Explicit,
        version: "1".into(),
        source: TruthSourceRef {
            kind: TruthSourceKind::Explicit,
            source_id: "EXPLICIT_01".into(),
        },
        matrix: None,
    };
    if let TruthRecord::HaplotypeCopy(copy) = &mut terminal.records[0] {
        copy.segments[0].modifications = vec![TruthModification::Insertion {
            event_id: "EVT_03".into(),
            point: 0,
            bases: "A".into(),
            provenance,
        }];
    }
    assert!(matches!(
        terminal.validate_structure(),
        Err(TruthContractError::InvalidInsertion { .. })
    ));

    let mut matrix = valid_document();
    matrix.header.source_digests.push(NamedDigest {
        source: TruthSourceRef {
            kind: TruthSourceKind::Matrix,
            source_id: "MATRIX_02".into(),
        },
        sha256: "1".repeat(64),
    });
    let matrix_provenance = ModificationProvenance {
        kind: ModificationProvenanceKind::Matrix,
        version: "1".into(),
        source: TruthSourceRef {
            kind: TruthSourceKind::Matrix,
            source_id: "MATRIX_02".into(),
        },
        matrix: Some(MatrixProvenance {
            matrix_source_id: "MATRIX_02".into(),
            population_id: "POP_02".into(),
            variant_id: "VAR_02".into(),
        }),
    };
    if let TruthRecord::HaplotypeCopy(copy) = &mut matrix.records[0] {
        copy.segments[0].modifications = vec![TruthModification::Snp {
            event_id: "EVT_04".into(),
            position: 2,
            reference: "A".into(),
            alternate: "G".into(),
            provenance: matrix_provenance,
        }];
    }
    assert!(matches!(
        matrix.validate_structure(),
        Err(TruthContractError::MissingDigest { .. })
    ));
    matrix.header.group_digests.push(NamedDigest {
        source: TruthSourceRef {
            kind: TruthSourceKind::Matrix,
            source_id: "POP_02".into(),
        },
        sha256: "2".repeat(64),
    });
    assert!(matrix.validate_structure().is_ok());
    matrix
        .header
        .group_digests
        .push(matrix.header.group_digests[0].clone());
    assert!(matches!(
        matrix.validate_structure(),
        Err(TruthContractError::DuplicateDigest { .. })
    ));
}

#[test]
fn modification_order_resets_per_segment_and_keeps_the_domain_boundary_rules() {
    let mut document = valid_document();
    document.header.source_digests.push(NamedDigest {
        source: TruthSourceRef {
            kind: TruthSourceKind::Explicit,
            source_id: "EXPLICIT_02".into(),
        },
        sha256: "3".repeat(64),
    });
    let provenance = || ModificationProvenance {
        kind: ModificationProvenanceKind::Explicit,
        version: "1".into(),
        source: TruthSourceRef {
            kind: TruthSourceKind::Explicit,
            source_id: "EXPLICIT_02".into(),
        },
        matrix: None,
    };
    if let TruthRecord::HaplotypeCopy(copy) = &mut document.records[0] {
        copy.segments[0].modifications = vec![
            TruthModification::Insertion {
                event_id: "EVT_05".into(),
                point: 2,
                bases: "A".into(),
                provenance: provenance(),
            },
            TruthModification::Snp {
                event_id: "EVT_06".into(),
                position: 2,
                reference: "A".into(),
                alternate: "G".into(),
                provenance: provenance(),
            },
        ];
        copy.segments.push(TruthSegment {
            source: copy.segments[0].source.clone(),
            modifications: vec![TruthModification::Snp {
                event_id: "EVT_07".into(),
                position: 1,
                reference: "A".into(),
                alternate: "C".into(),
                provenance: provenance(),
            }],
        });
    }
    assert!(
        document.validate_structure().is_ok(),
        "insertion then SNP at one coordinate is ordered, and the next segment resets ordering"
    );
    let mut terminal_deletion = document;
    if let TruthRecord::HaplotypeCopy(copy) = &mut terminal_deletion.records[0] {
        copy.segments[0].modifications = vec![TruthModification::Deletion {
            event_id: "EVT_08".into(),
            start: 1,
            end: 8,
            provenance: provenance(),
        }];
    }
    assert!(matches!(
        terminal_deletion.validate_structure(),
        Err(TruthContractError::InvalidDeletion { .. })
    ));
}

fn assert_schema_result(name: &str, value: &Value, expected_valid: bool) {
    assert_raw_schema_result(
        name,
        &serde_json::to_string(value).expect("JSON text"),
        expected_valid,
    );
}

fn assert_raw_schema_result(name: &str, json: &str, expected_valid: bool) {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let nonce = NONCE.fetch_add(1, Ordering::Relaxed);
    let instance = std::env::temp_dir().join(format!("gps-truth-v1-{name}-{nonce}.json"));
    fs::write(&instance, json).expect("instance write");
    let status = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas/truth-v1.schema.json"))
        .arg(&instance)
        .stderr(Stdio::null())
        .status()
        .expect("external JSON Schema validator starts");
    fs::remove_file(&instance).expect("instance cleanup");
    assert_eq!(
        status.success(),
        expected_valid,
        "unexpected schema result for {name}"
    );
}
