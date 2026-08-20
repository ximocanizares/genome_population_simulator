use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

#[test]
fn f1_4_closes_six_requirements_and_every_evidence_symbol_exists() {
    let matrix = include_str!("../evidence/f1-reimplementation/f1.4-test-req-matrix.tsv");
    let traceability = include_str!("../design/trazabilidad-requisitos-0.1.md");
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut requirements = BTreeSet::new();
    let mut tests = BTreeSet::new();
    let mut passed = 0;
    let mut partial = 0;
    for (line_number, line) in matrix.lines().enumerate().skip(1) {
        let columns = line.split('\t').collect::<Vec<_>>();
        assert_eq!(columns.len(), 5, "line {line_number}");
        assert!(tests.insert(columns[0]), "duplicate test {}", columns[0]);
        assert!(
            requirements.insert(columns[1]),
            "duplicate requirement {}",
            columns[1]
        );
        match columns[2] {
            "passed" => {
                passed += 1;
                assert_eq!(columns[4], "F1.4");
                let row = traceability
                    .lines()
                    .find(|line| line.starts_with(&format!("| {} |", columns[1])))
                    .unwrap_or_else(|| panic!("missing traceability row for {}", columns[1]));
                assert!(row.contains("| F1.4 |"));
                assert!(row.contains("Superada"));
            }
            "partial" => {
                partial += 1;
                assert!(columns[4].contains(':'));
                assert!(columns[4].contains("falta"));
            }
            status => panic!("unknown status {status}"),
        }
        let (file, symbol) = columns[3]
            .split_once("::")
            .unwrap_or_else(|| panic!("invalid evidence at line {line_number}"));
        let source = fs::read_to_string(root.join(file))
            .unwrap_or_else(|error| panic!("cannot read {file}: {error}"));
        assert!(
            source.contains(&format!("fn {symbol}(")),
            "missing evidence symbol {symbol} in {file}"
        );
    }
    assert_eq!(requirements.len(), 10);
    assert_eq!(tests.len(), 10);
    assert_eq!(passed, 6);
    assert_eq!(partial, 4);

    let man_32 = traceability
        .lines()
        .find(|line| line.starts_with("| REQ-MAN-32 |"))
        .expect("REQ-MAN-32 traceability row");
    assert!(man_32.contains("output_categories_cover_zero_one_and_multiple_reductions"));
    assert!(!man_32.contains("output_categories_and_reduction_units"));

    let sal_03 = traceability
        .lines()
        .find(|line| line.starts_with("| REQ-SAL-03 |"))
        .expect("REQ-SAL-03 traceability row");
    assert!(sal_03.contains("configuration_copy_is_byte_exact_linked_and_matches_external_sha256"));
}
