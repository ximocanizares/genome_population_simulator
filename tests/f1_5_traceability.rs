use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

#[test]
fn f1_5_closes_exactly_22_requirements_with_existing_test_symbols() {
    let matrix = include_str!("../evidence/f1-reimplementation/f1.5-test-req-matrix.tsv");
    let traceability = include_str!("../design/trazabilidad-requisitos-0.1.md");
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut test_ids = BTreeSet::new();
    let mut requirements = BTreeSet::new();

    for (line_number, line) in matrix.lines().enumerate().skip(1) {
        let columns = line.split('\t').collect::<Vec<_>>();
        assert_eq!(columns.len(), 5, "line {line_number}");
        assert!(
            test_ids.insert(columns[0]),
            "duplicate test ID {}",
            columns[0]
        );
        assert!(
            requirements.insert(columns[1]),
            "duplicate requirement {}",
            columns[1]
        );
        assert_eq!(columns[2], "passed");
        assert_eq!(columns[4], "F1.5");

        let row = traceability
            .lines()
            .find(|line| line.starts_with(&format!("| {} |", columns[1])))
            .unwrap_or_else(|| panic!("missing traceability row for {}", columns[1]));
        assert!(row.contains("| F1.5 |"), "wrong milestone: {row}");
        assert!(row.contains("Superada"), "requirement is not closed: {row}");

        for evidence in columns[3].split(';') {
            let (file, symbol) = evidence
                .split_once("::")
                .unwrap_or_else(|| panic!("invalid evidence at line {line_number}: {evidence}"));
            let source = fs::read_to_string(root.join(file))
                .unwrap_or_else(|error| panic!("cannot read {file}: {error}"));
            assert!(
                source.contains(&format!("fn {symbol}(")),
                "missing evidence symbol {symbol} in {file}"
            );
        }
    }

    assert_eq!(test_ids.len(), 22);
    assert_eq!(requirements.len(), 22);
    assert!(requirements.contains("REQ-REF-14"));
    assert!(requirements.contains("REQ-MAN-86"));

    let evidence_for = |requirement: &str| {
        matrix
            .lines()
            .skip(1)
            .find_map(|line| {
                let columns = line.split('\t').collect::<Vec<_>>();
                (columns.get(1) == Some(&requirement)).then_some(columns[3])
            })
            .unwrap_or_else(|| panic!("missing matrix row for {requirement}"))
    };
    assert!(
        evidence_for("REQ-CONF-26").contains(
            "tests/f1_5_acceptance.rs::declared_and_system_csprng_roots_are_effective_before_reference_generation"
        ),
        "REQ-CONF-26 must exercise both a declared and a system-CSPRNG root"
    );
    const DERIVED_REFERENCE_WITNESS: &str = "tests/f1_6_acceptance.rs::target_reduction_publishes_reconstructible_fasta_map_and_manifest";
    for requirement in ["REQ-REF-03", "REQ-REF-22", "REQ-SAL-01"] {
        assert!(
            evidence_for(requirement).contains(DERIVED_REFERENCE_WITNESS),
            "{requirement} must cover its derived-reference clause"
        );
    }

    let ref_18 = traceability
        .lines()
        .find(|line| line.starts_with("| REQ-REF-18 |"))
        .expect("REQ-REF-18 row");
    assert!(!ref_18.contains("event"));

    let man_86 = traceability
        .lines()
        .find(|line| line.starts_with("| REQ-MAN-86 |"))
        .expect("REQ-MAN-86 row");
    assert!(man_86.contains("cabecera/descripción reversibles"));
    assert!(!man_86.contains("ploid"));
}
