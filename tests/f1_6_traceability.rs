use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

const MATRIX: &str =
    include_str!("../evidence/f1-reimplementation/f1-complete-test-req-matrix.tsv");
const TRACEABILITY: &str = include_str!("../design/trazabilidad-requisitos-0.1.md");
const EXECUTED_TESTS: &str = include_str!("../evidence/f1.6/executed-test-list.txt");

#[test]
fn simplified_f1_expedient_closes_exactly_54_requirements_with_real_symbols() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let expected_counts = BTreeMap::from([
        ("F1.2", 4_usize),
        ("F1.3", 1),
        ("F1.4", 6),
        ("F1.5", 22),
        ("F1.6", 21),
    ]);
    let mut test_ids = BTreeSet::new();
    let mut requirements = BTreeSet::new();
    let mut observed_counts = BTreeMap::<&str, usize>::new();

    for (line_number, line) in MATRIX.lines().enumerate().skip(1) {
        let columns = line.split('\t').collect::<Vec<_>>();
        assert_eq!(columns.len(), 5, "line {line_number}: {line}");
        let [test_id, requirement, status, milestone, evidence] = columns.as_slice() else {
            unreachable!("column count already checked")
        };
        assert!(test_ids.insert(*test_id), "duplicate test ID {test_id}");
        assert!(
            requirements.insert(*requirement),
            "duplicate requirement {requirement}"
        );
        assert_eq!(*test_id, format!("TEST-{requirement}"));
        assert_eq!(*status, "passed");
        assert!(expected_counts.contains_key(milestone));
        *observed_counts.entry(milestone).or_default() += 1;

        let trace_row = TRACEABILITY
            .lines()
            .find(|row| row.starts_with(&format!("| {requirement} |")))
            .unwrap_or_else(|| panic!("missing central trace row for {requirement}"));
        assert!(trace_row.contains(&format!("| {milestone} |")));
        assert!(
            trace_row.contains("Superada"),
            "requirement is not closed: {trace_row}"
        );

        for witness in evidence.split(';') {
            let (file, symbol) = witness
                .split_once("::")
                .unwrap_or_else(|| panic!("invalid witness {witness}"));
            let source = fs::read_to_string(root.join(file))
                .unwrap_or_else(|error| panic!("cannot read {file}: {error}"));
            assert!(
                source.contains(&format!("fn {symbol}(")),
                "missing executable symbol {symbol} in {file}"
            );
        }
    }

    assert_eq!(test_ids.len(), 54);
    assert_eq!(requirements.len(), 54);
    assert_eq!(observed_counts, expected_counts);
    assert!(!requirements.contains("REQ-CONF-126"));

    assert!(EXECUTED_TESTS.starts_with(
        "evidence_version=2\nexecuted_at=2026-08-26\ncommand=cargo test --all-targets --all-features\nscope=all executable witnesses cited by the current F1.6 requirement matrix\n"
    ));
    let listed_entries = EXECUTED_TESTS
        .lines()
        .filter(|line| line.contains("::"))
        .collect::<Vec<_>>();
    let listed_symbols = listed_entries
        .iter()
        .map(|line| line.rsplit_once("::").expect("qualified test name").1)
        .collect::<BTreeSet<_>>();
    assert_eq!(
        listed_entries.len(),
        listed_symbols.len(),
        "the executed list must not repeat a witness"
    );
    let mut expected_listed_symbols = MATRIX
        .lines()
        .skip(1)
        .filter(|line| line.split('\t').nth(3) == Some("F1.6"))
        .flat_map(|line| line.split('\t').nth(4).expect("matrix evidence").split(';'))
        .map(|witness| witness.rsplit_once("::").expect("matrix witness").1)
        .collect::<BTreeSet<_>>();
    expected_listed_symbols
        .insert("simplified_f1_expedient_closes_exactly_54_requirements_with_real_symbols");
    assert_eq!(
        listed_symbols, expected_listed_symbols,
        "the executed list must contain exactly the live F1.6 witnesses"
    );

    let central_f1 = TRACEABILITY
        .lines()
        .filter(|line| line.starts_with("| REQ-"))
        .filter_map(|line| {
            let columns = line.split('|').map(str::trim).collect::<Vec<_>>();
            let requirement = *columns.get(1)?;
            let milestone = *columns.get(3)?;
            milestone.starts_with("F1.").then_some(requirement)
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(central_f1, requirements);

    let f1_6_rows = MATRIX
        .lines()
        .skip(1)
        .filter(|line| line.split('\t').nth(3) == Some("F1.6"))
        .collect::<Vec<_>>()
        .join("\n");
    for removed_claim in [
        "ResolutionDomain",
        "resolution_evidence",
        "resolution-evidence",
        "symlink",
        "hardlink",
        "TOCTOU",
    ] {
        assert!(
            !f1_6_rows.contains(removed_claim),
            "retired claim in current F1.6 expedient: {removed_claim}"
        );
    }
}
