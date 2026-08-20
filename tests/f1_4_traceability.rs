use std::collections::BTreeSet;

#[test]
fn f1_4_reaudits_exactly_25_requirements_without_claiming_missing_consumers() {
    let matrix = include_str!("../evidence/f1.4/test-req-matrix.tsv");
    let mut requirements = BTreeSet::new();
    let mut tests = BTreeSet::new();
    let mut passed = 0;
    let mut reassigned = 0;
    for (line_number, line) in matrix.lines().enumerate().skip(1) {
        let columns = line.split('\t').collect::<Vec<_>>();
        assert_eq!(columns.len(), 5, "line {line_number}");
        assert!(tests.insert(columns[0]), "duplicate test {}", columns[0]);
        assert!(
            requirements.insert(columns[1]),
            "duplicate REQ {}",
            columns[1]
        );
        match columns[2] {
            "passed" => {
                passed += 1;
                assert_eq!(columns[4], "F1.4");
            }
            "reassigned" => {
                reassigned += 1;
                assert!(columns[4].contains(':'));
                assert!(columns[4].contains("falta") || columns[4].contains("faltan"));
            }
            status => panic!("unknown status {status}"),
        }
        assert!(columns[0].starts_with("TEST-REQ-"));
        assert!(columns[3].starts_with("tests/"));
    }
    assert_eq!(requirements.len(), 25);
    assert_eq!(tests.len(), 25);
    assert_eq!(passed, 6);
    assert_eq!(reassigned, 19);
}
