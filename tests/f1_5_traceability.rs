use std::collections::BTreeSet;

#[test]
fn f1_5_reaudits_exactly_27_requirements_with_complete_consumers() {
    let matrix = include_str!("../evidence/f1.5/requirements.md");
    let rows = matrix
        .lines()
        .filter(|line| line.starts_with("| TEST-REQ-"))
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 27);
    assert_eq!(rows.iter().collect::<BTreeSet<_>>().len(), 27);
    assert_eq!(
        rows.iter()
            .filter(|row| row.contains("| Superado |"))
            .count(),
        27
    );
    assert_eq!(
        rows.iter()
            .filter(|row| row.contains("| Parcial |"))
            .count(),
        0
    );
    assert!(
        rows.iter()
            .any(|row| row.starts_with("| TEST-REQ-REF-05 | Superado |")
                && row.contains("DTF-048"))
    );
    assert!(
        rows.iter()
            .any(|row| row.starts_with("| TEST-REQ-EVO-63 | Superado |")
                && row.contains("REQ-EVO-78"))
    );
}
