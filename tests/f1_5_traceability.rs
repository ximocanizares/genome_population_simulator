use std::collections::BTreeSet;

#[test]
fn f1_5_reaudits_22_complete_requirements_and_6_partial_witnesses() {
    let matrix = include_str!("../evidence/f1.5/requirements.md");
    let rows = matrix
        .lines()
        .filter(|line| line.starts_with("| TEST-REQ-"))
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 28);
    assert_eq!(rows.iter().collect::<BTreeSet<_>>().len(), 28);
    assert_eq!(
        rows.iter()
            .filter(|row| row.contains("| Superado |"))
            .count(),
        22
    );
    assert_eq!(
        rows.iter()
            .filter(|row| row.contains("| Parcial |"))
            .count(),
        6
    );
}
