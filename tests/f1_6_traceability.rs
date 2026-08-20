use std::collections::BTreeSet;

#[test]
fn simplified_f1_assignment_contains_exactly_54_complete_requirements() {
    let traceability = include_str!("../design/trazabilidad-requisitos-0.1.md");
    let observed = traceability
        .lines()
        .filter(|line| line.starts_with("| REQ-"))
        .filter_map(|line| {
            let columns = line.split('|').map(str::trim).collect::<Vec<_>>();
            let requirement = *columns.get(1)?;
            let milestone = *columns.get(3)?;
            milestone
                .starts_with("F1.")
                .then(|| format!("{milestone}:{requirement}"))
        })
        .collect::<BTreeSet<_>>();
    let counts = [
        ("F1.2", 4_usize),
        ("F1.3", 1),
        ("F1.4", 6),
        ("F1.5", 22),
        ("F1.6", 21),
    ];
    for (milestone, expected) in counts {
        assert_eq!(
            observed
                .iter()
                .filter(|entry| entry.starts_with(&format!("{milestone}:")))
                .count(),
            expected,
        );
    }
    assert_eq!(observed.len(), 54);
}
