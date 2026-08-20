use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

const MATRIX: &str = include_str!("../evidence/f1.2/test-req-matrix.tsv");
const TRACEABILITY: &str = include_str!("../design/trazabilidad-requisitos-0.1.md");
const ROADMAP: &str = include_str!("../design/hoja-de-ruta-0.1.md");

#[test]
fn f1_2_audits_exactly_36_test_req_cases_without_claiming_future_biology() {
    let expected = [
        "REQ-CONF-01",
        "REQ-CONF-04",
        "REQ-CONF-05",
        "REQ-CONF-07",
        "REQ-CONF-09",
        "REQ-CONF-10",
        "REQ-CONF-13",
        "REQ-CONF-14",
        "REQ-CONF-16",
        "REQ-CONF-17",
        "REQ-CONF-18",
        "REQ-CONF-20",
        "REQ-CONF-22",
        "REQ-CONF-24",
        "REQ-CONF-29",
        "REQ-CONF-30",
        "REQ-CONF-31",
        "REQ-CONF-32",
        "REQ-CONF-43",
        "REQ-CONF-65",
        "REQ-CONF-83",
        "REQ-CONF-84",
        "REQ-CONF-106",
        "REQ-CONF-116",
        "REQ-CONF-117",
        "REQ-CONF-118",
        "REQ-CONF-119",
        "REQ-CONF-120",
        "REQ-CONF-121",
        "REQ-CONF-122",
        "REQ-IND-13",
        "REQ-IND-15",
        "REQ-IND-21",
        "REQ-REP-02",
        "REQ-VAL-40",
        "REQ-VAL-51",
    ]
    .into_iter()
    .collect::<BTreeSet<_>>();
    let expected_reassignments = BTreeMap::from([
        ("REQ-CONF-07", "F8.1"),
        ("REQ-CONF-09", "F8.1"),
        ("REQ-CONF-10", "F8.1"),
        ("REQ-CONF-13", "F9.7"),
        ("REQ-CONF-16", "F9.7"),
        ("REQ-CONF-17", "F9.7"),
        ("REQ-CONF-18", "F9.7"),
        ("REQ-CONF-20", "F2.3"),
        ("REQ-CONF-22", "F2.3"),
        ("REQ-CONF-24", "F2.2"),
        ("REQ-CONF-29", "F9.7"),
        ("REQ-CONF-30", "F9.7"),
        ("REQ-CONF-31", "F9.7"),
        ("REQ-CONF-32", "F6.2"),
        ("REQ-CONF-43", "F7.6"),
        ("REQ-CONF-65", "F5.3"),
        ("REQ-CONF-83", "F4.2"),
        ("REQ-CONF-84", "F4.2"),
        ("REQ-CONF-106", "F6.2"),
        ("REQ-CONF-116", "F9.7"),
        ("REQ-CONF-117", "F9.7"),
        ("REQ-CONF-118", "F9.7"),
        ("REQ-CONF-119", "F9.7"),
        ("REQ-CONF-120", "F9.7"),
        ("REQ-CONF-121", "F9.7"),
        ("REQ-CONF-122", "F9.7"),
        ("REQ-IND-13", "F2.2"),
        ("REQ-IND-15", "F2.3"),
        ("REQ-IND-21", "F2.2"),
        ("REQ-REP-02", "F7.6"),
        ("REQ-VAL-40", "F6.2"),
        ("REQ-VAL-51", "F9.7"),
    ]);

    let rows = MATRIX.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(rows.len(), 36);
    let mut test_ids = BTreeSet::new();
    let mut requirements = BTreeSet::new();
    let mut observed_reassignments = BTreeMap::new();
    let mut passed = 0;
    for row in rows {
        let columns = row.split('\t').collect::<Vec<_>>();
        assert_eq!(columns.len(), 6, "fila mal formada: {row}");
        let [test_id, requirement, status, milestone, evidence, behavior] = columns.as_slice()
        else {
            unreachable!("la longitud ya se comprobó")
        };
        assert_eq!(*test_id, format!("TEST-{requirement}"));
        assert!(test_ids.insert(*test_id), "caso repetido: {test_id}");
        assert!(requirements.insert(*requirement));
        assert!(
            !behavior.trim().is_empty(),
            "falta conducta auditada: {test_id}"
        );
        let traceability_prefix = format!("| {requirement} |");
        let traceability_row = TRACEABILITY
            .lines()
            .find(|line| line.starts_with(&traceability_prefix))
            .unwrap_or_else(|| panic!("falta {requirement} en la trazabilidad maestra"));
        assert!(
            traceability_row.contains(&format!("| {milestone} |")),
            "la matriz y la trazabilidad discrepan para {requirement}"
        );

        match *status {
            "pass" => {
                passed += 1;
                assert_eq!(*milestone, "F1.2");
                assert_ne!(*evidence, "-");
                for linked_test in evidence.split(';') {
                    let (file, function) = linked_test
                        .split_once("::")
                        .expect("cada evidencia superada debe ser fichero::función");
                    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(file);
                    let source =
                        fs::read_to_string(&path).expect("la prueba enlazada debe existir");
                    assert!(
                        source.contains(&format!("fn {function}")),
                        "no existe {linked_test}"
                    );
                }
            }
            "reassigned" => {
                assert_eq!(
                    *evidence, "-",
                    "un pendiente no debe fingir evidencia superada"
                );
                assert!(
                    observed_reassignments
                        .insert(*requirement, *milestone)
                        .is_none(),
                    "reasignación repetida: {requirement}"
                );
                assert!(
                    ROADMAP.contains(&format!("**{milestone}. ")),
                    "el destino {milestone} de {requirement} no existe en la hoja de ruta"
                );
            }
            other => panic!("estado no admitido {other}: {test_id}"),
        }
    }
    assert_eq!(passed, 4);
    assert_eq!(observed_reassignments, expected_reassignments);
    assert_eq!(requirements, expected);
}
