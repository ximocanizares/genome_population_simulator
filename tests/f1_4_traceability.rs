use std::collections::BTreeSet;
use std::fs;

use genome_population_simulator::application::load_and_validate_config;

#[test]
fn test_req_conf_38_validates_and_deduplicates_provided_individual_ids() {
    let root = std::env::temp_dir().join(format!("gps-f14-conf38-{}", std::process::id()));
    fs::create_dir_all(&root).expect("root");
    let base = include_str!("fixtures/config/general-minimal.toml");
    let valid = base.replacen(
        "[reference]",
        "[[individuals]]\nid = \"Sample.UPPER_1-ok\"\n\n[reference]",
        1,
    );
    let path = root.join("valid.toml");
    fs::write(&path, valid).expect("valid");
    let effective = load_and_validate_config(&path).expect("valid provided ID");
    assert_eq!(effective.individual_ids(), ["Sample.UPPER_1-ok"]);

    let invalid = base.replacen(
        "[reference]",
        "[[individuals]]\nid = \"../sample\"\n\n[reference]",
        1,
    );
    fs::write(&path, invalid).expect("invalid syntax");
    assert_eq!(
        load_and_validate_config(&path)
            .expect_err("path-like ID")
            .code(),
        "GPS-DOMAIN-001"
    );

    let duplicate = base.replacen(
        "[reference]",
        "[[individuals]]\nid = \"same\"\n[[individuals]]\nid = \"same\"\n\n[reference]",
        1,
    );
    fs::write(&path, duplicate).expect("duplicate");
    let error = load_and_validate_config(&path).expect_err("duplicate ID");
    assert_eq!(error.code(), "GPS-CONF-008");
    assert!(error.context()["detail"].contains("same"));
    fs::remove_dir_all(root).expect("cleanup");
}

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
    assert_eq!(passed, 17);
    assert_eq!(reassigned, 8);
}
