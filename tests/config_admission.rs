use std::path::Path;

use genome_population_simulator::adapters::admit_toml_header;
use genome_population_simulator::configuration::{ConfigKind, ExpectedConfigKind};

fn rejected(fixture: &str, expected_kind: ExpectedConfigKind) -> (String, u8) {
    let error = admit_toml_header(
        Path::new("fixture.toml"),
        fixture.as_bytes().to_vec(),
        expected_kind,
    )
    .expect_err("fixture must be rejected");
    (error.code().to_owned(), error.exit_code())
}

#[test]
fn admitted_document_keeps_exact_bytes_path_and_class() {
    let bytes = include_bytes!("fixtures/config/general-minimal.toml").to_vec();
    let admitted = admit_toml_header(
        Path::new("general.toml"),
        bytes.clone(),
        ExpectedConfigKind::General,
    )
    .expect("valid header");

    assert_eq!(admitted.bytes(), bytes);
    assert_eq!(admitted.path(), Path::new("general.toml"));
    assert_eq!(admitted.kind(), ConfigKind::General);
}

#[test]
fn class_and_version_rejections_are_stable_and_exit_two() {
    let cases = [
        (
            include_str!("fixtures/config/invalid-missing-kind.toml"),
            ExpectedConfigKind::General,
            "GPS-CONF-002",
        ),
        (
            include_str!("fixtures/config/invalid-unknown-kind.toml"),
            ExpectedConfigKind::General,
            "GPS-CONF-003",
        ),
        (
            include_str!("fixtures/config/invalid-crossed-kind.toml"),
            ExpectedConfigKind::General,
            "GPS-CONF-003",
        ),
        (
            include_str!("fixtures/config/invalid-missing-version.toml"),
            ExpectedConfigKind::General,
            "GPS-CONF-004",
        ),
        (
            include_str!("fixtures/config/invalid-version-zero.toml"),
            ExpectedConfigKind::General,
            "GPS-CONF-005",
        ),
        (
            include_str!("fixtures/config/invalid-version-two.toml"),
            ExpectedConfigKind::General,
            "GPS-CONF-005",
        ),
        (
            include_str!("fixtures/config/invalid-version-integer.toml"),
            ExpectedConfigKind::General,
            "GPS-CONF-005",
        ),
    ];

    for (fixture, kind, expected_code) in cases {
        assert_eq!(rejected(fixture, kind), (expected_code.to_owned(), 2));
    }
}

#[test]
fn malformed_and_duplicate_toml_are_localizable() {
    let cases = [
        (
            include_str!("fixtures/config/invalid-malformed.toml"),
            "GPS-CONF-001",
        ),
        (
            include_str!("fixtures/config/invalid-duplicate-key.toml"),
            "GPS-CONF-007",
        ),
    ];
    for (fixture, expected_code) in cases {
        let error = admit_toml_header(
            Path::new("broken.toml"),
            fixture.as_bytes().to_vec(),
            ExpectedConfigKind::Either,
        )
        .expect_err("fixture must fail");
        assert_eq!(error.code(), expected_code);
        assert!(error.context().contains_key("location"));
        if expected_code == "GPS-CONF-007" {
            assert_eq!(
                error.context().get("field").map(String::as_str),
                Some("route")
            );
            assert!(error.context().contains_key("locations"));
        }
        assert_eq!(error.exit_code(), 2);
    }
}

#[test]
fn invalid_version_wins_before_unknown_semantic_content() {
    let source = b"config_kind = \"general\"\nschema_version = \"2\"\nmystery = true\n";
    assert_eq!(
        rejected(
            std::str::from_utf8(source).expect("ASCII fixture"),
            ExpectedConfigKind::General,
        ),
        ("GPS-CONF-005".to_owned(), 2)
    );
}
