use genome_population_simulator::rules::{
    ParameterKey, RegionSelector, Rule, RuleProvenance, RuleValue, SelectorIdentityCatalog,
    SelectorSet, SpatialSelector, validate_rules,
};

fn rule(id: &str, value: f64, selectors: SelectorSet) -> Rule {
    Rule::new(
        id,
        ParameterKey::VariantRate,
        RuleValue::Float(value),
        selectors,
        RuleProvenance {
            path: "rules.toml".to_owned(),
            location: id.to_owned(),
        },
    )
    .expect("manual selector is allowed")
}

fn identities() -> SelectorIdentityCatalog {
    let mut catalog = SelectorIdentityCatalog::default();
    catalog.add_reference("ref");
    catalog.add_sequence("ref", "chr1");
    catalog.add_individual("founder");
    catalog.add_technical_event("event");
    catalog
}

fn region(start: u64, end: u64) -> SelectorSet {
    SelectorSet {
        spatial: Some(SpatialSelector::Region(RegionSelector {
            reference: "ref".to_owned(),
            sequence: "chr1".to_owned(),
            start,
            end,
        })),
        ..SelectorSet::default()
    }
}

#[test]
fn exact_duplicates_fail_with_equal_or_different_values() {
    for second_value in [1.0, 2.0] {
        let error = validate_rules(
            &[
                rule("first", 1.0, region(1, 10)),
                rule("second", second_value, region(1, 10)),
            ],
            &identities(),
        )
        .expect_err("exact selectors must be unique");
        assert_eq!(error.code(), "GPS-CONF-010");
        assert_eq!(error.context()["parameter"], "variant_rate");
        assert!(error.context()["selectors"].contains("Region"));
        assert!(error.context()["scope"].contains("Region"));
        assert!(error.context()["locations"].contains("first"));
        assert!(error.context()["locations"].contains("second"));
    }
}

#[test]
fn closed_regions_reject_partial_total_and_boundary_overlap() {
    for second in [region(5, 15), region(2, 9), region(10, 20)] {
        let error = validate_rules(
            &[
                rule("first", 1.0, region(1, 10)),
                rule("second", 1.0, second),
            ],
            &identities(),
        )
        .expect_err("closed overlap must fail");
        assert_eq!(error.code(), "GPS-CONF-011");
        assert_eq!(error.context()["parameter"], "variant_rate");
        assert_eq!(error.context()["reference"], "ref");
        assert_eq!(error.context()["sequence"], "chr1");
        assert!(error.context()["scope"].contains("Region"));
        assert!(error.context()["regions"].contains(';'));
    }
}

#[test]
fn adjacent_regions_are_valid_and_maximum_end_does_not_overflow() {
    validate_rules(
        &[
            rule("first", 1.0, region(1, 10)),
            rule("second", 1.0, region(11, u64::MAX)),
        ],
        &identities(),
    )
    .expect("adjacency is not overlap");
}

#[test]
fn future_missing_and_wrong_category_identities_are_rejected() {
    let missing = rule(
        "future",
        1.0,
        SelectorSet {
            individual: Some("future-child".to_owned()),
            ..SelectorSet::default()
        },
    );
    assert_eq!(
        validate_rules(&[missing], &identities())
            .expect_err("future individual must fail")
            .code(),
        "GPS-CONF-013"
    );

    let wrong_category = rule(
        "wrong-category",
        1.0,
        SelectorSet {
            individual: Some("event".to_owned()),
            ..SelectorSet::default()
        },
    );
    let error = validate_rules(&[wrong_category], &identities())
        .expect_err("event ID is not an individual ID");
    assert_eq!(error.code(), "GPS-CONF-013");
    assert_eq!(
        error.context().get("found_category").map(String::as_str),
        Some("technical_event")
    );
}

#[test]
fn unknown_reference_and_sequence_are_rejected() {
    for selectors in [
        SelectorSet {
            spatial: Some(SpatialSelector::Sequence {
                reference: "other".to_owned(),
                sequence: "chr1".to_owned(),
            }),
            ..SelectorSet::default()
        },
        SelectorSet {
            spatial: Some(SpatialSelector::Sequence {
                reference: "ref".to_owned(),
                sequence: "chr2".to_owned(),
            }),
            ..SelectorSet::default()
        },
    ] {
        assert_eq!(
            validate_rules(&[rule("invalid", 1.0, selectors)], &identities())
                .expect_err("unknown spatial identity must fail")
                .code(),
            "GPS-CONF-013"
        );
    }
}
