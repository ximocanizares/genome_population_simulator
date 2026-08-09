use genome_population_simulator::rules::{
    ParameterKey, PopulationSelector, RegionSelector, ResolutionContext, Rule, RuleProvenance,
    RuleSet, RuleValue, SelectorSet, SpatialSelector, TemporalSelector,
};

fn rule(id: &str, selectors: SelectorSet) -> Rule {
    Rule::new(
        id,
        ParameterKey::VariantRate,
        RuleValue::Float(id.len() as f64),
        selectors,
        RuleProvenance {
            path: "manual".to_owned(),
            location: id.to_owned(),
        },
    )
    .expect("manual selectors are valid")
}

#[test]
fn all_selectors_are_conjunctive_and_match_their_literal_scope() {
    let rules = RuleSet::new(
        ParameterKey::VariantRate,
        vec![
            rule("base", SelectorSet::default()),
            rule(
                "population",
                SelectorSet {
                    population: Some(PopulationSelector::Population("pop".to_owned())),
                    ..SelectorSet::default()
                },
            ),
            rule(
                "phase",
                SelectorSet {
                    temporal: Some(TemporalSelector::Phase("phase".to_owned())),
                    ..SelectorSet::default()
                },
            ),
            rule(
                "sequence",
                SelectorSet {
                    spatial: Some(SpatialSelector::Sequence {
                        reference: "ref".to_owned(),
                        sequence: "chr1".to_owned(),
                    }),
                    ..SelectorSet::default()
                },
            ),
            rule(
                "region",
                SelectorSet {
                    spatial: Some(SpatialSelector::Region(RegionSelector {
                        reference: "ref".to_owned(),
                        sequence: "chr1".to_owned(),
                        start: 10,
                        end: 20,
                    })),
                    ..SelectorSet::default()
                },
            ),
            rule(
                "locus",
                SelectorSet {
                    spatial: Some(SpatialSelector::Locus {
                        reference: "ref".to_owned(),
                        sequence: "chr1".to_owned(),
                        locus: "locus".to_owned(),
                    }),
                    ..SelectorSet::default()
                },
            ),
            rule(
                "combined",
                SelectorSet {
                    population: Some(PopulationSelector::Population("pop".to_owned())),
                    temporal: Some(TemporalSelector::Phase("phase".to_owned())),
                    ..SelectorSet::default()
                },
            ),
        ],
    )
    .expect("base exists");
    let context = ResolutionContext {
        population: Some("pop".to_owned()),
        phase: Some("phase".to_owned()),
        reference: Some("ref".to_owned()),
        sequence: Some("chr1".to_owned()),
        position: Some(10),
        locus: Some("locus".to_owned()),
        ..ResolutionContext::default()
    };
    let matching = rules
        .matching_rules(&context)
        .into_iter()
        .map(|rule| rule.id())
        .collect::<Vec<_>>();

    assert_eq!(
        matching,
        [
            "base",
            "population",
            "phase",
            "sequence",
            "region",
            "locus",
            "combined"
        ]
    );
}

#[test]
fn region_boundaries_are_closed_and_reference_qualified() {
    let region = rule(
        "region",
        SelectorSet {
            spatial: Some(SpatialSelector::Region(RegionSelector {
                reference: "ref".to_owned(),
                sequence: "chr1".to_owned(),
                start: 10,
                end: 20,
            })),
            ..SelectorSet::default()
        },
    );
    for position in [10, 20] {
        assert!(region.matches(&ResolutionContext {
            reference: Some("ref".to_owned()),
            sequence: Some("chr1".to_owned()),
            position: Some(position),
            ..ResolutionContext::default()
        }));
    }
    assert!(!region.matches(&ResolutionContext {
        reference: Some("other".to_owned()),
        sequence: Some("chr1".to_owned()),
        position: Some(15),
        ..ResolutionContext::default()
    }));
}

#[test]
fn base_rule_is_mandatory_for_a_resolvable_parameter() {
    let error = RuleSet::new(
        ParameterKey::VariantRate,
        vec![rule(
            "population",
            SelectorSet {
                population: Some(PopulationSelector::Population("pop".to_owned())),
                ..SelectorSet::default()
            },
        )],
    )
    .expect_err("missing base must fail");
    assert_eq!(error.code(), "GPS-CONF-008");
}
