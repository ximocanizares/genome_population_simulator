use genome_population_simulator::rules::{
    ParameterKey, PopulationSelector, RegionSelector, ResolutionContext, Rule, RuleProvenance,
    RuleSet, RuleValue, SelectorSet, SpatialSelector, TemporalSelector, dominates,
};

fn rule(id: &str, selectors: SelectorSet) -> Rule {
    parameter_rule(id, ParameterKey::VariantRate, selectors)
}

fn parameter_rule(id: &str, parameter: ParameterKey, selectors: SelectorSet) -> Rule {
    Rule::new(
        id,
        parameter,
        RuleValue::Float(1.0),
        selectors,
        RuleProvenance {
            path: "lattice".to_owned(),
            location: id.to_owned(),
        },
    )
    .expect("valid manual rule")
}

#[test]
fn every_internal_hierarchy_resolves_in_an_applicable_context_without_cross_axis_priority() {
    let temporal = RuleSet::new(
        ParameterKey::VariantRate,
        vec![
            rule("base", SelectorSet::default()),
            rule("phase", phase()),
            rule(
                "generation",
                SelectorSet {
                    temporal: Some(TemporalSelector::Generation(7)),
                    ..SelectorSet::default()
                },
            ),
        ],
    )
    .expect("temporal lattice");
    assert_eq!(
        temporal
            .resolve(&ResolutionContext {
                phase: Some("phase".to_owned()),
                generation: Some(7),
                ..ResolutionContext::default()
            })
            .expect("generation dominates phase")
            .id(),
        "generation"
    );

    let population_rules = RuleSet::new(
        ParameterKey::VariantRate,
        vec![
            rule("base", SelectorSet::default()),
            rule("population", population()),
            rule(
                "subpopulation",
                SelectorSet {
                    population: Some(PopulationSelector::Subpopulation("sub".to_owned())),
                    ..SelectorSet::default()
                },
            ),
        ],
    )
    .expect("population lattice");
    assert_eq!(
        population_rules
            .resolve(&ResolutionContext {
                population: Some("pop".to_owned()),
                subpopulation: Some("sub".to_owned()),
                ..ResolutionContext::default()
            })
            .expect("subpopulation dominates population")
            .id(),
        "subpopulation"
    );

    let spatial = RuleSet::new(
        ParameterKey::VariantRate,
        vec![
            rule("base", SelectorSet::default()),
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
                        locus: "locus-15".to_owned(),
                    }),
                    ..SelectorSet::default()
                },
            ),
        ],
    )
    .expect("spatial lattice");
    assert_eq!(
        spatial
            .resolve(&ResolutionContext {
                reference: Some("ref".to_owned()),
                sequence: Some("chr1".to_owned()),
                position: Some(15),
                locus: Some("locus-15".to_owned()),
                ..ResolutionContext::default()
            })
            .expect("locus dominates region and sequence")
            .id(),
        "locus"
    );

    for (id, selectors, context) in [
        (
            "individual",
            SelectorSet {
                individual: Some("founder".to_owned()),
                ..SelectorSet::default()
            },
            ResolutionContext {
                individual: Some("founder".to_owned()),
                ..ResolutionContext::default()
            },
        ),
        (
            "group",
            SelectorSet {
                group: Some("group".to_owned()),
                ..SelectorSet::default()
            },
            ResolutionContext {
                group: Some("group".to_owned()),
                ..ResolutionContext::default()
            },
        ),
        (
            "block",
            SelectorSet {
                block: Some("block".to_owned()),
                ..SelectorSet::default()
            },
            ResolutionContext {
                block: Some("block".to_owned()),
                ..ResolutionContext::default()
            },
        ),
    ] {
        let set = RuleSet::new(
            ParameterKey::VariantRate,
            vec![rule("base", SelectorSet::default()), rule(id, selectors)],
        )
        .expect("nominal lattice");
        assert_eq!(
            set.resolve(&context).expect("concrete beats absence").id(),
            id
        );
    }

    let technical = RuleSet::new(
        ParameterKey::Coverage,
        vec![
            parameter_rule("base", ParameterKey::Coverage, SelectorSet::default()),
            parameter_rule(
                "event",
                ParameterKey::Coverage,
                SelectorSet {
                    technical_event: Some("event".to_owned()),
                    ..SelectorSet::default()
                },
            ),
        ],
    )
    .expect("technical lattice");
    assert_eq!(
        technical
            .resolve(&ResolutionContext {
                technical_event: Some("event".to_owned()),
                ..ResolutionContext::default()
            })
            .expect("technical event beats absence")
            .id(),
        "event"
    );

    let incomparable = RuleSet::new(
        ParameterKey::VariantRate,
        vec![
            rule("base", SelectorSet::default()),
            rule("population", population()),
            rule("phase", phase()),
        ],
    )
    .expect("cross-axis lattice");
    assert_eq!(
        incomparable
            .resolve(&context())
            .expect_err("cross-axis rules have no implicit priority")
            .context()["rules"],
        "phase,population"
    );
}

#[test]
fn spatial_and_individual_rules_require_a_combined_rule_for_their_intersection() {
    let region = SelectorSet {
        spatial: Some(SpatialSelector::Region(RegionSelector {
            reference: "ref".to_owned(),
            sequence: "chr1".to_owned(),
            start: 10,
            end: 20,
        })),
        ..SelectorSet::default()
    };
    let individual = SelectorSet {
        individual: Some("founder".to_owned()),
        ..SelectorSet::default()
    };
    let combined = SelectorSet {
        individual: Some("founder".to_owned()),
        spatial: region.spatial.clone(),
        ..SelectorSet::default()
    };
    let context = ResolutionContext {
        individual: Some("founder".to_owned()),
        reference: Some("ref".to_owned()),
        sequence: Some("chr1".to_owned()),
        position: Some(15),
        ..ResolutionContext::default()
    };

    let conflict = RuleSet::new(
        ParameterKey::VariantRate,
        vec![
            rule("base", SelectorSet::default()),
            rule("individual", individual.clone()),
            rule("region", region.clone()),
        ],
    )
    .expect("base exists")
    .resolve(&context)
    .expect_err("spatial and individual are incomparable");
    assert_eq!(conflict.context()["rules"], "individual,region");

    for selectors in [
        vec![
            ("base", SelectorSet::default()),
            ("individual", individual.clone()),
            ("region", region.clone()),
            ("combined", combined.clone()),
        ],
        vec![
            ("combined", combined.clone()),
            ("region", region.clone()),
            ("individual", individual.clone()),
            ("base", SelectorSet::default()),
        ],
    ] {
        let set = RuleSet::new(
            ParameterKey::VariantRate,
            selectors
                .into_iter()
                .map(|(id, selectors)| rule(id, selectors))
                .collect(),
        )
        .expect("combined lattice");
        assert_eq!(
            set.resolve(&context).expect("unique maximum").id(),
            "combined"
        );
    }
}

fn population() -> SelectorSet {
    SelectorSet {
        population: Some(PopulationSelector::Population("pop".to_owned())),
        ..SelectorSet::default()
    }
}

fn phase() -> SelectorSet {
    SelectorSet {
        temporal: Some(TemporalSelector::Phase("phase".to_owned())),
        ..SelectorSet::default()
    }
}

fn combined() -> SelectorSet {
    SelectorSet {
        population: Some(PopulationSelector::Population("pop".to_owned())),
        temporal: Some(TemporalSelector::Phase("phase".to_owned())),
        ..SelectorSet::default()
    }
}

fn context() -> ResolutionContext {
    ResolutionContext {
        population: Some("pop".to_owned()),
        phase: Some("phase".to_owned()),
        ..ResolutionContext::default()
    }
}

#[test]
fn manual_lattice_has_only_the_declared_partial_order() {
    let base = SelectorSet::default();
    let population = population();
    let phase = phase();
    let combined = combined();

    assert!(dominates(&population, &base));
    assert!(dominates(&phase, &base));
    assert!(dominates(&combined, &population));
    assert!(dominates(&combined, &phase));
    assert!(!dominates(&population, &phase));
    assert!(!dominates(&phase, &population));
    assert!(!dominates(&base, &base));
}

#[test]
fn combined_rule_is_the_unique_maximum_independent_of_input_order() {
    let forward = vec![
        rule("base", SelectorSet::default()),
        rule("population", population()),
        rule("phase", phase()),
        rule("combined", combined()),
    ];
    let reverse = forward.iter().cloned().rev().collect::<Vec<_>>();
    for rules in [forward, reverse] {
        let set = RuleSet::new(ParameterKey::VariantRate, rules).expect("base exists");
        assert_eq!(
            set.resolve(&context()).expect("unique maximum").id(),
            "combined"
        );
    }
}

#[test]
fn incomparable_maxima_conflict_even_when_values_are_equal() {
    let set = RuleSet::new(
        ParameterKey::VariantRate,
        vec![
            rule("base", SelectorSet::default()),
            rule("population", population()),
            rule("phase", phase()),
        ],
    )
    .expect("base exists");
    let error = set
        .resolve(&context())
        .expect_err("incomparable maxima must be rejected");

    assert_eq!(error.code(), "GPS-CONF-009");
    assert_eq!(
        error.context().get("rules").map(String::as_str),
        Some("phase,population")
    );
}

#[test]
fn declared_individual_region_and_region_subpopulation_intersections_are_preflighted() {
    let base = rule("base", SelectorSet::default());
    let individual = rule(
        "individual",
        SelectorSet {
            individual: Some("founder".to_owned()),
            ..SelectorSet::default()
        },
    );
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
    let subpopulation = rule(
        "subpopulation",
        SelectorSet {
            population: Some(PopulationSelector::Subpopulation("sub".to_owned())),
            ..SelectorSet::default()
        },
    );

    for rules in [
        vec![base.clone(), individual, region.clone()],
        vec![base, region, subpopulation],
    ] {
        let error = RuleSet::new(ParameterKey::VariantRate, rules)
            .expect("base exists")
            .validate_declared_intersections()
            .expect_err("el cruce declarado debe fallar antes de simular");
        assert_eq!(error.code(), "GPS-CONF-009");
        assert_eq!(
            error.context().get("parameter").map(String::as_str),
            Some("variant_rate")
        );
        assert!(error.context().contains_key("scope"));
        assert_eq!(
            error.context()["rules"].split(',').count(),
            2,
            "debe enumerar todas las máximas"
        );
    }
}
