use std::collections::{BTreeMap, BTreeSet};

use genome_population_simulator::rules::{
    AxisKind, PARAMETER_CATALOG, ParameterKey, PopulationSelector, RegionSelector, Rule,
    RuleProvenance, RuleValue, SelectorSet, SpatialSelector, TemporalSelector, metadata,
    validate_selector_set,
};

#[test]
fn published_catalog_and_code_cover_the_same_parameter_families() {
    let catalog = include_str!("../catalogs/parameter-families-v1.tsv");
    let rows: BTreeMap<_, _> = catalog
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let fields: Vec<_> = line.split('|').collect();
            (fields[0], (fields[3], fields[5]))
        })
        .collect();

    assert_eq!(rows.len(), PARAMETER_CATALOG.len());
    for entry in PARAMETER_CATALOG {
        let axes = entry
            .allowed_axes
            .iter()
            .map(|axis| axis.as_str())
            .collect::<Vec<_>>()
            .join(";");
        assert_eq!(
            rows.get(entry.parameter.as_str()),
            Some(&(axes.as_str(), entry.value_domain.as_str()))
        );
    }
    assert_eq!(
        PARAMETER_CATALOG
            .into_iter()
            .map(|entry| entry.family)
            .collect::<BTreeSet<_>>()
            .len(),
        7
    );
}

#[test]
fn selector_acceptance_table_is_independent_of_resolution() {
    let individual = SelectorSet {
        individual: Some("founder-1".to_owned()),
        ..SelectorSet::default()
    };
    let population = SelectorSet {
        population: Some(PopulationSelector::Subpopulation("sub-1".to_owned())),
        ..SelectorSet::default()
    };
    let temporal = SelectorSet {
        temporal: Some(TemporalSelector::Generation(10)),
        ..SelectorSet::default()
    };
    let locus = SelectorSet {
        spatial: Some(SpatialSelector::Locus {
            reference: "ref".to_owned(),
            sequence: "chr1".to_owned(),
            locus: "locus-1".to_owned(),
        }),
        ..SelectorSet::default()
    };
    let region = SelectorSet {
        spatial: Some(SpatialSelector::Region(RegionSelector {
            reference: "ref".to_owned(),
            sequence: "chr1".to_owned(),
            start: 1,
            end: 10,
        })),
        ..SelectorSet::default()
    };
    let event = SelectorSet {
        technical_event: Some("event-1".to_owned()),
        ..SelectorSet::default()
    };

    let cases = [
        (ParameterKey::VariantRate, &individual, true),
        (ParameterKey::VariantRate, &population, true),
        (ParameterKey::VariantRate, &temporal, true),
        (ParameterKey::VariantRate, &locus, true),
        (ParameterKey::Ploidy, &region, false),
        (ParameterKey::SelfingRate, &locus, false),
        (ParameterKey::Coverage, &region, true),
        (ParameterKey::Coverage, &locus, false),
        (ParameterKey::Coverage, &event, true),
        (ParameterKey::VariantRate, &event, false),
    ];
    for (parameter, selectors, accepted) in cases {
        assert_eq!(
            validate_selector_set(parameter, selectors).is_ok(),
            accepted
        );
    }
}

#[test]
fn prohibited_selector_diagnostic_identifies_parameter_scope_and_domain_boundary() {
    let technical_event = SelectorSet {
        technical_event: Some("event-1".to_owned()),
        ..SelectorSet::default()
    };
    let error = validate_selector_set(ParameterKey::VariantRate, &technical_event)
        .expect_err("un evento técnico no puede seleccionar un parámetro biológico");
    assert_eq!(error.code(), "GPS-CONF-012");
    assert_eq!(error.context()["parameter"], "variant_rate");
    assert_eq!(error.context()["selector"], "technical_event");
    assert_eq!(error.context()["scope"], "declared_rule");

    validate_selector_set(ParameterKey::Coverage, &technical_event)
        .expect("el mismo eje sí está autorizado para un parámetro técnico");
}

#[test]
fn every_catalog_axis_and_spatial_level_has_an_acceptance_result() {
    let axes = [
        AxisKind::Individual,
        AxisKind::Population,
        AxisKind::Temporal,
        AxisKind::Spatial,
        AxisKind::Group,
        AxisKind::Block,
        AxisKind::TechnicalEvent,
    ];
    for entry in PARAMETER_CATALOG {
        for axis in axes {
            let selectors = selectors_for_axis(axis);
            let expected = entry.allowed_axes.contains(&axis)
                && (axis != AxisKind::Spatial
                    || entry
                        .allowed_spatial_levels
                        .contains(&genome_population_simulator::rules::SpatialLevel::Sequence));
            assert_eq!(
                validate_selector_set(entry.parameter, &selectors).is_ok(),
                expected,
                "{} / {}",
                entry.parameter.as_str(),
                axis.as_str()
            );
        }
        for spatial in [
            SpatialSelector::Sequence {
                reference: "ref".to_owned(),
                sequence: "chr1".to_owned(),
            },
            SpatialSelector::Region(RegionSelector {
                reference: "ref".to_owned(),
                sequence: "chr1".to_owned(),
                start: 1,
                end: 2,
            }),
            SpatialSelector::Locus {
                reference: "ref".to_owned(),
                sequence: "chr1".to_owned(),
                locus: "locus".to_owned(),
            },
        ] {
            let level = spatial.level();
            let selectors = SelectorSet {
                spatial: Some(spatial),
                ..SelectorSet::default()
            };
            assert_eq!(
                validate_selector_set(entry.parameter, &selectors).is_ok(),
                entry.allowed_spatial_levels.contains(&level)
            );
        }
    }
}

fn selectors_for_axis(axis: AxisKind) -> SelectorSet {
    match axis {
        AxisKind::Individual => SelectorSet {
            individual: Some("individual".to_owned()),
            ..SelectorSet::default()
        },
        AxisKind::Population => SelectorSet {
            population: Some(PopulationSelector::Population("population".to_owned())),
            ..SelectorSet::default()
        },
        AxisKind::Temporal => SelectorSet {
            temporal: Some(TemporalSelector::Phase("phase".to_owned())),
            ..SelectorSet::default()
        },
        AxisKind::Spatial => SelectorSet {
            spatial: Some(SpatialSelector::Sequence {
                reference: "ref".to_owned(),
                sequence: "chr1".to_owned(),
            }),
            ..SelectorSet::default()
        },
        AxisKind::Group => SelectorSet {
            group: Some("group".to_owned()),
            ..SelectorSet::default()
        },
        AxisKind::Block => SelectorSet {
            block: Some("block".to_owned()),
            ..SelectorSet::default()
        },
        AxisKind::TechnicalEvent => SelectorSet {
            technical_event: Some("event".to_owned()),
            ..SelectorSet::default()
        },
    }
}

#[test]
fn catalog_marks_technical_and_biological_families_explicitly() {
    assert!(metadata(ParameterKey::VariantRate).is_biological);
    assert!(!metadata(ParameterKey::Coverage).is_biological);
    assert_eq!(
        metadata(ParameterKey::Ploidy).allowed_axes,
        &[AxisKind::Spatial]
    );
}

#[test]
fn every_parameter_value_domain_accepts_its_inclusive_boundaries() {
    let provenance = RuleProvenance {
        path: "boundaries.toml".to_owned(),
        location: "line 1, column 1".to_owned(),
    };
    let cases = [
        (ParameterKey::Ploidy, RuleValue::Integer(1)),
        (ParameterKey::Ploidy, RuleValue::Integer(10)),
        (ParameterKey::SelfingRate, RuleValue::Float(0.0)),
        (ParameterKey::SelfingRate, RuleValue::Float(1.0)),
        (ParameterKey::ErrorRate, RuleValue::Float(0.0)),
        (ParameterKey::ErrorRate, RuleValue::Float(1.0)),
        (ParameterKey::ContaminationRate, RuleValue::Float(0.0)),
        (ParameterKey::ContaminationRate, RuleValue::Float(1.0)),
        (ParameterKey::VariantRate, RuleValue::Float(0.0)),
        (ParameterKey::RecombinationRate, RuleValue::Float(0.0)),
        (ParameterKey::Coverage, RuleValue::Float(0.0)),
        (ParameterKey::OutputEnabled, RuleValue::Boolean(false)),
        (ParameterKey::OutputEnabled, RuleValue::Boolean(true)),
    ];

    for (index, (parameter, value)) in cases.into_iter().enumerate() {
        Rule::new(
            format!("boundary-{index}"),
            parameter,
            value,
            SelectorSet::default(),
            provenance.clone(),
        )
        .unwrap_or_else(|error| panic!("{}: {error:?}", parameter.as_str()));
    }
}
