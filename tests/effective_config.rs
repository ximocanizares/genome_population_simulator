use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{
    DomainLocus, PhaseGenerationInterval, PopulationMembership, ResolutionDomain,
    ResolutionRequest, build_effective_config, build_effective_config_with_domain,
};
use genome_population_simulator::configuration::{DeclaredConfig, ExpectedConfigKind};
use genome_population_simulator::rules::{
    PARAMETER_CANONICAL_ORDER, PARAMETER_ORDER_VERSION, ParameterKey, ResolutionContext, RuleValue,
};

fn declared(
    bytes: Vec<u8>,
) -> Box<genome_population_simulator::configuration::DeclaredGeneralConfig> {
    let admitted = admit_toml_header(
        Path::new("author/general-effective.toml"),
        bytes,
        ExpectedConfigKind::General,
    )
    .expect("la cabecera debe admitirse");
    let DeclaredConfig::General(general) =
        parse_declared_config(admitted).expect("la semántica debe admitirse")
    else {
        panic!("se esperaba configuración general");
    };
    general
}

fn declared_scenario(
    path: &str,
    bytes: Vec<u8>,
) -> Box<genome_population_simulator::configuration::DeclaredScenarioConfig> {
    let admitted = admit_toml_header(Path::new(path), bytes, ExpectedConfigKind::Scenario)
        .expect("la cabecera de escenario debe admitirse");
    let DeclaredConfig::Scenario(scenario) =
        parse_declared_config(admitted).expect("la semántica de escenario debe admitirse")
    else {
        panic!("se esperaba configuración de escenario");
    };
    scenario
}

#[test]
fn effective_snapshot_owns_normalized_values_winner_selectors_and_provenance() {
    let bytes = include_bytes!("fixtures/config/general-effective.toml").to_vec();
    let general = declared(bytes);
    let effective = build_effective_config(
        &general,
        &[],
        &[ResolutionRequest {
            object_id: "ind-1".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                individual: Some("ind-1".to_owned()),
                ..ResolutionContext::default()
            },
        }],
    )
    .expect("la configuración debe quedar completamente resuelta");

    assert_eq!(effective.schema_version(), "1");
    assert_eq!(effective.reference_id(), "ref-main");
    assert_eq!(effective.total_genome_size(), 1000);
    assert_eq!(effective.sequences()[0].id(), "chr1");
    assert_eq!(effective.sequences()[0].ploidy(), 4);
    assert!(effective.explicit_event_ids().is_empty());
    assert_eq!(effective.normalized_collections().len(), 2);
    assert_eq!(effective.controlled_sets()[0].id(), "founders");
    assert_eq!(
        effective.controlled_sets()[0].counts()[0].identity(),
        ("north", "population")
    );
    assert_eq!(effective.controlled_sets()[0].counts()[0].count(), 2);

    let resolved = effective
        .values()
        .iter()
        .find(|value| value.object_id() == "ind-1")
        .expect("valor efectivo individual");
    assert_eq!(resolved.object_id(), "ind-1");
    assert_eq!(resolved.parameter(), ParameterKey::VariantRate);
    assert_eq!(
        resolved.resolution_context().individual.as_deref(),
        Some("ind-1")
    );
    assert_eq!(resolved.value(), &RuleValue::Float(0.02));
    assert_eq!(resolved.winning_rule_id(), "variant-ind-1");
    assert_eq!(resolved.selectors().individual.as_deref(), Some("ind-1"));
    assert_eq!(resolved.provenance().path, "author/general-effective.toml");
    assert!(resolved.provenance().location.contains("line 28"));
    let ploidy = effective
        .values()
        .iter()
        .find(|value| value.parameter() == ParameterKey::Ploidy)
        .expect("sustitución de ploidía efectiva");
    assert_eq!(ploidy.winning_rule_id(), "ploidy-chr1");
    assert_eq!(ploidy.value(), &RuleValue::Integer(4));
}

fn combined_rule_config(rule_blocks: &[&str]) -> Vec<u8> {
    config_with_entities_and_rules(
        r#"[[individuals]]
id = "ind-1""#,
        rule_blocks,
    )
}

fn config_with_entities_and_rules(entities: &str, rule_blocks: &[&str]) -> Vec<u8> {
    format!(
        r#"config_kind = "general"
schema_version = "1"
route = "direct"

[reference]
kind = "simulated"
id = "ref-main"
default_ploidy = 2
composition = [
  {{ id = "A", category = "base", value = 0.5 }},
  {{ id = "C", category = "base", value = 0.5 }},
]

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1000

{entities}

{}
"#,
        rule_blocks.join("\n\n")
    )
    .into_bytes()
}

const BASE_RULE: &str = r#"[[rules]]
id = "base"
parameter = "variant_rate"
value = 0.01"#;
const INDIVIDUAL_RULE: &str = r#"[[rules]]
id = "individual"
parameter = "variant_rate"
value = 0.02
individual = "ind-1""#;
const SEQUENCE_RULE: &str = r#"[[rules]]
id = "sequence"
parameter = "variant_rate"
value = 0.03
reference = "ref-main"
sequence = "chr1""#;
const COMBINED_RULE: &str = r#"[[rules]]
id = "combined"
parameter = "variant_rate"
value = 0.04
individual = "ind-1"
reference = "ref-main"
sequence = "chr1""#;

fn combined_request(object_id: &str) -> ResolutionRequest {
    ResolutionRequest {
        object_id: object_id.to_owned(),
        parameter: ParameterKey::VariantRate,
        context: ResolutionContext {
            individual: Some("ind-1".to_owned()),
            reference: Some("ref-main".to_owned()),
            sequence: Some("chr1".to_owned()),
            ..ResolutionContext::default()
        },
    }
}

#[test]
fn request_context_drives_combined_dominance_independently_of_rule_order() {
    for rules in [
        [BASE_RULE, INDIVIDUAL_RULE, SEQUENCE_RULE, COMBINED_RULE],
        [COMBINED_RULE, SEQUENCE_RULE, INDIVIDUAL_RULE, BASE_RULE],
    ] {
        let effective = build_effective_config(
            &declared(combined_rule_config(&rules)),
            &[],
            &[combined_request("ind-1")],
        )
        .expect("la regla combinada debe ser el máximo único");
        let value = effective
            .values()
            .iter()
            .find(|value| value.object_id() == "ind-1")
            .expect("valor combinado");
        assert_eq!(value.winning_rule_id(), "combined");
        assert_eq!(value.value(), &RuleValue::Float(0.04));
        assert_eq!(
            value.resolution_context(),
            &combined_request("ind-1").context
        );
    }
}

#[test]
fn foreign_object_and_incomparable_maxima_are_rejected_before_snapshot() {
    let complete = declared(combined_rule_config(&[
        BASE_RULE,
        INDIVIDUAL_RULE,
        SEQUENCE_RULE,
        COMBINED_RULE,
    ]));
    let foreign = build_effective_config(&complete, &[], &[combined_request("zzz-foreign")])
        .expect_err("un objeto ajeno no debe llegar a la instantánea");
    assert_eq!(foreign.code(), "GPS-CONF-008");

    let ambiguous = declared(combined_rule_config(&[
        BASE_RULE,
        INDIVIDUAL_RULE,
        SEQUENCE_RULE,
    ]));
    let conflict = build_effective_config(&ambiguous, &[], &[combined_request("ind-1")])
        .expect_err("dos máximas incomparables deben impedir la instantánea");
    assert_eq!(conflict.code(), "GPS-CONF-009");
}

#[test]
fn changing_or_reordering_input_after_admission_cannot_change_effective_snapshot() {
    let mut external = include_bytes!("fixtures/config/general-effective.toml").to_vec();
    let general = declared(external.clone());
    let effective = build_effective_config(&general, &[], &[])
        .expect("la regla base global debe resolverse automáticamente");
    let snapshot = effective.clone();

    external.reverse();
    external.fill(b'x');

    assert_eq!(effective, snapshot);
    let variant = effective
        .values()
        .iter()
        .find(|value| value.parameter() == ParameterKey::VariantRate)
        .expect("valor base global");
    assert_eq!(variant.winning_rule_id(), "variant-base");
    assert_eq!(variant.value(), &RuleValue::Float(0.01));
    assert_eq!(general.source().bytes()[0], b'c');
}

#[test]
fn changing_a_toml_parameter_changes_the_effective_value_without_code_changes() {
    let original = include_bytes!("fixtures/config/general-effective.toml").to_vec();
    let changed = String::from_utf8(original.clone())
        .expect("fixture UTF-8")
        .replacen("value = 0.01", "value = 0.03", 1)
        .into_bytes();

    let original = build_effective_config(&declared(original), &[], &[])
        .expect("la configuración original debe ser válida");
    let changed = build_effective_config(&declared(changed), &[], &[])
        .expect("la configuración editada debe ser válida");

    let base_value = |effective: &genome_population_simulator::application::EffectiveConfig| {
        effective
            .values()
            .iter()
            .find(|value| value.winning_rule_id() == "variant-base")
            .expect("valor global de variantes")
            .value()
            .clone()
    };
    assert_eq!(base_value(&original), RuleValue::Float(0.01));
    assert_eq!(base_value(&changed), RuleValue::Float(0.03));
}

#[test]
fn parameter_axes_reject_biological_rule_in_a_technical_event_context() {
    let config = declared(config_with_entities_and_rules(
        r#"[[events]]
id = "event-1""#,
        &[BASE_RULE],
    ));
    let error = build_effective_config(
        &config,
        &[],
        &[ResolutionRequest {
            object_id: "event-1".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                technical_event: Some("event-1".to_owned()),
                ..ResolutionContext::default()
            },
        }],
    )
    .expect_err("variant_rate no admite el eje de evento técnico");
    assert_eq!(error.code(), "GPS-CONF-012");
}

#[test]
fn complete_population_and_temporal_hierarchies_reach_ruleset_dominance() {
    let population_domain = ResolutionDomain::validated(
        "synthetic-precedence-test",
        [PopulationMembership {
            population: "pop-1".to_owned(),
            subpopulation: "sub-1".to_owned(),
        }],
        [],
        [],
    )
    .expect("dominio sintético explícito");
    let population_config = declared(config_with_entities_and_rules(
        r#"[[populations]]
id = "pop-1"
[[subpopulations]]
id = "sub-1""#,
        &[
            BASE_RULE,
            r#"[[rules]]
id = "population"
parameter = "variant_rate"
value = 0.02
population = "pop-1""#,
            r#"[[rules]]
id = "subpopulation"
parameter = "variant_rate"
value = 0.03
subpopulation = "sub-1""#,
        ],
    ));
    let population = build_effective_config_with_domain(
        &population_config,
        &[],
        &[ResolutionRequest {
            object_id: "sub-1".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                population: Some("pop-1".to_owned()),
                subpopulation: Some("sub-1".to_owned()),
                ..ResolutionContext::default()
            },
        }],
        &population_domain,
    )
    .expect("subpoblación debe dominar a población");
    assert_eq!(population.values()[0].winning_rule_id(), "subpopulation");

    let temporal_config = declared(config_with_entities_and_rules(
        r#"[[phases]]
id = "phase-1""#,
        &[
            BASE_RULE,
            r#"[[rules]]
id = "phase"
parameter = "variant_rate"
value = 0.02
phase = "phase-1""#,
            r#"[[rules]]
id = "generation"
parameter = "variant_rate"
value = 0.03
generation = 7"#,
        ],
    ));
    let temporal_domain = ResolutionDomain::validated(
        "synthetic-precedence-test",
        [],
        [PhaseGenerationInterval {
            phase: "phase-1".to_owned(),
            start: 1,
            end: 10,
        }],
        [],
    )
    .expect("intervalo sintético explícito");
    let temporal = build_effective_config_with_domain(
        &temporal_config,
        &[],
        &[ResolutionRequest {
            object_id: "phase-1".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                phase: Some("phase-1".to_owned()),
                generation: Some(7),
                ..ResolutionContext::default()
            },
        }],
        &temporal_domain,
    )
    .expect("generación debe dominar a fase");
    assert_eq!(temporal.values()[0].winning_rule_id(), "generation");
}

#[test]
fn complete_locus_position_hierarchy_reaches_ruleset_dominance() {
    let domain = ResolutionDomain::validated(
        "synthetic-precedence-test",
        [],
        [],
        [DomainLocus {
            id: "locus-15".to_owned(),
            reference: "ref-main".to_owned(),
            sequence: "chr1".to_owned(),
            position: 15,
        }],
    )
    .expect("catálogo sintético de loci");
    let config = declared(config_with_entities_and_rules(
        "",
        &[
            BASE_RULE,
            r#"[[rules]]
id = "sequence"
parameter = "variant_rate"
value = 0.02
reference = "ref-main"
sequence = "chr1""#,
            r#"[[rules]]
id = "region"
parameter = "variant_rate"
value = 0.03
reference = "ref-main"
sequence = "chr1"
region_start = 10
region_end = 20"#,
            r#"[[rules]]
id = "locus"
parameter = "variant_rate"
value = 0.04
reference = "ref-main"
sequence = "chr1"
locus = "locus-15""#,
        ],
    ));
    let effective = build_effective_config_with_domain(
        &config,
        &[],
        &[ResolutionRequest {
            object_id: "locus-15".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                position: Some(15),
                locus: Some("locus-15".to_owned()),
                ..ResolutionContext::default()
            },
        }],
        &domain,
    )
    .expect("locus debe dominar a región y secuencia");
    assert_eq!(effective.values()[0].winning_rule_id(), "locus");
}

#[test]
fn orphan_and_incoherent_spatial_coordinates_are_rejected() {
    let config = declared(include_bytes!("fixtures/config/general-effective.toml").to_vec());
    for (label, context) in [
        (
            "orphan-position",
            ResolutionContext {
                individual: Some("ind-1".to_owned()),
                position: Some(10),
                ..ResolutionContext::default()
            },
        ),
        (
            "locus-without-position",
            ResolutionContext {
                individual: Some("ind-1".to_owned()),
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                locus: Some("locus-15".to_owned()),
                ..ResolutionContext::default()
            },
        ),
        (
            "position-out-of-bounds",
            ResolutionContext {
                individual: Some("ind-1".to_owned()),
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                position: Some(1001),
                ..ResolutionContext::default()
            },
        ),
    ] {
        let error = build_effective_config(
            &config,
            &[],
            &[ResolutionRequest {
                object_id: "ind-1".to_owned(),
                parameter: ParameterKey::VariantRate,
                context,
            }],
        )
        .unwrap_err();
        assert_eq!(error.code(), "GPS-CONF-008", "case={label}");
    }
}

#[test]
fn spatial_level_and_primary_object_follow_parameter_metadata() {
    for (label, parameter, rule, object_id, context, expected_code) in [
        (
            "coverage-locus",
            ParameterKey::Coverage,
            r#"[[rules]]
id = "coverage-base"
parameter = "coverage"
value = 30.0"#,
            "locus-15",
            ResolutionContext {
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                position: Some(15),
                locus: Some("locus-15".to_owned()),
                ..ResolutionContext::default()
            },
            "GPS-CONF-012",
        ),
        (
            "recombination-region",
            ParameterKey::RecombinationRate,
            r#"[[rules]]
id = "recombination-base"
parameter = "recombination_rate"
value = 0.01"#,
            "chr1",
            ResolutionContext {
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                position: Some(15),
                ..ResolutionContext::default()
            },
            "GPS-CONF-012",
        ),
        (
            "reference-is-not-the-spatial-object",
            ParameterKey::VariantRate,
            BASE_RULE,
            "ref-main",
            ResolutionContext {
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                ..ResolutionContext::default()
            },
            "GPS-CONF-008",
        ),
    ] {
        let error = build_effective_config(
            &declared(config_with_entities_and_rules("", &[rule])),
            &[],
            &[ResolutionRequest {
                object_id: object_id.to_owned(),
                parameter,
                context,
            }],
        )
        .unwrap_err();
        assert_eq!(error.code(), expected_code, "case={label}");
    }
}

#[test]
fn hierarchical_contexts_are_not_inferred_from_flat_identifier_lists() {
    let config = declared(config_with_entities_and_rules(
        r#"[[populations]]
id = "pop-1"
[[subpopulations]]
id = "sub-1""#,
        &[
            BASE_RULE,
            r#"[[rules]]
id = "subpopulation"
parameter = "variant_rate"
value = 0.03
subpopulation = "sub-1""#,
        ],
    ));
    let error = build_effective_config(
        &config,
        &[],
        &[ResolutionRequest {
            object_id: "sub-1".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                population: Some("pop-1".to_owned()),
                subpopulation: Some("sub-1".to_owned()),
                ..ResolutionContext::default()
            },
        }],
    )
    .expect_err("las listas planas no prueban la relación jerárquica");
    assert_eq!(error.code(), "GPS-CONF-008");
}

#[test]
fn hierarchical_contexts_must_match_parent_interval_and_coordinate_evidence() {
    let population_config = declared(config_with_entities_and_rules(
        r#"[[populations]]
id = "pop-1"
[[populations]]
id = "pop-2"
[[subpopulations]]
id = "sub-1""#,
        &[BASE_RULE],
    ));
    let wrong_parent = ResolutionDomain::validated(
        "synthetic-negative-domain",
        [PopulationMembership {
            population: "pop-2".to_owned(),
            subpopulation: "sub-1".to_owned(),
        }],
        [],
        [],
    )
    .unwrap();
    let error = build_effective_config_with_domain(
        &population_config,
        &[],
        &[ResolutionRequest {
            object_id: "sub-1".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                population: Some("pop-1".to_owned()),
                subpopulation: Some("sub-1".to_owned()),
                ..ResolutionContext::default()
            },
        }],
        &wrong_parent,
    )
    .expect_err("una relación con otro padre no prueba este contexto");
    assert_eq!(error.code(), "GPS-CONF-008");

    let phase_config = declared(config_with_entities_and_rules(
        r#"[[phases]]
id = "phase-1""#,
        &[BASE_RULE],
    ));
    let wrong_interval = ResolutionDomain::validated(
        "synthetic-negative-domain",
        [],
        [PhaseGenerationInterval {
            phase: "phase-1".to_owned(),
            start: 1,
            end: 4,
        }],
        [],
    )
    .unwrap();
    let error = build_effective_config_with_domain(
        &phase_config,
        &[],
        &[ResolutionRequest {
            object_id: "phase-1".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                phase: Some("phase-1".to_owned()),
                generation: Some(5),
                ..ResolutionContext::default()
            },
        }],
        &wrong_interval,
    )
    .expect_err("una generación fuera del intervalo no está probada");
    assert_eq!(error.code(), "GPS-CONF-008");

    let locus_config = declared(config_with_entities_and_rules("", &[BASE_RULE]));
    let wrong_coordinate = ResolutionDomain::validated(
        "synthetic-negative-domain",
        [],
        [],
        [DomainLocus {
            id: "locus-15".to_owned(),
            reference: "ref-main".to_owned(),
            sequence: "chr1".to_owned(),
            position: 16,
        }],
    )
    .unwrap();
    let error = build_effective_config_with_domain(
        &locus_config,
        &[],
        &[ResolutionRequest {
            object_id: "locus-15".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                position: Some(15),
                locus: Some("locus-15".to_owned()),
                ..ResolutionContext::default()
            },
        }],
        &wrong_coordinate,
    )
    .expect_err("una coordenada distinta no prueba el locus solicitado");
    assert_eq!(error.code(), "GPS-CONF-008");
}

#[test]
fn every_parameter_is_checked_at_sequence_region_and_locus_resolution_levels() {
    let domain = ResolutionDomain::validated(
        "synthetic-spatial-matrix",
        [],
        [],
        [DomainLocus {
            id: "locus-15".to_owned(),
            reference: "ref-main".to_owned(),
            sequence: "chr1".to_owned(),
            position: 15,
        }],
    )
    .expect("catálogo espacial sintético");
    let cases = [
        (
            ParameterKey::Ploidy,
            r#"[[rules]]
id = "ploidy-chr1"
parameter = "ploidy"
value = 2
reference = "ref-main"
sequence = "chr1""#,
            [true, false, false],
        ),
        (ParameterKey::VariantRate, BASE_RULE, [true, true, true]),
        (
            ParameterKey::SelfingRate,
            r#"[[rules]]
id = "selfing-base"
parameter = "selfing_rate"
value = 0.1"#,
            [false, false, false],
        ),
        (
            ParameterKey::RecombinationRate,
            r#"[[rules]]
id = "recombination-base"
parameter = "recombination_rate"
value = 0.01"#,
            [true, false, false],
        ),
        (
            ParameterKey::Coverage,
            r#"[[rules]]
id = "coverage-base"
parameter = "coverage"
value = 30.0"#,
            [true, true, false],
        ),
        (
            ParameterKey::ErrorRate,
            r#"[[rules]]
id = "error-base"
parameter = "error_rate"
value = 0.01"#,
            [true, true, false],
        ),
        (
            ParameterKey::ContaminationRate,
            r#"[[rules]]
id = "contamination-base"
parameter = "contamination_rate"
value = 0.01"#,
            [false, false, false],
        ),
        (
            ParameterKey::OutputEnabled,
            r#"[[rules]]
id = "output-base"
parameter = "output_enabled"
value = true"#,
            [false, false, false],
        ),
    ];
    for (parameter, rule, expected) in cases {
        let config = declared(config_with_entities_and_rules("", &[rule]));
        for (index, context) in [
            ResolutionContext {
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                ..ResolutionContext::default()
            },
            ResolutionContext {
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                position: Some(15),
                ..ResolutionContext::default()
            },
            ResolutionContext {
                reference: Some("ref-main".to_owned()),
                sequence: Some("chr1".to_owned()),
                position: Some(15),
                locus: Some("locus-15".to_owned()),
                ..ResolutionContext::default()
            },
        ]
        .into_iter()
        .enumerate()
        {
            if parameter == ParameterKey::Ploidy && index == 0 {
                let effective = build_effective_config(&config, &[], &[])
                    .expect("la ploidía de secuencia es un positivo automático");
                assert!(effective.values().iter().any(|value| {
                    value.object_id() == "chr1" && value.parameter() == ParameterKey::Ploidy
                }));
                continue;
            }
            let object_id = if index == 2 { "locus-15" } else { "chr1" };
            let result = build_effective_config_with_domain(
                &config,
                &[],
                &[ResolutionRequest {
                    object_id: object_id.to_owned(),
                    parameter,
                    context,
                }],
                &domain,
            );
            assert_eq!(
                result.is_ok(),
                expected[index],
                "parameter={} level_index={index} error={:?}",
                parameter.as_str(),
                result.err().map(|error| error.render_human())
            );
        }
    }
}

#[test]
fn phase_ids_are_globally_unique_across_scenarios_independently_of_order() {
    let mut general =
        String::from_utf8(include_bytes!("fixtures/config/general-evolution.toml").to_vec())
            .expect("fixture general UTF-8");
    general.push_str(
        r#"

[[scenarios]]
id = "evolution-other"
path = "scenario-evolution-other.toml"
route = "evolution"
"#,
    );
    let general = declared(general.into_bytes());
    let first = declared_scenario(
        "author/scenario-evolution.toml",
        include_bytes!("fixtures/config/scenario-evolution.toml").to_vec(),
    );
    let second = declared_scenario(
        "author/scenario-evolution-other.toml",
        String::from_utf8(include_bytes!("fixtures/config/scenario-evolution.toml").to_vec())
            .expect("fixture scenario UTF-8")
            .replace("evolution-small", "evolution-other")
            .replace("sample-generation-10", "sample-generation-20")
            .into_bytes(),
    );

    for scenarios in [
        [first.as_ref(), second.as_ref()],
        [second.as_ref(), first.as_ref()],
    ] {
        let error = build_effective_config(&general, &scenarios, &[])
            .expect_err("phase-1 no puede pertenecer a dos escenarios v1");
        assert_eq!(error.code(), "GPS-CONF-008");
        let rendered = error.render_human();
        assert!(rendered.contains("phase-1"));
        assert!(rendered.contains("evolution-small"));
        assert!(rendered.contains("evolution-other"));
        assert!(rendered.contains("scenario-evolution.toml"));
        assert!(rendered.contains("scenario-evolution-other.toml"));
    }
}

#[test]
fn parameter_order_v1_is_literal_and_independent_of_rule_and_request_permutations() {
    assert_eq!(PARAMETER_ORDER_VERSION, "1");
    assert_eq!(
        PARAMETER_CANONICAL_ORDER.map(ParameterKey::as_str),
        [
            "ploidy",
            "variant_rate",
            "selfing_rate",
            "recombination_rate",
            "coverage",
            "error_rate",
            "contamination_rate",
            "output_enabled",
        ]
    );
    let variant = r#"[[rules]]
id = "variant-ind"
parameter = "variant_rate"
value = 0.02
individual = "ind-1""#;
    let variant_base = BASE_RULE;
    let output_base = r#"[[rules]]
id = "output-base"
parameter = "output_enabled"
value = false"#;
    let output = r#"[[rules]]
id = "output-ind"
parameter = "output_enabled"
value = true
individual = "ind-1""#;
    let requests = [
        ResolutionRequest {
            object_id: "ind-1".to_owned(),
            parameter: ParameterKey::VariantRate,
            context: ResolutionContext {
                individual: Some("ind-1".to_owned()),
                ..ResolutionContext::default()
            },
        },
        ResolutionRequest {
            object_id: "ind-1".to_owned(),
            parameter: ParameterKey::OutputEnabled,
            context: ResolutionContext {
                individual: Some("ind-1".to_owned()),
                ..ResolutionContext::default()
            },
        },
    ];
    for (rules, request_order) in [
        ([variant_base, variant, output_base, output], [0, 1]),
        ([output, output_base, variant, variant_base], [1, 0]),
    ] {
        let effective = build_effective_config(
            &declared(config_with_entities_and_rules(
                "[[individuals]]\nid = \"ind-1\"",
                &rules,
            )),
            &[],
            &[
                requests[request_order[0]].clone(),
                requests[request_order[1]].clone(),
            ],
        )
        .expect("ambos parámetros deben resolverse");
        let observed = effective
            .values()
            .iter()
            .filter(|value| value.object_id() == "ind-1")
            .map(|value| value.parameter().as_str())
            .collect::<Vec<_>>();
        assert_eq!(observed, ["variant_rate", "output_enabled"]);
    }
}
