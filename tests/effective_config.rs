use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{ResolutionRequest, build_effective_config};
use genome_population_simulator::configuration::{DeclaredConfig, ExpectedConfigKind};
use genome_population_simulator::rules::{ParameterKey, ResolutionContext, RuleValue};

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
