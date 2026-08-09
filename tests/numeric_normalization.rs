use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::configuration::{DeclaredConfig, ExpectedConfigKind};
use genome_population_simulator::numeric::{
    CANONICAL_ORDER, DEFAULT_NORMALIZATION_TOLERANCE, NORMALIZATION_ALGORITHM,
    NORMALIZATION_VERSION, neumaier, normalize_declared, roundtrip_decimal,
};

fn config(source: &[u8]) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    let admitted = admit_toml_header(
        Path::new("numeric.toml"),
        source.to_vec(),
        ExpectedConfigKind::General,
    )
    .expect("valid header");
    let DeclaredConfig::General(config) = parse_declared_config(admitted).expect("declared config")
    else {
        panic!("expected general")
    };
    *config
}

#[test]
fn toml_decimals_recover_independently_known_f64_bits() {
    let config = config(include_bytes!("fixtures/config/numeric-vectors.toml"));
    let set = &config.frequency_sets()[0];
    let normalized = normalize_declared("edge", &set.values, set.tolerance.as_ref())
        .expect("edge vector is accepted");

    assert_eq!(normalized.values[0].before.to_bits(), 0x3fe0_0000_0000_0000);
    assert_eq!(normalized.values[1].before.to_bits(), 0x3fe0_0000_0000_0002);
    assert_eq!(normalized.tolerance.to_bits(), f64::EPSILON.to_bits());
    assert_eq!(normalized.delta_before.to_bits(), f64::EPSILON.to_bits());
}

#[test]
fn inclusive_edge_inside_and_outside_are_distinct() {
    let edge = config(include_bytes!("fixtures/config/numeric-vectors.toml"));
    let set = &edge.frequency_sets()[0];
    normalize_declared("edge", &set.values, set.tolerance.as_ref()).expect("delta == tolerance");

    let inside = include_bytes!("fixtures/config/numeric-vectors.toml")
        .windows(b"2.220446049250313e-16".len())
        .position(|window| window == b"2.220446049250313e-16")
        .expect("tolerance literal exists");
    let mut inside_source = include_bytes!("fixtures/config/numeric-vectors.toml").to_vec();
    inside_source.splice(
        inside..inside + b"2.220446049250313e-16".len(),
        b"4.440892098500626e-16".iter().copied(),
    );
    let inside_config = config(&inside_source);
    let inside_set = &inside_config.frequency_sets()[0];
    normalize_declared("inside", &inside_set.values, inside_set.tolerance.as_ref())
        .expect("delta < tolerance");

    let mut outside_source = include_bytes!("fixtures/config/numeric-vectors.toml").to_vec();
    outside_source.splice(
        inside..inside + b"2.220446049250313e-16".len(),
        b"1.1102230246251565e-16".iter().copied(),
    );
    let outside_config = config(&outside_source);
    let outside_set = &outside_config.frequency_sets()[0];
    let error = normalize_declared(
        "outside",
        &outside_set.values,
        outside_set.tolerance.as_ref(),
    )
    .expect_err("delta > tolerance");
    assert_eq!(error.code(), "GPS-CONF-016");
}

#[test]
fn canonical_order_and_neumaier_are_stable_under_input_permutation() {
    let first = config(include_bytes!("fixtures/config/numeric-vectors.toml"));
    let mut permuted = first.frequency_sets()[0].values.clone();
    permuted.reverse();
    let first_result = normalize_declared("edge", &first.frequency_sets()[0].values, None)
        .expect("default tolerance");
    let second_result = normalize_declared("edge", &permuted, None).expect("permuted input");

    assert_eq!(first_result, second_result);
    assert_eq!(first_result.canonical_order, CANONICAL_ORDER);
    assert_eq!(first_result.algorithm, NORMALIZATION_ALGORITHM);
    assert_eq!(first_result.algorithm_version, NORMALIZATION_VERSION);
    assert_eq!(first_result.tolerance, DEFAULT_NORMALIZATION_TOLERANCE);
    assert_eq!(
        neumaier([1.0, 1e-16, 1e-16]).to_bits(),
        1.0000000000000002_f64.to_bits()
    );
}

#[test]
fn invalid_values_tolerances_and_zero_sum_are_rejected() {
    let template = include_bytes!("fixtures/config/numeric-vectors.toml");
    for (from, to, expected_code) in [
        ("0.5000000000000002", "nan", "GPS-CONF-015"),
        ("0.5000000000000002", "inf", "GPS-CONF-015"),
        ("0.5000000000000002", "1.1", "GPS-CONF-015"),
        ("0.5000000000000002", "-0.1", "GPS-CONF-015"),
        ("2.220446049250313e-16", "nan", "GPS-CONF-015"),
        ("2.220446049250313e-16", "-1e-9", "GPS-CONF-015"),
        ("2.220446049250313e-16", "1.1e-6", "GPS-CONF-015"),
    ] {
        let source = String::from_utf8(template.to_vec())
            .expect("UTF-8 fixture")
            .replacen(from, to, 1);
        let config = config(source.as_bytes());
        let set = &config.frequency_sets()[0];
        assert_eq!(
            normalize_declared("invalid", &set.values, set.tolerance.as_ref())
                .expect_err("numeric input must fail")
                .code(),
            expected_code
        );
    }

    let zero_source = String::from_utf8(template.to_vec())
        .expect("UTF-8 fixture")
        .replace("0.5000000000000002", "0.0")
        .replacen("value = 0.5", "value = 0.0", 1);
    let zero = config(zero_source.as_bytes());
    let set = &zero.frequency_sets()[0];
    assert_eq!(
        normalize_declared("zero", &set.values, set.tolerance.as_ref())
            .expect_err("zero sum cannot normalize")
            .code(),
        "GPS-CONF-016"
    );
}

#[test]
fn decimal_rendering_round_trips_every_recorded_f64_exactly() {
    let config = config(include_bytes!("fixtures/config/numeric-vectors.toml"));
    let set = &config.frequency_sets()[0];
    let normalized =
        normalize_declared("edge", &set.values, set.tolerance.as_ref()).expect("valid vector");
    let values = [
        normalized.tolerance,
        normalized.sum_before,
        normalized.delta_before,
        normalized.sum_after,
        normalized.delta_after,
        normalized.values[0].before,
        normalized.values[0].after,
        normalized.values[1].before,
        normalized.values[1].after,
    ];
    for value in values {
        let reparsed: f64 = roundtrip_decimal(value).parse().expect("valid decimal");
        assert_eq!(reparsed.to_bits(), value.to_bits());
    }
}
