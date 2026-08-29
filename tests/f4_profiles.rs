use sha2::{Digest, Sha256};

use genome_population_simulator::technical_profile::{
    F4TechnicalProfileError, GPS_HISEQX_PCRFREE_PE150_PROFILE_VERSION,
    GPS_HISEQX_PCRFREE_PE150_RESOURCE, GPS_HISEQX_PCRFREE_PE150_SHA256,
    GPS_HISEQX_PCRFREE_PE150_V1, LibraryPreparation, PairOrientation, ReadLayout,
    resolve_f4_technical_profile, verify_initial_profile_resource,
};

const DTPD_106_PROFILE_ID: &str = "gps-hiseqx-pcrfree-pe150-v1";
const DTPD_106_PROFILE_SHA256: &str =
    "3d8878ab920237dae03705276a26dbb0823941830d531374d35aa1c540f6deea";

#[test]
fn public_profile_constants_match_the_independent_dtpd_106_literals() {
    assert_eq!(GPS_HISEQX_PCRFREE_PE150_V1, DTPD_106_PROFILE_ID);
    assert_eq!(GPS_HISEQX_PCRFREE_PE150_SHA256, DTPD_106_PROFILE_SHA256);
}

#[test]
fn independent_resource_oracle_confirms_dtpd_106_literals_and_digest() {
    let bytes = include_bytes!("../resources/f4/gps-hiseqx-pcrfree-pe150-v1.json");
    assert_eq!(
        format!("{:x}", Sha256::digest(bytes)),
        DTPD_106_PROFILE_SHA256
    );

    let resource: serde_json::Value = serde_json::from_slice(bytes).expect("profile JSON");
    assert_eq!(resource["schema_version"], "1");
    assert_eq!(resource["profile_id"], DTPD_106_PROFILE_ID);
    assert_eq!(
        resource["classification"],
        "approximate_synthetic_not_empirical_calibration"
    );
    assert_eq!(
        resource["external_reference"],
        serde_json::json!({
            "name": "ART HSXn",
            "use": "comparative_reference_only_not_incorporated_or_redistributed",
        })
    );
    assert_eq!(resource["read"]["layout"], "paired_end");
    assert_eq!(resource["read"]["orientation"], "FR");
    assert_eq!(resource["read"]["minimum_supported_length_bases"], 1);
    assert_eq!(resource["read"]["maximum_supported_length_bases"], 150);
    assert_eq!(
        resource["insert"]["distribution"],
        "normal_truncated_box_muller_v1"
    );
    assert_eq!(resource["insert"]["mean_bases"], 350);
    assert_eq!(resource["insert"]["standard_deviation_bases"], 20);
    assert_eq!(resource["insert"]["minimum_bases"], 250);
    assert_eq!(resource["insert"]["maximum_bases"], 450);
    assert_eq!(
        resource["insert"]["integer_rounding"],
        "nearest_ties_away_from_zero"
    );
    assert_eq!(resource["insert"]["uniform"], "open_53_bit_v1");
    assert_eq!(resource["insert"]["no_cached_normal_variate"], true);
    assert_eq!(
        resource["quality"],
        serde_json::json!({
            "curve_kind": "piecewise_linear_phred_integer_round_half_up_v1",
            "r1": [
                {"first_cycle": 1, "last_cycle": 10, "first_phred": 36, "last_phred": 36},
                {"first_cycle": 11, "last_cycle": 100, "first_phred": 36, "last_phred": 32},
                {"first_cycle": 101, "last_cycle": 150, "first_phred": 32, "last_phred": 28},
            ],
            "r2": [
                {"first_cycle": 1, "last_cycle": 10, "first_phred": 35, "last_phred": 35},
                {"first_cycle": 11, "last_cycle": 100, "first_phred": 35, "last_phred": 30},
                {"first_cycle": 101, "last_cycle": 150, "first_phred": 30, "last_phred": 25},
            ],
        })
    );
    assert_eq!(
        resource["error"],
        serde_json::json!({
            "probability_from_phred": "10^(-Q/10)",
            "alternate_base": "uniform_among_other_acgt",
            "substitution_matrix": "deferred_dtf_053",
        })
    );
    assert_eq!(
        resource["reproducibility"]["rng_catalog"],
        "f4-technical-v3"
    );
    assert_eq!(
        resource["reproducibility"]["profile_algorithm_version"],
        "1"
    );
}

#[test]
fn only_the_approved_profile_resolves_to_typed_f4_facts() {
    verify_initial_profile_resource().expect("embedded resource validates");
    let profile = resolve_f4_technical_profile(DTPD_106_PROFILE_ID).expect("approved profile");

    assert_eq!(profile.id(), DTPD_106_PROFILE_ID);
    assert_eq!(profile.version(), GPS_HISEQX_PCRFREE_PE150_PROFILE_VERSION);
    assert_eq!(profile.resource_path(), GPS_HISEQX_PCRFREE_PE150_RESOURCE);
    assert_eq!(profile.resource_sha256(), DTPD_106_PROFILE_SHA256);
    assert_eq!(profile.preparation(), LibraryPreparation::PcrFree);
    assert_eq!(profile.layout(), ReadLayout::PairedEnd);
    assert_eq!(profile.orientation(), PairOrientation::Fr);
    assert_eq!(profile.common_read_length().minimum_bases(), 1);
    assert_eq!(profile.common_read_length().maximum_bases(), 150);
    assert!(profile.common_read_length().contains(1));
    assert!(profile.common_read_length().contains(150));
    assert!(!profile.common_read_length().contains(0));
    assert!(!profile.common_read_length().contains(151));
    assert_eq!(
        profile.insert_length().algorithm(),
        "normal_truncated_box_muller_v1"
    );
    assert_eq!(profile.insert_length().mean_bases(), 350);
    assert_eq!(profile.insert_length().standard_deviation_bases(), 20);
    assert_eq!(profile.insert_length().minimum_bases(), 250);
    assert_eq!(profile.insert_length().maximum_bases(), 450);
}

#[test]
fn single_end_long_read_and_unknown_profiles_are_rejected_without_output() {
    for profile_id in [
        "illumina-se-minimal-v1",
        "pacbio-hifi-v1",
        "unknown-profile-v1",
    ] {
        assert_eq!(
            resolve_f4_technical_profile(profile_id),
            Err(F4TechnicalProfileError::UnknownProfile {
                requested: profile_id.to_owned(),
            })
        );
    }
}
