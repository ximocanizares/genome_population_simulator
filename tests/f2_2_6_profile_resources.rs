use std::process::Command;

#[test]
fn independent_oracle_validates_f2_2_6_versioned_profile_resources() {
    let status = Command::new("python3")
        .arg("scripts/verify_f2_2_6_profiles.py")
        .status()
        .expect("python3 must be available for the independent resource oracle");
    assert!(status.success());
}
