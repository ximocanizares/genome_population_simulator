use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{
    build_effective_config_with_entropy, load_and_validate_config_with_entropy,
};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::randomness::{
    EntropySource, RootSeed, SeedOrigin, resolve_root_seed,
};

#[derive(Debug)]
struct FixedEntropy(Result<[u8; 32], &'static str>);

impl EntropySource for FixedEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        match self.0 {
            Ok(bytes) => {
                *destination = bytes;
                Ok(())
            }
            Err(detail) => Err(detail.to_owned()),
        }
    }
}

#[test]
fn root_seed_accepts_exactly_64_lowercase_hex_characters_and_round_trips_all_bytes() {
    let valid = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
    let seed = RootSeed::parse("seed", valid).expect("64 hexadecimales minúsculos");
    assert_eq!(
        seed.as_bytes(),
        &core::array::from_fn::<_, 32, _>(|index| index as u8)
    );
    assert_eq!(seed.to_hex(), valid);
    assert_eq!(seed.to_string(), valid);

    for invalid in [
        &valid[..63],
        &format!("{valid}0"),
        "000102030405060708090A0B0C0D0E0F101112131415161718191a1b1c1d1e1f",
        "g00102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
    ] {
        let error = RootSeed::parse("seed", invalid).expect_err("semilla no canónica");
        assert_eq!(error.code(), "GPS-CONF-017");
        assert_eq!(error.exit_code(), 2);
        assert_eq!(
            error.context().get("field").map(String::as_str),
            Some("seed")
        );
        assert_eq!(
            error.context().get("found").map(String::as_str),
            Some(invalid)
        );
        assert_eq!(
            error.context().get("expected").map(String::as_str),
            Some("64 hexadecimal lowercase characters")
        );
    }

    let high_bytes = RootSeed::from_bytes([0xff; 32]);
    assert_eq!(high_bytes.to_hex(), "ff".repeat(32));
}

#[test]
fn missing_seed_uses_injected_csprng_before_build_and_records_value_and_origin() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/config/general-minimal.toml");
    let expected = [0xa5; 32];
    let effective = load_and_validate_config_with_entropy(&path, &mut FixedEntropy(Ok(expected)))
        .expect("la entropía inyectada debe resolver la raíz");
    assert_eq!(
        effective.reproducibility().root_seed().as_bytes(),
        &expected
    );
    assert_eq!(
        effective.reproducibility().seed_origin(),
        SeedOrigin::SystemCsprng
    );

    let error = load_and_validate_config_with_entropy(
        &path,
        &mut FixedEntropy(Err("system entropy unavailable")),
    )
    .expect_err("no puede construirse configuración efectiva sin raíz");
    assert_eq!(error.code(), "GPS-RES-001");
    assert_eq!(error.exit_code(), 3);
    assert_eq!(error.category(), "RES");
    assert_eq!(
        error.context().get("operation").map(String::as_str),
        Some("generate_root_seed")
    );
}

#[test]
fn declared_seed_bypasses_entropy_and_enters_immutable_effective_config() {
    let bytes = include_bytes!("fixtures/config/general-seeded.toml").to_vec();
    let admitted = admit_toml_header(
        Path::new("general-seeded.toml"),
        bytes,
        ExpectedConfigKind::General,
    )
    .expect("cabecera válida");
    let general = parse_declared_config(admitted).expect("configuración declarada válida");
    let effective =
        build_effective_config_with_entropy(&general, &mut FixedEntropy(Err("must not be called")))
            .expect("la raíz declarada no consume entropía");
    assert_eq!(
        effective.reproducibility().seed_origin(),
        SeedOrigin::Declared
    );
    assert_eq!(
        effective.reproducibility().root_seed().to_hex(),
        "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
    );
    assert_eq!(effective.reproducibility().random_algorithm(), "ChaCha12");
    assert_eq!(effective.reproducibility().derivation_versions(), &[1, 2]);
    assert_eq!(effective.reproducibility().catalog_version(), "6");
}

#[test]
fn resolve_root_seed_reports_the_declared_origin_without_using_entropy() {
    let declared = "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff";
    let (seed, origin) =
        resolve_root_seed(Some(declared), &mut FixedEntropy(Err("must not be called")))
            .expect("raíz declarada");
    assert_eq!(origin, SeedOrigin::Declared);
    assert_eq!(seed.to_hex(), declared);
}
