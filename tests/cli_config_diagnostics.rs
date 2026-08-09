use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const BINARY: &str = env!("CARGO_BIN_EXE_genome-population-simulator");

fn run(path: &Path, mode: &str, output: &Path) -> Output {
    Command::new(BINARY)
        .args(["validate"])
        .arg(path)
        .args(["--diagnostics", mode, "--output"])
        .arg(output)
        .output()
        .expect("el binario CLI debe poder ejecutarse")
}

fn isolated_destination(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("el reloj debe ser posterior a UNIX_EPOCH")
        .as_nanos();
    std::env::temp_dir().join(format!("gps-f1-2-{}-{nonce}-{label}", std::process::id()))
}

fn code_from_human(stderr: &[u8]) -> &str {
    std::str::from_utf8(stderr)
        .expect("diagnóstico humano UTF-8")
        .split_whitespace()
        .next()
        .expect("código humano")
}

fn code_from_json(stdout: &[u8]) -> &str {
    let document = std::str::from_utf8(stdout)
        .expect("diagnóstico JSON UTF-8")
        .trim();
    let prefix = "\"code\":\"";
    let start = document.find(prefix).expect("campo code") + prefix.len();
    let end = document[start..].find('"').expect("fin de code") + start;
    &document[start..end]
}

fn general_with_rule(parameter: &str, value: &str, selectors: &str) -> String {
    format!(
        concat!(
            "config_kind = \"general\"\n",
            "schema_version = \"1\"\n",
            "route = \"direct\"\n",
            "[reference]\n",
            "kind = \"simulated\"\n",
            "id = \"ref-main\"\n",
            "default_ploidy = 2\n",
            "composition = [{{ id = \"A\", category = \"base\", value = 1.0 }}]\n",
            "[[reference.sequences]]\n",
            "id = \"chr1\"\n",
            "class = \"chromosome\"\n",
            "length = 1000\n",
            "[[rules]]\n",
            "id = \"rule-under-test\"\n",
            "parameter = \"{}\"\n",
            "value = {}\n",
            "{}\n"
        ),
        parameter,
        value,
        selectors.replace(';', "\n")
    )
}

#[test]
fn every_invalid_fixture_has_equivalent_human_and_single_json_diagnostics() {
    let cases = [
        ("invalid-missing-kind.toml", "GPS-CONF-002"),
        ("invalid-missing-version.toml", "GPS-CONF-004"),
        ("invalid-unknown-kind.toml", "GPS-CONF-003"),
        ("invalid-crossed-kind.toml", "GPS-CONF-003"),
        ("invalid-version-zero.toml", "GPS-CONF-005"),
        ("invalid-version-two.toml", "GPS-CONF-005"),
        ("invalid-version-integer.toml", "GPS-CONF-005"),
        ("invalid-malformed.toml", "GPS-CONF-001"),
        ("invalid-duplicate-key.toml", "GPS-CONF-007"),
        ("invalid-event-duplicate-id.toml", "GPS-CONF-008"),
        ("invalid-event-missing-id.toml", "GPS-CONF-008"),
        ("invalid-frequency-forms.toml", "GPS-CONF-008"),
        ("invalid-size-overflow.toml", "GPS-CONF-008"),
        ("invalid-supplied-ploidy-missing.toml", "GPS-CONF-014"),
        ("invalid-unknown-field.toml", "GPS-CONF-006"),
    ];
    for (fixture, expected_code) in cases {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/config")
            .join(fixture);
        let destination = isolated_destination(fixture);
        let human = run(&path, "human", &destination);
        let json = run(&path, "json", &destination);

        assert_eq!(human.status.code(), Some(2), "{fixture}");
        assert_eq!(json.status.code(), Some(2), "{fixture}");
        assert!(human.stdout.is_empty(), "{fixture}");
        assert!(json.stderr.is_empty(), "{fixture}");
        assert_eq!(code_from_human(&human.stderr), expected_code, "{fixture}");
        assert_eq!(code_from_json(&json.stdout), expected_code, "{fixture}");
        assert_eq!(
            code_from_human(&human.stderr),
            code_from_json(&json.stdout),
            "{fixture}"
        );

        let document = std::str::from_utf8(&json.stdout).expect("JSON UTF-8");
        assert_eq!(
            document.matches("\"diagnostic_kind\"").count(),
            1,
            "{fixture}"
        );
        assert!(document.starts_with('{'), "{fixture}");
        assert!(document.trim_end().ends_with('}'), "{fixture}");
        assert_eq!(document.trim_end().matches('\n').count(), 0, "{fixture}");
        assert!(document.contains("\"schema_version\":\"1\""), "{fixture}");
        assert!(document.contains("\"exit_code\":2"), "{fixture}");
        assert!(document.contains("\"context\":{"), "{fixture}");
        assert!(
            !destination.exists(),
            "{fixture} publicó un destino inválido"
        );
    }
}

#[test]
fn valid_validation_is_read_only_even_with_an_output_argument() {
    for fixture in [
        "general-minimal.toml",
        "general-supplied.toml",
        "general-evolution.toml",
        "general-families.toml",
        "general-combined-population-phase.toml",
    ] {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/config")
            .join(fixture);
        let destination = isolated_destination(fixture);
        let result = run(&path, "human", &destination);

        assert!(result.status.success(), "{fixture}");
        assert!(result.stdout.is_empty(), "{fixture}");
        assert!(result.stderr.is_empty(), "{fixture}");
        assert!(!destination.exists(), "{fixture}");
        assert!(fs::metadata(destination.parent().expect("parent")).is_ok());
    }
}

#[test]
fn incomparable_declared_intersections_fail_before_simulation_independent_of_rule_order() {
    for fixture in [
        "invalid-incomparable-population-phase.toml",
        "invalid-incomparable-population-phase-reordered.toml",
    ] {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/config")
            .join(fixture);
        let destination = isolated_destination(fixture);
        let human = run(&path, "human", &destination);
        let json = run(&path, "json", &destination);

        assert_eq!(human.status.code(), Some(2));
        assert_eq!(json.status.code(), Some(2));
        assert_eq!(code_from_human(&human.stderr), "GPS-CONF-009");
        assert_eq!(code_from_json(&json.stdout), "GPS-CONF-009");
        for rule in ["selfing-phase", "selfing-population"] {
            assert!(String::from_utf8_lossy(&human.stderr).contains(rule));
            assert!(String::from_utf8_lossy(&json.stdout).contains(rule));
        }
        assert!(String::from_utf8_lossy(&json.stdout).contains("scope"));
        assert!(!destination.exists());
    }
}

#[test]
fn rule_values_reject_wrong_types_non_finite_numbers_and_domain_boundaries_in_both_modes() {
    let cases = include_str!("fixtures/config/invalid-rule-values.tsv");
    for row in cases
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let columns = row.split('|').collect::<Vec<_>>();
        let [case, parameter, value, selectors, expected_code] = columns.as_slice() else {
            panic!("fixture mal formado: {row}");
        };
        let source = general_with_rule(parameter, value, selectors);
        let input = isolated_destination(&format!("{case}.toml"));
        fs::write(&input, source).expect("crear fixture TOML aislado");
        let destination = isolated_destination(&format!("{case}.output"));
        let human = run(&input, "human", &destination);
        let json = run(&input, "json", &destination);

        assert_eq!(human.status.code(), Some(2), "{case}");
        assert_eq!(json.status.code(), Some(2), "{case}");
        assert_eq!(code_from_human(&human.stderr), *expected_code, "{case}");
        assert_eq!(code_from_json(&json.stdout), *expected_code, "{case}");
        assert!(String::from_utf8_lossy(&human.stderr).contains(parameter));
        assert!(String::from_utf8_lossy(&json.stdout).contains(parameter));
        assert!(!destination.exists(), "{case}");
        fs::remove_file(&input).expect("retirar fixture aislado");
    }
}

#[test]
fn ploidy_rejections_cover_absent_non_integer_below_and_above_bounds_with_sequence_identity() {
    let cases = include_str!("fixtures/config/invalid-rule-values.tsv");
    for row in cases.lines().filter(|line| line.starts_with("ploidy_")) {
        let columns = row.split('|').collect::<Vec<_>>();
        let [case, parameter, value, selectors, expected_code] = columns.as_slice() else {
            panic!("fixture mal formado: {row}");
        };
        assert_eq!(*expected_code, "GPS-CONF-014");
        let input = isolated_destination(&format!("{case}.toml"));
        fs::write(&input, general_with_rule(parameter, value, selectors))
            .expect("crear fixture TOML de ploidía");
        let destination = isolated_destination(&format!("{case}.output"));
        for mode in ["human", "json"] {
            let result = run(&input, mode, &destination);
            assert_eq!(result.status.code(), Some(2), "{case}/{mode}");
            let output = if mode == "human" {
                &result.stderr
            } else {
                &result.stdout
            };
            let output = String::from_utf8_lossy(output);
            assert!(output.contains("GPS-CONF-014"), "{case}/{mode}");
            assert!(output.contains("ref-main"), "{case}/{mode}");
            assert!(output.contains("chr1"), "{case}/{mode}");
            assert!(!destination.exists(), "{case}/{mode}");
        }
        fs::remove_file(&input).expect("retirar fixture aislado");
    }

    let input = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/config/invalid-supplied-ploidy-missing.toml");
    let destination = isolated_destination("ploidy-missing.output");
    for mode in ["human", "json"] {
        let result = run(&input, mode, &destination);
        assert_eq!(result.status.code(), Some(2), "missing/{mode}");
        let output = if mode == "human" {
            &result.stderr
        } else {
            &result.stdout
        };
        let output = String::from_utf8_lossy(output);
        assert!(output.contains("GPS-CONF-014"), "missing/{mode}");
        assert!(output.contains("ref"), "missing/{mode}");
        assert!(output.contains("chr1"), "missing/{mode}");
        assert!(!destination.exists(), "missing/{mode}");
    }
}
