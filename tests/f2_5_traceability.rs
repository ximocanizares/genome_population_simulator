use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

const PLAN: &str = include_str!("../design/plan-implementacion-f2-0.1.md");
const TRACEABILITY: &str = include_str!("../design/trazabilidad-requisitos-0.1.md");
const RESUMPTION: &str = include_str!("../design/punto-reanudacion.md");
const CLOSURE: &str = include_str!("../evidence/f2/cierre.md");

#[test]
fn local_f2_5_closure_has_no_planned_steps_or_requirements() {
    for step in 1..=7 {
        let prefix = format!("- ✅ F2.5.{step}.");
        assert!(
            PLAN.lines().any(|line| line.starts_with(&prefix)),
            "{prefix}"
        );
    }
    assert!(!PLAN.contains("- ☐ F2.5."));
    assert!(PLAN.contains("Depende de: F2.5.1–F2.5.6."));
    assert!(PLAN.contains("DTC-066 externo con `bcftools`/`tabix` 1.24 sigue pendiente"));

    let expected = BTreeSet::from([
        "REQ-FMT-05",
        "REQ-FMT-06",
        "REQ-REC-26",
        "REQ-VER-02",
        "REQ-VER-06",
        "REQ-VER-12",
        "REQ-VER-15",
        "REQ-VER-16",
        "REQ-VER-17",
        "REQ-VER-30",
        "REQ-VER-32",
        "REQ-VER-33",
        "REQ-VER-34",
        "REQ-VER-36",
        "REQ-VER-37",
        "REQ-VER-39",
        "REQ-VER-40",
        "REQ-VER-41",
        "REQ-VER-43",
        "REQ-VER-44",
        "REQ-VER-47",
        "REQ-VAL-28",
        "REQ-MAN-27",
        "REQ-MAN-57",
        "REQ-SAL-21",
    ]);
    let rows = TRACEABILITY
        .lines()
        .filter(|line| line.contains("| F2.5 |"))
        .collect::<Vec<_>>();
    let observed = rows
        .iter()
        .map(|line| line.split('|').nth(1).expect("requirement column").trim())
        .collect::<BTreeSet<_>>();
    assert_eq!(observed, expected);

    for row in rows {
        assert!(row.contains("Implementada localmente F2.5"), "{row}");
        assert!(!row.contains("Planificada"), "{row}");
        assert!(row.contains("tests/") || row.contains("src/"), "{row}");
        let references = evidence_references(row);
        assert!(
            !references.is_empty(),
            "no parsed evidence references: {row}"
        );
        for (path, symbol) in references {
            assert_evidence_symbol_exists(&path, &symbol, row);
        }
    }

    for requirement in ["REQ-FMT-07", "REQ-SAL-12"] {
        let row = traceability_row(requirement);
        assert!(row.contains("DTC-066"), "{row}");
        assert!(row.contains("pendiente"), "{row}");
        assert!(row.contains("local"), "{row}");
        assert!(!row.contains("consulta 1.24"), "{row}");
        assert!(!row.contains("`bcftools norm` contra FASTA"), "{row}");
    }
}

#[test]
fn conditional_closure_keeps_external_acceptance_and_next_action_pending() {
    assert!(CLOSURE.contains("cerrada localmente de forma condicional"));
    assert!(CLOSURE.contains("DTC-066 externo pendiente: condición del cierre"));
    assert!(CLOSURE.contains("queda declarada como superada aquí."));
    for document in [PLAN, TRACEABILITY, RESUMPTION, CLOSURE] {
        assert_no_dtc_066_acceptance_claim(document);
    }

    for commit in ["f1f3205", "bb99a5d", "839f56e"] {
        assert!(RESUMPTION.contains(commit), "missing commit {commit}");
    }
    assert!(RESUMPTION.contains("## Siguiente acción vigente"));
    assert!(RESUMPTION.contains("Provisionar `bcftools 1.24`"));
    assert!(RESUMPTION.contains("F2.5.1–F2.5.7 están cerradas localmente"));
    assert!(RESUMPTION.contains("F2.5.1–F2.5.7 siguen sin\n  commit"));
    assert!(!RESUMPTION.contains("Revisar la pausa F2.4 y confirmar el cierre local"));
}

#[test]
fn dtc_066_claim_and_pending_are_mutually_exclusive() {
    let mutant = "DTC-066 está superada y pendiente.";
    assert!(has_affirmative_dtc_066_acceptance_claim(mutant));
    assert!(has_dtc_066_pending(mutant));
    assert!(
        !dtc_066_state_is_consistent(mutant),
        "the exact contradictory mutant must be rejected"
    );
}

#[test]
fn evidence_symbols_must_have_a_test_attribute() {
    let test_source = "#[test]\nfn cited_test() {}";
    let helper_source = "fn cited_test() {}";
    let prior_test_source = "#[test]\nfn previous() {}\nfn cited_test() {}";
    assert!(evidence_symbol_is_test(test_source, "cited_test"));
    assert!(!evidence_symbol_is_test(helper_source, "cited_test"));
    assert!(!evidence_symbol_is_test(prior_test_source, "cited_test"));
}

fn evidence_references(row: &str) -> Vec<(String, String)> {
    let mut references = Vec::new();
    let mut remainder = row;
    while let Some(start) = remainder.find('`') {
        let candidate = &remainder[start + 1..];
        let Some(end) = candidate.find('`') else {
            break;
        };
        let citation = &candidate[..end];
        if let Some((path, symbols)) = citation.split_once("::")
            && (path.starts_with("tests/") || path.starts_with("src/"))
        {
            let symbols = symbols.strip_prefix("tests::").unwrap_or(symbols);
            let symbols = symbols.trim_matches(['{', '}']);
            for symbol in symbols
                .split(',')
                .map(str::trim)
                .filter(|symbol| !symbol.is_empty())
            {
                references.push((path.to_owned(), symbol.to_owned()));
            }
        }
        remainder = &candidate[end + 1..];
    }
    references
}

fn traceability_row(requirement: &str) -> &str {
    TRACEABILITY
        .lines()
        .find(|line| line.starts_with(&format!("| {requirement} |")))
        .expect("traceability row")
}

fn assert_evidence_symbol_exists(path: &str, symbol: &str, row: &str) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let file = root.join(path);
    assert!(
        file.is_file(),
        "missing evidence file {} in {row}",
        file.display()
    );
    let source = fs::read_to_string(&file).expect("read evidence source");
    assert!(
        evidence_symbol_is_test(&source, symbol),
        "missing #[test] symbol {symbol} in {} cited by {row}",
        file.display()
    );
}

fn evidence_symbol_is_test(source: &str, symbol: &str) -> bool {
    let declaration = format!("fn {symbol}(");
    let Some(index) = source.find(&declaration) else {
        return false;
    };
    source[..index]
        .lines()
        .rev()
        .map(str::trim)
        .take_while(|line| line.is_empty() || line.starts_with("#["))
        .any(|attribute| {
            attribute == "#[test]"
                || attribute.starts_with("#[tokio::test")
                || attribute.starts_with("#[rstest")
                || attribute.starts_with("#[test_case")
        })
}

fn assert_no_dtc_066_acceptance_claim(document: &str) {
    assert!(
        dtc_066_state_is_consistent(document),
        "DTC-066 external acceptance cannot be both claimed and pending"
    );
}

fn dtc_066_state_is_consistent(document: &str) -> bool {
    !has_affirmative_dtc_066_acceptance_claim(document) && has_dtc_066_pending(document)
}

fn has_affirmative_dtc_066_acceptance_claim(document: &str) -> bool {
    document.replace('\n', " ").split(". ").any(|sentence| {
        let normalized = sentence.to_ascii_lowercase();
        let claims_external_acceptance = normalized.contains("aceptación externa superada")
            || (normalized.contains("dtc-066") && normalized.contains("superada"));
        let explicitly_negated = normalized.trim_start().starts_with("no ")
            || normalized.contains(" no ")
            || normalized.contains("ninguna");
        claims_external_acceptance && !explicitly_negated
    })
}

fn has_dtc_066_pending(document: &str) -> bool {
    document.replace('\n', " ").split(". ").any(|sentence| {
        let normalized = sentence.to_ascii_lowercase();
        normalized.contains("dtc-066") && normalized.contains("pendiente")
    })
}
