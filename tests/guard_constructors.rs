//! Los rechazos de formato de `src/results.rs` y `src/event_manifest.rs` nacen
//! en sus constructoras.
//!
//! *Por qué existe.* La puerta de cobertura de guardas cuenta y anota solo los
//! rechazos que pasan por las constructoras de error de cada módulo. Un
//! `GPS-FMT-006` construido a mano en otra función quedaría fuera de la puerta:
//! ni entraría en el denominador ni se anotaría al dispararse. Ya ocurrió en el
//! verificador del VCF, donde veinte guardas muertas eran invisibles.
//!
//! *Guarda léxica.* Lee el texto y busca cada `Diagnostic::format(` en la parte
//! que no es de pruebas. Las funciones permitidas son las constructoras y las
//! que **reetiquetan** un error ya rechazado (no deciden un rechazo nuevo).
//!
//! *Lo que NO cubre, sin adornos.* Solo ve `Diagnostic::format(`. Los
//! validadores de `results.rs` que devuelven `Err("texto")` y los rechazos que
//! pasan por `manifest_error` no están bajo la puerta (véase DT-001); esta
//! prueba no los vigila.

/// Nombre de la función de columna cero que contiene cada aparición de `needle`.
fn enclosing_functions(source: &str, needle: &str) -> Vec<String> {
    let mut found = Vec::new();
    let mut current = String::new();
    for line in source.lines() {
        let declaration = line
            .strip_prefix("pub fn ")
            .or_else(|| line.strip_prefix("fn "))
            .or_else(|| line.strip_prefix("pub(crate) fn "));
        if let Some(rest) = declaration {
            current = rest
                .split(|c: char| !(c.is_alphanumeric() || c == '_'))
                .next()
                .unwrap_or("")
                .to_owned();
        }
        if line.trim_start().starts_with("//") {
            continue;
        }
        if line.contains(needle) {
            found.push(current.clone());
        }
    }
    found.sort();
    found.dedup();
    found
}

fn non_test_part(source: &str) -> &str {
    source
        .split("\n#[cfg(test)]\n")
        .next()
        .expect("el módulo tiene una parte que no es de pruebas")
}

#[test]
fn format_rejections_of_results_and_event_manifest_are_built_only_in_their_constructors() {
    let results = non_test_part(include_str!("../src/results.rs"));
    let event = non_test_part(include_str!("../src/event_manifest.rs"));
    assert_eq!(
        enclosing_functions(results, "Diagnostic::format("),
        [
            "effective_error_for_schema",
            "manifest_error_for_schema",
            "retag_format_schema",
        ],
        "un Diagnostic::format fuera de las constructoras de results.rs queda fuera de la puerta de guardas"
    );
    assert_eq!(
        enclosing_functions(event, "Diagnostic::format("),
        ["event_error", "retag_event_manifest_schema"],
        "un Diagnostic::format fuera de la constructora de event_manifest.rs queda fuera de la puerta de guardas"
    );
    // Cada constructora vigilada por la puerta lleva la anotación que la ata
    // al barrido y llama a `record`.
    for (source, signature) in [
        (results, "fn effective_error("),
        (results, "fn effective_error_for_schema("),
        (event, "fn event_error("),
    ] {
        let start = source.find(signature).expect("la constructora existe");
        assert!(
            source[..start].ends_with("#[track_caller]\n"),
            "{signature} tiene que llevar #[track_caller] justo encima"
        );
    }
    for (source, signature) in [
        (results, "fn effective_error_for_schema("),
        (event, "fn event_error("),
    ] {
        let start = source.find(signature).expect("la constructora existe");
        let body = &source[start..];
        let end = body.find("\n}\n").expect("la constructora termina");
        assert!(
            body[..end].contains("crate::guard_trace::record();"),
            "{signature} tiene que anotar su disparo"
        );
    }
}

/// Control positivo: el instrumento ve una construcción a mano donde la hay.
#[test]
fn the_lexical_scan_sees_a_hand_built_rejection() {
    let sample = "fn honest() {}\nfn sneaky() {\n    Diagnostic::format(\"x\", \"y\")\n}\n";
    assert_eq!(
        enclosing_functions(sample, "Diagnostic::format("),
        ["sneaky"]
    );
    assert!(
        enclosing_functions("fn a() {}\n// Diagnostic::format(\n", "Diagnostic::format(")
            .is_empty()
    );
}
