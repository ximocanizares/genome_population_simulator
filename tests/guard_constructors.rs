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
//! *Lo que NO cubre, sin adornos.* Solo ve `Diagnostic::format(` y los
//! rechazos de texto suelto escritos como literal. No ve un rechazo de texto
//! que otra función construya a partir de una variable, ni un envoltorio que
//! llame a una constructora con un motivo ajeno (el cierre `reject` del juez del
//! abanico de deslizamiento, por ejemplo): ahí la garantía es la convención de
//! que cada rechazo pase por `reason`, no esta prueba.

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
            "effective_diagnostic",
            "manifest_diagnostic",
            "retag_format_schema",
        ],
        "un Diagnostic::format fuera de los constructores de results.rs queda fuera de la puerta de guardas"
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
        (results, "fn manifest_error("),
        (results, "fn manifest_error_for_schema("),
        (results, "fn physical_error("),
        (results, "fn repeat_section_error("),
        (results, "fn artifact_identity_error("),
        (results, "fn reason<"),
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
        (results, "fn manifest_error_for_schema("),
        (results, "fn artifact_identity_error("),
        (results, "fn reason<"),
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

/// El `GPS-IO-002` de `results.rs` nace solo en su constructora, y los
/// constructores sin anotar de `GPS-FMT-006` solo se llaman desde donde deben.
///
/// *Por qué existe.* `effective_diagnostic` y `manifest_diagnostic` no anotan
/// nada, y no pueden hacerlo: los llaman los envoltorios que reciben un motivo
/// ya anotado por `reason`. Pero eso los hace una vía abierta por la que
/// construir un rechazo que la puerta de guardas no ve (el mutante M09 de la
/// revisión del 04/10/2026 la usaba y sobrevivía a la batería). Esta prueba
/// limita **léxicamente** dónde pueden aparecer.
///
/// *Lo que NO cubre.* Un rechazo nuevo escrito dentro de una de las funciones
/// permitidas, o un `Diagnostic::io` con otro código.
#[test]
fn unrecorded_diagnostic_builders_are_called_only_from_their_wrappers() {
    let results = non_test_part(include_str!("../src/results.rs"));
    assert_eq!(
        enclosing_functions(results, "effective_diagnostic("),
        [
            "effective_diagnostic",
            "effective_error_for_schema",
            "effective_error_from_reason",
            "validate_published_slippage_fans",
        ]
    );
    assert_eq!(
        enclosing_functions(results, "manifest_diagnostic("),
        [
            "manifest_diagnostic",
            "manifest_error_for_schema",
            "manifest_error_from_reason",
            "repeat_section_error_from_reason",
        ]
    );
    assert_eq!(
        enclosing_functions(results, "\"GPS-IO-002\""),
        ["artifact_identity_error"]
    );
    assert_eq!(misused_from_reason_calls(results), Vec::<String>::new());
}

/// Líneas que llaman a un `*_from_reason(` fuera de la forma única
/// `.map_err(|detail| ..._from_reason(ruta, detail))`: con un literal como
/// motivo, el rechazo no lo anotaría nadie.
fn misused_from_reason_calls(source: &str) -> Vec<String> {
    source
        .lines()
        .map(str::trim)
        .filter(|line| !line.starts_with("//") && !line.starts_with("fn "))
        .filter(|line| line.contains("_from_reason("))
        .filter(|line| !(line.contains(".map_err(|detail| ") && line.contains(", detail))")))
        .map(str::to_owned)
        .collect()
}

#[test]
fn the_wrapper_scans_see_a_misuse() {
    let sample = "fn sneaky() {\n    Err(effective_diagnostic(p, \"s\", \"x\"))\n}\n";
    assert_eq!(
        enclosing_functions(sample, "effective_diagnostic("),
        ["sneaky"]
    );
    assert_eq!(
        misused_from_reason_calls("    Err(manifest_error_from_reason(path, \"x\"))\n"),
        ["Err(manifest_error_from_reason(path, \"x\"))"]
    );
    assert!(
        misused_from_reason_calls(
            "fn manifest_error_from_reason(p: &Path) {}\n  .map_err(|detail| manifest_error_from_reason(path, detail))?;\n"
        )
        .is_empty()
    );
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

/// Los rechazos de texto suelto de `results.rs` pasan por `reason`.
///
/// *Por qué existe.* Los validadores que devuelven `Result<_, &'static str>`
/// no construyen un diagnóstico: lo hace quien los llama, con una sola
/// llamada que recoge todos sus rechazos, y la puerta de cobertura de guardas
/// los contaba como uno. `reason` anota la línea de cada rechazo (DT-001).
///
/// *Lo que ve.* Un `Err("literal")` o un `.ok_or("literal")` sin pasar por
/// `reason`, el `Err(detail)` / `.ok_or(detail)` del auditor de la
/// normalización, un `reject("literal")` del cierre del juez del abanico, y
/// cualquier `reason(` que no vaya precedido de `Err(`, `||`, `reject(`,
/// `reject(&`, `|error|` o `|_|`: en particular `ok_or(reason(..))`, que
/// anotaría siempre, también cuando no hay rechazo.
///
/// *Lo que NO ve, y solo caza la puerta por recuento* (el número de
/// comprobaciones baja): un `map_err(|error| format!(..))` al que se le quita
/// el `reason`, porque ese patrón existe legítimamente fuera de los
/// validadores; y un texto que otra función construya y devuelva en
/// `Err(variable)` con otro nombre.
#[test]
fn text_rejections_of_results_go_through_reason() {
    let results = non_test_part(include_str!("../src/results.rs"));
    assert_eq!(bare_text_rejections(results), Vec::<String>::new());
}

fn bare_text_rejections(source: &str) -> Vec<String> {
    let code: String = source
        .lines()
        .filter(|line| !line.trim_start().starts_with("//"))
        .collect::<Vec<_>>()
        .join("\n");
    let squeezed: String = code.chars().filter(|c| !c.is_whitespace()).collect();
    let mut found: Vec<String> = [
        "Err(\"",
        ".ok_or(\"",
        "Err(detail)",
        ".ok_or(detail)",
        "reject(\"",
    ]
    .into_iter()
    .filter(|needle| squeezed.contains(needle))
    .map(str::to_owned)
    .collect();
    let allowed = ["Err(", "||", "reject(", "reject(&", "|error|", "|_|"];
    let bytes = squeezed.as_bytes();
    let mut from = 0;
    while let Some(offset) = squeezed[from..].find("reason(") {
        let at = from + offset;
        from = at + "reason(".len();
        let identifier_before =
            at > 0 && (bytes[at - 1].is_ascii_alphanumeric() || bytes[at - 1] == b'_');
        if identifier_before {
            continue; // `manifest_error_from_reason(` y similares
        }
        if !allowed
            .iter()
            .any(|prefix| squeezed[..at].ends_with(prefix))
        {
            let context_start = squeezed[..at]
                .char_indices()
                .rev()
                .nth(12)
                .map_or(0, |(i, _)| i);
            found.push(format!("{}reason(", &squeezed[context_start..at]));
        }
    }
    found
}

/// Control positivo: el instrumento ve un rechazo de texto sin `reason`.
#[test]
fn the_text_rejection_scan_sees_a_bare_one() {
    assert_eq!(
        bare_text_rejections("fn f() -> Result<(), &'static str> {\n    Err(\n \"x\"\n)\n}\n"),
        ["Err(\""]
    );
    assert_eq!(
        bare_text_rejections("fn f() { v.ok_or(\"x\")?; }"),
        [".ok_or(\""]
    );
    assert!(bare_text_rejections("fn f() { Err(reason(\"x\")) }\n// Err(\"y\")\n").is_empty());
    assert_eq!(
        bare_text_rejections("fn f() { v.ok_or(reason(\"x\"))?; }"),
        ["nf(){v.ok_or(reason("]
    );
    assert_eq!(
        bare_text_rejections("fn f() { v.ok_or_else(|| reject(\"x\"))?; }"),
        ["reject(\""]
    );
    assert!(
        bare_text_rejections(
            "fn f() { v.ok_or_else(|| reason(\"x\"))?; w.map_err(|error| reason(format!(\"{error}\")))?; reject(&reason(e)); }"
        )
        .is_empty()
    );
}
