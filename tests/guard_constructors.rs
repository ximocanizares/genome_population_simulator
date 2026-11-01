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
            "validate_published_detected_fans",
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

// ---------------------------------------------------------------------------
// Los nueve jueces que entraron en la puerta el 05/10/2026
// ---------------------------------------------------------------------------

/// Un juez que canaliza sus rechazos por `rejected` (la constructora genérica de
/// `src/guard_trace.rs`) o por constructoras propias anotadas.
struct Judge {
    /// Ruta del módulo, para los mensajes.
    path: &'static str,
    /// Texto completo del módulo.
    source: &'static str,
    /// Constructoras propias anotadas con `#[track_caller]` que llaman a
    /// `record`: cada llamada suya es un punto de rechazo, así que nunca pueden
    /// ir en `ok_or(..)`, que las llamaría también cuando no hay rechazo.
    constructors: &'static [&'static str],
    /// Constructores **sin anotar** (clausuras, que no admiten
    /// `#[track_caller]`): cada uso tiene que ir envuelto en `rejected(..)`.
    builders: &'static [&'static str],
    /// Prefijos de expresiones de error que **no** son un veredicto sobre lo
    /// juzgado (un fallo de entrada/salida, la propagación de un fallo ya
    /// producido, el generador de matrices) y que, por tanto, no pasan por
    /// `rejected`, con el **número exacto** de usos que se aceptan de cada
    /// uno. Cada exención se justifica en DT-001. Si el recuento no coincide
    /// la prueba falla: un rechazo nuevo escrito con un prefijo exento no
    /// puede colarse, y una exención que ya no se usa tiene que borrarse.
    exempt: &'static [(&'static str, usize)],
}

const NEW_JUDGES: &[Judge] = &[
    Judge {
        path: "src/truth_validation.rs",
        source: include_str!("../src/truth_validation.rs"),
        constructors: &["mismatch"],
        builders: &["missing"],
        // Reetiquetan un error que ya anotó su origen: el de la huella de la
        // referencia, y los de la reconstrucción y la semántica.
        exempt: &[
            ("TruthValidationError::BasesNotMaterialized", 1),
            ("InspectedValidationError::", 4),
        ],
    },
    Judge {
        path: "src/truth_schema.rs",
        source: include_str!("../src/truth_schema.rs"),
        constructors: &[],
        builders: &[],
        exempt: &[],
    },
    Judge {
        path: "src/fastq_validation.rs",
        source: include_str!("../src/fastq_validation.rs"),
        constructors: &["missing_canonical_eof"],
        builders: &["structural"],
        // Reetiqueta el texto que ya anotó `parse_final_header`.
        exempt: &[("FastqValidationError::Structural", 1)],
    },
    Judge {
        path: "src/population_snvs.rs",
        source: include_str!("../src/population_snvs.rs"),
        constructors: &["audit_error", "missing", "invalid"],
        builders: &[],
        // `AuthoritativeFasta` reetiqueta, en dos sitios, un error que ya
        // anotó `parse_metadata_json` o `parse_canonical_tsv`; el generador de
        // matrices no es un juez.
        exempt: &[
            ("PopulationSnvsPreparationError::AuthoritativeFasta", 2),
            ("PopulationSnvsGenerationError", 36),
        ],
    },
    Judge {
        path: "src/f5_publication.rs",
        source: include_str!("../src/f5_publication.rs"),
        constructors: &["document_error", "manifest_error", "failure_mismatch"],
        builders: &[],
        exempt: &[("Box::new(F5PublicationFailure", 2)],
    },
    Judge {
        path: "src/reference_fasta.rs",
        source: include_str!("../src/reference_fasta.rs"),
        constructors: &[
            "empty_header",
            "empty_fasta",
            "invalid_symbol",
            "configuration_mismatch",
            "reference_error",
        ],
        builders: &[],
        exempt: &[("Diagnostic::io(", 2)],
    },
    Judge {
        path: "src/reference_reduction_io.rs",
        source: include_str!("../src/reference_reduction_io.rs"),
        constructors: &["map_error", "out_of_range"],
        builders: &[],
        exempt: &[("Diagnostic::io(", 3)],
    },
    Judge {
        path: "src/integrity.rs",
        source: include_str!("../src/integrity.rs"),
        constructors: &[],
        builders: &[],
        exempt: &[],
    },
    Judge {
        path: "src/technical_profile.rs",
        source: include_str!("../src/technical_profile.rs"),
        constructors: &["invalid"],
        builders: &[],
        exempt: &[],
    },
    // El juez del catálogo de microsatélites, que entró el 01/11/2026 (entrega
    // 4 del plan de detección): un solo constructor anotado.
    Judge {
        path: "src/repeat_catalog_validation.rs",
        source: include_str!("../src/repeat_catalog_validation.rs"),
        constructors: &["catalog_reject"],
        builders: &[],
        exempt: &[],
    },
];

/// El texto sin los módulos `#[cfg(test)]`, pero **con** el código que venga
/// después de ellos: `fastq_validation.rs` tiene el juez entero detrás de un
/// módulo de pruebas, y cortar en el primero lo dejaba sin vigilar. Igual que
/// `scripts/guard-sites.awk`, un bloque de pruebas empieza en la anotación y
/// termina en la primera llave de cierre en la columna cero.
fn without_test_blocks(source: &str) -> String {
    let mut kept = String::new();
    let mut rest = source;
    while let Some(start) = rest.find("\n#[cfg(test)]\n") {
        kept.push_str(&rest[..start]);
        let block = &rest[start + 1..];
        let end = block.find("\n}\n").map_or(block.len(), |i| i + 3);
        rest = &block[end..];
    }
    kept.push_str(rest);
    kept
}

/// El código del módulo sin comentarios de línea y con cada tramo de espacios
/// en blanco reducido a uno solo.
fn squeezed_code(source: &str) -> String {
    without_test_blocks(source)
        .lines()
        .filter(|line| !line.trim_start().starts_with("//"))
        .collect::<Vec<_>>()
        .join(" ")
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

/// Un rechazo construido sin pasar por una constructora anotada.
///
/// Mira el comienzo de la expresión de error en `Err(..)`, `.ok_or_else(|| ..)`,
/// `.map_err(|e| ..)` y `.ok_or(..)`. Es un rechazo escondido si esa expresión
/// empieza por una mayúscula (una variante o un tipo de diagnóstico), por un
/// literal de texto, por `D::Error::custom(` o por un constructor sin anotar
/// del juez. En `.ok_or(..)` lo es además cualquier constructora, aunque esté
/// anotada: se evalúa siempre y anotaría también cuando no hay rechazo.
fn hidden_rejections(source: &str, judge: &Judge) -> Vec<String> {
    scan_rejections(source, judge).0
}

/// Cuántas veces aparece cada exención del juez, en el orden de `judge.exempt`.
fn exempt_uses(source: &str, judge: &Judge) -> Vec<usize> {
    scan_rejections(source, judge).1
}

/// Los rechazos escondidos y los usos de cada exención.
fn scan_rejections(source: &str, judge: &Judge) -> (Vec<String>, Vec<usize>) {
    let mut exempt_counts = vec![0; judge.exempt.len()];
    let code = squeezed_code(source);
    let bytes = code.as_bytes();
    let mut found = Vec::new();
    for opener in ["Err(", ".ok_or_else(||", ".ok_or(", ".map_err(|"] {
        let mut from = 0;
        while let Some(offset) = code[from..].find(opener) {
            let at = from + offset;
            from = at + opener.len();
            if opener == "Err("
                && at > 0
                && (bytes[at - 1].is_ascii_alphanumeric() || bytes[at - 1] == b'_')
            {
                continue; // `MapErr(` y similares
            }
            let mut head = &code[from..];
            if opener == ".map_err(|" {
                let Some(close) = head.find('|') else {
                    continue;
                };
                head = &head[close + 1..];
            }
            head = head.trim_start();
            head = head.strip_prefix('{').unwrap_or(head).trim_start();
            let is_error_expression = head.starts_with(|c: char| c.is_ascii_uppercase())
                || head.starts_with('"')
                || head.starts_with("D::Error::custom(")
                || judge
                    .builders
                    .iter()
                    .any(|name| head.starts_with(&format!("{name}(")));
            let always_evaluated = opener == ".ok_or("
                && (head.starts_with("rejected(")
                    || judge
                        .constructors
                        .iter()
                        .any(|name| head.starts_with(&format!("{name}("))));
            let exemption = judge
                .exempt
                .iter()
                .position(|(prefix, _)| head.starts_with(prefix));
            if let (true, Some(index)) = (is_error_expression, exemption) {
                exempt_counts[index] += 1;
            }
            if (is_error_expression && exemption.is_none()) || always_evaluated {
                found.push(format!(
                    "{opener}{}",
                    head.chars().take(60).collect::<String>()
                ));
            }
        }
    }
    (found, exempt_counts)
}

/// Ningún rechazo de los jueces añadidos se salta la constructora anotada.
///
/// *Por qué existe.* Es lo que mantiene honesta la puerta de cobertura en estos
/// módulos: un `Err(Variante { .. })` sin envolver en `rejected(..)`, o un
/// `ok_or(rejected(..))` que anota siempre, sería una comprobación que la
/// puerta no ve o ve disparada sin estarlo (DT-001).
///
/// *Lo que NO cubre, sin adornos.* Un rechazo cuya expresión de error empiece
/// por una minúscula que no sea un constructor listado (una variable
/// intermedia, una función de otro módulo); la propagación con `?`, que no
/// decide nada; los ayudantes que reciben el motivo como argumento y hacen un
/// único rechazo para muchos llamantes, que se cuentan una vez (`require` en
/// el juez de casos de `truth_validation.rs`); y **tres casos que solo caza la
/// puerta por recuento** (la cifra de comprobaciones baja, o la de sin
/// disparar sube): un rechazo escrito con una variable intermedia,
/// `<D::Error as serde::de::Error>::custom(..)` (la prueba solo conoce
/// `D::Error::custom(`) y un rechazo escrito con un prefijo exento, cuyo
/// recuento exacto sí se exige abajo pero que, si se cambia a la vez que el
/// recuento declarado, ya no se ve.
#[test]
fn no_rejection_of_the_added_judges_skips_its_constructor() {
    for judge in NEW_JUDGES {
        assert_eq!(
            hidden_rejections(judge.source, judge),
            Vec::<String>::new(),
            "{}: un rechazo que no pasa por rejected(..) ni por una constructora anotada queda fuera de la puerta de guardas",
            judge.path
        );
        let declared: Vec<usize> = judge.exempt.iter().map(|(_, n)| *n).collect();
        assert_eq!(
            exempt_uses(judge.source, judge),
            declared,
            "{}: el número de usos de cada exención ({:?}) no es el declarado; un rechazo nuevo con ese prefijo no pasaría por rejected(..)",
            judge.path,
            judge.exempt.iter().map(|(p, _)| *p).collect::<Vec<_>>()
        );
    }
}

/// Las constructoras de los jueces añadidos llevan `#[track_caller]` justo
/// encima y llaman a `record`, y `rejected` también.
#[test]
fn the_constructors_of_the_added_judges_are_annotated_and_record() {
    for judge in NEW_JUDGES {
        for name in judge.constructors {
            assert_eq!(
                constructor_defect(judge.source, name),
                None,
                "{}: {name}",
                judge.path
            );
        }
    }
    let trace = include_str!("../src/guard_trace.rs");
    assert_eq!(constructor_defect(trace, "rejected"), None);
}

/// Qué le falta a la constructora `name` de `source`, si algo.
fn constructor_defect(source: &str, name: &str) -> Option<&'static str> {
    let start = source
        .find(&format!("fn {name}("))
        .or_else(|| source.find(&format!("fn {name}<")))?;
    let line_start = source[..start].rfind('\n').map_or(0, |i| i + 1);
    if !source[..line_start].ends_with("#[track_caller]\n") {
        return Some("falta #[track_caller] justo encima");
    }
    let body = &source[start..];
    let end = body.find("\n}\n")?;
    let body = &body[..end];
    if !(body.contains("crate::guard_trace::record();") || body.contains("    record();")) {
        return Some("no llama a record");
    }
    None
}

/// La lista de jueces de esta prueba y el array `judges` de la puerta no
/// pueden diverger: un módulo vigilado aquí y ausente allí no se mediría.
#[test]
fn every_added_judge_is_in_the_coverage_gate() {
    let gate = include_str!("../scripts/check-guard-coverage.sh");
    let baseline = include_str!("../scripts/guard-coverage-baseline.tsv");
    for judge in NEW_JUDGES {
        assert!(
            gate.contains(&format!("\"{}|", judge.path)),
            "{} no está en el array judges de scripts/check-guard-coverage.sh",
            judge.path
        );
        assert!(
            baseline.contains(&format!("\n{}\t", judge.path)),
            "{} no tiene fila en scripts/guard-coverage-baseline.tsv",
            judge.path
        );
    }
}

fn sample_judge() -> Judge {
    Judge {
        path: "muestra",
        source: "",
        constructors: &["ctor"],
        builders: &["build"],
        exempt: &[("Diagnostic::io(", 1)],
    }
}

/// Controles positivos: el instrumento ve cada forma de saltarse la constructora.
#[test]
fn the_hidden_rejection_scan_sees_each_way_to_skip_the_constructor() {
    let judge = sample_judge();
    let sees = |sample: &str| hidden_rejections(sample, &judge);
    assert_eq!(
        sees("fn f() { return Err(MyError::Bad { x: 1 }); }").len(),
        1,
        "Err(Variante) sin envolver"
    );
    assert_eq!(sees("fn f() { Err(\"texto\") }").len(), 1, "Err(literal)");
    assert_eq!(
        sees("fn f() { v.ok_or(MyError::Bad)?; }").len(),
        1,
        "ok_or(Variante)"
    );
    assert_eq!(
        sees("fn f() { v.ok_or(rejected(MyError::Bad))?; }").len(),
        1,
        "ok_or(rejected(..)) anota siempre"
    );
    assert_eq!(
        sees("fn f() { v.ok_or(ctor(1))?; }").len(),
        1,
        "ok_or(constructora) anota siempre"
    );
    assert_eq!(
        sees("fn f() { v.ok_or_else(|| MyError::Bad)?; }").len(),
        1,
        "ok_or_else sin envolver"
    );
    assert_eq!(
        sees("fn f() { v.ok_or_else(|| {\n MyError::Bad { x: 1 }\n})?; }").len(),
        1,
        "ok_or_else con bloque sin envolver"
    );
    assert_eq!(
        sees("fn f() { v.map_err(|_| MyError::Bad)?; }").len(),
        1,
        "map_err sin envolver"
    );
    assert_eq!(
        sees("fn f() { v.map_err(|e| build(e))?; }").len(),
        1,
        "constructor sin anotar sin envolver"
    );
    assert_eq!(
        sees("fn f() { Err(D::Error::custom(\"x\")) }").len(),
        1,
        "error de serde sin envolver"
    );
    // Lo que no es un rechazo escondido.
    for clean in [
        "fn f() { return Err(rejected(MyError::Bad { x: 1 })); }",
        "fn f() { v.ok_or_else(|| rejected(MyError::Bad))?; }",
        "fn f() { v.ok_or_else(|| ctor(1))?; }",
        "fn f() { v.map_err(|_| rejected(MyError::Bad))?; }",
        "fn f() { v.map_err(|e| Diagnostic::io(\"GPS-IO-001\", \"x\"))?; }",
        "fn f() { match r { Err(_) => 1, Err(e) => g(e) } }",
        "// Err(MyError::Bad)\nfn f() {}\n",
    ] {
        assert_eq!(sees(clean), Vec::<String>::new(), "{clean}");
    }
}

/// Control positivo: el código que viene después de un módulo de pruebas se
/// sigue vigilando, y el del propio módulo no.
#[test]
fn the_scan_keeps_the_code_after_a_test_module() {
    let judge = sample_judge();
    let source = "fn a() {}\n#[cfg(test)]\nmod t {\n    fn x() { Err(MyError::InTest) }\n}\nfn b() { Err(MyError::AfterTests) }\n";
    let found = hidden_rejections(source, &judge);
    assert_eq!(found.len(), 1);
    assert!(found[0].contains("AfterTests"), "{found:?}");
}

/// Control positivo: el instrumento ve una constructora sin anotar, sin
/// `record`, y la que sí cumple.
#[test]
fn the_constructor_scan_sees_a_missing_annotation_or_record() {
    let good =
        "#[track_caller]\nfn ctor(x: u8) -> E {\n    crate::guard_trace::record();\n    E(x)\n}\n";
    assert_eq!(constructor_defect(good, "ctor"), None);
    let no_attribute = "fn ctor(x: u8) -> E {\n    crate::guard_trace::record();\n    E(x)\n}\n";
    assert_eq!(
        constructor_defect(no_attribute, "ctor"),
        Some("falta #[track_caller] justo encima")
    );
    let no_record = "#[track_caller]\nfn ctor(x: u8) -> E {\n    E(x)\n}\n";
    assert_eq!(
        constructor_defect(no_record, "ctor"),
        Some("no llama a record")
    );
    let generic = "#[track_caller]\npub(crate) fn rejected<E>(error: E) -> E {\n    record();\n    error\n}\n";
    assert_eq!(constructor_defect(generic, "rejected"), None);
}

/// Control positivo: el recuento de usos de una exención sube cuando se añade
/// un rechazo nuevo con ese prefijo, y no cuenta los envueltos en `rejected`.
#[test]
fn the_exemption_count_sees_a_new_use_of_an_exempt_prefix() {
    let judge = sample_judge();
    let one = "fn f() { v.map_err(|e| Diagnostic::io(\"GPS-IO-001\", \"x\"))?; }";
    assert_eq!(exempt_uses(one, &judge), [1]);
    let two = "fn f() { v.map_err(|e| Diagnostic::io(\"GPS-IO-001\", \"x\"))?;\n return Err(Diagnostic::io(\"GPS-IO-002\", \"y\")); }";
    assert_eq!(exempt_uses(two, &judge), [2]);
    let wrapped = "fn f() { return Err(rejected(Diagnostic::io(\"GPS-IO-002\", \"y\"))); }";
    assert_eq!(exempt_uses(wrapped, &judge), [0]);
}
