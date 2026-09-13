//! F5.6.3. Auditoría de cierre de F5.
//!
//! Este hito no añade contrato biológico: audita. Lo que fija es que la
//! proyección F5 de la matriz de trazabilidad sea **exactamente** la que el
//! plan declara, que cada una de sus filas cite pruebas que existen, que las
//! prohibiciones de alcance de F5 sigan siendo ciertas en el árbol, que los
//! auditores compuestos acepten también la mitad «fallo» sobre el destino ya
//! renombrado y que la atomicidad —hasta ahora demostrada en frontera de
//! biblioteca— se sostenga entre **procesos**, que es la deuda que F5.6.1 dejó
//! abierta al cablear el primer `run` productivo.
//!
//! Ninguna de las cuentas de este fichero se compara consigo misma: el
//! recuento sale de contar la tabla y se contrasta contra el literal del
//! contrato y contra las dos declaraciones en prosa que lo repiten —la de la
//! propia matriz y la de la entrada F5.6.3 del plan—, que se mantienen a mano
//! y por separado.

use std::collections::BTreeSet;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use genome_population_simulator::application::{
    CliValidatedConfig, EffectiveF5Config, execute_validated_f5_config,
    load_and_validate_cli_config,
};
use genome_population_simulator::f5_oracle::{
    F5DeclaredLink, F5FailureDeclaration, audit_f5_failure_package,
};
use serde_json::Value;
use sha2::{Digest, Sha256};

const TRACEABILITY: &str = include_str!("../design/trazabilidad-requisitos-0.1.md");
const PLAN: &str = include_str!("../design/plan-implementacion-f5-0.1.md");
const CLOSURE: &str = include_str!("../evidence/f5.6/cierre.md");

/// El literal del contrato. El plan de F5.6.3 pide «las 46 filas exactas».
const DECLARED_F5_ROWS: usize = 46;

/// Los tres REQ que el plan nombra uno a uno para que nadie los pierda.
const NAMED_REQUIREMENTS: [&str; 3] = ["REQ-VAL-54", "REQ-MAN-70", "REQ-SAL-02"];

static NONCE: AtomicU64 = AtomicU64::new(0);

// ---------------------------------------------------------------------------
// A. Recuento mecánico de la proyección F5
// ---------------------------------------------------------------------------

/// Una fila `REQ` de la matriz, con sus cuatro columnas.
#[derive(Clone, Copy, Debug)]
struct MatrixRow<'a> {
    requirement: &'a str,
    milestone: &'a str,
    evidence: &'a str,
    line: &'a str,
}

/// Extrae **todas** las filas `REQ` de la matriz, sin filtrar por hito.
fn matrix_rows(document: &str) -> Vec<MatrixRow<'_>> {
    document
        .lines()
        .filter(|line| line.starts_with("| REQ-"))
        .map(|line| {
            let mut cells = line.split('|').map(str::trim);
            cells
                .next()
                .expect("celda vacía antes del primer separador");
            let requirement = cells.next().expect("columna Requisito");
            cells.next().expect("columna ACE");
            let milestone = cells.next().expect("columna Hito");
            let evidence = cells.next().expect("columna Evidencia");
            MatrixRow {
                requirement,
                milestone,
                evidence,
                line,
            }
        })
        .collect()
}

/// `true` cuando **algún** término del campo Hito es F5 o un hito suyo.
///
/// El campo admite varios hitos separados por `/`, de modo que la pertenencia
/// se decide término a término y no por subcadena: `F50.1` no es F5 y
/// `F4.3/DTF-055` tampoco, mientras que `F2.0/F5.6` sí lo es.
fn projects_on_f5(milestone: &str) -> bool {
    milestone
        .split('/')
        .map(str::trim)
        .any(|term| term == "F5" || term.starts_with("F5."))
}

fn f5_rows(document: &str) -> Vec<MatrixRow<'_>> {
    matrix_rows(document)
        .into_iter()
        .filter(|row| projects_on_f5(row.milestone))
        .collect()
}

/// Lee un entero escrito en prosa entre dos marcas literales.
///
/// La prosa de los documentos va plegada a ochenta columnas, de modo que la
/// marca puede cruzar un salto de línea: la búsqueda trabaja sobre el texto con
/// los espacios en blanco colapsados.
fn declared_number(document: &str, prefix: &str, suffix: &str) -> usize {
    let document = &document.split_whitespace().collect::<Vec<_>>().join(" ");
    let start = document
        .find(prefix)
        .unwrap_or_else(|| panic!("el documento tiene que declarar «{prefix}…{suffix}»"))
        + prefix.len();
    let rest = &document[start..];
    let end = rest
        .find(suffix)
        .unwrap_or_else(|| panic!("la declaración «{prefix}» no cierra con «{suffix}»"));
    rest[..end]
        .trim()
        .parse()
        .expect("la cifra declarada tiene que ser un entero")
}

/// La proyección F5 son 46 filas, con 46 REQ distintos y los tres nombrados.
///
/// El recuento es **mecánico**: sale de contar la tabla, no de una lista
/// paralela escrita a mano. Se contrasta contra tres verdades independientes
/// entre sí: el literal del contrato, la prosa de la §1 de la propia matriz y
/// la entrada F5.6.3 del plan. Cambiar la tabla sin cambiar las dos prosas
/// —o al revés— tumba la prueba, que es lo que el hito exige.
#[test]
fn the_f5_projection_of_the_matrix_is_exactly_forty_six_unique_requirements() {
    let rows = f5_rows(TRACEABILITY);
    assert_eq!(
        rows.len(),
        DECLARED_F5_ROWS,
        "la proyección por hito F5 tiene que contener exactamente {DECLARED_F5_ROWS} filas"
    );

    let unique = rows
        .iter()
        .map(|row| row.requirement)
        .collect::<BTreeSet<_>>();
    assert_eq!(
        unique.len(),
        DECLARED_F5_ROWS,
        "los {DECLARED_F5_ROWS} REQ de la proyección tienen que ser distintos"
    );

    for requirement in NAMED_REQUIREMENTS {
        assert!(
            unique.contains(requirement),
            "{requirement} tiene que estar expresamente en la proyección F5"
        );
    }

    // Las dos declaraciones en prosa, que se mantienen aparte de la tabla.
    assert_eq!(
        declared_number(
            TRACEABILITY,
            "La proyección por hito contiene exactamente ",
            " filas"
        ),
        rows.len(),
        "la §1 de la matriz declara otro recuento que el que su tabla sostiene"
    );
    assert_eq!(
        declared_number(PLAN, "guardar evidencia y auditar las ", " filas exactas"),
        rows.len(),
        "la entrada F5.6.3 del plan declara otro recuento que el que la matriz sostiene"
    );

    // Ninguna fila ajena se cuela: todo hito proyectado nombra F5.
    for row in &rows {
        assert!(
            row.milestone.contains("F5"),
            "fila proyectada sin hito F5: {}",
            row.line
        );
    }
}

/// El proyector discrimina, y el recuento cambia cuando la matriz cambia.
///
/// Sin esta prueba el recuento anterior podría pasar con un proyector que
/// aceptase o rechazase de más: una búsqueda por subcadena contaría `F50.1`, y
/// una que exigiera el campo completo perdería `F2.0/F5.6`. Aquí se fija el
/// criterio término a término y, sobre una tabla sintética, que **añadir o
/// quitar una fila F5 mueve el recuento**.
#[test]
fn the_projection_discriminates_a_changed_matrix() {
    for accepted in ["F5", "F5.3", "F5.3/F5.4", "F2.0/F5.6", " F5.6 / F7.6 "] {
        assert!(projects_on_f5(accepted), "{accepted} sí proyecta sobre F5");
    }
    for rejected in [
        "F5X",
        "F50.1",
        "F6.5",
        "F4.3/DTF-055",
        "DTF-051",
        "F1.5",
        "",
    ] {
        assert!(!projects_on_f5(rejected), "{rejected} no proyecta sobre F5");
    }

    let table = "| Requisito | ACE | Hito | Evidencia |\n\
         | --- | --- | --- | --- |\n\
         | REQ-A-01 | ACE-01 | F5.3 | Implementada. |\n\
         | REQ-A-02 | ACE-02 | F6.5 | Planificada. |\n\
         | REQ-A-03 | ACE-03 | F2.0/F5.6 | Parcial. |\n";
    assert_eq!(matrix_rows(table).len(), 3);
    assert_eq!(f5_rows(table).len(), 2);

    let grown = format!("{table}| REQ-A-04 | ACE-04 | F5.4 | Implementada. |\n");
    assert_eq!(
        f5_rows(&grown).len(),
        3,
        "añadir una fila F5 tiene que mover el recuento"
    );
    let shrunk = table.replace("| REQ-A-01 | ACE-01 | F5.3 | Implementada. |\n", "");
    assert_ne!(shrunk, table, "la mutación tiene que mutar");
    assert_eq!(
        f5_rows(&shrunk).len(),
        1,
        "quitar una fila F5 tiene que mover el recuento"
    );
    let retargeted = table.replace(
        "| REQ-A-01 | ACE-01 | F5.3 |",
        "| REQ-A-01 | ACE-01 | F6.3 |",
    );
    assert_ne!(retargeted, table, "la mutación tiene que mutar");
    assert_eq!(
        f5_rows(&retargeted).len(),
        1,
        "cambiar el hito de una fila tiene que mover el recuento"
    );
}

// ---------------------------------------------------------------------------
// B. Auditoría fila a fila
// ---------------------------------------------------------------------------

/// Cada fila F5 dice qué hizo F5 y cita pruebas que existen.
///
/// «Planificada» es la palabra que la §1 de la matriz reserva para lo que aún
/// no está implementado: una fila cuyo hito es F5 y que sigue planificada
/// entera sería alcance F5 sin cerrar. Las filas compartidas con hitos
/// posteriores conservan su parte pendiente, pero tienen que declarar la parte
/// que F5 cierra y respaldarla con una prueba **que existe y es una prueba**.
#[test]
fn every_f5_row_states_what_f5_did_and_cites_tests_that_exist() {
    let rows = f5_rows(TRACEABILITY);
    assert_eq!(rows.len(), DECLARED_F5_ROWS);
    for row in &rows {
        let states_f5_work = (row.evidence.starts_with("Implementada ")
            || row.evidence.starts_with("Parcial "))
            && row.evidence.contains("F5");
        assert!(
            states_f5_work,
            "la fila de {} no declara qué cerró F5: {}",
            row.requirement, row.line
        );
        assert!(
            !row.evidence.starts_with("Planificada"),
            "la fila de {} sigue planificada entera: {}",
            row.requirement,
            row.line
        );
        let references = evidence_references(row.evidence);
        assert!(
            !references.is_empty(),
            "la fila de {} no cita ninguna prueba: {}",
            row.requirement,
            row.line
        );
        for (path, symbol) in references {
            assert_evidence_symbol_exists(&path, &symbol, row.requirement);
        }
    }
}

/// El acta de cierre cita pruebas reales y no esconde el estado del árbol.
///
/// La evidencia de cierre vale por lo que declara: si nombra una prueba, la
/// prueba tiene que existir; y si el árbol tiene fallos preexistentes ajenos a
/// F5, el acta tiene que nombrarlos con su fase.
#[test]
fn the_closure_record_cites_real_tests_and_declares_the_inherited_failures() {
    let references = evidence_references(CLOSURE);
    assert!(
        references.len() >= 8,
        "el acta tiene que respaldar sus afirmaciones con pruebas citadas"
    );
    for (path, symbol) in references {
        assert_evidence_symbol_exists(&path, &symbol, "evidence/f5.6/cierre.md");
    }
    for inherited in [
        "test_f2_5_vcf_ui",
        "test_ui_technical_event_v2_round_trips_to_rust_validation",
    ] {
        assert!(
            CLOSURE.contains(inherited),
            "el acta no puede esconder el fallo preexistente {inherited}"
        );
    }
    for debt in ["DTF-058", "f5-resource-accounting-v8"] {
        assert!(CLOSURE.contains(debt), "el acta tiene que recontar {debt}");
    }
}

/// Referencias `ruta::símbolo` citadas entre acentos graves.
///
/// Admite la forma abreviada `ruta::{a, b}` y las tres raíces que el proyecto
/// usa como evidencia: `tests/`, `src/` y `tests_ui/`. En Python una cita puede
/// interponer la clase —`fichero.py::Clase::test_x`—, de modo que el símbolo es
/// siempre el **último** segmento.
fn evidence_references(text: &str) -> Vec<(String, String)> {
    let mut references = Vec::new();
    let mut remainder = text;
    while let Some(start) = remainder.find('`') {
        let candidate = &remainder[start + 1..];
        let Some(end) = candidate.find('`') else {
            break;
        };
        let citation = &candidate[..end];
        if let Some((path, symbols)) = citation.split_once("::")
            && (path.starts_with("tests/")
                || path.starts_with("src/")
                || path.starts_with("tests_ui/"))
        {
            let symbols = symbols.trim_matches(['{', '}']);
            for symbol in symbols
                .split(',')
                .map(str::trim)
                .filter(|symbol| !symbol.is_empty())
            {
                let symbol = symbol.rsplit("::").next().unwrap_or(symbol);
                references.push((path.to_owned(), symbol.to_owned()));
            }
        }
        remainder = &candidate[end + 1..];
    }
    references
}

fn assert_evidence_symbol_exists(path: &str, symbol: &str, cited_by: &str) {
    let file = Path::new(env!("CARGO_MANIFEST_DIR")).join(path);
    assert!(
        file.is_file(),
        "{cited_by} cita un fichero que no existe: {}",
        file.display()
    );
    let source = fs::read_to_string(&file).expect("leer la fuente citada");
    let found = if path.ends_with(".py") {
        source.contains(&format!("def {symbol}("))
    } else {
        rust_symbol_is_test(&source, symbol)
    };
    assert!(
        found,
        "{cited_by} cita {symbol}, que no es una prueba de {}",
        file.display()
    );
}

/// `true` si `symbol` se declara en `source` con el atributo `#[test]`.
fn rust_symbol_is_test(source: &str, symbol: &str) -> bool {
    let declaration = format!("fn {symbol}(");
    let Some(index) = source.find(&declaration) else {
        return false;
    };
    source[..index]
        .lines()
        .rev()
        .map(str::trim)
        .take_while(|line| line.is_empty() || line.starts_with("#[") || line.starts_with("///"))
        .any(|attribute| attribute == "#[test]")
}

/// El lector de citas distingue una prueba de una función auxiliar.
///
/// Sin este control, una cita mal escrita pasaría por no encontrarse nada y la
/// auditoría anterior sería una aserción que no puede fallar.
#[test]
fn the_citation_reader_discriminates_tests_from_helpers() {
    assert!(rust_symbol_is_test("#[test]\nfn cited() {}", "cited"));
    assert!(rust_symbol_is_test(
        "#[test]\n/// prosa\nfn cited() {}",
        "cited"
    ));
    assert!(!rust_symbol_is_test("fn cited() {}", "cited"));
    assert!(!rust_symbol_is_test(
        "#[test]\nfn previous() {}\nfn cited() {}",
        "cited"
    ));
    assert!(!rust_symbol_is_test("#[test]\nfn other() {}", "cited"));

    let parsed = evidence_references(
        "prosa `tests/f5_6_closure.rs::{alpha, beta}` y `src/f5_layout.rs::gamma`, \
         con clase `tests_ui/caso.py::Clase::test_delta`, \
         más `otro/sitio.rs::epsilon` y `sin_separador`",
    );
    assert_eq!(
        parsed,
        vec![
            ("tests/f5_6_closure.rs".to_owned(), "alpha".to_owned()),
            ("tests/f5_6_closure.rs".to_owned(), "beta".to_owned()),
            ("src/f5_layout.rs".to_owned(), "gamma".to_owned()),
            ("tests_ui/caso.py".to_owned(), "test_delta".to_owned()),
        ]
    );
}

// ---------------------------------------------------------------------------
// C. El hueco que la auditoría destapa: la cota de ploidía de F5
// ---------------------------------------------------------------------------

/// REQ-CONF-65: la admisión F5 rechaza una ploidía fuera de `1..=10`.
///
/// La guarda vive en `src/adapters.rs` y **no tenía prueba**: la batería de
/// admisión de F5.0.1 solo cambiaba `ploidy = 2` por `ploidy = 1`, que es un
/// valor legal. Sin ella, relajar la comprobación a la conversión a `u8`
/// —que es lo que la construcción del efectivo hace después— dejaba pasar
/// `ploidy = 0` y `ploidy = 11` sin que ninguna prueba muriera, y la ploidía
/// efectiva de cada secuencia dejaba de estar acotada por el contrato.
#[test]
fn the_f5_admission_rejects_a_ploidy_outside_one_to_ten() {
    let root = temporary("ploidy-bound");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_pair(&root);
    let original = fs::read_to_string(&general).expect("general");
    assert!(original.contains("ploidy = 2"));

    // Los dos extremos legales se admiten: la cota es inclusiva por los dos
    // lados y el rechazo no puede venir de cualquier valor distinto de 2.
    for legal in ["1", "10"] {
        let candidate = original.replace("ploidy = 2", &format!("ploidy = {legal}"));
        assert_ne!(candidate, original, "la mutación tiene que mutar");
        fs::write(&general, &candidate).expect("general mutado");
        assert!(
            load_and_validate_cli_config(&general).is_ok(),
            "la ploidía {legal} pertenece a 1..=10"
        );
    }

    for illegal in ["0", "11", "250", "4096"] {
        let candidate = original.replace("ploidy = 2", &format!("ploidy = {illegal}"));
        assert_ne!(candidate, original, "la mutación tiene que mutar");
        fs::write(&general, &candidate).expect("general mutado");
        let rejection = load_and_validate_cli_config(&general)
            .err()
            .unwrap_or_else(|| panic!("la ploidía {illegal} no pertenece a 1..=10"));
        assert_eq!(rejection.code(), "GPS-CONF-101", "ploidía {illegal}");
        assert_eq!(
            rejection.context().get("fields").map(String::as_str),
            Some("ploidy"),
            "el diagnóstico tiene que identificar el campo, ploidía {illegal}"
        );
        assert_eq!(
            rejection.context().get("detail").map(String::as_str),
            Some("debe pertenecer a 1..=10"),
            "el diagnóstico tiene que identificar el problema, ploidía {illegal}"
        );
    }

    fs::remove_dir_all(&root).expect("limpieza");
}

// ---------------------------------------------------------------------------
// D. Búsquedas negativas de alcance, con control positivo
// ---------------------------------------------------------------------------

/// Las cinco prohibiciones de alcance de F5, buscadas sobre **código**.
///
/// Panel, variantes estructurales, repeticiones, alineamientos BAM/CRAM y
/// recombinación pertenecen a F6 y F7. La guarda es léxica y por eso necesita
/// dos controles: uno que demuestre que la búsqueda **encuentra** el término
/// cuando está en el código, y otro que demuestre que **no** lo encuentra
/// cuando solo aparece en la prosa que lo prohíbe —que es justo lo que los
/// módulos F5 hacen—, porque si no una búsqueda ingenua sería siempre roja.
#[test]
fn no_f5_source_or_schema_names_a_panel_sv_repeat_alignment_or_recombination() {
    let forbidden = [
        "panel",
        "structural_variant",
        "tandem_repeat",
        "microsatellite",
        "bam",
        "cram",
        "recombination",
        "recombinant",
        "crossover",
        "meiosis",
        "mosaic",
        "linkage_map",
    ];

    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut inspected = 0_usize;
    for entry in fs::read_dir(root.join("src")).expect("listar src") {
        let path = entry.expect("entrada de src").path();
        let name = path
            .file_name()
            .expect("nombre")
            .to_string_lossy()
            .into_owned();
        if !name.starts_with("f5_") || !name.ends_with(".rs") {
            continue;
        }
        inspected += 1;
        let code = strip_comments(&fs::read_to_string(&path).expect("leer módulo F5"));
        for needle in forbidden {
            assert!(
                !code.to_ascii_lowercase().contains(needle),
                "{name} nombra «{needle}» en código, y F5 no tiene ese alcance"
            );
        }
    }
    assert!(
        inspected >= 10,
        "la búsqueda tiene que recorrer los módulos F5, y solo vio {inspected}"
    );

    for schema in [
        "schemas/config-f5-scenario-v1.schema.json",
        "schemas/effective-config-v5.schema.json",
        "schemas/manifest-general-v5.schema.json",
        "schemas/f5-failure-v1.schema.json",
        "schemas/profile-metadata-f5-v1.schema.json",
        "schemas/truth-v2.schema.json",
    ] {
        let document = fs::read_to_string(root.join(schema))
            .unwrap_or_else(|error| panic!("leer {schema}: {error}"))
            .to_ascii_lowercase();
        for needle in forbidden {
            assert!(
                !document.contains(needle),
                "{schema} nombra «{needle}», y F5 no tiene ese alcance"
            );
        }
    }

    // Control positivo 1: la búsqueda encuentra el término en el código.
    let planted = strip_comments("fn build_panel() -> Panel { todo() }");
    assert!(planted.to_ascii_lowercase().contains("panel"));
    // Control positivo 2: y **no** lo encuentra cuando solo está en la prosa,
    // que es el modo real en que los módulos F5 nombran lo que prohíben.
    let prose = strip_comments("/// F5 no acepta ningún panel externo.\nfn build() {}");
    assert!(!prose.to_ascii_lowercase().contains("panel"));
    let inline = strip_comments("fn build() {} // sin panel ni mapa de recombinación");
    assert!(!inline.to_ascii_lowercase().contains("panel"));
    assert!(!inline.to_ascii_lowercase().contains("recombin"));
}

/// Quita los comentarios `//` de una fuente Rust sin tocar su código.
///
/// No interpreta literales de cadena: no hace falta, porque lo que la guarda
/// busca son identificadores y ningún módulo F5 construye cadenas con los
/// términos prohibidos. Lo que sí hace falta es no confundir una prohibición
/// escrita en prosa con una capacidad implementada.
fn strip_comments(source: &str) -> String {
    source
        .lines()
        .map(|line| match line.find("//") {
            Some(index) => &line[..index],
            None => line,
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// Ninguna vía de publicación de F5.6 llama a los publicadores desnudos.
///
/// F5.6.1 fijó la obligación sobre `application.rs` y `main.rs`, que son las
/// dos piezas que entonces publicaban. El cierre la extiende a **todo** el
/// crate menos el módulo que envuelve: cualquier módulo nuevo que publicase
/// por `publish_f5_success` o `publish_f5_failure` obtendría un destino sin
/// testigo y dejaría la clausura inaplicable sin ninguna señal.
///
/// La declaración de las dos funciones en `f5_publication.rs` no es una
/// llamada y no cuenta: la búsqueda descarta lo que va precedido de `fn `.
#[test]
fn no_f5_6_publication_path_calls_the_bare_publishers() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut inspected = 0_usize;
    for entry in fs::read_dir(root.join("src")).expect("listar src") {
        let path = entry.expect("entrada de src").path();
        let name = path
            .file_name()
            .expect("nombre")
            .to_string_lossy()
            .into_owned();
        if !name.ends_with(".rs") || name == "f5_closure.rs" {
            continue;
        }
        inspected += 1;
        let source = fs::read_to_string(&path).expect("leer módulo");
        for forbidden in ["publish_f5_success(", "publish_f5_failure("] {
            assert!(
                bare_call_positions(&source, forbidden).is_empty(),
                "{name} publica por {forbidden} fuera de los envoltorios que sellan"
            );
        }
    }
    assert!(
        inspected >= 40,
        "la búsqueda tiene que recorrer el crate, y solo vio {inspected} módulos"
    );

    // Control positivo 1: el orquestador sí publica, y por los envoltorios.
    let orchestrator = fs::read_to_string(root.join("src/application.rs")).expect("orquestador");
    assert!(orchestrator.contains("seal_and_publish_f5_success("));
    assert!(orchestrator.contains("seal_and_publish_f5_failure("));
    // Control positivo 2: la misma búsqueda encuentra la llamada desnuda
    // donde realmente vive, que es el módulo excluido.
    let wrappers = fs::read_to_string(root.join("src/f5_closure.rs")).expect("envoltorios");
    for forbidden in ["publish_f5_success(", "publish_f5_failure("] {
        assert!(
            !bare_call_positions(&wrappers, forbidden).is_empty(),
            "la búsqueda no discrimina si no encuentra {forbidden} en f5_closure.rs"
        );
    }
    // Control positivo 3: y distingue el envuelto del envoltorio, y la
    // declaración de la llamada.
    assert!(
        bare_call_positions("seal_and_publish_f5_success(x)", "publish_f5_success(").is_empty()
    );
    assert!(bare_call_positions("pub fn publish_f5_success(x)", "publish_f5_success(").is_empty());
    assert_eq!(
        bare_call_positions("publish_f5_success(x)", "publish_f5_success(").len(),
        1
    );
}

/// Posiciones de `needle` que son una **llamada** desnuda.
///
/// Descarta las precedidas por `seal_and_` —el nombre del envoltorio contiene
/// el del envuelto, de modo que una búsqueda ingenua daría siempre positivo— y
/// las precedidas por `fn `, que son la declaración y no una llamada.
fn bare_call_positions(source: &str, needle: &str) -> Vec<usize> {
    source
        .match_indices(needle)
        .filter(|(index, _)| {
            let before = &source[..*index];
            !before.ends_with("seal_and_") && !before.ends_with("fn ")
        })
        .map(|(index, _)| index)
        .collect()
}

// ---------------------------------------------------------------------------
// E. Los dos auditores sobre el destino ya renombrado: la mitad «fallo»
// ---------------------------------------------------------------------------

/// El auditor de paquete acepta el destino **renombrado** de un fallo.
///
/// F5.6.1 compuso `audit_f5_package` y `audit_f5_identity_closure` en la mitad
/// «éxito» y `audit_f5_failure_closure` en la de «fallo», y dejó dicho que
/// `audit_f5_failure_package` necesita una `F5FailureDeclaration` que el
/// orquestador no conserva. Aquí se cierra esa mitad: el fallo lo publica la
/// cadena productiva real —una cuota controlada que su tolerancia declarada no
/// sostiene, sin ninguna instrumentación—, la declaración se reconstruye desde
/// el **manifiesto** y el auditor relee `outputs/f5-failure.json`, que es otro
/// documento, sobre el destino ya renombrado. Las tres mutaciones demuestran
/// que discrimina.
#[test]
fn the_failure_package_auditor_accepts_the_renamed_destination_and_dies_of_each_mutation() {
    let root = temporary("failure-audit");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_infeasible_quota_pair(&root);

    let destination = root.join("package");
    let rejection = execute_validated_f5_config(&admitted(&general), &destination)
        .expect_err("una cuota imposible no puede publicar un éxito");
    assert_eq!(rejection.code(), "GPS-DOMAIN-103");

    // El destino está renombrado: no queda ningún provisional vivo.
    assert!(destination.is_dir());
    assert!(!destination.join("INCOMPLETE").exists());
    assert!(no_provisionals(&destination));

    let declaration = failure_declaration_from_manifest(&destination);
    assert_eq!(declaration.code, "GPS-DOMAIN-103");
    assert_eq!(declaration.stage, "founders");
    audit_f5_failure_package(&destination, &declaration)
        .expect("el auditor tiene que aceptar el paquete fallido recién renombrado");

    // Mutación 1: la etapa declarada deja de ser la publicada.
    let mut mutated = declaration.clone();
    mutated.stage = "layout".to_owned();
    assert_ne!(
        mutated.stage, declaration.stage,
        "la mutación tiene que mutar"
    );
    audit_f5_failure_package(&destination, &mutated)
        .expect_err("una etapa que el documento no sostiene tiene que rechazarse");

    // Mutación 2: el enlace deja de describir sus bytes.
    let mut mutated = declaration.clone();
    mutated.link.size_bytes += 1;
    assert_ne!(
        mutated.link.size_bytes, declaration.link.size_bytes,
        "la mutación tiene que mutar"
    );
    audit_f5_failure_package(&destination, &mutated)
        .expect_err("un enlace que no describe sus bytes tiene que rechazarse");

    // Mutación 3: un artefacto biológico aparece en el destino publicado.
    let intruder = destination.join("outputs/f5-founders.tsv");
    assert!(
        !intruder.exists(),
        "el fallo no publica artefactos biológicos"
    );
    fs::write(&intruder, b"copy_id\n").expect("intruso");
    audit_f5_failure_package(&destination, &declaration)
        .expect_err("un artefacto biológico en un paquete fallido tiene que rechazarse");
    fs::remove_file(&intruder).expect("retirar el intruso");
    audit_f5_failure_package(&destination, &declaration)
        .expect("retirado el intruso, el paquete vuelve a ser el que era");

    fs::remove_dir_all(&root).expect("limpieza");
}

/// Regresión del hallazgo F5.6.3: el auditor aceptaba la mitad de los fallos.
///
/// `audit_f5_failure_package` comparaba los tres contadores contra `limits`
/// **posición a posición**, de modo que a `swap_proposals` le aplicaba
/// `max_swap_proposals`. DTC-089 l. 1196 fija ese límite como presupuesto de
/// **un intento** y l. 1847-1848 enumera `swap_proposals` junto a
/// `haplotype_attempts_used`, que es acumulado: el contador vale
/// `haplotype_attempts_used * max_swap_proposals` en cuanto la búsqueda se
/// agota. Con el par de esta batería, `32 * 256 = 8192` frente a un límite de
/// `256`, y el auditor rechazaba con `GPS-RES-101` un paquete que la propia
/// cadena productiva acababa de publicar.
///
/// El defecto era de F5.4.2 y solo podía aparecer al componer este auditor
/// sobre un agotamiento **real**: ninguna revisión por hito lo vio porque las
/// baterías de F5.4 construían sus declaraciones con contadores pequeños.
///
/// La prueba fija las dos mitades. La cota tiene que ser el producto: ni la
/// vieja —`max_swap_proposals` a secas, que rechaza este paquete— ni una
/// demasiado laxa, porque una cota que no rechaza nada pasaría igual.
#[test]
fn the_failure_package_auditor_accepts_an_exhausted_search_package() {
    let root = temporary("exhausted-audit");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_pair(&root);
    make_ld_target_infeasible(&root);

    let destination = root.join("package");
    let rejection = execute_validated_f5_config(&admitted(&general), &destination)
        .expect_err("la búsqueda agotada no puede publicar un éxito");
    assert_eq!(rejection.code(), "GPS-DOMAIN-104");

    let declaration = failure_declaration_from_manifest(&destination);
    assert_eq!(declaration.stage, "haplotypes");
    let attempts = declaration.counters[1];
    let per_attempt = declaration.limits[2];
    let accumulated = declaration.counters[2];
    // El escenario **alcanza** la cota vieja y la supera: si no, la regresión
    // no discriminaría nada. Es el agotamiento completo, un intento tras otro
    // consumiendo su presupuesto entero.
    assert_eq!(attempts, declaration.limits[1]);
    assert_eq!(accumulated, attempts * per_attempt);
    assert!(
        accumulated > per_attempt,
        "sin superar la cota vieja esta prueba no sería una regresión"
    );
    audit_f5_failure_package(&destination, &declaration)
        .expect("el auditor tiene que aceptar el paquete que la cadena publicó");

    // Y la cota sigue **cerrando**: una propuesta más de las que los intentos
    // publicados pueden haber hecho se rechaza, y por su propia causa.
    //
    // La mutación va sobre los **bytes del paquete**, no sobre la declaración:
    // el auditor compara antes declaración y documento término a término, de
    // modo que tocar solo la declaración moriría de esa otra causa y no diría
    // nada de la cota. Aquí los dos documentos siguen de acuerdo entre sí y lo
    // único que falla es el contrato de recursos.
    let mutated = mutate_published_proposals(&destination, accumulated + 1);
    assert_eq!(mutated.counters[2], accumulated + 1);
    let rejection = audit_f5_failure_package(&destination, &mutated)
        .expect_err("una propuesta por encima de intentos × presupuesto tiene que rechazarse");
    let diagnostic = rejection.into_diagnostic();
    assert_eq!(diagnostic.code(), "GPS-RES-101");
    let detail = diagnostic
        .context()
        .get("detail")
        .expect("detalle del rechazo")
        .clone();
    assert!(detail.contains("swap_proposals"), "{detail}");
    assert!(detail.contains("haplotype_attempts_used"), "{detail}");

    fs::remove_dir_all(&root).expect("limpieza");
}

/// El orquestador compone `audit_f5_failure_package`, y se demuestra por el
/// diagnóstico que devuelve, no solo por el texto de la fuente.
///
/// Era la mitad que F5.6.1 no pudo componer. Con la cota corregida ya es
/// posible: el orquestador reconstruye la declaración desde el **manifiesto
/// publicado** —no desde su informe en memoria, que haría la auditoría
/// autocertificante— y se la da al oráculo sobre el destino ya renombrado.
///
/// La prueba discrimina porque el agotamiento de búsqueda es exactamente el
/// caso que la cota vieja rechazaba: si el auditor **no** estuviera compuesto,
/// devolver la cota a `max_swap_proposals` no cambiaría nada; compuesto,
/// convierte `GPS-DOMAIN-104` en `GPS-RES-101`. Comprobado con esa mutación.
#[test]
fn the_orchestrator_composes_the_failure_package_auditor() {
    let root = temporary("composed-failure-audit");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_pair(&root);
    make_ld_target_infeasible(&root);

    let destination = root.join("package");
    let returned = execute_validated_f5_config(&admitted(&general), &destination)
        .expect_err("la búsqueda agotada no puede publicar un éxito");
    // El diagnóstico es el del informe y **no** el de ningún auditor: los dos
    // gemelos corrieron sobre el destino renombrado y lo aceptaron.
    assert_eq!(returned.code(), "GPS-DOMAIN-104");
    assert_ne!(returned.code(), "GPS-RES-101");

    // Y el paquete que aceptaron es el que está en el disco.
    let manifest: Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifiesto"))
            .expect("JSON");
    assert_eq!(manifest["state"], "failed");
    assert_eq!(manifest["f5"]["f5_failure"]["stage"], "haplotypes");

    // Control léxico, con su control positivo: el orquestador nombra los dos
    // gemelos de la mitad «fallo» y sigue nombrando los dos de la de «éxito».
    let orchestrator =
        fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("src/application.rs"))
            .expect("orquestador");
    for auditor in [
        "audit_f5_failure_package(",
        "audit_f5_failure_closure(",
        "audit_f5_package(",
        "audit_f5_identity_closure(",
    ] {
        assert!(
            orchestrator.contains(auditor),
            "el orquestador tiene que componer {auditor}"
        );
    }
    assert!(!orchestrator.contains("audit_f5_failure_recombination("));

    fs::remove_dir_all(&root).expect("limpieza");
}

/// Reescribe `swap_proposals` en los **dos** documentos del paquete fallido.
///
/// Rehace también el tamaño y el SHA-256 del enlace, de modo que el paquete
/// mutado siga siendo coherente consigo mismo y la única regla que pueda
/// romper sea la cota de propuestas.
fn mutate_published_proposals(destination: &Path, proposals: u64) -> F5FailureDeclaration {
    let failure_path = destination.join("outputs/f5-failure.json");
    let mut failure: Value =
        serde_json::from_slice(&fs::read(&failure_path).expect("documento de fallo"))
            .expect("JSON del fallo");
    let previous = failure["counters"]["swap_proposals"].clone();
    failure["counters"]["swap_proposals"] = Value::from(proposals);
    assert_ne!(
        previous, failure["counters"]["swap_proposals"],
        "la mutación tiene que mutar"
    );
    let bytes = serde_json::to_vec(&failure).expect("serializar el fallo mutado");
    fs::write(&failure_path, &bytes).expect("documento de fallo mutado");

    let manifest_path = destination.join("manifest.json");
    let mut manifest: Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("manifiesto")).expect("JSON");
    let failure_node = &mut manifest["f5"]["f5_failure"];
    failure_node["counters"]["swap_proposals"] = Value::from(proposals);
    failure_node["artifacts"][0]["size_bytes"] = Value::from(bytes.len() as u64);
    failure_node["artifacts"][0]["sha256"] = Value::from(sha256(&bytes));
    fs::write(
        &manifest_path,
        serde_json::to_vec(&manifest).expect("serializar el manifiesto mutado"),
    )
    .expect("manifiesto mutado");

    failure_declaration_from_manifest(destination)
}

/// Reconstruye la declaración del fallo desde el manifiesto publicado.
fn failure_declaration_from_manifest(destination: &Path) -> F5FailureDeclaration {
    let manifest: Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifiesto"))
            .expect("manifiesto JSON");
    let failure = &manifest["f5"]["f5_failure"];
    let link = &failure["artifacts"][0];
    let text = |value: &Value| value.as_str().expect("cadena").to_owned();
    let optional = |value: &Value| value.as_str().map(str::to_owned);
    let numbers = |object: &Value, keys: &[&str]| {
        keys.iter()
            .map(|key| object[*key].as_u64().expect("entero"))
            .collect::<Vec<_>>()
    };
    F5FailureDeclaration {
        code: text(&failure["code"]),
        stage: text(&failure["stage"]),
        scenario_id: text(&failure["scenario_id"]),
        population_id: optional(&failure["population_id"]),
        block_id: optional(&failure["block_id"]),
        limits: numbers(
            &failure["limits"],
            &[
                "max_layout_attempts",
                "max_haplotype_attempts",
                "max_swap_proposals",
                "max_haplotype_cells",
                "max_operation_units",
                "max_working_memory_bytes",
                "max_output_bytes",
            ],
        )
        .try_into()
        .expect("siete límites"),
        counters: numbers(
            &failure["counters"],
            &[
                "layout_attempts_used",
                "haplotype_attempts_used",
                "swap_proposals",
                "swap_accepted",
            ],
        )
        .try_into()
        .expect("cuatro contadores"),
        best_layout: failure["best_layout"].clone(),
        best_attempt: failure["best_attempt"].clone(),
        link: F5DeclaredLink {
            artifact_id: text(&link["artifact_id"]),
            artifact_kind: text(&link["artifact_kind"]),
            schema_version: text(&link["schema_version"]),
            path: text(&link["path"]),
            size_bytes: link["size_bytes"].as_u64().expect("entero"),
            sha256: text(&link["sha256"]),
        },
    }
}

// ---------------------------------------------------------------------------
// F. Atomicidad en frontera de proceso — la deuda que F5.6.1 dejó abierta
// ---------------------------------------------------------------------------

/// Una ejecución F5 interrumpida de verdad no deja destino.
///
/// Hasta este hito la atomicidad de F5 estaba demostrada en frontera de
/// **biblioteca**: `an_interrupted_execution_leaves_the_provisional_marked_incomplete`
/// simula la interrupción dejando caer el provisional dentro del proceso de
/// prueba. DTC-045 exige la otra: un proceso real, muerto por `SIGKILL` a
/// mitad de la generación, que no puede ejecutar ningún código de limpieza.
/// Lo que queda tiene que ser un provisional marcado `INCOMPLETE` con permisos
/// `0700` y **ningún** destino.
#[test]
fn an_interrupted_f5_run_leaves_its_provisional_and_creates_no_destination() {
    let root = temporary("interrupted");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_slow_pair(&root);
    let destination = root.join("package");

    let mut child = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run"])
        .arg(&general)
        .arg("--output")
        .arg(&destination)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .expect("proceso F5 interrumpible");

    let provisional = wait_for_provisionals(&destination, 1).remove(0);
    assert!(
        child.try_wait().expect("estado del hijo").is_none(),
        "el escenario lento tiene que seguir vivo cuando se le interrumpe"
    );
    let kill = Command::new("kill")
        .args(["-KILL", &child.id().to_string()])
        .status()
        .expect("enviar SIGKILL");
    assert!(kill.success());
    let status = child.wait().expect("recoger el hijo muerto");
    assert!(!status.success(), "un proceso matado no termina bien");

    assert!(
        !destination.exists(),
        "una ejecución interrumpida no puede dejar destino publicado"
    );
    assert!(provisional.join("INCOMPLETE").is_file());
    assert_eq!(
        fs::metadata(&provisional)
            .expect("permisos del provisional")
            .permissions()
            .mode()
            & 0o777,
        0o700
    );
    // Y lo que quedó a medias vive **dentro** del provisional, no en el sitio
    // del destino: ningún byte biológico escapó del árbol marcado.
    for biological in ["outputs/truth.jsonl", "outputs/f5-founders.tsv"] {
        assert!(!destination.join(biological).exists());
    }

    fs::remove_dir_all(&root).expect("limpieza");
}

/// Dos procesos F5 concurrentes tienen un solo ganador.
///
/// La concurrencia **entre procesos** que DTC-045 exige solo recorría F1.
/// Aquí son dos `run` F5 reales sobre el mismo destino, solapados de verdad:
/// el escenario dura lo suficiente para que los dos abran su provisional antes
/// de que ninguno renombre. Exactamente uno publica; el otro se rechaza con
/// `GPS-IO-003` y conserva su provisional marcado, sin tocar el destino ajeno.
/// La primera mitad fija además que un destino que ya existe no se machaca.
#[test]
fn two_concurrent_f5_runs_have_one_winner_and_never_clobber_the_destination() {
    let root = temporary("race");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_slow_pair(&root);

    // Mitad 1: un destino que ya existe se rechaza sin tocar su contenido.
    let existing = root.join("existing");
    fs::create_dir_all(&existing).expect("destino existente");
    fs::write(existing.join("sentinel"), b"preserve").expect("centinela");
    let refused = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .args(["run"])
        .arg(&general)
        .arg("--output")
        .arg(&existing)
        .output()
        .expect("run F5 sobre un destino existente");
    assert!(!refused.status.success());
    assert!(
        String::from_utf8_lossy(&refused.stderr).contains("GPS-IO-003"),
        "{}",
        String::from_utf8_lossy(&refused.stderr)
    );
    assert_eq!(
        fs::read(existing.join("sentinel")).expect("centinela"),
        b"preserve"
    );

    // Mitad 2: dos procesos solapados sobre el mismo destino libre.
    let destination = root.join("winner");
    let spawn = || {
        Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
            .args(["run"])
            .arg(&general)
            .arg("--output")
            .arg(&destination)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("competidor F5")
    };
    let first = spawn();
    let second = spawn();
    // El solape es real y se comprueba: los dos provisionales coexisten antes
    // de que ninguno de los dos haya renombrado.
    let provisionals = wait_for_provisionals(&destination, 2);
    assert_eq!(provisionals.len(), 2);
    assert!(
        !destination.exists(),
        "los dos compiten antes de que exista el destino"
    );

    let outputs = [
        first.wait_with_output().expect("primer competidor"),
        second.wait_with_output().expect("segundo competidor"),
    ];
    assert_eq!(
        outputs
            .iter()
            .filter(|output| output.status.success())
            .count(),
        1,
        "exactamente un proceso publica"
    );
    let loser = outputs
        .iter()
        .find(|output| !output.status.success())
        .expect("un perdedor");
    assert!(
        String::from_utf8_lossy(&loser.stderr).contains("GPS-IO-003"),
        "{}",
        String::from_utf8_lossy(&loser.stderr)
    );

    // El ganador publicó un paquete completo y el perdedor no lo tocó.
    assert!(destination.join("manifest.json").is_file());
    assert!(destination.join("outputs/truth.jsonl").is_file());
    assert!(!destination.join("INCOMPLETE").exists());
    let surviving = provisionals
        .iter()
        .filter(|path| path.exists())
        .collect::<Vec<_>>();
    assert_eq!(
        surviving.len(),
        1,
        "el perdedor conserva su provisional y el ganador renombró el suyo"
    );
    assert!(surviving[0].join("INCOMPLETE").is_file());

    fs::remove_dir_all(&root).expect("limpieza");
}

/// Espera a que existan al menos `count` provisionales marcados del destino.
fn wait_for_provisionals(destination: &Path, count: usize) -> Vec<PathBuf> {
    let deadline = Instant::now() + Duration::from_secs(30);
    loop {
        let mut paths = provisionals(destination);
        paths.sort();
        if paths.len() >= count {
            return paths;
        }
        assert!(
            Instant::now() < deadline,
            "expiró la espera de {count} provisionales de {}",
            destination.display()
        );
        std::thread::sleep(Duration::from_millis(10));
    }
}

fn provisionals(destination: &Path) -> Vec<PathBuf> {
    let parent = destination.parent().expect("directorio padre del destino");
    let prefix = format!(
        ".{}.gps-incomplete-",
        destination
            .file_name()
            .expect("nombre del destino")
            .to_string_lossy()
    );
    fs::read_dir(parent)
        .expect("listar el directorio padre")
        .filter_map(Result::ok)
        .filter(|entry| entry.file_name().to_string_lossy().starts_with(&prefix))
        .map(|entry| entry.path())
        .filter(|path| path.join("INCOMPLETE").is_file())
        .collect()
}

fn no_provisionals(destination: &Path) -> bool {
    provisionals(destination).is_empty()
}

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------

fn sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-6-closure-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

fn admitted(general: &Path) -> Box<EffectiveF5Config> {
    match load_and_validate_cli_config(general).expect("par F5 admitido") {
        CliValidatedConfig::F5(effective) => effective,
        CliValidatedConfig::Legacy(_) => panic!("el dispatcher tenía que clasificar un general v3"),
    }
}

const SEED: &str = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8";

/// Escribe una autoría F5 factible mínima y devuelve la ruta del general.
fn write_pair(root: &Path) -> PathBuf {
    write_pair_sized(root, 2, 1, ">chr1\nACGTAC\n")
}

/// La misma autoría, dimensionada para que el `run` dure segundos.
///
/// La interrupción y la carrera necesitan un proceso que siga vivo cuando la
/// prueba actúa sobre él. Con el par mínimo las dos pasan —el `run` en
/// depuración tarda unos milisegundos y el provisional se abre antes del
/// layout, de modo que la ventana existe—, pero el margen depende de la carga
/// de la máquina. Con este par el `run` dura unos segundos y la ventana deja
/// de ser una coincidencia de reloj. Que la ventana se aprovecha de verdad no
/// se supone: la interrupción exige que el hijo siga vivo y la carrera exige
/// que los **dos** provisionales coexistan antes de que ninguno renombre.
fn write_slow_pair(root: &Path) -> PathBuf {
    let sequence = "ACGT".repeat(75_000);
    write_pair_sized(root, 192, 4, &format!(">chr1\n{sequence}\n"))
}

fn write_pair_sized(root: &Path, blocks: u32, individuals: u32, fasta: &str) -> PathBuf {
    let scenario_directory = root.join("scenarios");
    let profile_directory = scenario_directory.join("profiles");
    fs::create_dir_all(&profile_directory).expect("directorio de perfil");
    fs::create_dir_all(root.join("references")).expect("directorio de referencia");

    let files = [
        (
            "block_lengths.tsv",
            "f5_block_length_profile",
            "length_min_bp\tlength_max_bp\tweight\n3\t3\t1\n",
        ),
        (
            "snp_gaps.tsv",
            "f5_snp_gap_profile",
            "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n",
        ),
        (
            "snp_spectrum.tsv",
            "f5_snp_spectrum",
            "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t1\n",
        ),
        (
            "haplotype_weights.tsv",
            "f5_haplotype_weights",
            "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        ),
        (
            "ld_profile.tsv",
            "f5_ld_profile",
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t1\t0.05\nP1\t2\t2\t1\t0.05\n",
        ),
    ];
    let links = files
        .iter()
        .map(|(name, kind, contents)| {
            fs::write(profile_directory.join(name), contents).expect("TSV de perfil");
            format!(
                r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#,
                contents.len(),
                sha256(contents.as_bytes())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let metadata = format!(
        r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
    );
    fs::write(profile_directory.join("profile.metadata.json"), &metadata).expect("metadata");

    fs::write(root.join("references/ref.fa"), fasta).expect("referencia");
    let general_path = root.join("general.toml");
    fs::write(
        &general_path,
        format!(
            "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{SEED}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
            sha256(fasta.as_bytes())
        ),
    )
    .expect("general");
    let founders = (1..=individuals)
        .map(|index| format!("\"I{index}\""))
        .collect::<Vec<_>>()
        .join(", ");
    fs::write(
        scenario_directory.join("scenario.toml"),
        format!(
            "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = {individuals}\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = 32\nmax_swap_proposals = 256\nmax_haplotype_cells = 65536\nblock_count = {blocks}\n\n[[subpopulations]]\npopulation_id = \"P1\"\nindividual_count = {individuals}\nfounder_ids = [{founders}]\n",
            sha256(metadata.as_bytes())
        ),
    )
    .expect("escenario");
    general_path
}

/// Una autoría cuya **cuota controlada** no cabe en su tolerancia declarada.
///
/// Es la vía productiva que falla en la etapa `founders` sin instrumentar
/// nada: la búsqueda resuelve, los fundadores no. Con dos haplotipos de pesos
/// `0,6` y `0,4` y `C = 2` copias, el reparto entero es `[1, 1]` y el rango 0
/// se desvía `0,2` copias, que una tolerancia declarada de `0` no admite.
fn write_infeasible_quota_pair(root: &Path) -> PathBuf {
    let general = write_pair(root);
    let directory = root.join("scenarios/profiles");
    let weights = "haplotype_rank\tweight\n0\t0.6\n1\t0.4\n";
    let previous = fs::read_to_string(directory.join("haplotype_weights.tsv")).expect("pesos");
    assert_ne!(previous, weights, "la mutación tiene que mutar");
    fs::write(directory.join("haplotype_weights.tsv"), weights).expect("pesos asimétricos");
    rewrite_profile_metadata(root);

    let scenario_path = root.join("scenarios/scenario.toml");
    let scenario = fs::read_to_string(&scenario_path).expect("escenario");
    let rewritten = scenario.replace(
        "individual_count = 1\nfounder_ids",
        "individual_count = 1\nsampling = \"controlled\"\nfounder_ids",
    ) + "\n[[sampling_overrides]]\npopulation_id = \"P1\"\nblock_index = 0\nmode = \"controlled\"\ncontrolled_max_abs_error = 0.0\n";
    assert_ne!(rewritten, scenario, "la reescritura tiene que reescribir");
    fs::write(&scenario_path, rewritten).expect("escenario controlado");
    general
}

/// Deja el objetivo de LD fuera del alcance de la búsqueda.
fn make_ld_target_infeasible(root: &Path) {
    let profile = root.join("scenarios/profiles/ld_profile.tsv");
    let infeasible = "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\nP1\t2\t2\t0.5\t0.05\n";
    let previous = fs::read_to_string(&profile).expect("perfil LD");
    assert_ne!(previous, infeasible, "la mutación tiene que mutar");
    fs::write(&profile, infeasible).expect("perfil LD inviable");
    rewrite_profile_metadata(root);
}

/// Rehace el metadata del perfil y el hash que el escenario declara.
fn rewrite_profile_metadata(root: &Path) {
    let directory = root.join("scenarios/profiles");
    let kinds = [
        ("block_lengths.tsv", "f5_block_length_profile"),
        ("snp_gaps.tsv", "f5_snp_gap_profile"),
        ("snp_spectrum.tsv", "f5_snp_spectrum"),
        ("haplotype_weights.tsv", "f5_haplotype_weights"),
        ("ld_profile.tsv", "f5_ld_profile"),
    ];
    let links = kinds
        .iter()
        .map(|(name, kind)| {
            let contents = fs::read(directory.join(name)).expect("TSV de perfil");
            format!(
                r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#,
                contents.len(),
                sha256(&contents)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let metadata = format!(
        r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
    );
    fs::write(directory.join("profile.metadata.json"), &metadata).expect("metadata");

    let scenario_path = root.join("scenarios/scenario.toml");
    let scenario = fs::read_to_string(&scenario_path).expect("escenario");
    let previous = scenario
        .lines()
        .find(|line| line.starts_with("profile_metadata_sha256"))
        .expect("hash del perfil")
        .to_owned();
    let rewritten = scenario.replace(
        &previous,
        &format!(
            "profile_metadata_sha256 = \"{}\"",
            sha256(metadata.as_bytes())
        ),
    );
    assert_ne!(rewritten, scenario, "la reescritura tiene que reescribir");
    fs::write(&scenario_path, rewritten).expect("escenario reescrito");
}
