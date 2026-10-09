//! Cada rechazo del juez semántico de la verdad, con su motivo exacto.
//!
//! DT-001, tanda 5. El juez (`src/truth_validation.rs`) tenía 118
//! comprobaciones de las que 93 no las disparaba ninguna prueba. Aquí cada
//! caso parte de una verdad **publicada de verdad** por el programa (una
//! referencia de 600 bases con un tracto `(CA)50`, tres individuos de dos
//! poblaciones, una matriz de SNV con un locus multialélico, tres cambios de
//! base directos o de matriz, dos inserciones y dos deleciones directas, una
//! de ellas compuesta con una inserción anterior), la altera **en un solo
//! punto** y exige el motivo exacto del rechazo, con los dos lectores —el
//! completo y el de flujo— cuando los dos pueden verlo.
//!
//! # Lo que este fichero no es
//!
//! Las guardas que ningún `truth.jsonl` publicado y manipulado puede alcanzar
//! están anotadas **junto a la guarda** en `src/truth_validation.rs`, con el
//! porqué. Las pruebas de este fichero que solo «ejecutan el intento» sobre
//! ellas lo dicen en su docstring: **no discriminan** esas guardas.

use std::io::{BufRead, Cursor, Read};

use genome_population_simulator::reference::{Reference, ReferenceOrigin, SequenceClass};
use genome_population_simulator::truth_schema::{
    GenotypeCaseKind, TruthAllele, TruthContractError, TruthCurrentPiece, TruthHistoryEvent,
    TruthModification, TruthRepeatBlock, TruthSegment, TruthSourceKind,
};
use genome_population_simulator::truth_validation::{
    TruthValidationError, read_truth_jsonl, validate_and_reconstruct_truth,
    validate_and_reconstruct_truth_jsonl, validate_truth_jsonl_streaming,
};
use serde_json::{Value, json};

mod truth_guard_scenario;
use truth_guard_scenario::*;

// ---------------------------------------------------------------------------
// Control positivo
// ---------------------------------------------------------------------------

/// Sin alterar nada, la verdad publicada la aceptan los dos lectores. Sin este
/// control, un juez que rechazara todo satisfaría cada caso de abajo.
#[test]
fn the_published_truth_is_accepted_by_both_readers() {
    let text = jsonl(&lines());
    for result in both(text.as_bytes()) {
        result.expect("la verdad publicada intacta se acepta");
    }
}

// ---------------------------------------------------------------------------
// El lector de líneas: lo que se rechaza antes de mirar el contenido
// ---------------------------------------------------------------------------

fn schema_error(line: u64, detail: &str) -> TruthValidationError {
    TruthValidationError::Schema {
        line,
        detail: detail.to_owned(),
    }
}

/// Un byte que no es UTF-8 dentro de la línea 2 (la primera copia).
fn with_invalid_utf8_on_line_two() -> (Vec<u8>, usize) {
    let text = bytes_of(&lines());
    let first_newline = text
        .iter()
        .position(|byte| *byte == b'\n')
        .expect("línea 1");
    let marker = b"haplotype_copy";
    let start = text[first_newline + 1..]
        .windows(marker.len())
        .position(|window| window == marker)
        .expect("marca");
    let absolute = first_newline + 1 + start;
    let mut corrupt = text.clone();
    corrupt[absolute] = 0xFF;
    (corrupt, absolute - (first_newline + 1))
}

/// El texto que no es UTF-8 se rechaza con el desplazamiento exacto, en la
/// línea 2, y los dos lectores lo dicen igual.
#[test]
fn a_line_that_is_not_utf8_is_rejected_with_its_line_and_offset() {
    let (corrupt, offset) = with_invalid_utf8_on_line_two();
    assert_both(
        "UTF-8 inválido",
        &corrupt,
        &schema_error(
            2,
            &format!("invalid utf-8 sequence of 1 bytes from index {offset}"),
        ),
    );
}

fn with_blank_line_two(terminator: &str) -> Vec<u8> {
    let text = jsonl(&lines());
    let first_newline = text.find('\n').expect("línea 1");
    let mut corrupt = String::new();
    corrupt.push_str(&text[..=first_newline]);
    corrupt.push_str(terminator);
    corrupt.push_str(&text[first_newline + 1..]);
    corrupt.into_bytes()
}

/// Una línea vacía entre dos registros se rechaza, con `\n` y con `\r\n`: el
/// retorno de carro no cuenta como contenido.
#[test]
fn an_empty_line_is_rejected_with_either_line_terminator() {
    for terminator in ["\n", "\r\n"] {
        assert_both(
            &format!("línea vacía {terminator:?}"),
            &with_blank_line_two(terminator),
            &schema_error(2, "línea JSONL vacía"),
        );
    }
}

/// Un documento con terminadores `\r\n` en todas sus líneas se acepta.
#[test]
fn a_document_with_carriage_returns_is_accepted() {
    let text = jsonl(&lines()).replace('\n', "\r\n");
    for result in both(text.as_bytes()) {
        result.expect("`\\r\\n` es un terminador admitido");
    }
}

/// La primera línea es una copia, no la cabecera: los dos lectores lo rechazan.
#[test]
fn a_copy_before_the_header_is_rejected() {
    let mut all = lines();
    all.swap(0, 1);
    assert_both(
        "copia antes de la cabecera",
        &bytes_of(&all),
        &TruthValidationError::HeaderMustBeFirst,
    );
}

/// La primera línea es un caso genotípico, no la cabecera.
#[test]
fn a_genotype_before_the_header_is_rejected() {
    let mut all = lines();
    let genotype = all.pop().expect("caso");
    all.insert(0, genotype);
    assert_both(
        "caso antes de la cabecera",
        &bytes_of(&all),
        &TruthValidationError::HeaderMustBeFirst,
    );
}

/// Una segunda cabecera en la línea 2 no es una cabecera duplicada sino una
/// cabecera que no está la primera.
#[test]
fn a_header_on_line_two_is_rejected_as_not_first() {
    let mut all = lines();
    all.insert(1, all[0].clone());
    assert_both(
        "cabecera en la línea 2",
        &bytes_of(&all),
        &TruthValidationError::HeaderMustBeFirst,
    );
}

/// Sin ninguna línea no hay documento.
#[test]
fn an_empty_input_is_an_empty_document() {
    assert_both("vacío", b"", &TruthValidationError::EmptyDocument);
}

/// El lector por flujo exige el salto de línea final como el completo.
#[test]
fn a_missing_final_newline_is_rejected_by_both_readers() {
    let mut text = bytes_of(&lines());
    assert_eq!(text.pop(), Some(b'\n'));
    assert_both(
        "sin salto final",
        &text,
        &TruthValidationError::MissingFinalNewline,
    );
}

/// Un lector de bytes que falla tras la primera línea.
struct FailsAfterTheFirstLine {
    inner: Cursor<Vec<u8>>,
    limit: u64,
}

impl Read for FailsAfterTheFirstLine {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        self.inner.read(buffer)
    }
}

impl BufRead for FailsAfterTheFirstLine {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        if self.inner.position() >= self.limit {
            return Err(std::io::Error::other("disk gone"));
        }
        self.inner.fill_buf()
    }
    fn consume(&mut self, amount: usize) {
        self.inner.consume(amount);
    }
}

/// Un fallo de lectura se devuelve con el número de la línea que se leía, la 2.
#[test]
fn a_read_failure_is_reported_with_the_line_being_read() {
    let text = bytes_of(&lines());
    let limit = u64::try_from(
        text.iter()
            .position(|byte| *byte == b'\n')
            .expect("línea 1"),
    )
    .expect("cabe")
        + 1;
    let expected = TruthValidationError::Io {
        line: 2,
        detail: "disk gone".to_owned(),
    };
    let full = read_truth_jsonl(&mut FailsAfterTheFirstLine {
        inner: Cursor::new(text.clone()),
        limit,
    });
    assert_eq!(full.err(), Some(expected.clone()), "lector completo");
    let stream = validate_truth_jsonl_streaming(
        &mut FailsAfterTheFirstLine {
            inner: Cursor::new(text),
            limit,
        },
        &reference(),
        |_, _| {},
    );
    assert_eq!(stream.err(), Some(expected), "lector por flujo");
}

// ---------------------------------------------------------------------------
// Enteros que no caben
// ---------------------------------------------------------------------------

fn detail_of(error: Option<TruthValidationError>) -> String {
    match error {
        Some(TruthValidationError::Schema { detail, .. }) => detail,
        other => panic!("se esperaba un error de esquema, y se obtuvo {other:?}"),
    }
}

fn assert_schema_detail_prefix(label: &str, text: &str, prefix: &str) {
    let [full, stream] = both(text.as_bytes());
    for (reader, result) in [("completo", full), ("flujo", stream)] {
        let detail = detail_of(result.err());
        assert!(
            detail.starts_with(prefix),
            "{label} ({reader}): se esperaba «{prefix}» y se obtuvo «{detail}»"
        );
    }
}

/// La fase de una copia es de un byte: 256 se rechaza con su motivo.
#[test]
fn a_copy_index_that_does_not_fit_a_byte_is_rejected() {
    let mut all = lines();
    let at = copy_line(&all, "I1", 0);
    all[at]["identity"]["copy_index"] = json!(256);
    assert_schema_detail_prefix("copy_index 256", &jsonl(&all), "entero fuera de u8");
}

/// La dosis de un alelo es de un byte: 256 se rechaza con su motivo, y un
/// decimal que no es entero con el suyo.
#[test]
fn an_allele_dosage_that_is_not_a_byte_is_rejected_with_its_own_reason() {
    let mut all = lines();
    let at = case_line(&all, "I1", "snp", 49);
    all[at]["allele_dosages"] = json!([256, 1]);
    assert_schema_detail_prefix("dosis 256", &jsonl(&all), "dosis fuera de u8");
    let mut all = lines();
    all[at]["allele_dosages"] = json!([1.5, 1]);
    assert_schema_detail_prefix("dosis 1.5", &jsonl(&all), "dosis no es un entero exacto");
}

/// Con 255 la dosis cabe en el byte: el rechazo siguiente ya no es el del
/// formato sino el del recuento de portadores (lado aceptado del límite).
#[test]
fn a_dosage_of_two_hundred_fifty_five_is_a_byte_and_fails_later() {
    let mut all = lines();
    let at = case_line(&all, "I1", "snp", 49);
    all[at]["allele_dosages"] = json!([255, 1]);
    let [full, stream] = both(jsonl(&all).as_bytes());
    for result in [full, stream] {
        assert!(
            matches!(
                result.err(),
                Some(TruthValidationError::Contract(
                    TruthContractError::DosageMismatch { declared: 255, .. }
                ))
            ),
            "255 cabe en un byte y el rechazo es el de la dosis"
        );
    }
}

/// Con 255 cabe en un byte el índice de copia, y el rechazo ya es del contrato.
#[test]
fn a_copy_index_of_two_hundred_fifty_five_is_a_byte() {
    let mut all = lines();
    let at = copy_line(&all, "I1", 0);
    all[at]["identity"]["copy_index"] = json!(255);
    let [full, stream] = both(jsonl(&all).as_bytes());
    for result in [full, stream] {
        let detail = match result.err() {
            Some(TruthValidationError::Schema { detail, .. }) => detail,
            Some(_) => String::new(),
            None => panic!("debió rechazarse"),
        };
        assert!(
            !detail.starts_with("entero fuera de u8"),
            "255 cabe en un byte"
        );
    }
}

// ---------------------------------------------------------------------------
// Los bloques de un tracto (solo el lector completo: versión 5)
// ---------------------------------------------------------------------------

/// Los tres sitios donde la verdad publica bloques: el fragmento vigente de
/// una copia, su historial y el alelo de un caso. Devuelve la ruta hasta el
/// **primer bloque** de cada uno.
fn block_places(all: &[Value]) -> [(&'static str, Vec<String>); 3] {
    let copy = all
        .iter()
        .position(|line| {
            line["record_type"] == "haplotype_copy"
                && line["current_pieces"].as_array().is_some_and(|pieces| {
                    pieces.iter().any(|piece| piece["kind"] == "repeat_allele")
                })
        })
        .expect("una copia con un tracto");
    let position = |field: &str, kind: &str| {
        all[copy][field]
            .as_array()
            .expect("lista")
            .iter()
            .position(|item| item["kind"] == kind)
            .expect("el elemento existe")
            .to_string()
    };
    let case = all
        .iter()
        .position(|line| line["record_type"] == "genotype" && line["case_kind"] == "repeat")
        .expect("un caso de repetición");
    let allele = all[case]["alleles"]
        .as_array()
        .expect("alelos")
        .iter()
        .position(|allele| allele["kind"] == "repeat")
        .expect("un alelo de repetición")
        .to_string();
    [
        (
            "fragmento vigente",
            vec![
                copy.to_string(),
                "current_pieces".to_owned(),
                position("current_pieces", "repeat_allele"),
                "blocks".to_owned(),
                "0".to_owned(),
            ],
        ),
        (
            "historial",
            vec![
                copy.to_string(),
                "history".to_owned(),
                position("history", "repeat_allele"),
                "blocks".to_owned(),
                "0".to_owned(),
            ],
        ),
        (
            "alelo del caso",
            vec![
                case.to_string(),
                "alleles".to_owned(),
                allele,
                "blocks".to_owned(),
                "0".to_owned(),
            ],
        ),
    ]
}

fn node_mut<'a>(all: &'a mut [Value], path: &[String]) -> &'a mut Value {
    let mut node = &mut all[path[0].parse::<usize>().expect("línea")];
    for step in &path[1..] {
        node = match step.parse::<usize>() {
            Ok(index) => &mut node[index],
            Err(_) => &mut node[step.as_str()],
        };
    }
    node
}

/// Un número JSON con su grafía exacta (`serde_json` conserva el texto).
fn number(text: &str) -> Value {
    serde_json::from_str(text).expect("número JSON")
}

/// El número de unidades no es las bases entre la longitud del motivo, o la
/// longitud declarada del motivo no es la de su secuencia, o el bloque no
/// tiene bases, o el motivo está vacío, o las unidades no caben en un
/// decimal: cada una con su motivo, en los tres sitios donde se publica un
/// bloque (el mismo tipo lee los tres).
#[test]
fn a_malformed_block_is_rejected_with_its_own_reason_in_every_place() {
    type Edit = fn(&mut Value);
    let units = "el número de unidades no es las bases entre la longitud del motivo";
    let motif = "la longitud declarada del motivo no es la de su secuencia";
    let cases: [(&str, Edit, &str); 7] = [
        (
            "unidades de más",
            |block| block["units"] = json!(1000.5),
            units,
        ),
        (
            "unidades de menos",
            |block| {
                let derived = block["bases"].as_f64().expect("bases") / 2.0;
                block["units"] = json!(derived - 1.0);
            },
            units,
        ),
        (
            "longitud del motivo de más",
            |block| block["motif_length"] = json!(3),
            motif,
        ),
        (
            "motivo vacío con longitud cero",
            |block| {
                block["motif"] = json!("");
                block["motif_length"] = json!(0);
            },
            motif,
        ),
        (
            "bloque sin bases",
            |block| {
                block["bases"] = json!(0);
                block["units"] = number("0.0");
            },
            "un bloque sin bases no es un bloque",
        ),
        (
            "unidades que no caben en un decimal",
            |block| block["units"] = number("1e400"),
            "el número de unidades no es un decimal",
        ),
        (
            "unidades apenas fuera de la tolerancia (2e-9 relativo)",
            |block| {
                let derived = block["bases"].as_f64().expect("bases") / 2.0;
                block["units"] = json!(derived + derived * 2.0e-9);
            },
            units,
        ),
    ];
    for place in 0..3 {
        for (label, edit, reason) in cases {
            let mut all = lines_v5();
            let path = block_places(&all)[place].1.clone();
            let before = all[path[0].parse::<usize>().expect("línea")].clone();
            edit(node_mut(&mut all, &path));
            assert_ne!(
                before,
                all[path[0].parse::<usize>().expect("línea")],
                "la alteración alteró"
            );
            let detail = detail_of(full(jsonl(&all).as_bytes()).err());
            assert!(
                detail.starts_with(reason),
                "{} / {label}: se esperaba «{reason}» y se obtuvo «{detail}»",
                block_places(&all)[place].0
            );
        }
    }
}

/// Lado aceptado de los bloques: unas unidades dentro de la tolerancia, y un
/// número de unidades decimal con unidad parcial, se leen.
#[test]
fn a_block_within_tolerance_or_with_a_partial_unit_is_read() {
    for place in 0..3 {
        let mut all = lines_v5();
        let path = block_places(&all)[place].1.clone();
        let block = node_mut(&mut all, &path);
        let derived = block["bases"].as_f64().expect("bases") / 2.0;
        block["units"] = json!(derived + derived * 9.0e-10);
        read_truth_jsonl(&mut Cursor::new(jsonl(&all))).expect("dentro de la tolerancia");
    }
    // Una unidad parcial: 81 bases de un motivo de dos, 40,5 unidades, en el
    // fragmento vigente (el lector del contrato no cruza los tres sitios).
    let mut all = lines_v5();
    let path = block_places(&all)[0].1.clone();
    let line: usize = path[0].parse().expect("línea");
    let piece: usize = path[2].parse().expect("fragmento");
    let bases = all[line]["current_pieces"][piece]["bases"]
        .as_str()
        .expect("bases")
        .to_owned();
    all[line]["current_pieces"][piece]["bases"] = json!(format!("{bases}C"));
    let block = node_mut(&mut all, &path);
    block["bases"] = json!(block["bases"].as_u64().expect("bases") + 1);
    let units = (block["bases"].as_f64().expect("bases") / 2.0).to_string();
    block["units"] = number(&units);
    read_truth_jsonl(&mut Cursor::new(jsonl(&all))).expect("unidad parcial");
}

// ---------------------------------------------------------------------------
// Referencia, orden de los registros y recuentos de la cabecera
// ---------------------------------------------------------------------------

/// Una verdad escrita por completo con otro identificador de referencia es
/// una verdad de **otra** referencia: los dos lectores la rechazan con los dos
/// identificadores. El lector completo la rechaza al validar la semántica
/// (el contrato la acepta: es coherente consigo misma) y el de flujo en la
/// cabecera.
#[test]
fn a_truth_of_another_reference_is_rejected_by_both_readers() {
    let text = jsonl(&lines()).replace("ref-main", "ref-other");
    assert_both(
        "otra referencia",
        text.as_bytes(),
        &TruthValidationError::ReferenceIdMismatch {
            declared: "ref-other".to_owned(),
            actual: "ref-main".to_owned(),
        },
    );
}

fn contract(error: TruthContractError) -> TruthValidationError {
    TruthValidationError::Contract(error)
}

/// Una copia detrás de un caso genotípico: el orden canónico trae todas las
/// copias antes que los casos.
#[test]
fn a_copy_after_a_genotype_is_out_of_order() {
    // La última copia del documento, puesta tras el primer caso: su identidad
    // es mayor que la de todas las copias anteriores, de modo que solo el
    // orden entre copias y casos puede rechazarla.
    let mut all = lines();
    let last_copy = all
        .iter()
        .rposition(|line| line["record_type"] == "haplotype_copy")
        .expect("una copia");
    let copy = all.remove(last_copy);
    let first_genotype = all
        .iter()
        .position(|line| line["record_type"] == "genotype")
        .expect("un caso");
    all.insert(first_genotype + 1, copy);
    assert_both(
        "copia tras un caso",
        &bytes_of(&all),
        &contract(TruthContractError::RecordsOutOfOrder),
    );
}

/// Dos copias en orden decreciente.
#[test]
fn two_copies_in_decreasing_order_are_out_of_order() {
    let mut all = lines();
    let first = copy_line(&all, "I1", 0);
    all.swap(first, first + 1);
    assert_both(
        "copias permutadas",
        &bytes_of(&all),
        &contract(TruthContractError::RecordsOutOfOrder),
    );
}

/// Una copia repetida: el lector por flujo la toma por desordenada (no admite
/// dos con la misma identidad) y el completo, que acepta el orden no
/// decreciente, la rechaza como copia duplicada.
#[test]
fn a_repeated_copy_is_out_of_order_for_the_stream_and_duplicate_for_the_whole() {
    let mut all = lines();
    let at = copy_line(&all, "I1", 0);
    all.insert(at + 1, all[at].clone());
    let text = bytes_of(&all);
    let [whole, stream] = both(&text);
    assert_eq!(
        whole.err(),
        Some(contract(TruthContractError::DuplicateCopy {
            identity: identity("I1", 0)
        }))
    );
    assert_eq!(
        stream.err(),
        Some(contract(TruthContractError::RecordsOutOfOrder))
    );
}

/// Dos casos en orden decreciente.
#[test]
fn two_genotypes_in_decreasing_order_are_out_of_order() {
    let mut all = lines();
    let first = case_line(&all, "I1", "snp", 49);
    all.swap(first, first + 1);
    assert_both(
        "casos permutados",
        &bytes_of(&all),
        &contract(TruthContractError::RecordsOutOfOrder),
    );
}

/// Un caso repetido: el lector por flujo lo toma por desordenado y el
/// completo, por caso duplicado.
#[test]
fn a_repeated_genotype_is_out_of_order_for_the_stream_and_duplicate_for_the_whole() {
    let mut all = lines();
    let at = case_line(&all, "I1", "snp", 49);
    let case_id = all[at]["case_id"]
        .as_str()
        .expect("identificador")
        .to_owned();
    all.insert(at + 1, all[at].clone());
    let text = bytes_of(&all);
    let [whole, stream] = both(&text);
    assert_eq!(
        whole.err(),
        Some(contract(TruthContractError::DuplicateGenotype { case_id }))
    );
    assert_eq!(
        stream.err(),
        Some(contract(TruthContractError::RecordsOutOfOrder))
    );
}

/// El recuento de copias de la cabecera no es el de las copias, en los dos
/// lectores (el de flujo lo juzga al final).
#[test]
fn the_declared_copy_count_must_be_the_number_of_copies() {
    let mut all = lines();
    let actual = all[0]["counts"]["copies"].as_u64().expect("recuento");
    all[0]["counts"]["copies"] = json!(actual - 1);
    assert_both(
        "copias declaradas",
        &bytes_of(&all),
        &contract(TruthContractError::CopyCountMismatch {
            declared: actual - 1,
            actual,
        }),
    );
}

/// El recuento de casos de la cabecera no es el de los casos.
#[test]
fn the_declared_genotype_count_must_be_the_number_of_genotypes() {
    let mut all = lines();
    let actual = all[0]["counts"]["genotypes"].as_u64().expect("recuento");
    all[0]["counts"]["genotypes"] = json!(actual - 1);
    assert_both(
        "casos declarados",
        &bytes_of(&all),
        &contract(TruthContractError::GenotypeCountMismatch {
            declared: actual - 1,
            actual,
        }),
    );
}

// ---------------------------------------------------------------------------
// La réplica del historial
// ---------------------------------------------------------------------------

/// Un cambio de base del historial que cae en un tramo que una deleción
/// anterior ya quitó no tiene fragmento donde aplicarse.
#[test]
fn a_base_change_in_an_already_deleted_stretch_does_not_replay() {
    let mut all = lines();
    let at = copy_line(&all, "I1", 1);
    let base = reference_char(421);
    all[at]["segments"][0]["modifications"][0]["position"] = json!(421);
    all[at]["segments"][0]["modifications"][0]["reference"] = json!(base.to_string());
    all[at]["segments"][0]["modifications"][0]["alternate"] = json!(other_base(base).to_string());
    let mut history = all[at]["history"].as_array().expect("historial").clone();
    history.swap(0, 1);
    assert_eq!(history[1]["kind"], "snp");
    history[1]["position"] = json!(421);
    all[at]["history"] = Value::Array(history);
    assert_both(
        "cambio en un tramo borrado",
        &bytes_of(&all),
        &replay_failure("I1", 1),
    );
}

/// Dos fragmentos que cubren la posición de un cambio de base: el historial
/// no dice a cuál se aplica. El documento está armado como si se aplicara al
/// último, de modo que sin la guarda lo pararía otra comprobación distinta.
#[test]
fn a_base_change_covered_by_two_fragments_does_not_replay() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 0);
    let mut segments = all[at]["segments"].as_array().expect("segmentos").clone();
    segments.push(segment(0, 100, vec![]));
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    assert_eq!(snps.len(), 2, "el cambio de 49 y el de 119");
    let history = all[at]["history"].as_array().expect("historial").clone();
    set_copy(
        &mut all,
        "I3",
        0,
        segments,
        vec![
            source_piece(0, 600, vec![snps[1].clone()]),
            source_piece(0, 100, vec![snps[0].clone()]),
        ],
        history,
    );
    assert_both("dos fragmentos", &bytes_of(&all), &replay_failure("I3", 0));
}

/// Una inserción en el interior de unas bases insertadas que abren la copia:
/// no hay ninguna base de referencia delante donde anclar su locus. El
/// documento está armado con el fragmento vigente que resultaría de ignorar
/// la guarda, para que solo ella pueda pararlo.
#[test]
fn an_insertion_with_no_reference_base_before_it_does_not_replay() {
    let mut all = lines();
    let provenance = explicit("ins-1");
    set_copy(
        &mut all,
        "I2",
        1,
        vec![segment(
            10,
            600,
            vec![
                json!({"kind": "insertion", "event_id": "ins-seg", "point": 10,
                        "bases": "ACG", "provenance": provenance}),
            ],
        )],
        vec![
            inserted_piece("A", &provenance),
            inserted_piece("T", &provenance),
            inserted_piece("CG", &provenance),
            source_piece(10, 600, vec![]),
        ],
        vec![
            json!({"kind": "insertion", "event_id": "ins-hist", "point": 1,
                    "bases": "T", "provenance": provenance}),
        ],
    );
    assert_both(
        "inserción sin ancla",
        &bytes_of(&all),
        &replay_failure("I2", 1),
    );
}

/// Una deleción que solo quita bases insertadas y no tiene ninguna base de
/// referencia delante.
#[test]
fn a_deletion_of_inserted_bases_with_no_reference_base_before_it_does_not_replay() {
    let mut all = lines();
    let provenance = explicit("ins-1");
    set_copy(
        &mut all,
        "I2",
        1,
        vec![segment(
            10,
            600,
            vec![
                json!({"kind": "insertion", "event_id": "ins-seg", "point": 10,
                        "bases": "ACG", "provenance": provenance}),
            ],
        )],
        vec![
            inserted_piece("A", &provenance),
            source_piece(10, 600, vec![]),
        ],
        vec![
            json!({"kind": "deletion", "event_id": "del-hist", "start": 1, "end": 3,
                    "provenance": explicit("del-1")}),
        ],
    );
    assert_both(
        "deleción sin ancla",
        &bytes_of(&all),
        &replay_failure("I2", 1),
    );
}

/// Un cambio de base heredado y una deleción de segmento que lo cubre, en los
/// dos órdenes: **el contrato lo rechaza** (DT-047, cerrada), con los dos
/// lectores y su motivo exacto. La réplica, que ya no lo recibe de un
/// `truth.jsonl`, conserva su propia guarda para el documento en memoria que no
/// pasó el contrato: el cursor queda por delante del inicio de la deleción y se
/// rechaza; el documento está armado con los fragmentos que resultarían de
/// ignorar la guarda.
#[test]
fn an_inherited_base_change_covered_by_a_segment_deletion_is_rejected() {
    let mut all = lines();
    let provenance = explicit("snp-direct");
    let snp = snp_modification("snp-cov", 12, &provenance);
    let deletion = json!({"kind": "deletion", "event_id": "del-cov", "start": 10, "end": 20,
                          "provenance": explicit("del-1")});
    set_copy(
        &mut all,
        "I2",
        0,
        vec![segment(0, 600, vec![snp, deletion])],
        vec![
            source_piece(0, 12, vec![]),
            source_piece(12, 13, vec![snp_in_piece("snp-cov", 12, &provenance)]),
            source_piece(20, 600, vec![]),
        ],
        vec![],
    );
    assert_both(
        "cambio y después la deleción que lo cubre",
        &bytes_of(&all),
        &contract(TruthContractError::ModificationsNotOrdered {
            event_id: "del-cov".to_owned(),
        }),
    );
    let at = copy_line(&all, "I2", 0);
    all[at]["segments"][0]["modifications"]
        .as_array_mut()
        .expect("modificaciones")
        .swap(0, 1);
    assert_both(
        "deleción y después el cambio que cubre",
        &bytes_of(&all),
        &contract(TruthContractError::OverlappingModifications {
            event_id: "snp-cov".to_owned(),
        }),
    );
}

/// La guarda de la réplica ante un documento en memoria que no pasó el
/// contrato: el cambio heredado y la deleción que lo cubre.
#[test]
fn the_replay_rejects_a_covered_change_in_a_document_that_skipped_the_contract() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I2", 0);
    let snp = typed_snp("snp-cov", 12);
    copy.segments = vec![typed_segment(
        0,
        600,
        vec![
            typed_snp_modification("snp-cov", 12),
            TruthModification::Deletion {
                event_id: "del-cov".to_owned(),
                start: 10,
                end: 20,
                provenance: typed_provenance("del-1"),
            },
        ],
    )];
    copy.current_pieces = vec![
        typed_piece(0, 12, vec![]),
        typed_piece(12, 13, vec![snp]),
        typed_piece(20, 600, vec![]),
    ];
    copy.history = vec![];
    assert_eq!(judge(&all).err(), Some(replay_failure("I2", 0)));
}

/// El mismo identificador de evento heredado dos veces en los segmentos de
/// una copia (dos cambios de base de posiciones distintas). El documento está
/// armado con los fragmentos que saldrían de aceptarlo.
#[test]
fn a_repeated_inherited_event_id_does_not_replay() {
    let mut all = lines();
    let provenance = explicit("snp-direct");
    let first = snp_modification("EVT_REPEATED", 49, &provenance);
    let second = snp_modification("EVT_REPEATED", 60, &provenance);
    set_copy(
        &mut all,
        "I1",
        1,
        vec![segment(0, 600, vec![first, second])],
        vec![
            source_piece(0, 49, vec![]),
            source_piece(49, 50, vec![snp_in_piece("EVT_REPEATED", 49, &provenance)]),
            source_piece(50, 60, vec![]),
            source_piece(60, 61, vec![snp_in_piece("EVT_REPEATED", 60, &provenance)]),
            source_piece(61, 600, vec![]),
        ],
        vec![],
    );
    assert_both(
        "identificador repetido",
        &bytes_of(&all),
        &replay_failure("I1", 1),
    );
}

// ---------------------------------------------------------------------------
// La materialización de cada copia
// ---------------------------------------------------------------------------

/// Un cambio de base heredado (sin evento en el historial, por tanto no
/// contrastado en la réplica) con una base de referencia que no es la de la
/// referencia: lo rechaza la materialización, con su evento y su posición.
#[test]
fn an_inherited_base_change_with_the_wrong_reference_base_is_rejected() {
    let mut all = lines();
    let provenance = explicit("snp-direct");
    let real = reference_char(119);
    let wrong = other_base(real);
    let alternate = if wrong == 'C' { 'G' } else { 'C' };
    let modification = json!({
        "kind": "snp", "event_id": "EVT_WRONG", "position": 119,
        "reference": wrong.to_string(), "alternate": alternate.to_string(),
        "provenance": provenance,
    });
    let in_piece = json!({
        "event_id": "EVT_WRONG", "position": 119,
        "reference": wrong.to_string(), "alternate": alternate.to_string(),
        "provenance": provenance,
    });
    set_copy(
        &mut all,
        "I2",
        1,
        vec![segment(0, 600, vec![modification])],
        vec![
            source_piece(0, 119, vec![]),
            source_piece(119, 120, vec![in_piece]),
            source_piece(120, 600, vec![]),
        ],
        vec![],
    );
    assert_both(
        "referencia equivocada",
        &bytes_of(&all),
        &TruthValidationError::SnpReferenceMismatch {
            event_id: "EVT_WRONG".to_owned(),
            position: 119,
        },
    );
}

/// Dos fragmentos de fuente que se solapan: el segundo reclama una posición
/// que el primero ya ocupa. La posición es la mayor de los dos inicios.
#[test]
fn two_overlapping_source_fragments_are_rejected() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 0);
    let mut segments = all[at]["segments"].as_array().expect("segmentos").clone();
    segments.push(segment(300, 350, vec![]));
    let mut pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    pieces.push(source_piece(300, 350, vec![]));
    let history = all[at]["history"].as_array().expect("historial").clone();
    set_copy(&mut all, "I3", 0, segments, pieces, history);
    assert_both(
        "fragmentos solapados",
        &bytes_of(&all),
        &TruthValidationError::DuplicateCurrentSourcePosition {
            identity: identity("I3", 0),
            position: 300,
        },
    );
}

/// El mismo solapamiento si el que reclama la posición es un tracto de
/// repetición: el fragmento que lo ocupa antes está en la copia delante de él.
#[test]
fn a_repeat_tract_over_an_occupied_stretch_is_rejected() {
    let mut all = lines_v5();
    let at = copy_line(&all, "I3", 0);
    let mut segments = all[at]["segments"].as_array().expect("segmentos").clone();
    segments.insert(0, segment(250, 260, vec![]));
    let mut pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    pieces.insert(0, source_piece(250, 260, vec![]));
    let history = all[at]["history"].as_array().expect("historial").clone();
    set_copy(&mut all, "I3", 0, segments, pieces, history);
    let error = full(jsonl(&all).as_bytes()).err();
    assert_eq!(
        error,
        Some(TruthValidationError::DuplicateCurrentSourcePosition {
            identity: identity("I3", 0),
            position: 250,
        })
    );
}

// ---------------------------------------------------------------------------
// La fuente de cada fragmento
// ---------------------------------------------------------------------------

/// Una fuente del fragmento que es de la clase `founder` (que el contrato
/// admite) o que es de la clase `reference` pero con el identificador de otra
/// referencia: no nombra la referencia que se juzga.
#[test]
fn a_source_that_does_not_name_the_reference_is_rejected() {
    for (label, kind, source_id) in [
        (
            "fundador con el identificador de la referencia",
            "founder",
            "ref-main",
        ),
        (
            "referencia con otro identificador",
            "reference",
            "ref-other",
        ),
    ] {
        let mut all = lines();
        all[0]["source_digests"]
            .as_array_mut()
            .expect("huellas")
            .push(
                json!({"source": {"kind": kind, "source_id": source_id}, "sha256": "d".repeat(64)}),
            );
        let at = copy_line(&all, "I3", 1);
        let alien = json!({"kind": kind, "source_id": source_id});
        all[at]["segments"][0]["source"]["source"] = alien.clone();
        for piece in all[at]["current_pieces"]
            .as_array_mut()
            .expect("fragmentos")
        {
            piece["source"]["source"] = alien.clone();
        }
        assert_both(
            label,
            &bytes_of(&all),
            &TruthValidationError::SourceDoesNotNameReference {
                identity: identity("I3", 1),
            },
        );
    }
}

/// Una secuencia que la referencia no tiene.
#[test]
fn a_sequence_the_reference_lacks_is_rejected() {
    // «chr0» ordena antes que «chr2» como «chr1», de modo que el documento
    // sigue en orden canónico.
    let text = jsonl(&lines()).replace("chr1", "chr0");
    assert_both(
        "secuencia ajena",
        text.as_bytes(),
        &TruthValidationError::MissingReferenceSequence {
            sequence_id: "chr0".to_owned(),
        },
    );
}

/// Una longitud de secuencia declarada que no es la de la referencia.
#[test]
fn a_declared_sequence_length_that_is_not_the_references_is_rejected() {
    let text = jsonl(&lines()).replace("\"sequence_length\":600", "\"sequence_length\":601");
    assert_both(
        "longitud declarada",
        text.as_bytes(),
        &TruthValidationError::SourceLengthMismatch {
            sequence_id: "chr1".to_owned(),
            declared: 601,
            actual: 600,
        },
    );
}

// ---------------------------------------------------------------------------
// El contraste de cada caso con las copias
// ---------------------------------------------------------------------------

fn genotype_mismatch(all: &[Value], at: usize, detail: &'static str) -> TruthValidationError {
    TruthValidationError::GenotypeMismatch {
        case_id: all[at]["case_id"]
            .as_str()
            .expect("identificador")
            .to_owned(),
        detail,
    }
}

fn without_carriers(line: &mut Value) {
    line["carriers"] = json!([]);
    let count = line["alleles"].as_array().expect("alelos").len();
    line["allele_dosages"] = Value::Array(vec![json!(0); count]);
}

/// El último caso del documento que es de `chr1`: el siguiente es de `chr2`,
/// de modo que mover su locus dentro de `chr1` no rompe el orden.
fn last_case_on_chr1(all: &[Value]) -> usize {
    all.iter()
        .rposition(|line| {
            line["record_type"] == "genotype" && line["locus"]["sequence_id"] == "chr1"
        })
        .expect("un caso de chr1")
}

/// Un caso SNP en una secuencia que la referencia no tiene.
#[test]
fn a_snp_case_on_a_sequence_the_reference_lacks_is_rejected() {
    // Un caso SNP del último individuo, puesto el último y con la secuencia
    // «chrX», que ordena detrás de todas: el documento sigue en orden canónico.
    let mut all = lines();
    let line = all
        .iter()
        .rposition(|line| line["individual_id"] == "I3" && line["case_kind"] == "snp")
        .expect("un caso SNP de I3");
    let moved = all.remove(line);
    all.push(moved);
    let at = all.len() - 1;
    all[at]["locus"]["sequence_id"] = json!("chrX");
    without_carriers(&mut all[at]);
    let expected = genotype_mismatch(&all, at, "la secuencia del locus no existe");
    assert_both("SNP en secuencia ajena", &bytes_of(&all), &expected);
}

/// El locus de un caso SNP en la primera posición fuera de la referencia (600)
/// se rechaza; en la última dentro (599) el rechazo es otro, el de los
/// portadores: es el lado aceptado del límite.
#[test]
fn a_snp_locus_is_inside_the_reference_up_to_its_last_base() {
    for (start, detail) in [
        (600_u64, "el locus SNP sale de la referencia"),
        (
            599,
            "las copias portadoras no coinciden con la reconstrucción",
        ),
    ] {
        let mut all = lines();
        let at = last_case_on_chr1(&all);
        all[at]["locus"]["start"] = json!(start);
        all[at]["locus"]["end"] = json!(start + 1);
        without_carriers(&mut all[at]);
        let expected = genotype_mismatch(&all, at, detail);
        assert_both(&format!("locus {start}"), &bytes_of(&all), &expected);
    }
}

/// Un caso cuyas copias portadoras están cambiadas entre sí (con las mismas
/// dosis, de modo que el contrato no lo ve): lo rechaza la reconstrucción.
#[test]
fn carriers_swapped_between_copies_are_rejected_in_every_kind_of_case() {
    let reason = "las copias portadoras no coinciden con la reconstrucción";
    for (individual, kind, start) in [
        ("I1", "snp", 49_u64),
        ("I1", "indel", 350),
        ("I1", "indel", 420),
    ] {
        let mut all = lines();
        let at = case_line(&all, individual, kind, start);
        let carriers = all[at]["carriers"].as_array_mut().expect("portadores");
        let first = carriers[0]["allele_index"].clone();
        carriers[0]["allele_index"] = carriers[1]["allele_index"].clone();
        carriers[1]["allele_index"] = first;
        let expected = genotype_mismatch(&all, at, reason);
        assert_both(&format!("{kind} {start}"), &bytes_of(&all), &expected);
    }
}

/// Lo mismo en un caso de repetición (versión 5: solo el lector completo).
#[test]
fn carriers_swapped_between_copies_are_rejected_in_a_repeat_case() {
    let mut all = lines_v5();
    let at = all
        .iter()
        .position(|line| {
            line["case_kind"] == "repeat"
                && line["carriers"][0]["allele_index"] != line["carriers"][1]["allele_index"]
        })
        .expect("un caso con dos alelos distintos");
    let carriers = all[at]["carriers"].as_array_mut().expect("portadores");
    let first = carriers[0]["allele_index"].clone();
    carriers[0]["allele_index"] = carriers[1]["allele_index"].clone();
    carriers[1]["allele_index"] = first;
    let expected = genotype_mismatch(
        &all,
        at,
        "las copias portadoras no coinciden con la reconstrucción",
    );
    assert_eq!(full(jsonl(&all).as_bytes()).err(), Some(expected));
}

/// Un caso indel con dos alelos alternativos, o cuyo primer alelo no es el de
/// referencia: el juez solo admite REF y un único alelo de indel.
#[test]
fn an_indel_case_must_be_reference_and_one_alternate() {
    let reason = "F2 requiere REF y un único alelo indel alternativo";
    // Dos alternativos: a la inserción se le añade la deleción de otro caso.
    let mut all = lines();
    let at = case_line(&all, "I1", "indel", 350);
    let deletion = all[case_line(&all, "I1", "indel", 420)]["alleles"][1].clone();
    all[at]["alleles"]
        .as_array_mut()
        .expect("alelos")
        .push(deletion);
    all[at]["allele_dosages"] = json!([1, 1, 0]);
    let expected = genotype_mismatch(&all, at, reason);
    assert_both("dos alternativos", &bytes_of(&all), &expected);
    // El primer alelo no es el de referencia: se permutan los dos y sus dosis.
    let mut all = lines();
    all[at]["alleles"]
        .as_array_mut()
        .expect("alelos")
        .swap(0, 1);
    for carrier in all[at]["carriers"].as_array_mut().expect("portadores") {
        let index = carrier["allele_index"].as_u64().expect("índice");
        carrier["allele_index"] = json!(1 - index);
    }
    let expected = genotype_mismatch(&all, at, reason);
    assert_both("primero no es REF", &bytes_of(&all), &expected);
}

// ---------------------------------------------------------------------------
// Casos de repetición (versión 5: solo el lector completo)
// ---------------------------------------------------------------------------

fn repeat_case(all: &[Value], individual: &str) -> usize {
    case_line(all, individual, "repeat", 200)
}

/// Un caso de repetición con dos alelos de repetición de loci distintos.
#[test]
fn a_repeat_case_with_alleles_of_two_loci_is_rejected() {
    let mut all = lines_v5();
    let at = repeat_case(&all, "I2");
    let alleles = all[at]["alleles"].as_array_mut().expect("alelos");
    assert_eq!(alleles.len(), 3, "referencia y dos alternativos");
    alleles[2]["locus_id"] = json!("STR-2");
    let expected = genotype_mismatch(
        &all,
        at,
        "todos los alelos de un caso de repetición pertenecen a su locus",
    );
    assert_eq!(full(jsonl(&all).as_bytes()).err(), Some(expected));
}

/// Un caso de repetición en una secuencia que la referencia no tiene.
#[test]
fn a_repeat_case_on_a_sequence_the_reference_lacks_is_rejected() {
    // Tras el último caso del documento, que es el orden de «chrX».
    let mut all = lines_v5();
    let line = all.remove(repeat_case(&all, "I3"));
    all.push(line);
    let at = all.len() - 1;
    all[at]["locus"]["sequence_id"] = json!("chrX");
    without_carriers(&mut all[at]);
    let expected = genotype_mismatch(&all, at, "la secuencia del locus no existe");
    assert_eq!(full(jsonl(&all).as_bytes()).err(), Some(expected));
}

/// El locus de un caso de repetición tiene que ser un intervalo no vacío
/// dentro de la referencia: con el final en 601 o con inicio y final iguales
/// se rechaza; con el final en 600 (la última base) ya no, y el rechazo es el
/// de los portadores.
#[test]
fn a_repeat_locus_is_a_non_empty_interval_inside_the_reference() {
    let outside = "el locus de repetición sale de la referencia";
    for (label, start, end, detail) in [
        ("más allá del final", 200_u64, 601_u64, outside),
        ("vacío", 250, 250, outside),
        (
            "hasta la última base",
            200,
            600,
            "las copias portadoras no coinciden con la reconstrucción",
        ),
        (
            "con el inicio desplazado",
            201,
            300,
            "las copias portadoras no coinciden con la reconstrucción",
        ),
    ] {
        let mut all = lines_v5();
        let at = repeat_case(&all, "I3");
        all[at]["locus"]["start"] = json!(start);
        all[at]["locus"]["end"] = json!(end);
        let expected = genotype_mismatch(&all, at, detail);
        assert_eq!(
            full(jsonl(&all).as_bytes()).err(),
            Some(expected),
            "{label}"
        );
    }
}

/// Un alelo alternativo que repite el tracto de la referencia no es un alelo.
#[test]
fn a_repeat_allele_equal_to_the_reference_tract_is_rejected() {
    let mut all = lines_v5();
    let at = repeat_case(&all, "I3");
    let tract = std::str::from_utf8(&reference_bases()[200..300])
        .expect("ASCII")
        .to_owned();
    all[at]["alleles"][1]["bases"] = json!(tract);
    all[at]["alleles"][1]["blocks"] = json!([
        {"motif": "CA", "motif_length": 2, "bases": 100, "units": number("50.0")}
    ]);
    let expected = genotype_mismatch(
        &all,
        at,
        "un alelo alternativo repite el tracto de la referencia",
    );
    assert_eq!(full(jsonl(&all).as_bytes()).err(), Some(expected));
}

/// Una copia lleva un tracto (con su secuencia y sus bloques, coherentes en
/// su fragmento vigente y en su historial) que el caso no cataloga.
#[test]
fn a_copy_carrying_a_tract_the_case_does_not_catalogue_is_rejected() {
    let mut all = lines_v5();
    let at = copy_line(&all, "I3", 0);
    let interrupted = format!("CATT{}", "CA".repeat(43));
    assert_eq!(interrupted.len(), 90);
    let blocks = json!([
        {"motif": "CATT", "motif_length": 4, "bases": 4, "units": number("1.0")},
        {"motif": "CA", "motif_length": 2, "bases": 86, "units": number("43.0")},
    ]);
    for field in ["current_pieces", "history"] {
        for item in all[at][field].as_array_mut().expect("lista") {
            if item["kind"] == "repeat_allele" {
                item["bases"] = json!(interrupted);
                item["blocks"] = blocks.clone();
            }
        }
    }
    let case = repeat_case(&all, "I3");
    let expected = genotype_mismatch(
        &all,
        case,
        "la copia lleva un alelo de repetición que el caso no cataloga",
    );
    assert_eq!(full(jsonl(&all).as_bytes()).err(), Some(expected));
}

/// Una copia cuyo tracto es de **otro locus** que el del caso no cuenta como
/// portadora de este: el caso dice que la lleva y la reconstrucción que no.
#[test]
fn a_tract_of_another_locus_is_not_the_cases_allele() {
    let mut all = lines_v5();
    let at = copy_line(&all, "I3", 0);
    for field in ["current_pieces", "history"] {
        for item in all[at][field].as_array_mut().expect("lista") {
            if item["kind"] == "repeat_allele" {
                item["locus_id"] = json!("STR-2");
            }
        }
    }
    let case = repeat_case(&all, "I3");
    let expected = genotype_mismatch(
        &all,
        case,
        "las copias portadoras no coinciden con la reconstrucción",
    );
    assert_eq!(full(jsonl(&all).as_bytes()).err(), Some(expected));
}

/// Una inserción cuyo punto cae dentro de un tracto ya realizado: el tracto no
/// se parte.
#[test]
fn an_insertion_inside_a_realized_tract_does_not_replay() {
    let mut all = lines_v5();
    let at = copy_line(&all, "I3", 0);
    let mut history = all[at]["history"].as_array().expect("historial").clone();
    history.push(
        json!({"kind": "insertion", "event_id": "ins-x", "point": 250,
                        "bases": "T", "provenance": explicit("ins-1")}),
    );
    all[at]["history"] = Value::Array(history);
    assert_eq!(
        full(jsonl(&all).as_bytes()).err(),
        Some(replay_failure("I3", 0))
    );
}

/// Un tracto cuyo intervalo no cabe entero en un fragmento de fuente sin
/// cambios de base: lo cruza el borde de una inserción anterior, o contiene un
/// cambio de base.
#[test]
fn a_repeat_event_that_does_not_fit_one_clean_fragment_does_not_replay() {
    for (label, start, end) in [
        ("cruza un borde", 340_u64, 360_u64),
        ("contiene un cambio de base", 390, 410),
    ] {
        let mut all = lines_v5();
        let at = copy_line(&all, "I1", 0);
        for item in all[at]["history"].as_array_mut().expect("historial") {
            if item["kind"] == "repeat_allele" {
                item["start"] = json!(start);
                item["end"] = json!(end);
            }
        }
        assert_eq!(
            full(jsonl(&all).as_bytes()).err(),
            Some(replay_failure("I1", 0)),
            "{label}"
        );
    }
}

// ---------------------------------------------------------------------------
// Deleciones en el borde de la copia y longitudes que desbordan
// ---------------------------------------------------------------------------

/// La deleción `del-1` de la copia 1 de I1 con otro final, y con todo lo que
/// la verdad publica de ella (historial, fragmentos vigentes y caso) coherente.
fn with_deletion_ending_at(end: u64) -> Vec<Value> {
    let mut all = lines();
    let at = copy_line(&all, "I1", 1);
    for event in all[at]["history"].as_array_mut().expect("historial") {
        if event["kind"] == "deletion" {
            event["end"] = json!(end);
        }
    }
    let mut pieces = vec![all[at]["current_pieces"][0].clone()];
    if end < 600 {
        pieces.push(source_piece(end, 600, vec![]));
    }
    all[at]["current_pieces"] = Value::Array(pieces);
    let case = case_line(&all, "I1", "indel", 420);
    all[case]["locus"]["end"] = json!(end);
    all[case]["alleles"][1]["end"] = json!(end);
    all[case]["alleles"][1]["source_event"]["end"] = json!(end);
    // La deleción se come también la base 499, que esta copia ya no tiene: su
    // caso SNP solo lo lleva la otra copia.
    let snp = case_line(&all, "I1", "snp", 499);
    let carriers = all[snp]["carriers"].as_array_mut().expect("portadores");
    carriers.retain(|carrier| carrier["copy"]["copy_index"] == 0);
    all[snp]["allele_dosages"] = json!([0, 1]);
    all
}

/// Una deleción que llega a la última base de la copia no se admite (la
/// réplica exige que quede al menos una base de referencia detrás), y una que
/// termina una base antes sí: lado aceptado del límite, con el documento
/// coherente en los dos casos para que solo la guarda pueda distinguirlos.
#[test]
fn a_deletion_may_not_reach_the_end_of_the_copy() {
    for result in both(&bytes_of(&with_deletion_ending_at(599))) {
        result.expect("una deleción que deja la última base se acepta");
    }
    assert_both(
        "deleción hasta el final",
        &bytes_of(&with_deletion_ending_at(600)),
        &replay_failure("I1", 1),
    );
}

fn huge_reference_object(start: u64, end: u64) -> Value {
    reference_object(start, end, u64::MAX)
}

/// Dos fragmentos de fuente cuyas longitudes suman más de lo que cabe en 64
/// bits: la réplica de una inserción no puede medir la copia. El documento
/// está armado con los fragmentos que saldrían de ignorar el desborde.
#[test]
fn lengths_that_overflow_do_not_replay() {
    let mut all = lines();
    let huge = |start: u64, end: u64| json!({"kind": "source", "source": huge_reference_object(start, end), "snps": []});
    let at = copy_line(&all, "I3", 0);
    all[at]["segments"] = json!([
        {"source": huge_reference_object(0, u64::MAX), "modifications": []},
        {"source": huge_reference_object(0, u64::MAX), "modifications": []},
    ]);
    all[at]["current_pieces"] = json!([
        huge(0, 5),
        inserted_piece("T", &explicit("ins-1")),
        huge(5, u64::MAX),
        huge(0, u64::MAX),
    ]);
    all[at]["history"] = json!([{
        "kind": "insertion", "event_id": "ins-big", "point": 5, "bases": "T",
        "provenance": explicit("ins-1"),
    }]);
    assert_both(
        "longitudes que desbordan",
        &bytes_of(&all),
        &replay_failure("I3", 0),
    );
}

// ---------------------------------------------------------------------------
// Documentos en memoria que el contrato no deja pasar
// ---------------------------------------------------------------------------
//
// Lo que sigue **no se alcanza con un `truth.jsonl` leído**: el contrato lo
// rechaza antes. Sí se alcanza con un documento en memoria que no ha pasado
// por el contrato, que es lo que recibe `validate_and_reconstruct_truth`, y
// por eso el juez lo comprueba por su cuenta.

/// Sin la huella de la referencia, o con la huella repetida: ni cero ni dos
/// son una.
#[test]
fn a_header_without_exactly_one_reference_digest_is_rejected() {
    let mut missing = document(&lines());
    missing
        .header
        .source_digests
        .retain(|digest| digest.source.kind != TruthSourceKind::Reference);
    let mut repeated = document(&lines());
    let reference_digest = repeated
        .header
        .source_digests
        .iter()
        .find(|digest| digest.source.kind == TruthSourceKind::Reference)
        .expect("huella de la referencia")
        .clone();
    repeated.header.source_digests.push(reference_digest);
    for (label, document) in [("sin huella", missing), ("huella repetida", repeated)] {
        assert_eq!(
            judge(&document).err(),
            Some(TruthValidationError::MissingOrDuplicateReferenceDigest),
            "{label}"
        );
    }
}

/// Un cambio de base heredado sobre un fragmento que ya lleva un cambio de
/// base en esa posición. Los dos segmentos se solapan y el segundo quita la
/// posición con una deleción de segmento, de modo que solo el primero la
/// cubre. El documento está armado con el fragmento que resultaría de aceptar
/// los dos cambios, para que sin la guarda no lo pare el contraste final. El
/// contrato rechaza este documento (el cambio dentro de una deleción, DT-047),
/// así que solo se alcanza con un documento en memoria.
#[test]
fn two_base_changes_at_one_position_do_not_replay() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I3", 1);
    let history = TruthHistoryEvent::Snp {
        event_id: "EVT_SECOND".to_owned(),
        position: 60,
        provenance: typed_provenance("snp-direct"),
    };
    copy.segments = vec![
        typed_segment(0, 600, vec![typed_snp_modification("EVT_FIRST", 60)]),
        typed_segment(
            0,
            600,
            vec![
                typed_snp_modification("EVT_SECOND", 60),
                TruthModification::Deletion {
                    event_id: "DEL_COVER".to_owned(),
                    start: 55,
                    end: 65,
                    provenance: typed_provenance("del-1"),
                },
            ],
        ),
    ];
    copy.history = vec![history];
    copy.current_pieces = vec![
        typed_piece(0, 60, vec![]),
        typed_piece(
            60,
            61,
            vec![typed_snp("EVT_SECOND", 60), typed_snp("EVT_FIRST", 60)],
        ),
        typed_piece(61, 600, vec![]),
        typed_piece(0, 55, vec![]),
        typed_piece(65, 600, vec![]),
    ];
    assert_eq!(judge(&all).err(), Some(replay_failure("I3", 1)));
}

/// Una deleción del historial que empieza en la base cero. El documento está
/// armado con lo que quedaría de la copia si se aceptara.
#[test]
fn a_deletion_from_the_first_base_does_not_replay() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I1", 1);
    for event in &mut copy.history {
        if let TruthHistoryEvent::Deletion { start, end, .. } = event {
            *start = 0;
            *end = 424;
        }
    }
    copy.current_pieces = vec![typed_piece(424, 600, vec![])];
    assert_eq!(judge(&all).err(), Some(replay_failure("I1", 1)));
}

/// Una deleción del historial con el principio por detrás del final no quita
/// nada, y la réplica lo rechaza igual que las demás: el rechazo es el mismo
/// con o sin la comprobación de los extremos, y aquí solo se ejercita que el
/// documento se rechaza con el motivo de la réplica.
#[test]
fn a_deletion_with_its_start_after_its_end_does_not_replay() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I1", 1);
    for event in &mut copy.history {
        if let TruthHistoryEvent::Deletion { start, end, .. } = event {
            *start = 430;
            *end = 420;
        }
    }
    assert_eq!(judge(&all).err(), Some(replay_failure("I1", 1)));
}

/// Un cambio de base heredado cuya posición es la última representable: la
/// réplica no puede sumarle uno. El documento está armado con el único
/// fragmento que saldría de saturar la suma.
#[test]
fn an_inherited_base_change_at_the_last_position_does_not_replay() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I2", 1);
    let source = typed_source(0, u64::MAX, u64::MAX);
    copy.segments = vec![TruthSegment {
        source: source.clone(),
        modifications: vec![TruthModification::Snp {
            event_id: "EVT_LAST".to_owned(),
            position: u64::MAX,
            reference: "A".to_owned(),
            alternate: "C".to_owned(),
            provenance: typed_provenance("snp-direct"),
        }],
    }];
    copy.history = vec![];
    copy.current_pieces = vec![TruthCurrentPiece::Source {
        source,
        snps: vec![],
    }];
    assert_eq!(judge(&all).err(), Some(replay_failure("I2", 1)));
}

/// Una fuente que declara otra referencia que la que se juzga (el contrato
/// exige que coincida con la de la copia, y la copia con la cabecera).
#[test]
fn a_source_that_declares_another_reference_is_rejected() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I3", 1);
    copy.segments[0].source.reference_id = "ref-other".to_owned();
    for piece in &mut copy.current_pieces {
        if let TruthCurrentPiece::Source { source, .. } = piece {
            source.reference_id = "ref-other".to_owned();
        }
    }
    assert_eq!(
        judge(&all).err(),
        Some(TruthValidationError::SourceDoesNotNameReference {
            identity: identity("I3", 1),
        })
    );
}

/// Un fragmento que acaba una base más allá de la referencia, con la longitud
/// declarada que es la de la referencia (el contrato exige que el final no
/// pase de la longitud declarada).
#[test]
fn a_source_past_the_end_of_the_sequence_is_rejected() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I3", 1);
    copy.segments[0].source.end = 601;
    copy.current_pieces = vec![typed_piece(0, 601, vec![])];
    copy.history = vec![];
    copy.segments[0].modifications = vec![];
    assert_eq!(
        judge(&all).err(),
        Some(TruthValidationError::SourceIntervalOutOfBounds {
            identity: identity("I3", 1),
        })
    );
}

/// Un caso de repetición que no lleva ningún alelo de repetición.
#[test]
fn a_repeat_case_without_a_repeat_allele_is_rejected() {
    let mut all = document(&lines());
    let case = genotype_of(&mut all, "I1", GenotypeCaseKind::Snp, 49);
    case.case_kind = GenotypeCaseKind::Repeat;
    let case_id = case.case_id.clone();
    assert_eq!(
        judge(&all).err(),
        Some(TruthValidationError::GenotypeMismatch {
            case_id,
            detail: "un caso de repetición publica al menos un alelo alternativo",
        })
    );
}

/// Un caso indel sin ningún alelo alternativo: el juez lo rechaza en vez de
/// indexar el alelo que no hay.
#[test]
fn an_indel_case_with_no_alternate_allele_is_rejected() {
    let mut all = document(&lines());
    let case = genotype_of(&mut all, "I1", GenotypeCaseKind::Indel, 350);
    case.alleles = vec![TruthAllele::Reference];
    case.allele_dosages = vec![2];
    for carrier in &mut case.carriers {
        carrier.allele_index = 0;
    }
    let case_id = case.case_id.clone();
    assert_eq!(
        judge(&all).err(),
        Some(TruthValidationError::GenotypeMismatch {
            case_id,
            detail: "F2 requiere REF y un único alelo indel alternativo",
        })
    );
}

/// Un caso indel cuyo único alelo alternativo es de una base: no nombra su
/// evento estructural, y no hay ancla con que juzgar el conjunto de casos.
#[test]
fn an_indel_case_naming_no_structural_event_is_rejected() {
    let mut all = document(&lines());
    let case = genotype_of(&mut all, "I1", GenotypeCaseKind::Snp, 49);
    case.case_kind = GenotypeCaseKind::Indel;
    case.allele_dosages = vec![2, 0];
    for carrier in &mut case.carriers {
        carrier.allele_index = 0;
    }
    let case_id = case.case_id.clone();
    assert_eq!(
        judge(&all).err(),
        Some(TruthValidationError::GenotypeMismatch {
            case_id,
            detail: "un caso indel nombra su evento estructural",
        })
    );
}

/// Un tracto de repetición con el inicio detrás del final, o sin bases.
#[test]
fn a_repeat_event_without_an_interval_or_bases_does_not_replay() {
    for (label, start, end, empty) in [
        ("inicio igual al final", 250_u64, 250_u64, false),
        ("sin bases", 200, 300, true),
    ] {
        let mut all = document(&lines_v5());
        let copy = copy_of(&mut all, "I3", 0);
        let mut replaced = Vec::new();
        let mut event_provenance = typed_provenance("snp-direct");
        for event in &mut copy.history {
            if let TruthHistoryEvent::RepeatAllele {
                start: event_start,
                end: event_end,
                bases,
                blocks,
                provenance,
                ..
            } = event
            {
                event_provenance = provenance.clone();
                *event_start = start;
                *event_end = end;
                if empty {
                    bases.clear();
                    *blocks = vec![TruthRepeatBlock {
                        motif: "CA".to_owned(),
                        bases: 0,
                    }];
                }
                replaced.push((bases.clone(), blocks.clone()));
            }
        }
        // Los fragmentos vigentes que saldrían de aceptar el evento.
        let (bases, blocks) = replaced.pop().expect("un evento de repetición");
        let snps = match &copy.current_pieces[0] {
            TruthCurrentPiece::Source { snps, .. } => snps.clone(),
            _ => panic!("el primer fragmento es de fuente"),
        };
        copy.current_pieces = vec![
            typed_piece(0, start, snps),
            TruthCurrentPiece::RepeatAllele {
                locus_id: "STR-1".to_owned(),
                source: typed_source(start, end, 600),
                bases,
                blocks,
                provenance: event_provenance.clone(),
            },
            typed_piece(end, 600, vec![]),
        ];
        assert_eq!(judge(&all).err(), Some(replay_failure("I3", 0)), "{label}");
    }
}

// ---------------------------------------------------------------------------
// Una referencia sin bases, y las comprobaciones que ninguna entrada alcanza
// ---------------------------------------------------------------------------

/// La ficha de la referencia (con su huella calculada) sin las bases
/// materializadas: la que usa una ejecución que solo deriva referencias
/// reducidas.
fn reference_without_bases() -> Reference {
    use genome_population_simulator::reference_fasta::{
        SequenceExpectation, reference_from_summary, validate_supplied_fasta_summary,
    };
    let root = std::env::temp_dir().join(format!(
        "gps-truth-guards-summary-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("reloj")
            .as_nanos()
    ));
    std::fs::create_dir_all(&root).expect("directorio");
    let path = root.join("reference.fa");
    std::fs::write(
        &path,
        format!(
            ">chr1\n{}\n",
            std::str::from_utf8(&reference_bases()).expect("ASCII")
        ),
    )
    .expect("FASTA");
    let summary = validate_supplied_fasta_summary(
        &path,
        "ref-main",
        &[SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 600,
        }],
    )
    .expect("ficha de la referencia");
    let reference = reference_from_summary(&summary, "ref-main", ReferenceOrigin::Generated)
        .expect("referencia sin bases");
    std::fs::remove_dir_all(root).expect("limpieza");
    reference
}

/// Con una referencia sin bases, los tres lectores públicos rechazan la verdad
/// con `BasesNotMaterialized`, y lo hacen **al contrastar la huella de la
/// referencia** (que se calcula sobre las bases) antes de reconstruir
/// ninguna copia.
///
/// **No discrimina** los tres sitios del juez que leen bases de la referencia
/// (`bases_of`, la lectura de la base de un cambio durante la réplica y la de
/// la materialización): son inalcanzables precisamente por esto. Ninguna
/// entrada pública llega a ellos con una referencia sin bases, porque la huella
/// falla antes. Esta prueba ejecuta el intento y muestra que el rechazo es el
/// de la huella.
#[test]
fn a_reference_without_bases_is_rejected_before_any_copy_is_rebuilt() {
    let text = jsonl(&lines());
    let reference = reference_without_bases();
    assert_eq!(
        validate_and_reconstruct_truth_jsonl(&mut Cursor::new(text.clone()), &reference).err(),
        Some(TruthValidationError::BasesNotMaterialized)
    );
    assert_eq!(
        validate_truth_jsonl_streaming(&mut Cursor::new(text.clone()), &reference, |_, _| {
            panic!("ninguna copia llega a reconstruirse")
        })
        .err(),
        Some(TruthValidationError::BasesNotMaterialized)
    );
    assert_eq!(
        validate_and_reconstruct_truth(document(&lines()), &reference).err(),
        Some(TruthValidationError::BasesNotMaterialized)
    );
}

/// El intento de alcanzar las comprobaciones inalcanzables del juez, y lo que
/// gana la carrera en cada uno. Esta prueba **no discrimina** esas guardas;
/// solo muestra qué otra comprobación rechaza antes. Cada guarda lleva anotado
/// su porqué en `src/truth_validation.rs`.
#[test]
fn unreachable_guards_are_shadowed_by_other_checks() {
    // 1. Un identificador de evento heredado repetido, con su evento en el
    //    historial: el rechazo de `validate_structural_event_ids` llega antes
    //    que el de «dos modificaciones para el mismo evento» de la réplica de
    //    un cambio de base.
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    let mut repeated = all[at]["segments"][0]["modifications"][0].clone();
    repeated["position"] = json!(60);
    repeated["reference"] = json!(reference_char(60).to_string());
    repeated["alternate"] = json!(other_base(reference_char(60)).to_string());
    all[at]["segments"][0]["modifications"]
        .as_array_mut()
        .expect("modificaciones")
        .insert(1, repeated);
    let [whole, stream] = both(&bytes_of(&all));
    assert_eq!(whole.err(), Some(replay_failure("I3", 1)));
    assert_eq!(stream.err(), Some(replay_failure("I3", 1)));

    // 2. Un segmento con el inicio por detrás del final: la réplica lo rechaza
    //    al formar el fragmento inicial, de modo que ningún fragmento de fuente
    //    llega a tener longitud negativa y todas las sumas de longitudes
    //    funcionan.
    let mut document = document(&lines());
    let copy = copy_of(&mut document, "I3", 1);
    copy.segments[0].source.start = 300;
    copy.segments[0].source.end = 200;
    copy.history.clear();
    assert_eq!(judge(&document).err(), Some(replay_failure("I3", 1)));

    // 3. Una inserción en el último punto de la copia y una deleción hasta su
    //    última base: la réplica rechaza antes de partir ningún fragmento, de
    //    modo que la búsqueda del punto de corte siempre lo encuentra.
    for (label, history) in [
        (
            "inserción en el punto igual a la longitud",
            json!({"kind": "insertion", "event_id": "ins-x", "point": 600, "bases": "T",
                   "provenance": explicit("ins-1")}),
        ),
        (
            "deleción hasta la longitud",
            json!({"kind": "deletion", "event_id": "del-x", "start": 595, "end": 600,
                   "provenance": explicit("del-1")}),
        ),
    ] {
        let mut all = lines();
        let at = copy_line(&all, "I3", 1);
        let mut events = all[at]["history"].as_array().expect("historial").clone();
        events.push(history);
        all[at]["history"] = Value::Array(events);
        assert_eq!(
            full(jsonl(&all).as_bytes()).err(),
            Some(replay_failure("I3", 1)),
            "{label}"
        );
    }
}

// ---------------------------------------------------------------------------
// Los identificadores de evento heredados
// ---------------------------------------------------------------------------

/// Un evento del historial que comparte identificador con una modificación de
/// segmento sin ser su mismo cambio de base: con otra clase de evento, o con
/// el mismo identificador y otra posición, u otra procedencia. El único caso
/// admitido es el cambio de base de segmento y de historial que `with_snp`
/// conserva a propósito (la verdad publicada intacta lo lleva y se acepta).
#[test]
fn a_history_event_sharing_an_inherited_identifier_must_be_the_same_base_change() {
    let at = |all: &[Value]| copy_line(all, "I3", 1);
    let shared_id = |all: &[Value]| {
        all[at(all)]["segments"][0]["modifications"][0]["event_id"]
            .as_str()
            .expect("identificador")
            .to_owned()
    };
    // Una inserción del historial con el identificador de un cambio de segmento
    // (que sale del historial).
    let mut all = lines();
    let id = shared_id(&all);
    let line = at(&all);
    let history = all[line]["history"].as_array_mut().expect("historial");
    history.retain(|event| event["event_id"] != id.as_str());
    history.push(
        json!({"kind": "insertion", "event_id": id, "point": 5, "bases": "T",
                        "provenance": explicit("ins-1")}),
    );
    assert_both(
        "otra clase de evento",
        &bytes_of(&all),
        &replay_failure("I3", 1),
    );
    // El mismo identificador y otra posición.
    let mut all = lines();
    let line = at(&all);
    let id = shared_id(&all);
    for event in all[line]["history"].as_array_mut().expect("historial") {
        if event["event_id"] == id.as_str() {
            event["position"] = json!(60);
        }
    }
    assert_both("otra posición", &bytes_of(&all), &replay_failure("I3", 1));
    // El mismo identificador y otra procedencia.
    let mut all = lines();
    let line = at(&all);
    let id = shared_id(&all);
    for event in all[line]["history"].as_array_mut().expect("historial") {
        if event["event_id"] == id.as_str() {
            event["provenance"] = explicit("snp-direct");
        }
    }
    assert_both(
        "otra procedencia",
        &bytes_of(&all),
        &replay_failure("I3", 1),
    );
    // Una deleción de segmento con el identificador de un cambio del historial.
    let mut all = lines();
    let line = at(&all);
    let id = shared_id(&all);
    all[line]["segments"][0]["modifications"]
        .as_array_mut()
        .expect("modificaciones")
        .push(
            json!({"kind": "deletion", "event_id": id, "start": 590, "end": 595,
                     "provenance": explicit("del-1")}),
        );
    assert_both(
        "clase estructural",
        &bytes_of(&all),
        &replay_failure("I3", 1),
    );
}

// ---------------------------------------------------------------------------
// Documentos sin cabecera, y documentos con la copia en otro orden de fuentes
// ---------------------------------------------------------------------------

/// Un documento que solo trae copias, o solo casos, y ninguna cabecera: el
/// primer registro ya no es una cabecera (no es el caso de «la cabecera
/// llega después», que también se rechaza por no ser la primera).
#[test]
fn a_document_with_no_header_at_all_is_rejected_by_its_first_record() {
    let all = lines();
    let copies: Vec<Value> = all
        .iter()
        .filter(|l| l["record_type"] == "haplotype_copy")
        .cloned()
        .collect();
    let genotypes: Vec<Value> = all
        .iter()
        .filter(|l| l["record_type"] == "genotype")
        .cloned()
        .collect();
    assert_both(
        "solo copias",
        &bytes_of(&copies),
        &TruthValidationError::HeaderMustBeFirst,
    );
    assert_both(
        "solo casos",
        &bytes_of(&genotypes),
        &TruthValidationError::HeaderMustBeFirst,
    );
}

/// Un cambio de base heredado que empieza justo donde acaba el fragmento
/// anterior (el cursor coincide con el inicio): no hay hueco ni solapamiento y
/// el documento se acepta. Un segmento corta la copia en el 49 y el siguiente
/// empieza con el cambio.
#[test]
fn an_inherited_base_change_at_the_start_of_a_segment_is_accepted() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let history = all[at]["history"].as_array().expect("historial").clone();
    // El cambio del 49 deja el historial y pasa a heredado; el del 499 se queda.
    let first_id = mods[0]["event_id"].clone();
    let kept: Vec<Value> = history
        .iter()
        .filter(|e| e["event_id"] != first_id)
        .cloned()
        .collect();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    assert_eq!(snps[0]["position"], 49);
    set_copy(
        &mut all,
        "I3",
        1,
        vec![segment(0, 49, vec![]), segment(49, 600, mods)],
        vec![
            source_piece(0, 49, vec![]),
            source_piece(49, 50, vec![snps[0].clone()]),
            source_piece(50, 600, vec![snps[1].clone()]),
        ],
        kept,
    );
    for result in both(&bytes_of(&all)) {
        result.expect("el cambio al inicio del segmento se acepta");
    }
}

// ---------------------------------------------------------------------------
// Motivos exactos de rechazos que ya se disparaban
// ---------------------------------------------------------------------------

/// Una huella de referencia distinta (de forma válida): se rechaza con su
/// motivo, igual en los dos lectores.
#[test]
fn a_reference_digest_that_is_not_the_references_is_rejected() {
    let mut all = lines();
    all[0]["reference_content_sha256"] = json!("0".repeat(64));
    assert_both(
        "huella",
        &bytes_of(&all),
        &TruthValidationError::ReferenceDigestMismatch,
    );
}

/// Un cambio de base del historial sin su modificación de segmento, con una
/// base de referencia que no es la de la referencia, y con el fragmento
/// vigente que no es el que sale de la réplica.
#[test]
fn a_base_change_that_does_not_match_the_segments_or_the_reference_is_rejected() {
    // Sin modificación: el historial nombra otro evento.
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    all[at]["history"][0]["event_id"] = json!("EVT_DETACHED");
    assert_both(
        "sin modificación",
        &bytes_of(&all),
        &replay_failure("I3", 1),
    );
    // La modificación dice una referencia que no es la de la referencia.
    let mut all = lines();
    let modification = &mut all[at]["segments"][0]["modifications"][0];
    let real = modification["reference"]
        .as_str()
        .expect("referencia")
        .to_owned();
    let alternate = modification["alternate"]
        .as_str()
        .expect("alterna")
        .to_owned();
    let wrong = ["A", "C", "G", "T"]
        .into_iter()
        .find(|base| *base != real && *base != alternate)
        .expect("una base distinta de las dos");
    modification["reference"] = json!(wrong);
    let event_id = modification["event_id"]
        .as_str()
        .expect("evento")
        .to_owned();
    let position = modification["position"].as_u64().expect("posición");
    assert_both(
        "referencia equivocada",
        &bytes_of(&all),
        &TruthValidationError::SnpReferenceMismatch { event_id, position },
    );
    // El fragmento vigente lleva otra base alterna que la que resulta.
    let mut all = lines();
    let snp = &mut all[at]["current_pieces"][0]["snps"][0];
    let alternate = snp["alternate"].as_str().expect("alterna").to_owned();
    let reference_base = snp["reference"].as_str().expect("referencia").to_owned();
    let changed = ["A", "C", "G", "T"]
        .into_iter()
        .find(|base| *base != alternate && *base != reference_base)
        .expect("una base distinta de las dos");
    snp["alternate"] = json!(changed);
    assert_both(
        "fragmento vigente distinto",
        &bytes_of(&all),
        &replay_failure("I3", 1),
    );
}

/// Una inserción al final de la copia, justo en el punto igual a su longitud,
/// no se admite (tiene que quedar una base de referencia detrás). El documento
/// lleva su fragmento vigente, de modo que sin la guarda lo pararía otra cosa
/// distinta (el caso que falta).
#[test]
fn an_insertion_at_the_very_end_of_the_copy_does_not_replay() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    let mut history = all[at]["history"].as_array().expect("historial").clone();
    history.push(
        json!({"kind": "insertion", "event_id": "ins-end", "point": 600, "bases": "T",
                        "provenance": explicit("ins-1")}),
    );
    let mut pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    pieces.push(inserted_piece("T", &explicit("ins-1")));
    let segments = all[at]["segments"].as_array().expect("segmentos").clone();
    set_copy(&mut all, "I3", 1, segments, pieces, history);
    assert_both(
        "inserción al final",
        &bytes_of(&all),
        &replay_failure("I3", 1),
    );
}

/// Un evento del historial con el identificador de un cambio de base heredado
/// pero de otra clase: el documento está armado con los fragmentos que
/// resultarían de aceptarlo, y sin la guarda lo pararía el caso que falta.
#[test]
fn an_insertion_reusing_an_inherited_snp_identifier_does_not_replay() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    let shared = mods[0]["event_id"].clone();
    set_copy(
        &mut all,
        "I3",
        1,
        vec![segment(0, 600, mods)],
        vec![
            source_piece(0, 5, vec![]),
            inserted_piece("T", &explicit("ins-1")),
            source_piece(5, 49, vec![]),
            source_piece(49, 50, vec![snps[0].clone()]),
            source_piece(50, 499, vec![]),
            source_piece(499, 500, vec![snps[1].clone()]),
            source_piece(500, 600, vec![]),
        ],
        vec![
            json!({"kind": "insertion", "event_id": shared, "point": 5, "bases": "T",
                    "provenance": explicit("ins-1")}),
        ],
    );
    assert_both(
        "identificador compartido",
        &bytes_of(&all),
        &replay_failure("I3", 1),
    );
}

/// Un caso SNP cuyos alelos no incluyen el que una copia lleva.
#[test]
fn a_snp_case_missing_an_observed_allele_is_rejected() {
    let mut all = lines();
    let at = case_line(&all, "I1", "snp", 49);
    all[at]["alleles"] = json!([{"kind": "reference"}]);
    all[at]["allele_dosages"] = json!([2]);
    for carrier in all[at]["carriers"].as_array_mut().expect("portadores") {
        carrier["allele_index"] = json!(0);
    }
    let expected = genotype_mismatch(
        &all,
        at,
        "un alelo SNP reconstruido no aparece en el genotipo",
    );
    assert_both("alelo observado ausente", &bytes_of(&all), &expected);
}

/// El locus de un caso indel tiene que ser el que su evento ocupa en la
/// referencia, para una inserción (el ancla) y para una deleción (su intervalo).
#[test]
fn an_indel_locus_must_be_where_its_event_sits_in_the_reference() {
    let reason = "el locus del caso no es el que su evento ocupa en la referencia";
    let mut all = lines();
    let at = case_line(&all, "I1", "indel", 350);
    all[at]["locus"]["start"] = json!(349);
    all[at]["locus"]["end"] = json!(349);
    assert_both(
        "inserción en otro ancla",
        &bytes_of(&all),
        &genotype_mismatch(&all, at, reason),
    );
    let mut all = lines();
    let at = case_line(&all, "I1", "indel", 420);
    all[at]["locus"]["start"] = json!(421);
    assert_both(
        "deleción en otro intervalo",
        &bytes_of(&all),
        &genotype_mismatch(&all, at, reason),
    );
}

/// Dos casos del individuo para el mismo locus, y un caso indel de un evento
/// que ninguna copia lleva.
#[test]
fn a_case_published_twice_or_without_a_source_is_rejected() {
    // El mismo locus SNP dos veces, con identificadores distintos.
    let mut all = lines();
    let at = case_line(&all, "I1", "snp", 49);
    let mut twin = all[at].clone();
    // El mismo identificador con el último carácter cambiado por uno mayor.
    let original = twin["case_id"].as_str().expect("identificador").to_owned();
    twin["case_id"] = json!(format!("{}z", &original[..original.len() - 1]));
    all.insert(at + 1, twin);
    let count = all[0]["counts"]["genotypes"].as_u64().expect("recuento");
    all[0]["counts"]["genotypes"] = json!(count + 1);
    let expected = genotype_mismatch(
        &all,
        at + 1,
        "otro caso del individuo publica ya ese locus o ese evento",
    );
    assert_both("caso repetido", &bytes_of(&all), &expected);
    // Una inserción que ninguna copia lleva: el evento no está en ningún historial.
    let mut all = lines();
    let at = case_line(&all, "I2", "indel", 100);
    let mut orphan = all[at].clone();
    orphan["case_id"] = json!("CASE_zz");
    orphan["alleles"][1]["event_id"] = json!("ins-zz");
    orphan["alleles"][1]["source_event"]["event_id"] = json!("ins-zz");
    orphan["allele_dosages"] = json!([2, 0]);
    for carrier in orphan["carriers"].as_array_mut().expect("portadores") {
        carrier["allele_index"] = json!(0);
    }
    all.insert(at + 1, orphan);
    let count = all[0]["counts"]["genotypes"].as_u64().expect("recuento");
    all[0]["counts"]["genotypes"] = json!(count + 1);
    let expected = genotype_mismatch(
        &all,
        at + 1,
        "el caso no corresponde a ningún cambio de las copias, evento del historial ni fila de la matriz",
    );
    assert_both("caso sin fuente", &bytes_of(&all), &expected);
}

/// Una estructura escrita como lista (`"counts":[6,8]`) se rechaza con su
/// motivo, la lista donde se esperaba un objeto.
#[test]
fn a_structure_written_as_a_list_is_rejected_with_its_reason() {
    let mut all = lines();
    all[0]["counts"] = json!([6, 8]);
    assert_schema_detail_prefix(
        "counts como lista",
        &jsonl(&all),
        "invalid type: sequence, expected an object, not a list of its values",
    );
}

// ---------------------------------------------------------------------------
// Fronteras entre fragmentos y alelos que se distinguen por una sola mitad
// ---------------------------------------------------------------------------

/// Un cambio de base en la primera base del fragmento que lo recibe (el
/// fragmento anterior acaba justo donde ese empieza): se aplica al segundo, no
/// al primero ni a los dos. La copia 1 de I3 pasa a dos segmentos que se
/// cortan en la base 499, donde está su cambio.
#[test]
fn a_base_change_on_the_first_base_of_a_fragment_belongs_to_that_fragment() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    assert_eq!(snps[1]["position"], 499);
    let history = all[at]["history"].as_array().expect("historial").clone();
    set_copy(
        &mut all,
        "I3",
        1,
        vec![
            segment(0, 499, vec![mods[0].clone()]),
            segment(499, 600, vec![mods[1].clone()]),
        ],
        vec![
            source_piece(0, 499, vec![snps[0].clone()]),
            source_piece(499, 600, vec![snps[1].clone()]),
        ],
        history,
    );
    for result in both(&bytes_of(&all)) {
        result.expect("el cambio es del fragmento que empieza en él");
    }
}

/// Un tracto de repetición que acaba justo donde acaba su fragmento de fuente
/// (no hay base a la derecha): se aceptan el fragmento y el tracto. El documento
/// es el publicado con el fragmento izquierdo cortado en dos segmentos en la
/// base 300.
#[test]
fn a_repeat_tract_at_the_end_of_its_fragment_is_accepted() {
    let mut all = lines_v5();
    let at = copy_line(&all, "I3", 0);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let history = all[at]["history"].as_array().expect("historial").clone();
    assert_eq!(pieces.len(), 3);
    set_copy(
        &mut all,
        "I3",
        0,
        vec![segment(0, 300, mods), segment(300, 600, vec![])],
        pieces,
        history,
    );
    full(jsonl(&all).as_bytes()).expect("tracto al final de su fragmento");
}

/// Un tracto de repetición que empieza justo donde empieza su fragmento de
/// fuente (no hay base a la izquierda): sustituye el fragmento entero por la
/// izquierda y se acepta. El documento es el publicado con el fragmento
/// izquierdo partido en dos segmentos en la base 200.
#[test]
fn a_repeat_tract_at_the_start_of_its_fragment_is_accepted() {
    let mut all = lines_v5();
    let at = copy_line(&all, "I3", 0);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let history = all[at]["history"].as_array().expect("historial").clone();
    assert_eq!(pieces.len(), 3);
    set_copy(
        &mut all,
        "I3",
        0,
        vec![segment(0, 200, mods), segment(200, 600, vec![])],
        pieces,
        history,
    );
    full(jsonl(&all).as_bytes()).expect("tracto al inicio de su fragmento");
}

/// Un tracto cuyo intervalo contiene un cambio de base del fragmento que lo
/// recibiría: la réplica no sabe qué hacer con él y rechaza. El documento está
/// armado con los fragmentos que saldrían de ignorar el cambio (se pierde), de
/// modo que sin la comprobación lo pararía el caso del cambio.
#[test]
fn a_repeat_tract_over_a_base_change_does_not_replay() {
    let mut all = lines_v5();
    let at = copy_line(&all, "I3", 0);
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    assert_eq!(snps[1]["position"], 119);
    let mut tract = pieces[1].clone();
    tract["source"]["start"] = json!(100);
    tract["source"]["end"] = json!(150);
    let mut history = all[at]["history"].as_array().expect("historial").clone();
    for event in &mut history {
        if event["kind"] == "repeat_allele" {
            event["start"] = json!(100);
            event["end"] = json!(150);
        }
    }
    let segments = all[at]["segments"].as_array().expect("segmentos").clone();
    set_copy(
        &mut all,
        "I3",
        0,
        segments,
        vec![
            source_piece(0, 100, vec![snps[0].clone()]),
            tract,
            source_piece(150, 600, vec![]),
        ],
        history,
    );
    assert_eq!(
        full(jsonl(&all).as_bytes()).err(),
        Some(replay_failure("I3", 0))
    );
}

/// Un caso SNP cuyo locus cae dentro de un tracto realizado: ninguna copia
/// lleva ahí una base de referencia que leer, de modo que ninguna es portadora
/// y el caso, sin portadores, se acepta.
#[test]
fn a_snp_case_inside_a_realized_tract_has_no_carriers() {
    let mut all = lines_v5();
    let at = case_line(&all, "I1", "repeat", 200);
    let snp = json!({
        "record_type": "genotype",
        "case_id": "CASE_ZZ_INSIDE_TRACT",
        "reference_id": "ref-main",
        "individual_id": "I1",
        "case_kind": "snp",
        "locus": {"sequence_id": "chr1", "start": 250, "end": 251},
        "alleles": [{"kind": "reference"}, {"kind": "snp", "base": "A"}],
        "allele_dosages": [0, 0],
        "carriers": [],
    });
    all.insert(at + 1, snp);
    let count = all[0]["counts"]["genotypes"].as_u64().expect("recuento");
    all[0]["counts"]["genotypes"] = json!(count + 1);
    full(jsonl(&all).as_bytes()).expect("un caso sin portadores dentro de un tracto");
}

/// Dos alelos de repetición del mismo caso que comparten la secuencia y se
/// distinguen solo por su estructura de bloques (o al revés): cada copia se
/// identifica con el alelo que coincide en **las dos** cosas.
#[test]
fn a_repeat_allele_is_identified_by_its_sequence_and_its_blocks() {
    // (qué comparten los dos alelos alternativos, cómo se distingue el segundo)
    for shared_sequence in [true, false] {
        let mut all = lines_v5();
        let case = case_line(&all, "I2", "repeat", 200);
        let first = all[case]["alleles"][1].clone();
        let second = all[case]["alleles"][2].clone();
        assert_ne!(first["bases"], second["bases"]);
        assert_eq!(first["bases"].as_str().expect("bases").len() % 2, 0);
        // El segundo alelo conserva la secuencia del primero (o su estructura)
        // y cambia la otra mitad; la copia que lo lleva se ajusta igual.
        let mut distinct = second.clone();
        if shared_sequence {
            distinct["bases"] = first["bases"].clone();
            let total = first["bases"].as_str().expect("bases").len() as u64;
            distinct["blocks"] = json!([
                {"motif": "CA", "motif_length": 2, "bases": total - 2, "units": number(&format!("{}.0", (total - 2) / 2))},
                {"motif": "CA", "motif_length": 2, "bases": 2, "units": number("1.0")},
            ]);
        } else {
            // La misma estructura de bloques que el primero, con otra
            // secuencia de la misma longitud.
            let mut text = first["bases"].as_str().expect("bases").to_owned();
            text.replace_range(0..4, "CATT");
            distinct["bases"] = json!(text);
            distinct["blocks"] = first["blocks"].clone();
        }
        all[case]["alleles"][2] = distinct.clone();
        let copy = copy_line(&all, "I2", 1);
        for field in ["current_pieces", "history"] {
            for item in all[copy][field].as_array_mut().expect("lista") {
                if item["kind"] == "repeat_allele" {
                    item["bases"] = distinct["bases"].clone();
                    item["blocks"] = distinct["blocks"].clone();
                }
            }
        }
        full(jsonl(&all).as_bytes())
            .unwrap_or_else(|error| panic!("compartida la secuencia {shared_sequence}: {error:?}"));
    }
}

/// La tolerancia del número de unidades es **absoluta** cuando las unidades
/// son menos de una (un bloque más corto que su motivo): una milmillonésima.
/// Con unidades `2/3`, un error de 8e-10 se admite y uno de 1,2e-9 no.
#[test]
fn the_unit_tolerance_is_absolute_for_a_block_shorter_than_its_motif() {
    let blocks_with = |error: f64| {
        let mut all = lines_v5();
        let path = block_places(&all)[0].1.clone();
        let line: usize = path[0].parse().expect("línea");
        let piece: usize = path[2].parse().expect("fragmento");
        let total = all[line]["current_pieces"][piece]["bases"]
            .as_str()
            .expect("bases")
            .len() as u64;
        all[line]["current_pieces"][piece]["blocks"] = json!([
            {"motif": "CAG", "motif_length": 3, "bases": 2, "units": 2.0 / 3.0 + error},
            {"motif": "CA", "motif_length": 2, "bases": total - 2, "units": number(&format!("{}.0", (total - 2) / 2))},
        ]);
        jsonl(&all)
    };
    read_truth_jsonl(&mut Cursor::new(blocks_with(8.0e-10))).expect("8e-10 se admite");
    let rejected = read_truth_jsonl(&mut Cursor::new(blocks_with(1.2e-9))).err();
    assert!(
        matches!(&rejected, Some(TruthValidationError::Schema { detail, .. })
            if detail.starts_with("el número de unidades no es las bases entre la longitud del motivo")),
        "1,2e-9 no se admite: {rejected:?}"
    );
}

/// Una copia de **otra referencia** con la misma identidad salvo ese campo no
/// cuenta para ningún caso: no es candidata a portadora ni de un cambio de base,
/// ni de una inserción, ni de un tracto. Solo se alcanza con un documento en
/// memoria que no ha pasado por el contrato (que exige que toda copia sea de la
/// referencia de la cabecera).
#[test]
fn a_copy_of_another_reference_is_not_a_candidate_carrier() {
    let mut all = document(&lines_v5());
    let original = copy_of(&mut all, "I1", 0).clone();
    let mut foreign = original;
    foreign.identity.reference_id = "ref-other".to_owned();
    all.records.insert(
        0,
        genome_population_simulator::truth_schema::TruthRecord::HaplotypeCopy(foreign),
    );
    judge(&all).expect("la copia ajena no cambia ningún caso");
}

// ---------------------------------------------------------------------------
// El estado heredado recorre la referencia una vez y en orden (DT-048)
// ---------------------------------------------------------------------------
//
// Los segmentos de una copia, en el orden publicado, recorren su secuencia de
// la referencia exactamente una vez y en orden creciente. Los documentos de
// abajo están armados con los fragmentos vigentes que salen de replicar esos
// segmentos y su historial, de modo que la réplica, la materialización y el
// conjunto de casos los dejan pasar y solo la comprobación del recorrido puede
// pararlos (salvo el hueco, que antes paraba el conjunto de casos con otro
// motivo).

fn not_collinear(
    individual: &str,
    copy_index: u8,
    reached: u64,
    next_start: u64,
) -> TruthValidationError {
    TruthValidationError::SegmentsNotCollinear {
        identity: identity(individual, copy_index),
        reached,
        next_start,
    }
}

/// Dos fragmentos de fuente que se tocan sin solaparse, en el orden de la
/// copia **inverso** al de la referencia: la copia 1 de I3 pasa a los segmentos
/// `[300, 600)` y `[0, 300)` —el cromosoma rotado—, con sus fragmentos vigentes
/// en ese orden. Hasta el 09/10/2026 lo aceptaban los dos lectores.
#[test]
fn touching_source_fragments_in_reverse_order_are_rejected() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    assert_eq!(
        (snps[0]["position"].as_u64(), snps[1]["position"].as_u64()),
        (Some(49), Some(499))
    );
    let history = all[at]["history"].as_array().expect("historial").clone();
    set_copy(
        &mut all,
        "I3",
        1,
        vec![
            segment(300, 600, vec![mods[1].clone()]),
            segment(0, 300, vec![mods[0].clone()]),
        ],
        vec![
            source_piece(300, 600, vec![snps[1].clone()]),
            source_piece(0, 300, vec![snps[0].clone()]),
        ],
        history,
    );
    assert_both("rotación", &bytes_of(&all), &not_collinear("I3", 1, 0, 300));
}

/// Lo mismo con un tracto de repetición que acaba donde empieza el fragmento
/// anterior (versión 5: solo el lector completo). Hasta el 09/10/2026 se
/// aceptaba.
#[test]
fn a_repeat_tract_touching_the_previous_fragment_is_rejected() {
    let mut all = lines_v5();
    let at = copy_line(&all, "I3", 0);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let history = all[at]["history"].as_array().expect("historial").clone();
    assert_eq!(pieces.len(), 3, "izquierda, tracto, derecha");
    assert_eq!(pieces[1]["kind"], "repeat_allele");
    set_copy(
        &mut all,
        "I3",
        0,
        vec![segment(300, 600, vec![]), segment(0, 300, mods)],
        vec![pieces[2].clone(), pieces[0].clone(), pieces[1].clone()],
        history,
    );
    assert_eq!(
        full(jsonl(&all).as_bytes()).err(),
        Some(not_collinear("I3", 0, 0, 300))
    );
}

/// La copia 1 de I1 lleva la deleción `del-1` (420..424) y su caso. Con los
/// segmentos `[0, 424)`, `[420, 424)` y `[424, 600)` el tramo 420..424 está dos
/// veces y la deleción quita una: la copia final es la referencia con su cambio
/// en 49, y aun así el caso de la deleción la declara portadora. Hasta el
/// 09/10/2026 lo aceptaban los dos lectores.
#[test]
fn a_duplicated_stretch_that_the_history_deletes_again_is_rejected() {
    let mut all = lines();
    let at = copy_line(&all, "I1", 1);
    let segments = all[at]["segments"].as_array().expect("segmentos").clone();
    assert_eq!(segments.len(), 1);
    let mods = segments[0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    assert_eq!(
        (
            pieces[0]["source"]["end"].as_u64(),
            pieces[1]["source"]["start"].as_u64()
        ),
        (Some(420), Some(424)),
        "la deleción del-1 quita 420..424"
    );
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    let history = all[at]["history"].as_array().expect("historial").clone();
    assert_eq!(history[1]["event_id"], "del-1");
    set_copy(
        &mut all,
        "I1",
        1,
        vec![
            segment(0, 424, mods),
            segment(420, 424, vec![]),
            segment(424, 600, vec![]),
        ],
        vec![
            source_piece(0, 420, snps),
            source_piece(420, 424, vec![]),
            source_piece(424, 600, vec![]),
        ],
        history,
    );
    assert_both(
        "tramo duplicado y borrado",
        &bytes_of(&all),
        &not_collinear("I1", 1, 424, 420),
    );
}

/// Un solape parcial que la deleción del historial esconde: los segmentos
/// `[0, 424)` y `[422, 600)` se solapan en 422..424, y `del-1` (420..424) quita
/// justo el final del primero. Los fragmentos vigentes, `[0, 420)` y
/// `[422, 600)`, ya no se solapan, de modo que la comprobación de fragmentos
/// solapados no lo ve; la copia pierde 420..422 y el caso dice 420..424. Hasta
/// el 09/10/2026 lo aceptaban los dos lectores.
#[test]
fn a_partial_overlap_hidden_by_a_deletion_is_rejected() {
    let mut all = lines();
    let at = copy_line(&all, "I1", 1);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    let history = all[at]["history"].as_array().expect("historial").clone();
    set_copy(
        &mut all,
        "I1",
        1,
        vec![segment(0, 424, mods), segment(422, 600, vec![])],
        vec![source_piece(0, 420, snps), source_piece(422, 600, vec![])],
        history,
    );
    assert_both(
        "solape parcial escondido",
        &bytes_of(&all),
        &not_collinear("I1", 1, 424, 422),
    );
}

/// Un hueco entre dos segmentos (`[0, 300)` y `[310, 600)`): la copia pierde
/// 300..310 sin deleción ni caso. Hasta el 09/10/2026 lo rechazaba el conjunto
/// de casos (bases de referencia que faltan en el estado heredado); ahora lo
/// rechaza la réplica, con dónde se rompe el recorrido.
#[test]
fn a_gap_between_segments_is_rejected() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    let mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    let pieces = all[at]["current_pieces"]
        .as_array()
        .expect("fragmentos")
        .clone();
    let snps = pieces[0]["snps"].as_array().expect("cambios").clone();
    let history = all[at]["history"].as_array().expect("historial").clone();
    set_copy(
        &mut all,
        "I3",
        1,
        vec![
            segment(0, 300, vec![mods[0].clone()]),
            segment(310, 600, vec![mods[1].clone()]),
        ],
        vec![
            source_piece(0, 300, vec![snps[0].clone()]),
            source_piece(310, 600, vec![snps[1].clone()]),
        ],
        history,
    );
    assert_both("hueco", &bytes_of(&all), &not_collinear("I3", 1, 300, 310));
}

/// Tres segmentos que cubren la secuencia entera sin solaparse pero fuera de
/// orden, en la **segunda** secuencia (`chr2`): `[0, 100)`, `[200, 400)` y
/// `[100, 200)`. Hasta el 09/10/2026 lo aceptaban los dos lectores (la
/// búsqueda de huecos ordenaba los segmentos antes de mirarlos).
#[test]
fn segments_out_of_order_on_another_sequence_are_rejected() {
    let mut all = lines();
    let at = all
        .iter()
        .position(|line| {
            line["record_type"] == "haplotype_copy"
                && line["identity"]["individual_id"] == "I2"
                && line["identity"]["sequence_id"] == "chr2"
                && line["identity"]["copy_index"] == 1
        })
        .expect("la copia 1 de I2 en chr2");
    assert_eq!(all[at]["segments"].as_array().map(Vec::len), Some(1));
    assert_eq!(
        all[at]["segments"][0]["modifications"],
        json!([]),
        "sin cambios heredados"
    );
    assert_eq!(all[at]["history"], json!([]), "sin historial");
    let piece = |start: u64, end: u64| json!({"kind": "source", "source": reference_object_on("chr2", start, end, 400), "snps": []});
    let segment_on = |start: u64, end: u64| json!({"source": reference_object_on("chr2", start, end, 400), "modifications": []});
    all[at]["segments"] = json!([
        segment_on(0, 100),
        segment_on(200, 400),
        segment_on(100, 200)
    ]);
    all[at]["current_pieces"] = json!([piece(0, 100), piece(200, 400), piece(100, 200)]);
    let mut expected_identity = identity("I2", 1);
    expected_identity.sequence_id = "chr2".to_owned();
    assert_both(
        "desorden en chr2",
        &bytes_of(&all),
        &TruthValidationError::SegmentsNotCollinear {
            identity: expected_identity,
            reached: 100,
            next_start: 200,
        },
    );
}

fn reference_object_on(sequence: &str, start: u64, end: u64, length: u64) -> Value {
    let mut object = reference_object(start, end, length);
    object["sequence_id"] = json!(sequence);
    object
}

/// En memoria (sin el contrato, que exige segmentos): una copia sin ningún
/// segmento no recorre nada, y el recorrido se rompe en la longitud de la
/// secuencia. Si además su secuencia no está en la referencia, el rechazo es el
/// de la secuencia ausente: sin ningún fragmento, la materialización no llegó
/// a mirarla.
#[test]
fn a_copy_without_segments_does_not_walk_the_reference() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I3", 1);
    copy.segments.clear();
    copy.current_pieces.clear();
    copy.history.clear();
    assert_eq!(judge(&all).err(), Some(not_collinear("I3", 1, 0, 600)));

    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I3", 1);
    copy.segments.clear();
    copy.current_pieces.clear();
    copy.history.clear();
    copy.identity.sequence_id = "chr9".to_owned();
    assert_eq!(
        judge(&all).err(),
        Some(TruthValidationError::MissingReferenceSequence {
            sequence_id: "chr9".to_owned(),
        })
    );
}

/// Parte en dos el fragmento de fuente que contiene `point` (sin tocar sus
/// extremos), con cada cambio de base en su mitad. Devuelve también el índice
/// del fragmento partido.
fn split_source_piece_at(pieces: &[Value], point: u64) -> (Vec<Value>, u64) {
    let mut split = Vec::new();
    let mut found = Vec::new();
    for (index, piece) in pieces.iter().enumerate() {
        let start = piece["source"]["start"].as_u64();
        let end = piece["source"]["end"].as_u64();
        match (piece["kind"].as_str(), start, end) {
            (Some("source"), Some(start), Some(end)) if start < point && point < end => {
                found.push(u64::try_from(index).expect("índice"));
                let snps = piece["snps"].as_array().expect("cambios");
                let side = |left: bool| {
                    snps.iter()
                        .filter(|snp| (snp["position"].as_u64().expect("posición") < point) == left)
                        .cloned()
                        .collect::<Vec<_>>()
                };
                let mut left = piece.clone();
                left["source"]["end"] = json!(point);
                left["snps"] = Value::Array(side(true));
                let mut right = piece.clone();
                right["source"]["start"] = json!(point);
                right["snps"] = Value::Array(side(false));
                split.push(left);
                split.push(right);
            }
            _ => split.push(piece.clone()),
        }
    }
    assert_eq!(
        found.len(),
        1,
        "un único fragmento de fuente contiene {point}"
    );
    (split, found[0])
}

/// Tras partir el fragmento `split` de la copia, los enlaces de las inserciones
/// que lleva esa copia a fragmentos posteriores avanzan un puesto. Devuelve
/// cuántos enlaces movió.
fn shift_piece_links(all: &mut [Value], individual: &str, copy_index: u8, split: u64) -> usize {
    let mut moved = 0;
    for line in all.iter_mut() {
        if line["record_type"] != "genotype" {
            continue;
        }
        let carried = line["carriers"]
            .as_array()
            .expect("portadores")
            .iter()
            .filter(|carrier| {
                carrier["copy"]["individual_id"] == individual
                    && carrier["copy"]["sequence_id"] == "chr1"
                    && carrier["copy"]["copy_index"] == copy_index
            })
            .map(|carrier| carrier["allele_index"].as_u64().expect("alelo"))
            .collect::<Vec<_>>();
        for allele_index in carried {
            let carriers_of_allele = line["carriers"]
                .as_array()
                .expect("portadores")
                .iter()
                .filter(|carrier| carrier["allele_index"] == allele_index)
                .count();
            let allele = &mut line["alleles"][usize::try_from(allele_index).expect("índice")];
            let Some(links) = allele.get_mut("piece_links").and_then(Value::as_array_mut) else {
                continue;
            };
            if links.is_empty() {
                continue;
            }
            assert_eq!(
                carriers_of_allele, 1,
                "una inserción del escenario la lleva una copia"
            );
            for link in links {
                let index = link["piece_index"].as_u64().expect("índice");
                if index > split {
                    link["piece_index"] = json!(index + 1);
                    moved += 1;
                }
            }
        }
    }
    moved
}

/// El lado aceptado: segmentos colineales partidos en tres, sobre copias con
/// cambios de base heredados y del historial, inserciones, deleciones (también
/// una que se come bases insertadas y de referencia) y, en la versión 5, el
/// tracto. Cada copia de `chr1` pasa a los segmentos
/// `[0, 150)`, `[150, 375)` y `[375, 600)`, con sus modificaciones repartidas,
/// los fragmentos vigentes partidos en los mismos puntos y los enlaces de sus
/// inserciones renumerados. Los dos lectores (el completo en la versión 5) los
/// aceptan.
#[test]
fn collinear_segments_cut_anywhere_are_accepted() {
    let cut = |all: &mut Vec<Value>, individual: &str, copy_index: u8| -> usize {
        let at = copy_line(all, individual, copy_index);
        let mods = all[at]["segments"][0]["modifications"]
            .as_array()
            .expect("modificaciones")
            .clone();
        let in_range = |low: u64, high: u64| {
            mods.iter()
                .filter(|m| {
                    let place = m["position"]
                        .as_u64()
                        .or_else(|| m["point"].as_u64())
                        .or_else(|| m["start"].as_u64())
                        .expect("lugar");
                    low <= place && place < high
                })
                .cloned()
                .collect::<Vec<_>>()
        };
        let pieces = all[at]["current_pieces"]
            .as_array()
            .expect("fragmentos")
            .clone();
        let (pieces, first) = split_source_piece_at(&pieces, 150);
        let (pieces, second) = split_source_piece_at(&pieces, 375);
        all[at]["segments"] = json!([
            segment(0, 150, in_range(0, 150)),
            segment(150, 375, in_range(150, 375)),
            segment(375, 600, in_range(375, 600)),
        ]);
        all[at]["current_pieces"] = Value::Array(pieces);
        shift_piece_links(all, individual, copy_index, first)
            + shift_piece_links(all, individual, copy_index, second)
    };
    let copies = [
        ("I1", 0),
        ("I1", 1),
        ("I2", 0),
        ("I2", 1),
        ("I3", 0),
        ("I3", 1),
    ];
    for (version, published) in [("4", lines()), ("5", lines_v5())] {
        let mut all = published.clone();
        let moved: usize = copies
            .iter()
            .map(|(individual, copy_index)| cut(&mut all, individual, *copy_index))
            .sum();
        assert!(
            moved > 0,
            "versión {version}: algún enlace de inserción se movió"
        );
        assert_ne!(all, published, "versión {version}: la alteración alteró");
        let whole = full(jsonl(&all).as_bytes());
        whole.unwrap_or_else(|error| panic!("versión {version}: {error:?}"));
        if version == "4" {
            for result in both(&bytes_of(&all)) {
                result.expect("segmentos colineales partidos");
            }
        }
    }
}

/// La copia de ese individuo en `chr2`, del documento en memoria.
fn chr2_copy_of<'a>(
    document: &'a mut genome_population_simulator::truth_schema::TruthDocument,
    individual: &str,
    copy_index: u8,
) -> &'a mut genome_population_simulator::truth_schema::HaplotypeCopyRecord {
    document
        .records
        .iter_mut()
        .find_map(|record| match record {
            genome_population_simulator::truth_schema::TruthRecord::HaplotypeCopy(copy)
                if copy.identity.individual_id == individual
                    && copy.identity.sequence_id == "chr2"
                    && copy.identity.copy_index == copy_index =>
            {
                Some(copy)
            }
            _ => None,
        })
        .expect("la copia en chr2")
}

fn chr2_identity(
    individual: &str,
    copy_index: u8,
) -> genome_population_simulator::truth_schema::CopyIdentity {
    let mut identity = identity(individual, copy_index);
    identity.sequence_id = "chr2".to_owned();
    identity
}

fn source_on(
    sequence: &str,
    start: u64,
    end: u64,
    length: u64,
) -> genome_population_simulator::truth_schema::SegmentSource {
    let mut source = typed_source(start, end, length);
    sequence.clone_into(&mut source.sequence_id);
    source
}

/// En memoria (el contrato lo rechaza por bytes): un segmento **invertido**,
/// `[200, 100)`, entre `[0, 200)` y `[100, 400)` de la copia 1 de I2 en `chr2`,
/// con dos deleciones de segmento (200..100 y 100..200) que reponen el cursor
/// de la réplica, de modo que la réplica lo deja pasar con los fragmentos
/// `[0, 200)` y `[200, 400)`. La regla lo rechaza en el segmento invertido
/// aunque empiece donde iba el recorrido, sin depender de otra guarda. Antes de
/// endurecerla, lo paraba el conjunto de casos (una deleción de segmento sin
/// caso) o, en las fronteras del VCF, que no juzgan el conjunto, pasaba.
#[test]
fn an_inverted_segment_that_the_replay_lets_through_is_rejected() {
    let mut all = document(&lines());
    let copy = chr2_copy_of(&mut all, "I2", 1);
    assert!(copy.history.is_empty(), "sin historial");
    assert_eq!(copy.segments.len(), 1);
    assert!(
        copy.segments[0].modifications.is_empty(),
        "sin cambios heredados"
    );
    let deletion = |event_id: &str, start: u64, end: u64| TruthModification::Deletion {
        event_id: event_id.to_owned(),
        start,
        end,
        provenance: typed_provenance("del-1"),
    };
    copy.segments = vec![
        TruthSegment {
            source: source_on("chr2", 0, 200, 400),
            modifications: vec![],
        },
        TruthSegment {
            source: source_on("chr2", 200, 100, 400),
            modifications: vec![deletion("DEL_BACK", 200, 100)],
        },
        TruthSegment {
            source: source_on("chr2", 100, 400, 400),
            modifications: vec![deletion("DEL_FORWARD", 100, 200)],
        },
    ];
    copy.current_pieces = [(0, 200), (200, 400)]
        .into_iter()
        .map(|(start, end)| TruthCurrentPiece::Source {
            source: source_on("chr2", start, end, 400),
            snps: vec![],
        })
        .collect();
    assert_eq!(
        judge(&all).err(),
        Some(TruthValidationError::SegmentsNotCollinear {
            identity: chr2_identity("I2", 1),
            reached: 200,
            next_start: 200,
        })
    );
}

/// En memoria (el contrato lo rechaza por bytes): un segmento **vacío**,
/// `[300, 300)`, entre `[0, 300)` y `[300, 600)` de la copia 1 de I3. La
/// réplica lo ignora y los fragmentos son los de cortar la copia en 300, así
/// que solo la regla lo ve: un segmento que no avanza el recorrido no es un
/// segmento del recorrido.
#[test]
fn an_empty_segment_is_rejected() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I3", 1);
    let mods = copy.segments[0].modifications.clone();
    assert_eq!(mods.len(), 2, "cambios heredados en 49 y 499");
    let TruthCurrentPiece::Source { snps, .. } = copy.current_pieces[0].clone() else {
        panic!("un fragmento de fuente");
    };
    assert_eq!(snps.len(), 2);
    copy.segments = vec![
        typed_segment(0, 300, vec![mods[0].clone()]),
        typed_segment(300, 300, vec![]),
        typed_segment(300, 600, vec![mods[1].clone()]),
    ];
    copy.current_pieces = vec![
        typed_piece(0, 300, vec![snps[0].clone()]),
        typed_piece(300, 600, vec![snps[1].clone()]),
    ];
    assert_eq!(judge(&all).err(), Some(not_collinear("I3", 1, 300, 300)));
}

/// En memoria (por bytes lo rechaza el contrato con `SegmentSourceMismatch`):
/// una copia de `chr2` con un segmento de `chr1` en medio, `chr2[0, 100)`,
/// `chr1[100, 300)`, `chr2[300, 400)`. Las coordenadas encajan, pero el
/// segmento del medio no recorre la secuencia de la copia (DT-050): hasta el
/// 09/10/2026 se aceptaba, y la copia reconstruida no era `chr2`.
#[test]
fn a_segment_of_another_sequence_is_rejected() {
    let mut all = document(&lines());
    let copy = chr2_copy_of(&mut all, "I2", 1);
    assert!(copy.history.is_empty(), "sin historial");
    let parts = [
        ("chr2", 0, 100, 400),
        ("chr1", 100, 300, 600),
        ("chr2", 300, 400, 400),
    ];
    copy.segments = parts
        .iter()
        .map(|&(sequence, start, end, length)| TruthSegment {
            source: source_on(sequence, start, end, length),
            modifications: vec![],
        })
        .collect();
    copy.current_pieces = parts
        .iter()
        .map(
            |&(sequence, start, end, length)| TruthCurrentPiece::Source {
                source: source_on(sequence, start, end, length),
                snps: vec![],
            },
        )
        .collect();
    assert_eq!(
        judge(&all).err(),
        Some(TruthValidationError::SegmentsNotCollinear {
            identity: chr2_identity("I2", 1),
            reached: 100,
            next_start: 100,
        })
    );
}

/// Un segmento que **vuelve a empezar en 0**: la copia 1 de I2 en `chr2` pasa a
/// `[0, 50)` y `[0, 400)`, con una deleción del historial que quita el tramo
/// repetido (vigentes 50..100). Los fragmentos vigentes, `[0, 50)` y
/// `[50, 400)`, parecen la referencia entera. Distingue una regla que tolerara
/// reiniciar el recorrido en 0.
#[test]
fn a_segment_restarting_at_zero_is_rejected() {
    let mut all = lines();
    let at = all
        .iter()
        .position(|line| {
            line["record_type"] == "haplotype_copy"
                && line["identity"]["individual_id"] == "I2"
                && line["identity"]["sequence_id"] == "chr2"
                && line["identity"]["copy_index"] == 1
        })
        .expect("la copia 1 de I2 en chr2");
    assert_eq!(all[at]["history"], json!([]), "sin historial");
    let on_chr2 = |start: u64, end: u64| reference_object_on("chr2", start, end, 400);
    all[at]["segments"] = json!([
        {"source": on_chr2(0, 50), "modifications": []},
        {"source": on_chr2(0, 400), "modifications": []},
    ]);
    all[at]["current_pieces"] = json!([
        {"kind": "source", "source": on_chr2(0, 50), "snps": []},
        {"kind": "source", "source": on_chr2(50, 400), "snps": []},
    ]);
    all[at]["history"] = json!([{
        "kind": "deletion", "event_id": "del-dup", "start": 50, "end": 100,
        "provenance": explicit("del-1"),
    }]);
    assert_both(
        "reinicio en 0",
        &bytes_of(&all),
        &TruthValidationError::SegmentsNotCollinear {
            identity: chr2_identity("I2", 1),
            reached: 50,
            next_start: 0,
        },
    );
}

/// Una deleción de segmento que va **después** de otras modificaciones del
/// mismo segmento: la copia 1 de I3 conserva sus dos cambios heredados (49 y
/// 499) y añade la deleción de segmento 550..560, sin evento en el historial ni
/// caso. Los dos lectores la rechazan como bases de referencia que faltan en el
/// estado heredado. Distingue un juez que solo mirara la primera modificación de
/// cada segmento.
#[test]
fn a_segment_deletion_after_other_inherited_changes_has_no_case() {
    let mut all = lines();
    let at = copy_line(&all, "I3", 1);
    let mut mods = all[at]["segments"][0]["modifications"]
        .as_array()
        .expect("modificaciones")
        .clone();
    assert_eq!(mods.len(), 2, "cambios heredados en 49 y 499");
    mods.push(
        json!({"kind": "deletion", "event_id": "del-seg", "start": 550, "end": 560,
                     "provenance": explicit("del-1")}),
    );
    let snps = all[at]["current_pieces"][0]["snps"]
        .as_array()
        .expect("cambios")
        .clone();
    let history = all[at]["history"].as_array().expect("historial").clone();
    set_copy(
        &mut all,
        "I3",
        1,
        vec![segment(0, 600, mods)],
        vec![source_piece(0, 550, snps), source_piece(560, 600, vec![])],
        history,
    );
    assert_both(
        "deleción de segmento tras otros cambios",
        &bytes_of(&all),
        &TruthValidationError::MissingCase {
            individual_id: "I3".to_owned(),
            sequence_id: "chr1".to_owned(),
            start: 550,
            end: 560,
            detail: "faltan bases de referencia en el estado heredado de los segmentos, sin evento ni caso",
        },
    );
}

/// En memoria (por bytes el contrato exige que el final no pase de la
/// longitud): el último segmento de la copia 1 de I3 acaba en 610, diez bases
/// más allá de la secuencia, y una deleción de segmento 600..610 quita lo que
/// sobra, de modo que ningún fragmento vigente sale de la referencia y la
/// materialización no lo ve. El recorrido acaba **después** de la longitud y se
/// rechaza ahí; distingue una regla que solo mirara si se queda corto.
#[test]
fn a_walk_past_the_end_of_the_sequence_is_rejected() {
    let mut all = document(&lines());
    let copy = copy_of(&mut all, "I3", 1);
    let mut mods = copy.segments[0].modifications.clone();
    assert_eq!(mods.len(), 2, "cambios heredados en 49 y 499");
    mods.push(TruthModification::Deletion {
        event_id: "DEL_OVERHANG".to_owned(),
        start: 600,
        end: 610,
        provenance: typed_provenance("del-1"),
    });
    copy.segments = vec![typed_segment(0, 610, mods)];
    assert_eq!(judge(&all).err(), Some(not_collinear("I3", 1, 610, 600)));
}
