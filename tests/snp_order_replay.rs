//! DT-018 (3 de octubre de 2026): el juez de la verdad reproduce la historia
//! de una copia aunque sus cambios de base entren fuera de orden de posición.
//!
//! La verdad publica los cambios de base de cada tramo **ordenados por
//! posición**. El juez los iba añadiendo en el orden del historial y comparaba
//! las listas sin ordenar, de modo que dos cambios de base directos declarados
//! en orden descendente, o uno directo seguido de uno de la matriz en una
//! posición posterior, abortaban la ejecución con un error interno.
//!
//! Dos grupos de pruebas:
//!
//! - **de punta a punta**, con los dos escenarios del revisor
//!   (`tests/fixtures/snp-order/`): la ejecución publica, y la verdad
//!   publicada lleva los cambios en orden de posición;
//! - **del juez**, con una verdad escrita a mano: acepta el historial
//!   descendente y **sigue rechazando** un cambio de más, uno de menos, uno con
//!   otra base, uno en otra posición y uno repetido. Estas pruebas son las que
//!   impiden resolver el defecto relajando la comparación.

use std::fs;
use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::truth_schema::TruthContractError;
use genome_population_simulator::truth_validation::{
    TruthValidationError, reference_content_sha256, validate_and_reconstruct_truth_jsonl,
};
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);

// ---------------------------------------------------------------- de punta a punta

/// Publica el escenario y devuelve la verdad; aborta la prueba si la
/// ejecución falla.
fn publish(label: &str) -> Vec<Value> {
    let source = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/snp-order")
        .join(label);
    let root: PathBuf = std::env::temp_dir().join(format!(
        "gps-snp-order-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("root");
    for name in ["general.toml", "reference.fa", "snvs.tsv"] {
        fs::copy(source.join(name), root.join(name)).expect("fixture file");
    }
    let config = load_and_validate_config(&root.join("general.toml"))
        .unwrap_or_else(|error| panic!("{label}: {}", error.render_human()));
    let out = root.join("published");
    execute_validated_config(&config, &out)
        .unwrap_or_else(|error| panic!("{label}: la ejecución abortó: {}", error.render_human()));
    let truth = fs::read_to_string(out.join("outputs/truth.jsonl"))
        .expect("truth")
        .lines()
        .map(|line| serde_json::from_str(line).expect("JSON line"))
        .collect();
    fs::remove_dir_all(&root).ok();
    truth
}

/// Posiciones (base cero) de los cambios de base de la copia `phase` de
/// `ind-1`, en el orden en que la verdad los publica, y si ese orden es el de
/// posición dentro de cada tramo.
fn published_snp_positions(truth: &[Value], phase: u64) -> Vec<u64> {
    let record = truth
        .iter()
        .find(|record| {
            record["record_type"] == "haplotype_copy"
                && record["identity"]["individual_id"] == "ind-1"
                && record["phase"] == phase
        })
        .expect("copy ind-1");
    let mut positions = Vec::new();
    for piece in record["current_pieces"].as_array().expect("pieces") {
        let piece_positions = piece["snps"]
            .as_array()
            .map(|snps| {
                snps.iter()
                    .map(|snp| snp["position"].as_u64().expect("position"))
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        // Esta aserción no puede fallar por sí sola: la lectura estructural de
        // la verdad (`truth_schema.rs`) ya aborta la ejecución si un tramo
        // publica sus cambios fuera de orden. Se deja como control de lo que
        // exige el esquema, no como discriminante de DT-018.
        assert!(
            piece_positions.windows(2).all(|pair| pair[0] < pair[1]),
            "la verdad debe publicar los cambios de cada tramo por posición: {piece_positions:?}"
        );
        positions.extend(piece_positions);
    }
    positions
}

#[test]
fn two_direct_snps_declared_in_descending_order_publish() {
    // Base 180 y después 150, en la fase 0 de ind-1, y una sola fila de matriz
    // lejos de ambos (posición 20). El orden ascendente publicaba.
    let truth = publish("descending-direct");
    let positions = published_snp_positions(&truth, 0);
    assert!(
        positions.contains(&149) && positions.contains(&179),
        "los dos cambios directos deben estar en la copia: {positions:?}"
    );
}

#[test]
fn a_direct_snp_followed_by_a_later_matrix_snp_in_the_same_copy_publishes() {
    // La matriz trae las posiciones 20 y 200; el directo cae en la 180 de la
    // fase 0 de ind-1. El historial pone los de la matriz antes que el
    // directo, así que el de la 200 llega antes que el de la 180.
    let truth = publish("direct-then-matrix");
    let positions = published_snp_positions(&truth, 0);
    assert!(
        positions.contains(&179) && positions.contains(&199),
        "control: la copia debe llevar el directo (179) y el de la matriz (199): {positions:?}"
    );
}

// ---------------------------------------------------------------- del juez

fn reference() -> Reference {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "snp order test").expect("sequence ID"),
        FastaHeader::parse(b"chr1").expect("FASTA header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(b"AAAAAAAA").expect("DNA"),
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("REF_main", "snp order test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn provenance(variant: &str) -> Value {
    json!({"kind":"matrix","version":"1","source":{"kind":"matrix","source_id":"MATRIX_01"},
           "matrix":{"matrix_source_id":"MATRIX_01","population_id":"POP_01","variant_id":variant}})
}

fn snp(event: &str, position: u64, alternate: &str, variant: &str) -> Value {
    json!({"event_id":event,"position":position,"reference":"A","alternate":alternate,
           "provenance":provenance(variant)})
}

/// La modificación de segmento de un cambio de base: lleva `kind`, que las
/// piezas vigentes no llevan.
fn modification(snp: &Value) -> Value {
    let mut modification = snp.clone();
    modification["kind"] = json!("snp");
    modification
}

fn history(event: &str, position: u64, variant: &str) -> Value {
    json!({"kind":"snp","event_id":event,"position":position,"provenance":provenance(variant)})
}

/// Verdad de una copia con dos cambios de base, EVT_02 en la posición 5
/// (alternativa C) y EVT_01 en la 2 (alternativa G), **en el historial en ese
/// orden descendente**. Las piezas vigentes los llevan por posición.
fn two_snp_truth(mutate: impl FnOnce(&mut Vec<Value>)) -> String {
    let mut lines = include_str!("fixtures/truth-v4/valid.jsonl")
        .replacen(
            "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            &reference_content_sha256(&reference()).expect("huella"),
            1,
        )
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    let first = snp("EVT_01", 2, "G", "VAR_01");
    let second = snp("EVT_02", 5, "C", "VAR_02");
    lines[1]["segments"][0]["modifications"] = json!([modification(&first), modification(&second)]);
    lines[1]["history"] = json!([
        history("EVT_02", 5, "VAR_02"),
        history("EVT_01", 2, "VAR_01")
    ]);
    let source = lines[1]["current_pieces"][0]["source"].clone();
    lines[1]["current_pieces"] = json!([{"kind":"source","source":source,"snps":[first, second]}]);
    lines[0]["counts"]["genotypes"] = json!(2);
    lines[3]["locus"] = json!({"sequence_id":"chr1","start":2,"end":3});
    let mut second_case = lines[3].clone();
    second_case["case_id"] = json!("CASE_02");
    second_case["locus"] = json!({"sequence_id":"chr1","start":5,"end":6});
    second_case["alleles"] = json!([{"kind":"reference"},{"kind":"snp","base":"C"}]);
    lines.push(second_case);
    mutate(&mut lines);
    lines
        .iter()
        .map(|line| serde_json::to_string(line).expect("serialize"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

/// Tres cambios de base con historial no monótono (4, 1, 6) y piezas
/// publicadas [1, 4, 6]. Con dos cambios, insertar al principio o junto al
/// último da el mismo resultado por casualidad; con tres, no.
fn three_snp_truth() -> String {
    let mut lines = include_str!("fixtures/truth-v4/valid.jsonl")
        .replacen(
            "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            &reference_content_sha256(&reference()).expect("huella"),
            1,
        )
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("fixture JSON"))
        .collect::<Vec<_>>();
    let a = snp("EVT_01", 1, "G", "VAR_01");
    let b = snp("EVT_02", 4, "C", "VAR_02");
    let c = snp("EVT_03", 6, "T", "VAR_03");
    lines[1]["segments"][0]["modifications"] =
        json!([modification(&a), modification(&b), modification(&c)]);
    lines[1]["history"] = json!([
        history("EVT_02", 4, "VAR_02"),
        history("EVT_01", 1, "VAR_01"),
        history("EVT_03", 6, "VAR_03")
    ]);
    let source = lines[1]["current_pieces"][0]["source"].clone();
    lines[1]["current_pieces"] = json!([{"kind":"source","source":source,"snps":[a, b, c]}]);
    lines[0]["counts"]["genotypes"] = json!(3);
    lines[3]["locus"] = json!({"sequence_id":"chr1","start":1,"end":2});
    for (id, pos, base) in [("CASE_02", 4, "C"), ("CASE_03", 6, "T")] {
        let mut case = lines[3].clone();
        case["case_id"] = json!(id);
        case["locus"] = json!({"sequence_id":"chr1","start":pos,"end":pos + 1});
        case["alleles"] = json!([{"kind":"reference"},{"kind":"snp","base":base}]);
        lines.push(case);
    }
    lines
        .iter()
        .map(|line| serde_json::to_string(line).expect("serialize"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

fn judge(truth: String) -> Result<Vec<Vec<u8>>, TruthValidationError> {
    validate_and_reconstruct_truth_jsonl(&mut Cursor::new(truth), &reference())
        .map(|reconstructed| reconstructed.copies.into_values().collect())
}

fn assert_does_not_replay(label: &str, truth: String) {
    match judge(truth) {
        Err(TruthValidationError::HistoryDoesNotReplay { .. }) => {}
        other => panic!("{label}: el juez debía rechazar con HistoryDoesNotReplay, dio {other:?}"),
    }
}

#[test]
fn judge_accepts_a_descending_history_with_position_ordered_pieces() {
    let copies = judge(two_snp_truth(|_| {})).expect("descending history replays");
    assert_eq!(copies[0], b"AAGAACAA".to_vec());
}

#[test]
fn judge_still_rejects_an_extra_snp() {
    assert_does_not_replay(
        "SNP de más",
        two_snp_truth(|lines| {
            lines[1]["current_pieces"][0]["snps"]
                .as_array_mut()
                .expect("snps")
                .push(snp("EVT_03", 6, "T", "VAR_03"));
        }),
    );
}

#[test]
fn judge_still_rejects_a_missing_snp() {
    assert_does_not_replay(
        "SNP de menos",
        two_snp_truth(|lines| {
            lines[1]["current_pieces"][0]["snps"]
                .as_array_mut()
                .expect("snps")
                .pop();
        }),
    );
}

#[test]
fn judge_still_rejects_a_snp_with_another_base() {
    assert_does_not_replay(
        "SNP con otra base",
        two_snp_truth(|lines| {
            lines[1]["current_pieces"][0]["snps"][1]["alternate"] = json!("T");
        }),
    );
}

#[test]
fn judge_still_rejects_a_snp_at_another_position() {
    assert_does_not_replay(
        "SNP en otra posición",
        two_snp_truth(|lines| {
            lines[1]["current_pieces"][0]["snps"][1]["position"] = json!(6);
        }),
    );
}

#[test]
fn judge_accepts_a_non_monotone_history_of_three_snps() {
    let copies = judge(three_snp_truth()).expect("non-monotone history replays");
    assert_eq!(copies[0], b"AGAACATA".to_vec());
}

#[test]
fn two_snps_at_the_same_position_are_rejected_by_the_structural_reader() {
    // El rechazo es de la lectura estructural (`ModificationsNotOrdered`), no
    // de la guarda de posición repetida del reproductor, que ningún documento
    // válido alcanza.
    let truth = two_snp_truth(|lines| {
        let mut repeated = snp("EVT_02", 2, "C", "VAR_02");
        repeated["kind"] = json!("snp");
        lines[1]["segments"][0]["modifications"][1] = repeated;
    });
    match judge(truth) {
        Err(TruthValidationError::Contract(TruthContractError::ModificationsNotOrdered {
            ..
        })) => {}
        other => panic!("debía rechazar con ModificationsNotOrdered, dio {other:?}"),
    }
}
