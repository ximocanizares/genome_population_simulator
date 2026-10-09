//! Cada rechazo del contrato estructural de la verdad, con su motivo exacto.
//!
//! DT-001, tanda 5. El contrato (`TruthDocument::validate_structure`, en
//! `src/truth_schema.rs`) tenía 69 comprobaciones de las que 44 no las
//! disparaba ninguna prueba. Aquí cada caso parte de una verdad **publicada de
//! verdad** por el programa (la del escenario de `truth_guard_scenario`, en su
//! versión 5, con alelos de repetición), la altera **en un solo punto** y exige
//! el error exacto, con su campo y su valor. Cada tabla comprueba además que
//! el documento intacto se acepta y que la alteración alteró algo: sin lo
//! primero, un contrato que rechazara todo satisfaría la tabla entera.
//!
//! El contrato trabaja sobre el documento **en memoria**, no sobre bytes, y es
//! lo que recibe la publicación antes de escribir la verdad: por eso los
//! casos alteran el `TruthDocument` leído y no el texto.
//!
//! # Las comprobaciones inalcanzables
//!
//! Las que ningún documento puede alcanzar están anotadas **junto a la guarda**
//! en `src/truth_schema.rs`, con el porqué. `unreachable_guards_are_shadowed_by_earlier_checks`
//! ejecuta el intento y lo dice: esa prueba **no discrimina** esas guardas.

use genome_population_simulator::truth_schema::{
    CopyIdentity, GenotypeCaseKind, HaplotypeCopyRecord, ModificationProvenance,
    ModificationProvenanceKind, NamedDigest, SegmentSource, TruthAllele, TruthContractError,
    TruthCurrentPiece, TruthCurrentSnp, TruthDocument, TruthHistoryEvent, TruthModification,
    TruthPieceLink, TruthRecord, TruthRepeatBlock, TruthSegment, TruthSourceKind, TruthSourceRef,
    TruthStructuralEvent,
};

mod truth_guard_scenario;
use truth_guard_scenario::*;

fn base() -> TruthDocument {
    document(&lines_v5())
}

type Edit = Box<dyn Fn(&mut TruthDocument)>;

/// Una alteración. La función ayuda a la inferencia del tipo del cierre.
fn edit(change: impl Fn(&mut TruthDocument) + 'static) -> Edit {
    Box::new(change)
}

/// El documento intacto se acepta; cada alteración da **exactamente** el error
/// esperado, y de verdad alteró el documento.
#[track_caller]
fn assert_table(table: Vec<(&str, Edit, TruthContractError)>) {
    base()
        .validate_structure()
        .expect("el documento intacto se acepta");
    for (label, edit, expected) in table {
        let mut altered = base();
        edit(&mut altered);
        assert_ne!(altered, base(), "{label}: la alteración alteró");
        assert_eq!(
            altered.validate_structure().err(),
            Some(expected),
            "{label}"
        );
    }
}

/// Cada alteración se acepta.
#[track_caller]
fn assert_all_accepted(table: Vec<(&str, Edit)>) {
    for (label, edit) in table {
        let mut altered = base();
        edit(&mut altered);
        assert_ne!(altered, base(), "{label}: la alteración alteró");
        altered
            .validate_structure()
            .unwrap_or_else(|error| panic!("{label}: se esperaba aceptarlo y se obtuvo {error:?}"));
    }
}

fn lexical(field: &'static str, value: &str) -> TruthContractError {
    TruthContractError::InvalidLexicalField {
        field,
        value: value.to_owned(),
    }
}

const BAD_ID: &str = "bad id";

fn copy<'a>(
    document: &'a mut TruthDocument,
    individual: &str,
    index: u8,
) -> &'a mut HaplotypeCopyRecord {
    copy_of(document, individual, index)
}

fn case<'a>(
    document: &'a mut TruthDocument,
    individual: &str,
    kind: GenotypeCaseKind,
    start: u64,
) -> &'a mut genome_population_simulator::truth_schema::GenotypeRecord {
    genotype_of(document, individual, kind, start)
}

/// La copia 0 de I1 de la versión 5: tiene cambios de base, una inserción y un
/// tracto.
fn rich(document: &mut TruthDocument) -> &mut HaplotypeCopyRecord {
    copy(document, "I1", 0)
}

fn history_of(document: &mut TruthDocument, individual: &str, index: u8, kind: &str) -> usize {
    copy(document, individual, index)
        .history
        .iter()
        .position(|event| {
            matches!(
                (event, kind),
                (TruthHistoryEvent::Snp { .. }, "snp")
                    | (TruthHistoryEvent::Insertion { .. }, "insertion")
                    | (TruthHistoryEvent::Deletion { .. }, "deletion")
                    | (TruthHistoryEvent::RepeatAllele { .. }, "repeat")
            )
        })
        .expect("un evento de esa clase")
}

fn provenance_of_history(event: &mut TruthHistoryEvent) -> &mut ModificationProvenance {
    match event {
        TruthHistoryEvent::Snp { provenance, .. }
        | TruthHistoryEvent::Insertion { provenance, .. }
        | TruthHistoryEvent::Deletion { provenance, .. }
        | TruthHistoryEvent::RepeatAllele { provenance, .. } => provenance,
    }
}

fn event_id_of_history(event: &mut TruthHistoryEvent) -> &mut String {
    match event {
        TruthHistoryEvent::Snp { event_id, .. }
        | TruthHistoryEvent::Insertion { event_id, .. }
        | TruthHistoryEvent::Deletion { event_id, .. }
        | TruthHistoryEvent::RepeatAllele { event_id, .. } => event_id,
    }
}

fn first_source_piece(
    copy: &mut HaplotypeCopyRecord,
) -> (&mut SegmentSource, &mut Vec<TruthCurrentSnp>) {
    copy.current_pieces
        .iter_mut()
        .find_map(|piece| match piece {
            TruthCurrentPiece::Source { source, snps } => Some((source, snps)),
            _ => None,
        })
        .expect("un fragmento de fuente")
}

fn source_piece_with_snps(
    copy: &mut HaplotypeCopyRecord,
) -> (&mut SegmentSource, &mut Vec<TruthCurrentSnp>) {
    copy.current_pieces
        .iter_mut()
        .find_map(|piece| match piece {
            TruthCurrentPiece::Source { source, snps } if !snps.is_empty() => Some((source, snps)),
            _ => None,
        })
        .expect("un fragmento con cambios de base")
}

fn inserted_piece_of(copy: &mut HaplotypeCopyRecord) -> (&mut String, &mut ModificationProvenance) {
    copy.current_pieces
        .iter_mut()
        .find_map(|piece| match piece {
            TruthCurrentPiece::Inserted { bases, provenance } => Some((bases, provenance)),
            _ => None,
        })
        .expect("un fragmento insertado")
}

type RepeatPiece<'a> = (
    &'a mut String,
    &'a mut SegmentSource,
    &'a mut String,
    &'a mut Vec<TruthRepeatBlock>,
    &'a mut ModificationProvenance,
);

fn repeat_piece_of(copy: &mut HaplotypeCopyRecord) -> RepeatPiece<'_> {
    copy.current_pieces
        .iter_mut()
        .find_map(|piece| match piece {
            TruthCurrentPiece::RepeatAllele {
                locus_id,
                source,
                bases,
                blocks,
                provenance,
            } => Some((locus_id, source, bases, blocks, provenance)),
            _ => None,
        })
        .expect("un tracto vigente")
}

fn repeat_event_of(
    copy: &mut HaplotypeCopyRecord,
) -> (
    &mut String,
    &mut u64,
    &mut u64,
    &mut String,
    &mut Vec<TruthRepeatBlock>,
) {
    copy.history
        .iter_mut()
        .find_map(|event| match event {
            TruthHistoryEvent::RepeatAllele {
                locus_id,
                start,
                end,
                bases,
                blocks,
                ..
            } => Some((locus_id, start, end, bases, blocks)),
            _ => None,
        })
        .expect("un evento de tracto")
}

fn snp_modification_of(
    copy: &mut HaplotypeCopyRecord,
) -> (
    &mut String,
    &mut String,
    &mut String,
    &mut ModificationProvenance,
) {
    copy.segments
        .iter_mut()
        .flat_map(|segment| segment.modifications.iter_mut())
        .find_map(|modification| match modification {
            TruthModification::Snp {
                event_id,
                reference,
                alternate,
                provenance,
                ..
            } => Some((event_id, reference, alternate, provenance)),
            _ => None,
        })
        .expect("un cambio de base de segmento")
}

fn repeat_allele_of(
    genotype: &mut genome_population_simulator::truth_schema::GenotypeRecord,
) -> (&mut String, &mut String, &mut Vec<TruthRepeatBlock>) {
    genotype
        .alleles
        .iter_mut()
        .find_map(|allele| match allele {
            TruthAllele::Repeat {
                locus_id,
                bases,
                blocks,
            } => Some((locus_id, bases, blocks)),
            _ => None,
        })
        .expect("un alelo de repetición")
}

fn first_repeat_case(
    document: &mut TruthDocument,
) -> &mut genome_population_simulator::truth_schema::GenotypeRecord {
    document
        .records
        .iter_mut()
        .find_map(|record| match record {
            TruthRecord::Genotype(genotype) if genotype.case_kind == GenotypeCaseKind::Repeat => {
                Some(genotype)
            }
            _ => None,
        })
        .expect("un caso de repetición")
}

// ---------------------------------------------------------------------------
// Las condiciones de cada ayudante léxico, con un solo llamante
// ---------------------------------------------------------------------------

/// **`valid_id`**: hasta 64 caracteres, no vacío, que empiece por letra o
/// dígito y que el resto sean letras, dígitos, `.`, `_` o `-`. Cada condición
/// con su rechazo, y el lado aceptado de cada una (incluido el límite de 64).
#[test]
fn an_identifier_is_at_most_sixty_four_alphanumerics_dots_underscores_and_dashes() {
    let rejected = [
        ("sesenta y cinco caracteres", "a".repeat(65)),
        ("vacío", String::new()),
        ("empieza por guion bajo", "_a".to_owned()),
        ("empieza por punto", ".a".to_owned()),
        ("empieza por guion", "-a".to_owned()),
        ("contiene un espacio", "a b".to_owned()),
        ("contiene una barra", "a/b".to_owned()),
        ("contiene un carácter no ASCII", "aé".to_owned()),
    ];
    let table = rejected
        .iter()
        .map(|(label, value)| {
            let value = value.clone();
            let edit: Edit = Box::new(move |document| document.header.execution_id = value.clone());
            (
                *label,
                edit,
                lexical("execution_id", &rejected_value(label)),
            )
        })
        .collect();
    assert_table(table);
    assert_all_accepted(
        [
            ("sesenta y cuatro caracteres", "a".repeat(64)),
            ("empieza por dígito", "0a".to_owned()),
            ("con punto, guion bajo y guion", "a.b_c-d".to_owned()),
            ("un solo carácter", "z".to_owned()),
        ]
        .into_iter()
        .map(|(label, value)| {
            let edit: Edit = Box::new(move |document| document.header.execution_id = value.clone());
            (label, edit)
        })
        .collect(),
    );
}

fn rejected_value(label: &str) -> String {
    match label {
        "sesenta y cinco caracteres" => "a".repeat(65),
        "vacío" => String::new(),
        "empieza por guion bajo" => "_a".to_owned(),
        "empieza por punto" => ".a".to_owned(),
        "empieza por guion" => "-a".to_owned(),
        "contiene un espacio" => "a b".to_owned(),
        "contiene una barra" => "a/b".to_owned(),
        "contiene un carácter no ASCII" => "aé".to_owned(),
        other => panic!("etiqueta desconocida {other}"),
    }
}

/// **`valid_sha`**: sesenta y cuatro caracteres hexadecimales **en minúscula**.
#[test]
fn a_digest_is_sixty_four_lowercase_hexadecimal_characters() {
    let rejected = [
        ("sesenta y tres", "a".repeat(63)),
        ("sesenta y cinco", "a".repeat(65)),
        ("con mayúscula", format!("A{}", "a".repeat(63))),
        (
            "con una letra que no es hexadecimal",
            format!("g{}", "a".repeat(63)),
        ),
        ("con un signo", format!("-{}", "a".repeat(63))),
    ];
    let mut table = Vec::new();
    for (label, value) in rejected {
        let shown = value.clone();
        let edit: Edit = Box::new(move |document| {
            document.header.reference_content_sha256 = value.clone();
        });
        table.push((label, edit, lexical("sha256", &shown)));
    }
    assert_table(table);
    // El mismo rechazo por la huella de una fuente de la cabecera.
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "huella de una fuente con mayúscula",
            Box::new(|document| {
                document.header.source_digests[0].sha256 = "F".repeat(64);
            }),
            lexical("sha256", &"F".repeat(64)),
        ),
        (
            "huella de un grupo demasiado corta",
            Box::new(|document| {
                document.header.group_digests[0].sha256 = "a".repeat(63);
            }),
            lexical("sha256", &"a".repeat(63)),
        ),
    ];
    assert_table(table);
    // Lado aceptado: la huella con los dos extremos del alfabeto hexadecimal.
    assert_all_accepted(vec![(
        "huella 0123456789abcdef",
        edit(|document| {
            document.header.reference_content_sha256 = "0123456789abcdef".repeat(4);
        }),
    )]);
}

/// **`valid_base`**: una de `A`, `C`, `G`, `T`, en mayúscula.
#[test]
fn a_single_base_is_one_of_the_four_uppercase_letters() {
    for (label, value) in [
        ("N", "N"),
        ("minúscula", "a"),
        ("dos bases", "AC"),
        ("vacía", ""),
    ] {
        let owned = value.to_owned();
        let edit: Edit = Box::new(move |document| {
            snp_modification_of(rich(document)).1.clone_from(&owned);
        });
        assert_table(vec![(
            label,
            edit,
            lexical("historical_snp_reference", value),
        )]);
    }
    // Lado aceptado: cada una de las cuatro bases, como referencia del cambio.
    for (reference, alternate) in [("A", "C"), ("C", "G"), ("G", "T"), ("T", "A")] {
        let edit: Edit = Box::new(move |document| {
            let (_, base, other, _) = snp_modification_of(rich(document));
            *base = reference.to_owned();
            *other = alternate.to_owned();
        });
        assert_all_accepted(vec![(reference, edit)]);
    }
}

/// **`valid_bases`**: no vacía, de `A`, `C`, `G`, `T`, `N` en mayúscula.
#[test]
fn a_run_of_bases_is_non_empty_and_made_of_acgtn() {
    let table = ["", "ACX", "acgt", "A C"]
        .into_iter()
        .map(|value| {
            let owned = value.to_owned();
            let edit: Edit = Box::new(move |document| {
                inserted_piece_of(rich(document)).0.clone_from(&owned);
            });
            ("bases", edit, lexical("current_inserted_bases", value))
        })
        .collect();
    assert_table(table);
    assert_all_accepted(vec![(
        "con N",
        edit(|document| {
            *inserted_piece_of(rich(document)).0 = "ACGTN".to_owned();
        }),
    )]);
}

// ---------------------------------------------------------------------------
// Cada llamante de los ayudantes léxicos
// ---------------------------------------------------------------------------

fn bad_provenance(provenance: &mut ModificationProvenance) {
    provenance.source.source_id = BAD_ID.to_owned();
}

fn plain_provenance() -> ModificationProvenance {
    ModificationProvenance {
        kind: ModificationProvenanceKind::Explicit,
        version: "1".to_owned(),
        source: TruthSourceRef {
            kind: TruthSourceKind::Explicit,
            source_id: "ins-1".to_owned(),
        },
        matrix: None,
    }
}

fn add_insertion_modification(copy: &mut HaplotypeCopyRecord, event_id: &str, bases: &str) {
    copy.segments[0]
        .modifications
        .push(TruthModification::Insertion {
            event_id: event_id.to_owned(),
            point: 599,
            bases: bases.to_owned(),
            provenance: plain_provenance(),
        });
}

fn add_deletion_modification(copy: &mut HaplotypeCopyRecord, event_id: &str) {
    copy.segments[0]
        .modifications
        .push(TruthModification::Deletion {
            event_id: event_id.to_owned(),
            start: 590,
            end: 595,
            provenance: plain_provenance(),
        });
}

/// Un identificador que no cumple la sintaxis, en cada sitio donde el
/// documento lleva uno: cada llamante de `valid_id`, con su nombre de campo y
/// el valor que lo incumple. Con un solo llamante probado, borrar los demás
/// dejaría pasar la batería.
#[test]
fn every_identifier_the_document_carries_is_checked() {
    let id = BAD_ID;
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "execution_id",
            edit(|d| d.header.execution_id = BAD_ID.into()),
            lexical("execution_id", id),
        ),
        (
            "reference_id de la cabecera",
            edit(|d| d.header.reference_id = BAD_ID.into()),
            lexical("reference_id", id),
        ),
        (
            "fuente de una huella",
            edit(|d| d.header.source_digests[0].source.source_id = BAD_ID.into()),
            lexical("source_id", id),
        ),
        (
            "fuente de una huella de grupo",
            edit(|d| d.header.group_digests[0].source.source_id = BAD_ID.into()),
            lexical("source_id", id),
        ),
        (
            "reference_id de la identidad",
            edit(|d| rich(d).identity.reference_id = BAD_ID.into()),
            lexical("copy_identity", id),
        ),
        (
            "individual_id de la identidad",
            edit(|d| rich(d).identity.individual_id = BAD_ID.into()),
            lexical("copy_identity", id),
        ),
        (
            "sequence_id de la identidad",
            edit(|d| rich(d).identity.sequence_id = BAD_ID.into()),
            lexical("copy_identity", id),
        ),
        (
            "fuente de la procedencia de la copia",
            edit(|d| rich(d).provenance.sources[0].source_id = BAD_ID.into()),
            lexical("source_id", id),
        ),
        (
            "fuente de un segmento",
            edit(|d| rich(d).segments[0].source.source.source_id = BAD_ID.into()),
            lexical("source_id", id),
        ),
        (
            "referencia de un segmento",
            edit(|d| rich(d).segments[0].source.reference_id = BAD_ID.into()),
            lexical("segment_reference_id", id),
        ),
        (
            "secuencia de un segmento",
            edit(|d| rich(d).segments[0].source.sequence_id = BAD_ID.into()),
            lexical("segment_sequence_id", id),
        ),
        (
            "evento de un cambio de segmento",
            edit(|d| *snp_modification_of(rich(d)).0 = BAD_ID.into()),
            lexical("event_id", id),
        ),
        (
            "procedencia de un cambio de segmento",
            edit(|d| bad_provenance(snp_modification_of(rich(d)).3)),
            lexical("source_id", id),
        ),
        (
            "evento de una inserción de segmento",
            edit(|d| add_insertion_modification(rich(d), BAD_ID, "ACG")),
            lexical("event_id", id),
        ),
        (
            "evento de una deleción de segmento",
            edit(|d| add_deletion_modification(rich(d), BAD_ID)),
            lexical("event_id", id),
        ),
        (
            "fuente de un fragmento vigente",
            edit(|d| first_source_piece(rich(d)).0.source.source_id = BAD_ID.into()),
            lexical("source_id", id),
        ),
        (
            "referencia de un fragmento vigente",
            edit(|d| first_source_piece(rich(d)).0.reference_id = BAD_ID.into()),
            lexical("segment_reference_id", id),
        ),
        (
            "secuencia de un fragmento vigente",
            edit(|d| first_source_piece(rich(d)).0.sequence_id = BAD_ID.into()),
            lexical("segment_sequence_id", id),
        ),
        (
            "evento de un cambio vigente",
            edit(|d| source_piece_with_snps(rich(d)).1[0].event_id = BAD_ID.into()),
            lexical("event_id", id),
        ),
        (
            "procedencia de un cambio vigente",
            edit(|d| bad_provenance(&mut source_piece_with_snps(rich(d)).1[0].provenance)),
            lexical("source_id", id),
        ),
        (
            "procedencia de unas bases insertadas",
            edit(|d| bad_provenance(inserted_piece_of(rich(d)).1)),
            lexical("source_id", id),
        ),
        (
            "locus de un tracto vigente",
            edit(|d| *repeat_piece_of(rich(d)).0 = BAD_ID.into()),
            lexical("repeat_locus_id", id),
        ),
        (
            "fuente de un tracto vigente",
            edit(|d| repeat_piece_of(rich(d)).1.source.source_id = BAD_ID.into()),
            lexical("source_id", id),
        ),
        (
            "referencia de un tracto vigente",
            edit(|d| repeat_piece_of(rich(d)).1.reference_id = BAD_ID.into()),
            lexical("segment_reference_id", id),
        ),
        (
            "secuencia de un tracto vigente",
            edit(|d| repeat_piece_of(rich(d)).1.sequence_id = BAD_ID.into()),
            lexical("segment_sequence_id", id),
        ),
        (
            "procedencia de un tracto vigente",
            edit(|d| bad_provenance(repeat_piece_of(rich(d)).4)),
            lexical("source_id", id),
        ),
        (
            "evento de un cambio del historial",
            edit(|d| {
                let at = history_of(d, "I1", 0, "snp");
                *event_id_of_history(&mut rich(d).history[at]) = BAD_ID.into();
            }),
            lexical("event_id", id),
        ),
        (
            "procedencia de un cambio del historial",
            edit(|d| {
                let at = history_of(d, "I1", 0, "snp");
                bad_provenance(provenance_of_history(&mut rich(d).history[at]));
            }),
            lexical("source_id", id),
        ),
        (
            "evento de una inserción del historial",
            edit(|d| {
                let at = history_of(d, "I1", 0, "insertion");
                *event_id_of_history(&mut rich(d).history[at]) = BAD_ID.into();
            }),
            lexical("event_id", id),
        ),
        (
            "procedencia de una inserción del historial",
            edit(|d| {
                let at = history_of(d, "I1", 0, "insertion");
                bad_provenance(provenance_of_history(&mut rich(d).history[at]));
            }),
            lexical("source_id", id),
        ),
        (
            "evento de una deleción del historial",
            edit(|d| {
                let at = history_of(d, "I1", 1, "deletion");
                *event_id_of_history(&mut copy(d, "I1", 1).history[at]) = BAD_ID.into();
            }),
            lexical("event_id", id),
        ),
        (
            "procedencia de una deleción del historial",
            edit(|d| {
                let at = history_of(d, "I1", 1, "deletion");
                bad_provenance(provenance_of_history(&mut copy(d, "I1", 1).history[at]));
            }),
            lexical("source_id", id),
        ),
        (
            "evento de un tracto del historial",
            edit(|d| {
                let at = history_of(d, "I1", 0, "repeat");
                *event_id_of_history(&mut rich(d).history[at]) = BAD_ID.into();
            }),
            lexical("event_id", id),
        ),
        (
            "locus de un tracto del historial",
            edit(|d| *repeat_event_of(rich(d)).0 = BAD_ID.into()),
            lexical("repeat_locus_id", id),
        ),
        (
            "procedencia de un tracto del historial",
            edit(|d| {
                let at = history_of(d, "I1", 0, "repeat");
                bad_provenance(provenance_of_history(&mut rich(d).history[at]));
            }),
            lexical("source_id", id),
        ),
        (
            "case_id de un caso",
            edit(|d| case(d, "I1", GenotypeCaseKind::Snp, 49).case_id = BAD_ID.into()),
            lexical("genotype", id),
        ),
        (
            "reference_id de un caso",
            edit(|d| case(d, "I1", GenotypeCaseKind::Snp, 49).reference_id = BAD_ID.into()),
            lexical("genotype", id),
        ),
        (
            "individual_id de un caso",
            edit(|d| case(d, "I1", GenotypeCaseKind::Snp, 49).individual_id = BAD_ID.into()),
            lexical("genotype", id),
        ),
        (
            "secuencia del locus de un caso",
            edit(|d| case(d, "I1", GenotypeCaseKind::Snp, 49).locus.sequence_id = BAD_ID.into()),
            lexical("genotype", id),
        ),
        (
            "locus de un alelo de repetición",
            edit(|d| *repeat_allele_of(first_repeat_case(d)).0 = BAD_ID.into()),
            lexical("repeat_locus_id", id),
        ),
        (
            "evento de un alelo de inserción",
            edit(
                |d| match &mut case(d, "I1", GenotypeCaseKind::Indel, 350).alleles[1] {
                    TruthAllele::Insertion { event_id, .. } => *event_id = BAD_ID.into(),
                    _ => panic!("inserción"),
                },
            ),
            lexical("indel_event_id", id),
        ),
        (
            "evento de un alelo de deleción",
            edit(
                |d| match &mut case(d, "I1", GenotypeCaseKind::Indel, 420).alleles[1] {
                    TruthAllele::Deletion { event_id, .. } => *event_id = BAD_ID.into(),
                    _ => panic!("deleción"),
                },
            ),
            lexical("indel_event_id", id),
        ),
        (
            "reference_id de una copia portadora",
            edit(|d| {
                case(d, "I1", GenotypeCaseKind::Snp, 49).carriers[0]
                    .copy
                    .reference_id = BAD_ID.into()
            }),
            lexical("carrier_copy", id),
        ),
        (
            "individual_id de una copia portadora",
            edit(|d| {
                case(d, "I1", GenotypeCaseKind::Snp, 49).carriers[0]
                    .copy
                    .individual_id = BAD_ID.into()
            }),
            lexical("carrier_copy", id),
        ),
        (
            "sequence_id de una copia portadora",
            edit(|d| {
                case(d, "I1", GenotypeCaseKind::Snp, 49).carriers[0]
                    .copy
                    .sequence_id = BAD_ID.into()
            }),
            lexical("carrier_copy", id),
        ),
    ];
    assert_table(table);
}

fn matrix_provenance_of(
    copy: &mut HaplotypeCopyRecord,
) -> &mut genome_population_simulator::truth_schema::MatrixProvenance {
    copy.segments
        .iter_mut()
        .flat_map(|segment| segment.modifications.iter_mut())
        .find_map(|modification| match modification {
            TruthModification::Snp { provenance, .. } => provenance.matrix.as_mut(),
            _ => None,
        })
        .expect("un cambio de base de la matriz")
}

/// Los identificadores de la matriz de una procedencia: cada uno con su campo.
#[test]
fn the_identifiers_of_a_matrix_provenance_are_checked() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "matrix_source_id",
            edit(|d| matrix_provenance_of(rich(d)).matrix_source_id = BAD_ID.into()),
            lexical("matrix_source_id", BAD_ID),
        ),
        (
            "population_id",
            edit(|d| matrix_provenance_of(rich(d)).population_id = BAD_ID.into()),
            lexical("population_id", BAD_ID),
        ),
        (
            "variant_id",
            edit(|d| matrix_provenance_of(rich(d)).variant_id = BAD_ID.into()),
            lexical("variant_id", BAD_ID),
        ),
    ];
    assert_table(table);
}

/// Las bases de cada sitio donde el documento lleva una cadena de bases, con
/// una base que no es de `ACGTN`.
#[test]
fn every_run_of_bases_the_document_carries_is_checked() {
    let bad = "ACX";
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "bases de una inserción de segmento",
            edit(|d| add_insertion_modification(rich(d), "EVT_ADDED", "ACX")),
            lexical("historical_insertion_bases", bad),
        ),
        (
            "bases insertadas vigentes",
            edit(|d| *inserted_piece_of(rich(d)).0 = bad.into()),
            lexical("current_inserted_bases", bad),
        ),
        (
            "bases de un tracto vigente",
            edit(|d| *repeat_piece_of(rich(d)).2 = bad.into()),
            lexical("repeat_allele_bases", bad),
        ),
        (
            "bases de una inserción del historial",
            edit(|d| {
                let at = history_of(d, "I1", 0, "insertion");
                if let TruthHistoryEvent::Insertion { bases, .. } = &mut rich(d).history[at] {
                    *bases = bad.into();
                }
            }),
            lexical("history_insertion_bases", bad),
        ),
        (
            "bases de un tracto del historial",
            edit(|d| *repeat_event_of(rich(d)).3 = bad.into()),
            lexical("repeat_allele_bases", bad),
        ),
        (
            "bases de un alelo de repetición",
            edit(|d| *repeat_allele_of(first_repeat_case(d)).1 = bad.into()),
            lexical("repeat_allele_bases", bad),
        ),
        (
            "bases de un alelo de inserción",
            edit(|d| {
                if let TruthAllele::Insertion { bases, .. } =
                    &mut case(d, "I1", GenotypeCaseKind::Indel, 350).alleles[1]
                {
                    *bases = bad.into();
                }
            }),
            lexical("insertion_allele_bases", bad),
        ),
        (
            "motivo de un bloque",
            edit(|d| repeat_piece_of(rich(d)).3[0].motif = bad.into()),
            lexical("repeat_block_motif", bad),
        ),
    ];
    assert_table(table);
}

/// Las bases sueltas de cada sitio donde el documento lleva una: una base que
/// no es `A`, `C`, `G` ni `T`.
#[test]
fn every_single_base_the_document_carries_is_checked() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "alterna de un cambio de segmento",
            edit(|d| *snp_modification_of(rich(d)).2 = "N".into()),
            lexical("historical_snp_alternate", "N"),
        ),
        (
            "referencia de un cambio vigente",
            edit(|d| source_piece_with_snps(rich(d)).1[0].reference = "N".into()),
            lexical("current_snp_reference", "N"),
        ),
        (
            "alterna de un cambio vigente",
            edit(|d| source_piece_with_snps(rich(d)).1[0].alternate = "N".into()),
            lexical("current_snp_alternate", "N"),
        ),
        (
            "base de un alelo SNP",
            edit(|d| {
                if let TruthAllele::Snp { base } =
                    &mut case(d, "I1", GenotypeCaseKind::Snp, 49).alleles[1]
                {
                    *base = "N".into();
                }
            }),
            lexical("snp_allele_base", "N"),
        ),
    ];
    assert_table(table);
}

/// El índice y la fase de una copia, y el índice de una copia portadora, no
/// pasan de nueve; la procedencia de una copia no está vacía.
#[test]
fn copy_indices_phases_and_origin_are_bounded_and_not_empty() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "fase 10",
            edit(|d| rich(d).phase = 10),
            lexical("copy_index_or_phase", "0"),
        ),
        (
            "índice 10",
            edit(|d| {
                let copy = rich(d);
                copy.identity.copy_index = 10;
                copy.phase = 10;
            }),
            lexical("copy_index_or_phase", "10"),
        ),
        (
            "índice de una copia portadora 10",
            edit(|d| {
                case(d, "I1", GenotypeCaseKind::Snp, 49).carriers[0]
                    .copy
                    .copy_index = 10
            }),
            lexical("carrier_copy_index", "10"),
        ),
        (
            "procedencia vacía",
            edit(|d| rich(d).provenance.origin = String::new()),
            lexical("copy_provenance_origin", ""),
        ),
    ];
    assert_table(table);
    // Lado aceptado: el índice y la fase 9, con la copia y sus portadores renombrados.
    assert_all_accepted(vec![(
        "índice 9",
        edit(|d| {
            let renamed = copy(d, "I3", 1);
            renamed.identity.copy_index = 9;
            renamed.phase = 9;
            for record in &mut d.records {
                if let TruthRecord::Genotype(genotype) = record {
                    for carrier in &mut genotype.carriers {
                        if carrier.copy.individual_id == "I3"
                            && carrier.copy.sequence_id == "chr1"
                            && carrier.copy.copy_index == 1
                        {
                            carrier.copy.copy_index = 9;
                        }
                    }
                }
            }
        }),
    )]);
}

// ---------------------------------------------------------------------------
// Los bloques de un tracto, en cada sitio donde se publican
// ---------------------------------------------------------------------------

type Blocks = fn(&mut TruthDocument) -> &mut Vec<TruthRepeatBlock>;

fn current_blocks(document: &mut TruthDocument) -> &mut Vec<TruthRepeatBlock> {
    repeat_piece_of(rich(document)).3
}
fn history_blocks(document: &mut TruthDocument) -> &mut Vec<TruthRepeatBlock> {
    repeat_event_of(rich(document)).4
}
fn allele_blocks(document: &mut TruthDocument) -> &mut Vec<TruthRepeatBlock> {
    repeat_allele_of(first_repeat_case(document)).2
}

/// Las bases del tracto que acompaña a cada lista de bloques.
fn tract_lengths() -> [u64; 3] {
    let mut document = base();
    let current = repeat_piece_of(rich(&mut document)).2.len();
    let history = repeat_event_of(rich(&mut document)).3.len();
    let allele = repeat_allele_of(first_repeat_case(&mut document)).1.len();
    [current, history, allele].map(|length| u64::try_from(length).expect("cabe"))
}

/// Cada defecto de una lista de bloques, en los tres sitios donde se publica:
/// vacía; un bloque sin bases; un motivo vacío; una suma que desborda; una
/// suma que no es la longitud del tracto.
#[test]
fn a_malformed_list_of_blocks_is_rejected_in_every_place() {
    let places: [(&str, Blocks); 3] = [
        ("fragmento vigente", current_blocks),
        ("historial", history_blocks),
        ("alelo del caso", allele_blocks),
    ];
    let lengths = tract_lengths();
    for ((place, blocks), length) in places.into_iter().zip(lengths) {
        let big = u64::MAX / 2 + 1;
        let table: Vec<(String, Edit, TruthContractError)> = vec![
            (
                format!("{place}: sin bloques"),
                edit(move |d| blocks(d).clear()),
                lexical("repeat_blocks", ""),
            ),
            (
                format!("{place}: un bloque sin bases"),
                edit(move |d| blocks(d)[0].bases = 0),
                lexical("repeat_block_bases", "0"),
            ),
            (
                format!("{place}: un motivo vacío"),
                edit(move |d| blocks(d)[0].motif.clear()),
                lexical("repeat_block_motif", ""),
            ),
            (
                format!("{place}: una suma que desborda"),
                edit(move |d| {
                    let list = blocks(d);
                    list.clear();
                    for _ in 0..2 {
                        list.push(TruthRepeatBlock {
                            motif: "CA".to_owned(),
                            bases: big,
                        });
                    }
                }),
                TruthContractError::ArithmeticOverflow {
                    field: "repeat_block_bases",
                },
            ),
            (
                format!("{place}: una suma de más"),
                edit(move |d| blocks(d)[0].bases += 2),
                lexical("repeat_blocks_total_bases", &(length + 2).to_string()),
            ),
            (
                format!("{place}: una suma de menos"),
                edit(move |d| blocks(d)[0].bases -= 2),
                lexical("repeat_blocks_total_bases", &(length - 2).to_string()),
            ),
        ];
        base()
            .validate_structure()
            .expect("el documento intacto se acepta");
        for (label, change, expected) in table {
            let mut altered = base();
            change(&mut altered);
            assert_ne!(altered, base(), "{label}: la alteración alteró");
            assert_eq!(
                altered.validate_structure().err(),
                Some(expected),
                "{label}"
            );
        }
    }
}

/// Lado aceptado: una lista de dos bloques cuya suma es la longitud del
/// tracto, y el mayor bloque que cabe sin desbordar.
#[test]
fn a_list_of_blocks_adding_up_to_the_tract_is_accepted_in_every_place() {
    let places: [Blocks; 3] = [current_blocks, history_blocks, allele_blocks];
    for blocks in places {
        let mut document = base();
        let list = blocks(&mut document);
        let total: u64 = list.iter().map(|block| block.bases).sum();
        list.clear();
        list.push(TruthRepeatBlock {
            motif: "CA".to_owned(),
            bases: total - 1,
        });
        list.push(TruthRepeatBlock {
            motif: "C".to_owned(),
            bases: 1,
        });
        document
            .validate_structure()
            .expect("dos bloques que suman el tracto");
    }
}

// ---------------------------------------------------------------------------
// Cabecera, copias y segmentos
// ---------------------------------------------------------------------------

/// Sin huellas de fuentes no hay nada que contrastar.
#[test]
fn a_header_without_source_digests_is_rejected() {
    assert_table(vec![(
        "sin huellas",
        edit(|d| d.header.source_digests.clear()),
        TruthContractError::EmptySourceDigests,
    )]);
}

/// Una copia de otra referencia que la cabecera, con una fase que no es su
/// índice, con una procedencia repetida, sin procedencia o sin segmentos.
#[test]
fn a_copy_must_agree_with_the_header_and_have_provenance_and_segments() {
    let at = |index: u8, reference: &str| CopyIdentity {
        reference_id: reference.to_owned(),
        individual_id: "I1".to_owned(),
        sequence_id: "chr1".to_owned(),
        copy_index: index,
    };
    let source = rich(&mut base()).provenance.sources[0].clone();
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "referencia ajena",
            edit(|d| rich(d).identity.reference_id = "ref-other".to_owned()),
            TruthContractError::CopyReferenceMismatch {
                identity: at(0, "ref-other"),
            },
        ),
        (
            "fase que no es el índice",
            edit(|d| rich(d).phase = 1),
            TruthContractError::InvalidCopyPhase {
                identity: at(0, "ref-main"),
                phase: 1,
            },
        ),
        (
            "procedencia repetida",
            edit(|d| {
                let repeated = rich(d).provenance.sources[0].clone();
                rich(d).provenance.sources.push(repeated);
            }),
            lexical("copy_provenance_sources", &source.source_id),
        ),
        (
            "sin procedencia",
            edit(|d| rich(d).provenance.sources.clear()),
            TruthContractError::EmptySegment {
                identity: at(0, "ref-main"),
            },
        ),
        (
            "sin segmentos",
            edit(|d| rich(d).segments.clear()),
            TruthContractError::EmptySegment {
                identity: at(0, "ref-main"),
            },
        ),
        (
            "sin fragmentos vigentes",
            edit(|d| rich(d).current_pieces.clear()),
            TruthContractError::EmptyCurrentPieces {
                identity: at(0, "ref-main"),
            },
        ),
    ];
    assert_table(table);
}

fn founder_digest() -> NamedDigest {
    NamedDigest {
        source: TruthSourceRef {
            kind: TruthSourceKind::Founder,
            source_id: "FOUNDER_1".to_owned(),
        },
        sha256: "d".repeat(64),
    }
}

/// La fuente de un segmento tiene que ser de la referencia de la copia, de su
/// secuencia y de la clase `reference` o `founder`.
#[test]
fn a_segment_source_must_name_the_copys_reference_and_sequence() {
    let mismatch = TruthContractError::SegmentSourceMismatch {
        identity: identity("I1", 0),
    };
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "otra referencia",
            edit(|d| rich(d).segments[0].source.reference_id = "ref-other".to_owned()),
            mismatch.clone(),
        ),
        (
            "otra secuencia",
            edit(|d| rich(d).segments[0].source.sequence_id = "chr2".to_owned()),
            mismatch.clone(),
        ),
        (
            "una fuente explícita",
            edit(|d| {
                rich(d).segments[0].source.source = TruthSourceRef {
                    kind: TruthSourceKind::Explicit,
                    source_id: "ins-1".to_owned(),
                };
            }),
            mismatch,
        ),
    ];
    assert_table(table);
    // Lado aceptado: una fuente `founder` con su huella.
    assert_all_accepted(vec![(
        "fundador",
        edit(|d| {
            d.header.source_digests.push(founder_digest());
            rich(d).segments[0].source.source = founder_digest().source;
        }),
    )]);
}

/// El intervalo de un segmento: no vacío y dentro de la longitud declarada.
#[test]
fn a_segment_interval_is_non_empty_and_inside_the_declared_length() {
    let invalid = TruthContractError::InvalidSegmentInterval {
        identity: identity("I1", 0),
    };
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "vacío",
            edit(|d| {
                let source = &mut rich(d).segments[0].source;
                source.start = 10;
                source.end = 10;
            }),
            invalid.clone(),
        ),
        (
            "una base más allá de la longitud declarada",
            edit(|d| rich(d).segments[0].source.end = 601),
            invalid,
        ),
    ];
    assert_table(table);
    assert_all_accepted(vec![
        (
            "una sola base",
            edit(|d| {
                let copy = rich(d);
                copy.segments[0].modifications.clear();
                copy.segments[0].source.start = 10;
                copy.segments[0].source.end = 11;
            }),
        ),
        (
            "una base antes de la longitud declarada",
            edit(|d| rich(d).segments[0].source.end = 599),
        ),
    ]);
}

// ---------------------------------------------------------------------------
// Las modificaciones de un segmento
// ---------------------------------------------------------------------------

fn snp_at(event_id: &str, position: u64) -> TruthModification {
    typed_snp_modification(event_id, position)
}

fn insertion_at(event_id: &str, point: u64) -> TruthModification {
    TruthModification::Insertion {
        event_id: event_id.to_owned(),
        point,
        bases: "ACG".to_owned(),
        provenance: plain_provenance(),
    }
}

fn deletion_between(event_id: &str, start: u64, end: u64) -> TruthModification {
    TruthModification::Deletion {
        event_id: event_id.to_owned(),
        start,
        end,
        provenance: plain_provenance(),
    }
}

/// La copia 1 de I3 con un único segmento `[start, end)` y estas
/// modificaciones; sus fragmentos vigentes y su historial no se tocan (el
/// contrato no los cruza con los segmentos).
fn with_segment(start: u64, end: u64, modifications: Vec<TruthModification>) -> Edit {
    edit(move |d| {
        let copy = copy(d, "I3", 1);
        copy.segments = vec![typed_segment(start, end, modifications.clone())];
    })
}

/// Una modificación fuera del segmento `[10, 590)`, y las del borde que sí
/// están dentro. Una inserción puede estar en el final del segmento (el punto
/// es el límite entre dos bases); un cambio de base o una deleción no.
#[test]
fn a_modification_must_lie_inside_its_segment() {
    let outside = |event_id: &str| TruthContractError::ModificationOutsideSegment {
        event_id: event_id.to_owned(),
    };
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "cambio antes del segmento",
            with_segment(10, 590, vec![snp_at("E1", 9)]),
            outside("E1"),
        ),
        (
            "cambio en el final del segmento",
            with_segment(10, 590, vec![snp_at("E1", 590)]),
            outside("E1"),
        ),
        (
            "inserción antes del segmento",
            with_segment(10, 590, vec![insertion_at("E1", 9)]),
            outside("E1"),
        ),
        (
            "inserción después del segmento",
            with_segment(10, 590, vec![insertion_at("E1", 591)]),
            outside("E1"),
        ),
        (
            "deleción que empieza antes",
            with_segment(10, 590, vec![deletion_between("E1", 9, 20)]),
            outside("E1"),
        ),
        (
            "deleción que acaba después",
            with_segment(10, 590, vec![deletion_between("E1", 580, 591)]),
            outside("E1"),
        ),
    ];
    assert_table(table);
    let accepted: Vec<(&str, Edit)> = vec![
        (
            "cambio en el inicio",
            with_segment(10, 590, vec![snp_at("E1", 10)]),
        ),
        (
            "cambio en la última base",
            with_segment(10, 590, vec![snp_at("E1", 589)]),
        ),
        (
            "inserción en el inicio",
            with_segment(10, 590, vec![insertion_at("E1", 10)]),
        ),
        (
            "inserción en el final",
            with_segment(10, 590, vec![insertion_at("E1", 590)]),
        ),
        (
            "deleción desde el inicio",
            with_segment(10, 590, vec![deletion_between("E1", 10, 20)]),
        ),
        (
            "deleción hasta el final",
            with_segment(10, 590, vec![deletion_between("E1", 580, 590)]),
        ),
    ];
    assert_all_accepted(accepted);
}

/// Las modificaciones de un segmento están **ordenadas y no se solapan**, sin
/// excepción (DT-047, cerrada): ni dos deleciones solapadas, ni un cambio de
/// base o una inserción dentro de una deleción —en los dos órdenes—, ni a su
/// final. Cada rechazo con su motivo; y, al lado, lo que sí es válido
/// (contiguos, y el punto justo antes de una deleción o justo después de su
/// final).
#[test]
fn modifications_of_a_segment_are_ordered_and_never_overlap() {
    let overlapping = |id: &str| TruthContractError::OverlappingModifications {
        event_id: id.to_owned(),
    };
    let not_ordered = |id: &str| TruthContractError::ModificationsNotOrdered {
        event_id: id.to_owned(),
    };
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "dos deleciones solapadas",
            with_segment(
                0,
                600,
                vec![
                    deletion_between("E1", 100, 200),
                    deletion_between("E2", 150, 250),
                ],
            ),
            overlapping("E2"),
        ),
        (
            "una deleción dentro de otra",
            with_segment(
                0,
                600,
                vec![
                    deletion_between("E1", 100, 200),
                    deletion_between("E2", 120, 130),
                ],
            ),
            overlapping("E2"),
        ),
        (
            "cambio dentro de una deleción",
            with_segment(
                0,
                600,
                vec![deletion_between("E1", 100, 200), snp_at("E2", 150)],
            ),
            overlapping("E2"),
        ),
        (
            "cambio en la última base de la deleción",
            with_segment(
                0,
                600,
                vec![deletion_between("E1", 100, 200), snp_at("E2", 199)],
            ),
            overlapping("E2"),
        ),
        (
            "cambio en la primera base de la deleción, después",
            with_segment(
                0,
                600,
                vec![deletion_between("E1", 100, 200), snp_at("E2", 100)],
            ),
            not_ordered("E2"),
        ),
        (
            "cambio y después la deleción que lo cubre",
            with_segment(
                0,
                600,
                vec![snp_at("E1", 150), deletion_between("E2", 100, 200)],
            ),
            not_ordered("E2"),
        ),
        (
            "cambio en la última base y después la deleción",
            with_segment(
                0,
                600,
                vec![snp_at("E1", 199), deletion_between("E2", 100, 200)],
            ),
            not_ordered("E2"),
        ),
        (
            "cambio en la primera base y después la deleción",
            with_segment(
                0,
                600,
                vec![snp_at("E1", 100), deletion_between("E2", 100, 200)],
            ),
            overlapping("E2"),
        ),
        (
            "inserción dentro de una deleción",
            with_segment(
                0,
                600,
                vec![deletion_between("E1", 100, 200), insertion_at("E2", 150)],
            ),
            overlapping("E2"),
        ),
        (
            "inserción justo al final de la deleción",
            with_segment(
                0,
                600,
                vec![deletion_between("E1", 100, 200), insertion_at("E2", 200)],
            ),
            overlapping("E2"),
        ),
    ];
    assert_table(table);
    assert_all_accepted(vec![
        (
            "deleciones contiguas",
            with_segment(
                0,
                600,
                vec![
                    deletion_between("E1", 100, 200),
                    deletion_between("E2", 200, 250),
                ],
            ),
        ),
        (
            "dos inserciones en puntos distintos",
            with_segment(
                0,
                600,
                vec![insertion_at("E1", 100), insertion_at("E2", 101)],
            ),
        ),
        (
            "inserción justo en el inicio de una deleción",
            with_segment(
                0,
                600,
                vec![insertion_at("E1", 100), deletion_between("E2", 100, 200)],
            ),
        ),
        (
            "inserción estrictamente antes de una deleción",
            with_segment(
                0,
                600,
                vec![insertion_at("E1", 90), deletion_between("E2", 100, 200)],
            ),
        ),
        (
            "deleción seguida de una inserción tras su final",
            with_segment(
                0,
                600,
                vec![deletion_between("E1", 100, 200), insertion_at("E2", 201)],
            ),
        ),
        (
            "deleción seguida de un cambio tras su final",
            with_segment(
                0,
                600,
                vec![deletion_between("E1", 100, 200), snp_at("E2", 200)],
            ),
        ),
        (
            "dos cambios consecutivos",
            with_segment(0, 600, vec![snp_at("E1", 100), snp_at("E2", 101)]),
        ),
        (
            "cambio y deleción contigua",
            with_segment(
                0,
                600,
                vec![snp_at("E1", 99), deletion_between("E2", 100, 200)],
            ),
        ),
        (
            "inserción antes del cambio en la misma posición",
            with_segment(0, 600, vec![insertion_at("E1", 100), snp_at("E2", 100)]),
        ),
        (
            "cambio e inserción en la posición siguiente",
            with_segment(0, 600, vec![snp_at("E1", 99), insertion_at("E2", 100)]),
        ),
    ]);
}

// ---------------------------------------------------------------------------
// Los fragmentos vigentes y el historial
// ---------------------------------------------------------------------------

fn with_pieces(pieces: Vec<TruthCurrentPiece>) -> Edit {
    edit(move |d| copy(d, "I3", 1).current_pieces = pieces.clone())
}

fn bad_piece() -> TruthContractError {
    TruthContractError::InvalidCurrentPiece {
        identity: identity("I3", 1),
    }
}

fn piece_with_source(change: impl Fn(&mut SegmentSource)) -> TruthCurrentPiece {
    let mut source = typed_source(10, 100, 600);
    change(&mut source);
    TruthCurrentPiece::Source {
        source,
        snps: vec![],
    }
}

/// Cada defecto de un fragmento de fuente, con su lado aceptado.
#[test]
fn a_source_piece_must_be_a_non_empty_interval_of_the_copys_sequence() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "vacío",
            with_pieces(vec![piece_with_source(|s| s.end = s.start)]),
            bad_piece(),
        ),
        (
            "una base más allá de la longitud",
            with_pieces(vec![piece_with_source(|s| s.end = 601)]),
            bad_piece(),
        ),
        (
            "otra referencia",
            with_pieces(vec![piece_with_source(|s| {
                s.reference_id = "ref-other".into()
            })]),
            bad_piece(),
        ),
        (
            "otra secuencia",
            with_pieces(vec![piece_with_source(|s| s.sequence_id = "chr2".into())]),
            bad_piece(),
        ),
        (
            "fuente explícita",
            with_pieces(vec![piece_with_source(|s| {
                s.source = TruthSourceRef {
                    kind: TruthSourceKind::Explicit,
                    source_id: "ins-1".into(),
                };
            })]),
            bad_piece(),
        ),
    ];
    assert_table(table);
    assert_all_accepted(vec![
        (
            "una sola base",
            with_pieces(vec![piece_with_source(|s| s.end = s.start + 1)]),
        ),
        (
            "hasta la longitud declarada",
            with_pieces(vec![piece_with_source(|s| s.end = 600)]),
        ),
        (
            "fuente fundador",
            edit(|d| {
                d.header.source_digests.push(founder_digest());
                copy(d, "I3", 1).current_pieces =
                    vec![piece_with_source(|s| s.source = founder_digest().source)];
            }),
        ),
    ]);
}

fn piece_with_snps(positions: &[u64]) -> TruthCurrentPiece {
    TruthCurrentPiece::Source {
        source: typed_source(10, 100, 600),
        snps: positions
            .iter()
            .enumerate()
            .map(|(index, position)| typed_snp(&format!("S{index}"), *position))
            .collect(),
    }
}

/// Los cambios de base de un fragmento: dentro de su intervalo, distintos de la
/// referencia y en orden estrictamente creciente.
#[test]
fn the_base_changes_of_a_piece_are_inside_it_different_and_strictly_increasing() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "antes del fragmento",
            with_pieces(vec![piece_with_snps(&[9])]),
            bad_piece(),
        ),
        (
            "en el final del fragmento",
            with_pieces(vec![piece_with_snps(&[100])]),
            bad_piece(),
        ),
        (
            "alterna igual a la referencia",
            with_pieces(vec![TruthCurrentPiece::Source {
                source: typed_source(10, 100, 600),
                snps: vec![{
                    let mut snp = typed_snp("S0", 50);
                    snp.alternate = snp.reference.clone();
                    snp
                }],
            }]),
            bad_piece(),
        ),
        (
            "dos en la misma posición",
            with_pieces(vec![piece_with_snps(&[20, 20])]),
            bad_piece(),
        ),
        (
            "en orden decreciente",
            with_pieces(vec![piece_with_snps(&[30, 20])]),
            bad_piece(),
        ),
    ];
    assert_table(table);
    assert_all_accepted(vec![
        (
            "en el inicio y en la última base",
            with_pieces(vec![piece_with_snps(&[10, 99])]),
        ),
        (
            "consecutivos",
            with_pieces(vec![piece_with_snps(&[20, 21])]),
        ),
    ]);
}

fn repeat_piece_with(change: impl Fn(&mut SegmentSource)) -> TruthCurrentPiece {
    let mut source = typed_source(10, 14, 600);
    change(&mut source);
    TruthCurrentPiece::RepeatAllele {
        locus_id: "STR-1".to_owned(),
        source,
        bases: "CACA".to_owned(),
        blocks: vec![TruthRepeatBlock {
            motif: "CA".to_owned(),
            bases: 4,
        }],
        provenance: plain_provenance(),
    }
}

/// Cada defecto de un tracto vigente, con su lado aceptado.
#[test]
fn a_realized_tract_piece_must_be_a_non_empty_interval_of_the_copys_sequence() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "vacío",
            with_pieces(vec![repeat_piece_with(|s| s.end = s.start)]),
            bad_piece(),
        ),
        (
            "más allá de la longitud",
            with_pieces(vec![repeat_piece_with(|s| s.end = 601)]),
            bad_piece(),
        ),
        (
            "otra referencia",
            with_pieces(vec![repeat_piece_with(|s| {
                s.reference_id = "ref-other".into()
            })]),
            bad_piece(),
        ),
        (
            "otra secuencia",
            with_pieces(vec![repeat_piece_with(|s| s.sequence_id = "chr2".into())]),
            bad_piece(),
        ),
        (
            "fuente explícita",
            with_pieces(vec![repeat_piece_with(|s| {
                s.source = TruthSourceRef {
                    kind: TruthSourceKind::Explicit,
                    source_id: "ins-1".into(),
                };
            })]),
            bad_piece(),
        ),
    ];
    assert_table(table);
    assert_all_accepted(vec![
        (
            "intervalo de una base",
            with_pieces(vec![repeat_piece_with(|s| s.end = s.start + 1)]),
        ),
        (
            "hasta la longitud declarada",
            with_pieces(vec![repeat_piece_with(|s| {
                s.start = 596;
                s.end = 600;
            })]),
        ),
    ]);
}

fn with_history(events: Vec<TruthHistoryEvent>) -> Edit {
    edit(move |d| copy(d, "I3", 1).history = events.clone())
}

fn history_insertion(event_id: &str, point: u64) -> TruthHistoryEvent {
    TruthHistoryEvent::Insertion {
        event_id: event_id.to_owned(),
        point,
        bases: "ACG".to_owned(),
        provenance: plain_provenance(),
    }
}

fn history_deletion(event_id: &str, start: u64, end: u64) -> TruthHistoryEvent {
    TruthHistoryEvent::Deletion {
        event_id: event_id.to_owned(),
        start,
        end,
        provenance: plain_provenance(),
    }
}

fn history_repeat(event_id: &str, start: u64, end: u64) -> TruthHistoryEvent {
    TruthHistoryEvent::RepeatAllele {
        event_id: event_id.to_owned(),
        locus_id: "STR-1".to_owned(),
        start,
        end,
        bases: "CACA".to_owned(),
        blocks: vec![TruthRepeatBlock {
            motif: "CA".to_owned(),
            bases: 4,
        }],
        provenance: plain_provenance(),
    }
}

/// Un evento del historial que no es un evento: una inserción en el punto
/// cero, una deleción que empieza en la base cero o que no tiene bases, un
/// tracto sin intervalo; y los mismos eventos en su límite aceptado.
#[test]
fn a_history_event_must_be_a_real_event() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "inserción en el punto cero",
            with_history(vec![history_insertion("H1", 0)]),
            bad_piece(),
        ),
        (
            "deleción desde la base cero",
            with_history(vec![history_deletion("H1", 0, 5)]),
            bad_piece(),
        ),
        (
            "deleción vacía",
            with_history(vec![history_deletion("H1", 5, 5)]),
            bad_piece(),
        ),
        (
            "deleción invertida",
            with_history(vec![history_deletion("H1", 6, 5)]),
            bad_piece(),
        ),
        (
            "tracto vacío",
            with_history(vec![history_repeat("H1", 5, 5)]),
            bad_piece(),
        ),
        (
            "tracto invertido",
            with_history(vec![history_repeat("H1", 6, 5)]),
            bad_piece(),
        ),
    ];
    assert_table(table);
    assert_all_accepted(vec![
        (
            "inserción en el punto uno",
            with_history(vec![history_insertion("H1", 1)]),
        ),
        (
            "deleción de una base desde la uno",
            with_history(vec![history_deletion("H1", 1, 2)]),
        ),
        (
            "tracto de una base",
            with_history(vec![history_repeat("H1", 5, 6)]),
        ),
    ]);
}

/// Dos eventos con el mismo identificador, sean de la clase que sean.
#[test]
fn a_history_event_identifier_is_not_repeated() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "dos inserciones",
            with_history(vec![history_insertion("H1", 1), history_insertion("H1", 2)]),
            TruthContractError::DuplicateHistoryEvent {
                event_id: "H1".to_owned(),
            },
        ),
        (
            "una inserción y una deleción",
            with_history(vec![
                history_insertion("H1", 1),
                history_deletion("H1", 1, 2),
            ]),
            TruthContractError::DuplicateHistoryEvent {
                event_id: "H1".to_owned(),
            },
        ),
    ];
    assert_table(table);
    assert_all_accepted(vec![(
        "dos identificadores distintos",
        with_history(vec![history_insertion("H1", 1), history_insertion("H2", 1)]),
    )]);
}

// ---------------------------------------------------------------------------
// La procedencia de cada modificación
// ---------------------------------------------------------------------------

type Provenance = fn(&mut TruthDocument) -> &mut ModificationProvenance;

fn prov_segment_explicit(document: &mut TruthDocument) -> &mut ModificationProvenance {
    snp_modification_of(rich(document)).3
}
fn prov_segment_matrix(document: &mut TruthDocument) -> &mut ModificationProvenance {
    rich(document)
        .segments
        .iter_mut()
        .flat_map(|segment| segment.modifications.iter_mut())
        .find_map(|modification| match modification {
            TruthModification::Snp { provenance, .. } if provenance.matrix.is_some() => {
                Some(provenance)
            }
            _ => None,
        })
        .expect("un cambio de la matriz")
}
fn prov_current_snp(document: &mut TruthDocument) -> &mut ModificationProvenance {
    &mut source_piece_with_snps(rich(document)).1[0].provenance
}
fn prov_inserted(document: &mut TruthDocument) -> &mut ModificationProvenance {
    inserted_piece_of(rich(document)).1
}
fn prov_repeat_piece(document: &mut TruthDocument) -> &mut ModificationProvenance {
    repeat_piece_of(rich(document)).4
}
fn prov_history_snp(document: &mut TruthDocument) -> &mut ModificationProvenance {
    let at = history_of(document, "I1", 0, "snp");
    provenance_of_history(&mut rich(document).history[at])
}
fn prov_history_insertion(document: &mut TruthDocument) -> &mut ModificationProvenance {
    let at = history_of(document, "I1", 0, "insertion");
    provenance_of_history(&mut rich(document).history[at])
}
fn prov_history_repeat(document: &mut TruthDocument) -> &mut ModificationProvenance {
    let at = history_of(document, "I1", 0, "repeat");
    provenance_of_history(&mut rich(document).history[at])
}

fn matrix_snp_event_id(document: &mut TruthDocument) -> String {
    rich(document)
        .segments
        .iter()
        .flat_map(|segment| segment.modifications.iter())
        .find_map(|modification| match modification {
            TruthModification::Snp {
                event_id,
                provenance,
                ..
            } if provenance.matrix.is_some() => Some(event_id.clone()),
            _ => None,
        })
        .expect("un cambio de la matriz")
}

fn history_event_id(document: &mut TruthDocument, kind: &str) -> String {
    let at = history_of(document, "I1", 0, kind);
    event_id_of_history(&mut rich(document).history[at]).clone()
}

/// Cada incoherencia de una procedencia (versión, clase de la fuente respecto
/// de la clase de la procedencia, presencia de la matriz), en cada sitio donde
/// el documento lleva una. El error nombra el evento, o `inserted_piece` o
/// `repeat_allele_piece` para los fragmentos, de modo que también distingue
/// al llamante.
#[test]
fn an_incoherent_provenance_is_rejected_in_every_place() {
    type Owner = fn(&mut TruthDocument) -> String;
    let callers: [(&str, Provenance, Owner); 8] = [
        (
            "cambio de segmento explícito",
            prov_segment_explicit,
            |d| snp_modification_of(rich(d)).0.clone(),
        ),
        (
            "cambio de segmento de la matriz",
            prov_segment_matrix,
            matrix_snp_event_id,
        ),
        ("cambio vigente", prov_current_snp, |d| {
            source_piece_with_snps(rich(d)).1[0].event_id.clone()
        }),
        ("bases insertadas vigentes", prov_inserted, |_| {
            "inserted_piece".to_owned()
        }),
        ("tracto vigente", prov_repeat_piece, |_| {
            "repeat_allele_piece".to_owned()
        }),
        ("cambio del historial", prov_history_snp, |d| {
            history_event_id(d, "snp")
        }),
        ("inserción del historial", prov_history_insertion, |d| {
            history_event_id(d, "insertion")
        }),
        ("tracto del historial", prov_history_repeat, |d| {
            history_event_id(d, "repeat")
        }),
    ];
    for (place, prov, owner) in callers {
        let owner = owner(&mut base());
        let invalid = TruthContractError::InvalidModificationProvenance {
            event_id: owner.clone(),
        };
        let table: Vec<(String, Edit, TruthContractError)> = vec![
            (
                format!("{place}: versión 2"),
                edit(move |d| prov(d).version = "2".to_owned()),
                invalid.clone(),
            ),
            (
                format!("{place}: clase de la fuente que no es la de la procedencia"),
                edit(move |d| {
                    let provenance = prov(d);
                    provenance.source.kind = if provenance.source.kind == TruthSourceKind::Explicit
                    {
                        TruthSourceKind::Random
                    } else {
                        TruthSourceKind::Explicit
                    };
                }),
                invalid.clone(),
            ),
            (
                format!("{place}: la matriz no corresponde a la clase"),
                edit(move |d| {
                    let provenance = prov(d);
                    provenance.matrix = if provenance.matrix.is_some() {
                        None
                    } else {
                        Some(
                            genome_population_simulator::truth_schema::MatrixProvenance {
                                matrix_source_id: provenance.source.source_id.clone(),
                                population_id: "P1".to_owned(),
                                variant_id: "VAR_1".to_owned(),
                            },
                        )
                    };
                }),
                invalid.clone(),
            ),
        ];
        base()
            .validate_structure()
            .expect("el documento intacto se acepta");
        for (label, change, expected) in table {
            let mut altered = base();
            change(&mut altered);
            assert_ne!(altered, base(), "{label}: la alteración alteró");
            assert_eq!(
                altered.validate_structure().err(),
                Some(expected),
                "{label}"
            );
        }
    }
}

/// La matriz de una procedencia de clase matriz: su fuente es la de la
/// procedencia, y su población tiene huella de grupo en la cabecera.
#[test]
fn a_matrix_provenance_names_its_source_and_has_a_group_digest() {
    let owner = matrix_snp_event_id(&mut base());
    let population = {
        let mut document = base();
        prov_segment_matrix(&mut document)
            .matrix
            .clone()
            .expect("matriz")
            .population_id
    };
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "la matriz nombra otra fuente",
            edit(|d| {
                prov_segment_matrix(d)
                    .matrix
                    .as_mut()
                    .expect("matriz")
                    .matrix_source_id = "MATRIX_OTHER".to_owned();
            }),
            TruthContractError::InvalidModificationProvenance {
                event_id: owner.clone(),
            },
        ),
        (
            "sin huella de grupo para la población",
            edit(move |d| {
                let population = prov_segment_matrix(d)
                    .matrix
                    .clone()
                    .expect("matriz")
                    .population_id;
                d.header
                    .group_digests
                    .retain(|digest| digest.source.source_id != population);
            }),
            TruthContractError::MissingDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Matrix,
                    source_id: population,
                },
            },
        ),
    ];
    assert_table(table);
}

// ---------------------------------------------------------------------------
// Los casos genotípicos
// ---------------------------------------------------------------------------

fn snp_case(
    document: &mut TruthDocument,
) -> &mut genome_population_simulator::truth_schema::GenotypeRecord {
    case(document, "I1", GenotypeCaseKind::Snp, 49)
}
fn insertion_case(
    document: &mut TruthDocument,
) -> &mut genome_population_simulator::truth_schema::GenotypeRecord {
    case(document, "I1", GenotypeCaseKind::Indel, 350)
}
fn deletion_case(
    document: &mut TruthDocument,
) -> &mut genome_population_simulator::truth_schema::GenotypeRecord {
    case(document, "I1", GenotypeCaseKind::Indel, 420)
}

fn case_id_of(
    select: fn(
        &mut TruthDocument,
    ) -> &mut genome_population_simulator::truth_schema::GenotypeRecord,
) -> String {
    select(&mut base()).case_id.clone()
}

fn invalid_allele(
    select: fn(
        &mut TruthDocument,
    ) -> &mut genome_population_simulator::truth_schema::GenotypeRecord,
) -> TruthContractError {
    TruthContractError::InvalidGenotypeAllele {
        case_id: case_id_of(select),
    }
}

/// Un caso sin alelos, o con tantos alelos como dosis no hay.
#[test]
fn a_genotype_has_alleles_and_one_dosage_for_each() {
    let empty = TruthContractError::EmptyGenotype {
        case_id: case_id_of(snp_case),
    };
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "sin alelos ni dosis",
            edit(|d| {
                let case = snp_case(d);
                case.alleles.clear();
                case.allele_dosages.clear();
            }),
            empty.clone(),
        ),
        (
            "una dosis de más",
            edit(|d| snp_case(d).allele_dosages.push(0)),
            empty.clone(),
        ),
        (
            "una dosis de menos",
            edit(|d| {
                snp_case(d).allele_dosages.pop();
            }),
            empty,
        ),
    ];
    assert_table(table);
}

/// Un caso de otra referencia que la de la cabecera.
#[test]
fn a_genotype_of_another_reference_is_rejected() {
    assert_table(vec![(
        "otra referencia",
        edit(|d| snp_case(d).reference_id = "ref-other".to_owned()),
        TruthContractError::GenotypeReferenceMismatch {
            case_id: case_id_of(snp_case),
        },
    )]);
}

/// Cada forma de un caso que su clase no admite: un locus invertido, un SNP
/// cuyo locus no mide una base, un caso indel sin alelo de inserción ni de
/// deleción, un caso de repetición cuyo primer alelo no es el de referencia o
/// que no lleva ningún alelo de repetición.
#[test]
fn a_case_must_have_the_shape_of_its_kind() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "locus invertido",
            edit(|d| {
                let case = insertion_case(d);
                case.locus.start = 351;
                case.locus.end = 350;
            }),
            invalid_allele(insertion_case),
        ),
        (
            "SNP de dos bases",
            edit(|d| snp_case(d).locus.end += 1),
            invalid_allele(snp_case),
        ),
        (
            "SNP de cero bases",
            edit(|d| snp_case(d).locus.end -= 1),
            invalid_allele(snp_case),
        ),
        (
            "indel sin alelo de inserción ni de deleción",
            edit(|d| {
                let case = insertion_case(d);
                case.alleles[1] = TruthAllele::Snp {
                    base: "A".to_owned(),
                };
            }),
            invalid_allele(insertion_case),
        ),
        (
            "repetición cuyo primer alelo no es la referencia",
            edit(|d| {
                let case = first_repeat_case(d);
                case.alleles.swap(0, 1);
                for carrier in &mut case.carriers {
                    carrier.allele_index = 1 - carrier.allele_index;
                }
                case.allele_dosages.swap(0, 1);
            }),
            TruthContractError::InvalidGenotypeAllele {
                case_id: first_repeat_case(&mut base()).case_id.clone(),
            },
        ),
        (
            "repetición sin alelo de repetición",
            edit(|d| {
                let case = first_repeat_case(d);
                case.alleles.truncate(1);
                case.allele_dosages = vec![2];
                for carrier in &mut case.carriers {
                    carrier.allele_index = 0;
                }
            }),
            TruthContractError::InvalidGenotypeAllele {
                case_id: first_repeat_case(&mut base()).case_id.clone(),
            },
        ),
    ];
    assert_table(table);
}

/// Un locus de un solo punto (inicio igual al final) es el de una inserción y
/// el de un caso indel, y se acepta; el SNP mide exactamente una base.
#[test]
fn a_zero_width_locus_is_the_locus_of_an_insertion() {
    // La inserción publicada ya lo es (350..350); se comprueba sobre ella que
    // un locus de ancho no nulo para una inserción se rechaza.
    assert_table(vec![(
        "inserción con un locus de una base",
        edit(|d| insertion_case(d).locus.end = 351),
        invalid_allele(insertion_case),
    )]);
}

fn other_kind_allele(kind: &str) -> TruthAllele {
    match kind {
        "snp" => TruthAllele::Snp {
            base: "A".to_owned(),
        },
        "repetición" => TruthAllele::Repeat {
            locus_id: "STR-1".to_owned(),
            bases: "CACA".to_owned(),
            blocks: vec![TruthRepeatBlock {
                motif: "CA".to_owned(),
                bases: 4,
            }],
        },
        "inserción" => insertion_allele_of_the_document(),
        "deleción" => deletion_allele_of_the_document(),
        other => panic!("clase desconocida {other}"),
    }
}

fn insertion_allele_of_the_document() -> TruthAllele {
    insertion_case(&mut base()).alleles[1].clone()
}

fn deletion_allele_of_the_document() -> TruthAllele {
    deletion_case(&mut base()).alleles[1].clone()
}

/// Un alelo de una clase que el caso no admite: un caso SNP solo lleva
/// referencia y SNP; un caso indel, referencia, inserción y deleción; un caso
/// de repetición, referencia y repetición. Cada clase de alelo en cada clase
/// de caso que no es la suya.
#[test]
fn an_allele_must_be_of_a_kind_its_case_admits() {
    type Select =
        fn(&mut TruthDocument) -> &mut genome_population_simulator::truth_schema::GenotypeRecord;
    let cases: [(&str, Select, &[&str]); 3] = [
        ("SNP", snp_case, &["inserción", "deleción", "repetición"]),
        ("indel", insertion_case, &["snp", "repetición"]),
        (
            "repetición",
            first_repeat_case,
            &["snp", "inserción", "deleción"],
        ),
    ];
    for (case_label, select, foreign) in cases {
        for kind in foreign {
            let allele = other_kind_allele(kind);
            let expected = invalid_allele(select);
            let label = format!("caso {case_label} con un alelo de {kind}");
            let mut altered = base();
            let record = select(&mut altered);
            record.alleles.push(allele);
            record.allele_dosages.push(0);
            assert_eq!(
                altered.validate_structure().err(),
                Some(expected),
                "{label}"
            );
        }
    }
}

/// Dos alelos iguales en un caso.
#[test]
fn a_genotype_does_not_repeat_an_allele() {
    assert_table(vec![(
        "alelo repetido",
        edit(|d| {
            let case = snp_case(d);
            let repeated = case.alleles[1].clone();
            case.alleles.push(repeated);
            case.allele_dosages.push(0);
        }),
        TruthContractError::DuplicateAllele {
            case_id: case_id_of(snp_case),
        },
    )]);
}

/// Un portador cuyo alelo no existe: el límite es el número de alelos.
#[test]
fn a_carrier_allele_index_must_be_inside_the_alleles() {
    let case_id = case_id_of(snp_case);
    assert_table(vec![(
        "índice igual al número de alelos",
        edit(|d| snp_case(d).carriers[0].allele_index = 2),
        TruthContractError::CarrierAlleleOutOfRange {
            case_id,
            allele_index: 2,
        },
    )]);
    assert_all_accepted(vec![(
        "índice del último alelo",
        edit(|d| {
            // El portador 0 pasa al último alelo y las dosis acompañan.
            let case = snp_case(d);
            let from = case.carriers[0].allele_index;
            case.carriers[0].allele_index = 1;
            case.allele_dosages[from] -= 1;
            case.allele_dosages[1] += 1;
        }),
    )]);
}

/// El inicio de un SNP en el último valor representable: no se le puede sumar
/// uno para hallar el final esperado.
#[test]
fn a_snp_locus_at_the_last_position_overflows() {
    assert_table(vec![(
        "SNP en la última posición",
        edit(|d| {
            let case = snp_case(d);
            case.locus.start = u64::MAX;
            case.locus.end = u64::MAX;
        }),
        TruthContractError::ArithmeticOverflow {
            field: "snp_locus_start_plus_one",
        },
    )]);
}

/// Un caso repetido, y un recuento de casos que no es el declarado.
#[test]
fn genotype_records_are_unique_and_counted() {
    let case_id = case_id_of(snp_case);
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "caso repetido",
            edit(|d| {
                let at = d
                    .records
                    .iter()
                    .position(|record| matches!(record, TruthRecord::Genotype(g) if g.case_id == case_id_of(snp_case)))
                    .expect("el caso");
                let repeated = d.records[at].clone();
                d.records.insert(at + 1, repeated);
            }),
            TruthContractError::DuplicateGenotype {
                case_id: case_id.clone(),
            },
        ),
        (
            "un caso de menos declarado",
            edit(|d| d.header.genotype_count -= 1),
            TruthContractError::GenotypeCountMismatch {
                declared: base().header.genotype_count - 1,
                actual: base().header.genotype_count,
            },
        ),
        (
            "una copia de más declarada",
            edit(|d| d.header.copy_count += 1),
            TruthContractError::CopyCountMismatch {
                declared: base().header.copy_count + 1,
                actual: base().header.copy_count,
            },
        ),
    ];
    assert_table(table);
}

// ---------------------------------------------------------------------------
// Los alelos de inserción y de deleción
// ---------------------------------------------------------------------------

fn insertion_parts(
    document: &mut TruthDocument,
) -> (
    &mut String,
    &mut TruthStructuralEvent,
    &mut Vec<TruthPieceLink>,
) {
    match &mut insertion_case(document).alleles[1] {
        TruthAllele::Insertion {
            bases,
            source_event,
            piece_links,
            ..
        } => (bases, source_event, piece_links),
        _ => panic!("una inserción"),
    }
}

fn deletion_parts(
    document: &mut TruthDocument,
) -> (
    &mut u64,
    &mut u64,
    &mut TruthStructuralEvent,
    &mut Vec<TruthPieceLink>,
) {
    match &mut deletion_case(document).alleles[1] {
        TruthAllele::Deletion {
            start,
            end,
            source_event,
            piece_links,
            ..
        } => (start, end, source_event, piece_links),
        _ => panic!("una deleción"),
    }
}

fn the_insertion_event_id() -> String {
    match &insertion_case(&mut base()).alleles[1] {
        TruthAllele::Insertion { event_id, .. } => event_id.clone(),
        _ => panic!("una inserción"),
    }
}

fn the_deletion_event_id() -> String {
    match &deletion_case(&mut base()).alleles[1] {
        TruthAllele::Deletion { event_id, .. } => event_id.clone(),
        _ => panic!("una deleción"),
    }
}

/// Un alelo de inserción coherente con su evento fuente: su evento es una
/// inserción, con un punto positivo, las mismas bases, y enlaza fragmentos.
#[test]
fn an_insertion_allele_must_agree_with_its_source_event() {
    let id = the_insertion_event_id();
    let invalid = TruthContractError::InvalidGenotypeAllele {
        case_id: id.clone(),
    };
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "el evento fuente es una deleción",
            edit({
                let id = id.clone();
                move |d| {
                    *insertion_parts(d).1 = TruthStructuralEvent::Deletion {
                        event_id: id.clone(),
                        start: 5,
                        end: 9,
                    };
                }
            }),
            invalid.clone(),
        ),
        (
            "punto cero",
            edit(|d| {
                if let TruthStructuralEvent::Insertion { point, .. } = insertion_parts(d).1 {
                    *point = 0;
                }
            }),
            invalid.clone(),
        ),
        (
            "bases distintas",
            edit(|d| {
                if let TruthStructuralEvent::Insertion { bases, .. } = insertion_parts(d).1 {
                    *bases = "T".to_owned();
                }
            }),
            invalid.clone(),
        ),
        (
            "sin fragmentos enlazados",
            edit(|d| insertion_parts(d).2.clear()),
            invalid,
        ),
        (
            "evento fuente de otro identificador",
            edit(|d| {
                if let TruthStructuralEvent::Insertion { event_id, .. } = insertion_parts(d).1 {
                    *event_id = "other-event".to_owned();
                }
            }),
            lexical("source_event_id", "other-event"),
        ),
        (
            "evento fuente con un identificador inválido",
            edit(|d| {
                if let TruthStructuralEvent::Insertion { event_id, .. } = insertion_parts(d).1 {
                    *event_id = BAD_ID.to_owned();
                }
            }),
            lexical("source_event_id", BAD_ID),
        ),
    ];
    assert_table(table);
    let _ = id;
}

/// Un alelo de deleción coherente con su evento fuente: su evento es una
/// deleción del mismo intervalo, y su propio intervalo es de bases reales.
#[test]
fn a_deletion_allele_must_agree_with_its_source_event() {
    let id = the_deletion_event_id();
    let invalid = TruthContractError::InvalidGenotypeAllele {
        case_id: id.clone(),
    };
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "el evento fuente es una inserción",
            edit(|d| {
                *deletion_parts(d).2 = TruthStructuralEvent::Insertion {
                    event_id: the_deletion_event_id(),
                    point: 5,
                    bases: "A".to_owned(),
                };
            }),
            invalid.clone(),
        ),
        (
            "el inicio fuente no es el del alelo",
            edit(|d| {
                if let TruthStructuralEvent::Deletion { start, .. } = deletion_parts(d).2 {
                    *start += 1;
                }
            }),
            invalid.clone(),
        ),
        (
            "el final fuente no es el del alelo",
            edit(|d| {
                if let TruthStructuralEvent::Deletion { end, .. } = deletion_parts(d).2 {
                    *end += 1;
                }
            }),
            invalid,
        ),
        (
            "evento fuente de otro identificador",
            edit(|d| {
                if let TruthStructuralEvent::Deletion { event_id, .. } = deletion_parts(d).2 {
                    *event_id = "other-event".to_owned();
                }
            }),
            lexical("source_event_id", "other-event"),
        ),
        (
            "intervalo desde la base cero",
            edit(|d| *deletion_parts(d).0 = 0),
            lexical("deletion_allele_interval", &format!("0:{}", 424)),
        ),
        (
            "intervalo vacío",
            edit(|d| *deletion_parts(d).0 = 424),
            lexical("deletion_allele_interval", "424:424"),
        ),
        (
            "intervalo invertido",
            edit(|d| *deletion_parts(d).0 = 430),
            lexical("deletion_allele_interval", "430:424"),
        ),
    ];
    assert_table(table);
}

/// Los fragmentos enlazados de un alelo estructural: cada uno no vacío, sin
/// solapar el anterior y de un fragmento distinto. Con una inserción y con una
/// deleción (el mismo ayudante, dos llamantes).
#[test]
fn the_piece_links_of_a_structural_allele_are_non_empty_ordered_and_distinct() {
    let link = |piece_index: u64, start: u64, end: u64| TruthPieceLink {
        piece_index,
        start,
        end,
    };
    let cases: [(&str, Vec<TruthPieceLink>, TruthContractError); 3] = [
        (
            "vacío",
            vec![link(1, 7, 7)],
            lexical("piece_links", "1:7:7"),
        ),
        (
            "invertido",
            vec![link(1, 8, 7)],
            lexical("piece_links", "1:8:7"),
        ),
        (
            "solapado con el anterior",
            vec![link(1, 5, 10), link(2, 9, 12)],
            lexical("piece_links", "2:9:12"),
        ),
    ];
    for (label, links, expected) in cases {
        let links_in = links.clone();
        let table: Vec<(&str, Edit, TruthContractError)> = vec![(
            "inserción",
            edit(move |d| *insertion_parts(d).2 = links_in.clone()),
            expected.clone(),
        )];
        assert_table(table);
        let links_del = links;
        let table: Vec<(&str, Edit, TruthContractError)> = vec![(
            "deleción",
            edit(move |d| *deletion_parts(d).3 = links_del.clone()),
            expected,
        )];
        assert_table(table);
        let _ = label;
    }
    // Un fragmento repetido en dos enlaces, y el lado aceptado: enlaces
    // contiguos de fragmentos distintos.
    let repeated = vec![link(1, 5, 10), link(1, 10, 12)];
    assert_table(vec![(
        "fragmento repetido",
        edit(move |d| *insertion_parts(d).2 = repeated.clone()),
        lexical("piece_links", "1:10:12"),
    )]);
    assert_all_accepted(vec![(
        "enlaces contiguos de fragmentos distintos",
        edit(move |d| *insertion_parts(d).2 = vec![link(1, 5, 10), link(2, 10, 12)]),
    )]);
}

/// Las guardas del contrato que ninguna entrada puede alcanzar, con el intento
/// ejecutado: lo que gana la carrera en cada caso. Esta prueba **no discrimina**
/// esas guardas; solo muestra qué otra comprobación rechaza antes.
#[test]
fn unreachable_guards_are_shadowed_by_earlier_checks() {
    // Unas bases insertadas vacías en un fragmento vigente: el paso léxico las
    // rechaza antes que el del estado vigente.
    let mut document = base();
    *inserted_piece_of(rich(&mut document)).0 = String::new();
    assert_eq!(
        document.validate_structure().err(),
        Some(lexical("current_inserted_bases", ""))
    );
    // Un tracto vigente sin bases: lo mismo.
    let mut document = base();
    *repeat_piece_of(rich(&mut document)).2 = String::new();
    assert_eq!(
        document.validate_structure().err(),
        Some(lexical("repeat_allele_bases", ""))
    );
    // Una inserción del historial sin bases: lo mismo.
    let mut document = base();
    let at = history_of(&mut document, "I1", 0, "insertion");
    if let TruthHistoryEvent::Insertion { bases, .. } = &mut rich(&mut document).history[at] {
        bases.clear();
    }
    assert_eq!(
        document.validate_structure().err(),
        Some(lexical("history_insertion_bases", ""))
    );
    // Un alelo de inserción con el punto cero en su evento fuente: el control
    // del alelo (`point > 0`) rechaza antes que el del evento.
    let mut document = base();
    if let TruthStructuralEvent::Insertion { point, .. } = insertion_parts(&mut document).1 {
        *point = 0;
    }
    assert_eq!(
        document.validate_structure().err(),
        Some(TruthContractError::InvalidGenotypeAllele {
            case_id: the_insertion_event_id(),
        })
    );
    // Un alelo de deleción con el inicio cero en su evento fuente: el control
    // del intervalo del alelo rechaza antes.
    let mut document = base();
    *deletion_parts(&mut document).0 = 0;
    if let TruthStructuralEvent::Deletion { start, .. } = deletion_parts(&mut document).2 {
        *start = 0;
    }
    assert_eq!(
        document.validate_structure().err(),
        Some(lexical("deletion_allele_interval", "0:424"))
    );
}

// ---------------------------------------------------------------------------
// Cabecera: clase, versión, convención y huellas repetidas
// ---------------------------------------------------------------------------

/// La clase del documento, la versión (que se exige **según lo que el
/// documento publica**, en los dos sentidos) y la convención de coordenadas.
#[test]
fn the_header_names_the_document_kind_the_version_it_publishes_and_the_convention() {
    assert_table(vec![
        (
            "otra clase de documento",
            edit(|d| d.header.document_kind = "other_kind".to_owned()),
            TruthContractError::UnsupportedDocumentKind,
        ),
        (
            "versión 4 con alelos de repetición",
            edit(|d| d.header.schema_version = "4".to_owned()),
            TruthContractError::UnsupportedVersion,
        ),
        (
            "otra versión",
            edit(|d| d.header.schema_version = "9".to_owned()),
            TruthContractError::UnsupportedVersion,
        ),
        (
            "otra convención de coordenadas",
            edit(|d| d.header.coordinate_convention = "one_based_closed".to_owned()),
            TruthContractError::UnsupportedCoordinateConvention,
        ),
    ]);
    // La versión 5 sin ningún alelo de repetición tampoco se admite.
    let mut without_repeats = document(&lines());
    without_repeats
        .validate_structure()
        .expect("la versión 4 intacta se acepta");
    without_repeats.header.schema_version = "5".to_owned();
    assert_eq!(
        without_repeats.validate_structure().err(),
        Some(TruthContractError::UnsupportedVersion)
    );
}

/// Una huella de fuente o de grupo repetida.
#[test]
fn a_repeated_digest_is_rejected() {
    let repeated = |digests: &Vec<NamedDigest>| TruthContractError::DuplicateDigest {
        source: digests[0].source.clone(),
    };
    let sources = repeated(&base().header.source_digests);
    let groups = repeated(&base().header.group_digests);
    assert_table(vec![
        (
            "huella de fuente repetida",
            edit(|d| {
                let again = d.header.source_digests[0].clone();
                d.header.source_digests.push(again);
            }),
            sources,
        ),
        (
            "huella de grupo repetida",
            edit(|d| {
                let again = d.header.group_digests[0].clone();
                d.header.group_digests.push(again);
            }),
            groups,
        ),
    ]);
}

// ---------------------------------------------------------------------------
// Orden y validez de las modificaciones de un segmento
// ---------------------------------------------------------------------------

/// Las modificaciones de un segmento van en orden estrictamente creciente.
#[test]
fn the_modifications_of_a_segment_are_strictly_ordered() {
    let not_ordered = |id: &str| TruthContractError::ModificationsNotOrdered {
        event_id: id.to_owned(),
    };
    assert_table(vec![
        (
            "decrecientes",
            with_segment(0, 600, vec![snp_at("E1", 100), snp_at("E2", 50)]),
            not_ordered("E2"),
        ),
        (
            "la misma posición",
            with_segment(0, 600, vec![snp_at("E1", 100), snp_at("E2", 100)]),
            not_ordered("E2"),
        ),
    ]);
    assert_all_accepted(vec![(
        "crecientes",
        with_segment(0, 600, vec![snp_at("E1", 50), snp_at("E2", 100)]),
    )]);
}

/// Cada modificación válida por sí misma: un cambio de base que cambia la
/// base, una inserción con punto entre 1 y la longitud menos uno, una
/// deleción entre la base 1 y la longitud menos uno.
#[test]
fn every_modification_must_be_valid_on_its_own() {
    let table: Vec<(&str, Edit, TruthContractError)> = vec![
        (
            "un cambio que no cambia la base",
            edit(|d| {
                let (_, reference, alternate, _) = snp_modification_of(rich(d));
                *alternate = reference.clone();
            }),
            TruthContractError::SnpDoesNotChangeBase {
                event_id: snp_modification_of(rich(&mut base())).0.clone(),
            },
        ),
        (
            "inserción en el punto cero",
            with_segment(0, 600, vec![insertion_at("E1", 0)]),
            TruthContractError::InvalidInsertion {
                event_id: "E1".to_owned(),
            },
        ),
        (
            "inserción en el punto igual a la longitud",
            with_segment(0, 600, vec![insertion_at("E1", 600)]),
            TruthContractError::InvalidInsertion {
                event_id: "E1".to_owned(),
            },
        ),
        (
            "deleción desde la base cero",
            with_segment(0, 600, vec![deletion_between("E1", 0, 5)]),
            TruthContractError::InvalidDeletion {
                event_id: "E1".to_owned(),
            },
        ),
        (
            "deleción hasta la longitud",
            with_segment(0, 600, vec![deletion_between("E1", 590, 600)]),
            TruthContractError::InvalidDeletion {
                event_id: "E1".to_owned(),
            },
        ),
        (
            "cambio en la última posición representable",
            edit(|d| {
                let copy = copy(d, "I3", 1);
                copy.segments = vec![TruthSegment {
                    source: typed_source(0, u64::MAX, u64::MAX),
                    modifications: vec![TruthModification::Snp {
                        event_id: "E1".to_owned(),
                        position: u64::MAX,
                        reference: "A".to_owned(),
                        alternate: "C".to_owned(),
                        provenance: plain_provenance(),
                    }],
                }];
            }),
            TruthContractError::ArithmeticOverflow {
                field: "snp_position_plus_one",
            },
        ),
    ];
    assert_table(table);
    assert_all_accepted(vec![
        (
            "inserción en el punto uno",
            with_segment(0, 600, vec![insertion_at("E1", 1)]),
        ),
        (
            "inserción en el último punto",
            with_segment(0, 600, vec![insertion_at("E1", 599)]),
        ),
        (
            "deleción desde la base uno",
            with_segment(0, 600, vec![deletion_between("E1", 1, 5)]),
        ),
        (
            "deleción hasta la última base",
            with_segment(0, 600, vec![deletion_between("E1", 590, 599)]),
        ),
    ]);
}

// ---------------------------------------------------------------------------
// Portadores, locus de un caso y alelos de inserción
// ---------------------------------------------------------------------------

/// Un portador que no es de este caso (otra referencia, otro individuo, otra
/// secuencia), que no existe como copia, o repetido.
#[test]
fn a_carrier_must_be_a_distinct_copy_of_the_cases_individual() {
    let owner = |d: &mut TruthDocument| snp_case(d).carriers[0].copy.clone();
    let carrier = owner(&mut base());
    let case_id = case_id_of(snp_case);
    let altered = |change: fn(&mut CopyIdentity)| {
        let mut copy = carrier.clone();
        change(&mut copy);
        copy
    };
    let outside = |copy: CopyIdentity| TruthContractError::CarrierOutsideGenotype {
        case_id: case_id.clone(),
        copy,
    };
    assert_table(vec![
        (
            "otra referencia",
            edit(|d| snp_case(d).carriers[0].copy.reference_id = "ref-other".to_owned()),
            outside(altered(|c| c.reference_id = "ref-other".to_owned())),
        ),
        (
            "otro individuo",
            edit(|d| snp_case(d).carriers[0].copy.individual_id = "I2".to_owned()),
            outside(altered(|c| c.individual_id = "I2".to_owned())),
        ),
        (
            "otra secuencia",
            edit(|d| snp_case(d).carriers[0].copy.sequence_id = "chr2".to_owned()),
            outside(altered(|c| c.sequence_id = "chr2".to_owned())),
        ),
        (
            "una copia que no existe",
            edit(|d| snp_case(d).carriers[0].copy.copy_index = 5),
            TruthContractError::UnknownCarrierCopy {
                case_id: case_id.clone(),
                copy: altered(|c| c.copy_index = 5),
            },
        ),
        (
            "el mismo portador dos veces",
            edit(|d| {
                let first = snp_case(d).carriers[0].copy.clone();
                snp_case(d).carriers[1].copy = first;
            }),
            TruthContractError::DuplicateCarrier {
                case_id: case_id.clone(),
                copy: carrier.clone(),
            },
        ),
    ]);
}

/// Un índice de copia de diez, solo, con la fase en cero.
#[test]
fn a_copy_index_of_ten_is_rejected_even_with_a_valid_phase() {
    assert_table(vec![(
        "índice 10 y fase 0",
        edit(|d| rich(d).identity.copy_index = 10),
        lexical("copy_index_or_phase", "10"),
    )]);
}

/// Un caso indel con el locus invertido que no es una inserción, y un caso
/// indel con solo el alelo de referencia.
#[test]
fn an_indel_case_needs_a_forward_locus_and_an_indel_allele() {
    assert_table(vec![
        (
            "deleción con el locus invertido",
            edit(|d| {
                let case = deletion_case(d);
                case.locus.start = 425;
                case.locus.end = 424;
            }),
            invalid_allele(deletion_case),
        ),
        (
            "indel con solo la referencia",
            edit(|d| {
                let case = insertion_case(d);
                case.alleles.truncate(1);
                case.allele_dosages = vec![2];
                for carrier in &mut case.carriers {
                    carrier.allele_index = 0;
                }
            }),
            invalid_allele(insertion_case),
        ),
    ]);
    assert_all_accepted(vec![(
        "inserción en el punto uno",
        edit(|d| {
            if let TruthStructuralEvent::Insertion { point, .. } = insertion_parts(d).1 {
                *point = 1;
            }
        }),
    )]);
}

// ---------------------------------------------------------------------------
// Lo que una deleción cubre y lo que no, y las deleciones sin intervalo
// ---------------------------------------------------------------------------

/// Una deleción sin intervalo no es una deleción: ni vacía ni invertida.
#[test]
fn an_empty_deletion_is_not_a_deletion() {
    let outside = |id: &str| TruthContractError::ModificationOutsideSegment {
        event_id: id.to_owned(),
    };
    assert_table(vec![
        (
            "deleción vacía",
            with_segment(0, 600, vec![deletion_between("E1", 100, 100)]),
            outside("E1"),
        ),
        (
            "deleción invertida",
            with_segment(0, 600, vec![deletion_between("E1", 110, 100)]),
            outside("E1"),
        ),
    ]);
}
