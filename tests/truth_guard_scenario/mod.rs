//! Escenario compartido por las pruebas de las guardas de la verdad
//! (`truth_validation_guards.rs` y `truth_schema_guards.rs`): la verdad que el
//! programa publica de verdad para una referencia de 600 bases, y los
//! ayudantes para alterarla en un solo punto.
//!
//! Cada fichero de prueba que incluye este módulo usa solo una parte de él.
#![allow(dead_code)]

use std::fs;
use std::io::Cursor;
use std::path::PathBuf;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::truth_schema::{
    CopyIdentity, GenotypeCaseKind, ModificationProvenance, ModificationProvenanceKind,
    SegmentSource, TruthCurrentPiece, TruthCurrentSnp, TruthDocument, TruthModification,
    TruthRecord, TruthSegment, TruthSourceKind, TruthSourceRef,
};
use genome_population_simulator::truth_validation::{
    TruthValidationError, read_truth_jsonl, validate_and_reconstruct_truth,
    validate_and_reconstruct_truth_jsonl, validate_truth_jsonl_streaming,
};
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);

// ---------------------------------------------------------------------------
// La verdad publicada y su referencia
// ---------------------------------------------------------------------------

pub fn pseudo_random_bases(length: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[usize::try_from(state >> 62).expect("dos bits")]
        })
        .collect()
}

/// La primera secuencia, `chr1`, de 600 bases: aleatoria salvo un tracto
/// `(CA)50` en 200..300.
pub fn reference_bases() -> Vec<u8> {
    let mut bases = pseudo_random_bases(600, 0x2545_F491_4F6C_DD1D);
    for offset in 0..100 {
        bases[200 + offset] = if offset % 2 == 0 { b'C' } else { b'A' };
    }
    bases
}

/// La segunda secuencia, `chr2`, de 400 bases: aleatoria salvo un tracto
/// compuesto `(CAG)4(CA)9` en 100..130.
pub fn second_sequence_bases() -> Vec<u8> {
    let mut bases = pseudo_random_bases(400, 0x9E37_79B9_7F4A_7C15);
    let tract = format!("{}{}", "CAG".repeat(4), "CA".repeat(9));
    bases[100..130].copy_from_slice(tract.as_bytes());
    bases
}

pub fn reference_from(bases: &[u8]) -> Reference {
    let sequence = |id: &str, bases: &[u8]| {
        ReferenceSequence::new(
            SequenceId::parse(id, "truth guards").expect("sequence ID"),
            FastaHeader::parse(id.as_bytes()).expect("FASTA header"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(bases).expect("DNA"),
        )
        .expect("reference sequence")
    };
    Reference::new(
        ReferenceId::parse("ref-main", "truth guards").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![
            sequence("chr1", bases),
            sequence("chr2", &second_sequence_bases()),
        ],
    )
    .expect("reference")
}

pub fn reference() -> Reference {
    reference_from(&reference_bases())
}

/// Los dos loci de repetición: el tracto simple `(CA)50` de `chr1` con dos
/// alelos, y el tracto compuesto `(CAG)4(CA)9` de `chr2`, también con dos.
const REPEAT_LOCUS: &str = r#"[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 201
end = 300
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 45
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 40
frequency = 0.5

[[reference.repeat_loci]]
id = "STR-C"
sequence = "chr2"
start = 101
end = 130
motifs = ["CAG", "CA"]

[[reference.repeat_loci.alleles]]
units = [4, 5]
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = [4, 7]
frequency = 0.5
"#;

pub fn configuration(bases: &[u8], repeats: bool) -> (String, String) {
    let base = |position: usize| char::from(bases[position - 1]);
    let alternate = |value: char| if value == 'A' { 'C' } else { 'A' };
    let other = |value: char| if value == 'T' { 'G' } else { 'T' };
    let mut matrix = String::from("contig\tposition\tref\talt\tpopulation_id\talt_af\n");
    for population in ["P1", "P2"] {
        for position in [50, 120, 500] {
            matrix += &format!(
                "chr1\t{position}\t{}\t{}\t{population}\t0.5\n",
                base(position),
                alternate(base(position))
            );
        }
        matrix += &format!(
            "chr1\t120\t{}\t{}\t{population}\t0.25\n",
            base(120),
            other(base(120))
        );
    }
    let second = second_sequence_bases();
    let second_base = |position: usize| char::from(second[position - 1]);
    for population in ["P1", "P2"] {
        for position in [60, 250] {
            matrix += &format!(
                "chr2\t{position}\t{}\t{}\t{population}\t0.5\n",
                second_base(position),
                alternate(second_base(position))
            );
        }
    }
    let config = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{seed}"

[population_snvs]
mode = "supplied"
input_path = "population-snvs.tsv"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 600

[[reference.sequences]]
id = "chr2"
class = "chromosome"
length = 400

{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 2
ploidy = 2
genotype_model = "hwe"
individual_ids = ["I1", "I2"]

[[initial_populations]]
population_id = "P2"
individual_count = 1
ploidy = 2
genotype_model = "hwe"
individual_ids = ["I3"]

[[direct_events]]
event_id = "snp-direct"
individual_id = "I1"
phase = 0
kind = "snp"
contig = "chr1"
position_1based = 400
ref = "{r400}"
alt = "{a400}"

[[direct_events]]
event_id = "snp-boundary"
individual_id = "I1"
phase = 0
kind = "snp"
contig = "chr1"
position_1based = 351
ref = "{r351}"
alt = "{a351}"

[[direct_events]]
event_id = "ins-1"
individual_id = "I1"
phase = 0
kind = "insertion"
contig = "chr1"
point = 350
bases = "ACG"

[[direct_events]]
event_id = "ins-4"
individual_id = "I1"
phase = 0
kind = "insertion"
contig = "chr1"
point = 500
bases = "TT"

[[direct_events]]
event_id = "ins-5"
individual_id = "I2"
phase = 1
kind = "insertion"
contig = "chr1"
point = 450
bases = "CC"

[[direct_events]]
event_id = "del-span"
individual_id = "I2"
phase = 1
kind = "deletion"
contig = "chr1"
start = 448
end = 454

[[direct_events]]
event_id = "del-1"
individual_id = "I1"
phase = 1
kind = "deletion"
contig = "chr1"
start = 420
end = 424

[[direct_events]]
event_id = "ins-2"
individual_id = "I2"
phase = 0
kind = "insertion"
contig = "chr1"
point = 100
bases = "TTT"

[[direct_events]]
event_id = "del-2"
individual_id = "I2"
phase = 0
kind = "deletion"
contig = "chr1"
start = 101
end = 105

[[direct_events]]
event_id = "snp-direct-2"
individual_id = "I2"
phase = 0
kind = "snp"
contig = "chr2"
position_1based = 350
ref = "{r350}"
alt = "{a350}"

[[direct_events]]
event_id = "ins-first"
individual_id = "I1"
phase = 0
kind = "insertion"
contig = "chr2"
point = 1
bases = "A"

[[direct_events]]
event_id = "del-first"
individual_id = "I1"
phase = 1
kind = "deletion"
contig = "chr2"
start = 1
end = 3

[[direct_events]]
event_id = "ins-3"
individual_id = "I3"
phase = 1
kind = "insertion"
contig = "chr2"
point = 300
bases = "GT"

[[direct_events]]
event_id = "del-ins-only"
individual_id = "I3"
phase = 1
kind = "deletion"
contig = "chr2"
start = 301
end = 302
"#,
        loci = if repeats { REPEAT_LOCUS } else { "" },
        seed = "7".repeat(64),
        r400 = base(400),
        a400 = alternate(base(400)),
        r351 = base(351),
        a351 = alternate(base(351)),
        r350 = second_base(350),
        a350 = alternate(second_base(350)),
    );
    (config, matrix)
}

/// Las líneas de la verdad que el programa publica para ese escenario, con o
/// sin locus de repetición. **Sin él es la versión 4** y la leen los dos
/// lectores; **con él es la 5**, que el lector por flujo rechaza entera
/// (DT-015, abierta), de modo que esas pruebas usan solo el lector completo.
pub fn publish(repeats: bool) -> Vec<Value> {
    let root: PathBuf = std::env::temp_dir().join(format!(
        "gps-truth-guards-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("directorio");
    let bases = reference_bases();
    let (config, matrix) = configuration(&bases, repeats);
    fs::write(
        root.join("reference.fa"),
        format!(
            ">chr1\n{}\n>chr2\n{}\n",
            std::str::from_utf8(&bases).expect("ASCII"),
            std::str::from_utf8(&second_sequence_bases()).expect("ASCII")
        ),
    )
    .expect("FASTA");
    fs::write(root.join("population-snvs.tsv"), matrix).expect("matriz");
    fs::write(root.join("general.toml"), config).expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&root.join("general.toml")).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    let truth =
        fs::read_to_string(destination.join("outputs/truth.jsonl")).expect("verdad publicada");
    fs::remove_dir_all(root).expect("limpieza");
    truth
        .lines()
        .map(|line| serde_json::from_str(line).expect("línea JSON"))
        .collect()
}

/// La verdad de la versión 4 (sin repeticiones).
pub fn lines() -> Vec<Value> {
    static CELL: OnceLock<Vec<Value>> = OnceLock::new();
    CELL.get_or_init(|| publish(false)).clone()
}

/// La verdad de la versión 5 (con repeticiones).
pub fn lines_v5() -> Vec<Value> {
    static CELL: OnceLock<Vec<Value>> = OnceLock::new();
    CELL.get_or_init(|| publish(true)).clone()
}

pub fn jsonl(lines: &[Value]) -> String {
    lines
        .iter()
        .map(|line| serde_json::to_string(line).expect("serializar"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

pub fn identity(individual: &str, copy_index: u8) -> CopyIdentity {
    CopyIdentity {
        reference_id: "ref-main".to_owned(),
        individual_id: individual.to_owned(),
        sequence_id: "chr1".to_owned(),
        copy_index,
    }
}

/// Índice de la línea de la copia `copy_index` del individuo.
pub fn copy_line(lines: &[Value], individual: &str, copy_index: u8) -> usize {
    lines
        .iter()
        .position(|line| {
            line["record_type"] == "haplotype_copy"
                && line["identity"]["individual_id"] == individual
                && line["identity"]["copy_index"] == copy_index
        })
        .expect("la copia existe")
}

/// Índice del caso del individuo con ese tipo y ese inicio de locus.
pub fn case_line(lines: &[Value], individual: &str, kind: &str, start: u64) -> usize {
    lines
        .iter()
        .position(|line| {
            line["record_type"] == "genotype"
                && line["individual_id"] == individual
                && line["case_kind"] == kind
                && line["locus"]["start"] == start
        })
        .expect("el caso existe")
}

/// Lee con los dos lectores los mismos bytes: el completo y el de flujo.
pub fn both(bytes: &[u8]) -> [Result<(), TruthValidationError>; 2] {
    [
        validate_and_reconstruct_truth_jsonl(&mut Cursor::new(bytes.to_vec()), &reference())
            .map(|_| ()),
        validate_truth_jsonl_streaming(&mut Cursor::new(bytes.to_vec()), &reference(), |_, _| {})
            .map(|_| ()),
    ]
}

/// Solo el lector completo: la verdad de la versión 5 no la admite el de flujo.
pub fn full(bytes: &[u8]) -> Result<(), TruthValidationError> {
    validate_and_reconstruct_truth_jsonl(&mut Cursor::new(bytes.to_vec()), &reference()).map(|_| ())
}

/// Los dos lectores dan **exactamente** este error.
#[track_caller]
pub fn assert_both(label: &str, bytes: &[u8], expected: &TruthValidationError) {
    let [full, stream] = both(bytes);
    assert_eq!(
        full.as_ref().err(),
        Some(expected),
        "{label}: lector completo"
    );
    assert_eq!(
        stream.as_ref().err(),
        Some(expected),
        "{label}: lector por flujo"
    );
}

pub fn replay_failure(individual: &str, copy_index: u8) -> TruthValidationError {
    TruthValidationError::HistoryDoesNotReplay {
        identity: identity(individual, copy_index),
    }
}

// ---------------------------------------------------------------------------
// Ayudantes para armar copias a mano
// ---------------------------------------------------------------------------

pub fn reference_object(start: u64, end: u64, length: u64) -> Value {
    json!({
        "source": {"kind": "reference", "source_id": "ref-main"},
        "reference_id": "ref-main",
        "sequence_id": "chr1",
        "start": start,
        "end": end,
        "sequence_length": length,
    })
}

pub fn source_piece(start: u64, end: u64, snps: Vec<Value>) -> Value {
    json!({"kind": "source", "source": reference_object(start, end, 600), "snps": snps})
}

pub fn inserted_piece(bases: &str, provenance: &Value) -> Value {
    json!({"kind": "inserted", "bases": bases, "provenance": provenance})
}

pub fn segment(start: u64, end: u64, modifications: Vec<Value>) -> Value {
    json!({"source": reference_object(start, end, 600), "modifications": modifications})
}

/// Procedencia explícita de uno de los eventos directos del escenario, cuyas
/// huellas ya están en la cabecera.
pub fn explicit(event_id: &str) -> Value {
    json!({
        "kind": "explicit",
        "version": "1",
        "source": {"kind": "explicit", "source_id": event_id},
    })
}

pub fn reference_char(position: usize) -> char {
    char::from(reference_bases()[position])
}

pub fn other_base(base: char) -> char {
    if base == 'A' { 'C' } else { 'A' }
}

pub fn snp_modification(event_id: &str, position: usize, provenance: &Value) -> Value {
    let reference = reference_char(position);
    json!({
        "kind": "snp", "event_id": event_id, "position": position,
        "reference": reference.to_string(), "alternate": other_base(reference).to_string(),
        "provenance": provenance,
    })
}

pub fn snp_in_piece(event_id: &str, position: usize, provenance: &Value) -> Value {
    let reference = reference_char(position);
    json!({
        "event_id": event_id, "position": position,
        "reference": reference.to_string(), "alternate": other_base(reference).to_string(),
        "provenance": provenance,
    })
}

/// Sustituye los segmentos, los fragmentos vigentes y el historial de una copia.
pub fn set_copy(
    all: &mut [Value],
    individual: &str,
    copy_index: u8,
    segments: Vec<Value>,
    pieces: Vec<Value>,
    history: Vec<Value>,
) {
    let at = copy_line(all, individual, copy_index);
    all[at]["segments"] = Value::Array(segments);
    all[at]["current_pieces"] = Value::Array(pieces);
    all[at]["history"] = Value::Array(history);
}

pub fn document(all: &[Value]) -> TruthDocument {
    read_truth_jsonl(&mut Cursor::new(jsonl(all))).expect("el documento intacto se lee")
}

pub fn judge(document: &TruthDocument) -> Result<(), TruthValidationError> {
    validate_and_reconstruct_truth(document.clone(), &reference()).map(|_| ())
}

pub fn copy_of<'a>(
    document: &'a mut TruthDocument,
    individual: &str,
    copy_index: u8,
) -> &'a mut genome_population_simulator::truth_schema::HaplotypeCopyRecord {
    document
        .records
        .iter_mut()
        .find_map(|record| match record {
            TruthRecord::HaplotypeCopy(copy)
                if copy.identity.individual_id == individual
                    && copy.identity.copy_index == copy_index =>
            {
                Some(copy)
            }
            _ => None,
        })
        .expect("la copia existe")
}

pub fn genotype_of<'a>(
    document: &'a mut TruthDocument,
    individual: &str,
    kind: GenotypeCaseKind,
    start: u64,
) -> &'a mut genome_population_simulator::truth_schema::GenotypeRecord {
    document
        .records
        .iter_mut()
        .find_map(|record| match record {
            TruthRecord::Genotype(genotype)
                if genotype.individual_id == individual
                    && genotype.case_kind == kind
                    && genotype.locus.start == start =>
            {
                Some(genotype)
            }
            _ => None,
        })
        .expect("el caso existe")
}

pub fn source_ref(kind: TruthSourceKind, id: &str) -> TruthSourceRef {
    TruthSourceRef {
        kind,
        source_id: id.to_owned(),
    }
}

pub fn typed_provenance(id: &str) -> ModificationProvenance {
    ModificationProvenance {
        kind: ModificationProvenanceKind::Explicit,
        version: "1".to_owned(),
        source: source_ref(TruthSourceKind::Explicit, id),
        matrix: None,
    }
}

pub fn typed_source(start: u64, end: u64, length: u64) -> SegmentSource {
    SegmentSource {
        source: source_ref(TruthSourceKind::Reference, "ref-main"),
        reference_id: "ref-main".to_owned(),
        sequence_id: "chr1".to_owned(),
        start,
        end,
        sequence_length: length,
    }
}

pub fn typed_snp(event_id: &str, position: u64) -> TruthCurrentSnp {
    let reference = reference_char(usize::try_from(position).expect("cabe"));
    TruthCurrentSnp {
        event_id: event_id.to_owned(),
        position,
        reference: reference.to_string(),
        alternate: other_base(reference).to_string(),
        provenance: typed_provenance("snp-direct"),
    }
}

pub fn typed_snp_modification(event_id: &str, position: u64) -> TruthModification {
    let snp = typed_snp(event_id, position);
    TruthModification::Snp {
        event_id: snp.event_id,
        position,
        reference: snp.reference,
        alternate: snp.alternate,
        provenance: snp.provenance,
    }
}

pub fn typed_piece(start: u64, end: u64, snps: Vec<TruthCurrentSnp>) -> TruthCurrentPiece {
    TruthCurrentPiece::Source {
        source: typed_source(start, end, 600),
        snps,
    }
}

pub fn typed_segment(start: u64, end: u64, modifications: Vec<TruthModification>) -> TruthSegment {
    TruthSegment {
        source: typed_source(start, end, 600),
        modifications,
    }
}

pub fn bytes_of(lines: &[Value]) -> Vec<u8> {
    jsonl(lines).into_bytes()
}
