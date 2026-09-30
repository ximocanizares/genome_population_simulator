//! El tamaño de `truth.jsonl` no crece con la longitud del genoma por cada
//! *indel*.
//!
//! Hasta el 30 de septiembre de 2026, cada alelo de inserción o deleción de la
//! verdad copiaba la secuencia final entera de la copia portadora: con un
//! cromosoma de 500 kb y diez deleciones, el fichero medía 5.018.930 bytes y el
//! 99,6 % era esa copia. Estas comprobaciones **no miran el nombre de ningún
//! campo**: miden bytes, de modo que la copia no puede volver con otro nombre
//! sin romperlas.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

/// Longitud del cromosoma: lo bastante grande para que una sola copia de la
/// secuencia en un registro destaque sobre todo lo demás del fichero.
const CHROMOSOME_LENGTH: usize = 200_000;
/// Cota de cada registro de genotipo. Un registro de un *indel* sin la copia
/// mide unos cientos de bytes; con ella, más de `CHROMOSOME_LENGTH`.
const GENOTYPE_LINE_BOUND: usize = 4096;

/// **La verdad de un genoma de 200 kb con seis *indels* cabe en una cuarta
/// parte de su FASTA**, y ningún registro de genotipo se acerca a la longitud
/// de la copia.
///
/// Las dos cotas son independientes del nombre del campo:
///
/// - la del fichero entero falla con que **un solo** registro vuelva a llevar
///   la copia, porque esa copia ya mide lo que el FASTA;
/// - la de cada línea de genotipo señala cuál es el registro culpable.
///
/// Para que las cotas no pasen por vacío, se exige antes que la verdad
/// publique de verdad los seis casos de *indel*, en las dos copias.
#[test]
fn truth_size_does_not_grow_with_genome_length_per_indel() {
    let root = temporary("truth-size");
    fs::create_dir_all(&root).expect("fixture root");
    let config = write_fixture(&root);
    let destination = root.join("published");
    let effective = load_and_validate_config(&config).expect("valid configuration");
    execute_validated_config(&effective, &destination).expect("complete run");

    let truth_path = destination.join("outputs/truth.jsonl");
    let truth = fs::read_to_string(&truth_path).expect("truth.jsonl");
    let fasta_size = fs::metadata(root.join("reference.fa"))
        .expect("reference FASTA")
        .len();

    let records = truth
        .lines()
        .map(|line| {
            (
                line,
                serde_json::from_str::<Value>(line).expect("JSON line"),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(records[0].1["schema_version"], "4");
    let indel_cases = records
        .iter()
        .filter(|(_, record)| record["record_type"] == "genotype" && record["case_kind"] == "indel")
        .collect::<Vec<_>>();
    assert_eq!(
        indel_cases.len(),
        6,
        "the fixture must publish its six indels, or the bounds prove nothing"
    );
    for phase in [0, 1] {
        assert!(
            indel_cases.iter().any(|(_, record)| {
                record["carriers"]
                    .as_array()
                    .expect("carriers")
                    .iter()
                    .any(|carrier| {
                        carrier["allele_index"] == 1 && carrier["copy"]["copy_index"] == phase
                    })
            }),
            "copy {phase} must carry an indel"
        );
    }

    for (line, record) in &records {
        if record["record_type"] == "genotype" {
            assert!(
                line.len() < GENOTYPE_LINE_BOUND,
                "genotype {} measures {} bytes for a {CHROMOSOME_LENGTH}-base chromosome",
                record["case_id"],
                line.len()
            );
        }
    }
    let truth_size = u64::try_from(truth.len()).expect("size fits u64");
    assert!(
        truth_size * 4 <= fasta_size,
        "truth.jsonl measures {truth_size} bytes against a {fasta_size}-byte reference FASTA"
    );

    fs::remove_dir_all(&root).expect("fixture cleanup");
}

/// Escribe la referencia, la matriz y la configuración. Tres deleciones y tres
/// inserciones, repartidas entre las dos copias del único individuo diploide,
/// lejos del único sitio de la matriz.
fn write_fixture(root: &Path) -> PathBuf {
    let bases = pseudo_random_bases(CHROMOSOME_LENGTH);
    let mut fasta = String::with_capacity(CHROMOSOME_LENGTH + CHROMOSOME_LENGTH / 60 + 16);
    fasta.push_str(">chr1\n");
    for chunk in bases.chunks(60) {
        fasta.push_str(std::str::from_utf8(chunk).expect("ASCII bases"));
        fasta.push('\n');
    }
    fs::write(root.join("reference.fa"), fasta).expect("reference FASTA");

    let site = 1_000_usize;
    let reference_base = char::from(bases[site - 1]);
    let alternate_base = if reference_base == 'A' { 'C' } else { 'A' };
    fs::write(
        root.join("population-snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\n\
             chr1\t{site}\t{reference_base}\t{alternate_base}\tP1\t0.5\n"
        ),
    )
    .expect("population SNVs");

    let mut events = String::new();
    for (index, (phase, start)) in [(0, 20_000), (1, 60_000), (0, 100_000)]
        .into_iter()
        .enumerate()
    {
        events.push_str(&format!(
            "[[direct_events]]\nevent_id = \"del-{index}\"\nindividual_id = \"ind-1\"\n\
             phase = {phase}\nkind = \"deletion\"\ncontig = \"chr1\"\n\
             start = {start}\nend = {}\n\n",
            start + 5
        ));
    }
    for (index, (phase, point)) in [(1, 130_000), (0, 160_000), (1, 180_000)]
        .into_iter()
        .enumerate()
    {
        events.push_str(&format!(
            "[[direct_events]]\nevent_id = \"ins-{index}\"\nindividual_id = \"ind-1\"\n\
             phase = {phase}\nkind = \"insertion\"\ncontig = \"chr1\"\n\
             point = {point}\nbases = \"ACGTACGT\"\n\n"
        ));
    }
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{seed}\"\n\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"population-snvs.tsv\"\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = {CHROMOSOME_LENGTH}\n\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [\"ind-1\"]\n\n{events}",
            seed = "3".repeat(64),
        ),
    )
    .expect("configuration");
    config
}

/// Bases deterministas y sin estructura, para que ningún recorte o
/// alineamiento dependa de repeticiones accidentales.
fn pseudo_random_bases(length: usize) -> Vec<u8> {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[usize::try_from(state >> 62).expect("two bits")]
        })
        .collect()
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-truth-size-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
