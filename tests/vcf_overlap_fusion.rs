//! DTC-064: las variantes de copias distintas que se solapan salen en **un solo
//! registro complejo**, y ninguna copia queda anotada con un alelo que no tiene.
//!
//! # Qué hace que estas pruebas puedan fallar
//!
//! El defecto que cubren lo compartían el productor y el juez semántico: los
//! dos asignaban `0` a toda copia sin alelo en un registro, y por eso el juez no
//! podía verlo. De ahí que **ninguna** prueba de este fichero se conforme con
//! que el juez acepte. Cada escenario pasa además por dos comprobaciones que no
//! comparten nada con el programa:
//!
//! - **un oráculo propio**, escrito aquí, que aplica el `GT` de cada registro a
//!   la referencia para reconstruir cada copia y la compara con la copia que
//!   reconstruye, también aquí, desde las piezas de `truth.jsonl`. Si el VCF
//!   anotase a una copia con un alelo que no tiene, su reconstrucción no
//!   coincidiría. El oráculo exige además que ningún registro pise a otro;
//! - **`bcftools norm` 1.24**, el juez externo de DTC-066, que tiene que dejar
//!   intacto todo registro que no sea de repetición. Los de repetición los
//!   recorta a propósito (excepción de DTC-064) y para ellos se exige que
//!   recortarlos no cambie ninguna copia reconstruida.
//!
//! Las aserciones sobre la forma concreta de cada registro —cuántos hay, en qué
//! posición, con qué `GT`— están escritas a mano desde el escenario, no leídas
//! de nada que el programa produzca.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::{
    validate_effective_config_against_manifest_bytes, validate_general_manifest_bytes,
    validate_published_resources,
};
use genome_population_simulator::vcf_semantic_validation::{
    VcfSemanticValidationError, validate_published_vcf_semantics,
};
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "5555555555555555555555555555555555555555555555555555555555555555";
const SAMPLES: [&str; 3] = ["ind-1", "ind-2", "ind-3"];

/// La versión del algoritmo de normalización, **escrita aquí** y no importada:
/// importarla convertiría la comprobación de lo publicado en un espejo.
const NORMALIZATION_VERSION: &str = "vcf-small-allele-normalization-v2";

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-vcf-fusion-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

// ---------------------------------------------------------------------------
// Escenarios
// ---------------------------------------------------------------------------

/// Bases pseudoaleatorias reproducibles, sin depender del generador del
/// programa.
fn pseudo_random_bases(length: usize, mut state: u64) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[usize::try_from(state >> 62).expect("two bits")]
        })
        .collect()
}

/// Referencia de 3000 pb con `CCCTCTCT` en las posiciones 1246 a 1253 en base
/// uno, el tramo del ejemplo de la revisión.
fn plain_reference() -> Vec<u8> {
    let mut bases = pseudo_random_bases(3000, 0x2545_F491_4F6C_DD1D);
    bases[1245..1253].copy_from_slice(b"CCCTCTCT");
    bases
}

/// Cien bases al azar terminadas en `G`, el tracto `(CA)50` en las bases 101 a
/// 200 en base uno, y el flanco que se pida.
fn tract_reference(flank: &str) -> Vec<u8> {
    let mut bases = pseudo_random_bases(99, 0x9E37_79B9_7F4A_7C15);
    bases.push(b'G');
    bases.extend_from_slice("CA".repeat(50).as_bytes());
    bases.extend_from_slice(flank.as_bytes());
    bases
}

/// Locus `(CA)50` con dos alelos, el de la referencia entre ellos, de modo que
/// unas copias conservan el tracto de la referencia y otras lo sustituyen.
const TRACT_LOCUS: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 50
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 45
frequency = 0.5
"#;

struct Event(String);

fn snp(individual: &str, phase: u8, position_1based: usize, bases: &[u8]) -> Event {
    let reference = char::from(bases[position_1based - 1]);
    let alternate = if reference == 'G' { 'T' } else { 'G' };
    snp_to(individual, phase, position_1based, reference, alternate)
}

fn snp_to(
    individual: &str,
    phase: u8,
    position_1based: usize,
    reference: char,
    alternate: char,
) -> Event {
    Event(format!(
        "individual_id = \"{individual}\"\nphase = {phase}\nkind = \"snp\"\nposition_1based = {position_1based}\nref = \"{reference}\"\nalt = \"{alternate}\"\n"
    ))
}

/// Deleción de `[start, end)` en coordenadas vigentes de base cero.
fn deletion(individual: &str, phase: u8, start: usize, end: usize) -> Event {
    Event(format!(
        "individual_id = \"{individual}\"\nphase = {phase}\nkind = \"deletion\"\nstart = {start}\nend = {end}\n"
    ))
}

/// Inserción delante de la base de base cero `point`.
fn insertion(individual: &str, phase: u8, point: usize, inserted: &str) -> Event {
    Event(format!(
        "individual_id = \"{individual}\"\nphase = {phase}\nkind = \"insertion\"\npoint = {point}\nbases = \"{inserted}\"\n"
    ))
}

/// Escribe el escenario: FASTA, una fila de matriz poblacional en la base 20
/// —los cambios directos la exigen— y la configuración.
fn write_scenario(root: &Path, bases: &[u8], loci: &str, events: &[Event], vcf: bool) -> PathBuf {
    fs::create_dir_all(root).expect("fixture root");
    let mut fasta = String::from(">chr1\n");
    for chunk in bases.chunks(60) {
        fasta.push_str(std::str::from_utf8(chunk).expect("ASCII"));
        fasta.push('\n');
    }
    fs::write(root.join("reference.fa"), fasta).expect("FASTA");
    let reference = char::from(bases[19]);
    let alternate = if reference == 'A' { 'C' } else { 'A' };
    fs::write(
        root.join("snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\t{reference}\t{alternate}\tP1\t0.5\n"
        ),
    )
    .expect("matrix");
    let mut config = format!(
        "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
         [population_snvs]\nmode = \"supplied\"\ninput_path = \"snvs.tsv\"\n\n\
         [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
         [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = {}\n{loci}\n\
         [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 3\nploidy = 2\n\
         genotype_model = \"hwe\"\nindividual_ids = [\"ind-1\", \"ind-2\", \"ind-3\"]\n",
        bases.len()
    );
    if vcf {
        config.push_str("\n[vcf_output]\nmode = \"text\"\n");
    }
    for (index, event) in events.iter().enumerate() {
        config.push_str(&format!(
            "\n[[direct_events]]\nevent_id = \"e{index}\"\ncontig = \"chr1\"\n{}",
            event.0
        ));
    }
    let path = root.join("general.toml");
    fs::write(&path, config).expect("configuration");
    path
}

struct Published {
    root: PathBuf,
    out: PathBuf,
}

impl Published {
    fn vcf_path(&self) -> PathBuf {
        self.out.join("outputs/truth.vcf")
    }
    fn truth_path(&self) -> PathBuf {
        self.out.join("outputs/truth.jsonl")
    }
    fn fasta_path(&self) -> PathBuf {
        self.out.join("resources/references/ref-main.fa")
    }
    fn vcf(&self) -> String {
        fs::read_to_string(self.vcf_path()).expect("VCF")
    }
    fn records(&self) -> Vec<Record> {
        records(&self.vcf())
    }
    fn json(&self, name: &str) -> Value {
        serde_json::from_slice(&fs::read(self.out.join(name)).expect("JSON")).expect("JSON")
    }
}

impl Drop for Published {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.root).ok();
    }
}

fn run(config: &Path, out: &Path) -> Result<(), Diagnostic> {
    let effective = load_and_validate_config(config)?;
    execute_validated_config(&effective, out).map(|_| ())
}

fn publish(label: &str, bases: &[u8], loci: &str, events: &[Event]) -> Published {
    let root = temporary(label);
    let config = write_scenario(&root, bases, loci, events, true);
    let out = root.join("published");
    run(&config, &out).unwrap_or_else(|error| panic!("{label}: {}", error.render_human()));
    Published { root, out }
}

// ---------------------------------------------------------------------------
// Lectura del VCF y del oráculo propio
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Record {
    line: String,
    position: usize,
    reference: String,
    alternates: Vec<String>,
    info: String,
    /// Por muestra, en el orden de la cabecera, los índices de alelo por fase.
    genotypes: Vec<Vec<usize>>,
}

impl Record {
    fn interval(&self) -> (usize, usize) {
        (self.position - 1, self.position - 1 + self.reference.len())
    }
    fn is_repeat(&self) -> bool {
        self.info.split(';').any(|field| field.starts_with("RN="))
    }
    fn gt(&self, sample: usize) -> String {
        self.genotypes[sample]
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join("|")
    }
}

fn records(text: &str) -> Vec<Record> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            Record {
                line: line.to_owned(),
                position: fields[1].parse().expect("POS"),
                reference: fields[3].to_owned(),
                alternates: fields[4].split(',').map(str::to_owned).collect(),
                info: fields[7].to_owned(),
                genotypes: fields[9..]
                    .iter()
                    .map(|gt| {
                        gt.split('|')
                            .map(|allele| allele.parse().expect("GT index"))
                            .collect()
                    })
                    .collect(),
            }
        })
        .collect()
}

fn fasta_bases(path: &Path) -> Vec<u8> {
    fs::read_to_string(path)
        .expect("FASTA")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect()
}

/// Cada copia según `truth.jsonl`, por (individuo, fase), desde sus piezas.
fn copies_from_truth(truth: &Path, reference: &[u8]) -> BTreeMap<(String, usize), Vec<u8>> {
    let mut copies = BTreeMap::new();
    for line in fs::read_to_string(truth).expect("truth").lines() {
        let record: Value = serde_json::from_str(line).expect("JSON line");
        if record["record_type"] != "haplotype_copy" {
            continue;
        }
        let mut sequence = Vec::new();
        for piece in record["current_pieces"].as_array().expect("pieces") {
            match piece["kind"].as_str().expect("kind") {
                "source" => {
                    let start = piece["source"]["start"].as_u64().expect("start") as usize;
                    let end = piece["source"]["end"].as_u64().expect("end") as usize;
                    let mut bases = reference[start..end].to_vec();
                    for snp in piece["snps"].as_array().expect("snps") {
                        let position = snp["position"].as_u64().expect("position") as usize;
                        bases[position - start] =
                            snp["alternate"].as_str().expect("alternate").as_bytes()[0];
                    }
                    sequence.extend(bases);
                }
                "inserted" | "repeat_allele" => {
                    sequence.extend(piece["bases"].as_str().expect("bases").bytes());
                }
                other => panic!("pieza desconocida {other}"),
            }
        }
        copies.insert(
            (
                record["identity"]["individual_id"]
                    .as_str()
                    .expect("individual")
                    .to_owned(),
                usize::try_from(record["phase"].as_u64().expect("phase")).expect("phase"),
            ),
            sequence,
        );
    }
    copies
}

/// Cada copia según el VCF: la referencia con el alelo que el `GT` de cada
/// registro le da. Exige que ningún registro pise a otro.
fn copies_from_vcf(
    text: &str,
    reference: &[u8],
) -> Result<BTreeMap<(String, usize), Vec<u8>>, String> {
    let records = records(text);
    let samples = text
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("cabecera de columnas")
        .split('\t')
        .skip(9)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let mut copies = BTreeMap::new();
    for (sample_index, sample) in samples.iter().enumerate() {
        let ploidy = records
            .first()
            .map_or(2, |record| record.genotypes[sample_index].len());
        for phase in 0..ploidy {
            let mut sequence = Vec::new();
            let mut cursor = 0;
            for record in &records {
                let (start, end) = record.interval();
                if start < cursor {
                    return Err(format!(
                        "dos registros se pisan en {}: {}",
                        record.position, record.line
                    ));
                }
                assert_eq!(&reference[start..end], record.reference.as_bytes());
                sequence.extend_from_slice(&reference[cursor..start]);
                let allele = record.genotypes[sample_index][phase];
                if allele == 0 {
                    sequence.extend_from_slice(record.reference.as_bytes());
                } else {
                    sequence.extend_from_slice(record.alternates[allele - 1].as_bytes());
                }
                cursor = end;
            }
            sequence.extend_from_slice(&reference[cursor..]);
            copies.insert((sample.clone(), phase), sequence);
        }
    }
    Ok(copies)
}

/// El VCF describe **exactamente** las copias de la verdad, sin que ningún
/// registro pise a otro.
fn assert_vcf_describes_every_copy(published: &Published, text: &str) {
    let reference = fasta_bases(&published.fasta_path());
    let expected = copies_from_truth(&published.truth_path(), &reference);
    let observed = copies_from_vcf(text, &reference).unwrap_or_else(|problem| panic!("{problem}"));
    assert!(!expected.is_empty(), "control: hay copias que comparar");
    assert_eq!(
        observed.len(),
        expected.len(),
        "una copia por fase y muestra"
    );
    for (copy, sequence) in &expected {
        assert!(
            observed.get(copy) == Some(sequence),
            "la copia {copy:?} que describe el VCF no es la de la verdad"
        );
    }
}

/// `bcftools norm` 1.24 deja intacto todo registro que no sea de repetición,
/// y normalizar no cambia ninguna copia reconstruida.
fn assert_bcftools_norm_keeps_the_records(published: &Published) {
    let version = Command::new("bcftools")
        .arg("--version")
        .output()
        .expect("bcftools instalado");
    assert_eq!(
        String::from_utf8_lossy(&version.stdout).lines().next(),
        Some("bcftools 1.24"),
        "DTC-066 exige exactamente esta versión"
    );
    let normalized = Command::new("bcftools")
        .args(["norm", "--no-version", "--check-ref", "e", "--fasta-ref"])
        .arg(published.fasta_path())
        .arg(published.vcf_path())
        .output()
        .expect("bcftools norm");
    assert!(
        normalized.status.success(),
        "bcftools norm falló: {}",
        String::from_utf8_lossy(&normalized.stderr)
    );
    let text = String::from_utf8(normalized.stdout).expect("UTF-8");
    let after = records(&text);
    let before = published.records();
    assert_eq!(
        after.len(),
        before.len(),
        "norm no reparte ni junta registros"
    );
    for record in before.iter().filter(|record| !record.is_repeat()) {
        assert!(
            after.iter().any(|other| other.line == record.line),
            "bcftools norm alteró un registro que no es de repetición: {}",
            record.line
        );
    }
    // Tras recortar un registro de repetición, `bcftools norm` puede dejarlo
    // encima de un registro vecino —lo alinea a la izquierda sobre él—, y
    // entonces el fichero normalizado ya no se reconstruye registro a
    // registro. En ese caso, y solo en él, no se compara.
    let reference = fasta_bases(&published.fasta_path());
    if let Ok(observed) = copies_from_vcf(&text, &reference) {
        let expected = copies_from_truth(&published.truth_path(), &reference);
        for (copy, sequence) in &expected {
            assert!(
                observed.get(copy) == Some(sequence),
                "tras bcftools norm la copia {copy:?} cambia"
            );
        }
    } else {
        assert!(
            before.iter().any(Record::is_repeat),
            "solo el recorte de un registro de repetición puede llevarlo encima de otro"
        );
    }
}

/// El juez semántico, el oráculo propio y `bcftools norm`, sobre lo publicado.
fn assert_published_vcf_is_exact(published: &Published) {
    validate_published_vcf_semantics(
        &published.fasta_path(),
        &published.truth_path(),
        &published.vcf_path(),
    )
    .unwrap_or_else(|error| {
        panic!(
            "el juez semántico rechazó el VCF: {error}\n{}",
            published.vcf()
        )
    });
    assert_vcf_describes_every_copy(published, &published.vcf());
    assert_annotations_describe_the_tracts(published);
    assert_bcftools_norm_keeps_the_records(published);
}

/// Un tracto que alguna copia realiza, según `truth.jsonl`: su intervalo, su
/// motivo y la anotación declarada de cada copia que lo realiza.
struct TruthTract {
    start: usize,
    end: usize,
    motif: String,
    declared: BTreeMap<(String, usize), Vec<u64>>,
}

fn tracts_from_truth(truth: &Path) -> Vec<TruthTract> {
    let mut tracts = BTreeMap::<(usize, usize), TruthTract>::new();
    for line in fs::read_to_string(truth).expect("truth").lines() {
        let record: Value = serde_json::from_str(line).expect("JSON line");
        if record["record_type"] != "haplotype_copy" {
            continue;
        }
        let copy = (
            record["identity"]["individual_id"]
                .as_str()
                .expect("individual")
                .to_owned(),
            usize::try_from(record["phase"].as_u64().expect("phase")).expect("phase"),
        );
        for piece in record["current_pieces"].as_array().expect("pieces") {
            if piece["kind"] != "repeat_allele" {
                continue;
            }
            let start = piece["source"]["start"].as_u64().expect("start") as usize;
            let end = piece["source"]["end"].as_u64().expect("end") as usize;
            let blocks = piece["blocks"].as_array().expect("blocks");
            let tract = tracts.entry((start, end)).or_insert_with(|| TruthTract {
                start,
                end,
                motif: blocks[0]["motif"].as_str().expect("motif").to_owned(),
                declared: BTreeMap::new(),
            });
            tract.declared.insert(
                copy.clone(),
                blocks
                    .iter()
                    .map(|block| block["bases"].as_u64().expect("bases"))
                    .collect(),
            );
        }
    }
    tracts.into_values().collect()
}

/// Lo que la regla de conflicto de DTC-064 deriva de unas bases, escrito aquí
/// por su cuenta: de la primera a la última copia completa del motivo, sin
/// contar lo que el `ALT` comparte con el `REF` a los lados del tracto.
fn derived_tract_bases(
    alternate: &str,
    reference: &str,
    before: usize,
    after: usize,
    motif: &str,
) -> Vec<u64> {
    let alternate = alternate.as_bytes();
    let reference = reference.as_bytes();
    let lead = alternate
        .iter()
        .zip(reference)
        .take_while(|(left, right)| left == right)
        .count()
        .min(before);
    let trail = alternate[lead..]
        .iter()
        .rev()
        .zip(reference[lead..].iter().rev())
        .take_while(|(left, right)| left == right)
        .count()
        .min(after);
    let middle = &alternate[lead..alternate.len() - trail];
    let unit = motif.as_bytes();
    let (mut first, mut last, mut at) = (None, 0, 0);
    while at + unit.len() <= middle.len() {
        if &middle[at..at + unit.len()] == unit {
            first.get_or_insert(at);
            at += unit.len();
            last = at;
        } else {
            at += 1;
        }
    }
    first.map_or_else(Vec::new, |first| vec![(last - first) as u64])
}

/// **Oráculo de anotación propio.** Cada `ALT` de un registro de repetición
/// anota, tracto a tracto, lo que la copia lleva en ese sitio: la anotación
/// declarada si la realiza, y todas las bases del tracto de la referencia si
/// lo conserva. Si dos copias con el mismo `ALT` merecerían anotaciones
/// distintas, se exige la derivada de las bases. Mira solo las bases
/// publicadas en `RB` y su reparto por `RN`, y no usa nada del programa.
fn assert_annotations_describe_the_tracts(published: &Published) {
    let tracts = tracts_from_truth(&published.truth_path());
    let text = published.vcf();
    let samples = text
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("cabecera")
        .split('\t')
        .skip(9)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    for record in records(&text).iter().filter(|record| record.is_repeat()) {
        let (start, end) = record.interval();
        let inside = tracts
            .iter()
            .filter(|tract| start <= tract.start && tract.end <= end)
            .collect::<Vec<_>>();
        assert!(
            !inside.is_empty(),
            "registro de repetición sin tracto: {}",
            record.line
        );
        let counts = info_list(record, "RN")
            .iter()
            .map(|value| value.parse::<usize>().expect("RN"))
            .collect::<Vec<_>>();
        let bases = info_list(record, "RB")
            .iter()
            .map(|value| value.parse::<u64>().expect("RB"))
            .collect::<Vec<_>>();
        let mut published_by_alt = Vec::new();
        let mut cursor = 0;
        for count in &counts {
            published_by_alt.push(bases[cursor..cursor + count].to_vec());
            cursor += count;
        }
        // Una anotación no puede contar más bases de repetición de las que el
        // `ALT` lleva escritas.
        for (alternate, annotation) in record.alternates.iter().zip(&published_by_alt) {
            assert!(
                annotation.iter().sum::<u64>() <= alternate.len() as u64,
                "la anotación cuenta más bases de las que el ALT lleva: {}",
                record.line
            );
        }
        let mut expected_by_alt = BTreeMap::<usize, BTreeSet<Vec<u64>>>::new();
        for (sample_index, sample) in samples.iter().enumerate() {
            for (phase, allele) in record.genotypes[sample_index].iter().enumerate() {
                if *allele == 0 {
                    continue;
                }
                let copy = (sample.clone(), phase);
                let expected = inside
                    .iter()
                    .flat_map(|tract| {
                        tract
                            .declared
                            .get(&copy)
                            .cloned()
                            .unwrap_or_else(|| vec![(tract.end - tract.start) as u64])
                    })
                    .collect::<Vec<_>>();
                expected_by_alt.entry(*allele).or_default().insert(expected);
            }
        }
        for (allele, expected) in expected_by_alt {
            let observed = &published_by_alt[allele - 1];
            if expected.len() == 1 {
                assert_eq!(
                    Some(observed),
                    expected.iter().next(),
                    "la anotación del ALT {allele} no es el tracto que llevan sus copias: {}",
                    record.line
                );
            } else {
                assert_eq!(
                    inside.len(),
                    1,
                    "conflicto con varios tractos: {}",
                    record.line
                );
                let tract = inside[0];
                let derived = derived_tract_bases(
                    &record.alternates[allele - 1],
                    &record.reference,
                    tract.start - start,
                    end - tract.end,
                    &tract.motif,
                );
                assert_eq!(
                    observed, &derived,
                    "dos copias con el mismo ALT {allele}: se exige la anotación derivada: {}",
                    record.line
                );
            }
        }
    }
}

fn records_between(published: &Published, from: usize, to: usize) -> Vec<Record> {
    published
        .records()
        .into_iter()
        .filter(|record| from <= record.position && record.position <= to)
        .collect()
}

fn judge(published: &Published, text: &str) -> Result<(), VcfSemanticValidationError> {
    let path = published.root.join("candidate.vcf");
    fs::write(&path, text).expect("candidate VCF");
    validate_published_vcf_semantics(&published.fasta_path(), &published.truth_path(), &path)
}

// ---------------------------------------------------------------------------
// El ejemplo de la revisión
// ---------------------------------------------------------------------------

/// `ind-2` borra `[1246, 1253)` en su fase 0; `ind-1` en las dos fases e
/// `ind-3` en la fase 0 cambian la C de la posición 1250 por G.
fn example(label: &str) -> Published {
    let bases = plain_reference();
    publish(
        label,
        &bases,
        "",
        &[
            deletion("ind-2", 0, 1246, 1253),
            snp_to("ind-1", 0, 1250, 'C', 'G'),
            snp_to("ind-1", 1, 1250, 'C', 'G'),
            snp_to("ind-3", 0, 1250, 'C', 'G'),
        ],
    )
}

/// **Un cambio de base dentro de lo que otra copia borra sale en el mismo
/// registro que la deleción**, y la copia con la deleción no aparece con la C
/// de referencia en una base que no tiene.
///
/// Antes salían dos registros, `1246 CCCTCTCT>C` y `1250 C>G`, con `ind-2` en
/// `0|0` en el segundo. El registro esperado está escrito a mano desde el
/// escenario: `ALT` ordenados, la deleción en `ind-2` fase 0 y el cambio de
/// base en las tres copias que lo llevan.
#[test]
fn a_snp_inside_another_copys_deletion_is_one_complex_record() {
    let published = example("example");
    let region = records_between(&published, 1200, 1300);
    assert_eq!(region.len(), 1, "un solo registro: {region:?}");
    let record = &region[0];
    assert_eq!(record.position, 1246);
    assert_eq!(record.reference, "CCCTCTCT");
    assert_eq!(record.alternates, ["C", "CCCTGTCT"]);
    assert_eq!(
        (record.gt(0), record.gt(1), record.gt(2)),
        ("2|2".to_owned(), "1|0".to_owned(), "2|0".to_owned())
    );
    assert_eq!(
        record
            .info
            .trim_start_matches("GPS_CASE=")
            .split(',')
            .count(),
        3,
        "GPS_CASE lleva la deleción y los dos casos de cambio de base"
    );
    assert_published_vcf_is_exact(&published);
}

/// **El juez rechaza los dos registros de antes y un registro fundido que dé
/// la referencia a la copia borrada.**
///
/// Tres corrupciones del VCF correcto, cada una con su porqué:
///
/// - **los dos registros de antes**, tal como los publicaba el productor. El
///   primero ya no coincide con el registro esperado, de modo que el juez lo
///   rechaza en la primera línea por `REF`/`ALT`; esta variante demuestra que el
///   juez no acepta la forma vieja, no qué guarda la caza;
/// - **saturada para la guarda de solapamiento:** el registro fundido correcto
///   y, detrás, el registro `1250 C>G` de antes. La primera línea es correcta,
///   de modo que lo único que puede rechazar la segunda antes de buscarla en lo
///   esperado es la guarda que impide que un registro pise a otro; se exige su
///   detalle;
/// - **la copia borrada con la referencia:** el registro fundido con `ind-2`
///   en `0|0`.
///
/// Control: el VCF publicado se acepta.
#[test]
fn the_oracle_rejects_the_former_separate_records_and_a_reference_call_on_the_deleted_copy() {
    let published = example("oracle");
    let text = published.vcf();
    judge(&published, &text).expect("control: el VCF publicado se acepta");
    let fused = records(&text)
        .into_iter()
        .find(|record| record.position == 1246)
        .expect("registro fundido");
    let cases = fused.info.trim_start_matches("GPS_CASE=").to_owned();
    let replace = |lines: &str| text.replace(&format!("{}\n", fused.line), lines);

    let former = replace(&format!(
        "chr1\t1246\t.\tCCCTCTCT\tC\t.\tPASS\tGPS_CASE={cases}\tGT\t0|0\t1|0\t0|0\n\
         chr1\t1250\t.\tC\tG\t.\tPASS\tGPS_CASE={cases}\tGT\t1|1\t0|0\t1|0\n"
    ));
    assert!(
        judge(&published, &former).is_err(),
        "la forma vieja se rechaza"
    );

    let overlapping = replace(&format!(
        "{}\nchr1\t1250\t.\tC\tG\t.\tPASS\tGPS_CASE={cases}\tGT\t1|1\t0|0\t1|0\n",
        fused.line
    ));
    match judge(&published, &overlapping) {
        Err(VcfSemanticValidationError::Vcf { detail, .. }) => assert!(
            detail.contains("pisa"),
            "tiene que rechazarlo la guarda de solapamiento: {detail}"
        ),
        other => panic!("dos registros solapados no pueden aceptarse: {other:?}"),
    }

    // Un registro que comparte **una sola base** con el anterior —su `POS` es
    // la última base del registro fundido— también lo pisa.
    let last_base = replace(&format!(
        "{}\nchr1\t1253\t.\tT\tG\t.\tPASS\tGPS_CASE={cases}\tGT\t1|1\t0|0\t1|0\n",
        fused.line
    ));
    match judge(&published, &last_base) {
        Err(VcfSemanticValidationError::Vcf { detail, .. }) => assert!(
            detail.contains("pisa"),
            "una sola base compartida también es solaparse: {detail}"
        ),
        other => panic!("un registro que pisa la última base no puede aceptarse: {other:?}"),
    }

    // El `GT` de una copia se lee para reconstruirla antes de compararlo con
    // lo esperado, y lo que no puede leerse se rechaza con su propio motivo.
    for (gt, expected) in [
        (
            "\t2|2\t1\t2|0",
            "GT no trae un alelo para cada fase de la copia",
        ),
        ("\t2|2\tx|0\t2|0", "GT trae un alelo que no es un índice"),
        (
            "\t2|2\t5|0\t2|0",
            "GT nombra un ALT que el registro no tiene",
        ),
    ] {
        let unreadable = replace(&format!("{}\n", fused.line.replace("\t2|2\t1|0\t2|0", gt)));
        assert_ne!(unreadable, text, "la mutación cambió el VCF");
        match judge(&published, &unreadable) {
            Err(VcfSemanticValidationError::Vcf { detail, .. }) => {
                assert_eq!(detail, expected, "{gt}");
            }
            other => panic!("{gt}: {other:?}"),
        }
    }

    let reference_call = replace(&format!(
        "{}\n",
        fused.line.replace("\t2|2\t1|0\t2|0", "\t2|2\t0|0\t2|0")
    ));
    assert_ne!(reference_call, text, "la mutación cambió el VCF");
    assert!(
        judge(&published, &reference_call).is_err(),
        "la copia borrada no puede llevar la referencia"
    );
}

// ---------------------------------------------------------------------------
// La definición de solapamiento
// ---------------------------------------------------------------------------

/// **Se funde por transitividad.** A borra `[1500, 1505)` en `ind-1`, B borra
/// `[1503, 1510)` en `ind-2` y C cambia la base 1509 en base uno en `ind-3`.
/// A pisa a B, B pisa a C, y A y C no se tocan: `REF` de A es `[1499, 1505)` y
/// C está en la base 1508 de base cero. Sale un solo registro.
#[test]
fn chained_overlaps_fuse_even_when_the_ends_do_not_touch() {
    let bases = plain_reference();
    let published = publish(
        "chain",
        &bases,
        "",
        &[
            deletion("ind-1", 0, 1500, 1505),
            deletion("ind-2", 0, 1503, 1510),
            snp("ind-3", 0, 1509, &bases),
        ],
    );
    let region = records_between(&published, 1400, 1600);
    assert_eq!(region.len(), 1, "A, B y C en un registro: {region:?}");
    assert_eq!(region[0].alternates.len(), 3, "un ALT por copia distinta");
    assert_published_vcf_is_exact(&published);
}

/// **Lo que solo es adyacente no se funde.** Dos cambios de base contiguos en
/// copias distintas, y la deleción `[1246, 1253)` —su `REF` es `CCCTCTCT`,
/// `[1245, 1253)`, y no se desplaza— con un cambio de base en la base 1254 en
/// base uno, la siguiente a lo borrado: cuatro registros.
#[test]
fn adjacent_variants_stay_separate() {
    let bases = plain_reference();
    let published = publish(
        "adjacent",
        &bases,
        "",
        &[
            snp("ind-1", 0, 1300, &bases),
            snp("ind-2", 0, 1301, &bases),
            deletion("ind-1", 1, 1246, 1253),
            snp("ind-3", 1, 1254, &bases),
        ],
    );
    let positions = records_between(&published, 1200, 1400)
        .iter()
        .map(|record| record.position)
        .collect::<Vec<_>>();
    assert_eq!(positions, [1246, 1254, 1300, 1301]);
    assert_published_vcf_is_exact(&published);
}

/// **Una inserción dentro de lo borrado, o en su borde derecho, se funde con
/// la deleción; una delante de su ancla, no.**
///
/// La deleción de `ind-1` borra `[1246, 1253)` y su `REF` es `CCCTCTCT`,
/// `[1245, 1253)`. La inserción es `GA`, que no se desplaza al alinearla en
/// ninguno de los tres puntos.
///
/// - Dentro, en el punto 1250, se funde.
/// - En el borde derecho, en el punto 1253, también: su ancla es la base 1252,
///   la última borrada. Dejarla aparte anotaría a la copia de la deleción con
///   una base que no tiene, que es el defecto que esta regla corrige.
/// - En el punto 1245 su ancla es la base 1244, fuera del `REF` de la
///   deleción: solo es adyacente y sale en su propio registro.
#[test]
fn insertions_inside_or_at_the_right_border_of_a_deletion_fuse_and_one_before_its_anchor_does_not()
{
    let bases = plain_reference();
    for (point, fused) in [(1250, true), (1253, true), (1245, false)] {
        let published = publish(
            &format!("insertion-{point}"),
            &bases,
            "",
            &[
                deletion("ind-1", 0, 1246, 1253),
                insertion("ind-2", 0, point, "GA"),
            ],
        );
        let region = records_between(&published, 1200, 1290);
        assert_eq!(
            region.len(),
            if fused { 1 } else { 2 },
            "inserción en {point}: {region:?}"
        );
        assert_published_vcf_is_exact(&published);
    }
}

/// **Un cambio de base y una inserción anclada en la misma base siguen siendo
/// el registro multialélico de siempre**, con el `REF` de una base.
#[test]
fn a_snp_and_an_insertion_on_the_same_base_stay_one_multiallelic_record() {
    let bases = plain_reference();
    let published = publish(
        "same-base",
        &bases,
        "",
        &[
            snp("ind-1", 0, 1800, &bases),
            insertion("ind-2", 1, 1800, "GA"),
        ],
    );
    let region = records_between(&published, 1700, 1900);
    assert_eq!(region.len(), 1, "{region:?}");
    assert_eq!(region[0].position, 1800);
    assert_eq!(region[0].reference.len(), 1);
    assert_eq!(region[0].alternates.len(), 2);
    assert_published_vcf_is_exact(&published);
}

/// **Sin solapamiento entre copias, el VCF es idéntico byte a byte al del
/// algoritmo anterior.**
///
/// El fichero de referencia lo publicó el binario compilado desde el commit
/// `61b9177`, el último antes de la fusión, sobre este mismo escenario: una
/// deleción con un cambio de base y una inserción en su misma base
/// —multialélico por `CHROM+POS`—, cambios de base contiguos y una deleción
/// con un cambio de base justo detrás. Lo único que la versión nueva cambia en
/// la ejecución es el campo de versión de la configuración efectiva, que no
/// está en el VCF.
#[test]
fn without_cross_copy_overlaps_the_vcf_is_byte_identical_to_the_previous_algorithm() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/vcf-fusion");
    let published = publish_from_dir(fixture.join("no-cross-copy-overlap"), "identity");
    let expected =
        fs::read(fixture.join("no-cross-copy-overlap/truth-61b9177.vcf")).expect("VCF de 61b9177");
    assert_eq!(fs::read(published.vcf_path()).expect("VCF"), expected);
    assert_published_vcf_is_exact(&published);
}

// ---------------------------------------------------------------------------
// Microsatélites
// ---------------------------------------------------------------------------

/// Un cambio en cada una de las seis copias.
fn on_every_copy(event: impl Fn(&str, u8) -> Event) -> Vec<Event> {
    SAMPLES
        .iter()
        .flat_map(|sample| [0, 1].map(|phase| event(sample, phase)))
        .collect()
}

/// El registro de repetición del locus `(CA)50`, que empieza en la 101.
fn tract_record(published: &Published) -> Record {
    let region = records_between(published, 90, 260);
    assert_eq!(
        region.len(),
        1,
        "todo cae en el registro del tracto: {region:?}"
    );
    let record = region[0].clone();
    assert!(record.is_repeat(), "{}", record.line);
    record
}

fn info_list<'a>(record: &'a Record, key: &str) -> Vec<&'a str> {
    record
        .info
        .split(';')
        .find_map(|field| field.strip_prefix(&format!("{key}=")))
        .map(|value| value.split(',').collect())
        .unwrap_or_default()
}

/// La anotación de cada `ALT` del registro del tracto: (bases del `ALT`, `RB`).
fn alternate_tract_bases(record: &Record) -> Vec<(usize, u64)> {
    let bases = info_list(record, "RB");
    assert_eq!(bases.len(), record.alternates.len(), "un bloque por ALT");
    record
        .alternates
        .iter()
        .zip(bases)
        .map(|(alternate, bases)| (alternate.len(), bases.parse().expect("RB")))
        .collect()
}

/// **Una inserción en el borde derecho del tracto, en una copia que conserva
/// el de la referencia, se funde en el registro de repetición.**
///
/// La inserción va en las seis copias, de modo que la llevan copias que
/// conservan el tracto de la referencia y copias de 45 unidades sea cual sea
/// el sorteo; el control lo exige. El `ALT` de las primeras es el tracto de la
/// referencia más `TT`, 102 bases, y su anotación es **la del tracto que
/// conservan, por posición** (DTC-064): un bloque de `CA` con sus 100 bases;
/// el `TT` no es repetición. Las copias de 45 unidades llevan la anotación
/// declarada, 90.
/// El registro no se recorta: empieza en la 101 y su `REF` es el tracto
/// entero.
#[test]
fn an_insertion_at_the_right_border_of_a_tract_fuses_into_the_repeat_record() {
    let bases = tract_reference(&"G".repeat(100));
    let published = publish(
        "tract-insertion",
        &bases,
        TRACT_LOCUS,
        &on_every_copy(|sample, phase| insertion(sample, phase, 200, "TT")),
    );
    let record = tract_record(&published);
    assert_eq!(record.position, 101);
    assert_eq!(record.reference, "CA".repeat(50));
    let annotations = alternate_tract_bases(&record);
    assert!(
        annotations.contains(&(102, 100)),
        "el ALT del tracto de la referencia más TT anota el tracto: {annotations:?}"
    );
    assert!(
        annotations.contains(&(92, 90)),
        "control: el alelo de 45 unidades con TT anota sus 90 bases declaradas: {annotations:?}"
    );
    assert_published_vcf_is_exact(&published);
}

/// **Una deleción de `CA` justo detrás del tracto que, alineada a la izquierda,
/// lo cruza entero, se funde con el registro de repetición y conserva la
/// representación de siempre.**
///
/// El tracto `(CA)50` es puro y el flanco sigue con `CAT`, de modo que la
/// deleción `[200, 202)` se alinea hasta la `G` que precede al tracto: lo cruza
/// entero, y por la decisión del 2 de octubre de 2026 se le deja cruzarlo. Su
/// `REF` `[99, 102)` pisa el del tracto y se funde: el registro empieza en la
/// 100 y cubre `[99, 202)`. El `ALT` de la copia que conserva el tracto de la
/// referencia es `G` más cincuenta `CA`, y su anotación es la del tracto que
/// conserva, por posición: 100 bases. Las copias de 45 unidades llevan su
/// anotación declarada, 90, y ejercitan además la ventana que arranca en el
/// tracto.
#[test]
fn a_ca_deletion_that_aligns_across_the_whole_tract_fuses_into_the_repeat_record() {
    let bases = tract_reference(&format!("CAT{}", "G".repeat(97)));
    let published = publish(
        "tract-deletion",
        &bases,
        TRACT_LOCUS,
        &on_every_copy(|sample, phase| deletion(sample, phase, 200, 202)),
    );
    let record = tract_record(&published);
    assert_eq!(record.position, 100);
    assert_eq!(record.reference, format!("G{}", "CA".repeat(51)));
    let annotations = alternate_tract_bases(&record);
    assert!(
        annotations.contains(&(101, 100)),
        "la copia que conserva el tracto: {annotations:?}"
    );
    assert!(
        annotations.contains(&(91, 90)),
        "control: una copia de 45 unidades con la deleción: {annotations:?}"
    );
    assert_published_vcf_is_exact(&published);
}

/// **Regresión: una deleción justo detrás de un tracto que la propia copia
/// realiza ya no aborta.**
///
/// Hasta el 2 de octubre de 2026 la ejecución terminaba con un error interno:
/// el ancla de la deleción era la última base del tracto de 45 unidades de esa
/// copia, que no es la de la referencia, y la ventana partía el tracto. El
/// escenario pone la deleción en todas las copias para que la lleve alguna
/// copia de 45 unidades sea cual sea el sorteo. El flanco es de `G`, de modo
/// que, a diferencia de la prueba anterior, la deleción no podría alinearse
/// dentro del tracto: lo único que se ejercita es el ancla.
#[test]
fn a_deletion_right_after_a_realized_tract_no_longer_aborts() {
    let bases = tract_reference(&"G".repeat(100));
    let events = on_every_copy(|sample, phase| deletion(sample, phase, 200, 202));
    let published = publish("tract-anchor", &bases, TRACT_LOCUS, &events);
    let record = tract_record(&published);
    assert_eq!(record.reference, format!("{}GG", "CA".repeat(50)));
    assert!(
        alternate_tract_bases(&record).contains(&(90, 90)),
        "control: alguna copia de 45 unidades lleva la deleción"
    );
    assert_published_vcf_is_exact(&published);
}

// ---------------------------------------------------------------------------
// La cota de bases de tracto (DTC-097)
// ---------------------------------------------------------------------------

fn tract_output(published: &Published) -> Value {
    published.json("manifest.json")["repeat_loci"]["tract_output"].clone()
}

/// Bases de tracto que el VCF publicado lleva, contadas aquí sobre el texto:
/// `REF` más cada `ALT` de los registros que anotan repeticiones.
fn counted_tract_bases(published: &Published) -> u64 {
    published
        .records()
        .iter()
        .filter(|record| record.is_repeat())
        .map(|record| {
            (record.reference.len() + record.alternates.iter().map(String::len).sum::<usize>())
                as u64
        })
        .sum()
}

/// **La parte VCF de la cota es lo que el VCF lleva** (DTC-097).
///
/// Desde el 2 de octubre de 2026 no se estima: se mide al publicar, antes de
/// escribir. Se exige que la cifra del manifiesto, la de la configuración
/// efectiva y la suma por locus sean las bases que esta prueba cuenta en el
/// fichero. Los escenarios son los que la estimación dejaba por debajo de lo
/// escrito —y la ejecución abortaba al validar el manifiesto—: la inserción de
/// `TT` en el borde, la deleción de `CA` detrás del tracto, un único alelo de
/// 45 unidades con `TT` —ese ya abortaba antes de la fusión—, una deleción de
/// cinco bases con un cambio de base en cada una de las cinco bases que borra
/// en las otras copias, y el contraejemplo de la revisión (`w-bound-0`), que
/// al alinear recorre el flanco que continúa el motivo y publicaba 246 bases
/// frente a una cota de 185. Y `w-cap-2090`, con un tope declarado de 2090 y
/// una deleción de 2300 bases, que la estimación más holgada rechazaba y la
/// base admitía: publica.
#[test]
fn the_published_vcf_tract_figure_is_what_the_vcf_carries() {
    let single_allele = TRACT_LOCUS.replace(
        "[[reference.repeat_loci.alleles]]\nunits = 50\nfrequency = 0.5\n\n",
        "",
    );
    let single_allele =
        single_allele.replace("units = 45\nfrequency = 0.5", "units = 45\nfrequency = 1.0");
    assert_ne!(
        single_allele, TRACT_LOCUS,
        "control: el locus de un alelo cambió"
    );
    let snps_flank = tract_reference(&format!("CAT{}", "G".repeat(97)));
    let scenarios = [
        (
            "bound-insertion",
            tract_reference(&"G".repeat(100)),
            TRACT_LOCUS.to_owned(),
            vec![insertion("ind-1", 0, 200, "TT")],
        ),
        (
            "bound-deletion",
            tract_reference(&format!("CAT{}", "G".repeat(97))),
            TRACT_LOCUS.to_owned(),
            vec![deletion("ind-1", 0, 200, 202)],
        ),
        (
            "bound-single-allele",
            tract_reference(&"G".repeat(100)),
            single_allele,
            vec![insertion("ind-1", 0, 200, "TT")],
        ),
        (
            "bound-snps-inside-a-deletion",
            snps_flank.clone(),
            TRACT_LOCUS.to_owned(),
            vec![
                deletion("ind-1", 0, 200, 205),
                snp("ind-1", 1, 201, &snps_flank),
                snp("ind-2", 0, 202, &snps_flank),
                snp("ind-2", 1, 203, &snps_flank),
                snp("ind-3", 0, 204, &snps_flank),
                snp("ind-3", 1, 205, &snps_flank),
            ],
        ),
    ];
    let mut published = scenarios
        .into_iter()
        .map(|(label, bases, locus, events)| publish(label, &bases, &locus, &events))
        .collect::<Vec<_>>();
    published.push(publish_fixture("w-bound-0"));
    published.push(publish_fixture("w-cap-2090"));
    for published in &published {
        let counted = counted_tract_bases(published);
        assert!(counted > 0, "control: hay registro de repetición");
        let output = tract_output(published);
        assert_eq!(output["vcf_tract_bases_bound"].as_u64(), Some(counted));
        assert_eq!(output["realized_vcf_tract_bases"].as_u64(), Some(counted));
        let budget = &published.json("effective-config.json")["repeat_output_budget"];
        assert_eq!(budget["vcf_tract_bases_bound"].as_u64(), Some(counted));
        let per_locus = budget["per_locus"]
            .as_array()
            .expect("por locus")
            .iter()
            .map(|locus| locus["vcf_tract_bases_bound"].as_u64().expect("entero"))
            .sum::<u64>();
        assert_eq!(per_locus, counted, "la suma por locus es la del VCF");
        assert_published_vcf_is_exact(published);
    }
}

/// **Sin inserciones ni deleciones directas la cifra también se mide, y por
/// eso cambia respecto a `61b9177`.**
///
/// El escenario tiene un locus `(CA)50` con alelos de 50 y 45 unidades y un
/// cambio de base directo. El binario de `61b9177` publicaba la estimación:
/// cien bases del tracto más cien y noventa de los dos alelos declarados, 290,
/// y 2090 en total. Pero el alelo de 50 unidades es el de la referencia y no
/// produce `ALT`, de modo que el VCF lleva 190: el tracto en `REF` y el alelo
/// de 45 unidades. Eso es lo que ahora se publica, y el total pasa a 1990. La
/// cota de la verdad, 1800, no cambia.
#[test]
fn without_direct_indels_the_vcf_tract_figure_is_measured_too() {
    let published = publish_from_dir(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/vcf-fusion/tract-without-indels"),
        "bound-measured",
    );
    let output = tract_output(&published);
    assert_eq!(output["truth_tract_bases_bound"], 1800);
    assert_eq!(output["vcf_tract_bases_bound"], 190);
    assert_eq!(counted_tract_bases(&published), 190);
    assert_ne!(
        output["vcf_tract_bases_bound"], 290,
        "ya no es la estimación de antes"
    );
    assert_eq!(output["total_tract_bases_bound"], 1990);
}

/// Copia un escenario de `tests/fixtures` a un directorio temporal y lo
/// publica.
fn publish_from_dir(fixture: PathBuf, label: &str) -> Published {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("root");
    for name in ["reference.fa", "snvs.tsv", "general.toml"] {
        fs::copy(fixture.join(name), root.join(name)).expect("fixture");
    }
    let out = root.join("published");
    run(&root.join("general.toml"), &out)
        .unwrap_or_else(|error| panic!("{label}: {}", error.render_human()));
    Published { root, out }
}

/// Un escenario de la revisión adversarial del 2 de octubre de 2026.
fn publish_fixture(name: &str) -> Published {
    publish_from_dir(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/vcf-fusion/review")
            .join(name),
        name,
    )
}

fn review_fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/vcf-fusion/review")
        .join(name)
}

// ---------------------------------------------------------------------------
// Los escenarios de la revisión adversarial del 2 de octubre de 2026
// ---------------------------------------------------------------------------

/// **Los escenarios que la revisión encontró rotos publican un VCF exacto.**
///
/// Cada uno pasa por el juez semántico, el oráculo propio de esta prueba —que
/// reconstruye cada copia desde el `GT`— y `bcftools norm`:
///
/// - `w-atreg-4` y `w-flankdel-5`: abortaban con un error interno, dos copias
///   con las mismas bases y anotaciones distintas;
/// - `w-plain-cross` e `y-2813`: la normalización de una copia alineaba una
///   inserción a través de un cambio de base de la propia copia y el `ALT`
///   fundido perdía el cambio;
/// - `w-min-a`: una inserción en el borde de una ventana se contaba dos veces y
///   la copia salía con una base de más, con el juez dándolo por bueno;
/// - `y-1699` y `x-11063`: el juez rechazaba un VCF exacto porque sus reglas
///   de recomposición diferían de las del productor;
/// - `r03`: una deleción de `CA` dos unidades detrás del tracto;
/// - `v-136`, `v-163` y `v-175`, de la familia realista con microsatélite:
///   el algoritmo anterior abortaba con dos copias que daban el mismo `ALT`
///   con anotaciones distintas, y ahora las dos llevan la derivada de sus
///   bases (decisión del 2 de octubre de 2026);
/// - `w-bound-0` y `w-cap-2090` están en la prueba de la cota.
#[test]
fn the_scenarios_found_by_the_review_publish_an_exact_vcf() {
    for name in [
        "w-atreg-4",
        "w-flankdel-5",
        "w-plain-cross",
        "y-2813",
        "w-min-a",
        "y-1699",
        "x-11063",
        "r03",
        "v-136",
        "v-163",
        "v-175",
        "x3338",
        "40045",
        "20248",
    ] {
        let published = publish_fixture(name);
        assert_published_vcf_is_exact(&published);
    }
}

/// **La copia que conserva el tracto de la referencia se anota por posición:
/// el tracto entero, nunca lo que sobra del `ALT`.**
///
/// Tres escenarios de la segunda revisión, con la cifra esperada escrita a
/// mano desde el escenario:
///
/// - `x3338`: un `(A)18` y una `G` insertada delante en una copia que lo
///   conserva. Su `ALT` lleva las 18 `A` y la `G`; la anotación es 18, no 19.
/// - `40045`: un `(CA)6` con un cambio de base en el flanco de otra copia; la
///   anotación es 12, no 17.
/// - `v-136`: un `(CA)15` con un cambio de base; la anotación es 30, no 32.
///   Ahí dos copias tienen las mismas bases y anotaciones distintas, y la
///   derivada cuenta solo lo que repite el motivo.
///
/// Además el oráculo de anotación de esta prueba repasa todos los registros.
#[test]
fn a_copy_that_keeps_the_reference_tract_is_annotated_by_position() {
    for (name, position, alternate_index, expected) in [
        ("x3338", 98, 2, 18_u64),
        ("40045", 151, 3, 12),
        ("v-136", 200, 3, 30),
    ] {
        let published = publish_fixture(name);
        let record = published
            .records()
            .into_iter()
            .find(|record| record.position == position)
            .unwrap_or_else(|| panic!("{name}: registro en {position}"));
        let annotations = alternate_tract_bases(&record);
        assert_eq!(
            annotations[alternate_index - 1].1,
            expected,
            "{name}: {}",
            record.line
        );
        assert_published_vcf_is_exact(&published);
    }
}

/// **El límite de alineación también vale para un componente que empieza justo
/// donde acaba el tracto.**
///
/// En `20248`, una deleción de `ind-2` fase 1 empieza en la base que sigue al
/// tracto `[176, 184)`; alineada libre acabaría dentro de él. DTC-064 manda
/// repetir la alineación con el límite en la última base del tracto también
/// para ese componente, porque está detrás del tracto aunque no le quede ni una
/// base de separación. El registro de repetición sale entonces en la 177 con
/// un `REF` de 12 bases, `CACACACAAAAA`. Si el componente que empieza en el
/// borde quedara fuera del límite, la deleción se quedaría dentro del tracto y
/// el registro mediría 11.
#[test]
fn a_component_starting_right_at_the_end_of_a_tract_is_also_limited() {
    let published = publish_fixture("20248");
    let record = published
        .records()
        .into_iter()
        .find(|record| record.position == 177)
        .expect("registro del segundo tracto");
    assert_eq!(record.reference, "CACACACAAAAA", "{}", record.line);
    assert!(record.is_repeat());
    assert_published_vcf_is_exact(&published);
}

/// **Una deleción que cruza el tracto al alinearse lleva el registro hasta el
/// final de su tramo original.**
///
/// En `rb-over-alt-1831` —de la exploración final del 2 de octubre de 2026—
/// el tracto `(CA)7` ocupa `[163, 177)` y la referencia sigue con `CACA`.
/// `ind-3` fase 2 conserva el tracto de la referencia y borra esas cuatro
/// bases, `[177, 181)`. Alineada libre, la deleción cruza el tracto entero y
/// queda anclada en la `C` de la 162, delante de él; DTC-064 deja que lo cruce.
/// Antes, el registro de repetición acababa con el tracto, en la 177, y el
/// `ALT` de esa copia era `C(CA)5`, 11 bases, anotado por posición con las 14
/// del tracto que conserva: más repetición de la que el `ALT` lleva. Ahora el
/// registro llega hasta la 181, el final del tramo original de la deleción,
/// y el `ALT` de esa copia es `C(CA)7`, 15 bases con sus 14 de repetición.
/// Sin el ensanche, el `REF` mediría 15 bases y el oráculo de anotación de
/// esta batería —que exige que ninguna anotación cuente más bases de las que
/// su `ALT` lleva— rechazaría el registro.
#[test]
fn a_deletion_that_crosses_a_tract_extends_the_record_to_its_original_end() {
    let published = publish_fixture("rb-over-alt-1831");
    let record = published
        .records()
        .into_iter()
        .find(|record| record.position == 163)
        .expect("registro del tracto (CA)7");
    assert!(record.is_repeat(), "{}", record.line);
    assert_eq!(
        record.reference,
        format!("C{}", "CA".repeat(9)),
        "el registro llega hasta el final de la deleción: {}",
        record.line
    );
    let samples = published
        .vcf()
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("cabecera")
        .split('\t')
        .skip(9)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let ind3 = samples
        .iter()
        .position(|sample| sample == "ind-3")
        .expect("ind-3");
    let allele = record.genotypes[ind3][2];
    assert!(
        allele > 0,
        "control: la copia que borra lleva un ALT: {}",
        record.line
    );
    assert_eq!(
        record.alternates[allele - 1],
        format!("C{}", "CA".repeat(7)),
        "{}",
        record.line
    );
    assert_eq!(
        alternate_tract_bases(&record)[allele - 1],
        (15, 14),
        "{}",
        record.line
    );
    assert_published_vcf_is_exact(&published);
}

/// **Donde el algoritmo anterior ya era exacto, su VCF vuelve byte a byte.**
///
/// Con el límite de alineación de la primera versión de la fusión, estos tres
/// escenarios cambiaban de representación: una inserción que cruza el tracto
/// entero quedaba anclada en su última base. Desde el 2 de octubre de 2026 el
/// límite solo actúa si la alineación acabaría dentro del tracto, y la salida
/// de `61b9177`, publicada por ese binario sobre los mismos ficheros, vuelve.
#[test]
fn where_the_previous_algorithm_was_exact_its_vcf_comes_back() {
    for name in ["x-1392", "x-12168", "x-14771"] {
        let published = publish_fixture(name);
        let expected =
            fs::read(review_fixture(name).join("truth-61b9177.vcf")).expect("VCF de 61b9177");
        assert_eq!(
            fs::read(published.vcf_path()).expect("VCF"),
            expected,
            "{name}"
        );
        assert_published_vcf_is_exact(&published);
    }
}

/// **El juez rechaza el VCF que publicaban las variantes defectuosas.**
///
/// Son los VCF que el binario de la primera versión de la fusión publicó sobre
/// `w-min-a` —una inserción contada dos veces— y `w-plain-cross` —un cambio
/// de base perdido—; el de `w-min-a` el juez de entonces lo daba por bueno. La
/// verdad es la de la ejecución de ahora, que no cambia: lo que cambió es el
/// VCF.
#[test]
fn the_oracle_rejects_the_vcf_of_the_defective_variants() {
    for name in ["w-min-a", "w-plain-cross"] {
        let published = publish_fixture(name);
        let defective =
            fs::read_to_string(review_fixture(name).join("truth-defective.vcf")).expect("VCF");
        assert_ne!(defective, published.vcf(), "control: {name} es otro VCF");
        assert!(
            judge(&published, &defective).is_err(),
            "{name}: el juez tiene que rechazar el VCF defectuoso"
        );
        judge(&published, &published.vcf()).expect("control: el VCF de ahora se acepta");
    }
}

/// **Una alineación que acabaría dentro de un tracto realizado se queda en su
/// última base y se funde con el registro de repetición.**
///
/// El tracto de la referencia es impuro —`(CA)2 T (CA)8`— y lo sigue un
/// flanco `CACA`. Solo `ind-3` fase 0, que conserva el tracto de la referencia,
/// borra el `CA` de `[123, 125)`, dos unidades detrás del tracto. Alineada
/// libre, la deleción se detendría tras la `T`, dentro del tracto; por eso se
/// alinea de nuevo con el límite en la última base del tracto, la 120, y su
/// `REF` `[120, 123)` se funde con el registro de repetición. El registro
/// cubre entonces `[100, 123)`: 23 bases. Sin el límite, la deleción quedaría
/// dentro del tracto y el registro mediría 21; con el límite una base más a la
/// derecha, mediría 24. El `ALT` de esa copia son las 21 bases del tracto, y
/// su anotación es la de ese tracto, por posición: 21.
#[test]
fn an_alignment_that_would_end_inside_a_realized_tract_stops_at_its_last_base() {
    let published = publish_fixture("impure-tract");
    let records = records_between(&published, 90, 140);
    assert_eq!(records.len(), 1, "{records:?}");
    let record = &records[0];
    assert_eq!(record.position, 101);
    assert_eq!(record.reference.len(), 23, "{}", record.line);
    assert!(
        alternate_tract_bases(record).contains(&(21, 21)),
        "{}",
        record.line
    );
    assert_published_vcf_is_exact(&published);
}

// ---------------------------------------------------------------------------
// El locus compuesto
// ---------------------------------------------------------------------------

/// Locus compuesto `(AAAT)10(GGC)9` con el alelo de la referencia entre los
/// declarados, de modo que hay copias que conservan el tracto de la
/// referencia, e inserciones de `TT` en su borde derecho en todas las copias.
fn compound_scenario(root: &Path, vcf: bool) -> PathBuf {
    let mut bases = b"C".repeat(100);
    bases.extend_from_slice("AAAT".repeat(10).as_bytes());
    bases.extend_from_slice("GGC".repeat(9).as_bytes());
    bases.extend_from_slice(&b"T".repeat(100));
    let locus = r#"
[[reference.repeat_loci]]
id = "STR-COMP"
sequence = "chr1"
start = 101
end = 167
motifs = ["AAAT", "GGC"]

[[reference.repeat_loci.alleles]]
units = [10, 9]
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = [10, 6]
frequency = 0.5
"#;
    let events = SAMPLES
        .iter()
        .flat_map(|sample| [0, 1].map(|phase| insertion(sample, phase, 167, "TT")))
        .collect::<Vec<_>>();
    write_scenario(root, &bases, locus, &events, vcf)
}

/// **Una variante que se funde con un microsatélite compuesto, en una copia
/// que conserva el tracto de la referencia, se rechaza con un diagnóstico de
/// configuración.**
///
/// La estructura en bloques de ese tracto no se declara y con varios motivos
/// no se deriva sin ambigüedad. El rechazo llega al publicar el VCF, después
/// de construir las copias y antes de simular lecturas, y no publica nada.
/// La misma verdad, publicada sin VCF, la rechaza también el juez semántico.
#[test]
fn a_variant_fused_with_a_compound_tract_kept_from_the_reference_is_refused() {
    let root = temporary("compound");
    let config = compound_scenario(&root, true);
    let error = run(&config, &root.join("published")).expect_err("no hay regla de bloques");
    assert_eq!(error.code(), "GPS-CONF-008", "{}", error.render_human());
    assert!(
        error.context()["detail"].contains("compuesto")
            && error.context()["detail"].contains("STR-COMP"),
        "el diagnóstico nombra el locus: {}",
        error.render_human()
    );
    assert_eq!(error.context()["locus_id"], "STR-COMP");
    assert!(
        !root.join("published/outputs/truth.vcf").exists(),
        "no se publica ningún VCF"
    );

    let truth_only = root.join("truth-only");
    let config = compound_scenario(&truth_only, false);
    let out = truth_only.join("published");
    run(&config, &out).expect("sin VCF la verdad se publica");
    let empty = truth_only.join("empty.vcf");
    fs::write(&empty, "##fileformat=VCFv4.5\n").expect("VCF");
    match validate_published_vcf_semantics(
        &out.join("resources/references/ref-main.fa"),
        &out.join("outputs/truth.jsonl"),
        &empty,
    ) {
        Err(VcfSemanticValidationError::Vcf { detail, .. }) => {
            assert!(detail.contains("compuesto"), "{detail}");
        }
        other => panic!("el juez tampoco tiene regla de bloques: {other:?}"),
    }
    fs::remove_dir_all(root).ok();
}

// ---------------------------------------------------------------------------
// La versión del algoritmo
// ---------------------------------------------------------------------------

/// **La configuración efectiva v6 y el manifiesto general v6 publican la
/// versión del algoritmo, y sus lectores rechazan otra o su ausencia.**
#[test]
fn the_normalization_version_is_published_and_any_other_is_refused() {
    let bases = tract_reference(&"G".repeat(100));
    let published = publish("version", &bases, TRACT_LOCUS, &[]);
    let effective_path = published.out.join("effective-config.json");
    let manifest_path = published.out.join("manifest.json");
    let effective = published.json("effective-config.json");
    let manifest = published.json("manifest.json");
    assert_eq!(effective["schema_version"], "6");
    assert_eq!(manifest["schema_version"], "6");
    assert_eq!(
        effective["vcf_output"]["normalization_version"],
        NORMALIZATION_VERSION
    );
    assert_eq!(
        manifest["vcf_outputs"][0]["normalization_version"],
        NORMALIZATION_VERSION
    );
    let manifest_bytes = serde_json::to_vec(&manifest).expect("JSON");
    validate_general_manifest_bytes(&manifest_path, &manifest_bytes)
        .expect("control: el manifiesto publicado se acepta");

    let detail = |error: Diagnostic| error.context().get("detail").cloned().unwrap_or_default();
    for tampered in [
        Value::from("vcf-small-allele-normalization-v1"),
        Value::Null,
    ] {
        let mut altered = manifest.clone();
        if tampered.is_null() {
            altered["vcf_outputs"][0]
                .as_object_mut()
                .expect("enlace")
                .remove("normalization_version");
        } else {
            altered["vcf_outputs"][0]["normalization_version"] = tampered.clone();
        }
        let error = validate_general_manifest_bytes(
            &manifest_path,
            &serde_json::to_vec(&altered).expect("JSON"),
        )
        .expect_err("el manifiesto no puede declarar otra versión ni callarla");
        assert!(
            detail(error).contains("normalization_version"),
            "{tampered}"
        );

        let mut altered = effective.clone();
        if tampered.is_null() {
            altered["vcf_output"]
                .as_object_mut()
                .expect("vcf_output")
                .remove("normalization_version");
        } else {
            altered["vcf_output"]["normalization_version"] = tampered.clone();
        }
        let altered = serde_json::to_vec(&altered).expect("JSON");
        // El lector de la configuración efectiva frente al manifiesto es el
        // que juzga la versión; la comprobación completa de lo publicado la
        // rechaza antes, al comparar con el TOML resuelto, y se exige igual.
        let error = validate_effective_config_against_manifest_bytes(
            &effective_path,
            &altered,
            &manifest_path,
            &manifest_bytes,
        )
        .expect_err("la configuración efectiva no puede declarar otra versión ni callarla");
        let detail = detail(error);
        // Callarla la rechaza el lector cerrado, que exige el campo; declarar
        // otra, la comprobación de la versión.
        let expected = if tampered.is_null() {
            "campos ausentes o desconocidos en vcf_output"
        } else {
            "normalization_version"
        };
        assert!(detail.contains(expected), "{tampered}: {detail}");
        validate_published_resources(
            &published.out,
            &effective_path,
            &altered,
            &manifest_path,
            &manifest_bytes,
        )
        .expect_err("ni la comprobación completa de lo publicado");
    }
}
