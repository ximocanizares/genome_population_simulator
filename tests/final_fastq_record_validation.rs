//! DT-026: el juez de un mate del FASTQ final tiene que leer sus registros.
//!
//! Hasta el 05/10/2026, `validate_final_fastq_mate` comprobaba el marcador de
//! fin BGZF con un descriptor clonado que **compartía la posición de lectura**
//! con el original. La comprobación dejaba esa posición al final del fichero y
//! el lector BGZF que venía después leía cero líneas y devolvía `Ok`: ninguna
//! comprobación de registro se ejecutaba en producción.
//!
//! Cada prueba de rechazo de registros de este fichero **falla contra aquel
//! código**, porque allí cualquier contenido con el marcador canónico se
//! aceptaba. Las del final —falta el marcador de fin, mates que nombran parejas
//! distintas, artefacto alterado— no cubren aquel defecto, sino comprobaciones
//! anteriores que la revisión adversarial encontró sin proteger: borrarlas
//! dejaba pasar la batería entera. Cada prueba de rechazo compara el
//! **motivo** del rechazo, no solo que lo haya: varias guardas se solapan (una
//! cabecera sin `@` tampoco empieza por `GPS`; un byte no UTF-8 tampoco es un
//! identificador válido), de modo que sin el motivo una guarda borrada quedaría
//! tapada por la siguiente y la prueba seguiría pasando.
//!
//! Lo que estas pruebas **no** pueden ejercitar, dicho sin adornos: dentro de
//! `validate_final_record`, las tres comprobaciones de que una línea termina
//! en salto de línea (cabecera, secuencia y calidades). `read_final_line`
//! rechaza antes toda línea final sin salto («FASTQ final line is truncated»)
//! y el bucle rechaza antes todo registro con una línea vacía, así que cuando
//! `validate_final_record` recibe las cuatro líneas, todas terminan ya en
//! `\n`. Son defensivas e inalcanzables desde un fichero.

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::fastq_validation::{
    FastqValidationError, FinalFastqCounts, FinalFastqIdentity, validate_final_fastq_mate,
    validate_final_paired_fastq,
};
use genome_population_simulator::integrity::identify_file;
use genome_population_simulator::results::ArtifactRecord;
use genome_population_simulator::technical_identity::TechnicalMate;
use noodles_bgzf as bgzf;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
const R1_HEADER: &str = "@GPS:reads-main:LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZ:TP_ABCDEFGHIJKLMNOPQRSTUVWXYY 1:N:0:ACGTACGT+TGCATGCA\n";
const R2_HEADER: &str = "@GPS:reads-main:LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZ:TP_ABCDEFGHIJKLMNOPQRSTUVWXYY 2:N:0:ACGTACGT+TGCATGCA\n";
const LIBRARY: &str = "LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const PAIR_A: &str = "TP_ABCDEFGHIJKLMNOPQRSTUVWXYY";
const PAIR_B: &str = "TP_ZYXWVUTSRQPONMLKJIHGFEDCBA";
static NONCE: AtomicU64 = AtomicU64::new(0);

/// Identidad que esperan casi todas las pruebas: los registros de este fichero
/// miden cuatro bases salvo que la prueba diga otra cosa.
fn identity() -> FinalFastqIdentity {
    FinalFastqIdentity {
        event_id: "reads-main".to_owned(),
        library_id: LIBRARY.to_owned(),
        udi: "ACGTACGT+TGCATGCA".to_owned(),
        read_length: 4,
    }
}

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-final-fastq-record-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Escribe `content` como BGZF con el marcador de fin canónico, que es lo
/// que el escritor de `noodles-bgzf` añade al cerrar. Así el fichero supera
/// la primera comprobación del juez y lo único que puede rechazarlo es su
/// contenido.
fn bgzf_file(content: &[u8]) -> PathBuf {
    let path = temporary("fastq.gz");
    let mut writer = bgzf::io::Writer::new(fs::File::create(&path).expect("fixture"));
    writer.write_all(content).expect("BGZF content");
    writer.finish().expect("BGZF EOF marker");
    path
}

fn record(header: &str, sequence: &[u8], separator: &[u8], qualities: &[u8]) -> Vec<u8> {
    let mut bytes = header.as_bytes().to_vec();
    bytes.extend_from_slice(sequence);
    bytes.push(b'\n');
    bytes.extend_from_slice(separator);
    bytes.push(b'\n');
    bytes.extend_from_slice(qualities);
    bytes.push(b'\n');
    bytes
}

fn valid_record() -> Vec<u8> {
    record(R1_HEADER, b"ACNT", b"+", b"II!~")
}

/// Motivo estructural del rechazo de `content` como mate `mate`.
fn rejection(content: &[u8], mate: TechnicalMate) -> String {
    let path = bgzf_file(content);
    let result = validate_final_fastq_mate(&path, mate, &identity());
    fs::remove_file(&path).ok();
    match result {
        Err(FastqValidationError::Structural { detail, .. }) => detail.to_owned(),
        Err(other) => panic!("unexpected non-structural rejection: {other:?}"),
        Ok(_) => panic!(
            "the final FASTQ judge accepted invalid content:\n{}",
            String::from_utf8_lossy(content)
        ),
    }
}

fn assert_rejected(content: &[u8], detail: &str) {
    assert_eq!(rejection(content, TechnicalMate::R1), detail);
}

/// El caso de la reproducción de DT-026, tal cual.
#[test]
fn garbage_with_canonical_bgzf_eof_is_rejected() {
    assert_rejected(b"hello\nworld\n+\nxx\n", "FASTQ final header lacks @");
}

/// Control positivo. Sin él, un juez que lo rechazara todo pasaría todas las
/// pruebas de rechazo. Cubre además los dos extremos del alfabeto de
/// calidades (`!` es Phred 0 y `~` Phred 93, el máximo del contrato) y la `N`.
#[test]
fn canonical_records_are_accepted_for_both_mates() {
    let mut r1 = valid_record();
    r1.extend_from_slice(&record(
        &R1_HEADER.replace(PAIR_A, PAIR_B),
        b"TNNT",
        b"+",
        b"~!!~",
    ));
    let path = bgzf_file(&r1);
    let tally = validate_final_fastq_mate(&path, TechnicalMate::R1, &identity()).expect("R1");
    fs::remove_file(&path).ok();
    assert_eq!((tally.records, tally.bases, tally.n_bases), (2, 8, 3));
    let path = bgzf_file(&record(R2_HEADER, b"ACNT", b"+", b"II!~"));
    validate_final_fastq_mate(&path, TechnicalMate::R2, &identity()).expect("canonical R2");
    fs::remove_file(&path).ok();
}

#[test]
fn header_that_is_not_utf8_is_rejected() {
    let mut bytes = b"@GPS:reads-\xff".to_vec();
    bytes.extend_from_slice(&valid_record()[12..]);
    assert_rejected(&bytes, "FASTQ final header is not UTF-8");
}

#[test]
fn header_with_identity_or_mate_errors_is_rejected() {
    // Cabecera de R2 dentro del fichero de R1.
    assert_eq!(
        rejection(
            &record(R2_HEADER, b"ACGT", b"+", b"IIII"),
            TechnicalMate::R1
        ),
        "FASTQ final mate metadata is invalid"
    );
    let header = R1_HEADER.replace("ACGTACGT+", "ACGTACGN+");
    assert_rejected(
        &record(&header, b"ACGT", b"+", b"IIII"),
        "FASTQ final UDI is invalid",
    );
    let header = R1_HEADER.replace("@GPS:", "@SIM:");
    assert_rejected(
        &record(&header, b"ACGT", b"+", b"IIII"),
        "FASTQ final header lacks GPS prefix",
    );
}

/// Cada rechazo del analizador de la cabecera (REQ-DAT-68), con su motivo.
///
/// «La cabecera no tiene identidad» no figura porque no se alcanza: partir un
/// texto por espacios siempre da un primer trozo, aunque sea vacío, y ese
/// trozo vacío lo rechaza la comprobación del prefijo `GPS`.
#[test]
fn every_header_rejection_names_its_own_reason() {
    const MATE: &str = " 1:N:0:ACGTACGT+TGCATGCA\n";
    const PAIR: &str = "TP_ABCDEFGHIJKLMNOPQRSTUVWXYY";
    for (header, detail) in [
        (
            format!("@GPS:reads-main:LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZ:{PAIR}\n"),
            "FASTQ final header lacks mate metadata",
        ),
        (
            R1_HEADER.replace('\n', " extra\n"),
            "FASTQ final header has extra fields",
        ),
        (format!("@GPS{MATE}"), "FASTQ final header lacks event ID"),
        (
            format!("@GPS:reads-main{MATE}"),
            "FASTQ final header lacks library ID",
        ),
        (
            format!("@GPS:reads-main:LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZ{MATE}"),
            "FASTQ final header lacks technical pair ID",
        ),
        (
            format!("@GPS:reads-main:LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZ:{PAIR}:extra{MATE}"),
            "FASTQ final header has invalid identity fields",
        ),
        (
            format!("@GPS::LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZ:{PAIR}{MATE}"),
            "FASTQ final event ID is invalid",
        ),
        (
            format!("@GPS:reads-main:-library:{PAIR}{MATE}"),
            "FASTQ final library ID is invalid",
        ),
        (
            "@GPS:reads-main:LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZ:_TP 1:N:0:ACGTACGT+TGCATGCA\n"
                .to_owned(),
            "FASTQ final pair ID is invalid",
        ),
        (
            R1_HEADER.replace("+TGCATGCA", "+TGCATGCA+ACGTACGT"),
            "FASTQ final UDI is invalid",
        ),
    ] {
        assert_eq!(
            rejection(&record(&header, b"ACGT", b"+", b"IIII"), TechnicalMate::R1),
            detail,
            "{header}"
        );
    }
}

#[test]
fn separator_other_than_plus_is_rejected() {
    for separator in [&b"-"[..], b"+extra", b""] {
        assert_rejected(
            &record(R1_HEADER, b"ACGT", separator, b"IIII"),
            "FASTQ final record structure is invalid",
        );
    }
}

#[test]
fn empty_sequence_is_rejected() {
    assert_rejected(
        &record(R1_HEADER, b"", b"+", b""),
        "FASTQ final record structure is invalid",
    );
}

#[test]
fn quality_length_different_from_sequence_is_rejected() {
    for qualities in [&b"III"[..], b"IIIII"] {
        assert_rejected(
            &record(R1_HEADER, b"ACGT", b"+", qualities),
            "FASTQ final record structure is invalid",
        );
    }
}

/// Minúsculas, un código IUPAC y un símbolo de hueco: ninguno pertenece al
/// alfabeto del FASTQ final, que es `ACGTN` en mayúsculas porque la
/// referencia se normaliza a mayúsculas (DTC-041) y el error de secuenciación
/// solo sustituye entre `ACGT`.
#[test]
fn bases_outside_the_final_alphabet_are_rejected() {
    for sequence in [&b"ACgT"[..], b"ACRT", b"AC-T", b"AC T"] {
        assert_rejected(
            &record(R1_HEADER, sequence, b"+", b"IIII"),
            "FASTQ final sequence or quality alphabet is invalid",
        );
    }
}

/// Justo por debajo de `!` (el espacio, Phred -1) y justo por encima de `~`
/// (DEL, Phred 94): los dos vecinos inmediatos del intervalo `0..=93`.
#[test]
fn qualities_outside_phred_0_to_93_are_rejected() {
    for qualities in [&b"II I"[..], b"II\x7fI"] {
        assert_rejected(
            &record(R1_HEADER, b"ACGT", b"+", qualities),
            "FASTQ final sequence or quality alphabet is invalid",
        );
    }
}

#[test]
fn trailing_incomplete_record_is_rejected() {
    let mut bytes = valid_record();
    bytes.extend_from_slice(R1_HEADER.as_bytes());
    bytes.extend_from_slice(b"ACGT\n");
    assert_rejected(&bytes, "FASTQ final has a truncated record");
}

#[test]
fn final_line_without_newline_is_rejected() {
    let mut bytes = valid_record();
    bytes.pop();
    assert_rejected(&bytes, "FASTQ final line is truncated");
}

#[test]
fn line_beyond_the_profile_bound_is_rejected() {
    let sequence = vec![b'A'; 300];
    assert_rejected(
        &record(R1_HEADER, &sequence, b"+", &vec![b'I'; 300]),
        "FASTQ final line exceeds the profile bound",
    );
}

/// Publica un escenario real con `N` en la referencia y errores de
/// secuenciación, y devuelve la raíz y los dos mates publicados.
fn published_mates() -> (PathBuf, PathBuf, PathBuf) {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    // Una tirada de `N` dentro del alcance de las lecturas: produce bases `N`
    // con Phred 0, el extremo inferior del alfabeto de calidades.
    let bases = format!(
        "{}{}{}",
        "ACGTTGCAAGCTTCGA".repeat(50),
        "N".repeat(40),
        "GATCCTAGGCATTACG".repeat(50)
    );
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    fs::write(
        root.join("population-snvs.tsv"),
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tT\tG\tP1\t0.5\n",
    )
    .expect("matriz");
    let length = bases.len();
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "supplied"
input_path = "population-snvs.tsv"

[vcf_output]
mode = "text"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {length}

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 5
sequencing_errors = "enabled"

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        ),
    )
    .expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    let directory = destination.join("outputs/fastq/reads-main");
    (
        root,
        directory.join("sample-A_R1.fastq.gz"),
        directory.join("sample-A_R2.fastq.gz"),
    )
}

/// Identidad de lo publicado: la biblioteca sale del manifiesto del evento y el
/// UDI y la longitud, de la configuración que se ejecutó, no de los ficheros.
fn published_manifest_entry(root: &Path) -> serde_json::Value {
    let manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(root.join("published/events/reads-main/manifest.json")).expect("manifiesto"),
    )
    .expect("JSON");
    manifest["final_fastq_summary"][0].clone()
}

fn published_identity(root: &Path) -> FinalFastqIdentity {
    FinalFastqIdentity {
        event_id: "reads-main".to_owned(),
        library_id: published_manifest_entry(root)["library_id"]
            .as_str()
            .expect("biblioteca")
            .to_owned(),
        udi: "ACGTACGT+TGCATGCA".to_owned(),
        read_length: 150,
    }
}

fn decompressed(path: &Path) -> Vec<u8> {
    let mut bytes = Vec::new();
    bgzf::io::Reader::new(fs::File::open(path).expect("mate publicado"))
        .read_to_end(&mut bytes)
        .expect("BGZF publicado");
    bytes
}

/// Un FASTQ publicado de verdad se acepta, y la misma copia con un solo
/// byte de calidad cambiado se rechaza.
///
/// La segunda mitad es la que hace que la primera signifique algo: contra el
/// código de DT-026 la copia alterada también se aceptaba. Las comprobaciones
/// del fixture sobre `N` y `!` no son decorativas: sin ellas, un cambio del
/// fixture que dejara de producir `N` dejaría sin ejercitar el extremo
/// inferior del alfabeto sin que nadie se enterase.
#[test]
fn real_published_mates_are_accepted_and_a_corrupted_copy_is_rejected() {
    let (root, r1, r2) = published_mates();
    let published = published_identity(&root);
    validate_final_fastq_mate(&r1, TechnicalMate::R1, &published).expect("published R1");
    validate_final_fastq_mate(&r2, TechnicalMate::R2, &published).expect("published R2");

    // Los recuentos que publica el manifiesto del evento son los del contenido
    // real. Este escenario sí tiene `N` (la referencia lleva una tirada y las
    // lecturas la cruzan): con un fixture sin `N`, un manifiesto que publicara
    // cero `N` pasaría.
    let manifest = published_manifest_entry(&root);
    let published_counts = FinalFastqCounts {
        pair_count: manifest["pair_count"].as_u64().expect("parejas"),
        base_count: manifest["base_count"].as_u64().expect("bases"),
        n_base_count: manifest["n_base_count"].as_u64().expect("N"),
    };
    assert!(
        published_counts.n_base_count > 0,
        "the fixture must publish N bases"
    );
    validate_final_paired_fastq(
        &r1,
        &artifact_of(&r1, "outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
        &r2,
        &artifact_of(&r2, "outputs/fastq/reads-main/sample-A_R2.fastq.gz"),
        &published,
        &published_counts,
    )
    .expect("the published manifest counts are those of the files");
    for altered in [
        FinalFastqCounts {
            pair_count: published_counts.pair_count + 1,
            ..published_counts
        },
        FinalFastqCounts {
            base_count: published_counts.base_count + 1,
            ..published_counts
        },
        FinalFastqCounts {
            n_base_count: 0,
            ..published_counts
        },
    ] {
        assert!(
            validate_final_paired_fastq(
                &r1,
                &artifact_of(&r1, "outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
                &r2,
                &artifact_of(&r2, "outputs/fastq/reads-main/sample-A_R2.fastq.gz"),
                &published,
                &altered,
            )
            .is_err(),
            "a summary that differs from the files must be rejected: {altered:?}"
        );
    }

    let content = decompressed(&r1);
    let lines: Vec<&[u8]> = content.split(|byte| *byte == b'\n').collect();
    assert!(
        lines.len() > 40,
        "the published corpus must be representative"
    );
    assert!(
        lines
            .iter()
            .skip(1)
            .step_by(4)
            .any(|line| line.contains(&b'N')),
        "the fixture must publish template N bases"
    );
    assert!(
        lines
            .iter()
            .skip(3)
            .step_by(4)
            .any(|line| line.contains(&b'!')),
        "the fixture must publish Phred 0 qualities"
    );

    // Primer byte de la primera línea de calidades, desplazado fuera del
    // intervalo Sanger. Todo lo demás, incluida la estructura, queda igual.
    let quality_offset = lines[..3].iter().map(|line| line.len() + 1).sum::<usize>();
    let mut corrupted = content.clone();
    corrupted[quality_offset] = b' ';
    let path = bgzf_file(&corrupted);
    match validate_final_fastq_mate(&path, TechnicalMate::R1, &published) {
        Err(FastqValidationError::Structural { detail, .. }) => assert_eq!(
            detail,
            "FASTQ final sequence or quality alphabet is invalid"
        ),
        other => panic!("expected the alphabet rejection, got {other:?}"),
    }
    fs::remove_file(&path).ok();
    fs::remove_dir_all(root).ok();
}

/// Los registros válidos sin el marcador de fin BGZF se rechazan, y por ese
/// motivo. Es la única comprobación del marcador en el camino del FASTQ final:
/// el lector BGZF acepta un final truncado como fin de flujo, así que sin ella
/// un fichero cortado justo después de un bloque se daría por completo.
#[test]
fn valid_records_without_the_canonical_eof_marker_are_rejected() {
    let path = bgzf_file(&valid_record());
    let bytes = fs::read(&path).expect("fixture");
    let marker = 28;
    fs::write(&path, &bytes[..bytes.len() - marker]).expect("strip EOF marker");
    let result = validate_final_fastq_mate(&path, TechnicalMate::R1, &identity());
    fs::remove_file(&path).ok();
    match result {
        Err(FastqValidationError::Structural { detail, .. }) => assert_eq!(
            detail,
            "FASTQ BGZF stream is missing its canonical EOF marker"
        ),
        other => panic!("expected the missing-EOF rejection, got {other:?}"),
    }
}

/// Artefacto publicado que describe exactamente los bytes de `path`.
fn artifact_of(path: &Path, relative: &str) -> ArtifactRecord {
    let identity = identify_file(path, relative).expect("identity");
    ArtifactRecord::from_identity("fastq-test", "fastq", &identity).expect("artifact record")
}

/// Una pareja sintética: los mates válidos por separado, con los artefactos
/// que los describen.
struct Pair {
    r1: PathBuf,
    r2: PathBuf,
    r1_artifact: ArtifactRecord,
    r2_artifact: ArtifactRecord,
    /// Recuentos de una pareja de un solo registro, contados aquí a mano.
    counts: FinalFastqCounts,
}

/// Recuento independiente de un fichero de un solo registro.
fn counts_of(r1: &[u8], r2: &[u8]) -> FinalFastqCounts {
    let sequence = |content: &[u8]| {
        content
            .split(|byte| *byte == b'\n')
            .nth(1)
            .unwrap()
            .to_vec()
    };
    let (a, b) = (sequence(r1), sequence(r2));
    FinalFastqCounts {
        pair_count: 1,
        base_count: u64::try_from(a.len() + b.len()).unwrap(),
        n_base_count: u64::try_from(a.iter().chain(&b).filter(|base| **base == b'N').count())
            .unwrap(),
    }
}

impl Pair {
    fn new(r1: &[u8], r2: &[u8]) -> Self {
        let counts = counts_of(r1, r2);
        let r1 = bgzf_file(r1);
        let r2 = bgzf_file(r2);
        let r1_artifact = artifact_of(&r1, "outputs/fastq/reads-main/sample-A_R1.fastq.gz");
        let r2_artifact = artifact_of(&r2, "outputs/fastq/reads-main/sample-A_R2.fastq.gz");
        Self {
            r1,
            r2,
            r1_artifact,
            r2_artifact,
            counts,
        }
    }

    fn judge(&self) -> Result<(), FastqValidationError> {
        self.judge_expecting(&identity(), &self.counts)
    }

    fn judge_expecting(
        &self,
        identity: &FinalFastqIdentity,
        counts: &FinalFastqCounts,
    ) -> Result<(), FastqValidationError> {
        validate_final_paired_fastq(
            &self.r1,
            &self.r1_artifact,
            &self.r2,
            &self.r2_artifact,
            identity,
            counts,
        )
    }
}

impl Drop for Pair {
    fn drop(&mut self) {
        fs::remove_file(&self.r1).ok();
        fs::remove_file(&self.r2).ok();
    }
}

fn canonical_pair() -> Pair {
    Pair::new(
        &record(R1_HEADER, b"ACNT", b"+", b"II!~"),
        &record(R2_HEADER, b"TTGC", b"+", b"IIII"),
    )
}

/// Control positivo del juez de la pareja: sin él, uno que lo rechazara todo
/// pasaría las tres pruebas de rechazo que siguen.
#[test]
fn a_canonical_pair_with_its_artifacts_is_accepted() {
    canonical_pair().judge().expect("canonical pair");
}

/// Cada mate es válido por separado, pero R2 nombra otra pareja. Solo la
/// comparación de cabeceras del juez de la pareja puede verlo.
#[test]
fn mates_naming_different_pairs_are_rejected() {
    let other = R2_HEADER.replace(
        "TP_ABCDEFGHIJKLMNOPQRSTUVWXYY",
        "TP_ZYXWVUTSRQPONMLKJIHGFEDCBA",
    );
    let pair = Pair::new(
        &record(R1_HEADER, b"ACGT", b"+", b"IIII"),
        &record(&other, b"ACGT", b"+", b"IIII"),
    );
    match pair.judge() {
        Err(FastqValidationError::Structural { detail, .. }) => assert_eq!(
            detail,
            "FASTQ final mates are orphaned or have different pair IDs"
        ),
        other => panic!("expected the pair-ID rejection, got {other:?}"),
    }
}

/// Un artefacto cuyo SHA-256 o tamaño no son los de los bytes publicados se
/// rechaza, y por el campo que no coincide.
#[test]
fn an_artifact_that_does_not_describe_the_published_bytes_is_rejected() {
    let mut pair = canonical_pair();
    let original = pair.r2_artifact.digest.clone();
    let last = original.chars().last().expect("digest");
    let flipped = if last == '0' { '1' } else { '0' };
    pair.r2_artifact.digest = format!("{}{flipped}", &original[..original.len() - 1]);
    match pair.judge() {
        Err(FastqValidationError::Artifact { field, .. }) => assert_eq!(field, "sha256"),
        other => panic!("expected the SHA-256 rejection, got {other:?}"),
    }
    pair.r2_artifact.digest = original;
    pair.r2_artifact.size_bytes += 1;
    match pair.judge() {
        Err(FastqValidationError::Artifact { field, .. }) => assert_eq!(field, "size_bytes"),
        other => panic!("expected the size rejection, got {other:?}"),
    }
}

/// El juez de la pareja juzga también los registros de cada mate. Desde el
/// 05/10/2026 es la única pasada que los juzga al publicar: el escritor ya no
/// juzga cada mate por su cuenta antes de catalogarlo.
///
/// Aquí las cabeceras coinciden y los artefactos describen los bytes, de modo
/// que lo único que puede rechazar la pareja es el registro de R1.
#[test]
fn the_pair_judge_judges_the_records_of_each_mate() {
    let pair = Pair::new(
        &record(R1_HEADER, b"ACGT", b"+", b"II I"),
        &record(R2_HEADER, b"ACGT", b"+", b"IIII"),
    );
    match pair.judge() {
        Err(FastqValidationError::Structural { detail, .. }) => assert_eq!(
            detail,
            "FASTQ final sequence or quality alphabet is invalid"
        ),
        other => panic!("expected the record rejection, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// DT-029: lo que el juez exige además de la forma del registro. Cada prueba
// compara el motivo exacto, porque varias comprobaciones se solapan.
// ---------------------------------------------------------------------------

/// Motivo de rechazo de `content` como R1 con una identidad esperada dada.
fn rejection_expecting(content: &[u8], identity: &FinalFastqIdentity) -> String {
    let path = bgzf_file(content);
    let result = validate_final_fastq_mate(&path, TechnicalMate::R1, identity);
    fs::remove_file(&path).ok();
    match result {
        Err(FastqValidationError::Structural { detail, .. }) => detail.to_owned(),
        other => panic!("expected a structural rejection, got {other:?}"),
    }
}

/// Un fichero que mezcla registros de otro evento, otra biblioteca u otro UDI
/// se rechaza, cada cosa por su motivo (hueco 1; casos `undecl_mixed_*` de la
/// revisión). El registro ajeno va en segunda posición: el primero es el
/// correcto, así que solo la comparación con la identidad esperada, registro a
/// registro, puede verlo.
#[test]
fn a_file_mixing_events_libraries_or_udis_is_rejected() {
    for (foreign, detail) in [
        (
            R1_HEADER.replace("reads-main", "reads-other"),
            "FASTQ final record belongs to a different event than expected",
        ),
        (
            R1_HEADER.replace(LIBRARY, "LIB_ZYXWVUTSRQPONMLKJIHGFEDCBA"),
            "FASTQ final record belongs to a different library than expected",
        ),
        (
            R1_HEADER.replace("ACGTACGT+TGCATGCA", "CGTACGTA+GCATGCAT"),
            "FASTQ final record carries a different UDI than its library",
        ),
        // Mismo i7 y otro i5: la forma de un salto de índice real, que cambia
        // un solo índice. Comparar solo el i7 lo dejaría pasar.
        (
            R1_HEADER.replace("ACGTACGT+TGCATGCA", "ACGTACGT+GCATGCAT"),
            "FASTQ final record carries a different UDI than its library",
        ),
    ] {
        let mut content = valid_record();
        content.extend_from_slice(&record(
            &foreign.replace(PAIR_A, PAIR_B),
            b"ACGT",
            b"+",
            b"IIII",
        ));
        assert_eq!(rejection_expecting(&content, &identity()), detail);
    }
}

/// Hueco 2 (`undecl_duplicate_pair_id`): el mismo identificador de pareja dos
/// veces en un fichero.
#[test]
fn a_repeated_pair_id_within_one_file_is_rejected() {
    let mut content = valid_record();
    content.extend_from_slice(&valid_record());
    assert_eq!(
        rejection_expecting(&content, &identity()),
        "FASTQ final repeats a technical pair ID within the file"
    );
    // El repetido no es contiguo: A, B, A. La unicidad se comprueba ordenando
    // los identificadores al final, y sin ordenar este fichero pasaría.
    let mut separated = valid_record();
    separated.extend_from_slice(&record(
        &R1_HEADER.replace(PAIR_A, PAIR_B),
        b"ACNT",
        b"+",
        b"II!~",
    ));
    separated.extend_from_slice(&valid_record());
    assert_eq!(
        rejection_expecting(&separated, &identity()),
        "FASTQ final repeats a technical pair ID within the file"
    );
}

/// Hueco 3: la pareja es `TP_` y 26 símbolos Base32; la biblioteca, `LIB_` y
/// 26. Los identificadores de abajo son todos válidos como identificadores
/// genéricos, que es lo único que se exigía antes (las seis corrupciones del
/// corpus de DT-026 cambiaban `TP_` por `XX_`).
#[test]
fn pair_and_library_ids_must_have_the_derived_shape() {
    for (header, detail) in [
        (
            R1_HEADER.replace(PAIR_A, "XX_ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
            "FASTQ final pair ID is not TP_ followed by 26 canonical Base32 symbols",
        ),
        (
            R1_HEADER.replace(PAIR_A, "TP_ABCDEFGHIJKLMNOPQRSTUVWXY"),
            "FASTQ final pair ID is not TP_ followed by 26 canonical Base32 symbols",
        ),
        (
            R1_HEADER.replace(PAIR_A, "TP_ABCDEFGHIJKLMNOPQRSTUVWXY1"),
            "FASTQ final pair ID is not TP_ followed by 26 canonical Base32 symbols",
        ),
        // Un símbolo de más: sin la comprobación de longitud llegaría a
        // copiar 27 símbolos en un hueco de 26.
        (
            R1_HEADER.replace(PAIR_A, "TP_ABCDEFGHIJKLMNOPQRSTUVWXYYA"),
            "FASTQ final pair ID is not TP_ followed by 26 canonical Base32 symbols",
        ),
        // Veintiséis símbolos Base32, pero no canónico: el último lleva 01 en
        // los dos bits que sobran de los 128 (DTC-087).
        (
            R1_HEADER.replace(PAIR_A, "TP_ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
            "FASTQ final pair ID is not TP_ followed by 26 canonical Base32 symbols",
        ),
        (
            R1_HEADER.replace(LIBRARY, "LIB_ABCDEFGHIJKLMNOPQRSTUVWXYZA"),
            "FASTQ final library ID is not LIB_ followed by 26 Base32 symbols",
        ),
        (
            R1_HEADER.replace(PAIR_A, "TP_abcdefghijklmnopqrstuvwxyz"),
            "FASTQ final pair ID is not TP_ followed by 26 canonical Base32 symbols",
        ),
        (
            R1_HEADER.replace(LIBRARY, "XXX_ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
            "FASTQ final library ID is not LIB_ followed by 26 Base32 symbols",
        ),
        (
            R1_HEADER.replace(LIBRARY, "LIB_ABCDEFGHIJKLMNOPQRSTUVWX"),
            "FASTQ final library ID is not LIB_ followed by 26 Base32 symbols",
        ),
    ] {
        assert_eq!(
            rejection(&record(&header, b"ACGT", b"+", b"IIII"), TechnicalMate::R1),
            detail,
            "{header}"
        );
    }
}

/// Hueco 4: la longitud de cada lectura es la del perfil. Una más, una menos y
/// longitudes mezcladas dentro del fichero.
#[test]
fn read_length_must_equal_the_profile_length() {
    const DETAIL: &str = "FASTQ final read length differs from the profile read length";
    assert_eq!(
        rejection_expecting(&record(R1_HEADER, b"ACGTA", b"+", b"IIIII"), &identity()),
        DETAIL
    );
    assert_eq!(
        rejection_expecting(&record(R1_HEADER, b"ACG", b"+", b"III"), &identity()),
        DETAIL
    );
    let mut mixed = valid_record();
    mixed.extend_from_slice(&record(
        &R1_HEADER.replace(PAIR_A, PAIR_B),
        b"ACGTAC",
        b"+",
        b"IIIIII",
    ));
    assert_eq!(rejection_expecting(&mixed, &identity()), DETAIL);
}

/// La longitud esperada misma debe estar en 1..=150 (DTPD-106).
#[test]
fn an_expected_read_length_outside_one_to_150_is_rejected() {
    for read_length in [0_u16, 151] {
        let expected = FinalFastqIdentity {
            read_length,
            ..identity()
        };
        assert_eq!(
            rejection_expecting(&valid_record(), &expected),
            "FASTQ final expected read length is outside 1..=150"
        );
    }
}

/// Hueco 5 (REQ-DAT-72): una `N` con cualquier calidad que no sea Phred 0.
#[test]
fn an_n_base_must_carry_phred_zero() {
    for (sequence, qualities) in [
        (&b"ACNT"[..], &b"IIII"[..]),
        (b"NCGT", b"\"!!!"),
        (b"ACGN", b"!!!\""),
    ] {
        assert_eq!(
            rejection_expecting(&record(R1_HEADER, sequence, b"+", qualities), &identity()),
            "FASTQ final N base does not carry Phred 0",
        );
    }
}

/// Hueco 6: R1 y R2 de una pareja de longitud distinta. Lo rechaza el juez del
/// mate que se sale de la longitud común, por el motivo de la longitud: con
/// una longitud esperada igual para los dos mates, la comparación directa de
/// R1 con R2 no puede dispararse sola.
#[test]
fn mates_of_different_lengths_are_rejected() {
    let pair = Pair::new(
        &record(R1_HEADER, b"ACGT", b"+", b"IIII"),
        &record(R2_HEADER, b"ACGTA", b"+", b"IIIII"),
    );
    match pair.judge() {
        Err(FastqValidationError::Structural { detail, .. }) => assert_eq!(
            detail,
            "FASTQ final read length differs from the profile read length"
        ),
        other => panic!("expected the read-length rejection, got {other:?}"),
    }
}

/// Hueco 7: los recuentos del resumen frente al contenido, uno por uno. En
/// cada caso los otros dos son correctos y solo el alterado puede rechazar.
#[test]
fn summary_counts_that_differ_from_the_content_are_rejected() {
    let pair = canonical_pair();
    let right = pair.counts;
    assert_eq!(
        (right.pair_count, right.base_count, right.n_base_count),
        (1, 8, 1),
        "the fixture must contain a pair, eight bases and one N"
    );
    for (counts, detail) in [
        (
            FinalFastqCounts {
                pair_count: 2,
                ..right
            },
            "FASTQ final pair count differs from the published summary",
        ),
        (
            FinalFastqCounts {
                base_count: 9,
                ..right
            },
            "FASTQ final base count differs from the published summary",
        ),
        (
            FinalFastqCounts {
                n_base_count: 0,
                ..right
            },
            "FASTQ final N count differs from the published summary",
        ),
    ] {
        match pair.judge_expecting(&identity(), &counts) {
            Err(FastqValidationError::Structural { detail: found, .. }) => {
                assert_eq!(found, detail);
            }
            other => panic!("expected {detail}, got {other:?}"),
        }
    }
    pair.judge_expecting(&identity(), &right)
        .expect("the true counts are accepted");
}

/// Un mate con más registros que parejas publica el resumen, y el otro con las
/// justas. Se rechaza por el recuento de parejas y no por el de huérfanos, que
/// vendría después: sin esta prueba, comprobar solo uno de los dos mates no
/// la haría fallar, porque la comparación de cabeceras rechazaría igual aunque
/// por otro motivo.
#[test]
fn a_mate_with_more_records_than_published_pairs_is_rejected_by_count() {
    let one = |header: &str| record(header, b"ACNT", b"+", b"II!~");
    let two = |header: &str| {
        let mut bytes = one(header);
        bytes.extend_from_slice(&one(&header.replace(PAIR_A, PAIR_B)));
        bytes
    };
    for (r1, r2) in [
        (two(R1_HEADER), one(R2_HEADER)),
        (one(R1_HEADER), two(R2_HEADER)),
    ] {
        let pair = Pair::new(&r1, &r2);
        let counts = FinalFastqCounts {
            pair_count: 1,
            base_count: 8,
            n_base_count: 1,
        };
        match pair.judge_expecting(&identity(), &counts) {
            Err(FastqValidationError::Structural { detail, .. }) => assert_eq!(
                detail,
                "FASTQ final pair count differs from the published summary"
            ),
            other => panic!("expected the pair-count rejection, got {other:?}"),
        }
    }
}
