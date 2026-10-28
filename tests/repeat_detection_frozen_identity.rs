//! **Huella de cambios involuntarios en los escenarios sin microsatélites.**
//!
//! # Qué protege, dicho sin rodeos
//!
//! El trabajo de la detección de microsatélites (plan F6.4.a) toca la
//! construcción de la configuración efectiva, la preparación de la referencia,
//! el generador de la matriz de SNV y la publicación. Ningún escenario **sin
//! microsatélites** debería notarlo. Este fichero congela, como literales, el
//! SHA-256 de **todos los ficheros publicados** por cuatro escenarios pequeños
//! que no declaran un solo microsatélite, para que un cambio **involuntario**
//! de bytes mate la prueba en vez de pasar inadvertido.
//!
//! # Qué NO protege
//!
//! **No es una obligación de compatibilidad.** La etiqueta `v0.2.0-test` es una
//! versión alfa y no una publicación (decisión del responsable, 28/10/2026,
//! nota fechada de DTC-036): los esquemas y formatos se retocan en su versión
//! actual. Si un cambio **decidido** modifica estos bytes —por ejemplo, un campo
//! nuevo en la configuración efectiva que el plan apruebe—, se actualizan los
//! literales **tras decidirlo por escrito**, no para que la prueba calle. Y los
//! escenarios que **declaran** microsatélites quedan fuera a propósito: con el
//! modo único de detección cambian a propósito y no se congelan aquí.
//!
//! # De dónde salen los literales, y por qué no es circular
//!
//! Se obtuvieron el 28 de octubre de 2026 ejecutando **este mismo fichero**
//! contra un árbol aparte extraído con `git archive d101f8f` (el último commit
//! anterior al trabajo de detección), con su propio directorio de compilación y
//! su propio binario. No se calcularon con el árbol de trabajo. Los cuatro
//! escenarios ejecutan el binario `run`, no la biblioteca.
//!
//! # Qué tiene cada escenario, y qué ejercita
//!
//! 1. **Referencia simulada con matriz de SNV generada**: dos secuencias,
//!    individuos diploides de dos poblaciones, sin evento técnico. Cubre la
//!    construcción de la referencia al azar, el generador de la matriz, la
//!    verdad y el VCF.
//! 2. **Referencia suministrada, matriz suministrada y cambios directos**: un
//!    cambio de base, una inserción y una deleción escritos a mano.
//! 3. **Evento técnico sin PCR**: FASTQ con el perfil sin amplificación.
//! 4. **Evento técnico con PCR**: lo mismo con preparación amplificada.
//!
//! # Qué discrimina, comprobado el 28/10/2026
//!
//! Sobre una copia del árbol, fuera de git y con copia de seguridad aparte:
//!
//! - Consumir **una palabra de azar más** al elegir los loci de la matriz
//!   generada (`src/population_snvs.rs`, SHA-256 del fichero
//!   `ca91e37e…19b52` antes, `65a21c3a…f4335` con la mutación) mata **solo** el
//!   escenario 1, que es el único con matriz generada; los otros tres
//!   siguen en verde, como debe ser.
//! - Añadir **un campo** a la configuración efectiva (`src/results.rs`,
//!   `96df83de…2eb0` antes, `fe705d6a…f6a` con la mutación) mata **los cuatro**.
//!
//! Restaurado el original, los SHA vuelven a ser los de antes.
//!
//! **Lo que no hace:** no ejercita la detección, que no existe todavía; solo
//! vigila que el trabajo que la prepara no mueva lo que no debe mover.
//!
//! # Qué hacer cuando esta prueba falle
//!
//! El mensaje imprime **la tabla observada completa** de cada escenario. No se
//! copia sin más: una huella que cambia es o un defecto o una decisión que hay
//! que escribir. Actualizar el literal para que calle es lo que este fichero
//! existe para impedir.
//!
//! # Honestidad sobre el compresor
//!
//! Los `.gz` los comprime `noodles-bgzf`, que puede elegir rutina según el
//! procesador. Por eso cada fichero comprimido lleva **dos** huellas: la del
//! fichero y la de su contenido descomprimido (entrada `#content`). Si solo
//! fallara la primera, el culpable sería el compresor y no el productor.

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use noodles_bgzf as bgzf;
use serde_json::Value;
use sha2::{Digest, Sha256};

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

// ---------------------------------------------------------------------------
// Los literales congelados, obtenidos con el commit `d101f8f` compilado aparte.
// Cada fila es `(ruta publicada relativa, SHA-256)`.
// ---------------------------------------------------------------------------

const SIMULATED_GENERATED_MATRIX: &[(&str, &str)] = &[
    (
        "effective-config.json",
        "7603d6c6c86cb2d71b826d804f4a66c89139ff8fd2c5f2696839527662c33ba0",
    ),
    (
        "inputs/config/general.toml",
        "0cbcbd2628e83ed522d91e2a17a14a05d51f9804bbf8a6967acc86e83f67b37d",
    ),
    (
        "manifest.json",
        "0361eca167d548c4dbb8734a56bd496324b9cb6490e3e1f2db0b6eb0f245d0b5",
    ),
    (
        "outputs/truth.jsonl",
        "703fbbd00fc52b391853e54432dffd53e355ff0bc04bba14f8352bf88c222047",
    ),
    (
        "population_snvs.metadata.json",
        "b3431595e143a606ee0cc70c7364687600f2221ff1ebcc5e94983bbdf77f10b7",
    ),
    (
        "population_snvs.tsv",
        "db4d45285bf44d02aa70012f31bd4163ede5b611ca4b252d9ab1a735e843f611",
    ),
    (
        "resources/references/ref-main.fa",
        "44f3e02f87d360898065eb50863291d2ce131f142c558e756e25ce4152817e15",
    ),
];
const SUPPLIED_DIRECT_CHANGES: &[(&str, &str)] = &[
    (
        "effective-config.json",
        "5b3b119e9b8e6f60f5e7f2ba8b432f3749e5c1f46bc405be6f45446989404b76",
    ),
    (
        "inputs/config/general.toml",
        "94dd0daf1372dd86d8341975b0e4bc1bf70eefa526430612f492ff6f4a1d74ac",
    ),
    (
        "manifest.json",
        "5f1e0e56b6be224eb41a053db55f0ee209e5e091d630b5a99da53732b5ccffa9",
    ),
    (
        "outputs/truth.jsonl",
        "77e9d9971a6bf5494c1853f383ad3c9eca7eaabc259f9bcf924e2e0c3a3e06c1",
    ),
    (
        "population_snvs.metadata.json",
        "ae32e2a8822a764c13a3c612155640de4b22e92dd57768ad0a92fd94434e52f1",
    ),
    (
        "population_snvs.tsv",
        "8e67d323893dfa26c3e557cac1c1f9c6d8861c18ea06f7ab5e94618e45c547ab",
    ),
    (
        "resources/references/ref-main.fa",
        "1dc7783f375c11b5422f1e06513b7e793b358342c19ec879606b325462e665d5",
    ),
];
const TECHNICAL_WITHOUT_PCR: &[(&str, &str)] = &[
    (
        "effective-config.json",
        "07195cda2fce11a371b2c5722edd04f0d361260ec6457e15ba0123299e03dae1",
    ),
    (
        "events/reads-main/manifest.json",
        "a18d486a7c5e510727db105d78730cabd84ec547353bd1b66c6a2002032cd48a",
    ),
    (
        "inputs/config/general.toml",
        "45afa81f59557d56b01761783a55c95431bfa1eab2a851f41cefa8652335f261",
    ),
    (
        "manifest.json",
        "3e37e690a1485dba743ebd1ea2fdcb839e0fdf2851c1b1b31c63c5ae0550bfbc",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz",
        "044c22fe0ee24df20b2af833742f766cdc1a3be581d5c11ec126b81cf3a41d5d",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz#content",
        "5847b905e9cb8baba47dab3246fc7049002b8827d9b81be6e499a724e1744a62",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz",
        "43a605a7917d7941aaaf7c5cc6a163e76573b30c2fa5a348b976e85bb9308ca8",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz#content",
        "cb13209ca2591026d4052cdc16f5540ef720874488b62e3771a478fbf77fd23b",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R1.fastq.gz",
        "b7d2f688a545f112b07ea9cefccd8db77b91c9aaa85abf2dfac8b13a6762d609",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R1.fastq.gz#content",
        "33b811429e6e48166c3f282e09b66f3380cbe46c43b6525eb956a2bfa33017ff",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R2.fastq.gz",
        "3a32a86250cde242132b5da189505ddf0e069160e47207e18cba8f515db231b5",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R2.fastq.gz#content",
        "7045c794a4f10e1fbd1f89bab47bf567312d948e1216fbb5cdc9f8cbef79629c",
    ),
    (
        "outputs/truth.jsonl",
        "b3cf45c52833032ee387fde7243cf5c967e3c7aaa353b8b2d83ec490cd675777",
    ),
    (
        "population_snvs.metadata.json",
        "ae32e2a8822a764c13a3c612155640de4b22e92dd57768ad0a92fd94434e52f1",
    ),
    (
        "population_snvs.tsv",
        "8e67d323893dfa26c3e557cac1c1f9c6d8861c18ea06f7ab5e94618e45c547ab",
    ),
    (
        "resources/references/ref-main.fa",
        "1dc7783f375c11b5422f1e06513b7e793b358342c19ec879606b325462e665d5",
    ),
];
const TECHNICAL_WITH_PCR: &[(&str, &str)] = &[
    (
        "effective-config.json",
        "8eb46c3eb30ca7690d6eb2068d619cf3c2a99faa919d3f01a39d9d61c28e07dc",
    ),
    (
        "events/reads-main/manifest.json",
        "7025cfd5fe7011bcbeb2d88f2d580d980f49adc8b6fa8d409cc23055a94f731d",
    ),
    (
        "inputs/config/general.toml",
        "f330421687f142f33586318e97af36836f5b35b04e7e72c5974d7b31e86300a1",
    ),
    (
        "manifest.json",
        "831ddb439c512469b9aee0e3b8bc6857cd30fdbcbf8c84e96f9c796677703b66",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz",
        "07693f8de5e51d3e29a2e5e213c1617f0d45bf83614ec2ad4bb1e7fa6c3d7f2c",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz#content",
        "3d29046698c116be38feb7ebc772e8e66b918d5be3952b0858cf18d5584c7f35",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz",
        "6c7a20217b003a05a36297a44df1347934d638b5bd9d59454623ba301d749677",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz#content",
        "132378adebb31e3e401926dee4f9c7b17ac731ae14f3ae4c3b2c84eb76459447",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R1.fastq.gz",
        "0d61c8825b7d26e52ef307547286788eaa9b78ec1b7a70273f95516b01d458ef",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R1.fastq.gz#content",
        "e8c993796169134acb85f5b4bd4862ac4a10d7d76c1f5d2d5d5c1aadbc068350",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R2.fastq.gz",
        "2497f5025df98eb2a49912769965a8da5174593ebf4c8596ce17909949988050",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R2.fastq.gz#content",
        "285e1b16ed714399067b01d4d3c72946bb70740ddfce5328ef779748a14cc77a",
    ),
    (
        "outputs/truth.jsonl",
        "1f168f918bf9eff7f0750a8e7ed5158bf5d402cfb4cc08b48eb79620f5516776",
    ),
    (
        "population_snvs.metadata.json",
        "ae32e2a8822a764c13a3c612155640de4b22e92dd57768ad0a92fd94434e52f1",
    ),
    (
        "population_snvs.tsv",
        "8e67d323893dfa26c3e557cac1c1f9c6d8861c18ea06f7ab5e94618e45c547ab",
    ),
    (
        "resources/references/ref-main.fa",
        "1dc7783f375c11b5422f1e06513b7e793b358342c19ec879606b325462e665d5",
    ),
];

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-huella-deteccion-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("reloj")
            .as_nanos()
    ))
}

/// Directorio temporal de un escenario: se borra al soltarse.
struct TemporaryRoot(PathBuf);

impl Drop for TemporaryRoot {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Bases pseudoaleatorias reproducibles (xorshift64), escritas aquí y sin
/// depender del programa.
fn bases(seed: u64, length: usize) -> Vec<u8> {
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

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn fasta_of(sequences: &[(&str, &[u8])]) -> String {
    let mut text = String::new();
    for (id, bases) in sequences {
        text.push('>');
        text.push_str(id);
        text.push('\n');
        for chunk in bases.chunks(70) {
            text.push_str(std::str::from_utf8(chunk).expect("ASCII"));
            text.push('\n');
        }
    }
    text
}

fn declared_sequences(sequences: &[(&str, &[u8])]) -> String {
    sequences
        .iter()
        .map(|(id, bases)| {
            format!(
                "[[reference.sequences]]\nid = \"{id}\"\nclass = \"chromosome\"\nlength = {}\n\n",
                bases.len()
            )
        })
        .collect()
}

/// Una fila de la matriz suministrada. La base de referencia se **lee de la
/// secuencia**, no se escribe a mano.
fn matrix_row(
    sequences: &[(&str, &[u8])],
    contig: &str,
    position: usize,
    population: &str,
) -> String {
    let bases = sequences
        .iter()
        .find(|(id, _)| *id == contig)
        .expect("secuencia del SNV")
        .1;
    let reference = bases[position - 1];
    let alternate = *b"ACGT"
        .iter()
        .find(|base| **base != reference)
        .expect("alternativa");
    format!(
        "{contig}\t{position}\t{}\t{}\t{population}\t0.5\n",
        char::from(reference),
        char::from(alternate)
    )
}

const TECHNICAL_LIBRARIES: &str = "\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\n\
     i7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n\n\
     [[technical_event.libraries]]\nindividual_id = \"sample-B\"\n\
     i7 = \"TTTTAAAA\"\ni5 = \"CCCCGGGG\"\n";

/// Ejecuta `run` con el binario del proyecto y devuelve el directorio
/// publicado.
///
/// **La configuración se pasa por una ruta relativa, desde el propio
/// directorio del escenario.** La ruta tal como se escribe en la línea de
/// órdenes se publica en la configuración efectiva y en el manifiesto
/// (`source_path`, `original_path`): con una ruta absoluta, que cambia en cada
/// ejecución porque el directorio temporal es único, esos dos ficheros no
/// tendrían una huella estable.
fn run_binary(root: &Path, config: &Path) -> PathBuf {
    let destination = root.join("published");
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .current_dir(root)
        .arg("run")
        .arg(config.file_name().expect("nombre de la configuración"))
        .arg("--output")
        .arg("published")
        .output()
        .expect("binario del simulador");
    assert!(
        output.status.success(),
        "la ejecución falló\n{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    destination
}

/// La referencia suministrada de los escenarios 2 a 4, con una matriz
/// suministrada de dos poblaciones.
fn supplied_scenario(root: &Path, body: &str, populations: &str) -> PathBuf {
    let chr1 = bases(0x9E37_79B9_7F4A_7C15, 1_800);
    let chr2 = bases(0xD1B5_4A32_D192_ED03, 900);
    let sequences: [(&str, &[u8]); 2] = [("chr1", &chr1), ("chr2", &chr2)];
    fs::write(root.join("reference.fa"), fasta_of(&sequences)).expect("FASTA");
    let mut matrix = String::from("contig\tposition\tref\talt\tpopulation_id\talt_af\n");
    for (contig, position, population) in [
        ("chr1", 120, "P1"),
        ("chr1", 120, "P2"),
        ("chr1", 640, "P1"),
        ("chr1", 640, "P2"),
        ("chr1", 1_500, "P1"),
        ("chr1", 1_500, "P2"),
        ("chr2", 75, "P1"),
        ("chr2", 75, "P2"),
        ("chr2", 410, "P1"),
        ("chr2", 410, "P2"),
    ] {
        matrix.push_str(&matrix_row(&sequences, contig, position, population));
    }
    fs::write(root.join("population-snvs.tsv"), matrix).expect("matriz");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"population-snvs.tsv\"\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             {}{populations}{body}",
            declared_sequences(&sequences)
        ),
    )
    .expect("configuración");
    config
}

const TWO_DIPLOID_POPULATIONS: &str = "[[initial_populations]]\npopulation_id = \"P1\"\n\
     individual_count = 1\nploidy = 2\ngenotype_model = \"hwe\"\n\
     individual_ids = [\"sample-A\"]\n\n\
     [[initial_populations]]\npopulation_id = \"P2\"\nindividual_count = 1\nploidy = 2\n\
     genotype_model = \"hwe\"\nindividual_ids = [\"sample-B\"]\n\n";

/// Escenario 1.
fn simulated_generated_matrix(root: &Path) -> PathBuf {
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
             [population_snvs]\nmode = \"generate\"\n\n[population_snvs.generate]\n\
             count_mode = \"exact\"\nlocus_count = 12\nspatial_mode = \"uniform\"\n\
             substitution_profile = \"uniform_snv12_v1\"\n\n\
             [[population_snvs.generate.populations]]\npopulation_id = \"P1\"\n\
             frequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n\n\
             [[population_snvs.generate.populations]]\npopulation_id = \"P2\"\n\
             frequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n\n\
             [reference]\nkind = \"simulated\"\nid = \"ref-main\"\n\n\
             [reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1500\n\n\
             [[reference.sequences]]\nid = \"chr2\"\nclass = \"chromosome\"\nlength = 900\n\n\
             {TWO_DIPLOID_POPULATIONS}"
        ),
    )
    .expect("configuración");
    config
}

/// Escenario 2: tres cambios directos escritos a mano.
fn supplied_direct_changes(root: &Path) -> PathBuf {
    let chr1 = bases(0x9E37_79B9_7F4A_7C15, 1_800);
    let reference = chr1[299];
    let alternate = *b"ACGT"
        .iter()
        .find(|base| **base != reference)
        .expect("alternativa");
    let events = format!(
        "[[direct_events]]\nevent_id = \"by-hand-snp\"\nindividual_id = \"sample-A\"\nphase = 0\n\
         kind = \"snp\"\ncontig = \"chr1\"\nposition_1based = 300\nref = \"{}\"\nalt = \"{}\"\n\n\
         [[direct_events]]\nevent_id = \"by-hand-insertion\"\nindividual_id = \"sample-A\"\n\
         phase = 1\nkind = \"insertion\"\ncontig = \"chr1\"\npoint = 1000\nbases = \"ACGTA\"\n\n\
         [[direct_events]]\nevent_id = \"by-hand-deletion\"\nindividual_id = \"sample-B\"\n\
         phase = 0\nkind = \"deletion\"\ncontig = \"chr2\"\nstart = 600\nend = 605\n",
        char::from(reference),
        char::from(alternate)
    );
    supplied_scenario(root, &events, TWO_DIPLOID_POPULATIONS)
}

fn technical_event(preparation: &str) -> String {
    format!(
        "[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [\"sample-A\", \"sample-B\"]\n\
         target_coverage_wgs = 2\n{preparation}{TECHNICAL_LIBRARIES}"
    )
}

/// Escenario 3.
fn technical_without_pcr(root: &Path) -> PathBuf {
    supplied_scenario(root, &technical_event(""), TWO_DIPLOID_POPULATIONS)
}

/// Escenario 4.
fn technical_with_pcr(root: &Path) -> PathBuf {
    supplied_scenario(
        root,
        &technical_event(
            "library_preparation = \"pcr\"\n\n[technical_event.pcr]\n\
             duplicate_fraction = 0.2\ncycles = 10\nefficiency = 0.9\n",
        ),
        TWO_DIPLOID_POPULATIONS,
    )
}

fn walk(root: &Path, directory: &Path, found: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(directory).expect("directorio publicado") {
        let path = entry.expect("entrada").path();
        if path.is_dir() {
            walk(root, &path, found);
        } else {
            found.push(path.strip_prefix(root).expect("relativa").to_path_buf());
        }
    }
}

/// El contenido descomprimido de un fichero BGZF.
fn decompressed(path: &Path) -> Option<Vec<u8>> {
    let mut text = Vec::new();
    bgzf::io::Reader::new(fs::File::open(path).ok()?)
        .read_to_end(&mut text)
        .ok()?;
    Some(text)
}

/// Huella de **todos** los ficheros publicados: la del fichero y, para los
/// comprimidos, la de su contenido.
fn fingerprint(destination: &Path) -> BTreeMap<String, String> {
    let mut paths = Vec::new();
    walk(destination, destination, &mut paths);
    let mut table = BTreeMap::new();
    for relative in paths {
        let full = destination.join(&relative);
        let name = relative.to_string_lossy().into_owned();
        table.insert(name.clone(), digest(&fs::read(&full).expect("fichero")));
        if name.ends_with(".gz") {
            let content =
                decompressed(&full).unwrap_or_else(|| panic!("{name}: no es un BGZF legible"));
            table.insert(format!("{name}#content"), digest(&content));
        }
    }
    table
}

fn render(table: &BTreeMap<String, String>) -> String {
    table
        .iter()
        .map(|(path, sha)| format!("    (\"{path}\", \"{sha}\"),\n"))
        .collect()
}

/// Control positivo: el escenario es de verdad un escenario **sin
/// microsatélites**, de modo que la huella no es la de un caso que la
/// detección debería cambiar.
fn assert_without_microsatellites(destination: &Path) {
    let effective: Value = serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("configuración efectiva"),
    )
    .expect("JSON");
    let loci = effective.get("repeat_loci");
    assert!(
        loci.is_none_or(|value| value.is_null() || value.as_array().is_some_and(Vec::is_empty)),
        "el escenario congelado no puede declarar microsatélites: {loci:?}"
    );
}

fn check(name: &str, frozen: &[(&str, &str)], build: fn(&Path) -> PathBuf, required: &[&str]) {
    let root = temporary();
    fs::create_dir_all(&root).expect("raíz del escenario");
    let _cleanup = TemporaryRoot(root.clone());
    let config = build(&root);
    let destination = run_binary(&root, &config);
    assert_without_microsatellites(&destination);
    let observed = fingerprint(&destination);
    // Control positivo: el escenario publica de verdad lo que dice ejercitar,
    // de modo que la huella no es la de un paquete vacío o truncado.
    for needle in required {
        assert!(
            observed.keys().any(|path| path.contains(needle)),
            "{name}: el escenario tiene que publicar algo que contenga «{needle}»; publicó\n{}",
            render(&observed)
        );
    }
    let expected: BTreeMap<String, String> = frozen
        .iter()
        .map(|(path, sha)| ((*path).to_owned(), (*sha).to_owned()))
        .collect();
    assert_eq!(
        observed,
        expected,
        "{name}: los bytes publicados han cambiado. NO se actualiza el literal sin decidir por \
         escrito que el cambio es intencionado. Tabla observada:\n{}",
        render(&observed)
    );
}

/// **Referencia simulada con matriz generada: los bytes son los congelados.**
///
/// Muere si cambia la construcción de la referencia al azar, el consumo de
/// azar del generador de la matriz, la verdad, el VCF o cualquier fichero del
/// paquete, y también si la configuración efectiva gana un campo.
#[test]
fn the_simulated_reference_with_a_generated_matrix_keeps_its_frozen_bytes() {
    check(
        "simulada con matriz generada",
        SIMULATED_GENERATED_MATRIX,
        simulated_generated_matrix,
        &["truth", "ref-main.fa", "population_snvs"],
    );
}

/// **Referencia suministrada con matriz suministrada y cambios directos: los
/// bytes son los congelados.**
#[test]
fn the_supplied_reference_with_direct_changes_keeps_its_frozen_bytes() {
    check(
        "suministrada con cambios directos",
        SUPPLIED_DIRECT_CHANGES,
        supplied_direct_changes,
        &["truth", "ref-main.fa"],
    );
}

/// **Evento técnico sin PCR: los bytes, FASTQ incluidos, son los congelados.**
#[test]
fn the_technical_event_without_pcr_keeps_its_frozen_bytes() {
    check(
        "evento técnico sin PCR",
        TECHNICAL_WITHOUT_PCR,
        technical_without_pcr,
        &["sample-A_R1.fastq.gz", "sample-B_R2.fastq.gz"],
    );
}

/// **Evento técnico con PCR: los bytes, FASTQ incluidos, son los congelados.**
#[test]
fn the_technical_event_with_pcr_keeps_its_frozen_bytes() {
    check(
        "evento técnico con PCR",
        TECHNICAL_WITH_PCR,
        technical_with_pcr,
        &["sample-A_R1.fastq.gz", "sample-B_R2.fastq.gz"],
    );
}
