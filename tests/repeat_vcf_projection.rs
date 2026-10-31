//! F6.5.1.1: la proyección al VCF de los alelos de repetición en tándem.
//!
//! # Qué afirman estas pruebas y cómo están elegidos sus fixtures
//!
//! Lo que hay que demostrar es que el VCF publica los campos de repetición
//! **del propio estándar `VCFv4.5`**, que los alelos van siempre explícitos
//! —sin ningún umbral de tamaño— y que cada valor acompaña al alelo que le
//! corresponde. Para que las pruebas puedan fallar:
//!
//! - el fixture **compuesto** declara un locus de dos motivos y dos alelos que
//!   difieren **en los dos bloques**, con bases y copias todas distintas. Con
//!   un locus de un solo motivo las listas de motivo y de longitud de motivo
//!   serían iguales para todos los alelos y una lista desplazada sobreviviría;
//! - la correspondencia entre los valores y la longitud del alelo se calcula
//!   **desde el propio VCF** —`REF`, `ALT` y el intervalo declarado del
//!   locus—, no desde las listas que se están juzgando;
//! - la prueba del tracto largo usa **diez mil bases**, que es la cifra con la
//!   que el estándar ilustra la forma simbólica; y la búsqueda negativa lleva
//!   su control positivo, para que un patrón que no encontrara nunca nada no
//!   pueda pasar por una ausencia.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_fasta::{SequenceExpectation, validate_supplied_fasta};
use genome_population_simulator::truth_validation::read_truth_jsonl;
use genome_population_simulator::vcf_semantic_validation::{
    validate_published_vcf_semantics, validate_vcf_against_truth,
};

#[path = "support/detection_section.rs"]
mod detection_section;

static NONCE: AtomicU64 = AtomicU64::new(0);

const COMPOSITE_TRACT_BASES: i64 = 67;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

/// Las cinco declaraciones canónicas, **transcritas del estándar** y no
/// importadas del código que las escribe.
const CANONICAL_DECLARATIONS: [&str; 5] = [
    "##INFO=<ID=RN,Number=A,Type=Integer,Description=\"Total number of repeat sequences in this allele\">",
    "##INFO=<ID=RUS,Number=.,Type=String,Description=\"Repeat unit sequence of the corresponding repeat sequence\">",
    "##INFO=<ID=RUL,Number=.,Type=Integer,Description=\"Repeat unit length of the corresponding repeat sequence\">",
    "##INFO=<ID=RUC,Number=.,Type=Float,Description=\"Repeat unit count of corresponding repeat sequence\">",
    "##INFO=<ID=RB,Number=.,Type=Integer,Description=\"Total number of bases in the corresponding repeat sequence\">",
];

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f6-5-1-1-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Escenario con referencia suministrada, matriz mínima, dos individuos
/// diploides y VCF de texto. `bases` y `loci` se pegan tal cual.
fn scenario(bases: &str, loci: &str) -> (PathBuf, PathBuf) {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let length = bases.len();
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {length}
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"

[vcf_output]
mode = "text"
"#
        )),
    )
    .expect("configuración");
    (root, config)
}

fn publish(bases: &str, loci: &str) -> (PathBuf, PathBuf) {
    let (root, config) = scenario(bases, loci);
    let effective = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective, &destination).expect("publicación");
    (root, destination)
}

fn vcf_text(destination: &Path) -> String {
    fs::read_to_string(destination.join("outputs/truth.vcf")).expect("VCF publicado")
}

/// Un registro del VCF descompuesto en sus columnas.
struct Record {
    position: u64,
    reference: String,
    alternates: Vec<String>,
    info: BTreeMap<String, String>,
    genotypes: Vec<String>,
}

fn records(text: &str) -> Vec<Record> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            Record {
                position: fields[1].parse().expect("POS"),
                reference: fields[3].to_owned(),
                alternates: fields[4].split(',').map(str::to_owned).collect(),
                info: fields[7]
                    .split(';')
                    .filter_map(|field| field.split_once('='))
                    .map(|(key, value)| (key.to_owned(), value.to_owned()))
                    .collect(),
                genotypes: fields[9..]
                    .iter()
                    .map(|value| (*value).to_owned())
                    .collect(),
            }
        })
        .collect()
}

fn repeat_records(text: &str) -> Vec<Record> {
    let found = records(text)
        .into_iter()
        .filter(|record| record.info.contains_key("RN"))
        .collect::<Vec<_>>();
    assert!(
        !found.is_empty(),
        "el escenario tiene que producir algún registro de repetición"
    );
    found
}

fn repeat_record(text: &str) -> Record {
    let mut found = repeat_records(text);
    assert_eq!(
        found.len(),
        1,
        "el escenario tiene que producir exactamente un registro de repetición"
    );
    found.pop().expect("registro")
}

/// Comprueba que los valores de un registro se corresponden con los alelos que
/// acompañan, **en su mismo orden**.
///
/// La correspondencia se calcula desde `REF`, `ALT` y las bases que el tracto
/// ocupa en la referencia, que son datos ajenos a las listas juzgadas: si una
/// de ellas se desplaza respecto de la de alelos, las bases publicadas dejan de
/// ser las de la secuencia del alelo. Devuelve las cifras de bases vistas, para
/// que quien llama pueda exigir que sean distintas entre sí y la comparación no
/// sea vacua.
///
/// **El reparto de la lista entre los alelos solo se alimenta aquí con
/// recuentos iguales**, porque los alelos de un locus comparten su lista de
/// motivos y por tanto traen todos el mismo número de bloques. El caso de
/// recuentos distintos —dos bloques, ninguno y uno— lo cubre una prueba
/// unitaria del oráculo, que lo construye a mano.
fn check_allele_correspondence(record: &Record, reference_tract_bases: i64) -> Vec<u64> {
    let counts = integers(record, "RN");
    assert_eq!(
        counts.len(),
        record.alternates.len(),
        "el recuento de secuencias repetidas es uno por alelo alternativo, y el de referencia no lleva"
    );
    let total = counts.iter().sum::<u64>() as usize;
    let motifs = record.info["RUS"].split(',').collect::<Vec<_>>();
    let motif_lengths = integers(record, "RUL");
    let unit_counts = decimals(record, "RUC");
    let bases = integers(record, "RB");
    assert_eq!(motifs.len(), total, "RUS no suma los recuentos de RN");
    assert_eq!(
        motif_lengths.len(),
        total,
        "RUL no suma los recuentos de RN"
    );
    assert_eq!(unit_counts.len(), total, "RUC no suma los recuentos de RN");
    assert_eq!(bases.len(), total, "RB no suma los recuentos de RN");

    let mut cursor = 0_usize;
    for (index, count) in counts.iter().enumerate() {
        let count = *count as usize;
        let slice = cursor..cursor + count;
        cursor += count;
        let published = bases[slice.clone()].iter().sum::<u64>() as i64;
        // **La comparación es directa desde que el registro no se recorta**:
        // `ALT` es la secuencia completa del alelo, de modo que las bases que
        // los campos declaran tienen que ser exactamente las que `ALT` mide.
        // Con la forma recortada esto había que rodearlo sumando y restando la
        // longitud del tracto de la referencia, y esa vuelta era precisamente
        // la señal de que la anotación describía algo que el `ALT` no contenía.
        assert_eq!(
            published,
            record.alternates[index].len() as i64,
            "las bases publicadas del alelo {index} no son las que su ALT mide"
        );
        assert_eq!(
            record.reference.len() as i64,
            reference_tract_bases,
            "REF tiene que llevar el tracto de la referencia entero"
        );
        for position in slice {
            assert_eq!(
                motif_lengths[position],
                motifs[position].len() as u64,
                "la longitud del motivo no es la de su secuencia"
            );
            assert!(
                (unit_counts[position] - bases[position] as f64 / motif_lengths[position] as f64)
                    .abs()
                    < 1e-9,
                "el número de copias no son las bases entre la longitud del motivo"
            );
        }
    }
    bases
}

fn integers(record: &Record, key: &str) -> Vec<u64> {
    record.info[key]
        .split(',')
        .map(|value| value.parse().expect("entero"))
        .collect()
}

fn decimals(record: &Record, key: &str) -> Vec<f64> {
    record.info[key]
        .split(',')
        .map(|value| value.parse().expect("decimal"))
        .collect()
}

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// Trescientas bases: cien de flanco, cien de tracto `(CA)50` y cien más.
fn simple_reference() -> String {
    format!("{}{}{}", "G".repeat(100), "CA".repeat(50), "G".repeat(100))
}

/// Un locus de un motivo con dos alelos de longitud distinta.
const SIMPLE_POLYMORPHIC: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 45
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 40
frequency = 0.5
"#;

/// Un locus cuyo único alelo es el de la referencia: no publica ninguna
/// variante y, por tanto, ningún campo de repetición.
const SIMPLE_REFERENCE_ONLY: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 50
frequency = 1.0
"#;

/// Referencia con un microsatélite **compuesto** `(AAAT)10(GGC)9`, que es la
/// forma de muchos marcadores forenses y la que el estándar describe con sus
/// listas paralelas.
///
/// Los dos motivos miden **distinto** —cuatro y tres— a propósito: con motivos
/// de igual longitud la lista de longitudes de motivo sale toda igual y
/// desplazarla es indistinguible.
fn composite_reference() -> String {
    format!(
        "{}{}{}{}",
        "C".repeat(100),
        "AAAT".repeat(10),
        "GGC".repeat(9),
        "T".repeat(100)
    )
}

/// Dos alelos que comparten su primer bloque y difieren en el segundo:
/// `(AAAT)10(GGC)6` y `(AAAT)10(GGC)3`.
///
/// Compartir el primer bloque es lo que hace que los dos se alineen en el
/// **mismo sitio** y salgan como dos `ALT` de un único registro. Es el único
/// fixture que combina las dos cosas a la vez —varios alelos en un registro y
/// varios bloques por alelo—, que es justo donde las listas aplanadas pueden
/// desalinearse.
const COMPOSITE_ALLELES: &str = r#"
[[reference.repeat_loci]]
id = "STR-COMP"
sequence = "chr1"
start = 101
end = 167
motifs = ["AAAT", "GGC"]

[[reference.repeat_loci.alleles]]
units = [10, 6]
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = [10, 3]
frequency = 0.5
"#;

/// Referencia con un tracto `(CAG)30` para declarar un alelo cuya longitud no
/// es múltiplo del motivo.
fn partial_unit_reference() -> String {
    format!("{}{}{}", "G".repeat(100), "CAG".repeat(30), "G".repeat(110))
}

/// Un alelo de sesenta y dos bases sobre un motivo de tres: veinte unidades y
/// **dos tercios**.
fn partial_unit_locus() -> String {
    let sequence = format!("{}CA", "CAG".repeat(20));
    assert_eq!(sequence.len(), 62);
    format!(
        r#"
[[reference.repeat_loci]]
id = "STR-PART"
sequence = "chr1"
start = 101
end = 190
motif = "CAG"

[[reference.repeat_loci.alleles]]
sequence = "{sequence}"
frequency = 1.0
"#
    )
}

/// Un alelo de **diez mil bases**, la cifra con la que el estándar ilustra la
/// forma simbólica, sobre un tracto de referencia de cien.
const TEN_THOUSAND_BASE_ALLELE: &str = r#"
[[reference.repeat_loci]]
id = "STR-LONG"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 6000

[[reference.repeat_loci.alleles]]
units = 5000
frequency = 1.0
"#;

// ---------------------------------------------------------------------------
// Pruebas
// ---------------------------------------------------------------------------

/// **Un diploide con dos alelos produce un solo registro multialélico, faseado,
/// y con el tracto entero en `REF` y `ALT`.**
///
/// *Qué hace que pueda fallar:* se cuentan los registros del locus, de modo
/// que dos registros para el mismo sitio —el error que DTC-064 prohíbe— harían
/// caer la cuenta. Y se comprueba que el registro **no está recortado**: `REF`
/// mide el tracto de la referencia entero y empieza donde el locus declara,
/// que es justo lo que la forma recortada no cumplía.
#[test]
fn a_diploid_with_two_repeat_alleles_yields_one_untrimmed_multiallelic_record() {
    let (root, destination) = publish(&simple_reference(), SIMPLE_POLYMORPHIC);
    let text = vcf_text(&destination);
    let record = repeat_record(&text);

    assert_eq!(
        record.position, 101,
        "el registro empieza donde empieza el tracto, sin ancla previa"
    );
    assert_eq!(
        record.reference.len(),
        100,
        "REF lleva el tracto de la referencia entero, no una base de ancla"
    );
    assert_eq!(
        record.reference,
        "CA".repeat(50),
        "y es literalmente el tracto que hay en la referencia"
    );
    for alternate in &record.alternates {
        assert!(
            alternate.starts_with("CACA") && alternate.len() >= 80,
            "cada ALT lleva la secuencia completa de su alelo: {} bases",
            alternate.len()
        );
    }
    assert_eq!(
        record.alternates.len(),
        2,
        "los dos alelos van como dos ALT de un único registro: {:?}",
        record.alternates
    );
    for genotype in &record.genotypes {
        assert!(
            genotype.contains('|'),
            "el genotipo tiene que ir faseado: {genotype}"
        );
    }
    let at_position = records(&text)
        .into_iter()
        .filter(|other| other.position == record.position)
        .count();
    assert_eq!(at_position, 1, "un solo registro por CHROM+POS");
    fs::remove_dir_all(&root).ok();
}

/// **La cabecera declara cada campo usado, con su número y su tipo, y ninguno
/// más.**
///
/// *Qué hace que pueda fallar:* las declaraciones se comparan con un literal
/// **transcrito del estándar** en este fichero, no con la constante del
/// productor, de modo que cambiar el número o el tipo del productor rompe la
/// comparación. Y el escenario cuyo único alelo es el de la referencia exige
/// que **no** aparezca ninguna: declararlas siempre pasaría la primera mitad y
/// moriría en la segunda.
#[test]
fn the_header_declares_every_tandem_repeat_field_used_and_no_other() {
    let (root, destination) = publish(&simple_reference(), SIMPLE_POLYMORPHIC);
    let text = vcf_text(&destination);
    for declaration in CANONICAL_DECLARATIONS {
        assert_eq!(
            text.lines().filter(|line| *line == declaration).count(),
            1,
            "falta o sobra la declaración {declaration}"
        );
    }
    let info_declarations = text
        .lines()
        .filter(|line| line.starts_with("##INFO=<"))
        .count();
    assert_eq!(
        info_declarations,
        CANONICAL_DECLARATIONS.len() + 1,
        "solo GPS_CASE y los cinco campos del estándar"
    );
    fs::remove_dir_all(&root).ok();

    let (bare_root, bare) = publish(&simple_reference(), SIMPLE_REFERENCE_ONLY);
    let bare_text = vcf_text(&bare);
    for declaration in CANONICAL_DECLARATIONS {
        assert!(
            !bare_text.contains(declaration),
            "sin ningún alelo publicado no se declara {declaration}"
        );
    }
    assert!(
        !bare_text.contains("RN="),
        "y ningún registro publica los campos"
    );
    fs::remove_dir_all(&bare_root).ok();
}

/// **El número de copias se declara y se emite como decimal, y un tracto con
/// unidad parcial produce un valor no entero.**
///
/// *Qué hace que pueda fallar:* el alelo mide sesenta y dos bases sobre un
/// motivo de tres. Con un tracto de longitud múltiplo del motivo el valor
/// saldría entero y una implementación que lo emitiera como entero pasaría.
#[test]
fn the_repeat_unit_count_is_a_decimal_and_a_partial_unit_is_not_an_integer() {
    let (root, destination) = publish(&partial_unit_reference(), &partial_unit_locus());
    let text = vcf_text(&destination);
    assert!(
        text.contains("##INFO=<ID=RUC,Number=.,Type=Float,"),
        "el número de copias se declara decimal"
    );
    let record = repeat_record(&text);
    let counts = decimals(&record, "RUC");
    assert_eq!(counts.len(), 1);
    assert!(
        (counts[0] - 62.0 / 3.0).abs() < 1e-9,
        "sesenta y dos bases entre un motivo de tres: {}",
        counts[0]
    );
    assert!(
        counts[0].fract() > 0.0,
        "y el valor publicado no es entero: {}",
        record.info["RUC"]
    );
    assert_eq!(integers(&record, "RB"), vec![62]);
    fs::remove_dir_all(&root).ok();
}

/// **Cada alelo alternativo lleva su número de copias explícito, el de
/// referencia no, y cada valor se corresponde con la longitud del alelo al que
/// acompaña, en su mismo orden.**
///
/// *Qué hace que pueda fallar:* se comprueban dos escenarios complementarios.
/// El de un motivo produce un registro **multialélico** cuyos dos alelos miden
/// ochenta y noventa bases, de modo que intercambiar sus valores rompe la
/// correspondencia. El **compuesto** da dos bloques por alelo con motivos,
/// copias y bases todos distintos, de modo que desplazar cualquiera de las
/// cuatro listas dentro del alelo también rompe algo. La correspondencia se
/// calcula desde `REF`, `ALT` y las bases del tracto en la referencia, que no
/// salen de las listas juzgadas.
#[test]
fn every_alternate_allele_carries_repeat_counts_matching_its_own_length() {
    let (root, destination) = publish(&simple_reference(), SIMPLE_POLYMORPHIC);
    let text = vcf_text(&destination);
    let record = repeat_record(&text);
    assert_eq!(
        record.alternates.len(),
        2,
        "el fixture necesita un registro multialélico: {:?}",
        record.alternates
    );
    let bases = check_allele_correspondence(&record, 100);
    assert_eq!(
        bases
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        2,
        "los dos alelos miden distinto, que es lo que hace discriminante la comparación"
    );
    fs::remove_dir_all(&root).ok();

    let (composite_root, composite) = publish(&composite_reference(), COMPOSITE_ALLELES);
    let composite_text = vcf_text(&composite);
    let composite_record = repeat_record(&composite_text);
    assert_eq!(
        composite_record.alternates.len(),
        2,
        "el fixture compuesto también tiene que dar un registro multialélico: {:?}",
        composite_record.alternates
    );
    assert_eq!(
        composite_record.info["RUS"], "AAAT,GGC,AAAT,GGC",
        "los bloques de cada alelo van en el orden de los motivos declarados"
    );
    assert_eq!(
        composite_record.info["RN"], "2,2",
        "dos bloques por alelo, que es lo que permite desalinear la lista dentro del alelo"
    );
    let motif_lengths = integers(&composite_record, "RUL");
    assert_eq!(
        motif_lengths,
        vec![4, 3, 4, 3],
        "los motivos miden distinto, de modo que desplazar la lista se nota"
    );
    let composite_bases = check_allele_correspondence(&composite_record, COMPOSITE_TRACT_BASES);
    assert_eq!(
        composite_bases.len(),
        4,
        "dos alelos de dos bloques cada uno"
    );
    assert_ne!(
        composite_bases[0] + composite_bases[1],
        composite_bases[2] + composite_bases[3],
        "los dos alelos miden distinto: {composite_bases:?}"
    );
    fs::remove_dir_all(&composite_root).ok();
}

/// **El oráculo semántico independiente acepta el VCF publicado y rechaza sus
/// corrupciones características.**
///
/// *Qué hace que pueda fallar:* no se afirma solo que el oráculo acepta —eso lo
/// cumpliría también un oráculo que no mirase estos campos—. Se corrompe el
/// documento publicado de quince formas distintas y se exige el rechazo de cada
/// una. **Cada guarda del oráculo tiene aquí su corrupción**: si una desapareciera
/// entera, alguna de estas filas dejaría de fallar.
#[test]
fn the_independent_oracle_accepts_the_vcf_and_rejects_its_characteristic_corruptions() {
    let (root, destination) = publish(&simple_reference(), SIMPLE_POLYMORPHIC);
    let fasta = destination.join("resources/references/ref-main.fa");
    let truth = destination.join("outputs/truth.jsonl");
    let vcf = destination.join("outputs/truth.vcf");
    validate_published_vcf_semantics(&fasta, &truth, &vcf).expect("el VCF publicado es compatible");

    let text = vcf_text(&destination);
    let record = repeat_record(&text);
    let counts = integers(&record, "RN");
    assert_eq!(counts, vec![1, 1], "un bloque por alelo en este fixture");
    let corrupted = [
        (
            "bases intercambiadas entre los dos alelos",
            "RB no son las bases que la verdad declara",
            text.replace(
                &format!("RB={}", record.info["RB"]),
                &format!("RB={}", reversed(&record.info["RB"])),
            ),
        ),
        (
            "copias intercambiadas entre los dos alelos",
            "RUC no son las bases entre la longitud del motivo",
            text.replace(
                &format!("RUC={}", record.info["RUC"]),
                &format!("RUC={}", reversed(&record.info["RUC"])),
            ),
        ),
        (
            "número de copias declarado entero",
            "declaración INFO de repetición distinta de la del estándar",
            text.replace(
                "##INFO=<ID=RUC,Number=.,Type=Float,",
                "##INFO=<ID=RUC,Number=.,Type=Integer,",
            ),
        ),
        (
            "recuento por alelo incompleto",
            "RN no trae exactamente un valor por alelo alternativo",
            text.replace(
                &format!("RN={}", record.info["RN"]),
                &format!("RN={}", counts[0]),
            ),
        ),
        (
            "motivo cambiado",
            "RUS no es el motivo que la verdad declara",
            text.replace(
                &format!("RUS={}", record.info["RUS"]),
                &format!("RUS={}", record.info["RUS"].replacen("CA", "CG", 1)),
            ),
        ),
        (
            "longitud de motivo cambiada",
            "RUL no es la longitud del motivo que la acompaña",
            text.replace(
                &format!("RUL={}", record.info["RUL"]),
                &format!("RUL={}", record.info["RUL"].replacen('2', "3", 1)),
            ),
        ),
        (
            "un registro sin tracto al que se le añaden campos de repetición",
            "un registro sin tracto publica campos de repetición",
            with_appended_info(&text, |line| !line.contains("RN="), ";RN=0"),
        ),
        (
            "un registro con tracto al que se le quita uno de los cinco campos",
            "un registro con tracto no publica los cinco campos de repetición",
            text.replace(&format!(";RB={}", record.info["RB"]), ""),
        ),
        (
            "clave repetida en INFO",
            "INFO repite una clave",
            with_appended_info(&text, |line| line.contains("RN="), ";RN=1,1"),
        ),
        (
            "declaración de repetición duplicada en la cabecera",
            "declaración INFO de repetición duplicada en la cabecera",
            text.replace(
                CANONICAL_DECLARATIONS[0],
                &format!(
                    "{}\n{}",
                    CANONICAL_DECLARATIONS[0], CANONICAL_DECLARATIONS[0]
                ),
            ),
        ),
        (
            "recuento por alelo que no es un entero",
            "RN trae un valor que no es un entero",
            text.replace(&format!("RN={}", record.info["RN"]), "RN=uno,uno"),
        ),
        (
            "número de copias que no es un decimal",
            "RUC trae un valor que no es un decimal",
            text.replace(&format!("RUC={}", record.info["RUC"]), "RUC=mucho,poco"),
        ),
        (
            "lista más corta que la suma de los recuentos",
            "las listas de repetición no suman los recuentos declarados en RN",
            text.replace(
                &format!("RUS={}", record.info["RUS"]),
                &format!(
                    "RUS={}",
                    record.info["RUS"].split(',').next().expect("primer motivo")
                ),
            ),
        ),
    ];
    reject_each(&fasta, &truth, &destination, &text, &corrupted);
    fs::remove_dir_all(&root).ok();

    // La cabecera que declara los campos sin que ningún registro los use
    // necesita un escenario que **no** publique ningún tracto.
    let (bare_root, bare) = publish(&simple_reference(), SIMPLE_REFERENCE_ONLY);
    let bare_text = vcf_text(&bare);
    assert!(
        !bare_text.contains("RN="),
        "este escenario no publica ningún campo de repetición"
    );
    let declared = [(
        "cabecera que declara los campos sin que ningún registro los use",
        "la cabecera no declara exactamente los campos de repetición",
        bare_text.replace(
            "##FILTER=<",
            &format!("{}\n##FILTER=<", CANONICAL_DECLARATIONS.join("\n")),
        ),
    )];
    reject_each(
        &bare.join("resources/references/ref-main.fa"),
        &bare.join("outputs/truth.jsonl"),
        &bare,
        &bare_text,
        &declared,
    );
    fs::remove_dir_all(&bare_root).ok();

    // El compuesto añade la corrupción que un locus de un solo motivo no puede
    // expresar: **los dos bloques de un mismo alelo intercambiados**.
    let (composite_root, composite) = publish(&composite_reference(), COMPOSITE_ALLELES);
    let composite_text = vcf_text(&composite);
    let composite_record = &repeat_records(&composite_text)[0];
    let swapped = [(
        "bloques intercambiados dentro del alelo",
        "RB no son las bases que la verdad declara",
        composite_text.replace(
            &format!("RB={}", composite_record.info["RB"]),
            &format!("RB={}", reversed(&composite_record.info["RB"])),
        ),
    )];
    reject_each(
        &composite.join("resources/references/ref-main.fa"),
        &composite.join("outputs/truth.jsonl"),
        &composite,
        &composite_text,
        &swapped,
    );
    fs::remove_dir_all(&composite_root).ok();
}

/// Añade texto al campo `INFO` de los registros que el filtro acepta.
///
/// Se edita la **columna octava**, no la línea entera: insertar el texto entre
/// columnas añadiría una columna y el documento lo rechazaría por otra razón,
/// que es exactamente lo que no queremos de una corrupción dirigida.
fn with_appended_info(text: &str, keep: impl Fn(&str) -> bool, suffix: &str) -> String {
    let mut output = String::new();
    for line in text.lines() {
        if line.starts_with('#') || !keep(line) {
            output.push_str(line);
        } else {
            let mut fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
            assert!(fields.len() > 7, "un registro VCF tiene columna INFO");
            fields[7].push_str(suffix);
            output.push_str(&fields.join("\t"));
        }
        output.push('\n');
    }
    output
}

/// Invierte una lista separada por comas.
fn reversed(value: &str) -> String {
    let mut items = value.split(',').collect::<Vec<_>>();
    items.reverse();
    items.join(",")
}

/// Escribe cada corrupción y exige que el oráculo la rechace **por su motivo**.
///
/// Comprobar solo que rechaza no basta: una corrupción que cayera por una
/// guarda ajena —por ejemplo por haber añadido una columna de más— pasaría la
/// prueba sin proteger a la guarda que dice proteger. Por eso se compara
/// también el diagnóstico.
fn reject_each(
    fasta: &Path,
    truth: &Path,
    destination: &Path,
    original: &str,
    corrupted: &[(&str, &str, String)],
) {
    for (reason, expected_detail, content) in corrupted {
        assert_ne!(
            content, original,
            "la corrupción «{reason}» tiene que cambiar el documento"
        );
        let path = destination.join("outputs/corrupted.vcf");
        fs::write(&path, content).expect("VCF corrompido");
        let error = validate_published_vcf_semantics(fasta, truth, &path)
            .expect_err("el oráculo tiene que rechazar la corrupción");
        let detail = format!("{error:?}");
        assert!(
            detail.contains(expected_detail),
            "«{reason}» se rechazó por otra razón: {detail}"
        );
    }
}

/// **El juez externo y el productor describen los mismos alelos; solo difieren
/// en la codificación.**
///
/// *Por qué esta prueba existe.* Al dejar de recortar los registros de
/// repetición, `bcftools norm` **sí** los recorta: la comparación «el
/// normalizador externo no altera el fichero», que la puerta DTC-066 hace
/// sobre los demás registros, deja de valer para estos, porque la
/// representación buscada no es la canónica. En su lugar se comprueba que el
/// juez externo **acepta** el registro, que lo **recorta de verdad**, y que
/// recortarlo **no cambia ninguna secuencia de alelo**, reconstruyendo las
/// secuencias contra la referencia.
///
/// *Lo que esta prueba NO garantiza, dicho sin adornos.* No dice que los
/// alelos publicados sean los que la verdad declara. No puede: su única
/// entrada es el fichero publicado, de modo que una corrupción presente en él
/// aparecería por igual en las dos formas y las dos seguirían coincidiendo.
/// Esa propiedad la garantiza el **oráculo semántico**, que reconstruye los
/// alelos desde la verdad y los compara uno a uno.
///
/// *Qué hace que pueda fallar:* se exige que `bcftools norm` **realinee** el
/// registro —si no lo realineara, la premisa sería falsa y no se estaría
/// comprobando nada—, que el `REF` resultante sea más corto, y que las
/// secuencias reconstruidas coincidan como **multiconjunto**, con el control
/// de que los dos alelos comparados son distintos entre sí.
///
/// *Aserción que no puede ejercitarse desde aquí, dicho sin adornos.* La
/// comparación de secuencias tiene por sujeto al **programa externo**: sus dos
/// lados salen del mismo fichero publicado, de modo que ninguna mutación de
/// este código puede hacerla fallar. Solo fallaría si `bcftools norm` reparase
/// o descartase un alelo, que es justo lo que existe para impedir. Lo que sí
/// muere si se mutan estas fuentes es la premisa: volver a recortar los
/// registros de repetición deja al normalizador sin nada que acortar y la
/// prueba cae.
#[test]
fn the_external_normalizer_only_trims_the_repeat_record_without_changing_its_alleles() {
    let version = Command::new("bcftools")
        .arg("--version")
        .output()
        .expect("bcftools instalado");
    assert_eq!(
        String::from_utf8_lossy(&version.stdout).lines().next(),
        Some("bcftools 1.24"),
        "DTC-066 exige exactamente esta versión"
    );

    let (root, destination) = publish(&simple_reference(), SIMPLE_POLYMORPHIC);
    let text = vcf_text(&destination);
    let record = repeat_record(&text);
    let normalized = Command::new("bcftools")
        .args([
            "norm",
            "--no-version",
            "--check-ref",
            "e",
            "--fasta-ref",
            destination
                .join("resources/references/ref-main.fa")
                .to_str()
                .expect("ruta UTF-8"),
            destination
                .join("outputs/truth.vcf")
                .to_str()
                .expect("ruta UTF-8"),
        ])
        .output()
        .expect("bcftools norm");
    assert!(
        normalized.status.success(),
        "bcftools norm falló: {}",
        String::from_utf8_lossy(&normalized.stderr)
    );
    // La premisa se lee del recuento que `norm` imprime, no de que la palabra
    // «realigned» aparezca: esa palabra está siempre en la cabecera del
    // recuento y buscarla no comprobaría nada.
    let report = String::from_utf8_lossy(&normalized.stderr).into_owned();
    let realigned = report
        .lines()
        .find_map(|line| line.rsplit_once('\t'))
        .and_then(|(_, counts)| counts.split('/').nth(3))
        .and_then(|value| value.trim().parse::<u64>().ok())
        .expect("bcftools norm publica su recuento de realineados");
    assert!(
        realigned >= 1,
        "la premisa de esta prueba es que el juez externo sí recorta: {report}"
    );

    let normalized_text = String::from_utf8(normalized.stdout).expect("salida UTF-8");
    let normalized_record = records(&normalized_text)
        .into_iter()
        .find(|other| other.position <= record.position && other.alternates.len() == 2)
        .expect("el registro de repetición normalizado");

    assert!(
        normalized_record.reference.len() < record.reference.len(),
        "el juez externo tiene que producir un REF más corto: {} frente a {}",
        normalized_record.reference.len(),
        record.reference.len()
    );

    // **Se reconstruye la secuencia entera**, no se comparan longitudes.
    //
    // Comparar la diferencia de longitud entre cada `ALT` y su `REF` no dice
    // nada: recortar la conserva siempre, para cualquier entrada, correcta o
    // corrupta. Lo que sí dice algo es aplicar cada registro a la referencia y
    // exigir que las secuencias resultantes sean **las mismas**: si el
    // normalizador hubiera reparado o descartado un alelo —que es lo que esta
    // comprobación existe para impedir—, las secuencias diferirían.
    let bases = fs::read_to_string(destination.join("resources/references/ref-main.fa"))
        .expect("referencia publicada")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<String>();
    let haplotypes = |record: &Record| {
        let start = record.position as usize - 1;
        let end = start + record.reference.len();
        assert_eq!(
            &bases[start..end],
            record.reference,
            "el REF del registro tiene que ser el de la referencia"
        );
        let mut sequences = record
            .alternates
            .iter()
            .map(|alternate| format!("{}{alternate}{}", &bases[..start], &bases[end..]))
            .collect::<Vec<_>>();
        // Multiconjunto ordenado, no conjunto: dos alelos iguales tienen que
        // seguir siendo dos.
        sequences.sort();
        sequences
    };
    let published = haplotypes(&record);
    let normalized = haplotypes(&normalized_record);
    assert_eq!(
        published, normalized,
        "recortar no puede cambiar ninguna secuencia de alelo"
    );
    assert_eq!(published.len(), 2, "control: se compararon dos alelos");
    assert_ne!(
        published[0], published[1],
        "control: los dos alelos son distintos, de modo que la igualdad no es trivial"
    );
    fs::remove_dir_all(&root).ok();
}

/// **Dos copias que llevan el mismo tracto con descomposiciones en bloques
/// distintas reciben la anotación derivada de sus bases.**
///
/// *Por qué existe esta prueba:* hasta el 2 de octubre de 2026 el oráculo
/// rechazaba esa verdad y el productor abortaba; desde entonces la anotación
/// de ese ALT se deriva de sus bases, porque una anotación describe lo que
/// está escrito (DTC-064). La producción no construye esta verdad —los alelos
/// de un locus comparten su lista de motivos—, de modo que hay que operarla a
/// mano.
///
/// *Cómo se construye:* se publica un escenario real, se toma un individuo
/// **homocigoto** —sus dos copias llevan el mismo tracto— y a una de ellas se
/// le cambia la descomposición por otra equivalente: las mismas bases
/// repartidas en dos bloques en vez de uno. El caso genotípico cataloga el
/// alelo nuevo y su portador se repunta, de modo que la verdad sigue siendo
/// **reconstruible**.
///
/// *Qué la hace discriminante:* el VCF publicado se acepta con la verdad
/// operada porque su anotación es la derivada, y el mismo VCF con la
/// descomposición de la cirugía —dos secuencias repetidas para ese ALT— se
/// rechaza.
#[test]
fn two_copies_carrying_the_same_tract_with_different_blocks_are_annotated_from_their_bases() {
    use genome_population_simulator::truth_schema::{
        AlleleCarrier, GenotypeCaseKind, TruthAllele, TruthCurrentPiece, TruthHistoryEvent,
        TruthRecord, TruthRepeatBlock,
    };

    let (root, destination) = publish(&simple_reference(), SIMPLE_POLYMORPHIC);
    let reference = validate_supplied_fasta(
        &destination.join("resources/references/ref-main.fa"),
        "ref-main",
        &[SequenceExpectation {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 300,
        }],
    )
    .expect("referencia publicada");
    let file = fs::File::open(destination.join("outputs/truth.jsonl")).expect("verdad");
    let truth = read_truth_jsonl(&mut std::io::BufReader::new(file)).expect("verdad legible");

    // Control: la **misma llamada** sobre la verdad intacta no produce este
    // rechazo. Un VCF vacío falla igualmente —le falta la cabecera de
    // columnas—, y lo que se fija aquí es que el motivo es **otro**: sin eso,
    // el rechazo de abajo podría no venir de la cirugía.
    let intact = validate_vcf_against_truth(
        reference.reference(),
        &truth,
        &mut std::io::Cursor::new(Vec::new()),
    )
    .expect_err("un VCF vacío no tiene cabecera de columnas");
    assert!(
        !format!("{intact:?}").contains("dos copias dan el mismo ALT con tractos distintos"),
        "la verdad intacta no puede producir el rechazo de la guarda: {intact:?}"
    );

    // Un individuo homocigoto: sus dos copias llevan el mismo tracto.
    let homozygous = truth
        .records
        .iter()
        .filter_map(|record| match record {
            TruthRecord::Genotype(genotype)
                if genotype.case_kind == GenotypeCaseKind::Repeat
                    && genotype.carriers.len() == 2
                    && genotype.carriers[0].allele_index == genotype.carriers[1].allele_index
                    && genotype.carriers[0].allele_index != 0 =>
            {
                Some(genotype.individual_id.clone())
            }
            _ => None,
        })
        .next()
        .expect("el escenario tiene un individuo homocigoto para el tracto");

    let mut surgical = truth.clone();
    let mut split_blocks = None;
    for record in &mut surgical.records {
        let TruthRecord::HaplotypeCopy(copy) = record else {
            continue;
        };
        if copy.identity.individual_id != homozygous || copy.identity.copy_index != 1 {
            continue;
        }
        // Las mismas bases en dos bloques en vez de uno: una descomposición
        // equivalente y distinta, que es exactamente lo que la guarda prohíbe
        // que convivan bajo el mismo ALT.
        for piece in &mut copy.current_pieces {
            if let TruthCurrentPiece::RepeatAllele { bases, blocks, .. } = piece {
                let split = vec![
                    TruthRepeatBlock {
                        motif: "CA".into(),
                        bases: 2,
                    },
                    TruthRepeatBlock {
                        motif: "CA".into(),
                        bases: bases.len() as u64 - 2,
                    },
                ];
                assert_ne!(*blocks, split, "la cirugía tiene que cambiar los bloques");
                *blocks = split.clone();
                split_blocks = Some((bases.clone(), split));
            }
        }
        for event in &mut copy.history {
            if let TruthHistoryEvent::RepeatAllele { blocks, .. } = event {
                *blocks = split_blocks.as_ref().expect("bloques nuevos").1.clone();
            }
        }
    }
    let (bases, split) = split_blocks.expect("la copia homocigota publica su tracto");

    for record in &mut surgical.records {
        let TruthRecord::Genotype(genotype) = record else {
            continue;
        };
        if genotype.individual_id != homozygous || genotype.case_kind != GenotypeCaseKind::Repeat {
            continue;
        }
        let locus_id = genotype
            .alleles
            .iter()
            .find_map(|allele| match allele {
                TruthAllele::Repeat { locus_id, .. } => Some(locus_id.clone()),
                _ => None,
            })
            .expect("el caso cataloga un alelo de repetición");
        genotype.alleles.push(TruthAllele::Repeat {
            locus_id,
            bases: bases.clone(),
            blocks: split.clone(),
        });
        let added = genotype.alleles.len() - 1;
        let previous = genotype.carriers[1].allele_index;
        genotype.carriers[1] = AlleleCarrier {
            allele_index: added,
            copy: genotype.carriers[1].copy.clone(),
        };
        genotype.allele_dosages[previous] -= 1;
        genotype.allele_dosages.push(1);
    }

    // Desde el 2 de octubre de 2026 (DTC-064) la anotación de ese ALT se
    // deriva de sus bases: un bloque de `CA` con todas ellas, que es lo que el
    // VCF publicado ya lleva. El VCF se acepta con la verdad intacta y con la
    // operada, y uno que anote ese ALT con la descomposición de la cirugía se
    // rechaza: el oráculo no elige una de las dos anotaciones declaradas.
    let vcf = vcf_text(&destination);
    for (label, truth) in [("intacta", &truth), ("operada", &surgical)] {
        validate_vcf_against_truth(
            reference.reference(),
            truth,
            &mut std::io::Cursor::new(vcf.as_bytes()),
        )
        .unwrap_or_else(|error| panic!("verdad {label}: {error:?}"));
    }
    let record = repeat_record(&vcf);
    let line = vcf
        .lines()
        .find(|line| line.starts_with(&format!("chr1\t{}\t", record.position)))
        .expect("la línea del registro de repetición")
        .to_owned();
    let index = record
        .alternates
        .iter()
        .position(|alternate| *alternate == bases)
        .expect("el ALT de la copia homocigota");
    let mut counts = integers(&record, "RN");
    counts[index] = 2;
    let split_line = line.replacen(
        &format!("RN={}", record.info["RN"]),
        &format!(
            "RN={}",
            counts
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        ),
        1,
    );
    assert_ne!(split_line, line, "control: la anotación cambió");
    let split_vcf = vcf.replace(&line, &split_line);
    validate_vcf_against_truth(
        reference.reference(),
        &surgical,
        &mut std::io::Cursor::new(split_vcf.as_bytes()),
    )
    .expect_err("la descomposición declarada por una sola copia no es la anotación del ALT");
    fs::remove_dir_all(&root).ok();
}

/// **Ningún registro usa representación simbólica, ni siquiera con un tracto de
/// diez mil bases.**
///
/// *Qué hace que pueda fallar:* la búsqueda negativa lleva su **control
/// positivo**, de modo que un patrón incapaz de encontrar nada no puede pasar
/// por una ausencia. Y el tracto mide diez mil bases, que es la cifra con la
/// que el estándar ilustra la forma simbólica: una proyección que recayera en
/// ella por encima de algún tamaño moriría aquí.
#[test]
fn no_record_uses_a_symbolic_representation_not_even_with_a_ten_thousand_base_tract() {
    let symbolic = |text: &str| {
        text.lines()
            .filter(|line| !line.starts_with('#'))
            .any(|line| line.split('\t').nth(4).is_some_and(|alt| alt.contains('<')))
    };
    let control = "chr1\t100\t.\tT\t<CNV:TR>\t.\tPASS\t.\tGT\t1|0\n";
    assert!(
        symbolic(control),
        "control positivo: la búsqueda sí encuentra un alelo simbólico cuando lo hay"
    );

    let (root, destination) = publish(&simple_reference(), TEN_THOUSAND_BASE_ALLELE);
    let text = vcf_text(&destination);
    let record = repeat_record(&text);
    assert_eq!(
        integers(&record, "RB"),
        vec![10_000],
        "el tracto publicado mide diez mil bases"
    );
    assert!(
        record.alternates[0].len() > 9_000,
        "y su ALT lleva la secuencia completa: {} bases",
        record.alternates[0].len()
    );
    assert!(
        !symbolic(&text),
        "ningún registro publicado usa representación simbólica"
    );
    validate_published_vcf_semantics(
        &destination.join("resources/references/ref-main.fa"),
        &destination.join("outputs/truth.jsonl"),
        &destination.join("outputs/truth.vcf"),
    )
    .expect("el tracto largo sigue siendo compatible con la verdad");
    fs::remove_dir_all(&root).ok();
}
