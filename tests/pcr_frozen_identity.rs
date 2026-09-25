//! **Ancla absoluta de la reproducibilidad del recorrido técnico.**
//!
//! # Por qué existe este fichero, dicho sin rodeos
//!
//! Toda la batería de la preparación amplificada comparaba **dos ejecuciones
//! del mismo binario**: la amplificada neutra contra la sin PCR. Cualquier
//! cambio que mueva las dos por igual le es invisible. Se midió: cambiar una
//! sola línea del muestreador —hacer que el intento cero use las claves de
//! reintento en vez de las heredadas— cambia **los cuatro FASTQ y los dos
//! manifiestos** de un escenario que ni siquiera declara amplificación, y las
//! 1404 pruebas seguían en verde.
//!
//! Lo que faltaba era un ancla **absoluta**: una huella congelada que no
//! dependa de comparar el binario consigo mismo. Es lo que hay aquí.
//!
//! # Las dos anclas, y por qué las dos
//!
//! 1. **La huella del FASTQ publicado** de un escenario pequeño **sin PCR**.
//!    Cierra el contrato entero: cubre el sorteo del inserto, la elección de
//!    copia, el inicio, la emisión del par, la quimera, la calidad y la
//!    escritura. Si algo de eso se mueve, muere.
//! 2. **Las secuencias de la primera pareja** producidas por el muestreador
//!    **por su camino real**, congeladas aparte. Señala el defecto con el dedo
//!    en vez de decir solo que algo cambió, y corre en milisegundos.
//!
//! # Qué hacer cuando esta prueba falle
//!
//! **No se actualiza el literal sin más.** Una huella que cambia significa una
//! de dos cosas: o el cambio reescribe todas las ejecuciones históricas del
//! proyecto —y entonces hay que decidirlo por escrito, no taparlo—, o es un
//! defecto. Actualizar el número para que calle es exactamente lo que este
//! fichero existe para impedir.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::technical_profile::{
    GPS_HISEQX_PCRFREE_PE150_V1, resolve_f4_technical_profile,
};
use genome_population_simulator::wgs_sampling::WgsPairedEndSampler;
use noodles_bgzf as bgzf;
use sha2::{Digest, Sha256};

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

// ---------------------------------------------------------------------------
// Los literales congelados.
//
// Se obtuvieron el 25 de septiembre de 2026 del árbol de esta entrega, y se
// comprobó que coinciden con los del commit anterior a ella: la preparación
// amplificada **no cambia un solo byte** de una ejecución que no la declara.
// ---------------------------------------------------------------------------

/// SHA-256 del FASTQ R1 publicado por [`FROZEN_SCENARIO`].
const FROZEN_R1_SHA256: &str = "098b4bf54744c0b6cbe8f7888a3ee3e3c1be2b97b1db7c9d68a85ffcff481f6b";
/// SHA-256 del FASTQ R2 publicado por [`FROZEN_SCENARIO`].
const FROZEN_R2_SHA256: &str = "4564c492605722641aa243f12db87b91fddbdf3b681d0c639c3b198a044013b1";
/// SHA-256 de las secuencias de la primera pareja que el muestreador produce
/// por su camino real, sin pasar por la publicación.
const FROZEN_FIRST_PAIR_SHA256: &str =
    "4e1b056bc5c1bf0f150abe5ac4a0fc1312697882f8f97904ccfd7e19a392fa2e";

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-anclaje-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Bases equilibradas y sin periodo corto, fijadas por un generador
/// congruencial cuyo estado inicial es parte del ancla.
fn frozen_bases(length: usize) -> String {
    let mut state: u64 = 0x1234_5678_9ABC_DEF0;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            match (state >> 33) % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            }
        })
        .collect()
}

/// **El escenario congelado: sin amplificación, sin microsatélites y sin
/// artefactos.** Es deliberadamente lo más parecido a lo que el proyecto lleva
/// publicando desde antes de esta entrega.
const FROZEN_LENGTH: usize = 3_000;

fn publish_frozen_scenario() -> (PathBuf, PathBuf) {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let sequence = frozen_bases(FROZEN_LENGTH);
    let reference_base = sequence.as_bytes()[19] as char;
    let alternate = if reference_base == 'T' { 'G' } else { 'T' };
    fs::write(root.join("reference.fa"), format!(">chr1\n{sequence}\n")).expect("FASTA");
    fs::write(
        root.join("population-snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\n\
             chr1\t20\t{reference_base}\t{alternate}\tP1\t0.5\n"
        ),
    )
    .expect("matriz");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"
ploidy = 1

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
length = {FROZEN_LENGTH}

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 4

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
    (root, destination)
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn read_sequences(path: &Path) -> Vec<Vec<u8>> {
    let mut text = Vec::new();
    bgzf::io::Reader::new(fs::File::open(path).expect("FASTQ"))
        .read_to_end(&mut text)
        .expect("legible");
    String::from_utf8(text)
        .expect("ASCII")
        .lines()
        .skip(1)
        .step_by(4)
        .map(|line| line.as_bytes().to_vec())
        .collect()
}

/// **La huella del FASTQ publicado es la congelada.**
///
/// Esta es la prueba que muere si el muestreador cambia de clave aleatoria, si
/// el sorteo del inserto consume distinto, si la emisión del par se desplaza o
/// si la escritura cambia de orden. No compara el binario consigo mismo: el
/// literal es un ancla absoluta.
#[test]
fn the_unamplified_scenario_keeps_its_frozen_fastq_digest() {
    let (root, destination) = publish_frozen_scenario();
    let first =
        fs::read(destination.join("outputs/fastq/reads-main/sample-A_R1.fastq.gz")).expect("R1");
    let second =
        fs::read(destination.join("outputs/fastq/reads-main/sample-A_R2.fastq.gz")).expect("R2");
    // Control positivo: el escenario produjo lecturas de verdad, de modo que
    // la huella no es la de un fichero vacío.
    let sequences =
        read_sequences(&destination.join("outputs/fastq/reads-main/sample-A_R1.fastq.gz"));
    assert!(
        sequences.len() > 20,
        "el escenario congelado tiene que publicar lecturas: {}",
        sequences.len()
    );
    let observed_r1 = digest(&first);
    let observed_r2 = digest(&second);
    let _ = fs::remove_dir_all(&root);
    assert_eq!(
        observed_r1, FROZEN_R1_SHA256,
        "el R1 del escenario sin amplificar ha cambiado de bytes; NO se actualiza el literal \
         sin decidir por escrito que todas las ejecuciones históricas se reescriben"
    );
    assert_eq!(
        observed_r2, FROZEN_R2_SHA256,
        "el R2 del escenario sin amplificar ha cambiado de bytes"
    );
}

/// **La primera pareja que el muestreador produce por su camino real es la
/// congelada.**
///
/// Señala el defecto con el dedo: si el muestreador eligiera otra clave para el
/// intento cero, o consumiera distinto del sorteo del inserto, esta pareja
/// cambia. Corre en milisegundos y no depende de la publicación.
#[test]
fn the_first_pair_the_sampler_produces_keeps_its_frozen_digest() {
    let bases: Vec<u8> = frozen_bases(2_000).into_bytes();
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "test").expect("secuencia"),
        FastaHeader::parse(b"chr1").expect("cabecera"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(&bases).expect("ADN"),
    )
    .expect("secuencia de referencia");
    let reference = Reference::new(
        ReferenceId::parse("reference-main", "test").expect("referencia"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("referencia");
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individuo"),
        &reference,
        Ploidy::new(1).expect("ploidía"),
    );
    let identities = derive_technical_identity_plan(
        RootSeed::parse("root", SEED).expect("semilla"),
        reference.id(),
        EventId::parse("wgs-pe", "test").expect("evento"),
        &[individual.id().clone()],
    )
    .expect("identidades");
    let profile =
        resolve_f4_technical_profile(GPS_HISEQX_PCRFREE_PE150_V1).expect("perfil vigente");
    let mut sampler = WgsPairedEndSampler::new_for_pair_count(
        &reference,
        &individual,
        &identities,
        profile,
        100,
        1,
    )
    .expect("muestreador");
    let pair = sampler
        .next()
        .expect("una pareja")
        .expect("la pareja se produce");
    let mut material = Vec::new();
    material.extend_from_slice(pair.r1().sequence());
    material.push(b'|');
    material.extend_from_slice(pair.r2().sequence());
    // Control positivo: las dos lecturas tienen la longitud pedida, de modo que
    // la huella no es la de un material vacío o truncado.
    assert_eq!(pair.r1().sequence().len(), 100);
    assert_eq!(pair.r2().sequence().len(), 100);
    assert_eq!(
        digest(&material),
        FROZEN_FIRST_PAIR_SHA256,
        "la primera pareja del muestreador ha cambiado: o cambió la clave aleatoria del intento \
         cero, o cambió lo que se consume del sorteo del inserto"
    );
}
