//! DT-013 y S-03: **con referencia simulada, el tracto de cada microsatélite
//! declarado queda colocado en la referencia** antes del relleno al azar.
//!
//! Hasta el 01/10/2026 la preparación de la referencia llamaba al generador sin
//! regiones. El FASTA publicado llevaba bases al azar en el tramo del locus,
//! todas las copias recibían un alelo alternativo y el alelo de la referencia
//! no lo podía tener nadie. Estas pruebas fijan lo contrario, y cada una lo
//! hace **leyendo lo publicado**: el FASTA como texto, el VCF por columnas y la
//! verdad como JSON, sin llamar a ninguna función del productor.
//!
//! # Cómo están elegidos los fixtures, que es lo que hace que puedan fallar
//!
//! - **El tramo se compara con el motivo repetido escrito aquí**, no con lo que
//!   el productor diga que colocó. Con bases al azar la igualdad es imposible.
//! - **El tramo no múltiplo del motivo** mide 30 bases con un motivo de cuatro:
//!   siete copias enteras y media. Un relleno que se detuviera en la última
//!   copia entera, o que redondeara hacia arriba, deja otra cosa en las dos
//!   últimas bases o en la primera del flanco derecho.
//! - **El alelo que coincide con la referencia** es un alelo declarado de
//!   quince unidades sobre un tramo de treinta bases de `CA`. Con el tracto al
//!   azar ese alelo era alternativo, como todos; con el tracto colocado es la
//!   referencia. El fixture lo fuerza de dos maneras: con frecuencia uno —todas
//!   las copias lo llevan, y no puede quedar ni un registro del locus— y con un
//!   reparto mitad y mitad entre doce copias cuya mezcla se exige en la prueba,
//!   para que la propiedad no dependa de que el sorteo lo deje ver.
//! - **El generador por deslizamiento** usa `AGAT` sobre 48 bases, el caso de
//!   DT-013. Con un motivo de dos bases el defecto no siempre se veía como
//!   rechazo: unas bases al azar contienen casi siempre algún `CA`, y el
//!   generador partía entonces de un tramo al azar con pureza baja. Con cuatro
//!   bases, la ventana al azar de este fixture no contiene ninguna copia y el
//!   defecto se ve como el rechazo de DT-013. Las dos formas las caza la misma
//!   aserción: el alelo ancestral tiene que ser el tracto colocado.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::validate_published_resources;
use serde_json::Value;
use sha2::{Digest, Sha256};

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-sim-ref-tract-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Documento general v4 con referencia simulada de dos cromosomas de mil
/// bases, una variante de población confinada a las cincuenta primeras bases
/// de `chr1` —lejos de todo locus, porque una variante dentro de un tracto
/// rechaza el escenario por la decisión S-32— y VCF de texto.
fn scenario(loci: &str, individuals: usize) -> String {
    let ids = (1..=individuals)
        .map(|index| format!("\"I{index:02}\""))
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "weighted_regions"
outside_weight = 0.0
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.regions]]
contig = "chr1"
start = 1
end = 50
weight = 1.0

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "simulated"
id = "ref-main"

[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1000

[[reference.sequences]]
id = "chr2"
class = "chromosome"
length = 1000
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = {individuals}
individual_ids = [{ids}]
ploidy = 2
genotype_model = "hwe"

[vcf_output]
mode = "text"
"#
    )
}

/// Ejecuta el documento, valida los recursos publicados con el juez
/// independiente y devuelve el directorio publicado.
fn publish(document: &str) -> PathBuf {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    let config = root.join("general.toml");
    fs::write(&config, document).expect("configuración");
    let effective = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective, &destination).expect("publicación");
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    validate_published_resources(
        &destination,
        &effective_path,
        &fs::read(&effective_path).expect("configuración efectiva"),
        &manifest_path,
        &fs::read(&manifest_path).expect("manifiesto"),
        None,
    )
    .expect("los recursos publicados se validan juntos");
    destination
}

/// El FASTA publicado, leído como texto: cabecera y bases de cada secuencia.
fn published_fasta(destination: &Path) -> BTreeMap<String, String> {
    let text = fs::read_to_string(destination.join("resources/references/ref-main.fa"))
        .expect("FASTA publicado");
    let mut sequences = BTreeMap::new();
    let mut current: Option<String> = None;
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let id = header.split_whitespace().next().expect("id").to_owned();
            sequences.insert(id.clone(), String::new());
            current = Some(id);
        } else {
            sequences
                .get_mut(current.as_ref().expect("bases antes de la cabecera"))
                .expect("secuencia")
                .push_str(line);
        }
    }
    sequences
}

/// Las bases del tramo de base uno e inclusivo `start..=end`.
fn span(sequence: &str, start: usize, end: usize) -> &str {
    &sequence[start - 1..end]
}

/// Los registros del VCF con la marca de repetición, por columnas.
fn repeat_records(destination: &Path) -> Vec<Vec<String>> {
    fs::read_to_string(destination.join("outputs/truth.vcf"))
        .expect("VCF publicado")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_owned).collect::<Vec<_>>())
        .filter(|fields| fields[7].split(';').any(|field| field.starts_with("RN=")))
        .collect()
}

/// Los registros de genotipo de repetición de `truth.jsonl`.
fn repeat_genotypes(destination: &Path) -> Vec<Value> {
    fs::read_to_string(destination.join("outputs/truth.jsonl"))
        .expect("verdad publicada")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("línea JSON"))
        .filter(|record| record["record_type"] == "genotype" && record["case_kind"] == "repeat")
        .collect()
}

fn effective_locus(destination: &Path) -> Value {
    let effective: Value = serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("configuración efectiva"),
    )
    .expect("JSON");
    effective["repeat_loci"][0].clone()
}

const CA_LOCUS_HEADER: &str = r#"
[[reference.repeat_loci]]
id = "STR-2"
sequence = "chr2"
start = 701
end = 730
motif = "CA"
"#;

/// El caso de la revisión adversarial: `CA` en `chr2:701-730`.
///
/// Debe fallar contra: una preparación que no pase las regiones al generador
/// —el tramo sale al azar—, y una que las coloque una base más allá —el tramo
/// empieza por una base al azar y sigue `CACA…`—. Las dos mutaciones se
/// ejecutaron, cada una en el productor y en el juez a la vez: mutando solo uno
/// de los dos, lo que muere es la publicación, porque el juez que regenera la
/// referencia ya no coincide con el FASTA.
#[test]
fn the_published_fasta_carries_the_motif_over_the_declared_span() {
    let loci = format!(
        "{CA_LOCUS_HEADER}
[[reference.repeat_loci.alleles]]
units = 15
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 17
frequency = 0.5
"
    );
    let destination = publish(&scenario(&loci, 2));
    let fasta = published_fasta(&destination);
    assert_eq!(span(&fasta["chr2"], 701, 730), "CA".repeat(15));
}

/// Un tramo de 30 bases con un motivo de cuatro: siete copias enteras y la
/// mitad de la octava, que es `AG`.
///
/// Debe fallar contra: un relleno que se detenga en la última copia entera,
/// que deja las dos últimas bases al azar. Se ejecutó.
#[test]
fn a_span_that_is_not_a_multiple_of_the_motif_ends_in_a_partial_copy() {
    let loci = r#"
[[reference.repeat_loci]]
id = "STR-3"
sequence = "chr2"
start = 101
end = 130
motif = "AGAT"

[[reference.repeat_loci.alleles]]
units = 8
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 6
frequency = 0.5
"#;
    let destination = publish(&scenario(loci, 2));
    let fasta = published_fasta(&destination);
    assert_eq!(
        span(&fasta["chr2"], 101, 130),
        format!("{}AG", "AGAT".repeat(7))
    );
}

/// **Con frecuencia uno, todas las copias llevan el alelo de la referencia** y
/// no puede quedar nada del locus ni en el VCF ni en la verdad.
///
/// Debe fallar contra: el tracto al azar, con el que ese alelo es alternativo
/// para todas las copias. Se ejecutó.
#[test]
fn an_allele_equal_to_the_placed_tract_leaves_every_copy_on_the_reference() {
    let loci = format!(
        "{CA_LOCUS_HEADER}
[[reference.repeat_loci.alleles]]
units = 15
frequency = 1.0
"
    );
    let destination = publish(&scenario(&loci, 3));
    assert!(
        repeat_records(&destination).is_empty(),
        "un locus cuyo único alelo es el de la referencia no publica registro"
    );
    assert!(
        repeat_genotypes(&destination).is_empty(),
        "ni caso de genotipo en la verdad"
    );
}

/// **Una copia cuyo alelo coincide con el tracto se publica como referencia**:
/// `0` en el VCF y alelo `reference` en la verdad, mientras la otra mitad de
/// las copias lleva el alternativo.
///
/// Se exige la mezcla —alguna copia en `0` y alguna en `1`— para que la prueba
/// no pase en vacío si el sorteo dejara todas en el mismo alelo.
///
/// Debe fallar contra: el tracto al azar, con el que los dos alelos declarados
/// son alternativos, el registro lleva dos `ALT` y ninguna copia es `0`. Se
/// ejecutó.
#[test]
fn a_copy_carrying_the_reference_tract_is_published_as_reference() {
    let loci = format!(
        "{CA_LOCUS_HEADER}
[[reference.repeat_loci.alleles]]
units = 15
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 17
frequency = 0.5
"
    );
    let destination = publish(&scenario(&loci, 6));
    let records = repeat_records(&destination);
    assert_eq!(records.len(), 1, "un registro del locus");
    let record = &records[0];
    assert_eq!(record[1], "701");
    assert_eq!(record[3], "CA".repeat(15), "REF es el tracto colocado");
    assert_eq!(
        record[4],
        "CA".repeat(17),
        "un único ALT, el de 17 unidades"
    );
    let calls = record[9..]
        .iter()
        .flat_map(|genotype| genotype.split('|').map(str::to_owned))
        .collect::<Vec<_>>();
    assert_eq!(calls.len(), 12, "seis diploides");
    assert!(
        calls.iter().any(|call| call == "0"),
        "alguna copia en 0: {calls:?}"
    );
    assert!(
        calls.iter().any(|call| call == "1"),
        "alguna copia en 1: {calls:?}"
    );
    assert!(calls.iter().all(|call| call == "0" || call == "1"));

    // La verdad solo publica genotipo de repetición para el individuo que
    // lleva algún alternativo; se compara individuo a individuo con su columna
    // del VCF, y se exige algún heterocigoto para que la dosis de referencia
    // de la verdad no se compare en vacío.
    let samples = vcf_samples(&destination);
    let by_individual = samples
        .iter()
        .zip(&record[9..])
        .map(|(sample, genotype)| (sample.clone(), genotype.clone()))
        .collect::<BTreeMap<_, _>>();
    let genotypes = repeat_genotypes(&destination);
    let mut heterozygotes = 0;
    for genotype in &genotypes {
        let alleles = genotype["alleles"].as_array().expect("alelos");
        assert_eq!(alleles[0]["kind"], "reference");
        for allele in &alleles[1..] {
            assert_eq!(allele["kind"], "repeat");
            assert_ne!(
                allele["bases"],
                Value::from("CA".repeat(15)),
                "el alelo igual al tracto no puede publicarse como alternativo"
            );
        }
        let individual = genotype["individual_id"].as_str().expect("individuo");
        let zeros = by_individual[individual]
            .split('|')
            .filter(|call| *call == "0")
            .count() as u64;
        assert_eq!(
            genotype["allele_dosages"][0].as_u64().expect("dosis"),
            zeros,
            "la verdad de {individual} cuenta como referencia las mismas copias que el VCF"
        );
        if zeros == 1 {
            heterozygotes += 1;
        }
    }
    assert!(heterozygotes > 0, "hace falta algún heterocigoto");
    let carriers = by_individual
        .values()
        .filter(|genotype| genotype.contains('1'))
        .count();
    assert_eq!(genotypes.len(), carriers, "un genotipo por portador");
}

/// Los nombres de las muestras, de la cabecera `#CHROM` del VCF.
fn vcf_samples(destination: &Path) -> Vec<String> {
    let text = fs::read_to_string(destination.join("outputs/truth.vcf")).expect("VCF publicado");
    let header = text
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("cabecera de columnas");
    header.split('\t').skip(9).map(str::to_owned).collect()
}

const AGAT_SLIPPAGE: &str = r#"
[[reference.repeat_loci]]
id = "STR-2"
sequence = "chr2"
start = 501
end = 548
motif = "AGAT"

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = MEAN
"#;

/// **El generador por deslizamiento funciona sobre referencia simulada** y
/// parte del tracto colocado (DT-013).
///
/// Con media cero ningún linaje se mueve y el abanico es exactamente el alelo
/// ancestral: es la forma más directa de leer cuál fue. Con media tres el
/// abanico se reparte alrededor de él, y se exige que lo haga en unidades
/// enteras de `AGAT` puro, que es lo que sale de partir de `(AGAT)12`.
///
/// Debe fallar contra: el tracto al azar, con el que este escenario se rechaza
/// con «no hay alelo ancestral del que partir». Se ejecutó.
#[test]
fn the_slippage_generator_starts_from_the_placed_tract() {
    let destination = publish(&scenario(&AGAT_SLIPPAGE.replace("MEAN", "0.0"), 2));
    let fasta = published_fasta(&destination);
    assert_eq!(span(&fasta["chr2"], 501, 548), "AGAT".repeat(12));
    let alleles = effective_locus(&destination)["alleles"].clone();
    let alleles = alleles.as_array().expect("alelos");
    assert_eq!(
        alleles.len(),
        1,
        "con media cero el abanico es el ancestral"
    );
    assert_eq!(alleles[0]["sequence"], Value::from("AGAT".repeat(12)));
    assert_eq!(alleles[0]["origin"], "generated");

    let destination = publish(&scenario(&AGAT_SLIPPAGE.replace("MEAN", "3.0"), 2));
    let alleles = effective_locus(&destination)["alleles"].clone();
    let alleles = alleles.as_array().expect("alelos");
    assert!(alleles.len() > 1, "con media tres el abanico se abre");
    for allele in alleles {
        let units = allele["units"].as_f64().expect("unidades");
        assert_eq!(units.fract(), 0.0);
        assert_eq!(
            allele["sequence"],
            Value::from("AGAT".repeat(units as usize)),
            "cada alelo es el ancestral puro con unidades enteras de más o de menos"
        );
    }
}

/// **Una referencia simulada sin microsatélites sale idéntica byte a byte a la
/// de antes de la corrección.**
///
/// La huella se tomó con el binario del commit `7dd1fbe`, anterior a la
/// corrección, compilado aparte desde una copia extraída con `git archive`, y
/// ejecutando este mismo documento. No es una comparación entre dos
/// ejecuciones del binario de hoy, que sería una simetría y no una identidad.
/// El fixture lleva `N` y dos secuencias para recorrer también el sorteo de
/// bases desconocidas y el paso del flujo de azar de una secuencia a la
/// siguiente.
///
/// Debe fallar contra: una preparación que pase una región aunque no haya
/// loci —una sola base al principio de `chr1`—, que cambia el consumo de azar
/// de todo lo que viene detrás. Se ejecutó.
#[test]
fn without_loci_the_simulated_reference_matches_the_digest_frozen_at_7dd1fbe() {
    const FROZEN_DIGEST: &str = "8ac22b01ebfbe643e8239762fe1b23b8dab552cfa3127fdabab5ac96195c445e";
    let document = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "generate"

[population_snvs.generate]
count_mode = "exact"
locus_count = 2
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8

[reference]
kind = "simulated"
id = "ref-main"

[reference.composition]
A = 0.3
C = 0.2
G = 0.2
T = 0.3
N = 0.05

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1000

[[reference.sequences]]
id = "chr2"
class = "chromosome"
length = 1000

[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"
"#
    );
    let destination = publish(&document);
    let bytes =
        fs::read(destination.join("resources/references/ref-main.fa")).expect("FASTA publicado");
    let digest = Sha256::digest(&bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    assert_eq!(digest, FROZEN_DIGEST);
}

// ---------------------------------------------------------------------------
// DT-016: los microsatélites compuestos, con `reference_units`.
// ---------------------------------------------------------------------------

/// Locus compuesto **saturado** `(CAG)4(CA)9` sobre `chr2:101-130`: los dos
/// bloques tienen motivo, longitud de motivo, unidades y bases distintos, de
/// modo que invertir los bloques, cruzar las unidades o repetir el primer
/// motivo en todos deja otras bases en el tramo.
fn composite_locus(reference_units: &str) -> String {
    format!(
        r#"
[[reference.repeat_loci]]
id = "STR-C"
sequence = "chr2"
start = 101
end = 130
motifs = ["CAG", "CA"]
{reference_units}

[[reference.repeat_loci.alleles]]
units = [4, 9]
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = [4, 10]
frequency = 0.5
"#
    )
}

const COMPOSITE_TRACT: &str = "CAGCAGCAGCAGCACACACACACACACACA";

/// **Un compuesto se coloca bloque a bloque**, en el orden de sus motivos, y el
/// alelo que coincide con él es la referencia del registro.
///
/// Debe fallar contra: una preparación que ignore el campo —el tramo sale al
/// azar—, una que invierta los bloques —`CA` delante de `CAG`— y una que
/// cruce las unidades entre bloques. Las tres se ejecutaron en el productor y
/// en el juez a la vez. Las dos primeras mueren en la aserción del tramo; la
/// tercera muere antes, en la publicación, porque con este fixture cruzar las
/// unidades cambia la suma —`(CAG)9(CA)4` mide 35— y el juez lo rechaza por no
/// caber en el tramo. Mutando solo el productor, las tres mueren en la
/// publicación, porque el juez coloca los bloques por su cuenta.
#[test]
fn a_composite_locus_is_placed_block_by_block() {
    assert_eq!(COMPOSITE_TRACT, "CAG".repeat(4) + &"CA".repeat(9));
    let destination = publish(&scenario(&composite_locus("reference_units = [4, 9]"), 6));
    let fasta = published_fasta(&destination);
    assert_eq!(span(&fasta["chr2"], 101, 130), COMPOSITE_TRACT);
    let records = repeat_records(&destination);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0][3], COMPOSITE_TRACT, "REF es el tracto colocado");
    assert_eq!(
        records[0][4],
        "CAG".repeat(4) + &"CA".repeat(10),
        "el único ALT es el alelo que no coincide"
    );
    assert_eq!(
        effective_locus(&destination)["reference_units"],
        serde_json::json!([4, 9]),
        "la configuración efectiva publica las unidades por bloque"
    );
}

/// El rechazo de configuración, con su campo y su detalle.
fn rejection(document: &str) -> (String, String) {
    let root = temporary("rejected");
    fs::create_dir_all(&root).expect("directorio");
    fs::write(root.join("reference.fa"), supplied_fasta()).expect("FASTA");
    let config = root.join("general.toml");
    fs::write(&config, document).expect("configuración");
    let error = load_and_validate_config(&config).expect_err("la configuración se rechaza");
    (
        error.context().get("fields").cloned().unwrap_or_default(),
        error.context().get("detail").cloned().unwrap_or_default(),
    )
}

/// Dos cromosomas de mil bases con `(CAG)4(CA)9` en `chr2:101-130`, para que
/// la referencia suministrada lleve el tracto y el único motivo de rechazo
/// posible sea el campo.
fn supplied_fasta() -> String {
    let chr2 = "G".repeat(100) + COMPOSITE_TRACT + &"G".repeat(870);
    format!(">chr1\n{}\n>chr2\n{chr2}\n", "G".repeat(1000))
}

fn supplied(document: &str) -> String {
    document.replace(
        "kind = \"simulated\"\nid = \"ref-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n",
        "kind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n",
    )
}

const FIELD: &str = "reference.repeat_loci.reference_units";

/// **Unidades que no cuadran con el tramo se rechazan antes de simular**,
/// nombrando el locus, la suma y el tramo. Junto al rechazo va su pareja
/// válida, para que una lectura que rechazara siempre no pase por buena.
///
/// Debe fallar contra: una resolución que no compare la suma con el tramo, que
/// deja pasar 28 bases en un tramo de 30 —las dos últimas al azar—, y una que
/// solo rechace las sumas cortas, que deja pasar 32. Las dos se ejecutaron.
#[test]
fn reference_units_that_do_not_fill_the_span_are_rejected() {
    load_and_validate_config(&write(&scenario(
        &composite_locus("reference_units = [4, 9]"),
        2,
    )))
    .expect("la pareja válida se admite");
    let (fields, detail) = rejection(&scenario(&composite_locus("reference_units = [4, 8]"), 2));
    assert_eq!(fields, FIELD);
    assert!(detail.contains("locus STR-C"), "{detail}");
    assert!(detail.contains("suma 28 bases"), "{detail}");
    assert!(detail.contains("101..130 mide 30"), "{detail}");

    // Y la que se pasa: 32 bases en un tramo de 30. Sin ella, una resolución
    // que solo rechazara las sumas cortas pasaba.
    let (fields, detail) = rejection(&scenario(&composite_locus("reference_units = [4, 10]"), 2));
    assert_eq!(fields, FIELD);
    assert!(detail.contains("suma 32 bases"), "{detail}");
    assert!(detail.contains("101..130 mide 30"), "{detail}");

    let (fields, detail) = rejection(&scenario(&composite_locus("reference_units = [4]"), 2));
    assert_eq!(fields, FIELD);
    assert!(detail.contains("hace falta uno por bloque"), "{detail}");

    let (fields, detail) = rejection(&scenario(&composite_locus("reference_units = [0, 15]"), 2));
    assert_eq!(fields, FIELD);
    assert!(detail.contains("cero unidades en el bloque 0"), "{detail}");
}

/// **El campo en un locus de un solo motivo se rechaza**: su tramo ya fija el
/// tracto, y aceptarlo en silencio haría creer que cambia algo.
///
/// Debe fallar contra: una resolución que solo mire el campo en los
/// compuestos. Se ejecutó.
#[test]
fn reference_units_on_a_single_motif_locus_are_rejected() {
    let loci = format!(
        "{CA_LOCUS_HEADER}reference_units = [15]

[[reference.repeat_loci.alleles]]
units = 15
frequency = 1.0
"
    );
    let (fields, detail) = rejection(&scenario(&loci, 2));
    assert_eq!(fields, FIELD);
    assert!(detail.contains("locus STR-2"), "{detail}");
    assert!(
        detail.contains("no se declara en un locus de un solo motivo con tramo"),
        "{detail}"
    );
}

/// **El campo con referencia suministrada se rechaza**: allí el tracto se lee
/// del FASTA. El FASTA del fixture lleva el tracto, y sin el campo el mismo
/// documento se admite, de modo que el único motivo del rechazo es el campo.
///
/// Debe fallar contra: una resolución que no mire el tipo de referencia. Se
/// ejecutó.
#[test]
fn reference_units_with_a_supplied_reference_are_rejected() {
    let without = supplied(&scenario(&composite_locus(""), 2));
    assert!(without.contains("kind = \"supplied\""), "el fixture cambió");
    let root = temporary("supplied");
    fs::create_dir_all(&root).expect("directorio");
    fs::write(root.join("reference.fa"), supplied_fasta()).expect("FASTA");
    fs::write(root.join("general.toml"), &without).expect("configuración");
    load_and_validate_config(&root.join("general.toml"))
        .expect("sin el campo, el compuesto sobre referencia suministrada se admite");

    let (fields, detail) = rejection(&supplied(&scenario(
        &composite_locus("reference_units = [4, 9]"),
        2,
    )));
    assert_eq!(fields, FIELD);
    assert!(detail.contains("locus STR-C"), "{detail}");
    assert!(
        detail.contains("solo se declara con referencia simulada"),
        "{detail}"
    );
}

/// **Un compuesto sin el campo sobre referencia simulada se rechaza**: sin las
/// unidades por bloque no hay forma de construir su tracto, y hasta el
/// 01/10/2026 el tramo salía al azar sin avisar (DT-016).
///
/// Debe fallar contra: una resolución que admita el compuesto sin unidades. Se
/// ejecutó.
#[test]
fn a_composite_locus_on_a_simulated_reference_needs_reference_units() {
    let (fields, detail) = rejection(&scenario(&composite_locus(""), 2));
    assert_eq!(fields, FIELD);
    assert!(detail.contains("locus STR-C"), "{detail}");
    assert!(detail.contains("declare reference_units"), "{detail}");
}

/// **Ida y vuelta por la configuración efectiva**: el juez que regenera la
/// referencia simulada lee los bloques de `effective-config.json`, de modo que
/// un documento que publique otras unidades que las usadas, o que no las
/// publique, se rechaza.
///
/// La primera corrupción mantiene la suma —`(CAG)6(CA)6` también mide 30— para
/// que el rechazo no pueda venir de una comprobación de longitud: solo lo da la
/// comparación base a base con el FASTA publicado.
///
/// Debe fallar contra: un juez que coloque los bloques desde la configuración
/// reconstruida del TOML en vez de desde la efectiva publicada, que no ve
/// ninguna de las dos corrupciones. Se ejecutó.
#[test]
fn the_reference_judge_reads_the_blocks_from_the_effective_config() {
    let destination = publish(&scenario(&composite_locus("reference_units = [4, 9]"), 2));
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let manifest = fs::read(&manifest_path).expect("manifiesto");
    let original: Value =
        serde_json::from_slice(&fs::read(&effective_path).expect("configuración efectiva"))
            .expect("JSON");
    let judge = |document: &Value| {
        validate_published_resources(
            &destination,
            &effective_path,
            &serde_json::to_vec_pretty(document).expect("JSON"),
            &manifest_path,
            &manifest,
            None,
        )
    };
    judge(&original).expect("el documento reserializado se acepta");

    let mut moved = original.clone();
    moved["repeat_loci"][0]["reference_units"] = serde_json::json!([6, 6]);
    let error = judge(&moved).expect_err("otras unidades con la misma suma");
    let detail = error.context().get("detail").cloned().unwrap_or_default();
    assert!(detail.contains("no coincide base a base"), "{detail}");

    let mut dropped = original.clone();
    dropped["repeat_loci"][0]
        .as_object_mut()
        .expect("locus")
        .remove("reference_units");
    let error = judge(&dropped).expect_err("el campo borrado");
    let detail = error.context().get("detail").cloned().unwrap_or_default();
    assert!(detail.contains("locus STR-C"), "{detail}");
    assert!(detail.contains("no publica reference_units"), "{detail}");
}

fn write(document: &str) -> PathBuf {
    let root = temporary("valid");
    fs::create_dir_all(&root).expect("directorio");
    let path = root.join("general.toml");
    fs::write(&path, document).expect("configuración");
    path
}

// ---------------------------------------------------------------------------
// Correcciones de la revisión adversarial del 01/10/2026.
// ---------------------------------------------------------------------------

/// Dos loci en **secuencias distintas**, declarados en orden **inverso** al de
/// las secuencias: el compuesto de `chr1` después del simple de `chr2`.
const TWO_LOCI: &str = r#"
[[reference.repeat_loci]]
id = "STR-2"
sequence = "chr2"
start = 701
end = 730
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 15
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 17
frequency = 0.5

[[reference.repeat_loci]]
id = "STR-C"
sequence = "chr1"
start = 101
end = 130
motifs = ["CAG", "CA"]
reference_units = [4, 9]

[[reference.repeat_loci.alleles]]
units = [4, 9]
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = [4, 10]
frequency = 0.5
"#;

/// **Con dos loci se colocan los dos.** Ninguna otra prueba lee el FASTA con
/// más de un locus, y un productor que colocara solo el primero pasaba la
/// batería entera.
///
/// Debe fallar contra: colocar solo el primer locus, en el productor y en el
/// juez a la vez. Se ejecutó, y muere **antes** de esta aserción: la guarda de
/// admisión, que ahora corre también con referencia simulada, rechaza el tramo
/// al azar de `chr1`. Con la guarda apagada además, se ejecutó también y muere
/// en la aserción del tramo de `chr1`, de modo que la comprobación del FASTA
/// discrimina por sí sola.
#[test]
fn two_loci_on_different_sequences_are_both_placed() {
    let destination = publish(&scenario(TWO_LOCI, 2));
    let fasta = published_fasta(&destination);
    assert_eq!(span(&fasta["chr2"], 701, 730), "CA".repeat(15));
    assert_eq!(span(&fasta["chr1"], 101, 130), COMPOSITE_TRACT);
}

/// **Cada corrupción del tracto publicado muere en su propia rama del juez**, y
/// el rechazo cita el esquema de la versión del documento, la 6, no el v4 que
/// la constructora escribe fijo (DT-017).
///
/// Cada caso cambia **un solo** campo de un documento que, reserializado sin
/// tocar, se acepta, y exige el detalle de su rama: un juez que rechazara por
/// otro motivo no pasa.
///
/// **Una rama no se alcanza desde aquí, y se dice:** «repeat_loci no es una
/// lista». Antes de llegar al juez, la comprobación de que la configuración
/// efectiva y el manifiesto registran los mismos loci ya rechaza un documento
/// cuyo `repeat_loci` no es una lista; el caso se ejecuta abajo y se exige ese
/// otro rechazo, para que nadie lea la rama como cubierta.
///
/// Debe fallar contra: un juez que admita `reference_units` en un locus de un
/// motivo, y uno que no exija que los bloques cubran el tramo entero. Las dos
/// mutaciones sobrevivían a la batería entera; se ejecutaron y mueren aquí.
#[test]
fn every_corruption_of_a_published_tract_is_rejected_by_its_own_branch() {
    let destination = publish(&scenario(TWO_LOCI, 2));
    let effective_path = destination.join("effective-config.json");
    let manifest_path = destination.join("manifest.json");
    let manifest = fs::read(&manifest_path).expect("manifiesto");
    let original: Value =
        serde_json::from_slice(&fs::read(&effective_path).expect("configuración efectiva"))
            .expect("JSON");
    assert_eq!(original["repeat_loci"][0]["id"], "STR-2");
    assert_eq!(original["repeat_loci"][1]["id"], "STR-C");
    let judge = |document: &Value| {
        validate_published_resources(
            &destination,
            &effective_path,
            &serde_json::to_vec_pretty(document).expect("JSON"),
            &manifest_path,
            &manifest,
            None,
        )
    };
    judge(&original).expect("el documento reserializado se acepta");

    type Mutation = fn(&mut Value);
    let cases: [(&str, Mutation, &str); 11] = [
        (
            "sin secuencia",
            |document| {
                document["repeat_loci"][1]
                    .as_object_mut()
                    .expect("locus")
                    .remove("sequence_id");
            },
            "locus STR-C en la referencia simulada: sin secuencia",
        ),
        (
            "inicio cero",
            |document| document["repeat_loci"][1]["start"] = serde_json::json!(0),
            "el inicio no es una coordenada de base uno",
        ),
        (
            "final antes del inicio",
            |document| document["repeat_loci"][1]["end"] = serde_json::json!(100),
            "el final no es posterior al inicio",
        ),
        (
            "motivos vacíos",
            |document| document["repeat_loci"][1]["motifs"] = serde_json::json!([]),
            "los motivos no son una lista de textos no vacíos",
        ),
        (
            "unidades que no son enteros",
            |document| {
                document["repeat_loci"][1]["reference_units"] = serde_json::json!(["4", "9"]);
            },
            "reference_units no es una lista de enteros",
        ),
        (
            "unidades en un locus simple",
            |document| document["repeat_loci"][0]["reference_units"] = serde_json::json!([15]),
            "locus STR-2 en la referencia simulada: un locus de un solo motivo publica reference_units",
        ),
        (
            "compuesto sin unidades",
            |document| {
                document["repeat_loci"][1]
                    .as_object_mut()
                    .expect("locus")
                    .remove("reference_units");
            },
            "no publica reference_units",
        ),
        (
            "un valor para dos bloques",
            |document| document["repeat_loci"][1]["reference_units"] = serde_json::json!([4]),
            "reference_units no tiene un valor por motivo",
        ),
        (
            "bloque que se sale del tramo",
            |document| {
                document["repeat_loci"][1]["reference_units"] = serde_json::json!([4, 20]);
            },
            "un bloque de reference_units está vacío o se sale del tramo",
        ),
        (
            "bloque vacío",
            |document| {
                document["repeat_loci"][1]["reference_units"] = serde_json::json!([0, 15]);
            },
            "un bloque de reference_units está vacío o se sale del tramo",
        ),
        (
            "bloques que no cubren el tramo",
            |document| {
                document["repeat_loci"][1]["reference_units"] = serde_json::json!([4, 5]);
            },
            "los bloques de reference_units no cubren el tramo entero",
        ),
    ];
    for (name, mutate, expected) in cases {
        let mut document = original.clone();
        mutate(&mut document);
        assert_ne!(document, original, "{name}: la corrupción no cambió nada");
        let error = judge(&document).expect_err(name);
        let detail = error.context().get("detail").cloned().unwrap_or_default();
        assert!(detail.contains(expected), "{name}: {detail}");
        assert_eq!(
            error.context().get("schema").map(String::as_str),
            Some("effective-config-v6.schema.json"),
            "{name}: el rechazo cita el esquema de su versión"
        );
    }

    let mut not_a_list = original.clone();
    not_a_list["repeat_loci"] = serde_json::json!({"STR-2": {}});
    let error = judge(&not_a_list).expect_err("repeat_loci que no es una lista");
    let detail = error.context().get("detail").cloned().unwrap_or_default();
    assert!(
        detail.contains(
            "el manifiesto registra loci de repetición que la configuración efectiva no declara"
        ),
        "{detail}"
    );
}

/// **`motifs` con un solo motivo se rechaza**, en el lector y en el esquema: un
/// locus simple se declara con `motif`. Antes el lector lo admitía como simple
/// y el esquema exigía `reference_units` en cuanto veía `motifs`, de modo que
/// los dos se contradecían en los dos sentidos.
///
/// Debe fallar contra: un lector que admita la lista de un elemento. Se
/// ejecutó.
#[test]
fn a_one_element_motifs_list_is_rejected_in_favour_of_motif() {
    let one = CA_LOCUS_HEADER.replace("motif = \"CA\"", "motifs = [\"CA\"]");
    assert_ne!(one, CA_LOCUS_HEADER, "el fixture cambió");
    let loci = format!(
        "{one}
[[reference.repeat_loci.alleles]]
units = 15
frequency = 1.0
"
    );
    let (fields, detail) = rejection(&scenario(&loci, 2));
    assert_eq!(fields, "reference.repeat_loci.motifs");
    assert!(detail.contains("locus STR-2"), "{detail}");
    assert!(detail.contains("se declara con motif"), "{detail}");
}

/// **La guarda de admisión juzga también el tracto que la referencia simulada
/// coloca**, antes de simular. `ACG` sobre ocho bases coloca `ACGACGAC`: seis
/// de ocho bases en copias completas, pureza 0,75, por debajo del mínimo de
/// 0,8. Al lado, su pareja válida: el mismo motivo sobre nueve bases, pureza 1.
///
/// El alelo declarado tiene tres unidades enteras, de modo que la guarda de los
/// alelos, que corre al leer la configuración, lo admite: lo único que
/// rechaza es el tracto colocado.
///
/// Debe fallar contra: saltarse la guarda con referencia simulada, que es lo
/// que hacía el programa hasta la corrección. Se ejecutó.
#[test]
fn the_admission_guard_judges_the_tract_placed_in_a_simulated_reference() {
    let locus = |end: u64| {
        format!(
            r#"
[[reference.repeat_loci]]
id = "STR-A"
sequence = "chr2"
start = 101
end = {end}
motif = "ACG"

[[reference.repeat_loci.alleles]]
units = 3
frequency = 1.0
"#
        )
    };
    let destination = publish(&scenario(&locus(109), 2));
    assert_eq!(
        span(&published_fasta(&destination)["chr2"], 101, 109),
        "ACG".repeat(3)
    );

    let root = temporary("impure");
    fs::create_dir_all(&root).expect("directorio");
    let config = root.join("general.toml");
    fs::write(&config, scenario(&locus(108), 2)).expect("configuración");
    let effective =
        load_and_validate_config(&config).expect("la configuración se admite al leerla");
    let error = execute_validated_config(&effective, &root.join("published"))
        .expect_err("el tracto colocado no llega a la pureza mínima");
    let detail = error.context().get("detail").cloned().unwrap_or_default();
    assert_eq!(error.code(), "GPS-CONF-008");
    assert!(detail.contains("locus STR-A"), "{detail}");
    assert!(detail.contains("pureza 0.75"), "{detail}");
}

/// **Un compuesto que pide el generador se rechaza por el generador**, no
/// pidiéndole `reference_units`: con el campo declarado el rechazo seguiría
/// llegando, así que pedirlo era mandar al usuario por un camino sin salida. Y
/// el esquema no exige el campo cuando el locus declara el generador.
///
/// Debe fallar contra: resolver `reference_units` antes que el generador, que
/// es el orden de antes. Se ejecutó.
#[test]
fn a_composite_locus_with_the_generator_is_rejected_by_the_generator() {
    let loci = r#"
[[reference.repeat_loci]]
id = "STR-G"
sequence = "chr2"
start = 101
end = 130
motifs = ["CAG", "CA"]

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 1.0
"#;
    let (fields, detail) = rejection(&scenario(loci, 2));
    assert_ne!(fields, FIELD, "{detail}");
    assert!(detail.contains("el generador solo admite uno"), "{detail}");
}
