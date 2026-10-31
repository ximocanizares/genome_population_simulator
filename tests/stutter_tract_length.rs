//! La dependencia del tartamudeo con la longitud del tracto, sobre los
//! ficheros publicados.
//!
//! # Qué se comprueba aquí y qué no
//!
//! Las pruebas de la función pura viven en `src/stutter.rs` y son las que fijan
//! la tabla, el ancla, los extremos y la sustitución. Aquí se comprueba lo que
//! aquellas no pueden: que la modulación **llega de verdad a los bytes**, y que
//! el documento que la ejecución publica la lleva consigo con su marca de
//! tosquedad y su tabla entera.

use std::fs;
use std::io::Read;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::event_manifest::validate_event_manifest_bytes;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::validate_effective_config_against_manifest_bytes;
use noodles_bgzf as bgzf;
use serde_json::Value;

#[path = "support/detection_section.rs"]
mod detection_section;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

/// Unidades de los dos tractos del fixture.
///
/// **Diecisiete es la longitud de referencia** de los dinucleótidos —34 pares
/// de base—, donde el factor de modulación vale exactamente uno. **Cinco es el
/// extremo inferior** del tramo medido, diez pares de base, donde vale 0,0038.
/// Los dos comparten motivo, `CA`, y comparten la misma tasa de tabla: lo
/// único que los separa es su longitud.
const LONG_UNITS: usize = 17;
const SHORT_UNITS: usize = 5;

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-tract-length-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

/// La referencia del fixture, con los dos tractos dentro de los primeros 150
/// pares de base de cualquier molécula.
///
/// # La aritmética que garantiza que R1 los ve enteros, siempre
///
/// La molécula es `copia[s, s+I)` con inserto `I ∈ [250, 300]` y comienzo
/// `s ∈ [0, 300−I]`, de modo que `s ≤ 50`. R1 es `copia[s, s+150)`. Un tracto
/// que ocupa `[a, a+T)` cabe entero en R1 para **todo** comienzo posible
/// exactamente cuando `a ≥ 50` y `a + T ≤ 150`, contando `T` con la mayor
/// expansión posible, una unidad más.
///
/// El corto empieza en 60 y acaba, como mucho, en 72. El largo empieza en 76 y
/// acaba, como mucho, en 112. Los dos cumplen las dos condiciones con holgura.
fn reference_bases() -> String {
    let short = "CA".repeat(SHORT_UNITS);
    let long = "CA".repeat(LONG_UNITS);
    let head = format!("{}{short}GGGG{long}", "G".repeat(60));
    let bases = format!("{head}{}", "G".repeat(300 - head.len()));
    assert_eq!(bases.len(), 300);
    bases
}

struct Published {
    destination: PathBuf,
}

impl Published {
    fn json(&self, relative: &str) -> Value {
        serde_json::from_slice(&fs::read(self.destination.join(relative)).expect("fichero"))
            .expect("JSON")
    }

    fn effective(&self) -> Value {
        self.json("effective-config.json")
    }

    fn first_reads(&self) -> Vec<Vec<u8>> {
        let path = self
            .destination
            .join("outputs/fastq/reads-main/sample-A_R1.fastq.gz");
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
}

/// Las tiradas de `CA` de una lectura, **en orden de aparición**, con su
/// número de unidades.
///
/// No sirve buscar el separador `GGGG`: el relleno de la referencia también es
/// `G`, de modo que la primera aparición de `GGGG` está en el prefijo y no
/// entre los dos tractos. Lo que sí distingue a los dos tractos sin ambigüedad
/// es su **orden**: el corto va primero y el largo después, y entre ellos no
/// hay ninguna otra tirada de `CA`.
fn ca_runs(read: &[u8]) -> Vec<u64> {
    let mut runs = Vec::new();
    let mut offset = 0_usize;
    while offset + 2 <= read.len() {
        if &read[offset..offset + 2] == b"CA" {
            let mut units = 0_u64;
            while offset + 2 <= read.len() && &read[offset..offset + 2] == b"CA" {
                units += 1;
                offset += 2;
            }
            runs.push(units);
        } else {
            offset += 1;
        }
    }
    runs
}

/// Cuenta, entre las lecturas, cuántas muestran en ese tracto un número de
/// unidades distinto del verdadero.
///
/// `long` elige cuál de las dos tiradas se mira: la segunda o la primera.
fn stuttered(reads: &[Vec<u8>], long: bool, truth: u64) -> (u64, u64) {
    let mut covered = 0_u64;
    let mut moved = 0_u64;
    for read in reads {
        let runs = ca_runs(read);
        if runs.len() != 2 {
            // La lectura no cubre los dos tractos enteros. El control positivo
            // del denominador exige que esto no pase nunca en este fixture.
            continue;
        }
        let observed = if long { runs[1] } else { runs[0] };
        covered += 1;
        if observed != truth {
            moved += 1;
        }
    }
    (moved, covered)
}

fn run(coverage: u32) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    fs::write(
        root.join("reference.fa"),
        format!(">chr1\n{}\n", reference_bases()),
    )
    .expect("FASTA");
    fs::write(
        root.join("population-snvs.tsv"),
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tG\tT\tP1\t0.5\n",
    )
    .expect("matriz");
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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
length = 300

[[reference.repeat_loci]]
id = "STR-SHORT"
sequence = "chr1"
start = 61
end = 70
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{short}"
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-LONG"
sequence = "chr1"
start = 75
end = 108
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{long}"
frequency = 1.0

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = {coverage}
sequencing_errors = "disabled"

[technical_event.stutter]
enabled = true

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
rate = 1.0

[[technical_event.stutter.steps]]
units = -3
weight = 0.05

[[technical_event.stutter.steps]]
units = -2
weight = 0.15

[[technical_event.stutter.steps]]
units = -1
weight = 0.6

[[technical_event.stutter.steps]]
units = 1
weight = 0.2

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
            short = "CA".repeat(SHORT_UNITS),
            long = "CA".repeat(LONG_UNITS),
        )),
    )
    .expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Published { destination }
}

/// **Un tracto largo tartamudea medible y fuertemente más que uno corto del
/// mismo motivo, sobre los ficheros publicados.**
///
/// Los dos loci comparten motivo y comparten la tasa de tabla, que vale uno.
/// Lo único que los separa es la longitud: diecisiete unidades frente a cinco.
/// Con la modulación, el largo está en su ancla y tartamudea **siempre**,
/// mientras que el corto se queda en 0,0038 y no tartamudea **casi nunca**.
///
/// **Debe fallar contra** el modelo anterior, en el que la tasa de tabla valía
/// uno para los dos y los dos tartamudeaban siempre; y contra una modulación
/// aplicada al revés, que dejaría al corto por encima del largo.
#[test]
fn a_long_tract_stutters_far_more_than_a_short_one_of_the_same_motif() {
    let published = run(400);
    let reads = published.first_reads();
    assert!(reads.len() >= 100, "el corpus tiene que ser representativo");

    let (long_moved, long_covered) = stuttered(&reads, true, LONG_UNITS as u64);
    let (short_moved, short_covered) = stuttered(&reads, false, SHORT_UNITS as u64);

    // **Control positivo del denominador.** Sin él, un fixture en el que las
    // lecturas no cubrieran un tracto daría cero tartamudeos y la desigualdad
    // de abajo se cumpliría por la razón equivocada.
    assert_eq!(
        long_covered,
        reads.len() as u64,
        "todas las lecturas tienen que cubrir el tracto largo entero"
    );
    assert_eq!(
        short_covered,
        reads.len() as u64,
        "todas las lecturas tienen que cubrir el tracto corto entero"
    );

    let long_fraction = long_moved as f64 / long_covered as f64;
    let short_fraction = short_moved as f64 / short_covered as f64;
    assert!(
        long_fraction > 0.95,
        "en su ancla y con la tasa a uno el tracto largo tartamudea siempre: {long_fraction}"
    );
    assert!(
        short_fraction < 0.05,
        "con el factor a 0,0038 el tracto corto casi nunca tartamudea: {short_fraction}"
    );
    // **No hay aquí una tercera aserción de cociente**, y conviene decir por
    // qué: `long > short × 20` estaría **implicada** por las dos anteriores
    // —0,95 frente a 0,05 ya son diecinueve veces— y no añadiría nada. Quien
    // hace el trabajo son los dos umbrales absolutos, que es donde se separan
    // el modelo modulado y el plano.
}

/// **El documento publicado lleva la marca de tosquedad, la declaración de que
/// no se extrapola, las longitudes de referencia y la curva entera.**
///
/// Es el tercero de los tres caminos por los que la marca viaja —constante de
/// Rust, guarda que rechaza el documento y valor obligatorio en el esquema—, y
/// el único que un lector del fichero puede comprobar sin el programa.
#[test]
fn the_published_document_carries_the_curve_with_its_coarseness_mark() {
    let published = run(100);
    let effective = published.effective();
    let dependence = &effective["technical_event"]["stutter"]["tract_length_dependence"];

    assert_eq!(dependence["coarse"], true);
    assert_eq!(dependence["extrapolated"], false);
    let reason = dependence["coarse_reason"].as_str().expect("motivo");
    for fragment in [
        "read by eye",
        "normalised at the top of the measured range",
        "nothing is extrapolated",
        "carry no tract-length dependence",
    ] {
        assert!(
            reason.contains(fragment),
            "el motivo de la tosquedad tiene que decir «{fragment}»: {reason}"
        );
    }
    assert!(
        dependence["source"]
            .as_str()
            .expect("fuente")
            .contains("Fungtammasan"),
        "la curva tiene que viajar con su fuente"
    );

    let rows = dependence["by_motif_length"]
        .as_array()
        .expect("tabla por longitud de motivo");
    assert_eq!(rows.len(), 6);
    // Mononucleótidos y dinucleótidos tienen curva y ancla; los demás no, y
    // eso se publica dicho, no callado.
    assert_eq!(rows[0]["measured"], true);
    assert_eq!(rows[0]["reference_tract_bases"], 15);
    assert_eq!(rows[1]["measured"], true);
    assert_eq!(rows[1]["reference_tract_bases"], 34);
    for row in &rows[2..] {
        assert_eq!(row["measured"], false);
        assert!(
            row["reference_tract_bases"].is_null(),
            "un motivo sin curva no puede publicar longitud de referencia"
        );
        assert_eq!(
            row["factors"].as_array().expect("factores").len(),
            0,
            "un motivo sin curva publica la tabla vacía"
        );
    }
    // La curva entera, para que pueda rehacerse sin este programa: el último
    // factor de cada una vale exactamente uno, porque es el ancla.
    for (index, expected_len) in [(0_usize, 11_usize), (1, 13)] {
        let factors = rows[index]["factors"].as_array().expect("factores");
        assert_eq!(factors.len(), expected_len);
        assert_eq!(factors[expected_len - 1]["factor"], 1.0);
        assert!(
            factors[0]["factor"].as_f64().expect("factor") < 0.01,
            "el extremo corto tiene que estar muy por debajo del ancla"
        );
    }
}

/// **LAS DOS guardas —la de la configuración efectiva y la del manifiesto de
/// evento— rechazan el documento al que se le quita la marca de tosquedad, se
/// le mueve el tramo medido o se le cambia un factor de la curva.**
///
/// # Por qué los dos documentos y no uno
///
/// Porque son **dos guardas distintas, en dos ficheros distintos**, y que una
/// rechace no dice nada de la otra. La del manifiesto de evento llegó a tener
/// ciento treinta líneas que **ninguna prueba ejercitaba como rechazo**:
/// sustituirla entera por una línea vacía dejaba pasar la batería completa.
///
/// Es el segundo de los tres caminos. Cada mutación se aplica **sobre los
/// bytes publicados** y se vuelve a validar, de modo que la prueba distingue
/// una guarda que comprueba de una que solo declara.
///
/// **Debe fallar contra** una guarda que acepte cualquiera de las cuatro
/// mutaciones. El control positivo —que el documento intacto se acepta—
/// impide que un rechazo indiscriminado la satisfaga.
#[test]
fn the_guard_rejects_a_document_that_loses_the_mark_or_moves_the_curve() {
    let published = run(100);
    let effective_path = published.destination.join("effective-config.json");
    let manifest_path = published.destination.join("manifest.json");
    let event_path = published
        .destination
        .join("events/reads-main/manifest.json");
    let effective_bytes = fs::read(&effective_path).expect("configuración efectiva");
    let manifest_bytes = fs::read(&manifest_path).expect("manifiesto");
    let event_bytes = fs::read(&event_path).expect("manifiesto de evento");

    // Controles positivos: intactos, los dos documentos se aceptan. Sin ellos,
    // una guarda que rechazara cualquier cosa satisfaría el resto.
    validate_effective_config_against_manifest_bytes(
        &effective_path,
        &effective_bytes,
        &manifest_path,
        &manifest_bytes,
    )
    .expect("la configuración efectiva intacta tiene que aceptarse");
    validate_event_manifest_bytes(&event_path, &event_bytes)
        .expect("el manifiesto de evento intacto tiene que aceptarse");

    type Locate = for<'a> fn(&'a mut Value) -> &'a mut Value;
    type Mutation = (&'static str, fn(&mut Value));
    // Las seis mutaciones se aplican a **los dos** documentos, cada uno por su
    // raíz, porque cada uno tiene su propia guarda y son dos códigos
    // distintos. Que la de la configuración efectiva rechace no dice
    // absolutamente nada de la del manifiesto de evento.
    let mutations: [Mutation; 6] = [
        ("la marca de tosquedad en falso", |dependence| {
            dependence["coarse"] = Value::Bool(false);
        }),
        ("la declaración de no extrapolar en cierto", |dependence| {
            dependence["extrapolated"] = Value::Bool(true);
        }),
        ("la procedencia de la curva cambiada", |dependence| {
            dependence["source"] = Value::from("otra figura de otro artículo");
        }),
        ("una longitud de referencia cambiada", |dependence| {
            dependence["by_motif_length"][0]["reference_tract_bases"] = Value::from(13);
        }),
        ("el tramo medido cambiado", |dependence| {
            dependence["by_motif_length"][0]["measured_units_max"] = Value::from(14);
        }),
        ("un factor de la curva cambiado", |dependence| {
            dependence["by_motif_length"][0]["factors"][3]["factor"] = Value::from(0.5);
        }),
    ];

    let documents: [(&str, &[u8], Locate); 2] = [
        (
            "configuración efectiva",
            &effective_bytes,
            effective_dependence_of,
        ),
        ("manifiesto de evento", &event_bytes, event_dependence_of),
    ];

    for (document_name, bytes, locate) in documents {
        for (what, mutate) in mutations {
            let mut document: Value = serde_json::from_slice(bytes).expect("JSON");
            let before = serde_json::to_vec(&document).expect("antes");
            mutate(locate(&mut document));
            let mutated = serde_json::to_vec(&document).expect("después");
            // **Que la mutación mutó**, comprobado y no supuesto.
            assert_ne!(
                before, mutated,
                "en {document_name}, la mutación «{what}» no cambió nada"
            );
            let outcome = if document_name == "configuración efectiva" {
                validate_effective_config_against_manifest_bytes(
                    &effective_path,
                    &mutated,
                    &manifest_path,
                    &manifest_bytes,
                )
            } else {
                validate_event_manifest_bytes(&event_path, &mutated)
            };
            assert!(
                outcome.is_err(),
                "la guarda de {document_name} tiene que rechazar «{what}»"
            );
        }
    }
}

fn effective_dependence_of(document: &mut Value) -> &mut Value {
    &mut document["technical_event"]["stutter"]["tract_length_dependence"]
}

fn event_dependence_of(document: &mut Value) -> &mut Value {
    &mut document["artifacts"]["stutter"]["tract_length_dependence"]
}

// ---------------------------------------------------------------------------
// La modulación del tartamudeo CLONAL, sobre los bytes publicados
// ---------------------------------------------------------------------------

/// Unidades de los dos tractos del escenario amplificado.
///
/// **Diez está dentro del tramo medido y lejos del ancla**, donde el factor
/// vale 0,124; diecisiete es el ancla, donde vale uno. Es justo la
/// combinación que faltaba: todos los escenarios clonales anteriores se habían
/// quedado **donde la modulación es invisible** —en el ancla, por encima del
/// tramo, o con la probabilidad por ciclo declarada a mano, que prevalece sin
/// modular—, de modo que quitar la modulación clonal entera no mataba ninguna
/// prueba.
const CLONAL_SHORT_UNITS: usize = 10;
const CLONAL_LONG_UNITS: usize = 17;

/// Moléculas del escenario amplificado. Ver la nota de tamaño muestral en la
/// docstring de la prueba.
const CLONAL_MOLECULES: u32 = 4_000;

/// Referencia del escenario amplificado: dos tractos del **mismo motivo** y
/// distinta longitud, los dos dentro de los primeros 150 pares de base.
///
/// # La aritmética, hecha con la expansión en las DOS condiciones
///
/// La molécula es `copia[s, s+I)` con `I ∈ [250, 300]`, de modo que `s ≤ 50`,
/// y R1 es `copia[s, s+150)`. Un tracto en `[a, a+T)` cabe entero en R1 para
/// **todo** comienzo posible cuando `a ≥ 50` y `a + T ≤ 150`, **contando `T`
/// con la mayor extensión que el tracto pueda alcanzar**.
///
/// Aquí los pesos del paso dejan **solo la contracción**, de modo que los dos
/// tractos únicamente pueden encoger y `T` es su longitud de referencia. El
/// corto ocupa `[60, 80)` y el largo `[84, 118)`: **diez bases de margen por
/// la izquierda y treinta y dos por la derecha**. Encoger solo acerca los dos
/// límites, nunca los cruza.
fn clonal_reference_bases() -> String {
    let short = "CA".repeat(CLONAL_SHORT_UNITS);
    let long = "CA".repeat(CLONAL_LONG_UNITS);
    let head = format!("{}{short}GGGG{long}", "G".repeat(60));
    let bases = format!("{head}{}", "G".repeat(300 - head.len()));
    assert_eq!(bases.len(), 300);
    bases
}

fn run_amplified() -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    fs::write(
        root.join("reference.fa"),
        format!(">chr1\n{}\n", clonal_reference_bases()),
    )
    .expect("FASTA");
    fs::write(
        root.join("population-snvs.tsv"),
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tG\tT\tP1\t0.5\n",
    )
    .expect("matriz");
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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
length = 300

[[reference.repeat_loci]]
id = "STR-SHORT"
sequence = "chr1"
start = 61
end = 80
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{short}"
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-LONG"
sequence = "chr1"
start = 85
end = 118
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{long}"
frequency = 1.0

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = {coverage}
sequencing_errors = "disabled"
library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0
cycles = 10
efficiency = 0.9

[technical_event.stutter]
enabled = true

[[technical_event.stutter.locus_rates]]
locus_id = "STR-SHORT"
rate = 0.0

[[technical_event.stutter.locus_rates]]
locus_id = "STR-LONG"
rate = 0.0

[[technical_event.stutter.steps]]
units = -3
weight = 0.0

[[technical_event.stutter.steps]]
units = -2
weight = 0.0

[[technical_event.stutter.steps]]
units = -1
weight = 1.0

[[technical_event.stutter.steps]]
units = 1
weight = 0.0

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
            short = "CA".repeat(CLONAL_SHORT_UNITS),
            long = "CA".repeat(CLONAL_LONG_UNITS),
            coverage = CLONAL_MOLECULES,
        )),
    )
    .expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Published { destination }
}

/// **La probabilidad por ciclo del tartamudeo CLONAL depende de la longitud
/// del tracto, y se ve en los FASTQ publicados.**
///
/// # Qué agujero cierra, y por qué hacía falta esta prueba y no otra
///
/// Toda la razón de ser de la versión nueva del tartamudeo clonal es que su
/// probabilidad por ciclo se module con la longitud del tracto. Antes de esta
/// prueba **ninguna lo comprobaba**: deshacer la modulación en el punto donde
/// el sorteo la consume dejaba pasar la batería entera. Los escenarios
/// clonales existentes estaban todos donde la modulación es invisible —en el
/// ancla, por encima del tramo medido, o con la probabilidad declarada a mano,
/// que prevalece sin modular—.
///
/// # Cómo se aísla la parte clonal de la parte de base
///
/// Los dos loci llevan su tasa de base **sustituida a cero**, y la sustitución
/// prevalece sin modular. Con eso el tartamudeo de base no aporta nada y todo
/// lo que las lecturas muestran viene del linaje, cuya probabilidad por ciclo
/// es la **derivada** de las dos columnas —no una declarada— y por tanto la
/// que la modulación toca. Los pesos dejan solo la contracción de una unidad,
/// de modo que la dirección no añade ruido y la geometría del fixture es
/// exacta.
///
/// # El tamaño muestral, calculado y no elegido a ojo
///
/// Con diez ciclos, la fracción acumulada predicha es **0,0273** en el ancla y
/// **0,00336** a diez unidades: un cociente de **8,1**. Con 4000 moléculas, el
/// umbral de 0,012 que separa las dos queda a **9,5 desviaciones** por encima
/// de lo que produce el tracto corto modulado y a **5,9 desviaciones** por
/// debajo de lo que produciría sin modular.
///
/// **Debe fallar contra** una entrega que ignore la longitud del tracto en la
/// probabilidad por ciclo del linaje: el tracto corto subiría de 0,0034 a
/// 0,0273, la misma cifra que el largo, y cruzaría el umbral.
#[test]
fn the_clonal_per_cycle_probability_depends_on_the_tract_length() {
    let published = run_amplified();
    let reads = published.first_reads();
    assert!(
        reads.len() as u32 >= CLONAL_MOLECULES,
        "el plan de cobertura tiene que dar al menos las moléculas calculadas: {}",
        reads.len()
    );

    let (short_moved, short_covered) = stuttered(&reads, false, CLONAL_SHORT_UNITS as u64);
    let (long_moved, long_covered) = stuttered(&reads, true, CLONAL_LONG_UNITS as u64);

    // **Control positivo del denominador**, porque sin él un tracto que las
    // lecturas no cubrieran daría cero y la desigualdad se cumpliría por la
    // razón equivocada.
    assert_eq!(
        short_covered,
        reads.len() as u64,
        "todas las lecturas tienen que cubrir el tracto corto entero"
    );
    assert_eq!(
        long_covered,
        reads.len() as u64,
        "todas las lecturas tienen que cubrir el tracto largo entero"
    );

    let short_fraction = short_moved as f64 / short_covered as f64;
    let long_fraction = long_moved as f64 / long_covered as f64;

    assert!(
        long_fraction > 0.015,
        "en su ancla el tracto largo tiene que acumular cerca de 0,0273: {long_fraction:.5}"
    );
    assert!(
        short_fraction < 0.012,
        "a diez unidades la modulación tiene que dejarlo cerca de 0,0034, \
         muy por debajo del 0,0273 que daría sin modular: {short_fraction:.5}"
    );
    // Y la separación dicha como magnitud, que es lo que distingue esta prueba
    // de una desigualdad que cualquier ruido satisfaría.
    assert!(
        long_fraction > short_fraction * 3.0,
        "la separación predicha es de 8,1 veces: {long_fraction:.5} frente a {short_fraction:.5}"
    );
    // Control positivo del numerador: el tracto corto **sí** tartamudea alguna
    // vez. Sin esto, una entrega que apagara el linaje entero pasaría las tres
    // aserciones anteriores.
    assert!(
        short_moved > 0,
        "el tracto corto tiene que tartamudear alguna vez, o el linaje está apagado"
    );
}
