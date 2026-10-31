//! F6.4.4, F6.4.11 y F6.4.12: el **generador de alelos por deslizamiento
//! biológico** dentro del recorrido completo.
//!
//! # Cómo están elegidos los fixtures, que es lo que hace que puedan fallar
//!
//! El defecto número uno de este proyecto es una prueba aparcada donde el
//! efecto que dice medir es invisible. Aquí cada fixture está elegido para que
//! la propiedad **se vea**, y se dice cuál:
//!
//! - **El abanico corriente** usa una media de **tres pasos** sobre sesenta
//!   linajes y un tope muy por encima del tracto. Con una media mucho más baja
//!   ningún linaje se movería y todo el abanico sería el alelo de la
//!   referencia; con una mucho más alta, o con un tope alcanzable, lo que se
//!   mediría sería el tope. Con tres pasos salen nueve alelos repartidos
//!   alrededor del de la referencia, que es lo que hay que ver.
//! - **El agotamiento del tracto** usa un tracto de **cinco unidades** y sesgo
//!   de expansión **cero** con media seis. Así unos linajes se detienen en el
//!   mínimo y otros no: si todos se detuvieran, el abanico sería un solo alelo
//!   y no se distinguiría de un modelo que saturase en silencio.
//! - **El escenario mixto** declara **dos loci**, uno con alelos escritos a
//!   mano y otro generado. Un escenario de un solo tipo **no discrimina** un
//!   documento que marcara todos los alelos con el mismo origen, y así se dice
//!   en su docstring.
//! - **El locus que ordena antes** se inserta **en la misma secuencia**, a la
//!   izquierda del que ya había. Uno en otro cromosoma **no serviría**: dos
//!   loci de cromosomas distintos no comparten contenedor y los dos serían el
//!   índice cero, de modo que un ámbito indexado por posición sobreviviría.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::{
    validate_general_manifest_bytes, validate_published_resources,
};
use serde_json::Value;

#[path = "support/detection_section.rs"]
mod detection_section;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f6-4-11-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Trescientas bases con **dos tractos**, uno a la izquierda del otro y en la
/// misma secuencia.
///
/// El de la izquierda, `(GT)25` en 51..100, es el que permite comprobar que
/// insertar un locus que **ordena antes** no mueve el abanico del que ya había.
/// El de la derecha, `(CA)50` en 101..200, es el locus de trabajo.
fn reference_bases() -> String {
    format!(
        "{}{}{}{}",
        "G".repeat(50),
        "GT".repeat(25),
        "CA".repeat(50),
        "G".repeat(100)
    )
}

struct Published {
    root: PathBuf,
    destination: PathBuf,
    manifest: Value,
    effective: Value,
}

/// Escenario con referencia suministrada y `individuals` individuos diploides.
fn publish(loci: &str, individuals: usize) -> Published {
    let ids = (1..=individuals)
        .map(|index| format!("\"I{index:03}\""))
        .collect::<Vec<_>>()
        .join(", ");
    publish_with_populations(
        loci,
        &format!(
            r#"[[initial_populations]]
population_id = "P1"
individual_count = {individuals}
individual_ids = [{ids}]
ploidy = 2
genotype_model = "hwe"
"#
        ),
        &["P1"],
    )
}

/// El mismo escenario con las poblaciones iniciales que se den, ya escritas, y
/// la lista de sus identificadores para el generador de la matriz.
fn publish_with_populations(loci: &str, populations: &str, population_ids: &[&str]) -> Published {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    let bases = reference_bases();
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let snv_populations = population_ids
        .iter()
        .map(|id| {
            format!(
                r#"
[[population_snvs.generate.populations]]
population_id = "{id}"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
"#
            )
        })
        .collect::<String>();
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
# **La variante se confina a las primeras cincuenta bases**, que están fuera de
# los dos loci. Una posición sorteada dentro de un tracto declarado rechaza el
# escenario por la decisión S-32, y con el sorteo libre eso dependía de dónde
# cayera: la prueba habría fallado por el fixture y no por el producto.
spatial_mode = "weighted_regions"
outside_weight = 0.0
substitution_profile = "uniform_snv12_v1"

[[population_snvs.generate.regions]]
contig = "chr1"
start = 1
end = 50
weight = 1.0
{snv_populations}
[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {}
{loci}
{populations}"#,
            bases.len()
        )),
    )
    .expect("configuración");
    let effective_config = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective_config, &destination).expect("publicación");
    let manifest = serde_json::from_slice(
        &fs::read(destination.join("manifest.json")).expect("manifiesto publicado"),
    )
    .expect("JSON del manifiesto");
    let effective = serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("configuración efectiva"),
    )
    .expect("JSON de la configuración efectiva");
    Published {
        root,
        destination,
        manifest,
        effective,
    }
}

fn assert_external_schema(path: &Path, schema: &str) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(path)
        .output()
        .expect("validador externo de JSON Schema");
    assert!(
        output.status.success(),
        "validador externo para {path:?} contra {schema}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// El locus de trabajo: `(CA)50` generado con media de tres pasos y sesgo
/// simétrico, con el tope muy por encima.
const GENERATED_RIGHT: &str = r#"
[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5
"#;

/// El mismo locus con **otro delante**, en la misma secuencia y con un inicio
/// menor, que es lo que caza un ámbito indexado por posición o por orden.
const GENERATED_LEFT_AND_RIGHT: &str = r#"
[[reference.repeat_loci]]
id = "STR-LEFT"
sequence = "chr1"
start = 51
end = 100
motif = "GT"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5

[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5
"#;

/// Los mismos dos loci **declarados al revés**, para separar el orden de la
/// posición: el contenido es idéntico y solo cambia en qué orden se escriben.
const GENERATED_RIGHT_AND_LEFT: &str = r#"
[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5

[[reference.repeat_loci]]
id = "STR-LEFT"
sequence = "chr1"
start = 51
end = 100
motif = "GT"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5
"#;

/// El locus de trabajo **con otro nombre**, en el mismo sitio y con los mismos
/// parámetros. Es lo único que separa un ámbito por identificador de uno por
/// posición.
const GENERATED_RIGHT_RENAMED: &str = r#"
[[reference.repeat_loci]]
id = "STR-RENAMED"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5
"#;

/// Un locus **declarado a mano** y otro **generado** en el mismo documento.
const DECLARED_AND_GENERATED: &str = r#"
[[reference.repeat_loci]]
id = "STR-LEFT"
sequence = "chr1"
start = 51
end = 100
motif = "GT"

[[reference.repeat_loci.alleles]]
units = 25
frequency = 0.6

[[reference.repeat_loci.alleles]]
units = 24
frequency = 0.4

[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5
"#;

fn locus<'a>(document: &'a Value, section: &str, id: &str) -> &'a Value {
    let loci = if section == "manifest" {
        document["repeat_loci"]["loci"]
            .as_array()
            .expect("loci del manifiesto")
    } else {
        document["repeat_loci"].as_array().expect("loci efectivos")
    };
    loci.iter()
        .find(|locus| locus["id"] == id)
        .unwrap_or_else(|| panic!("el documento no publica el locus {id}"))
}

// ---------------------------------------------------------------------------
// F6.4.11 y F6.4.12 — el abanico y su publicación
// ---------------------------------------------------------------------------

/// **Los dos documentos publicados superan su validación externa, marcan cada
/// alelo como generado y registran los parámetros, las versiones y la
/// contabilidad del generador.**
///
/// El fixture produce **varios alelos** repartidos alrededor del de la
/// referencia; con un solo alelo la comprobación de que las frecuencias suman
/// uno se cumpliría sola y la del abanico no mediría nada.
#[test]
fn the_generated_fan_is_published_with_its_origin_parameters_and_accounting() {
    let published = publish(GENERATED_RIGHT, 30);
    assert_external_schema(
        &published.destination.join("effective-config.json"),
        "effective-config-v6.schema.json",
    );
    assert_external_schema(
        &published.destination.join("manifest.json"),
        "manifest-general-v6.schema.json",
    );

    let effective = locus(&published.effective, "effective", "STR-RIGHT");
    let alleles = effective["alleles"].as_array().expect("alelos");
    assert!(
        alleles.len() > 1,
        "el fixture tiene que producir varios alelos para que el abanico se vea"
    );
    assert!(
        alleles.iter().all(|allele| allele["origin"] == "generated"),
        "un alelo generado se publicó con otro origen"
    );
    // **Ningún alelo con fracciones de unidad.** Se parte de un tracto entero y
    // cada paso mueve una unidad completa, así que una fracción solo puede
    // venir de un paso que no sea de una unidad.
    for allele in alleles {
        let units = allele["units"].as_f64().expect("unidades");
        assert_eq!(units.fract(), 0.0, "alelo con fracción de unidad: {units}");
    }
    let sum: f64 = alleles
        .iter()
        .map(|allele| allele["effective_frequency"].as_f64().expect("frecuencia"))
        .sum();
    assert!((sum - 1.0).abs() < 1e-9, "las frecuencias suman {sum}");

    let slippage = &effective["slippage"];
    assert_eq!(slippage["model"], "stepwise");
    assert_eq!(slippage["step_units"], 1);
    assert_eq!(slippage["mean_events"], 3.0);
    assert_eq!(slippage["expansion_bias"], 0.5);
    assert_eq!(slippage["algorithm"], "biological-slippage-v1");
    assert_eq!(slippage["count_algorithm"], "poisson-count-v1");
    assert_eq!(slippage["bias_algorithm"], "bernoulli-v1");
    assert_eq!(slippage["rng_scope"], "initial_repeat_slippage");
    assert_eq!(slippage["rng_scope_catalog_version"], "10");
    assert_eq!(
        slippage["limitations"].as_array().map(Vec::len),
        Some(3),
        "las tres limitaciones declaradas tienen que viajar con el documento"
    );
    // **Un linaje por copia haplotípica**: treinta individuos diploides.
    assert_eq!(slippage["lineages"], 60);
    let drawn = slippage["drawn_events"].as_u64().expect("sorteados");
    let applied = slippage["applied_events"].as_u64().expect("aplicados");
    let unapplied = slippage["unapplied_events"].as_u64().expect("no aplicados");
    assert!(drawn > 0, "con media tres algún linaje tiene que moverse");
    assert_eq!(applied + unapplied, drawn);

    // El manifiesto publica **lo mismo** del generador y, además, lo realizado.
    let manifest = locus(&published.manifest, "manifest", "STR-RIGHT");
    assert_eq!(manifest["slippage"], effective["slippage"]);
    assert_eq!(manifest["drawn_copies"], 60);
    assert!(
        manifest["alleles"]
            .as_array()
            .expect("alelos")
            .iter()
            .all(|allele| allele["origin"] == "generated")
    );
    fs::remove_dir_all(published.root).ok();
}

/// **El abanico entero, congelado: qué alelo recibe cuántos linajes.**
///
/// # Por qué esta prueba existe, dicho sin rodeos
///
/// Es la afirmación central del generador —**la frecuencia de cada alelo es la
/// proporción de linajes que terminan en él**— y no la sostenía ninguna prueba
/// capaz de fallar. Se midió: emparejando cada alelo con la frecuencia del
/// simétrico de la lista, el espectro se publicaba **del revés** y la batería
/// entera seguía en verde. El conjunto de alelos no cambia, las frecuencias
/// siguen sumando uno y la contabilidad sigue cuadrando; lo único que cambia es
/// **a qué alelo le toca qué**, que es justo lo que nadie miraba.
///
/// # Por qué el fixture es este y no otro
///
/// Cinco individuos diploides dan **diez linajes**, de modo que cada frecuencia
/// es un múltiplo de un décimo y la tabla se puede escribir entera. Y el
/// espectro que sale es **asimétrico**: `0,1 · 0,3 · 0,3 · 0,2 · 0,1`. Eso es
/// deliberado y es la mitad de la prueba. Con cuatro individuos salía
/// `0,125 · 0,25 · 0,25 · 0,25 · 0,125`, que es **simétrico**, y la mutación del
/// espectro invertido habría sobrevivido a una tabla congelada.
///
/// Los valores se obtuvieron del árbol del 27 de septiembre de 2026. **No se
/// actualizan sin más si cambian**: una tabla que se mueve significa o que el
/// modelo cambió —y entonces hay que decidirlo por escrito— o que hay un
/// defecto.
#[test]
fn the_whole_fan_matches_its_frozen_table_allele_by_allele() {
    let published = publish(GENERATED_RIGHT, 5);
    let effective = locus(&published.effective, "effective", "STR-RIGHT");

    // unidades → frecuencia declarada, que es la proporción de linajes.
    const FROZEN: [(u64, f64); 5] = [(45, 0.1), (47, 0.3), (49, 0.3), (50, 0.2), (51, 0.1)];
    // El espectro **no** es simétrico, y sin eso la prueba no discriminaría un
    // abanico emparejado al revés.
    let frequencies: Vec<f64> = FROZEN.iter().map(|(_, frequency)| *frequency).collect();
    let mut reversed = frequencies.clone();
    reversed.reverse();
    assert_ne!(
        frequencies, reversed,
        "el fixture tiene que ser asimétrico para que la tabla discrimine"
    );

    let alleles = effective["alleles"].as_array().expect("alelos");
    assert_eq!(alleles.len(), FROZEN.len(), "el abanico cambió de tamaño");
    for (index, (units, frequency)) in FROZEN.into_iter().enumerate() {
        let allele = &alleles[index];
        assert_eq!(
            allele["units"].as_f64(),
            Some(units as f64),
            "alelo {index}: las unidades se movieron"
        );
        assert_eq!(
            allele["declared_frequency"].as_f64(),
            Some(frequency),
            "alelo {index} ({units} unidades): la proporción de linajes se movió"
        );
        assert_eq!(allele["effective_frequency"].as_f64(), Some(frequency));
        assert_eq!(allele["sequence"], "CA".repeat(units as usize));
    }
    let slippage = &effective["slippage"];
    assert_eq!(slippage["lineages"], 10);
    assert_eq!(slippage["drawn_events"], 28);
    assert_eq!(slippage["applied_events"], 28);
    assert_eq!(slippage["unapplied_events"], 0);

    // El manifiesto publica el **mismo** abanico, y no una segunda opinión.
    let manifest = locus(&published.manifest, "manifest", "STR-RIGHT");
    for (index, (_, frequency)) in FROZEN.into_iter().enumerate() {
        assert_eq!(
            manifest["alleles"][index]["declared_frequency"].as_f64(),
            Some(frequency)
        );
    }
    fs::remove_dir_all(published.root).ok();
}

/// **Las tres limitaciones se publican con su texto entero, palabra por
/// palabra.**
///
/// # Por qué los textos van escritos a mano aquí
///
/// El juez del manifiesto compara lo publicado contra la constante que lo
/// publica: por sí solo es un espejo. Se midió sustituyendo la primera
/// limitación entera por `"ancestry: fine"`, y la batería siguió verde. Cambiar
/// lo que el simulador le cuenta al usuario sobre las limitaciones de su propio
/// modelo era **invisible**.
///
/// Estos tres literales son la tercera copia —las otras dos son la constante de
/// Rust y DTC-095, la decisión que fija el generador— y se comparan contra los
/// **bytes publicados**. Tocar la constante obliga a tocar esta prueba, y
/// entonces el cambio deja de ser silencioso.
#[test]
fn the_three_limitations_are_published_word_for_word() {
    const LIMITATIONS: [&str; 3] = [
        "no shared ancestry: every haplotype copy derives independently from the reference allele, so related individuals share no history at this locus and the resulting frequency spectrum is flatter than a real population's, where drift concentrates a few alleles; it matters if the dataset is used to test population-structure tools",
        "the step count is not a process in time: each copy accumulates a random number of mutations without stating over how many generations, which is consistent with this being a single initialisation draw and not a simulation of generations",
        "the reference genome may carry an allele that no individual carries: it does happen in reality, but it must be known before chaining an analysis that rules it out",
    ];
    let published = publish(GENERATED_RIGHT, 4);
    for (document, section) in [
        (&published.effective, "effective"),
        (&published.manifest, "manifest"),
    ] {
        let slippage = &locus(document, section, "STR-RIGHT")["slippage"];
        let limitations = slippage["limitations"].as_array().expect("limitaciones");
        assert_eq!(limitations.len(), 3, "{section}: tienen que ser tres");
        for (index, expected) in LIMITATIONS.into_iter().enumerate() {
            assert_eq!(
                limitations[index].as_str(),
                Some(expected),
                "{section}: la limitación {index} no se publica entera"
            );
        }
    }
    // La primera lleva la frase que dice **por qué importa**, que es la que
    // faltaba en el texto publicado antes de esta prueba.
    assert!(LIMITATIONS[0].contains("population-structure tools"));
    fs::remove_dir_all(published.root).ok();
}

/// **El alelo de la referencia entra en el abanico solo si algún linaje termina
/// en él, y con media cero el abanico es exactamente ese alelo.**
///
/// Es el criterio del responsable: su presencia **no se fuerza**. Con media
/// cero ningún linaje se mueve y el abanico tiene que ser el alelo de la
/// referencia y nada más; que además **no se consuma aleatoriedad** lo
/// comprueba la prueba unitaria del módulo, que sí puede contar palabras.
#[test]
fn a_zero_mean_publishes_only_the_reference_allele() {
    let published = publish(
        r#"
[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 0.0
expansion_bias = 0.5
"#,
        4,
    );
    let effective = locus(&published.effective, "effective", "STR-RIGHT");
    let alleles = effective["alleles"].as_array().expect("alelos");
    assert_eq!(alleles.len(), 1);
    assert_eq!(alleles[0]["sequence"], "CA".repeat(50));
    assert_eq!(alleles[0]["origin"], "generated");
    assert_eq!(alleles[0]["effective_frequency"], 1.0);
    assert_eq!(effective["slippage"]["drawn_events"], 0);
    assert_eq!(effective["slippage"]["applied_events"], 0);
    assert_eq!(effective["slippage"]["unapplied_events"], 0);
    fs::remove_dir_all(published.root).ok();
}

/// **Un tracto que llega a su mínimo detiene al linaje, y los eventos que no se
/// aplicaron se publican aparte.**
///
/// El fixture es un tracto de **cinco unidades** con sesgo de expansión cero y
/// media seis: unos linajes agotan el tracto y otros no. Si todos lo agotaran,
/// el abanico sería un solo alelo y un modelo que saturase en silencio daría el
/// mismo tracto final; la cifra de eventos no aplicados es lo único que separa
/// a los dos.
#[test]
fn exhausting_the_tract_stops_the_lineage_and_the_lost_events_are_counted_apart() {
    let published = publish(
        r#"
[[reference.repeat_loci]]
id = "STR-SHORT"
sequence = "chr1"
start = 51
end = 60
motif = "GT"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 6.0
expansion_bias = 0.0

# El locus son las cinco primeras unidades del `(GT)×25` de la referencia, que
# el detector da entero (bases 52 a 100): la región lo contiene, y mandan las
# fronteras declaradas (REQ-REF-40). Suelo de una copia solo para los motivos
# de dos bases, que son los de la prueba.
[reference.repeat_detection]
min_copies = [1000, 1, 1000, 1000, 1000, 1000]

[reference.repeat_detection.slippage]
model = "stepwise"
mean_events = 1.0

[[reference.repeat_detection.regions]]
sequence = "chr1"
start = 51
end = 100
"#,
        30,
    );
    let effective = locus(&published.effective, "effective", "STR-SHORT");
    let slippage = &effective["slippage"];
    let drawn = slippage["drawn_events"].as_u64().expect("sorteados");
    let applied = slippage["applied_events"].as_u64().expect("aplicados");
    let unapplied = slippage["unapplied_events"].as_u64().expect("no aplicados");
    assert_eq!(applied + unapplied, drawn);
    assert!(
        unapplied > 0,
        "el fixture tiene que agotar el tracto en algún linaje"
    );
    assert!(
        applied > 0,
        "el fixture tiene que aplicar pasos en algún linaje: sin ellos no se \
distinguiría un generador que no hiciera nada"
    );
    let alleles = effective["alleles"].as_array().expect("alelos");
    assert!(
        alleles.len() > 1,
        "con sesgo cero y media seis tienen que quedar linajes a medio contraer"
    );
    // **Ningún alelo por debajo de una unidad completa**, que es el suelo de la
    // regla de agotamiento.
    for allele in alleles {
        assert!(
            allele["longest_run_units"].as_u64().expect("tirada") >= 1,
            "un alelo se quedó sin ninguna copia completa"
        );
        assert!(allele["units"].as_f64().expect("unidades") >= 1.0);
    }
    // Y **todos contraen**: con sesgo cero ninguno puede pasar de las cinco
    // unidades del tracto de la referencia.
    assert!(
        alleles
            .iter()
            .all(|allele| allele["units"].as_f64().expect("unidades") <= 5.0),
        "con sesgo cero ningún linaje puede expandirse"
    );
    fs::remove_dir_all(published.root).ok();
}

/// **Un documento con loci de los dos tipos marca cada alelo con su origen.**
///
/// *Lo que un escenario de un solo tipo no discrimina, dicho sin adornos:* un
/// productor que escribiera siempre el mismo origen pasaría cualquier prueba
/// hecha sobre un escenario con loci de una sola clase. Hace falta el documento
/// **mixto**, que es el de aquí.
#[test]
fn a_mixed_scenario_marks_declared_and_generated_alleles_in_the_same_document() {
    let published = publish(DECLARED_AND_GENERATED, 20);
    for (document, section) in [
        (&published.effective, "effective"),
        (&published.manifest, "manifest"),
    ] {
        let declared = locus(document, section, "STR-LEFT");
        let generated = locus(document, section, "STR-RIGHT");
        assert!(
            declared["alleles"]
                .as_array()
                .expect("alelos declarados")
                .iter()
                .all(|allele| allele["origin"] == "declared"),
            "{section}: un alelo escrito a mano se marcó como generado"
        );
        assert!(
            generated["alleles"]
                .as_array()
                .expect("alelos generados")
                .iter()
                .all(|allele| allele["origin"] == "generated"),
            "{section}: un alelo generado se marcó como declarado"
        );
        assert!(
            declared.get("slippage").is_none(),
            "{section}: un locus declarado no puede llevar bloque de generador"
        );
        assert!(generated.get("slippage").is_some());
    }
    assert_external_schema(
        &published.destination.join("manifest.json"),
        "manifest-general-v6.schema.json",
    );
    fs::remove_dir_all(published.root).ok();
}

/// **El abanico depende del identificador del locus y de nada más**: ni de su
/// posición en la referencia, ni del orden en que se declaró.
///
/// # Las tres comprobaciones, y por qué hacen falta las tres
///
/// Es la prueba del ámbito, y una sola de sus partes **no basta**. Se midió
/// mutando:
///
/// 1. **Insertar un locus que ordena antes** —en la misma secuencia y a la
///    izquierda— no mueve el abanico del que ya había. Uno en otro cromosoma
///    no serviría: dos loci de cromosomas distintos no comparten contenedor y
///    los dos serían el índice cero. Esta parte caza un ámbito indexado por el
///    **orden**.
/// 2. **Declarar los mismos dos loci al revés** tampoco los mueve, que es la
///    misma propiedad sin necesidad de añadir nada.
/// 3. **Renombrar el locus, dejándolo en el mismo sitio y con los mismos
///    parámetros, sí cambia su abanico.** Esta es la que caza un ámbito
///    indexado por la **posición**: con el inicio del tracto como clave, las
///    dos primeras comprobaciones pasan tan campantes —el inicio no se mueve
///    al insertar nada a la izquierda— y solo el cambio de nombre lo delata.
///    Se comprobó: una versión anterior de esta prueba, sin esta tercera
///    parte, **sobrevivía** a sustituir el identificador por el inicio.
#[test]
fn the_fan_follows_the_locus_identifier_and_neither_its_position_nor_its_order() {
    let alone = publish(GENERATED_RIGHT, 30);
    let with_neighbour = publish(GENERATED_LEFT_AND_RIGHT, 30);
    let reversed = publish(GENERATED_RIGHT_AND_LEFT, 30);
    let renamed = publish(GENERATED_RIGHT_RENAMED, 30);

    let before = locus(&alone.effective, "effective", "STR-RIGHT");
    let after = locus(&with_neighbour.effective, "effective", "STR-RIGHT");
    assert_eq!(
        before["alleles"], after["alleles"],
        "añadir un locus a la izquierda movió el abanico del de la derecha"
    );
    assert_eq!(before["slippage"], after["slippage"]);
    assert_eq!(
        locus(&reversed.effective, "effective", "STR-RIGHT")["alleles"],
        before["alleles"],
        "declarar los loci en otro orden movió el abanico"
    );
    assert_eq!(
        locus(&reversed.effective, "effective", "STR-LEFT")["alleles"],
        locus(&with_neighbour.effective, "effective", "STR-LEFT")["alleles"],
        "declarar los loci en otro orden movió el abanico del de la izquierda"
    );

    // **El mismo tracto, el mismo sitio, los mismos parámetros y otro nombre.**
    let other_name = locus(&renamed.effective, "effective", "STR-RENAMED");
    assert_eq!(
        other_name["start"], before["start"],
        "control: el locus renombrado tiene que estar en el mismo sitio"
    );
    assert_eq!(
        other_name["slippage"]["mean_events"],
        before["slippage"]["mean_events"]
    );
    assert_ne!(
        other_name["alleles"], before["alleles"],
        "renombrar el locus no cambió su abanico: el ámbito no depende de su \
identificador, sino de su posición o de su orden"
    );

    // **Control positivo**: el locus vecino existe de verdad y tiene su propio
    // abanico, distinto del otro. Sin él, una comparación que mirase un locus
    // ausente pasaría igual.
    let neighbour = locus(&with_neighbour.effective, "effective", "STR-LEFT");
    assert_ne!(neighbour["alleles"], after["alleles"]);
    fs::remove_dir_all(alone.root).ok();
    fs::remove_dir_all(with_neighbour.root).ok();
    fs::remove_dir_all(reversed.root).ok();
    fs::remove_dir_all(renamed.root).ok();
}

/// **La misma semilla publica el mismo abanico, byte a byte.**
///
/// *Lo que esta prueba mide y lo que no, dicho sin adornos.* Mide
/// **determinismo**: que el generador no dependa de nada fuera de sus entradas.
/// **No dice nada sobre corrección**, porque compara la función consigo misma:
/// un abanico equivocado sale igual de equivocado las dos veces. Quien sostiene
/// la corrección es la tabla congelada de
/// `the_whole_fan_matches_its_frozen_table_allele_by_allele`, el fixture
/// calculado a mano del módulo y el juez independiente del documento publicado.
#[test]
fn the_same_seed_publishes_the_same_fan() {
    let first = publish(GENERATED_RIGHT, 30);
    let second = publish(GENERATED_RIGHT, 30);
    assert_eq!(
        locus(&first.effective, "effective", "STR-RIGHT"),
        locus(&second.effective, "effective", "STR-RIGHT")
    );
    fs::remove_dir_all(first.root).ok();
    fs::remove_dir_all(second.root).ok();
}

// ---------------------------------------------------------------------------
// F6.4.4 — rechazos del servicio escalonado, antes de simular
// ---------------------------------------------------------------------------

fn rejection(loci: &str) -> String {
    let root = temporary("reject");
    fs::create_dir_all(&root).expect("directorio");
    let bases = reference_bases();
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {}
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I001", "I002"]
ploidy = 2
genotype_model = "hwe"
"#,
            bases.len()
        )),
    )
    .expect("configuración");
    let rendered = load_and_validate_config(&config)
        .err()
        .map(|diagnostic| diagnostic.render_human())
        .unwrap_or_else(|| "la configuración se aceptó".to_owned());
    fs::remove_dir_all(root).ok();
    rendered
}

/// **Un modelo distinto del escalonado y un salto de dos unidades se rechazan
/// antes de simular, nombrando el parámetro y el valor.**
///
/// El **control positivo** es la última comprobación: el mismo locus con el
/// modelo admitido se acepta. Sin él, un lector que rechazara siempre pasaría
/// la prueba entera.
#[test]
fn another_model_and_a_multi_unit_jump_are_rejected_naming_parameter_and_value() {
    let other_model = rejection(
        r#"
[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[reference.repeat_loci.slippage]
model = "infinite_alleles"
mean_events = 1.0
"#,
    );
    assert!(
        other_model.contains("slippage.model") && other_model.contains("infinite_alleles"),
        "el rechazo no nombra el parámetro y el valor: {other_model}"
    );

    let two_units = rejection(
        r#"
[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 1.0
step_units = 2
"#,
    );
    assert!(
        two_units.contains("slippage.step_units") && two_units.contains('2'),
        "el rechazo no nombra el parámetro y el valor: {two_units}"
    );

    let accepted = rejection(
        r#"
[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 1.0
step_units = 1
"#,
    );
    assert_eq!(
        accepted, "la configuración se aceptó",
        "el control positivo no se aceptó: {accepted}"
    );
}

/// **Declarar alelos y pedir además que se generen es excluyente**, y el
/// rechazo lo dice.
#[test]
fn declaring_alleles_and_asking_for_generation_at_once_is_rejected() {
    let rendered = rejection(
        r#"
[[reference.repeat_loci]]
id = "STR-RIGHT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 50
frequency = 1.0

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 1.0
"#,
    );
    assert!(
        rendered.contains("excluyentes"),
        "el rechazo no explica la exclusión: {rendered}"
    );
}

/// **Un microsatélite compuesto no puede generarse**, porque sin fronteras de
/// bloque declaradas no hay forma de partir el tracto de la referencia.
#[test]
fn a_composite_locus_cannot_be_generated() {
    let rendered = rejection(
        r#"
[[reference.repeat_loci]]
id = "STR-COMPOSITE"
sequence = "chr1"
start = 51
end = 200
motifs = ["GT", "CA"]

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 1.0
"#,
    );
    assert!(
        rendered.contains("2 motivos") && rendered.contains("solo admite uno"),
        "el rechazo no explica por qué: {rendered}"
    );
}

// ---------------------------------------------------------------------------
// F6.4.12 — las guardas del manifiesto, ejercitadas como rechazo
// ---------------------------------------------------------------------------

/// **Cada guarda nueva del manifiesto se ejercita mutando los bytes publicados
/// y exigiendo el rechazo exacto**, con control positivo del documento intacto.
///
/// Es la respuesta directa al hallazgo de la entrega anterior: ciento treinta
/// líneas de validación del manifiesto sobrevivieron a que se borraran enteras
/// porque ninguna prueba las provocaba. Aquí la primera corrupción es
/// justamente **borrar el bloque del generador entero**.
#[test]
fn the_rust_reader_rejects_every_corruption_of_the_generated_section() {
    let published = publish(GENERATED_RIGHT, 20);
    let path = published.destination.join("manifest.json");
    let bytes = fs::read(&path).expect("manifiesto publicado");
    validate_general_manifest_bytes(&path, &bytes).expect("el manifiesto intacto se acepta");

    let original: Value = serde_json::from_slice(&bytes).expect("JSON");
    // El control se hace contra el documento **reserializado**: si no, una
    // «corrupción» que no cambiara nada pasaría por el simple hecho de que el
    // formato difiere.
    let baseline = serde_json::to_vec(&original).expect("JSON reserializado");
    validate_general_manifest_bytes(&path, &baseline)
        .expect("el manifiesto reserializado se acepta");
    let corrupt = |mutate: &dyn Fn(&mut Value)| {
        let mut document = original.clone();
        mutate(&mut document);
        serde_json::to_vec(&document).expect("JSON corrompido")
    };
    let cases: Vec<(&str, &str, Vec<u8>)> = vec![
        (
            "bloque del generador borrado entero",
            "el origen del alelo contradice al del locus",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]
                    .as_object_mut()
                    .expect("locus")
                    .remove("slippage");
            }),
        ),
        (
            "origen inventado",
            "el origen del alelo no es declarado ni generado",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["alleles"][0]["origin"] =
                    serde_json::json!("inventado");
            }),
        ),
        (
            "alelo generado que se dice declarado",
            "el origen del alelo contradice al del locus",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["alleles"][0]["origin"] =
                    serde_json::json!("declared");
            }),
        ),
        (
            "modelo inventado",
            "el modelo del generador no es el escalonado",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["model"] =
                    serde_json::json!("infinite_alleles");
            }),
        ),
        (
            "salto de dos unidades",
            "el paso del generador no es de una unidad completa",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["step_units"] = serde_json::json!(2);
            }),
        ),
        (
            "media de pasos negativa",
            "la media de pasos no es un decimal no negativo",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["mean_events"] =
                    serde_json::json!(-1.0);
            }),
        ),
        (
            "sesgo fuera de su intervalo",
            "el sesgo de expansión no pertenece a [0,1]",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["expansion_bias"] =
                    serde_json::json!(1.5);
            }),
        ),
        (
            "versión de algoritmo cambiada",
            "una versión o un ámbito del generador no es el normativo",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["algorithm"] =
                    serde_json::json!("biological-slippage-v2");
            }),
        ),
        (
            "ámbito de azar cambiado por el del tartamudeo",
            "una versión o un ámbito del generador no es el normativo",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["rng_scope"] =
                    serde_json::json!("stutter");
            }),
        ),
        (
            "versión del catálogo de ámbitos cambiada",
            "una versión o un ámbito del generador no es el normativo",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["rng_scope_catalog_version"] =
                    serde_json::json!("9");
            }),
        ),
        (
            "una limitación recortada",
            "las limitaciones del generador no son las tres declaradas",
            corrupt(&|document| {
                let limitations = document["repeat_loci"]["loci"][0]["slippage"]["limitations"]
                    .as_array_mut()
                    .expect("limitaciones");
                limitations.pop();
            }),
        ),
        (
            "una limitación reescrita",
            "las limitaciones del generador no son las tres declaradas",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["limitations"][0] =
                    serde_json::json!("no limitations worth mentioning");
            }),
        ),
        (
            "contabilidad que no cuadra",
            "los pasos aplicados y los no aplicados no suman los sorteados",
            corrupt(&|document| {
                let slippage = &mut document["repeat_loci"]["loci"][0]["slippage"];
                let applied = slippage["applied_events"].as_u64().expect("aplicados");
                slippage["applied_events"] = serde_json::json!(applied + 1);
            }),
        ),
        (
            "sin un solo linaje",
            "los linajes del generador no son un entero positivo",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["lineages"] = serde_json::json!(0);
            }),
        ),
        (
            "campo desconocido en el bloque del generador",
            "el bloque del generador no tiene exactamente sus campos",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["slippage"]["sobra"] = serde_json::json!(1);
            }),
        ),
        (
            "más alelos que linajes simulados",
            "el abanico publica más alelos que linajes simulados",
            corrupt(&|document| {
                // Cada linaje termina en **un** alelo, así que no puede haber
                // más alelos que linajes. Se baja la cifra de linajes por
                // debajo del número de alelos y se recompone la contabilidad
                // para que lo único fuera de sitio sea esa relación.
                let slippage = &mut document["repeat_loci"]["loci"][0]["slippage"];
                slippage["lineages"] = serde_json::json!(1);
                slippage["drawn_events"] = serde_json::json!(0);
                slippage["applied_events"] = serde_json::json!(0);
                slippage["unapplied_events"] = serde_json::json!(0);
            }),
        ),
    ];
    for (name, expected, corrupted) in cases {
        assert_ne!(corrupted, baseline, "la corrupción «{name}» no corrompió");
        let error = validate_general_manifest_bytes(&path, &corrupted)
            .err()
            .unwrap_or_else(|| panic!("la corrupción «{name}» se aceptó"));
        let rendered = error.render_human();
        assert!(
            rendered.contains(expected),
            "la corrupción «{name}» se rechazó por otro motivo: {rendered}"
        );
    }
    fs::remove_dir_all(published.root).ok();
}

/// **Los alelos generados elevan la cota de bases de tracto, y el tope
/// declarado los rechaza.**
///
/// Es la rama que solo existe **después** de generar: al resolver la
/// configuración el locus no tiene ni un alelo, de modo que su cota es la del
/// tramo de la referencia y el tope no muerde. Cuando el abanico nace, la cota
/// sube y el tope tiene que volver a aplicarse; si no, el documento publicaría
/// una cota falsa y el juez de la contabilidad encontraría lo realizado por
/// encima de lo acotado.
///
/// **Control positivo**: el mismo escenario con un tope holgado se publica sin
/// problema, de modo que lo que rechaza es el tope y no el generador.
#[test]
fn the_generated_alleles_raise_the_tract_budget_and_a_tight_cap_rejects_them() {
    let rejected = publish_with_cap(30, Some(1_000))
        .err()
        .map(|diagnostic| diagnostic.render_human())
        .unwrap_or_else(|| "el escenario se aceptó".to_owned());
    assert!(
        rejected.contains("max_tract_output_bases") && rejected.contains("1000"),
        "el rechazo no nombra el tope ni su valor: {rejected}"
    );
    assert!(
        rejected.contains("deslizamiento"),
        "el rechazo no dice que la cota subió al generar: {rejected}"
    );
    publish_with_cap(30, Some(1_000_000)).expect("con un tope holgado se publica");
    publish_with_cap(30, None).expect("sin tope declarado no hay rechazo");
}

/// Publica el escenario del locus generado con un tope declarado de bases de
/// tracto, y devuelve el diagnóstico en vez de romper.
fn publish_with_cap(
    individuals: usize,
    cap: Option<u64>,
) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
    let root = temporary("cap");
    fs::create_dir_all(&root).expect("directorio");
    let bases = reference_bases();
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let ids = (1..=individuals)
        .map(|index| format!("\"I{index:03}\""))
        .collect::<Vec<_>>()
        .join(", ");
    let declared_cap = cap.map_or(String::new(), |cap| {
        format!("max_tract_output_bases = {cap}\n")
    });
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
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"
{declared_cap}
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {length}
{GENERATED_RIGHT}
[[initial_populations]]
population_id = "P1"
individual_count = {individuals}
individual_ids = [{ids}]
ploidy = 2
genotype_model = "hwe"
"#,
            length = bases.len()
        )),
    )
    .expect("configuración");
    let effective = load_and_validate_config(&config)?;
    let destination = root.join("published");
    let outcome = execute_validated_config(&effective, &destination).map(|_| ());
    fs::remove_dir_all(root).ok();
    outcome
}

// ---------------------------------------------------------------------------
// Prueba de sistema, conservada **con su aviso**
// ---------------------------------------------------------------------------

/// **Cambiar los ciclos de PCR no mueve ni un alelo de la verdad.**
///
/// # Esta prueba NO demuestra que los dos ámbitos de azar estén separados
///
/// Hay que decirlo sin rodeos, porque parece que sí y no lo es. En este
/// proyecto compartir ámbito significa **derivar la misma subsemilla**, y
/// DTC-081 exige que cada decisión arranque su propio generador desde su clave
/// y que no haya consumo acumulado, de modo que un consumidor **no desplaza** al
/// otro. Además la verdad se genera **antes** que las lecturas, así que ningún
/// parámetro de estas puede alcanzarla. La aserción se cumple sola.
///
/// La prueba del ámbito es comparar **las dos subsemillas**, y vive en
/// `tests/subseed_derivation.rs`, en
/// `biological_slippage_and_stutter_never_share_a_subseed`.
///
/// # Para qué sirve entonces, que no es poco
///
/// Es buena contra un generador que **leyera la configuración del tartamudeo**
/// por error: uno que tomara de ahí su sesgo, su tasa o su número de pasos
/// movería el abanico al cambiar los ciclos y moriría aquí. Por eso se
/// conserva.
#[test]
fn changing_the_pcr_cycles_does_not_move_a_single_generated_allele() {
    let ten = publish_amplified(10);
    let fourteen = publish_amplified(14);
    assert_ne!(
        ten.effective["technical_event"]["artifacts"]["pcr"]["cycles"],
        fourteen.effective["technical_event"]["artifacts"]["pcr"]["cycles"],
        "control positivo: los dos escenarios tienen que diferir en los ciclos"
    );
    assert_eq!(
        locus(&ten.effective, "effective", "STR-GEN"),
        locus(&fourteen.effective, "effective", "STR-GEN"),
        "los ciclos de PCR movieron el abanico de la verdad"
    );
    fs::remove_dir_all(ten.root).ok();
    fs::remove_dir_all(fourteen.root).ok();
}

/// Escenario con preparación amplificada por PCR y un locus **generado**.
///
/// La referencia lleva mil bases de flanco a cada lado del tracto porque el
/// perfil sortea insertos de unos trescientos cincuenta pares de base: sin
/// espacio alrededor, el muestreador no tendría de dónde sacarlos.
fn publish_amplified(cycles: u32) -> Published {
    let root = temporary("pcr");
    fs::create_dir_all(&root).expect("directorio");
    let bases = format!("{}{}{}", flank(1_000), "CA".repeat(50), flank(1_000));
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let config = root.join("general.toml");
    fs::write(
        &config,
        detection_section::with_detection_section(&format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"
ploidy = 2

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

[[reference.repeat_loci]]
id = "STR-GEN"
sequence = "chr1"
start = 1001
end = 1100
motif = "CA"
max_units = 90

[reference.repeat_loci.slippage]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5

[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I001", "I002"]
ploidy = 2
genotype_model = "hwe"

[technical_event]
event_id = "reads-main"
individual_ids = ["I001"]
target_coverage_wgs = 2.0
sequencing_errors = "disabled"
library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.2
cycles = {cycles}
efficiency = 0.9

[technical_event.stutter]
enabled = true

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
rate = 0.05

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
individual_id = "I001"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
            length = bases.len()
        )),
    )
    .expect("configuración");
    let effective_config = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective_config, &destination).expect("publicación");
    let manifest = serde_json::from_slice(
        &fs::read(destination.join("manifest.json")).expect("manifiesto publicado"),
    )
    .expect("JSON del manifiesto");
    let effective = serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("configuración efectiva"),
    )
    .expect("JSON de la configuración efectiva");
    Published {
        root,
        destination,
        manifest,
        effective,
    }
}

/// **El juez del abanico rechaza, y no cuelga, un número de linajes o una media
/// de pasos que el documento publicado no puede sostener** (DT-024).
///
/// El juez rehace el abanico con un paso por linaje y por evento, y los dos
/// números salían del documento sin acotar: con `lineages` de `u64::MAX` no
/// terminaba. Ahora los linajes tienen que ser las copias haplotípicas de las
/// poblaciones iniciales publicadas —lo que el productor usa— y la media, caer
/// en el tope del generador.
///
/// **Por qué el caso de dos millones va antes que el de `u64::MAX`:** la prueba
/// no puede colgar la batería si la cota desaparece. Con dos millones de linajes
/// el recorrido sin cota termina en segundos y rechaza **por otro motivo** —el
/// abanico rehecho no es el publicado—, de modo que la aserción sobre el motivo
/// muere ahí, antes de llegar al caso que no terminaría.
///
/// **Control positivo:** el documento sin tocar se acepta, de modo que lo que
/// rechaza es la corrupción y no el montaje.
#[test]
fn the_fan_judge_rejects_unsupportable_lineages_and_mean_events_without_hanging() {
    let published = publish(GENERATED_RIGHT, 30);
    let effective_path = published.destination.join("effective-config.json");
    let manifest_path = published.destination.join("manifest.json");
    let manifest = fs::read(&manifest_path).expect("manifiesto");
    let judge = |document: &Value| {
        validate_published_resources(
            &published.destination,
            &effective_path,
            &serde_json::to_vec_pretty(document).expect("JSON"),
            &manifest_path,
            &manifest,
            None,
        )
    };
    judge(&published.effective).expect("el documento sin tocar se acepta");
    assert_eq!(
        published.effective["repeat_loci"][0]["slippage"]["lineages"], 60,
        "el fixture son treinta diploides"
    );

    type Mutation = fn(&mut Value);
    let cases: [(&str, Mutation, &str); 3] = [
        (
            "dos millones de linajes",
            |document| {
                document["repeat_loci"][0]["slippage"]["lineages"] = serde_json::json!(2_000_000);
            },
            "los linajes no son las copias haplotípicas de las poblaciones iniciales",
        ),
        (
            "u64::MAX linajes",
            |document| {
                document["repeat_loci"][0]["slippage"]["lineages"] = serde_json::json!(u64::MAX);
            },
            "los linajes no son las copias haplotípicas de las poblaciones iniciales",
        ),
        (
            "una media de pasos sin tope",
            |document| {
                document["repeat_loci"][0]["slippage"]["mean_events"] = serde_json::json!(10_000.0);
            },
            "la media de pasos rebasa el tope que el generador admite",
        ),
    ];
    for (name, mutate, expected) in cases {
        let mut document = published.effective.clone();
        mutate(&mut document);
        assert_ne!(document, published.effective, "{name}: no cambió nada");
        let started = std::time::Instant::now();
        let error = judge(&document).expect_err(name);
        let detail = error.context().get("detail").cloned().unwrap_or_default();
        assert!(detail.contains(expected), "{name}: {detail}");
        assert!(
            started.elapsed() < std::time::Duration::from_secs(5),
            "{name}: el juez tardó {:?} en rechazar",
            started.elapsed()
        );
    }
    fs::remove_dir_all(published.root).ok();
}

/// **Controles positivos de las dos cotas del juez del abanico** (DT-024): lo
/// válido que una cota mal escrita rechazaría.
///
/// - **Dos poblaciones de ploidías distintas**: 3 individuos diploides y 2
///   tetraploides son 3×2 + 2×4 = 14 linajes. Una regla que sumara solo la
///   primera población (6), que fijara la ploidía en 2 (10) o que multiplicara
///   el total de individuos por la ploidía de la primera (10) rechazaría este
///   paquete, y ninguna otra prueba de este fichero combina deslizamiento con
///   varias poblaciones o con una ploidía distinta de 2.
/// - **Una media de pasos de exactamente 1000**, el tope: una cota exclusiva
///   (`0.0..1000.0`) la rechazaría.
///
/// Publicar ya pasa el juez, pero se vuelve a llamar explícitamente para que la
/// prueba no dependa de ello.
#[test]
fn the_fan_judge_accepts_mixed_ploidy_populations_and_the_maximum_mean() {
    let populations = r#"[[initial_populations]]
population_id = "P1"
individual_count = 3
individual_ids = ["A1", "A2", "A3"]
ploidy = 2
genotype_model = "hwe"

[[initial_populations]]
population_id = "P2"
individual_count = 2
individual_ids = ["B1", "B2"]
ploidy = 4
genotype_model = "hwe"
"#;
    let published = publish_with_populations(GENERATED_RIGHT, populations, &["P1", "P2"]);
    assert_eq!(
        published.effective["repeat_loci"][0]["slippage"]["lineages"], 14,
        "3×2 + 2×4"
    );
    let judge = |published: &Published| {
        validate_published_resources(
            &published.destination,
            &published.destination.join("effective-config.json"),
            &fs::read(published.destination.join("effective-config.json")).expect("efectiva"),
            &published.destination.join("manifest.json"),
            &fs::read(published.destination.join("manifest.json")).expect("manifiesto"),
            None,
        )
    };
    judge(&published).expect("dos poblaciones de ploidías distintas se aceptan");
    fs::remove_dir_all(published.root).ok();

    let maximum = GENERATED_RIGHT.replace("mean_events = 3.0", "mean_events = 1000.0");
    assert_ne!(maximum, GENERATED_RIGHT, "el fixture cambió");
    let published = publish(&maximum, 30);
    assert_eq!(
        published.effective["repeat_loci"][0]["slippage"]["mean_events"],
        1000.0
    );
    judge(&published).expect("la media máxima se acepta");
    fs::remove_dir_all(published.root).ok();
}

/// Bases de composición equilibrada y **sin periodo corto**, para que dos
/// posiciones distintas del genoma no produzcan la misma lectura.
fn flank(length: usize) -> String {
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
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
