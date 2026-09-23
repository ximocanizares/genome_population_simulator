//! F6.4.8: lo que la configuración efectiva y el manifiesto registran de las
//! repeticiones en tándem.
//!
//! # Cómo están elegidos los fixtures, que es lo que hace que puedan fallar
//!
//! - **La frecuencia realizada** se comprueba con un corpus de **cuatro
//!   copias** y frecuencias declaradas de 0,7 y 0,3. Con cuatro copias toda
//!   proporción realizada es un múltiplo de un cuarto, y ni 0,7 ni 0,3 lo son:
//!   la diferencia entre lo declarado y lo realizado está garantizada por
//!   aritmética y no depende del sorteo. Un manifiesto que copiara la
//!   frecuencia declarada bajo el nombre de la realizada muere aquí.
//! - **El corpus grande no discrimina** esa confusión, y así se dice en su
//!   docstring: lo que comprueba es la convergencia, que es otra cosa.
//! - **El tope y los mínimos** se comprueban con dos loci en el mismo
//!   escenario, uno que los sustituye y otro que no. Sin el segundo, un
//!   manifiesto que publicara siempre los valores de omisión podría pasar si
//!   los declarados coincidieran con ellos por casualidad.
//! - **La pureza** se compara con una cifra calculada a mano sobre la
//!   secuencia del alelo, escrita como la división que es.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::{
    validate_effective_config_against_manifest_bytes, validate_general_manifest_bytes,
};
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f6-4-8-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Trescientas bases: cien de flanco, cien de tracto `(CA)50` y cien más.
fn reference_bases() -> String {
    format!("{}{}{}", "G".repeat(100), "CA".repeat(50), "G".repeat(100))
}

/// Escenario con referencia suministrada y `individuals` individuos diploides.
fn scenario(loci: &str, individuals: usize) -> (PathBuf, PathBuf) {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    let bases = reference_bases();
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    // Cero individuos significa **ningún bloque de población inicial**, que es
    // un escenario admisible y el que destapó que la sección del manifiesto
    // dependiera del sorteo en vez de la configuración.
    let populations = if individuals == 0 {
        String::new()
    } else {
        let ids = (1..=individuals)
            .map(|index| format!("\"I{index:03}\""))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = {individuals}\nindividual_ids = [{ids}]\nploidy = 2\ngenotype_model = \"hwe\"\n"
        )
    };
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
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
length = {}
{loci}
{populations}"#,
            bases.len()
        ),
    )
    .expect("configuración");
    (root, config)
}

struct Published {
    root: PathBuf,
    manifest: Value,
    effective: Value,
    destination: PathBuf,
}

fn publish(loci: &str, individuals: usize) -> Published {
    let (root, config) = scenario(loci, individuals);
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
        manifest,
        effective,
        destination,
    }
}

/// Valida un documento publicado con el validador externo de JSON Schema, que
/// es un proceso ajeno al productor de Rust.
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

/// Un locus con dos alelos cuyas frecuencias declaradas **no son múltiplos de
/// un cuarto**, que es lo que garantiza la diferencia con lo realizado en un
/// corpus de cuatro copias.
const UNEVEN_FREQUENCIES: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 45
frequency = 0.7

[[reference.repeat_loci.alleles]]
units = 40
frequency = 0.3
"#;

/// Un alelo de cien bases con una interrupción `TT` tras la primera copia:
/// cuarenta y nueve copias completas, noventa y ocho bases cubiertas.
fn interrupted_locus() -> String {
    let sequence = format!("CATT{}", "CA".repeat(48));
    assert_eq!(sequence.len(), 100);
    format!(
        r#"
[[reference.repeat_loci]]
id = "STR-INT"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{sequence}"
frequency = 1.0
"#
    )
}

/// Dos loci en el mismo escenario: el primero **sustituye** el tope y los dos
/// mínimos; el segundo los deja en sus valores de omisión.
const OVERRIDDEN_AND_DEFAULT: &str = r#"
[[reference.repeat_loci]]
id = "STR-OVER"
sequence = "chr1"
start = 101
end = 150
motif = "CA"
max_units = 77
min_purity = 0.9
min_longest_run_units = 5

[[reference.repeat_loci.alleles]]
units = 25
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-PLAIN"
sequence = "chr1"
start = 151
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 25
frequency = 1.0
"#;

// ---------------------------------------------------------------------------
// Pruebas
// ---------------------------------------------------------------------------

/// **Los dos documentos publicados superan su validación externa y registran
/// los loci, la procedencia de la referencia, el motivo y los alelos con su
/// secuencia.**
#[test]
fn both_published_documents_validate_externally_and_carry_the_loci() {
    let published = publish(UNEVEN_FREQUENCIES, 2);
    assert_eq!(
        published.manifest["schema_version"], "6",
        "declarar loci promueve el manifiesto"
    );
    assert_eq!(published.effective["schema_version"], "6");
    assert_external_schema(
        &published.destination.join("manifest.json"),
        "manifest-general-v6.schema.json",
    );
    assert_external_schema(
        &published.destination.join("effective-config.json"),
        "effective-config-v6.schema.json",
    );

    let section = &published.manifest["repeat_loci"];
    assert_eq!(
        section["reference_source"], "supplied",
        "la procedencia de la referencia decide si la pureza es observada o construida"
    );
    assert_eq!(section["reference_id"], "ref-main");
    let locus = &section["loci"][0];
    assert_eq!(locus["id"], "STR-1");
    assert_eq!(locus["sequence_id"], "chr1");
    assert_eq!(locus["start"], 101);
    assert_eq!(locus["end"], 200);
    assert_eq!(locus["motifs"], serde_json::json!(["CA"]));
    assert_eq!(locus["motif_lengths"], serde_json::json!([2]));
    let sequences = locus["alleles"]
        .as_array()
        .expect("alelos")
        .iter()
        .map(|allele| allele["sequence"].as_str().expect("secuencia").to_owned())
        .collect::<Vec<_>>();
    assert_eq!(sequences, ["CA".repeat(45), "CA".repeat(40)]);
    fs::remove_dir_all(&published.root).ok();
}

/// **La frecuencia realizada difiere de la declarada en un corpus pequeño.**
///
/// *Qué hace que pueda fallar:* cuatro copias y frecuencias declaradas de 0,7 y
/// 0,3. Toda proporción realizada sobre cuatro copias es un múltiplo de un
/// cuarto, y ninguna de las dos declaradas lo es, de modo que la diferencia no
/// depende del sorteo. Un manifiesto que publicara la declarada bajo el nombre
/// de la realizada —la variante defectuosa que el paso nombra— no puede pasar.
#[test]
fn the_realized_frequency_differs_from_the_declared_one_in_a_small_corpus() {
    let published = publish(UNEVEN_FREQUENCIES, 2);
    let locus = &published.manifest["repeat_loci"]["loci"][0];
    assert_eq!(locus["drawn_copies"], 4, "dos individuos diploides");
    let alleles = locus["alleles"].as_array().expect("alelos");
    let mut realized_total = 0.0;
    let mut copies_total = 0;
    for allele in alleles {
        let declared = allele["declared_frequency"].as_f64().expect("declarada");
        let realized = allele["realized_frequency"].as_f64().expect("realizada");
        let copies = allele["realized_copies"].as_u64().expect("copias");
        assert!(
            (realized - declared).abs() > 1e-9,
            "la realizada ({realized}) no puede ser la declarada ({declared})"
        );
        assert!(
            (realized - copies as f64 / 4.0).abs() < 1e-12,
            "la realizada es el recuento entre las copias sorteadas"
        );
        realized_total += realized;
        copies_total += copies;
    }
    assert_eq!(copies_total, 4, "las copias sorteadas se reparten enteras");
    assert!((realized_total - 1.0).abs() < 1e-12);
    fs::remove_dir_all(&published.root).ok();
}

/// **Y converge a la declarada en un corpus grande.**
///
/// *Esta prueba no discrimina la confusión anterior*, y conviene decirlo sin
/// adornos: un manifiesto que copiara la frecuencia declarada bajo el nombre de
/// la realizada la pasaría con holgura, porque la diferencia sería exactamente
/// cero. Lo que comprueba es otra cosa —que el sorteo sigue las frecuencias que
/// se le dan—, y quien caza aquella confusión es la prueba del corpus pequeño.
#[test]
fn the_realized_frequency_converges_to_the_declared_one_in_a_large_corpus() {
    let published = publish(UNEVEN_FREQUENCIES, 150);
    let locus = &published.manifest["repeat_loci"]["loci"][0];
    assert_eq!(locus["drawn_copies"], 300);
    for allele in locus["alleles"].as_array().expect("alelos") {
        let declared = allele["declared_frequency"].as_f64().expect("declarada");
        let realized = allele["realized_frequency"].as_f64().expect("realizada");
        assert!(
            (realized - declared).abs() < 0.05,
            "con trescientas copias la realizada ({realized}) se acerca a la declarada ({declared})"
        );
    }
    fs::remove_dir_all(&published.root).ok();
}

/// **La pureza publicada por alelo coincide con la calculada a mano sobre su
/// secuencia, y las interrupciones son las suyas.**
///
/// *La cuenta, hecha a mano:* la secuencia es `CA` `TT` y luego `CA` cuarenta y
/// ocho veces. El barrido voraz de izquierda a derecha anota la primera copia,
/// tropieza con las dos `T` y anota las cuarenta y ocho restantes: cuarenta y
/// nueve copias completas, noventa y ocho bases cubiertas de cien. La pureza es
/// esa división y la interrupción son las dos bases del desplazamiento 2.
#[test]
fn the_published_purity_matches_the_one_computed_by_hand_over_the_sequence() {
    let published = publish(&interrupted_locus(), 2);
    let allele = &published.manifest["repeat_loci"]["loci"][0]["alleles"][0];
    assert_eq!(allele["sequence"], format!("CATT{}", "CA".repeat(48)));
    assert!(
        (allele["purity"].as_f64().expect("pureza") - 98.0 / 100.0).abs() < 1e-12,
        "pureza publicada: {}",
        allele["purity"]
    );
    assert_eq!(
        allele["interruptions"],
        serde_json::json!([{"offset": 2, "length": 2}]),
        "la interrupción es la pareja de bases que ninguna copia cubre"
    );
    assert_eq!(allele["longest_run_units"], 48);
    // La configuración efectiva dice lo mismo del mismo alelo: los dos
    // documentos no pueden contradecirse sobre un valor que no se declara.
    let effective_allele = &published.effective["repeat_loci"][0]["alleles"][0];
    assert_eq!(effective_allele["purity"], allele["purity"]);
    assert_eq!(effective_allele["interruptions"], allele["interruptions"]);
    fs::remove_dir_all(&published.root).ok();
}

/// **El tope y los mínimos publicados son los configurados, no los de omisión,
/// cuando la configuración los sustituye.**
///
/// *Qué hace que pueda fallar:* el escenario lleva **dos** loci, uno que los
/// sustituye y otro que no. El segundo fija cuáles son los valores de omisión,
/// de modo que la comparación del primero no puede salir bien por casualidad, y
/// un manifiesto que publicara siempre los de omisión moriría en el primero.
#[test]
fn the_published_cap_and_minimums_are_the_configured_ones_not_the_defaults() {
    let published = publish(OVERRIDDEN_AND_DEFAULT, 2);
    let section = &published.manifest["repeat_loci"];
    assert_eq!(
        section["admission_defaults"],
        serde_json::json!({"min_purity": 0.8, "min_longest_run_units": 3}),
        "los valores de omisión del escenario"
    );
    let loci = section["loci"].as_array().expect("loci");
    assert_eq!(loci.len(), 2);

    let overridden = loci
        .iter()
        .find(|locus| locus["id"] == "STR-OVER")
        .expect("locus con valores sustituidos");
    assert_eq!(overridden["max_units"], 77);
    assert_eq!(overridden["admission"]["min_purity"], 0.9);
    assert_eq!(overridden["admission"]["min_longest_run_units"], 5);

    let plain = loci
        .iter()
        .find(|locus| locus["id"] == "STR-PLAIN")
        .expect("locus sin sustituciones");
    assert_eq!(plain["max_units"], 200, "el tope de omisión");
    assert_eq!(plain["admission"]["min_purity"], 0.8);
    assert_eq!(plain["admission"]["min_longest_run_units"], 3);
    fs::remove_dir_all(&published.root).ok();
}

/// **Un escenario con loci y sin poblaciones iniciales registra los loci
/// igual, con cero copias sorteadas y sin frecuencia realizada.**
///
/// *Qué hace que pueda fallar:* la sección se construía desde el sorteo, de
/// modo que un escenario sin individuos materializados la dejaba vacía y
/// publicaba un manifiesto de la versión anterior frente a una configuración
/// efectiva que sí traía los loci. Aquí se exige el registro, la versión y el
/// emparejamiento de los dos documentos. Y de paso se recorre la única rama
/// que publica frecuencia realizada nula, que antes era inalcanzable.
#[test]
fn a_scenario_with_loci_and_no_initial_populations_still_registers_them() {
    let published = publish(UNEVEN_FREQUENCIES, 0);
    assert_eq!(
        published.manifest["schema_version"], "6",
        "declarar loci promueve el manifiesto aunque no haya poblaciones"
    );
    assert_eq!(published.effective["schema_version"], "6");
    assert_external_schema(
        &published.destination.join("manifest.json"),
        "manifest-general-v6.schema.json",
    );
    let locus = &published.manifest["repeat_loci"]["loci"][0];
    assert_eq!(locus["id"], "STR-1");
    assert_eq!(locus["max_units"], 200);
    assert_eq!(
        locus["drawn_copies"], 0,
        "sin poblaciones no se sortea ninguna copia"
    );
    for allele in locus["alleles"].as_array().expect("alelos") {
        assert_eq!(allele["realized_copies"], 0);
        assert!(
            allele["realized_frequency"].is_null(),
            "sin sorteo no hay proporción que publicar: {}",
            allele["realized_frequency"]
        );
        // Lo declarado sigue estando, que es lo que el requisito pide registrar.
        assert!(allele["declared_frequency"].as_f64().expect("declarada") > 0.0);
        assert!(!allele["sequence"].as_str().expect("secuencia").is_empty());
    }
    fs::remove_dir_all(&published.root).ok();
}

/// **El lector de Rust rechaza lo mismo que el esquema publicado, y algo más
/// que el esquema no puede expresar.**
///
/// *Qué hace que pueda fallar:* se toma el manifiesto **real** y se corrompe
/// de quince formas, cada una dirigida a una guarda distinta, exigiendo de
/// cada una que el lector la rechace **por su motivo exacto**: rechazar por
/// una guarda ajena no protege a la guarda que la corrupción dice proteger. La
/// primera aserción fija que el documento intacto se acepta: sin ella, un
/// lector que rechazara todo pasaría las demás. Seis de las corrupciones —el
/// par frecuencia/recuento, las listas paralelas, la suma de copias, la suma
/// de bases contra la secuencia, el tope de copias por alelo y la proporción
/// publicada sin sorteo— **el esquema no puede cazarlas**, porque dividir dos
/// campos, sumar una lista o comparar la longitud de propiedades hermanas no
/// es expresable en JSON Schema.
///
/// *Lo que esta prueba no cubre, dicho sin adornos.* El lector tiene muchas más
/// guardas que corrupciones hay aquí, y la mayoría de las restantes son
/// comprobaciones de tipo —«esto no es un objeto», «esto no es texto», «esta
/// cifra no es un entero»— que el esquema publicado ya rechaza antes y que
/// solo se alcanzarían fabricando un documento que el esquema no admite.
/// Forzarlas una a una añadiría ruido sin añadir discriminación, de modo que se
/// dejan deliberadamente sin corromper.
#[test]
fn the_rust_reader_rejects_a_repeat_section_as_strictly_as_its_schema() {
    let published = publish(UNEVEN_FREQUENCIES, 2);
    let path = published.destination.join("manifest.json");
    let bytes = fs::read(&path).expect("manifiesto publicado");
    validate_general_manifest_bytes(&path, &bytes).expect("el manifiesto intacto se acepta");

    // La comparación de control se hace contra el documento **reserializado**,
    // no contra los bytes publicados: si no, una «corrupción» que no cambiara
    // nada pasaría igualmente porque el formato difiere.
    let original: Value = serde_json::from_slice(&bytes).expect("JSON");
    let baseline = serde_json::to_vec(&original).expect("JSON reserializado");
    let corrupt = |mutate: &dyn Fn(&mut Value)| {
        let mut document = original.clone();
        mutate(&mut document);
        serde_json::to_vec(&document).expect("JSON corrompido")
    };
    let cases: Vec<(&str, &str, Vec<u8>)> = vec![
        (
            "sección vacía",
            "la sección repeat_loci no tiene exactamente sus campos",
            corrupt(&|document| document["repeat_loci"] = serde_json::json!({})),
        ),
        (
            "sección sin ningún locus",
            "la sección repeat_loci no publica ningún locus",
            corrupt(&|document| document["repeat_loci"]["loci"] = serde_json::json!([])),
        ),
        (
            "campo desconocido en la sección",
            "la sección repeat_loci no tiene exactamente sus campos",
            corrupt(&|document| document["repeat_loci"]["sobra"] = serde_json::json!(1)),
        ),
        (
            "procedencia de referencia inventada",
            "la procedencia de la referencia no es simulada ni suministrada",
            corrupt(&|document| {
                document["repeat_loci"]["reference_source"] = serde_json::json!("derivada");
            }),
        ),
        (
            "frecuencia realizada que contradice al recuento",
            "la frecuencia realizada no es las copias realizadas entre las sorteadas",
            corrupt(&|document| {
                let allele = &mut document["repeat_loci"]["loci"][0]["alleles"][0];
                // Un valor admisible y **distinto del correcto**, sea cual sea:
                // con cuatro copias la proporción correcta es múltiplo de un
                // cuarto, y estos dos la contradicen siempre.
                let correct = allele["realized_frequency"].as_f64().expect("realizada");
                allele["realized_frequency"] =
                    serde_json::json!(if correct <= 0.5 { 0.75 } else { 0.25 });
            }),
        ),
        (
            "copias realizadas que no suman las sorteadas",
            "las copias realizadas no suman las copias sorteadas del locus",
            corrupt(&|document| {
                // Se pone a cero el recuento del primer alelo **y su
                // frecuencia**, de modo que el par sigue cuadrando y lo único
                // que falla es la suma.
                let allele = &mut document["repeat_loci"]["loci"][0]["alleles"][0];
                allele["realized_copies"] = serde_json::json!(0);
                allele["realized_frequency"] = serde_json::json!(0.0);
            }),
        ),
        (
            "listas de bloque de distinta longitud",
            "las listas de bloque no miden lo que el recuento declara",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["alleles"][0]["block_bases"] =
                    serde_json::json!([80, 10]);
            }),
        ),
        (
            "longitud de motivo que no es la de su secuencia",
            "la longitud de un motivo de bloque no es la de su secuencia",
            corrupt(&|document| {
                document["repeat_loci"]["loci"][0]["alleles"][0]["block_motif_lengths"] =
                    serde_json::json!([3]);
            }),
        ),
        (
            "cota por componente por debajo de lo que ese componente realizó",
            "las bases de tracto realizadas superan su cota",
            corrupt(&|document| {
                // Se rebaja la cota de la verdad por debajo de lo que la verdad
                // realizó, y se recompone el total para que la suma siga
                // cuadrando: así lo único que puede fallar es la comparación.
                let accounting = &mut document["repeat_loci"]["tract_output"];
                let realized = accounting["realized_truth_tract_bases"]
                    .as_u64()
                    .expect("verdad realizada");
                let vcf_bound = accounting["vcf_tract_bases_bound"]
                    .as_u64()
                    .expect("cota vcf");
                accounting["truth_tract_bases_bound"] = serde_json::json!(realized - 1);
                accounting["total_tract_bases_bound"] = serde_json::json!(realized - 1 + vcf_bound);
            }),
        ),
        (
            "mínimos de admisión fuera de rango",
            "la pureza mínima no es finita ni pertenece a [0,1]",
            corrupt(&|document| {
                document["repeat_loci"]["admission_defaults"]["min_purity"] =
                    serde_json::json!(1.5);
            }),
        ),
        (
            "dos loci con el mismo identificador",
            "dos loci publicados comparten identificador",
            corrupt(&|document| {
                let first = document["repeat_loci"]["loci"][0].clone();
                document["repeat_loci"]["loci"]
                    .as_array_mut()
                    .expect("loci")
                    .push(first);
            }),
        ),
        (
            "un alelo con más copias que las sorteadas",
            "un alelo recibió más copias que las sorteadas del locus",
            corrupt(&|document| {
                let locus = &mut document["repeat_loci"]["loci"][0];
                let drawn = locus["drawn_copies"].as_u64().expect("sorteadas");
                locus["alleles"][0]["realized_copies"] = serde_json::json!(drawn + 1);
            }),
        ),
        (
            "frecuencia realizada sin copias sorteadas",
            "sin copias sorteadas la frecuencia realizada tiene que ser nula",
            corrupt(&|document| {
                // Se vacía el sorteo entero —ninguna copia en ningún alelo—, de
                // modo que la suma sigue cuadrando y lo único fuera de sitio es
                // una proporción publicada donde no hay nada que proporcionar.
                let locus = &mut document["repeat_loci"]["loci"][0];
                locus["drawn_copies"] = serde_json::json!(0);
                for allele in locus["alleles"].as_array_mut().expect("alelos") {
                    allele["realized_copies"] = serde_json::json!(0);
                    allele["realized_frequency"] = Value::Null;
                }
                locus["alleles"][0]["realized_frequency"] = serde_json::json!(0.5);
            }),
        ),
        (
            "bases de bloque que no suman la secuencia",
            "las bases de los bloques no suman las de la secuencia del alelo",
            corrupt(&|document| {
                // Se alarga el primer bloque en una unidad entera **y se
                // corrige su recuento de copias**, de modo que la coherencia
                // interna del bloque se mantiene y lo único que deja de cuadrar
                // es el total contra la secuencia.
                let allele = &mut document["repeat_loci"]["loci"][0]["alleles"][0];
                let motif = allele["block_motif_lengths"][0].as_u64().expect("longitud");
                let bases = allele["block_bases"][0].as_u64().expect("bases");
                let units = allele["block_unit_counts"][0].as_f64().expect("copias");
                allele["block_bases"][0] = serde_json::json!(bases + motif);
                allele["block_unit_counts"][0] = serde_json::json!(units + 1.0);
            }),
        ),
        (
            "motivo de bloque que no son bases",
            "un motivo publicado no son bases A, C, G o T",
            corrupt(&|document| {
                // Mismo número de caracteres, para que la comparación con la
                // longitud declarada siga cuadrando: lo único inadmisible son
                // las letras.
                let allele = &mut document["repeat_loci"]["loci"][0]["alleles"][0];
                let length = allele["block_motifs"][0].as_str().expect("motivo").len();
                allele["block_motifs"][0] = serde_json::json!("N".repeat(length));
            }),
        ),
        (
            "versión seis sin su sección",
            "la versión 6 exige la sección repeat_loci",
            corrupt(&|document| {
                document
                    .as_object_mut()
                    .expect("manifiesto")
                    .remove("repeat_loci");
            }),
        ),
    ];
    for (reason, expected_detail, document) in cases {
        assert_ne!(
            document, baseline,
            "«{reason}» tiene que cambiar el documento"
        );
        let error = validate_general_manifest_bytes(&path, &document)
            .expect_err("el lector tiene que rechazar la corrupción");
        // El motivo importa: rechazar por una guarda ajena no protege a la
        // guarda que la corrupción dice proteger.
        assert_eq!(
            error.context().get("detail").map(String::as_str),
            Some(expected_detail),
            "«{reason}» se rechazó por otra razón"
        );
    }
    fs::remove_dir_all(&published.root).ok();
}

/// **Una configuración efectiva con loci exige un manifiesto que los
/// registre.**
///
/// *Qué hace que pueda fallar:* se degrada el manifiesto real a su versión
/// anterior sin sección, que es exactamente el par incoherente que antes nadie
/// rechazaba, y se exige que la validación cruzada lo cace. El control es el
/// par intacto, que sí se acepta.
#[test]
fn an_effective_configuration_with_loci_demands_a_manifest_that_records_them() {
    let published = publish(UNEVEN_FREQUENCIES, 2);
    let manifest_path = published.destination.join("manifest.json");
    let effective_path = published.destination.join("effective-config.json");
    let manifest_bytes = fs::read(&manifest_path).expect("manifiesto");
    let effective_bytes = fs::read(&effective_path).expect("configuración efectiva");
    validate_effective_config_against_manifest_bytes(
        &effective_path,
        &effective_bytes,
        &manifest_path,
        &manifest_bytes,
    )
    .expect("el par intacto se acepta");

    let mut degraded: Value = serde_json::from_slice(&manifest_bytes).expect("JSON");
    degraded["schema_version"] = serde_json::json!("4");
    degraded
        .as_object_mut()
        .expect("raíz")
        .remove("repeat_loci");
    let degraded = serde_json::to_vec(&degraded).expect("JSON degradado");
    assert!(
        validate_effective_config_against_manifest_bytes(
            &effective_path,
            &effective_bytes,
            &manifest_path,
            &degraded,
        )
        .is_err(),
        "una configuración efectiva con loci no puede emparejarse con un manifiesto sin ellos"
    );
    fs::remove_dir_all(&published.root).ok();
}

/// **Un escenario sin loci publica el manifiesto de siempre, sin sección y en
/// su versión anterior.**
#[test]
fn a_scenario_without_loci_keeps_the_previous_manifest_version() {
    let published = publish("", 2);
    assert_eq!(published.manifest["schema_version"], "4");
    assert!(
        published.manifest.get("repeat_loci").is_none(),
        "sin loci no hay sección que publicar"
    );
    assert_external_schema(
        &published.destination.join("manifest.json"),
        "manifest-general-v4.schema.json",
    );
    fs::remove_dir_all(&published.root).ok();
}
