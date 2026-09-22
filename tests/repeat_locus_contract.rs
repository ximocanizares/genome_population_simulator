//! F6.1.1.1: contrato de entrada de los loci de repetición en tándem —
//! documento general **v4** y configuración efectiva **v6**.
//!
//! La forma es **declaración compacta, almacenamiento exacto**: por locus se
//! declara un puñado de alelos —motivo, número de unidades y sus interrupciones
//! si las tiene, o su secuencia exacta— y sus **frecuencias**, y el programa
//! materializa la secuencia. Es lo que resuelve la contradicción C-5: escribir
//! la cadena de bases de cada alelo de cada individuo reventaría el techo de
//! 65536 bytes del fichero general, mientras que un puñado de definiciones por
//! locus **no crece con el número de individuos**.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveConfig, load_and_validate_config};
use genome_population_simulator::execution::execute_validated_config;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn header(version: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "{version}"
seed = "{SEED}"

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
"#
    )
}

const MINIMAL_LOCUS: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 130
motif = "CAG"

[[reference.repeat_loci.alleles]]
units = 10
frequency = 0.6

[[reference.repeat_loci.alleles]]
units = 8
frequency = 0.4
"#;

fn write_config(body: &str) -> PathBuf {
    let path = temporary("toml");
    fs::write(&path, body).expect("documento general");
    path
}

fn load(
    body: &str,
) -> Result<EffectiveConfig, genome_population_simulator::diagnostics::Diagnostic> {
    let path = write_config(body);
    let result = load_and_validate_config(&path);
    fs::remove_file(&path).ok();
    result
}

/// La pareja mínima válida se admite, el locus se resuelve y la configuración
/// efectiva se publica en la versión **6**, con los mínimos de la guarda.
#[test]
fn a_minimal_v4_pair_is_admitted_and_publishes_effective_config_v6() {
    let effective = load(&(header("4") + MINIMAL_LOCUS)).expect("pareja mínima válida");
    assert_eq!(effective.repeat_loci().len(), 1);
    let locus = &effective.repeat_loci()[0];
    assert_eq!(locus.id(), "STR-1");
    assert_eq!(
        locus.single_motif().expect("un solo motivo").as_text(),
        "CAG"
    );
    assert_eq!(locus.max_units(), 200, "tope por omisión");
    assert_eq!(locus.alleles().len(), 2);
    assert_eq!(locus.alleles()[0].tract().as_text(), "CAG".repeat(10));
    assert_eq!(locus.alleles()[1].tract().as_text(), "CAG".repeat(8));
    assert!((locus.alleles()[0].frequency() - 0.6).abs() < 1e-12);

    let document = published_effective_config(&(header("4") + MINIMAL_LOCUS));
    assert_eq!(document["schema_version"], "6");
    assert_eq!(document["repeat_admission"]["min_purity"], 0.8);
    assert_eq!(document["repeat_admission"]["min_longest_run_units"], 3);
    assert_eq!(document["repeat_loci"][0]["alleles"][0]["units"], 10.0);
    assert_eq!(document["repeat_loci"][0]["alleles"][0]["purity"], 1.0);
    assert_eq!(
        document["repeat_loci"][0]["motifs"],
        serde_json::json!(["CAG"])
    );
    assert_eq!(
        document["repeat_loci"][0]["alleles"][0]["repeat_sequence_count"], 1,
        "un locus de un solo motivo publica exactamente un bloque"
    );
}

/// **Un fichero general v2 sin loci sigue publicando `effective-config-v4`, no
/// v6.** Es la regla que impide que la capacidad nueva mueva un byte de un
/// escenario que no la usa.
#[test]
fn a_general_v2_document_without_loci_still_publishes_effective_config_v4() {
    let effective = load(&header("2")).expect("documento v2");
    assert!(effective.repeat_loci().is_empty());
    let document = published_effective_config(&header("2"));
    assert_eq!(document["schema_version"], "4");
    assert!(
        document.get("repeat_loci").is_none(),
        "un documento v2 no publica la sección de repeticiones"
    );
    assert!(document.get("repeat_admission").is_none());
}

/// Y un documento **v4 sin loci** publica v6 igualmente: manda la versión del
/// documento, no la presencia de loci. Sin esta prueba, un lector que decidiera
/// por la presencia de loci pasaría la anterior sin despeinarse.
#[test]
fn a_general_v4_document_without_loci_still_publishes_effective_config_v6() {
    let effective = load(&header("4")).expect("documento v4 sin loci");
    assert!(effective.repeat_loci().is_empty());
    let document = published_effective_config(&header("4"));
    assert_eq!(document["schema_version"], "6");
    assert_eq!(document["repeat_loci"].as_array().map(Vec::len), Some(0));
}

/// Un locus declarado en un documento **v2** se rechaza como campo desconocido:
/// la capacidad pertenece al v4 y el v2 no la conoce.
#[test]
fn a_repeat_locus_in_a_v2_document_is_rejected_as_an_unknown_field() {
    let error = load(&(header("2") + MINIMAL_LOCUS)).expect_err("v2 no conoce las repeticiones");
    assert_eq!(error.code(), "GPS-CONF-006");
    assert_eq!(
        error.context().get("field").map(String::as_str),
        Some("repeat_loci")
    );
}

/// Una clave desconocida dentro del locus se rechaza: el objeto es cerrado.
#[test]
fn an_unknown_key_inside_a_repeat_locus_is_rejected() {
    let body = header("4") + MINIMAL_LOCUS;
    let error = load(&body.replace("motif = \"CAG\"", "motif = \"CAG\"\nrepeat_unit = \"CAG\""))
        .expect_err("clave desconocida");
    assert_eq!(error.code(), "GPS-CONF-006");
    assert_eq!(
        error.context().get("field").map(String::as_str),
        Some("repeat_unit")
    );
}

/// Un identificador de locus repetido se rechaza nombrándolo.
#[test]
fn a_repeated_locus_identifier_is_rejected() {
    let body = header("4") + MINIMAL_LOCUS + MINIMAL_LOCUS;
    let error = load(&body).expect_err("identificador repetido");
    assert!(
        format!("{error:?}").contains("identificador de locus repetido: STR-1"),
        "{error:?}"
    );
}

/// Frecuencias que no suman 1 dentro de la tolerancia se rechazan **nombrando
/// el locus**.
#[test]
fn frequencies_outside_the_tolerance_are_rejected_naming_the_locus() {
    let body = (header("4") + MINIMAL_LOCUS).replace("frequency = 0.4", "frequency = 0.5");
    let error = load(&body).expect_err("frecuencias que no suman 1");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-1"), "{rendered}");
    assert!(rendered.contains("suman 1.1"), "{rendered}");

    // Y la tolerancia **declarada** manda sobre la de omisión. El desajuste
    // elegido, 5e-8, cae **entre** las dos: la de omisión es 1e-9 y lo rechaza;
    // la declarada, 1e-6, lo admite.
    //
    // *Corregido tras la revisión adversarial:* el fixture anterior usaba un
    // desajuste de 1e-10, por debajo de la tolerancia de omisión, así que el
    // documento se habría admitido igual **sin declarar nada** y la declaración
    // no se ejercitaba.
    let skewed = (header("4") + MINIMAL_LOCUS).replace("frequency = 0.4", "frequency = 0.40000005");
    let with_default_tolerance =
        load(&skewed).expect_err("con la tolerancia de omisión, 5e-8 se aparta demasiado de 1");
    assert!(
        format!("{with_default_tolerance:?}").contains("locus STR-1"),
        "{with_default_tolerance:?}"
    );
    let with_declared_tolerance =
        skewed.replace("seed = ", "normalization_tolerance = 0.000001\nseed = ");
    load(&with_declared_tolerance).expect("con la tolerancia declarada, 5e-8 entra");
}

/// Una frecuencia negativa se rechaza, **nombrando el alelo que la declara**.
///
/// *Corregido tras la revisión adversarial.* El fixture anterior ponía 1,4 en el
/// primer alelo y −0,4 en el segundo, y el recorrido moría en el 1,4 sin llegar
/// nunca al negativo: retirando la cota inferior, la prueba seguía pasando. Aquí
/// el primer alelo es válido y **solo el segundo** es negativo, de modo que la
/// única comprobación que puede rechazarlo es la cota inferior.
#[test]
fn a_negative_frequency_is_rejected_naming_its_allele() {
    // El primer alelo queda **válido** a propósito: si también fuera inválido, el
    // recorrido moriría en él y el negativo no se llegaría a mirar, que es justo
    // el defecto que esta corrección repara.
    let body = (header("4") + MINIMAL_LOCUS)
        .replace("units = 8\nfrequency = 0.4", "units = 8\nfrequency = -0.4");
    let error = load(&body).expect_err("frecuencia negativa");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("fuera de [0,1]"), "{rendered}");
    assert!(
        rendered.contains("el alelo 1"),
        "el diagnóstico tiene que nombrar el alelo negativo, no el primero: {rendered}"
    );
}

/// Un alelo sin frecuencia se rechaza, y una frecuencia sin alelo —una entrada
/// que declara frecuencia y no dice qué alelo es— también.
#[test]
fn an_allele_without_frequency_and_a_frequency_without_allele_are_rejected() {
    let missing_frequency =
        (header("4") + MINIMAL_LOCUS).replace("units = 8\nfrequency = 0.4", "units = 8");
    let error = load(&missing_frequency).expect_err("alelo sin frecuencia");
    assert!(format!("{error:?}").contains("frequency"), "{error:?}");

    let frequency_without_allele =
        (header("4") + MINIMAL_LOCUS).replace("units = 8\nfrequency = 0.4", "frequency = 0.4");
    let error = load(&frequency_without_allele).expect_err("frecuencia sin alelo");
    assert!(
        format!("{error:?}").contains("debe declarar units o sequence"),
        "{error:?}"
    );

    // Y declarar las dos formas a la vez también se rechaza: son excluyentes.
    let both = (header("4") + MINIMAL_LOCUS).replace(
        "units = 8\nfrequency = 0.4",
        "units = 8\nsequence = \"CAGCAGCAG\"\nfrequency = 0.4",
    );
    let error = load(&both).expect_err("las dos formas a la vez");
    assert!(
        format!("{error:?}").contains("debe declarar units o sequence"),
        "{error:?}"
    );
}

/// Declarar una **pureza objetivo** se rechaza nombrando el locus y el
/// parámetro. La pureza es siempre una salida calculada.
#[test]
fn a_declared_target_purity_is_rejected_naming_the_locus_and_the_parameter() {
    for parameter in ["purity", "target_purity", "purity_target"] {
        let body = (header("4") + MINIMAL_LOCUS).replace(
            "motif = \"CAG\"",
            &format!("motif = \"CAG\"\n{parameter} = 0.9"),
        );
        let error = load(&body).expect_err("pureza declarada");
        let rendered = format!("{error:?}");
        assert!(rendered.contains("locus STR-1"), "{rendered}");
        assert!(rendered.contains(parameter), "{rendered}");
        assert!(rendered.contains("no se declara"), "{rendered}");
    }
}

/// El tope propio se lee de la configuración, no está codificado: con el tope
/// en 60 unidades un alelo de 60 se admite y uno de 61 se rechaza nombrando el
/// locus; **con el tope en 20, el de 60 pasa a rechazarse**.
#[test]
fn the_declared_cap_is_read_and_not_hardcoded() {
    let with_cap = |cap: u64, units: u64| {
        format!(
            r#"{}
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 130
motif = "CA"
max_units = {cap}

[[reference.repeat_loci.alleles]]
units = {units}
frequency = 1.0
"#,
            header("4")
        )
    };
    load(&with_cap(60, 60)).expect("sesenta unidades con el tope en sesenta");
    let error = load(&with_cap(60, 61)).expect_err("sesenta y una con el tope en sesenta");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-1"), "{rendered}");
    assert!(rendered.contains("rebasa el tope de 60"), "{rendered}");

    let error = load(&with_cap(20, 60)).expect_err("sesenta unidades con el tope en veinte");
    assert!(
        format!("{error:?}").contains("rebasa el tope de 20"),
        "el tope se lee de la configuración y no de una constante"
    );

    // Un alelo de 30 unidades de un dinucleótido son 60 pares de base, que el
    // umbral de 50 pb retirado habría rechazado. Se admite.
    load(&with_cap(200, 30)).expect("sesenta pares de base de un dinucleótido");
}

/// Una clave `max_bases` se rechaza **como campo desconocido**, porque el
/// objeto es cerrado.
///
/// *Corregido respecto de una redacción anterior:* no es realizable comprobar
/// que «un tope expresado en pares de base se rechaza», porque un entero no
/// lleva unidad consigo y el programa no puede distinguir 60 unidades de 60
/// pares de base. Lo único comprobable es el **nombre de la clave**.
#[test]
fn a_max_bases_key_is_rejected_as_an_unknown_field() {
    let body = (header("4") + MINIMAL_LOCUS)
        .replace("motif = \"CAG\"", "motif = \"CAG\"\nmax_bases = 400");
    let error = load(&body).expect_err("clave max_bases");
    assert_eq!(error.code(), "GPS-CONF-006");
    assert_eq!(
        error.context().get("field").map(String::as_str),
        Some("max_bases")
    );
}

/// **Contradicción C-5 resuelta.** Un escenario con 500 individuos y 5 loci
/// cabe holgadamente en el techo de 65536 bytes, porque el coste del fichero ya
/// no crece con el número de individuos: los alelos se declaran una vez por
/// locus y el simulador los sortea.
#[test]
fn five_hundred_individuals_and_five_loci_fit_well_inside_the_ceiling() {
    let mut body = header("4");
    let keys: Vec<String> = (1..=500).map(|index| format!("\"f{index:03}\"")).collect();
    body.push_str(&format!(
        r#"
[[initial_populations]]
population_id = "P1"
individual_count = 500
ploidy = 2
genotype_model = "hwe"
founder_keys = [{}]
"#,
        keys.join(", ")
    ));
    let without_loci = body.len();
    body.push_str(&five_loci_block());
    let loci_cost = body.len() - without_loci;
    assert!(
        body.len() < 65536,
        "el documento mide {} bytes y rebasa el techo de 65536",
        body.len()
    );
    let effective = load(&body).expect("quinientos individuos y cinco loci");
    assert_eq!(effective.repeat_loci().len(), 5);
    assert_eq!(effective.initial_populations()[0].individual_count(), 500);

    // Y lo que de verdad resuelve la contradicción C-5: el coste de los loci
    // **no crece con el número de individuos**.
    //
    // *Corregido tras la revisión adversarial.* La versión anterior pegaba el
    // mismo texto y comprobaba que medía lo que acababa de pegar —una identidad
    // escrita a mano junto al bucle— y además **nunca cargaba** el documento
    // grande. Ahora el bloque se genera de nuevo con la misma función y el
    // documento de mil individuos se carga y se comprueba.
    let mut larger = header("4");
    let more: Vec<String> = (1..=1000).map(|index| format!("\"f{index:04}\"")).collect();
    larger.push_str(&format!(
        r#"
[[initial_populations]]
population_id = "P1"
individual_count = 1000
ploidy = 2
genotype_model = "hwe"
founder_keys = [{}]
"#,
        more.join(", ")
    ));
    let larger_without_loci = larger.len();
    larger.push_str(&five_loci_block());
    assert_eq!(
        larger.len() - larger_without_loci,
        loci_cost,
        "el bloque de repeticiones no puede crecer con el número de individuos"
    );
    assert!(
        larger.len() < 65536,
        "con mil individuos el documento mide {} bytes y rebasa el techo",
        larger.len()
    );
    let larger_effective = load(&larger).expect("mil individuos y cinco loci");
    assert_eq!(larger_effective.repeat_loci().len(), 5);
    assert_eq!(
        larger_effective.initial_populations()[0].individual_count(),
        1000
    );
}

/// Bloque de cinco loci, generado por **una sola función**, de modo que el
/// documento pequeño y el grande no comparten un texto copiado sino su origen.
fn five_loci_block() -> String {
    let mut block = String::new();
    for locus in 1..=5 {
        block.push_str(&format!(
            r#"
[[reference.repeat_loci]]
id = "STR-{locus}"
sequence = "chr1"
start = {}
end = {}
motif = "CAG"

[[reference.repeat_loci.alleles]]
units = 10
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 12
frequency = 0.3

[[reference.repeat_loci.alleles]]
units = 8
frequency = 0.2
"#,
            100 * locus + 1,
            100 * locus + 30
        ));
    }
    block
}

/// Un documento legítimo de exactamente 65536 bytes se admite y el de 65537 se
/// rechaza. El relleno es un comentario TOML, de modo que el documento sigue
/// siendo válido en las dos medidas y lo único que cambia es su tamaño.
#[test]
fn a_legitimate_document_of_65536_bytes_is_admitted_and_65537_is_rejected() {
    let base = header("4") + MINIMAL_LOCUS;
    let pad = |target: usize| {
        let mut body = base.clone();
        body.push_str("\n# ");
        while body.len() < target {
            body.push('x');
        }
        body.truncate(target);
        body
    };
    let admitted = pad(65536);
    assert_eq!(admitted.len(), 65536);
    load(&admitted).expect("documento de 65536 bytes");

    let rejected = pad(65537);
    assert_eq!(rejected.len(), 65537);
    let error = load(&rejected).expect_err("documento de 65537 bytes");
    assert_eq!(error.code(), "GPS-INPUT-001");
}

/// La configuración efectiva publicada valida contra su esquema con el
/// validador externo, y la de un documento v2 contra el suyo.
#[test]
fn the_published_effective_configs_validate_against_their_schemas_externally() {
    for (body, schema) in [
        (
            header("4") + MINIMAL_LOCUS,
            "effective-config-v6.schema.json",
        ),
        (header("2"), "effective-config-v4.schema.json"),
    ] {
        let document = published_effective_config(&body);
        let path = temporary("json");
        fs::write(
            &path,
            serde_json::to_vec_pretty(&document).expect("serialización"),
        )
        .expect("documento");
        assert_external_schema(&path, schema);
        fs::remove_file(&path).ok();
    }
}

/// El documento general válido se comprueba también contra el esquema de
/// autoría v4, convertido a JSON, para que la vista de autoría y el lector no
/// puedan divergir en silencio.
#[test]
fn the_general_v4_authoring_schema_accepts_the_declared_document() {
    let document = serde_json::json!({
        "config_kind": "general",
        "schema_version": "4",
        "seed": SEED,
        "reference": {
            "kind": "simulated",
            "id": "ref-main",
            "composition": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            "sequences": [{"id": "chr1", "class": "chromosome", "length": 1000}],
            "repeat_loci": [{
                "id": "STR-1",
                "sequence": "chr1",
                "start": 101,
                "end": 130,
                "motif": "CAG",
                "alleles": [
                    {"units": 10, "frequency": 0.6},
                    {"units": 8, "frequency": 0.4, "interruptions": [{"offset": 6, "bases": "TT"}]}
                ]
            }]
        }
    });
    let path = temporary("json");
    fs::write(
        &path,
        serde_json::to_vec_pretty(&document).expect("serialización"),
    )
    .expect("documento");
    assert_external_schema(&path, "config-general-v4.schema.json");
    fs::remove_file(&path).ok();
}

/// Un `huge` compartido por las dos pruebas del tope: un locus con un único
/// alelo cuyo número de unidades se elige desde fuera.
fn huge_allele_document(units: &str) -> String {
    format!(
        r#"{}
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 130
motif = "CAG"

[[reference.repeat_loci.alleles]]
units = {units}
frequency = 1.0
"#,
        header("4")
    )
}

/// **El tope se decide sobre lo declarado, y el motivo del rechazo lo demuestra.**
///
/// Con mil trillones de unidades, materializar el tracto ni siquiera cabe en el
/// espacio de direcciones: la comparación **tardía** falla primero en el
/// constructor y devuelve *otro* diagnóstico. El motivo del rechazo es, por
/// tanto, un observable exacto que separa las dos implementaciones sin depender
/// de ningún reloj ni de cuánta memoria tenga la máquina.
#[test]
fn the_cap_is_decided_on_the_declaration_and_not_on_the_materialized_tract() {
    let error = load(&huge_allele_document("10000000000000000000")).expect_err("alelo desmesurado");
    let rendered = format!("{error:?}");
    assert!(
        rendered.contains("rebasa el tope de 200"),
        "el rechazo tiene que ser el del tope, no el del espacio de direcciones: {rendered}"
    );
    assert!(
        !rendered.contains("address space"),
        "materializar antes de comparar cambia el motivo del rechazo: {rendered}"
    );
}

/// **Y no se materializa nada por el camino.**
///
/// Trescientos millones de unidades sobre un motivo de tres bases son
/// novecientos megabytes de tracto. Es el caso que la revisión adversarial midió
/// en **26 segundos y 3,2 GB** con la comparación tardía; con la temprana se
/// resuelve sin asignar nada.
///
/// Esta prueba **usa el reloj**, y por eso va separada de la anterior: si
/// estuvieran juntas, la aserción determinista mataría al mutante antes de
/// llegar aquí y la cota de tiempo no se ejercitaría nunca. El margen es enorme
/// a propósito; no mide rendimiento, solo detecta que se ha materializado.
///
/// *Aviso:* si alguien devuelve la comparación a su sitio anterior, esta prueba
/// intentará reservar novecientos megabytes antes de fallar. Es el defecto, no
/// la prueba.
#[test]
fn rejecting_an_oversized_allele_allocates_nothing() {
    let started = std::time::Instant::now();
    let error = load(&huge_allele_document("300000000")).expect_err("alelo desmesurado");
    let elapsed = started.elapsed();
    assert!(
        format!("{error:?}").contains("rebasa el tope de 200"),
        "{error:?}"
    );
    assert!(
        elapsed < std::time::Duration::from_secs(10),
        "rechazar novecientos megabytes declarados tardó {elapsed:?}: se está materializando el tracto antes de compararlo"
    );
}

/// **Los mínimos de admisión se aplican al alelo declarado, no solo se
/// publican.**
///
/// Es el defecto que la revisión adversarial llamó «publicar el cartel de una
/// guarda que no existe»: la configuración efectiva anunciaba un umbral de
/// pureza de 0,8 y tres líneas más abajo publicaba alelos de pureza 0
/// presentados como unidades de una repetición. Quien leyera ese documento
/// concluiría que la guarda corrió y que esos alelos pasaron.
#[test]
fn the_admission_minima_reject_a_declared_allele_that_does_not_meet_them() {
    let with_sequence = |sequence: &str, minimum: &str| {
        format!(
            r#"{}
[reference.repeat_admission]
{minimum}

[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 130
motif = "CAG"

[[reference.repeat_loci.alleles]]
sequence = "{sequence}"
frequency = 1.0
"#,
            header("4")
        )
    };

    // Doce bases que no contienen ni una copia de `CAG`: pureza cero.
    let error = load(&with_sequence("TTTTTTTTTTTT", "min_purity = 0.8"))
        .expect_err("alelo sin ni una copia del motivo");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-1"), "{rendered}");
    assert!(rendered.contains("pureza 0"), "{rendered}");
    assert!(rendered.contains("mínimo declarado es 0.8"), "{rendered}");

    // Y con doce indeterminadas, que es el otro caso que la revisión reprodujo.
    let error = load(&with_sequence("NNNNNNNNNNNN", "min_purity = 0.8"))
        .expect_err("alelo de bases indeterminadas");
    assert!(format!("{error:?}").contains("pureza 0"), "{error:?}");

    // El umbral se **lee**: bajándolo a cero, el mismo alelo pasa la pureza y lo
    // que lo rechaza pasa a ser el mínimo de tirada, con su propio mensaje.
    let error = load(&with_sequence(
        "TTTTTTTTTTTT",
        "min_purity = 0.0\nmin_longest_run_units = 3",
    ))
    .expect_err("tirada más larga de cero unidades");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("tirada más larga de 0"), "{rendered}");
    assert!(rendered.contains("mínimo declarado es 3"), "{rendered}");

    // Y con los dos mínimos a su valor más permisivo el alelo entra: es lo que
    // demuestra que lo que rechazaba era el umbral y no otra comprobación.
    load(&with_sequence(
        "CAGCAGCAGTTT",
        "min_purity = 0.5\nmin_longest_run_units = 3",
    ))
    .expect("tres copias de CAG y una interrupción de tres bases");

    // Un tracto real impuro por encima del mínimo se admite. Es el caso que una
    // guarda de pureza 1 rechazaría por error.
    load(&with_sequence(
        "CAGCAGCAGTTCAGCAGCAG",
        "min_purity = 0.8\nmin_longest_run_units = 3",
    ))
    .expect("dieciocho de veinte bases en copias completas");
}

/// **El identificador del locus es tipado y se valida al leerlo.**
///
/// Sin esto, un identificador con espacios o con signos se admitía, se ejecutaba
/// y se publicaba, y solo el validador externo de la configuración efectiva lo
/// cazaba después contra su propio patrón: un artefacto que incumple su propio
/// contrato, descubierto únicamente si alguien pasa un validador que en
/// producción no corre.
#[test]
fn a_locus_identifier_with_invalid_syntax_is_rejected_when_it_is_read() {
    for bad in ["STR 1", "STR-1 -- ¡ojo!", "-STR1", ""] {
        let body =
            (header("4") + MINIMAL_LOCUS).replace("id = \"STR-1\"", &format!("id = \"{bad}\""));
        let error = load(&body).expect_err("identificador inválido");
        assert_eq!(
            error.code(),
            "GPS-DOMAIN-001",
            "identificador rechazado: {bad:?} — {error:?}"
        );
        assert_eq!(
            error.context().get("category").map(String::as_str),
            Some("repeat_locus")
        );
    }
    // Control positivo: el identificador legítimo del fixture sigue entrando.
    load(&(header("4") + MINIMAL_LOCUS)).expect("identificador válido");
}

/// Un intervalo más corto que su motivo se rechaza, y la guarda vive **en el
/// tipo del dominio**, no reescrita en la frontera de configuración.
///
/// Se añade porque al unificar las dos representaciones del locus la guarda se
/// movió, y una guarda que cambia de sitio sin prueba es una guarda que puede
/// desaparecer sin que nadie se entere.
#[test]
fn an_interval_shorter_than_its_motif_is_rejected() {
    let body = (header("4") + MINIMAL_LOCUS)
        .replace("end = 130", "end = 104")
        .replace("motif = \"CAG\"", "motif = \"CAGTT\"");
    let error = load(&body).expect_err("intervalo más corto que el motivo");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-1"), "{rendered}");
    assert!(
        rendered.contains("is shorter than the 5-base motif"),
        "{rendered}"
    );

    // Justo con la longitud del motivo se admite: es lo que demuestra que la
    // comparación es «más corto que», y no «más corto o igual».
    let exact = (header("4") + MINIMAL_LOCUS)
        .replace("end = 130", "end = 105")
        .replace("motif = \"CAG\"", "motif = \"CAGTT\"");
    load(&exact).expect("un intervalo de exactamente una copia del motivo");
}

/// Dos loci con intervalos solapados se rechazan nombrando el que ya estaba.
#[test]
fn two_loci_with_overlapping_intervals_are_rejected() {
    let body = header("4")
        + MINIMAL_LOCUS
        + &MINIMAL_LOCUS
            .replace("id = \"STR-1\"", "id = \"STR-2\"")
            .replace("start = 101", "start = 120")
            .replace("end = 130", "end = 150");
    let error = load(&body).expect_err("intervalos solapados");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-2"), "{rendered}");
    assert!(
        rendered.contains("se solapa con el del locus STR-1"),
        "{rendered}"
    );

    // Pegados pero sin solapar: se admiten. Es lo que demuestra que la
    // comprobación es de solapamiento y no de vecindad.
    let adjacent = header("4")
        + MINIMAL_LOCUS
        + &MINIMAL_LOCUS
            .replace("id = \"STR-1\"", "id = \"STR-2\"")
            .replace("start = 101", "start = 131")
            .replace("end = 130", "end = 160");
    load(&adjacent).expect("dos loci contiguos que no se solapan");
}

/// La auditoría de la normalización se publica **recomputada**, y las dos
/// frecuencias —la declarada y la efectiva— se distinguen.
///
/// Con frecuencias que no suman exactamente uno dentro de la tolerancia, la
/// normalizada **no es la que el usuario escribió**; publicar solo una bajo un
/// nombre que sugiere la otra es lo que DTC-049 obliga a evitar.
#[test]
fn the_frequency_normalization_audit_is_published_and_recomputed() {
    let skewed = (header("4") + MINIMAL_LOCUS)
        .replace("frequency = 0.4", "frequency = 0.40000005")
        .replace("seed = ", "normalization_tolerance = 0.000001\nseed = ");
    let document = published_effective_config(&skewed);
    let audit = &document["repeat_loci"][0]["frequency_normalization"];
    assert_eq!(audit["sum_algorithm"], "neumaier_sum_v1");
    assert_eq!(audit["tolerance"], 0.000001);

    let original_sum = audit["original_sum"].as_f64().expect("suma original");
    let original_delta = audit["original_delta"].as_f64().expect("delta original");
    let effective_sum = audit["effective_sum"].as_f64().expect("suma efectiva");
    let effective_delta = audit["effective_delta"].as_f64().expect("delta efectiva");
    assert!(
        (original_sum - 1.00000005).abs() < 1e-12,
        "la suma original es la de lo declarado: {original_sum}"
    );
    assert!(
        original_delta > 1e-9,
        "la desviación original tiene que ser la real, no cero: {original_delta}"
    );
    assert!(
        effective_delta < original_delta,
        "normalizar tiene que acercar la suma a uno: {effective_delta} frente a {original_delta}"
    );
    assert!((effective_sum - 1.0).abs() < 1e-12);

    let alleles = document["repeat_loci"][0]["alleles"]
        .as_array()
        .expect("alelos");
    let declared = alleles[1]["declared_frequency"]
        .as_f64()
        .expect("declarada");
    let effective = alleles[1]["effective_frequency"]
        .as_f64()
        .expect("efectiva");
    assert!(
        (declared - 0.40000005).abs() < 1e-15,
        "la declarada es literalmente la que se escribió: {declared}"
    );
    assert_ne!(
        declared, effective,
        "con una suma declarada distinta de uno, las dos frecuencias no coinciden"
    );
    let sum: f64 = alleles
        .iter()
        .map(|allele| allele["effective_frequency"].as_f64().expect("efectiva"))
        .sum();
    assert!((sum - 1.0).abs() < 1e-12);
}

/// El esquema de autoría v4 rechaza lo mismo que el lector.
///
/// La revisión adversarial encontró una divergencia real: un alelo que declaraba
/// **secuencia e interrupciones** lo admitía el esquema y lo rechazaba el lector.
/// Aquí se comprueban los tres casos, con su control positivo.
#[test]
fn the_authoring_schema_rejects_what_the_reader_rejects() {
    let document = |allele: serde_json::Value| {
        serde_json::json!({
            "config_kind": "general",
            "schema_version": "4",
            "reference": {
                "kind": "simulated",
                "id": "ref-main",
                "composition": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                "sequences": [{"id": "chr1", "class": "chromosome", "length": 1000}],
                "repeat_loci": [{
                    "id": "STR-1",
                    "sequence": "chr1",
                    "start": 101,
                    "end": 130,
                    "motif": "CAG",
                    "alleles": [allele]
                }]
            }
        })
    };
    let schema = "config-general-v4.schema.json";

    // Control positivo: sin él, un esquema que rechazase todo pasaría los tres
    // casos negativos sin discriminar nada.
    assert!(validates(
        &document(serde_json::json!({"units": 10, "frequency": 1.0})),
        schema
    ));
    assert!(validates(
        &document(serde_json::json!({
            "units": 10,
            "frequency": 1.0,
            "interruptions": [{"offset": 6, "bases": "TT"}]
        })),
        schema
    ));
    assert!(validates(
        &document(serde_json::json!({"sequence": "CAGCAGCAG", "frequency": 1.0})),
        schema
    ));

    assert!(
        !validates(
            &document(serde_json::json!({
                "sequence": "CAGCAGCAG",
                "frequency": 1.0,
                "interruptions": [{"offset": 6, "bases": "TT"}]
            })),
            schema
        ),
        "una secuencia exacta ya lleva sus interrupciones dentro: declarar las dos cosas es la divergencia que la revisión encontró"
    );
    assert!(
        !validates(
            &document(serde_json::json!({"units": 10, "sequence": "CAGCAGCAG", "frequency": 1.0})),
            schema
        ),
        "las dos formas son excluyentes"
    );
    assert!(
        !validates(&document(serde_json::json!({"frequency": 1.0})), schema),
        "una frecuencia sin alelo no describe ningún alelo"
    );
    assert!(
        !validates(
            &document(serde_json::json!({"units": 10, "frequency": 1.0, "purity": 0.9})),
            schema
        ),
        "la pureza no se declara ni en el esquema"
    );
}

/// El esquema de autoría también cierra las dos formas de declarar el motivo y
/// admite la lista de unidades por bloque.
#[test]
fn the_authoring_schema_closes_the_two_motif_forms() {
    let locus = |extra: serde_json::Value| {
        let mut value = serde_json::json!({
            "id": "STR-1",
            "sequence": "chr1",
            "start": 101,
            "end": 160,
            "alleles": [{"units": 10, "frequency": 1.0}]
        });
        for (key, item) in extra.as_object().expect("objeto") {
            value[key] = item.clone();
        }
        serde_json::json!({
            "config_kind": "general",
            "schema_version": "4",
            "reference": {
                "kind": "simulated",
                "id": "ref-main",
                "composition": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                "sequences": [{"id": "chr1", "class": "chromosome", "length": 1000}],
                "repeat_loci": [value]
            }
        })
    };
    let schema = "config-general-v4.schema.json";

    // Controles positivos: las dos formas legítimas.
    assert!(validates(
        &locus(serde_json::json!({"motif": "CAG"})),
        schema
    ));
    assert!(validates(
        &locus(serde_json::json!({
            "motifs": ["TCTA", "TCTG", "TCTA"],
            "alleles": [{"units": [5, 5, 5], "frequency": 1.0}]
        })),
        schema
    ));
    assert!(validates(
        &locus(serde_json::json!({"motif": "CAG", "min_purity": 0.6, "min_longest_run_units": 4})),
        schema
    ));

    assert!(
        !validates(
            &locus(serde_json::json!({"motif": "CAG", "motifs": ["CAG"]})),
            schema
        ),
        "motif y motifs son excluyentes"
    );
    assert!(
        !validates(&locus(serde_json::json!({})), schema),
        "un locus sin ningún motivo no describe nada"
    );
    assert!(
        !validates(&locus(serde_json::json!({"motifs": []})), schema),
        "una lista de motivos vacía no describe nada"
    );
    assert!(
        !validates(
            &locus(serde_json::json!({"motif": "CAG", "min_longest_run_units": 0})),
            schema
        ),
        "la tirada mínima es un entero positivo"
    );
}

/// El esquema de la configuración efectiva v6 exige la sección de repeticiones
/// y su auditoría, y rechaza un documento al que se las quiten.
#[test]
fn the_effective_config_v6_schema_requires_the_repeat_section() {
    let valid = published_effective_config(&(header("4") + MINIMAL_LOCUS));
    let schema = "effective-config-v6.schema.json";
    assert!(validates(&valid, schema));

    let mut without_admission = valid.clone();
    without_admission
        .as_object_mut()
        .expect("raíz")
        .remove("repeat_admission");
    assert!(
        !validates(&without_admission, schema),
        "los mínimos publicados son obligatorios en la versión 6"
    );

    let mut without_audit = valid.clone();
    without_audit["repeat_loci"][0]
        .as_object_mut()
        .expect("locus")
        .remove("frequency_normalization");
    assert!(
        !validates(&without_audit, schema),
        "la auditoría de la normalización es obligatoria"
    );

    let mut without_effective_frequency = valid.clone();
    without_effective_frequency["repeat_loci"][0]["alleles"][0]
        .as_object_mut()
        .expect("alelo")
        .remove("effective_frequency");
    assert!(
        !validates(&without_effective_frequency, schema),
        "las dos frecuencias son obligatorias"
    );
}

// ---------------------------------------------------------------------------
// Decisión 1 del responsable: varios motivos por locus.
// ---------------------------------------------------------------------------

/// Locus compuesto `(TCTA)5(TCTG)5(TCTA)5`, el caso biológico que motiva la
/// decisión. Sirve para la **admisión**, no para las listas: sus tres bloques
/// son indistinguibles entre sí.
fn composite_locus(extra: &str) -> String {
    format!(
        r#"{}
[[reference.repeat_loci]]
id = "STR-C"
sequence = "chr1"
start = 101
end = 160
motifs = ["TCTA", "TCTG", "TCTA"]
{extra}

[[reference.repeat_loci.alleles]]
units = [5, 5, 5]
frequency = 1.0
"#,
        header("4")
    )
}

/// Locus compuesto **saturado**: `(CAG)4(CA)5`.
///
/// *Añadido tras la segunda revisión adversarial.* El fixture de arriba es
/// **degenerado para las listas**: sus tres bloques tienen la misma longitud de
/// motivo, las mismas copias y las mismas bases, y la lista de motivos es
/// además simétrica, de modo que invertir, cruzar o colapsar las cuatro listas
/// **da exactamente el mismo documento**. El revisor lo demostró publicando la
/// longitud de motivo del primer bloque en todos: las 1198 pruebas pasaron.
///
/// Aquí las tres magnitudes son **todas distintas y asimétricas**:
///
/// | | bloque 0 | bloque 1 |
/// |---|---|---|
/// | motivo | `CAG` | `CA` |
/// | longitud del motivo | 3 | 2 |
/// | copias | 4 | 5 |
/// | bases | 12 | 10 |
///
/// Ninguna permutación ni colapso de las listas sobrevive a comprobarlas contra
/// estos valores literales.
fn saturated_locus() -> String {
    format!(
        r#"{}
[[reference.repeat_loci]]
id = "STR-S"
sequence = "chr1"
start = 101
end = 130
motifs = ["CAG", "CA"]

[[reference.repeat_loci.alleles]]
units = [4, 5]
frequency = 1.0
"#,
        header("4")
    )
}

/// **El caso que motiva la decisión, y que una implementación de un solo motivo
/// no puede pasar.**
///
/// El mismo tracto declarado con un solo motivo `TCTA` tiene pureza 40/60 y lo
/// rechaza el mínimo por omisión; declarado con sus tres bloques tiene pureza 1
/// y entra. La prueba comprueba **las dos mitades**: si solo comprobara que el
/// compuesto entra, un modelo que ignorase los bloques y midiera la pureza de
/// otra manera pasaría igual.
#[test]
fn a_composite_microsatellite_needs_its_blocks_to_be_admitted() {
    let effective = load(&composite_locus("")).expect("microsatélite compuesto");
    let locus = &effective.repeat_loci()[0];
    assert_eq!(locus.motifs().len(), 3);
    assert!(locus.single_motif().is_none());
    let annotation = locus.alleles()[0].tract().annotation();
    assert!((annotation.purity_as_f64() - 1.0).abs() < f64::EPSILON);
    assert_eq!(annotation.blocks().len(), 3);

    // El mismo tracto con un solo motivo: sesenta bases de las que solo cuarenta
    // caen en copias completas de `TCTA`.
    let as_one = format!(
        r#"{}
[[reference.repeat_loci]]
id = "STR-C"
sequence = "chr1"
start = 101
end = 160
motif = "TCTA"

[[reference.repeat_loci.alleles]]
sequence = "{}"
frequency = 1.0
"#,
        header("4"),
        "TCTA".repeat(5) + &"TCTG".repeat(5) + &"TCTA".repeat(5)
    );
    let error = load(&as_one).expect_err("leído con un solo motivo no llega al mínimo");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-C"), "{rendered}");
    assert!(rendered.contains("pureza 0.666"), "{rendered}");
}

/// **Las cuatro listas de bloque se comprueban contra sus valores literales.**
///
/// *Corregido tras la segunda revisión adversarial.* La versión anterior
/// comprobaba que «copias por longitud de motivo son las bases», y eso **es una
/// identidad circular**: el número de copias está *definido* como bases entre
/// longitud del motivo, así que se cumple siempre, con listas alineadas o no.
/// Aquí cada lista se compara con lo que tiene que valer, sobre un fixture donde
/// las tres magnitudes son distintas entre bloques.
#[test]
fn the_effective_config_publishes_the_blocks_as_parallel_lists() {
    let document = published_effective_config(&saturated_locus());
    let locus = &document["repeat_loci"][0];
    assert_eq!(locus["motifs"], serde_json::json!(["CAG", "CA"]));
    assert_eq!(locus["motif_lengths"], serde_json::json!([3, 2]));

    let allele = &locus["alleles"][0];
    assert_eq!(allele["repeat_sequence_count"], 2);
    assert_eq!(allele["block_motifs"], serde_json::json!(["CAG", "CA"]));
    assert_eq!(allele["block_motif_lengths"], serde_json::json!([3, 2]));
    assert_eq!(allele["block_unit_counts"], serde_json::json!([4.0, 5.0]));
    assert_eq!(allele["block_bases"], serde_json::json!([12, 10]));
    assert_eq!(allele["sequence"], "CAG".repeat(4) + &"CA".repeat(5));
    assert_eq!(allele["units"], 9.0);

    // Y el orden importa: invertir cualquiera de las listas produce un documento
    // distinto del que se publica, que es lo que el fixture degenerado no podía
    // demostrar.
    assert_ne!(allele["block_motif_lengths"], serde_json::json!([2, 3]));
    assert_ne!(allele["block_unit_counts"], serde_json::json!([5.0, 4.0]));
    assert_ne!(allele["block_bases"], serde_json::json!([10, 12]));
}

/// Las dos formas de declarar el motivo son excluyentes, y el número de
/// unidades tiene que tener un valor por bloque.
#[test]
fn the_two_motif_forms_are_exclusive_and_the_units_must_match_the_blocks() {
    let both = composite_locus("").replace(
        "motifs = [\"TCTA\", \"TCTG\", \"TCTA\"]",
        "motif = \"TCTA\"\nmotifs = [\"TCTA\", \"TCTG\", \"TCTA\"]",
    );
    assert!(
        format!("{:?}", load(&both).expect_err("motif y motifs a la vez"))
            .contains("son excluyentes")
    );

    let scalar = composite_locus("").replace("units = [5, 5, 5]", "units = 5");
    assert!(
        format!(
            "{:?}",
            load(&scalar).expect_err("un escalar con tres motivos")
        )
        .contains("un número de unidades por bloque"),
    );

    let short = composite_locus("").replace("units = [5, 5, 5]", "units = [5, 5]");
    assert!(
        format!(
            "{:?}",
            load(&short).expect_err("dos valores para tres bloques")
        )
        .contains("hace falta uno por bloque")
    );

    // Y la secuencia exacta no se admite con varios motivos: unas bases sueltas
    // no pueden decir dónde acaba un bloque y empieza el siguiente.
    let exact =
        composite_locus("").replace("units = [5, 5, 5]", "sequence = \"TCTATCTATCTATCTATCTA\"");
    assert!(
        format!(
            "{:?}",
            load(&exact).expect_err("secuencia exacta con tres motivos")
        )
        .contains("dónde acaba un bloque")
    );
}

/// Un locus de **un solo motivo** sigue admitiendo el escalar y se publica con
/// un único bloque. Es lo que garantiza que nada de lo anterior se movió.
#[test]
fn a_single_motif_locus_keeps_the_scalar_form() {
    let effective = load(&(header("4") + MINIMAL_LOCUS)).expect("un solo motivo");
    let locus = &effective.repeat_loci()[0];
    assert_eq!(locus.motifs().len(), 1);
    assert_eq!(locus.single_motif().expect("uno").as_text(), "CAG");
    assert_eq!(locus.alleles()[0].tract().blocks().len(), 1);
}

// ---------------------------------------------------------------------------
// Decisión 2 del responsable: los mínimos también por locus.
// ---------------------------------------------------------------------------

/// **Dos loci de la misma ejecución que necesitan umbrales distintos.**
///
/// El general exige pureza 0,9. El primer locus la cumple de sobra; el segundo
/// tiene pureza 0,75 y **solo entra porque declara el suyo**. Si los dos
/// pasaran con el umbral general la prueba no discriminaría nada, y por eso la
/// segunda mitad quita el umbral propio y exige que entonces sí se rechace.
#[test]
fn two_loci_in_the_same_run_need_different_minima() {
    let document = |with_own_minimum: bool| {
        format!(
            r#"{}
[reference.repeat_admission]
min_purity = 0.9

[[reference.repeat_loci]]
id = "STR-PURE"
sequence = "chr1"
start = 101
end = 140
motif = "GATA"

[[reference.repeat_loci.alleles]]
units = 10
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-IMPURE"
sequence = "chr1"
start = 201
end = 216
motif = "CA"
{}

[[reference.repeat_loci.alleles]]
sequence = "CACACACACACATTTT"
frequency = 1.0
"#,
            header("4"),
            if with_own_minimum {
                "min_purity = 0.6"
            } else {
                ""
            }
        )
    };

    let effective = load(&document(true)).expect("cada locus con su umbral");
    assert_eq!(effective.repeat_loci().len(), 2);
    assert!((effective.repeat_loci()[0].admission().min_purity() - 0.9).abs() < f64::EPSILON);
    assert!((effective.repeat_loci()[1].admission().min_purity() - 0.6).abs() < f64::EPSILON);
    assert!(
        (effective.repeat_admission().min_purity() - 0.9).abs() < f64::EPSILON,
        "el general sigue publicándose como valor de partida"
    );
    let impure = effective.repeat_loci()[1].alleles()[0].tract().annotation();
    assert!((impure.purity_as_f64() - 0.75).abs() < 1e-12);

    let error = load(&document(false)).expect_err("sin umbral propio, el general lo rechaza");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-IMPURE"), "{rendered}");
    assert!(rendered.contains("pureza 0.75"), "{rendered}");
    assert!(
        !rendered.contains("STR-PURE"),
        "el primer locus sí cumple el general: {rendered}"
    );
}

/// El umbral propio también puede ser **más estricto** que el general, que es
/// el caso de los motivos cortos: tres unidades de un mononucleótido son tres
/// bases, que salen por azar a cada paso.
#[test]
fn a_locus_can_declare_a_stricter_minimum_than_the_general_one() {
    let mononucleotide = |units: u64, minimum: u64| {
        format!(
            r#"{}
[[reference.repeat_loci]]
id = "POLY-A"
sequence = "chr1"
start = 101
end = 120
motif = "A"
min_longest_run_units = {minimum}

[[reference.repeat_loci.alleles]]
units = {units}
frequency = 1.0
"#,
            header("4")
        )
    };
    // Con el general —tres unidades— un poli-A de cinco entraría. Con doce, que
    // es lo que un mononucleótido necesita para no ser ruido, se rechaza.
    load(&mononucleotide(5, 3)).expect("cinco unidades con el mínimo general");
    let error = load(&mononucleotide(5, 12)).expect_err("cinco unidades con el mínimo propio");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus POLY-A"), "{rendered}");
    assert!(rendered.contains("tirada más larga de 5"), "{rendered}");
    assert!(rendered.contains("mínimo declarado es 12"), "{rendered}");
    load(&mononucleotide(12, 12)).expect("doce unidades con el mínimo propio");
}

/// Los mínimos de cada locus se publican en la configuración efectiva, y son
/// los suyos, no los generales.
#[test]
fn each_locus_publishes_its_own_minima() {
    let body = composite_locus("min_purity = 0.5\nmin_longest_run_units = 4");
    let document = published_effective_config(&body);
    let admission = &document["repeat_loci"][0]["admission"];
    assert_eq!(admission["min_purity"], 0.5);
    assert_eq!(admission["min_longest_run_units"], 4);
    assert_eq!(
        document["repeat_admission"]["min_purity"], 0.8,
        "el general sigue siendo el de omisión y se publica aparte"
    );
}

/// **El tramo decisivo de la cascada: el general pone la pureza y el locus solo
/// la tirada.**
///
/// *Añadido tras la segunda revisión adversarial.* Estaban cubiertos el locus
/// que no declara nada y el que declara el mismo campo que el general; faltaba
/// el **mixto**, que es donde la herencia se ejercita de verdad. Sustituyendo
/// los dos valores de partida por los de omisión, las 54 pruebas de los dos
/// ficheros pasaban.
///
/// Y el cambio va en la dirección peligrosa: con el general en 0,9 y un locus
/// que solo declara la tirada, el comportamiento correcto rechaza con 0,9 y el
/// defectuoso con 0,8. **La guarda se relajaría por debajo de lo que el usuario
/// escribió.**
#[test]
fn a_locus_that_declares_only_the_run_still_inherits_the_general_purity() {
    let body = format!(
        r#"{}
[reference.repeat_admission]
min_purity = 0.9

[[reference.repeat_loci]]
id = "STR-MIX"
sequence = "chr1"
start = 101
end = 116
motif = "CA"
min_longest_run_units = 5

[[reference.repeat_loci.alleles]]
sequence = "CACACACACACATTTT"
frequency = 1.0
"#,
        header("4")
    );
    // Pureza 0,75: por debajo del 0,9 general, que el locus NO declara y por
    // tanto hereda. El mensaje tiene que decir 0.9, no 0.8.
    let error = load(&body).expect_err("la pureza heredada del general lo rechaza");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("locus STR-MIX"), "{rendered}");
    assert!(rendered.contains("pureza 0.75"), "{rendered}");
    assert!(
        rendered.contains("mínimo declarado es 0.9"),
        "el locus hereda la pureza del general, no la de omisión: {rendered}"
    );

    // Con un alelo que sí llega al 0,9 heredado, el locus entra y se ve que la
    // tirada declarada por él sí manda: la suya es 5 y la general sería 3.
    let admitted = body.replace("sequence = \"CACACACACACATTTT\"", "units = 8");
    let effective = load(&admitted).expect("alelo puro de ocho unidades");
    let admission = effective.repeat_loci()[0].admission();
    assert!(
        (admission.min_purity() - 0.9).abs() < f64::EPSILON,
        "la pureza se hereda del general"
    );
    assert_eq!(
        admission.min_longest_run_units(),
        5,
        "la tirada es la que el locus declara"
    );

    // Y baja lo suficiente, la tirada propia del locus rechaza donde la general
    // habría admitido: cuatro unidades pasan el 3 general y no el 5 del locus.
    let short = admitted.replace("units = 8", "units = 4");
    let error = load(&short).expect_err("cuatro unidades con la tirada propia en cinco");
    assert!(
        format!("{error:?}").contains("mínimo declarado es 5"),
        "{error:?}"
    );
}

/// Un mínimo inválido se denuncia **con el campo que el usuario escribió** y
/// nombrando el locus.
///
/// *Añadido tras la segunda revisión adversarial.* `Diagnostic::with` sobrescribe
/// la clave, de modo que el envoltorio que añadía `fields` pisaba el campo que el
/// error interior había nombrado: con la tirada a cero, el diagnóstico mandaba a
/// mirar `min_purity`, una clave que el usuario no ha escrito.
#[test]
fn an_invalid_minimum_names_the_field_that_was_written_and_the_locus() {
    let locus_level = (header("4") + MINIMAL_LOCUS).replace(
        "motif = \"CAG\"",
        "motif = \"CAG\"\nmin_longest_run_units = 0",
    );
    let error = load(&locus_level).expect_err("tirada mínima cero en el locus");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("reference.repeat_loci.min_longest_run_units"),
        "el campo denunciado es el que se escribió: {error:?}"
    );
    assert!(
        format!("{error:?}").contains("locus STR-1"),
        "el diagnóstico nombra el locus, como todos los demás de repeticiones: {error:?}"
    );

    let locus_purity = (header("4") + MINIMAL_LOCUS)
        .replace("motif = \"CAG\"", "motif = \"CAG\"\nmin_purity = 1.5");
    let error = load(&locus_purity).expect_err("pureza mínima fuera de rango");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("reference.repeat_loci.min_purity")
    );
    assert!(format!("{error:?}").contains("locus STR-1"), "{error:?}");

    // Y el mismo error en los mínimos generales denuncia su propio campo y no
    // nombra ningún locus, porque no pertenece a ninguno.
    let general = (header("4") + MINIMAL_LOCUS).replace(
        "[reference]",
        "[reference.repeat_admission]\nmin_longest_run_units = 0\n\n[reference]",
    );
    let error = load(&general).expect_err("tirada mínima cero en los generales");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("reference.repeat_admission.min_longest_run_units")
    );
    assert!(!format!("{error:?}").contains("locus "), "{error:?}");
}

/// Las tres divergencias que la decisión de varios motivos abrió entre el
/// esquema de autoría y el lector.
///
/// *Añadido tras la segunda revisión adversarial.* Dos se cierran en el esquema;
/// la tercera —que la lista de unidades tenga la misma longitud que la de
/// motivos— **no es expresable** en el formato, y por eso el esquema lo dice en
/// su descripción y la fija esta prueba sobre el lector.
#[test]
fn the_authoring_schema_closes_what_it_can_of_the_multi_motif_frontier() {
    let with_motifs = |allele: serde_json::Value| {
        serde_json::json!({
            "config_kind": "general",
            "schema_version": "4",
            "reference": {
                "kind": "simulated",
                "id": "ref-main",
                "composition": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                "sequences": [{"id": "chr1", "class": "chromosome", "length": 1000}],
                "repeat_loci": [{
                    "id": "STR-1", "sequence": "chr1", "start": 101, "end": 160,
                    "motifs": ["CAG", "CA"], "alleles": [allele]
                }]
            }
        })
    };
    let schema = "config-general-v4.schema.json";

    // Control positivo.
    assert!(validates(
        &with_motifs(serde_json::json!({"units": [4, 5], "frequency": 1.0})),
        schema
    ));

    assert!(
        !validates(
            &with_motifs(serde_json::json!({"sequence": "CAGCAGCAG", "frequency": 1.0})),
            schema
        ),
        "con varios motivos la secuencia exacta no vale, y el esquema ya lo cierra"
    );
    assert!(
        !validates(
            &with_motifs(serde_json::json!({"units": 4, "frequency": 1.0})),
            schema
        ),
        "con varios motivos el escalar no vale, y el esquema ya lo cierra"
    );

    // La tercera: el esquema **admite** una lista de longitud distinta, porque
    // no puede compararla con la de motivos. Lo dice en su descripción, y quien
    // la rechaza es el lector, nombrando el locus.
    let mismatched = with_motifs(serde_json::json!({"units": [4, 5, 6], "frequency": 1.0}));
    assert!(
        validates(&mismatched, schema),
        "si algún día el formato permitiera expresarlo, este aviso sobra y hay que cerrarlo"
    );
    let body = format!(
        r#"{}
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 160
motifs = ["CAG", "CA"]

[[reference.repeat_loci.alleles]]
units = [4, 5, 6]
frequency = 1.0
"#,
        header("4")
    );
    let error = load(&body).expect_err("tres valores para dos bloques");
    assert!(
        format!("{error:?}").contains("hace falta uno por bloque"),
        "{error:?}"
    );
}

/// **Los esquemas heredados no cambian ni un byte.** La regla de la sección
/// 10.3 del plan: los esquemas publicados no se editan; la capacidad nueva va
/// en una versión nueva.
#[test]
fn the_inherited_schemas_are_byte_for_byte_unchanged() {
    let expected = [
        (
            "config-general-v1.schema.json",
            "dbe5d16813612bdede5fb65a97d8646129322b174cf1a637841312657cd81b57",
        ),
        (
            "config-general-v2.schema.json",
            "29bd4aeec6663f7a7f08095a075e4a71ffca216a803d8aa4f35c0d5130c84740",
        ),
        (
            "effective-config-v4.schema.json",
            "1647da1cda5c512182dd7d9565c36d7dd6497e208c08acf6760efc0a2c3b5d3a",
        ),
        (
            "effective-config-v5.schema.json",
            "71f3758fbae8603103aecb7c1851812e6571d7627ce07c50f56748db88e8d9e5",
        ),
        (
            "truth-v1.schema.json",
            "e3d49abd11447370d34e4acf29b755a288fb310641474218086ed60099233fd7",
        ),
        (
            "truth-v2.schema.json",
            "be53a6b32eb04689bcb4113731bb97f84e57f3566116e98fa2a2500313ba647b",
        ),
    ];
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for (name, digest) in expected {
        let identity = genome_population_simulator::integrity::identify_file(
            &root.join("schemas").join(name),
            name,
        )
        .expect("huella del esquema");
        assert_eq!(identity.digest(), digest, "el esquema {name} ha cambiado");
    }
}

/// Publica de verdad la ejecución y devuelve el `effective-config.json` que
/// queda escrito, en vez de reconstruir el documento a mano: el oráculo es el
/// artefacto publicado, no otra llamada a la misma función que lo produce.
fn published_effective_config(body: &str) -> serde_json::Value {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio de trabajo");
    let config = root.join("general.toml");
    fs::write(&config, body).expect("documento general");
    let effective = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective, &destination).expect("publicación");
    let document: serde_json::Value = serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("configuración efectiva"),
    )
    .expect("JSON");
    fs::remove_dir_all(&root).ok();
    document
}

/// Ejecuta el validador externo y devuelve si admitió el documento.
fn validates(document: &serde_json::Value, schema: &str) -> bool {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let path = temporary("json");
    fs::write(
        &path,
        serde_json::to_vec_pretty(document).expect("serialización"),
    )
    .expect("documento");
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(&path)
        .output()
        .expect("validador externo de JSON Schema");
    fs::remove_file(&path).ok();
    output.status.success()
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

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f6-1-1-1-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}
