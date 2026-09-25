//! F4.9.7 a F4.9.11: el recorrido completo con preparación amplificada.
//!
//! # Lo que estas pruebas miran de verdad
//!
//! Casi todas cuentan cosas **sobre los ficheros publicados**, no sobre el
//! estado interno del productor: los duplicados se cuentan releyendo los FASTQ,
//! el sesgo por composición se mide clasificando las lecturas escritas, y la
//! reducción al caso sin amplificar se comprueba **byte a byte**.
//!
//! # Cuidado con los errores de secuenciación
//!
//! Con los errores activados los duplicados dejan de ser idénticos y el oráculo
//! que los cuenta desde el FASTQ no discrimina nada. Todos los fixtures que
//! cuentan duplicados los **desactivan**, y aquí queda dicho.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use noodles_bgzf as bgzf;
use serde_json::Value;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-9-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Qué distingue una ejecución de otra. Todo lo demás es idéntico entre ellas,
/// que es lo que permite atribuir una diferencia al parámetro que cambió.
struct Scenario {
    bases: String,
    coverage: f64,
    /// Bloque completo del evento técnico que va tras la cobertura.
    event_extra: String,
    /// Individuos de la tanda. El padrón siempre declara los dos, de modo que
    /// cambiar esta lista **no cambia ningún genotipo**.
    selected: Vec<&'static str>,
    /// Modo de errores de secuenciación. Los fixtures que cuentan duplicados
    /// los desactivan, porque con ellos los duplicados dejan de ser idénticos.
    sequencing_errors: &'static str,
}

struct Published {
    root: PathBuf,
    destination: PathBuf,
}

impl Drop for Published {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

impl Published {
    fn effective(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("effective-config.json")).expect("configuración"),
        )
        .expect("JSON")
    }

    fn event(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("events/reads-main/manifest.json"))
                .expect("manifiesto de evento"),
        )
        .expect("JSON")
    }

    fn fastq_bytes(&self, individual: &str, mate: &str) -> Vec<u8> {
        fs::read(self.destination.join(format!(
            "outputs/fastq/reads-main/{individual}_{mate}.fastq.gz"
        )))
        .expect("FASTQ")
    }

    /// Las parejas publicadas de un individuo, en el orden del fichero.
    fn pairs(&self, individual: &str) -> Vec<(Vec<u8>, Vec<u8>)> {
        let first = read_sequences(
            &self
                .destination
                .join(format!("outputs/fastq/reads-main/{individual}_R1.fastq.gz")),
        );
        let second = read_sequences(
            &self
                .destination
                .join(format!("outputs/fastq/reads-main/{individual}_R2.fastq.gz")),
        );
        assert_eq!(first.len(), second.len());
        first.into_iter().zip(second).collect()
    }
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

fn run(scenario: &Scenario) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let length = scenario.bases.len();
    fs::write(
        root.join("reference.fa"),
        format!(">chr1\n{}\n", scenario.bases),
    )
    .expect("FASTA");
    // La base de referencia del sitio SNV se **lee de la secuencia**, no se
    // escribe a mano: un fixture con otra composición cambiaría esa base y el
    // lector rechazaría la matriz por una discrepancia que no tiene nada que
    // ver con lo que la prueba mira.
    let reference_base = scenario.bases.as_bytes()[19] as char;
    let alternate = if reference_base == 'T' { 'G' } else { 'T' };
    fs::write(
        root.join("population-snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\n\
             chr1\t20\t{reference_base}\t{alternate}\tP1\t0.5\n\
             chr1\t20\t{reference_base}\t{alternate}\tP2\t0.5\n"
        ),
    )
    .expect("matriz");
    let selected = scenario
        .selected
        .iter()
        .map(|id| format!("\"{id}\""))
        .collect::<Vec<_>>()
        .join(", ");
    // **El UDI va atado a la individua, no a su puesto en la lista.** Si se
    // repartieran por orden, quitar una biblioteca de la tanda cambiaría el
    // UDI de la otra y con él la cabecera de todas sus lecturas: la prueba de
    // independencia entre bibliotecas fallaría por culpa del fixture y no del
    // producto. Ocurrió, y por eso queda escrito.
    let libraries = scenario
        .selected
        .iter()
        .map(|id| {
            let (i7, i5) = match *id {
                "sample-A" => ("ACGTACGT", "TGCATGCA"),
                "sample-B" => ("TTTTAAAA", "CCCCGGGG"),
                other => panic!("individua sin UDI asignado: {other}"),
            };
            format!(
                "\n[[technical_event.libraries]]\nindividual_id = \"{id}\"\ni7 = \"{i7}\"\ni5 = \"{i5}\"\n"
            )
        })
        .collect::<String>();
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
length = {length}

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[[initial_populations]]
population_id = "P2"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-B"]

[technical_event]
event_id = "reads-main"
individual_ids = [{selected}]
target_coverage_wgs = {coverage}
sequencing_errors = "{errors}"
{event_extra}{libraries}"#,
            coverage = scenario.coverage,
            errors = scenario.sequencing_errors,
            event_extra = scenario.event_extra,
        ),
    )
    .expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Published { root, destination }
}

/// Bases de composición equilibrada y **sin periodo corto**.
///
/// El generador congruencial no es un adorno: con un patrón periódico, dos
/// posiciones distintas del genoma producen la misma lectura y el oráculo que
/// cuenta duplicados desde el FASTQ empieza a contar coincidencias que no son
/// duplicados de amplificación. Se midió: un patrón de periodo cuatro daba 194
/// coincidencias donde el productor había hecho 121 familias repetidas.
fn uniform_bases(length: usize) -> String {
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

/// Una mitad pobre en GC (30 %) y otra rica (70 %).
fn two_halves(half: usize) -> String {
    let poor: String = (0..half)
        .map(|index| if index % 10 < 3 { 'G' } else { 'A' })
        .collect();
    let rich: String = (0..half)
        .map(|index| if index % 10 < 7 { 'C' } else { 'T' })
        .collect();
    format!("{poor}{rich}")
}

fn amplified(extra: &str) -> String {
    format!("library_preparation = \"pcr\"\n\n[technical_event.pcr]\n{extra}")
}

fn pcr_section(event: &Value) -> &Value {
    &event["artifacts"]["pcr"]
}

// ---------------------------------------------------------------------------
// F4.9.11 — reducción al caso sin amplificar.
// ---------------------------------------------------------------------------

/// **Con fracción de duplicados cero, curva neutra y tasas a cero, los FASTQ
/// son byte a byte iguales a los de la preparación sin PCR.**
///
/// Es el contrato que protege todo lo publicado antes de esta entrega. Falla
/// contra una implementación que consuma aleatoriedad con parámetros neutros,
/// contra una que desplace los ordinales de plantilla y contra una inversión
/// numérica que devuelva un `m` distinto de cero para fracción cero.
#[test]
fn the_neutral_amplified_run_writes_the_same_bytes_as_the_unamplified_one() {
    let bases = uniform_bases(4_000);
    let plain = run(&Scenario {
        bases: bases.clone(),
        coverage: 8.0,
        event_extra: String::new(),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let neutral = run(&Scenario {
        bases,
        coverage: 8.0,
        event_extra: amplified("duplicate_fraction = 0\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    for mate in ["R1", "R2"] {
        assert_eq!(
            plain.fastq_bytes("sample-A", mate),
            neutral.fastq_bytes("sample-A", mate),
            "el {mate} amplificado con parámetros neutros tiene que salir byte a byte igual"
        );
    }
    // Control positivo del mismo par de ficheros: una fracción de duplicados
    // distinta de cero **sí** cambia los bytes, de modo que la igualdad de
    // arriba no pasa por comparar dos ficheros que nunca podrían diferir.
    let busy = run(&Scenario {
        bases: uniform_bases(4_000),
        coverage: 8.0,
        event_extra: amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    assert_ne!(
        plain.fastq_bytes("sample-A", "R1"),
        busy.fastq_bytes("sample-A", "R1"),
        "con duplicados de verdad los bytes tienen que cambiar"
    );

    // **Y con los errores de secuenciación encendidos también.** El plan pide
    // el fixture con las tasas a cero, pero los errores son el mayor consumidor
    // de aleatoriedad aguas abajo: si la amplificación neutra desplazara una
    // sola palabra, aquí se vería aunque en el fixture neutro no.
    let noisy_plain = run(&Scenario {
        bases: uniform_bases(4_000),
        coverage: 8.0,
        event_extra: String::new(),
        selected: vec!["sample-A"],
        sequencing_errors: "enabled",
    });
    let noisy_neutral = run(&Scenario {
        bases: uniform_bases(4_000),
        coverage: 8.0,
        event_extra: amplified("duplicate_fraction = 0\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "enabled",
    });
    for mate in ["R1", "R2"] {
        assert_eq!(
            noisy_plain.fastq_bytes("sample-A", mate),
            noisy_neutral.fastq_bytes("sample-A", mate),
            "con errores activados el {mate} tampoco puede cambiar"
        );
    }
    // Control positivo del par ruidoso: los errores se aplicaron de verdad, de
    // modo que la igualdad no sale de comparar dos ficheros sin errores.
    assert_ne!(
        noisy_plain.fastq_bytes("sample-A", "R1"),
        plain.fastq_bytes("sample-A", "R1"),
        "con los errores activados las lecturas tienen que cambiar"
    );
}

/// **El manifiesto de la ejecución neutra difiere solo en campos declarativos,
/// y sus recuentos confirman que la etapa de amplificación no consumió nada.**
///
/// El recuento de intentos igual al de familias demuestra que **ninguna clave
/// de reintento se abrió**; el histograma con una sola celda en el tamaño 1
/// demuestra que ningún tamaño de familia se sorteó; y la `m` exactamente cero
/// demuestra que la inversión salió por su retorno anticipado.
#[test]
fn the_neutral_amplified_manifest_shows_no_amplification_decisions() {
    let plain = run(&Scenario {
        bases: uniform_bases(4_000),
        coverage: 8.0,
        event_extra: String::new(),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let neutral = run(&Scenario {
        bases: uniform_bases(4_000),
        coverage: 8.0,
        event_extra: amplified("duplicate_fraction = 0\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let pcr = pcr_section(&neutral.event()).clone();
    let families = pcr["families"].as_u64().expect("familias");
    let pairs = pcr["pairs"].as_u64().expect("parejas");
    assert!(pairs > 0, "control positivo: tiene que haber parejas");
    assert_eq!(families, pairs, "cada molécula da exactamente una pareja");
    assert_eq!(
        pcr["gc_attempts_total"].as_u64(),
        Some(families),
        "con curva neutra no se abre ni un reintento"
    );
    assert_eq!(pcr["gc_attempts_exhausted_families"].as_u64(), Some(0));
    assert_eq!(pcr["duplicate_pairs"].as_u64(), Some(0));
    assert_eq!(pcr["realized_duplicate_fraction"].as_f64(), Some(0.0));
    assert_eq!(
        pcr["family_size_histogram"],
        serde_json::json!([{"family_size": 1, "families": families}])
    );
    assert_eq!(
        pcr["mean_pairs_per_molecule"].as_f64(),
        Some(0.0),
        "la inversión de la fracción cero tiene que dar cero exacto"
    );
    assert_eq!(
        pcr["duplicate_fraction_inversion"]["iterations"].as_u64(),
        Some(0),
        "la fracción cero sale por retorno anticipado, sin iterar"
    );

    // Los dos manifiestos coinciden en todo lo que no es declarativo.
    let mut plain_event = plain.event();
    let mut neutral_event = neutral.event();
    for document in [&mut plain_event, &mut neutral_event] {
        let object = document.as_object_mut().expect("objeto");
        object.remove("schema_version");
        object.remove("rng_catalog");
        object.remove("profile");
        let artifacts = object
            .get_mut("artifacts")
            .and_then(Value::as_object_mut)
            .expect("artefactos");
        artifacts.remove("pcr");
    }
    assert_eq!(
        plain_event, neutral_event,
        "los dos manifiestos solo pueden diferir en los campos declarativos"
    );
}

// ---------------------------------------------------------------------------
// F4.9.8 — familias, duplicados e independencia entre bibliotecas.
// ---------------------------------------------------------------------------

/// Cuenta, **leyendo solo el FASTQ**, cuántas parejas repiten la secuencia de
/// otra.
fn duplicate_pairs_from_reads(pairs: &[(Vec<u8>, Vec<u8>)]) -> u64 {
    let distinct = pairs
        .iter()
        .cloned()
        .collect::<std::collections::BTreeSet<_>>();
    (pairs.len() - distinct.len()) as u64
}

/// **Los duplicados que el manifiesto publica son los que hay escritos.**
///
/// *Cuidado:* con los errores de secuenciación activados los duplicados dejan
/// de ser idénticos y este oráculo no discriminaría nada. El fixture los
/// desactiva, la quimera está a cero y no hay loci declarados.
///
/// Falla contra un manifiesto que publique el recuento **esperado** en vez del
/// realizado, y contra una implementación que sortee el tamaño de familia pero
/// emita una sola pareja: la exigencia de que el recuento sea mayor que cero
/// evita que la igualdad se cumpla por el motivo equivocado.
#[test]
fn the_published_duplicate_count_is_the_one_written_in_the_reads() {
    let published = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 6.0,
        event_extra: amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let event = published.event();
    let pcr = pcr_section(&event);
    let declared = pcr["duplicate_pairs"].as_u64().expect("duplicados");
    let observed = duplicate_pairs_from_reads(&published.pairs("sample-A"));
    assert!(
        declared > 0,
        "con una fracción de duplicados alta tiene que haber duplicados"
    );
    assert_eq!(
        observed, declared,
        "el recuento del manifiesto tiene que ser el realizado, no el esperado"
    );
    // Y la fracción realizada **no** es la declarada: la desigualdad de los
    // primeros ciclos la sube.
    let realized = pcr["realized_duplicate_fraction"]
        .as_f64()
        .expect("fracción");
    let declared_fraction = pcr["declared_duplicate_fraction"]
        .as_f64()
        .expect("declarada");
    assert!(
        realized > declared_fraction,
        "la realizada ({realized}) tiene que superar a la declarada ({declared_fraction})"
    );
}

/// **Las parejas de una familia no salen adyacentes en el fichero publicado.**
///
/// La escritura ordena por identificador de instancia, no por orden de
/// generación, exactamente como en un FASTQ real.
#[test]
fn the_pairs_of_a_family_are_not_adjacent_in_the_published_file() {
    let published = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 6.0,
        event_extra: amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let pairs = published.pairs("sample-A");
    let mut positions: std::collections::BTreeMap<&(Vec<u8>, Vec<u8>), Vec<usize>> =
        std::collections::BTreeMap::new();
    for (index, pair) in pairs.iter().enumerate() {
        positions.entry(pair).or_default().push(index);
    }
    let families: Vec<&Vec<usize>> = positions
        .values()
        .filter(|indices| indices.len() > 1)
        .collect();
    // Control positivo: tiene que haber familias con más de una pareja, o la
    // comprobación de abajo no estaría mirando nada.
    assert!(
        !families.is_empty(),
        "el fixture tiene que producir familias de más de una pareja"
    );
    let separated = families
        .iter()
        .any(|indices| indices.windows(2).any(|window| window[1] != window[0] + 1));
    assert!(
        separated,
        "alguna familia tiene que quedar separada en el fichero: {families:?}"
    );
}

/// **Añadir una biblioteca no cambia las familias de las demás.**
///
/// Falla contra un ámbito aleatorio indexado por **orden de biblioteca**: con
/// él, `sample-B` sería la primera en una tanda y la segunda en la otra, y sus
/// familias cambiarían. El padrón declara las dos individuas en los dos casos,
/// de modo que sus genotipos son idénticos y la única diferencia es la tanda.
#[test]
fn adding_a_library_does_not_change_the_families_of_the_others() {
    let extra = amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n");
    let alone = run(&Scenario {
        bases: uniform_bases(8_000),
        coverage: 6.0,
        event_extra: extra.clone(),
        selected: vec!["sample-B"],
        sequencing_errors: "disabled",
    });
    let together = run(&Scenario {
        bases: uniform_bases(8_000),
        coverage: 6.0,
        event_extra: extra,
        selected: vec!["sample-A", "sample-B"],
        sequencing_errors: "disabled",
    });
    for mate in ["R1", "R2"] {
        assert_eq!(
            alone.fastq_bytes("sample-B", mate),
            together.fastq_bytes("sample-B", mate),
            "el {mate} de sample-B no puede depender de quién más esté en la tanda"
        );
    }
    // Control positivo: la tanda con dos bibliotecas publica de verdad la
    // segunda, de modo que la igualdad no pasa por no haber añadido nada.
    assert!(
        together
            .destination
            .join("outputs/fastq/reads-main/sample-A_R1.fastq.gz")
            .exists(),
        "la tanda conjunta tiene que publicar también sample-A"
    );
    let event = together.event();
    assert!(
        pcr_section(&event)["families"].as_u64().expect("familias")
            > pcr_section(&alone.event())["families"]
                .as_u64()
                .expect("familias"),
        "la tanda conjunta tiene más familias que la de una sola biblioteca"
    );
}

// ---------------------------------------------------------------------------
// F4.9.7 — la composición actúa al elegir la molécula.
// ---------------------------------------------------------------------------

/// Proporción de lecturas que vienen de la mitad rica en GC.
fn rich_half_share(pairs: &[(Vec<u8>, Vec<u8>)]) -> f64 {
    let rich = pairs
        .iter()
        .filter(|(first, _)| {
            let gc = first
                .iter()
                .filter(|base| matches!(base, b'G' | b'C'))
                .count();
            (gc as f64) / (first.len() as f64) > 0.5
        })
        .count();
    (rich as f64) / (pairs.len() as f64)
}

/// **Una curva que penaliza la mitad rica en GC la infrarrepresenta.**
///
/// Falla contra una implementación que aplique la composición solo al tamaño
/// de familia: con una fracción de duplicados pequeña casi todas las familias
/// tienen tamaño 1 en las dos mitades, de modo que esa variante dejaría la
/// proporción en aproximadamente 0,5. El fixture exige la fracción **pequeña**
/// a propósito: con una fracción grande la variante defectuosa también produce
/// sesgo y las dos dejan de distinguirse.
#[test]
fn the_composition_curve_underrepresents_the_penalised_half() {
    let curve = "\n[[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 0.0\n\
                 relative_efficiency = 1.0\n\n\
                 [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 0.5\n\
                 relative_efficiency = 1.0\n\n\
                 [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 0.55\n\
                 relative_efficiency = 0.2\n\n\
                 [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 1.0\n\
                 relative_efficiency = 0.2\n";
    let biased = run(&Scenario {
        bases: two_halves(5_000),
        coverage: 20.0,
        event_extra: amplified(&format!(
            "duplicate_fraction = 0.01\ncycles = 10\nefficiency = 0.9\n{curve}"
        )),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let share = rich_half_share(&biased.pairs("sample-A"));
    assert!(
        share < 0.15,
        "la mitad penalizada tiene que caer muy por debajo de 0,5, y vale {share}"
    );
    // Control positivo: con la curva neutra y todo lo demás igual, las dos
    // mitades se reparten las lecturas casi a partes iguales. Sin él, una
    // implementación que perdiera la mitad rica por cualquier otro motivo
    // pasaría la comprobación de arriba.
    let neutral = run(&Scenario {
        bases: two_halves(5_000),
        coverage: 20.0,
        event_extra: amplified("duplicate_fraction = 0.01\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let neutral_share = rich_half_share(&neutral.pairs("sample-A"));
    assert!(
        (0.4..=0.6).contains(&neutral_share),
        "con curva neutra las dos mitades se reparten las lecturas, y vale {neutral_share}"
    );
    // Y el modo publicado distingue las dos ejecuciones.
    assert_eq!(biased.event()["artifacts"]["gc"]["mode"], "curve");
    assert_eq!(neutral.event()["artifacts"]["gc"]["mode"], "neutral");
}

/// **Con un tope de dos intentos hay familias aceptadas al agotarlo, y con
/// dieciséis hay menos.**
///
/// El discriminador es el **recuento de aceptadas al agotar el tope**: un
/// rechazo sin tope lo dejaría necesariamente en cero. La primera redacción de
/// esta prueba decía que la variante sin tope «no termina», y no es cierto: con
/// una curva que rechaza casi siempre el bucle sí termina, solo que tarde, de
/// modo que la prueba quedaba dependiendo del reloj y era intermitente.
#[test]
fn exhausting_the_attempt_cap_is_counted_and_shrinks_with_a_larger_cap() {
    let curve = "\n[[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 0.0\n\
                 relative_efficiency = 0.76\n\n\
                 [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 1.0\n\
                 relative_efficiency = 0.76\n";
    let counts: Vec<(u64, u64)> = [2_u32, 16]
        .into_iter()
        .map(|attempts| {
            let published = run(&Scenario {
                bases: uniform_bases(8_000),
                coverage: 6.0,
                event_extra: amplified(&format!(
                    "duplicate_fraction = 0.05\ncycles = 10\nefficiency = 0.9\n\
                     gc_max_attempts = {attempts}\n{curve}"
                )),
                selected: vec!["sample-A"],
                sequencing_errors: "disabled",
            });
            let event = published.event();
            let pcr = pcr_section(&event);
            (
                pcr["gc_attempts_exhausted_families"]
                    .as_u64()
                    .expect("agotadas"),
                pcr["families"].as_u64().expect("familias"),
            )
        })
        .collect();
    let (exhausted_two, families_two) = counts[0];
    let (exhausted_sixteen, families_sixteen) = counts[1];
    assert!(
        exhausted_two > 0,
        "con tope de dos intentos tiene que haber familias aceptadas al agotarlo"
    );
    assert!(
        exhausted_sixteen < exhausted_two,
        "con dieciséis intentos tiene que haber menos ({exhausted_sixteen}) que con dos \
         ({exhausted_two})"
    );
    // Control positivo: las dos ejecuciones produjeron familias de verdad, de
    // modo que la desigualdad no sale de que una no llegara a ejecutarse.
    assert!(families_two > 0 && families_sixteen > 0);
}

// ---------------------------------------------------------------------------
// F4.9.9 — recuentos agregados del manifiesto.
// ---------------------------------------------------------------------------

/// **Las identidades internas del histograma y de la fracción realizada.**
#[test]
fn the_aggregate_counts_are_internally_consistent() {
    let published = run(&Scenario {
        bases: uniform_bases(8_000),
        coverage: 6.0,
        event_extra: amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let event = published.event();
    let pcr = pcr_section(&event);
    let families = pcr["families"].as_u64().expect("familias");
    let pairs = pcr["pairs"].as_u64().expect("parejas");
    let histogram = pcr["family_size_histogram"].as_array().expect("histograma");
    // **Ponderada por su índice**, que es lo que impide que una celda cuente
    // familias en vez de parejas.
    let weighted: u64 = histogram
        .iter()
        .map(|entry| {
            entry["family_size"].as_u64().expect("tamaño")
                * entry["families"].as_u64().expect("familias")
        })
        .sum();
    let total_families: u64 = histogram
        .iter()
        .map(|entry| entry["families"].as_u64().expect("familias"))
        .sum();
    assert_eq!(
        weighted, pairs,
        "la suma ponderada tiene que dar las parejas"
    );
    assert_eq!(total_families, families);
    assert_eq!(
        pcr["duplicate_pairs"].as_u64(),
        Some(pairs - families),
        "los duplicados son las parejas menos las familias"
    );
    assert_eq!(
        pcr["realized_duplicate_fraction"].as_f64(),
        Some(1.0 - (families as f64) / (pairs as f64))
    );
    // Control positivo: el histograma tiene que contener más de un tamaño, o
    // las identidades de arriba se cumplirían con una sola celda trivial.
    assert!(
        histogram.len() > 1,
        "el fixture tiene que producir familias de varios tamaños: {histogram:?}"
    );
    // **La cobertura de moléculas distintas, contada desde los bytes.**
    //
    // Rehacer aquí la fórmula con los mismos números que el manifiesto publica
    // sería un espejo: los dos lados se moverían juntos. El conteo
    // independiente son los **grupos de parejas idénticas** del FASTQ, que es
    // lo que un analista tiene tras eliminar duplicados.
    let basic: f64 = event["coverage"]
        .as_array()
        .expect("cobertura")
        .iter()
        .map(|entry| entry["basic_x"].as_f64().expect("básica"))
        .sum();
    let written = published.pairs("sample-A");
    let distinct = written
        .iter()
        .cloned()
        .collect::<std::collections::BTreeSet<_>>()
        .len() as u64;
    assert_eq!(
        distinct, families,
        "los grupos distintos escritos tienen que ser las familias publicadas"
    );
    assert_eq!(
        pcr["unique_molecule_coverage_x"].as_f64(),
        Some(basic * (distinct as f64) / (written.len() as f64)),
        "la cobertura de moléculas distintas tiene que salir de lo escrito"
    );
    // Control positivo: hubo duplicados de verdad, de modo que la cobertura de
    // moléculas distintas es **menor** que la básica y la igualdad de arriba no
    // se cumple por ser las dos la misma cifra.
    assert!(
        pcr["unique_molecule_coverage_x"].as_f64().expect("única") < basic,
        "con duplicados, la cobertura de moléculas distintas baja"
    );
}

/// **Con dos bibliotecas, la cobertura de moléculas distintas sigue siendo la
/// suma de las de cada una.**
///
/// El manifiesto usa la proporción **agregada**, y la duda razonable es si eso
/// es una media disfrazada. No lo es: todas las bibliotecas de un evento
/// comparten el denominador de la cobertura, de modo que las dos expresiones
/// son idénticas. Aquí se comprueba contra los grupos escritos de **cada**
/// biblioteca por separado.
#[test]
fn the_unique_molecule_coverage_adds_up_across_two_libraries() {
    let published = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 6.0,
        event_extra: amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A", "sample-B"],
        sequencing_errors: "disabled",
    });
    let event = published.event();
    let coverage: std::collections::BTreeMap<String, f64> = event["coverage"]
        .as_array()
        .expect("cobertura")
        .iter()
        .map(|entry| {
            (
                entry["individual_id"]
                    .as_str()
                    .expect("individua")
                    .to_owned(),
                entry["basic_x"].as_f64().expect("básica"),
            )
        })
        .collect();
    assert_eq!(coverage.len(), 2, "control positivo: hay dos bibliotecas");
    let mut expected = 0.0_f64;
    for individual in ["sample-A", "sample-B"] {
        let pairs = published.pairs(individual);
        let distinct = pairs
            .iter()
            .cloned()
            .collect::<std::collections::BTreeSet<_>>()
            .len() as f64;
        assert!(
            distinct < pairs.len() as f64,
            "cada biblioteca tiene que tener duplicados propios"
        );
        expected += coverage[individual] * distinct / (pairs.len() as f64);
    }
    let published_value = pcr_section(&event)["unique_molecule_coverage_x"]
        .as_f64()
        .expect("cobertura de moléculas distintas");
    assert!(
        (published_value - expected).abs() < 1e-12,
        "el agregado publicado ({published_value}) tiene que ser la suma por biblioteca \
         ({expected})"
    );
}

/// **El par de ejecuciones que aísla la desigualdad de los primeros ciclos.**
///
/// Con **cero** ciclos explícitos, curva neutra y sin clases, la fracción
/// realizada coincide con la declarada dentro de la tolerancia muestral. Con
/// **tres** ciclos explícitos y todo lo demás igual, es **estrictamente
/// mayor**.
///
/// *Cuidado, y por eso son dos ejecuciones:* una sola no discrimina, porque
/// declarada y realizada quedan cerca. Y el mecanismo que hay que aislar es el
/// de los primeros ciclos, que actúa **ya por omisión**, no las clases de
/// eficiencia, que son opcionales. Un manifiesto que copiara la declarada en el
/// campo de la realizada pasaría la primera mitad y moriría en la segunda,
/// porque la segunda tendría que subir y no subiría.
#[test]
fn the_early_cycles_inequality_raises_the_realized_fraction_above_the_declared_one() {
    const DECLARED: f64 = 0.2;
    let realized = |early: u32| {
        let published = run(&Scenario {
            bases: uniform_bases(20_000),
            coverage: 8.0,
            event_extra: amplified(&format!(
                "duplicate_fraction = {DECLARED}\ncycles = 10\nefficiency = 0.9\n\
                 early_cycles = {early}\n"
            )),
            selected: vec!["sample-A"],
            sequencing_errors: "disabled",
        });
        let event = published.event();
        let pcr = pcr_section(&event);
        assert_eq!(pcr["declared_duplicate_fraction"].as_f64(), Some(DECLARED));
        assert_eq!(pcr["early_cycles"].as_u64(), Some(u64::from(early)));
        (
            pcr["realized_duplicate_fraction"]
                .as_f64()
                .expect("fracción"),
            pcr["pairs"].as_u64().expect("parejas"),
        )
    };
    let (flat, pairs) = realized(0);
    assert!(pairs > 300, "hacen falta parejas para que la muestra valga");
    assert!(
        (flat - DECLARED).abs() < 0.05,
        "sin ciclos explícitos la realizada ({flat}) tiene que rondar la declarada ({DECLARED})"
    );
    let (uneven, _) = realized(3);
    assert!(
        uneven > flat,
        "con tres ciclos explícitos la realizada ({uneven}) tiene que superar a la de cero \
         ({flat})"
    );
    assert!(
        uneven > DECLARED,
        "y tiene que superar también a la declarada ({DECLARED})"
    );
}

/// El manifiesto y la configuración efectiva pasan su **validación externa**.
#[test]
fn the_published_documents_pass_their_external_schema_validation() {
    let published = run(&Scenario {
        bases: uniform_bases(8_000),
        coverage: 6.0,
        event_extra: amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for (schema, document) in [
        (
            "manifest-event-v3.schema.json",
            published
                .destination
                .join("events/reads-main/manifest.json"),
        ),
        (
            "effective-config-v6.schema.json",
            published.destination.join("effective-config.json"),
        ),
    ] {
        let output = std::process::Command::new("python3")
            .arg(root.join("tests/external_jsonschema_validator.py"))
            .arg(root.join("schemas").join(schema))
            .arg(&document)
            .output()
            .expect("el validador externo se ejecuta");
        assert!(
            output.status.success(),
            "{schema} rechaza el documento publicado: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    // Y el catálogo publicado es el de la etapa que el evento abre.
    assert_eq!(published.event()["rng_catalog"], "f4-technical-v5");
    assert_eq!(published.event()["schema_version"], "3");
    assert_eq!(
        published.effective()["technical_event"]["rng_catalog"],
        "f4-technical-v5"
    );
    assert_eq!(
        published.effective()["technical_event"]["profile"]["preparation"],
        "pcr"
    );
}

// ---------------------------------------------------------------------------
// F4.9.3 — las clases de eficiencia, ejercitadas de extremo a extremo.
// ---------------------------------------------------------------------------

fn two_classes(first_weight: &str, second_weight: &str, efficiency: &str) -> String {
    amplified(&format!(
        "duplicate_fraction = 0.3\ncycles = 10\nefficiency = {efficiency}\n\
         \n[[technical_event.pcr.efficiency_classes]]\nefficiency = 0.95\n\
         weight = {first_weight}\n\
         \n[[technical_event.pcr.efficiency_classes]]\nefficiency = 0.5\n\
         weight = {second_weight}\n"
    ))
}

/// **El sorteo de la clase de eficiencia se ejecuta de verdad y decide algo.**
///
/// # Por qué los pesos y no las clases
///
/// Dos ejecuciones con **las mismas dos clases** y **pesos intercambiados**
/// solo pueden diferir si el sorteo elige de verdad entre ellas. Comparar una
/// clase contra dos no valdría: la diferencia podría venir de cualquier otra
/// parte del documento.
///
/// Falla contra una implementación que anule el resultado del sorteo y tome
/// siempre la primera clase, que es la variante que sobrevivía a toda la
/// batería anterior.
#[test]
fn the_efficiency_class_draw_decides_which_class_a_family_gets() {
    let mostly_first = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 8.0,
        event_extra: two_classes("0.95", "0.05", "0.9"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let mostly_second = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 8.0,
        event_extra: two_classes("0.05", "0.95", "0.9"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    assert_ne!(
        mostly_first.fastq_bytes("sample-A", "R1"),
        mostly_second.fastq_bytes("sample-A", "R1"),
        "intercambiar los pesos de las clases tiene que cambiar las familias"
    );
    // Control positivo: la ejecución es reproducible, de modo que la
    // desigualdad no sale de que dos ejecuciones iguales ya difieran.
    let repeated = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 8.0,
        event_extra: two_classes("0.95", "0.05", "0.9"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    assert_eq!(
        mostly_first.fastq_bytes("sample-A", "R1"),
        repeated.fastq_bytes("sample-A", "R1"),
        "la misma configuración tiene que dar los mismos bytes"
    );
    // Y los dos documentos publican sus dos clases.
    let event = mostly_first.event();
    let classes = pcr_section(&event)["efficiency_classes"]
        .as_array()
        .expect("clases");
    assert_eq!(classes.len(), 2);
    assert_eq!(classes[0]["efficiency"].as_f64(), Some(0.95));
    assert_eq!(classes[1]["efficiency"].as_f64(), Some(0.5));
}

/// **Con clases declaradas, la eficiencia escalar deja de intervenir.**
///
/// Es una consecuencia del contrato que conviene dejar **probada** y no solo
/// escrita: la eficiencia de la tabla es obligatoria, pero cuando se declaran
/// clases el modelo usa la de cada clase y la escalar solo se publica. Un
/// usuario que declare eficiencia 0,9 y clases 0,95 y 0,5 verá publicado un 0,9
/// que no interviene en nada.
///
/// Dos ejecuciones que solo difieren en la escalar producen los **mismos
/// bytes**, y sus documentos sí difieren. Ese par es la prueba.
#[test]
fn the_scalar_efficiency_is_only_published_when_classes_are_declared() {
    let nine = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 8.0,
        event_extra: two_classes("0.7", "0.3", "0.9"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let four = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 8.0,
        event_extra: two_classes("0.7", "0.3", "0.4"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    for mate in ["R1", "R2"] {
        assert_eq!(
            nine.fastq_bytes("sample-A", mate),
            four.fastq_bytes("sample-A", mate),
            "con clases declaradas la eficiencia escalar no puede cambiar el {mate}"
        );
    }
    // Y sí cambia lo publicado, de modo que la igualdad de arriba no sale de
    // haber ejecutado dos veces la misma configuración.
    assert_eq!(pcr_section(&nine.event())["efficiency"].as_f64(), Some(0.9));
    assert_eq!(pcr_section(&four.event())["efficiency"].as_f64(), Some(0.4));
    // Control positivo del mecanismo: **sin** clases declaradas, esa misma
    // escalar sí manda, y cambiarla cambia los bytes. Sin esto, la igualdad de
    // arriba podría venir de que la eficiencia no se use nunca.
    let alone_nine = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 8.0,
        event_extra: amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.9\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    let alone_four = run(&Scenario {
        bases: uniform_bases(20_000),
        coverage: 8.0,
        event_extra: amplified("duplicate_fraction = 0.3\ncycles = 10\nefficiency = 0.4\n"),
        selected: vec!["sample-A"],
        sequencing_errors: "disabled",
    });
    assert_ne!(
        alone_nine.fastq_bytes("sample-A", "R1"),
        alone_four.fastq_bytes("sample-A", "R1"),
        "sin clases, la eficiencia escalar sí manda"
    );
}
