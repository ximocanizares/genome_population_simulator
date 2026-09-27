//! F4.10.1, F4.10.2, F4.10.4 y F4.10.6: la genealogía y su aplicación.
//!
//! # Cómo está partido este fichero
//!
//! - **La genealogía**, sobre un fixture saturado que llama directamente al
//!   resolutor: es la única forma de forzar la divergencia y la probabilidad,
//!   que es lo que el plan pide y lo que hace que la aserción discrimine.
//! - **El ámbito**, comparando ejecuciones completas: lo que se comprueba es
//!   que declarar un locus nuevo no mueve los resbalones de los anteriores.
//! - **La aplicación**, comparando lo escrito en los FASTQ contra una cuenta
//!   hecha a mano.

use std::collections::BTreeSet;
use std::fs;
use std::io::Read;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::clonal_stutter::{
    LineageKey, LineageMolecule, LineageTract, resolve_family_lineage,
};
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::identity::{EventId, LibraryId, RepeatLocusId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::stutter::StutterStepWeights;
use noodles_bgzf as bgzf;
use serde_json::Value;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

// ---------------------------------------------------------------------------
// F4.10.2 — la genealogía mínima, sobre el fixture saturado.
// ---------------------------------------------------------------------------

fn root() -> RootSeed {
    RootSeed::parse("root", SEED).expect("semilla")
}

fn identifiers() -> (EventId, LibraryId, RepeatLocusId) {
    (
        EventId::parse("reads-main".to_owned(), "test").expect("evento"),
        LibraryId::parse("LIB_0001".to_owned(), "test").expect("biblioteca"),
        RepeatLocusId::parse("STR-1".to_owned(), "test").expect("locus"),
    )
}

fn step_weights(weights: [f64; 4]) -> StutterStepWeights {
    StutterStepWeights::new(weights).expect("pesos")
}

/// Pesos de divergencia **forzados**: toda la masa en un ciclo.
///
/// Es la sustitución que el plan admite «solo en prueba», y aquí se hace
/// pasándolos como argumento en vez de tocar el código que se prueba.
fn forced_divergence(cycles: usize, at: usize) -> Vec<f64> {
    let mut weights = vec![0.0_f64; cycles + 1];
    weights[at] = 1.0;
    weights
}

/// Genealogía de una familia de ocho copias sobre un tracto saturado.
fn saturated_family(divergence: &[f64], family_size: u64, cycles: u32) -> Vec<Vec<i64>> {
    let (event_id, library_id, locus_id) = identifiers();
    // **Probabilidad 1 en todos los ciclos.** Con la tasa calibrada, la
    // inmensa mayoría de los ciclos no resbala y cualquier aserción sobre la
    // herencia pasaría por casualidad.
    let tracts = [LineageTract::new(&locus_id, 40, 2, 1.0)];
    let lineage = resolve_family_lineage(
        LineageKey {
            root: root(),
            event_id: &event_id,
            library_id: &library_id,
            family_ordinal: 17,
        },
        family_size,
        divergence,
        cycles,
        &tracts,
        &step_weights([0.05, 0.15, 0.6, 0.2]),
        LineageMolecule {
            length_bases: 400,
            read_length_bases: 150,
        },
    )
    .expect("genealogía");
    (0..family_size as usize)
        .map(|index| lineage.member(index).expect("copia").steps().to_vec())
        .collect()
}

/// **F4.10.2 (a). Con la divergencia en el último ciclo, todas las copias de
/// la familia muestran el mismo número de unidades, distinto del verdadero.**
///
/// # Cuidado: qué está forzado y por qué la aserción discrimina
///
/// Están forzadas **dos** cosas, y las dos hacen falta. La divergencia va toda
/// al último ciclo, de modo que cada copia hereda el linaje **entero** de su
/// progenitor; y la probabilidad de resbalón vale 1 en todos los ciclos, de
/// modo que todos resbalan. **Sin forzar la divergencia y sin probabilidad 1,
/// algunas discordancias dentro de la familia son legítimas y esta aserción no
/// discrimina nada.**
///
/// El plan pedía la probabilidad a 1 en el primer ciclo y 0 en los demás. Aquí
/// va a 1 en todos, que es **más exigente**: cualquier ciclo que no se
/// heredara bien rompería la igualdad, no solo el primero.
///
/// **Debe fallar contra:** un tartamudeo sorteado por pareja o por lectura, que
/// daría copias distintas; y una herencia que copie los resbalones
/// **posteriores** a la divergencia en vez de los anteriores, que con la
/// divergencia en el último ciclo no heredaría nada y también daría copias
/// distintas.
#[test]
fn with_divergence_at_the_last_cycle_every_copy_of_the_family_shows_the_same_units() {
    let cycles = 10_u32;
    let members = saturated_family(
        &forced_divergence(cycles as usize, cycles as usize),
        8,
        cycles,
    );
    let first = members[0][0];
    assert_ne!(
        first, 0,
        "control positivo: con probabilidad 1 el linaje tiene que desplazar el tracto"
    );
    for (index, member) in members.iter().enumerate() {
        assert_eq!(
            member[0], first,
            "la copia {index} tendría que heredar el linaje entero: {member:?}"
        );
    }
}

/// **F4.10.2 (b). Con la divergencia en el ciclo cero, las copias son
/// independientes.**
///
/// Es la variante que el plan exige junto a la anterior, y sin ella la primera
/// no significaría nada: una implementación que diera **siempre** el mismo
/// paso a toda la familia —por ejemplo una que ignorara el índice de copia en
/// la clave— pasaría la de arriba y muere aquí.
#[test]
fn with_divergence_at_cycle_zero_the_copies_are_independent() {
    let cycles = 10_u32;
    let members = saturated_family(&forced_divergence(cycles as usize, 0), 8, cycles);
    let distinct = members
        .iter()
        .map(|member| member[0])
        .collect::<BTreeSet<_>>();
    assert!(
        distinct.len() > 1,
        "ocho linajes independientes no pueden dar todos el mismo paso neto: {members:?}"
    );
}

/// **F4.10.6. Los resbalones heredados se cuentan aparte, y valen cero cuando
/// no hay de quién heredar.**
///
/// **Debe fallar contra** un recuento de heredados que cuente **todos** los
/// resbalones: con una familia de una sola copia daría diez en vez de cero.
#[test]
fn inherited_slips_are_counted_apart_from_the_rest() {
    let (event_id, library_id, locus_id) = identifiers();
    let tracts = [LineageTract::new(&locus_id, 40, 2, 1.0)];
    let weights = step_weights([0.05, 0.15, 0.6, 0.2]);
    let molecule = LineageMolecule {
        length_bases: 400,
        read_length_bases: 150,
    };
    let key = LineageKey {
        root: root(),
        event_id: &event_id,
        library_id: &library_id,
        family_ordinal: 17,
    };
    let alone = resolve_family_lineage(
        key,
        1,
        &forced_divergence(10, 10),
        10,
        &tracts,
        &weights,
        molecule,
    )
    .expect("genealogía");
    assert_eq!(
        alone.slips_by_locus()[&locus_id],
        10,
        "diez ciclos, diez resbalones"
    );
    assert_eq!(
        alone.inherited_slips_by_locus()[&locus_id],
        0,
        "una copia sola no hereda de nadie"
    );
    let family = resolve_family_lineage(
        key,
        8,
        &forced_divergence(10, 10),
        10,
        &tracts,
        &weights,
        molecule,
    )
    .expect("genealogía");
    assert_eq!(family.slips_by_locus()[&locus_id], 80);
    assert_eq!(
        family.inherited_slips_by_locus()[&locus_id],
        70,
        "siete copias heredan sus diez ciclos"
    );
    // Y con la divergencia en cero no se hereda nada aunque la familia sea
    // grande: es el control que separa «heredado» de «hay más de una copia».
    let independent = resolve_family_lineage(
        key,
        8,
        &forced_divergence(10, 0),
        10,
        &tracts,
        &weights,
        molecule,
    )
    .expect("genealogía");
    assert_eq!(independent.slips_by_locus()[&locus_id], 80);
    assert_eq!(independent.inherited_slips_by_locus()[&locus_id], 0);
}

/// **F4.10.4. La guarda se aplica al paso neto acumulado, no a cada paso.**
///
/// Treinta ciclos con probabilidad 1 y solo la contracción de una unidad con
/// peso: una implementación que comprobara cada paso por separado dejaría el
/// tracto en **menos veintiséis** unidades. Con la guarda sobre el acumulado,
/// la tirada nunca baja de una unidad completa.
#[test]
fn the_viability_guard_applies_to_the_accumulated_net_step() {
    let (event_id, library_id, locus_id) = identifiers();
    let run_units = 5_i64;
    let tracts = [LineageTract::new(
        &locus_id,
        run_units.unsigned_abs(),
        2,
        1.0,
    )];
    let lineage = resolve_family_lineage(
        LineageKey {
            root: root(),
            event_id: &event_id,
            library_id: &library_id,
            family_ordinal: 17,
        },
        1,
        &forced_divergence(30, 0),
        30,
        &tracts,
        // Solo la contracción de una unidad tiene peso, de modo que cada
        // resbalón que ocurra resta exactamente uno.
        &step_weights([0.0, 0.0, 1.0, 0.0]),
        LineageMolecule {
            length_bases: 400,
            read_length_bases: 150,
        },
    )
    .expect("genealogía");
    let net = lineage.member(0).expect("copia").steps()[0];
    assert!(
        net < 0,
        "control positivo: con probabilidad 1 y peso en la contracción tiene que contraerse"
    );
    assert_eq!(
        run_units + net,
        1,
        "la tirada se para en una unidad completa y no baja de ahí: neto {net}"
    );
}

// ---------------------------------------------------------------------------
// Las ejecuciones completas.
// ---------------------------------------------------------------------------

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-10-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
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
    fn event(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("events/reads-main/manifest.json"))
                .expect("manifiesto de evento"),
        )
        .expect("JSON")
    }

    fn first_reads(&self) -> Vec<Vec<u8>> {
        let mut text = Vec::new();
        bgzf::io::Reader::new(
            fs::File::open(
                self.destination
                    .join("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
            )
            .expect("FASTQ"),
        )
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

    /// Resbalones de linaje de un locus, tal como el manifiesto los publica.
    fn lineage_slips(&self, locus_id: &str) -> u64 {
        self.event()["artifacts"]["pcr"]["lineage_slips_by_locus"]
            .as_array()
            .expect("recuentos")
            .iter()
            .find(|entry| entry["locus_id"] == locus_id)
            .map(|entry| entry["slips"].as_u64().expect("recuento"))
            .unwrap_or_else(|| panic!("el locus {locus_id} no aparece en el manifiesto"))
    }

    /// Pasos netos de linaje que **de verdad se escribieron**, por locus.
    fn applied_steps(&self, locus_id: &str) -> u64 {
        self.event()["artifacts"]["pcr"]["lineage_applied_steps_by_locus"]
            .as_array()
            .expect("recuentos")
            .iter()
            .find(|entry| entry["locus_id"] == locus_id)
            .map(|entry| entry["steps"].as_u64().expect("recuento"))
            .unwrap_or_else(|| panic!("el locus {locus_id} no aparece en el manifiesto"))
    }

    /// Parejas quiméricas realizadas, tal como el manifiesto las publica.
    fn realized_chimera_pairs(&self) -> u64 {
        self.event()["artifacts"]["chimera"]["realized_pairs"]
            .as_u64()
            .expect("parejas quiméricas realizadas")
    }

    fn inherited_slips(&self, locus_id: &str) -> u64 {
        self.event()["artifacts"]["pcr"]["lineage_inherited_slips_by_locus"]
            .as_array()
            .expect("recuentos")
            .iter()
            .find(|entry| entry["locus_id"] == locus_id)
            .map(|entry| entry["slips"].as_u64().expect("recuento"))
            .unwrap_or_else(|| panic!("el locus {locus_id} no aparece en el manifiesto"))
    }
}

/// Lo que distingue una ejecución de otra dentro de este fichero.
struct Scenario<'a> {
    /// Secuencias de la referencia, en orden, como pares nombre/bases.
    sequences: &'a [(&'a str, String)],
    /// Bloques `[[reference.repeat_loci]]` ya escritos.
    loci: String,
    coverage: u64,
    ploidy: u8,
    /// Bloque del evento que va tras la cobertura.
    event_extra: String,
}

/// Copias completas contiguas del motivo en la tirada más larga.
///
/// Contabilidad escrita **aquí**, no la del productor: comparar la salida del
/// productor con su propio contador no distinguiría nada.
fn longest_run_units(sequence: &[u8], motif: &[u8]) -> u64 {
    let mut best = 0_u64;
    let mut current = 0_u64;
    let mut offset = 0_usize;
    while offset < sequence.len() {
        if offset + motif.len() <= sequence.len()
            && &sequence[offset..offset + motif.len()] == motif
        {
            current += 1;
            best = best.max(current);
            offset += motif.len();
        } else {
            current = 0;
            offset += 1;
        }
    }
    best
}

fn run(scenario: &Scenario<'_>) -> Published {
    let root = temporary();
    fs::create_dir_all(&root).expect("directorio");
    let fasta = scenario
        .sequences
        .iter()
        .map(|(name, bases)| format!(">{name}\n{bases}\n"))
        .collect::<String>();
    fs::write(root.join("reference.fa"), fasta).expect("FASTA");
    let first = scenario.sequences[0].0;
    fs::write(
        root.join("population-snvs.tsv"),
        format!("contig\tposition\tref\talt\tpopulation_id\talt_af\n{first}\t20\tG\tT\tP1\t0.5\n"),
    )
    .expect("matriz");
    let declarations = scenario
        .sequences
        .iter()
        .map(|(name, bases)| {
            format!(
                "\n[[reference.sequences]]\nid = \"{name}\"\nclass = \"chromosome\"\nlength = {}\n",
                bases.len()
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
ploidy = {ploidy}

[population_snvs]
mode = "supplied"
input_path = "population-snvs.tsv"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"
{declarations}{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = {ploidy}
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = {coverage}
sequencing_errors = "disabled"
{event_extra}
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#,
            ploidy = scenario.ploidy,
            loci = scenario.loci,
            coverage = scenario.coverage,
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

/// Un bloque de locus con su alelo exacto.
fn locus_block(id: &str, sequence: &str, start: usize, motif: &str, allele: &str) -> String {
    format!(
        "\n[[reference.repeat_loci]]\nid = \"{id}\"\nsequence = \"{sequence}\"\nstart = {start}\n\
         end = {end}\nmotif = \"{motif}\"\n\n[[reference.repeat_loci.alleles]]\n\
         sequence = \"{allele}\"\nfrequency = 1.0\n",
        start = start + 1,
        end = start + allele.len(),
    )
}

/// Bases con un tracto en el desplazamiento 60 y `G` en todo lo demás.
///
/// # La aritmética que garantiza que R1 ve el tracto entero, siempre
///
/// La molécula es `copia[s, s+I)` con `I ∈ [250, 300]` y `s ∈ [0, 300−I]`, de
/// modo que `s ≤ 50`. R1 es `copia[s, s+150)`. Un tracto en `[a, a+T)` cabe en
/// R1 para **todo** comienzo exactamente cuando `a ≥ 50` y `a + T ≤ 150`. Con
/// `a = 60` y `T ≤ 24` las dos se cumplen con holgura.
fn bases_with_tract(tract: &str) -> String {
    let head = format!("{}{tract}", "G".repeat(60));
    let bases = format!("{head}{}", "G".repeat(300 - head.len()));
    assert_eq!(bases.len(), 300);
    bases
}

/// Bloque del tartamudeo de base con la tasa del dinucleótido sustituida y un
/// solo paso con peso.
fn base_stutter(rate: f64, step: i64) -> String {
    let weights = [-3_i64, -2, -1, 1]
        .into_iter()
        .map(|units| {
            let weight = if units == step { 1.0 } else { 0.0 };
            format!("\n[[technical_event.stutter.steps]]\nunits = {units}\nweight = {weight}\n")
        })
        .collect::<String>();
    format!(
        "\n[technical_event.stutter]\nenabled = true\n\n\
         [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = {rate}\n\
         {weights}"
    )
}

/// Bloque de tartamudeo de base con la tasa del dinucleótido sustituida y los
/// pesos declarados uno a uno.
fn base_stutter_with_weights(rate: f64, weights: [f64; 4]) -> String {
    let declared = [-3_i64, -2, -1, 1]
        .into_iter()
        .zip(weights)
        .map(|(units, weight)| {
            format!("\n[[technical_event.stutter.steps]]\nunits = {units}\nweight = {weight}\n")
        })
        .collect::<String>();
    format!(
        "\n[technical_event.stutter]\nenabled = true\n\n\
         [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = {rate}\n\
         {declared}"
    )
}

/// Sustitución de la tasa de un locus, que **prevalece sin modular**.
///
/// Hace falta en los fixtures cuyo tracto es corto **a propósito** —el régimen
/// en el que la guarda de viabilidad muerde y las dos partes dejan de
/// conmutar—. Desde que el modelo depende de la longitud del tracto, poner la
/// tabla por longitud de motivo a uno ya no significa «siempre»: en un tracto
/// de cuatro o cinco unidades la curva vale unas 266 veces menos que en su
/// referencia. La sustitución por locus es la forma que el modelo tiene de
/// decir «este locus tartamudea con probabilidad `rate`, mida lo que mida», y
/// es justo lo que estos fixtures quieren decir.
fn locus_rate(locus_id: &str, rate: f64) -> String {
    format!("\n[[technical_event.stutter.locus_rates]]\nlocus_id = \"{locus_id}\"\nrate = {rate}\n")
}

/// Bloque de la amplificación con la probabilidad por ciclo sustituida.
fn amplified(duplicate_fraction: f64, cycles: u32, per_cycle: Option<f64>) -> String {
    let substitution = per_cycle.map_or_else(String::new, |probability| {
        format!(
            "\n[[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 2\n\
             probability = {probability}\n"
        )
    });
    format!(
        "library_preparation = \"pcr\"\n\n[technical_event.pcr]\n\
         duplicate_fraction = {duplicate_fraction}\ncycles = {cycles}\nefficiency = 0.9\n\
         {substitution}"
    )
}

// ---------------------------------------------------------------------------
// F4.10.1 — el ámbito se indexa por locus, no por posición.
// ---------------------------------------------------------------------------

/// **Añadir un locus en otro cromosoma no cambia los resbalones de los loci
/// previos.**
///
/// Es la comprobación que el plan pide con esas palabras. El segundo
/// cromosoma existe en **las dos** ejecuciones y lo único que cambia es que la
/// segunda lo **declara** como locus: así el tamaño del genoma, el número de
/// parejas y las moléculas son idénticos, y la única diferencia posible viene
/// del ámbito aleatorio.
///
/// # Lo que esta prueba **no** caza, dicho sin adornos
///
/// El plan afirma que esta es «la única que caza un ámbito indexado por la
/// posición del tracto dentro de la molécula». **No lo es, y conviene que
/// conste.** Un tracto en otro cromosoma nunca comparte molécula con el
/// primero, de modo que una implementación que indexara por la posición del
/// tracto dentro de su molécula le seguiría dando el índice cero al primero y
/// pasaría esta prueba sin inmutarse. Quien sí lo caza es la prueba siguiente,
/// que añade el locus **dentro de la misma molécula**.
#[test]
fn declaring_a_locus_on_another_chromosome_does_not_move_the_previous_slips() {
    let sequences = [
        ("chr1", bases_with_tract(&"CA".repeat(10))),
        ("chr2", bases_with_tract(&"GT".repeat(10))),
    ];
    let first = locus_block("STR-1", "chr1", 60, "CA", &"CA".repeat(10));
    let second = locus_block("STR-2", "chr2", 60, "GT", &"GT".repeat(10));
    let alone = run(&Scenario {
        sequences: &sequences,
        loci: first.clone(),
        coverage: 600,
        ploidy: 1,
        event_extra: format!("{}{}", amplified(0.0, 10, None), base_stutter(0.5, -1)),
    });
    let together = run(&Scenario {
        sequences: &sequences,
        loci: format!("{first}{second}"),
        coverage: 600,
        ploidy: 1,
        event_extra: format!("{}{}", amplified(0.0, 10, None), base_stutter(0.5, -1)),
    });
    assert!(
        alone.lineage_slips("STR-1") > 0,
        "control positivo: tiene que haber resbalones que comparar"
    );
    assert_eq!(
        alone.lineage_slips("STR-1"),
        together.lineage_slips("STR-1"),
        "declarar un locus en otro cromosoma ha movido los resbalones del primero"
    );
    // Y el locus nuevo sí produce los suyos, de modo que la igualdad de arriba
    // no pasa por no haber añadido nada.
    assert!(together.lineage_slips("STR-2") > 0);
}

/// **Añadir un locus a la izquierda, dentro de la misma molécula, tampoco
/// mueve los resbalones del que ya había.**
///
/// **Esta es la prueba que caza un ámbito indexado por la posición del tracto
/// dentro de la molécula.** Con el locus nuevo a la izquierda, el que ya había
/// pasa de ser el primer tracto de la molécula a ser el segundo: una
/// implementación que usara ese índice en la clave cambiaría todas sus
/// decisiones, y sus resbalones dejarían de coincidir.
#[test]
fn declaring_a_locus_to_the_left_of_another_does_not_move_its_slips() {
    // El homopolímero va a la izquierda del dinucleótido, dentro de la misma
    // molécula. Los flancos son `G`, de modo que ninguno de los dos extiende
    // la tirada del otro.
    let bases = {
        let head = format!(
            "{}{}GGGG{}",
            "G".repeat(40),
            "A".repeat(12),
            "CA".repeat(10)
        );
        format!("{head}{}", "G".repeat(300 - head.len()))
    };
    let sequences = [("chr1", bases)];
    let dinucleotide = locus_block("STR-DI", "chr1", 56, "CA", &"CA".repeat(10));
    let homopolymer = locus_block("STR-HOMO", "chr1", 40, "A", &"A".repeat(12));
    let event = format!("{}{}", amplified(0.0, 10, None), base_stutter(0.5, -1));
    let alone = run(&Scenario {
        sequences: &sequences,
        loci: dinucleotide.clone(),
        coverage: 600,
        ploidy: 1,
        event_extra: event.clone(),
    });
    let together = run(&Scenario {
        sequences: &sequences,
        loci: format!("{homopolymer}{dinucleotide}"),
        coverage: 600,
        ploidy: 1,
        event_extra: event,
    });
    assert!(
        alone.lineage_slips("STR-DI") > 0,
        "control positivo: tiene que haber resbalones que comparar"
    );
    assert_eq!(
        alone.lineage_slips("STR-DI"),
        together.lineage_slips("STR-DI"),
        "declarar un locus a la izquierda ha movido los resbalones del de la derecha"
    );
    // **El locus añadido tiene que entrar de verdad en la genealogía.** Sin
    // este control la igualdad de arriba no afirma nada: si el homopolímero se
    // quedara fuera —por caer donde no hay tirada, o por no llegar a la etapa—
    // las dos ejecuciones serían la misma y la igualdad se cumpliría sola.
    assert!(
        together.lineage_slips("STR-HOMO") > 0,
        "el locus añadido a la izquierda tiene que producir sus propios resbalones"
    );
}

// ---------------------------------------------------------------------------
// F4.10.4 — el orden de aplicación.
// ---------------------------------------------------------------------------

/// **Con el paso de PCR fijado y la tasa de base a 1 con su paso fijado, el
/// resultado es la suma de los dos.**
///
/// Un ciclo con probabilidad 1 da un resbalón de linaje; la tasa de base a 1
/// da otro; los dos pesos dejan una sola alternativa, `−1`. El tracto
/// verdadero tiene **dieciocho unidades** y **todas** las lecturas tienen que
/// mostrar dieciséis.
///
/// # Por qué dieciocho y no diecisiete
///
/// Porque el paso de base se decide **después** del de linaje y **sobre la
/// molécula ya actualizada**: cuando le toca, el tracto tiene una unidad
/// menos. Con dieciocho verdaderas, la decisión de base cae sobre diecisiete,
/// que son 34 pares de base, la **longitud de tracto de referencia** de los
/// dinucleótidos, donde el factor de modulación vale exactamente `1.0`. Así la
/// tasa de base a uno significa de verdad «siempre».
///
/// Con diecisiete verdaderas la decisión de base caía sobre dieciséis, donde
/// el factor vale 0,794, y una de cada cinco lecturas se quedaba sin el
/// segundo paso. Se midió: la prueba falló exactamente así. Es un recordatorio
/// útil de que la modulación se lee **en el momento de decidir**, no al
/// principio.
///
/// **Debe fallar contra** una entrega que aplique solo una de las dos partes,
/// que daría diecisiete.
#[test]
fn the_two_parts_add_up_on_every_published_read() {
    let sequences = [("chr1", bases_with_tract(&"CA".repeat(18)))];
    let published = run(&Scenario {
        sequences: &sequences,
        loci: locus_block("STR-1", "chr1", 60, "CA", &"CA".repeat(18)),
        coverage: 300,
        ploidy: 1,
        event_extra: format!("{}{}", amplified(0.0, 1, Some(1.0)), base_stutter(1.0, -1)),
    });
    let reads = published.first_reads();
    assert!(
        !reads.is_empty(),
        "control positivo: tiene que haber lecturas"
    );
    for read in &reads {
        assert_eq!(
            longest_run_units(read, b"CA"),
            16,
            "dieciocho unidades menos el paso de linaje menos el de base son dieciséis"
        );
    }
}

/// **Tras dos pasos sucesivos, la tirada más larga se recalcula entre uno y
/// otro.**
///
/// El tracto es **impuro** a propósito: cinco unidades, una interrupción `TT`
/// y cuatro unidades. El paso de linaje, de tres unidades, deja la primera
/// tirada en dos, de modo que **la más larga pasa a ser la segunda**; el paso
/// de base tiene que caer allí.
///
/// - Correcto: `CACA` + `TT` + `CA`.
/// - Una implementación que calcule la tirada más larga **una sola vez** y
///   aplique los dos pasos en la misma deja `CACACACA`, sin la interrupción.
///
/// Con un tracto **puro** las dos convenciones coinciden y la prueba no
/// discriminaría nada: de ahí la interrupción.
///
/// # Lo que esta prueba **no** ancla, dicho sin adornos
///
/// **No distingue el orden entre las dos partes**, y no puede: aquí los dos
/// pasos están forzados **al mismo valor**, `−3`, de modo que intercambiarlos
/// aplica el mismo par de pasos sobre el mismo par de tiradas y escribe los
/// mismos bytes. Se midió.
///
/// Quien sí distingue el orden es
/// [`the_lineage_step_goes_before_the_base_decision`], con un tracto corto y
/// los pesos cargados a la contracción: ahí la guarda de viabilidad del segundo
/// paso depende de lo que el primero deja, y las dos partes dejan de conmutar.
/// La primera versión de esta entrega llegó a declarar esa separación
/// imposible; **lo es solo mientras la guarda no muerda**, y queda dicho aquí
/// para que nadie vuelva a concluir lo mismo desde este fixture.
#[test]
fn the_longest_run_is_recomputed_between_the_two_steps() {
    let tract = format!("{}TT{}", "CA".repeat(5), "CA".repeat(4));
    let sequences = [("chr1", bases_with_tract(&tract))];
    let published = run(&Scenario {
        sequences: &sequences,
        loci: locus_block("STR-1", "chr1", 60, "CA", &tract),
        coverage: 300,
        ploidy: 1,
        event_extra: format!(
            "{}{}{}",
            amplified(0.0, 1, Some(1.0)),
            base_stutter(1.0, -3),
            locus_rate("STR-1", 1.0)
        ),
    });
    let reads = published.first_reads();
    assert!(!reads.is_empty());
    for read in &reads {
        let text = String::from_utf8(read.clone()).expect("ASCII");
        assert!(
            text.contains("GCACATTCAG"),
            "el segundo paso tiene que caer en la otra tirada: {text}"
        );
        assert!(
            !text.contains("GCACACACAG"),
            "una tirada calculada una sola vez borraría la interrupción: {text}"
        );
    }
}

/// **El paso de linaje va primero y la decisión de base después, y aquí se
/// nota.**
///
/// # Por qué este fixture sí separa los dos órdenes, y el anterior no
///
/// La guarda de viabilidad **no es fija: depende del estado que deja el primer
/// paso**. La parte de base calcula sus candidatos con la tirada que encuentra,
/// de modo que si el primer paso la deja corta, el segundo pierde
/// alternativas. Las dos partes **dejan de conmutar en cuanto la guarda
/// muerde**, y eso ocurre con tractos cortos y pesos cargados a la contracción.
///
/// El tracto tiene **cuatro** unidades, un ciclo con probabilidad de resbalón
/// 1, la tasa de base a 1 y los pesos cargados: `−3` con 0,9 y `+1` con 0,1.
///
/// - **Orden vigente.** El linaje sortea sobre una tirada de cuatro, donde los
///   cuatro pasos caben. Si sale `−3` —nueve de cada diez veces— la tirada
///   queda en **una** unidad, y ahí la parte de base ya solo puede expandir,
///   porque cualquier contracción la vaciaría: sale `+1` sin sortear siquiera,
///   y la lectura muestra **dos**. Si sale `+1`, la tirada queda en cinco y la
///   base contrae a dos o expande a seis. **Ninguna lectura puede mostrar una
///   unidad.**
/// - **Orden invertido.** La base sortea primero sobre cuatro y contrae a una;
///   entonces el paso neto `−3` del linaje ya no cabe y se descarta, y la
///   lectura muestra **una** unidad. Ocurre en ocho de cada diez parejas.
///
/// De ahí la aserción: **cero lecturas con una unidad**. Se midió que el orden
/// invertido las produce en 237 de 300.
///
/// # Por qué el fixture del tracto impuro no servía
///
/// Allí los dos pasos están forzados **al mismo valor**, `−3`, y con eso los
/// dos órdenes aplican el mismo par de pasos sobre el mismo par de tiradas y
/// escriben los mismos bytes. Fue el error de la primera versión de esta
/// entrega, que llegó a declarar la separación imposible: lo es solo mientras
/// la guarda no muerda.
#[test]
fn the_lineage_step_goes_before_the_base_decision() {
    let sequences = [("chr1", bases_with_tract(&"CA".repeat(4)))];
    let published = run(&Scenario {
        sequences: &sequences,
        loci: locus_block("STR-1", "chr1", 60, "CA", &"CA".repeat(4)),
        coverage: 300,
        ploidy: 1,
        event_extra: format!(
            "{}{}{}",
            amplified(0.0, 1, Some(1.0)),
            base_stutter_with_weights(1.0, [0.9, 0.0, 0.0, 0.1]),
            locus_rate("STR-1", 1.0)
        ),
    });
    let reads = published.first_reads();
    assert!(
        !reads.is_empty(),
        "control positivo: tiene que haber lecturas"
    );
    let mut counts: std::collections::BTreeMap<u64, u64> = std::collections::BTreeMap::new();
    for read in &reads {
        *counts.entry(longest_run_units(read, b"CA")).or_insert(0) += 1;
    }
    // **La aserción que discrimina, y por qué va sola.** El orden invertido
    // deja 237 de 300 lecturas en una unidad; el vigente, ninguna.
    assert_eq!(
        counts.get(&1).copied().unwrap_or(0),
        0,
        "con el linaje primero ninguna lectura puede quedarse en una unidad: {counts:?}"
    );
    // **Esta de aquí es más fuerte que la anterior y la implica**, porque un
    // corpus que solo muestre dos o seis unidades no puede tener ninguna de
    // una. No es un control positivo, aunque lo pareciera: queda dicho para que
    // nadie borre la de arriba creyendo que sobra. Se conservan las dos porque
    // la primera nombra **el** valor que el orden invertido produce, y un
    // fallo que la señale a ella se lee de un vistazo.
    assert!(
        counts.keys().all(|units| *units == 2 || *units == 6),
        "el fixture solo admite dos o seis unidades: {counts:?}"
    );
    // Controles positivos de verdad: la rama que estrecha la guarda dominó, de
    // modo que el régimen en el que las dos partes dejan de conmutar se
    // ejercitó, y el linaje llegó a escribir.
    let two = counts.get(&2).copied().unwrap_or(0);
    assert!(
        two * 10 >= (reads.len() as u64) * 9,
        "la rama que estrecha la guarda tiene que dominar: {counts:?}"
    );
    assert!(
        published.applied_steps("STR-1") > 0,
        "control positivo: el linaje tiene que haber escrito algún paso"
    );
}

/// **Un corpus con tasa de quimera alta no rompe la comprobación de
/// procedencia de la pareja.**
///
/// La quimera rematerializa el prefijo del fragmento desde la copia
/// haplotípica y lo compara byte a byte; el tartamudeo, clonal o de base, va
/// **después**. Esta prueba existe porque poner el paso de linaje antes de la
/// quimera devolvería `ChimeraError::PairSourceMismatch` y abortaría la
/// ejecución.
#[test]
fn a_high_chimera_rate_does_not_break_the_pair_source_check() {
    let sequences = [("chr1", bases_with_tract(&"CA".repeat(10)))];
    let published = run(&Scenario {
        sequences: &sequences,
        loci: locus_block("STR-1", "chr1", 60, "CA", &"CA".repeat(10)),
        coverage: 300,
        ploidy: 2,
        event_extra: format!(
            "chimera_rate = 0.5\n{}{}",
            amplified(0.2, 10, None),
            base_stutter(0.5, -1)
        ),
    });
    assert!(
        !published.first_reads().is_empty(),
        "la ejecución tiene que llegar a publicar"
    );
    // **Los dos controles positivos, sin los cuales esto no afirma nada.** Una
    // ejecución en la que no hubiera ocurrido ni una quimera, o en la que el
    // linaje no hubiera escrito ni un paso, publicaría lecturas igual y la
    // aserción de arriba pasaría sin haber ejercitado el cruce que la prueba
    // dice cubrir.
    assert!(
        published.realized_chimera_pairs() > 0,
        "con tasa 0,5 tiene que haber quimeras realizadas"
    );
    assert!(
        published.applied_steps("STR-1") > 0,
        "tiene que haberse escrito algún paso de linaje"
    );
}

/// **Un paso neto que no cabe en la molécula quimérica no se escribe.**
///
/// # Qué guarda ata esta prueba, y por qué no la ataba ninguna
///
/// El paso neto del linaje se sortea sobre **la molécula original de la
/// familia**; al aplicarlo se vuelve a comprobar contra **la molécula que se va
/// a emitir**, que con una quimera puede no ser la misma. Esa reverificación no
/// la ejercitaba nadie.
///
/// # Por qué este fixture la dispara, y por qué hacen falta menos condiciones
/// de las que la primera versión de esta entrega supuso
///
/// Se creyó que hacían falta tres cosas a la vez: un locus heterocigoto, un
/// punto de rotura dentro del tracto y una tirada resultante demasiado corta.
/// **Con dos basta.** El sufijo de una quimera viene de un desplazamiento **al
/// azar** sobre la copia hermana, de modo que en cuanto el punto de rotura cae
/// dentro del tracto la tirada queda cortada; no hace falta que los dos alelos
/// midan distinto.
///
/// Los números del fixture: dinucleótido de **cuatro** unidades —para que una
/// contracción de tres deje la tirada al borde—, un ciclo con probabilidad de
/// resbalón 1 —todas las copias resbalan—, **todo** el peso en la contracción
/// de tres unidades, tartamudeo de base apagado y tasa de quimera 1.
///
/// # El control es la mitad que hace que esto signifique algo
///
/// La misma ejecución **sin quimeras**: ahí los pasos escritos tienen que valer
/// exactamente los resbalones sorteados. Sin ese control, una implementación
/// que simplemente perdiera pasos por cualquier motivo satisfaría la
/// desigualdad de arriba y nadie sabría que la diferencia no era de la quimera.
///
/// # Por qué la desigualdad NO basta, y por qué hacen falta dos escenarios
///
/// Todo medido, nada razonado. La reverificación tiene **dos mitades**: que el
/// tracto siga estando en la ventana —`run_units > 0`— y que el paso neto
/// **quepa** en la tirada que queda. Se instrumentó el descarte y las dos
/// disparan por igual con contracciones: veintiséis parejas por pasada cada
/// una.
///
/// Una primera versión se quedaba en `escritos < sorteados` y **sobrevivía a
/// suprimir la mitad de la viabilidad**, porque los descartes de la otra mitad
/// bastaban para que la desigualdad siguiera cumpliéndose. Con las dos mitades
/// se escriben **248** de 300; suprimiendo la viabilidad, **274**. De ahí la
/// cifra congelada.
///
/// Y con contracciones la otra mitad **tampoco se distingue**, porque
/// `step_is_viable` ya excluye por su cuenta una contracción sobre una tirada
/// de cero unidades. La mitad del tracto ausente solo manda con **expansiones**,
/// donde la viabilidad no filtra nada y lo único que impide insertar copias de
/// un motivo que no está es esa comprobación. De ahí el **segundo escenario**,
/// con todo el peso en `+1`: **274** de 300 con quimera y **300** sin ella.
///
/// **Ninguno de los dos literales se actualiza sin más si algún día cambian**:
/// o el cambio mueve la aleatoriedad de estos escenarios —y entonces hay que
/// decidirlo por escrito— o es un defecto.
#[test]
fn a_net_step_that_does_not_fit_the_chimeric_molecule_is_not_written() {
    let tract = "CA".repeat(4);
    let sequences = [("chr1", bases_with_tract(&tract))];
    let scenario = |chimera_rate: f64| Scenario {
        sequences: &sequences,
        loci: locus_block("STR-1", "chr1", 60, "CA", &tract),
        coverage: 300,
        ploidy: 2,
        event_extra: format!(
            "chimera_rate = {chimera_rate}\n{}{}",
            amplified(0.0, 1, Some(1.0)),
            base_stutter_with_weights(0.0, [1.0, 0.0, 0.0, 0.0])
        ),
    };
    let chimeric = run(&scenario(1.0));
    let clean = run(&scenario(0.0));

    let chimeric_slips = chimeric.lineage_slips("STR-1");
    let chimeric_written = chimeric.applied_steps("STR-1");
    let clean_slips = clean.lineage_slips("STR-1");
    let clean_written = clean.applied_steps("STR-1");

    // Controles positivos: las quimeras ocurrieron de verdad en una ejecución y
    // en la otra no, y en las dos se sortearon resbalones.
    assert!(
        chimeric.realized_chimera_pairs() > 0,
        "con tasa 1 tiene que haber quimeras realizadas"
    );
    assert_eq!(
        clean.realized_chimera_pairs(),
        0,
        "el control no puede tener quimeras"
    );
    assert!(chimeric_slips > 0 && clean_slips > 0, "control positivo");

    // **Sin quimeras, todo paso sorteado se escribe.**
    assert_eq!(
        clean_written, clean_slips,
        "sin quimeras la reverificación no puede descartar ningún paso"
    );
    // **Con quimeras, algunos no caben y la reverificación los descarta.**
    assert!(
        chimeric_written < chimeric_slips,
        "la reverificación tiene que descartar algún paso: sorteados {chimeric_slips}, \
         escritos {chimeric_written}"
    );
    // **Y la cifra exacta**, que es lo que distingue la mitad de la viabilidad.
    // Ver la nota de la documentación: la desigualdad sola sobrevive a
    // suprimirla, y entonces da 274.
    assert_eq!(
        chimeric_written, 248,
        "pasos escritos congelados de este escenario; ver la documentación antes de tocar \
         este número"
    );

    // ------------------------------------------------------------------
    // Segundo escenario: **expansiones**, que es donde manda la otra mitad.
    // ------------------------------------------------------------------
    let expanding = |chimera_rate: f64| Scenario {
        sequences: &sequences,
        loci: locus_block("STR-1", "chr1", 60, "CA", &tract),
        coverage: 300,
        ploidy: 2,
        event_extra: format!(
            "chimera_rate = {chimera_rate}\n{}{}",
            amplified(0.0, 1, Some(1.0)),
            base_stutter_with_weights(0.0, [0.0, 0.0, 0.0, 1.0])
        ),
    };
    let grown = run(&expanding(1.0));
    let grown_clean = run(&expanding(0.0));
    assert_eq!(
        grown_clean.applied_steps("STR-1"),
        grown_clean.lineage_slips("STR-1"),
        "sin quimeras ninguna expansión puede descartarse"
    );
    assert_eq!(
        grown.applied_steps("STR-1"),
        274,
        "expansiones escritas congeladas: con la comprobación del tracto ausente suprimida \
         se escribirían las 300, porque la viabilidad no filtra una expansión"
    );
}

/// **El recuento de pasos de linaje escritos cuenta lo escrito, no lo
/// sorteado.**
///
/// Con un solo ciclo a probabilidad 1, la parte de base apagada y un único paso
/// con peso, **todas** las parejas reciben exactamente un paso de linaje de una
/// unidad: el recuento publicado tiene que valer el número de parejas, y todas
/// las lecturas tienen que mostrar nueve unidades de las diez verdaderas.
///
/// Lo que esto ata es el recuento **a los bytes**, y solo eso. **Lo que NO
/// caza, medido:** un contador que sumara los resbalones **sorteados** en vez
/// de los escritos pasa por aquí sin inmutarse, porque en este fixture
/// coinciden —un resbalón por copia y todos se escriben—. Quien lo caza es
/// `a_net_step_that_does_not_fit_the_chimeric_molecule_is_not_written`, donde
/// la quimera hace que dejen de coincidir. Y un contador que doblara la cuenta
/// muere por el tope del manifiesto, que exige que no haya más pasos escritos
/// que parejas.
#[test]
fn the_applied_lineage_count_matches_what_the_reads_show() {
    let sequences = [("chr1", bases_with_tract(&"CA".repeat(10)))];
    let published = run(&Scenario {
        sequences: &sequences,
        loci: locus_block("STR-1", "chr1", 60, "CA", &"CA".repeat(10)),
        coverage: 300,
        ploidy: 1,
        event_extra: format!("{}{}", amplified(0.0, 1, Some(1.0)), base_stutter(0.0, -1)),
    });
    let reads = published.first_reads();
    assert!(!reads.is_empty());
    let displaced = reads
        .iter()
        .filter(|read| longest_run_units(read, b"CA") == 9)
        .count() as u64;
    assert_eq!(
        displaced,
        reads.len() as u64,
        "todas las lecturas tienen que mostrar nueve unidades"
    );
    assert_eq!(
        published.applied_steps("STR-1"),
        reads.len() as u64,
        "el recuento publicado tiene que valer lo que las lecturas muestran"
    );
}

// ---------------------------------------------------------------------------
// F4.10.6 — los recuentos con familias grandes.
// ---------------------------------------------------------------------------

/// **Con familias grandes, los resbalones heredados son mayores que cero.**
///
/// Es la otra mitad del par que lo distingue de un recuento que contara todos
/// los resbalones; la mitad del cero vive en `tests/clonal_composition.rs`,
/// con familias de tamaño 1.
#[test]
fn large_families_inherit_slips() {
    // Diecisiete unidades: la longitud de tracto de referencia, donde el
    // factor de modulación vale 1. Esta prueba usa la probabilidad por ciclo
    // **derivada** de las dos columnas, no una sustituida, y con un tracto más
    // corto la modulación la reduciría hasta dejar sin resbalones el corpus.
    let sequences = [("chr1", bases_with_tract(&"CA".repeat(17)))];
    let published = run(&Scenario {
        sequences: &sequences,
        loci: locus_block("STR-1", "chr1", 60, "CA", &"CA".repeat(17)),
        coverage: 600,
        ploidy: 1,
        event_extra: format!("{}{}", amplified(0.5, 10, None), base_stutter(0.011, -1)),
    });
    let slips = published.lineage_slips("STR-1");
    let inherited = published.inherited_slips("STR-1");
    assert!(slips > 0, "control positivo: tiene que haber resbalones");
    assert!(
        inherited > 0,
        "con fracción de duplicados 0,5 tiene que haber familias que hereden"
    );
    assert!(
        inherited < slips,
        "no todos los resbalones pueden ser heredados: la copia 0 de cada familia sortea los suyos"
    );
}
