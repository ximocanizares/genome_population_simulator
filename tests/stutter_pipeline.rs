//! F4.8.4 a F4.8.7: el tartamudeo dentro del recorrido técnico.
//!
//! # El fixture saturado, y por qué su aritmética es la que es
//!
//! La comprobación de coherencia del par exige que **R1 y el complemento
//! inverso de R2 cubran los dos el tracto entero**, en todas las parejas. Si
//! solo R1 lo cubriera, la aserción se cumpliría sola y no discriminaría nada.
//!
//! El plan lo pedía con una copia haplotípica de 250 pb, que es el inserto
//! mínimo del perfil. **Eso no es ejecutable**: el perfil sortea el inserto de
//! una normal de media 350 y desviación 20, truncada **por rechazo**, de modo
//! que 250 queda a **cinco desviaciones** de la media. Una copia de 250 pb
//! obliga a acertar exactamente ese valor y cuesta del orden de trece millones
//! de sorteos por pareja; medido, una ejecución de dos parejas tardaba 111
//! segundos. El fixture se ha rehecho con la misma propiedad y coste
//! razonable, y su aritmética está en [`saturated_scenario`].
//!
//! # Las demás elecciones, y qué se perdería sin ellas
//!
//! - **La tasa a uno.** Con la tasa calibrada del dinucleótido, 0,011, un
//!   corpus de noventa moléculas produciría una o dos slips y casi cualquier
//!   aserción pasaría por casualidad. Con la tasa a uno tartamudean todas.
//! - **La tirada corta a la derecha.** Si la tirada larga estuviera a la
//!   derecha, la convención de «la tirada más larga» y la del «extremo
//!   derecho» coincidirían y la prueba no distinguiría entre las dos.
//! - **La tasa intermedia del manifiesto.** Con la tasa a uno el número
//!   esperado de moléculas tartamudeadas y el realizado coinciden, así que un
//!   manifiesto que publicara el esperado pasaría. Con 0,5 se separan.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use noodles_bgzf as bgzf;
use serde_json::Value;

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

/// Longitud de la copia haplotípica del fixture saturado.
///
/// Es el techo del inserto: el muestreador condiciona la normal truncada al
/// intervalo que alguna copia puede alojar. Con 280 pb el inserto cae siempre
/// en `[250, 280]`, y el solapamiento de los dos mates de 150 pb mide
/// `300 − inserto`, es decir **entre 20 y 50 bases**.
const COPY_LENGTH: usize = 280;

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-8-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

/// Los cuatro pesos, siempre declarados: son obligatorios cuando el tartamudeo
/// está activo y no tienen valor por omisión.
const CONTRACTION_BIASED_STEPS: &str = r#"
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
"#;

/// Solo la contracción de una unidad: hace determinista **el paso**, no la
/// decisión de tartamudear, que sigue saliendo del Bernoulli.
const ONLY_MINUS_ONE_STEPS: &str = r#"
[[technical_event.stutter.steps]]
units = -3
weight = 0

[[technical_event.stutter.steps]]
units = -2
weight = 0

[[technical_event.stutter.steps]]
units = -1
weight = 1

[[technical_event.stutter.steps]]
units = 1
weight = 0
"#;

struct Scenario {
    root: PathBuf,
    destination: PathBuf,
}

impl Scenario {
    fn manifest(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("manifest.json")).expect("manifiesto"),
        )
        .expect("JSON del manifiesto")
    }

    fn effective(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("effective-config.json")).expect("configuración"),
        )
        .expect("JSON de la configuración efectiva")
    }

    fn event_manifest(&self) -> Value {
        serde_json::from_slice(
            &fs::read(self.destination.join("events/reads-main/manifest.json"))
                .expect("manifiesto de evento"),
        )
        .expect("JSON del manifiesto de evento")
    }

    fn mates(&self) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
        (
            read_fastq_sequences(
                &self
                    .destination
                    .join("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
            ),
            read_fastq_sequences(
                &self
                    .destination
                    .join("outputs/fastq/reads-main/sample-A_R2.fastq.gz"),
            ),
        )
    }

    /// La verdad publicada **sin su identificador de ejecución**.
    ///
    /// Dos escenarios que difieren en una línea de configuración tienen
    /// identificadores de ejecución distintos por contrato, porque el
    /// identificador se deriva de los bytes del documento. Compararlos crudos
    /// haría fallar la prueba por un motivo que no tiene nada que ver con el
    /// tartamudeo, y taparía el que sí importa. Todo lo demás —los genotipos,
    /// las coordenadas, los alelos— se compara íntegro.
    fn truth_without_execution_id(&self) -> Vec<Value> {
        let text = fs::read_to_string(self.destination.join("outputs/truth.jsonl"))
            .expect("verdad publicada");
        text.lines()
            .map(|line| {
                let mut record: Value = serde_json::from_str(line).expect("registro JSON");
                if let Some(object) = record.as_object_mut() {
                    object.remove("execution_id");
                }
                record
            })
            .collect()
    }

    /// El VCF publicado **sin las líneas que llevan el identificador de
    /// ejecución**, por el mismo motivo que la verdad: dos escenarios que
    /// difieren en una línea de configuración tienen identificadores distintos
    /// por contrato, y compararlos crudos haría fallar la prueba por algo que
    /// no tiene nada que ver con el tartamudeo.
    fn vcf_without_execution_id(&self) -> Vec<String> {
        let text =
            fs::read_to_string(self.destination.join("outputs/truth.vcf")).expect("VCF publicado");
        text.lines()
            // Se descarta **la línea de cabecera concreta** que lleva el
            // identificador, no cualquiera que contenga cierta subcadena. Hoy
            // las dos formas caen sobre la misma línea, pero un filtro por
            // subcadena se lleva por delante cualquier registro futuro que la
            // mencione, y este documento puede crecer.
            .filter(|line| !line.starts_with("##gps_execution_id="))
            .map(str::to_owned)
            .collect()
    }

    fn digest(&self, relative: &str) -> String {
        use sha2::{Digest, Sha256};
        let bytes = fs::read(self.destination.join(relative)).expect("artefacto publicado");
        let mut hasher = Sha256::new();
        hasher.update(&bytes);
        format!("{:x}", hasher.finalize())
    }
}

impl Drop for Scenario {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn read_fastq_sequences(path: &Path) -> Vec<Vec<u8>> {
    let mut text = Vec::new();
    bgzf::io::Reader::new(fs::File::open(path).expect("FASTQ publicado"))
        .read_to_end(&mut text)
        .expect("FASTQ legible");
    String::from_utf8(text)
        .expect("FASTQ ASCII")
        .lines()
        .skip(1)
        .step_by(4)
        .map(|line| line.as_bytes().to_vec())
        .collect()
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect()
}

/// Copias completas **contiguas** del motivo en la tirada más larga.
///
/// Es una contabilidad escrita **aquí**, en la prueba, y no la del productor:
/// comparar la salida del productor con su propio contador sería un oráculo
/// circular y no distinguiría nada.
fn longest_run_units(sequence: &[u8], motif: &[u8]) -> u64 {
    let mut best = 0;
    let mut current = 0;
    let mut offset = 0;
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

/// Escenario saturado: una copia de 280 pb, un individuo haploide y un locus.
///
/// # La aritmética que garantiza que los dos mates ven el tracto entero
///
/// La molécula es `copia[s, s+I)` con inserto `I ∈ [250, 280]` y comienzo
/// `s ∈ [0, 280−I]`. Un tracto que ocupa `[a, a+T)` de la copia y que tras el
/// tartamudeo cambia de longitud en `d` bases queda dentro de los dos mates de
/// 150 pb, **para todo comienzo posible**, exactamente cuando
///
/// ```text
/// a ≥ 130 + d        (lo exige R2, en el peor comienzo s = 280 − I)
/// a ≤ 150 − T − d    (lo exige R1, en el peor comienzo s = 0)
/// ```
///
/// De ahí salen las dos posiciones que usan las pruebas: `a = 135` con `T = 10`
/// admite los cuatro pasos —`d` va de −6 a +2—, y `a = 132` con `T = 16` admite
/// las contracciones. **Ninguna de las dos es un número redondo elegido a ojo**:
/// fuera de esos márgenes hay parejas en las que un mate corta el tracto y la
/// comparación dejaría de tener sentido.
fn saturated_scenario(tract: &str, tract_offset: usize, stutter: &str, coverage: u32) -> Scenario {
    scenario_with_length(tract, tract_offset, stutter, coverage, COPY_LENGTH)
}

fn scenario_with_length(
    tract: &str,
    tract_offset: usize,
    stutter: &str,
    coverage: u32,
    copy_length: usize,
) -> Scenario {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    let bases = format!(
        "{}{tract}{}",
        "G".repeat(tract_offset),
        "G".repeat(copy_length - tract_offset - tract.len())
    );
    assert_eq!(bases.len(), copy_length);
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    // La matriz es obligatoria en cuanto hay población inicial. Su único sitio
    // cae **fuera del tracto**, en el flanco: dentro invadiría el locus y la
    // admisión lo rechazaría por un motivo ajeno al tartamudeo.
    fs::write(
        root.join("population-snvs.tsv"),
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tG\tT\tP1\t0.5\n",
    )
    .expect("matriz");
    let config = root.join("general.toml");
    let start = tract_offset + 1;
    let end = tract_offset + tract.len();
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "supplied"
input_path = "population-snvs.tsv"

[vcf_output]
mode = "text"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {copy_length}

[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = {start}
end = {end}
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{tract}"
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
{stutter}
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
    Scenario { root, destination }
}

/// Escenario con sitio de sobra: 1200 pb, insertos del perfil sin truncar y
/// espacio para una segunda fuente de quimera.
fn roomy_scenario(tract: &str, tract_offset: usize, stutter: &str, coverage: u32) -> Scenario {
    scenario_with_length(tract, tract_offset, stutter, coverage, 1200)
}

/// Bloque de tartamudeo con la tasa del dinucleótido sustituida.
fn stutter_block(rate: f64, steps: &str) -> String {
    format!(
        r#"
[technical_event.stutter]
enabled = true

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
rate = {rate}
{steps}"#
    )
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

/// **F4.8.4 (a). El tracto entero cabe en los dos mates, y el corpus varía.**
///
/// # Qué comprueba de verdad esta prueba, y qué no
///
/// La igualdad entre R1 y el complemento inverso de R2 **no puede fallar por
/// ningún defecto de `src/stutter.rs`**, y conviene decirlo aquí en vez de
/// dejar que el nombre prometa lo que no da. Los dos mates se cortan del
/// **mismo vector de bytes**: R1 es su prefijo y R2 el complemento inverso de
/// su sufijo. Al deshacer ese complemento inverso, la comparación vuelve a
/// mirar el mismo tramo del mismo vector, de modo que la igualdad es una
/// **identidad estructural** del emisor de parejas, no una comprobación del
/// tartamudeo. Se ha medido: aplicar el paso **dos veces** por molécula —una
/// edición gravemente incorrecta— deja esta igualdad intacta.
///
/// Lo que sí discrimina, y por eso la prueba se conserva:
///
/// - que el tracto siga **entero dentro de los dos mates** después del paso,
///   que es la aritmética del fixture y lo que un paso mal colocado rompe;
/// - el **control positivo** del final: con la tasa a uno el corpus tiene que
///   mostrar varios tamaños de tracto. Sin él, un tartamudeo que no ocurriera
///   nunca satisfaría todo lo anterior.
///
/// La afirmación «una decisión por molécula y no por lectura» la caza otra
/// cosa: construir R2 desde una molécula que haya resbalado por su cuenta. Esa
/// mutación sí mata esta prueba, pero hay que hacerla en el emisor de parejas,
/// no aquí.
#[test]
fn both_mates_of_a_pair_show_the_same_unit_count() {
    let scenario = saturated_scenario(
        &"CA".repeat(5),
        135,
        &stutter_block(1.0, CONTRACTION_BIASED_STEPS),
        100,
    );
    let (r1, r2) = scenario.mates();
    assert_eq!(r1.len(), r2.len());
    assert!(r1.len() >= 50, "el corpus tiene que ser representativo");
    let mut observed_counts = std::collections::BTreeSet::new();
    for (first, second) in r1.iter().zip(&r2) {
        let forward = longest_run_units(first, b"CA");
        let reverse = longest_run_units(&reverse_complement(second), b"CA");
        // Identidad estructural del emisor, no una comprobación del
        // tartamudeo: ver la docstring.
        assert_eq!(
            forward, reverse,
            "los dos mates de una pareja tienen que mostrar el mismo número de unidades"
        );
        assert!(
            forward > 0,
            "el fixture saturado deja el tracto entero en los dos mates"
        );
        observed_counts.insert(forward);
    }
    // Control positivo: con la tasa a uno y los cuatro pasos declarados, el
    // corpus **tiene** que contener más de un número de unidades. Sin esto, un
    // tartamudeo que nunca ocurriera satisfaría la aserción anterior.
    assert!(
        observed_counts.len() > 1,
        "la tasa a uno tiene que producir varios tamaños: {observed_counts:?}"
    );
}

/// **F4.8.4 (b). El tartamudeo va después de la quimera.**
///
/// La quimera vuelve a materializar el prefijo del fragmento desde la copia
/// haplotípica y lo compara byte a byte. Una implementación que modificara el
/// fragmento **antes** rompería esa comprobación con un
/// `ChimeraError::PairSourceMismatch`, que aborta la ejecución. Sin una tasa de
/// quimera alta esta prueba no ejercitaría el camino que dice cubrir, y por eso
/// la tasa es uno; el manifiesto lo confirma con el recuento realizado.
#[test]
fn a_high_chimera_rate_together_with_a_high_stutter_rate_publishes_without_error() {
    // **Este caso no usa el fixture saturado, y no puede usarlo.** Formar una
    // quimera exige una segunda fuente que no solape con la primera, y en una
    // copia de 280 pb con insertos de 250 a 280 no cabe ninguna: el
    // muestreador se detiene con `ChimeraError::NoSecondSource` antes de
    // llegar al camino que esta prueba quiere ejercitar. Con 1200 pb sí cabe, y
    // la saturación no hace falta aquí: lo que se comprueba es que los dos
    // artefactos conviven, no la coherencia entre mates.
    let scenario = roomy_scenario(
        &"CA".repeat(5),
        600,
        &format!(
            "chimera_rate = 1.0\n{}",
            stutter_block(1.0, CONTRACTION_BIASED_STEPS)
        ),
        40,
    );
    let event = scenario.event_manifest();
    assert!(
        event["artifacts"]["chimera"]["realized_pairs"]
            .as_u64()
            .expect("quimeras")
            > 0,
        "con la tasa a uno tiene que haber quimeras realizadas, o la prueba no ejercita el camino"
    );
    assert!(
        event["artifacts"]["stutter"]["stuttered_tracts"]
            .as_u64()
            .expect("tartamudeos")
            > 0,
        "y tartamudeo a la vez"
    );
}

/// **F4.8.4 (c). El paso cae en la tirada más larga, no en el extremo derecho.**
///
/// El tracto es impuro: una tirada de **cinco** unidades a la izquierda, una
/// interrupción `TT` y una tirada de **dos** a la derecha. Con la tasa a uno y
/// solo el paso de menos uno declarado, toda molécula pierde exactamente una
/// unidad. Si el paso se aplicara por el extremo derecho, la tirada de dos
/// bajaría a una; aquí se exige que **siga midiendo dos**, byte a byte.
///
/// *La posición importa:* con la tirada larga a la derecha las dos convenciones
/// coincidirían y esta prueba no discriminaría nada.
#[test]
fn the_step_lands_on_the_longest_run_and_not_on_the_rightmost_one() {
    // Cinco unidades, interrupción `TT` y dos unidades: 14 de 16 bases
    // cubiertas, pureza 0,875, por encima del mínimo de admisión de 0,8.
    let tract = format!("{}TT{}", "CA".repeat(5), "CA".repeat(2));
    assert_eq!(tract.len(), 16);
    let scenario = saturated_scenario(&tract, 132, &stutter_block(1.0, ONLY_MINUS_ONE_STEPS), 100);
    let (r1, _) = scenario.mates();
    let survivor = format!("{}TT{}", "CA".repeat(4), "CA".repeat(2));
    let untouched = format!("{}TT{}", "CA".repeat(5), "CA".repeat(2));
    assert!(!r1.is_empty());
    for read in &r1 {
        let text = String::from_utf8(read.clone()).expect("ASCII");
        assert!(
            text.contains(&survivor),
            "la tirada de dos de la derecha tiene que sobrevivir byte a byte: {text}"
        );
        assert!(
            !text.contains(&untouched),
            "la tirada de cinco tiene que haber perdido una unidad: {text}"
        );
    }
}

/// **F4.8.5. Las dos lecturas conservan la longitud del perfil.**
///
/// *Observar solo las lecturas no discrimina* el cambio de longitud de la
/// molécula, porque siguen midiendo 150 pase lo que pase: esa mitad del
/// contrato se comprueba **sobre la molécula**, en la prueba unitaria
/// `a_contraction_shortens_the_molecule_by_the_motif_length` de
/// `src/stutter.rs`. Lo que sí se comprueba aquí es que el cambio de inserto
/// **no** se compensa recortando ni rellenando las lecturas, y que ninguna
/// molécula se descarta: el número de parejas publicadas es el mismo con
/// tartamudeo y sin él.
#[test]
fn the_reads_keep_their_length_and_no_molecule_is_discarded() {
    let quiet = saturated_scenario(&"CA".repeat(5), 135, "", 100);
    let loud = saturated_scenario(
        &"CA".repeat(5),
        135,
        &stutter_block(1.0, ONLY_MINUS_ONE_STEPS),
        100,
    );
    let (quiet_r1, _) = quiet.mates();
    let (loud_r1, loud_r2) = loud.mates();
    assert_eq!(
        quiet_r1.len(),
        loud_r1.len(),
        "el tartamudeo no descarta moléculas"
    );
    for (first, second) in loud_r1.iter().zip(&loud_r2) {
        assert_eq!(first.len(), 150, "R1 conserva la longitud del perfil");
        assert_eq!(second.len(), 150, "R2 conserva la longitud del perfil");
        assert_eq!(
            longest_run_units(first, b"CA"),
            4,
            "cinco unidades menos una son cuatro"
        );
    }
    for read in &quiet_r1 {
        assert_eq!(
            longest_run_units(read, b"CA"),
            5,
            "sin tartamudeo el tracto sale entero"
        );
    }
}

/// **F4.8.6 (c). La verdad del individuo no cambia.**
///
/// El escenario declara un locus cuyo **único alelo es el de la referencia, a
/// frecuencia 1**: el tracto existe —y por tanto tartamudea— pero no hay
/// ninguna variante que publicar. Eso es lo que hace la comparación posible.
/// Con un locus polimórfico la verdad cambiaría legítimamente y exigir su
/// identidad sería pedir lo imposible.
///
/// # Qué mitad de esta prueba puede fallar, y cuál no
///
/// La igualdad de la verdad y del VCF **no puede fallar por ningún defecto de
/// `src/stutter.rs`**, y conviene decirlo aquí igual que se dice en la
/// coherencia del par. El transformador recibe el individuo **prestado e
/// inmutable**: nada dentro del módulo puede escribir en el estado genotípico
/// ni en lo que de él se publica. La invariante está sostenida por el tipo, no
/// por esta comparación. Medido: con el tartamudeo anulado, esta prueba muere
/// **solo** en la desigualdad de los ficheros de lecturas.
///
/// Lo que sí discrimina es esa desigualdad —que el tartamudeo ocurra de
/// verdad— y los dos controles positivos, que impiden que la igualdad se
/// cumpla sobre documentos vacíos.
#[test]
fn stutter_changes_the_reads_and_leaves_the_truth_untouched() {
    let without = saturated_scenario(&"CA".repeat(5), 135, "", 100);
    let with = saturated_scenario(
        &"CA".repeat(5),
        135,
        &stutter_block(1.0, CONTRACTION_BIASED_STEPS),
        100,
    );
    assert_eq!(
        without.truth_without_execution_id(),
        with.truth_without_execution_id(),
        "el tartamudeo vive en las moléculas, no en la verdad del individuo"
    );
    // Control positivo de que la verdad no está vacía: si lo estuviera, la
    // igualdad anterior se cumpliría sola.
    assert!(
        with.truth_without_execution_id().len() > 1,
        "la verdad tiene cabecera y al menos un registro"
    );
    // **El VCF, que el paso enumera uno a uno y ninguna prueba comparaba.**
    assert_eq!(
        without.vcf_without_execution_id(),
        with.vcf_without_execution_id(),
        "el VCF de verdad tampoco cambia: no hay ningún alelo de repetición que publicar"
    );
    assert!(
        with.vcf_without_execution_id()
            .iter()
            .any(|line| !line.starts_with('#')),
        "el VCF tiene que traer algún registro, o su igualdad no compara nada"
    );
    for mate in ["R1", "R2"] {
        assert_ne!(
            without.digest(&format!(
                "outputs/fastq/reads-main/sample-A_{mate}.fastq.gz"
            )),
            with.digest(&format!(
                "outputs/fastq/reads-main/sample-A_{mate}.fastq.gz"
            )),
            "con la tasa a uno los FASTQ tienen que cambiar"
        );
    }
}

/// **F4.8.6 (d). Dos ejecuciones que solo difieren en la tasa.**
///
/// Difieren en los FASTQ y en los recuentos técnicos agregados, y en nada de la
/// verdad. Con la tasa a cero `bernoulli-v1` no sortea: no se consume ni una
/// palabra y los FASTQ salen **byte a byte** como los del escenario que ni
/// siquiera declara el bloque.
#[test]
fn two_runs_that_differ_only_in_the_rate_differ_only_in_reads_and_aggregates() {
    let quiet = saturated_scenario(
        &"CA".repeat(5),
        135,
        &stutter_block(0.0, CONTRACTION_BIASED_STEPS),
        100,
    );
    let loud = saturated_scenario(
        &"CA".repeat(5),
        135,
        &stutter_block(1.0, CONTRACTION_BIASED_STEPS),
        100,
    );
    assert_eq!(
        quiet.truth_without_execution_id(),
        loud.truth_without_execution_id()
    );
    assert_ne!(
        quiet.digest("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
        loud.digest("outputs/fastq/reads-main/sample-A_R1.fastq.gz")
    );
    assert_eq!(
        quiet.event_manifest()["artifacts"]["stutter"]["stuttered_tracts"],
        Value::from(0_u64),
        "con la tasa a cero `bernoulli-v1` no sortea y nada tartamudea"
    );
    assert!(
        loud.event_manifest()["artifacts"]["stutter"]["stuttered_tracts"]
            .as_u64()
            .expect("recuento")
            > 0
    );
    let absent = saturated_scenario(&"CA".repeat(5), 135, "", 100);
    assert_eq!(
        quiet.digest("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
        absent.digest("outputs/fastq/reads-main/sample-A_R1.fastq.gz"),
        "una tasa de cero no puede desplazar ninguna decisión anterior"
    );
}

/// **F4.8.7. Lo que el manifiesto registra del tartamudeo.**
///
/// La suma de la distribución de pasos por locus tiene que valer exactamente lo
/// mismo que el número de moléculas que tartamudearon. Un manifiesto que
/// publicara el número **esperado** —la tasa por las moléculas evaluadas— en
/// vez del realizado rompe esa identidad, y con la tasa intermedia de 0,5 las
/// dos cifras se separan de verdad.
#[test]
fn the_event_manifest_publishes_realized_counts_with_their_provenance() {
    let scenario = saturated_scenario(
        &"CA".repeat(5),
        135,
        &stutter_block(0.5, CONTRACTION_BIASED_STEPS),
        100,
    );
    let event = scenario.event_manifest();
    assert_eq!(event["schema_version"], "3");
    assert_eq!(event["rng_catalog"], "f4-technical-v4");
    let stutter = &event["artifacts"]["stutter"];
    assert_eq!(stutter["algorithm"], "stutter-v1");
    assert_eq!(stutter["run_length_coefficient"]["calibrated"], false);
    assert_eq!(stutter["run_length_coefficient"]["value"], 0.0);
    assert_eq!(stutter["calibration"]["species"], "human");
    assert_eq!(stutter["calibration"]["technology"], "illumina");
    assert_eq!(stutter["calibration"]["column"], "pcr_free");
    assert!(
        stutter["calibration"]["source"]
            .as_str()
            .expect("fuente")
            .contains("Raz")
    );
    assert!(
        stutter["calibration"]["known_mismatch"]
            .as_str()
            .expect("desajuste")
            .contains("motif type")
    );
    let evaluated = stutter["evaluated_tracts"].as_u64().expect("evaluadas");
    let stuttered = stutter["stuttered_tracts"].as_u64().expect("realizadas");
    assert!(evaluated > 0);
    assert!(
        stuttered > 0 && stuttered < evaluated,
        "con la tasa a 0,5 el realizado tiene que separarse del total evaluado"
    );
    let total: u64 = stutter["steps_by_locus"]
        .as_array()
        .expect("distribución")
        .iter()
        .flat_map(|entry| entry["steps"].as_array().expect("pasos"))
        .map(|step| step["count"].as_u64().expect("recuento"))
        .sum();
    assert_eq!(
        total, stuttered,
        "la distribución de pasos suma las moléculas que tartamudearon"
    );
    // **Búsqueda negativa con control positivo:** ningún campo nombra una
    // lectura ni una pareja. El control positivo es que el texto sí contiene el
    // identificador del locus, de modo que la búsqueda mira donde debe.
    let text = serde_json::to_string(stutter).expect("texto");
    assert!(text.contains("STR-1"));
    for forbidden in ["read_id", "pair_id", "template_ordinal", "pre_hop"] {
        assert!(
            !text.contains(forbidden),
            "el manifiesto no puede nombrar {forbidden}"
        );
    }
    assert_external_schema(
        &scenario.destination.join("events/reads-main/manifest.json"),
        "manifest-event-v3.schema.json",
    );
    assert_external_schema(
        &scenario.destination.join("effective-config.json"),
        "effective-config-v6.schema.json",
    );
    let effective = scenario.effective();
    let published = &effective["technical_event"]["stutter"];
    assert_eq!(published["rng_catalog"], "f4-technical-v4");
    assert_eq!(
        published["motif_length_rates"].as_array().map(Vec::len),
        Some(6)
    );
    assert_eq!(published["motif_length_rates"][1]["rate"], 0.5);
    assert_eq!(
        published["motif_length_rates"][0]["rate"], 0.027,
        "las longitudes no sustituidas conservan la cifra calibrada"
    );
    assert_eq!(
        effective["technical_event"]["rng_catalog"],
        "f4-technical-v4"
    );
    assert_eq!(
        scenario.manifest()["technical_events"][0]["rng_catalog"],
        "f4-technical-v4"
    );
}

/// Copias completas contiguas del motivo **inmediatamente después** de un
/// marcador, contadas aquí.
///
/// Con dos tractos del mismo motivo en la misma molécula, contar «la tirada más
/// larga de la lectura» devolvería el máximo de los dos y no distinguiría
/// cuál cambió. El marcador es lo que ata cada cuenta a su locus.
fn units_after(read: &[u8], marker: &[u8], motif: &[u8]) -> Option<u64> {
    let start = read
        .windows(marker.len())
        .position(|window| window == marker)?
        + marker.len();
    let mut units = 0;
    let mut offset = start;
    while offset + motif.len() <= read.len() && &read[offset..offset + motif.len()] == motif {
        units += 1;
        offset += motif.len();
    }
    Some(units)
}

/// Escenario con **dos loci del mismo motivo**, cada uno tras su marcador.
///
/// Los dos caen dentro de R1 en todas las parejas: R1 cubre `copia[s, s+150)`
/// con `s ∈ [0, 30]`, de modo que cualquier tramo dentro de `copia[30, 148)`
/// está siempre entero en la primera lectura, tartamudee o no.
fn two_loci_scenario(stutter: &str) -> Scenario {
    let root = temporary("dos-loci");
    let config = write_two_loci(&root, stutter);
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .expect("publicación");
    Scenario { root, destination }
}

/// El mismo documento **escrito y no ejecutado**, para las pruebas que solo
/// miran el rechazo. Reutiliza `Scenario` como dueño del directorio temporal,
/// con un destino que nunca llega a existir.
fn two_loci_config(stutter: &str) -> (Scenario, PathBuf) {
    let root = temporary("dos-loci");
    let config = write_two_loci(&root, stutter);
    let destination = root.join("published");
    (Scenario { root, destination }, config)
}

fn write_two_loci(root: &Path, stutter: &str) -> PathBuf {
    fs::create_dir_all(root).expect("directorio");
    let tract = "CA".repeat(5);
    let bases = format!(
        "{}TTTT{tract}{}AAAA{tract}{}",
        "G".repeat(36),
        "G".repeat(46),
        "G".repeat(170)
    );
    assert_eq!(bases.len(), COPY_LENGTH);
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    fs::write(
        root.join("population-snvs.tsv"),
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\tG\tT\tP1\t0.5\n",
    )
    .expect("matriz");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[population_snvs]
mode = "supplied"
input_path = "population-snvs.tsv"

[vcf_output]
mode = "text"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {COPY_LENGTH}

[[reference.repeat_loci]]
id = "STR-A"
sequence = "chr1"
start = 41
end = 50
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{tract}"
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-B"
sequence = "chr1"
start = 101
end = 110
motif = "CA"

[[reference.repeat_loci.alleles]]
sequence = "{tract}"
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
target_coverage_wgs = 100
sequencing_errors = "disabled"
{stutter}
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        ),
    )
    .expect("configuración");
    config
}

/// **F4.8.3. La sustitución de tasa por locus llega de verdad a la simulación.**
///
/// # Por qué hacen falta dos loci y por qué del mismo motivo
///
/// Los dos tractos son `CA` de cinco unidades, así que la **tabla por longitud
/// de motivo les da exactamente la misma tasa**. Lo único que los separa es la
/// sustitución: `STR-A` a uno y `STR-B` a cero. Si el modelo ignorara el
/// identificador de locus y usara la tabla, los dos recibirían 0,5 y **los dos
/// mostrarían una mezcla de tamaños**; aquí se exige que uno tartamudee
/// siempre y el otro nunca.
///
/// La tasa de la tabla se deja en 0,5 a propósito, en un valor que **no** es
/// ninguna de las dos sustituciones: con la tabla a cero o a uno, un modelo que
/// ignorara la sustitución podría coincidir con uno de los dos loci por
/// casualidad y la prueba se volvería ambigua.
///
/// El paso es determinista —solo el de menos una unidad tiene peso— para que
/// lo único aleatorio sea **si** se resbala, que es justo lo que la tasa
/// gobierna.
#[test]
fn a_per_locus_rate_substitution_reaches_the_simulation() {
    let scenario = two_loci_scenario(&format!(
        r#"
[technical_event.stutter]
enabled = true

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
rate = 0.5

[[technical_event.stutter.locus_rates]]
locus_id = "STR-A"
rate = 1.0

[[technical_event.stutter.locus_rates]]
locus_id = "STR-B"
rate = 0.0
{ONLY_MINUS_ONE_STEPS}"#
    ));
    let (r1, _) = scenario.mates();
    assert!(r1.len() >= 50, "el corpus tiene que ser representativo");
    for read in &r1 {
        assert_eq!(
            units_after(read, b"TTTT", b"CA"),
            Some(4),
            "el locus sustituido a uno tartamudea en todas las moléculas"
        );
        assert_eq!(
            units_after(read, b"AAAA", b"CA"),
            Some(5),
            "el locus sustituido a cero no tartamudea en ninguna"
        );
    }

    // Y lo mismo en el manifiesto: la tasa efectiva por locus sale de la misma
    // función que gobierna la simulación, de modo que esta aserción muere con
    // la misma mutación.
    let event = scenario.event_manifest();
    let stutter = &event["artifacts"]["stutter"];
    let rate_of = |locus: &str| {
        stutter["locus_parameters"]
            .as_array()
            .expect("parámetros por locus")
            .iter()
            .find(|entry| entry["locus_id"] == locus)
            .map(|entry| {
                (
                    entry["base_rates"][0]["rate"].as_f64().expect("tasa"),
                    entry["rate_substituted"].as_bool().expect("marca"),
                )
            })
            .expect("el locus tiene que estar publicado")
    };
    assert_eq!(rate_of("STR-A"), (1.0, true));
    assert_eq!(rate_of("STR-B"), (0.0, true));

    // La distribución de pasos solo nombra al locus que resbaló. Un modelo que
    // ignorara la sustitución haría aparecer también al otro.
    let stuttering_loci = stutter["steps_by_locus"]
        .as_array()
        .expect("distribución")
        .iter()
        .map(|entry| entry["locus_id"].as_str().expect("locus").to_owned())
        .collect::<Vec<_>>();
    assert_eq!(stuttering_loci, vec!["STR-A".to_owned()]);

    // Control de la convención de recuento: con dos tractos por molécula, los
    // tractos evaluados son **el doble** de las parejas. Si el recuento fuera
    // por molécula, esta igualdad no se cumpliría.
    let evaluated = stutter["evaluated_tracts"].as_u64().expect("evaluados");
    assert_eq!(evaluated, 2 * r1.len() as u64);
    assert_eq!(
        stutter["stuttered_tracts"].as_u64().expect("resbalados"),
        r1.len() as u64,
        "resbala exactamente uno de los dos tractos de cada molécula"
    );
}

/// Una sustitución que nombra un locus **que no existe** se rechaza en vez de
/// tirarse en silencio. Un identificador mal escrito es la forma que toma ese
/// error, y sin esta guarda el manifiesto publicaría obedientemente una
/// sustitución que la simulación nunca usó.
#[test]
fn a_substitution_naming_an_undeclared_locus_is_rejected() {
    let (_owner, config) = two_loci_config(&format!(
        r#"
[technical_event.stutter]
enabled = true

[[technical_event.stutter.locus_rates]]
locus_id = "STR-Z"
rate = 0.5
{CONTRACTION_BIASED_STEPS}"#
    ));
    let rejection =
        load_and_validate_config(&config).expect_err("la configuración tenía que rechazarse");
    let rendered = format!("{rejection:?}");
    assert_eq!(rejection.code(), "GPS-CONF-008");
    assert!(
        rendered.contains("technical_event.stutter.locus_rates.locus_id"),
        "nombra el campo: {rendered}"
    );
    assert!(rendered.contains("STR-Z"), "nombra el locus: {rendered}");
    // Control positivo: el mismo escenario con el locus **bien escrito** se
    // admite. Sin él, una guarda que rechazara toda sustitución pasaría igual.
    let (_accepted_owner, accepted) = two_loci_config(&format!(
        r#"
[technical_event.stutter]
enabled = true

[[technical_event.stutter.locus_rates]]
locus_id = "STR-A"
rate = 0.5
{CONTRACTION_BIASED_STEPS}"#
    ));
    assert!(load_and_validate_config(&accepted).is_ok());
}
