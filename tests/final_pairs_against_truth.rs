//! REQ-VAL-03: cada pareja final del FASTQ es la que se espera de su clave
//! técnica, de la verdad biológica publicada y de las transformaciones
//! autorizadas.
//!
//! Decisión del responsable del 27 de octubre de 2026 (DTC-101): la comparación
//! pareja a pareja **no se repite en cada ejecución**; se garantiza con esta
//! batería obligatoria, que recorre escenarios que combinan las
//! transformaciones autorizadas y juzga cada uno con un oráculo escrito aparte
//! del productor (`final_pairs_oracle`). Lo que se acepta perder: una
//! combinación de configuración que la batería no contemple y falle no se
//! detectaría en la ejecución del usuario.
//!
//! Cada escenario ejecuta el **binario** —el mismo recorrido que usa el
//! usuario— sobre una configuración escrita aquí, y exige además que el
//! oráculo haya visto de verdad la transformación que el escenario dice
//! cubrir: un escenario que no la ejercitara no podría fallar por ella.

mod final_pairs_oracle;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use final_pairs_oracle::{OracleTally, check_final_pairs_against_truth};

#[path = "support/detection_section.rs"]
mod detection_section;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "6a1f0c3e9b2d47586e7f8091a2b3c4d5e6f708192a3b4c5d6e7f8091a2b3c4d5";

/// Una secuencia de referencia del escenario.
struct Sequence {
    id: &'static str,
    bases: Vec<u8>,
}

/// Escribe la referencia, la matriz de SNV y la configuración, ejecuta el
/// binario y juzga lo publicado con el oráculo.
struct Scenario {
    name: &'static str,
    sequences: Vec<Sequence>,
    /// `(secuencia, posición en base 1, población)`; la alternativa es la
    /// primera base de `ACGT` distinta de la de referencia.
    snvs: Vec<(&'static str, usize, &'static str)>,
    /// Tablas TOML adicionales del documento general: poblaciones, loci y
    /// cambios directos.
    body: String,
    /// Tabla `[technical_event]` completa.
    event: String,
}

impl Scenario {
    fn run(self) -> OracleTally {
        let root = std::env::temp_dir().join(format!(
            "gps-final-pairs-{}-{}-{}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed),
            self.name
        ));
        fs::create_dir_all(&root).expect("raíz del escenario");
        // Se borra al soltarse, también si una aserción falla a mitad.
        let _cleanup = TemporaryRoot(root.clone());
        let mut fasta = String::new();
        let mut declared = String::new();
        for sequence in &self.sequences {
            fasta.push('>');
            fasta.push_str(sequence.id);
            fasta.push('\n');
            for chunk in sequence.bases.chunks(70) {
                fasta.push_str(std::str::from_utf8(chunk).expect("ASCII"));
                fasta.push('\n');
            }
            declared.push_str(&format!(
                "[[reference.sequences]]\nid = \"{}\"\nclass = \"chromosome\"\nlength = {}\n\n",
                sequence.id,
                sequence.bases.len()
            ));
        }
        fs::write(root.join("reference.fa"), fasta).expect("FASTA");
        let mut matrix = String::from("contig\tposition\tref\talt\tpopulation_id\talt_af\n");
        for (contig, position, population) in &self.snvs {
            let bases = &self
                .sequences
                .iter()
                .find(|sequence| sequence.id == *contig)
                .expect("secuencia del SNV")
                .bases;
            let reference = bases[position - 1];
            let alternate = *b"ACGT"
                .iter()
                .find(|base| **base != reference)
                .expect("alternativa");
            matrix.push_str(&format!(
                "{contig}\t{position}\t{}\t{}\t{population}\t0.5\n",
                char::from(reference),
                char::from(alternate)
            ));
        }
        fs::write(root.join("population-snvs.tsv"), matrix).expect("matriz de SNV");
        let config = root.join("general.toml");
        fs::write(
            &config,
            detection_section::with_detection_section(&format!(
                "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
                 [population_snvs]\nmode = \"supplied\"\ninput_path = \"population-snvs.tsv\"\n\n\
                 [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
                 {declared}{}\n{}",
                self.body, self.event
            )),
        )
        .expect("configuración");
        let destination = root.join("published");
        let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
            .arg("run")
            .arg(&config)
            .arg("--output")
            .arg(&destination)
            .output()
            .expect("binario del simulador");
        assert!(
            output.status.success(),
            "{}: la ejecución falló\n{}\n{}",
            self.name,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        check_final_pairs_against_truth(&destination, Path::new(env!("CARGO_MANIFEST_DIR")))
    }
}

/// Directorio temporal de un escenario: se borra al soltarse, de modo que un
/// escenario que falla no deja restos.
struct TemporaryRoot(PathBuf);

impl Drop for TemporaryRoot {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Bases pseudoaleatorias reproducibles (xorshift64), sin depender del
/// programa.
fn bases(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[usize::try_from(state >> 62).expect("dos bits")]
        })
        .collect()
}

/// Sustituye `bases[start..]` (base 0) por `content`.
fn plant(bases: &mut [u8], start: usize, content: &[u8]) {
    bases[start..start + content.len()].copy_from_slice(content);
}

fn library(individual: &str, i7: &str, i5: &str) -> String {
    format!(
        "[[technical_event.libraries]]\nindividual_id = \"{individual}\"\ni7 = \"{i7}\"\ni5 = \"{i5}\"\n"
    )
}

const UDIS: [(&str, &str); 6] = [
    ("ACGTACGT", "TGCATGCA"),
    ("CGTACGTA", "GCATGCAT"),
    ("GTACGTAC", "CATGCATG"),
    ("TACGTACG", "ATGCATGC"),
    ("AACCGGTT", "TTGGCCAA"),
    ("CCAAGGTT", "GGTTCCAA"),
];

fn libraries(individuals: &[&str]) -> String {
    individuals
        .iter()
        .zip(UDIS)
        .map(|(individual, (i7, i5))| library(individual, i7, i5))
        .collect()
}

fn quoted(individuals: &[&str]) -> String {
    individuals
        .iter()
        .map(|individual| format!("\"{individual}\""))
        .collect::<Vec<_>>()
        .join(", ")
}

/// **Haploide, sin errores ni artefactos, sin PCR, con una tanda de `N`.**
///
/// Aísla el muestreo del fragmento (inserto, copia, comienzo), la orientación
/// FR y la curva de calidad, que con los errores apagados es la del perfil
/// literal. La tanda de `N` en la referencia ejercita su calidad cero.
#[test]
fn haploid_reads_without_errors_match_the_truth() {
    let mut chr1 = bases(0x9E37_79B9_7F4A_7C15, 4000);
    plant(&mut chr1, 2000, &[b'N'; 30]);
    let individuals = ["H1", "H2", "H3"];
    let tally = Scenario {
        name: "haploid",
        sequences: vec![Sequence {
            id: "chr1",
            bases: chr1,
        }],
        snvs: vec![
            ("chr1", 300, "P1"),
            ("chr1", 1500, "P1"),
            ("chr1", 3100, "P1"),
        ],
        body: format!(
            "[[initial_populations]]\npopulation_id = \"P1\"\ngenotype_model = \"hwe\"\n\
             individual_count = 3\nploidy = 1\nindividual_ids = [{}]\n",
            quoted(&individuals)
        ),
        event: format!(
            "[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [{}]\n\
             target_coverage_wgs = 8\nsequencing_errors = \"disabled\"\n{}",
            quoted(&individuals),
            libraries(&individuals)
        ),
    }
    .run();
    assert!(tally.final_pairs > 200, "{tally:?}");
    assert!(tally.n_bases > 0, "la tanda de N no se muestreó: {tally:?}");
    assert!(tally.snv_bases_sampled > 0, "{tally:?}");
}

/// Microsatélite `(CA)n` plantado en la referencia en `start` (base 0) y su
/// declaración con alelos generados por deslizamiento.
fn generated_locus(id: &str, sequence: &str, start: usize, units: usize) -> String {
    format!(
        "[[reference.repeat_loci]]\nid = \"{id}\"\nsequence = \"{sequence}\"\nstart = {}\nend = {}\n\
         motif = \"CA\"\nmax_units = {}\n\n[reference.repeat_loci.slippage]\nmodel = \"stepwise\"\n\
         mean_events = 3.0\nexpansion_bias = 0.5\n\n",
        start + 1,
        start + 2 * units,
        units + 15
    )
}

const STUTTER_STEPS: &str = "[[technical_event.stutter.steps]]\nunits = -3\nweight = 0.05\n\
     [[technical_event.stutter.steps]]\nunits = -2\nweight = 0.15\n\
     [[technical_event.stutter.steps]]\nunits = -1\nweight = 0.6\n\
     [[technical_event.stutter.steps]]\nunits = 1\nweight = 0.2\n";

/// **Escenario saturado: todo encendido a la vez.**
///
/// Diploide con SNV, una inserción y una deleción explícitas, dos
/// microsatélites con alelos generados, errores de secuenciación, quimeras,
/// contaminación biológica entre tres individuas, salto de índice entre tres
/// bibliotecas, preparación con PCR con curva de composición no neutra y dos
/// clases de eficiencia, tartamudeo de base y tartamudeo clonal por linaje.
#[test]
fn saturated_scenario_matches_the_truth() {
    let mut chr1 = bases(0xD1B5_4A32_D192_ED03, 5000);
    plant(&mut chr1, 1000, &b"CA".repeat(25));
    plant(&mut chr1, 3500, &b"CA".repeat(18));
    let individuals = ["S1", "S2", "S3"];
    let tally = Scenario {
        name: "saturated",
        sequences: vec![Sequence { id: "chr1", bases: chr1 }],
        snvs: vec![
            ("chr1", 200, "P1"),
            ("chr1", 900, "P1"),
            ("chr1", 1700, "P1"),
            ("chr1", 2500, "P1"),
            ("chr1", 4300, "P1"),
        ],
        body: format!(
            "{}{}[[initial_populations]]\npopulation_id = \"P1\"\ngenotype_model = \"hwe\"\n\
             individual_count = 3\nploidy = 2\nindividual_ids = [{}]\n\n\
             [[direct_events]]\nevent_id = \"ins-1\"\nindividual_id = \"S1\"\nphase = 0\n\
             kind = \"insertion\"\ncontig = \"chr1\"\npoint = 2200\nbases = \"TTGTACG\"\n\n\
             [[direct_events]]\nevent_id = \"del-1\"\nindividual_id = \"S2\"\nphase = 1\n\
             kind = \"deletion\"\ncontig = \"chr1\"\nstart = 2800\nend = 2812\n",
            generated_locus("STR-A", "chr1", 1000, 25),
            generated_locus("STR-B", "chr1", 3500, 18),
            quoted(&individuals)
        ),
        event: format!(
            "[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [{}]\n\
             target_coverage_wgs = 10\nchimera_rate = 0.3\nindex_hopping_rate = 0.2\n\
             contamination_min_percent = 10\ncontamination_max_percent = 30\n\
             library_preparation = \"pcr\"\n\n\
             [technical_event.pcr]\nduplicate_fraction = 0.4\ncycles = 10\nefficiency = 0.9\n\
             [[technical_event.pcr.efficiency_classes]]\nefficiency = 0.95\nweight = 0.6\n\
             [[technical_event.pcr.efficiency_classes]]\nefficiency = 0.7\nweight = 0.4\n\
             [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 0.0\nrelative_efficiency = 0.6\n\
             [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 0.5\nrelative_efficiency = 1.0\n\
             [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 1.0\nrelative_efficiency = 0.6\n\
             [technical_event.pcr.stutter]\n[[technical_event.pcr.stutter.per_cycle_slip]]\n\
             motif_length_bases = 2\nprobability = 0.05\n\n\
             [technical_event.stutter]\nenabled = true\n\
             [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = 0.5\n\
             [[technical_event.stutter.locus_rates]]\nlocus_id = \"STR-B\"\nrate = 0.6\n{STUTTER_STEPS}\n{}",
            quoted(&individuals),
            libraries(&individuals)
        ),
    }
    .run();
    assert!(tally.chimeric_pairs > 0, "{tally:?}");
    assert!(tally.base_stutter_steps > 0, "{tally:?}");
    assert!(tally.lineage_steps > 0, "{tally:?}");
    assert!(tally.sequencing_errors > 0, "{tally:?}");
    assert!(tally.contaminant_pairs > 0, "{tally:?}");
    assert!(tally.hopped_pairs > 0, "{tally:?}");
    assert!(tally.duplicate_pairs > 0, "{tally:?}");
    assert!(tally.gc_retries > 0, "{tally:?}");
    assert!(tally.inserted_bases_sampled > 0, "{tally:?}");
    assert!(tally.repeat_bases_sampled > 0, "{tally:?}");
}

/// **Reproducción mínima de DT-060: el tartamudeo de una molécula quimérica
/// solo cae en sus tractos declarados.**
///
/// Haploide, una individua, un solo microsatélite `(CA)20` con su alelo de
/// referencia, todas las parejas quiméricas (`chimera_rate = 1`) y el locus
/// resbalando siempre (`locus_rates = 1`), sin errores ni PCR. El contrato
/// (REQ-DAT-74 a 76) dice que el tartamudeo se aplica a la molécula quimérica
/// ya formada, a **cada tracto declarado contenido en ella**, y que una región
/// que no sea un tracto declarado **no tartamudea**.
///
/// Hasta el 27 de octubre de 2026 el programa buscaba el tracto en la copia
/// del prefijo y aplicaba su ventana a toda la molécula: más allá del punto de
/// corte expandía `CA` sueltos del segundo origen y no evaluaba nunca un tracto
/// que entrara con él. Esta prueba fallaba entonces con 10 parejas distintas;
/// vuelve a fallar si el defecto se reintroduce.
#[test]
fn chimeric_molecules_stutter_only_their_declared_tracts() {
    let mut chr1 = bases(0x2545_F491_4F6C_DD1D, 3000);
    plant(&mut chr1, 1500, &b"CA".repeat(20));
    let tally = Scenario {
        name: "chimera-stutter",
        sequences: vec![Sequence { id: "chr1", bases: chr1 }],
        snvs: vec![("chr1", 400, "P1")],
        body: "[[reference.repeat_loci]]\nid = \"STR-1\"\nsequence = \"chr1\"\nstart = 1501\nend = 1540\n\
               motif = \"CA\"\n\n[[reference.repeat_loci.alleles]]\nunits = 20\nfrequency = 1.0\n\n\
               [[initial_populations]]\npopulation_id = \"P1\"\ngenotype_model = \"hwe\"\n\
               individual_count = 1\nploidy = 1\nindividual_ids = [\"C1\"]\n"
            .to_owned(),
        event: format!(
            "[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [\"C1\"]\n\
             target_coverage_wgs = 10\nsequencing_errors = \"disabled\"\nchimera_rate = 1.0\n\n\
             [technical_event.stutter]\nenabled = true\n\
             [[technical_event.stutter.locus_rates]]\nlocus_id = \"STR-1\"\nrate = 1.0\n{STUTTER_STEPS}\n{}",
            libraries(&["C1"])
        ),
    }
    .run();
    assert!(
        tally.chimeric_pairs > 0 && tally.base_stutter_steps > 0,
        "{tally:?}"
    );
}

fn population(id: &str, ploidy: u8, individuals: &[&str]) -> String {
    format!(
        "[[initial_populations]]\npopulation_id = \"{id}\"\ngenotype_model = \"hwe\"\n\
         individual_count = {}\nploidy = {ploidy}\nindividual_ids = [{}]\n\n",
        individuals.len(),
        quoted(individuals)
    )
}

fn event_header(individuals: &[&str], coverage: u32, extra: &str) -> String {
    format!(
        "[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [{}]\n\
         target_coverage_wgs = {coverage}\n{extra}",
        quoted(individuals)
    )
}

/// **Diploide con SNV, una inserción y una deleción explícitas y errores de
/// secuenciación activos; sin artefactos ni PCR.**
///
/// Aísla la sustitución de errores (Bernoulli por ciclo y base alternativa) y
/// su calidad, sobre copias cuya longitud cambia con los *indels*. Una tanda
/// de `N` comprueba que, con los errores activos, una `N` sigue sin error y
/// con calidad cero.
#[test]
fn diploid_reads_with_indels_and_sequencing_errors_match_the_truth() {
    fn diploid_reference() -> Vec<u8> {
        let mut chr1 = bases(0xA076_1D64_78BD_642F, 4000);
        plant(&mut chr1, 3600, &[b'N'; 25]);
        chr1
    }
    let individuals = ["D1", "D2"];
    let tally = Scenario {
        name: "diploid-errors",
        sequences: vec![Sequence {
            id: "chr1",
            bases: diploid_reference(),
        }],
        snvs: vec![
            ("chr1", 250, "P1"),
            ("chr1", 1900, "P1"),
            ("chr1", 3300, "P1"),
        ],
        body: format!(
            "{}[[direct_events]]\nevent_id = \"ins-1\"\nindividual_id = \"D1\"\nphase = 1\n\
             kind = \"insertion\"\ncontig = \"chr1\"\npoint = 1200\nbases = \"GATTACAGATTACA\"\n\n\
             [[direct_events]]\nevent_id = \"del-1\"\nindividual_id = \"D2\"\nphase = 0\n\
             kind = \"deletion\"\ncontig = \"chr1\"\nstart = 2600\nend = 2640\n",
            population("P1", 2, &individuals)
        ),
        event: event_header(&individuals, 10, &libraries(&individuals)),
    }
    .run();
    assert!(tally.sequencing_errors > 20, "{tally:?}");
    assert!(tally.n_bases > 0, "la tanda de N no se muestreó: {tally:?}");
    assert!(tally.inserted_bases_sampled > 0, "{tally:?}");
    assert_eq!(tally.copy_indices_sampled.len(), 2, "{tally:?}");
}

/// **Tetraploide con un microsatélite de alelos generados, sin tartamudeo.**
///
/// Las cuatro copias llevan alelos de longitud distinta (la verdad cambia de
/// longitud), de modo que un error en la elección de copia o en el reparto de
/// comienzos entre copias de longitud desigual se ve en las bases.
#[test]
fn tetraploid_reads_over_generated_microsatellite_alleles_match_the_truth() {
    let mut chr1 = bases(0xE703_7ED1_A0B4_28DB, 3500);
    plant(&mut chr1, 1600, &b"CA".repeat(20));
    let individuals = ["T1", "T2"];
    let tally = Scenario {
        name: "tetraploid",
        sequences: vec![Sequence {
            id: "chr1",
            bases: chr1,
        }],
        snvs: vec![("chr1", 500, "P1"), ("chr1", 2900, "P1")],
        body: format!(
            "{}{}",
            generated_locus("STR-T", "chr1", 1600, 20),
            population("P1", 4, &individuals)
        ),
        event: event_header(
            &individuals,
            8,
            &format!(
                "sequencing_errors = \"disabled\"\n{}",
                libraries(&individuals)
            ),
        ),
    }
    .run();
    assert_eq!(tally.copy_indices_sampled.len(), 4, "{tally:?}");
    assert!(tally.repeat_bases_sampled > 0, "{tally:?}");
}

/// **Dos secuencias declaradas en orden distinto del alfabético, con quimeras.**
///
/// `chrB` se declara antes que `chrA`. El fragmento recorre las copias en el
/// orden de la referencia, y el segundo origen de la quimera ordena sus
/// candidatas como fija DTC-084 —corregida el 27 de octubre de 2026—: por
/// identificador de secuencia en orden de bytes de su nombre, después por
/// índice de copia y después por inicio. Aquí ese orden no coincide con el de
/// la referencia, y con dos copias por secuencia tampoco con uno que pusiera
/// el índice de copia primero: el escenario distingue los tres (comprobado
/// cambiando el orden en el oráculo: 120 parejas distintas con el índice de
/// copia primero, 239 con el orden de la referencia). Es la regresión de ese
/// orden.
#[test]
fn chimeras_across_two_sequences_match_the_truth() {
    let individuals = ["Q1", "Q2"];
    let tally = Scenario {
        name: "two-sequences",
        sequences: vec![
            Sequence {
                id: "chrB",
                bases: bases(0x8EBC_6AF0_9C88_C6E3, 2600),
            },
            Sequence {
                id: "chrA",
                bases: bases(0x5899_65CC_7537_4CC3, 1800),
            },
        ],
        snvs: vec![("chrB", 700, "P1"), ("chrA", 900, "P1")],
        body: population("P1", 2, &individuals),
        event: event_header(
            &individuals,
            10,
            &format!(
                "sequencing_errors = \"disabled\"\nchimera_rate = 0.6\n{}",
                libraries(&individuals)
            ),
        ),
    }
    .run();
    assert_eq!(tally.sequences_sampled.len(), 2, "{tally:?}");
    assert!(tally.chimeras_across_sequences > 20, "{tally:?}");
}

/// **Contaminación biológica con salto de índice, en dos poblaciones.**
///
/// Cuatro individuas en dos poblaciones: una receptora solo puede recibir de
/// su propia población. Receptoras declaradas explícitamente (no todas), salto
/// de índice entre las cuatro bibliotecas y errores activos, de modo que una
/// contaminante lleva los errores de su biblioteca de origen.
#[test]
fn contamination_with_index_hopping_matches_the_truth() {
    let individuals = ["K1", "K2", "K3", "K4"];
    let tally = Scenario {
        name: "contamination",
        sequences: vec![Sequence {
            id: "chr1",
            bases: bases(0x1D8E_4E27_C47D_124F, 3000),
        }],
        snvs: vec![
            ("chr1", 600, "PA"),
            ("chr1", 600, "PB"),
            ("chr1", 2100, "PA"),
            ("chr1", 2100, "PB"),
        ],
        body: format!(
            "{}{}",
            population("PA", 2, &["K1", "K2"]),
            population("PB", 2, &["K3", "K4"])
        ),
        event: event_header(
            &individuals,
            10,
            &format!(
                "contamination_min_percent = 20\ncontamination_max_percent = 40\n\
                 contamination_recipient_ids = [\"K1\", \"K3\", \"K4\"]\nindex_hopping_rate = 0.3\n{}",
                libraries(&individuals)
            ),
        ),
    }
    .run();
    assert!(tally.contaminant_pairs > 30, "{tally:?}");
    assert!(tally.hopped_pairs > 30, "{tally:?}");
    assert!(tally.sequencing_errors > 0, "{tally:?}");
}

/// **Preparación con PCR sola: familias, duplicados y aceptación por
/// composición, sin microsatélites ni artefactos.**
///
/// Dos clases de eficiencia y una curva de composición que penaliza los
/// extremos, de modo que hay reintentos del fragmento con sus claves propias.
#[test]
fn pcr_families_without_microsatellites_match_the_truth() {
    let mut chr1 = bases(0x2F2E_1A36_B5C9_AA4B, 4000);
    // Una región rica en AT y otra rica en GC, para que la curva rechace.
    plant(&mut chr1, 500, &b"AATTATATTAAT".repeat(40));
    plant(&mut chr1, 2600, &b"GCGGCCGCGCGG".repeat(40));
    let individuals = ["R1", "R2"];
    let tally = Scenario {
        name: "pcr",
        sequences: vec![Sequence { id: "chr1", bases: chr1 }],
        snvs: vec![("chr1", 1500, "P1")],
        body: population("P1", 2, &individuals),
        event: event_header(
            &individuals,
            10,
            &format!(
                "sequencing_errors = \"disabled\"\nlibrary_preparation = \"pcr\"\n\n\
                 [technical_event.pcr]\nduplicate_fraction = 0.5\ncycles = 12\nefficiency = 0.9\n\
                 early_cycles = 2\nmax_family_size = 6\ngc_max_attempts = 4\n\
                 [[technical_event.pcr.efficiency_classes]]\nefficiency = 0.95\nweight = 0.5\n\
                 [[technical_event.pcr.efficiency_classes]]\nefficiency = 0.6\nweight = 0.5\n\
                 [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 0.0\nrelative_efficiency = 0.4\n\
                 [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 0.5\nrelative_efficiency = 1.0\n\
                 [[technical_event.pcr.gc_relative_efficiency]]\ngc_fraction = 1.0\nrelative_efficiency = 0.4\n\n{}",
                libraries(&individuals)
            ),
        ),
    }
    .run();
    assert!(tally.duplicate_pairs > 50, "{tally:?}");
    assert!(tally.gc_retries > 0, "{tally:?}");
}

/// Tracto impuro con dos tiradas de seis `CA` separadas por `GG`: la tirada
/// más larga está empatada y el contrato elige la de la izquierda
/// (REQ-DAT-77).
const IMPURE_TRACT: &[u8] = b"CACACACACACAGGCACACACACACA";

/// **Tartamudeo de base solo, sin PCR ni quimeras.**
///
/// Un dinucleótido con alelos generados, un trinucleótido con dos alelos
/// declarados, un mononucleótido y un tracto impuro con dos tiradas
/// empatadas, con tasas por longitud de motivo, sustituciones por locus y
/// coeficiente de tirada positivo. El tracto impuro ejercita el desempate por
/// la izquierda: contraer la tirada de la derecha da otras bases.
#[test]
fn base_stutter_without_pcr_matches_the_truth() {
    let mut chr1 = bases(0x7A5D_3C1B_9E4F_2A61, 4000);
    plant(&mut chr1, 800, &b"CA".repeat(22));
    plant(&mut chr1, 2000, &b"GAT".repeat(14));
    plant(&mut chr1, 3000, &[b'A'; 16]);
    plant(&mut chr1, 3500, IMPURE_TRACT);
    let individuals = ["U1", "U2"];
    let tally = Scenario {
        name: "stutter",
        sequences: vec![Sequence { id: "chr1", bases: chr1 }],
        snvs: vec![("chr1", 400, "P1")],
        body: format!(
            "{}[[reference.repeat_loci]]\nid = \"STR-TRI\"\nsequence = \"chr1\"\nstart = 2001\nend = 2042\n\
             motif = \"GAT\"\n\n[[reference.repeat_loci.alleles]]\nunits = 14\nfrequency = 0.5\n\n\
             [[reference.repeat_loci.alleles]]\nunits = 12\nfrequency = 0.5\n\n\
             [[reference.repeat_loci]]\nid = \"STR-MONO\"\nsequence = \"chr1\"\nstart = 3001\nend = 3016\n\
             motif = \"A\"\n\n[[reference.repeat_loci.alleles]]\nunits = 16\nfrequency = 1.0\n\n\
             [[reference.repeat_loci]]\nid = \"STR-IMPURE\"\nsequence = \"chr1\"\nstart = 3501\nend = {}\n\
             motif = \"CA\"\n\n[[reference.repeat_loci.alleles]]\nsequence = \"{}\"\nfrequency = 1.0\n\n{}",
            generated_locus("STR-DI", "chr1", 800, 22),
            3500 + IMPURE_TRACT.len(),
            std::str::from_utf8(IMPURE_TRACT).expect("ASCII"),
            population("P1", 2, &individuals)
        ),
        event: event_header(
            &individuals,
            12,
            &format!(
                "sequencing_errors = \"disabled\"\n\n[technical_event.stutter]\nenabled = true\n\
                 run_length_coefficient = 0.01\n\
                 [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 1\nrate = 0.4\n\
                 [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = 0.3\n\
                 [[technical_event.stutter.locus_rates]]\nlocus_id = \"STR-TRI\"\nrate = 0.5\n\
                 [[technical_event.stutter.locus_rates]]\nlocus_id = \"STR-IMPURE\"\nrate = 1.0\n{STUTTER_STEPS}\n{}",
                libraries(&individuals)
            ),
        ),
    }
    .run();
    assert!(tally.base_stutter_steps > 30, "{tally:?}");
    assert!(tally.stutter_steps_on_tied_runs > 0, "{tally:?}");
}

/// **Tartamudeo clonal por linaje de PCR, sin tartamudeo de base efectivo y
/// sin quimeras.**
///
/// Las tasas de base son cero (la decisión de base se evalúa y nunca resbala),
/// y la probabilidad por ciclo del dinucleótido se declara alta, de modo que
/// lo que cambia las bases es solo la genealogía de cada familia: progenitor,
/// ciclo de divergencia, resbalón por ciclo y paso.
#[test]
fn clonal_stutter_lineage_matches_the_truth() {
    let mut chr1 = bases(0x3C6E_F372_FE94_F82B, 3500);
    plant(&mut chr1, 1500, &b"CA".repeat(20));
    let individuals = ["L1", "L2"];
    let tally = Scenario {
        name: "clonal",
        sequences: vec![Sequence { id: "chr1", bases: chr1 }],
        snvs: vec![("chr1", 700, "P1")],
        body: format!(
            "{}{}",
            generated_locus("STR-L", "chr1", 1500, 20),
            population("P1", 2, &individuals)
        ),
        event: event_header(
            &individuals,
            20,
            &format!(
                "sequencing_errors = \"disabled\"\nlibrary_preparation = \"pcr\"\n\n\
                 [technical_event.pcr]\nduplicate_fraction = 0.6\ncycles = 10\nefficiency = 0.9\n\
                 [technical_event.pcr.stutter]\n[[technical_event.pcr.stutter.per_cycle_slip]]\n\
                 motif_length_bases = 2\nprobability = 0.2\n\n\
                 [technical_event.stutter]\nenabled = true\n\
                 [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = 0.0\n{STUTTER_STEPS}\n{}",
                libraries(&individuals)
            ),
        ),
    }
    .run();
    assert!(tally.lineage_steps > 30, "{tally:?}");
    assert!(tally.duplicate_pairs > 50, "{tally:?}");
}

/// **Tartamudeo clonal con quimeras: el paso del linaje solo cae en el
/// primer tramo.**
///
/// Diploide con el mismo locus en las dos copias y casi todas las parejas
/// quiméricas, de modo que hay moléculas cuyo sufijo trae el locus desde la
/// otra copia. El paso heredado pertenece a la molécula original de la familia
/// y no se aplica al tracto del sufijo (DT-060); el escenario exige que ese
/// caso haya ocurrido.
#[test]
fn clonal_stutter_with_chimeras_matches_the_truth() {
    let mut chr1 = bases(0x4B4D_5F3A_8C1E_6D27, 2400);
    plant(&mut chr1, 1100, &b"CA".repeat(20));
    let individuals = ["M1"];
    let tally = Scenario {
        name: "clonal-chimera",
        sequences: vec![Sequence { id: "chr1", bases: chr1 }],
        snvs: vec![("chr1", 300, "P1")],
        body: format!(
            "{}{}",
            generated_locus("STR-M", "chr1", 1100, 20),
            population("P1", 2, &individuals)
        ),
        event: event_header(
            &individuals,
            40,
            &format!(
                "sequencing_errors = \"disabled\"\nchimera_rate = 0.9\nlibrary_preparation = \"pcr\"\n\n\
                 [technical_event.pcr]\nduplicate_fraction = 0.6\ncycles = 10\nefficiency = 0.9\n\
                 [technical_event.pcr.stutter]\n[[technical_event.pcr.stutter.per_cycle_slip]]\n\
                 motif_length_bases = 2\nprobability = 0.1\n\n\
                 [technical_event.stutter]\nenabled = true\n\
                 [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = 0.2\n{STUTTER_STEPS}\n{}",
                libraries(&individuals)
            ),
        ),
    }
    .run();
    assert!(tally.lineage_steps_withheld_in_suffix > 0, "{tally:?}");
    assert!(
        tally.lineage_steps > 0 && tally.base_stutter_steps > 0,
        "{tally:?}"
    );
}

/// **Quimeras entre dos secuencias con un microsatélite en cada una.**
///
/// `chrB` lleva `(CA)20` y `chrA` lleva `(GAT)12`; todas las parejas son
/// quiméricas y los dos loci resbalan siempre. Una molécula cuyo prefijo sale
/// de una secuencia y su sufijo de la otra lleva en el sufijo un tracto que
/// solo existe en la segunda: el tartamudeo tiene que buscarlo con la copia
/// del sufijo, no con la del prefijo (DT-060). El escenario exige que hayan
/// tartamudeado tractos así. Lo pidió la revisión del 27/10/2026: filtrar los
/// loci por la secuencia del primer tramo pasaba toda la batería anterior.
#[test]
fn chimeric_stutter_across_two_sequences_matches_the_truth() {
    let mut chr_b = bases(0x8EBC_6AF0_9C88_C6E3, 1400);
    plant(&mut chr_b, 600, &b"CA".repeat(20));
    let mut chr_a = bases(0x5899_65CC_7537_4CC3, 1200);
    plant(&mut chr_a, 500, &b"GAT".repeat(12));
    let tally = Scenario {
        name: "chimera-stutter-two-sequences",
        sequences: vec![
            Sequence {
                id: "chrB",
                bases: chr_b,
            },
            Sequence {
                id: "chrA",
                bases: chr_a,
            },
        ],
        snvs: vec![("chrB", 100, "P1")],
        body: format!(
            "[[reference.repeat_loci]]\nid = \"STR-B\"\nsequence = \"chrB\"\nstart = 601\nend = 640\n\
             motif = \"CA\"\n\n[[reference.repeat_loci.alleles]]\nunits = 20\nfrequency = 0.5\n\n\
             [[reference.repeat_loci.alleles]]\nunits = 17\nfrequency = 0.5\n\n\
             [[reference.repeat_loci]]\nid = \"STR-A\"\nsequence = \"chrA\"\nstart = 501\nend = 536\n\
             motif = \"GAT\"\n\n[[reference.repeat_loci.alleles]]\nunits = 12\nfrequency = 1.0\n\n{}",
            population("P1", 2, &["X1"])
        ),
        event: event_header(
            &["X1"],
            60,
            &format!(
                "sequencing_errors = \"disabled\"\nchimera_rate = 1.0\n\n[technical_event.stutter]\nenabled = true\n\
                 [[technical_event.stutter.locus_rates]]\nlocus_id = \"STR-B\"\nrate = 1.0\n\
                 [[technical_event.stutter.locus_rates]]\nlocus_id = \"STR-A\"\nrate = 1.0\n{STUTTER_STEPS}\n{}",
                libraries(&["X1"])
            ),
        ),
    }
    .run();
    assert!(tally.chimeras_across_sequences > 0, "{tally:?}");
    assert!(
        tally.suffix_stutter_steps_from_other_sequence > 0,
        "ningún tracto del sufijo de otra secuencia tartamudeó: {tally:?}"
    );
}
