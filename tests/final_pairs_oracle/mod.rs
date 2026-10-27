//! Oráculo independiente de las parejas finales del FASTQ (REQ-VAL-03).
//!
//! Lee **solo lo publicado** —la configuración efectiva, la referencia, la
//! verdad (`truth.jsonl`), el perfil técnico incorporado y los FASTQ finales—
//! y rehace cada pareja esperada desde su clave técnica, la verdad biológica y
//! las transformaciones autorizadas. Después la compara letra a letra con lo
//! que el programa escribió: bases, calidades, mate, cabecera, biblioteca tras
//! el salto de índice y bandera de contaminación agregada.
//!
//! # Qué comparte con el productor y qué no
//!
//! Comparte **únicamente** las primitivas normativas documentadas:
//!
//! - la derivación de subsemillas del catálogo técnico F4 (DTC-086 y sus
//!   ampliaciones) a través de `derive_f4_technical_subseed` y de los
//!   constructores de `F4TechnicalDecision`, que son la definición del
//!   catálogo;
//! - las primitivas de muestreo de DTC-088 (`f4_uniform_index`, `f4_bernoulli`,
//!   `f4_uniform_closed`, `f4_truncated_normal`, `f4_partial_fisher_yates`) y
//!   `categorical_f64` (DTC-072);
//! - la derivación de los identificadores de biblioteca del catálogo de
//!   identidades (`derive_technical_identity_plan`), cuyo resultado se exige
//!   además igual al publicado en la configuración efectiva.
//!
//! Todo lo demás está **reescrito aquí** sin llamar al productor: la
//! reconstrucción de cada copia haplotípica desde la verdad, el número de
//! parejas, el sorteo del inserto y su colocación, la familia de
//! amplificación y su aceptación por composición, la genealogía clonal del
//! tartamudeo, la quimera, el tartamudeo de base, la orientación FR, la
//! calidad y los errores, el plan de contaminación, el demultiplexado con
//! salto de índice y las identidades `PH_`/`TP_` (BLAKE3 de DTC-087, escrito de
//! nuevo con el crate `blake3`). Ninguna de esas funciones del programa se
//! invoca para calcular lo esperado: si una cambiara, el oráculo no cambiaría
//! con ella.
//!
//! Los parámetros numéricos que el programa **publica** (tasas, pesos, tablas
//! de la dependencia con la longitud del tracto, fracciones de la PCR, media de
//! parejas por molécula) se leen de la configuración efectiva publicada y no se
//! copian como constantes: así el oráculo juzga contra lo que el programa
//! declara haber usado.
//!
//! # Lo que no compara
//!
//! El **orden** de las parejas dentro de cada FASTQ no es contrato: se compara
//! el conjunto exacto de parejas de cada biblioteca final (ninguna falta,
//! ninguna sobra, ninguna repetida) y el contenido de cada una.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::num::NonZeroU64;
use std::path::{Path, PathBuf};
use std::process::Command;

use genome_population_simulator::identity::{
    EventId, IndividualId, LibraryId, ReferenceId, RepeatLocusId,
};
use genome_population_simulator::randomness::{
    F4Mate, F4TechnicalDecision, FragmentRetryDecision, LocalRng, RootSeed, categorical_f64,
    derive_f4_technical_subseed, f4_bernoulli, f4_partial_fisher_yates, f4_truncated_normal,
    f4_uniform_closed, f4_uniform_index,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use serde_json::Value;
use sha2::{Digest, Sha256};

/// Pasos de tartamudeo admitidos, en el orden en que el contrato publica sus
/// pesos (REQ-CONF-131; la estructura la fija el algoritmo).
const ADMITTED_STEPS: [i64; 4] = [-3, -2, -1, 1];

/// Recuentos de lo que el oráculo vio pasar, para que cada escenario pueda
/// exigir que de verdad ejercitó la transformación que dice cubrir.
#[derive(Clone, Debug, Default)]
pub struct OracleTally {
    pub final_pairs: u64,
    pub libraries: u64,
    pub chimeric_pairs: u64,
    pub base_stutter_steps: u64,
    pub lineage_steps: u64,
    pub sequencing_errors: u64,
    pub contaminant_pairs: u64,
    pub hopped_pairs: u64,
    pub duplicate_pairs: u64,
    pub gc_retries: u64,
    pub n_bases: u64,
    pub inserted_bases_sampled: u64,
    pub repeat_bases_sampled: u64,
    pub snv_bases_sampled: u64,
    pub copy_indices_sampled: BTreeSet<u8>,
    /// Tractos del segundo tramo de una quimera cuyo locus tenía paso de
    /// linaje en la familia: el paso no se les aplica.
    pub lineage_steps_withheld_in_suffix: u64,
    /// Quimeras cuyo segundo origen está en otra secuencia de referencia.
    pub chimeras_across_sequences: u64,
    /// Pasos de tartamudeo de base escritos en un tracto del segundo tramo de
    /// una quimera que viene de otra secuencia que el primero.
    pub suffix_stutter_steps_from_other_sequence: u64,
    /// Pasos de tartamudeo de base escritos en una ventana con dos tiradas
    /// empatadas en longitud: ejercitan el desempate por la izquierda.
    pub stutter_steps_on_tied_runs: u64,
    /// Moléculas quiméricas con el mismo locus en los dos tramos: sus dos
    /// tractos comparten la clave de la decisión de base (DTC-102).
    pub same_locus_in_both_segments: u64,
    pub sequences_sampled: BTreeSet<String>,
}

/// Compara el paquete publicado en `published` con la reconstrucción
/// independiente y devuelve los recuentos. Falla con un informe legible en
/// cuanto una sola pareja difiere.
pub fn check_final_pairs_against_truth(published: &Path, repository: &Path) -> OracleTally {
    let run = PublishedRun::load(published, repository);
    let mut tally = OracleTally::default();
    let expected = run.expected_final_pairs(&mut tally);
    let mut mismatches = Vec::new();
    for library in run.event.libraries.values() {
        tally.libraries += 1;
        let observed = run.read_final_library(library);
        let wanted = expected
            .get(&library.library_id)
            .cloned()
            .unwrap_or_default();
        compare_library(library, &wanted, &observed, &mut mismatches);
        // El resumen agregado del manifiesto es la única huella publicada de
        // qué parejas finales son contaminantes; no está en el FASTQ.
        let contaminant = wanted.values().filter(|pair| pair.is_contaminant).count() as u64;
        let own = wanted.len() as u64 - contaminant;
        let summary = run.final_summary(&library.individual_id);
        if summary.0 != own || summary.1 != contaminant {
            mismatches.push(format!(
                "biblioteca {}: el manifiesto publica {} propias y {} contaminantes; el oráculo espera {own} y {contaminant}",
                library.individual_id, summary.0, summary.1
            ));
        }
        tally.final_pairs += observed.len() as u64;
    }
    assert!(
        mismatches.is_empty(),
        "{} discrepancias entre el FASTQ final y la reconstrucción independiente; primeras:\n{}",
        mismatches.len(),
        mismatches
            .iter()
            .take(8)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    tally
}

// ---------------------------------------------------------------------------
// Lectura de lo publicado.
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Piece {
    offset: u64,
    length: u64,
    kind: PieceKind,
}

#[derive(Clone, Debug)]
enum PieceKind {
    Source { start: u64, end: u64 },
    Inserted,
    Repeat { locus_id: String },
}

/// Una copia haplotípica reconstruida desde `current_pieces` de la verdad.
#[derive(Clone, Debug)]
struct CopyView {
    sequence_id: String,
    copy_index: u8,
    bases: Vec<u8>,
    pieces: Vec<Piece>,
    /// Marca por base: 0 referencia, 1 SNV, 2 insertada, 3 tracto de
    /// repetición. Solo alimenta los recuentos de cobertura del oráculo.
    origin: Vec<u8>,
}

#[derive(Clone, Debug)]
struct IndividualView {
    id: String,
    population: String,
    /// Copias en orden de referencia (orden declarado de las secuencias) y,
    /// dentro de cada secuencia, por índice de copia.
    copies: Vec<CopyView>,
}

#[derive(Clone, Debug)]
struct Locus {
    id: String,
    typed: RepeatLocusId,
    sequence_id: String,
    start: u64,
    end: u64,
    motifs: Vec<Vec<u8>>,
}

#[derive(Clone, Copy, Debug)]
struct QualitySegment {
    first_cycle: u32,
    last_cycle: u32,
    first_phred: i64,
    last_phred: i64,
}

#[derive(Clone, Debug)]
struct Profile {
    insert_mean: f64,
    insert_sd: f64,
    insert_min: u64,
    insert_max: u64,
    quality_r1: Vec<QualitySegment>,
    quality_r2: Vec<QualitySegment>,
}

#[derive(Clone, Debug)]
struct StutterParams {
    enabled: bool,
    motif_rates: [f64; 6],
    coefficient: f64,
    step_weights: [f64; 4],
    locus_rates: BTreeMap<String, f64>,
    /// Tabla publicada de la dependencia con la longitud del tracto, por
    /// longitud de motivo: `(unidades, factor)` en orden creciente.
    factors: BTreeMap<u64, Vec<(u64, f64)>>,
}

#[derive(Clone, Debug)]
struct ClonalParams {
    cycles: u32,
    /// Por longitud de motivo 1..=6: `(declarada, probabilidad, con PCR, sin PCR)`.
    per_motif: [(bool, f64, f64, f64); 6],
}

#[derive(Clone, Debug)]
struct PcrParams {
    cycles: u32,
    classes: Vec<(f64, f64)>,
    early_cycles: u32,
    max_family_size: u64,
    gc_max_attempts: u64,
    gc_curve: Vec<(f64, f64)>,
    mean_pairs: f64,
    clonal: Option<ClonalParams>,
}

#[derive(Clone, Debug)]
struct LibraryView {
    individual_id: String,
    library_id: String,
    typed: LibraryId,
    i7: String,
    i5: String,
}

#[derive(Clone, Debug)]
struct EventView {
    event_id: String,
    typed: EventId,
    libraries: BTreeMap<String, LibraryView>,
    read_length: u64,
    target_coverage: f64,
    chimera_rate: f64,
    hopping_rate: f64,
    contamination_min: f64,
    contamination_max: f64,
    recipients: Vec<String>,
    errors_enabled: bool,
    stutter: StutterParams,
    pcr: Option<PcrParams>,
    profile: Profile,
}

struct PublishedRun {
    root: PathBuf,
    seed: RootSeed,
    reference_length: u64,
    individuals: BTreeMap<String, IndividualView>,
    loci: Vec<Locus>,
    event: EventView,
    event_manifest: Value,
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(
        &std::fs::read(path).unwrap_or_else(|error| panic!("{}: {error}", path.display())),
    )
    .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn number(value: &Value, what: &str) -> f64 {
    value
        .as_f64()
        .unwrap_or_else(|| panic!("{what}: se esperaba un número y hay {value}"))
}

fn integer(value: &Value, what: &str) -> u64 {
    value
        .as_u64()
        .unwrap_or_else(|| panic!("{what}: se esperaba un entero y hay {value}"))
}

fn text<'a>(value: &'a Value, what: &str) -> &'a str {
    value
        .as_str()
        .unwrap_or_else(|| panic!("{what}: se esperaba un texto y hay {value}"))
}

fn read_fasta(path: &Path) -> BTreeMap<String, Vec<u8>> {
    let content =
        std::fs::read_to_string(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let mut sequences = BTreeMap::new();
    let mut current: Option<(String, Vec<u8>)> = None;
    for line in content.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some((id, bases)) = current.take() {
                sequences.insert(id, bases);
            }
            let id = header.split_whitespace().next().unwrap_or("").to_owned();
            current = Some((id, Vec::new()));
        } else if let Some((_, bases)) = current.as_mut() {
            bases.extend(line.trim().bytes().map(|base| base.to_ascii_uppercase()));
        }
    }
    if let Some((id, bases)) = current {
        sequences.insert(id, bases);
    }
    sequences
}

impl PublishedRun {
    fn load(root: &Path, repository: &Path) -> Self {
        let effective = read_json(&root.join("effective-config.json"));
        let seed = RootSeed::parse("seed", text(&effective["seed"]["root"], "seed.root"))
            .expect("semilla raíz publicada");
        let reference_id = text(&effective["reference"]["id"], "reference.id").to_owned();
        let declared_order = effective["reference"]["sequences"]
            .as_array()
            .expect("secuencias de referencia")
            .iter()
            .map(|sequence| text(&sequence["id"], "reference.sequences.id").to_owned())
            .collect::<Vec<_>>();
        let fasta = read_fasta(&root.join(format!("resources/references/{reference_id}.fa")));
        let reference_length = declared_order
            .iter()
            .map(|id| {
                fasta
                    .get(id)
                    .expect("secuencia declarada en el FASTA")
                    .len() as u64
            })
            .sum();

        let mut population_of = BTreeMap::new();
        for population in effective["initial_populations"]
            .as_array()
            .expect("poblaciones iniciales")
        {
            let population_id = text(&population["population_id"], "population_id");
            for individual in population["individual_ids"]
                .as_array()
                .expect("individuos de la población")
            {
                population_of.insert(
                    text(individual, "individual_id").to_owned(),
                    population_id.to_owned(),
                );
            }
        }

        let individuals = read_truth_copies(
            &root.join("outputs/truth.jsonl"),
            &fasta,
            &declared_order,
            &population_of,
        );

        let loci = effective["repeat_loci"]
            .as_array()
            .map(|loci| {
                loci.iter()
                    .map(|locus| {
                        let id = text(&locus["id"], "repeat_loci.id").to_owned();
                        Locus {
                            typed: RepeatLocusId::parse(id.clone(), "oracle")
                                .expect("identificador de locus"),
                            id,
                            sequence_id: text(&locus["sequence_id"], "repeat_loci.sequence_id")
                                .to_owned(),
                            // Coordenadas declaradas en base 1 e inclusivas.
                            start: integer(&locus["start"], "repeat_loci.start") - 1,
                            end: integer(&locus["end"], "repeat_loci.end"),
                            motifs: locus["motifs"]
                                .as_array()
                                .expect("motivos del locus")
                                .iter()
                                .map(|motif| text(motif, "motif").as_bytes().to_vec())
                                .collect(),
                        }
                    })
                    .collect()
            })
            .unwrap_or_default();

        let event = read_event(&effective, seed, &reference_id, repository);
        let event_manifest =
            read_json(&root.join(format!("events/{}/manifest.json", event.event_id)));
        Self {
            root: root.to_path_buf(),
            seed,
            reference_length,
            individuals,
            loci,
            event,
            event_manifest,
        }
    }

    /// Propias y contaminantes que el manifiesto de evento publica para la
    /// biblioteca final de una individua.
    fn final_summary(&self, individual_id: &str) -> (u64, u64) {
        let summary = self.event_manifest["final_fastq_summary"]
            .as_array()
            .expect("resumen FASTQ final")
            .iter()
            .find(|entry| entry["individual_id"] == individual_id)
            .unwrap_or_else(|| panic!("resumen de {individual_id}"));
        (
            integer(&summary["own_pair_count"], "own_pair_count"),
            integer(&summary["contaminant_pair_count"], "contaminant_pair_count"),
        )
    }

    fn read_final_library(&self, library: &LibraryView) -> BTreeMap<String, ObservedRecordPair> {
        let directory = self
            .root
            .join(format!("outputs/fastq/{}", self.event.event_id));
        let r1 = read_gzip_fastq(&directory.join(format!("{}_R1.fastq.gz", library.individual_id)));
        let r2 = read_gzip_fastq(&directory.join(format!("{}_R2.fastq.gz", library.individual_id)));
        assert_eq!(
            r1.len(),
            r2.len(),
            "{}: R1 y R2 tienen distinto número de registros",
            library.individual_id
        );
        let mut pairs = BTreeMap::new();
        for (first, second) in r1.into_iter().zip(r2) {
            let key = first.pair_key();
            assert_eq!(
                key,
                second.pair_key(),
                "{}: R1 y R2 desemparejados ({} frente a {})",
                library.individual_id,
                first.header,
                second.header
            );
            assert!(
                pairs
                    .insert(
                        key.clone(),
                        ObservedRecordPair {
                            r1: first,
                            r2: second,
                        },
                    )
                    .is_none(),
                "{}: pareja repetida {key}",
                library.individual_id
            );
        }
        pairs
    }
}

/// Reconstruye cada copia concatenando sus piezas vigentes: tramos de la
/// referencia con sus SNV aplicados, bases insertadas y alelos de repetición.
fn read_truth_copies(
    path: &Path,
    fasta: &BTreeMap<String, Vec<u8>>,
    declared_order: &[String],
    population_of: &BTreeMap<String, String>,
) -> BTreeMap<String, IndividualView> {
    let content =
        std::fs::read_to_string(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let mut individuals: BTreeMap<String, IndividualView> = BTreeMap::new();
    for line in content.lines() {
        let record: Value = serde_json::from_str(line).expect("línea JSON de la verdad");
        if record["record_type"] != "haplotype_copy" {
            continue;
        }
        let identity = &record["identity"];
        let individual_id = text(&identity["individual_id"], "identity.individual_id").to_owned();
        let sequence_id = text(&identity["sequence_id"], "identity.sequence_id").to_owned();
        let copy_index =
            u8::try_from(integer(&identity["copy_index"], "copy_index")).expect("índice de copia");
        let mut bases = Vec::new();
        let mut origin = Vec::new();
        let mut pieces = Vec::new();
        for piece in record["current_pieces"]
            .as_array()
            .expect("piezas vigentes")
        {
            let offset = bases.len() as u64;
            match text(&piece["kind"], "piece.kind") {
                "source" => {
                    let source = &piece["source"];
                    let source_sequence = text(&source["sequence_id"], "source.sequence_id");
                    let start = integer(&source["start"], "source.start");
                    let end = integer(&source["end"], "source.end");
                    let reference = fasta
                        .get(source_sequence)
                        .expect("secuencia fuente en el FASTA");
                    let mut slice = reference[usize::try_from(start).expect("inicio")
                        ..usize::try_from(end).expect("fin")]
                        .to_vec();
                    let mut marks = vec![0_u8; slice.len()];
                    for snp in piece["snps"].as_array().expect("SNV de la pieza") {
                        let position = integer(&snp["position"], "snp.position");
                        let index = usize::try_from(position - start).expect("posición del SNV");
                        assert_eq!(
                            slice[index],
                            text(&snp["reference"], "snp.reference").as_bytes()[0],
                            "la base de referencia del SNV no coincide con el FASTA publicado"
                        );
                        slice[index] = text(&snp["alternate"], "snp.alternate").as_bytes()[0];
                        marks[index] = 1;
                    }
                    bases.extend_from_slice(&slice);
                    origin.extend_from_slice(&marks);
                    pieces.push(Piece {
                        offset,
                        length: end - start,
                        kind: PieceKind::Source { start, end },
                    });
                }
                "inserted" => {
                    let inserted = text(&piece["bases"], "inserted.bases").as_bytes();
                    bases.extend_from_slice(inserted);
                    origin.extend(std::iter::repeat_n(2_u8, inserted.len()));
                    pieces.push(Piece {
                        offset,
                        length: inserted.len() as u64,
                        kind: PieceKind::Inserted,
                    });
                }
                "repeat_allele" => {
                    let allele = text(&piece["bases"], "repeat_allele.bases").as_bytes();
                    bases.extend_from_slice(allele);
                    origin.extend(std::iter::repeat_n(3_u8, allele.len()));
                    pieces.push(Piece {
                        offset,
                        length: allele.len() as u64,
                        kind: PieceKind::Repeat {
                            locus_id: text(&piece["locus_id"], "repeat_allele.locus_id").to_owned(),
                        },
                    });
                }
                other => panic!("pieza vigente desconocida en la verdad: {other}"),
            }
        }
        let population = population_of
            .get(&individual_id)
            .cloned()
            .unwrap_or_else(|| panic!("población de {individual_id}"));
        individuals
            .entry(individual_id.clone())
            .or_insert_with(|| IndividualView {
                id: individual_id.clone(),
                population,
                copies: Vec::new(),
            })
            .copies
            .push(CopyView {
                sequence_id,
                copy_index,
                bases,
                pieces,
                origin,
            });
    }
    for individual in individuals.values_mut() {
        individual.copies.sort_by_key(|copy| {
            (
                declared_order
                    .iter()
                    .position(|id| *id == copy.sequence_id)
                    .expect("secuencia declarada"),
                copy.copy_index,
            )
        });
    }
    individuals
}

fn read_event(
    effective: &Value,
    seed: RootSeed,
    reference_id: &str,
    repository: &Path,
) -> EventView {
    let event = &effective["technical_event"];
    let event_id = text(&event["event_id"], "event_id").to_owned();
    let typed_event = EventId::parse(event_id.clone(), "oracle").expect("evento");
    let individual_ids = event["individual_ids"]
        .as_array()
        .expect("individuos del evento")
        .iter()
        .map(|value| {
            IndividualId::parse(text(value, "individual_id"), "oracle").expect("individuo")
        })
        .collect::<Vec<_>>();
    let plan = derive_technical_identity_plan(
        seed,
        &ReferenceId::parse(reference_id, "oracle").expect("referencia"),
        typed_event.clone(),
        &individual_ids,
    )
    .expect("identidades del evento");
    let mut libraries = BTreeMap::new();
    for declared in event["libraries"].as_array().expect("bibliotecas") {
        let individual_id = text(&declared["individual_id"], "libraries.individual_id");
        let typed = plan
            .library(&IndividualId::parse(individual_id, "oracle").expect("individuo"))
            .expect("biblioteca derivada")
            .library_id()
            .clone();
        assert_eq!(
            typed.as_str(),
            text(&declared["library_id"], "libraries.library_id"),
            "la biblioteca publicada no es la que deriva el catálogo de identidades"
        );
        libraries.insert(
            individual_id.to_owned(),
            LibraryView {
                individual_id: individual_id.to_owned(),
                library_id: typed.as_str().to_owned(),
                typed,
                i7: text(&declared["i7"], "i7").to_owned(),
                i5: text(&declared["i5"], "i5").to_owned(),
            },
        );
    }

    let profile_view = &event["profile"];
    let resource = repository.join(text(
        &profile_view["resource_path"],
        "profile.resource_path",
    ));
    let resource_bytes =
        std::fs::read(&resource).unwrap_or_else(|error| panic!("{}: {error}", resource.display()));
    assert_eq!(
        hex(&Sha256::digest(&resource_bytes)),
        text(&profile_view["sha256"], "profile.sha256"),
        "el recurso del perfil no es el que la configuración efectiva declara"
    );
    let profile_json: Value = serde_json::from_slice(&resource_bytes).expect("perfil JSON");
    let segments = |mate: &str| {
        profile_json["quality"][mate]
            .as_array()
            .expect("curva de calidad")
            .iter()
            .map(|segment| QualitySegment {
                first_cycle: u32::try_from(integer(&segment["first_cycle"], "first_cycle"))
                    .expect("ciclo"),
                last_cycle: u32::try_from(integer(&segment["last_cycle"], "last_cycle"))
                    .expect("ciclo"),
                first_phred: i64::try_from(integer(&segment["first_phred"], "first_phred"))
                    .expect("Phred"),
                last_phred: i64::try_from(integer(&segment["last_phred"], "last_phred"))
                    .expect("Phred"),
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(
        profile_json["quality"]["curve_kind"],
        "piecewise_linear_phred_integer_round_half_up_v1"
    );
    assert_eq!(profile_json["read"]["orientation"], "FR");
    let profile = Profile {
        insert_mean: number(&profile_json["insert"]["mean_bases"], "mean_bases"),
        insert_sd: number(
            &profile_json["insert"]["standard_deviation_bases"],
            "standard_deviation_bases",
        ),
        insert_min: integer(&profile_json["insert"]["minimum_bases"], "minimum_bases"),
        insert_max: integer(&profile_json["insert"]["maximum_bases"], "maximum_bases"),
        quality_r1: segments("r1"),
        quality_r2: segments("r2"),
    };

    let artifacts = &event["artifacts"];
    let contamination = &event["biological_contamination"];
    let (contamination_min, contamination_max, recipients) = if contamination.is_object() {
        (
            number(
                &contamination["percent_range"]["minimum"],
                "percent_range.minimum",
            ),
            number(
                &contamination["percent_range"]["maximum"],
                "percent_range.maximum",
            ),
            contamination["recipient_ids"]
                .as_array()
                .expect("receptoras")
                .iter()
                .map(|value| text(value, "recipient_id").to_owned())
                .collect(),
        )
    } else {
        (0.0, 0.0, libraries.keys().cloned().collect())
    };

    let stutter = read_stutter(&event["stutter"]);
    let pcr = read_pcr(
        &artifacts["pcr"],
        text(&profile_view["preparation"], "preparation"),
    );
    EventView {
        typed: typed_event,
        event_id,
        libraries,
        read_length: integer(&event["read_length"], "read_length"),
        target_coverage: number(&event["target"]["coverage_wgs"], "target.coverage_wgs"),
        chimera_rate: number(&artifacts["chimera"]["rate"], "chimera.rate"),
        hopping_rate: number(&artifacts["index_hopping"]["rate"], "index_hopping.rate"),
        contamination_min,
        contamination_max,
        recipients,
        errors_enabled: text(
            &event["sequencing_errors"]["effective_mode"],
            "sequencing_errors.effective_mode",
        ) == "enabled",
        stutter,
        pcr,
        profile,
    }
}

fn read_stutter(stutter: &Value) -> StutterParams {
    if !stutter.is_object() {
        return StutterParams {
            enabled: false,
            motif_rates: [0.0; 6],
            coefficient: 0.0,
            step_weights: [0.0; 4],
            locus_rates: BTreeMap::new(),
            factors: BTreeMap::new(),
        };
    }
    let mut motif_rates = [0.0; 6];
    for entry in stutter["motif_length_rates"].as_array().expect("tasas") {
        let length = integer(&entry["motif_length_bases"], "motif_length_bases");
        motif_rates[usize::try_from(length - 1).expect("longitud de motivo")] =
            number(&entry["rate"], "rate");
    }
    let mut step_weights = [0.0; 4];
    for entry in stutter["step_weights"].as_array().expect("pesos") {
        let units = entry["units"].as_i64().expect("unidades del paso");
        let slot = ADMITTED_STEPS
            .iter()
            .position(|step| *step == units)
            .expect("paso admitido");
        step_weights[slot] = number(&entry["weight"], "weight");
    }
    let locus_rates = stutter["locus_rates"]
        .as_array()
        .map(|entries| {
            entries
                .iter()
                .map(|entry| {
                    (
                        text(&entry["locus_id"], "locus_rates.locus_id").to_owned(),
                        number(&entry["rate"], "locus_rates.rate"),
                    )
                })
                .collect()
        })
        .unwrap_or_default();
    let mut factors = BTreeMap::new();
    for entry in stutter["tract_length_dependence"]["by_motif_length"]
        .as_array()
        .expect("dependencia con la longitud del tracto")
    {
        if entry["measured"] != true {
            continue;
        }
        let table = entry["factors"]
            .as_array()
            .expect("factores")
            .iter()
            .map(|factor| {
                (
                    integer(&factor["units"], "factor.units"),
                    number(&factor["factor"], "factor.factor"),
                )
            })
            .collect::<Vec<_>>();
        factors.insert(
            integer(&entry["motif_length_bases"], "motif_length_bases"),
            table,
        );
    }
    StutterParams {
        enabled: stutter["enabled"] == true,
        motif_rates,
        coefficient: number(
            &stutter["run_length_coefficient"]["value"],
            "run_length_coefficient",
        ),
        step_weights,
        locus_rates,
        factors,
    }
}

fn read_pcr(pcr: &Value, preparation: &str) -> Option<PcrParams> {
    if preparation != "pcr" {
        return None;
    }
    assert_eq!(
        pcr["applicable"], true,
        "preparación amplificada sin PCR aplicable"
    );
    let cycles = u32::try_from(integer(&pcr["cycles"], "pcr.cycles")).expect("ciclos");
    let clonal = pcr["stutter"].is_object().then(|| {
        let stutter = &pcr["stutter"];
        let mut per_motif = [(false, 0.0, 0.0, 0.0); 6];
        for entry in stutter["per_cycle_slip"].as_array().expect("por ciclo") {
            let length = integer(&entry["motif_length_bases"], "motif_length_bases");
            per_motif[usize::try_from(length - 1).expect("longitud")] = (
                text(&entry["source"], "source") == "declared",
                number(&entry["probability"], "probability"),
                number(&entry["with_pcr_fraction"], "with_pcr_fraction"),
                number(&entry["without_pcr_fraction"], "without_pcr_fraction"),
            );
        }
        ClonalParams {
            cycles: u32::try_from(integer(&stutter["cycles"], "pcr.stutter.cycles"))
                .expect("ciclos"),
            per_motif,
        }
    });
    Some(PcrParams {
        cycles,
        classes: pcr["efficiency_classes"]
            .as_array()
            .expect("clases de eficiencia")
            .iter()
            .map(|class| {
                (
                    number(&class["efficiency"], "efficiency"),
                    number(&class["weight"], "weight"),
                )
            })
            .collect(),
        early_cycles: u32::try_from(integer(&pcr["early_cycles"], "early_cycles"))
            .expect("ciclos tempranos"),
        max_family_size: integer(&pcr["max_family_size"], "max_family_size"),
        gc_max_attempts: integer(&pcr["gc_max_attempts"], "gc_max_attempts"),
        gc_curve: pcr["gc_relative_efficiency"]
            .as_array()
            .expect("curva de composición")
            .iter()
            .map(|knot| {
                (
                    number(&knot["gc_fraction"], "gc_fraction"),
                    number(&knot["relative_efficiency"], "relative_efficiency"),
                )
            })
            .collect(),
        mean_pairs: number(&pcr["mean_pairs_per_molecule"], "mean_pairs_per_molecule"),
        clonal,
    })
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().fold(String::new(), |mut text, byte| {
        let _ = write!(text, "{byte:02x}");
        text
    })
}

// ---------------------------------------------------------------------------
// FASTQ publicado.
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct FastqRecord {
    header: String,
    sequence: Vec<u8>,
    quality: Vec<u8>,
}

impl FastqRecord {
    /// Identidad de pareja: todo lo anterior al espacio, que R1 y R2 comparten.
    fn pair_key(&self) -> String {
        self.header.split(' ').next().unwrap_or_default().to_owned()
    }
}

#[derive(Clone, Debug)]
struct ObservedRecordPair {
    r1: FastqRecord,
    r2: FastqRecord,
}

/// Descomprime con el `gzip` del sistema: el lector BGZF del programa no
/// interviene en lo que el oráculo ve.
fn read_gzip_fastq(path: &Path) -> Vec<FastqRecord> {
    let output = Command::new("gzip")
        .arg("-dc")
        .arg(path)
        .output()
        .unwrap_or_else(|error| panic!("gzip sobre {}: {error}", path.display()));
    assert!(
        output.status.success(),
        "gzip falló sobre {}: {}",
        path.display(),
        String::from_utf8_lossy(&output.stderr)
    );
    let content = String::from_utf8(output.stdout).expect("FASTQ ASCII");
    let lines = content.lines().collect::<Vec<_>>();
    assert_eq!(lines.len() % 4, 0, "{}: registro truncado", path.display());
    lines
        .chunks(4)
        .map(|record| {
            assert_eq!(record[2], "+", "{}: separador", path.display());
            FastqRecord {
                header: record[0].to_owned(),
                sequence: record[1].as_bytes().to_vec(),
                quality: record[3].as_bytes().to_vec(),
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Identidades técnicas (DTC-087), reescritas con BLAKE3.
// ---------------------------------------------------------------------------

const ID_CONTEXT: &str = "genome_population_simulator/v1/id";
const TECHNICAL_ID_DERIVATION_VERSION: u32 = 2;

fn technical_id(seed: RootSeed, prefix: &str, components: &[(&str, &str)]) -> String {
    let mut material = Vec::new();
    material.extend_from_slice(seed.as_bytes());
    material.extend_from_slice(&TECHNICAL_ID_DERIVATION_VERSION.to_be_bytes());
    material.extend_from_slice(
        &u32::try_from(components.len())
            .expect("componentes")
            .to_be_bytes(),
    );
    for (name, value) in components {
        for part in [name.as_bytes(), value.as_bytes()] {
            material.extend_from_slice(&u32::try_from(part.len()).expect("longitud").to_be_bytes());
            material.extend_from_slice(part);
        }
    }
    let digest = blake3::derive_key(ID_CONTEXT, &material);
    format!("{prefix}{}", base32(&digest[..16]))
}

fn base32(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 32] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
    let mut bits = String::new();
    for byte in bytes {
        let _ = write!(bits, "{byte:08b}");
    }
    while !bits.len().is_multiple_of(5) {
        bits.push('0');
    }
    bits.as_bytes()
        .chunks(5)
        .map(|chunk| {
            let index = chunk
                .iter()
                .fold(0_usize, |value, bit| value * 2 + usize::from(*bit == b'1'));
            char::from(ALPHABET[index])
        })
        .collect()
}

fn own_pre_hop_instance(seed: RootSeed, event_id: &str, library_id: &str, ordinal: u64) -> String {
    let template_id = format!("TPL_{library_id}_{ordinal}");
    technical_id(
        seed,
        "PH_",
        &[
            ("class", "PreHopInstanceId"),
            ("event_id", event_id),
            ("library_id", library_id),
            ("origin", "own"),
            ("template_id", &template_id),
        ],
    )
}

fn contamination_pre_hop_instance(
    seed: RootSeed,
    event_id: &str,
    recipient_library_id: &str,
    source: &str,
) -> String {
    technical_id(
        seed,
        "PH_",
        &[
            ("class", "PreHopInstanceId"),
            ("event_id", event_id),
            ("library_id", recipient_library_id),
            ("origin", "contamination"),
            ("source_pre_hop_instance_id", source),
        ],
    )
}

fn technical_pair(seed: RootSeed, event_id: &str, final_library_id: &str, pre_hop: &str) -> String {
    technical_id(
        seed,
        "TP_",
        &[
            ("class", "TechnicalPairId"),
            ("event_id", event_id),
            ("final_library_id", final_library_id),
            ("pre_hop_instance_id", pre_hop),
        ],
    )
}

// ---------------------------------------------------------------------------
// Reconstrucción.
// ---------------------------------------------------------------------------

fn decision_rng(seed: RootSeed, decision: F4TechnicalDecision<'_>) -> LocalRng {
    derive_f4_technical_subseed(seed, decision)
        .expect("subsemilla del catálogo técnico")
        .local_rng()
}

fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            other => panic!("base imposible {other}"),
        })
        .collect()
}

/// Un tramo de la molécula: de qué copia sale y desde dónde.
#[derive(Clone, Debug)]
struct Segment {
    copy: usize,
    copy_start: u64,
    molecule_start: u64,
    length: u64,
}

#[derive(Clone, Debug)]
struct Molecule {
    copy: usize,
    start: u64,
    bases: Vec<u8>,
}

/// Lo que el oráculo espera de una pareja antes del demultiplexado.
#[derive(Clone, Debug)]
struct OwnPair {
    pre_hop: String,
    r1: Vec<u8>,
    q1: Vec<u8>,
    r2: Vec<u8>,
    q2: Vec<u8>,
    provenance: String,
}

#[derive(Clone, Debug)]
struct ExpectedPair {
    header_prefix: String,
    udi: String,
    r1: Vec<u8>,
    q1: Vec<u8>,
    r2: Vec<u8>,
    q2: Vec<u8>,
    is_contaminant: bool,
    provenance: String,
}

/// Genealogía clonal de una familia: los loci y el paso neto de cada copia.
#[derive(Clone, Debug, Default)]
struct Lineage {
    loci: Vec<String>,
    net: Vec<Vec<i64>>,
}

struct Family {
    molecule: Molecule,
    remaining: u64,
    next_member: usize,
    lineage: Lineage,
    family_ordinal: u64,
}

/// Tracto declarado situado sobre la molécula: `[start, end)` en coordenadas de
/// la molécula.
#[derive(Clone, Debug)]
struct PlacedTract {
    locus: usize,
    start: u64,
    end: u64,
    /// Sale del primer tramo, el de la molécula original de la familia: solo
    /// ahí cae el paso heredado del linaje.
    first_segment: bool,
    /// Sale de un tramo cuya copia es de otra secuencia que la del primero.
    other_sequence: bool,
}

/// Tirada más larga de un tracto: unidades, inicio relativo y longitud de su
/// motivo. Se recorre de izquierda a derecha probando los motivos del locus en
/// su orden declarado; el empate conserva la tirada de más a la izquierda
/// (REQ-DAT-77).
fn longest_run(bases: &[u8], motifs: &[Vec<u8>]) -> (u64, u64, u64) {
    let (mut best, mut best_start, mut best_motif) = (0_u64, 0_u64, 0_u64);
    let (mut run, mut run_start) = (0_u64, 0_u64);
    let mut run_motif: Option<usize> = None;
    let mut offset = 0_usize;
    while offset < bases.len() {
        let matched = motifs
            .iter()
            .position(|motif| bases[offset..].starts_with(motif));
        match matched {
            Some(index) => {
                if run_motif == Some(index) {
                    run += 1;
                } else {
                    run_motif = Some(index);
                    run = 1;
                    run_start = offset as u64;
                }
                if run > best {
                    best = run;
                    best_start = run_start;
                    best_motif = motifs[index].len() as u64;
                }
                offset += motifs[index].len();
            }
            None => {
                run_motif = None;
                run = 0;
                offset += 1;
            }
        }
    }
    (best, best_start, best_motif)
}

/// ¿Hay dos tiradas que alcanzan la longitud máxima? Recorre las tiradas
/// completas por su cuenta, sin reutilizar `longest_run`.
fn tied_longest_run(bases: &[u8], motifs: &[Vec<u8>]) -> bool {
    let mut runs = Vec::new();
    let (mut run, mut run_motif, mut offset) = (0_u64, None, 0_usize);
    while offset < bases.len() {
        match motifs
            .iter()
            .position(|motif| bases[offset..].starts_with(motif))
        {
            Some(index) => {
                if run_motif != Some(index) {
                    if run > 0 {
                        runs.push(run);
                    }
                    run = 0;
                    run_motif = Some(index);
                }
                run += 1;
                offset += motifs[index].len();
            }
            None => {
                if run > 0 {
                    runs.push(run);
                }
                run = 0;
                run_motif = None;
                offset += 1;
            }
        }
    }
    if run > 0 {
        runs.push(run);
    }
    let best = runs.iter().copied().max().unwrap_or(0);
    runs.iter().filter(|units| **units == best).count() > 1
}

fn viable(
    run_units: u64,
    motif_length: u64,
    molecule_length: u64,
    read_length: u64,
    step: i64,
) -> bool {
    if step > 0 {
        return true;
    }
    let removed = step.unsigned_abs();
    removed < run_units
        && molecule_length
            .checked_sub(removed * motif_length)
            .is_some_and(|length| length >= read_length)
}

impl PublishedRun {
    fn expected_final_pairs(
        &self,
        tally: &mut OracleTally,
    ) -> BTreeMap<String, BTreeMap<String, ExpectedPair>> {
        let event = &self.event;
        // 1. Parejas propias de cada biblioteca, ya con calidad y errores.
        let mut own: BTreeMap<String, Vec<OwnPair>> = BTreeMap::new();
        for library in event.libraries.values() {
            let individual = self
                .individuals
                .get(&library.individual_id)
                .unwrap_or_else(|| panic!("verdad de {}", library.individual_id));
            own.insert(
                library.individual_id.clone(),
                self.simulate_library(individual, library, tally),
            );
        }
        // 2. Contaminación biológica: qué parejas propias se retiran y qué
        //    parejas de donantes entran en cada receptora.
        let mut emitted: Vec<(String, String, OwnPair, bool)> = Vec::new();
        let recipients = event.recipients.iter().cloned().collect::<BTreeSet<_>>();
        for (individual_id, pairs) in &own {
            let library = &event.libraries[individual_id];
            if !recipients.contains(individual_id) {
                for pair in pairs {
                    emitted.push((
                        library.library_id.clone(),
                        pair.pre_hop.clone(),
                        pair.clone(),
                        false,
                    ));
                }
                continue;
            }
            let recipient = &self.individuals[individual_id];
            let typed_individual =
                IndividualId::parse(individual_id.clone(), "oracle").expect("id");
            let percent = f4_uniform_closed(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::contamination_percent(
                        &event.typed,
                        &library.typed,
                        &typed_individual,
                    ),
                ),
                event.contamination_min,
                event.contamination_max,
            )
            .expect("porcentaje");
            let clean = pairs.len() as u64;
            let raw = clean as f64 * percent / 100.0;
            let target = if raw - raw.floor() >= 0.5 {
                raw.floor() as u64 + 1
            } else {
                raw.floor() as u64
            };
            // Candidatas donantes: las parejas limpias de las demás individuas
            // de la misma población, en orden de individua y de instancia.
            let mut donors = own
                .iter()
                .filter(|(other, _)| {
                    *other != individual_id
                        && self.individuals[*other].population == recipient.population
                })
                .flat_map(|(other, pairs)| {
                    pairs.iter().map(move |pair| (other.clone(), pair.clone()))
                })
                .collect::<Vec<_>>();
            donors.sort_by(|left, right| {
                (&left.0, &left.1.pre_hop).cmp(&(&right.0, &right.1.pre_hop))
            });
            let count = usize::try_from(target).expect("recuento");
            f4_partial_fisher_yates(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::contamination_donors(
                        &event.typed,
                        &library.typed,
                        &typed_individual,
                    ),
                ),
                &mut donors,
                count,
            )
            .expect("donantes");
            let mut removable = pairs.clone();
            removable.sort_by(|left, right| left.pre_hop.cmp(&right.pre_hop));
            f4_partial_fisher_yates(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::contamination_recipients(
                        &event.typed,
                        &library.typed,
                        &typed_individual,
                    ),
                ),
                &mut removable,
                count,
            )
            .expect("retiradas");
            let removed = removable[..count]
                .iter()
                .map(|pair| pair.pre_hop.clone())
                .collect::<BTreeSet<_>>();
            for pair in pairs {
                if !removed.contains(&pair.pre_hop) {
                    emitted.push((
                        library.library_id.clone(),
                        pair.pre_hop.clone(),
                        pair.clone(),
                        false,
                    ));
                }
            }
            for (donor, pair) in &donors[..count] {
                let pre_hop = contamination_pre_hop_instance(
                    self.seed,
                    &event.event_id,
                    &library.library_id,
                    &pair.pre_hop,
                );
                let mut contaminant = pair.clone();
                contaminant.provenance = format!(
                    "contaminante de {donor} en {individual_id} ({})",
                    pair.provenance
                );
                emitted.push((library.library_id.clone(), pre_hop, contaminant, true));
                tally.contaminant_pairs += 1;
            }
        }
        // 3. Demultiplexado con salto de índice e identidad final.
        let by_library = event
            .libraries
            .values()
            .map(|library| (library.library_id.clone(), library))
            .collect::<BTreeMap<_, _>>();
        let mut expected: BTreeMap<String, BTreeMap<String, ExpectedPair>> = BTreeMap::new();
        for (pre_hop_library, pre_hop, pair, is_contaminant) in emitted {
            let source = by_library[&pre_hop_library];
            let mut final_library = pre_hop_library.clone();
            if event.hopping_rate > 0.0 {
                let instance =
                    genome_population_simulator::randomness::F4IndexHopInstance::parse(&pre_hop)
                        .expect("instancia PH_");
                let hops = f4_bernoulli(
                    &mut decision_rng(
                        self.seed,
                        F4TechnicalDecision::index_hop_bernoulli(
                            &event.typed,
                            &source.typed,
                            instance,
                        ),
                    ),
                    event.hopping_rate,
                )
                .expect("salto");
                if hops {
                    let others = by_library
                        .keys()
                        .filter(|library| **library != pre_hop_library)
                        .cloned()
                        .collect::<Vec<_>>();
                    let index = f4_uniform_index(
                        &mut decision_rng(
                            self.seed,
                            F4TechnicalDecision::index_hop_destination(
                                &event.typed,
                                &source.typed,
                                instance,
                            ),
                        ),
                        others.len() as u64,
                    )
                    .expect("destino");
                    final_library = others[usize::try_from(index).expect("índice")].clone();
                    tally.hopped_pairs += 1;
                }
            }
            let destination = by_library[&final_library];
            let technical = technical_pair(self.seed, &event.event_id, &final_library, &pre_hop);
            let header_prefix = format!("@GPS:{}:{final_library}:{technical}", event.event_id);
            let previous = expected.entry(final_library.clone()).or_default().insert(
                header_prefix.clone(),
                ExpectedPair {
                    header_prefix,
                    udi: format!("{}+{}", destination.i7, destination.i5),
                    r1: pair.r1,
                    q1: pair.q1,
                    r2: pair.r2,
                    q2: pair.q2,
                    is_contaminant,
                    provenance: format!(
                        "{} [biblioteca previa al salto {}]",
                        pair.provenance, source.individual_id
                    ),
                },
            );
            assert!(
                previous.is_none(),
                "el oráculo derivó dos veces la misma pareja final"
            );
        }
        expected
    }

    /// Parejas propias de una biblioteca, por ordinal de plantilla.
    fn simulate_library(
        &self,
        individual: &IndividualView,
        library: &LibraryView,
        tally: &mut OracleTally,
    ) -> Vec<OwnPair> {
        let event = &self.event;
        let profile = &event.profile;
        // Número de parejas: objetivo por la referencia haploide entera entre
        // las bases de una pareja, al entero más próximo con empate hacia
        // arriba.
        let raw =
            event.target_coverage * self.reference_length as f64 / (2 * event.read_length) as f64;
        let pair_count = raw.round() as u64;
        // El techo efectivo del inserto es la copia más larga capaz de alojar
        // el mínimo, sin pasar del máximo del perfil.
        let Some(longest) = individual
            .copies
            .iter()
            .map(|copy| copy.bases.len() as u64)
            .filter(|length| *length >= profile.insert_min)
            .max()
        else {
            return Vec::new();
        };
        let insert_max = longest.min(profile.insert_max);
        let mut pairs = Vec::new();
        let mut family: Option<Family> = None;
        for ordinal in 0..pair_count {
            // Molécula de esta pareja, y paso neto heredado si hay familia.
            let (molecule, lineage_steps, provenance_family) = match &event.pcr {
                None => (
                    self.sample_molecule(individual, library, insert_max, ordinal, 0),
                    None,
                    String::new(),
                ),
                Some(pcr) => {
                    if family.as_ref().is_none_or(|open| open.remaining == 0) {
                        family = Some(self.open_family(
                            individual,
                            library,
                            pcr,
                            insert_max,
                            ordinal,
                            pair_count - ordinal,
                            tally,
                        ));
                    }
                    let open = family.as_mut().expect("familia abierta");
                    let member = open.next_member;
                    open.next_member += 1;
                    open.remaining -= 1;
                    if member > 0 {
                        tally.duplicate_pairs += 1;
                    }
                    let steps = (!open.lineage.loci.is_empty()).then(|| {
                        open.lineage
                            .loci
                            .iter()
                            .cloned()
                            .zip(open.lineage.net[member].iter().copied())
                            .collect::<Vec<_>>()
                    });
                    (
                        open.molecule.clone(),
                        steps,
                        format!(" familia {} copia {member}", open.family_ordinal),
                    )
                }
            };
            let copy = &individual.copies[molecule.copy];
            let insert = molecule.bases.len() as u64;
            let mut provenance = format!(
                "{} {}#{} [{}, {}){provenance_family} ordinal {ordinal}",
                individual.id,
                copy.sequence_id,
                copy.copy_index,
                molecule.start,
                molecule.start + insert
            );
            tally.copy_indices_sampled.insert(copy.copy_index);
            tally.sequences_sampled.insert(copy.sequence_id.clone());
            for mark in &copy.origin[usize::try_from(molecule.start).expect("inicio")
                ..usize::try_from(molecule.start + insert).expect("fin")]
            {
                match mark {
                    1 => tally.snv_bases_sampled += 1,
                    2 => tally.inserted_bases_sampled += 1,
                    3 => tally.repeat_bases_sampled += 1,
                    _ => {}
                }
            }
            // Quimera (DTC-084), decidida por pareja sobre la molécula intacta.
            let (bases, segments) = self.chimera(
                individual,
                library,
                ordinal,
                &molecule,
                &mut provenance,
                tally,
            );
            // Tartamudeo: primero el paso heredado del linaje y después la
            // decisión de base, tracto a tracto (REQ-DAT-75 a 79).
            let bases = self.stutter(
                individual,
                library,
                ordinal,
                bases,
                &segments,
                lineage_steps.as_deref(),
                &mut provenance,
                tally,
            );
            // Orientación FR: R1 es la hebra directa del inserto y R2 el
            // complemento inverso de su extremo (REQ-DAT-12).
            let read_length = usize::try_from(event.read_length).expect("longitud de lectura");
            assert!(
                bases.len() >= read_length,
                "inserto más corto que la lectura"
            );
            let r1_true = bases[..read_length].to_vec();
            let r2_true = reverse_complement(&bases[bases.len() - read_length..]);
            let (r1, q1) = self.observe_mate(library, ordinal, F4Mate::R1, &r1_true, tally);
            let (r2, q2) = self.observe_mate(library, ordinal, F4Mate::R2, &r2_true, tally);
            pairs.push(OwnPair {
                pre_hop: own_pre_hop_instance(
                    self.seed,
                    &event.event_id,
                    &library.library_id,
                    ordinal,
                ),
                r1,
                q1,
                r2,
                q2,
                provenance,
            });
        }
        pairs
    }

    /// Sorteo del fragmento biológico: longitud del inserto, copia y comienzo,
    /// con las claves de `fragment` (el intento `0`) o sus reintentos.
    fn sample_molecule(
        &self,
        individual: &IndividualView,
        library: &LibraryView,
        insert_max: u64,
        ordinal: u64,
        attempt: u64,
    ) -> Molecule {
        let event = &self.event;
        let profile = &event.profile;
        let decision = |kind: FragmentRetryDecision| match NonZeroU64::new(attempt) {
            None => match kind {
                FragmentRetryDecision::InsertLength => F4TechnicalDecision::fragment_insert_length(
                    &event.typed,
                    &library.typed,
                    ordinal,
                ),
                FragmentRetryDecision::Copy => {
                    F4TechnicalDecision::fragment_copy(&event.typed, &library.typed, ordinal)
                }
                FragmentRetryDecision::Start => {
                    F4TechnicalDecision::fragment_start(&event.typed, &library.typed, ordinal)
                }
            },
            Some(attempt) => F4TechnicalDecision::fragment_retry(
                &event.typed,
                &library.typed,
                ordinal,
                kind,
                attempt,
            ),
        };
        let insert = u64::try_from(
            f4_truncated_normal(
                &mut decision_rng(self.seed, decision(FragmentRetryDecision::InsertLength)),
                profile.insert_mean,
                profile.insert_sd,
                i64::try_from(profile.insert_min).expect("mínimo"),
                i64::try_from(insert_max).expect("máximo"),
            )
            .expect("inserto"),
        )
        .expect("inserto no negativo");
        // Cada comienzo admisible de cada copia pesa lo mismo.
        let starts = individual
            .copies
            .iter()
            .map(|copy| (copy.bases.len() as u64 + 1).saturating_sub(insert))
            .collect::<Vec<_>>();
        let total = starts.iter().sum::<u64>();
        let mut selected = f4_uniform_index(
            &mut decision_rng(self.seed, decision(FragmentRetryDecision::Copy)),
            total,
        )
        .expect("copia");
        let copy = starts
            .iter()
            .position(|count| {
                if selected < *count {
                    true
                } else {
                    selected -= count;
                    false
                }
            })
            .expect("copia candidata");
        let start = f4_uniform_index(
            &mut decision_rng(self.seed, decision(FragmentRetryDecision::Start)),
            starts[copy],
        )
        .expect("comienzo");
        let bases = individual.copies[copy].bases[usize::try_from(start).expect("inicio")
            ..usize::try_from(start + insert).expect("fin")]
            .to_vec();
        Molecule { copy, start, bases }
    }

    #[allow(clippy::too_many_arguments)]
    fn open_family(
        &self,
        individual: &IndividualView,
        library: &LibraryView,
        pcr: &PcrParams,
        insert_max: u64,
        ordinal: u64,
        remaining_pairs: u64,
        tally: &mut OracleTally,
    ) -> Family {
        let event = &self.event;
        let class = if pcr.classes.len() == 1 {
            0
        } else {
            let weights = pcr.classes.iter().map(|class| class.1).collect::<Vec<_>>();
            categorical_f64(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::amplification_efficiency_class(
                        &event.typed,
                        &library.typed,
                        ordinal,
                    ),
                ),
                &weights,
            )
            .expect("clase de eficiencia")
        };
        let efficiency = pcr.classes[class].0;
        // Copias tras los primeros ciclos: cada copia se duplica con
        // probabilidad igual a la eficiencia, y la familia observada se elige
        // con probabilidad proporcional a su número de copias.
        let early = if pcr.early_cycles == 0 || pcr.mean_pairs == 0.0 {
            1
        } else {
            let distribution = early_copy_distribution(pcr.early_cycles, efficiency);
            let weighted = distribution
                .iter()
                .enumerate()
                .map(|(index, probability)| (index + 1) as f64 * probability)
                .collect::<Vec<_>>();
            let total = weighted.iter().sum::<f64>();
            let biased = weighted
                .iter()
                .map(|weight| weight / total)
                .collect::<Vec<_>>();
            categorical_f64(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::amplification_early_copies(
                        &event.typed,
                        &library.typed,
                        ordinal,
                    ),
                ),
                &biased,
            )
            .expect("copias tempranas") as u64
                + 1
        };
        // Molécula original con aceptación por composición.
        let mut attempt = 0_u64;
        let (molecule, composition) = loop {
            let molecule = self.sample_molecule(individual, library, insert_max, ordinal, attempt);
            if attempt > 0 {
                tally.gc_retries += 1;
            }
            let relative = evaluate_gc_curve(&pcr.gc_curve, gc_fraction(&molecule.bases));
            let factor = ((1.0 + efficiency * relative) / (1.0 + efficiency))
                .powi(i32::try_from(pcr.cycles).expect("ciclos"))
                .clamp(0.0, 1.0);
            attempt += 1;
            if factor == 1.0 {
                break (molecule, factor);
            }
            let accepted = f4_bernoulli(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::amplification_gc_accept(
                        &event.typed,
                        &library.typed,
                        ordinal,
                        attempt - 1,
                    ),
                ),
                factor,
            )
            .expect("aceptación");
            if accepted || attempt >= pcr.gc_max_attempts {
                break (molecule, factor);
            }
        };
        let mean_early = (1.0 + efficiency).powi(i32::try_from(pcr.early_cycles).expect("ciclos"));
        let lambda = pcr.mean_pairs * (early as f64 / mean_early) * composition;
        let drawn = if lambda == 0.0 {
            1
        } else {
            categorical_f64(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::amplification_family_size(
                        &event.typed,
                        &library.typed,
                        ordinal,
                    ),
                ),
                &zero_truncated_poisson(lambda, pcr.max_family_size),
            )
            .expect("tamaño de familia") as u64
                + 1
        };
        let size = drawn.min(remaining_pairs.max(1));
        let lineage = self.resolve_lineage(library, pcr, ordinal, efficiency, size, &molecule);
        Family {
            molecule,
            remaining: size,
            next_member: 0,
            lineage,
            family_ordinal: ordinal,
        }
    }

    /// Tractos declarados contenidos en la molécula, situados por los tramos
    /// de copia de los que la molécula se compone.
    fn place_tracts(&self, individual: &IndividualView, segments: &[Segment]) -> Vec<PlacedTract> {
        let mut placed = Vec::new();
        for (segment_index, segment) in segments.iter().enumerate() {
            let copy = &individual.copies[segment.copy];
            for (index, locus) in self.loci.iter().enumerate() {
                if locus.sequence_id != copy.sequence_id {
                    continue;
                }
                let Some((start, end)) = locate_tract(copy, locus) else {
                    continue;
                };
                let from = start.max(segment.copy_start);
                let to = end.min(segment.copy_start + segment.length);
                if from >= to {
                    continue;
                }
                placed.push(PlacedTract {
                    locus: index,
                    start: segment.molecule_start + from - segment.copy_start,
                    end: segment.molecule_start + to - segment.copy_start,
                    first_segment: segment_index == 0,
                    other_sequence: copy.sequence_id
                        != individual.copies[segments[0].copy].sequence_id,
                });
            }
        }
        placed.sort_by(|left, right| {
            (left.start, &self.loci[left.locus].id).cmp(&(right.start, &self.loci[right.locus].id))
        });
        placed
    }

    fn tract_factor(&self, motif_length: u64, run_units: u64) -> f64 {
        let Some(table) = self.event.stutter.factors.get(&motif_length) else {
            return 1.0;
        };
        let first = table.first().expect("tabla").0;
        let last = table.last().expect("tabla").0;
        let units = run_units.clamp(first, last);
        table
            .iter()
            .find(|(candidate, _)| *candidate == units)
            .expect("factor tabulado")
            .1
    }

    fn resolve_lineage(
        &self,
        library: &LibraryView,
        pcr: &PcrParams,
        family_ordinal: u64,
        efficiency: f64,
        size: u64,
        molecule: &Molecule,
    ) -> Lineage {
        let stutter = &self.event.stutter;
        let Some(clonal) = &pcr.clonal else {
            return Lineage::default();
        };
        if !stutter.enabled || self.loci.is_empty() || clonal.cycles == 0 {
            return Lineage::default();
        }
        let individual_id = self
            .event
            .libraries
            .values()
            .find(|candidate| candidate.library_id == library.library_id)
            .expect("biblioteca")
            .individual_id
            .clone();
        let individual = &self.individuals[&individual_id];
        let segments = [Segment {
            copy: molecule.copy,
            copy_start: molecule.start,
            molecule_start: 0,
            length: molecule.bases.len() as u64,
        }];
        // Tractos con su tirada más larga en la molécula original y su
        // probabilidad por ciclo.
        let mut tracts = Vec::new();
        for placed in self.place_tracts(individual, &segments) {
            let locus = &self.loci[placed.locus];
            let window = &molecule.bases[usize::try_from(placed.start).expect("inicio")
                ..usize::try_from(placed.end).expect("fin")];
            let (units, _, motif_length) = longest_run(window, &locus.motifs);
            if units == 0 || !(1..=6).contains(&motif_length) {
                continue;
            }
            let (declared, probability, with_pcr, without_pcr) =
                clonal.per_motif[usize::try_from(motif_length - 1).expect("longitud")];
            let probability = if declared {
                probability
            } else {
                let factor = self.tract_factor(motif_length, units);
                per_cycle_probability(with_pcr * factor, without_pcr * factor, clonal.cycles)
            };
            tracts.push((placed.locus, units, motif_length, probability));
        }
        if tracts.is_empty() {
            return Lineage::default();
        }
        let cycles = usize::try_from(clonal.cycles).expect("ciclos");
        // Ciclo de divergencia: se separa del progenitor en el ciclo k con
        // probabilidad (1 − r)·r^k, r = 1/(1 + eficiencia); el último recoge
        // la cola.
        let ratio = 1.0 / (1.0 + efficiency);
        let mut divergence = (0..cycles)
            .map(|cycle| (1.0 - ratio) * ratio.powi(i32::try_from(cycle).expect("ciclo")))
            .collect::<Vec<_>>();
        divergence.push(ratio.powi(i32::try_from(cycles).expect("ciclos")));
        let members = usize::try_from(size).expect("tamaño");
        let event = &self.event;
        let mut history = vec![vec![vec![0_i64; cycles]; tracts.len()]; members];
        let mut net_by_member = Vec::with_capacity(members);
        for member in 0..members {
            let copy_index = member as u64;
            let (parent, diverges_at) = if member == 0 {
                (0, 0)
            } else {
                let parent = f4_uniform_index(
                    &mut decision_rng(
                        self.seed,
                        F4TechnicalDecision::lineage_parent(
                            &event.typed,
                            &library.typed,
                            family_ordinal,
                            copy_index,
                        ),
                    ),
                    copy_index,
                )
                .expect("progenitor");
                let cycle = categorical_f64(
                    &mut decision_rng(
                        self.seed,
                        F4TechnicalDecision::lineage_divergence(
                            &event.typed,
                            &library.typed,
                            family_ordinal,
                            copy_index,
                        ),
                    ),
                    &divergence,
                )
                .expect("divergencia");
                (usize::try_from(parent).expect("progenitor"), cycle)
            };
            let mut net = vec![0_i64; tracts.len()];
            let mut length = molecule.bases.len() as i64;
            for cycle in 1..=cycles {
                for (index, (locus, units, motif_length, probability)) in tracts.iter().enumerate()
                {
                    let inherited = member > 0 && cycle <= diverges_at;
                    let step = if inherited {
                        history[parent][index][cycle - 1]
                    } else {
                        let typed = &self.loci[*locus].typed;
                        let slips = f4_bernoulli(
                            &mut decision_rng(
                                self.seed,
                                F4TechnicalDecision::lineage_slip(
                                    &event.typed,
                                    &library.typed,
                                    family_ordinal,
                                    copy_index,
                                    cycle as u64,
                                    typed,
                                ),
                            ),
                            *probability,
                        )
                        .expect("resbalón");
                        let run = *units as i64 + net[index];
                        if !slips || run <= 0 || length <= 0 {
                            0
                        } else {
                            let candidates = ADMITTED_STEPS
                                .iter()
                                .copied()
                                .filter(|step| {
                                    viable(
                                        run as u64,
                                        *motif_length,
                                        length as u64,
                                        event.read_length,
                                        *step,
                                    )
                                })
                                .collect::<Vec<_>>();
                            if candidates.is_empty() {
                                0
                            } else {
                                self.choose_step(
                                    &candidates,
                                    F4TechnicalDecision::lineage_step(
                                        &event.typed,
                                        &library.typed,
                                        family_ordinal,
                                        copy_index,
                                        cycle as u64,
                                        typed,
                                    ),
                                )
                                .unwrap_or(0)
                            }
                        }
                    };
                    if step == 0 {
                        continue;
                    }
                    history[member][index][cycle - 1] = step;
                    net[index] += step;
                    length += step * *motif_length as i64;
                }
            }
            net_by_member.push(net);
        }
        Lineage {
            loci: tracts
                .iter()
                .map(|tract| self.loci[tract.0].id.clone())
                .collect(),
            net: net_by_member,
        }
    }

    /// Elige un paso entre los viables con los pesos publicados; sin peso
    /// positivo entre ellos no hay paso y no se sortea.
    fn choose_step(&self, candidates: &[i64], decision: F4TechnicalDecision<'_>) -> Option<i64> {
        let weights = candidates
            .iter()
            .map(|step| {
                self.event.stutter.step_weights[ADMITTED_STEPS
                    .iter()
                    .position(|admitted| admitted == step)
                    .expect("paso admitido")]
            })
            .collect::<Vec<_>>();
        if !weights.iter().any(|weight| *weight > 0.0) {
            return None;
        }
        let chosen =
            categorical_f64(&mut decision_rng(self.seed, decision), &weights).expect("paso");
        Some(candidates[chosen])
    }

    fn chimera(
        &self,
        individual: &IndividualView,
        library: &LibraryView,
        ordinal: u64,
        molecule: &Molecule,
        provenance: &mut String,
        tally: &mut OracleTally,
    ) -> (Vec<u8>, Vec<Segment>) {
        let event = &self.event;
        let whole = vec![Segment {
            copy: molecule.copy,
            copy_start: molecule.start,
            molecule_start: 0,
            length: molecule.bases.len() as u64,
        }];
        let selected = f4_bernoulli(
            &mut decision_rng(
                self.seed,
                F4TechnicalDecision::chimera_bernoulli(&event.typed, &library.typed, ordinal),
            ),
            event.chimera_rate,
        )
        .expect("quimera");
        if !selected {
            return (molecule.bases.clone(), whole);
        }
        let insert = molecule.bases.len() as u64;
        let breakpoint = f4_uniform_index(
            &mut decision_rng(
                self.seed,
                F4TechnicalDecision::chimera_breakpoint(&event.typed, &library.typed, ordinal),
            ),
            insert - 1,
        )
        .expect("punto de corte")
            + 1;
        let suffix = insert - breakpoint;
        let prefix_end = molecule.start + breakpoint;
        // Candidatas: intervalos de longitud `suffix` de todas las copias del
        // mismo individuo, ordenadas por identidad de copia (secuencia y
        // después índice) y por inicio; en la copia fuente se excluyen los que
        // solapan o tocan el prefijo.
        let mut order = (0..individual.copies.len()).collect::<Vec<_>>();
        order.sort_by(|left, right| {
            let left = &individual.copies[*left];
            let right = &individual.copies[*right];
            (&left.sequence_id, left.copy_index).cmp(&(&right.sequence_id, right.copy_index))
        });
        let mut ranges = Vec::new();
        for index in order {
            let length = individual.copies[index].bases.len() as u64;
            let Some(last) = length.checked_sub(suffix) else {
                continue;
            };
            for first in 0..=last {
                let touches = index == molecule.copy
                    && first <= prefix_end
                    && first + suffix >= molecule.start;
                if !touches {
                    ranges.push((index, first));
                }
            }
        }
        let chosen = f4_uniform_index(
            &mut decision_rng(
                self.seed,
                F4TechnicalDecision::chimera_second_source(
                    &event.typed,
                    &library.typed,
                    ordinal,
                    NonZeroU64::new(breakpoint).expect("punto de corte positivo"),
                ),
            ),
            ranges.len() as u64,
        )
        .expect("segundo origen");
        let (second, second_start) = ranges[usize::try_from(chosen).expect("índice")];
        let mut bases = molecule.bases[..usize::try_from(breakpoint).expect("corte")].to_vec();
        bases.extend_from_slice(
            &individual.copies[second].bases[usize::try_from(second_start).expect("inicio")
                ..usize::try_from(second_start + suffix).expect("fin")],
        );
        tally.chimeric_pairs += 1;
        if individual.copies[second].sequence_id != individual.copies[molecule.copy].sequence_id {
            tally.chimeras_across_sequences += 1;
        }
        let _ = write!(
            provenance,
            " quimera corte {breakpoint} segundo {}#{} [{second_start}, {})",
            individual.copies[second].sequence_id,
            individual.copies[second].copy_index,
            second_start + suffix
        );
        (
            bases,
            vec![
                Segment {
                    copy: molecule.copy,
                    copy_start: molecule.start,
                    molecule_start: 0,
                    length: breakpoint,
                },
                Segment {
                    copy: second,
                    copy_start: second_start,
                    molecule_start: breakpoint,
                    length: suffix,
                },
            ],
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn stutter(
        &self,
        individual: &IndividualView,
        library: &LibraryView,
        ordinal: u64,
        mut bases: Vec<u8>,
        segments: &[Segment],
        lineage: Option<&[(String, i64)]>,
        provenance: &mut String,
        tally: &mut OracleTally,
    ) -> Vec<u8> {
        let event = &self.event;
        let stutter = &event.stutter;
        if !stutter.enabled || self.loci.is_empty() {
            return bases;
        }
        let mut shift = 0_i64;
        let placements = self.place_tracts(individual, segments);
        let loci = placements
            .iter()
            .map(|placed| placed.locus)
            .collect::<BTreeSet<_>>();
        if loci.len() < placements.len() {
            tally.same_locus_in_both_segments += 1;
        }
        for placed in placements {
            let locus = &self.loci[placed.locus];
            let window_start = usize::try_from(placed.start as i64 + shift).expect("ventana");
            let mut window_end = usize::try_from(placed.end as i64 + shift).expect("ventana");
            // 1. Paso neto heredado del linaje, si cabe en la molécula que se
            //    emite.
            let inherited = lineage
                .and_then(|steps| steps.iter().find(|(id, _)| *id == locus.id))
                .map_or(0, |(_, step)| *step);
            let net = if placed.first_segment {
                inherited
            } else {
                if inherited != 0 {
                    tally.lineage_steps_withheld_in_suffix += 1;
                }
                0
            };
            if net != 0 {
                let (units, run_start, motif_length) =
                    longest_run(&bases[window_start..window_end], &locus.motifs);
                if units > 0
                    && viable(
                        units,
                        motif_length,
                        bases.len() as u64,
                        event.read_length,
                        net,
                    )
                {
                    apply_step(
                        &mut bases,
                        window_start + run_start as usize,
                        units,
                        motif_length,
                        net,
                    );
                    let delta = net * motif_length as i64;
                    window_end = usize::try_from(window_end as i64 + delta).expect("ventana");
                    shift += delta;
                    tally.lineage_steps += 1;
                    let _ = write!(provenance, " linaje {}{net:+}", locus.id);
                }
            }
            // 2. Decisión de base por pareja sobre la molécula ya actualizada.
            let (units, run_start, motif_length) =
                longest_run(&bases[window_start..window_end], &locus.motifs);
            if units == 0 {
                continue;
            }
            let base_rate = match stutter.locus_rates.get(&locus.id) {
                Some(rate) => *rate,
                None => {
                    if !(1..=6).contains(&motif_length) {
                        continue;
                    }
                    stutter.motif_rates[usize::try_from(motif_length - 1).expect("longitud")]
                        * self.tract_factor(motif_length, units)
                }
            };
            let candidates = ADMITTED_STEPS
                .iter()
                .copied()
                .filter(|step| {
                    viable(
                        units,
                        motif_length,
                        bases.len() as u64,
                        event.read_length,
                        *step,
                    )
                })
                .collect::<Vec<_>>();
            if candidates.is_empty() {
                continue;
            }
            let probability = stutter
                .coefficient
                .mul_add(units as f64, base_rate)
                .clamp(0.0, 1.0);
            let slips = f4_bernoulli(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::stutter_bernoulli(
                        &event.typed,
                        &library.typed,
                        ordinal,
                        &locus.typed,
                    ),
                ),
                probability,
            )
            .expect("tartamudeo");
            if !slips {
                continue;
            }
            let Some(step) = self.choose_step(
                &candidates,
                F4TechnicalDecision::stutter_step(
                    &event.typed,
                    &library.typed,
                    ordinal,
                    &locus.typed,
                ),
            ) else {
                continue;
            };
            if tied_longest_run(&bases[window_start..window_end], &locus.motifs) {
                tally.stutter_steps_on_tied_runs += 1;
            }
            if placed.other_sequence {
                tally.suffix_stutter_steps_from_other_sequence += 1;
            }
            apply_step(
                &mut bases,
                window_start + run_start as usize,
                units,
                motif_length,
                step,
            );
            shift += step * motif_length as i64;
            tally.base_stutter_steps += 1;
            let _ = write!(provenance, " tartamudeo {}{step:+}", locus.id);
        }
        bases
    }

    fn observe_mate(
        &self,
        library: &LibraryView,
        ordinal: u64,
        mate: F4Mate,
        bases: &[u8],
        tally: &mut OracleTally,
    ) -> (Vec<u8>, Vec<u8>) {
        let event = &self.event;
        let curve = match mate {
            F4Mate::R1 => &event.profile.quality_r1,
            F4Mate::R2 => &event.profile.quality_r2,
        };
        let mut observed = Vec::with_capacity(bases.len());
        let mut quality = Vec::with_capacity(bases.len());
        for (index, base) in bases.iter().copied().enumerate() {
            let cycle = index as u64 + 1;
            if base == b'N' {
                observed.push(b'N');
                quality.push(b'!');
                tally.n_bases += 1;
                continue;
            }
            let phred = curve_phred(curve, cycle);
            // La calidad emitida es la probabilidad realizada: -10·log10 de
            // 10^(-Q/10), que para un Q entero devuelve el propio Q.
            quality.push(u8::try_from(phred + 33).expect("carácter de calidad"));
            if !event.errors_enabled {
                observed.push(base);
                continue;
            }
            let probability = 10_f64.powf(-(phred as f64) / 10.0);
            let errs = f4_bernoulli(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::error_bernoulli(
                        &event.typed,
                        &library.typed,
                        ordinal,
                        mate,
                        cycle,
                    ),
                ),
                probability,
            )
            .expect("error");
            if !errs {
                observed.push(base);
                continue;
            }
            let alternatives = b"ACGT"
                .iter()
                .copied()
                .filter(|candidate| *candidate != base)
                .collect::<Vec<_>>();
            let choice = f4_uniform_index(
                &mut decision_rng(
                    self.seed,
                    F4TechnicalDecision::error_substitution(
                        &event.typed,
                        &library.typed,
                        ordinal,
                        mate,
                        cycle,
                    ),
                ),
                3,
            )
            .expect("sustitución");
            observed.push(alternatives[usize::try_from(choice).expect("índice")]);
            tally.sequencing_errors += 1;
        }
        (observed, quality)
    }
}

/// Phred entero de un ciclo: interpolación lineal entre los extremos del
/// tramo, redondeada a la mitad hacia arriba.
fn curve_phred(segments: &[QualitySegment], cycle: u64) -> i64 {
    let segment = segments
        .iter()
        .find(|segment| {
            (u64::from(segment.first_cycle)..=u64::from(segment.last_cycle)).contains(&cycle)
        })
        .expect("ciclo cubierto por la curva");
    let span = i64::from(segment.last_cycle - segment.first_cycle);
    if span == 0 {
        return segment.first_phred;
    }
    let offset = cycle as i64 - i64::from(segment.first_cycle);
    // Q = primero + (último − primero)·offset/span, a la mitad hacia arriba:
    // floor(x + 1/2) con x racional.
    let numerator = 2
        * (segment.first_phred * span + (segment.last_phred - segment.first_phred) * offset)
        + span;
    numerator.div_euclid(2 * span)
}

fn per_cycle_probability(with_pcr: f64, without_pcr: f64, cycles: u32) -> f64 {
    if cycles == 0 || with_pcr <= without_pcr || !(0.0..1.0).contains(&without_pcr) {
        return 0.0;
    }
    let increment = 1.0 - (1.0 - with_pcr) / (1.0 - without_pcr);
    if increment <= 0.0 {
        return 0.0;
    }
    if increment >= 1.0 {
        return 1.0;
    }
    (1.0 - (1.0 - increment).powf(1.0 / f64::from(cycles))).clamp(0.0, 1.0)
}

fn locate_tract(copy: &CopyView, locus: &Locus) -> Option<(u64, u64)> {
    for piece in &copy.pieces {
        match &piece.kind {
            PieceKind::Repeat { locus_id } if *locus_id == locus.id => {
                return Some((piece.offset, piece.offset + piece.length));
            }
            PieceKind::Source { start, end } if *start <= locus.start && locus.end <= *end => {
                let from = piece.offset + locus.start - start;
                return Some((from, from + locus.end - locus.start));
            }
            _ => {}
        }
    }
    None
}

fn apply_step(bases: &mut Vec<u8>, run_start: usize, units: u64, motif_length: u64, step: i64) {
    let motif_length = usize::try_from(motif_length).expect("motivo");
    let run_end = run_start + usize::try_from(units).expect("unidades") * motif_length;
    if step > 0 {
        let motif = bases[run_start..run_start + motif_length].to_vec();
        let inserted = motif.repeat(usize::try_from(step).expect("paso"));
        bases.splice(run_end..run_end, inserted);
    } else {
        let removed = usize::try_from(step.unsigned_abs()).expect("paso") * motif_length;
        bases.drain(run_end - removed..run_end);
    }
}

/// Distribución del número de copias tras `cycles` ciclos tempranos: cada
/// copia se duplica con probabilidad `efficiency`. Índice `i` = `i + 1` copias.
fn early_copy_distribution(cycles: u32, efficiency: f64) -> Vec<f64> {
    let states = 1_usize << cycles;
    let mut distribution = vec![0.0; states];
    distribution[0] = 1.0;
    for _ in 0..cycles {
        let mut next = vec![0.0; states];
        for (index, probability) in distribution.iter().enumerate() {
            if *probability == 0.0 {
                continue;
            }
            let copies = index + 1;
            // Binomial(copies, efficiency) por la fórmula directa.
            let mut coefficient = 1.0_f64;
            for successes in 0..=copies {
                if successes > 0 {
                    coefficient = coefficient * (copies - successes + 1) as f64 / successes as f64;
                }
                let weight = coefficient
                    * efficiency.powi(i32::try_from(successes).expect("éxitos"))
                    * (1.0 - efficiency).powi(i32::try_from(copies - successes).expect("fracasos"));
                next[copies + successes - 1] += probability * weight;
            }
        }
        distribution = next;
    }
    distribution
}

/// Poisson truncado en cero y acotado: el último tamaño recoge toda la cola.
fn zero_truncated_poisson(lambda: f64, max_family_size: u64) -> Vec<f64> {
    let slots = usize::try_from(max_family_size.max(1)).expect("tamaño máximo");
    let mut weights = vec![0.0; slots];
    let tail = (-lambda).exp();
    if slots == 1 || lambda.is_nan() || lambda <= 0.0 || 1.0 - tail <= 0.0 {
        weights[0] = 1.0;
        return weights;
    }
    let mut term = lambda * tail / (1.0 - tail);
    for (index, slot) in weights.iter_mut().enumerate().take(slots - 1) {
        if index > 0 {
            term *= lambda / (index + 1) as f64;
        }
        *slot = term;
    }
    let head = weights[..slots - 1].iter().sum::<f64>();
    weights[slots - 1] = (1.0 - head).max(0.0);
    weights
}

/// Fracción GC sobre las bases A/C/G/T; sin ninguna, 1.
fn gc_fraction(bases: &[u8]) -> f64 {
    let acgt = bases.iter().filter(|base| b"ACGT".contains(base)).count();
    if acgt == 0 {
        return 1.0;
    }
    bases.iter().filter(|base| b"GC".contains(base)).count() as f64 / acgt as f64
}

fn evaluate_gc_curve(knots: &[(f64, f64)], gc: f64) -> f64 {
    let first = knots[0];
    let last = knots[knots.len() - 1];
    if gc <= first.0 {
        return first.1;
    }
    if gc >= last.0 {
        return last.1;
    }
    for window in knots.windows(2) {
        let (left, right) = (window[0], window[1]);
        if gc <= right.0 {
            return left.1 + (right.1 - left.1) * ((gc - left.0) / (right.0 - left.0));
        }
    }
    last.1
}

// ---------------------------------------------------------------------------
// Comparación.
// ---------------------------------------------------------------------------

fn compare_library(
    library: &LibraryView,
    expected: &BTreeMap<String, ExpectedPair>,
    observed: &BTreeMap<String, ObservedRecordPair>,
    mismatches: &mut Vec<String>,
) {
    for (key, pair) in expected {
        let Some(found) = observed.get(key) else {
            mismatches.push(format!(
                "biblioteca {}: falta la pareja {key} ({})",
                library.individual_id, pair.provenance
            ));
            continue;
        };
        for (mate, record, bases, quality) in [
            (1, &found.r1, &pair.r1, &pair.q1),
            (2, &found.r2, &pair.r2, &pair.q2),
        ] {
            let header = format!("{} {mate}:N:0:{}", pair.header_prefix, pair.udi);
            if record.header != header {
                mismatches.push(format!(
                    "biblioteca {}: cabecera R{mate} {} frente a {header} ({})",
                    library.individual_id, record.header, pair.provenance
                ));
            }
            if record.sequence != *bases {
                mismatches.push(format!(
                    "biblioteca {}: bases R{mate} de {key} difieren en {} posiciones ({})\n  publicada {}\n  esperada  {}",
                    library.individual_id,
                    record
                        .sequence
                        .iter()
                        .zip(bases.iter())
                        .filter(|(left, right)| left != right)
                        .count()
                        + record.sequence.len().abs_diff(bases.len()),
                    pair.provenance,
                    String::from_utf8_lossy(&record.sequence),
                    String::from_utf8_lossy(bases)
                ));
            }
            if record.quality != *quality {
                mismatches.push(format!(
                    "biblioteca {}: calidades R{mate} de {key} difieren ({})\n  publicada {}\n  esperada  {}",
                    library.individual_id,
                    pair.provenance,
                    String::from_utf8_lossy(&record.quality),
                    String::from_utf8_lossy(quality)
                ));
            }
        }
    }
    for key in observed.keys() {
        if !expected.contains_key(key) {
            mismatches.push(format!(
                "biblioteca {}: sobra la pareja {key}, que el oráculo no espera en esta biblioteca",
                library.individual_id
            ));
        }
    }
}
