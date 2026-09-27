//! Contrato y escritura del primer manifiesto de un evento técnico F3.
//!
//! La unidad no ensambla ni publica el árbol de una ejecución: recibe hechos
//! cerrados por sus productores. F3.5.3 será responsable de calcular estados,
//! integrarlo en el manifiesto general y publicar ambos documentos.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use serde_json::{Value, json};

use crate::application::EffectiveTechnicalEvent;
use crate::biological_contamination::{
    BiologicalContaminationManifestRecipient, contamination_target_pair_count,
};
use crate::demultiplexing::{DemultiplexingAggregate, DemultiplexingReport};
use crate::diagnostics::{Diagnostic, diagnostic_matches_catalog_version};
use crate::error_quality::{
    ErrorQualityReport, ErrorQualitySettings, MateErrorCounts, SequencingErrorsMode,
};
use crate::fastq_output::paired_fastq_relative_paths;
use crate::fastq_validation::F4CoverageReport;
use crate::identity::{EventId, IndividualId, LibraryId, ReadGroupId, ReferenceId};
use crate::results::{ArtifactRecord, UnitState, json_document_bytes, validate_relative_link};
use crate::technical_identity::TechnicalIdentityPlan;
use crate::technical_identity::TechnicalMate;

/// DTC-082 activa esta versión junto con effective-config v4 y manifest v4.
pub const EVENT_MANIFEST_SCHEMA_VERSION: &str = "2";
/// Versión que estrena la sección del tartamudeo.
///
/// Solo la publica un evento que **de verdad** abre la etapa `stutter`. Un
/// evento sin tartamudeo sigue publicando la v2, y eso es lo que mantiene byte
/// a byte los manifiestos anteriores a esta ampliación: subir la versión en
/// todos habría cambiado cada uno de ellos sin que nada cambiara dentro.
pub const EVENT_MANIFEST_SCHEMA_VERSION_V3: &str = "3";
/// F4 publica manifiestos de evento bajo el catálogo diagnóstico v3. F5
/// estrenará su propio documento v5 con el catálogo v4.
const EVENT_MANIFEST_DIAGNOSTIC_CATALOG_VERSION: &str = "3";
pub const EVENT_MANIFEST_KIND: &str = "gps_event_manifest";
pub const TECHNICAL_SEQUENCING_EVENT_TYPE: &str = "technical_sequencing";
/// Tecnología del único perfil técnico que F3 admite. El perfil conserva sus
/// parámetros efectivos; este identificador distingue la tecnología del evento.
pub const ILLUMINA_SHORT_READ_TECHNOLOGY_V1: &str = "illumina-short-read-pe-v1";
const REFERENCE_FASTA_ARTIFACT_KIND: &str = "reference_fasta";
const FASTQ_ARTIFACT_KIND: &str = "fastq";

/// Referencia común usada por el evento. Se enlaza, no se copia dentro del
/// directorio de evento.
#[derive(Clone, Copy, Debug)]
pub struct EventManifestReference<'a> {
    pub reference_id: &'a ReferenceId,
    pub fasta: &'a ArtifactRecord,
}

/// Estado de una salida propia del evento. Las salidas no completas no tienen
/// ruta, tamaño ni digest porque nunca se catalogan como artefactos válidos.
#[derive(Clone, Copy, Debug)]
pub struct EventManifestArtifact<'a> {
    pub id: &'a str,
    pub kind: &'a str,
    pub state: UnitState,
    pub individual_id: &'a IndividualId,
    pub library_id: &'a LibraryId,
    pub mate: TechnicalMate,
    pub artifact: Option<&'a ArtifactRecord>,
    pub diagnostic: Option<&'a Diagnostic>,
}

/// Resumen de cobertura final de una biblioteca, sin procedencia por lectura.
#[derive(Clone, Copy, Debug)]
pub struct EventManifestCoverage<'a> {
    pub individual_id: &'a IndividualId,
    pub library_id: &'a LibraryId,
    pub report: F4CoverageReport,
}

/// Resumen publicable de una pareja final. Sus campos son contadores de la
/// biblioteca completa y no permiten recuperar una pareja, donante o clave.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct EventManifestFinalFastqSummary<'a> {
    pub individual_id: &'a IndividualId,
    pub library_id: &'a LibraryId,
    pub pair_count: u64,
    pub read_count: u64,
    pub base_count: u64,
    pub own_pair_count: u64,
    pub contaminant_pair_count: u64,
    pub own_base_count: u64,
    pub contaminant_base_count: u64,
    pub n_base_count: u64,
}

/// Hechos efectivos de error y calidad por mate. Se agregan en toda la tanda:
/// no identifican una lectura, donante ni pareja concreta.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EventManifestMateErrorQuality {
    pub quality_shift: i16,
    pub error_scale: f64,
    pub opportunities: u64,
    pub errors: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EventManifestErrorQuality {
    pub requested_mode: SequencingErrorsMode,
    pub effective_mode: SequencingErrorsMode,
    pub qualities_represent_realized_error_probability: bool,
    pub r1: EventManifestMateErrorQuality,
    pub r2: EventManifestMateErrorQuality,
}

impl EventManifestErrorQuality {
    #[must_use]
    pub fn from_report(settings: ErrorQualitySettings, report: ErrorQualityReport) -> Self {
        Self {
            requested_mode: settings.sequencing_errors(),
            effective_mode: settings.sequencing_errors(),
            qualities_represent_realized_error_probability: settings
                .sequencing_errors()
                .qualities_represent_realized_error_probability(),
            r1: mate_error_quality(settings, TechnicalMate::R1, report.r1()),
            r2: mate_error_quality(settings, TechnicalMate::R2, report.r2()),
        }
    }

    #[must_use]
    pub fn neutral() -> Self {
        Self::from_report(
            ErrorQualitySettings::default(),
            ErrorQualityReport::default(),
        )
    }
}

fn mate_error_quality(
    settings: ErrorQualitySettings,
    mate: TechnicalMate,
    counts: MateErrorCounts,
) -> EventManifestMateErrorQuality {
    EventManifestMateErrorQuality {
        quality_shift: settings.quality_shift(mate),
        error_scale: settings.error_scale(mate),
        opportunities: counts.opportunities(),
        errors: counts.errors(),
    }
}

/// Hechos agregados del demultiplexado F4.5.1. Las rutas y los FASTQ son la
/// única salida durable; el informe no serializa instancias ni destinos por
/// pareja.
#[derive(Clone, Debug, PartialEq)]
pub struct EventManifestDemultiplexing {
    index_hopping_rate: f64,
    reassigned_pairs: u64,
    discarded_pairs: u64,
    unassigned_pairs: u64,
    reassigned_by_sample: Vec<DemultiplexingAggregate>,
    discarded_by_sample: Vec<DemultiplexingAggregate>,
    unassigned_by_sample: Vec<DemultiplexingAggregate>,
}

impl EventManifestDemultiplexing {
    #[must_use]
    pub fn from_report(report: DemultiplexingReport) -> Self {
        Self {
            index_hopping_rate: report.index_hopping_rate(),
            reassigned_pairs: report.reassigned_pairs(),
            discarded_pairs: report.discarded_pairs(),
            unassigned_pairs: report.unassigned_pairs(),
            reassigned_by_sample: report.reassigned_by_sample().to_vec(),
            discarded_by_sample: report.discarded_by_sample().to_vec(),
            unassigned_by_sample: report.unassigned_by_sample().to_vec(),
        }
    }

    #[must_use]
    pub fn disabled() -> Self {
        Self {
            index_hopping_rate: 0.0,
            reassigned_pairs: 0,
            discarded_pairs: 0,
            unassigned_pairs: 0,
            reassigned_by_sample: Vec::new(),
            discarded_by_sample: Vec::new(),
            unassigned_by_sample: Vec::new(),
        }
    }
}

fn demultiplexing_aggregates_document(entries: &[DemultiplexingAggregate]) -> Vec<Value> {
    entries
        .iter()
        .map(|entry| {
            let mut document = serde_json::Map::new();
            document.insert(
                "source_library_id".to_owned(),
                Value::String(entry.source_library_id().as_str().to_owned()),
            );
            if let Some(destination) = entry.destination_library_id() {
                document.insert(
                    "destination_library_id".to_owned(),
                    Value::String(destination.as_str().to_owned()),
                );
            }
            document.insert(
                "cause".to_owned(),
                Value::String(entry.cause().as_str().to_owned()),
            );
            document.insert("pair_count".to_owned(), Value::from(entry.pair_count()));
            Value::Object(document)
        })
        .collect()
}

/// Hechos del primer evento técnico que el escritor registra sin inferirlos de
/// rutas, nombres de fichero o mensajes.
pub struct EventManifest<'a> {
    pub event: &'a EffectiveTechnicalEvent,
    pub state: UnitState,
    pub identities: &'a TechnicalIdentityPlan,
    pub reference: EventManifestReference<'a>,
    pub artifacts: &'a [EventManifestArtifact<'a>],
    pub coverage: &'a [EventManifestCoverage<'a>],
    pub final_fastq_summary: &'a [EventManifestFinalFastqSummary<'a>],
    pub error_quality: EventManifestErrorQuality,
    pub demultiplexing: EventManifestDemultiplexing,
    /// Hecho agregado de F4.3. Nunca contiene origen, ID ni coordenadas de
    /// una pareja; esas piezas solo existen durante la transformación.
    pub chimera_rate: f64,
    pub realized_chimera_pairs: u64,
    /// Vista agregada del plan y la materialización F4.4. No contiene claves
    /// de pareja, IDs de instancia ni datos materializables individuales.
    pub biological_contamination: &'a [BiologicalContaminationManifestRecipient],
    /// Sección del tartamudeo. `None` cuando el evento no abre la etapa
    /// `stutter`, y entonces el documento no la menciona.
    pub stutter: Option<EventManifestStutter>,
    /// Sección de la amplificación. `None` con preparación sin PCR, y entonces
    /// `artifacts.pcr` conserva exactamente los bytes de antes de F4.9.
    pub pcr: Option<EventManifestPcr>,
    pub diagnostics: &'a [Diagnostic],
}

/// Lo que el manifiesto publica de la amplificación: sus parámetros efectivos
/// y sus **recuentos agregados**, más la cota de recursos calculada.
///
/// Nada por pareja. Lo más fino que aparece es un tamaño de familia, que es un
/// número y no una identidad: el manifiesto **no** publica qué parejas forman
/// familia, y un validador puede rederivarlo desde la semilla.
#[derive(Clone, Debug, PartialEq)]
pub struct EventManifestPcr {
    pub settings: crate::pcr_amplification::PcrSettings,
    pub report: crate::pcr_amplification::PcrReport,
    pub bound: crate::pcr_amplification::PcrResourceBound,
    /// Cobertura básica antes de artefactos, sumada sobre las bibliotecas. De
    /// ella y de la fracción realizada sale la cobertura de moléculas
    /// distintas.
    pub basic_wgs_coverage_x: f64,
}

impl EventManifestPcr {
    fn document(&self) -> Value {
        crate::pcr_amplification::manifest_document(
            &self.settings,
            &self.report,
            self.bound,
            self.basic_wgs_coverage_x,
        )
    }
}

/// Lo que el manifiesto publica del tartamudeo: sus parámetros efectivos con su
/// procedencia y sus **recuentos agregados**.
///
/// Nada por lectura ni por pareja. Lo más fino que aparece es un locus, que es
/// una entidad de la configuración, no de la salida.
#[derive(Clone, Debug, PartialEq)]
pub struct EventManifestStutter {
    pub settings: crate::stutter::StutterSettings,
    pub report: crate::stutter::StutterReport,
    /// Los loci declarados, en orden canónico, con las longitudes de motivo que
    /// cada uno usa. De aquí salen los parámetros efectivos por locus.
    pub loci: Vec<(String, Vec<u64>)>,
}

impl EventManifestStutter {
    fn document(&self) -> Value {
        let locus_parameters = self
            .loci
            .iter()
            .map(|(locus_id, motif_lengths)| {
                json!({
                    "locus_id": locus_id,
                    "base_rates": motif_lengths
                        .iter()
                        .map(|length| json!({
                            "motif_length_bases": length,
                            "rate": self
                                .settings
                                .base_rate_for(locus_id, *length)
                                .unwrap_or(0.0),
                        }))
                        .collect::<Vec<_>>(),
                    "rate_substituted": self.settings.locus_rates().contains_key(locus_id),
                })
            })
            .collect::<Vec<_>>();
        let steps_by_locus = self
            .report
            .steps_by_locus()
            .iter()
            .map(|(locus_id, steps)| {
                json!({
                    "locus_id": locus_id,
                    "steps": steps
                        .iter()
                        .map(|(units, count)| json!({"units": units, "count": count}))
                        .collect::<Vec<_>>(),
                })
            })
            .collect::<Vec<_>>();
        // El bloque de parámetros es **el mismo** que publica la configuración
        // efectiva; aquí solo se le añaden los loci y los recuentos realizados.
        let mut document = crate::stutter::effective_settings_document(&self.settings);
        if let Some(object) = document.as_object_mut() {
            // `locus_rates` es la sustitución declarada; `locus_parameters` es
            // el parámetro **efectivo** de cada locus. El manifiesto publica el
            // segundo, que es lo que de verdad se usó.
            object.remove("locus_rates");
            object.insert(
                "locus_parameters".to_owned(),
                Value::Array(locus_parameters),
            );
            object.insert(
                "evaluated_tracts".to_owned(),
                Value::from(self.report.evaluated_tracts()),
            );
            object.insert(
                "stuttered_tracts".to_owned(),
                Value::from(self.report.stuttered_tracts()),
            );
            object.insert("steps_by_locus".to_owned(), Value::Array(steps_by_locus));
        }
        document
    }
}

/// Ruta canónica del manifiesto dentro de la raíz de una ejecución.
#[must_use]
pub fn event_manifest_relative_path(event_id: &EventId) -> PathBuf {
    PathBuf::from("events")
        .join(event_id.as_str())
        .join("manifest.json")
}

/// Forma el documento JSON y verifica que los hechos recibidos sean coherentes
/// con la identidad técnica ya derivada.
pub fn event_manifest_document(manifest: &EventManifest<'_>) -> Result<Value, Diagnostic> {
    validate_writer_facts(manifest)?;
    let selected_individuals = manifest
        .identities
        .libraries()
        .map(|library| {
            let configured = manifest
                .event
                .libraries()
                .iter()
                .find(|configured| configured.individual_id() == library.individual_id())
                .expect("effective event invariants require one UDI per library");
            json!({
                "individual_id": library.individual_id().as_str(),
                "library_id": library.library_id().as_str(),
                "read_group_id": library.read_group_id().as_str(),
                "i7": configured.i7(),
                "i5": configured.i5(),
            })
        })
        .collect::<Vec<_>>();
    let mut artifacts = manifest.artifacts.to_vec();
    artifacts.sort_by_key(|artifact| {
        (
            artifact.individual_id.as_str(),
            artifact.library_id.as_str(),
            artifact.id,
        )
    });
    let artifacts = artifacts
        .into_iter()
        .map(|artifact| {
            json!({
                "id": artifact.id,
                "kind": artifact.kind,
                "state": artifact.state.as_str(),
                "individual_id": artifact.individual_id.as_str(),
                "library_id": artifact.library_id.as_str(),
                "mate": artifact.mate.number(),
                "artifact": artifact.artifact.map(ArtifactRecord::to_json),
                "diagnostic_code": artifact.diagnostic.map(Diagnostic::code),
            })
        })
        .collect::<Vec<_>>();
    let mut coverage = manifest.coverage.to_vec();
    coverage.sort_by_key(|entry| (entry.individual_id.as_str(), entry.library_id.as_str()));
    let coverage = coverage
        .into_iter()
        .map(|entry| {
            json!({
                "individual_id": entry.individual_id.as_str(),
                "library_id": entry.library_id.as_str(),
                "target_x": entry.report.target_x(),
                "basic_x": entry.report.basic_wgs_x(),
                "final_x": {
                    "own_x": entry.report.final_fastq_own_x(),
                    "contaminant_x": entry.report.final_fastq_contaminant_x(),
                },
                "denominator_bases": entry.report.denominator_bases(),
                "n_bases": entry.report.n_bases(),
            })
        })
        .collect::<Vec<_>>();
    let mut final_fastq_summary = manifest.final_fastq_summary.to_vec();
    final_fastq_summary
        .sort_by_key(|entry| (entry.individual_id.as_str(), entry.library_id.as_str()));
    let final_fastq_summary = final_fastq_summary
        .into_iter()
        .map(|entry| {
            json!({
                "individual_id": entry.individual_id.as_str(),
                "library_id": entry.library_id.as_str(),
                "pair_count": entry.pair_count,
                "read_count": entry.read_count,
                "base_count": entry.base_count,
                "own_pair_count": entry.own_pair_count,
                "contaminant_pair_count": entry.contaminant_pair_count,
                "own_base_count": entry.own_base_count,
                "contaminant_base_count": entry.contaminant_base_count,
                "n_base_count": entry.n_base_count,
            })
        })
        .collect::<Vec<_>>();
    let document = json!({
        "manifest_kind": EVENT_MANIFEST_KIND,
        "schema_version": event_manifest_schema_version(manifest),
        "event_id": manifest.event.event_id().as_str(),
        "event_type": TECHNICAL_SEQUENCING_EVENT_TYPE,
        "technology": ILLUMINA_SHORT_READ_TECHNOLOGY_V1,
        "rng_catalog": event_manifest_rng_catalog(manifest),
        "id_catalog": f4_id_catalog(),
        "state": manifest.state.as_str(),
        "profile": {
            "id": manifest.event.profile(),
            "resource_path": manifest_profile_resource_path(manifest),
            "sha256": manifest_profile_sha256(manifest),
            "preparation": manifest.event.preparation().as_str(),
            "read_length": manifest.event.read_length(),
            "quality": {
                "acgt_phred": manifest.event.acgt_phred(),
                "n_phred": manifest.event.n_phred(),
            },
        },
        "technical_batch_id": manifest.identities.batch_id().as_str(),
        "selected_individuals": selected_individuals,
        "event_subseed": manifest.identities.event_subseed().to_hex(),
        "reference": {
            "reference_id": manifest.reference.reference_id.as_str(),
            "fasta": manifest.reference.fasta.to_json(),
        },
        "coverage": coverage,
        "final_fastq_summary": final_fastq_summary,
        "error_quality": {
            "requested_mode": manifest.error_quality.requested_mode.as_str(),
            "effective_mode": manifest.error_quality.effective_mode.as_str(),
            "qualities_represent_realized_error_probability": manifest.error_quality.qualities_represent_realized_error_probability,
            "r1": error_quality_document(manifest.error_quality.r1),
            "r2": error_quality_document(manifest.error_quality.r2),
        },
        "demultiplexing": {
            "index_hopping_rate": manifest.demultiplexing.index_hopping_rate,
            "reassigned_pairs": manifest.demultiplexing.reassigned_pairs,
            "discarded_pairs": manifest.demultiplexing.discarded_pairs,
            "unassigned_pairs": manifest.demultiplexing.unassigned_pairs,
            "reassigned_by_sample": demultiplexing_aggregates_document(&manifest.demultiplexing.reassigned_by_sample),
            "discarded_by_sample": demultiplexing_aggregates_document(&manifest.demultiplexing.discarded_by_sample),
            "unassigned_by_sample": demultiplexing_aggregates_document(&manifest.demultiplexing.unassigned_by_sample),
            "clean_fastq_published": false,
            "individual_provenance_published": false,
        },
        "artifacts": {
            "chimera": {
                "rate": manifest.chimera_rate,
                "realized_pairs": manifest.realized_chimera_pairs,
            },
            "biological_contamination": manifest.biological_contamination.iter().map(|entry| json!({
                "recipient_id": entry.recipient_id().as_str(),
                "recipient_library_id": entry.recipient_library_id().as_str(),
                "percent_range": {
                    "minimum": entry.minimum_percent(),
                    "maximum": entry.maximum_percent(),
                },
                "sampled_percent": entry.sampled_percent(),
                "clean_pair_count": entry.clean_pair_count(),
                "target_pair_count": entry.target_pair_count(),
                "removed_own_pair_count": entry.removed_own_pair_count(),
                "inserted_contaminant_pair_count": entry.inserted_contaminant_pair_count(),
                "retained_own_pair_count": entry.retained_own_pair_count(),
                "final_pair_count": entry.final_pair_count(),
                "subseeds": {
                    "percent": entry.percent_subseed().to_hex(),
                    "donors": entry.donor_subseed().to_hex(),
                    "removal": entry.removal_subseed().to_hex(),
                },
                "donor_contributions": entry.donor_contributions().iter().map(|(donor, count)| json!({
                    "individual_id": donor.as_str(),
                    "pairs": count,
                })).collect::<Vec<_>>(),
            })).collect::<Vec<_>>(),
            "pcr": manifest.pcr.as_ref().map_or_else(
                || json!({"applicable": manifest.event.pcr_duplicates_applicable()}),
                EventManifestPcr::document,
            ),
            "optical_duplicates": {"realized_pairs": manifest.event.realized_optical_duplicate_pairs()},
            "gc": {
                "mode": manifest.pcr.as_ref().map_or("neutral", |pcr| pcr.settings.gc_mode()),
            },
            "adapters": {"enabled": false},
        },
        "outputs": artifacts,
        "diagnostics": canonical_diagnostics(manifest.diagnostics),
    });
    // La sección se **añade** en vez de escribirse siempre con un nulo: un
    // manifiesto de un evento sin tartamudeo tiene que salir con los mismos
    // bytes que antes de que el tartamudeo existiera, y una clave con `null`
    // ya sería un byte distinto.
    let mut document = document;
    if let Some(stutter) = manifest.stutter.as_ref() {
        let artifacts = document
            .get_mut("artifacts")
            .and_then(Value::as_object_mut)
            .ok_or_else(|| {
                event_error(
                    &event_manifest_relative_path(manifest.event.event_id()),
                    "el documento del manifiesto no conserva su sección de artefactos",
                )
            })?;
        artifacts.insert("stutter".to_owned(), stutter.document());
    }
    validate_event_manifest_bytes(
        &event_manifest_relative_path(manifest.event.event_id()),
        &json_document_bytes(&document)?,
    )?;
    Ok(document)
}

/// Serializa el manifiesto ya validado. F3.5.3 lo escribirá como artefacto
/// temporal del `ExecutionWorkspace` y lo finalizará con la política DTC-045.
pub fn event_manifest_bytes(manifest: &EventManifest<'_>) -> Result<Vec<u8>, Diagnostic> {
    json_document_bytes(&event_manifest_document(manifest)?)
}

/// Valida el esquema propio y las relaciones internas que JSON Schema no puede
/// expresar (identidades, orden estable y enlaces entre bibliotecas y FASTQ).
pub fn validate_event_manifest_bytes(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    // El esquema que el diagnóstico cita sale de la versión que el documento
    // **declara**, no de una constante. Denunciar siempre el de la versión 2
    // mandaba a mirar un fichero que no era el suyo en cuanto el evento abría
    // el tartamudeo o la amplificación. Se retiqueta **al salir**, que es el
    // único sitio donde la versión ya se conoce con seguridad y donde no hay
    // que tocar las ciento treinta y seis llamadas de dentro.
    validate_event_manifest_document(path, bytes).map_err(|error| {
        retag_event_manifest_schema(error, event_manifest_schema_file_for_bytes(bytes))
    })
}

/// El esquema que hay que citar, según lo que el documento declare.
///
/// Un documento cuya versión no se pueda leer —JSON ilegible, raíz que no es un
/// objeto, `schema_version` ausente— se queda con el de la versión 2, que es la
/// base: es la única respuesta honesta cuando el documento no dice de cuál
/// pretende ser.
fn event_manifest_schema_file_for(version: Option<&str>) -> &'static str {
    if version == Some(EVENT_MANIFEST_SCHEMA_VERSION_V3) {
        "manifest-event-v3.schema.json"
    } else {
        "manifest-event-v2.schema.json"
    }
}

/// Lee la versión declarada sin exigir que el documento sea válido.
///
/// Se hace **solo en el camino del error**, de modo que una validación que pasa
/// no vuelve a analizar el documento.
fn event_manifest_schema_file_for_bytes(bytes: &[u8]) -> &'static str {
    let declarada = serde_json::from_slice::<Value>(bytes)
        .ok()
        .and_then(|documento| {
            documento
                .get("schema_version")
                .and_then(Value::as_str)
                .map(str::to_owned)
        });
    event_manifest_schema_file_for(declarada.as_deref())
}

/// Reescribe el esquema citado por un diagnóstico de formato ya construido.
fn retag_event_manifest_schema(error: Diagnostic, schema: &'static str) -> Diagnostic {
    if error.code() != "GPS-FMT-006" {
        return error;
    }
    let context = error.context();
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with(
        "path",
        context
            .get("path")
            .cloned()
            .unwrap_or_else(|| "documento JSON".to_owned()),
    )
    .with("schema", schema)
    .with(
        "detail",
        context
            .get("detail")
            .cloned()
            .unwrap_or_else(|| "documento inválido".to_owned()),
    )
}

fn validate_event_manifest_document(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    let document: Value = serde_json::from_slice(bytes)
        .map_err(|error| event_error(path, format!("JSON no legible: {error}")))?;
    let root = required_object(&document, "raíz", path)?;
    validate_closed(root, ROOT_FIELDS, "raíz", path)?;
    require_text(root, "manifest_kind", path, EVENT_MANIFEST_KIND)?;
    // La versión del documento y el catálogo de ámbitos **van juntos**: el v4
    // solo existe porque la etapa `stutter` existe, y la sección del
    // tartamudeo solo cabe en la v3. Comprobarlos por separado dejaría pasar
    // un documento que dijera v3 y catálogo v3, o v2 con sección.
    let declared_version = required_string(root, "schema_version", path)?;
    let opens_new_stage = match declared_version {
        EVENT_MANIFEST_SCHEMA_VERSION => false,
        EVENT_MANIFEST_SCHEMA_VERSION_V3 => true,
        other => {
            return Err(event_error(
                path,
                format!(
                    "schema_version debe ser {EVENT_MANIFEST_SCHEMA_VERSION} o \
                     {EVENT_MANIFEST_SCHEMA_VERSION_V3}, se encontró {other}"
                ),
            ));
        }
    };
    // Las dos etapas nuevas se leen **del propio documento**, no de la versión:
    // la v3 cabe para el tartamudeo, para la amplificación o para las dos, y
    // cada una tiene su catálogo. Lo que la versión decide es si cabe alguna.
    let artifacts_object =
        required_object(required_value(root, "artifacts", path)?, "artifacts", path)?;
    let expects_stutter = artifacts_object.contains_key("stutter");
    let pcr_object = artifacts_object.get("pcr").and_then(Value::as_object);
    let expects_pcr = pcr_object.and_then(|pcr| pcr.get("applicable")) == Some(&Value::Bool(true));
    // La amplificación abre las decisiones de linaje **solo** cuando publica
    // su bloque de tartamudeo clonal, que es lo mismo que decir cuando el
    // escenario declara loci de repetición.
    let expects_lineage = expects_pcr && pcr_object.is_some_and(|pcr| pcr.contains_key("stutter"));
    if (expects_stutter || expects_pcr) != opens_new_stage {
        return Err(event_error(
            path,
            "la versión del documento no corresponde a las etapas que declara",
        ));
    }
    let event_id = parse_id::<EventId>(required_string(root, "event_id", path)?, "event", path)?;
    require_text(root, "event_type", path, TECHNICAL_SEQUENCING_EVENT_TYPE)?;
    require_text(root, "technology", path, ILLUMINA_SHORT_READ_TECHNOLOGY_V1)?;
    require_text(
        root,
        "rng_catalog",
        path,
        if expects_lineage {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V6
        } else if expects_pcr {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V5
        } else if expects_stutter {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
        } else {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
        },
    )?;
    validate_f4_id_catalog(required_value(root, "id_catalog", path)?, path)?;
    let state = parse_state(required_string(root, "state", path)?, path)?;
    let read_length = validate_profile(required_value(root, "profile", path)?, path)?;
    let _batch = parse_id::<crate::identity::BatchId>(
        required_string(root, "technical_batch_id", path)?,
        "batch",
        path,
    )?;
    let selected =
        validate_selected_individuals(required_value(root, "selected_individuals", path)?, path)?;
    validate_hex(
        required_string(root, "event_subseed", path)?,
        "event_subseed",
        path,
    )?;
    let reference = validate_reference(required_value(root, "reference", path)?, path)?;
    let coverage = validate_coverage(required_value(root, "coverage", path)?, &selected, path)?;
    let final_fastq_summary = validate_final_fastq_summary(
        required_value(root, "final_fastq_summary", path)?,
        &selected,
        read_length,
        path,
    )?;
    validate_error_quality(required_value(root, "error_quality", path)?, path)?;
    validate_final_fastq_summary_coverage(
        required_value(root, "coverage", path)?,
        required_value(root, "final_fastq_summary", path)?,
        path,
    )?;
    let final_pair_count = required_array(
        required_value(root, "final_fastq_summary", path)?,
        "final_fastq_summary",
        path,
    )?
    .iter()
    .try_fold(0_u64, |total, summary| {
        required_u64(
            required_object(summary, "final_fastq_summary[]", path)?,
            "pair_count",
            path,
        )?
        .checked_add(total)
        .ok_or_else(|| event_error(path, "pares FASTQ finales desbordan"))
    })?;
    validate_demultiplexing(
        required_value(root, "demultiplexing", path)?,
        &selected,
        final_pair_count,
        path,
    )?;
    let diagnostics = validate_diagnostics(required_value(root, "diagnostics", path)?, path)?;
    validate_technical_artifacts(
        required_value(root, "artifacts", path)?,
        expects_stutter,
        expects_pcr,
        path,
    )?;
    let artifact_states = validate_outputs(
        required_value(root, "outputs", path)?,
        &event_id,
        &selected,
        &reference,
        &diagnostics,
        path,
    )?;
    validate_event_state(state, &artifact_states, &diagnostics, path)?;
    if state == UnitState::Complete && final_fastq_summary != coverage {
        return Err(event_error(
            path,
            "resúmenes FASTQ y coberturas no cubren las mismas bibliotecas",
        ));
    }
    if state == UnitState::Complete && coverage.len() != selected.len() {
        return Err(event_error(
            path,
            "falta cobertura de una biblioteca completa",
        ));
    }
    if state == UnitState::Complete && final_fastq_summary.len() != selected.len() {
        return Err(event_error(
            path,
            "falta resumen FASTQ de una biblioteca completa",
        ));
    }
    if state != UnitState::Complete && !coverage.is_empty() {
        return Err(event_error(
            path,
            "un evento no completo no publica cobertura F4",
        ));
    }
    if state != UnitState::Complete && !final_fastq_summary.is_empty() {
        return Err(event_error(
            path,
            "un evento no completo no publica resúmenes FASTQ",
        ));
    }
    Ok(())
}

const ROOT_FIELDS: &[&str] = &[
    "manifest_kind",
    "schema_version",
    "event_id",
    "event_type",
    "technology",
    "rng_catalog",
    "id_catalog",
    "state",
    "profile",
    "technical_batch_id",
    "selected_individuals",
    "event_subseed",
    "reference",
    "coverage",
    "final_fastq_summary",
    "error_quality",
    "demultiplexing",
    "artifacts",
    "outputs",
    "diagnostics",
];

/// Ruta y huella del recurso del perfil que este evento usó.
fn manifest_profile_resource_path(manifest: &EventManifest<'_>) -> &'static str {
    match manifest.event.preparation() {
        crate::technical_profile::LibraryPreparation::PcrFree => {
            crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_RESOURCE
        }
        crate::technical_profile::LibraryPreparation::Pcr => {
            crate::technical_profile::GPS_HISEQX_PCR_PE150_RESOURCE
        }
    }
}

fn manifest_profile_sha256(manifest: &EventManifest<'_>) -> &'static str {
    match manifest.event.preparation() {
        crate::technical_profile::LibraryPreparation::PcrFree => {
            crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256
        }
        crate::technical_profile::LibraryPreparation::Pcr => {
            crate::technical_profile::GPS_HISEQX_PCR_PE150_SHA256
        }
    }
}

/// Versión del documento: la v3 cuando el evento abre la etapa `stutter` o la
/// etapa `amplification`. Las dos secciones nuevas solo caben en ella.
fn event_manifest_schema_version(manifest: &EventManifest<'_>) -> &'static str {
    if manifest.stutter.is_some() || manifest.pcr.is_some() {
        EVENT_MANIFEST_SCHEMA_VERSION_V3
    } else {
        EVENT_MANIFEST_SCHEMA_VERSION
    }
}

/// Catálogo de ámbitos aleatorios del que salen **las decisiones de este
/// evento**, no el más alto que el binario conoce.
fn event_manifest_rng_catalog(manifest: &EventManifest<'_>) -> &'static str {
    if let Some(pcr) = manifest.pcr.as_ref() {
        // El v6 solo cuando la amplificación abre de verdad las decisiones de
        // linaje. Una biblioteca amplificada sin loci de repetición sigue
        // publicando el v5, byte a byte como antes de F4.10.
        pcr.settings.rng_catalog()
    } else if manifest.stutter.is_some() {
        crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
    } else {
        crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
    }
}

fn validate_writer_facts(manifest: &EventManifest<'_>) -> Result<(), Diagnostic> {
    let path = event_manifest_relative_path(manifest.event.event_id());
    if manifest.error_quality.requested_mode != manifest.event.sequencing_errors()
        || manifest.error_quality.effective_mode != manifest.event.sequencing_errors()
        || manifest
            .error_quality
            .qualities_represent_realized_error_probability
            != manifest
                .event
                .sequencing_errors()
                .qualities_represent_realized_error_probability()
    {
        return Err(event_error(
            &path,
            "el modo de errores no coincide con la configuración efectiva",
        ));
    }
    if !has_fixed_error_quality_adjustments(manifest.error_quality.r1)
        || !has_fixed_error_quality_adjustments(manifest.error_quality.r2)
    {
        return Err(event_error(
            &path,
            "F4 no admite ajustes de error/calidad distintos de quality_shift=0 y error_scale=1",
        ));
    }
    if manifest.chimera_rate != manifest.event.chimera_rate()
        || !manifest.chimera_rate.is_finite()
        || !(0.0..=1.0).contains(&manifest.chimera_rate)
        || (manifest.chimera_rate == 0.0 && manifest.realized_chimera_pairs != 0)
    {
        return Err(event_error(
            &path,
            "la tasa de quimera no coincide con la configuración efectiva",
        ));
    }
    if manifest.identities.event_id() != manifest.event.event_id()
        || manifest.identities.reference_id() != manifest.reference.reference_id
    {
        return Err(event_error(
            &path,
            "el plan técnico no coincide con el evento o referencia",
        ));
    }
    if manifest.reference.fasta.kind != REFERENCE_FASTA_ARTIFACT_KIND {
        return Err(event_error(
            &path,
            "el enlace de referencia no apunta a un FASTA",
        ));
    }
    let configured = manifest
        .event
        .individual_ids()
        .iter()
        .collect::<BTreeSet<_>>();
    let planned = manifest
        .identities
        .libraries()
        .map(|library| library.individual_id())
        .collect::<BTreeSet<_>>();
    if configured.len() != manifest.event.individual_ids().len() || configured != planned {
        return Err(event_error(
            &path,
            "individuos técnicos y bibliotecas no coinciden",
        ));
    }
    let diagnostics = manifest
        .diagnostics
        .iter()
        .map(Diagnostic::code)
        .collect::<BTreeSet<_>>();
    if diagnostics.len() != manifest.diagnostics.len() {
        return Err(event_error(&path, "diagnósticos repetidos"));
    }
    let mut artifact_mates = BTreeSet::new();
    for artifact in manifest.artifacts {
        let library = manifest
            .identities
            .library(artifact.individual_id)
            .ok_or_else(|| {
                event_error(&path, "un artefacto refiere un individuo no seleccionado")
            })?;
        if library.library_id() != artifact.library_id {
            return Err(event_error(
                &path,
                "un artefacto refiere una biblioteca ajena",
            ));
        }
        if !artifact_mates.insert((artifact.individual_id.as_str(), artifact.mate)) {
            return Err(event_error(&path, "una biblioteca repite un mate FASTQ"));
        }
        match artifact.state {
            UnitState::Complete => {
                let record = artifact.artifact.ok_or_else(|| {
                    event_error(&path, "un artefacto completo no tiene identidad física")
                })?;
                if artifact.diagnostic.is_some()
                    || record.id != artifact.id
                    || record.kind != artifact.kind
                    || record.path
                        != paired_path(
                            manifest.event.event_id().as_str(),
                            artifact.individual_id.as_str(),
                            artifact.mate,
                        )
                {
                    return Err(event_error(&path, "artefacto completo incoherente"));
                }
            }
            UnitState::Incomplete | UnitState::Failed => {
                let diagnostic = artifact
                    .diagnostic
                    .ok_or_else(|| event_error(&path, "artefacto no completo sin diagnóstico"))?;
                if artifact.artifact.is_some() || !diagnostics.contains(diagnostic.code()) {
                    return Err(event_error(&path, "artefacto no completo incoherente"));
                }
            }
        }
    }
    if artifact_mates.len() != planned.len() * 2
        || !planned.iter().all(|individual| {
            TechnicalMate::all()
                .into_iter()
                .all(|mate| artifact_mates.contains(&(individual.as_str(), mate)))
        })
    {
        return Err(event_error(
            &path,
            "falta el resultado FASTQ de una biblioteca seleccionada",
        ));
    }
    if manifest.state == UnitState::Complete {
        let coverage_libraries = manifest
            .coverage
            .iter()
            .map(|entry| (entry.individual_id.as_str(), entry.library_id.as_str()))
            .collect::<BTreeSet<_>>();
        let planned_libraries = manifest
            .identities
            .libraries()
            .map(|library| {
                (
                    library.individual_id().as_str(),
                    library.library_id().as_str(),
                )
            })
            .collect::<BTreeSet<_>>();
        if coverage_libraries.len() != manifest.coverage.len()
            || coverage_libraries != planned_libraries
        {
            return Err(event_error(
                &path,
                "cobertura por biblioteca incompleta o repetida",
            ));
        }
        let summary_libraries = manifest
            .final_fastq_summary
            .iter()
            .map(|entry| (entry.individual_id.as_str(), entry.library_id.as_str()))
            .collect::<BTreeSet<_>>();
        if summary_libraries.len() != manifest.final_fastq_summary.len()
            || summary_libraries != planned_libraries
        {
            return Err(event_error(
                &path,
                "resumen FASTQ por biblioteca incompleto o repetido",
            ));
        }
    } else if !manifest.coverage.is_empty() {
        return Err(event_error(
            &path,
            "un evento no completo no publica cobertura F4",
        ));
    }
    Ok(())
}

fn f4_id_catalog() -> Value {
    json!({
        "version": crate::identity::ID_CATALOG_VERSION,
        "derivations": {
            "v1": ["execution", "reference", "individual", "event", "batch", "library", "read_group", "template"],
            "v2": ["PreHopInstanceId", "TechnicalPairId"],
        }
    })
}

fn validate_f4_id_catalog(value: &Value, path: &Path) -> Result<(), Diagnostic> {
    let catalog = required_object(value, "id_catalog", path)?;
    validate_closed(catalog, &["version", "derivations"], "id_catalog", path)?;
    require_text(
        catalog,
        "version",
        path,
        crate::identity::ID_CATALOG_VERSION,
    )?;
    let derivations = required_object(
        required_value(catalog, "derivations", path)?,
        "id_catalog.derivations",
        path,
    )?;
    validate_closed(derivations, &["v1", "v2"], "id_catalog.derivations", path)?;
    let expected_v1 = [
        "execution",
        "reference",
        "individual",
        "event",
        "batch",
        "library",
        "read_group",
        "template",
    ];
    let expected_v2 = ["PreHopInstanceId", "TechnicalPairId"];
    for (field, expected) in [
        ("v1", expected_v1.as_slice()),
        ("v2", expected_v2.as_slice()),
    ] {
        let values = required_array(required_value(derivations, field, path)?, field, path)?;
        if values.iter().filter_map(Value::as_str).collect::<Vec<_>>() != expected {
            return Err(event_error(path, "id_catalog incompatible"));
        }
    }
    Ok(())
}

fn validate_profile(value: &Value, path: &Path) -> Result<u64, Diagnostic> {
    let profile = required_object(value, "profile", path)?;
    validate_closed(
        profile,
        &[
            "id",
            "resource_path",
            "sha256",
            "preparation",
            "read_length",
            "quality",
        ],
        "profile",
        path,
    )?;
    // **La terna va junta.** Comprobar el identificador, la ruta, la huella y
    // la preparación por separado dejaría pasar un documento que publicara el
    // identificador de un perfil con la huella del otro.
    let (resource_path, sha256, preparation) = match required_string(profile, "id", path)? {
        crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_V1 => (
            crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_RESOURCE,
            crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256,
            crate::technical_profile::LibraryPreparation::PcrFree.as_str(),
        ),
        crate::technical_profile::GPS_HISEQX_PCR_PE150_V1 => (
            crate::technical_profile::GPS_HISEQX_PCR_PE150_RESOURCE,
            crate::technical_profile::GPS_HISEQX_PCR_PE150_SHA256,
            crate::technical_profile::LibraryPreparation::Pcr.as_str(),
        ),
        _ => return Err(event_error(path, "id incompatible")),
    };
    require_text(profile, "resource_path", path, resource_path)?;
    require_text(profile, "sha256", path, sha256)?;
    require_text(profile, "preparation", path, preparation)?;
    let read_length = profile
        .get("read_length")
        .and_then(Value::as_u64)
        .filter(|value| (1..=150).contains(value))
        .ok_or_else(|| event_error(path, "longitud de lectura de perfil inválida"))?;
    let quality = required_object(
        required_value(profile, "quality", path)?,
        "profile.quality",
        path,
    )?;
    validate_closed(quality, &["acgt_phred", "n_phred"], "profile.quality", path)?;
    if quality.get("acgt_phred").and_then(Value::as_u64)
        != Some(u64::from(
            crate::application::ILLUMINA_SE_MINIMAL_ACGT_PHRED,
        ))
        || quality.get("n_phred").and_then(Value::as_u64)
            != Some(u64::from(crate::application::ILLUMINA_SE_MINIMAL_N_PHRED))
    {
        return Err(event_error(path, "calidad de perfil inválida"));
    }
    Ok(read_length)
}

fn validate_selected_individuals(
    value: &Value,
    path: &Path,
) -> Result<BTreeMap<String, String>, Diagnostic> {
    let values = required_array(value, "selected_individuals", path)?;
    if values.is_empty() {
        return Err(event_error(path, "el evento no selecciona individuos"));
    }
    let mut selected = BTreeMap::new();
    let mut previous = None::<String>;
    let mut read_groups = BTreeSet::new();
    for value in values {
        let item = required_object(value, "selected_individuals[]", path)?;
        validate_closed(
            item,
            &["individual_id", "library_id", "read_group_id", "i7", "i5"],
            "selected_individuals[]",
            path,
        )?;
        let individual = parse_id::<IndividualId>(
            required_string(item, "individual_id", path)?,
            "individual",
            path,
        )?;
        let library =
            parse_id::<LibraryId>(required_string(item, "library_id", path)?, "library", path)?;
        let read_group = parse_id::<ReadGroupId>(
            required_string(item, "read_group_id", path)?,
            "read_group",
            path,
        )?;
        for field in ["i7", "i5"] {
            let udi = required_string(item, field, path)?;
            if udi.len() != 8
                || !udi
                    .bytes()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            {
                return Err(event_error(path, "UDI de biblioteca inválido"));
            }
        }
        if previous
            .as_deref()
            .is_some_and(|prior| prior >= individual.as_str())
            || selected
                .insert(individual.as_str().to_owned(), library.as_str().to_owned())
                .is_some()
            || !read_groups.insert(read_group.as_str().to_owned())
        {
            return Err(event_error(
                path,
                "individuos, bibliotecas o grupos no son únicos y ordenados",
            ));
        }
        previous = Some(individual.as_str().to_owned());
    }
    Ok(selected)
}

fn validate_coverage(
    value: &Value,
    selected: &BTreeMap<String, String>,
    path: &Path,
) -> Result<BTreeSet<(String, String)>, Diagnostic> {
    let values = required_array(value, "coverage", path)?;
    let mut coverage_libraries = BTreeSet::new();
    let mut previous = None::<(String, String)>;
    for value in values {
        let entry = required_object(value, "coverage[]", path)?;
        validate_closed(
            entry,
            &[
                "individual_id",
                "library_id",
                "target_x",
                "basic_x",
                "final_x",
                "denominator_bases",
                "n_bases",
            ],
            "coverage[]",
            path,
        )?;
        let individual = required_string(entry, "individual_id", path)?;
        let library = required_string(entry, "library_id", path)?;
        let ordering = (individual.to_owned(), library.to_owned());
        if selected.get(individual).map(String::as_str) != Some(library)
            || previous.as_ref().is_some_and(|prior| prior >= &ordering)
            || !coverage_libraries.insert(ordering.clone())
        {
            return Err(event_error(
                path,
                "coberturas no son únicas, ordenadas o compatibles con las bibliotecas",
            ));
        }
        previous = Some(ordering);
        for field in ["target_x", "basic_x"] {
            validate_nonnegative_finite(entry, field, path)?;
        }
        let final_x = required_object(
            required_value(entry, "final_x", path)?,
            "coverage[].final_x",
            path,
        )?;
        validate_closed(
            final_x,
            &["own_x", "contaminant_x"],
            "coverage[].final_x",
            path,
        )?;
        for field in ["own_x", "contaminant_x"] {
            validate_nonnegative_finite(final_x, field, path)?;
        }
        let _denominator = entry
            .get("denominator_bases")
            .and_then(Value::as_u64)
            .filter(|value| *value > 0)
            .ok_or_else(|| event_error(path, "denominador de cobertura inválido"))?;
        let _n_bases = entry
            .get("n_bases")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento N de cobertura inválido"))?;
    }
    Ok(coverage_libraries)
}

fn validate_final_fastq_summary(
    value: &Value,
    selected: &BTreeMap<String, String>,
    read_length: u64,
    path: &Path,
) -> Result<BTreeSet<(String, String)>, Diagnostic> {
    let values = required_array(value, "final_fastq_summary", path)?;
    let mut libraries = BTreeSet::new();
    let mut previous = None::<(String, String)>;
    for value in values {
        let entry = required_object(value, "final_fastq_summary[]", path)?;
        validate_closed(
            entry,
            &[
                "individual_id",
                "library_id",
                "pair_count",
                "read_count",
                "base_count",
                "own_pair_count",
                "contaminant_pair_count",
                "own_base_count",
                "contaminant_base_count",
                "n_base_count",
            ],
            "final_fastq_summary[]",
            path,
        )?;
        let individual = required_string(entry, "individual_id", path)?;
        let library = required_string(entry, "library_id", path)?;
        let ordering = (individual.to_owned(), library.to_owned());
        if selected.get(individual).map(String::as_str) != Some(library)
            || previous.as_ref().is_some_and(|prior| prior >= &ordering)
            || !libraries.insert(ordering.clone())
        {
            return Err(event_error(
                path,
                "resúmenes FASTQ no son únicos, ordenados o compatibles",
            ));
        }
        previous = Some(ordering);
        let pair_count = required_u64(entry, "pair_count", path)?;
        let read_count = required_u64(entry, "read_count", path)?;
        let base_count = required_u64(entry, "base_count", path)?;
        let own_pairs = required_u64(entry, "own_pair_count", path)?;
        let contaminant_pairs = required_u64(entry, "contaminant_pair_count", path)?;
        let own_bases = required_u64(entry, "own_base_count", path)?;
        let contaminant_bases = required_u64(entry, "contaminant_base_count", path)?;
        let n_bases = required_u64(entry, "n_base_count", path)?;
        let total_bases = own_bases
            .checked_add(contaminant_bases)
            .ok_or_else(|| event_error(path, "contadores agregados FASTQ desbordan"))?;
        let expected_read_count = pair_count
            .checked_mul(2)
            .ok_or_else(|| event_error(path, "contador de lecturas FASTQ desborda"))?;
        let expected_base_count = expected_read_count
            .checked_mul(read_length)
            .ok_or_else(|| event_error(path, "contador de bases FASTQ desborda"))?;
        if own_pairs.checked_add(contaminant_pairs) != Some(pair_count)
            || read_count != expected_read_count
            || base_count != expected_base_count
            || total_bases != expected_base_count
            || n_bases > base_count
        {
            return Err(event_error(path, "contadores agregados FASTQ incoherentes"));
        }
    }
    Ok(libraries)
}

fn validate_final_fastq_summary_coverage(
    coverage: &Value,
    summary: &Value,
    path: &Path,
) -> Result<(), Diagnostic> {
    let coverage = required_array(coverage, "coverage", path)?;
    let summary = required_array(summary, "final_fastq_summary", path)?;
    if summary.is_empty() {
        return Ok(());
    }
    let mut coverage_by_library = BTreeMap::new();
    for entry in coverage {
        let entry = required_object(entry, "coverage[]", path)?;
        let final_x = required_object(
            required_value(entry, "final_x", path)?,
            "coverage[].final_x",
            path,
        )?;
        coverage_by_library.insert(
            (
                required_string(entry, "individual_id", path)?,
                required_string(entry, "library_id", path)?,
            ),
            (
                required_u64(entry, "denominator_bases", path)?,
                required_u64(entry, "n_bases", path)?,
                final_x
                    .get("own_x")
                    .and_then(Value::as_f64)
                    .ok_or_else(|| event_error(path, "cobertura propia inválida"))?,
                final_x
                    .get("contaminant_x")
                    .and_then(Value::as_f64)
                    .ok_or_else(|| event_error(path, "cobertura contaminante inválida"))?,
            ),
        );
    }
    for entry in summary {
        let entry = required_object(entry, "final_fastq_summary[]", path)?;
        let key = (
            required_string(entry, "individual_id", path)?,
            required_string(entry, "library_id", path)?,
        );
        let Some((denominator, n_bases, own_x, contaminant_x)) = coverage_by_library.get(&key)
        else {
            return Err(event_error(
                path,
                "resumen FASTQ no tiene cobertura asociada",
            ));
        };
        let own_bases = required_u64(entry, "own_base_count", path)?;
        let contaminant_bases = required_u64(entry, "contaminant_base_count", path)?;
        if required_u64(entry, "n_base_count", path)? != *n_bases
            || (own_bases as f64 / *denominator as f64).to_bits() != own_x.to_bits()
            || (contaminant_bases as f64 / *denominator as f64).to_bits() != contaminant_x.to_bits()
        {
            return Err(event_error(
                path,
                "resumen FASTQ no coincide con su cobertura final",
            ));
        }
    }
    Ok(())
}

fn validate_demultiplexing(
    value: &Value,
    selected: &BTreeMap<String, String>,
    final_pair_count: u64,
    path: &Path,
) -> Result<(), Diagnostic> {
    let report = required_object(value, "demultiplexing", path)?;
    validate_closed(
        report,
        &[
            "index_hopping_rate",
            "reassigned_pairs",
            "discarded_pairs",
            "unassigned_pairs",
            "reassigned_by_sample",
            "discarded_by_sample",
            "unassigned_by_sample",
            "clean_fastq_published",
            "individual_provenance_published",
        ],
        "demultiplexing",
        path,
    )?;
    if !matches!(report.get("index_hopping_rate").and_then(Value::as_f64), Some(value) if value.is_finite() && (0.0..=1.0).contains(&value))
        || report
            .get("reassigned_pairs")
            .and_then(Value::as_u64)
            .is_none()
        || report
            .get("discarded_pairs")
            .and_then(Value::as_u64)
            .is_none()
        || report
            .get("unassigned_pairs")
            .and_then(Value::as_u64)
            .is_none()
        || report.get("clean_fastq_published") != Some(&Value::Bool(false))
        || report.get("individual_provenance_published") != Some(&Value::Bool(false))
    {
        return Err(event_error(
            path,
            "resumen de demultiplexado inválido o publica procedencia",
        ));
    }
    let discarded = required_u64(report, "discarded_pairs", path)?;
    let unassigned = required_u64(report, "unassigned_pairs", path)?;
    let reassigned = required_u64(report, "reassigned_pairs", path)?;
    let index_hopping_rate = report
        .get("index_hopping_rate")
        .and_then(Value::as_f64)
        .expect("validated index_hopping_rate is finite");
    if unassigned > discarded {
        return Err(event_error(
            path,
            "parejas sin asignar exceden los descartes",
        ));
    }
    if reassigned > final_pair_count {
        return Err(event_error(
            path,
            "reasignaciones exceden las parejas FASTQ finales",
        ));
    }
    if index_hopping_rate == 0.0 && reassigned != 0 {
        return Err(event_error(
            path,
            "un índice de salto nulo no puede reasignar parejas",
        ));
    }
    let selected_libraries = selected
        .values()
        .map(String::as_str)
        .collect::<BTreeSet<_>>();
    validate_demultiplexing_aggregates(
        report,
        "reassigned_by_sample",
        "index_hopping",
        true,
        reassigned,
        &selected_libraries,
        path,
    )?;
    let discarded_by_sample = validate_demultiplexing_aggregates(
        report,
        "discarded_by_sample",
        "unrecognized_udi",
        false,
        discarded,
        &selected_libraries,
        path,
    )?;
    let unassigned_by_sample = validate_demultiplexing_aggregates(
        report,
        "unassigned_by_sample",
        "unrecognized_udi",
        false,
        unassigned,
        &selected_libraries,
        path,
    )?;
    if discarded_by_sample != unassigned_by_sample {
        return Err(event_error(
            path,
            "descartes y parejas sin asignar no coinciden por biblioteca",
        ));
    }
    Ok(())
}

fn validate_demultiplexing_aggregates(
    report: &serde_json::Map<String, Value>,
    field: &str,
    expected_cause: &str,
    requires_destination: bool,
    expected_total: u64,
    selected_libraries: &BTreeSet<&str>,
    path: &Path,
) -> Result<BTreeMap<String, u64>, Diagnostic> {
    let entries = report
        .get(field)
        .and_then(Value::as_array)
        .ok_or_else(|| event_error(path, "agregados de demultiplexado inválidos"))?;
    let mut total = 0_u64;
    let mut previous = None;
    let mut counts_by_source = BTreeMap::new();
    for entry in entries {
        let entry = entry
            .as_object()
            .ok_or_else(|| event_error(path, "agregado de demultiplexado inválido"))?;
        let expected_fields: &[&str] = if requires_destination {
            &[
                "source_library_id",
                "destination_library_id",
                "cause",
                "pair_count",
            ]
        } else {
            &["source_library_id", "cause", "pair_count"]
        };
        validate_closed(entry, expected_fields, field, path)?;
        let source = entry
            .get("source_library_id")
            .and_then(Value::as_str)
            .filter(|value| valid_manifest_id(value) && selected_libraries.contains(value))
            .ok_or_else(|| event_error(path, "origen de demultiplexado inválido"))?;
        let destination = entry.get("destination_library_id").and_then(Value::as_str);
        if requires_destination {
            let destination = destination
                .filter(|value| valid_manifest_id(value) && selected_libraries.contains(value))
                .ok_or_else(|| event_error(path, "destino de demultiplexado inválido"))?;
            if destination == source {
                return Err(event_error(path, "reasignación de demultiplexado circular"));
            }
        } else if destination.is_some() {
            return Err(event_error(path, "descarte de demultiplexado con destino"));
        }
        if entry.get("cause").and_then(Value::as_str) != Some(expected_cause) {
            return Err(event_error(path, "causa de demultiplexado inválida"));
        }
        let count = required_u64(entry, "pair_count", path)?;
        if count == 0 {
            return Err(event_error(path, "agregado de demultiplexado vacío"));
        }
        total = total
            .checked_add(count)
            .ok_or_else(|| event_error(path, "recuento de demultiplexado desbordado"))?;
        counts_by_source.insert(source.to_owned(), count);
        let key = (source, expected_cause, destination.unwrap_or(""));
        if previous.is_some_and(|previous_key| previous_key >= key) {
            return Err(event_error(
                path,
                "agregados de demultiplexado desordenados",
            ));
        }
        previous = Some(key);
    }
    if total != expected_total {
        return Err(event_error(
            path,
            "agregados de demultiplexado no suman el total",
        ));
    }
    Ok(counts_by_source)
}

fn valid_manifest_id(value: &str) -> bool {
    let bytes = value.as_bytes();
    (1..=64).contains(&bytes.len())
        && matches!(bytes.first(), Some(first) if first.is_ascii_alphanumeric())
        && bytes
            .iter()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
}

fn error_quality_document(entry: EventManifestMateErrorQuality) -> Value {
    json!({
        "quality_shift": entry.quality_shift,
        "error_scale": entry.error_scale,
        "opportunities": entry.opportunities,
        "errors": entry.errors,
    })
}

fn has_fixed_error_quality_adjustments(entry: EventManifestMateErrorQuality) -> bool {
    entry.quality_shift == 0 && entry.error_scale == 1.0
}

fn validate_error_quality(value: &Value, path: &Path) -> Result<(), Diagnostic> {
    let error_quality = required_object(value, "error_quality", path)?;
    validate_closed(
        error_quality,
        &[
            "requested_mode",
            "effective_mode",
            "qualities_represent_realized_error_probability",
            "r1",
            "r2",
        ],
        "error_quality",
        path,
    )?;
    let requested_mode = parse_sequencing_errors_mode(
        required_string(error_quality, "requested_mode", path)?,
        path,
    )?;
    let effective_mode = parse_sequencing_errors_mode(
        required_string(error_quality, "effective_mode", path)?,
        path,
    )?;
    let qualities_represent_realized_error_probability = error_quality
        .get("qualities_represent_realized_error_probability")
        .and_then(Value::as_bool)
        .ok_or_else(|| event_error(path, "marca de calidad de error inválida"))?;
    if requested_mode != effective_mode
        || qualities_represent_realized_error_probability
            != effective_mode.qualities_represent_realized_error_probability()
    {
        return Err(event_error(
            path,
            "modo de error y marca de calidad incoherentes",
        ));
    }
    let mut realized_errors = [0_u64; 2];
    for mate in ["r1", "r2"] {
        let entry = required_object(
            required_value(error_quality, mate, path)?,
            "error_quality mate",
            path,
        )?;
        validate_closed(
            entry,
            &["quality_shift", "error_scale", "opportunities", "errors"],
            "error_quality mate",
            path,
        )?;
        if entry.get("quality_shift").and_then(Value::as_i64) != Some(0)
            || entry.get("error_scale").and_then(Value::as_f64) != Some(1.0)
        {
            return Err(event_error(
                path,
                "F4 no admite ajustes de error/calidad distintos de quality_shift=0 y error_scale=1",
            ));
        }
        let opportunities = required_u64(entry, "opportunities", path)?;
        let mate_errors = required_u64(entry, "errors", path)?;
        if mate_errors > opportunities {
            return Err(event_error(
                path,
                "errores realizados exceden oportunidades",
            ));
        }
        realized_errors[usize::from(mate == "r2")] = mate_errors;
    }
    if effective_mode == SequencingErrorsMode::Disabled && realized_errors != [0, 0] {
        return Err(event_error(
            path,
            "sequencing_errors disabled exige cero errores por mate",
        ));
    }
    Ok(())
}

fn parse_sequencing_errors_mode(
    value: &str,
    path: &Path,
) -> Result<SequencingErrorsMode, Diagnostic> {
    match value {
        "enabled" => Ok(SequencingErrorsMode::Enabled),
        "disabled" => Ok(SequencingErrorsMode::Disabled),
        _ => Err(event_error(
            path,
            "modo de errores de secuenciación inválido",
        )),
    }
}

fn required_u64(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<u64, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_u64)
        .ok_or_else(|| event_error(path, format!("{field} no es entero sin signo")))
}

fn validate_nonnegative_finite(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if !matches!(object.get(field).and_then(Value::as_f64), Some(value) if value.is_finite() && value >= 0.0)
    {
        return Err(event_error(path, format!("{field} de cobertura inválida")));
    }
    Ok(())
}

fn validate_reference(value: &Value, path: &Path) -> Result<ArtifactRecord, Diagnostic> {
    let reference = required_object(value, "reference", path)?;
    validate_closed(reference, &["reference_id", "fasta"], "reference", path)?;
    let _ = parse_id::<ReferenceId>(
        required_string(reference, "reference_id", path)?,
        "reference",
        path,
    )?;
    let fasta = validate_artifact(required_value(reference, "fasta", path)?, path)?;
    if fasta.kind != REFERENCE_FASTA_ARTIFACT_KIND {
        return Err(event_error(path, "la referencia no enlaza un FASTA"));
    }
    Ok(fasta)
}

fn validate_diagnostics(value: &Value, path: &Path) -> Result<BTreeSet<String>, Diagnostic> {
    let values = required_array(value, "diagnostics", path)?;
    let mut codes = BTreeSet::new();
    let mut previous = None::<String>;
    for value in values {
        let diagnostic = required_object(value, "diagnostics[]", path)?;
        validate_closed(diagnostic, &["code", "message"], "diagnostics[]", path)?;
        let code = required_string(diagnostic, "code", path)?;
        let message = required_string(diagnostic, "message", path)?;
        if previous.as_deref().is_some_and(|prior| prior >= code)
            || !diagnostic_matches_catalog_version(
                EVENT_MANIFEST_DIAGNOSTIC_CATALOG_VERSION,
                code,
                message,
            )
            || !codes.insert(code.to_owned())
        {
            return Err(event_error(
                path,
                "diagnóstico desconocido, alterado o repetido",
            ));
        }
        previous = Some(code.to_owned());
    }
    Ok(codes)
}

/// Valida los únicos artefactos técnicos agregados activos en F4.3. El
/// contrato es deliberadamente cerrado para que no se filtren listas de
/// parejas, identificadores o coordenadas al manifiesto durable.
fn validate_technical_artifacts(
    value: &Value,
    expects_stutter: bool,
    expects_pcr: bool,
    path: &Path,
) -> Result<(), Diagnostic> {
    let artifacts = required_object(value, "artifacts", path)?;
    validate_closed(
        artifacts,
        if expects_stutter {
            &[
                "chimera",
                "biological_contamination",
                "pcr",
                "optical_duplicates",
                "gc",
                "adapters",
                "stutter",
            ]
        } else {
            &[
                "chimera",
                "biological_contamination",
                "pcr",
                "optical_duplicates",
                "gc",
                "adapters",
            ]
        },
        "artifacts",
        path,
    )?;
    if expects_stutter {
        validate_stutter_manifest(required_value(artifacts, "stutter", path)?, path)?;
    }
    let chimera = required_object(
        required_value(artifacts, "chimera", path)?,
        "artifacts.chimera",
        path,
    )?;
    validate_biological_contamination_manifest(
        required_value(artifacts, "biological_contamination", path)?,
        path,
    )?;
    validate_closed(
        chimera,
        &["rate", "realized_pairs"],
        "artifacts.chimera",
        path,
    )?;
    let rate = chimera.get("rate").and_then(Value::as_f64);
    let realized_pairs = chimera.get("realized_pairs").and_then(Value::as_u64);
    if !matches!(rate, Some(rate) if rate.is_finite() && (0.0..=1.0).contains(&rate))
        || realized_pairs.is_none()
        || (rate == Some(0.0) && realized_pairs != Some(0))
    {
        return Err(event_error(path, "artefacto de quimera inválido"));
    }
    let pcr = required_object(
        required_value(artifacts, "pcr", path)?,
        "artifacts.pcr",
        path,
    )?;
    let optical = required_object(
        required_value(artifacts, "optical_duplicates", path)?,
        "artifacts.optical_duplicates",
        path,
    )?;
    let gc = required_object(required_value(artifacts, "gc", path)?, "artifacts.gc", path)?;
    let adapters = required_object(
        required_value(artifacts, "adapters", path)?,
        "artifacts.adapters",
        path,
    )?;
    if expects_pcr {
        validate_pcr_manifest(pcr, path)?;
    } else {
        validate_closed(pcr, &["applicable"], "artifacts.pcr", path)?;
        if pcr.get("applicable") != Some(&Value::Bool(false)) {
            return Err(event_error(path, "artefactos técnicos fijos inválidos"));
        }
    }
    validate_closed(
        optical,
        &["realized_pairs"],
        "artifacts.optical_duplicates",
        path,
    )?;
    validate_closed(gc, &["mode"], "artifacts.gc", path)?;
    validate_closed(adapters, &["enabled"], "artifacts.adapters", path)?;
    // El modo `curve` solo existe con amplificación: publicarlo sin ella
    // afirmaría un sesgo por composición que nadie aplicó.
    let gc_mode_valid = match gc.get("mode").and_then(Value::as_str) {
        Some("neutral") => true,
        Some("curve") => expects_pcr,
        _ => false,
    };
    if optical.get("realized_pairs").and_then(Value::as_u64) != Some(0)
        || !gc_mode_valid
        || adapters.get("enabled") != Some(&Value::Bool(false))
    {
        return Err(event_error(path, "artefactos técnicos fijos inválidos"));
    }
    Ok(())
}

/// Comprueba la sección de amplificación del manifiesto de evento.
///
/// Las relaciones que se comprueban aquí son las que JSON Schema no puede
/// expresar y las que de verdad discriminan un manifiesto inventado:
///
/// - la suma del histograma **ponderada por su índice** vale exactamente el
///   número de parejas, que es lo que impide que una celda cuente familias en
///   vez de parejas;
/// - las parejas duplicadas son `parejas − familias`;
/// - la fracción realizada es `1 − familias/parejas`, de modo que no puede ser
///   una copia de la declarada;
/// - la derivación de `m` viaja con el método y las iteraciones.
fn validate_pcr_manifest(
    pcr: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    // Los tres campos del linaje **van juntos o no van**: publicar la
    // derivación sin los recuentos, o al revés, dejaría un documento en el que
    // no se puede saber si el tartamudeo clonal llegó a ocurrir. Y solo
    // aparecen cuando el escenario declara loci de repetición, de modo que una
    // biblioteca amplificada sin microsatélites tiene que seguir publicando
    // **exactamente** los mismos campos que antes de F4.10.
    let lineage_fields = [
        "stutter",
        "lineage_slips_by_locus",
        "lineage_inherited_slips_by_locus",
        "lineage_applied_steps_by_locus",
    ];
    let present = lineage_fields
        .iter()
        .filter(|field| pcr.contains_key(**field))
        .count();
    if present != 0 && present != lineage_fields.len() {
        return Err(event_error(
            path,
            "la sección de tartamudeo clonal y sus tres recuentos de linaje van juntos o no van",
        ));
    }
    let opens_lineage = present == lineage_fields.len();
    let mut fields = vec![
        "applicable",
        "algorithm",
        "rng_catalog",
        "cycles",
        "efficiency",
        "early_cycles",
        "declared_duplicate_fraction",
        "expected_duplicate_fraction",
        "mean_pairs_per_molecule",
        "duplicate_fraction_inversion",
        "max_family_size",
        "gc_max_attempts",
        "gc_relative_efficiency",
        "efficiency_classes",
        "families",
        "pairs",
        "duplicate_pairs",
        "realized_duplicate_fraction",
        "family_size_histogram",
        "families_truncated_at_end",
        "gc_attempts_total",
        "gc_attempts_exhausted_families",
        "unique_molecule_coverage_x",
        "resource_bound",
    ];
    if opens_lineage {
        fields.extend_from_slice(&lineage_fields);
    }
    validate_closed(pcr, &fields, "artifacts.pcr", path)?;
    require_text(
        pcr,
        "algorithm",
        path,
        crate::pcr_amplification::PCR_AMPLIFICATION_ALGORITHM,
    )?;
    require_text(
        pcr,
        "rng_catalog",
        path,
        if opens_lineage {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V6
        } else {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V5
        },
    )?;
    if opens_lineage {
        validate_lineage_counts(pcr, path, required_u64(pcr, "pairs", path)?)?;
    }
    let inversion = required_object(
        required_value(pcr, "duplicate_fraction_inversion", path)?,
        "artifacts.pcr.duplicate_fraction_inversion",
        path,
    )?;
    validate_closed(
        inversion,
        &["method", "iterations"],
        "artifacts.pcr.duplicate_fraction_inversion",
        path,
    )?;
    require_text(
        inversion,
        "method",
        path,
        crate::pcr_amplification::PCR_DUPLICATE_FRACTION_INVERSION_METHOD,
    )?;
    let bound = required_object(
        required_value(pcr, "resource_bound", path)?,
        "artifacts.pcr.resource_bound",
        path,
    )?;
    validate_closed(
        bound,
        &[
            "operation_units",
            "memory_bytes",
            "regeneration_passes",
            "measured",
        ],
        "artifacts.pcr.resource_bound",
        path,
    )?;
    // La cota es un valor **calculado sin ejecutar**, no una medida. Decirlo en
    // el documento evita que alguien la lea como consumo observado.
    if bound.get("measured") != Some(&Value::Bool(false)) {
        return Err(event_error(
            path,
            "la cota de recursos se publica calculada, no medida",
        ));
    }
    let families = required_u64(pcr, "families", path)?;
    let pairs = required_u64(pcr, "pairs", path)?;
    let duplicate_pairs = required_u64(pcr, "duplicate_pairs", path)?;
    if families > pairs || duplicate_pairs != pairs - families {
        return Err(event_error(
            path,
            "los recuentos de familias y duplicados no cuadran",
        ));
    }
    let histogram = required_array(
        required_value(pcr, "family_size_histogram", path)?,
        "artifacts.pcr.family_size_histogram",
        path,
    )?;
    let mut histogram_families = 0_u64;
    let mut histogram_pairs = 0_u64;
    let mut previous: Option<u64> = None;
    for entry in histogram {
        let entry = required_object(entry, "artifacts.pcr.family_size_histogram[]", path)?;
        validate_closed(
            entry,
            &["family_size", "families"],
            "artifacts.pcr.family_size_histogram[]",
            path,
        )?;
        let size = required_u64(entry, "family_size", path)?;
        let count = required_u64(entry, "families", path)?;
        if size == 0 || count == 0 || previous.is_some_and(|previous| previous >= size) {
            return Err(event_error(
                path,
                "el histograma de tamaños de familia no está en orden estricto",
            ));
        }
        previous = Some(size);
        histogram_families = histogram_families
            .checked_add(count)
            .ok_or_else(|| event_error(path, "el histograma de familias desborda"))?;
        histogram_pairs = histogram_pairs
            .checked_add(
                size.checked_mul(count)
                    .ok_or_else(|| event_error(path, "el histograma de familias desborda"))?,
            )
            .ok_or_else(|| event_error(path, "el histograma de familias desborda"))?;
    }
    if histogram_families != families || histogram_pairs != pairs {
        return Err(event_error(
            path,
            "el histograma de tamaños de familia no suma las familias y parejas publicadas",
        ));
    }
    let realized = pcr
        .get("realized_duplicate_fraction")
        .and_then(Value::as_f64)
        .ok_or_else(|| event_error(path, "falta la fracción realizada de duplicados"))?;
    let expected_realized = if pairs == 0 {
        0.0
    } else {
        1.0 - (families as f64) / (pairs as f64)
    };
    if realized != expected_realized {
        return Err(event_error(
            path,
            "la fracción realizada de duplicados no sale de los recuentos publicados",
        ));
    }
    let exhausted = required_u64(pcr, "gc_attempts_exhausted_families", path)?;
    let attempts = required_u64(pcr, "gc_attempts_total", path)?;
    let truncated = required_u64(pcr, "families_truncated_at_end", path)?;
    if exhausted > families || attempts < families || truncated > families {
        return Err(event_error(
            path,
            "los recuentos de intentos de composición o de familias truncadas no cuadran",
        ));
    }
    Ok(())
}

/// Comprueba los tres recuentos de linaje del manifiesto.
///
/// # La relación que de verdad discrimina
///
/// Los heredados son un **subconjunto** de los resbalones, y los dos listados
/// nombran **los mismos loci**. Un recuento de heredados que contara todos los
/// resbalones daría las dos listas iguales locus a locus, cosa que aquí no se
/// prohíbe —puede ocurrir legítimamente si todo se heredó—, pero uno que los
/// contara **de más** rompe la desigualdad y muere aquí.
fn validate_lineage_counts(
    pcr: &serde_json::Map<String, Value>,
    path: &Path,
    pairs: u64,
) -> Result<(), Diagnostic> {
    let slips = locus_slip_counts(pcr, "lineage_slips_by_locus", "slips", path)?;
    let inherited = locus_slip_counts(pcr, "lineage_inherited_slips_by_locus", "slips", path)?;
    let applied = locus_slip_counts(pcr, "lineage_applied_steps_by_locus", "steps", path)?;
    if slips.keys().ne(inherited.keys()) || slips.keys().ne(applied.keys()) {
        return Err(event_error(
            path,
            "los tres recuentos de linaje no nombran los mismos loci",
        ));
    }
    for (locus_id, count) in &inherited {
        if slips.get(locus_id).copied().unwrap_or(0) < *count {
            return Err(event_error(
                path,
                "hay más resbalones heredados que resbalones en un locus",
            ));
        }
    }
    // Un paso neto escrito es **uno por molécula y locus**, de modo que no
    // puede haber más que parejas. Sin esta cota, un recuento que sumara
    // resbalones en vez de pasos aplicados pasaría desapercibido en cuanto una
    // sola copia resbalara dos veces.
    for count in applied.values() {
        if *count > pairs {
            return Err(event_error(
                path,
                "hay más pasos de linaje escritos que parejas en un locus",
            ));
        }
    }
    Ok(())
}

/// Lee una lista de recuentos por locus, comprobando que cada fila está
/// cerrada y que ningún locus aparece dos veces.
fn locus_slip_counts(
    pcr: &serde_json::Map<String, Value>,
    field: &str,
    count_field: &str,
    path: &Path,
) -> Result<BTreeMap<String, u64>, Diagnostic> {
    let entries = required_array(
        required_value(pcr, field, path)?,
        &format!("artifacts.pcr.{field}"),
        path,
    )?;
    let mut counts = BTreeMap::new();
    for entry in entries {
        let entry = required_object(entry, &format!("artifacts.pcr.{field}[]"), path)?;
        validate_closed(
            entry,
            &["locus_id", count_field],
            &format!("artifacts.pcr.{field}[]"),
            path,
        )?;
        let locus_id = required_string(entry, "locus_id", path)?.to_owned();
        let slips = required_u64(entry, count_field, path)?;
        if counts.insert(locus_id, slips).is_some() {
            return Err(event_error(
                path,
                "un locus aparece dos veces en un recuento de linaje",
            ));
        }
    }
    Ok(counts)
}

fn validate_biological_contamination_manifest(
    value: &Value,
    path: &Path,
) -> Result<(), Diagnostic> {
    let records = required_array(value, "artifacts.biological_contamination", path)?;
    let mut recipients = BTreeSet::new();
    for value in records {
        let record = required_object(value, "artifacts.biological_contamination[]", path)?;
        validate_closed(
            record,
            &[
                "recipient_id",
                "recipient_library_id",
                "percent_range",
                "sampled_percent",
                "clean_pair_count",
                "target_pair_count",
                "removed_own_pair_count",
                "inserted_contaminant_pair_count",
                "retained_own_pair_count",
                "final_pair_count",
                "subseeds",
                "donor_contributions",
            ],
            "artifacts.biological_contamination[]",
            path,
        )?;
        let recipient = required_string(record, "recipient_id", path)?;
        let _ = parse_id::<IndividualId>(recipient, "individual", path)?;
        let _ = parse_id::<LibraryId>(
            required_string(record, "recipient_library_id", path)?,
            "library",
            path,
        )?;
        if !recipients.insert(recipient.to_owned()) {
            return Err(event_error(
                path,
                "receptora de contaminación biológica repetida",
            ));
        }
        let range = required_object(
            required_value(record, "percent_range", path)?,
            "percent_range",
            path,
        )?;
        validate_closed(range, &["minimum", "maximum"], "percent_range", path)?;
        let minimum = range.get("minimum").and_then(Value::as_f64);
        let maximum = range.get("maximum").and_then(Value::as_f64);
        let sampled = record.get("sampled_percent").and_then(Value::as_f64);
        if !matches!((minimum, maximum, sampled), (Some(minimum), Some(maximum), Some(sampled)) if minimum.is_finite() && maximum.is_finite() && sampled.is_finite() && (0.0..=100.0).contains(&minimum) && (0.0..=100.0).contains(&maximum) && minimum <= sampled && sampled <= maximum)
        {
            return Err(event_error(
                path,
                "porcentaje de contaminación biológica inválido",
            ));
        }
        let sampled = sampled
            .ok_or_else(|| event_error(path, "porcentaje de contaminación biológica ausente"))?;
        let clean = record
            .get("clean_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento limpio de contaminación inválido"))?;
        let target = record
            .get("target_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "objetivo de contaminación inválido"))?;
        let removed = record
            .get("removed_own_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento retirado de contaminación inválido"))?;
        let inserted = record
            .get("inserted_contaminant_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento insertado de contaminación inválido"))?;
        let retained = record
            .get("retained_own_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento retenido de contaminación inválido"))?;
        let final_pairs = record
            .get("final_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| event_error(path, "recuento final de contaminación inválido"))?;
        if target > clean {
            return Err(event_error(
                path,
                "objetivo de contaminación supera parejas limpias",
            ));
        }
        let expected_target = contamination_target_pair_count(clean, sampled)
            .map_err(|_| event_error(path, "objetivo de contaminación no es recomputable"))?;
        if target != expected_target {
            return Err(event_error(
                path,
                "objetivo de contaminación no aplica el redondeo normativo",
            ));
        }
        if removed != target
            || inserted != target
            || retained.checked_add(removed) != Some(clean)
            || retained.checked_add(inserted) != Some(final_pairs)
            || final_pairs != clean
        {
            return Err(event_error(
                path,
                "los contadores realizados de contaminación no conservan las parejas",
            ));
        }
        let subseeds =
            required_object(required_value(record, "subseeds", path)?, "subseeds", path)?;
        validate_closed(
            subseeds,
            &["percent", "donors", "removal"],
            "subseeds",
            path,
        )?;
        for value in ["percent", "donors", "removal"] {
            validate_subseed_hex(required_string(subseeds, value, path)?, value, path)?;
        }
        let contributions = required_array(
            required_value(record, "donor_contributions", path)?,
            "donor_contributions",
            path,
        )?;
        let mut donors = BTreeSet::new();
        let mut total = 0_u64;
        for value in contributions {
            let contribution = required_object(value, "donor_contributions[]", path)?;
            validate_closed(
                contribution,
                &["individual_id", "pairs"],
                "donor_contributions[]",
                path,
            )?;
            let donor = required_string(contribution, "individual_id", path)?;
            let _ = parse_id::<IndividualId>(donor, "individual", path)?;
            let pairs = contribution
                .get("pairs")
                .and_then(Value::as_u64)
                .ok_or_else(|| event_error(path, "aporte de donante inválido"))?;
            if donor == recipient || pairs == 0 || !donors.insert(donor.to_owned()) {
                return Err(event_error(path, "aporte de donante inválido"));
            }
            total = total
                .checked_add(pairs)
                .ok_or_else(|| event_error(path, "desbordamiento de aportes"))?;
        }
        if total != target {
            return Err(event_error(
                path,
                "los aportes no suman el objetivo de contaminación",
            ));
        }
    }
    Ok(())
}

fn validate_subseed_hex(value: &str, field: &str, path: &Path) -> Result<(), Diagnostic> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        return Err(event_error(
            path,
            format!("{field} no es una subsemilla hexadecimal de 64 caracteres"),
        ));
    }
    Ok(())
}

fn validate_outputs(
    value: &Value,
    event_id: &EventId,
    selected: &BTreeMap<String, String>,
    reference: &ArtifactRecord,
    diagnostics: &BTreeSet<String>,
    path: &Path,
) -> Result<Vec<UnitState>, Diagnostic> {
    let values = required_array(value, "outputs", path)?;
    if values.is_empty() {
        return Err(event_error(path, "el evento no tiene resultados FASTQ"));
    }
    let mut ids = BTreeSet::new();
    let mut paths = BTreeSet::new();
    let mut artifact_mates = BTreeSet::new();
    let mut states = Vec::with_capacity(values.len());
    let mut previous = None::<(String, String, u8, String)>;
    for value in values {
        let artifact = required_object(value, "outputs[]", path)?;
        validate_closed(
            artifact,
            &[
                "id",
                "kind",
                "state",
                "individual_id",
                "library_id",
                "mate",
                "artifact",
                "diagnostic_code",
            ],
            "outputs[]",
            path,
        )?;
        let id = required_string(artifact, "id", path)?;
        let kind = required_string(artifact, "kind", path)?;
        if kind != FASTQ_ARTIFACT_KIND {
            return Err(event_error(path, "F3 solo admite artefactos FASTQ"));
        }
        let state = parse_state(required_string(artifact, "state", path)?, path)?;
        let individual = required_string(artifact, "individual_id", path)?;
        let library = required_string(artifact, "library_id", path)?;
        let mate = parse_mate(artifact, path)?;
        let ordering = (
            individual.to_owned(),
            library.to_owned(),
            mate.number(),
            id.to_owned(),
        );
        if selected.get(individual).map(String::as_str) != Some(library)
            || !ids.insert(id.to_owned())
            || !artifact_mates.insert((individual.to_owned(), mate))
            || previous.as_ref().is_some_and(|prior| prior >= &ordering)
        {
            return Err(event_error(
                path,
                "artefactos FASTQ no son únicos, completos por biblioteca u ordenados",
            ));
        }
        previous = Some(ordering);
        match state {
            UnitState::Complete => {
                let record =
                    validate_artifact(required_non_null(artifact, "artifact", path)?, path)?;
                if record.id != id
                    || record.kind != kind
                    || record.path != paired_path(event_id.as_str(), individual, mate)
                    || !paths.insert(record.path.clone())
                    || record.id == reference.id
                    || record.path == reference.path
                    || artifact.get("diagnostic_code") != Some(&Value::Null)
                {
                    return Err(event_error(path, "artefacto FASTQ completo incoherente"));
                }
            }
            UnitState::Incomplete | UnitState::Failed => {
                if artifact.get("artifact") != Some(&Value::Null) {
                    return Err(event_error(
                        path,
                        "artefacto no completo no puede tener identidad física",
                    ));
                }
                let code = required_string(artifact, "diagnostic_code", path)?;
                if !diagnostics.contains(code) {
                    return Err(event_error(path, "diagnóstico de artefacto ausente"));
                }
            }
        }
        states.push(state);
    }
    if artifact_mates.len() != selected.len() * 2
        || !selected.keys().all(|individual| {
            TechnicalMate::all()
                .into_iter()
                .all(|mate| artifact_mates.contains(&(individual.clone(), mate)))
        })
    {
        return Err(event_error(
            path,
            "falta el resultado FASTQ de una biblioteca seleccionada",
        ));
    }
    Ok(states)
}

fn paired_path(event_id: &str, individual_id: &str, mate: TechnicalMate) -> String {
    let (r1, r2) = paired_fastq_relative_paths(event_id, individual_id);
    match mate {
        TechnicalMate::R1 => r1,
        TechnicalMate::R2 => r2,
    }
    .to_string_lossy()
    .into_owned()
}

fn parse_mate(
    artifact: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<TechnicalMate, Diagnostic> {
    match artifact.get("mate").and_then(Value::as_u64) {
        Some(1) => Ok(TechnicalMate::R1),
        Some(2) => Ok(TechnicalMate::R2),
        _ => Err(event_error(path, "mate FASTQ inválido")),
    }
}

fn validate_event_state(
    state: UnitState,
    artifact_states: &[UnitState],
    diagnostics: &BTreeSet<String>,
    path: &Path,
) -> Result<(), Diagnostic> {
    match state {
        UnitState::Complete
            if !diagnostics.is_empty()
                || artifact_states
                    .iter()
                    .any(|artifact_state| *artifact_state != UnitState::Complete) =>
        {
            Err(event_error(
                path,
                "un evento completo exige FASTQ completos y ningún diagnóstico",
            ))
        }
        UnitState::Incomplete | UnitState::Failed
            if diagnostics.is_empty()
                || artifact_states
                    .iter()
                    .all(|artifact_state| *artifact_state == UnitState::Complete) =>
        {
            Err(event_error(
                path,
                "un evento no completo exige diagnóstico y un FASTQ no completo",
            ))
        }
        _ => Ok(()),
    }
}

fn canonical_diagnostics(diagnostics: &[Diagnostic]) -> Vec<Value> {
    let mut diagnostics = diagnostics.iter().collect::<Vec<_>>();
    diagnostics.sort_by_key(|diagnostic| diagnostic.code());
    diagnostics
        .into_iter()
        .map(|diagnostic| {
            json!({
                "code": diagnostic.code(), "message": diagnostic.message(),
            })
        })
        .collect()
}

fn validate_artifact(value: &Value, path: &Path) -> Result<ArtifactRecord, Diagnostic> {
    let artifact = required_object(value, "artifact", path)?;
    validate_closed(
        artifact,
        &["id", "kind", "path", "size_bytes", "checksum"],
        "artifact",
        path,
    )?;
    let id = required_string(artifact, "id", path)?;
    let kind = required_string(artifact, "kind", path)?;
    let artifact_path = required_string(artifact, "path", path)?;
    validate_relative_link(artifact_path)
        .map_err(|_| event_error(path, "ruta de artefacto no canónica"))?;
    let size_bytes = artifact
        .get("size_bytes")
        .and_then(Value::as_u64)
        .ok_or_else(|| event_error(path, "tamaño de artefacto inválido"))?;
    let checksum = required_object(
        required_value(artifact, "checksum", path)?,
        "checksum",
        path,
    )?;
    validate_closed(checksum, &["algorithm", "digest"], "checksum", path)?;
    require_text(checksum, "algorithm", path, "sha256")?;
    let digest = required_string(checksum, "digest", path)?;
    validate_hex(digest, "checksum.digest", path)?;
    Ok(ArtifactRecord {
        id: id.to_owned(),
        kind: kind.to_owned(),
        path: artifact_path.to_owned(),
        size_bytes,
        digest: digest.to_owned(),
    })
}

fn parse_state(value: &str, path: &Path) -> Result<UnitState, Diagnostic> {
    match value {
        "incomplete" => Ok(UnitState::Incomplete),
        "complete" => Ok(UnitState::Complete),
        "failed" => Ok(UnitState::Failed),
        _ => Err(event_error(path, "estado inválido")),
    }
}

fn parse_id<T>(value: &str, category: &str, path: &Path) -> Result<T, Diagnostic>
where
    T: ParseTypedId,
{
    T::parse_id(value).map_err(|_| event_error(path, format!("ID {category} inválido")))
}

trait ParseTypedId: Sized {
    fn parse_id(value: &str) -> Result<Self, Diagnostic>;
}
macro_rules! typed_id_parser { ($($type:ty),+ $(,)?) => { $(impl ParseTypedId for $type { fn parse_id(value: &str) -> Result<Self, Diagnostic> { Self::parse(value, "event_manifest") } })+ }; }
typed_id_parser!(
    EventId,
    IndividualId,
    LibraryId,
    ReadGroupId,
    ReferenceId,
    crate::identity::BatchId,
    crate::identity::TemplateId
);

fn validate_hex(value: &str, field: &str, path: &Path) -> Result<(), Diagnostic> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'a'..=b'f'))
    {
        return Err(event_error(
            path,
            format!("{field} debe ser SHA-256 hexadecimal minúsculo"),
        ));
    }
    Ok(())
}

/// Comprueba la sección del tartamudeo: lo que JSON Schema no puede expresar.
///
/// La comprobación que importa es la **identidad contable**: la suma de la
/// distribución de pasos por locus tiene que valer exactamente lo mismo que el
/// número de moléculas que tartamudearon. Un manifiesto que publicara el número
/// **esperado** en vez del realizado la rompe, y esta es la única comprobación
/// del documento que lo caza.
fn validate_stutter_manifest(value: &Value, path: &Path) -> Result<(), Diagnostic> {
    let stutter = required_object(value, "artifacts.stutter", path)?;
    validate_closed(
        stutter,
        &[
            "enabled",
            "algorithm",
            "rng_catalog",
            "calibration",
            "motif_length_rates",
            "tract_length_dependence",
            "run_length_coefficient",
            "step_weights",
            "locus_parameters",
            "evaluated_tracts",
            "stuttered_tracts",
            "steps_by_locus",
        ],
        "artifacts.stutter",
        path,
    )?;
    require_text(
        stutter,
        "algorithm",
        path,
        crate::stutter::STUTTER_ALGORITHM,
    )?;
    require_text(
        stutter,
        "rng_catalog",
        path,
        crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4,
    )?;
    if stutter.get("enabled").and_then(Value::as_bool) != Some(true) {
        return Err(event_error(
            path,
            "un manifiesto que publica la sección de tartamudeo la publica activa",
        ));
    }
    let coefficient = required_object(
        required_value(stutter, "run_length_coefficient", path)?,
        "artifacts.stutter.run_length_coefficient",
        path,
    )?;
    validate_closed(
        coefficient,
        &["value", "calibrated"],
        "artifacts.stutter.run_length_coefficient",
        path,
    )?;
    // **La marca de no calibrado no es decorativa.** Sin ella, un coeficiente
    // que nadie ha medido pasaría por medido.
    if coefficient.get("calibrated").and_then(Value::as_bool)
        != Some(crate::stutter::STUTTER_RUN_LENGTH_COEFFICIENT_CALIBRATED)
    {
        return Err(event_error(
            path,
            "el coeficiente de la tirada tiene que publicarse como no calibrado",
        ));
    }
    let calibration = required_object(
        required_value(stutter, "calibration", path)?,
        "artifacts.stutter.calibration",
        path,
    )?;
    validate_closed(
        calibration,
        &[
            "source",
            "species",
            "technology",
            "column",
            "known_mismatch",
        ],
        "artifacts.stutter.calibration",
        path,
    )?;
    require_text(
        calibration,
        "known_mismatch",
        path,
        crate::stutter::STUTTER_KNOWN_CALIBRATION_MISMATCH,
    )?;
    // **La marca de tosquedad de la curva tampoco es decorativa.** Una curva
    // leída a ojo de una figura, publicada sin ella, se leería como una
    // medida; y sin `extrapolated` en falso nadie sabría que por encima del
    // tramo medido se mantiene el extremo.
    let dependence = required_object(
        required_value(stutter, "tract_length_dependence", path)?,
        "artifacts.stutter.tract_length_dependence",
        path,
    )?;
    validate_closed(
        dependence,
        &[
            "source",
            "coarse",
            "coarse_reason",
            "extrapolated",
            "by_motif_length",
        ],
        "artifacts.stutter.tract_length_dependence",
        path,
    )?;
    require_text(
        dependence,
        "source",
        path,
        crate::stutter::STUTTER_TRACT_LENGTH_SOURCE,
    )?;
    require_text(
        dependence,
        "coarse_reason",
        path,
        crate::stutter::STUTTER_TRACT_LENGTH_COARSE_REASON,
    )?;
    if dependence.get("coarse").and_then(Value::as_bool)
        != Some(crate::stutter::STUTTER_TRACT_LENGTH_COARSE)
    {
        return Err(event_error(
            path,
            "la dependencia con la longitud del tracto se publica marcada como tosca",
        ));
    }
    if dependence.get("extrapolated").and_then(Value::as_bool) != Some(false) {
        return Err(event_error(
            path,
            "fuera del tramo medido se mantiene el extremo y no se extrapola",
        ));
    }
    let by_motif_length = required_array(
        required_value(dependence, "by_motif_length", path)?,
        "artifacts.stutter.tract_length_dependence.by_motif_length",
        path,
    )?;
    if by_motif_length.len() != crate::stutter::STUTTER_MOTIF_LENGTH_MAX {
        return Err(event_error(
            path,
            "la tabla de longitudes de referencia no cubre el tramo calibrado",
        ));
    }
    for (index, entry) in by_motif_length.iter().enumerate() {
        let entry = required_object(
            entry,
            "artifacts.stutter.tract_length_dependence.by_motif_length",
            path,
        )?;
        validate_closed(
            entry,
            &[
                "motif_length_bases",
                "measured",
                "reference_tract_bases",
                "measured_units_min",
                "measured_units_max",
                "factors",
            ],
            "artifacts.stutter.tract_length_dependence.by_motif_length",
            path,
        )?;
        let expected = crate::stutter::STUTTER_TRACT_LENGTH_REFERENCE_BASES[index];
        let published = entry.get("reference_tract_bases").and_then(|value| {
            if value.is_null() {
                None
            } else {
                value.as_u64()
            }
        });
        if published != expected
            || entry.get("measured").and_then(Value::as_bool) != Some(expected.is_some())
        {
            return Err(event_error(
                path,
                "una longitud de tracto de referencia no es la del algoritmo",
            ));
        }
        let motif_length = (index as u64) + 1;
        // **El tramo medido se compara contra el algoritmo, no solo se exige
        // que las claves estén.** Faltaba, y era una asimetría con la guarda
        // de la configuración efectiva: un manifiesto podía declarar un tramo
        // medido que no era el suyo y nadie lo distinguía de uno correcto.
        let range = crate::stutter::tract_length_measured_units(motif_length);
        if entry.get("measured_units_min").and_then(Value::as_u64) != range.map(|(min, _)| min)
            || entry.get("measured_units_max").and_then(Value::as_u64) != range.map(|(_, max)| max)
        {
            return Err(event_error(
                path,
                "el tramo medido publicado no es el del algoritmo",
            ));
        }
        let factors = required_array(
            required_value(entry, "factors", path)?,
            "artifacts.stutter.tract_length_dependence.by_motif_length.factors",
            path,
        )?;
        if factors.len() != range.map_or(0, |(min, max)| (max - min + 1) as usize) {
            return Err(event_error(
                path,
                "la tabla de factores no cubre el tramo medido",
            ));
        }
        for (offset, factor) in factors.iter().enumerate() {
            let factor = required_object(
                factor,
                "artifacts.stutter.tract_length_dependence.by_motif_length.factors",
                path,
            )?;
            validate_closed(
                factor,
                &["units", "factor"],
                "artifacts.stutter.tract_length_dependence.by_motif_length.factors",
                path,
            )?;
            let units = range.map_or(0, |(min, _)| min) + offset as u64;
            if factor.get("units").and_then(Value::as_u64) != Some(units)
                || factor.get("factor").and_then(Value::as_f64)
                    != Some(crate::stutter::tract_length_factor(motif_length, units))
            {
                return Err(event_error(
                    path,
                    "un factor de la curva no es el del algoritmo",
                ));
            }
        }
    }
    let step_weights = required_array(
        required_value(stutter, "step_weights", path)?,
        "artifacts.stutter.step_weights",
        path,
    )?;
    let declared_steps = step_weights
        .iter()
        .map(|entry| {
            let entry = required_object(entry, "artifacts.stutter.step_weights", path)?;
            validate_closed(
                entry,
                &["units", "weight"],
                "artifacts.stutter.step_weights",
                path,
            )?;
            entry
                .get("units")
                .and_then(Value::as_i64)
                .ok_or_else(|| event_error(path, "un paso publicado no declara sus unidades"))
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;
    if declared_steps != crate::stutter::STUTTER_ADMITTED_STEPS {
        return Err(event_error(
            path,
            "los pesos publicados no cubren exactamente la estructura del algoritmo",
        ));
    }
    let stuttered = required_u64(stutter, "stuttered_tracts", path)?;
    let evaluated = required_u64(stutter, "evaluated_tracts", path)?;
    if stuttered > evaluated {
        return Err(event_error(
            path,
            "no pueden tartamudear más moléculas de las evaluadas",
        ));
    }
    let mut total = 0_u64;
    for entry in required_array(
        required_value(stutter, "steps_by_locus", path)?,
        "artifacts.stutter.steps_by_locus",
        path,
    )? {
        let entry = required_object(entry, "artifacts.stutter.steps_by_locus", path)?;
        validate_closed(
            entry,
            &["locus_id", "steps"],
            "artifacts.stutter.steps_by_locus",
            path,
        )?;
        for step in required_array(
            required_value(entry, "steps", path)?,
            "artifacts.stutter.steps_by_locus.steps",
            path,
        )? {
            let step = required_object(step, "artifacts.stutter.steps_by_locus.steps", path)?;
            validate_closed(
                step,
                &["units", "count"],
                "artifacts.stutter.steps_by_locus.steps",
                path,
            )?;
            let units = step
                .get("units")
                .and_then(Value::as_i64)
                .ok_or_else(|| event_error(path, "un paso realizado no declara sus unidades"))?;
            if !crate::stutter::STUTTER_ADMITTED_STEPS.contains(&units) {
                return Err(event_error(
                    path,
                    format!("el paso realizado {units} no pertenece a la estructura"),
                ));
            }
            total = total
                .checked_add(required_u64(step, "count", path)?)
                .ok_or_else(|| event_error(path, "la distribución de pasos desborda"))?;
        }
    }
    if total != stuttered {
        return Err(event_error(
            path,
            format!(
                "la distribución de pasos suma {total} y las moléculas que tartamudearon son \
                 {stuttered}"
            ),
        ));
    }
    Ok(())
}

fn required_object<'a>(
    value: &'a Value,
    context: &str,
    path: &Path,
) -> Result<&'a serde_json::Map<String, Value>, Diagnostic> {
    value
        .as_object()
        .ok_or_else(|| event_error(path, format!("{context} no es objeto")))
}
fn required_array<'a>(
    value: &'a Value,
    context: &str,
    path: &Path,
) -> Result<&'a [Value], Diagnostic> {
    value
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| event_error(path, format!("{context} no es lista")))
}
fn required_value<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a Value, Diagnostic> {
    object
        .get(field)
        .ok_or_else(|| event_error(path, format!("falta {field}")))
}
fn required_non_null<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a Value, Diagnostic> {
    let value = required_value(object, field, path)?;
    if value.is_object() {
        Ok(value)
    } else {
        Err(event_error(path, format!("{field} no puede ser nulo")))
    }
}
fn required_string<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    required_value(object, field, path)?
        .as_str()
        .filter(|value| !value.is_empty())
        .ok_or_else(|| event_error(path, format!("{field} debe ser texto no vacío")))
}
fn require_text(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
    expected: &str,
) -> Result<(), Diagnostic> {
    if required_string(object, field, path)? == expected {
        Ok(())
    } else {
        Err(event_error(path, format!("{field} incompatible")))
    }
}
fn validate_closed(
    object: &serde_json::Map<String, Value>,
    fields: &[&str],
    context: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if object.len() == fields.len() && object.keys().all(|field| fields.contains(&field.as_str())) {
        Ok(())
    } else {
        Err(event_error(path, format!("campos inválidos en {context}")))
    }
}
fn event_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path.display().to_string())
    // Nombre **provisional**: `validate_event_manifest_bytes` lo retiqueta al
    // salir con el que corresponde a la versión que el documento declara. Aquí
    // no se conoce, y las ciento treinta y seis llamadas de este módulo no
    // tienen por qué conocerlo.
    .with("schema", "manifest-event-v2.schema.json")
    .with("detail", detail.into())
}
