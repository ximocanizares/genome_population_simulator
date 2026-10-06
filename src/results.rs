//! Documentos JSON mínimos producidos por F1.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Component, Path, PathBuf};

use serde_json::{Value, json};

use crate::adapters::{admit_toml_header, parse_declared_config};
use crate::application::{
    DeclaredDocumentVersion, EffectiveConfig, EffectiveDirectEvent, EffectiveInitialGenotypeModel,
    EffectiveReducedReference, EffectiveReferenceSource, EffectiveRepeatAllele,
    EffectiveRepeatSlippage, EffectiveTechnicalEvent, EffectiveTechnicalEventTarget,
    build_effective_config_with_entropy,
};
use crate::biological_contamination::contamination_target_pair_count;
use crate::configuration::ExpectedConfigKind;
use crate::diagnostics::{Diagnostic, diagnostic_matches_catalog_version};
use crate::event_manifest::{event_manifest_relative_path, validate_event_manifest_bytes};
use crate::execution::declared_small_event;
use crate::identity::{
    EventId, FounderKey, IndividualId, LibraryId, PopulationId, ReferenceId, SequenceId,
    derive_execution_id,
};
use crate::initial_individual_integration::IncompatibleGenotypeAuthorization;
use crate::initial_individual_materialization::{
    InitialIndividualsMaterializationError, materialize_initial_individuals,
};
use crate::integrity::{ArtifactIdentity, configuration_digest, identify_file};
use crate::population_snvs::{
    PreparedPopulationSnvs, validate_reusable_population_snvs,
    validate_reusable_population_snvs_bytes,
};
use crate::randomness::{
    EntropySource, F4TechnicalDecision, RootSeed, RuntimeSubseedScope, derive_f4_technical_subseed,
    derive_subseed, f4_uniform_closed,
};
use crate::reference::{
    AlphabetValidation, DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence,
    ReversibleText, SequenceClass,
};
use crate::reference_fasta::{
    CANONICAL_FASTA_VERSION, CANONICAL_LINE_BASE_COUNT, SequenceExpectation,
    reference_from_summary, validate_supplied_fasta, validate_supplied_fasta_summary,
};
use crate::reference_generation::COMPOSITION_TOLERANCE;
use crate::reference_generation::{BaseFrequency, GeneratedSequenceSpec};
use crate::reference_reduction::{
    EligibleSequence, MaterializedPlan, materialize_selection, select_reduced_reference,
};
use crate::reference_reduction_io::{
    REDUCTION_MAP_DOCUMENT_KIND, REDUCTION_MAP_SCHEMA_VERSION, ReductionMapParameters,
    build_reduced_reference, reduction_map_document, source_for, validate_reduction_map_document,
};
use crate::repeat_genotypes::RepeatAlleleDraw;
use crate::tandem_repeat::{BlockAnnotation, RepeatMotif, TandemRepeatTract};
use crate::technical_identity::derive_technical_identity_plan;
use crate::vcf_contract::{TBI_MAXIMUM_POSITION_EXCLUSIVE, VcfIndexKind, VcfOutputMode};
use crate::vcf_projection::VCF_NORMALIZATION_VERSION;
use crate::vcf_publication::{
    VCF_ARTIFACT_ID, VCF_ARTIFACT_KIND, VCF_BGZF_ARTIFACT_PATH, VCF_CSI_ARTIFACT_PATH,
    VCF_INDEX_ARTIFACT_ID, VCF_TBI_ARTIFACT_PATH, VCF_TEXT_ARTIFACT_PATH,
};

/// Versión vigente del manifiesto general. DTC-082 la activa junto con el
/// primer productor PE y no mantiene lectores F3.
pub const GENERAL_MANIFEST_SCHEMA_VERSION: &str = "4";

/// Versión del manifiesto general que abre este trabajo, para las ejecuciones
/// que publican repeticiones en tándem.
///
/// Es la **seis** y no la cinco porque la cinco ya la ocupa el manifiesto de
/// la fase de población fundadora, que es otro documento con la misma clase.
/// Una ejecución sin loci declarados sigue publicando la versión 4, byte a
/// byte como antes.
pub const GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS: &str = "6";

/// Versión vigente de la configuración efectiva para el productor PE F4.
pub const EFFECTIVE_CONFIG_SCHEMA_VERSION: &str = "4";

/// Versión de configuración efectiva que abre este trabajo, para los
/// documentos generales v4 que pueden declarar repeticiones en tándem. Un
/// documento general v2 sigue publicando la versión 4, byte a byte como antes.
pub const EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS: &str = "6";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum UnitState {
    Incomplete,
    Complete,
    Failed,
}

impl UnitState {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Incomplete => "incomplete",
            Self::Complete => "complete",
            Self::Failed => "failed",
        }
    }
}

/// Única descripción física de un fichero publicado por una ejecución.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArtifactRecord {
    pub id: String,
    pub kind: String,
    pub path: String,
    pub size_bytes: u64,
    pub digest: String,
}

impl ArtifactRecord {
    pub fn from_identity(
        id: impl Into<String>,
        kind: impl Into<String>,
        identity: &ArtifactIdentity,
    ) -> Result<Self, Diagnostic> {
        let path = identity
            .relative_path()
            .to_string_lossy()
            .replace('\\', "/");
        validate_relative_link(&path)?;
        Ok(Self {
            id: id.into(),
            kind: kind.into(),
            path,
            size_bytes: identity.size_bytes(),
            digest: identity.digest().to_owned(),
        })
    }

    #[must_use]
    pub fn to_json(&self) -> Value {
        json!({
            "id": self.id,
            "kind": self.kind,
            "path": self.path,
            "size_bytes": self.size_bytes,
            "checksum": {"algorithm": "sha256", "digest": self.digest},
        })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConfigurationInputRecord {
    pub kind: String,
    pub original_path: String,
    pub copy: ArtifactRecord,
}

#[derive(Clone, Copy)]
pub struct ReferenceOutput<'a> {
    pub reference: &'a Reference,
    pub fasta: &'a ArtifactRecord,
}

#[derive(Clone, Copy)]
pub struct ReducedReferenceOutput<'a> {
    pub effective: &'a EffectiveReducedReference,
    pub reference: &'a Reference,
    pub fasta: &'a ArtifactRecord,
    pub map: &'a ArtifactRecord,
    pub plan: &'a MaterializedPlan,
    pub source_fasta: &'a ArtifactRecord,
    pub state: UnitState,
}

/// Enlace semántico de una referencia al catálogo físico único.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ManifestReferenceLink {
    pub id: String,
    pub fasta_artifact_id: String,
}

/// Referencia F1 resuelta desde un manifiesto general ya validado.
///
/// Esta frontera conserva juntos el ID solicitado, el enlace semántico del
/// manifiesto y la identidad física comprobada del FASTA. Los consumidores F2
/// no deben aceptar un ID ni un FASTA alternativos fuera de esta resolución.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ResolvedManifestReference {
    reference_id: ReferenceId,
    fasta_path: PathBuf,
    fasta: ArtifactRecord,
}

impl ResolvedManifestReference {
    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub fn fasta_path(&self) -> &Path {
        &self.fasta_path
    }

    #[must_use]
    pub const fn fasta(&self) -> &ArtifactRecord {
        &self.fasta
    }
}

/// Estado alcanzado por una reducción solicitada, incluso si la ejecución falla.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ManifestReductionLink {
    pub id: String,
    pub source_reference_id: String,
    pub eligible_bases: u64,
    pub region_size_bases: u64,
    pub target_mode: String,
    pub target_declared: String,
    pub target_bases: u64,
    pub root_seed: String,
    pub selection_subseed: String,
    pub selection_order: Vec<u64>,
    pub state: UnitState,
    pub fasta_artifact_id: Option<String>,
    pub map_artifact_id: Option<String>,
    pub source_fasta_artifact_id: Option<String>,
    pub map_document_kind: String,
    pub map_schema_version: String,
}

/// Enlace verificable al metadato de la matriz SNV canónica. En `reuse` el
/// recurso permanece externo e inmutable, por eso no hay identificador de
/// artefacto publicado.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PopulationSnvsManifestLink {
    pub mode: String,
    pub metadata_sha256: String,
    pub metadata_artifact_id: Option<String>,
    pub normalization_audit: Value,
    pub initial_population_statistics: Option<Value>,
}

/// Enlace de un VCF opcional a su única referencia, muestras e inventario
/// físico. Las identidades SHA-256 y tamaños viven en el catálogo compartido.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfManifestLink {
    pub state: UnitState,
    pub reference_id: String,
    pub samples: Vec<String>,
    pub mode: VcfOutputMode,
    pub vcf_artifact_id: String,
    pub index_kind: Option<VcfIndexKind>,
    pub index_artifact_id: Option<String>,
}

/// Enlace del manifiesto general al único manifiesto de un evento real. El
/// detalle técnico continúa en el documento de evento; la raíz solo enlaza su
/// identidad, estado y artefacto físico validado.
#[derive(Clone, Debug, PartialEq)]
pub struct ManifestEventLink {
    pub event_id: String,
    pub state: UnitState,
    pub manifest_artifact_id: String,
    /// Resumen agregado duplicado en la raíz para inspección sin abrir el
    /// documento de evento. No transporta procedencia por pareja.
    pub chimera_rate: f64,
    pub realized_chimera_pairs: u64,
    /// Perfil incorporado que usó el evento, con su huella. El inventario del
    /// manifiesto general repite lo que el manifiesto del evento ya publicó,
    /// en vez de deducirlo otra vez.
    pub profile_id: &'static str,
    pub profile_sha256: &'static str,
    /// Catálogo del que salieron **las decisiones de este evento**. El
    /// inventario del manifiesto general repite el mismo nombre que el
    /// manifiesto del evento, en vez de deducirlo otra vez.
    pub rng_catalog: &'static str,
}

pub struct GeneralManifest<'a> {
    pub execution_id: &'a str,
    pub state: UnitState,
    pub configuration_input: Option<&'a ConfigurationInputRecord>,
    pub effective_config: Option<&'a ArtifactRecord>,
    pub references: &'a [ManifestReferenceLink],
    pub reduced_references: &'a [ManifestReductionLink],
    pub population_snvs: Option<&'a PopulationSnvsManifestLink>,
    pub vcf_outputs: &'a [VcfManifestLink],
    pub events: Option<&'a [ManifestEventLink]>,
    pub authorizations: crate::execution::ExecutionAuthorizations,
    pub artifacts: &'a [ArtifactRecord],
    pub diagnostics: &'a [Diagnostic],
    /// Sección de repeticiones en tándem, tal y como la construye
    /// [`repeat_loci_manifest_json`].
    ///
    /// Su presencia —y solo su presencia— promueve el documento a la versión
    /// 6. Una ejecución sin loci declarados la deja en `None` y publica
    /// exactamente el documento de antes.
    pub repeat_loci: Option<&'a Value>,
}

pub fn effective_config_document_with_reference(
    config: &EffectiveConfig,
    output: ReferenceOutput<'_>,
) -> Result<Value, Diagnostic> {
    effective_config_document_with_reductions(config, output, &[])
}

pub fn effective_config_document_with_reductions(
    config: &EffectiveConfig,
    output: ReferenceOutput<'_>,
    reductions: &[ReducedReferenceOutput<'_>],
) -> Result<Value, Diagnostic> {
    validate_reference_output(config, output)?;
    let composition = config.composition().map(|composition| {
        let values = composition.values();
        json!({"A": values[0], "C": values[1], "G": values[2], "T": values[3], "N": composition.unknown_base_fraction()})
    });
    let mut document = json!({
        "effective_config_kind": "gps_effective_config",
        "schema_version": config.document_version().effective_config_schema_version(),
        "source_path": config.source_path(),
        "ploidy": config.ploidy().copy_count(),
        "normalization": {
            "tolerance": config.normalization().tolerance(),
            "canonical_order": "population_id,contig,position,ref,alt,variant_id,source_id",
            "sum_algorithm": "neumaier_sum_v1",
        },
        "seed": {
            "root": config.reproducibility().root_seed().to_hex(),
            "origin": config.reproducibility().seed_origin().as_str(),
        },
        "reference": reference_json(config, output, composition),
        "reduced_references": reductions.iter().map(reduction_json).collect::<Vec<_>>(),
        "initial_populations": effective_initial_populations_json(config),
        "direct_events": config.direct_events().iter().map(direct_event_json).collect::<Vec<_>>(),
        "vcf_output": effective_vcf_output_json(config),
    });
    if let Some(event) = config.technical_event() {
        document["technical_event"] = technical_event_json(config, event);
    }
    if config.document_version() == DeclaredDocumentVersion::V4 {
        // Los mínimos de la guarda se publican siempre que el documento pueda
        // declarar repeticiones, se hayan declarado o no: quien lea el
        // documento tiene que ver contra qué umbral se admitió el tracto.
        document["repeat_loci"] = repeat_loci_json(config)?;
        document["repeat_admission"] = json!({
            "min_purity": config.repeat_admission().min_purity(),
            "min_longest_run_units": config.repeat_admission().min_longest_run_units(),
        });
        document["repeat_output_budget"] = repeat_output_budget_json(config);
    }
    Ok(document)
}

/// Publica cada locus con su motivo, su intervalo, su tope, sus alelos y la
/// **auditoría** de la normalización de sus frecuencias.
///
/// De cada alelo se publican la secuencia, la frecuencia **declarada tal y como
/// se escribió**, la frecuencia **efectiva** ya normalizada, y los valores
/// realizados que se derivan de la secuencia: unidades, pureza, interrupciones y
/// tirada más larga. La pureza no se puede declarar, solo leer.
///
/// Se publican las dos frecuencias porque, si la suma declarada no era
/// exactamente uno, la normalizada **no es la que el usuario escribió**, y
/// publicar solo una de las dos bajo un nombre que sugiere la otra es
/// exactamente el error que DTC-049 obliga a evitar.
fn repeat_loci_json(config: &EffectiveConfig) -> Result<Value, Diagnostic> {
    let mut loci = Vec::with_capacity(config.repeat_loci().len());
    for locus in config.repeat_loci() {
        let normalization = locus.normalization();
        loci.push(json!({
            "id": locus.id(),
            "sequence_id": locus.sequence_id(),
            "start": locus.start()?,
            "end": locus.end()?,
            "motifs": locus.motifs().iter().map(RepeatMotif::as_text).collect::<Vec<_>>(),
            "motif_lengths": locus.motifs().iter().map(RepeatMotif::length).collect::<Vec<_>>(),
            "max_units": locus.max_units(),
            "admission": {
                "min_purity": locus.admission().min_purity(),
                "min_longest_run_units": locus.admission().min_longest_run_units(),
            },
            "frequency_normalization": {
                "tolerance": normalization.tolerance(),
                "sum_algorithm": "neumaier_sum_v1",
                "canonical_order": "locus_id,allele_index",
                "original_sum": normalization.original_sum(),
                "original_delta": normalization.original_delta(),
                "effective_sum": normalization.effective_sum(),
                "effective_delta": normalization.effective_delta(),
            },
            "alleles": locus.alleles().iter()
                .map(|allele| Value::Object(repeat_allele_json(allele)))
                .collect::<Vec<_>>(),
        }));
        // **Con inserción, `start`/`end` son del genoma final** (DTC-098), y
        // lo que el usuario escribió en el de partida se publica aparte.
        if let Some(declared) = locus.declared_position() {
            loci.last_mut()
                .and_then(Value::as_object_mut)
                .ok_or_else(|| {
                    Diagnostic::internal_invariant("results", "repeat_locus_declared")
                        .with("detail", "el locus recién insertado no es un objeto")
                })?
                .insert(
                    "declared".to_owned(),
                    match declared {
                        crate::application::DeclaredLocusPosition::Interval { start, end } => {
                            json!({"start": start, "end": end})
                        }
                        crate::application::DeclaredLocusPosition::InsertAfter(point) => {
                            json!({"insert_after": point})
                        }
                    },
                );
        }
        if let Some(units) = locus.reference_units() {
            loci.last_mut()
                .and_then(Value::as_object_mut)
                .ok_or_else(|| {
                    Diagnostic::internal_invariant("results", "repeat_locus_reference_units")
                        .with("detail", "el locus recién insertado no es un objeto")
                })?
                .insert("reference_units".to_owned(), json!(units));
        }
        if let Some(slippage) = locus.slippage() {
            loci.last_mut()
                .ok_or_else(|| {
                    Diagnostic::internal_invariant("results", "repeat_locus_slippage")
                        .with("detail", "el locus recién insertado desapareció")
                })?
                .as_object_mut()
                .ok_or_else(|| {
                    Diagnostic::internal_invariant("results", "repeat_locus_slippage")
                        .with("detail", "el locus publicado no es un objeto")
                })?
                .insert("slippage".to_owned(), slippage_json(slippage));
        }
    }
    Ok(Value::Array(loci))
}

/// Juez independiente de la contabilidad de bases de tracto.
///
/// # Qué comprueba, y por qué no es un espejo
///
/// Vuelve a contar las bases de tracto **releyendo los bytes publicados**: el
/// `truth.jsonl` como JSON línea a línea y el VCF como texto, sin tocar ni una
/// estructura del productor. Con esas dos cifras exige tres cosas:
///
/// 1. que coincidan con las que el manifiesto publica —si el productor contó
///    mal o mintió, aquí se ve—;
/// 2. que la cota del manifiesto sea la de la configuración efectiva;
/// 3. que lo realizado **no supere** la cota.
///
/// Las tres comparaciones enfrentan caminos distintos: la cota sale de lo
/// declarado, la cifra del manifiesto del recorrido de las estructuras en
/// memoria, y la de aquí de los bytes de los ficheros. Si los tres fueran la
/// misma función, ninguna de las tres podría fallar.
pub fn validate_published_tract_output_accounting(
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
    truth_path: &Path,
    vcf_path: Option<&Path>,
) -> Result<(), Diagnostic> {
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let Some(section) = manifest.get("repeat_loci") else {
        // Sin loci no hay contabilidad de tracto que juzgar, y exigirla sería
        // inventar una obligación donde el manifiesto no promete nada.
        return Ok(());
    };
    let accounting = section
        .get("tract_output")
        .and_then(Value::as_object)
        .ok_or_else(|| repeat_section_error(manifest_path, "falta la contabilidad de tracto"))?;
    let read = |key: &str| -> Result<u64, Diagnostic> {
        accounting.get(key).and_then(Value::as_u64).ok_or_else(|| {
            repeat_section_error(manifest_path, "una cifra de tracto no es un entero")
        })
    };

    let effective: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let declared_bound = effective
        .get("repeat_output_budget")
        .and_then(|budget| budget.get("total_tract_bases_bound"))
        .and_then(Value::as_u64)
        .ok_or_else(|| {
            effective_error(
                effective_path,
                "la configuración efectiva no publica la cota de bases de tracto",
            )
        })?;
    if read("total_tract_bases_bound")? != declared_bound {
        return Err(repeat_section_error(
            manifest_path,
            "la cota del manifiesto no es la de la configuración efectiva",
        ));
    }

    let recomputed_truth = recount_truth_tract_bases(truth_path)?;
    // Los loci tal como el manifiesto los publica, en base uno e inclusivos:
    // a cada registro de repetición se le atribuye el primero, por posición,
    // cuyo intervalo contiene.
    let loci = section
        .get("loci")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|locus| {
            Some((
                locus.get("sequence_id")?.as_str()?.to_owned(),
                locus.get("start")?.as_u64()?,
                locus.get("end")?.as_u64()?,
                locus.get("id")?.as_str()?.to_owned(),
            ))
        })
        .collect::<BTreeSet<_>>();
    let (recomputed_vcf, recomputed_by_locus) = match vcf_path {
        None => (0, BTreeMap::new()),
        Some(path) => recount_vcf_tract_bases(path, &loci, manifest_path)?,
    };
    if recomputed_truth != read("realized_truth_tract_bases")?
        || recomputed_vcf != read("realized_vcf_tract_bases")?
    {
        return Err(repeat_section_error(
            manifest_path,
            "las bases de tracto realizadas no son las que los documentos publicados llevan",
        ));
    }
    // **La parte VCF de la cota es lo que el VCF lleva** (DTC-097): se mide al
    // publicar, de modo que tiene que coincidir con lo que este juez recuenta
    // sobre los bytes, y la configuración efectiva tiene que decir lo mismo,
    // en total y sumando por locus.
    if read("vcf_tract_bases_bound")? != recomputed_vcf {
        return Err(repeat_section_error(
            manifest_path,
            "la parte VCF de la cota no es lo que el VCF publicado lleva",
        ));
    }
    let effective_budget = effective.get("repeat_output_budget");
    if effective_budget
        .and_then(|budget| budget.get("vcf_tract_bases_bound"))
        .and_then(Value::as_u64)
        != Some(recomputed_vcf)
    {
        return Err(effective_error(
            effective_path,
            "la configuración efectiva no publica las bases de tracto que el VCF lleva",
        ));
    }
    let per_locus = effective_budget
        .and_then(|budget| budget.get("per_locus"))
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .map(|locus| locus.get("vcf_tract_bases_bound").and_then(Value::as_u64))
        .try_fold(0_u64, |total, bases| {
            bases.map(|bases| total.saturating_add(bases))
        });
    if per_locus != Some(recomputed_vcf) {
        return Err(effective_error(
            effective_path,
            "las bases de tracto del VCF por locus no suman las del VCF publicado",
        ));
    }
    // Y cada locus, lo que llevan **sus** registros: sin esto, atribuir todo
    // al primer locus de la secuencia conservaría la suma.
    for locus in effective_budget
        .and_then(|budget| budget.get("per_locus"))
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
    {
        let id = locus.get("id").and_then(Value::as_str).unwrap_or_default();
        let published = locus.get("vcf_tract_bases_bound").and_then(Value::as_u64);
        if published != Some(recomputed_by_locus.get(id).copied().unwrap_or(0)) {
            return Err(effective_error(
                effective_path,
                format!(
                    "las bases de tracto del VCF del locus {id} no son las que llevan sus registros"
                ),
            ));
        }
    }
    let realized = recomputed_truth.saturating_add(recomputed_vcf);
    if realized != read("realized_total_tract_bases")? {
        return Err(repeat_section_error(
            manifest_path,
            "el total realizado no es la suma de sus dos términos",
        ));
    }
    if realized > declared_bound {
        return Err(repeat_section_error(
            manifest_path,
            "las bases de tracto realizadas superan la cota declarada",
        ));
    }
    Ok(())
}

/// Recuenta las bases de tracto releyendo el `truth.jsonl` como JSON.
fn recount_truth_tract_bases(path: &Path) -> Result<u64, Diagnostic> {
    let text = std::fs::read_to_string(path)
        .map_err(|error| manifest_error(path, format!("verdad no legible: {error}")))?;
    let mut total = 0_u64;
    for line in text.lines().filter(|line| !line.trim().is_empty()) {
        let record: Value = serde_json::from_str(line)
            .map_err(|error| manifest_error(path, format!("línea de verdad ilegible: {error}")))?;
        let mut add = |value: Option<&Value>| {
            if let Some(bases) = value.and_then(Value::as_str) {
                total = total.saturating_add(bases.len() as u64);
            }
        };
        for piece in record
            .get("current_pieces")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .chain(
                record
                    .get("history")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten(),
            )
        {
            if piece.get("kind").and_then(Value::as_str) == Some("repeat_allele") {
                add(piece.get("bases"));
            }
        }
        for allele in record
            .get("alleles")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
        {
            if allele.get("kind").and_then(Value::as_str) == Some("repeat") {
                add(allele.get("bases"));
            }
        }
    }
    Ok(total)
}

/// Recuenta las bases de tracto releyendo el VCF como texto.
///
/// Un registro es de repetición si publica el recuento de secuencias
/// repetidas; sus bases son las de `REF` más las de cada `ALT`, que desde la
/// decisión S-09 llevan la secuencia completa del alelo.
fn recount_vcf_tract_bases(
    path: &Path,
    loci: &BTreeSet<(String, u64, u64, String)>,
    manifest_path: &Path,
) -> Result<(u64, BTreeMap<String, u64>), Diagnostic> {
    let bytes = std::fs::read(path)
        .map_err(|error| manifest_error(path, format!("VCF no legible: {error}")))?;
    let text = if path.extension().is_some_and(|extension| extension == "gz") {
        let mut bytes = std::io::Cursor::new(bytes);
        crate::bgzf_blocks::reject_oversized_bgzf_blocks(&mut bytes)
            .map_err(|fault| manifest_error(path, format!("VCF BGZF no legible: {fault}")))?;
        let mut reader = noodles_bgzf::io::Reader::new(bytes);
        let mut decoded = String::new();
        std::io::Read::read_to_string(&mut reader, &mut decoded)
            .map_err(|error| manifest_error(path, format!("VCF BGZF no legible: {error}")))?;
        decoded
    } else {
        String::from_utf8(bytes).map_err(|_| manifest_error(path, "VCF no es UTF-8"))?
    };
    let mut total = 0_u64;
    let mut by_locus = BTreeMap::<String, u64>::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let fields = line.split('\t').collect::<Vec<_>>();
        let Some(info) = fields.get(7) else { continue };
        if !info.split(';').any(|field| field.starts_with("RN=")) {
            continue;
        }
        let mut bases = fields[3].len() as u64;
        for alternate in fields[4].split(',') {
            bases = bases.saturating_add(alternate.len() as u64);
        }
        total = total.saturating_add(bases);
        let position = fields[1].parse::<u64>().unwrap_or(0);
        let last = position
            .saturating_add(fields[3].len() as u64)
            .saturating_sub(1);
        // Un registro de repetición contiene siempre el intervalo de su
        // locus; si no contiene ninguno, el fichero no es el que se publicó.
        let Some((_, _, _, id)) = loci.iter().find(|(sequence, start, end, _)| {
            sequence == fields[0] && position <= *start && *end <= last
        }) else {
            return Err(repeat_section_error(
                manifest_path,
                "un registro de repetición del VCF no contiene el intervalo de ningún locus",
            ));
        };
        let entry = by_locus.entry(id.clone()).or_insert(0);
        *entry = entry.saturating_add(bases);
    }
    Ok((total, by_locus))
}

/// Bases de tracto que los documentos publicados llevan de verdad.
///
/// Las dos cifras se cuentan **sobre lo publicado**, cada una donde se produce:
/// la de la verdad recorriendo el documento ya construido, la del VCF
/// recorriendo los registros ya agrupados. Ninguna de las dos sabe nada de la
/// cota.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RealizedTractOutput {
    pub truth_bases: u64,
    pub vcf_bases: u64,
    /// Las bases del VCF por locus. Desde el 2 de octubre de 2026 son también
    /// la parte VCF de la cota publicada: se miden, no se estiman (DTC-097).
    pub vcf_by_locus: Vec<(String, u64)>,
}

impl RealizedTractOutput {
    #[must_use]
    pub const fn total_bases(&self) -> u64 {
        self.truth_bases.saturating_add(self.vcf_bases)
    }
}

/// Escribe en la configuración efectiva ya construida la parte VCF de la cota
/// de bases de tracto, **medida** al publicar el VCF (DTC-097).
///
/// La configuración la deja a cero, porque antes de construir las copias no
/// se conoce; esta función la completa con lo que el VCF lleva, por locus y en
/// total, y rehace el total de la cota. Sin sección de repeticiones no hace
/// nada.
pub fn publish_measured_vcf_tract_bases(document: &mut Value, realized: &RealizedTractOutput) {
    let Some(budget) = document
        .get_mut("repeat_output_budget")
        .and_then(Value::as_object_mut)
    else {
        return;
    };
    let truth = budget
        .get("truth_tract_bases_bound")
        .and_then(Value::as_u64)
        .unwrap_or(0);
    budget.insert(
        "vcf_tract_bases_bound".to_owned(),
        json!(realized.vcf_bases),
    );
    budget.insert(
        "total_tract_bases_bound".to_owned(),
        json!(truth.saturating_add(realized.vcf_bases)),
    );
    if let Some(loci) = budget.get_mut("per_locus").and_then(Value::as_array_mut) {
        for locus in loci {
            let measured = locus
                .get("id")
                .and_then(Value::as_str)
                .and_then(|id| {
                    realized
                        .vcf_by_locus
                        .iter()
                        .find(|(locus, _)| locus == id)
                        .map(|(_, bases)| *bases)
                })
                .unwrap_or(0);
            locus["vcf_tract_bases_bound"] = json!(measured);
        }
    }
}

/// Cota superior de bases de tracto publicadas, tal y como la calculó la
/// frontera de configuración.
///
/// Se publica **siempre** que el documento pueda declarar repeticiones, aunque
/// no declare ninguna: quien lea el documento tiene que poder ver contra qué
/// techo se admitió el escenario, y un techo ausente no se distingue de un
/// techo infinito.
fn repeat_output_budget_json(config: &EffectiveConfig) -> Value {
    let budget = config.repeat_output_budget();
    json!({
        "truth_tract_bases_bound": budget.truth_bases(),
        // Se mide al publicar el VCF y lo completa
        // `publish_measured_vcf_tract_bases`; antes no se conoce.
        "vcf_tract_bases_bound": 0,
        "total_tract_bases_bound": budget.truth_bases(),
        "declared_max_tract_output_bases": budget.declared_cap(),
        "truth_copies_per_haplotype_copy":
            crate::application::EffectiveRepeatOutputBudget::TRUTH_COPIES_PER_HAPLOTYPE_COPY,
        "per_locus": budget.per_locus().iter().map(|locus| json!({
            "id": locus.locus_id(),
            "truth_tract_bases_bound": locus.truth_bases(),
            "vcf_tract_bases_bound": 0,
        })).collect::<Vec<_>>(),
    })
}

/// Bases de secuencia de tracto que un documento de verdad **lleva de verdad**.
///
/// # Por qué se cuenta aquí y no se deriva de la cota
///
/// Esta función recorre el documento ya construido y suma las bases de cada
/// tracto que publica: la vista vigente de cada copia, su historial y los
/// alelos de cada caso genotípico. No sabe nada de alelos declarados, de
/// individuos ni de ploidía, que es de lo único que se alimenta la cota. Son
/// dos caminos que no se tocan, y por eso compararlos dice algo.
#[must_use]
pub fn realized_truth_tract_bases(truth: &crate::truth_schema::TruthDocument) -> u64 {
    use crate::truth_schema::{TruthAllele, TruthCurrentPiece, TruthHistoryEvent, TruthRecord};

    let mut total = 0_u64;
    for record in &truth.records {
        match record {
            TruthRecord::HaplotypeCopy(copy) => {
                for piece in &copy.current_pieces {
                    if let TruthCurrentPiece::RepeatAllele { bases, .. } = piece {
                        total = total.saturating_add(bases.len() as u64);
                    }
                }
                for event in &copy.history {
                    if let TruthHistoryEvent::RepeatAllele { bases, .. } = event {
                        total = total.saturating_add(bases.len() as u64);
                    }
                }
            }
            TruthRecord::Genotype(genotype) => {
                for allele in &genotype.alleles {
                    if let TruthAllele::Repeat { bases, .. } = allele {
                        total = total.saturating_add(bases.len() as u64);
                    }
                }
            }
        }
    }
    total
}

/// Un alelo de repetición: su secuencia, sus frecuencias declarada y efectiva
/// y los **valores realizados** que se leen de su secuencia.
///
/// Los bloques se publican como **listas paralelas**, que es la forma que el
/// estándar VCF 4.5 define para una repetición compuesta: cuántas secuencias
/// repetidas tiene el alelo, y por cada una su motivo, la longitud del motivo,
/// el número de copias —decimal, porque puede acabar en unidad parcial— y las
/// bases totales. Publicarlas aquí con la misma forma hace que la proyección al
/// VCF sea un volcado y no una traducción.
///
/// La pureza, las interrupciones y la tirada más larga **no se pueden
/// declarar**: se calculan sobre la secuencia cada vez (decisión S-02). Por eso
/// el manifiesto y la configuración efectiva comparten esta función en vez de
/// escribirlas dos veces: dos copias podrían decir cosas distintas del mismo
/// alelo.
fn repeat_allele_json(allele: &EffectiveRepeatAllele) -> serde_json::Map<String, Value> {
    let annotation = allele.tract().annotation();
    let Value::Object(object) = json!({
        "sequence": allele.tract().as_text(),
        // **El origen va por alelo, no por locus.** Un documento con loci
        // declarados y loci generados tiene que poder distinguirlos alelo a
        // alelo sin volver a la configuración.
        "origin": allele.origin().as_str(),
        "declared_frequency": allele.declared_frequency(),
        "effective_frequency": allele.frequency(),
        "units": annotation.units_as_f64(),
        "purity": annotation.purity_as_f64(),
        "longest_run_units": annotation.longest_run_units(),
        "repeat_sequence_count": annotation.blocks().len(),
        "block_motifs": annotation.blocks().iter()
            .map(|block| block.motif().as_text()).collect::<Vec<_>>(),
        "block_motif_lengths": annotation.blocks().iter()
            .map(|block| block.motif().length()).collect::<Vec<_>>(),
        "block_unit_counts": annotation.blocks().iter()
            .map(BlockAnnotation::units_as_f64).collect::<Vec<_>>(),
        "block_bases": annotation.blocks().iter()
            .map(BlockAnnotation::length).collect::<Vec<_>>(),
        "interruptions": annotation.interruptions().iter().map(|interruption| {
            json!({"offset": interruption.start(), "length": interruption.length()})
        }).collect::<Vec<_>>(),
    }) else {
        unreachable!("json! de un objeto literal es un objeto");
    };
    object
}

/// El generador de alelos por deslizamiento de un locus, tal y como se publica.
///
/// # Qué lleva y por qué cada cosa
///
/// - **Los parámetros** del generador, que dicen qué se pidió.
/// - **La versión del algoritmo** y **la del catálogo de ámbitos de azar**, sin
///   las cuales dos ejecuciones con la misma semilla y resultados distintos
///   serían indistinguibles de un defecto.
/// - **La contabilidad de la generación**: linajes, pasos sorteados, pasos
///   aplicados y, sobre todo, **eventos no aplicados al agotar el tracto**. Sin
///   esa última cifra, un servicio que ignorase en silencio los eventos
///   sobrantes y otro que los contase producirían el mismo tracto final y nadie
///   notaría la diferencia.
/// - **Las tres limitaciones declaradas del modelo**, palabra por palabra.
///   Son decisiones biológicas, no detalles: quien encadene un análisis sobre
///   estos datos tiene que leerlas sin abrir el código.
fn slippage_json(slippage: &EffectiveRepeatSlippage) -> Value {
    let settings = slippage.settings();
    let mut document = json!({
        "model": crate::repeat_slippage::SLIPPAGE_MODEL_STEPWISE,
        "step_units": crate::repeat_slippage::SLIPPAGE_STEP_UNITS,
        "mean_events": settings.mean_events(),
        "expansion_bias": settings.expansion_bias(),
        "algorithm": crate::repeat_slippage::SLIPPAGE_ALGORITHM,
        "count_algorithm": crate::repeat_slippage::SLIPPAGE_COUNT_ALGORITHM,
        "bias_algorithm": crate::repeat_slippage::SLIPPAGE_BIAS_ALGORITHM,
        "rng_scope": crate::repeat_slippage::SLIPPAGE_SCOPE,
        "rng_scope_catalog_version": crate::repeat_slippage::SLIPPAGE_SCOPE_CATALOG_VERSION,
        "limitations": crate::repeat_slippage::SLIPPAGE_LIMITATIONS,
    });
    if let Some(generation) = slippage.generation() {
        document["lineages"] = json!(generation.lineages());
        document["drawn_events"] = json!(generation.drawn_events());
        document["applied_events"] = json!(generation.applied_events());
        document["unapplied_events"] = json!(generation.unapplied_events());
    }
    document
}

/// Sección de repeticiones en tándem del manifiesto general.
///
/// # Qué añade respecto de la configuración efectiva, y por qué no basta con
/// aquella
///
/// La configuración efectiva describe **lo que se pidió**: los loci, sus
/// parámetros, sus alelos y las frecuencias declarada y normalizada. Lo que
/// solo el manifiesto puede decir es **lo que salió**: cuántas copias recibió
/// cada alelo en el sorteo y qué proporción representa. Copiar ahí la
/// frecuencia declarada sería publicar una afirmación falsa con nombre de
/// medida, y en un corpus pequeño la diferencia entre las dos es grande.
///
/// # Por qué la procedencia de la referencia va aquí
///
/// De ella depende el significado de todo lo demás: con una referencia real el
/// tracto **ya estaba ahí** con las interrupciones que tuviera, y la pureza es
/// un valor observado; con una simulada, el tracto lo puso el propio simulador
/// (decisión S-06). Quien lea el manifiesto sin la configuración tiene que
/// poder distinguir los dos casos.
///
/// Devuelve `None` cuando la configuración no declara ningún locus: sin loci no
/// hay sección, el manifiesto conserva su versión anterior y sus bytes son los
/// de siempre.
pub fn repeat_loci_manifest_json(
    config: &EffectiveConfig,
    draws: &[RepeatAlleleDraw],
    realized: &RealizedTractOutput,
) -> Result<Option<Value>, Diagnostic> {
    if config.repeat_loci().is_empty() {
        return Ok(None);
    }
    let mut loci = Vec::with_capacity(config.repeat_loci().len());
    for locus in config.repeat_loci() {
        let drawn = draws
            .iter()
            .filter(|draw| draw.locus_id() == locus.id())
            .collect::<Vec<_>>();
        let drawn_copies = drawn.len() as u64;
        let alleles = locus
            .alleles()
            .iter()
            .enumerate()
            .map(|(index, allele)| {
                let realized_copies = drawn
                    .iter()
                    .filter(|draw| draw.allele_index() == index)
                    .count() as u64;
                let mut object = repeat_allele_json(allele);
                object.insert("realized_copies".to_owned(), json!(realized_copies));
                // Sin copias sorteadas no hay proporción que publicar, y un
                // cero sería tan falso como cualquier otra cifra inventada.
                object.insert(
                    "realized_frequency".to_owned(),
                    if drawn_copies == 0 {
                        Value::Null
                    } else {
                        json!(realized_copies as f64 / drawn_copies as f64)
                    },
                );
                Value::Object(object)
            })
            .collect::<Vec<_>>();
        loci.push(json!({
            "id": locus.id(),
            "sequence_id": locus.sequence_id(),
            "start": locus.start()?,
            "end": locus.end()?,
            "motifs": locus.motifs().iter().map(RepeatMotif::as_text).collect::<Vec<_>>(),
            "motif_lengths": locus.motifs().iter().map(RepeatMotif::length).collect::<Vec<_>>(),
            "max_units": locus.max_units(),
            "admission": {
                "min_purity": locus.admission().min_purity(),
                "min_longest_run_units": locus.admission().min_longest_run_units(),
            },
            "drawn_copies": drawn_copies,
            "alleles": alleles,
        }));
        if let Some(slippage) = locus.slippage() {
            loci.last_mut()
                .ok_or_else(|| {
                    Diagnostic::internal_invariant("results", "manifest_locus_slippage")
                        .with("detail", "el locus recién insertado desapareció")
                })?
                .as_object_mut()
                .ok_or_else(|| {
                    Diagnostic::internal_invariant("results", "manifest_locus_slippage")
                        .with("detail", "el locus publicado no es un objeto")
                })?
                .insert("slippage".to_owned(), slippage_json(slippage));
        }
    }
    Ok(Some(json!({
        "reference_id": config.simulation_reference_id(),
        "reference_source": match config.reference_source() {
            crate::application::EffectiveReferenceSource::Simulated => "simulated",
            crate::application::EffectiveReferenceSource::Supplied { .. } => "supplied",
        },
        "admission_defaults": {
            "min_purity": config.repeat_admission().min_purity(),
            "min_longest_run_units": config.repeat_admission().min_longest_run_units(),
        },
        // La cota se repite aquí **desde la configuración efectiva**, no se
        // recalcula: quien lea solo el manifiesto tiene que poder comparar las
        // dos cifras sin abrir el otro documento, y el juez independiente
        // comprueba que las dos coinciden.
        "tract_output": {
            "truth_tract_bases_bound": config.repeat_output_budget().truth_bases(),
            // La parte del VCF no se estima: es lo que el VCF lleva, medido al
            // publicarlo (DTC-097). Sin VCF, o si la ejecución falla antes de
            // escribirlo, es cero.
            "vcf_tract_bases_bound": realized.vcf_bases,
            "total_tract_bases_bound": config
                .repeat_output_budget()
                .truth_bases()
                .saturating_add(realized.vcf_bases),
            "realized_truth_tract_bases": realized.truth_bases,
            "realized_vcf_tract_bases": realized.vcf_bases,
            "realized_total_tract_bases": realized.total_bases(),
        },
        "loci": loci,
    })))
}

fn effective_initial_populations_json(config: &EffectiveConfig) -> Value {
    Value::Array(
        config
            .initial_populations()
            .iter()
            .map(|population| {
                let genotype_model = match population.genotype_model() {
                    EffectiveInitialGenotypeModel::HardyWeinberg => json!({"kind": "hwe"}),
                    EffectiveInitialGenotypeModel::FDiploid {
                        inbreeding_coefficient,
                    } => json!({"kind": "f_diploid", "inbreeding_coefficient": inbreeding_coefficient}),
                };
                let mut result = json!({
                    "population_id": population.population_id().as_str(),
                    "individual_count": population.individual_count(),
                    "individual_ids": population.individual_ids().iter().map(|id| id.as_str()).collect::<Vec<_>>(),
                    "ploidy": population.ploidy().copy_count(),
                    "genotype_model": genotype_model,
                });
                if let Some(founder_keys) = population.founder_keys() {
                    result["founder_keys"] =
                        json!(founder_keys.iter().map(|key| key.as_str()).collect::<Vec<_>>());
                }
                result
            })
            .collect(),
    )
}

fn direct_event_json(event: &EffectiveDirectEvent) -> Value {
    let common = |kind: &str| {
        json!({
            "event_id": event.event_id().as_str(),
            "individual_id": event.individual_id().as_str(),
            "phase": event.phase().copy_index(),
            "kind": kind,
            "contig": event.sequence_id().as_str(),
        })
    };
    match event {
        EffectiveDirectEvent::Snp {
            position_1based,
            reference,
            alternate,
            ..
        } => {
            let mut value = common("snp");
            let object = value.as_object_mut().expect("JSON object");
            object.insert("position_1based".to_owned(), json!(position_1based));
            object.insert("ref".to_owned(), json!(char::from(*reference)));
            object.insert("alt".to_owned(), json!(char::from(*alternate)));
            value
        }
        EffectiveDirectEvent::Insertion { point, bases, .. } => {
            let mut value = common("insertion");
            let object = value.as_object_mut().expect("JSON object");
            object.insert("point".to_owned(), json!(point));
            object.insert(
                "bases".to_owned(),
                json!(std::str::from_utf8(bases).expect("effective bases are canonical ASCII")),
            );
            value
        }
        EffectiveDirectEvent::Deletion { start, end, .. } => {
            let mut value = common("deletion");
            let object = value.as_object_mut().expect("JSON object");
            object.insert("start".to_owned(), json!(start));
            object.insert("end".to_owned(), json!(end));
            value
        }
    }
}

fn technical_event_json(config: &EffectiveConfig, event: &EffectiveTechnicalEvent) -> Value {
    let target = match event.target() {
        EffectiveTechnicalEventTarget::CoverageWgs(value) => json!({"coverage_wgs": value}),
        EffectiveTechnicalEventTarget::ReadCount(value) => json!({"read_count": value}),
    };
    let stutter_active = config.stutter_plan().is_active();
    // **Sin `expect` y sin recalcular la huella.** La escritura del documento
    // no puede entrar en pánico porque el recurso esté alterado —eso lo caza la
    // admisión, que sí lo comprueba—, y volver a digerir el recurso en cada
    // llamada era trabajo gratis: los hechos salen de las constantes.
    let profile = crate::technical_profile::f4_profile_facts(event.preparation());
    let mut document = json!({
        "event_id": event.event_id().as_str(),
        "individual_ids": event.individual_ids().iter().map(|id| id.as_str()).collect::<Vec<_>>(),
        "profile": {
            "id": profile.id(),
            "resource_path": profile.resource_path(),
            "sha256": profile.sha256(),
            "preparation": profile.preparation().as_str(),
            "insert_length": {
                "minimum_bases": profile.insert_length().minimum_bases(),
                "maximum_bases": profile.insert_length().maximum_bases(),
            },
        },
        "read_length": event.read_length(),
        "quality": {"acgt_phred": event.acgt_phred(), "n_phred": event.n_phred()},
        "sequencing_errors": {
            "requested_mode": event.sequencing_errors().as_str(),
            "effective_mode": event.sequencing_errors().as_str(),
            "qualities_represent_realized_error_probability": event.sequencing_errors().qualities_represent_realized_error_probability(),
        },
        "artifacts": {
            "chimera": {
                "rate": event.chimera_rate(),
            },
            "index_hopping": {
                "rate": event.index_hopping_rate(),
            },
            "pcr": event.pcr().map_or_else(
                || json!({"applicable": event.pcr_duplicates_applicable()}),
                crate::pcr_amplification::effective_settings_document,
            ),
            "optical_duplicates": {
                "realized_pairs": event.realized_optical_duplicate_pairs(),
            },
            "gc": {
                "mode": event
                    .pcr()
                    .map_or("neutral", crate::pcr_amplification::PcrSettings::gc_mode),
            },
            "adapters": {"enabled": false},
        },
        "biological_contamination": {
            "percent_range": {
                "minimum": event.biological_contamination().minimum_percent(),
                "maximum": event.biological_contamination().maximum_percent(),
            },
            "recipient_ids": resolved_contamination_recipients(event),
        },
        // El catálogo es el **más alto** cuya etapa abre este evento. La
        // amplificación estrena el v5 y no desplaza ninguna decisión: las
        // trece heredadas siguen en el v3 y las dos del tartamudeo en el v4.
        "rng_catalog": if event.pcr().is_some() {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V5
        } else if stutter_active {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
        } else {
            crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
        },
        "id_catalog": f4_id_catalog_value(),
        "libraries": resolved_technical_libraries(config, event),
        "target": target,
    });
    // La clave se **añade** solo cuando el evento abre de verdad la etapa
    // `stutter`. Escribirla siempre, aunque fuera con el interruptor apagado,
    // cambiaría los bytes de toda configuración efectiva ya publicada.
    if stutter_active && let Some(object) = document.as_object_mut() {
        object.insert(
            "stutter".to_owned(),
            crate::stutter::effective_settings_document(event.stutter()),
        );
    }
    document
}

fn resolved_technical_libraries(
    config: &EffectiveConfig,
    event: &EffectiveTechnicalEvent,
) -> Vec<Value> {
    let reference_id = ReferenceId::parse(config.simulation_reference_id(), "effective_config")
        .expect("effective configuration preserves a valid reference ID");
    let identities = derive_technical_identity_plan(
        config.reproducibility().root_seed(),
        &reference_id,
        event.event_id().clone(),
        event.individual_ids(),
    )
    .expect("effective technical event preserves derivable identities");
    let recipients = resolved_contamination_recipients(event)
        .into_iter()
        .filter_map(|value| value.as_str().map(str::to_owned))
        .collect::<std::collections::BTreeSet<_>>();
    identities
        .libraries()
        .map(|identity| {
            let configured = event
                .libraries()
                .iter()
                .find(|library| library.individual_id() == identity.individual_id())
                .expect("effective event has one configured library per individual");
            json!({
                "individual_id": identity.individual_id().as_str(),
                "library_id": identity.library_id().as_str(),
                "i7": configured.i7(),
                "i5": configured.i5(),
                "is_contamination_recipient": recipients.contains(identity.individual_id().as_str()),
            })
        })
        .collect()
}

fn resolved_contamination_recipients(event: &EffectiveTechnicalEvent) -> Vec<Value> {
    event
        .biological_contamination()
        .configured_recipients()
        .map(|recipients| recipients.to_vec())
        .unwrap_or_else(|| event.individual_ids().to_vec())
        .into_iter()
        .map(|id| json!(id.as_str()))
        .collect()
}

fn f4_id_catalog_value() -> Value {
    json!({
        "version": crate::identity::ID_CATALOG_VERSION,
        "derivations": {
            "v1": ["execution", "reference", "individual", "event", "batch", "library", "read_group", "template"],
            "v2": ["PreHopInstanceId", "TechnicalPairId"],
        }
    })
}

/// Bloque `profile` que corresponde a un perfil incorporado.
///
/// Se construye **desde el perfil resuelto**, no desde literales repetidos:
/// así la ruta, la huella y la preparación que el documento publica son las
/// mismas que el validador del recurso ya comprobó, y no pueden separarse.
fn f4_technical_profile_object(profile_id: &str) -> Option<Value> {
    let profile = crate::technical_profile::f4_profile_facts_by_id(profile_id)?;
    Some(json!({
        "id": profile.id(),
        "resource_path": profile.resource_path(),
        "sha256": profile.sha256(),
        "preparation": profile.preparation().as_str(),
        "insert_length": {
            "minimum_bases": profile.insert_length().minimum_bases(),
            "maximum_bases": profile.insert_length().maximum_bases(),
        },
    }))
}

pub fn effective_config_document_with_population_snvs(
    config: &EffectiveConfig,
    authorizations: crate::execution::ExecutionAuthorizations,
    output: ReferenceOutput<'_>,
    reductions: &[ReducedReferenceOutput<'_>],
    population_snvs: Option<&PopulationSnvsManifestLink>,
) -> Result<Value, Diagnostic> {
    let mut document = effective_config_document_with_reductions(config, output, reductions)?;
    let object = document.as_object_mut().ok_or_else(|| {
        Diagnostic::internal_invariant(
            "effective-config",
            "effective_config_document_is_not_object",
        )
    })?;
    object.insert(
        "authorizations".to_owned(),
        json!({
            "allow_empty_snv_matrix": authorizations.allow_empty_snv_matrix,
            "allow_incompatible_genotypes": authorizations.allow_incompatible_genotypes,
        }),
    );
    object.insert(
        "population_snvs".to_owned(),
        population_snvs.map_or(Value::Null, population_snvs_link_json),
    );
    Ok(document)
}

#[must_use]
pub fn general_manifest_document(manifest: &GeneralManifest<'_>) -> Value {
    let mut document = json!({
        "manifest_kind": "gps_execution_manifest",
        "schema_version": if manifest.repeat_loci.is_some() {
            GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS
        } else {
            GENERAL_MANIFEST_SCHEMA_VERSION
        },
        "execution_id": manifest.execution_id,
        "state": manifest.state.as_str(),
        "configuration_input": manifest.configuration_input.map(|input| json!({
            "kind": input.kind,
            "original_path": input.original_path,
            "artifact_id": input.copy.id,
        })),
        "effective_config_artifact_id": manifest.effective_config.map(|artifact| artifact.id.as_str()),
        "references": manifest.references.iter().map(manifest_reference_json).collect::<Vec<_>>(),
        "reference_reductions": manifest.reduced_references.iter().map(manifest_reduction_json).collect::<Vec<_>>(),
        "population_snvs": manifest.population_snvs.map(population_snvs_link_json),
        "vcf_outputs": manifest
            .vcf_outputs
            .iter()
            .map(|link| vcf_output_json(link, manifest.repeat_loci.is_some()))
            .collect::<Vec<_>>(),
        "authorizations": {
            "allow_empty_snv_matrix": manifest.authorizations.allow_empty_snv_matrix,
            "allow_incompatible_genotypes": manifest.authorizations.allow_incompatible_genotypes,
        },
        "output_categories": output_categories(manifest),
        "artifacts": manifest.artifacts.iter().map(ArtifactRecord::to_json).collect::<Vec<_>>(),
        "diagnostics": manifest.diagnostics.iter().map(|diagnostic| json!({
            "code": diagnostic.code(),
            "message": diagnostic.message(),
        })).collect::<Vec<_>>(),
        "technical_events": manifest.events.map_or_else(Vec::new, |events| {
            events.iter().map(|link| technical_event_manifest_json(link, manifest)).collect()
        }),
    });
    if let Some(events) = manifest.events {
        document["events"] = Value::Array(events.iter().map(manifest_event_json).collect());
    }
    if let Some(repeat_loci) = manifest.repeat_loci {
        document["repeat_loci"] = repeat_loci.clone();
    }
    document
}

fn technical_event_manifest_json(
    link: &ManifestEventLink,
    manifest: &GeneralManifest<'_>,
) -> Value {
    let event_manifest = manifest
        .artifacts
        .iter()
        .find(|artifact| artifact.id == link.manifest_artifact_id);
    json!({
        "event_id": link.event_id,
        "effective_config": manifest.effective_config.map(ArtifactRecord::to_json),
        "profile": {
            "id": link.profile_id,
            "sha256": link.profile_sha256,
        },
        "rng_catalog": link.rng_catalog,
        "id_catalog": f4_id_catalog_value(),
        "artifacts": {
            "chimera": {
                "rate": link.chimera_rate,
                "realized_pairs": link.realized_chimera_pairs,
            }
        },
        "event_manifest": event_manifest.map(ArtifactRecord::to_json),
    })
}

fn manifest_event_json(link: &ManifestEventLink) -> Value {
    json!({
        "event_id": link.event_id,
        "state": link.state.as_str(),
        "manifest_artifact_id": link.manifest_artifact_id,
    })
}

/// La activación VCF de la configuración efectiva.
///
/// La versión 6 publica además la del algoritmo de normalización, que DTC-064
/// exige registrar; la 4, ya cerrada, no la lleva. Es la única forma de esta
/// sección, y la usan tanto quien escribe el documento como quien lo
/// reconstruye desde el TOML para compararlo.
fn effective_vcf_output_json(config: &EffectiveConfig) -> Value {
    config.vcf_output().map_or(Value::Null, |mode| {
        let mode = match mode {
            VcfOutputMode::Text => "text",
            VcfOutputMode::Bgzf => "bgzf",
        };
        if config.document_version() == DeclaredDocumentVersion::V4 {
            json!({"mode": mode, "normalization_version": VCF_NORMALIZATION_VERSION})
        } else {
            json!({"mode": mode})
        }
    })
}

/// Enlace VCF del manifiesto general. La versión 6 registra la versión del
/// algoritmo de normalización (DTC-064); la 4, ya cerrada, no.
fn vcf_output_json(link: &VcfManifestLink, with_normalization_version: bool) -> Value {
    let mut value = json!({
        "state": link.state.as_str(),
        "reference_id": link.reference_id,
        "samples": link.samples,
        "mode": match link.mode { VcfOutputMode::Text => "text", VcfOutputMode::Bgzf => "bgzf" },
        "vcf_artifact_id": link.vcf_artifact_id,
        "index_kind": link.index_kind.map(|kind| match kind { VcfIndexKind::Tbi => "tbi", VcfIndexKind::Csi => "csi" }),
        "index_artifact_id": link.index_artifact_id,
    });
    if with_normalization_version {
        value["normalization_version"] = json!(VCF_NORMALIZATION_VERSION);
    }
    value
}

fn population_snvs_link_json(link: &PopulationSnvsManifestLink) -> Value {
    let mut object = serde_json::Map::from_iter([
        ("mode".to_owned(), json!(link.mode)),
        ("metadata_sha256".to_owned(), json!(link.metadata_sha256)),
        (
            "metadata_artifact_id".to_owned(),
            json!(link.metadata_artifact_id),
        ),
        (
            "normalization_audit".to_owned(),
            link.normalization_audit.clone(),
        ),
    ]);
    if let Some(statistics) = &link.initial_population_statistics {
        object.insert(
            "initial_population_statistics".to_owned(),
            statistics.clone(),
        );
    }
    Value::Object(object)
}

pub fn write_json_document(
    path: &Path,
    relative_path: &Path,
    document: &Value,
) -> Result<ArtifactIdentity, Diagnostic> {
    let bytes = json_document_bytes(document)?;
    let mut file = File::create(path).map_err(|error| io_error(path, "create", error))?;
    file.write_all(&bytes)
        .map_err(|error| io_error(path, "write", error))?;
    file.flush()
        .map_err(|error| io_error(path, "flush", error))?;
    drop(file);
    identify_file(path, relative_path)
}

pub fn json_document_bytes(document: &Value) -> Result<Vec<u8>, Diagnostic> {
    let mut bytes = serde_json::to_vec_pretty(document).map_err(json_error)?;
    bytes.push(b'\n');
    Ok(bytes)
}

/// Analiza y valida el contrato propio y la resolución de todos los enlaces.
pub fn validate_general_manifest_bytes(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    parse_general_manifest_bytes(path, bytes).map(|_| ())
}

/// Relee cada manifiesto de evento enlazado por la raíz y comprueba que sus
/// hechos físicos siguen siendo exactamente los del catálogo general. Esta
/// puerta es deliberadamente independiente de la configuración efectiva para
/// que también proteja un fallo controlado que ya no la publica.
pub fn validate_published_event_manifests(
    root: &Path,
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    let records = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let document: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "la raíz no es un objeto"))?;
    let catalog = records
        .iter()
        .map(|record| (record.id.as_str(), record))
        .collect::<BTreeMap<_, _>>();
    let Some(events) = object.get("events") else {
        return Ok(());
    };
    let mut linked_event_artifacts = BTreeSet::new();
    let mut linked_fastq_artifacts = BTreeSet::new();
    for event in events
        .as_array()
        .ok_or_else(|| manifest_error(manifest_path, "events no es una lista"))?
    {
        let link = event
            .as_object()
            .ok_or_else(|| manifest_error(manifest_path, "un enlace de evento no es objeto"))?;
        let event_id = required_non_empty_text(link, "event_id", manifest_path)?;
        let state = required_non_empty_text(link, "state", manifest_path)?;
        let manifest_id = required_non_empty_text(link, "manifest_artifact_id", manifest_path)?;
        linked_event_artifacts.insert(manifest_id.to_owned());
        let artifact = catalog.get(manifest_id).ok_or_else(|| {
            manifest_error(manifest_path, "el artefacto del manifiesto de evento falta")
        })?;
        let event_path = root.join(&artifact.path);
        let observed = identify_file(&event_path, Path::new(&artifact.path))?;
        if observed.size_bytes() != artifact.size_bytes || observed.digest() != artifact.digest {
            return Err(manifest_error(
                manifest_path,
                "el manifiesto de evento no coincide con su identidad publicada",
            ));
        }
        let bytes = fs::read(&event_path)
            .map_err(|error| io_error(&event_path, "read_event_manifest", error))?;
        validate_event_manifest_bytes(Path::new(&artifact.path), &bytes)?;
        let event_document: Value = serde_json::from_slice(&bytes).map_err(|error| {
            manifest_error(
                manifest_path,
                format!("el manifiesto de evento no es JSON legible: {error}"),
            )
        })?;
        if event_document.get("event_id").and_then(Value::as_str) != Some(event_id)
            || event_document.get("state").and_then(Value::as_str) != Some(state)
        {
            return Err(manifest_error(
                manifest_path,
                "el enlace general no coincide con el manifiesto de evento",
            ));
        }
        let reference = event_document
            .get("reference")
            .and_then(Value::as_object)
            .and_then(|reference| reference.get("fasta"))
            .ok_or_else(|| manifest_error(manifest_path, "falta FASTA de evento"))?;
        validate_event_artifact_against_catalog(reference, &catalog, manifest_path)?;
        for artifact in event_document
            .get("outputs")
            .and_then(Value::as_array)
            .ok_or_else(|| manifest_error(manifest_path, "faltan salidas de evento"))?
        {
            if artifact.get("state").and_then(Value::as_str) == Some("complete") {
                let physical = artifact.get("artifact").ok_or_else(|| {
                    manifest_error(manifest_path, "artefacto completo de evento ausente")
                })?;
                validate_event_artifact_against_catalog(physical, &catalog, manifest_path)?;
                let id = physical.get("id").and_then(Value::as_str).ok_or_else(|| {
                    manifest_error(manifest_path, "artefacto completo de evento sin ID")
                })?;
                linked_fastq_artifacts.insert(id.to_owned());
            }
        }
    }
    let cataloged_event_artifacts = records
        .iter()
        .filter(|artifact| artifact.kind == "event_manifest")
        .map(|artifact| artifact.id.clone())
        .collect::<BTreeSet<_>>();
    let cataloged_fastq_artifacts = records
        .iter()
        .filter(|artifact| artifact.kind == "fastq")
        .map(|artifact| artifact.id.clone())
        .collect::<BTreeSet<_>>();
    if cataloged_event_artifacts != linked_event_artifacts
        || cataloged_fastq_artifacts != linked_fastq_artifacts
    {
        return Err(manifest_error(
            manifest_path,
            "el catálogo contiene o omite artefactos de evento respecto de sus manifiestos enlazados",
        ));
    }
    Ok(())
}

fn validate_event_artifact_against_catalog(
    value: &Value,
    catalog: &BTreeMap<&str, &ArtifactRecord>,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let object = value
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "artefacto de evento inválido"))?;
    let id = object
        .get("id")
        .and_then(Value::as_str)
        .ok_or_else(|| manifest_error(manifest_path, "artefacto de evento sin ID"))?;
    let record = catalog
        .get(id)
        .ok_or_else(|| manifest_error(manifest_path, "artefacto de evento no catalogado"))?;
    if value != &record.to_json() {
        return Err(manifest_error(
            manifest_path,
            "artefacto de evento no coincide con el catálogo físico general",
        ));
    }
    Ok(())
}

/// Resuelve una referencia exclusivamente desde su manifiesto general F1 y
/// comprueba que el FASTA enlazado sigue teniendo los bytes declarados.
///
/// Es la frontera reutilizable de DTC-060: el llamador aporta conjuntamente el
/// manifiesto y el `reference_id`; ni una ruta FASTA libre ni un ID aislado
/// pueden llegar al consumidor posterior.
///
/// `reused_snv_metadata` es la ruta de los metadatos de la matriz de SNV si el
/// paquete la reutilizó (`mode = "reuse"`); conviene una ruta absoluta, porque
/// una relativa se resuelve contra el directorio de trabajo. Ver
/// [`validate_published_resources`].
///
/// Esta función juzga además la verdad entera del paquete (DT-023, cerrada por
/// aceptación): su coste es el de ese juicio, y debe reconsiderarse si una orden
/// de la CLI la usa con paquetes grandes.
pub fn resolve_manifest_reference(
    manifest_path: &Path,
    reference_id: &ReferenceId,
    reused_snv_metadata: Option<&Path>,
) -> Result<ResolvedManifestReference, Diagnostic> {
    let bytes = fs::read(manifest_path).map_err(|error| {
        Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
            .with("path", manifest_path.display().to_string())
            .with("operation", "read_reference_manifest")
            .with("detail", error.to_string())
    })?;
    let artifacts = parse_general_manifest_bytes(manifest_path, &bytes)?;
    let parent = manifest_path.parent().ok_or_else(|| {
        manifest_error(
            manifest_path,
            "el manifiesto no tiene un directorio padre para resolver recursos",
        )
    })?;
    let document: Value = serde_json::from_slice(&bytes).map_err(|error| {
        manifest_error(
            manifest_path,
            format!("JSON no legible tras validación: {error}"),
        )
    })?;
    let references = document
        .get("references")
        .and_then(Value::as_array)
        .ok_or_else(|| manifest_error(manifest_path, "falta references tras validación"))?;
    let fasta_artifact_id = references
        .iter()
        .filter_map(Value::as_object)
        .find(|reference| {
            reference.get("id").and_then(Value::as_str) == Some(reference_id.as_str())
        })
        .and_then(|reference| reference.get("fasta_artifact_id").and_then(Value::as_str))
        .ok_or_else(|| {
            manifest_error(
                manifest_path,
                format!(
                    "la referencia solicitada no existe: {}",
                    reference_id.as_str()
                ),
            )
        })?;
    let fasta = artifacts
        .iter()
        .find(|artifact| artifact.id == fasta_artifact_id && artifact.kind == "reference_fasta")
        .cloned()
        .ok_or_else(|| {
            manifest_error(
                manifest_path,
                format!("el FASTA de referencia no existe: {fasta_artifact_id}"),
            )
        })?;
    let fasta_path = parent.join(&fasta.path);
    let observed = identify_file(&fasta_path, &fasta.path)?;
    if observed.size_bytes() != fasta.size_bytes || observed.digest() != fasta.digest {
        return Err(artifact_identity_error()
            .with("path", fasta_path.display().to_string())
            .with("expected_size", fasta.size_bytes.to_string())
            .with("observed_size", observed.size_bytes().to_string())
            .with("expected_digest", &fasta.digest)
            .with("observed_digest", observed.digest()));
    }

    // El digest del FASTA seleccionado no autentica por sí solo una referencia
    // reducida: su mapa puede haber sido sustituido conservando intactos el
    // FASTA y sus metadatos. Reutilizamos la puerta física F1 para reconstruir
    // todas las relaciones de la publicación antes de entregar cualquier
    // enlace al consumidor F2.
    let effective_id = document
        .get("effective_config_artifact_id")
        .and_then(Value::as_str)
        .ok_or_else(|| manifest_error(manifest_path, "falta effective_config_artifact_id"))?;
    let effective = artifacts
        .iter()
        .find(|artifact| artifact.id == effective_id && artifact.kind == "effective_config")
        .ok_or_else(|| {
            manifest_error(
                manifest_path,
                format!("la configuración efectiva no existe: {effective_id}"),
            )
        })?;
    let effective_path = parent.join(&effective.path);
    let effective_bytes = fs::read(&effective_path)
        .map_err(|error| io_error(&effective_path, "read_effective_config", error))?;
    validate_published_resources(
        parent,
        &effective_path,
        &effective_bytes,
        manifest_path,
        &bytes,
        reused_snv_metadata,
    )?;
    Ok(ResolvedManifestReference {
        reference_id: reference_id.clone(),
        fasta_path,
        fasta,
    })
}

/// Valida la configuración efectiva releyéndola y cruza todos sus enlaces con
/// el catálogo y los enlaces semánticos del manifiesto ya validado.
pub fn validate_effective_config_against_manifest_bytes(
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    // El esquema que el diagnóstico cita sale de la versión que el documento
    // **declara**, no de una constante. Denunciar siempre el de la versión 4
    // mandaba a mirar un fichero que no era el suyo en cuanto el escenario
    // declaraba repeticiones.
    let version = effective_schema_version(effective_path, effective_bytes)?;
    match version {
        EFFECTIVE_CONFIG_SCHEMA_VERSION | EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS => {
            validate_effective_config_v3_against_manifest_bytes(
                effective_path,
                effective_bytes,
                manifest_path,
                manifest_bytes,
            )
            .map_err(|error| {
                retag_effective_config_schema(
                    error,
                    effective_config_schema_file_for(version),
                    effective_path,
                )
            })
        }
        _ => unreachable!("effective_schema_version only returns supported versions"),
    }
}

/// El esquema que hay que citar en un diagnóstico de la configuración
/// efectiva, según la versión que el documento declara.
fn effective_config_schema_file_for(version: &str) -> &'static str {
    if version == EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS {
        "effective-config-v6.schema.json"
    } else {
        "effective-config-v4.schema.json"
    }
}

/// Reescribe el esquema citado, **solo** si el diagnóstico habla de la
/// configuración efectiva.
///
/// # Por qué las dos condiciones y no solo una
///
/// Esta validación mira **dos** documentos, y sus dos familias de diagnóstico
/// viajan por el mismo `Result`. Retiquetar por el código solo convertiría un
/// rechazo del manifiesto general en uno que cita el esquema de la
/// configuración efectiva: cambiaría un destino equivocado por otro. Se exige
/// que coincidan la **ruta** y la **familia** del esquema ya citado.
fn retag_effective_config_schema(
    error: Diagnostic,
    schema: &'static str,
    effective_path: &Path,
) -> Diagnostic {
    if error.code() != "GPS-FMT-006" {
        return error;
    }
    let context = error.context();
    let habla_de_la_efectiva = context.get("path").map(String::as_str)
        == Some(effective_path.display().to_string().as_str())
        && context
            .get("schema")
            .is_some_and(|nombre| nombre.starts_with("effective-config-"));
    if habla_de_la_efectiva {
        retag_format_schema(error, schema)
    } else {
        error
    }
}

fn validate_effective_config_v3_against_manifest_bytes(
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    let artifacts = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let manifest = manifest
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "la raíz no es un objeto"))?;
    if !is_supported_general_manifest_version(
        manifest.get("schema_version").and_then(Value::as_str),
    ) {
        return Err(effective_error(
            effective_path,
            "el manifiesto no usa la versión general vigente",
        ));
    }
    let catalog = artifacts
        .iter()
        .map(|record| (record.id.as_str(), record.kind.as_str()))
        .collect::<BTreeMap<_, _>>();
    let effective_id = optional_text(manifest, "effective_config_artifact_id", manifest_path)?
        .ok_or_else(|| {
            effective_error(
                effective_path,
                "el manifiesto no enlaza la configuración efectiva",
            )
        })?;
    require_artifact_kind(&catalog, effective_id, "effective_config", manifest_path)?;
    let manifest_references = validate_references(
        required_array(manifest, "references", manifest_path)?,
        &catalog,
        manifest_path,
    )?;

    let document: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "la raíz no es un objeto"))?;
    const LEGACY_FIELDS: [&str; 11] = [
        "effective_config_kind",
        "schema_version",
        "source_path",
        "ploidy",
        "normalization",
        "seed",
        "reference",
        "population_snvs",
        "authorizations",
        "reduced_references",
        "initial_populations",
    ];
    let mut fields = LEGACY_FIELDS.to_vec();
    for field in ["direct_events", "vcf_output", "technical_event"] {
        if object.contains_key(field) {
            fields.push(field);
        }
    }
    // La sección de repeticiones **pertenece a la versión 6 y solo a ella**: su
    // presencia en un documento v4 y su ausencia en uno v6 son errores de
    // forma, no campos opcionales. Por eso se añaden a la lista cerrada según
    // la versión declarada y no según si aparecen.
    let declared_version = object.get("schema_version").and_then(Value::as_str);
    if declared_version == Some(EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS) {
        fields.push("repeat_loci");
        fields.push("repeat_admission");
        fields.push("repeat_output_budget");
        // Opcional dentro de la versión 6: solo cuando un locus inserta su
        // tracto (DTC-098). Esta lista es exacta, así que entra solo si está.
        if object.contains_key("inserted_reference") {
            fields.push("inserted_reference");
        }
    }
    // El nombre del esquema que se denuncia es el de la versión que el documento
    // declara. Antes se denunciaba siempre el de la versión 4, de modo que una
    // configuración efectiva v6 defectuosa mandaba a mirar un esquema que no era
    // el suyo.
    validate_effective_closed_object(object, &fields, "effective-config.json", effective_path)
        .map_err(|error| retag_format_schema(error, effective_schema_file_for(declared_version)))?;
    // **Los dos documentos tienen que decir lo mismo sobre las repeticiones.**
    //
    // La condición no es «configuración efectiva v6 exige manifiesto v6»: un
    // documento general que **puede** declarar loci publica configuración
    // efectiva v6 aunque no declare ninguno, y ese escenario no estrena
    // sección en el manifiesto ni debe hacerlo. Lo que no puede existir es una
    // configuración efectiva **con loci** frente a un manifiesto que no los
    // registra, ni al revés, y eso es exactamente lo que se comprueba.
    let effective_declares_loci = object
        .get("repeat_loci")
        .and_then(Value::as_array)
        .is_some_and(|loci| !loci.is_empty());
    let manifest_declares_loci = manifest.get("schema_version").and_then(Value::as_str)
        == Some(GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS);
    if effective_declares_loci != manifest_declares_loci {
        return Err(effective_error(
            effective_path,
            if effective_declares_loci {
                "la configuración efectiva declara loci de repetición y el manifiesto no los registra"
            } else {
                "el manifiesto registra loci de repetición que la configuración efectiva no declara"
            },
        ));
    }
    if object.get("effective_config_kind").and_then(Value::as_str) != Some("gps_effective_config")
        || !matches!(
            declared_version,
            Some(EFFECTIVE_CONFIG_SCHEMA_VERSION)
                | Some(EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS)
        )
    {
        return Err(effective_error(
            effective_path,
            "clase o versión incompatible",
        ));
    }
    let source_path = effective_non_empty_text(object, "source_path", effective_path)?;
    let effective_ploidy = effective_positive_u64(object, "ploidy", effective_path)?;
    if effective_ploidy > 10 {
        return Err(effective_error(effective_path, "ploidy fuera de rango"));
    }
    validate_effective_normalization(
        effective_object(object, "normalization", effective_path)?,
        effective_path,
    )?;
    let configuration_input = required_object(manifest, "configuration_input", manifest_path)?;
    if required_non_empty_text(configuration_input, "original_path", manifest_path)? != source_path
    {
        return Err(effective_error(
            effective_path,
            "source_path no coincide con la procedencia de la copia TOML",
        ));
    }
    let seed = effective_object(object, "seed", effective_path)?;
    validate_effective_seed(seed, effective_path)?;
    let effective_root_seed = effective_non_empty_text(seed, "root", effective_path)?;
    let effective_reference = effective_object(object, "reference", effective_path)?;
    let reference = validate_effective_reference(effective_reference, &catalog, effective_path)?;
    let eligible_bases =
        effective_positive_u64(effective_reference, "total_genome_size", effective_path)?;
    if manifest_references.get(&reference.0).map(String::as_str) != Some(reference.1.as_str()) {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no coincide con el enlace del manifiesto",
        ));
    }
    let mut effective_references = BTreeMap::from([(reference.0.clone(), reference.1.clone())]);
    let inserted = validate_effective_inserted_reference(
        object,
        manifest,
        effective_reference,
        &reference.1,
        effective_path,
    )?;
    if let Some((id, fasta)) = &inserted {
        effective_references.insert(id.clone(), fasta.clone());
    }
    // El VCF habla de la referencia sobre la que se simula: la final si hay
    // inserción.
    let simulation_reference = object
        .get("inserted_reference")
        .and_then(Value::as_object)
        .unwrap_or(effective_reference);
    validate_effective_population_snvs(object, manifest, &catalog, effective_path, manifest_path)?;
    validate_effective_authorizations(object, manifest, effective_path, manifest_path)?;
    let initial_individual_ids = validate_effective_initial_populations(object, effective_path)?;
    validate_effective_vcf_against_manifest(
        object,
        manifest,
        &catalog,
        simulation_reference,
        effective_path,
        manifest_path,
    )?;
    let direct_event_ids = if object.contains_key("direct_events") {
        validate_effective_direct_events(object, effective_path)?
    } else {
        BTreeSet::new()
    };
    if object.contains_key("technical_event") {
        validate_effective_technical_event(
            object,
            &initial_individual_ids,
            &direct_event_ids,
            object.get("schema_version").and_then(Value::as_str)
                == Some(EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS),
            effective_path,
        )?;
    }
    validate_effective_technical_event_link(object, manifest, effective_path, manifest_path)?;

    let manifest_reductions = required_array(manifest, "reference_reductions", manifest_path)?
        .iter()
        .map(|value| {
            let reduction = value
                .as_object()
                .ok_or_else(|| manifest_error(manifest_path, "una reducción no es objeto"))?;
            let id = required_non_empty_text(reduction, "id", manifest_path)?;
            Ok((id.to_owned(), reduction))
        })
        .collect::<Result<BTreeMap<_, _>, Diagnostic>>()?;
    let effective_reductions = effective_array(object, "reduced_references", effective_path)?;
    if effective_reductions.len() != manifest_reductions.len() {
        return Err(effective_error(
            effective_path,
            "las reducciones efectivas no coinciden con el manifiesto",
        ));
    }
    let mut effective_reduction_ids = BTreeSet::new();
    for value in effective_reductions {
        let reduction = validate_effective_reduction(value, &catalog, effective_path)?;
        if !effective_reduction_ids.insert(reduction.id.clone()) {
            return Err(effective_error(
                effective_path,
                "ID de reducción efectiva repetido",
            ));
        }
        let manifest_reduction = manifest_reductions.get(&reduction.id).ok_or_else(|| {
            effective_error(
                effective_path,
                format!("la reducción {} no existe en el manifiesto", reduction.id),
            )
        })?;
        let effective_reduction = value
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        if manifest_reduction.get("state").and_then(Value::as_str) != Some("complete")
            || manifest_reduction
                .get("eligible_bases")
                .and_then(Value::as_u64)
                != Some(eligible_bases)
            || manifest_reduction.get("root_seed").and_then(Value::as_str)
                != Some(effective_root_seed)
            || manifest_reduction
                .get("map_document_kind")
                .and_then(Value::as_str)
                != Some(REDUCTION_MAP_DOCUMENT_KIND)
            || manifest_reduction
                .get("map_schema_version")
                .and_then(Value::as_str)
                != Some(REDUCTION_MAP_SCHEMA_VERSION)
            || manifest_reduction
                .get("source_reference_id")
                .and_then(Value::as_str)
                != Some(reduction.source_reference_id.as_str())
            || manifest_reduction
                .get("fasta_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.fasta_artifact_id.as_str())
            || manifest_reduction
                .get("map_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.map_artifact_id.as_str())
            || manifest_reduction
                .get("source_fasta_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.source_fasta_artifact_id.as_str())
            || [
                "region_size_bases",
                "target",
                "selection_subseed",
                "selection_order",
            ]
            .iter()
            .any(|field| manifest_reduction.get(*field) != effective_reduction.get(*field))
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "la reducción efectiva {} no coincide con el manifiesto",
                    reduction.id
                ),
            ));
        }
        if effective_references
            .insert(reduction.id.clone(), reduction.fasta_artifact_id.clone())
            .is_some()
        {
            return Err(effective_error(
                effective_path,
                "una reducción reutiliza un ID de referencia",
            ));
        }
    }
    if effective_references != manifest_references {
        return Err(effective_error(
            effective_path,
            "el conjunto de referencias del manifiesto no coincide exactamente con la configuración efectiva",
        ));
    }
    Ok(())
}

/// **La sección `inserted_reference` de la configuración efectiva** frente a
/// la `reference_insertion` del manifiesto (DTC-098).
///
/// Las dos secciones existen a la vez o no existe ninguna; dicen lo mismo de
/// la referencia final, de su origen, de sus artefactos y de los cambios; el
/// origen es la referencia declarada con su FASTA; y el catálogo final suma
/// exactamente las bases de la suministrada más las insertadas. Devuelve el
/// identificador de la final y su FASTA, que entran en el conjunto de
/// referencias que el manifiesto tiene que registrar.
fn validate_effective_inserted_reference(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    reference: &serde_json::Map<String, Value>,
    reference_fasta: &str,
    effective_path: &Path,
) -> Result<Option<(String, String)>, Diagnostic> {
    const FIELDS: [&str; 7] = [
        "id",
        "source_reference_id",
        "sequences",
        "insertions",
        "fasta_artifact_id",
        "map_artifact_id",
        "source_fasta_artifact_id",
    ];
    let (value, section) = match (
        effective.get("inserted_reference"),
        manifest
            .get("reference_insertion")
            .and_then(Value::as_object),
    ) {
        (None, None) => return Ok(None),
        (Some(value), Some(section)) => (value, section),
        _ => {
            return Err(effective_error(
                effective_path,
                "la configuración efectiva y el manifiesto no coinciden en si hay referencia final por inserción",
            ));
        }
    };
    let inserted = value
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "inserted_reference no es objeto"))?;
    validate_effective_closed_object(inserted, &FIELDS, "inserted_reference", effective_path)?;
    let id = effective_non_empty_text(inserted, "id", effective_path)?;
    let fasta = effective_non_empty_text(inserted, "fasta_artifact_id", effective_path)?;
    if inserted.get("source_reference_id") != reference.get("id")
        || inserted
            .get("source_fasta_artifact_id")
            .and_then(Value::as_str)
            != Some(reference_fasta)
        || [
            ("id", "id"),
            ("source_reference_id", "source_reference_id"),
            ("fasta_artifact_id", "fasta_artifact_id"),
            ("map_artifact_id", "map_artifact_id"),
            ("source_fasta_artifact_id", "source_fasta_artifact_id"),
            ("insertions", "insertions"),
        ]
        .iter()
        .any(|(mine, theirs)| inserted.get(*mine) != section.get(*theirs))
    {
        return Err(effective_error(
            effective_path,
            "inserted_reference no coincide con la referencia declarada o con reference_insertion del manifiesto",
        ));
    }
    let inserted_bases = section
        .get("insertions")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|change| change.get("inserted_bases").and_then(Value::as_u64))
        .try_fold(0_u64, u64::checked_add);
    let expected_total = inserted_bases
        .zip(reference.get("total_genome_size").and_then(Value::as_u64))
        .and_then(|(added, total)| total.checked_add(added))
        // Un tamaño que no se puede calcular no lo alcanza ningún catálogo, y
        // lo rechaza la comprobación del catálogo de abajo.
        .unwrap_or(u64::MAX);
    validate_effective_derived_sequences(inserted, expected_total, effective_path)?;
    Ok(Some((id.to_owned(), fasta.to_owned())))
}

fn validate_effective_technical_event_link(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    if !effective.contains_key("technical_event") {
        return match manifest.get("events") {
            None => Ok(()),
            Some(Value::Array(events)) if events.is_empty() => Ok(()),
            _ => Err(effective_error(
                effective_path,
                "un manifiesto con eventos exige technical_event en la configuración efectiva",
            )),
        };
    }
    let event = effective_object(effective, "technical_event", effective_path)?;
    let event_id = effective_non_empty_text(event, "event_id", effective_path)?;
    let events = manifest
        .get("events")
        .and_then(Value::as_array)
        .ok_or_else(|| {
            effective_error(
                effective_path,
                "technical_event no tiene enlace en el manifiesto general",
            )
        })?;
    if events.len() != 1 {
        return Err(effective_error(
            effective_path,
            "F3 exige exactamente un enlace de evento técnico",
        ));
    }
    let link = events[0]
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "enlace de evento inválido"))?;
    if link.get("event_id").and_then(Value::as_str) != Some(event_id)
        || link.get("state").and_then(Value::as_str) != Some("complete")
    {
        return Err(effective_error(
            effective_path,
            "el enlace del evento no coincide con technical_event efectivo",
        ));
    }
    let manifest_id = link
        .get("manifest_artifact_id")
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(effective_path, "enlace de evento sin manifiesto"))?;
    let artifacts = required_array(manifest, "artifacts", manifest_path)?;
    if !artifacts.iter().any(|artifact| {
        artifact.get("id").and_then(Value::as_str) == Some(manifest_id)
            && artifact.get("kind").and_then(Value::as_str) == Some("event_manifest")
    }) {
        return Err(effective_error(
            effective_path,
            "el manifiesto de evento enlazado no está catalogado",
        ));
    }
    Ok(())
}

/// Relee los artefactos físicos enlazados y reconstruye sus relaciones sin
/// aceptar effective/manifest como prueba de sus propias declaraciones.
///
/// `reused_snv_metadata` es la ubicación de la matriz de SNV reutilizada: la
/// ruta de su `population_snvs.metadata.json`, con el TSV al lado. **Es
/// obligatorio si y solo si el paquete es de `mode = "reuse"`**, porque la
/// matriz reutilizada no viaja en el paquete y su ruta, escrita en el TOML
/// copiado, es relativa al directorio de trabajo de la ejecución original; el
/// resultado del juez no puede depender de dónde se lance. Sin reutilización el
/// parámetro se ignora. El juez acepta la matriz indicada solo si coinciden
/// las huellas: `metadata_sha256` del manifiesto y `matrix.sha256` de los
/// metadatos. **Conviene una ruta absoluta**: una relativa se resuelve contra
/// el directorio de trabajo del proceso que llama.
pub fn validate_published_resources(
    root: &Path,
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
    reused_snv_metadata: Option<&Path>,
) -> Result<(), Diagnostic> {
    match effective_schema_version(effective_path, effective_bytes)? {
        EFFECTIVE_CONFIG_SCHEMA_VERSION | EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS => {
            validate_published_resources_v3(
                root,
                effective_path,
                effective_bytes,
                manifest_path,
                manifest_bytes,
                reused_snv_metadata,
            )
        }
        _ => unreachable!("effective_schema_version only returns supported versions"),
    }
}

fn validate_published_resources_v3(
    root: &Path,
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
    reused_snv_metadata: Option<&Path>,
) -> Result<(), Diagnostic> {
    validate_published_event_manifests(root, manifest_path, manifest_bytes)?;
    validate_effective_config_against_manifest_bytes(
        effective_path,
        effective_bytes,
        manifest_path,
        manifest_bytes,
    )?;
    let artifacts = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let artifact_paths = artifacts
        .iter()
        .map(|artifact| (artifact.id.as_str(), artifact.path.as_str()))
        .collect::<BTreeMap<_, _>>();
    let effective: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let effective = effective
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "la raíz no es objeto"))?;
    let physical_config = validate_effective_targets_against_configuration(
        root,
        effective,
        &artifact_paths,
        manifest_path,
        manifest_bytes,
        effective_path,
    )?;
    let root_seed = RootSeed::parse(
        "effective-config.seed.root",
        effective["seed"]["root"]
            .as_str()
            .ok_or_else(|| effective_error(effective_path, "seed.root inválida"))?,
    )
    .map_err(|error| physical_error(effective_path, "root_seed", &error))?;
    let reference_value = effective_object(effective, "reference", effective_path)?;
    let declared_reference = read_physical_reference(
        root,
        reference_value,
        &artifact_paths,
        effective_path,
        if effective_needs_reference_bases(effective)
            || effective.contains_key("inserted_reference")
        {
            CatalogBases::Materialized
        } else {
            CatalogBases::SummaryOnly
        },
    )?;
    // **Con inserción, todo lo demás se juzga contra la referencia final**, y
    // solo después de haberla rehecho por cuenta propia (DTC-098).
    let inserted_reference = judge_inserted_reference(
        root,
        effective,
        manifest_bytes,
        &declared_reference,
        &artifact_paths,
        effective_path,
    )?;
    let source_reference = inserted_reference.as_ref().unwrap_or(&declared_reference);
    validate_published_reference_names(
        root,
        effective,
        manifest_bytes,
        &artifacts,
        source_reference,
        effective_path,
    )?;
    // **El generador de alelos por deslizamiento corre también aquí.** Este
    // juez reconstruye la configuración desde el TOML copiado, donde un locus
    // generado **no trae alelos**: son el resultado de aplicar el modelo a las
    // bases de la referencia. Sin regenerarlos, el juez compararía las
    // estadísticas publicadas contra una configuración sin un solo alelo y no
    // contrastaría nada.
    //
    // Y es la comparación correcta por sí misma: el abanico se rehace desde la
    // **referencia publicada** y la semilla publicada, de modo que un productor
    // que hubiera escrito otros alelos en el documento muere aquí.
    // **Primero el juez, y después la reconstrucción.** El orden importa: la
    // reconstrucción llama al productor, de modo que si fuera lo único que hay
    // el abanico se estaría comparando consigo mismo. El juez de abajo rehace
    // el abanico **por su cuenta**, desde los bytes publicados, y contrasta lo
    // escrito con lo rehecho.
    validate_published_slippage_fans(effective, source_reference, root_seed, effective_path)?;
    let physical_config = physical_config
        .with_generated_repeat_alleles(source_reference)?
        .unwrap_or(physical_config);
    validate_physical_technical_event(
        root,
        effective,
        manifest_path,
        manifest_bytes,
        &artifacts,
        &root_seed,
        source_reference,
        effective_path,
    )?;
    validate_physical_vcf_output(
        root,
        &physical_config,
        manifest_bytes,
        &artifacts,
        source_reference,
        effective_path,
        manifest_path,
    )?;
    let allow_empty_snv_matrix = effective
        .get("authorizations")
        .and_then(|value| value.get("allow_empty_snv_matrix"))
        .and_then(Value::as_bool)
        .ok_or_else(|| effective_error(effective_path, "authorizations inválido"))?;
    let matrix = validate_physical_population_snvs(
        root,
        effective,
        &physical_config,
        &artifacts,
        source_reference,
        effective_path,
        allow_empty_snv_matrix,
        reused_snv_metadata,
    )?;
    validate_physical_initial_population_statistics(
        effective,
        &physical_config,
        source_reference,
        matrix.as_ref(),
        effective_path,
    )?;
    validate_physical_truth(
        root,
        &artifacts,
        &physical_config,
        source_reference,
        matrix.as_ref(),
        effective_path,
    )?;
    if reference_value.get("source").and_then(Value::as_str) == Some("simulated") {
        validate_generated_reference(
            reference_value,
            root_seed,
            source_reference,
            effective_path,
            effective.get("repeat_loci"),
        )
        // DT-017: el rechazo cita el esquema de la versión que el documento
        // declara, no el v4 que la constructora escribe fijo.
        .map_err(|error| {
            retag_format_schema(
                error,
                effective_schema_file_for(effective.get("schema_version").and_then(Value::as_str)),
            )
        })?;
    }
    drop(inserted_reference);
    let mut references = BTreeMap::from([(
        declared_reference.id().as_str().to_owned(),
        declared_reference,
    )]);

    for reduction_value in effective_array(effective, "reduced_references", effective_path)? {
        let reduction = reduction_value
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        let id = effective_non_empty_text(reduction, "id", effective_path)?;
        let source_id = effective_non_empty_text(reduction, "source_reference_id", effective_path)?;
        let source = references.get(source_id).ok_or_else(|| {
            effective_error(
                effective_path,
                format!("la referencia fuente física {source_id} no existe"),
            )
        })?;
        let derived_subseed = derive_subseed(
            root_seed,
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: source_id,
                derived_reference_id: id,
            },
        )
        .map_err(|error| physical_error(effective_path, "derive_reduction_subseed", &error))?;
        if effective_non_empty_text(reduction, "selection_subseed", effective_path)?
            != derived_subseed.to_hex()
        {
            return Err(effective_error(
                effective_path,
                format!("selection_subseed de {id} no reproduce el ámbito canónico"),
            ));
        }
        let eligible = source
            .sequences()
            .iter()
            .map(|sequence| {
                EligibleSequence::new(sequence.id().clone(), sequence.class(), sequence.length())
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| {
                effective_error(
                    effective_path,
                    format!("catálogo elegible inválido: {error:?}"),
                )
            })?;
        let target_bases = effective_positive_u64(
            effective_object(reduction, "target", effective_path)?,
            "bases",
            effective_path,
        )?;
        let region_size_bases =
            effective_positive_u64(reduction, "region_size_bases", effective_path)?;
        let mut rng = derived_subseed.local_rng();
        let selection = select_reduced_reference(
            source.id(),
            &eligible,
            target_bases,
            region_size_bases,
            &mut rng,
        )
        .map_err(|error| {
            effective_error(
                effective_path,
                format!("selección física no reproducible: {error:?}"),
            )
        })?;
        let derived_id = ReferenceId::parse(id, "published-resource-validation")
            .map_err(|error| physical_error(effective_path, "derived_reference_id", &error))?;
        let plan = materialize_selection(&selection, &derived_id, &eligible).map_err(|error| {
            effective_error(
                effective_path,
                format!("plan físico no materializable: {error:?}"),
            )
        })?;
        let declared_order = effective_array(reduction, "selection_order", effective_path)?
            .iter()
            .map(Value::as_u64)
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| effective_error(effective_path, "selection_order inválido"))?;
        let expected_order = plan
            .segments()
            .iter()
            .map(|segment| segment.selection_order())
            .collect::<Vec<_>>();
        if declared_order != expected_order {
            return Err(effective_error(
                effective_path,
                format!("selection_order de {id} no coincide con la subsemilla"),
            ));
        }
        let expected_derived = build_reduced_reference(
            source_for(source, None)
                .map_err(|error| physical_error(effective_path, "reduction_source", &error))?
                .as_ref(),
            &plan,
        )
        .map_err(|error| physical_error(effective_path, "build_reduced_reference", &error))?;
        let fasta_id = effective_non_empty_text(reduction, "fasta_artifact_id", effective_path)?;
        let derived = read_reference_with_expected(
            root,
            id,
            fasta_id,
            effective_array(reduction, "sequences", effective_path)?,
            expected_derived.sequences(),
            &artifact_paths,
            effective_path,
        )?;
        compare_reference_bases(&expected_derived, &derived, effective_path, id)?;
        let map_id = effective_non_empty_text(reduction, "map_artifact_id", effective_path)?;
        let map_relative = artifact_paths.get(map_id).ok_or_else(|| {
            effective_error(effective_path, format!("mapa físico ausente: {map_id}"))
        })?;
        let map_path = root.join(map_relative);
        let map_bytes = std::fs::read(&map_path)
            .map_err(|error| io_error(&map_path, "read_reduction_map", error))?;
        let map: Value = serde_json::from_slice(&map_bytes).map_err(|error| {
            effective_error(
                effective_path,
                format!("mapa {map_id} no es JSON legible: {error}"),
            )
        })?;
        if map != reduction_map_document(&plan) {
            return Err(effective_error(
                effective_path,
                format!("mapa físico {map_id} no coincide con subsemilla y selección"),
            ));
        }
        validate_reduction_map_document(
            source_for(source, None)
                .map_err(|error| physical_error(effective_path, "reduction_source", &error))?
                .as_ref(),
            &derived,
            &map,
            ReductionMapParameters {
                region_size_bases,
                target_bases,
            },
        )
        .map_err(|error| physical_error(effective_path, "validate_reduction_map", &error))?;
        references.insert(id.to_owned(), derived);
    }
    Ok(())
}

/// **Todo lo publicado nombra la referencia de la simulación** (DTC-098): la
/// final si hay inserción y la suministrada si no.
///
/// Es una comprobación barata, de nombres y longitudes, no una relectura de
/// la verdad ni del VCF: la cabecera de la verdad —identificador, huella del
/// contenido y huella del FASTA—, el identificador y la longitud de cada
/// fragmento de la verdad que los lleva, las líneas `##reference` y
/// `##contig` del VCF, y el identificador de la sección de loci del
/// manifiesto. Con dos referencias en la misma ejecución, que una salida
/// nombre la otra es justo el error que nadie más vería.
///
/// La huella del contenido solo se compara si la referencia se leyó con sus
/// bases; sin ellas no puede calcularse, y entonces la ejecución no publica
/// nada que dependa de las bases.
fn validate_published_reference_names(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    manifest_bytes: &[u8],
    artifacts: &[ArtifactRecord],
    reference: &Reference,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let id = reference.id().as_str();
    let lengths = reference
        .sequences()
        .iter()
        .map(|sequence| (sequence.id().as_str(), sequence.length()))
        .collect::<Vec<_>>();
    let fasta_digest = simulation_reference_object(effective, effective_path)?
        .get("fasta_artifact_id")
        .and_then(Value::as_str)
        .and_then(|fasta| artifacts.iter().find(|artifact| artifact.id == fasta))
        .map(|artifact| artifact.digest.as_str());
    let read = |artifact_id: &str| {
        artifacts
            .iter()
            .find(|artifact| artifact.id == artifact_id)
            .and_then(|artifact| std::fs::read(root.join(&artifact.path)).ok())
    };
    if let Some(truth) = read("truth-jsonl") {
        let text = String::from_utf8_lossy(&truth);
        let mut lines = text.lines();
        let header: Value = lines
            .next()
            .and_then(|line| serde_json::from_str(line).ok())
            .unwrap_or(Value::Null);
        if header.get("reference_id").and_then(Value::as_str) != Some(id) {
            return Err(effective_error(
                effective_path,
                format!("la cabecera de la verdad no nombra la referencia de la simulación, {id}"),
            ));
        }
        if let Ok(content) = crate::truth_validation::reference_content_sha256(reference)
            && header
                .get("reference_content_sha256")
                .and_then(Value::as_str)
                != Some(content.as_str())
        {
            return Err(effective_error(
                effective_path,
                format!("la huella del contenido en la cabecera de la verdad no es la de {id}"),
            ));
        }
        let digests = header
            .get("source_digests")
            .and_then(Value::as_array)
            .map_or(&[][..], Vec::as_slice)
            .iter()
            .filter(|digest| digest["source"]["kind"] == "reference")
            .collect::<Vec<_>>();
        if digests.len() != 1
            || digests[0]["source"]["source_id"].as_str() != Some(id)
            || digests[0]["sha256"].as_str() != fasta_digest
        {
            return Err(effective_error(
                effective_path,
                format!("la huella de referencia de la verdad no es la del FASTA de {id}"),
            ));
        }
        fn foreign_fragment(value: &Value, id: &str, lengths: &[(&str, u64)]) -> bool {
            match value {
                Value::Object(object) => {
                    (object.contains_key("sequence_length")
                        && (object.get("reference_id").and_then(Value::as_str) != Some(id)
                            || object
                                .get("sequence_id")
                                .and_then(Value::as_str)
                                .and_then(|sequence| {
                                    lengths.iter().find(|(known, _)| *known == sequence)
                                })
                                .map(|(_, length)| *length)
                                != object.get("sequence_length").and_then(Value::as_u64)))
                        || object
                            .values()
                            .any(|child| foreign_fragment(child, id, lengths))
                }
                Value::Array(items) => items.iter().any(|item| foreign_fragment(item, id, lengths)),
                _ => false,
            }
        }
        if lines.any(|line| {
            serde_json::from_str::<Value>(line)
                .map_or(true, |record| foreign_fragment(&record, id, &lengths))
        }) {
            return Err(effective_error(
                effective_path,
                format!(
                    "un fragmento de la verdad no está en {id} o no tiene la longitud de su secuencia"
                ),
            ));
        }
    }
    let vcf = match read(VCF_ARTIFACT_ID) {
        None => None,
        Some(bytes) if bytes.starts_with(&[0x1f, 0x8b]) => {
            // Un bloque que declara más de 64 KiB haría entrar en pánico al
            // lector (DT-030): se rechaza antes de dárselo.
            let mut cursor = std::io::Cursor::new(bytes);
            crate::bgzf_blocks::reject_oversized_bgzf_blocks(&mut cursor).map_err(|fault| {
                effective_error(effective_path, format!("VCF BGZF no legible: {fault}"))
            })?;
            let mut text = String::new();
            let _ = std::io::Read::read_to_string(
                &mut noodles_bgzf::io::Reader::new(cursor),
                &mut text,
            );
            Some(text)
        }
        Some(bytes) => Some(String::from_utf8_lossy(&bytes).into_owned()),
    };
    if let Some(vcf) = vcf {
        let header = vcf
            .lines()
            .take_while(|line| line.starts_with("##"))
            .collect::<Vec<_>>();
        if !header.contains(&format!("##reference={id}").as_str()) {
            return Err(effective_error(
                effective_path,
                format!("la cabecera del VCF no declara ##reference={id}"),
            ));
        }
        let contigs = header
            .iter()
            .filter(|line| line.starts_with("##contig="))
            .map(|line| (*line).to_owned())
            .collect::<Vec<_>>();
        let expected = lengths
            .iter()
            .map(|(sequence, length)| format!("##contig=<ID={sequence},length={length}>"))
            .collect::<Vec<_>>();
        if contigs != expected {
            return Err(effective_error(
                effective_path,
                format!("los ##contig del VCF no son el catálogo de {id}"),
            ));
        }
    }
    let manifest: Value = serde_json::from_slice(manifest_bytes).unwrap_or(Value::Null);
    if let Some(section) = manifest.get("repeat_loci")
        && section.get("reference_id").and_then(Value::as_str) != Some(id)
    {
        return Err(effective_error(
            effective_path,
            format!(
                "la sección de loci del manifiesto no nombra la referencia de la simulación, {id}"
            ),
        ));
    }
    Ok(())
}

/// **Juez independiente de la referencia final por inserción** (DTC-098).
///
/// No llama a `insert_tracts` ni a `coordinate_map_document`: todo lo que
/// compara lo rehace aquí, desde la copia suministrada publicada y lo que la
/// configuración efectiva declara de cada locus.
///
/// 1. Los cambios publicados en `inserted_reference` son los que salen de los
///    loci con `declared.insert_after`, con el tracto de cada uno construido
///    desde sus motivos y sus `reference_units`.
/// 2. El FASTA final es la suministrada con esos tractos empalmados, base a
///    base.
/// 3. El mapa es exactamente el que describe ese empalme: segmentos
///    contiguos, conservados con su intervalo de origen e insertados en su
///    punto. El esquema del mapa admite huecos, intervalos invertidos e
///    identificadores iguales, así que no basta con él.
/// 4. El intervalo publicado de cada locus, en la configuración efectiva y en
///    el manifiesto, es la traducción de lo declarado.
///
/// Devuelve la referencia final leída del FASTA publicado, que es contra la
/// que se juzga todo lo demás, o `None` si la ejecución no inserta nada.
fn judge_inserted_reference(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    manifest_bytes: &[u8],
    supplied: &Reference,
    artifact_paths: &BTreeMap<&str, &str>,
    effective_path: &Path,
) -> Result<Option<Reference>, Diagnostic> {
    // Que la sección sea un objeto cerrado y concuerde con el manifiesto ya lo
    // exigió `validate_effective_inserted_reference`, antes que este juez.
    let Some(inserted) = effective
        .get("inserted_reference")
        .and_then(Value::as_object)
    else {
        return Ok(None);
    };
    let derived_id = effective_non_empty_text(inserted, "id", effective_path)?;
    let loci = effective
        .get("repeat_loci")
        .and_then(Value::as_array)
        .map_or(&[][..], Vec::as_slice);
    // (1) Los tractos, desde los motivos y las unidades declarados.
    let catalog_index = |sequence: &str| {
        supplied
            .sequences()
            .iter()
            .position(|candidate| candidate.id().as_str() == sequence)
    };
    let mut points: Vec<(usize, u64, String, String, Vec<u8>)> = Vec::new();
    // **Cota física de lo que el juez reserva** (DT-025): los tractos insertados
    // caben en el FASTA final publicado, que los contiene. El tope `max_units`
    // se lee del mismo documento que las unidades y no acota nada por sí solo.
    // Sin el fichero, la cota es cero y el locus se rechaza.
    let fasta_size = effective_non_empty_text(inserted, "fasta_artifact_id", effective_path)
        .ok()
        .and_then(|artifact| artifact_paths.get(artifact))
        .and_then(|relative| std::fs::metadata(root.join(relative)).ok())
        .map_or(0, |metadata| metadata.len());
    let mut reserved_bases = 0_u64;
    for locus in loci {
        let Some(insert_after) = locus
            .get("declared")
            .and_then(|declared| declared.get("insert_after"))
            .and_then(Value::as_u64)
        else {
            continue;
        };
        let id = locus.get("id").and_then(Value::as_str).unwrap_or_default();
        let sequence = locus
            .get("sequence_id")
            .and_then(Value::as_str)
            .unwrap_or_default();
        let motifs = locus
            .get("motifs")
            .and_then(Value::as_array)
            .map_or(&[][..], Vec::as_slice);
        let units = locus
            .get("reference_units")
            .and_then(Value::as_array)
            .map_or(&[][..], Vec::as_slice);
        let index = catalog_index(sequence);
        // El punto tiene que dejar una base a cada lado (REQ-IND-40): si no,
        // el empalme de abajo no tendría dónde cortar.
        let anchored = index
            .and_then(|index| supplied.sequences().get(index))
            .is_some_and(|host| insert_after >= 1 && insert_after < host.length());
        // Las unidades no pueden pasar del tope publicado del locus: un
        // documento corrupto con una cifra enorme no puede hacer reservar al
        // juez gigabytes antes de rechazarlo.
        let total_units = units
            .iter()
            .try_fold(0_u64, |total, count| total.checked_add(count.as_u64()?));
        let within_cap = total_units
            .zip(locus.get("max_units").and_then(Value::as_u64))
            .is_some_and(|(total, cap)| total <= cap);
        // Σ unidades × longitud del motivo, comprobada, y acumulada entre loci:
        // todo lo que el juez va a reservar tiene que caber en el FASTA final.
        let locus_bases = motifs
            .iter()
            .zip(units)
            .try_fold(0_u64, |total, (motif, count)| {
                total.checked_add((motif.as_str()?.len() as u64).checked_mul(count.as_u64()?)?)
            });
        let fits_in_fasta = locus_bases
            .and_then(|bases| reserved_bases.checked_add(bases))
            .filter(|total| *total <= fasta_size);
        if let Some(total) = fits_in_fasta {
            reserved_bases = total;
        }
        let mut bases = Vec::new();
        let mut well_formed = anchored
            && within_cap
            && fits_in_fasta.is_some()
            && !motifs.is_empty()
            && motifs.len() == units.len();
        let checked = if well_formed { motifs.len() } else { 0 };
        for (motif, count) in motifs.iter().zip(units).take(checked) {
            match (motif.as_str(), count.as_u64()) {
                (Some(motif), Some(count)) if !motif.is_empty() && count > 0 => {
                    for _ in 0..count {
                        bases.extend_from_slice(motif.as_bytes());
                    }
                }
                _ => well_formed = false,
            }
        }
        let (Some(index), true) = (index, well_formed) else {
            return Err(effective_error(
                effective_path,
                format!(
                    "referencia final: el locus {id} inserta un tracto sin secuencia conocida, sin una base a cada lado, sin motivos y unidades por bloque, o con más unidades que su tope o más bases que el FASTA final"
                ),
            ));
        };
        points.push((
            index,
            insert_after,
            id.to_owned(),
            sequence.to_owned(),
            bases,
        ));
    }
    points.sort_by_key(|(index, point, ..)| (*index, *point));
    let expected_changes = points
        .iter()
        .map(|(_, insert_after, locus, sequence, bases)| {
            json!({
                "locus_id": locus,
                "sequence_id": sequence,
                "insert_after": insert_after,
                "inserted_bases": bases.len(),
            })
        })
        .collect::<Vec<_>>();
    if points.is_empty() || inserted.get("insertions") != Some(&Value::Array(expected_changes)) {
        return Err(effective_error(
            effective_path,
            "referencia final: los cambios publicados no son los que declaran los loci",
        ));
    }
    // (2) y (3): el empalme, y el mapa que lo describe, rehechos aquí.
    let mut sequences = Vec::with_capacity(supplied.sequences().len());
    let mut segments = Vec::new();
    for (index, sequence) in supplied.sequences().iter().enumerate() {
        // La copia suministrada se lee con sus bases siempre que la efectiva
        // declara `inserted_reference` (ver la lectura en el llamador), así que
        // un fallo aquí es una invariante rota y no un documento rechazable.
        let source_bases = sequence.sequence().as_bytes().map_err(|_| {
            Diagnostic::internal_invariant("results", "judge_inserted_reference")
                .with("detail", "la copia suministrada se leyó sin sus bases")
        })?;
        let length = sequence.length();
        let class = sequence_class(sequence.class());
        let mut spliced = Vec::with_capacity(source_bases.len());
        let (mut source_next, mut derived_next) = (1_u64, 1_u64);
        let retained = |from: u64, to: u64, derived_next: &mut u64, segments: &mut Vec<Value>| {
            if from > to {
                return;
            }
            let span = to - from;
            segments.push(json!({
                "materialization_order": segments.len() + 1,
                "segment_class": "retained",
                "sequence_id": sequence.id().as_str(),
                "sequence_class": class,
                "source_interval": {"start": from, "end": to},
                "derived_interval": {"start": *derived_next, "end": *derived_next + span},
                "orientation": "forward",
            }));
            *derived_next += span + 1;
        };
        for (_, insert_after, _, _, bases) in points.iter().filter(|point| point.0 == index) {
            let anchor = *insert_after;
            spliced.extend_from_slice(&source_bases[(source_next - 1) as usize..anchor as usize]);
            retained(source_next, anchor, &mut derived_next, &mut segments);
            source_next = anchor + 1;
            spliced.extend_from_slice(bases);
            let span = bases.len() as u64;
            segments.push(json!({
                "materialization_order": segments.len() + 1,
                "segment_class": "inserted",
                "sequence_id": sequence.id().as_str(),
                "sequence_class": class,
                "derived_interval": {"start": derived_next, "end": derived_next + span - 1},
                "orientation": "forward",
            }));
            derived_next += span;
        }
        spliced.extend_from_slice(&source_bases[(source_next - 1) as usize..]);
        retained(source_next, length, &mut derived_next, &mut segments);
        // Todo lo que se construye aquí sale de bases ya validadas: un fallo
        // sería una invariante rota, no un documento rechazable.
        let invariant = |what: &str, error: &dyn std::fmt::Debug| {
            Diagnostic::internal_invariant("results", "judge_inserted_reference")
                .with("detail", format!("{what}: {error:?}"))
        };
        let sequence_id = SequenceId::parse(sequence.id().as_str(), derived_id)
            .map_err(|error| invariant("identificador", &error))?;
        let header = FastaHeader::parse(sequence.id().as_str().as_bytes())
            .map_err(|error| invariant("cabecera", &error))?;
        let dna =
            DnaSequence::try_from_bytes(&spliced).map_err(|error| invariant("bases", &error))?;
        sequences.push(
            ReferenceSequence::new(sequence_id, header, sequence.class(), &dna)
                .map_err(|error| invariant("secuencia", &error))?,
        );
    }
    let expected = Reference::new(
        ReferenceId::parse(derived_id, "published-resource-validation")
            .map_err(|error| physical_error(effective_path, "inserted_reference_id", &error))?,
        ReferenceOrigin::TractInsertion {
            source_reference_id: supplied.id().clone(),
        },
        sequences,
    )
    .map_err(|error| {
        Diagnostic::internal_invariant("results", "judge_inserted_reference")
            .with("detail", format!("{error:?}"))
    })?;
    let physical = read_physical_sequence_catalog(
        root,
        derived_id,
        effective_non_empty_text(inserted, "fasta_artifact_id", effective_path)?,
        effective_array(inserted, "sequences", effective_path)?,
        "n_count",
        artifact_paths,
        effective_path,
        CatalogBases::Materialized,
    )?;
    compare_reference_bases(&expected, &physical, effective_path, derived_id)?;
    let map_id = effective_non_empty_text(inserted, "map_artifact_id", effective_path)?;
    // Un mapa ausente o ilegible no es el esperado: lo rechaza la misma
    // comparación de abajo. Que el artefacto exista con su clase ya lo exige
    // el manifiesto.
    let map = artifact_paths
        .get(map_id)
        .and_then(|relative| std::fs::read(root.join(relative)).ok())
        .and_then(|bytes| serde_json::from_slice::<Value>(&bytes).ok())
        .unwrap_or(Value::Null);
    let expected_map = json!({
        "document_kind": "reference_coordinate_map",
        "schema_version": "2",
        "operation": "tract_insertion",
        "source_reference_id": supplied.id().as_str(),
        "derived_reference_id": derived_id,
        "segments": segments,
    });
    if map != expected_map {
        return Err(effective_error(
            effective_path,
            format!(
                "referencia final: el mapa {map_id} no describe el empalme de los tractos declarados"
            ),
        ));
    }
    // (4) El intervalo publicado de cada locus es la traducción de lo declarado.
    let shift = |sequence: &str, before: u64| -> u64 {
        points
            .iter()
            .filter(|point| point.3 == sequence && point.1 < before)
            .map(|point| point.4.len() as u64)
            .sum()
    };
    let manifest: Value = serde_json::from_slice(manifest_bytes).unwrap_or(Value::Null);
    let manifest_loci = manifest
        .get("repeat_loci")
        .and_then(|section| section.get("loci"))
        .and_then(Value::as_array)
        .map_or(&[][..], Vec::as_slice);
    for locus in loci {
        let id = locus.get("id").and_then(Value::as_str).unwrap_or_default();
        let sequence = locus
            .get("sequence_id")
            .and_then(Value::as_str)
            .unwrap_or_default();
        let declared = locus.get("declared");
        let expected_interval = match (
            declared
                .and_then(|value| value.get("insert_after"))
                .and_then(Value::as_u64),
            declared
                .and_then(|value| value.get("start"))
                .and_then(Value::as_u64),
            declared
                .and_then(|value| value.get("end"))
                .and_then(Value::as_u64),
        ) {
            (Some(anchor), None, None) => points.iter().find(|point| point.2 == id).map(|point| {
                let start = anchor + shift(sequence, anchor) + 1;
                (start, start + point.4.len() as u64 - 1)
            }),
            (None, Some(start), Some(end))
                if !points
                    .iter()
                    .any(|point| point.3 == sequence && start <= point.1 && point.1 < end) =>
            {
                Some((start + shift(sequence, start), end + shift(sequence, end)))
            }
            _ => None,
        };
        let published = |object: &Value| {
            (
                object.get("start").and_then(Value::as_u64),
                object.get("end").and_then(Value::as_u64),
            )
        };
        let in_manifest = manifest_loci
            .iter()
            .find(|candidate| candidate.get("id").and_then(Value::as_str) == Some(id))
            .map(published);
        // Lo declarado frente al TOML ya lo comprobó la resolución, antes de
        // este juez, así que aquí basta una rama: o el intervalo publicado es
        // la traducción de lo declarado, o no lo es.
        let translated = expected_interval.map(|(start, end)| (Some(start), Some(end)));
        if translated.is_none() || Some(published(locus)) != translated || in_manifest != translated
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "referencia final: el intervalo publicado del locus {id} no es la traducción al genoma final de lo declarado"
                ),
            ));
        }
    }
    Ok(Some(physical))
}

#[allow(clippy::too_many_arguments)]
fn validate_physical_technical_event(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    manifest_path: &Path,
    manifest_bytes: &[u8],
    artifacts: &[ArtifactRecord],
    root_seed: &RootSeed,
    reference: &Reference,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let Some(technical_event) = effective.get("technical_event") else {
        return Ok(());
    };
    let technical_event = technical_event
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "technical_event inválido"))?;
    let event_id = EventId::parse(
        effective_non_empty_text(technical_event, "event_id", effective_path)?,
        "effective_config",
    )
    .map_err(|_| effective_error(effective_path, "event_id técnico inválido"))?;
    let individuals = effective_array(technical_event, "individual_ids", effective_path)?
        .iter()
        .map(|value| {
            let value = value
                .as_str()
                .ok_or_else(|| effective_error(effective_path, "individual_id técnico inválido"))?;
            IndividualId::parse(value, "effective_config")
                .map_err(|_| effective_error(effective_path, "individual_id técnico inválido"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let identities =
        derive_technical_identity_plan(*root_seed, reference.id(), event_id.clone(), &individuals)
            .map_err(|error| {
                physical_error(effective_path, "technical_event_identities", &error)
            })?;
    let general: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let event_link = general
        .get("events")
        .and_then(Value::as_array)
        .and_then(|events| events.first())
        .and_then(Value::as_object)
        .ok_or_else(|| effective_error(effective_path, "falta enlace físico de technical_event"))?;
    let manifest_id = required_non_empty_text(event_link, "manifest_artifact_id", manifest_path)?;
    let artifact = artifacts
        .iter()
        .find(|artifact| artifact.id == manifest_id)
        .ok_or_else(|| effective_error(effective_path, "manifiesto técnico no catalogado"))?;
    let event_path = root.join(&artifact.path);
    let event: Value = serde_json::from_slice(
        &fs::read(&event_path)
            .map_err(|error| io_error(&event_path, "read_event_manifest", error))?,
    )
    .map_err(|error| {
        effective_error(
            effective_path,
            format!("manifiesto técnico ilegible: {error}"),
        )
    })?;
    let configured_libraries = effective_array(technical_event, "libraries", effective_path)?;
    let selected = identities
        .libraries()
        .map(|library| {
            let configured = configured_libraries
                .iter()
                .filter_map(Value::as_object)
                .find(|configured| {
                    configured.get("individual_id").and_then(Value::as_str)
                        == Some(library.individual_id().as_str())
                })
                .ok_or_else(|| {
                    effective_error(
                        effective_path,
                        "falta la UDI de una biblioteca técnica en la configuración efectiva",
                    )
                })?;
            Ok(json!({
                "individual_id": library.individual_id().as_str(),
                "library_id": library.library_id().as_str(),
                "read_group_id": library.read_group_id().as_str(),
                "i7": effective_non_empty_text(configured, "i7", effective_path)?,
                "i5": effective_non_empty_text(configured, "i5", effective_path)?,
            }))
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;
    let expected_event_subseed = identities.event_subseed().to_hex();
    let profile_matches = event
        .get("profile")
        .and_then(Value::as_object)
        .is_some_and(|profile| {
            technical_event
                .get("profile")
                .and_then(Value::as_object)
                .is_some_and(|technical_profile| {
                    profile.get("id") == technical_profile.get("id")
                        && profile.get("read_length") == technical_event.get("read_length")
                        && profile.get("quality") == technical_event.get("quality")
                })
        });
    let effective_chimera_rate = technical_event
        .get("artifacts")
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("chimera"))
        .and_then(Value::as_object)
        .and_then(|chimera| chimera.get("rate"));
    let effective_index_hopping_rate = technical_event
        .get("artifacts")
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("index_hopping"))
        .and_then(Value::as_object)
        .and_then(|index_hopping| index_hopping.get("rate"));
    let effective_sequencing_errors = technical_event.get("sequencing_errors");
    let event_sequencing_errors = event.get("error_quality").and_then(Value::as_object);
    let sequencing_errors_match = event_sequencing_errors.is_some_and(|error_quality| {
        effective_sequencing_errors.is_some_and(|effective_errors| {
            error_quality.get("requested_mode") == effective_errors.get("requested_mode")
                && error_quality.get("effective_mode") == effective_errors.get("effective_mode")
                && error_quality.get("qualities_represent_realized_error_probability")
                    == effective_errors.get("qualities_represent_realized_error_probability")
        })
    });
    let event_chimera = event
        .get("artifacts")
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("chimera"));
    let event_index_hopping_rate = event
        .get("demultiplexing")
        .and_then(Value::as_object)
        .and_then(|demultiplexing| demultiplexing.get("index_hopping_rate"));
    let general_chimera = general
        .get("technical_events")
        .and_then(Value::as_array)
        .and_then(|events| events.first())
        .and_then(Value::as_object)
        .and_then(|technical| technical.get("artifacts"))
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("chimera"));
    validate_f4_coverage_against_effective(&event, technical_event, reference, effective_path)?;
    validate_physical_biological_contamination_summary(
        &event,
        technical_event,
        effective,
        root_seed,
        &event_id,
        effective_path,
    )?;
    if event.get("event_id").and_then(Value::as_str) != Some(event_id.as_str())
        || event.get("event_subseed").and_then(Value::as_str)
            != Some(expected_event_subseed.as_str())
        || event.get("technical_batch_id").and_then(Value::as_str)
            != Some(identities.batch_id().as_str())
        || event
            .get("reference")
            .and_then(Value::as_object)
            .and_then(|reference| reference.get("reference_id"))
            .and_then(Value::as_str)
            != Some(reference.id().as_str())
        || event.get("selected_individuals").and_then(Value::as_array) != Some(&selected)
        || !profile_matches
        || event_chimera
            .and_then(Value::as_object)
            .and_then(|chimera| chimera.get("rate"))
            != effective_chimera_rate
        || event_index_hopping_rate != effective_index_hopping_rate
        || event_chimera != general_chimera
        || !sequencing_errors_match
    {
        return Err(effective_error(
            effective_path,
            "el manifiesto técnico físico no coincide con la configuración efectiva, semilla o referencia",
        ));
    }
    Ok(())
}

/// Contrasta el resumen agregado F4.4 con la configuración efectiva sin
/// retener ni reconstruir parejas. La materialización solo entrega los hechos
/// auditables del plan; aquí no se reintroduce procedencia individual.
fn validate_physical_biological_contamination_summary(
    event: &Value,
    technical_event: &serde_json::Map<String, Value>,
    effective: &serde_json::Map<String, Value>,
    root_seed: &RootSeed,
    event_id: &EventId,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let contamination =
        effective_object(technical_event, "biological_contamination", effective_path)?;
    let range = effective_object(contamination, "percent_range", effective_path)?;
    let minimum = range
        .get("minimum")
        .and_then(Value::as_f64)
        .ok_or_else(|| {
            effective_error(effective_path, "mínimo efectivo de contaminación inválido")
        })?;
    let maximum = range
        .get("maximum")
        .and_then(Value::as_f64)
        .ok_or_else(|| {
            effective_error(effective_path, "máximo efectivo de contaminación inválido")
        })?;
    let configured_recipients = effective_array(contamination, "recipient_ids", effective_path)?
        .iter()
        .map(|value| {
            value.as_str().map(str::to_owned).ok_or_else(|| {
                effective_error(
                    effective_path,
                    "receptora efectiva de contaminación inválida",
                )
            })
        })
        .collect::<Result<BTreeSet<_>, _>>()?;
    let selected_libraries = effective_array(technical_event, "libraries", effective_path)?
        .iter()
        .map(|value| {
            let library = value.as_object().ok_or_else(|| {
                effective_error(effective_path, "biblioteca técnica efectiva inválida")
            })?;
            Ok((
                effective_non_empty_text(library, "individual_id", effective_path)?.to_owned(),
                effective_non_empty_text(library, "library_id", effective_path)?.to_owned(),
            ))
        })
        .collect::<Result<BTreeMap<_, _>, Diagnostic>>()?;
    let populations = effective_array(effective, "initial_populations", effective_path)?
        .iter()
        .map(|value| {
            let population = value.as_object().ok_or_else(|| {
                effective_error(effective_path, "población inicial efectiva inválida")
            })?;
            let population_id =
                effective_non_empty_text(population, "population_id", effective_path)?;
            let individuals = effective_array(population, "individual_ids", effective_path)?;
            individuals
                .iter()
                .map(|value| {
                    value
                        .as_str()
                        .map(|individual| (individual.to_owned(), population_id.to_owned()))
                        .ok_or_else(|| {
                            effective_error(effective_path, "individual inicial efectivo inválido")
                        })
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?
        .into_iter()
        .flatten()
        .collect::<BTreeMap<_, _>>();
    let records = event
        .get("artifacts")
        .and_then(Value::as_object)
        .and_then(|artifacts| artifacts.get("biological_contamination"))
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(effective_path, "falta resumen físico de contaminación"))?;
    // Un resumen neutro puede omitirse. Cualquier intervalo positivo exige el
    // plan entero antes de que la publicación PE (F4.5) pueda abrir una ruta.
    if records.is_empty() {
        if maximum == 0.0 {
            return Ok(());
        }
        return Err(effective_error(
            effective_path,
            "un rango positivo de contaminación exige resumen por receptora",
        ));
    }
    let mut observed_recipients = BTreeSet::new();
    for value in records {
        let record = value.as_object().ok_or_else(|| {
            effective_error(effective_path, "resumen físico de contaminación inválido")
        })?;
        let recipient = effective_non_empty_text(record, "recipient_id", effective_path)?;
        let recipient_library =
            effective_non_empty_text(record, "recipient_library_id", effective_path)?;
        if !configured_recipients.contains(recipient)
            || selected_libraries.get(recipient).map(String::as_str) != Some(recipient_library)
            || !observed_recipients.insert(recipient.to_owned())
        {
            return Err(effective_error(
                effective_path,
                "receptora o biblioteca del resumen de contaminación no coincide con la configuración efectiva",
            ));
        }
        let record_range = effective_object(record, "percent_range", effective_path)?;
        if record_range.get("minimum").and_then(Value::as_f64) != Some(minimum)
            || record_range.get("maximum").and_then(Value::as_f64) != Some(maximum)
        {
            return Err(effective_error(
                effective_path,
                "rango del resumen de contaminación no coincide con la configuración efectiva",
            ));
        }
        let sampled = record
            .get("sampled_percent")
            .and_then(Value::as_f64)
            .ok_or_else(|| effective_error(effective_path, "porcentaje contaminante inválido"))?;
        let clean = record
            .get("clean_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento limpio contaminante inválido")
            })?;
        let target = record
            .get("target_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| effective_error(effective_path, "objetivo contaminante inválido"))?;
        let removed = record
            .get("removed_own_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento retirado contaminante inválido")
            })?;
        let inserted = record
            .get("inserted_contaminant_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento insertado contaminante inválido")
            })?;
        let retained = record
            .get("retained_own_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento retenido contaminante inválido")
            })?;
        let final_pairs = record
            .get("final_pair_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| {
                effective_error(effective_path, "recuento final contaminante inválido")
            })?;
        let recipient_id = IndividualId::parse(recipient, "effective_config")
            .map_err(|_| effective_error(effective_path, "receptora contaminante inválida"))?;
        let library_id = LibraryId::parse(recipient_library, "effective_config").map_err(|_| {
            effective_error(effective_path, "biblioteca receptora contaminante inválida")
        })?;
        let percent_subseed = derive_f4_technical_subseed(
            *root_seed,
            F4TechnicalDecision::contamination_percent(event_id, &library_id, &recipient_id),
        )
        .map_err(|_| {
            effective_error(
                effective_path,
                "no se pudo rederivar una subsemilla contaminante",
            )
        })?;
        let expected_subseeds = [
            ("percent", Ok(percent_subseed)),
            (
                "donors",
                derive_f4_technical_subseed(
                    *root_seed,
                    F4TechnicalDecision::contamination_donors(event_id, &library_id, &recipient_id),
                ),
            ),
            (
                "removal",
                derive_f4_technical_subseed(
                    *root_seed,
                    F4TechnicalDecision::contamination_recipients(
                        event_id,
                        &library_id,
                        &recipient_id,
                    ),
                ),
            ),
        ];
        let subseeds = effective_object(record, "subseeds", effective_path)?;
        for (name, expected) in expected_subseeds {
            let expected = expected.map_err(|_| {
                effective_error(
                    effective_path,
                    "no se pudo rederivar una subsemilla contaminante",
                )
            })?;
            if subseeds.get(name).and_then(Value::as_str) != Some(expected.to_hex().as_str()) {
                return Err(effective_error(
                    effective_path,
                    "subsemilla de contaminación no coincide con el catálogo F4",
                ));
            }
        }
        let expected_sampled =
            f4_uniform_closed(&mut percent_subseed.local_rng(), minimum, maximum).map_err(
                |_| {
                    effective_error(
                        effective_path,
                        "no se pudo remuestrear el porcentaje contaminante",
                    )
                },
            )?;
        if sampled.to_bits() != expected_sampled.to_bits() {
            return Err(effective_error(
                effective_path,
                "porcentaje contaminante no coincide con la subsemilla F4",
            ));
        }
        if contamination_target_pair_count(clean, sampled).ok() != Some(target) {
            return Err(effective_error(
                effective_path,
                "objetivo contaminante no coincide con el redondeo normativo",
            ));
        }
        if removed != target
            || inserted != target
            || retained.checked_add(removed) != Some(clean)
            || retained.checked_add(inserted) != Some(final_pairs)
            || final_pairs != clean
        {
            return Err(effective_error(
                effective_path,
                "los contadores físicos de contaminación no conservan las parejas",
            ));
        }
        let recipient_population = populations.get(recipient).ok_or_else(|| {
            effective_error(
                effective_path,
                "receptora contaminante fuera de población seleccionada",
            )
        })?;
        let contributions = effective_array(record, "donor_contributions", effective_path)?;
        let contribution_total = contributions.iter().try_fold(0_u64, |total, value| {
            let contribution = value
                .as_object()
                .ok_or_else(|| effective_error(effective_path, "aporte contaminante inválido"))?;
            let donor = effective_non_empty_text(contribution, "individual_id", effective_path)?;
            let pairs = contribution
                .get("pairs")
                .and_then(Value::as_u64)
                .ok_or_else(|| {
                    effective_error(effective_path, "pares de aporte contaminante inválidos")
                })?;
            if donor == recipient
                || !selected_libraries.contains_key(donor)
                || populations.get(donor) != Some(recipient_population)
            {
                return Err(effective_error(
                    effective_path,
                    "donante contaminante fuera de la población seleccionada",
                ));
            }
            total.checked_add(pairs).ok_or_else(|| {
                effective_error(effective_path, "desbordamiento de aportes contaminantes")
            })
        })?;
        if contribution_total != target {
            return Err(effective_error(
                effective_path,
                "los aportes físicos de contaminación no suman el objetivo",
            ));
        }
    }
    if observed_recipients != configured_recipients {
        return Err(effective_error(
            effective_path,
            "el resumen de contaminación no cubre exactamente las receptoras efectivas",
        ));
    }
    Ok(())
}

/// Comprueba los hechos que el manifiesto F4 puede contrastar sin convertir la
/// observación FASTQ final en una nueva configuración u objetivo.
fn validate_f4_coverage_against_effective(
    event: &Value,
    technical_event: &serde_json::Map<String, Value>,
    reference: &Reference,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let target = effective_object(technical_event, "target", effective_path)?;
    let target_x = target
        .get("coverage_wgs")
        .and_then(Value::as_f64)
        .filter(|value| value.is_finite() && *value >= 0.0)
        .ok_or_else(|| effective_error(effective_path, "objetivo WGS técnico inválido"))?;
    let expected_denominator = reference.total_length();
    let coverage = event
        .get("coverage")
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(effective_path, "falta cobertura F4 de evento"))?;
    for entry in coverage {
        let entry = entry
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "cobertura F4 de evento inválida"))?;
        if entry.get("target_x").and_then(Value::as_f64) != Some(target_x)
            || entry.get("denominator_bases").and_then(Value::as_u64) != Some(expected_denominator)
        {
            return Err(effective_error(
                effective_path,
                "la cobertura F4 no coincide con el objetivo o denominador estructural efectivo",
            ));
        }
    }
    Ok(())
}

#[allow(dead_code)]
fn validate_published_resources_v1(
    root: &Path,
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    validate_effective_config_against_manifest_bytes(
        effective_path,
        effective_bytes,
        manifest_path,
        manifest_bytes,
    )?;
    let artifacts = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let artifact_paths = artifacts
        .iter()
        .map(|artifact| (artifact.id.as_str(), artifact.path.as_str()))
        .collect::<BTreeMap<_, _>>();
    let effective: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let effective = effective
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "la raíz no es objeto"))?;
    validate_effective_targets_against_configuration_v1(
        root,
        effective,
        &artifact_paths,
        manifest_path,
        manifest_bytes,
        effective_path,
    )?;
    let root_seed = RootSeed::parse(
        "effective-config.seed.root",
        effective["seed"]["root"]
            .as_str()
            .ok_or_else(|| effective_error(effective_path, "seed.root inválida"))?,
    )
    .map_err(|error| physical_error(effective_path, "root_seed", &error))?;
    let reference_value = effective_object(effective, "reference", effective_path)?;
    let source_reference = read_physical_reference(
        root,
        reference_value,
        &artifact_paths,
        effective_path,
        if effective_needs_reference_bases(effective) {
            CatalogBases::Materialized
        } else {
            CatalogBases::SummaryOnly
        },
    )?;
    if reference_value.get("source").and_then(Value::as_str) == Some("simulated") {
        // Las versiones de este lector no conocen los microsatélites: no hay
        // tracto que colocar.
        validate_generated_reference(
            reference_value,
            root_seed,
            &source_reference,
            effective_path,
            None,
        )?;
    }
    let mut references =
        BTreeMap::from([(source_reference.id().as_str().to_owned(), source_reference)]);

    for reduction_value in effective_array(effective, "reduced_references", effective_path)? {
        let reduction = reduction_value
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        let id = effective_non_empty_text(reduction, "id", effective_path)?;
        let source_id = effective_non_empty_text(reduction, "source_reference_id", effective_path)?;
        let source = references.get(source_id).ok_or_else(|| {
            effective_error(
                effective_path,
                format!("la referencia fuente física {source_id} no existe"),
            )
        })?;
        let derived_subseed = derive_subseed(
            root_seed,
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: source_id,
                derived_reference_id: id,
            },
        )
        .map_err(|error| physical_error(effective_path, "derive_reduction_subseed", &error))?;
        if effective_non_empty_text(reduction, "selection_subseed", effective_path)?
            != derived_subseed.to_hex()
        {
            return Err(effective_error(
                effective_path,
                format!("selection_subseed de {id} no reproduce el ámbito canónico"),
            ));
        }
        let eligible = source
            .sequences()
            .iter()
            .map(|sequence| {
                EligibleSequence::new(sequence.id().clone(), sequence.class(), sequence.length())
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| {
                effective_error(
                    effective_path,
                    format!("catálogo elegible inválido: {error:?}"),
                )
            })?;
        let target_bases = effective_positive_u64(
            effective_object(reduction, "target", effective_path)?,
            "bases",
            effective_path,
        )?;
        let region_size_bases =
            effective_positive_u64(reduction, "region_size_bases", effective_path)?;
        let mut rng = derived_subseed.local_rng();
        let selection = select_reduced_reference(
            source.id(),
            &eligible,
            target_bases,
            region_size_bases,
            &mut rng,
        )
        .map_err(|error| {
            effective_error(
                effective_path,
                format!("selección física no reproducible: {error:?}"),
            )
        })?;
        let derived_id = ReferenceId::parse(id, "published-resource-validation")
            .map_err(|error| physical_error(effective_path, "derived_reference_id", &error))?;
        let plan = materialize_selection(&selection, &derived_id, &eligible).map_err(|error| {
            effective_error(
                effective_path,
                format!("plan físico no materializable: {error:?}"),
            )
        })?;
        let declared_order = effective_array(reduction, "selection_order", effective_path)?
            .iter()
            .map(Value::as_u64)
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| effective_error(effective_path, "selection_order inválido"))?;
        let expected_order = plan
            .segments()
            .iter()
            .map(|segment| segment.selection_order())
            .collect::<Vec<_>>();
        if declared_order != expected_order {
            return Err(effective_error(
                effective_path,
                format!("selection_order de {id} no coincide con la subsemilla"),
            ));
        }
        let expected_derived = build_reduced_reference(
            source_for(source, None)
                .map_err(|error| physical_error(effective_path, "reduction_source", &error))?
                .as_ref(),
            &plan,
        )
        .map_err(|error| physical_error(effective_path, "build_reduced_reference", &error))?;
        let fasta_id = effective_non_empty_text(reduction, "fasta_artifact_id", effective_path)?;
        let derived = read_reference_with_expected(
            root,
            id,
            fasta_id,
            effective_array(reduction, "sequences", effective_path)?,
            expected_derived.sequences(),
            &artifact_paths,
            effective_path,
        )?;
        compare_reference_bases(&expected_derived, &derived, effective_path, id)?;
        let map_id = effective_non_empty_text(reduction, "map_artifact_id", effective_path)?;
        let map_relative = artifact_paths.get(map_id).ok_or_else(|| {
            effective_error(effective_path, format!("mapa físico ausente: {map_id}"))
        })?;
        let map_path = root.join(map_relative);
        let map_bytes = std::fs::read(&map_path)
            .map_err(|error| io_error(&map_path, "read_reduction_map", error))?;
        let map: Value = serde_json::from_slice(&map_bytes).map_err(|error| {
            effective_error(
                effective_path,
                format!("mapa {map_id} no es JSON legible: {error}"),
            )
        })?;
        if map != reduction_map_document(&plan) {
            return Err(effective_error(
                effective_path,
                format!("mapa físico {map_id} no coincide con subsemilla y selección"),
            ));
        }
        validate_reduction_map_document(
            source_for(source, None)
                .map_err(|error| physical_error(effective_path, "reduction_source", &error))?
                .as_ref(),
            &derived,
            &map,
            ReductionMapParameters {
                region_size_bases,
                target_bases,
            },
        )
        .map_err(|error| physical_error(effective_path, "validate_reduction_map", &error))?;
        references.insert(id.to_owned(), derived);
    }
    Ok(())
}

#[allow(dead_code)]
fn validate_effective_targets_against_configuration_v1(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    artifact_paths: &BTreeMap<&str, &str>,
    manifest_path: &Path,
    manifest_bytes: &[u8],
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let configuration_id = manifest["configuration_input"]["artifact_id"]
        .as_str()
        .ok_or_else(|| effective_error(effective_path, "configuration_input no enlaza TOML"))?;
    let relative = artifact_paths.get(configuration_id).ok_or_else(|| {
        effective_error(
            effective_path,
            format!("TOML físico ausente: {configuration_id}"),
        )
    })?;
    let configuration_path = root.join(relative);
    let bytes = std::fs::read(&configuration_path)
        .map_err(|error| io_error(&configuration_path, "read_configuration", error))?;
    let digest = configuration_digest(&bytes)
        .map_err(|error| physical_error(effective_path, "configuration_digest", &error))?;
    let original_configuration_path = manifest["configuration_input"]["original_path"]
        .as_str()
        .ok_or_else(|| effective_error(effective_path, "configuration_input sin procedencia"))?;
    let admitted = admit_toml_header(
        Path::new(original_configuration_path),
        bytes,
        ExpectedConfigKind::General,
    )
    .map_err(|error| physical_error(effective_path, "admit_configuration", &error))?;
    let declared = parse_declared_config(admitted)
        .map_err(|error| physical_error(effective_path, "parse_configuration", &error))?;
    let seed = effective_object(effective, "seed", effective_path)?;
    let root_seed = RootSeed::parse(
        "effective-config.seed.root",
        effective_non_empty_text(seed, "root", effective_path)?,
    )
    .map_err(|error| physical_error(effective_path, "root_seed", &error))?;
    let mut entropy = ReplayEntropy(root_seed);
    let expected = build_effective_config_with_entropy(&declared, &mut entropy)
        .map_err(|error| physical_error(effective_path, "rebuild_effective_config", &error))?;
    let expected_execution = derive_execution_id(root_seed, &digest)
        .map_err(|error| physical_error(effective_path, "derive_execution_id", &error))?;
    if manifest.get("execution_id").and_then(Value::as_str) != Some(expected_execution.as_str()) {
        return Err(effective_error(
            effective_path,
            "execution_id no coincide con la semilla y los bytes TOML",
        ));
    }

    if seed.get("root").and_then(Value::as_str)
        != Some(expected.reproducibility().root_seed().to_hex().as_str())
        || seed.get("origin").and_then(Value::as_str)
            != Some(expected.reproducibility().seed_origin().as_str())
    {
        return Err(effective_error(
            effective_path,
            "la semilla efectiva no coincide con la declaración TOML",
        ));
    }

    let reference = effective_object(effective, "reference", effective_path)?;
    let expected_source = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => "simulated",
        EffectiveReferenceSource::Supplied { .. } => "supplied",
    };
    if effective_non_empty_text(reference, "id", effective_path)? != expected.reference_id()
        || reference.get("source").and_then(Value::as_str) != Some(expected_source)
        || reference.get("total_genome_size").and_then(Value::as_u64)
            != Some(expected.total_genome_size())
    {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no coincide con la declaración TOML",
        ));
    }
    let expected_original_path = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => Value::Null,
        EffectiveReferenceSource::Supplied { original_path } => {
            Value::String(original_path.to_string_lossy().into_owned())
        }
    };
    if reference.get("original_fasta_path") != Some(&expected_original_path) {
        return Err(effective_error(
            effective_path,
            "la procedencia FASTA efectiva no coincide con la declaración TOML",
        ));
    }
    let expected_canonical_format = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => canonical_fasta_json(),
        EffectiveReferenceSource::Supplied { .. } => Value::Null,
    };
    if reference.get("canonical_format") != Some(&expected_canonical_format) {
        return Err(effective_error(
            effective_path,
            "el formato FASTA efectivo no corresponde al origen declarado",
        ));
    }

    let expected_composition = expected.composition().map(|composition| {
        let values = composition.values();
        json!({"A": values[0], "C": values[1], "G": values[2], "T": values[3]})
    });
    if reference.get("composition") != Some(&expected_composition.unwrap_or(Value::Null)) {
        return Err(effective_error(
            effective_path,
            "la composición efectiva no coincide con la declaración TOML",
        ));
    }

    let effective_sequences = effective_array(reference, "sequences", effective_path)?;
    if effective_sequences.len() != expected.sequences().len() {
        return Err(effective_error(
            effective_path,
            "el número de secuencias efectivo no coincide con el TOML",
        ));
    }
    for (observed, expected) in effective_sequences.iter().zip(expected.sequences()) {
        let observed = observed
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        if effective_non_empty_text(observed, "id", effective_path)? != expected.id()
            || observed.get("class").and_then(Value::as_str)
                != Some(sequence_class(expected.class()))
            || observed.get("length").and_then(Value::as_u64) != Some(expected.length())
        {
            return Err(effective_error(
                effective_path,
                "orden, ID, clase o longitud de secuencia no coincide con el TOML",
            ));
        }
    }

    let effective_reductions = effective_array(effective, "reduced_references", effective_path)?;
    if effective_reductions.len() != expected.reduced_references().len() {
        return Err(effective_error(
            effective_path,
            "las reducciones efectivas no coinciden con las declaradas en TOML",
        ));
    }
    for (expected, effective) in expected
        .reduced_references()
        .iter()
        .zip(effective_reductions)
    {
        let effective = effective
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        let target = effective_object(effective, "target", effective_path)?;
        let expected_declared = expected.target().declared_value();
        if effective_non_empty_text(effective, "id", effective_path)? != expected.id()
            || effective_non_empty_text(effective, "source_reference_id", effective_path)?
                != expected.source_reference_id()
            || effective.get("region_size_bases").and_then(Value::as_u64)
                != Some(expected.region_size_bases())
            || target.get("mode").and_then(Value::as_str) != Some(expected.target().mode())
            || target.get("declared").and_then(Value::as_str) != Some(expected_declared.as_str())
            || target.get("bases").and_then(Value::as_u64) != Some(expected.target_bases())
            || effective.get("selection_subseed").and_then(Value::as_str)
                != Some(expected.selection_subseed().to_hex().as_str())
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "la reducción efectiva {} no coincide completamente con el TOML",
                    expected.id()
                ),
            ));
        }
    }
    Ok(())
}

/// El objeto de la configuración efectiva que describe la referencia **de la
/// simulación**: `inserted_reference` si la ejecución inserta tractos
/// (DTC-098) y `reference` en los demás casos.
fn simulation_reference_object<'a>(
    effective: &'a serde_json::Map<String, Value>,
    effective_path: &Path,
) -> Result<&'a serde_json::Map<String, Value>, Diagnostic> {
    match effective
        .get("inserted_reference")
        .and_then(Value::as_object)
    {
        Some(inserted) => Ok(inserted),
        None => effective_object(effective, "reference", effective_path),
    }
}

#[allow(clippy::too_many_arguments)]
fn validate_physical_population_snvs(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    physical_config: &EffectiveConfig,
    artifacts: &[ArtifactRecord],
    reference: &Reference,
    effective_path: &Path,
    allow_empty_snv_matrix: bool,
    reused_snv_metadata: Option<&Path>,
) -> Result<Option<PreparedPopulationSnvs>, Diagnostic> {
    let Some(link) = effective.get("population_snvs").and_then(Value::as_object) else {
        return Ok(None);
    };
    if link.get("mode").and_then(Value::as_str) == Some("reuse") {
        let crate::application::EffectivePopulationSnvs::Reuse { .. } =
            physical_config.population_snvs().ok_or_else(|| {
                effective_error(effective_path, "reutilización SNV sin fuente física")
            })?
        else {
            return Err(effective_error(
                effective_path,
                "modo SNV efectivo no coincide con reutilización",
            ));
        };
        // La matriz reutilizada no viaja en el paquete: se lee de la ubicación
        // que indica quien juzga, nunca de la ruta del TOML copiado, que es
        // relativa al directorio de trabajo de la ejecución original.
        //
        // Falta de un argumento del llamador, no un defecto del documento: el
        // catálogo no tiene un código de «argumento ausente», y el más cercano
        // es el de entrada declarada que no se puede leer.
        let metadata_path = reused_snv_metadata.ok_or_else(|| {
            Diagnostic::input(
                "GPS-INPUT-001",
                "No se puede abrir o leer una entrada declarada.",
            )
            .with("path", effective_path.display().to_string())
            .with(
                "detail",
                "el paquete reutiliza la matriz SNV (mode = \"reuse\") y esta no viaja en él: \
                 falta indicar la ruta de population_snvs.metadata.json",
            )
        })?;
        let reference_id = simulation_reference_object(effective, effective_path)?
            .get("fasta_artifact_id")
            .and_then(Value::as_str)
            .ok_or_else(|| effective_error(effective_path, "referencia sin artefacto FASTA"))?;
        let reference_artifact = artifacts
            .iter()
            .find(|artifact| artifact.id == reference_id)
            .ok_or_else(|| effective_error(effective_path, "FASTA SNV ausente del catálogo"))?;
        let (prepared, _, metadata_sha256, normalization_audit) =
            validate_reusable_population_snvs(
                metadata_path,
                reference,
                &reference_artifact.digest,
                physical_config.normalization().tolerance(),
                allow_empty_snv_matrix,
            )
            .map_err(|error| {
                effective_error(
                    effective_path,
                    format!("reuso SNV físico inválido: {error}"),
                )
            })?;
        if link.get("metadata_sha256").and_then(Value::as_str) != Some(metadata_sha256.as_str())
            || link.get("normalization_audit") != Some(&normalization_audit)
        {
            return Err(effective_error(
                effective_path,
                "enlace SNV reutilizado no coincide con el recurso físico",
            ));
        }
        return Ok(Some(prepared));
    }
    let metadata_id = link
        .get("metadata_artifact_id")
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(effective_path, "population_snvs publicado sin metadata"))?;
    let metadata = artifacts
        .iter()
        .find(|artifact| artifact.id == metadata_id)
        .ok_or_else(|| effective_error(effective_path, "metadata SNV ausente del catálogo"))?;
    let matrix = artifacts
        .iter()
        .find(|artifact| artifact.id == "population-snvs-tsv")
        .ok_or_else(|| effective_error(effective_path, "TSV SNV ausente del catálogo"))?;
    let reference_id = simulation_reference_object(effective, effective_path)?
        .get("fasta_artifact_id")
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(effective_path, "referencia sin artefacto FASTA"))?;
    let reference_artifact = artifacts
        .iter()
        .find(|artifact| artifact.id == reference_id)
        .ok_or_else(|| effective_error(effective_path, "FASTA SNV ausente del catálogo"))?;
    let metadata_path = root.join(&metadata.path);
    let matrix_path = root.join(&matrix.path);
    let metadata_bytes = read_verified_artifact_bytes(&metadata_path, metadata, effective_path)?;
    let matrix_bytes = read_verified_artifact_bytes(&matrix_path, matrix, effective_path)?;
    let (prepared, provenance, metadata_sha256, normalization_audit) =
        validate_reusable_population_snvs_bytes(
            &metadata_bytes,
            &matrix_bytes,
            reference,
            &reference_artifact.digest,
            effective
                .get("normalization")
                .and_then(|value| value.get("tolerance"))
                .and_then(Value::as_f64)
                .ok_or_else(|| {
                    effective_error(effective_path, "normalization.tolerance inválida")
                })?,
            allow_empty_snv_matrix,
        )
        .map_err(|error| {
            effective_error(effective_path, format!("par SNV físico inválido: {error}"))
        })?;
    if link.get("metadata_sha256").and_then(Value::as_str) != Some(metadata_sha256.as_str())
        || link.get("mode").and_then(Value::as_str) != Some(provenance.as_str())
        || link.get("normalization_audit") != Some(&normalization_audit)
    {
        return Err(effective_error(
            effective_path,
            "enlace population_snvs no coincide con el par físico",
        ));
    }
    Ok(Some(prepared))
}

fn validate_physical_initial_population_statistics(
    effective: &serde_json::Map<String, Value>,
    physical_config: &EffectiveConfig,
    reference: &Reference,
    matrix: Option<&PreparedPopulationSnvs>,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let declared = effective
        .get("population_snvs")
        .and_then(Value::as_object)
        .and_then(|link| link.get("initial_population_statistics"));
    if physical_config.initial_populations().is_empty() {
        if declared.is_some() {
            return Err(effective_error(
                effective_path,
                "hay estadísticas iniciales sin poblaciones iniciales",
            ));
        }
        return Ok(());
    }
    let matrix = matrix.ok_or_else(|| {
        effective_error(
            effective_path,
            "faltan matriz y estadísticas para poblaciones iniciales",
        )
    })?;
    let authorizations = effective
        .get("authorizations")
        .and_then(Value::as_object)
        .ok_or_else(|| effective_error(effective_path, "authorizations inválido"))?;
    let authorization = if authorizations
        .get("allow_incompatible_genotypes")
        .and_then(Value::as_bool)
        .ok_or_else(|| effective_error(effective_path, "allow_incompatible_genotypes inválido"))?
    {
        IncompatibleGenotypeAuthorization::allowed()
    } else {
        IncompatibleGenotypeAuthorization::denied()
    };
    let rebuilt = materialize_initial_individuals(
        physical_config,
        reference,
        matrix,
        authorization,
        |individual_id| {
            physical_config
                .direct_events()
                .iter()
                .filter(|event| event.individual_id() == individual_id)
                .map(|event| {
                    declared_small_event(event, reference).map_err(|error| {
                        InitialIndividualsMaterializationError::Sampling(error.render_human())
                    })
                })
                .collect()
        },
    )
    .map_err(|error| {
        effective_error(
            effective_path,
            format!("estadísticas iniciales no reconstruibles: {error}"),
        )
    })?;
    if declared != Some(rebuilt.statistics()) {
        return Err(effective_error(
            effective_path,
            "estadísticas iniciales no coinciden con matriz, semilla y eventos físicos",
        ));
    }
    Ok(())
}

/// **El juez completo de la verdad sobre los bytes publicados** (DT-020 y
/// DT-021), con la matriz publicada y la población de cada individuo según la
/// configuración efectiva: cada caso contra las copias que reconstruye por su
/// cuenta, y el conjunto de casos, todo-REF incluidos.
///
/// La publicación de la verdad ya juzgó el documento en memoria; esto juzga
/// lo que quedó escrito en el paquete, que es lo que leerá quien lo use. Hasta
/// el 3 de octubre de 2026 este juez solo miraba la cabecera de la verdad.
fn validate_physical_truth(
    root: &Path,
    artifacts: &[ArtifactRecord],
    physical_config: &EffectiveConfig,
    reference: &Reference,
    matrix: Option<&PreparedPopulationSnvs>,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    // Sin poblaciones iniciales la verdad no lleva copias ni casos, y la
    // referencia puede haberse leído sin bases: su cabecera la juzga ya
    // `validate_published_reference_names`.
    let Some(truth) = artifacts
        .iter()
        .find(|artifact| artifact.id == "truth-jsonl")
        .filter(|_| !physical_config.initial_populations().is_empty())
    else {
        return Ok(());
    };
    let populations = physical_config
        .initial_populations()
        .iter()
        .flat_map(|population| {
            population
                .individual_ids()
                .iter()
                .map(move |individual| (individual.as_str(), population.population_id().as_str()))
        });
    let case_loci = match matrix {
        Some(matrix) => crate::truth_validation::MatrixCaseLoci::from_matrix(matrix, populations),
        None => crate::truth_validation::MatrixCaseLoci::from_populations(populations),
    };
    let path = root.join(&truth.path);
    let document = std::fs::File::open(&path)
        .map_err(|error| format!("{error}"))
        .and_then(|file| {
            crate::truth_validation::read_truth_jsonl(&mut std::io::BufReader::new(file))
                .and_then(|document| {
                    crate::truth_validation::validate_published_truth(
                        &document, reference, &case_loci,
                    )
                    .map(|()| document)
                })
                .map_err(|error| format!("{error:?}"))
        })
        .map_err(|detail| {
            effective_error(
                effective_path,
                format!("la verdad publicada no supera el juez de la verdad: {detail}"),
            )
        })?;
    // **Qué copias hay.** El juez de la verdad recorre los individuos y las
    // copias que el documento trae; al proyectar, el contraste con lo simulado
    // exige además que estén todas, pero aquí no hay estado simulado. Las
    // copias tienen que ser exactamente cada individuo de la configuración
    // efectiva, por cada secuencia de la referencia, por cada copia de la
    // ploidía de su población.
    let expected = physical_config
        .initial_populations()
        .iter()
        .flat_map(|population| {
            population
                .individual_ids()
                .iter()
                .flat_map(move |individual| {
                    reference.sequences().iter().flat_map(move |sequence| {
                        (0..population.ploidy().copy_count()).map(move |copy_index| {
                            crate::truth_schema::CopyIdentity {
                                reference_id: reference.id().as_str().to_owned(),
                                individual_id: individual.as_str().to_owned(),
                                sequence_id: sequence.id().as_str().to_owned(),
                                copy_index,
                            }
                        })
                    })
                })
        })
        .collect::<std::collections::BTreeSet<_>>();
    let published = document
        .records
        .iter()
        .filter_map(|record| match record {
            crate::truth_schema::TruthRecord::HaplotypeCopy(copy) => Some(&copy.identity),
            crate::truth_schema::TruthRecord::Genotype(_) => None,
        })
        .collect::<std::collections::BTreeSet<_>>();
    if let Some(missing) = expected
        .iter()
        .find(|identity| !published.contains(identity))
    {
        return Err(effective_error(
            effective_path,
            format!("a la verdad publicada le falta una copia de la simulación: {missing:?}"),
        ));
    }
    if let Some(extra) = published
        .iter()
        .find(|identity| !expected.contains(**identity))
    {
        return Err(effective_error(
            effective_path,
            format!("la verdad publicada trae una copia que la simulación no tiene: {extra:?}"),
        ));
    }
    Ok(())
}

fn read_verified_artifact_bytes(
    path: &Path,
    artifact: &ArtifactRecord,
    effective_path: &Path,
) -> Result<Vec<u8>, Diagnostic> {
    let bytes = fs::read(path).map_err(|error| io_error(path, "read_population_snvs", error))?;
    let observed = identify_file(path, Path::new(&artifact.path))?;
    if observed.size_bytes() != artifact.size_bytes || observed.digest() != artifact.digest {
        return Err(effective_error(
            effective_path,
            format!(
                "artefacto SNV físico {} no coincide con su SHA-256",
                artifact.id
            ),
        ));
    }
    Ok(bytes)
}

fn validate_physical_vcf_output(
    root: &Path,
    effective: &EffectiveConfig,
    manifest_bytes: &[u8],
    artifacts: &[ArtifactRecord],
    reference: &Reference,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let manifest = manifest
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "la raíz no es un objeto"))?;
    if !is_supported_general_manifest_version(
        manifest.get("schema_version").and_then(Value::as_str),
    ) {
        return Err(effective_error(
            effective_path,
            "el manifiesto no usa la versión general vigente",
        ));
    }
    let outputs = required_array(manifest, "vcf_outputs", manifest_path)?;
    let expected = match effective.vcf_output() {
        None => Vec::new(),
        Some(VcfOutputMode::Text) => {
            vec![(VCF_ARTIFACT_ID, VCF_ARTIFACT_KIND, VCF_TEXT_ARTIFACT_PATH)]
        }
        Some(VcfOutputMode::Bgzf) => {
            let index = if reference
                .sequences()
                .iter()
                .all(|sequence| sequence.length() < TBI_MAXIMUM_POSITION_EXCLUSIVE)
            {
                ("vcf_tbi_index", VCF_TBI_ARTIFACT_PATH)
            } else {
                ("vcf_csi_index", VCF_CSI_ARTIFACT_PATH)
            };
            vec![
                (VCF_ARTIFACT_ID, VCF_ARTIFACT_KIND, VCF_BGZF_ARTIFACT_PATH),
                (VCF_INDEX_ARTIFACT_ID, index.0, index.1),
            ]
        }
    };
    if (effective.vcf_output().is_none() && !outputs.is_empty())
        || (effective.vcf_output().is_some() && outputs.len() != 1)
    {
        return Err(effective_error(
            effective_path,
            "vcf_outputs no coincide con la activación efectiva",
        ));
    }
    let actual = artifacts
        .iter()
        .filter(|artifact| {
            matches!(
                artifact.kind.as_str(),
                "truth_vcf" | "vcf_tbi_index" | "vcf_csi_index"
            )
        })
        .collect::<Vec<_>>();
    if actual.len() != expected.len() {
        return Err(effective_error(
            effective_path,
            "el catálogo físico VCF no contiene exactamente los artefactos esperados",
        ));
    }
    for (id, kind, relative_path) in expected {
        let artifact = artifacts
            .iter()
            .find(|artifact| artifact.id == id)
            .ok_or_else(|| {
                effective_error(
                    effective_path,
                    format!("falta el artefacto VCF esperado {id}"),
                )
            })?;
        if artifact.kind != kind || artifact.path != relative_path {
            return Err(effective_error(
                effective_path,
                format!("el artefacto {id} no tiene clase o ruta VCF canónica"),
            ));
        }
        let physical_path = root.join(relative_path);
        let identity = identify_file(&physical_path, Path::new(relative_path))?;
        if identity.size_bytes() != artifact.size_bytes || identity.digest() != artifact.digest {
            return Err(artifact_identity_error()
                .with("path", physical_path.display().to_string())
                .with("artifact_id", id));
        }
    }
    Ok(())
}

fn validate_effective_targets_against_configuration(
    root: &Path,
    effective: &serde_json::Map<String, Value>,
    artifact_paths: &BTreeMap<&str, &str>,
    manifest_path: &Path,
    manifest_bytes: &[u8],
    effective_path: &Path,
) -> Result<EffectiveConfig, Diagnostic> {
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let configuration_id = manifest["configuration_input"]["artifact_id"]
        .as_str()
        .ok_or_else(|| effective_error(effective_path, "configuration_input no enlaza TOML"))?;
    let relative = artifact_paths.get(configuration_id).ok_or_else(|| {
        effective_error(
            effective_path,
            format!("TOML físico ausente: {configuration_id}"),
        )
    })?;
    let configuration_path = root.join(relative);
    let bytes = std::fs::read(&configuration_path)
        .map_err(|error| io_error(&configuration_path, "read_configuration", error))?;
    let digest = configuration_digest(&bytes)
        .map_err(|error| physical_error(effective_path, "configuration_digest", &error))?;
    let original_configuration_path = manifest["configuration_input"]["original_path"]
        .as_str()
        .ok_or_else(|| effective_error(effective_path, "configuration_input sin procedencia"))?;
    let admitted = admit_toml_header(
        Path::new(original_configuration_path),
        bytes,
        ExpectedConfigKind::General,
    )
    .map_err(|error| physical_error(effective_path, "admit_configuration", &error))?;
    let declared = parse_declared_config(admitted)
        .map_err(|error| physical_error(effective_path, "parse_configuration", &error))?;
    let seed = effective_object(effective, "seed", effective_path)?;
    let root_seed = RootSeed::parse(
        "effective-config.seed.root",
        effective_non_empty_text(seed, "root", effective_path)?,
    )
    .map_err(|error| physical_error(effective_path, "root_seed", &error))?;
    let mut entropy = ReplayEntropy(root_seed);
    let expected = build_effective_config_with_entropy(&declared, &mut entropy)
        .map_err(|error| physical_error(effective_path, "rebuild_effective_config", &error))?;
    let expected_execution = derive_execution_id(root_seed, &digest)
        .map_err(|error| physical_error(effective_path, "derive_execution_id", &error))?;
    if manifest.get("execution_id").and_then(Value::as_str) != Some(expected_execution.as_str()) {
        return Err(effective_error(
            effective_path,
            "execution_id no coincide con la semilla y los bytes TOML",
        ));
    }

    if seed.get("root").and_then(Value::as_str)
        != Some(expected.reproducibility().root_seed().to_hex().as_str())
        || seed.get("origin").and_then(Value::as_str)
            != Some(expected.reproducibility().seed_origin().as_str())
    {
        return Err(effective_error(
            effective_path,
            "la semilla efectiva no coincide con la declaración TOML",
        ));
    }

    if effective.get("ploidy").and_then(Value::as_u64)
        != Some(u64::from(expected.ploidy().copy_count()))
    {
        return Err(effective_error(
            effective_path,
            "la ploidía efectiva no coincide con la declaración TOML",
        ));
    }
    if effective.get("initial_populations") != Some(&effective_initial_populations_json(&expected))
    {
        return Err(effective_error(
            effective_path,
            "las cohortes iniciales efectivas no coinciden exactamente con el TOML",
        ));
    }
    let normalization = effective_object(effective, "normalization", effective_path)?;
    if normalization
        .get("tolerance")
        .and_then(Value::as_f64)
        .map(f64::to_bits)
        != Some(expected.normalization().tolerance().to_bits())
    {
        return Err(effective_error(
            effective_path,
            "la tolerancia efectiva no coincide con la declaración TOML",
        ));
    }

    let reference = effective_object(effective, "reference", effective_path)?;
    let expected_source = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => "simulated",
        EffectiveReferenceSource::Supplied { .. } => "supplied",
    };
    if effective_non_empty_text(reference, "id", effective_path)? != expected.reference_id()
        || reference.get("source").and_then(Value::as_str) != Some(expected_source)
        || reference.get("total_genome_size").and_then(Value::as_u64)
            != Some(expected.total_genome_size())
    {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no coincide con la declaración TOML",
        ));
    }
    let expected_original_path = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => Value::Null,
        EffectiveReferenceSource::Supplied { original_path } => {
            Value::String(original_path.to_string_lossy().into_owned())
        }
    };
    if reference.get("original_fasta_path") != Some(&expected_original_path) {
        return Err(effective_error(
            effective_path,
            "la procedencia FASTA efectiva no coincide con la declaración TOML",
        ));
    }
    let expected_canonical_format = match expected.reference_source() {
        EffectiveReferenceSource::Simulated => canonical_fasta_json(),
        EffectiveReferenceSource::Supplied { .. } => Value::Null,
    };
    if reference.get("canonical_format") != Some(&expected_canonical_format) {
        return Err(effective_error(
            effective_path,
            "el formato FASTA efectivo no corresponde al origen declarado",
        ));
    }

    let expected_composition = expected.composition().map(|composition| {
        let values = composition.values();
        json!({"A": values[0], "C": values[1], "G": values[2], "T": values[3], "N": composition.unknown_base_fraction()})
    });
    if reference.get("composition") != Some(&expected_composition.unwrap_or(Value::Null)) {
        return Err(effective_error(
            effective_path,
            "la composición efectiva no coincide con la declaración TOML",
        ));
    }

    let effective_sequences = effective_array(reference, "sequences", effective_path)?;
    if effective_sequences.len() != expected.sequences().len() {
        return Err(effective_error(
            effective_path,
            "el número de secuencias efectivo no coincide con el TOML",
        ));
    }
    for (observed, expected) in effective_sequences.iter().zip(expected.sequences()) {
        let observed = observed
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        if effective_non_empty_text(observed, "id", effective_path)? != expected.id()
            || observed.get("class").and_then(Value::as_str)
                != Some(sequence_class(expected.class()))
            || observed.get("length").and_then(Value::as_u64) != Some(expected.length())
        {
            return Err(effective_error(
                effective_path,
                "orden, ID, clase o longitud de secuencia no coincide con el TOML",
            ));
        }
    }

    let effective_reductions = effective_array(effective, "reduced_references", effective_path)?;
    if effective_reductions.len() != expected.reduced_references().len() {
        return Err(effective_error(
            effective_path,
            "las reducciones efectivas no coinciden con las declaradas en TOML",
        ));
    }
    for (expected, effective) in expected
        .reduced_references()
        .iter()
        .zip(effective_reductions)
    {
        let effective = effective
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        let target = effective_object(effective, "target", effective_path)?;
        let expected_declared = expected.target().declared_value();
        if effective_non_empty_text(effective, "id", effective_path)? != expected.id()
            || effective_non_empty_text(effective, "source_reference_id", effective_path)?
                != expected.source_reference_id()
            || effective.get("region_size_bases").and_then(Value::as_u64)
                != Some(expected.region_size_bases())
            || target.get("mode").and_then(Value::as_str) != Some(expected.target().mode())
            || target.get("declared").and_then(Value::as_str) != Some(expected_declared.as_str())
            || target.get("bases").and_then(Value::as_u64) != Some(expected.target_bases())
            || effective.get("selection_subseed").and_then(Value::as_str)
                != Some(expected.selection_subseed().to_hex().as_str())
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "la reducción efectiva {} no coincide completamente con el TOML",
                    expected.id()
                ),
            ));
        }
    }
    // **Las regiones del generador de la matriz, traducidas aquí por cuenta
    // propia** (DTC-098). La resolución que se acaba de rehacer es la del
    // productor; si su traducción estuviera mal, compararla consigo misma no lo
    // vería. Por eso cada región del TOML se traduce con una suma escrita aquí,
    // desde los cambios que la configuración efectiva publica —que el juez de
    // la referencia final ata a los loci—, y se exige la misma región que la
    // resolución entregó al generador. Una región que un tracto partiría tiene
    // que haberse rechazado antes de publicar.
    if let (
        Some(changes),
        Some(crate::configuration::DeclaredPopulationSnvs::Generate { parameters }),
    ) = (
        effective
            .get("inserted_reference")
            .and_then(|inserted| inserted.get("insertions"))
            .and_then(Value::as_array),
        declared.population_snvs(),
    ) {
        let points = changes
            .iter()
            .filter_map(|change| {
                Some((
                    change.get("sequence_id")?.as_str()?,
                    change.get("insert_after")?.as_u64()?,
                    change.get("inserted_bases")?.as_u64()?,
                ))
            })
            .collect::<Vec<_>>();
        let translate = |contig: &str, position: u64| {
            position
                + points
                    .iter()
                    .filter(|(sequence, anchor, _)| *sequence == contig && *anchor < position)
                    .map(|(_, _, bases)| bases)
                    .sum::<u64>()
        };
        let translated = parameters
            .regions
            .iter()
            .map(|region| {
                let (contig, start, end) = (
                    region.contig.value().as_str(),
                    *region.start.value(),
                    *region.end.value(),
                );
                let split = points.iter().any(|(sequence, anchor, _)| {
                    *sequence == contig && start <= *anchor && *anchor < end
                });
                (!split).then(|| {
                    (
                        contig.to_owned(),
                        translate(contig, start),
                        translate(contig, end),
                    )
                })
            })
            .collect::<Option<Vec<_>>>();
        let resolved = match expected
            .population_snvs_generation()
            .map(|generation| &generation.spatial)
        {
            Some(crate::population_snvs::GeneratedSpatialMode::WeightedRegions {
                regions, ..
            }) => regions
                .iter()
                .map(|region| {
                    (
                        region.contig.clone(),
                        region.start_1based,
                        region.end_1based,
                    )
                })
                .collect::<Vec<_>>(),
            _ => Vec::new(),
        };
        if translated.as_ref() != Some(&resolved) {
            return Err(effective_error(
                effective_path,
                "las regiones del generador de la matriz no son la traducción al genoma final de las del TOML",
            ));
        }
    }
    // **La inserción publicada es la que el TOML declara** (DTC-098): la misma
    // referencia final, los mismos cambios y, en cada locus, lo mismo
    // declarado. La traducción de lo declarado al genoma final la juzga por su
    // cuenta `judge_inserted_reference`, sin pasar por esta resolución.
    let expected_inserted = expected
        .tract_insertion()
        .map(|plan| {
            Ok::<_, Diagnostic>((
                plan.derived_reference_id().to_owned(),
                insertion_changes_json(&expected)?,
            ))
        })
        .transpose()
        .map_err(|error| physical_error(effective_path, "insertion_changes", &error))?;
    let published_inserted = effective
        .get("inserted_reference")
        .and_then(Value::as_object)
        .map(|inserted| {
            (
                inserted
                    .get("id")
                    .and_then(Value::as_str)
                    .map(str::to_owned),
                inserted.get("insertions").cloned(),
            )
        });
    let expected_declared = expected
        .repeat_loci()
        .iter()
        .map(|locus| {
            locus.declared_position().map(|declared| match declared {
                crate::application::DeclaredLocusPosition::Interval { start, end } => {
                    json!({"start": start, "end": end})
                }
                crate::application::DeclaredLocusPosition::InsertAfter(point) => {
                    json!({"insert_after": point})
                }
            })
        })
        .collect::<Vec<_>>();
    let published_declared = effective
        .get("repeat_loci")
        .and_then(Value::as_array)
        .map_or_else(Vec::new, |loci| {
            loci.iter()
                .map(|locus| locus.get("declared").cloned())
                .collect()
        });
    if published_inserted != expected_inserted.map(|(id, changes)| (Some(id), Some(changes)))
        || (expected.document_version() == DeclaredDocumentVersion::V4
            && published_declared != expected_declared)
    {
        return Err(effective_error(
            effective_path,
            "la referencia final, sus cambios o lo declarado de cada locus no coinciden con el TOML",
        ));
    }
    let expected_direct_events = Value::Array(
        expected
            .direct_events()
            .iter()
            .map(direct_event_json)
            .collect(),
    );
    if effective.get("direct_events") != Some(&expected_direct_events) {
        return Err(effective_error(
            effective_path,
            "direct_events efectivo no coincide exactamente con el TOML resuelto",
        ));
    }
    let expected_technical_event = expected
        .technical_event()
        .map(|event| technical_event_json(&expected, event));
    if effective.get("technical_event") != expected_technical_event.as_ref() {
        return Err(effective_error(
            effective_path,
            "technical_event efectivo no coincide exactamente con el TOML resuelto",
        ));
    }
    let expected_vcf = effective_vcf_output_json(&expected);
    if effective.get("vcf_output") != Some(&expected_vcf) {
        return Err(effective_error(
            effective_path,
            "vcf_output efectivo no coincide exactamente con el TOML resuelto",
        ));
    }
    Ok(expected)
}

struct ReplayEntropy(RootSeed);

impl EntropySource for ReplayEntropy {
    fn fill_root_seed(&mut self, destination: &mut [u8; 32]) -> Result<(), String> {
        destination.copy_from_slice(self.0.as_bytes());
        Ok(())
    }
}

struct EffectiveReductionLinks {
    id: String,
    source_reference_id: String,
    fasta_artifact_id: String,
    map_artifact_id: String,
    source_fasta_artifact_id: String,
}

/// Indica si el documento efectivo declara alguna salida que trabaje sobre las
/// bases de la referencia principal.
///
/// Una ejecución que solo deriva referencias reducidas no las necesita: el
/// catálogo basta para comprobar ficha y alfabeto, y los tramos se piden por
/// índice. Cualquier otra salida sí las recorre.
fn effective_needs_reference_bases(effective: &serde_json::Map<String, Value>) -> bool {
    // Las claves ausentes aparecen como `null` en el documento efectivo, así que
    // no basta con comprobar que existan.
    let declared = |key: &str| effective.get(key).is_some_and(|value| !value.is_null());
    let non_empty = |key: &str| {
        effective
            .get(key)
            .and_then(Value::as_array)
            .is_some_and(|values| !values.is_empty())
    };
    // Una referencia simulada se comprueba regenerándola desde la semilla y
    // comparándola base a base: sin bases, esa comparación no podría hacerse.
    let simulated = effective
        .get("reference")
        .and_then(Value::as_object)
        .and_then(|reference| reference.get("source"))
        .and_then(Value::as_str)
        == Some("simulated");
    simulated
        || declared("population_snvs")
        || declared("technical_event")
        || declared("vcf_output")
        || non_empty("initial_populations")
        || non_empty("direct_events")
        || non_empty("vcf_outputs")
}

fn read_physical_reference(
    root: &Path,
    reference: &serde_json::Map<String, Value>,
    artifact_paths: &BTreeMap<&str, &str>,
    effective_path: &Path,
    bases: CatalogBases,
) -> Result<Reference, Diagnostic> {
    let id = effective_non_empty_text(reference, "id", effective_path)?;
    let sequences = effective_array(reference, "sequences", effective_path)?;
    let artifact_id = effective_non_empty_text(reference, "fasta_artifact_id", effective_path)?;
    read_physical_sequence_catalog(
        root,
        id,
        artifact_id,
        sequences,
        "unknown_base_count",
        artifact_paths,
        effective_path,
        bases,
    )
}

/// Modo de lectura del catálogo físico.
///
/// La referencia de partida puede ser un genoma completo del que solo se
/// comprueban ficha y alfabeto; una referencia derivada, en cambio, se compara
/// base a base contra su reconstrucción y necesita cargarse.
#[derive(Clone, Copy, Eq, PartialEq)]
enum CatalogBases {
    Materialized,
    SummaryOnly,
}

#[allow(clippy::too_many_arguments)] // cada argumento es un hecho distinto del catálogo publicado.
fn read_physical_sequence_catalog(
    root: &Path,
    reference_id: &str,
    artifact_id: &str,
    sequences: &[Value],
    n_count_field: &str,
    artifact_paths: &BTreeMap<&str, &str>,
    effective_path: &Path,
    bases: CatalogBases,
) -> Result<Reference, Diagnostic> {
    let mut expectations = Vec::with_capacity(sequences.len());
    for sequence in sequences {
        let sequence = sequence
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        expectations.push(SequenceExpectation {
            id: effective_non_empty_text(sequence, "id", effective_path)?,
            class: effective_sequence_class(sequence, effective_path)?,
            length: effective_positive_u64(sequence, "length", effective_path)?,
        });
    }
    let relative = artifact_paths.get(artifact_id).ok_or_else(|| {
        effective_error(
            effective_path,
            format!("FASTA físico ausente: {artifact_id}"),
        )
    })?;
    let path = root.join(relative);
    // Cuando solo se comprueban ficha y alfabeto basta validar sin conservar las
    // bases; quien las necesite las pide por índice. Una referencia que se
    // compara base a base sí se carga: es pequeña por construcción.
    let observed = match bases {
        CatalogBases::SummaryOnly => {
            let summary = validate_supplied_fasta_summary(&path, reference_id, &expectations)
                .map_err(|error| {
                    physical_error(effective_path, "validate_reference_fasta", &error)
                })?;
            reference_from_summary(
                &summary,
                reference_id,
                crate::reference::ReferenceOrigin::Supplied {
                    original_path: path.clone(),
                },
            )
            .map_err(|error| physical_error(effective_path, "validate_reference_fasta", &error))?
        }
        CatalogBases::Materialized => validate_supplied_fasta(&path, reference_id, &expectations)
            .map_err(|error| physical_error(effective_path, "validate_reference_fasta", &error))?
            .into_reference(),
    };
    for (declared, physical) in sequences.iter().zip(observed.sequences()) {
        let declared = declared
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        let unknown = declared
            .get(n_count_field)
            .and_then(Value::as_u64)
            .ok_or_else(|| effective_error(effective_path, format!("{n_count_field} inválido")))?;
        if unknown != physical.unknown_base_count() {
            return Err(effective_error(
                effective_path,
                format!(
                    "{n_count_field} físico de {} no coincide",
                    physical.id().as_str()
                ),
            ));
        }
        if declared.get("header") != Some(&reversible_bytes_json(physical.header().raw_bytes()))
            || declared.get("description")
                != Some(&reversible_description_json(
                    physical.header().description(),
                ))
            || declared.get("alphabet_validation").and_then(Value::as_str)
                != Some(alphabet_validation(physical.alphabet_validation()))
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "cabecera, descripción o validación de alfabeto física de {} no coincide",
                    physical.id().as_str()
                ),
            ));
        }
    }
    Ok(observed)
}

fn read_reference_with_expected(
    root: &Path,
    reference_id: &str,
    artifact_id: &str,
    declared_sequences: &[Value],
    expected: &[ReferenceSequence],
    artifact_paths: &BTreeMap<&str, &str>,
    effective_path: &Path,
) -> Result<Reference, Diagnostic> {
    let observed = read_physical_sequence_catalog(
        root,
        reference_id,
        artifact_id,
        declared_sequences,
        "n_count",
        artifact_paths,
        effective_path,
        CatalogBases::Materialized,
    )?;
    if observed.sequences().len() != expected.len() {
        return Err(effective_error(
            effective_path,
            "el catálogo físico derivado no coincide con el plan",
        ));
    }
    Ok(observed)
}

/// Un tracto de la referencia simulada tal y como lo publica la configuración
/// efectiva: secuencia, inicio en base cero y bases.
type PublishedTract<'a> = (&'a str, u64, Vec<u8>);

/// Los tractos de una referencia simulada tal y como los publica la
/// configuración efectiva: secuencia, inicio en base cero y bases.
///
/// Un locus de un motivo repite el motivo sobre el tramo, con la última copia
/// parcial si no es múltiplo; uno compuesto concatena cada bloque con su motivo
/// repetido tantas veces como diga `reference_units`, y la suma tiene que
/// cubrir el tramo exactamente.
///
/// **Cada rechazo llama a la constructora en su propia rama**, para que la
/// puerta de cobertura de los jueces los cuente uno a uno: con una sola
/// llamada compartida contaba nueve rechazos como uno y dos mutantes de esta
/// función sobrevivían a la batería entera. Las ramas que este juez puede
/// alcanzar desde un documento publicado las ejercita
/// `every_corruption_of_a_published_tract_is_rejected_by_its_own_branch`.
fn published_simulated_tracts<'a>(
    repeat_loci: Option<&'a Value>,
    physical: &Reference,
    effective_path: &Path,
) -> Result<Vec<PublishedTract<'a>>, Diagnostic> {
    let Some(loci) = repeat_loci else {
        return Ok(Vec::new());
    };
    // **No se alcanza desde un documento publicado**: la comprobación de que la
    // configuración efectiva y el manifiesto registran los mismos loci rechaza
    // antes un `repeat_loci` que no sea una lista. Se conserva porque este
    // lector no puede depender de que esa otra comprobación siga donde está.
    let loci = loci.as_array().ok_or_else(|| {
        effective_error(
            effective_path,
            "tractos de la referencia simulada: repeat_loci no es una lista",
        )
    })?;
    let mut tracts = Vec::with_capacity(loci.len());
    for locus in loci {
        let id = locus.get("id").and_then(Value::as_str).unwrap_or("sin id");
        let detail =
            |what: &str| format!("tracto del locus {id} en la referencia simulada: {what}");
        let sequence_id = locus
            .get("sequence_id")
            .and_then(Value::as_str)
            .ok_or_else(|| effective_error(effective_path, detail("sin secuencia")))?;
        let start = locus
            .get("start")
            .and_then(Value::as_u64)
            .filter(|start| *start >= 1)
            .ok_or_else(|| {
                effective_error(
                    effective_path,
                    detail("el inicio no es una coordenada de base uno"),
                )
            })?;
        let end = locus
            .get("end")
            .and_then(Value::as_u64)
            .filter(|end| *end >= start)
            .ok_or_else(|| {
                effective_error(effective_path, detail("el final no es posterior al inicio"))
            })?;
        let motifs = locus
            .get("motifs")
            .and_then(Value::as_array)
            .and_then(|motifs| {
                motifs
                    .iter()
                    .map(|motif| motif.as_str().filter(|motif| !motif.is_empty()))
                    .collect::<Option<Vec<_>>>()
            })
            .filter(|motifs| !motifs.is_empty())
            .ok_or_else(|| {
                effective_error(
                    effective_path,
                    detail("los motivos no son una lista de textos no vacíos"),
                )
            })?;
        let units = match locus.get("reference_units") {
            None => None,
            Some(value) => Some(
                value
                    .as_array()
                    .and_then(|units| units.iter().map(Value::as_u64).collect::<Option<Vec<_>>>())
                    .ok_or_else(|| {
                        effective_error(
                            effective_path,
                            detail("reference_units no es una lista de enteros"),
                        )
                    })?,
            ),
        };
        // **El tramo cabe en la secuencia publicada, antes de reservar nada con
        // él** (DT-025): `repeat_motif_over` reserva `length` bytes, y un `end`
        // leído del documento sin acotar abortaba el juez con `capacity
        // overflow` en vez de rechazarlo.
        if !physical
            .sequences()
            .iter()
            .any(|sequence| sequence.id().as_str() == sequence_id && end <= sequence.length())
        {
            return Err(effective_error(
                effective_path,
                detail("el tramo no cabe en la secuencia de la referencia publicada"),
            ));
        }
        let length = end - start + 1;
        let bases = match (motifs.as_slice(), units) {
            ([motif], None) => repeat_motif_over(motif.as_bytes(), length),
            ([_], Some(_)) => {
                return Err(effective_error(
                    effective_path,
                    detail("un locus de un solo motivo publica reference_units"),
                ));
            }
            (_, None) => {
                return Err(effective_error(
                    effective_path,
                    detail(
                        "un locus compuesto sobre referencia simulada no publica reference_units",
                    ),
                ));
            }
            (motifs, Some(units)) => {
                if units.len() != motifs.len() {
                    return Err(effective_error(
                        effective_path,
                        detail("reference_units no tiene un valor por motivo"),
                    ));
                }
                let mut bases = Vec::new();
                for (motif, units) in motifs.iter().zip(units) {
                    // Aritmética comprobada y **antes de reservar** (DT-025): un
                    // `reference_units` de `u64::MAX` desbordaba la suma y llegaba a
                    // `repeat_motif_over` con un tamaño imposible. `bases.len()` no
                    // pasa nunca de `length`, de modo que la resta no desborda.
                    let room = length - bases.len() as u64;
                    let Some(block) = units
                        .checked_mul(motif.len() as u64)
                        .filter(|block| *block != 0 && *block <= room)
                    else {
                        return Err(effective_error(
                            effective_path,
                            detail("un bloque de reference_units está vacío o se sale del tramo"),
                        ));
                    };
                    bases.extend(repeat_motif_over(motif.as_bytes(), block));
                }
                if bases.len() as u64 != length {
                    return Err(effective_error(
                        effective_path,
                        detail("los bloques de reference_units no cubren el tramo entero"),
                    ));
                }
                bases
            }
        };
        tracts.push((sequence_id, start - 1, bases));
    }
    Ok(tracts)
}

/// `length` bases del motivo repetido desde su primera base: cada una es la
/// del motivo en la posición que le toca por módulo.
fn repeat_motif_over(motif: &[u8], length: u64) -> Vec<u8> {
    (0..length)
        .map(|offset| motif[(offset % motif.len() as u64) as usize])
        .collect()
}

/// Regenera la referencia simulada desde lo publicado y la compara, base a
/// base, con el FASTA publicado.
///
/// # Los tractos se colocan aquí por cuenta propia
///
/// Desde el 01/10/2026 la referencia simulada lleva colocado el tracto de cada
/// microsatélite antes del relleno al azar (DT-013 y DT-016). Este juez lee la
/// geometría de cada tracto **de la configuración efectiva publicada** —la
/// secuencia, el tramo, los motivos y, en un compuesto, `reference_units`— y lo
/// construye **con su propio cálculo**, sin llamar a la función del productor.
/// Así muere aquí un productor que deje de colocar los tractos, que los
/// desplace, que coloque mal un bloque, que ignore `reference_units` o que lo
/// publique distinto de como lo usó.
///
/// **Lo que sí comparte, dicho:** el relleno lo hace el mismo generador de
/// referencia. Lo que este juez juzga es la colocación de los tractos, no el
/// relleno, que tiene sus propias pruebas y su huella congelada.
fn validate_generated_reference(
    reference: &serde_json::Map<String, Value>,
    root_seed: RootSeed,
    physical: &Reference,
    effective_path: &Path,
    repeat_loci: Option<&Value>,
) -> Result<(), Diagnostic> {
    let tracts = published_simulated_tracts(repeat_loci, physical, effective_path)?;
    let regions = tracts
        .iter()
        .map(
            |(sequence_id, start, bases)| crate::reference_generation::DeclaredRegion {
                sequence_id,
                start: *start,
                bases,
            },
        )
        .collect::<Vec<_>>();
    let sequences = effective_array(reference, "sequences", effective_path)?;
    let mut specs = Vec::with_capacity(sequences.len());
    for sequence in sequences {
        let sequence = sequence
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una secuencia no es objeto"))?;
        specs.push(GeneratedSequenceSpec {
            id: effective_non_empty_text(sequence, "id", effective_path)?,
            class: effective_sequence_class(sequence, effective_path)?,
            length: effective_positive_u64(sequence, "length", effective_path)?,
        });
    }
    let composition = effective_object(reference, "composition", effective_path)?;
    let frequencies = (*b"ACGT")
        .into_iter()
        .zip(["A", "C", "G", "T"])
        .map(|(base, label)| {
            composition
                .get(label)
                .and_then(Value::as_f64)
                .map(|probability| BaseFrequency { base, probability })
                .ok_or_else(|| effective_error(effective_path, "composición efectiva inválida"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let unknown_base_fraction = composition.get("N").map_or(Ok(0.0), |value| {
        value
            .as_f64()
            .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
            .ok_or_else(|| effective_error(effective_path, "fracción N efectiva inválida"))
    })?;
    let expected = crate::reference_generation::generate_reference_with_regions(
        effective_non_empty_text(reference, "id", effective_path)?,
        root_seed,
        &specs,
        &frequencies,
        unknown_base_fraction,
        &regions,
    )
    .map_err(|error| physical_error(effective_path, "regenerate_reference", &error))?;
    compare_reference_bases(&expected, physical, effective_path, expected.id().as_str())
}

fn compare_reference_bases(
    expected: &Reference,
    physical: &Reference,
    effective_path: &Path,
    reference_id: &str,
) -> Result<(), Diagnostic> {
    if expected.id() != physical.id()
        || expected.sequences().len() != physical.sequences().len()
        || expected
            .sequences()
            .iter()
            .zip(physical.sequences())
            .any(|(expected, physical)| {
                expected.id() != physical.id()
                    || expected.class() != physical.class()
                    || expected.length() != physical.length()
                    || expected.unknown_base_count() != physical.unknown_base_count()
                    // Las bases solo se comparan cuando ambas partes las
                    // conservan. Si alguna se validó sin ellas, la
                    // correspondencia con el fichero ya la garantizan su
                    // SHA-256 y el recuento de bases desconocidas, que sí se
                    // comprueban aquí arriba.
                    || match (
                        expected.sequence().as_bytes(),
                        physical.sequence().as_bytes(),
                    ) {
                        (Ok(expected_bases), Ok(physical_bases)) => {
                            expected_bases != physical_bases
                        }
                        _ => false,
                    }
            })
    {
        return Err(effective_error(
            effective_path,
            format!("FASTA físico de {reference_id} no coincide base a base"),
        ));
    }
    Ok(())
}

fn effective_sequence_class(
    sequence: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<SequenceClass, Diagnostic> {
    match sequence.get("class").and_then(Value::as_str) {
        Some("chromosome") => Ok(SequenceClass::Chromosome),
        Some("scaffold") => Ok(SequenceClass::Scaffold),
        _ => Err(effective_error(path, "clase de secuencia inválida")),
    }
}

/// Rechazo de la autoridad física F1. Lleva `#[track_caller]`: sin ella, sus
/// treinta y cuatro llamadores contaban como un solo punto de rechazo para la
/// puerta de cobertura de guardas.
#[track_caller]
fn physical_error(path: &Path, operation: &str, error: &Diagnostic) -> Diagnostic {
    effective_error(
        path,
        format!(
            "autoridad física rechazada en {operation}: {} {}",
            error.code(),
            error.message()
        ),
    )
}

#[track_caller]
fn repeat_section_error(path: &Path, detail: &str) -> Diagnostic {
    manifest_error_for_schema(path, "manifest-general-v6.schema.json", detail.to_owned())
}

/// Lo mismo que [`repeat_section_error`] para un motivo que ya anotó
/// [`reason`]: no anota la línea del envoltorio.
fn repeat_section_error_from_reason(path: &Path, detail: &str) -> Diagnostic {
    manifest_diagnostic(path, "manifest-general-v6.schema.json", detail)
}

/// Comprueba la sección de repeticiones con la misma severidad que el resto
/// del documento.
///
/// # Por qué no basta con que exista
///
/// La primera versión de esta guarda solo exigía que la sección fuera un
/// objeto, y con eso el lector de Rust admitía documentos que su propio
/// esquema publicado rechaza —una sección vacía, por ejemplo— y no cazaba una
/// frecuencia realizada que contradijera a su recuento de copias. El documento
/// vecino, el de las estadísticas de la población inicial, sí comprueba ese
/// par desde hace tiempo; esta sección tiene que ser igual de estricta y no
/// más débil que la que tiene al lado.
///
/// # Lo que el esquema no puede comprobar y aquí sí
///
/// Las cuatro listas de bloque son **paralelas** y miden lo que el recuento de
/// secuencias repetidas declara. Comparar la longitud de varias propiedades
/// hermanas entre sí y contra un entero no es expresable en JSON Schema, de
/// modo que esa correspondencia la sostiene esta función, y la frecuencia
/// realizada se compara **bit a bit** contra copias realizadas entre copias
/// sorteadas, igual que en el documento vecino.
fn validate_repeat_loci_section(value: &Value) -> Result<(), &'static str> {
    let section = value
        .as_object()
        .ok_or_else(|| reason("la sección repeat_loci no es un objeto"))?;
    const FIELDS: [&str; 5] = [
        "reference_id",
        "reference_source",
        "admission_defaults",
        "tract_output",
        "loci",
    ];
    if section.len() != FIELDS.len() || section.keys().any(|key| !FIELDS.contains(&key.as_str())) {
        return Err(reason(
            "la sección repeat_loci no tiene exactamente sus campos",
        ));
    }
    if !section
        .get("reference_id")
        .and_then(Value::as_str)
        .is_some_and(|value| !value.is_empty())
    {
        return Err(reason("la sección repeat_loci no nombra su referencia"));
    }
    if !matches!(
        section.get("reference_source").and_then(Value::as_str),
        Some("simulated" | "supplied")
    ) {
        return Err(reason(
            "la procedencia de la referencia no es simulada ni suministrada",
        ));
    }
    validate_repeat_admission(
        section
            .get("admission_defaults")
            .ok_or_else(|| reason("faltan los mínimos de omisión"))?,
    )?;
    validate_repeat_tract_output(
        section
            .get("tract_output")
            .ok_or_else(|| reason("falta la contabilidad de bases de tracto"))?,
    )?;
    let loci = section
        .get("loci")
        .and_then(Value::as_array)
        .filter(|loci| !loci.is_empty())
        .ok_or_else(|| reason("la sección repeat_loci no publica ningún locus"))?;
    let mut identifiers = BTreeSet::new();
    for locus in loci {
        let locus = locus
            .as_object()
            .ok_or_else(|| reason("un locus publicado no es un objeto"))?;
        const LOCUS_FIELDS: [&str; 10] = [
            "id",
            "sequence_id",
            "start",
            "end",
            "motifs",
            "motif_lengths",
            "max_units",
            "admission",
            "drawn_copies",
            "alleles",
        ];
        // `slippage` es el único campo **opcional** del locus: está cuando el
        // abanico lo generó el deslizamiento biológico y no está cuando los
        // alelos los escribió el usuario. Por eso el recuento se comprueba
        // contra los obligatorios más ese, y no contra un número fijo.
        if LOCUS_FIELDS.iter().any(|field| !locus.contains_key(*field))
            || locus
                .keys()
                .any(|key| !LOCUS_FIELDS.contains(&key.as_str()) && key != "slippage")
        {
            return Err(reason("un locus publicado no tiene exactamente sus campos"));
        }
        if let Some(slippage) = locus.get("slippage") {
            validate_repeat_slippage(slippage)?;
        }
        let identifier = locus
            .get("id")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .ok_or_else(|| reason("un locus publicado no tiene identificador"))?;
        if !identifiers.insert(identifier) {
            return Err(reason("dos loci publicados comparten identificador"));
        }
        if !locus
            .get("sequence_id")
            .and_then(Value::as_str)
            .is_some_and(|value| !value.is_empty())
        {
            return Err(reason("un locus publicado no nombra su secuencia"));
        }
        let start = locus
            .get("start")
            .and_then(Value::as_u64)
            .filter(|start| *start >= 1)
            .ok_or_else(|| reason("el inicio del locus no es una coordenada de base uno"))?;
        let end = locus
            .get("end")
            .and_then(Value::as_u64)
            .filter(|end| *end >= start)
            .ok_or_else(|| reason("el final del locus no llega a su inicio"))?;
        let _ = end;
        if locus.get("max_units").and_then(Value::as_u64).unwrap_or(0) == 0 {
            return Err(reason("el tope del locus no es un entero positivo"));
        }
        validate_repeat_admission(
            locus
                .get("admission")
                .ok_or_else(|| reason("faltan los mínimos del locus"))?,
        )?;
        let motifs = repeat_motif_list(locus.get("motifs"))?;
        let motif_lengths = repeat_u64_list(locus.get("motif_lengths"))?;
        if motif_lengths.len() != motifs.len()
            || motifs
                .iter()
                .zip(&motif_lengths)
                .any(|(motif, length)| motif.len() as u64 != *length)
        {
            return Err(reason(
                "las longitudes de motivo del locus no son las de sus motivos",
            ));
        }
        let drawn = locus
            .get("drawn_copies")
            .and_then(Value::as_u64)
            .ok_or_else(|| reason("las copias sorteadas del locus no son un entero"))?;
        let alleles = locus
            .get("alleles")
            .and_then(Value::as_array)
            .filter(|alleles| !alleles.is_empty())
            .ok_or_else(|| reason("el locus publicado no tiene alelos"))?;
        let generated = locus.contains_key("slippage");
        // **Cada linaje termina en exactamente un alelo**, de modo que no puede
        // haber más alelos que linajes. Es la única relación aritmética que ata
        // la cifra de linajes con la lista de alelos, y sin ella esa cifra se
        // publicaba sin que nada la usara.
        if let Some(lineages) = locus
            .get("slippage")
            .and_then(|slippage| slippage.get("lineages"))
            .and_then(Value::as_u64)
            && alleles.len() as u64 > lineages
        {
            return Err(reason(
                "el abanico publica más alelos que linajes simulados",
            ));
        }
        let mut realized_total = 0_u64;
        for allele in alleles {
            // El alelo se juzga **primero** por sí mismo —que su origen sea uno
            // de los dos que existen, entre lo demás— y solo después se compara
            // con el del locus. Al revés, la comparación cazaba antes cualquier
            // texto inventado y la guarda del dominio no llegaba a dispararse
            // nunca.
            let copies = validate_repeat_allele(allele, drawn)?;
            // **El origen del alelo y el del locus tienen que concordar.** Un
            // locus con generador no puede publicar alelos «declarados», ni al
            // revés: sería una procedencia inventada sobre datos reales.
            let origin = allele
                .get("origin")
                .and_then(Value::as_str)
                .ok_or_else(|| reason("un alelo publicado no declara su origen"))?;
            if generated != (origin == "generated") {
                return Err(reason("el origen del alelo contradice al del locus"));
            }
            realized_total = realized_total
                .checked_add(copies)
                .ok_or_else(|| reason("las copias realizadas del locus desbordan"))?;
        }
        if realized_total != drawn {
            return Err(reason(
                "las copias realizadas no suman las copias sorteadas del locus",
            ));
        }
    }
    Ok(())
}

/// **Juez independiente del abanico de alelos generado por deslizamiento.**
///
/// # Por qué existe, dicho sin rodeos
///
/// La afirmación central del generador es que **la frecuencia de cada alelo es
/// la proporción de linajes que terminan en él**. Antes de que esto existiera,
/// lo único que miraba el abanico publicado era la reconstrucción de la
/// configuración, que **llama al productor**: era el espejo que este proyecto
/// prohíbe por escrito. Se midió con una mutación que emparejaba cada alelo con
/// la frecuencia del simétrico de la lista —el conjunto de alelos no cambia, las
/// frecuencias siguen sumando uno, la contabilidad sigue cuadrando, solo cambia
/// a qué alelo le toca qué— y la batería entera seguía en verde con el espectro
/// publicado del revés.
///
/// # Qué rehace por su cuenta, y qué comparte
///
/// Rehace **el recorrido entero**: su filtro de viabilidad, su regla del borde,
/// su regla de peso, su recuento por alelo y su orden canónico están escritos
/// aquí, aparte del productor. Lee el tramo ancestral de la **referencia
/// publicada**, los parámetros del **documento publicado** y la semilla
/// **publicada**, y después compara.
///
/// **Lo que comparte se dice, en vez de fingir que no lo comparte**: las
/// primitivas de azar —la derivación de la subsemilla, el recuento de Poisson y
/// el ensayo binario— y la operación de paso del dominio con su anotación. Las
/// primitivas tienen sus propios vectores congelados y su propio oráculo
/// independiente; la operación de paso tiene su fichero de pruebas del dominio.
/// Lo que este juez cubre y aquello no es **qué alelo recibe qué linaje**, que
/// es justo lo que la mutación del espectro invertido atacaba.
fn validate_published_slippage_fans(
    effective: &serde_json::Map<String, Value>,
    reference: &Reference,
    root_seed: RootSeed,
    effective_path: &Path,
) -> Result<(), Diagnostic> {
    let Some(loci) = effective.get("repeat_loci").and_then(Value::as_array) else {
        return Ok(());
    };
    // Con inserción, el ámbito de azar es el de la referencia final
    // (DTC-098): la que el productor siembra y la que se juzga.
    let reference_id = effective
        .get("inserted_reference")
        .or_else(|| effective.get("reference"))
        .and_then(Value::as_object)
        .and_then(|reference| reference.get("id"))
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(effective_path, "la referencia publicada no tiene id"))?;
    for locus in loci {
        let Some(slippage) = locus.get("slippage") else {
            continue;
        };
        // Cada `reject(reason(..))` anota su propia línea; el cierre no.
        let reject = |detail: &str| {
            effective_diagnostic(
                effective_path,
                "effective-config-v4.schema.json",
                format!(
                    "abanico generado de {}: {detail}",
                    locus
                        .get("id")
                        .and_then(Value::as_str)
                        .unwrap_or("locus sin id")
                ),
            )
        };
        let locus_id = locus
            .get("id")
            .and_then(Value::as_str)
            .ok_or_else(|| reject(reason("el locus no tiene identificador")))?;
        let ancestral = published_ancestral_tract(locus, reference, &reject)?;
        let max_units = locus
            .get("max_units")
            .and_then(Value::as_u64)
            .ok_or_else(|| reject(reason("el tope del locus no es un entero")))?;
        let mean_events = slippage
            .get("mean_events")
            .and_then(Value::as_f64)
            .ok_or_else(|| reject(reason("la media de pasos no es un decimal")))?;
        let expansion_bias = slippage
            .get("expansion_bias")
            .and_then(Value::as_f64)
            .ok_or_else(|| reject(reason("el sesgo de expansión no es un decimal")))?;
        let lineages = slippage
            .get("lineages")
            .and_then(Value::as_u64)
            .ok_or_else(|| reject(reason("los linajes no son un entero")))?;
        // **Las dos cotas, antes de rehacer el abanico** (DT-024). El
        // recorrido cuesta un paso por linaje y por evento, y ambos números
        // salen del documento: con `lineages` de `u64::MAX` el juez no
        // terminaba. Los linajes son las copias haplotípicas de las
        // poblaciones iniciales, que es lo que el productor usa
        // (`haplotype_copy_count`) y que el documento publica con la lista de
        // sus individuos, de modo que el valor está acotado por lo publicado.
        // La media de pasos tiene el tope que el productor admite.
        if Some(lineages) != published_haplotype_copies(effective) {
            return Err(reject(reason(
                "los linajes no son las copias haplotípicas de las poblaciones iniciales",
            )));
        }
        if !(0.0..=crate::repeat_slippage::SLIPPAGE_MAX_MEAN_EVENTS).contains(&mean_events) {
            return Err(reject(reason(
                "la media de pasos rebasa el tope que el generador admite",
            )));
        }
        let mut rng = crate::repeat_slippage::slippage_subseed(root_seed, reference_id, locus_id)
            .map_err(|error| reject(&reason(error.detail())))?
            .local_rng();
        let recomputed = recompute_slippage_fan(
            &mut rng,
            &ancestral,
            max_units,
            mean_events,
            expansion_bias,
            lineages,
        )
        .map_err(|detail| reject(&detail))?;
        compare_published_fan_against_oracle(locus, slippage, &recomputed, lineages, &reject)?;
    }
    Ok(())
}

/// Las copias haplotípicas de las poblaciones iniciales **publicadas**:
/// individuos por ploidía, sumadas con aritmética comprobada. `None` si algo
/// falta o desborda. Es el mismo recuento que usa el productor
/// (`haplotype_copy_count`), leído del documento.
fn published_haplotype_copies(effective: &serde_json::Map<String, Value>) -> Option<u64> {
    effective
        .get("initial_populations")?
        .as_array()?
        .iter()
        .try_fold(0_u64, |total, population| {
            let individuals = population.get("individual_count")?.as_u64()?;
            let ploidy = population.get("ploidy")?.as_u64()?;
            total.checked_add(individuals.checked_mul(ploidy)?)
        })
}

/// El tramo que la **referencia publicada** lleva en el intervalo del locus, ya
/// estructurado como tracto de su motivo.
fn published_ancestral_tract(
    locus: &Value,
    reference: &Reference,
    reject: &dyn Fn(&str) -> Diagnostic,
) -> Result<TandemRepeatTract, Diagnostic> {
    let sequence_id = locus
        .get("sequence_id")
        .and_then(Value::as_str)
        .ok_or_else(|| reject(reason("el locus no nombra su secuencia")))?;
    let start = locus
        .get("start")
        .and_then(Value::as_u64)
        .filter(|start| *start >= 1)
        .ok_or_else(|| reject(reason("el inicio no es una coordenada de base uno")))?;
    let end = locus
        .get("end")
        .and_then(Value::as_u64)
        .ok_or_else(|| reject(reason("el final no es una coordenada")))?;
    let motifs = locus
        .get("motifs")
        .and_then(Value::as_array)
        .ok_or_else(|| reject(reason("el locus no publica sus motivos")))?;
    let [motif] = motifs.as_slice() else {
        return Err(reject(reason(
            "un locus generado tiene que tener un solo motivo",
        )));
    };
    let motif = motif
        .as_str()
        .ok_or_else(|| reject(reason("el motivo publicado no es texto")))?;
    let motif =
        RepeatMotif::new(motif.as_bytes()).map_err(|error| reject(&reason(error.detail())))?;
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id().as_str() == sequence_id)
        .ok_or_else(|| reject(reason("la referencia publicada no contiene esa secuencia")))?;
    let bases = sequence.sequence().as_bytes().map_err(|_| {
        reject(reason(
            "las bases de la referencia publicada no están materializadas",
        ))
    })?;
    // El intervalo publicado es de base uno cerrada; el corte, de base cero
    // semiabierta. La conversión se escribe aquí una sola vez y a la vista.
    let from = usize::try_from(start - 1).map_err(|_| reject(reason("el inicio desborda")))?;
    let to = usize::try_from(end).map_err(|_| reject(reason("el final desborda")))?;
    let slice = bases
        .get(from..to)
        .ok_or_else(|| reject(reason("el intervalo rebasa la secuencia publicada")))?;
    TandemRepeatTract::new(motif, slice).map_err(|error| reject(&reason(error.detail())))
}

/// El recorrido **escrito aparte del productor**: cada linaje parte del
/// ancestro, recibe su número de pasos y camina hasta detenerse en un borde.
///
/// Devuelve los alelos en orden canónico —por longitud del tracto y, a igual
/// longitud, por sus bases— con los linajes que terminaron en cada uno, y la
/// contabilidad de los eventos.
fn recompute_slippage_fan(
    rng: &mut crate::randomness::LocalRng,
    ancestral: &TandemRepeatTract,
    max_units: u64,
    mean_events: f64,
    expansion_bias: f64,
    lineages: u64,
) -> Result<RecomputedFan, String> {
    let mut counts: BTreeMap<(usize, Vec<u8>), u64> = BTreeMap::new();
    let mut drawn = 0_u64;
    let mut applied = 0_u64;
    let mut unapplied = 0_u64;
    for _ in 0..lineages {
        let events = crate::randomness::poisson_count(rng, mean_events)
            .map_err(|error| reason(format!("el recuento de pasos falló: {error:?}")))?;
        drawn = drawn
            .checked_add(events)
            .ok_or_else(|| reason("el recuento de pasos desborda".to_owned()))?;
        let mut current = ancestral.clone();
        let mut stopped = false;
        for _ in 0..events {
            if stopped {
                unapplied += 1;
                continue;
            }
            // Viabilidad, escrita aquí: la contracción exige dos unidades en la
            // tirada más larga y la expansión tiene que caber bajo el tope.
            let annotation = current.annotation();
            let can_contract = annotation.longest_run_units() >= 2;
            let can_expand = annotation.units_as_f64() + 1.0 <= max_units as f64;
            if !(can_contract && can_expand) {
                stopped = true;
                unapplied += 1;
                continue;
            }
            let step = if expansion_bias == 1.0 {
                1
            } else if expansion_bias == 0.0 {
                -1
            } else if crate::randomness::f4_bernoulli(rng, expansion_bias)
                .map_err(|error| reason(format!("el ensayo del sesgo falló: {error:?}")))?
            {
                1
            } else {
                -1
            };
            current = current
                .apply_unit_step(step)
                .map_err(|error| reason(format!("el paso falló: {}", error.detail())))?;
            applied += 1;
        }
        *counts
            .entry((current.bases().len(), current.bases().to_vec()))
            .or_insert(0) += 1;
    }
    Ok(RecomputedFan {
        alleles: counts
            .into_iter()
            .map(|((_, bases), count)| (bases, count))
            .collect(),
        drawn,
        applied,
        unapplied,
    })
}

struct RecomputedFan {
    alleles: Vec<(Vec<u8>, u64)>,
    drawn: u64,
    applied: u64,
    unapplied: u64,
}

/// Contrasta el abanico **escrito** contra el **rehecho**, alelo a alelo.
fn compare_published_fan_against_oracle(
    locus: &Value,
    slippage: &Value,
    recomputed: &RecomputedFan,
    lineages: u64,
    reject: &dyn Fn(&str) -> Diagnostic,
) -> Result<(), Diagnostic> {
    let published = locus
        .get("alleles")
        .and_then(Value::as_array)
        .ok_or_else(|| reject(reason("el locus no publica alelos")))?;
    if published.len() != recomputed.alleles.len() {
        return Err(reject(reason(
            "el número de alelos publicados no es el que el modelo produce",
        )));
    }
    for (index, (expected_bases, expected_lineages)) in recomputed.alleles.iter().enumerate() {
        let allele = &published[index];
        let sequence = allele
            .get("sequence")
            .and_then(Value::as_str)
            .ok_or_else(|| reject(reason("un alelo publicado no tiene secuencia")))?;
        if sequence.as_bytes() != expected_bases.as_slice() {
            return Err(reject(reason(
                "la secuencia de un alelo publicado no es la que el modelo produce, o el orden \
canónico no es el suyo",
            )));
        }
        // **Aquí muere el espectro emparejado al revés.** La frecuencia se
        // compara bit a bit contra los linajes que este juez contó terminando
        // en ESTE alelo, no contra el conjunto de frecuencias.
        let expected_frequency = *expected_lineages as f64 / lineages as f64;
        let published_frequency = allele.get("declared_frequency").and_then(Value::as_f64);
        if published_frequency != Some(expected_frequency) {
            // El diagnóstico lleva las dos cifras: sin ellas, quien lo lea sabe
            // que algo no cuadra y no **qué**.
            return Err(reject(&reason(format!(
                "la frecuencia de un alelo no es la proporción de linajes que terminan en él \
(alelo {index}: publicada {published_frequency:?}, rehecha {expected_frequency} = \
{expected_lineages}/{lineages})"
            ))));
        }
    }
    for (field, expected) in [
        ("drawn_events", recomputed.drawn),
        ("applied_events", recomputed.applied),
        ("unapplied_events", recomputed.unapplied),
    ] {
        if slippage.get(field).and_then(Value::as_u64) != Some(expected) {
            return Err(reject(reason(
                "la contabilidad publicada del generador no es la que el modelo produce",
            )));
        }
    }
    Ok(())
}

/// Comprueba el bloque del generador de alelos por deslizamiento.
///
/// # Qué exige, y por qué cada cosa
///
/// El modelo y el tamaño de paso son **literales fijos**: cualquier otro valor
/// describiría un modelo que el simulador no tiene. La media y el sesgo son
/// decimales en su dominio. Las dos versiones —la del algoritmo y la del
/// catálogo de ámbitos— y el nombre del ámbito tienen que estar, porque sin
/// ellas dos ejecuciones con la misma semilla y resultados distintos serían
/// indistinguibles de un defecto.
///
/// La contabilidad se comprueba **entre sí**: los pasos aplicados más los no
/// aplicados tienen que sumar exactamente los sorteados. Un servicio que
/// perdiera eventos por el camino, o que contase dos veces, muere aquí y no en
/// una cifra suelta que nadie puede contrastar.
///
/// Las tres limitaciones se exigen **una a una y con su texto**: son decisiones
/// biológicas del modelo, y un documento que las recortara dejaría a quien lo
/// lea creyendo que el conjunto sirve para lo que no sirve.
fn validate_repeat_slippage(value: &Value) -> Result<(), &'static str> {
    let slippage = value
        .as_object()
        .ok_or_else(|| reason("el bloque del generador no es un objeto"))?;
    const FIELDS: [&str; 14] = [
        "model",
        "step_units",
        "mean_events",
        "expansion_bias",
        "algorithm",
        "count_algorithm",
        "bias_algorithm",
        "rng_scope",
        "rng_scope_catalog_version",
        "limitations",
        "lineages",
        "drawn_events",
        "applied_events",
        "unapplied_events",
    ];
    if slippage.len() != FIELDS.len() || slippage.keys().any(|key| !FIELDS.contains(&key.as_str()))
    {
        return Err(reason(
            "el bloque del generador no tiene exactamente sus campos",
        ));
    }
    if slippage.get("model").and_then(Value::as_str)
        != Some(crate::repeat_slippage::SLIPPAGE_MODEL_STEPWISE)
    {
        return Err(reason("el modelo del generador no es el escalonado"));
    }
    if slippage.get("step_units").and_then(Value::as_u64)
        != Some(crate::repeat_slippage::SLIPPAGE_STEP_UNITS)
    {
        return Err(reason("el paso del generador no es de una unidad completa"));
    }
    let mean = slippage
        .get("mean_events")
        .and_then(Value::as_f64)
        .ok_or_else(|| reason("la media de pasos no es un decimal"))?;
    if !mean.is_finite() || mean < 0.0 {
        return Err(reason("la media de pasos no es un decimal no negativo"));
    }
    let bias = slippage
        .get("expansion_bias")
        .and_then(Value::as_f64)
        .ok_or_else(|| reason("el sesgo de expansión no es un decimal"))?;
    if !bias.is_finite() || !(0.0..=1.0).contains(&bias) {
        return Err(reason("el sesgo de expansión no pertenece a [0,1]"));
    }
    for (field, expected) in [
        ("algorithm", crate::repeat_slippage::SLIPPAGE_ALGORITHM),
        (
            "count_algorithm",
            crate::repeat_slippage::SLIPPAGE_COUNT_ALGORITHM,
        ),
        (
            "bias_algorithm",
            crate::repeat_slippage::SLIPPAGE_BIAS_ALGORITHM,
        ),
        ("rng_scope", crate::repeat_slippage::SLIPPAGE_SCOPE),
        (
            "rng_scope_catalog_version",
            crate::repeat_slippage::SLIPPAGE_SCOPE_CATALOG_VERSION,
        ),
    ] {
        if slippage.get(field).and_then(Value::as_str) != Some(expected) {
            return Err(reason(
                "una versión o un ámbito del generador no es el normativo",
            ));
        }
    }
    let limitations = slippage
        .get("limitations")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("las limitaciones del generador no son una lista"))?;
    if limitations.len() != crate::repeat_slippage::SLIPPAGE_LIMITATIONS.len()
        || limitations
            .iter()
            .zip(crate::repeat_slippage::SLIPPAGE_LIMITATIONS)
            .any(|(published, expected)| published.as_str() != Some(expected))
    {
        return Err(reason(
            "las limitaciones del generador no son las tres declaradas",
        ));
    }
    if slippage
        .get("lineages")
        .and_then(Value::as_u64)
        .is_none_or(|lineages| lineages == 0)
    {
        return Err(reason(
            "los linajes del generador no son un entero positivo",
        ));
    }
    let drawn = slippage
        .get("drawn_events")
        .and_then(Value::as_u64)
        .ok_or_else(|| reason("los pasos sorteados no son un entero"))?;
    let applied = slippage
        .get("applied_events")
        .and_then(Value::as_u64)
        .ok_or_else(|| reason("los pasos aplicados no son un entero"))?;
    let unapplied = slippage
        .get("unapplied_events")
        .and_then(Value::as_u64)
        .ok_or_else(|| reason("los pasos no aplicados no son un entero"))?;
    if applied.checked_add(unapplied) != Some(drawn) {
        return Err(reason(
            "los pasos aplicados y los no aplicados no suman los sorteados",
        ));
    }
    Ok(())
}

/// Comprueba la forma de la contabilidad de tracto.
///
/// Aquí solo se juzga lo que el documento puede decir de sí mismo: que las seis
/// cifras existan, que cada total sea la suma de sus dos términos y que lo
/// realizado no pase de la cota. **Que lo realizado sea lo que los ficheros
/// llevan no se puede comprobar desde el manifiesto**, y no se finge que sí:
/// eso lo hace el juez que relee el `truth.jsonl` y el VCF.
fn validate_repeat_tract_output(value: &Value) -> Result<(), &'static str> {
    let accounting = value
        .as_object()
        .ok_or_else(|| reason("la contabilidad de tracto no es un objeto"))?;
    const FIELDS: [&str; 6] = [
        "truth_tract_bases_bound",
        "vcf_tract_bases_bound",
        "total_tract_bases_bound",
        "realized_truth_tract_bases",
        "realized_vcf_tract_bases",
        "realized_total_tract_bases",
    ];
    if accounting.len() != FIELDS.len()
        || accounting.keys().any(|key| !FIELDS.contains(&key.as_str()))
    {
        return Err(reason(
            "la contabilidad de tracto no tiene exactamente sus campos",
        ));
    }
    let read = |key: &str| {
        accounting
            .get(key)
            .and_then(Value::as_u64)
            .ok_or_else(|| reason("una cifra de la contabilidad de tracto no es un entero"))
    };
    let truth_bound = read("truth_tract_bases_bound")?;
    let vcf_bound = read("vcf_tract_bases_bound")?;
    let total_bound = read("total_tract_bases_bound")?;
    let truth = read("realized_truth_tract_bases")?;
    let vcf = read("realized_vcf_tract_bases")?;
    let total = read("realized_total_tract_bases")?;
    if truth_bound.checked_add(vcf_bound) != Some(total_bound) {
        return Err(reason("la cota total no es la suma de sus dos términos"));
    }
    if truth.checked_add(vcf) != Some(total) {
        return Err(reason(
            "el total realizado no es la suma de sus dos términos",
        ));
    }
    if truth > truth_bound || vcf > vcf_bound {
        return Err(reason("las bases de tracto realizadas superan su cota"));
    }
    Ok(())
}

fn validate_repeat_admission(value: &Value) -> Result<(), &'static str> {
    let admission = value
        .as_object()
        .ok_or_else(|| reason("unos mínimos de admisión no son un objeto"))?;
    const FIELDS: [&str; 2] = ["min_purity", "min_longest_run_units"];
    if admission.len() != FIELDS.len()
        || admission.keys().any(|key| !FIELDS.contains(&key.as_str()))
    {
        return Err(reason(
            "unos mínimos de admisión no tienen exactamente sus campos",
        ));
    }
    let purity = admission
        .get("min_purity")
        .and_then(Value::as_f64)
        .ok_or_else(|| reason("la pureza mínima no es un decimal"))?;
    if !purity.is_finite() || !(0.0..=1.0).contains(&purity) {
        return Err(reason("la pureza mínima no es finita ni pertenece a [0,1]"));
    }
    if admission
        .get("min_longest_run_units")
        .and_then(Value::as_u64)
        .unwrap_or(0)
        == 0
    {
        return Err(reason("la tirada mínima no es un entero positivo"));
    }
    Ok(())
}

/// Valida un alelo publicado y devuelve sus copias realizadas.
fn validate_repeat_allele(value: &Value, drawn: u64) -> Result<u64, &'static str> {
    let allele = value
        .as_object()
        .ok_or_else(|| reason("un alelo publicado no es un objeto"))?;
    const FIELDS: [&str; 15] = [
        "sequence",
        "origin",
        "declared_frequency",
        "effective_frequency",
        "realized_frequency",
        "realized_copies",
        "units",
        "purity",
        "longest_run_units",
        "repeat_sequence_count",
        "block_motifs",
        "block_motif_lengths",
        "block_unit_counts",
        "block_bases",
        "interruptions",
    ];
    if allele.len() != FIELDS.len() || allele.keys().any(|key| !FIELDS.contains(&key.as_str())) {
        return Err(reason("un alelo publicado no tiene exactamente sus campos"));
    }
    if !allele
        .get("sequence")
        .and_then(Value::as_str)
        .is_some_and(|sequence| !sequence.is_empty() && sequence.bytes().all(is_repeat_allele_base))
    {
        return Err(reason("la secuencia del alelo no son bases"));
    }
    if !matches!(
        allele.get("origin").and_then(Value::as_str),
        Some("declared" | "generated")
    ) {
        return Err(reason("el origen del alelo no es declarado ni generado"));
    }
    for field in ["declared_frequency", "effective_frequency", "purity"] {
        let value = allele
            .get(field)
            .and_then(Value::as_f64)
            .ok_or_else(|| reason("una proporción del alelo no es un decimal"))?;
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(reason(
                "una proporción del alelo no es finita ni pertenece a [0,1]",
            ));
        }
    }
    let units = allele
        .get("units")
        .and_then(Value::as_f64)
        .ok_or_else(|| reason("las unidades del alelo no son un decimal"))?;
    if !units.is_finite() || units <= 0.0 {
        return Err(reason("las unidades del alelo no son un decimal positivo"));
    }
    if allele
        .get("longest_run_units")
        .and_then(Value::as_u64)
        .is_none()
    {
        return Err(reason("la tirada más larga del alelo no es un entero"));
    }
    let copies = allele
        .get("realized_copies")
        .and_then(Value::as_u64)
        .ok_or_else(|| reason("las copias realizadas del alelo no son un entero"))?;
    if copies > drawn {
        return Err(reason(
            "un alelo recibió más copias que las sorteadas del locus",
        ));
    }
    // **La misma comparación que el documento vecino**: bit a bit contra
    // copias realizadas entre copias sorteadas. Sin sorteo no hay proporción, y
    // publicar un cero afirmaría una medida que no se hizo.
    let realized = allele
        .get("realized_frequency")
        .ok_or_else(|| reason("falta la frecuencia realizada del alelo"))?;
    if drawn == 0 {
        if !realized.is_null() {
            return Err(reason(
                "sin copias sorteadas la frecuencia realizada tiene que ser nula",
            ));
        }
    } else {
        let realized = realized
            .as_f64()
            .ok_or_else(|| reason("la frecuencia realizada del alelo no es un decimal"))?;
        let expected = copies as f64 / drawn as f64;
        if realized.to_bits() != expected.to_bits() {
            return Err(reason(
                "la frecuencia realizada no es las copias realizadas entre las sorteadas",
            ));
        }
    }
    let count = allele
        .get("repeat_sequence_count")
        .and_then(Value::as_u64)
        .filter(|count| *count >= 1)
        .and_then(|count| usize::try_from(count).ok())
        .ok_or_else(|| reason("el recuento de secuencias repetidas no es un entero positivo"))?;
    let motifs = repeat_motif_list(allele.get("block_motifs"))?;
    let motif_lengths = repeat_u64_list(allele.get("block_motif_lengths"))?;
    let bases = repeat_u64_list(allele.get("block_bases"))?;
    let unit_counts = allele
        .get("block_unit_counts")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("las copias por bloque no son una lista"))?
        .iter()
        .map(|value| value.as_f64())
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| reason("una cifra de copias por bloque no es un decimal"))?;
    if motifs.len() != count
        || motif_lengths.len() != count
        || bases.len() != count
        || unit_counts.len() != count
    {
        return Err(reason(
            "las listas de bloque no miden lo que el recuento declara",
        ));
    }
    let mut total = 0_u64;
    for index in 0..count {
        if motifs[index].len() as u64 != motif_lengths[index] {
            return Err(reason(
                "la longitud de un motivo de bloque no es la de su secuencia",
            ));
        }
        if bases[index] == 0 {
            return Err(reason("un bloque sin bases no es un bloque"));
        }
        let derived = bases[index] as f64 / motif_lengths[index] as f64;
        if !unit_counts[index].is_finite()
            || (unit_counts[index] - derived).abs() > 1e-9 * derived.max(1.0)
        {
            return Err(reason(
                "las copias de un bloque no son sus bases entre la longitud del motivo",
            ));
        }
        total = total
            .checked_add(bases[index])
            .ok_or_else(|| reason("las bases de los bloques desbordan"))?;
    }
    if total
        != allele
            .get("sequence")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .len() as u64
    {
        return Err(reason(
            "las bases de los bloques no suman las de la secuencia del alelo",
        ));
    }
    for interruption in allele
        .get("interruptions")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("las interrupciones del alelo no son una lista"))?
    {
        let interruption = interruption
            .as_object()
            .ok_or_else(|| reason("una interrupción publicada no es un objeto"))?;
        const INTERRUPTION_FIELDS: [&str; 2] = ["offset", "length"];
        if interruption.len() != INTERRUPTION_FIELDS.len()
            || interruption
                .keys()
                .any(|key| !INTERRUPTION_FIELDS.contains(&key.as_str()))
        {
            return Err(reason(
                "una interrupción publicada no tiene exactamente sus campos",
            ));
        }
        if interruption.get("offset").and_then(Value::as_u64).is_none()
            || interruption
                .get("length")
                .and_then(Value::as_u64)
                .unwrap_or(0)
                == 0
        {
            return Err(reason(
                "una interrupción publicada no tiene desplazamiento y longitud válidos",
            ));
        }
    }
    Ok(copies)
}

fn repeat_motif_list(value: Option<&Value>) -> Result<Vec<&str>, &'static str> {
    let motifs = value
        .and_then(Value::as_array)
        .filter(|motifs| !motifs.is_empty())
        .ok_or_else(|| reason("una lista de motivos está ausente o vacía"))?
        .iter()
        .map(Value::as_str)
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| reason("un motivo publicado no es texto"))?;
    if motifs
        .iter()
        .any(|motif| motif.is_empty() || !motif.bytes().all(is_repeat_motif_base))
    {
        return Err(reason("un motivo publicado no son bases A, C, G o T"));
    }
    Ok(motifs)
}

fn repeat_u64_list(value: Option<&Value>) -> Result<Vec<u64>, &'static str> {
    value
        .and_then(Value::as_array)
        .filter(|values| !values.is_empty())
        .ok_or_else(|| reason("una lista de enteros está ausente o vacía"))?
        .iter()
        .map(Value::as_u64)
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| reason("un entero publicado no lo es"))
}

const fn is_repeat_motif_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

const fn is_repeat_allele_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T' | b'N')
}

pub(crate) fn parse_general_manifest_bytes(
    path: &Path,
    bytes: &[u8],
) -> Result<Vec<ArtifactRecord>, Diagnostic> {
    let schema_version = serde_json::from_slice::<Value>(bytes)
        .ok()
        .and_then(|document| {
            document
                .get("schema_version")
                .and_then(Value::as_str)
                .map(str::to_owned)
        });
    parse_general_manifest_bytes_unretagged(path, bytes).map_err(|error| {
        if is_supported_general_manifest_version(schema_version.as_deref()) {
            retag_format_schema(
                error,
                general_manifest_schema_file_for(schema_version.as_deref()),
            )
        } else {
            error
        }
    })
}

fn parse_general_manifest_bytes_unretagged(
    path: &Path,
    bytes: &[u8],
) -> Result<Vec<ArtifactRecord>, Diagnostic> {
    let document: Value = serde_json::from_slice(bytes)
        .map_err(|error| manifest_error(path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| manifest_error(path, "la raíz no es un objeto"))?;
    let schema_version = object.get("schema_version").and_then(Value::as_str);
    if !is_supported_general_manifest_version(schema_version) {
        return Err(manifest_error_for_schema(
            path,
            general_manifest_schema_file_for(schema_version),
            format!(
                "versión de manifiesto no admitida: {schema_version:?}; se admiten exclusivamente \"{GENERAL_MANIFEST_SCHEMA_VERSION}\" ({}) y \"{GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS}\" ({})",
                general_manifest_schema_file_for(Some(GENERAL_MANIFEST_SCHEMA_VERSION)),
                general_manifest_schema_file_for(Some(
                    GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS
                )),
            ),
        ));
    }
    const FIELDS: [&str; 16] = [
        "manifest_kind",
        "schema_version",
        "execution_id",
        "state",
        "configuration_input",
        "effective_config_artifact_id",
        "references",
        "reference_reductions",
        "population_snvs",
        "vcf_outputs",
        "authorizations",
        "output_categories",
        "artifacts",
        "diagnostics",
        "events",
        "technical_events",
    ];
    // La sección de repeticiones **pertenece a la versión 6 y solo a ella**: su
    // presencia en un documento v4 y su ausencia en uno v6 son errores de
    // forma, no un campo opcional. Se añade a la lista cerrada según la versión
    // que el documento declara, no según si aparece.
    let mut fields = FIELDS.to_vec();
    if schema_version == Some(GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS) {
        fields.push("repeat_loci");
        // La inserción de tractos exige loci, de modo que su sección solo
        // puede vivir en la versión 6; es opcional dentro de ella (DTC-098).
        fields.push("reference_insertion");
    }
    validate_closed_object(object, &fields, "manifest.json", path).map_err(|error| {
        retag_format_schema(error, general_manifest_schema_file_for(schema_version))
    })?;
    if object.get("manifest_kind").and_then(Value::as_str) != Some("gps_execution_manifest") {
        return Err(manifest_error(path, "clase o versión incompatible"));
    }
    if schema_version == Some(GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS) {
        let section = object.get("repeat_loci").ok_or_else(|| {
            repeat_section_error(path, "la versión 6 exige la sección repeat_loci")
        })?;
        validate_repeat_loci_section(section)
            .map_err(|detail| repeat_section_error_from_reason(path, detail))?;
    }
    required_non_empty_text(object, "execution_id", path)?;
    let state = match object.get("state").and_then(Value::as_str) {
        Some("complete") => UnitState::Complete,
        Some("failed") => UnitState::Failed,
        _ => return Err(manifest_error(path, "estado de ejecución inválido")),
    };

    let values = required_array(object, "artifacts", path)?;
    let mut ids = BTreeSet::new();
    let mut paths = BTreeSet::new();
    let mut records = Vec::with_capacity(values.len());
    for value in values {
        let record = validate_artifact(value, path)?;
        if !ids.insert(record.id.clone()) || !paths.insert(record.path.clone()) {
            return Err(manifest_error(path, "ID o ruta de artefacto repetidos"));
        }
        records.push(record);
    }
    let catalog = records
        .iter()
        .map(|record| (record.id.as_str(), record.kind.as_str()))
        .collect::<BTreeMap<_, _>>();
    if state == UnitState::Complete {
        require_truth_artifact(&records, path)?;
    }

    let configuration_input = object
        .get("configuration_input")
        .ok_or_else(|| manifest_error(path, "falta configuration_input"))?;
    let configuration_artifact_id = match configuration_input {
        Value::Null => None,
        Value::Object(value) => Some(validate_configuration_input(value, &catalog, path)?),
        _ => return Err(manifest_error(path, "configuration_input inválido")),
    };
    let effective_id = optional_text(object, "effective_config_artifact_id", path)?;
    if let Some(id) = effective_id {
        require_artifact_kind(&catalog, id, "effective_config", path)?;
    }
    let references =
        validate_references(required_array(object, "references", path)?, &catalog, path)?;
    let vcf_requested = validate_vcf_outputs(
        required_array(object, "vcf_outputs", path)?,
        &catalog,
        &references,
        schema_version == Some(GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS),
        path,
    )?;
    let reductions = validate_reference_reductions(
        required_array(object, "reference_reductions", path)?,
        &catalog,
        &references,
        path,
    )?;
    let insertion = object
        .get("reference_insertion")
        .map(|value| validate_manifest_reference_insertion(value, &catalog, &references, path))
        .transpose()?;
    let population_snvs =
        validate_manifest_population_snvs(object.get("population_snvs"), &catalog, path)?;
    let events = validate_manifest_events(object.get("events"), &records, &catalog, path)?;
    validate_manifest_technical_events(
        required_array(object, "technical_events", path)?,
        &events,
        effective_id,
        &records,
        path,
    )?;
    validate_manifest_authorizations(object.get("authorizations"), path)?;
    if let Some(id) = population_snvs.as_deref() {
        let declared = object
            .get("population_snvs")
            .and_then(Value::as_object)
            .and_then(|link| link.get("metadata_sha256"))
            .and_then(Value::as_str);
        let observed = records
            .iter()
            .find(|artifact| artifact.id == id)
            .map(|artifact| artifact.digest.as_str());
        if declared != observed {
            return Err(manifest_error(
                path,
                "population_snvs.metadata_sha256 no coincide con el artefacto",
            ));
        }
    }
    let mut expected_minimum = BTreeSet::new();
    if let Some(id) = configuration_artifact_id.as_deref() {
        expected_minimum.insert(id);
    }
    if let Some(id) = effective_id {
        expected_minimum.insert(id);
    }
    if let Some(id) = population_snvs.as_deref() {
        expected_minimum.insert(id);
        require_artifact_kind(&catalog, "population-snvs-tsv", "population_snvs_tsv", path)?;
        expected_minimum.insert("population-snvs-tsv");
    }
    if state == UnitState::Complete {
        expected_minimum.insert("truth-jsonl");
    }
    if !events.is_empty() {
        expected_minimum.extend(
            records
                .iter()
                .filter(|artifact| artifact.kind == "fastq" || artifact.kind == "event_manifest")
                .map(|artifact| artifact.id.as_str()),
        );
    }
    for (reference_id, artifact_id) in &references {
        if !reductions.derived_reference_ids.contains(reference_id) {
            expected_minimum.insert(artifact_id.as_str());
        }
    }
    // El mapa de la referencia final es parte del mínimo, como su FASTA: sin
    // él la simulación no puede relacionarse con el genoma suministrado.
    if let Some(insertion) = &insertion {
        expected_minimum.insert(insertion.map_artifact_id.as_str());
    }
    validate_output_categories(
        required_object(object, "output_categories", path)?,
        &catalog,
        &expected_minimum,
        &reductions
            .requested_artifact_ids
            .union(&vcf_requested)
            .cloned()
            .collect(),
        path,
    )?;
    let diagnostics = required_array(object, "diagnostics", path)?;
    validate_diagnostics(diagnostics, "3", path)?;

    match state {
        UnitState::Complete => {
            // Rama propia, y no una condición más de la de abajo: así la
            // puerta de cobertura la cuenta por separado.
            if insertion
                .as_ref()
                .is_some_and(|insertion| insertion.state != UnitState::Complete)
            {
                return Err(manifest_error(
                    path,
                    "una ejecución complete exige que su referencia final por inserción esté complete",
                ));
            }
            if configuration_input.is_null()
                || effective_id.is_none()
                || required_array(object, "references", path)?.is_empty()
                || reductions
                    .states
                    .iter()
                    .any(|state| *state != UnitState::Complete)
                || !diagnostics.is_empty()
                || events
                    .iter()
                    .any(|(_, state, _)| *state != UnitState::Complete)
            {
                return Err(manifest_error(
                    path,
                    "una ejecución complete exige entrada, configuración efectiva, referencia y ningún diagnóstico",
                ));
            }
        }
        UnitState::Failed => {
            if diagnostics.is_empty() {
                return Err(manifest_error(
                    path,
                    "una ejecución failed exige al menos un diagnóstico",
                ));
            }
        }
        UnitState::Incomplete => unreachable!("manifest parser never constructs incomplete"),
    }
    Ok(records)
}

/// Reads the frozen F1 manifest contract.  It deliberately does not route a
/// v1 document through F2 requirements such as `truth.jsonl`, population-SNV
/// metadata, authorizations, or VCF outputs.
#[allow(dead_code)]
fn parse_general_manifest_v1_object(
    path: &Path,
    object: &serde_json::Map<String, Value>,
) -> Result<Vec<ArtifactRecord>, Diagnostic> {
    const FIELDS: [&str; 11] = [
        "manifest_kind",
        "schema_version",
        "execution_id",
        "state",
        "configuration_input",
        "effective_config_artifact_id",
        "references",
        "reference_reductions",
        "output_categories",
        "artifacts",
        "diagnostics",
    ];
    validate_closed_object(object, &FIELDS, "manifest.json", path)?;
    if object.get("manifest_kind").and_then(Value::as_str) != Some("gps_execution_manifest")
        || object.get("schema_version").and_then(Value::as_str) != Some("1")
    {
        return Err(manifest_error(path, "clase o versión incompatible"));
    }
    required_non_empty_text(object, "execution_id", path)?;
    let state = match object.get("state").and_then(Value::as_str) {
        Some("complete") => UnitState::Complete,
        Some("failed") => UnitState::Failed,
        _ => return Err(manifest_error(path, "estado de ejecución inválido")),
    };

    let values = required_array(object, "artifacts", path)?;
    let mut ids = BTreeSet::new();
    let mut paths = BTreeSet::new();
    let mut records = Vec::with_capacity(values.len());
    for value in values {
        let record = validate_artifact(value, path)?;
        if !ids.insert(record.id.clone()) || !paths.insert(record.path.clone()) {
            return Err(manifest_error(path, "ID o ruta de artefacto repetidos"));
        }
        records.push(record);
    }
    let catalog = records
        .iter()
        .map(|record| (record.id.as_str(), record.kind.as_str()))
        .collect::<BTreeMap<_, _>>();
    let configuration_input = object
        .get("configuration_input")
        .ok_or_else(|| manifest_error(path, "falta configuration_input"))?;
    let configuration_artifact_id = match configuration_input {
        Value::Null => None,
        Value::Object(value) => Some(validate_configuration_input(value, &catalog, path)?),
        _ => return Err(manifest_error(path, "configuration_input inválido")),
    };
    let effective_id = optional_text(object, "effective_config_artifact_id", path)?;
    if let Some(id) = effective_id {
        require_artifact_kind(&catalog, id, "effective_config", path)?;
    }
    let references =
        validate_references(required_array(object, "references", path)?, &catalog, path)?;
    let reductions = validate_reference_reductions(
        required_array(object, "reference_reductions", path)?,
        &catalog,
        &references,
        path,
    )?;
    let mut expected_minimum = BTreeSet::new();
    if let Some(id) = configuration_artifact_id.as_deref() {
        expected_minimum.insert(id);
    }
    if let Some(id) = effective_id {
        expected_minimum.insert(id);
    }
    for (reference_id, artifact_id) in &references {
        if !reductions.derived_reference_ids.contains(reference_id) {
            expected_minimum.insert(artifact_id.as_str());
        }
    }
    validate_output_categories(
        required_object(object, "output_categories", path)?,
        &catalog,
        &expected_minimum,
        &reductions.requested_artifact_ids,
        path,
    )?;
    let diagnostics = required_array(object, "diagnostics", path)?;
    validate_diagnostics(diagnostics, "3", path)?;
    match state {
        UnitState::Complete
            if configuration_input.is_null()
                || effective_id.is_none()
                || required_array(object, "references", path)?.is_empty()
                || reductions
                    .states
                    .iter()
                    .any(|state| *state != UnitState::Complete)
                || !diagnostics.is_empty() =>
        {
            Err(manifest_error(
                path,
                "una ejecución complete exige entrada, configuración efectiva, referencia y ningún diagnóstico",
            ))
        }
        UnitState::Failed if diagnostics.is_empty() => Err(manifest_error(
            path,
            "una ejecución failed exige al menos un diagnóstico",
        )),
        UnitState::Incomplete => unreachable!("manifest parser never constructs incomplete"),
        _ => Ok(records),
    }
}

#[allow(dead_code)]
fn validate_effective_config_v1_against_manifest_bytes(
    effective_path: &Path,
    effective_bytes: &[u8],
    manifest_path: &Path,
    manifest_bytes: &[u8],
) -> Result<(), Diagnostic> {
    let artifacts = parse_general_manifest_bytes(manifest_path, manifest_bytes)?;
    let manifest: Value = serde_json::from_slice(manifest_bytes)
        .map_err(|error| manifest_error(manifest_path, format!("JSON no legible: {error}")))?;
    let manifest = manifest
        .as_object()
        .ok_or_else(|| manifest_error(manifest_path, "la raíz no es un objeto"))?;
    if manifest.get("schema_version").and_then(Value::as_str) != Some("1") {
        return Err(effective_error(
            effective_path,
            "el manifiesto no usa la misma versión de contrato que la configuración efectiva",
        ));
    }
    let catalog = artifacts
        .iter()
        .map(|record| (record.id.as_str(), record.kind.as_str()))
        .collect::<BTreeMap<_, _>>();
    let effective_id = optional_text(manifest, "effective_config_artifact_id", manifest_path)?
        .ok_or_else(|| {
            effective_error(
                effective_path,
                "el manifiesto no enlaza la configuración efectiva",
            )
        })?;
    require_artifact_kind(&catalog, effective_id, "effective_config", manifest_path)?;
    let manifest_references = validate_references(
        required_array(manifest, "references", manifest_path)?,
        &catalog,
        manifest_path,
    )?;
    let document: Value = serde_json::from_slice(effective_bytes)
        .map_err(|error| effective_error(effective_path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "la raíz no es un objeto"))?;
    const FIELDS: [&str; 6] = [
        "effective_config_kind",
        "schema_version",
        "source_path",
        "seed",
        "reference",
        "reduced_references",
    ];
    validate_effective_closed_object(object, &FIELDS, "effective-config.json", effective_path)?;
    if object.get("effective_config_kind").and_then(Value::as_str) != Some("gps_effective_config")
        || object.get("schema_version").and_then(Value::as_str) != Some("1")
    {
        return Err(effective_error(
            effective_path,
            "clase o versión incompatible",
        ));
    }
    let source_path = effective_non_empty_text(object, "source_path", effective_path)?;
    let configuration_input = required_object(manifest, "configuration_input", manifest_path)?;
    if required_non_empty_text(configuration_input, "original_path", manifest_path)? != source_path
    {
        return Err(effective_error(
            effective_path,
            "source_path no coincide con la procedencia de la copia TOML",
        ));
    }
    let seed = effective_object(object, "seed", effective_path)?;
    validate_effective_seed(seed, effective_path)?;
    let effective_root_seed = effective_non_empty_text(seed, "root", effective_path)?;
    let effective_reference = effective_object(object, "reference", effective_path)?;
    let reference = validate_effective_reference(effective_reference, &catalog, effective_path)?;
    let eligible_bases =
        effective_positive_u64(effective_reference, "total_genome_size", effective_path)?;
    if manifest_references.get(&reference.0).map(String::as_str) != Some(reference.1.as_str()) {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no coincide con el enlace del manifiesto",
        ));
    }
    let mut effective_references = BTreeMap::from([(reference.0.clone(), reference.1.clone())]);
    let manifest_reductions = required_array(manifest, "reference_reductions", manifest_path)?
        .iter()
        .map(|value| {
            let reduction = value
                .as_object()
                .ok_or_else(|| manifest_error(manifest_path, "una reducción no es objeto"))?;
            let id = required_non_empty_text(reduction, "id", manifest_path)?;
            Ok((id.to_owned(), reduction))
        })
        .collect::<Result<BTreeMap<_, _>, Diagnostic>>()?;
    let effective_reductions = effective_array(object, "reduced_references", effective_path)?;
    if effective_reductions.len() != manifest_reductions.len() {
        return Err(effective_error(
            effective_path,
            "las reducciones efectivas no coinciden con el manifiesto",
        ));
    }
    for value in effective_reductions {
        let reduction = validate_effective_reduction(value, &catalog, effective_path)?;
        let manifest_reduction = manifest_reductions.get(&reduction.id).ok_or_else(|| {
            effective_error(
                effective_path,
                format!("la reducción {} no existe en el manifiesto", reduction.id),
            )
        })?;
        let effective_reduction = value
            .as_object()
            .ok_or_else(|| effective_error(effective_path, "una reducción no es objeto"))?;
        if manifest_reduction.get("state").and_then(Value::as_str) != Some("complete")
            || manifest_reduction
                .get("eligible_bases")
                .and_then(Value::as_u64)
                != Some(eligible_bases)
            || manifest_reduction.get("root_seed").and_then(Value::as_str)
                != Some(effective_root_seed)
            || manifest_reduction
                .get("source_reference_id")
                .and_then(Value::as_str)
                != Some(reduction.source_reference_id.as_str())
            || manifest_reduction
                .get("fasta_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.fasta_artifact_id.as_str())
            || manifest_reduction
                .get("map_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.map_artifact_id.as_str())
            || manifest_reduction
                .get("source_fasta_artifact_id")
                .and_then(Value::as_str)
                != Some(reduction.source_fasta_artifact_id.as_str())
            || [
                "region_size_bases",
                "target",
                "selection_subseed",
                "selection_order",
            ]
            .iter()
            .any(|field| manifest_reduction.get(*field) != effective_reduction.get(*field))
        {
            return Err(effective_error(
                effective_path,
                format!(
                    "la reducción efectiva {} no coincide con el manifiesto",
                    reduction.id
                ),
            ));
        }
        if effective_references
            .insert(reduction.id.clone(), reduction.fasta_artifact_id.clone())
            .is_some()
        {
            return Err(effective_error(
                effective_path,
                "una reducción reutiliza un ID de referencia",
            ));
        }
    }
    if effective_references != manifest_references {
        return Err(effective_error(
            effective_path,
            "el conjunto de referencias del manifiesto no coincide exactamente con la configuración efectiva",
        ));
    }
    Ok(())
}

fn require_truth_artifact(records: &[ArtifactRecord], path: &Path) -> Result<(), Diagnostic> {
    let Some(truth) = records.iter().find(|artifact| artifact.id == "truth-jsonl") else {
        return Err(manifest_error(
            path,
            "falta el artefacto obligatorio truth-jsonl",
        ));
    };
    if truth.kind != "truth_jsonl" || truth.path != "outputs/truth.jsonl" {
        return Err(manifest_error(
            path,
            "truth-jsonl debe tener clase truth_jsonl y ruta outputs/truth.jsonl",
        ));
    }
    Ok(())
}

fn validate_output_categories(
    object: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    expected_minimum: &BTreeSet<&str>,
    expected_requested: &BTreeSet<String>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 4] = ["minimum", "defaulted", "requested", "disabled"];
    validate_closed_object(object, &FIELDS, "output_categories", path)?;
    let mut categories = BTreeMap::new();
    let mut globally_classified = BTreeSet::new();
    for field in FIELDS {
        let mut ids = BTreeSet::new();
        for value in required_array(object, field, path)? {
            let id = value
                .as_str()
                .filter(|id| !id.is_empty())
                .ok_or_else(|| manifest_error(path, format!("{field} contiene un ID inválido")))?;
            if !catalog.contains_key(id) {
                return Err(manifest_error(
                    path,
                    format!("la categoría {field} enlaza un artefacto ausente: {id}"),
                ));
            }
            if !ids.insert(id) {
                return Err(manifest_error(
                    path,
                    format!("el artefacto {id} se repite en {field}"),
                ));
            }
            if !globally_classified.insert(id) {
                return Err(manifest_error(
                    path,
                    format!("el artefacto {id} aparece en más de una categoría"),
                ));
            }
        }
        categories.insert(field, ids);
    }
    if !categories["defaulted"].is_empty() || !categories["disabled"].is_empty() {
        return Err(manifest_error(
            path,
            "defaulted y disabled deben estar vacías en F1",
        ));
    }
    let expected_requested = expected_requested
        .iter()
        .map(String::as_str)
        .collect::<BTreeSet<_>>();
    if !expected_minimum.is_disjoint(&expected_requested) {
        return Err(manifest_error(
            path,
            "un mismo artefacto no puede tener simultáneamente rol minimum y requested",
        ));
    }
    if categories["minimum"] != *expected_minimum || categories["requested"] != expected_requested {
        return Err(manifest_error(
            path,
            "minimum/requested no coinciden con los roles normativos de F1",
        ));
    }
    let expected_catalog = expected_minimum
        .iter()
        .copied()
        .chain(expected_requested.iter().copied())
        .collect::<BTreeSet<_>>();
    if expected_catalog != catalog.keys().copied().collect::<BTreeSet<_>>() {
        return Err(manifest_error(
            path,
            "el catálogo contiene un artefacto sin rol de salida permitido en F1",
        ));
    }
    Ok(())
}

fn validate_manifest_events(
    value: Option<&Value>,
    records: &[ArtifactRecord],
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<Vec<(String, UnitState, String)>, Diagnostic> {
    let Some(value) = value else {
        return Ok(Vec::new());
    };
    let values = value
        .as_array()
        .ok_or_else(|| manifest_error(path, "events no es una lista"))?;
    let mut previous = None;
    let mut links = Vec::with_capacity(values.len());
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "un enlace de evento no es objeto"))?;
        const FIELDS: [&str; 3] = ["event_id", "state", "manifest_artifact_id"];
        validate_closed_object(object, &FIELDS, "events[]", path)?;
        let event_id = required_non_empty_text(object, "event_id", path)?;
        EventId::parse(event_id.to_owned(), "manifest.events")
            .map_err(|error| manifest_error(path, error.render_human()))?;
        if previous.is_some_and(|previous: &str| previous >= event_id) {
            return Err(manifest_error(path, "events no está ordenado por event_id"));
        }
        previous = Some(event_id);
        let state = match required_non_empty_text(object, "state", path)? {
            "complete" => UnitState::Complete,
            "failed" => UnitState::Failed,
            _ => return Err(manifest_error(path, "estado de evento inválido")),
        };
        let artifact_id = required_non_empty_text(object, "manifest_artifact_id", path)?;
        require_artifact_kind(catalog, artifact_id, "event_manifest", path)?;
        let artifact = records
            .iter()
            .find(|artifact| artifact.id == artifact_id)
            .ok_or_else(|| manifest_error(path, "manifiesto de evento ausente"))?;
        let expected_path = event_manifest_relative_path(
            &EventId::parse(event_id.to_owned(), "manifest.events")
                .map_err(|error| manifest_error(path, error.render_human()))?,
        );
        if artifact.path != expected_path.to_string_lossy() {
            return Err(manifest_error(
                path,
                "la ruta del manifiesto de evento no es canónica",
            ));
        }
        links.push((event_id.to_owned(), state, artifact_id.to_owned()));
    }
    Ok(links)
}

fn validate_manifest_technical_events(
    values: &[Value],
    events: &[(String, UnitState, String)],
    effective_id: Option<&str>,
    records: &[ArtifactRecord],
    path: &Path,
) -> Result<(), Diagnostic> {
    if values.len() != events.len() {
        return Err(manifest_error(
            path,
            "technical_events no coincide con events",
        ));
    }
    for (value, (event_id, _, manifest_id)) in values.iter().zip(events) {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "technical_event inválido"))?;
        const FIELDS: [&str; 7] = [
            "event_id",
            "effective_config",
            "profile",
            "rng_catalog",
            "id_catalog",
            "artifacts",
            "event_manifest",
        ];
        validate_closed_object(object, &FIELDS, "technical_events[]", path)?;
        if required_non_empty_text(object, "event_id", path)? != event_id
            || !matches!(
                object.get("profile"),
                Some(profile) if *profile == json!({
                    "id": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_V1,
                    "sha256": crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_SHA256,
                }) || *profile == json!({
                    "id": crate::technical_profile::GPS_HISEQX_PCR_PE150_V1,
                    "sha256": crate::technical_profile::GPS_HISEQX_PCR_PE150_SHA256,
                })
            )
            || !matches!(
                object.get("rng_catalog").and_then(Value::as_str),
                Some(
                    crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
                        | crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
                        | crate::randomness::F4_TECHNICAL_RNG_CATALOG_V5
                )
            )
            || object.get("id_catalog") != Some(&f4_id_catalog_value())
        {
            return Err(manifest_error(path, "inventario técnico F4 incompatible"));
        }
        validate_embedded_artifact(object.get("effective_config"), effective_id, records, path)?;
        validate_embedded_artifact(
            object.get("event_manifest"),
            Some(manifest_id),
            records,
            path,
        )?;
        let artifacts = object
            .get("artifacts")
            .and_then(Value::as_object)
            .ok_or_else(|| manifest_error(path, "artefactos técnicos agregados inválidos"))?;
        validate_closed_object(
            artifacts,
            &["chimera"],
            "technical_events[].artifacts",
            path,
        )?;
        let chimera = artifacts
            .get("chimera")
            .and_then(Value::as_object)
            .ok_or_else(|| manifest_error(path, "resumen de quimera inválido"))?;
        validate_closed_object(
            chimera,
            &["rate", "realized_pairs"],
            "technical_events[].artifacts.chimera",
            path,
        )?;
        let rate = chimera.get("rate").and_then(Value::as_f64);
        let realized_pairs = chimera.get("realized_pairs").and_then(Value::as_u64);
        if !matches!(rate, Some(rate) if rate.is_finite() && (0.0..=1.0).contains(&rate))
            || realized_pairs.is_none()
            || (rate == Some(0.0) && realized_pairs != Some(0))
        {
            return Err(manifest_error(path, "resumen de quimera inválido"));
        }
    }
    Ok(())
}

fn validate_embedded_artifact(
    value: Option<&Value>,
    expected_id: Option<&str>,
    records: &[ArtifactRecord],
    path: &Path,
) -> Result<(), Diagnostic> {
    match (value, expected_id) {
        (Some(Value::Null), None) => Ok(()),
        (Some(value), Some(expected_id))
            if value
                == &records
                    .iter()
                    .find(|record| record.id == expected_id)
                    .map_or(Value::Null, ArtifactRecord::to_json) =>
        {
            Ok(())
        }
        _ => Err(manifest_error(
            path,
            "artefacto técnico embebido no coincide con el catálogo",
        )),
    }
}

fn validate_manifest_population_snvs(
    value: Option<&Value>,
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<Option<String>, Diagnostic> {
    let Some(value) = value else { return Ok(None) };
    if value.is_null() {
        return Ok(None);
    }
    let object = value
        .as_object()
        .ok_or_else(|| manifest_error(path, "population_snvs inválido"))?;
    const FIELDS: [&str; 5] = [
        "mode",
        "metadata_sha256",
        "metadata_artifact_id",
        "normalization_audit",
        "initial_population_statistics",
    ];
    let fields = if object.contains_key("initial_population_statistics") {
        &FIELDS[..]
    } else {
        &FIELDS[..4]
    };
    validate_closed_object(object, fields, "population_snvs", path)?;
    let mode = match object.get("mode").and_then(Value::as_str) {
        Some(mode @ ("supplied" | "generated" | "reuse")) => mode,
        _ => return Err(manifest_error(path, "population_snvs.mode inválido")),
    };
    let digest = required_non_empty_text(object, "metadata_sha256", path)?;
    if digest.len() != 64
        || !digest
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        return Err(manifest_error(
            path,
            "population_snvs.metadata_sha256 inválido",
        ));
    }
    validate_serialized_normalization_audit(
        object.get("normalization_audit"),
        "population_snvs.normalization_audit inválida",
    )
    .map_err(|detail| manifest_error_from_reason(path, detail))?;
    if let Some(statistics) = object.get("initial_population_statistics") {
        validate_initial_population_statistics(statistics)
            .map_err(|detail| manifest_error_from_reason(path, detail))?;
    }
    match (mode, object.get("metadata_artifact_id")) {
        ("reuse", None | Some(Value::Null)) => Ok(None),
        ("supplied" | "generated", Some(Value::String(id))) if !id.is_empty() => {
            require_artifact_kind(catalog, id, "population_snvs_metadata", path)?;
            Ok(Some(id.clone()))
        }
        _ => Err(manifest_error(
            path,
            "population_snvs.metadata_artifact_id inválido",
        )),
    }
}

fn validate_initial_population_statistics(value: &Value) -> Result<(), &'static str> {
    let object = value
        .as_object()
        .ok_or_else(|| reason("estadísticas iniciales inválidas"))?;
    const FIELDS: [&str; 4] = ["kind", "version", "provenance", "populations"];
    if object.len() != FIELDS.len() || object.keys().any(|field| !FIELDS.contains(&field.as_str()))
    {
        return Err(reason("estadísticas iniciales inválidas"));
    }
    if object.get("kind").and_then(Value::as_str) != Some("initial_population_statistics")
        || object.get("version").and_then(Value::as_str) != Some("1")
    {
        return Err(reason("estadísticas iniciales inválidas"));
    }
    let provenance = object
        .get("provenance")
        .and_then(Value::as_object)
        .ok_or_else(|| reason("procedencia de estadísticas iniciales inválida"))?;
    validate_initial_population_statistics_provenance(provenance)?;
    let populations = object
        .get("populations")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("estadísticas iniciales inválidas"))?;
    if populations.is_empty() {
        return Err(reason("estadísticas iniciales sin poblaciones"));
    }
    let mut population_ids = BTreeSet::new();
    for population in populations {
        validate_initial_population_statistics_population(population)?;
        let population_id = population
            .get("population_id")
            .and_then(Value::as_str)
            .ok_or_else(|| reason("población estadística inválida"))?;
        if !population_ids.insert(population_id) {
            return Err(reason("población estadística repetida"));
        }
    }
    Ok(())
}

fn validate_initial_population_statistics_provenance(
    provenance: &serde_json::Map<String, Value>,
) -> Result<(), &'static str> {
    const FIELDS: [&str; 6] = [
        "matrix_reference_id",
        "matrix_reference_fingerprint",
        "matrix_canonical_tsv_sha256",
        "genotype_algorithms",
        "missing_locus_policy",
        "explicit_snp_constraints",
    ];
    if provenance.len() != FIELDS.len()
        || provenance
            .keys()
            .any(|field| !FIELDS.contains(&field.as_str()))
    {
        return Err(reason("procedencia de estadísticas iniciales inválida"));
    }
    let reference_id = provenance
        .get("matrix_reference_id")
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty());
    if reference_id.is_none()
        || !valid_lower_sha256(provenance.get("matrix_reference_fingerprint"))
        || !valid_lower_sha256(provenance.get("matrix_canonical_tsv_sha256"))
        || provenance
            .get("missing_locus_policy")
            .and_then(Value::as_str)
            != Some("missing_by_deletion_excluded_from_allele_frequency_denominator_v1")
    {
        return Err(reason("procedencia de estadísticas iniciales inválida"));
    }
    let algorithms = provenance
        .get("genotype_algorithms")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("algoritmos de genotipo inválidos"))?;
    if algorithms.len() != 3
        || algorithms[0].as_str() != Some("hwe-composition-v1")
        || algorithms[1].as_str() != Some("f-diploid-genotype-v1")
        || algorithms[2].as_str() != Some("phase-permutation-v1")
    {
        return Err(reason("algoritmos de genotipo inválidos"));
    }
    let constraints = provenance
        .get("explicit_snp_constraints")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("restricciones SNV explícitas inválidas"))?;
    let mut constraint_ids = BTreeSet::new();
    for constraint in constraints {
        let constraint = constraint
            .as_object()
            .ok_or_else(|| reason("restricción SNV explícita inválida"))?;
        const CONSTRAINT_FIELDS: [&str; 8] = [
            "population_id",
            "individual_id",
            "event_id",
            "phase",
            "matrix_variant_id",
            "sampled_alternate",
            "explicit_alternate",
            "incompatible",
        ];
        if constraint.len() != CONSTRAINT_FIELDS.len()
            || constraint
                .keys()
                .any(|field| !CONSTRAINT_FIELDS.contains(&field.as_str()))
            || ![
                "population_id",
                "individual_id",
                "event_id",
                "matrix_variant_id",
            ]
            .iter()
            .all(|field| {
                constraint
                    .get(*field)
                    .and_then(Value::as_str)
                    .is_some_and(|value| !value.is_empty())
            })
            || !matches!(constraint.get("phase").and_then(Value::as_u64), Some(0..=9))
            || !constraint.get("sampled_alternate").is_some_and(|value| {
                value.is_null() || value.as_str().is_some_and(is_canonical_base_text)
            })
            || !constraint
                .get("explicit_alternate")
                .and_then(Value::as_str)
                .is_some_and(is_canonical_base_text)
            || !constraint
                .get("incompatible")
                .is_some_and(Value::is_boolean)
        {
            return Err(reason("restricción SNV explícita inválida"));
        }
        let key = (
            constraint["population_id"].as_str().unwrap_or_default(),
            constraint["individual_id"].as_str().unwrap_or_default(),
            constraint["event_id"].as_str().unwrap_or_default(),
        );
        if !constraint_ids.insert(key) {
            return Err(reason("restricción SNV explícita repetida"));
        }
    }
    Ok(())
}

fn validate_initial_population_statistics_population(value: &Value) -> Result<(), &'static str> {
    let population = value
        .as_object()
        .ok_or_else(|| reason("población estadística inválida"))?;
    const FIELDS: [&str; 5] = [
        "population_id",
        "genotype_model",
        "individual_count",
        "ploidy",
        "loci",
    ];
    if population.len() != FIELDS.len()
        || population.keys().any(|key| !FIELDS.contains(&key.as_str()))
    {
        return Err(reason("población estadística inválida"));
    }
    let population_id = population
        .get("population_id")
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty());
    let individual_count = population.get("individual_count").and_then(Value::as_u64);
    let ploidy = population.get("ploidy").and_then(Value::as_u64);
    if population_id.is_none()
        || individual_count == Some(0)
        || individual_count.is_none()
        || !matches!(ploidy, Some(1..=10))
    {
        return Err(reason("población estadística inválida"));
    }
    validate_initial_population_statistics_model(
        population
            .get("genotype_model")
            .and_then(Value::as_object)
            .ok_or_else(|| reason("modelo genotípico estadístico inválido"))?,
        ploidy.unwrap_or_default(),
    )?;
    let loci = population
        .get("loci")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("loci estadísticos inválidos"))?;
    let mut locus_keys = BTreeSet::new();
    for locus in loci {
        validate_initial_population_statistics_locus(
            locus,
            individual_count.unwrap_or_default(),
            ploidy.unwrap_or_default(),
        )?;
        let key = (
            locus
                .get("contig")
                .and_then(Value::as_str)
                .unwrap_or_default(),
            locus
                .get("position_1based")
                .and_then(Value::as_u64)
                .unwrap_or_default(),
            locus.get("ref").and_then(Value::as_str).unwrap_or_default(),
        );
        if !locus_keys.insert(key) {
            return Err(reason("locus estadístico repetido"));
        }
    }
    Ok(())
}

fn validate_initial_population_statistics_model(
    model: &serde_json::Map<String, Value>,
    ploidy: u64,
) -> Result<(), &'static str> {
    let kind = model.get("kind").and_then(Value::as_str);
    match kind {
        Some("hwe") if model.len() == 1 => Ok(()),
        Some("f_diploid")
            if ploidy == 2
                && model.len() == 2
                && model
                    .get("inbreeding_coefficient")
                    .and_then(Value::as_f64)
                    .is_some_and(|value| value.is_finite() && (0.0..=1.0).contains(&value)) =>
        {
            Ok(())
        }
        _ => Err(reason("modelo genotípico estadístico inválido")),
    }
}

fn validate_initial_population_statistics_locus(
    value: &Value,
    individual_count: u64,
    ploidy: u64,
) -> Result<(), &'static str> {
    let locus = value
        .as_object()
        .ok_or_else(|| reason("locus estadístico inválido"))?;
    const FIELDS: [&str; 9] = [
        "contig",
        "position_1based",
        "ref",
        "total_copies",
        "observed_copies",
        "missing_by_deletion",
        "alleles",
        "complete_genotype_observations",
        "genotypes",
    ];
    if locus.len() != FIELDS.len() || locus.keys().any(|key| !FIELDS.contains(&key.as_str())) {
        return Err(reason("locus estadístico inválido"));
    }
    if !locus
        .get("contig")
        .and_then(Value::as_str)
        .is_some_and(|value| !value.is_empty())
        || !matches!(
            locus.get("position_1based").and_then(Value::as_u64),
            Some(1..)
        )
        || !locus
            .get("ref")
            .and_then(Value::as_str)
            .is_some_and(is_canonical_base_text)
    {
        return Err(reason("locus estadístico inválido"));
    }
    let total = locus
        .get("total_copies")
        .and_then(Value::as_u64)
        .ok_or_else(|| reason("total de copias inválido"))?;
    let observed = locus
        .get("observed_copies")
        .and_then(Value::as_u64)
        .ok_or_else(|| reason("copias observadas inválidas"))?;
    let missing = locus
        .get("missing_by_deletion")
        .and_then(Value::as_u64)
        .ok_or_else(|| reason("copias ausentes inválidas"))?;
    if total
        != individual_count
            .checked_mul(ploidy)
            .ok_or_else(|| reason("total de copias inválido"))?
        || observed.checked_add(missing) != Some(total)
    {
        return Err(reason("las copias observadas y ausentes no suman el total"));
    }
    let alleles = locus
        .get("alleles")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("alelos estadísticos inválidos"))?;
    if alleles.len() < 2 {
        return Err(reason("alelos estadísticos inválidos"));
    }
    let reference = locus.get("ref").and_then(Value::as_str).unwrap_or_default();
    let mut allele_texts = BTreeSet::new();
    let mut count_sum = 0_u64;
    let mut requested_sum = 0.0;
    for (index, allele) in alleles.iter().enumerate() {
        let allele = allele
            .as_object()
            .ok_or_else(|| reason("alelo estadístico inválido"))?;
        const REF_FIELDS: [&str; 5] = [
            "allele",
            "kind",
            "requested_frequency",
            "realized_count",
            "realized_frequency",
        ];
        const ALT_FIELDS: [&str; 7] = [
            "allele",
            "kind",
            "variant_id",
            "source_id",
            "requested_frequency",
            "realized_count",
            "realized_frequency",
        ];
        let kind = allele.get("kind").and_then(Value::as_str);
        let fields = match kind {
            Some("ref") if index == 0 => &REF_FIELDS[..],
            Some("alt") if index > 0 => &ALT_FIELDS[..],
            _ => return Err(reason("clase u orden de alelo estadístico inválido")),
        };
        if allele.len() != fields.len()
            || allele.keys().any(|field| !fields.contains(&field.as_str()))
        {
            return Err(reason("alelo estadístico inválido"));
        }
        let allele_text = allele
            .get("allele")
            .and_then(Value::as_str)
            .filter(|value| is_canonical_base_text(value))
            .ok_or_else(|| reason("alelo estadístico inválido"))?;
        if !allele_texts.insert(allele_text) || (index == 0 && allele_text != reference) {
            return Err(reason(
                "alelo estadístico duplicado o referencia incoherente",
            ));
        }
        if index > 0
            && (!allele
                .get("variant_id")
                .and_then(Value::as_str)
                .is_some_and(|value| !value.is_empty())
                || allele.get("source_id").and_then(Value::as_str).is_none())
        {
            return Err(reason("procedencia de alelo alternativo inválida"));
        }
        let count = allele
            .get("realized_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| reason("recuento alélico inválido"))?;
        let requested = allele
            .get("requested_frequency")
            .and_then(Value::as_f64)
            .ok_or_else(|| reason("frecuencia objetivo inválida"))?;
        let realized = allele
            .get("realized_frequency")
            .and_then(Value::as_f64)
            .ok_or_else(|| reason("frecuencia realizada inválida"))?;
        if !requested.is_finite()
            || !realized.is_finite()
            || !(0.0..=1.0).contains(&requested)
            || !(0.0..=1.0).contains(&realized)
        {
            return Err(reason("frecuencia estadística no finita o fuera de rango"));
        }
        let expected_realized = if observed == 0 {
            0.0
        } else {
            count as f64 / observed as f64
        };
        if realized.to_bits() != expected_realized.to_bits() {
            return Err(reason(
                "frecuencia realizada no coincide con recuento y denominador",
            ));
        }
        requested_sum += requested;
        count_sum = count_sum
            .checked_add(count)
            .ok_or_else(|| reason("recuento alélico desbordado"))?;
    }
    if count_sum != observed {
        return Err(reason(
            "los recuentos alélicos no suman las copias observadas",
        ));
    }
    if !requested_sum.is_finite() || (requested_sum - 1.0).abs() > 1e-12 {
        return Err(reason("las frecuencias objetivo no suman uno"));
    }
    let complete = locus
        .get("complete_genotype_observations")
        .and_then(Value::as_u64)
        .ok_or_else(|| reason("observaciones genotípicas inválidas"))?;
    if complete > individual_count {
        return Err(reason("observaciones genotípicas fuera de rango"));
    }
    let genotypes = locus
        .get("genotypes")
        .and_then(Value::as_array)
        .ok_or_else(|| reason("genotipos estadísticos inválidos"))?;
    if genotypes.is_empty() {
        return Err(reason("genotipos estadísticos inválidos"));
    }
    let mut genotype_keys = BTreeSet::new();
    let mut expected_probability_sum = 0.0;
    let mut realized_genotype_sum = 0_u64;
    for genotype in genotypes {
        let genotype = genotype
            .as_object()
            .ok_or_else(|| reason("genotipo estadístico inválido"))?;
        const GENOTYPE_FIELDS: [&str; 4] = [
            "allele_dosages",
            "expected_probability",
            "expected_count",
            "realized_count",
        ];
        if genotype.len() != GENOTYPE_FIELDS.len()
            || genotype
                .keys()
                .any(|field| !GENOTYPE_FIELDS.contains(&field.as_str()))
        {
            return Err(reason("genotipo estadístico inválido"));
        }
        let dosages = genotype
            .get("allele_dosages")
            .and_then(Value::as_array)
            .ok_or_else(|| reason("dosis genotípicas inválidas"))?;
        let dosages = dosages
            .iter()
            .map(Value::as_u64)
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| reason("dosis genotípicas inválidas"))?;
        let dosage_sum = dosages.iter().try_fold(0_u64, |sum, dosage| {
            sum.checked_add(*dosage)
                .ok_or_else(|| reason("dosis genotípicas desbordadas"))
        })?;
        if dosages.len() != alleles.len() || dosage_sum != ploidy || !genotype_keys.insert(dosages)
        {
            return Err(reason("dosis genotípicas inválidas o repetidas"));
        }
        let probability = genotype
            .get("expected_probability")
            .and_then(Value::as_f64)
            .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
            .ok_or_else(|| reason("probabilidad genotípica inválida"))?;
        genotype
            .get("expected_count")
            .and_then(Value::as_f64)
            .filter(|value| value.is_finite() && *value >= 0.0)
            .ok_or_else(|| reason("recuento genotípico esperado inválido"))?;
        // El recuento esperado se conserva como una magnitud de salida. Su
        // correspondencia exacta con la probabilidad y la población se vuelve
        // a calcular desde TOML, semilla y matriz en la puerta física; aquí se
        // valida sólo su representación cerrada y finita para no confundir el
        // contrato de esquema con ese oráculo independiente.
        let realized = genotype
            .get("realized_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| reason("recuento genotípico realizado inválido"))?;
        expected_probability_sum += probability;
        realized_genotype_sum = realized_genotype_sum
            .checked_add(realized)
            .ok_or_else(|| reason("recuento genotípico realizado desbordado"))?;
    }
    if !expected_probability_sum.is_finite() || (expected_probability_sum - 1.0).abs() > 1e-12 {
        return Err(reason("probabilidades genotípicas no suman uno"));
    }
    if realized_genotype_sum != complete {
        return Err(reason("recuentos genotípicos realizados incoherentes"));
    }
    Ok(())
}

fn valid_lower_sha256(value: Option<&Value>) -> bool {
    value.and_then(Value::as_str).is_some_and(|digest| {
        digest.len() == 64
            && digest
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    })
}

fn is_canonical_base_text(value: &str) -> bool {
    matches!(value, "A" | "C" | "G" | "T")
}

fn validate_serialized_normalization_audit(
    value: Option<&Value>,
    detail: &'static str,
) -> Result<(), &'static str> {
    let groups = value
        .and_then(Value::as_array)
        .ok_or_else(|| reason(detail))?;
    let mut previous_group = None;
    for group in groups {
        let group = group.as_object().ok_or_else(|| reason(detail))?;
        const GROUP_FIELDS: [&str; 14] = [
            "population_id",
            "contig",
            "position",
            "ref",
            "tolerance",
            "canonical_order",
            "sum_algorithm",
            "original_sum",
            "original_delta",
            "effective_sum",
            "effective_delta",
            "reference_frequency",
            "decision",
            "values",
        ];
        if group.len() != GROUP_FIELDS.len()
            || group
                .keys()
                .any(|field| !GROUP_FIELDS.contains(&field.as_str()))
        {
            return Err(reason(detail));
        }
        let population = group
            .get("population_id")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .ok_or_else(|| reason(detail))?;
        let contig = group
            .get("contig")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .ok_or_else(|| reason(detail))?;
        let position = group
            .get("position")
            .and_then(Value::as_str)
            .and_then(|value| value.parse::<u64>().ok())
            .filter(|value| *value > 0)
            .ok_or_else(|| reason(detail))?;
        let reference = group
            .get("ref")
            .and_then(Value::as_str)
            .and_then(single_canonical_base)
            .ok_or_else(|| reason(detail))?;
        if group
            .get("tolerance")
            .and_then(Value::as_f64)
            .filter(|value| value.is_finite() && (0.0..=1e-6).contains(value))
            .is_none()
            || group.get("canonical_order").and_then(Value::as_str)
                != Some("population_id,contig,position,ref,alt,variant_id,source_id")
            || group.get("sum_algorithm").and_then(Value::as_str) != Some("neumaier_sum_v1")
            || [
                "original_sum",
                "original_delta",
                "effective_sum",
                "effective_delta",
            ]
            .iter()
            .any(|field| {
                group
                    .get(*field)
                    .and_then(Value::as_f64)
                    .filter(|value| value.is_finite() && *value >= 0.0)
                    .is_none()
            })
            || group
                .get("reference_frequency")
                .and_then(Value::as_f64)
                .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
                .is_none()
            || !matches!(
                group.get("decision").and_then(Value::as_str),
                Some("reference_remainder" | "normalized")
            )
        {
            return Err(reason(detail));
        }
        let group_key = (population, contig, position, reference);
        if previous_group
            .as_ref()
            .is_some_and(|last| last >= &group_key)
        {
            return Err(reason(detail));
        }
        previous_group = Some(group_key);
        let values = group
            .get("values")
            .and_then(Value::as_array)
            .filter(|values| !values.is_empty())
            .ok_or_else(|| reason(detail))?;
        let mut previous_value = None;
        for value in values {
            let value = value.as_object().ok_or_else(|| reason(detail))?;
            const VALUE_FIELDS: [&str; 5] =
                ["alt", "variant_id", "source_id", "original", "effective"];
            if value.len() != VALUE_FIELDS.len()
                || value
                    .keys()
                    .any(|field| !VALUE_FIELDS.contains(&field.as_str()))
            {
                return Err(reason(detail));
            }
            let alternate = value
                .get("alt")
                .and_then(Value::as_str)
                .and_then(single_canonical_base)
                .filter(|alternate| *alternate != reference)
                .ok_or_else(|| reason(detail))?;
            let variant_id = value
                .get("variant_id")
                .and_then(Value::as_str)
                .filter(|value| !value.is_empty())
                .ok_or_else(|| reason(detail))?;
            let source_id = value
                .get("source_id")
                .and_then(Value::as_str)
                .ok_or_else(|| reason(detail))?;
            if ["original", "effective"].iter().any(|field| {
                value
                    .get(*field)
                    .and_then(Value::as_f64)
                    .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
                    .is_none()
            }) {
                return Err(reason(detail));
            }
            let value_key = (alternate, variant_id, source_id);
            if previous_value
                .as_ref()
                .is_some_and(|last| last >= &value_key)
            {
                return Err(reason(detail));
            }
            previous_value = Some(value_key);
        }
    }
    Ok(())
}

fn single_canonical_base(value: &str) -> Option<u8> {
    (value.len() == 1)
        .then(|| value.as_bytes()[0])
        .filter(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

fn validate_effective_population_snvs(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let effective_value = effective.get("population_snvs").unwrap_or(&Value::Null);
    let manifest_value = manifest.get("population_snvs").unwrap_or(&Value::Null);
    if effective_value != manifest_value {
        return Err(effective_error(
            effective_path,
            "population_snvs no coincide con el manifiesto",
        ));
    }
    if effective_value.is_null() {
        return Ok(());
    }
    let object = effective_value
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "population_snvs inválido"))?;
    let id = validate_manifest_population_snvs(Some(effective_value), catalog, manifest_path)?;
    if object.get("mode").and_then(Value::as_str) != Some("reuse") && id.is_none() {
        return Err(effective_error(
            effective_path,
            "population_snvs publicado requiere metadata_artifact_id",
        ));
    }
    Ok(())
}

fn validate_effective_authorizations(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let effective_value = effective
        .get("authorizations")
        .ok_or_else(|| effective_error(effective_path, "falta authorizations"))?;
    let manifest_value = manifest
        .get("authorizations")
        .ok_or_else(|| manifest_error(manifest_path, "falta authorizations"))?;
    if effective_value != manifest_value {
        return Err(effective_error(
            effective_path,
            "authorizations no coincide con el manifiesto",
        ));
    }
    validate_authorizations_value(effective_value)
        .map_err(|detail| effective_error_from_reason(effective_path, detail))
}

fn validate_effective_vcf_output(
    effective: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<Option<VcfOutputMode>, Diagnostic> {
    let Some(value) = effective.get("vcf_output") else {
        // Los documentos F1/F2 previos a F2.5.5 no incluyen la clave. La
        // reconstrucción física contra el TOML actual sí exige `null` explícito.
        return Ok(None);
    };
    if value.is_null() {
        return Ok(None);
    }
    let output = value
        .as_object()
        .ok_or_else(|| effective_error(path, "vcf_output inválido"))?;
    // La versión 6 registra la del algoritmo de normalización (DTC-064), y
    // solo admite la que este programa aplica: un documento que declare otra
    // describe un VCF que no es el que se publicó.
    if effective.get("schema_version").and_then(Value::as_str)
        == Some(EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS)
    {
        const FIELDS: [&str; 2] = ["mode", "normalization_version"];
        validate_effective_closed_object(output, &FIELDS, "vcf_output", path)?;
        if output.get("normalization_version").and_then(Value::as_str)
            != Some(VCF_NORMALIZATION_VERSION)
        {
            return Err(effective_error(
                path,
                "vcf_output.normalization_version no es la versión del algoritmo de normalización VCF",
            ));
        }
    } else {
        const FIELDS: [&str; 1] = ["mode"];
        validate_effective_closed_object(output, &FIELDS, "vcf_output", path)?;
    }
    match output.get("mode").and_then(Value::as_str) {
        Some("text") => Ok(Some(VcfOutputMode::Text)),
        Some("bgzf") => Ok(Some(VcfOutputMode::Bgzf)),
        _ => Err(effective_error(path, "vcf_output.mode inválido")),
    }
}

fn validate_effective_vcf_against_manifest(
    effective: &serde_json::Map<String, Value>,
    manifest: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    reference: &serde_json::Map<String, Value>,
    effective_path: &Path,
    manifest_path: &Path,
) -> Result<(), Diagnostic> {
    let mode = validate_effective_vcf_output(effective, effective_path)?;
    let outputs = required_array(manifest, "vcf_outputs", manifest_path)?;
    let Some(mode) = mode else {
        if !outputs.is_empty() {
            return Err(effective_error(
                effective_path,
                "el manifiesto declara VCF aunque la configuración efectiva lo desactiva",
            ));
        }
        return Ok(());
    };
    if outputs.len() != 1 {
        return Err(effective_error(
            effective_path,
            "la configuración efectiva VCF exige exactamente un enlace en el manifiesto",
        ));
    }
    let output = outputs[0]
        .as_object()
        .ok_or_else(|| effective_error(effective_path, "enlace VCF no es objeto"))?;
    let effective_reference_id = effective_non_empty_text(reference, "id", effective_path)?;
    let expected_mode = match mode {
        VcfOutputMode::Text => "text",
        VcfOutputMode::Bgzf => "bgzf",
    };
    if output.get("reference_id").and_then(Value::as_str) != Some(effective_reference_id)
        || output.get("mode").and_then(Value::as_str) != Some(expected_mode)
        || output.get("vcf_artifact_id").and_then(Value::as_str) != Some(VCF_ARTIFACT_ID)
        || catalog.get(VCF_ARTIFACT_ID) != Some(&VCF_ARTIFACT_KIND)
    {
        return Err(effective_error(
            effective_path,
            "referencia, modo o artefacto VCF del manifiesto no coincide con la configuración efectiva",
        ));
    }
    match mode {
        VcfOutputMode::Text => {
            if !output.get("index_kind").is_some_and(Value::is_null)
                || !output.get("index_artifact_id").is_some_and(Value::is_null)
            {
                return Err(effective_error(
                    effective_path,
                    "VCF de texto efectivo no admite índice en el manifiesto",
                ));
            }
        }
        VcfOutputMode::Bgzf => {
            let expected_index = expected_vcf_index_kind(reference, effective_path)?;
            let (index_kind, artifact_kind) = match expected_index {
                VcfIndexKind::Tbi => ("tbi", "vcf_tbi_index"),
                VcfIndexKind::Csi => ("csi", "vcf_csi_index"),
            };
            if output.get("index_kind").and_then(Value::as_str) != Some(index_kind)
                || output.get("index_artifact_id").and_then(Value::as_str)
                    != Some(VCF_INDEX_ARTIFACT_ID)
                || catalog.get(VCF_INDEX_ARTIFACT_ID) != Some(&artifact_kind)
            {
                return Err(effective_error(
                    effective_path,
                    "índice VCF del manifiesto no coincide con referencia y configuración efectiva",
                ));
            }
        }
    }
    Ok(())
}

fn expected_vcf_index_kind(
    reference: &serde_json::Map<String, Value>,
    effective_path: &Path,
) -> Result<VcfIndexKind, Diagnostic> {
    let sequences = effective_array(reference, "sequences", effective_path)?;
    if sequences.is_empty() {
        return Err(effective_error(
            effective_path,
            "la referencia efectiva no tiene secuencias para seleccionar índice VCF",
        ));
    }
    if sequences.iter().any(|sequence| {
        sequence
            .as_object()
            .and_then(|sequence| sequence.get("length"))
            .and_then(Value::as_u64)
            .is_none_or(|length| length >= TBI_MAXIMUM_POSITION_EXCLUSIVE)
    }) {
        Ok(VcfIndexKind::Csi)
    } else {
        Ok(VcfIndexKind::Tbi)
    }
}

fn validate_effective_initial_populations(
    effective: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<BTreeSet<IndividualId>, Diagnostic> {
    let populations = effective_array(effective, "initial_populations", path)?;
    let mut population_ids = BTreeSet::new();
    let mut individual_ids = BTreeSet::new();
    for value in populations {
        let population = value
            .as_object()
            .ok_or_else(|| effective_error(path, "initial_populations[] inválida"))?;
        const FIELDS: [&str; 5] = [
            "population_id",
            "individual_count",
            "individual_ids",
            "ploidy",
            "genotype_model",
        ];
        const DERIVED_FIELDS: [&str; 6] = [
            "population_id",
            "individual_count",
            "individual_ids",
            "ploidy",
            "genotype_model",
            "founder_keys",
        ];
        let fields = if population.contains_key("founder_keys") {
            &DERIVED_FIELDS[..]
        } else {
            &FIELDS[..]
        };
        validate_effective_closed_object(population, fields, "initial_populations[]", path)?;
        let population_id = effective_non_empty_text(population, "population_id", path)?;
        PopulationId::parse(population_id, "effective_config")
            .map_err(|_| effective_error(path, "population_id inicial inválido"))?;
        if !population_ids.insert(population_id) {
            return Err(effective_error(path, "population_id inicial repetido"));
        }
        let count = effective_positive_u64(population, "individual_count", path)?;
        let ploidy = effective_positive_u64(population, "ploidy", path)?;
        if ploidy > 10 {
            return Err(effective_error(path, "ploidy inicial fuera de rango"));
        }
        let ids = effective_array(population, "individual_ids", path)?;
        if usize::try_from(count).ok() != Some(ids.len()) {
            return Err(effective_error(
                path,
                "individual_count inicial no coincide con individual_ids",
            ));
        }
        for id in ids {
            let id = id
                .as_str()
                .filter(|id| !id.is_empty())
                .ok_or_else(|| effective_error(path, "individual_id inicial inválido"))?;
            let individual_id = IndividualId::parse(id, "effective_config")
                .map_err(|_| effective_error(path, "individual_id inicial inválido"))?;
            if !individual_ids.insert(individual_id) {
                return Err(effective_error(path, "individual_id inicial repetido"));
            }
        }
        if let Some(founder_keys_value) = population.get("founder_keys") {
            let mut founder_keys = BTreeSet::new();
            let founder_key_values = founder_keys_value
                .as_array()
                .ok_or_else(|| effective_error(path, "founder_keys inicial inválida"))?;
            if founder_key_values.len() != ids.len() {
                return Err(effective_error(
                    path,
                    "founder_keys inicial no coincide con individual_ids",
                ));
            }
            for key in founder_key_values {
                let key = key
                    .as_str()
                    .ok_or_else(|| effective_error(path, "founder_key inicial inválida"))?;
                FounderKey::parse(key, "effective_config")
                    .map_err(|_| effective_error(path, "founder_key inicial inválida"))?;
                if !founder_keys.insert(key) {
                    return Err(effective_error(path, "founder_key inicial repetida"));
                }
            }
        }
        let model = effective_object(population, "genotype_model", path)?;
        match model.get("kind").and_then(Value::as_str) {
            Some("hwe") if model.len() == 1 => {}
            Some("f_diploid") if model.len() == 2 && ploidy == 2 => {
                let coefficient = model
                    .get("inbreeding_coefficient")
                    .and_then(Value::as_f64)
                    .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
                    .ok_or_else(|| {
                        effective_error(path, "inbreeding_coefficient inicial inválido")
                    })?;
                let _ = coefficient;
            }
            _ => return Err(effective_error(path, "genotype_model inicial inválido")),
        }
    }
    if !populations.is_empty()
        && !effective
            .get("population_snvs")
            .and_then(Value::as_object)
            .is_some_and(|link| link.contains_key("initial_population_statistics"))
    {
        return Err(effective_error(
            path,
            "faltan estadísticas iniciales para cohortes iniciales",
        ));
    }
    Ok(individual_ids)
}

fn validate_effective_direct_events(
    effective: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<BTreeSet<EventId>, Diagnostic> {
    let events = effective_array(effective, "direct_events", path)?;
    let mut event_ids = BTreeSet::new();
    for value in events {
        let event = value
            .as_object()
            .ok_or_else(|| effective_error(path, "direct_events[] inválido"))?;
        let kind = effective_non_empty_text(event, "kind", path)?;
        let fields: &[&str] = match kind {
            "snp" => &[
                "event_id",
                "individual_id",
                "phase",
                "kind",
                "contig",
                "position_1based",
                "ref",
                "alt",
            ],
            "insertion" => &[
                "event_id",
                "individual_id",
                "phase",
                "kind",
                "contig",
                "point",
                "bases",
            ],
            "deletion" => &[
                "event_id",
                "individual_id",
                "phase",
                "kind",
                "contig",
                "start",
                "end",
            ],
            _ => return Err(effective_error(path, "direct_events[].kind inválido")),
        };
        validate_effective_closed_object(event, fields, "direct_events[]", path)?;
        let event_id = effective_non_empty_text(event, "event_id", path)?;
        let event_id = EventId::parse(event_id, "effective_config")
            .map_err(|_| effective_error(path, "event_id directo inválido"))?;
        if !event_ids.insert(event_id) {
            return Err(effective_error(path, "event_id directo repetido"));
        }
        IndividualId::parse(
            effective_non_empty_text(event, "individual_id", path)?,
            "effective_config",
        )
        .map_err(|_| effective_error(path, "individual_id directo inválido"))?;
        SequenceId::parse(
            effective_non_empty_text(event, "contig", path)?,
            "effective_config",
        )
        .map_err(|_| effective_error(path, "contig directo inválido"))?;
        if effective_nonnegative_u64(event, "phase", path)? > 9 {
            return Err(effective_error(path, "fase directa fuera de rango"));
        }
        match kind {
            "snp" => {
                let _ = effective_positive_u64(event, "position_1based", path)?;
                for field in ["ref", "alt"] {
                    let base = effective_non_empty_text(event, field, path)?;
                    if base.len() != 1 || !matches!(base.as_bytes()[0], b'A' | b'C' | b'G' | b'T') {
                        return Err(effective_error(path, "alelo directo inválido"));
                    }
                }
            }
            "insertion" => {
                let _ = effective_nonnegative_u64(event, "point", path)?;
                if effective_non_empty_text(event, "bases", path)?
                    .bytes()
                    .any(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
                {
                    return Err(effective_error(path, "bases de inserción inválidas"));
                }
            }
            "deletion" => {
                let _ = effective_nonnegative_u64(event, "start", path)?;
                let _ = effective_nonnegative_u64(event, "end", path)?;
            }
            _ => unreachable!("kind was checked"),
        }
    }
    Ok(event_ids)
}

/// Comprueba la sección del tartamudeo de la configuración efectiva.
///
/// Lo que de verdad discrimina aquí es que la **procedencia** viaje con los
/// valores: unas cifras embebidas sin fuente, sin especie y sin tecnología no
/// se distinguen de unas inventadas, y un coeficiente sin la marca de no
/// calibrado pasaría por medido.
fn validate_effective_stutter(
    stutter: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    validate_effective_closed_object(
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
            "locus_rates",
        ],
        "technical_event.stutter",
        path,
    )?;
    if effective_non_empty_text(stutter, "algorithm", path)? != crate::stutter::STUTTER_ALGORITHM
        || effective_non_empty_text(stutter, "rng_catalog", path)?
            != crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
    {
        return Err(effective_error(path, "algoritmo de tartamudeo inválido"));
    }
    let calibration = effective_object(stutter, "calibration", path)?;
    validate_effective_closed_object(
        calibration,
        &[
            "source",
            "species",
            "technology",
            "column",
            "known_mismatch",
        ],
        "technical_event.stutter.calibration",
        path,
    )?;
    if effective_non_empty_text(calibration, "source", path)?
        != crate::stutter::STUTTER_CALIBRATION_SOURCE
        || effective_non_empty_text(calibration, "species", path)?
            != crate::stutter::STUTTER_CALIBRATION_SPECIES
        || effective_non_empty_text(calibration, "technology", path)?
            != crate::stutter::STUTTER_CALIBRATION_TECHNOLOGY
        || effective_non_empty_text(calibration, "column", path)?
            != crate::stutter::STUTTER_CALIBRATION_COLUMN
    {
        return Err(effective_error(
            path,
            "la procedencia de la calibración del tartamudeo no es la publicada",
        ));
    }
    validate_effective_tract_length_dependence(
        effective_object(stutter, "tract_length_dependence", path)?,
        path,
    )?;
    let coefficient = effective_object(stutter, "run_length_coefficient", path)?;
    validate_effective_closed_object(
        coefficient,
        &["value", "calibrated"],
        "technical_event.stutter.run_length_coefficient",
        path,
    )?;
    if coefficient.get("calibrated").and_then(Value::as_bool)
        != Some(crate::stutter::STUTTER_RUN_LENGTH_COEFFICIENT_CALIBRATED)
    {
        return Err(effective_error(
            path,
            "el coeficiente de la tirada tiene que publicarse como no calibrado",
        ));
    }
    let rates = stutter
        .get("motif_length_rates")
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(path, "faltan las tasas por longitud de motivo"))?;
    if rates.len() != crate::stutter::STUTTER_MOTIF_LENGTH_MAX {
        return Err(effective_error(
            path,
            "la calibración publica exactamente seis tasas por longitud de motivo",
        ));
    }
    let steps = stutter
        .get("step_weights")
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(path, "faltan los pesos de los pasos"))?;
    let declared = steps
        .iter()
        .map(|entry| entry.get("units").and_then(Value::as_i64))
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| effective_error(path, "un peso publicado no declara sus unidades"))?;
    if declared != crate::stutter::STUTTER_ADMITTED_STEPS {
        return Err(effective_error(
            path,
            "los pesos publicados no cubren exactamente la estructura del algoritmo",
        ));
    }
    Ok(())
}

/// Comprueba la dependencia con la longitud del tracto.
///
/// # Lo que de verdad discrimina aquí
///
/// Tres cosas, y ninguna es decorativa. Que la **marca de tosquedad** esté y
/// sea cierta: una curva leída a ojo de una figura, presentada sin ella, se
/// leería como una medida. Que se declare que **no se extrapola**, porque un
/// lector que no lo sepa creerá que un tracto de 30 unidades tartamudea más
/// que uno de 15, y no es así. Y que cada longitud de motivo diga **si tiene
/// curva medida o no**, con su longitud de referencia cuando la tiene, porque
/// los motivos de tres a seis pares de base no tienen dependencia ninguna y
/// callarlo invitaría a suponer que sí.
fn validate_effective_tract_length_dependence(
    dependence: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    validate_effective_closed_object(
        dependence,
        &[
            "source",
            "coarse",
            "coarse_reason",
            "extrapolated",
            "by_motif_length",
        ],
        "technical_event.stutter.tract_length_dependence",
        path,
    )?;
    if effective_non_empty_text(dependence, "source", path)?
        != crate::stutter::STUTTER_TRACT_LENGTH_SOURCE
    {
        return Err(effective_error(
            path,
            "la procedencia de la curva de longitud de tracto no es la publicada",
        ));
    }
    if dependence.get("coarse") != Some(&Value::Bool(crate::stutter::STUTTER_TRACT_LENGTH_COARSE)) {
        return Err(effective_error(
            path,
            "la dependencia con la longitud del tracto se publica marcada como tosca",
        ));
    }
    if effective_non_empty_text(dependence, "coarse_reason", path)?
        != crate::stutter::STUTTER_TRACT_LENGTH_COARSE_REASON
    {
        return Err(effective_error(
            path,
            "el motivo de la tosquedad de la curva no es el publicado",
        ));
    }
    if dependence.get("extrapolated") != Some(&Value::Bool(false)) {
        return Err(effective_error(
            path,
            "fuera del tramo medido se mantiene el extremo y no se extrapola",
        ));
    }
    let entries = dependence
        .get("by_motif_length")
        .and_then(Value::as_array)
        .ok_or_else(|| {
            effective_error(path, "falta la tabla de longitudes de tracto de referencia")
        })?;
    if entries.len() != crate::stutter::STUTTER_MOTIF_LENGTH_MAX {
        return Err(effective_error(
            path,
            "la tabla de longitudes de referencia no cubre el tramo calibrado",
        ));
    }
    for (index, entry) in entries.iter().enumerate() {
        let entry = entry.as_object().ok_or_else(|| {
            effective_error(path, "una fila de la tabla de longitudes no es objeto")
        })?;
        validate_effective_closed_object(
            entry,
            &[
                "motif_length_bases",
                "measured",
                "reference_tract_bases",
                "measured_units_min",
                "measured_units_max",
                "factors",
            ],
            "technical_event.stutter.tract_length_dependence.by_motif_length",
            path,
        )?;
        let motif_length = (index as u64) + 1;
        if effective_nonnegative_u64(entry, "motif_length_bases", path)? != motif_length {
            return Err(effective_error(
                path,
                "la tabla de longitudes no va en orden de longitud de motivo",
            ));
        }
        let measured = entry
            .get("measured")
            .and_then(Value::as_bool)
            .ok_or_else(|| effective_error(path, "una fila no declara si tiene curva medida"))?;
        let expected = crate::stutter::STUTTER_TRACT_LENGTH_REFERENCE_BASES[index];
        if measured != expected.is_some() {
            return Err(effective_error(
                path,
                "una fila dice tener curva medida y no la tiene, o al revés",
            ));
        }
        let published = entry.get("reference_tract_bases").and_then(|value| {
            if value.is_null() {
                None
            } else {
                value.as_u64()
            }
        });
        if published != expected {
            return Err(effective_error(
                path,
                "una longitud de tracto de referencia no es la del algoritmo",
            ));
        }
        let range = crate::stutter::tract_length_measured_units(motif_length);
        let published_min = entry.get("measured_units_min").and_then(Value::as_u64);
        let published_max = entry.get("measured_units_max").and_then(Value::as_u64);
        if (published_min, published_max) != (range.map(|(min, _)| min), range.map(|(_, max)| max))
        {
            return Err(effective_error(
                path,
                "el tramo medido publicado no es el del algoritmo",
            ));
        }
        // **Cada factor, uno a uno.** Publicar la tabla y no comprobarla
        // dejaría que un documento llevara una curva distinta de la que la
        // ejecución usó, y nadie podría distinguirlo de uno correcto.
        let factors = entry
            .get("factors")
            .and_then(Value::as_array)
            .ok_or_else(|| effective_error(path, "falta la tabla de factores"))?;
        let expected_len = range.map_or(0, |(min, max)| (max - min + 1) as usize);
        if factors.len() != expected_len {
            return Err(effective_error(
                path,
                "la tabla de factores no cubre el tramo medido",
            ));
        }
        for (offset, factor) in factors.iter().enumerate() {
            let factor = factor.as_object().ok_or_else(|| {
                effective_error(path, "una fila de la tabla de factores no es objeto")
            })?;
            validate_effective_closed_object(
                factor,
                &["units", "factor"],
                "technical_event.stutter.tract_length_dependence.by_motif_length.factors",
                path,
            )?;
            let units = effective_nonnegative_u64(factor, "units", path)?;
            let expected_units = range.map_or(0, |(min, _)| min) + offset as u64;
            let published = factor
                .get("factor")
                .and_then(Value::as_f64)
                .ok_or_else(|| effective_error(path, "falta un factor de la tabla"))?;
            if units != expected_units
                || published != crate::stutter::tract_length_factor(motif_length, units)
            {
                return Err(effective_error(
                    path,
                    "un factor de la curva no es el del algoritmo",
                ));
            }
        }
    }
    Ok(())
}

/// Comprueba la sección del tartamudeo clonal de la configuración efectiva.
///
/// # Lo que de verdad discrimina aquí
///
/// Que **la probabilidad publicada salga de las dos fracciones publicadas y
/// del número de ciclos publicado**. Sin esta comprobación, el documento
/// podría llevar tres números que no tienen nada que ver entre sí y nadie
/// podría distinguirlo de uno correcto. Y que la marca de tosquedad esté: una
/// tasa que arrastra efectos de individuo y de lote presentada sin ella se
/// leería como una medida.
fn validate_effective_clonal_stutter(
    stutter: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    validate_effective_closed_object(
        stutter,
        &[
            "algorithm",
            "cycles",
            "formula",
            "coarse",
            "coarse_reason",
            "calibration",
            "per_cycle_slip",
        ],
        "technical_event.artifacts.pcr.stutter",
        path,
    )?;
    if effective_non_empty_text(stutter, "algorithm", path)?
        != crate::clonal_stutter::CLONAL_STUTTER_ALGORITHM
        || effective_non_empty_text(stutter, "formula", path)?
            != crate::clonal_stutter::CLONAL_STUTTER_FORMULA
    {
        return Err(effective_error(
            path,
            "algoritmo o fórmula del tartamudeo clonal inválidos",
        ));
    }
    if stutter.get("coarse") != Some(&Value::Bool(true)) {
        return Err(effective_error(
            path,
            "las tasas por ciclo se publican marcadas como toscas",
        ));
    }
    let cycles = effective_nonnegative_u64(stutter, "cycles", path)?;
    let cycles =
        u32::try_from(cycles).map_err(|_| effective_error(path, "ciclos fuera de dominio"))?;
    let entries = stutter
        .get("per_cycle_slip")
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(path, "falta la tabla de probabilidades por ciclo"))?;
    if entries.len() != crate::stutter::STUTTER_MOTIF_LENGTH_MAX {
        return Err(effective_error(
            path,
            "la tabla de probabilidades por ciclo no cubre el tramo calibrado",
        ));
    }
    for (index, entry) in entries.iter().enumerate() {
        let entry = entry
            .as_object()
            .ok_or_else(|| effective_error(path, "una fila de la tabla por ciclo no es objeto"))?;
        validate_effective_closed_object(
            entry,
            &[
                "motif_length_bases",
                "with_pcr_fraction",
                "without_pcr_fraction",
                "probability",
                "source",
            ],
            "technical_event.artifacts.pcr.stutter.per_cycle_slip",
            path,
        )?;
        if effective_nonnegative_u64(entry, "motif_length_bases", path)? != (index as u64) + 1 {
            return Err(effective_error(
                path,
                "la tabla por ciclo no va en orden de longitud de motivo",
            ));
        }
        let source = effective_non_empty_text(entry, "source", path)?;
        let with_pcr = entry
            .get("with_pcr_fraction")
            .and_then(Value::as_f64)
            .ok_or_else(|| effective_error(path, "falta la fracción con PCR"))?;
        let without_pcr = entry
            .get("without_pcr_fraction")
            .and_then(Value::as_f64)
            .ok_or_else(|| effective_error(path, "falta la fracción sin PCR"))?;
        let probability = entry
            .get("probability")
            .and_then(Value::as_f64)
            .ok_or_else(|| effective_error(path, "falta la probabilidad por ciclo"))?;
        match source {
            "declared" => {}
            "derived" => {
                let expected = crate::clonal_stutter::per_cycle_slip_probability(
                    with_pcr,
                    without_pcr,
                    cycles,
                );
                if probability != expected {
                    return Err(effective_error(
                        path,
                        "una probabilidad por ciclo no sale de las fracciones publicadas",
                    ));
                }
            }
            _ => {
                return Err(effective_error(
                    path,
                    "la procedencia de una probabilidad por ciclo no es ni derivada ni declarada",
                ));
            }
        }
    }
    Ok(())
}

/// Comprueba la sección de amplificación de la configuración efectiva.
///
/// Lo que de verdad discrimina aquí es que la **derivación** viaje con los
/// valores: una `m` publicada sin el método ni las iteraciones que la
/// produjeron no se distingue de una inventada, y quien quiera repetir la
/// ejecución en otra máquina no tendría con qué comprobarla. Se comprueba
/// además que la `m` publicada reproduzca de verdad la fracción declarada, que
/// es lo que convierte la ida y vuelta en una afirmación verificable.
fn validate_effective_pcr(
    pcr: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    // El catálogo que se publica es el de las decisiones que esta
    // amplificación abre: el v6 en cuanto lleva tartamudeo clonal, el v5
    // mientras no lo lleve. Un v6 sin la sección, o un v5 con ella, sería un
    // documento en el que el catálogo y el contenido se contradicen. Y sin
    // loci de repetición la sección **no aparece**, de modo que una
    // biblioteca amplificada sin microsatélites publica exactamente el
    // documento de antes de F4.10.
    let opens_lineage = pcr.contains_key("stutter");
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
    ];
    if opens_lineage {
        fields.push("stutter");
    }
    validate_effective_closed_object(pcr, &fields, "technical_event.artifacts.pcr", path)?;
    let expected_catalog = if opens_lineage {
        crate::randomness::F4_TECHNICAL_RNG_CATALOG_V6
    } else {
        crate::randomness::F4_TECHNICAL_RNG_CATALOG_V5
    };
    if effective_non_empty_text(pcr, "algorithm", path)?
        != crate::pcr_amplification::PCR_AMPLIFICATION_ALGORITHM
        || effective_non_empty_text(pcr, "rng_catalog", path)? != expected_catalog
    {
        return Err(effective_error(path, "algoritmo de amplificación inválido"));
    }
    if opens_lineage {
        validate_effective_clonal_stutter(effective_object(pcr, "stutter", path)?, path)?;
    }
    let inversion = effective_object(pcr, "duplicate_fraction_inversion", path)?;
    validate_effective_closed_object(
        inversion,
        &["method", "iterations"],
        "technical_event.artifacts.pcr.duplicate_fraction_inversion",
        path,
    )?;
    if effective_non_empty_text(inversion, "method", path)?
        != crate::pcr_amplification::PCR_DUPLICATE_FRACTION_INVERSION_METHOD
    {
        return Err(effective_error(
            path,
            "método de inversión de la fracción de duplicados inválido",
        ));
    }
    let declared = pcr
        .get("declared_duplicate_fraction")
        .and_then(Value::as_f64)
        .ok_or_else(|| effective_error(path, "falta la fracción de duplicados declarada"))?;
    let mean_pairs = pcr
        .get("mean_pairs_per_molecule")
        .and_then(Value::as_f64)
        .ok_or_else(|| effective_error(path, "falta la m derivada"))?;
    let iterations = effective_nonnegative_u64(inversion, "iterations", path)?;
    // Con fracción cero la inversión **no itera**: devuelve cero exacto por
    // retorno anticipado, y un residuo minúsculo cambiaría todos los FASTQ.
    let expected_iterations = if declared == 0.0 {
        0
    } else {
        u64::from(crate::pcr_amplification::PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS)
    };
    if iterations != expected_iterations || (declared == 0.0 && mean_pairs != 0.0) {
        return Err(effective_error(
            path,
            "la inversión publicada no corresponde a la fracción declarada",
        ));
    }
    let expected = pcr
        .get("expected_duplicate_fraction")
        .and_then(Value::as_f64)
        .ok_or_else(|| effective_error(path, "falta la fracción esperada"))?;
    if expected != crate::pcr_amplification::expected_duplicate_fraction(mean_pairs) {
        return Err(effective_error(
            path,
            "la fracción esperada publicada no sale de la m publicada",
        ));
    }
    let classes = pcr
        .get("efficiency_classes")
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(path, "faltan las clases de eficiencia"))?;
    if classes.is_empty()
        || classes.len() > crate::technical_profile::PCR_MAXIMUM_EFFICIENCY_CLASSES as usize
    {
        return Err(effective_error(
            path,
            "las clases de eficiencia publicadas están fuera del contrato",
        ));
    }
    let knots = pcr
        .get("gc_relative_efficiency")
        .and_then(Value::as_array)
        .ok_or_else(|| effective_error(path, "falta la curva de composición"))?;
    if knots.len() < 2 {
        return Err(effective_error(
            path,
            "la curva de composición publicada necesita al menos dos nudos",
        ));
    }
    Ok(())
}

fn validate_effective_technical_event(
    effective: &serde_json::Map<String, Value>,
    initial_individual_ids: &BTreeSet<IndividualId>,
    direct_event_ids: &BTreeSet<EventId>,
    may_declare_repeats: bool,
    path: &Path,
) -> Result<(), Diagnostic> {
    let event = effective_object(effective, "technical_event", path)?;
    const FIELDS: [&str; 12] = [
        "event_id",
        "individual_ids",
        "profile",
        "read_length",
        "quality",
        "sequencing_errors",
        "artifacts",
        "biological_contamination",
        "libraries",
        "rng_catalog",
        "id_catalog",
        "target",
    ];
    const FIELDS_WITH_STUTTER: [&str; 13] = [
        "event_id",
        "individual_ids",
        "profile",
        "read_length",
        "quality",
        "sequencing_errors",
        "artifacts",
        "biological_contamination",
        "libraries",
        "rng_catalog",
        "id_catalog",
        "target",
        "stutter",
    ];
    // La sección del tartamudeo y el catálogo v4 **van juntos**, igual que en
    // el manifiesto del evento: uno sin el otro sería un documento que dice
    // haber tartamudeado bajo un catálogo que no tiene esa etapa, o al revés.
    let declares_stutter = event.contains_key("stutter");
    // La sección del tartamudeo **solo cabe en la versión que admite
    // repeticiones**: sin loci declarados no hay nada que tartamudear, y sin
    // repeticiones no hay loci. Admitirla en una versión anterior dejaría al
    // lector de Rust diciendo algo distinto del esquema publicado, que la
    // rechaza como campo desconocido.
    if declares_stutter && !may_declare_repeats {
        return Err(effective_error(
            path,
            "la sección de tartamudeo no pertenece a esta versión de configuración efectiva",
        ));
    }
    validate_effective_closed_object(
        event,
        if declares_stutter {
            &FIELDS_WITH_STUTTER
        } else {
            &FIELDS
        },
        "technical_event",
        path,
    )?;
    let event_id = EventId::parse(
        effective_non_empty_text(event, "event_id", path)?,
        "effective_config",
    )
    .map_err(|_| effective_error(path, "event_id técnico inválido"))?;
    if direct_event_ids.contains(&event_id) {
        return Err(effective_error(
            path,
            "event_id técnico colisiona con un evento directo",
        ));
    }
    let read_length = effective_nonnegative_u64(event, "read_length", path)?;
    let profile = effective_object(event, "profile", path)?;
    let expected_profile = profile
        .get("id")
        .and_then(Value::as_str)
        .and_then(f4_technical_profile_object);
    if expected_profile.as_ref() != Some(&Value::Object(profile.clone()))
        || !(1..=150).contains(&read_length)
    {
        return Err(effective_error(path, "perfil técnico mínimo inválido"));
    }
    // La amplificación estrena el v5, y **manda sobre el v4**: un evento que
    // amplifica y además tartamudea abre las dos etapas y publica el más alto.
    let declares_pcr = effective_object(event, "artifacts", path)
        .ok()
        .and_then(|artifacts| artifacts.get("pcr"))
        .and_then(Value::as_object)
        .and_then(|pcr| pcr.get("applicable"))
        == Some(&Value::Bool(true));
    let expected_catalog = if declares_pcr {
        crate::randomness::F4_TECHNICAL_RNG_CATALOG_V5
    } else if declares_stutter {
        crate::randomness::F4_TECHNICAL_RNG_CATALOG_V4
    } else {
        crate::randomness::F4_TECHNICAL_RNG_CATALOG_V3
    };
    if effective_non_empty_text(event, "rng_catalog", path)? != expected_catalog
        || event.get("id_catalog") != Some(&f4_id_catalog_value())
    {
        return Err(effective_error(path, "catálogo técnico F4 inválido"));
    }
    if declares_stutter {
        validate_effective_stutter(effective_object(event, "stutter", path)?, path)?;
    }
    let contamination = effective_object(event, "biological_contamination", path)?;
    validate_effective_closed_object(
        contamination,
        &["percent_range", "recipient_ids"],
        "technical_event.biological_contamination",
        path,
    )?;
    let percent_range = effective_object(contamination, "percent_range", path)?;
    validate_effective_closed_object(
        percent_range,
        &["minimum", "maximum"],
        "technical_event.biological_contamination.percent_range",
        path,
    )?;
    let minimum_percent = percent_range.get("minimum").and_then(Value::as_f64);
    let maximum_percent = percent_range.get("maximum").and_then(Value::as_f64);
    if !matches!((minimum_percent, maximum_percent), (Some(minimum), Some(maximum)) if minimum.is_finite() && maximum.is_finite() && (0.0..=100.0).contains(&minimum) && (0.0..=100.0).contains(&maximum) && minimum <= maximum)
    {
        return Err(effective_error(
            path,
            "rango de contaminación biológica inválido",
        ));
    }
    let quality = effective_object(event, "quality", path)?;
    validate_effective_closed_object(
        quality,
        &["acgt_phred", "n_phred"],
        "technical_event.quality",
        path,
    )?;
    if effective_nonnegative_u64(quality, "acgt_phred", path)? != 40
        || effective_nonnegative_u64(quality, "n_phred", path)? != 0
    {
        return Err(effective_error(path, "calidad técnica mínima inválida"));
    }
    let sequencing_errors = effective_object(event, "sequencing_errors", path)?;
    validate_effective_closed_object(
        sequencing_errors,
        &[
            "requested_mode",
            "effective_mode",
            "qualities_represent_realized_error_probability",
        ],
        "technical_event.sequencing_errors",
        path,
    )?;
    let requested_mode = effective_non_empty_text(sequencing_errors, "requested_mode", path)?;
    let effective_mode = effective_non_empty_text(sequencing_errors, "effective_mode", path)?;
    let quality_mark = sequencing_errors
        .get("qualities_represent_realized_error_probability")
        .and_then(Value::as_bool);
    if requested_mode != effective_mode
        || !matches!(requested_mode, "enabled" | "disabled")
        || quality_mark != Some(requested_mode == "enabled")
    {
        return Err(effective_error(
            path,
            "modo de errores de secuenciación efectivo inválido",
        ));
    }
    let artifacts = effective_object(event, "artifacts", path)?;
    validate_effective_closed_object(
        artifacts,
        &[
            "chimera",
            "index_hopping",
            "pcr",
            "optical_duplicates",
            "gc",
            "adapters",
        ],
        "technical_event.artifacts",
        path,
    )?;
    let chimera = effective_object(artifacts, "chimera", path)?;
    validate_effective_closed_object(
        chimera,
        &["rate"],
        "technical_event.artifacts.chimera",
        path,
    )?;
    if !matches!(chimera.get("rate").and_then(Value::as_f64), Some(value) if value.is_finite() && (0.0..=1.0).contains(&value))
    {
        return Err(effective_error(
            path,
            "artefacto de quimera efectivo inválido",
        ));
    }
    let index_hopping = effective_object(artifacts, "index_hopping", path)?;
    validate_effective_closed_object(
        index_hopping,
        &["rate"],
        "technical_event.artifacts.index_hopping",
        path,
    )?;
    if !matches!(index_hopping.get("rate").and_then(Value::as_f64), Some(value) if value.is_finite() && (0.0..=1.0).contains(&value))
    {
        return Err(effective_error(
            path,
            "artefacto de index hopping efectivo inválido",
        ));
    }
    let pcr = effective_object(artifacts, "pcr", path)?;
    if declares_pcr {
        validate_effective_pcr(pcr, path)?;
    } else {
        validate_effective_closed_object(
            pcr,
            &["applicable"],
            "technical_event.artifacts.pcr",
            path,
        )?;
        if pcr.get("applicable") != Some(&Value::Bool(false)) {
            return Err(effective_error(
                path,
                "artefactos de duplicación del perfil PCR-free inválidos",
            ));
        }
    }
    let optical_duplicates = effective_object(artifacts, "optical_duplicates", path)?;
    validate_effective_closed_object(
        optical_duplicates,
        &["realized_pairs"],
        "technical_event.artifacts.optical_duplicates",
        path,
    )?;
    if effective_nonnegative_u64(optical_duplicates, "realized_pairs", path)? != 0 {
        return Err(effective_error(
            path,
            "artefactos de duplicación del perfil PCR-free inválidos",
        ));
    }
    let gc = effective_object(artifacts, "gc", path)?;
    validate_effective_closed_object(gc, &["mode"], "technical_event.artifacts.gc", path)?;
    let adapters = effective_object(artifacts, "adapters", path)?;
    validate_effective_closed_object(
        adapters,
        &["enabled"],
        "technical_event.artifacts.adapters",
        path,
    )?;
    let gc_mode = gc.get("mode").and_then(Value::as_str);
    // El modo `curve` solo existe con amplificación: sin ella nadie declara
    // una curva y publicarlo sería afirmar un sesgo que no se aplicó.
    let gc_mode_valid = match gc_mode {
        Some("neutral") => true,
        Some("curve") => declares_pcr,
        _ => false,
    };
    if !gc_mode_valid || adapters.get("enabled") != Some(&Value::Bool(false)) {
        return Err(effective_error(
            path,
            "artefactos GC o adaptadores del perfil inicial inválidos",
        ));
    }
    let mut individuals = BTreeSet::new();
    let individual_ids = effective_array(event, "individual_ids", path)?;
    if individual_ids.is_empty() {
        return Err(effective_error(path, "technical_event sin individuos"));
    }
    for individual_id in individual_ids {
        let individual_id = individual_id
            .as_str()
            .ok_or_else(|| effective_error(path, "individual_id técnico inválido"))?;
        let individual_id = IndividualId::parse(individual_id, "effective_config")
            .map_err(|_| effective_error(path, "individual_id técnico inválido"))?;
        if !individuals.insert(individual_id.clone()) {
            return Err(effective_error(path, "individual_id técnico repetido"));
        }
        if !initial_individual_ids.contains(&individual_id) {
            return Err(effective_error(
                path,
                "individual_id técnico fuera de las poblaciones iniciales",
            ));
        }
    }
    let recipient_ids = effective_array(contamination, "recipient_ids", path)?;
    if recipient_ids.is_empty() {
        return Err(effective_error(
            path,
            "contaminación biológica sin receptoras",
        ));
    }
    let mut recipients = BTreeSet::new();
    for recipient_id in recipient_ids {
        let recipient_id = recipient_id.as_str().ok_or_else(|| {
            effective_error(path, "receptora de contaminación biológica inválida")
        })?;
        let recipient_id = IndividualId::parse(recipient_id, "effective_config")
            .map_err(|_| effective_error(path, "receptora de contaminación biológica inválida"))?;
        if !individuals.contains(&recipient_id) || !recipients.insert(recipient_id) {
            return Err(effective_error(
                path,
                "receptoras de contaminación inválidas o repetidas",
            ));
        }
    }
    let libraries = effective_array(event, "libraries", path)?;
    if libraries.len() != individuals.len() {
        return Err(effective_error(
            path,
            "libraries no coincide con individual_ids",
        ));
    }
    let mut library_individuals = BTreeSet::new();
    let mut i7 = BTreeSet::new();
    let mut i5 = BTreeSet::new();
    for library in libraries {
        let library = library
            .as_object()
            .ok_or_else(|| effective_error(path, "library inválida"))?;
        validate_effective_closed_object(
            library,
            &[
                "individual_id",
                "library_id",
                "i7",
                "i5",
                "is_contamination_recipient",
            ],
            "technical_event.libraries",
            path,
        )?;
        let individual_id = IndividualId::parse(
            effective_non_empty_text(library, "individual_id", path)?,
            "effective_config",
        )
        .map_err(|_| effective_error(path, "library individual_id inválido"))?;
        let expected_library = derive_technical_identity_plan(
            RootSeed::parse(
                "root",
                effective_non_empty_text(effective_object(effective, "seed", path)?, "root", path)?,
            )?,
            &ReferenceId::parse(
                effective_non_empty_text(
                    effective_object(effective, "reference", path)?,
                    "id",
                    path,
                )?,
                "effective_config",
            )?,
            event_id.clone(),
            &individuals.iter().cloned().collect::<Vec<_>>(),
        )
        .map_err(|_| effective_error(path, "biblioteca técnica no derivable"))?
        .library(&individual_id)
        .ok_or_else(|| effective_error(path, "library técnica no derivable"))?
        .library_id()
        .as_str()
        .to_owned();
        let valid_udi = |field: &str| {
            effective_non_empty_text(library, field, path).map(|value| {
                value.len() == 8
                    && value
                        .bytes()
                        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            })
        };
        if !individuals.contains(&individual_id)
            || !library_individuals.insert(individual_id.clone())
            || effective_non_empty_text(library, "library_id", path)? != expected_library
            || library.get("is_contamination_recipient")
                != Some(&Value::Bool(recipients.contains(&individual_id)))
            || !valid_udi("i7")?
            || !valid_udi("i5")?
            || !i7.insert(effective_non_empty_text(library, "i7", path)?.to_owned())
            || !i5.insert(effective_non_empty_text(library, "i5", path)?.to_owned())
        {
            return Err(effective_error(path, "biblioteca UDI inválida o repetida"));
        }
    }
    let target = effective_object(event, "target", path)?;
    validate_effective_closed_object(target, &["coverage_wgs"], "technical_event.target", path)?;
    if !matches!(target.get("coverage_wgs").and_then(Value::as_f64), Some(value) if value.is_finite() && value >= 0.0)
    {
        return Err(effective_error(path, "cobertura WGS técnica inválida"));
    }
    Ok(())
}

fn validate_manifest_authorizations(value: Option<&Value>, path: &Path) -> Result<(), Diagnostic> {
    let value = value.ok_or_else(|| manifest_error(path, "falta authorizations"))?;
    validate_authorizations_value(value).map_err(|detail| manifest_error_from_reason(path, detail))
}

fn validate_authorizations_value(value: &Value) -> Result<(), &'static str> {
    let object = value
        .as_object()
        .ok_or_else(|| reason("authorizations inválido"))?;
    if object.len() != 2
        || !matches!(object.get("allow_empty_snv_matrix"), Some(Value::Bool(_)))
        || !matches!(
            object.get("allow_incompatible_genotypes"),
            Some(Value::Bool(_))
        )
    {
        return Err(reason("authorizations inválido"));
    }
    Ok(())
}

fn validate_configuration_input(
    object: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<String, Diagnostic> {
    const FIELDS: [&str; 3] = ["kind", "original_path", "artifact_id"];
    validate_closed_object(object, &FIELDS, "configuration_input", path)?;
    if object.get("kind").and_then(Value::as_str) != Some("general") {
        return Err(manifest_error(path, "configuration_input.kind inválido"));
    }
    required_non_empty_text(object, "original_path", path)?;
    let id = required_non_empty_text(object, "artifact_id", path)?;
    require_artifact_kind(catalog, id, "config_general_toml", path)?;
    Ok(id.to_owned())
}

fn validate_references(
    values: &[Value],
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<BTreeMap<String, String>, Diagnostic> {
    const FIELDS: [&str; 2] = ["id", "fasta_artifact_id"];
    let mut references = BTreeMap::new();
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "una referencia no es objeto"))?;
        validate_closed_object(object, &FIELDS, "references[]", path)?;
        let id = required_non_empty_text(object, "id", path)?;
        let fasta_id = required_non_empty_text(object, "fasta_artifact_id", path)?;
        if references
            .insert(id.to_owned(), fasta_id.to_owned())
            .is_some()
        {
            return Err(manifest_error(path, "ID de referencia repetido"));
        }
        require_artifact_kind(catalog, fasta_id, "reference_fasta", path)?;
    }
    Ok(references)
}

/// Lo que el manifiesto publica de la referencia final (DTC-098), ya validado.
struct ValidatedManifestInsertion {
    state: UnitState,
    map_artifact_id: String,
}

/// **La sección `reference_insertion` del manifiesto general** (DTC-098,
/// REQ-MAN-07).
///
/// Comprueba la forma cerrada, que la referencia final y la de origen están
/// entre las referencias del manifiesto con los FASTA que dicen, que el mapa
/// es un artefacto de su clase y versión, y que la lista de cambios no está
/// vacía y cada cambio tiene sus cuatro campos. Que los cambios sean **los
/// declarados** y que el FASTA final sea el suministrado con esos cambios lo
/// juzga la validación física, que tiene la configuración y las bases.
fn validate_manifest_reference_insertion(
    value: &Value,
    catalog: &BTreeMap<&str, &str>,
    references: &BTreeMap<String, String>,
    path: &Path,
) -> Result<ValidatedManifestInsertion, Diagnostic> {
    const FIELDS: [&str; 9] = [
        "id",
        "source_reference_id",
        "state",
        "fasta_artifact_id",
        "map_artifact_id",
        "source_fasta_artifact_id",
        "map_document_kind",
        "map_schema_version",
        "insertions",
    ];
    const CHANGE_FIELDS: [&str; 4] = ["locus_id", "sequence_id", "insert_after", "inserted_bases"];
    let object = value
        .as_object()
        .ok_or_else(|| manifest_error(path, "reference_insertion no es objeto"))?;
    validate_closed_object(object, &FIELDS, "reference_insertion", path)?;
    let id = required_non_empty_text(object, "id", path)?;
    let source = required_non_empty_text(object, "source_reference_id", path)?;
    let fasta = required_non_empty_text(object, "fasta_artifact_id", path)?;
    let source_fasta = required_non_empty_text(object, "source_fasta_artifact_id", path)?;
    if id == source
        || references.get(id).map(String::as_str) != Some(fasta)
        || references.get(source).map(String::as_str) != Some(source_fasta)
    {
        return Err(manifest_error(
            path,
            "reference_insertion no enlaza una referencia final distinta de su origen, con los FASTA que el manifiesto registra",
        ));
    }
    let map = required_non_empty_text(object, "map_artifact_id", path)?;
    require_artifact_kind(
        catalog,
        map,
        crate::reference_tract_insertion::COORDINATE_MAP_DOCUMENT_KIND,
        path,
    )?;
    if object.get("map_document_kind").and_then(Value::as_str)
        != Some(crate::reference_tract_insertion::COORDINATE_MAP_DOCUMENT_KIND)
        || object.get("map_schema_version").and_then(Value::as_str)
            != Some(crate::reference_tract_insertion::COORDINATE_MAP_SCHEMA_VERSION)
    {
        return Err(manifest_error(
            path,
            "reference_insertion no declara la clase y la versión del mapa de coordenadas",
        ));
    }
    let state = match object.get("state").and_then(Value::as_str) {
        Some("complete") => UnitState::Complete,
        Some("failed") => UnitState::Failed,
        _ => {
            return Err(manifest_error(
                path,
                "reference_insertion no tiene un estado válido",
            ));
        }
    };
    let changes = object
        .get("insertions")
        .and_then(Value::as_array)
        .filter(|changes| !changes.is_empty())
        .ok_or_else(|| manifest_error(path, "reference_insertion no publica ningún cambio"))?;
    for change in changes {
        let well_formed = change.as_object().is_some_and(|change| {
            change.len() == CHANGE_FIELDS.len()
                && CHANGE_FIELDS
                    .iter()
                    .all(|field| change.contains_key(*field))
                && change
                    .get("locus_id")
                    .and_then(Value::as_str)
                    .is_some_and(|text| !text.is_empty())
                && change
                    .get("sequence_id")
                    .and_then(Value::as_str)
                    .is_some_and(|text| !text.is_empty())
                && change
                    .get("insert_after")
                    .and_then(Value::as_u64)
                    .is_some_and(|value| value >= 1)
                && change
                    .get("inserted_bases")
                    .and_then(Value::as_u64)
                    .is_some_and(|value| value >= 1)
        });
        if !well_formed {
            return Err(manifest_error(
                path,
                "un cambio de reference_insertion no tiene exactamente locus, secuencia, punto y bases",
            ));
        }
    }
    Ok(ValidatedManifestInsertion {
        state,
        map_artifact_id: map.to_owned(),
    })
}

struct ValidatedManifestReductions {
    states: Vec<UnitState>,
    derived_reference_ids: BTreeSet<String>,
    requested_artifact_ids: BTreeSet<String>,
}

fn validate_reference_reductions(
    values: &[Value],
    catalog: &BTreeMap<&str, &str>,
    references: &BTreeMap<String, String>,
    path: &Path,
) -> Result<ValidatedManifestReductions, Diagnostic> {
    const FIELDS: [&str; 14] = [
        "id",
        "source_reference_id",
        "eligible_bases",
        "region_size_bases",
        "target",
        "root_seed",
        "selection_subseed",
        "selection_order",
        "state",
        "fasta_artifact_id",
        "map_artifact_id",
        "source_fasta_artifact_id",
        "map_document_kind",
        "map_schema_version",
    ];
    const TARGET_FIELDS: [&str; 3] = ["mode", "declared", "bases"];
    let mut reduction_ids = BTreeSet::new();
    let mut states = Vec::with_capacity(values.len());
    let mut derived_reference_ids = BTreeSet::new();
    let mut requested_artifact_ids = BTreeSet::new();
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "una reducción no es objeto"))?;
        validate_closed_object(object, &FIELDS, "reference_reductions[]", path)?;
        let id = required_non_empty_text(object, "id", path)?;
        let source_id = required_non_empty_text(object, "source_reference_id", path)?;
        if !reduction_ids.insert(id) {
            return Err(manifest_error(path, "ID de reducción repetido"));
        }
        derived_reference_ids.insert(id.to_owned());
        required_positive_u64(object, "eligible_bases", path)?;
        required_positive_u64(object, "region_size_bases", path)?;
        let target = required_object(object, "target", path)?;
        validate_closed_object(
            target,
            &TARGET_FIELDS,
            "reference_reductions[].target",
            path,
        )?;
        match target.get("mode").and_then(Value::as_str) {
            Some("target_bases" | "reduction_percent") => {}
            _ => return Err(manifest_error(path, "modo de reducción inválido")),
        }
        required_non_empty_text(target, "declared", path)?;
        required_positive_u64(target, "bases", path)?;
        validate_lower_hex_64(
            required_non_empty_text(object, "root_seed", path)?,
            "root_seed",
            path,
        )?;
        if required_non_empty_text(object, "map_document_kind", path)?
            != REDUCTION_MAP_DOCUMENT_KIND
            || required_non_empty_text(object, "map_schema_version", path)?
                != REDUCTION_MAP_SCHEMA_VERSION
        {
            return Err(manifest_error(
                path,
                "algoritmo o contrato de mapa de reducción incompatible",
            ));
        }
        validate_lower_hex_64(
            required_non_empty_text(object, "selection_subseed", path)?,
            "selection_subseed",
            path,
        )?;
        let selection_order = required_array(object, "selection_order", path)?;
        if selection_order.is_empty() {
            return Err(manifest_error(path, "selection_order está vacío"));
        }
        for order in selection_order {
            if order.as_u64().is_none_or(|value| value == 0) {
                return Err(manifest_error(path, "selection_order inválido"));
            }
        }
        let state = match object.get("state").and_then(Value::as_str) {
            Some("incomplete") => UnitState::Incomplete,
            Some("complete") => UnitState::Complete,
            Some("failed") => UnitState::Failed,
            _ => return Err(manifest_error(path, "estado de reducción inválido")),
        };
        states.push(state);
        let mut linked = BTreeMap::new();
        for (field, kind) in [
            ("fasta_artifact_id", "reference_fasta"),
            ("map_artifact_id", "reference_reduction_map"),
            ("source_fasta_artifact_id", "reference_fasta"),
        ] {
            let id = optional_text(object, field, path)?;
            if let Some(id) = id {
                require_artifact_kind(catalog, id, kind, path)?;
                if field != "source_fasta_artifact_id" {
                    requested_artifact_ids.insert(id.to_owned());
                }
            }
            linked.insert(field, id);
        }
        if state == UnitState::Complete && linked.values().any(Option::is_none) {
            return Err(manifest_error(
                path,
                "una reducción complete exige sus tres enlaces de artefacto",
            ));
        }
        if linked["fasta_artifact_id"].is_some()
            && linked["fasta_artifact_id"] == linked["source_fasta_artifact_id"]
        {
            return Err(manifest_error(
                path,
                "una reducción no puede reutilizar el mismo FASTA como fuente y derivado",
            ));
        }
        if let Some(fasta_id) = linked["fasta_artifact_id"]
            && references.get(id).map(String::as_str) != Some(fasta_id)
        {
            return Err(manifest_error(
                path,
                "la reducción no coincide con el enlace de su referencia derivada",
            ));
        }
        if let Some(source_fasta_id) = linked["source_fasta_artifact_id"]
            && references.get(source_id).map(String::as_str) != Some(source_fasta_id)
        {
            return Err(manifest_error(
                path,
                "la reducción no coincide con el enlace de su referencia fuente",
            ));
        }
    }
    Ok(ValidatedManifestReductions {
        states,
        derived_reference_ids,
        requested_artifact_ids,
    })
}

fn validate_vcf_outputs(
    values: &[Value],
    catalog: &BTreeMap<&str, &str>,
    references: &BTreeMap<String, String>,
    with_normalization_version: bool,
    path: &Path,
) -> Result<BTreeSet<String>, Diagnostic> {
    if values.len() > 1 {
        return Err(manifest_error(
            path,
            "solo se admite un VCF final por ejecución",
        ));
    }
    let mut requested = BTreeSet::new();
    let mut reference_ids = BTreeSet::new();
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "vcf_outputs contiene un valor no objeto"))?;
        const FIELDS: [&str; 7] = [
            "state",
            "reference_id",
            "samples",
            "mode",
            "vcf_artifact_id",
            "index_kind",
            "index_artifact_id",
        ];
        if with_normalization_version {
            // La versión 6 registra la del algoritmo de normalización
            // (DTC-064), obligatoria y solo la que este programa aplica.
            let mut fields = FIELDS.to_vec();
            fields.push("normalization_version");
            validate_closed_object(object, &fields, "vcf_outputs[]", path)?;
            if object.get("normalization_version").and_then(Value::as_str)
                != Some(VCF_NORMALIZATION_VERSION)
            {
                return Err(manifest_error_for_schema(
                    path,
                    "manifest-general-v6.schema.json",
                    "vcf_outputs[].normalization_version no es la versión del algoritmo de normalización VCF",
                ));
            }
        } else {
            validate_closed_object(object, &FIELDS, "vcf_outputs[]", path)?;
        }
        if object.get("state").and_then(Value::as_str) != Some("complete") {
            return Err(manifest_error(path, "un VCF publicado debe estar complete"));
        }
        let reference_id = required_non_empty_text(object, "reference_id", path)?;
        if !references.contains_key(reference_id) || !reference_ids.insert(reference_id.to_owned())
        {
            return Err(manifest_error(
                path,
                "VCF con referencia ausente o duplicada",
            ));
        }
        let samples = required_array(object, "samples", path)?;
        let mut sample_ids = BTreeSet::new();
        for sample in samples {
            let sample = sample
                .as_str()
                .filter(|sample| !sample.is_empty())
                .ok_or_else(|| manifest_error(path, "muestra VCF inválida"))?;
            if !sample_ids.insert(sample) {
                return Err(manifest_error(path, "muestra VCF duplicada"));
            }
        }
        let vcf_id = required_non_empty_text(object, "vcf_artifact_id", path)?;
        require_artifact_kind(catalog, vcf_id, "truth_vcf", path)?;
        match object.get("mode").and_then(Value::as_str) {
            Some("text") => {
                if !object.get("index_kind").is_some_and(Value::is_null)
                    || !object.get("index_artifact_id").is_some_and(Value::is_null)
                {
                    return Err(manifest_error(path, "VCF de texto no admite índice"));
                }
            }
            Some("bgzf") => {
                let expected_kind = match object.get("index_kind").and_then(Value::as_str) {
                    Some("tbi") => "vcf_tbi_index",
                    Some("csi") => "vcf_csi_index",
                    _ => return Err(manifest_error(path, "VCF BGZF exige índice TBI o CSI")),
                };
                let index_id = optional_text(object, "index_artifact_id", path)?
                    .ok_or_else(|| manifest_error(path, "VCF BGZF exige artefacto de índice"))?;
                require_artifact_kind(catalog, index_id, expected_kind, path)?;
                requested.insert(index_id.to_owned());
            }
            _ => return Err(manifest_error(path, "modo VCF inválido")),
        }
        requested.insert(vcf_id.to_owned());
    }
    Ok(requested)
}

fn validate_artifact(value: &Value, path: &Path) -> Result<ArtifactRecord, Diagnostic> {
    const FIELDS: [&str; 5] = ["id", "kind", "path", "size_bytes", "checksum"];
    const CHECKSUM_FIELDS: [&str; 2] = ["algorithm", "digest"];
    let artifact = value
        .as_object()
        .ok_or_else(|| manifest_error(path, "un artefacto no es objeto"))?;
    validate_closed_object(artifact, &FIELDS, "artifact", path)?;
    let id = required_non_empty_text(artifact, "id", path)?;
    let kind = required_non_empty_text(artifact, "kind", path)?;
    let relative = required_non_empty_text(artifact, "path", path)?;
    validate_relative_link(relative)?;
    let checksum = required_object(artifact, "checksum", path)?;
    validate_closed_object(checksum, &CHECKSUM_FIELDS, "artifact.checksum", path)?;
    if checksum.get("algorithm").and_then(Value::as_str) != Some("sha256") {
        return Err(manifest_error(path, "checksum no usa sha256"));
    }
    let digest = required_non_empty_text(checksum, "digest", path)?;
    validate_lower_hex_64(digest, "digest", path)?;
    Ok(ArtifactRecord {
        id: id.to_owned(),
        kind: kind.to_owned(),
        path: relative.to_owned(),
        size_bytes: artifact
            .get("size_bytes")
            .and_then(Value::as_u64)
            .ok_or_else(|| manifest_error(path, "size_bytes inválido"))?,
        digest: digest.to_owned(),
    })
}

fn validate_diagnostics(
    values: &[Value],
    catalog_version: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 2] = ["code", "message"];
    for value in values {
        let object = value
            .as_object()
            .ok_or_else(|| manifest_error(path, "un diagnóstico no es objeto"))?;
        validate_closed_object(object, &FIELDS, "diagnostics[]", path)?;
        let code = required_non_empty_text(object, "code", path)?;
        validate_diagnostic_code(code, path)?;
        let message = required_non_empty_text(object, "message", path)?;
        if !diagnostic_matches_catalog_version(catalog_version, code, message) {
            return Err(manifest_error(
                path,
                format!("el diagnóstico {code} no coincide con diagnostics-v{catalog_version}.tsv"),
            ));
        }
    }
    Ok(())
}

fn require_artifact_kind(
    catalog: &BTreeMap<&str, &str>,
    id: &str,
    kind: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    match catalog.get(id) {
        Some(observed) if *observed == kind => Ok(()),
        Some(observed) => Err(manifest_error(
            path,
            format!("el enlace {id} exige {kind}, pero el catálogo declara {observed}"),
        )),
        None => Err(manifest_error(path, format!("el enlace {id} no existe"))),
    }
}

fn require_effective_artifact_kind(
    catalog: &BTreeMap<&str, &str>,
    id: &str,
    kind: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    match catalog.get(id) {
        Some(observed) if *observed == kind => Ok(()),
        Some(observed) => Err(effective_error(
            path,
            format!("el enlace {id} exige {kind}, pero el catálogo declara {observed}"),
        )),
        None => Err(effective_error(path, format!("el enlace {id} no existe"))),
    }
}

fn validate_closed_object(
    object: &serde_json::Map<String, Value>,
    fields: &[&str],
    context: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if object.keys().any(|key| !fields.contains(&key.as_str())) {
        return Err(manifest_error(
            path,
            format!("campo desconocido en {context}"),
        ));
    }
    Ok(())
}

fn required_object<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a serde_json::Map<String, Value>, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_object)
        .ok_or_else(|| manifest_error(path, format!("{field} no es objeto")))
}

fn required_array<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a [Value], Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .ok_or_else(|| manifest_error(path, format!("{field} no es lista")))
}

fn required_non_empty_text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    let value = object
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| manifest_error(path, format!("falta {field}")))?;
    if value.is_empty() {
        return Err(manifest_error(path, format!("{field} está vacío")));
    }
    Ok(value)
}

fn optional_text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<Option<&'a str>, Diagnostic> {
    match object.get(field) {
        Some(Value::Null) => Ok(None),
        Some(Value::String(value)) if !value.is_empty() => Ok(Some(value)),
        _ => Err(manifest_error(path, format!("{field} inválido"))),
    }
}

fn required_positive_u64(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<u64, Diagnostic> {
    match object.get(field).and_then(Value::as_u64) {
        Some(value) if value > 0 => Ok(value),
        _ => Err(manifest_error(path, format!("{field} inválido"))),
    }
}

fn validate_lower_hex_64(value: &str, field: &str, path: &Path) -> Result<(), Diagnostic> {
    if value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Ok(());
    }
    Err(manifest_error(path, format!("{field} SHA-256 inválido")))
}

fn validate_diagnostic_code(code: &str, path: &Path) -> Result<(), Diagnostic> {
    let mut parts = code.split('-');
    let valid_category = matches!(parts.next(), Some("GPS"))
        && matches!(
            parts.next(),
            Some("CONF" | "INPUT" | "FMT" | "DOMAIN" | "RES" | "TOOL" | "IO" | "INTERNAL")
        );
    let valid_number = parts.next().is_some_and(|number| {
        number.len() == 3 && number.bytes().all(|byte| byte.is_ascii_digit())
    });
    if valid_category && valid_number && parts.next().is_none() {
        return Ok(());
    }
    Err(manifest_error(path, "código de diagnóstico inválido"))
}

fn validate_effective_seed(
    object: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 2] = ["root", "origin"];
    validate_effective_closed_object(object, &FIELDS, "seed", path)?;
    validate_lower_hex_64(
        effective_non_empty_text(object, "root", path)?,
        "seed.root",
        path,
    )?;
    match object.get("origin").and_then(Value::as_str) {
        Some("declared" | "system_csprng") => Ok(()),
        _ => Err(effective_error(path, "seed.origin inválido")),
    }
}

fn validate_effective_normalization(
    object: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 3] = ["tolerance", "canonical_order", "sum_algorithm"];
    validate_effective_closed_object(object, &FIELDS, "normalization", path)?;
    let tolerance = object
        .get("tolerance")
        .and_then(Value::as_f64)
        .filter(|value| value.is_finite() && (0.0..=1e-6).contains(value))
        .ok_or_else(|| effective_error(path, "normalization.tolerance inválida"))?;
    let _ = tolerance;
    if object.get("canonical_order").and_then(Value::as_str)
        != Some("population_id,contig,position,ref,alt,variant_id,source_id")
        || object.get("sum_algorithm").and_then(Value::as_str) != Some("neumaier_sum_v1")
    {
        return Err(effective_error(
            path,
            "metadatos de normalización incompatibles",
        ));
    }
    Ok(())
}

fn validate_effective_reference(
    object: &serde_json::Map<String, Value>,
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<(String, String), Diagnostic> {
    const FIELDS: [&str; 8] = [
        "id",
        "source",
        "original_fasta_path",
        "canonical_format",
        "total_genome_size",
        "composition",
        "sequences",
        "fasta_artifact_id",
    ];
    validate_effective_closed_object(object, &FIELDS, "reference", path)?;
    let id = effective_non_empty_text(object, "id", path)?.to_owned();
    let source = match object.get("source").and_then(Value::as_str) {
        Some(source @ ("simulated" | "supplied")) => source,
        _ => return Err(effective_error(path, "reference.source inválido")),
    };
    let declared_total = effective_positive_u64(object, "total_genome_size", path)?;
    let composition_is_object = match object.get("composition") {
        Some(Value::Null) => false,
        Some(Value::Object(composition)) => {
            const BASES: [&str; 4] = ["A", "C", "G", "T"];
            const BASES_WITH_N: [&str; 5] = ["A", "C", "G", "T", "N"];
            if composition
                .keys()
                .any(|key| !BASES_WITH_N.contains(&key.as_str()))
            {
                return Err(effective_error(
                    path,
                    "campos desconocidos en reference.composition",
                ));
            }
            for base in BASES {
                let frequency = composition
                    .get(base)
                    .and_then(Value::as_f64)
                    .ok_or_else(|| {
                        effective_error(path, format!("reference.composition.{base} inválida"))
                    })?;
                if !(0.0..=1.0).contains(&frequency) {
                    return Err(effective_error(
                        path,
                        format!("reference.composition.{base} fuera de rango"),
                    ));
                }
            }
            let sum = BASES
                .iter()
                .map(|base| composition[*base].as_f64().expect("validated frequency"))
                .sum::<f64>();
            if !sum.is_finite() || (sum - 1.0).abs() > COMPOSITION_TOLERANCE {
                return Err(effective_error(
                    path,
                    "reference.composition no suma 1 dentro de la tolerancia",
                ));
            }
            if let Some(value) = composition.get("N") {
                let value = value
                    .as_f64()
                    .ok_or_else(|| effective_error(path, "reference.composition.N inválida"))?;
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    return Err(effective_error(
                        path,
                        "reference.composition.N fuera de rango",
                    ));
                }
            }
            true
        }
        _ => return Err(effective_error(path, "reference.composition inválida")),
    };
    if (source == "simulated") != composition_is_object {
        return Err(effective_error(
            path,
            "reference.composition no corresponde al origen declarado",
        ));
    }
    match (source, object.get("original_fasta_path")) {
        ("simulated", Some(Value::Null)) => {}
        ("supplied", Some(Value::String(value))) if !value.is_empty() => {}
        _ => {
            return Err(effective_error(
                path,
                "reference.original_fasta_path no corresponde al origen declarado",
            ));
        }
    }
    match (source, object.get("canonical_format")) {
        ("simulated", Some(Value::Object(format))) => {
            validate_canonical_fasta(format, path)?;
        }
        ("supplied", Some(Value::Null)) => {}
        _ => {
            return Err(effective_error(
                path,
                "reference.canonical_format no corresponde al origen declarado",
            ));
        }
    }
    let sequences = effective_array(object, "sequences", path)?;
    if sequences.is_empty() {
        return Err(effective_error(path, "reference.sequences está vacío"));
    }
    let mut sequence_ids = BTreeSet::new();
    let mut summed_length = 0_u64;
    for value in sequences {
        let sequence = value
            .as_object()
            .ok_or_else(|| effective_error(path, "una secuencia no es objeto"))?;
        const SEQUENCE_FIELDS: [&str; 7] = [
            "id",
            "class",
            "length",
            "header",
            "description",
            "alphabet_validation",
            "unknown_base_count",
        ];
        validate_effective_closed_object(
            sequence,
            &SEQUENCE_FIELDS,
            "reference.sequences[]",
            path,
        )?;
        let sequence_id = effective_non_empty_text(sequence, "id", path)?;
        if !sequence_ids.insert(sequence_id) {
            return Err(effective_error(path, "ID de secuencia repetido"));
        }
        match sequence.get("class").and_then(Value::as_str) {
            Some("chromosome" | "scaffold") => {}
            _ => return Err(effective_error(path, "clase de secuencia inválida")),
        }
        let length = effective_positive_u64(sequence, "length", path)?;
        summed_length = summed_length
            .checked_add(length)
            .ok_or_else(|| effective_error(path, "la suma de longitudes desborda u64"))?;
        let unknown = sequence
            .get("unknown_base_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| effective_error(path, "unknown_base_count inválido"))?;
        if unknown > length {
            return Err(effective_error(path, "unknown_base_count excede length"));
        }
        validate_reversible_text(
            sequence
                .get("header")
                .ok_or_else(|| effective_error(path, "header ausente"))?,
            "reference.sequences[].header",
            path,
        )?;
        match sequence.get("description") {
            Some(Value::Null) => {}
            Some(value) => {
                validate_reversible_text(value, "reference.sequences[].description", path)?
            }
            None => return Err(effective_error(path, "description ausente")),
        }
        if sequence.get("alphabet_validation").and_then(Value::as_str)
            != Some("ACGTN-case-insensitive")
        {
            return Err(effective_error(path, "alphabet_validation inválida"));
        }
    }
    if summed_length != declared_total {
        return Err(effective_error(
            path,
            "reference.total_genome_size no coincide con la suma de secuencias",
        ));
    }
    let fasta_id = effective_non_empty_text(object, "fasta_artifact_id", path)?.to_owned();
    require_effective_artifact_kind(catalog, &fasta_id, "reference_fasta", path)?;
    Ok((id, fasta_id))
}

fn validate_effective_reduction(
    value: &Value,
    catalog: &BTreeMap<&str, &str>,
    path: &Path,
) -> Result<EffectiveReductionLinks, Diagnostic> {
    let object = value
        .as_object()
        .ok_or_else(|| effective_error(path, "una reducción no es objeto"))?;
    const FIELDS: [&str; 11] = [
        "id",
        "source_reference_id",
        "region_size_bases",
        "target",
        "selection_subseed",
        "selection_order",
        "sequences",
        "state",
        "fasta_artifact_id",
        "map_artifact_id",
        "source_fasta_artifact_id",
    ];
    validate_effective_closed_object(object, &FIELDS, "reduced_references[]", path)?;
    let id = effective_non_empty_text(object, "id", path)?.to_owned();
    let source_reference_id =
        effective_non_empty_text(object, "source_reference_id", path)?.to_owned();
    effective_positive_u64(object, "region_size_bases", path)?;
    let target = effective_object(object, "target", path)?;
    const TARGET_FIELDS: [&str; 3] = ["mode", "declared", "bases"];
    validate_effective_closed_object(target, &TARGET_FIELDS, "reduced_references[].target", path)?;
    match target.get("mode").and_then(Value::as_str) {
        Some("target_bases" | "reduction_percent") => {}
        _ => return Err(effective_error(path, "target.mode inválido")),
    }
    effective_non_empty_text(target, "declared", path)?;
    let target_bases = effective_positive_u64(target, "bases", path)?;
    validate_lower_hex_64(
        effective_non_empty_text(object, "selection_subseed", path)?,
        "selection_subseed",
        path,
    )?;
    let selection_order = effective_array(object, "selection_order", path)?;
    if selection_order.is_empty()
        || selection_order
            .iter()
            .any(|value| value.as_u64().is_none_or(|value| value == 0))
    {
        return Err(effective_error(path, "selection_order inválido"));
    }
    validate_effective_derived_sequences(object, target_bases, path)?;
    if object.get("state").and_then(Value::as_str) != Some("complete") {
        return Err(effective_error(
            path,
            "estado de reducción efectiva inválido",
        ));
    }
    let fasta_artifact_id = effective_non_empty_text(object, "fasta_artifact_id", path)?.to_owned();
    let map_artifact_id = effective_non_empty_text(object, "map_artifact_id", path)?.to_owned();
    let source_fasta_artifact_id =
        effective_non_empty_text(object, "source_fasta_artifact_id", path)?.to_owned();
    require_effective_artifact_kind(catalog, &fasta_artifact_id, "reference_fasta", path)?;
    require_effective_artifact_kind(catalog, &map_artifact_id, "reference_reduction_map", path)?;
    require_effective_artifact_kind(catalog, &source_fasta_artifact_id, "reference_fasta", path)?;
    Ok(EffectiveReductionLinks {
        id,
        source_reference_id,
        fasta_artifact_id,
        map_artifact_id,
        source_fasta_artifact_id,
    })
}

fn validate_effective_derived_sequences(
    reduction: &serde_json::Map<String, Value>,
    target_bases: u64,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 7] = [
        "id",
        "class",
        "length",
        "header",
        "description",
        "alphabet_validation",
        "n_count",
    ];
    let sequences = effective_array(reduction, "sequences", path)?;
    if sequences.is_empty() {
        return Err(effective_error(path, "catálogo derivado vacío"));
    }
    let mut ids = BTreeSet::new();
    let mut total = 0_u64;
    for value in sequences {
        let sequence = value
            .as_object()
            .ok_or_else(|| effective_error(path, "una secuencia derivada no es objeto"))?;
        validate_effective_closed_object(
            sequence,
            &FIELDS,
            "reduced_references[].sequences[]",
            path,
        )?;
        let id = effective_non_empty_text(sequence, "id", path)?;
        if !ids.insert(id) {
            return Err(effective_error(path, "ID de secuencia derivada repetido"));
        }
        match sequence.get("class").and_then(Value::as_str) {
            Some("chromosome" | "scaffold") => {}
            _ => {
                return Err(effective_error(
                    path,
                    "clase de secuencia derivada inválida",
                ));
            }
        }
        let length = effective_positive_u64(sequence, "length", path)?;
        total = total
            .checked_add(length)
            .ok_or_else(|| effective_error(path, "la longitud derivada desborda u64"))?;
        let n_count = sequence
            .get("n_count")
            .and_then(Value::as_u64)
            .ok_or_else(|| effective_error(path, "n_count derivado inválido"))?;
        if n_count > length {
            return Err(effective_error(path, "n_count derivado excede length"));
        }
        validate_reversible_text(
            sequence
                .get("header")
                .ok_or_else(|| effective_error(path, "header derivado ausente"))?,
            "reduced_references[].sequences[].header",
            path,
        )?;
        match sequence.get("description") {
            Some(Value::Null) => {}
            Some(value) => validate_reversible_text(
                value,
                "reduced_references[].sequences[].description",
                path,
            )?,
            None => return Err(effective_error(path, "description derivada ausente")),
        }
        if sequence.get("alphabet_validation").and_then(Value::as_str)
            != Some("ACGTN-case-insensitive")
        {
            return Err(effective_error(
                path,
                "alphabet_validation derivada inválida",
            ));
        }
    }
    if total != target_bases {
        return Err(effective_error(
            path,
            "la suma del catálogo derivado no coincide con target.bases",
        ));
    }
    Ok(())
}

fn validate_effective_closed_object(
    object: &serde_json::Map<String, Value>,
    fields: &[&str],
    context: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if object.len() != fields.len() || object.keys().any(|key| !fields.contains(&key.as_str())) {
        return Err(effective_error(
            path,
            format!("campos ausentes o desconocidos en {context}"),
        ));
    }
    Ok(())
}

fn effective_object<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a serde_json::Map<String, Value>, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_object)
        .ok_or_else(|| effective_error(path, format!("{field} no es objeto")))
}

fn effective_array<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a [Value], Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .ok_or_else(|| effective_error(path, format!("{field} no es lista")))
}

fn effective_non_empty_text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty())
        .ok_or_else(|| effective_error(path, format!("{field} falta o está vacío")))
}

fn effective_positive_u64(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<u64, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_u64)
        .filter(|value| *value > 0)
        .ok_or_else(|| effective_error(path, format!("{field} inválido")))
}

fn effective_nonnegative_u64(
    object: &serde_json::Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<u64, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_u64)
        .ok_or_else(|| effective_error(path, format!("{field} inválido")))
}

fn reference_json(
    config: &EffectiveConfig,
    output: ReferenceOutput<'_>,
    composition: Option<Value>,
) -> Value {
    json!({
        "id": config.reference_id(),
        "source": match config.reference_source() {
            crate::application::EffectiveReferenceSource::Simulated => "simulated",
            crate::application::EffectiveReferenceSource::Supplied { .. } => "supplied",
        },
        "original_fasta_path": match config.reference_source() {
            crate::application::EffectiveReferenceSource::Simulated => Value::Null,
            crate::application::EffectiveReferenceSource::Supplied { original_path } => {
                Value::String(original_path.to_string_lossy().into_owned())
            },
        },
        "canonical_format": match config.reference_source() {
            crate::application::EffectiveReferenceSource::Simulated => canonical_fasta_json(),
            crate::application::EffectiveReferenceSource::Supplied { .. } => Value::Null,
        },
        "total_genome_size": config.total_genome_size(),
        "composition": composition,
        "sequences": output.reference.sequences().iter().map(|sequence| json!({
            "id": sequence.id().as_str(),
            "class": sequence_class(sequence.class()),
            "length": sequence.length(),
            "header": reversible_bytes_json(sequence.header().raw_bytes()),
            "description": reversible_description_json(sequence.header().description()),
            "alphabet_validation": alphabet_validation(sequence.alphabet_validation()),
            "unknown_base_count": sequence.unknown_base_count(),
        })).collect::<Vec<_>>(),
        "fasta_artifact_id": output.fasta.id,
    })
}

/// Lista de los cambios que la inserción aplicó, en el orden del plan.
fn insertion_changes_json(config: &EffectiveConfig) -> Result<Value, Diagnostic> {
    let plan = config.tract_insertion().ok_or_else(|| {
        Diagnostic::internal_invariant("results", "inserted_reference")
            .with("detail", "se publica una inserción sin plan de inserción")
    })?;
    Ok(Value::Array(
        plan.points()
            .iter()
            .map(|point| {
                json!({
                    "locus_id": point.locus_id(),
                    "sequence_id": point.sequence_id(),
                    "insert_after": point.insert_after(),
                    "inserted_bases": point.bases().len(),
                })
            })
            .collect(),
    ))
}

/// Sección `inserted_reference` de la configuración efectiva v6 (DTC-098): la
/// referencia **final**, sobre la que se simula, con su catálogo, los cambios
/// que la distinguen de la suministrada y sus artefactos.
pub fn inserted_reference_json(
    config: &EffectiveConfig,
    reference: &Reference,
    fasta: &ArtifactRecord,
    map: &ArtifactRecord,
    source_fasta: &ArtifactRecord,
) -> Result<Value, Diagnostic> {
    Ok(json!({
        "id": config.simulation_reference_id(),
        "source_reference_id": config.reference_id(),
        "sequences": reference.sequences().iter().map(|sequence| json!({
            "id": sequence.id().as_str(),
            "class": sequence_class(sequence.class()),
            "length": sequence.length(),
            "header": reversible_bytes_json(sequence.header().raw_bytes()),
            "description": reversible_description_json(sequence.header().description()),
            "alphabet_validation": alphabet_validation(sequence.alphabet_validation()),
            "n_count": sequence.unknown_base_count(),
        })).collect::<Vec<_>>(),
        "insertions": insertion_changes_json(config)?,
        "fasta_artifact_id": fasta.id,
        "map_artifact_id": map.id,
        "source_fasta_artifact_id": source_fasta.id,
    }))
}

/// Sección `reference_insertion` del manifiesto general v6 (DTC-098, REQ-MAN-07):
/// la referencia final, la de origen, los artefactos y los cambios aplicados.
pub fn manifest_insertion_json(
    config: &EffectiveConfig,
    state: UnitState,
    fasta: &ArtifactRecord,
    map: &ArtifactRecord,
    source_fasta: &ArtifactRecord,
) -> Result<Value, Diagnostic> {
    Ok(json!({
        "id": config.simulation_reference_id(),
        "source_reference_id": config.reference_id(),
        "state": state.as_str(),
        "fasta_artifact_id": fasta.id,
        "map_artifact_id": map.id,
        "source_fasta_artifact_id": source_fasta.id,
        "map_document_kind": crate::reference_tract_insertion::COORDINATE_MAP_DOCUMENT_KIND,
        "map_schema_version": crate::reference_tract_insertion::COORDINATE_MAP_SCHEMA_VERSION,
        "insertions": insertion_changes_json(config)?,
    }))
}

fn canonical_fasta_json() -> Value {
    json!({
        "version": CANONICAL_FASTA_VERSION,
        "line_base_count": CANONICAL_LINE_BASE_COUNT,
        "line_ending": "LF",
        "base_case": "uppercase",
    })
}

fn validate_canonical_fasta(
    object: &serde_json::Map<String, Value>,
    path: &Path,
) -> Result<(), Diagnostic> {
    const FIELDS: [&str; 4] = ["version", "line_base_count", "line_ending", "base_case"];
    validate_effective_closed_object(object, &FIELDS, "reference.canonical_format", path)?;
    if object.get("version").and_then(Value::as_str) != Some(CANONICAL_FASTA_VERSION)
        || object.get("line_base_count").and_then(Value::as_u64)
            != Some(CANONICAL_LINE_BASE_COUNT as u64)
        || object.get("line_ending").and_then(Value::as_str) != Some("LF")
        || object.get("base_case").and_then(Value::as_str) != Some("uppercase")
    {
        return Err(effective_error(path, "reference.canonical_format inválido"));
    }
    Ok(())
}

fn reversible_bytes_json(bytes: &[u8]) -> Value {
    match std::str::from_utf8(bytes) {
        Ok(text) => json!({"encoding": "utf-8", "value": text}),
        Err(_) => json!({"encoding": "bytes-hex", "value": encode_hex(bytes)}),
    }
}

fn reversible_description_json(description: Option<&ReversibleText>) -> Value {
    description.map_or(
        Value::Null,
        |text| json!({"encoding": text.encoding(), "value": text.encoded()}),
    )
}

fn validate_reversible_text(value: &Value, field: &str, path: &Path) -> Result<(), Diagnostic> {
    let object = value
        .as_object()
        .ok_or_else(|| effective_error(path, format!("{field} no es objeto")))?;
    const FIELDS: [&str; 2] = ["encoding", "value"];
    validate_effective_closed_object(object, &FIELDS, field, path)?;
    let encoded = object
        .get("value")
        .and_then(Value::as_str)
        .ok_or_else(|| effective_error(path, format!("{field}.value inválido")))?;
    if encoded.is_empty() {
        return Err(effective_error(path, format!("{field}.value está vacío")));
    }
    match object.get("encoding").and_then(Value::as_str) {
        Some("utf-8") => Ok(()),
        Some("bytes-hex")
            if encoded.len() % 2 == 0
                && encoded
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)) =>
        {
            Ok(())
        }
        _ => Err(effective_error(path, format!("{field}.encoding inválido"))),
    }
}

fn alphabet_validation(validation: AlphabetValidation) -> &'static str {
    match validation {
        AlphabetValidation::AcgtnCaseInsensitive => "ACGTN-case-insensitive",
    }
}

fn encode_hex(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(char::from(HEX[usize::from(byte >> 4)]));
        output.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    output
}

fn manifest_reference_json(output: &ManifestReferenceLink) -> Value {
    json!({
        "id": output.id,
        "fasta_artifact_id": output.fasta_artifact_id,
    })
}

fn manifest_reduction_json(output: &ManifestReductionLink) -> Value {
    json!({
        "id": output.id,
        "source_reference_id": output.source_reference_id,
        "eligible_bases": output.eligible_bases,
        "region_size_bases": output.region_size_bases,
        "target": {
            "mode": output.target_mode,
            "declared": output.target_declared,
            "bases": output.target_bases,
        },
        "root_seed": output.root_seed,
        "selection_subseed": output.selection_subseed,
        "selection_order": output.selection_order,
        "state": output.state.as_str(),
        "fasta_artifact_id": output.fasta_artifact_id,
        "map_artifact_id": output.map_artifact_id,
        "source_fasta_artifact_id": output.source_fasta_artifact_id,
        "map_document_kind": output.map_document_kind,
        "map_schema_version": output.map_schema_version,
    })
}

fn reduction_json(output: &ReducedReferenceOutput<'_>) -> Value {
    json!({
        "id": output.effective.id(),
        "source_reference_id": output.effective.source_reference_id(),
        "region_size_bases": output.effective.region_size_bases(),
        "target": {
            "mode": output.effective.target().mode(),
            "declared": output.effective.target().declared_value(),
            "bases": output.effective.target_bases(),
        },
        "selection_subseed": output.effective.selection_subseed().to_hex(),
        "selection_order": output.plan.segments().iter().map(|segment| segment.selection_order()).collect::<Vec<_>>(),
        "sequences": output.reference.sequences().iter().map(|sequence| json!({
            "id": sequence.id().as_str(),
            "class": sequence_class(sequence.class()),
            "length": sequence.length(),
            "header": reversible_bytes_json(sequence.header().raw_bytes()),
            "description": reversible_description_json(sequence.header().description()),
            "alphabet_validation": alphabet_validation(sequence.alphabet_validation()),
            "n_count": sequence.unknown_base_count(),
        })).collect::<Vec<_>>(),
        "state": output.state.as_str(),
        "fasta_artifact_id": output.fasta.id,
        "map_artifact_id": output.map.id,
        "source_fasta_artifact_id": output.source_fasta.id,
    })
}

fn output_categories(manifest: &GeneralManifest<'_>) -> Value {
    let requested = manifest
        .reduced_references
        .iter()
        .flat_map(|reduction| {
            [
                reduction.fasta_artifact_id.as_deref(),
                reduction.map_artifact_id.as_deref(),
            ]
        })
        .flatten()
        .chain(manifest.vcf_outputs.iter().flat_map(|output| {
            std::iter::once(output.vcf_artifact_id.as_str())
                .chain(output.index_artifact_id.as_deref())
        }))
        .collect::<BTreeSet<_>>();
    let minimum = manifest
        .artifacts
        .iter()
        .filter(|artifact| !requested.contains(artifact.id.as_str()))
        .map(|artifact| artifact.id.as_str())
        .collect::<Vec<_>>();
    json!({
        "minimum": minimum,
        "defaulted": [],
        "requested": requested,
        "disabled": [],
    })
}

fn validate_reference_output(
    config: &EffectiveConfig,
    output: ReferenceOutput<'_>,
) -> Result<(), Diagnostic> {
    if output.reference.authoritative_resource().is_none()
        || output
            .reference
            .authoritative_resource()
            .and_then(Path::to_str)
            != Some(&output.fasta.path)
    {
        return Err(Diagnostic::internal_invariant(
            "results",
            "reference artifact is not the authoritative resource",
        ));
    }
    let origin_matches = match (config.reference_source(), output.reference.origin()) {
        (EffectiveReferenceSource::Simulated, ReferenceOrigin::Generated) => true,
        (
            EffectiveReferenceSource::Supplied { original_path },
            ReferenceOrigin::Supplied {
                original_path: observed,
            },
        ) => original_path == observed,
        // **Aviso para el paso que conecte la inserción de tractos.** El origen
        // `ReferenceOrigin::TractInsertion` caerá aquí y la comparación dirá que
        // no coincide, porque la fuente efectiva sigue siendo la suministrada
        // mientras la referencia publicada ya es la derivada. Hoy no hay error
        // porque ningún camino produce todavía esa referencia; en cuanto lo
        // haya, esta rama necesita su propio caso y la fuente efectiva necesita
        // saber que la referencia se derivó.
        _ => false,
    };
    if !origin_matches
        || output.reference.id().as_str() != config.reference_id()
        || output.reference.total_length() != config.total_genome_size()
        || output.reference.sequences().len() != config.sequences().len()
        || output
            .reference
            .sequences()
            .iter()
            .zip(config.sequences())
            .any(|(observed, expected)| {
                observed.id().as_str() != expected.id()
                    || observed.class() != expected.class()
                    || observed.length() != expected.length()
            })
    {
        return Err(Diagnostic::internal_invariant(
            "results",
            "reference domain does not match effective configuration",
        ));
    }
    Ok(())
}

fn sequence_class(class: SequenceClass) -> &'static str {
    match class {
        SequenceClass::Chromosome => "chromosome",
        SequenceClass::Scaffold => "scaffold",
    }
}

pub(crate) fn validate_relative_link(path: &str) -> Result<(), Diagnostic> {
    let candidate = Path::new(path);
    if path.is_empty()
        || path == "INCOMPLETE"
        || candidate.is_absolute()
        || path.contains("//")
        || path.contains('\\')
        || path.ends_with('/')
        || path.bytes().any(|byte| byte <= b'\x1f' || byte == b'\x7f')
        || path
            .split('/')
            .any(|component| matches!(component, "." | ".."))
        || candidate
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(manifest_error(
            Path::new("manifest.json"),
            format!("la ruta de artefacto no es canónica: {path:?}"),
        ));
    }
    Ok(())
}

fn effective_schema_version(path: &Path, bytes: &[u8]) -> Result<&'static str, Diagnostic> {
    let document: Value = serde_json::from_slice(bytes)
        .map_err(|error| effective_error(path, format!("JSON no legible: {error}")))?;
    let object = document
        .as_object()
        .ok_or_else(|| effective_error(path, "la raíz no es un objeto"))?;
    match object.get("schema_version").and_then(Value::as_str) {
        Some(EFFECTIVE_CONFIG_SCHEMA_VERSION) => Ok(EFFECTIVE_CONFIG_SCHEMA_VERSION),
        Some(EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS) => {
            Ok(EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS)
        }
        observed => Err(effective_error_for_schema(
            path,
            "effective-config-v4.schema.json",
            format!(
                "versión de configuración efectiva no admitida: {:?}; se admiten \"{}\" y \"{}\"",
                observed,
                EFFECTIVE_CONFIG_SCHEMA_VERSION,
                EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS
            ),
        )),
    }
}

/// Nombre del esquema que le corresponde a una versión declarada.
///
/// Sale del propio tipo que decide qué versión se publica, de modo que el nombre
/// no pueda desviarse de la versión emitida.
fn effective_schema_file_for(version: Option<&str>) -> &'static str {
    if version == Some(EFFECTIVE_CONFIG_SCHEMA_VERSION_WITH_REPEATS) {
        DeclaredDocumentVersion::V4.effective_config_schema_file()
    } else {
        DeclaredDocumentVersion::V2.effective_config_schema_file()
    }
}

fn retag_format_schema(error: Diagnostic, schema: &'static str) -> Diagnostic {
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

#[allow(dead_code)]
fn retag_v1_pair_error(
    error: Diagnostic,
    effective_path: &Path,
    manifest_path: &Path,
) -> Diagnostic {
    if error.code() != "GPS-FMT-006" {
        return error;
    }
    let context = error.context();
    let error_path = context.get("path").map(String::as_str);
    let error_schema = context.get("schema").map(String::as_str);
    if error_path == Some(&manifest_path.display().to_string())
        || error_schema.is_some_and(|schema| schema.starts_with("manifest-general-"))
    {
        retag_format_schema(error, "manifest-general-v1.schema.json")
    } else if error_path == Some(&effective_path.display().to_string())
        || error_schema.is_some_and(|schema| schema.starts_with("effective-config-"))
    {
        retag_format_schema(error, "effective-config-v1.schema.json")
    } else {
        error
    }
}

/// Envoltorio puro sobre [`manifest_error_for_schema`]. Lleva
/// `#[track_caller]` para que el barrido de guardas ([`crate::guard_trace`])
/// anote el punto de rechazo real y no esta línea.
#[track_caller]
fn manifest_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    manifest_error_for_schema(path, "manifest-general-v4.schema.json", detail)
}

/// Las dos versiones que un manifiesto general de este recorrido puede
/// declarar: la de siempre y la que añade las repeticiones en tándem.
#[must_use]
pub fn is_supported_general_manifest_version(version: Option<&str>) -> bool {
    matches!(
        version,
        Some(GENERAL_MANIFEST_SCHEMA_VERSION) | Some(GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS)
    )
}

/// El esquema que hay que citar en un diagnóstico, según lo que el documento
/// declare. Denunciar siempre el de la versión 4 mandaba a mirar un fichero
/// que no era el suyo.
fn general_manifest_schema_file_for(version: Option<&str>) -> &'static str {
    if version == Some(GENERAL_MANIFEST_SCHEMA_VERSION_WITH_REPEATS) {
        "manifest-general-v6.schema.json"
    } else {
        "manifest-general-v4.schema.json"
    }
}

/// Constructora de los rechazos del **manifiesto general**, el documento que
/// resume la ejecución entera.
///
/// La anotación `#[track_caller]` y la llamada a `record` son el barrido de
/// guardas descrito en [`crate::guard_trace`]: apagado no cuesta nada, y
/// encendido anota **desde qué línea** se rechazó, no esta.
#[track_caller]
fn manifest_error_for_schema(
    path: &Path,
    schema: &'static str,
    detail: impl Into<String>,
) -> Diagnostic {
    crate::guard_trace::record();
    manifest_diagnostic(path, schema, detail)
}

/// El diagnóstico de un rechazo del manifiesto general, **sin anotar nada**.
///
/// Lo construyen [`manifest_error_for_schema`], que sí anota el punto de
/// rechazo, y los envoltorios que reciben un motivo ya anotado por
/// [`reason`]: ahí el rechazo nació en otra línea, y anotar también la del
/// envoltorio contaría como punto de rechazo una línea que no decide nada.
fn manifest_diagnostic(path: &Path, schema: &'static str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path.display().to_string())
    .with("schema", schema)
    .with("detail", detail.into())
}

/// Envoltorio puro sobre [`effective_error_for_schema`]. Lleva
/// `#[track_caller]` para que el barrido de guardas ([`crate::guard_trace`])
/// anote el punto de rechazo real y no esta línea.
#[track_caller]
fn effective_error(path: &Path, detail: impl Into<String>) -> Diagnostic {
    effective_error_for_schema(path, "effective-config-v4.schema.json", detail)
}

/// La anotación `#[track_caller]` y la llamada a `record` son el barrido de
/// guardas descrito en [`crate::guard_trace`]: apagado no cuesta nada, y
/// encendido anota **desde qué línea** se rechazó, no esta.
#[track_caller]
fn effective_error_for_schema(
    path: &Path,
    schema: &'static str,
    detail: impl Into<String>,
) -> Diagnostic {
    crate::guard_trace::record();
    effective_diagnostic(path, schema, detail)
}

/// El diagnóstico de un rechazo de la configuración efectiva, **sin anotar nada**.
///
/// Lo construyen [`effective_error_for_schema`], que sí anota el punto de
/// rechazo, y los envoltorios que reciben un motivo ya anotado por
/// [`reason`]: ahí el rechazo nació en otra línea, y anotar también la del
/// envoltorio contaría como punto de rechazo una línea que no decide nada.
fn effective_diagnostic(
    path: &Path,
    schema: &'static str,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-006",
        "Un documento JSON producido no satisface su esquema publicado.",
    )
    .with("path", path.display().to_string())
    .with("schema", schema)
    .with("detail", detail.into())
}

/// El motivo de un rechazo que viaja como **texto** hasta que un envoltorio lo
/// convierte en diagnóstico.
///
/// Los validadores de este módulo que devuelven `Result<_, &'static str>` o
/// `Result<_, String>` no construyen un [`Diagnostic`]: lo hace quien los llama,
/// con **una sola** llamada a [`manifest_error`] o [`effective_error`] que
/// recoge todos sus rechazos. Para la puerta de cobertura de guardas
/// ([`crate::guard_trace`]) eran entonces un único punto de rechazo, y bastaba
/// con que se disparase uno para dar por ejercitados los demás.
///
/// Cada rechazo de texto se escribe ahora `Err(reason("..."))`. `#[track_caller]`
/// hace que se anote **la línea del rechazo**, no esta; el texto no cambia ni
/// un byte. Hay que llamarla **solo en la rama que rechaza**: por eso se usa
/// `ok_or_else(|| reason(..))` y no `ok_or(reason(..))`, que anotaría siempre.
#[track_caller]
fn reason<T>(detail: T) -> T {
    crate::guard_trace::record();
    detail
}

/// [`manifest_error`] para un motivo que ya anotó [`reason`].
fn manifest_error_from_reason(path: &Path, detail: &str) -> Diagnostic {
    manifest_diagnostic(path, "manifest-general-v4.schema.json", detail)
}

/// [`effective_error`] para un motivo que ya anotó [`reason`].
fn effective_error_from_reason(path: &Path, detail: &str) -> Diagnostic {
    effective_diagnostic(path, "effective-config-v4.schema.json", detail)
}

fn json_error(error: serde_json::Error) -> Diagnostic {
    Diagnostic::internal_invariant("results", format!("JSON serialization failed: {error}"))
}

/// Constructora del rechazo `GPS-IO-002`: un artefacto publicado no tiene el
/// tamaño o el SHA-256 que su registro declara. Devuelve el diagnóstico base y
/// el llamante le añade el contexto que le toca. Lleva `#[track_caller]` para
/// que la puerta de guardas anote el rechazo real y no esta línea.
#[track_caller]
fn artifact_identity_error() -> Diagnostic {
    crate::guard_trace::record();
    Diagnostic::io(
        "GPS-IO-002",
        "El tamaño o SHA-256 de un artefacto no coincide con el valor esperado.",
    )
}

fn io_error(path: &Path, operation: &str, error: std::io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}
