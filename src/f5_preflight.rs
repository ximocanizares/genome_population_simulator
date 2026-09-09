//! Preflight puro de recursos F5 (DTC-089).
//!
//! Este módulo no abre un espacio de publicación ni inicializa un generador
//! aleatorio. Valida SHA-256, ID y sintaxis mientras recorre el FASTA una vez
//! junto con la máscara opcional; no construye un `Reference` ni conserva las
//! bases completas antes de admitir los presupuestos. Solo retiene metadatos,
//! el intervalo de máscara vigente y los tramos elegibles.

#[cfg(test)]
use std::cell::Cell;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;

use sha2::{Digest, Sha256};

use crate::application::EffectiveF5Config;
use crate::configuration::{DeclaredF5Target, DeclaredNumber};
use crate::diagnostics::Diagnostic;
use crate::exact_decimal::ExactDecimal;
pub use crate::f5_admission::F5_RESOURCE_ACCOUNTING_ALGORITHM;
use crate::f5_admission::{
    AdmissionPhase, AdmittedInputIdentity, F5AdmissionReceipt, F5ResourceLedger, InputRole,
    LedgerCharge, OperationCost, OutputReservation, PhysicalInputIdentity, WorkingMemoryPeak,
};
use crate::identity::{ReferenceId, SequenceId};
use crate::reference_fasta::{empty_fasta, empty_header, input_error, invalid_symbol};

pub const F5_HARD_OPERATION_LIMIT: u64 = 1 << 28;
pub const F5_HARD_WORKING_MEMORY_LIMIT: u64 = 1 << 30;
pub const F5_HARD_OUTPUT_LIMIT: u64 = 1 << 31;

pub const F5_DEFAULT_MAX_LAYOUT_ATTEMPTS: u64 = 1_000;
const DEFAULT_MAX_HAPLOTYPE_ATTEMPTS: u64 = 100;
const DEFAULT_MAX_SWAP_PROPOSALS: u64 = 100_000;
const DEFAULT_MAX_HAPLOTYPE_CELLS: u64 = 16_777_216;
const F5_SCAN_FIXED_MEMORY_BYTES: u64 = 65_536;
const F5_SCAN_RETAINED_ITEM_BYTES: u64 = 256;

// El observador es exclusivamente instrumental: no forma parte de la
// admisión ni existe en el binario de producción. F5.1.2 lo usa para probar
// que la visita de capacidad presupuestada por v4 ocurre en el recorrido real.
#[cfg(test)]
thread_local! {
    static CAPACITY_TRACT_VISITS: Cell<u64> = const { Cell::new(0) };
}

#[cfg(test)]
pub(crate) fn reset_test_capacity_tract_visits() {
    CAPACITY_TRACT_VISITS.with(|visits| visits.set(0));
}

#[cfg(test)]
pub(crate) fn test_capacity_tract_visits() -> u64 {
    CAPACITY_TRACT_VISITS.with(Cell::get)
}

fn record_capacity_tract_visit() {
    #[cfg(test)]
    CAPACITY_TRACT_VISITS.with(|visits| visits.set(visits.get().saturating_add(1)));
}

/// Aportaciones obligatorias o activadas por la misma ejecución pero ajenas a
/// F5. Se mantienen separadas para impedir mezclar unidades en la fórmula.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct F5OtherResources {
    operation_units: u64,
    working_memory_bytes: u64,
    output_bytes: u64,
}

impl F5OtherResources {
    #[must_use]
    pub const fn new(operation_units: u64, working_memory_bytes: u64, output_bytes: u64) -> Self {
        Self {
            operation_units,
            working_memory_bytes,
            output_bytes,
        }
    }
    #[must_use]
    pub const fn operation_units(self) -> u64 {
        self.operation_units
    }
    #[must_use]
    pub const fn working_memory_bytes(self) -> u64 {
        self.working_memory_bytes
    }
    #[must_use]
    pub const fn output_bytes(self) -> u64 {
        self.output_bytes
    }
}

/// Tramo maximal de A/C/G/T no excluido. Las coordenadas son inclusivas y
/// 1-based; la lista es transitoria y nunca forma parte de la salida F5.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F5EligibleTract {
    sequence_ordinal: u64,
    start_1based: u64,
    end_1based: u64,
}

impl F5EligibleTract {
    /// Constructor exclusivo de pruebas; producción solo recibe tramos del
    /// escáner F5.0.4, que conserva sus invariantes de coordenadas.
    #[cfg(test)]
    #[must_use]
    pub(crate) const fn new(sequence_ordinal: u64, start_1based: u64, end_1based: u64) -> Self {
        Self {
            sequence_ordinal,
            start_1based,
            end_1based,
        }
    }
    #[must_use]
    pub const fn sequence_ordinal(self) -> u64 {
        self.sequence_ordinal
    }
    #[must_use]
    pub const fn start_1based(self) -> u64 {
        self.start_1based
    }
    #[must_use]
    pub const fn end_1based(self) -> u64 {
        self.end_1based
    }
    #[must_use]
    pub const fn len(self) -> u64 {
        self.end_1based - self.start_1based + 1
    }
    #[must_use]
    pub const fn is_empty(self) -> bool {
        false
    }
}

/// Objetivo F5 declarado y resuelto que vincula el preflight con su consumidor.
///
/// La fracción se conserva como racional reducido, no como `f64`: dos
/// literales que denotan el mismo racional son el mismo objetivo, pero nunca
/// se confunden con `covered_bases` aunque den el mismo T en una referencia.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum F5ResolvedTarget {
    BlockCount { declared_block_count: u64 },
    CoveredBases { declared_covered_bases: u64 },
    CoverageFraction { numerator: u64, denominator: u64 },
}

/// Identidad canónica de todas las entradas efectivas que condicionan la
/// admisión o la ejecución F5. No depende de la ubicación de los documentos
/// de autoría, pero sí de la presencia, ruta resuelta e identidad declarada de
/// la máscara.
///
/// La codificación usa campos con longitud y órdenes canónicos de población y
/// override; no concatena literales TOML ni ubicaciones de autoría. Así dos
/// configuraciones semánticamente iguales comparten identidad, mientras que
/// semilla/origen, fundadores, F, modo o tolerancia no pueden cruzar un
/// preflight ya admitido.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct F5PreflightConfigurationIdentity {
    canonical_bytes: Vec<u8>,
}

impl F5PreflightConfigurationIdentity {
    fn from_config(config: &EffectiveF5Config) -> Result<Self, Diagnostic> {
        let scenario = config.declared_scenario();
        let mut canonical = CanonicalF5Identity::new();
        canonical.field_text("scenario_id", config.scenario_id());
        canonical.field_text("reference_id", config.reference_id());
        canonical.field_text("reference_sha256", config.reference_sha256());
        canonical.field_text("profile_id", config.profile().profile_id());
        canonical.field_text(
            "profile_metadata_sha256",
            scenario.profile_metadata_sha256.value(),
        );
        canonical.field_text("seed", config.seed());
        canonical.field_text("seed_origin", config.seed_origin());
        canonical.field_u64("ploidy", u64::from(config.ploidy()));
        canonical.field_u64("founder_count", *scenario.founder_count.value());
        canonical.field_f64_bits(
            "scenario_inbreeding_f",
            Some(
                canonical_f5_number(config, "f5.inbreeding_f", scenario.inbreeding_f.as_ref())?
                    .unwrap_or(0.0_f64.to_bits()),
            ),
        );

        match (
            config.excluded_regions_path(),
            scenario.excluded_regions_sha256.as_ref(),
        ) {
            (None, None) => canonical.field_none("excluded_regions"),
            (Some(path), Some(digest)) => {
                canonical.field_text("excluded_regions_path", &path.display().to_string());
                canonical.field_text("excluded_regions_sha256", digest.value());
            }
            _ => {
                return Err(f5_configuration(
                    config,
                    "excluded_regions",
                    "la presencia de la ruta y del SHA-256 de máscara debe coincidir",
                ));
            }
        }

        match resolved_target_for_config(config)? {
            F5ResolvedTarget::BlockCount {
                declared_block_count,
            } => {
                canonical.field_text("target_kind", "block_count");
                canonical.field_u64("target_block_count", declared_block_count);
            }
            F5ResolvedTarget::CoveredBases {
                declared_covered_bases,
            } => {
                canonical.field_text("target_kind", "covered_bases");
                canonical.field_u64("target_covered_bases", declared_covered_bases);
            }
            F5ResolvedTarget::CoverageFraction {
                numerator,
                denominator,
            } => {
                canonical.field_text("target_kind", "coverage_fraction");
                canonical.field_u64("target_fraction_numerator", numerator);
                canonical.field_u64("target_fraction_denominator", denominator);
            }
        }

        let mut populations = scenario.subpopulations.iter().collect::<Vec<_>>();
        populations.sort_unstable_by_key(|population| population.population_id.value());
        canonical.field_u64(
            "subpopulation_count",
            u64::try_from(populations.len()).map_err(|_| {
                f5_configuration(config, "subpopulations", "el recuento no cabe en u64")
            })?,
        );
        for population in populations {
            canonical.field_text(
                "subpopulation.population_id",
                population.population_id.value(),
            );
            canonical.field_u64(
                "subpopulation.individual_count",
                *population.individual_count.value(),
            );
            canonical.field_text(
                "subpopulation.sampling",
                population
                    .sampling
                    .as_ref()
                    .map_or("probabilistic", |mode| mode.value()),
            );
            canonical.field_f64_bits(
                "subpopulation.inbreeding_f",
                Some(
                    canonical_f5_number(
                        config,
                        "subpopulations.inbreeding_f",
                        population
                            .inbreeding_f
                            .as_ref()
                            .or(scenario.inbreeding_f.as_ref()),
                    )?
                    .unwrap_or(0.0_f64.to_bits()),
                ),
            );
            match (&population.founder_ids, &population.founder_keys) {
                (Some(founders), None) => {
                    canonical.field_text("subpopulation.founder_kind", "founder_ids");
                    canonical.field_u64(
                        "subpopulation.founder_count",
                        u64::try_from(founders.len()).map_err(|_| {
                            f5_configuration(
                                config,
                                "subpopulations.founder_ids",
                                "el recuento no cabe en u64",
                            )
                        })?,
                    );
                    for founder in founders {
                        canonical.field_text("subpopulation.founder", founder.value());
                    }
                }
                (None, Some(founders)) => {
                    canonical.field_text("subpopulation.founder_kind", "founder_keys");
                    canonical.field_u64(
                        "subpopulation.founder_count",
                        u64::try_from(founders.len()).map_err(|_| {
                            f5_configuration(
                                config,
                                "subpopulations.founder_keys",
                                "el recuento no cabe en u64",
                            )
                        })?,
                    );
                    for founder in founders {
                        canonical.field_text("subpopulation.founder", founder.value());
                    }
                }
                _ => {
                    return Err(f5_configuration(
                        config,
                        "subpopulations.founder_ids,subpopulations.founder_keys",
                        "la configuración efectiva no conserva exactamente un conjunto de fundadores",
                    ));
                }
            }
        }

        let mut overrides = scenario.sampling_overrides.iter().collect::<Vec<_>>();
        overrides.sort_unstable_by_key(|override_| {
            (
                override_.population_id.value(),
                *override_.block_index.value(),
            )
        });
        canonical.field_u64(
            "sampling_override_count",
            u64::try_from(overrides.len()).map_err(|_| {
                f5_configuration(config, "sampling_overrides", "el recuento no cabe en u64")
            })?,
        );
        for override_ in overrides {
            canonical.field_text(
                "sampling_override.population_id",
                override_.population_id.value(),
            );
            canonical.field_u64(
                "sampling_override.block_index",
                *override_.block_index.value(),
            );
            canonical.field_text("sampling_override.mode", override_.mode.value());
            let population = scenario
                .subpopulations
                .iter()
                .find(|population| {
                    population.population_id.value() == override_.population_id.value()
                })
                .ok_or_else(|| {
                    f5_configuration(
                        config,
                        "sampling_overrides.population_id",
                        "la configuración efectiva no conserva la población sustituida",
                    )
                })?;
            let default_tolerance = 0.5
                / (population
                    .individual_count
                    .value()
                    .checked_mul(u64::from(config.ploidy()))
                    .ok_or_else(|| {
                        f5_configuration(
                            config,
                            "sampling_overrides.controlled_max_abs_error",
                            "las copias efectivas desbordan u64",
                        )
                    })? as f64);
            let tolerance = if override_.mode.value() == "controlled" {
                canonical_f5_number(
                    config,
                    "sampling_overrides.controlled_max_abs_error",
                    override_.controlled_max_abs_error.as_ref(),
                )?
                .unwrap_or(default_tolerance.to_bits())
            } else {
                0_u64
            };
            canonical.field_f64_bits(
                "sampling_override.controlled_max_abs_error",
                Some(tolerance),
            );
        }

        canonical.field_u64(
            "max_layout_attempts",
            limit(
                scenario.max_layout_attempts.as_ref(),
                F5_DEFAULT_MAX_LAYOUT_ATTEMPTS,
            ),
        );
        canonical.field_u64(
            "max_haplotype_attempts",
            limit(
                scenario.max_haplotype_attempts.as_ref(),
                DEFAULT_MAX_HAPLOTYPE_ATTEMPTS,
            ),
        );
        canonical.field_u64(
            "max_swap_proposals",
            limit(
                scenario.max_swap_proposals.as_ref(),
                DEFAULT_MAX_SWAP_PROPOSALS,
            ),
        );
        canonical.field_u64(
            "max_haplotype_cells",
            limit(
                scenario.max_haplotype_cells.as_ref(),
                DEFAULT_MAX_HAPLOTYPE_CELLS,
            ),
        );
        canonical.field_u64(
            "operation_limit",
            limit(
                scenario.max_operation_units.as_ref(),
                F5_HARD_OPERATION_LIMIT,
            ),
        );
        canonical.field_u64(
            "working_memory_limit_bytes",
            limit(
                scenario.max_working_memory_bytes.as_ref(),
                F5_HARD_WORKING_MEMORY_LIMIT,
            ),
        );
        canonical.field_u64(
            "output_limit_bytes",
            limit(scenario.max_output_bytes.as_ref(), F5_HARD_OUTPUT_LIMIT),
        );
        Ok(Self {
            canonical_bytes: canonical.finish(),
        })
    }
}

/// Codificador sin ambigüedades: cada etiqueta y valor llevan su longitud
/// big-endian. Nunca se serializan posiciones, espacios TOML ni rutas de los
/// dos documentos de autoría.
struct CanonicalF5Identity {
    bytes: Vec<u8>,
}

impl CanonicalF5Identity {
    fn new() -> Self {
        Self {
            bytes: b"f5-preflight-configuration-v1".to_vec(),
        }
    }

    fn field(&mut self, tag: &str, value: &[u8]) {
        self.bytes.extend_from_slice(
            &u64::try_from(tag.len())
                .expect("las etiquetas de identidad F5 caben en u64")
                .to_be_bytes(),
        );
        self.bytes.extend_from_slice(tag.as_bytes());
        self.bytes.extend_from_slice(
            &u64::try_from(value.len())
                .expect("los valores de identidad F5 caben en u64")
                .to_be_bytes(),
        );
        self.bytes.extend_from_slice(value);
    }

    fn field_text(&mut self, tag: &str, value: &str) {
        self.field(tag, value.as_bytes());
    }

    fn field_u64(&mut self, tag: &str, value: u64) {
        self.field(tag, &value.to_be_bytes());
    }

    fn field_f64_bits(&mut self, tag: &str, value: Option<u64>) {
        match value {
            Some(bits) => self.field(tag, &bits.to_be_bytes()),
            None => self.field_none(tag),
        }
    }

    fn field_none(&mut self, tag: &str) {
        self.field(tag, b"none");
    }

    fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

fn canonical_f5_number(
    config: &EffectiveF5Config,
    field: &str,
    value: Option<&DeclaredNumber>,
) -> Result<Option<u64>, Diagnostic> {
    value
        .map(|value| {
            value
                .literal()
                .parse::<f64>()
                .map(f64::to_bits)
                .map_err(|_| {
                    f5_configuration(
                        config,
                        field,
                        "el valor efectivo no puede codificarse como binary64",
                    )
                })
        })
        .transpose()
}

/// Resultado admitido —todavía no sellado— del preflight
/// `f5-resource-accounting-v7`.
///
/// Es contabilidad interna: ningún consumidor puede derivar layout ni contrato
/// de él hasta que `seal` lo une a su `F5AdmissionReceipt`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct F5ResourcePreflight {
    configuration_identity: F5PreflightConfigurationIdentity,
    scenario_id: String,
    reference_id: String,
    reference_sha256: String,
    profile_id: String,
    profile_metadata_sha256: String,
    block_length_minimum: u64,
    block_length_maximum: u64,
    // Compartir el único backing `Vec` permite devolver el resultado público
    // sin una segunda asignación proporcional a la fragmentación.
    eligible_tracts: Arc<Vec<F5EligibleTract>>,
    eligible_bases: u64,
    reference_bases_processed: u64,
    fasta_bytes_processed: u64,
    mask_bytes_processed: u64,
    general_toml_bytes: u64,
    scenario_toml_bytes: u64,
    profile_metadata_bytes: u64,
    profile_tsv_bytes: u64,
    input_parser_peak_bytes: u64,
    tract_count: u64,
    retained_sequence_ids: u64,
    gap_support_cardinality: u64,
    gap_selector_operation_units: u64,
    site_bound_per_block: u64,
    bin_count_total: u64,
    profile_rows: u64,
    normalization_rows: u64,
    tsv_rows_bound: u64,
    json_rows_bound: u64,
    other: F5OtherResources,
    operation_units: u64,
    working_memory_bytes: u64,
    output_bytes: u64,
    operation_limit: u64,
    working_memory_limit_bytes: u64,
    output_limit_bytes: u64,
    max_layout_attempts: u64,
    block_count: u64,
    target_bases: Option<u64>,
    resolved_target: F5ResolvedTarget,
    /// Identidades de la primera pasada, todavía no selladas ni publicadas.
    physical_input_identities: Vec<AdmittedInputIdentity>,
}

impl F5ResourcePreflight {
    /// Comprueba el enlace semántico completo con la configuración efectiva,
    /// sin incluir aportaciones ajenas que el consumidor debe entregar por
    /// separado mediante `F5LayoutRequest::new_with_other_resources`.
    pub(crate) fn matches_configuration(&self, config: &EffectiveF5Config) -> bool {
        F5PreflightConfigurationIdentity::from_config(config)
            .is_ok_and(|identity| identity == self.configuration_identity)
    }

    /// Comprueba también las contribuciones externas incorporadas a las tres
    /// cotas de recursos del preflight.
    pub(crate) fn matches_admission(
        &self,
        config: &EffectiveF5Config,
        other: F5OtherResources,
    ) -> bool {
        self.matches_configuration(config) && self.other == other
    }

    /// Completa el binding semántico después de resolver la semilla. La
    /// entropía solo puede solicitarse tras admitir físicamente todas las
    /// entradas, por lo que este paso no recalcula ningún coste ni reabre E/S.
    pub(crate) fn bind_configuration(
        &mut self,
        config: &EffectiveF5Config,
    ) -> Result<(), Diagnostic> {
        self.configuration_identity = F5PreflightConfigurationIdentity::from_config(config)?;
        Ok(())
    }

    /// Sella la admisión uniendo la proyección v7 con su receipt. Es la única
    /// transición de admitido a sellado y no recalcula ningún coste.
    pub(crate) fn seal(self, receipt: F5AdmissionReceipt) -> SealedF5Admission {
        SealedF5Admission {
            receipt,
            preflight: self,
        }
    }
}

/// Admisión F5 sellada: la única forma de la que un consumidor puede derivar
/// layout, publicación o proyección de contrato.
///
/// Une el receipt —identidades físicas y contabilidad por fase— con la
/// proyección de recursos v7 que las mismas lecturas produjeron. No existe
/// constructor público: solo la máquina de admisión completa puede sellarla,
/// y ningún receipt v5/v6 puede aceptarse ni migrarse a esta forma.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SealedF5Admission {
    receipt: F5AdmissionReceipt,
    preflight: F5ResourcePreflight,
}

impl SealedF5Admission {
    #[must_use]
    pub const fn algorithm(&self) -> &'static str {
        F5_RESOURCE_ACCOUNTING_ALGORITHM
    }

    #[must_use]
    pub const fn receipt(&self) -> &F5AdmissionReceipt {
        &self.receipt
    }

    /// Proyección exacta de `resource_preflight` que efectivo y manifiesto v5
    /// copian escalar y posicionalmente: algoritmo v7, identidades físicas en
    /// orden de rol contractual, fases en el orden cerrado de la máquina y los
    /// veintisiete escalares del algoritmo.
    #[must_use]
    pub fn to_contract_json(&self) -> serde_json::Value {
        let preflight = &self.preflight;
        serde_json::json!({
            "algorithm": self.algorithm(),
            "input_identities": self
                .receipt
                .input_identities()
                .iter()
                .map(|identity| {
                    serde_json::json!({
                        "role": identity.role().as_str(),
                        "normalized_path": identity.normalized_path().display().to_string(),
                        "observed_size_bytes": identity.observed_size_bytes(),
                        "observed_sha256": identity.observed_sha256(),
                        "declared_size_bytes": identity.declared_size_bytes(),
                        "declared_sha256": identity.declared_sha256(),
                    })
                })
                .collect::<Vec<_>>(),
            "phases": self
                .receipt
                .phases()
                .iter()
                .map(|record| {
                    serde_json::json!({
                        "phase": record.phase().as_str(),
                        "operation_units_cumulative": record.operation_units_cumulative(),
                        "retained_before_bytes": record.retained_before_bytes(),
                        "transient_bytes": record.transient_bytes(),
                        "retained_after_bytes": record.retained_after_bytes(),
                        "phase_peak_bytes": record.phase_peak_bytes(),
                        "output_bytes_reserved": record.output_bytes_reserved(),
                    })
                })
                .collect::<Vec<_>>(),
            "reference_bases_processed": preflight.reference_bases_processed,
            "fasta_bytes_processed": preflight.fasta_bytes_processed,
            "mask_bytes_processed": preflight.mask_bytes_processed,
            "general_toml_bytes": preflight.general_toml_bytes,
            "scenario_toml_bytes": preflight.scenario_toml_bytes,
            "profile_metadata_bytes": preflight.profile_metadata_bytes,
            "profile_tsv_bytes": preflight.profile_tsv_bytes,
            "tract_count": preflight.tract_count,
            "retained_sequence_ids": preflight.retained_sequence_ids,
            "gap_support_cardinality": preflight.gap_support_cardinality,
            "gap_selector_operation_units": preflight.gap_selector_operation_units,
            "site_bound_per_block": preflight.site_bound_per_block,
            "bin_count_total": preflight.bin_count_total,
            "profile_rows": preflight.profile_rows,
            "normalization_rows": preflight.normalization_rows,
            "tsv_rows_bound": preflight.tsv_rows_bound,
            "json_rows_bound": preflight.json_rows_bound,
            "other_operation_units": preflight.other.operation_units(),
            "other_working_memory_bytes": preflight.other.working_memory_bytes(),
            "other_output_bytes": preflight.other.output_bytes(),
            "operation_units": preflight.operation_units,
            "working_memory_bytes": preflight.working_memory_bytes,
            "output_bytes": preflight.output_bytes,
            "operation_limit": preflight.operation_limit,
            "working_memory_limit_bytes": preflight.working_memory_limit_bytes,
            "output_limit_bytes": preflight.output_limit_bytes,
        })
    }

    /// Comprueba el binding completo antes de abrir cualquier ámbito RNG: la
    /// identidad semántica, las aportaciones externas y —lo que .5 añade— las
    /// identidades físicas de los ocho documentos que la configuración leyó.
    ///
    /// Dos TOML semánticamente iguales que difieran en un comentario o en su
    /// relleno tienen la misma identidad semántica y distinto SHA observado;
    /// solo esta segunda mitad del binding impide cruzarlos.
    pub(crate) fn matches_admission(
        &self,
        config: &EffectiveF5Config,
        other: F5OtherResources,
    ) -> bool {
        self.preflight.matches_admission(config, other) && self.matches_physical_documents(config)
    }

    fn matches_physical_documents(&self, config: &EffectiveF5Config) -> bool {
        let expected = config.admitted_document_identities();
        self.receipt
            .input_identities()
            .get(..expected.len())
            .is_some_and(|sealed| sealed == expected)
    }

    /// Solo la mitad semántica del binding. `F5BlockLengthPlan` la usa para
    /// comprobar objetivo, perfil y longitudes; la mitad física pertenece a la
    /// puerta de layout, que es la que abre RNG y publicación.
    pub(crate) fn matches_configuration(&self, config: &EffectiveF5Config) -> bool {
        self.preflight.matches_configuration(config)
    }

    #[must_use]
    pub fn scenario_id(&self) -> &str {
        &self.preflight.scenario_id
    }
    #[must_use]
    pub fn reference_id(&self) -> &str {
        &self.preflight.reference_id
    }
    #[must_use]
    pub fn reference_sha256(&self) -> &str {
        &self.preflight.reference_sha256
    }
    #[must_use]
    pub fn profile_id(&self) -> &str {
        &self.preflight.profile_id
    }
    #[must_use]
    pub fn profile_metadata_sha256(&self) -> &str {
        &self.preflight.profile_metadata_sha256
    }
    #[must_use]
    pub const fn block_length_minimum(&self) -> u64 {
        self.preflight.block_length_minimum
    }
    #[must_use]
    pub const fn block_length_maximum(&self) -> u64 {
        self.preflight.block_length_maximum
    }
    #[must_use]
    pub fn eligible_tracts(&self) -> &[F5EligibleTract] {
        &self.preflight.eligible_tracts
    }
    #[must_use]
    pub const fn eligible_bases(&self) -> u64 {
        self.preflight.eligible_bases
    }
    #[must_use]
    pub const fn reference_bases_processed(&self) -> u64 {
        self.preflight.reference_bases_processed
    }
    /// Bytes crudos del FASTA recorridos por cada una de las dos pasadas F5.
    #[must_use]
    pub const fn fasta_bytes_processed(&self) -> u64 {
        self.preflight.fasta_bytes_processed
    }
    /// Bytes crudos de la máscara recorridos durante el escáner de preflight.
    #[must_use]
    pub const fn mask_bytes_processed(&self) -> u64 {
        self.preflight.mask_bytes_processed
    }
    /// Bytes físicos del TOML general retenidos durante la admisión F5.
    #[must_use]
    pub const fn general_toml_bytes(&self) -> u64 {
        self.preflight.general_toml_bytes
    }
    /// Bytes físicos del TOML de escenario retenidos durante la admisión F5.
    #[must_use]
    pub const fn scenario_toml_bytes(&self) -> u64 {
        self.preflight.scenario_toml_bytes
    }
    /// Bytes físicos del metadata JSON del perfil.
    #[must_use]
    pub const fn profile_metadata_bytes(&self) -> u64 {
        self.preflight.profile_metadata_bytes
    }
    /// Suma de los cinco TSV físicos del perfil.
    #[must_use]
    pub const fn profile_tsv_bytes(&self) -> u64 {
        self.preflight.profile_tsv_bytes
    }
    /// Reserva conservadora de lectura, retención y parser de los documentos
    /// F5, calculada antes de abrir metadata o TSV.
    #[must_use]
    pub const fn input_parser_peak_bytes(&self) -> u64 {
        self.preflight.input_parser_peak_bytes
    }
    #[must_use]
    pub const fn tract_count(&self) -> u64 {
        self.preflight.tract_count
    }
    /// Pico de IDs FASTA retenidos simultáneamente para detectar duplicados.
    #[must_use]
    pub const fn retained_sequence_ids(&self) -> u64 {
        self.preflight.retained_sequence_ids
    }
    /// Número de gaps enteros con probabilidad positiva, sumado por los
    /// intervalos disjuntos del perfil sin rellenar huecos entre ellos.
    #[must_use]
    pub const fn gap_support_cardinality(&self) -> u64 {
        self.preflight.gap_support_cardinality
    }
    /// Cota conservadora de visitas del selector literal de gaps: dos pasadas
    /// completas por cada sorteo potencial de cada bloque e intento.
    #[must_use]
    pub const fn gap_selector_operation_units(&self) -> u64 {
        self.preflight.gap_selector_operation_units
    }
    #[must_use]
    pub const fn site_bound_per_block(&self) -> u64 {
        self.preflight.site_bound_per_block
    }
    #[must_use]
    pub const fn bin_count_total(&self) -> u64 {
        self.preflight.bin_count_total
    }
    #[must_use]
    pub const fn profile_rows(&self) -> u64 {
        self.preflight.profile_rows
    }
    #[must_use]
    pub const fn normalization_rows(&self) -> u64 {
        self.preflight.normalization_rows
    }
    #[must_use]
    pub const fn tsv_rows_bound(&self) -> u64 {
        self.preflight.tsv_rows_bound
    }
    #[must_use]
    pub const fn json_rows_bound(&self) -> u64 {
        self.preflight.json_rows_bound
    }
    #[must_use]
    pub const fn other_resources(&self) -> F5OtherResources {
        self.preflight.other
    }
    #[must_use]
    pub const fn operation_units(&self) -> u64 {
        self.preflight.operation_units
    }
    #[must_use]
    pub const fn working_memory_bytes(&self) -> u64 {
        self.preflight.working_memory_bytes
    }
    #[must_use]
    pub const fn output_bytes(&self) -> u64 {
        self.preflight.output_bytes
    }
    #[must_use]
    pub const fn operation_limit(&self) -> u64 {
        self.preflight.operation_limit
    }
    #[must_use]
    pub const fn working_memory_limit_bytes(&self) -> u64 {
        self.preflight.working_memory_limit_bytes
    }
    #[must_use]
    pub const fn output_limit_bytes(&self) -> u64 {
        self.preflight.output_limit_bytes
    }
    /// Límite efectivo que el preflight admitió para el consumidor de layout.
    #[must_use]
    pub const fn max_layout_attempts(&self) -> u64 {
        self.preflight.max_layout_attempts
    }
    #[must_use]
    pub const fn block_count(&self) -> u64 {
        self.preflight.block_count
    }
    #[must_use]
    pub const fn target_bases(&self) -> Option<u64> {
        self.preflight.target_bases
    }
    #[must_use]
    pub const fn resolved_target(&self) -> &F5ResolvedTarget {
        &self.preflight.resolved_target
    }
}

/// Recupera la única admisión sellada, calculada antes de RNG, layout o
/// publicación. No reabre entradas ni acepta aportaciones distintas a las
/// ligadas entonces.
/// Recupera la admisión ya sellada de una configuración y comprueba que las
/// aportaciones externas coinciden con las que se admitieron.
///
/// No sella: el sellado ocurre una sola vez, al materializar la configuración.
/// La garantía física de este hito la impone `F5LayoutRequest`, que exige el
/// receipt antes de abrir RNG o publicación.
pub fn recover_sealed_f5_admission(
    config: &EffectiveF5Config,
    other: F5OtherResources,
) -> Result<SealedF5Admission, Diagnostic> {
    let admitted = config.admitted().ok_or_else(|| {
        resource(
            config,
            "la configuración F5 no conserva su admisión física sellada",
        )
    })?;
    if !admitted.matches_admission(config, other) {
        return Err(f5_configuration(
            config,
            "preflight.other_resources",
            "las aportaciones externas deben coincidir con las admitidas antes de abrir referencia y máscara",
        ));
    }
    Ok(admitted.clone())
}

/// Reúne, en orden de rol contractual, las identidades de las diez entradas
/// posibles: los dos TOML de bootstrap, el metadata, los cinco TSV, la
/// referencia y la máscara cuando el escenario la declara.
///
/// Las identidades no se recalculan: cada una procede del mismo handle que
/// contó y hasheó sus bytes durante su fase.
pub(crate) fn sealed_input_identities(
    document_identities: &[AdmittedInputIdentity],
    preflight: &F5ResourcePreflight,
) -> Vec<AdmittedInputIdentity> {
    let mut identities =
        Vec::with_capacity(document_identities.len() + preflight.physical_input_identities.len());
    identities.extend_from_slice(document_identities);
    identities.extend_from_slice(&preflight.physical_input_identities);
    identities
}

/// Ejecuta .4 sobre el ledger procedente del perfil. El resultado se conserva
/// sin sellar para .5 y para llamadas posteriores al preflight sin nueva E/S.
pub(crate) fn admit_f5_resources_under_ledger(
    config: &EffectiveF5Config,
    other: F5OtherResources,
    ledger: &mut F5ResourceLedger,
) -> Result<F5ResourcePreflight, Diagnostic> {
    reserve_known_before_scan(config, other, ledger)?;
    let scanned = scan_reference_and_mask_under_ledger(config, other, ledger)?;
    preflight_resolved_f5_resources(config, scanned, other, Some(ledger))
}

/// Reserva antes de la primera apertura todos los costes futuros que ya son
/// resolubles sin observar la referencia. Los términos dependientes de
/// secuencias/tramos se completan en `FinalFormula` sobre este mismo ledger.
fn reserve_known_before_scan(
    config: &EffectiveF5Config,
    other: F5OtherResources,
    ledger: &mut F5ResourceLedger,
) -> Result<(), Diagnostic> {
    let profile = config.profile();
    let l_min = profile
        .block_lengths()
        .first()
        .ok_or_else(|| domain(config, "perfil sin longitudes"))?
        .minimum();
    let l_max = profile
        .block_lengths()
        .last()
        .ok_or_else(|| domain(config, "perfil sin longitudes"))?
        .maximum();
    let gap_min = profile
        .snp_gaps()
        .first()
        .ok_or_else(|| domain(config, "perfil sin gaps"))?
        .minimum();
    let site_bound = add(
        "S_bound",
        1,
        div(
            "S_bound",
            sub("S_bound", l_max, 1, config)?,
            gap_min,
            config,
        )?,
        config,
    )?;
    if site_bound > 4_096 {
        return Err(resource(config, "site_bound_per_block supera 4096"));
    }
    let h = checked_usize(
        "haplotype_weights",
        profile.haplotype_weights().len(),
        config,
    )?;
    let max_cells = limit(
        config.declared_scenario().max_haplotype_cells.as_ref(),
        DEFAULT_MAX_HAPLOTYPE_CELLS,
    );
    if mul("H*S_bound", h, site_bound, config)? > max_cells {
        return Err(resource(
            config,
            "H*site_bound_per_block supera max_haplotype_cells",
        ));
    }
    let n_pop = checked_usize(
        "subpopulations",
        config.declared_scenario().subpopulations.len(),
        config,
    )?;
    let bin_count_total = sum_bins(config)?;
    let normalization_rows = profile
        .normalizations()
        .iter()
        .try_fold(0_u64, |total, audit| {
            add(
                "normalization_rows",
                total,
                checked_usize("normalization_rows", audit.input_values().len(), config)?,
                config,
            )
        })?;
    let profile_rows = add_many(
        config,
        "profile_rows",
        [
            checked_usize("block_length_rows", profile.block_lengths().len(), config)?,
            checked_usize("snp_gap_rows", profile.snp_gaps().len(), config)?,
            checked_usize("snp_spectrum_rows", profile.snp_spectrum().len(), config)?,
            h,
            bin_count_total,
        ],
    )?;
    let block_count = known_block_count_before_scan(config, l_min, l_max)?.unwrap_or(0);
    let bn = mul("B*N_pop", block_count, n_pop, config)?;
    let founder_count = *config.declared_scenario().founder_count.value();
    let c_total = mul("C_total", founder_count, u64::from(config.ploidy()), config)?;
    let layout_attempts = limit(
        config.declared_scenario().max_layout_attempts.as_ref(),
        F5_DEFAULT_MAX_LAYOUT_ATTEMPTS,
    );
    let haplotype_attempts = limit(
        config.declared_scenario().max_haplotype_attempts.as_ref(),
        DEFAULT_MAX_HAPLOTYPE_ATTEMPTS,
    );
    let swap_proposals = limit(
        config.declared_scenario().max_swap_proposals.as_ref(),
        DEFAULT_MAX_SWAP_PROPOSALS,
    );
    let k_bound = mul("K_bound", block_count, site_bound, config)?;
    let length_support = sub("length_support", l_max, l_min, config)?
        .checked_add(1)
        .ok_or_else(|| resource(config, "overflow en length_support"))?;
    let gap_support = profile
        .snp_gaps()
        .iter()
        .try_fold(0_u64, |total, interval| {
            let width = sub(
                "gap_interval_width",
                interval.maximum(),
                interval.minimum(),
                config,
            )?
            .checked_add(1)
            .ok_or_else(|| resource(config, "overflow en gap_interval_width"))?;
            add("gap_support", total, width, config)
        })?;
    let known_layout_fragment_growth = if block_count == 0 {
        0
    } else {
        triangular("B*(B-1)/2", block_count, config)?
    };
    let pair_bound = div(
        "pair_bound",
        mul(
            "S^2",
            site_bound,
            sub("S-1", site_bound, 1, config)?,
            config,
        )?,
        2,
        config,
    )?;
    let tsv_rows = add_many(
        config,
        "known_tsv_rows",
        [
            block_count,
            k_bound,
            mul("B*N*H", bn, h, config)?,
            mul("B*N*S", bn, site_bound, config)?,
            mul("B*bins", block_count, bin_count_total, config)?,
            mul("B*C", block_count, c_total, config)?,
            mul("B*N*A", bn, haplotype_attempts, config)?,
        ],
    )?;
    let json_rows = add_many(
        config,
        "known_json_rows",
        [
            founder_count,
            n_pop,
            mul("3*B*N", 3, bn, config)?,
            mul("2*B*N*H", 2, mul("B*N*H", bn, h, config)?, config)?,
            mul("B*N*S", bn, site_bound, config)?,
            mul("2*normalization_rows", 2, normalization_rows, config)?,
            9,
        ],
    )?;
    let known_operations = add_many(
        config,
        "known_future_operation_units",
        [
            other.operation_units(),
            profile_rows,
            normalization_rows,
            mul(
                "2*L*B*length_support",
                2,
                mul(
                    "L*B*length_support",
                    mul("L*B", layout_attempts, block_count, config)?,
                    length_support,
                    config,
                )?,
                config,
            )?,
            mul(
                "2*L*B*S*gap_support",
                2,
                mul(
                    "L*B*S*gap_support",
                    mul(
                        "L*B*S",
                        mul("L*B", layout_attempts, block_count, config)?,
                        site_bound,
                        config,
                    )?,
                    gap_support,
                    config,
                )?,
                config,
            )?,
            mul(
                "L*B*(B+S)",
                mul("L*B", layout_attempts, block_count, config)?,
                add("B+S", block_count, site_bound, config)?,
                config,
            )?,
            // Partes de consulta y desplazamiento que solo dependen de B,
            // no del número de tramos observado: 2*L*tri + 3*L*tri.
            mul(
                "5*L*B*(B-1)/2",
                5,
                mul(
                    "L*B*(B-1)/2",
                    layout_attempts,
                    known_layout_fragment_growth,
                    config,
                )?,
                config,
            )?,
            mul(
                "B*N*A*(H*S^2+Q*S)",
                mul("B*N*A", bn, haplotype_attempts, config)?,
                add(
                    "H*S^2+Q*S",
                    mul(
                        "H*S^2",
                        h,
                        mul("S^2", site_bound, site_bound, config)?,
                        config,
                    )?,
                    mul("Q*S", swap_proposals, site_bound, config)?,
                    config,
                )?,
                config,
            )?,
            mul("B*C", block_count, c_total, config)?,
            mul("G", founder_count, k_bound, config)?,
            mul("B*N*S", bn, site_bound, config)?,
            tsv_rows,
            json_rows,
        ],
    )?;
    ledger.charge(
        InputRole::ReferenceFasta,
        LedgerCharge::Operation(OperationCost::new(known_operations)),
    )?;
    // Las sondas EOF no son operaciones publicables: reservan como respaldo
    // una parte de trabajo futuro ya cargado y la devuelven al observar cero
    // bytes. Si aparece un byte real, consume ese respaldo y la fórmula final
    // debe admitir su coste adicional.
    ledger.earmark_read_operation_escrow(
        InputRole::ReferenceFasta,
        2 + u64::from(config.excluded_regions_path().is_some()),
    )?;
    let known_output = add_many(
        config,
        "known_future_output_bytes",
        [
            other.output_bytes(),
            1_048_576,
            mul(
                "4096*(tsv_rows+json_rows)",
                4_096,
                add("tsv_rows+json_rows", tsv_rows, json_rows, config)?,
                config,
            )?,
            65_536,
            mul(
                "4096*G",
                4_096,
                mul("G", founder_count, k_bound, config)?,
                config,
            )?,
        ],
    )?;
    ledger.charge(
        InputRole::ReferenceFasta,
        LedgerCharge::Output(OutputReservation::new(known_output)),
    )?;

    // El pico final conocido no incluye secuencias/tramos todavía, pero sí
    // la retención del perfil y las estructuras derivables del escenario.
    let hs = mul("H*S", h, site_bound, config)?;
    let known_final_transient = add_many(
        config,
        "known_final_transient",
        [
            other.working_memory_bytes(),
            mul("1024*B", 1_024, block_count, config)?,
            mul(
                "B*N*population_memory",
                bn,
                add_many(
                    config,
                    "population_memory",
                    [
                        ceil_div("ceil(H*S/8)", hs, 8, config)?,
                        mul("64*H", 64, h, config)?,
                        mul("256*S", 256, site_bound, config)?,
                        mul("192*pair", 192, pair_bound, config)?,
                    ],
                )?,
                config,
            )?,
            mul(
                "256*B*bins",
                256,
                mul("B*bins", block_count, bin_count_total, config)?,
                config,
            )?,
            mul("128*C", 128, c_total, config)?,
            F5_SCAN_FIXED_MEMORY_BYTES,
        ],
    )?;
    ledger.reserve_transient_memory(InputRole::ReferenceFasta, known_final_transient)
}

/// Dimensión B que ya es autoridad antes de observar FASTA. `coverage_fraction`
/// depende de E y se completa después del scan; `covered_bases`, en cambio,
/// debe validar su representabilidad y reservar con B=ceil(T/Lmax) ahora.
fn known_block_count_before_scan(
    config: &EffectiveF5Config,
    l_min: u64,
    l_max: u64,
) -> Result<Option<u64>, Diagnostic> {
    match &config.declared_scenario().target {
        DeclaredF5Target::BlockCount(value) => Ok(Some(*value.value())),
        DeclaredF5Target::CoveredBases(value) => {
            let target = *value.value();
            if target == 0 {
                return Err(domain(config, "target_bases debe ser positivo"));
            }
            let block_count = ceil_div("B", target, l_max, config)?;
            if mul("B*Lmin", block_count, l_min, config)? > target {
                return Err(domain(
                    config,
                    "target_bases no es representable por el soporte de longitudes",
                ));
            }
            Ok(Some(block_count))
        }
        DeclaredF5Target::CoverageFraction(_) => Ok(None),
    }
}

fn preflight_resolved_f5_resources(
    config: &EffectiveF5Config,
    scanned: F5ReferenceScan,
    other: F5OtherResources,
    ledger: Option<&mut F5ResourceLedger>,
) -> Result<F5ResourcePreflight, Diagnostic> {
    let configuration_identity = F5PreflightConfigurationIdentity::from_config(config)?;
    let F5ReferenceScan {
        eligible_tracts,
        eligible_bases,
        reference_bases_processed,
        fasta_bytes_processed,
        mask_bytes_processed,
        sequence_count,
        retained_sequence_ids,
        physical_input_identities,
    } = scanned;
    let tract_count = checked_usize("tract_count", eligible_tracts.len(), config)?;
    let input_documents = config.input_documents();
    let general_toml_bytes = input_documents.general_toml_bytes();
    let scenario_toml_bytes = input_documents.scenario_toml_bytes();
    let profile_metadata_bytes = input_documents.profile_metadata_bytes();
    let profile_tsv_bytes = input_documents.profile_tsv_bytes();
    let input_parser_peak_bytes = input_documents.parser_peak_bytes(config.scenario_id())?;
    let input_document_bytes = add_many(
        config,
        "input_document_bytes",
        [
            general_toml_bytes,
            scenario_toml_bytes,
            profile_metadata_bytes,
            profile_tsv_bytes,
        ],
    )?
    .checked_mul(2)
    .ok_or_else(|| resource(config, "overflow en 2*input_document_bytes"))?;
    let profile = config.profile();
    let l_min = profile
        .block_lengths()
        .first()
        .ok_or_else(|| domain(config, "perfil sin longitudes"))?
        .minimum();
    let l_max = profile
        .block_lengths()
        .last()
        .ok_or_else(|| domain(config, "perfil sin longitudes"))?
        .maximum();
    let gap_min = profile
        .snp_gaps()
        .first()
        .ok_or_else(|| domain(config, "perfil sin gaps"))?
        .minimum();
    let h = checked_usize(
        "haplotype_weights",
        profile.haplotype_weights().len(),
        config,
    )?;
    let site_bound_per_block = add(
        "S_bound",
        1,
        div(
            "S_bound",
            sub("S_bound", l_max, 1, config)?,
            gap_min,
            config,
        )?,
        config,
    )?;
    if site_bound_per_block > 4_096 {
        return Err(resource(config, "site_bound_per_block supera 4096"));
    }
    let max_cells = limit(
        config.declared_scenario().max_haplotype_cells.as_ref(),
        DEFAULT_MAX_HAPLOTYPE_CELLS,
    );
    if mul("H*S_bound", h, site_bound_per_block, config)? > max_cells {
        return Err(resource(
            config,
            "H*site_bound_per_block supera max_haplotype_cells",
        ));
    }

    let (block_count, target_bases, resolved_target) =
        resolve_capacity(config, eligible_bases, &eligible_tracts, l_min, l_max)?;
    if config
        .declared_scenario()
        .sampling_overrides
        .iter()
        .any(|item| *item.block_index.value() >= block_count)
    {
        return Err(Diagnostic::f5_configuration(
            config.scenario_path().display().to_string(),
            "sampling_overrides.block_index",
            format!("debe pertenecer a 0..{}", block_count - 1),
        ));
    }
    let n_pop = checked_usize(
        "subpopulations",
        config.declared_scenario().subpopulations.len(),
        config,
    )?;
    let founder_count = *config.declared_scenario().founder_count.value();
    let ploidy = u64::from(config.ploidy());
    let c_total = mul("C_total", founder_count, ploidy, config)?;
    let layout_attempts = limit(
        config.declared_scenario().max_layout_attempts.as_ref(),
        F5_DEFAULT_MAX_LAYOUT_ATTEMPTS,
    );
    let haplotype_attempts = limit(
        config.declared_scenario().max_haplotype_attempts.as_ref(),
        DEFAULT_MAX_HAPLOTYPE_ATTEMPTS,
    );
    let swap_proposals = limit(
        config.declared_scenario().max_swap_proposals.as_ref(),
        DEFAULT_MAX_SWAP_PROPOSALS,
    );
    let k_bound = mul("K_bound", block_count, site_bound_per_block, config)?;
    let pair_bound = div(
        "pair_bound",
        mul(
            "pair_bound",
            site_bound_per_block,
            sub("pair_bound", site_bound_per_block, 1, config)?,
            config,
        )?,
        2,
        config,
    )?;
    let x = mul("X", c_total, sequence_count, config)?;
    let g_bound = mul("G_bound", founder_count, k_bound, config)?;
    let bin_count_total = sum_bins(config)?;
    // El selector literal recorre cada entero factible dos veces (total y
    // prefijo). La restricción residual nunca puede ampliar este soporte.
    let length_support_cardinality = sub("length_support_cardinality", l_max, l_min, config)?
        .checked_add(1)
        .ok_or_else(|| resource(config, "overflow en length_support_cardinality"))?;
    // Los gaps pueden tener soporte discontinuo. El selector literal recorre
    // todos los valores soportados para el total Neumaier y, conservadoramente,
    // todos de nuevo para el prefijo. Cada bloque puede pedir como máximo
    // S_bound gaps (S_bound-1 que caben y uno terminal que ya no cabe) en cada
    // intento. Esta cota se calcula sin enumerar el soporte y antes de RNG.
    let gap_support_cardinality =
        profile
            .snp_gaps()
            .iter()
            .try_fold(0_u64, |total, interval| {
                let width = sub(
                    "gap_interval_width",
                    interval.maximum(),
                    interval.minimum(),
                    config,
                )?
                .checked_add(1)
                .ok_or_else(|| resource(config, "overflow en gap_interval_width"))?;
                add("gap_support_cardinality", total, width, config)
            })?;
    let gap_selector_operation_units = mul(
        "2*L*B*S_bound*gap_support_cardinality",
        2,
        mul(
            "L*B*S_bound*gap_support_cardinality",
            mul(
                "L*B*S_bound",
                mul("L*B", layout_attempts, block_count, config)?,
                site_bound_per_block,
                config,
            )?,
            gap_support_cardinality,
            config,
        )?,
        config,
    )?;
    let profile_rows = add_many(
        config,
        "profile_rows",
        [
            checked_usize("block_length_rows", profile.block_lengths().len(), config)?,
            checked_usize("snp_gap_rows", profile.snp_gaps().len(), config)?,
            checked_usize("snp_spectrum_rows", profile.snp_spectrum().len(), config)?,
            h,
            bin_count_total,
        ],
    )?;
    let normalization_rows = profile
        .normalizations()
        .iter()
        .try_fold(0_u64, |total, audit| {
            add(
                "normalization_rows",
                total,
                checked_usize("normalization_rows", audit.input_values().len(), config)?,
                config,
            )
        })?;
    let bn = mul("B*N_pop", block_count, n_pop, config)?;
    let tsv_rows_bound = add_many(
        config,
        "tsv_rows_bound",
        [
            block_count,
            k_bound,
            mul("B*N_pop*H", bn, h, config)?,
            mul("B*N_pop*S_bound", bn, site_bound_per_block, config)?,
            mul("B*bin_count_total", block_count, bin_count_total, config)?,
            mul("B*C_total", block_count, c_total, config)?,
            mul("B*N_pop*A", bn, haplotype_attempts, config)?,
        ],
    )?;
    let json_rows_bound = add_many(
        config,
        "json_rows_bound",
        [
            founder_count,
            n_pop,
            mul("3*B*N_pop", 3, bn, config)?,
            mul("2*B*N_pop*H", 2, mul("B*N_pop*H", bn, h, config)?, config)?,
            mul("B*N_pop*S_bound", bn, site_bound_per_block, config)?,
            mul("2*normalization_rows", 2, normalization_rows, config)?,
            9,
        ],
    )?;

    let layout_fragment_growth = triangular("B*(B-1)/2", block_count, config)?;
    let operation_units = add_many(
        config,
        "operation_units",
        [
            other.operation_units(),
            // Los documentos de autoría y perfil se leen y parsean antes del
            // preflight geométrico: dos unidades por byte, observadas por sus
            // propios handles, no por metadata del sistema de ficheros.
            input_document_bytes,
            // El escáner de admisión y la recuperación REF posterior recorren
            // cada byte físico del FASTA, incluidas cabeceras, descripciones y
            // separadores de línea. El cursor de máscara avanza por base.
            mul("2*fasta_bytes_processed", 2, fasta_bytes_processed, config)?,
            reference_bases_processed,
            // La máscara se lee y se hashea una vez durante el preflight. Sus
            // IDs y filas no están acotados por la geometría de la referencia.
            mask_bytes_processed,
            retained_sequence_ids,
            tract_count,
            // v4: una visita para comprobar la capacidad de los tramos y una
            // para validar la lista enlazada antes de cualquier RNG.
            mul("2*tract_count", 2, tract_count, config)?,
            mul(
                "2*L*(B*tract_count+B*(B-1)/2)",
                2,
                mul(
                    "L*(B*tract_count+B*(B-1)/2)",
                    layout_attempts,
                    add(
                        "B*tract_count+B*(B-1)/2",
                        mul("B*tract_count", block_count, tract_count, config)?,
                        layout_fragment_growth,
                        config,
                    )?,
                    config,
                )?,
                config,
            )?,
            // Cada intento copia la lista inicial y puede desplazar hasta tres
            // veces los tramos vivos por bloque (remove + hasta dos insert).
            // El término triangular cubre que cada división aumenta la lista
            // como máximo en un tramo.
            mul(
                "L*(tract_count+3*B*tract_count+3*B*(B-1)/2)",
                layout_attempts,
                add_many(
                    config,
                    "tract_count+3*B*tract_count+3*B*(B-1)/2",
                    [
                        tract_count,
                        mul(
                            "3*B*tract_count",
                            3,
                            mul("B*tract_count", block_count, tract_count, config)?,
                            config,
                        )?,
                        mul("3*B*(B-1)/2", 3, layout_fragment_growth, config)?,
                    ],
                )?,
                config,
            )?,
            profile_rows,
            normalization_rows,
            mul(
                "2*L*B*length_support_cardinality",
                2,
                mul(
                    "L*B*length_support_cardinality",
                    mul("L*B", layout_attempts, block_count, config)?,
                    length_support_cardinality,
                    config,
                )?,
                config,
            )?,
            gap_selector_operation_units,
            mul(
                "L*B*(B+S_bound)",
                mul("L*B", layout_attempts, block_count, config)?,
                add("B+S_bound", block_count, site_bound_per_block, config)?,
                config,
            )?,
            mul(
                "B*N_pop*A*(H*S^2+Q*S)",
                mul("B*N_pop*A", bn, haplotype_attempts, config)?,
                add(
                    "H*S^2+Q*S",
                    mul(
                        "H*S^2",
                        h,
                        mul("S^2", site_bound_per_block, site_bound_per_block, config)?,
                        config,
                    )?,
                    mul("Q*S", swap_proposals, site_bound_per_block, config)?,
                    config,
                )?,
                config,
            )?,
            mul("B*C_total", block_count, c_total, config)?,
            mul("X*K_bound", x, k_bound, config)?,
            g_bound,
            mul("B*N_pop*S_bound", bn, site_bound_per_block, config)?,
            tsv_rows_bound,
            json_rows_bound,
        ],
    )?;
    let hs = mul("H*S_bound", h, site_bound_per_block, config)?;
    let final_formula_peak = add_many(
        config,
        "final_formula_peak",
        [
            other.working_memory_bytes(),
            mul(
                "256*(4*B+2*tract_count+retained_sequence_ids+profile_rows+normalization_rows)",
                F5_SCAN_RETAINED_ITEM_BYTES,
                add_many(
                    config,
                    "4*B+2*tract_count+retained_sequence_ids+profile_rows+normalization_rows",
                    [
                        mul("4*B", 4, block_count, config)?,
                        mul("2*tract_count", 2, tract_count, config)?,
                        retained_sequence_ids,
                        profile_rows,
                        normalization_rows,
                    ],
                )?,
                config,
            )?,
            mul(
                "B*N_pop*(ceil(H*S/8)+64*H+256*S+192*pair)",
                bn,
                add_many(
                    config,
                    "per_population_memory",
                    [
                        ceil_div("ceil(H*S/8)", hs, 8, config)?,
                        mul("64*H", 64, h, config)?,
                        mul("256*S", 256, site_bound_per_block, config)?,
                        mul("192*pair", 192, pair_bound, config)?,
                    ],
                )?,
                config,
            )?,
            mul(
                "256*B*bin_count_total",
                256,
                mul("B*bin_count_total", block_count, bin_count_total, config)?,
                config,
            )?,
            mul("128*C_total", 128, c_total, config)?,
            F5_SCAN_FIXED_MEMORY_BYTES,
        ],
    )?;
    // La memoria no es una cuota acumulativa: cada fase conserva solo sus
    // estructuras propias. El resultado contractual es el mayor pico físico,
    // no el presupuesto menos un máximo histórico ni una suma de parsers.
    let bootstrap_peak = add_many(
        config,
        "bootstrap_peak",
        [
            other.working_memory_bytes(),
            mul(
                "2*(general_toml_bytes+scenario_toml_bytes)",
                2,
                add(
                    "general_toml_bytes+scenario_toml_bytes",
                    general_toml_bytes,
                    scenario_toml_bytes,
                    config,
                )?,
                config,
            )?,
        ],
    )?;
    let profile_metadata_peak = add_many(
        config,
        "profile_metadata_peak",
        [
            other.working_memory_bytes(),
            mul(
                "2*profile_metadata_bytes",
                2,
                profile_metadata_bytes,
                config,
            )?,
        ],
    )?;
    let profile_table_peak = add_many(
        config,
        "profile_table_peak",
        [
            other.working_memory_bytes(),
            mul(
                "2*profile_tsv_peak_bytes",
                2,
                input_documents.profile_tsv_peak_bytes(),
                config,
            )?,
            mul(
                "256*(profile_rows+normalization_rows)",
                F5_SCAN_RETAINED_ITEM_BYTES,
                add(
                    "profile_rows+normalization_rows",
                    profile_rows,
                    normalization_rows,
                    config,
                )?,
                config,
            )?,
        ],
    )?;
    let reference_scan_peak = add_many(
        config,
        "reference_and_mask_scan_peak",
        [
            other.working_memory_bytes(),
            F5_SCAN_FIXED_MEMORY_BYTES,
            if config.excluded_regions_path().is_some() {
                MAX_MASK_LINE_BYTES as u64
            } else {
                0
            },
            mul(
                "256*(profile_rows+normalization_rows)",
                F5_SCAN_RETAINED_ITEM_BYTES,
                add(
                    "profile_rows+normalization_rows",
                    profile_rows,
                    normalization_rows,
                    config,
                )?,
                config,
            )?,
            mul(
                "256*(tract_count+retained_sequence_ids)",
                F5_SCAN_RETAINED_ITEM_BYTES,
                add(
                    "tract_count+retained_sequence_ids",
                    tract_count,
                    retained_sequence_ids,
                    config,
                )?,
                config,
            )?,
        ],
    )?;
    let working_memory_bytes = [
        bootstrap_peak,
        profile_metadata_peak,
        profile_table_peak,
        reference_scan_peak,
        final_formula_peak,
    ]
    .into_iter()
    .max()
    .expect("five non-empty F5 memory phases");
    let truth_bound = add_many(
        config,
        "truth_bound",
        [
            65_536,
            mul(
                "X*(4096+2048*K_bound)",
                x,
                add(
                    "4096+2048*K_bound",
                    4_096,
                    mul("2048*K_bound", 2_048, k_bound, config)?,
                    config,
                )?,
                config,
            )?,
            mul("4096*G_bound", 4_096, g_bound, config)?,
        ],
    )?;
    let output_bytes = add_many(
        config,
        "output_bytes",
        [
            other.output_bytes(),
            1_048_576,
            mul(
                "4096*(tsv_rows_bound+json_rows_bound)",
                4_096,
                add(
                    "tsv_rows_bound+json_rows_bound",
                    tsv_rows_bound,
                    json_rows_bound,
                    config,
                )?,
                config,
            )?,
            truth_bound,
        ],
    )?;

    let operation_limit = limit(
        config.declared_scenario().max_operation_units.as_ref(),
        F5_HARD_OPERATION_LIMIT,
    );
    let working_memory_limit_bytes = limit(
        config.declared_scenario().max_working_memory_bytes.as_ref(),
        F5_HARD_WORKING_MEMORY_LIMIT,
    );
    let output_limit_bytes = limit(
        config.declared_scenario().max_output_bytes.as_ref(),
        F5_HARD_OUTPUT_LIMIT,
    );
    if operation_units > operation_limit
        || working_memory_bytes > working_memory_limit_bytes
        || output_bytes > output_limit_bytes
    {
        return Err(resource(
            config,
            format!(
                "operation_units={operation_units}/{operation_limit}, working_memory_bytes={working_memory_bytes}/{working_memory_limit_bytes}, output_bytes={output_bytes}/{output_limit_bytes}, length_support_cardinality={length_support_cardinality}, gap_support_cardinality={gap_support_cardinality}, gap_selector_operation_units={gap_selector_operation_units}"
            ),
        ));
    }
    if let Some(ledger) = ledger {
        ledger.transition_to(AdmissionPhase::FinalFormula)?;
        let charged_operations = ledger.totals().operation_units().value();
        let remaining_operations =
            operation_units
                .checked_sub(charged_operations)
                .ok_or_else(|| {
                    resource(
                        config,
                        "el ledger incremental excede la fórmula final de operaciones",
                    )
                })?;
        ledger.charge(
            InputRole::ReferenceFasta,
            LedgerCharge::Operation(OperationCost::new(remaining_operations)),
        )?;
        let charged_output = ledger.totals().output_bytes().value();
        let remaining_output = output_bytes.checked_sub(charged_output).ok_or_else(|| {
            resource(
                config,
                "el ledger incremental excede la fórmula final de salida",
            )
        })?;
        ledger.charge(
            InputRole::ReferenceFasta,
            LedgerCharge::Output(OutputReservation::new(remaining_output)),
        )?;
        ledger.charge(
            InputRole::ReferenceFasta,
            LedgerCharge::WorkingMemoryPeak(WorkingMemoryPeak::new(working_memory_bytes)),
        )?;
    }
    Ok(F5ResourcePreflight {
        configuration_identity,
        scenario_id: config.scenario_id().to_owned(),
        reference_id: config.reference_id().to_owned(),
        reference_sha256: config.reference_sha256().to_owned(),
        profile_id: profile.profile_id().to_owned(),
        profile_metadata_sha256: config
            .declared_scenario()
            .profile_metadata_sha256
            .value()
            .to_owned(),
        block_length_minimum: l_min,
        block_length_maximum: l_max,
        eligible_tracts: Arc::new(eligible_tracts),
        eligible_bases,
        reference_bases_processed,
        fasta_bytes_processed,
        mask_bytes_processed,
        general_toml_bytes,
        scenario_toml_bytes,
        profile_metadata_bytes,
        profile_tsv_bytes,
        input_parser_peak_bytes,
        tract_count,
        retained_sequence_ids,
        gap_support_cardinality,
        gap_selector_operation_units,
        site_bound_per_block,
        bin_count_total,
        profile_rows,
        normalization_rows,
        tsv_rows_bound,
        json_rows_bound,
        other,
        operation_units,
        working_memory_bytes,
        output_bytes,
        operation_limit,
        working_memory_limit_bytes,
        output_limit_bytes,
        max_layout_attempts: layout_attempts,
        block_count,
        target_bases,
        resolved_target,
        physical_input_identities,
    })
}

fn resolve_capacity(
    config: &EffectiveF5Config,
    eligible: u64,
    tracts: &[F5EligibleTract],
    l_min: u64,
    l_max: u64,
) -> Result<(u64, Option<u64>, F5ResolvedTarget), Diagnostic> {
    match &config.declared_scenario().target {
        DeclaredF5Target::BlockCount(value) => {
            let b = *value.value();
            if mul("B*Lmin", b, l_min, config)? > eligible
                || tracts.iter().try_fold(0_u64, |sum, tract| {
                    record_capacity_tract_visit();
                    add("sum(floor(tract/Lmin))", sum, tract.len() / l_min, config)
                })? < b
            {
                return Err(domain(
                    config,
                    "block_count no cabe en los tramos elegibles",
                ));
            }
            Ok((
                b,
                None,
                F5ResolvedTarget::BlockCount {
                    declared_block_count: b,
                },
            ))
        }
        DeclaredF5Target::CoveredBases(value) => {
            let target = *value.value();
            target_capacity(config, target, eligible, tracts, l_min, l_max).map(|(b, t)| {
                (
                    b,
                    t,
                    F5ResolvedTarget::CoveredBases {
                        declared_covered_bases: target,
                    },
                )
            })
        }
        DeclaredF5Target::CoverageFraction(value) => {
            let target = fractional_target(eligible, value.literal(), config)?;
            let (numerator, denominator) = coverage_fraction_rational(value.literal(), config)?;
            target_capacity(config, target, eligible, tracts, l_min, l_max).map(|(b, t)| {
                (
                    b,
                    t,
                    F5ResolvedTarget::CoverageFraction {
                        numerator,
                        denominator,
                    },
                )
            })
        }
    }
}

/// Reconstruye el enlace del objetivo desde una configuración ya admitida.
/// El consumidor de F5.1.1 lo usa para rechazar un preflight cruzado antes de
/// abrir RNG o aceptar una geometría que correspondía a otro objetivo.
pub(crate) fn resolved_target_for_config(
    config: &EffectiveF5Config,
) -> Result<F5ResolvedTarget, Diagnostic> {
    match &config.declared_scenario().target {
        DeclaredF5Target::BlockCount(value) => Ok(F5ResolvedTarget::BlockCount {
            declared_block_count: *value.value(),
        }),
        DeclaredF5Target::CoveredBases(value) => Ok(F5ResolvedTarget::CoveredBases {
            declared_covered_bases: *value.value(),
        }),
        DeclaredF5Target::CoverageFraction(value) => {
            let (numerator, denominator) = coverage_fraction_rational(value.literal(), config)?;
            Ok(F5ResolvedTarget::CoverageFraction {
                numerator,
                denominator,
            })
        }
    }
}

fn target_capacity(
    config: &EffectiveF5Config,
    target: u64,
    eligible: u64,
    tracts: &[F5EligibleTract],
    l_min: u64,
    l_max: u64,
) -> Result<(u64, Option<u64>), Diagnostic> {
    if target == 0
        || target > eligible
        || tracts.iter().try_fold(0_u64, |sum, tract| {
            record_capacity_tract_visit();
            add("sum(tract_length)", sum, tract.len(), config)
        })? < target
    {
        return Err(domain(
            config,
            "target_bases no cabe en la referencia elegible",
        ));
    }
    let b = ceil_div("B", target, l_max, config)?;
    if mul("B*Lmin", b, l_min, config)? > target {
        return Err(domain(
            config,
            "target_bases no es representable por el soporte de longitudes",
        ));
    }
    Ok((b, Some(target)))
}

fn fractional_target(
    eligible: u64,
    literal: &str,
    config: &EffectiveF5Config,
) -> Result<u64, Diagnostic> {
    ExactDecimal::parse("f5.coverage_fraction", literal)
        .and_then(|fraction| fraction.floor_product(eligible))
        .map_err(|_| resource(config, "coverage_fraction exacta no cabe en u64"))
}

fn coverage_fraction_rational(
    literal: &str,
    config: &EffectiveF5Config,
) -> Result<(u64, u64), Diagnostic> {
    let (whole, fraction) = literal.split_once('.').ok_or_else(|| {
        f5_configuration(config, "f5.coverage_fraction", "literal decimal ausente")
    })?;
    let places = u32::try_from(fraction.len()).map_err(|_| {
        f5_configuration(
            config,
            "f5.coverage_fraction",
            "la escala decimal no cabe en u32",
        )
    })?;
    let denominator = 10_u64.checked_pow(places).ok_or_else(|| {
        f5_configuration(
            config,
            "f5.coverage_fraction",
            "el denominador decimal no cabe en u64",
        )
    })?;
    let whole = whole.parse::<u64>().map_err(|_| {
        f5_configuration(
            config,
            "f5.coverage_fraction",
            "la parte entera no cabe en u64",
        )
    })?;
    let fraction = fraction.parse::<u64>().map_err(|_| {
        f5_configuration(
            config,
            "f5.coverage_fraction",
            "la parte fraccionaria no cabe en u64",
        )
    })?;
    let numerator = whole
        .checked_mul(denominator)
        .and_then(|value| value.checked_add(fraction))
        .ok_or_else(|| {
            f5_configuration(
                config,
                "f5.coverage_fraction",
                "el numerador decimal no cabe en u64",
            )
        })?;
    let divisor = gcd(numerator, denominator);
    Ok((numerator / divisor, denominator / divisor))
}

const fn gcd(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

fn sum_bins(config: &EffectiveF5Config) -> Result<u64, Diagnostic> {
    config
        .declared_scenario()
        .subpopulations
        .iter()
        .try_fold(0_u64, |sum, population| {
            let bins = config
                .profile()
                .ld_bins(population.population_id.value())
                .ok_or_else(|| domain(config, "perfil LD sin población"))?;
            add(
                "bin_count_total",
                sum,
                checked_usize("bin_count", bins.len(), config)?,
                config,
            )
        })
}

#[derive(Clone, Debug)]
pub(crate) struct F5ReferenceScan {
    eligible_tracts: Vec<F5EligibleTract>,
    eligible_bases: u64,
    reference_bases_processed: u64,
    fasta_bytes_processed: u64,
    mask_bytes_processed: u64,
    sequence_count: u64,
    retained_sequence_ids: u64,
    physical_input_identities: Vec<AdmittedInputIdentity>,
}

struct EligibleScan<'a> {
    config: &'a EffectiveF5Config,
    masks: MaskStream<'a>,
    tracts: Vec<F5EligibleTract>,
    retained_id_count: usize,
    eligible: u64,
    processed: u64,
    ordinal: u64,
    sequence_id: String,
    sequence_length: u64,
    run_start: Option<u64>,
    ledger: &'a mut F5ResourceLedger,
    coexisting_transient_bytes: u64,
    tract_operation_units: u64,
    sequence_operation_units: u64,
}

impl<'a> EligibleScan<'a> {
    fn new(
        config: &'a EffectiveF5Config,
        other: F5OtherResources,
        ledger: &'a mut F5ResourceLedger,
    ) -> Result<Self, Diagnostic> {
        let mask_transient = if config.excluded_regions_path().is_some() {
            MAX_MASK_LINE_BYTES as u64
        } else {
            0
        };
        let coexisting_transient_bytes = F5_SCAN_FIXED_MEMORY_BYTES
            .checked_add(mask_transient)
            .and_then(|value| value.checked_add(other.working_memory_bytes()))
            .ok_or_else(|| resource(config, "overflow de transitorios simultáneos del scan"))?;
        let profile = config.profile();
        let l_min = profile
            .block_lengths()
            .first()
            .ok_or_else(|| domain(config, "perfil sin longitudes"))?
            .minimum();
        let l_max = profile
            .block_lengths()
            .last()
            .ok_or_else(|| domain(config, "perfil sin longitudes"))?
            .maximum();
        let gap_min = profile
            .snp_gaps()
            .first()
            .ok_or_else(|| domain(config, "perfil sin gaps"))?
            .minimum();
        let site_bound = add(
            "S_bound",
            1,
            div(
                "S_bound",
                sub("S_bound", l_max, 1, config)?,
                gap_min,
                config,
            )?,
            config,
        )?;
        let known_block_count = known_block_count_before_scan(config, l_min, l_max)?;
        let layout_attempts = limit(
            config.declared_scenario().max_layout_attempts.as_ref(),
            F5_DEFAULT_MAX_LAYOUT_ATTEMPTS,
        );
        let tract_operation_units = if let Some(block_count) = known_block_count {
            add(
                "incremental_tract_units",
                3,
                mul(
                    "L*(5*B+1)",
                    layout_attempts,
                    add("5*B+1", mul("5*B", 5, block_count, config)?, 1, config)?,
                    config,
                )?,
                config,
            )?
        } else {
            3
        };
        let sequence_operation_units = if let Some(block_count) = known_block_count {
            let founder_copies = mul(
                "C_total",
                *config.declared_scenario().founder_count.value(),
                u64::from(config.ploidy()),
                config,
            )?;
            add(
                "incremental_sequence_units",
                1,
                mul(
                    "C_total*B*S_bound",
                    founder_copies,
                    mul("B*S_bound", block_count, site_bound, config)?,
                    config,
                )?,
                config,
            )?
        } else {
            1
        };
        let scan = Self {
            config,
            masks: MaskStream::new(config, ledger)?,
            tracts: Vec::new(),
            retained_id_count: 0,
            eligible: 0,
            processed: 0,
            ordinal: 0,
            sequence_id: String::new(),
            sequence_length: 0,
            run_start: None,
            ledger,
            coexisting_transient_bytes,
            tract_operation_units,
            sequence_operation_units,
        };
        scan.ensure_memory_bound(0, 0)?;
        Ok(scan)
    }

    fn begin_sequence(&mut self, ordinal: u64, sequence_id: &str) -> Result<(), Diagnostic> {
        // La copia del ID activo y la posterior entrada del mapa de unicidad
        // se admiten antes de retener texto de cabecera.
        self.ordinal = ordinal;
        self.sequence_id.clear();
        self.sequence_id.push_str(sequence_id);
        self.sequence_length = 0;
        self.run_start = None;
        self.masks.begin_sequence(ordinal, sequence_id, self.ledger)
    }

    fn accept_base(&mut self, base: u8) -> Result<(), Diagnostic> {
        self.sequence_length = add(
            "reference_sequence_length",
            self.sequence_length,
            1,
            self.config,
        )?;
        self.processed = add("reference_bases_processed", self.processed, 1, self.config)?;
        let valid = matches!(base, b'A' | b'C' | b'G' | b'T')
            && !self
                .masks
                .excludes(&self.sequence_id, self.sequence_length, self.ledger)?;
        if valid {
            self.eligible = add("eligible_bases", self.eligible, 1, self.config)?;
            if self.run_start.is_none() {
                self.run_start = Some(self.sequence_length);
            }
        } else if let Some(start) = self.run_start.take() {
            self.push_tract(start, self.sequence_length - 1)?;
        }
        Ok(())
    }

    fn finish_sequence(&mut self) -> Result<(), Diagnostic> {
        if let Some(start) = self.run_start.take() {
            self.push_tract(start, self.sequence_length)?;
        }
        self.masks.finish_sequence(
            self.ordinal,
            &self.sequence_id,
            self.sequence_length,
            self.ledger,
        )
    }

    fn push_tract(&mut self, start: u64, end: u64) -> Result<(), Diagnostic> {
        // Tanto el trabajo aprendido como la memoria se cargan antes de que
        // el backing Vec pueda crecer y retener el tramo.
        self.ledger.charge(
            InputRole::ReferenceFasta,
            LedgerCharge::Operation(OperationCost::new(self.tract_operation_units)),
        )?;
        let prospective_tract_count = self
            .tracts
            .len()
            .checked_add(1)
            .ok_or_else(|| resource(self.config, "tract_count no cabe en usize"))?;
        self.ensure_memory_bound(prospective_tract_count, self.retained_id_count)?;
        self.ledger.retain_memory(
            InputRole::ReferenceFasta,
            F5_SCAN_RETAINED_ITEM_BYTES,
            self.coexisting_transient_bytes,
        )?;
        self.tracts.try_reserve_exact(1).map_err(|_| {
            resource(
                self.config,
                "la lista de tramos elegibles no cabe en memoria",
            )
        })?;
        self.tracts.push(F5EligibleTract {
            sequence_ordinal: self.ordinal,
            start_1based: start,
            end_1based: end,
        });
        Ok(())
    }

    fn ensure_retained_id_capacity(&self, prospective_id_count: usize) -> Result<(), Diagnostic> {
        self.ensure_memory_bound(self.tracts.len(), prospective_id_count)
    }

    fn record_retained_id(&mut self) -> Result<(), Diagnostic> {
        // Incluye el ID retenido y, cuando B ya es conocido, el trabajo futuro
        // X*K que esta nueva secuencia hace obligatorio.
        self.ledger.charge(
            InputRole::ReferenceFasta,
            LedgerCharge::Operation(OperationCost::new(self.sequence_operation_units)),
        )?;
        self.retained_id_count = self
            .retained_id_count
            .checked_add(1)
            .ok_or_else(|| resource(self.config, "sequence_count no cabe en usize"))?;
        self.ledger.retain_memory(
            InputRole::ReferenceFasta,
            F5_SCAN_RETAINED_ITEM_BYTES,
            self.coexisting_transient_bytes,
        )?;
        Ok(())
    }

    fn ensure_memory_bound(
        &self,
        prospective_tract_count: usize,
        prospective_id_count: usize,
    ) -> Result<(), Diagnostic> {
        let tracts = checked_usize("scan_tract_count", prospective_tract_count, self.config)?;
        let ids = checked_usize("scan_sequence_ids", prospective_id_count, self.config)?;
        let retained_items = add("scan_retained_items", tracts, ids, self.config)?;
        let working_memory_bytes = add_many(
            self.config,
            "scan_working_memory_bytes",
            [
                self.coexisting_transient_bytes,
                mul(
                    "256*scan_retained_items",
                    F5_SCAN_RETAINED_ITEM_BYTES,
                    retained_items,
                    self.config,
                )?,
            ],
        )?;
        let working_memory_limit_bytes = limit(
            self.config
                .declared_scenario()
                .max_working_memory_bytes
                .as_ref(),
            F5_HARD_WORKING_MEMORY_LIMIT,
        );
        if working_memory_bytes > working_memory_limit_bytes {
            // Deuda de F5.1.2.4 revisada al sellar: el ledger rechaza antes por
            // el mismo valor más la retención de perfil, así que esta rama es
            // hoy inalcanzable. Se conserva porque adelanta la comprobación a
            // antes de `SequenceId::parse` —la «comprobación previa a retener»
            // que exige el hito— y se completa su contexto para que, si alguna
            // vez se alcanza, diagnostique como el resto del ledger.
            return Err(resource(
                self.config,
                format!(
                    "working_memory_bytes={working_memory_bytes}/{working_memory_limit_bytes}, phase=reference_scan, tract_count={tracts}, retained_sequence_ids={ids}"
                ),
            )
            .with("algorithm", F5_RESOURCE_ACCOUNTING_ALGORITHM)
            .with("phase", AdmissionPhase::ReferenceAndMaskScan.as_str())
            .with("role", InputRole::ReferenceFasta.as_str())
            .with("dimension", "working_memory_bytes")
            .with("accumulated", working_memory_bytes.to_string())
            .with("limit", working_memory_limit_bytes.to_string()));
        }
        Ok(())
    }

    fn finish(mut self, sequence_count: u64) -> Result<F5ReferenceScan, Diagnostic> {
        self.masks.ensure_exhausted(self.ledger)?;
        let mask_bytes_processed = self.masks.bytes_processed();
        let physical_input_identities = self.masks.identity.into_iter().collect();
        let retained_sequence_ids =
            checked_usize("retained_sequence_ids", self.retained_id_count, self.config)?;
        Ok(F5ReferenceScan {
            eligible_tracts: self.tracts,
            eligible_bases: self.eligible,
            reference_bases_processed: self.processed,
            fasta_bytes_processed: 0,
            mask_bytes_processed,
            sequence_count,
            retained_sequence_ids,
            physical_input_identities,
        })
    }
}

enum FastaLine {
    Start,
    Header(FastaHeaderToken),
    Sequence { pending: Option<u8> },
}

#[derive(Default)]
struct FastaHeaderToken {
    id: Vec<u8>,
    in_description: bool,
}

impl FastaHeaderToken {
    fn accept(&mut self, byte: u8, config: &EffectiveF5Config) -> Result<(), Diagnostic> {
        if self.in_description {
            return Ok(());
        }
        if byte.is_ascii_whitespace() {
            if !self.id.is_empty() {
                self.in_description = true;
            }
            return Ok(());
        }
        if self.id.len() == 64 {
            let mut invalid = String::from_utf8_lossy(&self.id).into_owned();
            invalid.push(char::from(byte));
            return SequenceId::parse(invalid, config.reference_id()).map(|_| ());
        }
        self.id.push(byte);
        Ok(())
    }

    fn finish(self, config: &EffectiveF5Config, line_number: u64) -> Result<String, Diagnostic> {
        if self.id.is_empty() {
            return Err(empty_header(config.reference_path(), line_number));
        }
        String::from_utf8(self.id).map_err(|_| empty_header(config.reference_path(), line_number))
    }
}

struct CurrentFastaSequence {
    id: String,
    header_line: u64,
    length: u64,
}

fn scan_reference_and_mask_under_ledger(
    config: &EffectiveF5Config,
    other: F5OtherResources,
    ledger: &mut F5ResourceLedger,
) -> Result<F5ReferenceScan, Diagnostic> {
    let path = config.reference_path();
    // Reserva antes de abrir: el búfer de streaming y, si procede, una fila
    // TSV acotada coexisten con los datos tipados del perfil ya retenidos.
    let scan_transient = F5_SCAN_FIXED_MEMORY_BYTES
        .checked_add(if config.excluded_regions_path().is_some() {
            MAX_MASK_LINE_BYTES as u64
        } else {
            0
        })
        .and_then(|value| value.checked_add(other.working_memory_bytes()))
        .ok_or_else(|| resource(config, "overflow de transitorios simultáneos del scan"))?;
    reserve_scan_transient(config, ledger, InputRole::ReferenceFasta, scan_transient)?;
    let file = File::open(path).map_err(|error| input_error(path, error))?;
    let mut reader = HashingReader::new(file);

    let mut scan = EligibleScan::new(config, other, ledger)?;
    let mut locations = Vec::<(String, u64)>::new();
    let mut current = None::<CurrentFastaSequence>;
    let mut sequence_count = 0_u64;
    let mut line_number = 1_u64;
    // Un byte pendiente permite reconocer CRLF incluso entre dos lecturas sin
    // conservar una línea de bases completa. De la cabecera solo se retiene
    // el token ID acotado; la descripción se descarta mientras se recorre.
    let mut line = FastaLine::Start;
    let mut buffer = [0_u8; 65_536];

    loop {
        let read = charged_read(
            &mut reader,
            &mut buffer,
            scan.ledger,
            InputRole::ReferenceFasta,
            2,
            path,
        )?;
        if read == 0 {
            break;
        }
        if reader.bytes_processed() <= u64::try_from(read).unwrap_or(u64::MAX)
            && buffer[..read].starts_with(&[0x1f, 0x8b])
        {
            return Err(Diagnostic::format(
                "GPS-FMT-001",
                "El FASTA suministrado está comprimido y v1 solo admite texto sin comprimir.",
            )
            .with("path", path.display().to_string())
            .with("detected_format", "gzip"));
        }
        for byte in &buffer[..read] {
            match &mut line {
                FastaLine::Start if *byte == b'>' => {
                    finish_fasta_sequence(
                        config,
                        &mut scan,
                        &mut current,
                        &mut locations,
                        &mut sequence_count,
                    )?;
                    line = FastaLine::Header(FastaHeaderToken::default());
                }
                FastaLine::Start => {
                    if current.is_none() {
                        return Err(empty_header(path, line_number));
                    }
                    line = FastaLine::Sequence { pending: None };
                    accept_sequence_byte(
                        config,
                        &mut scan,
                        current
                            .as_mut()
                            .ok_or_else(|| empty_header(path, line_number))?,
                        *byte,
                        &mut line,
                    )?;
                    if *byte == b'\n' {
                        line_number = add("FASTA line number", line_number, 1, config)?;
                    }
                }
                FastaLine::Header(header) if *byte == b'\n' => {
                    let header = std::mem::take(header);
                    begin_fasta_sequence(
                        config,
                        &mut scan,
                        &mut current,
                        header,
                        line_number,
                        sequence_count,
                    )?;
                    line = FastaLine::Start;
                    line_number = add("FASTA line number", line_number, 1, config)?;
                }
                FastaLine::Header(header) => header.accept(*byte, config)?,
                FastaLine::Sequence { .. } => {
                    accept_sequence_byte(
                        config,
                        &mut scan,
                        current
                            .as_mut()
                            .ok_or_else(|| empty_header(path, line_number))?,
                        *byte,
                        &mut line,
                    )?;
                    if *byte == b'\n' {
                        line_number = add("FASTA line number", line_number, 1, config)?;
                    }
                }
            }
        }
    }

    match &mut line {
        FastaLine::Start => {}
        FastaLine::Header(header) => {
            let header = std::mem::take(header);
            begin_fasta_sequence(
                config,
                &mut scan,
                &mut current,
                header,
                line_number,
                sequence_count,
            )?;
        }
        FastaLine::Sequence { pending } => {
            if let Some(byte) = pending.take()
                && byte != b'\r'
            {
                accept_validated_base(config, &mut scan, current.as_mut(), byte, path)?;
            }
        }
    }
    finish_fasta_sequence(
        config,
        &mut scan,
        &mut current,
        &mut locations,
        &mut sequence_count,
    )?;
    if sequence_count == 0 {
        return Err(empty_fasta(path));
    }
    ReferenceId::parse(config.reference_id(), "execution")?;
    let (observed_digest, fasta_bytes_processed) = reader.finish();
    let observed_digest = format!("{observed_digest:x}");
    if observed_digest != config.reference_sha256() {
        return Err(Diagnostic::input(
            "GPS-INPUT-002",
            "Una entrada cambió durante su copia o no coincide con su identidad declarada.",
        )
        .with("path", config.reference_path().display().to_string())
        .with("expected_digest", config.reference_sha256())
        .with("observed_digest", &observed_digest));
    }
    let mut result = scan.finish(sequence_count)?;
    result.fasta_bytes_processed = fasta_bytes_processed;
    // El general declara el SHA de la referencia, nunca su tamaño.
    result.physical_input_identities.insert(
        0,
        AdmittedInputIdentity::new(
            PhysicalInputIdentity::from_observed(
                InputRole::ReferenceFasta,
                path.to_path_buf(),
                fasta_bytes_processed,
                observed_digest,
            ),
            None,
            Some(config.reference_sha256().to_owned()),
        )?,
    );
    Ok(result)
}

fn charged_read<R: Read>(
    reader: &mut R,
    buffer: &mut [u8],
    ledger: &mut F5ResourceLedger,
    role: InputRole,
    units_per_byte: u64,
    path: &Path,
) -> Result<usize, Diagnostic> {
    let affordable = ledger.available_read_operation_units() / units_per_byte;
    let affordable = usize::try_from(affordable).unwrap_or(usize::MAX);
    let prospective = buffer.len().min(affordable);
    if prospective == 0 {
        // Produce el diagnóstico normal del ledger sin tocar el handle.
        let _ = ledger.reserve_read(role, 1, units_per_byte)?;
        unreachable!("a zero-capacity reservation cannot succeed")
    }
    let reservation = ledger.reserve_read(role, prospective, units_per_byte)?;
    match reader.read(&mut buffer[..prospective]) {
        Ok(count) => {
            ledger.settle_read(reservation, count)?;
            Ok(count)
        }
        Err(error) => {
            ledger.settle_read(reservation, 0)?;
            Err(input_error(path, error))
        }
    }
}

fn reserve_scan_transient(
    config: &EffectiveF5Config,
    ledger: &mut F5ResourceLedger,
    role: InputRole,
    transient_bytes: u64,
) -> Result<(), Diagnostic> {
    let prospective = ledger
        .retained_memory_bytes()
        .checked_add(transient_bytes)
        .ok_or_else(|| resource(config, "working_memory_bytes desborda durante el escáner"))?;
    let limit = limit(
        config.declared_scenario().max_working_memory_bytes.as_ref(),
        F5_HARD_WORKING_MEMORY_LIMIT,
    );
    if prospective > limit {
        return Err(resource(
            config,
            format!("working_memory_bytes={prospective}/{limit}, phase=reference_and_mask_scan"),
        ));
    }
    ledger.reserve_transient_memory(role, transient_bytes)
}

fn begin_fasta_sequence(
    config: &EffectiveF5Config,
    scan: &mut EligibleScan<'_>,
    current: &mut Option<CurrentFastaSequence>,
    header: FastaHeaderToken,
    line_number: u64,
    sequence_ordinal: u64,
) -> Result<(), Diagnostic> {
    let id = header.finish(config, line_number)?;
    scan.begin_sequence(sequence_ordinal, &id)?;
    *current = Some(CurrentFastaSequence {
        id,
        header_line: line_number,
        length: 0,
    });
    Ok(())
}

fn accept_sequence_byte(
    config: &EffectiveF5Config,
    scan: &mut EligibleScan<'_>,
    current: &mut CurrentFastaSequence,
    byte: u8,
    line: &mut FastaLine,
) -> Result<(), Diagnostic> {
    let FastaLine::Sequence { pending } = line else {
        return Err(resource(config, "estado interno FASTA inválido"));
    };
    if byte == b'\n' {
        if let Some(previous) = pending.take()
            && previous != b'\r'
        {
            accept_validated_base(
                config,
                scan,
                Some(current),
                previous,
                config.reference_path(),
            )?;
        }
        *line = FastaLine::Start;
    } else if let Some(previous) = pending.replace(byte) {
        accept_validated_base(
            config,
            scan,
            Some(current),
            previous,
            config.reference_path(),
        )?;
    }
    Ok(())
}

fn accept_validated_base(
    config: &EffectiveF5Config,
    scan: &mut EligibleScan<'_>,
    current: Option<&mut CurrentFastaSequence>,
    raw: u8,
    path: &std::path::Path,
) -> Result<(), Diagnostic> {
    let current = current.ok_or_else(|| empty_header(path, 0))?;
    scan.ledger.charge(
        InputRole::ReferenceFasta,
        LedgerCharge::Operation(OperationCost::new(1)),
    )?;
    let base = raw.to_ascii_uppercase();
    let position = add("FASTA sequence position", current.length, 1, config)?;
    if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
        return Err(invalid_symbol(path, &current.id, position, raw));
    }
    current.length = position;
    scan.accept_base(base)
}

fn finish_fasta_sequence(
    config: &EffectiveF5Config,
    scan: &mut EligibleScan<'_>,
    current: &mut Option<CurrentFastaSequence>,
    locations: &mut Vec<(String, u64)>,
    sequence_count: &mut u64,
) -> Result<(), Diagnostic> {
    let Some(sequence) = current.take() else {
        return Ok(());
    };
    if sequence.length == 0 {
        return Err(empty_fasta(config.reference_path()).with("sequence", sequence.id));
    }
    if let Some((_, first_line)) = locations.iter().find(|(id, _)| id == &sequence.id) {
        return Err(Diagnostic::format(
            "GPS-FMT-004",
            "Un identificador de secuencia FASTA está repetido.",
        )
        .with("path", config.reference_path().display().to_string())
        .with("identifier", &sequence.id)
        .with(
            "locations",
            format!("line {first_line}, line {}", sequence.header_line),
        ));
    }
    let prospective_id_count = locations
        .len()
        .checked_add(1)
        .ok_or_else(|| resource(config, "sequence_count no cabe en usize"))?;
    scan.ensure_retained_id_capacity(prospective_id_count)?;
    SequenceId::parse(&sequence.id, config.reference_id())?;
    scan.record_retained_id()?;
    locations.try_reserve_exact(1).map_err(|_| {
        resource(
            config,
            "el catálogo de IDs FASTA no cabe en memoria admitida",
        )
    })?;
    locations.push((sequence.id, sequence.header_line));
    scan.finish_sequence()?;
    *sequence_count = add("sequence_count", *sequence_count, 1, config)?;
    Ok(())
}

#[derive(Debug)]
struct MaskInterval {
    sequence_id: String,
    start: u64,
    end: u64,
}

#[derive(Debug)]
struct PreviousMaskInterval {
    sequence_id: String,
    start: u64,
    end: u64,
}

struct MaskStream<'a> {
    reader: Option<HashingReader<File>>,
    pending: Option<MaskInterval>,
    previous: Option<PreviousMaskInterval>,
    next_sequence_ordinal: u64,
    line_buffer: Vec<u8>,
    bytes_processed: u64,
    identity: Option<AdmittedInputIdentity>,
    config: &'a EffectiveF5Config,
}

impl<'a> MaskStream<'a> {
    fn new(
        config: &'a EffectiveF5Config,
        ledger: &mut F5ResourceLedger,
    ) -> Result<Self, Diagnostic> {
        let Some(path) = config.excluded_regions_path() else {
            return Ok(Self {
                reader: None,
                pending: None,
                previous: None,
                next_sequence_ordinal: 0,
                line_buffer: Vec::new(),
                bytes_processed: 0,
                identity: None,
                config,
            });
        };
        let mut header = Vec::new();
        header
            .try_reserve_exact(MAX_MASK_LINE_BYTES)
            .map_err(|_| mask_error(config, "no se puede reservar la fila acotada de máscara"))?;
        let file = File::open(path).map_err(|error| {
            Diagnostic::input(
                "GPS-INPUT-001",
                "No se puede abrir o leer una entrada declarada.",
            )
            .with("path", path.display().to_string())
            .with("detail", error.to_string())
        })?;
        let mut reader = HashingReader::new(file);
        read_bounded_mask_line(&mut reader, &mut header, ledger, config)?;
        if header != b"region_id\tsequence_id\tstart_1based\tend_1based\n" {
            return Err(mask_error(config, "cabecera de máscara inválida"));
        }
        Ok(Self {
            reader: Some(reader),
            pending: None,
            previous: None,
            next_sequence_ordinal: 0,
            line_buffer: header,
            bytes_processed: 0,
            identity: None,
            config,
        })
    }
    fn begin_sequence(
        &mut self,
        ordinal: u64,
        _sequence_id: &str,
        ledger: &mut F5ResourceLedger,
    ) -> Result<(), Diagnostic> {
        self.fill(ledger)?;
        if ordinal != self.next_sequence_ordinal {
            return Err(mask_error(self.config, "ordinal FASTA no contiguo"));
        }
        Ok(())
    }
    fn excludes(
        &mut self,
        sequence_id: &str,
        position: u64,
        ledger: &mut F5ResourceLedger,
    ) -> Result<bool, Diagnostic> {
        self.fill(ledger)?;
        while self
            .pending
            .as_ref()
            .is_some_and(|item| item.sequence_id == sequence_id && item.end < position)
        {
            self.pending = None;
            self.fill(ledger)?;
        }
        Ok(self.pending.as_ref().is_some_and(|item| {
            item.sequence_id == sequence_id && item.start <= position && position <= item.end
        }))
    }
    fn ensure_exhausted(&mut self, ledger: &mut F5ResourceLedger) -> Result<(), Diagnostic> {
        self.fill(ledger)?;
        if self.pending.is_some() {
            Err(mask_error(self.config, "sequence_id desconocido"))
        } else {
            self.verify_digest()
        }
    }
    fn finish_sequence(
        &mut self,
        _sequence_ordinal: u64,
        sequence_id: &str,
        length: u64,
        ledger: &mut F5ResourceLedger,
    ) -> Result<(), Diagnostic> {
        self.fill(ledger)?;
        while self
            .pending
            .as_ref()
            .is_some_and(|item| item.sequence_id == sequence_id)
        {
            if self.pending.as_ref().is_some_and(|item| item.end > length) {
                return Err(mask_error(
                    self.config,
                    "coordenadas de máscara fuera de rango",
                ));
            }
            self.pending = None;
            self.fill(ledger)?;
        }
        self.next_sequence_ordinal = self
            .next_sequence_ordinal
            .checked_add(1)
            .ok_or_else(|| mask_error(self.config, "ordinal FASTA fuera de u64"))?;
        Ok(())
    }
    fn fill(&mut self, ledger: &mut F5ResourceLedger) -> Result<(), Diagnostic> {
        if self.pending.is_some() || self.reader.is_none() {
            return Ok(());
        }
        self.line_buffer.clear();
        let reader = self
            .reader
            .as_mut()
            .ok_or_else(|| mask_error(self.config, "lector de máscara ausente"))?;
        let count = read_bounded_mask_line(reader, &mut self.line_buffer, ledger, self.config)?;
        if count == 0 {
            return Ok(());
        }
        if !self.line_buffer.ends_with(b"\n") {
            return Err(mask_error(
                self.config,
                "la máscara debe terminar cada fila con LF",
            ));
        }
        let row = self
            .line_buffer
            .strip_suffix(b"\n")
            .ok_or_else(|| mask_error(self.config, "fila de máscara sin LF"))?;
        let mut fields = row.split(|byte| *byte == b'\t');
        let (Some(id), Some(sequence_id), Some(start), Some(end)) =
            (fields.next(), fields.next(), fields.next(), fields.next())
        else {
            return Err(mask_error(self.config, "fila de máscara inválida"));
        };
        if fields.next().is_some()
            || [id, sequence_id, start, end]
                .iter()
                .any(|field| field.is_empty())
        {
            return Err(mask_error(self.config, "fila de máscara inválida"));
        }
        // `region_id` es una clave de orden después de coordenadas; dos filas
        // con las mismas coordenadas ya se rechazan por solape, por lo que no
        // hace falta retener su texto potencialmente enorme.
        std::str::from_utf8(id).map_err(|_| mask_error(self.config, "region_id no es UTF-8"))?;
        let sequence_id = std::str::from_utf8(sequence_id)
            .map_err(|_| mask_error(self.config, "sequence_id no es UTF-8"))?
            .to_owned();
        let start = parse_mask_u64(start, self.config)?;
        let end = parse_mask_u64(end, self.config)?;
        if start == 0 || start > end {
            return Err(mask_error(
                self.config,
                "coordenadas de máscara fuera de rango",
            ));
        }
        let current = MaskInterval {
            sequence_id,
            start,
            end,
        };
        if let Some(previous) = &self.previous
            && current.sequence_id == previous.sequence_id
            && (current.start < previous.start
                || (current.start == previous.start && current.end < previous.end)
                || current.start <= previous.end)
        {
            return Err(mask_error(self.config, "máscara desordenada o solapada"));
        }
        self.previous = Some(PreviousMaskInterval {
            sequence_id: current.sequence_id.clone(),
            start,
            end,
        });
        self.pending = Some(current);
        Ok(())
    }
    fn verify_digest(&mut self) -> Result<(), Diagnostic> {
        let Some(reader) = self.reader.take() else {
            return Ok(());
        };
        let (observed, bytes_processed) = reader.finish();
        self.bytes_processed = bytes_processed;
        let expected = self
            .config
            .declared_scenario()
            .excluded_regions_sha256
            .as_ref()
            .map(|value| value.value())
            .ok_or_else(|| mask_error(self.config, "hash de máscara ausente"))?;
        let observed = format!("{observed:x}");
        if observed != *expected {
            return Err(Diagnostic::input(
                "GPS-INPUT-002",
                "Una entrada cambió durante su copia o no coincide con su identidad declarada.",
            )
            .with(
                "path",
                self.config
                    .excluded_regions_path()
                    .expect("mask reader requires a path")
                    .display()
                    .to_string(),
            )
            .with("expected_digest", expected)
            .with("observed_digest", &observed));
        }
        // El escenario declara el SHA de la máscara, nunca su tamaño.
        self.identity = Some(AdmittedInputIdentity::new(
            PhysicalInputIdentity::from_observed(
                InputRole::Mask,
                self.config
                    .excluded_regions_path()
                    .expect("mask reader requires a path")
                    .to_path_buf(),
                bytes_processed,
                observed,
            ),
            None,
            Some(expected.clone()),
        )?);
        Ok(())
    }

    const fn bytes_processed(&self) -> u64 {
        self.bytes_processed
    }
}

const MAX_MASK_LINE_BYTES: usize = 65_536;

fn read_bounded_mask_line<R: Read>(
    reader: &mut R,
    line: &mut Vec<u8>,
    ledger: &mut F5ResourceLedger,
    config: &EffectiveF5Config,
) -> Result<usize, Diagnostic> {
    let initial_len = line.len();
    loop {
        if line.len() == MAX_MASK_LINE_BYTES {
            return Err(mask_error(config, "fila de máscara supera 65536 bytes"));
        }
        let mut byte = [0_u8; 1];
        let count = charged_read(
            reader,
            &mut byte,
            ledger,
            InputRole::Mask,
            1,
            config
                .excluded_regions_path()
                .expect("mask reader requires a path"),
        )?;
        if count == 0 {
            return Ok(line.len() - initial_len);
        }
        line.push(byte[0]);
        if byte[0] == b'\n' {
            return Ok(line.len() - initial_len);
        }
    }
}

struct HashingReader<R> {
    inner: R,
    digest: Sha256,
    bytes_processed: u64,
}
impl<R> HashingReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            digest: Sha256::new(),
            bytes_processed: 0,
        }
    }
    fn finish(self) -> (sha2::digest::Output<Sha256>, u64) {
        (self.digest.finalize(), self.bytes_processed)
    }
    const fn bytes_processed(&self) -> u64 {
        self.bytes_processed
    }
}
impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let count = self.inner.read(buffer)?;
        self.digest.update(&buffer[..count]);
        let counted_bytes = u64::try_from(count)
            .map_err(|_| std::io::Error::other("contador de bytes no cabe en u64"))?;
        self.bytes_processed = self
            .bytes_processed
            .checked_add(counted_bytes)
            .ok_or_else(|| std::io::Error::other("contador de bytes excede u64"))?;
        Ok(count)
    }
}

fn parse_mask_u64(value: &[u8], config: &EffectiveF5Config) -> Result<u64, Diagnostic> {
    if value.is_empty()
        || (value.len() > 1 && value.starts_with(b"0"))
        || !value.iter().all(u8::is_ascii_digit)
    {
        return Err(mask_error(config, "entero canónico de máscara inválido"));
    }
    std::str::from_utf8(value)
        .ok()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| mask_error(config, "entero de máscara fuera de u64"))
}
fn mask_error(config: &EffectiveF5Config, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-101",
        "El paquete de perfiles F5 no satisface su formato, versión, orden o identidad declarada.",
    )
    .with(
        "path",
        config
            .excluded_regions_path()
            .map_or_else(|| "<none>".to_owned(), |path| path.display().to_string()),
    )
    .with("field", "excluded_regions.tsv")
    .with("detail", detail.into())
}
fn domain(config: &EffectiveF5Config, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::domain(
        "GPS-DOMAIN-101",
        "La geometría F5 es demostrablemente inviable con el perfil y recursos declarados.",
    )
    .with("scenario_id", config.scenario_id().to_owned())
    .with("detail", detail.into())
}
fn f5_configuration(
    config: &EffectiveF5Config,
    fields: impl Into<String>,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::f5_configuration(config.scenario_path().display().to_string(), fields, detail)
}
fn resource(config: &EffectiveF5Config, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::resource(
        "GPS-RES-101",
        "Las dimensiones o presupuestos F5 exceden el límite declarado o duro.",
    )
    .with("scenario_id", config.scenario_id().to_owned())
    .with("detail", detail.into())
}
fn limit(value: Option<&crate::configuration::Located<u64>>, default: u64) -> u64 {
    value.map_or(default, |item| *item.value())
}
fn checked_usize(term: &str, value: usize, config: &EffectiveF5Config) -> Result<u64, Diagnostic> {
    u64::try_from(value).map_err(|_| resource(config, format!("{term} no cabe en u64")))
}
fn add(term: &str, left: u64, right: u64, config: &EffectiveF5Config) -> Result<u64, Diagnostic> {
    left.checked_add(right)
        .ok_or_else(|| resource(config, format!("overflow en {term}")))
}
fn add_many<const N: usize>(
    config: &EffectiveF5Config,
    term: &str,
    values: [u64; N],
) -> Result<u64, Diagnostic> {
    values
        .into_iter()
        .try_fold(0, |sum, value| add(term, sum, value, config))
}
fn sub(term: &str, left: u64, right: u64, config: &EffectiveF5Config) -> Result<u64, Diagnostic> {
    left.checked_sub(right)
        .ok_or_else(|| resource(config, format!("underflow en {term}")))
}
fn mul(term: &str, left: u64, right: u64, config: &EffectiveF5Config) -> Result<u64, Diagnostic> {
    left.checked_mul(right)
        .ok_or_else(|| resource(config, format!("overflow en {term}")))
}

fn triangular(term: &str, value: u64, config: &EffectiveF5Config) -> Result<u64, Diagnostic> {
    let previous = value
        .checked_sub(1)
        .ok_or_else(|| resource(config, format!("{term} con B=0")))?;
    if value.is_multiple_of(2) {
        mul(term, value / 2, previous, config)
    } else {
        mul(term, value, previous / 2, config)
    }
}
fn div(term: &str, left: u64, right: u64, config: &EffectiveF5Config) -> Result<u64, Diagnostic> {
    left.checked_div(right)
        .ok_or_else(|| resource(config, format!("división inválida en {term}")))
}
fn ceil_div(
    term: &str,
    value: u64,
    divisor: u64,
    config: &EffectiveF5Config,
) -> Result<u64, Diagnostic> {
    add(
        term,
        value / divisor,
        u64::from(!value.is_multiple_of(divisor)),
        config,
    )
}
