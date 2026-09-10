//! Longitudes de bloque F5 anteriores a la colocación física.
//!
//! F5.1.1 solo resuelve la geometría de las longitudes. No elige tramos,
//! posiciones ni SNPs: esos contratos pertenecen a F5.1.2 y posteriores.

use sha2::{Digest, Sha256};
#[cfg(test)]
use std::cell::RefCell;
use std::path::Path;

use crate::application::EffectiveF5Config;
use crate::diagnostics::Diagnostic;
use crate::f5_preflight::{
    F5_DEFAULT_MAX_LAYOUT_ATTEMPTS, F5EligibleTract, F5OtherResources, SealedF5Admission,
    resolved_target_for_config,
};
use crate::f5_profile::F5Profile;
use crate::identity::{
    ExecutionId, F5BlockId, F5IdentityContext, F5IdentityRegistry, F5SiteId, ReferenceId,
    ScenarioId, SequenceId,
};
use crate::randomness::{
    F5RandomDecision, RandomWordSource, RootSeed, derive_f5_subseed, f5_placed_layout_bytes,
    f5_uniform_index,
};
use crate::reference_fasta::stream_fasta_site_bases;

/// Algoritmo versionado de muestreo de longitudes del layout F5.
pub const F5_BLOCK_LENGTH_ALGORITHM: &str = "f5-block-lengths-v1";

/// Una fila del soporte entero de longitudes de bloque.
///
/// La distribución valida después que las filas sean contiguas y que sus pesos
/// sean probabilidades positivas y finitas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F5BlockLengthBand {
    minimum: u64,
    maximum: u64,
    weight: f64,
}

impl F5BlockLengthBand {
    #[must_use]
    pub const fn new(minimum: u64, maximum: u64, weight: f64) -> Self {
        Self {
            minimum,
            maximum,
            weight,
        }
    }

    #[must_use]
    pub const fn minimum(self) -> u64 {
        self.minimum
    }

    #[must_use]
    pub const fn maximum(self) -> u64 {
        self.maximum
    }

    #[must_use]
    pub const fn weight(self) -> f64 {
        self.weight
    }
}

/// Soporte discreto, continuo y ponderado de las longitudes de bloque.
#[derive(Clone, Debug, PartialEq)]
pub struct F5BlockLengthDistribution {
    bands: Vec<F5BlockLengthBand>,
    minimum: u64,
    maximum: u64,
}

impl F5BlockLengthDistribution {
    /// Construye el soporte validado independientemente del lector de perfiles.
    pub fn new(
        scenario_id: &ScenarioId,
        bands: Vec<F5BlockLengthBand>,
    ) -> Result<Self, Diagnostic> {
        let first = bands
            .first()
            .copied()
            .ok_or_else(|| domain(scenario_id, "perfil sin longitudes de bloque"))?;
        if first.minimum < 2 {
            return Err(domain(
                scenario_id,
                "la longitud mínima de bloque debe ser al menos dos",
            ));
        }
        let mut previous: Option<u64> = None;
        let mut weight_sum = 0.0_f64;
        let mut weight_compensation = 0.0_f64;
        for band in &bands {
            if band.minimum < 2
                || band.minimum > band.maximum
                || !band.weight.is_finite()
                || band.weight <= 0.0
            {
                return Err(domain(scenario_id, "intervalo o peso de longitud inválido"));
            }
            if let Some(last) = previous
                && last.checked_add(1) != Some(band.minimum)
            {
                return Err(domain(
                    scenario_id,
                    "el soporte de longitudes no es continuo",
                ));
            }
            // Neumaier, igual que el contrato de perfiles, para no aceptar una
            // distribución distinta de una probabilidad restringida.
            let next = weight_sum + band.weight;
            if weight_sum.abs() >= band.weight.abs() {
                weight_compensation += (weight_sum - next) + band.weight;
            } else {
                weight_compensation += (band.weight - next) + weight_sum;
            }
            weight_sum = next;
            previous = Some(band.maximum);
        }
        let normalized_sum = weight_sum + weight_compensation;
        if (normalized_sum - 1.0).abs() > 1e-9 {
            return Err(domain(
                scenario_id,
                "los pesos de longitudes no suman uno dentro de tolerancia",
            ));
        }
        Ok(Self {
            minimum: first.minimum,
            maximum: previous.expect("el primer intervalo ya comprobó que existe"),
            bands,
        })
    }

    pub fn from_profile(scenario_id: &ScenarioId, profile: &F5Profile) -> Result<Self, Diagnostic> {
        Self::new(
            scenario_id,
            profile
                .block_lengths()
                .iter()
                .map(|band| F5BlockLengthBand::new(band.minimum(), band.maximum(), band.weight()))
                .collect(),
        )
    }

    #[must_use]
    pub fn bands(&self) -> &[F5BlockLengthBand] {
        &self.bands
    }

    #[must_use]
    pub const fn minimum(&self) -> u64 {
        self.minimum
    }

    #[must_use]
    pub const fn maximum(&self) -> u64 {
        self.maximum
    }
}

/// Prueba entera previa de que un objetivo T puede partirse en B longitudes.
///
/// No usa pesos ni RNG. Su resultado es un testigo reconstruible para separar
/// una imposibilidad matemática (`GPS-DOMAIN-101`) de un intento de layout que
/// después no encuentre posiciones físicas.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F5LengthFeasibility {
    target_bases: u64,
    block_count: u64,
    minimum_total: u64,
    maximum_total: u64,
}

impl F5LengthFeasibility {
    pub fn for_target(
        scenario_id: &ScenarioId,
        target_bases: u64,
        distribution: &F5BlockLengthDistribution,
    ) -> Result<Self, Diagnostic> {
        if target_bases == 0 {
            return Err(domain(scenario_id, "target_bases debe ser positivo"));
        }
        let block_count = ceil_div(target_bases, distribution.maximum)
            .ok_or_else(|| domain(scenario_id, "B no cabe en u64"))?;
        Self::with_block_count(scenario_id, target_bases, block_count, distribution)
    }

    pub fn with_block_count(
        scenario_id: &ScenarioId,
        target_bases: u64,
        block_count: u64,
        distribution: &F5BlockLengthDistribution,
    ) -> Result<Self, Diagnostic> {
        if target_bases == 0 || block_count == 0 {
            return Err(domain(
                scenario_id,
                "target_bases y block_count deben ser positivos",
            ));
        }
        let minimum_total = block_count
            .checked_mul(distribution.minimum)
            .ok_or_else(|| domain(scenario_id, "B*Lmin desborda u64"))?;
        let maximum_total = block_count
            .checked_mul(distribution.maximum)
            .ok_or_else(|| domain(scenario_id, "B*Lmax desborda u64"))?;
        if !(minimum_total..=maximum_total).contains(&target_bases) {
            return Err(domain(
                scenario_id,
                "target_bases no es representable por el soporte de longitudes",
            ));
        }
        Ok(Self {
            target_bases,
            block_count,
            minimum_total,
            maximum_total,
        })
    }

    #[must_use]
    pub const fn target_bases(self) -> u64 {
        self.target_bases
    }

    #[must_use]
    pub const fn block_count(self) -> u64 {
        self.block_count
    }

    #[must_use]
    pub const fn minimum_total(self) -> u64 {
        self.minimum_total
    }

    #[must_use]
    pub const fn maximum_total(self) -> u64 {
        self.maximum_total
    }
}

/// Plan de longitudes: condicionado por T o estocástico por número de bloques.
#[derive(Clone, Debug, PartialEq)]
pub struct F5BlockLengthPlan {
    scenario_id: ScenarioId,
    distribution: F5BlockLengthDistribution,
    block_count: u64,
    target_bases: Option<u64>,
}

impl F5BlockLengthPlan {
    pub fn from_preflight(
        config: &EffectiveF5Config,
        preflight: &SealedF5Admission,
    ) -> Result<Self, Diagnostic> {
        let scenario_id = ScenarioId::parse(config.scenario_id(), "scenario_id").map_err(|_| {
            configuration(
                config,
                "scenario_id",
                "no se pudo reconstruir como identidad tipada",
            )
        })?;
        let profile = config.profile();
        let distribution =
            F5BlockLengthDistribution::from_profile(&scenario_id, profile).map_err(|error| {
                configuration(
                    config,
                    "f5.profile_metadata_path",
                    format!(
                        "el perfil admitido no conserva su distribución de longitudes: {}",
                        error.context()["detail"]
                    ),
                )
            })?;
        let declared_profile_sha256 = config.declared_scenario().profile_metadata_sha256.value();
        let resolved_target = resolved_target_for_config(config)?;
        if !preflight.matches_configuration(config)
            || preflight.scenario_id() != config.scenario_id()
            || preflight.reference_id() != config.reference_id()
            || preflight.reference_sha256() != config.reference_sha256()
            || preflight.profile_id() != profile.profile_id()
            || preflight.profile_metadata_sha256() != declared_profile_sha256
            || preflight.block_length_minimum() != distribution.minimum
            || preflight.block_length_maximum() != distribution.maximum
            || preflight.resolved_target() != &resolved_target
        {
            return Err(configuration(
                config,
                "preflight.binding",
                "el preflight no pertenece a la configuración, objetivo y perfil de longitudes recibidos",
            ));
        }
        let block_count = preflight.block_count();
        let target_bases = preflight.target_bases();
        if let Some(target) = target_bases {
            let oracle = F5LengthFeasibility::for_target(&scenario_id, target, &distribution)?;
            if oracle.block_count != block_count {
                return Err(domain(
                    &scenario_id,
                    "el block_count del preflight no es el mínimo factible para target_bases",
                ));
            }
        } else if block_count == 0 {
            return Err(domain(
                &scenario_id,
                "el preflight admitió block_count cero",
            ));
        }
        Ok(Self {
            scenario_id,
            distribution,
            block_count,
            target_bases,
        })
    }

    #[must_use]
    pub const fn algorithm(&self) -> &'static str {
        F5_BLOCK_LENGTH_ALGORITHM
    }

    #[must_use]
    pub const fn scenario_id(&self) -> &ScenarioId {
        &self.scenario_id
    }

    #[must_use]
    pub const fn block_count(&self) -> u64 {
        self.block_count
    }

    #[must_use]
    pub const fn target_bases(&self) -> Option<u64> {
        self.target_bases
    }

    #[must_use]
    pub const fn distribution(&self) -> &F5BlockLengthDistribution {
        &self.distribution
    }

    /// Sortea una longitud por bloque para un intento de layout.
    ///
    /// En presencia de T, los B-1 primeros sorteos se restringen a valores que
    /// dejan un residual factible y el último bloque recibe ese residual sin
    /// abrir su ámbito RNG. En `block_count`, cada bloque usa la distribución
    /// completa e independiente.
    /// Sortea las longitudes de un intento.
    ///
    /// No es público fuera del crate: `from_preflight` solo exige la mitad
    /// semántica del binding, de modo que una cadena pública hasta aquí
    /// permitiría abrir un ámbito RNG con una configuración físicamente
    /// distinta de la que se admitió. La vía pública al sorteo es
    /// `F5LayoutRequest`, que exige el receipt sellado.
    pub(crate) fn sample(
        &self,
        root_seed: RootSeed,
        attempt_index: u64,
    ) -> Result<Vec<u64>, Diagnostic> {
        let scenario_id = &self.scenario_id;
        let capacity =
            usize::try_from(self.block_count).map_err(|_| allocation_resource(self.block_count))?;
        let mut lengths = Vec::new();
        lengths
            .try_reserve_exact(capacity)
            .map_err(|_| allocation_resource(self.block_count))?;
        let mut remaining = self.target_bases;
        for block_index in 0..self.block_count {
            if let Some(target) = remaining {
                let remaining_blocks = self.block_count - block_index - 1;
                if remaining_blocks == 0 {
                    if target < self.distribution.minimum || target > self.distribution.maximum {
                        return Err(domain(
                            scenario_id,
                            "residual final fuera del soporte factible",
                        ));
                    }
                    lengths.push(target);
                    continue;
                }
                let lower = target
                    .checked_sub(
                        remaining_blocks
                            .checked_mul(self.distribution.maximum)
                            .ok_or_else(|| {
                                internal(scenario_id.as_str(), "residual*Lmax desborda u64")
                            })?,
                    )
                    .unwrap_or(self.distribution.minimum)
                    .max(self.distribution.minimum);
                let upper = target
                    .saturating_sub(
                        remaining_blocks
                            .checked_mul(self.distribution.minimum)
                            .ok_or_else(|| {
                                internal(scenario_id.as_str(), "residual*Lmin desborda u64")
                            })?,
                    )
                    .min(self.distribution.maximum);
                let value = self.sample_between(
                    root_seed,
                    scenario_id,
                    attempt_index,
                    block_index,
                    lower,
                    upper,
                )?;
                remaining = Some(
                    target
                        .checked_sub(value)
                        .ok_or_else(|| domain(scenario_id, "residual negativo"))?,
                );
                lengths.push(value);
            } else {
                lengths.push(self.sample_between(
                    root_seed,
                    scenario_id,
                    attempt_index,
                    block_index,
                    self.distribution.minimum,
                    self.distribution.maximum,
                )?);
            }
        }
        if let Some(target) = self.target_bases {
            let sum = lengths
                .iter()
                .try_fold(0_u64, |total, value| total.checked_add(*value))
                .ok_or_else(|| internal(scenario_id.as_str(), "suma de longitudes desborda u64"))?;
            if sum != target {
                return Err(domain(
                    scenario_id,
                    "la suma condicionada no coincide con target_bases",
                ));
            }
        }
        Ok(lengths)
    }

    fn sample_between(
        &self,
        root_seed: RootSeed,
        scenario_id: &ScenarioId,
        attempt_index: u64,
        block_index: u64,
        lower: u64,
        upper: u64,
    ) -> Result<u64, Diagnostic> {
        if lower > upper {
            return Err(domain(
                scenario_id,
                "el oráculo permitió un intervalo residual vacío",
            ));
        }
        let mut rng = derive_f5_subseed(
            root_seed,
            F5RandomDecision::layout_length(scenario_id, attempt_index, block_index),
        )?
        .local_rng();
        self.sample_between_with_source(lower, upper, &mut rng)
    }

    fn sample_between_with_source(
        &self,
        lower: u64,
        upper: u64,
        rng: &mut impl RandomWordSource,
    ) -> Result<u64, Diagnostic> {
        choose_weighted_integer_with_source(&self.distribution.bands, lower, upper, rng).map_err(
            |error| match error {
                WeightedIntegerError::Domain(detail) => domain(&self.scenario_id, detail),
                WeightedIntegerError::Internal(detail) => {
                    internal(self.scenario_id.as_str(), detail)
                }
            },
        )
    }
}

/// Distribución discreta de distancias entre SNP consecutivos dentro de un
/// bloque. Comparte la semántica U53 de las longitudes, pero admite `1` como
/// distancia mínima.
#[derive(Clone, Debug, PartialEq)]
pub struct F5SnpGapDistribution {
    bands: Vec<F5BlockLengthBand>,
}

impl F5SnpGapDistribution {
    fn from_bands(
        scenario_id: &ScenarioId,
        bands: Vec<F5BlockLengthBand>,
    ) -> Result<Self, Diagnostic> {
        validate_weighted_integer_bands(scenario_id, &bands, 1, "gaps SNP")?;
        Ok(Self { bands })
    }

    /// Constructor reservado a las pruebas de la selección discreta.
    #[cfg(test)]
    pub fn new(
        scenario_id: &ScenarioId,
        bands: Vec<F5BlockLengthBand>,
    ) -> Result<Self, Diagnostic> {
        Self::from_bands(scenario_id, bands)
    }

    fn from_profile(scenario_id: &ScenarioId, profile: &F5Profile) -> Result<Self, Diagnostic> {
        Self::from_bands(
            scenario_id,
            profile
                .snp_gaps()
                .iter()
                .map(|gap| F5BlockLengthBand::new(gap.minimum(), gap.maximum(), gap.weight()))
                .collect(),
        )
    }

    #[must_use]
    pub fn bands(&self) -> &[F5BlockLengthBand] {
        &self.bands
    }

    fn sample(
        &self,
        scenario_id: &ScenarioId,
        rng: &mut impl RandomWordSource,
    ) -> Result<u64, Diagnostic> {
        choose_weighted_integer_with_source(&self.bands, 1, u64::MAX, rng).map_err(|error| {
            match error {
                WeightedIntegerError::Domain(detail) => {
                    Diagnostic::f5_internal(scenario_id.as_str(), detail)
                }
                WeightedIntegerError::Internal(detail) => {
                    Diagnostic::f5_internal(scenario_id.as_str(), detail)
                }
            }
        })
    }
}

/// Secuencia de referencia ya verificada que permite recuperar la base REF sin
/// volver a escanear elegibilidad ni máscara.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5LayoutSequence {
    sequence_ordinal: u64,
    sequence_id: SequenceId,
    bases: Vec<u8>,
}

impl F5LayoutSequence {
    /// Constructor de pruebas de la mecánica de layout. La frontera de
    /// producción recupera únicamente las bases REF elegidas mediante
    /// streaming; no reconstruye las secuencias completas.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn new(sequence_ordinal: u64, sequence_id: SequenceId, bases: Vec<u8>) -> Self {
        Self {
            sequence_ordinal,
            sequence_id,
            bases,
        }
    }

    #[must_use]
    pub const fn sequence_ordinal(&self) -> u64 {
        self.sequence_ordinal
    }

    #[must_use]
    pub fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }
}

/// Sitio SNP común a todas las subpoblaciones, todavía sin marginales ni
/// haplotipos (F5.1.3+).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5LayoutSite {
    site_index: u64,
    site_id: F5SiteId,
    position_1based: u64,
    reference_base: u8,
    alternate_base: u8,
}

impl F5LayoutSite {
    #[must_use]
    pub const fn site_index(&self) -> u64 {
        self.site_index
    }
    #[must_use]
    pub fn site_id(&self) -> &F5SiteId {
        &self.site_id
    }
    #[must_use]
    pub const fn position_1based(&self) -> u64 {
        self.position_1based
    }
    #[must_use]
    pub const fn reference_base(&self) -> u8 {
        self.reference_base
    }
    #[must_use]
    pub const fn alternate_base(&self) -> u8 {
        self.alternate_base
    }
}

/// Bloque colocado en orden canónico FASTA/inicio/fin.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5LayoutBlock {
    block_index: u64,
    block_id: F5BlockId,
    sequence_ordinal: u64,
    sequence_id: SequenceId,
    start_1based: u64,
    end_1based: u64,
    length_bp: u64,
    sites: Vec<F5LayoutSite>,
}

impl F5LayoutBlock {
    #[must_use]
    pub const fn block_index(&self) -> u64 {
        self.block_index
    }
    #[must_use]
    pub fn block_id(&self) -> &F5BlockId {
        &self.block_id
    }
    #[must_use]
    pub const fn sequence_ordinal(&self) -> u64 {
        self.sequence_ordinal
    }
    #[must_use]
    pub fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }
    #[must_use]
    pub const fn start_1based(&self) -> u64 {
        self.start_1based
    }
    #[must_use]
    pub const fn end_1based(&self) -> u64 {
        self.end_1based
    }
    #[must_use]
    pub const fn length_bp(&self) -> u64 {
        self.length_bp
    }
    #[must_use]
    pub fn sites(&self) -> &[F5LayoutSite] {
        &self.sites
    }
}

/// Resultado completo de F5.1.2; los siguientes hitos consumen estos bloques
/// y sitios, no vuelven a decidir su geometría.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5PhysicalLayout {
    blocks: Vec<F5LayoutBlock>,
    successful_attempt_index: u64,
}

/// Resumen canónico del mejor intento de layout agotado.
///
/// Es el valor que el publicador necesita para construir `f5-failure-v1`; el
/// hash se calcula sobre los bytes canónicos que también participan en el
/// desempate y nunca se recupera analizando texto diagnóstico.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5BestLayout {
    attempt_index: u64,
    placed_blocks: u64,
    unplaced_blocks: u64,
    placed_bases: u64,
    remaining_target_bases: u64,
    layout_bytes_sha256: String,
}

impl F5BestLayout {
    #[must_use]
    pub const fn attempt_index(&self) -> u64 {
        self.attempt_index
    }
    #[must_use]
    pub const fn placed_blocks(&self) -> u64 {
        self.placed_blocks
    }
    #[must_use]
    pub const fn unplaced_blocks(&self) -> u64 {
        self.unplaced_blocks
    }
    #[must_use]
    pub const fn placed_bases(&self) -> u64 {
        self.placed_bases
    }
    #[must_use]
    pub const fn remaining_target_bases(&self) -> u64 {
        self.remaining_target_bases
    }
    #[must_use]
    pub fn layout_bytes_sha256(&self) -> &str {
        &self.layout_bytes_sha256
    }
}

/// Error de layout que conserva el testigo tipado de un agotamiento.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5LayoutFailure {
    diagnostic: Box<Diagnostic>,
    best_layout: Option<F5BestLayout>,
}

impl F5LayoutFailure {
    #[must_use]
    pub const fn diagnostic(&self) -> &Diagnostic {
        &self.diagnostic
    }
    #[must_use]
    pub const fn best_layout(&self) -> Option<&F5BestLayout> {
        self.best_layout.as_ref()
    }
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.diagnostic.code()
    }
    #[must_use]
    pub fn context(&self) -> &std::collections::BTreeMap<String, String> {
        self.diagnostic.context()
    }
    #[must_use]
    pub fn into_diagnostic(self) -> Diagnostic {
        *self.diagnostic
    }

    fn exhausted(scenario_id: &str, best_layout: F5BestLayout) -> Self {
        let detail = format!(
            "best_attempt_index={};unplaced_blocks={};remaining_target_bases={}",
            best_layout.attempt_index,
            best_layout.unplaced_blocks,
            best_layout.remaining_target_bases,
        );
        Self {
            diagnostic: Box::new(Diagnostic::f5_layout_exhausted(scenario_id, detail)),
            best_layout: Some(best_layout),
        }
    }
}

impl From<Diagnostic> for F5LayoutFailure {
    fn from(diagnostic: Diagnostic) -> Self {
        Self {
            diagnostic: Box::new(diagnostic),
            best_layout: None,
        }
    }
}

impl F5PhysicalLayout {
    #[must_use]
    pub fn blocks(&self) -> &[F5LayoutBlock] {
        &self.blocks
    }
    #[must_use]
    pub const fn successful_attempt_index(&self) -> u64 {
        self.successful_attempt_index
    }
}

/// Frontera de producción para el layout F5.1.2.
///
/// Solo enlaza la configuración efectiva, el preflight exacto que la admitió y
/// la identidad de ejecución. El consumidor deriva de esas fuentes los tramos
/// maximales, plan de longitudes, gaps, límite de intentos, semilla e
/// identidades; por tanto no puede sustituirlos por fragmentos o valores
/// compatibles solo de forma accidental.
#[derive(Debug)]
pub struct F5LayoutRequest<'a> {
    config: Option<&'a EffectiveF5Config>,
    preflight: Option<&'a SealedF5Admission>,
    execution_id: Option<&'a ExecutionId>,
    #[cfg(test)]
    test_inputs: Option<F5LayoutInputs<'a>>,
}

impl<'a> F5LayoutRequest<'a> {
    /// Rechaza antes de abrir cualquier ámbito RNG si el preflight no es el
    /// del escenario, referencia, perfil, objetivo o límite de intentos.
    pub fn new(
        config: &'a EffectiveF5Config,
        preflight: &'a SealedF5Admission,
        execution_id: &'a ExecutionId,
    ) -> Result<Self, Diagnostic> {
        Self::new_with_other_resources(config, preflight, F5OtherResources::default(), execution_id)
    }

    /// Construye una solicitud de producción cuando la ejecución incorpora
    /// recursos ajenos a F5. El valor debe ser exactamente el que participó
    /// en el preflight; así ningún consumidor puede reutilizar una admisión
    /// de menor coste antes de abrir RNG.
    pub fn new_with_other_resources(
        config: &'a EffectiveF5Config,
        preflight: &'a SealedF5Admission,
        other: F5OtherResources,
        execution_id: &'a ExecutionId,
    ) -> Result<Self, Diagnostic> {
        if !preflight.matches_admission(config, other) {
            return Err(configuration(
                config,
                "preflight.binding",
                "el preflight no pertenece a la configuración efectiva y aportaciones de recursos recibidas",
            ));
        }
        F5BlockLengthPlan::from_preflight(config, preflight)?;
        let attempts = config
            .declared_scenario()
            .max_layout_attempts
            .as_ref()
            .map_or(F5_DEFAULT_MAX_LAYOUT_ATTEMPTS, |value| *value.value());
        if preflight.max_layout_attempts() != attempts {
            return Err(configuration(
                config,
                "preflight.max_layout_attempts",
                "el preflight no pertenece al límite de intentos efectivo",
            ));
        }
        RootSeed::parse("seed", config.seed()).map_err(|_| {
            configuration(
                config,
                "seed",
                "la semilla efectiva admitida no es canónica",
            )
        })?;
        ReferenceId::parse(config.reference_id(), "f5_layout").map_err(|_| {
            configuration(
                config,
                "reference_id",
                "la identidad de referencia efectiva no es válida",
            )
        })?;
        Ok(Self {
            config: Some(config),
            preflight: Some(preflight),
            execution_id: Some(execution_id),
            #[cfg(test)]
            test_inputs: None,
        })
    }

    #[cfg(test)]
    fn from_test_parts(inputs: F5LayoutInputs<'a>) -> Self {
        // Los tests de la mecánica aislada son la única vía que puede aportar
        // piezas sueltas. Esta variante no existe en artefactos de producción.
        Self {
            config: None,
            preflight: None,
            execution_id: None,
            test_inputs: Some(inputs),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct F5LayoutInputs<'a> {
    length_plan: &'a F5BlockLengthPlan,
    eligible_tracts: &'a [F5EligibleTract],
    gap_distribution: &'a F5SnpGapDistribution,
    sequences: Option<&'a [F5LayoutSequence]>,
    streamed_reference: Option<F5StreamedReference<'a>>,
    max_layout_attempts: u64,
    root_seed: RootSeed,
    identity_context: F5IdentityContext<'a>,
}

#[derive(Clone, Copy, Debug)]
struct F5StreamedReference<'a> {
    path: &'a Path,
    reference_id: &'a str,
    /// Tamaño y SHA observados por la primera pasada y sellados en el receipt.
    /// La segunda pasada exige ambos, no la copia declarada en el general.
    sealed_size_bytes: u64,
    sealed_sha256: &'a str,
    expected_bases_processed: u64,
}

struct ResolvedProductionLayout<'a> {
    length_plan: F5BlockLengthPlan,
    gap_distribution: F5SnpGapDistribution,
    root_seed: RootSeed,
    reference_id: ReferenceId,
    preflight: &'a SealedF5Admission,
    reference_path: &'a Path,
    sealed_reference_size_bytes: u64,
    sealed_reference_sha256: &'a str,
}

impl<'a> ResolvedProductionLayout<'a> {
    fn from_request(request: &F5LayoutRequest<'a>) -> Result<Self, Diagnostic> {
        let config = request.config.ok_or_else(|| {
            Diagnostic::internal_invariant("f5_layout", "production request without config")
        })?;
        let preflight = request.preflight.ok_or_else(|| {
            Diagnostic::internal_invariant("f5_layout", "production request without preflight")
        })?;
        let length_plan = F5BlockLengthPlan::from_preflight(config, preflight)?;
        let scenario_id = length_plan.scenario_id();
        let gap_distribution = F5SnpGapDistribution::from_profile(scenario_id, config.profile())?;
        let root_seed = RootSeed::parse("seed", config.seed()).map_err(|_| {
            configuration(
                config,
                "seed",
                "la semilla efectiva admitida no es canónica",
            )
        })?;
        let reference_id =
            ReferenceId::parse(config.reference_id(), "f5_layout").map_err(|_| {
                configuration(
                    config,
                    "reference_id",
                    "la identidad de referencia efectiva no es válida",
                )
            })?;
        let sealed_reference = preflight.receipt().reference_identity()?;
        Ok(Self {
            length_plan,
            gap_distribution,
            root_seed,
            reference_id,
            preflight,
            // La ruta también procede del receipt: la segunda pasada recorre
            // exactamente el fichero que la admisión abrió, no el que una
            // configuración posterior declare.
            reference_path: sealed_reference.normalized_path(),
            sealed_reference_size_bytes: sealed_reference.observed_size_bytes(),
            sealed_reference_sha256: sealed_reference.observed_sha256(),
        })
    }

    fn inputs<'b>(&'b self, execution_id: &'b ExecutionId) -> F5LayoutInputs<'b> {
        F5LayoutInputs {
            length_plan: &self.length_plan,
            eligible_tracts: self.preflight.eligible_tracts(),
            gap_distribution: &self.gap_distribution,
            sequences: None,
            streamed_reference: Some(F5StreamedReference {
                path: self.reference_path,
                reference_id: self.reference_id.as_str(),
                sealed_size_bytes: self.sealed_reference_size_bytes,
                sealed_sha256: self.sealed_reference_sha256,
                expected_bases_processed: self.preflight.reference_bases_processed(),
            }),
            max_layout_attempts: self.preflight.max_layout_attempts(),
            root_seed: self.root_seed,
            identity_context: F5IdentityContext::new(
                self.root_seed,
                execution_id,
                self.length_plan.scenario_id(),
                &self.reference_id,
            ),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct FreeTract {
    sequence_ordinal: u64,
    start_1based: u64,
    end_1based: u64,
}

#[cfg(test)]
#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct LayoutInstrumentation {
    request_tract_visits: u64,
    total_tract_visits: u64,
    prefix_tract_visits: u64,
    initial_tract_copies: u64,
    reservation_calls: u64,
    reserved_slots: u64,
    fragment_growth: u64,
    divisions: u64,
    remove_shifts: u64,
    insert_shifts: u64,
    attempts_started: u64,
    retries: u64,
    peak_live_slots: u64,
    peak_live_bytes: u64,
    placement_trace: Vec<(u64, usize, u64, bool, bool)>,
}

#[cfg(test)]
thread_local! {
    static LAYOUT_INSTRUMENTATION: RefCell<LayoutInstrumentation> = RefCell::new(LayoutInstrumentation::default());
}

fn record_request_tract_visit() {
    #[cfg(test)]
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.request_tract_visits =
            instrumentation.request_tract_visits.saturating_add(1);
    });
}

fn record_total_tract_visit() {
    #[cfg(test)]
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.total_tract_visits = instrumentation.total_tract_visits.saturating_add(1);
    });
}

fn record_prefix_tract_visit() {
    #[cfg(test)]
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.prefix_tract_visits = instrumentation.prefix_tract_visits.saturating_add(1);
    });
}

fn record_initial_tract_copy() {
    #[cfg(test)]
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.initial_tract_copies =
            instrumentation.initial_tract_copies.saturating_add(1);
    });
}

fn record_remove_shift(count: usize) {
    #[cfg(not(test))]
    let _ = count;
    #[cfg(test)]
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.remove_shifts = instrumentation
            .remove_shifts
            .saturating_add(u64::try_from(count).unwrap_or(u64::MAX));
    });
}

fn record_reserved_slots(count: usize) {
    #[cfg(not(test))]
    let _ = count;
    #[cfg(test)]
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.reservation_calls = instrumentation.reservation_calls.saturating_add(1);
        instrumentation.reserved_slots = instrumentation
            .reserved_slots
            .saturating_add(u64::try_from(count).unwrap_or(u64::MAX));
    });
}

fn record_insert_shift(count: usize) {
    #[cfg(not(test))]
    let _ = count;
    #[cfg(test)]
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.insert_shifts = instrumentation
            .insert_shifts
            .saturating_add(u64::try_from(count).unwrap_or(u64::MAX));
    });
}

#[cfg(test)]
fn record_fragment_growth() {
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.fragment_growth = instrumentation.fragment_growth.saturating_add(1);
    });
}

#[cfg(test)]
fn record_division() {
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.divisions = instrumentation.divisions.saturating_add(1);
    });
}

#[cfg(test)]
fn record_attempt_started() {
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.attempts_started = instrumentation.attempts_started.saturating_add(1);
    });
}

#[cfg(test)]
fn record_retry() {
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        instrumentation.retries = instrumentation.retries.saturating_add(1);
    });
}

#[cfg(test)]
fn record_live_slots(slots: usize) {
    LAYOUT_INSTRUMENTATION.with(|instrumentation| {
        let mut instrumentation = instrumentation.borrow_mut();
        let slots = u64::try_from(slots).unwrap_or(u64::MAX);
        instrumentation.peak_live_slots = instrumentation.peak_live_slots.max(slots);
        instrumentation.peak_live_bytes = instrumentation.peak_live_slots.saturating_mul(256);
    });
}

#[cfg(test)]
fn reset_layout_instrumentation() {
    LAYOUT_INSTRUMENTATION
        .with(|instrumentation| *instrumentation.borrow_mut() = LayoutInstrumentation::default());
}

#[cfg(test)]
fn layout_instrumentation() -> LayoutInstrumentation {
    LAYOUT_INSTRUMENTATION.with(|instrumentation| instrumentation.borrow().clone())
}

// Adaptadores para las regresiones unitarias preexistentes; la prueba
// instrumental de producción usa el desglose completo de `LayoutInstrumentation`.
#[cfg(test)]
fn reset_tract_query_count() {
    reset_layout_instrumentation();
}

#[cfg(test)]
fn tract_query_count() -> u64 {
    let observed = layout_instrumentation();
    observed.total_tract_visits + observed.prefix_tract_visits
}

#[cfg(test)]
fn layout_vector_moves() -> u64 {
    let observed = layout_instrumentation();
    observed.initial_tract_copies + observed.remove_shifts + observed.insert_shifts
}

#[cfg(test)]
fn layout_reserved_slots() -> u64 {
    layout_instrumentation().reserved_slots
}

#[cfg(test)]
fn layout_placement_trace() -> Vec<(u64, usize, u64, bool, bool)> {
    layout_instrumentation().placement_trace
}

impl FreeTract {
    const fn len(self) -> u64 {
        self.end_1based - self.start_1based + 1
    }
}

#[derive(Clone, Debug)]
struct PlacedBlockDraft {
    source_block_index: u64,
    length_bp: u64,
    sequence_ordinal: u64,
    start_1based: u64,
    end_1based: u64,
}

/// Coloca bloques y sitios de F5.1.2. Cada intento construye sus propios
/// tramos y borradores; ninguno alcanza el resultado hasta tener todos los
/// bloques y al menos dos SNP por bloque.
pub fn build_f5_physical_layout(
    request: F5LayoutRequest<'_>,
) -> Result<F5PhysicalLayout, F5LayoutFailure> {
    #[cfg(test)]
    let production = request
        .test_inputs
        .is_none()
        .then(|| ResolvedProductionLayout::from_request(&request))
        .transpose()?;
    #[cfg(not(test))]
    let production = ResolvedProductionLayout::from_request(&request)?;
    #[cfg(test)]
    let inputs = if let Some(inputs) = request.test_inputs {
        inputs
    } else {
        let production = production.as_ref().ok_or_else(|| {
            Diagnostic::internal_invariant("f5_layout", "production layout not resolved")
        })?;
        let execution_id = request.execution_id.ok_or_else(|| {
            Diagnostic::internal_invariant("f5_layout", "production request without execution")
        })?;
        production.inputs(execution_id)
    };
    #[cfg(not(test))]
    let inputs = {
        let execution_id = request.execution_id.ok_or_else(|| {
            Diagnostic::internal_invariant("f5_layout", "production request without execution")
        })?;
        production.inputs(execution_id)
    };
    validate_layout_request(&inputs)?;
    let scenario_id = inputs.length_plan.scenario_id();
    let mut best: Option<BestLayoutCandidate> = None;

    for attempt_index in 0..inputs.max_layout_attempts {
        #[cfg(test)]
        record_attempt_started();
        let lengths = inputs.length_plan.sample(inputs.root_seed, attempt_index)?;
        let planned_bases = lengths.iter().try_fold(0_u64, |sum, length| {
            sum.checked_add(*length).ok_or_else(|| {
                Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "la suma de longitudes del intento desborda u64",
                )
            })
        })?;
        let mut pending = lengths
            .into_iter()
            .enumerate()
            .map(|(index, length_bp)| {
                Ok((
                    u64::try_from(index).map_err(|_| {
                        Diagnostic::f5_internal(scenario_id.as_str(), "block_index no cabe en u64")
                    })?,
                    length_bp,
                ))
            })
            .collect::<Result<Vec<_>, Diagnostic>>()?;
        #[cfg(test)]
        record_live_slots(
            inputs
                .eligible_tracts
                .len()
                .saturating_add(pending.capacity()),
        );
        pending.sort_unstable_by(|left, right| {
            right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0))
        });
        #[cfg(test)]
        let pending_capacity = pending.capacity();
        let free_capacity = inputs
            .eligible_tracts
            .len()
            .checked_add(pending.len())
            .ok_or_else(|| {
                Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "la reserva de tramos libres desborda usize",
                )
            })?;
        record_reserved_slots(free_capacity);
        let mut free = Vec::with_capacity(free_capacity);
        for tract in inputs.eligible_tracts {
            free.push(FreeTract {
                sequence_ordinal: tract.sequence_ordinal(),
                start_1based: tract.start_1based(),
                end_1based: tract.end_1based(),
            });
            record_initial_tract_copy();
        }
        let mut placed = Vec::with_capacity(pending.len());
        record_reserved_slots(pending.len());
        #[cfg(test)]
        record_live_slots(
            inputs
                .eligible_tracts
                .len()
                .saturating_add(pending_capacity)
                .saturating_add(free.capacity())
                .saturating_add(placed.capacity()),
        );
        let mut placement_failed = false;

        for (source_block_index, length_bp) in pending {
            let total_starts = total_starts_for_length(&free, length_bp).ok_or_else(|| {
                Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "el total de inicios válidos desborda u64",
                )
            })?;
            if total_starts == 0 {
                placement_failed = true;
                break;
            }
            let mut rng = derive_f5_subseed(
                inputs.root_seed,
                F5RandomDecision::layout_position(scenario_id, attempt_index, source_block_index),
            )?
            .local_rng();
            let chosen = f5_uniform_index(&mut rng, total_starts).map_err(|_| {
                Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "el selector uniforme de posición recibió un límite inválido",
                )
            })?;
            let (tract_index, start_1based) = resolve_start_by_prefix(&free, length_bp, chosen)
                .ok_or_else(|| {
                    Diagnostic::f5_internal(
                        scenario_id.as_str(),
                        "el prefijo de inicios no resolvió una posición válida",
                    )
                })?;
            #[cfg(test)]
            LAYOUT_INSTRUMENTATION.with(|instrumentation| {
                let tract = free[tract_index];
                instrumentation.borrow_mut().placement_trace.push((
                    attempt_index,
                    tract_index,
                    start_1based,
                    tract.start_1based < start_1based,
                    start_1based
                        .checked_add(length_bp - 1)
                        .is_some_and(|end| end < tract.end_1based),
                ));
            });
            record_remove_shift(free.len().saturating_sub(tract_index + 1));
            let tract = free.remove(tract_index);
            let end_1based = start_1based.checked_add(length_bp - 1).ok_or_else(|| {
                Diagnostic::f5_internal(scenario_id.as_str(), "el final de bloque desborda u64")
            })?;
            let mut insert_at = tract_index;
            #[cfg(test)]
            let mut inserted_residues = 0_u64;
            if tract.start_1based < start_1based {
                record_insert_shift(free.len().saturating_sub(insert_at));
                free.insert(
                    insert_at,
                    FreeTract {
                        sequence_ordinal: tract.sequence_ordinal,
                        start_1based: tract.start_1based,
                        end_1based: start_1based - 1,
                    },
                );
                #[cfg(test)]
                {
                    inserted_residues += 1;
                }
                insert_at += 1;
            }
            if end_1based < tract.end_1based {
                record_insert_shift(free.len().saturating_sub(insert_at));
                free.insert(
                    insert_at,
                    FreeTract {
                        sequence_ordinal: tract.sequence_ordinal,
                        start_1based: end_1based + 1,
                        end_1based: tract.end_1based,
                    },
                );
                #[cfg(test)]
                {
                    inserted_residues += 1;
                }
            }
            #[cfg(test)]
            {
                if inserted_residues == 2 {
                    record_division();
                    record_fragment_growth();
                }
            }
            placed.push(PlacedBlockDraft {
                source_block_index,
                length_bp,
                sequence_ordinal: tract.sequence_ordinal,
                start_1based,
                end_1based,
            });
            #[cfg(test)]
            record_live_slots(
                inputs
                    .eligible_tracts
                    .len()
                    .saturating_add(pending_capacity)
                    .saturating_add(free.capacity())
                    .saturating_add(placed.capacity()),
            );
        }

        placed.sort_unstable_by(|left, right| {
            (left.sequence_ordinal, left.start_1based, left.end_1based).cmp(&(
                right.sequence_ordinal,
                right.start_1based,
                right.end_1based,
            ))
        });
        let placed_bases = placed.iter().try_fold(0_u64, |sum, block| {
            sum.checked_add(block.length_bp).ok_or_else(|| {
                Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "la suma de bloques colocados desborda u64",
                )
            })
        })?;
        let remaining_target_bases = planned_bases.checked_sub(placed_bases).ok_or_else(|| {
            Diagnostic::f5_internal(
                scenario_id.as_str(),
                "los bloques colocados exceden la longitud planificada",
            )
        })?;
        let placed_layout_bytes = placed_layout_bytes(&placed);
        let placed_count = u64::try_from(placed.len()).map_err(|_| {
            Diagnostic::f5_internal(
                scenario_id.as_str(),
                "el recuento de bloques colocados no cabe en u64",
            )
        })?;
        consider_best(
            &mut best,
            inputs.length_plan.block_count() - placed_count,
            placed_count,
            placed_bases,
            remaining_target_bases,
            placed_layout_bytes,
            attempt_index,
        );

        if placement_failed {
            #[cfg(test)]
            record_retry();
            continue;
        }
        match materialize_attempt(&inputs, attempt_index, &placed) {
            Ok(layout) => return Ok(layout),
            Err(AttemptMaterialization::TooFewSites) => {
                #[cfg(test)]
                record_retry();
                continue;
            }
            Err(AttemptMaterialization::Diagnostic(error)) => return Err(error.into()),
        }
    }
    let best = best.ok_or_else(|| {
        Diagnostic::f5_internal(
            scenario_id.as_str(),
            "no se registró ningún intento de layout",
        )
    })?;
    Err(F5LayoutFailure::exhausted(
        scenario_id.as_str(),
        best.into_summary(),
    ))
}

enum AttemptMaterialization {
    TooFewSites,
    Diagnostic(Diagnostic),
}

fn validate_layout_request(inputs: &F5LayoutInputs<'_>) -> Result<(), Diagnostic> {
    let scenario_id = inputs.length_plan.scenario_id();
    if inputs.max_layout_attempts == 0 {
        return Err(Diagnostic::f5_internal(
            scenario_id.as_str(),
            "max_layout_attempts cero alcanzó F5.1.2 después de la admisión",
        ));
    }
    if inputs.root_seed != inputs.identity_context.root_seed()
        || inputs.identity_context.scenario_id() != scenario_id
    {
        return Err(Diagnostic::f5_internal(
            scenario_id.as_str(),
            "las entradas de layout no comparten el preflight, la semilla y la identidad F5 admitidos",
        ));
    }
    if inputs.sequences.is_none() && inputs.streamed_reference.is_none() {
        return Err(Diagnostic::f5_internal(
            scenario_id.as_str(),
            "el layout no recibió ni secuencias de prueba ni referencia FASTA de producción",
        ));
    }
    if let Some(sequences) = inputs.sequences {
        for (index, sequence) in sequences.iter().enumerate() {
            if sequences[..index]
                .iter()
                .any(|other| other.sequence_ordinal == sequence.sequence_ordinal)
            {
                return Err(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "las secuencias de referencia repiten un ordinal FASTA",
                ));
            }
        }
    }
    let mut previous = None;
    for tract in inputs.eligible_tracts {
        record_request_tract_visit();
        if tract.start_1based() == 0 || tract.end_1based() < tract.start_1based() {
            return Err(Diagnostic::f5_internal(
                scenario_id.as_str(),
                "un tramo elegible tiene coordenadas cerradas inválidas",
            ));
        }
        if let Some(sequences) = inputs.sequences {
            let sequence = sequences
                .iter()
                .find(|sequence| sequence.sequence_ordinal == tract.sequence_ordinal())
                .ok_or_else(|| {
                    Diagnostic::f5_internal(
                        scenario_id.as_str(),
                        "falta una secuencia para un tramo elegible",
                    )
                })?;
            let start = usize::try_from(tract.start_1based() - 1).map_err(|_| {
                Diagnostic::f5_internal(scenario_id.as_str(), "el inicio de tramo no cabe en usize")
            })?;
            let end = usize::try_from(tract.end_1based()).map_err(|_| {
                Diagnostic::f5_internal(scenario_id.as_str(), "el final de tramo no cabe en usize")
            })?;
            if start >= end || end > sequence.bases.len() {
                return Err(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "un tramo elegible no pertenece a las secuencias de referencia recibidas",
                ));
            }
        }
        let key = (
            tract.sequence_ordinal(),
            tract.start_1based(),
            tract.end_1based(),
        );
        if previous.is_some_and(|last: (u64, u64, u64)| {
            last.0 > key.0 || (last.0 == key.0 && (last.1 >= key.1 || last.2 >= key.1))
        }) {
            return Err(Diagnostic::f5_internal(
                scenario_id.as_str(),
                "los tramos de preflight no conservan el orden FASTA/coordenada",
            ));
        }
        previous = Some(key);
    }
    Ok(())
}

fn resolve_start_by_prefix(
    free: &[FreeTract],
    length_bp: u64,
    chosen: u64,
) -> Option<(usize, u64)> {
    let mut prefix = 0_u64;
    for (index, tract) in free.iter().enumerate() {
        record_prefix_tract_visit();
        if tract.len() < length_bp {
            continue;
        }
        let starts = tract.len() - length_bp + 1;
        if chosen < prefix.checked_add(starts)? {
            return Some((index, tract.start_1based + (chosen - prefix)));
        }
        prefix = prefix.checked_add(starts)?;
    }
    None
}

fn total_starts_for_length(free: &[FreeTract], length_bp: u64) -> Option<u64> {
    let mut total_starts = 0_u64;
    for tract in free {
        record_total_tract_visit();
        if tract.len() >= length_bp {
            total_starts = total_starts.checked_add(tract.len() - length_bp + 1)?;
        }
    }
    Some(total_starts)
}

fn placed_layout_bytes(placed: &[PlacedBlockDraft]) -> Vec<u8> {
    let rows = placed
        .iter()
        .map(|block| {
            (
                block.sequence_ordinal,
                block.start_1based,
                block.end_1based,
                block.length_bp,
            )
        })
        .collect::<Vec<_>>();
    f5_placed_layout_bytes(&rows)
}

struct BestLayoutCandidate {
    unplaced_blocks: u64,
    placed_blocks: u64,
    placed_bases: u64,
    remaining_target_bases: u64,
    bytes: Vec<u8>,
    attempt_index: u64,
}

impl BestLayoutCandidate {
    fn into_summary(self) -> F5BestLayout {
        F5BestLayout {
            attempt_index: self.attempt_index,
            placed_blocks: self.placed_blocks,
            unplaced_blocks: self.unplaced_blocks,
            placed_bases: self.placed_bases,
            remaining_target_bases: self.remaining_target_bases,
            layout_bytes_sha256: format!("{:x}", Sha256::digest(&self.bytes)),
        }
    }
}

fn consider_best(
    best: &mut Option<BestLayoutCandidate>,
    unplaced_blocks: u64,
    placed_blocks: u64,
    placed_bases: u64,
    remaining_target_bases: u64,
    placed_layout_bytes: Vec<u8>,
    attempt_index: u64,
) {
    let candidate = BestLayoutCandidate {
        unplaced_blocks,
        placed_blocks,
        placed_bases,
        remaining_target_bases,
        bytes: placed_layout_bytes,
        attempt_index,
    };
    if best.as_ref().is_none_or(|current| {
        (
            candidate.unplaced_blocks,
            candidate.remaining_target_bases,
            &candidate.bytes,
            candidate.attempt_index,
        ) < (
            current.unplaced_blocks,
            current.remaining_target_bases,
            &current.bytes,
            current.attempt_index,
        )
    }) {
        *best = Some(candidate);
    }
}

fn materialize_attempt(
    request: &F5LayoutInputs<'_>,
    attempt_index: u64,
    placed: &[PlacedBlockDraft],
) -> Result<F5PhysicalLayout, AttemptMaterialization> {
    let scenario_id = request.length_plan.scenario_id();
    let mut positions = Vec::with_capacity(placed.len());
    for block in placed {
        let mut sites = vec![block.start_1based];
        let mut gap_rng = derive_f5_subseed(
            request.root_seed,
            F5RandomDecision::layout_gap(scenario_id, attempt_index, block.source_block_index),
        )
        .map_err(AttemptMaterialization::Diagnostic)?
        .local_rng();
        loop {
            let gap = request
                .gap_distribution
                .sample(scenario_id, &mut gap_rng)
                .map_err(AttemptMaterialization::Diagnostic)?;
            let current = sites.last().copied().ok_or_else(|| {
                AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "el bloque no conserva su primer SNP",
                ))
            })?;
            let remaining = block.end_1based.checked_sub(current).ok_or_else(|| {
                AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "la posición SNP quedó fuera del bloque",
                ))
            })?;
            if gap > remaining {
                break;
            }
            // La comparación anterior demuestra que esta suma cabe y no supera
            // el final del bloque, incluso para un gap u64::MAX válido.
            let next = current.checked_add(gap).ok_or_else(|| {
                AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "la posición SNP desborda pese a la cota restante",
                ))
            })?;
            sites.push(next);
            if sites.len() > 4096 {
                return Err(AttemptMaterialization::Diagnostic(Diagnostic::f5_resource(
                    scenario_id.as_str(),
                    "S excede 4096 durante la inserción de SNP",
                )));
            }
        }
        if sites.len() < 2 {
            return Err(AttemptMaterialization::TooFewSites);
        }
        positions.push((block.sequence_ordinal, sites));
    }

    let requested_sites = positions
        .iter()
        .flat_map(|(ordinal, sites)| sites.iter().map(move |position| (*ordinal, *position)))
        .collect::<Vec<_>>();
    let site_bases = if let Some(sequences) = request.sequences {
        requested_sites
            .iter()
            .map(|(ordinal, position)| {
                let sequence = sequences
                    .iter()
                    .find(|sequence| sequence.sequence_ordinal == *ordinal)
                    .ok_or_else(|| {
                        AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                            scenario_id.as_str(),
                            "un bloque colocado no tiene secuencia de referencia",
                        ))
                    })?;
                let offset = usize::try_from(*position - 1).map_err(|_| {
                    AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                        scenario_id.as_str(),
                        "la coordenada SNP no cabe en usize",
                    ))
                })?;
                let base = *sequence.bases.get(offset).ok_or_else(|| {
                    AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                        scenario_id.as_str(),
                        "el preflight marcó elegible una coordenada fuera de referencia",
                    ))
                })?;
                Ok(F5ResolvedSiteBase {
                    sequence_id: sequence.sequence_id.clone(),
                    base,
                })
            })
            .collect::<Result<Vec<_>, AttemptMaterialization>>()?
    } else {
        let reference = request.streamed_reference.ok_or_else(|| {
            AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                scenario_id.as_str(),
                "el layout de producción no recibió la referencia a recorrer",
            ))
        })?;
        let streamed = stream_fasta_site_bases(
            reference.path,
            reference.reference_id,
            reference.sealed_size_bytes,
            reference.sealed_sha256,
            &requested_sites,
        )
        .map_err(AttemptMaterialization::Diagnostic)?;
        if streamed.reference_bases_processed() != reference.expected_bases_processed {
            return Err(AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                scenario_id.as_str(),
                "la pasada streaming de REF no recorrió las mismas bases que el preflight admitido",
            )));
        }
        streamed
            .sites()
            .iter()
            .zip(&requested_sites)
            .map(|(site, (ordinal, position))| {
                if site.sequence_ordinal() != *ordinal || site.position_1based() != *position {
                    return Err(AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                        scenario_id.as_str(),
                        "la recuperación streaming de REF alteró el orden de los sitios",
                    )));
                }
                Ok(F5ResolvedSiteBase {
                    sequence_id: site.sequence_id().clone(),
                    base: site.base(),
                })
            })
            .collect::<Result<Vec<_>, AttemptMaterialization>>()?
    };
    if site_bases.len() != requested_sites.len() {
        return Err(AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
            scenario_id.as_str(),
            "la recuperación streaming de REF no devolvió todos los sitios",
        )));
    }
    if site_bases
        .iter()
        .any(|site| !matches!(site.base, b'A' | b'C' | b'G' | b'T'))
    {
        return Err(AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
            scenario_id.as_str(),
            "el preflight marcó elegible una base no A/C/G/T",
        )));
    }

    let mut registry = F5IdentityRegistry::default();
    let mut blocks = Vec::with_capacity(placed.len());
    let mut resolved_sites = site_bases.into_iter();
    for (block_index, (draft, (_, site_positions))) in placed.iter().zip(positions).enumerate() {
        let block_index = u64::try_from(block_index).map_err(|_| {
            AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                scenario_id.as_str(),
                "el índice canónico de bloque no cabe en u64",
            ))
        })?;
        let block_id = registry
            .derive_block(request.identity_context, block_index)
            .map_err(AttemptMaterialization::Diagnostic)?;
        let mut sites = Vec::with_capacity(site_positions.len());
        let mut sequence_id = None;
        for (site_index, position_1based) in site_positions.into_iter().enumerate() {
            let site_index = u64::try_from(site_index).map_err(|_| {
                AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "el índice de sitio no cabe en u64",
                ))
            })?;
            let site_id = registry
                .derive_site(request.identity_context, block_index, site_index)
                .map_err(AttemptMaterialization::Diagnostic)?;
            let resolved = resolved_sites.next().ok_or_else(|| {
                AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "faltó una base REF tras la recuperación streaming",
                ))
            })?;
            if let Some(previous) = &sequence_id
                && previous != &resolved.sequence_id
            {
                return Err(AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "un bloque F5 recuperó sitios de más de una secuencia",
                )));
            }
            sequence_id = Some(resolved.sequence_id);
            let reference_base = resolved.base;
            let mut rng = derive_f5_subseed(
                request.root_seed,
                F5RandomDecision::site_alt(scenario_id, &block_id, &site_id),
            )
            .map_err(AttemptMaterialization::Diagnostic)?
            .local_rng();
            let alternatives = (*b"ACGT")
                .into_iter()
                .filter(|base| *base != reference_base)
                .collect::<Vec<_>>();
            let choice = f5_uniform_index(&mut rng, 3).map_err(|_| {
                AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "el selector ALT recibió un límite inválido",
                ))
            })?;
            let alternate_base = alternatives[usize::try_from(choice).map_err(|_| {
                AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "el índice ALT no cabe en usize",
                ))
            })?];
            sites.push(F5LayoutSite {
                site_index,
                site_id,
                position_1based,
                reference_base,
                alternate_base,
            });
        }
        blocks.push(F5LayoutBlock {
            block_index,
            block_id,
            sequence_ordinal: draft.sequence_ordinal,
            sequence_id: sequence_id.ok_or_else(|| {
                AttemptMaterialization::Diagnostic(Diagnostic::f5_internal(
                    scenario_id.as_str(),
                    "un bloque F5 no recuperó ninguna base REF",
                ))
            })?,
            start_1based: draft.start_1based,
            end_1based: draft.end_1based,
            length_bp: draft.length_bp,
            sites,
        });
    }
    Ok(F5PhysicalLayout {
        blocks,
        successful_attempt_index: attempt_index,
    })
}

struct F5ResolvedSiteBase {
    sequence_id: SequenceId,
    base: u8,
}

fn validate_weighted_integer_bands(
    scenario_id: &ScenarioId,
    bands: &[F5BlockLengthBand],
    minimum_allowed: u64,
    name: &str,
) -> Result<(), Diagnostic> {
    let first = bands
        .first()
        .ok_or_else(|| Diagnostic::f5_domain(scenario_id.as_str(), format!("perfil sin {name}")))?;
    if first.minimum < minimum_allowed {
        return Err(Diagnostic::f5_domain(
            scenario_id.as_str(),
            format!("mínimo inválido de {name}"),
        ));
    }
    let mut previous = None;
    let mut sum = 0.0_f64;
    let mut compensation = 0.0_f64;
    for band in bands {
        if band.minimum < minimum_allowed
            || band.minimum > band.maximum
            || !band.weight.is_finite()
            || band.weight <= 0.0
        {
            return Err(Diagnostic::f5_domain(
                scenario_id.as_str(),
                format!("intervalo o peso inválido de {name}"),
            ));
        }
        // DTC-089 exige soporte continuo únicamente para longitudes de
        // bloque. Los gaps solo deben estar ordenados y ser disjuntos: los
        // valores ausentes no son una imposibilidad geométrica.
        if previous.is_some_and(|last: u64| band.minimum <= last) {
            return Err(Diagnostic::f5_domain(
                scenario_id.as_str(),
                format!("intervalos solapados o desordenados de {name}"),
            ));
        }
        let next = sum + band.weight;
        compensation += if sum.abs() >= band.weight.abs() {
            (sum - next) + band.weight
        } else {
            (band.weight - next) + sum
        };
        sum = next;
        previous = Some(band.maximum);
    }
    if !(sum + compensation).is_finite() || ((sum + compensation) - 1.0).abs() > 1e-9 {
        return Err(Diagnostic::f5_domain(
            scenario_id.as_str(),
            format!("los pesos de {name} no suman uno"),
        ));
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum WeightedIntegerError {
    Domain(&'static str),
    Internal(&'static str),
}

/// Elección categórica literal de DTC-089.
///
/// Cada entero factible aporta por separado `row_weight / row_width`, en orden
/// creciente. La suma total y cada prefijo usan el mismo acumulador de
/// Neumaier. Una única alternativa no abre el consumidor U53.
fn choose_weighted_integer_with_source(
    bands: &[F5BlockLengthBand],
    lower: u64,
    upper: u64,
    rng: &mut impl RandomWordSource,
) -> Result<u64, WeightedIntegerError> {
    let mut total = 0.0_f64;
    let mut compensation = 0.0_f64;
    let mut alternatives = 0_u64;
    let mut only = None;
    for band in bands {
        let start = band.minimum.max(lower);
        let end = band.maximum.min(upper);
        if start > end {
            continue;
        }
        let width = band
            .maximum
            .checked_sub(band.minimum)
            .and_then(|value| value.checked_add(1))
            .ok_or(WeightedIntegerError::Internal(
                "el ancho del soporte desborda u64",
            ))?;
        let per_value = band.weight / width as f64;
        if !per_value.is_finite() || per_value <= 0.0 {
            return Err(WeightedIntegerError::Internal(
                "el peso por longitud no es finito y positivo",
            ));
        }
        for value in start..=end {
            neumaier_add(&mut total, &mut compensation, per_value)?;
            alternatives = alternatives
                .checked_add(1)
                .ok_or(WeightedIntegerError::Internal(
                    "la cardinalidad factible desborda u64",
                ))?;
            only = Some(value);
        }
    }
    let normalized_total = corrected(total, compensation)?;
    if alternatives == 0 {
        return Err(WeightedIntegerError::Domain(
            "soporte restringido sin valores",
        ));
    }
    if alternatives == 1 {
        return only.ok_or(WeightedIntegerError::Domain(
            "soporte restringido sin valor único",
        ));
    }

    let bits53 = rng.next_word() >> 11;
    let unit = bits53 as f64 / (1_u64 << 53) as f64;
    let threshold = unit * normalized_total;
    if !threshold.is_finite() {
        return Err(WeightedIntegerError::Internal(
            "el umbral categórico no es finito",
        ));
    }
    let mut cumulative = 0.0_f64;
    let mut cumulative_compensation = 0.0_f64;
    let mut last = None;
    for band in bands {
        let start = band.minimum.max(lower);
        let end = band.maximum.min(upper);
        if start > end {
            continue;
        }
        let width = band
            .maximum
            .checked_sub(band.minimum)
            .and_then(|value| value.checked_add(1))
            .ok_or(WeightedIntegerError::Internal(
                "el ancho del soporte desborda u64",
            ))?;
        let per_value = band.weight / width as f64;
        if !per_value.is_finite() || per_value <= 0.0 {
            return Err(WeightedIntegerError::Internal(
                "el peso por longitud no es finito y positivo",
            ));
        }
        for value in start..=end {
            let corrected_next =
                neumaier_add(&mut cumulative, &mut cumulative_compensation, per_value)?;
            if threshold < corrected_next {
                return Ok(value);
            }
            last = Some(value);
        }
    }
    // El último entero absorbe exclusivamente el residuo de redondeo permitido.
    last.ok_or(WeightedIntegerError::Domain(
        "soporte restringido sin valores",
    ))
}

fn neumaier_add(
    sum: &mut f64,
    compensation: &mut f64,
    value: f64,
) -> Result<f64, WeightedIntegerError> {
    let next = *sum + value;
    if !next.is_finite() {
        return Err(WeightedIntegerError::Internal(
            "la suma de pesos no es finita",
        ));
    }
    let correction = if sum.abs() >= value.abs() {
        (*sum - next) + value
    } else {
        (value - next) + *sum
    };
    let next_compensation = *compensation + correction;
    if !next_compensation.is_finite() {
        return Err(WeightedIntegerError::Internal(
            "la compensación de pesos no es finita",
        ));
    }
    *sum = next;
    *compensation = next_compensation;
    corrected(next, next_compensation)
}

fn corrected(sum: f64, compensation: f64) -> Result<f64, WeightedIntegerError> {
    let value = sum + compensation;
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(WeightedIntegerError::Internal(
            "la suma compensada de pesos no es finita y positiva",
        ))
    }
}

fn ceil_div(value: u64, divisor: u64) -> Option<u64> {
    (divisor != 0).then(|| value / divisor + u64::from(!value.is_multiple_of(divisor)))
}

fn domain(scenario_id: &ScenarioId, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::f5_domain(scenario_id.as_str(), detail)
}

fn internal(scenario_id: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::f5_internal(scenario_id, detail)
}

fn configuration(
    config: &EffectiveF5Config,
    fields: impl Into<String>,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::f5_configuration(config.scenario_path().display().to_string(), fields, detail)
}

fn allocation_resource(block_count: u64) -> Diagnostic {
    let required = block_count
        .checked_mul(u64::try_from(std::mem::size_of::<u64>()).expect("u64 size fits u64"))
        .map_or_else(|| ">u64::MAX".to_owned(), |bytes| bytes.to_string());
    Diagnostic::resource(
        "GPS-RES-002",
        "Los recursos disponibles no permiten completar el plan sin degradar resultados obligatorios.",
    )
    .with("resource", "f5_block_length_vector")
    .with("required", required)
    .with("available", usize::MAX.to_string())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::{
        F5BlockLengthBand, F5BlockLengthDistribution, F5BlockLengthPlan, F5LayoutRequest,
        F5LayoutSequence, F5LengthFeasibility, F5SnpGapDistribution, build_f5_physical_layout,
    };
    use crate::application::load_and_validate_f5_config;
    use crate::f5_preflight::{F5EligibleTract, F5OtherResources, recover_sealed_f5_admission};
    use crate::identity::{ExecutionId, F5IdentityContext, ReferenceId, ScenarioId, SequenceId};
    use crate::randomness::{RandomWordSource, RootSeed};
    use sha2::{Digest, Sha256};

    static PRODUCTION_FIXTURE_NONCE: AtomicU64 = AtomicU64::new(0);

    fn sha256(bytes: &[u8]) -> String {
        format!("{:x}", Sha256::digest(bytes))
    }

    fn production_fixture_root() -> PathBuf {
        std::env::temp_dir().join(format!(
            "gps-f5-layout-production-{}-{}",
            std::process::id(),
            PRODUCTION_FIXTURE_NONCE.fetch_add(1, Ordering::Relaxed)
        ))
    }

    /// Escribe el documento base con la línea del límite **siempre presente**.
    ///
    /// Los bytes del escenario entran dos veces en la cota de operaciones, de
    /// modo que añadir o quitar la línea `max_operation_units` desplaza la
    /// propia cota. Sólo un literal de anchura constante mantiene el punto
    /// fijo `operation_units == max_operation_units` que este fixture exige,
    /// como ya hace `f5_bootstrap_dispatch::fixed_width_bootstrap_memory_…`
    /// para memoria.
    fn write_production_layout_fixture(
        root: &Path,
        operation_limit: u64,
        seed_byte: u8,
    ) -> PathBuf {
        let scenario_directory = root.join("authoring/scenarios");
        let profile_directory = scenario_directory.join("profiles");
        fs::create_dir_all(&profile_directory).expect("production profile directory");
        let files = [
            (
                "block_lengths.tsv",
                "f5_block_length_profile",
                "length_min_bp\tlength_max_bp\tweight\n3\t3\t1\n",
            ),
            (
                "snp_gaps.tsv",
                "f5_snp_gap_profile",
                "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n",
            ),
            (
                "snp_spectrum.tsv",
                "f5_snp_spectrum",
                "alt_af_min\talt_af_max\tweight\n0\t1\t1\n",
            ),
            (
                "haplotype_weights.tsv",
                "f5_haplotype_weights",
                "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
            ),
            (
                "ld_profile.tsv",
                "f5_ld_profile",
                "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t2\t0.5\t0.05\n",
            ),
        ];
        let links = files
            .iter()
            .map(|(name, kind, contents)| {
                fs::write(profile_directory.join(name), contents).expect("production profile file");
                format!(
                    r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#,
                    contents.len(),
                    sha256(contents.as_bytes())
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let metadata = format!(
            r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
        );
        fs::write(profile_directory.join("profile.metadata.json"), &metadata)
            .expect("production profile metadata");
        let fasta = b">contig-a\nAAAAAANNNNAAANNNNAAA\n";
        let reference_directory = root.join("authoring/references");
        fs::create_dir_all(&reference_directory).expect("production reference directory");
        fs::write(reference_directory.join("ref.fa"), fasta).expect("production FASTA");
        let general = root.join("authoring/general.toml");
        fs::write(
            &general,
            format!(
                "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{}\"\nploidy = 2\nreference_id = \"reference-main\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
                format!("{seed_byte:02x}").repeat(32),
                sha256(fasta)
            ),
        )
        .expect("production general");
        let operation_limit = format!("max_operation_units = {operation_limit}\n");
        fs::write(
            scenario_directory.join("scenario.toml"),
            format!(
                "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = 1\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 2\nmax_haplotype_attempts = 1\nmax_swap_proposals = 1\nmax_haplotype_cells = 16\n{operation_limit}block_count = 4\n\n[[subpopulations]]\npopulation_id = \"P1\"\nindividual_count = 1\nfounder_ids = [\"I1\"]\n",
                sha256(metadata.as_bytes()),
            ),
        )
        .expect("production scenario");
        general
    }

    #[derive(Debug)]
    struct ScriptedWords {
        words: Vec<u64>,
        consumed: usize,
    }

    impl ScriptedWords {
        fn new(words: &[u64]) -> Self {
            Self {
                words: words.to_vec(),
                consumed: 0,
            }
        }
    }

    impl RandomWordSource for ScriptedWords {
        fn next_word(&mut self) -> u64 {
            let word = self.words[self.consumed];
            self.consumed += 1;
            word
        }
    }

    fn scenario() -> ScenarioId {
        ScenarioId::parse("scenario", "test").expect("scenario")
    }

    fn distribution() -> F5BlockLengthDistribution {
        F5BlockLengthDistribution::new(
            &scenario(),
            vec![
                F5BlockLengthBand::new(2, 3, 0.5),
                F5BlockLengthBand::new(4, 5, 0.5),
            ],
        )
        .expect("distribution")
    }

    fn fixed_plan(block_count: u64, length: u64) -> F5BlockLengthPlan {
        F5BlockLengthPlan {
            scenario_id: scenario(),
            distribution: F5BlockLengthDistribution {
                bands: vec![F5BlockLengthBand::new(length, length, 1.0)],
                minimum: length,
                maximum: length,
            },
            block_count,
            target_bases: None,
        }
    }

    fn gap_distribution(gap: u64) -> F5SnpGapDistribution {
        F5SnpGapDistribution::new(&scenario(), vec![F5BlockLengthBand::new(gap, gap, 1.0)])
            .expect("gap distribution")
    }

    fn sequence(ordinal: u64, id: &str, bases: &[u8]) -> F5LayoutSequence {
        F5LayoutSequence::new(
            ordinal,
            SequenceId::parse(id, "test").expect("sequence id"),
            bases.to_vec(),
        )
    }

    fn root_for_layout() -> RootSeed {
        RootSeed::from_bytes([0x31; 32])
    }

    #[test]
    fn conditioned_lengths_have_exact_sum_and_final_residual() {
        let scenario = scenario();
        let plan = F5BlockLengthPlan {
            scenario_id: scenario.clone(),
            distribution: distribution(),
            block_count: 2,
            target_bases: Some(8),
        };
        let lengths = plan
            .sample(RootSeed::from_bytes([7; 32]), 0)
            .expect("lengths");
        assert_eq!(lengths.len(), 2);
        assert_eq!(lengths.iter().sum::<u64>(), 8);
        assert!((2..=5).contains(lengths.last().expect("residual")));
        assert_eq!(lengths[1], 8 - lengths[0]);
    }

    #[test]
    fn conditioned_sampling_restricts_early_values_to_a_feasible_residual() {
        let scenario = scenario();
        let plan = F5BlockLengthPlan {
            scenario_id: scenario.clone(),
            distribution: distribution(),
            block_count: 2,
            target_bases: Some(6),
        };
        for attempt in 0..32 {
            let lengths = plan
                .sample(RootSeed::from_bytes([attempt as u8; 32]), attempt)
                .expect("lengths");
            assert_eq!(lengths.iter().sum::<u64>(), 6);
            assert_eq!(lengths[1], 6 - lengths[0]);
            assert!((2..=4).contains(&lengths[0]));
        }
    }

    #[test]
    fn block_count_mode_is_stochastic_and_has_no_target_sum() {
        let scenario = scenario();
        let plan = F5BlockLengthPlan {
            scenario_id: scenario.clone(),
            distribution: distribution(),
            block_count: 4,
            target_bases: None,
        };
        let first = plan
            .sample(RootSeed::from_bytes([1; 32]), 0)
            .expect("first");
        let second = plan
            .sample(RootSeed::from_bytes([2; 32]), 0)
            .expect("second");
        assert!(first.iter().all(|value| (2..=5).contains(value)));
        assert!(second.iter().all(|value| (2..=5).contains(value)));
        assert_ne!(first, second);
    }

    #[test]
    fn literal_integer_categorical_has_frozen_u53_and_neumaier_vector() {
        let scenario = scenario();
        let normative_plan = F5BlockLengthPlan {
            scenario_id: scenario.clone(),
            distribution: F5BlockLengthDistribution {
                bands: vec![
                    F5BlockLengthBand::new(2, 4, 0.1),
                    F5BlockLengthBand::new(5, 7, 0.9),
                ],
                minimum: 2,
                maximum: 7,
            },
            block_count: 1,
            target_bases: None,
        };
        let bits53 = 0x16_6666_6666_6666_u64;
        let mut normative_u53 = ScriptedWords::new(&[bits53 << 11]);
        assert_eq!(
            normative_plan.sample_between_with_source(2, 7, &mut normative_u53),
            Ok(7),
            "el vector normativo acumula cada entero creciente, no cada banda"
        );
        assert_eq!(normative_u53.consumed, 1);

        let plan = F5BlockLengthPlan {
            scenario_id: scenario.clone(),
            distribution: F5BlockLengthDistribution {
                bands: vec![
                    F5BlockLengthBand::new(2, 2, 1.0),
                    F5BlockLengthBand::new(3, 3, 1e-16),
                    F5BlockLengthBand::new(4, 4, 1e-16),
                ],
                minimum: 2,
                maximum: 4,
            },
            block_count: 1,
            target_bases: None,
        };
        let mut maximum_u53 = ScriptedWords::new(&[u64::MAX]);
        assert_eq!(
            plan.sample_between_with_source(2, 4, &mut maximum_u53),
            Ok(4),
            "la suma Neumaier literal 1+1e-16+1e-16 cruza el U53 máximo"
        );
        assert_eq!(maximum_u53.consumed, 1);

        let mut zero_u53 = ScriptedWords::new(&[0]);
        assert_eq!(plan.sample_between_with_source(2, 4, &mut zero_u53), Ok(2));
        assert_eq!(zero_u53.consumed, 1);

        let mut no_words = ScriptedWords::new(&[]);
        assert_eq!(plan.sample_between_with_source(3, 3, &mut no_words), Ok(3));
        assert_eq!(no_words.consumed, 0, "una alternativa no consume U53");
    }

    #[test]
    fn post_preflight_numeric_failure_is_typed_f5_internal_with_scenario() {
        let scenario = scenario();
        let non_finite_plan = F5BlockLengthPlan {
            scenario_id: scenario.clone(),
            distribution: F5BlockLengthDistribution {
                bands: vec![
                    F5BlockLengthBand::new(2, 2, f64::MAX),
                    F5BlockLengthBand::new(3, 3, f64::MAX),
                ],
                minimum: 2,
                maximum: 3,
            },
            block_count: 1,
            target_bases: None,
        };
        let error = non_finite_plan
            .sample_between_with_source(2, 3, &mut ScriptedWords::new(&[]))
            .expect_err("la suma no finita posterior al preflight debe fallar");
        assert_eq!(error.code(), "GPS-INTERNAL-101");
        assert_eq!(error.context()["scenario_id"], "scenario");
        assert!(error.context()["detail"].contains("no es finita"));

        let overflow_plan = F5BlockLengthPlan {
            scenario_id: scenario.clone(),
            distribution: F5BlockLengthDistribution {
                bands: vec![F5BlockLengthBand::new(2, u64::MAX, 1.0)],
                minimum: 2,
                maximum: u64::MAX,
            },
            block_count: 3,
            target_bases: Some(u64::MAX),
        };
        let error = overflow_plan
            .sample(RootSeed::from_bytes([0; 32]), 0)
            .expect_err("el overflow posterior al preflight debe fallar");
        assert_eq!(error.code(), "GPS-INTERNAL-101");
        assert_eq!(error.context()["scenario_id"], "scenario");
        assert!(error.context()["detail"].contains("desborda"));
    }

    #[test]
    fn failed_length_vector_reservation_is_a_resource_failure() {
        let plan = F5BlockLengthPlan {
            scenario_id: scenario(),
            distribution: F5BlockLengthDistribution {
                bands: vec![F5BlockLengthBand::new(2, 2, 1.0)],
                minimum: 2,
                maximum: 2,
            },
            block_count: u64::MAX,
            target_bases: None,
        };
        let error = plan
            .sample(RootSeed::from_bytes([0; 32]), 0)
            .expect_err("the allocator must reject an unrepresentable vector");
        assert_eq!(error.code(), "GPS-RES-002");
        assert_eq!(error.context()["resource"], "f5_block_length_vector");
    }

    #[test]
    fn feasibility_oracle_uses_exact_integer_bounds_and_rejects_overflow() {
        let scenario = scenario();
        let distribution = distribution();
        let result =
            F5LengthFeasibility::for_target(&scenario, 10, &distribution).expect("feasible target");
        assert_eq!(result.block_count(), 2);
        assert_eq!(result.minimum_total(), 4);
        assert_eq!(result.maximum_total(), 10);
        assert_eq!(
            F5LengthFeasibility::with_block_count(&scenario, u64::MAX, 2, &distribution)
                .expect_err("unrepresentable target")
                .code(),
            "GPS-DOMAIN-101"
        );
        let huge = F5BlockLengthDistribution::new(
            &scenario,
            vec![F5BlockLengthBand::new(u64::MAX - 1, u64::MAX, 1.0)],
        )
        .expect("huge support");
        assert_eq!(
            F5LengthFeasibility::with_block_count(&scenario, u64::MAX, 2, &huge)
                .expect_err("overflow")
                .code(),
            "GPS-DOMAIN-101"
        );
    }

    #[test]
    fn distribution_rejects_non_probability_or_discontinuous_weights() {
        let scenario = scenario();
        for bands in [
            vec![F5BlockLengthBand::new(2, 2, 0.0)],
            vec![F5BlockLengthBand::new(2, 2, f64::NAN)],
            vec![
                F5BlockLengthBand::new(2, 2, 0.5),
                F5BlockLengthBand::new(4, 4, 0.5),
            ],
            vec![F5BlockLengthBand::new(2, 2, 0.8)],
        ] {
            assert_eq!(
                F5BlockLengthDistribution::new(&scenario, bands)
                    .expect_err("restricted weights")
                    .code(),
                "GPS-DOMAIN-101"
            );
        }
    }

    #[test]
    fn gap_distribution_accepts_disjoint_support_without_post_preflight_domain_error() {
        let gaps = F5SnpGapDistribution::new(
            &scenario(),
            vec![
                F5BlockLengthBand::new(1, 1, 0.5),
                F5BlockLengthBand::new(3, 3, 0.5),
            ],
        )
        .expect("DTC-089 permits disjoint SNP gap intervals");
        let mut first = ScriptedWords::new(&[0]);
        let mut second = ScriptedWords::new(&[u64::MAX]);
        assert_eq!(gaps.sample(&scenario(), &mut first), Ok(1));
        assert_eq!(gaps.sample(&scenario(), &mut second), Ok(3));
    }

    fn run_layout(
        plan: &F5BlockLengthPlan,
        gaps: &F5SnpGapDistribution,
        sequences: &[F5LayoutSequence],
        tracts: &[F5EligibleTract],
        max_attempts: u64,
        root: RootSeed,
    ) -> Result<super::F5PhysicalLayout, super::F5LayoutFailure> {
        let execution =
            ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
        let reference = ReferenceId::parse("reference-main", "test").expect("reference");
        let context = F5IdentityContext::new(root, &execution, plan.scenario_id(), &reference);
        build_f5_physical_layout(F5LayoutRequest::from_test_parts(super::F5LayoutInputs {
            length_plan: plan,
            eligible_tracts: tracts,
            gap_distribution: gaps,
            sequences: Some(sequences),
            streamed_reference: None,
            max_layout_attempts: max_attempts,
            root_seed: root,
            identity_context: context,
        }))
    }

    #[test]
    fn physical_layout_uses_preflight_tracts_and_emits_canonical_blocks_sites_and_ids() {
        let plan = fixed_plan(2, 2);
        let gaps = gap_distribution(1);
        let sequences = [
            sequence(0, "contig-a", b"ACGT"),
            sequence(1, "contig-b", b"AANCGT"),
        ];
        // The N at contig-b:3 and an excluded base at :2 are absent from the
        // preflight tracts; this consumer does not scan or recreate them.
        let tracts = [F5EligibleTract::new(0, 1, 2), F5EligibleTract::new(1, 5, 6)];
        let layout = run_layout(&plan, &gaps, &sequences, &tracts, 1, root_for_layout())
            .expect("two complete blocks");

        let blocks = layout.blocks();
        assert_eq!(layout.successful_attempt_index(), 0);
        assert_eq!(blocks.len(), 2);
        assert_eq!(
            blocks
                .iter()
                .map(|block| (
                    block.sequence_id().as_str(),
                    block.start_1based(),
                    block.end_1based()
                ))
                .collect::<Vec<_>>(),
            vec![("contig-a", 1, 2), ("contig-b", 5, 6)]
        );
        assert_eq!(
            blocks[0]
                .sites()
                .iter()
                .map(|site| (
                    site.site_index(),
                    site.position_1based(),
                    site.reference_base(),
                    site.alternate_base()
                ))
                .collect::<Vec<_>>(),
            vec![(0, 1, b'A', b'T'), (1, 2, b'C', b'T')]
        );
        assert_eq!(
            blocks[0].block_id().as_str(),
            "BLK_GQGZK2AJLMHH44IPFRHNGXYV54"
        );
        assert_eq!(
            blocks[0].sites()[0].site_id().as_str(),
            "SITE_OMJCJYTYSN3WRJJOPS2TQPRXDE"
        );
        for block in blocks {
            for site in block.sites() {
                assert_ne!(site.reference_base(), site.alternate_base());
                assert!(matches!(site.alternate_base(), b'A' | b'C' | b'G' | b'T'));
            }
        }
    }

    #[test]
    fn adjacent_and_fragmented_tracts_never_overlap_and_preserve_attempt_invariance() {
        let plan = fixed_plan(2, 2);
        let gaps = gap_distribution(1);
        let sequences = [sequence(0, "contig-a", b"ACGTNNAC")];
        let tracts = [F5EligibleTract::new(0, 1, 2), F5EligibleTract::new(0, 3, 4)];
        let one = run_layout(&plan, &gaps, &sequences, &tracts, 1, root_for_layout())
            .expect("adjacent complete tracts");
        let many = run_layout(&plan, &gaps, &sequences, &tracts, 7, root_for_layout())
            .expect("the same first attempt must remain unchanged");
        assert_eq!(one, many);
        assert_eq!(
            one.blocks()
                .iter()
                .map(|block| (block.start_1based(), block.end_1based()))
                .collect::<Vec<_>>(),
            vec![(1, 2), (3, 4)]
        );
    }

    #[test]
    fn fewer_than_two_sites_discards_the_complete_layout_and_reports_exhaustion() {
        let plan = fixed_plan(1, 2);
        let gaps = gap_distribution(2);
        let sequences = [sequence(0, "contig-a", b"AC")];
        let tracts = [F5EligibleTract::new(0, 1, 2)];
        let error = run_layout(&plan, &gaps, &sequences, &tracts, 2, root_for_layout())
            .expect_err("one SNP is not an F5 block");
        assert_eq!(error.code(), "GPS-DOMAIN-106");
        assert!(error.context()["detail"].contains("unplaced_blocks=0"));
        assert!(error.context()["detail"].contains("remaining_target_bases=0"));
    }

    #[test]
    fn exhausted_layout_keeps_the_normative_best_layout_value_without_text_parsing() {
        let plan = fixed_plan(2, 3);
        let gaps = gap_distribution(1);
        let sequences = [sequence(0, "contig-a", b"ACGTA")];
        let tracts = [F5EligibleTract::new(0, 1, 5)];
        let error = run_layout(&plan, &gaps, &sequences, &tracts, 3, root_for_layout())
            .expect_err("two length-three blocks cannot occupy five bases");
        assert_eq!(error.code(), "GPS-DOMAIN-106");
        let best = error.best_layout().expect("typed best layout");
        assert_eq!(best.attempt_index(), 0);
        assert_eq!(best.placed_blocks(), 1);
        assert_eq!(best.unplaced_blocks(), 1);
        assert_eq!(best.placed_bases(), 3);
        assert_eq!(best.remaining_target_bases(), 3);
        assert_eq!(
            best.layout_bytes_sha256(),
            "d2795e187f9944079317b40cdcff7d7f410b55918213587eff1acd2c6b1509bc"
        );
    }

    #[test]
    fn retries_are_transactional_and_can_succeed_after_a_failed_attempt() {
        let plan = fixed_plan(2, 3);
        let gaps = gap_distribution(1);
        let sequences = [sequence(0, "contig-a", b"ACGTAC")];
        let tracts = [F5EligibleTract::new(0, 1, 6)];
        let root = RootSeed::from_bytes([0x05; 32]);
        assert!(
            run_layout(&plan, &gaps, &sequences, &tracts, 1, root).is_err(),
            "the documented seed exhausts the first transactional attempt"
        );
        super::reset_tract_query_count();
        let layout = run_layout(&plan, &gaps, &sequences, &tracts, 2, root)
            .expect("the second transaction succeeds");
        assert_eq!(layout.successful_attempt_index(), 1);
        assert_eq!(layout.blocks().len(), 2);
        assert_eq!(
            layout
                .blocks()
                .iter()
                .map(|block| block.length_bp())
                .sum::<u64>(),
            6
        );
        assert!(
            super::tract_query_count() <= 2 * 2 * (2 + 1),
            "dos intentos, dos bloques y un tramo inicial respetan la cota v7 de recorridos, incluidos los residuos"
        );
        assert!(
            super::layout_vector_moves() <= 2 * (1 + 3 * 2 + 3),
            "las copias y desplazamientos de remove/insert respetan la cota v7"
        );
        assert_eq!(
            super::layout_reserved_slots(),
            2 * ((1 + 2) + 2),
            "cada intento reserva los tramos libres (tramo inicial+B) y borradores B"
        );
    }

    #[test]
    fn fragmented_full_layout_accounts_for_interior_last_residues_and_retry() {
        let plan = fixed_plan(4, 3);
        let gaps = gap_distribution(1);
        let sequences = [sequence(0, "contig-a", b"AAAAAANNNNAAANNNNAAA")];
        let tracts = [
            F5EligibleTract::new(0, 1, 6),
            F5EligibleTract::new(0, 11, 13),
            F5EligibleTract::new(0, 18, 20),
        ];
        let mut witnessed = None;
        for byte in 0_u8..=u8::MAX {
            super::reset_tract_query_count();
            let root = RootSeed::from_bytes([byte; 32]);
            if let Ok(layout) = run_layout(&plan, &gaps, &sequences, &tracts, 2, root) {
                let trace = super::layout_placement_trace();
                if layout.successful_attempt_index() == 1
                    && trace
                        .iter()
                        .any(|(attempt, _, _, left, right)| *attempt == 0 && *left && *right)
                    && trace
                        .iter()
                        .any(|(attempt, index, _, _, _)| *attempt == 0 && *index >= 2)
                {
                    witnessed = Some((
                        byte,
                        trace,
                        super::tract_query_count(),
                        super::layout_vector_moves(),
                        super::layout_reserved_slots(),
                    ));
                    break;
                }
            }
        }
        let (byte, trace, traversals, moves, reserved) = witnessed
            .expect("a deterministic subsemilla must witness interior split, last tract and retry");
        assert_eq!(
            byte, 2,
            "the witness seed is part of the regression contract"
        );
        assert!(
            trace
                .iter()
                .any(|(attempt, _, _, left, right)| *attempt == 0 && *left && *right)
        );
        assert!(
            trace
                .iter()
                .any(|(attempt, index, _, _, _)| *attempt == 0 && *index >= 2)
        );
        let attempts = 2_u64;
        let blocks = 4_u64;
        let tract_count = 3_u64;
        assert!(
            traversals <= 2 * attempts * (blocks * tract_count + blocks * (blocks - 1) / 2),
            "v7 accounts for every total/prefix traversal after fragment growth"
        );
        assert!(
            moves
                <= attempts
                    * (tract_count + 3 * blocks * tract_count + 3 * blocks * (blocks - 1) / 2),
            "v7 accounts for copies and remove/insert shifts"
        );
        assert_eq!(
            reserved,
            attempts * ((tract_count + blocks) + blocks),
            "the test-only representation reserves free tracts and placed drafts per attempt"
        );
    }

    #[test]
    fn prefix_resolution_is_uniform_over_all_tract_starts_not_over_tracts() {
        let free = [
            super::FreeTract {
                sequence_ordinal: 0,
                start_1based: 10,
                end_1based: 12,
            },
            super::FreeTract {
                sequence_ordinal: 1,
                start_1based: 20,
                end_1based: 21,
            },
        ];
        assert_eq!(super::resolve_start_by_prefix(&free, 2, 0), Some((0, 10)));
        assert_eq!(super::resolve_start_by_prefix(&free, 2, 1), Some((0, 11)));
        assert_eq!(super::resolve_start_by_prefix(&free, 2, 2), Some((1, 20)));
    }

    #[test]
    fn last_tract_winner_uses_the_two_preflighted_tract_queries() {
        let free = [
            super::FreeTract {
                sequence_ordinal: 0,
                start_1based: 10,
                end_1based: 12,
            },
            super::FreeTract {
                sequence_ordinal: 1,
                start_1based: 20,
                end_1based: 21,
            },
        ];
        super::reset_tract_query_count();
        let total = super::total_starts_for_length(&free, 2)
            .expect("the two short tracts fit in the u64 start count");
        assert_eq!(total, 3);
        assert_eq!(
            super::resolve_start_by_prefix(&free, 2, total - 1),
            Some((1, 20)),
            "the global integer selects the only start in the last tract"
        );
        let executed = super::tract_query_count();
        let layout_attempts = 1_u64;
        let block_count = 1_u64;
        let preflighted = 2
            * layout_attempts
            * block_count
            * u64::try_from(free.len()).expect("tract count fits");
        assert_eq!(executed, preflighted);
        assert!(
            executed <= preflighted,
            "the last-tract winner cannot exceed the v7 two-pass tract budget"
        );
    }

    #[test]
    fn production_v7_instrumentation_covers_fragmentation_interior_last_residues_and_retry() {
        let root = production_fixture_root();
        // Semilla testigo del contrato de regresión: es la misma que
        // `fragmented_full_layout_accounts_for_interior_last_residues_and_retry`
        // fija para el mismo escenario fragmentado, y la primera que fuerza a
        // la vez ganador interior con dos residuos, último tramo y reintento.
        let seed_byte = 2;
        // Sonda de anchura fija: cuatro dígitos holgados que el documento
        // admite y que ya ocupan exactamente los mismos bytes que la cota
        // buscada. Sin este marcador la cota se mediría sobre un escenario
        // más corto y el punto fijo sería insatisfacible.
        let probe_limit = 9_999_u64;
        let probe_general = write_production_layout_fixture(&root, probe_limit, seed_byte);
        let baseline_config = load_and_validate_f5_config(&probe_general)
            .expect("the production fixture configuration is valid");
        let baseline = recover_sealed_f5_admission(&baseline_config, F5OtherResources::default())
            .expect("the fixed-width probe admits the production fixture");
        let exact_limit = baseline.operation_units();
        assert_eq!(
            exact_limit.to_string().len(),
            probe_limit.to_string().len(),
            "el punto fijo exige que sonda y cota ocupen los mismos bytes"
        );
        assert_eq!(
            (exact_limit - 1).to_string().len(),
            probe_limit.to_string().len(),
            "la frontera exacta−1 también debe conservar la anchura"
        );

        let exact_general = write_production_layout_fixture(&root, exact_limit, seed_byte);
        // El escaneo de capacidad ocurre dentro de la carga, que es la que
        // abre la referencia; `recover_sealed_f5_admission` ya no reabre nada.
        crate::f5_preflight::reset_test_capacity_tract_visits();
        super::reset_layout_instrumentation();
        let config = load_and_validate_f5_config(&exact_general)
            .expect("the exact operation limit retains the production fixture");
        let preflight = recover_sealed_f5_admission(&config, F5OtherResources::default())
            .expect("v7 admits operation_limit exactly equal to its bound");
        assert_eq!(
            preflight.operation_units(),
            exact_limit,
            "la cota declarada es exactamente la que el documento consume"
        );
        assert_eq!(preflight.operation_limit(), exact_limit);
        let execution = ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test")
            .expect("execution identity");
        let layout = build_f5_physical_layout(
            F5LayoutRequest::new(&config, &preflight, &execution)
                .expect("production preflight binding"),
        )
        .expect("the admitted production layout completes");

        let observed = super::layout_instrumentation();
        let capacity_visits = crate::f5_preflight::test_capacity_tract_visits();
        let attempts = 2_u64;
        let blocks = 4_u64;
        let tracts = 3_u64;
        let fragment_growth = blocks * (blocks - 1) / 2;
        let per_attempt_tract_visits = blocks * tracts + fragment_growth;
        let per_attempt_shift_bound = blocks * tracts + fragment_growth;
        let copy_and_shift_bound = attempts * (tracts + 3 * per_attempt_shift_bound);
        let slot_bound = 2 * tracts + 4 * blocks;

        assert_eq!(
            layout.successful_attempt_index(),
            1,
            "fixture exercises a retry; observed={observed:#?}"
        );
        assert_eq!(
            capacity_visits, tracts,
            "v7 capacity scan visits every preflight tract once"
        );
        assert_eq!(
            observed.request_tract_visits, tracts,
            "v7 request validation visits every preflight tract once"
        );
        assert_eq!(
            observed.initial_tract_copies,
            attempts * tracts,
            "one initial free-tract copy occurs for every tract and attempt"
        );
        assert_eq!(
            observed.reservation_calls,
            2 * attempts,
            "each attempt reserves free and placed vectors"
        );
        assert_eq!(
            observed.reserved_slots,
            attempts * (tracts + 2 * blocks),
            "the observed reservation slots match the two production vectors per attempt"
        );
        assert!(
            observed.total_tract_visits <= attempts * per_attempt_tract_visits,
            "v7 total-start scans include at most T+B(B-1)/2 live tracts per block"
        );
        assert!(
            observed.prefix_tract_visits <= attempts * per_attempt_tract_visits,
            "v7 prefix scans have the same independent bound"
        );
        assert!(
            observed.fragment_growth <= attempts * blocks,
            "a placement can grow the free list by at most one tract"
        );
        assert!(
            observed.divisions <= attempts * blocks,
            "at most one two-residue division occurs per placed block"
        );
        assert!(
            observed.remove_shifts <= attempts * per_attempt_shift_bound,
            "remove shifts fit their v7 live-list bound"
        );
        assert!(
            observed.insert_shifts <= 2 * attempts * per_attempt_shift_bound,
            "the two possible inserts fit the remaining v7 shift bound"
        );
        assert!(
            observed.initial_tract_copies + observed.remove_shifts + observed.insert_shifts
                <= copy_and_shift_bound,
            "v7 separately budgets copies plus remove/insert shifts"
        );
        assert_eq!(observed.attempts_started, attempts);
        assert_eq!(
            observed.retries,
            attempts - 1,
            "the first attempt is discarded transactionally"
        );
        assert!(
            observed.peak_live_slots <= slot_bound,
            "the production vectors fit v7's 2*T+4*B slot peak"
        );
        assert_eq!(
            observed.peak_live_bytes,
            observed.peak_live_slots * 256,
            "the measured slot peak uses the v7 256-byte accounting unit"
        );
        assert!(
            observed.peak_live_bytes <= 256 * slot_bound,
            "the production slot memory fits v7's 256*(2*T+4*B) component"
        );
        // Picos de fase de la admisión que el mismo documento produjo: el
        // fixture instrumental cierra así el desglose que .6 exige, desde la
        // carga hasta la fórmula final, sobre el receipt sellado.
        let phases = preflight.receipt().phases();
        assert_eq!(
            phases
                .iter()
                .map(crate::f5_admission::AdmissionPhaseRecord::phase)
                .collect::<Vec<_>>(),
            crate::f5_admission::ACCOUNTED_ADMISSION_PHASES.to_vec(),
            "las nueve fases se publican en el orden cerrado de la máquina"
        );
        assert!(
            phases
                .windows(2)
                .all(|pair| pair[0].operation_units_cumulative()
                    <= pair[1].operation_units_cumulative()),
            "el acumulado de operaciones por fase es monótono"
        );
        assert_eq!(
            phases
                .last()
                .expect("nueve fases")
                .operation_units_cumulative(),
            preflight.operation_units(),
            "la última fase cierra exactamente en la cota admitida"
        );
        assert_eq!(
            phases
                .iter()
                .map(crate::f5_admission::AdmissionPhaseRecord::phase_peak_bytes)
                .max()
                .expect("nueve picos de fase"),
            preflight.working_memory_bytes(),
            "el pico publicado es el mayor pico de fase, no una suma"
        );
        assert!(
            observed
                .placement_trace
                .iter()
                .any(|(attempt, _, _, left, right)| *attempt == 0 && *left && *right),
            "the failed attempt contains an interior winner with two residues"
        );
        assert!(
            observed
                .placement_trace
                .iter()
                .any(|(attempt, index, _, _, _)| *attempt == 0 && *index >= 2),
            "the failed attempt selects the final initial tract"
        );

        // La frontera exacta−1 conserva la anchura del literal, de modo que el
        // documento consume la misma cota y sólo cambia el límite declarado.
        // En v7 la admisión ocurre dentro de la carga: el rechazo llega antes
        // de que exista una admisión sellada que recuperar.
        let rejected_general = write_production_layout_fixture(&root, exact_limit - 1, seed_byte);
        let rejected = load_and_validate_f5_config(&rejected_general)
            .expect_err("v7 rejects operation_limit one below the exact admitted bound");
        assert_eq!(rejected.code(), "GPS-RES-101");
        assert_eq!(rejected.context()["dimension"], "operation_units");
        assert_eq!(rejected.context()["accumulated"], exact_limit.to_string());
        assert_eq!(
            rejected.context()["limit"],
            (exact_limit - 1).to_string(),
            "el rechazo cita el límite declarado, no la cota"
        );
        fs::remove_dir_all(root).expect("production fixture cleanup");
    }

    #[test]
    fn site_bound_is_enforced_even_when_preflight_input_is_internally_inconsistent() {
        let plan = fixed_plan(1, 4097);
        let gaps = gap_distribution(1);
        let sequence = sequence(0, "contig-a", &vec![b'A'; 4097]);
        let error = run_layout(
            &plan,
            &gaps,
            &[sequence],
            &[F5EligibleTract::new(0, 1, 4097)],
            1,
            root_for_layout(),
        )
        .expect_err("the dynamic S bound is still a F5 resource rejection");
        assert_eq!(error.code(), "GPS-RES-101");
    }
}
