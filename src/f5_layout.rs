//! Longitudes de bloque F5 anteriores a la colocación física.
//!
//! F5.1.1 solo resuelve la geometría de las longitudes. No elige tramos,
//! posiciones ni SNPs: esos contratos pertenecen a F5.1.2 y posteriores.

use crate::application::EffectiveF5Config;
use crate::diagnostics::Diagnostic;
use crate::f5_preflight::{F5ResourcePreflight, resolved_target_for_config};
use crate::f5_profile::F5Profile;
use crate::identity::ScenarioId;
use crate::randomness::{F5RandomDecision, RandomWordSource, RootSeed, derive_f5_subseed};

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
        preflight: &F5ResourcePreflight,
    ) -> Result<Self, Diagnostic> {
        let scenario_id = ScenarioId::parse("scenario_id", config.scenario_id()).map_err(|_| {
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
        if preflight.scenario_id() != config.scenario_id()
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
    pub fn sample(&self, root_seed: RootSeed, attempt_index: u64) -> Result<Vec<u64>, Diagnostic> {
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
    use super::{
        F5BlockLengthBand, F5BlockLengthDistribution, F5BlockLengthPlan, F5LengthFeasibility,
    };
    use crate::identity::ScenarioId;
    use crate::randomness::{RandomWordSource, RootSeed};

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
}
