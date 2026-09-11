//! Evaluador LD de F5: `p_hat`, `D`, `r²`, bins y signos de un intento.
//!
//! DTC-089 «LD, búsqueda y terminación» fija la aritmética literal y su orden.
//! En orden canónico `h`, `i`, `j`, y con suma de Neumaier **sin reducción
//! paralela**:
//!
//! - `p_hat_i = Σ_h q_h x_hi`,
//! - `p_hat_ij = Σ_h q_h x_hi x_hj`,
//! - `D_ij = p_hat_ij - p_hat_i*p_hat_j`,
//! - `r²_ij = D_ij² / [p_hat_i(1-p_hat_i) * p_hat_j(1-p_hat_j)]`.
//!
//! Las expresiones se escriben exactamente en ese orden y no usan `mul_add`:
//! un cambio de asociatividad o una contracción FMA cambiarían el último bit y
//! romperían la reproducibilidad. Rust no contrae multiplicación y suma por su
//! cuenta, de modo que basta con no pedirlo.
//!
//! La agrupación del denominador es la que escriben los corchetes de DTC-089:
//! `[p_hat_i(1-p_hat_i)] * [p_hat_j(1-p_hat_j)]`, es decir `(a*b)*(c*d)`. No
//! es una elección de estilo: la alternativa `((a*b)*c)*d` devuelve otro
//! último bit en una tercera parte de las parejas.
//!
//! Las parejas son **solo intrabloque**: el evaluador recorre `i < j` dentro de
//! un mismo bloque×subpoblación y nunca cruza bloques. Cada bin agrega las
//! parejas evaluables cuya distancia física cae en él y su media se reduce en
//! el orden `(site_index_i, site_index_j)`.
//!
//! Dos estados que **no** se colapsan en un cero:
//!
//! - `not_evaluable` es de pareja: el denominador de `r²` es cero porque una
//!   marginal realizada vale 0 o 1 —ahí caen necesariamente los sitios
//!   fijados—. Nunca se publica `r² = 0` en su lugar.
//! - `not_applicable` es de bin: ninguna pareja evaluable cayó en él. No es un
//!   error ni se sustituye por otro bin; si **todos** los bins de un
//!   bloque×subpoblación lo son, el intento es inválido.
//!
//! Junto a los bins, cada bloque×subpoblación publica los tres agregados que
//! `f5-attempts.tsv` exige —`bins_outside`, `max_r2_error` y `sum_r2_error`—
//! porque son funciones puras de los bins: el error de un bin es
//! `|media - objetivo|` contra la tolerancia inclusiva del perfil y la suma se
//! reduce con Neumaier en orden de bin. La búsqueda de F5.2.3 los compara; no
//! los recalcula ni los define.
//!
//! El veredicto de invalidez se publica, no se lanza: quien recorre los
//! intentos —la búsqueda de F5.2.3— es quien decide si corresponde
//! `GPS-DOMAIN-102` o `GPS-DOMAIN-104`. Este módulo tampoco consume entropía:
//! es una función determinista de la matriz, los pesos, las polaridades, las
//! posiciones y el perfil LD.

use std::collections::BTreeMap;
use std::fmt;

use crate::application::EffectiveF5Config;
use crate::diagnostics::Diagnostic;
use crate::f5_haplotypes::{F5AlleleMatrix, F5BlockHaplotypes, F5HaplotypeAttempt};
use crate::f5_layout::{F5LayoutBlock, F5PhysicalLayout};
use crate::f5_marginals::{F5RowVisits, NeumaierSum};
use crate::f5_profile::LdBin;
use crate::identity::{F5BlockId, PopulationId};

/// Literales de `status` en `f5-ld-evaluation.tsv`.
const EVALUATED: &str = "evaluated";
const NOT_APPLICABLE: &str = "not_applicable";

/// Causas por las que una pareja no puede evaluarse **aritméticamente**.
///
/// No describen el estado `not_evaluable` del contrato —ese es un resultado
/// legítimo y viaja dentro de [`F5SitePair`]—, sino un desajuste entre la
/// matriz, los pesos y los índices recibidos. El llamador debe poder
/// distinguirlas: dos de ellas acusan a la geometría y la tercera a la
/// aritmética.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F5LdError {
    /// La pareja nombra un sitio que la matriz no tiene.
    SiteOutOfRange,
    /// El perfil declara un número de haplotipos distinto al de la matriz.
    WeightCountMismatch,
    /// La suma compensada dejó de ser finita.
    NonFiniteSum,
}

impl F5LdError {
    #[must_use]
    pub const fn detail(self) -> &'static str {
        match self {
            Self::SiteOutOfRange => "la pareja nombra un sitio fuera de la matriz del bloque",
            Self::WeightCountMismatch => {
                "el perfil declara un número de haplotipos distinto al de la matriz"
            }
            Self::NonFiniteSum => "la suma compensada de la pareja no es finita",
        }
    }
}

impl fmt::Display for F5LdError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.detail())
    }
}

/// Estadísticos de una pareja intrabloque `(i, j)` con `i < j`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F5SitePair {
    site_index_i: u64,
    site_index_j: u64,
    p_hat_i: f64,
    p_hat_j: f64,
    p_hat_ij: f64,
    d: f64,
    r_squared: Option<f64>,
}

impl F5SitePair {
    #[must_use]
    pub const fn site_index_i(&self) -> u64 {
        self.site_index_i
    }

    #[must_use]
    pub const fn site_index_j(&self) -> u64 {
        self.site_index_j
    }

    /// `p_hat_i = Σ_h q_h x_hi`.
    #[must_use]
    pub const fn p_hat_i(&self) -> f64 {
        self.p_hat_i
    }

    /// `p_hat_j = Σ_h q_h x_hj`.
    #[must_use]
    pub const fn p_hat_j(&self) -> f64 {
        self.p_hat_j
    }

    /// `p_hat_ij = Σ_h q_h x_hi x_hj`.
    #[must_use]
    pub const fn p_hat_ij(&self) -> f64 {
        self.p_hat_ij
    }

    /// `D_ij = p_hat_ij - p_hat_i*p_hat_j`. Existe siempre, incluso cuando la
    /// pareja no es evaluable.
    #[must_use]
    pub const fn d(&self) -> f64 {
        self.d
    }

    /// `r²_ij`, o `None` cuando el denominador es cero: `not_evaluable` nunca
    /// se representa como `r² = 0`.
    #[must_use]
    pub const fn r_squared(&self) -> Option<f64> {
        self.r_squared
    }

    #[must_use]
    pub const fn is_evaluable(&self) -> bool {
        self.r_squared.is_some()
    }

    /// Regla de signos de DTC-089 sobre una pareja evaluable: `D>0` para
    /// polaridades iguales y `D<0` para distintas. `D=0` incumple. Una pareja
    /// no evaluable no tiene veredicto.
    #[must_use]
    pub fn satisfies_sign_rule(&self, polarities_match: bool) -> Option<bool> {
        self.is_evaluable()
            .then(|| sign_rule_is_satisfied(self.d, polarities_match))
    }
}

/// `D>0` si las polaridades coinciden, `D<0` si difieren; `D=0` incumple
/// siempre. No hay cuota de signos.
pub(crate) const fn sign_rule_is_satisfied(d: f64, polarities_match: bool) -> bool {
    if polarities_match { d > 0.0 } else { d < 0.0 }
}

/// Evalúa una pareja intrabloque a partir de la matriz y los pesos.
///
/// Recalcula las dos marginales además de la conjunta, de modo que la pareja
/// queda definida solo por `(q, matriz)`. El agregador de un bloque usa la
/// misma aritmética pero reutiliza las S marginales que ya calculó una vez.
///
/// # Errors
///
/// Devuelve [`F5LdError`] si los índices no pertenecen a la matriz, si los
/// pesos no la acompañan o si alguna suma deja de ser finita o de pertenecer a
/// `[0,1]`.
pub fn evaluate_f5_site_pair(
    weights: &[f64],
    matrix: &F5AlleleMatrix,
    site_index_i: usize,
    site_index_j: usize,
) -> Result<F5SitePair, F5LdError> {
    let mut visits = F5RowVisits::default();
    let p_hat_i = site_marginal(weights, matrix, site_index_i, &mut visits)?;
    let p_hat_j = site_marginal(weights, matrix, site_index_j, &mut visits)?;
    let p_hat_ij = joint_marginal(weights, matrix, site_index_i, site_index_j, &mut visits)?;
    Ok(pair_from_components(
        site_index_i as u64,
        site_index_j as u64,
        p_hat_i,
        p_hat_j,
        p_hat_ij,
    ))
}

/// `D` y `r²` a partir de las tres marginales, en el orden escrito de DTC-089.
///
/// Dos decisiones locales registradas, ninguna de las cuales toca el orden ni
/// la agrupación: `factor_i * factor_j` es el mismo árbol de expresión que
/// `p_hat_i(1-p_hat_i) * p_hat_j(1-p_hat_j)`.
///
/// 1. **La evaluabilidad se decide por factor, no por el producto.** El
///    contrato dice «denominador cero da `not_evaluable`». Una marginal puede
///    exceder 1 por un ULP —`normalize_slice` solo exige `|suma-1| <= 1e-9`
///    después de dividir, así que la suma compensada de una columna
///    monomórfica puede valer `1.0000000000000002`— y entonces su factor es
///    negativo. Con dos columnas así, el **producto** de dos factores
///    negativos es positivo y declararía evaluable una pareja sin variación,
///    con `r²` cercano a 1 y un signo incumplido. Exigir `factor > 0` en cada
///    lado cubre el cero del contrato y ese caso degenerado con la misma
///    condición. Una marginal fuera de `[0,1]` no aborta: DTC-089 no contempla
///    abortar ahí.
/// 2. **`r²` se recorta a 1 donde nace.** En el acoplamiento perfecto
///    `D = fl(p - fl(p*p))` y el factor es `fl(p * fl(1-p))`: dos redondeos
///    distintos del mismo real, cuyo cociente supera 1 en una fracción
///    apreciable de los pesos decimales —máximo observado `1.0000000000442435`
///    con `w = 0,999999`—. El TSV `f5-ld-evaluation.tsv` exige medias y errores
///    en `[0,1]`, así que un solo `min` en el origen deja dentro del dominio
///    tanto `r²` como la media y el error del bin, sin recortarlos aparte.
pub(crate) fn pair_from_components(
    site_index_i: u64,
    site_index_j: u64,
    p_hat_i: f64,
    p_hat_j: f64,
    p_hat_ij: f64,
) -> F5SitePair {
    let d = p_hat_ij - p_hat_i * p_hat_j;
    let factor_i = f5_evaluability_factor(p_hat_i);
    let factor_j = f5_evaluability_factor(p_hat_j);
    let denominator = factor_i * factor_j;
    let r_squared = (factor_i > 0.0 && factor_j > 0.0).then(|| (d * d / denominator).min(1.0));
    F5SitePair {
        site_index_i,
        site_index_j,
        p_hat_i,
        p_hat_j,
        p_hat_ij,
        d,
        r_squared,
    }
}

/// Factor `p(1-p)` del denominador de `r²`, escrito una sola vez.
///
/// La evaluabilidad de una pareja es `factor_i > 0 && factor_j > 0`, de modo
/// que `factor == 0.0` es la **frontera exacta** de esa regla. La búsqueda de
/// F5.2.3 la necesita con la misma aritmética que el evaluador para poder
/// atribuir una divergencia a la frontera en lugar de a un defecto.
pub(crate) fn f5_evaluability_factor(marginal: f64) -> f64 {
    marginal * (1.0 - marginal)
}

/// Veredicto canónico de una pareja intrabloque.
///
/// Existe para que la búsqueda de F5.2.3 pueda contrastar su puntuación
/// incremental **pareja a pareja** contra el canon y atribuir cada divergencia
/// a una frontera exacta del contrato. Se rellena durante la evaluación
/// canónica que ya se hace de todos modos, de modo que no cuesta ninguna visita
/// de fila adicional.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct F5CanonicalPairVerdict {
    pub(crate) d: f64,
    pub(crate) factor_i: f64,
    pub(crate) factor_j: f64,
    pub(crate) evaluable: bool,
    pub(crate) sign_ok: bool,
}

/// `p_hat_i = Σ_h q_h x_hi` en orden canónico `h` y con Neumaier.
///
/// Un `x_hi` nulo aporta `0*q_h`, que ni mueve la suma ni la compensación, de
/// modo que omitirlo produce exactamente los mismos bits que sumarlo.
///
/// `visits` cuenta una visita por fila recorrida en **el medidor que recibe**.
/// Es la instrumentación con la que la búsqueda de F5.2.3 acredita que sus tres
/// sumas canónicas por intento son las que dice y ninguna más: si alguien
/// añadiera una llamada a esta función dentro del bucle de propuestas pasándole
/// el medidor del intento, el recuento crecería con el presupuesto y la
/// identidad publicada fallaría. No cubre un recorrido que declare su propio
/// contador.
pub(crate) fn site_marginal(
    weights: &[f64],
    matrix: &F5AlleleMatrix,
    site_index: usize,
    visits: &mut F5RowVisits,
) -> Result<f64, F5LdError> {
    ensure_shape(weights, matrix, site_index)?;
    let mut compensated = NeumaierSum::default();
    let mut total = 0.0_f64;
    for (haplotype, weight) in weights.iter().enumerate() {
        visits.add(1);
        if matrix
            .allele(haplotype, site_index)
            .ok_or(F5LdError::SiteOutOfRange)?
        {
            total = compensated
                .add(*weight)
                .map_err(|_| F5LdError::NonFiniteSum)?;
        }
    }
    finite_total(total)
}

/// `p_hat_ij = Σ_h q_h x_hi x_hj` en orden canónico `h` y con Neumaier.
///
/// `visits` cuenta una visita por fila recorrida, igual que en
/// [`site_marginal`].
pub(crate) fn joint_marginal(
    weights: &[f64],
    matrix: &F5AlleleMatrix,
    site_index_i: usize,
    site_index_j: usize,
    visits: &mut F5RowVisits,
) -> Result<f64, F5LdError> {
    ensure_shape(weights, matrix, site_index_i)?;
    ensure_shape(weights, matrix, site_index_j)?;
    let mut compensated = NeumaierSum::default();
    let mut total = 0.0_f64;
    for (haplotype, weight) in weights.iter().enumerate() {
        visits.add(1);
        let alt_i = matrix
            .allele(haplotype, site_index_i)
            .ok_or(F5LdError::SiteOutOfRange)?;
        let alt_j = matrix
            .allele(haplotype, site_index_j)
            .ok_or(F5LdError::SiteOutOfRange)?;
        if alt_i && alt_j {
            total = compensated
                .add(*weight)
                .map_err(|_| F5LdError::NonFiniteSum)?;
        }
    }
    finite_total(total)
}

fn ensure_shape(
    weights: &[f64],
    matrix: &F5AlleleMatrix,
    site_index: usize,
) -> Result<(), F5LdError> {
    if weights.len() != matrix.haplotype_count() {
        return Err(F5LdError::WeightCountMismatch);
    }
    if site_index >= matrix.site_count() {
        return Err(F5LdError::SiteOutOfRange);
    }
    Ok(())
}

/// Una marginal solo se rechaza si deja de ser finita. Que exceda `1` por un
/// ULP es un resultado legítimo de un perfil admitido y lo absorbe la
/// evaluabilidad por factor de [`pair_from_components`].
fn finite_total(value: f64) -> Result<f64, F5LdError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(F5LdError::NonFiniteSum)
    }
}

/// Un bin del perfil LD ya evaluado sobre un bloque×subpoblación.
///
/// Reproduce fila a fila `f5-ld-evaluation.tsv`: `pair_count` cuenta solo
/// parejas **evaluables** y vale cero si y solo si el estado es
/// `not_applicable`, caso en que `mean_r2` y `abs_error` no existen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F5LdBinEvaluation {
    distance_start_bp: u64,
    distance_end_bp: u64,
    target_mean_r2: f64,
    max_abs_error: f64,
    pair_count: u64,
    mean_r2: Option<f64>,
    abs_error: Option<f64>,
}

impl F5LdBinEvaluation {
    #[must_use]
    pub const fn distance_start_bp(&self) -> u64 {
        self.distance_start_bp
    }

    #[must_use]
    pub const fn distance_end_bp(&self) -> u64 {
        self.distance_end_bp
    }

    #[must_use]
    pub const fn target_mean_r2(&self) -> f64 {
        self.target_mean_r2
    }

    /// Tolerancia absoluta del perfil, inclusiva.
    #[must_use]
    pub const fn max_abs_error(&self) -> f64 {
        self.max_abs_error
    }

    #[must_use]
    pub const fn pair_count(&self) -> u64 {
        self.pair_count
    }

    /// Media de `r²` de las parejas evaluables del bin, reducida en orden
    /// `(site_index_i, site_index_j)`; `None` es `not_applicable`.
    #[must_use]
    pub const fn mean_r2(&self) -> Option<f64> {
        self.mean_r2
    }

    /// `|mean_r2 - target_mean_r2|`, o `None` si el bin no es aplicable.
    #[must_use]
    pub const fn abs_error(&self) -> Option<f64> {
        self.abs_error
    }

    #[must_use]
    pub const fn is_evaluated(&self) -> bool {
        self.mean_r2.is_some()
    }

    /// Literal `evaluated`/`not_applicable` del TSV normativo.
    #[must_use]
    pub const fn status(&self) -> &'static str {
        if self.is_evaluated() {
            EVALUATED
        } else {
            NOT_APPLICABLE
        }
    }

    /// La tolerancia es inclusiva. Un bin `not_applicable` no está fuera de
    /// tolerancia: no es un error.
    #[must_use]
    pub fn is_within_tolerance(&self) -> bool {
        self.abs_error
            .is_none_or(|error| error <= self.max_abs_error)
    }
}

/// Evaluación LD completa de un bloque×subpoblación.
#[derive(Clone, Debug, PartialEq)]
pub struct F5BlockLdEvaluation {
    population_id: PopulationId,
    block_index: u64,
    block_id: F5BlockId,
    bins: Vec<F5LdBinEvaluation>,
    evaluable_pairs: u64,
    not_evaluable_pairs: u64,
    signs_failed: u64,
    bins_outside: u64,
    max_r2_error: f64,
    sum_r2_error: f64,
}

impl F5BlockLdEvaluation {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }

    #[must_use]
    pub const fn block_index(&self) -> u64 {
        self.block_index
    }

    #[must_use]
    pub const fn block_id(&self) -> &F5BlockId {
        &self.block_id
    }

    /// Bins en el orden del perfil, uno por fila declarada.
    #[must_use]
    pub fn bins(&self) -> &[F5LdBinEvaluation] {
        &self.bins
    }

    /// Parejas intrabloque con denominador no nulo.
    #[must_use]
    pub const fn evaluable_pairs(&self) -> u64 {
        self.evaluable_pairs
    }

    /// Parejas intrabloque `not_evaluable`: ahí caen los sitios fijados.
    #[must_use]
    pub const fn not_evaluable_pairs(&self) -> u64 {
        self.not_evaluable_pairs
    }

    /// Parejas evaluables que incumplen la regla de signos.
    #[must_use]
    pub const fn signs_failed(&self) -> u64 {
        self.signs_failed
    }

    /// Bins evaluados cuyo error supera la tolerancia inclusiva del perfil.
    #[must_use]
    pub const fn bins_outside(&self) -> u64 {
        self.bins_outside
    }

    /// Mayor `abs_error` entre los bins evaluados; `0` si no hay ninguno.
    #[must_use]
    pub const fn max_r2_error(&self) -> f64 {
        self.max_r2_error
    }

    /// Suma compensada de `abs_error` en orden de bin; `0` si no hay ninguno.
    #[must_use]
    pub const fn sum_r2_error(&self) -> f64 {
        self.sum_r2_error
    }

    /// Un bloque×subpoblación sin ninguna pareja evaluable invalida el intento.
    #[must_use]
    pub const fn is_valid(&self) -> bool {
        self.evaluable_pairs > 0
    }
}

#[cfg(test)]
impl F5BlockLdEvaluation {
    /// Constructor literal para las pruebas de la guarda de divergencia de
    /// F5.2.3. La vía productiva sigue siendo [`evaluate_block`].
    /// Añade una fila de bin ya evaluada, para las pruebas de la guarda.
    pub(crate) fn push_bin_for_test(
        &mut self,
        definition: LdBin,
        pair_count: u64,
        mean_r2: f64,
        abs_error: f64,
    ) {
        self.bins.push(F5LdBinEvaluation {
            distance_start_bp: definition.distance_start_bp(),
            distance_end_bp: definition.distance_end_bp(),
            target_mean_r2: definition.target_mean_r2(),
            max_abs_error: definition.max_abs_error(),
            pair_count,
            mean_r2: Some(mean_r2),
            abs_error: Some(abs_error),
        });
    }

    pub(crate) fn for_test(
        population_id: PopulationId,
        block_id: F5BlockId,
        evaluable_pairs: u64,
        signs_failed: u64,
        bins_outside: u64,
    ) -> Self {
        Self {
            population_id,
            block_index: 0,
            block_id,
            bins: Vec::new(),
            evaluable_pairs,
            not_evaluable_pairs: 0,
            signs_failed,
            bins_outside,
            max_r2_error: 0.0,
            sum_r2_error: 0.0,
        }
    }
}

/// Evaluación LD de un intento completo, en el orden canónico de la tabla
/// haplotípica.
#[derive(Clone, Debug, PartialEq)]
pub struct F5LdEvaluation {
    attempt_index: u64,
    blocks: Vec<F5BlockLdEvaluation>,
}

impl F5LdEvaluation {
    #[must_use]
    pub const fn attempt_index(&self) -> u64 {
        self.attempt_index
    }

    #[must_use]
    pub fn blocks(&self) -> &[F5BlockLdEvaluation] {
        &self.blocks
    }

    #[must_use]
    pub fn block(
        &self,
        population_id: &PopulationId,
        block_index: u64,
    ) -> Option<&F5BlockLdEvaluation> {
        self.blocks
            .iter()
            .find(|row| row.population_id == *population_id && row.block_index == block_index)
    }

    /// El intento es inválido en cuanto un bloque×subpoblación no tiene
    /// ninguna pareja evaluable, es decir, cuando todos sus bins son
    /// `not_applicable`.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.blocks.iter().all(F5BlockLdEvaluation::is_valid)
    }
}

/// Evalúa el LD realizado de un intento haplotípico.
///
/// Recibe el layout porque las distancias de los bins son físicas y viven en
/// él, no en la matriz. No consume entropía ni reutiliza estado del
/// constructor: recalcula cada `p_hat_i` desde la matriz y exige que coincida
/// bit a bit con el que la columna publicó, de modo que una divergencia entre
/// constructor y evaluador se detiene aquí en lugar de propagarse.
///
/// # Errors
///
/// Devuelve `GPS-CONF-101` si el intento y el layout no describen la misma
/// geometría —identificadores de bloque o de sitio distintos, otro número de
/// sitios— o si el perfil no declara bins para una subpoblación, y
/// `GPS-INTERNAL-101` si la aritmética de DTC-089 deja de estar definida o una
/// distancia no cae en ningún bin.
pub fn evaluate_f5_linkage_disequilibrium(
    config: &EffectiveF5Config,
    layout: &F5PhysicalLayout,
    attempt: &F5HaplotypeAttempt,
) -> Result<F5LdEvaluation, Diagnostic> {
    let scenario_id = config.scenario_id();
    let scenario_path = config.scenario_path().display().to_string();
    let weights = config.profile().haplotype_weights();
    let foreign = |detail: &'static str| {
        Diagnostic::f5_configuration(scenario_path.clone(), "layout.identity", detail)
    };

    let mut layout_blocks = BTreeMap::new();
    for block in layout.blocks() {
        if layout_blocks.insert(block.block_index(), block).is_some() {
            return Err(foreign("el layout repite un índice de bloque"));
        }
    }

    let mut blocks = Vec::with_capacity(attempt.blocks().len());
    for block in attempt.blocks() {
        let layout_block = layout_blocks
            .get(&block.block_index())
            .ok_or_else(|| foreign("el intento nombra un bloque que este layout no contiene"))?;
        let bins = config
            .profile()
            .ld_bins(block.population_id().as_str())
            .ok_or_else(|| {
                Diagnostic::f5_configuration(
                    scenario_path.clone(),
                    "profile.ld_profile",
                    "el perfil no declara bins LD para la subpoblación del intento",
                )
            })?;
        blocks.push(evaluate_block(
            scenario_id,
            &foreign,
            weights,
            bins,
            layout_block,
            block,
            None,
            &mut F5RowVisits::default(),
        )?);
    }
    Ok(F5LdEvaluation {
        attempt_index: attempt.attempt_index(),
        blocks,
    })
}

/// Evalúa el LD realizado de **un solo** bloque×subpoblación.
///
/// Existe porque DTC-089 cuenta los intentos por bloque×subpoblación: la
/// búsqueda de F5.2.3 resuelve cada unidad por separado y necesita su veredicto
/// canónico sin evaluar las demás. Produce exactamente la misma fila que
/// [`evaluate_f5_linkage_disequilibrium`] para ese bloque.
///
/// # Errors
///
/// Devuelve `GPS-CONF-101` si el bloque y el fragmento de layout no describen
/// la misma geometría o si el perfil no declara bins para la subpoblación, y
/// `GPS-INTERNAL-101` si la aritmética de DTC-089 deja de estar definida o una
/// distancia no cae en ningún bin.
pub fn evaluate_f5_block_linkage_disequilibrium(
    config: &EffectiveF5Config,
    layout_block: &F5LayoutBlock,
    block: &F5BlockHaplotypes,
) -> Result<F5BlockLdEvaluation, Diagnostic> {
    evaluate_f5_block_with_pair_verdicts(
        config,
        layout_block,
        block,
        None,
        &mut F5RowVisits::default(),
    )
}

/// Igual que [`evaluate_f5_block_linkage_disequilibrium`], recogiendo además
/// el veredicto canónico de cada pareja en el orden `(i, j)` creciente.
///
/// Es la vía que usa la búsqueda de F5.2.3 para atribuir las divergencias de su
/// puntuación incremental. No repite ningún recorrido: los veredictos se
/// apuntan mientras la evaluación canónica recorre las parejas.
pub(crate) fn evaluate_f5_block_with_pair_verdicts(
    config: &EffectiveF5Config,
    layout_block: &F5LayoutBlock,
    block: &F5BlockHaplotypes,
    verdicts: Option<&mut Vec<F5CanonicalPairVerdict>>,
    visits: &mut F5RowVisits,
) -> Result<F5BlockLdEvaluation, Diagnostic> {
    let scenario_path = config.scenario_path().display().to_string();
    let foreign = |detail: &'static str| {
        Diagnostic::f5_configuration(scenario_path.clone(), "layout.identity", detail)
    };
    if layout_block.block_index() != block.block_index() {
        return Err(foreign(
            "el bloque del intento y el del layout no comparten índice",
        ));
    }
    let bins = config
        .profile()
        .ld_bins(block.population_id().as_str())
        .ok_or_else(|| {
            Diagnostic::f5_configuration(
                scenario_path.clone(),
                "profile.ld_profile",
                "el perfil no declara bins LD para la subpoblación del intento",
            )
        })?;
    evaluate_block(
        config.scenario_id(),
        &foreign,
        config.profile().haplotype_weights(),
        bins,
        layout_block,
        block,
        verdicts,
        visits,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "el evaluador necesita identidad, perfil, layout, bloque y, además, los dos canales que la búsqueda de F5.2.3 le pide: los veredictos por pareja y el medidor de visitas de fila del intento"
)]
fn evaluate_block(
    scenario_id: &str,
    foreign: &impl Fn(&'static str) -> Diagnostic,
    weights: &[f64],
    bins: &[LdBin],
    layout_block: &F5LayoutBlock,
    block: &F5BlockHaplotypes,
    mut verdicts: Option<&mut Vec<F5CanonicalPairVerdict>>,
    visits: &mut F5RowVisits,
) -> Result<F5BlockLdEvaluation, Diagnostic> {
    let internal = |detail: String| Diagnostic::f5_internal(scenario_id, detail);
    let arithmetic = |error: F5LdError| Diagnostic::f5_internal(scenario_id, error.detail());

    if layout_block.block_id() != block.block_id() {
        return Err(foreign(
            "el bloque del intento y el del layout no comparten identificador",
        ));
    }
    let site_count = block.matrix().site_count();
    if block.columns().len() != site_count || layout_block.sites().len() != site_count {
        return Err(foreign(
            "el número de sitios del intento y el del layout no coinciden",
        ));
    }

    // Marginales realizadas: una sola pasada por columna. Recalcularlas dentro
    // de cada pareja daría los mismos bits pero triplicaría el recorrido de H.
    let mut marginals = Vec::with_capacity(site_count);
    let mut positions = Vec::with_capacity(site_count);
    for (index, (column, site)) in block
        .columns()
        .iter()
        .zip(layout_block.sites().iter())
        .enumerate()
    {
        if column.site_index() != index as u64 || site.site_index() != index as u64 {
            return Err(foreign(
                "un sitio no ocupa la columna que su índice anuncia",
            ));
        }
        if column.site_id() != site.site_id() {
            return Err(foreign(
                "el sitio del intento y el del layout no comparten identificador",
            ));
        }
        let realized = site_marginal(weights, block.matrix(), index, visits).map_err(arithmetic)?;
        ensure_published_column_matches(index, column.realized_alt_af(), realized)
            .map_err(internal)?;
        marginals.push(realized);
        positions.push(site.position_1based());
    }

    let mut accumulators = vec![BinAccumulator::default(); bins.len()];
    let mut evaluable_pairs = 0_u64;
    let mut not_evaluable_pairs = 0_u64;
    let mut signs_failed = 0_u64;

    // Orden canónico de parejas y ausencia de pares interbloque: los dos
    // bucles viven dentro de un único bloque y recorren `i < j` creciente.
    for i in 0..site_count {
        for j in (i + 1)..site_count {
            let pair = pair_from_components(
                i as u64,
                j as u64,
                marginals[i],
                marginals[j],
                joint_marginal(weights, block.matrix(), i, j, visits).map_err(arithmetic)?,
            );
            let polarities_match = block.columns()[i].polarity() == block.columns()[j].polarity();
            if let Some(collected) = verdicts.as_deref_mut() {
                collected.push(F5CanonicalPairVerdict {
                    d: pair.d(),
                    factor_i: f5_evaluability_factor(marginals[i]),
                    factor_j: f5_evaluability_factor(marginals[j]),
                    evaluable: pair.r_squared().is_some(),
                    sign_ok: pair.r_squared().is_some()
                        && sign_rule_is_satisfied(pair.d(), polarities_match),
                });
            }
            let Some(r_squared) = pair.r_squared() else {
                not_evaluable_pairs += 1;
                continue;
            };
            evaluable_pairs += 1;
            if !sign_rule_is_satisfied(pair.d(), polarities_match) {
                signs_failed += 1;
            }
            let distance = positions[j]
                .checked_sub(positions[i])
                .filter(|gap| *gap > 0)
                .ok_or_else(|| {
                    internal(format!(
                        "los sitios {i} y {j} no están en posiciones físicas crecientes"
                    ))
                })?;
            let bin = bin_of(bins, distance).ok_or_else(|| {
                internal(format!(
                    "la distancia {distance} de la pareja ({i},{j}) no cae en ningún bin del perfil"
                ))
            })?;
            accumulators[bin]
                .add(r_squared)
                .map_err(|detail| internal(detail.to_owned()))?;
        }
    }

    let mut evaluated_bins = Vec::with_capacity(bins.len());
    let mut bins_outside = 0_u64;
    let mut max_r2_error = 0.0_f64;
    let mut error_sum = NeumaierSum::default();
    let mut sum_r2_error = 0.0_f64;
    for (definition, accumulator) in bins.iter().zip(accumulators) {
        let pair_count = accumulator.pair_count;
        let (mean_r2, abs_error) = if pair_count == 0 {
            (None, None)
        } else {
            let mean = accumulator.total / pair_count as f64;
            if !mean.is_finite() {
                return Err(internal(format!(
                    "la media de {pair_count} parejas del bin {} no es finita",
                    definition.distance_start_bp()
                )));
            }
            let error = f5_bin_abs_error(mean, definition);
            if f5_bin_is_outside(error, definition) {
                bins_outside += 1;
            }
            if error > max_r2_error {
                max_r2_error = error;
            }
            sum_r2_error = error_sum
                .add(error)
                .map_err(|detail| internal(detail.to_owned()))?;
            (Some(mean), Some(error))
        };
        evaluated_bins.push(F5LdBinEvaluation {
            distance_start_bp: definition.distance_start_bp(),
            distance_end_bp: definition.distance_end_bp(),
            target_mean_r2: definition.target_mean_r2(),
            max_abs_error: definition.max_abs_error(),
            pair_count,
            mean_r2,
            abs_error,
        });
    }

    Ok(F5BlockLdEvaluation {
        population_id: block.population_id().clone(),
        block_index: block.block_index(),
        block_id: block.block_id().clone(),
        bins: evaluated_bins,
        evaluable_pairs,
        not_evaluable_pairs,
        signs_failed,
        bins_outside,
        max_r2_error,
        sum_r2_error,
    })
}

/// `|media - objetivo|` de un bin evaluado.
///
/// La regla se escribe una sola vez: la usan el evaluador canónico y la
/// puntuación incremental de la búsqueda de F5.2.3, que si no la duplicaría.
pub(crate) fn f5_bin_abs_error(mean: f64, definition: &LdBin) -> f64 {
    (mean - definition.target_mean_r2()).abs()
}

/// La tolerancia del perfil es **inclusiva**: solo un error estrictamente
/// mayor deja el bin fuera.
pub(crate) fn f5_bin_is_outside(error: f64, definition: &LdBin) -> bool {
    error > definition.max_abs_error()
}

/// La columna publicada debe ser **bit a bit** la suma canónica de su matriz.
///
/// Es la barrera que separa las dos aritméticas de F5.2.3: la puntuación del
/// bucle es incremental, pero lo que se publica se recomputa, y si alguien
/// olvidase recomputarlo el evaluador lo vería aquí en lugar de dejar pasar una
/// columna que su matriz ya no realiza.
fn ensure_published_column_matches(
    index: usize,
    published: f64,
    recomputed: f64,
) -> Result<(), String> {
    if published == recomputed {
        return Ok(());
    }
    Err(format!(
        "la columna {index} publicó {published} y el evaluador recalculó {recomputed}"
    ))
}

/// Acumulador de un bin mientras se recorren las parejas.
///
/// `total` conserva el valor que devolvió la última llamada a
/// [`NeumaierSum::add`]; el acumulador no expone su total por separado y
/// releerlo desde `sum` perdería la compensación.
#[derive(Clone, Copy, Default)]
struct BinAccumulator {
    compensated: NeumaierSum,
    total: f64,
    pair_count: u64,
}

impl BinAccumulator {
    fn add(&mut self, r_squared: f64) -> Result<(), &'static str> {
        self.total = self.compensated.add(r_squared)?;
        self.pair_count += 1;
        Ok(())
    }
}

/// Índice del bin que contiene la distancia. Los bins del perfil son cerrados,
/// disjuntos, contiguos y ordenados, de modo que basta una búsqueda binaria.
pub(crate) fn bin_of(bins: &[LdBin], distance: u64) -> Option<usize> {
    let index = bins.partition_point(|bin| bin.distance_end_bp() < distance);
    bins.get(index)
        .filter(|bin| bin.distance_start_bp() <= distance)
        .map(|_| index)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// La columna publicada y su recálculo deben coincidir **bit a bit**.
    ///
    /// Un solo ULP de diferencia es un fallo: la puntuación de la búsqueda es
    /// incremental y solo la recomputación canónica puede publicarse. Si la
    /// comprobación desapareciera, una columna obsoleta pasaría inadvertida.
    #[test]
    fn a_published_column_must_match_its_recomputation_bit_for_bit() {
        assert!(ensure_published_column_matches(3, 0.25, 0.25).is_ok());
        assert!(ensure_published_column_matches(3, 0.0, 0.0).is_ok());
        let one_ulp_above = f64::from_bits(0.5_f64.to_bits() + 1);
        let error = ensure_published_column_matches(3, 0.5, one_ulp_above)
            .expect_err("un ULP de diferencia debe rechazarse");
        assert!(error.contains("columna 3"), "{error}");
        assert!(ensure_published_column_matches(0, 0.5, 0.4).is_err());
    }

    /// La tolerancia del perfil es inclusiva, y lo es en el borde exacto.
    ///
    /// La regla vive en una sola función, que usan el evaluador canónico y la
    /// puntuación incremental de la búsqueda. Un error igual a la tolerancia
    /// **no** deja el bin fuera; el siguiente `f64` representable sí.
    #[test]
    fn the_bin_tolerance_is_inclusive_at_the_exact_boundary() {
        let bin = LdBin::for_test(1, 4, 0.25, 0.05);
        assert_eq!(f5_bin_abs_error(0.3, &bin), 0.049_999_999_999_999_99);
        assert!(
            !f5_bin_is_outside(0.05, &bin),
            "el borde exacto quedó fuera"
        );
        assert!(!f5_bin_is_outside(0.0, &bin));
        assert!(
            f5_bin_is_outside(f64::from_bits(0.05_f64.to_bits() + 1), &bin),
            "un ULP por encima de la tolerancia no quedó fuera"
        );
    }

    /// Cuatro pesos uniformes: todo racional del oráculo cabe en `u128` sin
    /// perder un bit, porque `1/4` y sus productos son diádicos exactos.
    const UNIFORM: [f64; 4] = [0.25, 0.25, 0.25, 0.25];

    /// Cuatro pesos **no uniformes** sobre la escala `1/8`: `4/8`, `2/8`,
    /// `1/8` y `1/8`. Siguen siendo diádicos exactos, de modo que el mismo
    /// oráculo racional vale sin tolerancia alguna.
    const NON_UNIFORM: [f64; 4] = [0.5, 0.25, 0.125, 0.125];

    /// Matriz `4x2` cuyas ocho casillas enumeran el entero `mask`.
    fn matrix_from_mask(mask: u32) -> F5AlleleMatrix {
        let rows = (0..4)
            .map(|haplotype| {
                (0..2)
                    .map(|site| {
                        if (mask >> (haplotype * 2 + site)) & 1 == 1 {
                            '1'
                        } else {
                            '0'
                        }
                    })
                    .collect::<String>()
            })
            .collect::<Vec<_>>();
        let borrowed = rows.iter().map(String::as_str).collect::<Vec<_>>();
        F5AlleleMatrix::from_rows(&borrowed)
    }

    /// Pesos acumulados de la matriz sobre la escala entera común.
    fn counts(matrix: &F5AlleleMatrix, numerators: &[u128]) -> (u128, u128, u128) {
        let mut n_i = 0_u128;
        let mut n_j = 0_u128;
        let mut n_ij = 0_u128;
        for (haplotype, numerator) in numerators.iter().enumerate() {
            let alt_i = matrix.allele(haplotype, 0).expect("sitio 0");
            let alt_j = matrix.allele(haplotype, 1).expect("sitio 1");
            n_i += numerator * u128::from(alt_i);
            n_j += numerator * u128::from(alt_j);
            n_ij += numerator * u128::from(alt_i && alt_j);
        }
        (n_i, n_j, n_ij)
    }

    /// Oráculo racional exacto sobre pesos diádicos `n_h / T`.
    ///
    /// Sobre esa escala común todo es entero: `p_i=N_i/T`, `p_ij=N_ij/T`,
    /// `D=(T*N_ij-N_i*N_j)/T²` y
    /// `r²=(T*N_ij-N_i*N_j)² / [N_i(T-N_i)*N_j(T-N_j)]`, donde numerador y
    /// denominador comparten `T⁴` y se cancelan. Con `T<=8` cada cociente
    /// intermedio es un múltiplo exacto de `1/T⁴`, representable en binary64
    /// sin redondeo, de modo que la comparación puede ser una igualdad y no
    /// una tolerancia. No llama a ninguna función del evaluador.
    fn rational_oracle(
        matrix: &F5AlleleMatrix,
        numerators: &[u128],
        total: u128,
    ) -> (f64, f64, f64, f64, Option<f64>) {
        let (n_i, n_j, n_ij) = counts(matrix, numerators);
        let d_numerator = i128::try_from(total).expect("T") * i128::try_from(n_ij).expect("N_ij")
            - i128::try_from(n_i * n_j).expect("N_i*N_j");
        let denominator = n_i * (total - n_i) * n_j * (total - n_j);
        let r_squared = (denominator > 0).then(|| {
            let squared = d_numerator * d_numerator;
            squared as f64 / denominator as f64
        });
        let total = total as f64;
        (
            n_i as f64 / total,
            n_j as f64 / total,
            n_ij as f64 / total,
            d_numerator as f64 / (total * total),
            r_squared,
        )
    }

    /// Recorre las 256 matrices `4x2` posibles contra el oráculo racional y
    /// devuelve cuántas veces se vio cada categoría.
    ///
    /// Al recorrerlas todas aparecen por construcción `D` positivo, negativo y
    /// cero, `r²` 0, 1 e intermedio, y parejas `not_evaluable` —las de una
    /// columna constante, que es justo lo que produce un sitio fijado—. El
    /// llamador exige que cada categoría se haya visto: si el evaluador dejara
    /// de producir alguna, fallaría ahí.
    /// Devuelve cuántas parejas dieron `r²` exactamente cero, que solo es
    /// alcanzable cuando `T*N_ij = N_i*N_j` tiene solución en la escala.
    fn sweep_every_four_by_two_matrix(weights: &[f64; 4], numerators: &[u128], total: u128) -> u32 {
        let mut positive = 0_u32;
        let mut negative = 0_u32;
        let mut zero = 0_u32;
        let mut r_zero = 0_u32;
        let mut r_one = 0_u32;
        let mut r_middle = 0_u32;
        let mut not_evaluable = 0_u32;

        for mask in 0..256_u32 {
            let matrix = matrix_from_mask(mask);
            let pair = evaluate_f5_site_pair(weights, &matrix, 0, 1).expect("pareja evaluada");
            let (p_i, p_j, p_ij, d, r_squared) = rational_oracle(&matrix, numerators, total);

            assert_eq!(pair.p_hat_i(), p_i, "p_hat_i de la máscara {mask}");
            assert_eq!(pair.p_hat_j(), p_j, "p_hat_j de la máscara {mask}");
            assert_eq!(pair.p_hat_ij(), p_ij, "p_hat_ij de la máscara {mask}");
            assert_eq!(pair.d(), d, "D de la máscara {mask}");
            assert_eq!(pair.r_squared(), r_squared, "r² de la máscara {mask}");

            if d > 0.0 {
                positive += 1;
            } else if d < 0.0 {
                negative += 1;
            } else {
                zero += 1;
            }
            match r_squared {
                None => not_evaluable += 1,
                Some(0.0) => r_zero += 1,
                Some(1.0) => r_one += 1,
                Some(_) => r_middle += 1,
            }
        }

        assert!(
            positive > 0 && negative > 0 && zero > 0,
            "faltó un signo de D"
        );
        assert!(
            r_one > 0 && r_middle > 0,
            "faltó un valor de r²: {r_one}/{r_middle}"
        );
        assert!(not_evaluable > 0, "ninguna pareja resultó no evaluable");
        r_zero
    }

    /// Oráculo racional exacto con pesos uniformes `1/4`.
    #[test]
    fn an_exact_rational_oracle_reproduces_every_four_by_two_matrix() {
        let r_zero = sweep_every_four_by_two_matrix(&UNIFORM, &[1, 1, 1, 1], 4);
        assert!(r_zero > 0, "faltó la independencia exacta r²=0");
    }

    /// El mismo recorrido con pesos **no uniformes** `4/8, 2/8, 1/8, 1/8`.
    ///
    /// Con pesos iguales `p_hat_i` solo puede valer `k/4` y muchas
    /// combinaciones de `D` y `r²` no se alcanzan; aquí `N_i` recorre los ocho
    /// valores de la escala y el oráculo sigue siendo exacto, de modo que la
    /// prueba cubre `D` y `r²` —no solo las marginales— fuera del caso
    /// uniforme. La independencia exacta `r²=0` no es alcanzable en esta
    /// escala —exigiría `8*N_ij = N_i*N_j` con `N` en `0..=8`— y la acredita
    /// el recorrido uniforme.
    #[test]
    fn the_exact_rational_oracle_also_covers_non_uniform_weights() {
        assert_eq!(
            sweep_every_four_by_two_matrix(&NON_UNIFORM, &[4, 2, 1, 1], 8),
            0
        );
    }

    /// Una columna constante —lo que produce `p_i` fijado en 0 o en 1— nunca
    /// publica `r² = 0`: publica `not_evaluable`, incluso cuando `D` vale
    /// exactamente cero.
    #[test]
    fn a_fixed_column_is_not_evaluable_and_never_r_squared_zero() {
        for rows in [
            ["01", "01", "01", "01"],
            ["11", "11", "11", "11"],
            ["00", "00", "00", "00"],
            ["10", "11", "10", "11"],
        ] {
            let matrix = F5AlleleMatrix::from_rows(&rows);
            let pair = evaluate_f5_site_pair(&UNIFORM, &matrix, 0, 1).expect("pareja evaluada");
            assert_eq!(
                pair.d(),
                0.0,
                "la columna constante deja D en cero: {rows:?}"
            );
            assert_eq!(pair.r_squared(), None, "{rows:?} debería ser no evaluable");
            assert!(!pair.is_evaluable());
            assert_eq!(pair.satisfies_sign_rule(true), None);
        }
    }

    /// Asociación completa, repulsión completa e independencia, con el valor
    /// exacto que el oráculo racional predice.
    #[test]
    fn coupling_repulsion_and_independence_have_the_expected_d_and_r_squared() {
        let coupling = F5AlleleMatrix::from_rows(&["11", "11", "00", "00"]);
        let pair = evaluate_f5_site_pair(&UNIFORM, &coupling, 0, 1).expect("acoplamiento");
        assert_eq!(pair.d(), 0.25);
        assert_eq!(pair.r_squared(), Some(1.0));
        assert_eq!(pair.satisfies_sign_rule(true), Some(true));
        assert_eq!(pair.satisfies_sign_rule(false), Some(false));

        let repulsion = F5AlleleMatrix::from_rows(&["10", "10", "01", "01"]);
        let pair = evaluate_f5_site_pair(&UNIFORM, &repulsion, 0, 1).expect("repulsión");
        assert_eq!(pair.d(), -0.25);
        assert_eq!(pair.r_squared(), Some(1.0));
        assert_eq!(pair.satisfies_sign_rule(true), Some(false));
        assert_eq!(pair.satisfies_sign_rule(false), Some(true));

        let independent = F5AlleleMatrix::from_rows(&["11", "10", "01", "00"]);
        let pair = evaluate_f5_site_pair(&UNIFORM, &independent, 0, 1).expect("independencia");
        assert_eq!(pair.d(), 0.0);
        assert_eq!(pair.r_squared(), Some(0.0));
        assert_eq!(
            pair.satisfies_sign_rule(true),
            Some(false),
            "D=0 incumple con polaridades iguales"
        );
        assert_eq!(
            pair.satisfies_sign_rule(false),
            Some(false),
            "D=0 incumple también con polaridades distintas"
        );

        let intermediate = F5AlleleMatrix::from_rows(&["11", "10", "01", "01"]);
        let pair = evaluate_f5_site_pair(&UNIFORM, &intermediate, 0, 1).expect("intermedio");
        assert_eq!(pair.d(), -0.125);
        assert_eq!(pair.r_squared(), Some(1.0 / 3.0));
    }

    /// Pesos no diádicos: la compensación de Neumaier deja de ser cero y la
    /// prueba distingue el acumulador compensado de una suma ingenua.
    ///
    /// `0,1+0,2+0,3` ingenuo da `0.6000000000000001`; compensado da `0.6`. Con
    /// pesos `0,25` la compensación es siempre cero y esta diferencia no
    /// existiría.
    #[test]
    fn the_compensated_sum_differs_from_a_naive_one_with_non_dyadic_weights() {
        let weights = [0.1, 0.2, 0.3, 0.4];
        let naive = weights[..3].iter().fold(0.0_f64, |sum, value| sum + value);
        assert_eq!(
            naive, 0.6000000000000001,
            "la suma ingenua arrastra el error"
        );

        let matrix = F5AlleleMatrix::from_rows(&["11", "11", "11", "00"]);
        let pair = evaluate_f5_site_pair(&weights, &matrix, 0, 1).expect("pareja evaluada");
        assert_eq!(pair.p_hat_i(), 0.6, "p_hat_i compensado");
        assert_eq!(pair.p_hat_j(), 0.6, "p_hat_j compensado");
        assert_eq!(pair.p_hat_ij(), 0.6, "p_hat_ij compensado");
        assert_ne!(
            pair.p_hat_i(),
            naive,
            "la compensación cambia el último bit"
        );
    }

    /// Pesos cuya suma compensada normalizada no cae exactamente en 1.
    ///
    /// `normalize_slice` solo exige `|suma-1| <= 1e-9` **después** de dividir,
    /// así que un perfil legal y admitido puede dejar estos cuatro pesos, cuya
    /// suma de Neumaier vale `1.0000000000000002`. Una columna monomórfica
    /// —lo que produce `p_i = 1`— realiza entonces ese mismo valor.
    const OVER_UNIT: [f64; 4] = [
        0.019_251_000_000_000_004,
        0.383_455_000_000_000_05,
        0.511_446_000_000_000_1,
        0.085_848_000_000_000_01,
    ];

    /// Regresión G-1: una marginal que excede 1 por un ULP no aborta la
    /// ejecución, y su pareja tampoco se declara evaluable.
    ///
    /// El estado correcto es `not_evaluable`: DTC-089 no contempla abortar. Y
    /// la decisión debe tomarse **por factor**, porque con dos columnas
    /// monomórficas los dos factores son negativos y su producto sería
    /// positivo, de modo que mirar solo el denominador declararía evaluable
    /// una pareja sin variación y con `r²` cercano a 1.
    #[test]
    fn a_marginal_above_one_is_not_evaluable_instead_of_aborting() {
        let mut sum = NeumaierSum::default();
        let mut total = 0.0_f64;
        for weight in OVER_UNIT {
            total = sum.add(weight).expect("suma finita");
        }
        assert_eq!(total, 1.000_000_000_000_000_2, "el fixture debe exceder 1");

        let both = F5AlleleMatrix::from_rows(&["11", "11", "11", "11"]);
        let pair = evaluate_f5_site_pair(&OVER_UNIT, &both, 0, 1).expect("no puede abortar");
        assert_eq!(pair.p_hat_i(), total);
        assert_eq!(pair.p_hat_j(), total);
        assert_eq!(
            pair.r_squared(),
            None,
            "dos factores negativos no hacen evaluable la pareja"
        );

        let one_side = F5AlleleMatrix::from_rows(&["11", "11", "11", "10"]);
        let pair = evaluate_f5_site_pair(&OVER_UNIT, &one_side, 0, 1).expect("no puede abortar");
        assert!(pair.p_hat_i() > 1.0, "el factor i es negativo");
        assert!(pair.p_hat_j() < 1.0, "el factor j es positivo");
        assert_eq!(
            pair.r_squared(),
            None,
            "un solo factor negativo ya impide evaluar"
        );
    }

    /// Regresión G-2: el acoplamiento perfecto publica `r²` exactamente 1.
    ///
    /// Con dos columnas idénticas `D = fl(w - fl(w*w))` y el factor es
    /// `fl(w * fl(1-w))`: dos redondeos distintos del mismo real, cuyo
    /// cociente supera 1 en una fracción apreciable de los pesos decimales.
    #[test]
    fn perfect_coupling_never_publishes_r_squared_above_one() {
        let coupling = F5AlleleMatrix::from_rows(&["11", "00"]);
        let pair = evaluate_f5_site_pair(&[0.901, 0.099], &coupling, 0, 1).expect("acoplamiento");
        assert_eq!(pair.r_squared(), Some(1.0), "r² del acoplamiento perfecto");

        // Ningún peso decimal de tres cifras deja `r²` fuera de `[0,1]`, y el
        // recorte hace trabajo real: sin él, el cociente literal excede 1 en
        // una parte de ellos. El contrato no exige el 1 exacto —un cociente
        // ligeramente menor es correcto—, solo el dominio.
        let mut above_one = 0_u32;
        for thousandths in 1..1000_u32 {
            let weight = f64::from(thousandths) / 1000.0;
            let complement = 1.0 - weight;
            let pair = evaluate_f5_site_pair(&[weight, complement], &coupling, 0, 1)
                .expect("acoplamiento");
            let value = pair.r_squared().expect("pareja evaluable");
            assert!(
                value <= 1.0,
                "r² salió del dominio con peso {weight}: {value}"
            );
            assert!(
                (1.0 - value) < 1e-9,
                "el acoplamiento perfecto debe quedar en 1 con peso {weight}: {value}"
            );
            let factor = pair.p_hat_i() * (1.0 - pair.p_hat_i());
            if pair.d() * pair.d() / (factor * factor) > 1.0 {
                above_one += 1;
            }
        }
        assert!(
            above_one > 0,
            "el fixture debe contener casos cuyo cociente literal excede 1"
        );
    }

    /// La regla de signos es la literal de DTC-089 y no tiene cuota.
    #[test]
    fn the_sign_rule_is_strict_and_rejects_zero() {
        assert!(sign_rule_is_satisfied(1e-300, true));
        assert!(!sign_rule_is_satisfied(1e-300, false));
        assert!(sign_rule_is_satisfied(-1e-300, false));
        assert!(!sign_rule_is_satisfied(-1e-300, true));
        assert!(!sign_rule_is_satisfied(0.0, true));
        assert!(!sign_rule_is_satisfied(0.0, false));
        assert!(!sign_rule_is_satisfied(-0.0, true));
    }

    /// Errores estables del contrato interno.
    #[test]
    fn the_shape_errors_name_their_cause() {
        let matrix = F5AlleleMatrix::from_rows(&["11", "00", "10", "01"]);
        assert_eq!(
            evaluate_f5_site_pair(&UNIFORM, &matrix, 0, 2),
            Err(F5LdError::SiteOutOfRange)
        );
        assert_eq!(
            evaluate_f5_site_pair(&[0.5, 0.5], &matrix, 0, 1),
            Err(F5LdError::WeightCountMismatch)
        );
        assert!(!F5LdError::NonFiniteSum.detail().is_empty());
        assert_eq!(
            F5LdError::SiteOutOfRange.to_string(),
            F5LdError::SiteOutOfRange.detail()
        );
    }
}
