//! Búsqueda local F5.2.3: propuestas 2×2, terminación y mejor intento.
//!
//! Este módulo es el primer llamador productivo de la cadena completa de F5.2:
//! recibe el layout físico y la admisión sellada, deriva de ellos los
//! marginales, construye cada intento haplotípico, lo refina con las
//! propuestas 2×2 de DTC-089 y lo evalúa con el evaluador LD. Ni
//! `build_f5_marginals`, ni `build_f5_haplotype_attempt`, ni
//! `evaluate_f5_linkage_disequilibrium` tenían hasta ahora quien los
//! encadenara.
//!
//! # El alcance de los intentos es la unidad bloque×subpoblación
//!
//! DTC-089 ordena la aceptación con siete campos y en este orden literal:
//!
//! ```text
//! (filas_duplicadas, bins_fuera_de_tolerancia, signos_incumplidos,
//!  max_error_r2, suma_error_r2, max_error_marginal, matrix_bytes)
//! ```
//!
//! Los siete son magnitudes de **un bloque×subpoblación**: una matriz tiene
//! unas filas duplicadas y unos `matrix_bytes`, y los bins y los signos se
//! agregan dentro de un bloque. Y no solo la tupla: **los intentos tienen ese
//! mismo alcance**. El contrato escribe el coste «por bloque×subpoblación×
//! intento», `f5-attempts.tsv` se indexa por `(population_id, block_id,
//! attempt_index)` con «una sola fila `true` por bloque×población», el contexto
//! público de `GPS-DOMAIN-104` nombra población y bloque, REQ-EVO-75 dice que
//! «solo falla el bloque×subpoblación sin ninguna pareja evaluable», y
//! `attempt_summary` publica `haplotype_attempts_total` y `max_attempt_index`
//! como **escalares distintos**, cosa que solo tiene sentido si cada unidad
//! lleva su propia cuenta.
//!
//! Por eso este módulo resuelve cada unidad por separado: su presupuesto de
//! intentos, su presupuesto de propuestas, su ámbito RNG y su selección del
//! mejor intento son suyos, y para en cuanto esa unidad termina. Un intento
//! global todo-o-nada sería incorrecto además de caro: con `U` unidades que
//! terminan de forma independiente con probabilidad `p`, el intento global
//! triunfa con `p^U` y `max_haplotype_attempts` está acotado, de modo que un
//! escenario resoluble unidad a unidad se declararía agotado. Eso es justo lo
//! que DTC-089 prohíbe al decir que el agotamiento «no afirma imposibilidad».
//!
//! La independencia entre unidades no es un supuesto de conveniencia: DTC-089
//! prohíbe parejas LD interbloque, el layout las sitúa en intervalos disjuntos,
//! los ámbitos RNG `initial_column`, `polarity` y `swap_proposal` llevan ya
//! `population_id`, `block_id` y `attempt_index`, y `f5-founders.tsv` asigna un
//! haplotipo por bloque de forma separada. Ningún documento exige coherencia
//! entre bloques en esta fase.
//!
//! # Dos aritméticas, una sola verdad
//!
//! La revisión de F5.2.2 dejó escrito que el valor normativo de `p_hat_i` es
//! la suma canónica completa sobre las H filas y que el evaluador exige que la
//! columna publicada coincida con ella **bit a bit**. Mantener eso en cada
//! propuesta costaría `O(H*S)` por propuesta y rompería la cota `O(Q*S)` que
//! DTC-089 declara. La conciliación, que es la que el propio contrato
//! describe, tiene dos niveles:
//!
//! - **Puntuación.** Dentro del bucle de propuestas, marginales y conjuntas se
//!   mantienen por incrementos `O(S)`. Es una aritmética determinista pero no
//!   canónica: un `+δ` no reproduce los bits de una suma de Neumaier rehecha.
//!   Solo decide qué propuesta se acepta.
//! - **Veredicto.** Al cerrar la búsqueda de un intento se recomputan
//!   canónicamente las columnas publicadas —`refresh_realized_columns`— y se
//!   vuelve a evaluar la unidad con
//!   [`crate::f5_ld::evaluate_f5_block_linkage_disequilibrium`]. Los valores que se publican,
//!   se comparan entre intentos y deciden éxito o agotamiento son **siempre**
//!   los canónicos. Sin la recomputación previa, cada propuesta aceptada
//!   abortaría con `GPS-INTERNAL-101`.
//!
//! El estado incremental se siembra desde una evaluación canónica, de modo que
//! en la propuesta cero las dos aritméticas coinciden exactamente; la deriva
//! posterior solo se acumula sobre las propuestas **aceptadas**, porque una
//! rechazada restaura el estado anterior. Y no se deja a la confianza: al
//! cerrar cada intento, `ensure_score_did_not_diverge` contrasta las dos
//! aritméticas **pareja a pareja** y **bin a bin**.
//!
//! La igualdad bit a bit no puede exigirse: `m - q_a + q_b` no es la suma de
//! Neumaier sobre el nuevo conjunto, y hacer que lo fuera costaría `O(H)` por
//! propuesta, que es exactamente lo que DTC-089 prohíbe al fijar la
//! actualización en `O(S)` y lo que la fórmula de admisión contabiliza. Lo que
//! sí se exige es que **ninguna** divergencia quede sin explicar: las tres
//! reglas discretas del contrato son discontinuas en un punto exacto —la de
//! signos en `D == 0.0`, la de evaluabilidad en `factor == 0.0`, la de bins en
//! `error == max_abs_error`—, y una divergencia se exime solo si uno de sus dos
//! valores está justo ahí. Cualquier otra, y cualquier diferencia en
//! `filas_duplicadas`, que no cuelga de ninguna frontera de coma flotante, es
//! `GPS-INTERNAL-101`: la búsqueda no puede explorar un problema que no es el
//! suyo.
//!
//! # Validez antes que tupla
//!
//! Un bloque sin ninguna pareja evaluable publica el mejor valor posible en
//! cuatro de los siete campos —`bins_outside=0`, `signs_failed=0`,
//! `max_r2_error=0` y `sum_r2_error=0`—, y en cinco si además no repite
//! ninguna fila, de modo que un intento degenerado ganaría la selección y
//! podría incluso parecer terminado. Por eso
//! `is_valid()` se consulta **antes** que la tupla, tanto para terminar como
//! para seleccionar, y la clave de orden antepone dos campos a los siete del
//! contrato: si el bloque terminó y si es válido.
//!
//! # Coste medido
//!
//! Por bloque×subpoblación×intento, y con `B` bins del perfil:
//!
//! - **Siembra canónica:** `H*S` visitas de fila para las S marginales más
//!   `H*S*(S-1)/2` para las conjuntas, es decir `O(H*S²)`. El registro publica
//!   ese recuento en `row_visits`, que vale exactamente `H*(S + S*(S-1)/2)`
//!   **sea cual sea el presupuesto de propuestas**.
//! - **Cada propuesta:** cuatro enteros acotados y, si es candidata, dos
//!   marginales movidas en `O(1)`, `2S-3` parejas actualizadas, dos filas
//!   rehasheadas —`O(S/8)` bytes— y un recorrido de los `B` bins. **Ninguna
//!   recorre las H filas.** El registro publica `pair_updates`, que vale
//!   exactamente `candidate_proposals*(2S-3)` y se cuenta dentro de
//!   `refresh_pair`.
//!
//!   `row_visits` no es una identidad escrita a mano: el contador vive dentro
//!   de `site_marginal` y `joint_marginal`, de modo que **cualquier** suma
//!   canónica lo incrementa, la haga la siembra o la haga una propuesta. Que
//!   valga `H*(S + S*(S-1)/2)` para cualquier presupuesto es, por tanto, una
//!   medida y no una tautología.
//! - **Cierre:** `O(H*S)` para recomputar canónicamente las columnas,
//!   `O(H*S²)` para la evaluación LD de la unidad y `O(H*S/8)` para el
//!   desempate lexicográfico de `matrix_bytes` contra el mejor candidato.
//!
//! El total es `O(H*S² + Q*(S+B))`, con `Q = max_swap_proposals`. DTC-089
//! escribe `O(H*S² + Q*S)` porque `B` es una propiedad del perfil y no de la
//! búsqueda; el término se declara aquí en lugar de esconderlo.
//!
//! La memoria añadida es la caché de parejas: 32 bytes por pareja intrabloque,
//! más los veredictos canónicos que la guarda de cierre recoge.
//! `f5-resource-accounting-v7` presupuesta `192*S*(S-1)/2` bytes por
//! bloque×subpoblación para ese fin, de modo que la caché cabe con holgura en
//! lo ya admitido y no reabre la contabilidad sellada.
//!
//! # Imposibilidad demostrada frente a búsqueda agotada
//!
//! `GPS-RES-101` por un fallo exacto de dimensiones —`H*S` por encima de
//! `max_haplotype_cells`, o un lado por encima de 4096— no se comprueba aquí:
//! la admisión sellada, que esta frontera exige y verifica, ya lo rechaza con
//! ese mismo código sobre la cota `H*site_bound_per_block`, que domina el
//! producto exacto del layout realizado. Repetirlo en la búsqueda añadiría una
//! rama inalcanzable y no una red.
//!
//! `GPS-DOMAIN-102` se reserva para una **demostración** hecha sobre los
//! marginales y antes de abrir ningún intento: un bloque con menos de dos
//! columnas cuyo objetivo no sea 0 o 1 no puede producir ninguna pareja
//! evaluable, porque una columna fijada es constante y su factor
//! `p(1-p)` nunca es positivo. Agotar `max_haplotype_attempts` devuelve en
//! cambio `GPS-DOMAIN-104`, conserva el mejor intento y **no afirma
//! imposibilidad**. Ningún fallo reduce H o S, cambia el perfil o una
//! tolerancia, fusiona bins ni acepta menos bloques.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use sha2::{Digest, Sha256};

use crate::application::EffectiveF5Config;
use crate::diagnostics::Diagnostic;
use crate::f5_haplotypes::{
    F5AlleleMatrix, F5BlockHaplotypes, build_f5_block_attempt, ensure_f5_haplotype_cardinality,
    is_fixed_target,
};
use crate::f5_layout::{F5LayoutBlock, F5PhysicalLayout};
use crate::f5_ld::{
    F5BlockLdEvaluation, F5CanonicalPairVerdict, bin_of, evaluate_f5_block_with_pair_verdicts,
    f5_bin_abs_error, f5_bin_is_outside, f5_evaluability_factor, joint_marginal,
    pair_from_components, sign_rule_is_satisfied, site_marginal,
};
use crate::f5_marginals::{F5MarginalTable, NeumaierSum, build_f5_marginals};
use crate::f5_preflight::{
    F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS, F5_DEFAULT_MAX_SWAP_PROPOSALS, F5OtherResources,
    SealedF5Admission,
};
use crate::f5_profile::LdBin;
use crate::identity::{ExecutionId, F5BlockId, PopulationId, ScenarioId};
use crate::randomness::{
    F5RandomDecision, LocalRng, RandomWordSource, RootSeed, derive_f5_subseed,
    f5_uniform_index_with_source,
};

/// Bin sin asignar dentro de la caché de parejas.
const NO_BIN: u32 = u32::MAX;

// ---------------------------------------------------------------------------
// Resultado publicable
// ---------------------------------------------------------------------------

/// Una fila de `f5-attempts.tsv`: un bloque×subpoblación dentro de un intento.
///
/// Los doce campos de la cabecera normativa están aquí con su nombre. Junto a
/// ellos viajan cuatro magnitudes que **ningún artefacto publica** y que
/// existen para que la batería pueda medir el comportamiento en lugar de
/// suponerlo: `evaluable_pairs` y `terminated`, que explican el veredicto, y
/// `candidate_proposals`, `row_visits` y `pair_updates`, que miden el coste.
/// `selected_best` es `true` exactamente una vez por bloque×subpoblación.
#[derive(Clone, Debug, PartialEq)]
pub struct F5AttemptRecord {
    population_id: PopulationId,
    block_index: u64,
    block_id: F5BlockId,
    attempt_index: u64,
    proposals: u64,
    accepted: u64,
    duplicate_rows: u64,
    bins_outside: u64,
    signs_failed: u64,
    max_r2_error: f64,
    sum_r2_error: f64,
    max_marginal_error: f64,
    evaluable_pairs: u64,
    terminated: bool,
    matrix_bytes_sha256: String,
    selected_best: bool,
    candidate_proposals: u64,
    row_visits: u64,
    pair_updates: u64,
}

impl F5AttemptRecord {
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

    #[must_use]
    pub const fn attempt_index(&self) -> u64 {
        self.attempt_index
    }

    /// Propuestas consumidas del presupuesto `max_swap_proposals`.
    #[must_use]
    pub const fn proposals(&self) -> u64 {
        self.proposals
    }

    /// Propuestas que mejoraron estrictamente la tupla. Un empate no cuenta.
    #[must_use]
    pub const fn accepted(&self) -> u64 {
        self.accepted
    }

    #[must_use]
    pub const fn duplicate_rows(&self) -> u64 {
        self.duplicate_rows
    }

    #[must_use]
    pub const fn bins_outside(&self) -> u64 {
        self.bins_outside
    }

    #[must_use]
    pub const fn signs_failed(&self) -> u64 {
        self.signs_failed
    }

    #[must_use]
    pub const fn max_r2_error(&self) -> f64 {
        self.max_r2_error
    }

    #[must_use]
    pub const fn sum_r2_error(&self) -> f64 {
        self.sum_r2_error
    }

    #[must_use]
    pub const fn max_marginal_error(&self) -> f64 {
        self.max_marginal_error
    }

    /// Parejas intrabloque con denominador no nulo. Cero invalida el bloque.
    #[must_use]
    pub const fn evaluable_pairs(&self) -> u64 {
        self.evaluable_pairs
    }

    /// El bloque alcanzó la terminación de DTC-089: es válido, tiene cero en
    /// los tres campos de error discreto y respeta todas las tolerancias.
    #[must_use]
    pub const fn terminated(&self) -> bool {
        self.terminated
    }

    /// SHA-256 de la codificación canónica `matrix_bytes`, que el esquema de
    /// fallo exige. El desempate **no** lo usa: compara los bytes.
    #[must_use]
    pub fn matrix_bytes_sha256(&self) -> &str {
        &self.matrix_bytes_sha256
    }

    #[must_use]
    pub const fn selected_best(&self) -> bool {
        self.selected_best
    }

    /// Propuestas cuya configuración `10/01 ↔ 01/10` sí era candidata.
    #[must_use]
    pub const fn candidate_proposals(&self) -> u64 {
        self.candidate_proposals
    }

    /// Visitas de fila que la búsqueda pidió a las sumas canónicas de `f5_ld`.
    ///
    /// El contador vive dentro de `site_marginal` y `joint_marginal`, de modo
    /// que **toda** visita lo incrementa venga de donde venga. Vale
    /// exactamente `H*S + H*S*(S-1)/2` por intento: si una propuesta hiciera
    /// una suma canónica, el recuento crecería con el presupuesto.
    #[must_use]
    pub const fn row_visits(&self) -> u64 {
        self.row_visits
    }

    /// Parejas rehechas por incrementos, contadas dentro de `refresh_pair`.
    /// Vale `candidate_proposals*(2S-3)`.
    #[must_use]
    pub const fn pair_updates(&self) -> u64 {
        self.pair_updates
    }
}

/// Búsqueda con éxito: el bloque ganador de cada unidad, su evaluación
/// canónica y el historial completo de intentos.
///
/// No hay un «intento publicado» único: DTC-089 cuenta los intentos por
/// bloque×subpoblación y cada unidad para en cuanto termina, de modo que el
/// resultado compone los ganadores de cada una. `attempt_summary` del
/// manifiesto pide exactamente eso —`haplotype_attempts_total` y
/// `max_attempt_index` son escalares distintos— y aquí se publican con ese
/// nombre.
#[derive(Clone, Debug, PartialEq)]
pub struct F5HaplotypeSearch {
    blocks: Vec<F5BlockHaplotypes>,
    ld_evaluation: Vec<F5BlockLdEvaluation>,
    records: Vec<F5AttemptRecord>,
    haplotype_attempts_total: u64,
    max_attempt_index: u64,
}

impl F5HaplotypeSearch {
    /// Bloques ganadores en el orden canónico de la tabla de marginales.
    #[must_use]
    pub fn blocks(&self) -> &[F5BlockHaplotypes] {
        &self.blocks
    }

    #[must_use]
    pub fn block(
        &self,
        population_id: &PopulationId,
        block_index: u64,
    ) -> Option<&F5BlockHaplotypes> {
        self.blocks
            .iter()
            .find(|row| row.population_id() == population_id && row.block_index() == block_index)
    }

    /// Evaluación LD canónica de cada bloque ganador, en el mismo orden.
    #[must_use]
    pub fn ld_evaluation(&self) -> &[F5BlockLdEvaluation] {
        &self.ld_evaluation
    }

    /// `haplotype_attempts_total` de `attempt_summary`: suma de los intentos
    /// ejecutados por todas las unidades.
    #[must_use]
    pub const fn haplotype_attempts_total(&self) -> u64 {
        self.haplotype_attempts_total
    }

    /// `max_attempt_index` de `attempt_summary`: el mayor índice de intento
    /// que alguna unidad necesitó.
    #[must_use]
    pub const fn max_attempt_index(&self) -> u64 {
        self.max_attempt_index
    }

    /// Historial en orden canónico de unidad y, dentro de ella, de intento.
    #[must_use]
    pub fn records(&self) -> &[F5AttemptRecord] {
        &self.records
    }
}

/// Fallo de la búsqueda que conserva el historial de intentos.
///
/// El agotamiento no afirma imposibilidad: el mejor intento de cada
/// bloque×subpoblación recorrido queda marcado en `records` con
/// `selected_best`.
#[derive(Clone, Debug, PartialEq)]
pub struct F5SearchFailure {
    diagnostic: Box<Diagnostic>,
    records: Vec<F5AttemptRecord>,
}

impl F5SearchFailure {
    #[must_use]
    pub const fn diagnostic(&self) -> &Diagnostic {
        &self.diagnostic
    }

    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.diagnostic.code()
    }

    #[must_use]
    pub fn context(&self) -> &BTreeMap<String, String> {
        self.diagnostic.context()
    }

    #[must_use]
    pub fn records(&self) -> &[F5AttemptRecord] {
        &self.records
    }

    #[must_use]
    pub fn into_diagnostic(self) -> Diagnostic {
        *self.diagnostic
    }
}

impl From<Diagnostic> for F5SearchFailure {
    fn from(diagnostic: Diagnostic) -> Self {
        Self {
            diagnostic: Box::new(diagnostic),
            records: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Frontera pública
// ---------------------------------------------------------------------------

/// Ejecuta la búsqueda haplotípica completa de DTC-089.
///
/// Encadena marginales, construcción, propuestas 2×2 y evaluación LD. Exige la
/// admisión sellada de la misma configuración porque abre ámbitos RNG, y
/// recibe el layout ya materializado porque las distancias de los bins son
/// físicas.
///
/// Recorre las unidades bloque×subpoblación en el orden canónico de la tabla
/// de marginales y resuelve cada una por separado, con su propio presupuesto
/// de intentos. La primera que agota `max_haplotype_attempts` detiene la
/// búsqueda: el resultado ya no puede publicarse y gastar el presupuesto de
/// las unidades restantes no aportaría diagnóstico.
///
/// # Errors
///
/// - `GPS-DOMAIN-102` si algún bloque demuestra, antes de abrir intentos, que
///   no puede existir ninguna pareja LD evaluable.
/// - `GPS-DOMAIN-104` si una unidad agota `max_haplotype_attempts` sin un
///   intento aceptable.
/// - Los diagnósticos de las etapas encadenadas: `GPS-CONF-101` por piezas
///   ajenas, `GPS-DOMAIN-105` por cardinalidad y `GPS-INTERNAL-101` por una
///   aritmética indefinida o por una divergencia entre la puntuación
///   incremental y la evaluación canónica.
pub fn search_f5_haplotypes(
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
    execution_id: &ExecutionId,
    layout: &F5PhysicalLayout,
) -> Result<F5HaplotypeSearch, F5SearchFailure> {
    if !sealed.matches_admission(config, F5OtherResources::default()) {
        return Err(Diagnostic::f5_configuration(
            config.scenario_path().display().to_string(),
            "preflight.binding",
            "la admisión sellada no corresponde a esta configuración",
        )
        .into());
    }
    let scenario_id = ScenarioId::parse(config.scenario_id(), "scenario_id").map_err(|_| {
        Diagnostic::f5_configuration(
            config.scenario_path().display().to_string(),
            "scenario_id",
            "no se pudo reconstruir como identidad tipada",
        )
    })?;
    let root_seed = RootSeed::parse("seed", config.seed()).map_err(|_| {
        Diagnostic::f5_configuration(
            config.scenario_path().display().to_string(),
            "seed",
            "la semilla efectiva no es canónica",
        )
    })?;

    let marginals = build_f5_marginals(config, sealed, execution_id, layout)?;
    let weights = config.profile().haplotype_weights();

    // Las dos comprobaciones previas de DTC-089, antes de abrir ningún intento
    // y, por tanto, antes de fundadores. El orden importa: un bloque sin
    // columnas variables incumple también `H <= 2^V`, y de los dos
    // diagnósticos el que describe la causa es `GPS-DOMAIN-102` —no existe
    // ninguna pareja evaluable—, no `GPS-DOMAIN-105`, que hablaría de
    // cardinalidad cuando el problema es que no hay variación que medir.
    ensure_some_pair_can_be_evaluable(&scenario_id, &marginals)?;
    ensure_f5_haplotype_cardinality(config, &marginals)?;

    let scenario = config.declared_scenario();
    let max_attempts = scenario
        .max_haplotype_attempts
        .as_ref()
        .map_or(F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS, |value| *value.value());
    let max_proposals = scenario
        .max_swap_proposals
        .as_ref()
        .map_or(F5_DEFAULT_MAX_SWAP_PROPOSALS, |value| *value.value());

    let layout_blocks = index_layout_blocks(config, layout)?;
    let mut records: Vec<F5AttemptRecord> = Vec::new();
    let mut blocks = Vec::with_capacity(marginals.blocks().len());
    let mut ld_evaluation = Vec::with_capacity(marginals.blocks().len());
    let mut haplotype_attempts_total = 0_u64;
    let mut max_attempt_index = 0_u64;

    for (block_ordinal, marginal_block) in marginals.blocks().iter().enumerate() {
        let layout_block = layout_blocks
            .get(&marginal_block.block_index())
            .ok_or_else(|| {
                Diagnostic::f5_configuration(
                    config.scenario_path().display().to_string(),
                    "layout.identity",
                    "los marginales nombran un bloque que este layout no contiene",
                )
            })?;
        let unit = search_block(
            config,
            &scenario_id,
            root_seed,
            weights,
            &marginals,
            block_ordinal,
            layout_block,
            max_attempts,
            max_proposals,
        )
        .map_err(|diagnostic| F5SearchFailure {
            diagnostic: Box::new(diagnostic),
            records: records.clone(),
        })?;

        haplotype_attempts_total += unit.attempts_executed;
        max_attempt_index = max_attempt_index.max(unit.best_attempt_index);
        records.extend(unit.records);
        let Some((block, evaluation)) = unit.terminated else {
            return Err(exhausted_failure(&scenario_id, records, max_attempts));
        };
        blocks.push(block);
        ld_evaluation.push(evaluation);
    }

    Ok(F5HaplotypeSearch {
        blocks,
        ld_evaluation,
        records,
        haplotype_attempts_total,
        max_attempt_index,
    })
}

// ---------------------------------------------------------------------------
// Comprobaciones previas a cualquier intento
// ---------------------------------------------------------------------------

/// `GPS-DOMAIN-102` por ausencia **demostrada** de toda pareja LD evaluable.
///
/// Una columna cuyo objetivo es 0 o 1 es constante por construcción, de modo
/// que su `p_hat` vale 0 o la suma de los pesos y su factor `p(1-p)` nunca es
/// positivo: toda pareja que la contenga es `not_evaluable`. Con menos de dos
/// columnas no fijadas no existe ninguna pareja que pueda evaluarse, y eso no
/// depende del intento ni de la semilla. La demostración se hace sobre los
/// marginales, antes de abrir el primer intento y antes de fundadores.
///
/// Alcance: basta **un** bloque así para detener el escenario. REQ-EVO-75 lo
/// dice literalmente —«solo falla el bloque×subpoblación sin ninguna pareja
/// evaluable»— y el contexto público del diagnóstico nombra la subpoblación
/// del bloque que lo demuestra, no el escenario entero.
fn ensure_some_pair_can_be_evaluable(
    scenario_id: &ScenarioId,
    marginals: &F5MarginalTable,
) -> Result<(), Diagnostic> {
    for block in marginals.blocks() {
        let variable = block
            .sites()
            .iter()
            .filter(|site| !is_fixed_target(site.target_alt_af()))
            .count();
        // `variable < 2` cubre también un bloque de un solo sitio: el layout
        // ya descarta menos de dos SNP, de modo que ese caso no llega aquí.
        if variable < 2 {
            return Err(Diagnostic::f5_no_evaluable_ld_pair(
                scenario_id.as_str(),
                block.population_id().as_str(),
                format!(
                    "el bloque {} tiene {} sitios y solo {variable} columnas no fijadas: ninguna pareja intrabloque puede tener denominador no nulo",
                    block.block_id().as_str(),
                    block.sites().len()
                ),
            ));
        }
    }
    Ok(())
}

fn index_layout_blocks<'a>(
    config: &EffectiveF5Config,
    layout: &'a F5PhysicalLayout,
) -> Result<BTreeMap<u64, &'a F5LayoutBlock>, Diagnostic> {
    let mut blocks = BTreeMap::new();
    for block in layout.blocks() {
        if blocks.insert(block.block_index(), block).is_some() {
            return Err(Diagnostic::f5_configuration(
                config.scenario_path().display().to_string(),
                "layout.identity",
                "el layout repite un índice de bloque",
            ));
        }
    }
    Ok(blocks)
}

fn exhausted_failure(
    scenario_id: &ScenarioId,
    records: Vec<F5AttemptRecord>,
    max_attempts: u64,
) -> F5SearchFailure {
    let offender = records
        .iter()
        .find(|record| record.selected_best && !record.terminated);
    let diagnostic = match offender {
        Some(record) => Diagnostic::f5_haplotype_search_exhausted(
            scenario_id.as_str(),
            record.population_id.as_str(),
            record.block_id.as_str(),
            format!(
                "best_attempt_index={};attempts={max_attempts};duplicate_rows={};bins_outside={};signs_failed={};max_r2_error={};sum_r2_error={};max_marginal_error={};evaluable_pairs={}",
                record.attempt_index,
                record.duplicate_rows,
                record.bins_outside,
                record.signs_failed,
                record.max_r2_error,
                record.sum_r2_error,
                record.max_marginal_error,
                record.evaluable_pairs,
            ),
        ),
        None => Diagnostic::f5_haplotype_search_exhausted(
            scenario_id.as_str(),
            "<none>",
            "<none>",
            format!("max_haplotype_attempts={max_attempts} no permitió ejecutar ningún intento"),
        ),
    };
    F5SearchFailure {
        diagnostic: Box::new(diagnostic),
        records,
    }
}

// ---------------------------------------------------------------------------
// Búsqueda de una unidad bloque×subpoblación
// ---------------------------------------------------------------------------

/// Candidato de una unidad: su matriz, su evaluación canónica y la clave de
/// orden con la que compite contra los demás intentos de la misma unidad.
struct UnitCandidate {
    block: F5BlockHaplotypes,
    evaluation: F5BlockLdEvaluation,
    score: ScoreTuple,
    terminated: bool,
    attempt_index: u64,
}

/// Clave de orden de un candidato dentro de su unidad.
///
/// Es la tupla de DTC-089 con las dos claves que la revisión de F5.2.2 exigió
/// anteponer y con el `attempt_index` detrás. Se expone como valor propio para
/// que la regla pueda probarse sin construir un bloque entero.
struct UnitOrderKey<'a> {
    terminated: bool,
    valid: bool,
    score: &'a ScoreTuple,
    matrix: &'a F5AlleleMatrix,
    attempt_index: u64,
}

impl UnitCandidate {
    fn order_key(&self) -> UnitOrderKey<'_> {
        UnitOrderKey {
            terminated: self.terminated,
            valid: self.evaluation.is_valid(),
            score: &self.score,
            matrix: self.block.matrix(),
            attempt_index: self.attempt_index,
        }
    }
}

/// Orden de DTC-089 con las dos claves que la revisión de F5.2.2 exigió
/// anteponer.
///
/// La tupla del contrato premiaría a un intento degenerado —sin ninguna pareja
/// evaluable tiene su mejor valor posible en cuatro de sus campos—, así que
/// antes van la terminación y la validez del bloque. Después llegan los siete
/// campos en su orden literal y, como último desempate, el menor
/// `attempt_index`.
///
/// El séptimo campo es `matrix_bytes` y se compara **byte a byte** con
/// [`F5AlleleMatrix::canonical_cmp`]. No se compara su SHA-256: el contrato
/// pide el mínimo *lexicográfico* de la codificación, y un digest es una
/// permutación pseudoaleatoria que invierte ese orden en la mitad de las
/// parejas. El digest sigue publicándose en el registro porque el esquema de
/// fallo lo exige, pero no decide nada.
fn order_key_is_better(candidate: &UnitOrderKey<'_>, current: &UnitOrderKey<'_>) -> bool {
    let leading = (
        !candidate.terminated,
        !candidate.valid,
        candidate.score.duplicate_rows,
        candidate.score.bins_outside,
        candidate.score.signs_failed,
    );
    let current_leading = (
        !current.terminated,
        !current.valid,
        current.score.duplicate_rows,
        current.score.bins_outside,
        current.score.signs_failed,
    );
    if leading != current_leading {
        return leading < current_leading;
    }
    for (left, right) in [
        (candidate.score.max_r2_error, current.score.max_r2_error),
        (candidate.score.sum_r2_error, current.score.sum_r2_error),
        (
            candidate.score.max_marginal_error,
            current.score.max_marginal_error,
        ),
    ] {
        let ordering = left.total_cmp(&right);
        if ordering != Ordering::Equal {
            return ordering == Ordering::Less;
        }
    }
    match candidate.matrix.canonical_cmp(current.matrix) {
        Ordering::Equal => candidate.attempt_index < current.attempt_index,
        other => other == Ordering::Less,
    }
}

/// Resultado de resolver una unidad bloque×subpoblación.
struct UnitOutcome {
    records: Vec<F5AttemptRecord>,
    /// Bloque y evaluación ganadores, o `None` si la unidad agotó su límite.
    terminated: Option<(F5BlockHaplotypes, F5BlockLdEvaluation)>,
    attempts_executed: u64,
    best_attempt_index: u64,
}

/// Divergencia máxima tolerada entre la puntuación incremental y la evaluación
/// canónica en los tres campos `f64`.
///
/// Los cuatro enteros se exigen idénticos: un incremento mal aplicado los mueve
/// y ninguna aritmética legítima lo hace. Los `f64` sí difieren por unos pocos
/// ULP —la puntuación mantiene los totales de bin por incrementos y el canon
/// los rehace con Neumaier—, y la cota separa ese ruido de un defecto.
const SCORE_DIVERGENCE_TOLERANCE: f64 = 1e-9;

#[expect(
    clippy::too_many_arguments,
    reason = "la unidad necesita identidad, semilla, perfil, marginales, layout y los dos límites; agruparlos en un tipo intermedio no añadiría ninguna invariante"
)]
fn search_block(
    config: &EffectiveF5Config,
    scenario_id: &ScenarioId,
    root_seed: RootSeed,
    weights: &[f64],
    marginals: &F5MarginalTable,
    block_ordinal: usize,
    layout_block: &F5LayoutBlock,
    max_attempts: u64,
    max_proposals: u64,
) -> Result<UnitOutcome, Diagnostic> {
    let mut records = Vec::new();
    let mut best: Option<UnitCandidate> = None;
    let mut attempts_executed = 0_u64;

    for attempt_index in 0..max_attempts {
        attempts_executed += 1;
        let mut block = build_f5_block_attempt(config, marginals, block_ordinal, attempt_index)?;
        let (measure, state) = refine_block(
            scenario_id,
            root_seed,
            weights,
            config
                .profile()
                .ld_bins(block.population_id().as_str())
                .ok_or_else(|| {
                    Diagnostic::f5_configuration(
                        config.scenario_path().display().to_string(),
                        "profile.ld_profile",
                        "el perfil no declara bins LD para la subpoblación del intento",
                    )
                })?,
            layout_block,
            &mut block,
            attempt_index,
            max_proposals,
        )?;
        let mut verdicts = Vec::new();
        let evaluation = evaluate_f5_block_with_pair_verdicts(
            config,
            layout_block,
            &block,
            Some(&mut verdicts),
        )?;
        let score = ScoreTuple {
            duplicate_rows: block.duplicate_rows(),
            bins_outside: evaluation.bins_outside(),
            signs_failed: evaluation.signs_failed(),
            max_r2_error: evaluation.max_r2_error(),
            sum_r2_error: evaluation.sum_r2_error(),
            max_marginal_error: block.max_abs_error(),
        };
        ensure_score_did_not_diverge(scenario_id, &state, &verdicts, &score, &evaluation)?;
        let terminated = block_has_terminated(
            &score,
            evaluation.is_valid(),
            block.marginal_max_abs_error(),
        );

        records.push(F5AttemptRecord {
            population_id: block.population_id().clone(),
            block_index: block.block_index(),
            block_id: block.block_id().clone(),
            attempt_index,
            proposals: measure.proposals,
            accepted: measure.accepted,
            duplicate_rows: score.duplicate_rows,
            bins_outside: score.bins_outside,
            signs_failed: score.signs_failed,
            max_r2_error: score.max_r2_error,
            sum_r2_error: score.sum_r2_error,
            max_marginal_error: score.max_marginal_error,
            evaluable_pairs: evaluation.evaluable_pairs(),
            terminated,
            matrix_bytes_sha256: format!("{:x}", Sha256::digest(block.matrix().canonical_bytes())),
            selected_best: false,
            candidate_proposals: measure.candidate_proposals,
            row_visits: measure.row_visits,
            pair_updates: measure.pair_updates,
        });

        let candidate = UnitCandidate {
            block,
            evaluation,
            score,
            terminated,
            attempt_index,
        };
        if best
            .as_ref()
            .is_none_or(|current| order_key_is_better(&candidate.order_key(), &current.order_key()))
        {
            best = Some(candidate);
        }
        if terminated {
            break;
        }
    }

    // `selected_best` señala el mínimo de la clave de orden, que con el bucle
    // por unidad es también la matriz publicada: la terminación es la primera
    // clave y el bucle para en el primer intento que termina, de modo que un
    // intento terminado siempre es el mejor de su unidad.
    let Some(best) = best else {
        return Ok(UnitOutcome {
            records,
            terminated: None,
            attempts_executed,
            best_attempt_index: 0,
        });
    };
    for record in &mut records {
        record.selected_best = record.attempt_index == best.attempt_index;
    }
    Ok(UnitOutcome {
        records,
        best_attempt_index: best.attempt_index,
        terminated: best.terminated.then_some((best.block, best.evaluation)),
        attempts_executed,
    })
}

/// Exige que cada divergencia entre la puntuación incremental y el canon esté
/// **atribuida** a una frontera exacta del contrato.
///
/// El bucle decide con aritmética incremental y el veredicto se recalcula
/// canónicamente. Las dos no coinciden bit a bit —`m - q_a + q_b` no es la suma
/// de Neumaier sobre el nuevo conjunto— y el contrato no permite hacerlas
/// coincidir: DTC-089 fija la actualización de una propuesta en `O(S)` y la
/// fórmula de admisión contabiliza ese coste exacto.
///
/// Lo que sí puede exigirse es que **ninguna** divergencia quede sin explicar.
/// Las tres reglas discretas del contrato son discontinuas en un punto exacto:
///
/// - la de signos, en `D == 0.0` («D=0 incumple»);
/// - la de evaluabilidad, en `factor == 0.0` (`p` exactamente 0 o 1);
/// - la de bins, en `error == max_abs_error` (tolerancia inclusiva).
///
/// Una pareja o un bin que difieran **teniendo** uno de sus valores
/// exactamente en la frontera son ruido de un ULP cayendo a un lado o al otro.
/// Cualquier otra divergencia es `GPS-INTERNAL-101`. `duplicate_rows` no cuelga
/// de ninguna frontera de coma flotante: ahí la igualdad es exacta siempre.
fn ensure_score_did_not_diverge(
    scenario_id: &ScenarioId,
    state: &BlockSearchState<'_>,
    verdicts: &[F5CanonicalPairVerdict],
    canonical: &ScoreTuple,
    evaluation: &F5BlockLdEvaluation,
) -> Result<(), Diagnostic> {
    let internal = |detail: String| Diagnostic::f5_internal(scenario_id.as_str(), detail);
    if state.current.duplicate_rows != canonical.duplicate_rows {
        return Err(internal(format!(
            "la puntuación incremental cuenta {} filas duplicadas y la canónica {}: ese recuento no cuelga de ninguna frontera de coma flotante",
            state.current.duplicate_rows, canonical.duplicate_rows
        )));
    }
    if verdicts.len() != state.pairs.len() {
        return Err(internal(format!(
            "el canon publicó {} veredictos de pareja y el estado incremental tiene {}",
            verdicts.len(),
            state.pairs.len()
        )));
    }

    for i in 0..state.site_count {
        for j in (i + 1)..state.site_count {
            let index = pair_index(i, j, state.site_count);
            let incremental = state.pairs[index];
            let canonical_pair = verdicts[index];
            if incremental.evaluable != canonical_pair.evaluable {
                // Frontera de evaluabilidad: `factor > 0.0`.
                let at_the_boundary = canonical_pair.factor_i == 0.0
                    || canonical_pair.factor_j == 0.0
                    || f5_evaluability_factor(state.marginals[i]) == 0.0
                    || f5_evaluability_factor(state.marginals[j]) == 0.0;
                if !at_the_boundary {
                    return Err(internal(format!(
                        "la pareja ({i},{j}) es evaluable para la puntuación incremental ({}) y no para el canon ({}), sin ningún factor exactamente nulo que lo explique",
                        incremental.evaluable, canonical_pair.evaluable
                    )));
                }
                continue;
            }
            if incremental.evaluable
                && incremental.sign_ok != canonical_pair.sign_ok
                // Frontera de signos: `D == 0.0` incumple.
                && incremental.d != 0.0
                && canonical_pair.d != 0.0
            {
                return Err(internal(format!(
                    "la pareja ({i},{j}) cumple la regla de signos para la puntuación incremental ({}) y no para el canon ({}), con D incremental {} y canónica {}: ninguna de las dos es exactamente cero",
                    incremental.sign_ok, canonical_pair.sign_ok, incremental.d, canonical_pair.d
                )));
            }
        }
    }

    for (index, (definition, evaluated)) in state.bins.iter().zip(evaluation.bins()).enumerate() {
        let incremental = (state.bin_count[index] > 0).then(|| {
            f5_bin_abs_error(
                state.bin_total[index] / state.bin_count[index] as f64,
                definition,
            )
        });
        let (Some(incremental_error), Some(canonical_error)) = (incremental, evaluated.abs_error())
        else {
            if incremental.is_some() != evaluated.abs_error().is_some() {
                return Err(internal(format!(
                    "el bin {index} es aplicable para una aritmética y no para la otra"
                )));
            }
            continue;
        };
        let incremental_outside = f5_bin_is_outside(incremental_error, definition);
        let canonical_outside = f5_bin_is_outside(canonical_error, definition);
        if incremental_outside != canonical_outside
            // Frontera de bin: el error justo en la tolerancia inclusiva.
            && incremental_error != definition.max_abs_error()
            && canonical_error != definition.max_abs_error()
        {
            return Err(internal(format!(
                "el bin {index} queda fuera de tolerancia para una aritmética ({incremental_outside}) y dentro para la otra ({canonical_outside}), con errores {incremental_error} y {canonical_error} frente a {}",
                definition.max_abs_error()
            )));
        }
    }

    // Los tres campos `f64` no tienen frontera discreta, pero sí una cota de
    // ruido: la puntuación los mantiene por incrementos y el canon los rehace
    // con Neumaier, de modo que unos pocos ULP son legítimos y una diferencia
    // apreciable no lo es.
    for (field, left, right) in [
        (
            "max_r2_error",
            state.current.max_r2_error,
            canonical.max_r2_error,
        ),
        (
            "sum_r2_error",
            state.current.sum_r2_error,
            canonical.sum_r2_error,
        ),
        (
            "max_marginal_error",
            state.current.max_marginal_error,
            canonical.max_marginal_error,
        ),
    ] {
        // Un `NaN` también debe fallar, de modo que la condición se escribe
        // sobre el orden total y no como una comparación negada.
        if (left - right).abs().total_cmp(&SCORE_DIVERGENCE_TOLERANCE) == Ordering::Greater {
            return Err(internal(format!(
                "la puntuación incremental de {field} vale {left} y la canónica {right}: la deriva supera {SCORE_DIVERGENCE_TOLERANCE}"
            )));
        }
    }
    Ok(())
}

/// Coste, consumo y puntuación final de un intento de una unidad.
#[derive(Clone, Copy, Debug)]
struct BlockSearchMeasure {
    proposals: u64,
    candidate_proposals: u64,
    accepted: u64,
    row_visits: u64,
    pair_updates: u64,
}

#[expect(
    clippy::too_many_arguments,
    reason = "el bucle de propuestas necesita identidad, semilla, perfil, layout, bloque, índice de intento y presupuesto"
)]
fn refine_block<'a>(
    scenario_id: &ScenarioId,
    root_seed: RootSeed,
    weights: &'a [f64],
    bins: &'a [LdBin],
    layout_block: &F5LayoutBlock,
    block: &mut F5BlockHaplotypes,
    attempt_index: u64,
    max_proposals: u64,
) -> Result<(BlockSearchMeasure, BlockSearchState<'a>), Diagnostic> {
    let block_id = block.block_id().clone();
    let internal = |detail: String| -> Diagnostic {
        Diagnostic::f5_internal(
            scenario_id.as_str(),
            format!("bloque {}: {detail}", block_id.as_str()),
        )
    };
    let population_id = block.population_id().clone();
    let mut state = BlockSearchState::seed(weights, bins, layout_block, block).map_err(internal)?;

    let mut proposals = 0_u64;
    let mut candidate_proposals = 0_u64;
    let mut accepted = 0_u64;
    let haplotype_count = weights.len() as u64;
    let site_count = state.site_count as u64;

    // El ámbito `swap_proposal` solo se abre si va a plantearse alguna
    // propuesta: un bloque que ya cumple la terminación no deriva subsemilla.
    let mut rng: Option<LocalRng> = (!state.has_terminated() && max_proposals > 0)
        .then(|| {
            open_swap_scope(
                root_seed,
                scenario_id,
                &population_id,
                &block_id,
                attempt_index,
            )
        })
        .transpose()?;

    while let Some(source) = rng.as_mut() {
        if state.has_terminated() || proposals >= max_proposals {
            break;
        }
        let proposal = draw_proposal(source, haplotype_count, site_count)
            .map_err(|detail| internal(detail.to_owned()))?;
        proposals += 1;
        let outcome = state
            .consider(block.matrix_mut_for_search(), proposal)
            .map_err(internal)?;
        match outcome {
            ProposalOutcome::NotCandidate => {}
            ProposalOutcome::Rejected => candidate_proposals += 1,
            ProposalOutcome::Accepted => {
                candidate_proposals += 1;
                accepted += 1;
            }
        }
    }

    // El veredicto solo puede calcularse sobre columnas canónicas.
    block
        .refresh_realized_columns(weights)
        .map_err(|detail| internal(detail.to_owned()))?;
    let measure = BlockSearchMeasure {
        proposals,
        candidate_proposals,
        accepted,
        row_visits: state.row_visits,
        pair_updates: state.pair_updates,
    };
    Ok((measure, state))
}

/// Abre el ámbito `swap_proposal` de DTC-089 para una unidad y un intento.
///
/// Los cuatro componentes del ámbito —escenario, subpoblación, bloque e índice
/// de intento— forman parte del contrato reproducible: dos intentos de la misma
/// unidad deben sortear secuencias distintas, y dos unidades del mismo intento
/// también.
fn open_swap_scope(
    root_seed: RootSeed,
    scenario_id: &ScenarioId,
    population_id: &PopulationId,
    block_id: &F5BlockId,
    attempt_index: u64,
) -> Result<LocalRng, Diagnostic> {
    Ok(derive_f5_subseed(
        root_seed,
        F5RandomDecision::swap_proposal(scenario_id, population_id, block_id, attempt_index),
    )?
    .local_rng())
}

/// Los cuatro enteros de una propuesta, ya incrementados y ordenados.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SwapProposal {
    row_a: usize,
    row_b: usize,
    site_i: usize,
    site_j: usize,
}

/// Consumo RNG literal de DTC-089.
///
/// Cuatro muestreos enteros sin sesgo en este orden: `a in 0..H`, `b in
/// 0..H-1`, `i in 0..S`, `j in 0..S-1`. Después, y no antes, `b>=a` y `j>=i`
/// se incrementan, y por último se ordena cada pareja. El orden y el número de
/// extracciones son contrato reproducible: una propuesta no candidata también
/// los consume, porque la candidatura se decide después de sortear.
fn draw_proposal<R: RandomWordSource>(
    rng: &mut R,
    haplotype_count: u64,
    site_count: u64,
) -> Result<SwapProposal, &'static str> {
    let bound = |value: u64| -> Result<u64, &'static str> {
        value
            .checked_sub(1)
            .filter(|reduced| *reduced >= 1)
            .ok_or("una propuesta 2x2 exige al menos dos filas y dos sitios")
    };
    let reduced_haplotypes = bound(haplotype_count)?;
    let reduced_sites = bound(site_count)?;
    let draw = |rng: &mut R, n: u64| -> Result<u64, &'static str> {
        f5_uniform_index_with_source(rng, n).map_err(|_| "el entero acotado rechazó su cota")
    };
    let a = draw(rng, haplotype_count)?;
    let mut b = draw(rng, reduced_haplotypes)?;
    let i = draw(rng, site_count)?;
    let mut j = draw(rng, reduced_sites)?;
    if b >= a {
        b += 1;
    }
    if j >= i {
        j += 1;
    }
    let (row_a, row_b) = if a <= b { (a, b) } else { (b, a) };
    let (site_i, site_j) = if i <= j { (i, j) } else { (j, i) };
    Ok(SwapProposal {
        row_a: usize::try_from(row_a).map_err(|_| "la fila sorteada no cabe en usize")?,
        row_b: usize::try_from(row_b).map_err(|_| "la fila sorteada no cabe en usize")?,
        site_i: usize::try_from(site_i).map_err(|_| "el sitio sorteado no cabe en usize")?,
        site_j: usize::try_from(site_j).map_err(|_| "el sitio sorteado no cabe en usize")?,
    })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ProposalOutcome {
    NotCandidate,
    Rejected,
    Accepted,
}

/// Estado de una pareja intrabloque dentro de la caché de puntuación.
///
/// `f5-resource-accounting-v7` presupuesta `192*S*(S-1)/2` bytes de parejas por
/// bloque×subpoblación; estos 32, más los 32 del veredicto canónico que la
/// guarda de cierre recoge, caben holgadamente, y son los que permiten que una
/// propuesta cueste `O(S)` en lugar de `O(H*S)`.
///
/// `d` se guarda porque es la **frontera exacta** de la regla de signos: la
/// guarda de cierre necesita saber si una divergencia cae justo en `D == 0.0`.
#[derive(Clone, Copy, Debug, Default)]
struct PairState {
    joint: f64,
    d: f64,
    r_squared: f64,
    evaluable: bool,
    sign_ok: bool,
    bin: u32,
}

/// Los seis campos numéricos de la tupla; el séptimo es la propia matriz.
#[derive(Clone, Copy, Debug, PartialEq)]
struct ScoreTuple {
    duplicate_rows: u64,
    bins_outside: u64,
    signs_failed: u64,
    max_r2_error: f64,
    sum_r2_error: f64,
    max_marginal_error: f64,
}

impl ScoreTuple {
    /// Comparación estricta de los seis primeros campos.
    ///
    /// Devuelve `None` cuando los seis coinciden: entonces decide
    /// `matrix_bytes`, que el llamador resuelve en `O(1)`.
    fn strictly_improves(&self, current: &Self) -> Option<bool> {
        let discrete = (self.duplicate_rows, self.bins_outside, self.signs_failed);
        let current_discrete = (
            current.duplicate_rows,
            current.bins_outside,
            current.signs_failed,
        );
        if discrete != current_discrete {
            return Some(discrete < current_discrete);
        }
        for (left, right) in [
            (self.max_r2_error, current.max_r2_error),
            (self.sum_r2_error, current.sum_r2_error),
            (self.max_marginal_error, current.max_marginal_error),
        ] {
            let ordering = left.total_cmp(&right);
            if ordering != std::cmp::Ordering::Equal {
                return Some(ordering == std::cmp::Ordering::Less);
            }
        }
        None
    }
}

/// Estado incremental de puntuación de un bloque×subpoblación.
struct BlockSearchState<'a> {
    weights: &'a [f64],
    bins: &'a [LdBin],
    site_count: usize,
    polarity_positive: Vec<bool>,
    targets: Vec<f64>,
    marginal_tolerance: f64,
    marginals: Vec<f64>,
    pairs: Vec<PairState>,
    bin_total: Vec<f64>,
    bin_count: Vec<u64>,
    bin_total_backup: Vec<f64>,
    bin_count_backup: Vec<u64>,
    evaluable_pairs: u64,
    signs_failed: u64,
    duplicates: DuplicateRowIndex,
    current: ScoreTuple,
    row_visits: u64,
    pair_updates: u64,
    touched: Vec<(usize, PairState)>,
    row_a_backup: Vec<u8>,
    row_b_backup: Vec<u8>,
    row_scratch: Vec<u8>,
}

impl<'a> BlockSearchState<'a> {
    /// Siembra el estado con una evaluación **canónica** del bloque.
    ///
    /// Reproduce literalmente lo que hace el evaluador LD: `p_hat_i` y
    /// `p_hat_ij` con Neumaier en orden canónico `h`, `D` y `r²` con la
    /// agrupación de DTC-089, y los bins acumulados en el orden
    /// `(site_index_i, site_index_j)`. Por eso la tupla de la propuesta cero
    /// coincide bit a bit con la que publica `evaluate_block`.
    fn seed(
        weights: &'a [f64],
        bins: &'a [LdBin],
        layout_block: &F5LayoutBlock,
        block: &F5BlockHaplotypes,
    ) -> Result<Self, String> {
        let site_count = block.matrix().site_count();
        if block.columns().len() != site_count || layout_block.sites().len() != site_count {
            return Err("el bloque y el layout no declaran el mismo número de sitios".to_owned());
        }
        if weights.len() != block.matrix().haplotype_count() {
            return Err("los pesos no acompañan a la matriz del bloque".to_owned());
        }
        let mut row_visits = 0_u64;
        let mut marginals = Vec::with_capacity(site_count);
        let mut targets = Vec::with_capacity(site_count);
        let mut polarity_positive = Vec::with_capacity(site_count);
        let mut positions = Vec::with_capacity(site_count);
        for (index, (column, site)) in block
            .columns()
            .iter()
            .zip(layout_block.sites().iter())
            .enumerate()
        {
            // La búsqueda lee del layout las posiciones físicas con las que
            // asigna los bins, así que comprueba la identidad antes de usarlas
            // en lugar de esperar al evaluador.
            if column.site_index() != index as u64
                || site.site_index() != index as u64
                || column.site_id() != site.site_id()
            {
                return Err(format!(
                    "el sitio {index} del intento y el del layout no coinciden"
                ));
            }
            marginals.push(
                site_marginal(weights, block.matrix(), index, &mut row_visits)
                    .map_err(|error| error.detail().to_owned())?,
            );
            targets.push(column.target_alt_af());
            polarity_positive.push(column.polarity() == crate::f5_haplotypes::F5Polarity::Positive);
            positions.push(site.position_1based());
        }

        let mut pairs = vec![PairState::default(); site_count * site_count.saturating_sub(1) / 2];
        let mut accumulators = vec![BinTotal::default(); bins.len()];
        let mut evaluable_pairs = 0_u64;
        let mut signs_failed = 0_u64;
        for i in 0..site_count {
            for j in (i + 1)..site_count {
                let joint = joint_marginal(weights, block.matrix(), i, j, &mut row_visits)
                    .map_err(|error| error.detail().to_owned())?;
                let distance = positions[j]
                    .checked_sub(positions[i])
                    .filter(|gap| *gap > 0)
                    .ok_or_else(|| {
                        format!("los sitios {i} y {j} no están en posiciones crecientes")
                    })?;
                let bin = bin_of(bins, distance).map_or(NO_BIN, |index| index as u32);
                let mut pair = evaluated_pair(
                    i,
                    j,
                    marginals[i],
                    marginals[j],
                    joint,
                    polarity_positive[i] == polarity_positive[j],
                    bin,
                );
                if pair.evaluable {
                    if pair.bin == NO_BIN {
                        return Err(format!(
                            "la distancia {distance} de la pareja ({i},{j}) no cae en ningún bin del perfil"
                        ));
                    }
                    evaluable_pairs += 1;
                    if !pair.sign_ok {
                        signs_failed += 1;
                    }
                    accumulators[pair.bin as usize].add(pair.r_squared)?;
                } else {
                    pair.sign_ok = false;
                }
                pairs[pair_index(i, j, site_count)] = pair;
            }
        }

        let bin_total = accumulators.iter().map(|total| total.value).collect();
        let bin_count = accumulators.iter().map(|total| total.count).collect();
        let mut duplicates = DuplicateRowIndex::new(fnv1a_64);
        for haplotype in 0..block.matrix().haplotype_count() {
            let row = block
                .matrix()
                .row_bytes(haplotype)
                .ok_or_else(|| "la matriz no expone una de sus filas".to_owned())?;
            duplicates.insert(row);
        }

        let bytes_per_row = site_count.div_ceil(8);
        let mut state = Self {
            weights,
            bins,
            site_count,
            polarity_positive,
            targets,
            marginal_tolerance: block.marginal_max_abs_error(),
            marginals,
            pairs,
            bin_total_backup: vec![0.0; bins.len()],
            bin_count_backup: vec![0; bins.len()],
            bin_total,
            bin_count,
            evaluable_pairs,
            signs_failed,
            duplicates,
            current: ScoreTuple {
                duplicate_rows: 0,
                bins_outside: 0,
                signs_failed: 0,
                max_r2_error: 0.0,
                sum_r2_error: 0.0,
                max_marginal_error: 0.0,
            },
            row_visits,
            pair_updates: 0,
            touched: Vec::with_capacity(2 * site_count),
            row_a_backup: vec![0; bytes_per_row],
            row_b_backup: vec![0; bytes_per_row],
            row_scratch: vec![0; bytes_per_row],
        };
        state.current = state.tuple()?;
        Ok(state)
    }

    /// Los tres agregados de bin, escritos como los escribe `evaluate_block`:
    /// recorrido en orden de bin, tolerancia inclusiva y Neumaier para la suma.
    fn tuple(&self) -> Result<ScoreTuple, String> {
        let mut bins_outside = 0_u64;
        let mut max_r2_error = 0.0_f64;
        let mut compensated = NeumaierSum::default();
        let mut sum_r2_error = 0.0_f64;
        for (index, definition) in self.bins.iter().enumerate() {
            if self.bin_count[index] == 0 {
                continue;
            }
            let mean = self.bin_total[index] / self.bin_count[index] as f64;
            if !mean.is_finite() {
                return Err(format!("la media del bin {index} no es finita"));
            }
            let error = f5_bin_abs_error(mean, definition);
            if f5_bin_is_outside(error, definition) {
                bins_outside += 1;
            }
            if error > max_r2_error {
                max_r2_error = error;
            }
            sum_r2_error = compensated.add(error).map_err(str::to_owned)?;
        }
        let mut max_marginal_error = 0.0_f64;
        for (realized, target) in self.marginals.iter().zip(self.targets.iter()) {
            let error = (realized - target).abs();
            if error > max_marginal_error {
                max_marginal_error = error;
            }
        }
        Ok(ScoreTuple {
            duplicate_rows: self.duplicates.duplicates,
            bins_outside,
            signs_failed: self.signs_failed,
            max_r2_error,
            sum_r2_error,
            max_marginal_error,
        })
    }

    /// Terminación de DTC-089 sobre el estado incremental.
    fn has_terminated(&self) -> bool {
        block_has_terminated(
            &self.current,
            self.evaluable_pairs > 0,
            self.marginal_tolerance,
        )
    }

    /// Evalúa una propuesta y la acepta solo si mejora estrictamente la tupla.
    fn consider(
        &mut self,
        matrix: &mut F5AlleleMatrix,
        proposal: SwapProposal,
    ) -> Result<ProposalOutcome, String> {
        let SwapProposal {
            row_a,
            row_b,
            site_i,
            site_j,
        } = proposal;
        if !matrix.is_swap_candidate(row_a, row_b, site_i, site_j) {
            return Ok(ProposalOutcome::NotCandidate);
        }
        // `a < b` y `i < j`, de modo que `(a, i)` es el primer bit que el
        // intercambio cambia dentro de `matrix_bytes`. Ese bit decide por sí
        // solo el orden lexicográfico de la codificación canónica.
        let alt_ai = matrix
            .allele(row_a, site_i)
            .ok_or_else(|| "la propuesta nombra una casilla fuera de la matriz".to_owned())?;

        self.bin_total_backup.copy_from_slice(&self.bin_total);
        self.bin_count_backup.copy_from_slice(&self.bin_count);
        let previous = self.current;
        let previous_signs_failed = self.signs_failed;
        let previous_evaluable_pairs = self.evaluable_pairs;
        let previous_marginal_i = self.marginals[site_i];
        let previous_marginal_j = self.marginals[site_j];
        self.touched.clear();

        copy_row(matrix, row_a, &mut self.row_a_backup)?;
        copy_row(matrix, row_b, &mut self.row_b_backup)?;
        self.duplicates.remove(&self.row_a_backup)?;
        self.duplicates.remove(&self.row_b_backup)?;

        matrix.apply_2x2_swap(row_a, row_b, site_i, site_j);

        copy_row(matrix, row_a, &mut self.row_scratch)?;
        self.duplicates.insert(&self.row_scratch);
        copy_row(matrix, row_b, &mut self.row_scratch)?;
        self.duplicates.insert(&self.row_scratch);

        self.apply_increments(matrix, proposal, alt_ai)?;

        let candidate = self.tuple()?;
        let accepted = proposal_is_accepted(&candidate, &previous, alt_ai);

        if accepted {
            self.current = candidate;
            return Ok(ProposalOutcome::Accepted);
        }

        // Reversión: el intercambio es su propia inversa y el resto del estado
        // vuelve desde las copias, de modo que una propuesta rechazada no deja
        // deriva numérica alguna.
        // El índice se deshace **antes** de revertir la matriz: las filas que
        // hay que retirar son las que la propuesta insertó, y solo existen
        // mientras el intercambio sigue aplicado.
        copy_row(matrix, row_a, &mut self.row_scratch)?;
        self.duplicates.remove(&self.row_scratch)?;
        copy_row(matrix, row_b, &mut self.row_scratch)?;
        self.duplicates.remove(&self.row_scratch)?;
        self.duplicates.insert(&self.row_a_backup);
        self.duplicates.insert(&self.row_b_backup);
        matrix.apply_2x2_swap(row_a, row_b, site_i, site_j);
        for (index, pair) in self.touched.drain(..) {
            self.pairs[index] = pair;
        }
        self.marginals[site_i] = previous_marginal_i;
        self.marginals[site_j] = previous_marginal_j;
        self.bin_total.copy_from_slice(&self.bin_total_backup);
        self.bin_count.copy_from_slice(&self.bin_count_backup);
        self.signs_failed = previous_signs_failed;
        self.evaluable_pairs = previous_evaluable_pairs;
        self.current = previous;
        Ok(ProposalOutcome::Rejected)
    }

    /// Actualizaciones `O(S)` de DTC-089 sobre marginales, conjuntas y bins.
    ///
    /// El intercambio solo mueve los alelos de dos filas en dos columnas, de
    /// modo que únicamente cambian `p_hat_i`, `p_hat_j` y las conjuntas que
    /// contienen a `i` o a `j`. La conjunta `p_hat_ij` **no** cambia: el
    /// producto `x_hi x_hj` valía cero en las dos filas antes y después.
    fn apply_increments(
        &mut self,
        matrix: &F5AlleleMatrix,
        proposal: SwapProposal,
        alt_ai: bool,
    ) -> Result<(), String> {
        let SwapProposal {
            row_a,
            row_b,
            site_i,
            site_j,
        } = proposal;
        let weight_a = *self
            .weights
            .get(row_a)
            .ok_or_else(|| "la propuesta nombra una fila sin peso".to_owned())?;
        let weight_b = *self
            .weights
            .get(row_b)
            .ok_or_else(|| "la propuesta nombra una fila sin peso".to_owned())?;
        // `10/01 -> 01/10` retira `a` de `i` y `b` de `j`; `01/10 -> 10/01`
        // hace lo contrario. Un solo signo describe los cuatro cambios.
        let delta_ai = if alt_ai { -1.0_f64 } else { 1.0_f64 };
        let delta_bi = -delta_ai;
        let delta_aj = -delta_ai;
        let delta_bj = delta_ai;

        // Marginales por incremento, en `O(1)`. DTC-089 fija en `O(S)` la
        // actualización de una propuesta y la fórmula de admisión contabiliza
        // ese coste exacto, de modo que rehacer aquí las dos sumas canónicas
        // —`2H` por propuesta candidata— subcontabilizaría el presupuesto justo
        // en la esquina `H=4096, S=12`. La consecuencia se conoce y está
        // vigilada: `m - q_a + q_b` difiere de la suma de Neumaier en un ULP y
        // la regla de signos es discontinua en `D=0`, así que la guarda de
        // cierre contrasta cada pareja contra el canon y solo exime las
        // divergencias que caen en esa frontera exacta.
        self.marginals[site_i] = self.marginals[site_i] + delta_ai * weight_a + delta_bi * weight_b;
        self.marginals[site_j] = self.marginals[site_j] + delta_aj * weight_a + delta_bj * weight_b;

        for site_k in 0..self.site_count {
            if site_k == site_i || site_k == site_j {
                continue;
            }
            let alt_ak = f64::from(
                matrix
                    .allele(row_a, site_k)
                    .ok_or_else(|| "la matriz no expone una casilla del bloque".to_owned())?,
            );
            let alt_bk = f64::from(
                matrix
                    .allele(row_b, site_k)
                    .ok_or_else(|| "la matriz no expone una casilla del bloque".to_owned())?,
            );
            self.refresh_pair(
                site_i,
                site_k,
                delta_ai * weight_a * alt_ak + delta_bi * weight_b * alt_bk,
            )?;
            self.refresh_pair(
                site_j,
                site_k,
                delta_aj * weight_a * alt_ak + delta_bj * weight_b * alt_bk,
            )?;
        }
        // La pareja `(i, j)` conserva su conjunta pero no su `D`: las dos
        // marginales acaban de cambiar.
        self.refresh_pair(site_i, site_j, 0.0)?;
        Ok(())
    }

    /// Rehace una pareja con la conjunta desplazada y actualiza sus agregados.
    fn refresh_pair(&mut self, left: usize, right: usize, delta: f64) -> Result<(), String> {
        let (i, j) = if left <= right {
            (left, right)
        } else {
            (right, left)
        };
        let index = pair_index(i, j, self.site_count);
        self.pair_updates += 1;
        let previous = self.pairs[index];
        self.touched.push((index, previous));
        if previous.evaluable {
            // `checked_sub` y no `-=`: un desajuste futuro envolvería en
            // release y dejaría `evaluable_pairs` gigante, con lo que la
            // validez del bloque sería siempre cierta.
            self.evaluable_pairs = self.evaluable_pairs.checked_sub(1).ok_or_else(|| {
                "el recuento de parejas evaluables se quedó sin margen".to_owned()
            })?;
            if !previous.sign_ok {
                self.signs_failed = self.signs_failed.checked_sub(1).ok_or_else(|| {
                    "el recuento de signos incumplidos se quedó sin margen".to_owned()
                })?;
            }
            let bin = previous.bin as usize;
            self.bin_total[bin] -= previous.r_squared;
            self.bin_count[bin] = self.bin_count[bin]
                .checked_sub(1)
                .ok_or_else(|| format!("el recuento del bin {bin} se quedó sin margen"))?;
        }
        let updated = evaluated_pair(
            i,
            j,
            self.marginals[i],
            self.marginals[j],
            previous.joint + delta,
            self.polarity_positive[i] == self.polarity_positive[j],
            previous.bin,
        );
        if updated.evaluable {
            if updated.bin == NO_BIN {
                return Err(format!(
                    "la pareja ({i},{j}) se volvió evaluable y su distancia no cae en ningún bin"
                ));
            }
            self.evaluable_pairs += 1;
            if !updated.sign_ok {
                self.signs_failed += 1;
            }
            let bin = updated.bin as usize;
            self.bin_total[bin] += updated.r_squared;
            self.bin_count[bin] += 1;
        }
        self.pairs[index] = updated;
        Ok(())
    }
}

/// Terminación de un bloque×subpoblación según DTC-089.
///
/// «Cero en los tres campos de error discreto —duplicados, bins y signos— y
/// todas las tolerancias». Las tolerancias son dos: la de cada bin, que
/// `bins_outside` ya resume, y la marginal `max_h(q_h)`, única cota marginal
/// activa del contrato. Un intercambio puede empeorarla, porque
/// `max_error_marginal` es el último campo numérico de la tupla y un campo
/// anterior puede compensarlo; la cota sigue siendo condición de terminación.
///
/// La validez se consulta **antes** que la tupla y es el primer factor: un
/// bloque sin ninguna pareja evaluable tiene cero en los tres campos discretos
/// y en los dos errores de `r²`, y aun así no ha resuelto nada.
fn block_has_terminated(tuple: &ScoreTuple, is_valid: bool, marginal_tolerance: f64) -> bool {
    is_valid
        && tuple.duplicate_rows == 0
        && tuple.bins_outside == 0
        && tuple.signs_failed == 0
        && tuple.max_marginal_error <= marginal_tolerance
}

/// Veredicto de los **siete** campos de DTC-089 sobre una propuesta candidata.
///
/// Los seis primeros los resuelve [`ScoreTuple::strictly_improves`]. Cuando
/// empatan decide `matrix_bytes`, y decide en `O(1)`: con `a < b` e `i < j`,
/// el bit `(a, i)` es el primero que el intercambio cambia dentro de la
/// codificación canónica —los bytes anteriores son idénticos—, de modo que la
/// codificación nueva es menor exactamente cuando ese bit valía 1. Un empate
/// en los siete campos es imposible después de un intercambio, porque ese bit
/// siempre cambia; un empate en los seis primeros con `matrix_bytes` mayor se
/// **rechaza**, como se rechaza cualquier empeoramiento.
fn proposal_is_accepted(candidate: &ScoreTuple, current: &ScoreTuple, alt_ai: bool) -> bool {
    match candidate.strictly_improves(current) {
        Some(verdict) => verdict,
        None => alt_ai,
    }
}

/// `D`, `r²` y el veredicto de signo de una pareja, con la aritmética del
/// evaluador LD y no una copia suya.
fn evaluated_pair(
    i: usize,
    j: usize,
    marginal_i: f64,
    marginal_j: f64,
    joint: f64,
    polarities_match: bool,
    bin: u32,
) -> PairState {
    let pair = pair_from_components(i as u64, j as u64, marginal_i, marginal_j, joint);
    match pair.r_squared() {
        Some(r_squared) => PairState {
            joint,
            d: pair.d(),
            r_squared,
            evaluable: true,
            sign_ok: sign_rule_is_satisfied(pair.d(), polarities_match),
            bin,
        },
        None => PairState {
            joint,
            d: pair.d(),
            r_squared: 0.0,
            evaluable: false,
            sign_ok: false,
            bin,
        },
    }
}

/// Posición de la pareja `(i, j)` con `i < j` dentro del triángulo superior.
const fn pair_index(i: usize, j: usize, site_count: usize) -> usize {
    i * (2 * site_count - i - 1) / 2 + (j - i - 1)
}

fn copy_row(
    matrix: &F5AlleleMatrix,
    haplotype: usize,
    destination: &mut [u8],
) -> Result<(), String> {
    let row = matrix
        .row_bytes(haplotype)
        .ok_or_else(|| "la matriz no expone una de sus filas".to_owned())?;
    if row.len() != destination.len() {
        return Err("la fila copiada no mide lo que la matriz declara".to_owned());
    }
    destination.copy_from_slice(row);
    Ok(())
}

/// Total de un bin sembrado con Neumaier, como lo hace el evaluador LD.
#[derive(Clone, Copy, Default)]
struct BinTotal {
    compensated: NeumaierSum,
    value: f64,
    count: u64,
}

impl BinTotal {
    fn add(&mut self, r_squared: f64) -> Result<(), String> {
        self.value = self.compensated.add(r_squared).map_err(str::to_owned)?;
        self.count += 1;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Filas duplicadas en O(S) por propuesta
// ---------------------------------------------------------------------------

/// FNV-1a de 64 bits sobre los bytes empaquetados de una fila.
///
/// La función solo reparte cubetas: el recuento nunca depende de ella, porque
/// una colisión se resuelve comparando los bytes. Está fijada aquí, y no
/// tomada de `RandomState`, para que dos ejecuciones recorran exactamente las
/// mismas cubetas.
fn fnv1a_64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Mapa de conteos de hashes de fila de DTC-089.
///
/// Mantiene `filas_duplicadas = Σ_v max(0, count(v)-1)` con dos retiradas y
/// dos inserciones por propuesta —solo las filas `a` y `b`—, y **contrasta los
/// bytes ante una colisión de hash**: dos vectores alélicos distintos que
/// compartan hash conservan entradas separadas. La comprobación no es
/// opcional; sin ella el recuento fundiría filas que no son iguales.
struct DuplicateRowIndex {
    hash: fn(&[u8]) -> u64,
    buckets: BTreeMap<u64, Vec<(Vec<u8>, u64)>>,
    duplicates: u64,
    byte_comparisons: u64,
}

impl DuplicateRowIndex {
    fn new(hash: fn(&[u8]) -> u64) -> Self {
        Self {
            hash,
            buckets: BTreeMap::new(),
            duplicates: 0,
            byte_comparisons: 0,
        }
    }

    fn insert(&mut self, row: &[u8]) {
        let bucket = self.buckets.entry((self.hash)(row)).or_default();
        for (candidate, count) in bucket.iter_mut() {
            self.byte_comparisons += 1;
            if candidate.as_slice() == row {
                if *count >= 1 {
                    self.duplicates += 1;
                }
                *count += 1;
                return;
            }
        }
        bucket.push((row.to_vec(), 1));
    }

    fn remove(&mut self, row: &[u8]) -> Result<(), String> {
        let key = (self.hash)(row);
        let bucket = self
            .buckets
            .get_mut(&key)
            .ok_or_else(|| "se retiró una fila que el índice no contenía".to_owned())?;
        let mut found = None;
        for (position, (candidate, _)) in bucket.iter().enumerate() {
            self.byte_comparisons += 1;
            if candidate.as_slice() == row {
                found = Some(position);
                break;
            }
        }
        let position =
            found.ok_or_else(|| "se retiró una fila que el índice no contenía".to_owned())?;
        let count = &mut bucket[position].1;
        if *count >= 2 {
            self.duplicates -= 1;
        }
        *count -= 1;
        if *count == 0 {
            bucket.swap_remove(position);
            if bucket.is_empty() {
                self.buckets.remove(&key);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{
        DuplicateRowIndex, ScoreTuple, SwapProposal, UnitOrderKey, block_has_terminated,
        draw_proposal, fnv1a_64, open_swap_scope, order_key_is_better, pair_index,
        proposal_is_accepted,
    };
    use crate::f5_haplotypes::F5AlleleMatrix;
    use crate::identity::{F5BlockId, PopulationId, ScenarioId};
    use crate::randomness::{RandomWordSource, RootSeed};
    use sha2::{Digest, Sha256};
    use std::cmp::Ordering;

    const MIDDLING: ScoreTuple = ScoreTuple {
        duplicate_rows: 1,
        bins_outside: 1,
        signs_failed: 1,
        max_r2_error: 0.25,
        sum_r2_error: 0.25,
        max_marginal_error: 0.25,
    };

    const PERFECT: ScoreTuple = ScoreTuple {
        duplicate_rows: 0,
        bins_outside: 0,
        signs_failed: 0,
        max_r2_error: 0.0,
        sum_r2_error: 0.0,
        max_marginal_error: 0.0,
    };

    /// Fuente que entrega palabras fijadas y entra en pánico si se le pide una
    /// de más: una prueba de consumo no puede pasar por accidente.
    struct ScriptedWords {
        words: Vec<u64>,
        consumed: usize,
    }

    impl ScriptedWords {
        fn new(words: Vec<u64>) -> Self {
            Self { words, consumed: 0 }
        }
    }

    impl RandomWordSource for ScriptedWords {
        fn next_word(&mut self) -> u64 {
            let word = self.words[self.consumed];
            self.consumed += 1;
            word
        }
    }

    /// Índice de pareja: biyección con el triángulo superior.
    #[test]
    fn pair_index_enumerates_the_upper_triangle_without_collisions() {
        for site_count in 2..24_usize {
            let mut seen = vec![false; site_count * (site_count - 1) / 2];
            for i in 0..site_count {
                for j in (i + 1)..site_count {
                    let index = pair_index(i, j, site_count);
                    assert!(!seen[index], "({i},{j}) repite el índice {index}");
                    seen[index] = true;
                }
            }
            assert!(seen.iter().all(|visited| *visited));
        }
    }

    /// El orden de extracción, los incrementos y la ordenación son literales.
    #[test]
    fn a_proposal_consumes_four_bounded_integers_in_the_contract_order() {
        // H = 8, S = 8: `a in 0..8`, `b in 0..7`, `i in 0..8`, `j in 0..7`.
        // Con `a=5`, `b=5` -> `b>=a` incrementa a 6; con `i=3`, `j=1` no.
        let mut rng = ScriptedWords::new(vec![5, 5, 3, 1]);
        let proposal = draw_proposal(&mut rng, 8, 8).expect("propuesta");
        assert_eq!(
            proposal,
            SwapProposal {
                row_a: 5,
                row_b: 6,
                site_i: 1,
                site_j: 3
            }
        );
        assert_eq!(rng.consumed, 4, "una propuesta consume cuatro enteros");
    }

    /// `b < a` no se incrementa y la pareja se ordena después.
    #[test]
    fn a_proposal_orders_each_pair_after_incrementing() {
        let mut rng = ScriptedWords::new(vec![6, 2, 7, 7]);
        let proposal = draw_proposal(&mut rng, 8, 8).expect("propuesta");
        // `b=2 < a=6` no se incrementa; ordenado da (2, 6).
        // `j=7 >= i=7` se incrementa a 8... pero `j in 0..7` nunca vale 7.
        assert_eq!(proposal.row_a, 2);
        assert_eq!(proposal.row_b, 6);
    }

    /// El sorteo nunca produce `a == b` ni `i == j`.
    #[test]
    fn no_proposal_can_name_the_same_row_or_the_same_site_twice() {
        for a in 0..5_u64 {
            for b in 0..4_u64 {
                for i in 0..5_u64 {
                    for j in 0..4_u64 {
                        let mut rng = ScriptedWords::new(vec![a, b, i, j]);
                        let proposal = draw_proposal(&mut rng, 5, 5).expect("propuesta");
                        assert_ne!(proposal.row_a, proposal.row_b);
                        assert_ne!(proposal.site_i, proposal.site_j);
                        assert!(proposal.row_a < proposal.row_b);
                        assert!(proposal.site_i < proposal.site_j);
                    }
                }
            }
        }
    }

    /// Solo `10/01 ↔ 01/10` es candidata: las otras catorce no lo son.
    #[test]
    fn only_the_two_crossed_configurations_are_candidates() {
        for mask in 0..16_u32 {
            let rows = [
                format!("{}{}", (mask >> 3) & 1, (mask >> 2) & 1),
                format!("{}{}", (mask >> 1) & 1, mask & 1),
            ];
            let matrix = F5AlleleMatrix::from_rows(&[&rows[0], &rows[1]]);
            let candidate = matrix.is_swap_candidate(0, 1, 0, 1);
            assert_eq!(
                candidate,
                mask == 0b1001 || mask == 0b0110,
                "la configuración {rows:?} se clasificó mal"
            );
        }
    }

    /// El intercambio es su propia inversa y solo toca cuatro casillas.
    #[test]
    fn the_swap_is_its_own_inverse_and_touches_four_cells() {
        let original = F5AlleleMatrix::from_rows(&["1010", "0110", "1111"]);
        let mut matrix = original.clone();
        matrix.apply_2x2_swap(0, 1, 0, 1);
        assert_ne!(matrix, original);
        assert_eq!(matrix.allele(0, 0), Some(false));
        assert_eq!(matrix.allele(0, 1), Some(true));
        assert_eq!(matrix.allele(0, 2), Some(true), "la columna 2 no cambia");
        matrix.apply_2x2_swap(0, 1, 0, 1);
        assert_eq!(matrix, original);
    }

    /// El desempate por `matrix_bytes` se resuelve en `O(1)` mirando `(a, i)`.
    ///
    /// La prueba materializa las dos codificaciones canónicas y comprueba que
    /// su orden lexicográfico coincide con la regla del bit. Recorre todas las
    /// parejas candidatas de una matriz pequeña, incluidas las que cruzan un
    /// límite de byte.
    #[test]
    fn the_matrix_bytes_tiebreak_matches_the_materialised_encoding() {
        let base = F5AlleleMatrix::from_rows(&[
            "100101100110",
            "011010011001",
            "110011001100",
            "001100110011",
        ]);
        let mut compared = 0_u32;
        for a in 0..4 {
            for b in (a + 1)..4 {
                for i in 0..12 {
                    for j in (i + 1)..12 {
                        if !base.is_swap_candidate(a, b, i, j) {
                            continue;
                        }
                        let before = base.canonical_bytes();
                        let alt_ai = base.allele(a, i).expect("casilla");
                        let mut swapped = base.clone();
                        swapped.apply_2x2_swap(a, b, i, j);
                        let after = swapped.canonical_bytes();
                        assert_eq!(
                            after < before,
                            alt_ai,
                            "({a},{b},{i},{j}) rompe la regla del bit (a,i)"
                        );
                        compared += 1;
                    }
                }
            }
        }
        assert!(compared > 50, "la prueba recorrió solo {compared} parejas");
    }

    /// `canonical_cmp` es el orden lexicográfico de `matrix_bytes`.
    #[test]
    fn canonical_cmp_agrees_with_the_materialised_encoding() {
        let matrices = [
            F5AlleleMatrix::from_rows(&["0000", "0000"]),
            F5AlleleMatrix::from_rows(&["0001", "0000"]),
            F5AlleleMatrix::from_rows(&["1000", "0000"]),
            F5AlleleMatrix::from_rows(&["0000", "0001"]),
            F5AlleleMatrix::from_rows(&["1111", "1111"]),
            F5AlleleMatrix::from_rows(&["01", "10"]),
            F5AlleleMatrix::from_rows(&["010101010101", "101010101010"]),
        ];
        for left in &matrices {
            for right in &matrices {
                assert_eq!(
                    left.canonical_cmp(right),
                    left.canonical_bytes().cmp(&right.canonical_bytes()),
                    "canonical_cmp no reproduce el orden de los bytes"
                );
            }
        }
    }

    /// El séptimo campo se decide por bytes y **no** por su digest.
    ///
    /// La prueba busca por sí sola una pareja de matrices cuyo orden de bytes
    /// sea el contrario al de sus SHA-256 y exige que gane la de bytes
    /// menores. Si la selección volviese a comparar digests, fallaría.
    #[test]
    fn the_selection_orders_matrix_bytes_and_not_their_digest() {
        let mut opposed = None;
        for left_mask in 0..64_u32 {
            for right_mask in 0..64_u32 {
                let left = matrix_from_mask(left_mask);
                let right = matrix_from_mask(right_mask);
                let bytes = left.canonical_bytes().cmp(&right.canonical_bytes());
                let digest = format!("{:x}", Sha256::digest(left.canonical_bytes()))
                    .cmp(&format!("{:x}", Sha256::digest(right.canonical_bytes())));
                if bytes == Ordering::Less && digest == Ordering::Greater {
                    opposed = Some((left, right));
                    break;
                }
            }
            if opposed.is_some() {
                break;
            }
        }
        let (smaller_bytes, larger_bytes) =
            opposed.expect("no se halló una pareja con órdenes opuestos");
        fn key<'a>(matrix: &'a F5AlleleMatrix, attempt_index: u64) -> UnitOrderKey<'a> {
            UnitOrderKey {
                terminated: false,
                valid: true,
                score: &MIDDLING,
                matrix,
                attempt_index,
            }
        }
        // Mismo `attempt_index` para que solo decida el séptimo campo.
        assert!(order_key_is_better(
            &key(&smaller_bytes, 7),
            &key(&larger_bytes, 7)
        ));
        assert!(!order_key_is_better(
            &key(&larger_bytes, 7),
            &key(&smaller_bytes, 7)
        ));
    }

    /// Matriz `2x3` cuyas seis casillas enumeran el entero `mask`.
    fn matrix_from_mask(mask: u32) -> F5AlleleMatrix {
        let row = |offset: u32| {
            (0..3)
                .map(|bit| {
                    if (mask >> (offset + bit)) & 1 == 1 {
                        '1'
                    } else {
                        '0'
                    }
                })
                .collect::<String>()
        };
        F5AlleleMatrix::from_rows(&[&row(0), &row(3)])
    }

    /// El último desempate es el menor `attempt_index`, y solo se alcanza
    /// cuando los siete campos anteriores empatan.
    #[test]
    fn the_last_tiebreak_is_the_smallest_attempt_index() {
        let matrix = F5AlleleMatrix::from_rows(&["1010", "0101"]);
        let other = F5AlleleMatrix::from_rows(&["0101", "1010"]);
        fn key<'a>(
            terminated: bool,
            valid: bool,
            score: &'a ScoreTuple,
            matrix: &'a F5AlleleMatrix,
            attempt_index: u64,
        ) -> UnitOrderKey<'a> {
            UnitOrderKey {
                terminated,
                valid,
                score,
                matrix,
                attempt_index,
            }
        }
        assert!(order_key_is_better(
            &key(false, true, &MIDDLING, &matrix, 3),
            &key(false, true, &MIDDLING, &matrix, 9)
        ));
        assert!(!order_key_is_better(
            &key(false, true, &MIDDLING, &matrix, 9),
            &key(false, true, &MIDDLING, &matrix, 3)
        ));
        // `matrix_bytes` precede al índice: `0101…` es menor que `1010…`.
        assert!(order_key_is_better(
            &key(false, true, &MIDDLING, &other, 9),
            &key(false, true, &MIDDLING, &matrix, 3)
        ));
        // Un intento degenerado no gana pese a tener la tupla perfecta.
        assert!(order_key_is_better(
            &key(false, true, &MIDDLING, &matrix, 9),
            &key(false, false, &PERFECT, &matrix, 0)
        ));
        // Y uno terminado gana a uno más afilado que no lo está.
        assert!(order_key_is_better(
            &key(true, true, &MIDDLING, &matrix, 9),
            &key(false, true, &PERFECT, &matrix, 0)
        ));
    }

    /// El ámbito `swap_proposal` depende de sus cuatro componentes.
    ///
    /// DTC-089 lo fija con escenario, subpoblación, bloque e índice de intento.
    /// Si alguno no participara, dos ámbitos distintos compartirían secuencia y
    /// la búsqueda repetiría propuestas entre intentos o entre unidades.
    #[test]
    fn the_swap_scope_depends_on_all_four_components() {
        let seed = RootSeed::parse(
            "seed",
            "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8",
        )
        .expect("semilla");
        let first_word = |scenario: &str, population: &str, block: &str, attempt: u64| {
            let scenario_id = ScenarioId::parse(scenario, "test").expect("escenario");
            let population_id = PopulationId::parse(population, "test").expect("población");
            let block_id = F5BlockId::parse(block, "test").expect("bloque");
            let mut rng = open_swap_scope(seed, &scenario_id, &population_id, &block_id, attempt)
                .expect("ámbito");
            rng.next_word()
        };
        const SCENARIO: &str = "scenario";
        const POPULATION: &str = "P1";
        const BLOCK: &str = "BLK_Q3K5H3MRNH7HHOJZ4RZWGYYSTY";
        const OTHER_BLOCK: &str = "BLK_KEG5VVL3SEVTQNX26QLL6NH34U";
        let baseline = first_word(SCENARIO, POPULATION, BLOCK, 0);
        assert_ne!(
            baseline,
            first_word(SCENARIO, POPULATION, BLOCK, 1),
            "el ámbito ignora el índice de intento"
        );
        assert_ne!(
            baseline,
            first_word(SCENARIO, "P2", BLOCK, 0),
            "el ámbito ignora la subpoblación"
        );
        assert_ne!(
            baseline,
            first_word(SCENARIO, POPULATION, OTHER_BLOCK, 0),
            "el ámbito ignora el bloque"
        );
        assert_ne!(
            baseline,
            first_word("other", POPULATION, BLOCK, 0),
            "el ámbito ignora el escenario"
        );
    }

    /// La terminación exige validez y **todas** las tolerancias.
    #[test]
    fn termination_needs_validity_and_every_tolerance() {
        assert!(block_has_terminated(&PERFECT, true, 0.25));
        assert!(
            !block_has_terminated(&PERFECT, false, 0.25),
            "un bloque sin parejas evaluables se dio por terminado"
        );
        // La tolerancia marginal es inclusiva y decide por sí sola.
        let at_the_bound = ScoreTuple {
            max_marginal_error: 0.25,
            ..PERFECT
        };
        assert!(block_has_terminated(&at_the_bound, true, 0.25));
        let beyond_the_bound = ScoreTuple {
            max_marginal_error: 0.250_000_000_000_000_06,
            ..PERFECT
        };
        assert!(
            !block_has_terminated(&beyond_the_bound, true, 0.25),
            "un error marginal por encima de max_h(q_h) terminó igualmente"
        );
        // Cada campo discreto la impide por separado.
        for tuple in [
            ScoreTuple {
                duplicate_rows: 1,
                ..PERFECT
            },
            ScoreTuple {
                bins_outside: 1,
                ..PERFECT
            },
            ScoreTuple {
                signs_failed: 1,
                ..PERFECT
            },
        ] {
            assert!(!block_has_terminated(&tuple, true, 0.25));
        }
        // Los dos errores de `r²` no participan: los resume `bins_outside`.
        let sloppy = ScoreTuple {
            max_r2_error: 0.9,
            sum_r2_error: 0.9,
            ..PERFECT
        };
        assert!(block_has_terminated(&sloppy, true, 0.25));
    }

    /// El séptimo campo decide el empate y un empate sigue sin aceptarse.
    ///
    /// `strictly_improves` devuelve `None` cuando los seis primeros campos
    /// coinciden; ahí manda `matrix_bytes`. Un intercambio que agranda la
    /// codificación canónica —el bit `(a, i)` valía 0— empeora el séptimo
    /// campo y se rechaza; uno que la reduce la mejora estrictamente y se
    /// acepta. Ningún empate se acepta «por ser igual».
    #[test]
    fn the_seventh_field_decides_and_a_tie_is_never_accepted() {
        assert!(
            !proposal_is_accepted(&MIDDLING, &MIDDLING, false),
            "un empate con matrix_bytes mayor se aceptó"
        );
        assert!(
            proposal_is_accepted(&MIDDLING, &MIDDLING, true),
            "un empate con matrix_bytes menor no se aceptó"
        );
        // El séptimo campo no rescata un empeoramiento anterior.
        let worse = ScoreTuple {
            signs_failed: 2,
            ..MIDDLING
        };
        assert!(!proposal_is_accepted(&worse, &MIDDLING, true));
        // Ni estorba a una mejora anterior.
        let better = ScoreTuple {
            signs_failed: 0,
            ..MIDDLING
        };
        assert!(proposal_is_accepted(&better, &MIDDLING, false));
    }

    /// La tupla rechaza el empate y el empeoramiento en cada campo.
    #[test]
    fn the_tuple_accepts_only_a_strict_improvement() {
        let base = ScoreTuple {
            duplicate_rows: 2,
            bins_outside: 1,
            signs_failed: 3,
            max_r2_error: 0.25,
            sum_r2_error: 0.5,
            max_marginal_error: 0.125,
        };
        assert_eq!(base.strictly_improves(&base), None, "el empate no decide");
        let better = ScoreTuple {
            duplicate_rows: 1,
            ..base
        };
        assert_eq!(better.strictly_improves(&base), Some(true));
        assert_eq!(base.strictly_improves(&better), Some(false));
        // Un campo posterior no compensa el empeoramiento de uno anterior.
        let mixed = ScoreTuple {
            duplicate_rows: 3,
            bins_outside: 0,
            signs_failed: 0,
            max_r2_error: 0.0,
            sum_r2_error: 0.0,
            max_marginal_error: 0.0,
        };
        assert_eq!(mixed.strictly_improves(&base), Some(false));
        // El orden de los tres campos `f64` también es el del contrato.
        let sharper_sum = ScoreTuple {
            sum_r2_error: 0.25,
            ..base
        };
        assert_eq!(sharper_sum.strictly_improves(&base), Some(true));
        let sharper_marginal_worse_sum = ScoreTuple {
            sum_r2_error: 0.75,
            max_marginal_error: 0.0,
            ..base
        };
        assert_eq!(
            sharper_marginal_worse_sum.strictly_improves(&base),
            Some(false)
        );
    }

    /// El recuento de duplicados sobrevive a una función de hash degenerada.
    ///
    /// Con un hash constante todas las filas comparten cubeta, de modo que el
    /// índice solo puede distinguirlas comparando bytes. Si la comparación
    /// desapareciera, el recuento fundiría filas distintas.
    #[test]
    fn the_duplicate_index_contrasts_bytes_when_hashes_collide() {
        for hash in [fnv1a_64, |_: &[u8]| 0_u64] {
            let mut index = DuplicateRowIndex::new(hash);
            for row in [b"\x01", b"\x02", b"\x01", b"\x03", b"\x01"] {
                index.insert(row);
            }
            assert_eq!(
                index.duplicates, 2,
                "tres copias de 0x01 son dos duplicados"
            );
            index.remove(b"\x01").expect("retirada");
            assert_eq!(index.duplicates, 1);
            index.remove(b"\x01").expect("retirada");
            assert_eq!(index.duplicates, 0);
            index.remove(b"\x02").expect("retirada");
            assert_eq!(index.duplicates, 0);
            assert!(index.remove(b"\x02").is_err(), "no queda ninguna 0x02");
        }
    }

    /// El índice incremental coincide con el recuento directo de la matriz.
    #[test]
    fn the_duplicate_index_matches_the_direct_count_after_every_swap() {
        let mut matrix = F5AlleleMatrix::from_rows(&["1010", "0101", "1010", "0110"]);
        let mut index = DuplicateRowIndex::new(fnv1a_64);
        for haplotype in 0..matrix.haplotype_count() {
            index.insert(matrix.row_bytes(haplotype).expect("fila"));
        }
        assert_eq!(index.duplicates, matrix.duplicate_rows());
        for (a, b, i, j) in [
            (0_usize, 1_usize, 0_usize, 1_usize),
            (2, 3, 1, 2),
            (0, 3, 0, 3),
        ] {
            if !matrix.is_swap_candidate(a, b, i, j) {
                continue;
            }
            let row_a = matrix.row_bytes(a).expect("fila").to_vec();
            let row_b = matrix.row_bytes(b).expect("fila").to_vec();
            index.remove(&row_a).expect("retirada");
            index.remove(&row_b).expect("retirada");
            matrix.apply_2x2_swap(a, b, i, j);
            index.insert(matrix.row_bytes(a).expect("fila"));
            index.insert(matrix.row_bytes(b).expect("fila"));
            assert_eq!(
                index.duplicates,
                matrix.duplicate_rows(),
                "el índice se desvió tras ({a},{b},{i},{j})"
            );
        }
        assert!(
            index.byte_comparisons > 0,
            "la comparación de bytes no ocurrió"
        );
    }
}
