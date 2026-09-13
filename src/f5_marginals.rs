//! Marginales objetivo de F5: `p_i` por subpoblación×bloque×sitio.
//!
//! DTC-089 cierra el espectro en dos decisiones separadas: elegir una fila por
//! peso y, solo si el intervalo elegido no es degenerado, sortear dentro de él
//! un `f64` uniforme `min+(max-min)U` con `U` de 53 bits en `[0,1)`. Cada
//! decisión consume exactamente un U53 y una elección con una sola alternativa
//! no consume ninguno.
//!
//! Este módulo no decide geometría ni haplotipos: recibe el layout ya
//! producido bajo una admisión sellada y devuelve el objetivo marginal de cada
//! columna. Como cada sitio×subpoblación abre su propio ámbito
//! `marginal_spectrum`, el resultado de una subpoblación no depende de cuántas
//! otras existan ni del orden en que se declaren.
//!
//! La cota `H<=2^S`, que DTC-089 exige *antes* de sortear marginales, se
//! comprueba aquí. La segunda, `H<=2^V` sobre las columnas no fijadas, no
//! puede comprobarse hasta conocer los marginales y por eso vive en
//! `f5_haplotypes`, que es también quien llama a [`build_f5_marginals`] desde
//! F5.2.1.

use crate::application::EffectiveF5Config;
use crate::diagnostics::Diagnostic;
use crate::f5_layout::F5PhysicalLayout;
use crate::f5_preflight::{F5OtherResources, SealedF5Admission};
use crate::f5_profile::F5Profile;
use crate::identity::{
    ExecutionId, F5BlockId, F5IdentityContext, F5IdentityRegistry, F5SiteId, PopulationId,
    ReferenceId, ScenarioId,
};
use crate::randomness::{
    F5RandomDecision, RandomWordSource, RootSeed, derive_f5_subseed, f5_uniform_53,
};

/// Fila del espectro `f5_snp_spectrum`: un intervalo cerrado de `alt_af` y su
/// peso relativo.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F5SpectrumBand {
    minimum: f64,
    maximum: f64,
    weight: f64,
}

impl F5SpectrumBand {
    #[must_use]
    pub const fn new(minimum: f64, maximum: f64, weight: f64) -> Self {
        Self {
            minimum,
            maximum,
            weight,
        }
    }

    #[must_use]
    pub const fn minimum(&self) -> f64 {
        self.minimum
    }

    #[must_use]
    pub const fn maximum(&self) -> f64 {
        self.maximum
    }

    #[must_use]
    pub const fn weight(&self) -> f64 {
        self.weight
    }

    /// Un intervalo degenerado ya determina su valor y por tanto no abre el
    /// uniforme de valor.
    #[must_use]
    pub fn is_degenerate(&self) -> bool {
        self.minimum == self.maximum
    }
}

/// Espectro validado de frecuencias objetivo.
///
/// Conserva el orden del perfil, que DTC-089 exige creciente y disjunto, y
/// exige pesos normalizados, porque la selección compara U con el prefijo sin
/// reescalar, tal como fija DTC-089. El total se comprueba con Neumaier antes de
/// resolver el prefijo.
#[derive(Clone, Debug, PartialEq)]
pub struct F5MarginalSpectrum {
    bands: Vec<F5SpectrumBand>,
}

impl F5MarginalSpectrum {
    /// Construye el espectro validado con independencia del lector de
    /// perfiles, igual que hace el soporte de longitudes de bloque.
    pub fn new(scenario_id: &ScenarioId, bands: Vec<F5SpectrumBand>) -> Result<Self, Diagnostic> {
        if bands.is_empty() {
            return Err(domain(scenario_id, "perfil sin espectro de frecuencias"));
        }
        let mut previous: Option<f64> = None;
        for band in &bands {
            if !band.minimum.is_finite()
                || !band.maximum.is_finite()
                || !(0.0..=1.0).contains(&band.minimum)
                || !(0.0..=1.0).contains(&band.maximum)
                || band.minimum > band.maximum
                || !band.weight.is_finite()
                || band.weight <= 0.0
            {
                return Err(domain(
                    scenario_id,
                    "intervalo o peso del espectro inválido",
                ));
            }
            if previous.is_some_and(|last: f64| band.minimum <= last) {
                return Err(domain(
                    scenario_id,
                    "intervalos del espectro solapados o desordenados",
                ));
            }
            previous = Some(band.maximum);
        }
        // La selección compara U con el prefijo sin reescalar, así que los
        // pesos han de estar normalizados, con la misma tolerancia que exige
        // el soporte de longitudes de F5.1.
        let mut total = NeumaierSum::default();
        for band in &bands {
            total
                .add(band.weight)
                .map_err(|detail| domain(scenario_id, detail))?;
        }
        let total = total
            .positive_total()
            .map_err(|detail| domain(scenario_id, detail))?;
        if (total - 1.0).abs() > 1e-9 {
            return Err(domain(
                scenario_id,
                "los pesos del espectro no suman uno dentro de tolerancia",
            ));
        }
        Ok(Self { bands })
    }

    pub fn from_profile(scenario_id: &ScenarioId, profile: &F5Profile) -> Result<Self, Diagnostic> {
        Self::new(
            scenario_id,
            profile
                .snp_spectrum()
                .iter()
                .map(|row| F5SpectrumBand::new(row.minimum(), row.maximum(), row.weight()))
                .collect(),
        )
    }

    #[must_use]
    pub fn bands(&self) -> &[F5SpectrumBand] {
        &self.bands
    }

    /// Sortea el marginal de un sitio dentro de una subpoblación abriendo su
    /// propio ámbito `marginal_spectrum`.
    pub fn sample_site(
        &self,
        root_seed: RootSeed,
        scenario_id: &ScenarioId,
        population_id: &PopulationId,
        block_id: &F5BlockId,
        site_id: &F5SiteId,
    ) -> Result<f64, Diagnostic> {
        let mut rng = derive_f5_subseed(
            root_seed,
            F5RandomDecision::marginal_spectrum(scenario_id, population_id, block_id, site_id),
        )?
        .local_rng();
        self.sample_with_source(&mut rng)
            .map_err(|detail| Diagnostic::f5_internal(scenario_id.as_str(), detail))
    }

    /// Mecánica cerrada de DTC-089: fila y valor, en ese orden, cada una con su
    /// propio U53 y solo cuando tiene más de una alternativa.
    fn sample_with_source(&self, rng: &mut impl RandomWordSource) -> Result<f64, &'static str> {
        let band = self.choose_band(rng)?;
        if band.is_degenerate() {
            return Ok(band.minimum);
        }
        let unit = f5_uniform_53(rng);
        let value = band.minimum + (band.maximum - band.minimum) * unit;
        if !(0.0..=1.0).contains(&value) {
            return Err("el marginal sorteado no pertenece a [0,1]");
        }
        Ok(value)
    }

    fn choose_band(&self, rng: &mut impl RandomWordSource) -> Result<F5SpectrumBand, &'static str> {
        let (first, rest) = self
            .bands
            .split_first()
            .ok_or("el espectro validado perdió sus filas")?;
        if rest.is_empty() {
            // Una elección con una sola alternativa no consume RNG.
            return Ok(*first);
        }
        // DTC-089: se recorre `q` en orden de fila y se elige la primera suma
        // acumulada estrictamente mayor que U. Los pesos que el perfil entrega
        // son los efectivos que DTC-049 ya normalizó, así que U no se reescala
        // por el total: hacerlo divergía de la decisión cuando la suma
        // compensada no daba exactamente 1, que es lo habitual con pesos no
        // diádicos.
        let threshold = f5_uniform_53(rng);
        let mut prefix = NeumaierSum::default();
        for band in &self.bands {
            if threshold < prefix.add(band.weight)? {
                return Ok(*band);
            }
        }
        // La última fila absorbe en exclusiva el residuo de redondeo.
        Ok(*self.bands.last().unwrap_or(first))
    }
}

/// Marginal objetivo de un sitio dentro de un bloque×subpoblación.
#[derive(Clone, Debug, PartialEq)]
pub struct F5SiteMarginal {
    site_index: u64,
    site_id: F5SiteId,
    target_alt_af: f64,
}

impl F5SiteMarginal {
    #[must_use]
    pub const fn site_index(&self) -> u64 {
        self.site_index
    }

    #[must_use]
    pub const fn site_id(&self) -> &F5SiteId {
        &self.site_id
    }

    #[must_use]
    pub const fn target_alt_af(&self) -> f64 {
        self.target_alt_af
    }
}

/// Marginales de un bloque completo para una subpoblación.
#[derive(Clone, Debug, PartialEq)]
pub struct F5BlockMarginals {
    population_id: PopulationId,
    block_index: u64,
    block_id: F5BlockId,
    sites: Vec<F5SiteMarginal>,
}

impl F5BlockMarginals {
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
    pub fn sites(&self) -> &[F5SiteMarginal] {
        &self.sites
    }
}

/// Tabla completa de marginales en el orden canónico subpoblación, bloque,
/// sitio.
#[derive(Clone, Debug, PartialEq)]
pub struct F5MarginalTable {
    blocks: Vec<F5BlockMarginals>,
}

impl F5MarginalTable {
    #[must_use]
    pub fn blocks(&self) -> &[F5BlockMarginals] {
        &self.blocks
    }

    #[must_use]
    pub fn block(
        &self,
        population_id: &PopulationId,
        block_index: u64,
    ) -> Option<&F5BlockMarginals> {
        self.blocks
            .iter()
            .find(|row| row.population_id == *population_id && row.block_index == block_index)
    }
}

/// Exige que el layout recibido sea el que esta misma ejecución produjo.
///
/// El ámbito `marginal_spectrum` se deriva de `block_id` y `site_id`, que
/// vienen del layout y no de la configuración. Comparar solo el número de
/// bloques dejaría pasar el layout de otra ejecución con la misma geometría, y
/// los marginales saldrían de identificadores ajenos sin que ninguna capa
/// posterior pudiera notarlo. Rederivar cada identificador es la misma
/// operación con la que el layout los emitió.
fn ensure_layout_belongs_to(
    config: &EffectiveF5Config,
    scenario_id: &ScenarioId,
    root_seed: RootSeed,
    execution_id: &ExecutionId,
    layout: &F5PhysicalLayout,
    block_count: u64,
) -> Result<(), Diagnostic> {
    let reference_id = ReferenceId::parse(config.reference_id(), "reference_id").map_err(|_| {
        configuration(
            config,
            "reference_id",
            "no se pudo reconstruir como identidad tipada",
        )
    })?;
    let context = F5IdentityContext::new(root_seed, execution_id, scenario_id, &reference_id);
    let mut registry = F5IdentityRegistry::default();
    let foreign = |detail: &'static str| configuration(config, "layout.identity", detail);
    // El recuento y la identidad son complementarios: dos ejecuciones que solo
    // difieren en el número de bloques comparten el identificador del bloque
    // cero, porque se deriva del índice, así que la identidad sola no
    // distinguiría la geometría.
    if u64::try_from(layout.blocks().len()).ok() != Some(block_count) {
        return Err(configuration(
            config,
            "preflight.block_count",
            "el layout recibido no pertenece a la admisión sellada de esta configuración",
        ));
    }
    for block in layout.blocks() {
        let expected = registry
            .derive_block(context, block.block_index())
            .map_err(|_| foreign("no se pudo rederivar el identificador de bloque"))?;
        if &expected != block.block_id() {
            return Err(foreign(
                "el layout recibido no lo produjo esta ejecución: su identificador de bloque no coincide",
            ));
        }
        for site in block.sites() {
            let expected = registry
                .derive_site(context, block.block_index(), site.site_index())
                .map_err(|_| foreign("no se pudo rederivar el identificador de sitio"))?;
            if &expected != site.site_id() {
                return Err(foreign(
                    "el layout recibido no lo produjo esta ejecución: su identificador de sitio no coincide",
                ));
            }
        }
    }
    Ok(())
}

/// Genera los marginales de todas las subpoblaciones sobre un layout ya
/// materializado.
///
/// Deriva de la configuración efectiva el escenario, la semilla, el espectro y
/// las subpoblaciones, de modo que ningún consumidor pueda combinar piezas de
/// ejecuciones distintas. Exige además la admisión sellada de esa misma
/// configuración: esta frontera abre ámbitos RNG y F5.1.2.5 prohíbe abrir
/// ninguno antes del sellado.
pub fn build_f5_marginals(
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
    execution_id: &ExecutionId,
    layout: &F5PhysicalLayout,
) -> Result<F5MarginalTable, Diagnostic> {
    if !sealed.matches_admission(config, F5OtherResources::default()) {
        return Err(configuration(
            config,
            "preflight.binding",
            "la admisión sellada no corresponde a esta configuración",
        ));
    }
    let scenario_id = ScenarioId::parse(config.scenario_id(), "scenario_id").map_err(|_| {
        configuration(
            config,
            "scenario_id",
            "no se pudo reconstruir como identidad tipada",
        )
    })?;
    let root_seed = RootSeed::parse("seed", config.seed())
        .map_err(|_| configuration(config, "seed", "la semilla efectiva no es canónica"))?;
    ensure_layout_belongs_to(
        config,
        &scenario_id,
        root_seed,
        execution_id,
        layout,
        sealed.block_count(),
    )?;
    let spectrum = F5MarginalSpectrum::from_profile(&scenario_id, config.profile())?;
    let declared = &config.declared_scenario().subpopulations;
    let mut population_ids = Vec::with_capacity(declared.len());
    for subpopulation in declared {
        population_ids.push(
            PopulationId::parse(subpopulation.population_id.value().clone(), "population_id")
                .map_err(|_| {
                    configuration(
                        config,
                        "subpopulations.population_id",
                        "no se pudo reconstruir como identidad tipada",
                    )
                })?,
        );
    }
    if population_ids.is_empty() {
        return Err(configuration(
            config,
            "subpopulations",
            "no hay ninguna subpoblación cuyos marginales generar",
        ));
    }
    // Orden canónico por `population_id`: reordenar el escenario no reordena
    // la tabla, y un identificador repetido no puede duplicar una columna.
    population_ids.sort();
    if population_ids.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(configuration(
            config,
            "subpopulations.population_id",
            "una subpoblación aparece más de una vez",
        ));
    }

    // `h_max`, no el H sorteado de ningún bloque: esta cota se comprueba antes
    // de abrir el primer ámbito aleatorio de marginales y, por tanto, mucho
    // antes de que exista un recuento haplotípico por bloque. Usar la cota
    // superior la hace más exigente y, sobre todo, independiente de la
    // semilla.
    let haplotype_count = config.profile().haplotype_count_bound();
    let mut blocks = Vec::with_capacity(population_ids.len().saturating_mul(layout.blocks().len()));
    for population_id in &population_ids {
        for block in layout.blocks() {
            // DTC-089 exige `H <= 2^S` **antes** de sortear marginales. Con S
            // grande la cota es implícita, porque H nunca supera 4096; el
            // desbordamiento del turno se trata como bloque incumplidor.
            let capacity = u32::try_from(block.sites().len())
                .ok()
                .and_then(|exponent| 1_usize.checked_shl(exponent));
            if capacity.is_some_and(|capacity| haplotype_count > capacity) {
                return Err(Diagnostic::f5_haplotype_cardinality(
                    scenario_id.as_str(),
                    population_id.as_str(),
                    block.block_id().as_str(),
                    format!(
                        "el perfil admite hasta {haplotype_count} haplotipos y el bloque solo admite {} con {} sitios",
                        capacity.unwrap_or(0),
                        block.sites().len()
                    ),
                ));
            }
            let mut sites = Vec::with_capacity(block.sites().len());
            for site in block.sites() {
                sites.push(F5SiteMarginal {
                    site_index: site.site_index(),
                    site_id: site.site_id().clone(),
                    target_alt_af: spectrum.sample_site(
                        root_seed,
                        &scenario_id,
                        population_id,
                        block.block_id(),
                        site.site_id(),
                    )?,
                });
            }
            blocks.push(F5BlockMarginals {
                population_id: population_id.clone(),
                block_index: block.block_index(),
                block_id: block.block_id().clone(),
                sites,
            });
        }
    }
    Ok(F5MarginalTable { blocks })
}

/// Medidor de visitas de fila de las sumas canónicas de F5.
///
/// Vive aquí, junto al acumulador compensado, porque lo comparten el
/// constructor haplotípico, el evaluador LD y la búsqueda: las tres sumas
/// canónicas de un intento —la siembra de la puntuación, la recomputación de
/// las columnas publicadas y la evaluación LD de cierre— cuentan en **el
/// mismo** ejemplar, de modo que la cifra publicada describe el intento y no
/// una de sus fases.
///
/// Lo que el medidor **no** puede garantizar es que un recorrido nuevo lo use:
/// quien escriba un bucle sobre las H filas con su propio contador quedaría
/// fuera de la cifra. Es una medida de las sumas canónicas que pasan por
/// `site_marginal`, `joint_marginal` y `refresh_realized_columns`, no una
/// prueba de ausencia universal.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct F5RowVisits(u64);

impl F5RowVisits {
    pub(crate) const fn add(&mut self, rows: u64) {
        self.0 += rows;
    }

    #[must_use]
    pub(crate) const fn count(self) -> u64 {
        self.0
    }
}

/// Acumulador compensado de Neumaier.
///
/// El soporte discreto de longitudes tiene el suyo dentro de `f5_layout`, que
/// F5.1.2 cerró; aquí se repite en lugar de reabrir aquel módulo sellado. Los
/// consumidores posteriores de DTC-089 —la construcción haplotípica de F5.2—
/// reutilizan este, que por eso no lleva mensajes propios del espectro.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct NeumaierSum {
    sum: f64,
    compensation: f64,
}

impl NeumaierSum {
    /// Acumula un sumando y devuelve el total corregido hasta ese punto.
    pub(crate) fn add(&mut self, value: f64) -> Result<f64, &'static str> {
        let next = self.sum + value;
        if !next.is_finite() {
            return Err("la suma compensada de pesos no es finita");
        }
        let correction = if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        let next_compensation = self.compensation + correction;
        if !next_compensation.is_finite() {
            return Err("la compensación de pesos no es finita");
        }
        self.sum = next;
        self.compensation = next_compensation;
        Ok(next + next_compensation)
    }

    fn positive_total(self) -> Result<f64, &'static str> {
        let total = self.sum + self.compensation;
        if total.is_finite() && total > 0.0 {
            Ok(total)
        } else {
            Err("la suma compensada de pesos no es finita y positiva")
        }
    }
}

fn domain(scenario_id: &ScenarioId, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::f5_domain(scenario_id.as_str(), detail)
}

fn configuration(
    config: &EffectiveF5Config,
    fields: impl Into<String>,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::f5_configuration(config.scenario_path().display().to_string(), fields, detail)
}

#[cfg(test)]
mod tests {
    use super::{F5MarginalSpectrum, F5SpectrumBand, NeumaierSum};
    use crate::identity::ScenarioId;
    use crate::randomness::RandomWordSource;

    /// Fuente que entrega palabras fijadas y cuenta cuántas se consumieron.
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

    /// Palabra cuyos 53 bits altos son exactamente `bits53`.
    const fn word_for(bits53: u64) -> u64 {
        bits53 << 11
    }

    const TWO_POW_53: u64 = 1_u64 << 53;

    /// Oráculo independiente de la fila: aritmética entera exacta sobre pesos
    /// diádicos, sin `f64`, Neumaier ni prefijos productivos.
    ///
    /// Con pesos `w_k = numerators[k] / denominator` y `U = bits53/2^53`, el
    /// umbral `U*1` cae en la fila `k` menor tal que
    /// `bits53 * denominator < prefix_k * 2^53`.
    fn expected_row(numerators: &[u64], denominator: u64, bits53: u64) -> usize {
        let mut prefix = 0_u128;
        for (index, numerator) in numerators.iter().enumerate() {
            prefix += u128::from(*numerator);
            if u128::from(bits53) * u128::from(denominator) < prefix * u128::from(TWO_POW_53) {
                return index;
            }
        }
        numerators.len() - 1
    }

    #[test]
    fn spectrum_rejects_rows_outside_the_normative_domain() {
        let scenario = scenario();
        for (label, bands) in [
            ("vacío", vec![]),
            (
                "mínimo mayor que máximo",
                vec![F5SpectrumBand::new(0.5, 0.25, 1.0)],
            ),
            (
                "máximo mayor que uno",
                vec![F5SpectrumBand::new(0.5, 1.5, 1.0)],
            ),
            ("mínimo negativo", vec![F5SpectrumBand::new(-0.5, 0.5, 1.0)]),
            ("no finito", vec![F5SpectrumBand::new(0.0, f64::NAN, 1.0)]),
            ("peso cero", vec![F5SpectrumBand::new(0.0, 1.0, 0.0)]),
            (
                "peso no finito",
                vec![F5SpectrumBand::new(0.0, 1.0, f64::INFINITY)],
            ),
            (
                "solapado",
                vec![
                    F5SpectrumBand::new(0.0, 0.5, 0.5),
                    F5SpectrumBand::new(0.5, 1.0, 0.5),
                ],
            ),
            (
                "desordenado",
                vec![
                    F5SpectrumBand::new(0.6, 1.0, 0.5),
                    F5SpectrumBand::new(0.0, 0.5, 0.5),
                ],
            ),
        ] {
            let error = F5MarginalSpectrum::new(&scenario, bands)
                .expect_err(&format!("{label} debe rechazarse"));
            assert_eq!(error.code(), "GPS-DOMAIN-101", "{label}");
        }
        assert!(
            F5MarginalSpectrum::new(
                &scenario,
                vec![
                    F5SpectrumBand::new(0.0, 0.0, 0.25),
                    F5SpectrumBand::new(0.25, 0.75, 0.5),
                    F5SpectrumBand::new(1.0, 1.0, 0.25),
                ],
            )
            .is_ok(),
            "intervalos disjuntos y degenerados en los extremos son válidos"
        );
    }

    #[test]
    fn a_single_row_never_consumes_the_row_uniform() {
        let spectrum =
            F5MarginalSpectrum::new(&scenario(), vec![F5SpectrumBand::new(0.25, 0.75, 1.0)])
                .expect("espectro de una fila");
        // La única palabra disponible es la del valor: si el generador
        // sorteara también la fila, la fuente entraría en pánico por agotarse.
        let mut source = ScriptedWords::new(&[word_for(TWO_POW_53 / 2)]);
        assert_eq!(spectrum.sample_with_source(&mut source), Ok(0.5));
        assert_eq!(source.consumed, 1);
    }

    #[test]
    fn a_degenerate_interval_never_consumes_the_value_uniform() {
        let spectrum =
            F5MarginalSpectrum::new(&scenario(), vec![F5SpectrumBand::new(0.375, 0.375, 1.0)])
                .expect("espectro degenerado de una fila");
        let mut none = ScriptedWords::new(&[]);
        assert_eq!(spectrum.sample_with_source(&mut none), Ok(0.375));
        assert_eq!(
            none.consumed, 0,
            "fila única e intervalo degenerado no consumen entropía"
        );

        let mixed = F5MarginalSpectrum::new(
            &scenario(),
            vec![
                F5SpectrumBand::new(0.0, 0.0, 0.5),
                F5SpectrumBand::new(0.25, 0.75, 0.5),
            ],
        )
        .expect("espectro mixto");
        // Primera fila: solo se consume el U53 de fila.
        let mut degenerate = ScriptedWords::new(&[word_for(0)]);
        assert_eq!(mixed.sample_with_source(&mut degenerate), Ok(0.0));
        assert_eq!(degenerate.consumed, 1);
        // Segunda fila: fila y valor.
        let mut continuous = ScriptedWords::new(&[word_for(TWO_POW_53 - 1), word_for(0)]);
        assert_eq!(mixed.sample_with_source(&mut continuous), Ok(0.25));
        assert_eq!(continuous.consumed, 2);
    }

    #[test]
    fn row_choice_matches_an_exact_integer_oracle_over_every_boundary() {
        let numerators = [1_u64, 1, 2, 4];
        let denominator = 8_u64;
        let spectrum = F5MarginalSpectrum::new(
            &scenario(),
            vec![
                F5SpectrumBand::new(0.0, 0.0, 0.125),
                F5SpectrumBand::new(0.125, 0.125, 0.125),
                F5SpectrumBand::new(0.25, 0.25, 0.25),
                F5SpectrumBand::new(0.5, 0.5, 0.5),
            ],
        )
        .expect("espectro de cuatro filas degeneradas");
        let values = [0.0, 0.125, 0.25, 0.5];
        let mut boundaries = vec![0_u64, 1, TWO_POW_53 - 1];
        // Fronteras exactas de prefijo y sus vecinos inmediatos.
        for prefix in [1_u64, 2, 4] {
            let exact = prefix * TWO_POW_53 / u64::from(8_u8);
            boundaries.extend([exact - 1, exact, exact + 1]);
        }
        for bits53 in boundaries {
            let mut source = ScriptedWords::new(&[word_for(bits53)]);
            let observed = spectrum
                .sample_with_source(&mut source)
                .expect("fila degenerada");
            let expected = values[expected_row(&numerators, denominator, bits53)];
            assert_eq!(observed, expected, "bits53={bits53}");
            assert_eq!(source.consumed, 1, "bits53={bits53}");
        }
    }

    #[test]
    fn the_value_uniform_is_the_literal_dtc_089_affine_map() {
        let spectrum =
            F5MarginalSpectrum::new(&scenario(), vec![F5SpectrumBand::new(0.25, 0.75, 1.0)])
                .expect("espectro de una fila continua");
        // Oráculo exacto: p = 1/4 + (1/2)*k/2^53 = (2^52 + k)/2^54, elegido con
        // `k` par para que el racional sea representable en binary64.
        for bits53 in [0_u64, 2, 1 << 20, TWO_POW_53 / 2, TWO_POW_53 - 2] {
            let mut source = ScriptedWords::new(&[word_for(bits53)]);
            let observed = spectrum.sample_with_source(&mut source).expect("valor");
            let numerator = u128::from(1_u64 << 52) + u128::from(bits53);
            let expected = numerator as f64 / (1_u128 << 54) as f64;
            assert_eq!(observed, expected, "bits53={bits53}");
            assert_eq!(source.consumed, 1);
        }
    }

    #[test]
    fn zero_and_one_survive_the_continuous_draw() {
        let lower = F5MarginalSpectrum::new(&scenario(), vec![F5SpectrumBand::new(0.0, 0.5, 1.0)])
            .expect("intervalo con cero");
        let mut zero = ScriptedWords::new(&[word_for(0)]);
        assert_eq!(
            lower.sample_with_source(&mut zero),
            Ok(0.0),
            "U=0 conserva el extremo inferior sin recortarlo a un interior abierto"
        );

        let upper = F5MarginalSpectrum::new(&scenario(), vec![F5SpectrumBand::new(0.5, 1.0, 1.0)])
            .expect("intervalo con uno");
        let mut one = ScriptedWords::new(&[word_for(TWO_POW_53 - 1)]);
        assert_eq!(
            upper.sample_with_source(&mut one),
            Ok(1.0),
            "0,5+(0,5)(1-2^-53) redondea a uno bajo ties-to-even"
        );

        let fixed = F5MarginalSpectrum::new(
            &scenario(),
            vec![
                F5SpectrumBand::new(0.0, 0.0, 0.5),
                F5SpectrumBand::new(1.0, 1.0, 0.5),
            ],
        )
        .expect("espectro fijado en los extremos");
        let mut absent = ScriptedWords::new(&[word_for(TWO_POW_53 - 1)]);
        assert_eq!(fixed.sample_with_source(&mut absent), Ok(1.0));
        assert_eq!(absent.consumed, 1);
    }

    #[test]
    fn every_draw_stays_inside_its_row_and_inside_zero_one() {
        let spectrum = F5MarginalSpectrum::new(
            &scenario(),
            vec![
                F5SpectrumBand::new(0.0, 0.1, 0.3),
                F5SpectrumBand::new(0.100_000_000_000_000_1, 0.9, 0.4),
                F5SpectrumBand::new(0.900_000_000_000_000_1, 1.0, 0.3),
            ],
        )
        .expect("tres intervalos disjuntos");
        for bits53 in [
            0_u64,
            1,
            TWO_POW_53 / 3,
            TWO_POW_53 / 2,
            TWO_POW_53 - 2,
            TWO_POW_53 - 1,
        ] {
            for value_bits in [0_u64, TWO_POW_53 / 2, TWO_POW_53 - 1] {
                let mut source = ScriptedWords::new(&[word_for(bits53), word_for(value_bits)]);
                let observed = spectrum.sample_with_source(&mut source).expect("marginal");
                assert!(
                    (0.0..=1.0).contains(&observed),
                    "bits53={bits53} value_bits={value_bits} produjo {observed}"
                );
                let band = spectrum
                    .bands()
                    .iter()
                    .find(|band| observed >= band.minimum() && observed <= band.maximum());
                assert!(
                    band.is_some(),
                    "bits53={bits53} value_bits={value_bits} produjo {observed} fuera de todo intervalo"
                );
            }
        }
    }

    #[test]
    fn unnormalized_weights_are_rejected() {
        // La selección compara U con el prefijo sin reescalar, de modo que unos
        // pesos que no sumen uno dejarían tramos de [0,1) sin fila.
        let error = F5MarginalSpectrum::new(
            &scenario(),
            vec![
                F5SpectrumBand::new(0.0, 0.0, 3.0),
                F5SpectrumBand::new(1.0, 1.0, 1.0),
            ],
        )
        .expect_err("los pesos sin normalizar se rechazan");
        assert_eq!(error.code(), "GPS-DOMAIN-101");
    }

    /// F51-R02: con pesos normalizados cuya suma compensada no es exactamente
    /// uno, la regla de DTC-089 —«la primera suma acumulada estrictamente mayor
    /// que U»— y la de reescalar por el total eligen filas distintas. Ninguna
    /// prueba lo distinguía porque todos los fixtures usaban pesos diádicos.
    #[test]
    fn the_rule_compares_u_with_the_prefix_without_rescaling_by_the_total() {
        let weights = [0.5554, 0.486, 0.401, 0.539, 0.501, 0.801];
        let raw_total: f64 = {
            let mut sum = NeumaierSum::default();
            for weight in weights {
                sum.add(weight).expect("suma finita");
            }
            sum.positive_total().expect("total positivo")
        };
        let normalized: Vec<f64> = weights.iter().map(|weight| weight / raw_total).collect();
        let mut total = NeumaierSum::default();
        for weight in &normalized {
            total.add(*weight).expect("suma finita");
        }
        let total = total.positive_total().expect("total positivo");
        assert_ne!(total, 1.0, "el fixture necesita una suma distinta de uno");
        assert!((total - 1.0).abs() <= 1e-9, "y dentro de la tolerancia");

        let bands: Vec<F5SpectrumBand> = normalized
            .iter()
            .enumerate()
            .map(|(index, weight)| {
                let value = index as f64 / 10.0;
                F5SpectrumBand::new(value, value, *weight)
            })
            .collect();
        let spectrum = F5MarginalSpectrum::new(&scenario(), bands).expect("espectro normalizado");

        // Prefijo cuarto: por debajo de U con la regla literal, por encima si se
        // reescala por el total, que es menor que uno.
        let mut prefix = NeumaierSum::default();
        let mut fourth = 0.0;
        for weight in normalized.iter().take(4) {
            fourth = prefix.add(*weight).expect("suma finita");
        }
        let bits53 = (fourth * TWO_POW_53 as f64) as u64;
        let mut source = ScriptedWords::new(&[word_for(bits53)]);
        let observed = spectrum
            .sample_with_source(&mut source)
            .expect("fila seleccionada");
        let literal = if (bits53 as f64 / TWO_POW_53 as f64) < fourth {
            0.3
        } else {
            0.4
        };
        assert_eq!(observed, literal, "la regla es la literal de DTC-089");
    }

    /// H4: con pesos no diádicos la compensación deja de ser cero, así que
    /// esta prueba distingue el acumulador compensado de una suma ingenua.
    /// Todos los fixtures anteriores usaban pesos exactos en binario y no la
    /// discriminaban.
    #[test]
    fn the_compensated_sum_differs_from_the_naive_one_with_inexact_weights() {
        let tiny = f64::EPSILON / 2.0;
        let mut naive = 1.0_f64;
        for _ in 0..4 {
            naive += tiny;
        }
        assert_eq!(naive, 1.0, "la suma ingenua pierde los cuatro sumandos");

        let mut compensated = NeumaierSum::default();
        let mut last = 0.0;
        for value in [1.0, tiny, tiny, tiny, tiny] {
            last = compensated.add(value).expect("suma finita");
        }
        assert!(
            last > 1.0,
            "la compensación conserva los sumandos que la suma ingenua pierde: {last}"
        );
        assert_eq!(last, 1.0 + 4.0 * tiny);
        assert_eq!(compensated.positive_total().expect("total positivo"), last);
    }
}
