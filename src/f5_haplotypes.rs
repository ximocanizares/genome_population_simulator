//! Construcción haplotípica inicial de F5: la matriz binaria H×S de un intento.
//!
//! DTC-089 cierra cómo se crea esa matriz **sin enumerar los `2^H`**
//! subconjuntos de filas. Para cada columna cuyo objetivo `p_i` no está fijado
//! se permutan las H filas con una Fisher–Yates independiente, se recorren sus
//! pesos acumulados y se elige, entre el prefijo anterior y el que cruza
//! `p_i`, el que minimiza `|p_hat_i - p_i|`; un empate conserva el prefijo
//! anterior. Una columna con `p_i` igual a 0 o a 1 es constante y **no abre
//! ningún ámbito RNG**.
//!
//! Por construcción `marginal_max_abs_error = max_h(q_h)` —con pesos uniformes,
//! `1/H`— y esa es la **única** cota marginal activa: el prefijo elegido y su
//! vecino difieren exactamente en el peso de una fila, de modo que el mejor de
//! los dos nunca dista de `p_i` más de la mitad de ese peso. El módulo la
//! reverifica sobre el `p_hat_i` normativo, que se suma en orden canónico `h`
//! con el acumulador de Neumaier y no reutiliza el prefijo de la permutación.
//!
//! Este hito encadena por primera vez layout, marginales y haplotipos: la
//! frontera pública exige la admisión sellada, deriva de ella los marginales y
//! solo entonces abre los ámbitos `initial_column` y `polarity`. La segunda
//! cota de cardinalidad de DTC-089 —`H <= 2^V` sobre las V columnas cuyo
//! objetivo no es 0/1— se comprueba después de esos sorteos marginales y antes
//! de abrir el intento de cada bloque×subpoblación.
//!
//! Fuera de alcance en F5.2.1 y previstos en los dos hitos siguientes: el
//! evaluador LD (`D`, `r²`, bins y signos) de F5.2.2 y las propuestas 2×2 que
//! eliminan los duplicados de F5.2.3. Aquí las filas duplicadas se cuentan y se
//! publican, nunca se fusionan ni se suman en silencio.

use std::collections::BTreeMap;

use crate::application::EffectiveF5Config;
use crate::diagnostics::Diagnostic;
use crate::f5_marginals::{F5BlockMarginals, F5MarginalTable, NeumaierSum};
use crate::identity::{F5BlockId, F5SiteId, PopulationId, ScenarioId};
use crate::randomness::{
    F5RandomDecision, RandomWordSource, RootSeed, derive_f5_subseed, f5_bernoulli_with_source,
    f5_fisher_yates_with_source,
};

/// Polaridad de un sitio dentro de una subpoblación.
///
/// DTC-089 la fija en `positive` para el primer SNP del bloque y la sortea con
/// Bernoulli(0,5) para cada sitio posterior. Vive en `f5-frequencies.tsv`, no
/// en la geometría común, porque es propia de cada población.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F5Polarity {
    Positive,
    Negative,
}

impl F5Polarity {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Positive => "positive",
            Self::Negative => "negative",
        }
    }
}

/// Matriz binaria H×S empaquetada en bits.
///
/// El empaquetado es literalmente el de DTC-089 para `matrix_bytes`: filas por
/// rango de haplotipo, sitio creciente y bit más significativo primero, con el
/// relleno bajo de la última posición a cero. Así la representación en memoria
/// y la codificación canónica de desempate son la misma secuencia de bytes y
/// no pueden divergir.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5AlleleMatrix {
    haplotype_count: usize,
    site_count: usize,
    bytes_per_row: usize,
    data: Vec<u8>,
}

impl F5AlleleMatrix {
    fn zeroed(haplotype_count: usize, site_count: usize) -> Option<Self> {
        let bytes_per_row = site_count.div_ceil(8);
        let payload = haplotype_count.checked_mul(bytes_per_row)?;
        Some(Self {
            haplotype_count,
            site_count,
            bytes_per_row,
            data: vec![0_u8; payload],
        })
    }

    #[must_use]
    pub const fn haplotype_count(&self) -> usize {
        self.haplotype_count
    }

    #[must_use]
    pub const fn site_count(&self) -> usize {
        self.site_count
    }

    /// Alelo alternativo del haplotipo en el sitio, o `None` fuera de rango.
    #[must_use]
    pub fn allele(&self, haplotype: usize, site: usize) -> Option<bool> {
        (haplotype < self.haplotype_count && site < self.site_count)
            .then(|| self.bit(haplotype, site))
    }

    /// Bytes empaquetados de una fila, o `None` fuera de rango.
    #[must_use]
    pub fn row_bytes(&self, haplotype: usize) -> Option<&[u8]> {
        (haplotype < self.haplotype_count).then(|| self.row(haplotype))
    }

    /// Codificación canónica de DTC-089: `H` y `S` big-endian seguidos de las
    /// filas empaquetadas.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(16 + self.data.len());
        bytes.extend_from_slice(&(self.haplotype_count as u64).to_be_bytes());
        bytes.extend_from_slice(&(self.site_count as u64).to_be_bytes());
        bytes.extend_from_slice(&self.data);
        bytes
    }

    /// Filas repetidas: `Σ_v max(0, count(v)-1)` sobre los vectores alélicos.
    ///
    /// Dos filas con el mismo vector son un duplicado inválido; el recuento las
    /// conserva para que la puntuación del intento pueda verlas.
    #[must_use]
    pub fn duplicate_rows(&self) -> u64 {
        let mut counts: BTreeMap<&[u8], u64> = BTreeMap::new();
        for haplotype in 0..self.haplotype_count {
            *counts.entry(self.row(haplotype)).or_default() += 1;
        }
        counts.values().map(|count| count - 1).sum()
    }

    /// Índices calculados por bucles acotados a `haplotype_count`/`site_count`:
    /// ninguna entrada de la usuaria alcanza este acceso.
    fn row(&self, haplotype: usize) -> &[u8] {
        let start = haplotype * self.bytes_per_row;
        &self.data[start..start + self.bytes_per_row]
    }

    fn bit(&self, haplotype: usize, site: usize) -> bool {
        let byte = self.data[haplotype * self.bytes_per_row + site / 8];
        (byte >> (7 - site % 8)) & 1 == 1
    }

    fn set(&mut self, haplotype: usize, site: usize) {
        self.data[haplotype * self.bytes_per_row + site / 8] |= 1_u8 << (7 - site % 8);
    }

    /// Constructor literal para las pruebas unitarias del evaluador LD.
    ///
    /// Existe solo bajo `cfg(test)`: la vía productiva sigue siendo
    /// [`build_f5_haplotype_attempt`], que es quien sostiene la relación entre
    /// la matriz y los marginales. Aquí una fila es una cadena de `0`/`1` y la
    /// matriz no pretende cumplir ninguna cota marginal.
    #[cfg(test)]
    pub(crate) fn from_rows(rows: &[&str]) -> Self {
        let site_count = rows.first().map_or(0, |row| row.len());
        let mut matrix = Self::zeroed(rows.len(), site_count).expect("matriz de prueba");
        for (haplotype, row) in rows.iter().enumerate() {
            assert_eq!(row.len(), site_count, "las filas deben medir lo mismo");
            for (site, allele) in row.bytes().enumerate() {
                match allele {
                    b'0' => {}
                    b'1' => matrix.set(haplotype, site),
                    other => panic!("alelo {other} inesperado"),
                }
            }
        }
        matrix
    }
}

/// Columna realizada de un bloque×subpoblación.
#[derive(Clone, Debug, PartialEq)]
pub struct F5HaplotypeColumn {
    site_index: u64,
    site_id: F5SiteId,
    target_alt_af: f64,
    realized_alt_af: f64,
    abs_error: f64,
    alternate_haplotypes: u64,
    polarity: F5Polarity,
}

impl F5HaplotypeColumn {
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

    /// `p_hat_i = Σ_h q_h x_hi`, sumado en orden canónico `h` con Neumaier.
    #[must_use]
    pub const fn realized_alt_af(&self) -> f64 {
        self.realized_alt_af
    }

    #[must_use]
    pub const fn abs_error(&self) -> f64 {
        self.abs_error
    }

    /// Cuántas filas llevan el alelo alternativo; una columna constante tiene 0
    /// o H.
    #[must_use]
    pub const fn alternate_haplotypes(&self) -> u64 {
        self.alternate_haplotypes
    }

    #[must_use]
    pub const fn is_fixed(&self) -> bool {
        is_fixed_target(self.target_alt_af)
    }

    #[must_use]
    pub const fn polarity(&self) -> F5Polarity {
        self.polarity
    }
}

/// Matriz y columnas de un bloque×subpoblación dentro de un intento.
#[derive(Clone, Debug, PartialEq)]
pub struct F5BlockHaplotypes {
    population_id: PopulationId,
    block_index: u64,
    block_id: F5BlockId,
    matrix: F5AlleleMatrix,
    columns: Vec<F5HaplotypeColumn>,
    marginal_max_abs_error: f64,
    max_abs_error: f64,
    duplicate_rows: u64,
    construction_operations: u64,
}

impl F5BlockHaplotypes {
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
    pub const fn matrix(&self) -> &F5AlleleMatrix {
        &self.matrix
    }

    #[must_use]
    pub fn columns(&self) -> &[F5HaplotypeColumn] {
        &self.columns
    }

    /// `max_h(q_h)`: la única cota marginal activa de DTC-089.
    #[must_use]
    pub const fn marginal_max_abs_error(&self) -> f64 {
        self.marginal_max_abs_error
    }

    /// Mayor `|p_hat_i - p_i|` observado; nunca supera la cota anterior.
    #[must_use]
    pub const fn max_abs_error(&self) -> f64 {
        self.max_abs_error
    }

    #[must_use]
    pub const fn duplicate_rows(&self) -> u64 {
        self.duplicate_rows
    }

    /// Operaciones elementales que costó construir el bloque: los `H-1`
    /// sorteos de cada permutación no constante, las visitas del recorrido de
    /// prefijos y las `H` visitas de la suma canónica de cada columna.
    ///
    /// Es la medida instrumental que acredita el «sin enumerar `2^H`»: crece de
    /// forma lineal en H, no exponencial.
    #[must_use]
    pub const fn construction_operations(&self) -> u64 {
        self.construction_operations
    }
}

/// Un intento haplotípico completo: un bloque×subpoblación por fila, en el
/// orden canónico que fija la tabla de marginales.
#[derive(Clone, Debug, PartialEq)]
pub struct F5HaplotypeAttempt {
    attempt_index: u64,
    blocks: Vec<F5BlockHaplotypes>,
}

impl F5HaplotypeAttempt {
    #[must_use]
    pub const fn attempt_index(&self) -> u64 {
        self.attempt_index
    }

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
            .find(|row| row.population_id == *population_id && row.block_index == block_index)
    }
}

/// Construye la matriz inicial de un intento haplotípico para todas las
/// subpoblaciones y bloques.
///
/// Recibe la tabla de marginales ya construida en lugar de derivarla: `p_i` no
/// depende del intento —su ámbito no lleva `attempt_index`— y la admisión
/// presupuesta esos sorteos una sola vez, sin factor de intentos. Derivarlos
/// por intento gastaría entropía y operaciones fuera de lo admitido, en
/// silencio, porque en esta fase no hay contabilidad en tiempo de ejecución.
/// Quien recorra los intentos construye la tabla una vez con
/// [`build_f5_marginals`], que ya exige la admisión sellada y la pertenencia
/// del layout a esta ejecución.
pub fn build_f5_haplotype_attempt(
    config: &EffectiveF5Config,
    marginals: &F5MarginalTable,
    attempt_index: u64,
) -> Result<F5HaplotypeAttempt, Diagnostic> {
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
    let weights = config.profile().haplotype_weights();
    let marginal_max_abs_error = validated_maximum_weight(&scenario_id, weights)?;

    // Segunda cota de cardinalidad de DTC-089, para todos los bloques antes de
    // abrir ningún intento: así ninguna subpoblación consume entropía cuando
    // otra ya demuestra que su geometría no admite H vectores distintos.
    for block in marginals.blocks() {
        ensure_block_represents_every_haplotype(&scenario_id, weights.len(), block)?;
    }

    let mut permutation = vec![0_u32; weights.len()];
    let mut blocks = Vec::with_capacity(marginals.blocks().len());
    for block in marginals.blocks() {
        blocks.push(build_block(
            &scenario_id,
            root_seed,
            weights,
            marginal_max_abs_error,
            attempt_index,
            block,
            &mut permutation,
        )?);
    }
    Ok(F5HaplotypeAttempt {
        attempt_index,
        blocks,
    })
}

/// `max_h(q_h)` sobre unos pesos que el perfil ya normalizó.
///
/// La comprobación es local y barata: el resto del módulo asume pesos finitos
/// y positivos, y un peso inválido produciría una cota marginal sin sentido.
fn validated_maximum_weight(scenario_id: &ScenarioId, weights: &[f64]) -> Result<f64, Diagnostic> {
    if weights.len() < 2 {
        return Err(Diagnostic::f5_domain(
            scenario_id.as_str(),
            "el perfil debe declarar al menos dos haplotipos",
        ));
    }
    let mut maximum = 0.0_f64;
    for weight in weights {
        if !weight.is_finite() || *weight <= 0.0 || *weight > 1.0 {
            return Err(Diagnostic::f5_domain(
                scenario_id.as_str(),
                "peso de haplotipo fuera de (0,1]",
            ));
        }
        if *weight > maximum {
            maximum = *weight;
        }
    }
    Ok(maximum)
}

/// Segunda cota de cardinalidad de DTC-089: `H <= 2^V`.
///
/// Ocurre tras los sorteos marginales —los hizo `build_f5_marginals`— y antes
/// de abrir un intento haplotípico: con V columnas no fijadas solo existen
/// `2^V` vectores alélicos distintos, de modo que un H mayor obligaría a
/// duplicar filas y ningún intercambio podría deshacerlo.
fn ensure_block_represents_every_haplotype(
    scenario_id: &ScenarioId,
    haplotype_count: usize,
    block: &F5BlockMarginals,
) -> Result<(), Diagnostic> {
    let variable_sites = block
        .sites()
        .iter()
        .filter(|site| !is_fixed_target(site.target_alt_af()))
        .count();
    // Con V grande la cota es implícita, porque H nunca supera 4096; el
    // desbordamiento del turno se trata como bloque que sí cabe.
    let capacity = u32::try_from(variable_sites)
        .ok()
        .and_then(|exponent| 1_usize.checked_shl(exponent));
    if capacity.is_some_and(|capacity| haplotype_count > capacity) {
        return Err(Diagnostic::f5_haplotype_cardinality(
            scenario_id.as_str(),
            block.population_id().as_str(),
            block.block_id().as_str(),
            format!(
                "el perfil declara {haplotype_count} haplotipos y las {variable_sites} columnas no fijadas del bloque solo distinguen {}",
                capacity.unwrap_or(0)
            ),
        ));
    }
    Ok(())
}

fn build_block(
    scenario_id: &ScenarioId,
    root_seed: RootSeed,
    weights: &[f64],
    marginal_max_abs_error: f64,
    attempt_index: u64,
    block: &F5BlockMarginals,
    permutation: &mut [u32],
) -> Result<F5BlockHaplotypes, Diagnostic> {
    let haplotype_count = weights.len();
    let site_count = block.sites().len();
    let internal = |detail: String| Diagnostic::f5_internal(scenario_id.as_str(), detail);

    let mut matrix = F5AlleleMatrix::zeroed(haplotype_count, site_count).ok_or_else(|| {
        internal(format!(
            "la matriz {haplotype_count}x{site_count} desborda el espacio direccionable"
        ))
    })?;
    let mut columns = Vec::with_capacity(site_count);
    let mut operations = 0_u64;
    let mut max_abs_error = 0.0_f64;

    for (column, site) in block.sites().iter().enumerate() {
        if site.site_index() != column as u64 {
            return Err(internal(format!(
                "el sitio {} no ocupa la columna {column} que su índice anuncia",
                site.site_index()
            )));
        }
        let target = site.target_alt_af();
        if !target.is_finite() || !(0.0..=1.0).contains(&target) {
            return Err(internal(format!(
                "el marginal del sitio {column} no pertenece a [0,1]"
            )));
        }

        // El ámbito `initial_column` solo se abre para una columna no
        // constante: `p_i = 0/1` no consume RNG.
        let mut scope = (!is_fixed_target(target))
            .then(|| {
                derive_f5_subseed(
                    root_seed,
                    F5RandomDecision::initial_column(
                        scenario_id,
                        block.population_id(),
                        block.block_id(),
                        attempt_index,
                        site.site_index(),
                    ),
                )
            })
            .transpose()?
            .map(crate::randomness::Subseed::local_rng);
        let selection = select_column(permutation, weights, target, scope.as_mut())
            .map_err(|detail| internal(format!("columna {column}: {detail}")))?;
        operations = operations.saturating_add(selection.operations);
        for row in &permutation[..selection.alternate_rows] {
            matrix.set(*row as usize, column);
        }

        // `p_hat_i` normativo: orden canónico `h` y Neumaier, no el prefijo de
        // la permutación.
        let mut compensated = NeumaierSum::default();
        let mut realized = 0.0_f64;
        for (haplotype, weight) in weights.iter().enumerate() {
            operations = operations.saturating_add(1);
            if matrix.bit(haplotype, column) {
                realized = compensated
                    .add(*weight)
                    .map_err(|detail| internal(format!("columna {column}: {detail}")))?;
            }
        }
        let abs_error = (realized - target).abs();
        if !abs_error.is_finite() || abs_error > marginal_max_abs_error {
            return Err(internal(format!(
                "la columna {column} realizó {realized} frente a {target}: su error {abs_error} supera la cota {marginal_max_abs_error}"
            )));
        }
        if abs_error > max_abs_error {
            max_abs_error = abs_error;
        }

        let polarity = site_polarity(
            scenario_id,
            root_seed,
            block,
            attempt_index,
            site.site_index(),
        )?;
        columns.push(F5HaplotypeColumn {
            site_index: site.site_index(),
            site_id: site.site_id().clone(),
            target_alt_af: target,
            realized_alt_af: realized,
            abs_error,
            alternate_haplotypes: selection.alternate_rows as u64,
            polarity,
        });
    }

    Ok(F5BlockHaplotypes {
        population_id: block.population_id().clone(),
        block_index: block.block_index(),
        block_id: block.block_id().clone(),
        duplicate_rows: matrix.duplicate_rows(),
        matrix,
        columns,
        marginal_max_abs_error,
        max_abs_error,
        construction_operations: operations,
    })
}

/// Polaridad del sitio: el primero es `positive` sin RNG y cada posterior abre
/// su propio ámbito `polarity` para un único U53.
fn site_polarity(
    scenario_id: &ScenarioId,
    root_seed: RootSeed,
    block: &F5BlockMarginals,
    attempt_index: u64,
    site_index: u64,
) -> Result<F5Polarity, Diagnostic> {
    if site_index == 0 {
        return Ok(F5Polarity::Positive);
    }
    let mut rng = derive_f5_subseed(
        root_seed,
        F5RandomDecision::polarity(
            scenario_id,
            block.population_id(),
            block.block_id(),
            attempt_index,
            site_index,
        ),
    )?
    .local_rng();
    draw_polarity(site_index, &mut rng)
        .map_err(|detail| Diagnostic::f5_internal(scenario_id.as_str(), detail))
}

fn draw_polarity(
    site_index: u64,
    rng: &mut impl RandomWordSource,
) -> Result<F5Polarity, &'static str> {
    if site_index == 0 {
        return Ok(F5Polarity::Positive);
    }
    // Decisión local registrada: DTC-089 fija la Bernoulli(0,5) pero no cuál de
    // las dos etiquetas es el éxito; el éxito es `positive`.
    match f5_bernoulli_with_source(rng, 0.5) {
        Ok(true) => Ok(F5Polarity::Positive),
        Ok(false) => Ok(F5Polarity::Negative),
        Err(_) => Err("la Bernoulli de polaridad rechazó la probabilidad 0,5"),
    }
}

/// Prefijo elegido para una columna y coste elemental de elegirlo.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ColumnSelection {
    alternate_rows: usize,
    operations: u64,
}

/// Elige las filas alternativas de una columna sin enumerar subconjuntos.
///
/// Deja en `permutation` la permutación de filas usada —identidad si la columna
/// es constante— y devuelve cuántas de sus primeras posiciones reciben el alelo
/// alternativo. Con `p_i` igual a 0 o a 1 no toca la fuente de entropía; en
/// cualquier otro caso exige el ámbito abierto y consume exactamente `H-1`
/// enteros acotados.
fn select_column<R: RandomWordSource>(
    permutation: &mut [u32],
    weights: &[f64],
    target: f64,
    rng: Option<&mut R>,
) -> Result<ColumnSelection, &'static str> {
    if permutation.len() != weights.len() {
        return Err("la permutación no tiene una posición por haplotipo");
    }
    for (index, slot) in permutation.iter_mut().enumerate() {
        *slot = u32::try_from(index).map_err(|_| "H no cabe en la permutación")?;
    }
    if target == 0.0 {
        return Ok(ColumnSelection {
            alternate_rows: 0,
            operations: 0,
        });
    }
    if target == 1.0 {
        return Ok(ColumnSelection {
            alternate_rows: permutation.len(),
            operations: 0,
        });
    }
    let rng = rng.ok_or("una columna no constante exige el ámbito initial_column")?;
    f5_fisher_yates_with_source(rng, permutation)
        .map_err(|_| "la permutación Fisher–Yates no pudo completarse")?;
    let mut operations = permutation.len().saturating_sub(1) as u64;

    let mut compensated = NeumaierSum::default();
    let mut previous = 0.0_f64;
    // Sin cruce —solo alcanzable si la suma de pesos se queda por debajo de
    // `p_i`— el conjunto completo es el prefijo más próximo.
    let mut alternate_rows = permutation.len();
    for (visited, row) in permutation.iter().enumerate() {
        operations += 1;
        let weight = *weights
            .get(*row as usize)
            .ok_or("la permutación apunta fuera de los pesos")?;
        let current = compensated.add(weight)?;
        if current >= target {
            // Empate: se conserva el prefijo anterior.
            alternate_rows = if current - target < target - previous {
                visited + 1
            } else {
                visited
            };
            break;
        }
        previous = current;
    }
    Ok(ColumnSelection {
        alternate_rows,
        operations,
    })
}

const fn is_fixed_target(target: f64) -> bool {
    target == 0.0 || target == 1.0
}

#[cfg(test)]
mod tests {
    use super::{
        F5AlleleMatrix, F5Polarity, draw_polarity, is_fixed_target, select_column,
        validated_maximum_weight,
    };
    use crate::identity::ScenarioId;
    use crate::randomness::{RandomWordSource, f5_matrix_bytes};

    /// Fuente que entrega palabras fijadas y entra en pánico si se le pide una
    /// de más: así una prueba de consumo cero no puede pasar por accidente.
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

    const fn word_for(bits53: u64) -> u64 {
        bits53 << 11
    }

    /// Oráculo entero de la permutación Fisher–Yates de DTC-089.
    ///
    /// Reimplementa el recorrido `i = len-1..1`, el rechazo desde `next_u64` y
    /// el resto, con aritmética `u128` y sin tocar ninguna función productiva.
    /// Devuelve la permutación y cuántas palabras costó: los `len-1` enteros
    /// acotados consumen una palabra cada uno salvo cuando el rechazo obliga a
    /// repetir el sorteo.
    fn oracle_permutation(length: usize, words: &[u64]) -> (Vec<usize>, usize) {
        let mut permutation = (0..length).collect::<Vec<_>>();
        let mut cursor = 0_usize;
        for index in (1..length).rev() {
            let n = (index + 1) as u128;
            let threshold = ((1_u128 << 64) / n) * n;
            let other = loop {
                let value = u128::from(words[cursor]);
                cursor += 1;
                if value < threshold {
                    break (value % n) as usize;
                }
            };
            permutation.swap(index, other);
        }
        assert!(
            cursor >= length.saturating_sub(1),
            "consumo de la permutación"
        );
        (permutation, cursor)
    }

    /// Oráculo entero del prefijo ponderado.
    ///
    /// Trabaja con numeradores sobre un denominador común y con el objetivo
    /// como racional exacto, de modo que ni la suma ni la comparación pasan
    /// por `f64`. Devuelve cuántas filas iniciales de la permutación reciben
    /// el alelo alternativo.
    fn oracle_prefix(
        permutation: &[usize],
        numerators: &[u128],
        denominator: u128,
        target_numerator: u128,
        target_denominator: u128,
    ) -> usize {
        let mut prefix = 0_u128;
        let mut previous = 0_u128;
        // Escala común: prefijo*td*? se compara con tn*D.
        let target = target_numerator * denominator;
        for (visited, row) in permutation.iter().enumerate() {
            prefix += numerators[*row];
            let current = prefix * target_denominator;
            if current >= target {
                let crossing_error = current - target;
                let previous_error = target - previous * target_denominator;
                return if crossing_error < previous_error {
                    visited + 1
                } else {
                    visited
                };
            }
            previous = prefix;
        }
        permutation.len()
    }

    /// Palabras que dejan la permutación en la identidad: en cada paso
    /// `i = len-1..1` el valor `i` es menor que el umbral de rechazo y su resto
    /// módulo `i+1` vuelve a ser `i`.
    fn identity_words(length: usize) -> Vec<u64> {
        (1..length).rev().map(|index| index as u64).collect()
    }

    #[test]
    fn a_fixed_target_builds_a_constant_column_without_touching_the_source() {
        let weights = vec![0.25_f64; 4];
        let mut permutation = vec![0_u32; 4];
        for (target, expected) in [(0.0_f64, 0_usize), (1.0, 4)] {
            // La fuente no tiene ni una palabra: pedir cualquiera sería pánico.
            let mut spy = ScriptedWords::new(&[]);
            let selection = select_column(&mut permutation, &weights, target, Some(&mut spy))
                .expect("columna constante");
            assert_eq!(selection.alternate_rows, expected, "objetivo {target}");
            assert_eq!(selection.operations, 0, "objetivo {target}");
            assert_eq!(spy.consumed, 0, "objetivo {target}");
            assert_eq!(permutation, vec![0, 1, 2, 3], "objetivo {target}");
        }
        assert!(is_fixed_target(0.0) && is_fixed_target(1.0));
        assert!(!is_fixed_target(0.5));
    }

    #[test]
    fn a_variable_column_consumes_exactly_h_minus_one_bounded_integers() {
        let weights = vec![0.125_f64; 8];
        let mut permutation = vec![0_u32; 8];
        let words = identity_words(8);
        let mut spy = ScriptedWords::new(&words);
        let selection = select_column(&mut permutation, &weights, 0.5, Some(&mut spy))
            .expect("columna no constante");
        assert_eq!(spy.consumed, 7, "H-1 enteros por permutación");
        assert_eq!(selection.alternate_rows, 4);
        // Operaciones = 7 sorteos + 4 visitas de prefijo.
        assert_eq!(selection.operations, 11);
    }

    #[test]
    fn a_variable_column_without_an_open_scope_is_refused() {
        let weights = vec![0.5_f64; 2];
        let mut permutation = vec![0_u32; 2];
        assert_eq!(
            select_column::<ScriptedWords>(&mut permutation, &weights, 0.5, None),
            Err("una columna no constante exige el ámbito initial_column")
        );
        let mut short = vec![0_u32; 1];
        assert_eq!(
            select_column::<ScriptedWords>(&mut short, &weights, 0.5, None),
            Err("la permutación no tiene una posición por haplotipo")
        );
    }

    /// Pesos que no alcanzan el objetivo: sin prefijo que lo cruce, el
    /// conjunto completo es el más próximo.
    #[test]
    fn a_target_above_the_total_weight_takes_every_row() {
        let weights = [0.1_f64, 0.1];
        let mut permutation = vec![0_u32; 2];
        let mut spy = ScriptedWords::new(&identity_words(2));
        let selection = select_column(&mut permutation, &weights, 0.5, Some(&mut spy))
            .expect("sin cruce posible");
        assert_eq!(selection.alternate_rows, 2);
        assert_eq!(spy.consumed, 1, "la permutación sigue costando H-1 enteros");
    }

    /// Pesos uniformes contra el oráculo entero, incluida una palabra
    /// rechazada por el muestreo acotado.
    #[test]
    fn uniform_weights_match_the_integer_oracle_over_permutation_and_prefix() {
        let haplotype_count = 8_usize;
        let weights = vec![0.125_f64; haplotype_count];
        let numerators = vec![2_u128; haplotype_count];
        let denominator = 16_u128;
        // `n = 3` es el único paso con umbral de rechazo `2^64-1`, así que
        // `u64::MAX` obliga a repetir el sorteo exactamente ahí.
        let word_sets: [Vec<u64>; 3] = [
            identity_words(haplotype_count),
            vec![3, 11, 29, 4, 2, 1, 0],
            vec![7, 6, 5, 4, 3, u64::MAX, 2, 1],
        ];
        for words in &word_sets {
            for (target, target_numerator) in [
                (0.5_f64, 8_u128),
                (0.0625, 1),
                (0.1875, 3),
                (0.3125, 5),
                (0.9375, 15),
            ] {
                let mut permutation = vec![0_u32; haplotype_count];
                let mut spy = ScriptedWords::new(words);
                let selection = select_column(&mut permutation, &weights, target, Some(&mut spy))
                    .expect("columna no constante");
                let (expected_permutation, expected_words) =
                    oracle_permutation(haplotype_count, words);
                assert_eq!(
                    permutation
                        .iter()
                        .map(|row| *row as usize)
                        .collect::<Vec<_>>(),
                    expected_permutation,
                    "permutación con {words:?}"
                );
                assert_eq!(spy.consumed, expected_words, "palabras con {words:?}");
                assert!(
                    selection.operations - (haplotype_count as u64 - 1) <= haplotype_count as u64,
                    "el recorrido de prefijos no supera H visitas"
                );
                let expected = oracle_prefix(
                    &expected_permutation,
                    &numerators,
                    denominator,
                    target_numerator,
                    16,
                );
                assert_eq!(
                    selection.alternate_rows, expected,
                    "objetivo {target} con {words:?}"
                );
            }
        }
    }

    /// Pesos no uniformes: la fila que cruza `p_i` ya no aporta siempre lo
    /// mismo, de modo que el oráculo distingue el prefijo correcto del vecino.
    #[test]
    fn non_uniform_weights_match_the_integer_oracle() {
        let weights = [0.125_f64, 0.125, 0.25, 0.5];
        let numerators = [2_u128, 2, 4, 8];
        let words_sets: [Vec<u64>; 3] = [identity_words(4), vec![0, 0, 0], vec![1, 2, 1]];
        for words in &words_sets {
            for (target, target_numerator) in [
                (0.0625_f64, 1_u128),
                (0.1875, 3),
                (0.5, 8),
                (0.6875, 11),
                (0.9375, 15),
            ] {
                let mut permutation = vec![0_u32; 4];
                let mut spy = ScriptedWords::new(words);
                let selection = select_column(&mut permutation, &weights, target, Some(&mut spy))
                    .expect("columna no constante");
                let (expected_permutation, expected_words) = oracle_permutation(4, words);
                assert_eq!(spy.consumed, expected_words, "palabras con {words:?}");
                let expected =
                    oracle_prefix(&expected_permutation, &numerators, 16, target_numerator, 16);
                assert_eq!(
                    selection.alternate_rows, expected,
                    "objetivo {target} con {words:?}"
                );
            }
        }
    }

    /// Empate exacto: `p_i` equidista del prefijo anterior y del que lo cruza,
    /// y DTC-089 conserva el anterior.
    #[test]
    fn an_exact_tie_keeps_the_previous_prefix() {
        let weights = vec![0.25_f64; 4];
        let mut permutation = vec![0_u32; 4];
        // Prefijos 0, 1/4, 1/2, 3/4, 1. Con 3/8 el anterior (1/4) y el que
        // cruza (1/2) distan ambos 1/8.
        let mut spy = ScriptedWords::new(&identity_words(4));
        let selection =
            select_column(&mut permutation, &weights, 0.375, Some(&mut spy)).expect("empate");
        assert_eq!(
            selection.alternate_rows, 1,
            "el empate conserva el anterior"
        );

        // Un objetivo apenas superior rompe el empate hacia el que cruza.
        let mut spy = ScriptedWords::new(&identity_words(4));
        let broken =
            select_column(&mut permutation, &weights, 0.376, Some(&mut spy)).expect("empate roto");
        assert_eq!(broken.alternate_rows, 2);

        // Y el empate del primer prefijo conserva el conjunto vacío.
        let mut spy = ScriptedWords::new(&identity_words(4));
        let empty =
            select_column(&mut permutation, &weights, 0.125, Some(&mut spy)).expect("empate en 0");
        assert_eq!(empty.alternate_rows, 0);
    }

    /// El recorrido usa el acumulador compensado: con pesos que no son
    /// diádicos exactos, una suma ingenua elegiría otro prefijo.
    #[test]
    fn the_prefix_walk_is_compensated_and_a_naive_sum_would_choose_another_prefix() {
        let tiny = f64::EPSILON / 4.0; // 2^-54
        let weights = [0.5_f64, tiny, tiny, tiny, 0.5 - 4.0 * tiny];
        assert_eq!(0.5 + tiny, 0.5, "la suma ingenua pierde cada sumando");
        let target = 0.5 + 2.0 * tiny;

        let mut naive = 0.0_f64;
        let mut naive_previous = 0.0_f64;
        let mut naive_rows = weights.len();
        for (visited, weight) in weights.iter().enumerate() {
            naive += *weight;
            if naive >= target {
                naive_rows = if naive - target < target - naive_previous {
                    visited + 1
                } else {
                    visited
                };
                break;
            }
            naive_previous = naive;
        }
        assert_eq!(naive_rows, 4, "la suma ingenua cruza en la última fila");

        let mut permutation = vec![0_u32; weights.len()];
        let mut spy = ScriptedWords::new(&identity_words(weights.len()));
        let selection = select_column(&mut permutation, &weights, target, Some(&mut spy))
            .expect("columna compensada");
        assert_eq!(
            selection.alternate_rows, 3,
            "la suma compensada alcanza el objetivo exactamente en la tercera fila"
        );
    }

    /// Medida instrumental de la ausencia de `2^H`.
    ///
    /// Para cada H se cuentan las palabras consumidas y las operaciones
    /// elementales de la columna. Ambas quedan por debajo de `2H`, y al doblar
    /// H no llegan a triplicarse: el coste es lineal, no exponencial. Con
    /// H=4096 una enumeración de subconjuntos necesitaría `2^4096` pasos.
    #[test]
    fn the_column_cost_grows_linearly_in_h_and_never_enumerates_subsets() {
        let mut previous: Option<(usize, u64)> = None;
        for exponent in 3..=12_u32 {
            let haplotype_count = 1_usize << exponent;
            let weights = vec![1.0 / haplotype_count as f64; haplotype_count];
            // Palabras deterministas y variadas, sin depender del generador
            // productivo.
            let words = (0..haplotype_count as u64 - 1)
                .map(|index| index.wrapping_mul(0x9e37_79b9_7f4a_7c15))
                .collect::<Vec<_>>();
            let mut permutation = vec![0_u32; haplotype_count];
            let mut spy = ScriptedWords::new(&words);
            let selection = select_column(&mut permutation, &weights, 0.5, Some(&mut spy))
                .expect("columna no constante");
            assert_eq!(spy.consumed, haplotype_count - 1, "H={haplotype_count}");
            assert!(
                selection.operations < 2 * haplotype_count as u64,
                "H={haplotype_count} costó {} operaciones",
                selection.operations
            );
            eprintln!(
                "medida H={haplotype_count}: palabras={} operaciones={}",
                spy.consumed, selection.operations
            );
            if let Some((previous_count, previous_operations)) = previous {
                assert_eq!(haplotype_count, previous_count * 2);
                assert!(
                    selection.operations < previous_operations * 3,
                    "al doblar H de {previous_count} a {haplotype_count} el coste pasó de {previous_operations} a {}",
                    selection.operations
                );
            }
            previous = Some((haplotype_count, selection.operations));
        }
        let (largest, operations) = previous.expect("la escalera recorrió al menos un H");
        assert_eq!(largest, 4096);
        assert!(operations < 12_288, "H=4096 costó {operations} operaciones");
    }

    #[test]
    fn the_bit_matrix_packs_the_canonical_dtc_089_bytes() {
        let mut matrix = F5AlleleMatrix::zeroed(3, 11).expect("matriz 3x11");
        let expected_rows = vec![
            vec![
                true, false, false, false, false, false, false, false, true, false, true,
            ],
            vec![false; 11],
            vec![true; 11],
        ];
        for (haplotype, row) in expected_rows.iter().enumerate() {
            for (site, allele) in row.iter().enumerate() {
                if *allele {
                    matrix.set(haplotype, site);
                }
            }
        }
        for (haplotype, row) in expected_rows.iter().enumerate() {
            for (site, allele) in row.iter().enumerate() {
                assert_eq!(matrix.allele(haplotype, site), Some(*allele));
            }
        }
        assert_eq!(matrix.allele(3, 0), None);
        assert_eq!(matrix.allele(0, 11), None);
        assert_eq!(
            matrix.canonical_bytes(),
            f5_matrix_bytes(&expected_rows).expect("codificación normativa"),
            "el empaquetado interno y `matrix_bytes` son la misma secuencia"
        );
        assert_eq!(matrix.row_bytes(2), Some([0xff_u8, 0xe0].as_slice()));
        assert_eq!(matrix.row_bytes(3), None);
    }

    #[test]
    fn duplicate_rows_count_every_repetition_without_merging_them() {
        let mut matrix = F5AlleleMatrix::zeroed(5, 3).expect("matriz 5x3");
        // Filas: 101, 101, 010, 101, 000 → tres iguales y dos únicas.
        for (haplotype, sites) in [
            (0, vec![0, 2]),
            (1, vec![0, 2]),
            (2, vec![1]),
            (3, vec![0, 2]),
        ] {
            for site in sites {
                matrix.set(haplotype, site);
            }
        }
        assert_eq!(matrix.duplicate_rows(), 2);

        let empty = F5AlleleMatrix::zeroed(4, 2).expect("matriz 4x2");
        assert_eq!(
            empty.duplicate_rows(),
            3,
            "cuatro filas nulas son idénticas"
        );
    }

    #[test]
    fn the_first_site_is_positive_without_rng_and_the_rest_draw_one_u53() {
        let mut none = ScriptedWords::new(&[]);
        assert_eq!(draw_polarity(0, &mut none), Ok(F5Polarity::Positive));
        assert_eq!(none.consumed, 0);

        // El umbral de Bernoulli(0,5) es `floor(0,5*2^53) = 2^52`.
        let mut success = ScriptedWords::new(&[word_for((1_u64 << 52) - 1)]);
        assert_eq!(draw_polarity(1, &mut success), Ok(F5Polarity::Positive));
        assert_eq!(success.consumed, 1);

        let mut failure = ScriptedWords::new(&[word_for(1_u64 << 52)]);
        assert_eq!(draw_polarity(7, &mut failure), Ok(F5Polarity::Negative));
        assert_eq!(failure.consumed, 1);

        assert_eq!(F5Polarity::Positive.as_str(), "positive");
        assert_eq!(F5Polarity::Negative.as_str(), "negative");
    }

    /// Los dos ámbitos de este hito reproducen su vector publicado en
    /// `f5-rng-catalog-v1`, y el vector se cierra hasta el valor producido: la
    /// permutación que sale de esas palabras y la polaridad que sale de la
    /// primera.
    #[test]
    fn initial_column_and_polarity_reproduce_their_published_catalog_vectors() {
        use crate::identity::{F5BlockId, PopulationId};
        use crate::randomness::{F5RandomDecision, RootSeed, derive_f5_subseed};

        const CATALOG: &str = include_str!("../catalogs/f5-rng-catalog-v1.tsv");
        let root = RootSeed::from_bytes(std::array::from_fn(|index| index as u8));
        let scenario_id = ScenarioId::parse("scenario-f5", "test").expect("scenario");
        let population_id =
            PopulationId::parse("population-A".to_owned(), "test").expect("population");
        let block_id = F5BlockId::parse("BLK_blockA", "test").expect("block");

        let row = |decision: &str| {
            let line = CATALOG
                .lines()
                .find(|line| line.starts_with(&format!("{decision}|")))
                .expect("la línea publicada del ámbito");
            let fields = line.split('|').map(str::to_owned).collect::<Vec<_>>();
            assert_eq!(fields.len(), 9);
            fields
        };

        let column_row = row("initial_column");
        assert_eq!(
            column_row[1],
            "scenario_id=scenario-f5;population_id=population-A;block_id=BLK_blockA;stage=initial_column;attempt_index_u64=3;salt.site_index_u64=5"
        );
        let subseed = derive_f5_subseed(
            root,
            F5RandomDecision::initial_column(&scenario_id, &population_id, &block_id, 3, 5),
        )
        .expect("subsemilla");
        assert_eq!(subseed.to_hex(), column_row[4], "subsemilla publicada");
        let mut rng = subseed.local_rng();
        let words: [u64; 4] = std::array::from_fn(|_| rng.next_u64());
        for (index, word) in words.iter().enumerate() {
            assert_eq!(
                format!("{word:016x}"),
                column_row[index + 5],
                "palabra {index}"
            );
        }

        // Cierre del vector: esas palabras construyen una columna concreta.
        let weights = vec![0.25_f64; 4];
        let mut permutation = vec![0_u32; 4];
        let mut published = ScriptedWords::new(&words);
        let selection = select_column(&mut permutation, &weights, 0.5, Some(&mut published))
            .expect("columna del ámbito publicado");
        let (expected_permutation, expected_words) = oracle_permutation(4, &words);
        assert_eq!(published.consumed, expected_words);
        assert_eq!(
            permutation
                .iter()
                .map(|row| *row as usize)
                .collect::<Vec<_>>(),
            expected_permutation
        );
        assert_eq!(
            selection.alternate_rows,
            oracle_prefix(&expected_permutation, &[1, 1, 1, 1], 4, 1, 2)
        );

        let polarity_row = row("polarity");
        assert_eq!(
            polarity_row[1],
            "scenario_id=scenario-f5;population_id=population-A;block_id=BLK_blockA;stage=polarity;attempt_index_u64=3;salt.site_index_u64=5"
        );
        let subseed = derive_f5_subseed(
            root,
            F5RandomDecision::polarity(&scenario_id, &population_id, &block_id, 3, 5),
        )
        .expect("subsemilla");
        assert_eq!(subseed.to_hex(), polarity_row[4], "subsemilla publicada");
        let mut rng = subseed.local_rng();
        let polarity_words: [u64; 4] = std::array::from_fn(|_| rng.next_u64());
        for (index, word) in polarity_words.iter().enumerate() {
            assert_eq!(
                format!("{word:016x}"),
                polarity_row[index + 5],
                "palabra {index}"
            );
        }
        // Cierre del vector: la primera palabra decide la polaridad.
        let expected = if polarity_words[0] >> 11 < (1_u64 << 52) {
            F5Polarity::Positive
        } else {
            F5Polarity::Negative
        };
        let mut published = ScriptedWords::new(&polarity_words[..1]);
        assert_eq!(draw_polarity(5, &mut published), Ok(expected));
        assert_eq!(
            published.consumed, 1,
            "un U53 por sitio posterior al primero"
        );
    }

    #[test]
    fn the_marginal_bound_is_the_largest_weight() {
        let scenario = ScenarioId::parse("scenario", "test").expect("scenario");
        assert_eq!(
            validated_maximum_weight(&scenario, &[0.25, 0.5, 0.25]),
            Ok(0.5)
        );
        assert_eq!(
            validated_maximum_weight(&scenario, &[0.25_f64; 4]),
            Ok(0.25),
            "con pesos uniformes la cota es 1/H"
        );
        for invalid in [
            vec![1.0_f64],
            vec![0.5, 0.0],
            vec![0.5, -0.5],
            vec![0.5, f64::NAN],
        ] {
            assert_eq!(
                validated_maximum_weight(&scenario, &invalid)
                    .expect_err("peso inválido")
                    .code(),
                "GPS-DOMAIN-101"
            );
        }
    }
}
