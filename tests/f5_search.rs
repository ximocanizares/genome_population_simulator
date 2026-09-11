//! F5.2.3. Búsqueda local por propuestas 2×2 y selección del mejor intento.
//!
//! Las pruebas recorren la cadena productiva completa —configuración,
//! admisión sellada, layout, marginales, construcción, propuestas y evaluación
//! LD— y contrastan el resultado con un oráculo racional independiente escrito
//! en enteros de 128 bits. Con pesos uniformes `1/H` todo el cálculo de
//! DTC-089 es un cociente de enteros exactos, de modo que el oráculo no
//! reutiliza ninguna función del productor.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::f5_haplotypes::{
    F5BlockHaplotypes, build_f5_block_attempt, build_f5_haplotype_attempt,
};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_ld::evaluate_f5_block_linkage_disequilibrium;
use genome_population_simulator::f5_marginals::build_f5_marginals;
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::f5_search::{
    F5AttemptRecord, F5HaplotypeSearch, F5SearchFailure, search_f5_haplotypes,
};
use genome_population_simulator::identity::{ExecutionId, ScenarioId};
use genome_population_simulator::randomness::{
    F5RandomDecision, RootSeed, derive_f5_subseed, f5_uniform_index,
};
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

const EXECUTION: &str = "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4";
const SEED: &str = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8";

/// Espectro de una sola banda continua `[1/4, 3/4]`: ninguna columna queda
/// fijada y toda pareja puede ser evaluable.
const CONTINUOUS_SPECTRUM: &str = "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t1\n";

/// Espectro degenerado en 0: todas las columnas quedan fijadas y **ninguna**
/// pareja puede tener denominador no nulo.
const FIXED_SPECTRUM: &str = "alt_af_min\talt_af_max\tweight\n0\t0\t1\n";

fn sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-search-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

struct Fixture<'a> {
    populations: &'a [&'a str],
    block_count: u64,
    block_length: u64,
    weights: Vec<f64>,
    gap: u64,
    /// Máximo del intervalo de gaps; por omisión, igual a `gap`.
    gap_max: Option<u64>,
    spectrum: &'a str,
    /// Bins LD `(inicio, fin, media objetivo, tolerancia)`; cubren `1..=L-1`.
    bins: Vec<(u64, u64, f64, f64)>,
    max_haplotype_attempts: u64,
    max_swap_proposals: u64,
    max_haplotype_cells: u64,
    seed: &'a str,
    /// Cuántas veces cabe el conjunto de bloques en la referencia. `1` obliga
    /// a un teselado exacto; un valor mayor deja holgura al layout.
    slack: u64,
}

impl Default for Fixture<'_> {
    fn default() -> Self {
        Self {
            populations: &["P1"],
            block_count: 2,
            block_length: 4,
            weights: vec![0.25; 4],
            gap: 1,
            gap_max: None,
            spectrum: CONTINUOUS_SPECTRUM,
            bins: vec![(1, 3, 0.5, 0.05)],
            max_haplotype_attempts: 8,
            max_swap_proposals: 64,
            max_haplotype_cells: 65_536,
            seed: SEED,
            slack: 1,
        }
    }
}

impl Fixture<'_> {
    fn write(&self, root: &Path) -> PathBuf {
        let scenario_directory = root.join("authoring/scenarios");
        let profile_directory = scenario_directory.join("profiles");
        fs::create_dir_all(&profile_directory).expect("profile directory");
        fs::create_dir_all(root.join("authoring/references")).expect("reference directory");

        let mut ld_rows = String::from(
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\n",
        );
        let mut sorted = self.populations.to_vec();
        sorted.sort_unstable();
        for population in &sorted {
            for (start, end, target, tolerance) in &self.bins {
                ld_rows.push_str(&format!(
                    "{population}\t{start}\t{end}\t{target}\t{tolerance}\n"
                ));
            }
        }

        let length = self.block_length;
        let mut weight_rows = String::from("haplotype_rank\tweight\n");
        for (rank, weight) in self.weights.iter().enumerate() {
            weight_rows.push_str(&format!("{rank}\t{weight}\n"));
        }
        let files = [
            (
                "block_lengths.tsv",
                "f5_block_length_profile",
                format!("length_min_bp\tlength_max_bp\tweight\n{length}\t{length}\t1\n"),
            ),
            (
                "snp_gaps.tsv",
                "f5_snp_gap_profile",
                format!(
                    "gap_min_bp\tgap_max_bp\tweight\n{}\t{}\t1\n",
                    self.gap,
                    self.gap_max.unwrap_or(self.gap)
                ),
            ),
            (
                "snp_spectrum.tsv",
                "f5_snp_spectrum",
                self.spectrum.to_owned(),
            ),
            ("haplotype_weights.tsv", "f5_haplotype_weights", weight_rows),
            ("ld_profile.tsv", "f5_ld_profile", ld_rows),
        ];
        let links = files
            .iter()
            .map(|(name, kind, contents)| {
                fs::write(profile_directory.join(name), contents).expect("profile file");
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
        fs::write(profile_directory.join("profile.metadata.json"), &metadata).expect("metadata");

        let bases = b"ACGT";
        let span =
            usize::try_from(self.block_length * self.block_count * self.slack).expect("span");
        let sequence = (0..span)
            .map(|index| char::from(bases[index % 4]))
            .collect::<String>();
        let fasta = format!(">chr1\n{sequence}\n");
        fs::write(root.join("authoring/references/ref.fa"), &fasta).expect("reference");
        let general = root.join("authoring/general.toml");
        fs::write(
            &general,
            format!(
                "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
                self.seed,
                sha256(fasta.as_bytes())
            ),
        )
        .expect("general");

        let mut subpopulations = String::new();
        for (index, population) in self.populations.iter().enumerate() {
            subpopulations.push_str(&format!(
                "\n[[subpopulations]]\npopulation_id = \"{population}\"\nindividual_count = 1\nfounder_ids = [\"I{}\"]\n",
                index + 1
            ));
        }
        fs::write(
            scenario_directory.join("scenario.toml"),
            format!(
                "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = {}\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = {}\nmax_swap_proposals = {}\nmax_haplotype_cells = {}\nblock_count = {}\n{subpopulations}",
                self.populations.len(),
                sha256(metadata.as_bytes()),
                self.max_haplotype_attempts,
                self.max_swap_proposals,
                self.max_haplotype_cells,
                self.block_count,
            ),
        )
        .expect("scenario");
        general
    }
}

struct Chain {
    config: EffectiveF5Config,
    sealed: SealedF5Admission,
    execution: ExecutionId,
    layout: F5PhysicalLayout,
    root: PathBuf,
}

/// Cada cadena retira su fixture al terminar.
impl Drop for Chain {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn build_chain(root: &Path, fixture: &Fixture<'_>) -> Chain {
    let general = fixture.write(root);
    let config = load_and_validate_f5_config(&general).expect("configuración F5 admitida");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("admisión sellada");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("solicitud de layout"),
    )
    .expect("layout físico");
    Chain {
        config,
        sealed,
        execution,
        layout,
        root: root.to_path_buf(),
    }
}

fn search(chain: &Chain) -> Result<F5HaplotypeSearch, F5SearchFailure> {
    search_f5_haplotypes(
        &chain.config,
        &chain.sealed,
        &chain.execution,
        &chain.layout,
    )
}

// ---------------------------------------------------------------------------
// Oráculo racional exacto sobre la matriz publicada
// ---------------------------------------------------------------------------

/// Descompone un `f64` finito y positivo en `numerador / 2^exponente` exacto.
fn dyadic(value: f64) -> (u128, u32) {
    assert!(value.is_finite() && value > 0.0, "peso {value} inválido");
    let bits = value.to_bits();
    let raw_exponent = i32::try_from((bits >> 52) & 0x7ff).expect("exponente");
    let mantissa = u128::from(bits & ((1_u64 << 52) - 1));
    if raw_exponent == 0 {
        (mantissa, 1074)
    } else {
        let significand = mantissa | (1_u128 << 52);
        let exponent = 1075 - raw_exponent;
        assert!(exponent > 0, "peso demasiado grande para el oráculo");
        (significand, u32::try_from(exponent).expect("exponente"))
    }
}

/// Pesos como numeradores enteros sobre una escala común `2^scale`.
///
/// Los cuatro perfiles de la batería usan pesos diádicos —`1/H` uniforme o
/// `1/2, 1/4, 1/8, 1/8`—, de modo que cada peso es exactamente
/// `numerador / 2^scale` y todo el álgebra de DTC-089 cabe en `u128` sin
/// perder un bit. La escala se reduce a la mínima común para que los productos
/// de cuatro factores sigan cabiendo en 53 bits.
fn scaled_weights(weights: &[f64]) -> (Vec<u128>, u128) {
    let scale = weights
        .iter()
        .map(|weight| dyadic(*weight).1)
        .max()
        .expect("al menos un peso");
    let numerators: Vec<u128> = weights
        .iter()
        .map(|weight| {
            let (significand, exponent) = dyadic(*weight);
            significand << (scale - exponent)
        })
        .collect();
    let unit: u128 = numerators.iter().sum();
    // Reducir por el mayor común divisor mantiene la escala pequeña sin
    // cambiar ningún cociente.
    let mut divisor = unit;
    for numerator in &numerators {
        divisor = gcd(divisor, *numerator);
    }
    assert!(divisor > 0, "los pesos no suman nada");
    (
        numerators.iter().map(|value| value / divisor).collect(),
        unit / divisor,
    )
}

const fn gcd(mut left: u128, mut right: u128) -> u128 {
    while right != 0 {
        let next = left % right;
        left = right;
        right = next;
    }
    left
}

/// Sumas enteras `(W_i, W_j, W_ij)` de una pareja de columnas, en unidades de
/// `1/unit`.
fn weighted_counts(
    block: &F5BlockHaplotypes,
    numerators: &[u128],
    i: usize,
    j: usize,
) -> (u128, u128, u128) {
    let matrix = block.matrix();
    let mut w_i = 0_u128;
    let mut w_j = 0_u128;
    let mut w_ij = 0_u128;
    for (haplotype, weight) in numerators.iter().enumerate() {
        let alt_i = matrix.allele(haplotype, i).expect("alelo i");
        let alt_j = matrix.allele(haplotype, j).expect("alelo j");
        w_i += u128::from(alt_i) * weight;
        w_j += u128::from(alt_j) * weight;
        w_ij += u128::from(alt_i && alt_j) * weight;
    }
    (w_i, w_j, w_ij)
}

/// `(D, r²)` exactos para pesos **diádicos** cualesquiera.
///
/// Con `q_h = w_h / U` y `U = Σ w_h`: `D = (U*W_ij - W_i*W_j) / U²` y
/// `r² = (U*W_ij - W_i*W_j)² / [W_i(U-W_i) * W_j(U-W_j)]`, donde numerador y
/// denominador comparten `U⁴` y se cancelan. Todos los enteros caben en 53
/// bits con los perfiles de la batería, de modo que la comparación es una
/// igualdad y no una tolerancia. No llama a ninguna función del productor.
fn rational_pair(
    block: &F5BlockHaplotypes,
    numerators: &[u128],
    unit: u128,
    i: usize,
    j: usize,
) -> (f64, Option<f64>) {
    let (w_i, w_j, w_ij) = weighted_counts(block, numerators, i, j);
    let numerator =
        i128::try_from(unit * w_ij).expect("U*W_ij") - i128::try_from(w_i * w_j).expect("W_i*W_j");
    let denominator = w_i * (unit - w_i) * w_j * (unit - w_j);
    let r_squared = (denominator > 0).then(|| {
        let squared = numerator * numerator;
        (squared as f64 / denominator as f64).min(1.0)
    });
    (numerator as f64 / (unit * unit) as f64, r_squared)
}

/// Suma de Neumaier escrita desde su definición.
fn neumaier(values: &[f64]) -> f64 {
    let mut sum = 0.0_f64;
    let mut compensation = 0.0_f64;
    for value in values {
        let next = sum + value;
        compensation += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + compensation
}

/// Los siete campos de la tupla recalculados desde la matriz publicada, sin
/// llamar a ninguna función del productor.
struct OracleTuple {
    duplicate_rows: u64,
    bins_outside: u64,
    signs_failed: u64,
    max_r2_error: f64,
    sum_r2_error: f64,
    max_marginal_error: f64,
    evaluable_pairs: u64,
}

fn oracle_tuple(chain: &Chain, block: &F5BlockHaplotypes) -> OracleTuple {
    let site_count = block.matrix().site_count();
    let haplotypes = block.matrix().haplotype_count();
    let (numerators, unit) = scaled_weights(chain.config.profile().haplotype_weights());
    let layout_block = chain
        .layout
        .blocks()
        .iter()
        .find(|candidate| candidate.block_index() == block.block_index())
        .expect("bloque del layout");
    let bins = chain
        .config
        .profile()
        .ld_bins(block.population_id().as_str())
        .expect("bins del perfil");

    let mut rows: Vec<String> = Vec::with_capacity(haplotypes);
    for haplotype in 0..haplotypes {
        rows.push(
            (0..site_count)
                .map(|site| {
                    if block.matrix().allele(haplotype, site).expect("alelo") {
                        '1'
                    } else {
                        '0'
                    }
                })
                .collect(),
        );
    }
    let distinct: BTreeSet<&String> = rows.iter().collect();
    let duplicate_rows = (rows.len() - distinct.len()) as u64;

    let mut per_bin: Vec<Vec<f64>> = vec![Vec::new(); bins.len()];
    let mut signs_failed = 0_u64;
    let mut evaluable_pairs = 0_u64;
    for i in 0..site_count {
        for j in (i + 1)..site_count {
            let (d, r_squared) = rational_pair(block, &numerators, unit, i, j);
            let Some(value) = r_squared else {
                continue;
            };
            evaluable_pairs += 1;
            let matched = block.columns()[i].polarity() == block.columns()[j].polarity();
            if !(if matched { d > 0.0 } else { d < 0.0 }) {
                signs_failed += 1;
            }
            let distance = layout_block.sites()[j].position_1based()
                - layout_block.sites()[i].position_1based();
            let bin = bins
                .iter()
                .position(|bin| {
                    bin.distance_start_bp() <= distance && distance <= bin.distance_end_bp()
                })
                .expect("bin de la distancia");
            per_bin[bin].push(value);
        }
    }

    let mut bins_outside = 0_u64;
    let mut max_r2_error = 0.0_f64;
    let mut errors = Vec::new();
    for (definition, values) in bins.iter().zip(&per_bin) {
        if values.is_empty() {
            continue;
        }
        let mean = neumaier(values) / values.len() as f64;
        let error = (mean - definition.target_mean_r2()).abs();
        if error > definition.max_abs_error() {
            bins_outside += 1;
        }
        if error > max_r2_error {
            max_r2_error = error;
        }
        errors.push(error);
    }

    let mut max_marginal_error = 0.0_f64;
    for (site, column) in block.columns().iter().enumerate() {
        let mut realized_numerator = 0_u128;
        for (haplotype, weight) in numerators.iter().enumerate() {
            if block.matrix().allele(haplotype, site).expect("alelo") {
                realized_numerator += weight;
            }
        }
        let realized = realized_numerator as f64 / unit as f64;
        let error = (realized - column.target_alt_af()).abs();
        if error > max_marginal_error {
            max_marginal_error = error;
        }
    }

    OracleTuple {
        duplicate_rows,
        bins_outside,
        signs_failed,
        max_r2_error,
        sum_r2_error: neumaier(&errors),
        max_marginal_error,
        evaluable_pairs,
    }
}

/// Visitas de fila de las **tres** sumas canónicas de un intento.
///
/// - siembra de la puntuación: `H*S` marginales y `H*S*(S-1)/2` conjuntas;
/// - recomputación de las columnas publicadas: `H*S`, más `H` del recuento de
///   filas duplicadas;
/// - evaluación LD de cierre: otra vez `H*S + H*S*(S-1)/2`.
///
/// La suma se simplifica a `H*(S+1)²`, y la prueba escribe las dos formas para
/// que un cambio de fase se vea en el desglose y no solo en el total.
fn canonical_row_visits(haplotypes: u64, sites: u64) -> u64 {
    let pairs = sites * (sites - 1) / 2;
    let seed = haplotypes * (sites + pairs);
    let refresh = haplotypes * sites + haplotypes;
    let closing = haplotypes * (sites + pairs);
    let total = seed + refresh + closing;
    assert_eq!(total, haplotypes * (sites + 1) * (sites + 1));
    total
}

fn records_of<'a>(
    records: &'a [F5AttemptRecord],
    block: &F5BlockHaplotypes,
) -> Vec<&'a F5AttemptRecord> {
    records
        .iter()
        .filter(|record| {
            record.population_id() == block.population_id()
                && record.block_index() == block.block_index()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Fixtures nombrados
// ---------------------------------------------------------------------------

/// Búsqueda que termina: `H=2`, `S=3` y bins que piden acoplamiento perfecto.
///
/// Con dos haplotipos `r²` solo puede valer 0 o 1, de modo que el objetivo
/// `1,0` es alcanzable y el bucle acaba decidiéndose por los signos. Un bloque
/// termina antes de la primera propuesta y el otro necesita intercambios.
fn terminating_fixture<'a>() -> Fixture<'a> {
    Fixture {
        weights: vec![0.5, 0.5],
        block_length: 3,
        bins: vec![(1, 1, 1.0, 0.05), (2, 2, 1.0, 0.05)],
        max_haplotype_attempts: 32,
        max_swap_proposals: 256,
        ..Fixture::default()
    }
}

/// La misma búsqueda sobre ocho bloques×subpoblación.
fn wide_terminating_fixture<'a>() -> Fixture<'a> {
    Fixture {
        populations: &["P1", "P2"],
        block_count: 4,
        ..terminating_fixture()
    }
}

/// Dieciséis unidades bloque×subpoblación con la misma geometría.
///
/// Es el fixture de la regresión B-1: cada unidad termina por su cuenta en muy
/// pocos intentos, pero un intento global todo-o-nada exigiría que las
/// dieciséis coincidieran en el mismo índice y agotaría el límite.
fn sixteen_unit_fixture<'a>() -> Fixture<'a> {
    Fixture {
        populations: &["P1", "P2"],
        block_count: 8,
        slack: 3,
        ..terminating_fixture()
    }
}

/// Búsqueda que se agota: `H=4`, `S=4` y una media objetivo que la aritmética
/// exacta de cuatro haplotipos no puede alcanzar dentro de `0,05`.
fn exhausting_fixture<'a>() -> Fixture<'a> {
    Fixture::default()
}

/// Espectro de dos bandas degeneradas: `0,1` deja la columna constante con
/// `H=2` —el prefijo `0,5` dista más de `0,1` que el conjunto vacío— y `0,6`
/// la deja variable. Un bloque con una sola columna variable no tiene ninguna
/// pareja evaluable y, aun así, tiene cero filas duplicadas.
const MIXED_SPECTRUM: &str = "alt_af_min\talt_af_max\tweight\n0.1\t0.1\t0.5\n0.6\t0.6\t0.5\n";

fn degenerate_fixture<'a>() -> Fixture<'a> {
    Fixture {
        populations: &["P1", "P2"],
        block_count: 4,
        weights: vec![0.5, 0.5],
        block_length: 3,
        spectrum: MIXED_SPECTRUM,
        bins: vec![(1, 1, 1.0, 0.05), (2, 2, 1.0, 0.05)],
        max_haplotype_attempts: 8,
        max_swap_proposals: 64,
        ..Fixture::default()
    }
}

// ---------------------------------------------------------------------------
// B-1. El alcance de los intentos es la unidad bloque×subpoblación
// ---------------------------------------------------------------------------

/// Dieciséis unidades resolubles se resuelven.
///
/// Cada unidad termina por su cuenta en pocos intentos; un intento global
/// todo-o-nada exigiría que las dieciséis lo hicieran en el **mismo** índice,
/// cuya probabilidad decae exponencialmente, y declararía `GPS-DOMAIN-104`
/// sobre un problema resuelto. DTC-089 cuenta los intentos por
/// bloque×subpoblación y aquí se comprueba que así ocurre.
#[test]
fn sixteen_independent_units_are_all_resolved() {
    let root = temporary("sixteen");
    let chain = build_chain(&root, &sixteen_unit_fixture());
    let success = search(&chain).expect("dieciséis unidades resolubles");
    assert_eq!(
        success.blocks().len(),
        16,
        "no se resolvieron las 16 unidades"
    );
    assert_eq!(success.ld_evaluation().len(), 16);
    for evaluation in success.ld_evaluation() {
        assert!(evaluation.is_valid());
    }
    // Cada unidad para en su primer intento terminado y no ejecuta ninguno
    // después: `f5-attempts.tsv` no puede contener un intento posterior.
    for block in success.blocks() {
        let records = records_of(success.records(), block);
        let winner = records
            .iter()
            .find(|record| record.selected_best())
            .expect("una fila selected_best por unidad");
        assert!(winner.terminated(), "la unidad publicada no terminó");
        assert_eq!(
            records.len() as u64,
            winner.attempt_index() + 1,
            "la unidad ejecutó intentos posteriores a su primer éxito"
        );
        for record in &records {
            if record.attempt_index() < winner.attempt_index() {
                assert!(
                    !record.terminated(),
                    "hubo un intento terminado anterior al elegido"
                );
            }
        }
    }
    // `attempt_summary`: el total suma los intentos de todas las unidades y el
    // máximo índice es el de la unidad que más necesitó.
    let expected_total: u64 = success.records().len() as u64;
    assert_eq!(success.haplotype_attempts_total(), expected_total);
    let expected_max = success
        .records()
        .iter()
        .filter(|record| record.selected_best())
        .map(F5AttemptRecord::attempt_index)
        .max()
        .expect("registros");
    assert_eq!(success.max_attempt_index(), expected_max);
    assert!(
        expected_max > 0,
        "el fixture no ejercita ninguna unidad que necesite más de un intento"
    );
}

/// El constructor por bloque produce exactamente el bloque del intento entero.
///
/// Es la premisa que permite contar los intentos por unidad: los ámbitos
/// `initial_column` y `polarity` ya llevan bloque y subpoblación, de modo que
/// reconstruir una unidad no depende de las demás.
#[test]
fn the_per_block_builder_reproduces_the_whole_attempt() {
    let root = temporary("per-block");
    let chain = build_chain(&root, &wide_terminating_fixture());
    let marginals = build_f5_marginals(
        &chain.config,
        &chain.sealed,
        &chain.execution,
        &chain.layout,
    )
    .expect("tabla de marginales");
    for attempt_index in 0..3 {
        let whole = build_f5_haplotype_attempt(&chain.config, &marginals, attempt_index)
            .expect("intento completo");
        for (ordinal, block) in whole.blocks().iter().enumerate() {
            let alone = build_f5_block_attempt(&chain.config, &marginals, ordinal, attempt_index)
                .expect("bloque suelto");
            assert_eq!(
                &alone, block,
                "el bloque {ordinal} difiere del intento entero"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Oráculo independiente
// ---------------------------------------------------------------------------

/// La tupla publicada coincide **exactamente** con el oráculo racional.
///
/// Es la red contra la divergencia silenciosa entre la puntuación incremental
/// y la evaluación canónica. La puntuación decide qué propuestas se aceptan,
/// pero lo que se publica se recalcula desde la matriz final; el oráculo lo
/// verifica con enteros de 128 bits y sin llamar a ninguna función del
/// productor, incluidas las marginales realizadas de cada columna. Se ejecuta
/// con pesos uniformes y con pesos diádicos no uniformes.
#[test]
fn the_published_tuple_matches_an_independent_rational_oracle() {
    for (label, fixture) in [
        ("uniform", wide_terminating_fixture()),
        (
            "dyadic",
            Fixture {
                weights: vec![0.5, 0.25, 0.125, 0.125],
                block_length: 3,
                bins: vec![(1, 1, 1.0, 0.05), (2, 2, 1.0, 0.05)],
                max_haplotype_attempts: 64,
                max_swap_proposals: 256,
                ..Fixture::default()
            },
        ),
    ] {
        let root = temporary(label);
        let chain = build_chain(&root, &fixture);
        let outcome = search(&chain);
        let (blocks, records): (Vec<&F5BlockHaplotypes>, &[F5AttemptRecord]) = match &outcome {
            Ok(success) => (success.blocks().iter().collect(), success.records()),
            // Con pesos diádicos la geometría puede no resolverse; entonces el
            // oráculo no tiene matriz publicada y la prueba solo exige que el
            // desenlace sea un agotamiento legítimo.
            Err(failure) => {
                assert_eq!(failure.code(), "GPS-DOMAIN-104", "{label}");
                continue;
            }
        };
        let mut accepted_total = 0_u64;
        for block in blocks {
            let oracle = oracle_tuple(&chain, block);
            let selected = records_of(records, block)
                .into_iter()
                .find(|record| record.selected_best())
                .expect("una fila selected_best por bloque×población");
            accepted_total += selected.accepted();
            assert_eq!(
                selected.duplicate_rows(),
                oracle.duplicate_rows,
                "{label}: duplicados"
            );
            assert_eq!(
                selected.bins_outside(),
                oracle.bins_outside,
                "{label}: bins fuera"
            );
            assert_eq!(
                selected.signs_failed(),
                oracle.signs_failed,
                "{label}: signos"
            );
            assert_eq!(
                selected.max_r2_error(),
                oracle.max_r2_error,
                "{label}: max_r2_error"
            );
            assert_eq!(
                selected.sum_r2_error(),
                oracle.sum_r2_error,
                "{label}: sum_r2_error"
            );
            assert_eq!(
                selected.max_marginal_error(),
                oracle.max_marginal_error,
                "{label}: max_marginal_error"
            );
            assert_eq!(
                selected.evaluable_pairs(),
                oracle.evaluable_pairs,
                "{label}: parejas evaluables"
            );

            // Las columnas publicadas volvieron a ser la suma canónica sobre
            // las H filas de la matriz **final**, no la de la construida.
            let (numerators, unit) = scaled_weights(chain.config.profile().haplotype_weights());
            for (site, column) in block.columns().iter().enumerate() {
                let mut realized = 0_u128;
                for (haplotype, weight) in numerators.iter().enumerate() {
                    if block.matrix().allele(haplotype, site).expect("alelo") {
                        realized += weight;
                    }
                }
                assert_eq!(
                    column.realized_alt_af(),
                    realized as f64 / unit as f64,
                    "{label}: la columna {site} publica una marginal que su matriz no realiza"
                );
            }
        }
        if label == "uniform" {
            assert!(
                accepted_total > 0,
                "el oráculo no llegó a ver ninguna matriz intercambiada"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Propuesta no candidata, mejora y empate rechazado
// ---------------------------------------------------------------------------

/// Una propuesta no candidata consume RNG y no cambia nada.
#[test]
fn most_proposals_are_discarded_without_applying_the_swap() {
    let root = temporary("not-candidate");
    let chain = build_chain(&root, &exhausting_fixture());
    let failure = search(&chain).expect_err("la media objetivo no es alcanzable");
    assert!(
        failure
            .records()
            .iter()
            .any(|record| record.candidate_proposals() < record.proposals()),
        "ninguna propuesta se descartó por no candidata"
    );
    for record in failure.records() {
        assert_eq!(
            record.proposals(),
            64,
            "el presupuesto no se consumió entero"
        );
        assert!(record.candidate_proposals() <= record.proposals());
        assert!(
            record.accepted() <= record.candidate_proposals(),
            "se aceptó una propuesta que no era candidata"
        );
    }
}

/// La búsqueda acepta propuestas y ninguna deja la tupla peor que la encontró.
#[test]
fn accepted_proposals_never_leave_the_tuple_worse_than_it_started() {
    let root = temporary("improve");
    let chain = build_chain(&root, &exhausting_fixture());
    let failure = search(&chain).expect_err("la media objetivo no es alcanzable");
    let accepted: u64 = failure
        .records()
        .iter()
        .map(F5AttemptRecord::accepted)
        .sum();
    assert!(accepted > 0, "ninguna propuesta mejoró la tupla");
    assert!(
        failure
            .records()
            .iter()
            .any(|record| record.accepted() > 0 && record.duplicate_rows() == 0)
    );
}

/// Un bloque ya terminado no plantea ninguna propuesta.
#[test]
fn a_terminated_block_consumes_no_proposal() {
    let root = temporary("early");
    let chain = build_chain(&root, &terminating_fixture());
    let success = search(&chain).expect("búsqueda con éxito");
    assert!(
        success
            .records()
            .iter()
            .any(|record| record.terminated() && record.proposals() == 0),
        "ningún bloque terminó antes de la primera propuesta"
    );
}

/// Parar antes de agotar el presupuesto solo puede significar terminación.
///
/// Es la única señal observable desde fuera de que la puntuación incremental
/// se desvió hasta creerse terminada: si el bucle saliera por un veredicto que
/// el canon no confirma, quedaría un registro con `proposals` por debajo del
/// presupuesto y `terminated == false`.
#[test]
fn stopping_before_the_budget_can_only_mean_termination() {
    for (label, fixture, budget) in [
        ("wide", wide_terminating_fixture(), 256_u64),
        ("narrow", terminating_fixture(), 256),
        ("exhausting", exhausting_fixture(), 64),
        ("degenerate", degenerate_fixture(), 64),
    ] {
        let root = temporary(label);
        let chain = build_chain(&root, &fixture);
        let records = match search(&chain) {
            Ok(success) => success.records().to_vec(),
            Err(failure) => failure.records().to_vec(),
        };
        assert!(!records.is_empty(), "{label}: sin registros");
        for record in &records {
            if record.proposals() < budget {
                assert!(
                    record.terminated(),
                    "{label}: el bloque {} del intento {} paró con {} propuestas de {budget} sin terminar",
                    record.block_index(),
                    record.attempt_index(),
                    record.proposals()
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Éxito temprano y número exacto de intentos
// ---------------------------------------------------------------------------

/// Cada unidad ejecuta exactamente los intentos que necesita.
#[test]
fn every_unit_runs_exactly_the_attempts_it_needs() {
    let root = temporary("attempts");
    let chain = build_chain(&root, &wide_terminating_fixture());
    let success = search(&chain).expect("búsqueda con éxito");
    for block in success.blocks() {
        let records = records_of(success.records(), block);
        let selected: Vec<&&F5AttemptRecord> = records
            .iter()
            .filter(|record| record.selected_best())
            .collect();
        assert_eq!(selected.len(), 1, "selected_best no es único por unidad");
        assert!(selected[0].terminated());
        assert_eq!(
            records.len() as u64,
            selected[0].attempt_index() + 1,
            "la unidad siguió intentando tras terminar"
        );
        let highest = records
            .iter()
            .map(|record| record.attempt_index())
            .max()
            .expect("registros");
        assert_eq!(highest, selected[0].attempt_index());
    }
    assert_eq!(
        success.records().len() as u64,
        success.haplotype_attempts_total()
    );
}

// ---------------------------------------------------------------------------
// Agotamiento, sin relajación ni falsa imposibilidad
// ---------------------------------------------------------------------------

/// Agotar el presupuesto devuelve `GPS-DOMAIN-104`, conserva el mejor intento
/// y **no afirma imposibilidad**.
#[test]
fn exhaustion_keeps_the_best_attempt_without_claiming_impossibility() {
    let root = temporary("exhausted");
    let chain = build_chain(&root, &exhausting_fixture());
    let failure = search(&chain).expect_err("la media objetivo no es alcanzable");
    assert_eq!(failure.code(), "GPS-DOMAIN-104");
    let context = failure.context();
    assert!(context.contains_key("population_id"));
    assert!(context.contains_key("block_id"));
    assert!(
        context["detail"].contains("best_attempt_index="),
        "el diagnóstico no conserva el mejor intento: {}",
        context["detail"]
    );

    // La primera unidad agotó sus ocho intentos y la búsqueda se detuvo ahí:
    // gastar el presupuesto de las demás no añadiría diagnóstico.
    let attempts: BTreeSet<u64> = failure
        .records()
        .iter()
        .map(F5AttemptRecord::attempt_index)
        .collect();
    assert_eq!(attempts, BTreeSet::from([0, 1, 2, 3, 4, 5, 6, 7]));
    let units: BTreeSet<u64> = failure
        .records()
        .iter()
        .map(F5AttemptRecord::block_index)
        .collect();
    assert_eq!(units.len(), 1, "la búsqueda siguió tras la unidad agotada");
    // El problema **no** era imposible: cada intento tuvo parejas evaluables,
    // de modo que `GPS-DOMAIN-102` habría sido una declaración falsa.
    assert!(
        failure
            .records()
            .iter()
            .all(|record| record.evaluable_pairs() > 0)
    );
    assert_eq!(
        failure
            .records()
            .iter()
            .filter(|record| record.selected_best())
            .count(),
        1,
        "la unidad agotada no tiene un único mejor intento"
    );
}

/// Ningún fallo reduce H o S, cambia el perfil o la tolerancia, fusiona bins
/// ni acepta menos bloques.
#[test]
fn exhaustion_does_not_relax_the_problem() {
    let root = temporary("no-relax");
    let fixture = Fixture {
        populations: &["P1", "P2"],
        block_count: 3,
        bins: vec![(1, 1, 1.0, 0.05), (2, 2, 1.0, 0.05), (3, 3, 1.0, 0.05)],
        max_haplotype_attempts: 2,
        max_swap_proposals: 16,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let failure = search(&chain).expect_err("objetivo inalcanzable");
    assert_eq!(failure.code(), "GPS-DOMAIN-104");

    for population in ["P1", "P2"] {
        let bins = chain
            .config
            .profile()
            .ld_bins(population)
            .expect("bins del perfil");
        assert_eq!(bins.len(), 3, "los bins se fusionaron");
        for bin in bins {
            assert_eq!(bin.target_mean_r2(), 1.0, "el objetivo se relajó");
            assert_eq!(bin.max_abs_error(), 0.05, "la tolerancia se relajó");
        }
    }
    // H y S no se redujeron: `row_visits = H*(S + S*(S-1)/2)` lo acredita.
    let haplotypes = chain.config.profile().haplotype_weights().len() as u64;
    let sites = 4_u64;
    for record in failure.records() {
        assert_eq!(
            record.row_visits(),
            canonical_row_visits(haplotypes, sites),
            "H o S cambiaron en el intento {}",
            record.attempt_index()
        );
        assert!(record.bins_outside() <= 3);
    }
}

// ---------------------------------------------------------------------------
// Imposibilidad demostrada y límites de dimensión
// ---------------------------------------------------------------------------

/// `GPS-DOMAIN-102` se emite antes de abrir ningún intento y sobre una
/// demostración, no sobre una observación.
#[test]
fn a_fully_fixed_spectrum_proves_no_pair_can_be_evaluable() {
    let root = temporary("impossible");
    let fixture = Fixture {
        spectrum: FIXED_SPECTRUM,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let failure = search(&chain).expect_err("ninguna pareja puede evaluarse");
    assert_eq!(failure.code(), "GPS-DOMAIN-102");
    assert!(
        failure.records().is_empty(),
        "la demostración no debe abrir ningún intento"
    );
    assert!(failure.context().contains_key("population_id"));
}

/// Un solo bloque demostrablemente estéril detiene el escenario.
///
/// REQ-EVO-75 dice «solo falla el bloque×subpoblación sin ninguna pareja
/// evaluable»: el alcance del fallo es el bloque, y basta uno. El diagnóstico
/// nombra su subpoblación y el bloque aparece en el detalle, de modo que la
/// usuaria sabe cuál es sin recorrer el escenario entero.
#[test]
fn one_sterile_block_among_valid_ones_stops_the_scenario() {
    let root = temporary("one-sterile");
    // Una banda degenerada en 0 con peso 1/4: con cuatro sitios por bloque, hay
    // bloques cuyas columnas quedan casi todas fijadas y bloques que no.
    const MOSTLY_FIXED: &str = "alt_af_min\talt_af_max\tweight\n0\t0\t0.75\n0.5\t0.5\t0.25\n";
    let fixture = Fixture {
        block_count: 8,
        slack: 3,
        spectrum: MOSTLY_FIXED,
        bins: vec![(1, 3, 1.0, 0.05)],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let failure = search(&chain).expect_err("algún bloque queda sin columnas variables");
    assert_eq!(failure.code(), "GPS-DOMAIN-102");
    assert!(failure.records().is_empty());
    assert!(
        failure.context()["detail"].contains("columnas no fijadas"),
        "el detalle no nombra el bloque estéril: {}",
        failure.context()["detail"]
    );
}

/// Dos columnas no fijadas bastan: la demostración es la mínima que el
/// contrato permite y no rechaza un bloque que sí puede tener una pareja.
///
/// Con `S=2` la única pareja es `(0,1)` y `2S-3 = 1`: el bucle de `k` de
/// `apply_increments` no da ninguna vuelta y solo se rehace esa pareja. Es el
/// borde exacto de la actualización incremental.
#[test]
fn two_variable_columns_are_enough_to_open_the_search() {
    let root = temporary("two-columns");
    let fixture = Fixture {
        block_length: 2,
        bins: vec![(1, 1, 0.5, 0.05)],
        max_haplotype_attempts: 2,
        max_swap_proposals: 8,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let records = match search(&chain) {
        Ok(success) => success.records().to_vec(),
        Err(failure) => {
            assert_eq!(
                failure.code(),
                "GPS-DOMAIN-104",
                "dos columnas variables no son una imposibilidad demostrada"
            );
            failure.records().to_vec()
        }
    };
    assert!(
        !records.is_empty(),
        "la búsqueda no llegó a abrir un intento"
    );
    for record in &records {
        assert_eq!(
            record.row_visits(),
            canonical_row_visits(4, 2),
            "el bloque de dos sitios no recorrió las tres sumas canónicas"
        );
        assert_eq!(
            record.pair_updates(),
            record.candidate_proposals(),
            "con S=2 cada propuesta candidata rehace exactamente una pareja"
        );
    }
    assert!(
        records
            .iter()
            .any(|record| record.candidate_proposals() > 0),
        "ninguna propuesta llegó a ser candidata con S=2"
    );
}

/// El presupuesto de casillas se hace cumplir **antes** de la búsqueda.
#[test]
fn the_cell_budget_is_enforced_before_the_search_can_run() {
    let root = temporary("cells");
    let fixture = Fixture {
        max_haplotype_cells: 8,
        ..Fixture::default()
    };
    let general = fixture.write(&root);
    let error = load_and_validate_f5_config(&general)
        .expect_err("la admisión debe rechazar el presupuesto de casillas");
    assert_eq!(error.code(), "GPS-RES-101");
    assert!(error.context()["detail"].contains("max_haplotype_cells"));
    let _ = fs::remove_dir_all(&root);
}

// ---------------------------------------------------------------------------
// Complejidad medida
// ---------------------------------------------------------------------------

/// Ninguna propuesta recorre las H filas.
///
/// `row_visits` no es una identidad escrita a mano: el contador vive dentro de
/// `site_marginal` y `joint_marginal`, de modo que cualquier suma canónica lo
/// incrementa, la haga la siembra o la haga una propuesta. Vale exactamente
/// `H*(S + S*(S-1)/2)` sea cual sea el presupuesto. `pair_updates`, contado
/// dentro de `refresh_pair`, vale exactamente `candidate_proposals*(2S-3)` y
/// no depende de H.
#[test]
fn the_measured_cost_separates_the_canonical_pass_from_the_proposals() {
    let sites = 5_u64;
    let mut measured = Vec::new();
    for (label, haplotypes, proposals) in [
        ("h4-q64", 4_usize, 64_u64),
        ("h4-q256", 4, 256),
        ("h8-q256", 8, 256),
    ] {
        let root = temporary(label);
        let fixture = Fixture {
            weights: vec![1.0 / haplotypes as f64; haplotypes],
            block_length: 5,
            bins: vec![(1, 4, 1.0, 0.05)],
            max_haplotype_attempts: 1,
            max_swap_proposals: proposals,
            ..Fixture::default()
        };
        let chain = build_chain(&root, &fixture);
        let failure = search(&chain).expect_err("objetivo inalcanzable");
        let attempt_visits = canonical_row_visits(haplotypes as u64, sites);
        for record in failure.records() {
            assert_eq!(
                record.proposals(),
                proposals,
                "{label}: el presupuesto no se consumió entero"
            );
            assert_eq!(
                record.row_visits(),
                attempt_visits,
                "{label}: una propuesta recorrió las H filas"
            );
            assert_eq!(
                record.pair_updates(),
                record.candidate_proposals() * (2 * sites - 3),
                "{label}: una propuesta actualizó un número de parejas distinto de 2S-3"
            );
        }
        measured.push(attempt_visits);
    }
    // Cuadruplicar el presupuesto no cambia el recorrido canónico de siembra.
    assert_eq!(measured[0], measured[1]);
    // Duplicar H lo duplica.
    assert_eq!(measured[2], 2 * measured[1]);
}

// ---------------------------------------------------------------------------
// Puntuación incremental frente a canon
// ---------------------------------------------------------------------------

/// La puntuación incremental no puede separarse del canon.
///
/// Con pesos **no diádicos** cada intercambio mueve las dos marginales por una
/// cantidad que no es exacta en binario, de modo que el estado incremental es
/// el más expuesto a desviarse. Al cerrar cada intento la búsqueda exige que
/// sus cuatro campos enteros coincidan con la evaluación canónica y que los
/// tres `f64` no disten más de `1e-9`; una divergencia sería
/// `GPS-INTERNAL-101`. La prueba comprueba que el desenlace es legítimo y que
/// hubo propuestas aceptadas que ejercieran la comparación.
#[test]
fn the_incremental_score_never_diverges_from_the_canonical_one() {
    for (label, weights, block_length) in [
        ("non-dyadic", vec![0.3, 0.2, 0.3, 0.2], 6_u64),
        (
            "non-dyadic-eight",
            vec![0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            5,
        ),
        ("dyadic", vec![0.5, 0.25, 0.125, 0.125], 6),
    ] {
        let root = temporary(label);
        let fixture = Fixture {
            populations: &["P1", "P2"],
            weights,
            block_length,
            slack: 2,
            bins: vec![(1, block_length - 1, 1.0, 0.05)],
            max_haplotype_attempts: 6,
            max_swap_proposals: 1_024,
            ..Fixture::default()
        };
        let chain = build_chain(&root, &fixture);
        let records = match search(&chain) {
            Ok(success) => success.records().to_vec(),
            Err(failure) => {
                assert_eq!(
                    failure.code(),
                    "GPS-DOMAIN-104",
                    "{label}: la puntuación incremental se separó del canon"
                );
                failure.records().to_vec()
            }
        };
        let accepted: u64 = records.iter().map(F5AttemptRecord::accepted).sum();
        assert!(
            accepted > 0,
            "{label}: sin propuestas aceptadas la prueba no vigila nada"
        );
    }
}

/// Con pesos no uniformes, un intercambio mueve las marginales.
///
/// `p_hat_i` cambia en `q_b - q_a`, que solo es cero cuando los dos pesos
/// coinciden. Vigila la condición que la revisión de F5.2.2 dejó escrita: sin
/// recomputar canónicamente la columna publicada antes de volver a evaluar, la
/// primera propuesta aceptada haría abortar al evaluador con
/// `GPS-INTERNAL-101` en lugar de terminar en un agotamiento legítimo.
#[test]
fn non_uniform_weights_survive_every_accepted_swap() {
    let root = temporary("non-uniform");
    let fixture = Fixture {
        weights: vec![0.5, 0.25, 0.125, 0.125],
        block_length: 5,
        bins: vec![(1, 4, 1.0, 0.05)],
        max_haplotype_attempts: 4,
        max_swap_proposals: 256,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let failure = search(&chain).expect_err("objetivo inalcanzable");
    assert_eq!(
        failure.code(),
        "GPS-DOMAIN-104",
        "la búsqueda publicó una columna que su matriz ya no realiza"
    );
    let accepted: u64 = failure
        .records()
        .iter()
        .map(F5AttemptRecord::accepted)
        .sum();
    assert!(
        accepted > 0,
        "sin propuestas aceptadas la prueba no vigila nada"
    );
    let tolerance = chain
        .config
        .profile()
        .haplotype_weights()
        .iter()
        .copied()
        .fold(0.0_f64, f64::max);
    assert_eq!(tolerance, 0.5, "la cota marginal activa es max_h(q_h)");
    assert!(
        failure
            .records()
            .iter()
            .any(|record| record.accepted() > 0 && record.max_marginal_error() > 0.0),
        "ningún intento aceptado movió una marginal"
    );
}

// ---------------------------------------------------------------------------
// Invariantes del intercambio
// ---------------------------------------------------------------------------

/// El intercambio candidato conserva los recuentos de fila y de columna.
///
/// `10/01 ↔ 01/10` retira un alelo de una columna y lo devuelve en la otra
/// fila, de modo que ni el número de haplotipos alternativos de una columna ni
/// el peso de una fila cambian; las catorce configuraciones que **no** son
/// candidatas rompen esa invariante. Y una propuesta rechazada deshace su
/// intercambio, así que la matriz final solo puede diferir de la construida en
/// las cuatro casillas de cada propuesta aceptada.
#[test]
fn only_crossed_swaps_are_applied_so_every_row_and_column_count_survives() {
    let root = temporary("counts");
    let chain = build_chain(&root, &wide_terminating_fixture());
    let success = search(&chain).expect("búsqueda con éxito");
    let marginals = build_f5_marginals(
        &chain.config,
        &chain.sealed,
        &chain.execution,
        &chain.layout,
    )
    .expect("tabla de marginales");

    let mut changed_blocks = 0_u32;
    for (ordinal, searched) in success.blocks().iter().enumerate() {
        let winner = records_of(success.records(), searched)
            .into_iter()
            .find(|record| record.selected_best())
            .expect("registro publicado");
        let original =
            build_f5_block_attempt(&chain.config, &marginals, ordinal, winner.attempt_index())
                .expect("bloque reconstruido");
        assert_eq!(
            original.matrix().haplotype_count(),
            searched.matrix().haplotype_count()
        );
        assert_eq!(
            original.matrix().site_count(),
            searched.matrix().site_count()
        );

        let mut differing = 0_u64;
        for haplotype in 0..searched.matrix().haplotype_count() {
            for site in 0..searched.matrix().site_count() {
                if original.matrix().allele(haplotype, site)
                    != searched.matrix().allele(haplotype, site)
                {
                    differing += 1;
                }
            }
        }
        if differing > 0 {
            changed_blocks += 1;
        }
        assert!(
            differing <= 4 * winner.accepted(),
            "la matriz difiere en {differing} casillas y solo se aceptaron {} propuestas",
            winner.accepted()
        );
        for site in 0..searched.matrix().site_count() {
            assert_eq!(
                original.columns()[site].alternate_haplotypes(),
                searched.columns()[site].alternate_haplotypes(),
                "la columna {site} cambió su número de alelos alternativos"
            );
        }
        for haplotype in 0..searched.matrix().haplotype_count() {
            let weight_of = |block: &F5BlockHaplotypes| {
                (0..block.matrix().site_count())
                    .filter(|site| block.matrix().allele(haplotype, *site).expect("alelo"))
                    .count()
            };
            assert_eq!(
                weight_of(&original),
                weight_of(searched),
                "la fila {haplotype} cambió su número de alelos alternativos"
            );
        }
    }
    assert!(
        changed_blocks > 0,
        "ninguna propuesta llegó a modificar una matriz"
    );
}

/// La búsqueda nunca deja los tres campos discretos peor que los encontró.
#[test]
fn the_search_never_degrades_the_discrete_fields_it_started_from() {
    for (label, fixture) in [
        ("exhausting", exhausting_fixture()),
        (
            "non-uniform",
            Fixture {
                populations: &["P1", "P2"],
                weights: vec![0.5, 0.25, 0.125, 0.125],
                block_length: 6,
                slack: 2,
                bins: vec![(1, 5, 1.0, 0.05)],
                max_haplotype_attempts: 6,
                max_swap_proposals: 1_024,
                ..Fixture::default()
            },
        ),
    ] {
        let root = temporary(label);
        let chain = build_chain(&root, &fixture);
        let failure = search(&chain).expect_err("objetivo inalcanzable");
        let marginals = build_f5_marginals(
            &chain.config,
            &chain.sealed,
            &chain.execution,
            &chain.layout,
        )
        .expect("tabla de marginales");
        let mut compared = 0_u32;
        for record in failure.records() {
            if record.accepted() == 0 {
                continue;
            }
            let ordinal = marginals
                .blocks()
                .iter()
                .position(|candidate| {
                    candidate.population_id() == record.population_id()
                        && candidate.block_index() == record.block_index()
                })
                .expect("unidad en la tabla de marginales");
            let block =
                build_f5_block_attempt(&chain.config, &marginals, ordinal, record.attempt_index())
                    .expect("bloque reconstruido");
            let layout_block = chain
                .layout
                .blocks()
                .iter()
                .find(|candidate| candidate.block_index() == record.block_index())
                .expect("bloque del layout");
            let initial =
                evaluate_f5_block_linkage_disequilibrium(&chain.config, layout_block, &block)
                    .expect("evaluación inicial");
            let before = (
                block.duplicate_rows(),
                initial.bins_outside(),
                initial.signs_failed(),
            );
            let after = (
                record.duplicate_rows(),
                record.bins_outside(),
                record.signs_failed(),
            );
            assert!(
                after <= before,
                "{label}: el bloque {} del intento {} pasó de {before:?} a {after:?} tras {} aceptaciones",
                record.block_index(),
                record.attempt_index(),
                record.accepted()
            );
            compared += 1;
        }
        assert!(compared > 0, "{label}: ningún intento aceptó propuestas");
    }
}

/// Un bloque degenerado no termina y no gana la selección.
#[test]
fn a_degenerate_block_neither_terminates_nor_is_selected() {
    let root = temporary("degenerate");
    let chain = build_chain(&root, &degenerate_fixture());
    let records = match search(&chain) {
        Ok(success) => success.records().to_vec(),
        Err(failure) => failure.records().to_vec(),
    };
    let degenerate: Vec<&F5AttemptRecord> = records
        .iter()
        .filter(|record| record.evaluable_pairs() == 0)
        .collect();
    assert!(
        !degenerate.is_empty(),
        "el fixture no produjo ningún bloque sin parejas evaluables"
    );
    // Sin ninguna pareja evaluable, todos los bins son `not_applicable` y
    // cuatro de los siete campos quedan en su mejor valor posible. Es la tupla
    // que haría ganar a un intento degenerado si la validez no se consultara
    // antes.
    for record in &degenerate {
        assert_eq!(record.bins_outside(), 0);
        assert_eq!(record.signs_failed(), 0);
        assert_eq!(record.max_r2_error(), 0.0);
        assert_eq!(record.sum_r2_error(), 0.0);
    }
    for record in &degenerate {
        assert!(
            !record.terminated(),
            "un bloque sin parejas evaluables se dio por terminado"
        );
    }
    for record in &records {
        if !record.selected_best() || record.evaluable_pairs() > 0 {
            continue;
        }
        let has_valid_alternative = records.iter().any(|other| {
            other.population_id() == record.population_id()
                && other.block_index() == record.block_index()
                && other.evaluable_pairs() > 0
        });
        assert!(
            !has_valid_alternative,
            "se eligió un intento degenerado habiendo uno válido"
        );
    }
}

/// La búsqueda es reproducible bajo su semilla y depende de ella.
#[test]
fn the_search_is_reproducible_under_its_seed() {
    const OTHER_SEED: &str = "0f0e0d0c0b0a09080706050403020100f0e0d0c0b0a090807060504030201000";
    let first = {
        let root = temporary("seed-a");
        let chain = build_chain(&root, &exhausting_fixture());
        search(&chain).expect_err("agotamiento").records().to_vec()
    };
    let repeated = {
        let root = temporary("seed-a-bis");
        let chain = build_chain(&root, &exhausting_fixture());
        search(&chain).expect_err("agotamiento").records().to_vec()
    };
    assert_eq!(first, repeated, "la misma semilla dio otro historial");

    let root = temporary("seed-b");
    let chain = build_chain(
        &root,
        &Fixture {
            seed: OTHER_SEED,
            ..exhausting_fixture()
        },
    );
    let other = match search(&chain) {
        Ok(success) => success.records().to_vec(),
        Err(failure) => failure.records().to_vec(),
    };
    assert_ne!(
        first, other,
        "otra semilla recorrió exactamente el mismo camino"
    );
}

/// El ámbito `swap_proposal` de cada intento es el suyo, no el del primero.
///
/// El fixture da **una sola propuesta** por intento, de modo que lo que ocurre
/// en él queda determinado por los cuatro enteros que sortea. La prueba deriva
/// esos cuatro enteros por su cuenta —subsemilla del ámbito, `a in 0..H`,
/// `b in 0..H-1`, `i in 0..S`, `j in 0..S-1`, incrementos y ordenación—, los
/// aplica sobre la matriz que el constructor entrega para **ese**
/// `attempt_index` y comprueba que la candidatura coincide con la que el
/// registro publica. Si la búsqueda abriera el ámbito con un índice constante,
/// todos los intentos sortearían la propuesta del intento cero y la
/// candidatura dejaría de coincidir.
#[test]
fn every_attempt_draws_from_its_own_swap_scope() {
    let root = temporary("swap-scope");
    let fixture = Fixture {
        max_haplotype_attempts: 8,
        max_swap_proposals: 1,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let failure = search(&chain).expect_err("una propuesta por intento no alcanza el objetivo");
    let marginals = build_f5_marginals(
        &chain.config,
        &chain.sealed,
        &chain.execution,
        &chain.layout,
    )
    .expect("tabla de marginales");
    let scenario_id = ScenarioId::parse(chain.config.scenario_id(), "test").expect("escenario");
    let root_seed = RootSeed::parse("seed", chain.config.seed()).expect("semilla");
    let haplotypes = chain.config.profile().haplotype_weights().len() as u64;

    let mut checked = 0_u32;
    for record in failure.records() {
        let ordinal = marginals
            .blocks()
            .iter()
            .position(|candidate| {
                candidate.population_id() == record.population_id()
                    && candidate.block_index() == record.block_index()
            })
            .expect("unidad en la tabla de marginales");
        let block =
            build_f5_block_attempt(&chain.config, &marginals, ordinal, record.attempt_index())
                .expect("bloque reconstruido");
        let sites = block.matrix().site_count() as u64;

        // Los cuatro muestreos literales de DTC-089, derivados aquí.
        let mut rng = derive_f5_subseed(
            root_seed,
            F5RandomDecision::swap_proposal(
                &scenario_id,
                record.population_id(),
                record.block_id(),
                record.attempt_index(),
            ),
        )
        .expect("subsemilla")
        .local_rng();
        let a = f5_uniform_index(&mut rng, haplotypes).expect("a");
        let mut b = f5_uniform_index(&mut rng, haplotypes - 1).expect("b");
        let i = f5_uniform_index(&mut rng, sites).expect("i");
        let mut j = f5_uniform_index(&mut rng, sites - 1).expect("j");
        if b >= a {
            b += 1;
        }
        if j >= i {
            j += 1;
        }
        let (row_a, row_b) = (a.min(b) as usize, a.max(b) as usize);
        let (site_i, site_j) = (i.min(j) as usize, i.max(j) as usize);

        let allele =
            |haplotype: usize, site: usize| block.matrix().allele(haplotype, site).expect("alelo");
        let candidate = allele(row_a, site_i) != allele(row_a, site_j)
            && allele(row_b, site_i) != allele(row_b, site_j)
            && allele(row_a, site_i) != allele(row_b, site_i);
        assert_eq!(
            record.candidate_proposals(),
            u64::from(candidate),
            "el intento {} del bloque {} no sorteó la propuesta de su propio ámbito",
            record.attempt_index(),
            record.block_index()
        );
        checked += 1;
    }
    assert!(checked >= 8, "la prueba solo comprobó {checked} intentos");
}

/// El coste del bucle **no depende del número de bins del perfil**.
///
/// REQ-EVO-75 exige bins contiguos que cubran de 1 a `Lmax-1`, así que un
/// perfil de resolución por par de bases es legal y admitido; y un `gap_min`
/// grande frente a `Lmax` deja `S_bound` pequeño. Si cada propuesta recorriese
/// los `B` bins —o respaldase los dos vectores de bin enteros—, el coste
/// crecería como `B/S_bound` sin cota: con 255 bins y dos sitios por bloque
/// llegó a ser 62 veces el presupuesto sellado de la admisión.
///
/// La prueba fija la geometría y solo cambia el número de bins: con un bin y
/// con `L-1` bins de un par de bases, el recuento de visitas elementales del
/// bucle debe ser **el mismo**. Y comprueba además que se mantiene dentro de un
/// múltiplo pequeño y declarado del presupuesto `H*S_bound² + Q*S_bound` que
/// `f5-resource-accounting-v7` concede por bloque×subpoblación×intento.
#[test]
fn the_proposal_loop_cost_does_not_depend_on_the_bin_count() {
    /// Cota declarada de la **forma** del coste por propuesta: una constante
    /// —los cuatro sorteos, los cuatro bits del intercambio y las cuatro
    /// operaciones del índice de filas— más un término lineal en `S`. Es lo
    /// que DTC-089 fija al decir que las actualizaciones 2×2 mantienen las
    /// estadísticas en `O(S)`. Las «unidades de operación» de v7 son una
    /// medida más gruesa que estas visitas elementales, así que la cota se
    /// escribe con sus dos constantes a la vista en lugar de como un múltiplo
    /// de `Q*S_bound`.
    const CONSTANT_PER_PROPOSAL: u64 = 24;
    const PER_SITE_PER_PROPOSAL: u64 = 8;
    const PROPOSALS: u64 = 1_024;

    let mut measured = Vec::new();
    for (label, block_length, gap, one_bin_per_base) in [
        ("L256-g128-1bin", 256_u64, 128_u64, false),
        ("L256-g128-255bins", 256, 128, true),
        ("L64-g32-1bin", 64, 32, false),
        ("L64-g32-63bins", 64, 32, true),
    ] {
        let bins: Vec<(u64, u64, f64, f64)> = if one_bin_per_base {
            (1..block_length)
                .map(|distance| (distance, distance, 1.0, 0.05))
                .collect()
        } else {
            vec![(1, block_length - 1, 1.0, 0.05)]
        };
        let bin_count = bins.len();
        let root = temporary(label);
        let fixture = Fixture {
            block_count: 1,
            block_length,
            gap,
            bins,
            max_haplotype_attempts: 2,
            max_swap_proposals: PROPOSALS,
            ..Fixture::default()
        };
        let chain = build_chain(&root, &fixture);
        let records = match search(&chain) {
            Ok(success) => success.records().to_vec(),
            Err(failure) => failure.records().to_vec(),
        };
        assert!(!records.is_empty(), "{label}: sin registros");
        let site_bound = chain.sealed.site_bound_per_block();
        let haplotypes = chain.config.profile().haplotype_weights().len() as u64;
        let coarse_budget = haplotypes * site_bound * site_bound + PROPOSALS * site_bound;
        for record in &records {
            let bound =
                record.proposals() * (CONSTANT_PER_PROPOSAL + PER_SITE_PER_PROPOSAL * site_bound);
            assert!(
                record.loop_operations() <= bound,
                "{label}: el bucle costó {} visitas para {} propuestas, por encima de {bound} (presupuesto grueso de v7: {coarse_budget})",
                record.loop_operations(),
                record.proposals()
            );
        }
        let total: u64 = records.iter().map(F5AttemptRecord::loop_operations).sum();
        measured.push((label, bin_count, total, records.len()));
    }

    // Misma geometría, un bin frente a `L-1`: el coste debe coincidir.
    assert_eq!(
        measured[0].3, measured[1].3,
        "los intentos no son comparables"
    );
    assert_eq!(
        measured[0].2, measured[1].2,
        "con {} bins el bucle costó {} y con {} bins costó {}",
        measured[0].1, measured[0].2, measured[1].1, measured[1].2
    );
    assert_eq!(measured[2].3, measured[3].3);
    assert_eq!(
        measured[2].2, measured[3].2,
        "con {} bins el bucle costó {} y con {} bins costó {}",
        measured[2].1, measured[2].2, measured[3].1, measured[3].2
    );
    assert!(
        measured[1].1 >= 255 && measured[3].1 >= 63,
        "los fixtures de muchos bins no se construyeron"
    );
}

/// El desempate final por `attempt_index`, alcanzado por la vía productiva.
///
/// Con `H=2`, `S=2`, espectro degenerado en `0,5` y un objetivo de bin
/// inalcanzable, varios intentos producen **la misma matriz ganadora** y
/// empatan en los siete campos de la tupla: el único árbitro que queda es el
/// menor `attempt_index`, y es quien decide qué fila lleva `selected_best` en
/// `f5-attempts.tsv`. Hasta aquí solo había prueba unitaria del comparador.
#[test]
fn the_attempt_index_tiebreak_is_reachable_end_to_end() {
    const HALF_SPECTRUM: &str = "alt_af_min\talt_af_max\tweight\n0.5\t0.5\t1\n";
    let root = temporary("tiebreak");
    let fixture = Fixture {
        weights: vec![0.5, 0.5],
        block_count: 1,
        block_length: 2,
        spectrum: HALF_SPECTRUM,
        bins: vec![(1, 1, 0.25, 0.05)],
        max_haplotype_attempts: 40,
        max_swap_proposals: 8,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let failure = search(&chain).expect_err("con r² en {0,1} el objetivo 0,25±0,05 no se alcanza");
    assert_eq!(failure.code(), "GPS-DOMAIN-104");

    let selected = failure
        .records()
        .iter()
        .find(|record| record.selected_best())
        .expect("un mejor intento");
    // Intentos que empatan con el elegido en los siete campos de la tupla.
    let tied: Vec<&F5AttemptRecord> = failure
        .records()
        .iter()
        .filter(|record| {
            record.duplicate_rows() == selected.duplicate_rows()
                && record.bins_outside() == selected.bins_outside()
                && record.signs_failed() == selected.signs_failed()
                && record.max_r2_error() == selected.max_r2_error()
                && record.sum_r2_error() == selected.sum_r2_error()
                && record.max_marginal_error() == selected.max_marginal_error()
                && record.matrix_bytes_sha256() == selected.matrix_bytes_sha256()
                && record.evaluable_pairs() == selected.evaluable_pairs()
                && record.terminated() == selected.terminated()
        })
        .collect();
    assert!(
        tied.len() >= 2,
        "el fixture no produjo ningún empate completo: solo {} intento(s)",
        tied.len()
    );
    let lowest = tied
        .iter()
        .map(|record| record.attempt_index())
        .min()
        .expect("empates");
    assert_eq!(
        selected.attempt_index(),
        lowest,
        "el desempate no eligió el menor attempt_index entre {} intentos idénticos",
        tied.len()
    );
}

/// Un bin sin ninguna pareja evaluable llega a `not_applicable` a través de la
/// frontera productiva.
///
/// `f5-ld-evaluation.tsv` distingue `evaluated` de `not_applicable`, y el
/// segundo no es un error ni se sustituye por otro bin. Hasta aquí solo se
/// comprobaba sobre el evaluador suelto.
#[test]
fn an_empty_bin_reaches_not_applicable_through_the_search() {
    let root = temporary("not-applicable");
    // Bins de un par de bases sobre un bloque de cuatro bases cuyos dos sitios
    // distan tres: solo el último bin puede contener parejas, y los otros dos
    // quedan sin ninguna. El perfil sigue cubriendo `1..=Lmax-1`, como exige
    // REQ-EVO-75.
    let fixture = Fixture {
        block_count: 1,
        block_length: 4,
        gap: 3,
        bins: vec![(1, 1, 0.5, 0.05), (2, 2, 0.5, 0.05), (3, 3, 1.0, 0.05)],
        max_haplotype_attempts: 2,
        max_swap_proposals: 8,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let outcome = search(&chain);
    let evaluation = match &outcome {
        Ok(success) => success.ld_evaluation().to_vec(),
        Err(failure) => {
            assert_eq!(failure.code(), "GPS-DOMAIN-104");
            let _ = fs::remove_dir_all(&root);
            return;
        }
    };
    let block = &evaluation[0];
    assert_eq!(block.bins().len(), 3, "los bins del perfil se fusionaron");
    assert!(
        block.bins()[2].is_evaluated(),
        "el bin de distancia 3 no se evaluó"
    );
    for bin in &block.bins()[..2] {
        assert!(
            !bin.is_evaluated(),
            "un bin sin parejas se dio por evaluado"
        );
        assert_eq!(bin.status(), "not_applicable");
        assert_eq!(bin.pair_count(), 0);
        assert_eq!(bin.mean_r2(), None);
        assert_eq!(bin.abs_error(), None);
        assert!(
            bin.is_within_tolerance(),
            "un bin vacío no está fuera de tolerancia"
        );
    }
    assert_eq!(
        block.bins_outside(),
        u64::from(!block.bins()[2].is_within_tolerance()),
        "un bin no aplicable contó como fuera de tolerancia"
    );
}

/// Rehacer el máximo de `r²` recorre los bins **activos**, no los declarados.
///
/// El máximo no admite mantenimiento por incrementos cuando el bin que lo
/// ostentaba baja, así que se rehace; y ese rehacer es el último sitio donde el
/// coste podría volver a depender de `B`. La prueba elige una geometría con
/// pocos sitios y muchos bins —gaps en `[30,40]` sobre bloques de 120 bases,
/// con un bin por par de bases—, comprueba que el rehacer **dispara**, y exige
/// la identidad exacta
///
/// ```text
/// max_rebuild_visits = max_rebuilds * (bins con al menos una pareja)
/// ```
///
/// El número de bins activos lo calcula la prueba por su cuenta, desde el
/// layout y el perfil, sin preguntárselo al productor.
#[test]
fn rebuilding_the_maximum_visits_only_the_active_bins() {
    let root = temporary("rebuild");
    let block_length = 120_u64;
    let fixture = Fixture {
        block_count: 1,
        block_length,
        gap: 30,
        gap_max: Some(40),
        bins: (1..block_length)
            .map(|distance| (distance, distance, 1.0, 0.05))
            .collect(),
        max_haplotype_attempts: 2,
        max_swap_proposals: 1_024,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let records = match search(&chain) {
        Ok(success) => success.records().to_vec(),
        Err(failure) => failure.records().to_vec(),
    };
    assert!(!records.is_empty(), "sin registros");

    // Bins activos, calculados aquí: los que contienen alguna distancia
    // intrabloque. El conjunto no cambia durante la búsqueda porque las
    // posiciones son del layout.
    let layout_block = &chain.layout.blocks()[0];
    let bins = chain
        .config
        .profile()
        .ld_bins("P1")
        .expect("bins del perfil");
    let sites = layout_block.sites();
    let mut active: BTreeSet<usize> = BTreeSet::new();
    for i in 0..sites.len() {
        for j in (i + 1)..sites.len() {
            let distance = sites[j].position_1based() - sites[i].position_1based();
            let bin = bins
                .iter()
                .position(|bin| {
                    bin.distance_start_bp() <= distance && distance <= bin.distance_end_bp()
                })
                .expect("bin de la distancia");
            active.insert(bin);
        }
    }
    assert!(
        active.len() < bins.len(),
        "el fixture no distingue bins activos de declarados: {} de {}",
        active.len(),
        bins.len()
    );
    assert!(
        active.len() > 1,
        "el rehacer debe recorrer más de un bin para que la identidad discrimine"
    );
    assert!(
        sites.len() * (sites.len() - 1) / 2 >= active.len(),
        "los bins activos superan el número de parejas"
    );

    let rebuilds: u64 = records.iter().map(F5AttemptRecord::max_rebuilds).sum();
    assert!(rebuilds > 0, "el rehacer del máximo no llegó a dispararse");
    for record in &records {
        assert_eq!(
            record.max_rebuild_visits(),
            record.max_rebuilds() * active.len() as u64,
            "el rehacer visitó {} bins en {} disparos, y los activos son {}",
            record.max_rebuild_visits(),
            record.max_rebuilds(),
            active.len()
        );
    }
}
