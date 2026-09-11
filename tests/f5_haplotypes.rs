//! F5.2.1. Construcción haplotípica inicial: columnas por permutación y
//! prefijo ponderado, sin enumerar `2^H`.
//!
//! Las pruebas recorren la cadena productiva completa —configuración,
//! admisión sellada, layout, marginales y haplotipos— y contrastan cada
//! columna con un oráculo escrito en enteros `u128`: permutación Fisher–Yates
//! con su rechazo, prefijo ponderado exacto sobre racionales diádicos y regla
//! de empate. El oráculo no llama a ninguna función productiva de muestreo.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::f5_haplotypes::{
    F5BlockHaplotypes, F5HaplotypeAttempt, F5Polarity, build_f5_haplotype_attempt,
};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_marginals::{F5MarginalTable, build_f5_marginals};
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::identity::{ExecutionId, PopulationId, ScenarioId};
use genome_population_simulator::randomness::{F5RandomDecision, RootSeed, derive_f5_subseed};
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

const EXECUTION: &str = "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4";
const SEED: &str = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8";
/// Escala común del oráculo: todo racional se expresa sobre `2^54`.
const TARGET_SCALE: u32 = 54;

/// Espectro de una sola banda continua `[1/4, 3/4]`: la fila no se sortea y el
/// valor consume un único U53.
const CONTINUOUS_SPECTRUM: &str = "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t1\n";

fn sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-haplotypes-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

struct Fixture<'a> {
    populations: &'a [&'a str],
    block_count: u64,
    block_length: u64,
    weights: Vec<f64>,
    spectrum: &'a str,
    max_haplotype_cells: u64,
}

impl Default for Fixture<'_> {
    fn default() -> Self {
        Self {
            populations: &["P1"],
            block_count: 2,
            block_length: 3,
            weights: vec![0.25; 4],
            spectrum: CONTINUOUS_SPECTRUM,
            max_haplotype_cells: 65_536,
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
            // Cada población LD debe terminar exactamente en Lmax-1.
            ld_rows.push_str(&format!(
                "{population}\t1\t{}\t0.5\t0.05\n",
                self.block_length - 1
            ));
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
                "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n".to_owned(),
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
        let span = usize::try_from(self.block_length * self.block_count).expect("span");
        let sequence = (0..span)
            .map(|index| char::from(bases[index % 4]))
            .collect::<String>();
        let fasta = format!(">chr1\n{sequence}\n");
        fs::write(root.join("authoring/references/ref.fa"), &fasta).expect("reference");
        let general = root.join("authoring/general.toml");
        fs::write(
            &general,
            format!(
                "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{SEED}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
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
                "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = {}\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = 1\nmax_swap_proposals = 1\nmax_haplotype_cells = {}\nblock_count = {}\n{subpopulations}",
                self.populations.len(),
                sha256(metadata.as_bytes()),
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

/// Cada cadena retira su fixture al terminar. Sin esto, cada ejecución de la
/// batería dejaba una docena de directorios en el temporal del sistema, que
/// se acumulan sin límite.
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

fn marginals_of(chain: &Chain) -> F5MarginalTable {
    build_f5_marginals(
        &chain.config,
        &chain.sealed,
        &chain.execution,
        &chain.layout,
    )
    .expect("tabla de marginales")
}

fn attempt(chain: &Chain, attempt_index: u64) -> F5HaplotypeAttempt {
    build_f5_haplotype_attempt(&chain.config, &marginals_of(chain), attempt_index)
        .expect("intento haplotípico")
}

// ---------------------------------------------------------------------------
// Oráculo independiente en enteros
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
fn scaled_weights(weights: &[f64]) -> (Vec<u128>, u32) {
    let scale = weights
        .iter()
        .map(|weight| dyadic(*weight).1)
        .max()
        .expect("al menos un peso")
        .max(TARGET_SCALE);
    let numerators = weights
        .iter()
        .map(|weight| {
            let (significand, exponent) = dyadic(*weight);
            significand << (scale - exponent)
        })
        .collect();
    (numerators, scale)
}

/// Permutación Fisher–Yates de DTC-089 reconstruida en `u128`.
fn oracle_permutation(length: usize, words: &mut impl Iterator<Item = u64>) -> Vec<usize> {
    let mut permutation = (0..length).collect::<Vec<_>>();
    for index in (1..length).rev() {
        let n = (index + 1) as u128;
        let threshold = ((1_u128 << 64) / n) * n;
        let other = loop {
            let value = u128::from(words.next().expect("palabra de permutación"));
            if value < threshold {
                break usize::try_from(value % n).expect("índice");
            }
        };
        permutation.swap(index, other);
    }
    permutation
}

/// Prefijo elegido: el anterior o el que cruza, el que minimice la distancia;
/// el empate conserva el anterior. Todo en enteros sobre `2^scale`.
fn oracle_prefix(
    permutation: &[usize],
    numerators: &[u128],
    target: u128,
    scale: u32,
    target_scale: u32,
) -> usize {
    let target = target << (scale - target_scale);
    let mut prefix = 0_u128;
    let mut previous = 0_u128;
    for (visited, row) in permutation.iter().enumerate() {
        prefix += numerators[*row];
        if prefix >= target {
            return if prefix - target < target - previous {
                visited + 1
            } else {
                visited
            };
        }
        previous = prefix;
    }
    permutation.len()
}

/// Marginal de una banda única `[lo, hi]`, como `f64` y como numerador exacto
/// sobre `2^54`.
///
/// Reproduce el mapa afín de DTC-089 con redondeo explícito ties-to-even, sin
/// llamar al generador productivo.
fn oracle_marginal(low: f64, high: f64, words: &mut impl Iterator<Item = u64>) -> (f64, u128) {
    if low == high {
        let (significand, exponent) = if low == 0.0 {
            (0_u128, 0_u32)
        } else {
            dyadic(low)
        };
        assert!(exponent <= TARGET_SCALE, "banda degenerada fuera de escala");
        return (low, significand << (TARGET_SCALE - exponent));
    }
    assert!(low == 0.25 && high == 0.75, "el oráculo cubre [1/4, 3/4]");
    let bits53 = u128::from(words.next().expect("palabra de valor") >> 11);
    // Valor exacto antes de redondear: (2^52 + bits53) / 2^54.
    let numerator = (1_u128 << 52) + bits53;
    if numerator < (1_u128 << 53) {
        (numerator as f64 / (1_u128 << 54) as f64, numerator)
    } else {
        let mut half = numerator / 2;
        if numerator % 2 == 1 && half % 2 == 1 {
            half += 1;
        }
        (half as f64 / (1_u128 << 53) as f64, half * 2)
    }
}

fn scenario_id() -> ScenarioId {
    ScenarioId::parse("scenario", "test").expect("scenario")
}

fn root_seed() -> RootSeed {
    RootSeed::parse("seed", SEED).expect("semilla")
}

/// Recorre un bloque×subpoblación con el oráculo y devuelve, por columna, el
/// marginal exacto, las filas alternativas esperadas y la polaridad esperada.
fn oracle_block(
    block: &F5BlockHaplotypes,
    weights: &[f64],
    band: (f64, f64),
    attempt_index: u64,
) -> Vec<(f64, usize, F5Polarity)> {
    let (numerators, scale) = scaled_weights(weights);
    let scenario = scenario_id();
    let mut expected = Vec::with_capacity(block.columns().len());
    for column in block.columns() {
        let mut marginal_rng = derive_f5_subseed(
            root_seed(),
            F5RandomDecision::marginal_spectrum(
                &scenario,
                block.population_id(),
                block.block_id(),
                column.site_id(),
            ),
        )
        .expect("subsemilla marginal")
        .local_rng();
        let mut marginal_words = std::iter::from_fn(|| Some(marginal_rng.next_u64()));
        let (target, target_numerator) = oracle_marginal(band.0, band.1, &mut marginal_words);

        let alternate = if target == 0.0 {
            0
        } else if target == 1.0 {
            weights.len()
        } else {
            let mut column_rng = derive_f5_subseed(
                root_seed(),
                F5RandomDecision::initial_column(
                    &scenario,
                    block.population_id(),
                    block.block_id(),
                    attempt_index,
                    column.site_index(),
                ),
            )
            .expect("subsemilla de columna")
            .local_rng();
            let mut column_words = std::iter::from_fn(|| Some(column_rng.next_u64()));
            let permutation = oracle_permutation(weights.len(), &mut column_words);
            oracle_prefix(
                &permutation,
                &numerators,
                target_numerator,
                scale,
                TARGET_SCALE,
            )
        };

        let polarity = if column.site_index() == 0 {
            F5Polarity::Positive
        } else {
            let mut polarity_rng = derive_f5_subseed(
                root_seed(),
                F5RandomDecision::polarity(
                    &scenario,
                    block.population_id(),
                    block.block_id(),
                    attempt_index,
                    column.site_index(),
                ),
            )
            .expect("subsemilla de polaridad")
            .local_rng();
            if u128::from(polarity_rng.next_u64() >> 11) < (1_u128 << 52) {
                F5Polarity::Positive
            } else {
                F5Polarity::Negative
            }
        };
        expected.push((target, alternate, polarity));
    }
    expected
}

/// Filas alternativas observadas en una columna, leídas de la matriz.
fn observed_alternate_rows(block: &F5BlockHaplotypes, column: usize) -> usize {
    (0..block.matrix().haplotype_count())
        .filter(|haplotype| {
            block
                .matrix()
                .allele(*haplotype, column)
                .expect("celda dentro de rango")
        })
        .count()
}

// ---------------------------------------------------------------------------
// Pruebas
// ---------------------------------------------------------------------------

#[test]
fn every_column_matches_the_integer_oracle_with_uniform_weights() {
    let root = temporary("uniform");
    let fixture = Fixture {
        populations: &["P2", "P1"],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let produced = attempt(&chain, 0);

    assert_eq!(produced.attempt_index(), 0);
    assert_eq!(
        produced.blocks().len(),
        4,
        "dos subpoblaciones × dos bloques"
    );
    for block in produced.blocks() {
        assert_eq!(block.matrix().haplotype_count(), 4);
        assert_eq!(block.matrix().site_count(), 3);
        assert_eq!(
            block.marginal_max_abs_error(),
            0.25,
            "max(q) con H=4 uniforme"
        );
        let expected = oracle_block(block, &fixture.weights, (0.25, 0.75), 0);
        assert_eq!(expected.len(), block.columns().len());
        for (column, (target, alternate, polarity)) in expected.iter().enumerate() {
            let produced_column = &block.columns()[column];
            assert_eq!(produced_column.target_alt_af(), *target, "columna {column}");
            assert_eq!(
                produced_column.alternate_haplotypes(),
                *alternate as u64,
                "columna {column}"
            );
            assert_eq!(
                observed_alternate_rows(block, column),
                *alternate,
                "la matriz y el recuento coinciden en la columna {column}"
            );
            // Con pesos uniformes diádicos `p_hat = k/H` es exacto.
            assert_eq!(
                produced_column.realized_alt_af(),
                *alternate as f64 / 4.0,
                "columna {column}"
            );
            assert_eq!(
                produced_column.abs_error(),
                (produced_column.realized_alt_af() - target).abs()
            );
            assert!(
                produced_column.abs_error() <= block.marginal_max_abs_error(),
                "columna {column}: la única cota marginal activa es max(q)"
            );
            assert_eq!(produced_column.polarity(), *polarity, "columna {column}");
            assert!(!produced_column.is_fixed());
        }
        assert!(
            block.max_abs_error() <= block.marginal_max_abs_error(),
            "el peor error del bloque respeta max(q)"
        );
    }

    // Orden canónico: subpoblación y después bloque, con la subpoblación
    // reordenada en el escenario.
    let keys = produced
        .blocks()
        .iter()
        .map(|block| {
            (
                block.population_id().as_str().to_owned(),
                block.block_index(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        keys,
        vec![
            ("P1".to_owned(), 0),
            ("P1".to_owned(), 1),
            ("P2".to_owned(), 0),
            ("P2".to_owned(), 1),
        ]
    );
    let population = PopulationId::parse("P2".to_owned(), "test").expect("population");
    assert!(produced.block(&population, 1).is_some());
    assert!(produced.block(&population, 7).is_none());
}

#[test]
fn non_uniform_and_non_dyadic_weights_match_the_integer_oracle() {
    let root = temporary("non-uniform");
    let fixture = Fixture {
        // Pesos que no son diádicos exactos: la compensación de Neumaier deja
        // de ser idénticamente cero, así que el recorrido de prefijos no
        // coincide con una suma ingenua por casualidad.
        weights: vec![0.1, 0.2, 0.3, 0.4],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let produced = attempt(&chain, 0);
    // El perfil normaliza los pesos; el oráculo parte de los normalizados, que
    // son una entrada, no una salida del algoritmo.
    let weights = chain.config.profile().haplotype_weights().to_vec();
    assert!(
        weights.iter().any(|weight| dyadic(*weight).1 > 54),
        "los pesos normalizados no son diádicos de escala corta"
    );

    for block in produced.blocks() {
        assert_eq!(
            block.marginal_max_abs_error(),
            weights.iter().copied().fold(0.0_f64, f64::max)
        );
        let expected = oracle_block(block, &weights, (0.25, 0.75), 0);
        for (column, (target, alternate, _)) in expected.iter().enumerate() {
            let produced_column = &block.columns()[column];
            assert_eq!(produced_column.target_alt_af(), *target, "columna {column}");
            assert_eq!(
                produced_column.alternate_haplotypes(),
                *alternate as u64,
                "columna {column}: el prefijo exacto en enteros"
            );
            assert!(
                produced_column.abs_error() <= block.marginal_max_abs_error(),
                "columna {column}"
            );
        }
    }
}

#[test]
fn an_exact_tie_keeps_the_previous_prefix_through_the_whole_chain() {
    let root = temporary("tie");
    let fixture = Fixture {
        // Banda degenerada en 3/8: con H=4 uniforme el objetivo equidista de
        // los prefijos 1/4 y 1/2.
        spectrum: "alt_af_min\talt_af_max\tweight\n0.375\t0.375\t1\n",
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let produced = attempt(&chain, 0);
    for block in produced.blocks() {
        let expected = oracle_block(block, &fixture.weights, (0.375, 0.375), 0);
        for (column, (target, alternate, _)) in expected.iter().enumerate() {
            assert_eq!(*target, 0.375);
            assert_eq!(*alternate, 1, "el empate conserva el prefijo anterior");
            assert_eq!(block.columns()[column].alternate_haplotypes(), 1);
            assert_eq!(block.columns()[column].realized_alt_af(), 0.25);
            assert_eq!(block.columns()[column].abs_error(), 0.125);
        }
    }
}

#[test]
fn fixed_targets_build_constant_columns_without_error() {
    let root = temporary("fixed");
    let fixture = Fixture {
        // El espectro fija 0 y 1 con un cuarto de peso cada uno.
        spectrum: "alt_af_min\talt_af_max\tweight\n0\t0\t0.25\n0.25\t0.75\t0.5\n1\t1\t0.25\n",
        weights: vec![0.5, 0.5],
        block_count: 4,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let produced = attempt(&chain, 0);

    let mut zero = 0_usize;
    let mut one = 0_usize;
    let mut variable = 0_usize;
    for block in produced.blocks() {
        for (column, produced_column) in block.columns().iter().enumerate() {
            let rows = observed_alternate_rows(block, column);
            if produced_column.target_alt_af() == 0.0 {
                zero += 1;
                assert_eq!(rows, 0, "p=0 da una columna constante de referencia");
                assert_eq!(produced_column.realized_alt_af(), 0.0);
                assert_eq!(produced_column.abs_error(), 0.0);
                assert!(produced_column.is_fixed());
            } else if produced_column.target_alt_af() == 1.0 {
                one += 1;
                assert_eq!(rows, 2, "p=1 da una columna constante alternativa");
                assert_eq!(produced_column.realized_alt_af(), 1.0);
                assert_eq!(produced_column.abs_error(), 0.0);
                assert!(produced_column.is_fixed());
            } else {
                variable += 1;
                assert!(!produced_column.is_fixed());
            }
            assert!(produced_column.abs_error() <= block.marginal_max_abs_error());
        }
    }
    assert!(
        zero > 0 && one > 0,
        "el fixture debe fijar columnas en 0 y en 1"
    );
    assert!(variable > 0, "y conservar columnas no fijadas");
}

#[test]
fn more_haplotypes_than_the_block_can_represent_is_rejected_before_the_attempt() {
    // Primera cota: H > 2^S con S=3 y H=16.
    let root = temporary("beyond-two-pow-s");
    let chain = build_chain(
        &root,
        &Fixture {
            weights: vec![0.0625; 16],
            ..Fixture::default()
        },
    );
    // La primera cota vive en los marginales, que es donde DTC-089 la sitúa:
    // antes de sortear. Desde que el intento recibe la tabla ya construida, el
    // rechazo ocurre al construirla, no al abrir el intento.
    let error = build_f5_marginals(
        &chain.config,
        &chain.sealed,
        &chain.execution,
        &chain.layout,
    )
    .expect_err("H=16 no cabe en 2^3 vectores");
    assert_eq!(error.code(), "GPS-DOMAIN-105");
    assert!(
        error.context()["detail"].contains("sitios"),
        "la primera cota cita los sitios del bloque: {}",
        error.context()["detail"]
    );

    // Segunda cota: todas las columnas quedan fijadas, de modo que V=0 y solo
    // existe un vector alélico posible.
    let root = temporary("beyond-two-pow-v");
    let chain = build_chain(
        &root,
        &Fixture {
            spectrum: "alt_af_min\talt_af_max\tweight\n0\t0\t1\n",
            weights: vec![0.5, 0.5],
            ..Fixture::default()
        },
    );
    let error = build_f5_haplotype_attempt(&chain.config, &marginals_of(&chain), 0)
        .expect_err("H=2 no cabe en 2^0 vectores");
    assert_eq!(error.code(), "GPS-DOMAIN-105");
    assert!(
        error.context()["detail"].contains("columnas no fijadas"),
        "la segunda cota cita las columnas libres: {}",
        error.context()["detail"]
    );
    assert_eq!(error.context()["population_id"], "P1");
}

#[test]
fn duplicate_rows_are_counted_and_never_merged() {
    let root = temporary("duplicates");
    let fixture = Fixture {
        // Con H=2 y p=0,1 el prefijo vacío gana a 1/2, así que toda columna
        // sale constante y las dos filas coinciden. El intento respeta la cota
        // marginal y aun así es inválido: F5.2.3 deberá deshacerlo con swaps o
        // rechazarlo.
        spectrum: "alt_af_min\talt_af_max\tweight\n0.1\t0.1\t1\n",
        weights: vec![0.5, 0.5],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let produced = attempt(&chain, 0);
    for block in produced.blocks() {
        for column in block.columns() {
            assert_eq!(column.target_alt_af(), 0.1);
            assert_eq!(column.alternate_haplotypes(), 0);
            assert_eq!(column.abs_error(), 0.1);
            assert!(column.abs_error() <= block.marginal_max_abs_error());
        }
        assert_eq!(
            block.duplicate_rows(),
            1,
            "dos filas idénticas cuentan un duplicado"
        );
        assert_eq!(block.matrix().row_bytes(0), block.matrix().row_bytes(1));
    }
}

#[test]
fn the_attempt_index_changes_the_matrix_because_it_changes_the_scope() {
    let root = temporary("attempt-index");
    let chain = build_chain(&root, &Fixture::default());
    let first = attempt(&chain, 0);
    let second = attempt(&chain, 1);
    assert_ne!(
        first.blocks()[0].matrix(),
        second.blocks()[0].matrix(),
        "dos intentos abren subsemillas distintas"
    );
    // El objetivo no depende del intento: los marginales son previos.
    for (left, right) in first.blocks()[0]
        .columns()
        .iter()
        .zip(second.blocks()[0].columns())
    {
        assert_eq!(left.target_alt_af(), right.target_alt_af());
    }
    // Y el segundo intento sigue reproduciendo el oráculo con su índice.
    let expected = oracle_block(&second.blocks()[0], &[0.25; 4], (0.25, 0.75), 1);
    for (column, (_, alternate, polarity)) in expected.iter().enumerate() {
        assert_eq!(
            second.blocks()[0].columns()[column].alternate_haplotypes(),
            *alternate as u64
        );
        assert_eq!(second.blocks()[0].columns()[column].polarity(), *polarity);
    }
}

/// Medida instrumental sobre la ruta productiva: con H=4096 el coste declarado
/// sigue siendo lineal en H. Una enumeración de `2^H` necesitaría `2^4096`
/// pasos y esta prueba no terminaría nunca.
#[test]
fn four_thousand_ninety_six_haplotypes_stay_linear_in_h() {
    let root = temporary("h4096");
    let haplotype_count = 4096_usize;
    let fixture = Fixture {
        block_count: 1,
        block_length: 12,
        weights: vec![1.0 / haplotype_count as f64; haplotype_count],
        max_haplotype_cells: 1 << 20,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let started = std::time::Instant::now();
    let produced = attempt(&chain, 0);
    let elapsed = started.elapsed();

    let block = &produced.blocks()[0];
    let site_count = block.matrix().site_count();
    assert_eq!(site_count, 12, "S=12 es el mínimo que admite H=2^12");
    assert_eq!(block.matrix().haplotype_count(), haplotype_count);
    assert_eq!(block.marginal_max_abs_error(), 1.0 / haplotype_count as f64);
    // El error marginal publicado se contrasta contra el máximo
    // `|p_hat_i - p_i|` calculado con enteros: reafirmar la cota que
    // `build_block` ya impone abortando no comprobaría nada que el propio
    // constructor no garantice.
    let (numerators, scale) = scaled_weights(chain.config.profile().haplotype_weights());
    let unit: u128 = numerators.iter().sum();
    assert_eq!(unit, 1_u128 << scale, "los pesos uniformes suman la unidad");
    let mut expected_max_abs_error = 0.0_f64;
    for (site, column) in block.columns().iter().enumerate() {
        let mut realized = 0_u128;
        for (haplotype, weight) in numerators.iter().enumerate() {
            if block.matrix().allele(haplotype, site).expect("alelo") {
                realized += weight;
            }
        }
        let error = (realized as f64 / unit as f64 - column.target_alt_af()).abs();
        assert_eq!(column.abs_error(), error, "el error de la columna {site}");
        if error > expected_max_abs_error {
            expected_max_abs_error = error;
        }
    }
    assert_eq!(
        block.max_abs_error(),
        expected_max_abs_error,
        "el máximo publicado no es el máximo de los errores por columna"
    );

    let cells = (haplotype_count * site_count) as u64;
    eprintln!(
        "medida productiva H={haplotype_count} S={site_count}: operaciones={} tiempo={elapsed:?}",
        block.construction_operations()
    );
    assert!(
        block.construction_operations() < 4 * cells,
        "H=4096 costó {} operaciones para {cells} celdas en {elapsed:?}",
        block.construction_operations()
    );
    // La codificación canónica ocupa H*ceil(S/8) bytes más la cabecera.
    assert_eq!(
        block.matrix().canonical_bytes().len(),
        16 + haplotype_count * 2
    );
}

/// La medida instrumental frente a la enumeración: al doblar H, el coste
/// declarado no llega a triplicarse.
#[test]
fn the_declared_cost_doubles_with_h_instead_of_squaring_the_search_space() {
    let mut previous: Option<u64> = None;
    for exponent in 4..=8_u32 {
        let haplotype_count = 1_usize << exponent;
        let root = temporary(&format!("ladder-{haplotype_count}"));
        let fixture = Fixture {
            block_count: 1,
            block_length: 12,
            weights: vec![1.0 / haplotype_count as f64; haplotype_count],
            max_haplotype_cells: 1 << 20,
            ..Fixture::default()
        };
        let chain = build_chain(&root, &fixture);
        let produced = attempt(&chain, 0);
        let operations = produced.blocks()[0].construction_operations();
        assert!(
            operations < 4 * (haplotype_count * 12) as u64,
            "H={haplotype_count} costó {operations} operaciones"
        );
        eprintln!("medida productiva H={haplotype_count} S=12: operaciones={operations}");
        if let Some(previous) = previous {
            assert!(
                operations < previous * 3,
                "al doblar H el coste pasó de {previous} a {operations}"
            );
        }
        previous = Some(operations);
    }
}

/// M-1: los marginales no dependen del intento —su ámbito no lleva
/// `attempt_index`— y la admisión los presupuesta una sola vez, sin factor de
/// intentos. La firma recibe la tabla ya construida precisamente para que
/// recorrer los intentos no rederive `B·N_pop·S` subsemillas por vuelta.
#[test]
fn several_attempts_share_one_marginal_table() {
    let root = temporary("shared-marginals");
    let chain = build_chain(&root, &Fixture::default());
    let marginals = marginals_of(&chain);

    let first = build_f5_haplotype_attempt(&chain.config, &marginals, 0).expect("primer intento");
    let second = build_f5_haplotype_attempt(&chain.config, &marginals, 1).expect("segundo intento");

    // La misma tabla sirve a los dos intentos y los objetivos son idénticos.
    for block in marginals.blocks() {
        let from_first = first
            .block(block.population_id(), block.block_index())
            .expect("bloque del primer intento");
        let from_second = second
            .block(block.population_id(), block.block_index())
            .expect("bloque del segundo intento");
        for (index, site) in block.sites().iter().enumerate() {
            assert_eq!(
                from_first.columns()[index].target_alt_af(),
                site.target_alt_af()
            );
            assert_eq!(
                from_second.columns()[index].target_alt_af(),
                site.target_alt_af()
            );
        }
    }
    // Y sin embargo las matrices difieren, porque el ámbito del intento sí
    // lleva `attempt_index`.
    assert_ne!(
        first.blocks()[0].matrix().canonical_bytes(),
        second.blocks()[0].matrix().canonical_bytes()
    );
    fs::remove_dir_all(&root).expect("cleanup");
}
