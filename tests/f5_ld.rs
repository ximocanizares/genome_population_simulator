//! F5.2.2. Evaluador LD: `p_hat`, `D`, `r²`, bins y signos con el orden y el
//! acumulador de Neumaier normativos.
//!
//! Las pruebas recorren la cadena productiva completa —configuración, admisión
//! sellada, layout, marginales, haplotipos y evaluación LD— y contrastan cada
//! pareja con un oráculo racional escrito en enteros de 128 bits: con pesos
//! uniformes `1/H` todo el cálculo de DTC-089 es un cociente de enteros
//! exactos, de modo que el oráculo no reutiliza ninguna función del evaluador.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::f5_haplotypes::F5BlockHaplotypes;
use genome_population_simulator::f5_haplotypes::{F5HaplotypeAttempt, build_f5_haplotype_attempt};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_ld::{
    F5LdEvaluation, evaluate_f5_linkage_disequilibrium, evaluate_f5_site_pair,
};
use genome_population_simulator::f5_marginals::build_f5_marginals;
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::identity::ExecutionId;
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

const EXECUTION: &str = "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4";
const SEED: &str = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8";

/// Espectro de una sola banda continua `[1/4, 3/4]`.
const CONTINUOUS_SPECTRUM: &str = "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t1\n";

fn sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-ld-{}-{}-{label}",
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
    spectrum: &'a str,
    /// Bins LD `(inicio, fin, media objetivo)`; deben cubrir `1..=L-1`.
    bins: Vec<(u64, u64, f64)>,
}

impl Default for Fixture<'_> {
    fn default() -> Self {
        Self {
            populations: &["P1"],
            block_count: 2,
            block_length: 4,
            weights: vec![0.25; 4],
            gap: 1,
            spectrum: CONTINUOUS_SPECTRUM,
            bins: vec![(1, 3, 0.5)],
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
            for (start, end, target) in &self.bins {
                ld_rows.push_str(&format!("{population}\t{start}\t{end}\t{target}\t0.05\n"));
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
                format!("gap_min_bp\tgap_max_bp\tweight\n{0}\t{0}\t1\n", self.gap),
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
                "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = {}\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = 1\nmax_swap_proposals = 1\nmax_haplotype_cells = 65536\nblock_count = {}\n{subpopulations}",
                self.populations.len(),
                sha256(metadata.as_bytes()),
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

/// Cada cadena retira su fixture al terminar: la batería no deja directorios
/// temporales acumulándose.
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

fn attempt_of(chain: &Chain) -> F5HaplotypeAttempt {
    let marginals = build_f5_marginals(
        &chain.config,
        &chain.sealed,
        &chain.execution,
        &chain.layout,
    )
    .expect("tabla de marginales");
    build_f5_haplotype_attempt(&chain.config, &marginals, 0).expect("intento haplotípico")
}

fn evaluation_of(chain: &Chain, attempt: &F5HaplotypeAttempt) -> F5LdEvaluation {
    evaluate_f5_linkage_disequilibrium(&chain.config, &chain.layout, attempt)
        .expect("evaluación LD")
}

// ---------------------------------------------------------------------------
// Oráculo racional exacto
// ---------------------------------------------------------------------------

/// Recuentos enteros `(N_i, N_j, N_ij)` de una pareja de columnas.
fn counts(block: &F5BlockHaplotypes, i: usize, j: usize) -> (u128, u128, u128) {
    let matrix = block.matrix();
    let mut n_i = 0_u128;
    let mut n_j = 0_u128;
    let mut n_ij = 0_u128;
    for haplotype in 0..matrix.haplotype_count() {
        let alt_i = matrix.allele(haplotype, i).expect("alelo i");
        let alt_j = matrix.allele(haplotype, j).expect("alelo j");
        n_i += u128::from(alt_i);
        n_j += u128::from(alt_j);
        n_ij += u128::from(alt_i && alt_j);
    }
    (n_i, n_j, n_ij)
}

/// Oráculo racional exacto para `H` pesos uniformes `1/H`.
///
/// Sobre esa escala todo es entero: `p_i=N_i/H`, `p_ij=N_ij/H`,
/// `D=(H*N_ij-N_i*N_j)/H²` y `r²=(H*N_ij-N_i*N_j)²/[N_i(H-N_i)*N_j(H-N_j)]`,
/// donde numerador y denominador comparten `H⁴` y se cancelan. Con `H=4` cada
/// cociente intermedio es un múltiplo exacto de `1/256`, representable en
/// binary64 sin redondeo, de modo que la comparación puede ser una igualdad y
/// no una tolerancia. No llama a ninguna función del evaluador.
fn rational_oracle(
    block: &F5BlockHaplotypes,
    i: usize,
    j: usize,
) -> (f64, f64, f64, f64, Option<f64>) {
    let haplotypes = u128::try_from(block.matrix().haplotype_count()).expect("H");
    let (n_i, n_j, n_ij) = counts(block, i, j);
    let scale = i128::try_from(haplotypes).expect("H");
    let d_numerator =
        scale * i128::try_from(n_ij).expect("N_ij") - i128::try_from(n_i * n_j).expect("N_i*N_j");
    let denominator = n_i * (haplotypes - n_i) * n_j * (haplotypes - n_j);
    let r_squared = (denominator > 0).then(|| {
        let squared = d_numerator * d_numerator;
        squared as f64 / denominator as f64
    });
    let haplotypes = haplotypes as f64;
    (
        n_i as f64 / haplotypes,
        n_j as f64 / haplotypes,
        n_ij as f64 / haplotypes,
        d_numerator as f64 / (haplotypes * haplotypes),
        r_squared,
    )
}

/// Suma de Neumaier escrita desde su definición, para reducir en el orden
/// canónico de parejas las `r²` que produjo el oráculo.
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
        .expect("al menos un peso");
    let numerators = weights
        .iter()
        .map(|weight| {
            let (significand, exponent) = dyadic(*weight);
            significand << (scale - exponent)
        })
        .collect();
    (numerators, scale)
}

/// El `f64` correctamente redondeado de `numerator / 2^scale`.
///
/// Compararse con este valor es **más estricto** que lo que DTC-089 garantiza:
/// el contrato fija el orden y el acumulador, no que la suma compensada sea el
/// redondeo correcto de la fracción exacta. Con los pesos de esta prueba lo es;
/// con otros podría no serlo sin que hubiera defecto, y entonces la
/// comparación tendría que pasar a «a un ULP».
///
/// El casting de `u128` a `f64` redondea al más próximo y la división por una
/// potencia de dos es exacta, así que la composición es el redondeo único de
/// la fracción exacta.
fn rounded(numerator: u128, scale: u32) -> f64 {
    numerator as f64 / 2.0_f64.powi(i32::try_from(scale).expect("escala"))
}

// ---------------------------------------------------------------------------
// Pruebas
// ---------------------------------------------------------------------------

/// Bins de una unidad de ancho: cada distancia intrabloque tiene el suyo.
fn unit_bins() -> Vec<(u64, u64, f64)> {
    vec![(1, 1, 0.5), (2, 2, 0.4), (3, 3, 0.3)]
}

/// Cada pareja y cada bin coinciden con el oráculo racional exacto.
#[test]
fn pairs_and_bins_reproduce_an_exact_rational_oracle() {
    let root = temporary("oracle");
    let fixture = Fixture {
        bins: unit_bins(),
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);
    let evaluation = evaluation_of(&chain, &attempt);
    let weights = chain.config.profile().haplotype_weights();
    assert_eq!(evaluation.attempt_index(), attempt.attempt_index());
    assert_eq!(evaluation.blocks().len(), attempt.blocks().len());

    for block in attempt.blocks() {
        let evaluated = evaluation
            .block(block.population_id(), block.block_index())
            .expect("bloque evaluado");
        let site_count = block.matrix().site_count();
        let mut per_bin: Vec<Vec<f64>> = vec![Vec::new(); evaluated.bins().len()];
        let mut evaluable = 0_u64;
        let mut not_evaluable = 0_u64;
        let mut signs_failed = 0_u64;

        for i in 0..site_count {
            for j in (i + 1)..site_count {
                let pair =
                    evaluate_f5_site_pair(weights, block.matrix(), i, j).expect("pareja evaluada");
                let (p_i, p_j, p_ij, d, r_squared) = rational_oracle(block, i, j);
                assert_eq!(pair.p_hat_i(), p_i, "p_hat_i de ({i},{j})");
                assert_eq!(pair.p_hat_j(), p_j, "p_hat_j de ({i},{j})");
                assert_eq!(pair.p_hat_ij(), p_ij, "p_hat_ij de ({i},{j})");
                assert_eq!(pair.d(), d, "D de ({i},{j})");
                assert_eq!(pair.r_squared(), r_squared, "r² de ({i},{j})");

                let Some(value) = r_squared else {
                    not_evaluable += 1;
                    continue;
                };
                evaluable += 1;
                let matched = block.columns()[i].polarity() == block.columns()[j].polarity();
                if !(if matched { d > 0.0 } else { d < 0.0 }) {
                    signs_failed += 1;
                }
                // Las distancias del fixture son `j-i`: los sitios van de uno
                // en uno desde el inicio del bloque. La aserción lo comprueba
                // en lugar de suponerlo.
                let layout_sites = chain.layout.blocks()[block.block_index() as usize].sites();
                let distance =
                    layout_sites[j].position_1based() - layout_sites[i].position_1based();
                assert_eq!(
                    distance,
                    (j - i) as u64,
                    "geometría supuesta por el oráculo"
                );
                let bin = evaluated
                    .bins()
                    .iter()
                    .position(|bin| {
                        bin.distance_start_bp() <= distance && distance <= bin.distance_end_bp()
                    })
                    .expect("bin de la distancia");
                per_bin[bin].push(value);
            }
        }

        assert_eq!(evaluated.evaluable_pairs(), evaluable);
        assert_eq!(evaluated.not_evaluable_pairs(), not_evaluable);
        assert_eq!(evaluated.signs_failed(), signs_failed);
        assert!(evaluated.is_valid());

        let mut expected_outside = 0_u64;
        let mut expected_max = 0.0_f64;
        let mut expected_errors = Vec::new();
        for (bin, values) in evaluated.bins().iter().zip(&per_bin) {
            assert_eq!(bin.pair_count(), values.len() as u64, "parejas del bin");
            if values.is_empty() {
                assert_eq!(bin.mean_r2(), None);
                assert_eq!(bin.abs_error(), None);
                assert_eq!(bin.status(), "not_applicable");
                continue;
            }
            let mean = neumaier(values) / values.len() as f64;
            assert_eq!(bin.mean_r2(), Some(mean), "media del bin");
            let error = (mean - bin.target_mean_r2()).abs();
            assert_eq!(bin.abs_error(), Some(error), "error del bin");
            assert_eq!(bin.status(), "evaluated");
            assert_eq!(bin.is_within_tolerance(), error <= bin.max_abs_error());
            if error > bin.max_abs_error() {
                expected_outside += 1;
            }
            expected_max = expected_max.max(error);
            expected_errors.push(error);
        }
        assert_eq!(evaluated.bins_outside(), expected_outside);
        assert_eq!(evaluated.max_r2_error(), expected_max);
        assert_eq!(evaluated.sum_r2_error(), neumaier(&expected_errors));
    }
    assert!(evaluation.is_valid());
}

/// Los bins solo reciben parejas intrabloque y respetan sus extremos.
///
/// Con dos bloques contiguos de cuatro sitios consecutivos, las distancias
/// interbloque `1`, `2` y `3` existen físicamente —los sitios 4 y 5, 4 y 6,
/// 4 y 7— y caerían en los mismos bins si el evaluador las considerase. Los
/// recuentos exactos `3`, `2` y `1` por bloque demuestran que no lo hace.
#[test]
fn each_bin_receives_only_the_intrablock_pairs_at_its_distance() {
    let root = temporary("intrablock");
    let fixture = Fixture {
        bins: unit_bins(),
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);
    let evaluation = evaluation_of(&chain, &attempt);

    let positions = chain
        .layout
        .blocks()
        .iter()
        .flat_map(|block| block.sites().iter().map(|site| site.position_1based()))
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        vec![1, 2, 3, 4, 5, 6, 7, 8],
        "geometría del fixture"
    );

    assert_eq!(evaluation.blocks().len(), 2);
    for block in evaluation.blocks() {
        let counts = block
            .bins()
            .iter()
            .map(|bin| bin.pair_count())
            .collect::<Vec<_>>();
        assert_eq!(counts, vec![3, 2, 1], "solo las parejas del propio bloque");
        assert_eq!(counts.iter().sum::<u64>(), 6, "C(4,2) parejas intrabloque");
        assert_eq!(block.evaluable_pairs(), 6);
    }
}

/// Un bin sin ninguna pareja evaluable es `not_applicable`, nunca un cero.
#[test]
fn a_bin_without_pairs_is_not_applicable_and_never_zero() {
    let root = temporary("empty-bin");
    let fixture = Fixture {
        gap: 2,
        bins: unit_bins(),
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);
    let evaluation = evaluation_of(&chain, &attempt);

    for block in evaluation.blocks() {
        let bins = block.bins();
        assert_eq!(bins.len(), 3);
        for empty in [&bins[0], &bins[2]] {
            assert_eq!(empty.status(), "not_applicable");
            assert_eq!(empty.pair_count(), 0);
            assert_eq!(empty.mean_r2(), None, "un bin vacío no vale cero");
            assert_eq!(empty.abs_error(), None);
            assert!(empty.is_within_tolerance(), "un bin vacío no es un error");
            assert!(!empty.is_evaluated());
        }
        assert_eq!(bins[1].status(), "evaluated");
        assert_eq!(bins[1].pair_count(), 1);
        assert!(bins[1].mean_r2().is_some());
        assert!(block.is_valid(), "un bin vacío no invalida el bloque");
    }
    assert!(evaluation.is_valid());
}

/// Los sitios fijados no entran en ningún bin y pueden invalidar el intento.
///
/// El espectro mezcla una banda degenerada en 1 con una continua: el primer
/// bloque recibe tres sitios fijados y uno variable, de modo que sus seis
/// parejas son `not_evaluable` y todos sus bins `not_applicable`. El segundo
/// bloque no tiene ningún sitio fijado y evalúa sus seis parejas.
#[test]
fn fixed_sites_are_not_evaluable_and_can_invalidate_the_attempt() {
    let root = temporary("fixed");
    let fixture = Fixture {
        weights: vec![0.5; 2],
        spectrum: "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t0.5\n1\t1\t0.5\n",
        bins: unit_bins(),
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);
    let evaluation = evaluation_of(&chain, &attempt);

    let first = &attempt.blocks()[0];
    assert_eq!(
        first.columns().iter().filter(|c| c.is_fixed()).count(),
        3,
        "el fixture debe dejar tres sitios fijados en el primer bloque"
    );
    let evaluated = evaluation
        .block(first.population_id(), first.block_index())
        .expect("bloque evaluado");
    assert_eq!(evaluated.evaluable_pairs(), 0);
    assert_eq!(evaluated.not_evaluable_pairs(), 6);
    assert!(
        evaluated.bins().iter().all(|bin| !bin.is_evaluated()),
        "sin pareja evaluable ningún bin puede evaluarse"
    );
    assert!(!evaluated.is_valid());
    assert!(
        !evaluation.is_valid(),
        "un bloque sin pareja evaluable invalida el intento"
    );

    let second = &attempt.blocks()[1];
    assert!(second.columns().iter().all(|c| !c.is_fixed()));
    let evaluated = evaluation
        .block(second.population_id(), second.block_index())
        .expect("bloque evaluado");
    assert_eq!(evaluated.evaluable_pairs(), 6);
    assert!(evaluated.is_valid());
}

/// Una marginal realizada en 0 deja la pareja `not_evaluable` aunque su
/// objetivo no estuviera fijado: el estado es del denominador, no del perfil.
#[test]
fn a_realized_zero_marginal_is_not_evaluable_without_being_fixed() {
    let root = temporary("realized-zero");
    let fixture = Fixture {
        weights: vec![0.5; 2],
        spectrum: "alt_af_min\talt_af_max\tweight\n0.01\t0.1\t1\n",
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);
    let evaluation = evaluation_of(&chain, &attempt);

    for block in attempt.blocks() {
        assert!(
            block.columns().iter().all(|column| !column.is_fixed()),
            "ningún objetivo está fijado"
        );
        assert!(
            block
                .columns()
                .iter()
                .all(|column| column.realized_alt_af() == 0.0),
            "toda columna realiza 0"
        );
    }
    for block in evaluation.blocks() {
        assert_eq!(block.evaluable_pairs(), 0);
        assert_eq!(block.not_evaluable_pairs(), 6);
        assert_eq!(
            block.signs_failed(),
            0,
            "sin pareja evaluable no hay signos"
        );
        assert_eq!(block.bins()[0].status(), "not_applicable");
        assert!(!block.is_valid());
    }
    assert!(!evaluation.is_valid());
}

/// Con pesos no diádicos la compensación de Neumaier deja de ser cero: las
/// marginales simple y conjunta se comparan contra la fracción exacta sobre
/// `2^scale` y no contra una suma ingenua.
#[test]
fn non_dyadic_weights_match_the_exact_rational_marginals() {
    let root = temporary("non-dyadic");
    let fixture = Fixture {
        weights: vec![0.1; 10],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);
    let weights = chain.config.profile().haplotype_weights();
    let (numerators, scale) = scaled_weights(weights);
    assert!(
        numerators.iter().any(|numerator| numerator % (1 << 3) != 0),
        "los pesos deben ser inexactos en binario para que la compensación cuente"
    );

    let mut seen_compensation = false;
    for block in attempt.blocks() {
        let matrix = block.matrix();
        for i in 0..matrix.site_count() {
            for j in (i + 1)..matrix.site_count() {
                let pair = evaluate_f5_site_pair(weights, matrix, i, j).expect("pareja");
                let mut exact_i = 0_u128;
                let mut exact_j = 0_u128;
                let mut exact_ij = 0_u128;
                let mut naive_i = 0.0_f64;
                for (haplotype, numerator) in numerators.iter().enumerate() {
                    let alt_i = matrix.allele(haplotype, i).expect("alelo i");
                    let alt_j = matrix.allele(haplotype, j).expect("alelo j");
                    if alt_i {
                        exact_i += numerator;
                        naive_i += weights[haplotype];
                    }
                    if alt_j {
                        exact_j += numerator;
                    }
                    if alt_i && alt_j {
                        exact_ij += numerator;
                    }
                }
                assert_eq!(pair.p_hat_i(), rounded(exact_i, scale), "p_hat_i ({i},{j})");
                assert_eq!(pair.p_hat_j(), rounded(exact_j, scale), "p_hat_j ({i},{j})");
                assert_eq!(
                    pair.p_hat_ij(),
                    rounded(exact_ij, scale),
                    "p_hat_ij ({i},{j})"
                );
                seen_compensation |= naive_i != pair.p_hat_i();
            }
        }
    }
    assert!(
        seen_compensation,
        "alguna columna debe distinguir la suma compensada de la ingenua"
    );
}

/// Un layout que no describe la misma geometría que el intento se rechaza.
#[test]
fn a_layout_with_another_geometry_is_refused() {
    let mine = temporary("geometry-mine");
    let other = temporary("geometry-other");
    let chain = build_chain(&mine, &Fixture::default());
    let foreign = build_chain(
        &other,
        &Fixture {
            block_count: 1,
            block_length: 8,
            bins: vec![(1, 7, 0.5)],
            ..Fixture::default()
        },
    );
    let attempt = attempt_of(&chain);
    let diagnostic = evaluate_f5_linkage_disequilibrium(&chain.config, &foreign.layout, &attempt)
        .expect_err("el layout ajeno debe rechazarse");
    assert_eq!(diagnostic.code(), "GPS-CONF-101");
    assert_eq!(
        diagnostic.context().get("fields").map(String::as_str),
        Some("layout.identity")
    );
}

/// Regresión G-1 sobre la cadena productiva: un perfil legal cuyos pesos
/// normalizados suman `1.0000000000000002` no aborta la evaluación.
///
/// `normalize_slice` solo exige `|suma-1| <= 1e-9` después de dividir, así que
/// estos cuatro pesos —admitidos sin reproche— dejan una columna monomórfica
/// con `realized_alt_af` por encima de 1. F5.2.1 la acepta, porque su cota es
/// `|realizada - objetivo| <= max(q)` y el exceso es de un ULP; el evaluador
/// debe declarar sus parejas `not_evaluable`, no matar la ejecución.
#[test]
fn a_profile_whose_weights_overshoot_one_does_not_abort_the_evaluation() {
    let root = temporary("overshoot");
    let fixture = Fixture {
        block_count: 1,
        block_length: 8,
        weights: vec![0.019251, 0.383455, 0.511446, 0.085848],
        spectrum: "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t0.5\n1\t1\t0.5\n",
        bins: vec![(1, 7, 0.5)],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);

    let block = &attempt.blocks()[0];
    let overshooting = block
        .columns()
        .iter()
        .filter(|column| column.realized_alt_af() > 1.0)
        .count();
    assert!(
        overshooting > 0,
        "el fixture debe producir alguna marginal por encima de 1"
    );

    let evaluation = evaluation_of(&chain, &attempt);
    let evaluated = evaluation
        .block(block.population_id(), block.block_index())
        .expect("bloque evaluado");
    for bin in evaluated.bins() {
        if let Some(mean) = bin.mean_r2() {
            assert!(
                (0.0..=1.0).contains(&mean),
                "media fuera del dominio: {mean}"
            );
        }
    }
    // Cada pareja que toca una columna monomórfica es `not_evaluable`: sus dos
    // factores serían negativos y su producto, positivo.
    let expected_not_evaluable = block
        .columns()
        .iter()
        .enumerate()
        .flat_map(|(i, left)| {
            block
                .columns()
                .iter()
                .skip(i + 1)
                .map(move |right| (left, right))
        })
        .filter(|(left, right)| {
            !(0.0..1.0).contains(&left.realized_alt_af())
                || !(0.0..1.0).contains(&right.realized_alt_af())
        })
        .count() as u64;
    assert!(expected_not_evaluable > 0);
    assert_eq!(evaluated.not_evaluable_pairs(), expected_not_evaluable);
}

/// Dos bins anchos parten exactamente en su frontera `end_k` / `start_{k+1}`.
///
/// Con ocho sitios consecutivos, la distancia `k` aparece `8-k` veces: el bin
/// `1..3` recibe `7+6+5` parejas y el `4..7`, `4+3+2+1`. Un error de un paso en
/// `bin_of` movería una pareja de un lado al otro y rompería ambos recuentos.
#[test]
fn wide_bins_split_exactly_at_their_boundary() {
    let root = temporary("wide-bins");
    let fixture = Fixture {
        block_count: 1,
        block_length: 8,
        bins: vec![(1, 3, 0.5), (4, 7, 0.3)],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);
    let evaluation = evaluation_of(&chain, &attempt);

    let block = &evaluation.blocks()[0];
    let counts = block
        .bins()
        .iter()
        .map(|bin| bin.pair_count())
        .collect::<Vec<_>>();
    assert_eq!(counts, vec![18, 10], "reparto exacto en la frontera 3|4");
    assert_eq!(counts.iter().sum::<u64>(), 28, "C(8,2) parejas intrabloque");
}

/// Dos subpoblaciones se evalúan por separado, cada una con sus propios bins.
#[test]
fn two_subpopulations_are_evaluated_independently() {
    let root = temporary("two-populations");
    let fixture = Fixture {
        populations: &["P1", "P2"],
        bins: unit_bins(),
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let attempt = attempt_of(&chain);
    let evaluation = evaluation_of(&chain, &attempt);

    assert_eq!(evaluation.blocks().len(), 4, "dos bloques por subpoblación");
    let populations = evaluation
        .blocks()
        .iter()
        .map(|block| block.population_id().as_str().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(populations, vec!["P1", "P1", "P2", "P2"]);

    for block in evaluation.blocks() {
        assert_eq!(
            block
                .bins()
                .iter()
                .map(|bin| bin.pair_count())
                .collect::<Vec<_>>(),
            vec![3, 2, 1],
            "cada bloque×subpoblación agrega solo sus propias parejas"
        );
    }
    // Las matrices de las dos subpoblaciones son independientes: sus ámbitos
    // RNG llevan el `population_id`, así que sus medias no tienen por qué
    // coincidir y de hecho no coinciden.
    let first = evaluation.blocks()[0].bins()[0].mean_r2();
    let second = evaluation.blocks()[2].bins()[0].mean_r2();
    assert_ne!(first, second, "las subpoblaciones no comparten haplotipos");
}
